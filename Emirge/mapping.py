"""Read mapping module

This module abstracts from the various read mapping tools.
"""
import optparse
import os
import re
from multiprocessing import cpu_count
from subprocess import CalledProcessError

import Emirge.log as log
from Emirge.io import EnumerateReads, filter_fastq, decompressed, make_pipe, \
    check_call, TempDir, File, PIPE, command_avail, \
    fastq_count_reads, AlignmentFile
from Emirge.log import INFO, DEBUG, WARNING


# comment the two lines below to enable DEBUG info for this module
def DEBUG(*args):
    pass


class MappingFailure(Exception):
    """
    Something went wrong during read mapping
    """


class IndexBuildFailure(Exception):
    """
    Exception raised when something went wrong trying to build
    the mapping index
    """
    pass


# sam flags:
# 0x001 read paired
# 0x002 read mapped proper
# 0x004 read unmapped
# 0x008 mate unmapped
# 0x010 read reverse strand
# 0x020 mate reverse strand
# 0x040 first in pair
# 0x080 second in pair
# 0x100 not primary alignment
# 0x200 read fails platform/vendor quality checks
# 0x400 read is PCR or optical duplicate
# 0x800 supplementary alignment

Sam2Bam_aligned_only = make_pipe("Sam2Bam", [
    "samtools", "view",
    "-h",  # include header
    "-u",  # uncompressed output
    "-b",  # output bam
    "-F", "0x0004",  # exclude unmapped reads
    "-",  # output to stdout
])


class Mapper(object):
    """Base class for read mapping.

    Props:
        candidates -- fasta file containing reference sequences
        fwd_reads  -- forward reads file (may be changed by filtering)
        rev_reads  -- None or as fwd_reads
        phred33    -- use phred33 quality encoding if true, phred64 otherwise
        threads    -- number of threads mappers may use
        reindex    -- True if reads should be renamed to consecutive numbers
        workdir    -- directory to hold generated files
        insert_mean -- mean insert size
        insert_sd   -- std dev of insert size
        do_prefiltering -- if true, pre-filter reads using candidates
        n_reads -- number of reads
    """
    binary = None

    def __init__(self, candidates, candidates_index, fwd_reads, rev_reads=None,
                 phred33=False,
                 threads=cpu_count(), reindex=True, workdir=None,
                 insert_mean=1500, insert_sd=500, prefilter_reads=False):
        """Create new mapper.

        Args:
            candidates -- fasta file containing reference sequences
            fwd_reads  -- forward reads file
            rev_reads  -- reverse reads file, may be none for single ended
            phred33    -- true if quality encoded as phred33, false for 64
            threads    -- number of threads mappers should use
            reindex    -- change read names to consecutive numbers
            workdir    -- (optional) designed work directory
            insert_mean -- mean insert size
            insert_sd   -- std dev of insert size
            prefilter_reads -- if true, pre-filter reads using candidates
        """
        super(Mapper, self).__init__()

        if isinstance(fwd_reads, str):
            fwd_reads = decompressed(fwd_reads)
        if isinstance(rev_reads, str):
            rev_reads = decompressed(rev_reads)

        self.candidates = candidates
        self.candidates_index = candidates_index
        self.phred33 = phred33
        self.threads = threads
        self.insert_mean = insert_mean
        self.insert_sd = insert_sd
        self.do_prefiltering = prefilter_reads
        self.n_reads = None  # filled by prefilter reads

        if reindex:
            fwd_reads = EnumerateReads(fwd_reads)
            if rev_reads:
                rev_reads = EnumerateReads(rev_reads)

        self.fwd_reads, self.rev_reads = fwd_reads, rev_reads

        if workdir is None:
            self.tmp_workdir = TempDir()
            workdir = self.tmp_workdir.name
        self.workdir = workdir

    @classmethod
    def available(cls):
        """Check if pre-requisites for this mapper are met by system"""
        try:
            return cls.__available
        except AttributeError:
            cls.__available = command_avail(cls.binary)
            return cls.__available

    def prepare_reads(self):
        if self.do_prefiltering:
            self.prefilter_reads()
        else:
            self.n_reads = fastq_count_reads(self.fwd_reads)

    def prefilter_reads(self):
        raise NotImplementedError()


class Bowtie2(Mapper):
    """Handle read mapping with Bowtie2."""
    binary = "bowtie2"

    def __init__(self, *args, **kwargs):
        super(Bowtie2, self).__init__(*args, **kwargs)

    def make_basecmd(self):
        """Create array with common arguments to run bowtie2."""
        cmd = [
            self.binary,
            "-p", str(self.threads),
            "-t",  # print timings
        ]

        if self.rev_reads:
            cmd += [
                "-1", self.fwd_reads,
                "-2", self.rev_reads,
            ]
        else:
            cmd += [
                "-U", self.fwd_reads,
            ]

        if self.phred33:
            cmd += ["--phred33"]
        else:
            cmd += ["--phred64"]

        return cmd

    @log.timed("Pre-Filtering reads")
    def prefilter_reads(self, dirname="pre_mapping"):
        """Pre-filter reads using candidates.

        Switches the fwd/rev_reads property to a temporary file holding only
        reads mapping (loosely) to the candidates file
        """

        # prepare bowtie2 command
        cmd = self.make_basecmd()
        cmd += [
            "-x", self.prep_index(self.candidates, self.candidates_index),
            "--very-sensitive-local",
            "-k", "1",  # report up to 1 alignments per read
            "--no-unal",  # suppress SAM records for unaligned reads
        ]
        bowtie_cmd = make_pipe("prefilter", cmd)

        # run bowtie2 prefiltering command, use pysam to parse matching
        # read names from output into set 'keepers'
        keepers = set()
        i = 0
        with AlignmentFile(bowtie_cmd()) as sam:
            for read in sam.fetch(until_eof=True):
                if not read.is_unmapped:
                    # if we have reverse, query name is reported w/o "/1"
                    # or "/2" at the end, if we have only one read file,
                    # the full identifier is kept, so remove everything
                    # after a "/" before adding to keepers set:
                    keepers.add(read.query_name.split("/")[0])
                i += 1

        INFO("Pre-Mapping: mapper found %i matches" % len(keepers))

        for i, keeper in zip(range(10), keepers):
            DEBUG("keeper {}: '{}'".format(i, keeper))

        # replace forward read file with temporary file containing only
        # matching reads
        pathname = os.path.join(self.workdir, dirname)
        os.mkdir(pathname)
        pathname = os.path.join(pathname, "pre_mapped_reads.{}.fastq")

        with self.fwd_reads.reader() as reads, \
                File(pathname.format("1")).writer() as outfile:
            self.fwd_reads, fwd_count, fwd_matches = \
                filter_fastq(reads, keepers, outfile=outfile)

        # if we have a reverse read file, do the same for that
        if self.rev_reads is not None:
            with self.rev_reads.reader() as reads, \
                    File(pathname.format("2")).writer() as outfile:
                self.rev_reads, rev_count, rev_matches = \
                    filter_fastq(reads, keepers, outfile=outfile)

            # number of reads and must be identical for forward and reverse
            assert fwd_count == rev_count
            assert fwd_matches == rev_matches

        INFO("Pre-Mapping: {} out of {} reads remaining"
             .format(fwd_matches, fwd_count))
        self.n_reads = fwd_matches

    @log.timed('Mapping reads to "{fastafile}"')
    def map_reads(self, fastafile, workdir):
        """Create bam file containing mappings for fastafile in workdir.

        :param fastafile:   target sequences
        :param workdir:     directory to hold bamfile
        """
        minins = max((self.insert_mean - 3 * self.insert_sd), 0)
        maxins = self.insert_mean + 3 * self.insert_sd

        # use pre-computed index for candidates if available
        if fastafile == self.candidates:
            index = self.candidates_index
        else:
            index = None


        cmd = self.make_basecmd()
        cmd += [
            "-x", self.prep_index(fastafile, index),  # index
            "-D", "20",  # number of extension attempts (15)
            "-R", "3",  # try 3 sets of seeds (2)
            "-N", "0",  # max 0 mismatches, can be 0 or 1 (0)
            "-L", "20",  # length of seed substrings 4-31 (22)
            "-i", "S,1,0.50",  # interval between seed substrings (S,1,1.15)
            "-k", "20",  # report 20 alns per read
            "--no-unal",  # suppress sam output for unaligned reads
        ]

        if self.rev_reads is not None:
            cmd += [
                "--minins", str(minins),  # minimum fragment length
                "--maxins", str(maxins),  # maximum fragment length
                "--no-mixed",  # suppress unpaired alignments
                "--no-discordant",  # suppress discordant alignments
            ]

        bowtie_mapper = make_pipe("Bowtie Mapper", cmd)
        bam_file = File(os.path.join(workdir, "bt2.bam"))
        bam_writer = Sam2Bam_aligned_only(stdin=PIPE, stdout=bam_file)
        frags_mapped = 0

        with bowtie_mapper(stdout=bam_writer, stderr=PIPE) as stderr:
            for line in stderr:
                INFO("BOWTIE2: " + line.rstrip())
                match = re.search("\s([0-9]+) .* aligned .*1 time", line)
                if match is not None:
                    frags_mapped += int(match.group(1))
                match = re.search("Error", line)
                if match is not None:
                    raise MappingFailure("Bowtie2 failed: '{}'".format(line))

        INFO("Bowtie 2 reported {} alignments written to '{}'".format(
                frags_mapped, bam_file.name
        ))

        # TODO:delete index

        return frags_mapped, bam_file.reader()

    @staticmethod
    def have_index(indexname, fastafile=None):
        """Check if an index indexname exists and is newer than fastafile.

        :param indexname: name of bowtie2 index
        :param fastafile: name of corresponding fasta file
        """
        suffixes = ["1.bt2", "2.bt2", "3.bt2", "4.bt2",
                    "rev.1.bt2", "rev.2.bt2"]
        no_access = []
        too_old = []
        for suffix in suffixes:
            file_name = indexname + "." + suffix
            if not os.access(file_name, os.R_OK):
                no_access.append(file_name)
            elif fastafile is not None:
                if os.path.getctime(file_name) < os.path.getctime(fastafile):
                    too_old.append(file_name)

        if len(no_access) == 0 and len(too_old) == 0:
            return True
        if len(no_access) == len(suffixes):
            return False
        if len(too_old) > 0:
            WARNING('Bowtie2 index for "{}" is out of date'.format(fastafile))
            return False
        if len(no_access) > 0:
            WARNING('No (or incomplete) Bowtie2 index found with basename "{'
                    '}"'.format(indexname))
            return False

    @log.timed("Preparing index for {fastafile} ({indexname})")
    def prep_index(self, fastafile, indexname=None):
        """Prepare index for fastafile. If indexname given, use that name.

        Check if an index exists for fastafile, and if not, build one. If
        indexname is given, the index must be at that location.
        Otherwise, it is expected/build either next to the fastafile or in the
        workdir.

        :param fastafile: full path of fasta file to index
        :param indexname: name of index to create (optional)
        :returns valid indexname
        :raises  IndexBuildFailure if unable to create an index
        """

        if indexname is not None:
            # we have an index name, check if it exists or build a new one
            if not self.have_index(indexname, fastafile):
                WARNING("provided index not valid, trying to rebuild")
                self.build_index(fastafile, indexname)
        else:  # we have no indexname
            locations = [
                fastafile[:-6],  # fasta w/o ".fasta"
                os.path.join(self.workdir,
                             os.path.basename(fastafile))[:-6]
            ]

            # check for existing indices
            for name in locations:
                if self.have_index(name, fastafile):
                    return name

            # try building new index
            last_exception = None
            for name in locations:
                try:
                    self.build_index(fastafile, name)
                except IndexBuildFailure as e:
                    last_exception = e
                else:
                    indexname = name
                    break
            if last_exception is not None:
                raise last_exception
        return indexname

    @log.timed("Building Bowtie 2 index {indexname} from {fastafile}")
    def build_index(self, fastafile, indexname):
        """Build a index for fastafile named indexname.

        :param fastafile: fasta file to index
        :param indexname: name of index to create
        :raises IndexBuildFailure if index build command failed
        """
        if not os.access(fastafile, os.R_OK):
            raise IndexBuildFailure('Cannot read "{}"'.format(fastafile))
        cmd = [
            "bowtie2-build",
            "--offrate", "3",  # "SA is sampled every 2^<int> BWT chars"
            fastafile,
            indexname
        ]

        # TODO: this file will cause TempDir to remain. Should it be deleted?
        logfile = os.path.join(self.workdir,
                               os.path.basename(fastafile) + "_bt2_build.log")

        try:
            with File(logfile).writer() as out:
                check_call(cmd, stdout=out)
        except CalledProcessError as e:
            raise IndexBuildFailure(e.args)

        return indexname


class Bowtie(Mapper):
    binary = "bowtie"

    def __init__(self, *args, **kwargs):
        super(Bowtie, self).__init__(*args, **kwargs)
        self.binary = "bowtie"

    def prefilter_reads(self):
        raise NotImplementedError()

    def map_reads(self):
        raise NotImplementedError()


class Bwa(Mapper):
    binary = "bwa"

    def __init__(self, *args, **kwargs):
        super(Bwa, self).__init__(*args, **kwargs)
        self.binary = "bwa"

    def prefilter_reads(self):
        raise NotImplementedError()

    def map_reads(self):
        raise NotImplementedError()


class BBMap(Mapper):
    binary = "bbmap.sh"

    def __init__(self, *args, **kwargs):
        super(BBMap, self).__init__(*args, **kwargs)
        self.binary = "bbmap.sh"

    def prefilter_reads(self):
        raise NotImplementedError()

    def map_reads(self):
        raise NotImplementedError()


# module level helper functions for command line interface

_mappers = [Bowtie2, Bwa, BBMap]  # implemented mappers, in order of preference

_msg_no_mapper_found = """

No supported read mapper found in PATH!

  Please make sure that one of the read mappers supported by EMIRGE is installed
  on your system and that the directory it has been installed to is part of your
  PATH.

Supported read mappers: {}
""".format(", ".join([mapper.__name__ for mapper in _mappers]))


def get_options(parser):
    assert isinstance(parser, optparse.OptionParser)

    # check which mappers are available on system
    mappers_avail = ["'{}'".format(mpr.__name__) for mpr in _mappers
                     if mpr.available()]
    if len(mappers_avail) == 0:
        parser.error(_msg_no_mapper_found)

    # prepare help message for not available mappers
    mappers_supported = ["'{}'".format(mpr.__name__) for mpr in _mappers
                         if not mpr.available()]
    mappers_supported_msg = ""
    if len(mappers_supported) > 0:
        mappers_supported_msg = "; supported (but not installed): {}".format(
                ", ".join(mappers_supported)
        )

    opts = optparse.OptionGroup(parser, "Read-Mapping parameters")
    opts.add_option(
            "-m", "--mapper",
            choices=mappers_avail, default=mappers_avail[0],
            help="Choose read mapper to use. "
                 "(default: %default; available: {}{})"
                 .format(", ".join(mappers_avail), mappers_supported_msg)
    )
    opts.add_option(
            "-b", "--mapper-index",
            help="Alternative name or location for index over FASTA_DB. By "
                 "default, EMIRGE will try to find an appropriate index in the "
                 "same directory as the FASTA_DB."
    )
    opts.add_option(
            "-i", "--insert_mean",
            type="int", default=1500,
            help="Insert size distribution mean. (default: %default)")
    opts.add_option(
            "-s", "--insert_stddev",
            type="int", default=500,
            help="Insert size distribution standard deviation. (default: "
                 "%default")
    opts.add_option(
            "-P", "--prefilter-reads",
            action="store_true", default=False,
            help="Pre-filter reads using candidate db (FASTA_DB). Only reads "
                 "(loosely) matching the candidate db will be considered by "
                 "EMIRGE in all iterations. Recommended for metagenome data."
    )
    parser.add_option_group(opts)

    # TODO
    # -m, --mapping: path to pre-computed initial mapping file


def get_mapper(opts, workdir, candidates, fwd_reads, rev_reads, threads,
               phred33, ):
    mapperclass = None
    for mpr in _mappers:
        if mpr.__name__ == opts.mapper.strip("'"):
            mapperclass = mpr
            break
    if mapperclass is None:
        raise Exception("Could not find selected mapper?! (mapper={})"
                        .format(opts.mapper))
    assert issubclass(mapperclass, Mapper)
    INFO("Using read mapper '{}'".format(mapperclass.__name__))
    mpr = mapperclass(candidates=candidates,
                      candidates_index=opts.mapper_index,
                      workdir=workdir,
                      fwd_reads=fwd_reads,
                      rev_reads=rev_reads,
                      threads=threads,
                      phred33=phred33,
                      insert_mean=opts.insert_mean,
                      insert_sd=opts.insert_stddev,
                      prefilter_reads=opts.prefilter_reads)
    return mpr
