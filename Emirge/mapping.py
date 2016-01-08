"""Read mapping module

This module abstracts from the various read mapping tools.
"""
import os
from subprocess import CalledProcessError

import pysam

import Emirge.log as log
from Emirge.io import EnumerateReads, filter_fastq, decompressed, make_pipe, \
    check_call, TempDir, File
from Emirge.log import INFO, DEBUG, WARNING


class IndexBuildFailure(Exception):
    """
    Exception raised when something went wrong trying to build
    the mapping index
    """
    pass


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
    """

    def __init__(self, candidates, fwd_reads, rev_reads=None, phred33=False,
                 threads=8, reindex=True, workdir=None):
        """Create new mapper.

        Args:
            candidates -- fasta file containing reference sequences
            fwd_reads  -- forward reads file
            rev_reads  -- reverse reads file, may be none for single ended
            phred33    -- true if quality encoded as phred33, false for 64
            threads    -- number of threads mappers should use
            reindex   -- change read names to consecutive numbers
            workdir    -- (optional) designed work directory
        """
        super(Mapper, self).__init__()

        if isinstance(fwd_reads, str):
            fwd_reads = decompressed(fwd_reads)
        if isinstance(rev_reads, str):
            rev_reads = decompressed(rev_reads)

        self.candidates = candidates
        self.phred33 = phred33
        self.threads = threads

        if reindex:
            fwd_reads = EnumerateReads(fwd_reads)
            if rev_reads:
                rev_reads = EnumerateReads(rev_reads)

        self.fwd_reads, self.rev_reads = fwd_reads, rev_reads

        if workdir is None:
            workdir = TempDir()
        self.workdir = workdir


class Bowtie2(Mapper):
    """Handle read mapping with Bowtie2."""

    def __init__(self, *args, **kwargs):
        super(Bowtie2, self).__init__(*args, **kwargs)
        self.binary = "bowtie2"

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
    def prefilter_reads(self):
        """Pre-filter reads using candidates.

        Switches the fwd/rev_reads property to a temporary file holding only
        reads mapping (loosely) to the candidates file
        """

        # prepare bowtie2 command
        cmd = self.make_basecmd()
        cmd += [
            "-x", self.prep_index(self.candidates),  # index
            "--very-sensitive-local",
            "-k", "1",  # report up to 1 alignments per read
            "--no-unal",  # suppress SAM records for unaligned reads
        ]

        # run bowtie2 prefiltering command, use pysam to parse matching
        # read names from output into set 'keepers'
        keepers = set()
        i = 0
        with make_pipe("prefilter", cmd)(None) as p:
            with pysam.AlignmentFile(p, "r") as sam:
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
        with self.fwd_reads.reader() as reads:
            self.fwd_reads, fwd_count, fwd_matches = \
                filter_fastq(reads, keepers)

        # if we have a reverse read file, do the same for that
        if self.rev_reads is not None:
            with self.rev_reads.reader() as reads:
                self.rev_reads, rev_count, rev_matches = \
                    filter_fastq(reads, keepers)

            # number of reads and must be identical for forward and reverse
            assert fwd_count == rev_count
            assert fwd_matches == rev_matches

        INFO("Pre-Mapping: {} out of {} reads remaining"
             .format(fwd_matches, fwd_count))

    @log.timed('Mapping reads to "{fastafile}"')
    def map_reads(self, fastafile, workdir, insert_mean=1500, insert_sd=500):
        """Create bam file containing mappings for fastafile in workdir.

        :param fastafile:   target sequences
        :param workdir:     directory to hold bamfile
        :param insert_mean: mean insert size
        :param insert_sd:   stddev of insert size distribution
        """
        minins = max((insert_mean - 3 * insert_sd), 0)
        maxins = insert_mean + 3 * insert_sd

        log_file = os.path.join(workdir, "bt2.log")

        cmd = self.make_basecmd()
        cmd += [
            "-x", self.prep_index(fastafile),  # index
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

        assert False, "TODO: complete function"

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
            file_name = indexname + suffix
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
                self.build_index(fastafile, indexname)
        else:  # we have no indexname
            locations = [
                fastafile[:-6],  # fasta w/o ".fasta"
                os.path.join(self.workdir.name,
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

    @log.timed("Building Bowtie 2 index {fastafile} from {indexname}")
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
        logfile = os.path.join(self.workdir.name,
                               os.path.basename(fastafile) + "_bt2_build.log")

        try:
            with File(logfile).writer() as out:
                check_call(cmd, stdout=out)
        except CalledProcessError as e:
            raise IndexBuildFailure(e.args)

        return indexname


class Bowtie(Mapper):
    def __init__(self, *args, **kwargs):
        super(Bowtie, self).__init__(*args, **kwargs)
        self.binary = "bowtie"


class Bwa(Mapper):
    def __init__(self, *args, **kwargs):
        super(Bwa, self).__init__(*args, **kwargs)
        self.binary = "bwa"


class BBMap(Mapper):
    def __init__(self, *args, **kwargs):
        super(BBMap, self).__init__(*args, **kwargs)
        self.binary = "bbmap.sh"
