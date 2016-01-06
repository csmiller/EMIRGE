"""implement mappers"""

import pysam

from Emirge.io import EnumerateReads, filter_fastq, decompressed, make_pipe
from Emirge.log import INFO, DEBUG


class Mapper(object):
    def __init__(self,
                 candidates,
                 fwd_reads,
                 rev_reads=None,
                 phred33=False,
                 threads=8,
                 reindex=True):
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


class Bowtie2(Mapper):
    def __init__(self, *args, **kwargs):
        super(Bowtie2, self).__init__(*args, **kwargs)
        self.binary = "bowtie2"

    def make_basecmd(self):
        cmd = [
            self.binary,
            "-p", str(self.threads),
            "-x", self.candidates,
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

    def prefilter_reads(self):
        """Pre-filters read file to include only reads aligning to seed data"""

        # prepare bowtie2 command
        cmd = self.make_basecmd()
        cmd += [
            "--very-sensitive-local",
            "-k1",  # report up to 1 alignments per read
            "--no-unal",  # suppress SAM records for unaligned reads
        ]
        Prefilter = make_pipe("prefilter", cmd)

        # run bowtie2 prefiltering command, use pysam to parse matching
        # read names from output into set 'keepers'
        keepers = set()
        i = 0
        with Prefilter(None) as p:
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

    def map_reads(self):
        pass


class Bowtie(Mapper):
    def __init__(self, *args, **kwargs):
        super(Bowtie, self).__init__(*args, **kwargs)
        self.binary = "bowtie"


class Bwa(Mapper):
    def __init__(self, *args, **kwargs):
        super(Bwa, self).__init__(*args, **kwargs)
        self.binary = "bwa"


class Bbmap(Mapper):
    def __init__(self, *args, **kwargs):
        super(Bbmap, self).__init__(*args, **kwargs)
        self.binary = "bbmap.sh"
