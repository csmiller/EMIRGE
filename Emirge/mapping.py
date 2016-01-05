"""implement mappers"""

import pysam

from Emirge.io import EnumerateReads, filter_fastq, File, make_pipe
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
            fwd_reads = File(fwd_reads)
        if isinstance(rev_reads, str):
            rev_reads = File(rev_reads)
        # if isinstance(candidates, str):
        #     candidates = File(candidates)

        self.candidates = candidates
        self.phred33 = phred33
        self.threads = threads

        if reindex:
            self.fwd_reads = EnumerateReads(fwd_reads)
            if rev_reads:
                self.rev_reads = EnumerateReads(rev_reads)
        else:
            self.fwd_reads = open(fwd_reads, "r")
            if rev_reads:
                self.rev_reads = open(rev_reads, "r")


class Bowtie2(Mapper):
    def __init__(self, *args, **kwargs):
        super(Bowtie2, self).__init__(*args, **kwargs)
        self.binary = "bowtie2"

    def make_basecmd(self):
        cmd = [
            self.binary,
            "-p", str(self.threads)
        ]
        if self.rev_reads:
            cmd += [
                "-1", self.fwd_reads,
                "-2", self.rev_reads,
                "-x", self.candidates,
            ]
        else:
            cmd += [
                "-U", self.fwd_reads
            ]
        if self.phred33:
            cmd += ["--phred33"]
        else:
            cmd += ["--phred64"]

        return cmd

    def prefilter_reads(self):
        """Pre-filters read file to include only reads aligning to seed data"""
        cmd = self.make_basecmd()
        cmd += [
            "--very-sensitive-local",
            "-k1",  # report up to 1 alignments per read
            "--no-unal",  # suppress SAM records for unaligned reads
        ]
        Prefilter = make_pipe("prefilter", cmd)

        keepers = set()
        i = 0
        with Prefilter(None) as p:
            with pysam.AlignmentFile(p, "r") as sam:
                for read in sam.fetch(until_eof=True):
                    if not read.is_unmapped:
                        keepers.add(read.query_name)
                    i+=1

        INFO("Found %i sequences" % len(keepers))

        with self.fwd_reads.reader() as reads:
            self.fwd_reads = filter_fastq(reads, keepers)
        with self.rev_reads.reader() as reads:
            self.rev_reads = filter_fastq(reads, keepers)


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
