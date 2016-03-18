"""test functions in Emirge/amplicon.pyx"""
from timeit import Timer

import numpy
from numpy.lib.function_base import median
from numpy.testing.utils import assert_array_equal

from Emirge import amplicon, io
from nose.tools import assert_equal

from Emirge.log import INFO

read_file_1 = "tests/test_data/ten_seq_community_000_50K_L150_I350.2.fastq.xz"


sequence_sample = (
    "GTGCAAAGTTGTGTAGTGCGATCGGTGGATGCCTTGGCACCAAGAGCCGATGAAGGACGT"
    "TGTGACCTGCGATAAGCCCTGGGGAGTTGGTGAGCGAGCTGTGATCCGGGGGTGTCCGAA"
    "TGGGGAAACCTGGAATGTCCGGAGTAGTGTCCGGTGGCCCTGCCCTGAATGTATAGGGGT"
    "GTGGGTGGTAACGCGGGGAAGTGAAACATCTTAGTACCCGTAGGAAGAGAAAACAAGTGT"
)


def test_base_alpha2int():
    base2i = {"A": 0, "T": 1, "C": 2, "G": 3, "N": 4}
    for base in base2i:
        assert_equal(base2i[base], amplicon.base_alpha2int(ord(base)))


def test_seq_alpha2int():
    bases = "ATCGN"
    ibases = [0, 1, 2, 3, 4]
    itestbases = amplicon.seq_alpha2int(bases, 5)
    assert_equal(len(ibases), len(itestbases))
    for ibase, itestbase in zip(ibases, itestbases):
        assert_equal(ibase, itestbase)


def test_complement_numeric_base():
    # A:0, T:1, C:2, G:3, N:4 => 0<->1, 2<->3, 4<->4
    for complement, base in enumerate([1, 0, 3, 2, 4]):
        assert_equal(amplicon.complement_numeric_base(base), complement)


class Complement_test():
    @classmethod
    def setup_class(cls):
        cls.sequence = [amplicon.base_alpha2int(ord(x))
                        for x in sequence_sample]
        cls.sequence_comp = [amplicon.complement_numeric_base(x)
                             for x in cls.sequence]
        cls.nd_sequence = numpy.array(cls.sequence, dtype=numpy.uint8)
        cls.nd_sequence_comp = numpy.array(cls.sequence_comp, dtype=numpy.uint8)
        cls.result_buf = numpy.zeros_like(cls.nd_sequence)

    def test_complement_sequence(self):
        for start in 0, 20, 150:
            for stop in 0, 1, 22, 149, 150:
                if start+stop > 150 or stop < start:
                    continue
                amplicon.complement_sequence(
                    self.nd_sequence[start:stop],
                    self.result_buf
                )
                assert_array_equal(
                    self.nd_sequence_comp[start:stop],
                    self.result_buf[0:stop-start],
                    "start={} stop={}".format(start,stop)
                )



class EM_test():
    @classmethod
    def setup_class(cls):
        cls.reads1_filepath = io.decompressed(read_file_1)
        cls.reads2_filepath = None
        cls.n_reads = 50000
        cls.max_read_length = 150
        cls.reads = numpy.zeros((cls.n_reads, 2, cls.max_read_length),
                                dtype=numpy.uint8)
        cls.quals = numpy.zeros_like(cls.reads)
        cls.readlengths = numpy.zeros((cls.n_reads, 2),
                                      dtype=numpy.uint16)
        cls.reads_ascii_offset = 33

    def test_popuplate_reads_arrays(self):
        amplicon.populate_reads_arrays(self)
        total_bases = self.readlengths.sum(axis=0)[0]
        nuc_counts = [(self.reads[:,0,:] == b).sum() for b in range(4)]

        assert_equal(total_bases, self.n_reads * self.max_read_length)
        assert_equal(nuc_counts, [1679866, 1665236, 2052473, 2102425])
        assert_equal(self.quals.sum(), 271176838)

    def complement_reads(self):
        x = numpy.zeros_like(self.reads[0,0])
        for read in range(10000):
            amplicon.complement_sequence(self.reads[read,0], x)

    def complement_reads_mv(self):
        x = numpy.zeros_like(self.reads[0,0])
        for read in range(10000):
            amplicon.complement_sequence_mv(self.reads[read,0], x)

    def test_complement_sequence_speed(self):
        tm = Timer(stmt="x.complement_reads()",
                   setup="from tests.test_amplicon import EM_test;"
                         "x = EM_test(); ",
                   )
        x = tm.repeat(number=1, repeat=20)
        INFO(median(x) * 10)

    def test_complement_sequence_speed_mv(self):
        tm = Timer(stmt="x.complement_reads_mv()",
                   setup="from tests.test_amplicon import EM_test;"
                         "x = EM_test(); ",
                   )
        x = tm.repeat(number=1, repeat=20)
        INFO(median(x) * 10)


    def test_dumm2(self):
        pass
#\
# def test_complement():
#     max_readlen = 150
#     n_seqs = 100
#     x = numpy.ndarray([n_seqs, max_readlen], dtype=numpy.uint8)
#     with io.Kseq(read_file_1) as k:
#         for read in range(n_seqs):
#             k.read_next_sequence()
#             for base in range()
#             x[read, k.seq
