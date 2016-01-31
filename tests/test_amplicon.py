"""test functions in Emirge/amplicon.pyx"""
import numpy

from Emirge import amplicon, io
from nose.tools import assert_equal

from Emirge.log import INFO

read_file_1 = "tests/test_data/ten_seq_community_000_50K_L150_I350.2.fastq"


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


class EM_test():
    @classmethod
    def setup_class(cls):
        cls.reads1_filepath = io.File(read_file_1)
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