"""test functions in Emirge/amplicon.pyx"""

from Emirge import amplicon
from nose.tools import assert_equal


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
