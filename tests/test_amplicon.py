from Emirge import amplicon
from nose.tools import assert_equal

def test_complement_numeric_base():
    """
    A:0, T:1, C:2, G:3, N:4
    => 0<->1, 2<->3, 4<->4
    """

    for complement, base in enumerate([1,0,3,2,4]):
       assert_equal(amplicon._complement_numeric_base(base), complement)
