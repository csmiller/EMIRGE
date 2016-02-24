from nose.tools import assert_equal

from Emirge.clustering import vsearch


def test_vsearch_available():
    assert_equal(vsearch.available(), True)