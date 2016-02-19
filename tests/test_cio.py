from tempfile import NamedTemporaryFile
from timeit import Timer

from nose.tools import assert_equal
from numpy.lib.function_base import median

from Emirge.cio import fastq_process
from Emirge.io import File
from Emirge.log import INFO

read_file_1 = "tests/test_data/ten_seq_community_000_50K_L150_I350.1.fastq"


keepers = {
    "JQ186203.1.1363_864_1173_0:0:0_0:0:0_0",  # first sequence
    "FJ906931.1.1234_474_831_0:0:0_0:0:0_119b",  # last sequence
    "JQ186203.1.1363_361_721_0:0:0_0:0:0_4",
    "JQ186203.1.1363_82_418_0:0:0_0:0:0_3d7",
    "JQ186203.1.1363_171_531_0:0:0_0:0:0_121",
    "JQ186203.1.1363_788_1096_0:0:0_0:0:0_133",
    "JQ186203.1.1363_520_833_0:0:0_0:0:0_10",
    "AY154489.1.1330_293_612_0:0:0_0:0:0_1d5",
}


class TestFastqProcess(object):
    def __init__(self):
        self.out = None

    def setup(self):  # run before each test
        self.out = NamedTemporaryFile()

    def teardown(self):  # run after each test
        self.out = None

    def test_enumerate(self):
        with File(read_file_1) as f:
            # reads, read_len, phred33 = enumerate_fastq(f, out)
            reads, read_len, phred33 = fastq_process(f, self.out)

        assert_equal(reads, 50000)
        assert_equal(read_len, 150)
        assert_equal(phred33, True)

    def test_filter(self):
        with File(read_file_1) as f:
            reads, read_len, phred33 = fastq_process(f, self.out, keepers)

        assert_equal(reads, len(keepers))
        assert_equal(read_len, 150)
        assert_equal(phred33, True)

    @staticmethod
    def test_detect_len_and_encoding():
        with File(read_file_1) as f:
            reads, read_len, phred33 = fastq_process(f)
        assert_equal(reads, 50000)
        assert_equal(read_len, 150)
        assert_equal(phred33, True)

    @staticmethod
    def test_detect_len_and_encoding_fast():
        with File(read_file_1) as f:
            reads, read_len, phred33 = fastq_process(f, num=1000)
        assert_equal(reads, 1000)
        assert_equal(read_len, 150)
        assert_equal(phred33, True)


def test_timing_enumerate():
    tm = Timer(stmt="T.test_enumerate();",
               setup="from tests.test_cio import TestFastqProcess;"
                     "T=TestFastqProcess(); T.setup()")
    times = tm.repeat(number=1, repeat=10)
    INFO(median(times))


def test_timing_filter():
    tm = Timer(stmt="T.test_filter();",
               setup="from tests.test_cio import TestFastqProcess;"
                     "T=TestFastqProcess(); T.setup()")
    times = tm.repeat(number=1, repeat=10)
    INFO(median(times))
