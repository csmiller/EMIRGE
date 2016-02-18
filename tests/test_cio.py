from tempfile import TemporaryFile, NamedTemporaryFile
from timeit import Timer

from nose.tools import assert_equal
from numpy.lib.function_base import median

from Emirge.cio import enumerate_fastq, fastq_detect_len_and_encoding
from Emirge.io import File, check_call
from Emirge.log import INFO

read_file_1 = "tests/test_data/ten_seq_community_000_50K_L150_I350.1.fastq"

def test_enumerate_fastq():
    out = NamedTemporaryFile()

    with File(read_file_1) as f:
        reads, read_len, phred33 = enumerate_fastq(f, out)

    assert_equal(reads, 50000)
    assert_equal(read_len, 150)
    assert_equal(phred33, True)


def test_fastq_detect_len_and_encoding():
    with File(read_file_1) as f:
        reads, read_len, phred33 = fastq_detect_len_and_encoding(f)
    assert_equal(reads, 1000)
    assert_equal(read_len, 150)
    assert_equal(phred33, True)

    with File(read_file_1) as f:
        reads, read_len, phred33 = fastq_detect_len_and_encoding(f,num=0)
    assert_equal(reads, 50000)
    assert_equal(read_len, 150)
    assert_equal(phred33, True)


def notest_enumerate_fastq_timing():
    tm = Timer(stmt="test_enumerate_fastq()",
               setup="from tests.test_cio import test_enumerate_fastq")
    times = tm.repeat(number=1, repeat=10)
    INFO(median(times))
