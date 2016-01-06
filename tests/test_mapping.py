"""Test mapping code"""

from nose.tools import assert_equal

from Emirge import mapping, io, log

log.setup(debug=True)

read_file_1 = "tests/test_data/ten_seq_community_000_50K_L150_I350.1.fastq.xz"
read_file_2 = "tests/test_data/ten_seq_community_000_50K_L150_I350.2.fastq.xz"
cand_file = "tests/test_data/twenty_seq_database"


def test_bt2_prefilter():
    bt2 = mapping.Bowtie2(cand_file, read_file_1, read_file_2,
                          phred33=True, reindex=True)
    bt2.prefilter_reads()
    assert_equal(io.fastq_count_reads(bt2.fwd_reads), 50000)
    assert_equal(io.fastq_count_reads(bt2.rev_reads), 50000)


def test_bt2_prefilter_noreindex():
    bt2 = mapping.Bowtie2(cand_file, read_file_1, read_file_2,
                          phred33=True, reindex=False)
    bt2.prefilter_reads()
    assert_equal(io.fastq_count_reads(bt2.fwd_reads), 50000)
    assert_equal(io.fastq_count_reads(bt2.rev_reads), 50000)


def test_bt2_prefilter_onlyfwd():
    bt2 = mapping.Bowtie2(cand_file, read_file_1,
                          phred33=True, reindex=True)
    bt2.prefilter_reads()
    assert_equal(io.fastq_count_reads(bt2.fwd_reads), 49993)


def test_bt2_prefilter_onlyfwd_noreindex():
    bt2 = mapping.Bowtie2(cand_file, read_file_1,
                          phred33=True, reindex=False)
    bt2.prefilter_reads()
    assert_equal(io.fastq_count_reads(bt2.fwd_reads), 49993)
