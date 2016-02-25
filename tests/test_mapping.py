"""Test mapping code"""

import os

from nose.tools import assert_equal
from shutil import copy

from Emirge import mapping, io, log
from Emirge.io import TempDir

log.setup(debug=True)

read_file_1 = "tests/test_data/ten_seq_community_000_50K_L150_I350.1.fastq.xz"
read_file_2 = "tests/test_data/ten_seq_community_000_50K_L150_I350.2.fastq.xz"
cand_file = "tests/test_data/twenty_seq_database.fasta"


def test_bt2_prefilter():
    bt2 = mapping.Bowtie2(cand_file, read_file_1, read_file_2, phred33=True)
    bt2.prefilter_reads()
    assert_equal(io.fastq_count_reads(bt2.fwd_reads), 50000)
    assert_equal(io.fastq_count_reads(bt2.rev_reads), 50000)


def test_bt2_prefilter_onlyfwd():
    bt2 = mapping.Bowtie2(cand_file, read_file_1, phred33=True)
    bt2.prefilter_reads()
    assert_equal(io.fastq_count_reads(bt2.fwd_reads), 49993)


def test_bt2_build_index():
    tmpdir = TempDir()
    tmpfasta = os.path.join(tmpdir.name, "test.fasta")
    copy(cand_file, tmpfasta)
    bt2 = mapping.Bowtie2(cand_file, read_file_1, phred33=True)
    bt2.prep_index(tmpfasta)


def test_bt2_map_reads():
    tmpdir = TempDir()
    bt2 = mapping.Bowtie2(cand_file, read_file_1, read_file_2, phred33=True)
    frags_mapped, bam_file = bt2.map_reads(cand_file, tmpdir.name)
    assert_equal(frags_mapped, 49435)


if __name__ == '__main__':
    test_bt2_prefilter()
