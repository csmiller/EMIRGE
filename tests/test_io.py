"""Tests functions in Emirge.io"""

from StringIO import StringIO
from subprocess import Popen, PIPE
from tempfile import NamedTemporaryFile

from Emirge import io

# === FASTA test data ===

# sequence formatted at 60 cols
fasta_sample_60 = """>id123
GUGCAAAGUUGUGUAGUGCGAUCGGUGGAUGCCUUGGCACCAAGAGCCGAUGAAGGACGU
UGUGACCUGCGAUAAGCCCUGGGGAGUUGGUGAGCGAGCUGUGAUCCGGGGGUGUCCGAA
UGGGGAAACCUGGAAUGUCCGGAGUAGUGUCCGGUGGCCCUGCCCUGAAUGUAUAGGGGU
GUGGGUGGUAACGCGGGGAAGUGAAACAUCUUAGUACCCGUAGGAAGAGAAAACAAGUGU
"""

# same sequence formated at 77 cols
fasta_sample_77 = """>id123
GUGCAAAGUUGUGUAGUGCGAUCGGUGGAUGCCUUGGCACCAAGAGCCGAUGAAGGACGUUGUGACCUGCGAUAAGC
CCUGGGGAGUUGGUGAGCGAGCUGUGAUCCGGGGGUGUCCGAAUGGGGAAACCUGGAAUGUCCGGAGUAGUGUCCGG
UGGCCCUGCCCUGAAUGUAUAGGGGUGUGGGUGGUAACGCGGGGAAGUGAAACAUCUUAGUACCCGUAGGAAGAGAA
AACAAGUGU
"""

# same sequence, one line
sequence_sample = (
    "GUGCAAAGUUGUGUAGUGCGAUCGGUGGAUGCCUUGGCACCAAGAGCCGAUGAAGGACGU"
    "UGUGACCUGCGAUAAGCCCUGGGGAGUUGGUGAGCGAGCUGUGAUCCGGGGGUGUCCGAA"
    "UGGGGAAACCUGGAAUGUCCGGAGUAGUGUCCGGUGGCCCUGCCCUGAAUGUAUAGGGGU"
    "GUGGGUGGUAACGCGGGGAAGUGAAACAUCUUAGUACCCGUAGGAAGAGAAAACAAGUGU"
)

read_file_1 = "tests/test_data/ten_seq_community_000_50K_L150_I350.2.fastq"


# === test functions ===

def test_Record_empty():
    record = io.Record()
    assert record.title == ""
    assert record.sequence == ""
    assert str(record) == ">\n\n"


def test_Record_formatting():
    record = io.Record(title="id123", sequence=sequence_sample)
    print sequence_sample
    assert str(record) == fasta_sample_60


def test_FastIterator():
    n = 10
    fasta_file = StringIO(fasta_sample_77 * n)
    i = 0
    for record in io.FastIterator(fasta_file):
        assert str(record) == fasta_sample_60
        i = i+1
    assert i == n


def cmp_reindexed_fq_files(orig, reindexed, nseq):
    orig.seek(0)
    lineno = 0
    for sline, dline in zip(orig, reindexed):
        if lineno % 4 == 0:
            assert dline.rstrip() == "@" + str(lineno/4)
        else:
            assert dline == sline
        lineno += 1
    assert nseq == lineno / 4


def test_ReindexReads():
    nlines = 4 * 1200

    # prep input files
    src = NamedTemporaryFile()
    src_zipped = NamedTemporaryFile(suffix=".gz")
    zipper = Popen(["gzip", "-c"], stdin=PIPE, stdout=src_zipped).stdin

    with open(read_file_1) as f:
        for line, n in zip(f, range(0, nlines)):
            src.write(line)
            zipper.write(line)

    src.flush()
    zipper.close()
    src_zipped.flush()

    dst, n_reads = io.ReindexReads(src.name)
    assert n_reads == nlines / 4
    cmp_reindexed_fq_files(src, dst, n_reads)

    dst, n_reads = io.ReindexReads(src_zipped.name)
    assert n_reads == nlines / 4
    cmp_reindexed_fq_files(src, dst, n_reads)
