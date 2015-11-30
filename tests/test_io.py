"""Tests functions in Emirge.io"""

from StringIO import StringIO

from Emirge import io

## FASTA test data

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


def test_ReindexReads():
    pass
    # read_file = "test_data/ten_seq_community_000_50K_L150_I350.2.fastq"
    # new_read_file = ReindexReads(read_file)
