"""Tests functions in Emirge.io"""

from Emirge import io


def test_Record_empty():
    record = io.Record()
    assert record.title == ""
    assert record.sequence == ""
    assert str(record) == ">\n\n"


def test_Record_formatting():
    record = io.Record(
        title="id123",
        sequence="GUGCAAAGUUGUGUAGUGCGAUCGGUGGAUGCCUUGGCACCA"
                 "AGAGCCGAUGAAGGACGUUGUGACCUGCGAUAAGCCCU"
                 "GGGGAGUUGGUGAGCGAGCUGUGAUCCGGGGGUGUCCGAAUG"
                 "GGGAAACCUGGAAUGUCCGGAGUAGUGUCCGGUGGCCC"
                 "UGCCCUGAAUGUAUAGGGGUGUGGGUGGUAACGCGGGGAAGU"
                 "GAAACAUCUUAGUACCCGUAGGAAGAGAAAACAAGUGU"
    )
    assert str(record) == \
        """>id123
GUGCAAAGUUGUGUAGUGCGAUCGGUGGAUGCCUUGGCACCAAGAGCCGAUGAAGGACGU
UGUGACCUGCGAUAAGCCCUGGGGAGUUGGUGAGCGAGCUGUGAUCCGGGGGUGUCCGAA
UGGGGAAACCUGGAAUGUCCGGAGUAGUGUCCGGUGGCCCUGCCCUGAAUGUAUAGGGGU
GUGGGUGGUAACGCGGGGAAGUGAAACAUCUUAGUACCCGUAGGAAGAGAAAACAAGUGU
"""


def test_FastIterator():
    pass


def test_ReindexReads():
    pass
    # read_file = "test_data/ten_seq_community_000_50K_L150_I350.2.fastq"
    # new_read_file = ReindexReads(read_file)
