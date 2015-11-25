"""Implements IO functions"""

import re

class Record:
    """
    stripped down FASTA record class with same members as
    Biopython FASTA Record Class
        title    -- with > removed
        sequence -- as string
    """
    _colwidth = 60 # colwidth specifies the number of residues to put on each line when generating FASTA format.
    def __init__(self, title = "", sequence = ""):
        """
        Create a new Record.

        """
        self.title = title
        self.sequence = sequence
    def __str__(self):
        return ">%s\n%s\n"%(self.title, "\n".join(re.findall(".{1,%s}"%self._colwidth, self.sequence)))

def FastIterator(filehandle, dummyParser = None, record = None):
    """
    a generator to return records one at a time.  Maybe 160% faster on test case.
    MAY RAISE MemoryError ON LARGE FASTA FILES
    IN:  file object
         dummyParser is a placeholder for RecordParser from Biopython.  Not used.
         a record to use for yielding.  Otherwise create an empty one with standard init

    NOTE: this fasta iterator is fast, but it breaks easily with nonstandard input, for example
    if there are "\r" in endlines.
    """
    if record is None:
        record = Record()

    for recordstring in re.split('\n>', filehandle.read()[1:]):
        record.title, record.sequence = recordstring.split('\n',1)
        record.sequence = record.sequence.replace('\n','').replace(' ','')
        yield record
