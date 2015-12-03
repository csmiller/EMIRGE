"""Implements IO functions"""

import re
from tempfile import NamedTemporaryFile
from subprocess import CalledProcessError, check_call, Popen, PIPE

from Emirge import log


class Record:
    """
    Stripped down FASTA record class with same members as
    Biopython FASTA Record Class:

        title    -- with > removed
        sequence -- as string
                    assumed to be free of whitespace and other non-sequence
                    characters
    """

    _colwidth = 60  # number of residues per FASTA line

    def __init__(self, title="", sequence=""):
        """
        Create a new Record.
        """
        self.title = title
        self.sequence = sequence

    def __str__(self):
        return ">%s\n%s\n" % (
            self.title,
            "\n".join(re.findall(".{1,%s}" % self._colwidth, self.sequence))
        )


def FastIterator(filehandle, dummyParser=None, record=None):
    """
    Generator returning Records from FASTA. Maybe 160% faster on test case.
    MAY RAISE MemoryError ON LARGE FASTA FILES
    IN:  file object
         dummyParser is a placeholder for RecordParser from Biopython.  Unused.
         (optional) record to use as template

    NOTE: this iterator is fast, but breaks easily with nonstandard input,
    e.g. if there are "\r" in endlines.
    """
    if record is None:
        record = Record()

    for recordstring in re.split('\n>', filehandle.read()[1:]):
        record.title, record.sequence = recordstring.split('\n', 1)
        record.sequence = record.sequence.replace('\n', '').replace(' ', '')
        yield record


@log.timed("Rewriting reads with indices in headers")
def ReindexReads(reads_filepath):
    """
    Replaces sequence headers ("@...") with the sequence number.
    Although this requires an inefficient rewrite of the fastq file,
    it means that reading of bam files does not require a costly separate
    id lookup step on the read name.

    Returns tuple:
        new_reads_file  NamedTemporaryFile (will self-delete) of rewritten fq
        num_reads       int number of sequences in fq
    """

    tmp_n_reads_file_path = NamedTemporaryFile()
    new_reads_filepath = NamedTemporaryFile(suffix="reindexed_reads.fq")

    try:
        if reads_filepath.endswith('.gz'):
            src = Popen(["gunzip", "-c", reads_filepath],
                        stdout=PIPE).stdout
        else:
            src = open(reads_filepath)

        check_call(['awk',
                    '{ if ((NR-1) %% 4 == 0) { print "@"(NR-1)/4 }'
                    'else { print } }'
                    'END { print (NR)/4 > "%s" }'
                    % (tmp_n_reads_file_path.name)],
                   stdin=src, stdout=new_reads_filepath)

        n_reads = int(tmp_n_reads_file_path.readline().strip())
    except CalledProcessError:
        log.error("awk rewrite of reads failed! Is awk installed?")
        raise

    new_reads_filepath.seek(0)
    return (new_reads_filepath, n_reads)

@log.timed("Counting reads in input files")
def FastqCountReads(filename):
    cmd = "cat %s | wc -l" % (filename)
    if filename.endswith('.gz'):
        cmd = "z" + cmd
    p = Popen(cmd, shell=True, stdout=PIPE)
    stdoutdata, stderrdata = p.communicate()
    return int(stdoutdata.strip()) / 4