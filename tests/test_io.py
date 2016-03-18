"""Tests functions in Emirge.io"""

import os
import re
from StringIO import StringIO
from subprocess import Popen, PIPE
from tempfile import NamedTemporaryFile

from nose.tools import assert_true, assert_false, assert_equal, assert_raises

from Emirge import io
# === FASTA test data ===

# sequence formatted at 60 cols
from Emirge.io import decompressed
from Emirge.log import WARNING

fasta_sample_60 = """>id123
GUGCAAAGUUGUGUAGUGCGAUCGGUGGAUGCCUUGGCACCAAGAGCCGAUGAAGGACGU
UGUGACCUGCGAUAAGCCCUGGGGAGUUGGUGAGCGAGCUGUGAUCCGGGGGUGUCCGAA
UGGGGAAACCUGGAAUGUCCGGAGUAGUGUCCGGUGGCCCUGCCCUGAAUGUAUAGGGGU
GUGGGUGGUAACGCGGGGAAGUGAAACAUCUUAGUACCCGUAGGAAGAGAAAACAAGUGU
"""

# same sequence formatted at 77 cols
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

read_file_1 = decompressed("tests/test_data" \
                          "/ten_seq_community_000_50K_L150_I350.2.fastq.xz")

bam_file = "tests/test_data/test.bam"
sam_file = "tests/test_data/test.sam"


# === helper functions ===


def assert_re_match(regex, string):
    try:
        assert re.match(regex, string) is not None
    except AssertionError as e:
        # suppress PyCharm warning (yes, AssertionError has setter for args!)
        # noinspection PyPropertyAccess
        e.args += ('"{}" does not match regex "{}"'.format(string, regex),)
        raise


# === test functions ===

def test_Record_empty():
    record = io.Record()
    assert record.title == ""
    assert record.sequence == ""
    assert str(record) == ">\n\n"


def test_Record_formatting():
    record = io.Record(title="id123", sequence=sequence_sample)
    assert str(record) == fasta_sample_60


def test_FastIterator():
    n = 10
    fasta_file = StringIO(fasta_sample_77 * n)
    i = 0
    for record in io.FastIterator(fasta_file):
        assert str(record) == fasta_sample_60
        i += 1
    assert i == n


def cmp_reindexed_fq_files(orig, reindexed, nseq):
    # orig.seek(0)
    line_no = 0
    for src_line, dst_line in zip(orig, reindexed):
        if line_no % 4 == 0:
            assert_equal(dst_line.rstrip(), "@" + str(line_no / 4))
        else:
            assert_equal(dst_line, src_line)
        line_no += 1
    assert_equal(nseq, line_no / 4)


def test_reindex_reads():
    num_lines = 4 * 1200

    # prep input files
    src = NamedTemporaryFile()
    with read_file_1 as f:
        for line, n in zip(f, range(0, num_lines)):
            src.write(line)

    src.flush()

    dst, n_reads = io.reindex_reads(src.name)
    assert_equal(n_reads, num_lines / 4)
    src.seek(0)
    cmp_reindexed_fq_files(src, dst, n_reads)


def test_reindex_reads_zipped():
    num_lines = 4 * 1200
    # prep input files
    src = NamedTemporaryFile(suffix=".gz")
    zipper = Popen(["gzip", "-c"], stdin=PIPE, stdout=src)

    with read_file_1 as f:
        for line, n in zip(f, range(0, num_lines)):
            zipper.stdin.write(line)

    zipper.stdin.close()
    zipper.wait()
    src.flush()

    dst, n_reads = io.reindex_reads(src.name)

    assert_equal(n_reads, num_lines / 4)
    with read_file_1 as f:
        cmp_reindexed_fq_files(f, dst, n_reads)


# noinspection PyPep8Naming
def test_NamedPipe():
    pipe = io.NamedPipe()
    pipe_file = pipe.name
    assert_true(io.ispipe(pipe_file))
    del pipe
    assert_false(io.ispipe(pipe_file))


# noinspection PyPep8Naming
def notest_InputFileName():
    io.InputFileName(read_file_1.name, check=True)
    with assert_raises(Exception) as ex:
        io.InputFileName("/tmp", check=True)
    assert_re_match(".*is a dir.*", ex.exception.args[0])
    with assert_raises(Exception) as ex:
        io.InputFileName("/tmp/this_should_not_exist_d9js9d$HHx", check=True)
    assert_re_match(".*does not exist.*", ex.exception.args[0])
    tmpfile = NamedTemporaryFile()
    os.chmod(tmpfile.name, 0)
    with assert_raises(Exception) as ex:
        io.InputFileName(tmpfile.name, check=True)
    assert_re_match(".*cannot be read.*", ex.exception.args[0])


# noinspection PyPep8Naming
def notest_OutputFileName():
    io.OutputFileName("/tmp/valid_output_file_lkjad9k", check=True)
    with assert_raises(Exception) as ex:
        io.OutputFileName("/tmp")
    assert_re_match(".*is a dir.*", ex.exception.args[0])
    with assert_raises(Exception) as ex:
        io.OutputFileName("/sbin/cannot_write_here")
    assert_re_match(".*directory is not writable.*", ex.exception.args[0])
    with assert_raises(Exception) as ex:
        io.OutputFileName("/this_path_no_exists/filename")
    assert_re_match(".*directory does not exist.*", ex.exception.args[0])
    tmpfile = NamedTemporaryFile()
    with assert_raises(Exception) as ex:
        io.OutputFileName(tmpfile.name, overwrite=False)
    assert_re_match(".*cowardly refusing to overwrite.*", ex.exception.args[0])
    os.chmod(tmpfile.name, 0)
    with assert_raises(Exception) as ex:
        io.OutputFileName(tmpfile.name, overwrite=True)
    assert_re_match(".*write protected.*", ex.exception.args[0])


def test_decompressed():
    data = ["this is a test\n", "and a second line\n"]
    methods = [("gzip", "gz"), ("lz4", "lz4"),
               ("xz", "xz"), ("bzip2", "bz2")]

    for compressor, suffix in methods:
        src = NamedTemporaryFile(suffix="." + suffix)
        try:
            zipper = io.Popen([compressor, "-c"], stdin=PIPE, stdout=src).stdin
        except OSError as e:
            if e.errno == 2:  # ENOENT
                if io.command_avail(compressor):
                    raise
            WARNING('Command "{}" not found. Skipping test.'.format(compressor))
            continue

        zipper.writelines(data)
        zipper.close()
        src.flush()

        with io.decompressed(src.name).reader() as f:
            for orig, test in zip(data, f):
                assert_equal(orig, test)


# noinspection PyPep8Naming
def notest_EnumerateReads():
    enumerated_reads = io.EnumerateReads(read_file_1)
    i = 0
    with enumerated_reads as reads:
        for _ in reads:
            i += 1
    with enumerated_reads as reads:
        for _ in reads:
            i += 1
    with enumerated_reads as reads, read_file_1 as orig_reads:
        cmp_reindexed_fq_files(orig_reads, reads, 50000)


def test_FastqCountReads():
    n_reads = io.fastq_count_reads(read_file_1)
    assert_equal(n_reads, 50000)


def test_Pipe_simple():
    class EchoPipe(io.Pipe):
        cmd = ["echo", "test word"]

    with EchoPipe() as p:
        line = p.next()
        assert_equal(line.strip(), "test word")


def notest_Pipe_chained():
    with io.Gunzip(io.Gzip(read_file_1)) as f, read_file_1 as g:
        for f_line, g_line in zip(f, g):
            assert_equal(f_line, g_line)


# noinspection PyPep8Naming
def notest_Pipe_cmd_subst():
    cmd = ["cat", io.EnumerateReads(read_file_1),
           io.EnumerateReads(read_file_1)]
    pipe = io.make_pipe("test", cmd)
    i = 0
    with pipe(None) as f:
        for _ in f:
            i += 1

    assert_equal(i, 400000)


def test_Pipe_cmd_stderr():
    cmd = ["bash", "-c", "for n in 0 1 2 3; do echo $n >&2; done"]
    pipe = io.make_pipe("test", cmd)
    with pipe(stdout=None, stderr=PIPE) as p:
        for a, b in zip(p, range(4)):
            assert_equal(a.strip(), str(b))


def test_Pipe_outfile():
    tmp = NamedTemporaryFile()
    pipe = io.make_pipe("test", ["echo", "test123"])
    with pipe(stdout=io.File(tmp.name)):
        pass
    with open(tmp.name) as f:
        assert_equal(f.next().strip(), "test123")


def test_Cat():
    with io.Cat(sam_file) as f:
        n = sum(1 for _ in f)
        assert_equal(n, 24)


def test_AlignmentFile():
    def check_AlignmentFile(obj, mode):
        with io.AlignmentFile(obj, mode) as af:
            n = sum(1 for _ in af.fetch(until_eof=True))
            assert_equal(n, 2)

    for obj in str, io.File, io.Cat:
        yield check_AlignmentFile, obj(bam_file), "rb"
        yield check_AlignmentFile, obj(sam_file), "r"


def test_Kseq():
    reads = 0
    bp = 0
    x = io.Kseq(read_file_1)
    with io.Kseq(read_file_1) as f:
        for seq in f:
            reads +=1
            assert_equal(len(seq.seq), len(seq.qual))
            bp += len(seq.seq)

    assert_equal(reads, 50000)
    assert_equal(bp, 50000*150)


if __name__ == '__main__':
    test_Pipe_cmd_subst()
