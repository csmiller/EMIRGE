"""Implements IO functions"""

import re
import os
import weakref
import stat
from tempfile import NamedTemporaryFile, mkdtemp
import subprocess

from Emirge.log import ERROR, DEBUG, INFO, timed


PIPE = subprocess.PIPE


class FileError(Exception):
    """Raised when file operations fail"""
    pass


def ispipe(path):
    """Checks if file at @path is a pipe"""
    try:
        return stat.S_ISFIFO(os.stat(path).st_mode)
    except OSError:
        return False


class Record(object):
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


@timed("Rewriting reads with indices in headers")
def reindex_reads(reads_filepath):
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

    with decompressed(reads_filepath).reader() as src:
        try:
            check_call(['awk',
                        '{ if ((NR-1) %% 4 == 0) { print "@"(NR-1)/4;nr=nr+1 }'
                        'else { print } }'
                        'END { print nr > "%s" }'
                        % tmp_n_reads_file_path.name],
                       stdin=src, stdout=new_reads_filepath)
        except subprocess.CalledProcessError:
            ERROR("awk rewrite of reads failed! Is awk installed?")
            raise

    new_reads_filepath.file.seek(0)
    n_reads = int(tmp_n_reads_file_path.readline().strip())

    return new_reads_filepath, n_reads


class TempDir(object):
    """Creates a temporary directory
    The directory is cleaned automatically at object deletion."""
    def __init__(self, suffix="", prefix="tmp", dir=None):
        self.__suffix = suffix
        self.__prefix = prefix
        self.__dir = dir
        self.__name = None

    def __del__(self):
        if self.__name:
            DEBUG("Unlinking temporary directory \"%s\"" % self.__name)
            try:
                os.rmdir(self.__name)
            except OSError, e:
                ERROR("Unlinking failed: %s" % e.strerror)

    @property
    def name(self):
        if not self.__name:
            self.__name = mkdtemp(self.__suffix, self.__prefix, self.__dir)
            DEBUG("Created temporary directory \"%s\"" % self.__name)

        return self.__name


class FileName(object):
    def __init__(self, name):
        self.name = name

    def __str__(self):
        return self.name


class InputFileName(FileName):
    """Class to contain the path to an input file. Checks that file exists
    and is readable in constructor."""

    def __init__(self, name, check=True):
        super(InputFileName, self).__init__(name)
        self.mode = "rb"

        if not check:
            return

        # do some checking to tell user about potential problems at startup
        # (these errors still have to be checked later, as the accessibility
        # of the file may have changed until then)
        msg = 'Input file "{}" '.format(name)
        if not os.path.exists(name):
            raise Exception(msg + "does not exist")
        if os.path.isdir(name):
            raise Exception(msg + "is a directory (file expected)")
        if not os.path.isfile(name) and not ispipe(name):
            raise Exception(msg + "is neither file nor pipe")
        if not os.access(name, os.R_OK):
            raise Exception(msg + "cannot be read")


class OutputFileName(FileName):
    def __init__(self, name, check=True, overwrite=True):
        super(OutputFileName, self).__init__(name)
        self.mode = "wb+"

        if not check:
            return

        msg = 'Output file "{}" '.format(name)
        if os.path.isdir(name):
            raise Exception(msg + "is a directory (file expected)")
        if not os.path.exists(os.path.dirname(name)):
            raise Exception(msg + "cannot be created: directory does not exist")
        if not os.path.isfile(name) and not ispipe(name):
            if not os.access(os.path.dirname(name), os.W_OK):
                raise Exception(msg + "cannot be created: " +
                                      "directory is not writable")
        if os.path.isfile(name):
            if not overwrite:
                raise Exception(msg + "exists; cowardly refusing to overwrite")
            if not os.access(name, os.W_OK):
                raise Exception(msg + "is write protected")


class FileLike(object):
    """Base class for objects that behave like a file

    MUST OVERRIDE __enter__ and __exit__!

    I.e.:
     - can be opened, closed and have a closed property
     - can be iterated over, reading OR writing
     - have a file-descriptor / fileno() function

    """
    def __init__(self):
        self.__isReader = True
        self.__fileobj = None

    @property
    def isReader(self):
        return self.__isReader

    @property
    def closed(self):
        return self.__fileobj is None

    def reader(self):
        if self.__fileobj and not self.__isReader:
            raise FileError("Cannot change open writer to reader ")
        self.__isReader = True
        return self

    def writer(self):
        if self.__fileobj and self.__isReader:
            raise FileError("Cannot change open reader to writer")
        self.__isReader = False
        return self

    def close(self):
        if self.__fileobj is not None:
            DEBUG("closing {}".format(repr(self.__fileobj)))
            self.__fileobj.close()
            self.__fileobj = None

    def fileno(self):
        if self.__fileobj is None:
            raise FileError("No fd for closed File")
        return self.__fileobj.fileno()

    def __iter__(self):
        if self.__fileobj is None:
            raise FileError("Cannot read from closed File")
        if not self.__isReader:
            raise FileError("Cannot read from writer")
        return self.__fileobj

    def next(self):
        if self.__fileobj is None:
            raise FileError("Cannot read from closed File")
        if not self.__isReader:
            raise FileError("Cannot read from writer")
        return self.__fileobj.next()

    def write(self, string):
        if self.__fileobj is None:
            raise FileError("Cannot write to closed File")
        if self.__isReader:
            raise FileError("Cannot write to reader")
        return self.__fileobj.write(string)

    def writelines(self, lines):
        if self.__fileobj is None:
            raise FileError("Cannot write to closed File")
        if self.__isReader:
            raise FileError("Cannot write to reader")
        return self.__fileobj.writelines(lines)

    def _setFileObj(self, fileobj):
        self.__fileobj = fileobj

    def __enter__(self):
        raise NotImplementedError()

    def __exit__(self, exc_type, exc_val, exc_tb):
        raise NotImplementedError()


class File(FileLike):
    def __init__(self, filename):
        super(File, self).__init__()
        self.__filename = filename

    @property
    def name(self):
        return self.__filename

    def __enter__(self):
        if self.isReader:
            DEBUG("opening {} read".format(repr(self.__filename)))
            self._setFileObj(open(self.__filename, "rb"))
        else:
            DEBUG("opening {} write".format(repr(self.__filename)))
            self._setFileObj(open(self.__filename, "wb+"))
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()


class NamedPipe(File):
    """Creates a named pipe
    Named pipes are created in a shared temporary directory, pipes and
    directory are automatically removed once the NamedPipe object is destroyed

    - opening a named pipe for reading blocks until write open
    - opening a named pipe for writing blocks until read open
    - opening a named pipe for read/write (wb+) does not block
    - pipe will EOF once all writers are closed

    """
    numPipes = 0  # number of created named pipes
    pipes_dir = None  # holds WeakRef to TempDir object

    def __init__(self, suffix="", prefix="pipe", dir=None):
        self.__suffix = suffix
        self.__prefix = prefix
        self.__dir = dir
        self.filename = None

        self.__pipeNum = NamedPipe.numPipes
        NamedPipe.numPipes += 1

        if self.__dir is None:
            pipes_dir = NamedPipe.pipes_dir
            if pipes_dir is None or pipes_dir() is None:
                # no pipe dir created, or pipe dir deleted (weakref evaluates
                #  to None)
                # => create new TempDir and store weak reference in static
                #  pipes_dir
                self.__dir = TempDir("pipes")
                NamedPipe.pipes_dir = weakref.ref(self.__dir)
            else:
                # reference shared pipe directory
                self.__dir = pipes_dir()

        filename = "_".join([self.__prefix, str(self.__pipeNum),
                             self.__suffix])
        filename = os.path.join(self.__dir.name, filename)
        os.mkfifo(filename)
        DEBUG("Created named pipe \"%s\"" % filename)
        super(NamedPipe, self).__init__(filename)

    def __del__(self):
        if self.name:
            DEBUG("Unlinking named pipe \"%s\"" % self.name)
            os.unlink(self.name)


def _obj2str(obj):
    if hasattr(obj, "name"):
        return obj.name
    elif isinstance(obj, FileLike):
        raise Exception("what?")
    return obj

def _expand_args(args):
    return [_obj2str(x) for x in args]

class Popen(subprocess.Popen):
    """Convenience wrapper around subprocess.Popen
    - always has close_fds=True
    - can have File objects as argument
    - prints log when executing command
    """
    def __init__(self, args, *otherargs, **kwargs):
        args = _expand_args(args)
        cmdline = " ".join([str(x) for x in args])
        INFO("Executing '{}'".format(cmdline))
        kwargs["close_fds"] = True
        super(Popen, self).__init__(args, *otherargs, **kwargs)


def check_call(args, *otherargs, **kwargs):
    """Convenience wrapper around subprocess.check_call
    - can have File objects as argument
    - prints log when executing command
    """
    args = _expand_args(args)
    cmdline = " ".join([str(x) for x in args])
    INFO("Executing '{}'".format(cmdline))
    return subprocess.check_call(args, *otherargs, **kwargs)


class Pipe(FileLike):
    cmd = None
    name = "-"  # needed to make pysam accept this as a pipe

    def __init__(self, infile=PIPE, outfile=PIPE):
        super(Pipe, self).__init__()
        if isinstance(infile, FileLike):
            self.__input = infile.reader()
        else:
            self.__input = infile

        if isinstance(outfile, FileLike):
            self.__output = outfile.writer()
        else:
            self.__output = outfile

        self.__proc = None
        self.closed = True
        self.__input_entered = None
        self.__output_entered = None
        self.__tmppipe = None

    def __enter__(self):
        to_close = []
        # enter input
        if isinstance(self.__input, FileLike):
            self.__input_entered = self.__input.__enter__().fileno()
        elif hasattr(self.__input, "name"):
            self.__input_entered = open(self.__input.name, "rb+")
            to_close.append(self.__input_entered)
        else:
            self.__input_entered = self.__input

        # enter ouput
        if isinstance(self.__output, FileLike):
            self.__output_entered = self.__output.__enter__().fileno()
        elif isinstance(self.__output, NamedPipe):
            self.__output_entered = open(self.__output.name, "wb+")
            to_close.append(self.__output_entered)
        else:
            self.__output_entered = self.__output

        # replace Pipe objects with named pipes
        args = [x.enter_as_named_pipe() if isinstance(x, Pipe) else x
                for x in self.cmd]

        self.__proc = Popen(args, close_fds=True,
                            stdin=self.__input_entered,
                            stdout=self.__output_entered)

        for f in to_close:
            f.close()

        if self.__input == PIPE:
            self.infile = self.__proc.stdin
        if self.__output == PIPE:
            self.outfile = self.__proc.stdout
        self.closed = False

        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if isinstance(self.__output, FileLike):
            self.__output.__exit__(exc_type, exc_val, exc_tb)
        if isinstance(self.__input, FileLike):
            self.__input.__exit__(exc_type, exc_val, exc_tb)
        for x in self.cmd:
            if isinstance(x, Pipe):
                x.__exit__(exc_type, exc_val, exc_tb)
        self.closed = True
        self.__tmppipe = None
        self.__proc.wait()
        del self.__proc
        self.__proc = None
        del self.outfile
        self.outfile = None
        self.infile = None

    def enter_as_named_pipe(self):
        if self.__input == PIPE:
            self.__input = NamedPipe(suffix=type(self).__name__)
            self.__tmppipe = self.__input
            self.__enter__()
            self.__input == PIPE
        elif self.__output == PIPE:
            self.__output = NamedPipe(suffix=type(self).__name__)
            self.__tmppipe = self.__output
            self.__enter__()
            self.__output = PIPE
        else:
            raise Exception("named pipe with infile={} and outfile={}"
                            .format(self.infile, self.outfile))

        return self.__tmppipe.name


def make_pipe(name, args):
    def __init__(self, filename):
        Pipe.__init__(self, filename)
    return type(name, (Pipe,), {"__init__": __init__,
                                "cmd": args})


Gunzip = make_pipe("Gunzip", ["gzip", "-dc"])
Unxz = make_pipe("Unxz", ["xz", "-dc"])
Unlz4 = make_pipe("Unlz4", ["lz4", "-dc"])
Bunzip2 = make_pipe("Bunzip2", ["bzip2", "-dc"])

Gzip = make_pipe("Gzip", ["gzip", "-c"])


def decompressed(filething):
    if isinstance(filething, str):
        name = filething
        filething = File(filething)
    else:
        name = filething.name

    if name.endswith(".gz"):
        return Gunzip(filething)
    if name.endswith(".xz"):
        return Unxz(filething)
    if name.endswith(".lz4"):
        return Unlz4(filething)
    if name.endswith(".bz2"):
        return Bunzip2(filething)

    return filething


EnumerateReads = make_pipe(
    "EnumerateReads",
    ['awk',
     '{ if ((NR-1) % 4 == 0) { print "@"(NR-1)/4 }'
     'else { print } }']
)
LineCount = make_pipe("LineCount", ['wc', '-l'])


@timed("Counting reads in input files")
def fastq_count_reads(filename):
    with LineCount(decompressed(filename)) as f:
        return int(f.next().strip()) / 4


def filter_fastq(infile, readnames, outfile=None):
    if outfile is None:
        outfile = NamedTemporaryFile(suffix="filtered.fq")

    matches = reads = 0
    for line in infile:
        # fastq header lines are formatted as:
        # @ followed by alphanumeric+dots, optionally followed by /1 or /2 to
        # mark fwd/reverse reads, followed by line-end
        # => skip first character, strip whitespace, consider only part before
        #    first "/"
        if line[1:].strip().split("/")[0] in readnames:
            outfile.write(line)
            for n in range(3):
                outfile.write(infile.next())
            matches += 1
        else:
            for n in range(3):
                infile.next()
        reads += 1
        if reads < 10:
            DEBUG("read {}: '{}'".format(reads, line.strip()))

    outfile.seek(0)
    return outfile, reads, matches
