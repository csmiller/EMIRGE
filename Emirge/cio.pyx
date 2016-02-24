import os

from libc.string cimport strchr
from libc.stdio cimport fdopen, fwrite, getline, FILE, fflush, ferror, \
    sprintf
from libc.stdlib cimport free
from libc.errno cimport errno

from Emirge.log import timed


@timed("Pre-processing input FastQ file {self.name}")
def fastq_process(object filein, object fileout=None, set keep=None, int num=0):
    """
    Renumbers, filters, and inspects reads.

    Reads from :param filein are copied to :param fileout if their
    ID ("@<ID>[/1|/2]") is in :param keep or :param keep is None.
    Returns the number of reads copied, the maximum length of all reads
    copied and whether any phred33 encoding characters where encountered in
    the quality lines.

    :param filein: Input file object. Must have fileno() property
    :param fileout: Output file object. Must have fileno() property (optional)
    :param keep: Set of IDs to copy (optional)
    :param num: Max number of reads to copy / process (optional)
    :return: Tuple of n_reads_copied, max_readlen, phred33
    """
    cdef FILE *ptr_in = NULL # input FILE pointer
    cdef FILE *ptr_out = NULL # output FILE pointer
    cdef char *ptr_linebuf = <char *> 0  # getline buffer pointer
    cdef size_t linebuf_len = 0  # getline buffer length
    cdef char *ptr_towrite = <char *> 0  # fwrite buffer pointer
    cdef size_t towrite_len = 0  # fwrite length

    cdef const char *plus = "+\n"  # string to write for + line
    cdef int slash = ord('/')
    cdef bytes header  # for writing header line
    cdef char *ptr_header
    cdef char c  # for iterating over quality line
    cdef int i
    cdef Py_ssize_t res  # return code from fread/write
    cdef int lineno = -1  # line number
    cdef int fqline = 0  # line number modulo 4
    cdef bint phred33 = False  # have we found phred33 characters?
    cdef int max_readlen = 0  # max read length
    cdef int n_reads = 0
    cdef int skip = 0

    ptr_in = fdopen(filein.fileno(), "r")
    if ptr_in is NULL:
        raise IOError(os.strerror(errno))

    if fileout is not None:
        ptr_out = fdopen(fileout.fileno(), "w")
        if ptr_in is NULL:
            raise IOError(os.strerror(errno))

    while num == 0 or n_reads < num:
        lineno += 1
        res = getline(&ptr_linebuf, &linebuf_len, ptr_in)
        if res < 0:
            break
        if skip > 0:
            skip -= 1
            continue

        ptr_towrite = ptr_linebuf
        towrite_len = res

        fqline = lineno % 4
        if fqline == 0:  # @header
            # fastq header lines are formatted as:
            # @ followed by alphanumeric+dots, optionally followed by /1 or /2 to
            # mark fwd/reverse reads, followed by line-end
            # => skip first character, strip whitespace, consider only part before
            #    first "/"
            assert res > 2
            if keep is not None:
                ptr_header = strchr(ptr_linebuf, slash)
                if ptr_header:
                    header = ptr_linebuf[1:ptr_header-ptr_linebuf]
                else:
                    header = ptr_linebuf[1:res-1]
                if header not in keep:
                    skip = 3
                    continue
            # rewrite to @<read-number>
            if ptr_out is not NULL:
                towrite_len = sprintf(ptr_linebuf, "@%i\n", n_reads)
                # header = "@{}\n".format(n_reads)
                # ptr_towrite = header
                # towrite_len = len(header)
        elif fqline == 2:  # +
            # just write "+"
            # (legacy line, but some tools like fastq-dump still put
            #  the read id after, no need to waste that space)
            ptr_towrite = plus
            towrite_len = 2
        elif fqline == 3:  # quality
            if not phred33:
                # if we find a character <64, we must be parsing
                # phred 33 encoded fastq
                for i in range(res-1):
                    if ptr_linebuf[i] < 64:
                        phred33=True
                        break
            max_readlen = max(max_readlen, res-1)
            n_reads += 1

        if ptr_out is not NULL:
            res = fwrite(ptr_towrite, 1, towrite_len, ptr_out)
            if res < towrite_len:
                raise IOError(os.strerror(errno))

    if res < 0 and ferror(ptr_in) != 0:
        raise IOError(os.strerror(errno))

    if ptr_linebuf != NULL:
        # free buffer allocated by getline
        free(ptr_linebuf)
    if ptr_out is not NULL:
        # flush output, in case it's not closed before it's read externally
        res = fflush(ptr_out)
        if res != 0:
            raise IOError(os.strerror(errno))

    return n_reads, max_readlen, phred33
