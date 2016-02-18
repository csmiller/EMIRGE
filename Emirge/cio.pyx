cimport cython
from libc.stdio cimport FILE, fdopen, fclose, fread, fwrite, getline
from libc.stdlib cimport malloc, free

from Emirge.log import INFO, timed

@timed("Re-numbering reads")
def enumerate_fastq(object filein, object fileout):
    """
    Renames reads with consecutive numbers. Also gathers
    max read length, number of reads and phred type

    :param filein: Filelike object to read from
    :param fileout: Filelike object to write to
    :return: tuple (number of reads, max read length, is phred33)
    """
    cdef FILE *ptr_in  # input FILE pointer
    cdef FILE *ptr_out  # output FILE pointer
    cdef char *ptr_linebuf = <char *> 0  # getline buffer pointer
    cdef size_t linebuf_len = 0  # getline buffer length

    cdef const char *plus = "+\n"  # string to write for + line
    cdef bytes header  # for writing header line
    cdef char *ptr_header
    cdef char c  # for iterating over quality line
    cdef int i
    cdef Py_ssize_t res  # return code from fread/write
    cdef int lineno = 0  # line number
    cdef int fqline = 0  # line number modulo 4
    cdef bint phred33 = False  # have we found phred33 characters?
    cdef int max_readlen = 0  # max read length

    ptr_in = fdopen(filein.fileno(), "rb")
    ptr_out = fdopen(fileout.fileno(), "wb")

    while True:
        res = getline(&ptr_linebuf, &linebuf_len, ptr_in)
        if res <= 0:
            break

        fqline = lineno % 4
        if fqline == 0:  # @header
            header = "@{}\n".format(lineno/4)
            ptr_header = header
            res = fwrite(ptr_header, 1, len(header), ptr_out)
        elif fqline == 1:  # sequence
            res = fwrite(ptr_linebuf, 1, res, ptr_out)
        elif fqline == 2:  # +
            res = fwrite(plus, 1, 2, ptr_out)
        elif fqline == 3:  # quality
            if not phred33:
                for i in range(res-1):
                    if ptr_linebuf[i] < 64:
                        phred33=True
                        break
            max_readlen = max(max_readlen, res-1)
            res = fwrite(ptr_linebuf, 1, res, ptr_out)

        lineno += 1

    if ptr_linebuf != NULL:
        free(ptr_linebuf)

    # no closing of FDs. That is left to file object.

    return lineno / 4, max_readlen, phred33

@timed("Detecting read length and quality encoding")
def fastq_detect_len_and_encoding(object filein, int num=1000):
    """
        Parses first num lines of fastq file to guess max read length and quality
        encoding.

        :param filein: Input file object, needs to have fileno()
        :param num: Number of reads to check from beginning of file
        :return: (max_readlen, is_phred33)
        """
    cdef FILE *ptr_in  # input FILE pointer
    cdef char *ptr_linebuf = <char *> 0  # getline buffer pointer
    cdef size_t linebuf_len = 0  # getline buffer length
    cdef int lineno = 0  # line number
    cdef bint phred33 = False  # have we found phred33 characters?
    cdef int max_readlen = 0  # max read length

    ptr_in = fdopen(filein.fileno(), "rb")

    while True:
        res = getline(&ptr_linebuf, &linebuf_len, ptr_in)
        if res <= 0:
            break

        if lineno % 4 == 3:  # quality line
            if not phred33:
                for i in range(res-1):
                    if ptr_linebuf[i] < 64:
                        phred33=True
                        break
            max_readlen = max(max_readlen, res-1)

        lineno += 1
        if lineno > num * 4:
            break

    return max_readlen, phred33
