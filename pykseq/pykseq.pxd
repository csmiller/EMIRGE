import sys
cimport kseq
from libc.stdio cimport FILE, fopen, fclose, fread, const_char
from zlib cimport gzFile, gzopen, gzclose, gzrewind

ctypedef kseq.kseq_t *kseq_t_pointer

cdef class Kseq:
    cdef gzFile *_c_file
    cdef int l
    cdef kseq.kseq_t *seq

    cdef kseq.kseq_t* c_read_sequence(self)
    cpdef rewind(self)
