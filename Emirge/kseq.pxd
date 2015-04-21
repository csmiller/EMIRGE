"""
For exporting Heng Li's kseq.h FASTA/FASTQ parser.

see http://lh3lh3.users.sourceforge.net/parsefastq.shtml

only deals with uncompressed files, for the moment.

Chris Miller
csmiller@gmail.com

"""
from zlib cimport gzFile

cdef extern from "kseq.h":
    cdef struct __kstring_t:
        int l
        int m
        char *s
    ctypedef __kstring_t kstring_t
    cdef struct __kstream_t:
        int begin
        int end
        int is_eof
    ctypedef __kstream_t kstream_t

    ctypedef struct  kseq_t:
        kstring_t name
        kstring_t comment
        kstring_t seq
        kstring_t qual
        int last_char
        kstream_t *f
    void _KSEQ_INIT "KSEQ_INIT" ()
    
    kseq_t* kseq_init(gzFile *fp)

    void kseq_rewind(kseq_t *seq)
    int kseq_read(kseq_t *seq)



    
    
