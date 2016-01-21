cdef extern from "kseq.h":
    cdef struct __kstring_t:
        int l
        int m
        char *s
    ctypedef __kstring_t kstring_t

    cdef struct __kstream_t:
        char *buf
        int begin
        int end
        int is_eof
        int fd
    ctypedef __kstream_t kstream_t

    ctypedef struct kseq_t:
        kstring_t name
        kstring_t comment
        kstring_t seq
        kstring_t qual
        int last_char
        kstream_t *f

    void _KSEQ_INIT "KSEQ_INIT"()
    kseq_t*kseq_init(int fd)
    int kseq_read(kseq_t *seq)


cdef class Kseq:
    cdef file
    cdef int fd
    cdef int l
    cdef kseq_t *ks

    cdef Kseq open(self)
    cdef close(self, exc_type= *, exc_val= *, exc_tb= *)
    cdef int read_next_sequence(self)
