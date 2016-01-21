# distutils: libraries = z

"""
For wrapping Heng Li's kseq.h FASTA/FASTQ parser.

see http://lh3lh3.users.sourceforge.net/parsefastq.shtml

only deals with uncompressed files, for the moment.

Chris Miller
csmiller@gmail.com

"""
import os

from Emirge.io import FileLike

cdef class Kseq:
    def __cinit__(self, fileobj):
        self.file = fileobj

    def __enter__(self):
        self.open()
        return self

    cdef open(self):
        if isinstance(self.file, int):
            self.fd = self.file
        elif isinstance(self.file, FileLike):
            self.file.__enter__()
            self.fd = self.file.fileno()
        elif hasattr(self.file, "fileno"):
            self.fd = self.file.fileno()
        elif isinstance(self.file, str):
            self.fd = open(self.file, "r")

        self.ks = kseq_init(self.fd)

        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close(exc_type, exc_val, exc_tb)

    cdef close(self, exc_type=None, exc_val=None, exc_tb=None):
        # don't do anything if we were given an fd or a file object
        if isinstance(self.file, FileLike):
            self.file.__exit__(exc_type, exc_val, exc_tb)
        elif isinstance(self.file, str):
            os.close(self.fd)
        return False

    def __iter__(self):
        return self

    def __next__(self):
        cdef int l
        l = self.read_next_sequence()
        if l == -1:
            raise StopIteration()
        elif l == -2:
            raise ValueError("Malformed read '{}'".format(self.ks.name.s))
        return self

    cdef int read_next_sequence(self):
        return kseq_read(self.ks)

    def __getattr__(self, item):
        if item == "seq":
            return self.ks.seq.s
        elif item == "qual":
            return self.ks.qual.s
        elif item == "name":
            return self.ks.name.s
        elif item == "comment":
            return self.ks.comment.s
