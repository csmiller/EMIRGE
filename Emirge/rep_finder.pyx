from cpython.mem cimport PyMem_Malloc, PyMem_Realloc, PyMem_Free
from libc.string cimport memset

cdef extern from "rep_finder.h":
    void find_repeats(unsigned int*, int, int, char*, int, int*, int)

cdef class RepFinder:
    cdef int k
    cdef unsigned int* kmer_map

    def __cinit__(self, k=16):
        self.k = k
        # we need bits to store whether we've seen any of the
        # 4^k = 2^2k possible kmers. since a byte has 8 or 2^3 bits,
        # we need 2^(2k-3) bytes:
        map_size = 2**(2*k -3)
        self.kmer_map = <unsigned int*> PyMem_Malloc(map_size)
        if not self.kmer_map:
            raise MemoryError()
        memset(self.kmer_map, 0, map_size)

    def __dealloc__(self):
        PyMem_Free(self.kmer_map)

    def set_k(self, k):
        self.k = k
        map_size = 2**(2*k -3)
        cdef unsigned int* new_mem
        new_mem = <unsigned int*> PyMem_Realloc(self.kmer_map, map_size)
        if not new_mem:
            raise MemoryError()
        self.kmer_map = new_mem
        memset(self.kmer_map, 0, map_size)

    def check(self, bytes sequence, int minlen=50):
        cdef int maxresults = 91
        cdef int[91] results
        find_repeats(self.kmer_map, self.k, minlen,
                     sequence, len(sequence),
                     results, maxresults)

        res = []
        for a,b,c in zip(*[iter(results)]*3):
            if a == 0: break
            res.append([a,b,c])

        return res

