"""
_emirge.pyx contains helper cython functions for emirge.py

EMIRGE: Expectation-Maximization Iterative Reconstruction of Genes from the Environment
Copyright (C) 2010-2016 Christopher S. Miller  (christopher.s.miller@ucdenver.edu)

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>

https://github.com/csmiller/EMIRGE
"""

from libc.stdlib cimport malloc, free, atoi
from libc.math cimport M_E, log, pow as c_pow
from cpython.tuple cimport PyTuple_GetItem
import numpy as np
cimport numpy as np
import sys

# from time import time  # DEBUG

# for lookup table of qual values
cdef extern from *:
    ctypedef double* static_const_double_ptr "static const double*"
cdef extern from "_emirge_C.h":
    static_const_double_ptr qual2one_minus_p   # lookup tables to avoid too many double calcs.
    static_const_double_ptr qual2p_div_3

cpdef int base_alpha2int(char base_ascii):
    """
    base2i = {"A":0,"T":1,"C":2,"G":3}
    """
    if base_ascii == <int>('A'):
        return 0
    elif base_ascii == <int>('T'):
        return 1
    elif base_ascii == <int>('C'):
        return 2
    elif base_ascii == <int>('G'):
        return 3
    else:
        return 4
        
cpdef double calc_likelihood_cell(unsigned int seq_i, unsigned int read_i, unsigned int pair_i,
                                  unsigned int pos,
                                  np.ndarray[np.uint8_t, ndim=1] numeric_bases,
                                  np.ndarray[np.uint8_t, ndim=1] qualints,
                                  np.ndarray[np.float_t, ndim=2] probN):  # an individual probN numpy matrix
    """
    numeric_bases is an array of unsigned char (numpy.ubyte)  already converted to 0-4-based index
    
    return single entry for lik_data, keeping main looping in python still
    RETURNS: double to place in lik_data
    """
    cdef double p = 0.0
    cdef double s

    cdef int i
    cdef int j

    for i in range(numeric_bases.shape[0]):
    # for i from 0 <= i < numeric_bases.shape[0]:  # old syntax, range syntax above fast if i is cdef-ed.
        s = 0.0
        for j in range(4):
            if numeric_bases[i] == j:   # this is called base, set 1-P
                # s += (1. - error_p_i) * probN[pos + i, j]
                s += ( qual2one_minus_p[qualints[i]] * probN[pos + i, j] )    # lookup table
            else:                       # this is not the called base, so set to P/3
                # s += (error_p_i / 3 ) * probN[pos + i, j]
                s += ( qual2p_div_3[qualints[i]] * probN[pos + i, j] )          # lookup table
        p += log(s)
        
    # do at least this product in log space (0.94005726833168002 vs. 0.94005726833167991)
    # likelihood[seq_i, read_i] = e**(np_log(prob_b_i).sum())
    return c_pow(M_E, p)

def _calc_likelihood(np.ndarray[np.int_t, ndim=2] bamfile_data,
                     list reads,
                     list quals,
                     list probN,
                     np.ndarray[np.uint_t, ndim=1] lik_row_seqi,
                     np.ndarray[np.uint_t, ndim=1] lik_col_readi,
                     np.ndarray[np.float_t, ndim=1] lik_data):
    """
    do looping here in Cython

    return values are placed in the arrays passed as the last 3 args.
    """
    cdef int alignedread_i
    cdef int seq_i
    cdef int read_i
    cdef int pair_i
    cdef int rlen
    cdef int pos

    cdef double p
    cdef double s

    cdef int i
    cdef int j

    cdef np.ndarray[np.uint8_t, ndim=1] numeric_bases
    cdef np.ndarray[np.uint8_t, ndim=1] qualints
    cdef np.ndarray[np.float_t, ndim=2] probN_single  # an individual probN numpy matrix

    for alignedread_i in range(bamfile_data.shape[0]):
        seq_i, read_i, pair_i, rlen, pos = bamfile_data[alignedread_i]
        lik_row_seqi[alignedread_i] = seq_i
        lik_col_readi[alignedread_i] = read_i

        numeric_bases = reads[read_i]
        qualints      = quals[read_i]
        probN_single  = probN[seq_i]

        p = 0.0

        for i in range(rlen):   # numeric_bases.shape[0]):
            s = 0.0
            for j in range(4):
                if numeric_bases[i] == j:   # this is called base, set 1-P
                    # s += (1. - error_p_i) * probN[pos + i, j]
                    s += ( qual2one_minus_p[qualints[i]] * probN_single[pos + i, j] )    # lookup table
                else:                       # this is not the called base, so set to P/3
                    # s += (error_p_i / 3 ) * probN[pos + i, j]
                    s += ( qual2p_div_3[qualints[i]] * probN_single[pos + i, j] )          # lookup table
            p += log(s)

        lik_data[alignedread_i] = c_pow(M_E, p)

    return

def calc_probN_read(bint initial_iteration,
                    unsigned int seq_i, unsigned int read_i,
                    unsigned int pos,
                    np.ndarray[np.float_t, ndim=1] priors,
                    posteriors,  # pass as python object for now
                    np.ndarray[np.uint8_t, ndim=1] numeric_bases,
                    np.ndarray[np.uint8_t, ndim=1] qualints,
                    np.ndarray[np.float_t, ndim=2] probN):  # an individual probN numpy matrix. Passed by ref?
    """
    calculates effect of a single read on probN

    this is inside loop over bamfile_data
    """
    cdef int i
    cdef int j
    cdef double weight

    ## posteriors is a sparse matrix, so this presents a little bit of a tricky situation until
    ## I can think about how to pass this.  Have to choose correct data struct, manip internals?
    ## for now, just keep as python obj.
    
    if initial_iteration:
        weight = priors[seq_i]
    else:
        if seq_i < posteriors.shape[0] and read_i < posteriors.shape[1]:
            weight = posteriors[seq_i, read_i]
        else:
            weight = priors[seq_i]

    for i in range(numeric_bases.shape[0]):
    # for i from 0 <= i < numeric_bases.shape[0]:  # old syntax, range syntax above fast if i is cdef-ed.
        s = 0.0
        for j in range(4):
            if numeric_bases[i] == j:   # this is called base, add (1-P) * weight
                probN[pos + i, j] += qual2one_minus_p[qualints[i]] * weight
            else:                       # this is not the called base, so add P/3 * weight
                probN[pos + i, j] += qual2p_div_3[qualints[i]] * weight

    return

def _calc_probN(np.ndarray[np.int_t, ndim=2] bamfile_data,
                bint initial_iteration,
                np.ndarray[np.float_t, ndim=1] priors,
                posteriors,  # pass as python object for now
                list numeric_bases,
                list qualints,
                list probN):  # list of  probN numpy matrices
    """
    calculates effect of *all* reads on probN

    entire loop over bamfile_data.
    ON RETURN:  probN is modified for this iteration.

    ## posteriors is a sparse matrix, so this presents a little bit of a tricky situation until
    ## I can think about how to pass this.  Have to choose correct data struct, manip internals?
    ## for now, just keep as python obj.
    """
    cdef int alignedread_i
    cdef int seq_i
    cdef int read_i
    cdef int pair_i
    cdef int rlen
    cdef int pos

    cdef int i
    cdef int j
    cdef double weight

    cdef np.ndarray[np.uint8_t, ndim=1] numeric_bases_single # single numeric bases numpy array -- per read
    cdef np.ndarray[np.uint8_t, ndim=1] qualints_single      # # single quals numpy array -- per read
    cdef np.ndarray[np.float_t, ndim=2] probN_single  # an individual probN numpy matrix.

    cdef int posteriors_shape_0 = 0
    cdef int posteriors_shape_1 = 0
    if posteriors is not None:
        posteriors_shape_0 = posteriors.shape[0]
        posteriors_shape_1 = posteriors.shape[1]
    


    for alignedread_i in range(bamfile_data.shape[0]):
        seq_i, read_i, pair_i, rlen, pos = bamfile_data[alignedread_i]
        # about 1/3 of time is spent doing this weight step
        if initial_iteration:
            weight = priors[seq_i]
        else:
            if seq_i < posteriors_shape_0 and read_i < posteriors_shape_1:
                weight = posteriors[seq_i, read_i]
            else:
                weight = priors[seq_i]

        numeric_bases_single = numeric_bases[read_i]
        qualints_single      = qualints[read_i]
        probN_single         = probN[seq_i]
        
        for i in range(numeric_bases_single.shape[0]): # for i in range(rlen)  # ????
            s = 0.0
            for j in range(4):
                if numeric_bases_single[i] == j:   # this is called base, add (1-P) * weight
                    probN_single[pos + i, j] += qual2one_minus_p[qualints_single[i]] * weight
                else:                       # this is not the called base, so add P/3 * weight
                    probN_single[pos + i, j] += qual2p_div_3[qualints_single[i]] * weight
    return


def count_cigar_aln(char* query_seq, char* hit_seq,
                    np.ndarray[np.uint8_t, ndim=1] query_unmapped_bases,
                    np.ndarray[np.uint8_t, ndim=1] hit_unmapped_bases,
                    unsigned int query_i,  # query_start -- 0-based
                    unsigned int hit_i,    # hit_start -- 0-based
                    alncode_list):
    """
    alncode list comes from a cigar string (see call in emirge_amplicon.py)

    returns number of aligned columns and number of matches in those aligned columns as tuple
    only includes columns where both sequences have mapped bases
    """

    # cdef unsigned int query_i = query_start  # 0   # for global alignment, this used to always be zero.
    # cdef unsigned int hit_i   = hit_start    # 0
    cdef int matches = 0
    cdef int aln_columns = 0
    cdef int count
    # cdef char aln_code  # doesn't work.  Expects integer
    cdef int i
    cdef int j

    for i in range(len(alncode_list)):
        count_a, aln_code = alncode_list[i]
        # print >> sys.stderr, "***", count_a, aln_code, len(query_seq), len(query_unmapped_bases), len(hit_seq), len(hit_unmapped_bases)
        if not len(count_a):
            count = 1  # "As a special case, if n=1 then n is omitted"
        else:
            count = int(count_a)

        if aln_code == 'M':    # match (letter-letter column; could be mismatch)
            for j in range(count):
                # print >> sys.stderr, "***", query_i+j, hit_i+j, query_seq[query_i+j], hit_seq[hit_i+j]
                
                if query_unmapped_bases[query_i+j] == 0:
                    if hit_unmapped_bases[hit_i+j] == 0:
                        aln_columns += 1
                        if query_seq[query_i+j] == hit_seq[hit_i+j]:
                            matches += 1
            query_i += count
            hit_i += count
        elif aln_code == 'I':
            query_i += count
        elif aln_code == 'D':
            hit_i += count
        else:
            raise ValueError, "unknown alignment code: '%s'"%aln_code

    return aln_columns, matches

def process_bamfile_predata(int seq_i, int read_i,
                            list predata, tuple samfile_references,
                            dict sequence_name2sequence_i, dict sequence_i2sequence_name,
                            dict read_name2read_i, dict read_i2read_name,
                            list reads, list quals, list readlengths,
                            list coverage, int ascii_offset,
                            np.ndarray[np.int_t, ndim=2] bamfile_data):
    """
    returns updated:
    (seq_i, read_i)
    """
    
    cdef int alignedread_i
    cdef int pos
    cdef int tid

    cdef int seq_i_to_cache
    cdef int pair_i
    cdef int read_i_to_cache

    for alignedread_i in range(len(predata)):
        pos, tid, pair_i, qname, qual, seq = predata[alignedread_i]
        # refname = PyTuple_GetItem(samfile_references, tid)  # Cannot convert 'PyObject *' to Python object
        refname = samfile_references[tid]
        seq_i_to_cache = sequence_name2sequence_i.get(refname, seq_i)
        if refname not in sequence_name2sequence_i:  # new sequence we haven't seen before
            sequence_name2sequence_i[refname] = seq_i
            sequence_i2sequence_name[seq_i] = refname
            coverage.append(0)
            seq_i += 1
        readname = "%s/%d"%(qname, pair_i+1)
        read_i_to_cache = read_name2read_i.get(readname, read_i)
        if readname not in read_name2read_i: # new read we haven't seen before
            read_name2read_i[readname] = read_i
            read_i2read_name[read_i] = readname
            read_i += 1
            # add to self.reads and self.quals and self.readlengths
            readlengths.append(len(seq))
            reads.append(seq_alpha2int(seq.encode('utf8'), len(seq)))
            quals.append(quals_alpha2int(qual.encode('utf8'), len(qual), ascii_offset))
            # reads.append(array([base_alpha2int(x) for x in fromstring(seq, dtype=uint8)], dtype=uint8))
            # quals.append(np.fromstring(qual, dtype=np.uint8) - ascii_offset)

        coverage[seq_i_to_cache] += len(seq)
        bamfile_data[alignedread_i] = [seq_i_to_cache, read_i_to_cache, pair_i, len(seq), pos]
    return seq_i, read_i

cdef np.ndarray[np.uint8_t, ndim=1] seq_alpha2int(char* seq, int seqlen):
    cdef int i
    np_as_ints = np.empty(seqlen, dtype=np.uint8)
    for i in range(seqlen):
        np_as_ints[i] = base_alpha2int(seq[i])
    return np_as_ints

cdef np.ndarray[np.uint8_t, ndim=1] quals_alpha2int(char* qual, int quallen, int ascii_offset):
    cdef int i
    np_as_ints = np.empty(quallen, dtype=np.uint8)
    for i in range(quallen):
        np_as_ints[i] = qual[i] - ascii_offset
    return np_as_ints

# def adjust_posteriors_for_split(posteriors, int seq_i, double major_fraction_avg):
#     raise NotImplementedError
#     return minor_data, minor_row, minor_col
