"""
_emirge_amplicon.pyx contains helper cython functions for emirge_amplicon.py

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

cimport cython
from cpython.dict cimport PyDict_SetItemString, PyDict_GetItemString
from libc.stdlib cimport malloc, free, atoi
from libc.math cimport M_E, log, pow as c_pow
from cpython.tuple cimport PyTuple_GetItem
import numpy as np
cimport numpy as np
from scipy import sparse
from bisect import bisect_left
import sys
import os
from time import ctime, time
from datetime import timedelta
import pysam
# cimport ctrie
cimport pykseq

# for lookup table of qual values
cdef extern from *:
    ctypedef double* static_const_double_ptr "static const double*"
cdef extern from "_emirge_C.h":
    static_const_double_ptr qual2one_minus_p   # lookup tables to avoid too many double calcs.
    static_const_double_ptr qual2p_div_3

cdef inline int base_alpha2int(char base_ascii):
    """
    base2i = {"A":0,"T":1,"C":2,"G":3}
    """
    if base_ascii == 65:   # <int>('A'):
        return 0
    elif base_ascii == 84: # <int>('T'):
        return 1
    elif base_ascii == 67:  # <int>('C'):
        return 2
    elif base_ascii == 71:  # <int>('G'):
        return 3
    elif base_ascii == 78:  # <int>('N'):
        return 4
    else:
        print >> sys.stderr, "WARNING: invalid base in base_alpha2int: %s"%chr(base_ascii)
        return 4

cdef inline unsigned char complement_numeric_base(unsigned char c):
        if c == 0:
            return 1
        elif c == 1:
            return 0
        elif c == 2:
            return 3
        elif c == 3:
            return 2
        else:  # just set any non ACTG to N
            return 4

@cython.boundscheck(False)
def _calc_likelihood(em):
    """
    Cython helper function with typed data structures, fast looping.
    """
    t0 = time()
    cdef np.ndarray[np.uint32_t, ndim=2] bamfile_data = em.bamfile_data
    cdef np.ndarray[np.uint8_t, ndim=3] reads = em.reads
    cdef np.ndarray[np.uint8_t, ndim=3] quals = em.quals                    
    cdef list probN = em.probN

    # probN_array_of_pointers = 

    # result arrays for constructing coo_matrix for likelihood.
    cdef bint paired_reads
    if em.reads2_filepath is not None:
        l = len(em.bamfile_data) / 2
        paired_reads = True
    else:
        l = len(em.bamfile_data)
        paired_reads = False
    cdef np.ndarray[np.uint_t, ndim=1]  lik_row_seqi  = np.empty(l, dtype=np.uint)
    cdef np.ndarray[np.uint_t, ndim=1]  lik_col_readi = np.empty_like(lik_row_seqi)
    cdef np.ndarray[np.float_t, ndim=1] lik_data      = np.zeros(l, dtype=np.float)
                     
    cdef unsigned int alignedread_i
    cdef unsigned int seq_i
    cdef unsigned int read_i
    cdef unsigned int pair_i
    cdef unsigned int rlen
    cdef unsigned int pos
    cdef bint is_reverse

    cdef double p
    cdef double s
    cdef double result

    cdef unsigned int i
    cdef unsigned int ii
    cdef unsigned int j

    # keep track of whether we have collected 2 of 2 pairs or not.
    # start at 1 because of way I flip this flag every iteration
    # below. (a bit counter-intuitive, but fastest thing I could think
    # of)  Assumption is that every two reads is a mapped pair spit out
    # from bowtie in bamfile
    cdef bint pair_collected = 1        
    cdef unsigned int result_i       = 0
    cdef np.ndarray[np.float_t, ndim=2] probN_single  # an individual probN numpy matrix

    cdef np.ndarray[np.uint8_t, ndim=1] dead_seqs # these sequences have been removed previously
    dead_seqs = np.zeros(len(em.probN), dtype=np.uint8)
    for seq_i in range(len(em.probN)):
        if em.probN[seq_i] is None:
            dead_seqs[seq_i] = 1

    p = 0.0  # start off with no pair collected
    loop_t0 = time()
    for alignedread_i in range(bamfile_data.shape[0]):
        if paired_reads:
            pair_collected = not pair_collected  # flip flag to indicate we've got 2/2 reads accumulated in p

        #   0       1       2      3     4       5   # "efficient indexing only affects certain index operations, namely those with exactly ndim number of typed integer indices”
        # seq_i, read_i, pair_i, rlen, pos, is_reverse = bamfile_data[alignedread_i]
        seq_i = bamfile_data[alignedread_i, 0]
        read_i = bamfile_data[alignedread_i, 1]

        if dead_seqs[seq_i]:
            result = 0.0        # still store a result so that we can have th coo arrays allocated with right size.  These get ignored with sparse construction anyway
        else: # else this is a good seq_i
            pair_i = bamfile_data[alignedread_i, 2]
            rlen = bamfile_data[alignedread_i, 3]
            pos = bamfile_data[alignedread_i, 4]
            is_reverse= bamfile_data[alignedread_i, 5]

            probN_single  = probN[seq_i]  # probably slow!! (list access)
            # probN_single = PyList_GET_ITEM(probN, seq_i)

            if is_reverse:
                # This code is redundant, but by basically manually inlining everything here in the special reverse-complement case, I gain some speed
                for i in range(rlen):
                    ii = rlen-1-i
                    s = 0.0
                    for j in range(4):
                        if complement_numeric_base(reads[read_i, pair_i, ii]) == j:   # this is called base, set 1-P
                            # s += (1. - error_p_i) * probN[pos + i, j]
                            s += ( qual2one_minus_p[quals[read_i, pair_i, ii]] * probN_single[pos + i, j] )    # lookup table
                        else:                       # this is not the called base, so set to P/3
                            # s += (error_p_i / 3 ) * probN[pos + i, j]
                            s += ( qual2p_div_3[quals[read_i, pair_i, ii]] * probN_single[pos + i, j] )          # lookup table
                    p += log(s)

            else: # not reverse
                for i in range(rlen):
                    s = 0.0
                    for j in range(4):
                        if reads[read_i, pair_i, i] == j:   # this is called base, set 1-P
                            # s += (1. - error_p_i) * probN[pos + i, j]
                            s += ( qual2one_minus_p[quals[read_i, pair_i, i]] * probN_single[pos + i, j] )    # lookup table
                        else:                       # this is not the called base, so set to P/3
                            # s += (error_p_i / 3 ) * probN[pos + i, j]
                            s += ( qual2p_div_3[quals[read_i, pair_i, i]] * probN_single[pos + i, j] )          # lookup table
                    p += log(s)
            if pair_collected:
                result = c_pow(M_E, p)
        
        if pair_collected:
            # for either good seq or dead seq, add the values to the coo construction arrays
            lik_row_seqi[result_i] = seq_i
            lik_col_readi[result_i] = read_i
            lik_data[result_i] = result
            p = 0.0  # reset for next pair or single read
            result_i += 1
    td = timedelta(seconds = time()-loop_t0)
    # print >> sys.stderr, "DEBUG: iteration %02d _calc_likelihood time per 1000 bam entries: %.2e seconds"%(em.iteration_i, (((td.microseconds + (td.seconds + td.days * 24 * 3600) * 10**6) / float(10**6)) / alignedread_i) * 1000.0) # DEBUG
    # now actually construct sparse matrix.
    em.likelihoods = sparse.coo_matrix((lik_data, (lik_row_seqi, lik_col_readi)), em.likelihoods.shape, dtype=em.likelihoods.dtype).tocsr()
    # print >> sys.stderr, "DEBUG: iteration %02d _calc_likelihood total cython function call time: %s seconds"%(em.iteration_i, timedelta(seconds = time()-t0)) # DEBUG
    return

@cython.boundscheck(False)
def _calc_probN(em):
    """
    helper function to emirge_amplicon.py calc_probN
    
    calculates effect of *all* reads on probN
    i.e. loop over entire bamfile_data.
    
    ON RETURN:  probN contains sums of weighted prob vals
    """
    cdef bint initial_iteration = em.iteration_i < 1#  or em.resume   # MARK if resume case, no posteriors, so treat like initial iter and use priors.
    
    cdef np.ndarray[np.uint8_t, ndim=3] reads = em.reads  # shape = (self.n_reads, 2, self.max_read_length) (dim 1 = pair)
    cdef np.ndarray[np.uint8_t, ndim=3] quals = em.quals
    cdef np.ndarray[np.uint32_t, ndim=2] bamfile_data = em.bamfile_data

    # if no posteriors are present (or initial iter), falls back on priors from previous round
    cdef np.ndarray[np.float_t, ndim=1] priors = np.array(em.priors[-2], dtype = np.float)    

    ## posteriors is a sparse matrix, so this presents a little bit of a tricky situation until
    ## I can think about how to get to C level.  Have to choose correct data struct, manip internals?
    ## for now, just keep as python obj.

    # these are for manipulating internals of csr type sparse matrix
    cdef np.ndarray[np.float_t, ndim=1] posteriors_data
    cdef np.ndarray[np.int32_t, ndim=1] posteriors_indices
    cdef np.ndarray[np.int32_t, ndim=1] posteriors_indptr
    
    if not initial_iteration:
        convert_t0 = time()
        em.posteriors[-2] = em.posteriors[-2].tocsr()
        em.posteriors[-2].sort_indices()
        # print >> sys.stderr, "DEBUG: Time for conversion of posterios to csr sparse matrix: %s"%(timedelta(seconds = time() - convert_t0))
        posteriors = em.posteriors[-2]  # this depends on PREVIOUS iteration's posteriors (seq_n x read_n)
        posteriors_data = posteriors.data
        posteriors_indices = posteriors.indices
        posteriors_indptr = posteriors.indptr
    else:
        posteriors = None

    t_check = time()                    # DEBUG
    cdef list probN      = em.probN       # list of numpy arrays.  At this point, only None values will be due to just-culled sequences.

    # for bamfile_data
    cdef unsigned int alignedread_i
    cdef unsigned int seq_i
    cdef unsigned int read_i
    cdef unsigned int pair_i
    cdef unsigned int rlen
    cdef unsigned int pos
    cdef bint is_reverse

    cdef unsigned int i
    cdef unsigned int ii
    cdef unsigned int j
    cdef int start # for indexing sparse matrix data structures
    cdef int end
    cdef int low  # binary search in sparse matrix data struct
    cdef int mid
    cdef int hi
    cdef int midval
    
    cdef double weight

    cdef np.ndarray[np.float_t, ndim=2] probN_single  # an individual probN numpy matrix.
    # cdef np.ndarray[np.uint8_t, ndim=1] new_read      # see below
    cdef np.ndarray[np.uint8_t, ndim=1] reads_seen = em.reads_seen      # ...in previous rounds
    cdef np.ndarray[np.uint8_t, ndim=1] reads_seen_this_round = np.zeros_like(reads_seen)      # see below

    cdef np.ndarray[np.uint8_t, ndim=1] dead_seqs # these sequences have been removed just prior in reset_probN
    # dead_seqs = np.array([1 if p is None else 0 for p in em.probN], dtype=np.uint8)  # doesn't compile?
    dead_seqs = np.zeros(len(em.probN), dtype=np.uint8)  # this lookup faster than list lookup in main loop
    for seq_i in range(len(em.probN)):
        if em.probN[seq_i] is None:
            dead_seqs[seq_i] = 1

    # print >> sys.stderr, "DEBUG: t_check 1: %s"%(timedelta(seconds = time() - t_check))

    loop_t0 = time()
    for alignedread_i in range(bamfile_data.shape[0]):
        #   0       1       2      3     4       5   # "efficient indexing only affects certain index operations, namely those with exactly ndim number of typed integer indices”
        # seq_i, read_i, pair_i, rlen, pos, is_reverse = bamfile_data[alignedread_i]
        seq_i = bamfile_data[alignedread_i, 0]

        if dead_seqs[seq_i]:  # probN is None -- just culled.
            continue

        read_i = bamfile_data[alignedread_i, 1]
        pair_i = bamfile_data[alignedread_i, 2]
        rlen = bamfile_data[alignedread_i, 3]
        pos = bamfile_data[alignedread_i, 4]
        is_reverse= bamfile_data[alignedread_i, 5]

        reads_seen_this_round[read_i] = 1
        # find weight
        if initial_iteration:
            weight = priors[seq_i]
        else:
            if reads_seen[read_i] == 0:
                weight = priors[seq_i]
            else: # get weight from posteriors internal data structures
                # here: no bounds checking, try to do as much in cython as possible
                # weight = posteriors[seq_i, read_i]    # simple way to do this, but slower

                lo = posteriors_indptr[seq_i]
                hi   = posteriors_indptr[seq_i+1]
                # binary search (indices are sorted)
                weight = 0.0  # default value if read_i not found
                while lo < hi:
                    mid = (lo+hi)//2
                    midval = posteriors_indices[mid]
                    if midval < read_i:
                        lo = mid + 1
                    elif midval > read_i:
                        hi = mid
                    else: # found it
                        weight = posteriors_data[mid]
                        break
                # assert weight == posteriors[seq_i, read_i], "%s vs %s"%(weight, posteriors[seq_i, read_i])  # DEBUG.  Make sure I'm doing this right, since no error checking otherwise

        probN_single = probN[seq_i]
        
        if is_reverse:
            # manually inline reverse complement loop here.  It's redundant, but it's all in the name of speed.
            for i in range(rlen):  
                ii = rlen-1-i # the index to the base given that this seq is reversed
                for j in range(4):
                    if complement_numeric_base(reads[read_i, pair_i, ii]) == j:   # this is called base, add (1-P) * weight
                        probN_single[pos + i, j] += ( qual2one_minus_p[quals[read_i, pair_i, ii]] * weight )
                    else:                       # this is not the called base, so add P/3 * weight
                        probN_single[pos + i, j] += ( qual2p_div_3[quals[read_i, pair_i, ii]] * weight )

        else: # not reverse
            for i in range(rlen):
                for j in range(4):
                    if reads[read_i, pair_i, i] == j:   # this is called base, add (1-P) * weight
                        probN_single[pos + i, j] += ( qual2one_minus_p[quals[read_i, pair_i, i]] * weight )
                    else:                       # this is not the called base, so add P/3 * weight
                        probN_single[pos + i, j] += ( qual2p_div_3[quals[read_i, pair_i, i]] * weight )

    td = timedelta(seconds = time()-loop_t0)
    # print >> sys.stderr, "DEBUG: iteration %02d _calc_probN time per 1000 bam entries: %.2e seconds"%(em.iteration_i, (((td.microseconds + (td.seconds + td.days * 24 * 3600) * 10**6) / float(10**6)) / alignedread_i) * 1000.0) # DEBUG
    # print >> sys.stderr, "DEBUG: total loop time (%s bam entries): %s"%(alignedread_i, td)

    t_check = time()                    # DEBUG

    # at this point, probN contains sums.
    # Now sum probs for each base, divide cells by this sum for probs.
    # Only divide cells with at least one base mapped, and assign
    # DEFAULT_ERROR to bases with no mappings.
    cdef double base_sum
    cdef np.ndarray[np.uint8_t, ndim=1] seq_bases
    cdef np.ndarray[np.uint8_t, ndim=1] unmapped_bases
    cdef char* seq
    cdef unsigned short reference_base
    cdef double one_minus_p = 1.0 - em.DEFAULT_ERROR
    cdef double p_div_3     = em.DEFAULT_ERROR / 3.0

    sequence_i2sequence_name_array = np.array(em.sequence_i2sequence_name)  # faster slicing. can't cdef string arrays yet.
    
    for seq_i, probN_single in enumerate(em.probN):
        if probN_single is None:  # means this sequence no longer present
            em.unmapped_bases[seq_i] = None
            continue
        # get sequence for nonmapped bases.  Just get for all, even though only need for ones with nonmapped.
        # not sure if this is particularly slow.
        py_seq = em.fastafile.fetch(str(sequence_i2sequence_name_array[seq_i]))
        seq = py_seq  # necessary workaround for strings in Cython
        seq_bases = seq_alpha2int(seq, len(py_seq))

        unmapped_bases = np.zeros(probN_single.shape[0], dtype=np.uint8)
        # sum probs for each base, divide cells by this sum for probs.
        # only divide cells with at least one base mapped.

        for i in range(probN_single.shape[0]):
            base_sum = 0
            for j in range(4):
                base_sum = base_sum + probN_single[i, j]
            # if no mappings, give 1. - em.DEFAULT_ERROR to existing reference base and em.DEFAULT_ERROR / 3 to everyone else.
            if base_sum == 0:
                unmapped_bases[i] = 1
                reference_base = seq_bases[i]
                for j in range(4):
                    if j == reference_base:
                        probN_single[i, j] = one_minus_p
                    else:
                        probN_single[i, j] = p_div_3
            else: # has mappings.  adjust probs.
                for j in range(4):
                    probN_single[i, j] = probN_single[i, j] / base_sum
        em.unmapped_bases[seq_i] = unmapped_bases.copy()

    # print >> sys.stderr, "DEBUG: t_check 2: %s"%(timedelta(seconds = time() - t_check))
    t_check = time()                    # DEBUG

    reads_seen |= reads_seen_this_round  # logical or equals, i.e. update newly seen reads this round

    # print >> sys.stderr, "DEBUG: t_check 3: %s"%(timedelta(seconds = time() - t_check))
    return

# will I want to add @cython.boundscheck(False) ?
@cython.boundscheck(False)
def _calc_posteriors(em):
    """
    Calculates Pr(S|R) for all sequence read pairs
    requires that the likelihood and priors are already calculated.
    """
    # alternately:  use coo format.  pre-populate posterior size same as likelihood size.  pre-populate denom array with size of n_reads.
    # Go through coo lik, putting numerator in posterior and adding to denom sum as you go.
    # lik and post are (seq_n x read_n)
    lik_coo = em.likelihoods.tocoo()
    cdef np.ndarray[np.float_t, ndim=1] lik_data = lik_coo.data
    cdef np.ndarray[np.int32_t, ndim=1] lik_row = lik_coo.row
    cdef np.ndarray[np.int32_t, ndim=1] lik_col = lik_coo.col

    cdef np.ndarray[np.float_t, ndim=1] posteriors_data = np.zeros_like(lik_coo.data)
    cdef np.ndarray[np.float_t, ndim=1] priors = np.array(em.priors[-2], dtype = np.float)    
    cdef np.ndarray[np.float_t, ndim=1] denominator = np.zeros(em.n_reads, dtype=np.float)

    cdef unsigned int i
    cdef float product

    for i in range(posteriors_data.shape[0]):
        # lik_row[i] should be row in posteriors matrix, which should be seq_i
        product = lik_data[i] * priors[lik_row[i]]
        posteriors_data[i] = product
        # lik_col[i] should return read_i
        denominator[lik_col[i]] += product

    # make second pass through posterior coo, dividing data by denom sum for that read.
    for i in range(posteriors_data.shape[0]):
        if denominator[lik_col[i]] == 0:  # e.g. if prior was 0
            # assert posteriors_data[i] == 0
            posteriors_data[i] = 0
        else:
            posteriors_data[i] = posteriors_data[i] / denominator[lik_col[i]]

    # now replace most recent posterior data structure with this built one (in coo format)
    em.posteriors[-1] = sparse.coo_matrix((posteriors_data, (lik_row, lik_col)), shape=em.likelihoods.shape, dtype=em.posteriors[-1].dtype)
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
    cdef int aln_columns    = 0
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


def populate_reads_arrays(em):
    """
    2.  populates:
            em.reads
            em.quals
            em.readlengths
    """
    cdef np.ndarray[np.uint8_t, ndim=3] reads = em.reads  # n_reads X pair X max_read_length
    cdef np.ndarray[np.uint8_t, ndim=3] quals = em.quals
    cdef np.ndarray[np.uint16_t, ndim=2] readlengths = em.readlengths

    cdef unsigned int read_index
    cdef unsigned int readtype_index
    cdef int readlength
    cdef unsigned int i
    cdef int ascii_offset = em.reads_ascii_offset
    
    cdef pykseq.Kseq ks
    cdef pykseq.kseq_t_pointer ks_p

    # first do /1 reads, then do /2 reads if present
    for reads_filepath, readtype_index in [(em.reads1_filepath, 0),
                                           (em.reads2_filepath, 1)]:
        if reads_filepath is None:  # should only happen if no reads2_filepath present
            continue
        read_index = 0
        ks = pykseq.Kseq(reads_filepath.encode('utf8'))
        while 1:
            ks_p = ks.c_read_sequence()
            if ks_p == NULL:
                break
            else:
                # read_index = atoi(ks_p.name.s)
                readlengths[read_index] = ks_p.seq.l
                for i in range(0, ks_p.seq.l):
                    reads[read_index, readtype_index, i] = base_alpha2int(ks_p.seq.s[i])
                    quals[read_index, readtype_index, i] = ks_p.qual.s[i] - ascii_offset
                read_index += 1
        
    return
    
def process_bamfile(em, int ascii_offset):
    """
    read in bamfile and populate:
       bamfile_data    # [seq_i, read_i, pair_i, rlen, pos, is_reverse]
       sequence_name2sequence_i if new entries are here that we haven't seen in previous rounds.

    assume read name is an integer == read index.
    """
    # for samfile use mode "r" (ascii text, currently in use because
    # empirical tests show it is faster when bowtie creates this file
    # not to pipe it to samtools, even though CPU load (always < 100%
    # ) doesn't explain this)
    if em.current_bam_filename.endswith("bam"):
        mode = "rb"
    elif em.current_bam_filename.endswith("sam"):
        mode = "r"
    bamfile = pysam.Samfile(em.current_bam_filename, mode)   
    # TODO: any Cython tricks with string-->int dictionaries (looks like maybe in pxd)
    cdef dict sequence_name2sequence_i = em.sequence_name2sequence_i
    cdef list sequence_i2sequence_name = em.sequence_i2sequence_name

    cdef np.ndarray[np.uint32_t, ndim=2] bamfile_data
    bamfile_data = np.empty((em.n_alignments, 6), dtype=np.uint32)

    cdef int i                          # generic counter

    cdef unsigned int alignedread_i
    cdef unsigned int pos
    cdef unsigned int tid
    # cdef char *refname # ??
    cdef unsigned int read_i
    cdef unsigned int seq_i
    cdef unsigned int new_seq_i
    cdef unsigned int pair_i
    cdef unsigned int rlen
    cdef unsigned int read_i_to_cache
    cdef char *qname
    cdef bint is_reverse
    cdef unsigned int n_reads_mapped = 0  # num. of single reads or pairs
    cdef np.ndarray[np.uint8_t, ndim=1] reads_mapped = np.zeros(em.n_reads, dtype=np.uint8) # mapped this round = 1
    cdef np.ndarray[np.uint32_t, ndim=1] base_coverage # single base coverage
    
    samfile_references = np.array(bamfile.references)  # faster to slice numpy array
    samfile_lengths    = np.array(bamfile.lengths)  # faster to slice numpy array

    # go through header here and create mapping between reference tid
    # (in each alignedread) and my internal seq_i.  If a new seq_i is
    # seen, map name to seq_i here once, rather than doing costly
    # dictionary lookup for each and every alignedread
    cdef np.ndarray[np.uint32_t, ndim=1] tid2seq_i = np.zeros(len(bamfile.references), dtype=np.uint32)
    cdef np.ndarray[np.uint8_t, ndim=1] tid_mapped = np.zeros(len(bamfile.references), dtype=np.uint8)
    t_check = time()
    for alignedread in bamfile:
        tid_mapped[alignedread.tid] = 1
    # print >> sys.stderr, "DEBUG: first pass bamfile time: %s"%(timedelta(seconds = time()-t_check))

    bamfile = pysam.Samfile(em.current_bam_filename, mode)  # reopen file to reset.  seek(0) has problems? 

    new_seq_i = len(sequence_i2sequence_name)
    for tid, refname in enumerate(bamfile.references):
        # assert bamfile.getrname(tid) == refname # this is true
        if tid_mapped[tid] == 1:
            if not sequence_name2sequence_i.has_key(refname):  # new sequence we haven't seen before
                seq_i = new_seq_i
                new_seq_i = new_seq_i + 1
                sequence_name2sequence_i[refname] =  seq_i
                sequence_i2sequence_name.append(refname)
                em.base_coverages.append(np.zeros(samfile_lengths[tid], dtype=np.uint32))
            else:
                seq_i = sequence_name2sequence_i[refname]
            tid2seq_i[tid] = seq_i
        
    # reset this here to be size of n_sequences
    em.coverage = np.zeros(len(sequence_name2sequence_i), dtype=np.uint32)
    cdef np.ndarray[np.uint32_t, ndim=1] coverage = em.coverage

    # and now keep temporary base_coverages in numpy array for speed below.
    cdef list base_coverages = em.base_coverages  # list of arrays
    max_length = max([base_coverage.shape[0] for base_coverage in base_coverages])
    cdef np.ndarray[np.uint32_t, ndim=2] base_coverages_2d = np.zeros((len(base_coverages), max_length), dtype=np.uint32)
    
    for seq_i, base_coverage in enumerate(base_coverages):
        for i in range(base_coverage.shape[0]):
            base_coverages_2d[seq_i, i] = base_coverage[i]

    alignedread_i = 0
    for alignedread in bamfile:
        #TODO: decide if it's better to assign first to typed variable, or assign directly to bamfile_data
        # tid = alignedread.tid
        # pair_i = alignedread.is_read2
        qname_temp = alignedread.qname
        qname = qname_temp      # this extra assignment is necessary to go from python object to char*. (see Cython Language Basics)
        # is_reverse = alignedread.is_reverse
        # rlen = alignedread.rlen

        read_i = atoi(qname)
        # read_i = int(qname_temp)
        reads_mapped[read_i] = 1
        seq_i = tid2seq_i[alignedread.tid]
        coverage[seq_i] += <unsigned int>alignedread.rlen
        # base_coverage = em.base_coverages[seq_i]  # list lookup.  SLOW?
        # for i in range(alignedread.pos, alignedread.pos+alignedread.rlen):
        #     base_coverage[i] += 1
        for i in range(alignedread.pos, alignedread.pos+alignedread.rlen):
             base_coverages_2d[seq_i, i] += 1    # base_coverage[i] += 1
        
            
        bamfile_data[alignedread_i, 0] = seq_i
        bamfile_data[alignedread_i, 1] = read_i
        bamfile_data[alignedread_i, 2] = alignedread.is_read2
        bamfile_data[alignedread_i, 3] = alignedread.rlen
        bamfile_data[alignedread_i, 4] = alignedread.pos
        bamfile_data[alignedread_i, 5] = alignedread.is_reverse

        alignedread_i += 1
    assert alignedread_i == bamfile_data.shape[0]
    em.bamfile_data = bamfile_data.copy()
    bamfile.close()

    # get base_coverages back in list:
    for seq_i, base_coverage in enumerate(base_coverages):
        for i in range(base_coverage.shape[0]):
            base_coverage[i] = base_coverages_2d[seq_i, i]
    
    cdef int n_reads_total = reads_mapped.shape[0]
    for i in range(n_reads_total):
        n_reads_mapped = n_reads_mapped + reads_mapped[i]
    em.n_reads_mapped = n_reads_mapped
    return

def reset_probN(em):
    """
    sets probN back to zeros for all sequences
    gets coverage set for all sequences by dividing current values by sequence lengths
    does culling based on base_coverages
    
    """
    if em.current_bam_filename.endswith("bam"):
        mode = "rb"
    elif em.current_bam_filename.endswith("sam"):
        mode = "r"
    bamfile = pysam.Samfile(em.current_bam_filename, mode)
    
    cdef int i
    cdef int l
    cdef unsigned int cov_thresh    = em.min_length_coverage_def  # above what coverage
    cdef double length_thresh       = em.min_length_coverage      # percent length required
    cdef double percent_length_covered
    cdef int cullcount = 0
    
    references_array = np.array(bamfile.references) # slicing these tuples is stupid-slow, for some reason.
    references_lengths = np.array(bamfile.lengths)
    for i in range(len(bamfile.lengths)):
        seq_i = em.sequence_name2sequence_i.get(references_array[i])
        if seq_i is not None:  # since not all seqs in header will be ones with mapped reads
            # CULLING: check for coverage > 0 across a % threshold of bases
            percent_length_covered = float((em.base_coverages[seq_i] >= cov_thresh).sum()) / em.base_coverages[seq_i].shape[0]
            if percent_length_covered < length_thresh:
                em.probN[seq_i] = None
                em.coverage[seq_i] = 0
                cullcount += 1
                continue
            # else is a valid seq, so create empty probN matrix
            l = references_lengths[i]
            em.probN[seq_i] = np.zeros((l, 5), dtype=np.float)   #ATCG[other] --> 01234
            em.coverage[seq_i] = em.coverage[seq_i] / float(l)
    bamfile.close()

    # EXPERIMENTAL: if initial iteration, then redistribute Priors out over only sequences that weren't culled.
    # if em.iteration_i < 1:  # initial_iteration
        

    if em._VERBOSE:
        if cullcount > 0:
            print >> sys.stderr, "Culled %s sequences in iteration %02d due to low fraction of reference sequence bases covered by >= %s reads"%(cullcount, em.iteration_i, cov_thresh)

    return

@cython.boundscheck(False) # turn off bounds-checking for entire function
cdef np.ndarray[np.uint8_t, ndim=1] seq_alpha2int(char* seq, int seqlen):
    cdef unsigned int i
    np_as_ints = np.empty(seqlen, dtype=np.uint8)
    for i in range(seqlen):
        np_as_ints[i] = base_alpha2int(seq[i])
    return np_as_ints

