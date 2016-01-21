"""
_emirge_amplicon.pyx contains helper cython functions for emirge_amplicon.py

EMIRGE: Expectation-Maximization Iterative Reconstruction of Genes from the Environment
Copyright (C) 2010-2012 Christopher S. Miller  (christopher.s.miller@ucdenver.edu)

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
from _kseq cimport Kseq
from Emirge import log as logger


# for lookup table of qual values
from Emirge.log import INFO
cdef extern from *:
    ctypedef double* static_const_double_ptr "static const double*"


cdef extern from "tables.h":
    static_const_double_ptr qual2one_minus_p   # lookup tables to avoid too many double calcs.
    static_const_double_ptr qual2p_div_3


cpdef inline int base_alpha2int(char base_ascii):
    """
    Convert ASCII encoded base to numeric base.

    Equivalent to lookup in dict: base2i = {"A":0,"T":1,"C":2,"G":3,"N":4}
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


@cython.boundscheck(False) # turn off bounds-checking for entire function
cpdef np.ndarray[np.uint8_t, ndim=1] seq_alpha2int(char* seq, int seqlen):
    cdef unsigned int i
    np_as_ints = np.empty(seqlen, dtype=np.uint8)
    for i in range(seqlen):
        np_as_ints[i] = base_alpha2int(seq[i])
    return np_as_ints


cpdef inline unsigned char complement_numeric_base(unsigned char c):
    return (c ^ 1) ^ (c >> 2)
    # Explanation:
    # (c ^ 1) because:
    # A = 000  <=>  T = 001
    # C = 010  <=>  G = 011  ==>  xor c with 1 to complement
    # ... ^ (c >> 2) because:
    # N = 100  <=>  N = 100, but c ^ 1 gives us 101
    # ==> xor with c >> 2, which is 1 only for N (=100) to unset last bit


@cython.boundscheck(False)
def calc_likelihood(em):
    """
    Cython helper function with typed data structures, fast looping.
    """
    cdef np.ndarray[np.uint32_t, ndim=2] bamfile_data = em.bamfile_data
    cdef np.ndarray[np.uint8_t, ndim=3] reads = em.reads
    cdef np.ndarray[np.uint8_t, ndim=3] quals = em.quals                    
    cdef list probN = em.probN
    cdef list prob_indels = em.prob_indels
    cdef list cigars = em.cigars
    

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
    cdef unsigned int ref_i
    cdef unsigned int rlen
    cdef unsigned int pos
    cdef bint is_reverse

    cdef double p
    cdef float s
    cdef double result

    cdef unsigned int i
    cdef unsigned int ii
    cdef unsigned char j
    cdef unsigned int k
    cdef unsigned int matchlen
    cdef unsigned int matchtype

    # keep track of whether we have collected 2 of 2 pairs or not.
    # start at 1 because of way I flip this flag every iteration
    # below. (a bit counter-intuitive, but fastest thing I could think
    # of)  Assumption is that every two reads is a mapped pair spit out
    # from bowtie in bamfile
    cdef bint pair_collected = 1        
    cdef unsigned int result_i       = 0
    cdef np.ndarray[np.float_t, ndim=2] probN_single  # an individual probN numpy matrix
    cdef np.ndarray[np.float_t, ndim=2] prob_indels_single
    cdef np.ndarray[np.uint8_t, ndim=1] numeric_bases_single
    cdef np.ndarray[np.uint8_t, ndim=1] quals_single

    cdef np.ndarray[np.uint8_t, ndim=1] dead_seqs # these sequences have been removed previously
    dead_seqs = np.zeros(len(em.probN), dtype=np.uint8)
    for seq_i in range(len(em.probN)):
        if em.probN[seq_i] is None:
            dead_seqs[seq_i] = 1

    p = 0.0  # start off with no pair collected
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

            numeric_bases_single = reads[read_i,pair_i] # and why are we calling this numeric bases better as reads_single?
            quals_single = quals[read_i,pair_i]
            probN_single = probN[seq_i]
            prob_indels_single   = prob_indels[seq_i]
            cigarlist = cigars[alignedread_i] #Gives list representing cigar string for alignedread_i
        
            ref_i = 0
            if is_reverse:
                # put reverse case here, see below:
                #Iterate through the alignedread cigars, and calc likelihood for each read
                i = rlen-1
                # for each tuple in cigar string in alignedread_i, representing
                # a mapping in the cigar string and length
                for matchtype, matchlen in cigarlist:
                    s = 0.0  #here?
                    if matchtype == 0:  # Cigar String - Match
                        for k in range(matchlen):
                            # increment i, and ref_i once for each base,
                            # and calc prob(n) for each base
                            for j in range(4):
                                if complement_numeric_base(numeric_bases_single[i]) == j:   # this is called base, set 1-P
                                    # s += (1. - error_p_i) * probN[pos + i, j]
                                    s += ( qual2one_minus_p[quals_single[i]] * probN_single[pos + ref_i, j] )    # lookup table
                                else:                       # this is not the called base, so set to P/3
                                    # s += (error_p_i / 3 ) * probN[pos + i, j]
                                    s += ( qual2p_div_3[quals_single[i]] * probN_single[pos + ref_i, j] )
                            ref_i += 1
                            i -= 1
                            
                    elif matchtype == 1:  # Cigar String - Insertion to reference
                        for k in range(matchlen): #increment i for each inserted base
                            i -= 1
                            
                    elif matchtype == 2:  # Cigar String - Deletion from the reference
                        for k in range(matchlen): #increment ref_i for each deleted base
                            ref_i += 1   
                    else:
                        logger.critical("Calc Likelihood Failure: Invalid Cigar "
                                        "String with header value of %i" % (matchlen))
                        sys.exit("Calc Likelihood Failure: Invalid Cigar String Header Value")
                        
                    p += log(s)
            
            
            else:  #not reverse    
                #Iterate through the alignedread cigars, and calc likelihood for each read
                i = 0 # i is the distance from the pos of a mapping in a read.  Doesn't increment for deletions in the read, but does for insertions. In other words, i is the base(n) of a read sequence
                # for each tuple in cigar string in alignedread_i, representing
                # a mapping in the cigar string and length
                for matchtype, matchlen in cigarlist:
                    s = 0.0  #here?
                    if matchtype == 0:  # Cigar String - Match 
                        for k in range(matchlen):
                            # increment i, and ref_i once for each base,
                            # and calc prob(n) for each base
                            for j in range(4):
                                if numeric_bases_single[i] == j:   # this is called base, set 1-P
                                    # s += (1. - error_p_i) * probN[pos + i, j]
                                    s += ( qual2one_minus_p[quals_single[i]] * probN_single[pos + ref_i, j] )    # lookup table
                                else:                       # this is not the called base, so set to P/3
                                    # s += (error_p_i / 3 ) * probN[pos + i, j]
                                    s += ( qual2p_div_3[quals_single[i]] * probN_single[pos + ref_i, j] )
                            ref_i += 1
                            i += 1
                            
                    elif matchtype == 1:  # Cigar String - Insertion to reference
                        for k in range(matchlen): #increment i for each inserted base, and calc prob(n) for each base
                            i += 1
                    elif matchtype == 2:  # Cigar String - Deletion from the reference
                        for k in range(matchlen): #increment ref_i for each deleted base
                            ref_i += 1
                    else:
                        print >> sys.stderr, "Calc Likelihood Failure: Invalid Cigar String with header value of", matchlen
                        sys.exit("Calc Likelihood Failure: Invalid Cigar String Header Value")
                    
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
    
    em.likelihoods = sparse.coo_matrix((lik_data, (lik_row_seqi, lik_col_readi)), em.likelihoods.shape, dtype=em.likelihoods.dtype).tocsr()

@cython.boundscheck(False)
def calc_probN(em):
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
        em.posteriors[-2] = em.posteriors[-2].tocsr()
        em.posteriors[-2].sort_indices()
        posteriors = em.posteriors[-2]  # this depends on PREVIOUS iteration's posteriors (seq_n x read_n)
        posteriors_data = posteriors.data
        posteriors_indices = posteriors.indices
        posteriors_indptr = posteriors.indptr
    else:
        posteriors = None

    cdef list probN      = em.probN     # list of numpy arrays.  At this point, only None values will be due to just-culled sequences.
    cdef list cigars = em.cigars        # should be list of lists of pysam cigartuples read from bamfile in process_bamfile()
    cdef list prob_indels = em.prob_indels  
    
    # for bamfile_data
    cdef unsigned int alignedread_i
    cdef unsigned int seq_i
    cdef unsigned int read_i
    cdef unsigned int pair_i
    cdef unsigned int ref_i
    cdef unsigned int rlen
    cdef unsigned int pos
    cdef bint is_reverse

    cdef unsigned int i
    cdef unsigned int ii
    cdef unsigned int j
    cdef unsigned int k
    cdef int start # for indexing sparse matrix data structures
    cdef int end
    cdef int lo  # binary search in sparse matrix data struct
    cdef int mid
    cdef int hi
    cdef int midval

    cdef int matchlen
    cdef int matchtype
    
    cdef double weight

    cdef np.ndarray[np.float_t, ndim=2] probN_single        # an individual probN numpy matrix.
    cdef np.ndarray[np.float_t, ndim=2] prob_indels_single  # indel data structure.
    cdef np.ndarray[np.uint8_t, ndim=1] numeric_bases_single
    cdef np.ndarray[np.uint8_t, ndim=1] quals_single

    cdef np.ndarray[np.uint8_t, ndim=1] reads_seen = em.reads_seen      # ...in previous rounds
    cdef np.ndarray[np.uint8_t, ndim=1] reads_seen_this_round = np.zeros_like(reads_seen)      # see below

    cdef np.ndarray[np.uint8_t, ndim=1] dead_seqs # these sequences have been removed just prior in reset_probN
    # dead_seqs = np.array([1 if p is None else 0 for p in em.probN], dtype=np.uint8)  # doesn't compile?
    dead_seqs = np.zeros(len(em.probN), dtype=np.uint8)  # this lookup faster than list lookup in main loop
    for seq_i in range(len(em.probN)):
        if em.probN[seq_i] is None:
            dead_seqs[seq_i] = 1

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

        numeric_bases_single = reads[read_i, pair_i]
        quals_single = quals[read_i, pair_i]
        probN_single         = probN[seq_i]
        prob_indels_single   = prob_indels[seq_i]
        cigarlist = cigars[alignedread_i] #Gives list of cigartuples representing cigar string for alignedread_i

        ref_i = 0 # ref_i is the distance from the pos of a mapping in the reference.  Initialized to start at 0

        if is_reverse:
            i = rlen-1 # index into read, start at right end or read as stored
            # for each tuple in cigar string in alignedread_i, representing
            # a mapping in the cigar string and length
            for matchtype, matchlen in cigarlist:
                if matchtype == 0:  # Cigar String - Match
                    for k in range(matchlen):
                        # increment i, and ref_i once for each base,
                        # and calc prob(n) for each base
                        for j in range(4):
                            if complement_numeric_base(numeric_bases_single[i]) == j:   # this is called base, add (1-P) * weight
                                probN_single[pos + ref_i, j] += qual2one_minus_p[quals_single[i]] * weight
                            else:                       # this is not the called base, so add P/3 * weight
                                probN_single[pos + ref_i, j] += qual2p_div_3[quals_single[i]] * weight
                        prob_indels_single[pos + ref_i, 0] += weight   #Add indel weight for a match
                        ref_i += 1
                        i -= 1
                elif matchtype == 1:  # Cigar String - Insertion to reference
                    # (insertion(s) in read relative to reference when aligned)
                    #Add weight for an insertion after pos+ref_i in reference sequence, subtract 1 since ref_i was already incremented past
                    prob_indels_single[pos + ref_i - 1 , 1] += weight
                    if matchlen > prob_indels_single[pos + ref_i - 1 , 3]:
                         prob_indels_single[pos + ref_i - 1 , 3] = matchlen
                    i -= matchlen
                elif matchtype == 2:  # Cigar String - Deletion from the reference
                    # (gaps in read relative to reference when aligned)
                    for k in range(matchlen):
                        prob_indels_single[pos + ref_i, 2] += weight #Add weight for a deletion at pos+ref_i in reference sequence
                        probN_single[pos+ref_i,4]+= weight  #Put the same weight in the probN matrix for evalution against the other bases
                        ref_i += 1
                else:
                    logger.critical("ProbN Failure: Invalid Cigar String with "
                                    "header value of %i" % (matchlen))
                    sys.exit("Calc Likelihood Failure: Invalid Cigar String Header Value")
            
        else: # not reverse:
            i = 0     # index into read
            # Iterate through the alignedread cigars, and calc prob(n) for each position
            # for each tuple in cigar string in alignedread_i, representing
            # a mapping in the cigar string and length
            for matchtype, matchlen in cigarlist:
                if matchtype == 0:  # Cigar String - Match 
                    for k in range(matchlen):
                        # increment i, and ref_i once for each base,
                        # and calc prob(n) for each base
                        for j in range(4):
                            if numeric_bases_single[i] == j:   # this is called base, add (1-P) * weight
                                probN_single[pos + ref_i, j] += qual2one_minus_p[quals_single[i]] * weight
                            else:                       # this is not the called base, so add P/3 * weight
                                probN_single[pos + ref_i, j] += qual2p_div_3[quals_single[i]] * weight
                        prob_indels_single[pos + ref_i, 0] += weight   #Add indel weight for a match
                        ref_i += 1
                        i += 1
                elif matchtype == 1:  # Cigar String - Insertion to the reference
                    # (i.e. need to insert gap(s) in reference seq to align)
                    # Add weight for an insertion after pos+ref_i in reference sequence,
                    # subtract 1 since ref_i was already incremented past
                    prob_indels_single[pos + ref_i - 1, 1] += weight
                    if matchlen > prob_indels_single[pos + ref_i - 1 , 3]:
                         prob_indels_single[pos + ref_i - 1 , 3] = matchlen
                    i += matchlen    # move index in read 
                elif matchtype == 2:  # Cigar String - Deletion from the reference
                    # (i.e. need to insert gap(s) in read to align)
                    for k in range(matchlen):
                        prob_indels_single[pos + ref_i, 2] += weight #Add weight for a deletion at pos+ref_i in reference sequence
                        probN_single[pos+ref_i,4]+= weight  #Put the same weight in the probN matrix for evalution against the other bases
                        ref_i += 1
                else:
                    logger.critical("ProbN Failure: Invalid Cigar String with "
                                    "header value of %i" % (matchlen))
                    sys.exit("Calc Likelihood Failure: Invalid Cigar String Header Value")

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

    # faster slicing. can't cdef string arrays yet.
    sequence_i2sequence_name_array = np.array(em.sequence_i2sequence_name)

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
            for j in range(5): #upping this to 5 from 4 because we should consider weight for a deletion as a read mapping
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
                for j in range(5): #as above, increased this to include deletion mappings
                    probN_single[i, j] = probN_single[i, j] / base_sum
        em.unmapped_bases[seq_i] = unmapped_bases.copy()

    reads_seen |= reads_seen_this_round  # logical or equals, i.e. update newly seen reads this round

    return

# will I want to add @cython.boundscheck(False) ?
@cython.boundscheck(False)
def calc_posteriors(em):
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

    cdef int matches = 0
    cdef int aln_columns    = 0
    cdef int count
    cdef int i
    cdef int j

    for i in range(len(alncode_list)):
        count_a, aln_code = alncode_list[i]
        if not len(count_a):
            count = 1  # "As a special case, if n=1 then n is omitted"
        else:
            count = int(count_a)

        if aln_code == 'M':    # match (letter-letter column; could be mismatch)
            for j in range(count):
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

@logger.timed("Preallocating reads and quals in memory")
def populate_reads_arrays(em):
    """
    2.  populates:
            em.reads
            em.quals
            em.readlengths
    """
    # n_reads X pair X max_read_length
    cdef np.ndarray[np.uint8_t, ndim=3] reads = em.reads
    cdef np.ndarray[np.uint8_t, ndim=3] quals = em.quals
    cdef np.ndarray[np.uint16_t, ndim=2] readlengths = em.readlengths

    cdef unsigned int read_index
    cdef unsigned int readtype_index
    cdef int readlength
    cdef unsigned int i
    cdef int ascii_offset = em.reads_ascii_offset
    
    cdef Kseq ks

    # first do /1 reads, then do /2 reads if present
    for reads_filepath, readtype_index in [(em.reads1_filepath, 0),
                                           (em.reads2_filepath, 1)]:
        if reads_filepath is None:
            # should only happen if no reads2_filepath present
            continue
        read_file = Kseq(reads_filepath).open()
        read_index = 0
        while not read_file.read_next_sequence() < 0:
            readlengths[read_index] = read_file.ks.seq.l
            for i in range(0, read_file.ks.seq.l):
                reads[read_index, readtype_index, i] = \
                    base_alpha2int(read_file.ks.seq.s[i])
                quals[read_index, readtype_index, i] = \
                    read_file.ks.qual.s[i] - ascii_offset
            read_index += 1
        read_file.close()


def process_bamfile(em, bamfile, int ascii_offset):
    """
    read in bamfile and populate:
       bamfile_data    # [seq_i, read_i, pair_i, rlen, pos, is_reverse]
       sequence_name2sequence_i if new entries are here that we haven't seen in previous rounds.
       cigars  # list of cigar tuples per alignment.  Same length as bamfile_data

    assume read name is an integer == read index.
    """
    # TODO: any Cython tricks with string-->int dictionaries (looks like maybe in pxd)
    cdef dict sequence_name2sequence_i = em.sequence_name2sequence_i
    cdef list sequence_i2sequence_name = em.sequence_i2sequence_name
    cdef list refseq_lengths = em.refseq_lengths

    cdef np.ndarray[np.uint32_t, ndim=2] bamfile_data
    # TODO: replace em.n_alignments with bamfile.mapped
    bamfile_data = np.empty((em.n_alignments, 6), dtype=np.uint32)
    cdef list cigars = []  # it's an empty list stores cigartuples from pysam alignments
    cigars = range(bamfile_data.shape[0])  # make same size as bamfile_data
    cdef list base_coverages = em.base_coverages

    cdef int i                          # generic counter

    cdef unsigned int alignedread_i
    cdef unsigned int pos
    cdef unsigned int tid
    cdef unsigned int read_i
    cdef unsigned int seq_i
    cdef unsigned int new_seq_i
    cdef unsigned int pair_i
    cdef unsigned int rlen
    cdef unsigned int read_i_to_cache
    cdef unsigned int matchtype
    cdef unsigned int matchlen
    cdef char *qname
    cdef bint is_reverse
    cdef unsigned int n_reads_mapped = 0  # num. of single reads or pairs
    cdef np.ndarray[np.uint8_t, ndim=1] reads_mapped = np.zeros(em.n_reads, dtype=np.uint8) # mapped this round = 1
    cdef np.ndarray[np.uint32_t, ndim=1] base_coverage # single base coverage
    cdef np.ndarray[np.uint32_t, ndim=1] tid2seq_i
    cdef np.ndarray[np.uint8_t, ndim=1] tid_mapped

    with bamfile as bam:
        # copy some data to faster numpy arrays
        samfile_references = np.array(bam.references)
        samfile_lengths    = np.array(bam.lengths)
        max_length = max(bam.lengths)

        # go through header here and create mapping between reference tid
        # (in each alignedread) and my internal seq_i.  If a new seq_i is
        # seen, map name to seq_i here once, rather than doing costly
        # dictionary lookup for each and every alignedread
        tid2seq_i = np.zeros(len(bam.references), dtype=np.uint32)
        tid_mapped = np.zeros(len(bam.references), dtype=np.uint8)

        for alignedread in bam:
            tid_mapped[alignedread.tid] = 1

    new_seq_i = len(sequence_i2sequence_name)
    for tid, refname in enumerate(samfile_references):
        if tid_mapped[tid] != 1:
            continue

        if not sequence_name2sequence_i.has_key(refname):
            # new sequence we haven't seen before
            seq_i = new_seq_i
            new_seq_i += 1
            sequence_name2sequence_i[refname] = seq_i
            sequence_i2sequence_name.append(refname)
            em.base_coverages.append(np.zeros(samfile_lengths[tid], dtype=np.uint32))
            refseq_lengths.append(samfile_lengths[tid])
        else:
            seq_i = sequence_name2sequence_i[refname]

        tid2seq_i[tid] = seq_i
        refseq_lengths[seq_i] = samfile_lengths[tid]

    # and now keep temporary base_coverages in numpy array for speed below.
    cdef np.ndarray[np.uint32_t, ndim=2] base_coverages_2d
    base_coverages_2d = np.zeros((len(base_coverages), max_length),
                                 dtype=np.uint32)

    for seq_i, base_coverage in enumerate(base_coverages):
        em.base_coverages[seq_i] = np.zeros(max_length,dtype=np.uint32)

    with bamfile as bam: # re-open to start from beginning
        # Now actually iterate through alignments
        alignedread_i = 0
        for alignedread in bam:
            #TODO: decide if it's better to assign first to typed variable,
            # or assign directly to bamfile_data
            # tid = alignedread.tid
            # pair_i = alignedread.is_read2
            qname_temp = alignedread.qname
            qname = qname_temp      # this extra assignment is necessary to
            # go from python object to char*. (see Cython Language Basics)
            # is_reverse = alignedread.is_reverse
            # rlen = alignedread.rlen

            read_i = atoi(qname)
            # read_i = int(qname_temp)
            reads_mapped[read_i] = 1
            seq_i = tid2seq_i[alignedread.tid]
            # indexed by alignedread_i just as bamfile_data
            cigars[alignedread_i] = alignedread.cigartuples

            i = alignedread.pos
            for matchtype, matchlen in cigars[alignedread_i]:
                if matchtype == 0:  # match
                    for j in range(matchlen):
                        base_coverages_2d[seq_i, i] += 1
                        i+=1
                elif matchtype == 1:  # Cigar String - Insertion wrt reference
                    # (i.e. need to insert gap(s) in reference seq to align)
                    for j in range(matchlen):
                        i+=0
                elif matchtype == 2:  # Cigar String - Deletion wrt reference
                    # (i.e. need to insert gap(s) in read to align)
                    for j in range(matchlen):
                        i+=1

            bamfile_data[alignedread_i, 0] = seq_i
            bamfile_data[alignedread_i, 1] = read_i
            bamfile_data[alignedread_i, 2] = alignedread.is_read2
            bamfile_data[alignedread_i, 3] = alignedread.rlen
            bamfile_data[alignedread_i, 4] = alignedread.pos
            bamfile_data[alignedread_i, 5] = alignedread.is_reverse

            alignedread_i += 1

        assert alignedread_i == bamfile_data.shape[0]
        em.bamfile_data = bamfile_data.copy()
        em.cigars = cigars[:]  # is slicing necessary for safe copy?

    # get base_coverages back in list:
    for seq_i, base_coverage in enumerate(base_coverages):
        base_coverage[:] = base_coverages_2d[seq_i][:base_coverage.shape[0]]
    em.base_coverages = base_coverages[:]
    
    cdef int n_reads_total = reads_mapped.shape[0]
    for i in range(n_reads_total):
        n_reads_mapped = n_reads_mapped + reads_mapped[i]
    em.n_reads_mapped = n_reads_mapped


def reset_probN(em, bamfile):
    """
    sets probN back to zeros for all sequences
    NEW: also resets prob_indels to zeros - this also updates array size for new seq length 
    gets coverage set for all sequences by dividing current values by sequence lengths
    does culling based on base_coverages
    
    """
    cdef int i
    cdef int l
    cdef unsigned int cov_thresh    = em.min_length_coverage_def  # above what coverage
    cdef double length_thresh       = em.min_length_coverage      # percent length required
    cdef double percent_length_covered
    cdef int cullcount = 0

    with bamfile as bam:
        references_array = np.array(bam.references) # slicing these tuples is stupid-slow, for some reason.
        references_lengths = np.array(bam.lengths)
        for i in range(len(bam.lengths)):
            seq_i = em.sequence_name2sequence_i.get(references_array[i])
            if seq_i is not None:  # since not all seqs in header will be ones with mapped reads
                # CULLING: check for coverage > 0 across a % threshold of bases
                #percent_length_covered = float((em.base_coverages[seq_i] >= cov_thresh).sum()) / em.base_coverages[seq_i].shape[0]
                percent_length_covered = float((em.base_coverages[seq_i] >= cov_thresh).sum()) / em.refseq_lengths[seq_i]
                if percent_length_covered < length_thresh:
                    em.probN[seq_i] = None
                    em.prob_indels[seq_i] = None
                    #em.coverage[seq_i] = 0
                    cullcount += 1
                    continue
                # else is a valid seq, so create empty probN matrix
                l = references_lengths[i]
                em.probN[seq_i] = np.zeros((l, 5), dtype=np.float)   #ATCG[other] --> 01234
                em.prob_indels[seq_i] = np.zeros((l, 4), dtype=np.float)  #0 = match weight, 1 = insertion weight, 2 = deletion weight, 3= insertion length
                #em.coverage[seq_i] = em.coverage[seq_i] / float(l)  # Don't think we need this anymore

    if em._VERBOSE:
        if cullcount > 0:
            print >> sys.stderr, "Culled %s sequences in iteration %02d due to low fraction of reference sequence bases covered by >= %s reads"%(cullcount, em.iteration_i, cov_thresh)

    return
