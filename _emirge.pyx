"""
_emirge.pyx contains helper cython functions for emirge.py

EMIRGE: Expectation-Maximization Iterative Reconstruction of Genes from the Environment
Copyright (C) 2010 Christopher S. Miller  (csmiller@berkeley.edu)

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

from libc.stdlib cimport malloc, free
from libc.math cimport M_E, log, pow as c_pow
import numpy as np
cimport numpy as np
import sys

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

def calc_probN_read(bint initial_iteration,
                    unsigned int seq_i, unsigned int read_i,
                    unsigned int pos,
                    double weight,
                    # np.ndarray[np.float_t, ndim=1] priors,
                    # np.ndarray[np.float_t, ndim=2] posteriors,
                    np.ndarray[np.uint8_t, ndim=1] numeric_bases,
                    np.ndarray[np.uint8_t, ndim=1] qualints,
                    np.ndarray[np.float_t, ndim=2] probN):  # an individual probN numpy matrix. Passed by ref?
    """
    calculates effect of a single read on probN

    this is inside loop over bamfile_data
    """
    cdef int i
    cdef int j
    # cdef double weight

    ## posteriors is a sparse matrix, so this presents a little bit of a tricky situation until
    ## I can think about how to pass this.  Have to choose correct data struct, manip internals
    
    # if initial_iteration:
    #     weight = priors[seq_i]
    # else:
    #     if seq_i < posteriors.shape[0] and read_i < posteriors.shape[1]:
    #         weight = posteriors[seq_i, read_i]
    #     else:
    #         weight = priors[seq_i]

    for i in range(numeric_bases.shape[0]):
    # for i from 0 <= i < numeric_bases.shape[0]:  # old syntax, range syntax above fast if i is cdef-ed.
        s = 0.0
        for j in range(4):
            if numeric_bases[i] == j:   # this is called base, add (1-P) * weight
                probN[pos + i, j] += qual2one_minus_p[qualints[i]] * weight
            else:                       # this is not the called base, so add P/3 * weight
                probN[pos + i, j] += qual2p_div_3[qualints[i]] * weight

    return

