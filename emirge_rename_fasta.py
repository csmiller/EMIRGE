#!/usr/bin/env python
"""
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

for help, type:
emirge_rename_fasta.py --help

re-write the fasta file to stdout with renamed sequences:
>SEQ_ID|ORIG_SEQ_NAME prior_prob prior/length

sort by decending prior probability
bases with no read support are labeled 'N', and terminal N's are trimmed
"""

import sys
import os
from optparse import OptionParser
from Bio import SeqIO, Seq
import copy
from math import log, e
import cPickle
from gzip import GzipFile
import numpy

# from emirge.py: self.DEFAULT_ERROR = 0.05
DEFAULT_ERROR = 0.05

def replace_with_Ns(probN, seq_i, seq):
    """
    IN:  probN matrix, seq_i for sequence, seq (BioPython Seq object) 
    OUT: returns the sequence where bases with no read support are replaced with "N", and terminal N's are trimmed
    """
    default_probN = 1.0 - DEFAULT_ERROR

    try:
        this_probN = probN[seq_i]
    except:
        print >> sys.stderr, seq_i, len(probN)
        raise
        
    indices = numpy.where(numpy.max(this_probN, axis=1) == default_probN)
    newseq = numpy.array(str(seq), dtype='c')
    newseq[indices] = 'N'
    newseq = ''.join(newseq)
    newseq = newseq.strip("N")
    return Seq.Seq(newseq)

def rename(wd = os.getcwd(), prob_min = None, record_prefix = '', no_N = False):
    """
    wd is an iteration directory
    prob_min is minimum prior prob of a sequence to include in output.  If None, include all seqs.
    """
    if prob_min is None:
        prob_min = -1
    current_iter = int(wd.split('.')[-1])
    prior_file = file(os.path.join(wd, 'priors.iter.%02d.txt'%current_iter))
    probN_filename = os.path.join(wd, 'probN.pkl.gz')
    probN = cPickle.load(GzipFile(probN_filename))

    name2seq_i = {}
    name2prior = {}
    name2normed_prior = {}  # normed by length of sequences
    for line in prior_file:
        atoms = line.split()
        name2seq_i[atoms[1]] = int(atoms[0])
        name2prior[atoms[1]] = atoms[2]

    sorted_records = []
    for record in SeqIO.parse(file(os.path.join(wd, "iter.%02d.cons.fasta"%current_iter)), "fasta"):
        name = record.description.split()[0]
        record.id = "%s%d|%s"%(record_prefix, name2seq_i[name], name)
        if not no_N:
            record.seq = replace_with_Ns(probN, name2seq_i[name], record.seq)
        record.description = ""
        sorted_records.append((float(name2prior[name]), record)) 

    # normalize priors by length
    sum_lengths = sum(len(record.seq) for prior, record in sorted_records)
    normed_priors = [prior/ len(record.seq) for prior, record in sorted_records]
    sum_norm = float(sum(normed_priors))
    normed_priors = [x/sum_norm for x in normed_priors]
    for i, (prior, record) in enumerate(sorted_records):
        record.description = "Prior=%06f Length=%d NormPrior=%06f"%(prior, len(record.seq), normed_priors[i])

    for prior, record in sorted(sorted_records, reverse=True):
        if prior < prob_min:
            break
        sys.stdout.write(record.format('fasta'))

def main(argv = None):
    if argv is None:
        argv = sys.argv[1:]

    parser = OptionParser("usage: %prog [options] <iter.DIR>  > renamed.fasta\n\nRewrites emirge fasta file to include proper sequence names and prior probabilities (abundance estimates) in record headers.\nIf no iter.DIR is given, assumes cwd.\nBases with no read support are labeled 'N', and terminal N's are trimmed")
    parser.add_option('-p', '--prob_min',
                      type='float',
                      help='Only include sequences in output with prior probability above PROB_MIN (Default: include all sequences)')
    parser.add_option('-r', '--record_prefix',
                      type='string', default = '',
                      help='Add the specified prefix to each fasta record title')
    parser.add_option('-n', '--no_N',
                      action="store_true", 
                      help="Don't change bases with no read support to N.  Caution: these bases are not supported by reads in the input data, but will usually be from a closely related sequence.")

    
    (options, args) = parser.parse_args(argv)
    if len(args) == 0:
        wd = os.getcwd()
    elif len(args) == 1:
        wd = os.path.abspath(args[0])
    else:
        parser.error("Found more than one argument on command line (expects only iter.DIR): %s"%(args))
    assert os.path.exists(wd), "Directory not found: %s"%(wd)
    rename(wd, options.prob_min, options.record_prefix, options.no_N)
    


if __name__ == '__main__':
    main()
