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
"""

import sys
import os
import re
from optparse import OptionParser
from Bio import SeqIO
import copy
from math import log, e


def rename(wd = os.getcwd(), prob_min = None, record_prefix = ''):
    """
    wd is an iteration directory
    prob_min is minimum prior prob of a sequence to include in output.  If None, include all seqs.
    """
    if prob_min is None:
        prob_min = -1
    current_iter = int(wd.split('.')[-1])
    prior_file = file(os.path.join(wd, 'priors.iter.%02d.txt'%current_iter))

    name2seq_i = {}
    name2prior = {}
    name2normed_prior = {}  # normed by length of sequences
    for line in prior_file:
        atoms = line.split()
        name2seq_i[atoms[1]] = int(atoms[0])
        name2prior[atoms[1]] = atoms[2]

    clustermark_pat = re.compile(r'(\d+\|.?\|)?(.*)')
    name2length = dict((clustermark_pat.search(record.description).groups()[1], len(record.seq)) for record in SeqIO.parse(file(os.path.join(wd, "iter.%02d.cons.fasta"%current_iter)), 'fasta'))


    name2normed_prior = dict((k, float(name2prior[k]) / length) for k, length in name2length.iteritems())
    # name2normed_prior = dict((k, e**(log(float(name2prior[k])) - log(float(length)))) for k, length in name2length.iteritems())
    sumproducts = float(sum((float(x) for x in name2normed_prior.values())))
    name2normed_prior = dict((k, v/sumproducts) for k,v in name2normed_prior.iteritems())


    sorted_records = []
    for record in SeqIO.parse(file(os.path.join(wd, "iter.%02d.cons.fasta"%current_iter)), "fasta"):
        name = clustermark_pat.search(record.description).groups()[1]  # strip off beginning cluster marks
        record.id = "%s%d|%s"%(record_prefix, name2seq_i[name], name)
        record.description = "Prior=%s Length=%d NormPrior=%06f"%(name2prior[name], len(record.seq), name2normed_prior[name])
        sorted_records.append((float(name2prior[name]), record)) 


    for prior, record in sorted(sorted_records, reverse=True):
        if prior < prob_min:
            break
        sys.stdout.write(record.format('fasta'))

def main(argv = None):
    if argv is None:
        argv = sys.argv[1:]

    parser = OptionParser("usage: %prog [options] <iter.DIR>  > renamed.fasta\n\nRewrites emirge fasta file to include proper sequence names and prior probabilities (abundance estimates) in record headers.\nIf no iter.DIR is given, assumes cwd.")
    parser.add_option('-p', '--prob_min',
                      type='float',
                      help='Only include sequences in output with prior probability above PROB_MIN (Default: include all sequences)')
    parser.add_option('-r', '--record_prefix',
                      type='string', default = '',
                      help='Add the specified prefix to each fasta record title')
    
    (options, args) = parser.parse_args(argv)
    if len(args) == 0:
        wd = os.getcwd()
    elif len(args) == 1:
        wd = os.path.abspath(args[0])
    else:
        parser.error("Found more than one argument on command line (expects only iter.DIR): %s"%(args))
    assert os.path.exists(wd), "Directory not found: %s"%(wd)
    rename(wd, options.prob_min, options.record_prefix)
    


if __name__ == '__main__':
    main()
