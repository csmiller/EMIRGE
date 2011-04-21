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
python emirge.py --help
"""

import sys
import os
import re
import csv
from optparse import OptionParser, OptionGroup
# from MyFasta import FastIterator, Record  # moved this code into this file.
import pysam
import numpy
from scipy import sparse
from subprocess import Popen, PIPE, check_call
from time import ctime, time
from datetime import timedelta
import gzip
import cPickle

BOWTIE_l = 20
BOWTIE_e  = 300

BOWTIE_ASCII_OFFSET = 33   # currently, bowtie writes quals with an ascii offset of 33 
# ILLUMINA_ASCII_OFFSET = 64  # not used... right now emirge.py assumes Illumina1.3+ 64-offset and sends appropriate command to bowtie.

class Record:
    """
    stripped down FASTA record class with same members as
    Biopython FASTA Record Class
        title    -- with > removed
        sequence -- as string
    """
    _colwidth = 60 # colwidth specifies the number of residues to put on each line when generating FASTA format.
    def __init__(self, title = "", sequence = ""):
        """
        Create a new Record.  
        
        """
        self.title = title
        self.sequence = sequence
    def __str__(self):
        return ">%s\n%s\n"%(self.title, "\n".join(re.findall(".{1,%s}"%self._colwidth, self.sequence)))

def FastIterator(filehandle, dummyParser = None, record = None):
    """
    a generator to return records one at a time.  Maybe 160% faster on test case.
    MAY RAISE MemoryError ON LARGE FASTA FILES
    IN:  file object
         dummyParser is a placeholder for RecordParser from Biopython.  Not used.
         a record to use for yielding.  Otherwise create an empty one with standard init

    NOTE: this fasta iterator is fast, but it breaks if there are ">" characters in the title.
    """
    if record is None:
        record = Record()
        
    for recordstring in re.split('\n>', filehandle.read()[1:]):
        record.title, record.sequence= recordstring.split('\n',1)
        record.sequence = record.sequence.replace('\n','').replace(' ','')
        yield record

class EM(object):
    """
    driver class for EM algorithm
    """
    _VERBOSE = True
    base2i = {"A":0,"T":1,"C":2,"G":3}
    i2base = dict([(v,k) for k,v in base2i.iteritems()])
    clustermark_pat = re.compile(r'(\d+\|.?\|)?(.*)')  # cludgey code associated with this should go into a method: get_seq_i()
    DEFAULT_ERROR = 0.05
    def __init__(self, reads1_filepath, reads2_filepath,
                 insert_mean,
                 insert_sd,
                 n_cpus = 1,
                 cwd = os.getcwd(), max_read_length = 76,
                 iterdir_prefix = "iter.", cluster_thresh = 0.97,
                 mapping_nice = None):
        """

        n_cpus is how many processors to use for multithreaded steps (currently only the bowtie mapping)
        mapping_nice is nice value to add to mapping program 
        """
        assert not reads1_filepath.endswith('.gz')
        # assert reads2_filepath.endswith('.gz')
        self.reads1_filepath = reads1_filepath
        self.reads2_filepath = reads2_filepath
        self.insert_mean = insert_mean
        self.insert_sd = insert_sd
        self.n_cpus = n_cpus
        self.mapping_nice   = mapping_nice

        self.iteration_i = None    # keeps track of which iteration we are on.
        self.cwd = cwd
        self.max_read_length = max_read_length
        self.iterdir_prefix = iterdir_prefix
        self.cluster_thresh = cluster_thresh   # if two sequences evolve to be >= cluster_thresh identical (via usearch), then merge them. [0, 1.0]
        assert self.cluster_thresh >= 0 and self.cluster_thresh <= 1.0
        # Single numpy arrays.  Has the shape: (numsequences x numreads)
        self.likelihoods = None      # = Pr(R_i|S_i), the likelihood of generating read R_i given sequence S_i
        # list of numpy arrays.  list index is iteration number.  Each numpy array has the shape: (numsequences,)
        self.priors      = []      # = Pr(S_i), the prior probability that sequence S generated any read
        # list of numpy arrays.  list index is iteration number.  Each numpy array has the shape: (numsequences x numreads)
        self.posteriors      = []  # = Pr(S_i|R_i), the posterior probability that sequence S_i generated read R_i
        # list of dictionaries.  list index is iteration number.  each dict is a mapping from sequence header name and internal id, or vice-versa
        # index is stable between iterations.  If sequence_i2sequence value is None, means this sequence was abandoned in previous round
        self.sequence_name2sequence_i = []
        self.sequence_i2sequence_name = []
        # list of strings.  list index is iteration number.
        
        # list of dictionaries.  list index is iteration number.  similar to above except for reads
        self.read_name2read_i = []
        self.read_i2read_name = []

        # this one for mapping between names with and without cluster numbers.
        self.sequence_name2fasta_name = {}

        # lists of numpy arrays of shape (max_read_length,)  Paired reads, so 0 = read 1, 1 = read 2
        # self.read_lengths = [[], []]
        # self.reads = [[], []]  # not used in current implementation
        self.quals = [[], []]

        # other constants, potentially changeable, tunable later, or could incorporate into probabilistic model.
        self.min_depth = 5.0    # minimum depth to keep sequence around for next round
        self.min_prior = None   # minimum prior probability for a sequence to keep it around for next round (alternative to depth, which is
                                # a little weird when you allow mappings to more than one place.
        self.min_length_coverage = None
        self.snp_minor_prob_thresh = 0.10      # if prob(N) for minor allele base N is >= this threshold, call site a minor allele
        self.snp_percentage_thresh = 0.10      # if >= this percentage of bases are minor alleles (according to self.snp_minor_prob_thresh),
                                               # then split this sequence into two sequences.
        
        return

    def read_bam(self, bam_filename, reference_fasta_filename):
        """
        reads a bam file and...
        populates a new entry for (appends to list, removes t-2)
                self.sequence_i2sequence_name
                self.sequence_name2sequence_i
                self.read_i2read_name           
                self.read_name2read_i
                
        creates a new EMPTY entry for (appends to list, removes t-2
                self.priors
                self.posteriors

        creates new each iteration (overwrites):
                # self.read_lengths
                self.likelihoods
                # self.reads
                self.quals
                self.probN
                self.unmapped_bases
                self.coverage
                self.sequence_name2fasta_name
                self.bamfile_data
                
        paired reads are treated as separate, individual reads.

        This MUST maintain seq_i to name and read_i to name mappings between iterations, so that a single
        name always maintains the same index from one iteration to the next.  One result of this requirement
        is that the various matrices can always get larger in a later t, but never smaller (as reads or seqs are added)
        """
        if self._VERBOSE:
            sys.stderr.write("Reading bam file %s at %s...\n"%(bam_filename, ctime()))

        initial_iteration = self.iteration_i < 0  # this is initial iteration            
        self.current_bam_filename = bam_filename
        self.current_reference_fasta_filename = reference_fasta_filename
        self.fastafile = pysam.Fastafile(self.current_reference_fasta_filename)
        

        for d in [self.sequence_name2sequence_i, self.sequence_i2sequence_name,
                  self.read_name2read_i, self.read_i2read_name]:
                  
            if initial_iteration:
                d.append({})
            else:
                d.append(d[-1].copy())  # get same name mappings from previous round, add to them if new reads or seqs
                if len(d) > 2:
                    trash = d.pop(0)  # no longer care about t-2
                    del trash
        # start new seq_i's or read_i's at next integer if there is a dictionary from the previous iteration.
        if initial_iteration:
            seq_i = 0
            read_i = 0
        else:
            seq_i = max(self.sequence_i2sequence_name[-2].keys()) + 1
            read_i = max(self.read_i2read_name[-2].keys()) + 1

        # FIRST PASS through file just to get number of reads, id mappings between names and integer ids
        # open bam file to get all sequences that have >= 1 mapping, and all reads
        bamfile = pysam.Samfile(bam_filename, "rb")
        getrname = bamfile.getrname  # speed hack
        # clustermark_pat_search = self.clustermark_pat.search

        # this is data struct to avoid lookups, since I read through bam file several times in an iteration.
        self.bamfile_data = [] # (refname, readname, seq_i, read_i, pair_i)
        bamfile_data = self.bamfile_data  # speed hack
        # this loop is pretty fast still (6 seconds on slow iterations)
        for alignedread_i, (alignedread) in enumerate(bamfile):
            refname = getrname(alignedread.rname)
            seq_i_to_cache = self.sequence_name2sequence_i[-1].get(refname, seq_i)
            if refname not in self.sequence_name2sequence_i[-1]:
                self.sequence_name2sequence_i[-1][refname] = seq_i
                self.sequence_i2sequence_name[-1][seq_i] = refname
                seq_i += 1
            pair_i = int(alignedread.is_read2)
            readname = "%s/%d"%(alignedread.qname, pair_i+1)
            read_i_to_cache = self.read_name2read_i[-1].get(readname, read_i)
            if readname not in self.read_name2read_i[-1]:
                self.read_name2read_i[-1][readname] = read_i
                self.read_i2read_name[-1][read_i] = readname
                read_i += 1
            bamfile_data.append((refname, readname, seq_i_to_cache, read_i_to_cache, pair_i))
        # self.read_lengths = numpy.zeros((2, read_i), dtype=numpy.uint8)
        # self.reads = numpy.zeros((2, read_i, self.max_read_length), dtype='c')
        self.quals = numpy.zeros((2, read_i, self.max_read_length), dtype=numpy.uint8)

        # these are pretty big in memory.  1 GB for 42272 reads and 1584 ref. sequences .
        # whole script uses about 3.5 GB in practice.
        # Most reads will have only one or a few sequences, so this is a very
        # sparse table.
        # one alternative is to use sparse matrices for the 2D arrays.
        # These act mostly like numpy arrays but if you need to convert, then numpy.array(a.todense())
        # self.priors.append(numpy.zeros(seq_i+1, dtype = numpy.float))
        # self.likelihoods.append(sparse.lil_matrix((seq_i+1, read_i+1)))
        # self.posteriors.append(sparse.lil_matrix((seq_i+1, read_i+1)))
        # really only need to keep t-1 and t, so maybe not so bad
        self.priors.append(numpy.zeros(seq_i, dtype = numpy.float))
        # self.likelihoods = numpy.zeros((seq_i, read_i), dtype = numpy.float)
        self.likelihoods = sparse.coo_matrix((seq_i, read_i), dtype = numpy.float)  # init all to zero.
        # self.posteriors.append(numpy.zeros_like(self.likelihoods))
        self.posteriors.append(sparse.lil_matrix((seq_i+1, read_i+1), dtype=numpy.float))

        self.coverage = numpy.zeros_like(self.priors[-1])

        # now second pass to populate (self.read_lengths,) XXself.readsXX, self.quals, self.coverage
        bamfile.close()
        bamfile = pysam.Samfile(bam_filename, "rb")

        # speed hacks:
        getrname = bamfile.getrname  
        ord_local = ord
        # reads = self.reads
        quals = self.quals
        coverage = self.coverage
        sequence_name2sequence_i = self.sequence_name2sequence_i[-1]
        read_name2read_i = self.read_name2read_i[-1]
        bamfile_data = self.bamfile_data
        ascii_offset = BOWTIE_ASCII_OFFSET
        # clustermark_pat_search = self.clustermark_pat.search

        # slowish.  1m8s... now 48s... now 42s on slow iteration
        # print >> sys.stderr, "***", ctime()
        
        for alignedread_i, alignedread in enumerate(bamfile):
            refname, readname, seq_i, read_i, pair_i = bamfile_data[alignedread_i]

            # refname = getrname(alignedread.rname)
            # seq_i = sequence_name2sequence_i[refname]
            # pair_i = int(alignedread.is_read2)
            # readname = "%s/%d"%(alignedread.qname, int(alignedread.is_read2)+1)
            # read_i = read_name2read_i[readname]

            # self.read_lengths[pair_i, read_i] = alignedread.rlen
            # reads[pair_i, read_i, :alignedread.rlen] = alignedread.seq
            quals[pair_i, read_i, :alignedread.rlen] = [ord_local(l)-ascii_offset for l in  alignedread.qual]  # could speed this up with ctypes
            coverage[seq_i] += alignedread.rlen
        bamfile.close()

        self.probN = [None for x in range(max(self.sequence_name2sequence_i[-1].values())+1)]
        self.unmapped_bases = [None for x in self.probN]
        
        # self.sequence_lengths = numpy.zeros(len(self.sequence_name2sequence_i[-1]), dtype=numpy.uint32)
        # print >> sys.stderr, "***", ctime()

        self.sequence_name2fasta_name = {}
        for record in FastIterator(file(self.current_reference_fasta_filename)):
            # refname = clustermark_pat_search(record.title.split()[0]).groups()[1]  # strip off beginning cluster marks
            refname = record.title.split()[0]
            self.sequence_name2fasta_name[refname] = record.title.split()[0]
            seq_i = self.sequence_name2sequence_i[-1].get(refname)
            if seq_i is not None:
                self.probN[seq_i] = numpy.zeros((len(record.sequence), 5), dtype=numpy.float)   #ATCG[other] --> 01234
                self.coverage[seq_i] = self.coverage[seq_i] / float(len(record.sequence))
                # self.sequence_lengths[seq_i] = len(record.sequence)

        for d in [self.priors, self.posteriors,
                  self.sequence_name2sequence_i, self.sequence_i2sequence_name,
                  self.read_name2read_i, self.read_i2read_name]:
            if len(d) > 2:
                trash = d.pop(0)  # no longer care about t-2
                del trash

        if self._VERBOSE:
            sys.stderr.write("DONE Reading bam file %s at %s...\n"%(bam_filename, ctime()))
        return

    def initialize_EM(self, bam_filename, reference_fasta_filename):
        """
        Set up EM with two things so that first iteration can proceed:
           - Initial guesses of Pr(S) are made purely based on read counts, where each read is only allowed to
             map only once to a single best reference  (**if more than one alignment reported per read, raise exception!**).
           - Initial guess of Pr(N=n) (necessary for likelihood in Pr(S|R) is also calculated simply, with the assumption
             of 1 read (the best again) mapped to exactly 1 sequence.  Thus Pr(N=n) only takes the base call errors
             into account.  This is actually not done here, but rather the first time self.calc_probN is called.
        
           - bamfile for iteration 0 is assumed to have just one ("best") mapping per read.
           - there is no t-1 for t = 0, hence the need to set up Pr(S)
        """
        if self._VERBOSE:
            sys.stderr.write("Beginning initialization at %s...\n"%(ctime()))

        self.iteration_i = -1
        self.read_bam(bam_filename, reference_fasta_filename)
        # initialize priors.  Here just adding a count for each read mapped to each reference sequence
        # also initialize Pr(N=n), which is the only other thing besides Pr(S) in the M step that depends on iteration t-1.
        bamfile = pysam.Samfile(self.current_bam_filename, "rb")
        getrname = bamfile.getrname
        sequence_name2sequence_i = self.sequence_name2sequence_i[-1]
        bamfile_data = self.bamfile_data
        # clustermark_pat_search = self.clustermark_pat.search

        # PRIOR
        for alignedread_i, alignedread in enumerate(bamfile):
            refname, readname, seq_i, read_i, pair_i = bamfile_data[alignedread_i]

            # refname = getrname(alignedread.rname)
            # seq_i = sequence_name2sequence_i[refname]
            self.priors[-1][seq_i] += 1
        bamfile.close()
        nonzero_indices = numpy.nonzero(self.priors[-1])  # only divide cells with at least one count.  Set all others to Pr(S) = 0
        self.priors[-1] = self.priors[-1][nonzero_indices] / self.priors[-1][nonzero_indices].sum()  # turn these into probabilities
        self.priors.append(self.priors[-1].copy())  # push this back to t-1 (index == -2)
        if self._VERBOSE:
            sys.stderr.write("DONE with initialization at %s...\n"%(ctime()))
        return

    def save_state(self, filename = None):
        """
        save state
        """
        tup_to_save = (self.bamfile_data, self.base2i, self.cluster_thresh, self.coverage, self.current_bam_filename, self.current_reference_fasta_filename, self.cwd, self.i2base, self.insert_mean, self.insert_sd, self.iteration_i, self.iterdir, self.iterdir_prefix, self.k, self.likelihoods, self.mapping_nice, self.max_read_length, self.min_depth, self.min_length_coverage, self.min_prior, self.n_cpus, self.posteriors, self.priors, self.probN, self.quals, self.read_i2read_name, self.read_name2read_i, self.reads1_filepath, self.reads2_filepath, self.sequence_i2sequence_name, self.sequence_name2fasta_name, self.sequence_name2sequence_i, self.snp_minor_prob_thresh, self.snp_percentage_thresh, self.unmapped_bases, self.v)
        if filename is None:
            filename = os.path.join(self.iterdir, 'em.%02i.data.pkl'%self.iteration_i)
        try:
            cPickle.dump(tup_to_save, file(filename, 'w'), cPickle.HIGHEST_PROTOCOL)
        except SystemError:  # cPickle problem with numpy arrays in latest emacs???
            sys.stderr.write("oops!  cPickle error!  Falling back to pickle.\n")
            import pickle
            pickle.dump(tup_to_save, file(filename, 'w'), pickle.HIGHEST_PROTOCOL)
        return filename
    
    def load_state(self, filename = None):
        """
        load pickled data structures stored with self.save_state
        """
        if filename is None:
            filename = os.path.join(self.iterdir, 'em.%02i.data.pkl'%self.iteration_i)

        if filename.endswith('bz2'):
            infile = Popen("bzcat %s"%filename, shell=True, stdout=PIPE).stdout
        else:
            infile = file(filename)

        for name in ("bamfile_data", "base2i", "cluster_thresh", "coverage", "current_bam_filename", "current_reference_fasta_filename", "cwd", "i2base", "insert_mean", "insert_sd", "iteration_i", "iterdir", "iterdir_prefix", "k", "likelihoods", "mapping_nice", "max_read_length", "min_depth", "min_length_coverage", "min_prior", "n_cpus", "posteriors", "priors", "probN", "quals", "read_i2read_name", "read_name2read_i", "reads1_filepath", "reads2_filepath", "sequence_i2sequence_name", "sequence_name2fasta_name", "sequence_name2sequence_i", "snp_minor_prob_thresh", "snp_percentage_thresh", "unmapped_bases", "v"):
            if not hasattr(self, name):
                setattr(self, name, None)
        try:
            (self.bamfile_data, self.base2i, self.cluster_thresh, self.coverage, self.current_bam_filename, self.current_reference_fasta_filename, self.cwd, self.i2base, self.insert_mean, self.insert_sd, self.iteration_i, self.iterdir, self.iterdir_prefix, self.k, self.likelihoods, self.mapping_nice, self.max_read_length, self.min_depth, self.min_length_coverage, self.min_prior, self.n_cpus, self.posteriors, self.priors, self.probN, self.quals, self.read_i2read_name, self.read_name2read_i, self.reads1_filepath, self.reads2_filepath, self.sequence_i2sequence_name, self.sequence_name2fasta_name, self.sequence_name2sequence_i, self.snp_minor_prob_thresh, self.snp_percentage_thresh, self.unmapped_bases, self.v) = \
                                cPickle.load(infile)
        except ValueError:  # old version didn't have bamfile_data
            (self.base2i, self.cluster_thresh, self.coverage, self.current_bam_filename, self.current_reference_fasta_filename, self.cwd, self.i2base, self.insert_mean, self.insert_sd, self.iteration_i, self.iterdir_prefix, self.k, self.likelihoods, self.mapping_nice, self.max_read_length, self.min_depth, self.min_length_coverage, self.min_prior, self.n_cpus, self.posteriors, self.priors, self.probN, self.quals, self.read_i2read_name, self.read_name2read_i, self.reads1_filepath, self.reads2_filepath, self.sequence_i2sequence_name, self.sequence_name2fasta_name, self.sequence_name2sequence_i, self.snp_minor_prob_thresh, self.snp_percentage_thresh, self.unmapped_bases, self.v) = \
                                cPickle.load(infile)                                
        self.fastafile = pysam.Fastafile(self.current_reference_fasta_filename)

        # change self.posteriors to sparse matrix if we are loading old data type
        if type(self.posteriors[-1]) != type(sparse.lil_matrix([1])):
            if self._VERBOSE:
                sys.stderr.write("\tConverting old posteriors to sparse matrix format [%d items to ... "%(self.posteriors[-1].shape[0] * self.posteriors[-1].shape[1]))
            for i, p in enumerate(self.posteriors):
                self.posteriors[i] = sparse.lil_matrix(p)
            if self._VERBOSE:
                sys.stderr.write("%d items] Done.\n"%(self.posteriors[-1].nnz))
                
        return
    
    def do_iteration(self, bam_filename, reference_fasta_filename):
        """
        This starts with the M-step, so it requires that Pr(S) and Pr(N=n) from previous round are set.
        Pr(S) is used from the previous round's E-step.
        Pr(N=n) partially depends on the previous round's M-step.  Is this okay?
        Once M-step is done, then E-step calculates Pr(S) based upon the just-calculated M-step.
        
        
        """
        self.iteration_i += 1
        if self._VERBOSE:
            sys.stderr.write("Starting iteration %d at %s...\n"%(self.iteration_i, ctime()))
            start_time = time()

        self.iterdir = os.path.join(self.cwd, "%s%02d"%(self.iterdir_prefix, self.iteration_i))
        check_call("mkdir -p %s"%(self.iterdir), shell=True)
        self.read_bam(bam_filename, reference_fasta_filename)  # initializes all data structures.

        # m-step
        self.calc_likelihoods()  
        self.calc_posteriors()

        # now e-step
        self.calc_priors()

        # now write a new fasta file.  Cull sequences below self.min_depth
        consensus_filename = os.path.join(self.iterdir, "iter.%02d.cons.fasta"%(self.iteration_i))
        self.write_consensus(consensus_filename)    # culls and splits
        self.cluster_sequences(consensus_filename)  # merges sequences that have evolved to be the same (USEARCH)

        # leave a few things around for later.  Note that print_priors also leaves sequence_name2sequence_i mapping, basically.
        self.print_priors()
        cPickle.dump(self.probN, gzip.GzipFile(os.path.join(self.iterdir, 'probN.pkl.gz'), 'w'), cPickle.HIGHEST_PROTOCOL)

        # now do a new mapping run for next iteration
        self.do_mapping(consensus_filename, nice = self.mapping_nice)
        
        if self._VERBOSE:
            sys.stderr.write("Finished iteration %d at %s...\n"%(self.iteration_i, ctime()))
            sys.stderr.write("Total time for iteration %d: %s\n"%(self.iteration_i, timedelta(seconds = time()-start_time)))
        
        return
    def print_priors(self, ofname = None):
        """
        leave a file in directory with nonzero priors printed out.
        """
        if ofname is not None:
            of = file(ofname, 'w')
        of = file(os.path.join(self.iterdir, "priors.iter.%02d.txt"%(self.iteration_i)), 'w')
        for seq_i, prior in enumerate(self.priors[-1]):
            seqname = self.sequence_i2sequence_name[-1][seq_i]
            of.write("%d\t%s\t%f\n"%(seq_i, seqname, prior))

        of.close()
             
            
    def calc_priors(self):
        """
        calculates priors [ Pr(S) ] based on
            Pr(S|R) (current posteriors from previous M step, this iteration)
        """
        # here we do have column summing with the posteriors; will have to sort this out with sparse.
        # should be csc
        self.posteriors[-1] = self.posteriors[-1].tocsc()
        self.priors[-1] = numpy.asarray(self.posteriors[-1].sum(axis = 1)).flatten() / self.posteriors[-1].sum()
        
        return

    def write_consensus(self, outputfilename):
        """
        writes a consensus, taking the most probable base at each position, according to current
        values in Pr(N=n) (self.probN)

        only write sequences above with coverage above self.min_depth (culling)
        split sequences with many minor alleles:
             self.snp_minor_prob_thresh     # if prob(N) for minor allele base N is >= this threshold, call site a minor allele
             self.snp_percentage_thresh     # if >= this percentage of bases are minor alleles (according to self.snp_minor_prob_thresh),
                                            # then split this sequence into two sequences.
        """
        if self._VERBOSE:
            sys.stderr.write("Writing consensus for iteration %d at %s...\n"%(self.iteration_i, ctime()))
            sys.stderr.write("\tsnp_minor_prob_thresh = %.3f\n"%(self.snp_minor_prob_thresh))
            sys.stderr.write("\tsnp_percentage_thresh = %.3f\n"%(self.snp_percentage_thresh))

        splitcount = 0
        cullcount  = 0
        of = file(outputfilename, 'w')
        
        i2base = self.i2base
        for seq_i, probNarray in enumerate(self.probN):
            if probNarray is None: # means this sequence is no longer present in this iteration
                continue          
            # check if coverage passes self.min_depth, if not don't write it (culling happens here)
            if self.min_depth is not None and self.coverage[seq_i] < self.min_depth:  # could also do this only after self.iteration_i > 5 or something
                # could adjust priors and posteriors here, but because prior will already be low (b/c of low coverage)
                # and because next round will have 0 mappings (no sequence in reference file to map to), this seems
                # unneccesary.
                cullcount += 1
                continue # continue == don't write it to consensus.
            # check if prior is above threshold... otherwise cull it:
            if self.min_prior is not None and self.priors[-1][seq_i] < self.min_prior:
                cullcount += 1
                continue
            if self.min_length_coverage is not None:
                num_mapped_indices = numpy.argwhere(probNarray > 1-default_error).shape[0]
                if float(num_mapped_indices) / float(probNarray.shape[0]) <= self.min_length_coverage:
                    cullcount += 1
                    continue
            # passes coverage threshold
            title = self.sequence_i2sequence_name[-1][seq_i]
            consensus = numpy.array([i2base.get(x, "N") for x in numpy.argsort(probNarray)[:,-1]])

            # check for minor allele consensus, SPLIT sequence into two candidate sequences if passes thresholds.
            minor_indices = numpy.argwhere((probNarray >= self.snp_minor_prob_thresh).sum(axis=1) >= 2)[:,0]
            if minor_indices.shape[0] / float(probNarray.shape[0]) >= self.snp_percentage_thresh:
                splitcount += 1
                if self._VERBOSE:
                    sys.stderr.write("splitting sequence %d at %s ..."%(seq_i, ctime()))
                minor_bases   = numpy.array([i2base.get(x, "N") for x in numpy.argsort(probNarray[minor_indices])[:,-2]]) # -2 gets second most probably base
                minor_consensus = consensus.copy()               # get a copy of the consensus
                minor_consensus[minor_indices] = minor_bases     # replace the bases that pass minor threshold
                # now deal with naming.  Just carry over short name in bam file, stripped of cluster info.
                title_root = re.search(r'(.+)(_m(\d+))$', title)
                if title_root is None: # no _m00 on this name 
                    title_root = title[:]
                else:
                    title_root = title_root.groups()[0]
                # now check for any known name with same root and a _m on it.
                previous_m_max = max([0] + [int(x) for x in re.findall(r'%s_m(\d+)'%title_root, " ".join(self.sequence_i2sequence_name[-1].values()))])
                m_title = "%s_m%02d"%(title_root, previous_m_max+1)

                # also split out Priors and Posteriors (which will be used in next round), split with average ratio of major to minor alleles.
                # updating priors first:
                minor_fraction_avg = numpy.average(probNarray[(minor_indices, numpy.argsort(probNarray[minor_indices])[:, -2])])
                major_fraction_avg = 1.-minor_fraction_avg # if there's >=3 alleles, major allele keeps prob of other minors)
                old_prior = self.priors[-1][seq_i]
                self.priors[-1][seq_i] = old_prior * major_fraction_avg
                seq_i_minor = max(self.sequence_i2sequence_name[-1].keys()) + 1
                self.sequence_i2sequence_name[-1][seq_i_minor] = m_title
                self.sequence_name2sequence_i[-1][m_title] = seq_i_minor
                new_priors = numpy.zeros(seq_i_minor+1, dtype=self.priors[-1].dtype)
                new_priors[:-1] = self.priors[-1].copy()
                new_priors[seq_i_minor] = old_prior * minor_fraction_avg
                trash = self.priors.pop()
                del trash
                self.priors.append(new_priors)

                # updating posteriors. for each seq-read pair with prob > 0, split prob out to major and minor seq.
                read_indices = numpy.nonzero(self.posteriors[-1][seq_i, :])[1]
                new_read_probs  = self.posteriors[-1][seq_i, read_indices].toarray().flatten() * minor_fraction_avg  # already returns copy
                # construct new_posteriors via coo, then convert to csr.
                new_posteriors = self.posteriors[-1].tocoo()  # first make a copy in coo format
                # then create new coo matrix with new shape, using same row, col, data as before.
                # leaves an empty row as last row, and the matrix is then converted to lil for slicing further down.
                new_posteriors = sparse.coo_matrix((new_posteriors.data, (new_posteriors.row, new_posteriors.col)),
                                                   shape=(seq_i_minor+1, self.posteriors[-1].shape[1]),
                                                   dtype=new_posteriors.dtype).tolil()
                # adjust old read probs
                new_posteriors[seq_i, read_indices] = new_posteriors[seq_i, read_indices] * major_fraction_avg
                # add new read probs
                new_posteriors[seq_i_minor, read_indices] = new_read_probs
                trash = self.posteriors.pop()
                del trash
                self.posteriors.append(new_posteriors.tocsr())
                # adjust self.unmapped_bases (used in clustering).  For now give same pattern as parent
                self.unmapped_bases.append(self.unmapped_bases[seq_i].copy())

                of.write(">%s\n"%(m_title))
                of.write("%s\n"%("".join(minor_consensus)))
                if self._VERBOSE:
                    sys.stderr.write("Done.\n")

            # now write major allele consensus, regardless of whether there was a minor allele consensus
            of.write(">%s\n"%(title))
            of.write("%s\n"%("".join(consensus)))

        if self._VERBOSE:
            sys.stderr.write("\tSplit out %d new minor allele sequences.\n"%(splitcount))
            sys.stderr.write("\tCulled %d sequences\n"%(cullcount))
            sys.stderr.write("DONE Writing consensus for iteration %d at %s...\n"%(self.iteration_i, ctime()))

        return

    def cluster_sequences(self, fastafilename):
        """
        uses Edgar's USEARCH to sort and then merge sequences above self.cluster_thresh %ID over the
        length of the shorter sequence

        "Search and clustering orders of magnitude faster than BLAST"
        Robert C. Edgar
        Bioinformatics 2010

        also adjusts Pr(S) [prior] and Pr(S_t-1) [posteriors] as needed after merging.  
        """
        return self.cluster_sequences2(fastafilename)

        return

    def cluster_sequences2(self, fastafilename):
        """
        uses USEARCH  to globally align sequences.  Merge two sequences if the
        *NON-GAPPED* positions have % identity >= self.cluster_thresh

        also adjusts Pr(S) [prior] and Pr(S_t-1) [posteriors] as needed after merging.
        """
        if self._VERBOSE:
            sys.stderr.write("Clustering sequences for iteration %d at %s...\n"%(self.iteration_i, ctime()))
            sys.stderr.write("\tcluster threshold = %.3f\n"%(self.cluster_thresh))

        # get posteriors ready for slicing:
        self.posteriors[-1] = self.posteriors[-1].tolil()
            
        # sort fasta sequences longest to shortest
        tmp_fastafilename = fastafilename + ".sorted.tmp.fasta"
        # check_call("usearch --sort %s --output %s"%(fastafilename, tmp_fastafilename), shell=True, stdout = sys.stdout, stderr = sys.stderr) # NEW version 4
        check_call("uclust --sort %s --output %s"%(fastafilename, tmp_fastafilename), shell=True, stdout = sys.stdout, stderr = sys.stderr) 
        tmp_fastafile = pysam.Fastafile(tmp_fastafilename)
        # do global alignments with USEARCH/UCLUST.
        # turn off banding in later rounds to increase alignment accuracy (more merging) at expense of speed
        band_string = ""
        # if em.iteration_i > 10:
        num_seqs = len([x for x in self.probN if x is not None])
        if num_seqs < 400:
            band_string = "--band 128"
        # if really few seqs, then no use not doing smith-waterman alignments
        if num_seqs < 50:
            band_string = "--nofastalign"
        cmd = "uclust %s --query %s --db %s --id 0 --iddef 0 --uc %s.uc --self --maxaccepts 0 --maxrejects 0 --allhits"%(band_string, tmp_fastafilename, tmp_fastafilename, tmp_fastafilename)
        # cmd = "usearch %s --cluster %s --db %s --id 0 --iddef 0 --uc %s.uc --maxaccepts 0 --maxrejects 0"%(band_string, tmp_fastafilename, tmp_fastafilename, tmp_fastafilename) # NEW version 4

        if self._VERBOSE:
            sys.stderr.write("usearch command was:\n%s\n"%(cmd))
        
        check_call(cmd, shell=True, stdout = sys.stdout, stderr = sys.stderr)
        # read clustering file to adjust Priors and Posteriors, summing merged reference sequences
        # Tab-separated fields:
        # 1=Type, 2=ClusterNr, 3=SeqLength or ClusterSize, 4=PctId, 5=Strand, 6=QueryStart, 7=SeedStart, 8=Alignment, 9=QueryLabel, 10=TargetLabel
        nummerged = 0
        alnstring_pat = re.compile(r'(\d*)([MDI])')
        already_removed = set()  # seq_ids
        # this is a bit slow and almost certainly could be sped up with algorithmic improvements.
        for row in csv.reader(file("%s.uc"%tmp_fastafilename), delimiter='\t'):
            if row[0] == "H":  # here's an alignment
                # member == query
                member_name = self.clustermark_pat.search(row[8]).groups()[1]  # strip off beginning cluster marks
                member_seq_id = self.sequence_name2sequence_i[-1][member_name]
                seed_name = self.clustermark_pat.search(row[9]).groups()[1]  # strip off beginning cluster marks
                if member_name == seed_name:
                    continue # new version of usearch doesn't officially support --self?
                seed_seq_id = self.sequence_name2sequence_i[-1][seed_name]
                if member_seq_id in already_removed or seed_seq_id in already_removed:
                    continue

                # decide if these pass the cluster_thresh *over non-gapped columns*
                member_i = 0
                seed_i   = 0
                matches = 0
                aln_columns    = 0
                member_fasta_seq = tmp_fastafile.fetch(member_name)
                seed_fasta_seq   = tmp_fastafile.fetch(seed_name)
                member_unmapped = self.unmapped_bases[member_seq_id]  # unmapped positions (default prob)
                seed_unmapped = self.unmapped_bases[seed_seq_id]
                
                for count, aln_code in alnstring_pat.findall(row[7]):
                    if not len(count):
                        count = 1  # "As a special case, if n=1 then n is omitted"
                    else:
                        count = int(count)
                    if aln_code == 'M':    # match (letter-letter column; could be mismatch)
                        s1 = member_fasta_seq[member_i:member_i+count]
                        s1_unmapped = member_unmapped[member_i:member_i+count]
                        s2 = seed_fasta_seq[seed_i:seed_i+count]
                        s2_unmapped = seed_unmapped[seed_i:seed_i+count]

                        both_columns_mapped = sum((s1_unmapped==False) &  (s2_unmapped==False))
                        aln_columns += both_columns_mapped  # count

                        # pure python (actually faster than numpy version below):
                        matches += sum(ch1 == ch2 for ch1, ch2, s1_u, s2_u in zip(s1, s2, s1_unmapped, s2_unmapped) if \
                                       not s1_u and not s2_u ) # only count match if both bases have read support
                        # numpy (a touch slower):
                        # matches += (numpy.array(s1, dtype='c') == numpy.array(s2, dtype='c')).sum()
                        member_i += count
                        seed_i += count
                    elif aln_code == 'D':  # seems like doc is backwards on this
                        member_i += count  
                    elif aln_code == 'I':
                        seed_i += count
                    else:
                        raise ValueError, "unknown alignment code: '%s'"%aln_code
                

                # DEBUG PRINT:
                if self._VERBOSE and num_seqs < 50:
                    print >> sys.stderr, "\t\t%s|%s vs %s|%s %.3f over %s aligned columns"%(member_seq_id, member_name, seed_seq_id, seed_name,
                                                                    float(matches) / aln_columns, aln_columns)
                # if alignment is less that 500 bases, or identity over those 500+ bases is not above thresh, then continue                
                if (aln_columns < 500) or ((float(matches) / aln_columns) < self.cluster_thresh):
                    continue

                # OKAY TO MERGE AT THIS POINT
                # if above thresh, then first decide which sequence to keep, (one with higher prior probability).
                if self.priors[-1][seed_seq_id] > self.priors[-1][member_seq_id]:
                    keep_seq_id = seed_seq_id
                    remove_seq_id = member_seq_id
                    keep_name = seed_name
                    remove_name = member_name
                else:
                    keep_seq_id = member_seq_id
                    remove_seq_id = seed_seq_id
                    keep_name   = member_name
                    remove_name = seed_name
                
                # merge priors (add remove_seq_id probs to keep_seq_id probs).
                self.priors[-1][keep_seq_id] += self.priors[-1][remove_seq_id]
                self.priors[-1][remove_seq_id] = 0.0
                # now merge posteriors (all remove probs go to keep).
                self.posteriors[-1][keep_seq_id, :] = self.posteriors[-1][keep_seq_id, :] + self.posteriors[-1][remove_seq_id, :]
                self.posteriors[-1][remove_seq_id, :] = self.posteriors[-1][remove_seq_id, :] * 0.0
                already_removed.add(remove_seq_id)
                nummerged += 1
                if self._VERBOSE:
                    sys.stderr.write("\t...merging %d|%s into %d|%s\n"%(remove_seq_id, remove_name,
                                                                        keep_seq_id,   keep_name))
            else:  # not an alignment line in input .uc file
                continue

        # write new fasta file with only new sequences
        if self._VERBOSE:
            sys.stderr.write("Writing new fasta file for iteration %d at %s...\n"%(self.iteration_i, ctime()))
        tmp_fastafile.close()
        outfile = file(fastafilename, 'w')
        for record in FastIterator(file(tmp_fastafilename)): # read through file again, overwriting orig file if we keep the seq
            seqname = self.clustermark_pat.search(record.title.split()[0]).groups()[1]  # strip off beginning cluster marks
            seq_id = self.sequence_name2sequence_i[-1][seqname]
            if seq_id not in already_removed:
                outfile.write(str(record))
        outfile.close()
            
        # clean up.
        check_call("sed -i 's/^>[0-9]\+|.\?|/>/g' %s"%(fastafilename), shell=True)  # remove cluster marks left by USEARCH (still needed?)
        os.remove("%s.uc"%(tmp_fastafilename))
        os.remove(tmp_fastafilename)
        os.remove("%s.fai"%(tmp_fastafilename))

        if self._VERBOSE:
            sys.stderr.write("\tremoved %d sequences after merging\n"%(nummerged))
            sys.stderr.write("DONE Clustering sequences for iteration %d at %s...\n"%(self.iteration_i, ctime()))

        return

    def do_mapping(self, full_fasta_path, nice = None):
        """
        IN:  path of fasta file to map reads to 
        run external mapping program to produce bam file
        right now this is bowtie
        """
        if self._VERBOSE:
            sys.stderr.write("Starting read mapping for iteration %d at %s...\n"%(self.iteration_i, ctime()))

        self.do_mapping_bowtie(full_fasta_path, nice = nice)

        if self._VERBOSE:
            sys.stderr.write("DONE with read mapping for iteration %d at %s...\n"%(self.iteration_i, ctime()))
        
        return
    def do_mapping_bowtie(self, full_fasta_path, nice = None):
        """
        run bowtie to produce bam file for next iteration
        
        """
        bowtie_logfile = os.path.join(self.iterdir, "bowtie.iter.%02d.log"%(self.iteration_i))
        bowtie_index   = os.path.join(self.iterdir, "bowtie.index.iter.%02d"%(self.iteration_i))
        # 1. build index
        cmd = "bowtie-build -o 3 %s %s >& %s"%(full_fasta_path , bowtie_index, bowtie_logfile) # -o 3 for speed? magnitude of speedup untested!
        if self._VERBOSE:
            sys.stderr.write("\tbowtie-build command:\n")
            sys.stderr.write("\t%s\n"%cmd)
        check_call(cmd, shell=True, stdout = sys.stdout, stderr = sys.stderr)
        
        # 2. run bowtie
        nice_command = ""
        if nice is not None:
            nice_command = "nice -n %d"%(nice)
        # these in theory could be used for single reads too.
        shared_bowtie_params = "--solexa1.3-quals -t -p %s  -n 3 -l %s -e %s  --best --strata --all --sam --chunkmbs 128"%(self.n_cpus, BOWTIE_l, BOWTIE_e)
        
        minins = max((self.insert_mean - 3*self.insert_sd), self.max_read_length)
        maxins = self.insert_mean + 3*self.insert_sd
        output_prefix = os.path.join(self.iterdir, "bowtie.iter.%02d"%(self.iteration_i))
        bowtie_command = "%s | %s bowtie %s -1 %s -2 - --minins %d --maxins %d %s  | samtools view -b -S -u -F 0x0004 - | samtools sort - %s.sort.PE >> %s 2>&1"%(\
            self.reads2_filepath,
            nice_command,
            shared_bowtie_params, self.reads1_filepath,
            minins, maxins,
            bowtie_index,
            output_prefix,
            bowtie_logfile)

        if self.reads2_filepath.endswith('.gz'):
            bowtie_command = "gzip -dc " + bowtie_command
        else:
            bowtie_command = "cat " + bowtie_command
            

        if self._VERBOSE:
            sys.stderr.write("\tbowtie command:\n")
            sys.stderr.write("\t%s\n"%bowtie_command)
        
        check_call(bowtie_command, shell=True, stdout = sys.stdout, stderr = sys.stderr)

        if self._VERBOSE:
            sys.stderr.write("\tFinished Bowtie for iteration %02d at %s:\n"%(self.iteration_i, ctime()))
            sys.stderr.write("\tRunning samtools to index bam file, and cleaning up\n")

        # 3. clean up (compress and create sam file for grepping, remove index files)
        check_call("samtools index %s.sort.PE.bam"%(output_prefix), shell=True, stdout = sys.stdout, stderr = sys.stderr)
        # commented out next because it is redundant
        # check_call("samtools view -h  %s.sort.PE.bam | gzip -c -f > %s.sort.PE.matches.sam.gz"%(output_prefix, output_prefix), shell=True)
        check_call("gzip -f %s"%(bowtie_logfile), shell=True)

        assert self.iterdir != '/'
        for filename in os.listdir(self.iterdir):
            assert(len(os.path.basename(bowtie_index)) >= 20)  # weak check that I'm not doing anything dumb.
            if os.path.basename(bowtie_index) in filename:
                os.remove(os.path.join(self.iterdir, filename))
        
        
        return
    def calc_likelihoods(self):
        """
        sets self.likelihoods  (seq_n x read_n) for this round
        """
        if self._VERBOSE:
            sys.stderr.write("Calculating likelihood %s for iteration %d at %s...\n"%(self.likelihoods.shape, self.iteration_i, ctime()))
        # first calculate self.probN from mapped reads, previous round's posteriors
        self.calc_probN()   # (handles initial iteration differently within this method)
        bamfile = pysam.Samfile(self.current_bam_filename, "rb")

        # for speed:
        probN = self.probN
        numpy_float = numpy.float
        quals       = self.quals
        sequence_name2sequence_i = self.sequence_name2sequence_i[-1]
        read_name2read_i         = self.read_name2read_i[-1]
        getrname = bamfile.getrname
        base2i_get = self.base2i.get
        arange = numpy.arange
        numpy_log = numpy.log
        e = numpy.e
        lik_row_seqi = []  # these for constructing coo_matrix for likelihood.
        lik_col_readi = []
        lik_data = []
        bamfile_data = self.bamfile_data
        zeros = numpy.zeros
        # clustermark_pat_search = self.clustermark_pat.search # no longer needed
        # could multi-thread this part here.  Only data structure I write to is self.likelihoods (here: likelihood), so I'd have to be careful about concurrent writes to that.
        for alignedread_i, alignedread in enumerate(bamfile):
            refname, readname, seq_i, read_i, pair_i = bamfile_data[alignedread_i]

            # refname = getrname(alignedread.rname)
            # seq_i = sequence_name2sequence_i[refname]
            # readname = "%s/%d"%(alignedread.qname, int(alignedread.is_read2)+1)
            # read_i = read_name2read_i[readname]
            # pair_i = int(alignedread.is_read2)
            prob_bases = zeros((alignedread.rlen, 5), dtype=numpy_float)
            # add P/3 to all bases.
            error_P = 10**(quals[pair_i, read_i, :alignedread.rlen] / -10.)
            prob_bases[:, :4] += (error_P / 3.).reshape((alignedread.rlen, 1))
            # subtract P/3 from called base
            numeric_bases = [base2i_get(base, 4) for base in alignedread.seq]
            read_pos_index = arange(alignedread.pos, alignedread.pos + alignedread.qlen)
            prob_bases[(arange(len(numeric_bases)), numeric_bases)] -= (error_P / 3.)
            # add (1-P) to called base
            prob_bases[(arange(len(numeric_bases)), numeric_bases)] += (1. - error_P)
            prob_b_i = (prob_bases[:, :4] * probN[seq_i][read_pos_index, :4]).sum(axis=1)  # this is Pr(b_i|S), where axis=0 is i
            # do at least this product in log space (0.94005726833168002 vs. 0.94005726833167991)
            # likelihood[seq_i, read_i] = e**(numpy_log(prob_b_i).sum())
            lik_row_seqi.append(seq_i)
            lik_col_readi.append(read_i)
            lik_data.append(e**(numpy_log(prob_b_i).sum()))
            
        bamfile.close()
        # now actually construct sparse matrix.
        self.likelihoods = sparse.coo_matrix((lik_data, (lik_row_seqi, lik_col_readi)), self.likelihoods.shape, dtype=self.likelihoods.dtype).tocsr()
        if self._VERBOSE:
            sys.stderr.write("DONE Calculating likelihood for iteration %d at %s...\n"%(self.iteration_i, ctime()))
        return

    def calc_posteriors(self):
        """
        Calculates Pr(S|R) for all sequence read pairs
        requires that the likelihood and priors are already calculated.
        """
        if self._VERBOSE:
            sys.stderr.write("Calculating posteriors for iteration %d at %s...\n"%(self.iteration_i, ctime()))
            t_start = time()

        # first populate with numerator of Bayes' theorum.  Use t-1 for priors, which corresponds to previous item in list/queue
        # do row by row now that I'm using sparse matrix
        lik_csr = self.likelihoods.tocsr()  # for efficient row slicing
        data = []  # these for constructing coo matrix
        ii   = []
        jj   = []
        for seq_i, this_prior in enumerate(self.priors[-2]):
            # self.posteriors[-1][seq_i, :] = l[seq_i, :] * this_prior
            this_row_coo = lik_csr[seq_i, :].tocoo()
            if this_row_coo.nnz == 0:
                continue
            data.append(this_row_coo.data * this_prior)
            jj.append(this_row_coo.col)
            ii.append(numpy.ones_like(jj[-1]) * seq_i)
        # build new coo_matrix for posteriors, convert to csc_matrix (for column slicing in calc_priors)
        data = numpy.concatenate(data)
        ii   = numpy.concatenate(ii)
        jj   = numpy.concatenate(jj)
        self.posteriors[-1] = sparse.coo_matrix((data, (ii, jj)), shape=self.likelihoods.shape, dtype=self.posteriors[-1].dtype)
        # now divide by sum of likelihoods*priors -- normalization factor (summed for each read over all possible sequences)
        denom = numpy.asarray(self.posteriors[-1].tocsc().sum(axis=0)).flatten()
        # only divide by nonzero denom -- this is taken care of by coo format!
        self.posteriors[-1].data = self.posteriors[-1].data / denom[(self.posteriors[-1].col,)]  # index out denom with column indices from coo format. 

        # convert to csc format for storage and use in self.calc_prior later.
        self.posteriors[-1] = self.posteriors[-1].tocsc()
        if self._VERBOSE:
            sys.stderr.write("DONE Calculating posteriors for iteration %d at %s [%.2f seconds]...\n"%(self.iteration_i, ctime(), time() - t_start))
        
    def calc_probN(self):
        """
        Pr(N=n)
        If read or sequence is new this round (not seen at t-1), then there is no Pr(S|R) from previous round, so we substitute Pr(S), the unbiased prior
        If initial iteration, all reads and seqs are new, so all calcs for Pr(N=n) use the prior as weighting factor instead of previous round's posterior.
        """
        if self._VERBOSE:
            sys.stderr.write("\tCalculating Pr(N=n) for iteration %d at %s...\n"%(self.iteration_i, ctime()))

        bamfile = pysam.Samfile(self.current_bam_filename, "rb")
        # basecounts = [seqprobNarray.astype(numpy.uint32) for seqprobNarray in self.probN]
        initial_iteration = self.iteration_i < 1
        
        # for speed:
        getrname = bamfile.getrname 
        arange = numpy.arange
        probN = self.probN
        if not initial_iteration:
            self.posteriors[-2] = self.posteriors[-2].tolil()
            posteriors = self.posteriors[-2]  # this depends on PREVIOUS iteration's posteriors (seq_n x read_n)
        priors     = self.priors[-2]          # and if no posteriors are present (or initial iter), falls back on priors from previous round
        base2i_get = self.base2i.get
        sequence_name2sequence_i = self.sequence_name2sequence_i[-1]
        read_name2read_i = self.read_name2read_i[-1]
        bamfile_data = self.bamfile_data
        quals = self.quals
        # clustermark_pat_search = self.clustermark_pat.search

        # could multithread this too.  self.probN is what is being written to.
        for alignedread_i, alignedread in enumerate(bamfile):
            refname, readname, seq_i, read_i, pair_i = bamfile_data[alignedread_i]

            # refname = getrname(alignedread.rname)
            # seq_i = sequence_name2sequence_i[refname]
            # readname = "%s/%d"%(alignedread.qname, int(alignedread.is_read2)+1)
            # read_i = read_name2read_i[readname]
            if initial_iteration:
                weight = priors[seq_i]
            else:
                try:
                    weight = posteriors[seq_i, read_i]
                except IndexError:
                    weight = priors[seq_i]
                
            numeric_bases = [base2i_get(base, 4) for base in alignedread.seq]
            read_pos_index = arange(alignedread.pos, alignedread.pos + alignedread.qlen)
            # basecounts[seq_i][(read_pos_index, numeric_bases)] += 1
            # pair_i = int(alignedread.is_read2)
            # add P/3 * Pr(S_t-1 | R) to all bases.
            error_P = 10**(quals[pair_i, read_i, :alignedread.rlen] / -10.)
            weighted_vals = (error_P / 3.) * weight
            # try:
            probN[seq_i][alignedread.pos: alignedread.pos+alignedread.qlen, 0:4] += weighted_vals.reshape((alignedread.rlen, 1)) # TODO: broadcasting correct?
            # except:
            #     print alignedread
            #     print seq_i
            #     print self.sequence_i2sequence_name[-1][seq_i]
            #     raise
            # subtract P/3 from called base
            probN[seq_i][(read_pos_index, numeric_bases)] -= weighted_vals
            # add (1-P) to called base
            probN[seq_i][(read_pos_index, numeric_bases)] += ((1. - error_P) * weight)
            # TODO: figure out if all this can be kept in log space... rounding errors might be important
        bamfile.close()

        numpy_where = numpy.where
        numpy_nonzero = numpy.nonzero
        # Here is Pr(N=n) = 0.95 (set default error_P to 0.05); kind of arbitrary.
        # NOW setting so all other bases < snp_minor_prob_thresh, but ref base not as high as before.
        # default_error = 1 - (3. * self.snp_minor_prob_thresh * 0.95) # make sure don't create snps by splitting
        default_error = self.DEFAULT_ERROR
        for seq_i, probNarray in enumerate(probN):
            if probNarray is None:  # means this sequence is no longer present in this iteration
                self.unmapped_bases[seq_i] = None
                continue
            # only divide cells with at least one base mapped.
            nonzero_indices = numpy_nonzero(probNarray.sum(axis=1))
            nonzero_probsums = probNarray.sum(axis=1)[nonzero_indices[0]]
            nonzero_probsums = nonzero_probsums.reshape((nonzero_probsums.shape[0], 1))
            probN[seq_i][nonzero_indices[0]] = probNarray[nonzero_indices[0]] / nonzero_probsums

            # bases with no mappings now -- (arbitrarily) set Pr(N=n) to 0.95 where n = reference base
            # 
            zero_indices    = numpy_where(probNarray.sum(axis=1) == 0)
            self.unmapped_bases[seq_i] = numpy.zeros(probNarray.shape[0], dtype=numpy.bool)
            self.unmapped_bases[seq_i][zero_indices[0]] = True
            if zero_indices[0].shape[0] > 0:  # there are bases without mappings.
                fastaname = self.sequence_name2fasta_name[self.sequence_i2sequence_name[-1][seq_i]]
                bases = numpy.array(self.fastafile.fetch(fastaname), dtype='c')[zero_indices[0]]
                numeric_bases = [base2i_get(base, 4) for base in bases]
                error_P = numpy.zeros(zero_indices[0].shape) + default_error
                # add P/3 to all bases without mappings.
                probNarray[zero_indices[0], :4] += (error_P / 3.).reshape((len(zero_indices[0]), 1)) # TODO: check this broadcasting...
                # subtract P/3 from previous reference base
                probNarray[(zero_indices[0], numeric_bases)] -= (error_P / 3.)
                # add (1-P) to previous reference base
                probNarray[(zero_indices[0], numeric_bases)] += (1. - error_P)
                # TODO: figure out if all this can be kept in log space... rounding errors might be +/- 3e-12
        if self._VERBOSE:
            sys.stderr.write("\tDONE calculating Pr(N=n) for iteration %d at %s...\n"%(self.iteration_i, ctime()))
        return
        

    def update_prior(self):
        """
        Updates Pr(S) based on the calculated Pr(S|R) from previous M step.
        """

        return
    
    def write_fastq_for_seq_i(self, seq_i, output_prefix = None):
        """
        for a specific seq_i, write the reads that most of their
        probability (>50%) assigned to this sequence into output_prefix.fastq
        as a fastq file.
        """
        if output_prefix is None:
            output_prefix = os.path.join(em.iterdir, "%s.reads"%seq_i)
        of_fastq = file('%s.fastq'%(output_prefix), 'w')
        # go through bam file instead of original sequencing reads, as aligned reads only a fraction
        # of original file.
        self.posteriors[-1] = self.posteriors[-1].tolil()  # seq_i x read_i
        posteriors = self.posteriors[-1]
        bamfile = pysam.Samfile(self.current_bam_filename, "rb")
        bamfile_data = self.bamfile_data
        reads = 0
        for alignedread_i, alignedread in enumerate(bamfile):
            refname, readname, this_seq_i, read_i, pair_i = bamfile_data[alignedread_i]
            if this_seq_i != seq_i or posteriors[seq_i, read_i] < 0.5:
                continue
            # means we have a read with more than 50% prob assigned to this sequence
            of_fastq.write("@%s\n%s\n+\n%s\n"%(readname, alignedread.seq, alignedread.qual))
            reads += 1
        of_fastq.close()
        if self._VERBOSE:
            sys.stderr.write("Wrote %d sequences to %s\n"%(reads, of_fastq.name))
        return
    def write_sam_for_seq_i(self, seq_i, output_prefix = None):
        """
        for a specific seq_i, write the readmappings for reads with most of their
        probability (>50%) assigned to this sequence
        also write a reference fasta file for Velvet's Columbus module.

        PREFIX.mappings.sam
        PREFIX.ref.fasta
        
        """
        if output_prefix is None:
            output_prefix = os.path.join(em.iterdir, "%s.reads"%seq_i)
        
        # go through bam file instead of original sequencing reads, as aligned reads only a fraction
        # of original file.
        self.posteriors[-1] = self.posteriors[-1].tolil()  # seq_i x read_i
        posteriors = self.posteriors[-1]  # seq_i x read_i
        bamfile = pysam.Samfile(self.current_bam_filename, "rb")
        of_sam_name = '%s.mappings.sam'%(output_prefix)
        of_sam = pysam.Samfile(of_sam_name, 'w', template = bamfile)
        bamfile_data = self.bamfile_data
        reads = 0
        for alignedread_i, alignedread in enumerate(bamfile):
            refname, readname, this_seq_i, read_i, pair_i = bamfile_data[alignedread_i]
            if this_seq_i != seq_i or posteriors[seq_i, read_i] < 0.5:
                continue
            # means we have a read with more than 50% prob assigned to this sequence
            # MARK!!!
            of_sam.write(alignedread)
            reads += 1
        of_sam.close()
        of_fasta = file('%s.ref.fasta'%(output_prefix), 'w')
        fasta_seq   = self.fastafile.fetch(self.sequence_i2sequence_name[-1][seq_i])
        of_fasta.write("%s"%(str(Record(self.sequence_i2sequence_name[-1][seq_i], fasta_seq))))
        of_fasta.close()
        if self._VERBOSE:
            sys.stderr.write("Wrote %d sequences to %s\n"%(reads, of_sam_name))
        return

    def write_phrap_for_seq_i(self, seq_i, output_prefix = None):
        """
        for a specific seq_i, write the reads and quality scores for reads with most of their
        probability (>50%) assigned to this sequence
        also write a reference fasta file for Velvet's Columbus module.

        PREFIX.reads.fasta
        PREFIX.quals.fasta
        
        """
        if output_prefix is None:
            output_prefix = os.path.join(em.iterdir, "%s"%seq_i)
        
        # go through bam file instead of original sequencing reads, as aligned reads only a fraction
        # of original file.
        self.posteriors[-1] = self.posteriors[-1].tolil()  # seq_i x read_i
        posteriors = self.posteriors[-1]  # seq_i x read_i
        bamfile = pysam.Samfile(self.current_bam_filename, "rb")
        of_fasta_name = '%s.reads.fasta'%(output_prefix)
        of_fasta = file(of_fasta_name, 'w')
        bamfile_data = self.bamfile_data
        reads = 0
        for alignedread_i, alignedread in enumerate(bamfile):
            refname, readname, this_seq_i, read_i, pair_i = bamfile_data[alignedread_i]
            if this_seq_i != seq_i or posteriors[seq_i, read_i] < 0.5:
                continue
            # means we have a read with more than 50% prob assigned to this sequence
            header = "%s DIRECTION: fwd CHEM: unknown TEMPLATE: %s"%(readname, readname)
            # header = "%s"%(reads)
            of_fasta.write(">%s\n%s\n"%(header, alignedread.seq))
            #, alignedread.qual))
            reads += 1
        of_fasta.close()
        if self._VERBOSE:
            sys.stderr.write("Wrote %d sequences to %s\n"%(reads, of_fasta_name))
        return

    def iterations_done(self):
        """
        check if we are done iterating, i.e. are the current reference sequences the same as that from the last round

        returns True or False
        """

        return False


def test_initialize(n_cpus = 10, cwd = os.getcwd(),
                    bam_ref   = '/work/csmiller/16S/emess/SSURef_102_tax_silva.sorted.fixed.97.fasta.sort.PE.bam',
                    fasta_ref = '/work/csmiller/16S/emess/SSURef_102_tax_silva.sorted.fixed.97.fasta',
                    mapping_nice = None):
    """
    This one used for 5WS data.

    
    """
    # 
    # bam_ref   = '/work/csmiller/16S/emess/SSURef_102_tax_silva.sorted.fixed.97.fasta.sort.PE.bam'
    # fasta_ref = '/work/csmiller/16S/emess/SSURef_102_tax_silva.sorted.fixed.97.fasta'

    # pe1_file = "/work/csmiller/min44/887.937.951.PE.1.fastq"
    # pe2_file = "/work/csmiller/min44/887.937.951.PE.2.fastq.gz"

    pe1_file = "/work/csmiller/Btrim60/951.Btrim60.PE.1.fastq"
    pe2_file = "/work/csmiller/Btrim60/951.Btrim60.PE.2.fastq.gz"

    em = EM(reads1_filepath = pe1_file,
            reads2_filepath = pe2_file,
            insert_mean = 194,
            insert_sd   = 25,
            n_cpus = n_cpus,
            cwd = cwd,
            mapping_nice = mapping_nice)
    em.snp_percentage_thresh = 0.04
    em.min_depth = 1
    em.initialize_EM(bam_ref, fasta_ref)
    return em

def test_generic(cwd, bam_ref, fasta_ref, 
                 pe1_file, pe2_file,
                 insert_mean, insert_sd,
                 max_read_length,
                 snp_percentage_thresh = 0.04,
                 n_cpus = 10
                 ):
    """
    INITIAL BOWTIE COMMAND (right now run separately.  This reports 1 alignment per read.):

    [csmiller@hydra 16S]> pwd
    /localdisk1/work/csmiller/16S
    [csmiller@hydra 16S]> bowtie-build SSURef_102_tax_silva.sorted.fixed.97.fasta ./bowtie_indices/SSURef_102_tax_silva.sorted.fixed.97

    gzip -dc /work/csmiller/Singer/GTXX.Btrim60.PE.2.fastq.gz | bowtie -t -p 12 -n 3 -l 20 -e 300 --best --sam --chunkmbs 128 -1 /work/csmiller/Singer/GTXX.Btrim60.PE.1.fastq -2 - --minins 123 --maxins 339 /work/csmiller/16S/bowtie_indices/SSURef_102_tax_silva.sorted.fixed.97 | samtools view -b -S -u -F 0x0004 - | samtools sort - /work/csmiller/Singer/emess/initial_mapping/GTXX.Btrim60.PE.sort >> & /work/csmiller/Singer/emess/initial_mapping/bowtie.log


    # for singer GTXX: 
    pe1_file = "/work/csmiller/Singer/GTXX.Btrim60.PE.1.fastq"
    pe2_file = "/work/csmiller/Singer/GTXX.Btrim60.PE.2.fastq.gz"
    insert_mean = 231
    insert_sd = 36
    """
    em = EM(reads1_filepath = pe1_file,
            reads2_filepath = pe2_file,
            insert_mean = insert_mean,
            insert_sd   = insert_sd,
            max_read_length = max_read_length,
            n_cpus = n_cpus,
            cwd = cwd)
    em.snp_percentage_thresh = snp_percentage_thresh
    em.min_depth = 1
    em.initialize_EM(bam_ref, fasta_ref)

    return em

def hack_it_together(max_iter):
    """
    until I get that UI going... I've got to see what it does!

    ASSUMES A GLOBAL VARIABLE CALLED em
    """
    assert 'em' in globals()
    global em
    os.chdir(em.cwd)

    # em = test_initialize(n_cpus = 12, cwd = workdir)
    # em.snp_percentage_thresh = 0.04
    if em.iteration_i < 0:  # first run
        em.do_iteration(em.current_bam_filename, em.current_reference_fasta_filename)

    while em.iteration_i < max_iter:
        subdir = os.path.join(em.cwd, "iter.%02d"%(em.iteration_i))
        em.do_iteration(os.path.join(subdir, "bowtie.iter.%02d.sort.PE.bam"%(em.iteration_i)),
                        os.path.join(subdir, "iter.%02d.cons.fasta"%(em.iteration_i)))

def hack_all():
    assert 'em' in globals()
    global em
    os.chdir(em.cwd)

    for newmax in range(10,130, 10):
	hack_it_together(newmax)
 	filename = em.save_state()
        os.system("bzip2 -f %s &"%(filename))

def countem(thresh):
    i = 0
    indices = []
    for seq_i, probNarray in enumerate(em.probN):
        if probNarray is None:
            continue
        try:
            num_mapped_indices = numpy.argwhere(probNarray > (1-0.05)).shape[0]
        except:
            print seq_i
            raise
 	if float(num_mapped_indices) / float(probNarray.shape[0]) <= thresh:
            indices.append(seq_i)
    return indices

def do_Rifle():
    global em
    wd = "/localdisk1/work/csmiller/16S/emess_Rifle"
    em = test_generic(wd, os.path.join(wd, "s_7_trim_B.min60.PE.sort.bam"), "/work/csmiller/16S/good_reference_database/SSURef_102_tax_silva.sorted.fixed.97.fasta", os.path.join(wd, "s_7_trim_B.min60.PE.1.fastq"), os.path.join(wd, "s_7_trim_B.min60.PE.2.fastq.gz"), 479, 195, max_read_length = 105, snp_percentage_thresh = 0.04, n_cpus = 20)
    em.min_depth = 3
    em.mapping_nice = 4
    hack_all()

def do_iterations(em, max_iter, save_every):
    """
    an EM object is passed in, so that one could in theory start from a saved state
    """
    os.chdir(em.cwd)

    if em.iteration_i < 0:  # first run
        em.do_iteration(em.current_bam_filename, em.current_reference_fasta_filename)

    while em.iteration_i < max_iter:
        subdir = os.path.join(em.cwd, "iter.%02d"%(em.iteration_i))
        em.do_iteration(os.path.join(subdir, "bowtie.iter.%02d.sort.PE.bam"%(em.iteration_i)),
                        os.path.join(subdir, "iter.%02d.cons.fasta"%(em.iteration_i)))
        if em.iteration_i > 0 and (em.iteration_i % save_every == 0):
            filename = em.save_state()
            os.system("bzip2 -f %s &"%(filename))
    return

def do_initial_mapping(working_dir, options):
    """
    IN:  takes the working directory and an OptionParser options object

    does the initial 1-reference-per-read bowtie mapping to initialize the algorithm
    OUT:  path to the bam file from this initial mapping
    """

    initial_mapping_dir = os.path.join(working_dir, "initial_mapping")
    if not os.path.exists(initial_mapping_dir):
        os.mkdir(initial_mapping_dir)
    # Template:
    # gzip -dc /work/csmiller/Singer/GSGA.Btrim60.PE.2.fastq.gz | bowtie --solexa1.3-quals -t -p 12 -n 3 -l 20 -e 300 --best --sam --chunkmbs 128 -1 /work/csmiller/Singer/GSGA.Btrim60.PE.1.fastq -2 - --minins 78 --maxins 336 ./bowtie_indices/SSURef_102_tax_silva.sorted.fixed.97 | samtools view -b -S -u -F 0x0004 - | samtools sort - GSGA.Btrim60.PE.sort >> & initial_mapping.log &

    minins = max((options.insert_mean - 3*options.insert_stddev), options.max_read_length)
    maxins = options.insert_mean + 3*options.insert_stddev
    bampath_prefix = os.path.join(initial_mapping_dir, "initial_bowtie_mapping.sorted")

    nicestring = ""
    if options.nice_mapping is not None:
        nicestring = "nice -n %d"%(options.nice_mapping)  # TODO: fix this so it isn't such a hack and will work in bash.  Need to rewrite all subprocess code, really (shell=False)
    option_strings = [options.fastq_reads_2, nicestring, options.processors, BOWTIE_l, BOWTIE_e, options.fastq_reads_1, minins, maxins, options.bowtie_db, bampath_prefix]
    if options.fastq_reads_2.endswith(".gz"):
        option_strings = ["gzip -dc "] + option_strings
    else:
        option_strings = ["cat "] + option_strings

    cmd = "%s %s | %s bowtie --solexa1.3-quals -t -p %s -n 3 -l %s -e %s --best --sam --chunkmbs 128 -1 %s -2 - --minins %s --maxins %s %s | samtools view -b -S -u -F 0x0004 - | samtools sort - %s "%tuple(option_strings)    
    print "Performing initial mapping with command:\n%s"%cmd
    check_call(cmd, shell=True, stdout = sys.stdout, stderr = sys.stderr)
    sys.stdout.flush()
    sys.stderr.flush()
    os.fsync(sys.stdout.fileno())  # make sure it gets synced?
    os.fsync(sys.stderr.fileno())  # make sure it gets synced?
    
    return bampath_prefix+".bam"

def resume(working_dir, options):
    """
    resume from a previous run.

    Takes the emirge working dir, and an OptionParser options object
    """
    raise NotImplementedError, "This option is currently broken, and will be fixed in a later version."
    em = EM("", "", 0, 0)  # reads1_filepath, reads2_filepath, insert_mean, insert_sd
    data_path = os.path.join(working_dir, "iter.%02d"%(options.resume_from), 'em.%02i.data.pkl.bz2'%(options.resume_from))
    sys.stdout.write("Loading saved state from %s...\n"%data_path)
    em.load_state(data_path)
    sys.stdout.write("Done.\n")

    # if the current working dir has been copied or moved, the old state will no longer be valid,
    # so need to set this again:
    em.cwd = working_dir

    # secretly (not advertised) options allowed to change in a resume
    if options.fastq_reads_1 is not None:
        em.reads1_filepath = os.path.abspath(options.fastq_reads_1)
    if options.fastq_reads_2 is not None:
        em.reads2_filepath = os.path.abspath(options.fastq_reads_2)

    # now process any *relevant* options:
    # this is broken right now because it just reverts all to defaults.
    # if options.processors is not None:
    #     em.n_cpus = options.processors
    # if options.snp_fraction_thresh is not None:
    #     em.snp_percentage_thresh = options.snp_fraction_thresh
    # if options.variant_fraction_thresh is not None:
    #     em.snp_minor_prob_thresh = options.variant_fraction_thresh
    # if options.join_threshold is not None:
    #     em.cluster_thresh = options.join_threshold
    # if options.min_depth is not None:
    #     em.min_depth = options.min_depth
    # if options.nice_mapping is not None:
    #     em.mapping_nice = options.nice_mapping

    do_iterations(em, max_iter = options.iterations, save_every = options.save_every)
    return

def main(argv = sys.argv[1:]):
    """
    command line interface to emirge
    
    """
    parser = OptionParser("usage: %prog DIR <required_parameters> [options]\n\nEMIRGE attempts to reconstruct rRNA SSU genes from Illumina metagenomic data.\n DIR is the working directory to process data in.\nuse --help to see a list of required and optional arguments")


    # REQUIRED
    group_reqd = OptionGroup(parser, "Required flags",
                             "These flags are all required to run EMIRGE (and may be supplied in any order)")

    group_reqd.add_option("-1", dest="fastq_reads_1", metavar="reads_1.fastq",
                      type="string",
                      help="path to fastq file with \\1 (forward) reads from paired-end run.  Must be unzipped for mapper.  EMIRGE expects ASCII-ofset of 64 for quality scores.  Required.")
    group_reqd.add_option("-2", dest="fastq_reads_2", metavar="reads_2.fastq[.gz]",
                      type="string",
                      help="path to fastq file with \\2 (reverse) reads from paired-end run.  May be gzipped for mapper.  EMIRGE expects ASCII-ofset of 64 for quality scores.  Required.")
    group_reqd.add_option("-f", "--fasta_db",
                      type="string",
                      help="path to fasta file of candidate SSU sequences")
    group_reqd.add_option("-b", "--bowtie_db",
                      type="string",
                      help="precomputed bowtie index of candidate SSU sequences (path to appropriate prefix; see --fasta_db)")
    group_reqd.add_option("-i", "--insert_mean",
                      type="int",
                      help="insert size distribution mean.  Required.")
    group_reqd.add_option("-s", "--insert_stddev",
                      type="int",
                      help="insert size distribution standard deviation.  Required.")
    group_reqd.add_option("-l", "--max_read_length",
                      type="int",
                      help="""length of longest read in input data.  Required.""")
    parser.add_option_group(group_reqd)

    # OPTIONAL
    group_opt = OptionGroup(parser, "Optional parameters",
                             "Defaults should normally be fine for these options in order to run EMIRGE")
    group_opt.add_option("-n", "--iterations",
                      type="int", default=40,
                      help="""Number of iterations to perform.  (default=%default)""")
    group_opt.add_option("-a", "--processors",
                      type="int", default=1,
                      help="""Number of processors to use in the mapping steps.  You probably want to raise this if you have the processors. (default: %default)""")
    group_opt.add_option("-m", "--mapping",
                         type="string",
                         help="path to precomputed initial mapping (bam file).  If not provided, and initial mapping will be run for you.")
    group_opt.add_option("-p", "--snp_fraction_thresh",
                      type="float", default="0.04",
                      help="If fraction of variants in a candidate sequence exceeds this threhold, then split the candidate into two sequences for next iteration.  See also --variant_fraction_thresh. (default: %default)")
    group_opt.add_option("-v", "--variant_fraction_thresh",
                      type="float", default="0.1",
                      help="minimum probability of second most probable base at a site required in order to call site a variant.  See also --snp_fraction_thresh.  (default: %default)")
    group_opt.add_option("-j", "--join_threshold",
                      type="float", default="0.97",
                      help="If two candidate sequences share >= this fractional identity over their bases with mapped reads, then merge the two sequences into one for the next iteration.  (default: %default) ")
    group_opt.add_option("-c", "--min_depth",
                      type="float",
                      default = 3,
                      help = "minimum average read depth below which a candidate sequence is discarded for next iteration(default: %default)")
    group_opt.add_option("--nice_mapping",
                      type="int",
                      help="""If set, during mapping phase, the mapper will be "niced" by the Linux kernel with this value (default: no nice)""")
    group_opt.add_option("-e", "--save_every",
                      type="int", default=10,
                      help="""every SAVE_EVERY iterations, save the programs state.  This allows you to run further iterations later starting from these save points.  The program will always save its state after the final iteration.  (default=%default)""")

    parser.add_option_group(group_opt)
    # # RESUME
    # group_resume = OptionGroup(parser, "Resuming iterations",
    #                          "These options allow you to resume iterations from a previous saved state.  Other options set on the command line are ignored")
    #                            # With --resume_from, all optional flags except --mapping are also allowed.  You'll likely want to set --iterations as well.  If you do set other options, they overwrite the defaults otherwise saved in the iteration being resumed from")
    # group_resume.add_option("-r", "--resume_from",
    #                         type="int",
    #                         help="Resume iterations from saved state in iteration specified.  Requires that a em data pickle was saved with the --save_every flag for that iteration (em.N.data.pkl.bz2)")
    # parser.add_option_group(group_resume)

    # ACTUALLY PARSE ARGS
    (options, args) = parser.parse_args(argv)

    if len(args) !=1:
        parser.error("DIR is required, and all options except DIR should have a flag associated with them (options without flags: %s)"%args)
    working_dir = os.path.abspath(args[0])

    sys.stdout.write("Command:\n")
    sys.stdout.write(' '.join([__file__]+argv))
    sys.stdout.write('\n\n')
    sys.stdout.flush()
    os.fsync(sys.stdout.fileno())  # make sure it gets synced?

    # first handle RESUME case
    # if options.resume_from is not None:
    #     resume(working_dir, options)
    #     return # ends the program, as we resumed from previous run

    # below here, means that we are handling the NEW case (as opposed to resume)
    if sum([int(x is None) for x in [options.fastq_reads_1, options.fastq_reads_2, options.insert_mean, options.insert_stddev, options.max_read_length]]):
        parser.error("Some required arguments are missing (try --help)")

    if not os.path.exists(working_dir):
        os.mkdir(working_dir)

    # clean up options to be absolute paths
    for o in ["fastq_reads_1", "fastq_reads_2", "fasta_db", "bowtie_db", "mapping"]:
        current_o_value = getattr(options, o)
        if current_o_value is not None:
            setattr(options, o, os.path.abspath(current_o_value))

    # DO INITIAL MAPPING if not provided with --mapping
    if options.mapping is None:
        options.mapping = do_initial_mapping(working_dir, options)


    # finally, CREATE EM OBJECT
    em = EM(reads1_filepath = options.fastq_reads_1,
            reads2_filepath = options.fastq_reads_2,
            insert_mean = options.insert_mean,
            insert_sd   = options.insert_stddev,
            max_read_length = options.max_read_length,
            cluster_thresh = options.join_threshold,
            n_cpus = options.processors,
            cwd = working_dir)

    #  if >= this percentage of bases are minor alleles, split candidate sequence
    em.snp_percentage_thresh = options.snp_fraction_thresh
    # if prob(N) for minor allele base N is >= this threshold, call site a minor allele
    em.snp_minor_prob_thresh = options.variant_fraction_thresh 
    em.min_depth = options.min_depth
    if options.nice_mapping is not None:
        em.mapping_nice = options.nice_mapping

    em.initialize_EM(options.mapping, options.fasta_db)

    # BEGIN ITERATIONS
    do_iterations(em, max_iter = options.iterations, save_every = options.save_every)

    return

if __name__ == '__main__':
    main()

