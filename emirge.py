#!/usr/bin/env python2
"""
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

for help, type:
python emirge.py --help
"""

USAGE = \
"""usage: %prog DIR <required_parameters> [options]

This version of EMIRGE (%prog) attempts to reconstruct rRNA SSU genes from
Illumina metagenomic data.
DIR is the working directory to process data in.
Use --help to see a list of required and optional arguments

Additional information:
https://groups.google.com/group/emirge-users
https://github.com/csmiller/EMIRGE/wiki

If you use EMIRGE in your work, please cite these manuscripts, as appropriate.

Miller CS, Baker BJ, Thomas BC, Singer SW, Banfield JF (2011)
EMIRGE: reconstruction of full-length ribosomal genes from microbial community short read sequencing data.
Genome biology 12: R44. doi:10.1186/gb-2011-12-5-r44.

Miller CS, Handley KM, Wrighton KC, Frischkorn KR, Thomas BC, Banfield JF (2013)
Short-Read Assembly of Full-Length 16S Amplicons Reveals Bacterial Diversity in Subsurface Sediments.
PloS one 8: e56018. doi:10.1371/journal.pone.0056018.
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
import _emirge

BOWTIE_l = 20
BOWTIE_e  = 300

BOWTIE_ASCII_OFFSET = 33   # currently, bowtie writes quals with an ascii offset of 33


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
          all characters are read as upper-case
    """
    if record is None:
        record = Record()

    for recordstring in re.split('\n>', filehandle.read()[1:]):
        record.title, record.sequence= recordstring.split('\n',1)
        record.sequence = record.sequence.replace('\n','').replace(' ','').upper()
        yield record

class EM(object):
    """
    driver class for EM algorithm
    """
    _VERBOSE = True
    base2i = {"A":0,"T":1,"C":2,"G":3}
    i2base = dict([(v,k) for k,v in base2i.iteritems()])
    # asciibase2i = {65:0,84:1,67:2,71:3}

    DEFAULT_ERROR = 0.05

    def __init__(self, reads1_filepath, reads2_filepath,
                 insert_mean,
                 insert_sd,
                 n_cpus = 1,
                 cwd = os.getcwd(), max_read_length = 76,
                 iterdir_prefix = "iter.", cluster_thresh = 0.97,
                 mapping_nice = None,
                 reads_ascii_offset = 64,
                 expected_coverage_thresh = 20):
        """

        n_cpus is how many processors to use for multithreaded steps (currently only the bowtie mapping)
        mapping_nice is nice value to add to mapping program
        """
        self.reads1_filepath = reads1_filepath
        self.reads2_filepath = reads2_filepath
        self.insert_mean = insert_mean
        self.insert_sd = insert_sd
        self.n_cpus = n_cpus
        self.mapping_nice   = mapping_nice
        self.reads_ascii_offset = reads_ascii_offset

        self.iteration_i = None         # keeps track of which iteration we are on.
        self.resume_i    = None         # special case for resuming iterations
        self.cwd = cwd
        self.max_read_length = max_read_length
        self.iterdir_prefix = iterdir_prefix
        self.cluster_thresh = cluster_thresh   # if two sequences evolve to be >= cluster_thresh identical (via usearch), then merge them. [0, 1.0]
        assert self.cluster_thresh >= 0 and self.cluster_thresh <= 1.0
        self.expected_coverage_thresh = expected_coverage_thresh

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

        self.quals = []
        self.reads = []
        self.readlengths = []

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
                self.reads
                self.quals
                self.readlengths

        creates a new EMPTY entry for (appends to list, removes t-2
                self.priors
                self.posteriors

        creates new each iteration (overwrites):
                self.likelihoods
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
            start_time = time()

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
            read_i = max([-1] + self.read_i2read_name[-2].keys()) + 1           # [-1] added for resume case

        # reset this every iteration
        self.coverage = [0]*seq_i

        # FIRST PASS through file just to get values out of file in memory.
        # no actual processing of data here.
        # Goal is to use multiprocess.Pool, but farming out parsing of samfile
        # bogs down with concurrent disk access from worker threads.

        # speed hacks
        bamfile = pysam.Samfile(bam_filename, "rb")
        getrname = bamfile.getrname
        quals = self.quals
        readlengths = self.readlengths
        reads = self.reads
        coverage = self.coverage
        ascii_offset = BOWTIE_ASCII_OFFSET
        fromstring = numpy.fromstring
        array = numpy.array
        uint8 = numpy.uint8
        base_alpha2int = _emirge.base_alpha2int

        # this pass just to read from file, get disk access over with.  As little processing as possible.
        predata = [(alignedread.pos, alignedread.tid, alignedread.is_read2, alignedread.qname,
                    alignedread.qual, alignedread.seq) for alignedread in bamfile]

        # cache bamfile information in this data structure, bamfile_data:
        # [seq_i, read_i, pair_i, rlen, pos], dtype=int))
        self.bamfile_data = numpy.empty((len(predata), 5), dtype=int)
        bamfile_data = self.bamfile_data # speed hack to avoid dot lookups

        # set here:
        #       self.sequence_name2sequence_i
        #       self.sequence_i2sequence_name
        #       self.read_i2read_name
        #       self.read_name2read_i
        #       bamfile_data
        #       self.reads
        #       self.readlengths
        #       self.quals

        seq_i, read_i = _emirge.process_bamfile_predata(seq_i, read_i,
                                                        predata, bamfile.references,
                                                        self.sequence_name2sequence_i[-1], self.sequence_i2sequence_name[-1],
                                                        self.read_name2read_i[-1], self.read_i2read_name[-1],
                                                        self.reads, self.quals, self.readlengths,
                                                        self.coverage, BOWTIE_ASCII_OFFSET,
                                                        self.bamfile_data)

        bamfile.close()

        self.priors.append(numpy.zeros(seq_i, dtype = numpy.float))
        self.likelihoods = sparse.coo_matrix((seq_i, read_i), dtype = numpy.float)  # init all to zero.
        self.posteriors.append(sparse.lil_matrix((seq_i+1, read_i+1), dtype=numpy.float))


        self.probN = [None for x in range(max(self.sequence_name2sequence_i[-1].values())+1)]
        self.unmapped_bases = [None for x in self.probN]
        self.mean_read_length = numpy.mean(self.readlengths)

        self.sequence_name2fasta_name = {}
        for record in FastIterator(file(self.current_reference_fasta_filename)):
            refname = record.title.split()[0]
            self.sequence_name2fasta_name[refname] = record.title.split()[0]
            seq_i = self.sequence_name2sequence_i[-1].get(refname)
            if seq_i is not None:
                self.probN[seq_i] = numpy.zeros((len(record.sequence), 5), dtype=numpy.float)   #ATCG[other] --> 01234
                self.coverage[seq_i] = self.coverage[seq_i] / float(len(record.sequence))

        for d in [self.priors, self.posteriors,
                  self.sequence_name2sequence_i, self.sequence_i2sequence_name,
                  self.read_name2read_i, self.read_i2read_name]:
            if len(d) > 2:
                trash = d.pop(0)  # no longer care about t-2
                del trash

        if self._VERBOSE:
            sys.stderr.write("DONE Reading bam file %s at %s [%s]...\n"%(bam_filename, ctime(), timedelta(seconds = time()-start_time)))
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
        # PRIOR
        for (seq_i, read_i, pair_i, rlen, pos) in self.bamfile_data:
            self.priors[-1][seq_i] += 1

        nonzero_indices = numpy.nonzero(self.priors[-1])  # only divide cells with at least one count.  Set all others to Pr(S) = 0
        self.priors[-1] = self.priors[-1][nonzero_indices] / self.priors[-1][nonzero_indices].sum()  # turn these into probabilities
        self.priors.append(self.priors[-1].copy())  # push this back to t-1 (index == -2)

        # write initial priors as special case:
        self.print_priors(os.path.join(self.cwd, "priors.initialized.txt"))

        if self._VERBOSE:
            sys.stderr.write("DONE with initialization at %s...\n"%(ctime()))
        return
    def resume(self, resume_from):
        """
        INPUT:   iteration number (integer) of previously completed EMIRGE iteration.
                 Expects completed iteration directory to contain:
                   probN.pkl.gz
                 and PREVIOUS iteration directory to contain:
                   priors.iter.NN.txt
                   iter.NN.cons.fasta
                   bowtie.iter.NN.PE.bam
        ON EXIT: ready for do_iteration to be called with next iteration.
        """
        self.resume_i = resume_from
        resume_iterdir = os.path.join(self.cwd, "%s%02d"%(self.iterdir_prefix, self.resume_i))
        previous_iterdir = os.path.join(self.cwd, "%s%02d"%(self.iterdir_prefix, self.resume_i - 1))
        if not os.path.exists(resume_iterdir):
            raise OSError, "\n\nERROR: Cannot resume from non-existent directory %s"%(resume_iterdir)
        if not os.path.exists(previous_iterdir):
            raise OSError, "\n\nERROR: Resume requires the previous iteration directory (%02d) also be present."%(resume_iterdir - 1)
        if not os.path.exists(os.path.join(resume_iterdir, "bowtie.iter.%02d.log.gz"%(self.resume_i))):
            raise OSError, "\n\nERROR: directory %s appears to be an incomplete iteration (no bowtie log file).  You must resume from a *completed* iteration (perhaps try iteration %02d)."%(resume_iterdir, self.resume_i-1)
        if self._VERBOSE:
            sys.stderr.write("Resuming EMIRGE from iteration %02d at %s ...\nStarting from information in directory:\n%s\n"%(self.resume_i,
                                                                                                                             ctime(),
                                                                                                                             resume_iterdir))
        # before calling read_bam in do_iteration, need to reset
        # sequence_name2sequence_i, sequence_i2sequence_name, priors 
        # from priors file (read_bam will add new entries for all otherwise)
        # FROM PREVIOUS ITER

        seq_i_seen = []
        priors_seen = []

        self.sequence_name2sequence_i = [{}, {}]
        self.sequence_i2sequence_name = [{}, {}]
        self.read_name2read_i         = [{}, {}]
        self.read_i2read_name         = [{}, {}]
        
        with open(os.path.join(previous_iterdir, "priors.iter.%02d.txt"%(self.resume_i - 1))) as f:
            for line in f:
                seq_i, seq_name, prior = line.split()
                seq_i = int(seq_i)
                prior = float(prior)

                for dict_i in [-1, -2]:
                    self.sequence_name2sequence_i[dict_i][seq_name] = seq_i
                    self.sequence_i2sequence_name[dict_i][seq_i]    = seq_name

                seq_i_seen.append(seq_i)
                priors_seen.append(prior)

        self.priors = []
        self.priors.append(numpy.zeros(max(seq_i_seen)+1, dtype=numpy.float))
        self.priors[-1][seq_i_seen] = priors_seen
        self.priors.append(self.priors[-1].copy())
        # note that self.priors gets appended with empty priors list upon read_bam being called again in initial do_iteration call

        # this is tricky.  We are actually wanting to re-calc likelihood and posteriors for CURRENT ITER based on priors from previous iter, but probN calced for this iter. Since do_iteration immediately increments, need to set this.
        self.iteration_i = self.resume_i - 1 

        if self._VERBOSE:
            sys.stderr.write("DONE with resume initialization at %s...\n"%(ctime()))
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
        if not os.path.exists(self.iterdir):
            os.mkdir(self.iterdir)
        # remove any previous fai files in this directory (e.g. from previous run we are resuming from)
        for filename in os.listdir(self.iterdir):
            if filename.endswith('.fai'):
                os.remove(os.path.join(self.iterdir, filename))
        
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
        if self._VERBOSE:
            sys.stderr.write("Writing priors and probN to disk for iteration %d at %s...\n"%(self.iteration_i, ctime()))
        self.print_priors()
        # python gzip.GzipFile is slow.  Use system call instead
        # cPickle.dump(self.probN, gzip.GzipFile(os.path.join(self.iterdir, 'probN.pkl.gz'), 'w'), cPickle.HIGHEST_PROTOCOL)
        pickled_filename = os.path.join(self.iterdir, 'probN.pkl')
        cPickle.dump(self.probN, file(pickled_filename, 'w'), cPickle.HIGHEST_PROTOCOL)
        check_call("gzip -f %s"%(pickled_filename), shell=True, stdout = sys.stdout, stderr = sys.stderr)
        if self._VERBOSE:
            sys.stderr.write("DONE Writing priors and probN to disk for iteration %d at %s...\n"%(self.iteration_i, ctime()))

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
        else:
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
            t0 = time()

        splitcount = 0
        cullcount  = 0
        of = file(outputfilename, 'w')

        times_split   = []              # DEBUG
        times_posteriors   = []              # DEBUG
        seqs_to_process = len(self.probN) # DEBUG

        i2base = self.i2base
        updated_seq_i = max(self.sequence_i2sequence_name[-1].keys()) + 1
        rows_to_add = []                # these are for updating posteriors at end with new minor strains
        cols_to_add = []
        data_to_add = []
        probNtoadd  = []  # for newly split out sequences

        self.posteriors[-1] = self.posteriors[-1].tolil()  # just to make sure this is in row-access-friendly format


        loop_t0 = time()
        for seq_i, probNarray in enumerate(self.probN):
            seq_i_t0 = time()
            if probNarray is None: # means this sequence is no longer present in this iteration
                continue
            # check if coverage passes self.min_depth, if not don't write it (culling happens here)
            if self.min_depth is not None and self.coverage[seq_i] < self.min_depth:  # could also do this only after self.iteration_i > 5 or something
                # could adjust priors and posteriors here, but because prior will already be low (b/c of low coverage)
                # and because next round will have 0 mappings (no sequence in reference file to map to), this seems
                # unneccesary.
                cullcount += 1
                self.probN[seq_i] = None
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

            # passes coverage thresholds
            title = self.sequence_i2sequence_name[-1][seq_i]
            consensus = numpy.array([i2base.get(x, "N") for x in numpy.argsort(probNarray)[:,-1]])

            # check for minor allele consensus, SPLIT sequence into two candidate sequences if passes thresholds.

            minor_indices = numpy.argwhere((probNarray >= self.snp_minor_prob_thresh).sum(axis=1) >= 2)[:,0]
            if minor_indices.shape[0] > 0:
                minor_fraction_avg = numpy.mean(probNarray[(minor_indices, numpy.argsort(probNarray[minor_indices])[:, -2])])
            else:
                minor_fraction_avg = 0.0
            # NEW rule: only split sequence if *expected* coverage
            # of newly split minor sequence (assuming uniform read
            # coverage over reconstructed sequence) is > some
            # threshold.  Here, expected coverage is calculated
            # based on:
            # Prior(seq_i) * number of reads * avg read length
            expected_coverage_minor = ( self.priors[-1][seq_i] * minor_fraction_avg * len(self.reads) * self.mean_read_length ) / probNarray.shape[0]
            expected_coverage_major = ( self.priors[-1][seq_i] * (1-minor_fraction_avg) * len(self.reads) * self.mean_read_length ) / probNarray.shape[0]
            # print >> sys.stderr, "DEBUG: ", seq_i, expected_coverage_minor, expected_coverage_major, minor_indices.shape, self.priors[-1][seq_i] , minor_fraction_avg , len(self.reads) , self.mean_read_length , probNarray.shape[0]

            if minor_indices.shape[0] / float(probNarray.shape[0]) >= self.snp_percentage_thresh and \
                   expected_coverage_minor >= self.expected_coverage_thresh:
                splitcount += 1
                if self._VERBOSE:
                    t0_split = time()
                major_fraction_avg = 1.-minor_fraction_avg # if there's >=3 alleles, major allele keeps prob of other minors)
                minor_bases   = numpy.array([i2base.get(x, "N") for x in numpy.argsort(probNarray[minor_indices])[:,-2]]) # -2 gets second most probably base
                minor_consensus = consensus.copy()               # get a copy of the consensus
                minor_consensus[minor_indices] = minor_bases     # replace the bases that pass minor threshold
                # now deal with naming.
                title_root = re.search(r'(.+)(_m(\d+))$', title)
                if title_root is None: # no _m00 on this name
                    title_root = title[:]
                else:
                    title_root = title_root.groups()[0]
                # now check for any known name with same root and a _m on it.
                previous_m_max = max([0] + [int(x) for x in re.findall(r'%s_m(\d+)'%re.escape(title_root), " ".join(self.sequence_i2sequence_name[-1].values()))])
                m_title = "%s_m%02d"%(title_root, previous_m_max+1)

                # also split out Priors and Posteriors (which will be used in next round), split with average ratio of major to minor alleles.
                # updating priors first:
                old_prior = self.priors[-1][seq_i]
                self.priors[-1][seq_i] = old_prior * major_fraction_avg
                seq_i_minor = updated_seq_i  # grow data structs by 1
                updated_seq_i += 1
                self.sequence_i2sequence_name[-1][seq_i_minor] = m_title
                self.sequence_name2sequence_i[-1][m_title] = seq_i_minor
                self.sequence_name2fasta_name[m_title] = m_title  # need to depreciate this data structure...
                # how I adjust probN here for newly split seq doesn't really matter,
                # as it is re-calculated next iter.
                # this only matters for probN.pkl.gz file left behind for this iteration.
                # for now just set prob(major base) = 0 and redistribute prob to other bases for minor,
                # and set prob(minor base) = 0 and redistribute prob to other bases for major
                # MINOR
                major_base_i = numpy.argsort(probNarray[minor_indices])[:, -1]
                newprobNarray = probNarray.copy()
                newprobNarray[(minor_indices, major_base_i)] = 0
                newprobNarray = newprobNarray / numpy.sum(newprobNarray, axis=1).reshape(newprobNarray.shape[0], 1)
                probNtoadd.append(newprobNarray)
                # MAJOR
                minor_base_i = numpy.argsort(probNarray[minor_indices])[:, -2]
                probNarray[(minor_indices, minor_base_i)] = 0
                probNarray = probNarray / numpy.sum(probNarray, axis=1).reshape(probNarray.shape[0], 1)

                new_priors = numpy.zeros(seq_i_minor+1, dtype=self.priors[-1].dtype)
                new_priors[:-1] = self.priors[-1].copy()
                new_priors[seq_i_minor] = old_prior * minor_fraction_avg
                trash = self.priors.pop()
                del trash
                self.priors.append(new_priors)


                # --- THIS WAS SLOW STEP ---
                # keep track of all new minor data to add and add it
                # once at end for ALL split sequences with one coo
                # matrix construction, instead of each iteration.

                t_posterior = time()

                # new_read_probs, new_rows, new_cols = adjust_posteriors_for_split(AAAA, BBBB, CCCC) # TODO: could move to Cython

                # updating posteriors. for each seq-read pair with prob > 0, split prob out to major and minor seq.
                new_cols = self.posteriors[-1].rows[seq_i] # col in coo format
                new_read_probs  = [x * minor_fraction_avg for x in self.posteriors[-1].data[seq_i]]  # data in coo format
                new_rows = [seq_i_minor for x in new_cols]  # row in coo format

                # add new read probs to cache of new read probs to add at end of loop
                rows_to_add.extend(new_rows)
                cols_to_add.extend(new_cols)
                data_to_add.extend(new_read_probs)

                # adjust old read probs to reflect major strain
                self.posteriors[-1].data[seq_i] = [x  * major_fraction_avg for x in self.posteriors[-1].data[seq_i]]

                times_posteriors.append(time() - t_posterior)
                # --- END SLOW STEP ---

                # adjust self.unmapped_bases (used in clustering).  For now give same pattern as parent
                self.unmapped_bases.append(self.unmapped_bases[seq_i].copy())

                # write out minor strain consensus
                of.write(">%s\n"%(m_title))
                of.write("%s\n"%("".join(minor_consensus)))
                if self._VERBOSE:
                    sys.stderr.write("splitting sequence %d (%s) to %d (%s)...\n"%(seq_i, title,
                                                                                   seq_i_minor, m_title))


                times_split.append(time()-seq_i_t0)

            # now write major strain consensus, regardless of whether there was a minor strain consensus
            of.write(">%s\n"%(title))
            of.write("%s\n"%("".join(consensus)))

        # END LOOP
        loop_t_total = time() - loop_t0
        # update posteriors matrix with newly added minor sequences new_posteriors via coo, then convert to csr.
        new_posteriors = self.posteriors[-1].tocoo()  # first make a copy in coo format
        # then create new coo matrix with new shape, appending new row, col, data to old row, col, data

        new_posteriors = sparse.coo_matrix((numpy.concatenate((new_posteriors.data, data_to_add)),
                                           (numpy.concatenate((new_posteriors.row, rows_to_add)),
                                            numpy.concatenate((new_posteriors.col, cols_to_add)))),
                                           shape=(updated_seq_i, self.posteriors[-1].shape[1]),
                                           dtype=new_posteriors.dtype).tocsr()

        # finally, exchange in this new matrix
        trash = self.posteriors.pop()
        del trash
        self.posteriors.append(new_posteriors)

        # update probN array:
        self.probN.extend(probNtoadd)

        if self._VERBOSE:
            total_time = time()-t0
            sys.stderr.write("\tSplit out %d new minor strain sequences.\n"%(splitcount))
            if splitcount > 0:
                sys.stderr.write("\tAverage time for split sequence: [%.6f seconds]\n"%numpy.mean(times_split))
                sys.stderr.write("\tAverage time for posterior update: [%.6f seconds]\n"%numpy.mean(times_posteriors))
            sys.stderr.write("\tAverage time for non-split sequences: [%.6f seconds]\n"%((loop_t_total - sum(times_split)) / (seqs_to_process - len(times_split))))
            sys.stderr.write("\tCulled %d sequences\n"%(cullcount))
            sys.stderr.write("DONE Writing consensus for iteration %d at %s [%s]...\n"%(self.iteration_i, ctime(), timedelta(seconds = total_time)))

        return

    def write_consensus_with_mask(self, reference_fastafilename, output_fastafilename, mask):
        """
        write a consensus sequence to output_fastafilename for each
        sequence in probN where unmapped bases are replaced with:
          mask == "hard"  --> N
          mask == "soft"  --> lowercase letters
        If masking with soft bases, use reference_fastafilename for bases to use for unmapped bases.
        this is useful prior to usearch clustering.

        OUT: number of sequences processed
        """
        n_seqs = 0
        # i2base = self.i2base
        i2base_get = self.i2base.get # for speed
        of = file(output_fastafilename, 'w')
        reference_fastafile = pysam.Fastafile(reference_fastafilename)

        for seq_i in range(len(self.probN)):
            if self.probN[seq_i] is None:
                continue

            title = self.sequence_i2sequence_name[-1][seq_i]
            consensus = numpy.array([i2base_get(x, "N") for x in numpy.argsort(self.probN[seq_i])[:,-1]])
            # consensus = numpy.array([i2base[x] for x in numpy.argsort(self.probN[seq_i])[:,-1]]) # -- crashing, think because some bases unmapped so argsort defaults to index 4 (unknown) if no bases mapped?
            orig_bases = numpy.array(reference_fastafile.fetch(self.sequence_name2fasta_name[self.sequence_i2sequence_name[-1][seq_i]]).lower(), dtype='c')
            # now replace consensus bases with no read support with N
            # unmapped_indices = numpy.where(self.unmapped_bases[seq_i] == 1)
            unmapped_indices = numpy.where(consensus == "N")
            if mask == "hard":
                consensus[unmapped_indices] = 'N'
            elif mask == "soft":
                for unmapped_i in unmapped_indices[0]:
                    consensus[unmapped_i] = orig_bases[unmapped_i] # return to original base if unmapped.
                # consensus[unmapped_indices] = [letter.lower() for letter in consensus[unmapped_indices]]
            else:
                raise ValueError, "Invalid value for mask: %s (choose one of {soft, hard}"%mask
            of.write(">%s\n"%(title))
            of.write("%s\n"%("".join(consensus)))
            n_seqs += 1
        of.close()
        return n_seqs

    def cluster_sequences(self, fastafilename):
        """
        uses Edgar's USEARCH to merge sequences above self.cluster_thresh %ID over the
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

        only supports usearch version > 6, as the command line substantially changed in this version.
        """
        if self._VERBOSE:
            sys.stderr.write("Clustering sequences for iteration %d at %s...\n"%(self.iteration_i, ctime()))
            sys.stderr.write("\tcluster threshold = %.3f\n"%(self.cluster_thresh))
            start_time = time()
        tocleanup = []                  # list of temporary files to remove after done
        # get posteriors ready for slicing (just prior to this call, is csr matrix?):
        self.posteriors[-1] = self.posteriors[-1].tolil()

        tmp_fastafilename = fastafilename + ".tmp.fasta"
        num_seqs = self.write_consensus_with_mask(fastafilename, tmp_fastafilename, mask="soft")
        tocleanup.append(tmp_fastafilename)
        tmp_fastafile = pysam.Fastafile(tmp_fastafilename)
        tocleanup.append("%s.fai"%(tmp_fastafilename))
        # do global alignments with USEARCH/UCLUST.
        # I don't use --cluster because it doesn't report alignments
        # usearch is fast but will sometimes miss things -- I've tried to tune params as best as I can.
        # Also, I use a lower %ID thresh than specified for joining because I really calculate %ID over *mapped* sequence positions.

        sens_string = "--maxaccepts 8 --maxrejects 256"
        uclust_id = 0.80
        algorithm="-usearch_global"
        # uclust_id = self.cluster_thresh - 0.05

        # if em.iteration_i > 10:
        # num_seqs = len([x for x in self.probN if x is not None])
        assert num_seqs == len([x for x in self.probN if x is not None])
        if num_seqs < 1000:
            sens_string = "--maxaccepts 16 --maxrejects 256"
        if num_seqs < 500:
            sens_string = "--maxaccepts 32 --maxrejects 256"
        if num_seqs < 150:
            algorithm="-search_global"
            sens_string = "--maxaccepts 0 --maxrejects 0"  # slower, but more sensitive.
        # if really few seqs, then no use not doing smith-waterman or needleman wunsch alignments
        if num_seqs < 50:
            algorithm="-search_global"
            sens_string = "-fulldp"

        # there is a bug in usearch threads that I can't figure out (fails with many threads).  Currently limiting to max 6 threads
        usearch_threads = min(6, self.n_cpus)
        cmd = "usearch %s %s --db %s --id %.3f -quicksort -query_cov 0.5 -target_cov 0.5 -strand plus --userout %s.us.txt --userfields query+target+id+caln+qlo+qhi+tlo+thi -threads %d %s"%\
              (algorithm,
               tmp_fastafilename, tmp_fastafilename,
               uclust_id,
               tmp_fastafilename,
               usearch_threads,
               sens_string)

        if self._VERBOSE:
            sys.stderr.write("usearch command was:\n%s\n"%(cmd))

        check_call(cmd, shell=True, stdout = sys.stdout, stderr = sys.stderr)
        # read clustering file to adjust Priors and Posteriors, summing merged reference sequences
        tocleanup.append("%s.us.txt"%tmp_fastafilename)

        nummerged = 0
        alnstring_pat = re.compile(r'(\d*)([MDI])')
        already_removed = set()  # seq_ids
        # this is a bit slow and almost certainly could be sped up with algorithmic improvements.
        times = []  # DEBUG
        for row in csv.reader(file("%s.us.txt"%tmp_fastafilename), delimiter='\t'):
            # each row an alignment in userout file
            t0 = time()
            # member == query
            member_name = row[0]
            seed_name = row[1]
            if member_name == seed_name:
                continue # usearch allows self-hits, which we don't care about
            member_seq_id = self.sequence_name2sequence_i[-1].get(member_name)
            seed_seq_id = self.sequence_name2sequence_i[-1].get(seed_name)
            if member_seq_id in already_removed or seed_seq_id in already_removed:
                continue

            # decide if these pass the cluster_thresh *over non-gapped, mapped columns*
            member_fasta_seq = tmp_fastafile.fetch(member_name).upper()
            seed_fasta_seq   = tmp_fastafile.fetch(seed_name).upper()
            member_unmapped = self.unmapped_bases[member_seq_id]  # unmapped positions (default prob)
            seed_unmapped = self.unmapped_bases[seed_seq_id]
            # query+target+id+caln+qlo+qhi+tlo+thi %s"%\
            #   0     1     2   3   4   5  6    7
            member_start = int(row[4]) - 1   # printed as 1-based by usearch now
            seed_start   = int(row[6]) - 1

            t0 = time()
            # print >> sys.stderr, "DEBUG", alnstring_pat.findall(row[3])
            aln_columns, matches = _emirge.count_cigar_aln(tmp_fastafile.fetch(seed_name).upper(),
                                                           tmp_fastafile.fetch(member_name).upper(),
                                                           self.unmapped_bases[seed_seq_id].astype(numpy.uint8),
                                                           self.unmapped_bases[member_seq_id].astype(numpy.uint8),
                                                           seed_start,
                                                           member_start,
                                                           alnstring_pat.findall(row[3]))
            ## print >> sys.stderr, "DEBUG: %.6e seconds"%(time()-t0)# timedelta(seconds = time()-t0)

            # if alignment is less that 500 bases, or identity over those 500+ bases is not above thresh, then continue
            if (aln_columns < 500) or ((float(matches) / aln_columns) < self.cluster_thresh):
                continue

            # DEBUG PRINT:
            if self._VERBOSE and num_seqs < 50:
                # print >> sys.stderr, row
                # print >> sys.stderr, "%s %s %s %s  --  %s, %s"%(member_seq_id, member_name, seed_seq_id, seed_name, float(matches), aln_columns)
                print >> sys.stderr, "\t\t%s|%s vs %s|%s %.3f over %s aligned columns"%(member_seq_id, member_name, seed_seq_id, seed_name,
                                                                float(matches) / aln_columns, aln_columns)


            # if above thresh, then first decide which sequence to keep, (one with higher prior probability).
            percent_id = (float(matches) / aln_columns) * 100.
            t0 = time()
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

            # now merge posteriors (all removed probs from remove_seq_id go to keep_seq_id).
            # self.posteriors[-1] at this point is lil_matrix
            # some manipulations of underlying sparse matrix data structures for efficiency here.
            # 1st, do addition in csr format (fast), convert to lil format, and store result in temporary array.
            new_row = (self.posteriors[-1].getrow(keep_seq_id).tocsr() + self.posteriors[-1].getrow(remove_seq_id).tocsr()).tolil() # NEW 4
            # then change linked lists directly in the posteriors data structure -- this is very fast
            self.posteriors[-1].data[keep_seq_id] = new_row.data[0] # NEW 4
            self.posteriors[-1].rows[keep_seq_id] = new_row.rows[0] # NEW 4
            # these two lines remove the row from the linked list (or rather, make them empty rows), essentially setting all elements to 0
            self.posteriors[-1].rows[remove_seq_id] = []
            self.posteriors[-1].data[remove_seq_id] = []

            # set self.probN[removed] to be None -- note that this doesn't really matter, except for
            # writing out probN.pkl.gz every iteration, as probN is recalculated from bam file
            # with each iteration
            self.probN[remove_seq_id] = None

            already_removed.add(remove_seq_id)
            nummerged += 1
            if self._VERBOSE:
                times.append(time()-t0)
                sys.stderr.write("\t...merging %d|%s into %d|%s (%.2f%% ID over %d columns) in %.3f seconds\n"%\
                                 (remove_seq_id, remove_name,
                                  keep_seq_id,   keep_name,
                                  percent_id, aln_columns,
                                  times[-1]))

        # if len(times) and self._VERBOSE:  # DEBUG
        #     sys.stderr.write("merges: %d\n"%(len(times)))
        #     sys.stderr.write("total time for all merges: %.3f seconds\n"%(numpy.sum(times)))
        #     sys.stderr.write("average time per merge: %.3f seconds\n"%(numpy.mean(times)))
        #     sys.stderr.write("min time per merge: %.3f seconds\n"%(numpy.min(times)))
        #     sys.stderr.write("max time per merge: %.3f seconds\n"%(numpy.max(times)))

        # write new fasta file with only new sequences
        if self._VERBOSE:
            sys.stderr.write("Writing new fasta file for iteration %d at %s...\n"%(self.iteration_i, ctime()))
        tmp_fastafile.close()
        tocleanup.append("%s.fai"%(fastafilename))  # this file will change!  So must remove index file.  pysam should check timestamps of these!
        recordstrings=""
        num_seqs = 0
        for record in FastIterator(file(fastafilename)): # read through file again, overwriting orig file if we keep the seq
            seqname = record.title.split()[0]  # strip off beginning cluster marks
            seq_id = self.sequence_name2sequence_i[-1].get(seqname)
            if seq_id not in already_removed:
                recordstrings += str(record)  # could do a better job here of actually "merging" a new consensus, rather than just keeping one or the other.
                num_seqs += 1
        outfile = file(fastafilename, 'w')
        outfile.write(recordstrings)
        outfile.close()

        for fn in tocleanup:  # quite important, actually, to remove old fai index files.
            os.remove(fn)

        if self._VERBOSE:
            sys.stderr.write("\tremoved %d sequences after merging\n"%(nummerged))
            sys.stderr.write("\tsequences remaining for iteration %02d: %d\n"%(self.iteration_i, num_seqs))
            sys.stderr.write("DONE Clustering sequences for iteration %d at %s [%s]...\n"%(self.iteration_i, ctime(), timedelta(seconds = time()-start_time)))

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
        cmd = "bowtie-build -o 3 %s %s > %s 2>&1"%(full_fasta_path , bowtie_index, bowtie_logfile) # -o 3 for speed? magnitude of speedup untested!
        if self._VERBOSE:
            sys.stderr.write("\tbowtie-build command:\n")
            sys.stderr.write("\t%s\n"%cmd)
        check_call(cmd, shell=True, stdout = sys.stdout, stderr = sys.stderr)

        # 2. run bowtie
        nicestring = ""
        if nice is not None:
            nicestring = "nice -n %d"%(nice)

        if self.reads1_filepath.endswith(".gz"):
            cat_cmd = "gzip -dc "
        else:
            cat_cmd = "cat "

        # these are used for single reads too.
        shared_bowtie_params = "--phred%d-quals -t -p %s  -n 3 -l %s -e %s  --best --strata --all --sam --chunkmbs 128"%(self.reads_ascii_offset, self.n_cpus, BOWTIE_l, BOWTIE_e)

        minins = max((self.insert_mean - 3*self.insert_sd), self.max_read_length)
        maxins = self.insert_mean + 3*self.insert_sd
        output_prefix = os.path.join(self.iterdir, "bowtie.iter.%02d"%(self.iteration_i))

        if self.reads2_filepath is not None:
            bowtie_command = "%s %s | %s bowtie %s --minins %d --maxins %d %s -1 - -2 %s | samtools view -b -S -F 0x0004 - > %s.PE.bam 2> %s "%(\
                cat_cmd,
                self.reads1_filepath,
                nicestring,
                shared_bowtie_params,
                minins, maxins,
                bowtie_index,
                self.reads2_filepath,
                output_prefix,
                bowtie_logfile)
        else: # single reads
            bowtie_command = "%s %s | %s bowtie %s %s - | samtools view -b -S -F 0x0004 - > %s.PE.bam 2> %s "%(\
                cat_cmd,
                self.reads1_filepath,
                nicestring,
                shared_bowtie_params,
                bowtie_index,
                output_prefix,
                bowtie_logfile)

        if self._VERBOSE:
            sys.stderr.write("\tbowtie command:\n")
            sys.stderr.write("\t%s\n"%bowtie_command)

        check_call(bowtie_command, shell=True, stdout = sys.stdout, stderr = sys.stderr)

        if self._VERBOSE:
            sys.stderr.write("\tFinished Bowtie for iteration %02d at %s:\n"%(self.iteration_i, ctime()))

        # 3. clean up
        # check_call("samtools index %s.sort.PE.bam"%(output_prefix), shell=True, stdout = sys.stdout, stderr = sys.stderr)
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
            start_time = time()
        # first calculate self.probN from mapped reads, previous round's posteriors
        self.calc_probN()   # (handles initial iteration differently within this method)

        # for speed:
        probN = self.probN
        numpy_float = numpy.float
        sequence_name2sequence_i = self.sequence_name2sequence_i[-1]
        read_name2read_i         = self.read_name2read_i[-1]
        base2i_get = self.base2i.get
        arange = numpy.arange
        numpy_log = numpy.log
        e = numpy.e
        lik_row_seqi = numpy.empty(len(self.bamfile_data), dtype=numpy.uint) # these for constructing coo_matrix for likelihood.
        lik_col_readi = numpy.empty_like(lik_row_seqi)
        lik_data = numpy.empty(len(self.bamfile_data), dtype=numpy.float)
        bamfile_data = self.bamfile_data
        reads = self.reads
        quals = self.quals
        zeros = numpy.zeros

        # TODO: multiprocess -- worker pools would just return three arrays that get appended at end (extend instead of append) before coo matrix construction
        # data needed outside of bamfile:
        # optional: bamfile_data, quals
        # required: base2i, probN

        # keep looping here in python so that we can (in the future) use multiprocessing with an iterator (no cython looping).
        # calc_likelihood_cell = _emirge.calc_likelihood_cell
        # for alignedread_i, (seq_i, read_i, pair_i, rlen, pos) in enumerate(bamfile_data):
        #     lik_row_seqi[alignedread_i] = seq_i
        #     lik_col_readi[alignedread_i] = read_i
        #     lik_data[alignedread_i] = calc_likelihood_cell(seq_i, read_i, pair_i,
        #                                                    pos,
        #                                                    reads[read_i],
        #                                                    quals[read_i],
        #                                                    probN[seq_i])

        # Move looping to cython for small speed gain if we don't multiprocess.
        _emirge._calc_likelihood(bamfile_data,
                                 reads,
                                 quals,
                                 probN,
                                 lik_row_seqi,
                                 lik_col_readi,
                                 lik_data)

        # now actually construct sparse matrix.
        self.likelihoods = sparse.coo_matrix((lik_data, (lik_row_seqi, lik_col_readi)), self.likelihoods.shape, dtype=self.likelihoods.dtype).tocsr()
        if self._VERBOSE:
            sys.stderr.write("DONE Calculating likelihood for iteration %d at %s [%s]...\n"%(self.iteration_i, ctime(), timedelta(seconds = time()-start_time)))
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
            sys.stderr.write("DONE Calculating posteriors for iteration %d at %s [%.3f seconds]...\n"%(self.iteration_i, ctime(), time() - t_start))

    def calc_probN(self):
        """
        Pr(N=n)
        If read or sequence is new this round (not seen at t-1), then there is no Pr(S|R) from previous round,
        so we substitute Pr(S), the unbiased prior.
        If initial iteration, all reads and seqs are new, so all calcs for Pr(N=n) use the prior as weighting
        factor instead of previous round's posterior.
        """
        if self._VERBOSE:
            sys.stderr.write("\tCalculating Pr(N=n) for iteration %d at %s...\n"%(self.iteration_i, ctime()))
            start_time = time()

        # basecounts = [seqprobNarray.astype(numpy.uint32) for seqprobNarray in self.probN]
        initial_iteration = self.iteration_i < 1

        # for speed:
        probN = self.probN
        if initial_iteration or (self.resume_i is not None and self.iteration_i == self.resume_i):
            posteriors = None
        else:
            self.posteriors[-2] = self.posteriors[-2].tolil()
            posteriors = self.posteriors[-2]  # this depends on PREVIOUS iteration's posteriors (seq_n x read_n)

        priors     = self.priors[-2]          # and if no posteriors are present (or initial iter), falls back on priors from previous round
        base2i_get = self.base2i.get
        sequence_name2sequence_i = self.sequence_name2sequence_i[-1]
        read_name2read_i = self.read_name2read_i[-1]
        bamfile_data = self.bamfile_data
        reads = self.reads
        quals = self.quals
        np_int = numpy.int

        # could multithread this too.  self.probN is what is being written to.
        # calc_probN_read = _emirge.calc_probN_read
        # for alignedread_i in range(bamfile_data.shape[0]):
        #     seq_i, read_i, pair_i, rlen, pos = bamfile_data[alignedread_i]
        #     calc_probN_read(initial_iteration,
        #                     seq_i, read_i, pos,
        #                     priors,
        #                     posteriors,
        #                     reads[read_i],
        #                     quals[read_i],
        #                     probN[seq_i])

        if self.resume_i is not None and (self.iteration_i == self.resume_i):  # we have just resumed a run
            pickled_filename = os.path.join(self.iterdir, 'probN.pkl')
            sys.stderr.write("\tLoading probN for resume case from %s\n"%pickled_filename)
            check_call("gzip -d -c %s > %s"%(pickled_filename+'.gz', pickled_filename), shell=True, stdout = sys.stdout, stderr = sys.stderr)
            self.probN = cPickle.load(open(pickled_filename))     # and need to use already-calculated probN
            os.remove(pickled_filename) # git rid of decompressed file
        else:
            # here do looping in Cython (this loop is about 95% of the time in this method on test data):
            _emirge._calc_probN(self.bamfile_data,
                                initial_iteration,
                                priors,
                                posteriors,
                                self.reads,
                                self.quals,
                                self.probN)

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
                bases = numpy.array(self.fastafile.fetch(fastaname).upper(), dtype='c')[zero_indices[0]]
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
            sys.stderr.write("\tDONE calculating Pr(N=n) for iteration %d at %s [%s]...\n"%(self.iteration_i, ctime(), timedelta(seconds = time()-start_time)))

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
            this_seq_i, read_i, pair_i, rlen, pos = bamfile_data[alignedread_i]
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
            this_seq_i, read_i, pair_i, rlen, pos = bamfile_data[alignedread_i]
            if this_seq_i != seq_i or posteriors[seq_i, read_i] < 0.5:
                continue
            # means we have a read with more than 50% prob assigned to this sequence
            of_sam.write(alignedread)
            reads += 1
        of_sam.close()
        of_fasta = file('%s.ref.fasta'%(output_prefix), 'w')
        fasta_seq   = self.fastafile.fetch(self.sequence_i2sequence_name[-1][seq_i]).upper()
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
        raise NotImplementedError, "Broken with most recent revision of data structures (no more bamfile_readnames)"
        if output_prefix is None:
            output_prefix = os.path.join(em.iterdir, "%s"%seq_i)

        # go through bam file instead of original sequencing reads, as aligned reads only a fraction
        # of original file.
        self.posteriors[-1] = self.posteriors[-1].tolil()  # seq_i x read_i
        posteriors = self.posteriors[-1]  # seq_i x read_i
        of_fasta_name = '%s.reads.fasta'%(output_prefix)
        of_fasta = file(of_fasta_name, 'w')
        bamfile_data = self.bamfile_data
        reads = 0
        for alignedread_i, (this_seq_i, read_i, pair_i, rlen, pos) in enumerate(self.bamfile_data):
            if this_seq_i != seq_i or posteriors[seq_i, read_i] < 0.5:
                continue
            # means we have a read with more than 50% prob assigned to this sequence
            header = "%s DIRECTION: fwd CHEM: unknown TEMPLATE: %s"%(bamfile_readnames[alignedread_i],
                                                                     bamfile_readnames[alignedread_i])
            # header = "%s"%(reads)
            of_fasta.write(">%s\n%s\n"%(header, self.reads[alignedread_i, :rlen]))
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

def do_iterations(em, max_iter, save_every):
    """
    an EM object is passed in, so that one could in theory start from a saved state
    """
    os.chdir(em.cwd)

    if em.iteration_i < 0:  # first run
        em.do_iteration(em.current_bam_filename, em.current_reference_fasta_filename)

    while em.iteration_i < max_iter:
        subdir = os.path.join(em.cwd, "iter.%02d"%(em.iteration_i))
        em.do_iteration(os.path.join(subdir, "bowtie.iter.%02d.PE.bam"%(em.iteration_i)),
                        os.path.join(subdir, "iter.%02d.cons.fasta"%(em.iteration_i)))
        if save_every is not None and em.iteration_i > 0 and (em.iteration_i % save_every == 0):
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

    minins = max((options.insert_mean - 3*options.insert_stddev), options.max_read_length)
    maxins = options.insert_mean + 3*options.insert_stddev
    bampath_prefix = os.path.join(initial_mapping_dir, "initial_bowtie_mapping.PE")

    nicestring = ""
    if options.nice_mapping is not None:
        nicestring = "nice -n %d"%(options.nice_mapping)  # TODO: fix this so it isn't such a hack and will work in bash.  Need to rewrite all subprocess code, really (shell=False)
    reads_ascii_offset = {False: 64, True: 33}[options.phred33]
    if options.fastq_reads_1.endswith(".gz"):
        option_strings = ["gzip -dc "]
    else:
        option_strings = ["cat "]
    # shared regardless of whether paired mapping or not
    option_strings.extend([options.fastq_reads_1, nicestring, reads_ascii_offset, options.processors, BOWTIE_l, BOWTIE_e])

    # PAIRED END MAPPING
    if options.fastq_reads_2 is not None:
        option_strings.extend([minins, maxins, options.bowtie_db, options.fastq_reads_2, bampath_prefix])
        cmd = "%s %s | %s bowtie --phred%d-quals -t -p %s -n 3 -l %s -e %s --best --sam --chunkmbs 128 --minins %s --maxins %s %s -1 - -2 %s | samtools view -b -S -u -F 0x0004 - > %s.bam "%tuple(option_strings)    
    # SINGLE END MAPPING
    else:
        option_strings.extend([options.bowtie_db, bampath_prefix])
        cmd = "%s %s | %s bowtie --phred%d-quals -t -p %s -n 3 -l %s -e %s --best --sam --chunkmbs 128  %s - | samtools view -b -S -u -F 0x0004 - > %s.bam "%tuple(option_strings)    

    print "Performing initial mapping with command:\n%s"%cmd
    check_call(cmd, shell=True, stdout = sys.stdout, stderr = sys.stderr)
    sys.stdout.flush()
    sys.stderr.flush()

    return bampath_prefix+".bam"

def dependency_check():
    """
    check presense, versions of programs used in emirge
    TODO: right now just checking usearch, as the command line params
    and behavior are finicky and seem to change from version to
    version
    """
    # usearch
    working_maj = 6
    working_minor = 0
    working_minor_minor = 203
    match = re.search(r'usearch([^ ])* v([0-9]*)\.([0-9]*)\.([0-9]*)', Popen("usearch --version", shell=True, stdout=PIPE).stdout.read())
    if match is None:
        print >> sys.stderr, "FATAL: usearch not found in path!"
        exit(0)
    binary_name, usearch_major, usearch_minor, usearch_minor_minor = match.groups()
    usearch_major = int(usearch_major)
    usearch_minor = int(usearch_minor)
    usearch_minor_minor = int(usearch_minor_minor)
    if usearch_major < working_maj or \
       (usearch_major == working_maj and (usearch_minor < working_minor or \
                                          (usearch_minor == working_minor and usearch_minor_minor < working_minor_minor))):
        print >> sys.stderr, "FATAL: usearch version found was %s.%s.%s.\nemirge works with version >=  %s.%s.%s\nusearch has different command line arguments and minor bugs in previous versions that can cause problems."%(usearch_major, usearch_minor, usearch_minor_minor, working_maj, working_minor, working_minor_minor)
        exit(0)
    return

def main(argv = sys.argv[1:]):
    """
    command line interface to emirge

    """
    dependency_check()

    parser = OptionParser(USAGE)

    # REQUIRED
    group_reqd = OptionGroup(parser, "Required flags",
                             "These flags are all required to run EMIRGE, and may be supplied in any order.")

    group_reqd.add_option("-1", dest="fastq_reads_1", metavar="reads_1.fastq[.gz]",
                      type="string",
                      help="path to fastq file with \\1 (forward) reads from paired-end sequencing run, or all reads from single-end sequencing run.  File may optionally be gzipped.  EMIRGE expects ASCII-offset of 64 for quality scores.  (Note that running EMIRGE with single-end reads is largely untested.  Please let me know how it works for you.)")
    group_reqd.add_option("-f", "--fasta_db",
                      type="string",
                      help="path to fasta file of candidate SSU sequences")
    group_reqd.add_option("-b", "--bowtie_db",
                      type="string",
                      help="precomputed bowtie index of candidate SSU sequences (path to appropriate prefix; see --fasta_db)")
    group_reqd.add_option("-l", "--max_read_length",
                      type="int",
                      help="""length of longest read in input data.""")
    parser.add_option_group(group_reqd)

    # REQUIRED for paired end
    group_reqd_PE = OptionGroup(parser, "Required flags for paired-end reads",
                             "These flags are required to run EMIRGE when you have paired-end reads (the standard way of running EMIRGE), and may be supplied in any order.")
    group_reqd_PE.add_option("-2", dest="fastq_reads_2", metavar="reads_2.fastq",
                      type="string",
                      help="path to fastq file with \\2 (reverse) reads from paired-end run.  File must be unzipped for mapper.  EMIRGE expects ASCII-offset of 64 for quality scores.")
    group_reqd_PE.add_option("-i", "--insert_mean",
                      type="int", default=0,
                      help="insert size distribution mean.")
    group_reqd_PE.add_option("-s", "--insert_stddev",
                      type="int", default=0,
                      help="insert size distribution standard deviation.")
    parser.add_option_group(group_reqd_PE)

    # OPTIONAL
    group_opt = OptionGroup(parser, "Optional parameters",
                             "Defaults should normally be fine for these options in order to run EMIRGE")
    group_opt.add_option("-n", "--iterations",
                      type="int", default=40,
                      help="""Number of iterations to perform.  It may be necessary to use more iterations for more complex samples (default=%default)""")
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
                      help="If two candidate sequences share >= this fractional identity over their bases with mapped reads, then merge the two sequences into one for the next iteration.  (default: %default; valid range: [0.95, 1.0] ) ")
    group_opt.add_option("-c", "--min_depth",
                      type="float",
                      default = 3,
                      help = "minimum average read depth below which a candidate sequence is discarded for next iteration(default: %default)")
    group_opt.add_option("--nice_mapping",
                      type="int",
                      help="""If set, during mapping phase, the mapper will be "niced" by the Linux kernel with this value (default: no nice)""")
    group_opt.add_option("--phred33",
                         action="store_true", default=False,
                         help="Illumina quality values in fastq files are the (fastq standard) ascii offset of Phred+33.  This is the new default for Illumina pipeline >= 1.8. DEFAULT is still to assume that quality scores are Phred+64")
    group_opt.add_option("-e", "--save_every",
                      type="int", default=None,
                      help="""every SAVE_EVERY iterations, save some information about the program's state.  This is solely for debugging information, and is NOT required to resume a run (see --resume_from below).  (default=%default)""")

    parser.add_option_group(group_opt)
    # # RESUME
    group_resume = OptionGroup(parser, "Resuming iterations",
                               "These options allow you to resume iterations from a previously completed EMIRGE iteration.  This requires that directories for the iteration to resume from and the previous iteration both be present.  It is STRONGLY recommended that other options set on the command line be identical to the original run.  Note that EMIRGE does not check this for you!")
    group_resume.add_option("-r", "--resume_from",
                            type="int",
                            help="Resume iterations from COMPLETED iteration specified.  Requires that the iteration and previous iteration fully completed, i.e. a priors file, bam file, and fasta file are all present in the iteration directory.")
    parser.add_option_group(group_resume)

    # ACTUALLY PARSE ARGS
    (options, args) = parser.parse_args(argv)

    if len(args) !=1:
        parser.error("DIR is required, and all options except DIR should have a flag associated with them (options without flags: %s)"%args)
    if options.join_threshold < 0.95 or options.join_threshold > 1:
        parser.error("join_threshold must be between [0.95, 1.0].  You supplied %.3f"%options.join_threshold)

    working_dir = os.path.abspath(args[0])

    sys.stdout.write("""If you use EMIRGE in your work, please cite these manuscripts, as appropriate.

Miller CS, Baker BJ, Thomas BC, Singer SW, Banfield JF (2011)
EMIRGE: reconstruction of full-length ribosomal genes from microbial community short read sequencing data.
Genome biology 12: R44. doi:10.1186/gb-2011-12-5-r44.

Miller CS, Handley KM, Wrighton KC, Frischkorn KR, Thomas BC, Banfield JF (2013)
Short-Read Assembly of Full-Length 16S Amplicons Reveals Bacterial Diversity in Subsurface Sediments.
PloS one 8: e56018. doi:10.1371/journal.pone.0056018.\n\n""")

    sys.stdout.write("imported _emirge C functions from: %s\n"%(_emirge.__file__))
    sys.stdout.write("Command:\n")
    sys.stdout.write(' '.join([__file__]+argv))
    sys.stdout.write('\n\n')
    total_start_time = time()
    sys.stdout.write("EMIRGE started at %s\n"%(ctime()))
    sys.stdout.flush()

    if options.bowtie_db is None:
        if options.fasta_db:
            parser.error("Bowtie index is missing (--bowtie_db). You need to build it before running EMIRGE\nTry:\n\nbowtie-build %s bowtie_prefix" % options.fasta_db)
        else:
            parser.error("Bowtie index is missing (--bowtie_db). You need to build it before running EMIRGE\nTry:\n\nbowtie-build candidate_db.fasta bowtie_prefix")
    if options.fasta_db is None:
        parser.error("Fasta file for candidate database is missing. Specify --fasta_db. (try --help for more information)")

    if any((val is None) for val in [options.fastq_reads_1, options.max_read_length]): 
        parser.error("Some required arguments are missing (try --help)")

    if options.fastq_reads_2 is not None and any((val is 0) for val in [options.insert_mean, options.insert_stddev]): 
        parser.error("Some required arguments for paired-end mode are missing (try --help)")

    if options.fastq_reads_2 is not None and options.fastq_reads_2.endswith('.gz'):
        parser.error("Read 2 file cannot be gzipped (see --help)")

    if not os.path.exists(working_dir):
        os.mkdir(working_dir)

    # clean up options to be absolute paths
    for o in ["fastq_reads_1", "fastq_reads_2", "fasta_db", "bowtie_db", "mapping"]:
        current_o_value = getattr(options, o)
        if current_o_value is not None:
            setattr(options, o, os.path.abspath(current_o_value))

    # DO INITIAL MAPPING if not provided with --mapping
    if options.mapping is None and options.resume_from is None:
        options.mapping = do_initial_mapping(working_dir, options)

    # finally, CREATE EM OBJECT
    em = EM(reads1_filepath = options.fastq_reads_1,
            reads2_filepath = options.fastq_reads_2,
            insert_mean = options.insert_mean,
            insert_sd   = options.insert_stddev,
            max_read_length = options.max_read_length,
            cluster_thresh = options.join_threshold,
            n_cpus = options.processors,
            cwd = working_dir,
            reads_ascii_offset = {False: 64, True: 33}[options.phred33])

    #  if >= this percentage of bases are minor alleles, split candidate sequence
    em.snp_percentage_thresh = options.snp_fraction_thresh
    # if prob(N) for minor allele base N is >= this threshold, call site a minor allele
    em.snp_minor_prob_thresh = options.variant_fraction_thresh
    em.min_depth = options.min_depth
    if options.nice_mapping is not None:
        em.mapping_nice = options.nice_mapping

    if options.resume_from is None:                             # new EMIRGE run
        em.initialize_EM(options.mapping, options.fasta_db)
    else:                                                       # resumed EMIRGE run
        if options.resume_from < 1:
            parser.error("--resume_from must be >= 1")
        em.resume(options.resume_from)                          # This **assumes** that cmd line args are the same!

    # BEGIN ITERATIONS
    do_iterations(em, max_iter = options.iterations, save_every = options.save_every)

    sys.stdout.write("EMIRGE finished at %s.  Total time: %s\n"%(ctime(), timedelta(seconds = time()-total_start_time)))

    return



if __name__ == '__main__':
    main()

