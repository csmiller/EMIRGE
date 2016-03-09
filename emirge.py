#!/usr/bin/env python2
"""
EMIRGE:
Expectation-Maximization Iterative Reconstruction of Genes from the Environment

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

for help, type:
python emirge_amplicon.py --help
"""
import cPickle
import csv
import multiprocessing
import os
import re
import sys
from datetime import timedelta
from optparse import OptionParser, OptionGroup, SUPPRESS_HELP
from subprocess import Popen, PIPE, check_call
from time import ctime, time

import Emirge.amplicon as amplicon
import numpy as np
import pysam
from Bio import SeqIO
from scipy import sparse

from Emirge import io, log, mapping
from Emirge.clustering import vsearch
from Emirge.io import make_pipe
from Emirge.log import DEBUG, INFO
from emirge_rename_fasta import rename

USAGE = """usage: %prog DIR <required_parameters> [options]

This version of EMIRGE (%prog) attempts to reconstruct rRNA SSU genes
from Illumina amplicon data.  It can handle up to a few million rRNA
reads at a time.
DIR is the working directory to process data in.
Use --help to see a list of required and optional arguments

Additional information:
https://groups.google.com/group/emirge-users
https://github.com/csmiller/EMIRGE/wiki

If you use EMIRGE in your work, please cite these manuscripts, as appropriate.

Miller CS, Baker BJ, Thomas BC, Singer SW, Banfield JF (2011)
EMIRGE: reconstruction of full-length ribosomal genes from microbial community
short read sequencing data.
Genome biology 12: R44. doi:10.1186/gb-2011-12-5-r44.

Miller CS, Handley KM, Wrighton KC, Frischkorn KR, Thomas BC,
Banfield JF (2013)
Short-Read Assembly of Full-Length 16S Amplicons Reveals Bacterial Diversity in
Subsurface Sediments.
PloS one 8: e56018. doi:10.1371/journal.pone.0056018.
"""


BOWTIE_l = 20
BOWTIE_e = 300

# currently, bowtie writes quals with an ascii offset of 33
BOWTIE_ASCII_OFFSET = 33


class EM(object):
    """
    driver class for EM algorithm

    Props:
        snp_minor_prob_thresh:
            if prob(N) for minor allele base N is >= this threshold, call
            site a minor allele
        snp_percentage_thresh:
            if >= this percentage of bases are minor alleles (according to
            self.snp_minor_prob_thresh), then split this sequence into
            two sequences.
        reads:
            holds the bases for all reads.
            uint8 ndarray with 3 dimensions: read number, pair (0,1), base
            allocated to have n_reads x 2 x max_read_length size
        quals:
            see reads, this array holds the base qualities
        readlengths:
            holds the length of each read
            uint16 ndarray with 2 dimensions: read number, pair (0,1)
        sequence_name2sequence_i:
        sequence_i2sequence_name:
        refseq_lengths:
    """
    _VERBOSE = True
    base2i = {"A": 0, "T": 1, "C": 2, "G": 3, "D": 4}
    i2base = dict([(v, k) for k, v in base2i.iteritems()])
    # asciibase2i = {65:0,84:1,67:2,71:3}

    DEFAULT_ERROR = 0.05

    def __init__(self,
                 mapper,
                 candidate_db,
                 reads1_filepath,
                 reads2_filepath,
                 insert_mean,
                 insert_sd,
                 output_files_prefix,
                 cov_thresh,
                 snp_percentage_thresh,
                 snp_minor_prob_thresh,
                 min_length_coverage,
                 n_cpus=1,
                 cwd=os.getcwd(),
                 max_read_length=76,
                 iterdir_prefix="iter.",
                 cluster_thresh=0.97,
                 indel_thresh=0.3,
                 mapping_nice=None,
                 reads_ascii_offset=64,
                 expected_coverage_thresh=10
                 ):
        """
        n_cpus is how many processors to use for multithreaded steps
        (currently only the bowtie mapping)
        mapping_nice is nice value to add to mapping program
        """

        assert 0 <= cluster_thresh <= 1.0

        self.mapper = mapper
        self.reads1_filepath = reads1_filepath
        self.reads2_filepath = reads2_filepath
        self.insert_mean = insert_mean
        self.insert_sd = insert_sd
        self.output_files_prefix = output_files_prefix
        self.n_cpus = n_cpus
        self.mapping_nice = mapping_nice
        self.reads_ascii_offset = reads_ascii_offset
        self.cwd = cwd
        self.iterdir_prefix = iterdir_prefix
        self.max_read_length = max_read_length
        self.cluster_thresh = cluster_thresh
        # if two sequences evolve to be >= cluster_thresh identical
        # (via vsearch), then merge them. [0, 1.0]
        self.expected_coverage_thresh = expected_coverage_thresh
        self.indel_thresh = indel_thresh
        self.snp_minor_prob_thresh = snp_minor_prob_thresh
        self.snp_percentage_thresh = snp_percentage_thresh
        self.min_length_coverage = min_length_coverage
        # EXPERIMENTAL.  Fraction of length that has to be covered by
        #                >= min_length_cov_depth
        self.current_reference_fasta_filename = candidate_db

        self.max_iterations = None
        self.bamfile_data = None
        self.fragments_mapped = None
        self.current_bam_filename = None
        self.n_alignments = None
        self.fastafile = None
        self.cigars = None
        self.probN = None
        self.prob_indels = None
        self.unmapped_bases = None
        self.mean_read_length = None
        self.iterdir = None
        self.initial_bam_filename_to_remove = None
        self.initial_compress_process = None
        self.avg_emirge_seq_length = None
        self.prob_min = None

        # keeps track of which iteration we are on.
        self.iteration_i = None

        # Single numpy array.  Has the shape: (numsequences x numreads)
        # [numreads can be numpairs]
        self.likelihoods = None
        # = Pr(R_i|S_i),
        # the likelihood of generating read R_i given sequence S_i

        # list of numpy arrays.  list index is iteration number.
        # Each numpy array has the shape: (numsequences,)
        self.priors = []
        # = Pr(S_i),
        # the prior probability that sequence S generated any read

        # list of numpy arrays.  list index is iteration number.
        # Each numpy array has the shape: (numsequences x numreads)
        self.posteriors = []
        # = Pr(S_i|R_i),
        # the posterior probability that sequence S_i generated read R_i

        # dict's and list keeping id mappings between sequence names
        # and internal indices (seq_i)
        # index is stable between iterations.
        # If sequence_i2sequence value is None, means this sequence was
        # abandoned in a previous round
        self.sequence_name2sequence_i = {}
        # list index is iteration number:
        self.sequence_i2sequence_name = []
        self.refseq_lengths = []

        # seq_i --> iteration first seen.  Useful for keeping track of
        # when a sequence first appeared, and not allowing merging of
        # recently split out sequences
        self.split_seq_first_appeared = {}

        # in fastq input (number of reads **or number of read pairs**)
        self.n_reads = mapper.n_reads
        self.n_reads_mapped = 0
        self.n_sequences = 0

        self.num_seqs = 0  # this is used in sequence clustering - might need
        # to rename it for clarity

        # other constants, potentially changeable, tunable later,
        # or could incorporate into probabilistic model:
        #
        # minimum depth to keep sequence around for next round
        self.min_depth = 5.0
        # minimum prior probability for a sequence to keep it
        # around for next round (alternative to depth, which is
        # a little weird when you allow mappings to more than one place.
        # NOT YET IMPLEMENTED
        self.min_prior = None

        # list of numpy arrays -- per base coverage values.
        self.base_coverages = []
        # EXPERIMENTAL:  Minimum coverage in order to be counted in
        #                min_length_coverage
        self.min_length_coverage_def = 1

        # cov_thresh=options.min_coverage_threshold
        self.cov_thresh = cov_thresh

        if self.reads2_filepath is None:
            self.paired_end = False
        else:
            self.paired_end = True

        self.temporary_files = []  # to remove at end of run

        log.info("Number of reads (or read pairs) in input file(s): %d"
                 % self.n_reads)

        # bool vector indicating whether read n was ever seen mapped
        self.reads_seen = np.zeros(self.n_reads, dtype=np.uint8)

        # where 1st dimension is read index (from rewritten file headers)
        # and second dimension is read number (0 or 1 ==> read /1 or read /2)
        # 3rd dimension for reads and quals is max_readlen
        self.reads = np.empty((self.n_reads, 2, self.max_read_length),
                              dtype=np.uint8)
        self.quals = np.empty_like(self.reads)
        self.readlengths = np.empty((self.n_reads, 2), dtype=np.uint16)

        # read through reads file again, fill these.
        amplicon.populate_reads_arrays(self)

    def do_initial_mapping(self):
        workdir = os.path.join(self.cwd, "initial_mapping")
        os.mkdir(workdir)

        self.fragments_mapped, self.current_bam_filename = \
            self.mapper.map_reads(self.current_reference_fasta_filename,
                                  workdir)

    def do_mapping(self, filename):
        self.current_reference_fasta_filename = filename
        self.fragments_mapped, self.current_bam_filename = \
            self.mapper.map_reads(self.current_reference_fasta_filename,
                                  self.iterdir)

    @log.timed("Reading bam file")
    def read_bam(self):
        """
        reads a bam file and...
        updates:
                self.sequence_i2sequence_name   # a numpy array
                self.sequence_name2sequence_i   # a dict
                self.read_name2read_i           # a dict
                self.probN
                self.refseq_lengths

        doesn't do anything with these anymore, they should be populated and
        stable with _emirge.populate_reads_arrays
                self.reads
                self.quals
                self.readlengths

        creates a new EMPTY entry for (appends to list, removes t-2
                self.priors
                self.posteriors

        creates new each iteration (overwrites):
                self.likelihoods
                self.unmapped_bases
                self.coverage
                self.bamfile_data
                self.cigars

        This MUST maintain seq_i to name and read_i to name mappings between
        iterations, so that a single name always maintains the same index from
        one iteration to the next.  One result of this requirement is that the
        various matrices can always get larger in a later t, but never smaller
        (as reads or seqs are added)
        """

        self.n_alignments = self.get_n_alignments_from_bowtie()
        self.fastafile = pysam.Fastafile(self.current_reference_fasta_filename)
        self.cigars = []  # list of pysam cigartuples
        # populate in cython:
        #  self.sequence_name2sequence_i
        #  self.sequence_i2sequence_name
        #  self.bamfile_data  numpy array
        #     with (seq_i, read_i, pair_i, rlen, pos, is_reverse)
        #  self.cigars

        bamfile = io.AlignmentFile(self.current_bam_filename)
        amplicon.process_bamfile(self, bamfile, BOWTIE_ASCII_OFFSET)

        self.n_sequences = len(self.sequence_name2sequence_i)
        INFO("Number of references with mappings: %s" % self.n_sequences)

        self.priors.append(np.zeros(self.n_sequences, dtype=np.float))
        self.likelihoods = sparse.coo_matrix(
            (self.n_sequences, self.n_reads), dtype=np.float
        )  # init all to zero.
        self.posteriors.append(sparse.lil_matrix(
            (self.n_sequences + 1, self.n_reads + 1), dtype=np.float)
        )

        # TODO: is this necessary any more?
        # or is bookkeeping with probN good enough now.
        self.probN = [None] * self.n_sequences
        # adjusted initialization to be same as probN as was done in emirge.py
        self.prob_indels = [None] * self.n_sequences
        self.unmapped_bases = [None] * self.n_sequences
        # need to calculate this each time? can't we set this once as
        # self.readlengths doesn't change with iters?
        self.mean_read_length = np.mean(self.readlengths)

        # reset probN for valid sequences (from
        # current_reference_fasta_filename). is this still necessary?
        # Or do I keep probN bookkeeping in order already?
        amplicon.reset_probN(self, bamfile)  # also updates coverage values and
        # culls via fraction of length covered, NEW: resets prob_indels as well

        for d in [self.priors, self.posteriors]:
            if len(d) > 2:
                trash = d.pop(0)  # no longer care about t-2
                del trash

    @log.timed("Initializing EM")
    def initialize_EM(self, randomize_priors=False):
        """
        Set up EM with two things so that first iteration can proceed:
           - Initial guesses of Pr(S) are made purely based on read counts,
             where each read is only allowed to map only once to a single best
             reference  (**if more than one alignment reported per read, raise
             exception!**).
           - Initial guess of Pr(N=n) (necessary for likelihood in Pr(S|R) is
             also calculated simply, with the assumption of 1 read (the best
             again) mapped to exactly 1 sequence.  Thus Pr(N=n) only takes the
             base call errors into account.  This is actually not done here,
             but rather the first time self.calc_probN is called.

           - bamfile for iteration 0 is assumed to have just one ("best")
             mapping per read.
           - there is no t-1 for t = 0, hence the need to set up Pr(S)

           if randomize_priors == True, then after calculating priors,
           shuffle them randomly.  This is useful for debugging
           purposes, to test effect of initialization, robustness of
           results, and how often the algorithm gets stuck in local
           maxima.
        """

        self.iteration_i = -1
        self.read_bam()

        # initialize priors.  Here just adding a count for each read mapped to
        # each reference sequence
        # since bowtie run with --best and reporting just 1 alignment at
        # random, there is some stochasticity here.
        for (seq_i, read_i, pair_i, rlen, pos, is_reverse) in self.bamfile_data:
            # if self.probN[seq_i] is not None:
            self.priors[-1][seq_i] += 1

        # this shouldn't be necessary with way I do initial mapping right now
        # (all seq_i in priors should be nonzero initially)
        # only divide cells with at least one count.
        # Set all others to Pr(S) = 0
        nonzero_indices = np.nonzero(self.priors[-1])
        # turn these into probabilities
        self.priors[-1] = self.priors[-1][nonzero_indices] / \
                          self.priors[-1][nonzero_indices].sum()

        if randomize_priors:
            np.random.shuffle(self.priors[-1])

        # push this back to t-1 (index == -2)
        self.priors.append(self.priors[-1].copy())

        # write priors as special case:
        self.print_priors(os.path.join(self.cwd, "priors.initialized.txt"))

    @log.timed("Running iteration {self.iteration_i}+1")  # FIXME i+1
    def do_iteration(self):
        """
        This starts with the M-step, so it requires that Pr(S) and Pr(N=n)
        from previous round are set.
        Pr(S) is used from the previous round's E-step.
        Pr(N=n) partially depends on the previous round's M-step.
        Once M-step is done, then E-step calculates Pr(S) based upon the
        just-calculated M-step.
        """
        self.iteration_i += 1

        self.iterdir = os.path.join(self.cwd, "%s%02d"
                                    % (self.iterdir_prefix, self.iteration_i))
        check_call("mkdir -p " + self.iterdir, shell=True)

        # initialize all data structures.
        self.read_bam()

        # m-step

        # first calculate self.probN from mapped reads, previous
        # round's posteriors
        # (handles initial iteration differently within this method):
        amplicon.calc_probN(self)
        amplicon.calc_likelihood(self)
        amplicon.calc_posteriors(self)

        # now e-step
        self.calc_priors()

        # now write a new fasta file.  Cull sequences below self.min_depth
        consensus_filename = os.path.join(self.iterdir, "iter.%02d.cons.fasta"
                                          % self.iteration_i)
        # culls and splits
        self.write_consensus(consensus_filename,
                             self.current_reference_fasta_filename)
        # if self.iteration_i == self.max_iterations:
        # merge sequences that have evolved to be the same (VSEARCH)
        self.cluster_sequences(consensus_filename)

        # leave a few things around for later.  Note that print_priors also
        #  leaves sequence_name2sequence_i mapping, basically.
        self.print_priors()
        self.print_probN()

        # TODO - compress bam from iteration 0

        # now do a new mapping run for next iteration
        self.do_mapping(consensus_filename)

    @log.timed("Writing priors for iteration {self.iteration_i}")
    def print_priors(self, ofname=None):
        """
        leave a file in directory with nonzero priors printed out.
        """
        if ofname is not None:
            of = file(ofname, 'w')
        else:
            of = file(os.path.join(self.iterdir,
                                   "priors.iter.%02d.txt" % self.iteration_i),
                      'w')
        sequence_i2sequence_name_array = \
            np.array(self.sequence_i2sequence_name)  # faster slicing?
        for seq_i, prior in enumerate(self.priors[-1]):
            seqname = sequence_i2sequence_name_array[seq_i]
            of.write("%d\t%s\t%.10f\n" % (seq_i, seqname, prior))

        of.close()

    @log.timed("Writing probN for iteration {self.iteration_i}")
    def print_probN(self):
        # python gzip.GzipFile is slow.  Use system call to gzip instead
        pickled_filename = os.path.join(self.iterdir, 'probN.pkl')
        cPickle.dump(self.probN, file(pickled_filename, 'w'),
                     cPickle.HIGHEST_PROTOCOL)
        check_call("gzip -f %s" % pickled_filename,
                   shell=True, stdout=sys.stdout, stderr=sys.stderr)

    def calc_priors(self):
        """
        calculates priors [ Pr(S) ] based on
            Pr(S|R) (current posteriors from previous M step, this iteration)
        """
        # here we do have column summing with the posteriors
        # therefore, should be csc sparse type for efficient summing
        self.posteriors[-1] = self.posteriors[-1].tocsc()
        self.priors[-1] = \
            np.asarray(self.posteriors[-1].sum(axis=1)).flatten() / \
            self.posteriors[-1].sum()

    @log.timed("Writing consensus for iteration {self.iteration_i}")
    def write_consensus(self, outputfilename, reference_fasta_filename):
        """
        writes a consensus, taking the most probable base at each position,
        according to current values in Pr(N=n) (self.probN)

        only write sequences with coverage above self.min_depth (culling)
        split sequences with many minor alleles:
             self.snp_minor_prob_thresh
                if prob(N) for minor allele base N is >= this threshold,
                call site a minor allele
             self.snp_percentage_thresh
                if >= this percentage of bases are minor alleles (according to
                self.snp_minor_prob_thresh), then split this sequence into two
                sequences.

        """
        INFO("\tsnp_minor_prob_thresh = %.3f" % self.snp_minor_prob_thresh)
        INFO("\tsnp_percentage_thresh = %.3f" % self.snp_percentage_thresh)

        splitcount = 0
        of = file(outputfilename, 'w')
        tmp_fastafilename = outputfilename + ".tmp.fasta" 
        of_tmp = file(tmp_fastafilename, 'w')

        # DEBUG:
        times_split = []
        times_posteriors = []
        seqs_to_process = len(self.probN)

        i2base = self.i2base
        # these are for updating posteriors at end with new minor strains
        rows_to_add = []
        cols_to_add = []
        data_to_add = []
        probNtoadd = []  # for newly split out sequences

        # just to make sure this is in row-access-friendly format
        self.posteriors[-1] = self.posteriors[-1].tolil()
        reference_fastafile = pysam.Fastafile(reference_fasta_filename)
        self.num_seqs = 0
        loop_t0 = time()
        for seq_i in range(len(self.probN)):
            seq_i_t0 = time()
            if self.probN[seq_i] is None:
                # this sequence is no longer present in this iteration
                # or was culled in reset_probN
                continue
            # FOLLOWING CULLING RULES REMOVED in favor of length-coverage
            # culling in reset_probN()
            # check if coverage passes self.min_depth, if not don't write it
            # (culling happens here)
            # if self.min_depth is not None and self.coverage[seq_i] \
            #    < self.min_depth: #  and self.iteration_i > 5:
            #     # could adjust priors and posteriors here, but because
            #     # prior will already be low (b/c of low coverage) and
            #     # because next round will have 0 mappings (no sequence
            #     # in reference file to map to), this seems
            #     # unnecessary.

            #     # probNarray = None  # NOT PASSED BY REF, assignment is only
            #                            local?
            #     self.probN[seq_i] = None
            #     cullcount += 1
            #     continue # continue == don't write it to consensus.

            # else passes culling thresholds
            title = str(self.sequence_i2sequence_name[seq_i])
            consensus = np.array([i2base.get(x, "N")
                                  for x in np.argsort(
                                       self.probN[seq_i])[:, -1]])
            orig_bases = np.array(reference_fastafile.fetch(title).lower(),
                                  dtype='c')

            # check for deletion, collect deletion sites:
            # need to decide relative weight of competing insertions/deletions
            deletion_threshold = self.indel_thresh / 2
            # retrieve single numpy matrix of prob_indel values for seq_i from
            # prob_indels list of numpy matrices
            prob_indels_single = self.prob_indels[seq_i]
            del_hits = []
            insertion_hits = []
            for base_i in range(prob_indels_single.shape[0]):
                # Eval if insertion or deletion exists at base position
                # Divides weight of deletion, by the sum of both deletions
                # and matches
                denominator = (prob_indels_single[base_i, 0] +
                               prob_indels_single[base_i, 2])
                if denominator < 0:
                    raise ValueError(
                            "denominator should never be < 0 (actual: %s)"
                            % denominator
                    )
                else:
                    if (prob_indels_single[base_i, 2] / denominator) \
                            > deletion_threshold:
                        del_hits.append(base_i)
                    else:
                        self.probN[seq_i][base_i, 4] = 0
                    
                    if (prob_indels_single[base_i, 1] / denominator) \
                        > self.indel_thresh:
                            insertion_hits.append(base_i)
                         
            deletion_indices = np.array(del_hits)
            insertion_indices = np.array(insertion_hits)
        

            # check for minor allele consensus, SPLIT sequence into two
            # candidate sequences if passes thresholds.
            minor_indices = np.argwhere(
                    (
                        self.probN[seq_i] >= self.snp_minor_prob_thresh
                    ).sum(axis=1) >= 2
            )[:, 0]
            if minor_indices.shape[0] > 0:
                minor_fraction_avg = \
                    np.mean(
                            self.probN[seq_i][
                                (
                                    minor_indices,
                                    np.argsort(
                                            self.probN[seq_i][minor_indices]
                                    )[:, -2]
                                )
                            ]
                    )
            else:
                minor_fraction_avg = 0.0
            # NEW rule: only split sequence if *expected* coverage
            # of newly split minor sequence (assuming uniform read
            # coverage over reconstructed sequence) is > some
            # threshold.  Here, expected coverage is calculated
            # based on:
            # Prior(seq_i) * number of MAPPED reads * avg read length
            # * 2 seq per pair
            expected_coverage_minor = \
                (self.priors[-1][seq_i] * minor_fraction_avg *
                 self.n_reads_mapped * self.mean_read_length) / \
                self.probN[seq_i].shape[0]
            expected_coverage_major = \
                (self.priors[-1][seq_i] * (1-minor_fraction_avg) *
                 self.n_reads_mapped * self.mean_read_length) / \
                self.probN[seq_i].shape[0]

            # multiply by 2 because n_reads_mapped is actually number of
            # mapped pairs
            if self.reads2_filepath is not None:
                expected_coverage_minor *= 2.0
                expected_coverage_major *= 2.0


            #if deletion_indices.shape[0] or insertion_indices.shape[0] > 0 or (
            #    (deletion_indices.shape[0] + insertion_indices.shape[0] + \
            #    minor_indices.shape[0]) /
            #    float(self.probN[seq_i].shape[0]) >=
            #    self.snp_percentage_thresh and
            #    expected_coverage_minor >= self.expected_coverage_thresh
            #):
            if minor_indices.shape[0]+ insertion_indices.shape[0] /float(self.probN[seq_i].shape[0]) >= self.snp_percentage_thresh and expected_coverage_minor >= self.expected_coverage_thresh:
                # We split!
                splitcount += 1
                # if there's >=3 alleles, major allele keeps prob of other
                # minors)
                major_fraction_avg = 1.-minor_fraction_avg
                # -2 gets second most probable base:
                minor_bases = np.array(
                        [i2base.get(x, "N")
                         for x in np.argsort(
                                self.probN[seq_i][minor_indices]
                         )[:, -2]]  # will this work?
                )
                # get a copy of the consensus
                minor_consensus = consensus.copy()
                # replace the bases that pass minor threshold
                minor_consensus[minor_indices] = minor_bases

                # remove deletion bases from minor_consensus - deletions form
                # part of this new sequence
                # new_minor_c=[]
                # for i in range(minor_consensus.shape[0]):
                #    if not i in deletion_indices:
                #        new_minor_c.append(minor_consensus[i])
                # new_minor_c = np.array(new_minor_c)

                # now deal with naming.
                title_root = re.search(r'(.+)(_m(\d+))$', title)
                if title_root is None:  # no _m00 on this name
                    title_root = title[:]
                else:
                    title_root = title_root.groups()[0]
                # now check for any known name with same root and a _m on it.
                previous_m_max = max([0] + [
                    int(x) for x in re.findall(
                            r'%s_m(\d+)' % re.escape(title_root),
                            " ".join(self.sequence_i2sequence_name))
                    ])
                m_title = "%s_m%02d" % (title_root, previous_m_max+1)

                # also split out Priors and Posteriors (which will be used in
                # next round), split with average ratio of major to minor
                # alleles.
                # updating priors first:
                old_prior = self.priors[-1][seq_i]
                self.priors[-1][seq_i] = old_prior * major_fraction_avg
                seq_i_minor = self.n_sequences
                self.n_sequences += 1
                self.sequence_i2sequence_name.append(m_title)
                assert len(self.sequence_i2sequence_name) == self.n_sequences
                assert len(self.sequence_i2sequence_name) == seq_i_minor + 1
                self.sequence_name2sequence_i[m_title] = seq_i_minor
                self.split_seq_first_appeared[seq_i] = self.iteration_i
                self.refseq_lengths.append(consensus.shape[0])
                assert len(self.refseq_lengths) == self.n_sequences
                # how I adjust probN here for newly split seq doesn't really
                # matter, as it is re-calculated next iter.
                # this only matters for probN.pkl.gz file left behind for this
                # iteration.
                # for now just set prob(major base) = 0 and redistribute prob
                # to other bases for minor, and set prob(minor base) = 0 and
                # redistribute prob to other bases for major
                # MINOR
                major_base_i = np.argsort(
                        self.probN[seq_i][minor_indices]
                )[:, -1]
                newprobNarray = self.probN[seq_i].copy()
                newprobNarray[(minor_indices, major_base_i)] = 0
                newprobNarray = newprobNarray / \
                                np.sum(newprobNarray, axis=1).reshape(
                                        newprobNarray.shape[0], 1)
                probNtoadd.append(newprobNarray)

                # why are we doing this?!?!?
                self.base_coverages.append(np.zeros_like(
                        self.base_coverages[seq_i]))

                # MAJOR
                minor_base_i = np.argsort(
                        self.probN[seq_i][minor_indices]
                )[:, -2]
                self.probN[seq_i][(minor_indices, minor_base_i)] = 0
                self.probN[seq_i] = self.probN[seq_i] / \
                                    np.sum(self.probN[seq_i], axis=1
                                           ).reshape(self.probN[seq_i].shape[0],
                                                     1)

                new_priors = np.zeros(seq_i_minor + 1,
                                      dtype=self.priors[-1].dtype)
                new_priors[:-1] = self.priors[-1].copy()
                new_priors[seq_i_minor] = old_prior * minor_fraction_avg
                trash = self.priors.pop()
                del trash
                self.priors.append(new_priors)

                # keep track of all new minor data to add and add it
                # once at end for ALL split sequences with one coo
                # matrix construction, instead of each iteration.

                t_posterior = time()
                # new_read_probs, new_rows, new_cols =
                #   adjust_posteriors_for_split(AAAA, BBBB, CCCC) # TODO:
                # could move to Cython
                # updating posteriors. for each seq-read pair with prob > 0,
                # split prob out to major and minor seq.
                new_cols = self.posteriors[-1].rows[seq_i]  # col in coo format
                # data in coo format:
                new_read_probs = [x * minor_fraction_avg
                                   for x in self.posteriors[-1].data[seq_i]]
                new_rows = [seq_i_minor for _ in new_cols]  # row in coo

                # add new read probs to cache of new read probs to add at
                # end of loop
                rows_to_add.extend(new_rows)
                cols_to_add.extend(new_cols)
                data_to_add.extend(new_read_probs)

                # adjust old read probs to reflect major strain
                self.posteriors[-1].data[seq_i] = \
                    [x * major_fraction_avg
                     for x in self.posteriors[-1].data[seq_i]]
                times_posteriors.append(time() - t_posterior)

                # adjust self.unmapped_bases (used in clustering).
                # For now give same pattern as parent
                self.unmapped_bases.append(self.unmapped_bases[seq_i].copy())

                # write out minor strain consensus
                of.write(">%s\n" % m_title)
                of_tmp.write(">%s\n" % m_title)
                # added for indels:
                minor_consensus = self.eval_indels(seq_i, minor_consensus,
                                                   m_title, orig_bases,
                                                   mask="soft", major="minor")
                of.write("%s\n" % "".join(minor_consensus))
                of_tmp.write("%s\n" % "".join(minor_consensus))
                self.unmapped_bases.append(np.zeros(len(consensus),
                                                    dtype=np.uint8))
                log.info("new sequence %s is length %s" % (title, len(consensus)))
                #print "new sequence %s is length %s" % (title, len(consensus))
                self.num_seqs += 1
                log.info("splitting sequence %d (%s) to %d (%s)...\n"
                         % (seq_i, title, seq_i_minor, m_title))
                times_split.append(time() - seq_i_t0)

            # now write major strain consensus, regardless of whether there
            # was a minor strain consensus
            of.write(">%s\n" % title)
            of_tmp.write(">%s\n" % title)
            consensus = self.eval_indels(seq_i, consensus, title, orig_bases,
                                         mask="soft",major="major")  # added for indels
            of.write("%s\n" % "".join(consensus))
            of_tmp.write("%s\n" % "".join(consensus))
            self.num_seqs += 1

        # END LOOP
        loop_t_total = time() - loop_t0
        # update posteriors matrix with newly added minor sequences
        # new_posteriors via coo, then convert to csr.
        # first make a copy in coo format:
        new_posteriors = self.posteriors[-1].tocoo()
        # then create new coo matrix with new shape, appending new row, col,
        # data to old row, col, data

        new_posteriors = sparse.coo_matrix(
                (np.concatenate((new_posteriors.data, data_to_add)),
                 (np.concatenate((new_posteriors.row, rows_to_add)),
                  np.concatenate((new_posteriors.col, cols_to_add)))),
                shape=(self.n_sequences, self.posteriors[-1].shape[1]),
                dtype=new_posteriors.dtype
        ).tocsr()

        # finally, exchange in this new matrix
        trash = self.posteriors.pop()
        del trash
        self.posteriors.append(new_posteriors)

        # update probN array:
        self.probN.extend(probNtoadd)

        log.info("Split out %d new minor strain sequences." % splitcount)
        if splitcount > 0:
            log.info("Average time for split sequence: [%.6f seconds]"
                     % np.mean(times_split))
            log.info("Average time for posterior update: [%.6f seconds]"
                     % np.mean(times_posteriors))
            log.info("Average time for non-split sequences: [%.6f seconds]"
                     % ((loop_t_total - sum(times_split)) /
                        (seqs_to_process - len(times_split))))

    def eval_indels(self, seq_i, consensus, title, orig_bases, mask, major):
        # Evaluates consensus sequence for write outs against the prob_indels
        # array.  deletes or inserts bases as appropriate
        # OUT:   returns a list of single character bases as new consensus
        insertion_threshold = self.indel_thresh
        new_cons = []
        # retrieve single numpy matrix of prob_indel values for seq_i from
        # prob_indels list of numpy matrices
        prob_indels_single = self.prob_indels[seq_i]
        unmapped_bases_single = self.unmapped_bases[seq_i]
    
        for base_i in range(prob_indels_single.shape[0]):
            # Eval if deletion exists at base position
            # Divides weight of deletion, by the sum of both deletions and
            # matches
            del_count = 0
            denominator = (prob_indels_single[base_i, 0] +
                           prob_indels_single[base_i, 2])
            
            if denominator < 0:
                raise ValueError(
                        "denominator should never be < 0 (actual: %s)"
                        % denominator
                )
            elif denominator == 0:
                # nothing mapped to this base.
                # Use consensus base from reference
                new_cons.append(consensus[base_i])
            else:
                if unmapped_bases_single[base_i] == 1:
                    if mask == "hard":
                        new_cons.append("N")
                    elif mask == "soft":
                        # return to original base if unmapped / ambiguous
                        new_cons.append(orig_bases[base_i])
                elif consensus[base_i] != "D":
                    # not deleted
                    new_cons.append(consensus[base_i])
                    del_count += 1
                else:
                    # delete (add nothing to new consensus)
                    INFO("Modified reference sequence %d (%s) with a deletion "
                         "of base %d " % (seq_i, title, base_i))
                    continue
                    
                # Keep insertions in seqs without deletions - DEBUG
                #if del_count > 0:
                 #   continue
                # Eval if insertion exists after base position
                if (prob_indels_single[base_i, 1]) == 0:
                    # no evidence for insertion
                    continue
                # if summed weights of insertion is greater than sum of reads
                # mapped nearby (at base on left flank of proposed insertion
                #  (because it's easy), i-1)
                
                elif (prob_indels_single[base_i, 1] /denominator) < insertion_threshold:
                    continue # no insertion
                
                elif ((prob_indels_single[base_i, 1] / denominator) > 2*insertion_threshold) and major=="major":
                        for i in range(int(prob_indels_single[base_i, 3])):
                            new_cons.append('N')
                        log.info("Modified reference sequence %d (%s) with an "
                             "insertion of %s bases after base %d"
                             % (seq_i, title, str(int(prob_indels_single[
                                                          base_i, 3])), base_i))
                else:
                    if ((prob_indels_single[base_i, 1] / denominator) > insertion_threshold) and major=="minor":
                        for i in range(int(prob_indels_single[base_i, 3])):
                            new_cons.append('N')
                    log.info("Modified reference sequence %d (%s) with an "
                             "insertion of %s bases after base %d"
                             % (seq_i, title, str(int(prob_indels_single[
                                                          base_i, 3])), base_i))

        return new_cons

    def cluster_sequences(self, fastafilename):
        """
        Right now, this simply calls cluster_sequences_vsearch, which
        uses VSEARCH.  Could swap in other functions here if there
        were faster or just alternative clustering methods to try out

        called function should also adjust Pr(S) [prior] and Pr(S_t-1)
        [posteriors] as needed after merging.
        """
        return self.cluster_sequences_vsearch(fastafilename)

    @log.timed("Clustering sequences for iteration {self.iteration_i}")
    def cluster_sequences_vsearch(self, fastafilename):
        """
        uses VSEARCH to merge sequences above self.cluster_thresh %ID over the
        length of the shorter sequence

        Merge two sequences if the *NON-GAPPED* positions have %
        identity >= self.cluster_thresh

        also adjusts Pr(S) [prior] and Pr(S_t-1) [posteriors] as needed after
        merging.
        """
        # Hard threshold in place for EMIRGE2.  cluster_thresh from -j flag
        # is now applied only in post-processing steps
        log.info("cluster threshold = 100.00%")

        tocleanup = []  # list of temporary files to remove after done

        # get posteriors ready for slicing (just prior to this call,
        # is csr matrix?):
        self.posteriors[-1] = self.posteriors[-1].tolil()

        # NOTE that this fasta file now contains N's where there are
        # no mapped bases, so that vsearch with iddef 0 will not count
        # positions aligned to these bases in the identity calculation

        tmp_fastafilename = fastafilename + ".tmp.fasta"
        tocleanup.append(tmp_fastafilename)
        tmp_fastafile = pysam.Fastafile(tmp_fastafilename)
        tocleanup.append("%s.fai" % tmp_fastafilename)
        # do global alignments with VSEARCH
        # I don't use --cluster because it doesn't report alignments
        # vsearch is fast but will sometimes miss things -- I've tried to tune
        # params as best as I can.
        # and I use different parameters depending on how many input sequences
        #  there are
        # Also, I use a lower %ID thresh than specified for joining because I
        # really calculate %ID over *mapped* sequence positions.

        # sens_string = "--maxaccepts 8 --maxrejects 256"
        # if self.iteration_i == self.max_iterations:
        #        uclust_id=0.80
        # else:
        #        uclust_id=1.0
        sens_string = "--maxaccepts 8 --maxrejects 256"
        uclust_id = 0.80
        algorithm = "-usearch_global"

        # if em.iteration_i > 10:
        # num_seqs = len([x for x in self.probN if x is not None])
        assert self.num_seqs == len([x for x in self.probN if x is not None])
        if self.num_seqs < 1000:
            sens_string = "--maxaccepts 16 --maxrejects 256"
        if self.num_seqs < 500:
            sens_string = "--maxaccepts 32 --maxrejects 256"
        if self.num_seqs < 150:
            algorithm = "-usearch_global"
            # slower, but more sensitive.
            sens_string = "--maxaccepts 0 --maxrejects 0"
        # if really few seqs, then no use not doing smith-waterman or needleman
        # wunsch alignments
        if self.num_seqs < 50:
            algorithm = "-usearch_global"
            sens_string = "-fulldp"

        cmd = "vsearch %s %s --db %s --id %.3f --query_cov 0.5 " \
              "--target_cov 0.5 --strand plus --userout %s.us.txt " \
              "--userfields query+target+id+caln+qlo+qhi+tlo+thi " \
              "--threads %d %s --quiet" % \
              (algorithm,
               tmp_fastafilename, tmp_fastafilename,
               uclust_id,
               tmp_fastafilename,
               self.n_cpus,
               sens_string)

        log.info("vsearch command was:\n%s" % cmd)

        check_call(cmd, shell=True, stdout=sys.stdout, stderr=sys.stderr)
        # read clustering file to adjust Priors and Posteriors, summing merged
        # reference sequences
        tocleanup.append("%s.us.txt" % tmp_fastafilename)

        nummerged = 0
        alnstring_pat = re.compile(r'(\d*)([MDI])')
        already_removed = set()  # seq_ids
        # this is a bit slow and almost certainly could be sped up with
        # algorithmic improvements.
        times = []  # DEBUG
        for row in csv.reader(file("%s.us.txt" % tmp_fastafilename),
                              delimiter='\t'):
            # each row an alignment in userout file
            # member == query
            member_name = row[0]
            seed_name = row[1]
            if member_name == seed_name:
                continue  # vsearch allows self-hits, which we don't care about
            member_seq_id = self.sequence_name2sequence_i.get(member_name)
            seed_seq_id = self.sequence_name2sequence_i.get(seed_name)
            if member_seq_id in already_removed or \
                            seed_seq_id in already_removed:
                continue
            # if sequences are different lengths, can't be 100% ID right?
            # So continue
            if self.unmapped_bases[member_seq_id].shape[0] != \
                    self.unmapped_bases[seed_seq_id].shape[0]:
                continue
            if len(alnstring_pat.findall(row[3])) > 1:
                continue
                # problem dealing with array lengths in _emirge.count_cigar_aln,
                # so bypassing for now not sure if this is right approach, but
                # should be excluding non-identical sequences (where
                # cigarstring indicates D or I)
                # decide if these pass the cluster_thresh *over non-gapped,
                # mapped columns*
                # query+target+id+caln+qlo+qhi+tlo+thi %s"%\
                #   0     1     2   3   4   5  6    7
            #member_start = int(row[4]) - 1  # printed as 1-based by vsearch now
            #seed_start = int(row[6]) - 1
            percent_id = float(row[2])
            aln_columns = alnstring_pat.findall(row[3])[0][0]

            #aln_columns, matches = amplicon.count_cigar_aln(
            #        tmp_fastafile.fetch(seed_name),
            #        tmp_fastafile.fetch(member_name),
            #        self.unmapped_bases[seed_seq_id],
            #        self.unmapped_bases[member_seq_id],
            #        seed_start,
            #        member_start,
            #        alnstring_pat.findall(row[3])
            #)

            # if alignment is less than 1000 bases, or identity over those
            # 500+ bases is not above thresh, then continue

            #if aln_columns < 500 \
             #  or ((float(matches) / aln_columns) < 1.0):
                # Fixed merging threshold for EMIRGE2.  During iterations only
                #  merge seqs that are 100% identical. Rest are cleaned up in
                # post-processing steps.
                # or ((float(matches) / aln_columns) < self.cluster_thresh):
                # or (float(aln_columns) / min(seed_n_mapped_bases,
                # member_n_mapped_bases) < 0.9)
              #  continue

            # how many iters does a newly split out seq have to be around
            # before it's allowed to merge again.  -1 to turn this off.
            minimum_residence_time = -1
            member_first_appeared = \
                self.split_seq_first_appeared.get(member_seq_id)
            if member_first_appeared is not None and \
               self.iteration_i - member_first_appeared <= \
               minimum_residence_time:
                continue
            seed_first_appeared = self.split_seq_first_appeared.get(seed_seq_id)
            if seed_first_appeared is not None and \
               self.iteration_i - seed_first_appeared <= minimum_residence_time:
                continue

            #if self.num_seqs < 50:
            #    log.info("\t\t%s|%s vs %s|%s %.3f over %s aligned columns"
            #             "(vsearch %%ID: %s)"
            #             % (member_seq_id, member_name, seed_seq_id, seed_name,
            #                float(matches) / aln_columns, aln_columns, row[2]))

            # if above thresh, then first decide which sequence to keep,
            # (one with higher prior probability).
            #percent_id = (float(matches) / aln_columns) * 100. 
            t0 = time()
            if self.priors[-1][seed_seq_id] > self.priors[-1][member_seq_id]:
                keep_seq_id = seed_seq_id
                remove_seq_id = member_seq_id
                keep_name = seed_name
                remove_name = member_name
            else:
                keep_seq_id = member_seq_id
                remove_seq_id = seed_seq_id
                keep_name = member_name
                remove_name = seed_name

            # merge priors (add remove_seq_id probs to keep_seq_id probs).
            self.priors[-1][keep_seq_id] += self.priors[-1][remove_seq_id]
            self.priors[-1][remove_seq_id] = 0.0

            # now merge posteriors (all removed probs from remove_seq_id go to
            # keep_seq_id).
            # self.posteriors[-1] at this point is lil_matrix
            # some manipulations of underlying sparse matrix data structures
            # for efficiency here.
            # 1st, do addition in csr format (fast), convert to lil format, and
            # store result in temporary array.
            new_row = (self.posteriors[-1].getrow(keep_seq_id).tocsr() +
                       self.posteriors[-1].getrow(remove_seq_id).tocsr()
                       ).tolil()
            # then change linked lists directly in the posteriors data
            # structure -- this is very fast
            self.posteriors[-1].data[keep_seq_id] = new_row.data[0]
            self.posteriors[-1].rows[keep_seq_id] = new_row.rows[0]
            # these two lines remove the row from the linked list (or rather,
            # make them empty rows), essentially setting all elements to 0
            self.posteriors[-1].rows[remove_seq_id] = []
            self.posteriors[-1].data[remove_seq_id] = []

            # set self.probN[removed] to be None -- note that this doesn't
            # really matter, except for
            # writing out probN.pkl.gz every iteration, as probN is recalculated
            # from bam file
            # with each iteration
            self.probN[remove_seq_id] = None

            already_removed.add(remove_seq_id)
            nummerged += 1

            times.append(time() - t0)
            log.info("\t...merging %d|%s into %d|%s (%.2f%% ID over %s columns)"
                     "in %.3f seconds"
                     % (remove_seq_id, remove_name, keep_seq_id, keep_name,
                        percent_id, aln_columns, times[-1]))

        # write new fasta file with only new sequences
        log.info("Writing new fasta file for iteration %d" % self.iteration_i)
        tmp_fastafile.close()
        recordstrings = ""
        num_seqs = 0
        # read through file again, overwriting orig file if we keep the seq
        for record in io.FastIterator(file(fastafilename)):
            seqname = record.title.split()[0]
            seq_id = self.sequence_name2sequence_i.get(seqname)
            if seq_id not in already_removed:
                # could do a better job here of actually "merging" a new
                # consensus, rather than just keeping one or the other.
                recordstrings += str(record)
                num_seqs += 1
        outfile = file(fastafilename, 'w')
        outfile.write(recordstrings)
        outfile.close()

        # clean up.  quite important, actually, to remove old fai index files.
        for fn in tocleanup:
            os.remove(fn)
            
        log.info("\tremoved %d sequences after merging" % nummerged)
        log.info("\tsequences remaining for iteration %02d: %d"
                 % (self.iteration_i, num_seqs))

    def get_n_alignments_from_bowtie(self):
        """
        bowtie2 output has changed,
        use this to get number of aligned reads
        """
        get_n = make_pipe("get_n_alignments_from_sam",
                          ["samtools", "view", "-c", self.current_bam_filename])
        with get_n() as samtools:
            n = samtools.next()
        return int(n)


    def calc_cutoff_threshold(self):  # done at the end of the final iteration
        """
        calculate the minimum abundance that will represent average coverage
        specified in cov_thresh (default will be 20X) this value is used in
        post-processing steps to filter final EMIRGE fasta file
        """

        with io.AlignmentFile(self.current_bam_filename, "rb") as bamfile:
            self.avg_emirge_seq_length = np.mean(bamfile.lengths)

        log.warning("Average read length is %.4d"
                    % self.mean_read_length)
        log.warning("Average EMIRGE sequence length is %.5d"
                    % self.avg_emirge_seq_length)
        log.warning("Fragments mapped = %.9d"
                    % self.fragments_mapped)
        self.prob_min = (self.avg_emirge_seq_length * float(self.cov_thresh)) \
            / (self.fragments_mapped * ((int(self.paired_end)+1) *
                                        self.mean_read_length))


def do_iterations(em, max_iter):
    """
    an EM object is passed in, so that one could in theory start from a saved
    state this should be moved into the EM object.
    """
    os.chdir(em.cwd)

    em.max_iterations = max_iter

    while em.iteration_i < max_iter:
        em.do_iteration()

    em.calc_cutoff_threshold()

    # clean up any global temporary files, i.e. rewritten reads files
    for filename in em.temporary_files:
        DEBUG("unlinking '{}'".format(filename))
        os.remove(filename)


def post_process(em, working_dir):
    """Do all the postprocessing for the EMIRGE run, producing a fasta file of
    all EMIRGE sequences produced that meet the estimated abundance threshold
    to achieve an estimated average coverage of 20X (default) or other user
    specified value.  Also produces a raw file of all EMIRGE sequences produced,
    but not recommended to be used in later analyses.  Final clustering is
    performed at 97% ID (default) but can be changed with the '-j' flag.
    """

    # first need to do rename(em) with no probmin - keeping all sequences for
    # clustering and writing out the raw output file
    nomin_fasta_filename = os.path.join(working_dir,
                                        em.output_files_prefix +
                                        "_nomin_iter."'%02d'
                                        % em.max_iterations + ".RAW.fasta")
    if os.path.exists(nomin_fasta_filename):
        log.warning("WARNING: overwriting file %s" % nomin_fasta_filename)

    rename(em.iterdir, nomin_fasta_filename,  em.output_files_prefix,
           prob_min=None, no_N=False, no_trim_N=True)

    # next need to cluster using vsearch at the user specified threshold,
    # default =0.97
    centroids = os.path.join(working_dir, "centroids.tmp")
    uc = os.path.join(working_dir, "uc.tmp")
    cmd = "vsearch --cluster_smallmem %s -usersort -notrunclabels -id %.2f" \
          " -centroids %s -uc %s " % (nomin_fasta_filename, em.cluster_thresh,
                                      centroids, uc)
    log.info("vsearch command was:\n%s" % cmd)

    check_call(cmd, shell=True, stdout=sys.stdout, stderr=sys.stderr)

    # now merging cluster abundances and writing out the fasta file with merged
    # abundance over cutoff
    log.warning("Minimum abundance for estimated %.2dX coverage is %.6f"
                % (em.cov_thresh, em.prob_min))
    clustered_fasta = os.path.join(working_dir, em.output_files_prefix +
                                   "_iter."'%02d' % em.max_iterations +
                                   "_"'%.6f' % em.prob_min + ".fasta")
    outfile = open(clustered_fasta, 'w')
    size_pattern = re.compile(r'NormPrior=([0-9\.]+)')

    cluster_sizes = {}
    for line in open(uc, 'r'):
        atoms = line.strip().split('\t')
        if atoms[0] == 'S':  # seed
            # get seed size, add to current cluster size
            cluster_sizes[atoms[8]] = \
                cluster_sizes.get(atoms[8], 0) + \
                float(size_pattern.search(atoms[8]).groups()[0])
        elif atoms[0] == 'H':  # hit
            # get hit size, add to current cluster size
            cluster_sizes[atoms[9]] = \
                cluster_sizes.get(atoms[9], 0) + \
                float(size_pattern.search(atoms[8]).groups()[0])
        elif atoms[0] == 'C':
            break
        else:
            raise ValueError("Unknown code in uc file: %s" % atoms[0])

    good_ids = {}
    for line in open(centroids, 'r'):
        if line.startswith('>'):
            id = line.strip().split(">")[1]
            if cluster_sizes[id] >= em.prob_min:
                # collect IDs of those sequences whose combined cluster
                # abundance is over the calculated$
                good_ids[id] = ((size_pattern.sub(
                        'comb_abund=%f' % cluster_sizes[id], line
                )))

    # header should be just containing comb_abund.  As it stands, haven't
    # replaced/removed Prior=
    for record in SeqIO.parse(open(centroids, 'r'), 'fasta'):
        if record.description in good_ids:
            atoms = good_ids[record.description].strip(">").split(" ")
            record.id = atoms[0]
            record.description = " ".join(atoms[1:])
            SeqIO.write(record, outfile, 'fasta')

    outfile.close()
    os.remove(centroids)
    os.remove(uc)


def main(argv=sys.argv[1:]):
    """
    command line interface to emirge

    """
    parser = OptionParser(USAGE)

    # REQUIRED
    group_reqd = OptionGroup(
            parser, "Required flags",
            "These flags are all required to run EMIRGE, and may be supplied "
            "in any order.")
    group_reqd.add_option(
            "-1", dest="fastq_reads_1",
            metavar="reads_1.fastq[.gz]",
            type="string",
            help="path to fastq file with \\1 (forward) reads from "
                 "paired-end sequencing run, or all reads from single-end "
                 "sequencing run.  File may optionally be gzipped.  EMIRGE "
                 "expects ASCII-offset of 64 for quality scores (but see "
                 "--phred33).  (Note that running EMIRGE with single-end "
                 "reads is largely untested.  Please let me know how it "
                 "works for you.)")
    group_reqd.add_option(
            "-f", "--fasta_db",
            type="string",
            help="path to fasta file of candidate SSU sequences")
    group_reqd.add_option(
            "-o", "--output_files_prefix",
            type="string", default="EMIRGE",
            help="prefix to be used for final output fasta file and EMIRGE "
                 "sequences reconstructed")
    parser.add_option_group(group_reqd)

    # REQUIRED for paired end
    group_reqd_pe = OptionGroup(
            parser, "Required flags for paired-end reads",
            "These flags are required to run EMIRGE when you have paired-end "
            "reads (the standard way of running EMIRGE), and may be supplied in"
            " any order.")
    group_reqd_pe.add_option(
            "-2", dest="fastq_reads_2", metavar="reads_2.fastq",
            type="string",
            help="path to fastq file with \\2 (reverse) reads from paired-end "
                 "run.  File must be unzipped for mapper.  EMIRGE expects "
                 "ASCII-offset of 64 for quality scores (but see --phred33).")
    parser.add_option_group(group_reqd_pe)

    # OPTIONAL
    group_opt = OptionGroup(
            parser, "Optional parameters",
            "Defaults should normally be fine for these options in order to "
            "run EMIRGE")
    group_opt.add_option(
            "-n", "--iterations",
            type="int", default=40,
            help="Number of iterations to perform.  It may be necessary to "
                 "use more iterations for more complex samples "
                 "(default=%default)")
    group_opt.add_option(
            "-a", "--processors",
            type="int", default=multiprocessing.cpu_count(),
            help="""Number of processors to use in the mapping steps.  You
            probably want to raise this if you have the processors. (default:
            use all available processors)""")
    group_opt.add_option(
            "-p", "--snp_fraction_thresh",
            type="float", default="0.04",
            help="If fraction of variants in a candidate sequence "
                 "exceeds this threshold, then split the candidate "
                 "into two sequences for next iteration.  See also "
                 "--variant_fraction_thresh. (default: %default)")
    group_opt.add_option(
            "-v", "--variant_fraction_thresh",
            type="float", default="0.1",
            help="minimum probability of second most probable base at a site "
                 "required in order to call site a variant.  See also "
                 "--snp_fraction_thresh.  (default: %default)")
    group_opt.add_option(
            "-j", "--join_threshold",
            type="float", default="0.97",
            help="If two candidate sequences share >= this fractional identity "
                 "over their bases with mapped reads, then merge the two "
                 "sequences into one for the next iteration.  (default: "
                 "%default; valid range: [0.95, 1.0] ) ")
    group_opt.add_option(
            "-d", "--debug",
            action="store_true", default=False,
            help="print debug information")
    group_opt.add_option(
            "-q", "--quiet",
            action="store_true", default=False,
            help="be less verbose")
    group_opt.add_option(
            "--indel_thresh",
            type="float", default=0.3,
            help="temporary flag to test indel thresholds")
    group_opt.add_option(
            "-c", "--min_length_coverage",
            type="float",
            default=0.3,
            help="minimum fraction of the length of a candidate reference "
                 "sequence that must be covered by mapped reads.  If not "
                 "met, a candidate sequence is discarded for the next "
                 "iteration.  (default: %default; valid range: (0.0, 1.0])")
    group_opt.add_option(
            '-t', '--min_coverage_threshold',
            type='int', default='20',
            help="Expected minimum depth coverage.  Sequences are only "
                 "included in output FASTA files if their expected coverage is "
                 "greater than this value (calculated based on "
                 "EMIRGE-estimated abundance and total number of reads "
                 "mapped).  Default: %default")

    # --- HIDDEN --- for debugging or special use case

    # this option randomizes the priors calculated for algorithm
    # initialization.  Useful for testing how init affects final
    # sequences and how likely we are to get stuck in a local maxima.
    group_opt.add_option(
            "--randomize_init_priors",
            action="store_true", default=False,
            help=SUPPRESS_HELP)
    # --- END HIDDEN ---

    parser.add_option_group(group_opt)

    mapping.get_options(parser)

    # ACTUALLY PARSE ARGS
    (options, args) = parser.parse_args(argv)

    if not vsearch.available():
        parser.error("Please install vsearch")

    # configure log level:
    log.setup(quiet=options.quiet, debug=options.debug)

    # minimal sanity checking of input
    if len(args) != 1:
        parser.error("DIR is required, and all options except DIR should have "
                     "a flag associated with them (options without flags: %s)"
                     % args)
    if options.join_threshold < 0.95 or options.join_threshold > 1:
        parser.error("join_threshold must be between [0.95, 1.0].  You "
                     "supplied %.3f. (see --help)" % options.join_threshold)
    if options.min_length_coverage is not None:
        if options.min_length_coverage <= 0 or options.min_length_coverage >= 1:
            parser.error("--min_length_coverage (-c) must be between "
                         "(0.0, 1.0).  You supplied %.3f. (see --help)"
                         % options.min_length_coverage)

    for filename_opt_string in ["fastq_reads_1", "fastq_reads_2", "fasta_db"]:
        filename_option = getattr(options, filename_opt_string)
        if filename_option is not None:
            if not os.path.exists(filename_option):
                parser.error("file not found for --%s: %s"
                             % (filename_opt_string, filename_option))

    working_dir = os.path.abspath(args[0])

    sys.stdout.write("""\
If you use EMIRGE in your work, please cite these manuscripts, as appropriate.

Miller CS, Baker BJ, Thomas BC, Singer SW, Banfield JF (2011)
EMIRGE: reconstruction of full-length ribosomal genes from microbial community
        short read sequencing data.
Genome biology 12: R44. doi:10.1186/gb-2011-12-5-r44.

Miller CS, Handley KM, Wrighton KC, Frischkorn KR, Thomas BC, Banfield JF (2013)
Short-Read Assembly of Full-Length 16S Amplicons Reveals Bacterial Diversity in
Subsurface Sediments.
PloS one 8: e56018. doi:10.1371/journal.pone.0056018.\n\n""")

    sys.stdout.write("imported _emirge C functions from: %s\n"
                     % amplicon.__file__)
    sys.stdout.write("Command:\n")
    sys.stdout.write(' '.join([__file__] + argv))
    sys.stdout.write('\n\n')
    total_start_time = time()
    sys.stdout.write("EMIRGE started at %s\n" % (ctime()))
    sys.stdout.flush()

    for o in ["fastq_reads_1", "fasta_db"]:
        if getattr(options, o) is None or getattr(options, o) == 0:
            if o == 'fasta_db':
                parser.error("Fasta file for candidate database is "
                             "missing. Specify --fasta_db. (try --help for "
                             "more information)")
            else:
                parser.error("--%s is required, but is not specified "
                             "(try --help)" % o)

    if not os.path.exists(working_dir):
        os.mkdir(working_dir)
    elif len(os.listdir(working_dir)) > 1:
        # allow 1 file in case log file is redirected here.
        print >> sys.stderr, os.listdir(working_dir)
        parser.error(
            "Directory not empty: %s\n"
            "It is recommended you run emirge in a new directory each run; "
            "delete this directory or specify a new one."
            % working_dir)

    # clean up options to be absolute paths
    for o in ["fastq_reads_1", "fastq_reads_2", "fasta_db"]:
        current_o_value = getattr(options, o)
        if current_o_value is not None:
            setattr(options, o, os.path.abspath(current_o_value))

    # create mapper object
    mapper = mapping.get_mapper(opts=options,
                                workdir=working_dir,
                                candidates=options.fasta_db,
                                fwd_reads=options.fastq_reads_1,
                                rev_reads=options.fastq_reads_2,
                                threads=options.processors)

    mapper.prepare_reads()

    # finally, CREATE EM OBJECT
    em = EM(mapper=mapper,
            candidate_db=options.fasta_db,
            reads1_filepath=mapper.fwd_reads,
            reads2_filepath=mapper.rev_reads,
            insert_mean=options.insert_mean,
            insert_sd=options.insert_stddev,
            output_files_prefix=options.output_files_prefix,
            cov_thresh=options.min_coverage_threshold,
            max_read_length=mapper.max_rlen,
            cluster_thresh=options.join_threshold,
            indel_thresh=options.indel_thresh,
            n_cpus=options.processors,
            cwd=working_dir,
            reads_ascii_offset={False: 64, True: 33}[mapper.phred33],
            snp_percentage_thresh=options.snp_fraction_thresh,
            snp_minor_prob_thresh=options.variant_fraction_thresh,
            min_length_coverage=options.min_length_coverage)

    em.do_initial_mapping()

    if options.randomize_init_priors:
        print >> sys.stderr, "*" * 60
        print >> sys.stderr, "DEBUG: initialized priors will be randomized " \
                             "for testing purposes"

    em.initialize_EM(randomize_priors=options.randomize_init_priors)

    # BEGIN ITERATIONS
    do_iterations(em, max_iter=options.iterations)

    # WRITE THE OUTPUT FASTA FILES- renames and writes raw file and file with
    # only seqs having estimated coverage over specified threshold (default=20X)
    post_process(em, working_dir)

    # print some info to user about final files produced,
    # brief description, filename
    # TODO
    
    sys.stdout.write("EMIRGE finished at %s.  Total time: %s\n"
                     % (ctime(), timedelta(seconds=time() - total_start_time)))


if __name__ == '__main__':
    main()
