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

EMIRGE reconstructs full length ribosomal genes from short read
sequencing data.  In the process, it also provides estimates of the
sequences' abundances.

EMIRGE uses a modification of the EM algorithm to iterate between
estimating the expected value of the abundance of all SSU sequences
present in a sample and estimating the probabilities for each read
that a specific sequence generated that read.  At the end of each
iteration, those probabilities are used to re-calculate (correct) a
consensus sequence for each reference SSU sequence, and the mapping is
repeated, followed by the estimations of probabilities.  The
iterations should usually stop when the reference sequences no longer
change from one iteration to the next.  Practically, 40 iterations is
usually sufficient for many samples.  Right now EMIRGE uses Bowtie
alignments internally, though in theory a different mapper could be
used.

EMIRGE was designed for Illumina reads in FASTQ format, from pipeline
version >= 1.3 (ASCII-offset of qual vals == 64)

DEPENDENCIES
------------------------------

EMIRGE expects the following programs to be installed and available in your path:

 -python (tested with version 2.6), with the following packages installed:
  -BioPython
  -Cython
  -pysam
  -scipy / numpy
 -usearch (www.drive5.com/usearch/ -- tested with usearch version 4.0.38).
 -samtools (http://samtools.sourceforge.net/ -- tested with verison 0.1.7)
 -bowtie (http://bowtie-bio.sourceforge.net/index.shtml -- tested with version 0.12.5)

INSTALLATION
------------------------------

After installing the dependencies listed above, type the following to build emirge:

  $ python setup.py build

To install (you may skip straight to this step), type the following as root, or with sudo:

  $ python setup.py install

You can also type the following for more options:
  $ python setup.py --help install

For example, to install instead to a location in your home directory
where you have permission to write, you might type something like:

  $ python setup.py install --prefix=$HOME/software

HELP
------------------------------

At the moment, there is very little documentation aside from running:

emirge --help

Once a run is completed, it is helpful to run emirge_rename_fasta.py
on the final iterations directory, for example:

emirge_rename_fasta.py iter.40 > renamed.fasta

Also see:

emirge_rename_fasta.py --help

Please feel free to contact me (csmiller@berkeley.edu) with any problems,
bug reports, or questions
  
CANDIDATE SSU DATABASE
------------------------------

SSU_candidate_db.fasta.gz is included with this distribution.  This was
made using Silva release 102 (http://www.arb-silva.de/).  Sequences
were clustered using uclust at 97% sequence identity, short and long
sequences were removed, and non-standard characters were changed to be
within {ACTG} (using utils/fix_nonstandard_chars.py).

You can use any reference SSU database with emirge, though this one is
recommended.  No matter your choice, you will need to first build a
bowtie index, with something like:
  $ bowtie-build SSU_candidate_db.fasta SSU_candidate_db_btindex
You might also consider changing the offrate (see
http://bowtie-bio.sourceforge.net/manual.shtml)



