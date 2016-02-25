==================================================
EMIRGE -- :WARNING: Experimental Version :WARNING:
==================================================

WARNING
.. image:: https://travis-ci.org/epruesse/EMIRGE.svg?branch=master
   :target: https://travis-ci.org/epruesse/EMIRGE

-------------------------------------------------------------------------------
Expectation-Maximization Iterative Reconstruction of Genes from the Environment
-------------------------------------------------------------------------------


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
change from one iteration to the next.  Practically, 40-80 iterations
is usually sufficient for many samples.  Right now EMIRGE uses Bowtie
alignments internally, though in theory a different mapper could be
used.

EMIRGE is licensed_ under the GPLv3.

.. contents:: Table of Contents
.. .. section-numbering::



Publications
------------

If you use EMIRGE in your work, please cite these manuscripts as
appropriate.

Miller CS, Baker BJ, Thomas BC, Singer SW, Banfield JF (2011) *EMIRGE:
reconstruction of full-length ribosomal genes from microbial community
short read sequencing data.* Genome biology 12:R44.
`doi:10.1186/gb-2011-12-5-r44`__.

__ http://doi.org/10.1186/gb-2011-12-5-r44

Miller CS, Handley KM, Wrighton KC, Frischkorn KR, Thomas BC, Banfield
JF (2013) *Short-Read Assembly of Full-Length 16S Amplicons Reveals
Bacterial Diversity in Subsurface Sediments.* PloS one
8:e56018. `doi:10.1371/journal.pone.0056018`__.

__ http://doi.org/10.1371/journal.pone.0056018


Installation
------------

Before installing EMIRGE, please make sure the following tools are
installed:

 - Python (tested with 2.7)
 - Pip (https://pip.pypa.io/en/latest/installing.html)
 - usearch 6.0.203 or vsearch >= 1.1
 - samtools (http://samtools.sourceforge.net/ -- tested with verison 0.1.18)
 - bowtie (http://bowtie-bio.sourceforge.net/index.shtml -- tested
   with version 0.12.7 and 0.12.8)

EMIRGE itself is most easily installed using pip. For a system wide installation run::

  sudo pip install emirge

To install in your user account only, run::

  pip install --user emirge

(The above will typically install the emirge executables to
``~/.local/bin``, so you need to add that directory to your PATH or
call ``emirge.py`` by its full path).

You can also install directly from the github repository::

  pip install https://github.com/csmiller/EMIRGE


Running EMIRGE
--------------

At the moment, there is no manual, though running the following is
helpful::

  emirge.py --help


EMIRGE Output
~~~~~~~~~~~~~

Once an EMIRGE run is completed, run ``emirge_rename_fasta.py`` on the
final iterations directory, for example::

  emirge_rename_fasta.py iter.40 > renamed.fasta

Also see::

  emirge_rename_fasta.py --help

Running ``emirge_rename_fasta.py`` will provide you with a fasta file with
EMIRGE output.  Dissecting a single example header::

  >3326|AF427479.1.1480_m01 Prior=0.000367 Length=1480 NormPrior=0.000414
    1        2           3        4            5                6

1. The internal EMIRGE ID -- unique for each sequence
2. The accession number of the starting candidate sequence
3. an optional suffix indicating this sequence was split out from
   another due to evidence in the mapping reads of 2 or more
   "strains."
4. The Prior, or abundance estimate (used in original publication)
5. The length of the sequence
6. The length-normalized abundance estimate (anecdotally, this is
   sometimes more accurate if there are lots of different sequence
   lengths)


Candidate SSU database
----------------------

You can download a standard candidate SSU database by running the
following command::

  python emirge_download_candidate_db.py

This script is included with EMIRGE.  The current version of this
database was made using Silva_ release SSURef_NR_111. Sequences were
clustered using uclust at 97% sequence identity, short and long
sequences were removed, and non-standard characters were changed to be
within {ACTG} (using ``utils/fix_nonstandard_chars.py``).

.. _Silva: http://www.arb-silva.de/

You can use any reference SSU database with emirge, though this one is
recommended.  No matter your choice, you should run
``utils/fix_nonstandard_chars.py`` on your fasta file.  You will also need
to first build a bowtie index, with something like::

  bowtie-build SSU_candidate_db.fasta SSU_candidate_db_btindex

You might also consider changing the offrate (see
http://bowtie-bio.sourceforge.net/manual.shtml).


Getting Help
------------

There is a google group (similar to a mailing list) for asking
questions about EMIRGE: https://groups.google.com/group/emirge-users

Also, there is some additional information (including a Frequently
Asked Questions section) on the github wiki:
https://github.com/csmiller/EMIRGE/wiki

Although I encourage use of the google group due to increased volume
of support emails, please feel free to contact me directly
(christopher.s.miller@ucdenver.edu) with any problems, bug reports, or
questions.


Other
-----

** A note about single-end sequencing:

EMIRGE was designed for and tested on paired-end sequencing reads.
However, you can now use EMIRGE on single-end reads as well: simply
omit the -2 parameter.  Although I have done some basic testing on
single-end reads, runs with single reads have NOT been as extensively
tested as runs with paired reads.  Please let me know how it works for
you if you try EMIRGE with single-end reads.


..  _licensed:

COPYRIGHT AND LICENSE
---------------------

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
