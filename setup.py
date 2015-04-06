#!/usr/bin/env python

import sys

try:
    from setuptools import setup
except ImportError:
    sys.exit("EMIRGE requires \"setuptools\" for installation.")

from distutils.extension import Extension

#def clean_cython(mod):
#    torm = [ modname[:-3] + ext
#             for modname in ext.sources
#             for ext in ('c','cpp','so')
#             if modname.endswith('.pyx') ]
#    for f in torm:
#        if os.path.exists(f):
#            os.unlink(f)

# class to evaluate a list lazily (to defer numpy/cython)
class lazy_eval_list(list):
    def __init__(self, callback):
        self._list, self.callback = None, callback

    def c_list(self):
        if self._list is None:
            self._list = self.callback()
        return self._list

    def __iter__(self):
        for e in self.c_list():
            yield e

    def __getitem__(self, ii):
        return self.c_list()[ii]

    def __len__(self):
        return len(self.c_list())


def extensions():
    from Cython.Build import cythonize
    from numpy import get_include
    numpy_include_dir = get_include()
    ext_modules = [
        Extension("*", ["Emirge/*.pyx"],
                  include_dirs=[numpy_include_dir, "Emirge/"],
                  )
#        Extension("emirge.pykseq", ["emirge/pykseq.pyx"],
#                  include_dirs=['./emirge/']),
#        Extension("_emirge", ["_emirge.pyx"],
#                  include_dirs=[numpy_include_dir],
#                  extra_compile_args=["-O3"]),
#        Extension("_emirge_amplicon", ["_emirge_amplicon.pyx"],
#                  libraries=["z"],
#                  include_dirs=[numpy_include_dir, './pykseq'],
#                  library_dirs=['./pykseq'],
#                  extra_compile_args=["-O3"])
    ]
    return cythonize(ext_modules)

setup(
    name='EMIRGE',
    version="0.6.1",
    description="EMIRGE reconstructs full length sequences from short sequencing reads",
    long_description="""
    EMIRGE: Expectation-Maximization Iterative Reconstruction of Genes
            from the Environment

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
    change from one iteration to the next.  Practically, 40-80 iterations is
    usually sufficient for many samples.  Right now EMIRGE uses Bowtie
    alignments internally, though in theory a different mapper could be
    used.

    EMIRGE was designed for Illumina reads in FASTQ format, from pipeline
    version >= 1.3

    """,
    classifiers=[
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Operating System :: POSIX",
        "Programming Language :: Cython",
        "Programming Language :: Python :: 2.7",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    author="Christopher Miller",
    author_email="christopher.s.miller@ucdenver.edu",
    url="https://github.com/csmiller/EMIRGE",
    scripts=["emirge.py",
             "emirge_amplicon.py",
             "emirge_rename_fasta.py"],
    ext_modules=lazy_eval_list(extensions),
    license="GPLv3+",
    keywords=["rRNA", "EM"],
    install_requires=["Cython", "numpy", "pysam", "scipy", "biopython"],
    setup_requires=["Cython", "numpy"]
)

print ""
print "NOTE:"
print "To download a standard candidate SSU database to use with EMIRGE, run"
print "python emirge_download_candidate_db.py"


