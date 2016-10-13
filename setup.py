#!/usr/bin/env python

"""
Installation script for EMIRGE

Large parts of this have been borrowed from the rlpy project
(https://github.com/rlpy/rlpy) which have been permitted for use under the
BSD license.
"""

import sys

try:
    from setuptools import setup, Command, find_packages
except ImportError:
    sys.exit("EMIRGE requires \"setuptools\" for installation.")

from distutils.extension import Extension
from distutils.command.build import build
from distutils.command.sdist import sdist
from distutils.command.build_ext import build_ext as _build_ext
import os
import shutil
import pkg_resources
from os.path import join as pjoin


version = '0.60.4'


try:
    from Cython.Distutils import build_ext as _build_ext
    cython = True
    print("Found Cython.")
except ImportError:
    cython = False
    print("Cython not installed.")


class build_ext(_build_ext):
    """extends build_ext to add numpy include dir at build time
    Using numpy.get_includes() would require numpy to have already been
    installed when this script is run. By deferring the time at which the
    path is obtained, we should be able to install numpy with Emirge in one go.
    """
    def build_extensions(self):
        numpy_incl = pkg_resources.resource_filename('numpy', 'core/include')
        for ext in self.extensions:
            if hasattr(ext, 'include_dirs') and \
                            numpy_incl not in ext.include_dirs:
                ext.include_dirs.append(numpy_incl)
        _build_ext.build_extensions(self)


class CheckingBuildExt(build_ext):
    """extend build_ext to check that all Cython generated c sources are
    present before going on to compile them. This is mainly to get a clearer
    error message if neither cython nor the c files are available.
    """
    @staticmethod
    def check_cython_extensions(extensions):
        for ext in extensions:
            for src in ext.sources:
                if not os.path.exists(src):
                    raise Exception("""Cython-generated file '%s' not found.
                Cython is required to compile EMIRGE from a development branch.
                Please install Cython or download a release package of EMIRGE.
                """ % src)

    def build_extensions(self):
        self.check_cython_extensions(self.extensions)
        build_ext.build_extensions(self)


class CythonCommand(build_ext):
    """Custom distutils command subclassed from Cython.Distutils.build_ext
    to compile pyx->c, and stop there. All this does is override the
    C-compile method build_extension() with a no-op."""
    def build_extension(self, ext):
        pass


class DummyBuildSrc(Command):
    """ numpy's build_src command interferes with Cython's build_ext.
    """
    user_options = []

    def initialize_options(self):
        self.py_modules_dict = {}

    def finalize_options(self):
        pass

    def run(self):
        pass


class CheckSDist(sdist):
    """Custom sdist that ensures Cython has compiled all pyx files to c."""
    _pyxfiles = ["pykseq/pykseq.pyx",
                 "_emirge.pyx",
                 "_emirge_amplicon.pyx"]

    def initialize_options(self):
        sdist.initialize_options(self)

    def run(self):
        if 'cython' in cmdclass:
            self.run_command('cython')
        else:
            for pyxfile in self._pyxfiles:
                cfile = pyxfile[:-3] + 'c'
                cppfile = pyxfile[:-3] + 'cpp'
                msg = "C-source file '%s' not found. Run 'setup.py cython' " \
                      "before sdist." % cfile
                assert os.path.isfile(cfile) or os.path.isfile(cppfile), msg
        sdist.run(self)


class CleanCommand(Command):
    """Custom distutils command to clean the .so and .pyc files."""
    user_options = [("all", "a", "")]

    def initialize_options(self):
        self.all = True
        self._clean_me = []
        self._clean_trees = []
        self._clean_exclude = []

        for root, dirs, files in list(os.walk('Emirge')):
            for f in files:
                if f in self._clean_exclude:
                    continue
                if os.path.splitext(f)[-1] in ('.pyc', '.so', '.o',
                                               '.pyo',
                                               '.pyd', '.c', '.orig'):
                    self._clean_me.append(pjoin(root, f))
            for d in dirs:
                if d == '__pycache__':
                    self._clean_trees.append(pjoin(root, d))
        for d in ('build',):
            if os.path.exists(d):
                self._clean_trees.append(d)

    def finalize_options(self):
        pass

    def run(self):
        for clean_me in self._clean_me:
            try:
                os.unlink(clean_me)
            except Exception:
                pass
        for clean_tree in self._clean_trees:
            try:
                shutil.rmtree(clean_tree)
            except Exception:
                pass


cmdclass = {'clean': CleanCommand,
            'build': build,
            'sdist': CheckSDist,
            'build_ext': CheckingBuildExt,
            }

if cython:
    cmdclass['cython'] = CythonCommand
else:
    cmdclass['build_src'] = DummyBuildSrc


extensions = [
    Extension("pykseq", ["pykseq/pykseq.pyx"],
              libraries=["z"]),
    Extension("_emirge", ["_emirge.pyx"],
              extra_compile_args=["-O3"]),
    Extension("_emirge_amplicon", ["_emirge_amplicon.pyx"],
              extra_compile_args=["-O3"],
              include_dirs=['./pykseq'],
              libraries=["z"])
    ]


def no_cythonize(extensions, **_):
    for extension in extensions:
        sources = []
        for sfile in extension.sources:
            path, ext = os.path.splitext(sfile)
            if ext in ('.pyx', '.py'):
                if extension.language == 'c++':
                    ext = '.cpp'
                else:
                    ext = '.c'
                sfile = path + ext
            elif ext in '.pxd':
                continue
            sources.append(sfile)
        extension.sources[:] = sources
    return extensions


if cython and os.path.exists("_emirge_amplicon.pyx"):
    from Cython.Build import cythonize
    extensions = cythonize(extensions, include_path=['pykseq'])
else:
    extensions = no_cythonize(extensions)


setup(
    name='EMIRGE',
    version=version,
    description="EMIRGE reconstructs full length sequences from short "
                "sequencing reads",
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
        "License :: OSI Approved :: GNU General Public License v3 or later "
        "(GPLv3+)",
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
    data_files=["pykseq/kseq.h","_emirge_C.h"],
    ext_modules=extensions,
    cmdclass=cmdclass,
    packages=find_packages(exclude=['tests', 'tests.*']),
    license="GPLv3+",
    keywords=["rRNA", "EM"],
    install_requires=["numpy", "pysam", "scipy", "biopython"],
    setup_requires=["numpy"]
)

print ""
print "NOTE:"
print "To download a standard candidate SSU database to use with EMIRGE, run"
print "python emirge_download_candidate_db.py"


