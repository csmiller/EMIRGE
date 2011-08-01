# typical compilation: python setup.py build_ext --inplace

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

from numpy import get_include
numpy_include_dir = get_include()

ext_modules = [Extension("_emirge", ["_emirge.pyx", "_emirge_C.h"], include_dirs=[numpy_include_dir])]

setup(
    name = '_emirge',
    description = 'emirge helper functions implemented in Cython for speed',
    cmdclass = {'build_ext': build_ext},
    ext_modules = ext_modules,

    scripts = ["emirge.py"],

    author="Christopher Miller",
    author_email="csmiller@gmail.com",
    version="0.0.1",
    license="GPLv3"
)
