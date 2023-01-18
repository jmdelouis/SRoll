from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

#import os
import numpy # only used for retrieving the numpy header directory
import os
import re


# The 'PYTHON_BASE' env variable is define in Makefile
#print("[DEBUG] $PYTHON_BASE='%s'" % os.environ['PYTHON_BASE'])

# Construct the path to the required numpy headers
#numpy_include_dir = os.environ['PYTHON_BASE'] + "/lib/python2.7/site-packages/numpy/core/include"
numpy_include_dir = numpy.get_include()
print("numpy_include_dir='%s'" % numpy_include_dir)

lib_file_full_path = "nodmclib/nodmclib.pyx"

extensions = [
  Extension("nodmclib.nodmclib", [lib_file_full_path],
    include_dirs = ["../NO_DMC_LIB", numpy_include_dir],
    libraries = ["no_dmc"],
    library_dirs = ["../Linux-x86_64"]),
]

def read_file(*paths):
  here = os.path.dirname(os.path.abspath(__file__))
  with open(os.path.join(here, *paths)) as f:
    return f.read()

def get_version():
  version_file = read_file(lib_file_full_path)
  version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]",
  version_file, re.M)
  if version_match:
    return version_match.group(1)
  raise RuntimeError("Unable to find version string.")


setup(
  name = 'nodmclib',
  description = 'Python wrapper for libno_dmc.so',
  version = get_version(),
  author = 'Christian Madsen',
  author_email = 'madsen@iap.fr',
  download_url = 'http://cvs.planck.fr/cvs/Level2/Lib_pkg/NO_DMC_LIB/pynodmclib/',
  ext_modules = cythonize(extensions),
)
