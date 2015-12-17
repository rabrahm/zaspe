import numpy as np
from distutils.core import setup, Extension

module =Extension('Cfunctions', sources = ['Cfunctions.c'],libraries=['gsl','gslcblas','m'],include_dirs=[np.get_include(),'/usr/local/include/', '/anaconda/lib/python2.7/site-packages/numpy/core/include/'], library_dirs=['/usr/local/lib/'])
setup(name = 'Funciones creadas por Nestor: Extensiones de C/Python', version = '1.0', ext_modules = [module])
