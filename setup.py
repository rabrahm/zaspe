from distutils.core import setup, Extension

module =Extension('integration', sources = ['integration.c'],libraries=['gsl','gslcblas','m'],include_dirs=['/usr/local/include/', '/anaconda/lib/python2.7/site-packages/numpy/core/include/'], library_dirs=['/usr/local/lib/'])
setup(name = 'Funciones creadas por R. Brahm: Extensiones de C/Python', version = '1.0', ext_modules = [module])
