from distutils.core import setup, Extension

module =Extension('integration', sources = ['integration.c'],libraries=['gsl','gslcblas','m'])
setup(name = 'Funciones creadas por R. Brahm: Extensiones de C/Python', version = '1.0', ext_modules = [module])
