from distutils.core import setup, Extension

module =Extension('Cfunctions', sources = ['Cfunctions.c'],libraries=['gsl','gslcblas','m'])
setup(name = 'Funciones creadas por Nestor: Extensiones de C/Python', version = '1.0', ext_modules = [module])
