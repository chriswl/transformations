from distutils.core import setup, Extension
import numpy
setup(name='_transformations',
      author='Christoph Gohlke',
      modules=['transformations'],
      ext_modules=[Extension('transformations._transformations', ['transformations/transformations.c'],
                             include_dirs=[numpy.get_include()])])
