import distutils 
from distutils.core import setup, Extension 

monoids = Extension("_monoids",sources = ['monoidsCpp.i','monoidsCpp.cxx',],swig_opts = ['-Wall','-c++'],) 

setup( name = 'monoids', ext_modules = [monoids,],  ) 
