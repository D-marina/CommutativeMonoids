import distutils 
from distutils.core import setup, Extension 

monoidsCpp = Extension("_monoidsCpp",sources = ['monoidsCpp.i','monoidsCpp.cxx',],swig_opts = ['-Wall','-c++'],) 

setup( name = 'monoidsCpp', ext_modules = [monoidsCpp,],  ) 
