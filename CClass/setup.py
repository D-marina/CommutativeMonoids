import distutils 
from distutils.core import setup, Extension 

monoids = Extension("_monoids",sources = ['monoids.i','monoids.cxx',],swig_opts = ['-Wall','-c++'],) 

setup( name = 'monoids', ext_modules = [monoids,],  ) 
