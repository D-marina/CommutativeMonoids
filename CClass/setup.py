import distutils 
from distutils.core import setup, Extension 

Csemigroups = Extension("_Csemigroups",sources = ['CsemigroupsCpp.i','Csemigroups.cxx',],swig_opts = ['-Wall','-c++'],) 

setup( name = 'Csemigroups', ext_modules = [Csemigroups,],  ) 
