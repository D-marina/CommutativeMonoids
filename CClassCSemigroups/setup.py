import distutils 
from distutils.core import setup, Extension 

CsemigroupsCpp = Extension("_CsemigroupsCpp",sources = ['CsemigroupsCpp.i','CsemigroupsCpp.cxx',],swig_opts = ['-Wall','-c++'],) 

setup( name = 'CsemigroupsCpp', ext_modules = [CsemigroupsCpp,],  ) 
