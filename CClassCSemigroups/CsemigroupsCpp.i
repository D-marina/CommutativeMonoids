%module CsemigroupsCpp 
%{ 
#define SWIG_FILE_WITH_INIT 
#include "CsemigroupsCpp.h" 
%} 

%include "std_vector.i"; 
namespace std { 
    %template(vector_long) vector<long>; 
    %template(vector_vector_long) vector<std::vector<long> >;
    %template(VecDouble) vector<double>;
    %template(VecVecdouble) vector< vector<double> >;
} 

%include "CsemigroupsCpp.h"; 