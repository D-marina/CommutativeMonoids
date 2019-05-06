%module monoidsCpp 
%{ 
#define SWIG_FILE_WITH_INIT 
#include "monoidsCpp.h" 
%} 

%include "std_vector.i"; 
namespace std { 
    %template(vector_long) vector<long>; 
    %template(vector_vector_long) vector<std::vector<long> >;
} 

%include "monoidsCpp.h"; 

