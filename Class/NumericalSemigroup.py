#!/usr/bin/env python
# coding: utf-8

# In[4]:


import numpy as np
from numpy import *
import sympy
from sympy import *


# In[5]:


import integerSmithNormalFormAndApplications
from integerSmithNormalFormAndApplications import *


# In[6]:


import auxiliars
from auxiliars import *


# In[1]:


class NumericalSemigroup:
    '''
    This class contains all methods and properties needs
    for working with a Numerical Semigroup.
    
    In order to use it you have to introduce the generators. For example S = NumericalSemigroup([3,4,5]).
    
    Functions:
        * FrobeniusNumber() returns the frobenius number of S.
        * Factorizations(x) returns the factorizations of x in S.
        * Belongs(x) returns True or False if x is or is not in S.
        * ComputeNS() returns a bound for the periodicity of Delta(S)
    '''
    
    def __init__(self,generators):
        self.generators = smgS(generators)
        self.multiplicity = self.generators[0]
        self.eDimension = len(self.generators)
        
    # Frobenius Number
    fNumber = 0
    # Bound for Delta(S) periodic
    NS = 0
    # Bound for Delta_nu periodic
    N0 = 0
    
    # This function gives us the Frobenius Number of a semigroup
    def FrobeniusNumber(self):
        if self.fNumber != 0:
            return self.fNumber
        self.fNumber = FrobeniusNumber(self.generators,self.eDimension)
        return self.fNumber
    
    # This function returns the factorizations of an element
    def Factorizations(self,x):
        return FSolve(self.generators,x,self.eDimension,False)
    
    # This function check if a number is in the semigroup
    def Belongs(self,x):
        if self.fNumber == 0:
            FrobeniusNumber(self.generators,self.eDimension)
        return Belong(self.generators,x,self.multiplicity,self.fNumber)
        
    def ComputeNS(self):
        if self.NS != 0:
            return self.NS
        self.NS = ComputeNs(self.generators,self.eDimension)
        return self.NS
    
    def ComputeN0(self):
        if self.N0 != 0:
            return self.N0
        if self.NS == 0:
            self.NS = ComputeNs(self.generators,self.eDimension)
        self.N0 = ComputeN0(self.generators,self.eDimension,self.NS)
        return self.N0
    
    def DeltaNu(self,n):
        if self.NS == 0:
            self.NS = ComputeNs(self.generators,self.eDimension)
        if self.N0 == 0:
            self.N0 = ComputeN0(self.generators,self.eDimension,self.NS)
        return ComputeDeltaNu(self.generators,n,self.eDimension,self.NS,self.N0)
        
    def W(self,n):
        smg=self.generators
        dim=self.eDimension
        return list(set( [ sum([x[i]*smg[i] for i in range(dim) ]) for x in f1(dim,n) ]))
    
    def L(self,x):
        l1=self.Factorizations(x)
        l2=[sum(y) for y in l1]
        l2.sort()
        return l2
    
    def nu(self,n,debug=False):
        waux=self.W(n)
        if debug:
            print("W of ",n,":",waux)
        longAux=list(set.union( *[set(self.L(x)) for x in waux] ))
        longAux.sort()
        if debug:
            print([self.L(x) for x in waux])
        return (longAux)


# In[8]:


ns = NumericalSemigroup([3,7,9])


# In[9]:


#ns.ComputeNS()


# In[10]:


#ns.ComputeN0()


# In[13]:


#ns.DeltaNu(70)


# In[ ]:


#ns.FrobeniusNumber()


# In[ ]:


#ns.ComputeN0()


# In[ ]:


#ns.Factorizations(30)


# In[ ]:


#ns.Belongs(4)

