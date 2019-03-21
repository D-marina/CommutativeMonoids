
# coding: utf-8

# In[1]:

import numpy as np
from numpy import *
import sympy
from sympy import *


# In[2]:

import integerSmithNormalFormAndApplications
from integerSmithNormalFormAndApplications import *


# In[3]:

import auxiliars
from auxiliars import *


# In[4]:

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
    # Minimum value of Delta(S)
    d = 0
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
    
    # This function compute the minimum of Delta(S)
    def ComputeMinimumDeltaS(self):
        if self.d != 0:
            return self.d
        self.d = ComputeD(self.generators,self.eDimension)
        return self.d
        
    def ComputeNS(self):
        if self.NS != 0:
            return self.NS
        if self.d == 0:
            self.d = ComputeD(self.generators,self.eDimension)
        self.NS = ComputeNs(self.generators,self.d,self.eDimension)
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
        laux=list(set( [ sum([x[i]*smg[i] for i in range(dim) ]) for x in f1(dim,n) ]))
        laux.sort()
        return laux
    
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
    def SminusIthMinimalGenerator(self,i):
        '''
        Return the numerical semigroup S minus its ith minimal generator.
        '''
        return None
    def Children(self):
        '''
        This function returns the children of a numerial semigroup.
        If S is a numerical semigroup, its children are the numerical semigroups S' verifying that S\S' 
        has cardinality 1 and the element in this set is a minimal generator of S greater than the Frobenius 
        number of S.
        '''
        return []


# In[5]:

#ns = NumericalSemigroup([5,7])


# In[6]:

#ns.ComputeNS()


# In[7]:

ns2 = NumericalSemigroup([15,17,27,35])


# In[8]:

ns2.ComputeNS()


# In[ ]:

ns.NS


# In[ ]:

ns.ComputeNS()


# In[ ]:

ns = NumericalSemigroup([6,12,15])


# In[ ]:

#set.union(*[set(ns.L(x)) for x in ns.W(10)])


# In[ ]:

#ns.nu(10)


# In[ ]:

#ns.ComputeNS()


# In[ ]:

#ns.ComputeN0()


# In[ ]:

#ns.DeltaNu(70)


# In[ ]:

#ns.FrobeniusNumber()


# In[ ]:

#ns.ComputeN0()


# In[ ]:

#ns.Factorizations(30)


# In[ ]:

#ns.Belongs(4)

