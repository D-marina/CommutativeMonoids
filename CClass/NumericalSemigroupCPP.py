
# coding: utf-8

# In[2]:


import monoidsCpp
from monoidsCpp import *


# In[47]:


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
        self.fNumber = FrobeniusNumber(self.generators)
        return self.fNumber
    
    # This function returns the factorizations of an element
    def Factorizations(self,x):
        return FSolve(self.generators,x)
    
    # This function check if a number is in the semigroup
    def Belongs(self,x):
        if self.fNumber == 0:
            FrobeniusNumber(self.generators)
        return Belong(self.generators,x,self.fNumber)
    
    # This function compute the minimum of Delta(S)
    def ComputeMinimumDeltaS(self):
        if self.d != 0:
            return self.d
        self.d = ComputeD(self.generators)
        return self.d
        
    def ComputeNS(self):
        if self.NS != 0:
            return self.NS
        if self.d == 0:
            self.d = ComputeD(self.generators)
        self.NS = ComputeNs(self.generators)
        return self.NS
    
    def ComputeN0(self):
        if self.N0 != 0:
            return self.N0
        if self.NS == 0:
            self.NS = ComputeNs(self.generators)
        self.N0 = ComputeN0(self.generators,self.NS)
        return self.N0
    
    def DeltaNu(self,n):
        if self.NS == 0:
            self.NS = ComputeNs(self.generators)
        if self.N0 == 0:
            self.N0 = ComputeN0(self.generators,self.NS)
        if self.N0 > n:
            return Delta(Nu(self.generators,n))
        return ComputeDeltaNu(self.generators,n)
        
    def W(self,n):
        return W(self.generators,n)
    
    def L(self,x):
        return L(self.generators,x)
    
    def Nu(self,n):
        return Nu(self.generators,n)
    
    ###################################################################################
    
    def RemoveGenerator(self,i):
        return SminusIthMinimalGenerator(self.generators, i)
        
    def Descendants(self):
        if self.fNumber == 0:
            self.fNumber = FrobeniusNumber(self.generators)
        return Children(self.generators,self.fNumber)

# In[48]:


ns = NumericalSemigroup([4,5,6,7])


# In[49]:


ns.generators


# In[50]:


ns.Descendants()


# In[51]:


ns.fNumber


# In[52]:


ns.RemoveGenerator(0)

