#!/usr/bin/env python
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


# In[7]:


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
        self.fNumber = 0
        self.d = 0
        self.NS = 0
        self.gaps = []
        self.pFrobenius = []
        self.t = 0
        
    # Frobenius Number
    fNumber = 0
    # Minimum value of Delta(S)
    d = 0
    # Bound for Delta(S) periodic
    NS = 0
    # Bound for Delta_nu periodic
    N0 = 0
    # Gaps of the semigroup
    gaps = []
    # Pseudo-Frobenius
    pFrobenius = []
    # Number of pseudo-Frobenius Elements
    t = 0
    
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
        print(self.generators,self.eDimension,self.d)
        self.NS = ComputeNs(self.generators,self.eDimension,self.d)
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
        if self.N0 > n:
            return Delta(Nu(self.generators,n,self.eDimension))
        return ComputeDeltaNu(self.generators,n,self.eDimension,self.NS,self.N0)
        
    def W(self,n):
        return W(self.generators,n,self.eDimension)
    
    def L(self,x):
        return L(self.generators,x,self.eDimension)
    
    def Nu(self,n):
        return Nu(self.generators,n,self.eDimension)
    
    ###################################################################################
    
    def SminusIthMinimalGenerator(self,i):
        '''
        Return the numerical semigroup S minus its ith minimal generator.
        '''
        saux=list(self.generators)
        x=saux[i]
        if saux==[1] and x==1:
            return NumericalSemigroup([2,3])
        saux.remove(x)
        #print(saux)
        saux= saux + [ x+y for y in saux ]
        saux=saux + [2*x,3*x]
        #print(saux)
        return NumericalSemigroup(saux)
    def Children(self):
        '''
        This function returns the children of a numerial semigroup.
        If S is a numerical semigroup, its children are the numerical semigroups S' verifying that S\S' 
        has cardinality 1 and the element in this set is a minimal generator of S greater than the Frobenius 
        number of S.
        ns=NumericalSemigroup([2,3])
        ns.Children()
        ns.Children()
        '''
        sg,nf=self.generators,self.FrobeniusNumber()
        familia=[]
        n=len(sg)
        for i in range(n):
            if nf < sg[i]:
                SNaux=self.SminusIthMinimalGenerator(i)
                familia=familia+[SNaux]
        return familia
    
    def Gaps(self):
        if self.gaps != []:
            return self.gaps
        m = min(self.generators)
        control = 0
        i = 1
        while control < m:
            if self.Belongs(i):
                control = control +1
            else:
                self.gaps.append(i)
                control = 0
            i = i+1
        return self.gaps
    
    def PseudoFrobenius(self):
        if self.pFrobenius != []:
            return self.pFrobenius
        if self.gaps == []:
            self.Gaps()
        pf = []
        for x in self.gaps:
            isPF = True
            for y in self.generators:
                if not self.Belongs(x+y):
                    isPF = False
                    break
            if isPF:
                self.pFrobenius.append(x)
        self.t = len(self.pFrobenius)
        return self.pFrobenius
    
    # Función auxiliar para la descomposición de irreducibles.
    def ComputeBP(self):
        aux = []
        for x in self.pFrobenius:
            if x > max(self.pFrobenius)/2:
                aux.append(x)
        return aux
    
    def IsIrreducible(self):
        if self.gaps == []:
            self.Gaps()
        if self.pFrobenius == []:
            self.PseudoFrobenius()
        bp = self.ComputeBP()
        if len(bp) == 1:
            return True
        return False
    
    # Función auxiliar para la descomposición de irreducibles.
    def Decomposition(self):
        if self.IsIrreducible():
            return [self]
        if self.gaps == []:
            self.Gaps()
        if self.pFrobenius == []:
            self.PseudoFrobenius()
        bp = self.ComputeBP()
        return [NumericalSemigroup(self.generators+[x]) for x in bp]
    
    def DecomposeIrreducible(self):
        candidatos = [self]
        end = False
        control = 0
        while end == False and control<20:
            candidatos2 = []
            for x in candidatos:
                candidatos2 = candidatos2 + x.Decomposition()
            aux = DeleteDuplicates([x.generators for x in candidatos2])
            candidatos = [NumericalSemigroup(x) for x in aux]
            end = all([x.IsIrreducible() for x in candidatos])
            control = control+1
        # Ahora calculamos los minimales
        if self.pFrobenius == []:
            self.PseudoFrobenius()
        bp = self.ComputeBP()
        minimales = []
        for x in candidatos:
            if max(x.PseudoFrobenius()) in bp:
                minimales.append(x)
        return minimales
    
    def Type(self):
        if self.t != 0:
            return self.t
        if self.pFrobenius == []:
            self.PseudoFrobenius()
            return self.t
        self.t = len(self.pFrobenius)
        return self.t
        

