# ## We import the libraries used

import math
from math import gcd
from sympy import Matrix, eye, init_printing, pprint, zeros
import sympy
import itertools
from subprocess import call, PIPE, Popen
import json
from  numpy.random import randint
import numpy
from fractions import Fraction
from sympy import groebner
import numpy as np
from numpy import array

init_printing(use_latex='mathjax')

# Several simple functions, some of them not used. 

def lcm(a,b): return int(round(abs(a * b) / math.gcd(a,b))) if a and b else 0
def lcmL(l):
    '''
    >>> lcmL([5,12,13])
    '''
    if len(l)==1:
        return l[0]
    if len(l)==2:
        return lcm(l[0],l[1])
    else:
        aux=l[2:]+[lcm(l[0],l[1])]
        return lcmL( aux )

def egcd(a, b):
    x,y, u,v = 0,1, 1,0
    while a != 0:
        q, r = b//a, b%a
        m, n = x-u*q, y-v*q
        b,a, x,y, u,v = a,r, u,v, m,n
    g = b
    return g, x, y

def gcdL(l):
    '''
    >>> gcdL([12,34,22])
    '''
    if len(l)==1:
        return l[0]
    if len(l)==2:
        return gcd(l[0],l[1])
    else:
        aux=l[2:]+[gcd(l[0],l[1])]
        return gcdL(aux)

def gcdMatrix(A):
    return gcdL(A[:])

def sympyMatrix2numpyArray(m,tipo=np.int):
    '''
    >>> sympyMatrix2numpyArray(Matrix([[1,2],[2,3]]))
    '''
    nf,nc=m.shape
    return array(list(m),dtype=tipo).reshape(nf,nc)
def numpyArray2sympyMatrix(npM):
    '''
    >>> numpyArray2sympyMatrix( np.array([[1,2],[3,4]]) )
    '''
    return Matrix( list( map( lambda x:list(map(int,x)) , npM ) ) )


# ## Several functions used to check the returned values of the function `integerSmithNormalForm`.


def isDiagonal(A):
    nf,nc=A.shape
    for i in range(nf):
        for j in range(nc):
            if i!=j and A[i,j]!=0:
                return False
    return True

def isSeqDiagOfDivisible(A):
    nf,nc=A.shape
    k=min(nf,nc)
    if not isDiagonal(A):
        return False
    for i in range(1,k):
        if A[i,i] % A[i-1,i-1]!=0:
            return False
    return True


# ---
# ## Function `integerSmithNormalForm`

# Some auxilary functions



def posMinNonNullOfMatrix(A):
    nf,nc=A.shape
    l=A[:]
    m=min([x for x in l if x!=0])
    aux=l.index(m)   
    return (aux//nc,aux%nc)
def putAbsMinInCorner(A):
    nf,nc=A.shape
    R=eye(nf)
    C=eye(nc)
    Aux=Matrix(A)
    i,j=posMinNonNullOfMatrix(Aux.applyfunc(abs))
    R.row_swap(i,0)
    C.col_swap(j,0)
    Aux=R.multiply(Aux).multiply(C)
    if Aux[0,0]<0:
        Raux=eye(nf)
        Raux[0,0]=-1
        R=Raux.multiply(R)
    return (R,C)



def addRowIfNecesary(A):
    nf,nc=A.shape
    Aux=Matrix(A)
    c=Aux[0,0]
    k=None
    for i in range(1,nf):
        if Aux[i,:].applyfunc(lambda x:x%c)!=Matrix.zeros(1,nc):
            k=i
            break;
    Maux=Matrix.eye(nf)
    if k:
        Maux[0,:]=Maux[0,:]+Maux[k,:]
    return Maux


def makeZeroInFirstColumn(A,i):
    Aux=Matrix(A)
    nf,nc=Aux.shape
    m1=sympy.eye(nf)
    if Aux[0,0]<0:
        m1[0,:]=-m1[0,:]
    Aux=m1.multiply(Aux)
    if Aux[i,0]==0:
        return m1
    while Aux[i,0]!=0:
        q=int(Aux[i,0])//int(Aux[0,0])
        maux=sympy.eye(nf)
        maux[i,:]=maux[i,:]-q*maux[0,:]
        Aux=maux.multiply(Aux)
        m1=maux.multiply(m1)
        if Aux[i,0]!=0:
            maux=eye(nf)
            maux.row_swap(0,i)
            m1=maux.multiply(m1)
            Aux=maux.multiply(Aux)
    return m1

def makeZeroInFirstRow(A,i):
    RT=makeZeroInFirstColumn(A.T,i)
    return RT.T

def makeZeroFirstColumn(A):
    nf,nc=A.shape
    Aux=Matrix(A)
    R=Matrix.eye(nf)
    for i in range(1,nf):
        Raux=makeZeroInFirstColumn(Aux,i)
        Aux=Raux.multiply(Aux)
        R=Raux.multiply(R)
    return R

def makeZeroFirstRow(A):
    B=Matrix(A.T)
    CT=makeZeroFirstColumn(B)
    return CT.T
def makeZerosFirstRowColumn(A):
    Aux=Matrix(A)
    nf,nc=A.shape
    R1=eye(nf)
    C1=eye(nc)
    c0=Aux[1:,0]
    f0=Aux[0,1:]
    while(c0!=zeros(nf-1,1) or f0!=zeros(1,nc-1)):
        R=makeZeroFirstColumn(Aux)
        Aux=R.multiply(Aux)
        R1=R.multiply(R1)
        C=makeZeroFirstRow(Aux)
        Aux=Aux.multiply(C)
        C1=C.multiply(C1)
        c0=Aux[1:,0]
        f0=Aux[0,1:]
    return (R1,C1)


# The function `integerSmithNormalForm`


def integerSmithNormalForm(A):
    '''
    >>> A=Matrix([[-3,11,3],[-48,15,12]])
    >>> R,C=integerSmithNormalForm(A)
    >>> [R.multiply(A).multiply(C),R.det(),C.det()]
    '''
    nf,nc=A.shape
    Aux=Matrix(A)
    R=eye(nf)
    C=eye(nc)
    if A==zeros(nf,nc):
        return (R,C)
    else:
        t=True
        while(t):
            R1,C1=putAbsMinInCorner(Aux)
            Aux=R1.multiply(Aux).multiply(C1)
            R=R1.multiply(R)
            C=C.multiply(C1)
        
            R1,C1=makeZerosFirstRowColumn(Aux)
            Aux=R1.multiply(Aux).multiply(C1)
            R=R1.multiply(R)
            C=C.multiply(C1)

            R1=addRowIfNecesary(Aux)
            Aux=R1.multiply(Aux)
            R=R1.multiply(R)
            t= not (R1==eye(nf))
            
    Rm,Cm=integerSmithNormalForm(Aux[1:,1:])
    Raux=eye(nf)
    Caux=eye(nc)
    Raux[1:,1:]=Rm
    Caux[1:,1:]=Cm
    return (Raux.multiply(R),C.multiply(Caux))



# ## Function to compute the generators of a subgroup of $\mathbb Z^p$ from its defining equations

# #### Function to obtain a minimal system of generators from a set of homogeneus equations

def equationsToGeneratorsHomogeneusCase(A):
    '''
    >>> equationsToGeneratorsHomogeneusCase(Matrix([[5,-7,3,-2],[6,9,-10,1]]))
    '''
    nf,nc=A.shape
    R,C=integerSmithNormalForm(A)
    D=R.multiply(A).multiply(C)
    noNullOfD=[D[i,i] for i in range(min(D.shape)) if D[i,i]!=0]
    r=len(noNullOfD)
    nGen=nc-r
    Caux=C[:,-nGen:].T
    return Caux



# #### Function to obtain the minimal system from a set of equations


def equationsToGenerators(A,modulus):
    '''
    >>> A=Matrix([[1,-2,3,4],[6,8,-10,-16],[5,7,-9,2]])
    >>> mm=Matrix([2,3])
    >>> sG=equationsToGenerators(A,mm)
    >>> pprint(sG)
    >>> A.multiply(sG.T)
    '''
    Aux=Matrix(A)
    nfm,ncm=modulus.shape
    nf,nc=A.shape
    idaux=eye(nfm)
    for i in range(nfm):
        idaux[i,i]=-modulus[i]
    zaux=zeros(nf-nfm,nfm)
    idaux=idaux.col_join(zaux)
    Aux=Aux.row_join(idaux)
    #pprint(Aux)
    r=equationsToGeneratorsHomogeneusCase(Aux)
    #pprint(Aux.multiply(r.T))
    return r[:,0:nc]


def minimalSystemOfGenerators(sGen):
    '''
    Returns a minimal system of generators of the subgroup generated by 
    the rows of the sympy.matrix sGen.
    >>> minimalSystemOfGenerators(Matrix([[5,1,0],[1,2,-3],[6,3,-3]]))
    '''
    A=Matrix(sGen)
    R,C=integerSmithNormalForm(A)
    D=R.multiply(A).multiply(C)
    RA=R.multiply(A)
    noNullOfD=[D[i,i] for i in range(min(D.shape)) if D[i,i]!=0]
    r=len(noNullOfD)
    return RA[0:r,:]

# #### Function to compute the equations from a system of generators

def generatorsToEquations(sGen,lVars=None):
    '''
    >>> sGen=Matrix([[2,-4,8],[3,2,-1]])
    >>> s1=generatorsToEquations(sGen,[x,y,z])
    >>> s2=generatorsToEquations(sGen)
    >>> [s1,s2]
    '''
    A=Matrix(sGen)
    nf,nc=A.shape
    R,C=integerSmithNormalForm(A)
    D=R.multiply(A).multiply(C)
    noNullOfD=[D[i,i] for i in range(min(D.shape)) if D[i,i]!=0]
    lmod=noNullOfD+[0 for i in range(nc-len(noNullOfD))]
    n1=len([x for x in noNullOfD if x==1])
    lmod=lmod[n1:]
    if not lVars:
        return (C.T[n1:,:],Matrix(lmod))
    else:
        return (Matrix(lVars).T.multiply(C[:,n1:]).T,Matrix(lmod))
