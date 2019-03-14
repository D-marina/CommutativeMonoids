import numpy as np
from numpy import *
import sympy
from sympy import *

# Function for computing the Frobenius Number of a semigroup.

def FrobeniusNumber(generators,eDimension=0):
    if eDimension == 0:
        eDimension = len(generators)
    lGen =[generators[i] for i in range(eDimension)]
    h=t=a=generators[0]
    b=lGen
    Q = [0 for i in range(a)]
    S = [a*lGen[-1] for i in range(a*lGen[-1])]
    S[a-1]=0
    P = [i for i in range(a)]
    P[a-1] = len(lGen)
    while(h != 0):
        QHaux = Q[h-1]
        v = h
        Q[h-1] = 0
        if(h == t):
            h = 0
        else:
            h = QHaux
        for j in range(1,P[v-1]+1):
            e = (b[j-1]+v) % a
            w = b[j-1] + S[v-1]
            if(w < S[e-1] and e != 0 ):
                S[e-1] = w
                P[e-1] = j
                if(Q[e-1] == 0):
                    if(h == 0 ):
                        h = e
                        Q[e-1] = e
                        t = e
                    else:
                        Q[t-1] = e
                        Q[e-1] = e
                        t = e
    aux = [S[i] for i in range(0,a)]
    fNumber = max(aux)-a
    return(fNumber)

# Function to compute integer solutions of a Diophantine Equation

def FSolve(lgen,x,dim=0,onlyFirst=True):
    posAModificar=0
    sumando=True
    if dim == 0:
        dim = len(lgen)
    ceros=np.array([0 for i in range(dim)],dtype=np.int)
    xaux=x
    tuplaActual=array([0 for i in range(dim)],dtype=np.int)  
    if x==0:
        return tuplaActual
    soluciones=[]
    while( not( np.all(tuplaActual==ceros) and not(sumando) ) ):
        if sumando:
            if xaux>=lgen[0]:
                tuplaActual[posAModificar]=tuplaActual[posAModificar]+1
                xaux=xaux-lgen[posAModificar]
                if xaux==0:
                    if onlyFirst:
                        return tuplaActual
                    else:
                        soluciones.append(array(tuplaActual))
                        sumando=False
                        continue
            if xaux!=0 and xaux<lgen[0]:
                sumando=False
                continue
        else:
            if ( posAModificar == dim-1 ) and ( tuplaActual[posAModificar] > 0 ):
                xaux=tuplaActual[posAModificar]*lgen[posAModificar]+xaux
                tuplaActual[posAModificar]=0
                posAModificar=posAModificar-1
                continue
            if ( posAModificar < dim-1 ) and ( tuplaActual[posAModificar] > 0 ):
                xaux=xaux+lgen[posAModificar]
                tuplaActual[posAModificar]=tuplaActual[posAModificar]-1
                posAModificar=posAModificar+1
                sumando=True
                continue
            if ( posAModificar < dim-1 ) and ( tuplaActual[posAModificar] == 0 ):
                posAModificar=np.max( [i for i in range(dim) if tuplaActual[i]!=0] )
                xaux=lgen[posAModificar]+xaux
                tuplaActual[posAModificar]=tuplaActual[posAModificar]-1
                sumando=True
                posAModificar=posAModificar+1
                continue
    return [list(xx) for xx in soluciones]

# Function for computing the minimal set of generators of a numerical semigroup.

def smgS(generators):
    smgS=array([],dtype=np.int)
    generators.sort()
    sgS=np.unique(generators)
    if len(sgS)==1:
        return sgS
    smgS=np.append(smgS,sgS[0])
    for i in range(1,len(sgS)):
        if sgS[i] % smgS[0] != 0:
            smgS=np.append(smgS,sgS[i])
            break
    for j in range(i+1,len(sgS)):
        if len(FSolve(smgS,sgS[j]))==0:
            smgS=np.append(smgS,sgS[j])
    return list(smgS)