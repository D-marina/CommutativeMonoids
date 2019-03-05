import numpy as np
from numpy import *

def FrobeniusSolve(lgen,x,onlyFirst=True):
    '''
    FrobeniusSolve([55,71,99],1000,onlyFirst=False)
    '''
    posAModificar=0
    sumando=True
    dim = len(lgen)
    ceros=np.array([0 for i in range(dim)],dtype=np.int)
    xaux=x
    tuplaActual=array([0 for i in range(dim)],dtype=np.int)
    #print(xaux,tuplaActual)    
    if x==0:
        return tuplaActual
    soluciones=[]
    while( not( np.all(tuplaActual==ceros) and not(sumando) ) ):
        #print(tuplaActual)
        if sumando:
            if xaux>=lgen[0]:
                tuplaActual[posAModificar]=tuplaActual[posAModificar]+1
                xaux=xaux-lgen[posAModificar]
                if xaux==0:
                    if onlyFirst:
                        return tuplaActual
                    else:
                        #print(tuplaActual)
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
    #return array([],dtype=np.int)
    return soluciones