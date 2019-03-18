import numpy as np
from numpy import *
import sympy
from sympy import *

import auxiliars
from auxiliars import *

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
    tuplaActual=np.array([0 for i in range(dim)],dtype=np.int)  
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
                        soluciones.append(np.array(tuplaActual))
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
    smgS=np.array([],dtype=np.int)
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

# Function for knowing if an element is in a numerical semigroup.

def Belong(generators,x,multiplicity=0,fNumber=0):
    if multiplicity == 0:
        multiplicity = generators[0]
    if fNumber == 0:
        fNumber = FrobeniusNumber(generators,len(generators))
    if x==0:
        return True
    if 0<x and x<multiplicity:
        return False
    if x in generators:
        return True
    if fNumber != 0 and x>fNumber:
        return True
    expression = FSolve(generators,x,len(generators))
    if len(expression)>0:
        return True
    else:
        return False
    

# This function compute the Ns of a numerical semigroup.

def ComputeNs(generators,dimension=0):
    # Calculamos la dimension
    a = generators
    if dimension == 0:
        dimension = len(a)
    # Creamos los simbolos ei
    e = symbols('e0:%d'%dimension)
    # Creamos el simbolo s
    s = symbols('s')
    # Calculamos d
    m = equationsToGeneratorsHomogeneusCase(Matrix([generators]))
    m1 = sympyMatrix2numpyArray(m)
    l = numpy.sum(m1,axis=1) # Esto es lo que fallaba porque l es de la forma [a b c] y deberia ser [a, b, c] con las comas, asi que lo transformo aqui.
    lista = []
    for i in range(len(l)):
        lista.append(l[i])
    d = gcdL(lista)
    # Calculo un vector con los mcd
    mcd = []
    for i in range(dimension-2):
        mcd.append(gcdL([a[i+1]-a[-1],-a[0]+a[-1],a[0]-a[i+1]]))
    # Calculo h
    h=d/(a[-1]-a[0])*(a[-1]*e[0]-a[0]*e[-1])
    # Calculo P2
    P2=s*(a[1]-a[-1])/(a[1]*(a[0]-a[-1]))*e[0]+s*(a[0]-a[1])/(a[1]*(a[0]-a[-1]))*e[-1]
    # Calculo Pp-1
    Ppm1=s*(a[-2]-a[-1])/(a[-2]*(a[0]-a[-1]))*e[0]+s*(a[0]-a[-2])/(a[-2]*(a[0]-a[-1]))*e[-1]
    # Calculo un vector con los qi
    Q = []
    for i in range(dimension-2):
        Q.append(1/mcd[i]*((a[i+1]-a[-1])*e[0]+(a[-1]-a[0])*e[i+1]+(a[0]-a[i+1])*e[-1]))
    # Calculo la suma de los Qi
    sumaQi = 0
    for i in range(dimension-2):
        sumaQi += Q[i]
    # Calculamos la última coordenada P2+h+∑qi
    primeraEcuacion = P2+h+sumaQi
    for i in range(dimension-1):
        primeraEcuacion = primeraEcuacion.subs(e[i],0)
    primeraEcuacion = primeraEcuacion.subs(e[dimension-1],1)
    sol1 = sympy.solve(primeraEcuacion,s)
    # Calculamos la primera coordenada de P(p−1)−h+∑qi
    ultimaEcuacion = Ppm1-h+sumaQi
    for i in range(1,dimension):
        ultimaEcuacion = ultimaEcuacion.subs(e[i],0)
    ultimaEcuacion = ultimaEcuacion.subs(e[0],1)
    sol2 = sympy.solve(ultimaEcuacion,s)
    NS = int(ceil(int(max(sol1,sol2)[0])))
    return NS    
    
# This function compute Lambda_1 for computing the Delta_nu

def CalcularLambda1(lgen,dimension = 0,Ns=0):
    if Ns == 0:
        ComputeNs(lgen,dimension)
    c1 = (lgen[-1]-lgen[-2])/lgen[-2]*Ns
    c4 = (lgen[0]/lgen[-2]-lgen[0]/lgen[-1]-lgen[0]/lgen[1]+1)*Ns
    return max(c1,c4)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
def f1(e, n):
    '''
    Function to compute the way to distribute n iqual objects in e different positions
    >>> f1(3,5)
    Puede que pudiera mejorarse con alguna función de itertools
    '''
    return FSolve( [1 for i in range(e)] , n , onlyFirst=False)

def Delta(laux):
    l1=[laux[i]-laux[i-1] for i in range(1,len(laux))]
    l1=list(set(l1))
    l1.sort()
    return l1
