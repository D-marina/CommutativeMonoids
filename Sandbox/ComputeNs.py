# Esta función calcula Ns usando la suma de los qi.
# Input: Generadores del semigrupo
# Output: Ns

import numpy as np
from numpy import *
import sympy
from sympy import symbols
import integerSmithNormalFormAndApplications
from integerSmithNormalFormAndApplications import *
import math
from math import *

def CalcularNs(a):
    # Calculamos la dimension
    dimension = len(a)
    # Creamos los simbolos ei
    e = symbols('e0:%d'%dimension)
    # Creamos el simbolo s
    s = symbols('s')
    # Calculamos d
    m = equationsToGeneratorsHomogeneusCase(Matrix([a]))
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
    return ceil(max(sol1,sol2)[0])