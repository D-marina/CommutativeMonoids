#!/usr/bin/env python
# coding: utf-8

import PyNormaliz
from PyNormaliz import *
import numpy
from numpy import *
import itertools
from scipy.spatial import ConvexHull


import sys
sys.path.insert(0, '../Class/')
sys.path.insert(0,'../ClassAffine')
import integerSmithNormalFormAndApplications
from integerSmithNormalFormAndApplications import *
import AffineSemigroup
from AffineSemigroup import *
import auxiliars
from auxiliars import *


# Check if all the coordinates are positives.
def BelongQPositive(v):
    for i in range(len(v)):
        if v[i] < 0:
            return False
    return True


# Pertenece a un eje y si es así, devuelve el múltiplo.
# INPUT:
#   - x: Value for checking if it is in the ray.
#   - r: Minimal value in the ray.
# OUTPUT:
#   - 0: If not belongs to the ray.
#   - A value if belongs to the ray.
def BelongAxis(x,r):
    coef = x[0]/r[0]
    aux2 = [j/coef for j in x]
    if aux2 == r:
        if(int(coef)==coef):
            return int(coef)
        return coef
    else:
        return 0


# Vemos para un eje si se forma un semigrupo.
# INPUT:
#   - gen: Set of generators.
#   - r: Minimal value in the ray.
# OUTPUT:
#   - True/False.
def AxisIsSemigroup(gen,r):
    aux = []
    for x in gen:
        aux.append(BelongAxis(x,r))
    aux2 = [x for x in aux if x != 0]
    if(gcdL(aux2) == 1):
        return True
    else:
        return False


# Vemos si los ejes forman un semigrupo.
# INPUT:
#   - gen: Set of generators.
#   - setR: Set of rays.
# OUTPUT:
#   - True/False.
def AxisAreSemigroup(gen,setR):
    aux = []
    for x in setR:
        aux.append(AxisIsSemigroup(gen,x))
    return all(aux)


# En primer lugar tenemos que calcular el diamante.
# INPUT:
#   - a: Minimal elements of the ray.
# OUTPUT:
#   - Points of the diamond
def Diamond(a):
    aux = list(a)
    aux.append([0 for x in range(len(a))])
    for i in range(len(a)):
        for j in range(i+1,len(a)):
            aux.append(array(a[i])+array(a[j]))
    return [list(x) for x in  aux]


# Ahora veremos si un punto pertenece al diamante o no.
# INPUT:
#   - pt: Punto para comprobar la pertenencia.
#   - eq: Ecuaciones del diamante.
# OUTPUT:
#   - True/False.
def PointBelongsDiamond(pt,eq):
    dim = len(pt)
    for x in eq:
        sum = 0
        for i in range(dim):
            sum = sum + pt[i]*x[i]
        sum = round(sum+x[-1],10)
        if sum > 0:
            return False
    return True


# Calculo el menor cuboide que contenga al diamante.
# INPUT:
#   - d: Vértices del diamante.
# OUTPUT:
#   - Máximos en cada coordenada del cubo.
def Cube(d):
    dim = len(d[0])
    aux = []
    for i in range(dim):
        aux.append(sorted(d, key = lambda x: x[i])[-1][i])
    return(aux)


# Calculamos todos los puntos del diamante.
# INPUT:
#   - eq: Ecuaciones que definen el diamante.
# OUTPUT:
#   - Puntos enteros del diamante.
def IntegerDiamond(eq,cube):
    d = []
    it = itertools.product(*[range(i+1) for i in cube])
    for x in it:
        if PointBelongsDiamond(list(x),eq):
            d.append(list(x))
    return d


def ProdEsc(v1,v2):
    n = len(v1)
    suma = 0
    for i in range(n):
        suma = suma + v1[i] * v2[i]
    return suma


# En primer lugar calculamos las ecuaciones que definen un rayo.
# INPUT:
#   - ray: rayo del cono.
#   - x: valor interno del cono para calibrar las ecuaciones.
# OUTPUT:
#   - ecuaciones del rayo.
def EqRay(ray,x):
    n = len(ray)
    eq = []
    for j in range(1,n):
        aux = []
        aux = aux + [-ray[j]]                 # Valor -aj
        aux = aux + [0 for k in range(j-1)]   # Ceros intermediso
        aux = aux + [ray[0]]                  # Vaor ai
        aux = aux + [0 for k in range(n-j-1)] # Ceros finales
        eq.append(aux)
    # Calibramos las ecuaciones
    for i in range(len(eq)):
        calibre = ProdEsc(eq[i],x)
        if calibre < 0:
            for j in range(len(eq[i])):
                eq[i][j] = -eq[i][j]
    return eq


# Eliminamos los elementos que son todos ceros.
# INPUT:
#   - m: matriz.
# OUTPUT:
#   - Matriz m sin filas nulas.
def DeleteRowZero(m):
    aux = []
    for v in m:
        allzero = True
        for i in range(len(v)):
            if v[i] != 0:
                allzero = False
        if not allzero:
            aux.append(v)
    return aux


# Calculamos los valores afines del diamante con respecto a un rayo.
# INPUT:
#   - eq: Ecuaciones del rayo.
#   - d: Diamante entero.
# OUTPUT:
#   - Terminos afines del diamante sin los ceros.
def AffineTerm(eq,d):
    neq = len(eq)
    dim = len(eq[0])
    afin = []
    for x in d:
        aux = []
        for e in eq:
            aux = aux + [ProdEsc(x,e)]
        afin.append(aux)
    return DeleteRowZero(afin)


# Multiplicar matrices.
# INPUT:
#   - X: una matriz.
#   - Y: otra matriz.
# OUTPUT:
#   - Producto de ambas.
def MultiplyMatrix(X,Y):
    result = [[0 for y in Y[0]] for x in X]
    # iterate through rows of X
    for i in range(len(X)):
        # iterate through columns of Y
        for j in range(len(Y[0])):
            # iterate through rows of Y
            for k in range(len(Y)):
                result[i][j] += X[i][k] * Y[k][j]
    return result


# Calculamos si en cada recta afín paralela al rayo hay un generador.
# INPUT:
#   - eqray: ecuaciones de un rayo.
#   - afinset: valor afín de la recta.
#   - smg: generadores del semigrupo.
# OUTPUT:
#   True/False si en esa recta hay un generador.
def ExistGenerator(eqray, afin,smg):
    numeq = len(eqray)
    for gen in smg:
        aux = MultiplyMatrix(eqray,[[x] for x in gen])
        if aux == afin:
            return True
    return False


# Esta función dice si un conjunto de elementos genera o no un C-semigrupo.
# INPUT:
#   - smg: Sistema generador.
# OUTPUT:
#   - True/False.
def IsCsemigroup(smg):
    # En primer lugar calculamos los rayos del cono con Pynormaliz.
    rayos = Cone(cone=smg).ExtremeRays()
    for ray in rayos:
        # En primer lugar comprobamos que el rayo es un semigrupo en sí mismo.
        if not AxisIsSemigroup(smg,ray):
            return False
    # Calculamos el diamante.
    diamante = Diamond(rayos)
    # Calculamos sus ecuaciones gracias a scipy.spatial -> ConvexHull.
    hull = ConvexHull(diamante)
    eqDiamante = [list(x) for x in hull.equations]
    # Calculamos una cota para los puntos del diamante.
    cotaDiamante = Cube(diamante)
    # Calculamos el diamante entero.
    diamanteEntero = IntegerDiamond(eqDiamante,cotaDiamante)
    # Veamos ahora que por cada rayo las paralelas afines "generadoras" cortan a los generadores del semigrupo.
    # En primer lugar elegimos un punto dentro del diamante para calibrarlo.
    dimension = len(rayos[0])
    calibre  = [0 for x in range(dimension)]
    for ray in rayos:
        calibre = [calibre[i]+ray[i] for i in range(dimension)]
    for ray in rayos:
        # Calculamos las ecuaciones del rayo.
        eqrayo = EqRay(ray,calibre)
        # Calculamos los términos afines de los puntos enteros del diamantes.
        afinesDiamante = AffineTerm(eqrayo,diamanteEntero)
        # Nos quedamos con los generadores afines.
        genAfines = AffineSemigroup(afinesDiamante, "generators").getMSG()
        # Comprobamos que por cada afín pasa un generador
        for afin in genAfines:
            if not ExistGenerator(eqrayo, [afin],smg):
                return False
    return True


# Calculamos el conductor.
# INPUT:
#   - gen: Set of generators.
#   - r: Minimal value in the ray.
# OUTPUT:
#   - conductor
def ConductorAxis(gen,r):
    aux = []
    for x in gen:
        aux.append(BelongAxis(x,r))
    aux2 = [x for x in aux if x != 0]
    return [x*(FrobeniusNumber(aux2)+1) for x in r]