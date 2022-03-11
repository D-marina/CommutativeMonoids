#!/usr/bin/env python
# coding: utf-8

import PyNormaliz
from PyNormaliz import *
import numpy
from numpy import *
import itertools
from scipy.spatial import ConvexHull


import sys
sys.path.insert(0,'../Class/')
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
    coef = 0
    for i in range(len(x)):
        if x[i] != 0 and r[i] != 0:
            coef = x[i]/r[i]
            break
    if coef == 0:
        return 0
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
#   - hp: hiperplanos soportes.
# OUTPUT:
#   - ecuaciones del rayo.
def EqRay(ray,hp):
    eq = []
    for x in hp:
        if ProdEsc(ray,x) == 0:
            eq.append(x)
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
    aux = AffineTerm(eqray,smg)
    if afin in aux:
        return True
    else:
        return False


# Esta función dice si un conjunto de elementos genera o no un C-semigrupo.
# INPUT:
#   - smg: Sistema generador.
# OUTPUT:
#   - True/False.
def IsCsemigroup(smg):
    # En primer lugar calculamos los rayos del cono con Pynormaliz.
    cono = Cone(cone=smg)
    rayos = cono.ExtremeRays()
    # Calculamos también los hiperplanos soportes.
    hp = cono.SupportHyperplanes()
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
    for ray in rayos:
        # Calculamos las ecuaciones del rayo.
        eqrayo = EqRay(ray,hp)
        # Calculamos los términos afines de los puntos enteros del diamantes.
        afinesDiamante = AffineTerm(eqrayo,diamanteEntero)
        # Nos quedamos con los generadores afines.
        genAfines = AffineSemigroup(afinesDiamante, "generators").getMSG()
        # Comprobamos que por cada afín pasa un generador
        for afin in genAfines:
            if not ExistGenerator(eqrayo, afin,smg):
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


# Calculamos la multiplicidad en cada rayo.
# INPUT:
#   - ray: rayo del cono.
#   - smg: sistema minimal de generadores.
# OUTPUT:
#   - Elemento minimal en el rayo.
def MultiplicityAxis(ray,smg):
    multiplicity = []
    for x in smg:
        aux = BelongAxis(x,ray)
        if aux != 0:
            multiplicity.append([aux,x])
    return sorted(multiplicity)[0][1]


# Calculamos el diamante entero de las multiplicidades.
# INPUT:
#   - cone: rayos del cono.
#   - smg: sistema minimal de generadores.
# OUTPUT:
#   - diamante entero de las multiplicidades
def DiamondMultiplicity(cone,smg):
    extremos = []
    for ray in cone:
        extremos.append(MultiplicityAxis(ray,smg))
    diamante = Diamond(extremos)
    hull = ConvexHull(diamante)
    eq = [list(x) for x in hull.equations]
    cotas = Cube(diamante)
    return IntegerDiamond(eq,cotas)


# Borramos elementos repetidos.
# INPUT:
#   - v: vector de entrada.
# OUTPUT:
#   - vector v sin elementos de entrada.
def DeleteDuplicates(v):
    w = []
    for i in range(len(v)):
        if v[i] not in w:
            w.append(v[i])
    return w


# Creamos una función que calcula la uno norma.
# INPUT:
#   - v: un vector.
# OUTPUT:
#   - Norma-1 del vector v.
def NormOne(v):
    suma = 0
    for x in v:
        suma = suma + x
    return suma


# Calculamos el punto de menor norma del diamante para una recta afín.
# INPUT:
#   - diamond: diamante entero.
#   - afin: término afin de la recta.
#   - eqray: ecuaciones de la recta.
# OUTPUT:
#   - Punto de menor norma del diamante que pasa por la recta afín.
def MinimumPointAffineDiamond(diamond,afin,eqray):
    numeq = len(eqray)
    aux2 = []
    for d in diamond:
        aux = MultiplyMatrix(eqray,[[x] for x in d])
        aux = list(array(aux).flatten())
        if aux == afin:
            aux2.append([NormOne(d),d])
    return(sorted(aux2)[0][1])

# Calculamos los elementos de norma mínima por cada recta afín.
# INPUT:
#   - vdir: vector director de la recta.
#   - seed: punto para empezar a mirar la recta.
#   - smgen: sistema minimal de generadores.
# OUTPUT:
#   - Norma mínima de los elementos del semigrupo de la recta. 
def ComputeMinimumNormInSemigroupLine(vdir,seed,smgen):
    v = seed
    afseg = AffineSemigroup(smgen, "generators")
    while True:
        if afseg.belongs(v):
            return NormOne(v)
        v = [v[i] + vdir[i] for i in range(len(vdir))]
        
        
# Calcula un valor en el rayo y en el semigrupo con norma mayor que una dada.
# INPUT:
#   - n: cota de la norma.
#   - ray: rayo del cono.
#   - gen: generadores del semigrupo.
# OUTPUT:
#   - Valor en el rayo que está en el semigrupo y tiene norma mayor que n.
def ComputeBoundRay(n,ray,gen):
    afseg = AffineSemigroup(gen, "generators")
    v = ray
    i = 1
    while True:
        v = [i*x for x in ray]
        if NormOne(v) > n and afseg.belongs(v):
            return v
        i = i+1
        
def ComputeGaps(gen):
    # Calculamos el cono.
    C = Cone(cone=gen)
    # Calculamos los rayos del cono.
    rayos = C.ExtremeRays()
    # Calculamso los hiperplanos soportes.
    hp = C.SupportHyperplanes()
    # Calculamos del diamante de multiplicidades.
    diamanteM = DiamondMultiplicity(rayos,gen)
    x = []
    # Para cada rayo del cono.
    for ray in rayos:
        # Calculamos el conductor.
        conductor = ConductorAxis(gen,ray)
        # Calculamos las ecuaciones del rayo.
        eq = EqRay(ray,hp)
        # Calculamos los términos afines.
        afinesDiamante = DeleteDuplicates(AffineTerm(eq,diamanteM))
        # Calculamos la máxima norma mínima en cada recta afin.
        maxMinNorm = 0
        for afin in afinesDiamante:
            # En primer lugar miramos la norma minima en el diamante.
            minimumNormDiamond = MinimumPointAffineDiamond(diamanteM,afin,eq)
            # Ahora buscamos el elemento mínimo en el diamante en dicha recta.
            minimumNormSG = ComputeMinimumNormInSemigroupLine(ray,minimumNormDiamond,gen)
            if maxMinNorm < minimumNormSG:
                maxMinNorm = minimumNormSG
        n = NormOne(conductor) + maxMinNorm
        x.append(ComputeBoundRay(n,ray,gen))
    # Calculamos el diamante de esos valores X
    diamanteX = Diamond(x)
    # Creamos un cubo que contenga al diamante.
    cotaX = Cube(diamanteX)
    # Calculamos las ecuaciones del diamante.
    hullX = ConvexHull(diamanteX)
    eqDiamanteX = [list(x) for x in hullX.equations]
    # Calculamos el diamante entero X.
    diamanteEnteroX = IntegerDiamond(eqDiamanteX,cotaX)
    # Inicialmente el conjunto de huecos es vacío.
    gapset = []
    # Definimos el semigrupo afín.
    afseg = AffineSemigroup(gen, "generators")
    for x in diamanteEnteroX:
        if not afseg.belongs(x):
            gapset.append(x)
    return gapset

# Esta función elimina los elentos de un vector cuyo doble está en él.
# INPUT:
#   - v: vector.
# OUTPUT:
#   - vector v sin sus mitades.
def DeleteHalves(v):
    aux = []
    for x in v:
        if not [2*x[i] for i in range(len(x))] in v:
            aux.append(x)
    return aux