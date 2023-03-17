#!/usr/bin/env python
# coding: utf-8

# In[1]:


import CsemigroupsCpp
from CsemigroupsCpp import *
import PyNormaliz
from PyNormaliz import *
import itertools
from scipy.spatial import ConvexHull


# In[2]:


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


# In[3]:


class Csemigroup:
    def __init__(self,generators):
        self.gaps = None
        self.pf = None
        self.irreducible = None
        aux = self.__ComputeMinimalGenerators(generators)
        self.generators = list([list(x) for x in aux])
        self.generators.sort(key=lambda row: row[1:])
        coneSg = Cone(cone=self.generators)
        self.rays = coneSg.ExtremeRays()
        self.hyperplanes = coneSg.SupportHyperplanes()
        if not self.__IsCSemigroup():
            raise Exception("The set do not form a C-Semigroup")
        
    def __IsCSemigroup(self):
        if not axisAreSemigroup(self.generators,self.rays):
            raise Exception("This set does not generate a numerical semigroup.")
        diamondA = diamond(self.rays)
        hull = ConvexHull(diamondA)
        eqDiamond = [list(x) for x in hull.equations]
        boundDiamond = Cube(diamondA)
        it = itertools.product(*[range(i+1) for i in boundDiamond])
        candidates = [list(x) for x in it]
        integerDiamond = filterPoints(candidates,eqDiamond)
        return studyRays(self.rays,self.hyperplanes,integerDiamond,self.generators);

        
    def __ComputeMinimalGenerators(self,generators):
        return computeMSG(generators)
        
    def GetGenerators(self):
        return self.generators
    def GetRays(self):
        return self.rays
    def GetHyperplanes(self):
        return self.hyperplanes
    def GetGaps(self):
        if self.gaps != None:
            return self.gaps
        else:
            multiplicitiesInAxes = multiplicityAllAxes(self.generators,self.rays)
            diamondMult = diamond(multiplicitiesInAxes)
            hull = ConvexHull(diamondMult)
            eqDiamond = [list(x) for x in hull.equations]
            boundDiamond = Cube(diamondMult)
            it = itertools.product(*[range(i+1) for i in boundDiamond])
            candidates = [list(x) for x in it]
            integerDiamond = filterPoints(candidates,eqDiamond)
            #print(len(integerDiamond))
            #print(integerDiamond)
            # self.gaps = computeGaps(self.generators,self.rays,self.hyperplanes)
            diamondX = computeXDiamond(self.generators, self.rays, self.hyperplanes, integerDiamond)
            #print(diamondX)
            boundDiamondX = Cube(diamondX)
            #print(boundDiamondX)
            hullX = ConvexHull(diamondX)
            eqDiamondX = [list(x) for x in hullX.equations]
            itX = itertools.product(*[range(i+1) for i in boundDiamondX])
            candidatesX = [list(x) for x in itX]
            integerDiamondX = filterPoints(candidatesX,eqDiamondX)
            #print(len(integerDiamondX))
            aux = filterGaps(self.generators, integerDiamondX)
            self.gaps = [list(xx) for xx in aux]
            return self.gaps
    
    def GetPseudoFrobenius(self):
        if self.pf != None:
            return self.pf
        if self.gaps == None:
            self.GetGaps()
        else:
            aux = computePseudoFrobenius(self.generators,self.gaps)
            self.pf = [list(xx) for xx in aux]
            return self.pf

    def IsIrreducible(self):
        if self.irreducible != None:
            return self.irreducible
        if self.pf != None:
            self.GetPseudoFrobenius()
        self.irreducible = IsIrreducible()
        return self.irreducible