import Auxiliares as aux

import PyNormaliz
from PyNormaliz import *
import numpy
from numpy import *
import itertools
from scipy.spatial import ConvexHull


import sys
sys.path.insert(0,'./CommutativeMonoids/Class/')
sys.path.insert(0,'./CommutativeMonoids/ClassAffine')
import integerSmithNormalFormAndApplications
from integerSmithNormalFormAndApplications import *
import AffineSemigroup
from AffineSemigroup import *
import auxiliars
from auxiliars import *

def DeleteCSemigruposDuplicates(v):
    limpio = []
    aux = []
    for x in v:
        if x.generadores not in aux:
            aux.append(x.generadores)
            limpio.append(x)
    return limpio

def IsCsemigroup(smg):
    return aux.IsCsemigroup(smg)

def ComputeGaps(generadores):
    return aux.ComputeGaps(generadores)

class Csemigroup:
    gaps = None
    pf = None
    def __init__(self,generadores):
        self.generadores = AffineSemigroup(generadores,input_type="generators").getMSG()
        self.generadores.sort(key=lambda row: row[1:])
        if not aux.IsCsemigroup(self.generadores):
            raise Exception("The set do not form a C-Semigroup")
        self.gaps = aux.ComputeGaps(generadores)
    def GetGenerators(self):
        return self.generadores
    def GetGaps(self):
        return self.gaps
    # This function compute the Pseudo-Frobenius elements of a C-semigroup.
    # INPUT:
    #   - gen: generadores del C-semigrupo.
    #   - gap: conjunto de huecos del C-semigrupo.
    # OUTPUT:
    #   - conjunto de pseudo-frobenius.
    def ComputePseudoFrobenius(self):
        if self.pf != None:
            return self.pf
        cs = AffineSemigroup(self.generadores,input_type="generators")
        pf = []
        for x in self.gaps:
            ispf = True
            for y in self.generadores:
                if not cs.belongs([x[i]+y[i] for i in range(len(x))]):
                    ispf = False
            if ispf:
                pf.append(x)
        self.pf = pf
        return pf
    
    # Esta función clasifica en irreducible o no a un C-semigrupo.
    # INPUT:
    #   - pf: Pseudo-Frobenius.
    # OUTPUT:
    #   True/False
    def IsIrreducible(self):
        if self.pf == None:
            self.ComputePseudoFrobenius()
        if (len(self.pf) == 1) or (len(self.pf) == 2 and len(aux.DeleteHalves(self.pf)) == 1):
            return True
        return False
    
    def DecomposeCSemigroup(self):
        if self.IsIrreducible():
            return [self]
        candidatos = []
        pf = aux.DeleteHalves(self.pf)
        for x in pf:
            cs = AffineSemigroup(self.generadores+[x],input_type="generators")
            candidatos.append(Csemigroup(cs.getMSG()))
            #candidatos.append(cs.getMSG())
        return candidatos
    
    def DecomposeIrreducible(self):
        candidatos = [self]
        irreducibles = []
        control = 0
        while candidatos != [] and control <= 100:
            auxiliar = []
            for x in candidatos:
                auxiliar = auxiliar + x.DecomposeCSemigroup()
            auxiliar = DeleteCSemigruposDuplicates(auxiliar)
            for x in auxiliar:
                if x.IsIrreducible():
                    irreducibles.append(x)
                    auxiliar.remove(x)
            candidatos = auxiliar
            control = control+1
        print([x.generadores for x in irreducibles])
        print()
        for x in irreducibles:
            for y in irreducibles:
                if y.IsSubsemigroup(x) and x != y:
                    irreducibles.remove(y)
        return irreducibles
            
    def __belongsByGens(self,x, gens):
        '''
        x is a numpy.array
        gens is a list of list whose size is equal to x length
        Returns True in case x = a_1*gens[0]+...+a_·*gens[·]
        '''
        if gens==[]:
            return False
        if all([e==0 for e in x]):
            return True
        if any([e<0 for e in x]):
            return False
        if belongsByGens(x-np.array(gens[0]), gens):
            return True
        elif len(gens)==1: 
            return False
        else: 
            return belongsByGens(x, gens[-(len(gens)-1):])
        
    def Belongs(self,x):
        return self.__belongsByGens(x,self.generadores)
    
    # Nos devuelve True si other es subsemigrupo de self
    def IsSubsemigroup(self,otherSemigroup):
        if all([x in otherSemigroup.gaps for x in self.gaps]):
            return True
        else:
            return False