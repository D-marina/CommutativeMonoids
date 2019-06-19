#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
from numpy import array
from PyNormaliz import *


# In[2]:


C = Cone(cone = [[6,0],[2,3],[0,7]])
C.HilbertBasis()


# In[3]:


class AffineSemigroupError(Exception):
    pass

class AffineSemigroup:
    '''
    Class for representing a numerical semigroup:
    as=AffineSemigroup("generators", [[5,1,3,7], [....], ....])
    as=AffineSemigroup("equations", \
                        [[[....],[c0]], [[....],[c1]], ...) A X = 0 % C
    as=AffineSemigroup("inequations", [[....], [....], ...]) A X >= 0
    '''
    def __init__(self, input_data, input_type="equations"):
        if not(isinstance(input_data, list)                and all(isinstance(elem, list) for elem in input_data)):
            raise AffineSemigroupError("input_data should be a list of lists")
        if input_type == "equations":
            """
               [[[a11, a21, ..., an1], [i1]],
                          ...               ,
                [[a1m, a2m, ..., anm], [im]]]
            """
            if len(input_data) >= 1 and len(input_data[0]) > 1:
                self.__dim = len(input_data[0][0])
            else:
                raise AffineSemigroupError("incorrect list of equations")
            if any(len(eq[0])!=self.__dim                    or any(not(isinstance(coef, int)) for coef in eq[0])                    or not(isinstance(eq[1], int)) for eq in input_data):
                raise AffineSemigroupError("incorrect input equations")               
            self.equations = input_data
            self.inequations = None
            self.generators = None
        else:
            if len(input_data) >= 1:
                self.__dim = len(input_data[0])
            else:
                raise AffineSemigroupError("a list of elements is expected")
            if input_type == "inequations":
                if any(len(eq)!=self.__dim                       or any(not(isinstance(coef, int)) for coef in eq)                      for eq in input_data):
                    raise AffineSemigroupError("incorrect input inequations")               
                self.equations = None
                self.inequations = input_data
                self.generators = None
            elif input_type == "generators":
                if any(len(eq)!=self.__dim                       or any(not(isinstance(coef, int)) for coef in eq)                      for eq in input_data):
                    raise AffineSemigroupError("incorrect input generators")               
                self.equations = None
                self.inequations = None
                self.generators = input_data
        self.has_min_generators = False

    def __str__(self):
        return('{ '+str(self.getMSG()).replace('[','<').               replace(']','>')+' }')
               
    def belongs(self, x):
        '''
        Check if the element x belong to the affine semigroup:
        >>> whatever.belongs([2,111,17])
        '''
        if all(v==0 for v in x):
            return True
        if any(v<0 for v in x):
            return False
        if ((self.has_min_generators or self.generators!=None)                and any([e==x for e in self.generators])):
            return True
        npax = np.array(x)
        if self.equations != None                 and all(np.array(e[0])*npa % e[1] == 0                         for e in self.equations):
            return True
        else:
            return False
        if self.inequations != None                 and all(np.array(e[0])*npa >= 0                         for e in self.inequations):
            return True
        else:
            return False

        if self.has_min_generators or self.generators != None:
            return belongsByGens(np.array(x), self.generators)

    def _computeMSG(self):
        '''
        Computes a minimal system of generator of the affine semigroup.
        >>> whatever.getMSG()
        '''
        if self.has_min_generators:
            return self.generators
        if self.equations != None:
            cngs = [eq[0]+[eq[1]] for eq in self.equations if eq[1]!=0]
            if len(cngs)==0:
                cono=Cone(equations=[eq[0] for eq in self.equations])
            else:
                eqns = [eq[0] for eq in self.equations if eq[1]==0]
                cono=Cone(congruences=cngs, equations=eqns)
                
            self.generators=cono.HilbertBasis(DualMode=True)
            if self.generators==[]:
                raise AffineSemigroupError("unable to compute Hilbert basis")
        elif self.inequations != None:
            cono=Cone(inequalities=self.inequations,                       signs=[1 for s in range(self.__dim)])
            self.generators=cono.HilbertBasis(DualMode=True)
            if self.generators==[]:
                raise AffineSemigroupError("unable to compute Hilbert basis")
        else: # "generators" case
            setgen=set([tuple(e) for e in self.generators])
            self.generators=[list(e) for e in setgen                    if all([not belongsByGens(list(e), list(setgen - {e}))])]
        self.has_min_generators=True
    
    def getMSG(self):
        if not self.has_min_generators:
            self._computeMSG()
        return self.generators
    
    def getExpressions(self,x):
        '''
        Computes the factorization(s) of x. 
        >>> whatever.getExpressions([9,1,74,2])
        '''
        assert len(x)==self.__dim,                 "parameter length differs from affine semigroup's dimension"
        if not self.has_min_generators:
            self._computeMSG()
        iheq=[list(e) for e in zip(*self.generators)]
        nterms=len(iheq[0])
        for i in range(len(x)):
            iheq[i].append(-x[i])
        cono = Cone(inhom_equations=iheq)
        facts = cono.ModuleGenerators(DualMode=True)
        if facts == []:
            return []
        else:
            return [e[0:nterms] for e in facts]
    
def belongsByGens(x, gens):
    '''
    x is a numpy.array
    gens is a list of list whose size is equal to x length
    Returns True in case x = a_1*gens[0]+...+a_·*gens[·]
    '''
    if gens==[]:
        return False
    if all(e==0 for e in x):
        return True
    if any(e<0 for e in x):
        return False
    if belongsByGens(x-np.array(gens[0]), gens):
        return True
    elif len(gens)==1: 
        return False
    else: 
        return belongsByGens(x, gens[-(len(gens)-1):])


# In[10]:


afs=AffineSemigroup([[2,0],[0,4],[4,4]], "generators")
print(afs.getMSG())


# In[4]:


afs=AffineSemigroup([[[-2,1],3],[[9,1],8]], "equations")
print(afs.getMSG())
print(afs)
afs.getExpressions([445, 11])


# In[5]:


afs=AffineSemigroup([[[4,-2,1],3],[[2,9,1],5]], "equations")
print(afs.getMSG())
afs.getExpressions([1, 26, 24])


# In[ ]:


afs=AffineSemigroup([[-2,1,1,0],[-10,21,3,9]], "inequations")
print(afs.getMSG())
afs.getExpressions([10, 1, 6, 4])
#afs.getExpressions([8, 1, 16, 5]) # Puede tardar...
#afs.getExpressions([38, 5, 41, 28])


# In[ ]:


afs=AffineSemigroup([[-23,11,9,1],[-10,21,13,9],[1,12,0,1]], "inequations")
print(afs.getMSG())
#afs.getExpressions([10, 21, 1, 0])
afs.getExpressions([8, 1, 16, 5])


# In[ ]:


afs=AffineSemigroup([[13,190, 9,2],[15,18,14,3],[6,7,2,1],[5,11,1,4],[3,4,10,1]], "generators")
print(afs.getMSG())
afs.getExpressions([28, 208, 23, 5])

