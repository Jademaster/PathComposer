from functools import reduce
from typing import List
import math
from scipy.sparse import identity
import graphblas as gb
import time

class ProgSet(list):
# Programmer's Sets: they are lists of strings.
    def __init__(self, L: List[str]):
        super(ProgSet, self).__init__(L)
        
    def __add__(self, other):
        return ProgSet(super(ProgSet, self).__add__(other))

    def __mul__(self, other):
        prod_elt = lambda x, y: x + (";" if x and y else "") + y
        return ProgSet([prod_elt(x, y) for x in self for y in other])
    
    def one():
        return ProgSet([''])
        
    def zero():
        return ProgSet([])

class Bimodule():
    def __init__(self, inp: ProgSet, out: ProgSet, mat,s):
      """
      inp and output are lists
      m is a dict relating pairs of inp, out to elements of the enriching category s taken as its underlying semiring s
      """
      self.input = inp
      self.output = out
      self.mat = mat
      self.semiring=s
      
    def __repr__(self):
        return str(self.mat)
   
    def __add__(self, other):
        if len(self.input) != len(other.input) or len(self.output) != len(other.output):
            raise ValueError("Sum of bimodules must have same shape.")
        if len(self.mat)>=len(other.mat):
            indexingmat=self.mat
        else:
            indexingmat=other.mat
        return Bimodule(self.input, self.output, {k: self.mat.get(k,self.semiring.zero()) + other.mat.get(k,self.semiring.zero()) for k in indexingmat},self.semiring)

    def __mul__(self, other):
        if len(self.output) != len(other.input):
            raise ValueError("Bimodule shapes are incompatible.")
        mat = {i+";"+j: reduce(lambda x, y: x + y,[self.mat.get(i+';'+k,self.semiring.zero()) * other.mat.get(k+';'+j,self.semiring.zero()) for k in self.output]) for i in self.input for j in other.output}
        return Bimodule(self.input, other.output, mat,self.semiring)

    def zero(inp, out, semiring):
        return Bimodule(inp, out, {x: semiring.zero() for x in inp * out},semiring)
        
    def one(obj, semiring):
        zeros = Bimodule.zero(obj, obj, semiring)
        for x in obj:
            zeros.mat[x+';'+x] = semiring.one()
        return zeros
            
    def exp(self, n, semiring):
        if n == 0:
            return Bimodule.one(self.input, semiring)
        else:
            return Bimodule.one(self.input, semiring) + self * self.exp(n - 1, semiring)

x = ProgSet(['g', 'h','k'])
m = Bimodule(
  inp=x,
  out=x,
  mat={
    "g;g": ProgSet(['g']),
    "g;k": ProgSet(['a']),
    "k;g": ProgSet(['b']),
    "k;k": ProgSet(['k']),
    "k;h": ProgSet(['c']),
    "h;k": ProgSet(['d']),
    "h;h": ProgSet(['h'])
  },
  s=ProgSet
)   

bdarysize=3
pathgraph=m.exp(bdarysize,ProgSet)
verts=m.input

roads=gb.io.mmread('GAP-road.mtx')
nodenum=roads.ncols
nonbdary=nodenum-bdarysize
size1=int(nonbdary/2)
size2=nonbdary - size1

road1=gb.Matrix(float,size1,size1)
road2=gb.Matrix(float,size2,size2)
intersection=gb.Matrix(float,bdarysize,bdarysize)
a=gb.Matrix(float,size1,bdarysize)
b=gb.Matrix(float,bdarysize,size1)
c=gb.Matrix(float,bdarysize,size2)
d=gb.Matrix(float,size2,bdarysize)

road1 << roads[:size1,:size1]
road2 << roads[size1 + bdarysize:, size1 + bdarysize:]
intersection << roads[size1 : size1 + bdarysize,size1 : size1 + bdarysize]
a << roads[:size1,size1: size1+bdarysize]
b << roads[size1:size1+bdarysize,:size1]
c << roads[size1: size1+ bdarysize,size1+bdarysize:]
d << roads[size1+bdarysize:,size1:size1+bdarysize]

def identity(source,target):
    a=gb.Matrix(float,source,target)
    a[
mapping={'g' : road1, 'h' : road2, 'k' : intersection, 'a':a,'b':b,'c':c,'d':d}
source_comp='g'
target_comp='h'

pths=pathgraph.mat[source_comp + ';' + target_comp]
print(pths)

#dictionary relating symbols to their sparse matrices
def evaluate(matdict,expr):
    sourcesize=matdict[expr[0][0]].nrows
    targetsize=matdict[expr[0][-1]].ncols
    result=gb.Matrix(float,sourcesize,targetsize)
    for i in expr:
        terms=i.split(';')
        prods=[matdict[terms[0]]]
        for j in terms[1:]:
            jmat=matdict[j]
            nextprod=gb.Matrix(float,prods[-1].nrows,jmat.ncols)
            nextprod << gb.semiring.min_plus(prods[-1] @ matdict[j])
            prods.append(nextprod)
        result(accum=gb.binary.min) << prods[-1]
    return result

tic=time.perf_counter()
answer=evaluate(mapping,pths)
toc=time.perf_counter()
print(toc-tic)


       



