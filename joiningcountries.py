from functools import reduce
from typing import Set
import graphblas as gb
import time
import sys

class pSet(set):

    def __init__(self, L:Set[str]):
        super(pSet,self).__init__(L)
        
    def __add__(self,other):
        return pSet(super(pSet, self).union(other))
        
    def __mul__(self,other):
        prod_elt = lambda x, y: x + (";" if x and y else "") + y
        return pSet({prod_elt(x, y) for x in self for y in other})
    
    def one():
        return pSet([''])
        
    def zero():
        return pSet([])

    

class Bimodule():
    def __init__(self, inp: pSet =pSet({}), out: pSet=pSet({}), mat={},s=pSet):
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

x = pSet({'UK', 'P1','P2','CA'})  
G = Bimodule(
  inp=x,
  out=x,
  mat={
    'UK;UK': pSet({'UK-UK'}),
    'UK;P2': pSet({'UK-P2'}),
    'P1;P1': pSet({'P1-P1'}),
    'P2;UK': pSet({'P2-UK'}),
    'P1;P2': pSet({'P1-P2'}),
    'P2;P1': pSet({'P2-P1'}),
    'P2;P2': pSet({'P2-P2'}),
    'P1;CA': pSet({'P1-CA'}),
    'CA;P1': pSet({'CA-P1'}),
    'CA;CA': pSet({'CA-CA'})
    },
  s=pSet
)



class Decomposition():
    def __init__(self,sh=Bimodule(),m={},obmap={}):
        #shape is a pSet-graph, mapping sends each edge to a gb matrix and mappingob sends each vertex to the size of its bag
        self.shape=sh
        self.mapping=m
        self.mappingob=obmap
    
    def __add__(self,other):
        newshape=Bimodule(self.shape.input,self.shape.output,self.shape.mat | other.shape.mat,pSet)
        newmap=self.mapping | other.mapping
        for key in self.mapping:
            if key in other.mapping:
                newmap[key]=self.mapping[key]+other.mapping[key]
        return Decomposition(newshape,newmap,self.mappingob)
        
    def rpaths(self,r,source,target):
        paths=self.shape.exp(r,pSet).mat[source+';'+target]
        result=gb.Matrix(float,self.mappingob[source],self.mappingob[target])
        for i in paths:
            if i=='':
                prods=[identitygb(self.mappingob[source])]
            else:
                edges=i.split(';')
                prods=[self.mapping[edges[0]]]
                for j in edges[1:]:
                    jmat=self.mapping[j]
                    nextprod=gb.Matrix(float,prods[-1].nrows,jmat.ncols)
                    nextprod << gb.semiring.min_plus(prods[-1] @ jmat)
                    prods.append(nextprod)
            result(accum=gb.binary.min) << prods[-1]
        return result


    
    def add(self,other,source,target):
        newmat=dict()
        newmap=dict()
        key=source+'-'+target
        if key in self.mapping and key in other.mapping:
            newmat[source+';'+target]=pSet({key})
            newmap[key]=self.mapping[key]+other.mapping[key]
        if key in self.mapping and not key in other.mapping:
            newmap[key]=self.mapping[key]
            newmat[source+';'+target]=pSet({key})
        if not key in self.mapping and key in other.mapping:
            newmap[key]=other.mapping[key]
            newmat[source+';'+target]=pSet({key})
        return Decomposition(Bimodule(self.shape.input,self.shape.output,newmat,pSet),newmap,self.mappingob)
            
   
    def __mul__(self,other):
        paths=self.shape*other.shape
        newmat=dict()
        newmapping=dict()
        for x in self.shape.input:
            for y in other.shape.output:
                terms = paths.mat[x+';'+y]
                if terms:
                    newmat[x+';'+y]=pSet({x+'-'+y})
                    result=gb.Matrix(float,self.mappingob[x],other.mappingob[y])
                    for term in terms:
                        edges=term.split(';')
                        firstmat=self.mapping[edges[0]]
                        secondmat=other.mapping[edges[1]]
                        prod=gb.Matrix(float,firstmat.nrows,secondmat.ncols)
                        prod << gb.semiring.min_plus(firstmat @ secondmat)
                        result(accum=gb.binary.min) << prod
                    newmapping[x+'-'+y]=result
        return Decomposition(Bimodule(self.shape.input,other.shape.output,newmat,pSet),newmapping,self.mappingob)
    
    def mult(self,other,source,target):
        terms=(self.shape*other.shape).mat[source +';'+target]
        newmat=dict()
        newmapping=dict()
        if terms:
            newmat[x+';'+y]=pSet({x+'-'+y})
            result=gb.Matrix(float,self.mappingob[source],other.mappingob[target])
            for term in terms:
                edges=term.split(';')
                firstmat=self.mapping[edges[0]]
                secondmat=other.mapping[edges[1]]
                prod=gb.Matrix(float,firstmat.nrows,secondmat.ncols)
                prod << gb.semiring.min_plus(firstmat @ secondmat)
                result(accum=gb.binary.min) << prod
            newmapping[source+'-'+target]=result
        return Decomposition(Bimodule(self.shape.input,other.shape.output,newmat,pSet),newmapping,self.mappingob)
        
    def identity(obmap):
        mapping=dict()
        newmat=dict()
        for i in obmap.keys():
            for j in obmap.keys():
                if i==j:
                    newmat[i+';'+j]=pSet({i+'-'+j})
                    mapping[i+'-'+j]=identitygb(obmap[i])
        return Decomposition(Bimodule(obmap.keys(),obmap.keys(),newmat,pSet),mapping,obmap)
        
    
    def totalgraph(self):
        nnodes=sum(self.mappingob.values())
        total=gb.Matrix(float,nnodes,nnodes)
        bags=list(self.shape.input)
        interval=dict()
        runningsumi=0
        for i in bags:
            runningsumj=0
            irange=[runningsumi]
            for j in bags:
                jrange=[runningsumj]
                if i+'-'+j in mapping.keys():
                    total[runningsumi:runningsumi+self.mappingob[i],runningsumj:runningsumj+self.mappingob[j]] << self.mapping[i+'-'+j]
                runningsumj+=self.mappingob[j]
                jrange.append(runningsumj)
                irange.append(runningsumi+self.mappingob[i])
                interval[(i,j)]=[irange,jrange]    
            runningsumi+=self.mappingob[i]
        return [total,interval]
            
def identitygb(n):
    return gb.Matrix.from_coo([range(n)],[range(n)],[0.0 for i in range(n)])                

if __name__=="__main__":

    nports=10

    CA=gb.io.mmread('roadNet-CA.mtx')
    UK=gb.io.mmread('great-britain_osm.mtx')

    CAsize=CA.ncols-nports
    CAint=gb.Matrix(float,CAsize,CAsize)
    CAint << CA[:CAsize,:CAsize]

    A=gb.Matrix(float,CAsize,nports)
    A << CA[:CAsize,CAsize:CAsize+nports]

    B=gb.Matrix(float,nports,CAsize)
    B << CA[CAsize:CAsize+nports,:CAsize]

    P1=gb.Matrix(float,nports,nports)
    P1 << CA[CAsize:CAsize+nports,CAsize:CAsize+nports]

    C=gb.Matrix.from_dense(np.random.rand(nports,nports))
    D=gb.Matrix.from_dense(np.random.rand(nports,nports))

    UKsize=UK.ncols-nports
    P2=gb.Matrix(float,nports,nports)
    P2 << UK[UKsize:,UKsize:]

    E=gb.Matrix(float,nports,UKsize)
    E << UK[UKsize:,:UKsize]

    F=gb.Matrix(float,UKsize,nports)
    F << UK[:UKsize,UKsize:]

    UKint=gb.Matrix(float,UKsize,UKsize)
    UKint << UK[:UKsize,:UKsize]

    mapping={'UK-UK' : UKint, 'UK-P2' : F, 'P2-UK' : E,'P1-P1':P1, 'P1-P2':     C,'P2-P1':D,'P1-CA':B,'CA-P1':A,'P2-P2':P2,'CA-CA':CAint}

    mappingob={'CA':CAsize,'P1':nports,'P2':nports,'UK':UKsize}

    De=Decomposition(G,mapping,mappingob)
    tic=time.perf_counter()
    a=rnaive(De,3,'CA','CA')
    toc=time.perf_counter()
    print(a.to_coo()[2],"in",toc-tic)

    tic=time.perf_counter()
    b=De.rpaths(3,'CA','CA')
    toc=time.perf_counter()
    print(b.to_coo()[2],"in",toc-tic)
# this function finds a matrix whose entries are the shortest paths in r steps
def rshortestpath(M,r):
    n=M.ncols
    exp=identitygb(n)
    for i in range(r):
        exp(accum=gb.binary.min) << gb.semiring.min_plus(M @ exp)
    return exp

def rshortestpathdecomp(D,r):
    exp=Decomposition.identity(D.mappingob)
    for i in range(r):
        exp = Decomposition.identity(D.mappingob) + (D*exp)
    return exp

def rnaive(D,r,source,target):
    tot,interval=D.totalgraph()
    sourceint=interval[(source,target)][0]
    targetint=interval[(source,target)][1]
    v=gb.Matrix(float,D.mappingob[source],tot.ncols)
    v << tot[sourceint[0]:sourceint[1],:]
    w = gb.Matrix(float,tot.nrows,D.mappingob[target])
    w << tot[:,targetint[0]:targetint[1]]
    res = gb.Matrix(float,sourceint[1]-sourceint[0],targetint[1]-targetint[0])
    res << gb.semiring.min_plus(gb.semiring.min_plus(v @ rshortestpath(tot,r-2)) @ w)
    return res
    






       


