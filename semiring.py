#this file contains the semiring and matsemiring classes as well as the relevant examples
import functools
import math
import numpy as np
import sympy as simp

class Semiring:
    def __init__(self,zro,one,pls,tms):
        self.zero=zro
        self.plus=pls
        self.times=tms
        self.one=one
    
    def exp(self,x,n):
        result=self.one
        for i in range(n):
           result=self.times(x,result)
        return result

#Semiring [0,inf]
def plus(a,b):
    return a + b
def times(a,b):
    return a*b
minplus=Semiring(math.inf,0,min,plus)

#lists
def flatten(l):
    flat=[]
    for a in l:
        if type(a)==list:
            for b in a:
                flat.append(b)
        else:
            flat.append(a)
    return flat
    
def listtimes(l,m):
    return [ flatten([a,b]) for a in l for b in m]
lists=Semiring([0],['*'],plus,listtimes)
#Viberti Semiring
viberti=Semiring(0,max,1,times)

# Booleans
def orfunc(x,y):
    return x or y
def andfunc(x,y):
    return x and y
B=Semiring(False,True,orfunc,andfunc)

#NFAsemiring
def setplus(A,B):
    return A | B
def setprod(A,B):
    prod=set()
    for x in A:
        for y in B:
            prod.add(x + y)
    return prod

kleisli=Semiring(set(),set(['']),setplus,setprod)
#regex semiring INCOMPLETE
#def pipe(r,s):
#    return (r+'|' + s)
#def concat(r,s):
#    return r+s
#kleen=Semiring(pipe,concat)
#matrix semirings
#matrix multiplication parameterized by semiring
def matprod(m1,m2,s):
    assert m1.shape[1]==m2.shape[0]
    product= np.full((m1.shape[0],m2.shape[1]),s.zero)
    for i in range(m1.shape[0]):
        for j in range(m2.shape[1]):
            for k in range(m1.shape[1]):
                product[i,j]=s.plus(product[i,j], s.times(m1[i,k],m2[k,j]))  
    return product
    
#constructs semiring object for n by n matrices valued in a semiring s
class MatSemiring(Semiring):
    def __init__(self,s):
        self.semi=s
        self.plus=lambda m1, m2 : np.vectorize(s.plus)(m1,m2)
        self.zero=lambda n: np.full((n,n),s.zero)
        self.times=lambda m1, m2 : matprod(m1,m2,s)
    
    def one(self,n):
        res=np.full((n,n),self.semi.zero)
        for i in range(n):
            res[i,i]=self.semi.one
        return res


def evaluateexpr(expr,s,values):
    if expr.args==():
        if expr.func==simp.core.numbers.One:
            return s.semi.one
        if expr.func==simp.Symbol:
            return values[expr]
    else:
        if expr.func==simp.Add:
            evaluateda=[evaluateexpr(arg,s,values) for arg in expr.args]
            return functools.reduce((lambda a, b:s.plus(a,b)), evaluateda) 
        if expr.func==simp.Mul:
            evaluatedt=[evaluateexpr(arg,s,values) for arg in expr.args]
            return functools.reduce((lambda a, b:s.times(a,b)), evaluatedt)  
        if expr.func==simp.Pow:
            return s.exp(evaluateexpr(expr.args[0],s,values),expr.args[1])
    
listmats=MatSemiring(lists)

#x = np.array([[2,1,3],[4,.3,99],[1,1,2]])
#print(a.fw(x))

#a=MatSemiring(kleisli,3)
#e=set([''])
#x = np.array([[e,set('a'),set()],[set(),e,set('bc')],[set('d'),set(),e]])
#print(a.fw(x))
#it works! 4-26-22


#a=[[1],[2],[3]]
#b=[['a'],['b'],['c']]
#c=np.array([a,b,a])
#print(matprod(c,c,lists))
