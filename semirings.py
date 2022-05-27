import math
import numpy as np
import compositionpatterns as symbols
import time
import sympy as simp
import copy
import networkx as nx
import random
from functools import reduce

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
        
class number():
    def __init__(self,x,s):
        self.value=x
        self.semiring=s
    
    def __add__(self,y):
        self.semiring.plus(self.value,y.value)
    
    def __mul__(self,y):
        self.semiring.times(self.value,y.value)    
        
#Semiring [0,inf]
def plus(a,b):
    return a + b
def times(a,b):
    return a*b
minplus=Semiring(math.inf,0,min,plus)

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
              
    def fw(self,g):
        dim=g.shape[0]
        g=self.plus(g,self.one(dim))
        # stepping loop
        for k in range(dim):
            # outer loop
            for i in range(dim):
                # inner loop
                for j in range(dim):
                    # replace direct path with path through k if direct path is longer
                    g[i,j] = self.semi.plus(g[i,j], self.semi.times(g[i,k], g[k,j]))
        return g
    
    def dykestra(self,g,s,t):
        univistednodes=g.verts
        shortest_path = {}
        previous_nodes = {} 
        for node in unvisited_nodes:
            shortest_path[node] = self.semi.zero 
        shortest_path[start_node] = self.semi.one
        while unvisitednodes:
            current_min_node = None
            for node in unvisitednodes: # Iterate over the nodes
                if current_min_node == None:
                    current_min_node = node
                elif shortest_path[node] < shortest_path[current_min_node]:
                    current_min_node = node
        # The code block below retrieves the current node's neighbors and updates their distances
                    neighbors = graph.get_outgoing_edges(current_min_node)
                    for neighbor in neighbors:
                            tentative_value = shortest_path[current_min_node] + graph.value(current_min_node, neighbor)
                            if tentative_value < shortest_path[neighbor]:
                                shortest_path[neighbor] = tentative_value
                                previous_nodes[neighbor] = current_min_node
            unvisited_nodes.remove(current_min_node)
        return previous_nodes, shortest_path
               
#minplusmat=MatSemiring(minplus,3)
#x = np.array([[2,1,3],[4,.3,99],[1,1,2]])
#print(a.fw(x))

#a=MatSemiring(kleisli,3)
#e=set([''])
#x = np.array([[e,set('a'),set()],[set(),e,set('bc')],[set('d'),set(),e]])
#print(a.fw(x))
#it works! 4-26-22


    
class graph():
#vertexmap is an array sending a set of vertices to the integer indices of the matrix
    def __init__(self,v,o):
        #verts is list of vertices, outedges is a dictionary relating each vertex to a list of row vectors of out weights (they can be multivalued! just for convenience, don't worry too much, it just means that each outvector is wrapped in a tuple
        self.verts=v
        self.outedges=o
        
    def getmatrix(self):
        return np.stack([self.outedges[i] for i in self.verts])
    
    def getinedges(self):
        mat=self.getmatrix()
        return { i: mat[:,self.verts.index(i)] for i in self.verts}
        
    def getweight(self,i,j):
        return self.outedges[i][self.verts.index(j)]
    
def graphfrommatrix(mat,verts):
    #takes in a list of vertex names and a matrix of the same size
    outedges={ i : mat[verts.index(i),:] for i in verts }
    return graph(verts,outedges)

#warning, returns multivalued matrix, some methods will not work on result
def join(g1,g2,inter,s):
#move keys of inter to the end of g1.verts and move values to the beginning of g2.verts
    for i in inter:
        g1.verts.append(g1.verts.pop(g1.verts.index(i)))
        g2.verts.insert(0,g2.verts.pop(g2.verts.index(inter[i])))
    #pad old incidences
    g1new={i:np.concatenate((g1.outedges[i],np.full((len(g2.verts)-len(inter),),s.semi.zero))) for i in g1.verts}
    g2new={i:np.concatenate((np.full((len(g1.verts)-len(inter),),s.semi.zero),g2.outedges[i])) for i in g2.verts}
    newoutedges= {**g1new, **g2new}
    for a in inter.keys():
        newoutedges[a]=(newoutedges[a],newoutedges[inter[a]])
        del newoutedges[inter[a]]
    return graph(list(newoutedges),newoutedges)

def evaluateexpr(expr,s,values):
    if expr.args==():
        if expr.func==simp.core.numbers.One:
            return s.semi.one
        if expr.func==simp.Symbol:
            return values[expr]
    else:
        if expr.func==simp.Add:
            evaluateda=[evaluateexpr(arg,s,values) for arg in expr.args]
            return reduce((lambda a, b:s.plus(a,b)), evaluateda) 
        if expr.func==simp.Mul:
            evaluatedt=[evaluateexpr(arg,s,values) for arg in expr.args]
            return reduce((lambda a, b:s.times(a,b)), evaluatedt)  
        if expr.func==simp.Pow:
            return s.exp(evaluateexpr(expr.args[0],s,values),expr.args[1])

def pushforward(joined,inter):
    g=copy.deepcopy(joined)
    h=copy.deepcopy(joined)
    for i in inter.keys():
        g.outedges[i]=g.outedges[i][0]
        h.outedges[i]=h.outedges[i][1]
    return (g.getmatrix(),h.getmatrix())

def getlengths(x,y,inter,joined):
    return  (len(x.verts)-len(inter), len(inter), len(y.verts)-len(inter),len(joined.verts))
    
def getblocks(gstar,hstar,lengthg,lengthk,lengthh,total):
    #find blocks
    gg=gstar[0:lengthg,0:lengthg]
    gk=gstar[0:lengthg,lengthg: lengthg+lengthk]
    kg=hstar[lengthg:lengthg+lengthk,0:lengthg]
    gkk=gstar[lengthg:lengthg+lengthk,lengthg:lengthg+lengthk]
    hkk=hstar[lengthg:lengthg+lengthk,lengthg:lengthg+lengthk]
    hh=hstar[total-lengthh:total,total-lengthh:total]
    hk=hstar[total-lengthh:total,lengthg:lengthg+lengthk]
    kh=hstar[lengthg:lengthg + lengthk,total-lengthh:total]
    return {symbols.GG:gg,symbols.GK:gk,symbols.KG:kg,symbols.gKK:gkk,symbols.hKK:hkk,symbols.HH:hh,symbols.HK:hk,symbols.KH:kh}  
          
def shortestpaths(joined,gstar,hstar,s,t,symbolmat,semiring,values,lengthg,lengthk,lengthh,total):
    #find s and t
    i=list(joined.outedges).index(s)
    j=list(joined.outedges).index(t)
    valueg=np.array([gstar[i,j]])
    valueh=np.array([hstar[i,j]])
    #add start and end vectors to value list
    newvalues={
    #start symbols
    symbols.GGs:gstar[i,0:lengthg][None,:],
    symbols.GKs:gstar[i,lengthg:lengthg+lengthk][None,:],
    symbols.KGs:gstar[i,0:lengthg][None,:],
    symbols.gKKs:gstar[i,lengthg:lengthg+lengthk][None,:],
    symbols.hKKs:hstar[i,lengthg:lengthg+lengthk][None,:],
    symbols.KHs:hstar[i,total-lengthh:total][None,:],
    symbols.HKs:hstar[i,lengthg:lengthg+lengthk][None,:],
    symbols.HHs:hstar[i,total-lengthh:total][None,:],
    #final symbols
    symbols.GGt:gstar[0:lengthg,j][:,None],
    symbols.GKt:gstar[0:lengthg,j][:,None],
    symbols.KGt:gstar[lengthg:lengthg+lengthk,j][:,None],
    symbols.gKKt:gstar[lengthg:lengthg+lengthk,j][:,None],
    symbols.hKKt:hstar[lengthg:lengthg+lengthk,j][:,None],
    symbols.KHt:hstar[lengthg:lengthg+lengthk,j][:,None],
    symbols.HKt:hstar[total-lengthh:total,j][:,None],
    symbols.HHt:hstar[total-lengthh:total,j][:,None],
    #length1 symbols
    symbols.GGst:valueg,
    symbols.GKst:valueg,
    symbols.KGst:valueg,
    symbols.gKKst:valueg,
    symbols.HHst:valueh,
    symbols.HKst:valueh,
    symbols.KHst:valueh,
    symbols.hKKst:valueh,
    symbols.I:np.array([semiring.one(total)[i,j]])
    }
    values.update(newvalues)
    #find the right quadrant and evaluate appropriate expression there
    if 0<=i<lengthg and 0<=j<lengthg:
        expression=symbolmat[0,0]
        return evaluateexpr(expression,semiring,values)
    if 0<=i<lengthg and lengthg <= j < lengthg + lengthk:
        expression=symbolmat[0,1]
        return evaluateexpr(expression,semiring,values)
    if 0<=i<lengthg and total-lengthh<=j<total:
        expression=symbolmat[0,2]
        return evaluateexpr(expression,semiring,values)
    if lengthg<=i<lengthg+lengthk and 0<=j<lengthg:
        expression=symbolmat[1,0]
        return evaluateexpr(expression,semiring,values)
    if lengthg<=i<lengthg + lengthk and lengthg <= j < lengthg + lengthk:
        expression=symbolmat[1,1]
        return evaluateexpr(expression,semiring,values)
    if lengthg<=i<lengthg+lengthk and total-lengthh<=j<total:
        expression=symbolmat[1,2]
        return evaluateexpr(expression,semiring,values)
    if lengthg+lengthk<=i<total and 0<=j<lengthg:
        expression=symbolmat[2,0]
        return evaluateexpr(expression,semiring,values)
    if lengthg+ lengthk<=i<total and lengthg <= j < lengthg + lengthk:
        expression=symbolmat[2,1]
        return evaluateexpr(expression,semiring,values)
    if lengthg+lengthk<=i<total and total-lengthh<=j<total:
        expression=symbolmat[2,2]
        return evaluateexpr(expression,semiring,values)


#algo! the following code runs benchmarking, make sure to adjust parameters

# initialize semiring
minplusmat=MatSemiring(minplus)
#create graphs
graphsize=200
x=np.random.randint(1,10,(graphsize,graphsize))
y=np.random.randint(1,10,(graphsize,graphsize))
#compute symbols
maxboundary=5
matsymbols=symbols.paths(maxboundary)
#precompilation
tic=time.perf_counter()
Fxmat=nx.floyd_warshall_numpy(nx.from_numpy_matrix(x))
Fymat=nx.floyd_warshall_numpy(nx.from_numpy_matrix(y))
toc=time.perf_counter()
precomptime=toc-tic
print('precompilation time is '+str(precomptime))
Fx=graphfrommatrix(Fxmat,[i for i in range(graphsize)])
Fy=graphfrommatrix(Fymat,[i+graphsize for i in range(graphsize)])


import matplotlib.pyplot as plt
averageingnumber=20
stepsize=1
nxtimes=np.zeros((maxboundary-1,averageingnumber))
comptimes=np.zeros((maxboundary-1,averageingnumber))
ks=np.arange(1,maxboundary,stepsize)

#range of size of intersection
for k in ks:
    #define intersection of right size, join matrices
    intersection=dict(zip([i for i in range(graphsize-k,graphsize)],[j for j in range(graphsize,graphsize+k)]))
    j=join(Fx,Fy,intersection,minplusmat)
    pushout=copy.deepcopy(j)
    for a in intersection:
        pushout.outedges[a]=minplusmat.plus(pushout.outedges[a][0],pushout.outedges[a][1])
    for i in range(averageingnumber):
        #choose s and t randomly
        s=random.choice(j.verts)
        t=random.choice(j.verts)
        #print(s,t)
        tic=time.perf_counter()
        sp=nx.dijkstra_path_length(nx.from_numpy_matrix(pushout.getmatrix()),list(pushout.outedges).index(s),list(pushout.outedges).index(t))
        toc=time.perf_counter()
        nxtimes[k-1,i]=toc-tic
        #print('nx is ' +str(sp)+' with time ' + str(toc-tic))
        tic=time.perf_counter()
        lengthg,lengthk,lengthh,total=getlengths(Fx,Fy,intersection,j)
        gstar,hstar=pushforward(j,intersection)
        values=getblocks(gstar,hstar,lengthg,lengthk,lengthh,total)
        #find compositional shortest path
        res=shortestpaths(j,gstar,hstar,s,t,matsymbols[k],minplusmat,values,lengthg,lengthk,lengthh,total)
        toc=time.perf_counter()
        comptimes[k-1,i]=toc-tic
        #print('comp is ' + str(res) + ' with time ' + str(toc-tic))
    
nxmean = np.mean(nxtimes, axis=1)
nxstd = np.std(nxtimes, axis=1)
compmean = np.mean(comptimes, axis=1)
compstd = np.std(comptimes, axis=1)
print(nxmean,nxstd)
print(compmean,compstd)
#plt.plot(ks,nxmean,label='networkx',color='#55CDFC')
#plt.fill_between(ks, nxmean-nxstd,nxmean + nxstd,alpha=0.5,color='#55CDFC')
#plt.plot(ks,compmean,label='compositional',color='#F7A8B8')
#plt.fill_between(ks, compmean-compstd,compmean + compstd,alpha=0.5,color='#F7A8B8')
#plt.show()





