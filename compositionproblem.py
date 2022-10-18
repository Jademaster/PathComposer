#defines classes of graph and composition problem, as well as helper functions

import numpy as np
import networkx as nx
import compositionpatterns as symbols
import copy
import semiring as se

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
    
#change this to arbitrary structured decompositions
class Compositionproblem():

    def __init__(self,g1,g2,intersection,matsemi):
        self.lgraph=g1
        self.rgraph=g2
        self.intersection=intersection
        self.symbols=None
        self.matsemi=matsemi
        self.values=None
        self.joined=None
        self.pushforwards=None
        self.lengths=None
    
    def precompilesymbols(self):
        maxboundary=len(self.intersection)
        self.symbols=symbols.paths(maxboundary)[-1]
    
    def precompilematrices(self):
        Fxmat=nx.floyd_warshall_numpy(nx.from_numpy_matrix(self.lgraph.getmatrix()))
        Fymat=nx.floyd_warshall_numpy(nx.from_numpy_matrix(self.rgraph.getmatrix()))
        self.lgraph=graphfrommatrix(Fxmat,self.lgraph.verts)
        self.rgraph=graphfrommatrix(Fymat,self.rgraph.verts)
    
    #warning, returns multivalued matrix, some methods will not work on result
    def join(self):
        #move keys of inter to the end of self.lgraph.verts and move values to the beginning of g2.verts
        for i in self.intersection:
            self.lgraph.verts.append(self.lgraph.verts.pop(self.lgraph.verts.index(i)))
            self.rgraph.verts.insert(0,self.rgraph.verts.pop(self.rgraph.verts.index(self.intersection[i])))
        #pad old incidences
        g1new={i:np.concatenate((self.lgraph.outedges[i],np.full((len(self.rgraph.verts)-len(self.intersection),),self.matsemi.semi.zero))) for i in self.lgraph.verts}
        g2new={i:np.concatenate((np.full((len(self.lgraph.verts)-len(self.intersection),),self.matsemi.semi.zero),self.rgraph.outedges[i])) for i in self.rgraph.verts}
        newoutedges= {**g1new, **g2new}
        for a in self.intersection.keys():
            newoutedges[a]=(newoutedges[a],newoutedges[self.intersection[a]])
            del newoutedges[self.intersection[a]]
        self.joined=graph(list(newoutedges),newoutedges)
        
    def precompile(self):
        self.precompilematrices()
        self.precompilesymbols()
        self.join()
        self.pushforward()
        self.getlengths()
        self.getblocks()
        
    def totalgraph(self):
        t=self.joined
        for a in self.intersection:
            t.outedges[a]=self.matsemi.plus(t.outedges[a][0],t.outedges[a][1])
        return t
        
    def getlengths(self):
        x=len(self.lgraph.verts)
        k=len(self.intersection)
        y=len(self.rgraph.verts)
        self.lengths=(x-k, k, y-k,x+y-k)
        
    def pushforward(self):
        joined=self.joined
        g=copy.deepcopy(joined)
        h=copy.deepcopy(joined)
        for i in self.intersection.keys():
            g.outedges[i]=g.outedges[i][0]
            h.outedges[i]=h.outedges[i][1]
        self.pushforwards=(g.getmatrix(),h.getmatrix())
        
    
    def getblocks(self):
        gstar,hstar = self.pushforwards
        lengthg,lengthk,lengthh,total=self.lengths
        gg=gstar[0:lengthg,0:lengthg]
        gk=gstar[0:lengthg,lengthg: lengthg+lengthk]
        kg=hstar[lengthg:lengthg+lengthk,0:lengthg]
        gkk=gstar[lengthg:lengthg+lengthk,lengthg:lengthg+lengthk]
        hkk=hstar[lengthg:lengthg+lengthk,lengthg:lengthg+lengthk]
        hh=hstar[total-lengthh:total,total-lengthh:total]
        hk=hstar[total-lengthh:total,lengthg:lengthg+lengthk]
        kh=hstar[lengthg:lengthg + lengthk,total-lengthh:total]
        self.values={symbols.GG:gg,symbols.GK:gk,symbols.KG:kg,symbols.gKK:gkk,symbols.hKK:hkk,symbols.HH:hh,symbols.HK:hk,symbols.KH:kh}  
        
    def shortestpath(self,s,t):
        lengthg,lengthk,lengthh,total=self.lengths
        gstar,hstar = self.pushforwards
        values=self.values
        #find compositional shortest path
        #find s and t
        i=list(self.joined.outedges).index(s)
        j=list(self.joined.outedges).index(t)
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
        symbols.I:np.array([self.matsemi.one(total)[i,j]])
        }
        values.update(newvalues)
        #find the right quadrant and evaluate appropriate expression there
        if 0<=i<lengthg and 0<=j<lengthg:
            expression=self.symbols[0,0]
            return se.evaluateexpr(expression,self.matsemi,values)
        if 0<=i<lengthg and lengthg <= j < lengthg + lengthk:
            expression=self.symbols[0,1]
            return se.evaluateexpr(expression,self.matsemi,values)
        if 0<=i<lengthg and total-lengthh<=j<total:
            expression=self.symbols[0,2]
            return se.evaluateexpr(expression,self.matsemi,values)
        if lengthg<=i<lengthg+lengthk and 0<=j<lengthg:
            expression=self.symbols[1,0]
            return se.evaluateexpr(expression,self.matsemi,values)
        if lengthg<=i<lengthg + lengthk and lengthg <= j < lengthg + lengthk:
            expression=self.symbols[1,1]
            return se.evaluateexpr(expression,self.matsemi,values)
        if lengthg<=i<lengthg+lengthk and total-lengthh<=j<total:
            expression=self.symbols[1,2]
            return se.evaluateexpr(expression,self.matsemi,values)
        if lengthg+lengthk<=i<total and 0<=j<lengthg:
            expression=self.symbols[2,0]
            return se.evaluateexpr(expression,self.matsemi,values)
        if lengthg+ lengthk<=i<total and lengthg <= j < lengthg + lengthk:
            expression=self.symbols[2,1]
            return se.evaluateexpr(expression,self.matsemi,values)
        if lengthg+lengthk<=i<total and total-lengthh<=j<total:
            expression=self.symbols[2,2]
            return se.evaluateexpr(expression,self.matsemi,values)
            


def randomcompproblem(graphsize,boundarysize,matsemi):
    x=np.random.randint(1,10,(graphsize,graphsize))
    y=np.random.randint(1,10,(graphsize,graphsize))
    gx=graphfrommatrix(x,[i for i in range(graphsize)])
    gy=graphfrommatrix(y,[i+graphsize for i in range(graphsize)])
    i=dict(zip([i for i in range(graphsize-boundarysize,graphsize)],[j for j in range(graphsize,graphsize+boundarysize)]))
    return Compositionproblem(gx,gy,i,matsemi)


