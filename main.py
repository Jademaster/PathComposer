import semiring as s
import compositionproblem as c
import time
import networkx as nx
import random
import numpy as np
from functools import reduce


# initialize semiring
minplusmat=s.MatSemiring(s.minplus)       
    
#create graphs
a=c.randomcompproblem(200,5,minplusmat)
a.precompile()
j=a.totalgraph()
s=random.choice(j.verts)
t=random.choice(j.verts)
#nx's attempt
tic=time.perf_counter()
sp=nx.dijkstra_path_length(nx.from_numpy_matrix(j.getmatrix()),list(j.outedges).index(s),list(j.outedges).index(t))
toc=time.perf_counter()
nxtime=toc-tic
#compositional attempt
tic=time.perf_counter()
cp=a.shortestpath(s,t)
toc=time.perf_counter()
ctime=toc-tic
print(sp,cp,"in",nxtime,ctime)



