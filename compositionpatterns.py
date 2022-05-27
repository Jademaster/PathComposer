from sympy import symbols, pprint, eye, simplify, diag
from sympy.matrices import Matrix
import numpy as np
GG, GK, KG, gKK, HH, HK, KH, GGs, GKs, KGs, gKKs, HHs, HKs, KHs, GGt, GKt, KGt, gKKt, HHt, HKt, KHt, hKK, hKKs, hKKt = symbols('GG, GK, KG, gKK, HH, HK, KH, GGs, GKs, KGs, gKKs, HHs, HKs, KHs, GGt, GKt, KGt, gKKt, HHt, HKt, KHt, hKK, hKKs, hKKt', commutative=False)
G = Matrix([[GG, GK, 0],[KG, gKK, 0],[0,0,0]])
Gs = Matrix([[GGs, GKs, 0],[KGs, gKKs, 0],[0,0,0]])
Gt = Matrix([[GGt, GKt, 0],[KGt, gKKt, 0],[0,0,0]])
H = Matrix([[0,0,0],[0, hKK, KH],[0,HK, HH]])
Hs = Matrix([[0,0,0],[0, hKKs, KHs],[0,HKs, HHs]])
Ht = Matrix([[0,0,0],[0, hKKt, KHt],[0,HKt, HHt]])

GGst, GKst, KGst, gKKst, HHst, HKst, KHst, hKKst, I = symbols('GGst, GKst, KGst, gKKst, HHst, HKst, KHst, hKKst, I',commutative=False)
Gst=Matrix([[GGst, GKst, 0],[KGst, gKKst, 0],[0,0,0]])
Hst=Matrix([[0,0,0],[0, hKKst, KHst],[0,HKst, HHst]])
B=[True,False]
One=diag(I,I,I)
def boolvect(n):
    return [B[x%2] for x in range(n)]
    
def alt(n):
    forward=Gs
    backward=Hs
    for i in boolvect(n-1):
        if i:
            forward=forward*H
            backward=backward*G
        else:
            forward=forward*G
            backward=backward*H  
    if n%2==1:
        forward=forward*Ht
        backward=backward*Gt
    else:
        forward=forward*Gt
        backward=backward*Ht
    return (forward + backward)

def altnostart(n):
    forward=G
    backward=H
    for i in boolvect(n):
        if i:
            forward=forward*H
            backward=backward*G
        else:
            forward=forward*G
            backward=backward*H
    return (forward + backward)

def paths(k):
  pathlist=[One,One + Gst + Hst]
  for i in range(k-1):
    pathlist.append(pathlist[-1]+alt(i+1))
  return pathlist
  
def pathsnostart(k):
    paths=eye(3)
    for i in range(k):
        paths=paths+altnostart(i)
    return paths
  
  



  


 

