import pandas as pd
import numpy as np
import math

metab = ["x1", "x2","x3","x4", "P"]
fluxes = ["v1","v2","v3","v4","v5","v6"]

#Metabolic flux model 1
data = np.array([[-1,-1,-1,-1,0,0],[0,1,1,0,-1,0],[0,0,0,1,0,-1],[0,0,0,0,1,1],[1,0,0,0,0,0]])
V = np.array([0.1,0.1,0.1,0.2,0.2,0.3])

N = pd.DataFrame(data, metab, fluxes)
print(N)
print(V)
#Make N so fluxes with zero are deleted.


#Calculation H(Pij)(k=1)
H = np.sum(-1 * V * np.log2(V))

#Calculation X(Pij)
outf = np.array(data*V)
outf[outf >0] = 0
inf = np.array(data*V)
inf[inf <0] = 0

X = 0
i = 0
for v in V:
    d_out = np.array(np.transpose(data)[i])
    d_out[d_out >0] =0
    d_in = np.array(np.transpose(data)[i])
    d_in[d_in <0] = 0
    b = np.sum(d_out*np.transpose(outf))
    c = np.sum(d_in*np.transpose(inf))
    #print("X =", V[i], "log2 ", V[i], "/", b, "* ",c)
    a = V[i]* math.log2(V[i]/(b*c))
    X = X + (a)
    i = i + 1

#Calculation Psi(Pij)
Psi = H - X

print("H=",H,"\n X = ",X,"\n Psi =",Psi)


#Calculation H(Pi) and H(pj)
#Is this the same as H(pi.) and H(p.j)???


#Make different dummy models, with X(Pij)=0 and Psi(Pij) = 0
#See what this does to growth and resilience
