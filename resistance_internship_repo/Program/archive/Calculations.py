import pandas as pd
import numpy as np

Metab = ["x1", "x2","x3","x4", "P"]
fluxes = ["v1","v2","v3","v4","v5","v6"]

#Metabolic flux model 1
data = np.array([[-1,-1,-1,-1,0,0],[0,1,1,0,-1,0],[0,0,0,1,0,-1],[0,0,0,0,1,1],[1,0,0,0,0,0]])
V = np.array([0.1,0.1,0.1,0.2,0.2,0.3])

#Pij's determination
Pij = V

#Determine P..
P = np.sum(V)

#Determine Pi.
data_outflux = np.array(data)
data_outflux[data_outflux == 1] = 0
data_outflux[data_outflux == -1] = 1
Pi = np.array([])
for row in data_outflux:
    temp = np.array(row * V)
    Pi = np.append(Pi, np.sum(temp))

#Determine P.j
data_influx = np.array(data)
data_influx[data_influx == -1] = 0
Pj = np.array([])
for row in data_influx:
    temp = np.array(row * V)
    Pj = np.append(Pj, np.sum(temp))

#Calculation H(Pij),  and Psi(Pij) (k=1)
H = np.sum(-1 * V * np.log2(V))

#Calculation X(Pij)
#For each row in V, each column in data can tell you what i is and j.
#in each row of data.T, check which metabolite is i (-1) and j (1)
outf = np.array(data*V)
outf[outf >0] = 0
inf = np.array(data*V)
inf[inf <0] = 0

print(data, "is N \n" )
print(V, "is V \n")
print(outf, "is i \n")
print(inf, "is j \n")
print(outf * data, "is i*N \n")
print(inf * data, "is j*N \n")
#print(np.transpose(data)[0])

X = np.sum(1)



#Calculation H(Pi) and H(pj)
#Is this the same as H(pi.) and H(p.j)???


#Make different dummy models, with X(Pij)=0 and Psi(Pij) = 0
#See what this does to growth and resilience





























#Calculation H(Pij)(k=1)
H = np.sum(-1 * V * np.log2(V))

#Calculation X(Pij)
#Make matrixes with only outflow (Pi) and only inflow (Pj)
outf = np.array(N*V)
outf[outf >0] = 0
inf = np.array(N*V)
inf[inf <0] = 0
print(outf, "= outf\n")
print(inf, "= inf \n")
X = 0
i = 0
for v in V:
    d_out = np.array(np.transpose(N)[i])
    d_out[d_out >0] =0
    #d_out = d_out / np.sum(np.absolute(d_out))
    d_in = np.array(np.transpose(N)[i])
    d_in[d_in <0] = 0
    #d_in = d_in / np.sum(np.absolute(d_in))

    print(d_out, "is d_out \n")
    print(d_in, "is d_in \n")
    print(np.transpose(outf))
    print(np.transpose(inf))
    b = np.sum(d_out*np.transpose(outf))
    c = np.sum(d_in*np.transpose(inf))
    print(b,c, "= b (i) and c (j)\n")
    print("X =", V[i], "log2 ", V[i], "/", b, "* ",c)
    z = V[i]* math.log2(V[i]/(b*c))
    X = X + (z)
    i = i + 1

#Calculation Psi(Pij)
psi = H - X

#Calculating flux that remains internal and external
#P3 - sum of P1,2,4,5,6,7
#T_intern = flux[2] - np.sum(flux[0:2]) - np.sum(flux[6:8]) - np.sum(flux[11:14])
#T_extern = flux[2] - T_intern
#print("T_intern=",T_intern,"\nT_extern=",  T_extern)

#Calculating A, C and \Psi
A = 1#T_intern * X
C = 1#T_intern * H
Psi = 1#T_intern * psi

#Calculating a
a = A/C

#Results
print("\nPsi=", Psi,"\nA=",A,"\nC=",C)
print("\na=", a)
print("\npsi=", psi,"\nX=",X,"\nH=",H)

return a
