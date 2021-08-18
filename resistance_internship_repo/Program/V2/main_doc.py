import pandas as pd
import numpy as np
import math

#dict = {1:1 , 2:2 , 3:0 , 4:0 , 5:0,
#        6:0 , 7:0 , 8:0 , 9:0 , 10:0,
#        11:0 ,12:0, 13:0, 14:0, 15:0,}
dict = {1:1 , 2:0 , 3:0 , 4:0.1 , 5:0.9,
        6:0 , 7:0 , 8:0 , 9:0.1 , 10:0,
        11:0.9 ,12:0, 13:0, 14:0, 15:0,}


#Stoichiometry matrix (does not change)
#network can only exist of -1,0,1 !
network = np.array([[1,-1,-1,-1,-1,-1,1,0,0,0,0,0,0,0,0],
                   [0,0,0,1,1,0,1,-1,-1,-1,-1,0,0,0,0],
                   [0,0,0,0,0,1,0,0,1,-1,0,-1,-1,-1,0],
                   [0,0,0,0,0,0,0,0,0,1,0,0,0,0,-1],
                   [0,0,0,0,0,0,0,0,0,0,1,1,0,0,1],
                   [-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
                   [0,1,0,0,0,0,0,0,0,0,0,0,0,0,0],
                   [0,0,1,0,0,0,0,0,0,0,0,0,0,0,0],
                   [0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0],
                   [0,0,0,0,0,0,0,1,0,0,0,0,0,0,0],
                   [0,0,0,0,0,0,0,0,0,0,0,0,1,0,0],
                   [0,0,0,0,0,0,0,0,0,0,0,0,0,1,0]])
flux = list(dict.values())

#determine index of fluxes that are zero
V_pos = np.array(flux)
V_pos[V_pos>0]=1
N_pos = np.array(flux * network)
N_pos[N_pos<0]= 1
col_sum = np.sum(N_pos, axis = 0)
idx = np.flatnonzero(col_sum)
#Create N and V without the fluxes that are zero
N = np.array(network)
N = N[:, idx]
V = np.array(flux)
V = V[V != 0]
#Determine probabilities insteat of fluxes for V
V = V / np.sum(V)
print(N, "=N\n\n",V, "=V")

#Calculation H(Pij)(k=1)
H = np.sum(-1 * V * np.log2(V))

#Calculation X(Pij)
#Make matrixes with only outflow (Pi) and only inflow (Pj)
outf = np.array(N*V)
outf[outf >0] = 0
inf = np.array(N*V)
inf[inf <0] = 0
#print(outf, "= outf\n")
#print(inf, "= inf \n")
X = 0
i = 0
for v in V:
    d_out = np.array(np.transpose(N)[i])
    d_out[d_out >0] =0
    d_in = np.array(np.transpose(N)[i])
    d_in[d_in <0] = 0
    #print(d_out, "is d_out \n")
    #print(d_in, "is d_in \n")
    b = np.sum(d_out*np.transpose(outf))
    c = np.sum(d_in*np.transpose(inf))
    #print(b,c, "= b (i) and c (j)\n")
    #print("X =", V[i], "log2 ", V[i], "/", b, "* ",c)
    z = V[i]* math.log2(V[i]/(b*c))
    X = X + (z)
    i = i + 1

#Calculation Psi(Pij)
psi = H - X

#Calculating flux that remains internal and external
##P3 - sum of P1,2,4,5,6,7
#T_intern = flux[2] - np.sum(flux[0:2]) - np.sum(flux[6:8]) - np.sum(flux[11:14])
#T_extern = flux[2] - T_intern
#print("T_intern=",T_intern,"\nT_extern=",  T_extern)

#Calculating A, C and \Psi
#A = T_intern * X
#C = T_intern * H
#Psi = T_intern * psi

#Calculating Fitnes
a = X/H
F = -1*a*math.log2(a)

#Results
#print("\nPsi=", Psi,"\nA=",A,"\nC=",C)
print("\na=", a,"\nF=", F)
print("\npsi=", psi,"\nX=",X,"\nH=",H)
