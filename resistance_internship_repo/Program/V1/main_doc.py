#import pandas as pd
import numpy as np

#dict = {1:1 , 2:2 , 3:0 , 4:0 , 5:0,
#        6:0 , 7:0 , 8:0 , 9:0 , 10:0,
#        11:0 ,12:0, 13:0, 14:0, 15:0,}
dict = {1:2 , 2:0 , 3:0 , 4:1 , 5:0,
        6:0 , 7:0 , 8:1 , 9:0 , 10:1,
        11:0 ,12:0, 13:0, 14:0, 15:0,}




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

#Get N and V
from get_N_and_V import get_N_and_V
N_and_V = get_N_and_V(network, flux)
N = N_and_V[0]
V = N_and_V[1]
print(V, "= V \n\n", N, "= N \n")

#Get i , j and i,j
from get_i_j import get_i_j
i_j_list = get_i_j(N, V, flux)
i = i_j_list[0]
j = i_j_list[1]
i_j = i_j_list[2]


#calculate results
from calculate_results import calculate_results
results = calculate_results(i, j, i_j)
H = results[0]
X = results[1]
psi = results[2]

#Results
print("\npsi=", psi,"\nI(i,j)=",X,"\nH(i,j)=",H)
