import numpy as np
import math

def calculate_results(i, j, i_j):
    #Calculation H(Pij)(k=1)
    H_i_j = np.sum(-1 * i_j * np.log2(i_j))
    #Calculating H(Pi)
    H_i = np.sum(-1 * i * np.log2(i))
    #Calculating H(Pj)
    H_j = np.sum(-1 * j * np.log2(j))

    #Calculating I(Pij)
    I_i_j = H_i + H_j - H_i_j
    #Calculating X(Pij)
    X_i_j = I_i_j

    #Calculating H(Pi | Pj)
    H_i_c_j = H_i_j - H_j
    #Calculating H(Pj | Pi)
    H_j_c_i = H_i_j - H_i

    #Calculating psi
    psi = H_i_j - X_i_j
    print("\nH(i)=", H_i,"\nH(j)=",H_j,"\nH(i|j)=",H_i_c_j,"\nH(j|i)=",H_j_c_i)

    return H_i_j, X_i_j, psi
