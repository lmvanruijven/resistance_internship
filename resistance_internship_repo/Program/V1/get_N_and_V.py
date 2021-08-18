import numpy as np

def get_N_and_V(network,flux):
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
    return [N,V]
