import pandas as pd
import numpy as np
import math
import cplex as cp
from scipy import sparse

def get_v_old():
    #dict = {1:1 , 2:2 , 3:0 , 4:0 , 5:0,
    #        6:0 , 7:0 , 8:0 , 9:0 , 10:0,
    #        11:0 ,12:0, 13:0, 14:0, 15:0,}
#   flux_dict = {1:5 , 2:0 , 3:0 , 4:4 , 5:0,
#                6:1 , 7:0 , 8:1 , 9:1 , 10:1,
#                11:1 ,12:1, 13:0, 14:0, 15:1,}
    flux_dict = {1:2 , 2:0 , 3:0 , 4:0.1 , 5:1.9,
                6:0 , 7:0 , 8:0 , 9:0 , 10:0,
                11:2 ,12:0, 13:0, 14:0, 15:0,}
    flux_list = list(flux_dict.values())
    v_old = np.array([float(i) for i in flux_list])
    return v_old

def get_N():
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
    N = np.array(network)
    return N


def get_N_internal(N):
    #External metabolites row index: 4,5,6,7,8,9,10,11
    N_internal = np.array(N[0:4, :])
    return N_internal



def get_v_new(v_i_old, v_old, index, N, delta):
    if delta > v_i_old:
        delta = v_i_old
    v_old_ub = np.array(v_old)
    v_old_ub[index] = v_i_old - delta
    v_old_lb = np.array([0.0]*len(v_old))
    v_old_lb[index] = v_i_old - delta

    #
    my_obj     = [float(i) for i in [0.0,0,0,0,0,0,0,0,0,0,1,1,0,0,1]] #value of variables: 1*v11_new + 1*v12_new + 1*v15_new = BM
    my_var     = ["v1","v2","v3","v4","v5","v6","v7","v8",
                  "v9","v10","v11","v12","v13","v14","v15"] #names of the variables
    my_ub      =  v_old_ub #All vi's need to be 0<=vi_new<=vi_old &&&&&&  c1_1_at_index <= c1_value
    my_lb      =  v_old_lb
    lp = cp.Cplex()
    lp.objective.set_sense(lp.objective.sense.maximize)
    lp.variables.add(names = my_var, obj = my_obj,
                        ub = my_ub, lb = my_lb)

    #change format N so can be used in constraints
    sN = sparse.csr_matrix(N)
    sN = sN.tocoo()
    rows = list(int(i) for i in sN.row) #geen idee waarom, maar dit geeft een error
    cols = list(int(i) for i in sN.col) # dit ook, moet dus met de hand overgetyped worden
    vals = list(int(i) for i in sN.data) #error ligt niet aan dat dit neg waarde bevat
    #print(rows, "\n", cols, "\n", vals, "\n")

    #set constraints
    #number of constraints = number of rows in N
    my_rhs     = [0.0] * len(N)
    my_sense   = "E" * len(N)
    lp.linear_constraints.add(rhs = my_rhs , senses = my_sense)
    lp.linear_constraints.set_coefficients(zip(rows, cols, vals))
    #disable printing
    lp.set_log_stream(None)
    lp.set_error_stream(None)
    lp.set_warning_stream(None)
    lp.set_results_stream(None)
    #solve
    lp.write("test" + str(index) + ".lp")
    lp.solve()
    x = lp.solution.get_values()
    #dj = lp.solution.get_reduced_costs()
    return x


#####################___main__program___####################################
v_old = get_v_old()
N = get_N()
N_internal = get_N_internal(N)
vbm_old = v_old[10]+ v_old[11]+ v_old[14]
#print(N, "\n")
#print(N_internal, "\n")
#print(N_internal * v_old, "\n")
#print(np.sum(N_internal * v_old))
#print known info:
print(v_old, "= original v (v_old)\n")

#for each v_i_old in v_old, find v_i_new
delta = 0.01
index = 0
controls = [0.00] * len(v_old)
for v_i_old in v_old:
    #delta = delta_perc * v_i_old
    if v_i_old > 0:
        v_new = get_v_new(v_i_old, v_old, index, N_internal, delta)
        vbm_new = v_new[10]+ v_new[11]+ v_new[14]
        print("\n", "Control of v", index+1, "with flux of", v_i_old, "  with delta of",delta,"\n",
                v_new, "= v_new \n",
                vbm_new, "= vbm_new \n",
                vbm_old, "= vbm_old \n",
                v_new[index], "= v_i_new \n",
                v_i_old, "= v_i_old \n")
        control = (v_i_old/vbm_old)*((vbm_old-vbm_new)/(v_i_old-v_new[index]))
        controls[index] = round(control, 5)
    #+1 to count for next v_i_old
    index = index +1

for ind, c in enumerate(controls):
    print("C",ind+1, " = ", c)
