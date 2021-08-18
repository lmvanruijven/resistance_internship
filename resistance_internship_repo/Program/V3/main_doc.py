import pandas as pd
import numpy as np
import math
import cplex as cp

def get_v_old():
    #dict = {1:1 , 2:2 , 3:0 , 4:0 , 5:0,
    #        6:0 , 7:0 , 8:0 , 9:0 , 10:0,
    #        11:0 ,12:0, 13:0, 14:0, 15:0,}
    flux_dict = {1:3 , 2:0 , 3:0 , 4:2 , 5:0,
                6:1 , 7:0 , 8:1.5 , 9:0 , 10:0,
                11:0.5 ,12:1, 13:0, 14:0, 15:0,}
    flux_list = list(flux_dict.values())
    #Create V without the fluxes that are zero
    v_old = np.array(flux_list)
    #v_old = v_old[v_old != 0]
    #Determine probabilities insteat of fluxes for V
    #v_old = v_old / np.sum(v_old)
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
    N_internal = np.array(N[4:11, :])
    return N_internal



def get_v_i_new(v_i_old, v_old, index):
    delta = 0.01
    c1_value = v_i_old - delta
    print(c1_value)
    ##Make variables
    my_obj     = [0,0,0,0,0,0,0,0,0,0,1,1,0,0,1] #value of variables: 1*v11_new + 1*v12_new + 1*v15_new = BM
    my_var     = ["v1","v2","v3","v4","v5","v6","v7","v8",
                  "v9","v10","v11","v12","v13","v14","v15"] #names of the variables
    #make my_ub. All vi's need to be 0<=vi_new<=vi_old
    my_ub      = v_old #[40.0, cp.infinity, cp.infinity] #bounds of the variables so: 0<=x1<=40
    #make my_rhs. V_index_new = v_i_old - delta
    my_rhs     = [c1_value, cp.infinity] #Subject to values. So 2 values because 2 constraints.
    my_con     = ["c1", "c2"] #names of the constraints
    my_sense   = "LL"

    #Make lp problem
    lp = cp.Cplex()
    lp.objective.set_sense(lp.objective.sense.maximize)
    #lp.variables.add(obj = obj, ub = ub, names = var)
    lp.variables.add(names = my_var, obj = my_obj, ub = my_ub)
    #p.linear_constraints.add(names = my_con, rhs = my_rhs , senses = my_sense,lin_expr =  rows)
    c1_1_at_index = np.zeros(15)
    c1_1_at_index[index] = 1
    rows = [   [[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14],c1_1_at_index]  ,
               [[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]]]
    lp.linear_constraints.add(names = my_con, rhs = my_rhs , senses = my_sense,
                           lin_expr =  rows)

    #Get v_i_new from the new objective solution
    lp.write("test" + str(index) + ".lp")
    lp.solve()

    x = lp.solution.get_values()
    dj = lp.solution.get_reduced_costs()
    print(x,"\n")

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
print(v_old, "= original v (v_old)","\nvbm_old =",vbm_old, "\n")

#for each v_i_old in v_old, find v_i_new
index = 0
v_new = np.array([])
for v_i_old in v_old:
    if v_i_old >0:
        print(v_i_old)
        v_i_new = get_v_i_new(v_i_old, v_old, index)
        print("NEXT V_I_OLD !!!!!!!!!!!!!!!!!!!! \n\n")
    #add v_i_new to v_new
    #v_new.append(1)
    #+1 to count for next v_i_old
    index = index +1


#for the found v_new, calculate Ci
for v_i_old in v_old:
    vbm_new = np.sum(v_new[10:11:14])
    Ci = 1
    #add Ci to list with all_C
