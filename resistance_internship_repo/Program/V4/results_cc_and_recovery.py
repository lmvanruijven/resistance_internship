import pandas as pd
import numpy as np
import math
import cplex as cp
from scipy import sparse

def get_v_old(network):
    flux_list = list(network.values())
    v_old = np.array([float(i) for i in flux_list])
    return v_old

def get_N():
    #Stoichiometry matrix (does not change)
    #network can only exist of -1,0,1
    Stoich_network = np.array([[1,-1,-1,-1,-1,-1,1,0,0,0,0,0,0,0,0,0,-1,0,-1],
                       [0,0,0,1,1,0,1,-1,-1,-1,-1,0,0,0,0,0,0,1,1],
                       [0,0,0,0,0,1,0,0,1,-1,0,-1,-1,-1,0,1,0,0,0],
                       [0,0,0,0,0,0,0,0,0,1,0,0,0,0,-1,0,0,0,0],
                       [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,-1,0],
                       [0,0,0,0,0,0,0,0,0,0,1,1,0,0,1,0,0,0,0],
                       [-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
                       [0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
                       [0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
                       [0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0],
                       [0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0],
                       [0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0],
                       [0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0],
                       [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0]])
    N = np.array(Stoich_network)
    return N


def get_N_internal(N):
    #External metabolites row index: 5,6,7,8,9,10,11,12,13
    N_internal = np.array(N[0:5, :])
    return N_internal

def get_v_new(v_i_old, v_old, index, N, delta):
    if delta > v_i_old:
        delta = v_i_old
    v_old_ub = np.array(v_old)
    v_old_ub[index] = v_i_old - delta
    v_old_lb = np.array([0.0]*len(v_old))
    v_old_lb[index] = v_i_old - delta

    #
    my_obj     = [float(i) for i in [0.0,0,0,0,0,0,0,0,0,0,1,1,0,0,1,0,0,0,0]] #value of variables: 1*v11_new + 1*v12_new + 1*v15_new = BM
    my_var     = ["v1","v2","v3","v4","v5","v6","v7","v8",
                  "v9","v10","v11","v12","v13","v14","v15"
                  ,"v16","v17","v18","v19"] #names of the variables
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

    return x


#####################___main__program___####################################
network_1a = {1:2 , 2:0 , 3:0 , 4:2, 5:0,
                6:0 , 7:0 , 8:0 , 9:0 , 10:0,
                11:2 ,12:0, 13:0, 14:0, 15:0,
                16:0, 17:0, 18:0, 19:0}
network_1c = {1:2 , 2:0 , 3:0 , 4:2, 5:0,
                6:0 , 7:0 , 8:0.5 , 9:0 , 10:0,
                11:1.5 ,12:0, 13:0, 14:0, 15:0,
                16:0, 17:0, 18:0, 19:0}
network_1b = {1:2 , 2:0 , 3:0 , 4:1, 5:1,
                6:0 , 7:0 , 8:0 , 9:0 , 10:0,
                11:2 ,12:0, 13:0, 14:0, 15:0,
                16:0, 17:0, 18:0, 19:0}
network_7 = {1:2 , 2:0 , 3:0 , 4:1, 5:0,
                6:0 , 7:0 , 8:0 , 9:0 , 10:0,
                11:2 ,12:0, 13:0, 14:0, 15:0,
                16:0, 17:1, 18:1, 19:0}
network_8 = {1:2 , 2:0 , 3:0 , 4:2/3, 5:2/3,
                6:0 , 7:0 , 8:0 , 9:0 , 10:0,
                11:2 ,12:0, 13:0, 14:0, 15:0,
                16:0, 17:0, 18:0, 19:2/3}
network_13 = {1:2 , 2:0 , 3:0 , 4:1.1, 5:0.9,
                6:0 , 7:0 , 8:0 , 9:0 , 10:0,
                11:2 ,12:0, 13:0, 14:0, 15:0,
                16:0, 17:0, 18:0, 19:0}
network_4 = {1:2 , 2:0 , 3:0 , 4:0.5, 5:1.5,
                6:0 , 7:0 , 8:0 , 9:0 , 10:0,
                11:2 ,12:0, 13:0, 14:0, 15:0,
                16:0, 17:0, 18:0, 19:0}
network_12 = {1:1 , 2:0 , 3:0 , 4:1, 5:0,
                6:0 , 7:0 , 8:0 , 9:0 , 10:1,
                11:0 ,12:0, 13:1, 14:0, 15:1,
                16:2, 17:0, 18:0, 19:0 }

#set networks and delta's in list so iterating is possible
networks_list = [network_1b, network_4]# [network_1c, network_12, network_1a, network_4, network_1b, network_8] #[network_1a, network_4, network_13, network_1b, network_7, network_8] #[network_1c,network_12, network_1a, network_4, network_13, network_1b, network_7, network_8]
delta_perc_list = [0.1]#[0.99, 0.9, 0.8, 0.7, 0.6,0.5,0.4,0.3,0.2,0.1,0.01]#[0.99,0.6,0.5,0.1]
#iterate over desired delta percentages
for delta_perc in delta_perc_list:
    print("\n")
#iterate over desired networks
    for network in networks_list:
        v_old = get_v_old(network)
        N = get_N()
        N_internal = get_N_internal(N)
        vbm_old = v_old[10]+ v_old[11]+ v_old[14]

        index = 0
        #declair new arrays wich will contain results for all reactions of the network
        vbm_new__to_old_perc =  np.array([-1.00] * len(v_old) )
        CC = np.array([-1.00] * len(v_old))
        CC_quadr =  np.array([-1.00] * len(v_old) )
        P_i =  np.array([-1.00] * len(v_old) )
        F_i =  np.array([-1.00] * len(v_old) )
        recovery =  np.array([-1.00] * len(v_old) )
       #for each v_i_old in v_old, find v_i_new by solving LP problem
        for v_i_old in v_old:
            delta = delta_perc * v_i_old
            if v_i_old > 0:
                v_new = get_v_new(v_i_old, v_old, index, N_internal, delta)
                vbm_new = v_new[10]+ v_new[11]+ v_new[14]
                control = (v_i_old/vbm_old)*((vbm_old-vbm_new)/(v_i_old-v_new[index]))
                CC[index] = round(control, 5)
                CC_quadr[index] = round( (control* control), 5)
                P_i[index] = v_i_old/np.sum(v_old)
                recovery[index] = (v_i_old - v_new[index] )/vbm_new
                #print("vi_old:", v_i_old, "v_new[index]:", v_new[index], "vbm_new", vbm_new)
            #+1 to count for next v_i_old
            index = index + 1

        #Calculate results with P_i as weighting method
        P_i_CC_quadr = 0
        P_i_recovery = 0
        for ind, c in enumerate(CC):
            if c >= 0:
                '''print("C",ind+1, " = ", c, "    CC^2 = ",  CC_quadr[ind],
                    "pi = ", P_i[ind], "    Fi = ", F_i[ind] ,
                    "      recovery_i = ", recovery[ind] )'''
                P_i_CC_quadr = P_i_CC_quadr + P_i[ind] * CC_quadr[ind]
                P_i_recovery = P_i_recovery + P_i[ind] * recovery[ind]


        #print results without details
        print(np.average(CC[CC>=0])#, "\n",
            #np.average(CC_quadr[CC_quadr>=0])#, "\n"
            #,P_i_CC_quadr, "\n",
            #np.average(recovery[recovery >=0])#, " \n"
            #,P_i_recovery, " \n"
            )

        #print results with details
        '''print("\n", "Results for network:", network, "\n"
                ,np.average(CC[CC>=0]), "= average CC \n"
                ,np.average(CC_quadr[CC_quadr>=0]), "= average CC^2 \n"
                ,P_i_CC_quadr, "= sum of Vi/T.. * CC^2 \n"
                ,np.average(recovery[recovery >=0]), "= average recovery \n"
                ,P_i_recovery, "=Sum of (Vi/T..) * recovery  \n"
                )'''
