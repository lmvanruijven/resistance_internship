#import cbmpy as cbm
import cplex as cp

# data common to all populateby functions
my_obj     = [1.0, 2.0, 3.0] #value of variables: 1*x1 + 2*x2 + 3*x3
my_var     = ["x1", "x2", "x3"] #names of the variables
my_ub      = [40.0, cp.infinity, cp.infinity] #bounds of the variables so: 0<=x1<=40
my_rhs     = [20.0, 30.0] #Subject to values. So two constraints because 2 values.
my_con     = ["c1", "c2"] #names of the constraints
my_sense   = "LL"

lp = cp.Cplex()

lp.objective.set_sense(lp.objective.sense.maximize)
#lp.variables.add(obj = obj, ub = ub, names = var)
lp.variables.add(names = my_var, obj = my_obj, ub = my_ub )

#p.linear_constraints.add(names = my_con, rhs = my_rhs , senses = my_sense,lin_expr =  rows)
rows = [[[0,1,2],[-1.0, 1.0,1.0]],
        [[0,1,2],[ 1.0,-3.0,1.0]]]

lp.linear_constraints.add(names = my_con, rhs = my_rhs , senses = my_sense,
                       lin_expr =  rows)

#print(lbs, "\n", ub1, "\n", names, "\n", rows)
#lp.maximize(my_obj)
lp.write("test.lp")
lp.solve()

x = lp.solution.get_values()
dj = lp.solution.get_reduced_costs()
print(x,"\n",dj)


#slack = lp.solution.get_linear_slacks()
#pi = lp.solution.get_dual_values()
#total_result = lp.solution.get_objective_value()
#solution_value = lp.solution.get_objective_value()
#print(total_result,"\n", solution_value,"\n",slack,"\n",pi,"\n",x,"\n",dj)
