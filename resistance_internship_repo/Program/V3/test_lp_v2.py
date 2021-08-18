import cplex
#from cplex.exceptions import CplexError
import sys



prob = cplex.Cplex()

my_obj      = [1.0, 0.0, 1.0]
my_ub       = [1.0] * len(my_obj)
my_lb       = [0.0] * len(my_obj)
my_colnames = ["a", "b", "c"]
prob.objective.set_sense(prob.objective.sense.minimize)
prob.variables.add(obj = my_obj, ub = my_ub, lb = my_lb ,names = my_colnames)

#combined constraints
my_rhs      = [20.0, 30.0, 1.0, 1.0]
#my_rownames = ["c1", "c2", "e1", "e2"]
my_sense    = "E" * len(my_rhs)
rows = [0,1,1,2,2,3]
cols =  [0,1,2,0,1,2]
vals = [20,20,30,1,1,1]
print(list(zip(rows, cols, vals)))
prob.linear_constraints.add(rhs = my_rhs, senses = my_sense)
prob.linear_constraints.set_coefficients(zip(rows, cols, vals))

prob.write("lpex1.lp")

prob.solve()

x = prob.solution.get_values()
#dj = lp.solution.get_reduced_costs()
print(x,"\n")
