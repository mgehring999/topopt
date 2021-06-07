from simple import model
from pyomo.opt import SolverFactory

model.display()
opt = SolverFactory("scip")
result = opt.solve(model)
print(result)