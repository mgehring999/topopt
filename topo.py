from model import model
from pyomo.opt import SolverFactory

model.display()
opt = SolverFactory("scip")
result = opt.solve(model,tee=True)
print(result)