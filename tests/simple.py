import pyomo.environ as pyo

model = pyo.ConcreteModel()

model.x = pyo.Var([1,2],domain=pyo.NonNegativeReals)
model.OBJ = pyo.Objective(expr=model.x[1]**2+model.x[2]**2)

model.Constraint1 = pyo.Constraint(expr=model.x[1] + 2 * model.x[2] == 2)

