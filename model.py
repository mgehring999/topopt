import solidspy.assemutil as ass
import solidspy.solutil as sol
import solidspy.postprocesor as pos
import solidspy.preprocesor as pre
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import pyomo.environ as pyo

# local packages
from mesh import rect_mesh

"""
construct rectangular mesh
"""
nodes,elements,loads = rect_mesh(20)
nelem = elements.shape[0]
nnodes = nodes.shape[0]

print("++++++ node array +++++++\n",nodes)
print("++++++ load array +++++++\n",loads)
print("++++++ element array +++++++\n",elements)

# assembly operator
DME , IBC , neq = ass.DME(nodes, elements)

"""
initial parameters for pyomo model
"""
volfrac = 1
xmin = 1e-3
nu = 0.3
penal = 3

# construct inital material array for topology optimiziation
x = np.ones((nelem,1))*volfrac

# Stiffness Matrix
elem_nu = np.ones((nelem,1))*nu
mats = np.concatenate((x,elem_nu),axis=1)
KG = ass.assembler(elements, mats, nodes, neq, DME)

# load vector
F = ass.loadasem(loads, IBC, neq)

"""
Pyomo Rules
"""
# volume fraction
def vol_rule(m,volfrac):
	return sum([m.x[j] for j in m.elems])==volfrac

# FE equation
def FKU_rule(m, i):
    return sum([m.K[(i,j)]*m.u[j] for j in m.eq]) == m.F[i]

# compliance rule for objective
def comp_rule(m,elements,nodes,neq,DME,penal):

	# update material definition
	nelem = len(m.elems)
	x = np.array([m.x[i].value**penal for i in m.elems]).reshape((nelem,1))
	elem_nu = np.ones((nelem,1))*0.3

	# update stiffness matrix
	mats = np.concatenate((x,elem_nu),axis=1)
	KG = ass.assembler(elements, mats, nodes, neq, DME)

	# assign valus to model
	for i,val in KG.todok().items():
		m.K[i] = val

	return sum([sum([m.K[(i,j)]*m.u[j] for j in m.eq])*m.u[i] for i in m.eq])

"""
Pyomo Topo Model
"""

model = pyo.ConcreteModel(name="topo")

model.elems = pyo.Set(initialize=elements[:,0],domain=pyo.NonNegativeIntegers)
model.eq = pyo.Set(initialize=range(neq),domain=pyo.NonNegativeIntegers)

model.x = pyo.Var(model.elems,bounds=(xmin,1),initialize=volfrac)
model.K = pyo.Param(model.eq,model.eq,initialize=dict(KG.todok()),default=0,mutable=True)
model.F = pyo.Param(model.eq,initialize=dict(enumerate(F)))
model.u = pyo.Var(model.eq,initialize=dict(enumerate(sp.sparse.linalg.spsolve(KG,F))))

model.FKU_con = pyo.Constraint(model.eq, rule=FKU_rule)
model.vol_con = pyo.Constraint(rule=vol_rule(model,volfrac))
model.obj = pyo.Objective(expr=comp_rule(model,elements,nodes,neq,DME,penal),sense=1)

if __name__ == "__main__":
	if nnodes < 30:
		model.pprint()