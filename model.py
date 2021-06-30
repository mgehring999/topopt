import solidspy.assemutil as ass
import solidspy.solutil as sol
import solidspy.postprocesor as pos
import solidspy.preprocesor as pre
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp

# local packages
from mesh import rect_mesh

# pyomo import
from pyomo.environ import *

"""
construct rectangular mesh
"""
nodes,elements,loads = rect_mesh(2)
nelem = elements.shape[0]
nnodes = nodes.shape[0]


# assembly operator
DME , IBC , neq = ass.DME(nodes, elements)

"""
initial parameters for pyomo model
"""
volfrac = 0.5
xmin = 1e-5
nu = 0.3
penal = 3

# construct inital material array for topology optimiziation
x = np.ones((nelem,1))#/nelem#*volfrac
V0 = x.sum()
vol = volfrac*V0

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
def vol_rule(m,vol): 
	return sum([m.x[j] for j in m.elems]) == vol 

# FE equation
def FKU_rule(m, i):
	return sum([m.K[(i,j)]*m.u[j] for j in m.eq]) == m.F[i]

# compliance rule for objective
def comp_rule(m,elements,nodes,neq,DME,penal):

	# update material definition
	nelem = len(m.elems)
	x = np.array([value(m.x[elem]**penal) for elem in m.elems]).reshape((nelem,1))
	elem_nu = np.ones((nelem,1))*0.3

	# update stiffness matrix
	mats = np.concatenate((x,elem_nu),axis=1)
	KG = ass.assembler(elements, mats, nodes, neq, DME)

	# assign values to model
	for elem,val in KG.todok().items():
		m.K[elem] = 0#val

	return sum([sum([m.K[(row,col)]*m.u[col] for col in m.eq])*m.u[row] for row in m.eq])

"""
Pyomo Topo Model
"""
model = ConcreteModel(name="topo")

model.elems = Set(initialize=elements[:,0],domain=NonNegativeIntegers)
model.eq = Set(initialize=range(neq),domain=NonNegativeIntegers)

model.x = Var(model.elems,bounds=(xmin,1),initialize=1)
model.K = Param(model.eq,model.eq,initialize=0,default=0,mutable=True) #dict(KG.todok())
model.F = Param(model.eq,initialize=dict(enumerate(F)))
model.u = Var(model.eq,initialize=dict(enumerate(sp.sparse.linalg.spsolve(KG,F))))

model.FKU_con = Constraint(model.eq, rule=FKU_rule)
model.vol_con = Constraint(rule=vol_rule(model,vol))
model.obj = Objective(expr=comp_rule(model,elements,nodes,neq,DME,penal),sense=1)

if __name__ == "__main__":

	if nnodes < 30:
		print("++++++ node array +++++++\n",nodes)
		print("++++++ load array +++++++\n",loads)
		print("++++++ element array +++++++\n",elements)
		model.pprint()