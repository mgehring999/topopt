import solidspy.assemutil as ass
import numpy as np
import scipy as sp

# local packages
from assembler import assemble

# pyomo import
from pyomo.environ import *

def assemble_model(nodes,elements,loads,volfrac):
	nelem = elements.shape[0]
	nnodes = nodes.shape[0]

	# assembly operator
	DME , IBC , neq = ass.DME(nodes, elements)

	# construct inital array for optimiziation var
	x = [1]*nelem
	V0 = sum(x)
	vol = volfrac*V0

	# material array (stiffness,nu)
	mats = np.ones((nelem,2))
	mats[:,1] *= 0.3

	# Stiffness Matrix
	Kglob_init = assemble(elements, mats, nodes, neq, DME,x)

	# load vector
	F_init = ass.loadasem(loads, IBC, neq)

	"""
	Pyomo Rules
	"""
	# volume fraction
	def vol_rule(m,vol): 
		return sum([m.x[j] for j in m.elems]) == vol 

	# FE equation
	def FKU_rule(m, i,penal=3):

		# update material definition
		x = [m.x[elem]**penal for elem in m.elems]

		# update stiffness matrix and assign to model
		m.K = assemble(elements, mats, nodes, neq, DME,x)

		return sum([m.K.value[i][j]*m.u[j] for j in m.eq]) == m.F[i]

	# compliance rule for objective
	def comp_rule(m,elements,mats,nodes,neq,DME,penal=3):

		# update material definition
		x = [m.x[elem]**penal for elem in m.elems]

		# update stiffness matrix and assign to model
		m.K = assemble(elements, mats, nodes, neq, DME,x)

		return sum([sum([m.K.value[row][col]*m.u[col] for col in m.eq])*m.u[row] for row in m.eq])

	"""
	Pyomo Topo Model
	"""
	model = ConcreteModel(name="topo")

	model.elems = Set(initialize=elements[:,0],domain=NonNegativeIntegers)
	model.eq = Set(initialize=range(neq),domain=NonNegativeIntegers)

	model.x = Var(model.elems,bounds=(1e-5,1),initialize=1)
	model.K = Param(initialize=Kglob_init,default=0,mutable=True,within=Any) 
	model.F = Param(model.eq,initialize=dict(enumerate(F_init)),within=Any)
	model.u = Var(model.eq,initialize=dict(enumerate(sp.sparse.linalg.spsolve(Kglob_init,F_init))))

	model.FKU_con = Constraint(model.eq, rule=FKU_rule)
	model.vol_con = Constraint(rule=vol_rule(model,vol))
	model.obj = Objective(expr=comp_rule(model,elements,mats,nodes,neq,DME),sense=minimize)

	return model