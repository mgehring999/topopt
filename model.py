import solidspy.assemutil as ass
import solidspy.solutil as sol
import solidspy.postprocesor as pos
import solidspy.preprocesor as pre
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import pyomo.environ as pyo

"""
load meshed domain
"""
folder = "./examples/mesh/"
nodes, _, elements, loads = pre.readin(folder=folder)
nelem = elements.shape[0]
nnodes = nodes.shape[0]

# assembly operator
DME , IBC , neq = ass.DME(nodes, elements)

"""
Variables and parameters
"""
# construct material array for topology optimiziation
x = np.ones((nelem,1))

# Stiffness Matrix
nu = 0.3
elem_nu = np.ones((nelem,1))*nu
mats = np.concatenate((x,elem_nu),axis=1)
KG = ass.assembler(elements, mats, nodes, neq, DME)
print(elements)

# load vector
RHSG = ass.loadasem(loads, IBC, neq)

# deformations
U = sp.sparse.linalg.spsolve(KG,RHSG)


V0 = nelem
V = sum(x)

# compliance
r = np.matmul(KG.todense(),U)
c = np.matmul(U,r.T)

"""
initial parameters for pyomo model
"""
volfrac = 0.5
xmin = 0.2
"""
Pyomo Rules
"""
# volume fraction
def vol_rule(m,volfrac):
	return sum([m.x[j] for j in m.elems])==volfrac

# FE equation
def FKU_rule(m, i):
    return sum([m.K[(i,j)]*m.u[j] for j in m.eq]) == m.F[i]


"""
Pyomo Topo Model
"""

model = pyo.ConcreteModel(name="topo")

model.elems = pyo.Set(initialize=elements[:,0],domain=pyo.NonNegativeIntegers)
model.eq = pyo.Set(initialize=range(neq),domain=pyo.NonNegativeIntegers)

model.x = pyo.Var(model.elems,bounds=(xmin,1),initialize=volfrac/nelem)
model.K = pyo.Param(model.eq,model.eq,initialize=dict(KG.todok()),default=0,mutable=True)
model.F = pyo.Param(model.eq,initialize=dict(enumerate(RHSG)))
model.u = pyo.Var(model.eq,initialize=1)

model.FKU_con = pyo.Constraint(model.eq, rule=FKU_rule)
model.vol_con = pyo.Constraint(rule=vol_rule(model,volfrac))

model.pprint()