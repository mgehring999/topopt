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

# load vector
RHSG = ass.loadasem(loads, IBC, neq)

# deformations
U = sp.sparse.linalg.spsolve(KG,RHSG)

# volume fraction
V0 = nelem
V = sum(x)

# compliance
r = np.matmul(KG.todense(),U)
c = np.matmul(U,r.T)

"""
Pyomo Topo Model
"""
model = pyo.ConcreteModel()