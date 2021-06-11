import solidspy.assemutil as ass
import solidspy.solutil as sol
import solidspy.postprocesor as pos
import solidspy.preprocesor as pre
import matplotlib.pyplot as plt
import numpy as np

folder = "./examples/mesh/"
nodes, _, elements, loads = pre.readin(folder=folder)

# assembly operator
DME , IBC , neq = ass.DME(nodes, elements)
nele = len(elements[:,0])

# construct material array for topology optimiziation
volfrac = 0.5
nu = 0.3
elem_E0 = np.ones((nele,1))*volfrac
elem_nu = np.ones((nele,1))*nu
mats = np.concatenate((elem_E0,elem_nu),axis=1)
print(mats)

# System assembly
KG = ass.assembler(elements, mats, nodes, neq, DME)
RHSG = ass.loadasem(loads, IBC, neq)

# solution
UG = sol.static_sol(KG, RHSG)

# calc compliance
KGdense = KG.todense()
r = np.matmul(KGdense,UG)
r = r.T
print("rechts",r)
c = np.matmul(UG,r)
print("compliance",c)

# post processing
UC = pos.complete_disp(IBC, nodes, UG)
print(UC)
print(UC.flatten())
pos.fields_plot(elements, nodes, UC)

# output
print(UC)
plt.show()