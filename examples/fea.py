#!/usr/bin/env python3

import solidspy.assemutil as ass
import solidspy.solutil as sol
import solidspy.postprocesor as pos
import solidspy.preprocesor as pre
import matplotlib.pyplot as plt
import numpy as np

def rect_mesh(ndiv):
	"""
	ndiv: int, number of division along each side of the domain
	load_coord: (x,y), coordinate of single force load
	load_direct: (fx,fy), load vector
	"""

	# generate rectangular grid
	x,y,elements=pre.rect_grid(2,2,ndiv,ndiv)

	# mesh properties
	nelem = elements.shape[0]
	nnodes = len(x)

	# reshape to column vector
	x = x.reshape((nnodes,1))
	y = y.reshape((nnodes,1))


	# node numbering
	nnum = np.arange(nnodes)
	nnum = nnum.reshape((nnodes,1))

	# init bc array, no bcs initially
	bc_init = np.zeros((nnodes,2))

	# concat to SolidsPy format
	nodes = np.concatenate((nnum,x,y,bc_init),axis=1)

	# write bcs and loads in array
	for bcrow in nodes:
		xcoord = bcrow[1]
		ycoord = bcrow[2]
		if np.isclose(xcoord,-1):
			bcrow[3:] = [-1,-1]

	return nodes,elements

def add_point_force(nodes,load_coord,load_direction):
	
	# init loads array, no loads initially
	nnodes = len(nodes)
	load_init = np.zeros((nnodes,2))
	nnum = np.arange(nnodes)

	loads = np.concatenate((nnum.reshape((nnodes,1)),load_init),axis=1)

	for node,lrow in zip(nodes,loads):
		xcoord = node[1]
		ycoord = node[2]
		if np.isclose(xcoord,load_coord[0]) and np.isclose(ycoord,load_coord[1]):
			lrow[1] = load_direction[0]
			lrow[2] = load_direction[1]

	# test if loads were added
	if not loads[:,1:3].any():
		raise ValueError("No loads were applied, check ndivs or load coordinates")

	return loads

nodes,elements = rect_mesh(6)
nelem = elements.shape[0]
nnodes = nodes.shape[0]

# add loads 
loads = add_point_force(nodes,(1,0),(0,-1))

print("++++++ node array +++++++\n",nodes)
#print("++++++ load array +++++++\n",loads)
#print("++++++ element array +++++++\n",elements)

# assembly operator
DME , IBC , neq = ass.DME(nodes, elements)

print(DME[1,:8])
print("++++++++++++++++++++++++++++",IBC)
# construct material array for topology optimiziation
volfrac = 1
nu = 0.3
elem_E0 = np.ones((nelem,1))*volfrac
elem_nu = np.ones((nelem,1))*nu
mats = np.concatenate((elem_E0,elem_nu),axis=1)

# System assembly
KG = ass.assembler(elements, mats, nodes, neq, DME)
RHSG = ass.loadasem(loads, IBC, neq)

# solution
UG = sol.static_sol(KG, RHSG)

# calc compliance
KGdense = KG.todense()
r = np.matmul(KGdense,UG)
r = r.T
c = np.matmul(UG,r)

# post processing
UC = pos.complete_disp(IBC, nodes, UG)
# pos.fields_plot(elements, nodes, UC)

# output
print(UC)

# stresses
strains,stresses=pos.strain_nodes(nodes,elements,mats,UC)
seqv = np.zeros((nnodes,1))
for idx,(sx,sy,txy) in enumerate(stresses):
    vMises = np.sqrt(sx**2+sy**2-sx*sy+3*txy**2)
    seqv[idx]=vMises
pos.fields_plot(elements,nodes,UC,S_nodes=stresses)

plt.show()