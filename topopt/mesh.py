import solidspy.preprocesor as pre
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

if __name__ == "__main__":
	nodes,elems = rect_mesh(3)
	print(nodes)

	loads = add_point_force(nodes,(0,1),(0,-1))
	print(loads)