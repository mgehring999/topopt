import solidspy.preprocesor as pre
import numpy as np

def rect_mesh(ndiv,load_coord=(1,0),load_direct=(0,-1)):
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

	# init loads array, no loads initially
	load_init = np.zeros((nnodes,2))

	# concat to SolidsPy format
	nodes = np.concatenate((nnum,x,y,bc_init),axis=1)
	loads = np.concatenate((nnum,load_init),axis=1)

	# write bcs and loads in array
	for bcrow,lrow in zip(nodes,loads):
		xcoord = bcrow[1]
		ycoord = bcrow[2]
		if np.isclose(xcoord,-1):
			bcrow[3:] = [-1,-1]
		if np.isclose(xcoord,load_coord[0]) and np.isclose(ycoord,load_coord[1]):
			lrow[1] = load_direct[0]
			lrow[2] = load_direct[1]

	return nodes,elements,loads
