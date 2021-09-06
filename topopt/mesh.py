from solidspy.preprocesor import rect_grid
import numpy as np

class Mesh:
    """Mesh creates or loads and saves the information about the meshed domain
    """
    def __init__(self):
        self.x = None
        self.y = None
        self.elements = None
        self.nodes = None
        self.nnum = None
        self.nnodes = None
    
    def rect_mesh(self,ndiv: int):
        """
        defines quadratic domain meshed with QUAD4 elements
        
        Paramters:
        ndiv: number of elements on each side of the domain
        """
        if ndiv < 1:
            ValueError("ndiv must be greater or equal to one")
        x,y,self.elements = rect_grid(2,2,ndiv,ndiv)

        self.nnodes = len(x)
        self.nelem = len(self.elements[:,0])

        x = x.reshape((self.nnodes,1))
        y = y.reshape((self.nnodes,1))
        
        # node numbering
        nnum = np.arange(self.nnodes)
        nnum = nnum.reshape((self.nnodes,1))

        # init bc array, no bcs initially
        bc_init = np.zeros((self.nnodes,2))

        # concat nodes array
        self.nodes = np.concatenate((nnum,x,y,bc_init),axis=1)

    def _set_bc(self,nodes,coord,value,dof):
        
        # write bcs and loads in array
        for row in nodes:
            if coord =="x":
                current_coord = row[1]
            elif coord == "y":
                current_coord = row[2]
            else:
                ValueError("unknown coordinate")

            # this assumes a 2x2 domain 
            if np.isclose(current_coord,value):
                if dof == 1:
                    row[3] = -1
                elif dof == 2:
                    row[4] = -1
                else:
                    ValueError("unknown DOF")
        return nodes

    def add_support(self,edge = "left",dof = "all"):
        # copied code ..
        # TODO: separate Mesh and BoundaryCondition
        if edge == "left":
            if dof is not int or list:
                ValueError("DOFs are specified as int or list of ints")

            if dof == "all":
                dof = [1,2]
                for i in dof:
                    self.nodes = self._set_bc(self.nodes,"x",-1,i)
            else:
                self.nodes = self._set_bc(self.nodes,"x",-1,dof)

        if edge == "right":
            if dof is not int or list:
                ValueError("DOFs are specified as int or list of ints")

            if dof == "all":
                dof = [1,2]
                for i in dof:
                    self.nodes = self._set_bc(self.nodes,"x",1,i)
            else:
                self.nodes = self._set_bc(self.nodes,"x",1,dof)

        if edge == "top":
            if dof is not int or list:
                ValueError("DOFs are specified as int or list of ints")

            if dof == "all":
                dof = [1,2]
                for i in dof:
                    self.nodes = self._set_bc(self.nodes,"y",1,i)
            else:
                self.nodes = self._set_bc(self.nodes,"y",1,dof)

        if edge == "bot":
            if dof is not int or list:
                ValueError("DOFs are specified as int or list of ints")

            if dof == "all":
                dof = [1,2]
                for i in dof:
                    self.nodes = self._set_bc(self.nodes,"y",-1,i)
            else:
                self.nodes = self._set_bc(self.nodes,"y",-1,dof)

    def load_mesh(self):
        # TODO
        pass

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

