from solidspy.preprocesor import rect_grid
from topopt.logging import TopOptLogger
import numpy as np

import logging
logger = logging.getLogger('topopt')

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
    
        self.logger_handler = TopOptLogger("debug")

    def rect_mesh(self,ndiv: int):
        """
        defines quadratic domain meshed with QUAD4 elements
        
        Paramters:
        ndiv: number of elements on each side of the domain
        """
        logger.info("started meshing")
        if ndiv < 1:
            ValueError("ndiv must be greater or equal to one")
        x,y,self.elements = rect_grid(2,2,ndiv,ndiv)

        self.nnodes = len(x)
        self.nelem = len(self.elements[:,0])

        x = x.reshape((self.nnodes,1))
        y = y.reshape((self.nnodes,1))
        
        # node numbering
        nnum = np.arange(self.nnodes)
        self.nnum = nnum.reshape((self.nnodes,1))

        # init bc array, no bcs initially
        bc_init = np.zeros((self.nnodes,2))

        # concat nodes array
        self.nodes = np.concatenate((self.nnum,x,y,bc_init),axis=1)

        # log
        logger.info("initialized mesh")


    def load_mesh(self):
        # TODO
        pass

class BoundaryCondition:
    def __init__(self):
        self.bc_by_coord = []
        self.bc_by_comp = []
        self.bc_by_edge = []

    def add_by_coord(self,coordinate,direction):
        """adds boundary condition to single node via cooridnate
        ## Parameters
        coordinate: tuple (x,y) 
        direction: tuple (xmag,ymag)
        """
        self.bc_by_coord.append((coordinate,direction))
	
    def add_by_component(self):
        # TODO: implement when mesh interface is built
        pass
    def add_by_edge(self,edge,dof):
        """adds boundary condition to whole edge via keywords

        ## Parameters
        edge: `str`: "left", "right", "top", "bot"
        dof: `int` or `str`: 1,2 or "all"
        """
        self.bc_by_edge.append((edge,dof))

    def apply(self,mesh):
        # Dummy method that gets overwritten
        pass

class Load(BoundaryCondition):
    def __init__(self, mesh):
        super().__init__()

        load_columns = np.zeros((mesh.nnodes,2))
        node_numbering_column = mesh.nnum.reshape((mesh.nnodes,1))
        self.loads = np.concatenate((node_numbering_column,load_columns),axis=1)

    def apply(self,mesh):
        # iterate over all loads to be applied
        for load_coord,load_direction in self.bc_by_coord:
            # iterate over loads array 
            for node,lrow in zip(mesh.nodes,self.loads):
                xcoord = node[1]
                ycoord = node[2]
                if np.isclose(xcoord,load_coord[0]) and np.isclose(ycoord,load_coord[1]):
                    lrow[1] = load_direction[0]
                    lrow[2] = load_direction[1]


class Displacement(BoundaryCondition):
    def __init__(self):
        super().__init__()

    def apply(self,mesh):
        for edge,dof in self.bc_by_edge:
            mesh = self._add_support(mesh,edge=edge,dof=dof)
        return mesh

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

    def _add_support(self,mesh,edge = "left",dof = "all"):
        # copied code ..
        if edge == "left":
            if dof is not int or list:
                ValueError("DOFs are specified as int or list of ints")

            if dof == "all":
                dof = [1,2]
                for i in dof:
                    mesh.nodes = self._set_bc(mesh.nodes,"x",-1,i)
            else:
                mesh.nodes = self._set_bc(mesh.nodes,"x",-1,dof)

        if edge == "right":
            if dof is not int or list:
                ValueError("DOFs are specified as int or list of ints")

            if dof == "all":
                dof = [1,2]
                for i in dof:
                    mesh.nodes = self._set_bc(mesh.nodes,"x",1,i)
            else:
                mesh.nodes = self._set_bc(mesh.nodes,"x",1,dof)

        if edge == "top":
            if dof is not int or list:
                ValueError("DOFs are specified as int or list of ints")

            if dof == "all":
                dof = [1,2]
                for i in dof:
                    mesh.nodes = self._set_bc(mesh.nodes,"y",1,i)
            else:
                mesh.nodes = self._set_bc(mesh.nodes,"y",1,dof)

        if edge == "bot":
            if dof is not int or list:
                ValueError("DOFs are specified as int or list of ints")

            if dof == "all":
                dof = [1,2]
                for i in dof:
                    mesh.nodes = self._set_bc(mesh.nodes,"y",-1,i)
            else:
                mesh.nodes = self._set_bc(mesh.nodes,"y",-1,dof)
        return mesh