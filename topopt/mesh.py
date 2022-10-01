from typing import Iterable
from solidspy.preprocesor import rect_grid
from topopt.logging import TopOptLogger
import numpy as np

import logging
logger = logging.getLogger('topopt')

class Mesh:
    """Mesh creates or loads and saves the information about the meshed domain
    """
    def __init__(self):
        self.nodal_coords = None
        self.elements = None
        self.nodes = None
        self.nnum = None
        self.nnodes = None
    
        self.logger_handler = TopOptLogger("debug")

    def _make_nodes_array(self):
        # node numbering
        nnum = np.arange(self.nnodes)
        self.nnum = nnum.reshape((self.nnodes,1))

        # concat nodes array
        return np.concatenate((self.nnum,self.nodal_coords),axis=1)

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
        #z = np.zeros_like(x)
        self.nodal_coords = np.concatenate((x,y),axis=1) #,z

        self.nodes = self._make_nodes_array()

        # log
        logger.info("initialized mesh")

    def load_mesh(self):
        # TODO
        pass

class BoundaryCondition:
    def __init__(self,mesh,dofs):
        self.mesh = mesh
        self.dofs = dofs
        self.values = np.array([])
    
    # methods to add boundary conditions
    def add_by_point(self,coord,value,dof=None):
        
        # find node next to provided point
        distance = (self.mesh.nodal_coords[:,0] - coord[0])**2 + (self.mesh.nodal_coords[:,1] - coord[1])**2 #+ (self.mesh.nodal_coords[:,2] - coord[1])**2
        next_node = np.argmin(distance).reshape((1,1))

        if not np.isclose(min(distance),0):
            logger.info("no node at location {0}. applied to node at {1} instead. min distance {2}".format(coord,(self.mesh.nodal_coords[next_node,0],self.mesh.nodal_coords[next_node,1]),min(distance)))
        
        self.values = self._make_bc_array(next_node,value,dof,self.values)

    def add_by_component(self,name,value):
        pass

    def add_by_plane(self,normal,location,value,dof=None):
        # find nodes to apply bcs to
        if normal == [1,0]:
            distance = self.mesh.nodal_coords[:,0] - location
        elif normal == [0,1]:
            distance = self.mesh.nodal_coords[:,1] - location
        else:
            raise ValueError

        select_idx = np.isclose(distance,0)
        node_numbers = self.mesh.nnum[select_idx].reshape((-1,1))

        self.values = self._make_bc_array(node_numbers,value,dof,self.values)
        
    def _make_bc_array(self,node_numbers,value,dof,old_bc_array):
        number_of_bc_nodes = node_numbers.shape[0]
        # put nodes and given bc value in datastructure
        if isinstance(value,Iterable):
            nodal_values = np.tile(value,(number_of_bc_nodes,1))
            is_constraint = np.tile([True]*self.dofs,(number_of_bc_nodes,1))
        elif dof and self.dofs>=1:
            nodal_values = np.zeros((number_of_bc_nodes,self.dofs))
            nodal_values[:,dof-1] = value
            is_constraint = np.zeros((number_of_bc_nodes,self.dofs))
            is_constraint[:,dof-1] = True
        elif self.dofs>=1:
            nodal_values = np.zeros((number_of_bc_nodes,self.dofs))
            nodal_values[:,:] = value
            is_constraint = np.zeros((number_of_bc_nodes,self.dofs))
            is_constraint[:,:] = True
        else:
            raise ValueError("provide dof as optional argument")

        new_bcs = np.concatenate((node_numbers,nodal_values,is_constraint),axis=1)

        if not len(old_bc_array)==0:
            return np.concatenate((old_bc_array,new_bcs),axis=0)
        else:
            return new_bcs
    @property
    def nodes(self):
        return self.values[:,0].astype(int)
    def _pick_idx_constrained_dofs(self):
        return self.values[:,self.dofs+1:].astype(bool)
    def get_constrained_dofs(self,node_to_dof_map):
        dofs = node_to_dof_map[self.values[:,0].astype(int)] 
        idx_dofs = self._pick_idx_constrained_dofs()
        return dofs[idx_dofs]
    def get_constrained_values(self):
        values = self.values[:,1:self.dofs+1]
        idx_dofs = self._pick_idx_constrained_dofs()
        return values[idx_dofs]

class Displacement(BoundaryCondition):
    def __init__(self,mesh):
        _domain_specific_dofs = 2
        super().__init__(mesh,_domain_specific_dofs)

class Force(BoundaryCondition):
    def __init__(self,mesh):
        _domain_specific_dofs = 2
        super().__init__(mesh,_domain_specific_dofs)

class Temperature(BoundaryCondition):
    def __init__(self,mesh):
        _domain_specific_dofs = 1
        super().__init__(mesh,_domain_specific_dofs)

class Heat(BoundaryCondition):
    def __init__(self,mesh):
        _domain_specific_dofs = 2
        super().__init__(mesh,_domain_specific_dofs)