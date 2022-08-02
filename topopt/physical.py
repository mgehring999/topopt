from topopt.mesh import Displacement, Force, Temperature, Heat

import solidspy.assemutil as ass
import numpy as np 
import logging

logger = logging.getLogger('topopt')

class PhysicalModel:
    """
    A `PhysicalModel` combines the information about the discretized domain, its physical properties, the boundaries and all system matrices. 

    Parameters
    ==========
    mesh: `Mesh` Provides information about the mesh
    mat: `Materal` Provides the material constants

    """
    def __init__(self,mesh,mat,bcs):
        self.mesh = mesh
        self.material = mat
        self.bcs = bcs
        #self.loads = loads

        # no initial constraints
        self.constrained_dofs = []

        # set instance variables 
        self.x = [1]*self.mesh.nelem

        # assemble system matrices needed
        load_in_bcs = any([isinstance(bc,Force) for bc in self.bcs])
        displacement_in_bcs = any([isinstance(bc,Displacement) for bc in self.bcs])
        if load_in_bcs or displacement_in_bcs:
            # initiallize load and displacement vector
            self.Kglob = None
            self.Fglob = None
            self.Uglob = None

            # set material paramters
            self.mats = self._set_materials()

            # assemble stiffness matrix
            nodes = np.concatenate((self.mesh.nodes,np.zeros((self.mesh.nnodes,2))),axis=1)
            self.element_to_dof_map , self.node_to_dof_map , self.ndof = ass.DME(nodes, self.mesh.elements)
            self.Kglob = np.array(self.update_system_matrix())
            print("Dims, first assem",self.Kglob.shape)

        # apply boundary conditions and assemble load vectors
        for bc in self.bcs:
            if isinstance(bc,Force):
                # assemble load 
                if self.Fglob:
                    self.Fglob += ass.loadasem(bc.values[:,0:3],self.node_to_dof_map,self.ndof)
                else:
                    self.Fglob = ass.loadasem(bc.values[:,0:3],self.node_to_dof_map,self.ndof)

            if isinstance(bc,Displacement):
                # write DOF values to solution vector
                self.Uglob = self.solve_system_eq()
                #self.Uglob += bc.values

        # constraints
        for bc in self.bcs:
            if isinstance(bc,Displacement):
                bc_nodes = bc.values[:,0].astype(int) # nodes with boundary conditions
                bc_node_dofs = self.node_to_dof_map[bc_nodes,:] # dofs for each node with bc
                idx_constrained_dofs = bc.values[:,bc.dofs+1:].astype(bool)
                
                self.constrained_dofs = bc_node_dofs[idx_constrained_dofs]
                u_constrained = bc.values[:,1:bc.dofs+1]
                self.U_constrained = u_constrained[bc.values[:,bc.dofs+1:].astype(bool)]
                print(self.constrained_dofs)
                # apply constraints by deleting equations 
                self.F_constrained = self.Fglob[self.constrained_dofs]
                self.K_constrained = self.Kglob[np.ix_(self.constrained_dofs,self.constrained_dofs)]
                self.Fglob = np.delete(self.Fglob,self.constrained_dofs)
                self.Uglob = np.delete(self.Uglob,self.constrained_dofs)
                self.Kglob = np.delete(self.Kglob,self.constrained_dofs,axis=0)
                self.Kglob = np.delete(self.Kglob,self.constrained_dofs,axis=1)
                print("Dims, after bc application",self.Kglob.shape)
                
                self.neq = len(self.Uglob)

        # log
        logger.info("started application of boundary conditions")

    def _set_materials(self):
        mats = np.ones((self.mesh.nelem,2))
        mats[:,0] *= self.material.youngs
        mats[:,1] *= self.material.nu
        return mats

    def update_system_matrix(self):
        # initialize empty system matrix as list of lists
        logger.debug("initialize stiffness matrix ...")
        Kglob = [ [0]*self.ndof for _ in range(self.ndof)]

        # update system matrix with optimiziation variable
        for el in range(self.mesh.nelem):
            kloc,ndof,_=ass.retriever(self.mesh.elements,self.mats,self.mesh.nodes,el)
            kloc = kloc.tolist()
            dme = self.element_to_dof_map[el,:ndof]
            #print(el,dme)
            for row in range(ndof):
                if not dme[row] in self.constrained_dofs:
                    for col in range(ndof):
                        if not dme[col] in self.constrained_dofs:
                            #print(dme[row],dme[col])
                            Kglob[dme[row]][dme[col]] += kloc[row][col]*self.x[el]

        return Kglob

    def solve_system_eq(self):
        return np.linalg.solve(self.Kglob,self.Fglob)

class Material:
    def __init__(self):
        self.nu = None
        self.youngs = None
        
    def set_structural_params(self,youngs,nu):
        self.youngs = youngs
        self.nu = nu
