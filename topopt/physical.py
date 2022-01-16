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
    def __init__(self,mesh,mat,bcs,loads):
        self.mesh = mesh
        self.material = mat
        self.boundary_conditions = bcs
        self.loads = loads

        # set instance variables 
        self.Kglob = None
        self.Fglob = None
        self.x = [1]*self.mesh.nelem
        
        

        # log
        logger.info("started application of boundary conditions")

        # apply boundary conditions
        self.mesh = self.boundary_conditions.apply(self.mesh)
        
        # apply loads
        self.loads.apply(self.mesh)

        # assemble everything 
        logger.info("starting assembly of physical system")
        self.mats = self._set_materials()
        self.assemble_system()

    def _set_materials(self):
        mats = np.ones((self.mesh.nelem,2))
        mats[:,0] *= self.material.youngs
        mats[:,1] *= self.material.nu
        return mats

    def assemble_system(self):
        self.DME , self.IBC , self.neq = ass.DME(self.mesh.nodes, self.mesh.elements)
        self.Kglob = self.update_system_matrix()
        self.Fglob = ass.loadasem(self.loads.loads,self.IBC,self.neq)

    def update_system_matrix(self):
        # initialize empty system matrix as list of lists
        logger.debug("initialize stiffness matrix ...")
        Kglob = [ [0]*self.neq for _ in range(self.neq)]

        # update system matrix with optimiziation variable
        for el in range(self.mesh.nelem):
            kloc,ndof,_=ass.retriever(self.mesh.elements,self.mats,self.mesh.nodes,el)
            kloc = kloc.tolist()
            dme = self.DME[el,:ndof]
            
            for row in range(ndof):
                if dme[row] != -1:
                    for col in range(ndof):
                        if dme[col] != -1:
                            Kglob[dme[row]][dme[col]] += kloc[row][col]*self.x[el]

        return Kglob

class Material:
    def __init__(self):
        self.nu = None
        self.youngs = None
        
    def set_structural_params(self,youngs,nu):
        self.youngs = youngs
        self.nu = nu