import solidspy.assemutil as ass
import numpy as np 

class PhysicalModel:
    """
    A `PhysicalModel` combines the information about the discretized domain, its physical properties, the boundaries and all system matrices. 

    Parameters
    ==========
    mesh: `Mesh` Provides information about the mesh
    mat: `Materal` Provides the material constants

    """
    def __init__(self,mesh,mat):
        self.mesh = mesh
        self.material = mat

        # set instance variables 
        self.Kglob = None
        self.Fglob = None
        self.x = [1]*self.mesh.nelem

        self.mats = self._set_materials()
        self.assemble_system()

        # TODO: implement loads
        self.Fglob = [1]*self.neq

    def _set_materials(self):
        mats = np.ones((self.mesh.nelem,2))
        mats[:,0] *= self.material.youngs
        mats[:,1] *= self.material.nu
        return mats

    def assemble_system(self):
        self.DME , self.IBC , self.neq = ass.DME(self.mesh.nodes, self.mesh.elements)
        self.Kglob = self.update_system_matrix()

    def update_system_matrix(self):

        # initialize empty system matrix as list of lists
        Kglob = [0]*self.neq
        for eq in range(self.neq):
            Kglob[eq] = [0]*self.neq

        # update system matrix with optimiziation variable
        for el in range(self.mesh.nelem):
            kloc,ndof,_=ass.retriever(self.mesh.elements,self.mats,self.mesh.nodes,el)
            kloc = kloc.tolist()
            dme = self.DME[el,:ndof]
            for row in range(ndof):
                glob_row=dme[row]
                if glob_row != -1:
                    for col in range(ndof):
                        glob_col = dme[col]
                        if glob_col != -1:
                            Kglob[glob_row][glob_col] = Kglob[glob_row][glob_col] +\
                                                        kloc[row][col]*self.x[el]
        return Kglob

class Material:
    def __init__(self):
        self.nu = None
        self.youngs = None
        
    def set_structural_params(self,youngs,nu):
        self.youngs = youngs
        self.nu = nu