from solidspy.preprocesor import rect_grid
import numpy as np

class Mesh:
    """Mesh creates or loads and saves the information about the meshed domain
    """
    def __init__(self):
        self.x = None
        self.y = None
        self.elements = None
    
    def rect_mesh(self,ndiv):
        x,y,self.elements = rect_grid(2,2,ndiv,ndiv)

    	nnodes = len(x)
        self.x = x.reshape((nnodes,1))
        self.y = y.reshape((nnodes,1))

    def _construct_nodes_array(self):

        # node numbering
        nnum = np.arange(nnodes)
        nnum = nnum.reshape((nnodes,1))

        # init bc array, no bcs initially
        bc_init = np.zeros((nnodes,2))

        # concat to SolidsPy format
        nodes = np.concatenate((nnum,x,y,bc_init),axis=1)

        pass

    def load_mesh(self):
        pass
    



class StructuralModel:
    """this holds all data necessary for the finite element model"""
    def __init__(self,mesh):
        self.mesh = mesh
        self.Kglob = None
        self.Fglob = None
    



class OptimModel:
    """idk"""

    def __init__(self):
        pass
        #self.FEmodel = 

if __name__ == "__main__":
    mesh = Mesh()
    mesh.rect_mesh(2)

    FEModel = StructuralModel(mesh)