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
        if edge == "left":
            if dof == "all":
                dof = [1,2]
            elif dof is not int or list:
                ValueError("DOFs are specified as int or list of ints")

            for i in dof:
                print(i)
                self.nodes = self._set_bc(self.nodes,"x",-1,i)
        
        elif edge == "right":
            if dof == "all":
                dof = [1,2]
            elif dof is not int or list:
                ValueError("DOFs are specified as int or list of ints")

            for i in dof:
                self.nodes = self._set_bc(self.nodes,"x",1,i)
            
        elif edge == "top":
            if dof == "all":
                dof = [1,2]
            elif dof is not int or list:
                ValueError("DOFs are specified as int or list of ints")

            for i in dof:
                self.nodes = self._set_bc(self.nodes,"y",1,i)
            
        elif edge == "bottom" or edge=="bot":
            if dof == "all":
                dof = [1,2]
            elif dof is not int or list:
                ValueError("DOFs are specified as int or list of ints")

            for i in dof:
                self.nodes = self._set_bc(self.nodes,"y",-1,i)

    def load_mesh(self):
        # TODO
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
    mesh.nodes = mesh._set_bc(mesh.nodes,"x",-1,1)
    mesh.nodes = mesh._set_bc(mesh.nodes,"x",-1,2)
    print(mesh.nodes)

    mesh = Mesh()
    mesh.rect_mesh(2)
    mesh.add_support()
    print(mesh.nodes)