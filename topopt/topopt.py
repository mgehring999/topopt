from solidspy.preprocesor import rect_grid
import solidspy.assemutil as ass

import numpy as np
import scipy as sp

from pyomo.environ import *

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

class Material:
    def __init__(self):
        self.nu = None
        self.youngs = None
        
    def set_structural_params(self,youngs,nu):
        self.youngs = youngs
        self.nu = nu

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

class OptimModel:
    def __init__(self,physical_model,volfrac):
        self.pmodel = physical_model
        self.vol = volfrac*sum(self.pmodel.x)

        self.model = ConcreteModel(name="topo")

        self.model.elems = Set(initialize=self.pmodel.mesh.elements[:,0],domain=NonNegativeIntegers)
        self.model.eq = Set(initialize=range(self.pmodel.neq),domain=NonNegativeIntegers)


class StructuralOptim(OptimModel):
    def __init__(self, physical_model, volfrac, penal):
        super().__init__(physical_model, volfrac)
        self.penal = penal
        self._init_system_matrices()
        self._make_constraints()
        self._make_objective()

        self.model.pprint()

    def _init_system_matrices(self):
        self.model.x = Var(self.model.elems,bounds=(1e-5,1),initialize=1)
        self.model.K = Param(initialize=self.pmodel.Kglob,default=0,mutable=True,within=Any) 
        self.model.F = Param(self.model.eq,initialize=dict(enumerate(self.pmodel.Fglob)),within=Any)
        self.model.u = Var(self.model.eq,initialize=dict(enumerate(sp.sparse.linalg.spsolve(self.pmodel.Kglob,self.pmodel.Fglob))))

    def _make_constraints(self):
        self.model.FKU_con = Constraint(self.model.eq, rule=self._FKU_rule)
        self.model.vol_con = Constraint(rule=self._vol_rule)
        
    def _make_objective(self):
        self.model.obj = Objective(expr=self._comp_rule(self.model),sense=minimize)

    def _comp_rule(self,m):
        # update material definition
        self.pmodel.x = [m.x[elem]**self.penal for elem in m.elems]

        # update stiffness matrix and assign to model
        m.K = self.pmodel.update_system_matrix()
        return sum([sum([m.K.value[row][col]*m.u[col] for col in m.eq])*m.u[row] for row in m.eq])

    # volume fraction
    def _vol_rule(self,m):
        return sum([m.x[j] for j in m.elems]) == self.vol 

    def _FKU_rule(self,m, i):

        # update material definition
        self.pmodel.x = [m.x[elem]**self.penal for elem in m.elems]

        # update stiffness matrix and assign to model
        m.K = self.pmodel.update_system_matrix()

        return sum([m.K.value[i][j]*m.u[j] for j in m.eq]) == m.F[i]
