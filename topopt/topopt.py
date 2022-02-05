import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

from pyomo.environ import *

import logging
logger = logging.getLogger('topopt')

class OptimModel:
    def __init__(self,physical_model,volfrac):
        self.pmodel = physical_model
        self.vol = volfrac*sum(self.pmodel.x)
        self.result = None

        logger.info("started assembly of optimization model")

        self.model = ConcreteModel(name="topo")

        self.model.elems = Set(initialize=self.pmodel.mesh.elements[:,0],domain=NonNegativeIntegers)
        self.model.eq = Set(initialize=range(self.pmodel.neq),domain=NonNegativeIntegers)

    def run(self):
        self.solver = SolverFactory("ipopt")
        logger.info("started solution process")
        self.result = self.solver.solve(self.model,tee=True)
        logger.info("finished solution process")

        # write the optimiziation variable back to the physical model
        self.pmodel.x = np.array([getattr(self.model,"x{}".format(i)).value for i in self.model.elems])

class StructuralOptim(OptimModel):
    def __init__(self, physical_model, volfrac, penal):
        super().__init__(physical_model, volfrac)
        
        # solve for the deformations to initialize optimiziation variables
        self.u_init = physical_model.solve_system_eq()
        
        # build pyomo model
        self.penal = penal
        self._init_system_matrices()
        self._init_optim_vars()
        logger.info("initialized system matrices")
        self._make_constraints()
        logger.info("made constraint functions")
        self._make_objective()
        logger.info("made objective functions")

    def _init_system_matrices(self):
        self.model.K = Param(initialize=self.pmodel.Kglob,default=0,mutable=True,within=Any) 
        self.model.F = Param(self.model.eq,initialize=dict(enumerate(self.pmodel.Fglob)),within=Any)

    def _init_optim_vars(self):
        # add deformation as optim vars to 
        for dof,u_init in zip(range(self.pmodel.neq),self.u_init):
            setattr(self.model,"u"+str(dof),Var(initialize=u_init))

        for elem in range(self.pmodel.mesh.nelem):
            setattr(self.model,"x"+str(elem),Var(bounds=(1e-5,1),initialize=1))

        # update material definition with pyomo reference for new stiffness matrix
        self.pmodel.x = [getattr(self.model,"x"+str(elem))**self.penal for elem in self.model.elems]

        # update stiffness matrix and assign to model
        self.model.K = self.pmodel.update_system_matrix()

    def _make_constraints(self):
        logger.info("making volume constraint ...")
        self.model.vol_con = Constraint(rule=self._vol_rule)

        logger.info("making finite element equation constraint ...")
        self.model.FKU_con = Constraint(self.model.eq, rule=self._FKU_rule)
        
    def _make_objective(self):
        self.model.obj = Objective(expr=self._comp_rule(self.model),sense=minimize)
        
    def _comp_rule(self,m):
        return sum([sum([m.K.value[row][col]*getattr(m,"u{}".format(col)) for col in m.eq])*getattr(m,"u{}".format(row)) for row in m.eq])

    # volume fraction
    def _vol_rule(self,m):
        return sum([getattr(m,"x"+str(j)) for j in m.elems]) == self.vol 

    def _FKU_rule(self,m, i):
        # logging
        if i % 100 == 0:
            logger.debug("making CE for DOF number {}".format(i))
        return sum([m.K.value[i][j]*getattr(m,"u{}".format(j)) for j in m.eq]) == m.F[i]

class Visualizer:
    def __init__(self,pmodel):
        self.model = pmodel
        self.result = None
        
        # default filename for read/write of optimization results
        self.filename = "optim.rst"

        # plot controls
        colors = ["white", "grey","grey","blue"]
        nodes = [0.0, 0.4, 0.6,1.0]
        self.cmap = LinearSegmentedColormap.from_list("mycmap", list(zip(nodes, colors)))

    def _load_result(self):
        ndiv = int(np.sqrt(len(self.model.x)))
        self.result = self.model.x.reshape((ndiv,ndiv))

    def write_result(self,filename=None):
        """
        writes results from PhysicalModel to disk
        """
        # overwrite default filename if custom filename is provided
        if filename:
            self.filename = filename+".rst"

        # write to text file
        np.savetxt(self.filename,self.model.x)

    def read_result(self,filename=None):
        """
        reads results from text file to self.result variable
        """

        if filename:
            self.filename = filename
        
        self.result = np.loadtxt(self.filename)
        print(self.result.shape)

    def show_result(self):
        self._load_result()
        self.fig = plt.figure()
        self.axs = plt.gca()
        self.axs.imshow(self.result,cmap=self.cmap,origin="lower")
        plt.show()
