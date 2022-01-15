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
        self.pmodel.x = np.array([self.model.x[i].value for i in self.model.elems])

class StructuralOptim(OptimModel):
    def __init__(self, physical_model, volfrac, penal):
        super().__init__(physical_model, volfrac)
        self.penal = penal
        self._init_system_matrices()
        logger.info("initialized system matrices")
        self._make_constraints()
        logger.info("made constraint functions")
        self._make_objective()
        logger.info("made objective functions")

    def _init_system_matrices(self):
        self.model.x = Var(self.model.elems,bounds=(1e-5,1),initialize=1)
        self.model.K = Param(initialize=self.pmodel.Kglob,default=0,mutable=True,within=Any) 
        self.model.F = Param(self.model.eq,initialize=dict(enumerate(self.pmodel.Fglob)),within=Any)
        self.model.u = Var(self.model.eq,initialize=dict(enumerate(sp.sparse.linalg.spsolve(self.pmodel.Kglob,self.pmodel.Fglob))))

    def _make_constraints(self):
        self.model.vol_con = Constraint(rule=self._vol_rule)
        logger.info("made volume rule")
        self.model.FKU_con = Constraint(self.model.eq, rule=self._FKU_rule)
        
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
        # logging
        logger.debug("making CE for DOF number {}".format(i))

        # update material definition
        self.pmodel.x = [m.x[elem]**self.penal for elem in m.elems]
        #logger.debug("updated x for DOF number {}".format(i))

        # update stiffness matrix and assign to model
        m.K = self.pmodel.update_system_matrix()
        #logger.debug("updated K for DOF number {}".format(i))

        return sum([m.K.value[i][j]*m.u[j] for j in m.eq]) == m.F[i]

class Visualizer:
    def __init__(self,pmodel):
        self.model = pmodel
        self.result = None
        colors = ["white", "grey","grey","blue"]
        nodes = [0.0, 0.4, 0.6,1.0]
        self.cmap = LinearSegmentedColormap.from_list("mycmap", list(zip(nodes, colors)))

    def _load_result(self):
        ndiv = int(np.sqrt(len(self.model.x)))
        self.result = self.model.x.reshape((ndiv,ndiv))

    def show_result(self):
        self._load_result()
        self.fig = plt.figure()
        self.axs = plt.gca()
        self.axs.imshow(self.result,cmap=self.cmap,origin="lower")
        plt.show()
