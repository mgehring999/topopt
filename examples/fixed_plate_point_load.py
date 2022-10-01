#!/usr/bin/env python3

from topopt.physical import PhysicalModel, Material
from topopt.mesh import Mesh, Displacement, Force
from topopt.topopt import StructuralOptim, Visualizer

import matplotlib.pyplot as plt

volfrac = 0.5
ndiv = 20

mesh = Mesh()
mesh.rect_mesh(ndiv)

force = Force(mesh)
force.add_by_point((1,0),-1,dof=2)

support = Displacement(mesh)
support.add_by_plane([1,0],-1,0)

mat = Material()
mat.set_structural_params(2.1e5,.3)

bcs = [force,support]

pmodel = PhysicalModel(mesh,mat,bcs)

optimizer = StructuralOptim(pmodel,volfrac,5)
optimizer.run()

visu = Visualizer(pmodel)
visu.plot_result()
plt.show()