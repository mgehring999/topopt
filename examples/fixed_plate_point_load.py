#!/usr/bin/env python3

from topopt.physical import PhysicalModel, Material
from topopt.mesh import Mesh, Displacement, Force
from topopt.topopt import StructuralOptim, Visualizer

import matplotlib.pyplot as plt

volfrac = 0.5
ndiv = 40

mesh = Mesh()
mesh.rect_mesh(ndiv)
print(mesh.nodal_coords)

force = Force(mesh)
force.add_by_point((1,0),-1,dof=2)

support = Displacement(mesh)
support.add_by_plane([1,0],-1,0)
print(support.values)

mat = Material()
mat.set_structural_params(2.1e5,.3)

bcs = [force,support]

pmodel = PhysicalModel(mesh,mat,bcs)

optimizer = StructuralOptim(pmodel,volfrac,5)
optimizer.run()

visu = Visualizer(pmodel)
visu.plot_result()
plt.show()


"""
loads = Load(mesh)
loads.add_by_coord((1,0),(0,-1))

support = Displacement(mesh)
support.add_by_edge("left","all")

mat = Material()
mat.set_structural_params(2.1e5,.3)

bcs = [loads,support]
pmodel = PhysicalModel(mesh,mat,bcs)

optimizer = StructuralOptim(pmodel,volfrac,5)
optimizer.run()

visu = Visualizer(pmodel)
visu.write_result()

visu.plot_result()
plt.show()

"""