#!/usr/bin/env python3

from topopt.physical import PhysicalModel, Material
from topopt.mesh import Mesh, Displacement, Force
from topopt.topopt import Visualizer
from fem.fem import FEModel,StructuralElement

import matplotlib.pyplot as plt

ndiv = 6

mesh = Mesh()
mesh.rect_mesh(ndiv)

mat = Material()
mat.set_structural_params(2.1e5,0.3)

support = Displacement(mesh)
support.add_by_plane([1,0],-1,0)

load = Force(mesh)
load.add_by_point((1,0),(0,-1))

fem = FEModel(mesh,mat,StructuralElement)
u = fem.solve(support,load=load)

print(u)