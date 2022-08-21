#!/usr/bin/env python3

from topopt.physical import PhysicalModel, Material
from topopt.mesh import Mesh, Displacement, Force
from fem.fem import FEModel,StructuralElement

import matplotlib.pyplot as plt

ndiv = 20

mesh = Mesh()
mesh.rect_mesh(ndiv)

mat = Material()
mat.set_structural_params(2.1e5,0.3)

fem = FEModel(mesh,mat,StructuralElement)
fem.assemble()