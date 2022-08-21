#!/usr/bin/env python3

from topopt.physical import PhysicalModel, Material
from topopt.mesh import Mesh, Temperature
from fem.fem import FEModel,ThermalElement

import matplotlib.pyplot as plt

ndiv = 20

mesh = Mesh()
mesh.rect_mesh(ndiv)

mat = Material()
mat.set_thermal_params(0.05)

fem = FEModel(mesh,mat,ThermalElement)
fem.assemble()