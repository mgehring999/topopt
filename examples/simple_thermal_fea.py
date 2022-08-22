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

temp = Temperature(mesh)
temp.add_by_plane([1,0],-1,100)
temp.add_by_plane([1,0],1,20)

fem = FEModel(mesh,mat,ThermalElement)
t = fem.solve(temp)
print(t)