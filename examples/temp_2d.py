#!/usr/bin/env python3

from topopt.topopt import StructuralOptim, Visualizer
from topopt.physical import PhysicalModel, Material
from topopt.mesh import Displacement, Heat, Mesh, Temperature
import matplotlib.pyplot as plt


volfrac = 0.5
ndiv = 10

mesh = Mesh()
mesh.rect_mesh(ndiv)

temp_loads = Temperature(mesh)
temp_loads.add_by_plane([1,0],-1,100)
temp_loads.add_by_plane([1,0],1,20)

mat = Material()
mat.set_structural_params(2.1e5,.3) # thermal conductivity, specific heat, density

bcs = [temp_loads]
pmodel = PhysicalModel(mesh,mat,bcs)
pmodel = PhysicalModel(mesh,temp_loads,mat)