#!/usr/bin/env python3

import numpy as np

from topopt.physical import PhysicalModel, Material 
from topopt.mesh import Mesh, Temperature

import matplotlib.pyplot as plt

volfrac = 0.5
ndiv = 2

mesh = Mesh()
mesh.rect_mesh(ndiv)

temp_loads = Temperature(mesh)
temp_loads.add_by_plane([1,0],-1,100)
temp_loads.add_by_plane([1,0],1,20)

mat = Material()
mat.set_thermal_params(0.05) # conductivity W/mmK

bcs = [temp_loads]
pmodel = PhysicalModel(mesh,mat,bcs)
