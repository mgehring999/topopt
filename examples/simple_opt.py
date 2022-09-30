#!/usr/bin/env python3

from topopt.physical import Material
from topopt.mesh import Mesh, Displacement, Force
from fem.fem import FEModel,StructuralElement

import numpy as np
import scipy.optimize as opt

ndiv = 20

mesh = Mesh()
mesh.rect_mesh(ndiv)

mat = Material()
mat.set_structural_params(2.1e5,0.3)

support = Displacement(mesh)
support.add_by_plane([1,0],-1,0)

load = Force(mesh)
load.add_by_point((1,0),(0,-1))

fem = FEModel(mesh,mat,StructuralElement)

def obj_func(fem,support,load):
    mat = fem.K
    vec = fem.solve(support,load)
    return np.dot(vec.T,np.matmul(mat,vec))

u = fem.solve(support,load)
K = fem.K

objective = obj_func(K,u)
print(objective)

res = opt.minimize(obj_func,x0,method="Nelder-Mead")
print(res.x)

