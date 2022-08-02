import pytest

from topopt.physical import PhysicalModel, Material
from topopt.mesh import Mesh, Displacement, Force, Temperature

volfrac = 0.5
ndiv = 20

mesh = Mesh()
mesh.rect_mesh(ndiv)

temp1 = Temperature(mesh)
temp1.add_by_plane([1,0],-1,20)

temp2 = Temperature(mesh)
temp2.add_by_plane([1,0],1,100)

force = Force(mesh)
force.add_by_point((1,0),(0,-1))

support = Displacement(mesh)
support.add_by_plane([1,0],-1,0)

mat = Material()
mat.set_structural_params(2.1e5,.3)


def test_structural_model():
    bcs = [force,support]
    pmodel = PhysicalModel(mesh,mat,bcs)

def test_thermal_model():
    bcs = [temp1,temp2]
    pmodel = PhysicalModel(mesh,mat,bcs)

def test_thermo_mechanical_model():
    bcs = [temp1,temp2,force,support]
    pmodel = PhysicalModel(mesh,mat,bcs)

