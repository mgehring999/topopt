from topopt.topopt import Mesh, OptimModel, PhysicalModel, Material
import pytest

def test_2x2():
    mesh = Mesh()
    mesh.rect_mesh(2)
    mesh.add_support()

    mat = Material()
    mat.set_structural_params(2.1e5,.3)

    pmodel = PhysicalModel(mesh,mat)
    
    optmizer = OptimModel(pmodel,0.5)