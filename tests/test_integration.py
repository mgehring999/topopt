from topopt.topopt import StructuralOptim
from topopt.physical import PhysicalModel, Material
from topopt.mesh import Displacement, Mesh, Force
import pytest

@pytest.mark.parametrize("ndiv, volfrac",[
    (2,0.5),
    (3,0.5),
    (5,1),
    (5,0.5),
    (8,0.5),
])
def test_parameters(ndiv,volfrac):
    mesh = Mesh()
    mesh.rect_mesh(ndiv)

    loads = Force(mesh)
    loads.add_by_point((1,1),(0,-1))

    support = Displacement(mesh)
    support.add_by_plane([1,0],-1,0)

    mat = Material()
    mat.set_structural_params(2.1e5,.3)

    pmodel = PhysicalModel(mesh,mat,[loads,support])
    
    optimizer = StructuralOptim(pmodel,volfrac,3)
    optimizer.run()

    assert optimizer.result