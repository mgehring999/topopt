from topopt.topopt import StructuralOptim
from topopt.physical import PhysicalModel, Material
from topopt.mesh import Displacement, Load, Mesh
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

    loads = Load(mesh)
    loads.add_by_coord((1,1),(0,-1))

    support = Displacement()
    support.add_by_edge("left","all")

    mat = Material()
    mat.set_structural_params(2.1e5,.3)

    pmodel = PhysicalModel(mesh,mat,support,loads)
    
    optimizer = StructuralOptim(pmodel,volfrac,3)
    optimizer.run()

    assert optimizer.result