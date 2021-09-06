from topopt.topopt import Mesh, PhysicalModel, Material, StructuralOptim
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
    mesh.add_support()

    mat = Material()
    mat.set_structural_params(2.1e5,.3)

    pmodel = PhysicalModel(mesh,mat)
    
    optimizer = StructuralOptim(pmodel,volfrac,3)
    optimizer.run()

    assert optimizer.result