from topopt import Mesh
import numpy as np
import pytest

@pytest.mark.parametrize("ndiv, coord, value, dof, result",[
    (2,"x",-1,1,[[ 0, -1, -1,  -1,  0],
                [ 1,  0, -1,  0,  0],
                [ 2,  1, -1,  0,  0],
                [ 3, -1,  0,  -1,  0],
                [ 4,  0,  0,  0,  0],
                [ 5,  1,  0,  0,  0],
                [ 6, -1,  1,  -1,  0],
                [ 7,  0,  1,  0,  0],
                [ 8,  1,  1,  0,  0]]),

    (2,"x",-1,2,[[ 0, -1, -1,  0,  -1],
                [ 1,  0, -1,  0,  0],
                [ 2,  1, -1,  0,  0],
                [ 3, -1,  0,  0,  -1],
                [ 4,  0,  0,  0,  0],
                [ 5,  1,  0,  0,  0],
                [ 6, -1,  1,  0,  -1],
                [ 7,  0,  1,  0,  0],
                [ 8,  1,  1,  0,  0]]),
])
def test_bc(ndiv,coord,value,dof,result):
    mesh = Mesh()
    mesh.rect_mesh(ndiv)
    test = mesh._set_bc(mesh.nodes,coord,value,dof)
    assert np.allclose(test,result)

@pytest.mark.parametrize("edge, dof, result",[
                ("left","all",
                [[ 0, -1, -1,  -1,  -1],
                [ 1,  0, -1,  0,  0],
                [ 2,  1, -1,  0,  0],
                [ 3, -1,  0,  -1,  -1],
                [ 4,  0,  0,  0,  0],
                [ 5,  1,  0,  0,  0],
                [ 6, -1,  1,  -1,  -1],
                [ 7,  0,  1,  0,  0],
                [ 8,  1,  1,  0,  0]]),
                
                ("left",1,
                [[ 0, -1, -1, -1,  0],
                [ 1,  0, -1,  0,  0],
                [ 2,  1, -1,  0,  0],
                [ 3, -1,  0,  -1,  0],
                [ 4,  0,  0,  0,  0],
                [ 5,  1,  0,  0,  0],
                [ 6, -1,  1,  -1,  0],
                [ 7,  0,  1,  0,  0],
                [ 8,  1,  1,  0,  0]]),
                
                ("left",2,
                [[ 0, -1, -1, 0,  -1],
                [ 1,  0, -1,  0,  0],
                [ 2,  1, -1,  0,  0],
                [ 3, -1,  0,  0,  -1],
                [ 4,  0,  0,  0,  0],
                [ 5,  1,  0,  0,  0],
                [ 6, -1,  1,  0,  -1],
                [ 7,  0,  1,  0,  0],
                [ 8,  1,  1,  0,  0]]),
                
                ("right","all",
                [[ 0, -1, -1, 0,  0],
                [ 1,  0, -1,  0,  0],
                [ 2,  1, -1, -1, -1],
                [ 3, -1,  0,  0,  0],
                [ 4,  0,  0,  0,  0],
                [ 5,  1,  0, -1, -1],
                [ 6, -1,  1,  0,  0],
                [ 7,  0,  1,  0,  0],
                [ 8,  1,  1, -1, -1]]),

                ("top","all",
                [[ 0, -1, -1, 0,  0],
                [ 1,  0, -1,  0,  0],
                [ 2,  1, -1,  0,  0],
                [ 3, -1,  0,  0,  0],
                [ 4,  0,  0,  0,  0],
                [ 5,  1,  0,  0,  0],
                [ 6, -1,  1, -1, -1],
                [ 7,  0,  1, -1, -1],
                [ 8,  1,  1, -1, -1]]),

                ("bot","all",
                [[ 0, -1, -1,-1, -1],
                [ 1,  0, -1, -1, -1],
                [ 2,  1, -1, -1, -1],
                [ 3, -1,  0,  0,  0],
                [ 4,  0,  0,  0,  0],
                [ 5,  1,  0,  0,  0],
                [ 6, -1,  1,  0,  0],
                [ 7,  0,  1,  0,  0],
                [ 8,  1,  1,  0,  0]])

])
def test_support(edge,dof,result):
    mesh=Mesh()
    mesh.rect_mesh(2)
    mesh.add_support(edge=edge,dof=dof)
    test = mesh.nodes
    assert np.allclose(test,result)