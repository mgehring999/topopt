from topopt.mesh import Mesh
import numpy as np

def test_DomianWdith():
    mesh = Mesh()
    mesh.rect_mesh(5)
    xmin = min(mesh.nodal_coords[:,0])
    xmax = max(mesh.nodal_coords[:,0])
    
    assert np.isclose(xmax - xmin,2)

def test_DomianHeigth():
    mesh = Mesh()
    mesh.rect_mesh(5)
    ymin = min(mesh.nodal_coords[:,1])
    ymax = max(mesh.nodal_coords[:,1])
    
    assert np.isclose(ymax - ymin,2)

def test_NodesArrayDims():
    mesh = Mesh()
    mesh.rect_mesh(5)
    
    shape = mesh.nodes.shape

    assert shape[0] == 36 and shape[1] == 3