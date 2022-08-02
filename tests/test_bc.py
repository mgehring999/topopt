from topopt.mesh import Displacement, Force, Mesh, Temperature
from topopt.physical import PhysicalModel
import numpy as np
import pytest
ndiv = 2

mesh = Mesh()
mesh.rect_mesh(ndiv)

@pytest.mark.parametrize("dof_values, dof, result",[
    ([1,0],None,
    np.array([[0,1,0,1,1],
              [3,1,0,1,1],
              [6,1,0,1,1]])),
    ([0,1],None,
    np.array([[0,0,1,1,1],
              [3,0,1,1,1],
              [6,0,1,1,1]])),
    (0,None,
    np.array([[0,0,0,1,1],
              [3,0,0,1,1],
              [6,0,0,1,1]])),
    (1,None,
    np.array([[0,1,1,1,1],
              [3,1,1,1,1],
              [6,1,1,1,1]])),
    (0,1,    
    np.array([[0,0,0,1,0],
              [3,0,0,1,0],
              [6,0,0,1,0]])),
    (0,2,    
    np.array([[0,0,0,0,1],
              [3,0,0,0,1],
              [6,0,0,0,1]])),    
    (1,1,    
    np.array([[0,1,0,1,0],
              [3,1,0,1,0],
              [6,1,0,1,0]])),
    (1,2,    
    np.array([[0,0,1,0,1],
              [3,0,1,0,1],
              [6,0,1,0,1]])),
])
def test_support_by_plane(dof_values, dof,result):
    support = Displacement(mesh)
    support.add_by_plane([1,0],-1,dof_values,dof=dof)
    
    assert np.allclose(support.values,result)

@pytest.mark.parametrize("coords, load, dof, result",[
    ([1,0],[0,-1],None,np.array([5,0,-1,1,1])),
    ([1,1],[0,-1],None,np.array([8,0,-1,1,1])),
    ([0,1],[0,-1],None,np.array([7,0,-1,1,1])),
    ([0,1],1,1,np.array([7,1,0,1,0])),
    ([0,1],0,2,np.array([7,0,0,0,1])),
])
def test_force_by_point(coords,load,dof,result):
    fload = Force(mesh)
    fload.add_by_point(coords,load,dof=dof)
    assert np.allclose(fload.values,result)

@pytest.mark.parametrize("dof_values, dof, result",[
    ([1,0],None,
    np.array([[0,1,0,1,1],
              [3,1,0,1,1],
              [6,1,0,1,1]])),
    ([0,1],None,
    np.array([[0,0,1,1,1],
              [3,0,1,1,1],
              [6,0,1,1,1]])),
    (0,None,
    np.array([[0,0,0,1,1],
              [3,0,0,1,1],
              [6,0,0,1,1]])),
    (1,None,
    np.array([[0,1,1,1,1],
              [3,1,1,1,1],
              [6,1,1,1,1]])),
    (0,1,    
    np.array([[0,0,0,1,0],
              [3,0,0,1,0],
              [6,0,0,1,0]])),
    (0,2,    
    np.array([[0,0,0,0,1],
              [3,0,0,0,1],
              [6,0,0,0,1]])),    
    (1,1,    
    np.array([[0,1,0,1,0],
              [3,1,0,1,0],
              [6,1,0,1,0]])),
    (1,2,    
    np.array([[0,0,1,0,1],
              [3,0,1,0,1],
              [6,0,1,0,1]])),
])
def test_force_by_plane(dof_values, dof,result):
    fload = Force(mesh)
    fload.add_by_plane([1,0],-1,dof_values,dof=dof)
    assert np.allclose(fload.values,result)

def test_add_multiple_bcs():
    fload = Force(mesh)
    fload.add_by_plane([1,0],-1,5)
    fload.add_by_plane([1,0],1,15)

    assert np.allclose(fload.values,np.array([
              [0,5,5,1,1],
              [3,5,5,1,1],
              [6,5,5,1,1],
              [2,15,15,1,1],
              [5,15,15,1,1],
              [8,15,15,1,1]]))

"""
def test_support_by_point():
    support = Displacement(mesh)
    support.add_by_point([1,0],-1)

def test_force_by_plane():
    fload = Force(mesh)
    fload.add_by_plane([1,0],-1,[0,-1])

#test_support_by_plane()
#test_support_by_point()

print("make displacement constraint")
#support.add_by_plane([1,0],-1,5,dof=1)
#support.add_by_plane([1,0],-1,5,dof=2)
print(support.values)

print("make force load")
#fload.add_by_point((1,0),[0,10]) #231
print(fload.values)

print("make displacement load")
dload = Displacement(mesh)
dload.add_by_point((1,0),[0,10]) #231
print(dload.values)

inlet_temp = Temperature(mesh)
inlet_temp.add_by_point((1,1),25)
# print(support.values)
# print(load.values)

bcs= [support, load]

# pmodel= PhysicalModel(mesh,mats,bcs)

"""
