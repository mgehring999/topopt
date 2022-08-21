from solidspy.assemutil import retriever
from solidspy.preprocesor import rect_grid

from topopt.physical import Material
from fem.fem import FEModel, StructuralElement
from topopt.mesh import Mesh
import numpy as np

x,y,element = rect_grid(1,-1,1,1)
nnodes = len(x)
nodes = np.concatenate((np.arange(nnodes).reshape((nnodes,1)),x.reshape((nnodes,1)),y.reshape((nnodes,1))),axis=1)
mats = np.array([[2.1e5,0.3]])

mesh = Mesh()
mesh.rect_mesh(10)

mat = Material()
mat.set_structural_params(2.1e5,0.3)

fem = FEModel(mesh,mat,StructuralElement)

def structural_quad4():
    kloc,_,_ = retriever(element,mats,nodes,0)
