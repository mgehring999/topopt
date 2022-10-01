from topopt.mesh import Mesh
from topopt.mesh import Temperature, Displacement, Force
from topopt.physical import Material, PhysicalModel

volfrac = 0.5
ndiv = 10

mesh = Mesh()
mesh.rect_mesh(ndiv)

temp_loads = Temperature(mesh)
temp_loads.add_by_plane([1,0],-1,100)
temp_loads.add_by_plane([1,0],1,20)

force_loads = Force(mesh)
force_loads.add_by_point([1,0],[0,-1])

support = Displacement(mesh)
support.add_by_plane([1,0],-1,0)

mat = Material()
mat.set_structural_params(2.1e5,0.3)

bcs = [support,temp_loads,force_loads]

pmodel = PhysicalModel(mesh,mat,bcs)