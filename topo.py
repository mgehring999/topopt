from model import assemble_model
from pyomo.opt import SolverFactory
import solidspy.postprocesor as pos
import matplotlib.pyplot as plt
import numpy as np
import sys

# local packages
from mesh import add_point_force, rect_mesh

# input arguments
if len(sys.argv) < 5:
    sys.exit("no arguments passed, <volfrac> <ndivisions> <load_coordinate> <load_direction> needed")

# plausibility
for argin in sys.argv[1:3]:
    try:
        float(argin)
    except ValueError:
        sys.exit("arguments must be float or int")

for argin in sys.argv[3:5]:
    try:
        assert "," in argin
    except AssertionError:
        sys.exit("pass load coordinates and directions in the format x,y" )

# assign to variables
volfrac = float(sys.argv[1])
ndiv = int(sys.argv[2])
load_coords = tuple([float(comp) for comp in sys.argv[3].split(",")])
load_direction = tuple([float(comp) for comp in sys.argv[4].split(",")])

# construct rectangular mesh
nodes,elements = rect_mesh(ndiv)

# add loads to mesh
loads = add_point_force(nodes,load_coords,load_direction)

# assemble model
model = assemble_model(nodes,elements,loads,volfrac)

# solve model with scip
opt = SolverFactory("scip")
result = opt.solve(model,tee=True)

# post process results
x = np.array([model.x[i].value for i in model.elems])
u = np.array([model.u[i].value for i in model.eq])
print(x)

# plot results
fig,axs = plt.subplots()
axs.imshow(x.reshape((ndiv,ndiv)))

plt.show()