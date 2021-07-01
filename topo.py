from model import assemble_model
from pyomo.opt import SolverFactory
import solidspy.postprocesor as pos
import matplotlib.pyplot as plt
import numpy as np
import sys

# local packages
from mesh import rect_mesh

# input arguments
if len(sys.argv) < 3:
    sys.exit("no arguments passed, <volfrac> and <ndivisions> needed")
for argin in sys.argv[1:]:
    try:
        float(argin)
    except ValueError:
        sys.exit("arguments must be float or int")

volfrac = float(sys.argv[1])
ndiv = int(sys.argv[2])

# construct rectangular mesh
nodes,elements,loads = rect_mesh(ndiv)

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