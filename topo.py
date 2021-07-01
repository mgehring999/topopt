from model import model
from pyomo.opt import SolverFactory
import solidspy.postprocesor as pos
import matplotlib.pyplot as plt
import numpy as np

opt = SolverFactory("scip")
result = opt.solve(model,tee=True)
print(result)

x = np.array([model.x[i].value for i in model.elems])
u = np.array([model.u[i].value for i in model.eq])
print(x)

fig,axs = plt.subplots()
axs.imshow(x.reshape((2,2)))

plt.show()