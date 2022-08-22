import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
import numpy as np

def plot_dof(fem,field,dof):

    triangs = []
    for elem in fem.elements:
        triangs.append(elem[[3, 4, 5]])
        triangs.append(elem[[5, 6, 3]])

    triangs = np.array(triangs)

    tri = Triangulation(fem.nodes[:,0], fem.nodes[:,1],triangs)
    nnodes,ndofs = fem.nodes.shape[0],len(fem.node_to_dof_map[0])
    u_node = np.zeros((nnodes,ndofs))
    for node,dofs in enumerate(fem.node_to_dof_map):
        u_node[node] = field[dofs]
        
    _,axs = plt.subplots()
    tcf=axs.tricontourf(tri,u_node[:,dof])
    plt.colorbar(tcf)