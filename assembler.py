from mesh import rect_mesh
from solidspy import assemutil as ass
import numpy as np

def assemble(elements,mats,nodes,neq,DME,x,uel=None):
    Kglob = [0]*neq
    for eq in range(neq):
        Kglob[eq] = [0]*neq
    nels = len(elements)
    print(nels)
    for el in range(nels):
        kloc,ndof,_=ass.retriever(elements,mats,nodes,el,uel=uel)
        kloc = kloc.tolist()
        dme = DME[el,:ndof]
        for row in range(ndof):
            glob_row=dme[row]
            if glob_row != -1:
                for col in range(ndof):
                    glob_col = dme[col]
                    if glob_col != -1:
                        print(type(kloc[row][col]))
                        Kglob[glob_row][glob_col] = Kglob[glob_row][glob_col] +\
                                                    kloc[row][col]*x[el]

    return Kglob

if __name__ == "__main__":
    # load 2x2 mesh
    nodes,elements,loads = rect_mesh(2)
    nelem = elements.shape[0]
    nnodes = nodes.shape[0]

    # assembly operator
    DME , IBC , neq = ass.DME(nodes, elements)

    # material
    mats = np.ones((nelem,2))
    mats[:,1] *= 0.3

    # optimiziation variable
    x = np.ones((nelem,)).tolist()

    # assemble
    Kglob = assemble(elements,mats,nodes,neq,DME,x)
    Kglob_ref = ass.assembler(elements,mats,nodes,neq,DME,sparse=False)



