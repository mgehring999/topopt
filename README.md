# ``topopt`` - Topology Optimization in Python 
This is a python implementation of a two-dimensional structural optimizitation framework started during summer term 2021 as a class project and continued since then.

## Installation
Install python requirements with pip
```
pip install -r requirements.txt
```
Download the IPOPT solver from [here](https://ampl.com/products/solvers/open-source/)
Linux Users place it in the /bin directory and on Windows you put the path to the binary in the PATH system variable

## Example script
This script models a quadratic plate clamped on the left side and loaded with a point force in the middle of the right edge. The plate is discretized with 400 linear plate elements. The optimization goal is to delete 50% of the volume and minimize the compliance. 

```python
    from topopt.physical import PhysicalModel, Material
    from topopt.mesh import Mesh, Displacement, Force
    from topopt.topopt import StructuralOptim, Visualizer

    import matplotlib.pyplot as plt

    volfrac = 0.5
    ndiv = 20

    mesh = Mesh()
    mesh.rect_mesh(ndiv)

    force = Force(mesh)
    force.add_by_point((1,0),-1,dof=2)

    support = Displacement(mesh)
    support.add_by_plane([1,0],-1,0)

    mat = Material()
    mat.set_structural_params(2.1e5,.3)

    bcs = [force,support]

    pmodel = PhysicalModel(mesh,mat,bcs)

    optimizer = StructuralOptim(pmodel,volfrac,5)
    optimizer.run()

    visu = Visualizer(pmodel)
    visu.write_result()

    visu.plot_result()
    plt.show()
```

The result is the follwing structure
![400 Element Example as MINLP](/docs/images/400elements_global.png)

The same setup in ANSYS leads to a different result.
![400 Element Example in ANSYS](/docs/images/400elements_local.png)


## ToDo's
- implement interface for arbitrary meshes and geometries
- implement post processing of results (FEA for final design)
- implement stress minimization 
- implement [decogo](https://github.com/ouyang-w-19/decogo) and scip

![TODO Graph](/docs/images/todos.png)