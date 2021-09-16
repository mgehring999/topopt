# ``topopt`` - Topology Optimization in Python 
This is a python implementation of a two-dimensional structural optimizitation framework started during summer term 2021 as a class project and continued since then.

## Installation
Install python requirements with pip
```
pip install -r requirements.txt
py setup.py install
```
To setup the solver do the following steps
- Download SCIPOptSuite 7.0.0 and scipampl 7.0.0 from [scipopt.org](https://scipopt.org/)
- install scipOptSuite and put scipampl binary in scipOptSuite\bin
- add \bin to PATH

## Example script
This script models a quadratic plate clamped on the left side and loaded with a point force in the middle of the right edge. The plate is discretized with 100 linear plate elements. The optimization goal is to delete 50% of the volume and minimize the compliance. 

```python
    from topopt.topopt import StructuralOptim, Visualizer
    from topopt.physical import PhysicalModel, Material
    from topopt.mesh import Displacement, Load, Mesh

    volfrac = 0.5
    ndiv = 20

    mesh = Mesh()
    mesh.rect_mesh(ndiv)

    loads = Load(mesh)
    loads.add_by_coord((1,0),(0,-1))

    support = Displacement()
    support.add_by_edge("left","all")

    mat = Material()
    mat.set_structural_params(2.1e5,.3)

    pmodel = PhysicalModel(mesh,mat,support,loads)

    optimizer = StructuralOptim(pmodel,volfrac,5)
    optimizer.run()

    visu = Visualizer(pmodel)
    visu.show_result()
```

The result is the follwing structure
![400 Element Example as MINLP](/docs/images/400elements_global.png)

The same setup in ANSYS leads to a different result.
![400 Element Example in ANSYS](/docs/images/400elements_local.png)


## ToDo's
- implement interface for arbitrary meshes and geometries
- save and load resultfile
- implement post processing of results (FEA for final design)
- implement stress minimization 
- implement [decogo](https://github.com/ouyang-w-19/decogo) and scip

![TODO Graph](/docs/images/todos.png)