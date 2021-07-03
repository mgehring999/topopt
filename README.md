# global-topopt
This is a python implementation of a two-dimensional structural optimizitation problem as a MINLP (mixed-integer non-linear programm) written during summer term 2021 as a class project. 

## installation
Install python requirements with pip
```
pip install -r requirements.txt
```

Download SCIPOptSuite 7.0.0 and scipampl 7.0.0 from scipopt.org
install scipOptSuite and put scipampl binary in scipOptSuite\bin
add \bin to PATH

## usage
Start the main script
```
py topo.py <volfrac> <ndivisions> <load_xcoord,load_ycoord> <load_xmagnitude,load_ymagnitude>
```

## example
This command starts the solver with a volume fraction constraint of 50% and 400 Elements (20x20 Elements). A single force load of -1 in y-direction is applied on the node at the coordinate (1,0)
```
py topo.py 0.5 20 1,0 0,-1
```

The result is the follwing structure
![400 Element Example as MINLP](/images/400elements_global.png)

The same setup in ANSYS leads to a different result. ANSYS uses a local optimiziation scheme.
![400 Element Example in ANSYS](/images/400elements_local.png)


## ToDo's
- implement interface for arbitrary meshes and geometries
- logging 
- save and load resultfile