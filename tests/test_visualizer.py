from pdb import pm
import numpy as np
from topopt.topopt import Visualizer
import glob, os

def test_write_file_exist():
    class PhysicalModelStub():
        def __init__(self,nelem):
            self.x = np.arange(nelem**2)

    pmodel = PhysicalModelStub(10)

    vis = Visualizer(pmodel)

    # delete file if it exists
    filename = vis.filename
    files = glob.glob(filename)
    print(filename,files)
    if filename in files:
        os.remove(filename)

    vis.write_result()

    assert vis.filename in glob.glob(vis.filename)

def test_data_in_file():
    class PhysicalModelStub():
        def __init__(self,nelem):
            self.x = np.arange(nelem**2)
    
    pmodel = PhysicalModelStub(10)

    vis = Visualizer(pmodel)
    vis.write_result()

    test_data = np.loadtxt(vis.filename)
    os.remove(vis.filename)
    
    assert np.allclose(pmodel.x,test_data)