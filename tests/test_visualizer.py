import numpy as np
from topopt.topopt import Visualizer
import glob, os, pytest

class PhysicalModelStub():
    def __init__(self,nelem):
        self.x = np.arange(nelem**2)

def test_write_file_exist():

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
    
    pmodel = PhysicalModelStub(10)

    vis = Visualizer(pmodel)
    vis.write_result()

    test_data = np.loadtxt(vis.filename)
    os.remove(vis.filename)
    
    assert np.allclose(pmodel.x,test_data)

@pytest.mark.parametrize("ndiv",range(1,100,10))
def test_vis_data(ndiv):
    class PhysicalModelStub():
        def __init__(self,nelem):
            self.x = np.arange(nelem**2)
    pmodel = PhysicalModelStub(ndiv)
    
    vis = Visualizer(pmodel)
    plot_data = vis.plot_result()
    
    assert list(pmodel.x) == list(plot_data.get_array().flatten())
