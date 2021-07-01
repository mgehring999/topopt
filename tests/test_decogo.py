from simple import model
import sys, os

# add Decogo source files to path
# very important to do it before importing Decogo
path = os.path.dirname(os.getcwd())
sys.path.insert(0, path+"\decogo\decogo")
from solver.decogo import Decogo
input("enter")

# overwrite solver default settings and put the file in the working directory
with open('decogo.set', 'w') as file:
    file.write('strategy = OA\n')
    file.write('decomp_estimate_var_bounds = True')
    file.close()

# create solver instance
solver = Decogo()

# solve the model
solver.optimize(model)

os.remove('decogo.set') # remove settings file