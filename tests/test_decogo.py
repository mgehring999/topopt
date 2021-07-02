from simple import from_decogo_docs, simple_scalar, simple_scalar_nonIndexed
import sys, os

from decogo.solver.decogo import DecogoSolver

if __name__ == "__main__":
    
    #model = from_decogo_docs()
    model = simple_scalar()                 # this isnt running
    #model = simple_scalar_nonIndexed()     # this is fine
    
    # overwrite solver default settings and put the file in the working directory
    with open('decogo.set', 'w') as file:
        file.write('strategy = OA\n')
        file.write('decomp_estimate_var_bounds = True')
        file.close()

    # create solver instance
    solver = DecogoSolver()

    # solve the model
    solver.optimize(model)

    os.remove('decogo.set') # remove settings file