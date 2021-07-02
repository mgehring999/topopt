from simple import model
import sys, os

from decogo.solver.decogo import DecogoSolver

if __name__ == "__main__":
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