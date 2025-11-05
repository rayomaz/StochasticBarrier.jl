# IMPORTS FROM INSTALLS
import time
import sympy as sp
import numpy as np

# IMPORTS FROM TOOL
from src.functions.parallel_dt_SS import parallel_dt_SS
from src.functions.dt_SS import dt_SS
# ========================= Parameters =========================

if __name__ == '__main__':

    dim = 2  # dimension of state space

    # Initial set
    L_initial = np.array([-5,-5])
    U_initial = np.array([5,5])

    # Unsafe set1
    L_unsafe1 = np.array([-7,-7])
    U_unsafe1 = np.array([-6,-6])
    
    # Unsafe set2
    L_unsafe2 = np.array([6,6])
    U_unsafe2 = np.array([7,7])

    # combine unsafe regions
    L_unsafe = np.array([L_unsafe1,L_unsafe2])
    U_unsafe = np.array([U_unsafe1,U_unsafe2])

    # State space
    L_space = np.array([-6,-6])
    U_space = np.array([6,6])

    # ========================= Symbolic Variables =========================
    x = sp.symbols(f'x0:{dim}')  # Create x1, x2, ..., x_degree symbols
    varsigma = sp.symbols(f'varsigma0:{dim}')
    # ========================= Dynamics =========================

    #noise term
    NoiseType = "normal"
    sigma = np.array([0.02, 0.02])
    mean = np.array([0, 0])
    
    tau = 0.1
    
    f1 = x[0] + tau*x[1] + varsigma[0]
    f2 = x[1] + (-x[0] + (1-x[0]**2)*x[1])*tau + varsigma[1]
    
    # Define the vector field
    f = np.array([f1,f2])
    
    #time horizon
    t = 10

    fixed_params = {
        'dim': dim,
        'L_initial': L_initial,
        'U_initial': U_initial,
        'L_unsafe': L_unsafe,
        'U_unsafe': U_unsafe,
        'L_space': L_space,
        'U_space': U_space,
        'x': x,
        'varsigma': varsigma,
        'f': f,
        't': t,
        'noise_type': NoiseType,
        'optimize': True,
        'solver': "mosek",
        'confidence': None,
        'gam': None,
        'lam': 1000,
        'c_val': None,
        'sigma': sigma,
        'mean' : mean,
        'rate': None,
        'a': None,
        'b': None,
        # Add other fixed parameters here
    }

    # ------------------ Run for multiple degrees -----------------
    degrees = [2, 4, 6, 8, 12]

    for degree in degrees:
        start = time.time()
        print("\n>>> Running dt_SS() Oscillator Model with degree =", degree)
        result = dt_SS(degree, **fixed_params)
        end = time.time()
        print("Elapsed time:", end - start)

        if not result:
            print("Results dictionary is empty.")
        else:
            print("Result:", result)