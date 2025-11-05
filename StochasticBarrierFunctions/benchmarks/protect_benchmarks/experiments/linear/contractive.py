# ===============================================================
# Two-tank stochastic system test
# ===============================================================

import time
import sympy as sp
import numpy as np

from src.functions.parallel_dt_SS import parallel_dt_SS
from src.functions.dt_SS import dt_SS

# ========================= Parameters =========================
if __name__ == '__main__':
    dim = 2  # dimension of state space

    # Initial set
    L_initial = np.array([-0.05, -0.05])
    U_initial = np.array([0.05, 0.05])

    # Unsafe set1
    L_unsafe1 = np.array([-2, -2])
    U_unsafe1 = np.array([-1, -1])

    # Unsafe set2
    L_unsafe2 = np.array([2, 2])
    U_unsafe2 = np.array([3, 3])

    # Combine unsafe regions
    L_unsafe = np.array([L_unsafe1, L_unsafe2])
    U_unsafe = np.array([U_unsafe1, U_unsafe2])

    # State space
    L_space = np.array([-1, -1])
    U_space = np.array([2, 2])

    # ========================= Symbolic Variables =========================
    x = sp.symbols(f'x0:{dim}')       # Create x0, x1
    varsigma = sp.symbols(f'varsigma0:{dim}')  # noise symbols

    print("State variables:", x)
    print("Noise variables:", varsigma)

    # ========================= Dynamics =========================
    NoiseType = "normal"
    sigma = np.array([0.1, 0.1])
    mean = np.array([0, 0])

    f1 = 0.95 * x[0] + varsigma[0]
    f2 = 0.95 * x[1] + varsigma[1]
    f = np.array([f1, f2])

    # Time horizon
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
        'solver': "cvxpy",
        'confidence': None,
        'gam': None,
        'lam': 10,
        'c_val': None,
        'sigma': sigma,
        'mean': mean,
        'rate': None,
        'a': None,
        'b': None,
    }

    # ========================= Run for multiple degrees =========================
    degrees = [2, 4, 8, 12, 24, 30]

    for degree in degrees:
        start = time.time()
        print("\n>>> Running dt_SS() Contractive Model with degree =", degree)
        result = dt_SS(degree, **fixed_params)
        end = time.time()
        print("Elapsed time:", end - start)

        if not result:
            print("Results dictionary is empty.")
        else:
            print("Result:", result)
