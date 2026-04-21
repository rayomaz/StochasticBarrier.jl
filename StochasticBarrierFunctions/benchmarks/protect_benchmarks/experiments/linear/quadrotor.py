# ===============================================================
# Linear 6D system test (adapted from room temperature model)
# ===============================================================

import sympy as sp
import numpy as np

from src.functions.dt_SS import dt_SS
from experiments._record import record

if __name__ == '__main__':
    # ------------------ System dimensions ------------------
    dim = 6

    # ------------------ System matrices -------------------
    A = np.array([
        [-4.84088e-16, -3.48202e-15, 0.274208, 0.701275, -0.00237413, -0.943452],
        [-0.00020332, -0.000245833, 3.61414e-17, 7.80912e-15, 2.58941e-15, 1.6826e-15],
        [0.00056152, 0.000488716, -0.000195302, -0.00038964, -0.000108086, -0.000212995],
        [-0.000245933, -0.000233841, -0.000165083, -0.000359714, 0.000206617, 0.000410682],
        [0.000186356, 0.000165408, 0.000319849, 0.000651892, -0.000123525, -0.000245799],
        [-0.000398372, -0.000373963, -0.000282334, -0.000600204, 0.000372937, 0.000757686]
    ])
    b = np.zeros(dim)
    sigma = np.array([0.001]*dim)

    # ------------------ Unsafe regions --------------------
    L_unsafe1 = np.array([-1.0, -2.0, -0.2, -0.2, -1.0, -1.0])
    U_unsafe1 = np.array([-0.5, 2.0, 0.2, 0.2, 4.0, 2.0])

    L_unsafe2 = np.array([2.0, -2.0, -0.2, -0.2, -1.0, -1.0])
    U_unsafe2 = np.array([3.0, 2.0, 0.2, 0.2, 4.0, 2.0])

    L_unsafe3 = np.array([-1.0, -2.0, -0.2, -0.2, -1.0, -1.0])
    U_unsafe3 = np.array([3.0, -1.0, 0.2, 0.2, 4.0, 2.0])

    L_unsafe4 = np.array([-1.0, 1.0, -0.2, -0.2, -1.0, -1.0])
    U_unsafe4 = np.array([3.0, 2.0, 0.2, 0.2, 4.0, 2.0])

    L_unsafe5 = np.array([-1.0, -2.0, -0.2, -0.2, -1.0, -1.0])
    U_unsafe5 = np.array([3.0, 2.0, -0.1, 0.2, 4.0, 2.0])

    L_unsafe6 = np.array([-1.0, -2.0, 0.1, -0.2, -1.0, -1.0])
    U_unsafe6 = np.array([3.0, 2.0, 0.2, 0.2, 4.0, 2.0])

    L_unsafe7 = np.array([-1.0, -2.0, -0.2, -0.2, -1.0, -1.0])
    U_unsafe7 = np.array([3.0, 2.0, 0.2, -0.1, 4.0, 2.0])

    L_unsafe8 = np.array([-1.0, -2.0, -0.2, 0.1, -1.0, -1.0])
    U_unsafe8 = np.array([3.0, 2.0, 0.2, 0.2, 4.0, 2.0])

    L_unsafe9 = np.array([-1.0, -2.0, -0.2, -0.2, -1.0, -1.0])
    U_unsafe9 = np.array([3.0, 2.0, 0.2, 0.2, -0.5, 2.0])

    L_unsafe10 = np.array([-1.0, -2.0, -0.2, -0.2, 3.0, -1.0])
    U_unsafe10 = np.array([3.0, 2.0, 0.2, 0.2, 4.0, 2.0])

    L_unsafe11 = np.array([-1.0, -2.0, -0.2, -0.2, -1.0, -0.5])
    U_unsafe11 = np.array([3.0, 2.0, 0.2, 0.2, 4.0, 2.0])
    
    L_unsafe12 = np.array([-1.0, -2.0, -0.2, -0.2, -1.0, 1.5])
    U_unsafe12 = np.array([3.0, 2.0, 0.2, 0.2, 4.0, 2.0])

    # Combine unsafe regions
    L_unsafe = np.array([L_unsafe1, L_unsafe2, L_unsafe3, L_unsafe4, L_unsafe5, L_unsafe6, L_unsafe7, L_unsafe8, L_unsafe9, L_unsafe10, L_unsafe11, L_unsafe12])
    U_unsafe = np.array([U_unsafe1, U_unsafe2, U_unsafe3, U_unsafe4, U_unsafe5, U_unsafe6, U_unsafe7, U_unsafe8, U_unsafe9, U_unsafe10, U_unsafe11, U_unsafe12])

    # ------------------ State space -----------------------
    L_space = np.array([-1.0, -2.0, -0.2, -0.2, -1.0, -1.0])
    U_space = np.array([3.0, 2.0, 0.2, 0.2, 4.0, 2.0])

    # ------------------ Initial region --------------------
    L_initial = np.array([-0.01, 0.01, -0.01, 0.01, -0.01, 0.01])
    U_initial = np.array([0.01, 0.03, 0.01, 0.03, 0.01, 0.03])

    # ------------------ Symbolic variables ----------------
    x = sp.symbols(f'x0:{dim}')
    varsigma = sp.symbols(f'varsigma0:{dim}')

    print("State variables:", x)
    print("Noise variables:", varsigma)

    # ------------------ Dynamics -------------------------
    # Linear system: f(x) = A*x + b + noise
    f = A @ sp.Matrix(x) + sp.Matrix(b) + sp.Matrix(varsigma)

    # ------------------ Fixed parameters -----------------
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
        't': 10,
        'noise_type': "normal",
        'optimize': True,
        'solver': "mosek",
        'confidence': None,
        'gam': None,
        'lam': 10,
        'c_val': None,
        'sigma': sigma,
        'mean': np.zeros(dim),
        'rate': None,
        'a': None,
        'b': None
    }

    # ------------------ Run for multiple degrees -----------------
    degrees = [2]

    for degree in degrees:
        print("\n>>> Running dt_SS() Quadrotor Model with degree =", degree)
        result = record(dt_SS,
                        system="6D Quadrotor",
                        cls="Linear", table=3, degree=degree,
                        **fixed_params)
        if not result:
            print("Results dictionary is empty.")
        else:
            print("Result:", result)
