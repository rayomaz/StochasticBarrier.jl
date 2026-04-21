# IMPORTS FROM INSTALLS
import sympy as sp
import numpy as np

# IMPORTS FROM TOOL
from src.functions.dt_SS import dt_SS
from experiments._record import record
# ========================= Parameters =========================

if __name__ == '__main__':

    dim = 1  # dimension of state space
    
    # Initial set
    L_initial = np.array([19.5])
    U_initial = np.array([20])

    # Unsafe set1
    L_unsafe1 = np.array([1])
    U_unsafe1 = np.array([17])

    # Unsafe set2
    L_unsafe2 = np.array([23])
    U_unsafe2 = np.array([50])

    # combine unsafe regions
    L_unsafe = np.array([L_unsafe1, L_unsafe2])
    U_unsafe = np.array([U_unsafe1, U_unsafe2])

    # State space
    L_space = np.array([1])
    U_space = np.array([50])

    # ========================= Symbolic Variables =========================
    x = sp.symbols(f'x0:{dim}')  # Create x1, x2, ..., x_degree symbols
    varsigma = sp.symbols(f'varsigma0:{dim}')
    # ========================= Dynamics =========================
    f1 = 4.32 + 0.7455988625*x[0] + 0.0017422474999999999*x[0]**2 + varsigma[0]

    # Define the vector field
    f = np.array([f1])

    # Noise terms
    NoiseType = "normal"
    sigma = np.array([0.1])
    mean = np.array([0])
    
    # Ttime horizon
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
        'lam': 10,
        'c_val': None,
        'sigma': sigma,
        'mean': mean,
        'rate': None,
        'a': None,
        'b': None,
        # Add other fixed parameters here
    }

    # ------------------ Run for multiple degrees -----------------
    degrees = [6, 8, 10]

    for degree in degrees:
        print("\n>>> Running dt_SS() Thermostat Model with degree =", degree)
        result = record(dt_SS,
                        system="1D Thermostat",
                        cls="Polynomial", table=3, degree=degree,
                        **fixed_params)
        if not result:
            print("Results dictionary is empty.")
        else:
            print("Result:", result)
