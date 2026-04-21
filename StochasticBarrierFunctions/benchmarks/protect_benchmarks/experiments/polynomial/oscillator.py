# IMPORTS FROM INSTALLS
import sympy as sp
import numpy as np

# IMPORTS FROM TOOL
from src.functions.dt_SS import dt_SS
from experiments._record import record
# ========================= Parameters =========================

if __name__ == '__main__':

    dim = 2  # dimension of state space

    # Initial set
    L_initial = np.array([-5,-5])
    U_initial = np.array([5,5])

    # Unsafe set1
    L_unsafe1 = np.array([-7,-7])
    U_unsafe1 = np.array([-6,7])

    # Unsafe set2
    L_unsafe2 = np.array([-7,-7])
    U_unsafe2 = np.array([7,-6])
    
    # Unsafe set3
    L_unsafe3 = np.array([6,-7])
    U_unsafe3 = np.array([7,7])

    # Unsafe set4
    L_unsafe4 = np.array([-7,6])
    U_unsafe4 = np.array([7,7])

    # combine unsafe regions
    L_unsafe = np.array([L_unsafe1,L_unsafe2,L_unsafe3,L_unsafe4])
    U_unsafe = np.array([U_unsafe1,U_unsafe2,U_unsafe3,U_unsafe4])

    # State space
    L_space = np.array([-7,-7])
    U_space = np.array([7,7])

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
        'solver': "cvxpy",
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
    degrees = [6, 8, 12]

    for degree in degrees:
        print("\n>>> Running dt_SS() Oscillator Model with degree =", degree)
        result = record(dt_SS,
                        system="2D Oscillator",
                        cls="Polynomial", table=3, degree=degree,
                        **fixed_params)
        if not result:
            print("Results dictionary is empty.")
        else:
            print("Result:", result)