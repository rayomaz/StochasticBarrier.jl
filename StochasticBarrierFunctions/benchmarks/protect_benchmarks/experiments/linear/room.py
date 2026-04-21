# IMPORTS FROM INSTALLS
import sympy as sp
import numpy as np

# IMPORTS FROM TOOL
from src.functions.dt_SS import dt_SS
from experiments._record import record
# ========================= Parameters =========================

if __name__ == '__main__':

    dim = 3  # dimension of state space

    # Initial set
    L_initial = np.array([21,21,21])
    U_initial = np.array([22,22,22])

    # Unsafe set1
    L_unsafe1 = np.array([16,16,16])
    U_unsafe1 = np.array([17,30,30])

    # Unsafe set2
    L_unsafe2 = np.array([29,16,16])
    U_unsafe2 = np.array([30,30,30])

    # Unsafe set3
    L_unsafe3 = np.array([16,16,16])
    U_unsafe3 = np.array([30,17,30])

    # Unsafe set4
    L_unsafe4 = np.array([16,29,16])
    U_unsafe4 = np.array([30,30,30])

    # Unsafe set5
    L_unsafe5 = np.array([16,16,16])
    U_unsafe5 = np.array([30,30,17])

    # Unsafe set6
    L_unsafe6 = np.array([16,16,29])
    U_unsafe6 = np.array([30,30,30])

    # combine unsafe regions
    L_unsafe = np.array([L_unsafe1, L_unsafe2, L_unsafe3, L_unsafe4, L_unsafe5, L_unsafe6])
    U_unsafe = np.array([U_unsafe1, U_unsafe2, U_unsafe3, U_unsafe4, U_unsafe5, U_unsafe6])

    # State space
    L_space = np.array([16,16,16])
    U_space = np.array([30,30,30])

    # ========================= Symbolic Variables =========================
    x = sp.symbols(f'x0:{dim}')  # Create x1, x2, ..., x_degree symbols
    varsigma = sp.symbols(f'varsigma0:{dim}')
    # ========================= Dynamics =========================

    #noise terms
    NoiseType = "normal"
    sigma = np.array([0.1,0.1,0.1])
    mean = np.array([0,0,0])
    
    T_e = 10
    alpha_e = 8e-3
    alpha = 6.2e-3
    tau = 5
    
    f1 = (1-tau*(alpha+alpha_e))*x[0] + tau*alpha*x[1] + tau*alpha_e*T_e + varsigma[0]
    f2 = (1-tau*(2*alpha+alpha_e))*x[1] + tau*alpha*(x[0]+x[2]) + tau*alpha_e*T_e + varsigma[1]
    f3 = (1-tau*(alpha+alpha_e))*x[2] + tau*alpha*x[1] + tau*alpha_e*T_e + varsigma[2]

    # Define the vector field
    f = np.array([f1,f2,f3])
    
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
        'lam': 10,
        'c_val': None,
        'sigma': sigma,
        'mean' : mean,
        'rate': None,
        'a': None,
        'b': None,
        # Add other fixed parameters here
    }

    # ========================= Run for multiple degrees =========================
    degrees = [6, 8, 14]

    for degree in degrees:
        print("\n>>> Running dt_SS() Room Temperature Model with degree =", degree)
        result = record(dt_SS,
                        system="3D Room Temp.",
                        cls="Linear", table=3, degree=degree,
                        **fixed_params)
        if not result:
            print("Results dictionary is empty.")
        else:
            print("Result:", result)