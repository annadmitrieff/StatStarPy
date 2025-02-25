# ode_integrator.py

# runge-kutta ODE integrator. Same approach as in the book.

import numpy as np
from .constants import dp

def rk_4(n, x0, h, y0, f0, f):
    """
    Fourth-order Runge-Kutta integration step
    
    Parameters:
    -----------
    n : int
        Number of ODEs
    x0 : float
        Initial value of independent variable
    h : float
        Step size
    y0 : ndarray
        Initial values of dependent variables
    f0 : ndarray
        Initial derivatives
    f : callable
        Function returning derivatives
        
    Returns:
    --------
    y4 : ndarray
        Updated y values after one step
    ok : bool
        True if integration successful, False otherwise
    """
    ok = True
    
    # Calculate intermediate derivatives
    k1 = h * f0
    
    k2 = np.zeros(n, dtype=dp)
    for i in range(n):
        derivative, ok = f(i, x0 + h/2, y0 + k1/2)
        if not ok:
            return None, False
        k2[i] = h * derivative
    
    k3 = np.zeros(n, dtype=dp)
    for i in range(n):
        derivative, ok = f(i, x0 + h/2, y0 + k2/2)
        if not ok:
            return None, False
        k3[i] = h * derivative
    
    k4 = np.zeros(n, dtype=dp)
    for i in range(n):
        derivative, ok = f(i, x0 + h, y0 + k3)
        if not ok:
            return None, False
        k4[i] = h * derivative
    
    # Compute the variables for the next shell using the 4th order Runge-Kutta formula
    y4 = y0 + k1/6 + k2/3 + k3/3 + k4/6
    
    return y4, ok