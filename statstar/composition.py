# composition.py

import numpy as np
from .constants import dp

def mean_molecular_weight(X, Y, Z):
    """
    Calculate the mean molecular weight of the gas mixture
    
    Parameters:
    -----------
    X : float
        Hydrogen mass fraction
    Y : float
        Helium mass fraction
    Z : float
        Metals mass fraction
        
    Returns:
    --------
    mu : float
        Mean molecular weight
    """
    # Here we assume complete ionization, Eq. (10.16)
    return 1.0 / (2 * X + 3 * Y / 4 + Z / 2)

def helium(X, Z):
    """
    Calculate the amount of Helium-4 in the mixture
    
    Parameters:
    -----------
    X : float
        Hydrogen mass fraction
    Z : float
        Metals mass fraction
        
    Returns:
    --------
    Y : float
        Helium mass fraction
    """
    return 1.0 - X - Z

def cno(Z):
    """
    Calculate the mass fraction of C, N, and O in the mixture
    
    Parameters:
    -----------
    Z : float
        Metals mass fraction
        
    Returns:
    --------
    XCNO : float
        CNO mass fraction
    """
    return Z / 2.0