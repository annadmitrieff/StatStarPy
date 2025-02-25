# user_io.py

import numpy as np
from .constants import dp, M_Sun, L_Sun, R_Sun, four_pi, sigma
from . import composition as comp

# simple, just takes in the SI values and spits out the solar values

def calculate_radius(Ls, Teff):
    """
    Calculate stellar radius from luminosity and effective temperature
    
    Parameters:
    -----------
    Ls : float
        Luminosity in SI units
    Teff : float
        Effective temperature in K
        
    Returns:
    --------
    Rs : float
        Stellar radius in SI units
    """
    # Eq. (3.17)
    return np.sqrt(Ls / (four_pi * sigma * Teff**4))

def calculate_solar_units(Ms, Ls, Rs):
    """
    Convert SI units to solar units
    
    Parameters:
    -----------
    Ms : float
        Mass in kg
    Ls : float
        Luminosity in W
    Rs : float
        Radius in m
        
    Returns:
    --------
    Msolar : float
        Mass in solar units
    Lsolar : float
        Luminosity in solar units
    Rsolar : float
        Radius in solar units
    """
    Msolar = Ms / M_Sun
    Lsolar = Ls / L_Sun
    Rsolar = Rs / R_Sun
    return Msolar, Lsolar, Rsolar

def si_to_solar_units(value, unit_type):
    """
    Convert a value from SI to solar units
    
    Parameters:
    -----------
    value : float
        Value in SI units
    unit_type : str
        Type of unit ('mass', 'luminosity', or 'radius')
        
    Returns:
    --------
    float
        Value in solar units
    """
    if unit_type.lower() == 'mass':
        return value / M_Sun
    elif unit_type.lower() == 'luminosity':
        return value / L_Sun
    elif unit_type.lower() == 'radius':
        return value / R_Sun
    else:
        raise ValueError(f"Unknown unit type: {unit_type}")

def solar_to_si_units(value, unit_type):
    """
    Convert a value from solar to SI units
    
    Parameters:
    -----------
    value : float
        Value in solar units
    unit_type : str
        Type of unit ('mass', 'luminosity', or 'radius')
        
    Returns:
    --------
    float
        Value in SI units
    """
    if unit_type.lower() == 'mass':
        return value * M_Sun
    elif unit_type.lower() == 'luminosity':
        return value * L_Sun
    elif unit_type.lower() == 'radius':
        return value * R_Sun
    else:
        raise ValueError(f"Unknown unit type: {unit_type}")