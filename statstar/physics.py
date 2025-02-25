# physics.py

"""
This subroutine (file?) computes:
    1. the pressure gradient with respect to temperature
    2. the specific heat ratio C_P/C_V
    3. the density of the gas assuming the ideal gas law and radiation pressure
    4. the opacity
    5. the change in optical depth across the zone
    6. the nuclear energy generation rates

"""

import numpy as np
from .constants import dp, a_rad_o3, k_B, m_H
from . import zone_quantities as zq
from . import composition as comp

def pt_gradient(Pm, P, Tm, T):
    """
    Compute the pressure gradient with respect to temperature 
    to determine whether convection is required.
    
    Parameters:
    -----------
    Pm : float
        Previous zone pressure
    P : float
        Current zone pressure
    Tm : float
        Previous zone temperature
    T : float
        Current zone temperature
        
    Returns:
    --------
    dlnPdlnT : float
        Pressure gradient with respect to temperature
    """
    dlnPdlnT = ((Tm + T)/(Pm + P))*((Pm - P)/(Tm - T))
    if dlnPdlnT > 99.9:
        dlnPdlnT = 99.9
    return dlnPdlnT

def specific_heat_ratio():
    """
    Compute the ratio C_P/C_V
    
    Returns:
    --------
    gamma : float
        Ratio of specific heats
    """
    # Assume a purely monatomic gas, Eq. (10.80)
    monatomic = 5.0/3.0
    return monatomic

def density(T, P, mu):
    """
    Compute the density of the gas, assuming the ideal gas law and radiation pressure
    
    Parameters:
    -----------
    T : float
        Temperature
    P : float
        Pressure
    mu : float
        Mean molecular weight
        
    Returns:
    --------
    rho : float
        Gas density (negative if error)
    """
    # Eq. (10.20)
    P_gas = P - a_rad_o3 * T**4
    
    if P_gas <= 0.0 and T > 0.0:
        # "Do something desperate"
        if zq.step_size_condition == 0:
            P_gas = P
        elif zq.step_size_condition == 1:
            P_gas = 0.001 * P
        elif zq.step_size_condition == 2:
            P_gas = 0.0001 * P
    
    if T > 0.0 and P_gas > 0.0:
        # Eq. (10.11)
        rho = P_gas * mu * m_H / (k_B * T)
    else:
        rho = -1.0
    
    if rho < 0.0:
        print("A negative density was computed!")
        print("Sorry but I am not programmed to handle this new physics :-)") # adds character
        print(f"Terminating calculation with:")
        print(f"         T     = {T:.6e}")
        print(f"         P     = {P:.6e}")
        print(f"         P_gas = {P_gas:.6e}")
    
    return rho

def opacity(T, rho, X, Z):
    """
    Compute an approximation of the Rosseland Mean Opacity
    
    Parameters:
    -----------
    T : float
        Temperature
    rho : float
        Density
    X : float
        Hydrogen mass fraction
    Z : float
        Metals mass fraction
        
    Returns:
    --------
    kappa : float
        Opacity
    """
    # The free-free Gaunt factor is on the order of unity
    g_ff = 1.0
    # Coefficients for bound-free, free-free, electron scattering, and H- opacity
    A_bf = 4.34e21
    A_ff = 3.68e18
    A_es = 0.02
    A_Hm = 7.9e-34
    
    # From Novotny (1973), p. 469
    tog_bf = 0.708 * (rho * (1.0 + X))**0.2
    
    # Eq. (9.22)
    kappa_bf = (A_bf / tog_bf) * Z * (1.0 + X) * rho / T**3.5
    
    # Eq. (9.23)
    kappa_ff = A_ff * g_ff * (1.0 - Z) * (1.0 + X) * rho / T**3.5
    
    # Eq. (9.27)
    kappa_es = A_es * (1.0 + X)
    
    # Check if H- opacity is applicable
    if ((T > 3000.0 and T < 6000.0) and 
        (rho > 1.0e-10 and rho < 1.0e-5) and 
        (Z > 0.001 and Z < 0.03)):
        # Eq. (9.28)
        kappa_Hminus = A_Hm * (Z / 0.02) * np.sqrt(rho) * T**9
    else:
        kappa_Hminus = 0.0
    
    # Total opacity
    kappa = kappa_bf + kappa_ff + kappa_es + kappa_Hminus
    
    return kappa

def optical_depth_change(kappa, kappam, rho, rhom, r, rm):
    """
    Compute the change in optical depth across the zone
    
    Parameters:
    -----------
    kappa : float
        Current opacity
    kappam : float
        Previous opacity
    rho : float
        Current density
    rhom : float
        Previous density
    r : float
        Current radius
    rm : float
        Previous radius
        
    Returns:
    --------
    dtau : float
        Change in optical depth
    """
    # Eq. (9.15)
    dtau = -(kappa * rho + kappam * rhom) * (r - rm) / 2.0
    
    return dtau

def nuclear(T, rho, X, Z):
    """
    Compute the nuclear energy generation rates
    
    Parameters:
    -----------
    T : float
        Temperature
    rho : float
        Density
    X : float
        Hydrogen mass fraction
    Z : float
        Metals mass fraction
        
    Returns:
    --------
    epsilon : float
        Total energy generation rate
    """
    # Screening factors
    fpp = 1.0
    f3a = 1.0
    
    # Reaction rate coefficients
    A_pp = 0.241
    A_CNO = 8.67e20
    A_He = 50.9
    
    # Temperature in millions and hundred millions of Kelvin (a lot!)
    T6 = T * 1.0e-6
    T8 = T * 1.0e-8
    
    # Fractional powers
    onethird = 1.0/3.0
    twothirds = 2.0/3.0
    fourthirds = 4.0/3.0
    fivethirds = 5.0/3.0
    
    # PP chains (see Hansen and Kawaler, Eq. 6.65, 6.73, and 6.74)
    psipp = 1.0 + 1.412e8 * (1.0/X - 1.0) * np.exp(-49.98 * T6**(-onethird))
    Cpp = 1.0 + 0.0123 * T6**onethird + 0.0109 * T6**twothirds + 0.000938 * T6
    # Eq. (10.46)
    eps_pp = A_pp * rho * X * X * fpp * psipp * Cpp * T6**(-twothirds) * np.exp(-33.80 * T6**(-onethird))
    
    # CNO cycle (Kippenhahn and Weigert, Eq. 18.65). Not our Dr. Weigert--different one :-(
    XCNO = comp.cno(Z)
    CCNO = 1.0 + 0.0027 * T6**onethird - 0.00778 * T6**twothirds - 0.000149 * T6
    # Eq. (10.58)
    eps_CNO = A_CNO * rho * X * XCNO * CCNO * T6**(-twothirds) * np.exp(-152.28 * T6**(-onethird))
    
    # Helium burning (Kippenhahn and Weigert, Eq. 18.67)
    Y = comp.helium(X, Z)
    # Eq. (10.62)
    eps_He = A_He * rho**2 * Y**3 / T8**3 * f3a * np.exp(-44.027 / T8)
    
    # Combined energy generation rate
    epsilon = eps_pp + eps_CNO + eps_He
    
    return epsilon