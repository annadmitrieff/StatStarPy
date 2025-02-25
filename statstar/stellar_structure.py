# stellar_structure.py

import numpy as np
from .constants import dp, G, four_pi, a_rad, four_ac_o3, c, m_H, k_B
from . import zone_quantities as zq
from . import composition as comp
from . import physics

"""
This subroutine/file:
    1. computes the pressure gradient with respect to radius
    2. computes the mass gradient with respect to radius
    3. computes the luminosity gradient with respect to radius
    4. computes the temperature gradient with respect to radius
"""

def structure_eqns(i, r, S):
    """
    Driver for stellar structure equations
    
    Parameters:
    -----------
    i : int
        Equation index
    r : float
        Radius
    S : ndarray
        State vector [P, M_r, L_r, T]
        
    Returns:
    --------
    dfdr : float
        Derivative value
    ok : bool
        True if derivative successful, False otherwise
    """
    ok = True
    
    P = S[0]
    M_r = S[1]
    L_r = S[2]
    T = S[3]
    
    X = zq.Xm
    Z = zq.Zm
    
    Y = comp.helium(X, Z)
    mu = comp.mean_molecular_weight(X, Y, Z)
    zq.rho = physics.density(T, P, mu)
    
    if zq.rho < 0:
        print("Density calculation error in FUNCTION Structure_Eqns")
        ok = False
    
    if i == 0:
        if ok:
            dfdr = dp_dr(M_r, zq.rho, r)
        else:
            dfdr = 0
        zq.dfdr0[0] = dfdr
    
    elif i == 1:
        if ok:
            dfdr = dm_dr(r, zq.rho)
        else:
            dfdr = 0
        zq.dfdr0[1] = dfdr
    
    elif i == 2:
        if ok:
            zq.epsilon = physics.nuclear(T, zq.rho, X, Z)
            dfdr = dl_dr(r, zq.rho, zq.epsilon)
        else:
            dfdr = 0
        zq.dfdr0[2] = dfdr
    
    elif i == 3:
        if ok:
            zq.kappa = physics.opacity(T, zq.rho, X, Z)
            zq.gamma = physics.specific_heat_ratio()
            zq.dlnPdlnT = physics.pt_gradient(zq.Pm, P, zq.Tm, T)
            dfdr = dt_dr(zq.kappa, zq.rho, T, L_r, r, mu, M_r, zq.gamma, zq.dlnPdlnT)
        else:
            dfdr = 0
        zq.dfdr0[3] = dfdr
    
    return dfdr, ok

def dp_dr(M_r, rho, r):
    """
    Hydrostatic Equilibrium
    
    Parameters:
    -----------
    M_r : float
        Mass enclosed within radius r
    rho : float
        Density
    r : float
        Radius
        
    Returns:
    --------
    dPdr : float
        Pressure gradient
    """
    # Eq. (10.6)
    dPdr = -G * M_r * rho / r**2
    return dPdr

def dm_dr(r, rho):
    """
    Mass Conservation
    
    Parameters:
    -----------
    r : float
        Radius
    rho : float
        Density
        
    Returns:
    --------
    dMdr : float
        Mass gradient
    """
    # Eq. (10.7)
    dMdr = four_pi * r**2 * rho
    return dMdr

def dl_dr(r, rho, epsilon):
    """
    Luminosity Gradient
    
    Parameters:
    -----------
    r : float
        Radius
    rho : float
        Density
    epsilon : float
        Energy generation rate
        
    Returns:
    --------
    dLdr : float
        Luminosity gradient
    """
    # Eq. (10.36)
    dLdr = four_pi * r**2 * rho * epsilon
    return dLdr

def dt_dr(kappa, rho, T, L_r, r, mu, M_r, gamma, dlnPdlnT):
    """
    Temperature Gradient
    
    Parameters:
    -----------
    kappa : float
        Opacity
    rho : float
        Density
    T : float
        Temperature
    L_r : float
        Luminosity at radius r
    r : float
        Radius
    mu : float
        Mean molecular weight
    M_r : float
        Mass enclosed within radius r
    gamma : float
        Ratio of specific heats
    dlnPdlnT : float
        Pressure gradient with respect to temperature
        
    Returns:
    --------
    dTdr : float
        Temperature gradient
    """
    gamma_ratio = gamma / (gamma - 1)
    
    # Radiation criterion, Eq. (10.95)
    if dlnPdlnT > gamma_ratio:
        # Radiation, Eq. (10.68)
        dTdr = -(kappa * rho / T**3) * (L_r / (four_pi * r**2)) / four_ac_o3
        zq.rc_flag = 'r'
    else:
        # Adiabatic convection, Eq. (10.89)
        dTdr = -(1 / gamma_ratio) * (mu * m_H / k_B) * (G * M_r / r**2)
        zq.rc_flag = 'c'
    
    return dTdr