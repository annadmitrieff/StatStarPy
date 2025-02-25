# boundary_conditions.py

import numpy as np
from .constants import dp, pi, two_pi, four_pi_o3, G, a_rad, a_rad_o3, c, k_B, m_H
from . import composition as comp
from . import physics
from . import stellar_structure as ss
from . import zone_quantities as zq

def surface(i, Ms, Ls, rm, X, Z, dr):
    """
    This sequence estimates surface boundary conditions.
    
    Parameters:
    -----------
    i : int
        Zone number
    Ms : float
        Total mass
    Ls : float
        Total luminosity
    rm : float
        Previous radius
    X : float
        Hydrogen mass fraction
    Z : float
        Metals mass fraction
    dr : float
        Radius step size
        
    Returns:
    --------
    r : float
        Radius
    P : float
        Pressure
    T : float
        Temperature
    M_r : float
        Mass enclosed within radius r
    L_r : float
        Luminosity at radius r
    rho : float
        Density
    kappa : float
        Opacity
    epsilon : float
        Energy generation rate
    good_surface : bool
        True if surface conditions are valid
    """
    # Maximum change in Ms and Ls over surface zone
    maximum = 1.0e-8
    j_max = 50
    
    r = rm + dr
    
    Y = comp.helium(X, Z)
    mu = comp.mean_molecular_weight(X, Y, Z)
    zq.gamma = physics.specific_heat_ratio()
    gamma_ratio = zq.gamma / (zq.gamma - 1)
    
    j = 0
    while True:
        # Compute the temperature and pressure for the radiative boundary condition
        zq.rc_flag = 'r'

        # Eq. (L.2); radiative assumption
        T = G * Ms * (mu * m_H / (4.25 * k_B)) * (1/r - 1/rm)
        
        if i < 2:
            # Assume small value for surfacea
            tog_bf = 0.01
        else:
            # From Novotny (1973), p. 469
            tog_bf = 2.82 * (zq.rho * (1 + X))**0.2
            
        # From Eq. (9.22) and (9.23)
        g_ff = 1.0
        A_bf = 4.34e21
        A_ff = 3.68e18
        Aop = (A_bf / tog_bf) * Z * (1 + X) + A_ff * g_ff * (1 - Z) * (1 + X)
        
        # Eq. (L.1)
        P = np.sqrt((1/4.25) * (16 * pi/3) * (G * Ms / Ls) * (a_rad * c * k_B / (Aop * mu * m_H))) * T**4.25
        
        # If the zone is convective, recompute the adiabatic temperature and pressure...
        zq.dlnPdlnT = physics.pt_gradient(zq.Pm, P, zq.Tm, T)
        if zq.dlnPdlnT < gamma_ratio and i > 2:
            zq.rc_flag = 'c'
            kPadiabatic = zq.Pm / zq.Tm**gamma_ratio
            # Eq. (L.3)
            T = G * Ms * (mu * m_H / (k_B * gamma_ratio)) * (1/r - 1/rm)
            # Eq. (10.83)
            P = kPadiabatic * T**gamma_ratio
        
        # Compute remaining surface quantities
        rho = physics.density(T, P, mu)
        if rho < 0:
            good_surface = False
            break
            
        kappa = physics.opacity(T, rho, X, Z)
        XCNO = comp.cno(Z)
        epsilon = physics.nuclear(T, rho, X, Z)
        
        # Test to be sure that variations in M_r and L_r are not too large
        M_r = zq.Mm + ss.dm_dr(r, rho) * dr
        L_r = zq.Lm + ss.dl_dr(r, rho, epsilon) * dr
        
        if abs((Ms - M_r) / Ms) < maximum and abs((Ls - L_r) / Ls) < maximum:
            good_surface = True
            break
        
        # If changes in M_r and L_r were too large, repeat with one-half the step size
        j += 1
        if j > j_max:
            print("Unable to converge in SUBROUTINE Surface --- Exiting")
            good_surface = False
            break
            
        dr /= 2
        r = rm + dr
    
    if not good_surface:
        print(f"The last values obtained by SUBROUTINE Surface were:")
        print(f"     M_r = {M_r:.6e}   dM_r/Ms = {(Ms - M_r)/Ms:.6e}")
        print(f"     L_r = {L_r:.6e}   dL_r/Ls = {(Ls - L_r)/Ls:.6e}")
    
    return r, P, T, M_r, L_r, rho, kappa, epsilon, good_surface

def core(M_r, L_r, P, T, X, Z, r):
    """
    Extrapolate core conditions based on the current zone conditions.
    
    Parameters:
    -----------
    M_r : float
        Mass enclosed within radius r
    L_r : float
        Luminosity at radius r
    P : float
        Pressure
    T : float
        Temperature
    X : float
        Hydrogen mass fraction
    Z : float
        Metals mass fraction
    r : float
        Radius
        
    Returns:
    --------
    P_0 : float
        Central pressure
    T_0 : float
        Central temperature
    rho_0 : float
        Central density
    kappa_0 : float
        Central opacity
    epsilon_0 : float
        Central energy generation rate
    rc_flag : str
        Radiative/convective flag
    dlnPdlnT : float
        Pressure gradient with respect to temperature
    good_T : bool
        True if temperature calculation successful
    """

    # Average density of the central core/sphere
    rho_0 = M_r / (four_pi_o3 * r**3)
    # Central pressure, Eq. (L.4)
    P_0 = P + (two_pi / 3) * G * rho_0**2 * r**2
    # Average energy generation rate of the central sphere
    epsilon_0 = L_r / M_r
    
    # Find core temperature by Newton-Raphson method (including radiation pressure)
    Y = comp.helium(X, Z)
    mu = comp.mean_molecular_weight(X, Y, Z)
    
    if rho_0 > 0:
        i = 0
        T_0 = T
        good_T = True
        
        # Convergence criteria
        converged = 1.0e-8
        i_max = 50

        # Newton-Raphson iteration to find core temperature
        while True:
            i += 1
            # Define f(T) = rho_0*k*T/(mu*m_H) + a_o3*T**4 - P_0
            f_T_0 = rho_0 * k_B * T_0 / (mu * m_H) + a_rad_o3 * T_0**4 - P_0
            # Define df/dT = rho_0*k/(mu*m_H) + 4*a_o3*T**3
            df_dT_0 = rho_0 * k_B / (mu * m_H) + 4 * a_rad_o3 * T_0**3
            
            dT = -f_T_0 / df_dT_0
            if abs(dT / T_0) < converged:
                break
                
            T_0 += dT
            if i > i_max:
                print("Unable to converge on core temperature in SUBROUTINE Core --- Exiting")
                good_T = False
                break
    else:
        T_0 = -T
        good_T = False
    
    if good_T:
        kappa_0 = physics.opacity(T_0, rho_0, X, Z)
        dlnPdlnT = physics.pt_gradient(P, P_0, T, T_0)
        gamma = physics.specific_heat_ratio()
        if dlnPdlnT < (gamma / (gamma - 1)):
            rc_flag = 'c'
        else:
            rc_flag = 'r'
    else:
        kappa_0 = -99.9
        dlnPdlnT = -99.9
        rc_flag = '*'
    
    return P_0, T_0, rho_0, kappa_0, epsilon_0, rc_flag, dlnPdlnT, good_T