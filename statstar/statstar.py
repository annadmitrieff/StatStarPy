# statstar.py

import numpy as np
import matplotlib.pyplot as plt
from .constants import dp, M_Sun, L_Sun, R_Sun
from . import user_io
from . import boundary_conditions as bc
from . import composition as comp
from . import physics
from . import zone_quantities as zq
from . import ode_integrator as ode
from . import stellar_structure as ss

"""
'main' file that puts it all together.
defines the class StellarModel and the methods that are used to compute the model given the inputs, which are defined in the class below.
"""


class StellarModel:
    """
    A class to generate and analyze a Zero Age Main Sequence (ZAMS) stellar model!
    """
    
    def __init__(self, Msolar=1.0, Lsolar=1.0, Teff=5778.0, X=0.73, Z=0.02):
        """
        Initialize the stellar model
        
        Parameters:
        -----------
        Msolar : float, optional
            Mass in solar units (default: 1.0)
        Lsolar : float, optional
            Luminosity in solar units (default: 1.0)
        Teff : float, optional
            Effective temperature in K (default: 5778.0)
        X : float, optional
            Hydrogen mass fraction (default: 0.73)
        Z : float, optional
            Metals mass fraction (default: 0.02)
        """
        # Store input parameters
        self.Msolar = Msolar
        self.Lsolar = Lsolar
        self.Teff = Teff
        self.X = X
        self.Z = Z
        
        # Calculate derived parameters
        self.Y = comp.helium(X, Z)
        self.Ms = Msolar * M_Sun
        self.Ls = Lsolar * L_Sun
        self.Rs = user_io.calculate_radius(self.Ls, self.Teff)
        self.Rsolar = self.Rs / R_Sun
        
        # Initialize results storage
        self.results = {
            'zone': [],
            'r': [],
            'tau': [],
            'M_r': [],
            'L_r': [],
            'T': [],
            'P': [],
            'rho': [],
            'kappa': [],
            'epsilon': [],
            'rc_flag': [],
            'dlnPdlnT': []
        }
        
        # Flags for model status
        self.model_computed = False
        self.ok_surface = False
        self.ok_core = False
    
    def compute_model(self):
        """
        Compute the stellar model by integrating inward from the surface
        
        Returns:
        --------
        bool
            True if model computation was successful
        """
        # Integration parameters
        dr_over_r = 1.0e-3  # Initial fractional step size
        M_fraction_limit = 0.01  # Mass fraction stop condition
        L_fraction_limit = 0.10  # Luminosity stop condition
        r_fraction_limit = 0.02  # radius stop condition
        maximum_zones = 10000  # Maximum number of zones allowed
        n_surface = 1  # Number of surface boundary zones
        
        # Reset zone quantities
        zq.Pm = 0.0
        zq.Tm = 0.0
        zq.Xm = self.X
        zq.Zm = self.Z
        zq.rm = self.Rs
        zq.taum = 0.0
        zq.rhom = 0.0
        zq.kappam = 0.0
        zq.epsilonm = 0.0
        zq.dlnPdlnT = 99.9
        zq.rc_flag = 'r'
        
        # Store initial conditions in results
        i = 0
        self.store_results(i, zq.rm, 0.0, 0.0, self.Ls, zq.Tm, zq.Pm, zq.rhom, 
                          zq.kappam, zq.epsilonm, zq.rc_flag, zq.dlnPdlnT)
        
        # Set up surface boundary conditions
        zq.Mm = self.Ms
        zq.Lm = self.Ls
        zq.rm = self.Rs
        dr = -dr_over_r * self.Rs
        zq.step_size_condition = 0
        
        # Compute surface zones
        for i in range(1, n_surface + 1):
            # Update last zone values if needed
            if i > 1:
                zq.Mm = M_r
                zq.Lm = L_r
                zq.rm = r
                zq.Pm = P
                zq.Tm = T
                zq.Xm = self.X
                zq.Zm = self.Z
                zq.taum = tau
                zq.rhom = rho
                zq.kappam = kappa
                zq.epsilonm = epsilon
            
            # Compute surface boundary values
            r, P, T, M_r, L_r, rho, kappa, epsilon, self.ok_surface = bc.surface(
                i, self.Ms, self.Ls, zq.rm, self.X, self.Z, dr)
            
            if not self.ok_surface:
                return False
            
            # Calculate optical depth
            tau = zq.taum + physics.optical_depth_change(kappa, zq.kappam, rho, zq.rhom, r, zq.rm)
            
            # Store results
            self.store_results(i, r, tau, 1 - M_r/self.Ms, L_r, T, P, rho, 
                              kappa, epsilon, zq.rc_flag, zq.dlnPdlnT)
        
        # If surface is ok, start main inward integration
        if self.ok_surface:
            # Load array of first derivatives to start the general inward integration
            Y = comp.helium(self.X, self.Z)
            mu = comp.mean_molecular_weight(self.X, Y, self.Z)
            zq.gamma = physics.specific_heat_ratio()
            zq.dlnPdlnT = physics.pt_gradient(zq.Pm, P, zq.Tm, T)
            
            zq.dfdr0[0] = ss.dp_dr(M_r, rho, r)
            zq.dfdr0[1] = ss.dm_dr(r, rho)
            zq.dfdr0[2] = ss.dl_dr(r, rho, epsilon)
            zq.dfdr0[3] = ss.dt_dr(kappa, rho, T, L_r, r, mu, M_r, zq.gamma, zq.dlnPdlnT)
            
            # Main inward integration loop
            ok_runge = True
            i_max = maximum_zones
            
            for i in range(i + 1, i_max + 1):
                # Update last zone values
                zq.Mm = M_r
                zq.Lm = L_r
                zq.Pm = P
                zq.Tm = T
                zq.Xm = self.X
                zq.Zm = self.Z
                zq.rm = r
                zq.taum = tau
                zq.rhom = rho
                zq.kappam = kappa
                zq.epsilonm = epsilon
                
                # Perform Runge-Kutta integration
                PMLT0 = np.array([zq.Pm, zq.Mm, zq.Lm, zq.Tm])
                
                def structure_function(j, x, y):
                    return ss.structure_eqns(j, x, y)
                
                PMLT, ok_runge = ode.rk_4(zq.n, zq.rm, dr, PMLT0, zq.dfdr0, structure_function)
                
                if not ok_runge:
                    break
                
                # Results from the current step
                P = PMLT[0]
                M_r = PMLT[1]
                L_r = PMLT[2]
                T = PMLT[3]
                
                # Calculate optical depth
                tau = zq.taum + physics.optical_depth_change(
                    zq.kappa, zq.kappam, zq.rho, zq.rhom, zq.rm + dr, zq.rm)
                
                # Store results
                self.store_results(i, r, tau, 1 - M_r/self.Ms, L_r, T, P, zq.rho,
                                  zq.kappa, zq.epsilon, zq.rc_flag, zq.dlnPdlnT)
                
                # Check termination conditions
                if ((M_r/self.Ms < M_fraction_limit and L_r/self.Ls < L_fraction_limit and r/self.Rs < r_fraction_limit)
                    or T < 0 or P < 0):
                    break
                
                # Is it time to change step size?
                if zq.step_size_condition == 0:
                    if M_r < 0.99 * self.Ms:
                        dr = -self.Rs / 100
                        zq.step_size_condition = 1
                elif zq.step_size_condition == 1:
                    if abs(dr) > 5 * r:
                        dr = dr / 10
                        zq.step_size_condition = 2
                
                # Update radius for next step
                r = r + dr
            
            # If RK integration was successful, extrapolate to the core
            if ok_runge:
                i += 1
                
                # Determine core conditions
                P_0, T_0, rho_0, kappa_0, epsilon_0, rc_flag, dlnPdlnT, self.ok_core = bc.core(
                    M_r, L_r, P, T, self.X, self.Z, r)
                
                if not self.ok_core:
                    print("WARNING: There was a problem with the core extrapolation routine") # :-(
                
                # Calculate optical depth to center
                tau = tau + physics.optical_depth_change(
                    kappa_0, zq.kappa, rho_0, zq.rho, 0.0, r)
                
                # Store core results
                self.store_results(i, 0.0, tau, 1 - M_r/self.Ms, L_r, T_0, P_0, rho_0,
                                  kappa_0, epsilon_0, rc_flag, dlnPdlnT)
            
            self.model_computed = True
            return True
        
        return False
    
    def store_results(self, i, r, tau, mass_frac, L_r, T, P, rho, kappa, epsilon, rc_flag, dlnPdlnT):
        """
        Store results from a single zone
        
        Parameters:
        -----------
        i : int
            Zone number
        r : float
            Radius
        tau : float
            Optical depth
        mass_frac : float
            1 - M_r/Ms
        L_r : float
            Luminosity at radius r
        T : float
            Temperature
        P : float
            Pressure
        rho : float
            Density
        kappa : float
            Opacity
        epsilon : float
            Energy generation rate
        rc_flag : str
            Radiative/convective flag
        dlnPdlnT : float
            Pressure gradient with respect to temperature
        """
        self.results['zone'].append(i)
        self.results['r'].append(r)
        self.results['tau'].append(tau)
        self.results['M_r'].append(self.Ms * (1.0 - mass_frac))
        self.results['L_r'].append(L_r)
        self.results['T'].append(T)
        self.results['P'].append(P)
        self.results['rho'].append(rho)
        self.results['kappa'].append(kappa)
        self.results['epsilon'].append(epsilon)
        self.results['rc_flag'].append(rc_flag)
        self.results['dlnPdlnT'].append(dlnPdlnT)
    
    def plot_structure(self, save_path=None):
        """
        Create a summary plot of the stellar structure
        
        Parameters:
        -----------
        save_path : str, optional
            Path to save the plot (if None, plot is displayed)
            
        Returns:
        --------
        matplotlib.figure.Figure
            The figure object
        """
        if not self.model_computed:
            print("Model has not been computed yet. Call compute_model() first.")
            return None
        
        # Convert to numpy arrays and normalize
        r_norm = np.array(self.results['r']) / self.Rs
        M_r_norm = np.array(self.results['M_r']) / self.Ms
        L_r_norm = np.array(self.results['L_r']) / self.Ls
        T_norm = np.array(self.results['T']) / np.max(self.results['T'])
        P_norm = np.array(self.results['P']) / np.max(self.results['P'])
        rho_norm = np.array(self.results['rho']) / np.max(self.results['rho'])
        
        # Create figure with 2 rows and 2 columns
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        fig.suptitle(f'Stellar Model: M = {self.Msolar:.2f} M$_\\odot$, '
                    f'L = {self.Lsolar:.2f} L$_\\odot$, '
                    f'T$_{{eff}}$ = {self.Teff:.0f} K, '
                    f'X = {self.X:.2f}, Z = {self.Z:.4f}',
                    fontsize=14)
        
        # Mass and Luminosity vs radius
        ax = axes[0, 0]
        ax.plot(r_norm, M_r_norm, 'b-', label='M$_r$/M')
        ax.plot(r_norm, L_r_norm, 'r-', label='L$_r$/L')
        ax.set_xlabel('r/R')
        ax.set_ylabel('Fractional Mass and Luminosity')
        ax.legend()
        ax.grid(True)
        
        # Temperature and Pressure vs radius
        ax = axes[0, 1]
        ax.plot(r_norm, T_norm, 'r-', label='T/T$_{max}$')
        ax.plot(r_norm, P_norm, 'b-', label='P/P$_{max}$')
        ax.set_xlabel('r/R')
        ax.set_ylabel('Normalized Temperature and Pressure')
        ax.legend()
        ax.grid(True)
        
        # Density vs radius
        ax = axes[1, 0]
        ax.plot(r_norm, rho_norm, 'k-')
        ax.set_xlabel('r/R')
        ax.set_ylabel('$\\rho$/$\\rho_{max}$')
        ax.grid(True)
        
        # Energy generation rate vs radius
        ax = axes[1, 1]
        eps_norm = np.array(self.results['epsilon']) / np.max(self.results['epsilon'])
        ax.plot(r_norm, eps_norm, 'g-')
        ax.set_xlabel('r/R')
        ax.set_ylabel('$\\epsilon$/$\\epsilon_{max}$')
        ax.grid(True)
        
        # Mark radiative and convective regions
        rc_flag = np.array(self.results['rc_flag'])
        r_rad = r_norm[rc_flag == 'r']
        r_conv = r_norm[rc_flag == 'c']
        
        for ax in axes.flatten():
            ymin, ymax = ax.get_ylim()
            ax.fill_between(r_rad, ymin, ymax, color='yellow', alpha=0.2, label='Radiative')
            ax.fill_between(r_conv, ymin, ymax, color='orange', alpha=0.2, label='Convective')
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path)
        
        return fig
    
    def save_model(self, filename='stellar_model.txt'):
        """
        Save the stellar model data to a text file
        
        Parameters:
        -----------
        filename : str, optional
            Filename to save results (default: 'stellar_model.txt')
        """
        if not self.model_computed:
            print("Model has not been computed yet. Call compute_model() first.")
            return
        
        with open(filename, 'w') as f:
            # Write header
            f.write(f"{'':<45}A ZAMS Stellar Model\n")
            f.write(f"{'':<45}--------------------\n\n")
            f.write(f"{'':<45}M    = {self.Msolar:11.5f} solar\n")
            f.write(f"{'':<45}L    = {self.Lsolar:11.5f} solar\n")
            f.write(f"{'':<45}R    = {self.Rsolar:11.5f} solar\n")
            f.write(f"{'':<45}Teff = {self.Teff:11.5f} K\n")
            f.write(f"{'':<45}X    = {self.X:11.5f}\n")
            f.write(f"{'':<45}Y    = {self.Y:11.5f}\n")
            f.write(f"{'':<45}Z    = {self.Z:11.5f}\n\n\n")
            
            # Write column headers
            f.write(" zone      r         tau     1-M_r/Ms      L_r         T          P         rho        kap        eps    dlnPdlnT\n")
            
            # Write data
            for i in range(len(self.results['zone'])):
                zone = self.results['zone'][i]
                r = self.results['r'][i]
                tau = self.results['tau'][i]
                mass_frac = 1.0 - self.results['M_r'][i] / self.Ms
                L_r = self.results['L_r'][i]
                T = self.results['T'][i]
                P = self.results['P'][i]
                rho = self.results['rho'][i]
                kappa = self.results['kappa'][i]
                epsilon = self.results['epsilon'][i]
                rc_flag = self.results['rc_flag'][i]
                dlnPdlnT = self.results['dlnPdlnT'][i]
                
                f.write(f"{zone:5d} {r:11.3e} {tau:11.3e} {mass_frac:11.3e} {L_r:11.3e} "
                       f"{T:11.3e} {P:11.3e} {rho:11.3e} {kappa:11.3e} {epsilon:11.3e} "
                       f"  {rc_flag} {dlnPdlnT:5.1f}\n")