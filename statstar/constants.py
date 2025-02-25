# constants.py

# All constants are either universal physical constants or taken from statstar's `constants.f90`.
# Instead of manual inputs, we will be using scipy and numpy :-) (21st century style)

import numpy as np
import scipy.constants as sc

# Precision parameters
dp = np.float64

# Pi and related constants
pi = np.pi
two_pi = 2 * pi
four_pi = 4 * pi
four_pi_o3 = four_pi / 3
pi_over_2 = pi / 2

# Physical constants
G = sc.gravitational_constant  # 6.673e-11
c = sc.speed_of_light  # 2.99792458e8
mu_0 = sc.mu_0
epsilon_0 = sc.epsilon_0

# Electron charge and energy units
e_C = sc.elementary_charge  # 1.602176462e-19
eV = e_C
keV = eV * 1.0e3
MeV = eV * 1.0e6
GeV = eV * 1.0e9

# Planck constant
h = sc.Planck  # 6.62606876e-34
hbar = h / two_pi

# Boltzmann constant
k_B = sc.Boltzmann  # 1.3806503e-23

# Stefan-Boltzmann and radiation constants
sigma = sc.Stefan_Boltzmann  # 5.670367e-8
a_rad = 4 * sigma / c
a_rad_o3 = a_rad / 3
four_ac_o3 = 4 * a_rad_o3 * c

# Masses
m_e = sc.electron_mass  # 9.10938188e-31
m_p = sc.proton_mass  # 1.67262158e-27
m_n = sc.neutron_mass  # 1.67492716e-27
m_H = 1.673532499e-27  # Hydrogen mass
u = sc.atomic_mass  # 1.66053873e-27

# Avogadro number and gas constant
N_A = sc.Avogadro  # 6.02214199e23
R_gas = sc.gas_constant  # 8.314472

# Atomic constants
a_0 = 4 * pi * epsilon_0 * hbar**2 / (m_e * e_C**2)  # Bohr radius
R_infty = m_e * e_C**4 / (64 * pi**3 * epsilon_0**2 * hbar**3 * c)  # Rydberg constant
R_H = m_p / (m_e + m_p) * R_infty

# Time constants
hr = 3600
day = 24 * hr
J_yr = 365.25 * day
yr = 3.15581450e7
T_yr = 3.155692519e7
G_yr = 3.1556952e7

# Astronomical length constants
AU = 1.4959787066e11
pc = 206264.806 * AU
ly = c * J_yr

# Solar constants
M_Sun = 1.9891e30
S_Sun = 1.365e3
L_Sun = four_pi * AU**2 * S_Sun
R_Sun = 6.95508e8
Te_Sun = (L_Sun / (four_pi * R_Sun**2 * sigma))**0.25

# Solar magnitudes
Mbol_Sun = 4.74
MU_Sun = 5.67
MB_Sun = 5.47
MV_Sun = 4.82
Mbol_Sun_ap = -26.83
MU_Sun_ap = -25.91
MB_Sun_ap = -26.10
MV_Sun_ap = -26.75
BC_Sun = -0.08

# Earth constants
M_Earth = 5.9736e24
R_Earth = 6.378136e6

# Unit conversions
cm = 1e-2
gram = 1e-3
erg = 1e-7
dyne = 1e-5
esu = 3.335640952e-10
statvolt = 2.997924580e2
gauss = 1e-4
angstrom = 1e-10
jansky = 1e-26