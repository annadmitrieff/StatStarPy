# zone_quantities.py

import numpy as np
from .constants import dp

# Number of stellar structure equations
n = 4

# Previous zone data
Mm = 0.0
Lm = 0.0
rm = 0.0
Pm = 0.0
Tm = 0.0
Xm = 0.0
Zm = 0.0
rhom = 0.0
kappam = 0.0
taum = 0.0
epsilonm = 0.0

# Current zone data
rho = 0.0
kappa = 0.0
tau = 0.0
epsilon = 0.0
gamma = 0.0
dlnPdlnT = 0.0
rc_flag = 'r'

# Current step size flag
step_size_condition = 0

# The first derivatives from the stellar structure equations to be used by Runge Kutta routines
dfdr0 = np.zeros(n, dtype=dp)