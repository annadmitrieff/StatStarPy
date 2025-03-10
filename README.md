# StatStar

A Python package for computing the internal structure of a zero-age main sequence (ZAMS) star.

## Description

StatStar is a Python translation of the Fortran code from "An Introduction to Modern Astrophysics" by Carroll & Ostlie, translated by me (Annie Dmitrieff). It performs inward integration of the stellar structure equations to create a model of a star's interior, given basic parameters like mass, luminosity, and composition.

The code includes:
- Hydrostatic equilibrium calculations
- Energy transport (radiative and convective)
- Nuclear energy generation (PP-chain, CNO cycle, and helium burning)
- Opacity calculations with various contributions
- Inward integration with adaptive step sizes

## Installation

Will be supported shortly...
```
pip install statstar
```

## Usage

Basic example:

```python
from statstar import StellarModel
import matplotlib.pyplot as plt
import numpy as np

# Example: A 1 solar mass model
model = StellarModel(Msolar=1.0, Lsolar=11.0, Teff=5800, X=0.73, Z=0.02)

if model.compute_model():
    # Extract normalized radius and temperature profile...
    r_norm = np.array(model.results['r']) / model.Rs
    T = np.array(model.results['T'])
    rho = np.array(model.results['rho'])
    P = np.array(model.results['P'])
    mass_fraction = 1 - np.array(model.results['M_r'])
    
    # Plot temperature profile!
    plt.figure(figsize=(10, 6))
    plt.plot(r_norm, T, 'r-', linewidth=2)
    plt.xlabel('r/R', fontsize=12)
    plt.ylabel('Temperature (K)', fontsize=12)
    plt.title('Temperature Profile for a Solar Mass Star', fontsize=14)
    plt.grid(True)
    plt.yscale('log')
    plt.show()

    # Plot density profile!
    plt.figure(figsize=(10, 6))
    plt.plot(r_norm, rho, 'b-', linewidth=2)
    plt.xlabel('r/R', fontsize=12)
    plt.ylabel('Density (kg/m^3)', fontsize=12)
    plt.title('Density Profile for a Solar Mass Star', fontsize=14)
    plt.grid(True)
    plt.yscale('log')
    plt.show()

    # Plot pressure profile!
    plt.figure(figsize=(10, 6))
    plt.plot(r_norm, P, 'g-', linewidth=2)
    plt.xlabel('r/R', fontsize=12)
    plt.ylabel('Pressure (N/m^2)', fontsize=12)
    plt.title('Pressure Profile for a Solar Mass Star', fontsize=14)
    plt.grid(True)
    plt.yscale('log')
    plt.show()
    
    # Plot mass profile!
    plt.figure(figsize=(10, 6))
    plt.plot(r_norm, mass_fraction, 'm-', linewidth=2)
    plt.xlabel('r/R', fontsize=12)
    plt.ylabel('Mass Fraction (1 - Mr/Ms)', fontsize=12)
    plt.title('Mass Profile for a Solar Mass Star', fontsize=14)
    plt.grid(True)
    plt.yscale('log')
    plt.show()
```
![DensProf](https://github.com/user-attachments/assets/7010c8c7-a6fc-40cd-a187-bf707c51c394)
![MassProf](https://github.com/user-attachments/assets/4cf5b0ad-e7dc-4473-8d55-afe90fad4b11)
![PressureProf](https://github.com/user-attachments/assets/84dace1c-b14b-45ed-bf51-9b6a637927a9)
![TempProf](https://github.com/user-attachments/assets/b9f81676-bfbf-4f37-b13e-42bc3c773ddc)

## Requirements

- numpy
- scipy
- matplotlib

## Limitations

This is an educational code with simplified physics. It uses:
- Kramers opacity laws
- Adiabatic convection
- Basic boundary conditions

It's not a research-grade stellar structure code but helps understand the basic concepts of stellar modeling.

## License

MIT License

## Author
Annie Dmitrieff (annadmitrieff at uga dot edu)
