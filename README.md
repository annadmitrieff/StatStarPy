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

# Create a model of the Sun
sun = StellarModel(Msolar=1.0, Lsolar=1.0, Teff=5778, X=0.73, Z=0.02)

# Compute the model
sun.compute_model()

# Plot the structure
sun.plot_structure()
plt.show()

# Save the model data to a file
sun.save_model('sun_model.txt')
```

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
