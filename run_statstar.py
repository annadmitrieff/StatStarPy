from statstar import StellarModel
import matplotlib.pyplot as plt
import numpy as np

# Example: A 2 solar mass model
model = StellarModel(Msolar=1.0, Lsolar=11.0, Teff=5800, X=0.73, Z=0.02)

if model.compute_model():
    # Extract normalized radius and temperature profile...
    r_norm = np.array(model.results['r']) / model.Rs
    T = np.array(model.results['T'])
    
    # Plot temperature profile!
    plt.figure(figsize=(10, 6))
    plt.plot(r_norm, T, 'r-', linewidth=2)
    plt.xlabel('r/R', fontsize=12)
    plt.ylabel('Temperature (K)', fontsize=12)
    plt.title('Temperature Profile for a Solar Mass Star', fontsize=14)
    plt.grid(True)
    plt.yscale('log')
    plt.show()