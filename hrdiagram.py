from statstar import StellarModel
import matplotlib.pyplot as plt
import numpy as np
import time

masses = np.logspace(-0.3, 1.3, 15)  # Range from about 0.5 to 20 solar masses
effective_temps = []
luminosities = []
radii = []
successful_masses = []

colors = plt.cm.seismic_r(np.linspace(0, 1, len(masses)))

print("Computing stellar models...")
for i, mass in enumerate(masses):
    print(f"Calculating model for {mass:.2f} solar masses...")
    
    # For each mass, estimate appropriate starting values for luminosity and temperature
    # These are approximate relations for main sequence stars
    if mass < 1.0:
        lum_estimate = mass**3.5
        teff_estimate = 5800 * mass**0.55  # Cooler (lower mass)
    else:
        lum_estimate = mass**3.9
        teff_estimate = 5800 * mass**0.60  # Hotter (higher mass)
    
    # Create model
    model = StellarModel(Msolar=mass, Lsolar=lum_estimate, Teff=teff_estimate, X=0.73, Z=0.02)
    
    # Some models fail to converge, so use try/except
    try:
        if model.compute_model():
            effective_temps.append(model.Teff)
            luminosities.append(model.Lsolar)
            radii.append(model.Rsolar)
            successful_masses.append(mass)
            print(f"  Success! Teff={model.Teff:.0f}K, L={model.Lsolar:.2f}Lsun")
        else:
            print(f"  Model computation failed for {mass:.2f} solar masses")
    except Exception as e:
        print(f"  Error for {mass:.2f} solar masses: {e}")
    
    time.sleep(0.1)

# Diagram setup:
plt.figure(figsize=(10, 8))

scatter = plt.scatter(effective_temps, luminosities, 
                     c=successful_masses, 
                     cmap='seismic_r', # tried to choose one similar to HR diag.
                     s=[r*20 for r in radii],  # Size points by radius
                     alpha=0.8,
                     edgecolor='black',
                     linewidth=0.5)

cbar = plt.colorbar(scatter)
cbar.set_label('Stellar Mass (solar masses)', fontsize=12)

# Note: temperature axis is traditionally reversed
plt.xscale('log')
plt.yscale('log')
plt.gca().invert_xaxis()

# Labels etc.
plt.xlabel('Effective Temperature (K)', fontsize=14)
plt.ylabel('Luminosity (L$_\\odot$)', fontsize=14)
plt.title('Hertzsprung-Russell Diagram', fontsize=16)

# Set axis limits to match traditional H-R diagram ranges
plt.xlim(40000, 2500)  # From hot O stars to cool M stars
plt.ylim(0.01, 1E6)    # From dim to very bright stars

plt.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('hr_diagram.png', dpi=300)
plt.show()

