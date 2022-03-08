import numpy as np

def thermal_energy(energies):
    """
    partition function
    """
    boltzmann_eV = 8.617333262145e-5 # eV / K
    hr_to_eV = 27.2114 # eV / hr

    boltzmann_hr = boltzmann_eV / hr_to_eV

    T = 300 # K 

    beta = 1 / (boltzmann_hr * T)

    E = energies

    boltz_factor = np.exp(-beta * E)

    Z = np.sum(boltz_factor, axis=1)

    w = boltz_factor / Z

    U = np.sum(E * w, axis=0)
    
    return U