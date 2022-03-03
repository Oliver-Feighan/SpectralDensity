import numpy as np
import glob
import re

def read_splits(prop, shape, run, screened=False):
    if not screened:
        files = glob.glob(f"../data/chlorophyll_xtb/{run}_LHII_states_energies*9.npy")
    else:
        files = glob.glob(f"../data/chlorophyll_xtb/{run}_LHII_states_energies*9_screened.npy")

    starting_frames = [re.findall(r"\d+",x)[2] for x in files]
    starting_frames.sort()
    
    res = None

    for i in starting_frames:
        if not screened:
            split = np.load(f"../data/chlorophyll_xtb/{run}_LHII_{prop}_{int(i)}_{int(i)+999}.npy")
        else:
            split = np.load(f"../data/chlorophyll_xtb/{run}_LHII_{prop}_{int(i)}_{int(i)+999}_screened.npy")

        assert(split.shape == shape)

        if res is not None:
            res = np.concatenate((res, split))
        else:
            res = split
            
    return res


def thermal_energy(energies):
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