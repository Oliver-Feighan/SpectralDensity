import os

import matplotlib.pyplot as plt
import numpy as np
import scipy
import json

import hessian
import system
import chl_xtb

def q_coordinate(diff, symbols, mode, freq, disp):
    """
    unitless coordinate for moving along a normal mode. Given by equation:
    $ q_i = \sqrt{\frac{\omega_i}{\hbar}} x^m_i $
    
    where $\omega_i$ is the angular frequency of the oscillator (normal mode) and 
    $x_^m_i$ is the normal mode in mass weighted coordinates 
    """
    
    hbar = scipy.constants.hbar
    au_to_kg = scipy.constants.physical_constants["atomic mass constant"][0] 
    
    omega = freq * 1e15 * 2 * np.pi # fs^-1 * s^-1/fs^-1 * 2\pi
    
    masses = [system.masses()[symbols[i]] * au_to_kg for i in range(len(mode))]
    mass_weighted = np.array([diff[i] * 1e-10 * np.sqrt(masses[i]) for i in range(len(mode))]) #1e-10 for A to m conversion
    mass_weighted_disp = np.sign(disp) * np.linalg.norm(mass_weighted, axis=1)
    
    q = np.sum(np.sqrt(omega/hbar) * mass_weighted_disp)
    
    return q


def move_along_mode(origin, symbols, mode, freq, disp):
    #move along mode
    mode = mode * 0.529177 # bohr to angstrom
    norm = np.sum(np.linalg.norm(mode, axis=1))
    res = origin + disp * mode/norm
    
    
    #assert amount displacement is correct
    diff = origin - res
    total_diff = np.sum(np.linalg.norm(diff, axis=1))
    assert(abs(total_diff - abs(disp)) < 1e-6)
    
    
    #calculate q
    q = q_coordinate(diff, symbols, mode, freq, disp)
    
    return res, q

def write_xyz(xyz, symbols):
    coord_list = ([f"['{s}', {c[0]}, {c[1]}, {c[2]}]" for s, c in zip(symbols, xyz)])
    
    return f"[{', '.join(coord_list)}]"


def run_displacement(origin, symbols, mode, freq, disp):
    
    moved, q = move_along_mode(origin, symbols, mode, freq, disp)
    qcore_str = f"\"res := xtb(model='chlorophyll' structure(xyz = {write_xyz(moved, symbols)}))\""

    res = chl_xtb.run_qcore(qcore_str)
    
    ground_energy = res["res"]["energy"]
    transition_energy = res["res"]["excitation_1_energy"]
    excited_energy = ground_energy + transition_energy
    
    return ground_energy, transition_energy, excited_energy, q

def predict_min(domain, features):
    
    coeffs = np.polyfit(domain, features, 2)
    poly = np.poly1d(coeffs)
    crit = poly.deriv().r[0]
    
    return crit


def run_mode(origin, symbols, mode, freq):
    disps = np.linspace(-1, 1, 21)
    
    res = [run_displacement(xyz, symbols, mode, freq, d) for d in disps]
    
    ground = [r[0] for r in res]
    excited = [r[2] for r in res]
    
    qs = [r[3] for r in res]

    g_min = predict_min(qs, ground)
    e_min = predict_min(qs, excited)

    hrf = (e_min - g_min)**2/2
    
    return {
        "displacements" : disps.tolist(),
        "ground_energies" : ground,
        "transition_energies" : [r[1] for r in res],
        "excited_energies" : excited,
        "q_coords" : qs,
        "g_min" : g_min,
        "e_min" : e_min,
        "hrf" : hrf
    }
    

if __name__ == "__main__":
    print("Calculating Huang-Rhys parameters")
    
    molden_file = os.environ["MOLDEN_FILE"]
    xyz_file = os.environ["XYZ_FILE"]
    
    print(f"molden file : {molden_file}")
    print(f"xyz file : {xyz_file}")
    
    symbols, xyz = system.read_xyz(xyz_file)
    
    n_atoms = len(symbols)
    
    print(f"no. atoms : {n_atoms}")
    
    modes = hessian.get_modes(molden_file, 3*n_atoms-6, n_atoms)
    modes = modes.transpose(2, 0, 1)
    
    wavenumbers = hessian.get_wavenumbers(molden_file)
    frequencies = hessian.wavenumber_to_frequency(wavenumbers)
    
    results = {}
    
    for m, mode in enumerate(modes):
        if m > 5:
            break
            
        print(f"Mode {m}, v : {wavenumbers[m]:4.1f} cm^-1")

        if wavenumbers[m] < 0:
            print("    imaginary mode -> ignored")
            print("-"*5)
            continue
        
        mode_props = run_mode(xyz, symbols, mode, frequencies[m])
        
        mode_props["wavenumber"] = wavenumbers[m]
        mode_props["frequency"] = frequencies[m]
        
        results[f"mode_{m}"] = mode_props
        
        print(f"HRF : {mode_props['hrf']:1.4f}")
        print("-"*5)
        
    with open('huang_rhys.json', 'w') as fp:
        json.dump(results, fp)

    exit(0)
    
