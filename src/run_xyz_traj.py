import os
import numpy as np
import re
import time

import chl_xtb


def extract_xyz(lines):

    symbols = [re.findall(r"[A-Za-z]+", l)[0] for l in lines]

    coords = np.array([[float(x) for x in re.findall(r"-?\d+\.\d+", l)] for l in lines])

    string_lines = [f"[\'{s}\', {c[0]}, {c[1]}, {c[2]}]" for s,c in zip(symbols, coords)]

    return f"[{','.join(string_lines)}]"


def run_frame(xyz):

    res = chl_xtb.run_qcore(f"\"res := xtb(structure(xyz = {xyz}) model='chlorophyll')\"")

    dipole = res["res"]["excitation_1_transition_dipole"]
    energy = res["res"]["excitation_1_energy"]

    return energy, dipole

def run_trajectory(file, n_atoms, start, end):
    with open(file) as line_file:
        lines = line_file.readlines()

        frames = list(range(start, end))
        
        energies = np.zeros(len(frames))
        dipoles = np.zeros((len(frames), 3))
        
        for enum, f in enumerate(frames):
            frame_lines = lines[(n_atoms+2)*(f):(n_atoms+2)*(f+1)]

            xyz = extract_xyz(frame_lines[2:])
     
            energy, dipole = run_frame(xyz)
        
            energies[enum] = energy
            dipoles[enum] = dipole

    return energies, dipoles


if __name__ == "__main__":
    print("Run monomer MD chl-xtb")

    start = int(os.environ["FRAME_START"])
    end = int(os.environ["FRAME_END"])
    
    traj_file = os.environ["TRAJ_FILE"]
    n_atoms = int(os.environ["N_ATOMS"])
    
    start_time = time.time()
    energies, dipoles = run_trajectory("monomer_trajectory.xyz", 140, start, end+1)
    
    print(f"frames {start} to {end} run in : {time.time() - start_time}")

    energies_name = f"excitation_energies_{start}_{end}.npy"
    dipoles_name = f"transition_dipoles_energies_{start}_{end}.npy"
    
    np.save(energies_name, energies)
    np.save(dipoles_name, dipoles)
    
    exit(0)
