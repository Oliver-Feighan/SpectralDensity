import os
import numpy as np
import re
import time
import mdtraj

import chl_xtb
import read_files


def extract_xyz(symbols, coords):

    string_lines = [f"[\'{s}\', {c[0]}, {c[1]}, {c[2]}]" for s,c in zip(symbols, coords)]

    return f"[{','.join(string_lines)}]"


def run_frame(xyz):

    res = chl_xtb.run_qcore(f"\"res := xtb(structure(xyz = {xyz}) model='chlorophyll')\"")

    dipole = res["res"]["excitation_1_transition_dipole"]
    energy = res["res"]["excitation_1_energy"]

    return energy, dipole

def run_trajectory(dcd_file, top_file, start, end):
    
    traj = read_files.load_dcd_file(dcd_file, top_file)
    top = traj.top
    
    bcl_indices = top.select("resname =~ 'BCL'")

    bcl_atoms = traj.atom_slice(bcl_indices)
    
    symbols = [a.element.symbol for a in bcl_atoms.top.atoms]
    
    frames = list(range(start, end))

    energies = np.zeros(len(frames))
    dipoles = np.zeros((len(frames), 3))

    for enum, f in enumerate(frames):
        xyz = extract_xyz(symbols, bcl_atoms.xyz[f] * 10) #x10 for nm to A conversion

        print(xyz)
        
        energy, dipole = run_frame(xyz)

        energies[enum] = energy
        dipoles[enum] = dipole

    return energies, dipoles


if __name__ == "__main__":
    print("Run Ether MD chl-xtb")
    
    dcd_file = os.environ["DCD_FILE"]
    top_file = os.environ["PRMTOP_FILE"]
    
    start = int(os.environ["FRAME_START"])
    end = int(os.environ["FRAME_END"])
    
    start_time = time.time()
    energies, dipoles = run_trajectory(dcd_file, top_file, start, end+1)
    
    print(f"frames {start} to {end} run in : {time.time() - start_time}")

    energies_name = f"excitation_energies_{start}_{end}.npy"
    dipoles_name = f"transition_dipoles_energies_{start}_{end}.npy"
    
    np.save(energies_name, energies)
    np.save(dipoles_name, dipoles)
    
    exit(0)
