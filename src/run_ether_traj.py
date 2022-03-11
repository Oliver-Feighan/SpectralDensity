import os
import numpy as np
import re
import time
import mdtraj

import chl_xtb


def extract_xyz(symbols, coords):

    string_lines = [f"[\'{s}\', {c[0]}, {c[1]}, {c[2]}]" for s,c in zip(symbols, coords)]

    return f"[{','.join(string_lines)}]"


def run_frame(xyz):

    res = chl_xtb.run_qcore(f"\"res := xtb(structure(xyz = {xyz}) model='chlorophyll')\"")

    dipole = res["res"]["excitation_1_transition_dipole"]
    energy = res["res"]["excitation_1_energy"]

    return energy, dipole

def run_trajectory(dcd_file, pdb_file, start, end):
    
    pdb = mdtraj.load_pdb(pdb_file)
    traj = mdtraj.load_dcd(dcd_file, pdb.top)
    
    cla_indices = pdb.top.select("resname =~ 'CLA'")

    cla_atoms = traj.atom_slice(cla_indices)
    
    symbols = [a.element.symbol for a in cla_atoms.top.atoms]
    
    frames = list(range(start, end))

    energies = np.zeros(len(frames))
    dipoles = np.zeros((len(frames), 3))

    for enum, f in enumerate(frames):
        xyz = extract_xyz(symbols, cla_atoms.xyz[f] * 10) #x10 for nm to A conversion

        energy, dipole = run_frame(xyz)

        energies[enum] = energy
        dipoles[enum] = dipole

    return energies, dipoles


if __name__ == "__main__":
    print("Run Ether MD chl-xtb")
    
    dcd_file = os.environ["DCD_FILE"]
    pdb_file = os.environ["PDB_FILE"]
    
    start = int(os.environ["FRAME_START"])
    end = int(os.environ["FRAME_END"])
    
    start_time = time.time()
    energies, dipoles = run_trajectory(dcd_file, pdb_file, start, end+1)
    
    print(f"frames {start} to {end} run in : {time.time() - start_time}")

    energies_name = f"excitation_energies_{start}_{end}.npy"
    dipoles_name = f"transition_dipoles_energies_{start}_{end}.npy"
    
    np.save(energies_name, energies)
    np.save(dipoles_name, dipoles)
    
    exit(0)
