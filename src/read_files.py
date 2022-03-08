import numpy as np
import glob
import re
import mdtraj

import utils

@utils.timer
def load_dcd_file(dcd_file, top_file):
    return mdtraj.load_dcd(dcd_file, top=top_file)


def read_LHII_data(prop, shape, run, screened=False):
    if not screened:
        files = glob.glob(f"../data/chlorophyll_xtb/{run}_LHII_states_energies*9.npy")
    else:
        files = glob.glob(f"../data/chlorophyll_xtb/{run}_LHII_states_energies*9_screened.npy")

    starting_frames = [re.findall(r"\d+",x)[2] for x in files]
    starting_frames.sort()
    
    res = np.zeros(0)

    for i in starting_frames:
        if not screened:
            split = np.load(f"../data/chlorophyll_xtb/{run}_LHII_{prop}_{int(i)}_{int(i)+999}.npy")
        else:
            split = np.load(f"../data/chlorophyll_xtb/{run}_LHII_{prop}_{int(i)}_{int(i)+999}_screened.npy")

        assert(split.shape == shape)

        res = np.concatenate((res, split))
            
    return res


def read_monomer_data():
    files = glob.glob("../data/monomer_MD/excitation_energies*npy")
    
    starts = [int(re.findall("\d+", x)[0]) for x in files]
    starts.sort()
    
    all_data = np.zeros(0)
    
    for i in starts:
        data = np.load(f"../data/monomer_MD/excitation_energies_{i}_{i+2499}.npy")
        
        all_data = np.concatenate((all_data, data))
        
    return all_data


def all_LHII_data(run):
    hamils = read_LHII_data("hamiltonians", (1000, 28, 28), run)
    distances = read_LHII_data("distances", (1000, 27, 27), run)
    dipoles = read_LHII_data("transition_dipoles", (1000, 27, 3), run)

    eigvec = read_LHII_data("eigenvectors", (1000, 28, 28), run)

    eigvec = np.transpose(eigvec, (0, 2, 1))

    density = np.square(eigvec)

    eigval = read_LHII_data("states_energies", (1000, 28), run)

    for d, H, vec, val in zip(distances, hamils, eigvec, eigval):
        assert(np.max(val - utils.reconstruct_val(vec, H)) < 1e-10)
        assert(utils.check_symmetric(H))

        assert(utils.check_symmetric(d))

    site_e = read_LHII_data("transition_energies", (1000, 27), run)
    
    return {
        "hamils" : hamils,
        "distances" : distances,
        "dipoles" : dipoles,
        "eigvec" : eigvec,
        "eigval" : eigval,
        "site_e" : site_e,
    }


def load_mg_xyz():
    traj = load_dcd_file("../LHII_MD/output/1ps_2fs_LHII.dcd", "../LHII_MD/LH2.prmtop")
    top = traj.top
    
    Mg_atom_indices = top.select("name =~ 'Mg*'")
    Mg_atoms = traj.atom_slice(Mg_atom_indices)
    
    return 10 * Mg_atoms.xyz # x10 for nm to A conversion