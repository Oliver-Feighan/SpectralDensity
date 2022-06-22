import numpy as np
import glob
import re
import mdtraj
import tqdm

import utils

@utils.timer
def load_dcd_file(dcd_file, top_file):
    return mdtraj.load_dcd(dcd_file, top=top_file)


def read_LHII_data(prop, shape, run, only_chl=False, screened=False):
    prefix = "chl_" if only_chl else ""
    
    if not screened:
        files = glob.glob(f"../data/chlorophyll_xtb/{run}/{prefix}{run}_LHII_states_energies*9.npy")
    else:
        files = glob.glob(f"../data/chlorophyll_xtb/{run}/{prefix}{run}_LHII_states_energies*9_screened.npy")

    starting_frames = [re.findall(r"\d+",x)[4] for x in files]
    starting_frames.sort()
    
    res = None

    for i in tqdm.tqdm(starting_frames):
        if not screened:
            split = np.load(f"../data/chlorophyll_xtb/{run}/{prefix}{run}_LHII_{prop}_{int(i)}_{int(i)+shape[0]-1}.npy")
        else:
            split = np.load(f"../data/chlorophyll_xtb/{run}/{prefix}{run}_LHII_{prop}_{int(i)}_{int(i)+shape[0]-1}_screened.npy")

        assert(split.shape == shape)

        if res is None:
            res = split
        else:
            res = np.concatenate((res, split))
            
    return res


def read_monomer_data():
    files = glob.glob("../data/monomer_MD/truncated_excitation_energies*npy")
    
    starts = [int(re.findall("\d+", x)[0]) for x in files]
    starts.sort()
    
    all_data = np.zeros(0)
    
    for i in starts:
        data = np.load(f"../data/monomer_MD/truncated_excitation_energies_{i}_{i+2499}.npy")
        
        all_data = np.concatenate((all_data, data))
        
    return all_data


def all_LHII_data(run, n, only_chl):
    hamils = read_LHII_data("hamiltonians", (n, 28, 28), run, only_chl)
    distances = read_LHII_data("distances", (n, 27, 27), run, only_chl)
    dipoles = read_LHII_data("transition_dipoles", (n, 27, 3), run, only_chl)

    eigvec = read_LHII_data("eigenvectors", (n, 28, 28), run, only_chl)
    
    eigvec = np.transpose(eigvec, (0, 2, 1))

    density = np.square(eigvec)

    eigval = read_LHII_data("states_energies", (n, 28), run, only_chl)

    for d, H, vec, val in tqdm.tqdm(zip(distances, hamils, eigvec, eigval)):
        assert(np.max(val - utils.reconstruct_val(vec, H)) < 1e-10)
        assert(utils.check_symmetric(H))

        assert(utils.check_symmetric(d))

    site_e = read_LHII_data("transition_energies", (n, 27), run, only_chl)
    
    exciton_energies = np.array([eigval[:,i] - eigval[:,0] for i in range(1, eigval.shape[1])]).T
    
    return {
        "hamils" : hamils,
        "distances" : distances,
        "dipoles" : dipoles,
        "eigvec" : eigvec,
        "eigval" : eigval,
        "site_e" : site_e,
        "exciton_energies" : exciton_energies
    }


def load_mg_xyz():
    traj = load_dcd_file("../LHII_MD/output/1ps_2fs_LHII.dcd", "../LHII_MD/LH2.prmtop")
    top = traj.top
    
    Mg_atom_indices = top.select("name =~ 'Mg*'")
    Mg_atoms = traj.atom_slice(Mg_atom_indices)
    
    return 10 * Mg_atoms.xyz # x10 for nm to A conversion