import numpy as np
import matplotlib.pyplot as plt


import read_files
import os



if __name__ == "__main__":
    print("running Na-Nc/Nb-Nd axis...")
    
    dcd_file, top_file = os.environ["DCD_FILE"], os.environ["TOP_FILE"]
    
    traj = read_files.load_dcd_file(dcd_file, top_file)
    top = traj.top
    
    bcl_atoms = top.select("resname =~ 'BCL'")
    bcl_atoms_list = bcl_atoms.tolist()
    
    bcl_traj = traj.atom_slice(bcl_atoms_list)
    
    N_indices = [5, 16, 25, 33]

    Nas = np.zeros((bcl_traj.n_frames, 3))
    Nbs = np.zeros((bcl_traj.n_frames, 3))
    Ncs = np.zeros((bcl_traj.n_frames, 3))
    Nds = np.zeros((bcl_traj.n_frames, 3))

    for f in range(bcl_traj.n_frames):
        for i in range(0, 3780, 140):

            Nas[f] = bcl_traj.xyz[f, N_indices[0]+i, :] * 10
            Nbs[f] = bcl_traj.xyz[f, N_indices[1]+i, :] * 10
            Ncs[f] = bcl_traj.xyz[f, N_indices[2]+i, :] * 10
            Nds[f] = bcl_traj.xyz[f, N_indices[3]+i, :] * 10
    
    
    plt.plot(np.linalg.norm(Nas-Nbs, axis=1))
    plt.show()