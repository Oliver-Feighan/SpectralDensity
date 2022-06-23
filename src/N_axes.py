import numpy as np

import os
import openmm
import openmm.app

import mdtraj

import read_files

if __name__ == "__main__":
    print("running Na-Nc/Nb-Nd axis...")
    
    dcd_file, top_file = os.environ["DCD_FILE"], os.environ["TOP_FILE"]
   
    prmtop_top = openmm.app.AmberPrmtopFile(top_file)
    full_top = mdtraj.Topology.from_openmm(prmtop_top.topology)
    chl_indices = full_top.select("resname =~ 'BCL*'")
    chl_top = full_top.subset(chl_indices)

    traj = mdtraj.load_dcd(dcd_file, top=chl_top)
    top = traj.top
 
    bcl_atoms = top.select("resname =~ 'BCL'")
    bcl_atoms_list = bcl_atoms.tolist()
    
    bcl_traj = traj.atom_slice(bcl_atoms_list)
    
    N_indices = [5, 16, 25, 33]

    Nas = np.zeros((bcl_traj.n_frames, 27, 3))
    Nbs = np.zeros((bcl_traj.n_frames, 27, 3))
    Ncs = np.zeros((bcl_traj.n_frames, 27, 3))
    Nds = np.zeros((bcl_traj.n_frames, 27, 3))

    for f in range(bcl_traj.n_frames):
        for i in range(0, 3780, 140):

            Nas[f][i//140] = bcl_traj.xyz[f, N_indices[0]+i, :] * 10
            Nbs[f][i//140] = bcl_traj.xyz[f, N_indices[1]+i, :] * 10
            Ncs[f][i//140] = bcl_traj.xyz[f, N_indices[2]+i, :] * 10
            Nds[f][i//140] = bcl_traj.xyz[f, N_indices[3]+i, :] * 10
    
   
    np.save("Na_300ps_2fs", Nas)
    np.save("Nb_300ps_2fs", Nbs)
    np.save("Nc_300ps_2fs", Ncs)
    np.save("Nd_300ps_2fs", Nds)



 
