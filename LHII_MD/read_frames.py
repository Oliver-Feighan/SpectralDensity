import numpy as np
import mdtraj


def load_dcd_file(dcd_file, top_file):
    traj = mdtraj.load_dcd(dcd_file, top=top_file)
    top = traj.top
    
    return traj, top

if __name__ == "__main__":
    traj, top = load_dcd_file("output/1ps_2fs_LHII.dcd", "LH2.prmtop")

    assert(len(traj) == 500)
    
    extracted_frames = traj[250::25]
    
    assert(len(extracted_frames) == 10)
    