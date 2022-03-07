import numpy as np
import mdtraj


def load_dcd_file(dcd_file, top_file):
    traj = mdtraj.load_dcd(dcd_file, top=top_file)
    top = traj.top
    
    return traj, top

if __name__ == "__main__":
    traj, top = load_dcd_file("output/1000ps_100fs_LHII.dcd", "LH2.prmtop")

    assert(len(traj) == 10000)

    selected_frames = traj[5000::10]

    assert(len(selected_frames) == 500)

    selected_frames.save("LHII_frames.pdb")
