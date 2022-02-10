import numpy as np
import mdtraj


def load_dcd_file(dcd_file, top_file):
    traj = mdtraj.load_dcd(dcd_file, top=top_file)
    top = traj.top
    
    return traj, top

if __name__ == "__main__":
    traj, top = load_dcd_file("output/1ps_2fs_LHII.dcd", "LH2.prmtop")

    traj[400:].save("all_frames.pdb")
    traj[400::10].save("last_hundred_every_ten.pdb")
