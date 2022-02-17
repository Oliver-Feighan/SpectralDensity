import analysis

if __name__ == "__main__":
    traj = analysis.load_dcd_file("../../LHII_MD/output/1ps_2fs_LHII.dcd", "../../LHII_MD/LH2.prmtop")
    top = traj.top
    
    chl_indices = top.select("resname =~ 'BCL*'")
    chl_atoms = traj.atom_slice(chl_indices)
    
    chl_atoms.save("chl_1ps_2fs.pdb")