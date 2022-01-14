import analysis

import numpy as np
import mdtraj

import subprocess
import json
import os

def run_qcore(qcore_str):
    #qcore_path = os.environ["QCORE_PATH"]
    qcore_path="~/qcore/cmake-build-release/bin/qcore"
    json_str = " -n 1 -f json --schema none -s "
    
    json_run = subprocess.run(qcore_path + json_str + qcore_str,
                              shell=True,
                              stdout=subprocess.PIPE,
                              executable="/bin/bash",
                              universal_newlines=True)

    json_results = json.loads(json_run.stdout)

    return json_results


def read_bcls(traj, top, frames):
    atoms = list(top.atoms)
    
    bcl_atoms = top.select("resname =~ 'BCL'")
    bcl_atoms_list = bcl_atoms.tolist()
    
    symbols = [atoms[i].element.symbol for i in bcl_atoms_list]

    res_strs = []
    
    for frame in frames:
    
        bcl_xyz = traj.xyz[frame][bcl_atoms]

        n_atoms = len(bcl_xyz)

        frame_strs = ["[" for i in range(len(top.select("name =~ 'Mg*'")))]

        for i, coord in enumerate(bcl_xyz):
            bcl = i // 140

            if (i+1) % 140 != 0:
                # * 10 for conversion to angstrom from nanometers
                frame_strs[bcl] += f"['{symbols[i]}', {coord[0] * 10:3.6f}, {coord[1] * 10:3.6f}, {coord[2] * 10:3.6f}],"
            else:
                frame_strs[bcl] += f"['{symbols[i]}', {coord[0] * 10:3.6f}, {coord[1] * 10:3.6f}, {coord[2] * 10:3.6f}]"

        for i in range(len(frame_strs)):
            frame_strs[i] += "]"
            
        res_strs.append(frame_strs)

    return res_strs

@analysis.timer
def run_trajectory(dcd_file, top_file, frames):
    
    traj = analysis.load_dcd_file(dcd_file, top_file)
    top = traj.top
    
    bcls = read_bcls(traj, top, frames)
    
    results = np.zeros((len(bcls), len(bcls[0])))
    
    for f, frame in enumerate(bcls):
        for b, bcl in enumerate(frame):
            qcore_str = f"\"res := xtb(structure(xyz = {bcl}) model='chlorophyll')\""    

            res = run_qcore(qcore_str)

            results[f][b] = res["res"]["excitation_1_energy"]

    return results

if __name__ == "__main__":
    run_trajectory("../LHII_MD/output/1ps_2fs_LHII.dcd", "../LHII_MD/LH2.prmtop")