import functools
import numpy as np
import mdtraj

import subprocess
import time
import json
import os

from multiprocessing import Pool

def timer(func):
        @functools.wraps(func)
        def wrapper_timer(*args, **kwargs):
            start_time = time.perf_counter()
            value = func(*args, **kwargs)
            end_time = time.perf_counter()
            run_time = end_time - start_time
            print(f"Finished {func.__name__!r} in {run_time:4f} secs")
            return value
        return wrapper_timer
    
@timer
def load_dcd_file(dcd_file, top_file):
    traj = mdtraj.load_dcd(dcd_file, top=top_file)
    top = traj.top
    
    return traj, top

def run_qcore(qcore_str):
    qcore_path = os.environ["QCORE_PATH"]
    #qcore_path="~/qcore/cmake-build-release/bin/qcore"
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

def run_single_chl(bcl):
    single_chl_str = f"\"res := xtb(structure(xyz = {bcl}) model='chlorophyll')\""    
        
    res = run_qcore(single_chl_str)

    return res["res"]["excitation_1_energy"]

@timer
def run_trajectory(dcd_file, top_file, frames):
    
    traj, top = load_dcd_file(dcd_file, top_file)    

    bcls = read_bcls(traj, top, frames)
    
    results = np.zeros((len(bcls), len(bcls[0])))
    
    for f, frame in enumerate(bcls):
        #single chlorophylls
        
        start = time.time()
        with Pool(20) as p:
            frame_energies = p.map(run_single_chl, frame)      

            results[f] = frame_energies
            
        end = time.time()
        print("Single ran in: ", end - start)
        
        #exciton system
        structures_str = " ".join([f"structure(xyz={bcl})" for bcl in frame])
        
        exciton_str = f"\"res := excitons({structures_str} use_chlorophyll = true hamiltonian = 'states')\""
        
        start = time.time()
        res = run_qcore(exciton_str)
        print("Exciton ran in: ", time.time() - start)
        
        print(res["res"].keys())
        print(len(res["res"]["eigenvalues"]))
        print([x - res["res"]["eigenvalues"][0] for x in res["res"]["eigenvalues"]])
            
        
            
    return results

if __name__ == "__main__":
    print(os.environ["DCD_FILE"])
    print(os.environ["PRMTOP_FILE"])

    run_trajectory("../LHII_MD/output/1ps_2fs_LHII.dcd", "../LHII_MD/LH2.prmtop", [0])
