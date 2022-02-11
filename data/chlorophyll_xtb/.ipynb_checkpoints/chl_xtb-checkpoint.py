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
    json_str = " -f json --schema none -s "
    
    json_run = subprocess.run(qcore_path + json_str + qcore_str,
                              shell=True,
                              stdout=subprocess.PIPE,
                              executable="/bin/bash",
                              universal_newlines=True)

    json_results = json.loads(json_run.stdout)

    return json_results


def run_single_chl(bcl):
    single_chl_str = f"\"res := xtb(structure(xyz = {bcl}) model='chlorophyll')\""    
        
    res = run_qcore(single_chl_str)

    return res["res"]["excitation_1_energy"]

def write_xyz(file_name, geom, symbols):
    assert(len(geom) == len(symbols))
    
    n_atoms = len(symbols)
    
    with open(file_name, 'w') as f:
        print(f"{n_atoms} \n", file=f)
        
        for i in range(n_atoms):
            print(f"{symbols[i]}\t{geom[i][0]}\t{geom[i][1]}\t{geom[i][2]}", file=f)

@timer
def run_trajectory(dcd_file, top_file, frames):
    
    traj, top = load_dcd_file(dcd_file, top_file)    

    atoms = list(top.atoms)
    
    bcl_atoms = top.select("resname =~ 'BCL'")
    bcl_atoms_list = bcl_atoms.tolist()
    
    symbols = [atoms[i].element.symbol for i in range(len(atoms))]

    mg_atoms = top.select("name =~ 'Mg*'").tolist()
    
    frames_xyz = traj.xyz[frames]
    
    for f, frame in enumerate(frames_xyz):
        for enum, mg_line in enumerate(mg_atoms):
            bcl = 10 * frame[mg_line:mg_line+140]
            
            write_xyz(f"bcl_{enum+1}.xyz", bcl, symbols[mg_line:mg_line+140])

        #exciton system
        structure_file_str = " ".join([f"structure(file = 'bcl_{enum+1}.xyz')" for enum, x in enumerate(mg_atoms)])
        
        exciton_str = f"\"res := excitons({structure_file_str} use_chlorophyll = true hamiltonian = 'states')\""
        
        start = time.time()
        res = run_qcore(exciton_str)
        print(f"Exciton, frame {frames[f]} ran in: ", time.time() - start)
        
        print(res["res"].keys())
        
    return 0

if __name__ == "__main__":
    print(os.environ["DCD_FILE"])
    print(os.environ["PRMTOP_FILE"])

    run_trajectory("../LHII_MD/output/1ps_2fs_LHII.dcd", "../LHII_MD/LH2.prmtop", [0])
