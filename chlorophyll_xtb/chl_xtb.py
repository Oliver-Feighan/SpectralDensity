import functools
import numpy as np
import mdtraj

import pathlib
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
    json_str = " -f json -s "
 
    json_run = subprocess.run(qcore_path + json_str + qcore_str,
                              shell=True,
                              stdout=subprocess.PIPE,
                              executable="/bin/bash",
                              universal_newlines=True)

    try:
        json_results = json.loads(json_run.stdout)
    except:
        debug_run = subprocess.run(qcore_path + " -s " + qcore_str,
                              shell=True,
                              stdout=subprocess.PIPE,
                              executable="/bin/bash",
                              universal_newlines=True)

        print(debug_run.stdout)
        exit(1) 
    else:
        return json_results

def write_xyz(file_name, geom, symbols):
    assert(len(geom) == len(symbols))
    
    n_atoms = len(symbols)
    
    open_file = open(file_name, 'w')
    print(f"{n_atoms} \n", file=open_file)
 
    for i in range(n_atoms):
        print(f"{symbols[i]}\t{geom[i][0]}\t{geom[i][1]}\t{geom[i][2]}", file=open_file)

    open_file.close()

def read_hex_data(data, key):
	
	hex_data = data[key]["data"]["bytes"]
	byte_data = bytes.fromhex(" ".join([format(n, "02x") for n in hex_data]))

	dtype = data[key]["dtype"]
	shape = data[key]["shape"]

	return np.frombuffer(byte_data, dtype).reshape(shape)

@timer
def run_trajectory(dcd_file, top_file, frames):
    
    traj, top = load_dcd_file(dcd_file, top_file)    

    atoms = list(top.atoms)
    
    bcl_atoms = top.select("resname =~ 'BCL'")
    bcl_atoms_list = bcl_atoms.tolist()
    
    symbols = [atoms[i].element.symbol for i in range(len(atoms))]

    mg_atoms = top.select("name =~ 'Mg*'").tolist()
    
    frames_xyz = traj.xyz[frames]
    
    transition_energies = np.empty((len(frames), len(mg_atoms)))
    state_energies = np.empty((len(frames), len(mg_atoms)+1))

    for f, frame in enumerate(frames_xyz):
        for enum, mg_line in enumerate(mg_atoms):
            bcl = 10 * frame[mg_line:mg_line+140]
            
            write_xyz(f"bcl_{enum+1}.xyz", bcl, symbols[mg_line:mg_line+140])

        #exciton system
        structure_file_str = " ".join([f"structure(file = 'bcl_{enum+1}.xyz')" for enum, x in enumerate(mg_atoms)])
    
        print(structure_file_str)
    
        exciton_str = f"\"res := excitons({structure_file_str} use_chlorophyll = true hamiltonian = 'states')\""
        
        start = time.time()
        res = run_qcore(exciton_str)
        print(f"Exciton, frame {frames[f]} ran in: ", time.time() - start)
        
        transition_energies[f] = read_hex_data(res["res"], "transition_energies")
        state_energies[f] = read_hex_data(res["res"], "eigenvalues")

    basename = pathlib.Path(dcd_file).name

    transition_energies_name = basename.replace(".dcd", "_transition_energies.npy")
    states_energies_name = basename.replace(".dcd", "_states_energies.npy")

    np.save(transition_energies_name, transition_energies)
    np.save(states_energies_name, state_energies)
        
    return 0

if __name__ == "__main__":
    dcd_file = os.environ["DCD_FILE"]
    prmtop_file = os.environ["PRMTOP_FILE"]

    run_trajectory(dcd_file, prmtop_file, list(range(200)))
