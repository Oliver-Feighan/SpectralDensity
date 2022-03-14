import functools
import openmm
import openmm.app
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
def load_dcd_file(dcd_file, top_file, only_chl):
    if only_chl:
        prmtop_top = openmm.app.AmberPrmtopFile(top_file)
        full_top = mdtraj.Topology.from_openmm(prmtop_top.topology)
        chl_indices = full_top.select("resname =~ 'BCL*'")
        chl_top = full_top.subset(self._chl_indices)

        traj = mdtraj.load_dcd(dcd_file, top=chl_top)
        top = traj.top
        
        return traj, top
    
    else:
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
def run_trajectory(dcd_file, top_file, frames, only_chl):
    
    traj, top = load_dcd_file(dcd_file, top_file, only_chl)    

    basename = pathlib.Path(dcd_file).name
    run = basename.replace(".dcd", "")

    atoms = list(top.atoms)
    
    bcl_atoms = top.select("resname =~ 'BCL'")
    bcl_atoms_list = bcl_atoms.tolist()
    
    symbols = [atoms[i].element.symbol for i in range(len(atoms))]

    mg_atoms = top.select("name =~ 'Mg*'").tolist()
    
    frames_xyz = traj.xyz[frames]
    
    transition_energies = np.empty((len(frames), len(mg_atoms)))
    transition_dipoles = np.empty((len(frames), len(mg_atoms), 3))
    state_energies = np.empty((len(frames), len(mg_atoms)+1))

    distances = np.empty((len(frames), len(mg_atoms), len(mg_atoms)))
    eigenvectors = np.empty((len(frames), len(mg_atoms)+1, len(mg_atoms)+1))
    hamiltonians = np.empty((len(frames), len(mg_atoms)+1, len(mg_atoms)+1))

    for f, frame in enumerate(frames_xyz):
        for enum, mg_line in enumerate(mg_atoms):
            bcl = 10 * frame[mg_line:mg_line+140]
            
            write_xyz(f"{run}_bcl_{enum+1}_{frames[0]}.xyz", bcl, symbols[mg_line:mg_line+140])

        #exciton system
        structure_file_str = " ".join([f"structure(file = '{run}_bcl_{enum+1}_{frames[0]}.xyz')" for enum, x in enumerate(mg_atoms)])
    
        exciton_str = f"\"res := excitons({structure_file_str} use_chlorophyll = true hamiltonian = 'states')\""
        
        start = time.time()
        res = run_qcore(exciton_str)
        print(f"Exciton, frame {frames[f]} ran in: ", time.time() - start)
        
        transition_energies[f] = read_hex_data(res["res"], "transition_energies")
        transition_dipoles[f] = read_hex_data(res["res"], "transition_dipoles").T
        state_energies[f] = read_hex_data(res["res"], "eigenvalues")
        distances[f] = read_hex_data(res["res"], "distances")
        eigenvectors[f] = read_hex_data(res["res"], "eigenvectors")
        hamiltonians[f] = read_hex_data(res["res"], "hamiltonian")

    frame_start=frames[0]
    frame_end  =frames[-1]

    transition_energies_name = basename.replace(".dcd", f"_transition_energies_{frame_start}_{frame_end}.npy")
    transition_dipoles_name = basename.replace(".dcd", f"_transition_dipoles_{frame_start}_{frame_end}.npy")
    states_energies_name = basename.replace(".dcd", f"_states_energies_{frame_start}_{frame_end}.npy")
    distances_name = basename.replace(".dcd", f"_distances_{frame_start}_{frame_end}.npy")
    eigenvectors_name = basename.replace(".dcd", f"_eigenvectors_{frame_start}_{frame_end}.npy")
    hamiltonians_name = basename.replace(".dcd", f"_hamiltonians_{frame_start}_{frame_end}.npy")

    np.save(transition_energies_name, transition_energies)
    np.save(transition_dipoles_name, transition_dipoles)
    np.save(states_energies_name, state_energies)
    np.save(distances_name, distances)
    np.save(eigenvectors_name, eigenvectors) 
    np.save(hamiltonians_name, hamiltonians)
    return 0

if __name__ == "__main__":
    dcd_file = os.environ["DCD_FILE"]
    prmtop_file = os.environ["PRMTOP_FILE"]
    only_chl = os.environ["ONLY_CHL"] if "ONLY_CHL" in os.inviron else False

    frame_start = int(os.environ["FRAME_START"])
    frame_end = int(os.environ["FRAME_END"])

    run_trajectory(dcd_file, prmtop_file, list(range(frame_start, frame_end)), only_chl)
