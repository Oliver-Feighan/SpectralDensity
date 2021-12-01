import numpy as np
import re

def read_hex_data(result, results_name, variable):
    
    hex_data = result[results_name][variable]["data"]["bytes"]
    byte_data = bytes.fromhex(" ".join([format(n, "02x") for n in hex_data]))
                                       
    dtype = result[results_name][variable]["dtype"]
    shape = result[results_name][variable]["shape"]
    
    return np.frombuffer(byte_data, dtype=dtype).reshape(shape)

def read_indexed_values(result, results_name, variable):
    
    data = result[results_name][variable]
    
    return data[0]

def get_wavenumbers(file_name):
    lines = list(open(file_name))
    
    break_points = []
    
    for enum, line in enumerate(lines):
        if "[" in line:
            break_points.append(enum)
            
    frequencies = [float(re.findall(r'-?\d+\.\d+', freq_line)[0]) for freq_line in lines[break_points[1]+1:break_points[2]-1]]
    
    return np.array(frequencies)
    

def get_mode(file_name, mode, n_atoms):
    norm_coord = []
    
    lines = list(open(file_name))
        
    for enum, line in enumerate(lines):
        if f"vibration {mode}\n" in line:
            for mode_line in range(enum+1, enum+1+n_atoms):
                norm_coord.append([float(x) for x in re.findall(r'-?\d+\.\d+', lines[mode_line])])
    
    return np.array(norm_coord)
                
def wavenumber_to_frequency(wavenumber):
    c = 2.99792458e10 # cm s-1
    
    frequency = wavenumber * c # cm-1 * cm s-1
    
    frequency_fs = frequency * 1e-15
    
    return frequency_fs

def read_symbols(file_name):
    symbols = []
    
    with open(file_name, 'r') as f:
        lines = list(f.readlines())
        
        for line in lines[2:]:
            symbol = re.findall(r'[A-Za-z]+', line)
            
            if len(symbol) != 0:
                symbols.append(symbol[0])
                
    return symbols

def get_all_displacements(molden_file, n_atoms, n_modes):
    assert(n_modes == 3 * n_atoms - 6)

    all_displacements = np.zeros((n_atoms, n_modes))

    for i in range(n_modes):
        mode_i = get_mode(molden_file, i+1, n_atoms)

        all_displacements[:,i] = np.linalg.norm(mode_i, axis=1)

    return all_displacements

def get_indices(symbols, symbol):
    return np.array([enum for enum, x in enumerate(symbols) if x == f"{symbol}"])