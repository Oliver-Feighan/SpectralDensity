import numpy as np

import functools
import time


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

def first_half(arr):
    return arr[1:len(arr)//2]


def random_unit_vectors_polar(n):
    phi = np.random.uniform(0, 2*np.pi, n)
    psi = np.random.uniform(0, np.pi, n)
    
    x = np.cos(phi) * np.sin(psi)
    y = np.sin(phi) * np.sin(psi)
    z = np.cos(psi)
    
    return np.vstack((x, y, z)).T


def random_unit_vectors(n):
    rand = np.random.normal(0, 0.01, size=(n, 3))
    
    norms = np.linalg.norm(rand, axis=1)
    norm_stk = np.vstack((norms, norms, norms)).T
    
    return rand/norm_stk


def random_unit_vectors_uniform(n):
    rand = np.random.uniform(-1, 1, size=(n, 3))
    
    norms = np.linalg.norm(rand, axis=1)
    norm_stk = np.vstack((norms, norms, norms)).T
    
    return rand/norm_stk


def hartree_to_nm(x):
    c = 299792458 #speed of light, m/s
    h = 6.626e-34 #Planck constant
    J_per_Eh = 4.35974e-18 # J/Eh
    
    wavelength = (h*c) / (x * J_per_Eh)
    
    return 1e9 * wavelength # m to nm conversion


def hartree_to_pcm(x):
    c = 299792458 #speed of light, m/s
    h = 6.626e-34 #Planck constant
    J_per_Eh = 4.35974e-18 # J/Eh
    
    wavenumber = (x * J_per_Eh) / (h*c)
    
    return 1e-2 * wavenumber # m-1 to cm-1 conversion


def check_symmetric(a, rtol=1e-05, atol=1e-08):
    return np.allclose(a, a.T, rtol=rtol, atol=atol)


def reconstruct_H(vec, val):
    val_mat = np.empty((len(val), len(val)))
    np.fill_diagonal(val_mat, val)
                       
    return vec @ val_mat @ np.linalg.inv(vec)


def reconstruct_val(vec, H):
    return (np.linalg.inv(vec) @ H @ vec).diagonal()