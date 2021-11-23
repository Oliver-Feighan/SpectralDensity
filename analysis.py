import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy
import mdtraj
import itertools
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
    
@timer
def load_dcd_file(dcd_file, top_file):
    traj = mdtraj.load_dcd(dcd_file, top=top_file)
    top = traj.top
    
    return traj, top

@timer
def acf(x, cutoff):
    """
    hand-rolled auto-correlation function as described in https://doi.org/10.1016/j.chemphys.2018.08.013 , eq 2.
    """
    N = x.size
        
    correlations = np.zeros(N-cutoff)
    
    x_bar = np.mean(x)
    
    for j in range(0, len(correlations)):
        c_j = np.sum([(x[k] - x_bar) * (x[k+j] - x_bar) for k in range(0, N-j)])
        
        correlations[j] = c_j
                
    return correlations

def autocorr(x):
    """
    autocorrelation function that uses scipy correlate. Faster than hand-rolled version
    """
    result = scipy.signal.correlate(x, x, mode='full', method='fft')
    return result

def spectrum_and_domain(x, dt):
    auto = autocorr(x)

    spectrum = scipy.fft.fft(auto) / len(auto)

    spectrum_domain = 2 * np.pi * np.fft.fftfreq(len(auto), dt) # 2 \pi for conversion between normal to angular frequency, 20 for fs sample spacing
    
    return spectrum, spectrum_domain
    
    
def Mg_distances(traj, top):
    Mg_atom_indices = top.select("name =~ 'Mg*'")

    pairs_iterator = itertools.combinations(Mg_atom_indices, 2)

    pairs = np.array(list(pairs_iterator))

    return mdtraj.compute_distances(traj, pairs)

def Mg_pair_spectral_densities(traj, top):
    Mg_d = Mg_distances(traj, top)
    avg_d = np.mean(Mg_d, axis=0)
    
    closest_neighbours = np.where(avg_d < 1.5)[0]
    
    spectra = []
    domain = []
    
    for i in closest_neighbours:
        x = Mg_d[:,i] - np.mean(Mg_d[:,i])
        spectrum, domain = spectrum_and_domain(x, 20)
        
        spectra.append(spectrum)
        
    return np.array(spectra), domain
        
    
