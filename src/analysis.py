import argparse
import numpy as np
import pandas as pd
import pathlib
import scipy
import scipy.signal
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

def first_half(arr):
    return arr[1:len(arr)//2]

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
    
    n=5
    tiled = np.tile(auto, n)

    spectrum = scipy.fft.fft(tiled)

    #spectrum_frequency_domain = 2 * np.pi * np.fft.fftfreq(len(auto), dt) # 2 \pi for conversion between normal to angular frequency
    spectrum_normal_domain = np.fft.fftfreq(len(tiled), dt) # 2 \pi for conversion between normal to angular frequency
    
    return tiled, spectrum, spectrum_normal_domain
    
    
def Mg_distances(traj, top):
    Mg_atom_indices = top.select("name =~ 'Mg*'")

    pairs_iterator = itertools.combinations(Mg_atom_indices, 2)

    pairs = np.array(list(pairs_iterator))

    return 10 * mdtraj.compute_distances(traj, pairs) # factor of 10 as mdtraj is in nanometers

def Mg_pair_spectral_densities(traj, top, dt):
    Mg_d = Mg_distances(traj, top)
    
    autocorrs = []
    spectra = []
    domain = []
    
    for i in range(Mg_d.shape[1]):
        x = Mg_d[:,i] - np.mean(Mg_d[:,i])
        auto, spectrum, domain = spectrum_and_domain(x, dt)
        
        autocorrs.append(auto)
        spectra.append(spectrum)
        
    return Mg_d, np.array(autocorrs), np.array(spectra), domain

@timer
def load_dcd_file(dcd_file, top_file):
    return mdtraj.load_dcd(dcd_file, top=top_file)

def make_parser():
    parser = argparse.ArgumentParser(description="""Calculate and save the spectral density and frequency domain for
 a given MD trajectory. Currently only does spectral density of inter-chromophore Mg distances.""")
    
    parser.add_argument("--inpcrd", action="store", type=pathlib.Path, help="Path to the inpcrd file.")
    parser.add_argument("--prmtop", action="store", type=pathlib.Path, help="Path to the prmtop file.")
    parser.add_argument("--dcd", action="store", type=pathlib.Path, help="Path to the dcd file.")
    parser.add_argument("--sample_period", action="store", type=int, help="Number of femtoseconds between frame sampling - necessary for fourier transform domain")
    
    return parser
    
    
if __name__ == "__main__":
    parser = make_parser()
    args = parser.parse_args()
    
    prmtop_path = args.prmtop.resolve()
    inpcrd_path = args.inpcrd.resolve()
    dcd_path = args.dcd.resolve()
    
    
    print("loading trajectory...")
    traj = load_dcd_file(str(dcd_path), str(prmtop_path))
    top = traj.top

    print(f"n. frames : {traj.n_frames}")   
    
    print("calculating spectral density...")
    dt = args.sample_period
    Mg_d, autocorrs, spectra, domain = Mg_pair_spectral_densities(traj, traj.top, dt)
    
    
    print("saving...")
    basename = dcd_path.name.replace(".dcd", "")

    np.save(f"{basename}_distances", Mg_d)
    np.save(f"{basename}_autocorrs", autocorrs)
    np.save(f"{basename}_spectra", spectra)
    np.save(f"{basename}_domain", domain)
    
        
