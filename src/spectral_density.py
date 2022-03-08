import argparse
import numpy as np
import scipy
import scipy.signal
import mdtraj
import functools
import time

def timer(func):
    """
    decorator for timing functions.
    """
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
    """
    Calculates spectral density from a given series. Does *NOT* take the series as deviations from the mean, which can sometimes
    be necessary for correct behaviour of the autocorrelation function.
    
        Parameters:
            x - the series
            dt - the 'frame rate' of the series. For example, a set of distances from an MD traj. taken every 2fs would have a frame
                rate of 2fs.
                
        Returns:
            tiled - repitition of the autocorrelation function (using np.tile()). The number of repeats is 5
            spectrum - the spectral density, windowed by the scipy.signal.flattop() function
            spectrum_normal_domain - the frequency domain of the spectrum. This can be reported as angular or normal frequency, and
                this function returns normal frequency. Multiply by 2\pi for angular frequency.
    
    """
    auto = autocorr(x)
    
    n=5
    tiled = np.tile(auto, n)

    window = scipy.signal.flattop(len(tiled), sym=False)
    windowed_tiled = window * tiled
    
    spectrum = scipy.fft.fft(windowed_tiled)

    #spectrum_frequency_domain = 2 * np.pi * np.fft.fftfreq(len(auto), dt) # 2 \pi for conversion between normal to angular frequency
    spectrum_normal_domain = np.fft.fftfreq(len(tiled), dt) # 2 \pi for conversion between normal to angular frequency
    
    return tiled, spectrum, spectrum_normal_domain


def gaussian(x, mu, s, a):
    """
    Gaussian function
        
        Parameters:
            x - position
            mu - mean
            s - standard deviation
            a - amplitude
            
        Returns:
            if x is within 3 standard deviations:
                a * np.exp(-0.5 * ((x-mu) / s)**2)
            else:
                0
            
    """
    if abs(x-s) > 3 * s:
        return 0
    
    return a * np.exp(-0.5 * ((x-mu) / s)**2)


def make_gaussian(mu, s, a):
    """
    returns wrapper for a gaussian function - i.e. saves mu, s and a. See gaussian() for more details.
    """
    return lambda x : gaussian(x, mu, s, a)


def gaussians_static_broadening(series, domain, broadening, cutoff):
    
    funcs = [make_gaussian(d0, broadening, amp) for amp, d0 in zip(series, domain) if amp > cutoff]
    return lambda x : np.sum([g(x) for g in funcs])


def gaussians_variable_broadening(series, domain, broadening, cutoff):
    
    funcs = [make_gaussian(d0, diff * broadening, amp) for amp, d0, diff in zip(series, domain, np.diff(domain)) if amp > cutoff]
    return lambda x : np.sum([g(x) for g in funcs])


def plot_gaussians(broadening, max_amp, ax, color='black'):
    #plotting
    fs = np.linspace(0, 5e-2, 1000)
    b = np.array([broadening(x) for x in fs])
    
    b *= max_amp / max(b)
    
    ax.plot(fs, b, color=color)

    
def broadened_spectral_density(prop, frame_rate, ax, color='black'):
    #load energies
    prop_rel = prop - np.mean(prop)

    #spectral density
    autocorr, spectrum, spectrum_normal_domain = analysis.spectrum_and_domain(prop_rel, frame_rate)
    
    max_amp = max(analysis.first_half(np.abs(spectrum)))
    
    #gaussian broadening
    cutoff = 0.1 * max_amp
    broadening = gaussians_variable_broadening(analysis.first_half(np.abs(spectrum)), analysis.first_half(spectrum_normal_domain), broadening=50, cutoff=cutoff)
    
    
    plot_gaussians(broadening, max_amp, ax)
    
    return autocorr, spectrum, spectrum_normal_domain