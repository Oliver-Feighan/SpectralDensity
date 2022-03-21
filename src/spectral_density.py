import argparse
import numpy as np
import scipy
import scipy.signal
import mdtraj
import functools
import time
import utils

@utils.timer
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
    if abs(x-mu) > 3 * s:
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
    
    funcs = [make_gaussian(d0, d0 * broadening, amp) for amp, d0, diff in zip(series, domain, np.diff(domain)) if amp > cutoff]
    return lambda x : np.sum([g(x) for g in funcs])


def plot_gaussians(broadening, max_amp, ax, p_range, color='black'):
    #plotting
    start, end, number = p_range
    fs = np.logspace(np.log10(start), np.log10(end), number)
    #fs = np.linspace(start, end, number)
    b = np.array([broadening(x) for x in fs])
    
    b *= max_amp / max(b)
    
    ax.plot(fs, b, color=color)
    
    return b

    
def broadened_spectral_density(prop, frame_rate, ax, broadening, p_range, color='black'):
    #load energies
    prop_rel = prop - np.mean(prop)

    #spectral density
    autocorr, spectrum, spectrum_normal_domain = spectrum_and_domain(prop_rel, frame_rate)
    
    max_amp = max(utils.first_half(np.abs(spectrum)))
    
    #gaussian broadening
    cutoff = 0.05 * max_amp
    broadening = gaussians_variable_broadening(utils.first_half(np.abs(spectrum)), utils.first_half(spectrum_normal_domain), broadening, cutoff=cutoff)
    
    
    broadened = plot_gaussians(broadening, max_amp, ax, color=color, p_range=p_range)
    
    return autocorr, spectrum, spectrum_normal_domain, broadened

def average_spectral_density(prop, dt, indices):
    all_autocorr = [[] for i in indices]
    all_spectra = [[] for i in indices]
    all_domain = [[] for i in indices]

    for enum, i in enumerate(indices):
        prop_rel = prop[:, i] - np.mean(prop[:, i])

        #spectral density
        autocorr, spectrum, spectrum_normal_domain = spectrum_and_domain(prop_rel, dt)
   
        all_autocorr[enum] = utils.first_half(autocorr)
        all_spectra[enum] = utils.first_half(np.abs(spectrum))
        all_domain[enum] = utils.first_half(spectrum_normal_domain)
        
    all_autocorr = np.array(all_autocorr)
    all_spectra = np.array(all_spectra)
    all_domain = np.array(all_domain)
        
    return all_autocorr, all_spectra, all_domain, np.average(all_autocorr, axis=0), np.average(all_spectra, axis=0), np.average(all_domain, axis=0)


def broadened_average_spectral_density(prop, dt, ax, indices, broadening, p_range, color='black'):
    #load energies
    spectral_objects = average_spectral_density(prop, dt, indices)
    
    all_autocorr, all_spectra, all_domain = spectral_objects[:3]
    avg_autocorr, avg_spectrum, avg_domain = spectral_objects[3:]
    
    max_amp = max(utils.first_half(avg_spectrum))
    
    #gaussian broadening
    cutoff = 0.05 * max_amp
    broadening = gaussians_variable_broadening(avg_spectrum, avg_domain, broadening=broadening, cutoff=cutoff)
    
    broadened = plot_gaussians(broadening, max_amp, ax, p_range, color=color)
    
    return all_autocorr, all_spectra, all_domain, broadened
