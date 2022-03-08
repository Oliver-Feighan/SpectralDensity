import numpy as np

def absorption(dipoles, eigvec):
    intensities = np.empty((len(dipoles), 27))
    
    for f in range(len(dipoles)):
        dip_f = dipoles[f]
        eig_f = eigvec[f]
        
        e = random_unit_vectors(1)[0]        
        overlap = dip_f @ e
        
        overlap = np.repeat(overlap, 27).reshape(27, 27)
        
        contrib = eig_f[1:, 1:] * overlap

        intensities[f] = np.square(np.sum(contrib, axis=0))
        
    return intensities

def make_screening(distances, a=2.68, b=0.14, f0=0.54):
    screen = np.ones((distances.shape[0], 28, 28))
    
    for i in range(len(distances)):
        screen[i, 1:, 1:] = a * np.exp(-b * distances[i]) + f0
    
        np.fill_diagonal(screen[i], 1)
    
    return screen


def line_shape(intensities, energies):
    assert(intensities.shape == energies.shape)
    
    intensities = intensities.flatten()
    energies = energies.flatten()
    
    args = np.argsort(energies)
    
    energies, intensities = energies[args], intensities[args]
    
    bins = np.linspace(600, 900, 151)

    digits = np.digitize(energies, bins)
    
    summed_intensities = np.zeros(len(bins))
    
    for i in range(len(bins)):
        indices = np.where(digits == i)[0]
        
        if len(indices) == 0:
            continue
        
        summed_intensities[i] = np.sum(intensities[indices])
    
    return bins, summed_intensities/np.max(summed_intensities)


def screen_hamiltonians(distances, hamiltonians, a, b, f0):
    states, eigvec = np.linalg.eigh(make_screening(distan, a, b, f0) * hamils)

    return states, eigvec

    
def trial_screening(distances, hamiltonians, dipoles,  a, b, f0, ax, peak_max):
    states, eigvec = screen_hamiltonians(distances, hamiltonians, a, b, f0)
    
    exciton_excitations = np.array([states[:, i] - states[:, 0] for i in range(1, 28)]).T
    exciton_excitations_nm = hartree_to_nm(exciton_excitations)

    intensities = absorption(dipoles, eigvec)

    bins, summed_ints = line_shape(intensities, exciton_excitations_nm)

    shift = peak_max - bins[np.argmax(summed_ints)]
    
    ax.plot(bins+shift, summed_ints);
    

def reference(states, dipoles, eigvecs, ax, peak_max):
    exciton_excitations = np.array([states[:, i] - states[:, 0] for i in range(1, 28)]).T
    exciton_excitations_nm = hartree_to_nm(exciton_excitations)

    intensities = absorption(dipoles, eigvecs)

    bins, summed_ints = line_shape(intensities, exciton_excitations_nm)

    shift = peak_max - bins[np.argmax(summed_ints)]
    
    ax.plot(bins+shift, summed_ints);

    
def gaussian_value_generator(peak, FWHM, x):
    gamma = FWHM/2 #see https://en.wikipedia.org/wiki/Full_width_at_half_maximum
    return 1/((np.pi * gamma) * (1 + ((x - peak)/gamma)**2))


def plot_experimental(ax):
    """
    data got from "G. Absorption and CD Spectroscopy and Modelling of Various LH2 Complexes from Purple Bacteria. Biophys. J. 2002, 82, 2184−2197"
    using Rps. acidophila (10050), RT row of data from Table 1 for values
    (LH2 complex was taken from this species of bacteria)
    peaks 801 nm, 856 nm, full width at half maximum 22 and 29 respectively.
    B850 peak maximum is 1.49 times greater than B800
    """

    B850_peak = 856
    B800_peak = 801

    B850_FWHM = 29
    B800_FWHM = 22
    
    nanometer_range = np.arange(650, 950)

    B850_line = gaussian_value_generator(B850_peak, B850_FWHM, nanometer_range)
    B850_line /= max(B850_line)

    B800_line = gaussian_value_generator(B800_peak, B800_FWHM, nanometer_range)
    B800_line /= (max(B800_line) * 1.49)

    spectrum = [sum(x) for x in zip(B800_line, B850_line)]
    
    ax.plot(nanometer_range, spectrum)

    
def transition_dipole_matrix_element(eigvec, dipoles):
    res = np.zeros(3)
    
    for i, c in enumerate(eigvec[1:]):
        res += c * dipoles[i]
    
    return res


def probability(eigvec, dipoles):
    tdme  = transition_dipole_matrix_element(eigvec, dipoles)

    e = random_unit_vectors(1000)

    return np.square(np.mean(e @ tdme)), np.sum(np.square(tdme))


def probabilities(eigvec, dipoles):
    assert(len(eigvec) == len(dipoles))
    N = len(eigvec)
    
    lights = np.zeros((N, 27))
    spheres = np.zeros((N, 27))

    for f in range(N):
        for i in range(1, 28):
            light, sphere = probability(eigvec[f][:, i], dipoles[f])

            lights[f][i-1] = light
            spheres[f][i-1] = sphere
            
    return lights, spheres