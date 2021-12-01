import numpy as np
import re
import matplotlib.pyplot as plt

from hessian import *
from system import *

def phytol_indices():
    return np.array([64, 134, 135, 136, 63, 65, 139, 137, 133, 138, 62,
                     132, 131, 61, 129, 60, 128, 127, 58, 123, 124, 59,
                     126, 125, 121, 122, 57, 56, 120, 119, 118, 117, 55,
                     53, 113, 114, 54, 116, 115, 52, 112, 111, 51, 109,
                     110, 107, 50, 108, 48, 104, 49, 106, 105, 47, 103,
                     46, 101, 102, 15, 13, 14, 130, 13])

def read_xyz(file_name):   
    lines = list(open(file_name))
    
    symbols = [] 
    coords = []
    
    lines = lines[2:]

    for line in lines:
        symbol = re.findall(r'[a-zA-Z]+', line)
        coord  = np.array([float(y) for y in re.findall(r'-?\d+.\d+', line)])
            
        if len(symbol) == 0 or len(coord) == 0:
            continue
                
        symbols.append(symbol[0])
        coords.append(coord)
    
    return symbols, np.array(coords)
            
def write_xyz(file_name, geom, symbols):
    assert(len(geom) == len(symbols))
    
    n_atoms = len(symbols)
    
    with open(file_name, 'w') as f:
        print(f"{n_atoms} \n", file=f)
        
        for i in range(n_atoms):
            print(f"{symbols[i]}\t{geom[i][0]}\t{geom[i][1]}\t{geom[i][2]}", file=f)
        
    
    
def rotation_matrix(axis, angle):
    cos = np.cos(angle)
    sin = np.sin(angle)
    
    cross_product_matrix = np.array([
                                    [ 0,         -axis[2],  axis[1]],
                                    [ axis[2],  0,         -axis[0]],
                                    [-axis[1],  axis[0],  0        ]
                                   ])
    
    outer_product = np.outer(axis, axis)
    
    identity = np.identity(3)
    
    R = cos * identity + sin * cross_product_matrix + (1-cos) * outer_product
    
    return R

def get_Nabcd(symbols):
    N_indices = get_indices(symbols, "N")

    Na, Nb, Nc, Nd = N_indices[1], N_indices[2], N_indices[3], N_indices[0]

    return Na, Nb, Nc, Nd

def normal_to_porphoryn(geom, symbols):
    Na, Nb, Nc, Nd = get_Nabcd(symbols)
    
    Qy = geom[Na]-geom[Nc]
    Qx = geom[Nb]-geom[Nd]
    
    normal = np.cross(Qx, Qy)
    
    return normal/np.linalg.norm(normal)

def translate(geom, vector):
    return geom + vector
    
def translate_Mg_to_origin(geom, symbols):
    Mg_index = get_indices(symbols, "Mg")[0]
    
    return translate(geom, 0-geom[Mg_index])

def rotate(geom, angle, axis):
    R = rotation_matrix(angle, axis)
    
    return np.matmul(geom, R)

def angle_between(v1, v2):
    v1_u = v1 / np.linalg.norm(v1)
    v2_u = v2 / np.linalg.norm(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

def rotate_QyQx_to_xy(geom, symbols):
    normal = normal_to_porphoryn(geom, symbols)
    
    unit_normal = normal/np.linalg.norm(normal)
    
    assert(np.linalg.norm(unit_normal) == 1.)
    
    z_axis = np.array([0.,0.,1.])
    
    axis = np.cross(z_axis, unit_normal)
    axis /= np.linalg.norm(axis)
    
    angle = angle_between(normal, z_axis)
    
    return rotate(geom, axis, angle)

def translate_and_orient(geom, symbols):
    translated = translate_Mg_to_origin(geom, symbols)
    
    oriented = rotate_QyQx_to_xy(translated, symbols)
    
    assert(np.allclose(distance_matrix(geom), distance_matrix(translated)))
    assert(np.allclose(distance_matrix(geom), distance_matrix(oriented)))
    assert(np.allclose(distance_matrix(translated), distance_matrix(oriented)))
    
    return oriented

def distance_matrix(geom):
    sqr = np.square(geom)
    d = -2 * np.matmul(geom, geom.T)
        
    addition = lambda x : x + np.sum(sqr, 1)

    d = np.apply_along_axis(addition, 0, arr=d)
    d = np.apply_along_axis(addition, 1, arr=d)

    np.fill_diagonal(d,0)
    
    return np.sqrt(d)

def plot_molecule(geom, symbols, ax, with_phytol=False):
    symbol_colors = {"Mg" : "gold", "H" : "white", "O" : "red", "N" : "blue", "C" : "grey"}

    atom_types = set(symbols)

    for symbol in atom_types:
        indices = get_indices(symbols, symbol)
        
        if not with_phytol:
            indices = [x for x in indices if x not in phytol_indices()]
    
        ax.scatter(geom[:,0][indices], geom[:,1][indices], color=symbol_colors[symbol])

    distances = distance_matrix(geom)
    
    for i in range(len(geom)):
        for j in range(0, i):
            if (i in phytol_indices() or j in phytol_indices()) and not with_phytol:
                continue
            
            if distances[i][j] < 1.7:
                x = [geom[:, 0][i], geom[:, 0][j]]
                y = [geom[:, 1][i], geom[:, 1][j]]

                ax.plot(x, y, color='black', zorder=0)

    ax.set_aspect('equal')
    ax.set_facecolor('lightgrey')
    
def get_bchla_in_xy_plane():
    symbols, geom = read_xyz("opt_bchla.xyz")
    bchla_in_xy_plane = translate_and_orient(geom, symbols)
    
    return bchla_in_xy_plane
    