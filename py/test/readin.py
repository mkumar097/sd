import numpy as np
from periodic import element

# read in the coordinates for your system
num_atoms = len(open('coords.dat', 'r').readlines())
coords = np.zeros((num_atoms, 3))
atm_mass = np.zeros(num_atoms)

fn = open('coords.dat', 'r')
for i in range(num_atoms):
    c_info = fn.readline().split()
    atm = c_info[0]
    coords[i, :] = map(float, c_info[1:4])
    atm_mass[i] = element(atm).mass
fn.close()

# read in the excited state forces on your atoms
num_eforces = len(open('eforces.dat', 'r').readlines())
eforces = np.zeros((num_eforces, 4))

fn = open('eforces.dat', 'r')
for i in range(num_eforces-1):
    e_info = fn.readline().split()
    eforces[i, :] = map(float, e_info[1:5])
fn.close()

# read in the ground state forces on your atoms
num_gforces = len(open('gforces.dat', 'r').readlines())
gforces = np.zeros((num_gforces, 4))

fn = open('gforces.dat', 'r')
for i in range(num_gforces-1):
    g_info = fn.readline().split()
    gforces[i, :] = map(float, g_info[1:5])
fn.close()
