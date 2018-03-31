import numpy as np
import math
from periodic import element

amu2kg = 1.66053892e-27                                                        # atomic mass unit to kg

def readin(coords_dat='coords.dat', eforces_dat='eforces.dat',
           gforces_dat='gforces.dat', modes_dat='modes.dat'):

    # read in the coordinates for your system
    with open(coords_dat, 'r') as file_coords:
        lines_coords = file_coords.readlines()
    num_atoms = len(lines_coords)
    coords = np.zeros((num_atoms, 3))
    atm_mass = np.zeros(num_atoms)

    for i, line in enumerate(lines_coords):
        c_info = line.split()
        coords[i, :] = map(float, c_info[1:4])
        atm_mass[i] = element(c_info[0]).mass

    atm_mass *= amu2kg

    # read in the modes for your system
    num_modes = (3*num_atoms)
    red_mass = np.zeros(num_modes)
    modes = np.zeros((num_modes, num_modes))
    with open(modes_dat, 'r') as file_modes:
        for i in range(num_modes-6):
            nm_info = file_modes.readline().split()
            red_mass[i+6] = float(nm_info[5])
            if nm_info[1] == str(i+1):
                norm = 0.0
                for j in range(3*num_atoms):
                    mode_info = file_modes.readline().split()
                    modes[i+6, j] = float(mode_info[3])*math.sqrt(atm_mass[int(math.floor(j/3))])
                    norm += modes[i+6, j]*modes[i+6, j]
                modes[i+6, :] /= math.sqrt(norm)
            else:
                print ("Unable to correctly read in the normal modes for your system.")

    # add a reduced mass of 1.0 for the first 6 normal modes which are rotational and translational modes
    red_mass[0:6] = 1.0

    # read in the excited state forces on your atoms
    eforces = read_forces(eforces_dat)

    # read in the ground state forces on your atoms
    gforces = read_forces(gforces_dat)

    return coords, eforces, gforces, atm_mass, modes


def read_forces(forces_dat):
    with open(forces_dat, 'r') as file_forces:
        lines_forces = file_forces.readlines()
    num_forces = len(lines_forces)
    forces = np.zeros((num_forces, 4))

    for i, line in enumerate(lines_forces[:-2]):  # skipping last line
        e_info = line.split()
        forces[i, :] = map(float, e_info[1:5])
    return forces

