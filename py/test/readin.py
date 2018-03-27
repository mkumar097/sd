import numpy as np
from periodic import element


def readin(coords_dat='coords.dat', eforces_dat='eforces.dat',
           gforces_dat='gforces.dat'):
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

    # read in the modes for your system
    with open(modes.dat, 'r') as file_modes:
        lines_modes = file_modes.readlines()

    modes = np.zeros((num_atoms*3,num_atoms))

    for i, line in enumerate((num_atoms*3)-6):
        m_info = line.split()
        if m_info[0] == i:
            for j in enumerate(num_atoms):
                m_info = line.split()
                modes[i,j] = map(float, m_info[3])
        else:
            print "Can't properly read in the normal modes for your system."
            break

    # read in the excited state forces on your atoms
    eforces = read_forces(eforces_dat)

    # read in the ground state forces on your atoms
    gforces = read_forces(gforces_dat)

    return coords, eforces, gforces


def read_forces(forces_dat):
    with open(forces_dat, 'r') as file_forces:
        lines_forces = file_forces.readlines()
    num_forces = len(lines_forces)
    forces = np.zeros((num_forces, 4))

    for i, line in enumerate(lines_forces[:-2]):  # skipping last line
        e_info = line.split()
        forces[i, :] = map(float, e_info[1:5])
    return forces

