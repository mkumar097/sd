from __future__ import absolute_import
from readin import readin
from constants_py import *
import numpy as np


coords, eforces, gforces, atm_mass, nm_modes = readin()                        # read in the coordinates, modes, forces


eforces *=  au2joules / (au2ang * ang2m) / np.sqrt(atm_mass)                   # excited state forces are in m/s**2
gforces *=  au2joules / (au2ang * ang2m) / np.sqrt(atm_mass)                   # excited state forces are in m/s**2



