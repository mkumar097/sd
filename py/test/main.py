import numpy as np
from readin import *
import math

# unit conversion factors
planck = 6.626068e-34                                                          # Planck's constant in J*s
avog = 6.0221415e23                                                            # Avogadro's number
light = 2.99792458e10                                                          # Speed of light (SI) in cm/s
pi = math.pi
# derived conversion factors
#**************** ENERGY ***************#
joules2wvnbr = 1.e0/planck/light                                               # Joules (SI) to wavenumber
au2joules = 4.35974434e-18                                                     # Hartrees to Joules
wvnbr2joules = 1.e0/joules2wvnbr                                               # wavenumber to Joules
joules2au = 1.e0/au2joules                                                     # Joules to Hartrees
#**************** LENGTH ***************#
au2ang = 0.52917725e0                                                          # atomic units of length to Angstrom
ang2au = 1.0e0/au2ang                                                          # Angstrom to atomic units of length
ang2m = 1.e-10                                                                 # Angstrom to meter
#***************** TIME ****************#
amu2kg = 1.66053892e-27                                                        # atomic mass unit to kg
#************** FREQUENCY **************#
wvnbr2Hz = 3.0e10                                                              # wavenumbers to Hertz
Hz2wvnbr = 1.0e0/wvnbr2Hz                                                      # Hertz to wavenumbers
Hz2aufreq = 2.418884e-17                                                       # Hertz to atomic units of frequency

