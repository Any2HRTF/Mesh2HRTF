# Read the data simulated by NumCalc and save to the folder
# Output2HRTF inside project folder.

import numpy
import os
import mesh2hrtf as m2h

Mesh2HRTF_version = '1.0.0'

# source information
sourceCenter = numpy.zeros((2, 3))
sourceArea = numpy.zeros((2, 1))

sourceType = 'Both ears'
numSources = 2
sourceCenter[0, :] = [0.074648, -0.056753, 0.033037]
sourceArea[0, 0] = 7.13566e-06
sourceCenter[1, :] = [-0.008213, 0.098506, -0.014679]
sourceArea[1, 0] = 5.83617e-06

# Reference to a point source in the origin
# accoring to the classical HRTF definition
# (https://doi.org/10.1016/0003-682X(92)90046-U)
reference = False

# Compute HRIRs via the inverse Fourier transfrom.
# This will add data at 0 Hz, mirror the single sided spectrum, and
# shift the HRIRs in time. Requires reference = true.
computeHRIRs = False

# Constants
speedOfSound = 343  # [m/s]
densityOfAir = 1.1839  # [kg/m^3]

# Collect the data simulated by NumCalc
m2h.Output2HRTF_Main(Mesh2HRTF_version, sourceType, 
                     numSources, sourceCenter, sourceArea,
                     reference, computeHRIRs,
                     speedOfSound, densityOfAir)
