# Read the data simulated by NumCalc and save to the folder
# Output2HRTF inside project folder.

import numpy
import os
import mesh2hrtf as m2h

projectPath = os.getcwd()

Mesh2HRTF_version = '1.0.0'

# source information
sourceCenter = numpy.zeros((2, 3))
sourceArea = numpy.zeros((2, 1))

sourceType = 'Both ears'
numSources = 2
sourceCenter[0, :] = [-0.002530, 0.087218, 0.000632]
sourceArea[0, 0] = 3.26585e-05
sourceCenter[1, :] = [0.004010, -0.087165, -0.000470]
sourceArea[1, 0] = 3.68826e-05

# Reference to a point source in the origin
# accoring to the classical HRTF definition
# (https://doi.org/10.1016/0003-682X(92)90046-U)
reference = True

# Compute HRIRs via the inverse Fourier transfrom.
# This will add data at 0 Hz, mirror the single sided spectrum, and
# shift the HRIRs in time. Requires reference = true.
computeHRIRs = True

# Constants
speedOfSound = 343  # [m/s]
densityOfAir = 1.1839  # [kg/m^3]

# Collect the data simulated by NumCalc
m2h.Output2HRTF_Main(Mesh2HRTF_version, sourceType,
                     numSources, sourceCenter, sourceArea,
                     reference, computeHRIRs,
                     speedOfSound, densityOfAir)
