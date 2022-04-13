# Read the data simulated by NumCalc and save to the folder
# Output2HRTF inside project folder.

import numpy
import mesh2hrtf as m2h

Mesh2HRTF_version = '1.0.0'

# source information
sourceCenter = numpy.zeros((1, 3))
sourceArea = numpy.zeros((1, 1))

sourceType = 'Right ear'
numSources = 1
sourceCenter[0, :] = [0.087107, -0.003103, 0.003888]
sourceArea[0, 0] = 3.68754e-05

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
m2h.output_two_hrtf(Mesh2HRTF_version, sourceType,
                    numSources, sourceCenter, sourceArea,
                    reference, computeHRIRs,
                    speedOfSound, densityOfAir)
