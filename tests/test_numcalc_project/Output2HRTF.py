# Read the data simulated by NumCalc and save to the folder
# Output2HRTF inside project folder.

import numpy
import os
import mesh2hrtf as m2h

projectPath = os.getcwd()

Mesh2HRTF_version = '1.0.0'

# source information
sourceCenter = numpy.zeros((1, 3))
sourceArea = numpy.zeros((1, 1))

sourceType = 'pointSource';
sourceCenter[0, :] = [0.0, 0.20000000298023224, 0.0]
sourceArea[0, 0]     = 1
# Reference to a point source in the origin
# accoring to the classical HRTF definition
# (https://doi.org/10.1016/0003-682X(92)90046-U)
reference = False

# Compute HRIRs via the inverse Fourier transfrom.
# This will add data at 0 Hz, mirror the single sided spectrum, and
# shift the HRIRs in time. Requires reference = true.
computeHRIRs = False

# Constants
speedOfSound = 343.18  # [m/s]
densityOfAir = 1.1839  # [kg/m^3]

# Distribution of ears across CPUs and cores.
# (Matrix of size [numCPUs x numCores])
cpusAndCores = numpy.array([
    [1]])

# Collect the data simulated by NumCalc
m2h.Output2HRTF_Main(projectPath, Mesh2HRTF_version, cpusAndCores,
                 sourceType, sourceCenter, sourceArea,
                 reference, computeHRIRs,
                 speedOfSound, densityOfAir)
