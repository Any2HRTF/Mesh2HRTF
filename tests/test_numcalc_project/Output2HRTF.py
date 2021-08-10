# Collect the data simulated by NumCalc and save to project folder.
import numpy
import sys
sys.path.insert(1, '/home/matheson/Apps/mesh2hrtf-git/mesh2hrtf/Output2HRTF/Python')
import Output2HRTF_Main as o2hrtfm

projectPath = '/home/matheson/Documents/jthomsen/Test_material/porting_test'

Mesh2HRTF_version = '1.0.0'

sourceCenter = numpy.zeros((1, 3))
sourceArea = numpy.zeros((1, 1))

# source information
sourceType = 'vibratingElement'
sourceCenter[0, :] = [0.000000, 0.000100, 0.000002]
sourceArea[0, 0] = 6.10114e-12

# Reference to a point source in the origin
# accoring to the classical HRTF definition
# (https://doi.org/10.1016/0003-682X(92)90046-U)
reference = True

# Compute HRIRs via the inverse Fourier transfrom.
# This will add data at 0 Hz, mirror the single sided spectrum, and
# shift the HRIRs in time. Requires reference = true.
computeHRIRs = True

# Constants
speedOfSound = 346.18  # [m/s]
densityOfAir = 1.1839  # [kg/m^3]

# Distribution of ears across CPUs and cores.
# (Matrix of size [numCPUs x numCores])
cpusAndCores = numpy.array([
    [1, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0]])

# Collect the data simulated by NumCalc
o2hrtfm.Output2HRTF_Main(projectPath, Mesh2HRTF_version, cpusAndCores,
                         sourceType, sourceCenter, sourceArea, reference,
                         computeHRIRs, speedOfSound, densityOfAir)
