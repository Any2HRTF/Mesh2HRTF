# TODO HEADER
# adapted from Andrés Pérez-López
# plotListenHRTF.py
# could still be extended to read HRIRs with help of his sofainfo.py script

import pysofaconventions as pysofa
import numpy


def hrtf_sofa_to_numpy(path):

    # path = '/home/matheson/Documents/jthomsen/Test_material/porting_test/Output2HRTF/HRTF_ARI.sofa'
    sofa = pysofa.SOFAFile(path, 'r')

    # File is actually not valid, but we can forgive them
    print("\n")
    print("File is valid:", sofa.isValid())

    # Convention is SimpleFreeFieldHRIR
    print("\n")
    print("SOFA Convention:", sofa.getGlobalAttributeValue('SOFAConventions'))

    # Let's see the dimensions:
    print("\n")
    print("Dimensions:")
    sofa.printSOFADimensions()

    # Read the data
    dataReal = sofa.getVariableValue('Data.Real')
    dataReal = numpy.asarray(dataReal)
    dataImag = sofa.getVariableValue('Data.Imag')
    dataImag = numpy.asarray(dataImag)
    # and get the HRTF associated with m=0
    hrtf = numpy.array(dataReal+1j*dataImag, dtype=complex)

    return hrtf
