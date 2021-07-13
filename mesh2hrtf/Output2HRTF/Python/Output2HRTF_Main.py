#                                Mesh2HRTF
#                Copyright (C) 2015 by Harald Ziegelwanger,
#        Acoustics Research Institute, Austrian Academy of Sciences
#                        mesh2hrtf.sourceforge.net
#
# Mesh2HRTF is licensed under the GNU Lesser General Public License as
# published by the Free Software Foundation, either version 3 of the License,
# or (at your option) any later version. Mesh2HRTF is distributed in the hope
# that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details. You should have received a
# copy of the GNU LesserGeneral Public License along with Mesh2HRTF. If not,
# see <http://www.gnu.org/licenses/lgpl.html>.
#
# If you use Mesh2HRTF:
# - Provide credits:
#   "Mesh2HRTF, H. Ziegelwanger, ARI, OEAW (mesh2hrtf.sourceforge.net)"
# - In your publication, cite both articles:
#   [1] Ziegelwanger, H., Kreuzer, W., and Majdak, P. (2015). "Mesh2HRTF:
#       Open-source software package for the numerical calculation of
#       head-related transfer functions," in Proceedings of the 22nd ICSV,
#       Florence, IT.
#   [2] Ziegelwanger, H., Majdak, P., and Kreuzer, W. (2015). "Numerical
#       calculation of listener-specific head-related transfer functions and
#       sound localization: Microphone model and mesh discretization," The
#       Journal of the Acoustical Society of America, 138, 208-222.
#
# Author: Harald Ziegelwanger
#        (Acoustics Research Institute, Austrian Academy of Sciences)
#        Fabian Brinkmann, Robert Pelzer, Jeffrey Thomsen
#        (Audio Communication Group, Technical University Berlin)

# TODO header
import os
import csv
import numpy
from netCDF4 import Dataset
import time


def Output2HRTF_Main(
        projectPath, Mesh2HRTF_version, cpusAndCores, sourceType, sourceCenter,
        sourceArea, reference, computeHRIRs, speedOfSound, densityOfAir):
    """
    Process NumCalc output and write data to disk.

    All parameters are written to Output2HRTF.py upon exporting a Mesh2HRTF
    project from Blender. This function will thus usually be called from an
    Output2HRTF.py file.

    Parameters
    ----------
    projectPath : string
        Path to the Mesh2HRTF folder, i.e., the folder containing the
        subfolders `EvaluationGrids`, `NumCal` etc.
    Mesh2HRTF_version : string
        Mesh2HRTF version as string
    cpusAndCores : array
        Array indicating which CPUs and Cores were used for the BEM
    sourceType : string
        The source type (vibrating element, point source)
    sourceCenter : array
        The source position
    sourceArea : array
        The area of the source(s) if using vibrating elements as sourceType.
        This is required for referencing the HRTFs. Use an area of 1 for a
        point source.
    reference : boolean
        Indicate if the HRTF are referenced to the pressure in the center of
        the head with the head absent.
    computeHRIRs : boolean
        Indicate if the HRIRs should be calculated by means of the inverse
        Fourier transform.
    speedOfSound : float
        The speed of sound in m/s
    densityOfAir : float
        The density of air in kg/m^3
    """

    # load meta data
    # output directory
    os.chdir(projectPath)
    if not os.path.exists(os.path.join(os.getcwd(), 'Output2HRTF')):
        os.makedirs(os.path.join(os.getcwd(), 'Output2HRTF'))

    # number of ears
    ears = numpy.amax(cpusAndCores)

    # get the evaluation grids
    evaluationGrids, evaluationGridsNumNodes = \
        read_nodes_and_elements(data='EvaluationGrids')

    # get the object mesh
    objectMeshes, objectMeshesNumNodes = \
        read_nodes_and_elements(data='ObjectMeshes')

    # Read computational effort
    print('\n Loading computational effort data ...')
    for ch in range(ears):
        computationTime = []
        for ii in range(cpusAndCores.shape[0]):
            for jj in range(cpusAndCores.shape[1]):
                if (cpusAndCores[ii, jj] == ch+1):
                    tmpFilename = os.path.join(
                        'NumCalc', 'CPU_%d_Core_%d' % ((ii+1), (jj+1)),
                        'NC.out')
                    tmp = read_computation_time(tmpFilename)
                    computationTime.append(tmp)
                    del tmp

    description = ['Frequency index', 'Frequency', 'Building', 'Solving',
                   'Postprocessing', 'Total']

    # save csv file for comparing computation times across simulations
    numpy.savetxt(os.path.join('Output2HRTF', 'computationTime.csv'),
                  computationTime[0],  fmt='%1.5e', delimiter=", ",
                  header=", ".join(description), comments="")

    del ch, ii, jj, description, computationTime

    # Load ObjectMesh data
    pressure, frequencies = read_pressure(ears, cpusAndCores,
                                          data='pBoundary')

    print('\nSave ObjectMesh data ...')
    cnt = 0
    for ii in range(len(objectMeshes)):
        nodes = objectMeshes[ii]["nodes"]
        elements = objectMeshes[ii]["elements"]
        element_data = pressure

        for jj in range(pressure.shape[1]):
            element_data[:, jj, :] = element_data[
                cnt:cnt+elements.shape[0], jj, :]

        file = open("ObjectMesh_"+objectMeshes[ii]["name"]+".npz", "w")
        numpy.savez_compressed(file.name, nodes=nodes, elements=elements,
                               frequencies=frequencies,
                               element_data=element_data)
        file.close()

        cnt = cnt + elements.shape[0]

    del pressure, nodes, elements, frequencies, ii, jj, cnt, element_data

    # Load EvaluationGrid data
    if not len(evaluationGrids) == 0:
        print('\nLoading data for the evaluation grids ...')

        pressure, frequencies = read_pressure(ears, cpusAndCores,
                                              data='pEvalGrid')

    # save to struct
    cnt = 0
    for ii in range(len(evaluationGrids)):
        evaluationGrids[ii]["pressure"] = pressure[
            cnt:cnt+evaluationGrids[ii]["num_nodes"], :, :]

        cnt = cnt + evaluationGrids[ii]["num_nodes"]

    del ii, pressure, cnt

    # reference to pressure in the middle of the head with the head absent
    # according to the HRTF definition.
    if reference:
        evaluationGrids = reference_HRTF(evaluationGrids, frequencies,
                                         sourceType, sourceArea, speedOfSound,
                                         densityOfAir, refMode=1)

    # Save complex pressure as SOFA file
    print('\nSaving complex pressure to SOFA file ...\n')

    for ii in range(len(evaluationGrids)):
        write_to_sofa(ii, evaluationGrids, Mesh2HRTF_version,
                      evaluationGridsNumNodes, frequencies,
                      ears, sourceCenter, type='HRTF')

    # Save impulse responses as SOFA file
    if computeHRIRs:

        print('\nSaving time data to SOFA file ...\n')

        for ii in range(len(evaluationGrids)):
            hrir, fs = compute_HRIR(ii, evaluationGrids, frequencies,
                                    reference, speedOfSound)

            write_to_sofa(ii, evaluationGrids, Mesh2HRTF_version,
                          evaluationGridsNumNodes, frequencies, ears,
                          sourceCenter, fs, hrir, type='HRIR')

    print('Done\n')


def read_nodes_and_elements(data):
    """
    Read the nodes and elements of the evaluation grids or object meshes.
    """
    # if data not in ['EvaluationGrids', 'ObjectMeshes']
    if not (data == 'EvaluationGrids' or data == 'ObjectMeshes'):
        raise ValueError('data must be EvaluationGrids or ObjectMeshes!')

    print('\n Loading the %s ...' % data)
    grids = []
    gridsList = os.listdir(data)
    gridsNumNodes = 0
    for ii in range(len(gridsList)):
        tmpNodes = numpy.loadtxt(os.path.join(
            data, gridsList[ii], 'Nodes.txt'),
            delimiter=' ', skiprows=1)
        tmpElements = numpy.loadtxt(os.path.join(
            data, gridsList[ii], 'Elements.txt'),
            delimiter=' ', skiprows=1)
        grids.append({"name": gridsList[ii],
                                "nodes": tmpNodes,
                                "elements": tmpElements,
                                "num_nodes": tmpNodes.shape[0]})
        gridsNumNodes += grids[ii]['num_nodes']

    return grids, gridsNumNodes


def read_pressure(ears, cpusAndCores, data):
    """Read the sound pressure on the object meshes or evaluation grid."""
    tmpPressure = []
    frequencies = []
    pressure = []

    if not (data == 'pBoundary' or data == 'pEvalGrid'):
        raise ValueError('data must be pBoundary or pEvalGrid!')

    print('\n Loading %s data ...' % data)
    for ch in range(ears):
        print('\n    Ear %d ...' % (ch+1))
        for ii in range(cpusAndCores.shape[0]):
            print('\n        CPU %d: ' % (ii+1))
            for jj in range(cpusAndCores.shape[1]):
                if (cpusAndCores[ii, jj] == (ch+1)):
                    print('%d, ' % (jj+1))
                    tmpFilename = os.path.join(
                        'NumCalc', 'CPU_%d_Core_%d' % ((ii+1), (jj+1)),
                        'be.out')
                    tmpData, tmpFrequencies = Output2HRTF_Load(tmpFilename,
                                                               data)
                    if tmpPressure:
                        tmpPressure.append(tmpData)
                        frequencies.append(tmpFrequencies)
                    else:
                        tmpPressure = tmpData
                        frequencies = tmpFrequencies
                    del tmpData, tmpFrequencies, tmpFilename
            print('...')

        idx = sorted(range(len(frequencies)), key=lambda k: frequencies[k])
        frequencies = sorted(frequencies)
        pressure.append(tmpPressure[idx, :])

        del tmpPressure

    pressure = numpy.transpose(numpy.array(pressure), (2, 0, 1))

    return pressure, frequencies


def reference_HRTF(evaluationGrids, frequencies, sourceType, sourceArea,
                   speedOfSound, densityOfAir, refMode):
    """Reference HRTF to the sound pressure in the center of the head."""
    # refmode
    # 1: reference to only one radius (the smallest found)
    # 2: reference to all individual radii
    for ii in range(len(evaluationGrids)):

        xyz = evaluationGrids[ii]["nodes"]
        pressure = evaluationGrids[ii]["pressure"]
        freqMatrix = numpy.tile(
            frequencies, (pressure.shape[0], pressure.shape[1], 1))

        # distance of source positions from the origin
        if refMode == 1:
            r = min(numpy.sqrt(xyz[:, 1]**2 + xyz[:, 2]**2 + xyz[:, 3]**2))
            r = numpy.tile(
                r,
                (pressure.shape[0], pressure.shape[1], pressure.shape[2]))
        else:
            r = numpy.sqrt(xyz[:, 1]**2 + xyz[:, 2]**2 + xyz[:, 3]**2)
            r = numpy.tile(numpy.transpose(r),
                           (pressure.shape[0], pressure.shape[1], 1))

        if sourceType == 'vibratingElement':

            volumeFlow = 0.1 * numpy.ones(
                (pressure.shape[0],
                    pressure.shape[1],
                    pressure.shape[2]))
            if 'sourceArea':
                # has to be fixed for both ears....
                for nn in range(len(sourceArea)):
                    volumeFlow[:, nn, :] = \
                        volumeFlow[:, nn, :] * sourceArea[nn]

            # point source in the origin evaluated at r
            # eq. (6.71) in: Williams, E. G. (1999). Fourier Acoustics.
            ps = -1j * densityOfAir * 2 * numpy.pi * freqMatrix * \
                volumeFlow / (4 * numpy.pi) * \
                numpy.exp(1j * 2 * numpy.pi *
                          freqMatrix / speedOfSound * r) / r

        elif sourceType == 'pointSource':

            amplitude = 0.1  # hard coded in Mesh2HRTF
            ps = amplitude * \
                numpy.exp(1j * 2 * numpy.pi * freqMatrix /
                          speedOfSound * r) / (4 * numpy.pi * r)

        else:
            raise ValueError(
                ("Referencing is currently only implemented for "
                    "sourceType 'vibratingElement' and 'pointSource'."))

        # here we go...
        evaluationGrids[ii]["pressure"] = pressure / ps

    return evaluationGrids


def compute_HRIR(ii, evaluationGrids, frequencies, reference, speedOfSound):
    """Compute HRIR from HRTF by means of the inverse Fourier transform"""
    # check if the frequency vector has the correct format
    if not all(numpy.abs(frequencies[0] - numpy.diff(frequencies)) < .1):
        raise ValueError(('The frequency vector must be if the format '
                          'a:a:fs/2, with a>0 and fs the sampling rate.'))

    if not reference:
        raise ValueError('HRIRs can only be computet if reference=true')

    pressure = evaluationGrids[ii]["pressure"]

    fs = 2*frequencies[-1]

    # add 0 Hz bin
    pressure = numpy.concatenate((numpy.ones((pressure.shape[0],
                                  pressure.shape[1], 1)), pressure),
                                 axis=2)
    # make fs/2 real
    pressure[:, :, -1] = numpy.abs(pressure[:, :, -1])
    # ifft (take complex conjugate because sign conventions differ)
    hrir = numpy.fft.irfft(numpy.conj(pressure))

    # shift 30 cm to make causal
    # (path differences between the origin and the ear are usually
    # smaller than 30 cm but numerical HRIRs show stringer pre-ringing)
    n_shift = int(numpy.round(.30 / (1/fs * speedOfSound)))
    hrir = numpy.roll(hrir, n_shift, axis=2)

    return hrir, fs


def write_to_sofa(ii, evaluationGrids, Mesh2HRTF_version,
                  evaluationGridsNumNodes, frequencies, ears,
                  sourceCenter, fs=None, hrir=None, type='HRTF'):
    """Write complex pressure or impulse responses as SOFA file."""

    xyz = evaluationGrids[ii]["nodes"]

    # Save as SOFA file
    path = os.path.join('Output2HRTF',
                        f'{type}_{evaluationGrids[ii]["name"]}%.sofa')

    # Need to delete it first if file already exists
    if os.path.exists(path):
        os.remove(path)
    Obj = Dataset(path, 'w', format='NETCDF4')

    # Required Attributes
    Obj.Conventions = 'SOFA'
    Obj.Version = '1.0'

    if type == 'HRTF':
        Obj.SOFAConventions = 'GeneralTF'
        Obj.DataType = 'TF'
    elif type == 'HRIR':
        Obj.SOFAConventions = 'GeneralFIR'
        Obj.DataType = 'FIR'

    Obj.SOFAConventionsVersion = '1.0'
    Obj.APIName = 'pysofaconventions'
    Obj.APIVersion = '0.1'
    Obj.AuthorContact = 'andres.perez@eurecat.org'
    Obj.Comment = ''

    Obj.License = 'No license provided, ask the author for permission'
    Obj.RoomType = 'free field'
    Obj.DateCreated = time.ctime(time.time())
    Obj.DateModified = time.ctime(time.time())
    Obj.Title = ''
    Obj.Organization = ''

    Obj.ApplicationName = 'Mesh2HRTF'
    Obj.ApplicationVersoin = Mesh2HRTF_version

    # Required Dimensions

    m = evaluationGridsNumNodes
    if type == 'HRTF':
        n = len(frequencies)
    elif type == 'HRIR':
        n = 2*len(frequencies)
    r = ears
    e = 1
    i = 1
    c = 3
    Obj.createDimension('M', m)
    Obj.createDimension('N', n)
    Obj.createDimension('E', e)
    Obj.createDimension('R', r)
    Obj.createDimension('I', i)
    Obj.createDimension('C', c)

    # Required Variables

    listenerPositionVar = Obj.createVariable('ListenerPosition', 'f8',
                                             ('I', 'C'))
    listenerPositionVar.Units = 'metre'
    listenerPositionVar.Type = 'cartesian'
    listenerPositionVar[:] = numpy.asarray([0, 0, 0])

    emitterPositionVar = Obj.createVariable('EmitterPosition', 'f8',
                                            ('E', 'C', 'I'))
    emitterPositionVar.Units = 'metre'
    emitterPositionVar.Type = 'cartesian'
    emitterPositionVar[:] = numpy.asarray([0, 0, 0])

    sourcePositionVar = Obj.createVariable('SourcePosition', 'f8', ('M', 'C'))
    sourcePositionVar.Units = 'metre'
    sourcePositionVar.Type = 'cartesian'
    sourcePositionVar[:] = xyz[:, 1:4]

    receiverPositionVar = Obj.createVariable('ReceiverPosition', 'f8',
                                             ('R', 'C', 'I'))
    receiverPositionVar.Units = 'metre'
    receiverPositionVar.Type = 'cartesian'
    receiverPositionVar[:] = sourceCenter

    if type == 'HRTF':
        pressure = evaluationGrids[ii]["pressure"]
        dataReal = Obj.createVariable('Data.Real', 'f8', ('M', 'R', 'N'))
        dataReal[:, :, :] = numpy.real(pressure)
        dataReal.LongName = 'pressure'
        dataReal.Units = 'pascal'
        dataImag = Obj.createVariable('Data.Imag', 'f8', ('M', 'R', 'N'))
        dataImag[:, :, :] = numpy.imag(pressure)
        dataImag.LongName = 'pressure'
        dataImag.Units = 'pascal'
        dataN = Obj.createVariable('N', 'f8', ('N'))
        dataN[:] = numpy.asarray(frequencies)
        dataN.LongName = 'frequency'
        dataN.Units = 'hertz'
    elif type == 'HRIR':
        samplingRateVar = Obj.createVariable('Data.SamplingRate', 'f8', ('I'))
        samplingRateVar.Units = 'hertz'
        samplingRateVar[:] = fs

        # delayVar = Obj.createVariable('Data.Delay', 'f8', ('I','R'))
        # delay = np.zeros((i,r))
        # delayVar[:,:] = delay

        dataIRVar = Obj.createVariable('Data.IR', 'f8', ('M', 'R', 'N'))
        dataIRVar[:, :, :] = hrir

    # Close it
    Obj.close()

    return


def Output2HRTF_Load(foldername, filename):
    """
    Load results of the BEM calculation.

    Parameters
    ----------
    foldername : string
        The folder from which the data is loaded. The data to be read is
        located in the folder be.out inside NumCalc/CPU_*_Core_*
    filename : string
        The kind of data that is loaded

        pBoundary
            The sound pressure on the object mesh
        vBoundary
            The sound velocity on the object mesh
        pEvalGrid
            The sound pressure on the evaulation grid
        vEvalGrid
            The sound velocity on the evaluation grid

    Returns
    -------
    data : numpy array
        Pressure or velocity values of shape (numFrequencies, numEntries)
    frequencies : numpy array
        The frequencies in Hz
    """

    # -----------------------check number of freq bins-------------------------
    for ii in range(1, 1_000_000):
        current_folder = os.path.join(foldername, 'be.%d' % ii, filename)
        if not os.path.exists(current_folder):
            numFrequencies = ii-1
            break

    # ---------------------check number of header and data lines---------------

    idx1 = 0
    idx2 = 0
    ii = 0
    with open(os.path.join(foldername, 'be.1', filename)) as file:
        line = csv.reader(file, delimiter=' ', skipinitialspace=True)
        for li in line:
            ii += 1
            if li not in (None, ""):
                if len(li) == 3:
                    idx1 = idx2
                    idx2 = int(li[0])
                    if idx2 - idx1 == 1:
                        numHeaderlines_BE = ii - 2
                        numDatalines = sum(1 for l in line) + 2
                        break

    # ------------------------------load data----------------------------------
    data = numpy.zeros((numFrequencies, numDatalines), dtype=complex)
    frequency = numpy.zeros(numFrequencies)

    for ii in range(numFrequencies):
        tmpData = []
        current_file = os.path.join(foldername, 'be.%d' % (ii+1), filename)
        with open(current_file) as file:
            line = csv.reader(file, delimiter=' ', skipinitialspace=True)
            for li in line:
                if line.line_num > numHeaderlines_BE:
                    tmpData.append(complex(float(li[1]), float(li[2])))

        if tmpData:
            current_file = os.path.join(
                foldername, '..', 'fe.out', 'fe.%d' % (ii+1), 'load')
            with open(current_file) as file:
                line = csv.reader(file, delimiter=' ', skipinitialspace=True)
                for li in line:
                    if line.line_num > 2:
                        tmpFrequency = float(li[0])
            data[ii, :] = tmpData
            frequency[ii] = tmpFrequency
        else:
            data[ii, :] = numpy.nan
            frequency[ii] = numpy.nan

    return data, frequency


def read_computation_time(filename):
    """
    Read compuation time

    Parameters
    ----------
    filename : string

    Returns
    -------
    computation_time : numpy array
    """

    f = open(filename, "r", encoding="utf8", newline="\n")
    count = -1

    for line in f:
        if line.find('Number of frequency steps') != -1:
            idx = line.find('=')
            nSteps = int(line[idx+2])
            data = numpy.zeros((nSteps, 6), dtype=int)

        if line.find('Frequency') != -1:
            count += 1
            idx1 = line.find('Step')
            idx2 = line.find('Frequency')
            data[count, 0] = int(line[idx1+5:idx2-2])
            idx1 = line.find('=')
            idx2 = line.find('Hz')
            data[count, 1] = int(line[idx1+2:idx2-1])

        if line.find('Assembling the equation system  ') != -1:
            idx = line.find(':')
            data[count, 2] = int(line[idx+2:-1])

        if line.find('Solving the equation system  ') != -1:
            idx = line.find(':')
            data[count, 3] = int(line[idx+2:-1])

        if line.find('Post processing  ') != -1:
            idx = line.find(':')
            data[count, 4] = int(line[idx+2:-1])

        if line.find('Total  ') != -1:
            idx = line.find(':')
            data[count, 5] = int(line[idx+2:-1])

        if line.find('Address computation ') != -1:
            return data
