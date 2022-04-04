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
import numpy as np
import sofar as sf


def Output2HRTF_Main(
        Mesh2HRTF_version, sourceType, numSources, sourceCenter,
        sourceArea, reference, computeHRIRs, speedOfSound, densityOfAir):
    """
    Process NumCalc output and write data to disk.

    All parameters are written to Output2HRTF.py upon exporting a Mesh2HRTF
    project from Blender. This function will thus usually be called from an
    Output2HRTF.py file.

    Parameters
    ----------
    Mesh2HRTF_version : string
        Mesh2HRTF version as string
    numSources : integer
        Number of sound sources in the simulation
    sourceType : string
        The source type ('Both ears', 'Left ear', 'Right ear',
                         'Point source', 'Plane wave')
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
    if not os.path.exists(os.path.join(os.getcwd(), 'Output2HRTF')):
        os.makedirs(os.path.join(os.getcwd(), 'Output2HRTF'))

    # get the number of frequency steps
    numFrequencies = get_number_of_frequencies('Info.txt')

    # get the evaluation grids
    evaluationGrids, evaluationGridsNumNodes = \
        read_nodes_and_elements(data='EvaluationGrids')

    # get the object mesh
    objectMeshes, objectMeshesNumNodes = \
        read_nodes_and_elements(data='ObjectMeshes')

    # # Read computational effort
    # print('\n Loading computational effort data ...')
    # computationTime = []
    # for source in range(numSources):
    #     for file in Path(os.path.join('NumCalc',
    #                                   f'source_{source+1}')).glob('NC*.out'):
    #         tmpFilename = os.path.join(file)
    #         tmp = read_computation_time(tmpFilename)
    #         computationTime.append(tmp)
    #         del tmp

    # description = ['Frequency index', 'Frequency', 'Building', 'Solving',
    #                'Postprocessing', 'Total']

    # # save csv file for comparing computation times across simulations
    # np.savetxt(os.path.join('Output2HRTF', 'computationTime.csv'),
    #               computationTime[0],  fmt='%1.5e', delimiter=", ",
    #               header=", ".join(description), comments="")

    # del source, file, description, computationTime

    # Load ObjectMesh data
    pressure, frequencies = read_pressure(
        numSources, numFrequencies, data='pBoundary')

    print('\nSaving ObjectMesh data ...')
    cnt = 0
    for ii in range(len(objectMeshes)):
        nodes = objectMeshes[ii]["nodes"]
        elements = objectMeshes[ii]["elements"]
        element_data = pressure

        for jj in range(pressure.shape[1]):
            element_data[:, jj, :] = element_data[
                cnt:cnt+elements.shape[0], jj, :]

        file = open("ObjectMesh_"+objectMeshes[ii]["name"]+".npz", "w")
        np.savez_compressed(file.name, nodes=nodes, elements=elements,
                            frequencies=frequencies,
                            element_data=element_data)
        file.close()

        cnt = cnt + elements.shape[0]

    del pressure, nodes, elements, frequencies, ii, jj, cnt, element_data

    # Load EvaluationGrid data
    if not len(evaluationGrids) == 0:
        print('\nLoading data for the evaluation grids ...')

        pressure, frequencies = read_pressure(
            numSources, numFrequencies, data='pEvalGrid')

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
                      frequencies, numSources, sourceCenter, type='HRTF')

    # Save impulse responses as SOFA file
    if computeHRIRs:

        print('\nSaving time data to SOFA file ...\n')

        for ii in range(len(evaluationGrids)):
            hrir, fs = compute_HRIR(ii, evaluationGrids, frequencies,
                                    reference, speedOfSound)

            write_to_sofa(ii, evaluationGrids, Mesh2HRTF_version, frequencies,
                          numSources, sourceCenter, fs, hrir, type='HRIR')

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
        tmpNodes = np.loadtxt(os.path.join(
            data, gridsList[ii], 'Nodes.txt'),
            delimiter=' ', skiprows=1)
        tmpElements = np.loadtxt(os.path.join(
            data, gridsList[ii], 'Elements.txt'),
            delimiter=' ', skiprows=1)
        grids.append({"name": gridsList[ii], "nodes": tmpNodes,
                      "elements": tmpElements, "num_nodes": tmpNodes.shape[0]})
        gridsNumNodes += grids[ii]['num_nodes']

    return grids, gridsNumNodes


def read_pressure(numSources, numFrequencies, data):
    """Read the sound pressure on the object meshes or evaluation grid."""
    pressure = []

    if not (data == 'pBoundary' or data == 'pEvalGrid'):
        raise ValueError('data must be pBoundary or pEvalGrid!')

    print('\n Loading %s data ...' % data)
    for source in range(numSources):
        print('\n    Source %d ...' % (source+1))

        tmpFilename = os.path.join('NumCalc', f'source_{source+1}', 'be.out')
        tmpPressure, frequencies = Output2HRTF_Load(
            tmpFilename, data, numFrequencies)

        pressure.append(tmpPressure)

    pressure = np.transpose(np.array(pressure), (2, 0, 1))

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
        freqMatrix = np.tile(
            frequencies, (pressure.shape[0], pressure.shape[1], 1))

        # distance of source positions from the origin
        if refMode == 1:
            r = min(np.sqrt(xyz[:, 1]**2 + xyz[:, 2]**2 + xyz[:, 3]**2))
            r = np.tile(
                r,
                (pressure.shape[0], pressure.shape[1], pressure.shape[2]))
        else:
            r = np.sqrt(xyz[:, 1]**2 + xyz[:, 2]**2 + xyz[:, 3]**2)
            r = np.tile(np.transpose(r),
                        (pressure.shape[0], pressure.shape[1], 1))

        if sourceType in {'Both ears', 'Left ear', 'Right ear'}:

            volumeFlow = 0.1 * np.ones(
                (pressure.shape[0], pressure.shape[1], pressure.shape[2]))
            if 'sourceArea':
                # has to be fixed for both ears....
                for nn in range(len(sourceArea)):
                    volumeFlow[:, nn, :] = \
                        volumeFlow[:, nn, :] * sourceArea[nn]

            # point source in the origin evaluated at r
            # eq. (6.71) in: Williams, E. G. (1999). Fourier Acoustics.
            ps = -1j * densityOfAir * 2 * np.pi * freqMatrix * \
                volumeFlow / (4 * np.pi) * \
                np.exp(1j * 2 * np.pi * freqMatrix / speedOfSound * r) / r

        elif sourceType == 'Point source':

            amplitude = 0.1  # hard coded in Mesh2HRTF
            ps = amplitude * \
                np.exp(1j * 2 * np.pi * freqMatrix /
                       speedOfSound * r) / (4 * np.pi * r)

        elif sourceType == 'Plane wave':
            raise ValueError(
                ("Plane wave not implemented yet."))

        else:
            raise ValueError(
                ("Referencing is currently only implemented for "
                    "sourceType 'Both ears', 'Left ear', 'Right ear',"
                    "'Point source' and 'Plane wave'."))

        # here we go...
        evaluationGrids[ii]["pressure"] = pressure / ps

    return evaluationGrids


def compute_HRIR(ii, evaluationGrids, frequencies, reference, speedOfSound):
    """Compute HRIR from HRTF by means of the inverse Fourier transform"""
    # check if the frequency vector has the correct format
    if any(np.abs(np.diff(frequencies, 2)) > .1) or frequencies[0] < .1:
        raise ValueError(
            ('The frequency vector must go from f_1 > 0 to'
             'f_2 (half the sampling rate) in equidistant steps.'))

    if not reference:
        raise ValueError('HRIRs can only be computet if reference=true')

    pressure = evaluationGrids[ii]["pressure"]

    fs = round(2*frequencies[-1])

    # add 0 Hz bin
    pressure = np.concatenate((np.ones((pressure.shape[0],
                               pressure.shape[1], 1)), pressure), axis=2)
    # make fs/2 real
    pressure[:, :, -1] = np.abs(pressure[:, :, -1])
    # ifft (take complex conjugate because sign conventions differ)
    hrir = np.fft.irfft(np.conj(pressure))

    # shift 30 cm to make causal
    # (path differences between the origin and the ear are usually
    # smaller than 30 cm but numerical HRIRs show stringer pre-ringing)
    n_shift = int(np.round(.30 / (1/fs * speedOfSound)))
    hrir = np.roll(hrir, n_shift, axis=2)

    return hrir, fs


def write_to_sofa(ii, evaluationGrids, Mesh2HRTF_version,
                  frequencies, numSources, sourceCenter, fs=None, hrir=None,
                  type='HRTF'):
    """Write complex pressure or impulse responses as SOFA file."""

    # get source coordinates in spherical convention
    xyz = evaluationGrids[ii]["nodes"][:, 1:4]

    radius = np.sqrt(xyz[:, 0]**2 + xyz[:, 1]**2 + xyz[:, 2]**2)
    z_div_r = np.where(radius != 0, xyz[:, 2] / radius, 0)
    elevation = 90 - np.arccos(z_div_r) / np.pi * 180
    azimuth = np.mod(np.arctan2(xyz[:, 1], xyz[:, 0]), 2 * np.pi) / np.pi * 180

    source_position = np.concatenate(
        (azimuth[..., np.newaxis],
         elevation[..., np.newaxis],
         radius[..., np.newaxis]), axis=-1)

    # Save as SOFA file
    path = os.path.join('Output2HRTF',
                        f'{type}_{evaluationGrids[ii]["name"]}.sofa')

    # Need to delete it first if file already exists
    if os.path.exists(path):
        os.remove(path)

    # create empty SOFA object
    if type == 'HRTF':
        convention = "SimpleFreeFieldHRTF" if numSources == 2 else "GeneralTF"
    elif type == 'HRIR':
        convention = "SimpleFreeFieldHRIR" if numSources == 2 else "GeneralFIR"

    sofa = sf.Sofa(convention)

    # write meta data
    sofa.GLOBAL_ApplicationName = 'Mesh2HRTF'
    sofa.GLOBAL_ApplicationVersion = Mesh2HRTF_version
    sofa.GLOBAL_History = "numerically simulated data"

    # Source and receiver data
    sofa.SourcePosition = source_position
    sofa.SourcePosition_Units = "degree, degree, meter"
    sofa.SourcePosition_Type = "spherical"

    sofa.ReceiverPosition = sourceCenter
    sofa.ReceiverPosition_Units = "meter"
    sofa.ReceiverPosition_Type = "cartesian"

    # HRTF/HRIR data
    if type == 'HRTF':
        sofa.Data_Real = np.real(evaluationGrids[ii]["pressure"])
        sofa.Data_Imag = np.imag(evaluationGrids[ii]["pressure"])
        sofa.N = frequencies
    elif type == 'HRIR':
        sofa.Data_IR = hrir
        sofa.Data_SamplingRate = fs
        sofa.Data_Delay = np.zeros((1, hrir.shape[1]))

    # Save
    sf.write_sofa(path, sofa)


def Output2HRTF_Load(foldername, filename, numFrequencies):
    """
    Load results of the BEM calculation.

    Parameters
    ----------
    foldername : string
        The folder from which the data is loaded. The data to be read is
        located in the folder be.out inside NumCalc/source_*
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
    numFrequencies : int
        the number of simulated frequencies

    Returns
    -------
    data : numpy array
        Pressure or velocity values of shape (numFrequencies, numEntries)
    frequencies : numpy array
        The frequencies in Hz
    """

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
                        numDatalines = sum(1 for ll in line) + 2
                        break

    # ------------------------------load data----------------------------------
    data = np.zeros((numFrequencies, numDatalines), dtype=complex)
    frequency = np.zeros(numFrequencies)

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
            data[ii, :] = np.nan
            frequency[ii] = np.nan

    return data, frequency


def get_number_of_frequencies(path):
    """
    Read number of simulated frequency steps from Info.txt

    Parameters
    ----------
    path : str
        path of Info.txt

    Returns
    -------
    frequency_steps : int
        number of simulated frequency steps
    """

    key = 'Frequency Steps: '

    # read info file
    with open(path) as f:
        lines = f.readlines()

    # find line containing the number of frequency steps
    for line in lines:
        if line.startswith(key):
            break

    return int(line.strip()[len(key):])


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
            data = np.zeros((nSteps, 6), dtype=int)

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
