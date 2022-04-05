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
import glob
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
    numFrequencies = _get_number_of_frequencies('Info.txt')

    # get the evaluation grids
    evaluationGrids, evaluationGridsNumNodes = \
        _read_nodes_and_elements(data='EvaluationGrids')

    # get the object mesh
    objectMeshes, objectMeshesNumNodes = \
        _read_nodes_and_elements(data='ObjectMeshes')

    # # Read computational effort
    # print('\n Loading computational effort data ...')
    # computationTime = []
    # for source in range(numSources):
    #     for file in Path(os.path.join('NumCalc',
    #                                   f'source_{source+1}')).glob('NC*.out'):
    #         tmpFilename = os.path.join(file)
    #         tmp = _read_computation_time(tmpFilename)
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
    pressure, frequencies = _read_pressure(
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

        pressure, frequencies = _read_pressure(
            numSources, numFrequencies, data='pEvalGrid')

    # save to struct
    cnt = 0
    for ii in range(len(evaluationGrids)):
        evaluationGrids[ii]["pressure"] = pressure[
            cnt:cnt+evaluationGrids[ii]["num_nodes"], :, :]

        cnt = cnt + evaluationGrids[ii]["num_nodes"]

    del ii, pressure, cnt

    # process BEM data for writing HRTFs and HRIRs to SOFA files
    for ii in range(len(evaluationGrids)):
        eval_grid_name = evaluationGrids[ii]["name"]
        eval_grid_pressure = evaluationGrids[ii]["pressure"]
        eval_grid_xzy = evaluationGrids[ii]["nodes"][:, 1:4]

        print(f'\nSaving HRTFs for EvaluationGrid {eval_grid_name} ...\n')

        # get pressure as SOFA object (all following steps are run on SOFA
        # objects. This way they are available to other users as well)
        sofa = _get_sofa_object(
            eval_grid_pressure, eval_grid_xzy, "cartesian", sourceCenter,
            "HRTF", Mesh2HRTF_version, frequencies=frequencies)

        # reference to sound pressure at the center of the head
        if reference:
            sofa = reference_HRTF(sofa, sourceType, sourceArea, speedOfSound,
                                  densityOfAir, mode="min")

        # write HRTF data to SOFA file
        sf.write_sofa(os.path.join(
            'Output2HRTF', f'HRTF_{eval_grid_name}.sofa'), sofa)

        # calculate and write HRIRs
        if computeHRIRs:
            if not reference:
                raise ValueError("Computing HRIRs requires prior referencing")

            # calculate shift value (equivalent to a 30 cm shift)
            fs = round(2*sofa.frequencies[-1])
            n_shift = int(np.round(.30 / (1/fs * speedOfSound)))

            sofa = compute_HRIR(sofa, n_shift)
            sf.write_sofa(os.path.join(
                'Output2HRTF', f'HRIR_{eval_grid_name}.sofa'), sofa)

    print('Done\n')


def reference_HRTF(sofa, sourceType, sourceArea, speedOfSound, densityOfAir,
                   mode="min"):
    """
    Reference HRTF to the sound pressure in the center of the head. After
    referencing the sound pressure approaches 1 (0 dB) for low frequencies.

    Parameters
    ----------
    sofa : sofar Sofa object, str
       Sofa object containing the sound pressure or filename of a SOFA file to
       be loaded. SOFA object/file must be of the convention
       SimpleFreeFieldHRTF or GeneralTF
    sourceType : str
        The referencing depends on the source type used for simulating the
        sound pressure. Can be "Both ears", "Left ear", "Right ear",
        "Point source", or "Plane wave"
    sourceArea : array like
        The area of the source is required if `sourceType` is "Both ears" in
        which case `sourceAre` is a list with two values or if `sourceType` is
        "Left ear" or "Right ear" in which case `sourceArea` is a list with one
        value.
    speedOfSound : number
        The speed of sound in m/s
    densityOfAir : number
        The density of air in kg / m**3
    mode : str, optional
        Pass "min", "max", or "mean" to reference to the minmum, maximum, or
        mean radius in `sofa.SourcePosition`. Pass "all" to normalize each
        HRTF to the radius of the corresponding source.

    Returns
    -------
    sofa : sofar Sofa.object
        A copy of the input data with referenced sound pressure
    """

    if isinstance(sofa, str):
        sofa = sf.read_sofa(sofa)
    else:
        sofa = sofa.copy()

    if sofa.GLOBAL_SOFAConventions not in ["SimpleFreeFieldHRTF", "GeneralTF"]:
        raise ValueError(("Sofa object must have the conventions "
                          "SimpleFreeFieldHRTF or GeneralTF"))

    if sofa.SourcePosition_Type == "spherical":
        radius = sofa.SourcePosition[:, 2]
    else:
        radius = np.sqrt(sofa.SourcePosition[:, 0]**2 +
                         sofa.SourcePosition[:, 1]**2 +
                         sofa.SourcePosition[:, 2]**2)

    pressure = sofa.Data_Real + 1j * sofa.Data_Imag
    frequencies = sofa.N

    # distance of source positions from the origin
    if mode == "min":
        r = np.min(radius)
    elif mode == "mean":
        r = np.mean(radius)
    elif mode == "max":
        r = np.max(radius)
    else:
        r = radius[..., np.newaxis, np.newaxis]

    if sourceType in {'Both ears', 'Left ear', 'Right ear'}:

        volumeFlow = 0.1 * np.ones(pressure.shape)
        if 'sourceArea':
            # has to be fixed for both ears....
            for nn in range(len(sourceArea)):
                volumeFlow[:, nn, :] = \
                    volumeFlow[:, nn, :] * sourceArea[nn]

        # point source in the origin evaluated at r
        # eq. (6.71) in: Williams, E. G. (1999). Fourier Acoustics.
        ps = -1j * densityOfAir * 2 * np.pi * frequencies * \
            volumeFlow / (4 * np.pi) * \
            np.exp(1j * 2 * np.pi * frequencies / speedOfSound * r) / r

    elif sourceType == 'Point source':

        amplitude = 0.1  # hard coded in Mesh2HRTF
        ps = amplitude * \
            np.exp(1j * 2 * np.pi * frequencies /
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
    pressure /= ps
    sofa.Data_Real = np.real(pressure)
    sofa.Data_Imag = np.imag(pressure)

    return sofa


def compute_HRIR(sofa, n_shift):
    """
    Compute HRIR from HRTF by means of the inverse Fourier transform.

    This requires the following:

    1. The HRTFs contained in `sofa` have been referenced using
       `reference_HRTFs()`.
    2. HRTF must be available for frequencies f_0 to f_1 in constant steps.
       f_0 must be > 0 Hz and f_1 is assumed to be half the sampling rate.

    HRIRs are computed with the following steps

    1. Add data for 0 Hz. The HRTF at 0 Hz is 1 (0 dB) by definition because
       the HRTF describes the filtering of incoming sound by the human
       anthropometry (which does not change the sound at 0 Hz).
    2. Only the absolute value for the data at half the sampling rate is used.
       Otherwise the next step would produce complex output
    3. The HRTF spectrum is mirrored and the HRIR is obtained through the
       inverse Fourier transform
    4. The HRIRs are circularly shifted by `n_shift` samples to enforce a
       causal system. A good shift value for a sampling rate of 44.1 kHz might
       be between 20 and 40 samples.

    Parameters
    ----------
    sofa : sofar Sofa object, str
       Sofa object containing the sound pressure or filename of a SOFA file to
       be loaded. SOFA object/file must be of the convention
       SimpleFreeFieldHRTF or GeneralTF
    n_shift : int
        Amount the HRIRs are shifted to enforce a causal system.

    Returns
    -------
    sofa : sofar Sofa.object
        HRIRs

    Notes
    -----
    HRIRs for different sampling rates can be generated from a single SOFA file
    if discarding or adding some data.
    """

    if isinstance(sofa, str):
        sofa = sf.read_sofa(sofa)
    else:
        sofa = sofa.copy()

    if sofa.GLOBAL_SOFAConventions not in ["SimpleFreeFieldHRTF", "GeneralTF"]:
        raise ValueError(("Sofa object must have the conventions "
                          "SimpleFreeFieldHRTF or GeneralTF"))

    # check if the frequency vector has the correct format
    frequencies = sofa.N
    if any(np.abs(np.diff(frequencies, 2)) > .1) or frequencies[0] < .1:
        raise ValueError(
            ('The frequency vector must go from f_1 > 0 to'
             'f_2 (half the sampling rate) in equidistant steps.'))

    pressure = sofa.Data_Real + 1j * sofa.Data_Imag

    fs = round(2*frequencies[-1])

    # add 0 Hz bin
    pressure = np.concatenate((np.ones((pressure.shape[0],
                               pressure.shape[1], 1)), pressure), axis=2)
    # make fs/2 real
    pressure[:, :, -1] = np.abs(pressure[:, :, -1])
    # ifft (take complex conjugate because sign conventions differ)
    hrir = np.fft.irfft(np.conj(pressure))

    # shift to make causal
    # (path differences between the origin and the ear are usually
    # smaller than 30 cm but numerical HRIRs show stringer pre-ringing)
    hrir = np.roll(hrir, n_shift, axis=-1)

    sofa = _get_sofa_object(
        hrir, sofa.SourcePosition, sofa.SourcePosition_Type,
        sofa.ReceiverPosition, "HRIR", sofa.GLOBAL_ApplicationVersion,
        sampling_rate=fs)

    return sofa


def merge_sofa_files(left, right, join="both"):
    """
    Join HRTFs and HRIRs from separate SOFA-files containing the left and right
    ear data. The joined data is to SOFA files to the directory/directories
    containing the left ear data. The file names is extended by ``"_joined"``.

    Parameters
    ----------
    left, right : str
        Two strings that specify what SOFA files to load for the left and right
        ear. There are two options

        1. `left` and `right` are full paths to SOFA files, i.e., the end with
           ".sofa"
        2. `left` and `right` give the folders that contain the SOFA files in
           the subfolders "Output2HRTF". In this case `left` and `right` can
           contain an asterisk to process data in multiple folders. E.g., if
           `left` is ``"some/path/*_left"`` and `right` is
           ``"some/path/*_right"`` all SOFA files in the matching folders
           will be joined. Note that the Output2HRTF folder pairs given be
           `left` and `right` must contain SOFA files with identical names.
    join : str
        Specifies what data is joined. ``"HRTF"`` joins only the HRTFs,
        ``"HRIR"`` joins only the HRIRs. ``"both"`` joins all data. This is
        only required if `left` and `right` contain folder names.
    """

    # join single SOFA files
    if left.lower().endswith(".sofa") and right.lower().endswith(".sofa"):

        # filename of joined data
        head, tail = os.path.split(left)
        tail = tail[:-len(".sofa")] + "_joined.sofa"

        # join data
        _merge_sofa_files(left, right, os.path.join(head, tail))

    # join SOFA files contained in one or more folders
    else:
        # get all directories containing SOFA files
        lefts = glob.glob(left)
        rights = glob.glob(right)

        if len(lefts) != len(rights):
            raise ValueError(("The umber of directories found with glob.glob()"
                              f" does not match for {left} and {right}"))

        # loop directories
        for L, R in zip(lefts, rights):

            # check if Output2HRTF folder exists
            L = os.path.join(L, "Output2HRTF")
            R = os.path.join(R, "Output2HRTF")

            if not os.path.isdir(L) or not os.path.isdir(R):
                raise ValueError(f"Directory {L} and/or {R} does not exist")

            # check which data to join
            join = ["HRIR", "HRTF"] if join == "both" else [join]

            for data in join:
                # get and check all SOFA files in Output2HRTF folder
                Ls = glob.glob(os.path.join(L, data + "*.sofa"))
                Rs = glob.glob(os.path.join(R, data + "*.sofa"))

                if len(Ls) != len(Rs):
                    raise ValueError((
                        "The umber of sofa files found with glob.glob()"
                        f" does not match for {L} and {R}"))

                # loop all SOFA files
                for left_file, right_file in zip(Ls, Rs):

                    # check file names
                    if (os.path.basename(left_file)
                            != os.path.basename(right_file)):
                        raise ValueError((
                            "Found mismatching. Each Output2HRTF folder must "
                            "contain SOFA files with the same names. Error for"
                            f" {left_file} and {right_file}"))

                    # join
                    merge_sofa_files(left_file, right_file)


def _merge_sofa_files(left, right, savename):
    """read two sofa files, join the data, and save the result"""

    left = sf.read_sofa(left)
    right = sf.read_sofa(right)

    # join Data
    if left.GLOBAL_DataType.startswith("TF"):

        # check if data can be joined
        if left.N.size != right.N.size or \
                np.any(np.abs(left.N - right.N)) > 1e-6:
            raise ValueError(("Number of frequencies or frequencies do not "
                              f"agree for {left} and {right}"))

        # join data
        left.Data_Real = np.concatenate(
            (left.Data_Real, right.Data_Real), axis=1)
        left.Data_Imag = np.concatenate(
            (left.Data_Imag, right.Data_Imag), axis=1)

    elif left.GLOBAL_DataType.startswith("FIR"):

        # check if data can be joined
        if left.get_dimension("N") != right.get_dimension("N") or \
                left.Data_SamplingRate != right.Data_SamplingRate:
            raise ValueError(("Number of samples or sampling rates do not"
                              f"agree for {left} and {right}"))

        # join data
        left.Data_IR = np.concatenate(
            (left.Data_IR, right.Data_IR), axis=1)
        left.Data_Delay = np.zeros((1, left.Data_IR.shape[1]))
    else:
        raise ValueError("Joining only works for DataTypes 'TF' and 'FIR'")

    left.ReceiverPosition = np.concatenate(
        (np.atleast_2d(left.ReceiverPosition),
         np.atleast_2d(right.ReceiverPosition)), axis=0)

    # write joined data to disk
    sf.write_sofa(savename, left)


def _read_nodes_and_elements(data):
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


def _read_pressure(numSources, numFrequencies, data):
    """Read the sound pressure on the object meshes or evaluation grid."""
    pressure = []

    if not (data == 'pBoundary' or data == 'pEvalGrid'):
        raise ValueError('data must be pBoundary or pEvalGrid!')

    print('\n Loading %s data ...' % data)
    for source in range(numSources):
        print('\n    Source %d ...' % (source+1))

        tmpFilename = os.path.join('NumCalc', f'source_{source+1}', 'be.out')
        tmpPressure, frequencies = _output_to_hrtf_load(
            tmpFilename, data, numFrequencies)

        pressure.append(tmpPressure)

    pressure = np.transpose(np.array(pressure), (2, 0, 1))

    return pressure, frequencies


def _get_sofa_object(data, source_position, source_position_type,
                     receiver_position, mode, Mesh2HRTF_version,
                     frequencies=None, sampling_rate=None):
    """
    Write complex pressure or impulse responses to a SOFA object.

    Parameters
    ----------
    data : numpy array
        The data as an array of shape (MRE)
    evaluation_grid : numpy array
        The evaluation grid in Cartesian coordinates as an array of shape (MC)
    receiver_position : numpy array
        The position of the receivers (ears) in Cartesian coordinates
    mode : str
        "HRTF" to save HRTFs, "HRIR" to save HRIRs
    Mesh2HRTF_version : str
    frequencies : numpy array
        The frequencies at which the HRTFs were calculated. Required if mode is
        "HRTF"
    sampling_rate :
        The sampling rate. Required if mode is "HRIR"

    Returns
    -------
    sofa : sofar.Sofa object
    """

    # get source coordinates in spherical convention
    if source_position_type == "cartesian":
        xyz = source_position

        radius = np.sqrt(xyz[:, 0]**2 + xyz[:, 1]**2 + xyz[:, 2]**2)
        z_div_r = np.where(radius != 0, xyz[:, 2] / radius, 0)
        elevation = 90 - np.arccos(z_div_r) / np.pi * 180
        azimuth = np.mod(np.arctan2(xyz[:, 1], xyz[:, 0]), 2 * np.pi) \
            / np.pi * 180

        source_position = np.concatenate(
            (azimuth[..., np.newaxis],
             elevation[..., np.newaxis],
             radius[..., np.newaxis]), axis=-1)

    # get number of sources from data
    numSources = data.shape[1]

    # create empty SOFA object
    if mode == "HRTF":
        convention = "SimpleFreeFieldHRTF" if numSources == 2 else "GeneralTF"
    else:
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

    sofa.ReceiverPosition = receiver_position
    sofa.ReceiverPosition_Units = "meter"
    sofa.ReceiverPosition_Type = "cartesian"

    # HRTF/HRIR data
    if mode == "HRTF":
        sofa.Data_Real = np.real(data)
        sofa.Data_Imag = np.imag(data)
        sofa.N = frequencies
    else:
        sofa.Data_IR = data
        sofa.Data_SamplingRate = sampling_rate
        sofa.Data_Delay = np.zeros((1, data.shape[1]))

    return sofa


def _output_to_hrtf_load(foldername, filename, numFrequencies):
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


def _get_number_of_frequencies(path):
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


def _read_computation_time(filename):
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
