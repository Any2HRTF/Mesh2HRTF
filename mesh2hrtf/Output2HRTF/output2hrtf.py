"""
Python tools for Mesh2HRTF including functions to generate SOFA files
containing the HRTF/HRIR data, merge SOFA files containing data for the left
and right ear and generate evaluation grids.
"""
import os
import csv
import warnings
import numpy as np
import json
import sofar as sf
import mesh2hrtf as m2h


def output2hrtf(folder=None):
    """
    Process NumCalc output and write data to disk.

    Processing the data is done in the following steps

    1. Read project parameter `from parameters.json`
    2. use :py:func:`~write_output_report` to parse files in
       project_folder/NumCalc/source_*/NC*.out, write project report to
       project_folder/Output2HRTF/report_source_*.csv. Raise a warning if any
       issues were detected and write report_issues.txt to the same folder
    3. Read simulated pressures from project_folder/NumCalc/source_*/be.out.
       This and the following steps are done, even if an issue was detected in
       the previous step
    4. use :py:func:`~reference_hrtf` and :py:func:`~compute_hrir` to save the
       results to SOFA files

    Parameters
    ----------
    folder : str, optional
        The path of the Mesh2HRTF project folder, i.e., the folder containing
        the subfolders EvaluationsGrids, NumCalc, and ObjectMeshes. The
        default, ``None`` uses the current working directory.
    """

    # check input
    if folder is None:
        folder = os.getcwd()

    # check and load parameters, required parameters are:
    # Mesh2HRTF_version, reference, computeHRIRs, speedOfSound, densityOfAir,
    # numSources, sourceType, sourceCenter, sourceArea,
    # numFrequencies, frequencies
    params = os.path.join(folder, "parameters.json")
    if not os.path.isfile(params):
        raise ValueError((
            f"The folder {folder} is not a valid Mesh2HRTF project. "
            "It must contain the file 'parameters.json'"))

    with open(params, "r") as file:
        params = json.load(file)

    # output directory
    if not os.path.exists(os.path.join(folder, 'Output2HRTF')):
        os.makedirs(os.path.join(folder, 'Output2HRTF'))

    # write the project report and check for issues
    print('\n Writing the project report ...')
    found_issues, report = m2h.write_output_report(folder)

    if found_issues:
        warnings.warn(report)

    # get the evaluation grids
    print('\n Loading the EvaluationGrids ...')
    evaluationGrids, _ = _read_nodes_and_elements(
        os.path.join(folder, 'EvaluationGrids'))

    # get the object mesh
    print('\n Loading the ObjectMeshes ...')
    objectMeshes, _ = _read_nodes_and_elements(
        os.path.join(folder, 'ObjectMeshes'))

    # Load ObjectMesh data
    pressure, _ = _read_pressure(
        params["numSources"], params["numFrequencies"], folder, 'pBoundary')

    print('\nSaving ObjectMesh data ...')
    cnt = 0
    for mesh in objectMeshes:
        elements = objectMeshes[mesh]["elements"]
        element_data = pressure

        for jj in range(pressure.shape[1]):
            element_data[:, jj, :] = element_data[
                cnt:cnt+elements.shape[0], jj, :]

        file = open(os.path.join(
            folder, "Output2HRTF",
            "ObjectMesh_" + mesh + ".npz"), "w")
        np.savez_compressed(
            file.name, frequencies=params["frequencies"],
            pressure=element_data)
        file.close()

        cnt = cnt + elements.shape[0]

    del pressure, elements, mesh, jj, cnt, element_data

    # Load EvaluationGrid data
    if not len(evaluationGrids) == 0:
        print('\nLoading data for the evaluation grids ...')

        pressure, _ = _read_pressure(
            params["numSources"], params["numFrequencies"],
            folder, 'pEvalGrid')

    # save to struct
    cnt = 0
    for grid in evaluationGrids:
        evaluationGrids[grid]["pressure"] = pressure[
            cnt:cnt+evaluationGrids[grid]["num_nodes"], :, :]

        cnt = cnt + evaluationGrids[grid]["num_nodes"]

    del grid, pressure, cnt

    # process BEM data for writing HRTFs and HRIRs to SOFA files
    for grid in evaluationGrids:

        print(f'\nSaving HRTFs for EvaluationGrid "{grid}" ...\n')

        # get pressure as SOFA object (all following steps are run on SOFA
        # objects. This way they are available to other users as well)
        sofa = _get_sofa_object(
            evaluationGrids[grid]["pressure"],
            evaluationGrids[grid]["nodes"][:, 1:4], "cartesian",
            params["sourceCenter"], "HRTF", params["Mesh2HRTF_Version"],
            frequencies=params["frequencies"])

        # reference to sound pressure at the center of the head
        if params["reference"]:
            sofa = reference_hrtf(
                sofa, params["sourceType"], params["sourceArea"],
                params["speedOfSound"], params["densityOfMedium"], mode="min")

        # write HRTF data to SOFA file
        sf.write_sofa(os.path.join(
            folder, 'Output2HRTF', f'HRTF_{grid}.sofa'), sofa)

        # calculate and write HRIRs
        if params["computeHRIRs"]:
            if not params["reference"]:
                raise ValueError("Computing HRIRs requires prior referencing")

            # calculate shift value (equivalent to a 30 cm shift)
            sampling_rate = round(2*sofa.N[-1])
            n_shift = int(np.round(
                .30 / (1/sampling_rate * params["speedOfSound"])))

            sofa = compute_hrir(sofa, n_shift)
            sf.write_sofa(os.path.join(
                folder, 'Output2HRTF', f'HRIR_{grid}.sofa'), sofa)

    print('Done\n')


def reference_hrtf(sofa, sourceType, sourceArea, speedOfSound, densityOfAir,
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


def compute_hrir(sofa, n_shift, sampling_rate=None):
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
    sampling_rate : int
        The sampling rate in Hz. The sampling rate can two times any frequency
        for which the HRTF was computed. The default ``None`` assumes the
        sampling rate to be two times the highest frequency for which the HRTF
        was computed.

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

    # get the HRTF
    pressure = sofa.Data_Real + 1j * sofa.Data_Imag

    if sampling_rate is None:
        # detect the sampling rate
        sampling_rate = round(2*frequencies[-1])
    else:
        # check if the specified sampling rate is valid
        idx = np.argmin(np.abs(2 * frequencies - sampling_rate))
        error = np.abs(2 * frequencies[idx] - sampling_rate)
        if np.abs(error) > 1e-6:
            raise ValueError((
                "The specified sampling rate is invalid. It must be two times "
                "any frequency for which the HRTF was computed."))

        # discard frequencies above sampling_rate / 2
        pressure = pressure[..., :idx + 1]
        frequencies = frequencies[:idx + 1]

    # add 0 Hz bin
    pressure = np.concatenate((np.ones((pressure.shape[0],
                               pressure.shape[1], 1)), pressure), axis=2)

    # make sampling_rate/2 real
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
        sampling_rate=sampling_rate)

    return sofa


def _read_nodes_and_elements(data):
    """
    Read the nodes and elements of the evaluation grids or object meshes.
    """
    if os.path.basename(data) not in ['EvaluationGrids', 'ObjectMeshes']:
        raise ValueError('data must be EvaluationGrids or ObjectMeshes!')

    grids = {}
    gridsList = os.listdir(data)
    gridsNumNodes = 0

    for grid in gridsList:
        tmpNodes = np.loadtxt(os.path.join(
            data, grid, 'Nodes.txt'),
            delimiter=' ', skiprows=1)

        tmpElements = np.loadtxt(os.path.join(
            data, grid, 'Elements.txt'),
            delimiter=' ', skiprows=1)

        grids[grid] = {
            "nodes": tmpNodes,
            "elements": tmpElements,
            "num_nodes": tmpNodes.shape[0]}

        gridsNumNodes += grids[grid]['num_nodes']

    return grids, gridsNumNodes


def _read_pressure(numSources, numFrequencies, folder, data):
    """Read the sound pressure on the object meshes or evaluation grid."""
    pressure = []

    if data not in ['pBoundary', 'pEvalGrid', 'vBoundary', 'vEvalGrid']:
        raise ValueError(
            'data must be pBoundary, pEvalGrid, vBoundary, or vEvalGrid')

    print('\n Loading %s data ...' % data)
    for source in range(numSources):
        print('\n    Source %d ...' % (source+1))

        tmpFilename = os.path.join(
            folder, 'NumCalc', f'source_{source+1}', 'be.out')
        tmpPressure, indices = _output_to_hrtf_load(
            tmpFilename, data, numFrequencies)

        pressure.append(tmpPressure)

    pressure = np.transpose(np.array(pressure), (2, 0, 1))

    return pressure, indices


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
        sofa.N = np.array(frequencies).flatten()
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
            The sound pressure on the evaluation grid
        vEvalGrid
            The sound velocity on the evaluation grid
    numFrequencies : int
        the number of simulated frequencies

    Returns
    -------
    data : numpy array
        Pressure or abs velocity values of shape (numFrequencies, numEntries)
    """

    # ---------------------check number of header and data lines---------------
    current_file = os.path.join(foldername, 'be.1', filename)
    with open(current_file) as file:
        line = csv.reader(file, delimiter=' ', skipinitialspace=True)
        for idx, li in enumerate(line):
            # read number of data points and head lines
            if idx == 2:
                if len(li) == 2:
                    numHeaderlines = 3
                    numDatalines = int(li[1])
                else:
                    raise ValueError((
                        f"error reading {current_file}. Expected three header "
                        "lines and the number of data points as second entry "
                        "in the third line."))
            # read starting index
            if idx == 3:
                start_index = int(li[0])

    # ------------------------------load data----------------------------------
    data = np.zeros((numFrequencies, numDatalines), dtype=complex)

    for ii in range(numFrequencies):
        tmpData = []
        current_file = os.path.join(foldername, 'be.%d' % (ii+1), filename)
        with open(current_file) as file:

            line = csv.reader(file, delimiter=' ', skipinitialspace=True)

            for li in line:

                if line.line_num <= numHeaderlines:
                    continue

                if filename.startswith("p"):
                    tmpData.append(complex(float(li[1]), float(li[2])))
                elif filename == "vBoundary":
                    tmpData.append(np.abs(complex(float(li[1]), float(li[2]))))
                elif filename == "vEvalGrid":
                    tmpData.append(np.sqrt(
                        np.abs(complex(float(li[1]), float(li[2])))**2 +
                        np.abs(complex(float(li[3]), float(li[4])))**2 +
                        np.abs(complex(float(li[5]), float(li[6])))**2
                        ))

        data[ii, :] = tmpData if tmpData else np.nan

    # make indices ------------------------------------------------------------
    indices = np.arange(start_index, numDatalines + start_index)

    return data, indices
