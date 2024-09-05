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
    4. use :py:func:`~mesh2hrtf.reference_hrtfs` and
       :py:func:`~mesh2hrtf.compute_hrirs` to save the results to SOFA files

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

    # Load EvaluationGrid data
    if not len(evaluationGrids) == 0:
        print('\nLoading data for the evaluation grids ...')

        pressure, _ = _read_numcalc_data(
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
            sofa = m2h.reference_hrtfs(
                sofa, params["sourceType"], params["sourceArea"],
                params["speedOfSound"], params["densityOfMedium"], mode="min")

        # Resampling spectrum to regular sampling, if necessary
        if num_octaves_decim > 0:
            pressure = sofa.Data_Real + 1j*sofa.Data_Imag
            frequencies = sofa.N
            sampling_rate = round(2*sofa.N[-1])

            num_octaves_decim = params["num_octaves_decim"]

            pressure, frequencies, inds_lin_freq = _convert_to_regular_samp_freq(
                                                    pressure, frequencies)

            inds_correct_phase = inds_lin_freq

            pressure = _correct_phase(pressure, inds_correct_phase)

            sofa.N = frequencies.flatten()
            sofa.Data_Real = np.real(pressure)
            sofa.Data_Imag = np.imag(pressure)

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

            sofa = m2h.compute_hrirs(sofa, n_shift)
            sf.write_sofa(os.path.join(
                folder, 'Output2HRTF', f'HRIR_{grid}.sofa'), sofa)

    print('Done\n')


def _read_nodes_and_elements(folder, objects=None):
    """
    Read the nodes and elements of the evaluation grids or object meshes.

    Parameters
    ----------
    folder : str
        Folder containing the object. Must end with EvaluationGrids or
        Object Meshes
    objects : str, options
        Name of the object. The default ``None`` reads all objects in folder

    Returns
    -------
    grids : dict
        One item per object (with the item name being the object name). Each
        item has the sub-items `nodes`, `elements`, `num_nodes`, `num_elements`
    gridsNumNodes : int
        Number of nodes in all grids
    """
    # check input
    if os.path.basename(folder) not in ['EvaluationGrids', 'ObjectMeshes']:
        raise ValueError('folder must be EvaluationGrids or ObjectMeshes!')

    if objects is None:
        objects = os.listdir(folder)
        # discard hidden folders that might occur on Mac OS
        objects = [o for o in objects if not o.startswith('.')]
    elif isinstance(objects, str):
        objects = [objects]

    grids = {}
    gridsNumNodes = 0

    for grid in objects:
        tmpNodes = np.loadtxt(os.path.join(
            folder, grid, 'Nodes.txt'),
            delimiter=' ', skiprows=1)

        tmpElements = np.loadtxt(os.path.join(
            folder, grid, 'Elements.txt'),
            delimiter=' ', skiprows=1)

        grids[grid] = {
            "nodes": tmpNodes,
            "elements": tmpElements,
            "num_nodes": tmpNodes.shape[0],
            "num_elements": tmpElements.shape[0]}

        gridsNumNodes += grids[grid]['num_nodes']

    return grids, gridsNumNodes


def _read_numcalc_data(numSources, numFrequencies, folder, data):
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
    numDatalines = None
    with open(current_file) as file:
        line = csv.reader(file, delimiter=' ', skipinitialspace=True)
        for idx, li in enumerate(line):
            # read number of data points and head lines
            if len(li) == 2 and not li[0].startswith("Mesh"):
                numDatalines = int(li[1])

            # read starting index
            elif numDatalines and len(li) > 2:
                start_index = int(li[0])
                break

    # ------------------------------load data----------------------------------
    dtype = complex if filename.startswith("p") else float
    data = np.zeros((numFrequencies, numDatalines), dtype=dtype)

    for ii in range(numFrequencies):
        tmpData = []
        current_file = os.path.join(foldername, 'be.%d' % (ii+1), filename)
        with open(current_file) as file:

            line = csv.reader(file, delimiter=' ', skipinitialspace=True)

            for li in line:

                # data lines have 3 ore more entries
                if len(li) < 3 or li[0].startswith("Mesh"):
                    continue

                if filename.startswith("p"):
                    tmpData.append(complex(float(li[1]), float(li[2])))
                elif filename == "vBoundary":
                    tmpData.append(np.abs(complex(float(li[1]), float(li[2]))))
                elif filename == "vEvalGrid":
                    tmpData.append(np.sqrt(
                        np.abs(complex(float(li[1]), float(li[2])))**2 +
                        np.abs(complex(float(li[3]), float(li[4])))**2 +
                        np.abs(complex(float(li[5]), float(li[6])))**2))

        data[ii, :] = tmpData if tmpData else np.nan

    return data, np.arange(start_index, numDatalines + start_index)



def _convert_to_regular_samp_freq(pressure, frequencies):
    """
    Converts irregularly sampled frequency axis into regular intervals via 
    interpolation.
    """
    print('\n \n Irregular frequency sampling identified. \n \n'
          'Converting to regular frequency sampling... \n \n')

    freqs_decim = np.array(frequencies).flatten()
    minFrequency = freqs_decim[0]
    frequencyStepSize = freqs_decim[1] - freqs_decim[0]
    frequencySteps = freqs_decim[-1]/frequencyStepSize

    # get all frequencies to be calculated
    frequencies = np.array([ff*frequencyStepSize+minFrequency
                            for ff in range(int(frequencySteps))])
    pressure_decim = np.squeeze(pressure)

    # Interpolate spectra back to regular sampling in frequency
    pressure_decim_mag = np.array([np.interp(frequencies, freqs_decim, np.abs(
                    pressure_decim)[i]) for i in range(pressure_decim.shape[0])])

    pressure_decim_phase = np.array([np.interp(frequencies, freqs_decim, 
                                    np.angle(pressure_decim)[i]) 
                                    for i in range(pressure_decim.shape[0])])

    # # Zeroing phase for high frequencies in the interpolated region
    inds_lin_freq = ((np.diff(np.diff(np.concatenate((freqs_decim, 
                        np.array([freqs_decim[-1]*2, freqs_decim[-1]*4])))))) == 0)
    
    pressure_decim_phase[:, np.sum(np.array(inds_lin_freq)) : ] = 0

    pressure = np.zeros([pressure.shape[0], pressure.shape[1], 
                         pressure_decim_mag.shape[1]], dtype = "complex_")

    pressure[:,0,:] = np.multiply(pressure_decim_mag, np.exp(1j*pressure_decim_phase))  

    inds_lin_freq = (frequencies < frequencies[np.sum(np.array(inds_lin_freq))])

    return pressure, frequencies, inds_lin_freq


def _correct_phase(pressure, inds_correct_phase):
    """
    Corrects the phase for high frequencies that were simulated using irregular
    frequency sampling using the overall group delay of the linearly sampled region.
    """
    new_phase_all = np.zeros(shape=pressure.shape)

    for Direction_ind in range(pressure.shape[0]):
        for ear_ind in range(pressure.shape[1]):

            phase = np.unwrap(np.angle(pressure[Direction_ind, ear_ind, :]))

            phase_diff = np.diff(phase[inds_correct_phase])
            phase_diff_mean = np.average(phase_diff)

            phase[~inds_correct_phase] = np.linspace(
                        start = (phase[inds_correct_phase])[-1] + phase_diff_mean, 
                        stop = (phase[inds_correct_phase])[-1] + phase_diff_mean*np.sum(~inds_correct_phase), 
                        num = np.sum(~inds_correct_phase))

            new_phase_all[Direction_ind, ear_ind, :] = phase
            
    pressure = np.multiply(np.abs(pressure), np.exp(1j*new_phase_all))
    
    return pressure