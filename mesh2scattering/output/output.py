import os
import warnings
import json
import mesh2scattering as m2s
import numpy as np
import pyfar as pf
import numpy as np
import glob
import sofar as sf


def angles2coords(
        azimuth, colatitude, radius=1., unit='rad') -> pf.Coordinates:
    """ ``data.cshape`` fits the cshape of ```coords``. Data get shifted
    through the ``coords`` Object around azimuth by ``shift_azimuth``.


    Parameters
    ----------
    azimuth : ndarray
        1D array or list of azimuth angles
    colatitude : ndarray
        1D array or list of incident angles
    radius : float, optional
        radius of the coordinate object, by default 1.
    unit : str, optional
        defines the unit of azimuth and colatitude, by default 'rad'

    Returns
    -------
    pf.Coordinates
        coordinate object ob chape (#azimuth, #colatitude) with radius radius.
    """
    azimuth = np.array(azimuth)
    colatitude = np.array(colatitude)
    if unit == 'deg':
        azimuth = azimuth * np.pi / 180.
        colatitude = colatitude * np.pi / 180.
    elif unit != 'rad':
        raise TypeError("Unknown Unit")
    phi, theta = np.meshgrid(azimuth, colatitude, indexing='ij')
    return pf.Coordinates(
        phi, theta, np.ones(phi.shape)*radius, 'sph')


def shift_data_coords(
        data, coords, shift_azimuth) -> pf.FrequencyData:
    """``data`` get shifted through the ``coords`` Object around azimuth by
    ``shift_azimuth``.

    Parameters
    ----------
    data : pf.FrequencyData
        data which need to be shifted. I has the cshape of (..., coords.cshape)
    coords : pf.Coordinates
        coordinates object of ``data```. It has the chape (#azimuth,
        #colatitude)
    shift_azimuth : float
        data get shifted by this azimuth angle in degree.

    Returns
    -------
    pf.FrequencyData
        shifted data object
    """
    # test input
    if not isinstance(data, pf.FrequencyData):
        raise TypeError(
            f'Data should be of type FrequencyData not {type(data)}')
    if not isinstance(coords, pf.Coordinates):
        raise TypeError(
            f'coords should be of type Coordinates not {type(coords)}')
    if not isinstance(shift_azimuth, (float, int)):
        raise TypeError(
            f'shift_azimuth should be of type float not {type(shift_azimuth)}')

    if shift_azimuth == 0:
        return data.copy()
    coords_ref = coords.copy()
    coords_cp = coords.copy()
    sph = coords_cp.get_sph(unit='deg')
    # shift azimuth by shift_azimuth in deg
    azimuth = np.remainder(sph[..., 0] + shift_azimuth, 360)
    coords_cp.set_sph(azimuth, sph[..., 1], sph[..., 2], unit='deg')
    xyz = coords_ref.get_cart()
    data_mask, _ = coords_cp.find_nearest_k(
        xyz[..., 0], xyz[..., 1], xyz[..., 2])
    data_mask = data_mask.flatten()
    freq = data.freq.copy()
    freq_new = freq[:, data_mask, ...]
    data_out = pf.FrequencyData(
        freq_new, data.frequencies)
    return data_out


def reshape_to_az_by_el(
        data: pf.FrequencyData, coords_in: pf.Coordinates,
        coords_out: pf.Coordinates, cdim: int = 0) -> (pf.FrequencyData):
    if cdim > 0:
        data.freq = np.moveaxis(data.freq, cdim, 0)
    freq_shape = list(coords_out.cshape)
    if len(data.cshape) > 1:
        for dim in data.cshape[1:]:
            freq_shape.append(dim)
    freq_shape.append(data.n_bins)
    freq = np.zeros(freq_shape, dtype=complex)
    data_in = data.freq
    xyz = coords_out.get_cart()
    index, _ = coords_in.find_nearest_k(xyz[..., 0], xyz[..., 1], xyz[..., 2])
    for iaz in range(coords_out.cshape[0]):
        res_data = data_in[index[iaz, :], ...]
        freq[iaz, ...] = res_data
    if cdim > 0:
        freq = np.moveaxis(freq, 0, cdim+1)
        freq = np.moveaxis(freq, 0, cdim+1)
    data_out = pf.FrequencyData(freq, data.frequencies)
    return data_out


def apply_symmetry_circular(
        data: pf.FrequencyData, coords_mic: pf.Coordinates,
        coords_inc: pf.Coordinates, coords_inc_out: pf.Coordinates):
    """apply symmetry for circular symmetrical surfaces.

    Parameters
    ----------
    data : pf.FrequencyData
        data which is rotated, cshape need to be (#theta_coords_inc,
        #coords_mic)
    coords_mic : pf.Coordinates
        Coordinate object from the receiver positions of the current reference
        plate of cshape (#theta_coords_inc)
    coords_inc : pf.Coordinates
        Coordinate object from the source positions of the reference of cshape
        (#coords_inc_reference.
    coords_inc_out : pf.Coordinates
        Coordinate object from the source positions of the sample of cshape
        (#coords_inc_sample).

    Returns
    -------
    pf.FrequencyData
        _description_
    """
    shape = [coords_inc_out.cshape[0], data.cshape[1], len(data.frequencies)]
    freq = np.empty(shape, dtype=complex)
    thetas = np.sort(np.array(list(set(
        np.round(coords_inc.get_sph(unit='deg')[..., 1], 5)))))
    for ii in range(coords_inc_out.cshape[0]):
        az = coords_inc_out[ii].get_sph(unit='deg')[0, 0]
        theta = coords_inc_out[ii].get_sph(unit='deg')[0, 1]
        data_in = data[np.abs(thetas-theta) < 1e-5, :]
        freq[ii, ...] = shift_data_coords(
            data_in, coords_mic, float(az)).freq.copy()
    data_out = pf.FrequencyData(freq, data.frequencies)
    return data_out


def apply_symmetry_mirror(data, coords_mic, incident_coords, mirror_axe=None):
    shape = list(data.cshape)
    shape.append(data.n_bins)
    shape[mirror_axe] = 2*shape[mirror_axe]-1
    index_min = int(shape[mirror_axe]/2)
    index_max = int(shape[mirror_axe])
    freq = np.empty(shape, dtype=complex)
    freq[:] = np.nan
    freq = np.moveaxis(freq, mirror_axe, -1)
    freq_in = np.moveaxis(data.freq, mirror_axe, -1)
    azimuths = incident_coords.get_sph()[:, 1, 0]
    max_aimuth = np.max(azimuths)
    elevations = incident_coords.get_sph()[0, :, 1]
    radius = np.median(incident_coords.get_sph()[:, :, 2])
    azimuths_new = []
    max_index = -1
    for iaz in range(index_max):
        if iaz > index_min:
            idx = index_max-iaz-1
            az = (max_aimuth-azimuths[idx]) * 2
            if azimuths[idx] + az > 2 * np.pi:
                max_index = iaz
                break
            data_in = shift_data_coords(data, coords_mic, az/np.pi*180).freq
            azimuths_new.append(azimuths[idx] + az)
        else:
            data_in = data.freq
            idx = iaz
            azimuths_new.append(azimuths[iaz])
        freq_in = np.moveaxis(data_in, mirror_axe, -1)
        freq[..., iaz] = freq_in[..., idx]
    if max_index > 0:
        freq = freq[..., :max_index]
    freq = np.moveaxis(freq, -1, mirror_axe)
    data_out = pf.FrequencyData(freq, data.frequencies)
    new_inc_coords = angles2coords(np.array(azimuths_new), elevations, radius)
    return data_out, new_inc_coords


def write_pattern(folder):
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

    if (not os.path.exists(os.path.join(folder, 'reference'))) \
            or (not os.path.exists(os.path.join(folder, 'sample'))):
        raise ValueError(
            "Folder need to contain reference and sample folders.")

    # read sample data
    evaluationGrids, params = read_numcalc(
        os.path.join(folder, 'sample'), False)

    # process BEM data for writing HRTFs and HRIRs to SOFA files
    for grid in evaluationGrids:
        print(f'\nWrite sample data "{grid}" ...\n')
        # get pressure as SOFA object (all following steps are run on SOFA
        # objects. This way they are available to other users as well)
        source_position = np.array(params["sourceCenter"])
        if source_position.shape[1] != 3:
            source_position = np.transpose(source_position)
        receiver_position = np.array(evaluationGrids[grid]["nodes"][:, 1:4])
        if receiver_position.shape[1] != 3:
            receiver_position = np.transpose(receiver_position)
        sofa = m2s.utils._get_sofa_object(
            evaluationGrids[grid]["pressure"],
            source_position,
            receiver_position,
            params["Mesh2HRTF_Version"],
            frequencies=params["frequencies"])

        sofa.GLOBAL_Title = folder.split(os.sep)[-1]

        # write HRTF data to SOFA file
        sf.write_sofa(os.path.join(
            folder, 'sample.pattern.sofa'), sofa)

    evaluationGrids, params = read_numcalc(
        os.path.join(folder, 'reference'), True)

    # process BEM data for writing HRTFs and HRIRs to SOFA files
    for grid in evaluationGrids:
        print(f'\nWrite sample data "{grid}" ...\n')
        # get pressure as SOFA object (all following steps are run on SOFA
        # objects. This way they are available to other users as well)
        # read source and receiver positions
        source_position_ref = np.array(params["sourceCenter"])
        if source_position_ref.shape[1] != 3:
            source_position_ref = np.transpose(source_position_ref)
        receiver_position_ref = np.array(
            evaluationGrids[grid]["nodes"][:, 1:4])
        if receiver_position_ref.shape[1] != 3:
            receiver_position_ref = np.transpose(receiver_position_ref)

        # apply symmetry of reference sample
        data = evaluationGrids[grid]["pressure"]
        data = np.swapaxes(data, 0, 1)
        data_out = apply_symmetry_circular(
            pf.FrequencyData(data, params["frequencies"]),
            _cart_coordinates(receiver_position_ref),
            _cart_coordinates(source_position_ref),
            _cart_coordinates(source_position))

        # create sofa file
        sofa = m2s.utils._get_sofa_object(
            data_out.freq,
            source_position,
            receiver_position_ref,
            params["Mesh2HRTF_Version"],
            frequencies=params["frequencies"])

        sofa.GLOBAL_Title = folder.split(os.sep)[-1]

        # write HRTF data to SOFA file
        sf.write_sofa(os.path.join(
            folder, 'reference.pattern.sofa'), sofa)

    print('Done\n')


def _cart_coordinates(xyz):
    return pf.Coordinates(xyz[:, 0], xyz[:, 1], xyz[:, 2])


def check_project(folder=None):
    r"""
    Generate project report from NumCalc output files.

    NumCalc (Mesh2HRTF's numerical core) writes information about the
    simulations to the files `NC*.out` located under `NumCalc/source_*`. The
    file `NC.out` exists if NumCalc was ran without the additional command line
    parameters ``-istart`` and ``-iend``. If these parameters were used, there
    is at least one `NC\*-\*.out`. If this is the case, information from
    `NC\*-\*.out` overwrites information from NC.out in the project report.

    .. note::

        The project reports are written to the files
        `Output2HRTF/report_source_*.csv`. If issues were detected, they are
        listed in `Output2HRTF/report_issues.csv`.

    The report contain the following information

    Frequency step
        The index of the frequency.
    Frequency in Hz
        The frequency in Hz.
    NC input
        Name of the input file from which the information was taken.
    Input check passed
        Contains a 1 if the check of the input data passed and a 0 otherwise.
        If the check failed for one frequency, the following frequencies might
        be affected as well.
    Converged
        Contains a 1 if the simulation converged and a 0 otherwise. If the
        simulation did not converge, the relative error might be high.
    Num. iterations
        The number of iterations that were required to converge
    relative error
        The relative error of the final simulation
    Comp. time total
        The total computation time in seconds
    Comp. time assembling
        The computation time for assembling the matrices in seconds
    Comp. time solving
        The computation time for solving the matrices in seconds
    Comp. time post-proc
        The computation time for post-processing the results in seconds


    Parameters
    ----------
    folder : str, optional
        The path of the Mesh2HRTF project folder, i.e., the folder containing
        the subfolders EvaluationsGrids, NumCalc, and ObjectMeshes. The
        default, ``None`` uses the current working directory.

    Returns
    -------
    found_issues : bool
        ``True`` if issues were found, ``False`` otherwise
    report : str
        The report or an empty string if no issues were found
    """

    if folder is None:
        folder = os.getcwd()

    # get sources and number of sources and frequencies
    sources = glob.glob(os.path.join(folder, "NumCalc", "source_*"))
    num_sources = len(sources)

    with open(os.path.join(folder, "parameters.json"), "r") as file:
        params = json.load(file)

    # sort source files (not read in correct order in some cases)
    nums = [int(source.split("_")[-1]) for source in sources]
    sources = np.array(sources)
    sources = sources[np.argsort(nums)]

    # parse all NC*.out files for all sources
    all_files, fundamentals, out, out_names = m2s.utils._parse_nc_out_files(
        sources, num_sources, params["numFrequencies"])

    return all_files, fundamentals, out, out_names


def merge_frequency_data(data_list):
    data_out = data_list[0].copy()
    frequencies = data_out.frequencies.copy()
    for idx in range(1, len(data_list)):
        data = data_list[idx]
        assert data_out.cshape == data.cshape
        frequencies = np.append(frequencies, data.frequencies)
        frequencies = np.array([i for i in set(frequencies)])
        frequencies = np.sort(frequencies)

        data_new = []
        for f in frequencies:
            if any(data_out.frequencies == f):
                freq_index = np.where(data_out.frequencies == f)
                data_new.append(data_out.freq[..., freq_index[0][0]])
            elif any(data.frequencies == f):
                freq_index = np.where(data.frequencies == f)
                data_new.append(data.freq[..., freq_index[0][0]])

        data_new = np.moveaxis(np.array(data_new), 0, -1)
        data_out = pf.FrequencyData(data_new, frequencies)
    return data_out


def read_numcalc(folder=None, is_ref=False):
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
    params = os.path.join(folder, '..', 'parameters.json')
    if not os.path.isfile(params):
        raise ValueError((
            f"The folder {folder} is not a valid Mesh2scattering project. "
            "It must contain the file 'parameters.json'"))

    with open(params, "r") as file:
        params = json.load(file)

    # get source positions
    source_coords = np.transpose(np.array(params['sourceCenter']))
    source_coords = pf.Coordinates(
        source_coords[..., 0], source_coords[..., 1], source_coords[..., 2])

    # output directory
    if not os.path.exists(os.path.join(folder, 'Output2HRTF')):
        os.makedirs(os.path.join(folder, 'Output2HRTF'))

    # write the project report and check for issues
    print('\n Writing the project report ...')
    found_issues, report = write_output_report(folder)

    if found_issues:
        warnings.warn(report)

    # get the evaluation grids
    evaluationGrids, _ = m2s.utils._read_nodes_and_elements(
        os.path.join(folder, 'EvaluationGrids'))

    # Load EvaluationGrid data
    if is_ref:
        xyz = np.array(params["sourceCenter"])
        coords = pf.Coordinates(xyz[..., 0], xyz[..., 1], xyz[..., 2])
        num_sources = np.sum(np.abs(coords.get_sph()[..., 0]) < 1e-12)
    else:
        num_sources = params["numSources"]

    if not len(evaluationGrids) == 0:
        pressure, _ = m2s.utils._read_numcalc_data(
            num_sources, params["numFrequencies"],
            folder, 'pEvalGrid')

    # save to struct
    cnt = 0
    for grid in evaluationGrids:
        evaluationGrids[grid]["pressure"] = pressure[
            cnt:cnt+evaluationGrids[grid]["num_nodes"], :, :]

        cnt = cnt + evaluationGrids[grid]["num_nodes"]

    receiver_coords = evaluationGrids[grid]["nodes"][:, 1:4]
    receiver_coords = pf.Coordinates(
        receiver_coords[..., 0], receiver_coords[..., 1],
        receiver_coords[..., 2])

    return evaluationGrids, params


def read_ram_estimates(folder: str):
    """
    Read estimated RAM consumption from Memory.txt.

    Note that the RAM consumption per frequency step can be estimated and
    written to `Memory.txt` by calling ``NumCalc -estimate_ram``. This must
    be done before calling this function.

    Parameters
    ----------
    folder : str
        full path to the source folder containing the `Memory.txt` file from
        which the estimates are read

    Returns
    -------
    estimates : numpy array
        An array of shape ``(N, 3)`` where ``N`` is the number of frequency
        steps. The first column contains the frequency step, the second the
        frequency in Hz, and the third the estimated RAM consumption in GB.
    """

    # check if file exists
    if not os.path.isfile(os.path.join(folder, "Memory.txt")):
        raise ValueError(f"{folder} does not contain a Memory.txt file")

    # read content of file
    with open(os.path.join(folder, "Memory.txt"), "r") as ff:
        content = ff.readlines()

    # parse data to nested list
    estimates = []
    for line in content:
        estimate = []
        for ee in line.strip().split(" "):
            estimate.append(float(ee))

        estimates.append(estimate)

    return np.asarray(estimates)


def write_output_report(folder=None):
    r"""
    Generate project report from NumCalc output files.

    NumCalc (Mesh2HRTF's numerical core) writes information about the
    simulations to the files `NC*.out` located under `NumCalc/source_*`. The
    file `NC.out` exists if NumCalc was ran without the additional command line
    parameters ``-istart`` and ``-iend``. If these parameters were used, there
    is at least one `NC\*-\*.out`. If this is the case, information from
    `NC\*-\*.out` overwrites information from NC.out in the project report.

    .. note::

        The project reports are written to the files
        `Output2HRTF/report_source_*.csv`. If issues were detected, they are
        listed in `Output2HRTF/report_issues.csv`.

    The report contain the following information

    Frequency step
        The index of the frequency.
    Frequency in Hz
        The frequency in Hz.
    NC input
        Name of the input file from which the information was taken.
    Input check passed
        Contains a 1 if the check of the input data passed and a 0 otherwise.
        If the check failed for one frequency, the following frequencies might
        be affected as well.
    Converged
        Contains a 1 if the simulation converged and a 0 otherwise. If the
        simulation did not converge, the relative error might be high.
    Num. iterations
        The number of iterations that were required to converge
    relative error
        The relative error of the final simulation
    Comp. time total
        The total computation time in seconds
    Comp. time assembling
        The computation time for assembling the matrices in seconds
    Comp. time solving
        The computation time for solving the matrices in seconds
    Comp. time post-proc
        The computation time for post-processing the results in seconds


    Parameters
    ----------
    folder : str, optional
        The path of the Mesh2HRTF project folder, i.e., the folder containing
        the subfolders EvaluationsGrids, NumCalc, and ObjectMeshes. The
        default, ``None`` uses the current working directory.

    Returns
    -------
    found_issues : bool
        ``True`` if issues were found, ``False`` otherwise
    report : str
        The report or an empty string if no issues were found
    """

    # get sources and number of sources and frequencies
    sources = glob.glob(os.path.join(folder, "NumCalc", "source_*"))
    num_sources = len(sources)

    if os.path.exists(os.path.join(folder, '..', "parameters.json")):
        with open(os.path.join(folder, '..', "parameters.json"), "r") as file:
            params = json.load(file)
    else:
        with open(os.path.join(folder, "parameters.json"), "r") as file:
            params = json.load(file)

    # sort source files (not read in correct order in some cases)
    nums = [int(source.split("_")[-1]) for source in sources]
    sources = np.array(sources)
    sources = sources[np.argsort(nums)]

    # parse all NC*.out files for all sources
    all_files, fundamentals, out, out_names = m2s.utils._parse_nc_out_files(
        sources, num_sources, params["numFrequencies"])

    # write report as csv file
    m2s.utils._write_project_reports(folder, all_files, out, out_names)

    # look for errors
    report = m2s.utils._check_project_report(folder, fundamentals, out)

    found_issues = True if report else False

    return found_issues, report
