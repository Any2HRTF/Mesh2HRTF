import os
import warnings
import json
import numpy as np
import pyfar as pf
import glob
import sofar as sf
from mesh2scattering import utils
import csv
import re


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
    freq_new = freq[..., data_mask, :]
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
        data_in = data[np.abs(thetas-theta) < 1e-3, :]
        freq[ii, ...] = shift_data_coords(
            data_in, coords_mic, float(az)).freq.copy()
    data_out = pf.FrequencyData(freq, data.frequencies)
    return data_out


def apply_symmetry_mirror(
        data_in, coords_mic, incident_coords, mirror_axe=None):
    all_az = np.sort(np.array(list(set(
        np.round(incident_coords.get_sph(unit='deg')[..., 0], 5)))))
    all_el = np.sort(np.array(list(set(
        np.round(incident_coords.get_sph(unit='deg')[..., 1], 5)))))
    radius = np.median(incident_coords.get_sph()[..., 2])
    source_coords_out = pf.Coordinates(
        *np.meshgrid(all_az, all_el, indexing='ij'), radius, 'sph', unit='deg')
    data = reshape_to_az_by_el(data_in, incident_coords, source_coords_out, 1)

    shape = list(data.cshape)
    shape.append(data.n_bins)
    azimuths = np.sort(np.array(list(set(
        np.round(incident_coords.get_sph()[..., 0], 5)))))
    index_min = len(azimuths)-1
    index_max = 2 * index_min + 1
    shape[mirror_axe] = index_max

    freq = np.empty(shape, dtype=complex)
    freq[:] = np.nan
    freq = np.moveaxis(freq, mirror_axe, -1)
    freq_in = np.moveaxis(data.freq.copy(), mirror_axe, -1)
    max_aimuth = np.max(azimuths)
    radius = np.median(incident_coords.get_sph()[:, 2])
    azimuths_new = []
    for iaz in range(index_max):
        if iaz > index_min:
            idx = index_max-iaz-1
            az = (max_aimuth-azimuths[idx]) * 2
            if azimuths[idx] + az > 2 * np.pi:
                az = 2 * np.pi - azimuths[idx]
            data_swap = pf.FrequencyData(
                np.moveaxis(data.freq, 0, -2), data.frequencies)
            data_in = shift_data_coords(
                data_swap, coords_mic, az/np.pi*180).freq
            data_in = np.moveaxis(data_in, -2, 0)

            azimuths_new.append(azimuths[idx] + az)
        else:
            data_in = data.freq
            idx = iaz
            azimuths_new.append(azimuths[iaz])
        freq_in = np.moveaxis(data_in, mirror_axe, 0)
        freq[..., iaz] = freq_in[idx, ...]
    # if max_index > 0:
    #     freq = freq[..., :max_index]
    freq = np.moveaxis(freq, -1, mirror_axe)
    data_out = pf.FrequencyData(freq, data.frequencies)
    elevations = np.sort(np.array(list(set(
        np.round(incident_coords.get_sph()[..., 1], 5)))))
    new_inc_coords = angles2coords(np.array(azimuths_new), elevations, radius)
    shape = data_out.cshape
    data_out.freq = np.reshape(
        data_out.freq, (shape[0], shape[1]*shape[2], data_out.n_bins))
    xyz = new_inc_coords.get_cart().reshape((new_inc_coords.csize, 3))
    new_inc_coords = pf.Coordinates(xyz[..., 0], xyz[..., 1], xyz[..., 2])
    mask_not_double = new_inc_coords.get_sph()[..., 1] != 0
    mask_not_double[np.argmax(~mask_not_double)] = True
    return data_out[:, mask_not_double], new_inc_coords[mask_not_double]


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
        source_position = np.array(params["sources"])
        if source_position.shape[1] != 3:
            source_position = np.transpose(source_position)
        receiver_position = np.array(evaluationGrids[grid]["nodes"][:, 1:4])
        if receiver_position.shape[1] != 3:
            receiver_position = np.transpose(receiver_position)

        # apply symmetry of reference sample
        data = pf.FrequencyData(
            evaluationGrids[grid]["pressure"], params["frequencies"])
        receiver_coords = _cart_coordinates(receiver_position)
        source_coords = _cart_coordinates(source_position)
        # data = np.swapaxes(data, 0, 1)
        data_out, source_coords = apply_symmetry_mirror(
            data, receiver_coords, source_coords, 1)
        data_out, source_coords = apply_symmetry_mirror(
            data_out, receiver_coords, source_coords, 1)

        # write data
        sofa = utils._get_sofa_object(
            data_out.freq,
            source_coords.get_cart(),
            receiver_position,
            params["mesh2scattering_version"],
            frequencies=params["frequencies"])

        sofa.GLOBAL_Title = folder.split(os.sep)[-1]

        # write scattered sound pressure to SOFA file
        sf.write_sofa(os.path.join(
            folder, 'sample.pattern.sofa'), sofa)

    evaluationGrids, params = read_numcalc(
        os.path.join(folder, 'reference'), True)

    # process BEM data for writing scattered sound pressure to SOFA files
    for grid in evaluationGrids:
        print(f'\nWrite sample data "{grid}" ...\n')
        # get pressure as SOFA object (all following steps are run on SOFA
        # objects. This way they are available to other users as well)
        # read source and receiver positions
        source_position_ref = np.array(params["sources"])
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
            source_coords)

        # create sofa file
        sofa = utils._get_sofa_object(
            data_out.freq,
            source_coords.get_cart(),
            receiver_position_ref,
            params["mesh2scattering_version"],
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
    all_files, fundamentals, out, out_names = _parse_nc_out_files(
        sources, num_sources, params["num_frequencies"])

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
    # sources_num, sourceType, sources, sourceArea,
    # num_frequencies, frequencies
    params = os.path.join(folder, '..', 'parameters.json')
    if not os.path.isfile(params):
        raise ValueError((
            f"The folder {folder} is not a valid Mesh2scattering project. "
            "It must contain the file 'parameters.json'"))

    with open(params, "r") as file:
        params = json.load(file)

    # get source positions
    source_coords = np.transpose(np.array(params['sources']))
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
    evaluationGrids, _ = _read_nodes_and_elements(
        os.path.join(folder, 'EvaluationGrids'))

    # Load EvaluationGrid data
    if is_ref:
        xyz = np.array(params["sources"])
        coords = pf.Coordinates(xyz[..., 0], xyz[..., 1], xyz[..., 2])
        num_sources = np.sum(np.abs(coords.get_sph()[..., 0]) < 1e-12)
    else:
        num_sources = params["sources_num"]

    if not len(evaluationGrids) == 0:
        pressure, _ = _read_numcalc_data(
            num_sources, params["num_frequencies"],
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
    all_files, fundamentals, out, out_names = _parse_nc_out_files(
        sources, num_sources, params["num_frequencies"])

    # write report as csv file
    _write_project_reports(folder, all_files, out, out_names)

    # look for errors
    report = _check_project_report(folder, fundamentals, out)

    found_issues = True if report else False

    return found_issues, report


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
            delimiter=' ', skiprows=1, dtype=np.float64)

        tmpElements = np.loadtxt(os.path.join(
            folder, grid, 'Elements.txt'),
            delimiter=' ', skiprows=1, dtype=np.float64)

        grids[grid] = {
            "nodes": tmpNodes,
            "elements": tmpElements,
            "num_nodes": tmpNodes.shape[0],
            "num_elements": tmpElements.shape[0]}

        gridsNumNodes += grids[grid]['num_nodes']

    return grids, gridsNumNodes


def _read_numcalc_data(sources_num, num_frequencies, folder, data):
    """Read the sound pressure on the object meshes or evaluation grid."""
    pressure = []

    if data not in ['pBoundary', 'pEvalGrid', 'vBoundary', 'vEvalGrid']:
        raise ValueError(
            'data must be pBoundary, pEvalGrid, vBoundary, or vEvalGrid')

    for source in range(sources_num):

        tmpFilename = os.path.join(
            folder, 'NumCalc', f'source_{source+1}', 'be.out')
        tmpPressure, indices = _load_results(
            tmpFilename, data, num_frequencies)

        pressure.append(tmpPressure)

    pressure = np.transpose(np.array(pressure), (2, 0, 1))

    return pressure, indices


def _load_results(foldername, filename, num_frequencies):
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
    num_frequencies : int
        the number of simulated frequencies

    Returns
    -------
    data : numpy array
        Pressure or abs velocity values of shape (num_frequencies, numEntries)
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
    data = np.zeros((num_frequencies, numDatalines), dtype=dtype)

    for ii in range(num_frequencies):
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


def _check_project_report(folder, fundamentals, out):

    # return if there are no fundamental errors or other issues
    if not all([all(f) for f in fundamentals]) and not np.any(out == -1) \
            and np.all(out[:, 3:5]):
        return

    # report detailed errors
    report = ""

    for ss in range(out.shape[2]):

        # currently we detect frequencies that were not calculated and
        # frequencies with convergence issues
        missing = "Frequency steps that were not calculated:\n"
        input_test = "Frequency steps with bad input:\n"
        convergence = "Frequency steps that did not converge:\n"

        any_missing = False
        any_input_failed = False
        any_convergence = False

        # loop frequencies
        for ff in range(out.shape[0]):

            f = out[ff, :, ss]

            # no value for frequency
            if f[1] == -1:
                any_missing = True
                missing += f"{int(f[0])}, "
                continue

            # input data failed
            if f[3] == 0:
                any_input_failed = True
                input_test += f"{int(f[0])}, "

            # convergence value is zero
            if f[4] == 0:
                any_convergence = True
                convergence += f"{int(f[0])}, "

        if any_missing or any_input_failed or any_convergence:
            report += f"Detected issues for source {ss+1}\n"
            report += "----------------------------\n"
            if any_missing:
                report += missing[:-2] + "\n\n"
            if any_input_failed:
                report += input_test[:-2] + "\n\n"
            if any_convergence:
                report += convergence[:-2] + "\n\n"

    if not report:
        report = ("\nDetected unknown issues\n"
                  "-----------------------\n"
                  "Check the project reports in Output2HRTF,\n"
                  "and the NC*.out files in NumCalc/source_*\n\n")

    report += ("For more information check Output2HRTF/report_source_*.csv "
               "and the NC*.out files located at NumCalc/source_*")

    # write to disk
    report_name = os.path.join(
        folder, "Output2HRTF", "report_issues.txt")
    with open(report_name, "w") as f_id:
        f_id.write(report)

    return report


def _write_project_reports(folder, all_files, out, out_names):
    """
    Write project report to disk at folder/Output2HRTF/report_source_*.csv

    For description of input parameter refer to write_output_report and
    _parse_nc_out_files
    """

    # loop sources
    for ss in range(out.shape[2]):

        report = ", ".join(out_names) + "\n"

        # loop frequencies
        for ff in range(out.shape[0]):
            f = out[ff, :, ss]
            report += (
                f"{int(f[0])}, "                # frequency step
                f"{float(f[1])}, "              # frequency in Hz
                f"{all_files[ss][int(f[2])]},"  # NC*.out file
                f"{int(f[3])}, "                # input check
                f"{int(f[4])}, "                # convergence
                f"{int(f[5])}, "                # number of iterations
                f"{float(f[6])}, "              # relative error
                f"{int(f[7])}, "                # total computation time
                f"{int(f[8])}, "                # assembling equations time
                f"{int(f[9])}, "                # solving equations time
                f"{int(f[10])}\n"               # post-processing time
                )

        # write to disk
        report_name = os.path.join(
            folder, "Output2HRTF", f"report_source_{ss + 1}.csv")
        with open(report_name, "w") as f_id:
            f_id.write(report)


def _parse_nc_out_files(sources, num_sources, num_frequencies):
    """
    Parse all NC*.out files for all sources.

    This function should never raise a value error, regardless of how mess
    NC*.out files are. Looking for error is done at a later step.

    Parameters
    ----------
    sources : list of strings
        full path to the source folders
    num_sources : int
        number of sources - len(num_sources)
    num_frequencies : int
        number of frequency steps

    Returns
    -------
    out : numpy array
        containing the extracted information for each frequency step
    out_names : list of string
        verbal information about the columns of `out`
    """

    # array for reporting fundamental errors
    fundamentals = []
    all_files = []

    # array for saving the detailed report
    out_names = ["frequency step",         # 0
                 "frequency in Hz",        # 1
                 "NC input file",          # 2
                 "Input check passed",     # 3
                 "Converged",              # 4
                 "Num. iterations",        # 5
                 "relative error",         # 6
                 "Comp. time total",       # 7
                 "Comp. time assembling",  # 8
                 "Comp. time solving",     # 9
                 "Comp. time post-proc."]  # 10
    out = -np.ones((num_frequencies, 11, num_sources))
    # values for steps
    out[:, 0] = np.arange(1, num_frequencies + 1)[..., np.newaxis]
    # values indicating failed input check and non-convergence
    out[:, 3] = 0
    out[:, 4] = 0

    # regular expression for finding a number that can be int or float
    re_number = r"(\d+(?:\.\d+)?)"

    # loop sources
    for ss, source in enumerate(sources):

        # list of NC*.out files for parsing
        files = glob.glob(os.path.join(source, "NC*.out"))

        # make sure that NC.out is first
        nc_out = os.path.join(source, "NC.out")
        if nc_out in files and files.index(nc_out):
            files = [files.pop(files.index(nc_out))] + files

        # update fundamentals
        fundamentals.append([0 for f in range(len(files))])
        all_files.append([os.path.basename(f) for f in files])

        # get content from all NC*.out
        for ff, file in enumerate(files):

            # read the file and join all lines
            with open(file, "r") as f_id:
                lines = f_id.readlines()
            lines = "".join(lines)

            # split header and steps
            lines = lines.split(
                ">> S T E P   N U M B E R   A N D   F R E Q U E N C Y <<")

            # look for fundamental errors
            if len(lines) == 1:
                fundamentals[ss][ff] = 1
                continue

            # parse frequencies (skip header)
            for line in lines[1:]:

                # find frequency step
                idx = re.search(r'Step \d+,', line)
                if idx:
                    step = int(line[idx.start()+5:idx.end()-1])

                # write number of input file (replaced by string later)
                out[step-1, 2, ss] = ff

                # find frequency
                idx = re.search(f'Frequency = {re_number} Hz', line)
                if idx:
                    out[step-1, 1, ss] = float(
                        line[idx.start()+12:idx.end()-3])

                # check if the input data was ok
                if "Too many integral points in the theta" not in line:
                    out[step-1, 3, ss] = 1

                # check and write convergence
                if 'Maximum number of iterations is reached!' not in line:
                    out[step-1, 4, ss] = 1

                # check iterations
                idx = re.search(r'number of iterations = \d+,', line)
                if idx:
                    out[step-1, 5, ss] = int(line[idx.start()+23:idx.end()-1])

                # check relative error
                idx = re.search('relative error = .+', line)
                if idx:
                    out[step-1, 6, ss] = float(line[idx.start()+17:idx.end()])

                # check time stats
                # -- assembling
                idx = re.search(
                    r'Assembling the equation system         : \d+',
                    line)
                if idx:
                    out[step-1, 8, ss] = float(line[idx.start()+41:idx.end()])

                # -- solving
                idx = re.search(
                    r'Solving the equation system            : \d+',
                    line)
                if idx:
                    out[step-1, 9, ss] = float(line[idx.start()+41:idx.end()])

                # -- post-pro
                idx = re.search(
                    r'Post processing                        : \d+',
                    line)
                if idx:
                    out[step-1, 10, ss] = float(line[idx.start()+41:idx.end()])

                # -- total
                idx = re.search(
                    r'Total                                  : \d+',
                    line)
                if idx:
                    out[step-1, 7, ss] = float(line[idx.start()+41:idx.end()])

    return all_files, fundamentals, out, out_names
