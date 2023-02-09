"""
Python tools for Mesh2HRTF including functions to generate SOFA files
containing the HRTF/HRIR data, merge SOFA files containing data for the left
and right ear and generate evaluation grids.
"""
import os
import warnings
import numpy as np
import json
import mesh2scattering as m2s
import pyfar as pf


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
    found_issues, report = m2s.NumCalc.write_output_report(folder)

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
