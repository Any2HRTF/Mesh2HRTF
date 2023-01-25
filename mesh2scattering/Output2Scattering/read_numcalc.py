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
from . import _utils


def read_numcalc(folder=None):
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

    # get source positions
    source_coords = np.transpose(np.array(params['sourceCenter']))
    source_coords = pf.Coordinates(
        source_coords[..., 0], source_coords[..., 1], source_coords[..., 2])

    # output directory
    if not os.path.exists(os.path.join(folder, 'Output2HRTF')):
        os.makedirs(os.path.join(folder, 'Output2HRTF'))

    # write the project report and check for issues
    print('\n Writing the project report ...')
    found_issues, report = m2s.write_output_report(folder)

    if found_issues:
        warnings.warn(report)

    # get the evaluation grids
    evaluationGrids, _ = _utils._read_nodes_and_elements(
        os.path.join(folder, 'EvaluationGrids'))

    # Load EvaluationGrid data
    if not len(evaluationGrids) == 0:
        pressure, _ = _utils._read_numcalc_data(
            params["numSources"], params["numFrequencies"],
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
