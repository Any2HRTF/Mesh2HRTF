"""
Python tools for Mesh2HRTF including functions to generate SOFA files
containing the HRTF/HRIR data, merge SOFA files containing data for the left
and right ear and generate evaluation grids.
"""
import os
import numpy as np
import sofar as sf
import mesh2scattering as m2s
import pyfar as pf
from . import _utils


def output2scattering(folder, strutural_wavelength):
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
    evaluationGrids, params = m2s.read_numcalc(os.path.join(folder, 'sample'))

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
        sofa = _utils._get_sofa_object(
            evaluationGrids[grid]["pressure"],
            source_position,
            receiver_position,
            params["Mesh2HRTF_Version"],
            frequencies=params["frequencies"])

        sofa.GLOBAL_Title = folder.split(os.sep)[-1]
        sofa.GLOBAL_References = f'{strutural_wavelength}'

        # write HRTF data to SOFA file
        sf.write_sofa(os.path.join(
            folder, f'sample_{grid}.pattern.sofa'), sofa)

    evaluationGrids, params = m2s.read_numcalc(
        os.path.join(folder, 'reference'))

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
        data_out = m2s.apply_symmetry_circular(
            pf.FrequencyData(data, params["frequencies"]),
            _cart_coordiantes(receiver_position_ref),
            _cart_coordiantes(source_position_ref),
            _cart_coordiantes(source_position))

        # create sofa file
        sofa = _utils._get_sofa_object(
            data_out.freq,
            source_position,
            receiver_position_ref,
            params["Mesh2HRTF_Version"],
            frequencies=params["frequencies"])

        sofa.GLOBAL_Title = folder.split(os.sep)[-1]

        # write HRTF data to SOFA file
        sf.write_sofa(os.path.join(
            folder, f'reference_{grid}.pattern.sofa'), sofa)

    print('Done\n')


def _cart_coordiantes(xyz):
    return pf.Coordinates(xyz[:, 0], xyz[:, 1], xyz[:, 2])
