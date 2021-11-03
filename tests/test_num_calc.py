import pytest
import subprocess
import tempfile
import shutil
import os
import hrtf_sofa_to_numpy as hstn
import scipy.io
import numpy
import tests_utils as tu


# def test_build():
#     """ test if make for NumCalc works """

#     create a temporary directory
#     tmp = tempfile.TemporaryDirectory(dir=os.getcwd())

#     shutil.copytree("../mesh2hrtf/NumCalc/Source", tmp.name+"/NumCalc")
#     tmp_path = os.path.join(tmp.name, "NumCalc")
#     subprocess.run(["make"], cwd=tmp_path, check=True)

@pytest.mark.parametrize("boundary_condition", [("rigid"), ("soft")])
@pytest.mark.parametrize("source", [("plane"), ("point")])
@pytest.mark.parametrize("bem_method", [("bem"), ("fmm-bem"), ("ml-fmm-bem")])
def test_numcalc(boundary_condition, source, bem_method):
    """
    test if NumCalc and Output2HRTF.py generate correct output by comparing to
    analytical solution
    """
    # Setup

    # create temporary directory
    tmp = tempfile.TemporaryDirectory(dir=os.getcwd())

    # copy test directory
    shutil.copytree(os.path.join(os.path.dirname(__file__),
                    'test_numcalc_project'), os.path.join(tmp.name, 'project'))
    # copy correct input file and rename it to NC.inp
    shutil.copyfile(os.path.join(os.path.dirname(__file__),
                    'test_numcalc_input_files', 'NC_'+boundary_condition+'_' +
                                 source+'_'+bem_method+'.inp'),
                    os.path.join(tmp.name, 'project', 'NumCalc',
                                 'CPU_1_Core_1', 'NC.inp'))

    # Exercise

    # run NumCalc with subprocess
    tmp_path = os.path.join(tmp.name, "project", "NumCalc", "CPU_1_Core_1")
    subprocess.run(["NumCalc"], cwd=tmp_path, check=True)
    # run Output2HRTF.py
    tmp_path = os.path.join(tmp.name, "project")
    subprocess.run(["python", "Output2HRTF.py"], cwd=tmp_path, check=True)

    # Verify

    # load HRTF data from simulation as numpy
    hrtf_sim = hstn.hrtf_sofa_to_numpy(os.path.join(tmp_path, "Output2HRTF",
                                                    "HRTF_HorPlane.sofa"))
    # normalize because only relative differences of interest
    hrtf_sim = hrtf_sim[:, :, 0]/numpy.mean(
                numpy.abs(hrtf_sim[numpy.isfinite(hrtf_sim)]))

    # load HRTF data from analytical comparison
    ana_path = os.path.join(os.path.dirname(__file__),
                            'test_numcalc_analytical_references',
                            'ref_'+boundary_condition+'_'+source+'.mat')
    mat_ana = scipy.io.loadmat(ana_path)
    hrtf_ana = mat_ana['p_total']
    # normalize because only relative differences of interest
    hrtf_ana = hrtf_ana/numpy.mean(
                numpy.abs(hrtf_ana[numpy.isfinite(hrtf_ana)]))

    # compare
    numpy.testing.assert_allclose(
        numpy.abs(hrtf_sim[numpy.isfinite(hrtf_sim)]),
        numpy.abs(hrtf_ana[numpy.isfinite(hrtf_ana)]), rtol=11.1)

    xyz = mat_ana["XYZ"]
    tu.scatter_reference_vs_analytic(hrtf_sim, hrtf_ana, xyz[:, 0], xyz[:, 1],
                                     (0, -10), (5, -5),
                                     boundary_condition, source, bem_method)


# test_bem("soft", "plane")
