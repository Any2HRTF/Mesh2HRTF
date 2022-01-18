import pytest
import subprocess
import tempfile
import shutil
import os
import hrtf_sofa_to_numpy as hstn
import scipy.io
import numpy
import test_utils as tu


def test_blender_export():
    # create a temporary directory
    tmp = tempfile.TemporaryDirectory(dir=os.getcwd())

    # copy test directory
    shutil.copytree(os.path.join(os.path.dirname(__file__),
                    'script_blender_export_test'), os.path.join(tmp.name, 'project'))
    blender_path = os.path.join('/home', 'matheson', 'Apps', 'blender-2.91.0', 'blender')
    tmp_path = os.path.join(tmp.name, 'project')
    # subprocess.run("pwd", cwd=tmp_path, check=True)
    # subprocess.run("ls", cwd=tmp_path, check=True)
    # subprocess.run([blender_path, "3dModel.blend",
    #                 "--background", "--python",
    #                 "blender_script.py"], cwd=tmp_path, check=True)

# subprocess.run(["/home/matheson/Apps/blender-2.91.0/blender 3dModel.blend --background --python blender_script.py"], cwd=tmp_path)


def test_build():
    """ test if make for NumCalc works """

    # create a temporary directory
    tmp = tempfile.TemporaryDirectory(dir=os.getcwd())

    shutil.copytree("../mesh2hrtf/NumCalc/Source", tmp.name+"/NumCalc")
    tmp_path = os.path.join(tmp.name, "NumCalc")
    subprocess.run(["make"], cwd=tmp_path, check=True)


@pytest.mark.parametrize("boundary_condition", [("rigid"), ("soft")])
@pytest.mark.parametrize("source,range_a", [("plane", (10, -20)),
                                            ("point", (40, -45))])
@pytest.mark.parametrize("bem_method", [("ml-fmm-bem"), ("fmm-bem"), ("bem")])
def test_numcalc(boundary_condition, source, bem_method, range_a, range_b=(-1, 1)):
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
                                 'source_1', 'NC.inp'))

    # Exercise

    # run NumCalc with subprocess
    tmp_path = os.path.join(tmp.name, "project", "NumCalc", "source_1")
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
                                     range_a, range_b,
                                     boundary_condition, source, bem_method)


def test_two_sources():
    """
    test Output2HRTF.py handling two sound sources ("Both ears" condition)
    """
    # Setup

    # create temporary directory
    tmp = tempfile.TemporaryDirectory(dir=os.getcwd())

    # copy test directory
    shutil.copytree(os.path.join(os.path.dirname(__file__),
                    'test_2sources_project'), os.path.join(tmp.name, 'project'))

    # Exercise

    # # run NumCalc with subprocess
    # for iSource in (1, 2):
    #     tmp_path = os.path.join(tmp.name, "project", "NumCalc", f"source_{iSource+1}")
    #     subprocess.run(["NumCalc"], cwd=tmp_path, check=True)
    
    # run Output2HRTF.py
    tmp_path = os.path.join(tmp.name, "project")
    subprocess.run(["python", "Output2HRTF.py"], cwd=tmp_path, check=True)

# test_blender_export()
# test_numcalc("rigid", "point", "ml-fmm-bem", (10, -20), range_b=(-1, 1))
