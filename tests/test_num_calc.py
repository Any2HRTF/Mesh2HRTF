import pytest
import subprocess
import tempfile
import shutil
import os
import scipy.io
import numpy
import utils

create_baseline = False

def test_build():
    """ test if make for NumCalc works """

    # create a temporary directory
    tmp = tempfile.TemporaryDirectory(dir=os.getcwd())

    shutil.copytree(os.path.join(os.path.dirname(__file__),
                                 "..", "mesh2hrtf", "NumCalc", "Source"),
                    tmp.name+"/NumCalc")
    tmp_path = os.path.join(tmp.name, "NumCalc")
    subprocess.run(["make"], cwd=tmp_path, check=True)


@pytest.mark.parametrize("boundary_condition", [("rigid"), ("soft")])
@pytest.mark.parametrize("source,range_a", [("plane", (10, -20)),
                                            ("point", (40, -45))])
@pytest.mark.parametrize("bem_method", [("ml-fmm-bem"), ("fmm-bem"), ("bem")])
def test_against_reference(boundary_condition, source, bem_method,
                           range_a, range_b=(-1, 1)):
    """
    test if NumCalc and Output2HRTF.py generate correct output by comparing to
    analytical solution
    """
    # Setup

    # create temporary directory
    tmp = tempfile.TemporaryDirectory(dir=os.getcwd())

    # copy test directory
    shutil.copytree(os.path.join(os.path.dirname(__file__),
                                 'test_numcalc_project'),
                    os.path.join(tmp.name, 'project'))
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
    hrtf_sim = utils.hrtf_sofa_to_numpy(
        os.path.join(tmp_path, "Output2HRTF", "HRTF_HorPlane.sofa"))
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
    utils.scatter_reference_vs_analytic(
        hrtf_sim, hrtf_ana, xyz[:, 0], xyz[:, 1],
        range_a, range_b, boundary_condition, source, bem_method)


@pytest.mark.parametrize("boundary_condition", [("rigid")])
@pytest.mark.parametrize("source,range_a", [("rightear", (40, -40))])
@pytest.mark.parametrize("bem_method", [("ml-fmm-bem")])
def test_ears(boundary_condition, source, bem_method,
                        range_a, range_b=(-1, 1)):
    """
    test if NumCalc and Output2HRTF.py generate correct output by comparing to
    analytical solution
    """
    # Setup

    # create temporary directory
    tmp = tempfile.TemporaryDirectory(dir=os.getcwd())

    # copy test directory
    shutil.copytree(os.path.join(os.path.dirname(__file__),
                    'test_numcalc_ear_projects',source),
                    os.path.join(tmp.name, 'project'))
    # shutil.copytree(os.path.join(os.path.dirname(__file__),
    #                             'test_numcalc_ear_projects','refproj'),
    #                 os.path.join(tmp.name, 'project'))
    
    # # copy correct input file and rename it to NC.inp
    # if source == "bothears":
    #     shutil.copyfile(os.path.join(os.path.dirname(__file__),
    #                     'test_numcalc_input_files', 'NC_'+
    #                                 'leftear'+'.inp'),
    #                     os.path.join(tmp.name, 'project', 'NumCalc',
    #                                 'source_1', 'NC.inp'))
    #     shutil.copyfile(os.path.join(os.path.dirname(__file__),
    #                     'test_numcalc_input_files', 'NC_'+
    #                                 'rightear'+'.inp'),
    #                     os.path.join(tmp.name, 'project', 'NumCalc',
    #                                 'source_2', 'NC.inp'))
    # else:
    #     shutil.copyfile(os.path.join(os.path.dirname(__file__),
    #                     'test_numcalc_input_files', 'NC_'+
    #                                 source+'.inp'),
    #                     os.path.join(tmp.name, 'project', 'NumCalc',
    #                                 'source_1', 'NC.inp'))

    # Exercise

    # run NumCalc with subprocess
    tmp_path = os.path.join(tmp.name, "project", "NumCalc", "source_1")
    subprocess.run(["NumCalc"], cwd=tmp_path, check=True)
    if source == "bothears":
        tmp_path = os.path.join(tmp.name, "project", "NumCalc", "source_2")
        subprocess.run(["NumCalc"], cwd=tmp_path, check=True)
    # run Output2HRTF.py
    tmp_path = os.path.join(tmp.name, "project")
    subprocess.run(["python", "Output2HRTF.py"], cwd=tmp_path, check=True)

    # Verify

    # load HRTF data from simulation as numpy
    hrtf_sim = utils.hrtf_sofa_to_numpy(
        os.path.join(tmp_path, "Output2HRTF", "HRTF_HorPlane.sofa"))
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
    utils.scatter_reference_vs_analytic(
        hrtf_sim, hrtf_ana, xyz[:, 0], xyz[:, 1],
        range_a, range_b, boundary_condition, source, bem_method)
