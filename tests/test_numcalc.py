# %%
import pytest
import subprocess
import tempfile
import shutil
import os
import scipy.io
import numpy
import utils
import mesh2scattering as m2s

# directory of this file
base_dir = os.path.dirname(__file__)

# Build NumCalc locally to use for testing
tmp = tempfile.TemporaryDirectory()
numcalc = os.path.join(tmp.name, "NumCalc", "bin", "NumCalc")

shutil.copytree(
    os.path.join(base_dir, "..", "mesh2scattering", "NumCalc"),
    os.path.join(tmp.name, "NumCalc"))

if os.path.isfile(numcalc):
    os.remove(numcalc)

subprocess.run(
    ["make"], cwd=os.path.join(tmp.name, "NumCalc", "src"), check=True)


def test_numcalc_invalid_parameter(capfd):
    """
    Test if NumCalc throws an error in case of invalid command line parameter.
    """

    try:
        subprocess.run([f'{numcalc} -invalid_parameter'],
                       check=True, shell=True)
    except subprocess.CalledProcessError:
        _, err = capfd.readouterr()
        assert "NumCalc was called with an unknown parameter or flag." in err
    else:
        ValueError("Num calc did not throw an error")


@pytest.mark.parametrize("nitermax, use", [(0, True), (1, True), (2, True),
                                           ([], False)])
def test_numcalc_commandline_nitermax(nitermax, use):
    """Test if command line parameter nitermax behaves as expected"""
    # Setup

    # create temporary directory
    tmp = tempfile.TemporaryDirectory()

    # copy test directory
    shutil.copytree(
        os.path.join(base_dir, 'resources', 'test_numcalc',
                     'project_folder_pspw'),
        os.path.join(tmp.name, 'project'))
    # copy correct input file and rename it to NC.inp
    os.mkdir(os.path.join(tmp.name, 'project', 'NumCalc'))
    os.mkdir(os.path.join(tmp.name, 'project', 'NumCalc', 'source_1'))
    shutil.copyfile(
        os.path.join(base_dir, 'resources', 'test_numcalc',
                     'ncinp_files', 'NC_commandline_parameters.inp'),
        os.path.join(tmp.name, 'project', 'NumCalc', 'source_1', 'NC.inp'))

    if use:
        commandLineArgument = f' -nitermax {nitermax}'
    else:
        commandLineArgument = ''

    # Exercise

    # run NumCalc with subprocess
    tmp_path = os.path.join(tmp.name, "project", "NumCalc", "source_1")
    subprocess.run([f'{numcalc}{commandLineArgument}'],
                   cwd=tmp_path, check=True, shell=True)

    # Verify
    out_filename = 'NC.out'
    out_filepath = os.path.join(tmp.name, "project", "NumCalc",
                                "source_1", out_filename)

    out_file = open(out_filepath)
    out_text = out_file.read()

    if use:
        assert f'CGS solver: number of iterations = {nitermax}' in out_text
        assert 'Warning: Maximum number of iterations is reached!' in out_text
    else:
        assert 'Warning: Maximum number of iterations is reached!' \
            not in out_text


@pytest.mark.parametrize("istart, iend", [(False, False), (3, False),
                                          (False, 3), (2, 3)])
def test_numcalc_commandline_istart_iend(istart, iend):
    """Test if command line parameters istart and iend behave as expected"""
    # Setup

    # create temporary directory
    tmp = tempfile.TemporaryDirectory()

    # copy test directory
    shutil.copytree(
        os.path.join(base_dir, 'resources', 'test_numcalc',
                     'project_folder_pspw'),
        os.path.join(tmp.name, 'project'))
    # copy correct input file and rename it to NC.inp
    os.mkdir(os.path.join(tmp.name, 'project', 'NumCalc'))
    os.mkdir(os.path.join(tmp.name, 'project', 'NumCalc', 'source_1'))
    shutil.copyfile(
        os.path.join(base_dir, 'resources', 'test_numcalc',
                     'ncinp_files', 'NC_commandline_parameters.inp'),
        os.path.join(tmp.name, 'project', 'NumCalc', 'source_1', 'NC.inp'))

    commandLineArgument = ''
    if istart > 0:
        commandLineArgument += f' -istart {istart}'
    if iend > 0:
        commandLineArgument += f' -iend {iend}'

    # Exercise

    # run NumCalc with subprocess
    tmp_path = os.path.join(tmp.name, "project", "NumCalc", "source_1")
    subprocess.run([f'{numcalc}{commandLineArgument}'],
                   cwd=tmp_path, check=True, shell=True)

    # Verify
    if (not istart and not iend):
        out_filename = 'NC.out'
    elif ((istart > 0) and not iend):
        out_filename = f'NCfrom{istart}.out'
    elif (not istart and (iend > 0)):
        out_filename = f'NCuntil{iend}.out'
    elif ((istart > 0) and (iend > 0)):
        out_filename = f'NC{istart}-{iend}.out'
    else:
        raise Exception("Wrong istart and/or iend parameters chosen")

    out_filepath = os.path.join(tmp.name, "project", "NumCalc",
                                "source_1", out_filename)

    with open(out_filepath) as out_file:
        out_text = out_file.read()

    if istart > 0:
        assert f'Step {istart-1}' not in out_text
        assert f'Step {istart}' in out_text
    else:
        assert 'Step 1' in out_text

    if iend > 0:
        assert f'Step {iend}' in out_text
        assert f'Step {iend+1}' not in out_text

    if istart > 0 and iend > 0:
        nStepsActual = out_text.count((
            '>> S T E P   N U M B E R   A N D   F R E Q U E N C Y <<'))
        nStepsExpected = iend - istart + 1
        assert nStepsActual == nStepsExpected


def test_numcalc_commandline_estimate_ram():
    """Test NumCalc's RAM estimation using -estimate_ram"""

    # copy test data
    cwd = tempfile.TemporaryDirectory()
    data_cwd = os.path.join(cwd.name, 'SHTF', 'NumCalc', 'source_1')
    data_shtf = os.path.join(os.path.dirname(__file__), 'resources', 'SHTF')
    shutil.copytree(data_shtf, os.path.join(cwd.name, 'SHTF'))

    # run NumCalc ram estimation
    subprocess.run([f'{numcalc} -estimate_ram'],
                   cwd=data_cwd, check=True, shell=True,
                   stdout=subprocess.DEVNULL)

    # check if Memory.txt exists
    assert os.path.isfile(os.path.join(data_cwd, 'Memory.txt'))
    # check if output files still exist
    assert os.path.isfile(os.path.join(
        data_cwd, 'be.out', 'be.1', 'pBoundary'))

    # check Memory.txt against reference
    with open(os.path.join(data_cwd, 'Memory.txt'), 'r') as file:
        current = file.readlines()

    with open(os.path.join(
            data_shtf, 'NumCalc', 'source_1', 'Memory.txt'), 'r') as file:
        reference = file.readlines()

    assert current == reference


@pytest.mark.parametrize("boundary_condition", [("rigid"), ("soft")])
@pytest.mark.parametrize("source,range_a", [("plane", (10, -20)),
                                            ("point", (40, -45))])
@pytest.mark.parametrize("bem_method", [("ml-fmm-bem"), ("fmm-bem"), ("bem")])
def test_numcalc_boundary_conditions_sources_types_numerical_methods(
        boundary_condition, source, bem_method, range_a, range_b=(-1, 1)):
    """
    Test if NumCalc and output2hrtf generate correct output by comparing to
    analytical solutions. Tests different single source types, boundary
    conditions and BEM methods.
    """

    # --- Setup ---
    # create temporary directory
    tmp = tempfile.TemporaryDirectory()

    # copy test directory
    shutil.copytree(
        os.path.join(base_dir, 'resources', 'test_numcalc',
                     'project_folder_pspw'),
        os.path.join(tmp.name, 'project'))
    # copy correct input file and rename it to NC.inp
    os.mkdir(os.path.join(tmp.name, 'project', 'NumCalc'))
    os.mkdir(os.path.join(tmp.name, 'project', 'NumCalc', 'source_1'))
    shutil.copyfile(
        os.path.join(base_dir, 'resources', 'test_numcalc', 'ncinp_files',
                     f'NC_{boundary_condition}_{source}_{bem_method}.inp'),
        os.path.join(tmp.name, 'project', 'NumCalc', 'source_1', 'NC.inp'))
    shutil.copyfile(
        os.path.join(base_dir, 'resources', 'test_numcalc', 'parameters',
                     f'parameters_{source}.json'),
        os.path.join(tmp.name, 'project', 'parameters.json'))

    # --- Exercise ---
    # run NumCalc with subprocess
    tmp_path = os.path.join(tmp.name, "project", "NumCalc", "source_1")
    subprocess.run([f'{numcalc}'], cwd=tmp_path, check=True)
    # run output2hrtf
    tmp_path = os.path.join(tmp.name, "project")
    m2s.output2hrtf(tmp_path)

    # --- Verify ---
    # load HRTF data from simulation
    hrtf_sim = utils.hrtf_sofa_to_numpy(
        os.path.join(tmp_path, "Output2HRTF", "HRTF_HorPlane.sofa"))
    # normalize because only relative differences of interest
    hrtf_sim = hrtf_sim[:, :, 0]/numpy.mean(
                numpy.abs(hrtf_sim[numpy.isfinite(hrtf_sim)]))

    # load HRTF data from analytical comparison
    ana_path = os.path.join(
        base_dir, 'resources', 'test_numcalc', 'analytical_references',
        f'ref_{boundary_condition}_{source}.mat')
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
@pytest.mark.parametrize("source,range_a", [
    ("leftear", (40, -40)),
    ("rightear", (40, -40)),
    ("bothears", (40, -40))
])
@pytest.mark.parametrize("bem_method", [("ml-fmm-bem")])
def test_numcalc_ear_source_types(boundary_condition, source, bem_method,
                                  range_a, range_b=(-1, 1)):
    """
    Test if NumCalc and output2hrtf generate correct output by comparing to
    analytical solution. Tests the simulation of HRTF for left, right and both
    ears.
    """

    # --- Setup ---
    # create temporary directory
    tmp = tempfile.TemporaryDirectory()

    # copy basic test directory
    shutil.copytree(
        os.path.join(base_dir, 'resources', 'test_numcalc',
                     'project_folder_ears', 'ears_basic_project'),
        os.path.join(tmp.name, 'project'))

    # copy correct input files for the source type
    shutil.copyfile(
        os.path.join(base_dir, 'resources', 'test_numcalc',
                     'parameters', f'parameters_{source}.json'),
        os.path.join(tmp.name, 'project', 'parameters.json'))
    shutil.copytree(
        os.path.join(base_dir, 'resources', 'test_numcalc',
                     'project_folder_ears', source, 'NumCalc'),
        os.path.join(tmp.name, 'project', 'NumCalc'))

    # --- Exercise ---
    # run NumCalc with subprocess
    tmp_path = os.path.join(tmp.name, "project", "NumCalc", "source_1")
    subprocess.run([f'{numcalc}'], cwd=tmp_path, check=True)
    if source == "bothears":
        tmp_path = os.path.join(tmp.name, "project", "NumCalc", "source_2")
        subprocess.run(
            [f'{numcalc}'], cwd=tmp_path, check=True)
    # run Output2HRTF.py
    tmp_path = os.path.join(tmp.name, "project")
    m2s.output2hrtf(tmp_path)

    # --- Verify ---
    # load HRTF data from simulation as numpy
    hrtf_sim = utils.hrtf_sofa_to_numpy(
        os.path.join(tmp_path, "Output2HRTF", "HRTF_HorPlane.sofa"))
    # normalize because only relative differences of interest
    hrtf_sim = hrtf_sim[:, :, 0]/numpy.mean(
                numpy.abs(hrtf_sim[numpy.isfinite(hrtf_sim)]))

    hrtf_sim = numpy.squeeze(hrtf_sim)

    # load HRTF data from analytical comparison
    ana_path = os.path.join(
        base_dir, 'resources', 'test_numcalc', 'analytical_references',
        f'ref_{boundary_condition}_{source}.mat')
    mat_ana = scipy.io.loadmat(ana_path)
    hrtf_ana = mat_ana['p_total_'+source]
    hrtf_ana = numpy.squeeze(hrtf_ana)
    # normalize because only relative differences of interest
    hrtf_ana = hrtf_ana/numpy.mean(
                numpy.abs(hrtf_ana[numpy.isfinite(hrtf_ana)]))

    # compare
    xyz = mat_ana['XYZ']
    utils.scatter_reference_vs_analytic(
        hrtf_sim, hrtf_ana, xyz[:, 0], xyz[:, 1],
        range_a, range_b, boundary_condition, source, bem_method)

    numpy.testing.assert_allclose(
        numpy.abs(hrtf_sim[numpy.isfinite(hrtf_sim)]),
        numpy.abs(hrtf_ana[numpy.isfinite(hrtf_ana)]), rtol=11.1)
