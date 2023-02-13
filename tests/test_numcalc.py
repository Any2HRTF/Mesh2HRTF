import pytest
import subprocess
import shutil
import os
import mesh2scattering as m2s
import glob
import warnings
import numpy.testing as npt
import numpy as np

# directory of this file
base_dir = os.path.dirname(__file__)

# ignore tests for windows since its difficult to build the exe
if os.name == 'nt':
    numcalc = os.path.join(
        m2s.utils.program_root(), "numcalc", "bin", "NumCalc.exe")
    numcalc_path = os.path.dirname(numcalc)
    warnings.warn(
        ('Under Windows the code is not compiling but an executable is '
         f'expected in {numcalc}.'), UserWarning)

else:
    # Build NumCalc locally to use for testing
    numcalc = os.path.join(
        m2s.utils.program_root(), "numcalc", "bin", "NumCalc")
    numcalc_path = numcalc

    if os.path.isfile(numcalc):
        os.remove(numcalc)

    subprocess.run(
        ["make"], cwd=os.path.join(
            m2s.utils.program_root(), "numcalc", "src"), check=True)


def test_import():
    from mesh2scattering import numcalc
    assert numcalc


def test_numcalc_invalid_parameter(capfd):
    """
    Test if NumCalc throws an error in case of invalid command line
    parameter.
    """

    try:
        # run NumCalc with subprocess
        if os.name == 'nt':  # Windows detected
            # run NumCalc and route all printouts to a log file
            subprocess.run(
                f'{numcalc} -invalid_parameter',
                stdout=subprocess.DEVNULL, check=True)
        else:  # elif os.name == 'posix': Linux or Mac detected
            # run NumCalc and route all printouts to a log file
            subprocess.run(
                [f'{numcalc} -invalid_parameter'],
                shell=True, stdout=subprocess.DEVNULL, check=True)
    except subprocess.CalledProcessError:
        _, err = capfd.readouterr()
        assert "NumCalc was called with an unknown parameter or flag." \
            in err
    else:
        ValueError("Num calc did not throw an error")


@pytest.mark.parametrize("nitermax, use", [
    (0, True), (1, True), (2, True), ([], False)])
def test_numcalc_commandline_nitermax(nitermax, use, tmpdir):
    """Test if command line parameter nitermax behaves as expected"""
    # Setup

    # copy test directory
    shutil.copytree(
        os.path.join(
            base_dir, 'resources', 'test_numcalc', 'project_folder_pspw'),
        os.path.join(tmpdir, 'project'))
    # copy correct input file and rename it to NC.inp
    os.mkdir(os.path.join(tmpdir, 'project', 'NumCalc'))
    os.mkdir(os.path.join(tmpdir, 'project', 'NumCalc', 'source_1'))
    shutil.copyfile(
        os.path.join(
            base_dir, 'resources', 'test_numcalc',
            'ncinp_files', 'NC_commandline_parameters.inp'),
        os.path.join(tmpdir, 'project', 'NumCalc', 'source_1', 'NC.inp'))

    if use:
        commandLineArgument = f' -nitermax {nitermax}'
    else:
        commandLineArgument = ''

    # Exercise

    # run NumCalc with subprocess
    tmp_path = os.path.join(tmpdir, "project", "NumCalc", "source_1")
    if os.name == 'nt':  # Windows detected
        # run NumCalc and route all printouts to a log file
        subprocess.run(
            f'{numcalc}{commandLineArgument}',
            stdout=subprocess.DEVNULL, cwd=tmp_path, check=True)
    else:  # elif os.name == 'posix': Linux or Mac detected
        # run NumCalc and route all printouts to a log file
        subprocess.run(
            [f'{numcalc}{commandLineArgument}'],
            shell=True, stdout=subprocess.DEVNULL, cwd=tmp_path, check=True)

    # Verify
    out_filename = 'NC.out'
    out_filepath = os.path.join(tmpdir, "project", "NumCalc",
                                "source_1", out_filename)

    out_file = open(out_filepath)
    out_text = out_file.read()

    if use:
        assert f'CGS solver: number of iterations = {nitermax}' in out_text
        assert 'Warning: Maximum number of iterations is reached!' \
            in out_text
    else:
        assert 'Warning: Maximum number of iterations is reached!' \
            not in out_text


@pytest.mark.parametrize("istart, iend", [
    (False, False), (3, False), (False, 3), (2, 3)])
def test_numcalc_commandline_istart_iend(istart, iend, tmpdir):
    """Test if command line parameters istart and iend behave as expected
    """
    # copy test directory
    shutil.copytree(
        os.path.join(
            base_dir, 'resources', 'test_numcalc', 'project_folder_pspw'),
        os.path.join(tmpdir, 'project'))
    # copy correct input file and rename it to NC.inp
    os.mkdir(os.path.join(tmpdir, 'project', 'NumCalc'))
    os.mkdir(os.path.join(tmpdir, 'project', 'NumCalc', 'source_1'))
    shutil.copyfile(
        os.path.join(
            base_dir, 'resources', 'test_numcalc',
            'ncinp_files', 'NC_commandline_parameters.inp'),
        os.path.join(tmpdir, 'project', 'NumCalc', 'source_1', 'NC.inp'))

    commandLineArgument = ''
    if istart > 0:
        commandLineArgument += f' -istart {istart}'
    if iend > 0:
        commandLineArgument += f' -iend {iend}'

    # Exercise
    # run NumCalc with subprocess
    tmp_path = os.path.join(tmpdir, "project", "NumCalc", "source_1")
    if os.name == 'nt':  # Windows detected
        # run NumCalc and route all printouts to a log file
        subprocess.run(
            f'{numcalc}{commandLineArgument}',
            stdout=subprocess.DEVNULL, cwd=tmp_path, check=True)
    else:  # elif os.name == 'posix': Linux or Mac detected
        # run NumCalc and route all printouts to a log file
        subprocess.run(
            [f'{numcalc}{commandLineArgument}'],
            shell=True, stdout=subprocess.DEVNULL, cwd=tmp_path, check=True)

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

    out_filepath = os.path.join(tmpdir, "project", "NumCalc",
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


def test_numcalc_commandline_estimate_ram(tmpdir):
    """Test NumCalc's RAM estimation using -estimate_ram"""

    # copy test data
    data_cwd = os.path.join(tmpdir, 'SHTF', 'NumCalc', 'source_1')
    data_shtf = os.path.join(
        os.path.dirname(__file__), 'resources', 'SHTF')
    shutil.copytree(data_shtf, os.path.join(tmpdir, 'SHTF'))

    if os.name == 'nt':  # Windows detected
        # run NumCalc and route all printouts to a log file
        subprocess.run(
            f"{numcalc} -estimate_ram",
            stdout=subprocess.DEVNULL, cwd=data_cwd, check=True)
    else:  # elif os.name == 'posix': Linux or Mac detected
        # run NumCalc and route all printouts to a log file
        subprocess.run(
            [f"{numcalc} -estimate_ram"],
            shell=True, stdout=subprocess.DEVNULL, cwd=data_cwd, check=True)

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


def test_defaults(tmpdir):
    """
    Test numcalc manager with default parameters by
    - directly calling the functions
    - running the script that calls the function
    """
    cwd = os.path.dirname(__file__)
    data_shtf = os.path.join(cwd, 'resources', 'SHTF')

    # copy test data to temporary directory and remove test critical data
    shutil.copytree(data_shtf, os.path.join(tmpdir, "SHTF"))
    os.remove(os.path.join(
        tmpdir, "SHTF", "NumCalc", "source_1", "Memory.txt"))
    shutil.rmtree(
        os.path.join(tmpdir, "SHTF", "NumCalc", "source_1", "be.out"))
    shutil.rmtree(os.path.join(tmpdir, "SHTF", "NumCalc", "source_2"))

    # run as function
    m2s.numcalc.manage_numcalc(
        tmpdir, numcalc_path=numcalc_path, wait_time=0)
    # check if files exist
    assert len(glob.glob(os.path.join(tmpdir, "manage_numcalc_*txt")))

    base = os.path.join(tmpdir, "SHTF", "NumCalc", "source_1")
    assert os.path.isfile(os.path.join(base, "Memory.txt"))
    for step in range(1, 61):
        assert os.path.isfile(os.path.join(base, f"NC{step}-{step}.out"))


@pytest.mark.parametrize("folders,issue,errors,nots", (
    # no issues single NC.out filejoin
    [["case_0"], False, [], []],
    # issues in NC.out that are corrected by second file NC1-1.out
    [["case_4"], False, [], []],
    # missing frequencies
    [["case_1"], True,
     ["Frequency steps that were not calculated:\n59, 60"], []],
    # convergence issues
    [["case_2"], True,
     ["Frequency steps that did not converge:\n18, 42"], []],
    # input/mesh issues
    [["case_3"], True,
     ["Frequency steps that were not calculated:\n59, 60",
      "Frequency steps with bad input:\n58"], []],
    # no isses in source 1 but issues in source 2
    [["case_0", "case_1"], True,
     ["Detected issues for source 2",
      "Frequency steps that were not calculated:\n59, 60"],
     ["Detected issues for source 1"]]
))
def test_project_report(folders, issue, errors, nots, tmpdir):
    """Test issues found by the project report"""

    cwd = os.path.dirname(__file__)
    data_nc = os.path.join(cwd, 'resources', 'nc.out')
    # create fake project structure
    os.mkdir(os.path.join(tmpdir, "NumCalc"))
    os.mkdir(os.path.join(tmpdir, "Output2HRTF"))
    shutil.copyfile(os.path.join(data_nc, "parameters.json"),
                    os.path.join(tmpdir, "parameters.json"))
    for ff, folder in enumerate(folders):
        shutil.copytree(os.path.join(data_nc, folder),
                        os.path.join(tmpdir, "NumCalc", f"source_{ff + 1}"))

    # run the project report
    issues, report = m2s.output.write_output_report(tmpdir)

    # test the output
    assert issues is issue
    for error in errors:
        assert error in report
    for no in nots:
        assert no not in report
    if issue:
        assert os.path.isfile(os.path.join(
            tmpdir, "Output2HRTF", "report_issues.txt"))
        assert ("For more information check Output2HRTF/report_source_*.csv "
                "and the NC*.out files located at NumCalc/source_*") in report
    else:
        assert not os.path.isfile(os.path.join(
            tmpdir, "Output2HRTF", "report_issues.txt"))


@pytest.mark.parametrize("boundary,grid", [
    (True, True), (True, False), (False, True)])
def test_purge_outputs_numcalc_data(boundary, grid, tmpdir):
    """Test purging the raw NumCalc output"""
    cwd = os.path.dirname(__file__)
    data_shtf = os.path.join(cwd, 'resources', 'SHTF')

    # copy required data to temporary directory
    shutil.copytree(data_shtf, os.path.join(tmpdir, "SHTF"))

    m2s.numcalc.remove_outputs(os.path.join(tmpdir, "*"), boundary, grid)

    for source in glob.glob(
            os.path.join(tmpdir, "SHTF", "NumCalc", "source_*")):
        if boundary and grid:
            assert not os.path.isdir(os.path.join(source, "be.out"))
        elif boundary:
            assert os.path.isdir(os.path.join(source, "be.out"))
            for be in glob.glob(os.path.join(source, "be.out", "be.*")):
                assert glob.glob(os.path.join(be, "*Boundary")) == []
        elif grid:
            assert os.path.isdir(os.path.join(source, "be.out"))
            for be in glob.glob(os.path.join(source, "be.out", "be.*")):
                assert glob.glob(os.path.join(be, "*EvalGrid")) == []


@pytest.mark.parametrize("hrtf,vtk,reports", [
    (False, True, False), (True, False, True)])
def test_purge_outputs_output_data(hrtf, vtk, reports, tmpdir):
    """Test purging the processed data in Output2HRTF"""
    cwd = os.path.dirname(__file__)
    data_shtf = os.path.join(cwd, 'resources', 'SHTF')
    shutil.copytree(data_shtf, os.path.join(tmpdir, "SHTF"))
    folder = os.path.join(tmpdir, "SHTF", "Output2HRTF")

    m2s.numcalc.remove_outputs(
        os.path.join(tmpdir, "*"), hrtf=hrtf, vtk=vtk, reports=reports)

    assert os.path.isfile(
        os.path.join(folder, "HRTF_FourPointHorPlane_r100cm.sofa")) \
        == (not hrtf)

    assert os.path.isdir(os.path.join(folder, "vtk")) \
        == (not vtk)

    assert os.path.isfile(os.path.join(folder, "report_source_1.csv")) == \
        (not reports)

    assert os.path.isfile(os.path.join(folder, "report_source_2.csv")) == \
        (not reports)


def test_read_ram_estimates():

    estimates = m2s.numcalc.read_ram_estimates(os.path.join(
        os.path.dirname(__file__), "resources", "SHTF", "NumCalc", "source_1"))

    assert isinstance(estimates, np.ndarray)
    assert estimates.shape == (60, 3)
    npt.assert_allclose([1.00000e+00, 1.00000e+02, 4.16414e-02], estimates[0])
    npt.assert_allclose([6.00000e+01, 6.00000e+03, 7.22010e-02], estimates[-1])


def test_read_ram_estimates_assertions():
    """test assertions for read_ram_estimates"""

    with pytest.raises(ValueError, match="does not contain a Memory.txt"):
        m2s.numcalc.read_ram_estimates(os.getcwd())
