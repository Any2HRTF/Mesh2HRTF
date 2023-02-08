# %%
import pytest
import subprocess
import tempfile
import shutil
import os
import mesh2scattering as m2s


# directory of this file
base_dir = os.path.dirname(__file__)

# ignore tests for wondows since its difficult to build the exe
if os.name == 'nt':
    numcalc = os.path.join(
        m2s.repository_root(), "NumCalc", "bin", "NumCalc.exe")
else:
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
def test_numcalc_commandline_nitermax(nitermax, use):
    """Test if command line parameter nitermax behaves as expected"""
    # Setup

    # create temporary directory
    tmp = tempfile.TemporaryDirectory()

    # copy test directory
    shutil.copytree(
        os.path.join(
            base_dir, 'resources', 'test_numcalc', 'project_folder_pspw'),
        os.path.join(tmp.name, 'project'))
    # copy correct input file and rename it to NC.inp
    os.mkdir(os.path.join(tmp.name, 'project', 'NumCalc'))
    os.mkdir(os.path.join(tmp.name, 'project', 'NumCalc', 'source_1'))
    shutil.copyfile(
        os.path.join(
            base_dir, 'resources', 'test_numcalc',
            'ncinp_files', 'NC_commandline_parameters.inp'),
        os.path.join(tmp.name, 'project', 'NumCalc', 'source_1', 'NC.inp'))

    if use:
        commandLineArgument = f' -nitermax {nitermax}'
    else:
        commandLineArgument = ''

    # Exercise

    # run NumCalc with subprocess
    tmp_path = os.path.join(tmp.name, "project", "NumCalc", "source_1")
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
    out_filepath = os.path.join(tmp.name, "project", "NumCalc",
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
def test_numcalc_commandline_istart_iend(istart, iend):
    """Test if command line parameters istart and iend behave as expected
    """
    # Setup

    # create temporary directory
    tmp = tempfile.TemporaryDirectory()

    # copy test directory
    shutil.copytree(
        os.path.join(
            base_dir, 'resources', 'test_numcalc', 'project_folder_pspw'),
        os.path.join(tmp.name, 'project'))
    # copy correct input file and rename it to NC.inp
    os.mkdir(os.path.join(tmp.name, 'project', 'NumCalc'))
    os.mkdir(os.path.join(tmp.name, 'project', 'NumCalc', 'source_1'))
    shutil.copyfile(
        os.path.join(
            base_dir, 'resources', 'test_numcalc',
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
    data_shtf = os.path.join(
        os.path.dirname(__file__), 'resources', 'SHTF')
    shutil.copytree(data_shtf, os.path.join(cwd.name, 'SHTF'))

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
