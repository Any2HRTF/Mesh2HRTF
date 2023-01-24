"""
Test manage_numcalc with default parameters. Testing parameter combinations is
not done. Be careful when changing manage_numcalc()!
"""
import pytest
import subprocess
from tempfile import TemporaryDirectory
import shutil
import os
import glob
import mesh2scattering as m2s

cwd = os.path.dirname(__file__)
data_shtf = os.path.join(cwd, 'resources', 'SHTF')

# Build NumCalc locally to use for testing
tmp = TemporaryDirectory()
numcalc = os.path.join(tmp.name, "NumCalc", "bin", "NumCalc")

shutil.copytree(
    os.path.join(cwd, "..", "mesh2scattering", "NumCalc"),
    os.path.join(tmp.name, "NumCalc"))

if os.path.isfile(numcalc):
    os.remove(numcalc)

subprocess.run(
    ["make"], cwd=os.path.join(tmp.name, "NumCalc", "src"), check=True)


@pytest.mark.parametrize("mode", ("function", "script"))
def test_defaults(mode):
    """
    Test numcalc manager with default parameters by
    - directly calling the functions
    - running the script that calls the function
    """

    # copy test data to temporary directory and remove test critical data
    temp = TemporaryDirectory()
    shutil.copytree(data_shtf, os.path.join(temp.name, "SHTF"))
    os.remove(os.path.join(
        temp.name, "SHTF", "NumCalc", "source_1", "Memory.txt"))
    shutil.rmtree(
        os.path.join(temp.name, "SHTF", "NumCalc", "source_1", "be.out"))
    shutil.rmtree(os.path.join(temp.name, "SHTF", "NumCalc", "source_2"))

    if mode == "function":
        # run as function
        m2s.manage_numcalc(temp.name, numcalc_path=numcalc, wait_time=0)
    elif mode == "script":
        # run as script
        script_path = os.path.join(cwd, "..", "mesh2scattering", "NumCalc")
        subprocess.run([
            (f'python manage_numcalc_script.py --project_path {temp.name} '
             f'--numcalc_path {numcalc} --wait_time 0 --confirm_errors False')
            ],
             cwd=script_path, check=True, shell=True)

    # check if files exist
    assert len(glob.glob(os.path.join(temp.name, "manage_numcalc_*txt")))

    base = os.path.join(temp.name, "SHTF", "NumCalc", "source_1")
    assert os.path.isfile(os.path.join(base, "Memory.txt"))
    for step in range(1, 61):
        assert os.path.isfile(os.path.join(base, f"NC{step}-{step}.out"))
