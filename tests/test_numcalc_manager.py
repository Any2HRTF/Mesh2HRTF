"""
Test numcalc_manger with default parameters. Testing parameter combinations is
not done. Be carefull when changing numcalc_manger()!
"""
import pytest
import subprocess
from tempfile import TemporaryDirectory
import shutil
import os
import glob
import mesh2hrtf as m2h

cwd = os.path.dirname(__file__)
data_shtf = os.path.join(cwd, 'resources', 'SHTF')


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
        m2h.numcalc_manager(temp.name, wait_time=0)
    elif mode == "script":
        # run as script
        script_path = os.path.join(cwd, "..", "mesh2hrtf", "NumCalc")
        subprocess.run([
            (f'python numcalc_manager.py --project_path {temp.name} '
             '--wait_time 0 --confirm_errors False')],
             cwd=script_path, check=True, shell=True)

    # check if files exist
    assert len(glob.glob(os.path.join(temp.name, "numcalc_manager_*txt")))

    base = os.path.join(temp.name, "SHTF", "NumCalc", "source_1")
    assert os.path.isfile(os.path.join(base, "Memory.txt"))
    for step in range(1, 61):
        assert os.path.isfile(os.path.join(base, f"NC{step}-{step}.out"))
