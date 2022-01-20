import subprocess
import tempfile
import shutil
import os


def test_two_sources():
    """
    test Output2HRTF.py handling two sound sources ("Both ears" condition)
    """
    # Setup

    # create temporary directory
    tmp = tempfile.TemporaryDirectory(dir=os.getcwd())

    # copy test directory
    shutil.copytree(os.path.join(os.path.dirname(__file__),
                                 'test_2sources_project'),
                    os.path.join(tmp.name, 'project'))

    # Exercise

    # # run NumCalc with subprocess
    # for iSource in (1, 2):
    #     tmp_path = os.path.join(tmp.name, "project", "NumCalc", f"source_{iSource+1}")
    #     subprocess.run(["NumCalc"], cwd=tmp_path, check=True)

    # run Output2HRTF.py
    tmp_path = os.path.join(tmp.name, "project")
    subprocess.run(["python", "Output2HRTF.py"], cwd=tmp_path, check=True)

# Only used for debugging
# test_blender_export()
# test_build()
# test_numcalc("rigid", "point", "ml-fmm-bem", (10, -20), range_b=(-1, 1))
# test_two_sources()
