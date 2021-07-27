import subprocess
import tempfile
import shutil
import os

# create a temporary directory
tmp = tempfile.TemporaryDirectory(dir=os.getcwd())


# def test_build():
#     """ test if make for NumCalc works """

#     shutil.copytree("../mesh2hrtf/NumCalc/Source", tmp.name+"/NumCalc")
#     tmp_path = os.path.join(tmp.name, "NumCalc")
#     subprocess.run(["make"], cwd=tmp_path, check=True)


def test_1ear_rigid_planewave_bem():
    """ 
    test if NumCalc and Output2HRTF.py generate correct output by comparing to
    analytical solution
    """

    # Setup

    # copy test directory
    shutil.copytree("/home/matheson/Documents/jthomsen/Test_material/test_prototype_bone", tmp.name+"/project")
    # copy correct input file and rename it to NC.inp
    shutil.copyfile("/home/matheson/Documents/jthomsen/Test_material/test_prototype/NumCalc/CPU_1_Core_1/NC_rigid_plane_bem.inp",
                    tmp.name+"/project/NumCalc/CPU_1_Core_1/NC.inp")

    # Exercise

    # run NumCalc with subprocess
    tmp_path = os.path.join(tmp.name, "project", "NumCalc", "CPU_1_Core_1")
    subprocess.run(["NumCalc"], cwd=tmp_path, check=True)
    # run Output2HRTF.py
    tmp_path = os.path.join(tmp.name, "project")
    subprocess.run(["python", "Output2HRTF.py"], cwd=tmp_path, check=True)

    # Verify

#test_1ear_rigid_planewave_bem()