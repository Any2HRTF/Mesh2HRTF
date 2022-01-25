import subprocess
import tempfile
import shutil
import os

blender_path = [os.path.join('/home', 'matheson', 'Apps', 'blender-2.91.0',
                             'blender')]


def test_blender_export():
    """ test the exportMesh2HRTF Blender plugin """

    for i in range(len(blender_path)):
        # create a temporary directory
        tmp = tempfile.TemporaryDirectory(dir=os.getcwd())

        # copy test directory
        shutil.copytree(os.path.join(os.path.dirname(__file__),
                                    'test_blender_export_project'),
                        os.path.join(tmp.name, 'project'))

        tmp_path = os.path.join(tmp.name, 'project')
        blender_file_path = os.path.join(tmp_path, '3dModel.blend')
        python_file_path = os.path.join(tmp_path, 'blender_script.py')

        # run exportMesh2HRTF from Blender with subprocess
        subprocess.run([blender_path[i], blender_file_path,
                        "--background", "--python",
                        python_file_path],
                    cwd=tmp_path, check=True, capture_output=True)

        # run NumCalc with subprocess
        tmp_path = os.path.join(tmp.name, "project", "NumCalc", "source_1")
        subprocess.run(["NumCalc"], cwd=tmp_path, check=True, capture_output=True)

        # run Output2HRTF.py
        tmp_path = os.path.join(tmp.name, "project")
        subprocess.run(["python", "Output2HRTF.py"], cwd=tmp_path, check=True,
                    capture_output=True)

# subprocess.run(["/home/matheson/Apps/blender-2.91.0/blender 3dModel.blend
# --background --python blender_script.py"], cwd=tmp_path)
