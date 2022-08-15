# %%
import pytest
import os
import subprocess
import tempfile
import utils
import mesh2hrtf as m2h

base_dir = os.path.dirname(__file__)
install_script = os.path.join(base_dir, 'resources', 'install_addons.py')

# define and check paths to your Blender versions (only use one blender)
blender_path = utils.blender_paths(2)[0]

# addons to be installed
addons = [
    (os.path.join(base_dir, 'mesh2hrtf', 'Mesh2Input', 'mesh2input.py'),
     'mesh2input'),
    (os.path.join(base_dir, 'mesh2hrtf', 'Mesh2Input', 'Meshes',
                  'AssignMaterials', 'AssignMaterials.py'),
     None)]

# generate script for installing addons
utils.blender_addons_installer(
    blender_path[0], os.path.join(blender_path[0], blender_path[1]),
    addons, install_script)

# install addons
subprocess.run(
        [os.path.join(blender_path[0], 'blender'), '--background',
         '--python', install_script],
        cwd=base_dir, check=True, capture_output=True)


@pytest.mark.parametrize('tutorial', ['rigid_sphere_scattering.py'])
def test_test(tutorial):

    # directory for testing the tutorial
    tmp = tempfile.TemporaryDirectory()
    print(f"preparing {tutorial} in {tmp.name}")

    # read tutorial file
    tutorial_file = os.path.join(
        base_dir, '..', 'mesh2hrtf', 'Mesh2Input', 'Tutorials', tutorial)
    with open(tutorial_file, 'r') as file:
        script = ''.join(file.readlines())

    # change paths
    script = script.replace(
        "path/to/your/project_folder",
        os.path.join(tmp.name, tutorial[:-3]))
    script = script.replace(
        "path/to/your/Mesh2HRTF/mesh2hrtf",
        os.path.join(base_dir, '..', 'mesh2hrtf'))

    # save tutorial to temp dir
    with open(os.path.join(tmp.name, tutorial), 'w') as file:
        file.write(script)

    # export the project folder
    print("exporting")
    subprocess.run(
        [os.path.join(blender_path[0], 'blender'), '--background',
         '--python', os.path.join(tmp.name, tutorial)],
        cwd=tmp.name, check=True, capture_output=True)

    # run manage_numcalc
    print("running NumCalc")
    m2h.manage_numcalc(os.path.join(tmp.name, tutorial[-3]))

    # run manage_numcalc
    print("running output2hrtf")
    m2h.output2hrtf(os.path.join(tmp.name, tutorial[-3]))
