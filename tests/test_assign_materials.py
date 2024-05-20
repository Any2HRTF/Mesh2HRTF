import pytest
import os
from tempfile import TemporaryDirectory
import subprocess
import utils

# blender path (uses only the latest)
blender_path = utils.blender_paths()[-1][0]
# location of assign materials script and reference data
base_dir = os.path.dirname(__file__)
assign_script = os.path.join(
    base_dir, '..', 'mesh2hrtf', 'Mesh2Input', 'Meshes',
    'AssignMaterials', 'AssignMaterials.py')
ref_dir = os.path.join(base_dir, 'resources', 'assign_materials')


@pytest.mark.parametrize('ear', ['Left ear', 'Right ear', 'Both ears'])
def test_assign_materials(ear):

    with TemporaryDirectory() as tmp_dir:

        # script for assigning materials to an ico sphere object
        savename = 'assign_materials ' + ear
        assign_materials = (
            "import bpy\n"
            "import os\n\n"
            "bpy.ops.object.delete()\n"
            "bpy.ops.mesh.primitive_ico_sphere_add()\n\n"
            f"bpy.ops.object.assignmaterials(ear='{ear}')\n\n"
            "bpy.ops.wm.obj_export(\n"
            f"    filepath=os.path.join('{tmp_dir}', '{savename}' + '.obj'),\n"
            "    export_selected_objects=True,\n"
            "    export_uv=False, export_normals=False,\n"
            "    export_materials=True)\n")

        with open(os.path.join(tmp_dir, savename + '.py'), 'w') as file:
            file.writelines(assign_materials)

        # run blender and assign materials scripts
        result = subprocess.run(  # noqa (result can be used for debugging)
            [os.path.join(blender_path, 'blender'), '--background',
             '--python', assign_script,
             '--python', os.path.join(tmp_dir, savename + '.py')],
             stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)

        # compare to generated data to reference
        for file_type in ['.obj', '.mtl']:
            with open(os.path.join(tmp_dir, savename + file_type)) as file:
                test = file.readlines()
            with open(os.path.join(ref_dir, savename + file_type)) as file:
                reference = file.readlines()

            # skip first to lines with blender version specific data
            assert test[2:] == reference[2:]
