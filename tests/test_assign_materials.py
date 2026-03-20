import pytest
import os
import numpy as np
import numpy.testing as npt
from tempfile import TemporaryDirectory
import subprocess
import warnings
import utils

# blender path (uses only the latest)
blender_paths = utils.blender_paths()
# location of assign materials script and reference data
base_dir = os.path.dirname(__file__)
assign_script = os.path.join(
    base_dir, '..', 'mesh2hrtf', 'Mesh2Input', 'Meshes',
    'AssignMaterials', 'AssignMaterials.py')
ref_dir = os.path.join(base_dir, 'resources', 'assign_materials')


@pytest.mark.parametrize('blender_path', blender_paths)
@pytest.mark.parametrize('ear', ['Left ear', 'Right ear', 'Both ears'])
def test_assign_materials(blender_path, ear):

    blender_path = blender_path[0]

    with TemporaryDirectory() as tmp_dir:

        # script for assigning materials to an ico sphere object
        savename = 'assign_materials ' + ear
        assign_materials = (
            "import bpy\n"
            "import numpy as np\n"
            "import os\n\n"
            # required variables
            f"tmp_dir = '{tmp_dir}'\n"
            f"savename = '{savename}'\n"
            # create ico sphere and assign materials
            "bpy.ops.object.delete()\n"
            "bpy.ops.mesh.primitive_ico_sphere_add()\n\n"
            f"bpy.ops.object.assignmaterials(ear='{ear}')\n\n"
            # find and save indices of faces assigned to each material
            "obj = bpy.context.object\n"
            "for mat_name in ['Skin', 'Left ear', 'Right ear']:\n"
            "    mat_idx = next(\n"
            "        (i for i, slot in enumerate(obj.material_slots)\n"
            "         if slot.material and slot.material.name == mat_name),\n"
            "        None)\n"
            "    face_indices = [p.index for p in obj.data.polygons if\n"
            "                    p.material_index == mat_idx]\n"
            "    np.savetxt(os.path.join(\n"
            "        tmp_dir, f'{savename}_face_indices_of_{mat_name}.csv'),\n"
            "        face_indices, fmt='%d', delimiter=',')\n"
            )

        with open(os.path.join(tmp_dir, savename + '.py'), 'w') as file:
            file.writelines(assign_materials)

        # run blender and assign materials scripts
        result = subprocess.run(  # noqa (result can be used for debugging)
            [os.path.join(blender_path, 'blender'), '--background',
             '--python', assign_script,
             '--python', os.path.join(tmp_dir, savename + '.py')],
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)

        for mat_name in ['Skin', 'Left ear', 'Right ear']:
            filename = f'{savename}_face_indices_of_{mat_name}.csv'

            with warnings.catch_warnings():
                warnings.simplefilter("ignore", UserWarning)
                test = np.loadtxt(
                    os.path.join(tmp_dir, filename), delimiter=',')
                reference = np.loadtxt(
                    os.path.join(ref_dir, filename), delimiter=',')

            npt.assert_equal(test, reference)
