# Export the 'HRTF' tutorial
import bpy
import os

# user parameters -------------------------------------------------------------
# this is the folder to which the project is exported
file_path = 'path/to/your/project_folder'
# this is the folder mesh2hrtf inside the Mesh2HRTF git repository
program_path = 'path/to/mesh2scattering/mesh2scattering'


# prepare the scene -----------------------------------------------------------
# force object mode
bpy.ops.object.mode_set(mode='OBJECT', toggle=False)

# set cursor to origin
bpy.context.scene.cursor.location = (0.0, 0, 0.0)

# remove cube
bpy.data.objects['Cube'].select_set(True)
bpy.ops.object.delete()

# load head mesh
mesh_file = os.path.join(program_path, 'Mesh2Input', 'Meshes', 'Data',
                         'example_head_max_f_16_kHz.ply')
bpy.ops.import_mesh.ply(filepath=mesh_file)

# assign materials Skin, Left ear, and Right ear
bpy.ops.object.assignmaterials()

# export the project ----------------------------------------------------------
bpy.ops.mesh2input.inp(
    filepath=file_path,
    programPath=program_path,
    title='Uniform HRTF calculation',
    method='ML-FMM BEM',
    sourceType='Both ears',
    pictures=False,
    reference=True,
    computeHRIRs=True,
    speedOfSound='343.18',
    evaluationGrids='PlaneHorizontal',
    minFrequency=100,
    maxFrequency=16000,
    frequencyVectorValue=100
)
