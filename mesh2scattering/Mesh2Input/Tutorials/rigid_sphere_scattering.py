# Export the 'Scattering from a rigid sphere' tutorial
import bpy

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

# add icosphere
bpy.ops.mesh.primitive_ico_sphere_add(subdivisions=6, radius=0.1)
bpy.data.objects['Icosphere'].name = 'Reference'

# add point source
bpy.ops.object.light_add(location=(0.0, 0.2, 0.0))
bpy.data.objects['Point'].name = 'Point source'


# export the project ----------------------------------------------------------
bpy.ops.mesh2input.inp(
    filepath=file_path,
    programPath=program_path,
    title='Scattering from a rigid sphere by a point source',
    method='ML-FMM BEM',
    sourceType='Point source',
    pictures=False,
    unit='m',
    speedOfSound='343.18',
    evaluationGrids='PlaneHorizontal',
    minFrequency=1000,
    maxFrequency=16000,
    frequencyVectorValue=1000
)
