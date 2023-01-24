# Export the 'Scattering from a rigid sphere' tutorial

import numpy as np
import os
import bpy

# user parameters -------------------------------------------------------------
# this is the folder to which the project is exported
data_path = '/home/anne/sciebo/2021_DFG-Projekt/data'
project_name_out = 'sine_test'
file_path = os.path.join(data_path, 'mesh2hrtf_results', project_name_out)

# this is the folder mesh2scattering inside the mesh2scattering git repository
program_path = '/home/anne/git/Mesh2scattering/mesh2scattering'

# this is the folder where the meshes are located, one called sample.stl the
# other called reference.stl
project_name_in = 'sine_10k'
meshes_path = os.path.join(data_path, 'meshes')

# this defines the source positions
source_distance = 10
source_phi_deg = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90]
source_theta_deg = [10, 20, 30, 40, 50, 60, 70, 80]

# this defines the receiver points, you can find them in the 
# EvaluationGrids/Data folder or create your own
evaluation_grid = 'coords5_'

# to create the project run the line below two times in the consol
# blender --background --python Mesh2Input/Tutorials/mesh2scattering.py


# prepare the scene -----------------------------------------------------------
def create_scene_with_stl(
        stl_path, stl_name, source_distance, source_theta_deg, source_phi_deg):
    name = stl_name.lower()
    folder_out = os.path.join(file_path, name)
    if os.path.exists(folder_out):
        return False
    # prepare the scene -------------------------------------------------------
    # force object mode
    bpy.ops.object.mode_set(mode='OBJECT', toggle=False)

    # set cursor to origin
    bpy.context.scene.cursor.location = (0.0, 0, 0.0)

    # remove cube
    bpy.data.objects['Cube'].select_set(True)
    bpy.ops.object.delete()

    # add mesh
    bpy.ops.import_mesh.stl(
        filepath=stl_path, axis_forward='Y', axis_up='Z')
    bpy.data.objects[stl_name].name = 'Reference'

    # add point source
    id = 0
    for theta in source_theta_deg:
        for phi in source_phi_deg:
            theta_rad = theta * np.pi / 180.
            phi_rad = phi * np.pi / 180.
            x = source_distance * np.sin(theta_rad) * np.cos(phi_rad)
            y = source_distance * np.sin(theta_rad) * np.sin(phi_rad)
            z = source_distance * np.cos(theta_rad)

            bpy.ops.object.light_add(location=(x, y, z))
            point_name = 'Point source' if id == 0 else f'Point source{id}'
            bpy.data.objects['Point'].name = point_name

    # export the project ------------------------------------------------------
    bpy.ops.mesh2input.inp(
        filepath=folder_out,
        programPath=program_path,
        title=name,
        method='ML-FMM BEM',
        sourceType='Point source',
        pictures=False,
        unit='m',
        speedOfSound='343.18',
        evaluationGrids=evaluation_grid,
        frequencyVectorType='Step size',
        minFrequency=1000,
        maxFrequency=16000,
        frequencyVectorValue=1000
    )
    bpy.ops.wm.quit_blender()
    return True


sample_path = os.path.join(meshes_path, project_name_in, 'sample.stl')
ref_path = os.path.join(meshes_path, project_name_in, 'reference.stl')
if not os.path.exists(file_path):
    os.mkdir(file_path)

for itype in [0, 1]:
    if itype == 0:
        created = create_scene_with_stl(
            sample_path, 'Sample',
            source_distance, source_theta_deg, source_phi_deg)
        if created:
            break

    if itype == 1:
        created = create_scene_with_stl(
            ref_path, 'Reference', source_distance, source_theta_deg, [0])
        if created:
            break
