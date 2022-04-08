import bpy
import os
import sys

if len(sys.argv) > 5:
    print('You have specified too many arguments in Blender command line interface!')
    sys.exit()

if len(sys.argv) < 5:
    print('You have specified too few arguments in Blender command line interface!')
    sys.exit()

# # remove old data
# if os.path.exists(dir_export):
#     shutil.rmtree(dir_export)
# os.mkdir(dir_export)

# Define paths
exportPath = os.getcwd()
mesh2HRTFpath = os.path.join(os.getcwd(), "..", "..", "mesh2hrtf")
eval_grid = 'HorPlane'
# os.path.join(os.getcwd(), "..", "export_project_test", "Evaluation Grids", "HorPlane")
addonPath = 
if not os.path.isdir(addonPath):
    print('The Blender addon path specified does not exist')
    sys.exit()
# '/home/matheson/Apps/blender-2.91.0/2.91/scripts/addons'
addonFile = os.path.join(os.getcwd(), "..", "..", "mesh2hrtf", "Mesh2Input",
                         "exportMesh2HRTF.py")
# '/home/matheson/Apps/mesh2hrtf-git/mesh2hrtf/Mesh2Input/exportMesh2HRTF.py'

# re-install export addon
bpy.context.preferences.filepaths.script_directory = addonPath
bpy.utils.refresh_script_paths()
bpy.ops.preferences.addon_install(overwrite=True, filepath=addonFile)
bpy.ops.preferences.addon_enable(module='exportMesh2HRTF')

# switch to object mode
bpy.ops.object.mode_set(mode='OBJECT', toggle=False)
# select mesh
obj = bpy.data.objects['Reference']
obj.select_set(True)

# purge unused data ---------------------------------------
ret = bpy.ops.outliner.orphans_purge()
while ret != {'CANCELLED'}:
    ret = bpy.ops.outliner.orphans_purge()

# apply transforms
bpy.ops.object.transform_apply(location=True)
bpy.ops.object.transform_apply(rotation=True)
bpy.ops.object.transform_apply(scale=True)

# save Mesh2HRTF project ----------------------------------
bpy.ops.export_mesh2hrtf.inp(
    filepath=exportPath,
export_args
    # title="head-related transfer functions",
    # filepath=exportPath,
    # programPath=mesh2HRTFpath,
    # materialSearchPaths='None',
    # method='ML-FMM BEM',
    # sourceType='Point source',
    # pictures=False,
    # unit='m',
    # evaluationGrids=eval_grid,
    # minFrequency=1000,
    # maxFrequency=1000,
    # frequencyVectorType='Step size',
    # frequencyVectorValue=1000,
    # reference=False,
    # computeHRIRs=False,
    # speedOfSound="343",
    # densityOfMedium='1.1839')

# purge unused data ---------------------------------------
ret = bpy.ops.outliner.orphans_purge()
while ret != {'CANCELLED'}:
    ret = bpy.ops.outliner.orphans_purge()
