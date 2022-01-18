import bpy
import os

# Mesh2HRTF export paths
exportPath = os.getcwd()
mesh2HRTFpath = os.path.join(os.getcwd(), "..", "..", "mesh2hrtf")
eval_grid = 'HorPlane'
# os.path.join(os.getcwd(), "..", "export_project_test", "Evaluation Grids", "HorPlane")

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
    title="head-related transfer functions",
    filepath=exportPath,
    programPath=mesh2HRTFpath,
    materialSearchPaths='None',
    method='ML-FMM BEM',
    sourceType='Point source',
    pictures=False,
    unit='m',
    evaluationGrids=eval_grid,
    minFrequency=1000,
    maxFrequency=1000,
    frequencyVectorType='Step size',
    frequencyVectorValue=1000,
    reference=False,
    computeHRIRs=False,
    speedOfSound="343",
    densityOfMedium='1.1839')

# purge unused data ---------------------------------------
ret = bpy.ops.outliner.orphans_purge()
while ret != {'CANCELLED'}:
    ret = bpy.ops.outliner.orphans_purge()
