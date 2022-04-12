import bpy
import sys

if len(sys.argv) > 5:
    print('You have specified too many arguments in Blender command line interface!')
    sys.exit()

if len(sys.argv) < 5:
    print('You have specified too few arguments in Blender command line interface!')
    sys.exit()

# Define paths (done in testing by substituting the strings)
exportPath = None
programPath = None
addonFile = None
addonPath = None

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
    programPath=programPath,
    # additional kwargs added in testing
)

# purge unused data ---------------------------------------
ret = bpy.ops.outliner.orphans_purge()
while ret != {'CANCELLED'}:
    ret = bpy.ops.outliner.orphans_purge()
