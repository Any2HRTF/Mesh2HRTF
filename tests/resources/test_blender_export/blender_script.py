import bpy

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

# save Mesh2HRTF project ----------------------------------
bpy.ops.export_mesh2hrtf.inp(
    filepath=exportPath,
    programPath=programPath,
    # additional kwargs added in testing
)
