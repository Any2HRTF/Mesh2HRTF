
import os
import bpy
import bmesh
import datetime
import math
import shutil
from math import pi
import json
from bpy.props import StringProperty, BoolProperty, EnumProperty, \
    IntProperty, FloatProperty
from bpy_extras.io_utils import ExportHelper
import numpy as np

bl_info = {
    "name": "Mesh2HRTF export add-on",
    "author": "The Mesh2HRTF developers",
    "version": (1, 0, 0),
    "blender": (2, 80, 0),
    "location": "File > Export",
    "description": "Export Blender scene as Mesh2HRTF project",
    "warning": "",
    "wiki_url": "",
    "tracker_url": "",
    "support": "COMMUNITY",
    "category": "Import-Export"}


class ExportMesh2HRTF(bpy.types.Operator, ExportHelper):
    '''Export Blender scene as Mesh2HRTF project'''
    bl_idname = "mesh2input.inp"
    bl_label = "Export Mesh2HRTF"

    filename_ext = ""
    filter_glob: StringProperty(default="", options={'HIDDEN'})

    # general settings --------------------------------------------------------
    title: StringProperty(
        name="Title",
        description="Title of the project",
        default="Head-Related Transfer Functions",
        )
    method: EnumProperty(
        name="BEM method",
        description="Method for numerical simulation",
        items=[('BEM', 'BEM', 'Traditional BEM'),
               ('SL-FMM BEM', 'SL-FMM BEM',
                'Singlelevel fast-multipole method BEM'),
               ('ML-FMM BEM', 'ML-FMM BEM',
                'Multilevel fast-multipole method BEM')],
        default='ML-FMM BEM',
        )
    sourceType: EnumProperty(
        name="Source type",
        description="Method for numerical simulation",
        items=[('Both ears', 'Both ears',
                    ("Mesh elements with user assigned material 'Left ear' "
                     "and 'Right ear' act as the source")),
               ('Left ear', 'Left ear',
                    ("Mesh elements with user assigned material 'Left ear' "
                     "act as the source")),
               ('Right ear', 'Right ear',
                    ("Mesh elements with user assigned material 'Right ear' "
                     "act as the source")),
               ('Point source', 'Point source',
                    ("Analytical point source. Coordinates taken from user "
                     "placed point light named 'Point source'")),
               ('Plane wave', 'Plane wave',
                    ("Analytical plane wave. Coordinates taken from the "
                    "location (not rotation) of the user placed "
                     "area light named 'Plane wave'"))],
        default='Both ears',
        )
    programPath: StringProperty(
        name="Mesh2HRTF-path",
        description=("Path to folder containing 'NumCalc' and other folders"
                     "(used to copy files to project folder during export)"),
        default=r"/Users/anne/git/Mesh2scattering/mesh2scattering",
        )
    pictures: BoolProperty(
        name="Pictures",
        description="Render pictures of the 3D mesh",
        default=True,
        )
    # post-processing ---------------------------------------------------------
    reference: BoolProperty(
        name="Reference",
        description=("Reference HRTF to the center of the head according to "
                     "the classic HRTF definition. For the HRTF definition "
                     "see https://doi.org/10.1016/0003-682X(92)90046-U"),
        default=False,
        )
    computeHRIRs: BoolProperty(
        name="Compute HRIRs",
        description=("Compute HRIRs by inverse Fourier Transform. This "
                     "requires referencing (see above) and frequencies "
                     "between f and fs/2 in steps of f, with f>0 and fs the "
                     "sampling frequency"),
        default=False,
        )
    # constants ---------------------------------------------------------------
    unit: EnumProperty(
        name="Unit",
        description="Unit of the 3D scene.",
        items=[('m', 'm', 'Meter'), ('mm', 'mm', 'Millimeter')],
        default='mm',
        )
    speedOfSound: StringProperty(
        name="c (m/s)",
        description=("Speed of sound in m/s "
                     "(Pass value as a string to avoid rounding errors)"),
        default="346.18",
        )
    densityOfMedium: StringProperty(
        name="rho (kg/m^3)",
        description=("Density of air in kg/m^3 "
                     "(Pass value as a string to avoid rounding errors)"),
        default="1.1839",
        )
    # evaluation grids --------------------------------------------------------
    evaluationGrids: StringProperty(
        name="Name",
        description=("Name of evaluation grid inside "
            "Mesh2Input/EvaluationsGrids/Data or absolute path to user grid. "
            "Multiple grids can be separated by semicolons (;)"),
        default='Default',
        )
    # material seach paths ----------------------------------------------------
    materialSearchPaths: StringProperty(
        name="Path(s)",
        description=("Absolute path(s) to folders that contain material data. "
            "Multiple paths can be separated by semicolons (;). The "
            "Mesh2HRTF path is added by default to the end of the list"),
        default='None',
        )
    # Frequency selection -----------------------------------------------------
    minFrequency: FloatProperty(
        name="Min. frequency",
        description=("Minimum frequency in Hz to be simulated. Can be 0 for "
            "constructing the frequency vector. But the 0 will not be "
            "simulated"),
        default=100,
        min=0,
        max=24000,
        )
    maxFrequency: FloatProperty(
        name="Max. frequency",
        description="Maximum frequency in Hz to be simulated",
        default=20000,
        min=1,
        max=96000,
        )
    frequencyVectorType: EnumProperty(
        name="Frequencies",
        description=("Select how the frequencies are distributed between min. "
                     "and max. frequency with the paramter 'Value' (below)."),
        items=[('Step size', 'Step size',
                    ("Simulate frequencies between the min. and max. "
                     "with a fixed step size.")),
               ('Num steps', 'Num steps',
                    "Simulate N frequencies the min. and max. frequency."),
               ('Nominal n-th octave', 'Nominal n-th octave',
                    "Simulate the nominal (rounded) center frequencies "
                    "specified in the standard. Nominal frequencies are only "
                    "returned for octave bands and third octave bands"),
               ('Exact n-th octave', 'Exact n-th octave',
                    "The exact center frequencies, resulting in a uniform "
                    "distribution of frequency bands over the frequency range.")],
        default='Step size',
        )
    frequencyVectorValue: FloatProperty(
        name="Value",
        description=("Value for distributing the frequencies according to the "
                     "parameter 'Frequencies' (above)."),
        default=100,
        min=0,
        max=24000,
        )

    @classmethod
    def poll(cls, context):
        return context.active_object is not None

    def execute(self, context):
        filepath = self.filepath
        filepath = bpy.path.ensure_ext(filepath, self.filename_ext)
        keywords = self.as_keywords(ignore=("check_existing", "filter_glob"))
        return ExportMesh2HRTF.save(self, context, **keywords)

    def draw(self, context):
        layout = self.layout
        layout.use_property_split = True
        layout.use_property_decorate = False

        # general settings
        layout.label(text="General:")
        row = layout.row()
        row.prop(self, "title")
        row = layout.row()
        row.prop(self, "method")
        row = layout.row()
        row.prop(self, "sourceType")
        row = layout.row()
        row.prop(self, "programPath")
        row = layout.row()
        row.prop(self, "pictures")
        # post-processing
        layout.label(text="Post-processing:")
        row = layout.row()
        row.prop(self, "reference")
        row = layout.row()
        row.prop(self, "computeHRIRs")
        row = layout.row()
        # constants
        layout.label(text="Constants:")
        row = layout.row()
        row.prop(self, "unit")
        row = layout.row()
        row.prop(self, "speedOfSound")
        row = layout.row()
        row.prop(self, "densityOfMedium")
        row = layout.row()
        # evaluation grids
        layout.label(text="Evaluation Grids:")
        row = layout.row()
        row.prop(self, "evaluationGrids")
        # material search paths
        layout.label(text="Materials:")
        row = layout.row()
        row.prop(self, "materialSearchPaths")
        # frequency distribution
        layout.label(text="Frequency distribution:")
        row = layout.row()
        row.prop(self, "minFrequency")
        row = layout.row()
        row.prop(self, "maxFrequency")
        row = layout.row()
        row.prop(self, "frequencyVectorType")
        row = layout.row()
        row.prop(self, "frequencyVectorValue")

    def save(operator,
             context,
             filepath="",
             title="head-related transfer functions",
             minFrequency=100,
             maxFrequency=20000,
             frequencyVectorType='Step size',
             frequencyVectorValue=100,
             pictures=True,
             evaluationGrids='ARI',
             materialSearchPaths='None',
             method='ML-FMM BEM',
             reference=False,
             computeHRIRs=False,
             speedOfSound='346.18',
             densityOfMedium='1.1839',
             unit='mm',
             programPath="path/to/mesh2hrtf",
             sourceType='Both ears'
             ):
        """Export Mesh2HRTF project."""

# General handling and constants ----------------------------------------------
        # purge unused data
        ret = bpy.ops.outliner.orphans_purge()
        while ret != {'CANCELLED'}:
            ret = bpy.ops.outliner.orphans_purge()

        # Switch to object mode to avoid export errors
        bpy.ops.object.mode_set(mode='OBJECT', toggle=False)

        # de-select all objects
        bpy.ops.object.select_all(action='DESELECT')

        # check if 'Reference' object exists and assign to `obj`
        referenceExist = False
        for obj in bpy.context.scene.objects[:]:
            if obj.type == 'MESH' and obj.name == 'Reference':
                referenceExist = True
                break

        # select and activate 'Reference'
        if referenceExist:
            bpy.data.objects['Reference'].select_set(True)
            bpy.context.view_layer.objects.active = \
                bpy.data.objects['Reference']
        else:
            raise ValueError("Did not find the 3D Mesh. It must be named"
                             "'Reference' (case sensitive).")

        # check if 'Reference' object is a triangular mesh
        has_error_message = False

        for p in bpy.data.objects['Reference'].data.polygons:
            # Select non quad face (polygons)
            p.select = len(p.vertices) != 3
            if p.select:
                has_error_message = True

        if has_error_message:
            # Go in edit mode to show the result
            bpy.ops.object.mode_set(mode='EDIT')
            raise TypeError(
                ('Not all faces in the Reference mesh are triangular! '
                 'Non triangular faces are selected'))

        del has_error_message

        # get Mesh2HRTF version
        with open(os.path.join(programPath, "..", "VERSION")) as read_version:
            version = read_version.readline()

        # Export path and export directory handling
        if not filepath.endswith(os.path.sep):
            filepath += os.path.sep
        filepath1, _ = os.path.split(filepath)
        # create project folder
        if not os.path.exists(filepath1):
            os.mkdir(filepath1)
        # check for existing projects
        if os.path.exists(os.path.join(filepath1, "NumCalc")):
            raise ValueError((f"Project folder {filepath1} already exists. "
                              "Choose another folder or delete files."))
        # create sub-folders
        subfolders = ["ObjectMeshes", "EvaluationGrids", "NumCalc"]
        for subfolder in subfolders:
            temp = os.path.join(filepath1, subfolder)
            # delete subfolder if it exists to ensure data consistency
            if os.path.exists(temp):
                shutil.rmtree(temp)
            # make subfolder
            os.mkdir(temp)

        # number of sources used to create NC.inp files
        numSources = 1 if sourceType != "Both ears" else 2

        # check input and assign unitFactor
        if unit not in ["mm", "m"]:
            raise ValueError("`unit` must be 'mm', 'm' "
                             "(case sensitive).")
        unitFactor = 1 if unit == 'm' else .001

        # apply transforms
        bpy.ops.object.transform_apply(location=True)
        bpy.ops.object.transform_apply(rotation=True)
        bpy.ops.object.transform_apply(scale=True)

        # empty list for saving objects
        objects = ([])


# Get point source position ---------------------------------------------------
        if sourceType in ["Point source", "Plane wave"]:
            sourceXPosition, sourceYPosition, sourceZPosition = \
                _get_source_position(sourceType, unitFactor)
            numSources = len(sourceXPosition)
        else:
            sourceXPosition = None
            sourceYPosition = None
            sourceZPosition = None


# Sort faces in the object according to materials -----------------------------
        _sort_faces_according_to_materials(obj)


# Write object data -----------------------------------------------------------
        for obj in bpy.context.scene.objects[:]:
            if obj.type == 'MESH' and obj.name == 'Reference':
                objects = _write_object_data(obj, objects, unitFactor,
                                             context, filepath1)


# save Blender project for documentation --------------------------------------
        bpy.ops.wm.save_as_mainfile(
            filepath=os.path.join(filepath1, "3d Model.blend"),
            check_existing=False, filter_blender=True, filter_image=False,
            filter_movie=False, filter_python=False, filter_font=False,
            filter_sound=False, filter_text=False, filter_btx=False,
            filter_collada=False, filter_folder=True, filemode=8,
            compress=False, relative_remap=True, copy=False)


# Write evaluation grid data --------------------------------------------------
        # split and sort evaluation grids and get absolute paths
        evaluationGrids, evalGridPaths = \
            _split_and_sort_evaluation_grids(evaluationGrids, programPath)

        # check if the evaluation grids exists
        _check_evaluation_grid_exists(evaluationGrids, evalGridPaths)

        # copy data to export directory
        for n, _ in enumerate(evaluationGrids):

            # create target directory
            savepath = os.path.join(
                filepath1, "EvaluationGrids", evaluationGrids[n])
            if not os.path.exists(savepath):
                os.mkdir(savepath)

            # copy the data
            shutil.copyfile(os.path.join(evalGridPaths[n], "Nodes.txt"),
                            os.path.join(savepath, "Nodes.txt"))
            shutil.copyfile(os.path.join(evalGridPaths[n], "Elements.txt"),
                            os.path.join(savepath, "Elements.txt"))

# Read material data ----------------------------------------------------------

        # add the default mesh2hrtf path to the material search path
        defaultPath = os.path.join(
            programPath, 'Mesh2Input', 'Materials', 'Data')
        if materialSearchPaths == "None":
            materialSearchPaths = defaultPath
        else:
            materialSearchPaths += f";  {defaultPath}"

        # get used materials as dictionary
        materials = _get_materials(bpy.data.objects['Reference'])

        if materials is None:
            if 'ear' in sourceType:
                raise ValueError(
                    ("Material 'Left ear' and/or 'Right ear' "
                     "must be defined for reciprocal simulations."))
        else:
            # add full path of material files to the dictionary
            materials = _get_material_files(materials, materialSearchPaths)

            # read the material data and add it to the dictionary
            materials = _read_material_data(materials)


# Calculate frequency information ---------------------------------------------

        frequencies, frequencyStepSize, numFrequencySteps = \
            _calc_frequencies_from_input(
                frequencyVectorType, frequencyVectorValue, 
                minFrequency, maxFrequency)

# Write parameters.json (feedback for user, not used by NumCalc) --------------
        _write_parameters_json(
            filepath1, title, programPath, version, method, pictures,
            evaluationGrids, materialSearchPaths, materials,
            speedOfSound, densityOfMedium, unit, unitFactor,
            reference, computeHRIRs,
            sourceType, numSources, sourceXPosition,
            sourceYPosition, sourceZPosition,
            frequencies, frequencyStepSize, numFrequencySteps)

# Render pictures of the model ------------------------------------------------
        if pictures:
            _render_pictures(filepath1, unitFactor)


# Write NumCalc input files for all sources (NC.inp) --------------------------
        _write_nc_inp(filepath1, version, title, speedOfSound,
                      densityOfMedium, frequencies, evaluationGrids, materials,
                      method, sourceType, numSources,
                      sourceXPosition, sourceYPosition, sourceZPosition)

# Finish ----------------------------------------------------------------------
        for obj in bpy.context.scene.objects[:]:
            bpy.data.objects[obj.name].select_set(False)

        return {'FINISHED'}


def _get_source_position(sourceType, unitFactor):
    # check if 'Reference' object exists and select it
    sourceExist = False
    allSourceTypes = list()
    for obj in bpy.context.scene.objects[:]:
        if obj.type == 'LIGHT' and obj.name.startswith(sourceType):
            allSourceTypes.append(obj.name)
            sourceExist = True

    if not sourceExist:
        light_type = "a 'Point'" if sourceType == "Point source" \
            else "an 'Area'"
        raise ValueError((
            f"Did not find the {sourceType.lower()}. It must be {light_type} "
            f"light object named '{sourceType}' (case sensitive)."))

    # check for multiple sources
    sourceXPosition = []
    sourceYPosition = []
    sourceZPosition = []
    # allSourceTypes = allSourceTypes.sort()
    for source in allSourceTypes:
        sourceXPosition.append(
            bpy.data.objects[source].location[0] * unitFactor)
        sourceYPosition.append(
            bpy.data.objects[source].location[1] * unitFactor)
        sourceZPosition.append(
            bpy.data.objects[source].location[2] * unitFactor)

    if sourceType == "Plane wave":
        for idx in range(len(sourceXPosition)):
            abs_source_position = math.sqrt(
                sourceXPosition[idx]**2 + sourceYPosition[idx]**2
                + sourceZPosition[idx]**2)

            sourceXPosition[idx] /= abs_source_position
            sourceYPosition[idx] /= abs_source_position
            sourceZPosition[idx] /= abs_source_position

    return sourceXPosition, sourceYPosition, sourceZPosition



def _sort_faces_according_to_materials(obj):
    """
    Sort faces in an object according to the materials. This makes the NC.inp
    files shorter in case boundary conditions are used."""

    # enforce object mode and select the object
    bpy.ops.object.mode_set(mode='OBJECT', toggle=False)
    obj.select_set(True)

    # sort mesh if required
    bpy.ops.object.mode_set(mode='EDIT', toggle=False)
    bpy.ops.mesh.select_all(action='SELECT')
    bpy.ops.mesh.sort_elements(type='MATERIAL', elements={'FACE'})
    bpy.ops.object.mode_set(mode='OBJECT', toggle=False)


def _write_object_data(obj, objects, unitFactor, context, filepath1):
    """Write object information to Nodes.txt and Elements.txt.

    Returns:
    objects: list
        Original list with name of the current object appended.

    """
    # apply transformations and activate
    bpy.context.view_layer.objects.active = obj
    bpy.ops.object.transform_apply(location=True)
    bpy.ops.object.transform_apply(rotation=True)
    bpy.ops.object.transform_apply(scale=True)
    obj = context.active_object
    obj.hide_render = False
    obj_data = obj.data

    # create target directory
    dirObject = os.path.join(filepath1, "ObjectMeshes", obj.name)
    if not os.path.exists(dirObject):
        os.mkdir(dirObject)

    # write object nodes
    file = open(os.path.join(dirObject, "Nodes.txt"), "w",
                encoding="utf8", newline="\n")
    fw = file.write
    # number of nodes
    fw("%i\n" % len(obj_data.vertices[:]))
    # coordinate of nodes (vertices)
    for ii in range(len(obj_data.vertices[:])):
        fw("%i " % ii)
        fw("%.6f %.6f %.6f\n" %
            (obj_data.vertices[ii].co[0] * unitFactor,
             obj_data.vertices[ii].co[1] * unitFactor,
             obj_data.vertices[ii].co[2] * unitFactor))
    file.close

    # write object elements
    file = open(os.path.join(dirObject, "Elements.txt"), "w",
                encoding="utf8", newline="\n")
    fw = file.write
    # number of faces (polygons)
    fw("%i\n" % len(obj_data.polygons[:]))
    # node (vertice) indicees
    for ii in range(len(obj_data.polygons[:])):
        fw("%i " % ii)
        if len(obj_data.polygons[ii].vertices[:]) == 3:
            fw("%d %d %d 0 0 0\n" %
               tuple(obj_data.polygons[ii].vertices[:]))
        else:
            fw("%d %d %d %d 0 0 0\n" %
               tuple(obj_data.polygons[ii].vertices[:]))
    file.close

    # add object to list
    objects.append(obj.name)

    return objects


def _split_and_sort_evaluation_grids(evaluationGrids, programPath):
    """Split and sort evaluation grids.

    Returns
    -------
    evaluationGrids: list
        Alphabetically sorted list of the evaluation grid names
        (case sensitive).
    evaluationGridPaths: list
        Full paths of the folders containing the evaluation grids in order of
        `evaluationGrids`.

    """
    if not len(evaluationGrids):
        raise ValueError("At least one evaluation grid mast be defines.")

    # split evaluationGrids and remove leading/trailing white spaces
    grids = evaluationGrids.split(';')
    grids = [eg.strip() for eg in grids]

    # default path for evaluation grids
    evaluationGridPath = os.path.join(
        programPath, "Mesh2Input", "EvaluationGrids", "Data")

    # construct full paths and names
    evalGridPaths = []
    evalGrids = []
    for eg in grids:
        # complete grids given by name to default path
        if os.path.sep not in eg:
            evalGridPaths.append(os.path.join(evaluationGridPath, eg))
            evalGrids.append(eg)
        else:
            evalGridPaths.append(eg)
            evalGrids.append(os.path.basename(os.path.normpath(eg)))

    # alphabetically sort according to evalGrids
    idx = sorted(range(len(evalGrids)), key=evalGrids.__getitem__)
    unsrt = evalGrids
    evaluationGrids = [unsrt[i] for i in idx]
    unsrt = evalGridPaths
    evalGridPaths = [unsrt[i] for i in idx]

    return evaluationGrids, evalGridPaths


def _check_evaluation_grid_exists(evaluationGrids, evalGridPaths):
    """Check if evaluation grid exists.

    Raises
    ------
    ValueError if the folder or neccessary files do not exist.
    """
    for n, _ in enumerate(evaluationGrids):
        if not os.path.isfile(
                os.path.join(evalGridPaths[n], 'Nodes.txt')) \
                or not \
                os.path.isfile(
                os.path.join(evalGridPaths[n], 'Elements.txt')):
            raise ValueError(
                f"Evalution grid folder {evalGridPaths[n]} does not exist "
                "or one of the files 'Nodes.txt' and 'Elements.txt' is "
                "missing.")


def _get_materials(obj):

    # enforce object mode and select the object
    bpy.ops.object.mode_set(mode='OBJECT', toggle=False)
    obj.select_set(True)

    # loop faces and log assigned materials
    materials = dict()

    for face in obj.data.polygons:
        # name of material of current face
        try:
            material = obj.material_slots[face.material_index].material.name
        except IndexError:
            return None

        # add material if it does not exist
        if material not in materials.keys():
            materials[material] = dict()
            materials[material]["is_used"] = True
            materials[material]["index_start"] = face.index
            materials[material]["index_end"] = face.index
        # update index_end if it does exist
        else:
            materials[material]["index_end"] = face.index

    # add default materials if missing
    for material in ["Skin", "Left ear", "Right ear"]:
        if material not in materials.keys():
            materials[material] = dict()
            materials[material]["is_used"] = False

    return materials


def _get_material_files(materials, materialSearchPaths):
    """Check if material files exists and return full paths as a list."""

    # split search paths and remove leading/trailing white spaces
    paths = materialSearchPaths.split(';')
    paths = [path.strip() for path in paths]

    # check if paths exist
    for path in paths:
        if not os.path.isdir(path):
            raise ValueError(f"Material search path '{path}' does not exist.")

    # search material in path
    for material in materials:
        materialFound = False
        for path in paths:
            if os.path.isfile(os.path.join(path, f"{material}.csv")):
                materialFound = True
                break

        # save full path of material
        if materialFound:
            materials[material]["path"] = os.path.join(path, f"{material}.csv")
        else:
            materials[material]["path"] = None
        # throw error if the material is not found
        if not materialFound \
           and material not in ['Skin', 'Left ear', 'Right ear']:
            raise ValueError(f"Material file '{material}.csv' not found.")

    return materials


def _read_material_data(materials):

    for material in materials:
        # current material file
        file = materials[material]["path"]
        # check if the file exists
        if file is None:
            continue

        # initilize data
        boundary = None
        freqs = []
        real = []
        imag = []

        # read the csv material file
        with open(file, 'r') as m:
            lines = m.readlines()

        # parse the file
        for line in lines:
            line = line.strip('\n')
            # skip empty lines and comments
            if not len(line):
                continue
            if line[0] == '#':
                continue

            # detect boundary keyword
            if line in ['ADMI', 'IMPE', 'VELO', 'PRES']:
                boundary = line
            # read curve value
            else:
                line = line.split(',')
                if not len(line) == 3:
                    raise ValueError(
                        (f'Expected three values in {file} '
                         f'definition but found {len(line)}'))
                freqs.append(line[0].strip())
                real.append(line[1].strip())
                imag.append(line[2].strip())

        # check if boundary keyword was found
        if boundary is None:
            raise ValueError(
                (f"No boundary definition found in {file}. "
                 "Must be 'ADMI', 'IMPE', 'VELO', or 'PRES'"))
        # check if frequency vector is valud
        for i in range(len(freqs)-1):
            if float(freqs[i+1]) <= float(freqs[i]):
                raise ValueError((f'Frequencies in {file} '
                                  'do not increase monotonously'))

        # create output
        materials[material]['boundary'] = boundary
        materials[material]['freqs'] = freqs
        materials[material]['real'] = real
        materials[material]['imag'] = imag

    return materials


def _write_parameters_json(
        filepath1, title, programPath, version, method, pictures,
        evaluationGrids, materialSearchPaths, materials,
        speedOfSound, densityOfMedium, unit, unitFactor,
        reference, computeHRIRs,
        sourceType, numSources, sourceXPosition,
        sourceYPosition, sourceZPosition,
        frequencies, frequencyStepSize, numFrequencySteps):

    # calculate missing parameters
    if "ear" in sourceType:
        # get the receiver/ear centers and areas
        obj = bpy.data.objects['Reference']
        obj_data = obj.data
        sourceCenter, sourceArea = _calculateReceiverProperties(
            obj, obj_data, unitFactor)

        # select ear specific data
        if sourceType == "Left ear":
            sourceCenter = sourceCenter[0]
            sourceArea = [sourceArea[0]]
        if sourceType == "Right ear":
            sourceCenter = sourceCenter[1]
            sourceArea = [sourceArea[1]]

    else:
        sourceCenter = np.array([sourceXPosition, sourceYPosition, sourceZPosition])
        sourceCenter = list(np.transpose(sourceCenter))
        sourceCenter = [list(x) for x in sourceCenter]
        sourceArea = [1]

    # write parameters to dict
    parameters = {
        # project Info
        "projectTitle": title,
        "Mesh2HRTF_Path": programPath,
        "Mesh2HRTF_Version": version,
        "BEM_Type": method,
        "exportPictures": pictures,
        # Constants
        "speedOfSound": float(speedOfSound),
        "densityOfMedium": float(densityOfMedium),
        "3D_SceneUnit": unit,
        # Grids and materials
        "evaluationGrids": evaluationGrids,
        "materialSearchPaths": materialSearchPaths,
        "materials": materials,
        # Source definition
        "sourceType": sourceType,
        "numSources": numSources,
        "sourceCenter": sourceCenter,
        "sourceArea": sourceArea,
        # post processing
        "reference": reference,
        "computeHRIRs": computeHRIRs,
        # frequencies
        "numFrequencies": numFrequencySteps,
        "frequencyStepSize": frequencyStepSize,
        "minFrequency": frequencies[0],
        "maxFrequency": frequencies[-1],
        "frequencies": frequencies
    }

    with open(os.path.join(filepath1, "parameters.json"), 'w') as file:
        json.dump(parameters, file, indent=4)


def _calculateReceiverProperties(obj, obj_data, unitFactor):
    """
    Calculate center and area of receiver elements in reciprocal simulation.

    """

    # init ear area
    earArea = [0., 0.]
    # init min and max xyz coordinates for left and right ear
    earCenter = [[[1e6, -1e6], [1e6, -1e6], [1e6, -1e6]],
                 [[1e6, -1e6], [1e6, -1e6], [1e6, -1e6]]]

    # loop elements in obj_data
    for ii in range(len(obj_data.polygons[:])):
        if obj.material_slots[obj_data.polygons[ii].material_index].name \
                == obj.material_slots['Left ear'].name \
                or \
                obj.material_slots[obj_data.polygons[ii].material_index].name \
                == obj.material_slots['Right ear'].name:

            # select the ear
            if obj.material_slots[obj_data.polygons[ii].material_index].name \
                    == obj.material_slots['Left ear'].name:
                ear = 0
            else:
                ear = 1

            # update min and max x,y,z-values
            for vertex in range(3):
                for coord in range(3):
                    value = obj_data.vertices[obj_data.polygons[ii].
                                              vertices[vertex]].co[coord]
                    if value < earCenter[ear][coord][0]:
                        earCenter[ear][coord][0] = value
                    if value > earCenter[ear][coord][1]:
                        earCenter[ear][coord][1] = value

            # to calculate the polygons (triangles) area first calculate the
            # side lengths using euclidean distance and then use Heron´s
            # formula to calculate the area
            # side_a = corner 0 to corner 1
            side_a = math.sqrt(((
                obj_data.vertices[obj_data.polygons[ii].vertices[0]].co[0] -
                obj_data.vertices[obj_data.polygons[ii].vertices[1]].co[0])**2)
                * unitFactor**2
                + ((
                obj_data.vertices[obj_data.polygons[ii].vertices[0]].co[1] -
                obj_data.vertices[obj_data.polygons[ii].vertices[1]].co[1])**2)
                * unitFactor**2
                + ((
                obj_data.vertices[obj_data.polygons[ii].vertices[0]].co[2] -
                obj_data.vertices[obj_data.polygons[ii].vertices[1]].co[2])**2)
                * unitFactor**2)

            # side_b = corner 1 to corner 2
            side_b = math.sqrt(((
                obj_data.vertices[obj_data.polygons[ii].vertices[1]].co[0] -
                obj_data.vertices[obj_data.polygons[ii].vertices[2]].co[0])**2)
                * unitFactor**2
                + ((
                obj_data.vertices[obj_data.polygons[ii].vertices[1]].co[1] -
                obj_data.vertices[obj_data.polygons[ii].vertices[2]].co[1])**2)
                * unitFactor**2
                + ((
                obj_data.vertices[obj_data.polygons[ii].vertices[1]].co[2] -
                obj_data.vertices[obj_data.polygons[ii].vertices[2]].co[2])**2)
                * unitFactor**2)

            # side_c = corner 2 to corner 0
            side_c = math.sqrt(((
                obj_data.vertices[obj_data.polygons[ii].vertices[2]].co[0] -
                obj_data.vertices[obj_data.polygons[ii].vertices[0]].co[0])**2)
                * unitFactor**2
                + ((
                obj_data.vertices[obj_data.polygons[ii].vertices[2]].co[1] -
                obj_data.vertices[obj_data.polygons[ii].vertices[0]].co[1])**2)
                * unitFactor**2
                + ((
                obj_data.vertices[obj_data.polygons[ii].vertices[2]].co[2] -
                obj_data.vertices[obj_data.polygons[ii].vertices[0]].co[2])**2)
                *unitFactor**2)

            # increment area using Heron´s formula
            earArea[ear] += 0.25 * math.sqrt(
                             (side_a+side_b+side_c)
                             * (-side_a+side_b+side_c)
                             * (side_a-side_b+side_c)
                             * (side_a+side_b-side_c))

    # estimate the center from min and max x,y,z-values
    earCenter = [
        [  # left ear center
         (earCenter[0][0][0] + earCenter[0][0][1]) / 2 * unitFactor,
         (earCenter[0][1][0] + earCenter[0][1][1]) / 2 * unitFactor,
         (earCenter[0][2][0] + earCenter[0][2][1]) / 2 * unitFactor],
        [  # right ear center
         (earCenter[1][0][0] + earCenter[1][0][1]) / 2 * unitFactor,
         (earCenter[1][1][0] + earCenter[1][1][1]) / 2 * unitFactor,
         (earCenter[1][2][0] + earCenter[1][2][1]) / 2 * unitFactor]]

    return earCenter, earArea


def _distribute_frequencies(minFrequency, maxFrequency,
                            frequencyStepSize, numFrequencySteps):
    """Calculate list of frequencies.

    Returns
    -------

    frequencies: list
        list that holds the frequencies in Hz that are calculated.
    frequencyStepSize: float
        Step size between successive frequencies (written to parameters.json).
        This is returned because it might be zero during the function call.
    numFrequencySteps: int
        Number of frequencies to be simulated (written to parameters.json).
        This is returned because it might be zero during the function call.

    """

    # check input
    if (numFrequencySteps == 0 and frequencyStepSize == 0) \
            or (numFrequencySteps != 0 and frequencyStepSize != 0):
        raise Exception(("Either 'frequencyStepSize' or 'numFrequencySteps' "
                         "must be zero while the other must not."))

    # Calculate Number of frequencies and frequency step size
    if minFrequency == maxFrequency:
        frequencySteps = (1, 0)
        frequencyStepSize = 0
    elif frequencyStepSize:
        frequencySteps = divmod(
            maxFrequency - minFrequency + frequencyStepSize, frequencyStepSize)
    else:
        if numFrequencySteps < 1:
            raise Exception("'numFrequencySteps' must be at least 1.")
        frequencySteps = (numFrequencySteps, 0)
        frequencyStepSize = (maxFrequency-minFrequency)/(numFrequencySteps-1)

    if not frequencySteps[1] == 0:
        raise Exception(("Error, frequencyStepSize is not a divisor of "
                         "maxFrequency-minFrequency"))

    # get all frequencies to be calculated
    frequencies = [ff*frequencyStepSize+minFrequency
                   for ff in range(int(frequencySteps[0]))]

    # remove 0 Hz if included in the list
    if frequencies[0] == 0:
        frequencies.pop(0)
        frequencySteps = (frequencySteps[0]-1, 0)
        print(('Warning: 0 Hz can not be calculated and was removed from the '
               'list of frequencies.'))

    numFrequencySteps = len(frequencies)

    if numFrequencySteps == 0:
        raise ValueError(("No frequencies to be calculated. "
                          "Check the input parameters."))
    elif numFrequencySteps == 1:
        frequencyStepSize = 0
    else:
        frequencyStepSize = frequencies[1] - frequencies[0]

    return frequencies, frequencyStepSize, numFrequencySteps


def _render_pictures(filepath1, unitFactor):
    """Render pictures of the 3D mesh and save to export folder."""

    # create directoy
    dirRender = os.path.join(filepath1, "Pictures")
    if not os.path.exists(dirRender):
        os.mkdir(dirRender)

    # general camera settings
    cam = bpy.data.objects['Camera']
    camdistance = .6 / unitFactor
    cam.data.clip_start = 0.001 / unitFactor
    cam.data.clip_end = 1 / unitFactor
    bpy.data.cameras['Camera'].lens_unit = 'MILLIMETERS'
    bpy.data.cameras['Camera'].lens = 50
    bpy.data.cameras['Camera'].clip_start = .01 / unitFactor
    bpy.data.cameras['Camera'].clip_end = 1 / unitFactor

    # general light settings
    light = bpy.data.objects['Light']
    lightradius = 5 / unitFactor
    bpy.data.lights['Light'].use_custom_distance = False
    bpy.data.lights['Light'].energy = 800 / unitFactor ** 2
    bpy.data.lights['Light'].shadow_soft_size = .1 / unitFactor
    bpy.data.lights['Light'].cutoff_distance = 10 / unitFactor
    bpy.data.lights['Light'].shadow_buffer_clip_start = .05 / unitFactor
    bpy.data.lights['Light'].shadow_buffer_bias = 1 / unitFactor

    # camera positions, rotations, and azimuth angles
    renderloc = [[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0],
                 [0.707, 0.707, 0], [-0.707, 0.707, 0], [0.707, -0.707, 0],
                 [-0.707, -0.707, 0]]
    renderrot = [[pi/2, 0, pi/2], [pi/2, 0, 3*pi/2], [pi/2, 0, pi],
                 [3/2*pi, pi, pi], [pi/2, 0, 3/4*pi], [pi/2, 0, 5/4*pi],
                 [pi/2, 0, pi/4], [pi/2, 0, -pi/4]]
    azim = [0, 180, 90, 270, 45, 135, 315, 225]

    # general rendering settings
    bpy.data.scenes['Scene'].render.pixel_aspect_x = 1
    bpy.data.scenes['Scene'].render.pixel_aspect_y = 1
    bpy.data.scenes['Scene'].render.resolution_x = 1440
    bpy.data.scenes['Scene'].render.resolution_y = 1920

    # render the images
    for ii in range(len(renderloc)):
        # set the camera and light
        cam.location = (renderloc[ii][0] * camdistance,
                        renderloc[ii][1] * camdistance,
                        renderloc[ii][2] * camdistance)
        cam.rotation_euler = (renderrot[ii][0],
                              renderrot[ii][1],
                              renderrot[ii][2])
        light.location = (renderloc[ii][0] * lightradius,
                          renderloc[ii][1] * lightradius,
                          renderloc[ii][2] * lightradius)
        # render
        bpy.ops.render.render()
        # save
        nameRender = ("%d_deg_azimuth" % azim[ii])
        bpy.data.images['Render Result'].save_render(
            os.path.join(dirRender, "%s.png" % nameRender))


def _write_nc_inp(filepath1, version, title,
                  speedOfSound, densityOfMedium, frequencies,
                  evaluationGrids, materials, method, sourceType, numSources,
                  sourceXPosition, sourceYPosition, sourceZPosition):
    """Write NC.inp file that is read by NumCalc to start the simulation.

    The file format is documented at:
    https://github.com/Any2HRTF/Mesh2HRTF/wiki/Structure_of_NC.inp
    """

    # check the BEM method
    if method == 'BEM':
        method_id = 0
    elif method == 'SL-FMM BEM':
        method_id = 1
    elif method == 'ML-FMM BEM':
        method_id = 4
    else:
        ValueError(
            f"Method must be BEM, SL-FMM BEM or ML-FMM BEM but is {method}")

    for source in range(numSources):

        # create directory
        filepath2 = os.path.join(
            filepath1, "NumCalc", f"source_{source+1}")
        if not os.path.exists(filepath2):
            os.mkdir(filepath2)

        # write NC.inp
        file = open(os.path.join(filepath2, "NC.inp"), "w",
                    encoding="utf8", newline="\n")
        fw = file.write

        obj_name = "Reference"

        obj = bpy.data.objects[obj_name]
        obj_data = obj.data

        # header --------------------------------------------------------------
        fw("##-------------------------------------------\n")
        fw("## This file was created by mesh2input\n")
        fw("## Date: %s\n" % datetime.date.today())
        fw("##-------------------------------------------\n")
        fw("Mesh2HRTF %s\n" % version)
        fw("##\n")
        fw("%s\n" % title)
        fw("##\n")

        # control parameter I (hard coded, not documented) --------------------
        fw("## Controlparameter I\n")
        fw("0 0 0 0 7 0\n")
        fw("##\n")

        # control parameter II ------------------------------------------------
        fw("## Controlparameter II\n")
        fw("1 %d 0.000001 0.00e+00 1 0 0\n" % (
            len(frequencies)))
        fw("##\n")
        fw("## Load Frequency Curve \n")
        fw("0 %d\n" % (len(frequencies)+1))
        fw("0.000000 0.000000e+00 0.0\n")
        for ii in range(len(frequencies)):
            fw("%f %fe+04 0.0\n" % (
                0.000001*(ii+1),
                frequencies[ii] / 10000))
        fw("##\n")

        # main parameters I ---------------------------------------------------
        fw("## 1. Main Parameters I\n")
        numNodes = 0
        numElements = 0
        for evaluationGrid in evaluationGrids:
            # read number of nodes
            nodes = open(os.path.join(
                filepath1, "EvaluationGrids", evaluationGrid,
                "Nodes.txt"))
            line = nodes.readline()
            numNodes = numNodes+int(line)
            # read number of elements
            elements = open(os.path.join(
                filepath1, "EvaluationGrids", evaluationGrid,
                "Elements.txt"))
            line = elements.readline()
            numElements = numElements+int(line)
        fw("2 %d " % (len(obj_data.polygons[:])+numElements))
        fw("%d 0 " % (len(obj_data.vertices[:])+numNodes))
        fw("0")
        fw(" 2 1 %s 0\n" % (method_id))
        fw("##\n")

        # main parameters II --------------------------------------------------
        fw("## 2. Main Parameters II\n")
        fw("0 ")
        if "ear" in sourceType:
            fw("0 ")
        else:
            fw("1 ")
        fw("0 0.0000e+00 0 0 0\n")
        fw("##\n")

        # main parameters III -------------------------------------------------
        fw("## 3. Main Parameters III\n")
        fw("0 0 0 0\n")
        fw("##\n")

        # main parameters IV --------------------------------------------------
        fw("## 4. Main Parameters IV\n")
        fw("%s %se+00 1.0 0.0e+00 0.0 e+00 0.0e+00 0.0e+00\n" % (
            speedOfSound, densityOfMedium))
        fw("##\n")

        # nodes ---------------------------------------------------------------
        fw("NODES\n")
        fw("../../ObjectMeshes/Reference/Nodes.txt\n")
        # write file path of nodes to input file
        for grid in evaluationGrids:
            fw("../../EvaluationGrids/%s/Nodes.txt\n" % grid)
        fw("##\n")
        fw("ELEMENTS\n")
        fw("../../ObjectMeshes/Reference/Elements.txt\n")
        # write file path of elements to input file
        for grid in evaluationGrids:
            fw("../../EvaluationGrids/%s/Elements.txt\n" % grid)
        fw("##\n")

        # SYMMETRY ------------------------------------------------------------
        fw("# SYMMETRY\n")
        fw("# 0 0 0\n")
        fw("# 0.0000e+00 0.0000e+00 0.0000e+00\n")
        fw("##\n")

        # assign mesh elements to boundary conditions -------------------------
        # (including both, left, right ear)
        fw("BOUNDARY\n")
        # write velocity condition for the ears if using vibrating
        # elements as the sound source
        if "ear" in sourceType:
            if source == 0 and \
                    sourceType in ['Both ears', 'Left ear']:
                tmpEar = 'Left ear'
            else:
                tmpEar = 'Right ear'
            fw(f"# {tmpEar} velocity source\n")
            fw("ELEM %i TO %i VELO 0.1 -1 0.0 -1\n" % (
                materials[tmpEar]["index_start"],
                materials[tmpEar]["index_end"]))
        # remaining conditions defined by frequency curves
        curves = 0
        steps = 0
        if materials is not None:
            for m in materials:
                if materials[m]["path"] is None:
                    continue
                # write information
                fw(f"# Material: {m}\n")
                fw("ELEM %i TO %i %s 1.0 %i 1.0 %i\n" % (
                    materials[m]["index_start"],
                    materials[m]["index_end"],
                    materials[m]["boundary"],
                    curves + 1, curves + 2))
                # update metadata
                steps = max(steps, len(materials[m]["freqs"]))
                curves += 2

        fw("RETU\n")
        fw("##\n")

        # source information: point source and plane wave ---------------------
        if sourceType == "Point source":
            fw("POINT SOURCES\n")
        elif sourceType == "Plane wave":
            fw("PLANE WAVES\n")
        if sourceType in ["Point source", "Plane wave"]:
            fw("0 %s %s %s 0.1 -1 0.0 -1\n" % (
                sourceXPosition[source], sourceYPosition[source], 
                sourceZPosition[source]))
        fw("##\n")

        # curves defining boundary conditions of the mesh ---------------------
        if curves > 0:
            fw("CURVES\n")
            # number of curves and maximum number of steps
            fw(f"{curves} {steps}\n")
            curves = 0
            for m in materials:
                if materials[m]["path"] is None:
                    continue
                # write curve for real values
                curves += 1
                fw(f"{curves} {len(materials[m]['freqs'])}\n")
                for f, v in zip(materials[m]['freqs'],
                                materials[m]['real']):
                    fw(f"{f} {v} 0.0\n")
                # write curve for imaginary values
                curves += 1
                fw(f"{curves} {len(materials[m]['freqs'])}\n")
                for f, v in zip(materials[m]['freqs'],
                                materials[m]['imag']):
                    fw(f"{f} {v} 0.0\n")

        else:
            fw("# CURVES\n")
        fw("##\n")

        # post process --------------------------------------------------------
        fw("POST PROCESS\n")
        fw("##\n")
        fw("END\n")
        file.close()


# ----------------------- Blender add-on registration -------------------------
def menu_func_export(self, context):
    self.layout.operator(ExportMesh2HRTF.bl_idname, text="Mesh2HRTF")


def register():
    bpy.utils.register_class(ExportMesh2HRTF)
    bpy.types.TOPBAR_MT_file_export.append(menu_func_export)


def unregister():
    bpy.utils.unregister_class(ExportMesh2HRTF)
    bpy.types.TOPBAR_MT_file_export.remove(menu_func_export)


def _calc_frequencies_from_input(
        frequencyVectorType, frequencyVectorValue, minFrequency, maxFrequency):
    # check how the frequency vector is defined
    if frequencyVectorType == 'Step size':
        frequencyStepSize = frequencyVectorValue
        numFrequencySteps = 0
    elif frequencyVectorType == 'Num steps':
        frequencyStepSize = 0
        numFrequencySteps = int(frequencyVectorValue)

    if frequencyVectorType in ('Step size', 'Num steps'):
        frequencies, frequencyStepSize, numFrequencySteps = \
            _distribute_frequencies(minFrequency, maxFrequency,
                                    frequencyStepSize, numFrequencySteps)
    if frequencyVectorType in ('Nominal n-th octave', 'Exact n-th octave'):
        frequencyStepSize = 0
        if frequencyVectorType == 'Nominal n-th octave':
            if int(frequencyVectorValue) == 1:
                all_frequencies = [
                    31.5, 63, 125, 250, 500, 1000,
                    2000, 4000, 8000, 16000, 32000,
                    64000, 128000]
            elif int(frequencyVectorValue) == 3:
                all_frequencies = [
                    25, 31.5, 40, 50, 63, 80, 100, 125, 160,
                    200, 250, 315, 400, 500, 630, 800, 1000,
                    1250, 1600, 2000, 2500, 3150, 4000, 5000,
                    6300, 8000, 10000, 12500, 16000, 20000,
                    25000, 31500, 40000, 50000, 63000, 80000,
                    100000, 125000, 160000, 200000]
            else:
                raise ValueError(
                    'Just 1st and 3rd ocatve are allowed in nominal, use exact'
                    ' instead.')
        elif frequencyVectorType == 'Exact n-th octave':
            if frequencyVectorValue < 0:
                raise ValueError('order must be larger than 0.')
            all_frequencies = 10**3 * np.power(
                2, np.arange(-18., 20., 1.)/frequencyVectorValue)
        frequencies = []
        for f in all_frequencies:
            if f >= minFrequency and f <= maxFrequency:
                frequencies.append(f)
        numFrequencySteps = len(frequencies)
    return frequencies, frequencyStepSize, numFrequencySteps
