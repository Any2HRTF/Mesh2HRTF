# Authors: The Mesh2HRTF developers
#
#                                Mesh2HRTF
#                Copyright (C) 2015 by Harald Ziegelwanger,
#        Acoustics Research Institute, Austrian Academy of Sciences
#                        mesh2hrtf.sourceforge.net
#
# Mesh2HRTF is licensed under the GNU Lesser General Public License as
# published by the Free Software Foundation, either version 3 of the License,
# or (at your option) any later version. Mesh2HRTF is distributed in the hope
# that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details. You should have received a
# copy of the GNU LesserGeneral Public License along with Mesh2HRTF. If not,
# see <http://www.gnu.org/licenses/lgpl.html>.
#
# If you use Mesh2HRTF:
# - Provide credits:
#   "Mesh2HRTF, H. Ziegelwanger, ARI, OEAW (mesh2hrtf.sourceforge.net)"
# - In your publication, cite both articles:
#   [1] Ziegelwanger, H., Kreuzer, W., and Majdak, P. (2015). "Mesh2HRTF:
#       Open-source software package for the numerical calculation of head-
#       related transfer functions," in Proceedings of the 22nd ICSV,
#       Florence, IT.
#   [2] Ziegelwanger, H., Majdak, P., and Kreuzer, W. (2015). "Numerical
#       calculation of listener-specific head-related transfer functions and
#       sound localization: Microphone model and mesh discretization," The
#       Journal of the Acoustical Society of America, 138, 208-222.

import os
import bpy
import bmesh
import datetime
import math
import shutil
from math import pi
from bpy.props import StringProperty, BoolProperty, EnumProperty, \
    IntProperty, FloatProperty
from bpy_extras.io_utils import ExportHelper

bl_info = {
    "name": "Mesh2HRTF input format",
    "author": "The Mesh2HRTF developers",
    "version": (1, 0, 0),
    "blender": (2, 80, 0),
    "location": "File > Export",
    "description": "Export Mesh2HRTF input files",
    "warning": "",
    "wiki_url": "",
    "tracker_url": "",
    "support": "COMMUNITY",
    "category": "Import-Export"}


class ExportMesh2HRTF(bpy.types.Operator, ExportHelper):
    '''Export an object as Mesh2HRTF input files'''
    bl_idname = "export_mesh2hrtf.inp"
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
        items=[('Vibrating element', 'Vibrating element',
                    ("Mesh elements with user assigned materials 'Left ear' "
                     "and 'Right ear' act as the source")),
               ('Point source', 'Point source',
                    ("Analytical point source. Coordinates taken from user "
                     "placed point light named 'Point'"))],
        default='Vibrating element',
        )
    ear: EnumProperty(
        name="Ear",
        description=("Selected ear(s) for simulation. Only required if using "
                     "vibrating elements as the source type"),
        items=[('Left ear', 'Left ear', 'Calculate left ear HRTFs only'),
               ('Right ear', 'Right ear', 'Calculate right ear HRTFs only'),
               ('Both ears', 'Both ears', 'Calculate HRTFs for both ears')],
        default='Both ears',
        )
    programPath: StringProperty(
        name="Mesh2HRTF-path",
        description=("Path to folder containing 'Mesh2Input' and other folders"
                     "(used to copy files to project folder during export)"),
        default=r"/home/matheson/Apps/mesh2hrtf-git/mesh2hrtf",
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
        description=("Name of evalation grid inside "
            "Mesh2Input/EvaluationsGrids or absolute path to user grid. "
            "Multiple grids can be separated by semicolons (;)"),
        default='ARI',
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
        max=24000,
        )
    frequencyVectorType: EnumProperty(
        name="Frequencies",
        description=("Select how the frequencies are distributed between min. "
                     "and max. frequency with the paramter 'Value' (below)."),
        items=[('Step size', 'Step size',
                    ("Simulate frequencies between the min. and max. "
                     "with a fixed step size.")),
               ('Num steps', 'Num steps',
                    "Simulate N frequencies the min. and max. frequency.")],
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
    # Job distribution --------------------------------------------------------
    cpuFirst: IntProperty(
        name="CPU (first)",
        description=("The simulation can be distributed to separate instances. "
                     "A separate folder is created created For each CPU_x_Core_x. "
                     "(This setting does NOT need to match your actual CPU and "
                     "some users might only care about the total nr. of instances)."),
        default=1,
        min=1,
        max=10000,
        )
    cpuLast: IntProperty(
        name="CPU (last)",
        description=("The simulation can be distributed to separate instances. "
                     "A separate folder is created created For each CPU_x_Core_x. "
                     "(This setting does NOT need to match your actual CPU and "
                     "some users might only care about the total nr. of instances)."),
        default=1,
        min=1,
        max=10000,
        )
    numCoresPerCPU: IntProperty(
        name="Cores per CPU",
        description=("The simulation can be distributed to separate instances. "
                     "A separate folder is created created For each CPU_x_Core_x. "
                     "(This setting does NOT need to match your actual CPU and "
                     "some users might only care about the total nr. of instances)."),
        default=1,
        min=1,
        max=10000,
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
        row.prop(self, "ear")
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
        # job distribution
        layout.label(text="Job distribution:")
        row = layout.row()
        row.prop(self, "cpuFirst")
        row = layout.row()
        row.prop(self, "cpuLast")
        row = layout.row()
        row.prop(self, "numCoresPerCPU")

    def save(operator,
             context,
             filepath="",
             title="head-related transfer functions",
             minFrequency=100,
             maxFrequency=20000,
             frequencyVectorType='Step size',
             frequencyVectorValue=100,
             cpuFirst=1,
             cpuLast=1,
             numCoresPerCPU=1,
             pictures=True,
             ear='Both ears',
             evaluationGrids='ARI',
             materialSearchPaths='None',
             method='ML-FMM BEM',
             reference=False,
             computeHRIRs=False,
             speedOfSound='346.18',
             densityOfMedium='1.1839',
             unit='mm',
             programPath="/home/matheson/Apps/mesh2hrtf-git/mesh2hrtf",
             sourceType='Vibrating element'
             ):
        """Export Mesh2HRTF project."""


# General handling and constants ----------------------------------------------
        # purge unused data
        ret = bpy.ops.outliner.orphans_purge()
        while ret != {'CANCELLED'}:
            ret = bpy.ops.outliner.orphans_purge()

        # Switch to object mode to avoid export errors
        bpy.ops.object.mode_set(mode='OBJECT', toggle=False)

        # check if 'Reference' object exists and select it
        referenceExist = False
        for obj in bpy.context.scene.objects[:]:
            if obj.type == 'MESH' and obj.name == 'Reference':
                referenceExist = True

        if referenceExist:
            bpy.data.objects['Reference'].select_set(True)
        else:
            raise ValueError("Did not find the 3D Mesh. It must be named"
                             "'Reference' (case sensitive).")

        # check if 'Reference' object is a triangular mesh
        obj2 = bpy.data.objects['Reference']
        has_error_message = False

        for p in obj2.data.polygons:
            # Select non quad face (polygons)
            p.select = len(p.vertices) != 3
            if p.select:
                has_error_message = True

        if has_error_message:
            # Go in edit mode to show the result
            bpy.ops.object.mode_set(mode = 'EDIT')
            raise TypeError(
                'Not all faces in the Reference mesh are triangular!')

        del obj2, has_error_message

        # get Mesh2HRTF version
        with open(os.path.join(programPath, "..", "VERSION")) as read_version:
            version = read_version.readline()

        # Export path and export directory handling
        if not filepath.endswith(os.path.sep):
            filepath += os.path.sep
        filepath1, _ = os.path.split(filepath)

        if not os.path.exists(filepath1): # create project folder with "Filename"
            os.mkdir(filepath1)
        elif os.path.exists( os.path.join(filepath1, "NumCalc") ):
            raise ValueError("Project folder ", filepath1,
                             " already exists. Please choose another Project "
                             "name or manually delete existing files.")

        subfolders = ["ObjectMeshes", "EvaluationGrids", "NumCalc"]
        for subfolder in subfolders:
            temp = os.path.join(filepath1, subfolder)
            # delete subfolder if it exists to ensure data consistency
            if os.path.exists(temp):
                shutil.rmtree(temp)
            # make subfolder
            os.mkdir(temp)

        # convert to int
        if sourceType == 'Vibrating element':
            sourceType_id = int(0)
        elif sourceType == 'Point source':
            sourceType_id = int(1)
        else:
            ValueError(("Source type must be 'Vibrating element' or "
                        f"'Point source' but is {sourceType}"))

        # assign numSources is used during assigning frequencies to cores. If
        # the simulation does not use a velocity source, numSources must be one
        # otherwise it is determined by the ears to be simulated
        if ear not in ["Left ear", "Right ear", "Both ears"]:
            raise ValueError("`ear` must be 'Left ear', 'Right ear' or "
                             "'Both ears' (case sensitive).")

        if sourceType_id != 0:
            numSources = 1
        else:
            numSources = 2 if ear == 'Both ears' else 1

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
        if sourceType_id == 1:
            sourceXPosition, sourceYPosition, sourceZPosition = \
                _get_point_source_position(unitFactor)
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
        defaultPath = os.path.join(programPath, 'Mesh2Input', 'Materials')
        materialSearchPaths += f";  {defaultPath}"

        # get used materials as dictionary
        materials = _get_materials(bpy.data.objects['Reference'])

        if materials is None:
            if sourceType_id == 0:
                raise ValueError(
                    ("Material 'Left ear' and/or 'Right ear' "
                     "must be defined for reciprocal simulations."))
        else:
            # add full path of material files to the dictionary
            materials = _get_material_files(materials, materialSearchPaths)

            # read the material data and add it to the dictionary
            materials = _read_material_data(materials)


# Calculate frequency information ---------------------------------------------
        # maximum numbers are used for initialization & status print-outs:
        maxCPUs = cpuLast             # used to be max(10, cpuLast)
        maxCores = numCoresPerCPU     # used to be max(8, numCoresPerCPU)


        # check how the frequency vector is defined
        if frequencyVectorType == 'Step size':
            frequencyStepSize = frequencyVectorValue
            numFrequencySteps = 0
        else:
            frequencyStepSize = 0
            numFrequencySteps = int(frequencyVectorValue)

        cpusAndCores, frequencies, frequencyStepsPerCore, numCoresAvailable, \
            freqs, frequencyStepSize, numFrequencySteps = \
            _distribute_frequencies(cpuFirst, cpuLast, maxCPUs,
                                    numCoresPerCPU, maxCores,
                                    numSources,
                                    minFrequency, maxFrequency,
                                    frequencyStepSize, numFrequencySteps)


# Write Info.txt (feedback for user, not used by NumCalc) ---------------------
        _write_info_txt(evalGridPaths, title, ear, filepath1, version,
                        cpuFirst, cpuLast, numCoresAvailable,
                        frequencyStepsPerCore, frequencies, freqs,
                        frequencyStepSize, numFrequencySteps, maxCores, maxCPUs)


# Write Output2HRTF.m function ------------------------------------------------
        _write_output2HRTF_m(filepath1, version, sourceType_id, ear, unitFactor,
                             reference, computeHRIRs,
                             speedOfSound, densityOfMedium,
                             cpusAndCores, maxCPUs, maxCores,
                             sourceXPosition, sourceYPosition, sourceZPosition)

# Write Output2VTK.m function ------------------------------------------------
        _write_output2VTK_m(filepath1, version)

# Write Output2HRTF.py function ------------------------------------------------
        _write_output2HRTF_py(filepath1, version, sourceType_id, ear, unitFactor,
                             reference, computeHRIRs,
                             speedOfSound, densityOfMedium,
                             cpusAndCores, maxCPUs, maxCores,
                             sourceXPosition, sourceYPosition, sourceZPosition,
                             programPath)


# Render pictures of the model ------------------------------------------------
        if pictures:
            _render_pictures(filepath1, unitFactor)


# Write NumCalc input files for all CPUs and Cores (NC.inp) -------------------
        _write_nc_inp(filepath1, version, title, ear, speedOfSound,
                      densityOfMedium, frequencies, cpusAndCores,
                      evaluationGrids, materials, method, sourceType_id,
                      sourceXPosition, sourceYPosition, sourceZPosition)

# Finish ----------------------------------------------------------------------
        for obj in bpy.context.scene.objects[:]:
            bpy.data.objects[obj.name].select_set(False)

        return {'FINISHED'}


def _get_point_source_position(unitFactor):
    # check if 'Reference' object exists and select it
    pointSourceExist = False
    for obj in bpy.context.scene.objects[:]:
        if obj.type == 'LIGHT' and obj.name == 'Point':
            pointSourceExist = True

    if not pointSourceExist:
        raise ValueError("Did not find the point source. It must be 'Point'"
                         "light object named 'Point' (case sensitive).")

    sourceXPosition = bpy.data.objects["Point"].location[0] * unitFactor
    sourceYPosition = bpy.data.objects["Point"].location[1] * unitFactor
    sourceZPosition = bpy.data.objects["Point"].location[2] * unitFactor

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
        programPath, "Mesh2Input", "EvaluationGrids")

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

    # remove default value if user path is not passed
    if paths[0] == "None":
        paths.pop(0)

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


def _write_info_txt(evalGridPaths, title, ear, filepath1, version,
                    cpuFirst, cpuLast, numCoresAvailable,
                    frequencyStepsPerCore, frequencies, freqs,
                    frequencyStepSize, numFrequencySteps, maxCores, maxCPUs):
    file = open(os.path.join(filepath1, "Info.txt"), "w",
                encoding="utf8", newline="\n")
    fw = file.write
    fw("#####################################\n")
    fw("######## General information ########\n")
    fw("#####################################\n\n")
    fw(f"Program: Mesh2HRTF {version}\n")
    fw("Title: %s\n" % title)
    fw("Ear: %s\n" % ear)
    fw("Evaluation Grids:\n")
    for evaluationGrid in evalGridPaths:
        fw("    %s\n" % evaluationGrid)
    fw("\n")
    fw("#####################################\n")
    fw("####### Frequency information #######\n")
    fw("#####################################\n\n")
    fw("Minimum evaluated Frequency: %f\n" % freqs[0])
    fw("Highest evaluated Frequency: %f\n" % freqs[-1])
    fw("Frequency Stepsize: %f\n" % frequencyStepSize)
    fw("Frequency Steps: %d\n" % numFrequencySteps)
    fw("Frequency steps per Core: %d\n\n" % frequencyStepsPerCore)
    fw("#####################################\n")
    fw("######## Cluster information ########\n")
    fw("#####################################\n\n")
    fw("Number of CPUs: %d\n" % (cpuLast-cpuFirst+1))
    fw("First CPU: 'CPU_%d'\n" % cpuFirst)
    fw("Last CPU: 'CPU_%d'\n" % cpuLast)
    fw("Number of Cores (available): %d\n" % numCoresAvailable)
    for cpu in range(1, maxCPUs+1):
        for core in range(1, maxCores+1):
            fw("CPU_%d (Core %d):\n" % (cpu, core))
            for ii in range(0, len(frequencies[cpu-1][core-1])):
                fw("    %f\n" % frequencies[cpu-1][core-1][ii])
        fw("\n")
    file.close


def _write_output2HRTF_m(filepath1, version,
                         sourceType_id, ear, unitFactor,
                         reference, computeHRIRs,
                         speedOfSound, densityOfMedium,
                         cpusAndCores, maxCPUs, maxCores,
                         sourceXPosition, sourceYPosition, sourceZPosition):

    # file handling
    file = open(os.path.join(filepath1, "Output2HRTF.m"), "w",
                encoding="utf8", newline="\n")
    fw = file.write

    # header
    fw("% Collect the data simulated by NumCalc and save to project folder.\n")
    fw("close all; clear\n\n")

    # Mesh2HRTF version
    fw(f"Mesh2HRTF_version = '{version}';\n\n")

    # add information about the source
    if sourceType_id == 0:
        fw("% source information\n")
        fw("sourceType = 'vibratingElement';\n")

        # get the receiver/ear centers and areas
        obj = bpy.data.objects['Reference']
        obj_data = obj.data
        earCenter, earArea = _calculateReceiverProperties(
            obj, obj_data, unitFactor)

        # write left ear data
        if ear == 'Left ear' or ear == 'Both ears':
            fw("sourceCenter(1,1:3) = [%f %f %f];\n" % (earCenter[0][0],
                                                        earCenter[0][1],
                                                        earCenter[0][2]))
            fw("sourceArea(1,1) = %g;\n" % earArea[0])

        # write right ear data
        if ear == 'Right ear' or ear == 'Both ears':
            if ear == 'Right ear':
                nn = 1
            if ear == 'Both ears':
                nn = 2

            fw("sourceCenter(%d,1:3) = [%f %f %f];\n" % (nn,
                                                         earCenter[1][0],
                                                         earCenter[1][1],
                                                         earCenter[1][2]))
            fw("sourceArea(%d,1) = %g;\n" % (nn, earArea[1]))

        fw("\n")
    else:

        fw("% source information\n")
        fw("sourceType = 'pointSource';\n")
        fw("sourceCenter(1,1:3) = [%s %s %s];\n" % (sourceXPosition,
                                                    sourceYPosition,
                                                    sourceZPosition))
        fw("sourceArea(1,1)     = 1;\n")

    # referencing
    fw("% Reference to a point source in the origin\n")
    fw("% accoring to the classical HRTF definition\n")
    fw("% (https://doi.org/10.1016/0003-682X(92)90046-U)\n")
    fw("reference = ")
    if reference:
        fw("true;\n\n")
    else:
        fw("false;\n\n")

    # compute HRIRs
    fw("% Compute HRIRs via the inverse Fourier transfrom.\n")
    fw("% This will add data at 0 Hz, mirror the single sided spectrum, and\n")
    fw("% shift the HRIRs in time. Requires reference = true.\n")
    fw("computeHRIRs = ")
    if computeHRIRs:
        fw("true;\n\n")
    else:
        fw("false;\n\n")

    # constants
    fw("% Constants\n")
    fw("speedOfSound = " + speedOfSound + "; % [m/s]\n")
    fw("densityOfAir = " + densityOfMedium + "; % [kg/m^3]\n\n")

    # write CPUs and Cores
    fw("% Distribution of ears across CPUs and cores.\n")
    fw("% (Matrix of size [numCPUs x numCores])\n")
    fw("cpusAndCores = [\n")
    for cpu in range(1, maxCPUs + 1):
        fw("    ")
        for core in range(1, maxCores + 1):
            fw("%i" % cpusAndCores[cpu-1][core-1])
            if core < maxCores:
                fw(" ")
        if cpu < maxCPUs:
            fw("; ...\n")
    fw("];\n\n")

    fw("% Collect the data simulated by NumCalc\n")
    fw("Output2HRTF_Main(Mesh2HRTF_version, cpusAndCores, ...\n")
    fw("                 sourceType, sourceCenter, sourceArea, ...\n")
    fw("                 reference, computeHRIRs, ...\n")
    fw("                 speedOfSound,densityOfAir);\n")
    file.close

def _write_output2VTK_m(filepath1, version):

    # file handling
    file = open(os.path.join(filepath1, "Output2VTK.m"), "w",
                encoding="utf8", newline="\n")
    fw = file.write

    # header
    fw("% Export the sound pressure files genereated by Output2HRTF.m as SPL to VTK-Files for visualization in Paraview.\n\n")

    fw("close all; clear\n\n")

    # Mesh2HRTF version
    fw(f"Mesh2HRTF_version = '{version}';\n\n")

    # create missing directories
    fw("if ~exist(fullfile(pwd, 'Visualization'),'dir')\n")
    fw("    mkdir(fullfile(pwd, 'Visualization'));\n")
    fw("end\n")
    fw("if ~exist(fullfile(pwd, 'Visualization', 'ObjectMesh'),'dir')\n")
    fw("    mkdir(fullfile(pwd, 'Visualization', 'ObjectMesh'))\n")
    fw("end\n\n")

    # load data struct created by Output2HRTF.m
    fw("load(fullfile(pwd, 'Output2HRTF', 'ObjectMesh_Reference.mat'))\n\n")

    # export the sound pressure distribution as dB SPL in VTK-files
    fw("EvalToolsExport2VTK(['Visualization' filesep 'ObjectMesh' filesep],nodes(:,2:end),elements(:,2:end),20*log10(abs(element_data{1})/0.00002),'amp')\n")
    file.close

def _write_output2HRTF_py(filepath1, version,
                         sourceType_id, ear, unitFactor,
                         reference, computeHRIRs,
                         speedOfSound, densityOfMedium,
                         cpusAndCores, maxCPUs, maxCores,
                         sourceXPosition, sourceYPosition, sourceZPosition,
                         programPath):

    # file handling
    file = open(os.path.join(filepath1, "Output2HRTF.py"), "w",
                encoding="utf8", newline="\n")
    fw = file.write

    # header
    fw("# Read the data simulated by NumCalc and save to the folder\n")
    fw("# Output2HRTF inside project folder.\n\n")

    fw("import numpy\n")
    fw("import os\n")
    fw("import mesh2hrtf as m2h\n\n")

    fw("projectPath = os.getcwd()\n\n")

    # Mesh2HRTF version
    fw(f"Mesh2HRTF_version = '{version}'\n\n")

    fw("# source information\n")
    # initialize arrays for source information
    if sourceType_id == 0 and ear == 'Both ears':
        fw("sourceCenter = numpy.zeros((2, 3))\n")
        fw("sourceArea = numpy.zeros((2, 1))\n\n")
    else:
        fw("sourceCenter = numpy.zeros((1, 3))\n")
        fw("sourceArea = numpy.zeros((1, 1))\n\n")

    # add information about the source
    if sourceType_id == 0:
        fw("sourceType = 'vibratingElement'\n")

        # get the receiver/ear centers and areas
        obj = bpy.data.objects['Reference']
        obj_data = obj.data
        earCenter, earArea = _calculateReceiverProperties(
            obj, obj_data, unitFactor)

        # write left ear data
        if ear == 'Left ear' or ear == 'Both ears':
            fw("sourceCenter[0, :] = [%f, %f, %f]\n"
                                            % (earCenter[0][0],
                                               earCenter[0][1],
                                               earCenter[0][2]))
            fw("sourceArea[0, 0] = %g\n" % earArea[0])

        # write right ear data
        if ear == 'Right ear' or ear == 'Both ears':
            if ear == 'Right ear':
                nn = 0
            if ear == 'Both ears':
                nn = 1

            fw("sourceCenter[%d, :] = [%f, %f, %f]\n"
                                            % (nn,
                                               earCenter[1][0],
                                               earCenter[1][1],
                                               earCenter[1][2]))
            fw("sourceArea[%d, 0] = %g\n" % (nn, earArea[1]))

        fw("\n")
    else:
        fw("sourceType = 'pointSource';\n")
        fw("sourceCenter[0, :] = [%s, %s, %s]\n"
                                            % (sourceXPosition,
                                               sourceYPosition,
                                               sourceZPosition))
        fw("sourceArea[0, 0]     = 1\n")

    # referencing
    fw("# Reference to a point source in the origin\n")
    fw("# accoring to the classical HRTF definition\n")
    fw("# (https://doi.org/10.1016/0003-682X(92)90046-U)\n")
    fw("reference = ")
    if reference:
        fw("True\n\n")
    else:
        fw("False\n\n")

    # compute HRIRs
    fw("# Compute HRIRs via the inverse Fourier transfrom.\n")
    fw("# This will add data at 0 Hz, mirror the single sided spectrum, and\n")
    fw("# shift the HRIRs in time. Requires reference = true.\n")
    fw("computeHRIRs = ")
    if computeHRIRs:
        fw("True\n\n")
    else:
        fw("False\n\n")

    # constants
    fw("# Constants\n")
    fw("speedOfSound = " + speedOfSound + "  # [m/s]\n")
    fw("densityOfAir = " + densityOfMedium + "  # [kg/m^3]\n\n")

    # write CPUs and Cores
    fw("# Distribution of ears across CPUs and cores.\n")
    fw("# (Matrix of size [numCPUs x numCores])\n")
    fw("cpusAndCores = numpy.array([\n")
    for cpu in range(1, maxCPUs + 1):
        fw("    [")
        for core in range(1, maxCores + 1):
            fw("%i" % cpusAndCores[cpu-1][core-1])
            if core < maxCores:
                fw(", ")
        if cpu < maxCPUs:
            fw("],\n")
    fw("]])\n\n")

    fw("# Collect the data simulated by NumCalc\n")
    fw("m2h.Output2HRTF_Main(projectPath, Mesh2HRTF_version, cpusAndCores,\n")
    fw("                 sourceType, sourceCenter, sourceArea,\n")
    fw("                 reference, computeHRIRs,\n")
    fw("                 speedOfSound, densityOfAir)\n")
    file.close

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


def _distribute_frequencies(cpuFirst, cpuLast, maxCPUs,
                            numCoresPerCPU, maxCores,
                            numSources,
                            minFrequency, maxFrequency,
                            frequencyStepSize, numFrequencySteps):
    """Calculate list of frequencies and distribute across cores and CPUs.

    Returns
    -------

    cpusAndCores: list
        Nested list that indicates which cpu and core is used to calculate data
        for which ear (e.g. cpuAndCores[0][1] holds the entry for the second
        core on the first cpu). Entries: 0=idle, 1=leftEar/right ear if
        calculating one ear, 2=rightEar if calculating two ears.
    frequencies: list
        Nested list that holds the frequencies that are calculated by each
        cpu/core (e.g. frequencies[0][1] holds a list of frequencies calculated
        by the second core on the first cpu.
    frequencyStepsPerCore: int
        The number of frequency steps calculated per core (written to
        Info.txt).
    numCoresAvailable: int
        The number of cores used for the computation (written to Info.txt).
    f: list
        A simple list of the frequencies to be simulated (for debugging).
    frequencyStepSize: float
        Step size between successive frequencies (written to Info.txt). This is
        returned because it might be zero during the function call.
    numFrequencySteps: int
        Number of frequencies to be simulated (written to Info.txt). This is
        returned because it might be zero during the function call.

    """

    # number of CPUs used
    numCPUs = cpuLast-cpuFirst+1

    # number of cores per ear
    numCoresUsedPerEar = numCPUs*numCoresPerCPU//numSources
    if not numCoresUsedPerEar:
        raise Exception("At least two cores must be available for calculating "
                        "both ears, i.e., two CPUs with one core each or one "
                        "CPU with two cores.")

    # check input
    if (numFrequencySteps == 0 and frequencyStepSize == 0) \
            or (numFrequencySteps != 0 and frequencyStepSize != 0):
        raise Exception("Either 'frequencyStepSize' or 'numFrequencySteps' \
                        must be zero while the other must not.")

    # Calculate Number of frequencies and frequency step size
    if minFrequency == maxFrequency:
        frequencySteps = (1, 0)
        frequencyStepSize = 0
    elif frequencyStepSize:
        frequencySteps = divmod(
            maxFrequency - minFrequency + frequencyStepSize, frequencyStepSize)
    else:
        if numFrequencySteps < 2:
            raise Exception("'numFrequencySteps' must be at least 2.")
        frequencySteps = (numFrequencySteps, 0)
        frequencyStepSize = (maxFrequency-minFrequency)/(numFrequencySteps-1)

    if not frequencySteps[1] == 0:
        raise Exception("Error, frequencyStepSize is not a divisor of \
                        maxFrequency-minFrequency")

    # get all frequencies to be calculated
    f = [ff*frequencyStepSize+minFrequency
            for ff in range(int(frequencySteps[0]))]

    # remove 0 Hz if included in the list
    if f[0] == 0:
        f.pop(0)
        frequencySteps = (frequencySteps[0]-1, 0)
        print('Warning: 0 Hz can not be calculated and was removed from the \
              list of frequencies.')

    if not len(f):
        raise ValueError("No frequencies to be calculated. \
                         Check the input parameters.")

    # check number of cores and frequencies
    if len(f) < numCoresUsedPerEar:
        raise Exception("More cores than frequencies, i.e., \
                        numCPUs*numCoresPerCPU//numSources < numFrequencies.")

    # distribution of frequencies across numCoresUsedPerEar
    F = [[] for ff in range(numCoresUsedPerEar)]
    for nn, ff in enumerate(f):
        F[nn % numCoresUsedPerEar].append(ff)

    # Initialize cpusAndCores:
    cpusAndCores = [[0]*maxCores for cc in range(maxCPUs)]

    # Initialize frequencies:
    freqs = [[] for cc in range(maxCores)]
    frequencies = [freqs.copy() for cc in range(maxCPUs)]

    # distribute ears and frequencies across cpus and cores.
    # Left ear is calculated on cpus 0 to numCoresUsedPerEar-1
    # Right ear is calculated on cpus numCoresUsedPerEar to numCoresAvailable
    for count in range(numCoresUsedPerEar*numSources):
        cpu, core = divmod(count + (cpuFirst-1)*numCoresPerCPU, numCoresPerCPU)
        cpusAndCores[cpu][core] = count//numCoresUsedPerEar + 1
        frequencies[cpu][core] = F[count % len(F)]
        # output for debugging
        # print(f"CPU {cpu+1:2d}, core {core+1}, \
        #       ear {count//numCoresUsedPerEar + 1}, \
        #       freqList {count%len(F)}")

    # meta data for Info.txt
    frequencyStepsPerCore = len(f)//numCoresUsedPerEar
    numCoresAvailable = numCoresUsedPerEar*numSources

    return cpusAndCores, frequencies, \
        frequencyStepsPerCore, numCoresAvailable, \
        f, frequencyStepSize, numFrequencySteps


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
    bpy.data.lights['Light'].energy = 800 / unitFactor** 2
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


def _write_nc_inp(filepath1, version, title, ear,
                  speedOfSound, densityOfMedium, frequencies, cpusAndCores,
                  evaluationGrids, materials, method, sourceType_id,
                  sourceXPosition, sourceYPosition, sourceZPosition):
    """Write NC.inp file that is read by NumCalc to start the simulation.

    The file format is documented at:
    https://sourceforge.net/p/mesh2hrtf/wiki/Structure%20of%20NC.inp_0.4.0/
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

    for core in range(len(cpusAndCores[0])):
        for cpu in range(len(cpusAndCores)):
            if not cpusAndCores[cpu][core] == 0:

                # create directoy
                filepath2 = os.path.join(
                    filepath1, "NumCalc", "CPU_%i_Core_%i" % (cpu+1, core+1))
                if not os.path.exists(filepath2):
                    os.mkdir(filepath2)

                # write NC.inp
                file = open(os.path.join(filepath2, "NC.inp"), "w",
                            encoding="utf8", newline="\n")
                fw = file.write

                obj_name = "Reference"

                obj = bpy.data.objects[obj_name]
                obj_data = obj.data

                # header ------------------------------------------------------
                fw("##-------------------------------------------\n")
                fw("## This file was created by export_mesh2hrtf\n")
                fw("## Date: %s\n" % datetime.date.today())
                fw("##-------------------------------------------\n")
                fw("Mesh2HRTF %s\n" % version)
                fw("##\n")
                fw("%s\n" % title)
                fw("##\n")

                # control parameter I (hard coded, not documented) ------------
                fw("## Controlparameter I\n")
                fw("0 0 0 0 7 0\n")
                fw("##\n")

                # control parameter II ----------------------------------------
                fw("## Controlparameter II\n")
                fw("1 %d 0.000001 0.00e+00 1 0 0\n" % (
                    len(frequencies[cpu][core])))
                fw("##\n")
                fw("## Load Frequency Curve \n")
                fw("0 %d\n" % (len(frequencies[cpu][core])+1))
                fw("0.000000 0.000000e+00 0.0\n")
                for ii in range(len(frequencies[cpu][core])):
                    fw("%f %fe+04 0.0\n" % (
                        0.000001*(ii+1),
                        frequencies[cpu][core][ii] / 10000))
                fw("##\n")

                # main parameters I -------------------------------------------
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

                # main parameters II ------------------------------------------
                fw("## 2. Main Parameters II\n")
                fw("0 ")
                if sourceType_id==0:
                    fw("0 ")
                else:
                    fw("1 ")
                fw("0 0.0000e+00 0 0 0\n")
                fw("##\n")

                # main parameters III -----------------------------------------
                fw("## 3. Main Parameters III\n")
                fw("0 0 0 0\n")
                fw("##\n")

                # main parameters IV ------------------------------------------
                fw("## 4. Main Parameters IV\n")
                fw("%s %se+00 1.0 0.0e+00 0.0 e+00 0.0e+00 0.0e+00\n" % (
                    speedOfSound, densityOfMedium))
                fw("##\n")

                # nodes -------------------------------------------------------
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

                # SYMMETRY ----------------------------------------------------
                fw("# SYMMETRY\n")
                fw("# 0 0 0\n")
                fw("# 0.0000e+00 0.0000e+00 0.0000e+00\n")
                fw("##\n")

                # boundary information ----------------------------------------
                fw("BOUNDARY\n")
                # write velocity condition for the ears if using vibrating
                # elements as the sound source
                if sourceType_id==0:
                    if cpusAndCores[cpu][core]==1 and ear!='Right ear':
                        tmpEar='Left ear'
                    else:
                        tmpEar='Right ear'
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

                # source information: plane wave ------------------------------
                fw("# PLANE WAVES\n")
                fw("# 0 0.0000e+00 -1.0000e+00 0.0000e+00 1.0000e-6 -1 0.0000e+00 -1\n")
                fw("##\n")

                # source information: point source ----------------------------
                if sourceType_id==0:
                    fw("# POINT SOURCES\n")
                    if cpusAndCores[cpu][core] == 1:
                        fw("# 0 0.0 0.101 0.0 0.1 -1 0.0 -1\n")
                    if cpusAndCores[cpu][core] == 2:
                        fw("# 0 0.0 -0.101 0.0 0.1 -1 0.0 -1\n")
                else:
                    fw("POINT SOURCES\n")
                    fw("0 %s %s %s 0.1 -1 0.0 -1\n" % (
                        sourceXPosition, sourceYPosition, sourceZPosition))
                fw("##\n")

                # curves ------------------------------------------------------
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

                # post process ------------------------------------------------
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
