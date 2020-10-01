#                                Mesh2HRTF
#                Copyright (C) 2015 by Harald Ziegelwanger,
#        Acoustics Research Institute, Austrian Academy of Sciences
#                        mesh2hrtf.sourceforge.net
#
# Mesh2HRTF is licensed under the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# Mesh2HRTF is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.
# You should have received a copy of the GNU LesserGeneral Public License along with Mesh2HRTF. If not, see <http://www.gnu.org/licenses/lgpl.html>.
#
# If you use Mesh2HRTF:
# - Provide credits:
#   "Mesh2HRTF, H. Ziegelwanger, ARI, OEAW (mesh2hrtf.sourceforge.net)"
# - In your publication, cite both articles:
#   [1] Ziegelwanger, H., Kreuzer, W., and Majdak, P. (2015). "Mesh2HRTF: Open-source software package for the numerical calculation of head-related transfer functions," in Proceedings of the 22nd ICSV, Florence, IT.
#   [2] Ziegelwanger, H., Majdak, P., and Kreuzer, W. (2015). "Numerical calculation of listener-specific head-related transfer functions and sound localization: Microphone model and mesh discretization," The Journal of the Acoustical Society of America, 138, 208-222.
#
# Authors: The Mesh2HRTF developers

import os
import bpy
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
    "version": (0, 2, 0),
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
    ear: EnumProperty(
        name="Ear",
        description="Selected ear(s) for simulation",
        items=[('Left ear', 'left', 'Left ear'),
               ('Right ear', 'right', 'Right ear'),
               ('Both ears', 'both', 'Both ears'),
               ('None', 'none', 'None')],
        default='Both ears',
        )
    method: EnumProperty(
        name="BEM method",
        description="Method for numerical simulation",
        items=[('0', 'BEM', 'Traditional BEM'),
               ('1', 'SL-FMM BEM', 'Singlelevel fast-multipole method BEM'),
               ('4', 'ML-FMM BEM', 'Multilevel fast-multipole method BEM')],
        default='4',
        )
    programPath: StringProperty(
        name="Mesh2HRTF-path",
        description="Path to folder containing 'Mesh2Input', 'NumCalc', etc.",
        default=r"path/to/mesh2hrtf",
        )
    pictures: BoolProperty(
        name="Pictures",
        description="Render pictures of the 3D mesh",
        default=True,
        )
    # source type -------------------------------------------------------------
    reciprocity: BoolProperty(
        name="Reciprocity",
        description=("Calculation with reciprocity. i.e., with vibrating "
            "elements as the source. (This disables the point source mode, "
            "i.e. the x, y, and z coordinate below"),
        default=True,
        )
    reference: BoolProperty(
        name="Reference",
        description=("Reference HRTF to the center of the head according to "
            "the classic HRTF definition. Only applicable if using "
            "reciprocity. For the HRTF definition see"
            "https://doi.org/10.1016/0003-682X(92)90046-U"),
        default=False,
        )
    sourceXPosition: StringProperty(
        name="Source (x)",
        description="Point source Position x-coordinate (according to units)",
        default="0",
        )
    sourceYPosition: StringProperty(
        name="Source (y)",
        description="Point source Position y-coordinate (according to units)",
        default="101",
        )
    sourceZPosition: StringProperty(
        name="Source (z)",
        description="Point source Position z-coordinate (according to units)",
        default="0",
        )
    # constants ---------------------------------------------------------------
    unit: EnumProperty(
        name="Unit",
        description="Unit of the 3D mesh",
        items=[('m', 'm', 'Meter'), ('mm', 'mm', 'Millimeter')],
        default='mm',
        )
    speedOfSound: StringProperty(
        name="c (m/s)",
        description="Speed of sound (m/s)",
        default="346.18",
        )
    densityOfMedium: StringProperty(
        name="rho ()",
        description="Density of air (kg/m^3)",
        default="1.1839",
        )
    # evaluation grids --------------------------------------------------------
    evaluationGrids: StringProperty(
        name="Name",
        description=("Name of evalation grid inside "
            "Mesh2Input/EvaluationsGrids or absolute path to user grid. "
            "Multiple grids can be separated by semicolons (;)."),
        default='ARI',
        )
    # Frequency selection -----------------------------------------------------
    minFrequency: FloatProperty(
        name="Min. frequency",
        description=("Minimum frequency in Hz to be simulated. Can be 0 for "
            "constructing single sided spectra. But the 0 will not be"
            "simulated."),
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
    frequencyStepSize: FloatProperty(
        name="Step size",
        description=("Simulate frequencies between the minimum and maximum "
            "frequency with this step size in Hz. Either this or the 'Steps' "
            "below must be zero."),
        default=100,
        min=0,
        max=24000,
        )
    numFrequencySteps: IntProperty(
        name="Steps",
        description=("Simulate N frequencies between the minimum and maximum "
            "frequency. Either this or the 'Step size' above must be zero."),
        default=0,
        min=0,
        max=24000,
        )
    # Job distribution --------------------------------------------------------
    cpuFirst: IntProperty(
        name="CPU (first)",
        description="First 'CPU' used",
        default=1,
        min=1,
        max=100,
        )
    cpuLast: IntProperty(
        name="CPU (last)",
        description="Last 'CPU' used",
        default=10,
        min=1,
        max=100,
        )
    numCoresPerCPU: IntProperty(
        name="Num. of used cores",
        description="Number of used cores per CPU",
        default=8,
        min=1,
        max=8,
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
        row.prop(self, "ear")
        row = layout.row()
        row.prop(self, "method")
        row = layout.row()
        row.prop(self, "programPath")
        row = layout.row()
        row.prop(self, "pictures")
        # source type
        layout.label(text="Source type:")
        row = layout.row()
        row.prop(self, "reciprocity")
        row = layout.row()
        row.prop(self, "reference")
        row = layout.row()
        row.prop(self, "sourceXPosition")
        row = layout.row()
        row.prop(self, "sourceYPosition")
        row = layout.row()
        row.prop(self, "sourceZPosition")
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
        # frequency distribution
        layout.label(text="Frequency distribution:")
        row = layout.row()
        row.prop(self, "minFrequency")
        row = layout.row()
        row.prop(self, "maxFrequency")
        row = layout.row()
        row.prop(self, "frequencyStepSize")
        row = layout.row()
        row.prop(self, "numFrequencySteps")
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
             frequencyStepSize=100,
             numFrequencySteps=0,
             cpuFirst=1,
             cpuLast=10,
             numCoresPerCPU=8,
             pictures=True,
             ear='Both ears',
             evaluationGrids='ARI',
             method='4',
             reciprocity=True,
             reference=False,
             sourceXPosition='0',
             sourceYPosition='101',
             sourceZPosition='0',
             speedOfSound='346.18',
             densityOfMedium='1.1839',
             unit='mm',
             programPath="",
             ):
        """Export Mesh2HRTF project."""


# General handling and constants ----------------------------------------------
        # Switch to object mode to avoid export errors
        bpy.ops.object.mode_set(mode='OBJECT', toggle=False)

        # get Mesh2HRTF version
        with open(os.path.join(programPath, "..", "VERSION")) as read_version:
            version = read_version.readline()

        # Export path and export directory handling
        (filepath1, filename1) = os.path.split(filepath)
        filename1 = "NC.inp"

        temp = os.path.join(filepath1, "ObjectMeshes")
        if not os.path.exists(temp):
            os.mkdir(temp)

        temp = os.path.join(filepath1, "EvaluationGrids")
        if not os.path.exists(temp):
            os.mkdir(temp)

        temp = os.path.join(filepath1, "NumCalc")
        if not os.path.exists(temp):
            os.mkdir(temp)

        # check input and assign numEars
        if ear not in ["Left ear", "Right ear", "Both ears"]:
            raise ValueError("`ear` must be 'Left ear', 'Right ear' or "
                             "'Both ears' (case sensitive).")
        numEars = 1 if ear == 'Both ears' else 2

        # check input and assign unitFactor
        if unit not in ["mm", "m"]:
            raise ValueError("`unit` must be 'mm', 'm' "
                             "(case sensitive).")
        unitFactor = 1 if unit == 'm' else .001

        # camera settings for rending pictures
        bpy.ops.object.transform_apply(location=True)
        bpy.ops.object.transform_apply(rotation=True)
        bpy.ops.object.transform_apply(scale=True)
        cam = bpy.data.objects['Camera']
        camradius = 400
        cam.location = (0, camradius, 0)
        cam.rotation_euler = (pi/2, 0, pi)
        cam.data.clip_end = 0.1
        cam.data.clip_end = 1000
        light = bpy.data.objects['Light']
        lightradius = 300
        light.location = (0, lightradius, 0)
        bpy.data.lights['Light'].energy = 800
        bpy.data.lights['Light'].distance = 100
        renderloc = [[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0],
                     [0.707, 0.707, 0], [-0.707, 0.707, 0], [0.707, -0.707, 0],
                     [-0.707, -0.707, 0]]
        renderrot = [[pi/2, 0, pi/2], [pi/2, 0, 3*pi/2], [pi/2, 0, pi],
                     [3/2*pi, pi, pi], [pi/2, 0, 3/4*pi], [pi/2, 0, 5/4*pi],
                     [pi/2, 0, pi/4], [pi/2, 0, -pi/4]]
        rendernam = [[0, 0], [180, 0], [90, 0], [270, 0],
                     [45, 0], [135, 0], [315, 0], [225, 0]]

        bpy.data.scenes['Scene'].render.pixel_aspect_x = 1
        bpy.data.scenes['Scene'].render.pixel_aspect_y = 1
        bpy.data.scenes['Scene'].render.resolution_x = 1440
        bpy.data.scenes['Scene'].render.resolution_y = 1920

        # empty list for saving objects
        objects = ([])


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


# Calculate frequency information ---------------------------------------------
        # maximum number of cpus. Still hard coded but might be:
        # `maxCPUs = max(10, cpuLast)` if it works with the rest of Mesh2HRTF
        maxCPUs = 10

        # maximum number of cores. Still hard coded but might be:
        # maxCores = max(8, numCoresPerCPU) if it works with the rest of
        # Mesh2HRTF
        maxCores = 8

        cpusAndCores, frequencies, frequencyStepsPerCore, numCoresAvailable, \
            freqs, frequencyStepSize, numFrequencySteps = \
            _distribute_frequencies(cpuFirst, cpuLast, maxCPUs,
                                    numCoresPerCPU, maxCores,
                                    numEars,
                                    minFrequency, maxFrequency,
                                    frequencyStepSize, numFrequencySteps)


# Write Info.txt (feedback for user, not used by NumCalc) ---------------------
        _write_info_txt(evalGridPaths, title, ear, filepath1,
                        cpuFirst, cpuLast, numCoresAvailable,
                        frequencyStepsPerCore, frequencies, freqs,
                        frequencyStepSize, numFrequencySteps)


# Write Output2HRTF.m function ------------------------------------------------
        _write_output2HRTF_m(filepath1, version, reciprocity, ear, unitFactor,
                             reference, speedOfSound, densityOfMedium,
                             cpusAndCores, maxCPUs, maxCores,
                             sourceXPosition, sourceYPosition, sourceZPosition)


# ----------------------- Render pictures of the model -------------------------
        if pictures:
            for ii in range(0, len(renderloc)):
                cam.location = (renderloc[ii][0]*camradius, renderloc[ii][1]*camradius, renderloc[ii][2]*camradius)
                cam.rotation_euler = (renderrot[ii][0], renderrot[ii][1], renderrot[ii][2])
                light.location = (renderloc[ii][0]*lightradius, renderloc[ii][1]*lightradius, renderloc[ii][2]*lightradius)
                bpy.ops.render.render()
                temp = os.path.join(filepath1, "Pictures")
                if not os.path.exists(temp):
                    os.mkdir(temp)
                temp = ("%d-%d" % (rendernam[ii][0], rendernam[ii][1]))
                bpy.data.images['Render Result'].save_render(os.path.join(filepath1, "Pictures", "%s.png" % temp))

# ----------------------- Write NumCalc input files for all CPUs and Cores -----
        for core in range(1, 9):
            for cpu in range(1, 11):
                if not cpusAndCores[cpu-1][core-1] == 0:

                    filepath2 = os.path.join(filepath1, "NumCalc", "CPU_%i_Core_%i" % (cpu, core))
                    if not os.path.exists(filepath2):
                        os.mkdir(filepath2)

                    file = open(os.path.join(filepath2, filename1), "w", encoding="utf8", newline="\n")
                    fw = file.write

                    obj_name = "Reference"

                    obj = bpy.data.objects[obj_name]
                    obj_data = obj.data

                    fw("##-------------------------------------------\n")
                    fw("## This file was created by export_mesh2hrtf\n")
                    fw("## Date: %s\n" % datetime.date.today())
                    fw("##-------------------------------------------\n")
                    fw("Mesh2HRTF %s\n" % version)
                    fw("##\n")
                    fw("%s\n" % title)
                    fw("##\n")
                    fw("## Controlparameter I\n")
                    fw("0 0 0 0 7 0\n")
                    fw("##\n")
                    fw("## Controlparameter II\n")
                    fw("1 %d %fe+00 0.00e+00 1 0 0\n" % (len(frequencies[cpu-1][core-1]), 1/(len(frequencies[cpu-1][core-1]))))
                    fw("##\n")
                    fw("## Load Frequency Curve \n")
                    fw("0 %d\n" % (len(frequencies[cpu-1][core-1])+1))
                    fw("0.000000e+00 0.000000e+00 0.0\n")
                    for ii in range(0, len(frequencies[cpu-1][core-1])):
                        fw("%fe+00 %fe+04 0.0\n" % (1/(len(frequencies[cpu-1][core-1]))*(ii+1), frequencies[cpu-1][core-1][ii]/10000))
                    fw("##\n")
                    fw("## 1. Main Parameters I\n")
                    numNodes = 0
                    numElements = 0
                    for evaluationGrid in evaluationGrids:
                        # read number of nodes
                        nodes = open(os.path.join(filepath1, "EvaluationGrids", evaluationGrid, "Nodes.txt"))
                        line = nodes.readline()
                        numNodes = numNodes+int(line)
                        # read number of elements
                        elements = open(os.path.join(filepath1, "EvaluationGrids", evaluationGrid, "Elements.txt"))
                        line = elements.readline()
                        numElements = numElements+int(line)
                    fw("2 %d " % (len(obj_data.polygons[:])+numElements))
                    fw("%d 0 " % (len(obj_data.vertices[:])+numNodes))
                    fw("0")
                    fw(" 2 1 %s 0\n" % (method))
                    fw("##\n")
                    fw("## 2. Main Parameters II\n")
                    fw("0 ")
                    if reciprocity:
                        fw("0 ")
                    else:
                        fw("1 ")
                    fw("0 0.0000e+00 0 0 0\n")
                    fw("##\n")
                    fw("## 3. Main Parameters III\n")
                    fw("0 0 0 0\n")
                    fw("##\n")
                    fw("## 4. Main Parameters IV\n")
                    fw("%s %se+00 1.0 0.0e+00 0.0 e+00 0.0e+00 0.0e+00\n" % (speedOfSound, densityOfMedium))
                    fw("##\n")
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
                    fw("# SYMMETRY\n")
                    fw("# 0 0 0\n")
                    fw("# 0.0000e+00 0.0000e+00 0.0000e+00\n")
                    fw("##\n")
                    if reciprocity:
                        if cpusAndCores[cpu-1][core-1]==1 and ear!='Right ear':
                            tmpEar='Left ear'
                        else:
                            tmpEar='Right ear'

                        fw("BOUNDARY\n")
                        for ii in range(len(obj_data.polygons[:])):
                            if obj.material_slots[obj_data.polygons[ii].material_index].name == obj.material_slots[tmpEar].name:
                                fw("ELEM %i TO %i VELO 0.1 IMAG 0.0\n" % (ii, ii))
                        fw("RETU\n")
                    else:
                        fw("BOUNDARY\n")
                        fw("# ELEM 0 TO 0 VELO 0.1 IMAG 0.0\n")
                        fw("RETU\n")
                    fw("##\n")
                    fw("# PLANE WAVES\n")
                    fw("# 0 0.0000e+00 -1.0000e+00 0.0000e+00 1.0000e-6 -1 0.0000e+00 -1\n")
                    fw("##\n")
                    if reciprocity:
                        fw("# POINT SOURCES\n")
                        if cpusAndCores[cpu-1][core-1] == 1:
                            fw("# 0 0.0 0.101 0.0 0.1 -1 0.0 -1\n")
                        if cpusAndCores[cpu-1][core-1] == 2:
                            fw("# 0 0.0 -0.101 0.0 0.1 -1 0.0 -1\n")
                    else:
                        fw("POINT SOURCES\n")
                        fw("0 %s %s %s 0.1 -1 0.0 -1\n" % (sourceXPosition, sourceYPosition, sourceZPosition))
                    fw("##\n")
                    fw("# CURVES\n")
                    fw("# Frequency Factor 0.0\n")
                    fw("##\n")
                    fw("POST PROCESS\n")
                    fw("##\n")
                    fw("END\n")
                    file.close()

        for obj in bpy.context.scene.objects[:]:
            bpy.data.objects[obj.name].select_set(False)

        return {'FINISHED'}


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
    evalGrids =[]
    for eg in grids:
        # complete grids given by name to default path
        if not os.path.sep in eg:
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


def _write_info_txt(evalGridPaths, title, ear, filepath1, cpuFirst, cpuLast,
                    numCoresAvailable, frequencyStepsPerCore,
                    frequencies, freqs, frequencyStepSize, numFrequencySteps):
    file = open(os.path.join(filepath1, "Info.txt"), "w",
                encoding="utf8", newline="\n")
    fw = file.write
    fw("#####################################\n")
    fw("######## General information ########\n")
    fw("#####################################\n\n")
    fw("Program: Mesh2HRTF\n")
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
    for core in range(1, 9):
        for cpu in range(1, 11):
            fw("CPU_%d (Core %d):\n" % (cpu, core))
            for ii in range(0, len(frequencies[cpu-1][core-1])):
                fw("    %f\n" % frequencies[cpu-1][core-1][ii])
        fw("\n")
    file.close


def _write_output2HRTF_m(filepath1, version,
                         reciprocity, ear, unitFactor, reference,
                         speedOfSound, densityOfMedium,
                         cpusAndCores, maxCPUs, maxCores,
                         sourceXPosition, sourceYPosition, sourceZPosition):

    # file handling
    file = open(os.path.join(filepath1, "Output2HRTF.m"), "w",
                encoding="utf8", newline="\n")
    fw = file.write

    # header
    fw("% Collect the data simulated by NumCalc and save to project folder.\n")
    fw(f"% Mesh2HRTF Version {version}\n")
    fw("close all; clear\n\n")

    # general information
    fw("% Constants\n")
    fw("speedOfSound = " + speedOfSound + "; % [m/s]\n")
    fw("densityOfAir = " + densityOfMedium + "; % [kg/m^3]\n\n")

    fw("% Reciprocal simulation\n")
    fw("reciprocity = ")
    if reciprocity:
        fw("1")
    else:
        fw("0")
    fw(";\n\n")

    fw("% Reference to a point source in the origin\n")
    fw("% accoring to the classical HRTF definition\n")
    fw("% (https://doi.org/10.1016/0003-682X(92)90046-U)\n")
    fw("reference = ")
    if reference:
        fw("true;\n\n")
    else:
        fw("false;\n\n")

    # add information about the reciever/point source
    if reciprocity:

        # get the receiver/ear centers and areas
        obj = bpy.data.objects['Reference']
        obj_data = obj.data
        earCenter, earArea = _calculateReceiverProperties(
            obj,obj_data,unitFactor)

        # write left ear data
        if ear=='Left ear' or ear=='Both ears':
            fw("% left ear / receiver\n")
            fw("receiverCenter(1,1:3)=[%f %f %f];\n" % (earCenter[0][0],
                                                        earCenter[0][1],
                                                        earCenter[0][2]))
            fw("receiverArea(1,1)    =%g;\n" % earArea[0])

        # write right ear data
        if ear=='Right ear' or ear=='Both ears':
            if ear=='Right ear':
                nn = 1
            if ear=='Both ears':
                nn = 2

            fw("% right ear / receiver\n")
            fw("receiverCenter(%d,1:3) = [%f %f %f];\n" % (nn,
                                                           earCenter[1][0],
                                                           earCenter[1][1],
                                                           earCenter[1][2]))
            fw("receiverArea(%d,1) = %g;\n" % (nn, earArea[1]))

        fw("\n")
    else:

        fw("% point source / receiver\n")
        fw("receiverCenter(1,1:3) = [%s %s %s];\n" % (sourceXPosition,
                                                      sourceYPosition,
                                                      sourceZPosition))
        fw("receiverArea(1,1)     = 1;\n")
    fw("\n")

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
    fw("Output2HRTF_Main(cpusAndCores,reciprocity, ...\n")
    fw("                 receiverCenter,receiverArea,reference, ...\n")
    fw("                 speedOfSound,densityOfAir);\n")
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
                or obj.material_slots[obj_data.polygons[ii].material_index].name \
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
                    value = obj_data.vertices[obj_data.polygons[ii].vertices[vertex]].co[coord]
                    if value < earCenter[ear][coord][0]:
                        earCenter[ear][coord][0] = value
                    if value > earCenter[ear][coord][1]:
                        earCenter[ear][coord][1] = value

            # to calculate the polygons (triangles) area first calculate the side lengths using euclidean distance and then use Heron´s formula to calculate the area
            # side_a = corner 0 to corner 1
            side_a = math.sqrt(((
                obj_data.vertices[obj_data.polygons[ii].vertices[0]].co[0] - \
                obj_data.vertices[obj_data.polygons[ii].vertices[1]].co[0])**2) \
                * unitFactor**2 + ((
                obj_data.vertices[obj_data.polygons[ii].vertices[0]].co[1] - \
                obj_data.vertices[obj_data.polygons[ii].vertices[1]].co[1])**2) \
                * unitFactor**2 + ((
                obj_data.vertices[obj_data.polygons[ii].vertices[0]].co[2] - \
                obj_data.vertices[obj_data.polygons[ii].vertices[1]].co[2])**2) \
                * unitFactor**2)

            # side_b = corner 1 to corner 2
            side_b = math.sqrt(((
                obj_data.vertices[obj_data.polygons[ii].vertices[1]].co[0] - \
                obj_data.vertices[obj_data.polygons[ii].vertices[2]].co[0])**2) \
                * unitFactor**2 + ((
                obj_data.vertices[obj_data.polygons[ii].vertices[1]].co[1] - \
                obj_data.vertices[obj_data.polygons[ii].vertices[2]].co[1])**2) \
                * unitFactor**2 + ((
                obj_data.vertices[obj_data.polygons[ii].vertices[1]].co[2] -
                obj_data.vertices[obj_data.polygons[ii].vertices[2]].co[2])**2) \
                * unitFactor**2)

            # side_c = corner 2 to corner 0
            side_c = math.sqrt(((
                obj_data.vertices[obj_data.polygons[ii].vertices[2]].co[0] -
                obj_data.vertices[obj_data.polygons[ii].vertices[0]].co[0])**2) \
                * unitFactor**2 + ((
                obj_data.vertices[obj_data.polygons[ii].vertices[2]].co[1] - \
                obj_data.vertices[obj_data.polygons[ii].vertices[0]].co[1])**2) \
                * unitFactor**2 + ((
                obj_data.vertices[obj_data.polygons[ii].vertices[2]].co[2] - \
                obj_data.vertices[obj_data.polygons[ii].vertices[0]].co[2])**2) \
                *unitFactor**2)

            # increment area using Heron´s formula
            earArea[ear] += 0.25 * math.sqrt(
                             (side_a+side_b+side_c)
                            *(-side_a+side_b+side_c)
                            *(side_a-side_b+side_c)
                            *(side_a+side_b-side_c))

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
                            numEars,
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
        Step size between sucsessive frequncies (written to Info.txt). This is
        returned because it might be zero during the function call.
    numFrequencySteps: int
        Number of frequncies to be simulated (written to Info.txt). This is
        returned because it might be zero during the function call.

    """

    # number of CPUs used
    numCPUs = cpuLast-cpuFirst+1

    # number of cores per ear
    numCoresUsedPerEar = numCPUs*numCoresPerCPU//numEars
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
    f = [ff*frequencyStepSize+minFrequency for ff in range(frequencySteps[0])]

    # remove 0 Hz if included in the list
    if f[0] == 0:
        f.pop(0)
        frequencySteps = (frequencySteps[0]-1, 0)
        print('Warning: 0 Hz can not be calculated and was removed from the \
              list of frequencies.')

    if not len(f):
        raise ValueError("No frequncies to be calulated. \
                         Check the input parameters.")

    # check number of cores and frequencies
    if len(f) < numCoresUsedPerEar:
        raise Exception("More cores than frequencies, i.e., \
                        numCPUs*numCoresPerCPU//numEars < numFrequencies.")

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
    for count in range(numCoresUsedPerEar*numEars):
        cpu, core = divmod(count + (cpuFirst-1)*numCoresPerCPU, numCoresPerCPU)
        cpusAndCores[cpu][core] = count//numCoresUsedPerEar + 1
        frequencies[cpu][core] = F[count % len(F)]
        # output for debugging
        # print(f"CPU {cpu+1:2d}, core {core+1}, \
        #       ear {count//numCoresUsedPerEar + 1}, \
        #       freqList {count%len(F)}")

    # meta data for Info.txt
    frequencyStepsPerCore = len(f)//numCoresUsedPerEar
    numCoresAvailable = numCoresUsedPerEar*numEars

    return cpusAndCores, frequencies, \
        frequencyStepsPerCore, numCoresAvailable, \
        f, frequencyStepSize, numFrequencySteps

# ----------------------- Blender add-on registration -------------------------
def menu_func_export(self, context):
    self.layout.operator(ExportMesh2HRTF.bl_idname, text="Mesh2HRTF")

def register():
    bpy.utils.register_class(ExportMesh2HRTF)
    bpy.types.TOPBAR_MT_file_export.append(menu_func_export)

def unregister():
    bpy.utils.unregister_class(ExportMesh2HRTF)
    bpy.types.TOPBAR_MT_file_export.remove(menu_func_export)
