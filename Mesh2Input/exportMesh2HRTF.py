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
# Author: Harald Ziegelwanger (Acoustics Research Institute, Austrian Academy of Sciences)
# Co-Authors: Fabian Brinkmann, Robert Pelzer (Audio Communication Group, Technical University Berlin)	

import os
import bpy
import datetime
import math
import shutil
from math import pi
from bpy.props import StringProperty, BoolProperty, EnumProperty, IntProperty
from bpy_extras.io_utils import ExportHelper

bl_info = {
    "name": "Mesh2HRTF input format",
    "author": "Harald Ziegelwanger",
    "version": (0, 1, 3),
    "blender": (2, 76),
    "location": "File > Import-Export",
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
    filter_glob = StringProperty(default="", options={'HIDDEN'})

    title = StringProperty(
        name="Title",
        description="Title",
        default="Head-Related Transfer Functions",
        )
    frequencyStepSize = IntProperty(
        name="Freq step",
        description="Lowest frequency and frequency step-size",
        default=100,
        min=10,
        max=24000,
        )
    maxFrequency = IntProperty(
        name="Freq max",
        description="Highest evaluated frequency",
        default=20000,
        min=10,
        max=24000,
        )
    cpuFirst = IntProperty(
        name="CPU (first)",
        description="First 'CPU' used",
        default=1,
        min=1,
        max=100,
        )
    cpuLast = IntProperty(
        name="CPU (last)",
        description="Last 'CPU' used",
        default=10,
        min=1,
        max=100,
        )
    numCoresPerCPU = IntProperty(
        name="Num. of used cores",
        description="Number of used cores per CPU",
        default=8,
        min=1,
        max=8,
        )
    pictures = BoolProperty(
        name="Pictures",
        description="Render pictures",
        default=True,
        )
    ear = EnumProperty(
        name="Ear",
        description="Selected ear",
        items=[('Left ear', 'left', 'Left ear'),
               ('Right ear', 'right', 'Right ear'),
               ('Both ears', 'both', 'Both ears'),
               ('None', 'none', 'None')],
        default='Both ears',
        )
    evaluationGrid1 = EnumProperty(
        name="Ev.Grid 1",
        description="Selected evaluation grid",
        items=[('21_NF', 'NF', 'NF HYPER (N=46)'),
               ('22_FF', 'FF', 'FF HYPER (N=46)'),
               ('3_ARI', 'ARI', 'Same as acoustical ARI HRTF measurement (1550)'),
               ('4_Low_ICO', 'Low ICO', 'Low spatial resolution (ICO ~2000)'),
               ('4_Low_UV', 'Low UV', 'Low spatial resolution (UV ~5000)'),
               ('5_High_ICO', 'High ICO', 'High spatial resolution (2°x2°)'),
               ('5_High_UV', 'High UV', 'High spatial resolution (2°x2°)'),
               ('7_SPlane', 'SPlane', 'Sagittal plane (???)'),
               ('8_HPlane', 'HPlane', 'Horizontal plane (???)'),
               ('9_FPlane', 'FPlane', 'Frontal plane (???)'),
               ('Custom', 'Custom', 'User defined evaluation grid'),
               ('None', 'None', 'None')],
        default='3_ARI',
        )
    evaluationGrid2 = EnumProperty(
        name="Ev.Grid 2",
        description="Selected evaluation grid",
        items=[('21_NF', 'NF', 'NF HYPER (N=46)'),
               ('22_FF', 'FF', 'FF HYPER (N=46)'),
               ('3_ARI', 'ARI', 'Same as acoustical ARI HRTF measurement (1550)'),
               ('4_Low_ICO', 'Low ICO', 'Low spatial resolution (ICO ~2000)'),
               ('4_Low_UV', 'Low UV', 'Low spatial resolution (UV ~5000)'),
               ('5_High_ICO', 'High ICO', 'High spatial resolution (2°x2°)'),
               ('5_High_UV', 'High UV', 'High spatial resolution (2°x2°)'),
               ('7_SPlane', 'SPlane', 'Sagittal plane (???)'),
               ('8_HPlane', 'HPlane', 'Horizontal plane (???)'),
               ('9_FPlane', 'FPlane', 'Frontal plane (???)'),
               ('Custom', 'Custom', 'User defined evaluation grid'),
               ('None', 'None', 'None')],
        default='None',
        )
    evaluationGrid3 = EnumProperty(
        name="Ev.Grid 3",
        description="Selected evaluation grid",
        items=[('21_NF', 'NF', 'NF HYPER (N=46)'),
               ('22_FF', 'FF', 'FF HYPER (N=46)'),
               ('3_ARI', 'ARI', 'Same as acoustical ARI HRTF measurement (1550)'),
               ('4_Low_ICO', 'Low ICO', 'Low spatial resolution (ICO ~2000)'),
               ('4_Low_UV', 'Low UV', 'Low spatial resolution (UV ~5000)'),
               ('5_High_ICO', 'High ICO', 'High spatial resolution (2°x2°)'),
               ('5_High_UV', 'High UV', 'High spatial resolution (2°x2°)'),
               ('7_SPlane', 'SPlane', 'Sagittal plane (???)'),
               ('8_HPlane', 'HPlane', 'Horizontal plane (???)'),
               ('9_FPlane', 'FPlane', 'Frontal plane (???)'),
               ('Custom', 'Custom', 'User defined evaluation grid'),
               ('None', 'None', 'None')],
        default='None',
        )
    evaluationGrid4 = EnumProperty(
        name="Ev.Grid 4",
        description="Selected evaluation grid",
        items=[('21_NF', 'NF', 'NF HYPER (N=46)'),
               ('22_FF', 'FF', 'FF HYPER (N=46)'),
               ('3_ARI', 'ARI', 'Same as acoustical ARI HRTF measurement (1550)'),
               ('4_Low_ICO', 'Low ICO', 'Low spatial resolution (ICO ~2000)'),
               ('4_Low_UV', 'Low UV', 'Low spatial resolution (UV ~5000)'),
               ('5_High_ICO', 'High ICO', 'High spatial resolution (2°x2°)'),
               ('5_High_UV', 'High UV', 'High spatial resolution (2°x2°)'),
               ('7_SPlane', 'SPlane', 'Sagittal plane (???)'),
               ('8_HPlane', 'HPlane', 'Horizontal plane (???)'),
               ('9_FPlane', 'FPlane', 'Frontal plane (???)'),
               ('Custom', 'Custom', 'User defined evaluation grid'),
               ('None', 'None', 'None')],
        default='None',
        )
    evaluationGrid5 = EnumProperty(
        name="Ev.Grid 5",
        description="Selected evaluation grid",
        items=[('21_NF', 'NF', 'NF HYPER (N=46)'),
               ('22_FF', 'FF', 'FF HYPER (N=46)'),
               ('3_ARI', 'ARI', 'Same as acoustical ARI HRTF measurement (1550)'),
               ('4_Low_ICO', 'Low ICO', 'Low spatial resolution (ICO ~2000)'),
               ('4_Low_UV', 'Low UV', 'Low spatial resolution (UV ~5000)'),
               ('5_High_ICO', 'High ICO', 'High spatial resolution (2°x2°)'),
               ('5_High_UV', 'High UV', 'High spatial resolution (2°x2°)'),
               ('7_SPlane', 'SPlane', 'Sagittal plane (???)'),
               ('8_HPlane', 'HPlane', 'Horizontal plane (???)'),
               ('9_FPlane', 'FPlane', 'Frontal plane (???)'),
               ('Custom', 'Custom', 'User defined evaluation grid'),
               ('None', 'None', 'None')],
        default='None',
        )
    method = EnumProperty(
        name="Method",
        description="Choose the calculation method",
        items=[('0', 'BEM', 'Traditional BEM'),
               ('1', 'SL-FMM BEM', 'Singlelevel fast-multipole method'),
               ('4', 'ML-FMM BEM', 'Multilevel fast-multipole method')],
        default='4',
        )
    reciprocity = BoolProperty(
        name="Recip.",
        description="Calculation with reciprocity",
        default=True,
        )
    sourceXPosition = StringProperty(
        name="Source (x)",
        description="Source Position (X-Coordinate)",
        default="0",
        )
    sourceYPosition = StringProperty(
        name="Source (y)",
        description="Source Position (Y-Coordinate)",
        default="101",
        )
    sourceZPosition = StringProperty(
        name="Source (z)",
        description="Source Position (Z-Coordinate)",
        default="0",
        )
    speedOfSound = StringProperty(
        name="c (m/s)",
        description="Speed of sound (m/s)",
        default="346.18",
        )
    densityOfMedium = StringProperty(
        name="rho ()",
        description="Density of air (kg/m^3)",
        default="1.1839",
        )
    unit = EnumProperty(
        name="Unit",
        description="Unit of the object",
        items=[('m', 'm', 'Meter'), ('mm', 'mm', 'Millimeter')],
        default='mm',
        )
    frequencyDependency = BoolProperty(
        name="Freq.-dep.",
        description="Use frequency-dependent meshes",
        default=False,
        )
    nearFieldCalculation = BoolProperty(
        name="NF-Calc.",
        description="Calculate near-field HRTFs",
        default=False,
        )
    programPath = StringProperty(
        name="Mesh2HRTF-path",
        description="Path to mesh2HRTF",
        default=r"C:\Users\jkhan\Documents\Mesh2HRTF - Kopie\trunk",
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

        row = layout.row()
        row.prop(self, "title")
        row = layout.row()
        row.prop(self, "ear")
        row = layout.row()
        row.prop(self, "pictures")
        layout.label("Point Source:")
        row = layout.row()
        row.prop(self, "sourceXPosition")
        row = layout.row()
        row.prop(self, "sourceYPosition")
        row = layout.row()
        row.prop(self, "sourceZPosition")
        row = layout.row()
        row.prop(self, "reciprocity")
        layout.label("Constants:")
        row = layout.row()
        row.prop(self, "speedOfSound")
        row = layout.row()
        row.prop(self, "densityOfMedium")
        layout.label("ObjectMeshes:")
        row = layout.row()
        row.prop(self, "unit")
        row = layout.row()
        layout.label("Evaluation Grids:")
        row = layout.row()
        row.prop(self, "evaluationGrid1")
        row = layout.row()
        row.prop(self, "evaluationGrid2")
        row = layout.row()
        row.prop(self, "evaluationGrid3")
        row = layout.row()
        row.prop(self, "evaluationGrid4")
        row = layout.row()
        row.prop(self, "evaluationGrid5")
        row = layout.row()
        row.prop(self, "nearFieldCalculation")
        layout.label("Frequencies:")
        row = layout.row()
        row.prop(self, "frequencyStepSize")
        row = layout.row()
        row.prop(self, "maxFrequency")
        row = layout.row()
        row.prop(self, "frequencyDependency")
        row = layout.row()
        row.prop(self, "method")
        layout.label("Cluster:")
        row = layout.row()
        row.prop(self, "cpuFirst")
        row = layout.row()
        row.prop(self, "cpuLast")
        row = layout.row()
        row.prop(self, "numCoresPerCPU")
        layout.label("Mesh2HRTF:")
        row = layout.row()
        row.prop(self, "programPath")

    def save(operator,
             context,
             filepath="",
             title="head-related transfer functions",
             frequencyStepSize=100,
             maxFrequency=20000,
             cpuFirst=1,
             cpuLast=10,
             numCoresPerCPU=8,
             pictures=True,
             ear='Both ears',
             evaluationGrid1='3_ARI',
             evaluationGrid2='None',
             evaluationGrid3='None',
             evaluationGrid4='None',
             evaluationGrid5='None',
             method='4',
             reciprocity=True,
             sourceXPosition='0',
             sourceYPosition='101',
             sourceZPosition='0',
             speedOfSound='346.18',
             densityOfMedium='1.1839',
             unit='mm',
             frequencyDependency=False,
             nearFieldCalculation=False,
             programPath="",
             ):

        def rvec3d(v):
            return round(v[0], 6), round(v[1], 6), round(v[2], 6)

        def rvec2d(v):
            return round(v[0], 6), round(v[1], 6)

# ----------------------- Initialize constants ---------------------------------
        bpy.ops.object.transform_apply(location=True)
        bpy.ops.object.transform_apply(rotation=True)
        bpy.ops.object.transform_apply(scale=True)
        cam = bpy.data.objects['Camera']
        camradius = 400
        cam.location = (0, camradius, 0)
        cam.rotation_euler = (pi/2, 0, pi)
#         bpy.data.scenes['Scene'].camera = cam
        cam.data.clip_end = 0.1
        cam.data.clip_end = 1000
        lamp = bpy.data.objects['Lamp']
        lampradius = 300
        lamp.location = (0, lampradius, 0)
        bpy.data.lamps['Lamp'].energy = 800
        bpy.data.lamps['Lamp'].distance = 100
        renderloc = [[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0.707, 0.707, 0], [-0.707, 0.707, 0], [0.707, -0.707, 0], [-0.707, -0.707, 0]]  # , [0.707, 0, 0.707], [-0.707, 0, 0.707], [0.707, 0,-0.707], [-0.707, 0,-0.707], [0, 0.707, 0.707], [0,-0.707, 0.707], [0, 0.707,-0.707], [0,-0.707,-0.707]]
        renderrot = [[pi/2, 0, pi/2], [pi/2, 0, 3*pi/2], [pi/2, 0, pi], [3/2*pi, pi, pi], [pi/2, 0, 3/4*pi], [pi/2, 0, 5/4*pi], [pi/2, 0, pi/4], [pi/2, 0, -pi/4]]  # , [pi/4, 0, pi/2], [5/4*pi, pi, pi/2], [3/4*pi, 0, pi/2], [0, 4*5/pi, 0], [pi/4, 0, pi], [pi/4, 0, 0], [3/4*pi, 0, pi], [3/4*pi, 0, 0]]
        rendernam = [[0, 0], [180, 0], [90, 0], [270, 0], [45, 0], [135, 0], [315, 0], [225, 0]]

        bpy.data.scenes['Scene'].render.pixel_aspect_x = 1
        bpy.data.scenes['Scene'].render.pixel_aspect_y = 1
        bpy.data.scenes['Scene'].render.resolution_x = 1440
        bpy.data.scenes['Scene'].render.resolution_y = 1920
        
        tmp = open("%s/VERSION" % programPath)
        version = tmp.readline()

        objects = ([])

        (filepath1, filename1) = os.path.split(filepath)
        filename1 = "NC.inp"

        temp = ("%s/ObjectMeshes/" % filepath1)
        if not os.path.exists(temp):
            os.mkdir(temp)

        temp = ("%s/EvaluationGrids/" % filepath1)
        if not os.path.exists(temp):
            os.mkdir(temp)

        temp = ("%s/NumCalc/" % filepath1)
        if not os.path.exists(temp):
            os.mkdir(temp)

        numCPUs = cpuLast-cpuFirst+1

        numEars = 1
        if ear == 'Both ears':
            numEars = 2

        unitFactor = 1
        if unit == 'mm':
            unitFactor = 0.001

        lowFrequency = 0
        lowFrequencyCores = 0
        if not frequencyDependency:
            obj = bpy.data.objects["Reference"]
            obj.hide_render = False
            obj_data = obj.data
            if len(obj_data.vertices[:]) > 40000 and numCPUs/numEars > 1:
                lowFrequency = frequencyStepSize*10
                lowFrequencyCores = 1
            else:
                lowFrequency = 0
                lowFrequencyCores = 0

        evaluationGridPath = ("%s/Mesh2Input/EvaluationGrids" % programPath)

# ------------------------ Write object data -----------------------------------
        for obj in bpy.context.scene.objects[:]:
            if obj.type == 'MESH' and not obj.name == 'User':
                bpy.context.scene.objects.active = obj
                bpy.ops.object.transform_apply(location=True)
                bpy.ops.object.transform_apply(rotation=True)
                bpy.ops.object.transform_apply(scale=True)
                obj = context.active_object
                obj.hide_render = False
                obj_data = obj.data

                temp = ("%s/ObjectMeshes/%s/" % (filepath1, obj.name))
                if not os.path.exists(temp):
                    os.mkdir(temp)

                file = open(("%s/ObjectMeshes/%s/Nodes.txt" % (filepath1, obj.name)), "w", encoding="utf8", newline="\n")
                fw = file.write
                fw("%i\n" % len(obj_data.vertices[:]))
                for ii in range(len(obj_data.vertices[:])):
                    fw("%i " % ii)
                    fw("%.6f %.6f %.6f\n" % (obj_data.vertices[ii].co[0]*unitFactor, obj_data.vertices[ii].co[1]*unitFactor, obj_data.vertices[ii].co[2]*unitFactor))
                file.close

                file = open(("%s/ObjectMeshes/%s/Elements.txt" % (filepath1, obj.name)), "w", encoding="utf8", newline="\n")
                fw = file.write
                fw("%i\n" % len(obj_data.polygons[:]))

                for ii in range(len(obj_data.polygons[:])):
                    fw("%i " % ii)
                    if len(obj_data.polygons[ii].vertices[:]) == 3:
                        fw("%d %d %d 0 0 0\n" % tuple(obj_data.polygons[ii].vertices[:]))
                    else:
                        fw("%d %d %d %d 0 0 0\n" % tuple(obj_data.polygons[ii].vertices[:]))
                file.close
                
                objects.append(obj.name)

        maxObjectFrequency = ([])
        for ii in range(0, len(objects)):
            if not objects[ii] == 'Reference' and not objects[ii] == 'User':
                try:
                    if maxObjectFrequency.count(int(objects[ii][1:len(objects[ii]):1])) == 0:
                        maxObjectFrequency.append(int(objects[ii][1:len(objects[ii]):1]))
                except:
                    print('No maximum object frequency found.\nPlease change object names to L{maxobjfq}/R{maxobjfq} e.g. L20000 or R20000.')
        maxObjectFrequency.sort()
        
        bpy.ops.wm.save_as_mainfile(filepath=("%s/3d Model.blend" % filepath1), check_existing=False, filter_blender=True, filter_image=False, filter_movie=False, filter_python=False, filter_font=False, filter_sound=False, filter_text=False, filter_btx=False, filter_collada=False, filter_folder=True, filemode=8, compress=False, relative_remap=True, copy=False)

# ------------------------ Write evaluation grid data --------------------------
        if nearFieldCalculation:
            radius = 0.0
            for ii in range(len(obj_data.vertices[:])):
                tmp = math.sqrt(math.pow(obj_data.vertices[ii].co[0], 2)+math.pow(obj_data.vertices[ii].co[1], 2)+math.pow(obj_data.vertices[ii].co[2], 2))
                if tmp > radius:
                    radius = tmp
            # bpy.ops.mesh.primitive_ico_sphere_add(subdivisions=5, size=radius+5, location=(0,0,0), rotation=(0, 0, 0))
            bpy.context.scene.cursor_location = (0.0, 0.0, 0.0)
            bpy.ops.mesh.primitive_hyper_add(orderN=46, size=(radius*unitFactor)+0.005)
            NFGrid = bpy.context.active_object
            NFGrid.name = 'NFGrid'
            NFGrid_data = NFGrid.data
            temp = ("%s/EvaluationGrids/NF_Sphere" % (filepath1))
            if not os.path.exists(temp):
                os.mkdir(temp)
           
            file = open(("%s/EvaluationGrids/NF_Sphere/Nodes.txt" % (filepath1)), "w", encoding="utf8", newline="\n")
            fw = file.write
            fw("%i\n" % len(NFGrid_data.vertices[:]))
            for ii in range(len(NFGrid_data.vertices[:])):
                fw("%i " % (ii+200000))
                fw("%.6f %.6f %.6f\n" % (NFGrid_data.vertices[ii].co[0]*unitFactor, NFGrid_data.vertices[ii].co[1]*unitFactor, NFGrid_data.vertices[ii].co[2]*unitFactor))
            file.close

            file = open(("%s/EvaluationGrids/NF_Sphere/Elements.txt" % (filepath1)), "w", encoding="utf8", newline="\n")
            fw = file.write
            fw("%i\n" % len(NFGrid_data.polygons[:]))
            for ii in range(len(NFGrid_data.polygons[:])):
                fw("%i " % (ii+200000))
                fw("%d %d %d" % (NFGrid_data.polygons[ii].vertices[0]+200000, NFGrid_data.polygons[ii].vertices[1]+200000, NFGrid_data.polygons[ii].vertices[2]+200000))
                fw(" 2 0 1\n")
            file.close

        else:
            if not evaluationGrid1 == 'None':
                if not evaluationGrid1 == 'User':
                    temp = ("%s/EvaluationGrids/%s" % (filepath1, evaluationGrid1))
                    if not os.path.exists(temp):
                        os.mkdir(temp)

                    shutil.copyfile(("%s/%s/Nodes.txt" % (evaluationGridPath, evaluationGrid1)), ("%s/EvaluationGrids/%s/Nodes.txt" % (filepath1, evaluationGrid1)))
                    shutil.copyfile(("%s/%s/Elements.txt" % (evaluationGridPath, evaluationGrid1)), ("%s/EvaluationGrids/%s/Elements.txt" % (filepath1, evaluationGrid1)))
                else:
                    obj = bpy.data.objects['User']
                    if obj.type == 'MESH':
                        bpy.context.scene.objects.active = obj
                        bpy.ops.object.transform_apply(location=True)
                        bpy.ops.object.transform_apply(rotation=True)
                        bpy.ops.object.transform_apply(scale=True)
                        obj = context.active_object
                        obj.hide_render = False
                        obj_data = obj.data

                        temp = ("%s/EvaluationGrids/User/" % filepath1)
                        if not os.path.exists(temp):
                            os.mkdir(temp)
              
                        file = open(("%s/EvaluationGrids/User/Nodes.txt" % filepath1), "w", encoding="utf8", newline="\n")
                        fw = file.write
                        fw("%i\n" % len(obj_data.vertices[:]))
                        for ii in range(len(obj_data.vertices[:])):
                            fw("%i " % (ii+350000))
                            fw("%.6f %.6f %.6f\n" % (obj_data.vertices[ii].co[0]*unitFactor, obj_data.vertices[ii].co[1]*unitFactor, obj_data.vertices[ii].co[2]*unitFactor))
                        file.close
            
                        file = open(("%s/EvaluationGrids/User/Elements.txt" % filepath1), "w", encoding="utf8", newline="\n")
                        fw = file.write
                        fw("%i\n" % len(obj_data.polygons[:]))
                        if len(obj_data.polygons[0].vertices[:]) == 3:
                            for ii in range(len(obj_data.polygons[:])):
                                fw("%i " % (ii+350000))
                                fw("%d %d %d 2 0 1\n" % (obj_data.polygons[ii].vertices[0]+350000, obj_data.polygons[ii].vertices[1]+350000, obj_data.polygons[ii].vertices[2]+350000))
                        else:
                            for ii in range(len(obj_data.polygons[:])):
                                fw("%i " % (ii+350000))
                                fw("%d %d %d %d 2 0 1\n" % (obj_data.polygons[ii].vertices[0]+350000, obj_data.polygons[ii].vertices[1]+350000, obj_data.polygons[ii].vertices[2]+350000, obj_data.polygons[ii].vertices[3]+350000))
                        file.close

            if not evaluationGrid2 == 'None':
                temp = ("%s/EvaluationGrids/%s" % (filepath1, evaluationGrid2))
                if not os.path.exists(temp):
                    os.mkdir(temp)

                shutil.copyfile(("%s/%s/Nodes.txt" % (evaluationGridPath, evaluationGrid2)), ("%s/EvaluationGrids/%s/Nodes.txt" % (filepath1, evaluationGrid2)))
                shutil.copyfile(("%s/%s/Elements.txt" % (evaluationGridPath, evaluationGrid2)), ("%s/EvaluationGrids/%s/Elements.txt" % (filepath1, evaluationGrid2)))

            if not evaluationGrid3 == 'None':
                temp = ("%s/EvaluationGrids/%s" % (filepath1, evaluationGrid3))
                if not os.path.exists(temp):
                    os.mkdir(temp)

                shutil.copyfile(("%s/%s/Nodes.txt" % (evaluationGridPath, evaluationGrid3)), ("%s/EvaluationGrids/%s/Nodes.txt" % (filepath1, evaluationGrid3)))
                shutil.copyfile(("%s/%s/Elements.txt" % (evaluationGridPath, evaluationGrid3)), ("%s/EvaluationGrids/%s/Elements.txt" % (filepath1, evaluationGrid3)))

            if not evaluationGrid4 == 'None':
                temp = ("%s/EvaluationGrids/%s" % (filepath1, evaluationGrid4))
                if not os.path.exists(temp):
                    os.mkdir(temp)

                shutil.copyfile(("%s/%s/Nodes.txt" % (evaluationGridPath, evaluationGrid4)), ("%s/EvaluationGrids/%s/Nodes.txt" % (filepath1, evaluationGrid4)))
                shutil.copyfile(("%s/%s/Elements.txt" % (evaluationGridPath, evaluationGrid4)), ("%s/EvaluationGrids/%s/Elements.txt" % (filepath1, evaluationGrid4)))

            if not evaluationGrid5 == 'None':
                temp = ("%s/EvaluationGrids/%s" % (filepath1, evaluationGrid5))
                if not os.path.exists(temp):
                    os.mkdir(temp)

                shutil.copyfile(("%s/%s/Nodes.txt" % (evaluationGridPath, evaluationGrid5)), ("%s/EvaluationGrids/%s/Nodes.txt" % (filepath1, evaluationGrid5)))
                shutil.copyfile(("%s/%s/Elements.txt" % (evaluationGridPath, evaluationGrid5)), ("%s/EvaluationGrids/%s/Elements.txt" % (filepath1, evaluationGrid5)))

# ------------------------ Calculate frequency information ---------------------
        frequencySteps = divmod(maxFrequency-lowFrequency, frequencyStepSize)
        if not frequencySteps[1] == 0:
            raise Exception("Error, frequencyStepSize is not a divisor of maxFrequency-lowFrequency")
            
        numCoresAvailable = (cpuLast-cpuFirst+1-lowFrequencyCores)*numCoresPerCPU
        numCoresUsedPerEar = int((cpuLast-cpuFirst+1-lowFrequencyCores*numEars)*numCoresPerCPU/numEars)
        frequencyStepsPerCore = divmod(frequencySteps[0], numCoresUsedPerEar)

        cpusAndCores = ([])
        tmp = ([])
        for core in range(1, 9):
            tmp.append(0)
        for cpu in range(1, 11):
            cpusAndCores.append(tmp[:])

        frequencies = ([])
        tmp = ([])
        for core in range(1, 9):
            tmp.append([])
        for cpu in range(1, 11):
            frequencies.append(tmp[:])
            
        if not frequencyDependency:
            coresteps = 0
            for tmpEar in range(1, numEars+1):
                count = 0
                for core in range(1, numCoresPerCPU+1):
                    if tmpEar == 2 and numCPUs == 1:
                        core += numCoresUsedPerEar
                    for cpu in range(cpuFirst+(int(numCPUs/2)*(tmpEar-1))*(numEars-1), cpuLast-(int(numCPUs/2)*(-(tmpEar-2)))*(numEars-1)+1):
                        tmp = ([])
                        if lowFrequencyCores > 0 and cpu == cpuFirst+(int(numCPUs/2)*(tmpEar-1))*(numEars-1):
                            if core < 3:
                                for ii in range(1, 6):
                                    tmp.append(frequencyStepSize*ii+lowFrequency/2*(core-1))
                        else:
                            count = count + 1
                            for ii in range(0, frequencyStepsPerCore[0]+1):
                                if ((frequencyStepSize*count+frequencyStepSize*numCoresUsedPerEar*ii)+lowFrequency) <= maxFrequency:
                                    tmp.append((frequencyStepSize*count+frequencyStepSize*numCoresUsedPerEar*ii)+lowFrequency)
                                else:
                                    break
                        frequencies[cpu-1][core-1] = tmp
                        if not frequencies[cpu-1][core-1] == ([]):
							#case: both ears on one cpu
                            if numCPUs == 1 and numEars == 2:
                                if numCoresPerCPU < 4:
                                    raise Exception('Please use at least 4 cores for calculation of 2 ears')
                                for temp_ear in range(0, numCoresUsedPerEar):
                                    if int((temp_ear+coresteps)/2) >= numCoresUsedPerEar:
                                        break
                                    else:
                                        cpusAndCores[cpu-1][temp_ear+coresteps] = tmpEar
                                coresteps += 1
							
                            else:
                                #general case
                                cpusAndCores[cpu-1][core-1] = tmpEar
				
                        if count == numCoresUsedPerEar:
                            break
                    if count == numCoresUsedPerEar:
                        break
        else:
            for tmpEar in range(1, numEars+1):
                countFreq = 0
                countCores = 0
                for core in range(1, numCoresPerCPU+1):
                    for cpu in range(cpuFirst+(int(numCPUs/2)*(tmpEar-1))*(numEars-1), cpuLast-(int(numCPUs/2)*(-(tmpEar-2)))*(numEars-1)+1):
                        countCores = countCores+1
                        tmp = ([])
                        if countCores <= frequencyStepsPerCore[1]:
                            tmpNumFrequencies = frequencyStepsPerCore[0]+1
                        else:
                            tmpNumFrequencies = frequencyStepsPerCore[0]
                        for ii in range(0, tmpNumFrequencies):
                            if (frequencyStepSize*countFreq) <= maxFrequency:
                                countFreq = countFreq + 1
                                tmp.append(frequencyStepSize*countFreq)
                            else:
                                break
                        frequencies[cpu-1][core-1] = tmp
                        cpusAndCores[cpu-1][core-1] = tmpEar

# ----------------------- Write general information ----------------------------
        file = open(("%s/Info.txt" % filepath1), "w", encoding="utf8", newline="\n")
        fw = file.write
        fw("#####################################\n")
        fw("######## General information ########\n")
        fw("#####################################\n\n")
        fw("Program: Mesh2HRTF\n")
        fw("Title: %s\n" % title)
        fw("Ear: %s\n" % ear)
        fw("Evaluation Grids:\n")
        fw("    %s\n" % evaluationGrid1)
        fw("    %s\n" % evaluationGrid2)
        fw("    %s\n" % evaluationGrid3)
        fw("    %s\n" % evaluationGrid4)
        fw("    %s\n\n" % evaluationGrid5)
        fw("#####################################\n")
        fw("####### Frequency information #######\n")
        fw("#####################################\n\n")
        fw("Highest evaluated Frequency: %d\n" % maxFrequency)
        fw("Frequency Stepsize: %d\n" % frequencyStepSize)
        fw("Frequency Steps: %d\n" % frequencySteps[0])
        fw("Frequency steps per Core: %d\n\n" % frequencyStepsPerCore[0])
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
                    fw("    %d\n" % frequencies[cpu-1][core-1][ii])
            fw("\n")
        file.close

# ----------------------- Write Output2HRTF.m function -------------------------
        file = open(("%s/Output2HRTF.m" % filepath1), "w", encoding="utf8", newline="\n")
        fw = file.write
        fw("close all\n")
        fw("clear\n")
        fw("\n")
        
        fw("cpusAndCores=[")
        for cpu in range(1, 11):
            for core in range(1, 9):
                fw("%i" % cpusAndCores[cpu-1][core-1])
                if core < 8:
                    fw(" ")
            if cpu < 10:
                fw("; ...\n")
        fw("];\n")
        fw("\n")
        
        fw("objectMeshes={")
        for cpu in range(1, 11):
            for core in range(1, 9):
                if not cpusAndCores[cpu-1][core-1] == 0:
                    if frequencyDependency:
                        if cpusAndCores[cpu-1][core-1] == 1:
                            tmpEar = "L"
                        if cpusAndCores[cpu-1][core-1] == 2:
                            tmpEar = "R"
                        tmpfmax = max(frequencies[cpu-1][core-1])
                        tmpfdiff = 24000
                        for ff in maxObjectFrequency:
                            if (ff-tmpfmax) >= 0 and (ff-tmpfmax) < tmpfdiff:
                                obj_name = ("%s%i" % (tmpEar, ff))
                                tmpfdiff = ff-tmpfmax
                    else:
                        obj_name = "Reference"
                else:
                    obj_name = ""
                fw("'%s'" % obj_name)
                if core < 8:
                    fw(" ")
            if cpu < 10:
                fw("; ...\n")
        fw("};\n")
        fw("\n")
                
        fw("reciprocity=")
        if reciprocity:
            fw("1")
        else:
            fw("0")
        fw(";\n")
        fw("\n")

        if reciprocity:
            obj = bpy.data.objects['Reference']
            obj_data = obj.data
            if ear=='Left ear' or ear=='Both ears':
                fw("receiverPositions(1,1:3)=[")
                for ii in range(len(obj_data.polygons[:])):
                    if obj.material_slots[obj_data.polygons[ii].material_index].name == obj.material_slots['Left ear'].name:
                        fw("%f " % ((obj_data.vertices[obj_data.polygons[ii].vertices[0]].co[0]+obj_data.vertices[obj_data.polygons[ii].vertices[1]].co[0]+obj_data.vertices[obj_data.polygons[ii].vertices[2]].co[0])/3*unitFactor))
                        fw("%f " % ((obj_data.vertices[obj_data.polygons[ii].vertices[0]].co[1]+obj_data.vertices[obj_data.polygons[ii].vertices[1]].co[1]+obj_data.vertices[obj_data.polygons[ii].vertices[2]].co[1])/3*unitFactor))
                        fw("%f];\n" % ((obj_data.vertices[obj_data.polygons[ii].vertices[0]].co[2]+obj_data.vertices[obj_data.polygons[ii].vertices[1]].co[2]+obj_data.vertices[obj_data.polygons[ii].vertices[2]].co[2])/3*unitFactor))
                        break
                        

                fw("microphone_area(1,1)=[")
                for ii in range(len(obj_data.polygons[:])):
                    if obj.material_slots[obj_data.polygons[ii].material_index].name == obj.material_slots['Left ear'].name:

                        ## to calculate the polygons (triangles) area first calculate the side lengths  using euclidean distance and then use Heron´s formula to calculate the area
                        #left_right_distance=math.sqrt((right[0]-left[0])**2 + (right[1]-left[1])**2 + (right[2]-left[2])**2)

                        # side_a = corner 0 to corner 1
                        side_a=math.sqrt(((obj_data.vertices[obj_data.polygons[ii].vertices[0]].co[0]-obj_data.vertices[obj_data.polygons[ii].vertices[1]].co[0])**2)*unitFactor**2 + ((obj_data.vertices[obj_data.polygons[ii].vertices[0]].co[1] - obj_data.vertices[obj_data.polygons[ii].vertices[1]].co[1])**2)*unitFactor**2 + ((obj_data.vertices[obj_data.polygons[ii].vertices[0]].co[2]-obj_data.vertices[obj_data.polygons[ii].vertices[1]].co[2])**2)*unitFactor**2)

                        # side_b = corner 1 to corner 2
                        side_b=math.sqrt(((obj_data.vertices[obj_data.polygons[ii].vertices[1]].co[0]-obj_data.vertices[obj_data.polygons[ii].vertices[2]].co[0])**2)*unitFactor**2 + ((obj_data.vertices[obj_data.polygons[ii].vertices[1]].co[1] - obj_data.vertices[obj_data.polygons[ii].vertices[2]].co[1])**2)*unitFactor**2 + ((obj_data.vertices[obj_data.polygons[ii].vertices[1]].co[2]-obj_data.vertices[obj_data.polygons[ii].vertices[2]].co[2])**2)*unitFactor**2)
                        
                        # side_c = corner 2 to corner 0
                        side_c=math.sqrt(((obj_data.vertices[obj_data.polygons[ii].vertices[2]].co[0]-obj_data.vertices[obj_data.polygons[ii].vertices[0]].co[0])**2)*unitFactor**2 + ((obj_data.vertices[obj_data.polygons[ii].vertices[2]].co[1] - obj_data.vertices[obj_data.polygons[ii].vertices[0]].co[1])**2)*unitFactor**2 + ((obj_data.vertices[obj_data.polygons[ii].vertices[2]].co[2]-obj_data.vertices[obj_data.polygons[ii].vertices[0]].co[2])**2)*unitFactor**2)
                        
                        # Heron´s formula
                        microphone_area=0.25 * math.sqrt((side_a+side_b+side_c)*(-side_a+side_b+side_c)*(side_a-side_b+side_c)*(side_a+side_b-side_c))
                        
                        
                        fw("%g];\n" % (microphone_area))
                        break
                    

                          
            if ear=='Right ear':
                nn = 1
            if ear=='Both ears':
                nn = 2
            

            if ear=='Right ear' or ear=='Both ears':
                fw("receiverPositions(%d,1:3)=[" %(nn))
               
                for ii in range(len(obj_data.polygons[:])):
                    if obj.material_slots[obj_data.polygons[ii].material_index].name == obj.material_slots['Right ear'].name:
                        fw("%f " % ((obj_data.vertices[obj_data.polygons[ii].vertices[0]].co[0]+obj_data.vertices[obj_data.polygons[ii].vertices[1]].co[0]+obj_data.vertices[obj_data.polygons[ii].vertices[2]].co[0])/3*unitFactor))
                        fw("%f " % ((obj_data.vertices[obj_data.polygons[ii].vertices[0]].co[1]+obj_data.vertices[obj_data.polygons[ii].vertices[1]].co[1]+obj_data.vertices[obj_data.polygons[ii].vertices[2]].co[1])/3*unitFactor))
                        fw("%f];\n" % ((obj_data.vertices[obj_data.polygons[ii].vertices[0]].co[2]+obj_data.vertices[obj_data.polygons[ii].vertices[1]].co[2]+obj_data.vertices[obj_data.polygons[ii].vertices[2]].co[2])/3*unitFactor))
                        break


                fw("microphone_area(%d,1)=[" %(nn))
                
                for ii in range(len(obj_data.polygons[:])):
                    if obj.material_slots[obj_data.polygons[ii].material_index].name == obj.material_slots['Right ear'].name:

                         ## to calculate the polygons (triangles) area first calculate the side lengths  using euclidean distance and then use Heron´s formula to calculate the area
                        #left_right_distance=math.sqrt((right[0]-left[0])**2 + (right[1]-left[1])**2 + (right[2]-left[2])**2)

                         # side_a = corner 0 to corner 1
                        side_a=math.sqrt(((obj_data.vertices[obj_data.polygons[ii].vertices[0]].co[0]-obj_data.vertices[obj_data.polygons[ii].vertices[1]].co[0])**2)*unitFactor**2 + ((obj_data.vertices[obj_data.polygons[ii].vertices[0]].co[1] - obj_data.vertices[obj_data.polygons[ii].vertices[1]].co[1])**2)*unitFactor**2 + ((obj_data.vertices[obj_data.polygons[ii].vertices[0]].co[2]-obj_data.vertices[obj_data.polygons[ii].vertices[1]].co[2])**2)*unitFactor**2)

                        # side_b = corner 1 to corner 2
                        side_b=math.sqrt(((obj_data.vertices[obj_data.polygons[ii].vertices[1]].co[0]-obj_data.vertices[obj_data.polygons[ii].vertices[2]].co[0])**2)*unitFactor**2 + ((obj_data.vertices[obj_data.polygons[ii].vertices[1]].co[1] - obj_data.vertices[obj_data.polygons[ii].vertices[2]].co[1])**2)*unitFactor**2 + ((obj_data.vertices[obj_data.polygons[ii].vertices[1]].co[2]-obj_data.vertices[obj_data.polygons[ii].vertices[2]].co[2])**2)*unitFactor**2)
                        
                        # side_c = corner 2 to corner 0
                        side_c=math.sqrt(((obj_data.vertices[obj_data.polygons[ii].vertices[2]].co[0]-obj_data.vertices[obj_data.polygons[ii].vertices[0]].co[0])**2)*unitFactor**2 + ((obj_data.vertices[obj_data.polygons[ii].vertices[2]].co[1] - obj_data.vertices[obj_data.polygons[ii].vertices[0]].co[1])**2)*unitFactor**2 + ((obj_data.vertices[obj_data.polygons[ii].vertices[2]].co[2]-obj_data.vertices[obj_data.polygons[ii].vertices[0]].co[2])**2)*unitFactor**2)
                        
                        # Heron´s formula
                        microphone_area=0.25 * math.sqrt((side_a+side_b+side_c)*(-side_a+side_b+side_c)*(side_a-side_b+side_c)*(side_a+side_b-side_c))
                        
                        
                        fw("%g];\n" % (microphone_area))
                        break
            fw("\n")
                
        fw("frequencyDependency=")
        if frequencyDependency:
            fw("1")
        else:
            fw("0")
        fw(";\n")
        fw("\n")
                
        fw("nearFieldCalculation=")
        if nearFieldCalculation:
            fw("1")
        else:
            fw("0")
        fw(";\n")
        fw("\n")

        fw("% Reference to a point source in the origin\n")
        fw("% accoring to the classical HRTF definition\n")
        fw("reference = false;\n\n")

        fw("Output2HRTF_Main(cpusAndCores,objectMeshes,reciprocity,")
        if reciprocity:
            fw("receiverPositions")
        else:
            fw("[0 0 0; 0 0 0]")
        fw(",frequencyDependency,nearFieldCalculation,reference);")
        file.close

# ----------------------- Render pictures of the model -------------------------
        if pictures:
            for ii in range(0, len(renderloc)):
                cam.location = (renderloc[ii][0]*camradius, renderloc[ii][1]*camradius, renderloc[ii][2]*camradius)
                cam.rotation_euler = (renderrot[ii][0], renderrot[ii][1], renderrot[ii][2])
                lamp.location = (renderloc[ii][0]*lampradius, renderloc[ii][1]*lampradius, renderloc[ii][2]*lampradius)
                bpy.ops.render.render()
                temp = ("%s/Pictures/" % filepath1)
                if not os.path.exists(temp):
                    os.mkdir(temp)
                temp = ("%d-%d" % (rendernam[ii][0], rendernam[ii][1]))
                bpy.data.images['Render Result'].save_render("%s/Pictures/%s.png" % (filepath1, temp))

# ----------------------- Write NumCalc input files for all CPUs and Cores -----
        for core in range(1, 9):
            for cpu in range(1, 11):
                if not cpusAndCores[cpu-1][core-1] == 0:

                    filepath2 = ("%s/NumCalc/CPU_%i_Core_%i/" % (filepath1, cpu, core))
                    if not os.path.exists(filepath2):
                        os.mkdir(filepath2)

                    file = open(("%s%s" % (filepath2, filename1)), "w", encoding="utf8", newline="\n")
                    fw = file.write
                    
                    if frequencyDependency:
                        if cpusAndCores[cpu-1][core-1] == 1:
                            tmpEar = "L"
                        if cpusAndCores[cpu-1][core-1] == 2:
                            tmpEar = "R"
                        tmpfmax = max(frequencies[cpu-1][core-1])
                        tmpfdiff = 24000
                        for ff in maxObjectFrequency:
                            if (ff-tmpfmax) >= 0 and (ff-tmpfmax) < tmpfdiff:
                                obj_name = ("%s%i" % (tmpEar, ff))
                                tmpfdiff = ff-tmpfmax
                    else:
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
                    if nearFieldCalculation:
                        numNodes = numNodes+len(NFGrid_data.vertices[:])
                    else:
                        if not evaluationGrid1 == 'None':
                            nodes = open("%s/EvaluationGrids/%s/Nodes.txt" % (filepath1, evaluationGrid1))
                            line = nodes.readline()
                            numNodes = numNodes+int(line)
                        if not evaluationGrid2 == 'None':
                            nodes = open("%s/EvaluationGrids/%s/Nodes.txt" % (filepath1, evaluationGrid2))
                            line = nodes.readline()
                            numNodes = numNodes+int(line)
                        if not evaluationGrid3 == 'None':
                            nodes = open("%s/EvaluationGrids/%s/Nodes.txt" % (filepath1, evaluationGrid3))
                            line = nodes.readline()
                            numNodes = numNodes+int(line)
                        if not evaluationGrid4 == 'None':
                            nodes = open("%s/EvaluationGrids/%s/Nodes.txt" % (filepath1, evaluationGrid4))
                            line = nodes.readline()
                            numNodes = numNodes+int(line)
                        if not evaluationGrid5 == 'None':
                            nodes = open("%s/EvaluationGrids/%s/Nodes.txt" % (filepath1, evaluationGrid5))
                            line = nodes.readline()
                            numNodes = numNodes+int(line)
                    numElements = 0
                    if nearFieldCalculation:
                        numElements = numElements+len(NFGrid_data.polygons[:])
                    else:
                        if not evaluationGrid1 == 'None':
                            elements = open("%s/EvaluationGrids/%s/Elements.txt" % (filepath1, evaluationGrid1))
                            line = elements.readline()
                            numElements = numElements+int(line)
                        if not evaluationGrid2 == 'None':
                            elements = open("%s/EvaluationGrids/%s/Elements.txt" % (filepath1, evaluationGrid2))
                            line = elements.readline()
                            numElements = numElements+int(line)
                        if not evaluationGrid3 == 'None':
                            elements = open("%s/EvaluationGrids/%s/Elements.txt" % (filepath1, evaluationGrid3))
                            line = elements.readline()
                            numElements = numElements+int(line)
                        if not evaluationGrid4 == 'None':
                            elements = open("%s/EvaluationGrids/%s/Elements.txt" % (filepath1, evaluationGrid4))
                            line = elements.readline()
                            numElements = numElements+int(line)
                        if not evaluationGrid5 == 'None':
                            elements = open("%s/EvaluationGrids/%s/Elements.txt" % (filepath1, evaluationGrid5))
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
                    if not frequencyDependency:
                        fw("../../ObjectMeshes/Reference/Nodes.txt\n")
                    else:
                        fw("../../ObjectMeshes/%s/Nodes.txt\n" % obj_name)
                    if nearFieldCalculation:
                        fw("../../EvaluationGrids/NF_Sphere/Nodes.txt\n")
                    else:
                        if not evaluationGrid1 == 'None':
                            fw("../../EvaluationGrids/%s/Nodes.txt\n" % evaluationGrid1)
                        if not evaluationGrid2 == 'None':
                            fw("../../EvaluationGrids/%s/Nodes.txt\n" % evaluationGrid2)
                        if not evaluationGrid3 == 'None':
                            fw("../../EvaluationGrids/%s/Nodes.txt\n" % evaluationGrid3)
                        if not evaluationGrid4 == 'None':
                            fw("../../EvaluationGrids/%s/Nodes.txt\n" % evaluationGrid4)
                        if not evaluationGrid5 == 'None':
                            fw("../../EvaluationGrids/%s/Nodes.txt\n" % evaluationGrid5)
                    fw("##\n")
                    fw("ELEMENTS\n")
                    if not frequencyDependency:
                        fw("../../ObjectMeshes/Reference/Elements.txt\n")
                    else:
                        fw("../../ObjectMeshes/%s/Elements.txt\n" % obj_name)
                    if nearFieldCalculation:
                        fw("../../EvaluationGrids/NF_Sphere/Elements.txt\n")
                    else:
                        if not evaluationGrid1 == 'None':
                            fw("../../EvaluationGrids/%s/Elements.txt\n" % evaluationGrid1)
                        if not evaluationGrid2 == 'None':
                            fw("../../EvaluationGrids/%s/Elements.txt\n" % evaluationGrid2)
                        if not evaluationGrid3 == 'None':
                            fw("../../EvaluationGrids/%s/Elements.txt\n" % evaluationGrid3)
                        if not evaluationGrid4 == 'None':
                            fw("../../EvaluationGrids/%s/Elements.txt\n" % evaluationGrid4)
                        if not evaluationGrid5 == 'None':
                            fw("../../EvaluationGrids/%s/Elements.txt\n" % evaluationGrid5)
                    fw("##\n")
                    fw("# SYMMETRY\n")
                    fw("# 0 0 0\n")
                    fw("# 0.0000e+00 0.0000e+00 0.0000e+00\n")
                    fw("##\n")
                    if reciprocity:
                        '''
                        if cpusAndCores[cpu-1][core-1]==1:
                            tmpEar='Left ear'
                        if cpusAndCores[cpu-1][core-1]==2:
                            tmpEar='Right ear'
                        '''
                        #changed tmpEar to ear 
                        if ear=='Right ear':
                            tmpEar='Right ear'
                        if ear=='Left ear':
                            tmpEar='Left ear'
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
            bpy.data.objects[obj.name].select = False
        if nearFieldCalculation:
            bpy.data.objects['NFGrid'].select = True
            bpy.ops.object.delete()

        return {'FINISHED'}


# ----------------------- Blender add-on registration --------------------------
def menu_func_export(self, context):
    self.layout.operator(ExportMesh2HRTF.bl_idname, text="Mesh2HRTF")


def register():
    bpy.utils.register_class(ExportMesh2HRTF)
    bpy.types.INFO_MT_file_export.append(menu_func_export)


def unregister():
    bpy.utils.unregister_class(ExportMesh2HRTF)
    bpy.types.INFO_MT_file_export.remove(menu_func_export)
