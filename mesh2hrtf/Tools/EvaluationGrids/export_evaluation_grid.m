# Blender Plugin for exporting evaluation grids - undocumented.
import os
import bpy
from bpy.props import StringProperty, EnumProperty, IntProperty
from bpy_extras.io_utils import ExportHelper

bl_info = {
    "name": "Export Evaluation Grid",
    "author": "The Mesh2HRTF developers",
    "version": (0, 2, 0),
    "blender": (2, 80, 0),
    "location": "File > Export",
    "description": "Export evaluation grid",
    "warning": "",
    "wiki_url": "",
    "tracker_url": "",
    "support": "COMMUNITY",
    "category": "Import-Export"}


class ExportEvaluationgrid(bpy.types.Operator, ExportHelper):
    '''Export a single object as an evaluation grid'''

    bl_idname = "export_evaluationgrid.inp"
    bl_label = "Export Evaluation Grid"

    filename_ext = ".txt"
    filter_glob: StringProperty(default="*.txt", options={'HIDDEN'})

    offset: IntProperty(
            name="Offset",
            description="Node and element index offset",
            default=200000,
            min=0,
            max=10000000,
            )
    suffix: StringProperty(
            name="Element suffix",
            description="Element suffix",
            default=" 2 0 1",
            )
    unit: EnumProperty(
            name="Unit",
            description="Unit of the evaluation grid",
            items=[('m', 'm', 'Meter'), ('mm', 'mm', 'Millimeter')],
            default='mm',
            )

    @classmethod
    def poll(cls, context):
        return context.active_object is not None

    def execute(self, context):
        filepath = self.filepath
        filepath = bpy.path.ensure_ext(filepath, self.filename_ext)
        keywords = self.as_keywords(ignore=("check_existing", "filter_glob"))
        return ExportEvaluationgrid.save(self, context, **keywords)

    def draw(self, context):
        layout = self.layout
        layout.use_property_split = True
        layout.use_property_decorate = False

        row = layout.row()
        row.prop(self, "offset")
        row = layout.row()
        row.prop(self, "suffix")
        row = layout.row()
        row.prop(self, "unit")

    def save(operator,
             context,
             filepath="",
             offset=0,
             suffix=" 2 0 1",
             unit='mm',
             ):

        # ----------------------- Initialize constants ------------------------
        obj = context.active_object

        if not obj:
            raise Exception("Error, Select 1 active object")

        obj_data = obj.data

        unitFactor = 1
        if unit == 'mm':
            unitFactor = 0.001

        (filepath, filename) = os.path.split(filepath)

        # ----------------------- Write object data ---------------------------
        obj_data = obj.data
        file = open(os.path.join(filepath, "Nodes.txt"), "w",
                    encoding="utf8", newline="\n")
        fw = file.write
        fw("%i\n" % len(obj_data.vertices[:]))
        for ii in range(len(obj_data.vertices[:])):
            fw("%i " % (ii+offset))
            fw("%.6f %.6f %.6f\n" % (obj_data.vertices[ii].co[0] * unitFactor,
                                     obj_data.vertices[ii].co[1] * unitFactor,
                                     obj_data.vertices[ii].co[2] * unitFactor))
        file.close

        file = open(os.path.join(filepath, "Elements.txt"), "w",
                    encoding="utf8", newline="\n")
        fw = file.write
        fw("%i\n" % len(obj_data.polygons[:]))
        if len(obj_data.polygons[0].vertices[:]) == 3:
            for ii in range(len(obj_data.polygons[:])):
                fw("%i " % (ii+offset))
                fw("%d %d %d" % (obj_data.polygons[ii].vertices[0] + offset,
                                 obj_data.polygons[ii].vertices[1] + offset,
                                 obj_data.polygons[ii].vertices[2] + offset))
                fw("%s\n" % suffix)
        else:
            for ii in range(len(obj_data.polygons[:])):
                fw("%i " % (ii+offset))
                fw("%d %d %d %d" % (obj_data.polygons[ii].vertices[0] + offset,
                                    obj_data.polygons[ii].vertices[1] + offset,
                                    obj_data.polygons[ii].vertices[2] + offset,
                                    obj_data.polygons[ii].vertices[3] + offset)
                   )
                fw("%s\n" % suffix)
        file.close

        return {'FINISHED'}


def menu_func_export(self, context):
    self.layout.operator(ExportEvaluationgrid.bl_idname,
                         text="Evaluation grid (.txt)")


def register():
    bpy.utils.register_class(ExportEvaluationgrid)
    bpy.types.TOPBAR_MT_file_export.append(menu_func_export)


def unregister():
    bpy.utils.unregister_class(ExportEvaluationgrid)
    bpy.types.TOPBAR_MT_file_export.append(menu_func_export)
