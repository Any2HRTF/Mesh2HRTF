"""
See documentaion of the function AssignMaterials for more information.

Robert Pelzer and Fabian Brinkmann, Audio Communication Group, Technical
University of Berlin, Germany
"""

import bpy
import bmesh
import math
from bpy.props import FloatProperty, StringProperty


class AssignMaterials(bpy.types.Operator):
    '''Create and assign skin and ear materials for Mesh2HRTF export.

    Actions on currently active object
    ----------------------------------
    0. remove all existing materials from the object
    1. assign the material `Skin` to the entire object
    2. assign the materials `Left ear` to the face that is closest to the
       y-axis and has a positive y-value
    3. assign the materials `Right ear` to the face that is closest to the
       y-axis and has a negative y-value
    4. rename the object to `Reference`

    Requirements
    ------------
    the head is already alligned and centered, i.e., the ear canals are on the
    the y-axis. This is important because this script defines the ear elements
    by finding the faces closest to the y-axis.

    Usage
    -----
    1. Load the script into blender, e.g., by copying it to scripts/startup
       (https://docs.blender.org/api/current/info_overview.html).
    2. execute `bpy.ops.object.assignmaterials()` in Blenders python
       terminal or inside a script. Can also called with custom tolerance,
       e.g., `bpy.ops.object.assignmaterials(tolerance=3)`

    Input
    -----
    tolerance : FloatProperty
        maximum distance from y-axis in mesh units ('mm' or 'm'). The default
        is 2.
    ear : StringProperty
        The ear canals to assign the materials to. Can be ``'Both ears'``,
        ``'Left ear'`` or ``'Right ear'``. The default is 'Both ears'.

    '''
    bl_idname = "object.assignmaterials"
    bl_label = "AssignMaterials"

    tolerance: FloatProperty(
        name = "tolerance",
        description = "maximum distance from y-axis in mesh units \
                       ('mm' or 'm'). The default is 2.",
        default = 2,
        )
    ear: StringProperty(
        name = "ear",
        description = "The ear canals to assign the materials to. \
                       The default is 'Both ears'.",
        default = "Both ears",
        )

    def execute(self, context):
        # unit = 'mm'
        keywords = self.as_keywords(ignore=("check_existing", "filter_glob"))
        assign_material(bpy.context.active_object, **keywords)
        return {'FINISHED'}


def setup_materials(obj):
    '''Remove existing and create new materials required by Mesh2HRTF.'''

    # remove existing materials
    for i in range(len(bpy.context.object.material_slots)):
        bpy.context.object.active_material_index = 0
        bpy.ops.object.material_slot_remove({'object': obj})

    # assign materials required for Mesh2HRTF
    for idx, (name, color) in enumerate(zip(
            ["Skin", "Left ear", "Right ear"],
            [(0, 1, 0, 1), (0, 0, 1, 1), (1, 0, 0, 1)])):

        # Create the material
        material = bpy.data.materials.new(name=name)
        # Assign Colour
        material.diffuse_color = color
        material.specular_intensity = .1
        # Create new amterial slot
        obj.data.materials.append(material)
        # override false names, such as Skin.001 etc
        bpy.context.object.active_material_index = idx
        bpy.context.object.active_material.name = name


def get_ear_indices(bm, obj, tolerance, ear):
    '''Return indicees of faces at the entrance to the ear channels'''

    bm.verts.ensure_lookup_table()
    bm.faces.ensure_lookup_table()

    # initialize helper variables
    left_indices = []
    left_xy_distance = []
    left_y = []
    right_indices = []
    right_xy_distance = []
    right_y = []
    min_y = [1000, 1000]

    # world matrix for obtaining vertex positions in world coorindates
    world = obj.matrix_world

    # find possible faces at the ear channel entrances
    for face in bm.faces:

        # current face location in world coordinates
        xyz = world @ face.calc_center_median()
        # distance from y axis, and y-value
        xy_distance = abs(xyz[0]) + abs(xyz[2])
        y = xyz[1]

        # potential left ear
        if y > 0 and abs(xyz[0]) < tolerance and abs(xyz[2]) < tolerance:
            left_indices.append(face.index)
            left_xy_distance.append(xy_distance)
            left_y.append(abs(y))
            if abs(y) < min_y[0]:
                min_y[0] = abs(y)
        # potential right ear
        elif y < 0 and abs(xyz[0]) < tolerance and abs(xyz[2]) < tolerance:
            right_indices.append(face.index)
            right_xy_distance.append(xy_distance)
            right_y.append(abs(y))
            if abs(y) < min_y[1]:
                min_y[1] = abs(y)

    # find left ear element and try to exclude faces at the tragus
    left_index = None
    if ear in ("Both ears", "Left ear"):
        min_xy_dist = 1000
        for n in range(len(left_indices)):
            if left_y[n] < min_y[0] + tolerance and left_xy_distance[n] < \
                    min_xy_dist:
                min_xy_dist = left_xy_distance[n]
                left_index = left_indices[n]
    # find left ear element and try to exclude faces at the tragus
    right_index = None
    if ear in ("Both ears", "Right ear"):
        min_xy_dist = 1000
        for n in range(len(right_indices)):
            if right_y[n] < min_y[1] + tolerance and right_xy_distance[n] < \
                    min_xy_dist:
                min_xy_dist = right_xy_distance[n]
                right_index = right_indices[n]

    return left_index, right_index


def assign_material(obj, tolerance, ear):

    # Switch to object mode to avoid export errors
    bpy.ops.object.mode_set(mode='OBJECT', toggle=False)

    # create materials
    setup_materials(obj)

    # get indicees
    print(f"Search tolerance: {tolerance} blender units")

    bpy.ops.object.mode_set(mode='EDIT')
    bm = bmesh.from_edit_mesh(obj.data)
    bm.faces.ensure_lookup_table()

    def point_in_triangle_xz_plane(p, a, b, c):
        v0 = (c[0] - a[0], c[1] - a[1])
        v1 = (b[0] - a[0], b[1] - a[1])
        v2 = (p[0] - a[0], p[1] - a[1])

        dot00 = v0[0]*v0[0] + v0[1]*v0[1]
        dot01 = v0[0]*v1[0] + v0[1]*v1[1]
        dot02 = v0[0]*v2[0] + v0[1]*v2[1]
        dot11 = v1[0]*v1[0] + v1[1]*v1[1]
        dot12 = v1[0]*v2[0] + v1[1]*v2[1]

        denom = dot00 * dot11 - dot01 * dot01
        if denom == 0:
            return False
        inv_denom = 1 / denom
        u = (dot11 * dot02 - dot01 * dot12) * inv_denom
        v = (dot00 * dot12 - dot01 * dot02) * inv_denom

        return (u >= 0) and (v >= 0) and (u + v <= 1)

    def find_face_intersect_y_axis(obj, bm, side):
        intersecting_faces = []
        origin = (0.0, 0.0)

        for face in bm.faces:
            verts = face.verts
            if len(verts) != 3:
                continue

            a = (verts[0].co.x, verts[0].co.z)
            b = (verts[1].co.x, verts[1].co.z)
            c = (verts[2].co.x, verts[2].co.z)

            if point_in_triangle_xz_plane(origin, a, b, c):
                intersecting_faces.append(face)

        best_y_face = None
        if side == 'left':
            best_avg_y = float('inf')
        elif side == 'right':
            best_avg_y = -float('inf')

        for face in intersecting_faces:
            avg_y = sum(v.co.y for v in face.verts) / len(face.verts)
            if side == 'left':
                if avg_y < best_avg_y and avg_y > 0:
                    best_avg_y = avg_y
                    best_y_face = face
            elif side == 'right':
                if avg_y > best_avg_y and avg_y < 0 :
                    best_avg_y = avg_y
                    best_y_face = face

        return best_y_face

    if ear in ("Both ears", "Left ear"):
        left_face = find_face_intersect_y_axis(obj, bm, side='left')
        left_y_coord = get_vertex_y(left_face)
        left_face.material_index = 1
        print(f"Left ear index: {left_face.index}, left y-coordinate: {left_y_coord}")

    if ear in ("Both ears", "Right ear"):
        right_face = find_face_intersect_y_axis(obj, bm, side='right')
        right_y_coord = get_vertex_y(right_face)
        right_face.material_index = 2
        print(f"Right ear index: {right_face.index}, rigjt y-coordinate: {right_y_coord}")

    # bm.to_mesh(obj.data)
    bm.free()
    bpy.ops.object.mode_set(mode='OBJECT')

    # renaming the object
    obj.name = "Reference"

def get_vertex_y(face):
    vertex=face.verts[0]
    return vertex.co.y



def register():
    bpy.utils.register_class(AssignMaterials)


def unregister():
    bpy.utils.unregister_class(AssignMaterials)


if __name__ == "__main__":
    register()


# for debugging only
# import AssignMaterials
# import importlib
# import bpy
#
# importlib.reload(AssignMaterials)
# AssignMaterials.assign_material(bpy.context.active_object, tolerance=2)
