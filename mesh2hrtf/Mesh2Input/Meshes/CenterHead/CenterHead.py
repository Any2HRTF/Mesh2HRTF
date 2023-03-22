"""
See documentaion of the function CenterHead for more information.

Robert Pelzer, Fabian Brinkmann and Tim Wennemann,
Audio Communication Group, Technical University of Berlin, Germany
"""

import bpy
import bmesh
import numpy as np


class CenterHead(bpy.types.Operator):
    """
    Center the head as required for numerical HRTF simulations using Mesh2HRTF.

    Usage
    ------

    1. Load the script into blender, e.g., by copying it to scripts/startup
       (https://docs.blender.org/api/current/info_overview.html) or run this
       file inside Blender.
    2. Roughly bring the head into the correct position. The center of the head
       should be close to the origin of coordinates. The nose should be
       pointing in positive x-direction, and the top of the head in positive
       z-direction.
    3. In edit mode, select two vertices to mark the center of the left and
       right ear channel entrances. These vertices define the interaural axis
       and their midpoint the interaural center. Based on this, the head is
       translated to make sure that the interaural center is at the origin of
       coordinates and rotated to make sure that the interaural axis coincides
       with the y-axis.
    4. An optional third point can be selected to establish a natural hearing
       position. If this is the case the mesh is rotated around the interaural
       axis (y-axis) to make sure the third selected vertex has a z-coordinate
       of 0, i.e., the vertex is on the x-y plane. We found a points on the
       nose to work good for establishing the natural hearing position.
    5. Run head centering by typing `bpy.ops.object.centerhead()`
       in Blenders python terminal or use the debugging code at the end
       of the script.

    """
    bl_idname = "object.centerhead"
    bl_label = "CenteringHead"

    def execute(self, context):
        center_head()
        return {'FINISHED'}


def center_head():
    # Set Blender into Edit Mode
    bpy.ops.object.mode_set(mode='EDIT')
    # Get mesh data
    obj = bpy.context.active_object
    bm = bmesh.from_edit_mesh(obj.data)
    # Get selected verts
    if any(v.select for v in bm.verts):
        verts = [v for v in bm.verts if v.select]
    else:
        raise RuntimeError("No vertex is selected.")
    if len(verts) in (2, 3):
        # Get index of left and right ear in verts
        point_idx = np.argsort([v.co.y for v in verts])
    else:
        raise RuntimeError(f"In total {len(verts)} points are selected. "
                           "Please select 2 or 3.")
    left = obj.matrix_world @ verts[point_idx[-1]].co
    right = obj.matrix_world @ verts[point_idx[0]].co
    if len(verts) == 3:
        center = obj.matrix_world @ verts[point_idx[1]].co
    else:
        center = None

    bpy.ops.object.mode_set(mode='OBJECT')

    # Distance between ear channles
    left_right_distance = np.sqrt((right[0] - left[0])**2
                                    + (right[1] - left[1])**2
                                    + (right[2] - left[2])**2)
    #print("distance between ears: ", left_right_distance)

    # Rotate around x-axis:
    # calculating rotation in degree
    delta_z = left[2] - right[2]
    #print("height difference between ear channels: ", delta_z)
    alpha = np.arcsin(delta_z/left_right_distance)

    # rotate around the left ear
    bpy.context.scene.cursor.location = left
    bpy.ops.object.origin_set(type='ORIGIN_CURSOR')
    bpy.ops.transform.rotate(value=alpha, orient_axis='X')

    # rotate around z-axis
    # calculating rotation degree
    delta_x = right[0]-left[0]
    #print("Front/back difference between ear channels: ", delta_x)
    gamma = np.arcsin(delta_x/left_right_distance)

    # rotate around the left ear
    bpy.context.scene.cursor.location = left
    bpy.ops.object.origin_set(type='ORIGIN_CURSOR')
    bpy.ops.transform.rotate(value=gamma, orient_axis='Z')

    # correcting x/y/z-offset
    translate = (-left[0], -left[1]+left_right_distance/2, -left[2])
    #print("head offset in x/y/z-direction: ", translate)
    bpy.ops.transform.translate(value=translate)

    # rotate around y-axis
    if center is not None:
        beta = np.arctan((center[2]+translate[2])/(center[0]+translate[0]))
        bpy.context.scene.cursor.location = (0.0, 0.0, 0.0)
        bpy.ops.object.origin_set(type='ORIGIN_CURSOR')
        bpy.ops.transform.rotate(value=-beta, orient_axis='Y')
        #print("head rotation around y-axis in deg:", np.rad2deg(-beta))


def register():
    bpy.utils.register_class(CenterHead)


def unregister():
    bpy.utils.unregister_class(CenterHead)


if __name__ == "__main__":
    register()


# for debugging only
# import CenterHead
# import importlib
# import bpy

# importlib.reload(CenterHead)
# CenterHead.center_head()
