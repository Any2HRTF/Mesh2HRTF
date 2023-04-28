"""
See documentaion of the function CenterHead for more information.

Robert Pelzer, Fabian Brinkmann and Tim Wennemann,
Audio Communication Group, Technical University of Berlin, Germany
"""

import bpy
from bpy.props import BoolProperty, IntProperty
import bmesh
import numpy as np


class CenterHead(bpy.types.Operator):
    """
    Center a head or head and torso mesh as required for numerical HRTF
    simulations using Mesh2HRTF.

    Parameters
    ----------
    verbose : bool
        Output information to console if True. The default is False.
    precision : int
        Precision in decimals that is used to check the position of the
        selected vertices after the alignment. A precision of 0.1 mm is
        recommended, wich refers to ``precision=4`` if the mesh is in meters
        and ``precision=1`` if the mesh is in mm. The default is 4.

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

    verbose: BoolProperty(
        name = "verbose",
        description = "Output debugging information. False by default.",
        default = False,
        )
    precision: IntProperty(
        name = "precision",
        description = ("Precison for checking if the alignment worked."
                       "Given by the number of decimals. The default is 4, which"
                       "equals 0.1 mmif the mesh unit is m"),
        default = 4,
        )

    def execute(self, context):
        keywords = self.as_keywords(ignore=("check_existing", "filter_glob"))
        center_head(**keywords)
        return {'FINISHED'}


def center_head(verbose, precision):

    # Rotate around x-axis to bring ears to same height
    left, right, center = get_vertices()
    bpy.context.scene.cursor.location = left
    bpy.ops.object.origin_set(type='ORIGIN_CURSOR')
    alpha = np.arctan2(left[2]-right[2], left[1]-right[1])
    bpy.ops.transform.rotate(value=alpha, orient_axis='X')

    # rotate around z-axis to bring ears to same azimuth angle
    left, right, center = get_vertices()
    bpy.context.scene.cursor.location = left
    bpy.ops.object.origin_set(type='ORIGIN_CURSOR')
    gamma = -np.arctan2(left[0]-right[0], left[1]-right[1])
    bpy.ops.transform.rotate(value=gamma, orient_axis='Z')

    # align interaural axis and center
    left, right, center = get_vertices()
    translate = (-(left[0] + right[0]) / 2,
                 -(left[1] + right[1]) / 2,
                 -(left[2] + right[2]) / 2)
    bpy.ops.transform.translate(value=translate)

    # rotate around y-axis
    if center is not None:
        left, right, center = get_vertices()
        beta = np.arctan2(center[2], center[0])
        bpy.context.scene.cursor.location = (0.0, 0.0, 0.0)
        bpy.ops.object.origin_set(type='ORIGIN_CURSOR')
        bpy.ops.transform.rotate(value=-beta, orient_axis='Y')
    else:
        beta = 0

    # output to console
    if verbose:
        print(("Moved in x/y/z-direction by\n"
               f"{translate[0]:.2f}, {translate[1]:.2f}, {translate[2]:.2f} units\n"
               "Rotated around x/y/z-axis by\n"
               f"{alpha/np.pi*180:.2f}, {beta/np.pi*180:.2f}, "
               f"{gamma/np.pi*180:.2f} degrees\n"))

    # check success
    left, right, center = get_vertices()

    atol = 10**(-precision)
    err = ""
    if np.abs(left[0]) > atol or np.abs(left[2]) > atol:
        err += ("Left ear channel x or z coordinate are not zero\n"
                f"position is {left[0]:.6f}, {left[1]:.6f}, {left[2]:.6f}\n")
    if np.abs(right[0]) > atol or np.abs(right[2]) > atol:
        err += ("Left ear channel x or z coordinate are not zero:\n"
                f"position is {right[0]:.6f}, {right[1]:.6f}, {right[2]:.6f}\n")
    if np.abs(left[1] + right[1]) > atol:
        err += ("y coordinate of left and right ear is not identical\n"
                f"positions are {left[1]:.6f}, {right[1]:.6f}\n")
    if center is not None and np.abs(center[2]) > atol:
        err += ("Natural head orientation not established.\n"
                f"z-coordinate of reference point is {center[2]:.6f}\n")
    if err:
        raise ValueError(err)


def get_vertices():
    """ get curent positions of selected vertices"""
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

    return left, right, center



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
