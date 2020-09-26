# Semi-automatic head centering
#
# Author: The Mesh2HRTF developers
#
# See class MaterialAssignment(bpy.types.Operator) for information on how to
# use this file. For debugging use the commeted code at the end insude Blenders
# text editor.
#
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


import bpy
import math

class CenterHead(bpy.types.Operator):
    """
    Center the head using three user selected Points to align the interaural
    center with the center of coordinates, align the interaural axis with the
    y-axis, and make the head look in positive x-direction.

    NOTE: The head mesh is not rotated above the y-axis (up/down). It is
    assumed that the natuarl hearing position is already given. Manually
    rotate if this is not the case.

    Usage
    ------

    1. Load the script into blender, e.g., by copying it to scripts/startup
       (https://docs.blender.org/api/current/info_overview.html).
    2. Pre-rotate the head until it is in generally in the right position
       (looking in x-direction and z-axis points upwards)
    3. Mark the center of the left ear channel and add a lamp point.
       (This will automatically be called "Point")
    4. Mark the center of the right ear channel and add a lamp point.
       (This will automatically be calles "Point.001")
    5. Run head centering by typing `bpy.ops.object.centerhead()`
       in Blenders python terminal or use the debuggin code at the end
       of the script

    A good way to set a lamp point is
    ---------------------------------

    1. Select the vertex where the point should be placed in edit mode
    2. Snap the 3D cursor to the vertex (Shift+s, 3)
    3. Add the lamp point in object mode

    """
    bl_idname = "object.centerhead"
    bl_label = "CenteringHead"

    def execute(self, context):
        center_head()
        return {'FINISHED'}


def center_head():
	# Get location of markers
	left_ear_canal = bpy.data.objects['Point']
	left = left_ear_canal.location
	print("Left ear channel is at ", left)

	right_ear_canal = bpy.data.objects['Point.001']
	right = right_ear_canal.location
	print("Right ear channel is at ", right)


	# Distance between ear channles
	left_right_distance=math.sqrt((right[0] - left[0])**2
	                            + (right[1] - left[1])**2
	                            + (right[2] - left[2])**2)
	print("distance between ears: ", left_right_distance)

	# Rotate around x-axis:
	# calculating rotation in degree
	delta_z=left[2]-right[2]
	print("height difference between ear channels: ", delta_z)
	alpha=math.asin(delta_z/left_right_distance)

	# rotate around the left ear
	bpy.context.scene.cursor.location = left
	bpy.ops.object.origin_set(type='ORIGIN_CURSOR')
	bpy.ops.transform.rotate(value=alpha, orient_axis='X')


	# rotate around z-axis
	# calculating rotation degree
	delta_x=right[0]-left[0]
	print("Front/back difference between ear channels: ", delta_x)
	beta=math.asin(delta_x/left_right_distance)

	# rotate around the left ear
	bpy.context.scene.cursor.location = left
	bpy.ops.object.origin_set(type='ORIGIN_CURSOR')
	bpy.ops.transform.rotate(value=beta, orient_axis='Z')


	# correcting x/y/z-offset
	translate=(-left[0], -left[1]+left_right_distance/2, -left[2])
	print("head offset in x/y/z-direction: ", translate)
	bpy.ops.transform.translate(value=translate)


	# remove the added lamp points
	bpy.data.objects.remove(bpy.data.objects['Point'])
	bpy.data.objects.remove(bpy.data.objects['Point.001'])


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
