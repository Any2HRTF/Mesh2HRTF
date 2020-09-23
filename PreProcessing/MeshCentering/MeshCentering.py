# head_centering.py
#
# Created by Robert Pelzer, Audio Communication Group, Technical University Berlin on 15.06.2017
#################################################
#
# This script will center the head mesh in blender after the user selected three points (the left and righ ear canal, and a point on the forehead, e.g. the tip of the nose or a point between the eyes)
# After the centering, the interaural center will coincide with the center of coordinates, the interaural axsi will coincide with the y-axis, and the mesh will look in positive x-direction.
# NOTE: The head mesh is not rotated above the y-axis (up/down) because it is assued that the natruarl hearing position is already given. If this is not the case, manually rotate abvoe the y-axis after running this script.
#
# The following steps must be fowllowed:
#
# 0. pre-rotate the head with 90 degree steps until it is in generally in the right position (looking in x-direction and z-axis points upwards)
#
# 1. mark the position of the left ear canal and add a lamp point. This point will automatically be called "Point"
# 2. mark the position of the right ear canal and add another lamp point which will be called "Point.001"
# 4. mark a point on the forehead with another lamp point which will be called "Point.002"
# 3. select the head mesh by righ clicking on it and run the script

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

# get location of markers
left_ear_canal = bpy.data.objects['Point']
left = left_ear_canal.location
print("left: ", left)
right_ear_canal = bpy.data.objects['Point.001']
right = right_ear_canal.location


nose_tip = bpy.data.objects['Point.002']
nose = nose_tip.location
print("nose: ", nose)

# preperation to delete Points later
objs = bpy.data.objects


#print(nose[2])
#cube = bpy.data.objects['Cube']
#cube.location[1]+=1
#print(cube.location)

## Distance between points on ear canals by calculating the euclidean distance
left_right_distance=math.sqrt((right[0]-left[0])**2 + (right[1]-left[1])**2 + (right[2]-left[2])**2)
print("distance",left_right_distance)

###############################################
# First Step: rotating tilted around on x-axis:
###############################################

# how much elevation in z-direction?
delta_z=left[2]-right[2]
print("deltaz ", delta_z)
#calculating rotation degree
alpha=math.asin(delta_z/left_right_distance)


# set curser on left ear canal and rotate for alpha degrees around x-axis in respect to left ear
bpy.context.scene.cursor_location = left
bpy.ops.object.origin_set(type='ORIGIN_CURSOR')
bpy.ops.transform.rotate(value=alpha, axis=(1, 0, 0))

#calculate the new position of the nose tip via matrix multiplication
###################

# place vector into origin first
zp=[0,0,0]
zp[1]=nose[1]-left[1]
zp[2]=nose[2]-left[2]

# rotate and add center point(left) to translate back parallel
new_nose=[0,0,0]
new_nose[0]=nose[0]
new_nose[1]=math.cos(alpha)*zp[1] - math.sin(alpha)*zp[2] + left[1]
new_nose[2]=math.sin(alpha)*zp[1] + math.cos(alpha)*zp[2] + left[2]

nose=[new_nose[0], new_nose[1], new_nose[2]]
nose_tip.location=nose


###############################################
# Second Step: rotating turned head on z-axis
###############################################

# delta x on x-axis
delta_x=right[0]-left[0]
print("deltax ", delta_x)
#calculating rotation degree
beta=math.asin(delta_x/left_right_distance)


# set curser on left ear canal and rotate for alpha degrees in respect to that position
bpy.context.scene.cursor_location = left
bpy.ops.object.origin_set(type='ORIGIN_CURSOR')
bpy.ops.transform.rotate(value=beta, axis=(0, 0, 1))

#again calculate the new position of the nose tip via matrix multiplication
###################

# place vector into origin first
#zp=[0,0,0]
zp[0]=nose[0]-left[0]
zp[1]=nose[1]-left[1]

# rotate and add center point(left) to translate back parallel
#new_nose=[0,0,0]
new_nose[2]=nose[2]
new_nose[0]=math.cos(beta)*zp[0] - math.sin(beta)*zp[1] + left[0]
new_nose[1]=math.sin(beta)*zp[0] + math.cos(beta)*zp[1] + left[1]

nose=[new_nose[0], new_nose[1], new_nose[2]]
nose_tip.location=nose

###############################################
# Third Step: correcting offset on x-axis
###############################################
translate_x=(-left[0], 0, 0)
print("transx",translate_x)
bpy.ops.transform.translate(value=translate_x)

nose_tip.location[0]=nose[0]-left[0] 

###############################################
# Fourth Step: centering on y-axis
###############################################

## this centers the head to the middle of the distance of the two ear canals
#translate_y=(0, -(left_right_distance/2 + left[1]), 0)


## this centers the head to the nose tip !!
translate_y=(0, -nose[1], 0)
print("transy",translate_y)
bpy.ops.transform.translate(value=translate_y)

nose_tip.location[1]=0.0
###############################################
# Fifth Step: correcting offset on z-axis
###############################################

translate_z=(0, 0, -left[2])
print("transz",translate_z)
bpy.ops.transform.translate(value=translate_z)
nose_tip.location[2]=nose[2]-left[2]

##########
# remove the added lamp points
objs.remove(objs["Point"], True)
objs.remove(objs["Point.001"], True)
objs.remove(objs["Point.002"], True)

