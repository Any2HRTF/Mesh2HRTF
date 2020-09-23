# material_and_assignment.py

# Created by Robert Pelzer, Audio Communication Group, Technical University Berlin on 27.06.2017
#################################################

# This script will create the materials Skin, Left and Right ear and will automatically assign the materials accordingly. 
# Two conditions have to be considered: 
#   1: the head is already alligned and centered so that the ear canals are on the y-axis. This is important because this script defines the ear elements by finding the elements closest to the y-axis.
#   2: if the head is remeshed the name of the ear (left or right) has to appear in the name. When neither of the expressions occur it is assumed that no re-mashing has taken place yet
#
# To start the script load the mesh in object mode, select it and run the script.

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
import bmesh
import math

head = bpy.context.active_object
name=head.name



#################################
#material part
################################

# Create Material Skin
skin = bpy.data.materials.new(name="Skin")
# Assign Colour
skin.diffuse_color=(0,1,0)
#Create new amterial slot
head.data.materials.append(skin)

# override false names, such as Skin.001 etc
bpy.context.object.active_material_index = 0
bpy.context.object.active_material.name = "Skin"

# Create Material for left ear
left_ear= bpy.data.materials.new(name="Left ear")
# Assign Colour
left_ear.diffuse_color=(0,0,1)
#Create new amterial slot
head.data.materials.append(left_ear)
# override false names
bpy.context.object.active_material_index = 1
bpy.context.object.active_material.name = "Left ear"

# Create Material for right ear
right_ear= bpy.data.materials.new(name="Right ear")
# Assign Colour
right_ear.diffuse_color=(1,0,0)
#Create new amterial slot
head.data.materials.append(right_ear)
# override false names
bpy.context.object.active_material_index = 2
bpy.context.object.active_material.name = "Right ear"

###switch into edit mode
bpy.ops.object.editmode_toggle()


##############################################
#finding ears and assigning material
##############################################




obj = bpy.context.edit_object
me = obj.data
bm = bmesh.from_edit_mesh(me)
bpy.ops.mesh.select_all(action = 'DESELECT')

# use before accessing bm.verts[] with blender 2.73
if hasattr(bm.verts, "ensure_lookup_table"): 
    bm.verts.ensure_lookup_table()
    bm.faces.ensure_lookup_table()

# initialize helper variables
left_index=0
right_index=0
dist_l_old=1
dist_r_old=1
y_l_old=1
y_r_old=1
y_delta=0.002

# notice in Bmesh polygons are called faces

# iterate through all faces
for face in bm.faces:
    #if face.select:
        #get the location
        face_location = face.calc_center_median()
        

        # because of deformed and graded meshes cases have to be handled differently for each side!
        if "left" in name:
    	
	        # finding left ear
	        if face_location[1]>0: 
	            
	            # when it´ getting close
	            if abs(face_location[0])<0.001 and abs(face_location[2])<0.001:
	                #y-location
	                y_loc_l=abs(face_location[1])
	                #x-z distance sum 
	                dist_l=abs(face_location[0]) + abs(face_location[2])
	                
	                #when it´ getting closer - y-location must be monitored for cases where the y-axis passes through the tragus
	                if dist_l < dist_l_old and not (y_loc_l > (y_l_old+ y_delta)) or y_loc_l <  (y_l_old - y_delta):
	                    left_index=face.index
	                    dist_l_old=dist_l
	                    y_l_old= y_loc_l
	   
	                
	        # finding graded right ear
	        if face_location[1]<0: 
	            
	            # when it´ getting close
	            if abs(face_location[0])<0.01 and abs(face_location[2])<0.01:
	                #y-location
	                y_loc_r=abs(face_location[1])
	                #x-z distance sum 
	                dist_r=abs(face_location[0]) + abs(face_location[2])
	                
	                if dist_r < dist_r_old:
	            
	                    right_index=face.index
	                    dist_r_old=dist_r
	                    y_r_old= y_loc_r

        elif "right" in name:
    	
	        # finding graded left  ear 
	        if face_location[1]>0: 
	            
	            # when it´ getting close
	            if abs(face_location[0])<0.01 and abs(face_location[2])<0.01:
	                #y-location
	                y_loc_l=abs(face_location[1])
	                #x-z distance sum 
	                dist_l=abs(face_location[0]) + abs(face_location[2])
	                
	                
	                if dist_l < dist_l_old:
	                    left_index=face.index
	                    dist_l_old=dist_l
	                    y_l_old= y_loc_l
	   
	                
	        # finding right ear
	        if face_location[1]<0: 
	            
	            # when it´ getting close
	            if abs(face_location[0])<0.001 and abs(face_location[2])<0.001:
	                #y-location
	                y_loc_r=abs(face_location[1])
	                #x-z distance sum 
	                dist_r=abs(face_location[0]) + abs(face_location[2])
	                
	                #when it´ getting closer - y-location must be monitored for cases where the y-axis passes through the tragus
	                if dist_r < dist_r_old and not(y_loc_r > (y_r_old+ y_delta)) or y_loc_r <  (y_r_old - y_delta):
	            
	                    right_index=face.index
	                    dist_r_old=dist_r
	                    y_r_old= y_loc_r


        else:
                
            # finding left ear
	        if face_location[1]>0: 
	            
	            # when it´ getting close
	            if abs(face_location[0])<0.001 and abs(face_location[2])<0.001:
	                #y-location
	                y_loc_l=abs(face_location[1])
	                #x-z distance sum 
	                dist_l=abs(face_location[0]) + abs(face_location[2])
	                
	                #when it´ getting closer - y-location must be monitored for cases where the y-axis passes through the tragus
	                if dist_l < dist_l_old and not (y_loc_l > (y_l_old+ y_delta)) or y_loc_l <  (y_l_old - y_delta):
	                    left_index=face.index
	                    dist_l_old=dist_l
	                    y_l_old= y_loc_l    
    	

	   
	                
	        # finding right ear
	        if face_location[1]<0: 
	            
	            # when it´ getting close
	            if abs(face_location[0])<0.001 and abs(face_location[2])<0.001:
	                #y-location
	                y_loc_r=abs(face_location[1])
	                #x-z distance sum 
	                dist_r=abs(face_location[0]) + abs(face_location[2])
	                
	                #when it´ getting closer - y-location must be monitored for cases where the y-axis passes through the tragus
	                if dist_r < dist_r_old and not(y_loc_r > (y_r_old+ y_delta)) or y_loc_r <  (y_r_old - y_delta):
	            
	                    right_index=face.index
	                    dist_r_old=dist_r
	                    y_r_old= y_loc_r

#print(f.index)

print(face_location, left_index, right_index)

bm.faces[left_index].select = True
bpy.context.object.active_material_index = 1
bpy.ops.object.material_slot_assign()
bm.faces[left_index].select = False
  
bm.faces[right_index].select = True 
bpy.context.object.active_material_index = 2
bpy.ops.object.material_slot_assign()
bm.faces[right_index].select = False

# Show the updates in the viewport
# and recalculate n-gon tessellation.
#bmesh.update_edit_mesh(me, True)

bpy.ops.object.editmode_toggle()
#renaming the object
head.name="Reference"