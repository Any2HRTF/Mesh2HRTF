import trimesh
import numpy as np
import os
from scipy.spatial import Delaunay, ConvexHull
import mesh2scattering as m2s


def write_mesh(vertices, faces, name='Reference', start=200000, discard=None):
    if vertices.ndim != 2 or vertices.shape[0] < 3 \
            or vertices.shape[1] != 3:
        raise ValueError(
            "vertices must be a 2D array of shape (N, 3) with N > 2")

    # discard dimension
    if discard == "x":
        mask = (1, 2)
    elif discard == "y":
        mask = (0, 2)
    elif discard == "z":
        mask = (0, 1)
    else:
        mask = (0, 1, 2)

    # triangulate
    if discard is None:
        tri = ConvexHull(vertices[:, mask])
    else:
        tri = Delaunay(vertices[:, mask])

    # check output directory
    if not os.path.isdir(name):
        os.mkdir(name)

    # write nodes
    N = int(vertices.shape[0])
    start = int(start)

    nodes = f"{N}\n"
    for nn in range(N):
        nodes += (f"{int(start + nn)} "
                  f"{vertices[nn, 0]} "
                  f"{vertices[nn, 1]} "
                  f"{vertices[nn, 2]}\n")

    with open(os.path.join(name, "Nodes.txt"), "w") as f_id:
        f_id.write(nodes)

    # write elements
    N = int(faces.shape[0])
    elems = f"{N}\n"
    for nn in range(N):
        elems += (f"{int(start + nn)} "
                  f"{tri.simplices[nn, 0] + start} "
                  f"{tri.simplices[nn, 1] + start} "
                  f"{tri.simplices[nn, 2] + start} "
                  "2 0 1\n")

    with open(os.path.join(name, "Elements.txt"), "w") as f_id:
        f_id.write(elems)


def write_stl(mesh_path, project_path):
    mesh = trimesh.load(r'D:\sciebo\2021_DFG-Projekt\data\meshes\ita_50k\sample.stl')
    write_mesh(mesh.vertices, mesh.faces)

    
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
