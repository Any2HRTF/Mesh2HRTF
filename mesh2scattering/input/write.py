import trimesh
import numpy as np
import os


def write_mesh(vertices, faces, path, start=200000, discard=None):
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

    # check output directory
    if not os.path.isdir(path):
        os.mkdir(path)

    # write nodes
    N = int(vertices.shape[0])
    start = int(start)

    nodes = f"{N}\n"
    for nn in range(N):
        nodes += (f"{int(start + nn)} "
                  f"{vertices[nn, 0]} "
                  f"{vertices[nn, 1]} "
                  f"{vertices[nn, 2]}\n")

    with open(os.path.join(path, "Nodes.txt"), "w") as f_id:
        f_id.write(nodes)

    # write elements
    N = int(faces.shape[0])
    elems = f"{N}\n"
    for nn in range(N):
        elems += (f"{int(start + nn)} "
                  f"{faces[nn, 0] + start} "
                  f"{faces[nn, 1] + start} "
                  f"{faces[nn, 2] + start} "
                  "0 0 0\n")

    with open(os.path.join(path, "Elements.txt"), "w") as f_id:
        f_id.write(elems)


def write_stl(mesh_path, project_path):
    mesh = trimesh.load(mesh_path)
    path = os.path.join(project_path, 'ObjectMeshes', 'Reference')
    write_mesh(mesh.vertices, mesh.faces, path, start=0)
