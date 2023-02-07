import mesh2scattering as m2s
import pytest
import os
import trimesh
import filecmp
import pathlib

def test_write_mesh(tmpdir):
    path = os.path.join(
        m2s.repository_root(), '..', 'tests', 'references', 'Mesh')
    mesh_path = os.path.join(path, 'sample.stl')
    mesh = trimesh.load(mesh_path)
    m2s.input.write_mesh(mesh.vertices, mesh.faces, tmpdir, start=0)
    assert filecmp.cmp(
        os.path.join(path, 'Elements.txt'),
        os.path.join(tmpdir, 'Elements.txt')
    )
    assert filecmp.cmp(
        os.path.join(path, 'Nodes.txt'),
        os.path.join(tmpdir, 'Nodes.txt')
    )
