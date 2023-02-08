import mesh2scattering as m2s
import os
import trimesh
import filecmp
import pytest
import numpy as np
import numpy.testing as npt
from tempfile import TemporaryDirectory
import pyfar as pf


def test_import():
    from mesh2scattering import input
    assert input


def test_write_mesh(tmpdir):
    path = os.path.join(
        m2s.utils.repository_root(), '..', 'tests', 'references', 'Mesh')
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


@pytest.mark.parametrize("n_dim", [3, 2])
@pytest.mark.parametrize("coordinates,show", [[False, True], [True, False]])
def test_read_and_write_evaluation_grid(n_dim, coordinates, show):
    cwd = os.path.dirname(__file__)
    data_grids = os.path.join(cwd, 'resources', 'evaluation_grids')

    tmp = TemporaryDirectory()

    # sampling grids
    if n_dim == 3:
        # 3D sampling grid (Lebedev, first order)
        points = np.array([
            [1., 0., 0.],
            [-1., 0., 0.],
            [0, 1., 0.],
            [0, -1., 0.],
            [0, 0., 1.],
            [0, 0., -1.]])
        discard = None
    else:
        # 2D sampling grid (all z = 0)
        points = np.array([
            [1., 0., 0.],
            [-1., 0., 0.],
            [0, 1., 0.],
            [0, -1., 0.]])
        discard = "z"

    # pass as Coordinates object
    if coordinates:
        points = pf.Coordinates(points[:, 0], points[:, 1], points[:, 2])

    # write grid
    m2s.input.write_evaluation_grid(
        points, os.path.join(tmp.name, "test"), discard=discard, show=show)

    # check if the plot exists
    if show:
        assert os.path.isfile(
            os.path.join(tmp.name, "test", "evaluation_grid.png"))
    else:
        assert not os.path.isfile(
            os.path.join(tmp.name, "test", "evaluation_grid.png"))

    # check the nodes and elements
    for file in ["Nodes.txt", "Elements.txt"]:
        with open(os.path.join(data_grids, f"{n_dim}D", file), "r") as f:
            ref = "".join(f.readlines())
        with open(os.path.join(tmp.name, "test", file), "r") as f:
            test = "".join(f.readlines())

        assert test == ref

    # read the grid
    coordinates = m2s.input.read_evaluation_grid(
        os.path.join(tmp.name, "test"))

    # check grid
    assert isinstance(coordinates, pf.Coordinates)
    npt.assert_equal(coordinates.get_cart(), points)


def test_write_material():
    # test write boundary condition with default values

    tmp = TemporaryDirectory()
    filename = os.path.join(tmp.name, "test_material.csv")

    # write data
    m2s.input.write_material(
        filename, "admittance", [100, 200], [1 + 0j, 1.5 + 0.5j])

    # read and check data
    with open(filename, "r") as f_id:
        file = f_id.readlines()
    file = "".join(file)

    assert file.startswith("# Keyword to define the boundary condition:\n")
    assert file.endswith("100, 1.0, 0.0\n200, 1.5, 0.5\n")


@pytest.mark.parametrize("kind,check_kind", (
    ["admittance", ["ADMI", "PRES", "VELO"]],
    ["pressure", ["PRES", "ADMI", "VELO"]],
    ["velocity", ["VELO", "ADMI", "PRES"]]))
def test_write_material_kind(kind, check_kind):
    # test if the kind of boundary condition is written correctly

    tmp = TemporaryDirectory()
    filename = os.path.join(tmp.name, "test_material.csv")

    # write data
    m2s.input.write_material(
        filename, kind, [100, 200], [1 + 0j, 1.5 + 0.5j])

    # read and check data
    with open(filename, "r") as f_id:
        file = f_id.readlines()

    assert f"{check_kind[0]}\n" in file
    assert f"{check_kind[1]}\n" not in file
    assert f"{check_kind[2]}\n" not in file


def test_write_material_comment():
    # test if the comment is written

    tmp = TemporaryDirectory()
    filename = os.path.join(tmp.name, "test_material.csv")
    comment = "Weird, random data"

    # write data
    m2s.input.write_material(
        filename, "pressure", [100, 200], [1 + 0j, 1.5 + 0.5j], comment)

    # read and check data
    with open(filename, "r") as f_id:
        file = f_id.readlines()

    assert file[0] == "# " + comment + "\n"
    assert file[1] == "#\n"
    assert file[2] == "# Keyword to define the boundary condition:\n"
