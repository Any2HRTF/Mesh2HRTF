import mesh2scattering as m2s
import os
import trimesh
import filecmp
import pytest
import numpy as np
import numpy.testing as npt
from tempfile import TemporaryDirectory
import pyfar as pf
import json


def test_import():
    from mesh2scattering import input
    assert input


def test_write_mesh(tmpdir):
    path = os.path.join(
        m2s.utils.program_root(), '..', 'tests', 'references', 'Mesh')
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
        points, os.path.join(tmp.name, "test"), discard=discard)

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


def test_create_source_positions():
    source_azimuth_deg = np.arange(0, 95, 10)
    source_colatitude_deg = np.arange(10, 85, 10)
    source_radius = 10

    sourcePositions = m2s.input.create_source_positions(
        source_azimuth_deg, source_colatitude_deg, source_radius)

    npt.assert_almost_equal(
        np.max(sourcePositions.get_sph()[..., 0]),
        np.max(source_azimuth_deg)/180*np.pi)
    npt.assert_almost_equal(
        np.min(sourcePositions.get_sph()[..., 0]),
        np.min(source_azimuth_deg)/180*np.pi)
    npt.assert_almost_equal(
        np.max(sourcePositions.get_sph()[..., 1]),
        np.max(source_colatitude_deg)/180*np.pi)
    npt.assert_almost_equal(
        np.min(sourcePositions.get_sph()[..., 1]),
        np.min(source_colatitude_deg)/180*np.pi)
    npt.assert_almost_equal(
        sourcePositions.get_sph()[..., 2], source_radius)


def test_write_scattering_parameter(source_coords_10deg, tmpdir):
    frequencies = pf.dsp.filter.fractional_octave_frequencies(
        3, (500, 5000))[0]
    path = os.path.join(
        m2s.utils.program_root(), '..',
        'tests', 'resources', 'mesh', 'sine_5k')
    sample_path = os.path.join(path, 'sample.stl')
    reference_path = os.path.join(path, 'reference.stl')
    receiver_delta_deg = 1
    receiver_radius = 5

    structural_wavelength = 0
    sample_diameter = 0.8
    model_scale = 2.5
    symmetry_azimuth = [90, 180]
    symmetry_rotational = False

    receiverPoints = pf.samplings.sph_equal_angle(
        receiver_delta_deg, receiver_radius)
    receiverPoints = receiverPoints[receiverPoints.get_sph()[..., 1] < np.pi/2]

    # execute
    m2s.input.write_scattering_project(
        project_path=tmpdir,
        frequencies=frequencies,
        sample_path=sample_path,
        reference_path=reference_path,
        receiver_coords=receiverPoints,
        source_coords=source_coords_10deg,
        structural_wavelength=structural_wavelength,
        model_scale=model_scale,
        sample_diameter=sample_diameter,
        symmetry_azimuth=symmetry_azimuth,
        symmetry_rotational=symmetry_rotational,
        )

    # test parameters
    f = open(os.path.join(tmpdir, 'parameters.json'))
    paras = json.load(f)
    source_list = [list(i) for i in list(source_coords_10deg.get_cart())]
    receiver_list = [list(i) for i in list(receiverPoints.get_cart())]
    parameters = {
        # project Info
        "project_title": 'scattering pattern',
        "mesh2scattering_path": m2s.utils.program_root(),
        "mesh2scattering_version": m2s.__version__,
        "bem_version": 'ML-FMM BEM',
        # Constants
        "speed_of_sound": float(346.18),
        "density_of_medium": float(1.1839),
        # Sample Information, post processing
        "structural_wavelength": structural_wavelength,
        "model_scale": model_scale,
        "sample_diameter": sample_diameter,
        # symmetry information
        "symmetry_azimuth": symmetry_azimuth,
        "symmetry_rotational": symmetry_rotational,
        # frequencies
        "num_frequencies": len(frequencies),
        "min_frequency": frequencies[0],
        "max_frequency": frequencies[-1],
        "frequencies": list(frequencies),
        # Source definition
        "source_type": 'Point source',
        "sources_num": len(source_list),
        "sources": source_list,
        # Receiver definition
        "receivers_num": len(receiver_list),
        "receivers": receiver_list,
    }
    npt.assert_array_almost_equal(paras['receivers'], parameters['receivers'])
    paras['receivers'] = parameters['receivers']
    npt.assert_equal(paras, parameters)
    # test folder structure
    assert os.path.isdir(os.path.join(tmpdir, 'sample'))
    assert os.path.isdir(os.path.join(tmpdir, 'reference'))
    assert os.path.isdir(os.path.join(tmpdir, 'sample', 'EvaluationGrids'))
    assert os.path.isdir(os.path.join(tmpdir, 'reference', 'EvaluationGrids'))
    assert os.path.isdir(os.path.join(tmpdir, 'sample', 'NumCalc'))
    assert os.path.isdir(os.path.join(tmpdir, 'reference', 'NumCalc'))
    assert os.path.isdir(os.path.join(tmpdir, 'sample', 'ObjectMeshes'))
    assert os.path.isdir(os.path.join(tmpdir, 'reference', 'ObjectMeshes'))
    assert os.path.isfile(os.path.join(
        tmpdir, 'sample', 'ObjectMeshes', 'sample.stl'))
    assert os.path.isfile(os.path.join(
        tmpdir, 'reference', 'ObjectMeshes', 'reference.stl'))

    # test sources
    for i in range(80):
        assert os.path.isdir(
            os.path.join(tmpdir, 'sample', 'NumCalc', f'source_{i+1}'))
    assert not os.path.isdir(
        os.path.join(tmpdir, 'sample', 'NumCalc', f'source_{81}'))
    for i in range(8):
        assert os.path.isdir(
            os.path.join(tmpdir, 'reference', 'NumCalc', f'source_{i+1}'))
    assert not os.path.isdir(
        os.path.join(tmpdir, 'reference', 'NumCalc', f'source_{9}'))
