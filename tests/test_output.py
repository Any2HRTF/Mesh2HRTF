import os
import mesh2scattering as m2s
import numpy.testing as npt
import pyfar as pf


def test_import():
    from mesh2scattering import output
    assert output


def test_write_pattern():
    project_path = os.path.join(
        m2s.utils.repository_root(), "examples", "project")
    m2s.output.write_pattern(project_path)
    reference, source_coords_ref, receiver_coords_ref = pf.io.read_sofa(
        os.path.join(project_path, 'reference.pattern.sofa'))
    sample, source_coords, receiver_coords = pf.io.read_sofa(
        os.path.join(project_path, 'sample.pattern.sofa'))
    assert sample.cshape[0] == source_coords.csize
    assert sample.cshape[1] == receiver_coords.csize
    assert sample.cshape[0] == source_coords.csize
    assert sample.cshape[1] == receiver_coords.csize
    assert reference.cshape[0] == source_coords_ref.csize
    assert reference.cshape[1] == receiver_coords_ref.csize
    assert reference.cshape[0] == source_coords_ref.csize
    assert reference.cshape[1] == receiver_coords_ref.csize
    assert reference.cshape == sample.cshape
    npt.assert_equal(reference.frequencies, sample.frequencies)
