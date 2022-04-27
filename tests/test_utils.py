import mesh2hrtf as m2h
import pytest
import os
from tempfile import TemporaryDirectory


def test_write_boundary_condition():
    # test write boundary condition with default values

    tmp = TemporaryDirectory()
    filename = os.path.join(tmp.name, "test_material.csv")

    # write data
    m2h.write_boundary_condition(
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
def test_write_boundary_condition_kind(kind, check_kind):
    # test if the kind of boundary condition is written correctly

    tmp = TemporaryDirectory()
    filename = os.path.join(tmp.name, "test_material.csv")

    # write data
    m2h.write_boundary_condition(
        filename, kind, [100, 200], [1 + 0j, 1.5 + 0.5j])

    # read and check data
    with open(filename, "r") as f_id:
        file = f_id.readlines()

    assert f"{check_kind[0]}\n" in file
    assert f"{check_kind[1]}\n" not in file
    assert f"{check_kind[2]}\n" not in file


def test_write_boundary_condition_comment():
    # test if the comment is written

    tmp = TemporaryDirectory()
    filename = os.path.join(tmp.name, "test_material.csv")
    comment = "Weird, random data"

    # write data
    m2h.write_boundary_condition(
        filename, "pressure", [100, 200], [1 + 0j, 1.5 + 0.5j], comment)

    # read and check data
    with open(filename, "r") as f_id:
        file = f_id.readlines()

    assert file[0] == "# " + comment + "\n"
    assert file[1] == "#\n"
    assert file[2] == "# Keyword to define the boundary condition:\n"
