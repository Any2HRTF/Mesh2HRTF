import pytest
import shutil
import os
from os.path import join
from glob import glob
from tempfile import TemporaryDirectory
import mesh2scattering as m2s

cwd = os.path.dirname(__file__)
data_shtf = join(cwd, 'resources', 'SHTF')


@pytest.mark.parametrize("boundary,grid", [
    (True, True), (True, False), (False, True)])
def test_purge_outputs_numcalc_data(boundary, grid):
    """Test purging the raw NumCalc output"""

    # copy required data to temporary directory
    tmp = TemporaryDirectory()
    shutil.copytree(data_shtf, join(tmp.name, "SHTF"))

    m2s.remove_outputs(join(tmp.name, "*"), boundary, grid)

    for source in glob(join(tmp.name, "SHTF", "NumCalc", "source_*")):
        if boundary and grid:
            assert not os.path.isdir(join(source, "be.out"))
        elif boundary:
            assert os.path.isdir(join(source, "be.out"))
            for be in glob(join(source, "be.out", "be.*")):
                assert glob(join(be, "*Boundary")) == []
        elif grid:
            assert os.path.isdir(join(source, "be.out"))
            for be in glob(join(source, "be.out", "be.*")):
                assert glob(join(be, "*EvalGrid")) == []


@pytest.mark.parametrize("hrtf,vtk,reports", [
    (False, True, False), (True, False, True)])
def test_purge_outputs_output_data(hrtf, vtk, reports):
    """Test purging the processed data in Output2HRTF"""

    # copy required data to temporary directory
    tmp = TemporaryDirectory()
    shutil.copytree(data_shtf, join(tmp.name, "SHTF"))
    folder = join(tmp.name, "SHTF", "Output2HRTF")

    m2s.remove_outputs(
        join(tmp.name, "*"), hrtf=hrtf, vtk=vtk, reports=reports)

    assert os.path.isfile(join(folder, "HRTF_FourPointHorPlane_r100cm.sofa")) \
        == (not hrtf)

    assert os.path.isdir(join(folder, "vtk")) \
        == (not vtk)

    assert os.path.isfile(join(folder, "report_source_1.csv")) == \
        (not reports)

    assert os.path.isfile(join(folder, "report_source_2.csv")) == \
        (not reports)
