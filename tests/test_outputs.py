import pytest
import shutil
import os
from os.path import join
from glob import glob
from tempfile import TemporaryDirectory
import mesh2hrtf as m2h

cwd = os.path.dirname(__file__)
data_shtf = join(cwd, 'resources', 'SHTF')


@pytest.mark.parametrize("savedir, generate_error", (
    [False, False], [False, True], [True, True]))
def test_outputs_to_hrtfs_minimum(savedir, generate_error):
    """
    Test without merging and inspecting files; using a single path with and
    without savedir
    """

    # copy required data to temporary directory
    tmp = TemporaryDirectory()
    shutil.copytree(data_shtf, join(tmp.name, "SHTF"))
    shutil.rmtree(join(tmp.name, "SHTF", "Output2HRTF"))
    shutil.copytree(join(tmp.name, "SHTF"), join(tmp.name, "HRTF"))

    if generate_error:
        # clear one output file to generate the error
        file = join(tmp.name, "HRTF", "NumCalc", "source_1", "NC1-20.out")
        with open(file, "w") as file:
            file.write("\n")

    # process outputs
    savedir = None if not savedir else join(tmp.name, "output")
    if not generate_error:
        m2h.process_multiple_outputs2hrtf(join(tmp.name, "*"))
    else:
        match = ("Detected issues in NumCalc output. Check report files in "
                 ".*\n.*HRTF")
        with pytest.raises(ValueError, match=match):
            m2h.process_multiple_outputs2hrtf(
                join(tmp.name, "*"), savedir=savedir)

    # check output directories (only if not moved to savedir)
    if not savedir:
        for folder in ["HRTF", "SHTF"]:

            files = ["HRIR_FourPointHorPlane_r100cm.sofa",
                     "HRTF_FourPointHorPlane_r100cm.sofa",
                     "report_source_1.csv",
                     "report_source_2.csv"]
            if generate_error and folder == "HRTF":
                files += ["report_issues.txt"]

            output = glob(join(tmp.name, folder, "Output2HRTF", "*"))
            output = [os.path.basename(o) for o in output]

            assert len(output) == len(files)
            for file in files:
                assert file in output
    else:
        files = ["HRTF_HRIR_FourPointHorPlane_r100cm.sofa",
                 "HRTF_HRTF_FourPointHorPlane_r100cm.sofa",
                 "SHTF_HRIR_FourPointHorPlane_r100cm.sofa",
                 "SHTF_HRTF_FourPointHorPlane_r100cm.sofa"]
        if generate_error:
            files += ["report_issues.txt"]

        output = glob(join(tmp.name, "output", "*"))
        output = [os.path.basename(o) for o in output]

        assert len(output) == len(files)
        for file in files:
            assert file in output

    # check issue report in savedir
    if savedir and generate_error:
        with open(join(tmp.name, "output", "report_issues.txt"), "r") as file:
            report = "\n".join(file.readlines())

        assert "Detected issues in NumCalc output" in report
        assert "HRTF" in report
        assert "SHTF" not in report


def test_outputs_to_hrtfs_full():
    """Test all options and passing multiple paths in tuple"""

    # copy required data to temporary directory (folder left, and right contain
    # two Mesh2HRTF projects each)
    tmp = TemporaryDirectory()
    os.mkdir(join(tmp.name, "left"))
    shutil.copytree(data_shtf, join(tmp.name, "left", "SHTF"))
    shutil.rmtree(join(tmp.name, "left", "SHTF", "Output2HRTF"))
    shutil.copytree(join(tmp.name, "left", "SHTF"),
                    join(tmp.name, "left", "HRTF"))
    shutil.copytree(join(tmp.name, "left"), join(tmp.name, "right"))

    # process outputs
    m2h.process_multiple_outputs2hrtf(
        (join(tmp.name, "left", "*"), join(tmp.name, "right", "*")),
        merge=True, inspect=True, pattern="HRIR",
        savedir=join(tmp.name, "output"))

    # check files in savedir
    output = glob(join(tmp.name, "output", "*"))
    output = [os.path.basename(o) for o in output]
    files = ["HRTF_HRIR_FourPointHorPlane_r100cm.sofa",
             "HRTF_HRIR_FourPointHorPlane_r100cm_2D.pdf",
             "HRTF_HRIR_FourPointHorPlane_r100cm_3D.jpeg",
             "SHTF_HRIR_FourPointHorPlane_r100cm.sofa",
             "SHTF_HRIR_FourPointHorPlane_r100cm_2D.pdf",
             "SHTF_HRIR_FourPointHorPlane_r100cm_3D.jpeg"]
    assert len(output) == len(files)
    for file in files:
        assert file in output

    # check files in project directories
    for folder_1 in ["left", "right"]:
        for folder_2 in ["HRTF", "SHTF"]:
            files = glob(join(
                tmp.name, folder_1, folder_2, "Output2HRTF", "*"))
            files = [os.path.basename(f) for f in files]

            assert len(files) == 2
            assert "report_source_1.csv" in files
            assert "report_source_2.csv" in files


@pytest.mark.parametrize("boundary,grid", [
    (True, True), (True, False), (False, True)])
def test_purge_outputs_numcalc_data(boundary, grid):
    """Test purging the raw NumCalc output"""

    # copy required data to temporary directory
    tmp = TemporaryDirectory()
    shutil.copytree(data_shtf, join(tmp.name, "SHTF"))

    m2h.remove_outputs(join(tmp.name, "*"), boundary, grid)

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

    m2h.remove_outputs(
        join(tmp.name, "*"), hrtf=hrtf, vtk=vtk, reports=reports)

    assert os.path.isfile(join(folder, "HRTF_FourPointHorPlane_r100cm.sofa")) \
        == (not hrtf)

    assert os.path.isdir(join(folder, "vtk")) \
        == (not vtk)

    assert os.path.isfile(join(folder, "report_source_1.csv")) == \
        (not reports)

    assert os.path.isfile(join(folder, "report_source_2.csv")) == \
        (not reports)
