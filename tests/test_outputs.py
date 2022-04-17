import pytest
import shutil
import os
import glob
from tempfile import TemporaryDirectory
import mesh2hrtf as m2h

cwd = os.path.dirname(__file__)
data_shtf = os.path.join(cwd, 'resources', 'SHTF')


@pytest.mark.parametrize("savedir, generate_error", (
    [False, False], [False, True], [True, True]))
def test_outputs_to_hrtfs_minimum(savedir, generate_error):
    """
    Test without merging and inspecting files; using a single path with and
    without savedir
    """

    # copy required data to temporary directory (folder left, and right contain
    # two Mesh2HRTF projects each)
    tmp = TemporaryDirectory()
    shutil.copytree(data_shtf, os.path.join(tmp.name, "SHTF"))
    shutil.rmtree(os.path.join(tmp.name, "SHTF", "Output2HRTF"))
    shutil.copytree(os.path.join(tmp.name, "SHTF"),
                    os.path.join(tmp.name, "HRTF"))

    if generate_error:
        # clear one output file to generate the error
        file = os.path.join(
            tmp.name, "HRTF", "NumCalc", "source_1", "NC1-20.out")
        with open(file, "w") as file:
            file.write("\n")

    # process outputs
    savedir = None if not savedir else os.path.join(tmp.name, "output")
    if not generate_error:
        m2h.outputs_to_hrtfs(os.path.join(tmp.name, "*"))
    else:
        match = ("Detected issues in NumCalc output. Check report files in "
                 ".*\n.*HRTF")
        with pytest.raises(ValueError, match=match):
            m2h.outputs_to_hrtfs(os.path.join(tmp.name, "*"), savedir=savedir)

    # check output directories (only if not moved to savedir)
    if not savedir:
        for folder in ["HRTF", "SHTF"]:

            files = ["HRIR_FourPointHorPlane_r100cm.sofa",
                     "HRTF_FourPointHorPlane_r100cm.sofa",
                     "ObjectMesh_Reference.npz",
                     "report_source_1.csv",
                     "report_source_2.csv"]
            if generate_error and folder == "HRTF":
                files += ["report_issues.txt"]

            output = glob.glob(os.path.join(
                tmp.name, folder, "Output2HRTF", "*"))
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

        output = glob.glob(os.path.join(tmp.name, "output", "*"))
        output = [os.path.basename(o) for o in output]

        assert len(output) == len(files)
        for file in files:
            assert file in output

    # check issue report in savedir
    if savedir and generate_error:
        with open(os.path.join(
                tmp.name, "output", "report_issues.txt"), "r") as file:
            report = "\n".join(file.readlines())

        assert "Detected issues in NumCalc output" in report
        assert "HRTF" in report
        assert "SHTF" not in report


def test_outputs_to_hrtfs_full():
    """Test all options and passing multiple paths in tuple"""

    # copy required data to temporary directory (folder left, and right contain
    # two Mesh2HRTF projects each)
    tmp = TemporaryDirectory()
    os.mkdir(os.path.join(tmp.name, "left"))
    shutil.copytree(data_shtf, os.path.join(tmp.name, "left", "SHTF"))
    shutil.rmtree(os.path.join(tmp.name, "left", "SHTF", "Output2HRTF"))
    shutil.copytree(os.path.join(tmp.name, "left", "SHTF"),
                    os.path.join(tmp.name, "left", "HRTF"))
    shutil.copytree(os.path.join(tmp.name, "left"),
                    os.path.join(tmp.name, "right"))

    # process outputs
    m2h.outputs_to_hrtfs(
        (os.path.join(tmp.name, "left", "*"),
         os.path.join(tmp.name, "right", "*")),
        merge=True, inspect=True, pattern="HRIR",
        savedir=os.path.join(tmp.name, "output"))

    # check files in savedir
    output = glob.glob(os.path.join(tmp.name, "output", "*"))
    output = [os.path.basename(o) for o in output]
    files = ["HRTF_HRIR_FourPointHorPlane_r100cm.sofa",
             "HRTF_HRIR_FourPointHorPlane_r100cm_2D.pdf",
             "HRTF_HRIR_FourPointHorPlane_r100cm_3D.pdf",
             "SHTF_HRIR_FourPointHorPlane_r100cm.sofa",
             "SHTF_HRIR_FourPointHorPlane_r100cm_2D.pdf",
             "SHTF_HRIR_FourPointHorPlane_r100cm_3D.pdf"]
    assert len(output) == len(files)
    for file in files:
        assert file in output

    # check files in project directories
    for folder_1 in ["left", "right"]:
        for folder_2 in ["HRTF", "SHTF"]:
            files = glob.glob(os.path.join(
                tmp.name, folder_1, folder_2, "Output2HRTF", "*"))
            files = [os.path.basename(f) for f in files]

            assert len(files) == 2
            assert "report_source_1.csv" in files
            assert "report_source_2.csv" in files
