import pytest
import numpy.testing as npt
import subprocess
from tempfile import TemporaryDirectory
import shutil
import os
import glob
import sofar as sf
import mesh2hrtf as m2h

cwd = os.path.dirname(__file__)
data_shtf = os.path.join(cwd, 'resources', 'SHTF')
data_nc = os.path.join(cwd, 'resources', 'nc.out')


@pytest.mark.parametrize("num_sources", ([1], [2]))
def test_Output2HRTF_Main_and_Output2HRTF(num_sources):
    """
    Run Output2HRTF.py script to do a round trip test:

    - does Output2HRTF_Main.py run without errors for projects with 1 and 2
      sources
    - are the report_source_*.csv files written correctly (2 sources only)
    - are the SOFA files written correctly (2 sources only)
    """

    # copy test data to new directory and delete output data
    tmp = TemporaryDirectory()
    tmp_shtf = os.path.join(tmp.name, "SHTF")
    shutil.copytree(data_shtf, tmp_shtf)
    shutil.rmtree(os.path.join(tmp_shtf, "Output2HRTF"))

    # manipulate output script
    if num_sources == 1:
        with open(os.path.join(tmp_shtf, "Output2HRTF.py"), "r") as f:
            script = "".join(f.readlines())

        script.replace("numSources = 2", "numSources = 1")
        assert "numSources = 1" in script

        with open(os.path.join(tmp_shtf, "Output2HRTF.py"), "w") as f:
            f.write(script)

    # run Output2HRTF.py
    subprocess.run(["python", "Output2HRTF.py"], cwd=tmp_shtf, check=True)

    if num_sources == 1:
        return

    # compare reports
    reports = ["report_source_1.csv", "report_source_2.csv"]
    for report in reports:

        with open(os.path.join(data_shtf, "Output2HRTF", report), "r") as r:
            ref = r.readlines()
        with open(os.path.join(tmp_shtf, "Output2HRTF", report), "r") as r:
            test = r.readlines()

        assert "".join(test) == "".join(ref)

    # compare sofa files
    sofas = ["HRTF_FourPointHorPlane_r100cm.sofa",
             "HRIR_FourPointHorPlane_r100cm.sofa"]
    for sofa in sofas:

        ref = sf.read_sofa(os.path.join(data_shtf, "Output2HRTF", sofa))
        test = sf.read_sofa(os.path.join(tmp_shtf, "Output2HRTF", sofa))

        assert sf.equals(test, ref, exclude="DATE")


@pytest.mark.parametrize("pattern,plot,created,not_created", (
    [None, None, ["HRIR_*_2D", "HRIR_*_3D", "HRTF_*_2D", "HRTF_*_3D"], []],
    ["HRIR", None, ["HRIR_*_2D", "HRIR_*_3D"], ["HRTF_*_2D", "HRTF_*_3D"]],
    ["HRIR", "2D", ["HRIR_*_2D"], ["HRIR_*_3D", "HRTF_*_2D", "HRTF_*_3D"]],
    ["HRIR", "3D", ["HRIR_*_3D"], ["HRIR_*_2D", "HRTF_*_2D", "HRTF_*_3D"]]
))
def test_inspect_sofa_files_single_project(
        pattern, plot, created, not_created):
    """
    Test if inspect_sofa_files creates the correct plots for a single project.
    Note: Not all options for reading from and saving to different directories
          are tested at the moment.
    """

    # copy test data to new directory and delete all plots
    tmp = TemporaryDirectory()
    tmp_shtf = os.path.join(tmp.name, "SHTF")
    shutil.copytree(data_shtf, tmp_shtf)

    for pdf in glob.glob(os.path.join(tmp_shtf, "Output2HRTF", "*.pdf")):
        os.remove(pdf)

    # create plots
    m2h.inspect_sofa_files(tmp_shtf, pattern, plot=plot)

    grid = "FourPointHorPlane_r100cm"

    # check if the correct files exist and are missing
    for pdf in created:
        pdf = os.path.join(tmp_shtf, "Output2HRTF", pdf.replace("*", grid))
        assert os.path.isfile(pdf + ".pdf")

    for pdf in not_created:
        pdf = os.path.join(tmp_shtf, "Output2HRTF", pdf.replace("*", grid))
        assert not os.path.isfile(pdf + ".pdf")


@pytest.mark.parametrize("pattern", (None, "HRIR", "HRTF"))
def test_merge_sofa_files(pattern):
    """
    Test if merge_sofa_files creates the correct files.
    Note: Not all options for reading from and saving to different directories
          are tested at the moment.
    """

    grid = "FourPointHorPlane_r100cm"
    tmp = TemporaryDirectory()

    # merge two identical files
    m2h.merge_sofa_files(data_shtf, data_shtf, pattern, tmp.name)

    # check merged files
    pattern = ["HRTF", "HRIR"] if not pattern else [pattern]

    for p in pattern:
        ref = sf.read_sofa(os.path.join(
            data_shtf, "Output2HRTF", f"{p}_{grid}.sofa"))
        test = sf.read_sofa(os.path.join(tmp.name, f"{p}_{grid}_merged.sofa"))

        assert test.get_dimension("R") == 2 * ref.get_dimension("R")

        # test receiver positions
        npt.assert_equal(test.ReceiverPosition[:2], ref.ReceiverPosition)
        npt.assert_equal(test.ReceiverPosition[2:], ref.ReceiverPosition)

        # check frequency data
        if p == "HRTF":
            npt.assert_equal(test.Data_Real[:, :2], ref.Data_Real)
            npt.assert_equal(test.Data_Real[:, 2:], ref.Data_Real)

            npt.assert_equal(test.Data_Imag[:, :2], ref.Data_Imag)
            npt.assert_equal(test.Data_Imag[:, 2:], ref.Data_Imag)
        # check time data
        if p == "HRIR":
            npt.assert_equal(test.Data_IR[:, :2], ref.Data_IR)
            npt.assert_equal(test.Data_IR[:, 2:], ref.Data_IR)


@pytest.mark.parametrize("folders,issue,errors,nots", (
    # no issues single NC.inp file
    [["case_0"], False, [], []],
    # issues in NC.inp that are corrected by second file NC1-1.inp
    [["case_4"], False, [], []],
    # missing frequencies
    [["case_1"], True,
     ["Frequency steps that were not calculated:\n59, 60"], []],
    # convergence issues
    [["case_2"], True,
     ["Frequency steps that did not converge:\n18, 42"], []],
    # input/mesh issues
    [["case_3"], True,
     ["Frequency steps that were not calculated:\n59, 60",
      "Frequency steps with bad input:\n58"], []],
    # no isses in source 1 but issues in source 2
    [["case_0", "case_1"], True,
     ["Detected issues for source 2",
      "Frequency steps that were not calculated:\n59, 60"],
     ["Detected issues for source 1"]]
))
def test_project_report(folders, issue, errors, nots):
    """Test issues found by the project report"""

    # create fake project structure
    tmp = TemporaryDirectory()
    os.mkdir(os.path.join(tmp.name, "NumCalc"))
    os.mkdir(os.path.join(tmp.name, "Output2HRTF"))
    shutil.copyfile(os.path.join(data_nc, "Info.txt"),
                    os.path.join(tmp.name, "Info.txt"))
    for ff, folder in enumerate(folders):
        shutil.copytree(os.path.join(data_nc, folder),
                        os.path.join(tmp.name, "NumCalc", f"source_{ff + 1}"))

    # run the project report
    issues, report = m2h.project_report(tmp.name)

    # test the output
    assert issues is issue
    for error in errors:
        assert error in report
    for no in nots:
        assert no not in report
    if issue:
        assert os.path.isfile(os.path.join(
            tmp.name, "Output2HRTF", "report_issues.txt"))
        assert ("For more information check Output2HRTF/report_source_*.csv "
                "and the NC*.inp files located at NumCalc/source_*") in report
    else:
        assert not os.path.isfile(os.path.join(
            tmp.name, "Output2HRTF", "report_issues.txt"))
