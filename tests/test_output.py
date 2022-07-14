import pytest
import numpy as np
import numpy.testing as npt
from tempfile import TemporaryDirectory
import shutil
import os
import glob
import json
import pyfar as pf
import sofar as sf
import mesh2hrtf as m2h

cwd = os.path.dirname(__file__)
data_shtf = os.path.join(cwd, 'resources', 'SHTF')
data_nc = os.path.join(cwd, 'resources', 'nc.out')
data_grids = os.path.join(cwd, 'resources', 'evaluation_grids')
data_sofa = os.path.join(cwd, 'resources', 'SOFA_files')


@pytest.mark.parametrize("num_sources", ([1], [2]))
def test_output_two_hrtf_and_Output2HRTF(num_sources):
    """
    Run output_to_hrtf to do a round trip test:

    - does output_to_hrtf run without errors for projects with 1 and 2
      sources
    - are the report_source_*.csv files written correctly (2 sources only)
    - are the SOFA files written correctly (2 sources only)
    """

    # copy test data to new directory and delete output data
    tmp = TemporaryDirectory()
    tmp_shtf = os.path.join(tmp.name, "SHTF")
    shutil.copytree(data_shtf, tmp_shtf)
    shutil.rmtree(os.path.join(tmp_shtf, "Output2HRTF"))

    # manipulate parameters
    if num_sources == 1:
        with open(os.path.join(tmp_shtf, "parameters.json"), "r") as f:
            params = json.load(f)

        params["sourceType"] = "Left ear"
        params["numSources"] = 1
        params["sourceCenter"] = params["sourceCenter"][0]
        params["sourceArea"] = [params["sourceArea"][0]]

        with open(os.path.join(tmp_shtf, "parameters.json"), "w") as f:
            json.dump(params, f, indent=4)

    # run output_to_hrtf
    m2h.output_to_hrtf(tmp_shtf)

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

        # test data entries with tolerance
        # (results differ across operating systems)
        if sofa.startswith("HRTF"):
            npt.assert_allclose(test.Data_Real, ref.Data_Real, rtol=1e-5)
            npt.assert_allclose(test.Data_Imag, ref.Data_Imag, rtol=1e-5)
        else:
            npt.assert_allclose(test.Data_IR, ref.Data_IR, rtol=1e-1)

        # test remaining entries
        ignore = ["Data_Real", "Data_Imag", "Data_IR", "GLOBAL_APIVersion"]
        for key, value in test.__dict__.items():
            if key.startswith("_") or "Date" in key or key in ignore:
                continue

            print(f"{sofa}: {key}")

            if isinstance(value, np.ndarray):
                npt.assert_allclose(
                    value, getattr(ref, key), atol=1e-6, rtol=1e-3)
            else:
                assert value == getattr(ref, key)


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

    for file in glob.glob(os.path.join(tmp_shtf, "Output2HRTF", "*.pdf")):
        os.remove(file)
    for file in glob.glob(os.path.join(tmp_shtf, "Output2HRTF", "*.jpeg")):
        os.remove(file)

    # create plots
    m2h.inspect_sofa_files(tmp_shtf, pattern, plot=plot)

    grid = "FourPointHorPlane_r100cm"

    # check if the correct files exist and are missing
    for file in created:
        file = os.path.join(tmp_shtf, "Output2HRTF", file.replace("*", grid))
        extension = ".pdf" if "2D" in file else ".jpeg"
        assert os.path.isfile(file + extension)

    for file in not_created:
        file = os.path.join(tmp_shtf, "Output2HRTF", file.replace("*", grid))
        extension = ".pdf" if "2D" in file else ".jpeg"
        assert not os.path.isfile(file + extension)


def test_compute_hrir_custom_sampling_rate():
    """Test compute HRIR with custom sampling rate"""

    # test with default (test file with constant spectrum of ones)
    sofa = m2h.compute_hrir(
        os.path.join(data_sofa, "HRTF_test_max_freq_24k.sofa"), 40)
    hrir = pf.Signal(sofa.Data_IR, sofa.Data_SamplingRate)

    assert sofa.GLOBAL_SOFAConventions == "SimpleFreeFieldHRIR"
    assert hrir.frequencies[-1] == 24000
    assert hrir.sampling_rate == 48000
    npt.assert_almost_equal(np.abs(hrir.freq_raw), np.ones_like(hrir.freq_raw))

    # test with valid sampling rate (test file with constant spectrum of ones)
    sofa = m2h.compute_hrir(
        os.path.join(data_sofa, "HRTF_test_max_freq_24k.sofa"), 40, 44100)
    hrir = pf.Signal(sofa.Data_IR, sofa.Data_SamplingRate)

    assert sofa.GLOBAL_SOFAConventions == "SimpleFreeFieldHRIR"
    assert hrir.frequencies[-1] == 22050
    assert hrir.sampling_rate == 44100
    npt.assert_almost_equal(np.abs(hrir.freq_raw), np.ones_like(hrir.freq_raw))

    # test with invalid sampling rate
    with pytest.raises(ValueError, match="sampling rate is invalid"):
        sofa = m2h.compute_hrir(
            os.path.join(data_sofa, "HRTF_test_max_freq_24k.sofa"), 40, 44110)


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
    m2h.merge_sofa_files((data_shtf, data_shtf), pattern, tmp.name)

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
    shutil.copyfile(os.path.join(data_nc, "parameters.json"),
                    os.path.join(tmp.name, "parameters.json"))
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


@pytest.mark.parametrize("n_dim", [3, 2])
@pytest.mark.parametrize("coordinates,show", [[False, True], [True, False]])
def test_read_and_write_evaluation_grid(n_dim, coordinates, show):

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
    m2h.write_evaluation_grid(points, os.path.join(tmp.name, "test"),
                              discard=discard, show=show)

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
    coordinates = m2h.read_evaluation_grid(os.path.join(tmp.name, "test"))

    # check grid
    assert isinstance(coordinates, pf.Coordinates)
    npt.assert_equal(coordinates.get_cart(), points)


@pytest.mark.parametrize("frequency_steps,dB", (
    [None, True], [[1, 2], True], [[1, 1], False]
))
def test_export_to_vtk(frequency_steps, dB):

    # copy test data
    tmp = TemporaryDirectory()
    cwd = os.path.join(tmp.name, "SHTF")
    shutil.copytree(data_shtf, cwd)

    # export to vtk
    m2h.export_to_vtk(cwd, frequency_steps=frequency_steps, dB=dB)

    # check results
    if frequency_steps is None:
        frequency_steps = [1, 60]

    prefix = "db_frequency_step_" if dB else "lin_frequency_step_"

    for ff in range(frequency_steps[0], frequency_steps[1]+1):

        file = os.path.join(
            "Output2HRTF", "Reference_vtk", f"{prefix}{ff}.vtk")

        # check if all files are there
        assert os.path.isfile(os.path.join(cwd, file))

        # test file content against reference
        # (references only exist for steps 1 and 2 to save space)
        if ff == 1 or (ff == 2 and 2 in frequency_steps):
            with open(os.path.join(data_shtf, file), "r") as f:
                ref = "".join(f.readlines())
            with open(os.path.join(cwd, file), "r") as f:
                test = "".join(f.readlines())
            assert test == ref
