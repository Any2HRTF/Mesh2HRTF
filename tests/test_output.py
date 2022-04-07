# %%
import subprocess
from tempfile import TemporaryDirectory
import shutil
import os
import sofar as sf


data_shtf = os.path.join(
    os.path.dirname(__file__), 'resources', 'SHTF')


def test_Output2HRTF_Main_and_Output2HRTF():
    """
    Run Output2HRTF.py script to do a round trip test
    """

    # copy test data to new directory and delete output data
    tmp = TemporaryDirectory()
    tmp_shtf = os.path.join(tmp.name, "SHTF")
    shutil.copytree(data_shtf, tmp_shtf)
    shutil.rmtree(os.path.join(tmp_shtf, "Output2HRTF"))

    # run Output2HRTF.py
    subprocess.run(["python", "Output2HRTF.py"], cwd=tmp_shtf, check=True)

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
