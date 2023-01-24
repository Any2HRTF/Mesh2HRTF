import pytest
import numpy.testing as npt
import subprocess
import tempfile
import os
import utils
import json
import numpy as np

# define and check paths to your Blender versions
blender_paths = utils.blender_paths(2)

# directory of this file and test data
base_dir = os.path.dirname(__file__)
data_dir = os.path.join(base_dir, 'resources', 'test_blender_export')


@pytest.mark.parametrize(
    "blender_path, addon_path, script_path", blender_paths)
@pytest.mark.parametrize("blender_file_name, params, match_nc, match_params", [
    # test default paramters - pictures disabled due to long rendering time
    (os.path.join(data_dir, 'test_export.blend'),
     {"pictures": False},
     [["##\nHead-Related Transfer Functions",
       "Controlparameter II\n1 200 0.000001 0.00e+00 1 0 0",
       "Main Parameters I\n2 24176 12092 0 0 2 1 4",
       "Main Parameters IV\n346.18 1.1839e+00",
       "../../EvaluationGrids/Default/Nodes.txt",
       "BOUNDARY\n"
       "# Left ear velocity source\n"
       "ELEM 20479 TO 20479 VELO 0.1 -1 0.0 -1\n"
       "RETU"],
      ["##\nHead-Related Transfer Functions",
       "Controlparameter II\n1 200 0.000001 0.00e+00 1 0 0",
       "Main Parameters I\n2 24176 12092 0 0 2 1 4",
       "Main Parameters IV\n346.18 1.1839e+00",
       "../../EvaluationGrids/Default/Nodes.txt",
       "BOUNDARY\n"
       "# Right ear velocity source\n"
       "ELEM 20478 TO 20478 VELO 0.1 -1 0.0 -1\n"
       "RETU"]],
     {"projectTitle": "Head-Related Transfer Functions",
      "BEM_Type": "ML-FMM BEM",
      "sourceType": "Both ears",
      "numSources": 2,
      "sourceCenter":
        [[-8.746444061399e-05, -3.090826794505e-11, 1.520398097352e-06],
         [8.746536821127e-05, 7.858034223318e-12, 1.500441343524e-06]],
      "sourceArea": [5.43588e-12, 5.29446e-12],
      "reference": False,
      "computeHRIRs": False,
      "speedOfSound": 346.18,
      "densityOfMedium": 1.1839}),

    # test varying most parameters
    (os.path.join(data_dir, 'test_export.blend'),
     {"title": "test title", "sourceType": "Point source",
      "method": "SL-FMM BEM", "evaluationGrids": "PlaneHorizontal",
      "speedOfSound": "300", "densityOfMedium": "1", "pictures": False},
     [["##\ntest title", "POINT SOURCES",
       "Main Parameters I\n2 45313 23132 0 0 2 1 1 0",
       "/EvaluationGrids/PlaneHorizontal/Nodes.txt",
       "Main Parameters IV\n300 1e+00"]],
     {"3D_SceneUnit": "mm",
      "sourceType": 'Point source',
      "numSources": 1,
      "sourceCenter": [0.00020000000298023225, 0.0, 0.0],
      "sourceArea": [1],
      "speedOfSound": 300,
      "densityOfMedium": 1}),

    # test unit m (default is mm)
    (os.path.join(data_dir, 'test_export.blend'),
     {"sourceType": "Point source", "unit": "m", "pictures": False},
     [["POINT SOURCES\n0 0.20000000298023224 0.0 0.0 0.1 -1 0.0 -1"]],
     {"sourceCenter": [0.20000000298023224, 0.0, 0.0],
      "3D_SceneUnit": "m"}),

    # test remaining BEM method
    (os.path.join(data_dir, 'test_export.blend'),
     {"sourceType": "Point source", "method": "BEM", "pictures": False},
     [["POINT SOURCES",
      "Main Parameters I\n2 24176 12092 0 0 2 1 0 0"]],
     {"BEM_Type": "BEM"}),

    # test point source
    (os.path.join(data_dir, 'test_export.blend'),
     {"sourceType": "Point source", "pictures": False},
     [["POINT SOURCES"]],
     {"sourceType": 'Point source'}),

    # test plane wave
    (os.path.join(data_dir, 'test_export.blend'),
     {"sourceType": "Plane wave", "pictures": False},
     [["PLANE WAVE"]],
     {"sourceType": 'Plane wave'}),

    # test left ear
    (os.path.join(data_dir, 'test_export.blend'),
     {"sourceType": "Left ear", "pictures": False},
     [["BOUNDARY\n"
       "# Left ear velocity source\n"
       "ELEM 20479 TO 20479 VELO 0.1 -1 0.0 -1\n"
       "RETU\n"]],
     {"sourceType": 'Left ear'}),

    # test right ear
    (os.path.join(data_dir, 'test_export.blend'),
     {"sourceType": "Right ear", "pictures": False},
     [["BOUNDARY\n"
       "# Right ear velocity source\n"
       "ELEM 20478 TO 20478 VELO 0.1 -1 0.0 -1\n"
       "RETU"]],
     {"sourceType": 'Right ear'}),

    # test both ears
    (os.path.join(data_dir, 'test_export.blend'),
     {"sourceType": "Both ears", "pictures": False},
     [["BOUNDARY\n"
       "# Left ear velocity source\n"
       "ELEM 20479 TO 20479 VELO 0.1 -1 0.0 -1\n"
       "RETU\n"],
      ["BOUNDARY\n"
       "# Right ear velocity source\n"
       "ELEM 20478 TO 20478 VELO 0.1 -1 0.0 -1\n"
       "RETU"]],
     {"sourceType": 'Both ears'}),

    # test built-in sound soft material
    (os.path.join(data_dir, 'test_export_soundsoft.blend'),
     {"sourceType": "Point source", "pictures": False},
     [["BOUNDARY\n"
       "# Material: SoundSoft\n"
       "ELEM 0 TO 20479 PRES 1.0 1 1.0 2\n"
       "RETU"]], {}),

    # test multiple built-in grids
    (os.path.join(data_dir, 'test_export.blend'),
     {"evaluationGrids": "ARI;PlaneFrontal;PlaneMedian", "pictures": False},
     [["../../EvaluationGrids/ARI/Nodes.txt\n"
       "../../EvaluationGrids/PlaneFrontal/Nodes.txt\n"
       "../../EvaluationGrids/PlaneMedian/Nodes.txt"]],
     {"evaluationGrids": ["ARI", "PlaneFrontal", "PlaneMedian"]}),

    # test multiple external grids
    (os.path.join(data_dir, 'test_export.blend'),
     {"evaluationGrids":
      (f"{os.path.join(data_dir, 'test_grids', 'test_grid_a')}; "
       f"{os.path.join(data_dir, 'test_grids', 'test_grid_b')}; "
       f"{os.path.join(data_dir, 'test_grids', 'test_grid_c')}"),
      "pictures": False},
     [["../../EvaluationGrids/test_grid_a/Nodes.txt\n"
       "../../EvaluationGrids/test_grid_b/Nodes.txt\n"
       "../../EvaluationGrids/test_grid_c/Nodes.txt"]],
     {"evaluationGrids": ["test_grid_a", "test_grid_b", "test_grid_c"]}),

    # test external materials
    (os.path.join(data_dir, 'test_export_external_materials.blend'),
     {"sourceType": "Point source",
      "materialSearchPaths": os.path.join(data_dir, 'test_materials'),
      "pictures": False},
     [["BOUNDARY\n"
       "# Material: test_material_a\n"
       "ELEM 0 TO 12743 ADMI 1.0 1 1.0 2\n"
       "# Material: test_material_b\n"
       "ELEM 12744 TO 20479 PRES 1.0 3 1.0 4\n"
       "RETU\n"]], {}),

    # test frequency vector option 'Num steps'
    (os.path.join(data_dir, 'test_export.blend'),
     {"sourceType": "Point source", "minFrequency": 2000, "maxFrequency": 4000,
      "frequencyVectorType": "Num steps", "frequencyVectorValue": 10,
      "pictures": False},
     [["## Controlparameter II\n"
       "1 10 0.000001 0.00e+00 1 0 0\n"
       "##\n"
       "## Load Frequency Curve \n"
       "0 11\n"
       "0.000000 0.000000e+00 0.0\n"
       "0.000001 0.200000e+04 0.0\n"
       "0.000002 0.222222e+04 0.0\n"
       "0.000003 0.244444e+04 0.0\n"
       "0.000004 0.266667e+04 0.0\n"
       "0.000005 0.288889e+04 0.0\n"
       "0.000006 0.311111e+04 0.0\n"
       "0.000007 0.333333e+04 0.0\n"
       "0.000008 0.355556e+04 0.0\n"
       "0.000009 0.377778e+04 0.0\n"
       "0.000010 0.400000e+04 0.0\n"]],
     {"numFrequencies": 10}),

    # test frequency vector option 'Step size'
    (os.path.join(data_dir, 'test_export.blend'),
     {"sourceType": "Point source", "minFrequency": 2000, "maxFrequency": 4000,
      "frequencyVectorType": "Step size", "frequencyVectorValue": 1000,
      "pictures": False},
     [["## Controlparameter II\n"
       "1 3 0.000001 0.00e+00 1 0 0\n"
       "##\n"
       "## Load Frequency Curve \n"
       "0 4\n"
       "0.000000 0.000000e+00 0.0\n"
       "0.000001 0.200000e+04 0.0\n"
       "0.000002 0.300000e+04 0.0\n"
       "0.000003 0.400000e+04 0.0\n"]],
     {"frequencyStepSize": 1000}),

    # test frequency vector option 'Numinal Ocatave'
    (os.path.join(data_dir, 'test_export.blend'),
     {"sourceType": "Point source", "minFrequency": 1250, "maxFrequency": 5000,
      "frequencyVectorType": "Nominal n-th octave", "frequencyVectorValue": 3,
      "pictures": False},
     [["## Controlparameter II\n"
       "1 7 0.000001 0.00e+00 1 0 0\n"
       "##\n"
       "## Load Frequency Curve \n"
       "0 8\n"
       "0.000000 0.000000e+00 0.0\n"
       "0.000001 0.125000e+04 0.0\n"
       "0.000002 0.160000e+04 0.0\n"
       "0.000003 0.200000e+04 0.0\n"
       "0.000004 0.250000e+04 0.0\n"
       "0.000005 0.315000e+04 0.0\n"
       "0.000006 0.400000e+04 0.0\n"
       "0.000007 0.500000e+04 0.0\n"]],
     {"numFrequencies": 7}),

    # test Output2HRTF flags reference and computeHRIRs
    (os.path.join(data_dir, 'test_export.blend'),
     {"reference": True, "computeHRIRs": True, "pictures": False},
     [[]],
     {"reference": True,
      "computeHRIRs": True}),

    # test picture generation (look for .png)
    (os.path.join(data_dir, 'test_export.blend'),
     {"pictures": True},
     [],
     {"exportPictures": True})
])
def test_blender_export(
        blender_path, addon_path, script_path, blender_file_name, params,
        match_nc, match_params):
    """
    Test the mesh2input Blender plugin

    1. Copy test data
    2. Do scripted Mesh2HRTF exports with varying parameters
    3.
    """

    # --- check path ---
    if not os.path.isdir(blender_path):
        raise ValueError("Blender path does not exist and must be configured")
    if not os.path.isdir(os.path.join(blender_path, addon_path)):
        raise ValueError("Addon path does not exist and must be configured")

    # --- Setup ---
    # create a temporary directory and write export script
    tmp = tempfile.TemporaryDirectory()

    utils.write_blender_export_script(
        os.path.join(tmp.name, 'blender_script.py'),
        tmp.name,
        os.path.join(base_dir, "..", "mesh2scattering"),
        os.path.join(base_dir, "..", "mesh2scattering", "Mesh2Input",
                     "mesh2input.py"),
        os.path.join(blender_path, addon_path),
        params
        )

    # --- Exercise ---
    # run mesh2input from Blender command line interface w/ Python script
    subprocess.run(
        [os.path.join(blender_path, 'blender'), '--background',
         blender_file_name, '--python',
         os.path.join(tmp.name, 'blender_script.py')],
        cwd=tmp.name, check=True, capture_output=True)

    # --- Verify ---
    # compare NC.inp against reference strings
    # (nested list to check NC.inp of multiple sources)
    for s in range(len(match_nc)):
        NCinp_filepath = os.path.join(
            tmp.name, "NumCalc", f'source_{s+1}', "NC.inp")

        with open(NCinp_filepath) as NCinp_file:
            NCinp_text = NCinp_file.read()

        for m in range(len(match_nc[s])):
            print(f"Searching {match_nc[s][m]} in {NCinp_filepath}")
            assert match_nc[s][m] in NCinp_text

    # compare parameters.json against reference
    with open(os.path.join(tmp.name, "parameters.json")) as params_json:
        act_params = json.load(params_json)

    for key, val in match_params.items():
        print(f"testing {key}")
        if isinstance(val, list):
            if isinstance(val[0], str):
                for v, a in zip(val, act_params[key]):
                    assert v == a
            else:
                npt.assert_array_almost_equal(val, act_params[key], decimal=12)
        else:
            assert val == act_params[key]

    # check if pictures were created or not
    if "pictures" in params:
        if params["pictures"]:  # pictures:True
            for az in [0, 45, 90, 135, 180, 225, 270, 315]:
                assert os.path.exists(os.path.join(
                    tmp.name, "Pictures", f'{az}_deg_azimuth.png'))
        else:  # pictures: False
            for az in [0, 45, 90, 135, 180, 225, 270, 315]:
                assert not os.path.exists(os.path.join(
                    tmp.name, "Pictures", f'{az}_deg_azimuth.png'))


@pytest.mark.parametrize(
  "frequencyType, minFrequency, maxFrequency, frequencyValue, frequencies", [
    ('Nominal n-th octave', 1250, 5000, 3,
     [1250, 1600, 2000, 2500, 3150, 4000, 5000]),
    ('Nominal n-th octave', 1240, 5004, 3,
     [1250, 1600, 2000, 2500, 3150, 4000, 5000]),
    ('Nominal n-th octave', 1240, 5004, 1,
     [2000, 4000]),
    ('Nominal n-th octave', 50, 20000, 1,
     [63, 125, 250, 500, 1000, 2000, 4000, 8000, 16000]),
    ('Exact n-th octave', 1000, 5000, 1,
     [1000, 2000, 4000]),
    ('Exact n-th octave', 20, 5000, 1,
     [31.25, 62.5, 125, 250, 500, 1000, 2000, 4000]),
    ('Exact n-th octave', 50, 20000, 1,
     [62.5, 125, 250, 500, 1000, 2000, 4000, 8000, 16000]),
    ('Exact n-th octave', 1000, 2000, 3,
     [1000, 1259.921, 1587.401, 2000]),
    ('Step size', 100, 200, 100, [100, 200]),
    ('Step size', 1000, 5000, 1000, [1000, 2000, 3000, 4000, 5000]),
    ('Num steps', 1000, 5000, 5, [1000, 2000, 3000, 4000, 5000]),
    ('Num steps', 37, 115, 5, [37., 56.5, 76., 95.5, 115.]),
  ])
@pytest.mark.parametrize(
    "blender_path, addon_path, script_path", blender_paths)
@pytest.mark.parametrize(
    "blender_file_name", [os.path.join(data_dir, "test_export.blend")])
def test_blender_export_frequencies(
        frequencyType, minFrequency, maxFrequency, frequencyValue,
        frequencies, blender_path, addon_path, script_path, blender_file_name):
    """
    Test the mesh2input Blender plugin

    1. Copy test data
    2. Do scripted Mesh2HRTF exports with varying parameters
    3.
    """

    # --- check path ---
    if not os.path.isdir(blender_path):
        raise ValueError("Blender path does not exist and must be configured")
    if not os.path.isdir(os.path.join(blender_path, addon_path)):
        raise ValueError("Addon path does not exist and must be configured")

    # --- Setup ---
    # create a temporary directory and write export script
    params = {
      "sourceType": "Point source",
      "minFrequency": minFrequency,
      "maxFrequency": maxFrequency,
      "frequencyVectorType": frequencyType,
      "frequencyVectorValue": frequencyValue,
      "pictures": False}

    tmp = tempfile.TemporaryDirectory()

    utils.write_blender_export_script(
        os.path.join(tmp.name, 'blender_script.py'),
        tmp.name,
        os.path.join(base_dir, "..", "mesh2scattering"),
        os.path.join(base_dir, "..", "mesh2scattering", "Mesh2Input",
                     "mesh2input.py"),
        os.path.join(blender_path, addon_path),
        params
        )

    # --- Exercise ---
    # run mesh2input from Blender command line interface w/ Python script
    subprocess.run(
        [os.path.join(blender_path, 'blender'), '--background',
         blender_file_name, '--python',
         os.path.join(tmp.name, 'blender_script.py')],
        cwd=tmp.name, check=True, capture_output=True)

    # --- Verify ---
    # compare parameters.json against reference
    with open(os.path.join(tmp.name, "parameters.json")) as params_json:
        act_params = json.load(params_json)

    if 'Exact' in frequencyType:
        npt.assert_array_almost_equal(
            np.array(act_params['frequencies']), np.array(frequencies),
            decimal=2)
    else:
        assert act_params['frequencies'] == frequencies
    assert act_params['minFrequency'] == frequencies[0]
    assert act_params['maxFrequency'] == frequencies[-1]
    assert act_params['numFrequencies'] == len(frequencies)
    if frequencyType in ('Step size', 'Num steps'):
        assert act_params['frequencyStepSize'] == \
          frequencies[1] - frequencies[0]
    else:
        assert act_params['frequencyStepSize'] == 0
