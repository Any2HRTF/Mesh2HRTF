# %%
import pytest
import subprocess
import tempfile
import shutil
import os
import utils

# define and check paths to your Blender versions
blender_paths = utils.blender_paths(2)

# directory of this file and test data
base_dir = os.path.dirname(__file__)
grid_dir = os.path.join(
    base_dir, 'resources', 'test_blender_export', 'test_grids')
material_dir = os.path.join(
    base_dir, 'resources', 'test_blender_export', 'test_materials')


@pytest.mark.parametrize("blender_path, addon_path", blender_paths)
@pytest.mark.parametrize("blender_file_name, params, match_nc, match_o2hrtf", [
    # test default paramters - pictures disabled due to long rendering time
    ('test_export.blend',
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
     ["sourceType = 'Both ears'\n"
      "numSources = 2\n"
      "sourceCenter[0, :] = [-0.000087, -0.000000, 0.000002]\n"
      "sourceArea[0, 0] = 5.43588e-12\n"
      "sourceCenter[1, :] = [0.000087, 0.000000, 0.000002]\n"
      "sourceArea[1, 0] = 5.29446e-12",
      "reference = False", "computeHRIRs = False",
      "speedOfSound = 346.18", "densityOfAir = 1.1839"]),

    # test varying most parameters
    ('test_export.blend',
     {"title": "test title", "sourceType": "Point source",
      "method": "SL-FMM BEM", "evaluationGrids": "PlaneHorizontal",
      "speedOfSound": "300", "densityOfMedium": "1", "pictures": False},
     [["##\ntest title", "POINT SOURCES",
       "Main Parameters I\n2 45313 23132 0 0 2 1 1 0",
       "/EvaluationGrids/PlaneHorizontal/Nodes.txt",
       "Main Parameters IV\n300 1e+00"]],
     ["sourceType = 'Point source'\n"
      "numSources = 1\n"
      "sourceCenter[0, :] = [0.00020000000298023225, 0.0, 0.0]\n"
      "sourceArea[0, 0] = 1\n",
      "speedOfSound = 300", "densityOfAir = 1"]),

    # test unit m (default is mm)
    ('test_export.blend',
     {"sourceType": "Point source", "unit": "m", "pictures": False},
     [["POINT SOURCES\n0 0.20000000298023224 0.0 0.0 0.1 -1 0.0 -1"]],
     ["sourceCenter[0, :] = [0.20000000298023224, 0.0, 0.0]\n"]),

    # test remaining BEM method
    ('test_export.blend',
     {"sourceType": "Point source", "method": "BEM", "pictures": False},
     [["POINT SOURCES",
      "Main Parameters I\n2 24176 12092 0 0 2 1 0 0"]], []),

    # test point source
    ('test_export.blend',
     {"sourceType": "Point source", "pictures": False},
     [["POINT SOURCES"]],
     ["sourceType = 'Point source'"]),
    # test plane wave
    ('test_export.blend',
     {"sourceType": "Plane wave", "pictures": False},
     [["PLANE WAVE"]],
     ["sourceType = 'Plane wave'"]),
    # test left ear
    ('test_export.blend',
     {"sourceType": "Left ear", "pictures": False},
     [["BOUNDARY\n"
       "# Left ear velocity source\n"
       "ELEM 20479 TO 20479 VELO 0.1 -1 0.0 -1\n"
       "RETU\n"]],
     ["sourceType = 'Left ear'"]),
    # test right ear
    ('test_export.blend',
     {"sourceType": "Right ear", "pictures": False},
     [["BOUNDARY\n"
       "# Right ear velocity source\n"
       "ELEM 20478 TO 20478 VELO 0.1 -1 0.0 -1\n"
       "RETU"]],
     ["sourceType = 'Right ear'"]),
    # test both ears
    ('test_export.blend',
     {"sourceType": "Both ears", "pictures": False},
     [["BOUNDARY\n"
       "# Left ear velocity source\n"
       "ELEM 20479 TO 20479 VELO 0.1 -1 0.0 -1\n"
       "RETU\n"],
      ["BOUNDARY\n"
       "# Right ear velocity source\n"
       "ELEM 20478 TO 20478 VELO 0.1 -1 0.0 -1\n"
       "RETU"]],
     ["sourceType = 'Both ears'"]),

    # test built-in sound soft material
    ('test_export_soundsoft.blend',
     {"sourceType": "Point source", "pictures": False},
     [["BOUNDARY\n"
       "# Material: SoundSoft\n"
       "ELEM 0 TO 20479 PRES 1.0 1 1.0 2\n"
       "RETU"]], []),

    # test multiple built-in grids
    ('test_export.blend',
     {"evaluationGrids": "ARI;PlaneFrontal;PlaneMedian", "pictures": False},
     [["../../EvaluationGrids/ARI/Nodes.txt\n"
       "../../EvaluationGrids/PlaneFrontal/Nodes.txt\n"
       "../../EvaluationGrids/PlaneMedian/Nodes.txt"]], []),

    # test multiple external grids
    ('test_export.blend',
     {"evaluationGrids":
      (f"{os.path.join(grid_dir, 'test_grid_a')}; "
       f"{os.path.join(grid_dir, 'test_grid_b')}; "
       f"{os.path.join(grid_dir, 'test_grid_c')}"),
      "pictures": False},
     [["../../EvaluationGrids/test_grid_a/Nodes.txt\n"
       "../../EvaluationGrids/test_grid_b/Nodes.txt\n"
       "../../EvaluationGrids/test_grid_c/Nodes.txt"]], []),

    # test external materials
    ('test_export_external_materials.blend',
     {"sourceType": "Point source",
      "materialSearchPaths": material_dir,
      "pictures": False},
     [["BOUNDARY\n"
       "# Material: test_material_a\n"
       "ELEM 0 TO 12743 ADMI 1.0 1 1.0 2\n"
       "# Material: test_material_b\n"
       "ELEM 12744 TO 20479 PRES 1.0 3 1.0 4\n"
       "RETU\n"]], []),

    # test frequency vector option 'Num steps'
    ('test_export.blend',
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
       "0.000010 0.400000e+04 0.0\n"]], []),

    # test frequency vector option 'Step size'
    ('test_export.blend',
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
       "0.000003 0.400000e+04 0.0\n"]], []),

    # test Output2HRTF flags reference and computeHRIRs
    ('test_export.blend',
     {"reference": True, "computeHRIRs": True, "pictures": False},
     [[]],
     ["reference = True", "computeHRIRs = True"]),

    # test picture generation (look for .png)
    ('test_export.blend',
     {"pictures": True},
     [], [])
    ])
def test_blender_export(blender_path, addon_path, blender_file_name, params,
                        match_nc, match_o2hrtf):
    """
    Test the exportMesh2HRTF Blender plugin

    1. Copy test data
    2. Modify python script to test different export parameters
    3. Run the python script inside blender
    4. Test the export by comparing the NC.inp and Output2HRTF.py files against
       given references
    """

    # --- check path ---
    if not os.path.isdir(blender_path):
        raise ValueError("Blender path does not exist and must be configured")
    if not os.path.isdir(os.path.join(blender_path, addon_path)):
        raise ValueError("Addon path does not exist and must be configured")

    # --- Setup ---
    # create a temporary directory
    tmp = tempfile.TemporaryDirectory()

    # copy test directory
    shutil.copytree(
        os.path.join(base_dir, 'resources', 'test_blender_export'),
        os.path.join(tmp.name, 'project'))

    tmp_path = os.path.join(tmp.name, 'project')
    program_path = os.path.join(base_dir, "..", "mesh2hrtf")
    export_path = os.path.join(tmp.name, "export")
    blender_file_path = os.path.join(tmp_path, blender_file_name)
    python_file_path = os.path.join(tmp_path, 'blender_script.py')

    utils.write_blender_export_script(
        python_file_path, export_path,
        os.path.join(base_dir, "..", "mesh2hrtf"),
        os.path.join(program_path, "Mesh2Input", "exportMesh2HRTF.py"),
        os.path.join(blender_path, addon_path),
        params
        )

    # --- Exercise ---
    # run exportMesh2HRTF from Blender command line interface w/ Python script
    subprocess.run(
        [os.path.join(blender_path, 'blender'), '--background',
         blender_file_path, '--python', python_file_path],
        cwd=tmp_path, check=True, capture_output=True)

    # --- Verify ---
    # compare NC.inp against refrerence strings
    # (nested list to check NC.inp of multiple sources)
    for s in range(len(match_nc)):
        NCinp_filepath = os.path.join(
            export_path, "NumCalc", f'source_{s+1}', "NC.inp")

        with open(NCinp_filepath) as NCinp_file:
            NCinp_text = NCinp_file.read()

        for m in range(len(match_nc[s])):
            print(f"Searching {match_nc[s][m]} in {NCinp_filepath}")
            assert match_nc[s][m] in NCinp_text

    # compare Output2HRTF.py against reference strings
    with open(os.path.join(export_path, "Output2HRTF.py")) as o2hrtf_file:
        o2hrtf_text = o2hrtf_file.read()

    for m in range(len(match_o2hrtf)):
        print((f"Searching {match_o2hrtf[m]} in "
               f"{os.path.join(export_path, 'Output2HRTF.py')}"))
        assert match_o2hrtf[m] in o2hrtf_text

    # check if pictures were created or not
    if "pictures" in params:
        if params["pictures"]:  # pictures:True
            for az in [0, 45, 90, 135, 180, 225, 270, 315]:
                assert os.path.exists(os.path.join(
                    export_path, "Pictures", f'{az}_deg_azimuth.png'))
        else:  # pictures: False
            for az in [0, 45, 90, 135, 180, 225, 270, 315]:
                assert not os.path.exists(os.path.join(
                    export_path, "Pictures", f'{az}_deg_azimuth.png'))
