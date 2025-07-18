"""Utilities to be used in testing"""
import os
import shutil
import numpy as np
import matplotlib.pyplot as plt
import sofar as sf

# set to correctly obtain paths from `blender_paths()` below:
COMPUTER_ID = 1


def blender_paths():
    """
    Define Blender paths for different computers for easier switching of
    test environments. First entry of tuple contains the path under which the
    `blender` executable is located. Second and third entry are relative to the
    first entry and contains the Blender addon and startup directories.
    """

    if COMPUTER_ID == 1:
        # bruel @ audio communication group
        blender_paths = [
            # earliest supported LTS version
            ('/home/bruel/Daten/Applications/blender-2.83.20-linux-x64/',
             '2.83/scripts/addons',
             '2.83/scripts/startup'),
            # latest LTS 3 version
            ('/home/bruel/Daten/Applications/blender-3.6.12-linux-x64/',
             '3.6/scripts/addons',
             '3.6/scripts/startup'),
            # latest LTS 4 version
            ('/home/bruel/Daten/Applications/blender-4.2.11-linux-x64/',
             '4.2/scripts/addons_core',
             '4.2/scripts/startup'),
            # latest version
            ('/home/bruel/Daten/Applications/blender-4.4.3-linux-x64/',
             '4.4/scripts/addons_core',
             '4.4/scripts/startup'),]
    elif COMPUTER_ID == 2:
        # panik macbook
        blender_paths = [
            # blender 4.1
            ('/Applications/Blender.app/Contents/MacOS/',
             '../Resources/4.1/scripts/addons',
             '../Resources/4.1/scripts/startup')
        ]
    else:
        raise ValueError("Invalid computer id")

    # check paths
    for path in blender_paths:
        for relative in [path[1], path[2]]:
            for p in [path[0], os.path.join(path[0], relative)]:
                if not os.path.isdir(p):
                    raise ValueError((
                        f"path {path} does not exist. Insert correct paths in"
                        " function blender_paths before testing."))

    return blender_paths


def install_blender_addons_and_scripts(
        blender_path, addon_path, addons, script_path, scripts,
        install_script):
    """
    Install Blender scripts and generate Python script to install Blender
    add-ons

    Parameters
    ----------
    blender_binary : str
        Blender binary for which the addons are refreshed
    addon_path : str
        Addon path under which the addons are saved
    addons : list
        List of tuples for installing the addons. The first entry of the tuple
        contains the path to the addon that is installed. The second entry
        contains the name of the addon for enabling it. The second entry is
        optional. Passing
    script_path : str
        Full path for saving the script.
    scripts : list
        List containing the full path to scripts. Will be copied to
        script_path. Path an empty list to install any scripts.
    install_script : str
        Full path and file name under which the generated python script to
        install addons will be saved.

        The script can be used to install Blender Add-Ons
        Option 1: run the script from within blender
        Option 2: Run blender in the background
        >>> subprocess.run(
        >>>     [os.path.join(blender_path, 'blender'), '--background',
        >>>      '--python', install_script,
        >>>      '--python', any_other_script.py],
        >>>     check=True, capture_output=True)
    """

    # generate script for installing the addons -------------------------------
    # initialize
    script = (
        "import bpy\n\n"
        "# set addon directory\n"
        '# Explicitly setting the addon dir worked in older Blender versions\n'
        '# But since the default addon dir is used, we dont really need this\n'
        "# bpy.context.preferences.filepaths.script_directory = "
        f"'{os.path.join(blender_path, addon_path)}'\n"
        "bpy.utils.refresh_script_paths()\n")
    # add code for installing addons
    for path, name in addons:
        script += (
            f"\n# install {os.path.basename(path)}\n"
            "bpy.ops.preferences.addon_install("
            f"overwrite=True, filepath='{path}')\n"
            f"bpy.ops.preferences.addon_enable(module='{name}')\n")

    with open(install_script, 'w') as file:
        file.writelines(script)

    # copy scripts to startup folder ------------------------------------------
    for script in scripts:
        shutil.copyfile(
            script,
            os.path.join(blender_path, script_path, os.path.basename(script)))


def write_blender_export_script(
        scriptPath, projectPath, programPath, addonFile, addonPath, params):
    """
    Write python script that exports a Mesh2HRTF project from Blender. The
    script can be used to export projects from the command line using

    ``blender --background "path/to/project.blend --python "scriptPath"``

    where project.blend is a blender Project ready for Mesh2HRTF export, i.e.,
    containing a mesh named Reference, with Left/Right ear materials, a point
    source, or a plane wave.


    Parameters
    ----------
    scriptPath : str
        Complete path and name for saving the export script to disk.
    projectPath : str
        Path to the folder into which the Mesh2HRTF project is exported
    programPath : str
        Path to the mesh2hrtf foler inside the Mesh2HRTF repository
    addonFile : str
        Path to the `mesh2input.py` Python addon for Blender
    addonPath : str
        Path under which the addon is to be installed (e.g. the folder
        3.1/scripts/addons for Blender 3.1 on Linux systems)
    params : dict
        Optional parameters for exporting the project, e.g.,
        ``{speedOfSound: "343", "Pictures": False}``. For a full list of
        parameters see `mesh2input.py` or the Mesh2HRTF export menue in
        Blender.
    """

    # parse export parameters to string
    export_args = ''
    for key, value in params.items():

        # check special parameter
        if key in ["speedOfSound", "densityOfMedium"] \
                and not isinstance(value, str):
            value = str(value)

        if type(value) is str:
            export_args += f'    {key}="{value}",\n'
        elif type(value) is bool:
            export_args += f'    {key}={value!s},\n'
        else:  # int or float
            export_args += f'    {key}={value},\n'

    # write export script to string
    script = (
        "import bpy\n\n"
        "# Define paths (done in testing by substituting the strings)\n\n"
        f"addonFile = '{addonFile}'\n"
        f"addonPath = '{addonPath}'\n\n"
        "# re-install export addon\n"
        '# Explicitly setting the addon dir worked in older Blender versions\n'
        '# But since the default addon dir is used, we dont really need this\n'
        "# bpy.context.preferences.filepaths.script_directory = addonPath\n"
        "bpy.utils.refresh_script_paths()\n"
        "bpy.ops.preferences.addon_install("
        "overwrite=True, filepath=addonFile)\n"
        "bpy.ops.preferences.addon_enable(module='mesh2input')\n\n"
        "# save Mesh2HRTF project\n"
        "bpy.ops.mesh2input.inp(\n"
        f"    filepath='{projectPath}',\n"
        f"    programPath='{programPath}',\n")

    if len(export_args):
        script += export_args[:-2]

    script += ")\n"

    # write export script to file
    with open(scriptPath, "w") as file:
        file.write(script)


def hrtf_sofa_to_numpy(path):
    """
    Read SOFA file with data type TF and return complex spectrum
    """

    sofa = sf.read_sofa(path)
    hrtf = sofa.Data_Real + 1j * sofa.Data_Imag

    return hrtf


def scatter_reference_vs_analytic(p_num, p_ana, x, y, range_a, range_b,
                                  boundary_condition, source, bem_method):
    """
    Plot numerical vs. analytical solution.

    Parameters
    ----------
    p_num, p_ana : numpy array
        complex pressure, normalized to max(abs(p)) = 1
    x, y : numpy array
        x and y coordinates
    range_a : (v_max, v_min)
        Range of the colormap in dB for plotting `p_num` and `p_ana`
    range_b : (v_max, v_min)
        Range of the colormap in dB for plotting `p_num / p_ana`
    boundary_condition, source: strings
        specifiying simulation test case
    """

    # Plot simulated and analytical data data
    if source == 'bothears':
        for iEar in range(2):
            plot_subfun(p_num[:, iEar], p_ana[:, iEar], range_a, range_b, x, y,
                        boundary_condition, source, bem_method, iEar)
    else:
        plot_subfun(p_num, p_ana, range_a, range_b, x, y, boundary_condition,
                    source, bem_method)


def plot_subfun(p_num, p_ana, range_a, range_b, x, y, boundary_condition,
                source, bem_method, iEar=0):

    _, axes = plt.subplots(1, 3, figsize=(30/2.54, 8/2.54))

    for n, (p, title, clabel, crange) in enumerate(zip(
            [p_num, p_ana, p_num / p_ana],          # pressure arrays
            ["Numerical solution",                  # titles
             "Analytical solution",
             "Difference: abs. mean (max) {:.2f} ({:.2f}) dB".format(
                np.nanmean(np.abs(20*np.log10(np.abs(p_num / p_ana)))),
                np.nanmax(np.abs(20*np.log10(np.abs(p_num / p_ana)))))],
            ["pressure", "pressure", "difference"],  # colorbar labels
            [range_a, range_a, range_b]              # range of colormap
            )):
        ax = axes[n]
        sc = ax.scatter(x, y, c=20*np.log10(np.abs(p)),
                        edgecolors=None, s=1, vmax=crange[0],
                        vmin=crange[1], cmap="RdBu_r")
        cb = plt.colorbar(sc)
        ax.axis("off")
        cb.set_label(clabel + ' in dB')
        ax.set_aspect("equal")
        ax.set_title(title)

    plt.tight_layout()
    if source == 'bothears':
        plt.savefig(os.path.dirname(__file__) +
                    "/resources/test_numcalc/analytical_references/" +
                    "comparisonplot_"+boundary_condition+"_"+source+"_" +
                    str(iEar+1)+"_"+bem_method+".jpg")
    else:
        plt.savefig(os.path.dirname(__file__) +
                    "/resources/test_numcalc/analytical_references/" +
                    "comparisonplot_"+boundary_condition+"_"+source+"_" +
                    bem_method+".jpg")

    plt.close()
