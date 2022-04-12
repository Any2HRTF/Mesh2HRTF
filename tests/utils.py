"""Utilities to be used in testing"""
import numpy as np
import matplotlib.pyplot as plt
import os
import sofar as sf

# remove this if used in test file
# %matplotlib qt


def blender_paths(computer_id):
    """
    Define Blender paths for different computers for easier switching of
    test environments
    """

    if computer_id == 1:
        # matheson @ audio communication group
        base_dir = os.path.join(os.sep, 'home', 'matheson', 'Apps')
        blender_paths = [
            (os.path.join(base_dir, 'blender-2.83.10'),
             os.path.join('2.83', 'scripts', 'addons')),
            (os.path.join(base_dir, 'blender-2.91.0'),
             os.path.join('2.91', 'scripts', 'addons')),
            (os.path.join(base_dir, 'blender-2.93.8'),
             os.path.join('2.93', 'scripts', 'addons')),
            (os.path.join(base_dir, 'blender-3.1.2'),
             os.path.join('3.1', 'scripts', 'addons'))]
    elif computer_id == 2:
        # bruel @ audio communication group
        blender_paths = [
            ('/home/chef/Desktop/Applications/blender-2.93.8-linux-x64',
             '2.93/scripts/addons')]
    else:
        raise ValueError("Invalid computer id")

    return blender_paths


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
    plt.figure(figsize=(30/2.54, 8/2.54))

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

    for n, (p, title, clabel, crange) in enumerate(zip(
            [p_num, p_ana, p_num / p_ana],          # pressure arrays
            ["Numerical solution",                  # titles
             "Analytical solution",
             "Difference (max. deviation {:.2f} dB)".format(
                np.nanmax(np.abs(20*np.log10(np.abs(p_num / p_ana)))))],
            ["pressure", "pressure", "difference"],  # colorbar labels
            [range_a, range_a, range_b]              # range of colormap
            )):
        ax = plt.subplot(1, 3, n+1)
        sc = plt.scatter(x, y, c=20*np.log10(np.abs(p)),
                         edgecolors=None, s=1, vmax=crange[0],
                         vmin=crange[1], cmap="RdBu_r")
        cb = plt.colorbar(sc)
        plt.axis("off")
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
