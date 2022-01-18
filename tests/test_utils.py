# %%
from matplotlib import cm
import scipy.io
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import os

# remove this if used in test file
# %matplotlib qt


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
    for n, (p, title, clabel, crange) in enumerate(zip(
            [p_num, p_ana, p_num / p_ana],           # pressure arrays
            ["Numerical solution",                   # titles
             "Analytical solution",
             "Difference (max. deviation {:.2f} dB)".format(
                 np.nanmax(np.abs(20*np.log10(np.abs(p_num / p_ana)))))],
            ["pressure", "pressure", "difference"],  # colorbar labels
            [range_a, range_a, range_b]              # range of the colormap
            )):
        ax = plt.subplot(1, 3, n+1)
        sc = plt.scatter(x, y, c=20*np.log10(np.abs(p)),
                         edgecolors=None, s=1, vmax=crange[0], vmin=crange[1],
                         cmap="RdBu_r")
        cb = plt.colorbar(sc)
        plt.axis("off")
        cb.set_label(clabel + ' in dB')
        ax.set_aspect("equal")
        ax.set_title(title)

    plt.tight_layout()
    plt.savefig(os.path.dirname(__file__) + "/test_numcalc_analytical_references/comparisonplot_" +
                boundary_condition+"_"+source+"_"+bem_method+".jpg")


# load data
# ref = (r"C:\Users\panik\Documents\Code\Third_party\mesh2hrtf-git\tests"
#        r"\test_numcalc_analytical_references\ref_rigid_plane.mat")

# ref = scipy.io.loadmat(ref)

# xyz = ref["XYZ"]
# p = ref["p_total"]

# p = p/np.max(np.abs(p[np.isfinite(p)]))

# # plot data
# scatter_reference_vs_analytic(p, p, xyz[:, 0], xyz[:, 1], (0, -10), (5, -5))

# save, e.g. to test_numcalc_analytical_references
# plt.savefig("something.pdf")

# %%
