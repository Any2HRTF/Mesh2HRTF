# %%
import mesh2scattering as m2s
import pyfar as pf
# import matplotlib.pyplot as plt
import os

parameters = (
    [None, "minimum", "equal"],                 # default parameter
    [1, "minimum", "equal"],                    # test smoothing
    [None, "linear", "equal"],                  # test linear phase filters
    [None, "zero", "equal"],                    # test zero phase filters
    [None, "minimum", 'voronoi'],               # test voronoi weights
    [None, "minimum", [1, 0, 0, 0, 0, 0]])      # test custom weights

for smooth_fractions, phase, weights in parameters:

    if isinstance(weights, list):
        name = f"compute_dtfs_{smooth_fractions}_{phase}_custom"
    else:
        name = f"compute_dtfs_{smooth_fractions}_{phase}_{weights}"

    print(name)

    dtf, dftf_inverse = m2s.compute_dtfs(
        os.path.join("..", "resources", "SOFA_files", "HRIR_6_points.sofa"),
        smooth_fractions, phase, weights)

    # plots used for manual verification of the results

    # pf.plot.time_freq(dftf_inverse)
    # plt.savefig(name + "_filter.png")
    # plt.close("all")

    # pf.plot.time_freq(pf.io.convert_sofa(dtf)[0])
    # plt.savefig(name + "_dtfs.png")
    # plt.close("all")

    pf.io.write(
        name, compress=True, dtf=dtf.Data_IR, dftf_inverse=dftf_inverse.time)
