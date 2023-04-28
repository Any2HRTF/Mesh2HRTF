import os
import glob
import warnings
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import pyfar as pf


def inspect_sofa_files(path, pattern=None, plot=None, plane="horizontal",
                       atol=0.1, savedir=None, dB_time=False, dB_freq=True,
                       freq_scale='log'):
    """
    Inspect SOFA files through plots.

    Generate and save plots for horizontal plane HRIRs (time domain) and HRTFs
    (frequency domain) for one or multiple SOFA files.

    Parameters
    ----------
    path : str
        Path to a folder containing Mesh2HRTF projects. SOFA files are searched
        in `path/Output2HRTF` if it exist and directly in `path` otherwise.

        The name may contain an asterisk to process data in multiple folders.
        E.g., if `path` is ``"some/path/HRIRs_*"`` files in all folder
        starting with "HRIRs" would be scanned for SOFA files.
    pattern : str
        Plot only files that contain `pattern` in their filename. The default
        ``None`` plots all SOFA files.
    plot : str, optional
        ``"2D"``
            generate line plots of four sources on the horizontal plane
            (front, back, left, right). The closest sources to the ideal
            positoins are used for plotting.
        ``"3D"``
            generate color coded plots of all sources on the horizontal
            plane. See also parameter `atol` below.

        The default ``None`` generate both plots.
    plane : str, optional
        Select the plane for which is shown in the 3D plot. Can be
        ``"horizontal"`` (default), ``"median"``, or ``"frontal"``
    atol : float, optional
        Sources on the `plane` are searched within a range +/- `atol` degree.
        The default is ``0.1``.
    savedir : str
        Directory for saving the merged SOFA files. The default ``None`` saves
        the files to the directory given by `path`.
    dB_time : bool, optional
        Plot the logarithmic time data. The default is ``'False'``.
    dB_freq : bool, optional
        Plot the logarithmic magnitude data. The default is ``'True'``.
    freq_scale : str, optional
        ``'log'`` to plot on a logarithmic frequency axis and ``'linear'`` to
        plot on a linear frequency axis. The default is ``'log'``.
    """

    # check input
    if plot is None:
        plot = ["2D", "3D"]
    elif plot not in ["2D", "3D"]:
        raise ValueError(
            f"plot is {plot} but must be 2D, 3D, or all")
    if plane not in ["horizontal", "median", "frontal"]:
        raise ValueError(
            f"plane is {plane} but must be horizontal, median, or frontal")

    # check which data to merge
    if pattern is None:
        pattern = "*.sofa"
    elif not pattern.endswith("sofa"):
        pattern = f"{pattern}*.sofa"

    # get all directories containing SOFA files
    folders = glob.glob(path)

    # loop directories
    for folder in folders:

        # check if Output2HRTF folder exists
        if os.path.isdir(os.path.join(folder, "Output2HRTF")):
            folder = os.path.join(folder, "Output2HRTF")

        # find matching SOFA files
        files = glob.glob(os.path.join(folder, pattern))
        if not files:
            raise ValueError((f"Did not find any SOFA files in {folder} "
                              f"that are matching {pattern}"))

        # loop and inspect all SOFA files
        for file in files:

            # inspect data
            save_to = folder if savedir is None else None
            _inspect_sofa_files(file, save_to, atol, plot, plane,
                                dB_time, dB_freq, freq_scale)


def _inspect_sofa_files(file, savedir, atol, plot, plane,
                        dB_time, dB_freq, freq_scale):

    with Dataset(file, "r", format="NETCDF4") as sofa_file:
        data_type = getattr(sofa_file, "DataType")

    if data_type == "FIR":
        mode = "hrir"
    elif data_type == "TF":
        mode = "hrtf"
    else:
        raise ValueError(
            f"The DataType of {file} is {data_type} but must 'FIR' or 'TF'")

    # load sofa file and source positions
    signal, sources, _ = pf.io.read_sofa(file)
    tail = os.path.basename(file)

    if "2D" in plot:
        with pf.plot.context():

            # generate plot layout and axes
            if mode == "hrir":
                _, ax = plt.subplots(2, 4, figsize=(16, 6), sharey="row")
                ax_time = ax[0]
                ax_freq = ax[1]
            else:
                _, ax_freq = plt.subplots(1, 4, figsize=(16, 3), sharey="row")

            # loop sources
            max_db = -300
            for nn, (az, name) in enumerate(zip(
                    [0, 180, 90, 270],
                    ["front", "back", "left", "right"])):

                # find current source
                idx, _ = sources.find_nearest_k(
                    az, 0, 1, 1, 'sph', 'top_elev', 'deg')

                # exact position for plotting
                source = sources.get_sph('top_elev', 'deg')[idx]
                name += (f" (az. {np.round(source[0])}, "
                         f"el. {np.round(source[1])} deg.)")

                # plot
                if mode == "hrir":
                    pf.plot.time_freq(
                        signal[idx], ax=[ax_time[nn], ax_freq[nn]],
                        dB_time=dB_time, dB_freq=dB_freq,
                        freq_scale=freq_scale)
                    ax_time[nn].set_title(name)
                else:
                    pf.plot.freq(signal[idx], ax=ax_freq[nn], dB=dB_freq,
                                 freq_scale=freq_scale)
                    ax_freq[nn].set_title(name)

                max_db = np.max([ax_freq[nn].get_ylim()[1], max_db])

            # format axis and legend
            ax_freq[nn].set_ylim(max_db - 60, max_db)

            if signal.cshape[-1] == 2:
                ax_freq[nn].legend(["left ear", "right ear"], loc=3)
            elif signal.cshape[-1] > 2:
                ax_freq[nn].legend(
                    [f"ch. {cc+1}" for cc in range(signal.cshape[-1])], loc=3)

            # save
            plt.tight_layout()
            plt.savefig(os.path.join(savedir, tail[:-5] + "_2D.pdf"),
                        bbox_inches="tight")

    if "3D" in plot:

        num_chanel = signal.cshape[-1]

        with pf.plot.context():

            # generate plot layout and axes
            if mode == "hrir":
                _, ax = plt.subplots(2, num_chanel, sharex=True, sharey="row",
                                     figsize=(4*num_chanel, 6),
                                     )
                ax_time = np.atleast_1d(ax[0])
                ax_freq = np.atleast_1d(ax[1])
            else:
                _, ax_freq = plt.subplots(
                    1, num_chanel, figsize=(4*num_chanel, 3),
                    sharex=True, sharey="row")
                ax_freq = np.atleast_1d(ax_freq)

            # find sources on the desired plane
            if plane == "horizontal":
                _, mask = sources.find_slice("elevation", "deg", 0, atol)
                angles = sources.get_sph('top_elev', 'deg')[mask, 0]
                angle = "azimuth"
            elif plane == "median":
                _, mask = sources.find_slice("lateral", "deg", 0, atol)
                angles = sources.get_sph('side', 'deg')[mask, 1]
                angle = "polar angle"
            else:
                _, mask = sources.find_slice("theta", "deg", 90, atol)
                angles = sources.get_sph('front', 'deg')[mask, 0]
                angle = "theta"

            if not np.any(mask):
                warnings.warn((
                    "Did not find and sources on the horizontal plane for "
                    f"within +/-{atol} deg. for {file}"))
                return

            # plot titles
            names = ["left ear", "right ear"] if signal.cshape[-1] == 2 \
                else [f"ch. {cc+1}" for cc in range(signal.cshape[-1])]

            # loop sources
            for cc, name in enumerate(names):

                # plot time data
                if mode == "hrir":
                    _, qm, _ = pf.plot.time_2d(
                        signal[mask, cc], indices=angles, dB=dB_time,
                        ax=ax_time[cc], cmap="coolwarm")
                    ax_time[cc].set_title(name)

                    # set limits of time plot (symmetric for nice colors, clip
                    # to improve visibility)
                    c_lim = qm.get_clim()
                    c_lim = np.round(.6 * np.max(np.abs(c_lim)), 1)
                    qm.set_clim(-c_lim, c_lim)

                # plot frequency data
                _, qm, _ = pf.plot.freq_2d(
                    signal[mask, cc], ax=ax_freq[cc], indices=angles,
                    dB=dB_freq, freq_scale=freq_scale, cmap="Reds")
                if mode == "hrir":
                    ax_freq[cc].set_xlabel(f"{angle} in degree")
                    ax_time[cc].set_xlabel("")
                else:
                    ax_freq[cc].set_title(name)
                    ax_freq[cc].set_xlabel(f"{angle} in degree")

                c_lim = qm.get_clim()
                c_lim = np.round(np.max(c_lim))
                qm.set_clim(c_lim - 60, c_lim)

            # save
            plt.tight_layout()
            plt.savefig(os.path.join(savedir, tail[:-5] + "_3D.jpeg"),
                        dpi=300, bbox_inches="tight")
