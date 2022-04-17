import os
import glob
import warnings
import numpy as np
from scipy.spatial import Delaunay, ConvexHull
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import sofar as sf
import pyfar as pf
from .output_to_hrtf import _read_nodes_and_elements


def inspect_sofa_files(path, pattern=None, plot=None, plane="horizontal",
                       atol=0.1, savedir=None):
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
        Merge only files that contain `pattern` in their filename. The default
        ``None`` merges all SOFA files.
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
            _inspect_sofa_files(file, save_to, atol, plot, plane)


def merge_sofa_files(paths, pattern=None, savedir=None):
    """
    Merge HRTFs and HRIRs from SOFA-files containing left and right ear data.

    The names of the merged SOFA files and with the "merged.sofa".

    Parameters
    ----------
    paths : tuple
        A tuple containing paths to folders. SOFA files are searched in
        `paths/Output2HRTF` if it exist and directly in `paths` otherwise.

        The names may contain an asterisk to process data in multiple folders.
        E.g., if ``paths[0]`` is ``"some/path/left/*"`` and ``paths[1]`` is
        ``"some/path/right/*"`` all SOFA files in the matching folders will be
        merged. Note the SOFA files contained in the folders must have the same
        names to be merged. Currently, `paths` must contain exactly two
        paths.
    pattern : str
        Merge only files that contain `pattern` in their filename. The default
        ``None`` merges all SOFA files.
    savedir : str
        Directory for saving the merged SOFA files. The default ``None`` saves
        the files to the directory given by `left`.
    """

    if savedir is not None and not os.path.isdir(savedir):
        raise ValueError(f"savedir {savedir} is not a directory")

    if not isinstance(paths, (tuple, list)) or len(paths) != 2:
        raise ValueError("paths, must be a tuple or list of length two")

    # check which data to merge
    if pattern is None:
        pattern = "*.sofa"
    elif not pattern.endswith("sofa"):
        pattern = f"{pattern}*.sofa"

    left = paths[0]
    right = paths[1]

    # get all search directories
    left_dirs = glob.glob(left)
    right_dirs = glob.glob(right)

    if len(left_dirs) != len(right_dirs):
        raise ValueError(("The number of directories found with glob.glob()"
                          f" does not match for {left} and {right}"))

    # loop directories
    for left_dir, right_dir in zip(left_dirs, right_dirs):

        # check if Output2HRTF folder exists
        if os.path.isdir(os.path.join(left_dir, "Output2HRTF")) \
                and os.path.isdir(os.path.join(right_dir, "Output2HRTF")):
            left_dir = os.path.join(left_dir, "Output2HRTF")
            right_dir = os.path.join(right_dir, "Output2HRTF")

        # get and check all SOFA files in Output2HRTF folder
        left_files = glob.glob(os.path.join(left_dir, pattern))
        right_files = glob.glob(os.path.join(right_dir, pattern))

        if len(left_files) != len(right_files):
            raise ValueError((
                "The umber of sofa files found with glob.glob()"
                f" does not match for {left_dir} and {right_dir}"))

        # loop all SOFA files
        for left_file, right_file in zip(left_files, right_files):

            # check file names
            if (os.path.basename(left_file)
                    != os.path.basename(right_file)):
                raise ValueError((
                    "Found mismatching. Each Output2HRTF folder must "
                    "contain SOFA files with the same names. Error for"
                    f" {left_file} and {right_file}"))

            # filename of merged data
            head, tail = os.path.split(left_file)
            tail = tail[:-len(".sofa")] + "_merged.sofa"

            if savedir is not None:
                head = savedir

            # merge data
            _merge_sofa_files(
                (left_file, right_file), os.path.join(head, tail))


def write_evaluation_grid(
        points, name, start=200000, discard=None, show=False):
    """
    Write evaluation grid for use in Mesh2HRTF.

    Mesh2HRTF evaluation grids consist of the two text files Nodes.txt and
    Elements.txt. Evaluations grids are always triangularized.

    Parameters
    ----------
    points : pyfar Coordinates, numpy array
        pyfar Coordinates object or 2D numpy array containing the cartesian
        points of the evaluation grid in meter. The array must be of shape
        (N, 3) with N > 2.
    name : str
        Name of the folder under which the evaluation grid is saved. If the
        folder does not exist, it is created.
    start : int, optional
        The nodes and elements of the evaluation grid are numbered and the
        first element will have the number `start`. In Mesh2HRTF, each Node
        must have a unique number. The nodes/elements of the mesh for which the
        HRTFs are simulated start at 1. Thus `start` must at least be greater
        than the number of nodes/elements in the evaluation grid.
    discard : "x", "y", "z", optional
        In case all values of the evaluation grid are constant for one
        dimension, this dimension has to be discarded during the
        triangularization. E.g. if all points have a z-value of 0 (or any other
        constant), discarded must be "z". The default ``None`` does not discard
        any dimension.
    show : bool, optional
        If ``True`` the evaluation grid is plotted and the plot is saved to
        the folder given by `name`

    Examples
    --------

    Generate a spherical sampling grid with pyfar and write it to the current
    working directory

    .. plot::

        >>> import mesh2hrtf as m2h
        >>> import pyfar as pf
        >>>
        >>> points = pf.samplings.sph_lebedev(sh_order=10)
        >>> m2h.write_evaluation_grid(
        ...     points, "Lebedev_N10", discard=None, show=True)
    """

    if isinstance(points, pf.Coordinates):
        points = points.get_cart()

    if points.ndim != 2 or points.shape[0] < 3 \
            or points.shape[1] != 3:
        raise ValueError(
            "points must be a 2D array of shape (N, 3) with N > 2")

    # discard dimension
    if discard == "x":
        mask = (1, 2)
    elif discard == "y":
        mask = (0, 2)
    elif discard == "z":
        mask = (0, 1)
    else:
        mask = (0, 1, 2)

    # triangulate
    if discard is None:
        tri = ConvexHull(points[:, mask])
    else:
        tri = Delaunay(points[:, mask])

    # check output directory
    if not os.path.isdir(name):
        os.mkdir(name)

    # write nodes
    N = int(points.shape[0])
    start = int(start)

    nodes = f"{N}\n"
    for nn in range(N):
        nodes += (f"{int(start + nn)} "
                  f"{points[nn, 0]} "
                  f"{points[nn, 1]} "
                  f"{points[nn, 2]}\n")

    with open(os.path.join(name, "Nodes.txt"), "w") as f_id:
        f_id.write(nodes)

    # write elements
    N = int(tri.simplices.shape[0])
    elems = f"{N}\n"
    for nn in range(N):
        elems += (f"{int(start + nn)} "
                  f"{tri.simplices[nn, 0] + start} "
                  f"{tri.simplices[nn, 1] + start} "
                  f"{tri.simplices[nn, 2] + start} "
                  "2 0 1\n")

    with open(os.path.join(name, "Elements.txt"), "w") as f_id:
        f_id.write(elems)

    # plot the evaluation grid
    if show:
        points = pf.Coordinates(
            points[:, 0], points[:, 1], points[:, 2]).show()
        plt.savefig(os.path.join(name, "evaluation_grid.png"), dpi=300)


def read_evaluation_grid(name, show=False):
    """
    Read Mesh2HRTF evaluation grid.

    Parameters
    ----------
    name : str
        Name of the folder containing the nodes of the evaluation grid in
        Nodes.txt
    show : bool, optional
        If ``True`` the points of the evaluation grid are plotted. The default
        is ``False``.

    Returns
    -------
    coordinates : pyfar Coordinates
        The points of the evaluation grid as a pyfar Coordinates object
    """

    # check if the grid exists
    if not os.path.isfile(os.path.join(name, "Nodes.txt")):
        raise ValueError(f"{os.path.join(name, 'Nodes.txt')} does not exist")

    # read the nodes
    with open(os.path.join(name, "Nodes.txt"), "r") as f_id:
        nodes = f_id.readlines()

    # get number of nodes
    N = int(nodes[0].strip())
    points = np.zeros((N, 3))

    # get points (first entry is node number)
    for nn in range(N):
        node = nodes[nn+1].strip().split(" ")
        points[nn, 0] = float(node[1])
        points[nn, 1] = float(node[2])
        points[nn, 2] = float(node[3])

    # make coordinates object
    coordinates = pf.Coordinates(points[:, 0], points[:, 1], points[:, 2])

    # plot and return
    if show:
        coordinates.show()

    return coordinates


def export_to_vtk(folder=None, object_mesh=None, frequency_steps=None,
                  dB=True, log_prefix=20, log_reference=1):
    """
    Export pressure on the (head) mesh to vtk files for importing in ParaView

    The exported vtk files are written to folder/Output2HRTF/Reference_vtk
    where 'Reference' is given by `object_mesh`.

    Parameters
    ----------
    folder : str, optional
        The Mesh2HRTF project folder. The default ``None`` uses the current
        folder
    object_mesh : str, optional
        The name of the mesh for which the pressure is exported. The default
        ``None`` selects the `Reference` mesh.
    frequency_steps : list, optional
        List with first and last frequency step for which the data is
        exported. The default ``None`` exports the data for all frequency
        steps.
    dB : bool, optional
        Save pressure in dB if ``True`` (default) and linear otherwise.
    log_prefix, log_reference : number, optional
        The pressure in dB is calculated as
        ``log_prefix * log_10(pressure/log_reference)``. The default values
        are ``20`` and ``1``.
    """

    # check input data
    if folder is None:
        folder = os.getcwd()
    if object_mesh is None:
        object_mesh = "Reference"

    # load object mesh data
    object_name = os.path.join(
            folder, "Output2HRTF", f"ObjectMesh_{object_mesh}.npz")
    if os.path.isfile(object_name):
        # contains: frequencies and pressure
        data = np.load(object_name, allow_pickle=False)
        pressure = data["pressure"]
        frequencies = data["frequencies"]
        del data
    else:
        raise ValueError((f"{object_name} does not exist. "
                          "Run output_to_hrtf to create it"))

    # convert pressure to dB
    if dB:
        eps = np.finfo(float).eps
        pressure = log_prefix*np.log10(np.abs(pressure)/log_reference + eps)

        amp_str = f"{log_prefix}log(pressure/{log_reference})"
        file_str = "db"
    else:
        amp_str = "pressure"
        file_str = "lin"

    # load object mesh
    grid, _ = _read_nodes_and_elements(os.path.join(folder, 'ObjectMeshes'))
    nodes = grid[object_mesh]["nodes"]
    elements = grid[object_mesh]["elements"]
    num_nodes = grid[object_mesh]["num_nodes"]
    del grid

    # create output folder
    savedir = os.path.join(folder, "Output2HRTF", f"{object_mesh}_vtk")
    if not os.path.isdir(savedir):
        os.mkdir(savedir)

    # parse frequency steps
    if frequency_steps is None:
        frequency_steps = [1, frequencies.size]
    if len(frequency_steps) != 2 or np.any(np.array(frequency_steps) < 1) \
            or any(np.array(frequency_steps) > frequencies.size):
        raise ValueError(("frequency_steps must contain two values between 1 "
                          f"and {frequencies.size}"))

    # write constant part of the vtk file
    vtk = ("# vtk DataFile Version 3.0\n"
           "Mesh2HRTF Files\n"
           "ASCII\n"
           "DATASET POLYDATA\n")

    # write nodes
    vtk += f"POINTS {num_nodes} float\n"
    nodes_txt = ""
    for nn in range(num_nodes):
        nodes_txt += (f"{nodes[nn, 1]: .4f} "
                      f"{nodes[nn, 2]: .4f} "
                      f"{nodes[nn, 3]: .4f}\n")
    vtk += nodes_txt

    # write elements
    vtk += f"POLYGONS {elements.shape[0]} {elements.shape[0]*4}\n"
    elements_txt = ""
    for nn in range(elements.shape[0]):
        elements_txt += (f"3 {int(elements[nn, 1])} "
                         f"{int(elements[nn, 2])} "
                         f"{int(elements[nn, 3])}\n")
    vtk += elements_txt

    vtk += "CELL_DATA 2412\n\n"

    # write vtk files
    for ff in range(frequency_steps[0]-1, frequency_steps[1]):

        pressure_txt = ""

        for ss in range(pressure.shape[1]):
            pressure_txt += f"SCALARS {amp_str}-source_{ss + 1} float 1\n"
            pressure_txt += "LOOKUP_TABLE default\n"

            pp = np.round(pressure[:, ss, ff], 5)
            for p in pp:
                pressure_txt += str(p) + "\n"

            pressure_txt += "\n"

        vtk_file = os.path.join(
            savedir, f"{file_str}_frequency_step_{ff + 1}.vtk")
        with open(vtk_file, "w") as f:
            f.write(vtk + pressure_txt)


def _inspect_sofa_files(file, savedir, atol, plot, plane):

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
                        signal[idx], ax=[ax_time[nn], ax_freq[nn]])
                    ax_time[nn].set_title(name)
                else:
                    pf.plot.freq(signal[idx], ax=ax_freq[nn])
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
                        signal[mask, cc], indices=angles,
                        ax=ax_time[cc], cmap="coolwarm")
                    ax_time[cc].set_title(name)

                    # set limits of time plot symmetric
                    c_lim = qm.get_clim()
                    c_lim = np.round(np.max(np.abs(c_lim)), 1)
                    qm.set_clim(-c_lim, c_lim)

                # plot frequency data
                _, qm, _ = pf.plot.freq_2d(signal[mask, cc], ax=ax_freq[cc],
                                           indices=angles, cmap="Reds")
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
            plt.savefig(os.path.join(savedir, tail[:-5] + "_3D.pdf"),
                        bbox_inches="tight")


def _merge_sofa_files(files, savename):
    """read two sofa files, join the data, and save the result"""

    left = sf.read_sofa(files[0])
    right = sf.read_sofa(files[1])

    # join Data
    if left.GLOBAL_DataType.startswith("TF"):

        # check if data can be joined
        if left.N.size != right.N.size or \
                np.any(np.abs(left.N - right.N)) > 1e-6:
            raise ValueError(("Number of frequencies or frequencies do not "
                              f"agree for {left} and {right}"))

        # join data
        left.Data_Real = np.concatenate(
            (left.Data_Real, right.Data_Real), axis=1)
        left.Data_Imag = np.concatenate(
            (left.Data_Imag, right.Data_Imag), axis=1)

    elif left.GLOBAL_DataType.startswith("FIR"):

        # check if data can be joined
        if left.get_dimension("N") != right.get_dimension("N") or \
                left.Data_SamplingRate != right.Data_SamplingRate:
            raise ValueError(("Number of samples or sampling rates do not"
                              f"agree for {left} and {right}"))

        # join data
        left.Data_IR = np.concatenate(
            (left.Data_IR, right.Data_IR), axis=1)
        left.Data_Delay = np.zeros((1, left.Data_IR.shape[1]))
    else:
        raise ValueError("Joining only works for DataTypes 'TF' and 'FIR'")

    left.ReceiverPosition = np.concatenate(
        (np.atleast_2d(left.ReceiverPosition),
         np.atleast_2d(right.ReceiverPosition)), axis=0)

    # write joined data to disk
    sf.write_sofa(savename, left)
