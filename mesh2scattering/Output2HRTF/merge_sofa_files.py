import os
import glob
import numpy as np
import sofar as sf


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
