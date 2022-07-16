import os
import shutil
import glob
import mesh2hrtf as m2h
from . import utils


def outputs2hrtfs(paths, merge=False, inspect=False, pattern=None,
                  plot=None, plane="horizontal", atol=.1, savedir=None):
    """
    Process NumCalc outputs from multiple projects and write data to disk.

    This is a convenience function that process multiple Mesh2HRTF projects and
    offers to merge and plot the results in one call:

    1.
        Run :py:func:`~mesh2hrtf.output2hrtf` in folders specified by
        `paths`. This also calls
        :py:func:`~mesh2hrtf.write_output_report`,
        :py:func:`~mesh2hrtf.reference_hrtf`, and
        :py:func:`~mesh2hrtf.compute_hrir`.
    2.
        Run :py:func:`~mesh2hrtf.merge_sofa_files` (optional).
    3.
        Run :py:func:`~mesh2hrtf.inspect_sofa_files` (optional). If merge is
        ``True``, only the merged data is inspected.
    4.
        Move data to `savedir` (optional). If merge is ``True`` only the merged
        data will be moved and the remaining data is deleted. The file names
        of the moved data stat with the name of the folder from which they were
        moved to `savedir`.

    Parameters
    ----------
    paths : str, tuple of strings
        Paths to search for Mesh2HRTF projects. If `paths` is a tuple, the
        projects found in the paths will be merged if `merge` is ``True``.
        E.g., if ``paths[0]`` is ``"some/path/left/*"`` and ``paths[1]`` is
        ``"some/path/right/*"`` all SOFA files in the matching folders will be
        merged. If `merge` is ``False`` all projects will be processed
        independently.
    merge : bool, optional
        ``True`` to merge data from different projets using
        :py:func:`~mesh2hrtf.merge_sofa_files`  (see `paths` above). The
        default ``False`` does not merge data.
    inspect : bool, optional
        ``True`` to plot data the results from different projets using
        :py:func:`~mesh2hrtf.inspect_sofa_files`. If `merge` is ``True`` this
        is only done for the merged data. The default ``False`` does not plot
        the data.
    pattern : str, optional
        Merge, inspect and move only SOFA files that have `pattern` in their
        file names. E.g., if `pattern` is ``"HRIR"`` only the SOFA files
        containing HRIR in their name will be merged, inspected, and coped to
        `savedir` (see below).
    plot : str, optional
        This parameter is passed to :py:func:`~mesh2hrtf.inspect_sofa_files`.
    plane : str, optional
        This parameter is passed to :py:func:`~mesh2hrtf.inspect_sofa_files`.
    atol : float, optional
        This parameter is passed to :py:func:`~mesh2hrtf.inspect_sofa_files`.
    savedir : str, optional
        Path under which the results are saved. The default ``None`` will
        leave the data in their project folders, any other value will move
        the data to `savedir` and delete all data in
        `project_folder/Output2HRTF` except for the project reports.

    Raises
    ------
    ValueError
        If issues were detected in any Mesh2HRTF project. The value error is
        raised after all data was processed.
    """

    # check input -------------------------------------------------------------
    if isinstance(paths, str):
        paths = (paths, )
    if not isinstance(paths, (tuple, list)) or len(paths) > 2:
        raise ValueError(
            "paths must be a string or a tuple with a maximum length of two")
    if merge and len(paths) != 2:
        raise ValueError("path must be a tuple of length two if merge is True")

    # create save directories if it does not exist
    # (already done here to throw early errors if savedir is invalid)
    if savedir is not None and not os.path.isdir(savedir):
        os.mkdir(savedir)

    # loop paths for running output_2_hrtf in folders -------------------------
    print("Generating SOFA files")
    print("---------------------")

    # for tracking of issues
    issues = []

    for pp, path in enumerate(paths):

        # loop folders in path
        folders = glob.glob(path)

        for ff, folder in enumerate(folders):

            # check
            name = os.path.basename(folder)
            if not os.path.isfile(os.path.join(folder, "parameters.json")):
                print(f"Skipping: {name}")
                continue

            print((f"{name} (path {pp+1}/{len(paths)}, "
                   f"folder {ff+1}/{len(folders)})"))

            # run output2hrtf
            m2h.output2hrtf(folder)

            # track issues
            if os.path.isfile(os.path.join(
                    folder, "Output2HRTF", "report_issues.txt")):
                issues.append(folder)

    # merge SOFA files --------------------------------------------------------
    if merge:
        print("\nMerging SOFA files")
        print("------------------")
        utils.merge_sofa_files(paths, pattern)

        # if SOFA files were merged, only the merged files are used in the
        # following (merged files saved under paths[0])
        paths_inspect = (paths[0], )
        pattern = "*merged*" if pattern is None else pattern + "*merged*"
    else:
        paths_inspect = paths

    # inspect SOFA files ------------------------------------------------------
    if inspect:

        print("\nInspecting SOFA files")
        print("---------------------")
        for pp, path in enumerate(paths_inspect):
            print(f"{path} ({pp+1}/{len(paths_inspect)})")
            utils.inspect_sofa_files(path, pattern, plot, plane, atol)

    # return if data should be left in place ----------------------------------
    if savedir is None:
        if issues:
            msg = ("Detected issues in NumCalc output. Check report files in "
                   "the following Output2HRTF folders:\n")
            msg += "\n".join(issues)
            raise ValueError(msg)

        print("\nfinished")
        return

    # move data to savedir ----------------------------------------------------
    print("\nMoving results to savedir")
    print("-------------------------")
    for pp, path in enumerate(paths):

        # loop folders in path
        folders = glob.glob(path)
        if pattern is None:
            pattern = "*"

        for ff, folder in enumerate(folders):

            name = os.path.basename(folder)

            # copy data
            if pp == 0 or not merge:
                # sofa files and plots
                for file in glob.glob(os.path.join(
                        folder, "Output2HRTF", pattern)):

                    if file.endswith((".sofa", ".pdf", ".jpeg")):
                        copy_name = name + "_" + os.path.basename(
                            file).replace("_merged", "")
                        shutil.copyfile(file,
                                        os.path.join(savedir, copy_name))

            # remove data
            for file in glob.glob(os.path.join(folder, "Output2HRTF", "*")):
                if "report" not in file:
                    os.remove(file)

    if issues:
        msg = ("Detected issues in NumCalc output. Check report files in "
               f"{savedir} and the following Output2HRTF folders:\n")
        msg += "\n".join(issues)

        with open(os.path.join(savedir, "report_issues.txt"), "w") as issues:
            issues.write(msg)

        raise ValueError(msg)

    print("\nfinished")
    return


def outputs2trash(paths, boundary=False, grid=False,
                  boundary_compressed=False, hrtf=False, vtk=False,
                  reports=False):
    """
    Remove output data from Mesh2HRTF project folder.

    Use this function to delete output that is no longer needed and is taking
    too much disk space.

    Parameters
    ----------
    paths : str, tuple of strings
        Paths under which Mesh2HRTF project folers are searched. Can contain
        `*` remove data from multiple project folders, e.g., `path/*left` will
        remove data from all folders in `path` that end with `left`.
    boundary : bool, optional
        Remove raw pressure and velocity simulated on the boundary, i.e., the
        mesh. This data is saved in
        `project_folder/NumCalc/source_*/be.out/be.*/*Boundary`
    grid : bool, optional
        Remove raw pressure and velocity simulated on the evaluation grid.This
        data is saved in
        `project_folder/NumCalc/source_*/be.out/be.*/*EvalGrid`
    boundary_compressed : bool, optional
        Remove compressed pressure and velocity simulated on the boundary,
        i.e., the mesh. This data is saved in
        `project_folder/Output2HRTF/ObjectMesh_*.npz`.
    hrtf : bool, optional
        Remove HRTFs saved in SOFA files saved in
        `project_folder/Output2HRTF/HRTF_*.sofa`.
    vtk : bool, optional
        Remove vtk exports generated with :py:func:`~mesh2hrtf.output2vtk`
        and saved in `project_folder/Output2HRTF/*_vtk`.
    reports : bool, optional
        Remove project reports saved in `project_folder/Output2HRTF/report_*`.
    """

    # check input
    if isinstance(paths, str):
        paths = (paths, )
    if not isinstance(paths, (tuple, list)):
        raise ValueError(
            "paths must be a string or a tuple of strings")

    # loop paths and contained folders
    for pp, path in enumerate(paths):

        folders = glob.glob(path)

        for ff, folder in enumerate(folders):

            print(f"Purging path {pp+1}/{len(paths)} folder {ff+1}/{folders}")
            print(os.path.basename(folder))

            # check if the directories exist ----------------------------------
            has_numcalc = os.path.isdir(os.path.join(folder, "NumCalc"))
            if has_numcalc:
                numcalc = glob.glob(os.path.join(
                    folder, "NumCalc", "source_*"))

            has_output = os.path.isdir(os.path.join(folder, "Output2HRTF"))
            if has_output:
                output = glob.glob(os.path.join(folder, "Output2HRTF", "*"))

            # data in source*/be.out/be.* folders -----------------------------
            # delete entire be.out folders
            if boundary and grid and has_numcalc:
                for nc in numcalc:
                    shutil.rmtree(os.path.join(nc, "be.out"))
            # delete only the boundary data
            elif boundary and has_numcalc:
                for nc in numcalc:
                    for be in glob.glob(os.path.join(nc, "be.out", "be.*")):
                        os.remove(os.path.join(be, "pBoundary"))
                        os.remove(os.path.join(be, "vBoundary"))
            # delete only the grid data
            elif grid and has_numcalc:
                for nc in numcalc:
                    for be in glob.glob(os.path.join(nc, "be.out", "be.*")):
                        os.remove(os.path.join(be, "pEvalGrid"))
                        os.remove(os.path.join(be, "vEvalGrid"))

            # data in Output2HRTF ---------------------------------------------
            if has_output:
                for oo in output:
                    base = os.path.basename(oo)
                    # remove compressed boundary data
                    if base.endswith(".npz") and boundary_compressed:
                        os.remove(oo)
                    # remove compressed boundary data
                    if base.startswith("HRTF_") and base.endswith(".sofa") \
                            and hrtf:
                        os.remove(oo)
                    # remove compressed boundary data
                    if base.endswith("vtk") and os.path.isdir(oo) \
                            and vtk:
                        shutil.rmtree(oo)
                    # remove compressed boundary data
                    if base.startswith("report_") and reports:
                        os.remove(oo)
