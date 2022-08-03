import os
import shutil
import glob


def remove_outputs(paths, boundary=False, grid=False,
                   hrtf=False, vtk=False,
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
    hrtf : bool, optional
        Remove HRTFs saved in SOFA files saved in
        `project_folder/Output2HRTF/HRTF_*.sofa`.
    vtk : bool, optional
        Remove vtk exports generated with :py:func:`~mesh2hrtf.export_vtk`
        and saved in `project_folder/Output2HRTF/vtk`.
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
