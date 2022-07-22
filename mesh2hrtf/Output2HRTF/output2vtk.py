import os
import glob
import numpy as np
from .output2hrtf import _read_nodes_and_elements


def output2vtk(folder=None, object=None, mode="magnitude",
               frequency_steps=None, source=None,
               dB=True, log_prefix=20, log_reference=1, unwrap_phase=False):
    """
    Export pressure on the (head) mesh to vtk files for importing in ParaView

    The exported vtk files are written to folder/Output2HRTF/Reference_vtk
    where 'Reference' is given by `object`.

    Parameters
    ----------
    folder : str, optional
        The Mesh2HRTF project folder. The default ``None`` uses the current
        folder
    object : str, optional
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

    if object is None:
        object = "Reference"
        object_type = "ObjectMeshes"
    else:
        object_type = None

    # search object, if required
    if object_type is None:

        for obj in ["ObjectMeshes", "EvaluationGrids"]:
            f = object in glob.glob(os.path.join(folder, obj, "*"))
            f = [os.path.basename(f) for f in f if os.path.isdir(f)]
            if object in f:
                object_type = obj
                break

    if object_type is None:
        raise ValueError(f"object '{object}' not found in {folder}")

    # look up

    # load object mesh data
    object_name = os.path.join(
            folder, "Output2HRTF", f"ObjectMesh_{object}.npz")
    if os.path.isfile(object_name):
        # contains: frequencies and pressure
        data = np.load(object_name, allow_pickle=False)
        pressure = data["pressure"]
        frequencies = data["frequencies"]
        del data
    else:
        raise ValueError((f"{object_name} does not exist. "
                          "Run output2hrtf to create it"))

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
    nodes = grid[object]["nodes"]
    elements = grid[object]["elements"]
    num_nodes = grid[object]["num_nodes"]
    del grid

    # create output folder
    savedir = os.path.join(folder, "Output2HRTF", f"{object}_vtk")
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

    vtk += f"\nCELL_DATA {pressure.shape[0]}\n\n"

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
