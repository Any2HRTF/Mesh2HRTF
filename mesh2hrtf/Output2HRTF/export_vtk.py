import os
import glob
import json
import numpy as np
from .output2hrtf import _read_nodes_and_elements, _read_numcalc_data


def export_vtk(folder=None, object=None, mode="pressure",
               frequency_steps=None, dB=True, log_prefix=20, log_reference=1,
               deg=False, unwrap=False):
    """
    Export pressure on the (head) mesh to vtk files for importing in ParaView

    The exported vtk files are written to folder/Output2HRTF/vtk
    in a subfolder that contains the name of the `object` and the `mode`.

    Parameters
    ----------
    folder : str, optional
        The Mesh2HRTF project folder. The default ``None`` uses the current
        folder
    object : str, optional
        The name of the mesh or evaluation grid for which the pressure is
        exported. The object is searched in the folders `ObjectMehshes` and
        `EvaluationGrids`. The default ``None`` selects the `Reference` mesh.
    mode : str
        Specify which data should be exported to VTK

        ``'pressure'``
            export the absolute sound pressure
        ``'phase'``
            export the phase
        ``'velocity'``
            export the RMS of the particle velocity

        The default is ``'pressure'``.
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
    deg : bool, optional
        Save the phase in degree. The default ``False`` saves the phase in
        radians.
    unwrap : bool optional
        Unwrap the phase. The default ``False`` does not unwrap the phase.
    """

    # check input data
    if folder is None:
        folder = os.getcwd()
    if not os.path.isfile(os.path.join(folder, "parameters.json")):
        raise ValueError((
            f"The folder {folder} is not a valid Mesh2HRTF project. "
            "It must contain the file 'parameters.json'"))

    if mode not in ["pressure", "phase", "velocity"]:
        raise ValueError((f"mode is '{mode}' but must be 'pressure', "
                          "'phase', or 'velocity'"))

    if object is None:
        object = "Reference"
        object_type = "ObjectMeshes"
    else:
        object_type = None

    # load project parameters
    with open(os.path.join(folder, "parameters.json"), "r") as file:
        params = json.load(file)

    # search object, if required
    if object_type is None:

        for obj in ["ObjectMeshes", "EvaluationGrids"]:
            f = glob.glob(os.path.join(folder, obj, "*"))
            # discard hidden folders that might occur on Mac OS
            f = [os.path.basename(f) for f in f
                 if os.path.isdir(f) and not f.startswith(".")]
            if object in f:
                object_type = obj
                break

    if object_type is None:
        raise ValueError(f"object '{object}' not found in {folder}")

    # get the filename
    file = "Boundary" if object_type == "ObjectMeshes" else "EvalGrid"
    file = "v" + file if mode == "velocity" else "p" + file

    # Load ObjectMesh data
    data, indices = _read_numcalc_data(
        params["numSources"], params["numFrequencies"], folder, file)

    # get the mesh
    grid, _ = _read_nodes_and_elements(
        os.path.join(folder, object_type), object)
    nodes = grid[object]["nodes"]
    elements = grid[object]["elements"]
    num_nodes = grid[object]["num_nodes"]
    num_elems = grid[object]["num_elements"]
    # remove offset in element ids (must start with 0)
    if object_type == "EvaluationGrids":
        elements -= indices[0]
    del grid

    # format data
    if mode == "pressure":
        # convert pressure to dB
        if dB:
            eps = np.finfo(float).eps
            data = log_prefix*np.log10(np.abs(data)/log_reference + eps)

            amp_str = f"{log_prefix}log(pressure/{log_reference})"
            folder_save = "pressure_db"
        # keep pressure linear
        else:
            data = np.abs(data)
            amp_str = "pressure"
            folder_save = "pressure_lin"
    elif mode == "phase":
        data = np.unwrap(np.angle(data))

        deg = "degree" if deg else "radians"

        amp_str = "phase_in_" + deg
        folder_save = "phase_" + deg
        if unwrap:
            amp_str += "_unwrapped"
            folder_save += "_unwrapped"
    else:
        amp_str = "RMS_velocity"
        folder_save = "velocity"

    # get values at element centers if saving data for an evaluation grid
    # this is done by taking the average of each node of an element
    if object_type == "EvaluationGrids":

        # init array of new size
        d = data.copy()
        data = np.zeros(
            (num_elems, params["numSources"], params["numFrequencies"]))

        # average data per element
        for ee in range(num_elems):
            elem_ids = elements[ee, 1:4].astype(int)
            data[ee] = np.mean(d[elem_ids], axis=0)

        del d

    # wrap phase and convert to correct unit
    # (must be done after the averaging above)
    if mode == "phase":
        if not unwrap:
            data = (data + np.pi) % (2 * np.pi) - np.pi
        if deg:
            data *= 180 / np.pi

    # create output folder
    folder_save = os.path.join(
        folder, "Output2HRTF", "vtk", object + "_" + folder_save)
    if not os.path.isdir(folder_save):
        os.makedirs(folder_save)

    # parse frequency steps
    if frequency_steps is None:
        frequency_steps = [1, params["numFrequencies"]]
    if isinstance(frequency_steps, (int, float)) \
            or len(frequency_steps) != 2 \
            or np.any(np.array(frequency_steps) < 1) \
            or any(np.array(frequency_steps) > params["numFrequencies"]):
        raise ValueError(("frequency_steps must contain two values between 1 "
                          f"and {params['numFrequencies']}"))

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
    vtk += f"\nPOLYGONS {elements.shape[0]} {elements.shape[0]*4}\n"
    elements_txt = ""
    for nn in range(elements.shape[0]):
        elements_txt += (f"3 {int(elements[nn, 1])} "
                         f"{int(elements[nn, 2])} "
                         f"{int(elements[nn, 3])}\n")
    vtk += elements_txt

    vtk += f"\nCELL_DATA {data.shape[0]}\n\n"

    # write vtk files
    for ff in range(frequency_steps[0]-1, frequency_steps[1]):

        pressure_txt = ""

        for ss in range(data.shape[1]):
            pressure_txt += f"SCALARS {amp_str}-source_{ss + 1} float 1\n"
            pressure_txt += "LOOKUP_TABLE default\n"

            pp = np.round(data[:, ss, ff], 5)
            for p in pp:
                pressure_txt += str(p) + "\n"

            pressure_txt += "\n"

        vtk_file = os.path.join(folder_save, f"frequency_step_{ff + 1}.vtk")
        with open(vtk_file, "w") as f:
            f.write(vtk + pressure_txt)
