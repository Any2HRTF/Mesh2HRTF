import os
import numpy as np
import pyfar as pf
from scipy.spatial import Delaunay, ConvexHull
import trimesh
import json
import datetime
from mesh2scattering import utils
import shutil
import mesh2scattering as m2s
from packaging import version


def create_source_positions(phi_deg, theta_deg, radius):
    theta_rad = theta_deg * np.pi / 180.
    phi_rad = phi_deg * np.pi / 180.
    theta, phi = np.meshgrid(theta_rad, phi_rad)
    theta = theta.flatten()
    phi = phi.flatten()
    # create coordinates
    x = radius * np.sin(theta) * np.cos(phi)
    y = radius * np.sin(theta) * np.sin(phi)
    z = radius * np.cos(theta)
    return pf.Coordinates(x, y, z)


def write_scattering_project(
        project_path, frequencies, sample_path, reference_path,
        receiver_coords, source_coords,
        structural_wavelength=0, model_scale=1, symmetry_azimuth=[90, 180],
        symmetry_rotational=False, sample_diameter=0.8,
        speed_of_sound='346.18',
        density_of_medium='1.1839'):

    if not os.path.isdir(project_path):
        os.mkdir(project_path)

    frequencyStepSize = 0
    title = 'scattering coefficient Sample'
    method = 'ML-FMM BEM'
    project_path_sample = os.path.join(project_path, 'sample')
    write_project(
        project_path_sample, title, frequencies, frequencyStepSize,
        sample_path,
        receiver_coords, source_coords, sourceType='Point source',
        method=method, materialSearchPaths=None,
        speedOfSound=speed_of_sound,
        densityOfMedium=density_of_medium, materials=None)

    title = 'scattering coefficient Reference'
    sourcePositions_ref = source_coords[
        np.abs(source_coords.get_sph()[..., 0]) < 1e-14]
    project_path_ref = os.path.join(project_path, 'reference')

    write_project(
        project_path_ref, title, frequencies, frequencyStepSize,
        reference_path, receiver_coords, sourcePositions_ref,
        sourceType='Point source',
        method=method,  materialSearchPaths=None,
        speedOfSound=speed_of_sound,
        densityOfMedium=density_of_medium, materials=None)

    source_list = [list(i) for i in list(source_coords.get_cart())]
    receiver_list = [list(i) for i in list(receiver_coords.get_cart())]
    title = 'scattering pattern'
    frequencies = np.array(frequencies, dtype=float)
    parameters = {
        # project Info
        "project_title": 'scattering pattern',
        "mesh2scattering_path": utils.program_root(),
        "mesh2scattering_version": m2s.__version__,
        "bem_version": 'ML-FMM BEM',
        # Constants
        "speed_of_sound": float(346.18),
        "density_of_medium": float(1.1839),
        # Sample Information, post processing
        "structural_wavelength": structural_wavelength,
        "model_scale": model_scale,
        "sample_diameter": sample_diameter,
        # symmetry information
        "symmetry_azimuth": symmetry_azimuth,
        "symmetry_rotational": symmetry_rotational,
        # frequencies
        "num_frequencies": len(frequencies),
        "min_frequency": frequencies[0],
        "max_frequency": frequencies[-1],
        "frequencies": list(frequencies),
        # Source definition
        "source_type": 'Point source',
        "sources_num": len(source_list),
        "sources": source_list,
        # Receiver definition
        "receivers_num": len(receiver_list),
        "receivers": receiver_list,

    }
    with open(os.path.join(project_path, "parameters.json"), 'w') as file:
        json.dump(parameters, file, indent=4)


def write_project(
        project_path, title, frequencies, frequencyStepSize, mesh_path,
        evaluationPoints, sourcePositions,
        sourceType='Point source', method='ML-FMM BEM',
        materialSearchPaths=None, speedOfSound='346.18',
        densityOfMedium='1.1839', materials=None):

    programPath = utils.program_root()
    defaultPath = os.path.join(
        programPath, 'Mesh2Input', 'Materials', 'Data')
    if materialSearchPaths is None:
        materialSearchPaths = defaultPath
    else:
        materialSearchPaths += f";  {defaultPath}"

    # create folders
    if not os.path.isdir(project_path):
        os.mkdir(project_path)
    if not os.path.isdir(os.path.join(project_path, 'ObjectMeshes')):
        os.mkdir(os.path.join(project_path, 'ObjectMeshes'))
    if not os.path.isdir(os.path.join(project_path, 'NumCalc')):
        os.mkdir(os.path.join(project_path, 'NumCalc'))
    if not os.path.isdir(os.path.join(project_path, 'EvaluationGrids')):
        os.mkdir(os.path.join(project_path, 'EvaluationGrids'))

    # write stl file
    mesh = trimesh.load(mesh_path)
    path = os.path.join(project_path, 'ObjectMeshes', 'Reference')
    write_mesh(mesh.vertices, mesh.faces, path, start=0)

    # copy stl file
    mesh = trimesh.load(mesh_path)
    path = shutil.copyfile(mesh_path, os.path.join(
        project_path, 'ObjectMeshes', mesh_path.split(os.sep)[-1]))

    # write evaluation grid
    write_evaluation_grid(
        evaluationPoints,
        os.path.join(project_path, 'EvaluationGrids', 'grid'))

    # Write NumCalc input files for all sources (NC.inp) ----------------------
    _write_nc_inp(
        project_path, version.parse(m2s.__version__), title, speedOfSound,
        densityOfMedium, frequencies, ['grid'], materials, method, sourceType,
        sourcePositions, len(mesh.faces), len(mesh.vertices))


def write_mesh(vertices, faces, path, start=200000, discard=None):
    if vertices.ndim != 2 or vertices.shape[0] < 3 \
            or vertices.shape[1] != 3:
        raise ValueError(
            "vertices must be a 2D array of shape (N, 3) with N > 2")

    # check output directory
    if not os.path.isdir(path):
        os.mkdir(path)

    # write nodes
    N = int(vertices.shape[0])
    start = int(start)

    nodes = f"{N}\n"
    for nn in range(N):
        nodes += (f"{int(start + nn)} "
                  f"{vertices[nn, 0]} "
                  f"{vertices[nn, 1]} "
                  f"{vertices[nn, 2]}\n")

    with open(os.path.join(path, "Nodes.txt"), "w") as f_id:
        f_id.write(nodes)

    # write elements
    N = int(faces.shape[0])
    elements = f"{N}\n"
    for nn in range(N):
        elements += (
            f"{int(start + nn)} "
            f"{faces[nn, 0] + start} "
            f"{faces[nn, 1] + start} "
            f"{faces[nn, 2] + start} "
            "0 0 0\n")

    with open(os.path.join(path, "Elements.txt"), "w") as f_id:
        f_id.write(elements)


def _write_nc_inp(
        filepath1, version, title, speedOfSound, densityOfMedium, frequencies,
        evaluationGrids, materials, method, sourceType, sourcePositions,
        numElementsMesh, numNodesMesh):
    """Write NC.inp file that is read by NumCalc to start the simulation.

    The file format is documented at:
    https://github.com/Any2HRTF/Mesh2HRTF/wiki/Structure_of_NC.inp
    """
    if isinstance(sourcePositions, pf.Coordinates):
        sourcePositions = sourcePositions.get_cart()

    # check the BEM method
    if method == 'BEM':
        method_id = 0
    elif method == 'SL-FMM BEM':
        method_id = 1
    elif method == 'ML-FMM BEM':
        method_id = 4
    else:
        ValueError(
            f"Method must be BEM, SL-FMM BEM or ML-FMM BEM but is {method}")

    for source in range(sourcePositions.shape[0]):

        # create directory
        filepath2 = os.path.join(
            filepath1, "NumCalc", f"source_{source+1}")
        if not os.path.exists(filepath2):
            os.mkdir(filepath2)

        # write NC.inp
        file = open(os.path.join(filepath2, "NC.inp"), "w",
                    encoding="utf8", newline="\n")
        fw = file.write

        # header --------------------------------------------------------------
        fw("##-------------------------------------------\n")
        fw("## This file was created by mesh2input\n")
        fw("## Date: %s\n" % datetime.date.today())
        fw("##-------------------------------------------\n")
        fw("Mesh2HRTF %s\n" % version)
        fw("##\n")
        fw("%s\n" % title)
        fw("##\n")

        # control parameter I (hard coded, not documented) --------------------
        fw("## Controlparameter I\n")
        fw("0 0 0 0 7 0\n")
        fw("##\n")

        # control parameter II ------------------------------------------------
        fw("## Controlparameter II\n")
        fw("1 %d 0.000001 0.00e+00 1 0 0\n" % (
            len(frequencies)))
        fw("##\n")
        fw("## Load Frequency Curve \n")
        fw("0 %d\n" % (len(frequencies)+1))
        fw("0.000000 0.000000e+00 0.0\n")
        for ii in range(len(frequencies)):
            fw("%f %fe+04 0.0\n" % (
                0.000001*(ii+1),
                frequencies[ii] / 10000))
        fw("##\n")

        # main parameters I ---------------------------------------------------
        fw("## 1. Main Parameters I\n")
        numNodes = 0
        numElements = 0
        for evaluationGrid in evaluationGrids:
            # read number of nodes
            nodes = open(os.path.join(
                filepath1, "EvaluationGrids", evaluationGrid,
                "Nodes.txt"))
            line = nodes.readline()
            numNodes = numNodes+int(line)
            # read number of elements
            elements = open(os.path.join(
                filepath1, "EvaluationGrids", evaluationGrid,
                "Elements.txt"))
            line = elements.readline()
            numElements = numElements+int(line)
        fw("2 %d " % (numElementsMesh+numElements))
        fw("%d 0 " % (numNodesMesh+numNodes))
        fw("0")
        fw(" 2 1 %s 0\n" % (method_id))
        fw("##\n")

        # main parameters II --------------------------------------------------
        fw("## 2. Main Parameters II\n")
        if "plane" in sourceType:
            fw("1 ")
        else:
            fw("0 ")
        if "ear" in sourceType:
            fw("0 ")
        else:
            fw("1 ")
        fw("0 0.0000e+00 0 0 0\n")
        fw("##\n")

        # main parameters III -------------------------------------------------
        fw("## 3. Main Parameters III\n")
        fw("0 0 0 0\n")
        fw("##\n")

        # main parameters IV --------------------------------------------------
        fw("## 4. Main Parameters IV\n")
        fw("%s %se+00 1.0 0.0e+00 0.0 e+00 0.0e+00 0.0e+00\n" % (
            speedOfSound, densityOfMedium))
        fw("##\n")

        # nodes ---------------------------------------------------------------
        fw("NODES\n")
        fw("../../ObjectMeshes/Reference/Nodes.txt\n")
        # write file path of nodes to input file
        for grid in evaluationGrids:
            fw("../../EvaluationGrids/%s/Nodes.txt\n" % grid)
        fw("##\n")
        fw("ELEMENTS\n")
        fw("../../ObjectMeshes/Reference/Elements.txt\n")
        # write file path of elements to input file
        for grid in evaluationGrids:
            fw("../../EvaluationGrids/%s/Elements.txt\n" % grid)
        fw("##\n")

        # SYMMETRY ------------------------------------------------------------
        fw("# SYMMETRY\n")
        fw("# 0 0 0\n")
        fw("# 0.0000e+00 0.0000e+00 0.0000e+00\n")
        fw("##\n")

        # assign mesh elements to boundary conditions -------------------------
        # (including both, left, right ear)
        fw("BOUNDARY\n")
        # write velocity condition for the ears if using vibrating
        # elements as the sound source
        if "ear" in sourceType:
            if source == 0 and \
                    sourceType in ['Both ears', 'Left ear']:
                tmpEar = 'Left ear'
            else:
                tmpEar = 'Right ear'
            fw(f"# {tmpEar} velocity source\n")
            fw("ELEM %i TO %i VELO 0.1 -1 0.0 -1\n" % (
                materials[tmpEar]["index_start"],
                materials[tmpEar]["index_end"]))
        # remaining conditions defined by frequency curves
        curves = 0
        steps = 0
        if materials is not None:
            for m in materials:
                if materials[m]["path"] is None:
                    continue
                # write information
                fw(f"# Material: {m}\n")
                fw("ELEM %i TO %i %s 1.0 %i 1.0 %i\n" % (
                    materials[m]["index_start"],
                    materials[m]["index_end"],
                    materials[m]["boundary"],
                    curves + 1, curves + 2))
                # update metadata
                steps = max(steps, len(materials[m]["freqs"]))
                curves += 2

        fw("RETU\n")
        fw("##\n")

        # source information: point source and plane wave ---------------------
        if sourceType == "Point source":
            fw("POINT SOURCES\n")
        elif sourceType == "Plane wave":
            fw("PLANE WAVES\n")
        if sourceType in ["Point source", "Plane wave"]:
            fw("0 %s %s %s 0.1 -1 0.0 -1\n" % (
                sourcePositions[source, 0], sourcePositions[source, 1],
                sourcePositions[source, 2]))
        fw("##\n")

        # curves defining boundary conditions of the mesh ---------------------
        if curves > 0:
            fw("CURVES\n")
            # number of curves and maximum number of steps
            fw(f"{curves} {steps}\n")
            curves = 0
            for m in materials:
                if materials[m]["path"] is None:
                    continue
                # write curve for real values
                curves += 1
                fw(f"{curves} {len(materials[m]['freqs'])}\n")
                for f, v in zip(materials[m]['freqs'],
                                materials[m]['real']):
                    fw(f"{f} {v} 0.0\n")
                # write curve for imaginary values
                curves += 1
                fw(f"{curves} {len(materials[m]['freqs'])}\n")
                for f, v in zip(materials[m]['freqs'],
                                materials[m]['imag']):
                    fw(f"{f} {v} 0.0\n")

        else:
            fw("# CURVES\n")
        fw("##\n")

        # post process --------------------------------------------------------
        fw("POST PROCESS\n")
        fw("##\n")
        fw("END\n")
        file.close()


def read_material_data(materials):

    for material in materials:
        # current material file
        file = materials[material]["path"]
        # check if the file exists
        if file is None:
            continue

        # initialize data
        boundary = None
        freqs = []
        real = []
        imag = []

        # read the csv material file
        with open(file, 'r') as m:
            lines = m.readlines()

        # parse the file
        for line in lines:
            line = line.strip('\n')
            # skip empty lines and comments
            if not len(line):
                continue
            if line[0] == '#':
                continue

            # detect boundary keyword
            if line in ['ADMI', 'IMPE', 'VELO', 'PRES']:
                boundary = line
            # read curve value
            else:
                line = line.split(',')
                if not len(line) == 3:
                    raise ValueError(
                        (f'Expected three values in {file} '
                         f'definition but found {len(line)}'))
                freqs.append(line[0].strip())
                real.append(line[1].strip())
                imag.append(line[2].strip())

        # check if boundary keyword was found
        if boundary is None:
            raise ValueError(
                (f"No boundary definition found in {file}. "
                 "Must be 'ADMI', 'IMPE', 'VELO', or 'PRES'"))
        # check if frequency vector is value
        for i in range(len(freqs)-1):
            if float(freqs[i+1]) <= float(freqs[i]):
                raise ValueError((f'Frequencies in {file} '
                                  'do not increase monotonously'))

        # create output
        materials[material]['boundary'] = boundary
        materials[material]['freqs'] = freqs
        materials[material]['real'] = real
        materials[material]['imag'] = imag

    return materials


def write_material(filename, kind, frequencies, data, comment=None):
    """
    Write boundary condition to file.

    Mesh2HRTF supports non-rigid boundary conditions in the form of text files.
    Such files can be written with this function.

    Parameters
    ----------
    filename : str
        Name of the material file that is written to disk. Must end with ".csv"
    kind : str
        Defines the kind of boundary condition

        ``"pressure"``
            A pressure boundary condition can be used to force a certain
            pressure on the boundary of the mesh. E.g., a pressure of 0 would
            define a sound soft boundary.
        ``"velocity"``
            A velocity boundary condition can be used to force a certain
            velocity on the boundary of the mesh. E.g., a velocity of 0 would
            define a sound hard boundary.
        ``admittance``
            A normalized admittance boundary condition can be used to define
            arbitrary boundaries. The admittance must be normalized, i.e.,
            admittance/(rho*c) has to be provided, which rho being the density
            of air in kg/m**3 and c the speed of sound in m/s.
    frequencies : array like
        The frequencies for which the boundary condition is given
    data : array like
        The values of the boundary condition at the frequencies given above.
    comment : str, optional
        A comment that is written to the beginning of the material file. The
        default ``None`` does omit the comment.

    Notes
    -----
    Mesh2HRTF performs an interpolation in case the boundary condition is
    required at frequencies that are not specified. The interpolation is linear
    between the lowest and highest provided frequency and uses the nearest
    neighbor outside this range.
    """

    # check input
    if not filename.endswith(".csv"):
        raise ValueError("The filename must end with .csv")

    if len(frequencies) != len(data):
        raise ValueError("frequencies and data must have the same lengths")

    # write the comment
    file = ""
    if comment is not None:
        file += "# " + comment + "\n#\n"

    # write the kind of boundary condition
    file += ("# Keyword to define the boundary condition:\n"
             "# ADMI: Normalized admittance boundary condition\n"
             "# PRES: Pressure boundary condition\n"
             "# VELO: Velocity boundary condition\n"
             "# NOTE: Mesh2HRTF expects normalized admittances, i.e., "
             "admittance/(rho*c).\n"
             "#       rho is the density of air and c the speed of sound. "
             "The normalization is\n"
             "#       beneficial because a single material file can be used "
             "for simulations\n"
             "#       with differing speed of sound and density of air.\n")

    if kind == "admittance":
        file += "ADMI\n"
    elif kind == "pressure":
        file += "PRES\n"
    elif kind == "velocity":
        file += "VELO\n"
    else:
        raise ValueError("kind must be admittance, pressure, or velocity")

    file += ("#\n"
             "# Frequency curve:\n"
             "# Mesh2HRTF performs an interpolation in case the boundary "
             "condition is required\n"
             "# at frequencies that are not specified. The interpolation is "
             "linear between the\n"
             "# lowest and highest provided frequency and uses the nearest "
             "neighbor outside\n"
             "# this range.\n"
             "#\n"
             "# Frequency in Hz, real value, imaginary value\n")

    # write data
    for f, d in zip(frequencies, data):
        file += f"{f}, {np.real(d)}, {np.imag(d)}\n"

    # write to file
    with open(filename, "w") as f_id:
        f_id.write(file)


def write_evaluation_grid(
        points, name, start=200000, discard=None):
    """
    Write evaluation grid for use in Mesh2HRTF.

    Mesh2HRTF evaluation grids consist of the two text files Nodes.txt and
    Elements.txt. Evaluations grids are always triangulated.

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

    Examples
    --------

    Generate a spherical sampling grid with pyfar and write it to the current
    working directory

    .. plot::

        >>> import mesh2scattering as m2s
        >>> import pyfar as pf
        >>>
        >>> points = pf.samplings.sph_lebedev(sh_order=10)
        >>> m2s.input.write_evaluation_grid(
        ...     points, "Lebedev_N10", discard=None)
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


def read_evaluation_grid(name):
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

    return coordinates
