import numpy as np


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
