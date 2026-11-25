import numpy as np
import sofar as sf


def reference_hrtfs(sofa, sourceType, sourceArea, speedOfSound, densityOfAir,
                    mode="min"):
    """
    Reference HRTFs to the sound pressure in the center of the head. After
    referencing the sound pressure approaches 1 (0 dB) for low frequencies.

    Parameters
    ----------
    sofa : sofar Sofa object, str
       Sofa object containing the sound pressure or filename of a SOFA file to
       be loaded. SOFA object/file must be of the convention
       SimpleFreeFieldHRTF or GeneralTF
    sourceType : str
        The referencing depends on the source type used for simulating the
        sound pressure. Can be "Both ears", "Left ear", "Right ear",
        "Vibrating elements", "Point source", or "Plane wave"
    sourceArea : array like
        The area of each vibrating source (one value per source). Required
        if vibrating elements are used as source.
    speedOfSound : number
        The speed of sound in m/s
    densityOfAir : number
        The density of air in kg / m**3
    mode : str, optional
        Pass "min", "max", or "mean" to reference to the minmum, maximum, or
        mean radius in `sofa.SourcePosition`. Pass "all" to normalize each
        HRTF to the radius of the corresponding source.

    Returns
    -------
    sofa : sofar Sofa.object
        A copy of the input data with referenced sound pressure
    """

    if isinstance(sofa, str):
        sofa = sf.read_sofa(sofa)
    else:
        sofa = sofa.copy()

    if sofa.GLOBAL_SOFAConventions not in ["SimpleFreeFieldHRTF", "GeneralTF"]:
        raise ValueError(("Sofa object must have the conventions "
                          "SimpleFreeFieldHRTF or GeneralTF"))

    if sofa.SourcePosition_Type == "spherical":
        radius = sofa.SourcePosition[:, 2]
    else:
        radius = np.sqrt(sofa.SourcePosition[:, 0]**2 +
                         sofa.SourcePosition[:, 1]**2 +
                         sofa.SourcePosition[:, 2]**2)

    pressure = sofa.Data_Real + 1j * sofa.Data_Imag
    frequencies = sofa.N

    # distance of source positions from the origin
    if mode == "min":
        r = np.min(radius)
    elif mode == "mean":
        r = np.mean(radius)
    elif mode == "max":
        r = np.max(radius)
    else:
        r = radius[..., np.newaxis, np.newaxis]

    if sourceType in {'Both ears', 'Left ear', 'Right ear',
                      'Vibrating elements'}:

        if len(sourceArea) != pressure.shape[1]:
            raise ValueError(
                (f"sourceArea has {len(sourceArea)} entries but the data "
                 f"contains {pressure.shape[1]} sources."))

        volumeFlow = 0.1 * np.ones(pressure.shape)
        for nn in range(len(sourceArea)):
            volumeFlow[:, nn, :] = \
                volumeFlow[:, nn, :] * sourceArea[nn]

        # point source in the origin evaluated at r
        # eq. (6.71) in: Williams, E. G. (1999). Fourier Acoustics.
        ps = -1j * densityOfAir * 2 * np.pi * frequencies * \
            volumeFlow / (4 * np.pi) * \
            np.exp(1j * 2 * np.pi * frequencies / speedOfSound * r) / r

    elif sourceType == 'Point source':

        amplitude = 0.1  # hard coded in Mesh2HRTF
        ps = amplitude * \
            np.exp(1j * 2 * np.pi * frequencies /
                   speedOfSound * r) / (4 * np.pi * r)

    elif sourceType == 'Plane wave':
        raise ValueError(
            ("Plane wave not implemented yet."))

    else:
        raise ValueError(
            ("Referencing is currently only implemented for "
                "sourceType 'Both ears', 'Left ear', 'Right ear', "
                "'Vibrating elements', 'Point source' and 'Plane wave'."))

    # here we go...
    pressure /= ps
    sofa.Data_Real = np.real(pressure)
    sofa.Data_Imag = np.imag(pressure)

    return sofa