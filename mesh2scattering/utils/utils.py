import os
import numpy as np
import sofar as sf


def program_root():
    """The root directory of the repository as absolute path. This function
    relies on the correct setting of the environment variable `REPOSITORY_ROOT`
    which is set during the setup of the utils module.

    Returns
    -------
    root : str
        String containing the root directory
    """
    environ = os.path.dirname(os.path.abspath(__file__))
    root = os.path.abspath(os.path.join(environ, os.pardir))
    return root


def repository_root():
    """The root directory of the repository as absolute path. This function
    relies on the correct setting of the environment variable `REPOSITORY_ROOT`
    which is set during the setup of the utils module.

    Returns
    -------
    root : str
        String containing the root directory
    """
    environ = os.path.dirname(os.path.abspath(__file__))
    root = os.path.abspath(os.path.join(environ, os.pardir, os.pardir))
    return root


def _get_sofa_object(data, source_position,
                     receiver_position, Mesh2HRTF_version,
                     frequencies=None, sampling_rate=None):
    """
    Write complex pressure or impulse responses to a SOFA object.

    Parameters
    ----------
    data : numpy array
        The data as an array of shape (MRE)
    evaluation_grid : numpy array
        The evaluation grid in Cartesian coordinates as an array of shape (MC)
    receiver_position : numpy array
        The position of the receivers (ears) in Cartesian coordinates
    mode : str
        "HRTF" to save HRTFs, "HRIR" to save HRIRs
    Mesh2HRTF_version : str
    frequencies : numpy array
        The frequencies at which the HRTFs were calculated. Required if mode is
        "HRTF"
    sampling_rate :
        The sampling rate. Required if mode is "HRIR"

    Returns
    -------
    sofa : sofar.Sofa object
    """

    # create empty SOFA object
    if True:
        convention = "GeneralTF"
    else:
        convention = "GeneralFIR"

    sofa = sf.Sofa(convention)

    # write meta data
    sofa.GLOBAL_ApplicationName = 'Mesh2scattering'
    sofa.GLOBAL_ApplicationVersion = Mesh2HRTF_version
    sofa.GLOBAL_History = "numerically simulated data"

    # Source and receiver data
    sofa.SourcePosition = source_position
    sofa.SourcePosition_Units = "meter"
    sofa.SourcePosition_Type = "cartesian"

    sofa.ReceiverPosition = receiver_position
    sofa.ReceiverPosition_Units = "meter"
    sofa.ReceiverPosition_Type = "cartesian"

    # HRTF/HRIR data
    if data.shape[0] != source_position.shape[0]:
        data = np.swapaxes(data, 0, 1)
    if True:
        sofa.Data_Real = np.real(data)
        sofa.Data_Imag = np.imag(data)
        sofa.N = np.array(frequencies).flatten()
    else:
        sofa.Data_IR = data
        sofa.Data_SamplingRate = sampling_rate
        sofa.Data_Delay = np.zeros((1, data.shape[1]))

    return sofa
