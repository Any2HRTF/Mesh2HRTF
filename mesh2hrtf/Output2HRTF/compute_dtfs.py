import numpy as np
import sofar as sf
import pyfar as pf
from pyfar.io.io import _sofa_pos


def compute_dtfs(sofa, smooth_fractions=None, phase="minimum"):
    r"""
    Compute directional transfer functions from HRIRs.

    The directional transfer functions (DTFs) are computed by a spectral
    division of the HRTFs by the diffuse field transfer function (DFTF)

    .. math::

        \mathrm{DTF} = \frac{\mathrm{HRTF}}{\mathrm{DTFT}}

    where the CTF is computed as the energetic average across source positions
    and ears

    .. math::

        \mathrm{DFTF} = \sqrt{\sum_{n=0}^{N-1} |\mathrm{HRTF}|^2_{l,n}\,w_n +
        \sum_{n=0}^{N-1} |\mathrm{HRTF}|^2_{r,n}\,w_n}

    where the index :math:`n` denotes the source position, :math:`|\cdot|` the
    absolute spectrum, and :math:`\sum_n w_n` normalized area weights for
    numerical integration. `Spherical Voronoi weights
    <https://pyfar.readthedocs.io/en/latest/modules/pyfar.samplings.html\
    ?highlight=voronoi#pyfar.samplings.calculate_sph_voronoi_weights>`_ are
    used for this purpose.

    Parameters
    ----------
    sofa : sofar Sofa object, str
       Sofa object containing the sound pressure or filename of a SOFA file to
       be loaded. SOFA object/file must be of the convention
       SimpleFreeFieldHRTF or GeneralTF
    smooth_fractions : number, None, optional
        Apply `fractional octave smoothing
        <https://pyfar.readthedocs.io/en/latest/modules/pyfar.dsp.html\
        ?highlight=smooth#pyfar.dsp.smooth_fractional_octave>`_ to theDFTF.
        E.g. a value of ``3`` applies third octave smoothing and a value of
        ``1`` applies octave smoothing. The default ``None`` does not apply any
        smoothing.
    phase : string, optional
        Define the phase of the DTFT. ``'minimum'`` generates a `minimum phase
        <https://pyfar.readthedocs.io/en/latest/modules/pyfar.dsp.html?\
        highlight=minimum%20phase#pyfar.dsp.minimum_phase>`_ and ``'linear'`` a
        `linear phase <https://pyfar.readthedocs.io/en/latest/modules/pyfar.\
        dsp.html?highlight=linear%20phase#pyfar.dsp.linear_phase>`_ response.
        The default is ``'minimum'``.

    Returns
    -------
    sofa : sofar Sofa.object
        The DTFs as Sofa object. Can be written to disk with `sofar.write_sofa
        <https://sofar.readthedocs.io/en/latest/sofar.html#\
        sofar.sofar.write_sofa>`_.
    DTFT_inverse : pyfar Signal object
        The inverse directional transfer function as a `pyfar signal object
        <https://pyfar.readthedocs.io/en/latest/classes/pyfar.audio.html#\
        pyfar.classes.audio.Signal>`_
    """

    if isinstance(sofa, str):
        sofa = sf.read_sofa(sofa)
    else:
        sofa = sofa.copy()

    if sofa.GLOBAL_SOFAConventions not in [
            "SimpleFreeFieldHRIR", "GeneralFIR"]:
        raise ValueError(("Sofa object must have the conventions "
                          "SimpleFreeFieldHRIR or GeneralFIR"))

    # get pyfar objects (use pf.io.convert_sofa once released)
    hrir = pf.Signal(sofa.Data_IR, sofa.Data_SamplingRate)
    s_domain, s_convention, s_unit = _sofa_pos(sofa.SourcePosition_Type)
    coordinates = pf.Coordinates(
        sofa.SourcePosition[:, 0],
        sofa.SourcePosition[:, 1],
        sofa.SourcePosition[:, 2],
        domain=s_domain,
        convention=s_convention,
        unit=s_unit)

    # compute DTFT
    dtft = hrir.copy()
    weights = pf.samplings.calculate_sph_voronoi_weights(coordinates)
    if dtft.cshape[1] == 2:
        weights /= 2
    dtft.freq_raw = np.sqrt(np.sum(
        np.abs(dtft.freq_raw)**2 * weights[..., None, None], (0, 1)))

    if smooth_fractions is not None:
        raise NotImplementedError(
            "Will be available with the next pyfar release")

    # get inverse DTFT
    if phase == "minimum":
        raise NotImplementedError(
            "Will be available with the next pyfar release")
    elif phase == "linear":
        dtft_inv = pf.dsp.linear_phase(1 / dtft, dtft.n_samples / 2)
    else:
        raise ValueError(
            f"phase is '{phase}' but must be 'minimum' or 'linear'")

    # compute DTF
    dtf = hrir * dtft_inv

    return dtf, dtft_inv
