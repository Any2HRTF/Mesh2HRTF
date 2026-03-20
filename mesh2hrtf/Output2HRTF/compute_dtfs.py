import numpy as np
import sofar as sf
import pyfar as pf
import spharpy


def compute_dtfs(
        sofa, smooth_fractions=None, phase="minimum", weights="equal"):
    r"""
    Compute directional transfer functions from HRIRs.

    The directional transfer functions (DTFs) are computed by a spectral
    division of the HRTFs by the diffuse field transfer function (DFTF)

    .. math::

        \mathrm{DTF} = \frac{\mathrm{HRTF}}{\mathrm{DFTF}}

    where the DFTF is computed as the energetic average across source positions
    and ears

    .. math::

        \mathrm{DFTF} = \sqrt{\sum_{n=0}^{N-1} |\mathrm{HRTF}|^2_{l,n}\,w_n +
        \sum_{n=0}^{N-1} |\mathrm{HRTF}|^2_{r,n}\,w_n}

    The index :math:`n` denotes the source position, :math:`|\cdot|` the
    absolute spectrum, and :math:`\sum_n w_n=1` normalized area weights for
    numerical integration (see below). The average across ears is made to not
    alter binaural cues.

    Parameters
    ----------
    sofa : sofar Sofa object, str
       Sofa object containing the sound pressure or filename of a SOFA file to
       be loaded. SOFA object/file must be of the convention
       SimpleFreeFieldHRIR or GeneralFIR
    smooth_fractions : number, None, optional
        Apply fractional octave smoothing to the DFTF using
        :py:func:`pyfar.dsp.smooth_fractional_octave` to the DFTF.
        E.g. a value of ``3`` applies third octave smoothing and a value of
        ``1`` applies octave smoothing. The default ``None`` does not apply any
        smoothing.
    phase : string, optional
        Define the phase of the inverse DFTF.

        ``'minimum'``
            generate a minimum phase response using
            :py:func:`pyfar.dsp.minimum_phase`
        ``'linear'``
            generate a linear phase response using
            :py:func:`pyfar.dsp.linear_phase`
        ``'zero'``
            generates a zero phase response.

        The default is ``'minimum'``.
    weights : optional
        Define the weights used for the numerical integration

        ``'equal'``
            Uses equal weights across source positions
        ``'voronoi'``
            Computes spherical Voronoi weights using
            :py:func:`spharpy.samplings.calculate_sampling_weights`. This
            requires that all points defined by `sofa.SourcePosition` have
            the same radius.
        array like
            Uses the weights provided in a list or numpy array. The size of the
            array like must agree with the number of HRTFs

        The default is ``'equal'``

    Returns
    -------
    sofa : sofar Sofa.object
        The DTFs as :py:class:`sofar.Sofa` object. Can be written to disk with
        :py:func:`sofar.write_sofa`.
    DTFT_inverse : pyfar Signal object
        The inverse diffuse field transfer function as a
        :py:class:`pyfar.Signal`.
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
    hrir, coordinates, _ = pf.io.convert_sofa(sofa)

    # get or check the weights
    if weights == "equal":
        weights = None
    elif weights == "voronoi":
        weights = spharpy.samplings.calculate_sampling_weights(
            spharpy.SamplingSphere.from_coordinates(coordinates))
        weights = weights[..., None]
    elif isinstance(weights, (list, np.ndarray)):
        weights = np.asarray(weights).flatten()
        if weights.size != hrir.cshape[0]:
            raise ValueError((
                f"{weights.size} provided but {hrir.cshape[0]} expected "
                "(number of HRIRs)"))
        weights = weights[..., None]
    else:
        raise ValueError("weights must be 'equal', 'voronoi' or an array like")

    # compute DFTF
    dftf = pf.dsp.average(
        hrir, "power", caxis=(0, 1), weights=weights)

    if smooth_fractions is not None:
        dftf, _ = pf.dsp.smooth_fractional_octave(dftf, smooth_fractions)

    # get inverse DFTF
    if phase == "minimum":
        dftf_inv = pf.dsp.minimum_phase(1 / dftf, truncate=False)
    elif phase == "linear":
        dftf_inv = pf.dsp.linear_phase(1 / dftf, dftf.n_samples / 2)
    elif phase == "zero":
        dftf_inv = 1 / dftf
    else:
        raise ValueError(
            f"phase is '{phase}' but must be 'minimum' or 'linear'")

    # compute DTF
    dtf = hrir * dftf_inv
    sofa.Data_IR = dtf.time

    return sofa, dftf_inv
