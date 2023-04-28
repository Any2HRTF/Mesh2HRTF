import sofar as sf
import pyfar as pf


def resample_sofa_file(sofa, sampling_rate, frac_limit=None, post_filter=True):
    """Resample Sofa file to new sampling rate.

    The SciPy function ``scipy.signal.resample_poly`` is used for resampling.
    The resampling ratio ``L = sampling_rate/signal.sampling_rate``
    is approximated by a fraction of two integer numbers `up/down` to first
    upsample the signal by `up` and then downsample by `down`. This way `up`
    and `down` are smaller than the respective new and old sampling rates.

    .. note ::

        `sampling_rate` should be divisible by 10, otherwise it can cause an
        infinite loop in the ``resample_poly`` function.

        The amplitudes of the resampled signal can match the amplitude of the
        input signal in the time or frequency domain. See the parameter
        `match_amplitude` and the examples for more information.

    Parameters
    ----------
    sofa : Sofa object, str
        Input data to be resampled as a sofar Sofa object or filename of the
        Sofa file on disk.
    sampling_rate : number
        The new sampling rate in Hz
    frac_limit : int
        Limit the denominator for approximating the resampling factor `L`
        (see above). This can be used in case the resampling gets stuck in an
        infinite loop (see note above) at the potenital cost of not exactly
        realizing the target sampling rate.

        The default is ``None``, which uses ``frac_limit = 1e6``.
    post_filter : bool, optional
        In some cases the up-sampling causes artifacts above the Nyquist
        frequency of the input signal, i.e., ``signal.sampling_rate/2``. If
        ``True`` the artifacts are suppressed by applying a zero-phase Elliptic
        filter with a pass band ripple of 0.1 dB, a stop band attenuation of 60
        dB. The pass band edge frequency is ``signal.sampling_rate/2``. The
        stop band edge frequency is the minimum of 1.05 times the pass band
        frequency and the new Nyquist frequency (``sampling_rate/2``). The
        default is ``True``. Note that this is only applied in case of
        up-sampling.

    Returns
    -------
    sofa_resampled : Sofa object
        The resampled version of the input data with a length of
        ``up/down * sofa.get_dimension("N)`` samples. The object can be written
        to disk with ``sofar.write_sofa()``.
    """

    # check input
    if isinstance(sofa, str):
        sofa = sf.read_sofa(sofa)
    elif not isinstance(sofa, sf.Sofa):
        raise TypeError(
            "sofa must be a sofar Sofa object or path to a sofa file")

    if not sofa.GLOBAL_DataType.startswith("FIR"):
        raise TypeError(("The DataType of the sofa file must be FIR "
                         f"but is {sofa.GLOBAL_DataType}"))

    # resample
    data, *_ = pf.io.convert_sofa(sofa)
    data = pf.dsp.resample(data, sampling_rate, 'freq', frac_limit=frac_limit,
                           post_filter=post_filter)

    # return as Sofa object
    sofa = sofa.copy()
    sofa.Data_IR = data.time
    sofa.Data_SamplingRate = sampling_rate

    return sofa
