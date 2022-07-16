"""
Python tools for Mesh2HRTF including functions to generate SOFA files
containing the HRTF/HRIR data, merge SOFA files containing data for the left
and right ear and generate evaluation grids.
"""
import os
import re
import csv
import warnings
import numpy as np
import glob
import json
import sofar as sf


def output2hrtf(folder=None):
    """
    Process NumCalc output and write data to disk.

    Processing the data is done in the following steps

    1. Read project parameter `from parameters.json`
    2. use :py:func:`~write_output_report` to parse files in
       project_folder/NumCalc/source_*/NC*.out, write project report to
       project_folder/Output2HRTF/report_source_*.csv. Raise a warning if any
       issues were detected and write report_issues.txt to the same folder
    3. Read simulated pressures from project_folder/NumCalc/source_*/be.out.
       This and the following steps are done, even if an issue was detected in
       the previous step
    4. use :py:func:`~reference_hrtf` and :py:func:`~compute_hrir` to save the
       results to SOFA files

    Parameters
    ----------
    folder : str, optional
        The path of the Mesh2HRTF project folder, i.e., the folder containing
        the subfolders EvaluationsGrids, NumCalc, and ObjectMeshes. The
        default, ``None`` uses the current working directory.
    """

    # check input
    if folder is None:
        folder = os.getcwd()

    # check and load parameters, required parameters are:
    # Mesh2HRTF_version, reference, computeHRIRs, speedOfSound, densityOfAir,
    # numSources, sourceType, sourceCenter, sourceArea,
    # numFrequencies, frequencies
    params = os.path.join(folder, "parameters.json")
    if not os.path.isfile(params):
        raise ValueError((
            f"The folder {folder} is not a valid Mesh2HRTF project. "
            "It must contain the file 'parameters.json'"))

    with open(params, "r") as file:
        params = json.load(file)

    # output directory
    if not os.path.exists(os.path.join(folder, 'Output2HRTF')):
        os.makedirs(os.path.join(folder, 'Output2HRTF'))

    # write the project report and check for issues
    print('\n Writing the project report ...')
    found_issues, report = write_output_report(folder)

    if found_issues:
        warnings.warn(report)

    # get the evaluation grids
    print('\n Loading the EvaluationGrids ...')
    evaluationGrids, _ = _read_nodes_and_elements(
        os.path.join(folder, 'EvaluationGrids'))

    # get the object mesh
    print('\n Loading the ObjectMeshes ...')
    objectMeshes, _ = _read_nodes_and_elements(
        os.path.join(folder, 'ObjectMeshes'))

    # Load ObjectMesh data
    pressure = _read_pressure(
        params["numSources"], params["numFrequencies"], folder, 'pBoundary')

    print('\nSaving ObjectMesh data ...')
    cnt = 0
    for mesh in objectMeshes:
        elements = objectMeshes[mesh]["elements"]
        element_data = pressure

        for jj in range(pressure.shape[1]):
            element_data[:, jj, :] = element_data[
                cnt:cnt+elements.shape[0], jj, :]

        file = open(os.path.join(
            folder, "Output2HRTF",
            "ObjectMesh_" + mesh + ".npz"), "w")
        np.savez_compressed(
            file.name, frequencies=params["frequencies"],
            pressure=element_data)
        file.close()

        cnt = cnt + elements.shape[0]

    del pressure, elements, mesh, jj, cnt, element_data

    # Load EvaluationGrid data
    if not len(evaluationGrids) == 0:
        print('\nLoading data for the evaluation grids ...')

        pressure = _read_pressure(
            params["numSources"], params["numFrequencies"],
            folder, 'pEvalGrid')

    # save to struct
    cnt = 0
    for grid in evaluationGrids:
        evaluationGrids[grid]["pressure"] = pressure[
            cnt:cnt+evaluationGrids[grid]["num_nodes"], :, :]

        cnt = cnt + evaluationGrids[grid]["num_nodes"]

    del grid, pressure, cnt

    # process BEM data for writing HRTFs and HRIRs to SOFA files
    for grid in evaluationGrids:

        print(f'\nSaving HRTFs for EvaluationGrid "{grid}" ...\n')

        # get pressure as SOFA object (all following steps are run on SOFA
        # objects. This way they are available to other users as well)
        sofa = _get_sofa_object(
            evaluationGrids[grid]["pressure"],
            evaluationGrids[grid]["nodes"][:, 1:4], "cartesian",
            params["sourceCenter"], "HRTF", params["Mesh2HRTF_Version"],
            frequencies=params["frequencies"])

        # reference to sound pressure at the center of the head
        if params["reference"]:
            sofa = reference_hrtf(
                sofa, params["sourceType"], params["sourceArea"],
                params["speedOfSound"], params["densityOfMedium"], mode="min")

        # write HRTF data to SOFA file
        sf.write_sofa(os.path.join(
            folder, 'Output2HRTF', f'HRTF_{grid}.sofa'), sofa)

        # calculate and write HRIRs
        if params["computeHRIRs"]:
            if not params["reference"]:
                raise ValueError("Computing HRIRs requires prior referencing")

            # calculate shift value (equivalent to a 30 cm shift)
            sampling_rate = round(2*sofa.N[-1])
            n_shift = int(np.round(
                .30 / (1/sampling_rate * params["speedOfSound"])))

            sofa = compute_hrir(sofa, n_shift)
            sf.write_sofa(os.path.join(
                folder, 'Output2HRTF', f'HRIR_{grid}.sofa'), sofa)

    print('Done\n')


def reference_hrtf(sofa, sourceType, sourceArea, speedOfSound, densityOfAir,
                   mode="min"):
    """
    Reference HRTF to the sound pressure in the center of the head. After
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
        "Point source", or "Plane wave"
    sourceArea : array like
        The area of the source is required if `sourceType` is "Both ears" in
        which case `sourceAre` is a list with two values or if `sourceType` is
        "Left ear" or "Right ear" in which case `sourceArea` is a list with one
        value.
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

    if sourceType in {'Both ears', 'Left ear', 'Right ear'}:

        volumeFlow = 0.1 * np.ones(pressure.shape)
        if 'sourceArea':
            # has to be fixed for both ears....
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
                "sourceType 'Both ears', 'Left ear', 'Right ear',"
                "'Point source' and 'Plane wave'."))

    # here we go...
    pressure /= ps
    sofa.Data_Real = np.real(pressure)
    sofa.Data_Imag = np.imag(pressure)

    return sofa


def compute_hrir(sofa, n_shift, sampling_rate=None):
    """
    Compute HRIR from HRTF by means of the inverse Fourier transform.

    This requires the following:

    1. The HRTFs contained in `sofa` have been referenced using
       `reference_HRTFs()`.
    2. HRTF must be available for frequencies f_0 to f_1 in constant steps.
       f_0 must be > 0 Hz and f_1 is assumed to be half the sampling rate.

    HRIRs are computed with the following steps

    1. Add data for 0 Hz. The HRTF at 0 Hz is 1 (0 dB) by definition because
       the HRTF describes the filtering of incoming sound by the human
       anthropometry (which does not change the sound at 0 Hz).
    2. Only the absolute value for the data at half the sampling rate is used.
       Otherwise the next step would produce complex output
    3. The HRTF spectrum is mirrored and the HRIR is obtained through the
       inverse Fourier transform
    4. The HRIRs are circularly shifted by `n_shift` samples to enforce a
       causal system. A good shift value for a sampling rate of 44.1 kHz might
       be between 20 and 40 samples.

    Parameters
    ----------
    sofa : sofar Sofa object, str
       Sofa object containing the sound pressure or filename of a SOFA file to
       be loaded. SOFA object/file must be of the convention
       SimpleFreeFieldHRTF or GeneralTF
    n_shift : int
        Amount the HRIRs are shifted to enforce a causal system.
    sampling_rate : int
        The sampling rate in Hz. The sampling rate can two times any frequency
        for which the HRTF was computed. The default ``None`` assumes the
        sampling rate to be two times the highest frequency for which the HRTF
        was computed.

    Returns
    -------
    sofa : sofar Sofa.object
        HRIRs

    Notes
    -----
    HRIRs for different sampling rates can be generated from a single SOFA file
    if discarding or adding some data.
    """

    if isinstance(sofa, str):
        sofa = sf.read_sofa(sofa)
    else:
        sofa = sofa.copy()

    if sofa.GLOBAL_SOFAConventions not in ["SimpleFreeFieldHRTF", "GeneralTF"]:
        raise ValueError(("Sofa object must have the conventions "
                          "SimpleFreeFieldHRTF or GeneralTF"))

    # check if the frequency vector has the correct format
    frequencies = sofa.N
    if any(np.abs(np.diff(frequencies, 2)) > .1) or frequencies[0] < .1:
        raise ValueError(
            ('The frequency vector must go from f_1 > 0 to'
             'f_2 (half the sampling rate) in equidistant steps.'))

    # get the HRTF
    pressure = sofa.Data_Real + 1j * sofa.Data_Imag

    if sampling_rate is None:
        # detect the sampling rate
        sampling_rate = round(2*frequencies[-1])
    else:
        # check if the specified sampling rate is valid
        idx = np.argmin(np.abs(2 * frequencies - sampling_rate))
        error = np.abs(2 * frequencies[idx] - sampling_rate)
        if np.abs(error) > 1e-6:
            raise ValueError((
                "The specified sampling rate is invalid. It must be two times "
                "any frequency for which the HRTF was computed."))

        # discard frequencies above sampling_rate / 2
        pressure = pressure[..., :idx + 1]
        frequencies = frequencies[:idx + 1]

    # add 0 Hz bin
    pressure = np.concatenate((np.ones((pressure.shape[0],
                               pressure.shape[1], 1)), pressure), axis=2)

    # make sampling_rate/2 real
    pressure[:, :, -1] = np.abs(pressure[:, :, -1])
    # ifft (take complex conjugate because sign conventions differ)
    hrir = np.fft.irfft(np.conj(pressure))

    # shift to make causal
    # (path differences between the origin and the ear are usually
    # smaller than 30 cm but numerical HRIRs show stringer pre-ringing)
    hrir = np.roll(hrir, n_shift, axis=-1)

    sofa = _get_sofa_object(
        hrir, sofa.SourcePosition, sofa.SourcePosition_Type,
        sofa.ReceiverPosition, "HRIR", sofa.GLOBAL_ApplicationVersion,
        sampling_rate=sampling_rate)

    return sofa


def write_output_report(folder=None):
    r"""
    Generate project report from NumCalc output files.

    NumCalc (Mesh2HRTF's numerical core) writes information about the
    simulations to the files `NC*.inp` located under `NumCalc/source_*`. The
    file `NC.inp` exists if NumCalc was ran without the additional command line
    parameters ``-istart`` and ``-iend``. If these parameters were used, there
    is at least one `NC\*-\*.inp`. If this is the case, information from
    `NC\*-\*.inp` overwrites information from NC.inp in the project report.

    .. note::

        The project reports are written to the files
        `Output2HRTF/report_source_*.csv`. If issues were detected, they are
        listed in `Output2HRTF/report_issues.csv`.

    The report contain the following information

    Frequency step
        The index of the frequency.
    Frequency in Hz
        The frequency in Hz.
    NC input
        Name of the input file from which the information was taken.
    Input check passed
        Contains a 1 if the check of the input data passed and a 0 otherwise.
        If the check failed for one frequency, the following frequencies might
        be affected as well.
    Converged
        Contains a 1 if the simulation converged and a 0 otherwise. If the
        simulation did not converge, the relative error might be high.
    Num. iterations
        The number of iterations that were required to converge
    relative error
        The relative error of the final simulation
    Comp. time total
        The total computation time in seconds
    Comp. time assembling
        The computation time for assembling the matrices in seconds
    Comp. time solving
        The computation time for solving the matrices in seconds
    Comp. time post-proc
        The computation time for post-processing the results in seconds


    Parameters
    ----------
    folder : str, optional
        The path of the Mesh2HRTF project folder, i.e., the folder containing
        the subfolders EvaluationsGrids, NumCalc, and ObjectMeshes. The
        default, ``None`` uses the current working directory.

    Returns
    -------
    found_issues : bool
        ``True`` if issues were found, ``False`` otherwise
    report : str
        The report or an empty string if no issues were found
    """

    if folder is None:
        folder = os.getcwd()

    # get sources and number of sources and frequencies
    sources = glob.glob(os.path.join(folder, "NumCalc", "source_*"))
    num_sources = len(sources)

    with open(os.path.join(folder, "parameters.json"), "r") as file:
        params = json.load(file)

    # sort source files (not read in correct order in some cases)
    nums = [int(source.split("_")[-1]) for source in sources]
    sources = np.array(sources)
    sources = sources[np.argsort(nums)]

    # parse all NC*.out files for all sources
    all_files, fundamentals, out, out_names = _parse_nc_out_files(
        sources, num_sources, params["numFrequencies"])

    # write report as csv file
    _write_project_reports(folder, all_files, out, out_names)

    # look for errors
    report = _check_project_report(folder, fundamentals, out)

    found_issues = True if report else False

    return found_issues, report


def _parse_nc_out_files(sources, num_sources, num_frequencies):
    """
    Parse all NC*.out files for all sources.

    This function should never raise a value error, regardless of how mess
    NC*.out files are. Looking for error is done at a later step.

    Parameters
    ----------
    sources : list of strings
        full path to the source folders
    num_sources : int
        number of sources - len(num_sources)
    num_frequencies : int
        number of frequency steps

    Returns
    -------
    out : numpy array
        containing the extracted information for each frequency step
    out_names : list of string
        verbal information about the columns of `out`
    """

    # array for reporting fundamental errors
    fundamentals = []
    all_files = []

    # array for saving the detailed report
    out_names = ["frequency step",         # 0
                 "frequency in Hz",        # 1
                 "NC input file",          # 2
                 "Input check passed",     # 3
                 "Converged",              # 4
                 "Num. iterations",        # 5
                 "relative error",         # 6
                 "Comp. time total",       # 7
                 "Comp. time assembling",  # 8
                 "Comp. time solving",     # 9
                 "Comp. time post-proc."]  # 10
    out = -np.ones((num_frequencies, 11, num_sources))
    # values for steps
    out[:, 0] = np.arange(1, num_frequencies + 1)[..., np.newaxis]
    # values indicating failed input check and non-convergence
    out[:, 3] = 0
    out[:, 4] = 0

    # regular expression for finding a number that can be int or float
    re_number = r"(\d+(?:\.\d+)?)"

    # loop sources
    for ss, source in enumerate(sources):

        # list of NC*.inp files for parsing
        files = glob.glob(os.path.join(source, "NC*.out"))

        # make sure that NC.out is first
        nc_out = os.path.join(source, "NC.out")
        if nc_out in files and files.index(nc_out):
            files = [files.pop(files.index(nc_out))] + files

        # update fundamentals
        fundamentals.append([0 for f in range(len(files))])
        all_files.append([os.path.basename(f) for f in files])

        # frequency steps simulated in each file
        # (code not needed and does not yet catch NCuntil*.inp and NCfrom*.inp)
        # steps = []
        # for file in files:
        #     parts = os.path.basename(file).split('-')
        #     if len(parts) == 1:
        #         steps.append([1, num_frequencies])
        #     else:
        #         steps.append([int(parts[0][2:]), int(parts[1][:-4])])

        # get content from all NC*.inp
        for ff, file in enumerate(files):

            # read the file and join all lines
            with open(file, "r") as f_id:
                lines = f_id.readlines()
            lines = "".join(lines)

            # split header and steps
            lines = lines.split(
                ">> S T E P   N U M B E R   A N D   F R E Q U E N C Y <<")

            # look for fundamental errors
            if len(lines) == 1:
                fundamentals[ss][ff] = 1
                continue

            # parse frequencies (skip header)
            for line in lines[1:]:

                # find frequency step
                idx = re.search(r'Step \d+,', line)
                if idx:
                    step = int(line[idx.start()+5:idx.end()-1])

                # write number of input file (replaced by string later)
                out[step-1, 2, ss] = ff

                # find frequency
                idx = re.search(f'Frequency = {re_number} Hz', line)
                if idx:
                    out[step-1, 1, ss] = float(
                        line[idx.start()+12:idx.end()-3])

                # check if the input data was ok
                if "Too many integral points in the theta" not in line:
                    out[step-1, 3, ss] = 1

                # check and write convergence
                if 'Maximum number of iterations is reached!' not in line:
                    out[step-1, 4, ss] = 1

                # check iterations
                idx = re.search(r'number of iterations = \d+,', line)
                if idx:
                    out[step-1, 5, ss] = int(line[idx.start()+23:idx.end()-1])

                # check relative error
                idx = re.search('relative error = .+', line)
                if idx:
                    out[step-1, 6, ss] = float(line[idx.start()+17:idx.end()])

                # check time stats
                # -- assembling
                idx = re.search(
                    r'Assembling the equation system         : \d+',
                    line)
                if idx:
                    out[step-1, 8, ss] = float(line[idx.start()+41:idx.end()])

                # -- solving
                idx = re.search(
                    r'Solving the equation system            : \d+',
                    line)
                if idx:
                    out[step-1, 9, ss] = float(line[idx.start()+41:idx.end()])

                # -- post-pro
                idx = re.search(
                    r'Post processing                        : \d+',
                    line)
                if idx:
                    out[step-1, 10, ss] = float(line[idx.start()+41:idx.end()])

                # -- total
                idx = re.search(
                    r'Total                                  : \d+',
                    line)
                if idx:
                    out[step-1, 7, ss] = float(line[idx.start()+41:idx.end()])

    return all_files, fundamentals, out, out_names


def _write_project_reports(folder, all_files, out, out_names):
    """
    Write project report to disk at folder/Output2HRTF/report_source_*.csv

    For description of input parameter refer to write_output_report and
    _parse_nc_out_files
    """

    # loop sources
    for ss in range(out.shape[2]):

        report = ", ".join(out_names) + "\n"

        # loop frequencies
        for ff in range(out.shape[0]):
            f = out[ff, :, ss]
            report += (
                f"{int(f[0])}, "                # frequency step
                f"{float(f[1])}, "              # frequency in Hz
                f"{all_files[ss][int(f[2])]},"  # NC*.inp file
                f"{int(f[3])}, "                # input check
                f"{int(f[4])}, "                # convergence
                f"{int(f[5])}, "                # number of iterations
                f"{float(f[6])}, "              # relative error
                f"{int(f[7])}, "                # total computation time
                f"{int(f[8])}, "                # assembling equations time
                f"{int(f[9])}, "                # solving equations time
                f"{int(f[10])}\n"               # post-processing time
                )

        # write to disk
        report_name = os.path.join(
            folder, "Output2HRTF", f"report_source_{ss + 1}.csv")
        with open(report_name, "w") as f_id:
            f_id.write(report)


def _check_project_report(folder, fundamentals, out):

    # return if there are no fundamental errors or other issues
    if not all([all(f) for f in fundamentals]) and not np.any(out == -1) \
            and np.all(out[:, 3:5]):
        return

    # report detailed errors
    report = ""

    for ss in range(out.shape[2]):

        # currently we detect frequencies that were not calculated and
        # frequencies with convergence issues
        missing = "Frequency steps that were not calculated:\n"
        input_test = "Frequency steps with bad input:\n"
        convergence = "Frequency steps that did not converge:\n"

        any_missing = False
        any_input_failed = False
        any_convergence = False

        # loop frequencies
        for ff in range(out.shape[0]):

            f = out[ff, :, ss]

            # no value for frequency
            if f[1] == -1:
                any_missing = True
                missing += f"{int(f[0])}, "
                continue

            # input data failed
            if f[3] == 0:
                any_input_failed = True
                input_test += f"{int(f[0])}, "

            # convergence value is zero
            if f[4] == 0:
                any_convergence = True
                convergence += f"{int(f[0])}, "

        if any_missing or any_input_failed or any_convergence:
            report += f"Detected issues for source {ss+1}\n"
            report += "----------------------------\n"
            if any_missing:
                report += missing[:-2] + "\n\n"
            if any_input_failed:
                report += input_test[:-2] + "\n\n"
            if any_convergence:
                report += convergence[:-2] + "\n\n"

    if not report:
        report = ("\nDetected unknown issues\n"
                  "-----------------------\n"
                  "Check the project reports in Output2HRTF,\n"
                  "and the NC*.inp files in NumCalc/source_*\n\n")

    report += ("For more information check Output2HRTF/report_source_*.csv "
               "and the NC*.inp files located at NumCalc/source_*")

    # write to disk
    report_name = os.path.join(
        folder, "Output2HRTF", "report_issues.txt")
    with open(report_name, "w") as f_id:
        f_id.write(report)

    return report


def _read_nodes_and_elements(data):
    """
    Read the nodes and elements of the evaluation grids or object meshes.
    """
    if os.path.basename(data) not in ['EvaluationGrids', 'ObjectMeshes']:
        raise ValueError('data must be EvaluationGrids or ObjectMeshes!')

    grids = {}
    gridsList = os.listdir(data)
    gridsNumNodes = 0

    for grid in gridsList:
        tmpNodes = np.loadtxt(os.path.join(
            data, grid, 'Nodes.txt'),
            delimiter=' ', skiprows=1)

        tmpElements = np.loadtxt(os.path.join(
            data, grid, 'Elements.txt'),
            delimiter=' ', skiprows=1)

        grids[grid] = {
            "nodes": tmpNodes,
            "elements": tmpElements,
            "num_nodes": tmpNodes.shape[0]}

        gridsNumNodes += grids[grid]['num_nodes']

    return grids, gridsNumNodes


def _read_pressure(numSources, numFrequencies, folder, data):
    """Read the sound pressure on the object meshes or evaluation grid."""
    pressure = []

    if not (data == 'pBoundary' or data == 'pEvalGrid'):
        raise ValueError('data must be pBoundary or pEvalGrid!')

    print('\n Loading %s data ...' % data)
    for source in range(numSources):
        print('\n    Source %d ...' % (source+1))

        tmpFilename = os.path.join(
            folder, 'NumCalc', f'source_{source+1}', 'be.out')
        tmpPressure = _output_to_hrtf_load(
            tmpFilename, data, numFrequencies)

        pressure.append(tmpPressure)

    pressure = np.transpose(np.array(pressure), (2, 0, 1))

    return pressure


def _get_sofa_object(data, source_position, source_position_type,
                     receiver_position, mode, Mesh2HRTF_version,
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

    # get source coordinates in spherical convention
    if source_position_type == "cartesian":
        xyz = source_position

        radius = np.sqrt(xyz[:, 0]**2 + xyz[:, 1]**2 + xyz[:, 2]**2)
        z_div_r = np.where(radius != 0, xyz[:, 2] / radius, 0)
        elevation = 90 - np.arccos(z_div_r) / np.pi * 180
        azimuth = np.mod(np.arctan2(xyz[:, 1], xyz[:, 0]), 2 * np.pi) \
            / np.pi * 180

        source_position = np.concatenate(
            (azimuth[..., np.newaxis],
             elevation[..., np.newaxis],
             radius[..., np.newaxis]), axis=-1)

    # get number of sources from data
    numSources = data.shape[1]

    # create empty SOFA object
    if mode == "HRTF":
        convention = "SimpleFreeFieldHRTF" if numSources == 2 else "GeneralTF"
    else:
        convention = "SimpleFreeFieldHRIR" if numSources == 2 else "GeneralFIR"

    sofa = sf.Sofa(convention)

    # write meta data
    sofa.GLOBAL_ApplicationName = 'Mesh2HRTF'
    sofa.GLOBAL_ApplicationVersion = Mesh2HRTF_version
    sofa.GLOBAL_History = "numerically simulated data"

    # Source and receiver data
    sofa.SourcePosition = source_position
    sofa.SourcePosition_Units = "degree, degree, meter"
    sofa.SourcePosition_Type = "spherical"

    sofa.ReceiverPosition = receiver_position
    sofa.ReceiverPosition_Units = "meter"
    sofa.ReceiverPosition_Type = "cartesian"

    # HRTF/HRIR data
    if mode == "HRTF":
        sofa.Data_Real = np.real(data)
        sofa.Data_Imag = np.imag(data)
        sofa.N = np.array(frequencies).flatten()
    else:
        sofa.Data_IR = data
        sofa.Data_SamplingRate = sampling_rate
        sofa.Data_Delay = np.zeros((1, data.shape[1]))

    return sofa


def _output_to_hrtf_load(foldername, filename, numFrequencies):
    """
    Load results of the BEM calculation.

    Parameters
    ----------
    foldername : string
        The folder from which the data is loaded. The data to be read is
        located in the folder be.out inside NumCalc/source_*
    filename : string
        The kind of data that is loaded

        pBoundary
            The sound pressure on the object mesh
        vBoundary
            The sound velocity on the object mesh
        pEvalGrid
            The sound pressure on the evaulation grid
        vEvalGrid
            The sound velocity on the evaluation grid
    numFrequencies : int
        the number of simulated frequencies

    Returns
    -------
    data : numpy array
        Pressure or velocity values of shape (numFrequencies, numEntries)
    """

    # ---------------------check number of header and data lines---------------

    idx1 = 0
    idx2 = 0
    ii = 0
    with open(os.path.join(foldername, 'be.1', filename)) as file:
        line = csv.reader(file, delimiter=' ', skipinitialspace=True)
        for li in line:
            ii += 1
            if li not in (None, ""):
                if len(li) == 3:
                    idx1 = idx2
                    idx2 = int(li[0])
                    if idx2 - idx1 == 1:
                        numHeaderlines_BE = ii - 2
                        numDatalines = sum(1 for ll in line) + 2
                        break

    # ------------------------------load data----------------------------------
    data = np.zeros((numFrequencies, numDatalines), dtype=complex)

    for ii in range(numFrequencies):
        tmpData = []
        current_file = os.path.join(foldername, 'be.%d' % (ii+1), filename)
        with open(current_file) as file:
            line = csv.reader(file, delimiter=' ', skipinitialspace=True)
            for li in line:
                if line.line_num > numHeaderlines_BE:
                    tmpData.append(complex(float(li[1]), float(li[2])))

        data[ii, :] = tmpData if tmpData else np.nan

    return data
