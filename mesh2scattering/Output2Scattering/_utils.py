import os
import re
import numpy as np
import glob
import csv
import sofar as sf


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
    root = os.path.abspath(os.path.join(environ, os.pardir))
    return root


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

        # list of NC*.out files for parsing
        files = glob.glob(os.path.join(source, "NC*.out"))

        # make sure that NC.out is first
        nc_out = os.path.join(source, "NC.out")
        if nc_out in files and files.index(nc_out):
            files = [files.pop(files.index(nc_out))] + files

        # update fundamentals
        fundamentals.append([0 for f in range(len(files))])
        all_files.append([os.path.basename(f) for f in files])

        # get content from all NC*.out
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
                f"{all_files[ss][int(f[2])]},"  # NC*.out file
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
                  "and the NC*.out files in NumCalc/source_*\n\n")

    report += ("For more information check Output2HRTF/report_source_*.csv "
               "and the NC*.out files located at NumCalc/source_*")

    # write to disk
    report_name = os.path.join(
        folder, "Output2HRTF", "report_issues.txt")
    with open(report_name, "w") as f_id:
        f_id.write(report)

    return report


def _read_nodes_and_elements(folder, objects=None):
    """
    Read the nodes and elements of the evaluation grids or object meshes.

    Parameters
    ----------
    folder : str
        Folder containing the object. Must end with EvaluationGrids or
        Object Meshes
    objects : str, options
        Name of the object. The default ``None`` reads all objects in folder

    Returns
    -------
    grids : dict
        One item per object (with the item name being the object name). Each
        item has the sub-items `nodes`, `elements`, `num_nodes`, `num_elements`
    gridsNumNodes : int
        Number of nodes in all grids
    """
    # check input
    if os.path.basename(folder) not in ['EvaluationGrids', 'ObjectMeshes']:
        raise ValueError('folder must be EvaluationGrids or ObjectMeshes!')

    if objects is None:
        objects = os.listdir(folder)
        # discard hidden folders that might occur on Mac OS
        objects = [o for o in objects if not o.startswith('.')]
    elif isinstance(objects, str):
        objects = [objects]

    grids = {}
    gridsNumNodes = 0

    for grid in objects:
        tmpNodes = np.loadtxt(os.path.join(
            folder, grid, 'Nodes.txt'),
            delimiter=' ', skiprows=1, dtype=np.float64)

        tmpElements = np.loadtxt(os.path.join(
            folder, grid, 'Elements.txt'),
            delimiter=' ', skiprows=1, dtype=np.float64)

        grids[grid] = {
            "nodes": tmpNodes,
            "elements": tmpElements,
            "num_nodes": tmpNodes.shape[0],
            "num_elements": tmpElements.shape[0]}

        gridsNumNodes += grids[grid]['num_nodes']

    return grids, gridsNumNodes


def _read_numcalc_data(numSources, numFrequencies, folder, data):
    """Read the sound pressure on the object meshes or evaluation grid."""
    pressure = []

    if data not in ['pBoundary', 'pEvalGrid', 'vBoundary', 'vEvalGrid']:
        raise ValueError(
            'data must be pBoundary, pEvalGrid, vBoundary, or vEvalGrid')

    for source in range(numSources):

        tmpFilename = os.path.join(
            folder, 'NumCalc', f'source_{source+1}', 'be.out')
        tmpPressure, indices = _load_results(
            tmpFilename, data, numFrequencies)

        pressure.append(tmpPressure)

    pressure = np.transpose(np.array(pressure), (2, 0, 1))

    return pressure, indices


def _load_results(foldername, filename, numFrequencies):
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
            The sound pressure on the evaluation grid
        vEvalGrid
            The sound velocity on the evaluation grid
    numFrequencies : int
        the number of simulated frequencies

    Returns
    -------
    data : numpy array
        Pressure or abs velocity values of shape (numFrequencies, numEntries)
    """

    # ---------------------check number of header and data lines---------------
    current_file = os.path.join(foldername, 'be.1', filename)
    numDatalines = None
    with open(current_file) as file:
        line = csv.reader(file, delimiter=' ', skipinitialspace=True)
        for idx, li in enumerate(line):
            # read number of data points and head lines
            if len(li) == 2 and not li[0].startswith("Mesh"):
                numDatalines = int(li[1])

            # read starting index
            elif numDatalines and len(li) > 2:
                start_index = int(li[0])
                break

    # ------------------------------load data----------------------------------
    dtype = complex if filename.startswith("p") else float
    data = np.zeros((numFrequencies, numDatalines), dtype=dtype)

    for ii in range(numFrequencies):
        tmpData = []
        current_file = os.path.join(foldername, 'be.%d' % (ii+1), filename)
        with open(current_file) as file:

            line = csv.reader(file, delimiter=' ', skipinitialspace=True)

            for li in line:

                # data lines have 3 ore more entries
                if len(li) < 3 or li[0].startswith("Mesh"):
                    continue

                if filename.startswith("p"):
                    tmpData.append(complex(float(li[1]), float(li[2])))
                elif filename == "vBoundary":
                    tmpData.append(np.abs(complex(float(li[1]), float(li[2]))))
                elif filename == "vEvalGrid":
                    tmpData.append(np.sqrt(
                        np.abs(complex(float(li[1]), float(li[2])))**2 +
                        np.abs(complex(float(li[3]), float(li[4])))**2 +
                        np.abs(complex(float(li[5]), float(li[6])))**2))

        data[ii, :] = tmpData if tmpData else np.nan

    return data, np.arange(start_index, numDatalines + start_index)


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
