import os
import glob
import time
import psutil
import subprocess
import numpy as np
import mesh2scattering as m2s


def manage_numcalc(project_path=os.getcwd(), numcalc_path=None,
                   max_ram_load=None, ram_safety_factor=1.05, max_cpu_load=90,
                   max_instances=psutil.cpu_count(), wait_time=15,
                   starting_order='alternate', confirm_errors=False):
    """
    Run NumCalc on one or multiple Mesh2HRTF project folders.

    This script monitors the RAM and CPU usage and starts a new NumCalc
    instance whenever enough resources are available. The required RAM for each
    frequency step is estimated using NumCalc's `estimate_ram` option. A log
    file is written to the `project_path` containing detailed information on
    the launched frequency steps, available resources, and detected errors.

    .. note ::

        `manage_numcalc` can also be launched by running the python script
        `manage_numcalc_script.py` contained in the subfolder
        `mesh2hrtf/NumCalc` of the Mesh2HRTF Git repository.

    Parameters
    ----------
    project_path : str, optional
        The directory to simulate: It can be path to either
        1- directory that contains multiple Mesh2HRTF project folders or
        2- one Mesh2HRTF project folder (folder containing "parameters.json").
        The default is os.getcwd()
    numcalc_path : str, optional
        On Unix, this is the path to the NumCalc binary (by default 'NumCalc'
        is used). On Windows, this is the path to the folder
        'NumCalc_WindowsExe' from
        https://sourceforge.net/projects/mesh2hrtf-tools/ (by default the
        `project_path` is searched for this folder)
    max_ram_load : number, optional
        The RAM that can maximally be used in GB. New NumCalc instances are
        only started if enough RAM is available. The default ``None`` uses all
        available RAM will be used.
    ram_safety_factor : number, optional
        A safety factor that is applied to the estimated RAM consumption. The
        estimate is obtained using NumCalc -estimate_ram. The default of
        ``1.05`` would for example assume that 10.5 GB ram are needed if a RAM
        consumption of 10 GB was estimated by NumCalc.
    max_cpu_load : number, optional
        Maximum allowed CPU load in percent. New instances are only launched if
        the current CPU load is below this value. The default is 90 percent.
    max_instances : int, optional
        The maximum numbers of parallel NumCalc instances. By default a new
        instance is launched until the number of available CPU cores given by
        ``psutil.cpu_count()`` is reached.
    wait_time : int, optional
        Delay in seconds for waiting until the RAM and CPU usage is checked
        after launching a NumCalc instance. This has to be sufficiently large
        for the RAM and CPU to be fully used by the started NumCalc instance.
        The default is 15. After this initial wait time, the resources are
        checked every second. And the next instance is started, once enough
        resources are available.
    starting_order : str, optional
        Control the order in which the frequency steps are launched.

        ``'high'``
            Always launches the step with the highest possible memory
            consumption.
        ``'low'``
            Alsways launches the step with the lowest estimated memory
            consumption
        ``'alternate'`` (default)
            mixes the two approaches above.

    confirm_errors : bool, optional
        If True, manage_numcalc waits for user input in case an error "
        occurs. The default false exits the function immediately if an error
        occurs.
    """

    # log_file initialization -------------------------------------------------
    current_time = time.strftime("%Y_%m_%d_%H-%M-%S", time.localtime())
    log_file = os.path.join(
        project_path, f"manage_numcalc_{current_time}.txt")

    # remove old log-file
    if os.path.isfile(log_file):
        os.remove(log_file)

    # default values ----------------------------------------------------------
    if os.name == "nt":
        numcalc_path = "Searching for NumCalc_WindowsExe"\
             if numcalc_path is None else numcalc_path
    else:
        numcalc_path = "NumCalc" if numcalc_path is None else numcalc_path

    ram_info = psutil.virtual_memory()
    max_ram_load = ram_info.total / 1073741824 \
        if max_ram_load is None else max_ram_load

    # helping variables -------------------------------------------------------

    # RAM that should not be used
    ram_offset = max([0, ram_info.total / 1073741824 - max_ram_load])

    # trick to get colored print-outs   https://stackoverflow.com/a/54955094
    text_color_red = '\033[31m'
    text_color_green = '\033[32m'
    text_color_reset = '\033[0m'

    # wait time in seconds before checking resources again if we are busy
    # (this is not the wait time directly after a NumCalc instance was
    # launched. That is given by wait time!)
    wait_time_busy = 1

    # check input -------------------------------------------------------------
    if max_instances > psutil.cpu_count():
        _raise_error(
            (f"max_instances is {max_instances} but can not be larger than "
             f"{psutil.cpu_count()} (The number of logical CPUs)"),
            text_color_red, log_file, confirm_errors)

    # Detect what the project_path or "getcwd()" is pointing to:
    if os.path.isfile(os.path.join(project_path, 'parameters.json')):
        # project_path is a Mesh2HRTF project folder
        all_projects = [project_path]
        log_file = os.path.join(project_path, log_file)
    else:
        # project_path contains multiple Mesh2HRTF project folders
        all_projects = []  # list of project folders to execute
        for subdir in os.listdir(project_path):
            if os.path.isdir(os.path.join(project_path, subdir,
                                          'ObjectMeshes', 'Reference')):
                all_projects.append(os.path.join(project_path, subdir))

        log_file = os.path.join(project_path, log_file)

        # stop if no project folders were detected
        if len(all_projects) == 0:
            message = ("manage_numcalc could not detect any Mesh2HRTF "
                       f"projects at project_path={project_path}")
            _raise_error(message, text_color_red, log_file, confirm_errors)

    # echo input parameters ---------------------------------------------------
    current_time = time.strftime("%b %d %Y, %H:%M:%S", time.localtime())
    message = ("\nStarting manage_numcalc with the following arguments "
               f"[{current_time}]\n")
    message += "-" * (len(message) - 2) + "\n"
    message += (
        f"project_path: {project_path}\n"
        f"numcalc_path: {numcalc_path}\n"
        f"max_ram_load: {max_ram_load}\n"
        f"ram_safety_factor: {ram_safety_factor}\n"
        f"max_cpu_load: {max_cpu_load}\n"
        f"max_instances: {max_instances}\n"
        f"wait_time: {wait_time}\n"
        f"starting_order: {starting_order}\n"
        f"confirm_errors: {confirm_errors}\n")

    _print_message(message, text_color_reset, log_file)

    # Check for NumCalc executable --------------------------------------------
    if os.name == 'nt':  # Windows detected

        # files that are needed to execute NumCalc
        NumCalc_runtime_files = ['NumCalc.exe', 'libgcc_s_seh-1.dll',
                                 'libstdc++-6.dll', 'libwinpthread-1.dll']

        if os.path.isdir(os.path.join(all_projects[0], 'NumCalc_WindowsExe')):
            # located inside the project folder
            numcalc_path = os.path.join(all_projects[0], 'NumCalc_WindowsExe')
        elif os.path.isdir(os.path.join(os.path.dirname(all_projects[0]),
                                        'NumCalc_WindowsExe')):
            # located inside the folder that contains all Mesh2HRTF projects
            numcalc_path = os.path.join(
                os.path.dirname(all_projects[0]), 'NumCalc_WindowsExe')
        elif os.path.isfile(os.path.join(all_projects[0],
                                         NumCalc_runtime_files[0])):
            # located directly in the project folder.
            numcalc_path = os.path.join(all_projects[0])
        else:
            # try path provided as it is
            pass

        # Check that each required runtime file is present:
        for calc_file in NumCalc_runtime_files:
            if not os.path.isfile(os.path.join(numcalc_path, calc_file)):
                message = (
                    f"The file {calc_file} is missing or manage_numcalc "
                    f"did not find the containing folder 'NumCalc_WindowsExe'")
                _raise_error(message, text_color_red, log_file, confirm_errors)

        # full path to the NumCalc executable
        numcalc_executable = os.path.join(numcalc_path, "NumCalc.exe")

        del calc_file, NumCalc_runtime_files
    else:
        if not numcalc_path.endswith("NumCalc"):
            _raise_error(
                "numcalc_path must end with 'NumCalc'", text_color_red,
                log_file, confirm_errors)
        p = subprocess.Popen(
            f"command -v {numcalc_path}", stdout=subprocess.PIPE, shell=True)
        if not len(p.stdout.read()):
            _raise_error(
                f"NumCalc executable does not exist at {numcalc_path}",
                text_color_red, log_file, confirm_errors)
        numcalc_executable = numcalc_path
        numcalc_path = os.path.dirname(numcalc_path)

    # echo the used NumCalc executable
    _print_message(f"NumCalc executable: {numcalc_executable}\n",
                   text_color_reset, log_file)

    # Check all projects that may need to be executed -------------------------
    projects_to_run = []
    message = ("\nPer project summary of instances that will be run\n"
               "-------------------------------------------------\n")

    message += f"Detected {len(all_projects)} Mesh2HRTF projects in\n"
    message += f"{os.path.dirname(log_file)}\n\n"

    for project in all_projects:
        all_instances, instances_to_run, *_ = _check_project(
            project, numcalc_executable, log_file)

        if instances_to_run is not None:
            projects_to_run.append(project)
            message += (
                f"{len(instances_to_run)}/{len(all_instances)} frequency "
                f"steps to run in {os.path.basename(project)}\n")
        else:
            message += f"{os.path.basename(project)} is already complete\n"

    _print_message(message, text_color_reset, log_file)

    # loop to process all projects --------------------------------------------
    for pp, project in enumerate(projects_to_run):

        current_time = time.strftime("%b %d %Y, %H:%M:%S", time.localtime())

        # Get number of instances in project and estimate their RAM consumption
        root_NumCalc = os.path.join(project, 'NumCalc')
        all_instances, instances_to_run, source_counter = \
            _check_project(project, numcalc_executable, log_file)
        total_nr_to_run = instances_to_run.shape[0]

        # Status printouts:
        message = (f"Started {os.path.basename(project)} "
                   f"({pp + 1}/{len(projects_to_run)}, {current_time})")
        message = "\n" + message + "\n" + "-" * len(message) + "\n"
        if total_nr_to_run:
            message += (
                f"Running {total_nr_to_run}/{len(all_instances)} unfinished "
                "frequency steps in the project\n")
        else:
            message += (
                "All NumCalc simulations in this project are complete")
            _print_message(message, text_color_reset, log_file)
            continue

        _print_message(message, text_color_reset, log_file)

        # sort instances according to RAM consumption (lowest first)
        instances_to_run = instances_to_run[np.argsort(instances_to_run[:, 3])]

        # check if available memory is enough for running the instance with the
        # highest memory consumption without ever exceeding 100% of RAM.
        ram_available, ram_used = _get_current_ram(ram_offset)
        if ram_available < instances_to_run[-1, 3] * ram_safety_factor:
            # note: it IS possible to run simulations that use even more than
            # 100% of available system RAM - only the performance will be poor.
            _raise_error((
                f"Stop - not sufficient free RAM for this simulation project: "
                f"Available RAM is {round(ram_available, 2)} GB, but frequency"
                f" step {int(instances_to_run[-1, 1])} of source "
                f"{int(instances_to_run[-1, 0])} requires "
                f"{round(instances_to_run[-1, 3] * ram_safety_factor, 2)} "
                "GB."), text_color_red, log_file, confirm_errors)

        # assure highest first if demanded
        if starting_order != "low":
            instances_to_run = np.flip(instances_to_run, axis=0)

        # main loop for starting instances
        started_instance = False  # init

        while instances_to_run.shape[0]:

            ram_required = np.min(instances_to_run[:, 3]) * ram_safety_factor

            # current time and resources
            current_time = time.strftime(
                "%b %d %Y, %H:%M:%S", time.localtime())
            ram_available, ram_used = _get_current_ram(ram_offset)
            cpu_load = psutil.cpu_percent(.1)
            running_instances = _numcalc_instances()

            # wait if
            # - CPU usage too high
            # - number of running instances is too large
            # - not enough RAM available
            if cpu_load > max_cpu_load \
                    or running_instances >= max_instances \
                    or ram_available < ram_required:

                # print message (only done once between launching instances)
                if started_instance:
                    _print_message(
                        (f"\n... waiting for resources (checking every "
                         f"second, {current_time}):\n"
                         f" {running_instances} NumCalc instances running ("
                         f"{cpu_load}% CPU load)\n"
                         f" {round(ram_available, 2)} GB RAM available ("
                         f"{round(ram_required, 2)} GB RAM needed next)\n"),
                        text_color_reset, log_file)
                    started_instance = False

                # wait and continue
                time.sleep(wait_time_busy)
                continue

            # find frequency step with the highest possible RAM consumption
            for idx, ram_required in enumerate(instances_to_run[:, 3]):
                if ram_required <= ram_available:
                    break

            # start new NumCalc instance
            source = int(instances_to_run[idx, 0])
            step = int(instances_to_run[idx, 1])
            progress = total_nr_to_run - instances_to_run.shape[0] + 1
            message = (
                f"{progress}/{total_nr_to_run} starting instance from: "
                f"{os.path.basename(project)} (source {source}, step {step}, "
                f"{current_time})")
            _print_message(message, text_color_reset, log_file)

            # new working directory
            cwd = os.path.join(root_NumCalc, "source_" + str(source))

            if os.name == 'nt':  # Windows detected
                # create a log file for all print-outs
                LogFileHandle = open(
                    os.path.join(cwd, "NC{step}-{step}_log.txt"), "w")
                # run NumCalc and route all printouts to a log file
                subprocess.Popen(
                    f"{numcalc_executable} -istart {step} -iend {step}",
                    stdout=LogFileHandle, cwd=cwd)

            else:  # elif os.name == 'posix': Linux or Mac detected
                # run NumCalc and route all printouts to a log file
                subprocess.Popen((
                    f"{numcalc_executable} -istart {step} -iend {step}"
                    f" >NC{step}-{step}_log.txt"), shell=True, cwd=cwd)

            # prepare instances for next loop
            instances_to_run = np.delete(instances_to_run, idx, 0)
            if starting_order == "alternate":
                instances_to_run = np.flip(instances_to_run, axis=0)

            started_instance = True
            time.sleep(wait_time)  # long wait to initialize RAM
        #  END of per project loop --------------------------------------------
    #  END of all projects loop -----------------------------------------------

    # wait for last NumCalc instances to finish
    current_time = time.strftime("%b %d %Y, %H:%M:%S", time.localtime())
    message = (f"\n... waiting for the last NumCalc instances to finish "
               f"(checking every second, {current_time})")
    _print_message(message, text_color_reset, log_file)
    while True:

        if _numcalc_instances() == 0:
            break

        time.sleep(wait_time_busy)

    # Check all projects that may need to be executed -------------------------
    current_time = time.strftime("%b %d %Y, %H:%M:%S", time.localtime())

    message = ("\nThe following instances did not finish\n"
               "--------------------------------------\n")

    for project in all_projects:
        all_instances, instances_to_run, *_ = _check_project(
            project, numcalc_executable, log_file)

        if instances_to_run is None:
            continue

        if instances_to_run.shape[0] > 0:
            message += f"{os.path.basename(project)}: "
            unfinished = [f"source {int(p[0])} step {int(p[1])}"
                          for p in instances_to_run]
            message += "; ".join(unfinished) + "\n"

    if message.count("\n") > 3:
        message += f"Finished at {current_time}"
        _raise_error(message, text_color_reset, log_file, confirm_errors)
    else:
        message = f"\nAll NumCalc projects finished at {current_time}"
        _print_message(message, text_color_reset, log_file)

        if confirm_errors:
            input(text_color_green + 'DONE. Hit Enter to exit')
            print(text_color_reset)


def _raise_error(message, text_color, log_file, confirm_errors):
    """Two different ways of error handling depending on `confirm_errors`"""

    # error to logfile
    with open(log_file, "a", encoding="utf8", newline="\n") as f:
        f.write("\n\n" + message + "\n")

    # error to console
    if confirm_errors:
        if os.name == 'nt':  # Windows detected
            print(message)
            input("Press Enter to exit manage_numcalc")
        else:  # elif os.name == 'posix': Linux or Mac detected
            print(text_color + message)
            input(text_color + "Press Enter to exit manage_numcalc\033[0m")
        raise Exception("manage_numcalc was stopped due to an error")
    else:
        raise ValueError(message)


def _print_message(message, text_color, log_file):
    """Print message to console and log file"""

    if os.name == 'nt':  # Windows detected
        text_color = ''  # color codes do not work as intended on Win10
    print(text_color + message)

    with open(log_file, "a", encoding="utf8", newline="\n") as f:
        f.write(message + "\n")


def _get_current_ram(ram_offset):
    """Get the available RAM = free RAM - ram_offset"""
    ram_info = psutil.virtual_memory()
    ram_available = max([0, ram_info.available / 1073741824 - ram_offset])
    ram_used = ram_info.used / 1073741824
    return ram_available, ram_used


def _numcalc_instances():
    """Return the number of currently running NumCalc instances"""

    numcalc_executable = 'NumCalc' if os.name != 'nt' else 'NumCalc.exe'

    num_instances = 0
    for p in psutil.process_iter(['name', 'memory_info']):
        if p.info['name'].endswith(numcalc_executable):
            num_instances += 1

    return num_instances


def _check_project(project, numcalc_executable, log_file):
    """
    Find unfinished instances (frequency steps) in a Mesh2HRTF project folder

    Parameters
    ----------
    project : str
        Full path of the Mesh2HRTF project folder

    Returns
    -------
    all_instances : numpy array
        Array of shape (N, 4) where N is the number of detected frequency
        steps in all source_* folders in the project. The first column contains
        the source number, the second the frequency step, the third the
        frequency in Hz, and the fourth the estimated RAM consumption in GB.
    instances_to_run : numpy array, None
        Array of size (M, 4) if any instances need to be run (in this case M
        gives the unfinished instances). ``None``, if all instances are
        finished.
    source_counter : int
        Number of sources in the project
    """

    # get source folders and number of sources
    sources = glob.glob(os.path.join(project, 'NumCalc', "source_*"))
    source_counter = len(sources)
    sources = [os.path.join(project, 'NumCalc', f"source_{s+1}")
               for s in range(source_counter)]

    # loop source_* folders
    for source_id, ff in enumerate(sources):

        # estimate RAM consumption if required
        if not os.path.isfile(os.path.join(ff, "Memory.txt")):

            _print_message(f"Obtaining RAM estimates for {ff}",
                           '\033[0m', log_file)

            if os.name == 'nt':  # Windows detected
                # run NumCalc and route all printouts to a log file
                subprocess.run(
                    f"{numcalc_executable} -estimate_ram",
                    stdout=subprocess.DEVNULL, cwd=ff, check=True)

            else:  # elif os.name == 'posix': Linux or Mac detected
                # run NumCalc and route all printouts to a log file
                subprocess.run(
                    [f"{numcalc_executable} -estimate_ram"],
                    shell=True, stdout=subprocess.DEVNULL, cwd=ff, check=True)

        # get RAM estimates and prepend source number
        estimates = m2s.read_ram_estimates(ff)
        estimates = np.concatenate(
            ((source_id + 1) * np.ones((estimates.shape[0], 1)), estimates),
            axis=1)

        if source_id == 0:
            all_instances = estimates
            instances_to_run = None
        else:
            all_instances = np.append(all_instances, estimates, axis=0)

        # loop frequency steps
        for step in range(estimates.shape[0]):

            if not os.path.isfile(os.path.join(
                    ff, "be.out", f"be.{1 + step}", "pEvalGrid")):

                # there are no output files, process this
                if instances_to_run is None:
                    instances_to_run = np.atleast_2d(estimates[step])
                else:
                    instances_to_run = np.append(
                        instances_to_run, np.atleast_2d(estimates[step]),
                        axis=0)

            elif os.path.isfile(os.path.join(
                    ff, f'NC{1 + step}-{1 + step}.out')):

                # check if "NCx-x.out" contains "End time:" to confirm that
                # the simulation was completed.
                nc_out = os.path.join(
                    ff, f'NC{1 + step}-{1 + step}.out')
                with open(nc_out, "r", encoding="utf8", newline="\n") as f:
                    nc_out = "".join(f.readlines())

                if 'End time:' not in nc_out:
                    # instance did not finish
                    if instances_to_run is None:
                        instances_to_run = np.atleast_2d(estimates[step])
                    else:
                        instances_to_run = np.append(
                            instances_to_run, np.atleast_2d(estimates[step]),
                            axis=0)

    return all_instances, instances_to_run, source_counter
