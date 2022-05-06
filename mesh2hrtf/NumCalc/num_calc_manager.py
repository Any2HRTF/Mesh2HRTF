"""
Run NumCalc on one or multiple Mesh2HRTF project folders.

This script monitors the RAM and CPU usage and starts a new NumCalc instance
whenever enough resources are available. It detects interrupted and incomplete
results and therefore can be interrupted and started again without losing
major progress.

run ``python num_calc_manager -h`` for help.

HOW TO USE
----------
 1- copy this .py file into the project folder OR the folder that contains multiple projects that you will simulate.
 2- Windows only: Make sure that NumCalc runtime files are located in a folder next to "start_path" or next to this
     "NumCalcManager.py" file (that is "NumCalc_WindowsExe" folder with compiled NumCalc.exe, .dll and .bat files.)
 3- "open/run" this "NumCalcManager.py" file with Python.      (Python3 and "psutil" must be installed)
     By default script automatically runs projects next to the script BUT you can specify custom "start_path".
 4- The script will manage all simulations and print out the progress until it is finished
 5- It is highly recommended to check NC _log.txt files of the first completed instances for issues before letting
     the simulation continue till the end.   DONE.

Extra tips:
  1- It is not recommended to run multiple NumCalcManager.py scripts on the same computer (but it would work)
  2- avoid folder names with spaces. Code is not tested for that.
  3- If there is a problem with some instance result - just delete its output folder "be.X" and the Manager
          will automatically re-run that instance on the next run.
  + It is recommended to re-launch NumCalcManager after the simulation is finished. This way it will quickly
      re-check the status and either confirm 100% ready or attempt to re-launch some instances that did not complete.
"""

# made for Python 3.7+
#     v1.00    2021-09-18     First versions. Tested on Windows 10.
#                             (Sergejs Dombrovskis)
#     v1.01    2021-09-21     Small fixes & improvements. (by S.D.)
#     v1.02    2021-09-23     More improvements. (by S.D.)
#     v2.00    2021-10-10     Almost-rewrite with major improvements including
#                             support for running multiple projects parsing
#                             text files and sorting of instances by actual
#                             frequency (by S.D.)
#     v2.01    2021-11-16     Critical bugfix. Minor improvements. Tested some
#                             more. (by S.D.)
#     v3.00    2022-01-20     Adapted code to the new Mesh2HRTF project
#                             structure. (by S.D.)
#     v3.05    2022-01-22     Added optimization to launch new instances faster
#                             (can save 15+ min of processing time when running
#                             very large number of instances simultaneously,
#                             by S.D.)
#     v3.10    2022-02-02     LINUX & MAC support added, note NumCalc must be
#                             already compiled (by S.D.)
#     v3.11    2022-02-05     minor fix for cases when re-running 100% complete
#                             projects (by S.D.)
#     v4.00    2022-05-05     added parser for command line arguments (by J.T.)
#     v4.01    2022-05-06     - flake8 and code styling and styling
#                             - unify parameters auto_set_max_instances and
#                               max_instances

# ToDo: future improvement ideas
#   - figure out a good method to split simulation projects over multiple
#     computers (so far it is only possible to split projects by
#     moving/deleting "../NumCalc/step_x/" folders - assuming there are more
#     than the usual one).
#   - add start_path input to the script (so-far start_path can only be changed
#     by editing the script itself)
#   - add automatic compilation of NumCalc on Linux/Mac (to automate-away one
#     more hassle)
#   - add nice printout for total processing time.
#   - Long-Term - launch "Output2HRTF" Octave/Matlab/Python processing at the
#     end of each simulation.
#   - Long-Term - try to run some low frequency instances in parallel with
#     high-frequency instances to better utilize RAM in the beginning of the
#     simulation project (this could be tricky and require new RAM monitoring
#     code).

import os
import shutil
import time
import psutil
import subprocess
import argparse

# TODO:
# - Format DocString
# - Add argument to pass NumCalc binary (default "NumCalc")
# - add argument for `confirm``

# Done
# - removed else case in check_project
# - removed `if len(all_projects) > 1: ... else ...`
# - removed double check for existing results (variable `pathToCheck`)


# helping functions -----------------------------------------------------------
def raise_error(message, text_color, confirm_errors):
    """Two different ways of error handling depending on `confirm_errors`"""
    if confirm_errors:
        print(text_color + message)
        input(text_color + "Press Enter to exit num_calc_manager")
        raise Exception("num_calc_manager was stopped due to an error")
    else:
        raise ValueError(message)


def get_num_calc_processes():
    """Return a list with the pid, names, and bytes of each NumCalc process"""
    pid_names_bytes = [
            (p.pid, p.info['name'], p.info['memory_info'].rss)
            for p in psutil.process_iter(['name', 'memory_info'])
            if p.info['name'] == NumCalc_filename]
    return pid_names_bytes


def check_project(project):
    """
    Find unfinished instances (frequency steps) in a Mesh2HRTF project folder

    Parameters
    ----------
    project : str
        Full path of the Mesh2HRTF project folder

    Returns
    -------
    _type_
        _description_
    """

    # initialize variables
    all_instances = []           # all folder for each instance
    instances_to_run = []        # folders for each instance that must be run
    source_counter: int = 0      # but up to sources 99999 are supported
    frequency_steps_nr: int = 0  # init
    frequency_step = 0           # init
    min_frequency = 0            # init

    project_numcalc = os.path.join(project, 'NumCalc')

    # parse "Info.txt" to get info about frequencies and instances
    with open(os.path.join(project, "Info.txt"), "r", encoding="utf8",
              newline="\n") as f:
        for line in f:
            if line.find('Minimum evaluated Frequency') != -1:
                idl = line.find(':')
                min_frequency = float(line[idl + 2:-1])
            elif line.find('Frequency Stepsize') != -1:
                idl = line.find(':')
                frequency_step = float(line[idl + 2:-1])
            elif line.find('Frequency Steps') != -1:
                idl = line.find(':')
                # NOTE: total instances = frequency_steps_nr * source_counter!!
                frequency_steps_nr = int(line[idl + 2:-1])

        del line, idl  # remove no longer needed variables

    # loop source_* folders
    for ff in os.listdir(project_numcalc):

        if not ff.startswith("source_"):
            continue

        source_counter += 1              # count this source
        source_id = int(ff[7:])
        base_nr_str = f'{source_id:05}'  # source_id as fixed length string

        # loop frequency steps
        for step in range(frequency_steps_nr):  # counting from zero

            # update list of all instances
            all_instances.append(
                'source' + base_nr_str + '_step' + f'{1 + step:05}')

            if not os.path.isfile(os.path.join(
                    project_numcalc, ff, "be.out", f"be.{1 + step}",
                    "pEvalGrid")):
                # there are no output files, process this
                instances_to_run.append(all_instances[-1])

            elif os.path.isfile(os.path.join(
                    project_numcalc, ff, f'NC{1 + step}-{1 + step}.out')):

                # check if "NCx-x.out" contains "End time:" to confirm that
                # the simulation was completed.
                nc_out = os.path.join(
                    project_numcalc, ff, f'NC{1 + step}-{1 + step}.out')
                with open(nc_out, "r", encoding="utf8", newline="\n") as f:
                    nc_out = "".join(f.readlines())

                if 'End time:' not in nc_out:
                    # instance did not finish
                    instances_to_run.append(all_instances[-1])

    return all_instances, instances_to_run, min_frequency, frequency_step, \
        frequency_steps_nr, source_counter


# parse command line input ----------------------------------------------------
parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument(
    "--path", default=False, type=str,
    help=("The working directory. This can be a directory that contains "
          "multiple Mesh2HRTF project folders, a Mesh2HRTF project folder or "
          "a NumCalc folder inside a Mesh2HRTF project folder. The default "
          "uses the current working directory"))
parser.add_argument(
    "--wait_time", default=15, type=int,
    help=("Delay in seconds for waiting until the RAM and CPU usage is checked"
          " after a NumCalc instance was started."))
parser.add_argument(
    "--ram_safety_factor", default=0.95, type=float,
    help=("Detect RAM usage of highest frequency and wait until free RAM "
          "exceeds `safety factor` times this value before starting the next "
          "NumCalc instance."))
parser.add_argument(
    "--cleanup_after_finish", default=True, type=bool,
    help="Delete NumCalc executables after completion")
parser.add_argument(
    "--max_cpu_load", default=80, type=int,
    help="Maximum allowed CPU load in percent")
parser.add_argument(
    "--max_instances", default=0, type=int,
    help=("The maximum numbers of parallel NumCalc instances. If this is 0 "
          "(default) the maximum number is estimated based on the CPU usage "
          "of the instance calculating the highest frequency and "
          "`max_cpu_load` (see above)"))
parser.add_argument(
    "--confirm_errors", default=True, type=bool,
    help=("If True, num_calc_manager waits for user input in case an error "
          "occurs."))

args = vars(parser.parse_args())

start_path = args["path"] if args["path"] \
    else os.path.dirname(os.path.realpath(__file__))
seconds_to_initialize = args["wait_time"]
ram_safety_factor = args["ram_safety_factor"]
clean_up_after_finish = args["cleanup_after_finish"]
max_cpu_load_percent = args["max_cpu_load"]
max_instances = args["max_instances"]
confirm_errors = args["confirm_errors"]

# initialization --------------------------------------------------------------
if os.name == 'nt':  # Windows detected
    NumCalc_filename = "NumCalc.exe"
    # folder where NumCalc_runtime_files are found
    NumCalc_runtime_folder = 'NumCalc_WindowsExe'
    # files that are needed to execute NumCalc
    NumCalc_runtime_files = ['NumCalc.exe', 'libgcc_s_seh-1.dll',
                             'libstdc++-6.dll', 'libwinpthread-1.dll']
else:  # elif os.name == 'posix': # Linux or Mac detected
    NumCalc_filename = "NumCalc"

# trick to get colored print-outs   https://stackoverflow.com/a/54955094
os.system("")
text_color_red = '\033[31m'
text_color_green = '\033[32m'
text_color_cyan = '\033[36m'
text_color_reset = '\033[0m'

# Detect what the start_path or "getcwd()" is pointing to:
print('NumCalcManager started with start_path = "' + start_path + '"')
if os.path.basename(start_path) == 'NumCalc':
    # start_path is a NumCalc folder
    all_projects = [os.path.dirname(start_path)]
elif os.path.isfile(os.path.join(start_path, 'Info.txt')):
    # start_path is a Mesh2HRTF project folder
    all_projects = [start_path]
else:
    # start_path contains multiple Mesh2HRTF project folders
    all_projects = []  # list of project folders to execute
    for subdir in os.listdir(start_path):
        if os.path.isfile(os.path.join(start_path, subdir, 'Info.txt')):
            all_projects.append(os.path.join(start_path, subdir))

# stop if no project folders were detected
if len(all_projects) == 0:
    message = ("num_calc_manager could not detect any Mesh2HRTF projects at "
               f"start_path={start_path}")
    raise_error(message, text_color_red, confirm_errors)

# Find 'NumCalc_WindowsExe' folder if we are on Windows
if os.name == 'nt':  # Windows detected
    if os.path.isfile(os.path.join(all_projects[0], NumCalc_runtime_folder,
                                   NumCalc_runtime_files[0])):
        # located inside the project folder
        NumCalc_runtime_location = os.path.join(
            all_projects[0], NumCalc_runtime_folder)
    elif os.path.isfile(os.path.join(os.path.dirname(all_projects[0]),
                                     NumCalc_runtime_folder,
                                     NumCalc_runtime_files[0])):
        # located is inside the folder that contains all Mesh2HRTF projects
        NumCalc_runtime_location = os.path.join(
            os.path.dirname(all_projects[0]), NumCalc_runtime_folder)
    elif os.path.isfile(os.path.join(all_projects[0],
                                     NumCalc_runtime_files[0])):
        # located directly in the project folder.
        NumCalc_runtime_location = os.path.join(all_projects[0])
    else:
        NumCalc_runtime_location = ""

    # Check that each required runtime file is present:
    for calc_file in NumCalc_runtime_files:
        if not os.path.isfile(os.path.join(NumCalc_runtime_location,
                                           calc_file)):
            message = (
                f"The file {calc_file} is missing in {NumCalc_runtime_folder}"
                "or num_calc_manager could not find the folder "
                f"{NumCalc_runtime_folder} that is usually located one level "
                "above the Mesh2HRTF project folder")
            raise_error(message, text_color_red, confirm_errors)

# else:  # elif os.name == 'posix': # Linux or Mac detected
    #  ToDo add check that NumCalc is compiled and working in Linux

# Check all projects that may need to be executed
projects_to_run = []

for project in all_projects:
    all_instances, instances_to_run, min_frequency, frequency_step, \
        frequency_steps_nr, source_counter = check_project(project)
    if len(instances_to_run) > 0:
        projects_to_run.append(project)
        print((f'Project "{os.path.basename(project)}" has '
               f'{len(instances_to_run)} out of {len(all_instances)} instances'
               ' to run'))
    else:
        print((f'Project "{os.path.basename(project)}" is '
               'already complete'))

del all_projects

# loop to process all projects
for pp, project in enumerate(projects_to_run):

    # Check how many instances are in this Project:
    root_NumCalc = os.path.join(project, 'NumCalc')
    all_instances, instances_to_run, min_frequency, frequency_step, \
        frequency_steps_nr, source_counter = check_project(project)
    total_nr_to_run = len(instances_to_run)

    # Status printouts:
    print(text_color_reset + "\n\n")
    print(f'{text_color_cyan}Started Project {pp + 1}/{len(projects_to_run)}')
    print(text_color_reset + "\n\n")
    print(f"--- {len(all_instances)} frequency steps contained in the project")
    if total_nr_to_run == 0:
        print("--- All NumCalc simulations in this project are complete")
        continue
    print(f"--- {total_nr_to_run} NumCalc simulations will be run")

    # ADVANCED sorting
    # (Build matching list of frequencies for each "instances_to_run")
    matched_freq_of_inst = [0] * len(instances_to_run)
    for inst in range(total_nr_to_run):
        idx = int(instances_to_run[inst][16:])
        # frequency in Hz
        matched_freq_of_inst[inst] = int(
            min_frequency + frequency_step * (idx - 1))

    # Sort list to run the largest frequencies that consume the most RAM first
    # (needed for Both ears!!!)
    Sorting_List = sorted(zip(matched_freq_of_inst, instances_to_run),
                          reverse=True)
    # sort "instances_to_run" according to decreasing frequency
    instances_to_run = [x for _, x in Sorting_List]
    # update Matched list to correspond to "instances_to_run"
    matched_freq_of_inst = [y for y, _ in Sorting_List]
    del Sorting_List, idx

    # main loop for each instance
    start_time = time.localtime()
    print(time.strftime("%d %b - %H:%M:%S", start_time))

    for NC_ins in range(total_nr_to_run):
        # current source and frequency step
        source = int(instances_to_run[NC_ins][6:11])
        step = int(instances_to_run[NC_ins][16:])

        print((f"--- {NC_ins + 1}/{total_nr_to_run} preparing >>> "
               f"{matched_freq_of_inst[NC_ins]}Hz <<< instance from "
               f"source_{source}, step{str(step)}"))

        if os.name == 'nt':  # Windows detected
            # copy the NumCalc.exe files into the source_ folder (if necessary)
            for calc_file in NumCalc_runtime_files:
                if not os.path.isfile(os.path.join(
                        root_NumCalc, f"source_{source}", calc_file)):
                    shutil.copyfile(
                        os.path.join(NumCalc_runtime_location, calc_file),
                        os.path.join(root_NumCalc, f"source_{source}",
                                     calc_file))

        # Check the RAM & run instance if feasible
        RAM_info = psutil.virtual_memory()
        print((f"{round((RAM_info.available / 1073741824), 2)} GB free RAM, "
               f"{RAM_info.percent}% used ("
               f"{time.strftime('%d %b - %H:%M:%S', time.localtime())})"))

        # use this to autodetect how many instances can at most be executed
        if NC_ins > 0 and not max_instances:
            # noinspection PyBroadException
            try:
                # it is better to get fresh pid (hopefully at least one NumCalc
                # process is still running)
                pid_names_bytes = get_num_calc_processes()

                PrcInfo = psutil.Process(pid_names_bytes[0][0])

                # wait until CPU usage is realistic
                Instance_CPU_usageNow = 0
                while Instance_CPU_usageNow < 0.001:
                    Instance_CPU_usageNow = \
                        PrcInfo.cpu_percent() / psutil.cpu_count()
                    time.sleep(0.01)

                # calculate optimal maximum number of NumCalc processes
                max_instances = round(
                    max_cpu_load_percent / Instance_CPU_usageNow)

                print(("One NumCalc instance requires "
                       f"{round(Instance_CPU_usageNow, 1)}% of the CPU."
                       "max_instances is now automatically set to "
                       f"{max_instances}"))

            except BaseException:
                message = (
                    "Failed to automatically set the maximum number of "
                    "parallel NumCalc instances. This can happen if a NumCalc "
                    "process finished very fast. Try to lower wait_time or "
                    "manually set max_instances")
                raise_error(message, text_color_red, confirm_errors)

        #  Main checks before launching the next instance
        # (to avoid system resource overload)
        wait_for_resources = False if NC_ins == 0 else True

        while wait_for_resources:
            # Find all NumCalc Processes
            pid_names_bytes = get_num_calc_processes()

            # DEBUGGING --- Find Processes consuming more than 250MB of memory:
            # pid_names_bytes = [
            #     (p.pid, p.info['name'], p.info['memory_info'].rss)
            #     for p in psutil.process_iter(['name', 'memory_info'])
            #     if p.info['memory_info'].rss > 250 * 1048576]

            # start NumCalc instance if none is running
            if len(pid_names_bytes) == 0:
                break

            # if the maximum number of instances to launch is not exceeded
            elif len(pid_names_bytes) < max_instances:

                # find out how much RAM is consumed by any NumCalc Instance
                Max_NumCalc_RAM = pid_names_bytes[0][2]
                if len(pid_names_bytes) > 1:
                    for prcNr in range(1, len(pid_names_bytes)):
                        if pid_names_bytes[prcNr][2] > Max_NumCalc_RAM:
                            Max_NumCalc_RAM = pid_names_bytes[prcNr][2]

                # check if we can run more:
                # IF free RAM is greater than RAM consumption of the biggest
                # NumCalc instance x ram_safety_factor
                RAM_info = psutil.virtual_memory()
                current_time = \
                    time.strftime('%d %b - %H:%M:%S', time.localtime())
                if RAM_info.available > Max_NumCalc_RAM * ram_safety_factor:
                    print(("   Enough RAM to run one more: "
                           f"{round((RAM_info.available / 1073741824), 1)} GB "
                           f"free [{current_time}]"))
                    break

                else:
                    print("   Waiting for more free RAM:   "
                          f"{round((RAM_info.available / 1073741824), 1)} GB "
                          f"free "
                          f"({round((Max_NumCalc_RAM * ram_safety_factor / 1073741824), 1)} "  # noqa
                          f"GB needed) [{current_time}]")

            else:
                print(("   No more instances allowed. Waiting for 1/"
                       f"{max_instances} instances to finish"))

            # delay before trying the while loop again
            time.sleep(seconds_to_initialize)

            # END of the while loop -------------------------------------------

        # start next instance -------------------------------------------------
        print((f"- {NC_ins + 1}/{total_nr_to_run} STARTING instance from: "
               f"{os.path.basename(project)} >>>    source_{source}, step "
               f"{step}"))

        # change working directory
        os.chdir(os.path.join(root_NumCalc, "source_" + str(source)))

        if os.name == 'nt':  # Windows detected
            # create a log file for all print-outs
            LogFileHandle = open(f"NC{step}-{step}_log.txt", "w")
            # run NumCalc and route all printouts to a log file
            subprocess.Popen(f"{NumCalc_filename} -istart {step} -iend {step}",
                             stdout=LogFileHandle)

        else:  # elif os.name == 'posix': Linux or Mac detected
            # run NumCalc and route all printouts to a log file
            subprocess.Popen((f"{NumCalc_filename} -istart {step} -iend {step}"
                              f" >NC{step}-{step}_log.txt"), shell=True)

        # optimize waiting time (important if available RAM >>> needed RAM)
        if NC_ins > 0:
            # noinspection PyUnboundLocalVariable
            if RAM_info.available > Max_NumCalc_RAM * 3:
                # wait less, if RAM is enough for three new instances
                waitTime = 0.5
            elif RAM_info.available > Max_NumCalc_RAM * 2:
                # wait less if RAM is enough for two new instances
                waitTime = seconds_to_initialize / 2
            else:
                # wait longer to assess how much RAM is available
                waitTime = seconds_to_initialize
            print((f"   ... waiting {waitTime} s for current instance to "
                   "initialize RAM"))

        else:
            # always wait for the 1st instance to initialize to get worst-case
            # RAM use estimate
            waitTime = seconds_to_initialize
            print((f"   ... waiting {seconds_to_initialize} s for the first "
                   "instance to initialize RAM"))

        # Wait for instance to initialize before attempting to start the next
        time.sleep(waitTime)

        # Check if all NumCalc processes crashed: Find all NumCalc Processes
        pid_names_bytes = get_num_calc_processes()

        if len(pid_names_bytes) == 0:
            message = (("No NumCalc processes running. Likely the last "
                        "launched instance crashed. Read NC.out log inside "
                        f"{instances_to_run[NC_ins]}"))
            raise_error(message, text_color_green, confirm_errors)

    #  END of the main project loop -------------------------------------------

    while True:
        # Find all NumCalc Processes
        pid_names_bytes = get_num_calc_processes()

        # no NumCalc processes are running, so Finish
        if len(pid_names_bytes) == 0:
            break

        print((f"... waiting {2 * seconds_to_initialize} s for last processes "
               "to finish"))
        time.sleep(2 * seconds_to_initialize)

    # Windows only: Delete all NumCalc files from individual source_ folders
    if os.name == 'nt' and clean_up_after_finish:
        print("Deleting local copies of NumCalc executables")
        for subdir in os.listdir(root_NumCalc):
            for calc_file in NumCalc_runtime_files:
                if os.path.isfile(os.path.join(
                        root_NumCalc, subdir, calc_file)):
                    os.remove(os.path.join(root_NumCalc, subdir, calc_file))

#  END of all_projects loop.

print(("All NumCalc projects finished. ["
       f"{time.strftime('%d %b - %H:%M:%S', time.localtime())}]"))
if confirm_errors:
    input(text_color_green + 'DONE. Hit Enter to exit')
