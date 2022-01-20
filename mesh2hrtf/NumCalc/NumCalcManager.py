import os
import shutil
import time
import psutil
import subprocess
#       ---  "psutil" is required !  https://github.com/giampaolo/psutil/blob/master/INSTALL.rst  ---
#

# -------------------------------------------  User Settings:  --------------------------------------------

StartPath = os.path.dirname(os.path.realpath(__file__))  # Automatically detect the working directory
#                            (usually in the main project folder (where "Info.txt" file is)
#                            OR NumCalcManager can start from a folder that contains multiple Mesh2HRTF projects to run)
# StartPath = r"C:\hrtf_NumCalc"  # alternatively specify which folder to execute (use r"..." raw strings!)
# StartPath = r"C:\hrtf_NumCalc\My_project"  # alternatively specify which folder to execute (use r"..." raw strings!)

SecondsToInitialize = 15  # delay in seconds during which new NumCalc instance should initialize its full RAM usage.
RAM_safetyFactor = 0.95  # extra instances are launched when:
#                           free RAM is greater than RAM consumption of the biggest NumCalc instance * RAM_safetyFactor.
CleanUp_afterFinish = True  # Delete all NumCalc executable files from individual instance folders after completion?
Auto_set_MaxInstances = True  # use this to autodetect how many instances can be at most executed.
Max_CPU_load_percent = 80  # target % for how much CPU can be used by Mesh2HRTF, assuming that RAM allows for so many
#                              instances. Recommended to NOT set more than 80% because the real CPU usage can be higher.
Max_Instances = 8  # default instance limit is used in case "Auto_set_MaxInstances = False"
# ----------------------------------------  End of User Settings  -----------------------------------------

#  "NumCalcManager.py" executable script automatically runs one or multiple Mesh2HRTF simulation projects while
#  maximizing the computer resource utilization (max RAM, max CPU without overloading).
#  It automatically detects interrupted/incomplete simulation and therefore NumCalc simulation can be interrupted and
#  started again without any user input of major loss of progress. It is recommended
#
#                               HOW TO USE:
#
#   1- copy this .py file into the project folder OR the folder that contains multiple projects that you will simulate.
#   2- Make sure that NumCalc runtime files are located in a folder next to "StartPath" or this "NumCalcManager.py" file
#       (on Windows that would be "NumCalc_WindowsExe" folder with compiled NumCalc.exe, .dll and .bat files.
#        On Linux/Mac it is enough to have NumCalc source files and the manager can compile the runtime)
#   3- "open/run" this "NumCalcManager.py" file with Python.      (Python3 and "psutil" must be installed)
#       By default script automatically runs projects next to the script BUT you can specify custom "StartPath".
#   4- The script will manage all simulations and print out the progress until it is finished
#   5- It is highly recommended to check NC _log.txt files of the first completed instances for issues before letting
#       the simulation continue till the end.   DONE.
#
# Extra tips:
#   1- It is not recommended to run multiple NumCalcManager.py scripts on the same computer (but it would work)
#   2- avoid folder names with spaces. Code is not tested for that.
#   3- If there is a problem with some instance result - just delete its output folder "be.X" and the Manager
#           will automatically re-run that instance on the next run.
#   + It is recommended to re-launch NumCalcManager after the simulation is finished. This way it will quickly
#       re-check the status and either confirm 100% ready or attempt to re-launch some instances that did not complete.
#
#
#   made for Python 3.7+                                            see license details at the end of the script.
#       v1.00    2021-09-18     First versions. Tested on Windows 10. (by Sergejs Dombrovskis)
#       v1.01    2021-09-21     Small fixes & improvements. (by S.D.)
#       v1.02    2021-09-23     More improvements. (by S.D.)
#       v2.00    2021-10-10     Almost-rewrite with major improvements including support for running multiple projects
#                                   parsing text files and sorting of instances by actual frequency (by S.D.)
#       v2.01    2021-11-16     Critical bugfix. Minor improvements. Tested some more. (by S.D.)
#       v3.00    2022-01-20     Adapted code to the new Mesh2HRTF project structure. (by S.D.)
#
# The logic is:
#   NumCalcManager checks RAM usage and launches more than 1 NumCalc instance IF it detects enough free RAM in the
#   system. It keeps launching new instances as resources free up. In addition, it makes sure to not overload the CPU
#   in case there is more than enough RAM. (all of this resource monitoring works reasonably well)
#   To have full control over RAM usage NumCalcManager runs NumCalc instances with a single frequency step at a time.
#
#


# ToDo: future improvement ideas
#   - figure out a good method to split simulation projects over multiple computers (so far it is only possible to
#        split projects by moving/deleting "../NumCalc/step_x/" folders - assuming there are more than the usual one).
#   - add proper support for Linux & Mac (changes are needed mostly regarding launching of NumCalc executable).
#   - fix Known minor bug - “The system can not find the file specified.” printouts - It does not cause any issues.
#   - add StartPath input to the script (so-far StartPath can only be changed by editing the script itself)
#   - add automatic compilation of NumCalc on Linux/Mac (to automate-away one more hassle)
#   - add nice printout for total processing time.
#   - Long-Term - launch "Output2HRTF" Octave/Matlab/Python processing at the end of each simulation.
#   - Long-Term - try to run some low frequency instances in parallel with high-frequency instances to better utilize
#       RAM in the beginning of the simulation project (this could be tricky and require new RAM monitoring code).
#


# init various stuff:

if os.name == 'nt':  # Windows detected
    NumCalc_filename = "NumCalc.exe"
    NumCalcRuntime_folder = 'NumCalc_WindowsExe'  # folder where NumCalcRuntime_files are found
    # list all files that are needed to execute NumCalc (all these files are mandatory):
    NumCalcRuntime_files = ['NumCalc.exe',
                            'libgcc_s_seh-1.dll', 'libstdc++-6.dll', 'libwinpthread-1.dll']
else:  # elif os.name == 'posix': # Linux or Mac detected
    NumCalc_filename = "This script needs to be upgraded to run on Linux or Mac"
    NumCalcRuntime_folder = 'not needed?'  # folder where NumCalcRuntime_files are found
    # list all files that are needed to execute NumCalc (all these files are mandatory):
    NumCalcRuntime_files = ['NumCalc.exe',
                            'libgcc_s_seh-1.dll', 'libstdc++-6.dll', 'libwinpthread-1.dll']
    print("ERROR - Sorry, this NumCalcManager version is not upgraded to run on this OS.")
    input("  -press Enter- to exit NumCalcManager.")
    raise Exception("not ok to continue")

os.system("")  # trick to get colored print-outs   https://stackoverflow.com/a/54955094
Text_Color_RED = '\033[31m'
Text_Color_GREEN = '\033[32m'
Text_Color_CYAN = '\033[36m'
Text_Color_RESET = '\033[0m'

# Detect what the StartPath or "getcwd()" is pointing to:
print('NumCalcManager started with StartPath = "' + StartPath + '"')
if os.path.basename(StartPath) == 'NumCalc':  # - is it a NumCalc Folder?
    all_projects = [os.path.dirname(StartPath)]  # correct project folder
elif os.path.isfile(os.path.join(StartPath, 'Info.txt')):  # - is it Project folder?
    all_projects = [StartPath]  # correct project folder
else:  # - is it a folder with multiple projects?
    all_projects = []  # list of project folders to execute
    for subdir in os.listdir(StartPath):
        if os.path.isfile(os.path.join(StartPath, subdir, 'Info.txt')):  # - is it Project folder?
            all_projects.append(os.path.join(StartPath, subdir))  # generate list of project folders to execute
    if len(all_projects) == 0:  # not good
        print(Text_Color_RED + 'ERROR - StartPath = "' + StartPath + '" does not contain any Mesh2HRTF projects to run')
        print(Text_Color_RED + ' Please make sure that "StartPath" variable is set correctly')
        input(Text_Color_RED + "  -press Enter- to exit NumCalcManager.")
        raise Exception("not ok to continue")

# Find NumCalc.exe and related files:
if os.path.isfile(os.path.join(all_projects[0], NumCalcRuntime_folder, NumCalcRuntime_files[0])):
    # basic location where 'NumCalc_WindowsExe' folder is inside the project folder.
    NumCalcRuntime_Location = os.path.join(all_projects[0], NumCalcRuntime_folder)
elif os.path.isfile(os.path.join(os.path.dirname(all_projects[0]), NumCalcRuntime_folder, NumCalcRuntime_files[0])):
    # recommended location in the folder that contains all Mesh2HRTF projects
    NumCalcRuntime_Location = os.path.join(os.path.dirname(all_projects[0]), NumCalcRuntime_folder)
elif os.path.isfile(os.path.join(all_projects[0], NumCalcRuntime_files[0])):
    # fall-back location of NumCalcRuntime_files found directly in the project folder.
    NumCalcRuntime_Location = os.path.join(all_projects[0])
else:
    NumCalcRuntime_Location = 'Failed to locate the "' + NumCalcRuntime_folder + \
                              '" (usually it is next to the Mesh2HRTF project folder - one level above)'
# Check that each required runtime file is present:
for CalcFile in NumCalcRuntime_files:
    if not os.path.isfile(os.path.join(NumCalcRuntime_Location, CalcFile)):
        print(Text_Color_RED + "ERROR - Required file " + CalcFile + " from compiled NumCalc is not found in folder:")
        print(Text_Color_RED + "         " + NumCalcRuntime_Location)
        input(Text_Color_RED + "  -press Enter- to exit NumCalcManager.")
        raise Exception("not ok to continue")


def check_project(project):  # declare function to find all unfinished instances in a given project folder
    all_instances = []  # list of folder names for each instance
    instances_to_run = []  # list of folder names for each instance that needs to be executed
    source_counter: int = 0  # usually can be 2 sources if both ears simulated at once, but up to 99999 are supported
    frequency_steps_nr: int = 0  # init
    frequency_step = 0  # init
    min_frequency = 0  # init
    project_numcalc = os.path.join(project, 'NumCalc')  # just to make life easier

    # parse "Info.txt" to get info about frequencies and instances
    f = open(os.path.join(project, "Info.txt"), "r", encoding="utf8", newline="\n")
    for line in f:
        if line.find('Minimum evaluated Frequency') != -1:
            idl = line.find(':')
            min_frequency = float(line[idl + 2:-1])
        elif line.find('Frequency Stepsize') != -1:
            idl = line.find(':')
            frequency_step = float(line[idl + 2:-1])
        elif line.find('Frequency Steps') != -1:
            idl = line.find(':')
            frequency_steps_nr = int(line[idl + 2:-1])  # NOTE: total instances = frequency_steps_nr * source_counter!!
        # elif line.find('Highest evaluated Frequency') != -1:
        #     idl = line.find(':')
        #     max_frequency = int(line[idl + 2:-1])
    del line, idl  # remove no longer needed variables
    f.close()

    for subFldr in os.listdir(project_numcalc):
        if subFldr.startswith("source_"):
            source_counter += 1  # count this source
            source_id = int(subFldr[7:])
            base_nr_str = f'{source_id:05}'  # storing source_id into fixed length string

            for step in range(frequency_steps_nr):  # counting from zero
                all_instances.append('source' + base_nr_str + '_step' + f'{1 + step:05}')  # generate list of all inst.

                if not os.path.isfile(os.path.join(project_numcalc, subFldr, "be.out", "be." + str(1 + step),
                                                   "vEvalGrid")):
                    instances_to_run.append(all_instances[-1])  # there are no output files, process this

                elif os.path.isfile(os.path.join(project_numcalc, subFldr,
                                                 "NC" + str(1 + step) + '-' + str(1 + step) + '.out')):
                    # check if "NCx-x.out" contains "End time:" to confirm that simulation was completed.
                    did_not_find_end_time = True  # re-init
                    fil = open(os.path.join(project_numcalc, subFldr, "NC" + str(1 + step) + '-' +
                                            str(1 + step) + '.out'), "r", encoding="utf8", newline="\n")
                    for f_line in fil:
                        if f_line.find('En') != -1:  # maybe performance improvement? :)
                            if f_line.find('End time:') != -1:
                                did_not_find_end_time = False
                    fil.close()

                    if did_not_find_end_time:
                        instances_to_run.append(all_instances[-1])  # this instance did not finish

                else:
                    print(Text_Color_RED + 'EXCEPTION - NumCalcManager will not process project:')
                    print(' "' + project_numcalc + '" because')
                    print(" NumCalc was previously executed outside NumCalcManager in ")
                    print("     " + os.path.join(project_numcalc, subFldr))
                    print(' (clean out "be.xx" folders that were not created by NumCalcManager to continue) ')
                    print(Text_Color_RESET + " ")
                    instances_to_run = []  # reset
                    all_instances = []  # reset
                    break

    return all_instances, instances_to_run, min_frequency, frequency_step, frequency_steps_nr, source_counter


# Check all projects that may need to be executed:
if len(all_projects) > 1:
    Projects_to_run = []  # init boolean
    for proj in range(len(all_projects)):
        All_Instances, Instances_to_Run, Min_Frequency, Frequency_step, Frequency_Steps_Nr, Source_Counter = \
            check_project(all_projects[proj])
        if len(Instances_to_Run) > 0:
            Projects_to_run.append(all_projects[proj])  # mark to run this project
            print('Project "' + os.path.basename(all_projects[proj]) + '" has ' + str(len(Instances_to_Run)) +
                  ' out of ' + str(len(All_Instances)) + ' instances to run')
        else:
            print('Project "' + os.path.basename(all_projects[proj]) + '" is already complete')
else:
    Projects_to_run = all_projects
    # #    if not all(Project_to_run):  # all projects are finished  (no problem)
del all_projects  # just to avoid bugs: removing variables that are no longer relevant

# loop to process all projects
for proj in range(len(Projects_to_run)):

    # Check how many instances are in this Project:
    root_NumCalc = os.path.join(Projects_to_run[proj], 'NumCalc')  # just to make life easier
    All_Instances, Instances_to_Run, Min_Frequency, Frequency_step, Frequency_Steps_Nr, Source_Counter = \
        check_project(Projects_to_run[proj])
    TotalNr_toRun = len(Instances_to_Run)

    # Status printouts:
    if len(Projects_to_run) > 1:
        print(Text_Color_RESET + " ")
        print(Text_Color_RESET + " ")
        print(Text_Color_CYAN + " Started Project", str(proj + 1), " out of ", str(len(Projects_to_run)),
              "detected projects to run")
        print(Text_Color_RESET + " ")
        print(Text_Color_RESET + " ")
    print("--- ", str(len(All_Instances)), "NumCalc Simulations defined in this Mesh2HRTF project")
    print("--- ", str(TotalNr_toRun), "NumCalc Simulations are not yet completed. (Starting from the max frequency)")

    # ADVANCED sorting (Build matching list of frequencies for each "Instances_to_Run")
    Matched_Freq_of_inst = [0] * len(Instances_to_Run)  # init
    for inst in range(0, len(Instances_to_Run)):
        idx = int(Instances_to_Run[inst][16:])
        Matched_Freq_of_inst[inst] = int(
            Min_Frequency + Frequency_step * (idx - 1))  # calculate and fill in frequency [Hz]

    # Sort list to run the largest frequencies that consume the most RAM first (needed for Both ears!!!)
    Sorting_List = sorted(zip(Matched_Freq_of_inst, Instances_to_Run), reverse=True)
    Instances_to_Run = [x for y, x in Sorting_List]  # sort "Instances_to_Run" according to decreasing frequency
    Matched_Freq_of_inst = [y for y, x in Sorting_List]  # update Matched list to correspond to "Instances_to_Run"
    del Sorting_List, idx  # remove no longer needed variables

    # main loop for each instance
    StartTime = time.localtime()
    print(time.strftime("%d %b - %H:%M:%S", StartTime))
    for NC_ins in range(0, TotalNr_toRun):
        Source = int(Instances_to_Run[NC_ins][6:11])    # decode to number
        Step = int(Instances_to_Run[NC_ins][16:])       # decode to number
        print("- ", str(NC_ins + 1), "/", str(TotalNr_toRun), " preparing >>>", str(Matched_Freq_of_inst[NC_ins]),
              "Hz <<< instance from:  source_" + str(Source), " Step", str(Step))

        # double check (roughly) that this instance does not have output
        pathToCheck = os.path.join(root_NumCalc, "source_" + str(Source), "be.out", "be." + str(Step), "vEvalGrid")
        if os.path.isfile(pathToCheck):                                                        
            if 500 < os.path.getsize(pathToCheck):
                print(Text_Color_CYAN, "instance --- source_" + str(Source), " Step", str(Step),
                      "--- already has output data! Skipping")
                print(Text_Color_RESET + " ")
                continue  # jump over this instance

        # copy the NumCalc.exe files into the source_ folder (if necessary)
        if not os.path.isfile(os.path.join(root_NumCalc, "source_" + str(Source), NumCalcRuntime_files[0])):
            for CalcFile in NumCalcRuntime_files:
                if not os.path.isfile(os.path.join(root_NumCalc, "source_" + str(Source), CalcFile)):
                    shutil.copyfile(os.path.join(NumCalcRuntime_Location, CalcFile),
                                    os.path.join(root_NumCalc, "source_" + str(Source), CalcFile))  # copy

        # Check the RAM & run instance if feasible
        RAM_info = psutil.virtual_memory()
        print(str(round((RAM_info.available / 1073741824), 2)), "GB free", "    ---    RAM memory used:",
              str(RAM_info.percent), "%     [", time.strftime("%d %b - %H:%M:%S", time.localtime()), "]")
        # # psutil.cpu_percent()

        # Run this once - normally before launching the 2nd instance IF "Auto_set_MaxInstances == True"
        if NC_ins > 0 and Auto_set_MaxInstances:  # use this to autodetect how many instances can at most be executed.
            # noinspection PyBroadException
            try:
                # it is better to get fresh pid (hopefully at least one NumCalc process is still running)
                Pid_Name_Bytes = [(p.pid, p.info['name'], p.info['memory_info'].rss) for p in
                                  psutil.process_iter(['name', 'memory_info']) if p.info['name'] == NumCalc_filename]
                PrcInfo = psutil.Process(Pid_Name_Bytes[0][0])
                Instance_CPU_usageNow = 0  # init
                while Instance_CPU_usageNow < 0.001:  # wait until CPU usage is realistic
                    Instance_CPU_usageNow = PrcInfo.cpu_percent() / psutil.cpu_count()
                    time.sleep(0.01)  # wait to retry

                # calculate optimal maximum number of NumCalc processes for this system:
                Max_Instances = round(Max_CPU_load_percent / Instance_CPU_usageNow)
                print("One instance loads CPU to", str(round(Instance_CPU_usageNow, 1)), "% on this machine, therefore",
                      "Max_Instances is now automatically set =", str(Max_Instances))
                Auto_set_MaxInstances = False  # mark that max instances does not need to be checked again.
            except BaseException:
                print("!!! Failed to Auto_set_MaxInstances - this can happen if NumCalc process finished very fast -",
                      "you could try to lower your ""SecondsToInitialize"" setting.")

        #
        #  Main checks before launching the next instance (to avoid system resource overload)
        Wait_for_resources = True  # re-init
        if NC_ins == 0:
            Wait_for_resources = False  # always run 1st instance.

        while Wait_for_resources:
            # Find all NumCalc Processes
            Pid_Name_Bytes = [(p.pid, p.info['name'], p.info['memory_info'].rss) for p in
                              psutil.process_iter(['name', 'memory_info']) if p.info['name'] == NumCalc_filename]
            # # FOR DEBUGGING --- Find Processes consuming more than 250MB of memory:
            # # Pid_Name_Bytes = [(p.pid, p.info['name'], p.info['memory_info'].rss) for p in psutil.process_iter([
            #                      'name', 'memory_info']) if p.info['memory_info'].rss > 250 * 1048576]

            if len(Pid_Name_Bytes) == 0:  # if no NumCalc processes are running, so Go start one instance.
                break

            elif len(Pid_Name_Bytes) < Max_Instances:  # if the maximum number of instances to launch is not exceeded

                # find out how much RAM consumes the biggest NumCalc Instance
                if len(Pid_Name_Bytes) == 1:  # only one process
                    Max_NumCalc_RAM = Pid_Name_Bytes[0][2]
                else:
                    Max_NumCalc_RAM = Pid_Name_Bytes[0][2]  # init
                    for prcNr in range(1, len(Pid_Name_Bytes)):
                        if Pid_Name_Bytes[prcNr][2] > Max_NumCalc_RAM:  # finding NumCalc process that consumes most RAM
                            Max_NumCalc_RAM = Pid_Name_Bytes[prcNr][2]

                # check if we can run more:
                # IF free RAM is greater than RAM consumption of the biggest NumCalc instance x RAM_safetyFactor
                RAM_info = psutil.virtual_memory()
                if RAM_info.available > Max_NumCalc_RAM * RAM_safetyFactor:
                    print("   enough RAM to run one more:     ", str(round((RAM_info.available / 1073741824), 1)),
                          "GB free", "     [", time.strftime("%d %b - %H:%M:%S", time.localtime()), "]")
                    break

                else:
                    print("   Waiting for more free RAM:     ", str(round((RAM_info.available / 1073741824), 1)),
                          "GB free    (", str(round((Max_NumCalc_RAM * RAM_safetyFactor / 1073741824), 1)),
                          "GB needed)", "     [", time.strftime("%d %b - %H:%M:%S", time.localtime()), "]")

                    if len(Pid_Name_Bytes) == 1:  # only one process
                        # extra delay before trying the while loop again for very large processes
                        time.sleep(4 * SecondsToInitialize)

            else:
                print("   No more instances allowed - waiting for 1 out of", str(Max_Instances), "instances to finish")

            # delay before trying the while loop again
            time.sleep(SecondsToInitialize)

        # END of Wait_for_resources while loop

        pathToCheck = os.path.join(root_NumCalc, "source_" + str(Source), "NC" + str(Step)+"-"+str(Step) + "_log.txt")

        #
        # START one more instance
        print("- ", str(NC_ins + 1), "/", str(TotalNr_toRun), " STARTING instance from:",
              os.path.basename(Projects_to_run[proj]), ">>>    source_" + str(Source), " Step", str(Step))
        os.chdir(os.path.join(root_NumCalc, "source_" + str(Source)))  # change to the correct directory
        LogFileHandle = open("NC" + str(Step)+"-"+str(Step) + "_log.txt", "w")  # create a log file for all print-outs
        # run NumCalc with "-istart x -iend x" parameters & route all printouts to a log file
        subprocess.Popen(NumCalc_filename + " -istart " + str(Step) + " -iend " + str(Step), stdout=LogFileHandle)
        # # (stderr=LogFileHandle, )
        # #     pid = Popen(["/bin/mycmd", "myarg"]).pid

        # Wait for current instance to initialize before attempting to start the next instance
        print("   ... waiting", str(SecondsToInitialize), "s for current instance to initialize RAM")
        time.sleep(SecondsToInitialize)

        # Check if all NumCalc processes crashed:    Find all NumCalc Processes
        Pid_Name_Bytes = [(p.pid, p.info['name'], p.info['memory_info'].rss) for p in
                          psutil.process_iter(['name', 'memory_info']) if p.info['name'] == NumCalc_filename]
        if len(Pid_Name_Bytes) == 0:
            print(Text_Color_RED +
                  "ERROR - NumCalc processes are NOT running! likely the last launched instance crashed.")
            print(Text_Color_RED + " Read NC.out log inside " + Instances_to_Run[NC_ins] + " folder")
            input(Text_Color_RED + "  -press Enter- to exit NumCalcManager.")
            raise Exception("not ok to continue")

    #  END of the main project loop.

    print("waiting for last processes to finish (in this project)")
    while True:
        # Find all NumCalc Processes
        Pid_Name_Bytes = [(p.pid, p.info['name'], p.info['memory_info'].rss) for p in
                          psutil.process_iter(['name', 'memory_info']) if p.info['name'] == NumCalc_filename]

        if len(Pid_Name_Bytes) == 0:
            break  # no NumCalc processes are running, so Finish.

        print("... waiting", str(2 * SecondsToInitialize), "s for last processes to finish")
        time.sleep(2 * SecondsToInitialize)

    if CleanUp_afterFinish:  # if Delete all NumCalc executable files from individual source_ folders after completion
        print("Cleaning away not-needed NumCalc executables because CleanUp_afterFinish == True")
        for subdir in os.listdir(root_NumCalc):
            if os.path.isfile(os.path.join(root_NumCalc, subdir, NumCalcRuntime_files[0])):
                # noinspection PyBroadException
                try:  # delete the NumCalc.exe files in the source_ folders
                    for CalcFile in NumCalcRuntime_files:  # delete all NumCalc executable files
                        os.remove(os.path.join(root_NumCalc, subdir, CalcFile))
                # except:
                #     print("    could not delete some file during cleanup")
                finally:
                    CleanUp_afterFinish = True  # do nothing! just to shut up PEP

#  END of all_projects loop.

# finalize
# # print("Total processing time with NumCalcManager:   ", time.strftime("%H:%M:%S", time.time() - StartTime))
print("NumCalc project is COMPLETE.   ", time.strftime("%d %b - %H:%M:%S", time.localtime()))

input(Text_Color_GREEN + 'DONE.  Hit Enter to exit   :)')

# Authors:    The Mesh2HRTF developers  &  Sergejs Dombrovskis
#               Part of the Mesh2HRTF codebase
#
#                        mesh2hrtf.sourceforge.net
#
# Mesh2HRTF is licensed under the GNU Lesser General Public License as
# published by the Free Software Foundation, either version 3 of the License,
# or (at your option) any later version. Mesh2HRTF is distributed in the hope
# that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details. You should have received a
# copy of the GNU LesserGeneral Public License along with Mesh2HRTF. If not,
# see <http://www.gnu.org/licenses/lgpl.html>.
#
# If you use Mesh2HRTF:
# - Provide credits:
#   "Mesh2HRTF, H. Ziegelwanger, ARI, OEAW (mesh2hrtf.sourceforge.net)"
# - In your publication, cite both articles:
#   [1] Ziegelwanger, H., Kreuzer, W., and Majdak, P. (2015). "Mesh2HRTF:
#       Open-source software package for the numerical calculation of head-
#       related transfer functions," in Proceedings of the 22nd ICSV,
#       Florence, IT.
#   [2] Ziegelwanger, H., Majdak, P., and Kreuzer, W. (2015). "Numerical
#       calculation of listener-specific head-related transfer functions and
#       sound localization: Microphone model and mesh discretization," The
#       Journal of the Acoustical Society of America, 138, 208-222.
