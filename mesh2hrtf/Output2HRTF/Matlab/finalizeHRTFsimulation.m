function finalizeHRTFsimulation(varargin)
%  "finalize_hrtf_simulation.py" merges SOFA data from two separate Mesh2HRTF simulations for LEFT & RIGHT ear.
%  it uses a lot of example code from Mesh2HRTF "Output2HRTF_Main.py" & "NumCalcManager.py"
%
%                               HOW TO USE:
% Main tutorial is here:   https://sourceforge.net/p/mesh2hrtf/wiki/Basic_HRTF_post_processing
%
% A - no inputs = scan start folder for 2 projects to merge.
% B - 1 input  = scan given folder for 2 projects to merge.
% C - 2 inputs = Merge 2 projects that were given as input.
%
% mode A - no inputs = (recommended usage mode) scans start folder for 2 projects to merge,
%                       + usually executes "Output2HRTF.py" to run full pre-processing in one go.
%   1- run this .py file from a dedicated folder (for example use folder "Finalize_HRTF_simulation")
%   2- Move into the "Finalize_HRTF_simulation" folder exactly the 2 project folders that need to be merged.
%   3- "open/run" this "finalize_hrtf_simulation.py" with Python.
%   4- Merged SOFA files will be saved in a folder next to the project that was found.   DONE.
%
% mode B - Input1 only = scan Input1 folder for 2 projects to merge.
%   1- Move into any folder exactly the 2 project folders that need to be merged.
%   2- "run" this "finalize_HRTF_simulation.py" file & specify "input_1" folder that contains projects to merge.
%       This script searches and merges the SOFA files from the 2 projects it finds inside "input_1" path.
%   3- Merged SOFA files will be saved in a folder next to the project that was found.   DONE.
%
% mode C - 2 inputs = Merge the 2 projects that were given as input_1 and input_2.
%
%   made for Python 3.7+                    see license details at the end of the script.
%       v0.50    2022-01-23     First working version. Tested on Windows 10. (by Sergejs Dombrovskis)
%       v0.90    2022-01-24     Mode-A1 is ready + most descriptions in place. (by S.D.)
%       v0.93    2022-02-06     added Mode-A2 and Mode-C. (by S.D.)
%       v1.00    2022-03-15     added Mode-A0 with automatic execution of "Output2HRTF.py" + small improvements +
%                                added multi-SamplingRate HRIR outputs + extra sanity checks (by S.D.)
%       v1.10    2022-03-17     added diagnostic plotting of the merged HRIR.sofa file (by S.D.)
%       v1.20    2022-03-20     added effective workaround to produce usable data even if simulation encountered the
%                                "Non-Convergence issue" - only in A0 mode (by S.D.)
%       v1.30    2022-03-26     added more and easier-to-use troubleshooting messages (more beginner-friendly),
%                                added universal handling of all possible simulation errors (by S.D.)
%       v1.31    2022-03-22     improved compatibility with Mac-OS + minor bugfix (by S.D.)
%       v1.32    2022-03-31     disabled "multi-SamplingRate" HRIR saving because the current code is wrong (by S.D.)
%       v1.50    2022-04-02     rewrite to make "multi-SamplingRate" HRIR saving work correctly (by S.D.)
%       v1.60    2022-04-20     Compatibility fixes to support latest version of "output_to_hrtf.py". Plus
%                                   renamed the script from the initial "join_sofa_files.py" name because of much
%                                   broader scope than in the initial v0.93 (by S.D.)
%       v1.70    2022-05-02     Improved repeated executions - now it is safe to re-run the script, added setting
%                                   "Re_Process_Output2HRTF", plus minor improvements for exception handling (by S.D.)
%       v2.00    2022-05-12     rewrite to make "multi-SamplingRate" HRIR use the new "output_to_hrtf" functionality
%                                   and avoid the crash due to other changes in recent "output_to_hrtf".
%                                   Added command line arguments, and simplified too much confusing modes.
%                                   Bugfix in HRTF merging based on code from "utils._merge_sofa_files" (by S.D.)
%       v2.05    2022-05-16     Various useful improvements & fixes (by S.D.)

% check input and input mode
% # different input scenarios:
% # A - no inputs = scan start folder for 2 projects to merge.
% # B - 1 input  = scan given folder for 2 projects to merge.
% # C - 2 inputs = Merge 2 projects that were given as input.


% def merge_write_sofa(sofa_left, sofa_right, basepath, filename, sofa_type='HRTF'):
% """Write complex pressure or impulse responses as SOFA file."""
%    adapted from "write_to_sofa()" in "Output2HRTF_Main.py"

% cart2sph conversion with DEGREE output
% def cart2sph_in_deg(xyz):

% def scan_for_errors(project_path):  # find any broken data in this project
%     # NOTE: made only to work with 1 source!!!

% def plot_hrir_pair(sofa_data_path, noise_floor=-70):
%     """
%     Plots both left and right ear HRIRs in an easy to diagnose format
% 
%     # based on tutorial: https://sofar.readthedocs.io/en/latest/working_with_sofa_files.html#working-with-sofa
% 
%     Params:
%     sofa_data_path : full path to the sofa HRIR file to load and plot
% 
%     noise_floor : were to limit colormap to provide easy to read image
%         = -50  # noise floor value used in "SOFAplotHRTF.m" from SOFA Matlab API
%         = -70  # value that visually matches the look of "SOFAplotHRTF.m" from SOFA MatlabAPI
% 
%     """



end
% EOF