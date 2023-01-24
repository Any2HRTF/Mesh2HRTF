% process_multiple_outputs2hrtf(
%    simulations_dir, overwrite,
%    purge_object_mesh_data, purge_eval_grid_data, purge_frequencies
%    purge_compressed_object_mesh_data, purge_hrft
%    assume_yes)
%
% Batch run output2hrtf.m in all subfolders in `simulations_dir`
% and delete temporary data if desired.
%
% output2hrtf.m reads simulation results and saves them as SOFA and mat
% files.
%
% Input:
% simulation_dir: str
%   The directory containing Msh2HRTF projects. The default is the current
%   working directory.
% overwrite: boolean
%   Run output2hrtf.m even if data already exists. The default is false
% purge_object_mesh_data: boolen
%   Delete the uncompressed output from NumCalc after running output2hrtf
%   (simulations_dir/*/NumCalc/CPU_*_Core_*/be.out/be.*/*_Boundary).
%   output2hrtf stores information contained in these files in
%   simulations_dir/*/output2hrtf/ObjectMesh_*.mat. The default is false.
% purge_eval_grid_data: boolen
%   Delete the uncompressed output from NumCalc after running output2hrtf
%   (simulations_dir/*/NumCalc/CPU_*_Core_*/be.out/be.*/*_Eval).
%   output2hrtf stores information contained in these files in
%   simulations_dir/*/output2hrtf/HRTF*.sofa. The default is false.
% purge_frequencies: boolean
%   Delete the frequencies written to the text files inside the fe.out
%   folder. This information is redundant and not required for anything.
%   The default is false.
% purge_compressed_object_mesh_data: boolen
%   Delete simulations_dir/*/output2hrtf/ObjectMesh_*.mat. This file might
%   not always be required but can be large. The default is false.
% purge_hrtf : boolean
%   Delete the HRTF SOFA file. This file might not be needed if the HRIRs
%   were calculated as well. The default is false.
% assume_yes: boolen
%   Ask the user if data should be purged if assume_yes=false. The default
%   is false.

function process_multiple_outputs2hrtf(...
    simulations_dir, overwrite_data, ...
    purge_object_mesh_data, purge_eval_grid_data, purge_frequencies, ...
    purge_compressed_object_mesh_data, purge_hrtf, ...
    assume_yes)

% globar vars are needed because output2hrtf.m contains a `clear`
global result_dirs overwrite current_dir nn purge_obj purge_eval purge_freq purge_compressed purge_HRTF folder

% default parameters
if ~exist('simulations_dir', 'var')
    simulations_dir = pwd;
end
if ~exist('overwrite_data', 'var')
    overwrite_data = false;
end
if ~exist('purge_object_mesh_data', 'var')
    purge_object_mesh_data = false;
end
if ~exist('purge_eval_grid_data', 'var')
    purge_eval_grid_data = false;
end
if ~exist('purge_frequencies', 'var')
    purge_frequencies = false;
end
if ~exist('purge_compressed_object_mesh_data', 'var')
    purge_compressed_object_mesh_data = false;
end
if ~exist('purge_hrtf', 'var')
    purge_hrtf = false;
end
if ~exist('assume_yes', 'var')
    assume_yes = false;
end

% initalize globals
current_dir = pwd;
overwrite = overwrite_data;
purge_obj = purge_object_mesh_data;
purge_eval = purge_eval_grid_data;
purge_freq = purge_frequencies;
purge_compressed = purge_compressed_object_mesh_data;
purge_HRTF = purge_hrtf;
result_dirs = dir(simulations_dir);

clear simulations_dir overwrite_data purge_object_mesh_data ...
    purge_eval_grid_data purge_hrtf

% ask if data should really, really, be deleted
if (purge_obj || purge_eval || purge_compressed || purge_HRTF) && ...
        ~assume_yes
    disp('Delete selected files after calculating the results:')
    if purge_obj; disp("- Raw pressure on the mesh"); end
    if purge_eval; disp("- Raw pressure on evaluation grid"); end
    if purge_freq; disp("- Frequencies"); end
    if purge_compressed; disp("- Compressed pressure on the mesh"); end
    if purge_HRTF; disp("- HRTF SOFA file"); end
    yes = input('Continue (y/n)?: ', 's');
    if ~strcmpi(yes, 'y')
        error('Aborted purging data')
    end
end

% loop directories
for nn = 1:numel(result_dirs)
    % current directory
    base = result_dirs(nn).folder;
    name = result_dirs(nn).name;
    folder = fullfile(base, name);

    % skip if it is not a folder
    if ~isfolder(folder) || strcmp(name(1), '.')
        continue
    end

    % skip if it does not contain output2hrtf.m
    if ~exist(fullfile(folder, 'parameters.json'), 'file')
        continue
    end

    disp('---------------------------------------------------------------')
    disp(['processing: ' name])
    disp('---------------------------------------------------------------')

    cd(folder)

    % run output2hrtf in subfolder
    if (exist(fullfile(folder, 'output2hrtf'), 'dir') && overwrite) || ...
       ~exist(fullfile(folder, 'output2hrtf'), 'dir')

        output2hrtf
    else
        disp('Output data already exists')
    end

    % recover global variables
    global result_dirs overwrite current_dir nn purge_obj purge_eval purge_compressed purge_HRTF folder%#ok<REDEFGG,TLEV>

    % purge uncompressed simulation results
    if purge_obj || purge_eval
        fprintf('\npurging uncompressed simulation results ... ')

        % NumCalc CPU* folder
        cores = dir(fullfile(folder, 'NumCalc', 'CPU*'));
        for cc = 1:numel(cores)
            % current core folder
            core = fullfile(cores(cc).folder, cores(cc).name);
            % skip any non folders
            if ~isfolder(core)
                continue
            end

            % remove data in be.out
            if isfolder(fullfile(core, 'be.out'))
                % remove entire folder ...
                if purge_obj && purge_eval
                    rmdir(fullfile(core, 'be.out'), 's')
                    continue
                end

                % or ... enter results directory and remove separate files
                cd(fullfile(core, 'be.out'))

                % loop over frequencies
                freqs = dir(fullfile(core, 'be.out', 'be.*'));
                for ff = 1:numel(freqs)
                    freq = fullfile(freqs(ff).folder, freqs(ff).name);
                    % skip any non folders
                    if ~isfolder(freq)
                        continue
                    end

                    if purge_obj
                        delete(fullfile(freq, '*Boundary'))
                    end
                    if purge_eval
                        delete(fullfile(freq, '*EvalGrid'))
                    end

                end
            end

            % purge frequencies
            if isfolder(fullfile(core, 'fe.out')) && purge_freq
                rmdir(fullfile(core, 'fe.out'), 's')
            end
        end
    end

    % purge compressed simulation results
    if purge_compressed
        fprintf('\npurging compressed simulation results ... ')
        delete(fullfile(folder, 'output2hrtf', 'ObjectMesh_*.mat'))
    end

    % purge compressed simulation results
    if purge_HRTF
        fprintf('\npurging HRTF SOFA file ... ')
        delete(fullfile(folder, 'output2hrtf', 'HRTF_*.sofa'))
    end

    fprintf('done\n\n')

end

% return to start directory
cd(current_dir)
% get rid of global variables
clear global