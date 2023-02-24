% process_multiple_outputs2hrtf( simulations_dir, overwrite,
%    purge_object_mesh_data, purge_eval_grid_data, purge_frequencies,
%    purge_compressed_object_mesh_data, purge_hrft, assume_yes)
%
% Batch run of `output2hrtf.m` in all subfolders of `simulations_dir`
% and optionally delete temporary data.
%
% Input parameters:
%
% simulation_dir: Directory containing the Meh2HRTF project. Default: current
%   working directory.
%
% overwrite: If true, output2hrtf will be executed even if data already exists. 
%   Default: false
%
% purge_object_mesh_data: If true, the uncompressed object mesh data from 
%   NumCalc will be deleted after running output2hrtf. The deleted files are 
%   usually files named simulations_dir/*/NumCalc/CPU_*_Core_*/be.out/be.*/*_Boundary.
%   Note that output2hrtf stores this information in
%   simulations_dir/*/output2hrtf/ObjectMesh_*.mat. Default: false
%
% purge_eval_grid_data: If true, uncompressed evaluation-grid data from 
%   NumCalc will be deleted after running output2hrtf. The deleted files
%   are usually named simulations_dir/*/NumCalc/CPU_*_Core_*/be.out/be.*/*_Eval.
%   Note that output2hrtf stores this information in
%   simulations_dir/*/output2hrtf/HRTF*.sofa. Default: false
%
% purge_frequencies: If true, frequencies written to the text files inside the 
%   fe.out folder will be deleted. This information is redundant and not 
%   required for anything. Default: false
%
% purge_compressed_object_mesh_data: if true, the compressed object mesh data
%   will be deleted after running output2hrtf. The deleted files are
%   simulations_dir/*/output2hrtf/ObjectMesh_*.mat. These files are not always
%   required but can be large. Default: false
%
% purge_hrtf: If true, the HRTF SOFA file will be deleted. This file might be 
%   not required if the HRIR SOFA files was calculated as well. Default: false.
%
% assume_yes: If false, the user will be asked to purge the data. Default: false

% This file is part of the Mesh2HRTF software package developed by the
% Mesh2HRTF Developer Team (https://mesh2hrtf.org) and licensed under the 
% EUPL, Version 1.2, or, as soon as approved by the European Commission, 
% subsequent versions of the EUPL. Details on the license can be found 
% in the file "license.txt" provided with Mesh2HRTF package
% or at https://joinup.ec.europa.eu/collection/eupl/eupl-text-eupl-12
%
% You may not use this work except in compliance with the license.
% Unless required by applicable law or agreed to in writing, software 
% distributed under the license is distributed on an "AS IS" basis,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.

% #Author: Harald Ziegelwanger (ARI, ÖAW): 2015, original implementation
% #Author: Fabian Brinkmann (TU-Berlin): 2020, integration in Mesh2HRTF 1.x
% #Author: Fabian Brinkmann (TU-Berlin): 2022, various improvements
% #Author: Piotr Majdak (ARI, ÖAW): 2023, help text, license boiler plate



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