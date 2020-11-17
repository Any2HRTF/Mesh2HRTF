% Output2HRTF_RunAll(simulations_dir, overwrite_data,
%                    purge_object_mesh_data, purge_eval_grid_data,
%                    purge_compressed_object_mesh_data,
%                    assume_yes)
%
% Batch run Output2HRTF.m files in all subfolder in `simulations_dir`.
%
% Output2HRTF.m reads simulation results and saves them as SOFA and mat
% files.
%
% Input:
% simulation_dir: str
%   The directory containing Msh2HRTF projects. The default is the current
%   working directory.
% overwrite_data: boolean
%   Run Output2HRTF.m even if data already exists. The default is false
% purge_object_mesh_data: boolen
%   Deletes the uncompressed output from NumCalc after running Output2HRTF
%   (simulations_dir/*/NumCalc/CPU_*_Core_*/be.out/be.*/*_Boundary).
%   Output2HRTF stores information contained in these files in
%   simulations_dir/*/Output2HRTF/ObjectMesh_*.mat. The default is false.
% purge_eval_grid_data: boolen
%   Deletes the uncompressed output from NumCalc after running Output2HRTF
%   (simulations_dir/*/NumCalc/CPU_*_Core_*/be.out/be.*/*_Eval).
%   Output2HRTF stores information contained in these files in
%   simulations_dir/*/Output2HRTF/HRTF*.sofa. The default is false.
% purge_compressed_object_mesh_data: boolen
%   Deletes simulations_dir/*/Output2HRTF/ObjectMesh_*.mat. This file might
%   not always be required but can be large. The default is false.
% assume_yes: boolen
%   Ask the user if data should be purged if assume_yes=false. The default
%   is false.

function Output2HRTF_RunAll(simulations_dir, overwrite_data, ...
                            purge_object_mesh_data, purge_eval_grid_data, ...
                            purge_compressed_object_mesh_data, ...
                            assume_yes)

% globar vars are needed because Output2HRTF.m contains a `clear`
global result_dirs overwrite current_dir nn purge_obj purge_eval purge_compressed

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
if ~exist('purge_compressed_object_mesh_data', 'var')
    purge_compressed_object_mesh_data = false;
end
if ~exist('assume_yes', 'var')
    assume_yes = false;
end

% initalize globals
current_dir = pwd;
overwrite = overwrite_data;
purge_obj = purge_object_mesh_data;
purge_eval = purge_eval_grid_data;
purge_compressed = purge_compressed_object_mesh_data;
result_dirs = dir(simulations_dir);

clear simulations_dir overwrite_data purge_object_mesh_data purge_eval_grid_data

% ask if data should really, really, be deleted
if (purge_obj || purge_eval) && ~assume_yes
    disp('You are about to delete uncompressed simulates results.')
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
    
    % skip if it does not contain Output2HRTF.m
    if ~exist(fullfile(folder, 'Output2HRTF.m'), 'file')
        continue
    end
    
    disp('---------------------------------------------------------------')
    disp(['processing: ' name])
    disp('---------------------------------------------------------------')
    
    cd(folder)
    
    % run Output2HRTF in subfolder
    if (exist(fullfile(folder, 'Output2HRTF'), 'dir') && overwrite) || ...
       ~exist(fullfile(folder, 'Output2HRTF'), 'dir')
   
        Output2HRTF
    else
        disp('Output data already exists')
    end
    
    % recover global variables
    global result_dirs overwrite current_dir nn purge_obj purge_eval purge_compressed %#ok<REDEFGG,TLEV>
    
    % purge uncompressed simulation results
    if ~purge_obj && ~purge_eval && ~purge_compressed
        continue
    end
    
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
        
        % enter results directory
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
                delete(fullfile(freq, '*Eval'))
            end
            
        end
    end
    
    if purge_compressed
        fprintf('\npurging compressed simulation results ... ')
        delete(fullfile(folder, 'Output2HRTF', 'ObjectMesh_*.mat'))
    end
    
    fprintf('done\n\n')
    
end

% return to start directory
cd(current_dir)
% get rid of global variables
clear global