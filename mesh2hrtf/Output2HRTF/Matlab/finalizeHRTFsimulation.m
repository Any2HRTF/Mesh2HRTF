function finalizeHRTFsimulation(varargin)
%  FINALIZEHRTFSIMULATION(varargin)
%  merge SOFA data from two separate Mesh2HRTF simulations for left & right ear.
%  For usage details, please check the documentation!
%
% A - no inputs = scan start folder for 2 projects to merge.
% B - 1 input  = scan given folder for 2 projects to merge.
% C - 2 inputs = Merge 2 projects that were given as input.
%
% mode A - no inputs = (recommended usage mode) scans start folder for 2 projects to merge,
%                       + usually executes "Output2HRTF.m" to run full pre-processing in one go.
%   run this .m file from a dedicated folder (for example use folder "Finalize_HRTF_simulation")
%   Move into the "Finalize_HRTF_simulation" folder exactly the 2 project folders that need to be merged.
%   run this "finalizeHRTFsimulation.m" with Matlab.
%   Merged SOFA files will be saved in a folder next to the project that was found.
%
% mode B - Input1 only = scan Input1 folder for 2 projects to merge.
%   Move into any folder exactly the 2 project folders that need to be merged.
%   run this "finalizeHRTFsimulation.m" file & specify "input_1" folder that contains projects to merge.
%       This script searches and merges the SOFA files from the 2 projects it finds inside "input_1" path.
%   Merged SOFA file will be saved in pwd.
%
% mode C - 2 inputs = Merge the 2 projects that were given as input_1 and input_2.
%   Merged SOFA file will be written in pwd.
%
% migrated from Python API from Sergejs Dombrovskis
% author(s): Katharina Pollack, June 2022

% check input and input mode
switch numel(varargin)
    case 0 % mode 'no inputs' = scan working directory for 2 projects to merge.
        start_folder = dir;
        for ii = 3:numel(start_folder)
            % check for folder name
            [fpath, fname, fext] = fileparts(start_folder(ii).name);
            if isempty(fext)
                % save folder name to variable
                folder_names{end+1} = fullfile(fpath, fname);
            end
        end
        if numel(folder_names) ~= 2
            error('More than two folders were found, only two expected. Program aborted.');
        end
    case 1 % mode '1 input' = scan given folder for 2 projects to merge.
        if ~ischar(varargin{1})
            error('Allowed input parameter is string.');
        elseif ~isfolder(varargin{1})
            error('Input parameter has to be valid folder name.');
        end
        start_folder = varargin{1};
        
        % pretty sure this can be solved in a more elegant way (than
        % copypasting the code from above ...)
        for ii = 3:numel(start_folder)
            % check for folder name
            [fpath, fname, fext] = fileparts(start_folder(ii).name);
            if isempty(fext)
                % save folder name to variable
                folder_names{end+1} = fullfile(fpath, fname);
            end
        end
        if numel(folder_names) ~= 2
            error('More than two folders were found, only two expected. Program aborted.');
        end
    case 2 % mode '2 inputs' = Merge 2 projects that were given as input.
        if ~ischar(varargin{1}) || ~ischar(varargin{2})
            error('Allowed input parameter is string.');
        elseif ~isfolder(varargin{1}) || ~isfolder(varargin{2})
            error('Both input parameters have to be valid folder names.');
        end
        folder_names{1} = varargin{1};
        folder_names{2} = varargin{2};
    otherwise
        error(['Wrong number of inputs given. Valid number of input arguments are:\n', ...
            '%s\n%s\n%s'], ...
            'A - no inputs = scan start folder for 2 projects to merge.', ...
            'B - 1 input  = scan given folder for 2 projects to merge.', ...
            'C - 2 inputs = Merge 2 projects that were given as input.');
end
clear ii

% def scan_for_errors(project_path):  # find any broken data in this project
%     NOTE: made only to work with 1 source!!!
% scan project for errors
% run Output2HRTF_Main in both folders
for ii = 1:2
    cd(fullfile(folder_names{ii}, 'NumCalc'))
    % check for number of sources
    NCdir_tmp = dir;
    if numel(NCdir_tmp) > 3 % including '.' and '..' entries
        error(['%s consists of two source folders.\n', ...
            'This script only merges SOFA files of two separately calculated sources.\n', ...
            'See documentation for details.'], folder_names{ii});
    end
    cd(NCdir_tmp{3}) % ignoring '.' and '..' entries
    % check for NC*.out files
    sourcedir_tmp = dir('NC*.out');
    for jj = 1:numel(sourcedir_tmp)
        NCout_files{end+1} = fullfile(sourcedir_tmp(jj).folder, sourcedir_tmp(jj).name);
    end
    % set keys for regexps
    frequency_key = ', Frequency = ';
    non_convergence_key = 'Warning: Maximum number of iterations is reached!';
    gaussean_lim_key = 'Too many integral points in the theta-direction!';
    happy_end_key = 'End time: ';
    for kk = 1:numel(NCout_files)
        % open each NC*.out file and scan for errors (i.e., search for keys in files)
        fid=fopen(NCout_files{kk});
        counter = 0; happy_end = 0; corrupt_freq = 0; lowest_corrupt_frq = 0; freq = [];
        while ~feof(fid)
            counter = counter+1;
            line=fgets(fid);
            switch line
                case frequency_key
                    line_nr = counter; % save line number as well
                    freq = [freq; line(length(frequency_key)+1:end-3)]; % omit last 3 chars (' Hz')
                case non_convergence_key
                    if freq(end) > 24e3 % affects very high frequencies only
                        warning(['Non-Convergence detected at over 24kHz  (at ', num2str(freq(end)), ' Hz']);
                    else
                        warning(['PROBLEM!! - Non-Convergence issue in important range at ', num2str(freq(end)), ' Hz']);
                        corrupt_freq = [corrupt_freq; freq(end)]; % note which at frequency the calculation failed
                    end
                case gaussean_lim_key
                    if freq(end) > 24e3 % affects very high frequencies only
                        warning(['Gaussean points limit issue detected at over 24kHz  (at ', num2str(freq(end)), ' Hz']);
                    else
                        warning(['PROBLEM!! - Gaussean points limit issue in important range at ', num2str(freq(end)), ' Hz']);
                        corrupt_freq = [corrupt_freq; freq(end)]; % note which at frequency the calculation failed
                    end
                case happy_end_key
                    happy_end = 1;
                otherwise
            end
        end
        fclose(fid);
        if happy_end == 1
            disp(['Great! The calculation of', NCout_files{kk}, 'finished successfully.']);
        else
            disp('FYI: The calculation of', NCout_files{kk}, 'did not finish successfully.');
        end

        if numel(corrupt_freq) > 1 % we have corrupt data
            if lowest_corrupt_frq == 0
                lowest_corrupt_frq = min(corrupt_freq);
            elseif lowest_corrupt_frq > min(corrupt_freq)
                lowest_corrupt_frq = min(corrupt_freq); % keep only the lowest problematic frequency
            end
        end
    end
    
    % run Output2HRTF script in project folder
    cd(folder_names{ii})
    Output2HRTF;
end


% determine left and right ear data
% source position along interaural axis. Positive = left, negative = right

% merge to one SOFA file
% merge_write_sofa(sofa_left, sofa_right, basepath, filename, sofa_type='HRTF'):
% merge_write_sofa(sofa_left, sofa_right, basepath, filename, sofa_type='HRIR'):

% write merged SOFA file to pwd

end
% EOF