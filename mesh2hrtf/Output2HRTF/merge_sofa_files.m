function merge_sofa_files(varargin)
%  merge_sofa_files(varargin)
%  merge SOFA data from two separate Mesh2HRTF simulations for left & right ear.
%  For usage details, please check the documentation!
%
%  A - no inputs ... scan pwd for 2 projects to merge
%  B - 1 input ..... scan given folder for 2 projects to merge
%  C - 2 inputs .... merge 2 projects that were given as input 1 and 2
%  D - 3 inputs .... merge 2 projects specified in input 1 and 2 and save to folder in input 3
%
%  MODE A - no inputs = scan start folder for 2 projects to merge,
%                       + usually executes "output2hrtf.m" to run full pre-processing in one go.
%   run this .m file from a dedicated folder (for example use folder "Finalize_HRTF_simulation")
%   Move into the "Finalize_HRTF_simulation" folder exactly the 2 project folders that need to be merged.
%   run this "merge_sofa_files.m" with Matlab.
%   Merged SOFA files will be saved in a folder next to the project that was found.
%
%  MODE B - Input1 only = scan Input1 folder for 2 projects to merge.
%   Move into any folder exactly the 2 project folders that need to be merged.
%   run this "merge_sofa_files.m" file & specify "input_1" folder that contains projects to merge.
%       This script searches and merges the SOFA files from the 2 projects it finds inside "input_1" path.
%   Merged SOFA file will be saved in pwd.
%
%  MODE C - 2 inputs = merge the 2 projects that were given as input_1 and input_2.
%   Merged SOFA file will be written in pwd.
%
%  MODE D - 3 inputs = (recommended usage mode) merge 2 projects and save to folder in input 3
%   Merged SOFA file will be written to input 3
%
%  migrated from finalize_hrtf_simulation.py from Sergejs Dombrovskis
%  author(s): Katharina Pollack, June 2022

pwdirectory = pwd;
targetdir = {};

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
    case 3 % mode '3 inputs' = Merge 2 projects and safe merged SOFA file to input 3
        if ~ischar(varargin{1}) || ~ischar(varargin{2}) || ~ischar(varargin{3})
            error('Allowed input parameter is string.');
        elseif ~isfolder(varargin{1}) || ~isfolder(varargin{2}) || ~isfolder(varargin{3})
            error('All input parameters have to be valid folder names.');
        end
        folder_names{1} = varargin{1};
        folder_names{2} = varargin{2};
        targetdir = varargin{3};
    otherwise
        error(['Wrong number of inputs given. Valid number of input arguments are:\n', ...
            '%s\n%s\n%s\n%s'], ...
            'A - no inputs = scan start folder for 2 projects to merge', ...
            'B - 1 input  = scan given folder for 2 projects to merge', ...
            'C - 2 inputs = merge 2 projects that were given as input', ...
            'D - 3 inputs = merge 2 projects specified in input 1 and 2 and save to folder in input 3');
end
clear ii

% create folder name for merged SOFA
if isempty(targetdir)
    if strcmp(folder_names{1}(end-1:end), '_L') || strcmp(folder_names{1}(end-1:end), '_R')
        targetdir = strcat(folder_names{1}(1:end-2), '_merged');
    elseif strcmpi(folder_names{1}, pwdirectory) || strcmpi(folder_names{2}, pwdirectory)
        error(['Please specify a target directory for the merged SOFA file to be stored as a third input. ', ...
            'It should not be one of the input directories.']);
    else
        targetdir = pwdirectory;
    end
end

% find any broken data in this project
for ii = 1:2
    NCout_files = {};
    disp(['Checking for errors in ', folder_names{ii}, ' ...']);
    cd(fullfile(folder_names{ii}, 'NumCalc'))
    % check for number of sources
    NCdir_tmp = dir;
    if numel(NCdir_tmp) > 3 % including '.' and '..' entries
        error(['%s consists of two source folders.\n', ...
            'This script only merges SOFA files of two separately calculated sources.\n', ...
            'See documentation for details.'], folder_names{ii});
    end
    cd(NCdir_tmp(3).name) % ignoring '.' and '..' entries
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
    lowest_corrupt_frq = [0, 0];
    for kk = 1:numel(NCout_files)
        % open each NC*.out file and scan for errors (i.e., search for keys in files)
        fid=fopen(NCout_files{kk});
        counter = 0; happy_end = 0; corrupt_freq = 0; freq = [];
        while ~feof(fid)
            counter = counter+1;
            line=fgetl(fid); % read line without newline character
            if contains(line, frequency_key)
                line_nr = counter; % save line number as well
                startIdx=regexp(line, frequency_key);
                freq = [freq; str2double(line(startIdx+length(frequency_key):end-3))]; % omit last 3 chars (' Hz')
            end
            if contains(line, non_convergence_key)
                if freq(end) > 24e3 % affects very high frequencies only
                    warning(['Non-Convergence detected at over 24kHz  (at ', num2str(freq(end)), ' Hz']);
                else
                    warning(['PROBLEM!! - Non-Convergence issue in important range at ', num2str(freq(end)), ' Hz']);
                    corrupt_freq = [corrupt_freq; freq(end)]; % note which at frequency the calculation failed
                end
            end
            if contains(line, gaussean_lim_key)
                if freq(end) > 24e3 % affects very high frequencies only
                    warning(['Gaussean points limit issue detected at over 24kHz  (at ', num2str(freq(end)), ' Hz']);
                else
                    warning(['PROBLEM!! - Gaussean points limit issue in important range at ', num2str(freq(end)), ' Hz']);
                    corrupt_freq = [corrupt_freq; freq(end)]; % note which at frequency the calculation failed
                end
            end
            if contains(line, happy_end_key)
                happy_end = 1;
            end
        end
        if happy_end == 1
            disp(['Great! The calculation of', NCout_files{kk}, ' finished successfully.']);
        else
            warning(['FYI: The calculation of', NCout_files{kk}, ' did not finish successfully.']);
        end
    fclose(fid);
    end

    if numel(corrupt_freq) > 1 % we have corrupt data
        if lowest_corrupt_frq(ii) == 0
            lowest_corrupt_frq(ii) = min(corrupt_freq);
        elseif lowest_corrupt_frq(ii) > min(corrupt_freq)
            lowest_corrupt_frq(ii) = min(corrupt_freq); % keep only the lowest problematic frequency
        end
    end
end
clear ii

if (lowest_corrupt_frq(1) > 1 && lowest_corrupt_frq(1) < 250) || ...
        (lowest_corrupt_frq(2) > 1 && lowest_corrupt_frq(2) < 250)
    error(['Practically all simulation instances failed. There is no usable data here. ', ...
        'See documentation for troubleshooting or open an issue in the repository.']);
end

% overwrite lowest_corrupt_frq with unified information
if min(lowest_corrupt_frq) == 0
    lowest_corrupt_frq = max(lowest_corrupt_frq);
else
    lowest_corrupt_frq = min(lowest_corrupt_frq);
end

for ii = 1:2
    % run output2hrtf script in project folder
    cd(folder_names{ii})
    if isfolder(fullfile(folder_names{ii}, 'Output2HRTF')) && ~isempty(dir(fullfile(folder_names{ii}, 'Output2HRTF', '*.sofa')))
        warning(['There are already SOFA files in ', fullfile(folder_names{ii}, 'Output2HRTF'), ...
            '. If you want to overwrite them, please delete the .sofa files and rerun merge_sofa_files.m']);
    else
        disp(['Running output2hrtf in folder ', folder_names{ii}, ' ...']);
        output2hrtf;
    end
    % determine left and right ear data
    % source position along interaural axis. Positive = left, negative = right
    disp(['Determine L/R ear of project ', folder_names{ii}, ' ...']);

    params = fullfile(folder_names{ii}, 'parameters.json');
    fid = fopen(params);
    params = fread(fid,inf);
    params = char(params');
    fclose(fid);
    params = jsondecode(params);

    switch params.sourceType
        case 'Left ear'
            sofa_file_L = folder_names{ii};
            sofa_files_tmp = dir(fullfile(folder_names{ii}, 'Output2HRTF', '*.sofa'));
            for jj = 1:numel(sofa_files_tmp)
                sofa_files_L{jj} = fullfile(sofa_files_tmp(jj).folder, sofa_files_tmp(jj).name);
            end
        case 'Right ear'
            sofa_file_R = folder_names{ii};
            sofa_files_tmp = dir(fullfile(folder_names{ii}, 'Output2HRTF', '*.sofa'));
            for jj = 1:numel(sofa_files_tmp)
                sofa_files_R{jj} = fullfile(sofa_files_tmp(jj).folder, sofa_files_tmp(jj).name);
            end
        otherwise
            error('This function only works for sourceTypes Left ear and Right ear')
    end
end
clear ii

% merge to one SOFA file
% merge_write_sofa(sofa_left, sofa_right, basepath, filename, sofa_type='HRTF'):
% merge_write_sofa(sofa_left, sofa_right, basepath, filename, sofa_type='HRIR'):
disp('Merge SOFA file(s) ...');
for ii = numel(sofa_files_L):-1:1  % loop for every SOFA file in one of the projects
    sofa_read_L = SOFAload(sofa_files_L{ii});
    sofa_read_R = SOFAload(sofa_files_R{ii});
    % no need to check whether there are the same SOFA files in L/R folders
    % since we created them from one file

    % detect data type
    if strcmp(sofa_read_L.GLOBAL_DataType, 'FIR')
        SOFA_type = 'HRIR';
        % sanity check whether files can be merged
        if sofa_read_L.Data.SamplingRate ~= sofa_read_R.Data.SamplingRate
            error('Sampling rates of the two input SOFA files do not match! Merging aborted.');
        elseif size(sofa_read_L.Data.IR, 2) ~= 1 || size(sofa_read_R.Data.IR, 2) ~= 1
            error('Input SOFA files contain data for more than one receiver! Merging aborted.');
        elseif size(sofa_read_L.Data.IR, 3) ~= size(sofa_read_R.Data.IR, 3)
            error('Sampling frequencies of the two input SOFA files do not match! Merging aborted.')
        end
    else
        SOFA_type = 'HRTF';
        % sanity check whether files can be merged
        if unique((sofa_read_L.N - sofa_read_R.N) > 1e-10)
            error('Frequency vectors of the two input SOFA files do not match! Merging aborted.');
        end
    end

    % create empty SOFA object
    if strcmp(SOFA_type, 'HRTF')
        Obj = SOFAgetConventions('SimpleFreeFieldHRTF');
        % write data to SOFA object
        Obj.Data.Real = sofa_read_L.Data.Real;
        Obj.Data.Real(:,2,:) = sofa_read_R.Data.Real;
        Obj.Data.Imag = sofa_read_L.Data.Imag;
        Obj.Data.Imag(:,2,:) = sofa_read_R.Data.Imag;
        Obj.N = sofa_read_L.N;
        fs_out = 2*Obj.N(end); % specify sampling frequency according to Nyquist
    elseif strcmp(SOFA_type, 'HRIR')
        Obj = SOFAgetConventions('SimpleFreeFieldHRIR');
        % write data to SOFA object
        Obj.Data.IR = sofa_read_L.Data.IR;
        Obj.Data.IR(:,2,:) = sofa_read_R.Data.IR;
        Obj.Data.SamplingRate = sofa_read_L.Data.SamplingRate;
        Obj.Data.Delay = [sofa_read_L.Data.Delay, sofa_read_R.Data.Delay];
        fs_out = Obj.Data.SamplingRate;
    else
        Obj = SOFAgetConventions('General'); % general SOFA object, no restrictions at all
    end

    % write meta data
    Obj.GLOBAL_ApplicationName = sofa_read_L.GLOBAL_ApplicationName;
    Obj.GLOBAL_ApplicationVersion = sofa_read_L.GLOBAL_ApplicationVersion;
    Obj.GLOBAL_History = sofa_read_L.GLOBAL_History;

    if isempty(sofa_read_L.GLOBAL_Title)
        Obj.GLOBAL_Title = ['untitled ', SOFA_type, ' SOFA data'];
    else
        Obj.GLOBAL_Title = sofa_read_L.GLOBAL_Title;
    end

    if isfield(sofa_read_L, 'GLOBAL_Origin') && isempty(sofa_read_L.GLOBAL_Origin) % added for spat5
        Obj.GLOBAL_Origin = 'Mesh2HRTF_simulation';
    elseif isfield(sofa_read_L, 'GLOBAL_Origin')
        Obj.GLOBAL_Origin = sofa_read_L.GLOBAL_Origin;
    end

    % source and receiver data
    if strcmp(sofa_read_L.SourcePosition_Type, 'cartesian') && strcmp(sofa_read_L.SourcePosition_Units, 'metre')
        [pos_tmp(:,1), pos_tmp(:,2), pos_tmp(:,3)] = cart2sph(sofa_read_L.SourcePosition(:,1), sofa_read_L.SourcePosition(:,2), sofa_read_L.SourcePosition(:,3));
        Obj.SourcePosition = [rad2deg(pos_tmp(:,1:2)), pos_tmp(:,3)];
        Obj.SourcePosition_Units = 'degree, degree, metre';
        Obj.SourcePosition_Type = 'spherical';
        clear pos_tmp
    else
        Obj.SourcePosition = sofa_read_L.SourcePosition;
        Obj.SourcePosition_Units = sofa_read_L.SourcePosition_Units;
        Obj.SourcePosition_Type = sofa_read_L.SourcePosition_Type;
    end

    Obj.ReceiverPosition = [sofa_read_L.ReceiverPosition; sofa_read_R.ReceiverPosition];
    Obj.ReceiverPosition_Units = sofa_read_L.ReceiverPosition_Units;
    Obj.ReceiverPosition_Type = sofa_read_L.ReceiverPosition_Type;

    % write merged SOFA file to targetdir
    disp(['Write merged SOFA file to ', targetdir, SOFA_type, '_', num2str(fs_out), 'fs.sofa ...']);
    % create targetdir if it does not exist yet
    if ~isfolder(targetdir)
        mkdir(targetdir)
    end
    cd(targetdir);

    % save as SOFA file
    Obj=SOFAupdateDimensions(Obj);
    SOFAsave(fullfile(targetdir, [SOFA_type, '_', num2str(fs_out), 'fs.sofa']), Obj, 9); % choose highest compression of 9
end

disp(['Done! Merged SOFA files saved in ', targetdir]);

end
% EOF