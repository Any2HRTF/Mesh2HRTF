function output2hrtf(folder)
%   [] = OUTPUT2HRTF_MAIN(folder)
%
%   Process NumCalc output and write data to disk. All parameters are read
%   from parameters.json.
%
%   Input:
%       folder ... The Mesh2HRTF project folder [string, the default is
%                  the current working directory]

%% ----------------------------load meta data------------------------------

if ~exist('folder', 'var')
    folder = pwd;
end

% load parameters
if ~isfile(fullfile(folder, 'parameters.json'))
    error('parameters.json not found.')
end
params = fullfile(folder, 'parameters.json');
fid = fopen(params);
params = fread(fid,inf);
params = char(params');
fclose(fid);
params = jsondecode(params);
if all(size(params.sourceCenter) == [3 1])
    params.sourceCenter = params.sourceCenter';
end

% output directory
if ~exist(fullfile(folder, 'Output2HRTF'), 'dir')
    mkdir(fullfile(folder, 'Output2HRTF'))
end

fprintf('Load meta data ...\n');

% get the evaluation grids
fprintf('Get evaluation grids ...\n');
evaluationGrids = dir(fullfile(folder, 'EvaluationGrids'));
evaluationGrids = evaluationGrids(~cellfun(@(x) strncmp(x, '.', 1), {evaluationGrids.name}));

evaluationGridsNumNodes = 0;
for ii=1:length(evaluationGrids)
    tmpNodes=importdata(fullfile(folder, 'EvaluationGrids', evaluationGrids(ii).name, 'Nodes.txt'),' ',1);
    tmpElements=importdata(fullfile(folder, 'EvaluationGrids', evaluationGrids(ii).name, 'Elements.txt'),' ',1);
    evaluationGrids(ii).nodes = tmpNodes.data;
    evaluationGrids(ii).elements = tmpElements.data;
    evaluationGrids(ii).num_nodes = size(tmpNodes.data, 1);
    evaluationGridsNumNodes = evaluationGridsNumNodes + evaluationGrids(ii).num_nodes;
end

% get the object mesh
fprintf('Get object meshes ...\n');
objectMeshes = dir(fullfile(folder, 'ObjectMeshes'));
objectMeshes = objectMeshes(~cellfun(@(x) strncmp(x, '.', 1), {objectMeshes.name}));

objectMeshesNumNodes = 0;
for ii=1:length(objectMeshes)
    tmpNodes=importdata(fullfile(folder, 'ObjectMeshes', objectMeshes(ii).name, 'Nodes.txt'),' ',1);
    tmpElements=importdata(fullfile(folder, 'ObjectMeshes', objectMeshes(ii).name, 'Elements.txt'),' ',1);
    objectMeshes(ii).nodes = tmpNodes.data;
    objectMeshes(ii).elements = tmpElements.data;
    objectMeshes(ii).num_nodes = size(tmpNodes.data, 1);
    objectMeshesNumNodes = objectMeshesNumNodes + objectMeshes(ii).num_nodes;
end

% build frequency vector
frequencies = params.frequencies;

clear ii tmpNodes tmpElements

%% Read computational effort
fprintf('Loading computational effort data ...\n');
% check for folder be.out and read computation time from every NC*.out file
for ii = 1:params.numSources
    computationTime{ii} = [];
    boundaryElements = dir([folder, filesep, 'NumCalc', filesep, 'source_', num2str(ii), filesep, 'NC*.out']);

    NC_all_flag = 0; NC_all_idx = []; NC_all_date = [];
    NC_single_flag = 0; NC_single_idx = []; NC_single_date = [];

    % read NC*.out configurations
    for jj = 1:size(boundaryElements, 1)
        if ~isempty(regexp(boundaryElements(jj).name, '^(NC.out)$')) % NC.out exists
            NC_all_flag = 1;
            NC_all_idx = jj;
            NC_all_date = datenum(boundaryElements(jj).date);
        elseif ~isempty(regexp(boundaryElements(jj).name, '(NC.)\S+(.out)')) % NC*.out exists
            NC_single_flag = 1;
            NC_single_idx = [NC_single_idx, jj];
            NC_single_date = [NC_single_date, datenum(boundaryElements(jj).date)];
        end
    end
    clear jj

    fprintf(['Reading computation time from source ', num2str(ii), '\n']);
    % read from .out files according to latest calculation
    if NC_all_flag && ~NC_single_flag % only NC.out exists
        % read computation time
        tmp=write_output_report([folder, filesep, 'NumCalc', filesep, 'source_', ...
            num2str(ii), filesep, boundaryElements(NC_all_idx).name]);
        computationTime{ii}=[computationTime{ii}; tmp];
    elseif NC_single_flag && ~NC_all_flag % only NC*.out exist
        % read single files
        for jj = 1:length(NC_single_idx)
            % read computation time
            tmp=write_output_report([folder, filesep, 'NumCalc', filesep, 'source_', ...
                num2str(ii), filesep, boundaryElements(NC_single_idx(jj)).name]);
            computationTime{ii}=[computationTime{ii}; tmp];
        end
    elseif (NC_all_flag && NC_single_flag) && any(NC_all_date < NC_single_date) % both exist, at least one NC*.out exists whose frequencies replace the ones in NC.out
        % start with NC.out ...
        tmp=write_output_report([folder, filesep, 'NumCalc', filesep, 'source_', ...
            num2str(ii), filesep, boundaryElements(NC_all_idx).name]);
        computationTime{ii}=[computationTime{ii}; tmp];
        % ... then read single files that are newer than NC.out
        for jj = 1:length(NC_single_idx)
            if NC_all_date < NC_single_date(jj)
                % read computation time
                tmp=write_output_report([folder, filesep, 'NumCalc', filesep, 'source_', ...
                    num2str(ii), filesep, boundaryElements(NC_single_idx(jj)).name]);
                computationTime{ii}(tmp(:,1),:) = tmp; % replace data only for current frequencies
            end
        end
    elseif (NC_all_flag && NC_single_flag) && all(NC_all_date > NC_single_date) % both exist, NC.out newest file
        % same as only NC.out exists
        % read computation time
        tmp=write_output_report([folder, filesep, 'NumCalc', filesep, 'source_', ...
            num2str(ii), filesep, boundaryElements(NC_all_idx).name]);
        computationTime{ii}=[computationTime{ii}; tmp];
    else % what possible case is this?
        error('This case is not yet implemented. Please open an issue at the project page: https://github.com/Any2HRTF/Mesh2HRTF/issues');
    end
    clear jj

    % check for error potential in numerical calculation
    warning('off', 'backtrace');
    iter_error_idx = find(computationTime{ii}(:,9));
    rel_error_idx = find(computationTime{ii}(:,7) > 1e-9);
    if ~isempty(iter_error_idx)
        for jj=1:length(iter_error_idx)
            warning(['Number of iterations for frequency ', num2str(computationTime{ii}(iter_error_idx(jj), 2)), ...
                ' Hz has reached the maximum.'], '\n');
        end
        clear jj
    end
    if ~isempty(rel_error_idx)
        for jj=1:length(rel_error_idx)
            warning(['Relative error for frequency ', num2str(computationTime{ii}(rel_error_idx(jj), 2)), ...
                ' Hz is greater than 1e-9.'], '\n');
        end
    end
    warning('on', 'backtrace');
    % sort computation time after frequency steps
    [~, sort_idx] = sort(computationTime{ii}(:,1));
    computationTime{ii} = computationTime{ii}(sort_idx,:);
end

fprintf('Write computation time to .mat file ...\n');
description={'Frequency index','Frequency','Building','Solving','Postprocessing','Total', 'relative error', 'number of iterations', 'maximum number of iterations reached'};
save(fullfile(folder, 'Output2HRTF', 'computationTime.mat'), 'description', 'computationTime', '-v6');
clear ii jj description computationTime tmp iter_error_idx rel_error_idx

%% Load ObjectMesh data
fprintf('Loading ObjectMesh data for ');
for ii=1:params.numSources
    fprintf(['source ', num2str(ii), '\n']);
    pressure{ii}=Output2HRTF_Load([folder, filesep, 'NumCalc', filesep, 'source_', num2str(ii), ...
        filesep, 'be.out', filesep], 'pBoundary', params.numFrequencies);
    clear data
end

clear ii

fprintf('Save ObjectMesh data ...\n');
cnt = 0;
for ii = 1:numel(objectMeshes)
    nodes=objectMeshes(ii).nodes;
    elements=objectMeshes(ii).elements;
    element_data=pressure;
    for jj = 1:numel(pressure)
        element_data{jj} = element_data{jj}(:, cnt+1:cnt+size(elements,1));
    end
    save(fullfile(folder, 'Output2HRTF', ['ObjectMesh_' objectMeshes(ii).name '.mat']), ...
        'nodes','elements','frequencies','element_data');
    cnt = cnt + size(elements,1);
end

clear pressure nodes elements ii jj cnt idx element_data

%% Load EvaluationGrid data
if ~isempty(evaluationGrids)
    fprintf('Loading evaluation grids for data in ');
    for ii=1:params.numSources
        fprintf(['source ', num2str(ii), '\n']);
        pressure(:,:,ii)=Output2HRTF_Load([folder, filesep, 'NumCalc', filesep, 'source_', num2str(ii), ...
            filesep, 'be.out', filesep], 'pEvalGrid', params.numFrequencies);
    end
end
clear ii

% save to struct
fprintf('Save EvaluationGrid data ...\n');
cnt = 0;
for ii = 1:numel(evaluationGrids)
    evaluationGrids(ii).pressure = pressure(:, cnt+1:cnt+evaluationGrids(ii).num_nodes, :);
    cnt = cnt + evaluationGrids(ii).num_nodes;
end


clear ii pressure cnt idx


% reference to pressure in the middle of the head with the head absent
% according to the HRTF definition.
if params.reference
    fprintf('Divide pressure by reference according to HRTF definition ...\n');
    % this might be a parameter in the function call
    refMode = 1;    % 1: reference to only one radius (the smallest found)
    % 2: reference to all indivudal radii

    for ii = 1:numel(evaluationGrids)

        xyz = evaluationGrids(ii).nodes;
        pressure = evaluationGrids(ii).pressure;
        freqMatrix = repmat(frequencies, [1 size(pressure,2) size(pressure,3)]);

        % distance of source positions from the origin
        if refMode == 1
            r = min(sqrt(xyz(:,2).^2 + xyz(:,3).^2 + xyz(:,4).^2));
            r = repmat(r, [size(pressure,1) size(pressure, 2) size(pressure, 3)]);
        else
            r = sqrt(xyz(:,2).^2 + xyz(:,3).^2 + xyz(:,4).^2);
            r = repmat(r', [size(pressure,1) 1 size(pressure, 3)]);
        end

        if strcmp(params.sourceType, 'Left ear') || ...
                strcmp(params.sourceType, 'Right ear') || ...
                strcmp(params.sourceType, 'Both ears')

            volumeFlow = .1 * ones(size(pressure));
            if isfield(params, 'sourceArea')
                for jj = 1:numel(params.sourceArea)
                    volumeFlow(:,:,jj) = volumeFlow(:,:,jj) * params.sourceArea(jj);
                end
            end

            % point source in the origin evaluated at r
            % eq. (6.71) in: Williams, E. G. (1999). Fourier Acoustics.
            ps   = -1j * params.densityOfMedium * 2*pi*freqMatrix .* volumeFlow ./ (4*pi) .* ...
                exp(1j * 2*pi*freqMatrix/params.speedOfSound .* r) ./ r;

        elseif strcmp(params.sourceType, 'Point source')

            amplitude = .1; % hard coded in Mesh2HRTF
            ps = amplitude * exp(1j * 2*pi*freqMatrix/params.speedOfSound .*r) ./ (4 * pi * r);

        elseif strcmp(params.sourceType, 'Plane wave')
            error('Referencing for plane wave source type not yet implemented.\n');
        else
            error('Referencing is currently only implemented for sourceType ''vibratingElement'' and ''pointSource''.\n')
        end

        evaluationGrids(ii).pressure = pressure ./ ps;

    end

    clear r freqMatrix ps ii jj
end

%% Save data as SOFA file
fprintf('Saving complex pressure to SOFA file ...\n')
SOFAstart;

for ii = 1:numel(evaluationGrids)

    xyz = evaluationGrids(ii).nodes;

    % check if .sofa file already exists
    filename = fullfile('Output2HRTF', ['HRTF_' evaluationGrids(ii).name '.sofa']);
    if exist(filename, 'file')
        % if yes, delete and write new object, but ask user
        prompt = [filename, ' already exists. Replace object? [y/n]\n'];
        replace_obj = input(prompt, 's');
        if strcmp(replace_obj, 'y') || strcmp(replace_obj, 'yes')
            delete(filename)
        elseif strcmp(replace_obj, 'n') || strcmp(replace_obj, 'no')
            error('Object was not replaced. Output2HRTF aborted.\n');
        else
            error('Invalid input. Output2HRTF aborted.\n');
        end
    end

    % get SOFA template according to number of sources
    if params.numSources == 2
        Obj = SOFAgetConventions('SimpleFreeFieldHRTF');
    else
        % Save as GeneralTF
        Obj = SOFAgetConventions('GeneralTF');
    end

    % pressure has size N (frequencies) x M (measurements) x R (sources)
    pressure = evaluationGrids(ii).pressure;
    % get dimenson 3 explicitly, because it would be dropped if R=1
    NMR = [size(pressure, 1), size(pressure, 2), size(pressure, 3)];

    % shift dimensions to MRN as expected by SOFA
    % (shifting loses leading and trailing dimensons of size 1)
    pressure = shiftdim(pressure, 1);
    % force dimensions of 1
    pressure = reshape(pressure, NMR(2), NMR(3), NMR(1));

    Obj.GLOBAL_ApplicationName = 'Mesh2HRTF';
    Obj.GLOBAL_ApplicationVersion = params.Mesh2HRTF_Version;
    Obj.GLOBAL_History = 'numerically calculated data';

    Obj.N=frequencies;
    Obj.Data.Real=real(pressure);
    Obj.Data.Imag=imag(pressure);
    Obj.Data.Real_LongName='pressure';
    Obj.Data.Real_Units='pascal';
    Obj.Data.Imag_LongName='pressure';
    Obj.Data.Imag_Units='pascal';

    Obj.ListenerPosition = [0 0 0];
    Obj.SourcePosition=xyz(:,2:4);
    Obj.SourcePosition_Type='cartesian';
    Obj.SourcePosition_Units='metre';
    Obj.ReceiverPosition=params.sourceCenter;
    Obj.ReceiverPosition_Type='cartesian';
    Obj.ReceiverPosition_Units='metre';

    Obj=SOFAupdateDimensions(Obj);
    SOFAsave(fullfile(folder, 'Output2HRTF', ['HRTF_' evaluationGrids(ii).name '.sofa']),Obj);
    fprintf(['HRTF_' evaluationGrids(ii).name '.sofa saved!\n']);
end

clear Obj ii xyz pressure prompt replace_Obj


%% Save time data data as SOFA file
if params.computeHRIRs
    fprintf('Saving HRIR data to SOFA file ...\n')
    for ii = 1:numel(evaluationGrids)
        % check if the frequency vector has the correct format
        if any(abs(diff(frequencies,2)) > .1) || frequencies(1) < .1
            error('The frequency vector must start at a frequency > 0.1 and continue in equidistant steps to the end.\n')
        end

        % check if reference exists
        if ~params.reference
            error('HRIRs can only be computed if refernce=true\n')
        end

        xyz = evaluationGrids(ii).nodes;
        pressure = evaluationGrids(ii).pressure;

        % check if .sofa file already exists
        filename = fullfile('Output2HRTF', ['HRIR_' evaluationGrids(ii).name '.sofa']);
        if exist(filename, 'file')
            % if yes, delete and write new object, but ask user
            prompt = [filename, ' already exists. Replace object? [y/n]\n'];
            replace_obj = input(prompt, 's');
            if strcmp(replace_obj, 'y') || strcmp(replace_obj, 'yes')
                delete(filename)
            elseif strcmp(replace_obj, 'n') || strcmp(replace_obj, 'no')
                error('Object was not replaced. Output2HRTF aborted.\n');
            else
                error('Invalid input. Output2HRTF aborted.\n');
            end
        end
        clear prompt

%         prompt = ['The default sampling frequency is twice the highest occuring frequency. ', ...
%             'In this case fs = ', num2str(2*frequencies(end)), '. \n', ...
%             'Do you want to specify a different sampling frequency? [y/n]\n'];
%         replace_fs = input(prompt, 's');
%         if strcmp(replace_fs, 'y') || strcmp(replace_fs, 'yes')
%             clear prompt
%             prompt = 'Please specify the sampling frequency in Hz: ';
%             fs = input(prompt);
%             if fs < 0
%                 error('The sampling frequency has to be positive.\n');
%             elseif fs < frequencies(end)
%                 warning('The sampling frequency is lower than the highest frequency in the signal. Further processing will potentially lead to aliasing.');
%             end
%             clear prompt
%         elseif strcmp(replace_fs, 'n') || strcmp(replace_fs, 'no')
%             fs = 2*frequencies(end);
%         else
%             error('Invalid input. Writing HRIR to SOFA object aborted.\n');
%         end

        fs = 2*frequencies(end);

        % add 0 Hz bin
        pressure = [ones(1, size(pressure,2), size(pressure,3));
            pressure];
        % make fs/2 real
        pressure(end,:,:) = abs(pressure(end,:,:));
        % mirror the spectrum
        pressure = [pressure; flipud(conj(pressure(2:end-1,:,:)))];
        % ifft (take complex conjugate because sign conventions differ)
        hrir = ifft(conj(pressure), 'symmetric');

        % shift 30 cm to make causal
        % (path differences between the origin and the ear are usually
        % smaller than 30 cm but numerical HRIRs show stringer pre-ringing)
        n_shift = round(.30 / (1/fs * params.speedOfSound));
        hrir = circshift(hrir, n_shift);

        % hrir has size N (samples) x M (measurements) x R (sources)
        % get dimenson 3 explicitly, because it would be dropped if R=1
        NMR = [size(hrir, 1), size(hrir, 2), size(hrir, 3)];

        % shift dimensions to MRN as expected by SOFA
        % (shifting loses leading and trailing dimensons of size 1)
        hrir = shiftdim(hrir, 1);
        % force dimensions of 1
        hrir = reshape(hrir, NMR(2), NMR(3), NMR(1));

        % get SOFA template according to number of sources
        if params.numSources == 2
            Obj = SOFAgetConventions('SimpleFreeFieldHRIR');
        else
            % Save as GeneralFIR
            Obj = SOFAgetConventions('GeneralFIR');
        end

        Obj.GLOBAL_ApplicationName = 'Mesh2HRTF';
        Obj.GLOBAL_ApplicationVersion = params.Mesh2HRTF_Version;
        Obj.GLOBAL_History = 'numerically calculated data';

        Obj.Data.IR=hrir;
        Obj.Data.SamplingRate=fs;
        Obj.Data.Delay=zeros(1, size(hrir, 2));

        Obj.ListenerPosition = [0 0 0];
        Obj.SourcePosition=xyz(:,2:4);
        Obj.SourcePosition_Type='cartesian';
        Obj.SourcePosition_Units='metre';
        Obj.ReceiverPosition=params.sourceCenter;
        Obj.ReceiverPosition_Type='cartesian';
        Obj.ReceiverPosition_Units='metre';

        Obj=SOFAupdateDimensions(Obj);
        SOFAsave(fullfile(folder, 'Output2HRTF', ['HRIR_' evaluationGrids(ii).name '.sofa']),Obj);
        fprintf(['HRIR_' evaluationGrids(ii).name '.sofa saved!\n']);
    end
end

clear Obj ii xyz pressure prompt replace_fs

fprintf('Done\n')

end % end of main function ------------------------------------------------


function data = Output2HRTF_Load(foldername,filename, numFrequencies)
%   data = OUTPUT2HRTF_LOAD(foldername,filename)
%
%   Loads results of the BEM-HRTF calculation from the NumCalc folder.
%
%   Input:
%       foldername ... path to NumCalc output
%       filename ..... either 'pBoundary' or 'pEvalGrid'
%
%   Output:
%       data ......... Matrix of complex values
%           dim1: frequencies
%           dim2: datapoints (complex pressure values)

%% --------------------check number of header lines------------------------
for ii=1:1000
    temp1=importdata([foldername, 'be.1', filesep, filename], ' ', ii);
    if ~isempty(temp1)
        if size(temp1.data,2)==3
            if temp1.data(2,1)-temp1.data(1,1)==1
                numHeaderlines_BE=ii;
                break
            end
        end
    end
end

%% -----------------------------load data----------------------------------
temp1 = importdata([foldername, 'be.1', filesep, filename], ' ', numHeaderlines_BE);

data=zeros(numFrequencies,size(temp1.data,1));

for ii=1:numFrequencies
    tmpData = importdata([foldername, 'be.', num2str(ii), filesep, filename], ' ', numHeaderlines_BE);
    if ~isempty(tmpData)
        data(ii,:)=transpose(tmpData.data(:,2)+1i*tmpData.data(:,3));
    else
        data(ii,:)      = NaN;
    end
end

end %of function
