%                                Mesh2HRTF
%                Copyright (C) 2015 by Harald Ziegelwanger,
%        Acoustics Research Institute, Austrian Academy of Sciences
%                        mesh2hrtf.sourceforge.net
% 
% Mesh2HRTF is licensed under the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
% Mesh2HRTF is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.
% You should have received a copy of the GNU LesserGeneral Public License along with Mesh2HRTF. If not, see <http://www.gnu.org/licenses/lgpl.html>.
% 
% If you use Mesh2HRTF:
% - Provide credits:
%   "Mesh2HRTF, H. Ziegelwanger, ARI, OEAW (mesh2hrtf.sourceforge.net)"
% - In your publication, cite both articles:
%   [1] Ziegelwanger, H., Kreuzer, W., and Majdak, P. (2015). "Mesh2HRTF: Open-source software package for the numerical calculation of head-related transfer functions," in Proceedings of the 22nd ICSV, Florence, IT.
%   [2] Ziegelwanger, H., Majdak, P., and Kreuzer, W. (2015). "Numerical calculation of listener-specific head-related transfer functions and sound localization: Microphone model and mesh discretization," The Journal of the Acoustical Society of America, 138, 208-222.
%
% Author: Harald Ziegelwanger (Acoustics Research Institute, Austrian Academy of Sciences)
% Co-Authors: Fabian Brinkmann, Robert Pelzer (Audio Communication Group, Technical University Berlin)	

function Output2HRTF_Main(Mesh2HRTF_version, cpusAndCores, ...
                          sourceType, sourceCenter, sourceArea, ...
                          reference, computeHRIRs, ...
                          speedOfSound, densityOfAir)
%OUTPUT2HRTF_MAIN
%   []=Output2HRTF_Main(cpusAndCores,objectMeshesUsed,reciprocity,
%   receiverPositions,receiverArea,reference,speedOfSound,densityOfAir) calculates
%   all relevant data after the NumCalc calculation and saves the results.
%
%   Input:
%       cpusAndCores:
% 
%       reciprocity:
%
%       receiverPositions:
%
%       microphoneArea: are of the mesh elements that were used as sound
%       source in reciprocal calculation in square meters
%
%       reference: 'true' complex pressure at the evaluation grid is
%       devided by a point source in the origin of the coordinate system.
%       This is the calssical HRTF definition (pressure st the ear divided
%       by pressure at the center of the head with the head being absent)

%% ----------------------------load meta data------------------------------
% output directory
if ~exist(fullfile(pwd, 'Output2HRTF'), 'dir')
    mkdir(fullfile(pwd, 'Output2HRTF'))
end

% number of ears
ears=max(unique(cpusAndCores));

% get the evaluation grids
evaluationGrids = dir('EvaluationGrids');
evaluationGrids = evaluationGrids(~ismember({evaluationGrids.name},{'.','..'}));

evaluationGridsNumNodes = 0;
for ii=1:length(evaluationGrids)
    tmpNodes=importdata(fullfile('EvaluationGrids', evaluationGrids(ii).name, 'Nodes.txt'),' ',1);
    tmpElements=importdata(fullfile('EvaluationGrids', evaluationGrids(ii).name, 'Elements.txt'),' ',1);
    evaluationGrids(ii).nodes = tmpNodes.data;
    evaluationGrids(ii).elements = tmpElements.data;
    evaluationGrids(ii).num_nodes = size(tmpNodes.data, 1);
    evaluationGridsNumNodes = evaluationGridsNumNodes + evaluationGrids(ii).num_nodes;
end

% get the object mesh
objectMeshes = dir('ObjectMeshes');
objectMeshes = objectMeshes(~ismember({objectMeshes.name},{'.','..'}));

% evaluationGridNodes=[];
% evaluationGridElements=[];

objectMeshesNumNodes = 0;
for ii=1:length(objectMeshes)
    tmpNodes=importdata(fullfile('ObjectMeshes', objectMeshes(ii).name, 'Nodes.txt'),' ',1);
    tmpElements=importdata(fullfile('ObjectMeshes', objectMeshes(ii).name, 'Elements.txt'),' ',1);
    objectMeshes(ii).nodes = tmpNodes.data;
    objectMeshes(ii).elements = tmpElements.data;
    objectMeshes(ii).num_nodes = size(tmpNodes.data, 1);
    objectMeshesNumNodes = objectMeshesNumNodes + objectMeshes(ii).num_nodes;
end

clear ii tmpNodes tmpElements

%% Read computational effort
fprintf('\nLoading computational effort data ...');
for ch=1:ears
    computationTime{ch}=[];
    fprintf(['\n    Ear ' num2str(ch) ' ...'])
    for ii=1:size(cpusAndCores,1)
        fprintf(['\n        CPU ' num2str(ii) ': ']);
        for jj=1:size(cpusAndCores,2)
            if cpusAndCores(ii,jj)==ch
                fprintf([num2str(jj) ', ']);
                tmp=Output2HRTF_ReadComputationTime(['NumCalc' filesep 'CPU_' num2str(ii) '_Core_' num2str(jj) filesep 'NC.out']);
                computationTime{ch}=[computationTime{ch}; tmp];
            end
            clear tmp
        end
        fprintf('...');
    end
end

description={'Frequency index','Frequency','Building','Solving','Postprocessing','Total'};
save(fullfile('Output2HRTF', 'computationTime.mat'), 'description', 'computationTime', '-v6');
clear ch ii jj description computationTime

%% Load ObjectMesh data
fprintf('\nLoading ObjectMesh data ...');
for ch=1:ears
    fprintf(['\n    Ear ' num2str(ch) ' ...'])
    for ii=1:size(cpusAndCores,1)
        fprintf(['\n        CPU ' num2str(ii) ': ']);
        for jj=1:size(cpusAndCores,2)
            if cpusAndCores(ii,jj)==ch
                fprintf([num2str(jj) ', ']);
                [tmpData,tmpFrequencies]=Output2HRTF_Load(['NumCalc' filesep 'CPU_' num2str(ii) '_Core_' num2str(jj) filesep 'be.out' filesep],'pBoundary');
                if exist('tmpPressure','var')
                    tmpPressure=[tmpPressure; tmpData];
                    frequencies=[frequencies; tmpFrequencies];
                else
                    tmpPressure=tmpData;
                    frequencies=tmpFrequencies;
                end
                clear tmpData tmpFrequencies
            end
        end
        fprintf('...');
    end
    [frequencies,idx]=sort(frequencies);
    pressure{ch}=tmpPressure(idx,:);
    
    clear tmpPressure tmpPhase
end

fprintf('\nSave ObjectMesh data ...');
cnt = 0;
for ii = 1:numel(objectMeshes)
    nodes=objectMeshes(ii).nodes;
    elements=objectMeshes(ii).elements;
    element_data=pressure;
    for jj = 1:numel(pressure)
        element_data{jj} = element_data{jj}(:, cnt+1:cnt+size(elements,1));
    end
    save(fullfile('Output2HRTF', ['ObjectMesh_' objectMeshes(ii).name '.mat']), ...
        'nodes','elements','frequencies','element_data');
    cnt = cnt + size(elements,1);
end

clear pressure nodes elements frequencies ii jj cnt ch idx element_data

%% Load EvaluationGrid data
if ~isempty(evaluationGrids)
    fprintf('\nLoading data for the evaluation grids ...');
    for ch=1:ears
        fprintf(['\n    Ear ' num2str(ch) ' ...'])
        for ii=1:size(cpusAndCores,1)
            fprintf(['\n        CPU ' num2str(ii) ': ']);
            for jj=1:size(cpusAndCores,2)
                if cpusAndCores(ii,jj)==ch
                    fprintf([num2str(jj) ', ']);
                    [tmpData,tmpFrequencies]=Output2HRTF_Load(['NumCalc' filesep 'CPU_' num2str(ii) '_Core_' num2str(jj) filesep 'be.out' filesep],'pEvalGrid');
                    if exist('tmpPressure','var')
                        tmpPressure=[tmpPressure; tmpData];
                        frequencies=[frequencies; tmpFrequencies];
                    else
                        tmpPressure=tmpData;
                        frequencies=tmpFrequencies;
                    end
                    clear tmpData tmpFrequencies
                end
            end
            fprintf('...');
        end
        pressure(:,:,ch)=tmpPressure;
        clear tmpPressure
    end
    [frequencies,idx]=sort(frequencies);
    pressure=pressure(idx,:,:);
end

% save to struct
cnt = 0;
for ii = 1:numel(evaluationGrids)
    evaluationGrids(ii).pressure = pressure(:, cnt+1:cnt+evaluationGrids(ii).num_nodes, :);
    cnt = cnt + evaluationGrids(ii).num_nodes;
end


clear ch ii jj tmpData tmpFrequencies tmpPressure pressure cnt idx


% reference to pressure in the middle of the head with the head absent
% according to the HRTF definition.
if reference
    
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
        
        if strcmp(sourceType, 'vibratingElement')
            
            volumeFlow = .1 * ones(size(pressure));
            if exist('sourceArea', 'var')
                %has to be fixed for both ears....
                for nn = 1:numel(sourceArea)
                    volumeFlow(:,:,nn) = volumeFlow(:,:,nn) * sourceArea(nn);
                end
            end
            
            % point source in the origin evaluated at r
            % eq. (6.71) in: Williams, E. G. (1999). Fourier Acoustics.
            ps   = -1j * densityOfAir * 2*pi*freqMatrix .* volumeFlow ./ (4*pi) .* ...
                exp(1j * 2*pi*freqMatrix/speedOfSound .* r) ./ r;
            
        elseif strcmp(sourceType, 'pointSource')
            
            amplitude = .1; % hard coded in Mesh2HRTF
            ps = amplitude * exp(1j * 2*pi*freqMatrix/speedOfSound .*r) ./ (4 * pi * r);
            
        else
            error('Referencing is currently only implemented for sourceType ''vibratingElement'' and ''pointSource''.')
        end
        
        % here we go...
        evaluationGrids(ii).pressure = pressure ./ ps;
        
    end
    
    clear Areceiver r freqMatrix ps ii
end

%% Save data as SOFA file
fprintf('\nSaving complex pressure to SOFA file ...\n')
SOFAstart;

for ii = 1:numel(evaluationGrids)
    
    xyz = evaluationGrids(ii).nodes;
    pressure = evaluationGrids(ii).pressure;
    
    % prepare pressure to be size of MRN when R=2
    pressure = shiftdim(pressure, 1);
    % prepare pressure to be size of MRN when R=1
    if size(pressure, 3) == 1
        pressure = reshape(pressure, size(pressure,1), 1, size(pressure, 2));
    end
    
    % Save as GeneralTF
    Obj = SOFAgetConventions('GeneralTF');
    
    Obj.GLOBAL_ApplicationName = 'Mesh2HRTF';
    Obj.GLOBAL_ApplicationVersion = Mesh2HRTF_version;
    Obj.GLOBAL_Organization = '';
    Obj.GLOBAL_Title = '';
    Obj.GLOBAL_DateCreated = date;
    Obj.GLOBAL_AuthorContact = '';
    
    Obj.ReceiverPosition=sourceCenter;
    Obj.N=frequencies;
    
    Obj.Data.Real=real(pressure);
    Obj.Data.Imag=imag(pressure);
    Obj.Data.Real_LongName='pressure';
    Obj.Data.Real_Units='pascal';
    Obj.Data.Imag_LongName='pressure';
    Obj.Data.Imag_Units='pascal';
    
    Obj.ListenerPosition = [0 0 0];
    Obj.SourcePosition_Type='cartesian';
    Obj.SourcePosition_Units='meter';
    Obj.SourcePosition=xyz(:,2:4);
    
    Obj=SOFAupdateDimensions(Obj);
    SOFAsave(fullfile('Output2HRTF', ['HRTF_' evaluationGrids(ii).name '.sofa']),Obj);
end

clear Obj ii xyz pressure


%% Save time data data as SOFA file
if computeHRIRs
    
    fprintf('\nSaving time data to SOFA file ...\n')
    
    for ii = 1:numel(evaluationGrids)
        
        % check if the frequency vector has the correct format
        if ~all(abs(frequencies(1) - diff(frequencies)) < .1)
            error('The frequency vector must be if the format a:a:fs/2, with a>0 and fs the sampling rate.')
        end
        
        if ~reference
            error('HRIRs can only be computet if refernce=true')
        end
        
        xyz = evaluationGrids(ii).nodes;
        pressure = evaluationGrids(ii).pressure;
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
        n_shift = round(.30 / (1/fs * speedOfSound));
        hrir = circshift(hrir, n_shift);
        
        
        % prepare pressure to be size of MRN when R=2
        hrir = shiftdim(hrir, 1);
        % prepare pressure to be size of MRN when R=1
        if size(hrir, 3) == 1
            hrir = reshape(hrir, size(hrir,1), 1, size(hrir, 2));
        end
        
        % Save as GeneralFIR
        Obj = SOFAgetConventions('GeneralFIR', 'm');
        
        Obj.GLOBAL_ApplicationName = 'Mesh2HRTF';
        Obj.GLOBAL_ApplicationVersion = Mesh2HRTF_version;
        Obj.GLOBAL_Organization = '';
        Obj.GLOBAL_Title = '';
        Obj.GLOBAL_DateCreated = date;
        Obj.GLOBAL_AuthorContact = '';
        
        Obj.ReceiverPosition=sourceCenter;
        
        Obj.Data.IR=hrir;
        Obj.Data.SamplingRate=fs;
        
        Obj.ListenerPosition = [0 0 0];
        Obj.SourcePosition_Type='cartesian';
        Obj.SourcePosition_Units='meter';
        Obj.SourcePosition=xyz(:,2:4);
        
        Obj=SOFAupdateDimensions(Obj);
        SOFAsave(fullfile('Output2HRTF', ['HRIR_' evaluationGrids(ii).name '.sofa']),Obj);
    end
end

clear Obj ii xyz pressure

fprintf('Done\n')