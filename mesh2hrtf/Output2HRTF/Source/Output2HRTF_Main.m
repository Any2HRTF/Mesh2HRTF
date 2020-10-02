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

function Output2HRTF_Main(Mesh2HRTF_version, cpusAndCores,receiverCenter,receiverArea,reference,speedOfSound,densityOfAir)
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
if ~exist('Output2HRTF', 'dir')
    mkdir('Output2HRTF')
end

% number of ears
ears=max(unique(cpusAndCores));

temp=dir('EvaluationGrids');
evaluationGrids={};
for ii=3:length(temp)
    evaluationGrids{ii-2}=temp(ii).name;
end
evaluationGridNodes=[];
evaluationGridElements=[];
for ii=1:length(evaluationGrids)
    if exist(['EvaluationGrids' filesep evaluationGrids{ii} filesep 'Nodes.txt'],'file')
        tmpNodes=importdata(['EvaluationGrids' filesep evaluationGrids{ii} filesep 'Nodes.txt'],' ',1);
        tmpElements=importdata(['EvaluationGrids' filesep evaluationGrids{ii} filesep 'Elements.txt'],' ',1);
        evaluationGridNodes=[evaluationGridNodes; tmpNodes.data];
        evaluationGridElements=[evaluationGridElements; tmpElements.data(:,1:end-3)];
    end
end

temp=dir('ObjectMeshes');
objectMeshes={};
for ii=3:length(temp)
    objectMeshes{ii-2}=temp(ii).name;
end
objectMeshNodes={};
objectMeshElements={};
objectMeshMaxFrequency={};
for ii=1:length(objectMeshes)
    if exist(['ObjectMeshes' filesep objectMeshes{ii} filesep 'Nodes.txt'],'file')
        tmpNodes=importdata(['ObjectMeshes' filesep objectMeshes{ii} filesep 'Nodes.txt'],' ',1);
        tmpElements=importdata(['ObjectMeshes' filesep objectMeshes{ii} filesep 'Elements.txt'],' ',1);
        objectMeshNodes{1}=tmpNodes.data;
        objectMeshElements{1}=tmpElements.data(:,1:end-3);
    end
end

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
clear description computationTime

%% Calculate ObjectMesh data
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
element_data=pressure;
nodes=objectMeshNodes;
elements=objectMeshElements;
save(fullfile('Output2HRTF', 'ObjectMesh.mat'),'nodes','elements','frequencies','element_data');

clear pressure nodes elements frequencies

%% Calculate EvaluationGrid data
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

% added Fabian Brinkmann
if reference
    
    % this might be a parameter in the function call
    refMode = 1;    % 1: reference to only one radius (the smallest found)
                    % 2: reference to all indivudal radii
    
    % reference to pressure in the middle of the head with the head absent
    % (HRTF definition). We do it reciprocal for ease of computation.
    
    volumeFlow = .1 * ones(size(pressure));
    if exist('receiverArea', 'var')
        %has to be fixed for both ears....
        for nn = 1:numel(receiverArea)
            volumeFlow(:,:,nn) = volumeFlow(:,:,nn) * receiverArea(nn);
        end       
    end
    
    % distance of source positions from the origin
    if refMode == 1
        r = min(sqrt(evaluationGridNodes(:,2).^2 + evaluationGridNodes(:,3).^2 + evaluationGridNodes(:,4).^2));
        r = repmat(r, [size(pressure,1) size(pressure, 2) size(pressure, 3)]);
    else
        r = sqrt(evaluationGridNodes(:,2).^2 + evaluationGridNodes(:,3).^2 + evaluationGridNodes(:,4).^2);
        r = repmat(r', [size(pressure,1) 1 size(pressure, 3)]);
    end

    % point source in the origin evaluated at r
    % eq. (6.71) in: Williams, E. G. (1999). Fourier Acoustics.
    freqMatrix = repmat(frequencies, [1 size(pressure,2) size(pressure,3)]);
    ps   = -1j * densityOfAir * 2*pi*freqMatrix .* volumeFlow ./ (4*pi) .* ...
           exp(1j * 2*pi*freqMatrix/speedOfSound .* r) ./ r;
    % here we go...
    pressure = pressure ./ ps;
    
    clear Areceiver r freqMatrix ps
end
% end of added Fabian Brinkmann

%% Save data as SOFA file
fprintf('\nSaving complex pressure to SOFA file ...\n')
SOFAstart;

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

Obj.ReceiverPosition=receiverCenter;
Obj.N=frequencies;

%add division G([0,0,0],evaluationGrid)
Obj.Data.Real=real(pressure);
Obj.Data.Imag=imag(pressure);
Obj.Data.Real_LongName='pressure';
Obj.Data.Real_Units='pascal';
Obj.Data.Imag_LongName='pressure';
Obj.Data.Imag_Units='pascal';

Obj.ListenerPosition = [0 0 0];
Obj.SourcePosition_Type='cartesian';
Obj.SourcePosition_Units='meter';
Obj.SourcePosition=evaluationGridNodes(:,2:4);

Obj=SOFAupdateDimensions(Obj);
SOFAsave(fullfile('Output2HRTF', 'EvaluationGrid_GeneralTF.sofa'),Obj);

fprintf('Done\n')