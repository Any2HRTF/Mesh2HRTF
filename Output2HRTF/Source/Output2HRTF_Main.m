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

function Output2HRTF_Main(cpusAndCores,objectMeshesUsed,reciprocity,receiverPositions,frequencyDependency,nearFieldCalculation)
%OUTPUT2HRTF_MAIN
%   []=Output2HRTF_Main(cpusAndCores,objectMeshesUsed,reciprocity,
%   receiverPositions,frequencyDependency,nearFieldCalculation) calculates
%   all relevant data after the NumCalc calculation and saves the results.
%
%   Input:
%       cpusAndCores:
%
%       objectMeshesUsed:
% 
%       reciprocity:
%
%       receiverPositions:
%
%       frequencyDependency:
% 
%       nearFieldCalculation:

ears=max(unique(cpusAndCores));

%% ----------------------------load meta data------------------------------
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
countL=0; countR=0; count=0;
for ii=1:length(objectMeshes)
    if exist(['ObjectMeshes' filesep objectMeshes{ii} filesep 'Nodes.txt'],'file')
        switch objectMeshes{ii}(1)
            case 'L'
                idxEar=1;
                countL=countL+1;
                count=countL;
            case 'R'
                switch objectMeshes{ii}(2)
                    case 'e'
                        idxEar=3;
                        count=1;
                    otherwise
                        idxEar=2;
                        countR=countR+1;
                        count=countR;
                    end
            otherwise
                idxEar=0;
        end
        if frequencyDependency==0 && objectMeshes{ii}(1)~='R'
            idxEar=0;
        end
        if idxEar~=0
            tmpNodes=importdata(['ObjectMeshes' filesep objectMeshes{ii} filesep 'Nodes.txt'],' ',1);
            tmpElements=importdata(['ObjectMeshes' filesep objectMeshes{ii} filesep 'Elements.txt'],' ',1);
            objectMeshNodes{count,idxEar}=tmpNodes.data;
            objectMeshElements{count,idxEar}=tmpElements.data(:,1:end-3);
            if idxEar~=3
                objectMeshMaxFrequency{count,idxEar}=str2double(objectMeshes{ii}(2:end));
            end
        end
    end
end
if isempty(objectMeshNodes{1,1}) && isempty(objectMeshElements{1,1}) && isempty(objectMeshNodes{1,2}) && isempty(objectMeshElements{1,2})
    objectMeshNodes{1,1}=objectMeshNodes{1,3};
    objectMeshElements{1,1}=objectMeshElements{1,3};
    objectMeshNodes{1,2}=objectMeshNodes{1,3};
    objectMeshElements{1,2}=objectMeshElements{1,3};
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
save('computationTime.mat','description','computationTime','-v6');
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
save('ObjectMesh.mat','nodes','elements','frequencies','element_data');

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

%% save EvalulaitonGrid data for visualization
fprintf('\nSave EvaluationGrid data ...');
if ~isempty(evaluationGrids)
    %save results in the sagittal plane
    if sum(evaluationGridNodes(:,1)>=700000 & evaluationGridNodes(:,1)<=799999)~=0
        idx=evaluationGridElements(:,1)>=700000 & evaluationGridElements(:,1)<=799999;
        elements=evaluationGridElements(idx,:)-min(min(evaluationGridElements(idx,:)));
        idx=evaluationGridNodes(:,1)>=700000 & evaluationGridNodes(:,1)<=799999;
        nodes=evaluationGridNodes(idx,:);
        nodes(:,1)=nodes(:,1)-min(nodes(:,1));
        node_data=pressure(:,idx,:);
        save('EvaluationGrid_SAGPLANE.mat','nodes','elements','frequencies','node_data');
        clear nodes elements idx node_data
    end

    %save results in the horizontal plane
    if sum(evaluationGridNodes(:,1)>=800000 & evaluationGridNodes(:,1)<=899999)~=0
        idx=evaluationGridElements(:,1)>=800000 & evaluationGridElements(:,1)<=899999;
        elements=evaluationGridElements(idx,:)-min(min(evaluationGridElements(idx,:)));
        idx=evaluationGridNodes(:,1)>=800000 & evaluationGridNodes(:,1)<=899999;
        nodes=evaluationGridNodes(idx,:);
        nodes(:,1)=nodes(:,1)-min(nodes(:,1));
        node_data=pressure(:,idx,:);
        save('EvaluationGrid_HORPLANE.mat','nodes','elements','frequencies','node_data');
        clear nodes elements idx node_data
    end

    %save results in the frontal plane
    if sum(evaluationGridNodes(:,1)>=900000 & evaluationGridNodes(:,1)<=999999)~=0
        idx=evaluationGridElements(:,1)>=900000 & evaluationGridElements(:,1)<=999999;
        elements=evaluationGridElements(idx,:)-min(min(evaluationGridElements(idx,:)));
        idx=evaluationGridNodes(:,1)>=900000 & evaluationGridNodes(:,1)<=999999;
        nodes=evaluationGridNodes(idx,:);
        nodes(:,1)=nodes(:,1)-min(nodes(:,1));
        node_data=pressure(:,idx,:);
        save('EvaluationGrid_FRTPLANE.mat','nodes','elements','frequencies','node_data');
        clear nodes elements idx node_data
    end
end
    
if ~isempty(evaluationGrids)
    if sum(evaluationGridNodes(:,1)>=210000 & evaluationGridNodes(:,1)<=219999)~=0
        idx=evaluationGridElements(:,1)>=210000 & evaluationGridElements(:,1)<=219999;
        elements=evaluationGridElements(idx,:)-min(min(evaluationGridElements(idx,:)));
        idx=evaluationGridNodes(:,1)>=210000 & evaluationGridNodes(:,1)<=219999;
        nodes=evaluationGridNodes(idx,:);
        nodes(:,1)=nodes(:,1)-min(nodes(:,1));
        node_data=pressure(:,idx,:);
        save('EvaluationGrid_NF.mat','nodes','elements','frequencies','node_data');
        clear nodes elements idx node_data
    end
    if sum(evaluationGridNodes(:,1)>=220000 & evaluationGridNodes(:,1)<=229999)~=0
        idx=evaluationGridElements(:,1)>=220000 & evaluationGridElements(:,1)<=229999;
        elements=evaluationGridElements(idx,:)-min(min(evaluationGridElements(idx,:)));
        idx=evaluationGridNodes(:,1)>=220000 & evaluationGridNodes(:,1)<=229999;
        nodes=evaluationGridNodes(idx,:);
        nodes(:,1)=nodes(:,1)-min(nodes(:,1));
        node_data=pressure(:,idx,:);
        save('EvaluationGrid_FF.mat','nodes','elements','frequencies','node_data');
        clear nodes elements idx node_data
    end
    if sum(evaluationGridNodes(:,1)>=300000 & evaluationGridNodes(:,1)<=349999)~=0
        idx=evaluationGridElements(:,1)>=300000 & evaluationGridElements(:,1)<=349999;
        elements=evaluationGridElements(idx,:)-min(min(evaluationGridElements(idx,:)));
        idx=evaluationGridNodes(:,1)>=300000 & evaluationGridNodes(:,1)<=349999;
        nodes=evaluationGridNodes(idx,:);
        nodes(:,1)=nodes(:,1)-min(nodes(:,1));
        node_data=pressure(:,idx,:);
        save('EvaluationGrid_ARI.mat','nodes','elements','frequencies','node_data');
        clear nodes elements idx node_data
    end
    if sum(evaluationGridNodes(:,1)>=350000 & evaluationGridNodes(:,1)<=399999)~=0
        idx=evaluationGridElements(:,1)>=350000 & evaluationGridElements(:,1)<=399999;
        elements=evaluationGridElements(idx,:)-min(min(evaluationGridElements(idx,:)));
        idx=evaluationGridNodes(:,1)>=350000 & evaluationGridNodes(:,1)<=399999;
        nodes=evaluationGridNodes(idx,:);
        nodes(:,1)=nodes(:,1)-min(nodes(:,1));
        node_data=pressure(:,idx,:);
        save('EvaluationGrid_User.mat','nodes','elements','frequencies','node_data');
        clear nodes elements idx node_data
    end
    if sum(evaluationGridNodes(:,1)>=400000 & evaluationGridNodes(:,1)<=499999)~=0
        idx=evaluationGridElements(:,1)>=400000 & evaluationGridElements(:,1)<=499999;
        elements=evaluationGridElements(idx,:)-min(min(evaluationGridElements(idx,:)));
        idx=evaluationGridNodes(:,1)>=400000 & evaluationGridNodes(:,1)<=499999;
        nodes=evaluationGridNodes(idx,:);
        nodes(:,1)=nodes(:,1)-min(nodes(:,1));
        node_data=pressure(:,idx,:);
        save('EvaluationGrid_LOW.mat','nodes','elements','frequencies','node_data');
        clear nodes elements idx node_data
    end
    if sum(evaluationGridNodes(:,1)>=500000 & evaluationGridNodes(:,1)<=599999)~=0
        idx=evaluationGridElements(:,1)>=500000 & evaluationGridElements(:,1)<=599999;
        elements=evaluationGridElements(idx,:)-min(min(evaluationGridElements(idx,:)));
        idx=evaluationGridNodes(:,1)>=500000 & evaluationGridNodes(:,1)<=599999;
        nodes=evaluationGridNodes(idx,:);
        nodes(:,1)=nodes(:,1)-min(nodes(:,1));
        node_data=pressure(:,idx,:);
        save('EvaluationGrid_HIGH.mat','nodes','elements','frequencies','node_data');
        clear nodes elements idx node_data
    end
end
    
%% Save HRTFs as SOFA file
if reciprocity

		SOFAstart;
		
    Obj = SOFAgetConventions('GeneralTF');

    Obj.GLOBAL_ApplicationName = 'Mesh2HRTF';
    Obj.GLOBAL_ApplicationVersion = '0.1.1';
    Obj.GLOBAL_Organization = '';
    Obj.GLOBAL_Title = '';
    Obj.GLOBAL_DateCreated = date;
    Obj.GLOBAL_AuthorContact = '';

    Obj.ReceiverPosition=receiverPositions;
    Obj.N=frequencies;

    %add division G([0,0,0],evaluationGrid)
    Obj.Data.Real=real(shiftdim(pressure,1));
    Obj.Data.Imag=imag(shiftdim(pressure,1));
    Obj.Data.Real_LongName='pressure';
    Obj.Data.Real_Units='pascal';
    Obj.Data.Imag_LongName='pressure';
    Obj.Data.Imag_Units='pascal';

    Obj.ListenerPosition = [0 0 0];
    Obj.SourcePosition_Type='cartesian';
    Obj.SourcePosition_Units='meter';
    Obj.SourcePosition=evaluationGridNodes(:,2:4);

    Obj=SOFAupdateDimensions(Obj);

    SOFAsave('EvaluationGrid.sofa',Obj);
end

end %of function
