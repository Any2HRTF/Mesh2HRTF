function export_vtk(path)
%   []=export_vtk(path)
%   exports data to VTK file for visualization in Paraview
%
%   Input:
%       path... The Mesh2HRTF project folder [string, the default is
%               the current working directory]

%% ----------------------check and initialize variables--------------------
if ~exist('path','var')
    path = pwd;
end

% create output directory
if ~exist(fullfile(path, 'Output2HRTF', 'ObjectMesh_vtk'),'dir')
    mkdir(fullfile(path, 'Output2HRTF', 'ObjectMesh_vtk'))
end

% load data struct created by output2hrtf.m
load(fullfile(path, 'Output2HRTF', 'ObjectMesh_Reference.mat'), 'elements', 'nodes', 'element_data')
nodes = nodes(:,2:end);
elements = elements(:,2:end);
data = 20*log10(abs(element_data{1})/0.00002);
datatype = 'amp';

% reset path to output directoy
path = fullfile(path, 'Output2HRTF', 'ObjectMesh_vtk');

temp=zeros(size(data,2),size(data,1),size(data,3));
for ii=1:size(data,3)
    temp(:,:,ii)=data(:,:,ii)';
end
data=temp;
clear temp


%% ----------------------------export data---------------------------------
for ii=1:size(data,2)
    file=fopen([path, filesep, datatype '_' num2str(ii) '.vtk'],'w');

    fprintf(file,'# vtk DataFile Version 3.0\n');
    fprintf(file,'Mesh2HRTF Files\n');
    fprintf(file,'ASCII\n');
    fprintf(file,'DATASET POLYDATA\n');

    fprintf(file,['POINTS ' num2str(size(nodes,1)) ' float\n']);
    fprintf(file,'%4.4f %4.4f %4.4f\n',transpose(nodes));

    fprintf(file,['POLYGONS ' num2str(size(elements,1)) ' ' num2str(size(elements,1)*(size(elements,2)+1)) '\n']);
    fprintf(file,'%i %i %i %i\n',transpose([size(elements,2)*ones(size(elements,1),1) elements]));

    fprintf(file,'\n');

	if size(data,1)==size(nodes,1)
    	fprintf(file,['POINT_DATA ' num2str(size(data,1)) '\n']);
	end
	if size(data,1)==size(elements,1)
		fprintf(file,['CELL_DATA ' num2str(size(data,1)) '\n']);
	end

    if size(data,3)>1
        fprintf(file,['SCALARS ' datatype '_l float 1\n']);
        fprintf(file,'LOOKUP_TABLE default\n');
        fprintf(file,'%5.5f\n',transpose(data(:,ii,1)));
        fprintf(file,['SCALARS ' datatype '_r float 1\n']);
        fprintf(file,'LOOKUP_TABLE default\n');
        fprintf(file,'%5.5f\n',transpose(data(:,ii,2)));
    else
        fprintf(file,['SCALARS ' datatype ' float 1\n']);
        fprintf(file,'LOOKUP_TABLE default\n');
        fprintf(file,'%5.5f\n',transpose(data(:,ii)));
    end

    fclose(file);
end

end %of function