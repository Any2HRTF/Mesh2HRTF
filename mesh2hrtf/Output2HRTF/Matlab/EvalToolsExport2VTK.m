function EvalToolsExport2VTK(path,nodes,elements,data,datatype)
%EvalTools_EXPORT2VTK
%   []=EvalTools_export2VTK(path,nodes,elements,data,datatype) exports
%   data to VTK file for visualization in Paraview
%
%   Input:
%       path:
% 
%       nodes:
% 
%       elements:
% 
%       data:

%% ----------------------check and initialize variables--------------------
if ~exist('path','var') || ~exist('nodes','var') || ~exist('elements','var') || ~exist('data','var') || ~exist('datatype','var')
    error('Not enough input arguments.')
end

temp=zeros(size(data,2),size(data,1),size(data,3));
for ii=1:size(data,3)
    temp(:,:,ii)=data(:,:,ii)';
end
data=temp;
clear temp


%% ----------------------------export data---------------------------------
for ii=1:size(data,2)
    file=fopen([path datatype '_' num2str(ii) '.vtk'],'w');

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