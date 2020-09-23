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