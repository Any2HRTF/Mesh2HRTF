% D O C U M E N T A T I O N
% function WriteMesh2HRTFevalGrid(nodes, path, doPlot)
%
% writes evaluation grid in the format required by mesh2HRTF in the
% location specified by path. nodes specify the spatial sampling points by
% x, y, and z coordinates in meter. nodes of dimension Qx3, where Q are the
% number of sampling points.
% doPlot (true, false) specifies weather or not to plot the grid passed to
% the function.
%
% Author: fabian.brinkmann@tu-berlin.de
% Audio Communication Group @ TU Berlin
	
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
function WriteMesh2HRTFevalGrid(nodes, path, doPlot)

% triangulate the nodes
elements  = delaunayTriangulation(nodes(:,1), nodes(:,2), nodes(:,3)); 
% get only the surface of the triangulation
elements  = freeBoundary(elements);

% start count for numbering the nodes, and elements in Nodes.txt, and
% Elements.txt
start_count = 200000;

% check and make the path
if ~exist(path, 'dir')
    mkdir(path);
end

% write the nodes
fileID = fopen(fullfile(path, 'Nodes.txt'), 'w+');
fprintf(fileID, '%u\n', size(nodes,1));

for n = 1:size(nodes,1)
    fprintf(fileID, '%u %.6f %.6f %.6f\n', start_count+n-1, nodes(n,1), nodes(n,2), nodes(n,3));
end

fclose(fileID);

% write the elements
fileID = fopen(fullfile(path, 'Elements.txt'), 'w+');
fprintf(fileID, '%u\n', size(elements,1));

for n = 1:size(elements,1)
    fprintf(fileID, '%u %u %u %u %u %u %u\n', ...
            start_count+n-1, ...
            elements(n,1)+start_count-1, ...
            elements(n,2)+start_count-1, ...
            elements(n,3)+start_count-1, ...
            2, 0, 1);
end

fclose(fileID);

if doPlot
    
    figure
    set(gcf,'Color',[1 1 1]);
    trisurf(elements,nodes(:,1),nodes(:,2),nodes(:,3), 'FaceColor', [.8 .8 .8], 'faceAlpha', 1, 'EdgeColor', 'none')
    hold on
    scatter3(nodes(:,1), nodes(:,2), nodes(:,3), 75, '.k')
    axis equal
    xlabel x; ylabel y; zlabel z;
    
end