%                                Mesh2HRTF
%                Copyright (C) 2015 by Harald Ziegelwanger,
%        Acoustics Research Institute, Austrian Academy of Sciences
%                        mesh2hrtf.sourceforge.net
%
% If you use Mesh2HRTF:
% - Provide credits:
%   "Mesh2HRTF, H. Ziegelwanger, ARI, OEAW (mesh2hrtf.sourceforge.net)"
% - In your publication, cite both articles:
%   [1] Ziegelwanger, H., Kreuzer, W., and Majdak, P. (2015). "Mesh2HRTF: Open-source software package for the numerical calculation of head-related transfer functions," in Proceedings of the 22nd ICSV, Florence, IT.
%   [2] Ziegelwanger, H., Majdak, P., and Kreuzer, W. (2015). "Numerical calculation of listener-specific head-related transfer functions and sound localization: Microphone model and mesh discretization," The Journal of the Acoustical Society of America, 138, 208-222.
function write_evaluation_grid(nodes, path, doPlot, startCount)
% write_evaluation_grid(nodes, path, doPlot, startCount)
%
% writes evaluation grid in the format required by mesh2HRTF in the
% location specified by path.
%
% nodes      - a [Q x 3] matrix that specifies the x/y/z coordinates of the
%              evaluation grid in meter. Q is the number of sampling points
% path       - string that specifies the path where the data are saved.
% doPlot     - boolean that specifies weather or not to plot the sampling
%              grid (default = true)
% startCount - Mesh2HRTF saves the sampling nodes in Nodes.txt. This file
%              assigns an integer ID to each node. startCount specifies the
%              ID of the first node. The second note will hav the ID
%              startCount+1, etc (default = 200000).
%
%              IMPORTANT: If you use more than one evaluation grid in a
%              single simulation, each node in each grid must have a unique
%              ID. For example use a startCount of 200000 for the first
%              grid and a startCount of 300000 for the second.
%
% Author: fabian.brinkmann@tu-berlin.de
% Audio Communication Group @ TU Berlin

if ~exist('doPlot', 'var')
    doPlot = true;
end

% start count for numbering the nodes, and elements in Nodes.txt, and
% Elements.txt
if ~exist('startCount', 'var')
    startCount = 200000;
end

% prepare nodes for triangulation in case they are on a plane
constCoords = [all(nodes(1,1) == nodes(:,1)) ...
               all(nodes(1,2) == nodes(:,2)) ...
               all(nodes(1,3) == nodes(:,3))];

% triangulate the nodes
elements = delaunayTriangulation(nodes(:,~constCoords));

if any(constCoords)
    elements = elements.ConnectivityList;
else
    % get only the surface of the triangulation
    elements = freeBoundary(elements);
end

% check and make the path
if ~exist(path, 'dir')
    mkdir(path);
end

% write the nodes
fileID = fopen(fullfile(path, 'Nodes.txt'), 'w+');
fprintf(fileID, '%u\n', size(nodes,1));

for n = 1:size(nodes,1)
    fprintf(fileID, '%u %.6f %.6f %.6f\n', startCount+n-1, nodes(n,1), nodes(n,2), nodes(n,3));
end

fclose(fileID);

% write the elements
fileID = fopen(fullfile(path, 'Elements.txt'), 'w+');
fprintf(fileID, '%u\n', size(elements,1));

for n = 1:size(elements,1)
    fprintf(fileID, '%u %u %u %u %u %u %u\n', ...
            startCount+n-1, ...
            elements(n,1)+startCount-1, ...
            elements(n,2)+startCount-1, ...
            elements(n,3)+startCount-1, ...
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
