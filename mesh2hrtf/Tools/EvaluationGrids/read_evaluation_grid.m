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
%
function [nodes, numNodes, elems, numElems] = read_evaluation_grid(path, doPlot)
% function [nodes, numNodes, elems, numElems] = read_evaluation_grid(path, doPlot)
%
% reads the mesh2HRTF evaluation grid from the location specified by path.
% doPlot (true, false) specifies weather or not to plot the grid.
%
% OUTPUT:
% nodes    - nodes specified by numberm and x, y, z coordinates
% numNodes - number of nodes given by size(nodes, 1)
% elems    - triangular element specified by number and indices of nodes
% numElems - number of elements given by size(elems, 1)
%
% Author: fabian.brinkmann@tu-berlin.de
% Audio Communication Group @ TU Berlin

% get nodes and number of nodes
nodesID = fopen(fullfile(path, 'Nodes.txt'), 'r');

newLine = fgetl(nodesID);
numNodes = str2double(newLine);

nodes = zeros(10^6, 4);
n     = 0;

newLine = fgetl(nodesID);
while ischar(newLine)

    n = n+1;
    id = find(newLine == ' ');

    nodes(n,1) = str2double(newLine(1:id(1)-1));
    nodes(n,2) = str2double(newLine(id(1)+1:id(2)-1));
    nodes(n,3) = str2double(newLine(id(2)+1:id(3)-1));
    nodes(n,4) = str2double(newLine(id(3)+1:end));

    newLine = fgetl(nodesID);
end

fclose(nodesID);
nodes = nodes(1:n,:);

% get elements and number of elements
elemsID = fopen(fullfile(path, 'Elements.txt'), 'r');

newLine = fgetl(elemsID);
numElems = str2double(newLine);

elems = zeros(10^6, 4);
n     = 0;

newLine = fgetl(elemsID);
while ischar(newLine)

    n = n+1;
    id = find(newLine == ' ');

    elems(n,1) = str2double(newLine(1:id(1)-1));
    elems(n,2) = str2double(newLine(id(1)+1:id(2)-1));
    elems(n,3) = str2double(newLine(id(2)+1:id(3)-1));
    elems(n,4) = str2double(newLine(id(3)+1:id(4)-1));

    newLine = fgetl(elemsID);
end

fclose(elemsID);
elems = elems(1:n,:);

elems = elems - nodes(1,1)+1;
nodes(:,1) = nodes(:,1) - nodes(1,1)+1;

% plot if desired
if doPlot

    set(gcf,'Color',[1 1 1]);
    trisurf(elems(:,2:4),nodes(:,2),nodes(:,3),nodes(:,4), 'FaceColor', [.8 .8 .8], 'faceAlpha', 1, 'EdgeColor', 'none')
    hold on
    scatter3(nodes(:,2), nodes(:,3), nodes(:,4), 75, '.k')
    axis equal
    xlabel x; ylabel y; zlabel z;

end
