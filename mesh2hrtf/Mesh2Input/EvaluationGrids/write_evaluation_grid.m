function write_evaluation_grid(nodes, path, doPlot, startCount)
% write_evaluation_grid(nodes, path, doPlot, startCount)
%
% write_evaluation_grid writes the evaluation grid to the file Nodes.txt
% which can be then used to calculate HRTFs for sources specified at that grid. 
% Each node is represented by an integer ID, with startCount specifying the
% ID of the first node.
%
% Input parameters:
%   nodes:      Q-by-3 matrix specifying the x, y, and z coordinates (in meter)
%               of the evaluation grid, with Q being the number of nodes
%               in the evaluation grid.
%                
%   path        Path to save the data to be saved.
%
%   doPlot      If true, the evaluation grid will be plotted for inspection. 
%               Default: true.
%
%   startCount: ID of the first node. Default: 200000.
%
% IMPORTANT: If multiple evaluation grids are used in a single simulation, 
% each node ID must be unique within the joint evaluation grid. 

% This file is part of the Mesh2HRTF software package developed by the
% Mesh2HRTF Developer Team (https://mesh2hrtf.org) and licensed under the 
% EUPL, Version 1.2, or, as soon as approved by the European Commission, 
% subsequent versions of the EUPL. Details on the license can be found 
% in the file "license.txt" provided with Mesh2HRTF package
% or at https://joinup.ec.europa.eu/collection/eupl/eupl-text-eupl-12
%
% You may not use this work except in compliance with the license.
% Unless required by applicable law or agreed to in writing, software 
% distributed under the license is distributed on an "AS IS" basis,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.

% #Author: Fabian Brinkmann (TU-Berlin): 2018, original implementation
% #Author: Fabian Brinkmann (TU-Berlin): 2022, integration in Mesh2HRTF 1.x
% #Author: Piotr Majdak (ARI, ÖAW): 2023, help text, license boiler plate


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
