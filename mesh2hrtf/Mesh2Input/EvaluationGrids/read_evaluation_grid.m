function [nodes, numNodes, elems, numElems] = read_evaluation_grid(path, doPlot)
% [nodes, numNodes, elems, numElems] = read_evaluation_grid(path, doPlot)
%
% read_evaluation_grid reads the evaluation grid from the file Nodes.txt
% located in a specific path and outputs its nodes and elements.
%
% Input parameters: 
%
%   path:   Path of the file Nodes.txt
%
%   doPlot: If true, the grid will be plotted for inspection. 
%
% Output parameters:
%
%   nodes:    Nodes of the grid specified by a node index and x, y, z coordinates
%
%   numNodes: total node number
%
%   elems:    triangular element specified by a number and 3 node indicies
%
%   numElems: total element number

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
