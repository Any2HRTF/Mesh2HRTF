% demo
%  This script demonstrates the handling of evaluation grid in Mesh2HRTF. 
%  -----------------------------------------------------------------------------
%   NOTE: The evaluation grid has to be in meters regardless of the parameter
%         'unit' in Blender.
%  -----------------------------------------------------------------------------

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

% #Author: Katharina Pollack (ARI, ÖAW): 2022, original implementation
% #Author: Fabian Brinkmann (TU-Berlin): 2022, various improvements
% #Author: Piotr Majdak (ARI, ÖAW): 2023, help text, license boiler plate



%% --------------------------------------------- read a grid from Mesh2HRTF
close all; clear

% path of the grid
evaluationGrid = fullfile('Data', 'ARI');
% show a figure
doPlot         = true;

% get the evaluation grid
[nodes, numNodes, elems, numElems] = read_evaluation_grid(evaluationGrid, doPlot);


%% --------------------------------- write an equal angle grid to Mesh2HRTF
close all; clear

% path of the grid for saving
evaluationGrid = fullfile('Custom');

% equal angle grid in spherical coordinates
sg        = [repmat( (0:10:350)', [17 1] ) repelem( (80:-10:-80)', 36)];
% transform to Carthesian coordinates
[x, y, z] = sph2cart( sg(:,1)/180*pi, sg(:,2)/180*pi, ones(size(sg,1),1) );


write_evaluation_grid([x y z], evaluationGrid, true)

% save additional data - if wanted
% print('-dpdf', fullfile(evaluationGrid, 'EvaluationGrid'))
% save(fullfile(evaluationGrid, 'EvaluationGrid'), 'x', 'y', 'z')
