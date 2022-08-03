% Author: fabian.brinkmann@tu-berlin.de
% Audio Communication Group @ TU Berlin

%% -----------------------------------------------------------------------------
%  NOTE: The evaluation grid has to be in meters regardless of the parameter
%        'unit' in Blender.
%  -----------------------------------------------------------------------------

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
