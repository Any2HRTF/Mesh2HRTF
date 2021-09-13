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
%
% Author: fabian.brinkmann@tu-berlin.de
% Audio Communication Group @ TU Berlin

%% -----------------------------------------------------------------------------
%  NOTE: The evaluation grid has to be in meters regardless of the parameter
%        'unit' in Blender.
%  -----------------------------------------------------------------------------

%% --------------------------------------------- read a grid from Mesh2HRTF
close all; clear

% path of the grid
evaluationGrid = fullfile('..', '..', 'Mesh2Input', 'EvaluationGrids', 'ARI');
% show a figure
doPlot         = true;

% get the evaluation grid
[nodes, numNodes, elems, numElems] = getMesh2HRTFevalGrid(evaluationGrid, doPlot);


%% --------------------------------- write an equal angle grid to Mesh2HRTF
close all; clear

% path of the grid for saving
evaluationGrid = fullfile('..', '..', 'Mesh2Input', 'EvaluationGrids', 'Custom');

% equal angle grid in spherical coordinates
sg        = [repmat( (0:10:350)', [17 1] ) repelem( (80:-10:-80)', 36)];
% transform to Carthesian coordinates
[x, y, z] = sph2cart( sg(:,1)/180*pi, sg(:,2)/180*pi, ones(size(sg,1),1) );


writeMesh2HRTFevalGrid([x y z], evaluationGrid, true)

% save additional data - if wanted
% print('-dpdf', fullfile(evaluationGrid, 'EvaluationGrid'))
% save(fullfile(evaluationGrid, 'EvaluationGrid'), 'x', 'y', 'z')
