function write_material(filename, kind, frequencies, data, comment)
% write_material(filename, kind, frequencies, data, comment)
%
% Mesh2HRTF supports non-rigid boundary conditions. The function write_material 
% writes the boundary conditions to a file used later by NumCalc.
%
% Input parameters:
%
%   filename:   Name of the material file to be written to disk. Must end with ".csv"
%   kind:       Defines the type of the boundary condition: 
%               `pressure`: The pressure condition can be used to force to a certain
%                  pressure on the boundary of the mesh, e.g., a pressure of 0 defines
%                  a sound soft boundary.
%               `velocity`: The velocity condition can be used to force to a certain
%                 velocity on the boundary of the mesh, e.g., a velocity of 0 defines
%                 a sound hard boundary.
%               `admittance`: A normalized admittance condition can be used to define
%                 arbitrary boundaries. The admittance must be normalized, i.e.,
%                 admittance/(rho*c), with rho being the air density (kg/m^3) and 
%                 c being the speed of sound (m/s).
%   frequencies: The frequencies for which the boundary conditions are given.
%   data:        Vector with the boundary condition for each frequency.
%   comment:     Optional comment to be written to the beginning of the material file.
%
% Note that Mesh2HRTF performs an interpolation in case the boundary condition is
% required at frequencies that are not specified. The interpolation is linear
% between the lowest and highest provided frequency and uses the nearest
% neighbor outside this range.

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
% #Author: Fabian Brinkmann (TU-Berlin): 2022, integration in Mesh2HRTF 1.x
% #Author: Piotr Majdak (ARI, ÖAW): 2023, help text, license boiler plate



% check input
if isempty(filename)
    error('Error: filename invalid.')
end

[~, fname, fext] = fileparts(filename);
if ~strcmp(fext, '.csv')
    warning(['Warning: filename %s does not end with .csv\n', ...
        'Changed file extension from %s to .csv'], fname, fext)
    fext = 'csv';
    filename = strcat(fname, '.csv');
end

if length(frequencies) ~= length(data)
    error('Frequencies and data must have the same length.');
end

fileID = fopen(filename, 'w+');
if fileID == -1
    error('Error: file with filename %s cannot be opened', filename);
end

% write comment
if nargin == 5
    fprintf(fileID, comment);
    fprintf(fileID, '\n#\n');
end

% write kind of boundary condition
fprintf(fileID, '# Keyword to define the boundary condition:\n');
fprintf(fileID, '# ADMI: Normalized admittance boundary condition\n');
fprintf(fileID, '# PRES: Pressure boundary condition\n');
fprintf(fileID, '# VELO: Velocity boundary condition\n');
fprintf(fileID, '# NOTE: Mesh2HRTF expects normalized admittances, i.e., ');
fprintf(fileID, 'admittance/(rho*c).\n');
fprintf(fileID, '#       rho is the density of air and c the speed of sound.');
fprintf(fileID, 'The normalization is\n');
fprintf(fileID, '#       beneficial because a single material file can be used ');
fprintf(fileID, 'for simulations\n');
fprintf(fileID, '#       with differing speed of sound and density of air.\n');

switch kind
    case 'admittance'
        fprintf(fileID, 'ADMI\n');
    case 'pressure'
        fprintf(fileID, 'PRES\n');
    case 'velocity'
        fprintf(fileID, 'VELO\n');
    otherwise
        error('Kind of boundary condition must be ''admittance'', ''pressure'', or ''velocity.''');
end

fprintf(fileID, '#\n');
fprintf(fileID, '# Frequency curve:\n');
fprintf(fileID, '# Mesh2HRTF performs an interpolation in case the boundary ');
fprintf(fileID, 'condition is required\n');
fprintf(fileID, '# at frequencies that are not specified. The interpolation is ');
fprintf(fileID, 'linear between the\n');
fprintf(fileID, '# lowest and highest provided frequency and uses the nearest ');
fprintf(fileID, 'neighbor outside\n');
fprintf(fileID, '# this range.\n');
fprintf(fileID, '#\n');
fprintf(fileID, '# Frequency in Hz, real value, imaginary value\n');

% write data
for ii = 1:length(frequencies)
    fprintf(fileID, '%f, %f, %f\n', frequencies(ii), real(data(ii)), imag(data(ii)));
end

fclose(fileID);

end
% EOF