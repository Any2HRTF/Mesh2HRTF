function write_material(filename, kind, frequencies, data, comment)
% write_material(filename, kind, frequencies, data, comment)
%
% Write boundary condition to file.
% Mesh2HRTF supports non-rigid boundary conditions in the form of text
% files. Such files can be written with this function.
%
% INPUT:
%
% filename ... Name of the material file that is written to disk. Must end with ".csv"
% kind ....... Defines the kind of boundary condition
%         pressure
%             A pressure boundary condition can be used to force a certain
%             pressure on the boundary of the mesh. E.g., a pressure of 0 would
%             define a sound soft boundary.
%         velocity
%             A velocity boundary condition can be used to force a certain
%             velocity on the boundary of the mesh. E.g., a velocity of 0 would
%             define a sound hard boundary.
%         admittance
%             A normalized admittance boundary condition can be used to define
%             arbitrary boundaries. The admittance must be normalized, i.e.,
%             admittance/(rho*c) has to be provided, which rho being the density
%             of air in kg/m**3 and c the speed of sound in m/s.
% frequencies ... The frequencies for which the boundary condition is given
% data .......... The values of the boundary condition at the frequencies given above.
% comment ....... A comment that is written to the beginning of the material file.
%                 The default ``None`` does omit the comment.
%
% Mesh2HRTF performs an interpolation in case the boundary condition is
% required at frequencies that are not specified. The interpolation is linear
% between the lowest and highest provided frequency and uses the nearest
% neighbor outside this range.

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