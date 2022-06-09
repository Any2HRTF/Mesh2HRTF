function writeBoundaryCondition(filename, kind, frequencies, data, comment)
% WRITEBOUNDARYCONDITION(filename, kind, frequencies, data, comment)
% 
% Write boundary condition to file.
% Mesh2HRTF supports non-rigid boundary conditions in the form of text
% files. Such files can be written wir this function.
% 
% INPUT:
% 
% filename : str
%         Name of the material file that is written to disk. Must end with ".csv"
%     kind : str
%         Defines the kind of boundary condition
%         ``"pressure"``
%             A pressure boundary condition can be used to force a certain
%             pressure on the boundary of the mesh. E.g., a pressure of 0 would
%             define a sound soft boundary.
%         ``"velocity"``
%             A velocity boundary condition can be used to force a certain
%             velocity on the boundary of the mesh. E.g., a velocity of 0 would
%             define a sound hard boundary.
%         ``admittance``
%             A normalized admittance boundary condition can be used to define
%             arbitrary boundaries. The admittance must be normalized, i.e.,
%             admittance/(rho*c) has to be provided, which rho being the density
%             of air in kg/m**3 and c the speed of sound in m/s.
%     frequencies : array like
%         The frequencies for which the boundary condition is given
%     data : array like
%         The values of the boundary condition at the frequencies given above.
%     comment : str, optional
%         A comment that is written to the beginning of the material file. The
%         default ``None`` does omit the comment.
%     Notes
%     -----
%     Mesh2HRTF performs an interpolation in case the boundary condition is
%     required at frequencies that are not specified. The interpolation is linear
%     between the lowest and highest provided frequency and uses the nearest
%     neighbor outside this range.




end
% EOF