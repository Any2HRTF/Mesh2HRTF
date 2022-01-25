% Collect the data simulated by NumCalc and save to project folder.
close all; clear

Mesh2HRTF_version = '1.0.0';

% source information
sourceType = 'Right ear';
numSources = 1;
sourceCenter(1,1:3) = [0.087107 -0.003103 0.003888];
sourceArea(1,1) = 3.68754e-05;

% Reference to a point source in the origin
% accoring to the classical HRTF definition
% (https://doi.org/10.1016/0003-682X(92)90046-U)
reference = false;

% Compute HRIRs via the inverse Fourier transfrom.
% This will add data at 0 Hz, mirror the single sided spectrum, and
% shift the HRIRs in time. Requires reference = true.
computeHRIRs = false;

% Constants
speedOfSound = 343; % [m/s]
densityOfAir = 1.1839; % [kg/m^3]

% Collect the data simulated by NumCalc
Output2HRTF_Main(Mesh2HRTF_version, ...
                 sourceType, numSources, sourceCenter, sourceArea, ...
                 reference, computeHRIRs, ...
                 speedOfSound, densityOfAir);
