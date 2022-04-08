% Collect the data simulated by NumCalc and save to project folder.
close all; clear

Mesh2HRTF_version = '1.0.0';

% source information
sourceType = 'pointSource';
sourceCenter(1,1:3) = [0.0 0.20000000298023224 0.0];
sourceArea(1,1)     = 1;
% Reference to a point source in the origin
% accoring to the classical HRTF definition
% (https://doi.org/10.1016/0003-682X(92)90046-U)
reference = false;

% Compute HRIRs via the inverse Fourier transfrom.
% This will add data at 0 Hz, mirror the single sided spectrum, and
% shift the HRIRs in time. Requires reference = true.
computeHRIRs = false;

% Constants
speedOfSound = 343.18; % [m/s]
densityOfAir = 1.1839; % [kg/m^3]

% Distribution of ears across CPUs and cores.
% (Matrix of size [numCPUs x numCores])
cpusAndCores = [
    1];

% Collect the data simulated by NumCalc
Output2HRTF_Main(Mesh2HRTF_version, cpusAndCores, ...
                 sourceType, sourceCenter, sourceArea, ...
                 reference, computeHRIRs, ...
                 speedOfSound,densityOfAir);
