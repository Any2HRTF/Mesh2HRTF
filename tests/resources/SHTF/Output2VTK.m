% Export the sound pressure files genereated by Output2HRTF.m as SPL to VTK-Files for visualization in Paraview.

close all; clear

Mesh2HRTF_version = '1.0.0';

if ~exist(fullfile(pwd, 'Visualization'),'dir')
    mkdir(fullfile(pwd, 'Visualization'));
end
if ~exist(fullfile(pwd, 'Visualization', 'ObjectMesh'),'dir')
    mkdir(fullfile(pwd, 'Visualization', 'ObjectMesh'))
end

load(fullfile(pwd, 'Output2HRTF', 'ObjectMesh_Reference.mat'))

EvalToolsExport2VTK(['Visualization' filesep 'ObjectMesh' filesep],nodes(:,2:end),elements(:,2:end),20*log10(abs(element_data{1})/0.00002),'amp')
