# %%
import mesh2scattering as m2s
import pyfar as pf
import numpy as np
import os

# %% 
# this is the project path which we want to simualte
# it should contain a reference and a sample folder
project_path = os.path.join(
    m2s.utils.repository_root(), '..', 'examples', 'project')
frequencies = pf.dsp.filter.fractional_octave_frequencies(3, (500, 5000))[0]
path = os.path.join(
    m2s.utils.repository_root(), '..', 
    'tests', 'resources', 'mesh', 'sine_5k')
sample_path = os.path.join(path, 'sample.stl')
reference_path = os.path.join(path, 'reference.stl')
receiver_delta_deg = 1
receiver_radius = 5
source_azimuth_deg = np.arange(0, 95, 10)
source_colatitude_deg = np.arange(10, 85, 10)
source_radius = 10

strcutal_wavelength = 0
sample_diameter = 0.8
modelScale = 2.5
symmetry_azimuth = [90, 180]
symmetry_rotational = False

numcalc_path = os.path.join(
    m2s.utils.repository_root(), 'NumCalc', 'bin', 'NumCalc')

# %% 
# create project
receiverPoints = pf.samplings.sph_equal_angle(
    receiver_delta_deg, receiver_radius)
receiverPoints = receiverPoints[receiverPoints.get_sph()[..., 1]<np.pi/2]
sourcePositions = m2s.input.create_source_positions(
    source_azimuth_deg, source_colatitude_deg, source_radius)

m2s.input.write_scattering_project(
    project_path=project_path,
    frequencies=frequencies,
    sample_path=sample_path,
    reference_path=reference_path,
    receiverPoints=receiverPoints, 
    sourcePositions=sourcePositions,
    structualWavelength=0,
    modelScale=1, 
    sample_diameter=sample_diameter,
    symmetry_azimuth=symmetry_azimuth,
    symmetry_rotational=symmetry_rotational,
    )

# %%
# rund simulation
m2s.NumCalc.manage_numcalc(
    os.path.join(project_path, 'reference'),
    numcalc_path)

m2s.NumCalc.manage_numcalc(
    os.path.join(project_path, 'sample'),
    numcalc_path)
# %%
