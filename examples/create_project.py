# %%
import mesh2scattering as m2s
import pyfar as pf
import numpy as np
import os

# %%
# this is the project path which we want to simulate
# it should contain a reference and a sample folder
project_path = os.path.join(
    m2s.utils.repository_root(), '..', 'examples', 'project')
frequencies = np.array([1250, 2500, 5000])
path = os.path.join(
    m2s.utils.repository_root(), '..',
    'tests', 'resources', 'mesh', 'sine_5k')
sample_path = os.path.join(path, 'sample.stl')
reference_path = os.path.join(path, 'reference.stl')
receiver_delta_deg = 5
receiver_radius = 5
source_azimuth_deg = np.arange(0, 95, 30)
source_colatitude_deg = np.arange(10, 85, 30)
source_radius = 10

structural_wavelength = 0
sample_diameter = 0.8
modelScale = 2.5
symmetry_azimuth = [90, 180]
symmetry_rotational = False

if os.name == "nt":
    numcalc_path = os.path.join(
        m2s.utils.repository_root(), 'NumCalc', 'bin')
else:
    numcalc_path = os.path.join(
        m2s.utils.repository_root(), 'NumCalc', 'bin', 'NumCalc')

# %%
# create project
receiverCoords = pf.samplings.sph_equal_angle(
    receiver_delta_deg, receiver_radius)
receiverCoords = receiverCoords[receiverCoords.get_sph()[..., 1] < np.pi/2]
sourceCoords = m2s.input.create_source_positions(
    source_azimuth_deg, source_colatitude_deg, source_radius)

m2s.input.write_scattering_project(
    project_path=project_path,
    frequencies=frequencies,
    sample_path=sample_path,
    reference_path=reference_path,
    receiver_coords=receiverCoords,
    source_coords=sourceCoords,
    structuralWavelength=structural_wavelength,
    modelScale=modelScale,
    sample_diameter=sample_diameter,
    symmetry_azimuth=symmetry_azimuth,
    symmetry_rotational=symmetry_rotational,
    )

# %%
# run simulation
m2s.numcalc.manage_numcalc(
    os.path.join(project_path, 'reference'),
    numcalc_path)

m2s.numcalc.manage_numcalc(
    os.path.join(project_path, 'sample'),
    numcalc_path)

# %%
