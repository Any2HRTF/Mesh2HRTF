# %%
import mesh2scattering as m2s
import pyfar as pf
import os

# this is the project path which we want to simualte
# it should contain a reference and a sample folder
project_path = '/Users/anne/sciebo/2021_DFG-Projekt/data/mesh2hrtf/ita'

# this is the structural wavelength of the test sample, it its not clear use 0
structural_wavelength = 1.

# this is the evaluation grid, means the microphone grid, which was defined
# in the mesh2scattering.py
evaluation_grid = '5m_1deg'

# %%
# export the scattering pattern for test sample and the reference sample to a
# sofa file since the reference sample is always rotational symmetric, the
# data for the missing angles are rotated in this was sample and reference
# data will have the same dimensions and coordinates
m2s.output2scattering(project_path, structural_wavelength)

# %%
# calcualtes the scattering coeffient for each incident angle and the random
# one from the scattering pattern
m2s.scattering.calc_coefficient(project_path, evaluation_grid)

# %%
# example of plotting the random scattering coefficient
project_name = os.path.split(project_path)[-1]
s_rand_path = os.path.join(
    project_path, f'{project_name}_{evaluation_grid}.scattering_rand.sofa')

s_rand, _, _ = pf.io.read_sofa(s_rand_path)
pf.plot.freq(s_rand, dB=False)

# %%
s_path = os.path.join(
    project_path, f'{project_name}_{evaluation_grid}.scattering.sofa')

s, source_pos, _ = pf.io.read_sofa(s_path)
pf.plot.freq(s, dB=False)

# %%
