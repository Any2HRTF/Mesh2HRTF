# %%
import mesh2scattering as m2s
import pyfar as pf
import os
path = r"D:\sciebo\2021_DFG-Projekt\data\mesh2hrtf\sine"

# %%
m2s.output2scattering(path, 0.177/2.5)

# %%
m2s.scattering.calc_coefficient(path, '5m_1deg')

# %%
#%matplotlib qt
s_rand_path = os.path.join(path, 'sample_coords.scattering_rand.sofa')

s_rand, _, _ = pf.io.read_sofa(s_rand_path)
pf.plot.freq(s_rand, dB=False)

# %%
s_path = os.path.join(path, 'sample_coords.scattering.sofa')

s, _, _ = pf.io.read_sofa(s_path)
pf.plot.freq(s, dB=False)

# %%
