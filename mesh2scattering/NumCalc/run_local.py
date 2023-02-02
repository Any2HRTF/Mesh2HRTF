# %%
import mesh2scattering as m2s
import os

# this is the project path which we want to simualte
# it should contain a reference and a sample folder
project_path = '/Users/anne/sciebo/2021_DFG-Projekt/data/mesh2scattering/sine_10k'

# this is the path to the NumCalc executable, on Windows it will be
# a path like C:⁄NumCalc.exe
numcalc_path = 'NumCalc'

# %%

m2s.manage_numcalc(
    os.path.join(project_path, 'reference'),
    numcalc_path)

m2s.manage_numcalc(
    os.path.join(project_path, 'sample'),
    numcalc_path)


# %%
