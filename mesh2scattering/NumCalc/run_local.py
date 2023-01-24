# %%
import mesh2hrtf as m2h
from utils import folders
import os
project_name = '01_kunsthaus_zuerich'
project_name = 'sine'
data_path = folders.sciebo_data_path()

# %%

m2h.manage_numcalc(
    os.path.join(
        data_path, 'mesh2hrtf_results', project_name, 'reference'),
    max_instances=8)

m2h.manage_numcalc(
    os.path.join(
        data_path, 'mesh2hrtf_results', project_name, 'sample'),
    starting_order='high',
    max_instances=8)


# %%
