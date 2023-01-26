# %%
import os
import subprocess
import mesh2scattering as m2s

# %% here you can make changes
blender_executable = 'blender'

# %%
script = 'mesh2scattering/Mesh2Input/Tutorials/mesh2scattering.py'
path = m2s.repository_root()
if os.name == 'nt':  # Windows detected
    # run NumCalc and route all printouts to a log file
    subprocess.run(
        f"{blender_executable} --background --python {script}",
        stdout=subprocess.DEVNULL, cwd=path, check=True)

else:  # elif os.name == 'posix': Linux or Mac detected
    # run NumCalc and route all printouts to a log file
    subprocess.run(
        [f"{blender_executable} --background --python {script}"],
        shell=True, stdout=subprocess.DEVNULL, cwd=path, check=True)