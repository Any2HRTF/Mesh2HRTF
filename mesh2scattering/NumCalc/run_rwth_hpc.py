# %%

import os
import subprocess
import numpy as np
import mesh2scattering as m2s


# %%
# create
project_name_out = '01_kunsthaus_zuerich'
numcalc_executable = 'NumCalc'
program_path = '/home/anne/git/Mesh2scattering/mesh2scattering'
data_path = '/home/anne/sciebo/2021_DFG-Projekt/data'

# %%
file_path = os.path.join(data_path, 'mesh2hrtf', project_name_out)

sample_source = os.path.join(file_path, 'sample', 'NumCalc', 'source_2')
ref_source = os.path.join(file_path, 'reference', 'NumCalc', 'source_2')
paths = [sample_source, ref_source]
for path in paths:
    if not os.path.isfile(os.path.join(path, "Memory.txt")):
        if os.name == 'nt':  # Windows detected
            # run NumCalc and route all printouts to a log file
            subprocess.run(
                f"{numcalc_executable} -estimate_ram",
                stdout=subprocess.DEVNULL, cwd=path, check=True)

        else:  # elif os.name == 'posix': Linux or Mac detected
            # run NumCalc and route all printouts to a log file
            subprocess.run(
                [f"{numcalc_executable} -estimate_ram"],
                shell=True, stdout=subprocess.DEVNULL, cwd=path, check=True)


# %%
ram = []
for idx in range(len(paths)):
    data = m2s.read_ram_estimates(paths[idx])
    data = np.append(data, np.zeros((data.shape[0], 1))+idx, axis=1)
    ram.append(data)

ram = np.vstack(ram)
cores = (np.array(ram[:, 2]*1.1/4, dtype=int)+1)
ram = np.append(ram, cores.reshape((len(cores), 1)), axis=1)
ram

# %%
cores_str = '$$CORES$$'
times_str = '$$TIME$$'
name_str = '$$NAME$$'
array_str = '$$ARRAY$$'
path_str = '$$PATH$$'
folder_str = '$$TYPE$$'
index_str = '$$INDEX$$'

times = '00-03:00:00'
path = f'$HOME/Dokumente/comsol_hpc/{project_name_out}'

# read draft
cwd = os.path.join(program_path, 'NumCalc')
with open(os.path.join(cwd, 'rwth_hpc_draft.sh')) as f:
    lines = f.read()
shell_scripte = []
for idx in range(ram.shape[0]):
    cores = int(ram[idx, 4])
    if ram[idx, 3] == 0:  # sample
        array = '1-80'
        folder = 'sample'
    else:
        array = '1-8'
        folder = 'reference'
    index = int(ram[idx, 0])

    all_files, fundamentals, out, out_names = m2s.check_project(
        os.path.join(file_path, folder))

    array_list = []
    is_error = False
    for ss in range(out.shape[2]):
        f = out[index-1, :, ss]
        if any(f < 0):
            array_list.append(f'{ss+1}')
            is_error = True
    if not is_error:
        continue
    array = ','.join(array_list)

    name = f'{project_name_out}_{folder}_{index}'
    # fill in form
    shell = lines.replace(cores_str, f'{cores}')
    shell = shell.replace(times_str, times)
    shell = shell.replace(name_str, name)
    shell = shell.replace(array_str, array)
    shell = shell.replace(path_str, path)
    shell = shell.replace(folder_str, folder)
    shell = shell.replace(index_str, f'{index}')

    file_out = os.path.join(file_path, f'{name}.sh')
    shell_scripte.append(f'{name}.sh')
    with open(file_out, "w") as f:
        f.write(shell)
        f.close()

for script in shell_scripte:
    print(f'sbatch < {script}')
print('')

# %%

# sftp ah664066@login18-x-1.hpc.itc.rwth-aachen.de

