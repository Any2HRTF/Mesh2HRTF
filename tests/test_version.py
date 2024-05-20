import mesh2hrtf as m2h
from packaging import version
import re
import os


def test_version():
    """Test if the version string is identical everywhere"""

    # version from top level version file
    with open('VERSION') as f:
        version_1 = f.readline().strip()

    # version from mesh2hrtf Python API
    version_2 = m2h.__version__

    # version from Blender export plug-in
    with open(os.path.join('mesh2hrtf', 'Mesh2Input', 'mesh2input.py')) as f:
        mesh2input = f.readlines()
    for line in mesh2input:
        if '"version":' in line:
            break
    # extract version from Blender export plug-in
    version_3 = re.search(r'\((.*?)\)', line).group(1).replace(', ', '.')

    # compare versions
    assert version.parse(version_1) == version.parse(version_2)
    assert version.parse(version_1) == version.parse(version_3)
