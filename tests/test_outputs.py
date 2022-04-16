import pytest
import numpy as np
import numpy.testing as npt
import subprocess
from tempfile import TemporaryDirectory
import shutil
import os
import glob
import pyfar as pf
import sofar as sf
import mesh2hrtf as m2h

cwd = os.path.dirname(__file__)
data_shtf = os.path.join(cwd, 'resources', 'SHTF')

def test_outputs_to_hrtfs():
    pass
