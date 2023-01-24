import os
import numpy as np


def read_ram_estimates(folder: str):
    """
    Read estimated RAM consumption from Memory.txt.

    Note that the RAM consumption per frequency step can be estimated and
    written to `Memory.txt` by calling ``NumCalc -estimate_ram``. This must
    be done before calling this function.

    Parameters
    ----------
    folder : str
        full path to the source folder containing the `Memory.txt` file from
        which the estimates are read

    Returns
    -------
    estimates : numpy array
        An array of shape ``(N, 3)`` where ``N`` is the number of frequency
        steps. The first column contains the frequency step, the second the
        frequency in Hz, and the third the estimated RAM consumption in GB.
    """

    # check if file exists
    if not os.path.isfile(os.path.join(folder, "Memory.txt")):
        raise ValueError(f"{folder} does not contain a Memory.txt file")

    # read content of file
    with open(os.path.join(folder, "Memory.txt"), "r") as ff:
        content = ff.readlines()

    # parse data to nested list
    estimates = []
    for line in content:
        estimate = []
        for ee in line.strip().split(" "):
            estimate.append(float(ee))

        estimates.append(estimate)

    return np.asarray(estimates)
