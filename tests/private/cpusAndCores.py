#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  4 16:12:36 2020

@author: root
"""

# %% test


def _distribute_frequencies(cpuFirst, cpuLast, maxCPUs,
                            numCoresPerCPU, maxCores,
                            numEars,
                            minFrequency, maxFrequency,
                            frequencyStepSize, numFrequencySteps):
    """Calculate list of frequencies and distribute across cores and CPUs.

    Returns
    -------

    cpusAndCores: list
        Nested list that indicates which cpu and core is used to calculate data
        for which ear (e.g. cpuAndCores[0][1] holds the entry for the second
        core on the first cpu). Entries: 0=idle, 1=leftEar/right ear if
        calculating one ear, 2=rightEar if calculating two ears.
    frequencies: list
        Nested list that holds the frequencies that are calculated by each
        cpu/core (e.g. frequencies[0][1] holds a list of frequencies calculated
        by the second core on the first cpu.
    frequencyStepsPerCore: int
        The number of frequency steps calculated per core (written to
        Info.txt).
    numCoresAvailable: int
        The number of cores used for the computation (written to Info.txt).
    f: list
        A simple list of the frequencies to be simulated (for debugging).
    frequencyStepSize: float
        Step size between sucsessive frequncies (written to Info.txt). This is
        returned because it might be zero during the function call.
    numFrequencySteps: int
        Number of frequncies to be simulated (written to Info.txt). This is
        returned because it might be zero during the function call.

    """

    # number of CPUs used
    numCPUs = cpuLast-cpuFirst+1

    # number of cores per ear
    numCoresUsedPerEar = numCPUs*numCoresPerCPU//numEars
    if not numCoresUsedPerEar:
        raise Exception("At least two cores must be available for calculating "
                        "both ears, i.e., two CPUs with one core each or one "
                        "CPU with two cores.")

    # check input
    if (numFrequencySteps == 0 and frequencyStepSize == 0) \
            or (numFrequencySteps != 0 and frequencyStepSize != 0):
        raise Exception("Either 'frequencyStepSize' or 'numFrequencySteps' \
                        must be zero while the other must not.")

    # Calculate Number of frequencies and frequency step size
    if minFrequency == maxFrequency:
        frequencySteps = (1, 0)
        frequencyStepSize = 0
    elif frequencyStepSize:
        frequencySteps = divmod(
            maxFrequency - minFrequency + frequencyStepSize, frequencyStepSize)
    else:
        if numFrequencySteps < 2:
            raise Exception("'numFrequencySteps' must be at least 2.")
        frequencySteps = (numFrequencySteps, 0)
        frequencyStepSize = (maxFrequency-minFrequency)/(numFrequencySteps-1)

    if not frequencySteps[1] == 0:
        raise Exception("Error, frequencyStepSize is not a divisor of \
                        maxFrequency-minFrequency")

    # get all frequencies to be calculated
    f = [ff*frequencyStepSize+minFrequency for ff in range(frequencySteps[0])]

    # remove 0 Hz if included in the list
    if f[0] == 0:
        f.pop(0)
        frequencySteps = (frequencySteps[0]-1, 0)
        print('Warning: 0 Hz can not be calculated and was removed from the \
              list of frequencies.')

    if not len(f):
        raise ValueError("No frequncies to be calulated. \
                         Check the input parameters.")

    # check number of cores and frequencies
    if len(f) < numCoresUsedPerEar:
        raise Exception("More cores than frequencies, i.e., \
                        numCPUs*numCoresPerCPU//numEars < numFrequencies.")

    # distribution of frequencies across numCoresUsedPerEar
    F = [[] for ff in range(numCoresUsedPerEar)]
    for nn, ff in enumerate(f):
        F[nn % numCoresUsedPerEar].append(ff)

    # Initialize cpusAndCores:
    cpusAndCores = [[0]*maxCores for cc in range(maxCPUs)]

    # Initialize frequencies:
    freqs = [[] for cc in range(maxCores)]
    frequencies = [freqs.copy() for cc in range(maxCPUs)]

    # distribute ears and frequencies across cpus and cores.
    # Left ear is calculated on cpus 0 to numCoresUsedPerEar-1
    # Right ear is calculated on cpus numCoresUsedPerEar to numCoresAvailable
    for count in range(numCoresUsedPerEar*numEars):
        cpu, core = divmod(count + (cpuFirst-1)*numCoresPerCPU, numCoresPerCPU)
        cpusAndCores[cpu][core] = count//numCoresUsedPerEar + 1
        frequencies[cpu][core] = F[count % len(F)]
        # output for debugging
        # print(f"CPU {cpu+1:2d}, core {core+1}, \
        #       ear {count//numCoresUsedPerEar + 1}, \
        #       freqList {count%len(F)}")

    # meta data for Info.txt
    frequencyStepsPerCore = len(f)//numCoresUsedPerEar
    numCoresAvailable = numCoresUsedPerEar*numEars

    return cpusAndCores, frequencies, \
        frequencyStepsPerCore, numCoresAvailable, \
        f, frequencyStepSize, numFrequencySteps


# %%
# export variables - try different parameters and see if everything works as
# expected
cpuFirst = 1
cpuLast = 1
numCoresPerCPU = 2
numEars = 2

minFrequency = 0       # new Variable must be added to Export-GUI
maxFrequency = 22_050

# either 'frequencyStepSize' or 'frequencySteps' must be 0
frequencyStepSize = 0
# new variabel must be added to Export-GUI
numFrequencySteps = 129

# maximum number of cpus. Still hard coded but might be:
# `maxCPUs = max(10, cpuLast)` if it works with the rest of Mesh2HRTF
maxCPUs = 10

# maximum number of cores. Still hard coded but might be:
# maxCores = max(8, numCoresPerCPU) if it works with the rest of Mesh2HRTF
maxCores = 8

cpusAndCores, frequencies, frequencyStepsPerCore, numCoresAvailable, f, *_ = \
    _distribute_frequencies(cpuFirst, cpuLast, maxCPUs,
                            numCoresPerCPU, maxCores,
                            numEars,
                            minFrequency, maxFrequency,
                            frequencyStepSize, numFrequencySteps)

# %% tests
maxCPUs = 10
maxCores = 8

# %% test single frequency calculation
cpusAndCores, frequencies, frequencyStepsPerCore, numCoresAvailable, f, *_ = \
    _distribute_frequencies(
        1, 1, maxCPUs,  # cpuFirst, cpuLast, maxCPUs
        1, maxCores,    # numCoresPerCPU, maxCores
        1,              # numEars
        100, 100,       # minFrequency, maxFrequency
        0, 1)           # frequencyStepSize, numFrequencySteps


# %% test single sided spectrum by numFrequencySteps
cpusAndCores, frequencies, frequencyStepsPerCore, numCoresAvailable, f, *_ = \
    _distribute_frequencies(
        1, 1, maxCPUs,  # cpuFirst, cpuLast, maxCPUs
        1, maxCores,    # numCoresPerCPU, maxCores
        1,              # numEars
        0, 22050,       # minFrequency, maxFrequency
        0, 129)         # frequencyStepSize, numFrequencySteps

# %% test single sided spectrum by frequencyStepSize
cpusAndCores, frequencies, frequencyStepsPerCore, numCoresAvailable, f, *_ = \
    _distribute_frequencies(
        1, 1, maxCPUs,  # cpuFirst, cpuLast, maxCPUs
        1, maxCores,    # numCoresPerCPU, maxCores
        1,              # numEars
        0, 22000,       # minFrequency, maxFrequency
        100, 0)         # frequencyStepSize, numFrequencySteps


# %% random debugging
cpusAndCores, frequencies, frequencyStepsPerCore, numCoresAvailable, f, *_ = \
    _distribute_frequencies(
        1, 1, maxCPUs,  # cpuFirst, cpuLast, maxCPUs
        1, maxCores,    # numCoresPerCPU, maxCores
        1,              # numEars
        100, 100,       # minFrequency, maxFrequency
        0, 129)         # frequencyStepSize, numFrequencySteps

print(f[0], f[-1])
