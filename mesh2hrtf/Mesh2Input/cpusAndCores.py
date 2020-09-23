#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  4 16:12:36 2020

@author: root
"""

# export variables - try different parameters and see if everything works as expected
cpuFirst          = 2
cpuLast           = 3
numCoresPerCPU    = 8
numEars           = 2
frequencyStepSize = 200
maxFrequency      = 20_000


# new variable:
# maximum number of allowed cpus - previously hard coded as number
# could be maxCPUs = cpuLast if the rests of Mesh2HRTF works with that
maxCPUs  = 10

# new variable:
# maximum nuber of allowed cores - previously hard coded as number
# could be numCoresPerCPU if the rest of Mesh2HRTF works with that
maxCores = 8
minFrequency = frequencyStepSize


# %% code at line 454
numCPUs = cpuLast-cpuFirst+1


# %% original script lines 464-475
lowFrequency = 0


# %% original block starting at 639

# number of cores per ear and in total
numCoresUsedPerEar    = numCPUs*numCoresPerCPU//numEars
if not numCoresUsedPerEar:
    raise Exception("At least two cores must be available for calculating both ears, i.e., two CPUs with one core each or one CPU with two cores.")

# Calculate number of frequencies to be calculated (MUST be integer)
frequencySteps = divmod(maxFrequency-lowFrequency, frequencyStepSize)
if not frequencySteps[1] == 0:
    raise Exception("Error, frequencyStepSize is not a divisor of maxFrequency-lowFrequency")

# all frequencies to be calculated
f = [ff*frequencyStepSize+minFrequency for ff in range(frequencySteps[0])]
if len(f) < numCoresUsedPerEar:
    raise Exception("More cores than frequencies, i.e., numCPUs*numCoresPerCPU//numEars < numFrequencies.")

# distribution across numCoresUsedPerEar
F = [[] for ff in range(numCoresUsedPerEar)]
for nn, ff in enumerate(f):
    F[nn%numCoresUsedPerEar].append(ff)

# Initialize cpusAndCores:
# Nested list that indicates which cpu and core is used to calculate data for
# which ear (e.g. cpuAndCores[0][1] holds the entry for the second core on the
# first cpu). Entries: 0=idle, 1=leftEar, 2=rightEar.
cpusAndCores = [[0]*maxCores for cc in range(maxCPUs)]

# Initialize frequencies:
# Nested list that holds the frequencies that are calculated by each cpu/core
# (e.g. frequencies[0][1] holds a list of frequencies calculated by the second
# core on the first cpu.
tmp         = [[] for cc in range(maxCores)]
frequencies = [tmp.copy() for cc in range(maxCPUs)]

# distribute ears and frequencies across cpus and cores.
# Left ear is calculated on cpus 0 to numCoresUsedPerEar-1
# Right ear is calculated on cpus numCoresUsedPerEar to numCoresAvailable
for count in range(numCoresUsedPerEar*numEars):
    cpu, core = divmod(count + (cpuFirst-1)*numCoresPerCPU, numCoresPerCPU)
    cpusAndCores[cpu][core] = count//numCoresUsedPerEar + 1
    frequencies[cpu][core]  = F[count%len(F)]
    # output for debugging
    print(f"CPU {cpu+1:2d}, core {core+1}, ear {count//numCoresUsedPerEar + 1}, freqList {count%len(F)}")


# %% deprecated variables

# only used for documentation replace with 'len(f)//numCoresUsedPerEar'
# frequencyStepsPerCore = divmod(frequencySteps[0], numCoresUsedPerEar)
# only used for documentation replace with 'numCoresUsedPerEar*numEars'
# numCoresAvailable = numCoresUsedPerEar*numEars

# remove all occurences
# lowFrequencyCores = 0