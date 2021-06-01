# -*- coding: utf-8 -*-

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
#   createTestFile.py
#
#   Reference implementation for SOFA file creation with pysofaconventions
#
#   (C) Andrés Pérez-López - Eurecat / UPF
#   24/08/2018
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

from netCDF4 import Dataset
import time
import numpy as np
import os

#----------Create it----------#

filePath = "testpysofaconventions.sofa"
# Need to delete it first if file already exists
if os.path.exists(filePath):
    os.remove(filePath)
rootgrp = Dataset(filePath, 'w', format='NETCDF4')


#----------Required Attributes----------#

rootgrp.Conventions = 'SOFA'
rootgrp.Version = '1.0'
rootgrp.SOFAConventions = 'SimpleFreeFieldHRIR'
rootgrp.SOFAConventionsVersion = '0.1'
rootgrp.APIName = 'pysofaconventions'
rootgrp.APIVersion = '0.1'
rootgrp.APIVersion = '0.1'
rootgrp.AuthorContact = 'andres.perez@eurecat.org'
rootgrp.Organization = 'Eurecat - UPF'
rootgrp.License = 'WTFPL - Do What the Fuck You Want to Public License'
rootgrp.DataType = 'FIR'
rootgrp.RoomType = 'reverberant'
rootgrp.DateCreated = time.ctime(time.time())
rootgrp.DateModified = time.ctime(time.time())
rootgrp.Title = 'testpysofaconventions'
rootgrp.RoomType = 'free field'
rootgrp.DatabaseName = 'CoolDatabase'
rootgrp.ListenerShortName = '001'



#----------Required Dimensions----------#

m = 3
n = 48000
r = 2
e = 1
i = 1
c = 3
rootgrp.createDimension('M', m)
rootgrp.createDimension('N', n)
rootgrp.createDimension('E', e)
rootgrp.createDimension('R', r)
rootgrp.createDimension('I', i)
rootgrp.createDimension('C', c)


#----------Required Variables----------#
listenerPositionVar = rootgrp.createVariable('ListenerPosition',    'f8',   ('I','C'))
listenerPositionVar.Units   = 'metre'
listenerPositionVar.Type    = 'cartesian'
listenerPositionVar[:] = np.zeros(c)

listenerUpVar       = rootgrp.createVariable('ListenerUp',          'f8',   ('I','C'))
listenerUpVar.Units         = 'metre'
listenerUpVar.Type          = 'cartesian'
listenerUpVar[:]    = np.asarray([0,0,1])

# Listener looking to the left (+Y axis)
listenerViewVar     = rootgrp.createVariable('ListenerView',        'f8',   ('I','C'))
listenerViewVar.Units       = 'metre'
listenerViewVar.Type        = 'cartesian'
listenerViewVar[:]  = np.asarray([0,1,0])

emitterPositionVar  = rootgrp.createVariable('EmitterPosition',     'f8',   ('E','C','I'))
emitterPositionVar.Units   = 'metre'
emitterPositionVar.Type    = 'spherical'
# Equidistributed speakers in circle
emitterPositionVar[:] = np.zeros((e,c,i))
# for idx in range(e):
#     azi = idx*360/float(e)
#     ele = 0.
#     dis = 1.
#     emitterPositionVar[idx,:,:] = np.asarray([azi,ele,dis])

sourcePositionVar = rootgrp.createVariable('SourcePosition',        'f8',   ('I','C'))
sourcePositionVar.Units   = 'metre'
sourcePositionVar.Type    = 'cartesian'
sourcePositionVar[:]      = np.zeros(c)

sourceUpVar       = rootgrp.createVariable('SourceUp',              'f8',   ('I','C'))
sourceUpVar.Units         = 'metre'
sourceUpVar.Type          = 'cartesian'
sourceUpVar[:]    = np.asarray([0,0,1])

sourceViewVar     = rootgrp.createVariable('SourceView',            'f8',   ('I','C'))
sourceViewVar.Units       = 'metre'
sourceViewVar.Type        = 'cartesian'
sourceViewVar[:]  = np.asarray([1,0,0])

receiverPositionVar = rootgrp.createVariable('ReceiverPosition',  'f8',   ('R','C','I'))
receiverPositionVar.Units   = 'metre'
receiverPositionVar.Type    = 'cartesian'
receiverPositionVar[:]      = np.zeros((r,c,i))

samplingRateVar =   rootgrp.createVariable('Data.SamplingRate', 'f8',   ('I'))
samplingRateVar.Units = 'hertz'
samplingRateVar[:] = 48000

delayVar        =   rootgrp.createVariable('Data.Delay',        'f8',   ('I','R'))
delay = np.zeros((i,r))
delayVar[:,:] = delay

dataIRVar =         rootgrp.createVariable('Data.IR', 'f8', ('M','R','N'))
dataIRVar.ChannelOrdering   = 'acn'
dataIRVar.Normalization     = 'sn3d'
dataIRVar[:] = np.random.rand(m,r,n)

#----------Close it----------#

rootgrp.close()