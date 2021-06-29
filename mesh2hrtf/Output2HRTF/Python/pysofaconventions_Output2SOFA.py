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
rootgrp.SOFAConventions = 'GeneralTF'
rootgrp.SOFAConventionsVersion = '1.0'
rootgrp.APIName = 'pysofaconventions'
rootgrp.APIVersion = '0.1'
rootgrp.AuthorContact = 'andres.perez@eurecat.org'
rootgrp.Comment = ''
rootgrp.DataType = 'TF'
rootgrp.License = 'No license provided, ask the author for permission'
rootgrp.RoomType = 'free field'
rootgrp.DateCreated = time.ctime(time.time())
rootgrp.DateModified = time.ctime(time.time())
rootgrp.Title = ''

#rootgrp.DatabaseName = ''
#rootgrp.ListenerShortName = ''

rootgrp.ApplicationName = 'Mesh2HRTF'
rootgrp.ApplicationVersoin = Mesh2HRTF_version



#----------Required Dimensions----------#

m = evaluationGridsNumNodes
n = len(frequencies)
r = ears
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
listenerPositionVar = rootgrp.createVariable('ListenerPosition', 'f8', ('I', 'C'))
listenerPositionVar.Units = 'metre'
listenerPositionVar.Type = 'cartesian'
listenerPositionVar[:] = numpy.asarray([0, 0, 0])

emitterPositionVar  = rootgrp.createVariable('EmitterPosition', 'f8', ('E','C','I'))
emitterPositionVar.Units = 'metre'
emitterPositionVar.Type = 'cartesian'

sourcePositionVar = rootgrp.createVariable('SourcePosition', 'f8', ('M','C'))
sourcePositionVar.Units = 'degree, degree, metre'
sourcePositionVar.Type = 'spherical'
sourcePositionVar[:] = xyz[:, 1:4]

receiverPositionVar = rootgrp.createVariable('ReceiverPosition', 'f8', ('R','C','I'))
receiverPositionVar.Units = 'metre'
receiverPositionVar.Type = 'cartesian'
receiverPositionVar[:] = sourceCenter

# samplingRateVar = rootgrp.createVariable('Data.SamplingRate', 'f8', ('I'))
# samplingRateVar.Units = 'hertz'
# samplingRateVar[:] = 48000

# delayVar = rootgrp.createVariable('Data.Delay', 'f8', ('I','R'))
# delay = np.zeros((i,r))
# delayVar[:,:] = delay

# dataIRVar = rootgrp.createVariable('Data.IR', 'f8', ('M','R','N'))
# dataIRVar.ChannelOrdering = 'acn'
# dataIRVar.Normalization = 'sn3d'
# dataIRVar[:] = np.random.rand(m,r,n)

dataReal = rootgrp.createVariable('Data.Real', 'f8', ('M','R','N'))
dataReal[:,:,:] = numpy.real(pressure)
dataReal.LongName = 'pressure'
dataReal.Units = 'pascal'
dataImag = rootgrp.createVariable('Data.Imag', 'f8', ('M','R','N'))
dataImag[:,:,:] = numpy.imag(pressure)
dataImag.LongName = 'pressure'
dataImag.Units = 'pascal'
dataN = rootgrp.createVariable('N', 'f8', ('N'))
dataN[:] = numpy.asarray(frequencies)
dataN.Units = 'hertz'

#----------Close it----------#

rootgrp.close()