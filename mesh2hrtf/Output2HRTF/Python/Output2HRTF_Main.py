#                                Mesh2HRTF
#                Copyright (C) 2015 by Harald Ziegelwanger,
#        Acoustics Research Institute, Austrian Academy of Sciences
#                        mesh2hrtf.sourceforge.net
#
# Mesh2HRTF is licensed under the GNU Lesser General Public License as
# published by the Free Software Foundation, either version 3 of the License,
# or (at your option) any later version. Mesh2HRTF is distributed in the hope
# that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details. You should have received a
# copy of the GNU LesserGeneral Public License along with Mesh2HRTF. If not,
# see <http://www.gnu.org/licenses/lgpl.html>.
#
# If you use Mesh2HRTF:
# - Provide credits:
#   "Mesh2HRTF, H. Ziegelwanger, ARI, OEAW (mesh2hrtf.sourceforge.net)"
# - In your publication, cite both articles:
#   [1] Ziegelwanger, H., Kreuzer, W., and Majdak, P. (2015). "Mesh2HRTF:
#       Open-source software package for the numerical calculation of
#       head-related transfer functions," in Proceedings of the 22nd ICSV,
#       Florence, IT.
#   [2] Ziegelwanger, H., Majdak, P., and Kreuzer, W. (2015). "Numerical
#       calculation of listener-specific head-related transfer functions and
#       sound localization: Microphone model and mesh discretization," The
#       Journal of the Acoustical Society of America, 138, 208-222.
#
# Author: Harald Ziegelwanger
#        (Acoustics Research Institute, Austrian Academy of Sciences)
#        Fabian Brinkmann, Robert Pelzer, Jeffrey Thomsen
#        (Audio Communication Group, Technical University Berlin)

# header
import os
import numpy
# TODO: I would suggest to move the functions into this module
import Output2HRTF_ReadComputationTime as o2hrtf_rct
import Output2HRTF_Load as o2hrtf_l
import sofa
from datetime import date
today = date.today()

# TODO: The module would be easier to test and read if we would separate it
#       into smaller functions. I've added empty function definitions at the
#       end for discussion.


def Output2HRTF_Main(
        projectPath, Mesh2HRTF_version, cpusAndCores, sourceType, sourceCenter,
        sourceArea, reference, computeHRIRs, speedOfSound, densityOfAir):
    """
    Process NumcCalc output and write data as SOFA files.

    All parameters are written upon exporting a Mesh2HRTF project from Blender.

    Parameters
    ----------
    projectPath : string
        Path to the Mesh2HRTF folder.
    Mesh2HRTF_version : string
        Mesh2HRTF version as string
    cpusAndCores : array
        Array indicating which CPUs and Cores were used for the BEM
    sourceType : string
        The source type (vibrating element, point source)
    sourceCenter : array
        The source position
    sourceArea : array
        The area of the source(s) if using vibrating elements as sourceType.
        This is required for referencing the HRTFs
    reference : boolean
        Indicate if the HRTF are referenced to the pressure in the center of
        the head with the head absent.
    computeHRIRs : boolean
        Indicate if the HRIRs should be calculated by means of the inverse
        Fourier transform.
    speedOfSound : float
        The speed of sound in m/s
    densityOfAir : float
        The density of air in kg/m^3
    """

    # ---------------------------- load meta data -----------------------------
    # output directory
    os.chdir(projectPath)
    if not os.path.exists(os.path.join(os.getcwd(), 'Output2HRTF')):
        os.makedirs(os.path.join(os.getcwd(), 'Output2HRTF'))

    # number of ears
    ears = numpy.amax(cpusAndCores)

    # get the evaluation grids
    print('\n Loading the EvaluationGrids ...')
    evaluationGrids = []
    evalGridsList = os.listdir('EvaluationGrids')
    evaluationGridsNumNodes = 0
    for ii in range(len(evalGridsList)):
        tmpNodes = numpy.loadtxt(os.path.join(
            'EvaluationGrids', evalGridsList[ii], 'Nodes.txt'),
            delimiter=' ', skiprows=1)
        tmpElements = numpy.loadtxt(os.path.join(
            'EvaluationGrids', evalGridsList[ii], 'Elements.txt'),
            delimiter=' ', skiprows=1)
        evaluationGrids.append({"name": evalGridsList[ii],
                                "nodes": tmpNodes,
                                "elements": tmpElements,
                                "num_nodes": tmpNodes.shape[0]})
        evaluationGridsNumNodes += evaluationGrids[ii]['num_nodes']

    # get the object mesh
    print('\n Loading the ObjectMeshes ...')
    objectMeshes = []
    objMeshesList = os.listdir('ObjectMeshes')
    objectMeshesNumNodes = 0
    for ii in range(len(objMeshesList)):
        tmpNodes = numpy.loadtxt(os.path.join(
            'ObjectMeshes', objMeshesList[ii], 'Nodes.txt'),
            delimiter=' ', skiprows=1)
        tmpElements = numpy.loadtxt(os.path.join(
            'ObjectMeshes', objMeshesList[ii], 'Elements.txt'),
            delimiter=' ', skiprows=1)
        objectMeshes.append({"name": objMeshesList[ii],
                             "nodes": tmpNodes,
                             "elements": tmpElements,
                             "num_nodes": tmpNodes.shape[0]})
        objectMeshesNumNodes += objectMeshes[ii]['num_nodes']

    del ii, tmpNodes, tmpElements

    # Read computational effort
    print('\n Loading computational effort data ...')
    # TODO: remove prints in the for loop. Loading this is fast and printing
    #       only clutters the command line in this case...
    for ch in range(ears):
        computationTime = []
        print('\n    Ear %d ...' % (ch+1))
        for ii in range(cpusAndCores.shape[0]):
            print('\n        CPU %d: ' % (ii+1))
            for jj in range(cpusAndCores.shape[1]):
                if (cpusAndCores[ii, jj] == ch+1):
                    print('%d, ' % (jj+1))
                    tmpFilename = os.path.join(
                        'NumCalc', 'CPU_%d_Core_%d' % ((ii+1), (jj+1)),
                        'NC.out')
                    tmp = o2hrtf_rct.Output2HRTF_ReadComputationTime(
                        tmpFilename)
                    computationTime.append(tmp)
                    del tmp
            print('...')

    description = ['Frequency index', 'Frequency', 'Building', 'Solving',
                   'Postprocessing', 'Total']

# TODO: This is just informal to compare or report computation times. Saving
#       it as a text file is fine and the most generic format. I would suggest
#       a csv file:
#       numpy.savetext("filename.csv", computationTime, delimiter=", ",
#                      header=", ".join(description), comments="")
# this following file save operation needs to be refined, perhaps turn it into
# a dictionary or save as a .npy array, let's see - I don't even know what this
#  file is used for
#    file = open(os.path.join('Output2HRTF', 'computationTime.txt'), "w")
#    file.write(repr(description)+'\n')
#    file.write(numpy.array2string(
#        computationTime, separator=',', formatter='int')+'\n')
#    file.close()

    del ch, ii, jj, description, computationTime

    # Load ObjectMesh data
    tmpPressure = []
    frequencies = []
    pressure = []
    print('\n Loading ObjectMesh data ...')
    for ch in range(ears):
        print('\n    Ear %d ...' % (ch+1))
        for ii in range(cpusAndCores.shape[0]):
            print('\n        CPU %d: ' % (ii+1))
            for jj in range(cpusAndCores.shape[1]):
                if (cpusAndCores[ii, jj] == (ch+1)):
                    print('%d, ' % (jj+1))
                    tmpFilename = os.path.join(
                        'NumCalc', 'CPU_%d_Core_%d' % ((ii+1), (jj+1)),
                        'be.out')
                    tmpData, tmpFrequencies = o2hrtf_l.Output2HRTF_Load(
                        tmpFilename, 'pBoundary')
                    if tmpPressure:
                        tmpPressure.append(tmpData)
                        frequencies.append(tmpFrequencies)
                    else:
                        tmpPressure = tmpData
                        frequencies = tmpFrequencies
                    del tmpData, tmpFrequencies, tmpFilename
            print('...')
        idx = sorted(range(len(frequencies)), key=lambda k: frequencies[k])
        frequencies = sorted(frequencies)
        pressure.append(tmpPressure[idx, :])

        del tmpPressure
    pressure = numpy.transpose(numpy.array(pressure), (2, 0, 1))

    print('\nSave ObjectMesh data ...')
    cnt = 0
    for ii in range(len(objectMeshes)):
        nodes = objectMeshes[ii]["nodes"]
        elements = objectMeshes[ii]["elements"]
        element_data = pressure

        # old version using list representation of pressure
        # for jj in range(len(pressure)):
        #     element_data[jj] = element_data[jj][:, cnt:cnt+elements.shape[0]]
        # new version using MRN-transposed array representation of pressure
        for jj in range(pressure.shape[1]):
            element_data[:, jj, :] = element_data[
                cnt:cnt+elements.shape[0], jj, :]

        # probably still need a better way of storing these values
        # TODO: Saving in numpy format is fine, but I would suggest to use
        #       savez_compressed :)
        file = open("ObjectMesh_"+objectMeshes[ii]["name"]+".npz", "w")
        numpy.savez(file.name, nodes=nodes, elements=elements,
                    frequencies=frequencies, element_data=element_data)
        file.close()

        # file = open("ObjectMeshAppend_"+objectMeshes[ii]["name"]+".txt", "a")
        # numpy.savetxt(file, nodes)
        # file.write('\n')
        # numpy.savetxt(file, elements)
        # file.write('\n')
        # numpy.savetxt(file, frequencies)
        # file.write('\n')
        # numpy.savetxt(file, element_data)
        # file.write('\n')
        # file.close()

        # file = open("ObjectMesh_"+objectMeshes[ii]["name"]+".txt", "w")
        # file.write(repr(nodes)+'\n')
        # file.write(repr(elements)+'\n')
        # file.write(repr(frequencies)+'\n')
        # file.write(repr(element_data)+'\n')
        # file.close()

        cnt = cnt + elements.shape[0]

    del pressure, nodes, elements, frequencies, ii, jj, cnt, ch, idx, \
        element_data

    # Load EvaluationGrid data
    tmpPressure = []
    frequencies = []
    pressure = []
    if not len(evaluationGrids) == 0:
        print('\nLoading data for the evaluation grids ...')
        for ch in range(ears):
            print('\n    Ear %d ...' % (ch+1))
            for ii in range(cpusAndCores.shape[0]):
                print('\n        CPU %d: ' % (ii+1))
                for jj in range(cpusAndCores.shape[1]):
                    if (cpusAndCores[ii, jj] == (ch+1)):
                        print('%d, ' % (jj+1))
                        tmpFilename = os.path.join(
                            'NumCalc', 'CPU_%d_Core_%d' % ((ii+1), (jj+1)),
                            'be.out')
                        [tmpData, tmpFrequencies] = o2hrtf_l.Output2HRTF_Load(
                            tmpFilename, 'pEvalGrid')
                        if tmpPressure:
                            tmpPressure.append(tmpData)
                            frequencies.append(tmpFrequencies)
                        else:
                            tmpPressure = tmpData
                            frequencies = tmpFrequencies
                        del tmpData, tmpFrequencies
                print('...')
            pressure.append(tmpPressure)
            del tmpPressure
        idx = sorted(range(len(frequencies)), key=lambda k: frequencies[k])
        frequencies = sorted(frequencies)
        pressure = [i[idx, :] for i in pressure]
    pressure = numpy.transpose(numpy.array(pressure), (2, 0, 1))

    # save to struct
    cnt = 0
    for ii in range(len(evaluationGrids)):
        # old version using list representation of pressure
        # evaluationGrids[ii]["pressure"] = \
        #     [i[:, cnt:cnt+evaluationGrids[ii]["num_nodes"],] \
        #     for i in pressure]
        # new version using MRN-transposed array representation of pressure
        evaluationGrids[ii]["pressure"] = pressure[
            cnt:cnt+evaluationGrids[ii]["num_nodes"], :, :]

        cnt = cnt + evaluationGrids[ii]["num_nodes"]

    del ch, ii, jj, pressure, cnt, idx

    # reference to pressure in the middle of the head with the head absent
    # according to the HRTF definition.
    if reference:

        # this might be a parameter in the function call
        # 1: reference to only one radius (the smallest found)
        # 2: reference to all individual radii
        refMode = 1

        for ii in range(len(evaluationGrids)):

            xyz = evaluationGrids[ii]["nodes"]
            pressure = evaluationGrids[ii]["pressure"]
            freqMatrix = numpy.tile(
                frequencies, (pressure.shape[0], pressure.shape[1], 1))

            # distance of source positions from the origin
            if refMode == 1:
                r = min(numpy.sqrt(xyz[:, 1]**2 + xyz[:, 2]**2 + xyz[:, 3]**2))
                r = numpy.tile(
                    r,
                    (pressure.shape[0], pressure.shape[1], pressure.shape[2]))
            else:
                r = numpy.sqrt(xyz[:, 1]**2 + xyz[:, 2]**2 + xyz[:, 3]**2)
                r = numpy.tile(numpy.transpose(r),
                               (pressure.shape[0], pressure.shape[1], 1))

            if sourceType == 'vibratingElement':

                volumeFlow = 0.1 * numpy.ones(
                    (pressure.shape[0],
                     pressure.shape[1],
                     pressure.shape[2]))
                if 'sourceArea':
                    # has to be fixed for both ears....
                    for nn in range(len(sourceArea)):
                        volumeFlow[:, nn, :] = \
                            volumeFlow[:, nn, :] * sourceArea[nn]

                # point source in the origin evaluated at r
                # eq. (6.71) in: Williams, E. G. (1999). Fourier Acoustics.
                ps = -1j * densityOfAir * 2 * numpy.pi * freqMatrix * \
                    volumeFlow / (4 * numpy.pi) * \
                    numpy.exp(1j * 2 * numpy.pi *
                              freqMatrix / speedOfSound * r) / \
                    r

            elif sourceType == 'pointSource':

                amplitude = 0.1  # hard coded in Mesh2HRTF
                ps = amplitude * \
                    numpy.exp(1j * 2 * numpy.pi * freqMatrix /
                              speedOfSound * r) / \
                    (4 * numpy.pi * r)

            else:
                raise ValueError(
                    ("Referencing is currently only implemented for "
                     "sourceType 'vibratingElement' and 'pointSource'."))

            # here we go...
            evaluationGrids[ii]["pressure"] = pressure / ps

        del r, freqMatrix, ps, ii

    # Save data as SOFA file
    print('\nSaving complex pressure to SOFA file ...\n')

    for ii in range(len(evaluationGrids)):

        xyz = evaluationGrids[ii]["nodes"]
        pressure = evaluationGrids[ii]["pressure"]

        # Save as GeneralTF
        TF_path = os.path.join(
            'Output2HRTF', 'HRTF_%s.sofa' % evaluationGrids[ii]["name"])

        Obj = sofa.Database.create(TF_path, "GeneralTF",
            dimensions={"M": evaluationGridsNumNodes, "R": ears,
            "N": len(frequencies)})

        Obj.Metadata.set_attribute('GLOBAL_ApplicationName', 'Mesh2HRTF')
        Obj.Metadata.set_attribute(
            'GLOBAL_ApplicationVersion', Mesh2HRTF_version)
        Obj.Metadata.set_attribute('GLOBAL_Organization', '')
        Obj.Metadata.set_attribute('GLOBAL_Title', '')
        Obj.Metadata.set_attribute(
            'GLOBAL_DateCreated', today.strftime("%b-%d-%Y"))
        Obj.Metadata.set_attribute('GLOBAL_AuthorContact', '')

        Obj.Listener.initialize(fixed=["Position", "View", "Up"])
        Obj.Listener.Position = [0, 0, 0]

        Obj.Receiver.initialize(fixed=["Position"], count=ears)
        Obj.Receiver.Position = sourceCenter

        Obj.Source.initialize(variances=["Position"], fixed=["View", "Up"])
        Obj.Source.Position.set_values(xyz[:, 1:4])
        
        Obj.Emitter.initialize(fixed=["Position", "View", "Up"], count=1)

        Obj.Data.Type = 'TF'
        Obj.Data.initialize()
        Obj.Data.Real.set_values(numpy.real(pressure))
        Obj.Data.Imag.set_values(numpy.imag(pressure))
        Obj.Data.create_attribute('Real_LongName', 'pressure')
        Obj.Data.create_attribute('Real_Units', 'pascal')
        Obj.Data.create_attribute('Imag_LongName', 'pressure')
        Obj.Data.create_attribute('Imag_Units', 'pascal')
        Obj.Data.N.set_values(frequencies)

    #     Obj.SourcePosition_Type='cartesian'
    #     Obj.SourcePosition_Units='meter'
        Obj.close()

    del Obj, ii, xyz, pressure


    #%% Save time data data as SOFA file
    if computeHRIRs:

        print('\nSaving time data to SOFA file ...\n')

        for ii in range(len(evaluationGrids)):

            # check if the frequency vector has the correct format
            if not all(numpy.abs(frequencies[0] - numpy.diff(frequencies)) < .1):
                raise ValueError('The frequency vector must be if the format a:a:fs/2, with a>0 and fs the sampling rate.')

            if not reference:
                raise ValueError('HRIRs can only be computet if reference=true')

            xyz = evaluationGrids[ii]["nodes"]
            pressure = evaluationGrids[ii]["pressure"]

            fs = 2*frequencies[-1]

            # add 0 Hz bin
            pressure = numpy.concatenate((numpy.ones((pressure.shape[0],
                pressure.shape[1], 1)), pressure), axis=2)
            # make fs/2 real
            pressure[:, :, -1] = numpy.abs(pressure[:, :, -1])
            # ifft (take complex conjugate because sign conventions differ)
            hrir = numpy.fft.irfft(numpy.conj(pressure))

            # shift 30 cm to make causal
            # (path differences between the origin and the ear are usually
            # smaller than 30 cm but numerical HRIRs show stringer pre-ringing)
            n_shift = int(numpy.round(.30 / (1/fs * speedOfSound)))
            hrir = numpy.roll(hrir, n_shift, axis=2)

            # Save as GeneralFIR
            FIR_path = os.path.join(
                'Output2HRTF', 'HRIR_%s.sofa' % evaluationGrids[ii]["name"])

            Obj = sofa.Database.create(FIR_path, "GeneralFIR",
                dimensions={"M": evaluationGridsNumNodes, "R": ears,
                "N": len(frequencies)})

            Obj.Metadata.set_attribute('GLOBAL_ApplicationName', 'Mesh2HRTF')
            Obj.Metadata.set_attribute(
                'GLOBAL_ApplicationVersion', Mesh2HRTF_version)
            Obj.Metadata.set_attribute('GLOBAL_Organization', '')
            Obj.Metadata.set_attribute('GLOBAL_Title', '')
            Obj.Metadata.set_attribute(
                'GLOBAL_DateCreated', today.strftime("%b-%d-%Y"))
            Obj.Metadata.set_attribute('GLOBAL_AuthorContact', '')

            Obj.Listener.initialize(fixed=["Position", "View", "Up"])
            Obj.Listener.Position = [0, 0, 0]

            Obj.Receiver.initialize(fixed=["Position"], count=ears)
            Obj.Receiver.Position = sourceCenter

            Obj.Source.initialize(variances=["Position"], fixed=["View", "Up"])
            Obj.Source.Position.set_values(xyz[:, 1:4])
            
            Obj.Emitter.initialize(fixed=["Position", "View", "Up"], count=1)

            Obj.Data.Type = 'FIR'
            Obj.Data.initialize()
            Obj.Data.IR = hrir
            Obj.Data.SamplingRate = fs

        #     Obj.SourcePosition_Type='cartesian'
        #     Obj.SourcePosition_Units='meter'
            Obj.close()

        del Obj, ii, xyz, pressure

    print('Done\n')


def read_nodes_and_elements():
    """
    Read the nodes and elements of the evaluation grids or object meshes.
    """
    # TODO: The code for reading the data above is very redundant. If this
    #       function uses the code with an additional parameter
    #       data='ObjectMeshes' or data='EvaluationGrids' one function might
    #       read both :)
    pass


def read_computation_time():
    """Read the computation time for each frequency."""
    pass


def read_pressure():
    """Read the sound pressure on the object meshes or evaluation grid."""
    # TODO: The code for reading the data above is very redundant. If this
    #       function uses the code with an additional parameter
    #       data='pEvalGrid' or data='pBoundary' one function might read
    #       both :)
    pass


def read_pressure_on_evaluation_grids():
    """read the sound pressure on the evaluation grids."""
    pass


def reference_HRTF():
    """Reference HRTF to the sound pressure in the center of the head."""
    pass


def compute_HRIR():
    """Compute HRIR from HRTF by means of the inverse Fourier transform"""
    pass


def write_to_sofa():
    # TODO: It might be possible to use this function to write HRTFs and HRIRs
    #       The two things that might differ are the convention and the data
    #       which could be passed as parameters
    pass
