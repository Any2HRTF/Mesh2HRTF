#                                Mesh2HRTF
#                Copyright (C) 2015 by Harald Ziegelwanger,
#        Acoustics Research Institute, Austrian Academy of Sciences
#                        mesh2hrtf.sourceforge.net
#
# Mesh2HRTF is licensed under the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# Mesh2HRTF is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.
# You should have received a copy of the GNU LesserGeneral Public License along with Mesh2HRTF. If not, see <http://www.gnu.org/licenses/lgpl.html>.
#
# If you use Mesh2HRTF:
# - Provide credits:
#   "Mesh2HRTF, H. Ziegelwanger, ARI, OEAW (mesh2hrtf.sourceforge.net)"
# - In your publication, cite both articles:
#   [1] Ziegelwanger, H., Kreuzer, W., and Majdak, P. (2015). "Mesh2HRTF: Open-source software package for the numerical calculation of head-related transfer functions," in Proceedings of the 22nd ICSV, Florence, IT.
#   [2] Ziegelwanger, H., Majdak, P., and Kreuzer, W. (2015). "Numerical calculation of listener-specific head-related transfer functions and sound localization: Microphone model and mesh discretization," The Journal of the Acoustical Society of America, 138, 208-222.
import os
import csv
import numpy


def Output2HRTF_Load(foldername, filename):
    # OUTPUT2HRTF_LOAD
    #   [data,frequency]=Output2HRTF_Load(foldername,filename) loads results
    #   of the BEM-HRTF calculation.
    #
    #   Input:
    #       foldername:
    #       filename:
    #
    #   Output:
    #       data: Matrix of complex values
    #           dim1: frequency
    #           dim2: datapoints

    #%% -------------------check and initialize variables----------------------
    if 'foldername' not in locals() or 'filename' not in locals():
        error('Not enough input arguments.')


    #%% ----------------------check number of freq bins------------------------
    for ii in range(1, 1000000):
        if not os.path.exists(os.path.join(foldername, 'be.%d' % ii, filename)):
            numFrequencies = ii-1
            break

    #%% --------------------check number of header and data lines--------------
    
    idx1 = 0
    idx2 = 0
    ii = 0
    with open(os.path.join(foldername, 'be.1', filename)) as file:
        line = csv.reader(file, delimiter=' ', skipinitialspace=True)
        for li in line:
            ii += 1
            if li not in (None, ""):
                if len(li) == 3:
                    idx1 = idx2
                    idx2 = int(li[0])
                    if idx2 - idx1 == 1:
                        numHeaderlines_BE = ii - 2
                        numDatalines = sum(1 for l in line) + 2
                        break

    #%% -----------------------------load data---------------------------------
    data = numpy.zeros((numFrequencies, numDatalines), dtype=complex)
    frequency = numpy.zeros(numFrequencies)

    for ii in range(numFrequencies):
        tmpData = []
        with open(os.path.join(foldername, 'be.%d' % (ii+1), filename)) as file:
            line = csv.reader(file, delimiter=' ', skipinitialspace=True)
            for li in line:
                if line.line_num > numHeaderlines_BE:
                    tmpData.append(complex(float(li[1]), float(li[2])))

        if tmpData:
            with open(os.path.join(foldername, '..', 'fe.out', 'fe.%d' % (ii+1), 'load')) as file:
                line = csv.reader(file, delimiter=' ', skipinitialspace=True)
                for li in line:
                    if line.line_num > 2:
                        tmpFrequency = float(li[0])
            data[ii, :] = tmpData
            frequency[ii] = tmpFrequency
        else:
            data[ii, :] = numpy.nan
            frequency[ii] = numpy.nan

    return data, frequency
