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
import numpy


def Output2HRTF_ReadComputationTime(filename):
    """
    Read compuation time

    Parameters
    ----------
    filename : string

    Returns
    -------
    computation_time : numpy array
    """

    f = open(filename, "r", encoding="utf8", newline="\n")
    count = -1

    for line in f:
        if line.find('Number of frequency steps') != -1:
            idx = line.find('=')
            nSteps = int(line[idx+2])
            data = numpy.zeros((nSteps, 6), dtype=int)

        if line.find('Frequency') != -1:
            count += 1
            idx1 = line.find('Step')
            idx2 = line.find('Frequency')
            data[count, 0] = int(line[idx1+5:idx2-2])
            idx1 = line.find('=')
            idx2 = line.find('Hz')
            data[count, 1] = int(line[idx1+2:idx2-1])

        if line.find('Assembling the equation system  ') != -1:
            idx = line.find(':')
            data[count, 2] = int(line[idx+2:-1])

        if line.find('Solving the equation system  ') != -1:
            idx = line.find(':')
            data[count, 3] = int(line[idx+2:-1])

        if line.find('Post processing  ') != -1:
            idx = line.find(':')
            data[count, 4] = int(line[idx+2:-1])

        if line.find('Total  ') != -1:
            idx = line.find(':')
            data[count, 5] = int(line[idx+2:-1])

        if line.find('Address computation ') != -1:
            return data
