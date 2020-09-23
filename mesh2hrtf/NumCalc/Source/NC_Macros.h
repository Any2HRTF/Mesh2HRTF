/*===========================================================================*\
 *                                                                            *
 *                                 Mesh2HRTF                                  *
 *                Copyright (C) 2015 by Harald Ziegelwanger,                  *
 *        Acoustics Research Institute, Austrian Academy of Sciences          *
 *                        mesh2hrtf.sourceforge.net                           *
 *                                                                            *
 *--------------------------------------------------------------------------- *
 *                                                                            *
 *  Mesh2HRTF is licensed under the GNU Lesser General Public License as      *
 *  published by the Free Software Foundation, either version 3 of            *
 *  the License, or (at your option) any later version.                       *
 *                                                                            *
 *  Mesh2HRTF is distributed in the hope that it will be useful,              *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of            *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the              *
 *  GNU Lesser General Public License for more details.                       *
 *                                                                            *
 *  You should have received a copy of the GNU LesserGeneral Public           *
 *  License along with Mesh2HRTF. If not, see                                 *
 *  <http://www.gnu.org/licenses/lgpl.html>.                                  *
 *                                                                            *
 *  If you use Mesh2HRTF:                                                     *
 *  - Provide credits:                                                        *
 *    "Mesh2HRTF, H. Ziegelwanger, ARI, OEAW (mesh2hrtf.sourceforge.net)"     *
 *  - In your publication, cite both articles:                                *
 *    [1] Ziegelwanger, H., Kreuzer, W., and Majdak, P. (2015). "Mesh2HRTF:   *
 *        Open-source software package for the numerical calculation of       *
 *        head-related transfer functions," in Proceedings of the 22nd        *
 *        ICSV, Florence, IT.                                                 *
 *    [2] Ziegelwanger, H., Majdak, P., and Kreuzer, W. (2015). "Numerical    *
 *        calculation of listener-specific head-related transfer functions    *
 *        and sound localization: Microphone model and mesh discretization,"  *
 *        The Journal of the Acoustical Society of America, 138, 208-222.     *
 *                                                                            *
 \*===========================================================================*/



/*===========================================================================*\
 *                                                                            *
 *  File name:      NC_Macros.h                                               *
 *  Description:    Definitions of macros                                     *
 *  Author:         H. Ziegelwanger, W. Kreuzer and Z.S. Chen                 *
 *                                                                            *
 \*===========================================================================*/



#ifndef NC_Macros_h
#define NC_Macros_h



//================================= I N C L U D E S ==================================
// local includes                                                                   //
//                                                                                  //
// system includes                                                                  //
#include<iostream>                                                                  //
#include<iomanip>                                                                   //
#include<math.h>                                                                    //
using namespace std;                                                                //
//====================================================================================



// maximum
#define max_(a, b) a > b ? a : b

// minimum
#define min_(a, b) a < b ? a : b

// norm of a real vector of dimension 3, normization of it
#define Vcnorm_dim3_(V, L) L = sqrt(V[0]*V[0] + V[1]*V[1] + V[2]*V[2]); V[0] /= L; V[1] /= L; V[2] /= L

// distance between two points
#define Dispoi_dim3_(P1, P2) sqrt((P1[0] - P2[0])*(P1[0] - P2[0]) + (P1[1] - P2[1])*(P1[1] - P2[1]) + (P1[2] - P2[2])*(P1[2] - P2[2]))

// scalar product of two vectors of dimension 3
#define Scprod_dim3_(V1, V2) (V1[0]*V2[0] + V1[1]*V2[1] + V1[2]*V2[2])

// error of integral of the kern function G over an element
#define Err_G_(m, x, i, err) err = 32.0; for(i=0; i<2*m+1; i++) {err *= x;}

// error of integral of the kern function H or H(T) over an element
#define Err_H_(m, x, i, err) err = 64.0*(2.0*(double)m + 1.0); for(i=0; i<2*m+2; i++) {err *= x;}

// error of integral of the kern function E over an element
#define Err_E_(m, x, i, err) err = 128.0*((double)m + 1.0)*(2.0*(double)m + 1.0); for(i=0; i<2*m+3; i++) {err *= x;}



#endif