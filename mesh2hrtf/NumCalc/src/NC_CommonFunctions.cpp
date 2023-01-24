/*===========================================================================*\
 *                                                                            *
 *  File name:      NC_CommonFunctions.cpp                                    *
 *  Description:    common public functions                                   *
 *  Author:         H. Ziegelwanger, W. Kreuzer and Z.S. Chen                 *
 *                                                                            *
 \*===========================================================================*/



//================================= I N C L U D E S ==================================
// local includes                                                                   //
#include"NC_ConstantsVariables.h"                                                   //
#include"NC_TypeDefinition.h"                                                       //
#include"NC_Arrays.h"                                                               //
#include"NC_IntegrationConstants.h"                                                 //
//                                                                                  //
// system includes                                                                  //
#include<iostream>                                                                  //
#include<iomanip>                                                                   //
#include<string>                                                                    //
#include<fstream>  // file class functions                                          //
#include<cstdlib> // c++ p625                                                       //
#include<math.h>                                                                    //
//====================================================================================



//  Local functions
void NC_MatrixVectorMultiplication(Vector<Complex>&, Vector<Complex>&);
void NC_ComputeScalingVector(ofstream&);
void NC_IncompleteLUDecomposition(ofstream&);
void NC_DeleteIncompleteLUMatrices();
void NC_IncompleteLUForBack(ofstream&, Vector<Complex>&, Vector<Complex>&);
void NC_IncompleteLUForBack();



double ErroIterSols = 1.0e-9;
//int niter_max = 1500;
double *dscaling;// scaling array scaling preconditiner
Complex *zL_incom, *zU_incom;// arrays for incomplete LU-decompostion preconditioner
int *nL_colnu, *nL_firnu, *nU_rownu, *nU_firnu;



// print erro message without output parameters
void NC_Error_Exit_0 // variant 0_0
(
	ofstream& NCout,
	const string& errorMessage
)
{
	cout << errorMessage << endl;
	NCout << errorMessage << endl;
	exit(1);
}

void NC_Error_Warning_0 // variant 1_0
(
	ofstream& NCout,
	const string& errorMessage
)
{
	cout << errorMessage << "\n" << endl;
	NCout << errorMessage << "\n" << endl;
}

void NC_Error_Exit_0 // variant 0_1
(
	const string& errorMessage
)
{
	cout << errorMessage << endl;
	exit(1);
}

void NC_Error_Warning_0 // variant 1_1
(
	const string& errorMessage
)
{
	cout << errorMessage << "\n" << endl;
}

//print erro message with one parameter with 1 out put parameter
void NC_Error_Exit_1 // variant 0_0
(
	ofstream& NCout,
	const string& errorMessage,
	const string& str1,
	const int& int1
)
{
	cout << errorMessage << endl;
	cout << str1 << int1 << endl;
	NCout << errorMessage << endl;
	NCout << str1 << int1 << endl;
	exit(1);
}

void NC_Error_Exit_1 // variant 1_0
(
	ofstream& NCout,
	const string& errorMessage,
	const string& str1,
	const double& d1
)
{
	cout << errorMessage << endl;
	cout << str1 << d1 << endl;
	NCout << errorMessage << endl;
	NCout << str1 << d1 << endl;
	exit(1);
}

void NC_Error_Exit_1 // variant 2_0
(
    ofstream& NCout,
    const string& errorMessage,
    const string& str1
)
{
    cout << errorMessage << " " << str1 << endl;
    NCout << errorMessage << " " << str1 << endl;
	exit(1);
}

void NC_Error_Warning_1 // variant 3_0
(
	ofstream& NCout,
	const string& errorMessage,
	const string& str1,
	const int& int1
)
{
	cout << "\nWarning: " << errorMessage << endl;
	cout << str1 << int1 << "\n" << endl;
	NCout << "\nWarning: " << errorMessage << endl;
	NCout << str1 << int1 << "\n" << endl;
}

void NC_Error_Warning_1 // variant 4_0
(
	ofstream& NCout,
	const string& errorMessage,
	const string& str1,
	const double& d1
)
{
	cout << "\nWarning: " << errorMessage << endl;
	cout << str1 << d1 << "\n" << endl;
	NCout << "\nWarning: " << errorMessage << endl;
	NCout << str1 << d1 << "\n" << endl;
}

void NC_Error_Exit_1 // variant 0_1
(
	const string& errorMessage,
	const string& str1,
	const int& int1
)
{
	cout << errorMessage << endl;
	cout << str1 << int1 << endl;
	exit(1);
}

void NC_Error_Exit_1 // variant 1_1
(
	const string& errorMessage,
	const string& str1,
	const double& d1
)
{
	cout << errorMessage << endl;
	cout << str1 << d1 << endl;
	exit(1);
}

void NC_Error_Exit_1 // variant 2_1
(
    const string& errorMessage,
    const string& str1
)
{
	cout << errorMessage << " " << str1 << endl;
	exit(1);
}

void NC_Error_Warning_1 // variant 3_1
(
	const string& errorMessage,
	const string& str1,
	const int& int1
)
{
	cout << "\nWarning: " << errorMessage << endl;
	cout << str1 << int1 << "\n" << endl;
}

void NC_Error_Warning_1 // variant 4_1
(
	const string& errorMessage,
	const string& str1,
	const double& d1
)
{
	cout << "\nWarning: " << errorMessage << endl;
	cout << str1 << d1 << "\n" << endl;
}

// print erro message with two out put parameters
void NC_Error_Exit_2 // variant 0_0
(
	ofstream& NCout,
	const string& errorMessage,
	const string& str1,
	const int& int1,
	const string& str2,
	const int& int2
)
{
	cout << errorMessage << endl;
	cout << str1 << int1 << endl;
	cout << str2 << int2 << endl;
	NCout << errorMessage << endl;
	NCout << str1 << int1 << endl;
	NCout << str2 << int2 << endl;
	exit(1);
}

void NC_Error_Exit_2 // variant 1_0
(
	ofstream& NCout,
	const string& errorMessage,
	const string& str1,
	const double& d1,
	const string& str2,
	const double& d2
)
{
	cout << errorMessage << endl;
	cout << str1 << d1 << endl;
	cout << str2 << d2 << endl;
	NCout << errorMessage << endl;
	NCout << str1 << d1 << endl;
	NCout << str2 << d2 << endl;
	exit(1);
}

void NC_Error_Warning_2 // variant 2_0
(
	ofstream& NCout,
	const string& errorMessage,
	const string& str1,
	const int& int1,
	const string& str2,
	const int& int2
)
{
	cout << "\nWarning: " << errorMessage << endl;
	cout << str1 << int1 << endl;
	cout << str2 << int2 << "\n" << endl;
	NCout << "\nWarning: " << errorMessage << endl;
	NCout << str1 << int1 << endl;
	NCout << str2 << int2 << "\n" << endl;
}

void NC_Error_Warning_2 // variant 3_0
(
	ofstream& NCout,
	const string& errorMessage,
	const string& str1,
	const double& d1,
	const string& str2,
	const double& d2
)
{
	cout << "\nWarning: " << errorMessage << endl;
	cout << str1 << d1 << endl;
	cout << str2 << d2 << "\n" << endl;
	NCout << "\nWarning: " << errorMessage << endl;
	NCout << str1 << d1 << endl;
	NCout << str2 << d2 << "\n" << endl;
}

void NC_Error_Warning_2 // variant 4_0
(
	ofstream& NCout,
	const string& errorMessage,
	const string& str1,
	const int& i1,
	const string& str2,
	const double& d2
)
{
	cout << "\nWarning: " << errorMessage << endl;
	cout << str1 << i1 << endl;
	cout << str2 << d2 << "\n" << endl;
	NCout << "\nWarning: " << errorMessage << endl;
	NCout << str1 << i1 << endl;
	NCout << str2 << d2 << "\n" << endl;
}


void NC_Error_Exit_2 // variant 0_1
(
	const string& errorMessage,
	const string& str1,
	const int& int1,
	const string& str2,
	const int& int2
)
{
	cout << errorMessage << endl;
	cout << str1 << int1 << endl;
	cout << str2 << int2 << endl;
	exit(1);
}

void NC_Error_Exit_2 // variant 1_1
(
	const string& errorMessage,
	const string& str1,
	const double& d1,
	const string& str2,
	const double& d2
)
{
	cout << errorMessage << endl;
	cout << str1 << d1 << endl;
	cout << str2 << d2 << endl;
	exit(1);
}

void NC_Error_Warning_2 // variant 2_1
(
	const string& errorMessage,
	const string& str1,
	const int& int1,
	const string& str2,
	const int& int2
)
{
	cout << "\nWarning: " << errorMessage << endl;
	cout << str1 << int1 << endl;
	cout << str2 << int2 << "\n" << endl;
}

void NC_Error_Warning_2 // variant 3_1
(
	const string& errorMessage,
	const string& str1,
	const double& d1,
	const string& str2,
	const double& d2
)
{
	cout << "\nWarning: " << errorMessage << endl;
	cout << str1 << d1 << endl;
	cout << str2 << d2 << "\n" << endl;
}

// Operator << for complex values
ostream& operator<<(ostream &s, const Complex &u)
{
    s << "(" << u.re() << ", " << u.im() << ")";
    return s;
};

// compute distance between two points, return: the result
double BLdispoi(Vector<double>& pointa, Vector<double>& pointb)
{
	double dist = 0.0, dif;
	int i;

	if(pointa.size() < NDIM || pointb.size() < NDIM)
		NC_Error_Exit_0(NCout, "Dimension of the vectors must be not less than NDIM!");

	for(i=0; i<NDIM; i++)
	{
		dif = pointa[i] - pointb[i];
		dist += dif*dif;
	}
	dist = sqrt(dist);

	return dist;
}

// compute square root of a complex number
Complex BLzsqrt(const Complex& z)
{
	Complex zsq;
	double re1, ph1;

	if(z.norm() == 0.0)
	{
		zsq.set(0.0, 0.0);
	}
	else
	{
		re1 = sqrt(z.norm());
		ph1 = atan2(z.im(), z.re())/2.0;
		zsq.set(re1*cos(ph1), re1*sin(ph1));
	}

	return zsq;
}

// compute the coordinates and weights of a given number of the Gauss integration points
void BLGauCooWe
(
	ofstream& NCout,
	const int& npt,          // number of Gaussean points
	Vector<double>& crdgp,   // O: coordinates of the Gaussean points
	Vector<double>& weigp    // O: weights of the Gaussean points
)
{
	int i, hpt = npt/2, hp1 = hpt - 1;

	switch(npt)
	{
	case 1:
		crdgp[0] = GCoorWei_1[0][0];
		weigp[0] = GCoorWei_1[0][1];
		break;
	case 2:
		for(i=0; i<hpt; i++) {
			crdgp[hpt + i] = GCoorWei_2[i][0];
			crdgp[hp1 - i] = -crdgp[hpt + i];
			weigp[hpt + i] = weigp[hp1 - i] = GCoorWei_2[i][1];
		}
		break;
	case 3:
		crdgp[hpt] = GCoorWei_3[0][0];
		weigp[hpt] = GCoorWei_3[0][1];
		for(i=1; i<=hpt; i++) {
			crdgp[hpt + i] = GCoorWei_3[i][0];
			crdgp[hpt - i] = -crdgp[hpt + i];
			weigp[hpt + i] = weigp[hpt - i] = GCoorWei_3[i][1];
		}
		break;
	case 4:
		for(i=0; i<hpt; i++) {
			crdgp[hpt + i] = GCoorWei_4[i][0];
			crdgp[hp1 - i] = -crdgp[hpt + i];
			weigp[hpt + i] = weigp[hp1 - i] = GCoorWei_4[i][1];
		}
		break;
	case 5:
		crdgp[hpt] = GCoorWei_5[0][0];
		weigp[hpt] = GCoorWei_5[0][1];
		for(i=1; i<=hpt; i++) {
			crdgp[hpt + i] = GCoorWei_5[i][0];
			crdgp[hpt - i] = -crdgp[hpt + i];
			weigp[hpt + i] = weigp[hpt - i] = GCoorWei_5[i][1];
		}
		break;
	case 6:
		for(i=0; i<hpt; i++) {
			crdgp[hpt + i] = GCoorWei_6[i][0];
			crdgp[hp1 - i] = -crdgp[hpt + i];
			weigp[hpt + i] = weigp[hp1 - i] = GCoorWei_6[i][1];
		}
		break;
	case 7:
		crdgp[hpt] = GCoorWei_7[0][0];
		weigp[hpt] = GCoorWei_7[0][1];
		for(i=1; i<=hpt; i++) {
			crdgp[hpt + i] = GCoorWei_7[i][0];
			crdgp[hpt - i] = -crdgp[hpt + i];
			weigp[hpt + i] = weigp[hpt - i] = GCoorWei_7[i][1];
		}
		break;
	case 8:
		for(i=0; i<hpt; i++) {
			crdgp[hpt + i] = GCoorWei_8[i][0];
			crdgp[hp1 - i] = -crdgp[hpt + i];
			weigp[hpt + i] = weigp[hp1 - i] = GCoorWei_8[i][1];
		}
		break;
	case 9:
		crdgp[hpt] = GCoorWei_9[0][0];
		weigp[hpt] = GCoorWei_9[0][1];
		for(i=1; i<=hpt; i++) {
			crdgp[hpt + i] = GCoorWei_9[i][0];
			crdgp[hpt - i] = -crdgp[hpt + i];
			weigp[hpt + i] = weigp[hpt - i] = GCoorWei_9[i][1];
		}
		break;
	case 10:
		for(i=0; i<hpt; i++) {
			crdgp[hpt + i] = GCoorWei_10[i][0];
			crdgp[hp1 - i] = -crdgp[hpt + i];
			weigp[hpt + i] = weigp[hp1 - i] = GCoorWei_10[i][1];
		}
		break;
	case 11:
		crdgp[hpt] = GCoorWei_11[0][0];
		weigp[hpt] = GCoorWei_11[0][1];
		for(i=1; i<=hpt; i++) {
			crdgp[hpt + i] = GCoorWei_11[i][0];
			crdgp[hpt - i] = -crdgp[hpt + i];
			weigp[hpt + i] = weigp[hpt - i] = GCoorWei_11[i][1];
		}
		break;
	case 12:
		for(i=0; i<hpt; i++) {
			crdgp[hpt + i] = GCoorWei_12[i][0];
			crdgp[hp1 - i] = -crdgp[hpt + i];
			weigp[hpt + i] = weigp[hp1 - i] = GCoorWei_12[i][1];
		}
		break;
	case 13:
		crdgp[hpt] = GCoorWei_13[0][0];
		weigp[hpt] = GCoorWei_13[0][1];
		for(i=1; i<=hpt; i++) {
			crdgp[hpt + i] = GCoorWei_13[i][0];
			crdgp[hpt - i] = -crdgp[hpt + i];
			weigp[hpt + i] = weigp[hpt - i] = GCoorWei_13[i][1];
		}
		break;
	case 14:
		for(i=0; i<hpt; i++) {
			crdgp[hpt + i] = GCoorWei_14[i][0];
			crdgp[hp1 - i] = -crdgp[hpt + i];
			weigp[hpt + i] = weigp[hp1 - i] = GCoorWei_14[i][1];
		}
		break;
	case 15:
		crdgp[hpt] = GCoorWei_15[0][0];
		weigp[hpt] = GCoorWei_15[0][1];
		for(i=1; i<=hpt; i++) {
			crdgp[hpt + i] = GCoorWei_15[i][0];
			crdgp[hpt - i] = -crdgp[hpt + i];
			weigp[hpt + i] = weigp[hpt - i] = GCoorWei_15[i][1];
		}
		break;
	case 16:
		for(i=0; i<hpt; i++) {
			crdgp[hpt + i] = GCoorWei_16[i][0];
			crdgp[hp1 - i] = -crdgp[hpt + i];
			weigp[hpt + i] = weigp[hp1 - i] = GCoorWei_16[i][1];
		}
		break;
	case 17:
		crdgp[hpt] = GCoorWei_17[0][0];
		weigp[hpt] = GCoorWei_17[0][1];
		for(i=1; i<=hpt; i++) {
			crdgp[hpt + i] = GCoorWei_17[i][0];
			crdgp[hpt - i] = -crdgp[hpt + i];
			weigp[hpt + i] = weigp[hpt - i] = GCoorWei_17[i][1];
		}
		break;
	case 18:
		for(i=0; i<hpt; i++) {
			crdgp[hpt + i] = GCoorWei_18[i][0];
			crdgp[hp1 - i] = -crdgp[hpt + i];
			weigp[hpt + i] = weigp[hp1 - i] = GCoorWei_18[i][1];
		}
		break;
	case 19:
		crdgp[hpt] = GCoorWei_19[0][0];
		weigp[hpt] = GCoorWei_19[0][1];
		for(i=1; i<=hpt; i++) {
			crdgp[hpt + i] = GCoorWei_19[i][0];
			crdgp[hpt - i] = -crdgp[hpt + i];
			weigp[hpt + i] = weigp[hpt - i] = GCoorWei_19[i][1];
		}
		break;
	case 20:
		for(i=0; i<hpt; i++) {
			crdgp[hpt + i] = GCoorWei_20[i][0];
			crdgp[hp1 - i] = -crdgp[hpt + i];
			weigp[hpt + i] = weigp[hp1 - i] = GCoorWei_20[i][1];
		}
		break;
	case 21:
		crdgp[hpt] = GCoorWei_21[0][0];
		weigp[hpt] = GCoorWei_21[0][1];
		for(i=1; i<=hpt; i++) {
			crdgp[hpt + i] = GCoorWei_21[i][0];
			crdgp[hpt - i] = -crdgp[hpt + i];
			weigp[hpt + i] = weigp[hpt - i] = GCoorWei_21[i][1];
		}
		break;
	case 22:
		for(i=0; i<hpt; i++) {
			crdgp[hpt + i] = GCoorWei_22[i][0];
			crdgp[hp1 - i] = -crdgp[hpt + i];
			weigp[hpt + i] = weigp[hp1 - i] = GCoorWei_22[i][1];
		}
		break;
	case 23:
		crdgp[hpt] = GCoorWei_23[0][0];
		weigp[hpt] = GCoorWei_23[0][1];
		for(i=1; i<=hpt; i++) {
			crdgp[hpt + i] = GCoorWei_23[i][0];
			crdgp[hpt - i] = -crdgp[hpt + i];
			weigp[hpt + i] = weigp[hpt - i] = GCoorWei_23[i][1];
		}
		break;
	case 24:
		for(i=0; i<hpt; i++) {
			crdgp[hpt + i] = GCoorWei_24[i][0];
			crdgp[hp1 - i] = -crdgp[hpt + i];
			weigp[hpt + i] = weigp[hp1 - i] = GCoorWei_24[i][1];
		}
		break;
	case 25:
		crdgp[hpt] = GCoorWei_25[0][0];
		weigp[hpt] = GCoorWei_25[0][1];
		for(i=1; i<=hpt; i++) {
			crdgp[hpt + i] = GCoorWei_25[i][0];
			crdgp[hpt - i] = -crdgp[hpt + i];
			weigp[hpt + i] = weigp[hpt - i] = GCoorWei_25[i][1];
		}
		break;
	case 26:
		for(i=0; i<hpt; i++) {
			crdgp[hpt + i] = GCoorWei_26[i][0];
			crdgp[hp1 - i] = -crdgp[hpt + i];
			weigp[hpt + i] = weigp[hp1 - i] = GCoorWei_26[i][1];
		}
		break;
	case 27:
		crdgp[hpt] = GCoorWei_27[0][0];
		weigp[hpt] = GCoorWei_27[0][1];
		for(i=1; i<=hpt; i++) {
			crdgp[hpt + i] = GCoorWei_27[i][0];
			crdgp[hpt - i] = -crdgp[hpt + i];
			weigp[hpt + i] = weigp[hpt - i] = GCoorWei_27[i][1];
		}
		break;
	case 28:
		for(i=0; i<hpt; i++) {
			crdgp[hpt + i] = GCoorWei_28[i][0];
			crdgp[hp1 - i] = -crdgp[hpt + i];
			weigp[hpt + i] = weigp[hp1 - i] = GCoorWei_28[i][1];
		}
		break;
	case 29:
		crdgp[hpt] = GCoorWei_29[0][0];
		weigp[hpt] = GCoorWei_29[0][1];
		for(i=1; i<=hpt; i++) {
			crdgp[hpt + i] = GCoorWei_29[i][0];
			crdgp[hpt - i] = -crdgp[hpt + i];
			weigp[hpt + i] = weigp[hpt - i] = GCoorWei_29[i][1];
		}
		break;
	case 30:
		for(i=0; i<hpt; i++) {
			crdgp[hpt + i] = GCoorWei_30[i][0];
			crdgp[hp1 - i] = -crdgp[hpt + i];
			weigp[hpt + i] = weigp[hp1 - i] = GCoorWei_30[i][1];
		}
		break;
	case 31:
		crdgp[hpt] = GCoorWei_31[0][0];
		weigp[hpt] = GCoorWei_31[0][1];
		for(i=1; i<=hpt; i++) {
			crdgp[hpt + i] = GCoorWei_31[i][0];
			crdgp[hpt - i] = -crdgp[hpt + i];
			weigp[hpt + i] = weigp[hpt - i] = GCoorWei_31[i][1];
		}
		break;
	case 32:
		for(i=0; i<hpt; i++) {
			crdgp[hpt + i] = GCoorWei_32[i][0];
			crdgp[hp1 - i] = -crdgp[hpt + i];
			weigp[hpt + i] = weigp[hp1 - i] = GCoorWei_32[i][1];
		}
		break;
	case 33:
		crdgp[hpt] = GCoorWei_33[0][0];
		weigp[hpt] = GCoorWei_33[0][1];
		for(i=1; i<=hpt; i++) {
			crdgp[hpt + i] = GCoorWei_33[i][0];
			crdgp[hpt - i] = -crdgp[hpt + i];
			weigp[hpt + i] = weigp[hpt - i] = GCoorWei_33[i][1];
		}
		break;
	case 34:
		for(i=0; i<hpt; i++) {
			crdgp[hpt + i] = GCoorWei_34[i][0];
			crdgp[hp1 - i] = -crdgp[hpt + i];
			weigp[hpt + i] = weigp[hp1 - i] = GCoorWei_34[i][1];
		}
		break;
	case 35:
		crdgp[hpt] = GCoorWei_35[0][0];
		weigp[hpt] = GCoorWei_35[0][1];
		for(i=1; i<=hpt; i++) {
			crdgp[hpt + i] = GCoorWei_35[i][0];
			crdgp[hpt - i] = -crdgp[hpt + i];
			weigp[hpt + i] = weigp[hpt - i] = GCoorWei_35[i][1];
		}
		break;
	case 36:
		for(i=0; i<hpt; i++) {
			crdgp[hpt + i] = GCoorWei_36[i][0];
			crdgp[hp1 - i] = -crdgp[hpt + i];
			weigp[hpt + i] = weigp[hp1 - i] = GCoorWei_36[i][1];
		}
		break;
	case 37:
		crdgp[hpt] = GCoorWei_37[0][0];
		weigp[hpt] = GCoorWei_37[0][1];
		for(i=1; i<=hpt; i++) {
			crdgp[hpt + i] = GCoorWei_37[i][0];
			crdgp[hpt - i] = -crdgp[hpt + i];
			weigp[hpt + i] = weigp[hpt - i] = GCoorWei_37[i][1];
		}
		break;
	case 38:
		for(i=0; i<hpt; i++) {
			crdgp[hpt + i] = GCoorWei_38[i][0];
			crdgp[hp1 - i] = -crdgp[hpt + i];
			weigp[hpt + i] = weigp[hp1 - i] = GCoorWei_38[i][1];
		}
		break;
	case 39:
		crdgp[hpt] = GCoorWei_39[0][0];
		weigp[hpt] = GCoorWei_39[0][1];
		for(i=1; i<=hpt; i++) {
			crdgp[hpt + i] = GCoorWei_39[i][0];
			crdgp[hpt - i] = -crdgp[hpt + i];
			weigp[hpt + i] = weigp[hpt - i] = GCoorWei_39[i][1];
		}
		break;
	case 40:
		for(i=0; i<hpt; i++) {
			crdgp[hpt + i] = GCoorWei_40[i][0];
			crdgp[hp1 - i] = -crdgp[hpt + i];
			weigp[hpt + i] = weigp[hp1 - i] = GCoorWei_40[i][1];
		}
		break;
	case 41:
		crdgp[hpt] = GCoorWei_41[0][0];
		weigp[hpt] = GCoorWei_41[0][1];
		for(i=1; i<=hpt; i++) {
			crdgp[hpt + i] = GCoorWei_41[i][0];
			crdgp[hpt - i] = -crdgp[hpt + i];
			weigp[hpt + i] = weigp[hpt - i] = GCoorWei_41[i][1];
		}
		break;
	case 42:
		for(i=0; i<hpt; i++) {
			crdgp[hpt + i] = GCoorWei_42[i][0];
			crdgp[hp1 - i] = -crdgp[hpt + i];
			weigp[hpt + i] = weigp[hp1 - i] = GCoorWei_42[i][1];
		}
		break;
	case 43:
		crdgp[hpt] = GCoorWei_43[0][0];
		weigp[hpt] = GCoorWei_43[0][1];
		for(i=1; i<=hpt; i++) {
			crdgp[hpt + i] = GCoorWei_43[i][0];
			crdgp[hpt - i] = -crdgp[hpt + i];
			weigp[hpt + i] = weigp[hpt - i] = GCoorWei_43[i][1];
		}
		break;
	case 44:
		for(i=0; i<hpt; i++) {
			crdgp[hpt + i] = GCoorWei_44[i][0];
			crdgp[hp1 - i] = -crdgp[hpt + i];
			weigp[hpt + i] = weigp[hp1 - i] = GCoorWei_44[i][1];
		}
		break;
	case 45:
		crdgp[hpt] = GCoorWei_45[0][0];
		weigp[hpt] = GCoorWei_45[0][1];
		for(i=1; i<=hpt; i++) {
			crdgp[hpt + i] = GCoorWei_45[i][0];
			crdgp[hpt - i] = -crdgp[hpt + i];
			weigp[hpt + i] = weigp[hpt - i] = GCoorWei_45[i][1];
		}
		break;
	case 46:
		for(i=0; i<hpt; i++) {
			crdgp[hpt + i] = GCoorWei_46[i][0];
			crdgp[hp1 - i] = -crdgp[hpt + i];
			weigp[hpt + i] = weigp[hp1 - i] = GCoorWei_46[i][1];
		}
		break;
	case 47:
		crdgp[hpt] = GCoorWei_47[0][0];
		weigp[hpt] = GCoorWei_47[0][1];
		for(i=1; i<=hpt; i++) {
			crdgp[hpt + i] = GCoorWei_47[i][0];
			crdgp[hpt - i] = -crdgp[hpt + i];
			weigp[hpt + i] = weigp[hpt - i] = GCoorWei_47[i][1];
		}
		break;
	case 48:
		for(i=0; i<hpt; i++) {
			crdgp[hpt + i] = GCoorWei_48[i][0];
			crdgp[hp1 - i] = -crdgp[hpt + i];
			weigp[hpt + i] = weigp[hp1 - i] = GCoorWei_48[i][1];
		}
		break;
	case 49:
		crdgp[hpt] = GCoorWei_49[0][0];
		weigp[hpt] = GCoorWei_49[0][1];
		for(i=1; i<=hpt; i++) {
			crdgp[hpt + i] = GCoorWei_49[i][0];
			crdgp[hpt - i] = -crdgp[hpt + i];
			weigp[hpt + i] = weigp[hpt - i] = GCoorWei_49[i][1];
		}
		break;
	case 50:
		for(i=0; i<hpt; i++) {
			crdgp[hpt + i] = GCoorWei_50[i][0];
			crdgp[hp1 - i] = -crdgp[hpt + i];
			weigp[hpt + i] = weigp[hp1 - i] = GCoorWei_50[i][1];
		}
		break;
	case 51:
		crdgp[hpt] = GCoorWei_51[0][0];
		weigp[hpt] = GCoorWei_51[0][1];
		for(i=1; i<=hpt; i++) {
			crdgp[hpt + i] = GCoorWei_51[i][0];
			crdgp[hpt - i] = -crdgp[hpt + i];
			weigp[hpt + i] = weigp[hpt - i] = GCoorWei_51[i][1];
		}
		break;
	case 52:
		for(i=0; i<hpt; i++) {
			crdgp[hpt + i] = GCoorWei_52[i][0];
			crdgp[hp1 - i] = -crdgp[hpt + i];
			weigp[hpt + i] = weigp[hp1 - i] = GCoorWei_52[i][1];
		}
		break;
	case 53:
		crdgp[hpt] = GCoorWei_53[0][0];
		weigp[hpt] = GCoorWei_53[0][1];
		for(i=1; i<=hpt; i++) {
			crdgp[hpt + i] = GCoorWei_53[i][0];
			crdgp[hpt - i] = -crdgp[hpt + i];
			weigp[hpt + i] = weigp[hpt - i] = GCoorWei_53[i][1];
		}
		break;
	case 54:
		for(i=0; i<hpt; i++) {
			crdgp[hpt + i] = GCoorWei_54[i][0];
			crdgp[hp1 - i] = -crdgp[hpt + i];
			weigp[hpt + i] = weigp[hp1 - i] = GCoorWei_54[i][1];
		}
		break;
	case 55:
		crdgp[hpt] = GCoorWei_55[0][0];
		weigp[hpt] = GCoorWei_55[0][1];
		for(i=1; i<=hpt; i++) {
			crdgp[hpt + i] = GCoorWei_55[i][0];
			crdgp[hpt - i] = -crdgp[hpt + i];
			weigp[hpt + i] = weigp[hpt - i] = GCoorWei_55[i][1];
		}
		break;
	case 56:
		for(i=0; i<hpt; i++) {
			crdgp[hpt + i] = GCoorWei_56[i][0];
			crdgp[hp1 - i] = -crdgp[hpt + i];
			weigp[hpt + i] = weigp[hp1 - i] = GCoorWei_56[i][1];
		}
		break;
	case 57:
		crdgp[hpt] = GCoorWei_57[0][0];
		weigp[hpt] = GCoorWei_57[0][1];
		for(i=1; i<=hpt; i++) {
			crdgp[hpt + i] = GCoorWei_57[i][0];
			crdgp[hpt - i] = -crdgp[hpt + i];
			weigp[hpt + i] = weigp[hpt - i] = GCoorWei_57[i][1];
		}
		break;
	case 58:
		for(i=0; i<hpt; i++) {
			crdgp[hpt + i] = GCoorWei_58[i][0];
			crdgp[hp1 - i] = -crdgp[hpt + i];
			weigp[hpt + i] = weigp[hp1 - i] = GCoorWei_58[i][1];
		}
		break;
	case 59:
		crdgp[hpt] = GCoorWei_59[0][0];
		weigp[hpt] = GCoorWei_59[0][1];
		for(i=1; i<=hpt; i++) {
			crdgp[hpt + i] = GCoorWei_59[i][0];
			crdgp[hpt - i] = -crdgp[hpt + i];
			weigp[hpt + i] = weigp[hpt - i] = GCoorWei_59[i][1];
		}
		break;
	case 60:
		for(i=0; i<hpt; i++) {
			crdgp[hpt + i] = GCoorWei_60[i][0];
			crdgp[hp1 - i] = -crdgp[hpt + i];
			weigp[hpt + i] = weigp[hp1 - i] = GCoorWei_60[i][1];
		}
		break;
	case 61:
		crdgp[hpt] = GCoorWei_61[0][0];
		weigp[hpt] = GCoorWei_61[0][1];
		for(i=1; i<=hpt; i++) {
			crdgp[hpt + i] = GCoorWei_61[i][0];
			crdgp[hpt - i] = -crdgp[hpt + i];
			weigp[hpt + i] = weigp[hpt - i] = GCoorWei_61[i][1];
		}
		break;
	case 62:
		for(i=0; i<hpt; i++) {
			crdgp[hpt + i] = GCoorWei_62[i][0];
			crdgp[hp1 - i] = -crdgp[hpt + i];
			weigp[hpt + i] = weigp[hp1 - i] = GCoorWei_62[i][1];
		}
		break;
	case 63:
		crdgp[hpt] = GCoorWei_63[0][0];
		weigp[hpt] = GCoorWei_63[0][1];
		for(i=1; i<=hpt; i++) {
			crdgp[hpt + i] = GCoorWei_63[i][0];
			crdgp[hpt - i] = -crdgp[hpt + i];
			weigp[hpt + i] = weigp[hpt - i] = GCoorWei_63[i][1];
		}
		break;
	case 64:
		for(i=0; i<hpt; i++) {
			crdgp[hpt + i] = GCoorWei_64[i][0];
			crdgp[hp1 - i] = -crdgp[hpt + i];
			weigp[hpt + i] = weigp[hp1 - i] = GCoorWei_64[i][1];
		}
		break;
	case 65:
		crdgp[hpt] = GCoorWei_65[0][0];
		weigp[hpt] = GCoorWei_65[0][1];
		for(i=1; i<=hpt; i++) {
			crdgp[hpt + i] = GCoorWei_65[i][0];
			crdgp[hpt - i] = -crdgp[hpt + i];
			weigp[hpt + i] = weigp[hpt - i] = GCoorWei_65[i][1];
		}
		break;
	case 66:
		for(i=0; i<hpt; i++) {
			crdgp[hpt + i] = GCoorWei_66[i][0];
			crdgp[hp1 - i] = -crdgp[hpt + i];
			weigp[hpt + i] = weigp[hp1 - i] = GCoorWei_66[i][1];
		}
		break;
	case 67:
		crdgp[hpt] = GCoorWei_67[0][0];
		weigp[hpt] = GCoorWei_67[0][1];
		for(i=1; i<=hpt; i++) {
			crdgp[hpt + i] = GCoorWei_67[i][0];
			crdgp[hpt - i] = -crdgp[hpt + i];
			weigp[hpt + i] = weigp[hpt - i] = GCoorWei_67[i][1];
		}
		break;
	case 68:
		for(i=0; i<hpt; i++) {
			crdgp[hpt + i] = GCoorWei_68[i][0];
			crdgp[hp1 - i] = -crdgp[hpt + i];
			weigp[hpt + i] = weigp[hp1 - i] = GCoorWei_68[i][1];
		}
		break;
	case 69:
		crdgp[hpt] = GCoorWei_69[0][0];
		weigp[hpt] = GCoorWei_69[0][1];
		for(i=1; i<=hpt; i++) {
			crdgp[hpt + i] = GCoorWei_69[i][0];
			crdgp[hpt - i] = -crdgp[hpt + i];
			weigp[hpt + i] = weigp[hpt - i] = GCoorWei_69[i][1];
		}
		break;
	case 70:
		for(i=0; i<hpt; i++) {
			crdgp[hpt + i] = GCoorWei_70[i][0];
			crdgp[hp1 - i] = -crdgp[hpt + i];
			weigp[hpt + i] = weigp[hp1 - i] = GCoorWei_70[i][1];
		}
		break;
	case 71:
		crdgp[hpt] = GCoorWei_71[0][0];
		weigp[hpt] = GCoorWei_71[0][1];
		for(i=1; i<=hpt; i++) {
			crdgp[hpt + i] = GCoorWei_71[i][0];
			crdgp[hpt - i] = -crdgp[hpt + i];
			weigp[hpt + i] = weigp[hpt - i] = GCoorWei_71[i][1];
		}
		break;
	case 72:
		for(i=0; i<hpt; i++) {
			crdgp[hpt + i] = GCoorWei_72[i][0];
			crdgp[hp1 - i] = -crdgp[hpt + i];
			weigp[hpt + i] = weigp[hp1 - i] = GCoorWei_72[i][1];
		}
		break;
	case 73:
		crdgp[hpt] = GCoorWei_73[0][0];
		weigp[hpt] = GCoorWei_73[0][1];
		for(i=1; i<=hpt; i++) {
			crdgp[hpt + i] = GCoorWei_73[i][0];
			crdgp[hpt - i] = -crdgp[hpt + i];
			weigp[hpt + i] = weigp[hpt - i] = GCoorWei_73[i][1];
		}
		break;
	case 74:
		for(i=0; i<hpt; i++) {
			crdgp[hpt + i] = GCoorWei_74[i][0];
			crdgp[hp1 - i] = -crdgp[hpt + i];
			weigp[hpt + i] = weigp[hp1 - i] = GCoorWei_74[i][1];
		}
		break;
	case 75:
		crdgp[hpt] = GCoorWei_75[0][0];
		weigp[hpt] = GCoorWei_75[0][1];
		for(i=1; i<=hpt; i++) {
			crdgp[hpt + i] = GCoorWei_75[i][0];
			crdgp[hpt - i] = -crdgp[hpt + i];
			weigp[hpt + i] = weigp[hpt - i] = GCoorWei_75[i][1];
		}
		break;
	case 76:
		for(i=0; i<hpt; i++) {
			crdgp[hpt + i] = GCoorWei_76[i][0];
			crdgp[hp1 - i] = -crdgp[hpt + i];
			weigp[hpt + i] = weigp[hp1 - i] = GCoorWei_76[i][1];
		}
		break;
	case 77:
		crdgp[hpt] = GCoorWei_77[0][0];
		weigp[hpt] = GCoorWei_77[0][1];
		for(i=1; i<=hpt; i++) {
			crdgp[hpt + i] = GCoorWei_77[i][0];
			crdgp[hpt - i] = -crdgp[hpt + i];
			weigp[hpt + i] = weigp[hpt - i] = GCoorWei_77[i][1];
		}
		break;
	case 78:
		for(i=0; i<hpt; i++) {
			crdgp[hpt + i] = GCoorWei_78[i][0];
			crdgp[hp1 - i] = -crdgp[hpt + i];
			weigp[hpt + i] = weigp[hp1 - i] = GCoorWei_78[i][1];
		}
		break;
	case 79:
		crdgp[hpt] = GCoorWei_79[0][0];
		weigp[hpt] = GCoorWei_79[0][1];
		for(i=1; i<=hpt; i++) {
			crdgp[hpt + i] = GCoorWei_79[i][0];
			crdgp[hpt - i] = -crdgp[hpt + i];
			weigp[hpt + i] = weigp[hpt - i] = GCoorWei_79[i][1];
		}
		break;
	case 80:
		for(i=0; i<hpt; i++) {
			crdgp[hpt + i] = GCoorWei_80[i][0];
			crdgp[hp1 - i] = -crdgp[hpt + i];
			weigp[hpt + i] = weigp[hp1 - i] = GCoorWei_80[i][1];
		}
		break;
	default:
		NC_Error_Exit_1(NCout, "Number of Gaussean points must be in the interval [1, 80]!",
			"Number of Gaussean points = ", npt);
	}
}

// CGS solver for the BEM equation system (Andreas Meister, p 168)
void NC_IterativeSolverCGS
(
	ofstream& NCout
)
{
	Vector<Complex> zX_j(numRowsOfCoefficientMatrix_), zR_0(numRowsOfCoefficientMatrix_), zR_j(numRowsOfCoefficientMatrix_), zP_j(numRowsOfCoefficientMatrix_), zV_j(numRowsOfCoefficientMatrix_),
		zQ_j(numRowsOfCoefficientMatrix_), zU_j(numRowsOfCoefficientMatrix_), zUQ_j(numRowsOfCoefficientMatrix_), zAUQ_j(numRowsOfCoefficientMatrix_);
	int i, j;
	double err_ori = 0.0, err_rel = 1.0, dwk1;
	Complex zalph, zbet, zrjr0, zrjr1;
	bool ifmodyprecond = false;

BeginCGS:
	// compute the arrays used by the preconditioners
	switch(methodPreconditioner_)
	{
	case 0: // ILU
		NC_IncompleteLUDecomposition(NCout);
		break;
	case 1: // row scanning
		NC_ComputeScalingVector(NCout);
		break;
	}


	// initialize the solution vector by using the random function
	for(j=0; j<numRowsOfCoefficientMatrix_; j++)
	{
	  //        zX_j[j].set(((double)(rand())/(double)(RAND_MAX) - 0.5)*2.0,
          //              ((double)(rand())/(double)(RAND_MAX) - 0.5)*2.0);
	  zX_j[j].set(0.0,0.0);
	}

	NC_MatrixVectorMultiplication(zX_j, zV_j);

	//see if the ILU-preconditioning matrix is almost singular
	if(!methodPreconditioner_)
	{
		for(i=0; i<numRowsOfCoefficientMatrix_; i++)
		{
			if(zV_j[i].norm() > 1.0e8)
			{
				if(methodFMM_)
				{
					ifmodyprecond = true;
					methodPreconditioner_ = 1;
					NC_Error_Warning_0(NCout,
					"ILU preconditioning failed, the row scanning method is used!");
					goto BeginCGS;
				}
				else
				{
					goto DirectCGS;
				}
			}
		}
	}

	for(i=0; i<numRowsOfCoefficientMatrix_; i++) zU_j[i] = zR_j[i] = zP_j[i] = zR_0[i] = zrhs[i] - zV_j[i];

	for(j=0; j<=niter_max_; j++)
	{
		NC_MatrixVectorMultiplication(zP_j, zV_j);

		if(j==0) {zrjr0 = zR_j*zR_0;} else {zrjr0 = zrjr1;}
		zalph = zrjr0/(zV_j*zR_0);

		for(i=0; i<numRowsOfCoefficientMatrix_; i++) zQ_j[i] = zU_j[i] - zalph*zV_j[i];
		zUQ_j = zU_j + zQ_j;

		NC_MatrixVectorMultiplication(zUQ_j, zAUQ_j);

		dwk1 = 0;
		for(i=0; i<numRowsOfCoefficientMatrix_; i++)
		{
		  zX_j[i] += zalph*zUQ_j[i];
		  zR_j[i] -= zalph*zAUQ_j[i];
		  
		  dwk1 += zR_j[i].qnorm();
		}

		if(j == 0)
		{
			err_ori = sqrt(dwk1);
			cout << "\nCGS: err_ori = " << err_ori << endl;
		}
		else
		{
            err_rel = sqrt(dwk1)/err_ori;
		}

		if(j > 0 && j/10*10 == j) {
		  cout << j << " " << err_rel << endl;
		  cout << j << " Abs Error:" << sqrt(dwk1) << endl;
		}

		if(err_rel < ErroIterSols || j == niter_max_)
		{
			if(j/10*10 != j) cout << j << " " << err_rel << "\n" << endl;
			for(i=0; i<numRowsOfCoefficientMatrix_; i++) zrhs[i] = zX_j[i];
			if(j == niter_max_) {
				NC_Error_Warning_0(NCout, "Warning: Maximum number of iterations is reached!");
			}
			break;
		}

		zrjr1 = zR_j*zR_0;
		zbet = zrjr1/zrjr0;
		for(i=0; i<numRowsOfCoefficientMatrix_; i++) zU_j[i] = zR_j[i] + zbet*zQ_j[i];
		for(i=0; i<numRowsOfCoefficientMatrix_; i++) zP_j[i] = zU_j[i] + zbet*(zQ_j[i] + zbet*zP_j[i]);
	} // end of loop j

	NCout << "\nCGS solver: number of iterations = "
		<< j << ", relative error = " << err_rel << endl;

	if(j >= niter_max_ && methodFMM_ == 0)
	{
DirectCGS:
		NC_Error_Warning_0(NCout, "Iteration method CGS failed, direct method is used!");

		// Gauss elimination method
		Tfactor_usy(zcoefl, numRowsOfCoefficientMatrix_);
		Tfbelim(zcoefl, zrhs, numRowsOfCoefficientMatrix_);

		return;
	}

	// if a preconditioner is used, modify the result
	switch(methodPreconditioner_)
	{
	case 0: // ILU
		NC_IncompleteLUForBack();
		NC_DeleteIncompleteLUMatrices();
		break;
	case 1: // row scanning
        for(i=0; i<numRowsOfCoefficientMatrix_; i++) zrhs[i] *= dscaling[i];
		delete [] dscaling;
		break;
	}

	if(ifmodyprecond) methodPreconditioner_ = 0;

}

// compute the product of the coefficient matrix and a vector
void NC_MatrixVectorMultiplication
(
	Vector<Complex>& zmultvect,		//I: multiplier vector
	Vector<Complex>& zresuvect      //O: result vector
)
{
	int i, j, k, ilv;
	Vector<Complex> ztmp(numRowsOfCoefficientMatrix_);

	// prconditioning
	switch(methodPreconditioner_)
	{
	case 0: // ILU
		NC_IncompleteLUForBack(NCout, zmultvect, ztmp);
		break;
	case 1: // scaling
		for(i=0; i<numRowsOfCoefficientMatrix_; i++) ztmp[i] = zmultvect[i]*dscaling[i];
		break;
	case 2: // no preconditioning
		ztmp = zmultvect;
		break;
	}

	// compute the product
	switch(methodFMM_)
	{
	case 0: // TBEM
		k = 0;
		for(i=0; i<numRowsOfCoefficientMatrix_; i++)
		{
			zresuvect[i].set(0.0, 0.0);
			for(j=0; j<numRowsOfCoefficientMatrix_; j++) zresuvect[i] += zcoefl[k + j]*ztmp[j];
			k += numRowsOfCoefficientMatrix_;
		}
		break;

	case 1: // SLFMBEM
		{
			Vector<Complex> zvc_t(numIntegrationPointsUnitSphere_*numOriginalReflectedClusters_), zvc_d(dmtxlev[0].nRowsD);

			for(i=0; i<numIntegrationPointsUnitSphere_*numOriginalReflectedClusters_; i++)
			{
				zvc_t[i].set(0.0, 0.0);
				for(j=irowtmtx[i]; j<irowtmtx[i + 1]; j++)
					zvc_t[i] += ztmtx[j]*ztmp[jcoltmtx[j]];
			}

			for(i=0; i<dmtxlev[0].nRowsD; i++)
			{
				zvc_d[i].set(0.0, 0.0);
				for(j=dmtxlev[0].irowDmxLv[i]; j<dmtxlev[0].irowDmxLv[i + 1]; j++)
					zvc_d[i] += dmtxlev[0].zDmxLv[j]*zvc_t[dmtxlev[0].jcolDmxLv[j]];
			}

			for(i=0; i<numRowsOfCoefficientMatrix_; i++)
			{
				zresuvect[i].set(0.0, 0.0);

				// contribution of the near field matrix
				for(j=irownea[i]; j<irownea[i + 1]; j++)
					zresuvect[i] += zcoefl[j]*ztmp[jcolnea[j]];

				// contribution of the far field matrix
				for(j=irowsmtx[i]; j<irowsmtx[i + 1]; j++)
					zresuvect[i] += zsmtx[j]*zvc_d[jcolsmtx[j]];
			}
		}
		break;
	case 3: // DMLFMBEM
		// contribution of the near fields
		for(i=0; i<numRowsOfCoefficientMatrix_; i++)
		{
			zresuvect[i].set(0.0, 0.0);
			for(j=irownea[i]; j<irownea[i + 1]; j++)
				zresuvect[i] += zcoefl[j]*ztmp[jcolnea[j]];
		}

		// contributions of the far fields
		for(ilv=0; ilv<numClusterLevels_; ilv++)
		{
			// compute the T-vector t = T*vect
			for(i=0; i<clulevarry[ilv].nPoinSpheLv*clulevarry[ilv].nClustSLv; i++)
			{
				clulevarry[ilv].zwkT[i].set(0.0, 0.0);
				for(j=tmtxlev[ilv].irowTmxLv[i]; j<tmtxlev[ilv].irowTmxLv[i + 1]; j++)
					clulevarry[ilv].zwkT[i] +=
					tmtxlev[ilv].zTmxLv[j]*ztmp[tmtxlev[ilv].jcolTmxLv[j]];
			}

			// compute the S-vectors s = D*t
			for(i=0; i<dmtxlev[ilv].nRowsD; i++)
			{
				clulevarry[ilv].zwkS[i].set(0.0, 0.0);
				for(j=dmtxlev[ilv].irowDmxLv[i]; j<dmtxlev[ilv].irowDmxLv[i+1]; j++)
					clulevarry[ilv].zwkS[i] +=
					dmtxlev[ilv].zDmxLv[j]*clulevarry[ilv].zwkT[dmtxlev[ilv].jcolDmxLv[j]];
			}

			// contributions to the result vector: resu += S*s
			for(i=0; i<numRowsOfCoefficientMatrix_; i++)
			{
				for(j=smtxlev[ilv].irowSmxLv[i]; j<smtxlev[ilv].irowSmxLv[i + 1]; j++)
					zresuvect[i] +=
					smtxlev[ilv].zSmxLv[j]*clulevarry[ilv].zwkS[smtxlev[ilv].jcolSmxLv[j]];
			}
		} // end of loop ILV
		break;
	}
}

// compute the scaling vector
void NC_ComputeScalingVector
(
	ofstream& NCout
)
{
	int i, j, k;
	double dwk;

	dscaling = new double[numRowsOfCoefficientMatrix_];

	switch(methodFMM_)
	{
	case 0: // TBEM
		k = 0;
		for(i=0; i<numRowsOfCoefficientMatrix_; i++)
		{
			dwk = 0.0;
			for(j=0; j<numRowsOfCoefficientMatrix_; j++) dwk += zcoefl[k + j].qnorm();
			dscaling[i] = sqrt((double)(numRowsOfCoefficientMatrix_)/dwk);
			k += numRowsOfCoefficientMatrix_;
		}
		break;
	case 1: // SLFMBEM
	case 3: // DMLFMBEM
		for(i=0; i<numRowsOfCoefficientMatrix_; i++)
		{
			dwk = 0.0;
			for(j=irownea[i]; j<irownea[i+1]; j++) dwk += zcoefl[j].qnorm();
			dscaling[i] = sqrt((double)(irownea[i+1] - irownea[i])/dwk);
		}
		break;
	}
}

// incomplete decomposition of the coefficient matrix by using three bool vectors to respresent the sparse muster (smaller storage, quick) (see Andreas Meister, p 197)
void NC_IncompleteLUDecomposition
(
	ofstream& NCout
)
{
	int i, j, k, kl, ku, m, j1 = 0, ml, mu, n_trues = 0, ilv;
	int nunzL, nunzU;
	double threshfac; // threshold factor
	double dwk, dw1;
	Vector<bool> Mirow(numRowsOfCoefficientMatrix_), Micol(numRowsOfCoefficientMatrix_), Mkvct(numRowsOfCoefficientMatrix_);

	if(methodFMM_ == 0) { // TRBEM
		if(scanningDegreeLU_ == 0) {
			threshfac = 1.2;
		} else if(scanningDegreeLU_ == 1) {
			threshfac = 1.0;
		} else if(scanningDegreeLU_ == 2) {
			threshfac = 0.8;
		} else {
			threshfac = 0.6;
		}
	} else if(methodFMM_ == 1) { // SLFMBEM
		if(scanningDegreeLU_ == 0) {
			//threshfac = 1.1;
			threshfac = 0.9;
		} else if(scanningDegreeLU_ == 1) {
			//threshfac = 0.5;
			threshfac = 0.35;
		} else if(scanningDegreeLU_ == 2) {
			//threshfac = 0.1;
			threshfac = 0.07;
		} else {
			threshfac = 0.01;
		}
	} else { // MLFMBEM
		if(scanningDegreeLU_ == 0) {
			//threshfac = 0.8;
			threshfac = 0.65;
		} else if(scanningDegreeLU_ == 1) {
			//threshfac = 0.3;
			threshfac = 0.15;
		} else if(scanningDegreeLU_ == 2) {
			//threshfac = 0.1;
			threshfac = 0.05;
		} else {
			threshfac = 0.005;
		}
	}

	switch(methodFMM_)
	{
	case 0: // TBEM
		// scaling the coefficient matrix and the right hand side vector, compute number of "TRUES"
		k = 0;
		for(i=0; i<numRowsOfCoefficientMatrix_; i++)
		{
			dw1 = 0.0;
			for(j=0; j<numRowsOfCoefficientMatrix_; j++) dw1 += zcoefl[k + j].qnorm();
			dwk = sqrt((double)(numRowsOfCoefficientMatrix_)/dw1);
			for(j=0; j<numRowsOfCoefficientMatrix_; j++)
			{
				zcoefl[k + j] *= dwk;
				if(zcoefl[k + j].norm() > threshfac || i == j) n_trues++;
			}
			zrhs[i] *= dwk;
			k += numRowsOfCoefficientMatrix_;
		} // end of loop I
		break;
	case 1: // SLFMBEM
		// scaling the near fied matrix, the S-matrix and the right hand side vector, compute number of "TRUES"
		for(i=0; i<numRowsOfCoefficientMatrix_; i++)
		{
			dw1 = 0.0;
			for(j=irownea[i]; j<irownea[i + 1]; j++) dw1 += zcoefl[j].qnorm();
			dwk = sqrt((double)(irownea[i + 1] - irownea[i])/dw1);
			for(j=irownea[i]; j<irownea[i + 1]; j++)
			{
				zcoefl[j] *= dwk;
				if(zcoefl[j].norm() > threshfac || jcolnea[j] == i) n_trues++;
			}
			for(j=irowsmtx[i]; j<irowsmtx[i + 1]; j++) zsmtx[j] *= dwk;
			zrhs[i] *= dwk;
		}
		break;
	case 3: // DMLFMBEM
		// scaling the near fied matrix, the S-matrices and the right hand side vector, compute number of "TRUES"
		for(i=0; i<numRowsOfCoefficientMatrix_; i++)
		{
			dw1 = 0.0;
			for(j=irownea[i]; j<irownea[i + 1]; j++) dw1 += zcoefl[j].qnorm();
			dwk = sqrt((double)(irownea[i + 1] - irownea[i])/dw1);
			for(j=irownea[i]; j<irownea[i + 1]; j++)
			{
				zcoefl[j] *= dwk;
				if(zcoefl[j].norm() > threshfac || jcolnea[j] == i) n_trues++;
			}

			for(ilv=0; ilv<numClusterLevels_; ilv++)
			{
				for(j=smtxlev[ilv].irowSmxLv[i]; j<smtxlev[ilv].irowSmxLv[i + 1]; j++)
					smtxlev[ilv].zSmxLv[j] *= dwk;
			}

			zrhs[i] *= dwk;
		}
		break;
	}

	Vector<int> jcol_tru(n_trues), nrow_tru(numRowsOfCoefficientMatrix_ + 1), irow_tru(n_trues),
		ncol_tru(numRowsOfCoefficientMatrix_ + 1);

	// compute the vectors to store the column numbers and row numbers of the "TRUES"
	switch(methodFMM_)
	{
	case 0: // TBEM
		kl = k = 0;
		for(i=0; i<numRowsOfCoefficientMatrix_; i++) // loop over rows
		{
			nrow_tru[i] = kl;
			for(j=0; j<numRowsOfCoefficientMatrix_; j++) // loop over columns
			{
				if(zcoefl[k + j].norm() > threshfac || i == j) jcol_tru[kl++] = j;
			}
			k += numRowsOfCoefficientMatrix_;
		}
		nrow_tru[numRowsOfCoefficientMatrix_] = kl;
		break;
	case 1: // SLFMBEM
	case 3: // DMLFMBEM
		kl = 0;
		for(i=0; i<numRowsOfCoefficientMatrix_; i++) // loop over rows
		{
			nrow_tru[i] = kl;
			for(j=irownea[i]; j<irownea[i + 1]; j++) // loop over nonzeros of the row
			{
				if(zcoefl[j].norm() > threshfac || jcolnea[j] == i)
					jcol_tru[kl++] = jcolnea[j];
			}
		}
		nrow_tru[numRowsOfCoefficientMatrix_] = kl;
		break;
	}

	// compute the arrays NCOL_TRU and IROW_TRUE
	kl = 0;
	for(j=0; j<numRowsOfCoefficientMatrix_; j++) // loop over columns
	{
		ncol_tru[j] = kl;
		for(i=0; i<numRowsOfCoefficientMatrix_; i++) // loop over rows
		{
			for(k=nrow_tru[i]; k<nrow_tru[i+1]; k++)
			{
				if(jcol_tru[k] == j)
				{
					irow_tru[kl++] = i;
					break;
				}
			}
		}
	}
	ncol_tru[numRowsOfCoefficientMatrix_] = kl;

	// generate the L-matrix (stored by rows) and the U-matrix (stored by columns)
	nunzL = nunzU = 0;
	for(i=0; i<numRowsOfCoefficientMatrix_; i++)
	{
		for(j=0; j<=i; j++) Mirow[j] = Micol[j] = false;
		for(j=nrow_tru[i]; j<nrow_tru[i+1]; j++) Mirow[jcol_tru[j]] = true;
		for(j=ncol_tru[i]; j<ncol_tru[i+1]; j++) Micol[irow_tru[j]] = true;

		for(j=0; j<=i; j++) if(Mirow[j]) nunzL++;
		for(j=0; j<i; j++) if(Micol[j]) nunzU++;
	}
	nL_colnu = new int[nunzL];
	nL_firnu = new int[numRowsOfCoefficientMatrix_ + 1];
	nU_rownu = new int[nunzU];
	nU_firnu = new int[numRowsOfCoefficientMatrix_ + 1];
	zL_incom = new Complex[nunzL];
	zU_incom = new Complex[nunzU];

	kl = ku = 0;
	for(i=0; i<numRowsOfCoefficientMatrix_; i++)
	{
		for(j=0; j<=i; j++) Mirow[j] = Micol[j] = false;
		for(j=nrow_tru[i]; j<nrow_tru[i+1]; j++) Mirow[jcol_tru[j]] = true;
		for(j=ncol_tru[i]; j<ncol_tru[i+1]; j++) Micol[irow_tru[j]] = true;

		nL_firnu[i] = kl;
		for(j=0; j<=i; j++) if(Mirow[j]) nL_colnu[kl++] = j;

		nU_firnu[i] = ku;
		for(j=0; j<i; j++) if(Micol[j]) nU_rownu[ku++] = j;
	}
	nL_firnu[numRowsOfCoefficientMatrix_] = kl;
	nU_firnu[numRowsOfCoefficientMatrix_] = ku;

	// compute the L- and U-matrices
	switch(methodFMM_)
	{
	case 0: // TBEM
		k = 0;
		for(i=0; i<numRowsOfCoefficientMatrix_; i++)
		{
			for(j=nL_firnu[i]; j<nL_firnu[i + 1]; j++) zL_incom[j] = zcoefl[k + nL_colnu[j]];

			for(j=nU_firnu[i]; j<nU_firnu[i + 1]; j++)
				zU_incom[j] = zcoefl[nU_rownu[j]*numRowsOfCoefficientMatrix_ + i];
			k += numRowsOfCoefficientMatrix_;
		}
		break;
	case 1: // SLFMBEM
	case 3: // DMLFMBEM
		for(i=0; i<numRowsOfCoefficientMatrix_; i++)
		{
			for(j=nL_firnu[i]; j<nL_firnu[i + 1]; j++)
			{
				for(k=irownea[i]; k<irownea[i + 1]; k++) if(jcolnea[k] == nL_colnu[j])
				{
					kl = k;
					break;
				}
				zL_incom[j] = zcoefl[kl];
			}

			for(j=nU_firnu[i]; j<nU_firnu[i + 1]; j++)
			{
				for(k=irownea[nU_rownu[j]]; k<irownea[nU_rownu[j] + 1]; k++)
				if(jcolnea[k] == i)
				{
					ku = k;
					break;
				}
				zU_incom[j] = zcoefl[ku];
			}
		} // end of loop I
		break;
	}

	// incomplete LU decomposition
	for(i=0; i<numRowsOfCoefficientMatrix_; i++)
	{
		for(j=0; j<numRowsOfCoefficientMatrix_; j++) Mirow[j] = Micol[j] = false;
		for(j=nrow_tru[i]; j<nrow_tru[i+1]; j++) Mirow[jcol_tru[j]] = true;
		for(j=ncol_tru[i]; j<ncol_tru[i+1]; j++) Micol[irow_tru[j]] = true;

		for(k=i; k<numRowsOfCoefficientMatrix_; k++)
		{
			if(!Micol[k]) continue;

			for(j=nL_firnu[k]; j<nL_firnu[k + 1]; j++)
			{
				if(nL_colnu[j] == i)
				{
					j1 = j;
					break;
				}
			}

			for(j=0; j<numRowsOfCoefficientMatrix_; j++) Mkvct[j] = false;
			for(j=nrow_tru[k]; j<nrow_tru[k+1]; j++) Mkvct[jcol_tru[j]] = true;

			ml = mu = 0;
			for(m=0; m<i; m++)
			{
				if(Mkvct[m] && Micol[m])
				{
					zL_incom[j1] -= zL_incom[nL_firnu[k] + ml]*zU_incom[nU_firnu[i] + mu];
				}
				if(Mkvct[m]) ml++;
				if(Micol[m]) mu++;
			}
		}

		for(k=i+1; k<numRowsOfCoefficientMatrix_; k++)
		{
			if(!Mirow[k]) continue;

			for(j=nU_firnu[k]; j<nU_firnu[k + 1]; j++)
			{
				if(nU_rownu[j] == i)
				{
					j1 = j;
					break;
				}
			}

			for(j=0; j<numRowsOfCoefficientMatrix_; j++) Mkvct[j] = false;
			for(j=ncol_tru[k]; j<ncol_tru[k+1]; j++) Mkvct[irow_tru[j]] = true;

			ml = mu = 0;
			for(m=0; m<i; m++)
			{
				if(Mirow[m] && Mkvct[m])
				{
					zU_incom[j1] -= zL_incom[nL_firnu[i] + ml]*zU_incom[nU_firnu[k] + mu];
				}
				if(Mirow[m]) ml++;
				if(Mkvct[m]) mu++;
			}
			zU_incom[j1] /= zL_incom[nL_firnu[i + 1] - 1];
		}
	} // end of loop I
}

// destroy the incomplete L- and U-matrices
void NC_DeleteIncompleteLUMatrices()
{
	delete [] nL_colnu;
	delete [] nL_firnu;
	delete [] nU_rownu;
	delete [] nU_firnu;
	delete [] zL_incom;
	delete [] zU_incom;
}

// multiplication of a vector with the inverse of the incomplete LU-matrices
void NC_IncompleteLUForBack
(
	ofstream& NCout,
	Vector<Complex>& zmultvect,		//I: multiplier vector
	Vector<Complex>& zprodvect      //O: product vector
)
{
	int i, j, k;

	zprodvect = zmultvect;

	// forward elemination by using the L-matrix (stored by rows)
	zprodvect[0] /= zL_incom[0];
	for(i=1; i<numRowsOfCoefficientMatrix_; i++)
	{
		k = nL_firnu[i + 1] - 1;
		for(j=nL_firnu[i]; j<k; j++) zprodvect[i] -= zL_incom[j]*zprodvect[nL_colnu[j]];
		zprodvect[i] /= zL_incom[k];
	}

	// Backward substitution by using the U-matrix (stored by columns)
	for(i=numRowsOfCoefficientMatrix_-1; i>0; i--)
	{
		for(j=nU_firnu[i]; j<nU_firnu[i+1]; j++)
			zprodvect[nU_rownu[j]] -= zU_incom[j]*zprodvect[i];
	}
}
void NC_IncompleteLUForBack
(
)
{
	int i, j, k;

	// forward elemination by using the L-matrix (stored by rows)
	zrhs[0] /= zL_incom[0];
	for(i=1; i<numRowsOfCoefficientMatrix_; i++)
	{
		k = nL_firnu[i + 1] - 1;
		for(j=nL_firnu[i]; j<k; j++) zrhs[i] -= zL_incom[j]*zrhs[nL_colnu[j]];
		zrhs[i] /= zL_incom[k];
	}

	// Backward substitution by using the U-matrix (stored by columns)
	for(i=numRowsOfCoefficientMatrix_-1; i>0; i--)
	{
		for(j=nU_firnu[i]; j<nU_firnu[i+1]; j++)
			zrhs[nU_rownu[j]] -= zU_incom[j]*zrhs[i];
	}
}
