/*===========================================================================*\
 *                                                                            *
 *  File name:      NC_3dFunctions.cpp                                        *
 *  Description:    public functions for 3D problems                          *
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



// reflect an element used by the Fast Multipol BEM (FMBEM)
void NC_ReflectElementFMM
(
    int& rfa_cl,             //  I: factor = 1 or -1
	bool& mir_cl,            //  I: if is a mirror image of the original element
	Vector<bool>& rfdi_cl,   //  I: if reflected in a given coordinate direction
	const int& inodj,        //  I: number of nodes of the element
	const int& ibvj_0,       //  I: = 0: velocity boundary condition
			                 //       1: pressure boundary condition
							 //       2: transfer boundary condition
	Vector<double>& centerj, // IO: center of the element
	Vector<double>& norvecj, // IO: normal vector at the center and nodes
	Vector<Complex>& zbvj_0, // IO: boundary condition for the middle face
	Matrix<double>& crdelj   // IO: nodal coordinates
)
{
	int i, nrdir, nod1, nod2, ntime;
	double xi;
    Complex zbl;

	// update the velocity boundary condition
	if(ibvj_0 <= 1 && rfa_cl == -1) for(i=0; i<inodj; i++) zbvj_0[i].nega();

	// update the center, the normal vectors and the nodal coordinates of the element
	for(nrdir=0; nrdir<NDIM; nrdir++)
	{
		if(!rfdi_cl[nrdir]) continue;
		centerj[nrdir] = 2.0*listSymmetryPlanesCoordinates[nrdir] - centerj[nrdir];
		norvecj[nrdir] = -norvecj[nrdir];
		for(i=0; i<inodj; i++)
		{
			norvecj[NDIM*(i + 1) + nrdir] = -norvecj[NDIM*(i + 1) + nrdir];
			crdelj(i, nrdir) = 2.0*listSymmetryPlanesCoordinates[nrdir] - crdelj(i, nrdir);
		}
	}

	// update the sequence of the nodes in the local coordinate system of the element
	if(mir_cl)
	{
		nod1 = 0;
		nod2 = 1;
		ntime = 0;
IterBegin_rfFMB:
		ntime++;
		for(i=0; i<NDIM; i++) 
		{
			xi = crdelj(nod1, i);
			crdelj(nod1, i) = crdelj(nod2, i);
			crdelj(nod2, i) = xi;
			xi = norvecj[NDIM*(nod1 + 1) + i];
			norvecj[NDIM*(nod1 + 1) + i] = norvecj[NDIM*(nod2 + 1) + i];
			norvecj[NDIM*(nod2 + 1) + i] = xi;
		}
		zbl = zbvj_0[nod1];
		zbvj_0[nod1] = zbvj_0[nod2];
		zbvj_0[nod2] = zbl;
		if(inodj == NETYP4 && ntime == 1) 
		{
			nod1 = 2;
			nod2 = 3;
			goto IterBegin_rfFMB;
		}
	}
}

// compute the global coordinates of a point in an element
void NC_ComputeGlobalCoordinates
(
    Vector<double>& crdpi, Vector<double>& shfunx,
    Matrix<double>& crdel, const double& s, const double& t,
    const int& inode)
{
	int i, j;
	double s1, s2, t1, t2;

	switch(inode) {
	case 3:
		shfunx[0] = s;
		shfunx[1] = t;
		shfunx[2] = 1.0 - s - t;
		break;
	case 4:
		s1 = 0.25*(s + 1.0);
		s2 = 0.25*(s - 1.0);
		t1 = t + 1.0;
		t2 = t - 1.0;
		shfunx[0] = s1*t1;
		shfunx[1] = -s2*t1;
		shfunx[2] = s2*t2;
		shfunx[3] = -s1*t2;
		break;
	default:
		NC_Error_Exit_0(NCout, "number of nodes of an element must be 3 or 4!");
	}

	for(i=0; i<NDIM; i++)
	{
		crdpi[i] = 0.0;
		for(j=0; j<inode; j++) crdpi[i] += shfunx[j]*crdel(j, i);
	}
}

// compute the unit normal vector to an element at a given poit, return: Jacobean at the point
double NC_ComputeUnitVectorOnElement
(
    Vector<double>& elnorv, Vector<double>& shnds,
    Vector<double>& shndt, Vector<double>& dxds,
    Vector<double>& dxdt, Matrix<double>& crdel,
    const double& s, const double& t, const int& inode)
{
	int i, j;
	double length;

	switch(inode) {
	case 3:
		shnds[0] = shndt[1] = 1.0;
		shnds[1] = shndt[0] = 0.0;
		shnds[2] = shndt[2] = -1.0;
		break;
	case 4:
		shnds[0] = 0.25*(t + 1.0);
		shnds[2] = 0.25*(t - 1.0);
		shnds[1] = -shnds[0];
		shnds[3] = -shnds[2];

		shndt[0] = 0.25*(s + 1.0);
		shndt[2] = 0.25*(s - 1.0);
		shndt[1] = -shndt[2];
		shndt[3] = -shndt[0];
		break;
	default:
		NC_Error_Exit_0(NCout, "number of nodes of an element must be 3 or 4!");
	}

	// derivatives of the global coordinates with respect to the intrinsic coordinates
	for(i=0; i<NDIM; i++) 
	{
		dxds[i] = dxdt[i] = 0.0;
		for(j=0; j<inode; j++)
		{
			dxds[i] += shnds[j]*crdel(j, i);
			dxdt[i] += shndt[j]*crdel(j, i);
		}
	}

	elnorv.mul_3(dxds, dxdt);
	Vcnorm_dim3_(elnorv, length);

	return length;
}

// reflect a point or a vector in a given direction, return: factor of the reflection
double NC_ReflectPointVector
(
	Vector<double>& po_ve, // IO: a vector of dimension NDIM
	const int& kref,       //  I: number of the current reflection
	double& reffac,        // IO: factor of the reflection
	const int &ifpoi       //  I: = 0: a vector is reflected
				           //     = 1: a point is reflected
)
{
	int i, nrdir=0, i0=0;
	int nrefdi[] = {0, 0, 1, 0, 2, 0, 1, 0};

	if(kref == 0)
	{
		return(reffac);
	}
	else if(kref < 0 || kref > 7)
		NC_Error_Exit_0("kref must be equal to 0, 1, ... 7!");

	// direction of the current reflection
	switch(numReflectionsOfElements_)
	{
	case 8:  // there are 3 symmetric planes
		nrdir = nrefdi[kref];
		break;

	case 4:  // there are 2 symmetric planes
		for(i=0; i<NDIM; i++) if(listSymmetryPlanes[i] == 0)
		{
			i0 = i + 1;
			break;
		}
		nrdir = i0 + nrefdi[kref];
		if(nrdir >= NDIM) nrdir -= NDIM;
		break;
	case 2:  // there is 1 symmetric plane
		for(i=0; i<NDIM; i++) if(listSymmetryPlanes[i] != 0)
		{
			nrdir = i;
			break;
		}
		break;
	case 1:  // there is not symmetric plane
		return(reffac);
		break;
	default:
			NC_Error_Exit_0(NCout, "numReflectionsOfElements_ must be equal to 8, 4, 2 or 1!");
	}

	// update the point coordinate (or vector direction) and the reflection factor
	if(ifpoi) // for a point
	{
		po_ve[nrdir] = 2.0*listSymmetryPlanesCoordinates[nrdir] - po_ve[nrdir];
	}
	else      // for a vector
	{
		po_ve[nrdir] = - po_ve[nrdir];
	}
	if(listSymmetryPlanes[nrdir] < 0) reffac = - reffac;

	return reffac;
}

// compute data of an element
void NC_ComputeElementData
(
	const int& jel,          // I: number of the element
	const int& icase,        // I: switch for different applications
			                 //    = 0: ibvj_0 = kind of boundary condition
							 //         admj   = admittance boundary condition
							 //         zbvj_0 = value of boundary condition
							 //    = 1: ibvj_0 = flag for nonzero velocity
							 //         admj   = value of potentials
							 //         zbvj_0 = value of velocity or 
							 //                  velocity difference
	int& inodj,              // O: number of nodes
	int& ibvj_0,             // O: kind of boundary condition for the
			                 //    under face/if velocity is not zero
	int& ifadmij,            // O: flag for admittance boundary condition
	Vector<int>& jdofaddre,  // O: address of the unknown DOFs of the element
	Vector<double>& centerj, // O: coordinate of the center
	Vector<double>& norvecj, // O: normal vector at the center
	Vector<Complex>& admj,   // O: two admittances/potential
	Vector<Complex>& zbvj_0, // O: boundary conditions of the middle face/
			                 //    velocity or difference of velocities
	Matrix<double>& crdelj   // O: nodal coordinates
)
{
	int i, j;

	inodj = listNumberNodesPerElement[jel];

	// boundary conditions
	switch(ibval[jel]) 
	{
	case 0:         // velocity prescribed
		ibvj_0 = 0;  
		ifadmij = 0;
		break;
	case 1:         // pressure prescribed
		ibvj_0 = 1;
		ifadmij = 0;
		break;
	case 2:         // velocity and surface admittance prescribed
		ibvj_0 = 0;  
		ifadmij = 1;
		break;
	case 3:         // transfer admittance prescribed
		ibvj_0 = 2;
		ifadmij = 0;
		break;
	case 4:         // transfer and surface admittances prescribed
		ibvj_0 = 2;
		ifadmij = 1;
		break;
	}

	// addresse of the unknown DOFs
	jdofaddre[0] = jelist[jel][0];
	jdofaddre[1] = jelist[jel][1];

	// center point and normal vectors
	for(i=0; i<NDIM; i++) centerj[i] = centel[jel][i];
	for(j=0; j<inodj + 1; j++) for(i=0; i<NDIM; i++) 
		norvecj[j*NDIM + i] = elenor[jel][j*NDIM + i];

	// boundary condition for the under face
	for(i=0; i<inodj; i++) zbvj_0[i] = zbval0[jel][i];

	// coordinates of the nodes
	for(i=0; i<inodj; i++) for(j=0; j<NDIM; j++) 
		crdelj(i, j) = nodesCoordinates[elementsConnectivity[jel][i]][j];

    if(ifadmij)
    {
        admj[0] = zbval1[jel];
    }
    else
    {
        admj[0] = 0.0;
    }

	if(icase == 1) 
	{
		Complex ztmp;
        if(ibvj_0 == 1) // pressure boundary condition is prescribed
        {
            for(i=0; i<inodj; i++) zbvj_0[i] = zrhs[jdofaddre[0]];
        }
        else // velocity boundary condition is prescribed
        {
            if(ifadmij)
            {
                ztmp = (admj[0]*zrhs[jdofaddre[0]])*Tao_;
                for(i=0; i<inodj; i++) zbvj_0[i] -= ztmp;
                ibvj_0 = 1;
            }
            else
            {
                ibvj_0 = 0;
                for(i=0; i<inodj; i++)
                {
                    if(zbvj_0[i].norm() > EPSY)
                    {
                        ibvj_0 = 1;
                        break;
                    }
                }
            }
            admj[0] = zrhs[jdofaddre[0]];
        }

	} // end of IF(ICASE == 1)
}

// reflect an element in a given direction
void NC_ReflectElement
(
	const int& kref,         //  I: number of the current reflection
	const int& inodj,        //  I: number of nodes of the element
	const int& ibvj_0,       //  I: = 0: velocity boundary condition
			                 //       1: pressure boundary condition
							 //       2: transfer boundary condition
	int& ireffac,            // IO: factor of the reflection (= 1 or -1)
	Vector<double>& centerj, // IO: center of the element
	Vector<double>& norvecj, // IO: normal vector at the center and nodes
	Vector<Complex>& zbvj_0, // IO: boundary condition for the middle face
	Matrix<double>& crdelj   // IO: nodal coordinates
)
{
	int i, nrdir=0, i0=0, nod1, nod2, ntime;
	double xi;
    Complex zbl;
	int nrefdi[] = {0, 0, 1, 0, 2, 0, 1, 0};

	if(kref == 0)
	{
		return;
	}
	else if(kref < 0 || kref > 7)
			NC_Error_Exit_0("kref must be equal to 0, 1, ... 7!");

	// direction of the current reflection
	switch(numReflectionsOfElements_)
	{
	case 8:  // there are 3 symmetric planes
		nrdir = nrefdi[kref];
		break;
	case 4:  // there are 2 symmetric planes
		for(i=0; i<NDIM; i++) if(listSymmetryPlanes[i] == 0) 
		{
			i0 = i + 1;
			break;
		}
		nrdir = i0 + nrefdi[kref];
		if(nrdir >= NDIM) nrdir -= NDIM;
		break;
	case 2:  // there is 1 symmetric plane
		for(i=0; i<NDIM; i++) if(listSymmetryPlanes[i] != 0)
		{
			nrdir = i;
			break;
		}
		break;
	case 1:  // there is not symmetric plane
		return;
		break;
	default:
		NC_Error_Exit_0(NCout, "numReflectionsOfElements_ must be equal to 8, 4, 2 or 1!");
	}

	// update the reflection factor and the velocity boundary condition in case of an antisymmetric reflection
	if(listSymmetryPlanes[nrdir] < 0)
	{
		ireffac = -ireffac;
		if(ibvj_0 == 0) for(i=0; i<inodj; i++) zbvj_0[i].nega();
	}

	// update the center, the normal vectors and the nodal coordinates of the element
	centerj[nrdir] = 2.0*listSymmetryPlanesCoordinates[nrdir] - centerj[nrdir];
	norvecj[nrdir] = -norvecj[nrdir];
	for(i=0; i<inodj; i++)
	{
		norvecj[NDIM*(i + 1) + nrdir] = -norvecj[NDIM*(i + 1) + nrdir];
		crdelj(i, nrdir) = 2.0*listSymmetryPlanesCoordinates[nrdir] - crdelj(i, nrdir);
	}

	// update the sequence of the nodes in the local coordinate system of the element
	nod1 = 0;
	nod2 = 1;
	ntime = 0;
IterBegin_rfe:
	ntime++;
	for(i=0; i<NDIM; i++) 
	{
		xi = crdelj(nod1, i);
		crdelj(nod1, i) = crdelj(nod2, i);
		crdelj(nod2, i) = xi;
		xi = norvecj[NDIM*(nod1 + 1) + i];
		norvecj[NDIM*(nod1 + 1) + i] = norvecj[NDIM*(nod2 + 1) + i];
		norvecj[NDIM*(nod2 + 1) + i] = xi;
	}
	zbl = zbvj_0[nod1];
	zbvj_0[nod1] = zbvj_0[nod2];
	zbvj_0[nod2] = zbl;
	if(inodj == NETYP4 && ntime == 1) 
	{
		nod1 = 2;
		nod2 = 3;
		goto IterBegin_rfe;
	}
}

// function for computing the number of Gaussean points,their local coordinates and weights, return: number of Gaussean point
int NC_ComputeGausseanPoints
(
    double* csigau,    // O: r-coordinate of the Gaussean points
    double* etagau,    // O: s-coordinate of the Gaussean points
    double* weigau,    // O: weight of the Gaussean points
    const int& inodj,  // I: number of nodes of the element
    const int& norder  // I: order of the integration
)
{
	int ngp=0, i, j;
	Vector<double> coordgau(norder), weitgau(norder);

	switch(inodj) 
	{
	case 3:
		switch(norder)
		{
		case 1:
			ngp = 1;
			csigau[0] = GAUCORWEI_TR1[0][0];
			etagau[0] = GAUCORWEI_TR1[0][1];
			weigau[0] = GAUCORWEI_TR1[0][2];
			break;
		case 2:
			ngp = 4;
			for(i=0; i<ngp; i++) 
			{
				csigau[i] = GAUCORWEI_TR4[i][0];
				etagau[i] = GAUCORWEI_TR4[i][1];
				weigau[i] = GAUCORWEI_TR4[i][2];
			}
			break;
		case 3:
			ngp = 7;
			for(i=0; i<ngp; i++) 
			{
				csigau[i] = GAUCORWEI_TR7[i][0];
				etagau[i] = GAUCORWEI_TR7[i][1];
				weigau[i] = GAUCORWEI_TR7[i][2];
			}
			break;
		case 4:
		case 5:
		case 6:
		case 7:
			ngp = 13;
			for(i=0; i<ngp; i++) 
			{
				csigau[i] = GAUCORWEI_TR13[i][0];
				etagau[i] = GAUCORWEI_TR13[i][1];
				weigau[i] = GAUCORWEI_TR13[i][2];
			}
			break;
		default:
			NC_Error_Exit_0(NCout, "norder must be between 1 and 7!");
		} // end id SWITCH NORDER
		break;
	case 4:

		BLGauCooWe(NCout, norder, coordgau, weitgau);

		ngp = 0;
		for(i=0; i<norder; i++) for(j=0; j<norder; j++)
		{
			csigau[ngp] = coordgau[i];
			etagau[ngp] = coordgau[j];
			weigau[ngp] = weitgau[i]*weitgau[j];
			ngp++;
		}
		break;
	default:
		NC_Error_Exit_0(NCout, "inodj must be equal to 3 or 4!");
	} // end of SWITCH INODJ

	return ngp;
}

// compute parameters at a integration point, return: Jacobian
double NC_ComputeParameters
(
    Vector<double>& shapfn,  // O: shape functions
    Vector<double>& shapds,  // O: dh/ds
    Vector<double>& shapdt,  // O: dh/dt
    Vector<double>& dxdxi,   // O: dx/ds
    Vector<double>& dxdet,   // O: dx/dt
    Vector<double>& crdgpi,  // O: global coordinates of the point
    Vector<double>& elnorm,  // O: unit normal vector to the el.
    Complex& zbvgau_0,       // O: boundary condition of under face at the point
    Matrix<double>& crden,   // I: nodal coordinates of the el.
    const int& ifcrhs,       // I: flag for computing velocity boundary condition
    Vector<Complex>& zbvj_0, // I: value of boundary condition of under face
    const int& inode,        // I: number of nodes of the el.
    const double& ss,        // I: local coordinate s
    const double& tt         // I: local coordinate t
)
{
	int i;
	double vjaco;

	// compute the global coordinate of the integration point
	NC_ComputeGlobalCoordinates(crdgpi, shapfn, crden, ss, tt, inode);

	// compute the normal vector and the Jacobian
    vjaco = NC_ComputeUnitVectorOnElement(elnorm, shapds, shapdt, dxdxi, dxdet, crden,
				ss, tt, inode);

	// boundary conditions at the point
	zbvgau_0.set(0.0, 0.0);
	if(ifcrhs)
	{
		for(i=0; i<inode; i++)
		{
			zbvgau_0 += zbvj_0[i]*shapfn[i]; 
		}
	}

	return vjaco;
}

// generate the subelements, return: number of subelements
int NC_GenerateSubelements
(
	ofstream& NCout,
	const int& ndip,        // I: number of the sorce node
	Vector<double>& crdip,  // I: coordinates of the source point //x
	const int& je,          // I: internal number of the element
	Matrix<double>& crdelj, // I: coordinates of nodes of the element //y
	const double& arelj,    // I: area of the element 
	const int& nvertj,      // I: number of vertexes of the element
	Vector<double>& csij,   // I: local coordinates of nodes of the el. 
	Vector<double>& etaj,   // I: local coordinates of nodes of the el.
	const int& jnode,       // I: number of nodes of the el.
	Vector<double>& xisbp,  // O: the first local coourdinate of the //center points of subelements
			                //    original point of each subelement
	Vector<double>& etsbp,  // O: the second local coourdinate of the
			                //    original point of each subelement 
	Vector<double>& fasbp,  // O: dimension factor of each subelement
	Vector<int>& ngsbp      // O: Gaussean order of each subelement
)
{
	const int MNSE = 60, MVX = 4, NSE = 4; //MNSE maximale Anzahl an unterteilungen //NSE factor of subdivision level

	int nsbe = -1, nsfl = 1, i, j, ndie, nsel, idi, nsu, j1, j2, jr = 0, nsf0,
		ngaumax = 6, ngaumin = 2;
	double tolf = 1.3,//toleranz f�r verh�ltnis Abstand Fl�che(subelement)
    faclin = 2.0, arels, scent, tcent, ratdis, xi, et, dxi = 0.0, det,
		disfac, err_g, err_h, err_e, gauaccu = 0.001;
	Vector<double> crdpoip(NDIM), xisp(MVX*2), etsp(MVX*2);
	Matrix<double> xisfp(MNSE, MVX), etsfp(MNSE, MVX), xisep(MNSE, MVX),
		etsep(MNSE, MVX);
	Vector<double> shfunx(NNODPE);

	for(j=0; j<nvertj; j++)
	{
		xisfp(0, j) = csij[j];
		etsfp(0, j) = etaj[j];
	}

	// loop using GOTO
Begin_Subel:
	ndie = 0;      // number of elments to be subdivided in this loop //subdivision level (kantesb = kante 2^-j)
	faclin *= 0.5; // length of the subelements / length of the original el. //relative Kantenl�nge der subelements 2^-j
	arels = arelj*faclin*faclin; // area of the current subelements
	nsel = nsfl;   // number of subelements to be treated in the current loop

	// local coordinates of the vertices of the current subelements
	for(idi=0; idi<nsel; idi++) for(j=0; j<nvertj; j++)
	{
		xisep(idi, j) = xisfp(idi, j);
		etsep(idi, j) = etsfp(idi, j);
	}

	// loop over subelements of the current depth of divides
	idi = 0; //index of divisions
Begin_Idi:

	// global coordinates of the center of the current subelement
	scent = tcent = 0.0;
	for(j=0; j<nvertj; j++)
	{
		scent += xisep(idi, j); //local coordinates of subelement center
		tcent += etsep(idi, j);
	}
	scent /= (double)nvertj;
	tcent /= (double)nvertj;
	NC_ComputeGlobalCoordinates(crdpoip, shfunx, crdelj, scent, tcent, jnode); //s t sind gewichte der Shapefunctions

	// ratio of the distace between the source point and the center of the current subelement to the size of the subelement
	ratdis = BLdispoi(crdpoip, crdip)/sqrt(arels); //l2 norm der vektoren

	if(ratdis < tolf) // the current subelement must be further subdivided
	{
		ndie++;
		if(ndie > 15)
		{
			cout << endl;
			if(ndip >= 0)
				cout << "number of the source node    = " << extNumbersOfNodes[ndip] << endl;
			cout << "coordinates of the source point: " << 
				crdip[0] << " " << crdip[1] << " " << crdip[2] << endl;
			NCout << endl;
			if(ndip >= 0)
				NCout << "number of the source node    = " << extNumbersOfNodes[ndip] << endl;
			NCout << "coordinates of the source point: " << 
				crdip[0] << " " << crdip[1] << " " << crdip[2] << endl;
			NC_Error_Exit_1(NCout, 
				"number of subels. which are subdivided in a loop must <= 15!",
				"number of the element = ", extNumbersOfElements[je]);
		}
		nsfl = ndie*NSE;
		nsf0 = nsfl - NSE;
		
		// generate the local coordinates of the midside points
		for(j=0; j<nvertj; j++)
		{
			j1 = j + 1;
			if(j1 == nvertj) j1 = 0;
			xisp[j] = xisep(idi, j);
			xisp[j + nvertj] = (xisep(idi, j) + xisep(idi, j1))/2.0;
			etsp[j] = etsep(idi, j);
			etsp[j + nvertj] = (etsep(idi, j) + etsep(idi, j1))/2.0;
		}

		// generate 4 subelements
		for(j=0; j<nvertj; j++)
		{
			nsu = nsf0 + j;
			j1 = j + nvertj;
			j2 = j1 - 1;
			if(j2 < nvertj) j2 += nvertj;
			xisfp(nsu, 0) = xisp[j];
			xisfp(nsu, 1) = xisp[j1];
			etsfp(nsu, 0) = etsp[j];
			etsfp(nsu, 1) = etsp[j1];
			if(nvertj == NETYP4) // for a 4 side element
			{
				xisfp(nsu, 2) = scent;
				xisfp(nsu, 3) = xisp[j2];
				etsfp(nsu, 2) = tcent;
				etsfp(nsu, 3) = etsp[j2];
			}
			else // for triangle element
			{
				xisfp(nsu, 2) = xisp[j2];
				etsfp(nsu, 2) = etsp[j2];
				xisfp(nsf0 + NSE - 1, j) = xisp[j1];
				etsfp(nsf0 + NSE - 1, j) = etsp[j1];
			}
		} // end of loop J
	}
	else // the current subelement need not be further subdivided, store the results
	{
		nsbe++;
		if(nsbe >= MSBE)
		{
			if(ndip >= 0)
			{
				cout << "number of the source node    = " << extNumbersOfNodes[ndip] << endl;
				NCout << "number of the source node    = " << extNumbersOfNodes[ndip] << endl;
			}
			cout << "coordinates of the source point: " << 
				crdip[0] << " " << crdip[1] << " " << crdip[2] << endl;
			NCout << "coordinates of the source point: " << 
				crdip[0] << " " << crdip[1] << " " << crdip[2] << endl;
			NC_Error_Exit_1(NCout, "nsbe must be < MSBE(= 110)", 
				"number of the element = ", extNumbersOfElements[je]);
		}

		if(nvertj == NETYP4) // for a 4 side el.
		{
			// store the local coordinates of the center of the subelement
			xisbp[nsbe] = (xisep(idi, 0) + xisep(idi, 1) + xisep(idi, 2) + 
				xisep(idi, 3))/4.0;
			etsbp[nsbe] = (etsep(idi, 0) + etsep(idi, 1) + etsep(idi, 2) + 
				etsep(idi, 3))/4.0;

			// store the ratio of the dimension 
			fasbp[nsbe] = faclin;
		}
		else // for triangle el.
		{
			// local coordinates of the center of the subel.
			xi = (xisep(idi, 0) + xisep(idi, 1) + xisep(idi, 2))/3.0;
			et = (etsep(idi, 0) + etsep(idi, 1) + etsep(idi, 2))/3.0; 

			// seek the vertex with the right angel in the local coordinate system
			for(j=0; j<nvertj; j++)
			{
				dxi = xi - xisep(idi, j);
				det = et - etsep(idi, j);
				if(fabs(dxi - det) < 1.0e-8)
				{
					jr = j;
					break;
				}
			}

			// store the local coordinates of the vertex with the right angle
			xisbp[nsbe] = xisep(idi, jr);
			etsbp[nsbe] = etsep(idi, jr);

			// store the ratio of the dimension
			if(dxi >= 0.0) fasbp[nsbe] = faclin; else fasbp[nsbe] = -faclin; //2^-j
		}

		// computing the Gaussean order
		disfac = 0.5/ratdis;//fl�che(subelement)/entfernung/2
		ngsbp[nsbe] = ngaumax;
		for(i=ngaumin; i<=ngaumax; i++)
		{
			Err_G_(i, disfac, j, err_g);
			Err_H_(i, disfac, j, err_h);
			Err_E_(i, disfac, j, err_e);
			if(err_e < gauaccu && err_h < gauaccu && err_g < gauaccu)
			{
				ngsbp[nsbe] = i;
				break;
			}
		}
	} // end of store results

	idi++;
	if(idi < nsel) goto Begin_Idi;

	if(ndie > 0) goto Begin_Subel;


	return (nsbe + 1);
}

// compute the unit vector and weight at each integration point of a unit sphere surface
void NC_ComputeUnitVectorOnUnitSphere
(
	ofstream& NCout,
	const int& np_thet,  // I: number of integration points in the theta-direction
	const int& np_phi,   // I: number of integration points in the phi-direction
	double* weisphe,	 // O: weight of each integration point
	double** uvcsphe     // O: coodinates of each integration point
)
{
	int i, j, k;
	double phi, cothet, sithet, delphi;
	Vector<double> crdgpth(np_thet), wegpth(np_thet);

	// compute the coordinates and weights of the Gaussean points
	BLGauCooWe(NCout, np_thet, crdgpth, wegpth);
	delphi = 2.0*PI/(double)(np_phi); 

	// compute the unit vector at each Gaussean point of the unit sphere
	k = 0;
	for(i=0; i<np_thet; i++)
	{
		cothet = crdgpth[i];
		sithet = sqrt(1.0 - cothet*cothet);
		for(j=0; j<np_phi; j++)
		{
			weisphe[k] = wegpth[i]*delphi/PI4;
			phi = delphi*(double)(j);
			uvcsphe[k][0] = sithet*cos(phi);
			uvcsphe[k][1] = sithet*sin(phi);
			uvcsphe[k++][2] = cothet;
		}
	}
}

// compute the spherical Hankel funktions of the first kind under a given order (see [1], eq. (19))
void NC_SphericalHankel
(
	ofstream& NCout,
	Vector<Complex>& zH_n,	// O: spherical Hankel function of the first kind and of orders from 0 to norder-1
	const int& norder,      // I: highst order + 1
	const double& argu      // I: argument
)
{
    if(norder < 2) NC_Error_Exit_1(NCout,
                              "Highst order the spherical Hankel funktions must greater than 1",
                              "Highst order = ", norder);
    
    int i, j=1;
    double dnu, gnu, errgn=1.0, di, Dj, Cj, aj, bj;
    double epsign = 1.0e-9, dp;
    Vector<double> dgn(norder), ggn(norder), dY_n(norder);
    
    // compute the imaginal part of the spherical Hankel function
    dY_n[0] = -cos(argu)/argu;
    dY_n[1] = -(cos(argu)/argu + sin(argu))/argu;
    
    for(i=2; i<norder; i++) dY_n[i] = dY_n[i - 1]*(double)(2*i - 1)/argu -
        dY_n[i - 2];
    
    // compute the real part of the spherical Hankel function (see K. Giebermann, "Schnelle Sumationsverfahren zur numerischen L�sung von Integralgleichungen f�r Streuprobleme im R^(3)", pp 87-88)
    dnu = (double)(norder - 1);
    
    di=(double)(2*(dnu+1)+1)/argu;
    Cj=di;
    Dj=0.0;
    while(errgn > epsign)
    {
        aj=-1.0;
        bj=(double)(2*(dnu+j+1)+1)/argu;
        Dj=bj+aj*Dj;
        if(Dj==0){Dj=1e-30;}
        Dj=1/Dj;
        Cj=bj+aj/Cj;
        if(Cj==0){Cj=1e-30;}
        di=di*Cj*Dj;
        errgn = fabs(Cj*Dj-1.0);
        j=j+1;
    }
    
    gnu=dnu/argu-1/di;
    
    ggn[norder - 1] = 1.0;
    dgn[norder - 1] = gnu;
    for(i=norder-2; i>=0; i--)
    {
        di = (double)(i);
        ggn[i] = (di + 2.0)/argu*ggn[i + 1] + dgn[i + 1];
        dgn[i] = di/argu*ggn[i] - ggn[i + 1];
    }
    
    if(fabs(ggn[0]) > 1.0e-5)
    {
        dp = sin(argu)/argu/ggn[0];
    }
    else
    {
        dp = (cos(argu) - sin(argu)/argu)/argu/dgn[0];
    }
    
    for(i=0; i<norder; i++) zH_n[i].set(dp*ggn[i], harmonicTimeFactor_*dY_n[i]);
}

// compute the Legendre polynomes under a given order
void NC_LegendrePolynomes
(
	ofstream& NCout,
	Vector<double>& PLe,	// O: values of the Legendre polynomes
	const int& norder,      // I: highst order + 1
	const double& argu      // I: argument
)
{
	if(norder < 2) NC_Error_Exit_1(NCout, 
		"Highst order of the Legendre polynomes must greater than 1",
		"Highst order = ", norder);

	PLe[0] = 1.0;
	PLe[1] = argu;

	for(int i=2; i<norder; i++) PLe[i] = (PLe[i - 1]*argu*(double)(2*i - 1) -
		PLe[i - 2]*(double)(i - 1))/(double)(i);
}