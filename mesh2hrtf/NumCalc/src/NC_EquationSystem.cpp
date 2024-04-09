/*===========================================================================*\
 *                                                                            *
 *  File name:      NC_EquationSystem.cpp                                     *
 *  Description:    set up equation system for BEM                            *
 *  Author:         H. Ziegelwanger, W. Kreuzer and Z.S. Chen                 *
 *                                                                            *
 \*===========================================================================*/



//================================= I N C L U D E S ==================================
// local includes                                                                   //
#include"NC_ConstantsVariables.h"                                                   //
#include"NC_TypeDefinition.h"                                                       //
#include"NC_Arrays.h"                                                               //
#include"NC_Macros.h"                                                               //
#include"NC_IntegrationConstants.h"                                                 //
//                                                                                  //
// system includes                                                                  //
#include<iostream>                                                                  //
#include<iomanip>                                                                   //
#include<string>                                                                    //
#include<fstream>  // file class functions                                          //
#include<math.h>                                                                    //
using namespace std;                                                                //
//====================================================================================



// local functions
void NC_BuildSystemTBEM(ofstream&);
void NC_AssembleTBEM(ofstream&,const int&,Vector<int>&,const int&);
void NC_ComputeEntriesForTBEM(ofstream&,const int&);
void NC_BuildSystemFMBEM(ofstream&);
void NC_BuildTforSLFMM(ofstream&);
void NC_BuildDforSLFMM(ofstream&);
void NC_BuildSforSLFMM(ofstream&);
void NC_BuildNforSLFMM(ofstream&);
void NC_BuildTforMLFMM(ofstream&);
void NC_BuildDforMLFMM(ofstream&);
void NC_BuildSforMLFMM(ofstream&);
void NC_BuildNforMLFMM(ofstream&);
void NC_RHSvecTofSLFMM(ofstream&);
void NC_RHSvecTofMLFMM(ofstream&);
void NC_ComputeEntriesTforFMM(ofstream&,const int&,const int&,int&,int&,double*,double**,Vector<Complex>&,Vector<Complex>&);
void NC_ComputeEntriesNforFMM(ofstream&,const int&,const int&,Vector<int>&);
void NC_AssembleNforFMM(ofstream&,const int&,const int&,Vector<int>&,int&);
void NC_SearchAddressNforFMM(ofstream&,int&,int&,int&,const char&);
void NC_SingularIntegration(ofstream&,Vector<Complex>&,const int&,const int&,Vector<Complex>&,Matrix<double>&);
void NC_RegularIntegration(ofstream&,Vector<Complex>&,const int&,const int&,const int&,const int&,Vector<Complex>&,Vector<double>&,Matrix<double>&);
Complex NC_IncidentWaveRHS(ofstream&);
double NC_FreqCurveInterpolation(ofstream&, int&);



// global variables or constants
int Ifadmii_3;
Complex Admia3, Admib3, Zbvi03;
Vector<double> Sourpoi3(NDIM);
Vector<int> Souradr3(2);
int Ibvj_03, Ifadmj3, Inoj3;
Vector<int> Jdofaddr3(20); //MORDER2D*2
Vector<double> Centej3(NDIM), Norvcj3(NDIM*(NNODPE + 1));
Vector<Complex> Admj3(2);
Vector<Complex> Zbvj03(NNODPE);
Matrix<double> Crdej3(NNODPE, NDIM);
int Ibvi03;
Vector<double> Norvci3(NDIM*(NNODPE + 1));
Vector<Complex> Zrsintel(6);
int Ifcrh3;
double Gama3;
Complex zBta3;
const int GauFMBOrd3 = 2;
double csiFMgau[GauFMBOrd3*GauFMBOrd3], etaFMgau[GauFMBOrd3*GauFMBOrd3], weiFMgau[GauFMBOrd3*GauFMBOrd3];



// the main program for seting up the equation system
void NC_SetupEquationSystem
(
	ofstream& NCout
)
{
	// set the coupling constant
	Gama3 = 1.0;
	if(isInternalProblem_) {zBta3.set(0.0, 0.0);} else {zBta3.set(0.0, harmonicTimeFactor_/waveNumbers_);}

	switch(methodFMM_)
	{
	case 0: // TBEM
		NC_BuildSystemTBEM(NCout);
		break;
	case 1: // SLFMBEM
	case 3: // DMLFMBEM
		NC_BuildSystemFMBEM(NCout);
	}
}

// generate the equation system of the traditional BEM (TBEM)
void NC_BuildSystemTBEM
(
	ofstream& NCout
)
{
	int i, iel;

	// initialize the near field coefficient matrix
	for(i=0; i<numComponentsOfCoefficientMatrix_; i++) zcoefl[i].set(0.0, 0.0);

	// loop over source elements
	for(iel=0; iel<numElements_; iel++)
	{
		if(listElementProperty[iel] == 2) continue;

		// center point and normal vector
		for(i=0; i<NDIM; i++)
		{
			Sourpoi3[i] = centel[iel][i];
			Norvci3[i] = elenor[iel][i];
		}

		// addresse of the unknown DOFs of the source element
		Souradr3[0] = jelist[iel][0];
		Souradr3[1] = jelist[iel][1];

		// boundary condition
		switch(ibval[iel]) 
		{
		case 0:         // velocity prescribed
			Ibvi03 = 0;  
			Ifadmii_3 = 0;
			break;
		case 1:         // pressure prescribed
			Ibvi03 = 1;
			Ifadmii_3 = 0;
			break;
		case 2:         // velocity and surface admittance prescribed
			Ibvi03 = 0;  
			Ifadmii_3 = 1;
			break;
		case 3:         // transfer admittance prescribed
			Ibvi03 = 2;
			Ifadmii_3 = 0;
			break;
		case 4:         // transfer and surface admittances prescribed
			Ibvi03 = 2;
			Ifadmii_3 = 1;
			break;
		}
		Zbvi03 = 0.0;
		for(i=0; i<listNumberNodesPerElement[iel]; i++) Zbvi03 += zbval0[iel][i];
		Zbvi03 /= (double)listNumberNodesPerElement[iel];

		// three cases
		if(listElementProperty[iel] == 0) // the element is a surface element
		{
			if(Ifadmii_3) {Admia3 = zbval1[iel];} else {Admia3.set(0.0, 0.0);}
			NC_ComputeEntriesForTBEM(NCout, iel);
		}

	} // end of loop IEL
}

// assemble the contributions of an element to the coefficient matrix and the right hand side vector (TMEM)
void NC_AssembleTBEM
(
	ofstream& NCout,
	const int& jel,			// I: number of the element
	Vector<int>& Souradr3,	// I: address of the DOF of the source el.
	const int& irefel
 )
{
    int ijadr00, ijadr01, ijadr10, ijadr11;
    Complex zsum0, zsum1;
    double reffact;
    
    ijadr00 = Souradr3[0]*numRowsOfCoefficientMatrix_ + Jdofaddr3[0];
    ijadr01 = Souradr3[0]*numRowsOfCoefficientMatrix_ + Jdofaddr3[1];
    ijadr10 = Souradr3[1]*numRowsOfCoefficientMatrix_ + Jdofaddr3[0];
    ijadr11 = Souradr3[1]*numRowsOfCoefficientMatrix_ + Jdofaddr3[1];
    
    if(irefel == -1) {reffact = -1.0;} else {reffact = 1.0;}
    
    if(Ibvj_03 == 0)
    {
        zsum0 = Zrsintel[1]*(Gama3*Tao_) + Zrsintel[3]*zBta3;
        if(Ifadmj3) zsum0 += (Zrsintel[0]*Gama3 + Zrsintel[2]*(zBta3*Tao_))*Admj3[0];
        
        if(Ifcrh3)
        {
            zrhs[Souradr3[0]] += Zrsintel[4];
        }
    }
    else if(Ibvj_03 == 1)
    {
        zsum0 = Zrsintel[0]*(-Gama3*Tao_) - Zrsintel[2]*(zBta3);
        if(Ifcrh3)
            zrhs[Souradr3[0]] += Zrsintel[4];
    }
    else
    {
        NC_Error_Exit_2(NCout, "Ibvj_03 must be equal to 0 or 1!", "Ibvj_03 = ", Ibvj_03,
                        "Element number = ", extNumbersOfElements[jel]);
    }
    
    zcoefl[ijadr00] += zsum0*reffact;
}

// compute a row of the equation system in case that the source point is located in a surface element (see eq. (13) in [1])
void NC_ComputeEntriesForTBEM
(
	ofstream& NCout, 
	const int& iel           // number of the source element
)
{
	int jel, kref, i, ireffac;
	Complex z0, Admib3(0.0, 0.0);

	// compute of the free terms of the coefficient matrix
	switch(Ibvi03)
	{
	case 0:
		zcoefl[Souradr3[0]*numRowsOfCoefficientMatrix_ + Souradr3[0]] += (Admia3*zBta3 - Gama3)*0.5;
		break;
	case 1:
		zcoefl[Souradr3[0]*numRowsOfCoefficientMatrix_ + Souradr3[0]] -= zBta3*(0.5*Tao_);
		break;
	default:
		NC_Error_Exit_2(NCout, "Ibvi03 must be equal to 0 or 1!", "Ibvi03 = ", Ibvi03,
			"element number = ", extNumbersOfElements[iel]);
	}

	// compute the free terms of the right hand side vector
	z0 = 0.0;

	// compute the contribution of the velocity boundary condition
	if(Ibvi03 == 0) // velocity B. C. prescribed
	{
		z0 += Zbvi03*zBta3*Tao_*0.5;
    }
    if(Ibvi03 == 1) // pres. B.C. prescribed
    {
        z0 += Zbvi03*Tao_*0.5;
    }

	// compute the contributions of the incident waves
	if(numIncidentPlaneWaves_ + numPointSources_ > 0)
	{
		z0 += NC_IncidentWaveRHS(NCout);
	}
	zrhs[Souradr3[0]] += z0;

	// loop over the field elements
	for(jel=0; jel<numElements_; jel++)
	{
		if(listElementProperty[jel] == 2) continue;

		// compute the element data 
		NC_ComputeElementData(jel, 0, Inoj3, Ibvj_03,
			Ifadmj3, Jdofaddr3, Centej3, Norvcj3, Admj3,
			Zbvj03, Crdej3);

		// set the control parameter for computing the contributions to the r.h.s.
		Ifcrh3 = 0;
		if((Ibvj_03 == 0 || Ibvj_03 == 1) && listElementProperty[jel] == 0)
		{
			for(i=0; i<Inoj3; i++)
			{
				if(Zbvj03[i].norm() > EPSY)
				{
					Ifcrh3 = 1;
					break;
				}
			}
		}

		// loop over the reflections of the element with respect to the symmetry/antisymmetry planes
		ireffac = 1;
		for(kref=0; kref<numReflectionsOfElements_; kref++)
		{
			// reflection of the element
			if(kref > 0) NC_ReflectElement(kref, Inoj3, Ibvj_03, ireffac, Centej3,
				Norvcj3, Zbvj03, Crdej3);

			if(kref == 0 && jel == iel)
			{
				// singular integration
				NC_SingularIntegration(NCout, Zrsintel, Inoj3, Ifcrh3, Zbvj03,
					Crdej3);
			}
			else
			{
				// regular integration
				NC_RegularIntegration(NCout, Zrsintel, iel, jel, Inoj3,
					Ifcrh3, Zbvj03, Centej3, Crdej3);
			}

			// assemble the contributions of the element to the equation system
			NC_AssembleTBEM(NCout, jel, Souradr3, ireffac);
		} // end of loop KREF
	} // end of loop JEL
}

// generate the fast multipole BEM equation system
void NC_BuildSystemFMBEM
(
	ofstream& NCout
 )
{
    int nlv;
    
    switch(methodFMM_)
    {
        case 1: // SLFMBEM
            // compute the near field coefficient matrix
            NC_BuildNforSLFMM(NCout);
            
            // compute the coordinates and weights of integration points on the unit sphere
            NC_ComputeUnitVectorOnUnitSphere(NCout, numIntegrationPointsThetaDirection_, numIntegrationPointsPhiDirection_, weisphe, uvcsphe);
            
            // compute the D-matrix or D-matrices
            NC_BuildDforSLFMM(NCout);
            
            // compute the T-matrix and the T-vector
            NC_BuildTforSLFMM(NCout);
            
            // compute the S-matrix
            NC_BuildSforSLFMM(NCout);
            
            break;
        case 3: // DMLFMBEM
            // compute the near field coefficient matrix
            NC_BuildNforMLFMM(NCout);
            
            // compute the coordinates and weights of integration points on the unit sphere
            for(nlv=0; nlv<numClusterLevels_; nlv++)
            {
                NC_ComputeUnitVectorOnUnitSphere(NCout, clulevarry[nlv].nPoinThetLv, clulevarry[nlv].nPoinPhiLv,
                                                 clulevarry[nlv].weisphe, clulevarry[nlv].uvcsphe);
            }
            
            // compute the D-matrices
            NC_BuildDforMLFMM(NCout);
            
            //weisphe = clulevarry[nlevtop_].weisphe;
            //uvcsphe = clulevarry[nlevtop_].uvcsphe;
            
            // compute the T-matrices the T-vectors
            NC_BuildTforMLFMM(NCout);
            
            // compute the S-matrices
            NC_BuildSforMLFMM(NCout);
    }
    
    // compute the contribution of the T-vector to the right hand side vector ({r} += [S]*[D]*{t})
    if(boolComputeTVector_) 
    {
        switch(methodFMM_)
        {
            case 1: // SLFMBEM
                NC_RHSvecTofSLFMM(NCout);
                break;
            case 3: // DMLFMBEM
                NC_RHSvecTofMLFMM(NCout);
                break;
        }
    }
}

// compute the T-matrix and the T-vector for SL-FMM
void NC_BuildTforSLFMM
(
	ofstream& NCout
 )
{
    int i, j, k, nelj, irownu, ientry, iadrek, ndofi, inod0 = 0, idfs, ncldofi, jdifdar;
    double *cencli;
    int ifvcj, nfmgp;
    
    int nrf_cl, rfa_cl;
    double dfa_cl;
    bool mir_cl;
    Vector<bool> rfdi_cl(NDIM);
    Complex *zt_mtx;
    
    // vectors for storing contributions of the element to the T-vector and the T-matrix
    Vector<Complex> TVele(numIntegrationPointsUnitSphere_), TMele(numIntegrationPointsUnitSphere_);
    
    // initialize the T-vector
    if(boolComputeTVector_) for(i=0; i<numIntegrationPointsUnitSphere_*numOriginalReflectedClusters_; i++) ztvct[i].set(0.0, 0.0);
    
    // define the array used to store the t-matrix
    zt_mtx = ztmtx;
    
    // initialize the T-matrix
    for(i=0; i<irowtmtx[numIntegrationPointsUnitSphere_*numOriginalReflectedClusters_]; i++) zt_mtx[i].set(0.0, 0.0);
    
    // loop over the clusters
    irownu = ientry = 0;
    for(i=0; i<numOriginalReflectedClusters_; i++)
    {
        // cluster parameters
        idfs = ClustArray[i].NDOFsPeEl;
        nrf_cl = ClustArray[i].nuref;
        rfa_cl = ClustArray[i].rffac;
        if(nrf_cl)
        {
            mir_cl = ClustArray[i].ifmirro;
            for(j=0; j<NDIM; j++) rfdi_cl[j] = ClustArray[i].ifrfdi[j];
        }
        
        // compute number of gaussean points and their coordinates and weights
        if(ClustArray[i].IfMonoEl)
        {
            inod0 = listNumberNodesPerElement[ClustArray[i].NumsOfEl[0]];
            nfmgp = NC_ComputeGausseanPoints(csiFMgau, etaFMgau, weiFMgau, inod0, GauFMBOrd3);
        }
        
        ndofi = ClustArray[i].NDOFsPeEl;   // number of unknowns of per element
        ncldofi = ClustArray[i].NumOfDOFs; // total number of unknowns of the cluster
        cencli = ClustArray[i].CoorCent;   // center of the cluster
        
        // loop over elements of the cluster
        for(j=0; j<ClustArray[i].NumOfEl; j++)
        {
            nelj = ClustArray[i].NumsOfEl[j];
            
            // compute number of gaussean points and their coordinates and weights
            if(!ClustArray[i].IfMonoEl && listNumberNodesPerElement[nelj] != inod0)
            {
                inod0 = listNumberNodesPerElement[nelj];
                nfmgp = NC_ComputeGausseanPoints(csiFMgau, etaFMgau, weiFMgau, inod0, GauFMBOrd3);
            }
            
            // compute the element data
            NC_ComputeElementData(nelj, 0, Inoj3, Ibvj_03,
                                  Ifadmj3, Jdofaddr3, Centej3, Norvcj3, Admj3,
                                  Zbvj03, Crdej3);
            
            // reflection of the element
            if(nrf_cl) NC_ReflectElementFMM(rfa_cl, mir_cl, rfdi_cl, Inoj3, Ibvj_03,
                                            Centej3, Norvcj3, Zbvj03, Crdej3);
            
            // see if the velocity boundary condition of the current element should be
            // considered
            ifvcj = 0;
            if(boolComputeTVector_ && (Ibvj_03 == 0 || Ibvj_03 == 1))
            {
                for(k=0; k<Inoj3; k++) if(Zbvj03[k].norm() > 0.0)
                {
                    ifvcj = 1;
                    break;
                }
            }
            
            // initialize the output vectors
            for(k=0; k<numIntegrationPointsUnitSphere_; k++)
            {
                TMele[k].set(0.0, 0.0);
                TVele[k].set(0.0, 0.0);
            }
            
            if(Ifadmj3) {Admia3 = zbval1[nelj];} else {Admia3.set(0.0, 0.0);}
            NC_ComputeEntriesTforFMM(NCout, nelj, ifvcj, numIntegrationPointsUnitSphere_, nfmgp,
                          cencli, uvcsphe, TMele, TVele);
            
            // assemble the contributions of the element to the T-matrix and the T-vector
            jdifdar = j*ndofi;
            if(rfa_cl == -1) {dfa_cl = -1.0;} else {dfa_cl = 1.0;}
            
            for(k=0; k<numIntegrationPointsUnitSphere_; k++)
            {
                iadrek = ientry + k*ncldofi + jdifdar;
                
                zt_mtx[iadrek] = TMele[k]*dfa_cl;
                if(idfs == 2)  zt_mtx[++iadrek] = TVele[k]*dfa_cl;
                
                if(ifvcj) ztvct[irownu + k] += TVele[k];
            }
        } // end of loop J
        irownu += numIntegrationPointsUnitSphere_;
        ientry += numIntegrationPointsUnitSphere_*ncldofi;
    } // end of loop I
}


// compute the D-matrix for SL-FMMM  (see eq. (18) in [1])
void NC_BuildDforSLFMM
(
	ofstream& NCout
)
{
	int i, j, jcl, k, jgp, identry, ikjadre;
	double dij, dkij, scprd;
	Complex zmun, zi;
	Vector<Complex> zH1m(numExpansionTerms_), zdiat(numIntegrationPointsUnitSphere_);
	Vector<double> Pm(numExpansionTerms_), disij(NDIM);
	Complex *zwork_d;

	// array to be used to store the D-matrices
	zwork_d = dmtxlev[0].zDmxLv;

	// loop over the original clusters
	identry = 0; // address of the first entry of the submatrices of the D-matrix that correspond to the cluster i

	for(i=0; i<numOriginalClusters_; i++) 
	{
		// compute and store the nonzeros corresponding to the current cluster
		for(j=0; j<ClustArray[i].NumFarClus; j++)
		{
			jcl = ClustArray[i].NumsFarClus[j];
			for(k=0; k<NDIM; k++) disij[k] = ClustArray[i].CoorCent[k] - 
				ClustArray[jcl].CoorCent[k];
			Vcnorm_dim3_(disij, dij);
			dkij = dij*waveNumbers_;

			// compute the Hankel functions
			NC_SphericalHankel(NCout, zH1m, numExpansionTerms_, dkij);

			for(jgp=0; jgp<numIntegrationPointsUnitSphere_; jgp++)
			{
				scprd = Scprod_dim3_(disij, uvcsphe[jgp]);

				// compute the Legendre plynomes
				NC_LegendrePolynomes(NCout, Pm, numExpansionTerms_, scprd);

				// compute the function Mu_n
				zmun.set(0.0, 0.0);
				zi.set(1.0, 0.0);
				for(k=0; k<numExpansionTerms_; k++)
				{
					zmun += zi*zH1m[k]*Pm[k]*(double)(2*k + 1);
					zi.mul_i(harmonicTimeFactor_); // _ruimag! necessary!
				}
				zdiat[jgp] = zmun;
			} // end of loop JGP

			// store the nonzeros into the D-matrix
			for(k=0; k<numIntegrationPointsUnitSphere_; k++)
			{
				ikjadre = identry + k*ClustArray[i].NumFarClus + j;
				zwork_d[ikjadre] = zdiat[k];
			}
		} // end of loop J
		identry += numIntegrationPointsUnitSphere_*ClustArray[i].NumFarClus;
	} // end of loop I
}

// compute the S-matrix for SL-FMM  (see eq. (21) in [1])
void NC_BuildSforSLFMM
(
	ofstream& NCout
)
{
	int icl, iel, nel, i, jip, iterm0;
	double scpr_dv, d0, wavruim = waveNumbers_*harmonicTimeFactor_;
	Complex z0;
	Vector<double> dcenel(NDIM);
	double *norel;
	bool ifadm;
	Vector<Complex> zvcun(numIntegrationPointsUnitSphere_), zvcdun(numIntegrationPointsUnitSphere_);
	Complex *zs_mtx;

	// define the array to store the S-matrix
	zs_mtx = zsmtx;

	// loop over the original clusters
	for(icl=0; icl<numOriginalClusters_; icl++)
	{
		ifadm = ClustArray[icl].IfAdmiBc;

		// loop over elements of the cluster
		for(iel=0; iel<ClustArray[icl].NumOfEl; iel++)
		{
			nel = ClustArray[icl].NumsOfEl[iel]; // number of the element

			// coordinate difference of the element center and the cluster center
			for(i=0; i<NDIM; i++) dcenel[i] = centel[nel][i] - ClustArray[icl].CoorCent[i];

			// unit normal vector at the element center
			norel = elenor[nel];

			// loop over the integral points on the unit sphere
			for(jip=0; jip<numIntegrationPointsUnitSphere_; jip++)
			{
				// scalar product of two vectors
				scpr_dv = Scprod_dim3_(dcenel, uvcsphe[jip]);

				// compute the result vectors
				d0 = scpr_dv*wavruim;
				z0.set(cos(d0), sin(d0));
				z0 *= weisphe[jip];
				zvcun[jip] = z0;

				d0 = Scprod_dim3_(uvcsphe[jip], norel);
				d0 *= wavruim;
				z0.mul_i(d0);
				zvcdun[jip] = z0;
			} // end of loop over JIP

			// store the results
			iterm0 = irowsmtx[jelist[nel][0]];
        
            for(i=0; i<numIntegrationPointsUnitSphere_; i++)
            {zs_mtx[iterm0 + i] = zvcun[i]*Gama3 + zvcdun[i]*(zBta3*Tao_);}

		} // end of loop over IEL
	} // end of loop ICL
}

// compute the near field equation system for SL-FMM
void NC_BuildNforSLFMM
(
	ofstream& NCout
 )
{
    int i, ic, ie, iel, k, k1, nori_clu;
    Vector<bool> bneacluso(numOriginalClusters_);
    Vector<int> iadreclu(numOriginalClusters_);
    
    // initialize the near field coefficient matrix
    for(i=0; i<irownea[numRowsOfCoefficientMatrix_]; i++) zcoefl[i].set(0.0, 0.0);
    
    // loop over source clusters
    for(ic=0; ic<numOriginalClusters_; ic++)
    {
        // compute the array indicating if an original cluster is a near field cluster
        bneacluso = false;
        for(k=0; k<ClustArray[ic].NumNeaClus; k++)
        {
            nori_clu = ClustArray[ClustArray[ic].NumsNeaClus[k]].OriClust;
            if(!bneacluso[nori_clu]) bneacluso[nori_clu] = true;
        }
        
        // compute the adress of the first entry of each near field cluster
        iadreclu = -1;
        k1 = 0;
        for(k=0; k<numOriginalClusters_; k++)
        {
            if(!bneacluso[k]) continue;
            iadreclu[k] = k1;
            k1 += ClustArray[k].NumOfDOFs;
        }
        
        // loop over elements of the cluster
        for(ie=0; ie<ClustArray[ic].NumOfEl; ie++)
        {
            // number of the source element
            iel = ClustArray[ic].NumsOfEl[ie];
            
            // center point and normal vector
            for(i=0; i<NDIM; i++)
            {
                Sourpoi3[i] = centel[iel][i];
                Norvci3[i] = elenor[iel][i];
            }
            
            // row numbers of the unknown DOFs of the source element
            Souradr3[0] = jelist[iel][0];
            Souradr3[1] = jelist[iel][1];
            
            // boundary condition
            switch(ibval[iel])
            {
                case 0:         // velocity prescribed
                    Ibvi03 = 0;
                    Ifadmii_3 = 0;
                    break;
                case 1:         // pressure prescribed
                    Ibvi03 = 1;
                    Ifadmii_3 = 0;
                    break;
                case 2:         // velocity and surface admittance prescribed
                    Ibvi03 = 0;
                    Ifadmii_3 = 1;
                    break;
                case 3:         // transfer admittance prescribed
                    Ibvi03 = 2;
                    Ifadmii_3 = 0;
                    break;
                case 4:         // transfer and surface admittances prescribed
                    Ibvi03 = 2;
                    Ifadmii_3 = 1;
                    break;
            }
            Zbvi03 = 0.0;
            for(i=0; i<listNumberNodesPerElement[iel]; i++) Zbvi03 += zbval0[iel][i];
            Zbvi03 /= (double)listNumberNodesPerElement[iel];
            
            // three cases
            if(listElementProperty[iel] == 0) // the source element is a surface element
            {
                if(Ifadmii_3) {Admia3 = zbval1[iel];} else {Admia3.set(0.0, 0.0);}
                NC_ComputeEntriesNforFMM(NCout, ic, iel, iadreclu);
            }
        } // end of loop IE
    } // end of loop IC
}

// compute the T-matrix and the T-vector for ML-FMM
void NC_BuildTforMLFMM
(
	ofstream& NCout
 )
{
    int i, j, k, nelj, irownu, ientry, iadrek = 0, ndofi, inod0 = 0, idfs, ncldofi, jdifdar;
    double *cencli;
    int ifvcj, nfmgp;
    
    int nrf_cl, rfa_cl, nlv, nPoSphe, nClustL, nrowT;
    double dfa_cl;
    bool mir_cl;
    Vector<bool> rfdi_cl(NDIM);
    
    ElCluster *clustArraLn;
    Complex *zTvcLn, *zTmxLn;
    int *irowTmxLn;
    
    // contribution of the element to the T-vector and the T-matrix
    Vector<Complex> TVele(clulevarry[0].nPoinSpheLv), TMele(clulevarry[0].nPoinSpheLv);
    
    // loop over cluster levels
    for(nlv=0; nlv<numClusterLevels_; nlv++)
    {
        // assign the parameters and arrays
        nPoSphe = clulevarry[nlv].nPoinSpheLv;
        nClustL = clulevarry[nlv].nClustSLv;
        nrowT = nPoSphe*nClustL;
        clustArraLn = clulevarry[nlv].ClastArLv;
        zTvcLn = tmtxlev[nlv].zTvcLv;
        zTmxLn = tmtxlev[nlv].zTmxLv;
        irowTmxLn = tmtxlev[nlv].irowTmxLv;
        
        // initialize the T-vector
        if(boolComputeTVector_) for(i=0; i<nrowT; i++) zTvcLn[i].set(0.0, 0.0);
        
        // initialize the T-matrix
        for(i=0; i<irowTmxLn[nrowT]; i++) zTmxLn[i].set(0.0, 0.0);
        
        // loop over clusters of the current level
        irownu = ientry = 0;
        for(i=0; i<nClustL; i++)
        {
            // cluster parameters
            idfs = clustArraLn[i].NDOFsPeEl;
            nrf_cl = clustArraLn[i].nuref;
            rfa_cl = clustArraLn[i].rffac;
            if(nrf_cl)
            {
                mir_cl = clustArraLn[i].ifmirro;
                for(j=0; j<NDIM; j++) rfdi_cl[j] = clustArraLn[i].ifrfdi[j];
            }
            
            // compute the number of gaussean points and the their local coordinates and  weights
            if(clustArraLn[i].IfMonoEl)
            {
                inod0 = listNumberNodesPerElement[clustArraLn[i].NumsOfEl[0]];
                nfmgp = NC_ComputeGausseanPoints(csiFMgau, etaFMgau, weiFMgau, inod0, GauFMBOrd3);
            }
            
            // number of the unknown DOFs of each element of the cluster and of the whole cluster
            ndofi = clustArraLn[i].NDOFsPeEl;
            ncldofi = clustArraLn[i].NumOfDOFs;
            // center of the cluster
            cencli = clustArraLn[i].CoorCent;
            
            // loop over elements of the cluster
            for(j=0; j<clustArraLn[i].NumOfEl; j++)
            {
                nelj = clustArraLn[i].NumsOfEl[j];
                
                // compute the number of gaussean points and the their local coordinates and weights
                if(!clustArraLn[i].IfMonoEl && listNumberNodesPerElement[nelj] != inod0)
                {
                    inod0 = listNumberNodesPerElement[nelj];
                    nfmgp = NC_ComputeGausseanPoints(csiFMgau, etaFMgau, weiFMgau, inod0, GauFMBOrd3);
                }
                
                // compute the element data
                NC_ComputeElementData(nelj, 0, Inoj3, Ibvj_03,
                                      Ifadmj3, Jdofaddr3, Centej3, Norvcj3, Admj3,
                                      Zbvj03, Crdej3);
                
                // reflection of the element
                if(nrf_cl) NC_ReflectElementFMM(rfa_cl, mir_cl, rfdi_cl, Inoj3, Ibvj_03,
                                                Centej3, Norvcj3, Zbvj03, Crdej3);
                
                // see if the velocity boundary condition of the current element should be considered
                ifvcj = 0;
                if(boolComputeTVector_ && (Ibvj_03 == 0 || Ibvj_03 == 1))
                {
                    for(k=0; k<Inoj3; k++) if(Zbvj03[k].norm() > 0.0)
                    {
                        ifvcj = 1;
                        break;
                    }
                }
                
                // initialize the output vectors
                for(k=0; k<nPoSphe; k++)
                {
                    TMele[k].set(0.0, 0.0);
                    TVele[k].set(0.0, 0.0);
                }
                
                if(Ifadmj3) {Admia3 = zbval1[nelj];} else {Admia3.set(0.0, 0.0);}
                NC_ComputeEntriesTforFMM(NCout, nelj, ifvcj, clulevarry[nlv].nPoinSpheLv, nfmgp,
                              cencli, clulevarry[nlv].uvcsphe, TMele, TVele);
                
                // assemble the contributions of the element to the T-matrix and the T-vector
                jdifdar = j*ndofi;
                if(rfa_cl == -1) {dfa_cl = -1.0;} else {dfa_cl = 1.0;}
                
                for(k=0; k<nPoSphe; k++)
                {
                    iadrek = ientry + k*ncldofi + jdifdar;
                    
                    zTmxLn[iadrek] = TMele[k]*dfa_cl;
                    if(idfs == 2)  zTmxLn[++iadrek] = TVele[k]*dfa_cl;
                    
                    if(ifvcj) zTvcLn[irownu + k] += TVele[k];
                }
            } // end of loop J
            irownu += nPoSphe;
            ientry += nPoSphe*ncldofi;
        } // end of loop I
        
        if(irowTmxLn[nrowT] != iadrek + 1)
        {
            NC_Error_Exit_2(NCout, "iadrek + 1 must be equal to irowTmxLn[nrowT]!",
                            "irowTmxLn[nrowT] = ", irowTmxLn[nrowT],
                            "iadrek = ", iadrek);
        }
    } // end of loop NLV
}

// compute the D-matrices for ML-FMM  (see eq. (18), (26) in [1])
void NC_BuildDforMLFMM
(
	ofstream& NCout
 )
{
    int i, j, jcl, k, jgp, identry, ikjadre;
    double dij, dkij, scprd;
    Complex zmun, zi;
    Vector<Complex> zH1m(clulevarry[0].nExpaTermLv), zdiat(clulevarry[0].nPoinSpheLv);
    Vector<double> Pm(clulevarry[0].nExpaTermLv), disij(NDIM);
    Complex *zwork_d;
    
    // loop over levels
    for(int nlv=0; nlv<numClusterLevels_; nlv++)
    {
        // array to be used to store the D-matrices
        zwork_d = dmtxlev[nlv].zDmxLv;
        
        // loop over the original clusters
        identry = 0; // address of first entry of the submatrices of the D-matrix that correspond to the cluster i
        for(i=0; i<clulevarry[nlv].nClustOLv; i++)
        {
            // compute and store the nonzeros corresponding to the current cluster
            for(j=0; j<clulevarry[nlv].ClastArLv[i].NumFanClus; j++)
            {
                jcl = clulevarry[nlv].ClastArLv[i].NumsFanClus[j];
                for(k=0; k<NDIM; k++) disij[k] = clulevarry[nlv].ClastArLv[i].CoorCent[k] -
                    clulevarry[nlv].ClastArLv[jcl].CoorCent[k];
                Vcnorm_dim3_(disij, dij);
                dkij = dij*waveNumbers_;
                
                // compute the Hankel functions
                NC_SphericalHankel(NCout, zH1m, clulevarry[nlv].nExpaTermLv, dkij);
                
                for(jgp=0; jgp<clulevarry[nlv].nPoinSpheLv; jgp++)
                {
                    scprd = Scprod_dim3_(disij, clulevarry[nlv].uvcsphe[jgp]);
                    
                    // compute the Legendre plynomes
                    NC_LegendrePolynomes(NCout, Pm, clulevarry[nlv].nExpaTermLv, scprd);
                    
                    // compute the function Mu_n
                    zmun.set(0.0, 0.0);
                    zi.set(1.0, 0.0);
                    for(k=0; k<clulevarry[nlv].nExpaTermLv; k++)
                    {
                        zmun += zi*zH1m[k]*Pm[k]*(double)(2*k + 1);
                        zi.mul_i(harmonicTimeFactor_); // _ruimag!
                    }
                    zdiat[jgp] = zmun;
                } // end of loop JGP
                
                // store the nonzeros into the D-matrix
                for(k=0; k<clulevarry[nlv].nPoinSpheLv; k++)
                {
                    ikjadre = identry + k*clulevarry[nlv].ClastArLv[i].NumFanClus + j;
                    zwork_d[ikjadre] = zdiat[k];
                }
            } // end of loop J
            identry += clulevarry[nlv].nPoinSpheLv*clulevarry[nlv].ClastArLv[i].NumFanClus;
        } // end of loop I
    } // end of loop NLV
}

// compute the S-matrices for ML-FMM  (see eq. (21), (26) in [1])
void NC_BuildSforMLFMM
(
	ofstream& NCout
)
{
	int icl, iel, nel, i, jip, iterm0, nlv;
	double scpr_dv, d0, wavruim = waveNumbers_*harmonicTimeFactor_;
	Complex z0;
	Vector<double> dcenel(NDIM);
	bool ifadm;
	Vector<Complex> zvcun(clulevarry[0].nPoinSpheLv), zvcdun(clulevarry[0].nPoinSpheLv);
	double *norel;

	int nPoSphe;
	ElCluster *clustArraLn;
	Complex *zSmxLn;
	int *irowSmxLn;
	double **uvshln;
	double *wshln;

	// loop over cluster levels
	for(nlv=0; nlv<numClusterLevels_; nlv++)
	{
		// assign the parameters and arrays
		nPoSphe = clulevarry[nlv].nPoinSpheLv;
		clustArraLn = clulevarry[nlv].ClastArLv;
		zSmxLn = smtxlev[nlv].zSmxLv;
		irowSmxLn = smtxlev[nlv].irowSmxLv;
		uvshln = clulevarry[nlv].uvcsphe;
		wshln = clulevarry[nlv].weisphe;

		// loop over clusters of the top level
		for(icl=0; icl<clulevarry[nlv].nClustOLv; icl++)
		{
			ifadm = clustArraLn[icl].IfAdmiBc;

			// loop over elements of the cluster
			for(iel=0; iel<clustArraLn[icl].NumOfEl; iel++)
			{
				nel = clustArraLn[icl].NumsOfEl[iel]; // number of the element

				// coordinate difference of the element center and the cluster center
				for(i=0; i<NDIM; i++) dcenel[i] = centel[nel][i] - clustArraLn[icl].CoorCent[i];

				// unit normal vector at the element center
				norel = elenor[nel];

				// loop over the integral points on the unit sphere
				for(jip=0; jip<nPoSphe; jip++)
				{
					// scalar product of two vectors
					scpr_dv = Scprod_dim3_(dcenel, uvshln[jip]);

					// compute the result vectors
					d0 = scpr_dv*wavruim;
					z0.set(cos(d0), sin(d0));
					z0 *= wshln[jip];
					zvcun[jip] = z0;

					d0 = Scprod_dim3_(uvshln[jip], norel);
					d0 *= wavruim;
					z0.mul_i(d0);
					zvcdun[jip] = z0;
				} // end of loop over JIP

				// store the results
				iterm0 = irowSmxLn[jelist[nel][0]];

                for(i=0; i<nPoSphe; i++)
                {zSmxLn[iterm0 + i] = zvcun[i]*Gama3 + zvcdun[i]*(zBta3*Tao_);}

			} // end of loop IEL
		} // end of loop ICL
	} // end of loop NLV
}

// compute the near field equation system for ML-FMM
void NC_BuildNforMLFMM
(
	ofstream& NCout
 )
{
    int i, ic, ie, iel, k, k1, nori_clu;
    Vector<bool> bneacluso(numOriginalClusters_);
    Vector<int> iadreclu(numOriginalClusters_);
    
    // initialize the near field coefficient matrix
    for(i=0; i<irownea[numRowsOfCoefficientMatrix_]; i++) zcoefl[i].set(0.0, 0.0);
    
    // loop over source clusters of the current level
    for(ic=0; ic<clulevarry[nlevtop_].nClustOLv; ic++)
    {
        // compute the array indicating if an original cluster is a near field cluster
        bneacluso = false;
        for(k=0; k<ClustArray[ic].NumNeaClus; k++)
        {
            nori_clu = ClustArray[ClustArray[ic].NumsNeaClus[k]].OriClust;
            if(!bneacluso[nori_clu]) bneacluso[nori_clu] = true;
        }
        
        // compute the adress of the first entry of each near field cluster
        iadreclu = -1;
        k1 = 0;
        for(k=0; k<numOriginalClusters_; k++)
        {
            if(!bneacluso[k]) continue;
            iadreclu[k] = k1;
            k1 += ClustArray[k].NumOfDOFs;
        }
        
        // loop over elements of the cluster
        for(ie=0; ie<ClustArray[ic].NumOfEl; ie++)
        {
            iel = ClustArray[ic].NumsOfEl[ie]; // number of the source element
            
            // center point and normal vector
            for(i=0; i<NDIM; i++)
            {
                Sourpoi3[i] = centel[iel][i];
                Norvci3[i] = elenor[iel][i];
            }
            
            // row numbers of the unknown DOFs of the source element
            Souradr3[0] = jelist[iel][0];
            Souradr3[1] = jelist[iel][1];
            
            // boundary condition
            switch(ibval[iel])
            {
                case 0:         // velocity prescribed
                    Ibvi03 = 0;
                    Ifadmii_3 = 0;
                    break;
                case 1:         // pressure prescribed
                    Ibvi03 = 1;
                    Ifadmii_3 = 0;
                    break;
                case 2:         // velocity and surface admittance prescribed
                    Ibvi03 = 0;
                    Ifadmii_3 = 1;
                    break;
                case 3:         // transfer admittance prescribed
                    Ibvi03 = 2;
                    Ifadmii_3 = 0;
                    break;
                case 4:         // velocity and transfer admittance prescribed
                    Ibvi03 = 2;
                    Ifadmii_3 = 1;
                    break;
            }
            Zbvi03 = 0.0;
            for(i=0; i<listNumberNodesPerElement[iel]; i++) Zbvi03 += zbval0[iel][i];
            Zbvi03 /= (double)listNumberNodesPerElement[iel];
            
            // three cases
            if(listElementProperty[iel] == 0) // the source element is a surface element
            {
                if(Ifadmii_3) {Admia3 = zbval1[iel];} else {Admia3.set(0.0, 0.0);}
                NC_ComputeEntriesNforFMM(NCout, ic, iel, iadreclu);
            }
        } // end of loop IE
    } // end of loop IC
}

// compute the contribution of the T-Vector to the right hand side vector (SLFMBEM)  (see eq. (25) in [1])
void NC_RHSvecTofSLFMM
(
	ofstream& NCout
)
{
	int i, j;
	Vector<Complex> zvc_d(dmtxlev[0].nRowsD);

	// D*T_vect
	for(i=0; i<dmtxlev[0].nRowsD; i++)
	{
		zvc_d[i].set(0.0, 0.0);
		for(j=dmtxlev[0].irowDmxLv[i]; j<dmtxlev[0].irowDmxLv[i + 1]; j++)
			zvc_d[i] += dmtxlev[0].zDmxLv[j]*ztvct[dmtxlev[0].jcolDmxLv[j]];
	}

	// ZRHS += S*D*T_vect
	for(i=0; i<numRowsOfCoefficientMatrix_; i++)
	{	
		for(j=irowsmtx[i]; j<irowsmtx[i + 1]; j++)
				zrhs[i] += zsmtx[j]*zvc_d[jcolsmtx[j]];
	}
}

// compute the contribution of the t-vector to the right hand side vector (DMLFMBEM)  (see eq. (26) in [1])
void NC_RHSvecTofMLFMM
(
	ofstream& NCout
)
{
	int i, j, ilv;
	Vector<Complex> zvc_d(maxRowNumberD_);

	// compute the S-vectors at each level S = D*T
	for(ilv=0; ilv<numClusterLevels_; ilv++)
	{
		for(i=0; i<dmtxlev[ilv].nRowsD; i++)
		{
			clulevarry[ilv].zwkS[i].set(0.0, 0.0);
			for(j=dmtxlev[ilv].irowDmxLv[i]; j<dmtxlev[ilv].irowDmxLv[i+1]; j++)
				clulevarry[ilv].zwkS[i] += 
				dmtxlev[ilv].zDmxLv[j]*tmtxlev[ilv].zTvcLv[dmtxlev[ilv].jcolDmxLv[j]];
		}

		for(i=0; i<numRowsOfCoefficientMatrix_; i++)
		{	
			for(j=smtxlev[ilv].irowSmxLv[i]; j<smtxlev[ilv].irowSmxLv[i+1]; j++) 
				zrhs[i] += smtxlev[ilv].zSmxLv[j]*clulevarry[ilv].zwkS[smtxlev[ilv].jcolSmxLv[j]];
		}
	} // end of loop ILV
}

// compute the contribution of a surface element to the T-matrix and the T-vector
/* kreiza: okay this is a bit tricky, because 
   a) Chen uses a negative sign for everything, so he looks at
      -(u+\beta v)/2 + \int (H + \beta E)u - \int (G + \beta H')v = -u_i - \beta v_i

   b)  G = (k I)/(4\pi) e_y M_L e_x, with e_y = exp(I k (z_1 - y) and e_x = exp(Ik (x - z_2)

   The T matrix covers everything with respect to y, thus the integral over \Gamma and the matrix parts
   For a velo b.c. we need 
         H + beta E as the matrix 
         G + beta H' as rhs vector
   thus
         ik/(4pi) (-ik (n_y.s) e_y M_L (1 + beta (n_x.s) e_x
	 ik/(4pi) z0(ey M_L (1 + beta ik (n_x.s))e_x)

   For a pres b.c. we need 
         -(G + beta H') as a matrix and 
         - (H + beta E) as rhs
   thus
        ik/(4pi) ( - e_y M_L (1 + beta ik n_x.s)e_x)
        ik/(4pi) (e_y (n_y.s) ik (1 + beta ik nx.s)e_x

   chen has written the code for VELO b.c. and forgot about the rest so the
   the I is never explicitely used because H and E have an I factor and the 
   kI/4pi has an I factor, thus we need only to multiply with -1, and additional
   -1 comes from the derivative with respect to y
   for the pres condition this does not work anymore, thus we currently multiply  with -I in the PRES part, which looks weird, but provides the right matrix
   
   for the vector Chen multiplies with ik/4pi thus we have to rethink again
 */
void NC_ComputeEntriesTforFMM
(
	ofstream& NCout,
	const int& nelj,	    // I: number of the element
	const int& ifvcj, 		// I: if the contribution to the T-vector should be computed
	int& nptshe,			// I: number of the integral points at the unit shere
	int& nfmgp,
	double* cencli,			// I: coordinates of the center of the cluster
	double** uvsph,			// I: coordinates of the integration points on the unit sphere
	Vector<Complex>& TMele, // O: contribution of the element to the T-matrix
	Vector<Complex>& TVele  // O: contribution of the element to the T-vector
)
{
    int i, inodj, igp, jgp;
    double vljaco, wga, dw1, wavruim = waveNumbers_*harmonicTimeFactor_;
    Vector<double> shfunx(NNODPE), shnds(NNODPE), shndt(NNODPE), dxds(NDIM),
    dxdt(NDIM), crdpoi(NDIM), elnorv(NDIM), vwk1(NDIM);
    Complex zbgao_0, zw1, zresu, zrev;
    
    
    inodj = listNumberNodesPerElement[nelj]; // number of nodes of the element
    
    // loop over the Gaussean points of the current element
    for(igp=0; igp<nfmgp; igp++)
    {
        
        // compute the parameters at the current Gaussean point
        vljaco = NC_ComputeParameters(shfunx, shnds, shndt, dxds, dxdt,
                                      crdpoi, elnorv, zbgao_0, Crdej3, ifvcj,
                                      Zbvj03, inodj, csiFMgau[igp], etaFMgau[igp]);
        wga = weiFMgau[igp]*vljaco;
        
        // loop over the Gaussean points of the unit sphere
        for(jgp=0; jgp<nptshe; jgp++)
        {
            for(i=0; i<NDIM; i++)
            {
                vwk1[i] = cencli[i] - crdpoi[i];
            }
            dw1 = (Scprod_dim3_(vwk1, uvsph[jgp]))*wavruim;
            zresu.set(cos(dw1), sin(dw1));
            if(ifvcj)
            {
                if(Ibvj_03 == 0)  {// veloc b.c.
                    zrev = zresu;
                    zrev *= (zbgao_0*wga);
                    TVele[jgp] += zrev;
                }
            }
            if(Ibvj_03 == 0) {// no pressure boundary condition
                dw1 = (Scprod_dim3_(elnorv, uvsph[jgp]))*waveNumbers_*Tao_;
                if(Admia3.norm() == 0.0) zresu *= (dw1*wga);
                else
                {
                    zw1.set(dw1 - harmonicTimeFactor_*Admia3.im(), harmonicTimeFactor_*Admia3.re());
                    zresu *= (zw1*wga);
                }
		// so zresu is now: k (n_y.s) e_y, so the -i factor is not used yet
            }
            if(Ibvj_03 == 1) {  // pressure b.c.
                zresu *= wga;
                /* set the tvector part, remember for the pressure
		   b.c. the rhs contribution is
		   -\int H - E, thus we need the derivative w.r.t
		   to n_y, the rest is done in the calculation of
		   S, the sign should be okay, because for the S
		   we also have a minus, and minus times minus is
		   plus
		   
		   below is a multiplication with
		   k/4pi, thus the i is missing, additionally we
		   the matrix is -(G+H')v, thus the negative sign

		*/
		zresu.mul_i(-1.0);

		
                if(ifvcj) {
                    dw1 = (Scprod_dim3_(elnorv, uvsph[jgp]))*waveNumbers_*Tao_;
                    TVele[jgp] -= zresu*zbgao_0*dw1;
                }
            }
            TMele[jgp] += zresu;
        } // end of loop JGP
    } // end of loop IGP
    
    dw1 = waveNumbers_/PI4;
    for(i=0; i<nptshe; i++) TMele[i] *= dw1;
    dw1 *= Tao_*harmonicTimeFactor_;
    if(ifvcj) for(i=0; i<nptshe; i++) TVele[i].mul_i(dw1);
}

// compute a near field equation if the source point is located in a surface element  (see eq. (13) [1])
void NC_ComputeEntriesNforFMM
(
	ofstream& NCout,
	const int& ic,           // number of the source cluster
	const int& iel,          // number of the source element
	Vector<int>& iadreclu    // address of the first entry of each near field cluster
)
{
    int jel, i, j, j1;
    Complex z0, Admib3(0.0, 0.0);
    
    int nu_cl, nrf_cl, rfa_cl, nuo_cl, jadr0_cl, jadrf_cl, jadre_cl;
    bool mir_cl;
    Vector<bool> rfdi_cl(NDIM);
    
    // compute of the free terms of the coefficient matrix
    NC_SearchAddressNforFMM(NCout, j1, Souradr3[0], Souradr3[0], 'o');
    switch(Ibvi03)
    {
        case 0:
	  /* that is slightly unnecessary if Admia3 is zero */
            zcoefl[j1] += (Admia3*zBta3 - Gama3)*0.5;
            break;
        case 1:
            zcoefl[j1] -= zBta3*(0.5*Tao_);
            break;
        default:
            NC_Error_Exit_2(NCout, "Ibvi03 must be equal to 0 or 1!", "Ibvi03 = ", Ibvi03,
                            "element number = ", extNumbersOfElements[iel]);
    }
    
    // compute of the free terms of the right hand side vector
    z0 = 0.0;
    // contributions of the velocity boundary condition
    if(Ibvi03 == 0) // velocity B. C. prescribed
    {
        z0 += Zbvi03*(zBta3*Tao_*0.5);
    }
    if(Ibvi03 == 1) // pres. B.C. prescribed
    {
        z0 += Zbvi03*Tao_*0.5;
    }
    // contributions of the incident waves
    if(numIncidentPlaneWaves_ + numPointSources_ > 0)
    {
        z0 += NC_IncidentWaveRHS(NCout);
    }
    zrhs[Souradr3[0]] += z0;
    
    // loop over the near clusters
    for(j=0; j<ClustArray[ic].NumNeaClus; j++)
    {
        // cluster parameters
        nu_cl = ClustArray[ic].NumsNeaClus[j];
        nuo_cl = ClustArray[nu_cl].OriClust;
        jadr0_cl = iadreclu[nuo_cl];
        jadrf_cl = ClustArray[nu_cl].NDOFsPeEl;
        nrf_cl = ClustArray[nu_cl].nuref;
        rfa_cl = ClustArray[nu_cl].rffac;
        if(nrf_cl)
        {
            mir_cl = ClustArray[nu_cl].ifmirro;
            for(i=0; i<NDIM; i++) rfdi_cl[i] = ClustArray[nu_cl].ifrfdi[i];
        }
        
        // loop over elements of the cluster
        for(j1=0; j1<ClustArray[nu_cl].NumOfEl; j1++)
        {
            // number of the field element
            jel = ClustArray[nu_cl].NumsOfEl[j1];
            
            // compute the element data
            NC_ComputeElementData(jel, 0, Inoj3, Ibvj_03,
                                  Ifadmj3, Jdofaddr3, Centej3, Norvcj3, Admj3,
                                  Zbvj03, Crdej3);
            
            // set the control parameter for computing the contributions to the r.h.s.
            Ifcrh3 = 0;
            if((Ibvj_03 == 0 || Ibvj_03 == 1) && listElementProperty[jel] == 0) for(i=0; i<Inoj3; i++)
                if(Zbvj03[i].norm() > EPSY) {Ifcrh3 = 1; break;}
            
            // reflection of the element
            if(nrf_cl) NC_ReflectElementFMM(rfa_cl, mir_cl, rfdi_cl, Inoj3, Ibvj_03,
                                            Centej3, Norvcj3, Zbvj03, Crdej3);
            
            // integration over the current field element
            if(nrf_cl == 0 && jel == iel) // singular integration
            {
                NC_SingularIntegration(NCout, Zrsintel, Inoj3, Ifcrh3, Zbvj03,
                                       Crdej3);
            }
            else // regular integration
            {
                NC_RegularIntegration(NCout, Zrsintel, iel, jel, Inoj3,
                                      Ifcrh3, Zbvj03, Centej3, Crdej3);
            }
            
            // assemble the contributions of the element to the equation system
            jadre_cl = jadr0_cl + j1*jadrf_cl;
            NC_AssembleNforFMM(NCout, jel, jadre_cl, Souradr3, rfa_cl);
        } // end of loop J1
    } // end of loop J
}

// assemble the contributions of an element to the near field coefficient matrix and right hand side vector
void NC_AssembleNforFMM
(
	ofstream& NCout,
	const int& jel,			 // I: number of the element
	const int& jadre_cl,     // adrress of the first entry of the current element
	Vector<int>& Souradr3,	 // I: address of the DOF of the source el.
	int& rfa_cl             // I: factor of reflection = 1 or -1
)
{
    Complex zsum0, zsum1;
    int iadr, kdadr;
    
    kdadr = jadre_cl;
    
    if(Ibvj_03 == 0) { // velo b.c.
        zsum0 = Zrsintel[1]*(Gama3*Tao_) + Zrsintel[3]*zBta3;
        if(Ifadmj3)
            zsum0 += (Zrsintel[0]*Gama3 + Zrsintel[2]*(zBta3*Tao_))*Admj3[0];
    }
    if(Ibvj_03 == 1) { // pres b.c.
        zsum0 = Zrsintel[0]*(-Gama3*Tao_) - Zrsintel[2]*(zBta3);
        if(Ifadmj3) {
            fprintf(stderr,"Sorry not implemented yet\n");
            exit(-1);
        }
    }
    
    iadr = irownea[Souradr3[0]] + kdadr;
    if(rfa_cl == -1) zcoefl[iadr] -= zsum0;
    else zcoefl[iadr] += zsum0;
    
    if(Ifcrh3) zrhs[Souradr3[0]] += Zrsintel[4];
}

// search address of an entry of the near fielf coefficient matrix
void NC_SearchAddressNforFMM
(
	ofstream& NCout,
	int& iadr,         //O: address
	int& nrow,         //I: row number
	int& ncol,         //I: column number
	const char& c	   //I: error output identifier
)
{
    int i;
    string strout("Address of an entry not found _ ");
    
    strout += c;
    
    iadr = -1;
    for(i=irownea[nrow]; i<irownea[nrow + 1]; i++) if(jcolnea[i] == ncol)
    {
        iadr = i;
        break;
    }
    if(iadr == -1) NC_Error_Exit_2(NCout, strout, "Row number = ", nrow,
                                   "Column number = ", ncol);
}

// function for performing the singular integration over an element (see eq. (14) in [1]), the idea is taken from the paper of Kirshnasamyetal where the hypersgl operator could be split into an integral over the edges \int_{edge} dG/dy_j dy_k times the normal vector and the Levi-Civita-symbol. Here the integral over the edge is one-dimensional because y_j = y_0 + t edge_j and everything can be written in terms of a cross product
void NC_SingularIntegration
(
	ofstream& NCout,
	Vector<Complex>& zrsintel, // O: result of the integration
 //    [0] integral of G
 //    [1] integral of H
 //    [2] integral of H(T)
 //    [3] integral of E
 //    [4] [5] contribution to RHS vector
	const int& inodj,          // I: number of nodes of the field el.
	const int& ifcrhs,         // I: flag for computing the r.s.h
	Vector<Complex>& zbvj_0,   // I: boundary condition for the under face
	Matrix<double>& crdelj    // I: nodal coordinates
)
{
    int ieg, i, ig1, ngpoi, ig, ngpo1 = 3, nsec1 = 4, nsec2 = 2, ig2, isec, ngausin = 4;
    double  sga, tga, wga, leneg, disfsp, re1, re2, sgg, tgg, vljaco, wnu2, secmid,
    delsec, wavruim= harmonicTimeFactor_*waveNumbers_;
    Vector<double> crdgp(NDIM), difpoi(NDIM), difpoo(NDIM), diffsp(NDIM);
    Complex zg, zhh, zht, z1, zre, zbgao_0, zbv;
    Vector<Complex> zdgdy(NDIM), zwk(NDIM);
    double ssub[3], tsub[3], aresub = 0.0, csigau[36], etagau[36], weigau[36];
    Vector<double> shfunx(NNODPE), shnds(NNODPE), shndt(NNODPE), dxds(NDIM),
    dxdt(NDIM), crdpoi(NDIM), elnorv(NDIM), coordgau(ngpo1), weitgau(ngpo1);
    
    // initialization
    for(i=0; i<6; i++) zrsintel[i].set(0.0, 0.0);
    wnu2 = waveNumbers_*waveNumbers_;
    
    // Gauss coordinates and weights for line integrations
    BLGauCooWe(NCout, ngpo1, coordgau, weitgau);
    
    // loop over edges of the element
    for(ieg=0; ieg<inodj; ieg++)
    {
        // set the nodal coordinates for the edge and subelement
        ig1 = ieg + 1;
        if(ig1 == inodj) ig1 = 0;
        ig2 = ieg + inodj;
        
        // perform the integration along the eadge
        leneg = 0.0;
        for(i=0; i<NDIM; i++)
        {
            difpoi[i] = crdelj(ig1, i) - crdelj(ieg, i);
            leneg += difpoi[i]*difpoi[i];
        }
        leneg = sqrt(leneg); // length of the edge
        for(i=0; i<NDIM; i++) difpoo[i] = difpoi[i]/leneg; // unit vector along the edge
        leneg /= (2.0*(double)nsec1);
        
        zre.set(0.0, 0.0);
        
        // each edge is divided into NSEC1 sections, loop over the sections
        delsec = 2.0/(double)nsec1;
        secmid = -1.0 - delsec/2.0;
        for(isec=0; isec<nsec1; isec++) {
            secmid += delsec;
            
            // loop over Gaussean points
            for(ig=0; ig<ngpo1; ig++)
            {
                sga = secmid + coordgau[ig]/(double)nsec1;
                wga = weitgau[ig]*leneg;
                for(i=0; i<NDIM; i++)
                {
                    crdgp[i] = crdelj(ieg, i) + difpoi[i]*(sga + 1.0)/2.0;
                    diffsp[i] = crdgp[i] - Sourpoi3[i];
                }
		// warning Vcnorm returns a vector of length 1
                Vcnorm_dim3_(diffsp, disfsp); // (y-x)/r and r
                re1 = wavruim*disfsp;
                re2 = PI4*disfsp;
                zg.set(cos(re1)/re2, sin(re1)/re2);
                
                z1.set(-1.0/disfsp, wavruim);
                zg.mul_c(z1);
                for(i=0; i<NDIM; i++) zdgdy[i] = zg*diffsp[i];
                zwk[0] = zdgdy[1]*difpoo[2] - zdgdy[2]*difpoo[1];
                zwk[1] = zdgdy[2]*difpoo[0] - zdgdy[0]*difpoo[2];
                zwk[2] = zdgdy[0]*difpoo[1] - zdgdy[1]*difpoo[0];
                zre += (zwk[0]*Norvci3[0] + zwk[1]*Norvci3[1] + zwk[2]*Norvci3[2])*wga;
                
            } // end of loop IG
        } // end of loop ISEC
        zrsintel[3] += zre;
        
        // each edge is divided into NSEC2 sections, loop over the sections
        for(isec=0; isec<nsec2; isec++) {
            
            // compute the area(/4) of the sub element in the local coordinate system and the local coordinates of nodes of the sub elements
            switch(inodj)
            {
	      // we construct 6 subtriangles, based on the midpoint of the element and the vertices + midpoints of edges
                case 3: // triangle element
                    aresub = 1.0/24.0/(double)nsec2;
                    ssub[0] = tsub[0] = 1.0/3.0;
                    if(isec == 0) {
		      // the difference is in the ordering/indices
                        ssub[1] = CSI6[ieg];
                        ssub[2] = CSI6[ig2];
                        tsub[1] = ETA6[ieg];
                        tsub[2] = ETA6[ig2];
                    }
                    else
                    {
                        ssub[1] = CSI6[ig2];
                        ssub[2] = CSI6[ig1];
                        tsub[1] = ETA6[ig2];
                        tsub[2] = ETA6[ig1];
                    }
                    break;
                case 4: // 4-sided element
                    aresub = 0.25/(double)nsec2;
                    ssub[0] = tsub[0] = 0.0;
                    if(isec == 0) {
                        ssub[1] = CSI8[ieg];
                        ssub[2] = CSI8[ig2];
                        tsub[1] = ETA8[ieg];
                        tsub[2] = ETA8[ig2];
                    }
                    else
                    {
                        ssub[1] = CSI8[ig2];
                        ssub[2] = CSI8[ig1];
                        tsub[1] = ETA8[ig2];
                        tsub[2] = ETA8[ig1];
                    }
                    break;
                default:
                    NC_Error_Exit_1(NCout, "inodj must be 3 or 4!", "inodj = ", inodj);
            }
            
            // compute the number of gaussean points, their local coordinates and weights
            ngpoi = NC_ComputeGausseanPoints(csigau, etagau, weigau, 4, ngausin);
            
            // loop over the gaussean points
            for(ig=0; ig<ngpoi; ig++)
            {
                sga = csigau[ig];
                tga = etagau[ig];
                
                // compute the local coordinates of the current integration point in the orginal element
                sgg = 0.5*(1.0 - sga)*ssub[0] +
                0.25*(1.0 + sga)*((1.0 - tga)*ssub[1] + (1.0 + tga)*ssub[2]);
                tgg = 0.5*(1.0 - sga)*tsub[0] +
                0.25*(1.0 + sga)*((1.0 - tga)*tsub[1] + (1.0 + tga)*tsub[2]);
                
                // compute the parameters at the current integration point
                vljaco = NC_ComputeParameters(shfunx, shnds, shndt, dxds, dxdt,
                                              crdpoi, elnorv, zbgao_0, crdelj, ifcrhs,
                                              zbvj_0, inodj, sgg, tgg);
                
                wga = weigau[ig]*(1.0 + sga)*aresub*vljaco;
                
                // compute the kernal functions
                for(i=0; i<NDIM; i++) diffsp[i] = crdpoi[i] - Sourpoi3[i];
                Vcnorm_dim3_(diffsp, disfsp); // "r" and the unit vector along "r"
                
                re1 = wavruim*disfsp;
                re2 = wga/PI4/disfsp;
                zg.set(cos(re1)*re2, sin(re1)*re2);
                
                z1.set(-1.0/disfsp, wavruim);
                zhh = zht = zg*z1;
                re1 = diffsp*elnorv;
                re2 = -(diffsp*Norvci3);
                zhh.mul_r(re1);
                zht.mul_r(re2);
                
                zrsintel[0] += zg;
                zrsintel[1] += zhh;
                zrsintel[2] += zht;
                
                zrsintel[3] += zg*(wnu2*Scprod_dim3_(Norvci3, elnorv));
                
                if(ifcrhs && Ibvj_03==0)
                {
                    zrsintel[4] += (zg*(Gama3*Tao_) + zht*zBta3)*zbgao_0;
                }
            } // end of loop IG
        } // end of loop ISEC
    } // end of loop IEG
    
    if(ifcrhs && Ibvj_03==1)
    {// pressure b.c.
        zrsintel[4] = ( zrsintel[1]*(Gama3*Tao_) + zrsintel[3]*zBta3 )*zbgao_0;
        zrsintel[4] *= -1;
    }
}

// function to perform the regular integration over an element
void NC_RegularIntegration
(
	ofstream& NCout,
	Vector<Complex>& zrsintel, // O: result of the integration
 //    [0] integral of G
 //    [1] integral of H
 //    [2] integral of H(T)
 //    [3] integral of E
 //    [4] [5] contribution to RHS vector
	const int& iel,            // I: number of the source element //x
	const int& jel,            // I: number of the field element
	const int& inodj,          // I: number of nodes of the field el.
	const int& ifcrhs,         // I: flag for computing the r.s.h
	Vector<Complex>& zbvj_0,   // I: boundary condition for the under face
	Vector<double>& centerj,   // I: center point of the field el.
	Matrix<double>& crdelj    // I: nodal coordinates
)
{
    int i, ngpoi, jgp;
    double vljaco, wga, disfsp, re1, re2, rq, dq, wavruim = harmonicTimeFactor_*waveNumbers_;
    double csigau[36], etagau[36], weigau[36];
    Vector<double> diffsp(NDIM);
    Complex zg, zhh, zht, ze, z1, zbgao_0, zbv;
    Vector<double> shfunx(NNODPE), shnds(NNODPE), shndt(NNODPE), dxds(NDIM),
    dxdt(NDIM), crdpoi(NDIM), elnorv(NDIM);
    
    // variables and arrays of subelements
    int nsbe, jsb, iforie;
    double fase, fase2, xio, eto, weih2, xice, etce;
    Vector<double> xisbp(MSBE), etsbp(MSBE), fasbp(MSBE), csip(NNODPE), etap(NNODPE);
    Vector<int> ngsbp(MSBE);
    
    // initialization
    for(i=0; i<6; i++) zrsintel[i].set(0.0, 0.0);
    
    // local coordinates of the nodes
    if(inodj == NETYP3)	for(i=0; i<inodj; i++)
    {
        csip[i] = CSI6[i];
        etap[i] = ETA6[i];
    }
    else for(i=0; i<inodj; i++)
    {
        csip[i] = CSI8[i];
        etap[i] = ETA8[i];
    }
    
    // generate the subelements
    nsbe = NC_GenerateSubelements(NCout, -1, Sourpoi3, jel, crdelj,
                                  areael[jel], inodj, csip, etap, inodj,
                                  xisbp, etsbp, fasbp, ngsbp); //jel field element, crdelj coordinates of field element, csi eta local parameters
    //nsbe ist anzahl der subelements
    
    // loop over the subelements
    for(jsb=0; jsb<nsbe; jsb++)
    {
        // assign values to some variables
        xice = xisbp[jsb]; // the first local coordinate //coordinates of subelement centers
        etce = etsbp[jsb]; // the second local cooedinate
        fase = fasbp[jsb]; // dimension factor //2^-j
        fase2 = fase*fase; // area factor
        
        if(fabs(fase - 1.0) < EPSY) iforie = 1; // the element is not subdivided
        else iforie = 0; // the element is subdivided
        
        // compute the number of gaussean points, their local coordinates and weights
        ngpoi = NC_ComputeGausseanPoints(csigau, etagau, weigau, inodj, ngsbp[jsb]); //lokale koordinaten vom subelement auf dem Einheitsdreieck
        
        // loop over the Gaussean points
        for(jgp=0; jgp<ngpoi; jgp++)
        {
            // compute the local coordinates and weight of the Gauss point in the original element
            if(iforie)
            {
                xio = csigau[jgp];
                eto = etagau[jgp];
                weih2 = weigau[jgp];
            }
            else
            {
                xio = xice + csigau[jgp]*fase;
                eto = etce + etagau[jgp]*fase;
                weih2 = weigau[jgp]*fase2;
            }
            
            // compute the parameters at the current Gaussean point
            vljaco = NC_ComputeParameters(shfunx, shnds, shndt, dxds, dxdt,
                                          crdpoi, elnorv, zbgao_0, crdelj, ifcrhs,
                                          zbvj_0, inodj, xio, eto);
            wga = weih2*vljaco;
            
            //G
            // compute the kernal functions //exp(ikr)/r/4/pi
            for(i=0; i<NDIM; i++) diffsp[i] = crdpoi[i] - Sourpoi3[i];
            Vcnorm_dim3_(diffsp, disfsp); // "r" and unit vector along "r"
            re1 = wavruim*disfsp;//ikr
            re2 = wga/PI4/disfsp;//gaußgewicht/4/pi/r
            zg.set(cos(re1)*re2, sin(re1)*re2);//exp(ikr)/r/4/pi
            
            //H
            z1.set(-1.0/disfsp, wavruim);
            zhh = zht = zg*z1;
            re1 = diffsp*elnorv;
            re2 = -(diffsp*Norvci3);
            zhh.mul_r(re1);
            zht.mul_r(re2);
            
            //E
            rq = re1*re2;
            re1 = Norvci3*elnorv; 
            dq = disfsp*disfsp;
            z1.set((3.0/dq - waveNumbers_*waveNumbers_)*rq + re1/dq, -wavruim/disfsp*(3.0*rq + re1));
            ze = zg*z1;
            
            zrsintel[0] += zg;
            zrsintel[1] += zhh;
            zrsintel[2] += zht;
            zrsintel[3] += ze;
            
            if(ifcrhs) 
            {
                if(Ibvj_03 == 0)  // velo b.c
                    zrsintel[4] += (zg*(Gama3*Tao_) + zht*zBta3)*zbgao_0;
                if(Ibvj_03 == 1)  // pres. b.c.
                    zrsintel[4] -= (zhh*(Gama3*Tao_) + ze*zBta3)*zbgao_0;
            } // end of IFCRHS
        } // end of loop JGP
    } // end of loop JSB
}

// compute the contributions of the incident waves to the right hand side vector
Complex NC_IncidentWaveRHS
(
	ofstream& NCout
)
{
    int i, j, kref;
    double re1, reffac, wavruim = harmonicTimeFactor_*waveNumbers_;
    Complex zrc(0.0, 0.0), zg, zv, zini;
    Vector<double> vect1(NDIM), vect2(NDIM);
    
    // compute the contribution of the plane waves
    if(numIncidentPlaneWaves_ > 0)
    {
        Vector<double> sourposym(Sourpoi3, NDIM);
        if(numSymmetricPlanes_ > 0)
        {
            for(i=0; i<NDIM; i++) if(listSymmetryPlanes[i] != 0) sourposym[i] -= listSymmetryPlanesCoordinates[i];
        }
        
        // loop over incident plane waves
        for(i=0; i<numIncidentPlaneWaves_; i++)
        {
            // direction of the cuurent incident plane wave
            for(j=0; j<NDIM; j++) vect1[j] = dirinw[i][j];
            
            // reflection of the incident plane wave with respect to the symmetry/antisymmetry planes
            zini = zinwav[i];
            reffac = 1.0;
            for(kref=0; kref<numReflectionsOfElements_; kref++)
            {
                if(kref > 0) reffac = NC_ReflectPointVector(vect1, kref, reffac, 0);
                
                re1 = vect1*sourposym*wavruim;
                zg.set(cos(re1)*reffac, sin(re1)*reffac);
                zg.mul_c(zini);
                
                zv = zg;
                re1 = wavruim*(vect1*Norvci3);
                zv.mul_i(re1);
                
                zrc -= zg*Gama3 + zv*(zBta3*Tao_);
            }
        } // end of loop I
    }
    
    // compute the contribution of the point sources
    if(numPointSources_ > 0)
    {
        double aglen, re2;
        
        // loop over the point sources
        for(i=0; i<numPointSources_; i++)
        {
            // coordinates of the point source
            for(j=0; j<NDIM; j++) vect1[j] = coorps[i][j];
            
            // reflection of the point source with respect to the symmetry/antisymmetry planes
            zini = zposoc[i];
            reffac = 1.0;
            for(kref=0; kref<numReflectionsOfElements_; kref++)
            {
                if(kref > 0) reffac = NC_ReflectPointVector(vect1, kref, reffac, 1);
                
                for(j=0; j<NDIM; j++) vect2[j] = Sourpoi3[j] - vect1[j];
                Vcnorm_dim3_(vect2, aglen);
                
                re1 = wavruim*aglen;
                re2 = reffac/PI4/aglen;
                zg.set(cos(re1)*re2, sin(re1)*re2);
                zg.mul_c(zini);
                
                re1 = vect2*Norvci3;
                zv.set(-re1/aglen, wavruim*re1);
                zv *= zg;
                
                zrc -= zg*Gama3 + zv*(zBta3*Tao_);
            }
        } // end of loop I
    }
    
    return zrc;
}

// update the admittance boundary conditions and the incident waves by using the curves describing their dependance on the frequency
void NC_UpdateFreqCurves
(
	ofstream& NCout,
	double* Freqs
 )
{
    int i, j, k, ibg = 0, ied = 0, ibc;
    double facf, gammfact = waveNumbers_*harmonicTimeFactor_;
    Complex zv;
    
    // loop over elements
    for(i=0; i<numElements_; i++)
    {
        if(listElementProperty[i] == 2) continue;
        
        // copy the boundary conditions
        for(j=0; j<listNumberNodesPerElement[i]; j++) zbval0[i][j] = zbvao0[i][j];
        zbval1[i] = zbvao1[i];
        zbval2[i] = zbvao2[i];
        
        // interpolation of the admittance boundary conditions
        if(ibval[i] < 2) continue;
        switch(ibval[i])
        {
            case 2:       // surface admittances are prescribed
                ibg = 2;
                ied = 5;
                break;
            case 3:       // transfer admittance is prescribed
                ibg = 0;
                ied = 1;
                break;
            case 4:       // transfer and surface admittances are prescribed
                ibg = 0;
                ied = 5;
                break;
            default:
                NC_Error_Exit_0(NCout, "ibval[i] must be not greater than 4!");
        }
        for(j=ibg; j<=ied; j++)
        {
            ibc = ibvcurv[i][j];
            if(ibc < 0) continue;
            facf = NC_FreqCurveInterpolation(NCout, ibc);
            switch(j)
            {
                case 0:
                    for(k=0; k<listNumberNodesPerElement[i]; k++)
                    {
                        zbval0[i][k].set(zbval0[i][k].re()*facf, zbval0[i][k].im());
                    }
                    break;
                case 1:
                    for(k=0; k<listNumberNodesPerElement[i]; k++)
                    {
                        zbval0[i][k].set(zbval0[i][k].re(), zbval0[i][k].im()*facf);
                    }
                    break;
                case 2:
                    zbval1[i].set(zbval1[i].re()*facf, zbval1[i].im());
                    break;
                case 3:
                    zbval1[i].set(zbval1[i].re(), zbval1[i].im()*facf);
                    break;
                case 4:
                    zbval2[i].set(zbval2[i].re()*facf, zbval2[i].im());
                    break;
                case 5:
                    zbval2[i].set(zbval2[i].re(), zbval2[i].im()*facf);
                    break;
            }
        } // end of loop J
    } // end of loop I
    
    // modyfy the admittance, transfer admittance and pressure boundary conditions because the velocity potential, instead of the sound pressure, is used as the unknown variable
    for(i=0; i<numElements_; i++)
    {
        if(listElementProperty[i] == 2) continue;
        switch(ibval[i])
        {
            case 1:  // for pressure boundary condition
                for(j=0; j<listNumberNodesPerElement[i]; j++)
                {
                    zbval0[i][j].div_i(rpfact_);
                }
                break;
            case 2: // for surface (normalized) admittance boundary condition
                zbval1[i].mul_i(gammfact);
                break;
            case 3: // for (normalized) transfer admittance
                for(j=0; j<listNumberNodesPerElement[i]; j++)
                {
                    zbval0[i][j].mul_i(gammfact);
                }
                break;
            case 4: // for (normalized) transfer admittance and (normalized) surface admittance
                for(j=0; j<listNumberNodesPerElement[i]; j++)
                {
                    zbval0[i][j].mul_i(gammfact);
                }
                zbval1[i].mul_i(gammfact);
                break;
        }
    } // end of loop I
    
    // update the incident plane waves
    if(numIncidentPlaneWaves_ > 0)
    {
        for(i=0; i<numIncidentPlaneWaves_; i++)
        {
            zv = zinwao[i];
            for(j=0; j<2; j++)
            {
                ibc = inwacurv[i][j];
                if(ibc >= 0)
                {
                    facf = NC_FreqCurveInterpolation(NCout, ibc);
                    if(j == 0) {zv.set(zv.re()*facf, zv.im());}
                    else {zv.set(zv.re(), zv.im()*facf);}
                }
            }
            zv.div_i(rpfact_);
            zinwav[i] = zv;
        }
    } // end of IF
    
    // update the point sources
    if(numPointSources_ > 0)
    {
        for(i=0; i<numPointSources_; i++)
        {
            zv = zposoo[i];
            for(j=0; j<2; j++)
            {
                ibc = posocurv[i][j];
                if(ibc >= 0)
                {
                    facf = NC_FreqCurveInterpolation(NCout, ibc);
                    if(j == 0) {zv.set(zv.re()*facf, zv.im());}
                    else {zv.set(zv.re(), zv.im()*facf);}
                }
            }
            zv.div_i(rpfact_);
            zposoc[i] = zv;
        }
    } // end of IF
}

// perform the interpolation along a curve
double NC_FreqCurveInterpolation
(
	ofstream& NCout, 
	int& nucv
 )
{
    double factor;
    int l, l1 = 0, l2, l3;
    
    // seek the relevant curve
    for(l=0; l<numCurvesFrequency_; l++)
    {
        if(nucv == nucurv[l]) 
        {
            l1 = l;
            goto Interpo_Curve;
        }
    }
    NC_Error_Exit_1(NCout, "A curve can not be found!", "Number of the curve = ", nucv);
    
Interpo_Curve:
    // perform the interpolation according to the curve
    if(frequency_ <= frqcurv[l1][0])
    {
        factor = faccurv[l1][0];
        return factor;
    }
    
    for(l=1; l<npcurv[l1]; l++)
    {
        if(frequency_ < frqcurv[l1][l])
        {
            l2 = l;
            goto Interpo_Point;
        }
    }
    factor = faccurv[l1][npcurv[l1] - 1];
    return factor;
    
Interpo_Point:
    l3 = l2 - 1;
    factor = faccurv[l1][l3] + (faccurv[l1][l2] - faccurv[l1][l3])/
    (frqcurv[l1][l2] - frqcurv[l1][l3])*(frequency_ - frqcurv[l1][l3]);
    
    return factor;
}
