/*===========================================================================*\
 *                                                                            *
 *  File name:      NC_PostProcessing.cpp                                     *
 *  Description:    computing and writing the results of the job              *
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
#ifdef isWindows                                                                    //
#include<direct.h> // for mkdir WINDOWS                                             //
#else                                                                               //
#include<sys/stat.h> // for mkdir UNIX                                              //
#endif                                                                              //
using namespace std;                                                                //
//====================================================================================



// local functions
void NC_ComputeTvecSLFMM(ofstream&);
void NC_ComputeTvecMLFMM(ofstream&);
void NC_ContributionTvecFMM(ofstream&,const int&,int&,Vector<double>&,Matrix<double>&,Complex&,Vector<Complex>&, double**,Vector<Complex>&);
void NC_WriteResultsObjectMesh(ofstream&,Matrix<Complex>&,Matrix<Complex>&,Matrix<Complex>&,Matrix<double>&);
void NC_WriteResultsEvaluationGrid(ofstream&, Matrix<Complex>&, Matrix<Complex>&);
Complex NC_ContributionIncidentWaves(ofstream&,Vector<double>&,Vector<Complex>&,const int&);
void NC_ContributionTBEM(ofstream&,Matrix<Complex>&,Matrix<Complex>&,int&,int&,Vector<double>&,Vector<Complex>&,Matrix<Complex>&);
void NC_ContributionFMM(ofstream&,Matrix<Complex>&,Matrix<Complex>&,Vector<int>&,Vector<Complex>&,Matrix<Complex>&);
Complex NC_IntegrationTBEM(ofstream&,Vector<double>&,const double&,const double&,const int&,const double&,const int&,Matrix<double>,const Complex,const int&,Vector<Complex>&,Vector<Complex>&,const int&, const int&);
void NC_Magnitude2dBdeg(Vector<double>&, const Complex&, const int&);



// cantant values
const double PREREF = 2.0e-5,
			 VELREF = 5.0e-8;



// files
FILE *Lu_pBoundary, *Lu_vBoundary, *Lu_pBoundaryVTK;



// global variables
char Dfname_3d[30];
Vector<double> Vmagdbph_pp_3d(3);
Vector<double> Vmagdbph_v_3d(3);
Vector<Complex> Zvelopoi_3d(NDIM);
Complex Zpotpoi_3d;
const int GaupFMBOrde = 2;
int npFMGp;
double csipFMgau[GaupFMBOrde*GaupFMBOrde], etapFMgau[GaupFMBOrde*GaupFMBOrde], weipFMgau[GaupFMBOrde*GaupFMBOrde];
Complex *zT_vc;
static bool IfTvc = false;



// the main program for postprecessing
void NC_PostProcessing
(
	ofstream& NCout
)
{
	Matrix<Complex> cVpotele(numElements_, 2), cVeloele(numElements_, NNODPE*2), cVeloelo(numElements_, NNODPE);
	Matrix<double> rEninele(numElements_, NNODPE*2);

	// create the BE output directories
	sprintf(Dfname_3d, "be.out/be.%d",  currentFrequency_ + 1);
#ifdef isWindows
	int ifmkd;
    ifmkd = _mkdir(Dfname_3d); // WINDOWS
#else
	mkdir(Dfname_3d, S_IRUSR | S_IWUSR | S_IXUSR | S_IRGRP | S_IWGRP | S_IXGRP | S_IROTH | S_IXOTH);  // UNIX
#endif

	// create the FE output directories
	//	sprintf(Dfname_3d, "fe.out/fe.%d",  currentFrequency_ + 1);
#ifdef isWindows
	//    ifmkd = _mkdir(Dfname_3d); // WINDOWS
#else
	//mkdir(Dfname_3d, S_IRUSR | S_IWUSR | S_IXUSR | S_IRGRP | S_IWGRP | S_IXGRP | S_IROTH | S_IXOTH);  // UNIX
#endif

	// compute the bool variable indicating if the T-vector should be computed
	if(numInternalPointsClusters_) IfTvc = true;

	// compute the T-vector
	if(IfTvc)
	{
		switch(methodFMM_)
		{
		case 1: // SLFMBEM
			zT_vc = new Complex[numOriginalReflectedClusters_*numIntegrationPointsUnitSphere_];
			NC_ComputeTvecSLFMM(NCout);
			break;
		case 3: // DMLFMBEM
			NC_ComputeTvecMLFMM(NCout);
			break;
		}
	}

	// compute and output the results at the boundary elements
	NC_WriteResultsObjectMesh(NCout, cVpotele, cVeloele, cVeloelo, rEninele);

	// compute and output the results at nodes of the evaluation mesh
	if(numNodesOfEvaluationMesh_ > 0)
    {
        NC_WriteResultsEvaluationGrid(NCout, cVpotele, cVeloele);
    }

	// destroy the T-vector
	if(IfTvc && (methodFMM_ == 1 || methodFMM_ == 2)) delete [] zT_vc;
}

// compute the T-vector for the SL-FMM
void NC_ComputeTvecSLFMM
(
	ofstream& NCout
 )
{
    int i, j, inod0 = 0, k, ibg, jel, iveloj, ifadmij, inodj;
    Vector<double> center_i(NDIM);
    Vector<int> jdofaddre(2);
    Vector<double> centerj(NDIM), norvecj(NDIM*(NNODPE + 1));
    Matrix<double> crdelj(NNODPE, NDIM);
    Vector<Complex> zpotej(2);
    Vector<Complex> zveloj(NNODPE);
    
    int nrf_cl, rfa_cl;
    bool mir_cl;
    Vector<bool> rfdi_cl(NDIM);
    
    
    // contribution of the element to the T-vector
    Vector<Complex> zTele(numIntegrationPointsUnitSphere_);
    
    // compute the T-vector
    for(i=0; i<numIntegrationPointsUnitSphere_*numOriginalReflectedClusters_; i++) zT_vc[i].set(0.0, 0.0);
    ibg = 0;
    for(i=0; i<numOriginalReflectedClusters_; i++)
    {
        // cluster parameters
        nrf_cl = ClustArray[i].nuref;
        rfa_cl = ClustArray[i].rffac;
        if(nrf_cl)
        {
            mir_cl = ClustArray[i].ifmirro;
            for(j=0; j<NDIM; j++) rfdi_cl[j] = ClustArray[i].ifrfdi[j];
        }
        
        // compute the number of gaussean points and their local coordinates and weights
        if(ClustArray[i].IfMonoEl)
        {
            inod0 = listNumberNodesPerElement[ClustArray[i].NumsOfEl[0]];
            npFMGp = NC_ComputeGausseanPoints(csipFMgau, etapFMgau, weipFMgau, inod0, GaupFMBOrde);
        }
        
        // center of the cluster
        for(k=0; k<NDIM; k++) center_i[k] = ClustArray[i].CoorCent[k];
        
        // loop over elements of the cluster
        for(j=0; j<ClustArray[i].NumOfEl; j++)
        {
            jel = ClustArray[i].NumsOfEl[j];
            
            // compute the number of gaussean points and their local coordinates and weights
            if(!ClustArray[i].IfMonoEl && listNumberNodesPerElement[jel] != inod0)
            {
                inod0 = listNumberNodesPerElement[jel];
                npFMGp = NC_ComputeGausseanPoints(csipFMgau, etapFMgau, weipFMgau, inod0, GaupFMBOrde);
            }
            
            // compute the element data
            NC_ComputeElementData(jel, 1, inodj, iveloj,
                                  ifadmij, jdofaddre, centerj, norvecj, zpotej,
                                  zveloj, crdelj);
            
            // reflection of the element
            if(nrf_cl) NC_ReflectElementFMM(rfa_cl, mir_cl, rfdi_cl, inodj, 0,
                                            centerj, norvecj, zveloj, crdelj);
            if(nrf_cl && rfa_cl == -1) zpotej[0].nega();
            
            // compute the contribution of the element to the T-vectors
            for(k=0; k<numIntegrationPointsUnitSphere_; k++) zTele[k].set(0.0, 0.0);
            NC_ContributionTvecFMM(NCout, jel, numIntegrationPointsUnitSphere_, center_i, crdelj, zpotej[0],
                                   zveloj, uvcsphe, zTele);
            
            for(k=0; k<numIntegrationPointsUnitSphere_; k++) zT_vc[ibg + k] += zTele[k];
        } // end of loop J
        ibg += numIntegrationPointsUnitSphere_;
    } // end of loop I
}

// function to compute the T-vectors for the ML-FMM
void NC_ComputeTvecMLFMM
(
	ofstream& NCout
 )
{
    int i, j, inod0 = 0, k, ibg, jel, iveloj, ifadmij, inodj, ilv;
    Vector<double> center_i(NDIM);
    Vector<int> jdofaddre(2);
    Vector<double> centerj(NDIM), norvecj(NDIM*(NNODPE + 1));
    Matrix<double> crdelj(NNODPE, NDIM);
    Vector<Complex> zpotej(2);
    Vector<Complex> zveloj(NNODPE);
    
    int nrf_cl, rfa_cl;
    bool mir_cl;
    Vector<bool> rfdi_cl(NDIM);
    
    ElCluster *clusArray;
    
    
    // contribution of the element to the T-vector
    Vector<Complex> zTele(clulevarry[0].nPoinSpheLv);
    
    // compute the T-vectors
    for(ilv=0; ilv<numClusterLevels_; ilv++)
    {
        // defines
        clusArray = clulevarry[ilv].ClastArLv;
        
        // initialize
        for(i=0; i<clulevarry[ilv].nPoinSpheLv*clulevarry[ilv].nClustSLv; i++)
            clulevarry[ilv].zwkT[i].set(0.0, 0.0);
        
        ibg = 0;
        for(i=0; i<clulevarry[ilv].nClustSLv; i++)
        {
            // cluster parameters
            nrf_cl = clusArray[i].nuref;
            rfa_cl = clusArray[i].rffac;
            if(nrf_cl)
            {
                mir_cl = clusArray[i].ifmirro;
                for(j=0; j<NDIM; j++) rfdi_cl[j] = clusArray[i].ifrfdi[j];
            }
            
            // compute the number of gaussean points and their local coordinates and weights
            if(clusArray[i].IfMonoEl)
            {
                inod0 = listNumberNodesPerElement[clusArray[i].NumsOfEl[0]];
                npFMGp = NC_ComputeGausseanPoints(csipFMgau, etapFMgau, weipFMgau, inod0, GaupFMBOrde);
            }
            
            // center of the cluster
            for(k=0; k<NDIM; k++) center_i[k] = clusArray[i].CoorCent[k];
            
            // loop over elements of the cluster
            for(j=0; j<clusArray[i].NumOfEl; j++)
            {
                jel = clusArray[i].NumsOfEl[j];
                
                // compute the number of gaussean points and their local coordinates and  weights
                if(!clusArray[i].IfMonoEl && listNumberNodesPerElement[jel] != inod0)
                {
                    inod0 = listNumberNodesPerElement[jel];
                    npFMGp = NC_ComputeGausseanPoints(csipFMgau, etapFMgau, weipFMgau, inod0, GaupFMBOrde);
                }
                
                // compute the element data 
                NC_ComputeElementData(jel, 1, inodj, iveloj,
                                      ifadmij, jdofaddre, centerj, norvecj, zpotej,
                                      zveloj, crdelj);
                
                // reflection of the element
                if(nrf_cl) NC_ReflectElementFMM(rfa_cl, mir_cl, rfdi_cl, inodj, 0,
                                                centerj, norvecj, zveloj, crdelj);
                if(nrf_cl && rfa_cl == -1) zpotej[0].nega();
                
                // compute the contribution of the element to the T-vectors
                for(k=0; k<clulevarry[ilv].nPoinSpheLv; k++) zTele[k].set(0.0, 0.0);
                NC_ContributionTvecFMM(NCout, jel, clulevarry[ilv].nPoinSpheLv, center_i, crdelj,
                                       zpotej[0], zveloj, clulevarry[ilv].uvcsphe, zTele);
                
                for(k=0; k<clulevarry[ilv].nPoinSpheLv; k++) 
                    clulevarry[ilv].zwkT[ibg + k] += zTele[k];
            } // end of loop J
            ibg += clulevarry[ilv].nPoinSpheLv;
        } // end of loop I
    } // end of loop ILV
}

// compute the contribution of an element to the T-vectors
void NC_ContributionTvecFMM
(
	ofstream& NCout,
	const int& nelj,		 // I: number of the element
	int& npshe,				 // I: number of points on the unit sphere
	Vector<double>& cencli,  // I: coordinates of the center of the cluster
	Matrix<double>& crdelj,  // I: corrdinates of nodes of the element
	Complex& zpote0,
	Vector<Complex>& zveloj, // I: nodal velocities or difference of nodal velocities
	double **uvshp,			 // I: coordinates of the integration points on the unit sphere
	Vector<Complex>& zTele   // O: Contribution of the element to the T-vector
)
{
    int i, inodj, igp, jgp;
    double vljaco, wga, dw1;
    Vector<double> shfunx(NNODPE), shnds(NNODPE), shndt(NNODPE), dxds(NDIM),
    dxdt(NDIM), crdpoi(NDIM), elnorv(NDIM), vwk1(NDIM);
    Complex zbgao_0, zw1, zG, zH;
    
    
    inodj = listNumberNodesPerElement[nelj]; // number of nodes of the element
    
    // loop over the Gaussean points of the current element
    for(igp=0; igp<npFMGp; igp++)
    {
        
        // compute the parameters at the current Gaussean point
        vljaco = NC_ComputeParameters(shfunx, shnds, shndt, dxds, dxdt,
                                      crdpoi, elnorv, zbgao_0, crdelj, 1,
                                      zveloj, inodj, csipFMgau[igp], etapFMgau[igp]);
        wga = weipFMgau[igp]*vljaco;
        
        // loop over the Gaussean points of the unit sphere
        for(jgp=0; jgp<npshe; jgp++)
        {
            vwk1 = cencli - crdpoi;
            dw1 = (Scprod_dim3_(vwk1, uvshp[jgp]))*waveNumbers_*harmonicTimeFactor_;
            zG.set(cos(dw1), sin(dw1));
            
            dw1 = (Scprod_dim3_(elnorv, uvshp[jgp]))*waveNumbers_;
            zH = zG*(dw1*wga);
            zG.mul_i(wga*harmonicTimeFactor_);

	    /* kreiza
	       here we have 2 problems:
	       a) Chen forgot the PRES version
	       b) BLeldat gives either the velocity potential or the particle velocity but not both
	       and I don't like to change that, because it might be needed somewhere for the admittance
	    */
	    if(ibval[nelj] == 1) {// PRES  kreiza that may be a bit tricky
	      zbgao_0 = zbvao0[nelj][0];
	      zbgao_0.div_i(rpfact_);
	      zTele[jgp] += zH*zbgao_0 - zG*zveloj[0];  
	    }
	    else
	      zTele[jgp] += zH*zpote0 - zG*zbgao_0;
        } // end of loop JGP
    } // end of loop IGP
    
    dw1 = waveNumbers_/PI4;
    for(i=0; i<npshe; i++) zTele[i].mul_r(dw1);
    if(!listElementProperty[nelj] && Tao_ != 1.0) for(i=0; i<npshe; i++) zTele[i].mul_r(Tao_);
}

// compute and write the results at the boundary
void NC_WriteResultsObjectMesh 
(
	ofstream& NCout,
	Matrix<Complex>& cVpotele,
	Matrix<Complex>& cVeloele,
	Matrix<Complex>& cVeloelo,  
	Matrix<double>& rEninele
)
{
	int i, j, k, k1, iel, igr;
	int ibvi_0, ifadmii, inodi;
	double re1;
	Complex z0, z1, zpotsum;
	Vector<int> ndgrp(numElementGroups_), thingrp(numElementGroups_), mdgrp(numNodes_), nundgrp(numNodes_);
	Vector<int> idofaddre(2);
	Vector<double> centeri(NDIM), norveci(NDIM*(NNODPE + 1)), wkpoi(NDIM);
	Vector<Complex> admi(2);
	Vector<Complex> zbvi_0(NNODPE);
	Matrix<double> crdeli(NNODPE, NDIM);
	/* fe.out will be removed
	   FILE *lu_load;
	*/

	// loop over the elements to compute the pressures, velocities and energy intensities in each element
	for(iel=0; iel<numElements_; iel++)
	{
		if(listElementProperty[iel] == 2) continue;

		// compute the element data
		NC_ComputeElementData(iel, 0, inodi, ibvi_0,
            ifadmii, idofaddre, centeri, norveci, admi,
            zbvi_0, crdeli);

        switch(ibvi_0)
        {
        case 0: // velocity boundary condition
            cVpotele(iel, 0) = zrhs[idofaddre[0]];
            for(i=0; i<inodi; i++) cVeloelo(iel, i) = cVeloele(iel, i) = zbvi_0[i];
            if(ifadmii) 
            {
                z1 = (cVpotele(iel, 0)*admi[0])*Tao_;
                for(i=0; i<inodi; i++) cVeloele(iel, i) -= z1;
            }
            break;
        case 1: // pressure boundary condition
            cVpotele(iel, 0) = zbvi_0[0];
            for(i=0; i<inodi; i++) cVeloelo(iel, i) = cVeloele(iel, i) =
                zrhs[idofaddre[0]];
            break;
        case 2: // transfer admittance boundary condition
            NC_Error_Exit_0(NCout, "transfer admittance boundary condition is not supported!");
        }

        for(i=0; i<inodi; i++) rEninele(iel, i) = 
            (cVpotele(iel, 0).re()*cVeloele(iel, i).im() - 
             cVpotele(iel, 0).im()*cVeloele(iel, i).re())*rpfact_*Tao_;
	} // end if loop IEL

	// open the files
    sprintf(Dfname_3d, "be.out/be.%d/pBoundary",  currentFrequency_ + 1);
	Lu_pBoundary = fopen(Dfname_3d, "w");
    if(!Lu_pBoundary) NC_Error_Exit_1(NCout,  "Can not open the file ", Dfname_3d);
    sprintf(Dfname_3d, "be.out/be.%d/vBoundary",  currentFrequency_ + 1);
    Lu_vBoundary = fopen(Dfname_3d, "w");
    if(!Lu_vBoundary) NC_Error_Exit_1(NCout,  "Can not open the file ", Dfname_3d);
    
    // open the VTK file
//    sprintf(Dfname_3d, "be.out/be.%d/pBoundaryVTK.vtk",  currentFrequency_ + 1);
//    Lu_pBoundaryVTK = fopen(Dfname_3d, "w");
//    if(!Lu_pBoundaryVTK) NC_Error_Exit_1(NCout,  "Can not open the file ", Dfname_3d);
//    fprintf(Lu_pBoundaryVTK, "# vtk DataFile Version 3.0\n");
//    fprintf(Lu_pBoundaryVTK, "%s\n", jobTitle_);
//    fprintf(Lu_pBoundaryVTK, "ASCII\n");
//    fprintf(Lu_pBoundaryVTK, "DATASET POLYDATA\n");

    /* fe.out removed kreiza
	sprintf(Dfname_3d, "fe.out/fe.%d/load",  currentFrequency_ + 1);
	lu_load = fopen(Dfname_3d, "w");
        if(!lu_load) NC_Error_Exit_1(NCout,  "Can not open the file ", Dfname_3d);
    */
    fprintf(Lu_pBoundary, "%s\n", versionNumber_.c_str());
    fprintf(Lu_vBoundary, "%s\n", versionNumber_.c_str());
    /* fe.out stuff removed for now kreiza
	fprintf(lu_load, "%s\n", versionNumber_);
	fprintf(lu_load, "1 0.0 %E %d\n", frequency_, currentFrequency_ + 1);
	fprintf(lu_load, "%E\n", frequency_);
	fclose(lu_load);
    */

	// compute the vector NDGRP (number of nodes of each element group) and THINGRP
	for(igr=0; igr<numElementGroups_; igr++)
	{
		thingrp[igr] = -1;
		for(j=0; j<numNodes_; j++) mdgrp[j] = -1;
		for(iel=0; iel<numElements_; iel++)
		{
			if(listElementProperty[iel] == 2 || indexOfElementGroup[igr] != listElementsElementGroup[iel]) continue;
			if(thingrp[igr] == -1) thingrp[igr] = listElementProperty[iel];
			for(j=0; j<listNumberNodesPerElement[iel]; j++) mdgrp[elementsConnectivity[iel][j]] = 1;			
		}
		k1 = 0;
		for(j=0; j<numNodes_; j++) if(mdgrp[j] >= 0) k1++;
		ndgrp[igr] = k1;
	}

	// write the number of element gouups into the output files
	int nelgrbou = 0;  // number of element groups having nodes at the boundary
	for(igr=0; igr<numElementGroups_; igr++) if(ndgrp[igr] > 0) nelgrbou++;

    fprintf(Lu_pBoundary, "%5d\n", numElementGroups_);
    fprintf(Lu_vBoundary, "%5d\n", numElementGroups_);

	// compute and output the velocity potentials, normal velocities and energy ientensitie at the nodes for each element group, the way is already paved for middle face elements, but forget the 2 sides for now
	Matrix<Complex> zpregrp(numNodes_, 2); // pressures at both sides
	Vector<Complex> zvelgrp(numNodes_);    // velocities at the middle face in normal direction
	Matrix<double> enegrp(numNodes_, 2);   // energy intensities at both sides

	for(igr=0; igr<numElementGroups_; igr++)
	{
	  /* thingrp denotes the element type of the group if it is -1 we have a eval element */
		if(ndgrp[igr] == 0 || thingrp[igr] == -1) continue;
		for(j=0; j<numNodes_; j++)
		{
			mdgrp[j] = 0;    // number of elements including the node
			nundgrp[j] = -1; // global nodal number of the node if it belongs to
			                 // the current element group
			zvelgrp[j].set(0.0, 0.0);
			for(k=0; k<2; k++)
			{
				zpregrp(j, k).set(0.0, 0.0);
				enegrp(j, k) = 0.0;
			}
		}
		for(iel=0; iel<numElements_; iel++)
		{
		  /* again if eval element ignore the element*/
			if(listElementProperty[iel] == 2 || indexOfElementGroup[igr] != listElementsElementGroup[iel]) continue;
			/* we are here at the node level for each element, 
			   before we looked at the collocation node level = elem midpoint level
			*/
			for(j=0; j<listNumberNodesPerElement[iel]; j++)
			{
				k = elementsConnectivity[iel][j];
				mdgrp[k]++;
				if(nundgrp[k] == -1) nundgrp[k] = extNumbersOfNodes[k];
				zpregrp(k, 0) += cVpotele(iel, 0);
				zvelgrp[k] += cVeloelo(iel, j);
				enegrp(k, 0) += rEninele(iel, j);
			}
		} // end of loop IEL

		for(j=0; j<numNodes_; j++)
		{
			if(mdgrp[j] > 0)
			{
				re1 = (double)mdgrp[j];
				zpregrp(j, 0).div_r(re1);
				zpregrp(j, 0).mul_i(rpfact_);
				zvelgrp[j].div_r(re1);
				enegrp(j, 0) /= re1;
			}
		}

        fprintf(Lu_pBoundary, "%5d %5d\n", indexOfElementGroup[igr], numberOfElementsInGroup[igr]);
        fprintf(Lu_vBoundary, "%5d %5d\n", indexOfElementGroup[igr], numberOfElementsInGroup[igr]);
        
        for(iel=0; iel<numElements_; iel++)
        {
            if(listElementProperty[iel] == 2 || indexOfElementGroup[igr] != listElementsElementGroup[iel]) continue;
	    /* 
	       currently we just write the output at the collocnodes, if you want the values for the
	       element nodes use the values calculated above
	    */
            Complex tmp;
            tmp=cVpotele(iel, 0);
            tmp.mul_i(rpfact_);
            fprintf(Lu_pBoundary, "%5d % E % E\n", extNumbersOfElements[iel], tmp.re(), tmp.im());
            tmp=cVeloele(iel, 0);
            fprintf(Lu_vBoundary, "%5d % E % E\n", extNumbersOfElements[iel], tmp.re(), tmp.im());
        }
	/* if you want dB values write Vmagdbph to a file */

		int j0 = -1;
		for(j=0; j<numNodes_; j++)
		{
			if(nundgrp[j] == -1) continue;
			k1 = nundgrp[j];
			j0++;

			// compute the magnitude, dB value and phase of the pressure, velocity
			// and energy intensity
			NC_Magnitude2dBdeg(Vmagdbph_pp_3d, zpregrp(j, 0), 0);
			NC_Magnitude2dBdeg(Vmagdbph_v_3d, zvelgrp[j], 1);
		} // end of loop J
	} // end of loop IGR
    
//    fprintf(Lu_pBoundaryVTK, "POINTS %d float\n", ndgrp[igr]);
//    for(i=0; i<numNodes_; i++)
//    {
//        if(nundgrp[j] == -1) continue;
//        fprintf(Lu_pBoundaryVTK, "%E %E %E\n", nodesCoordinates[i][0], nodesCoordinates[i][1], nodesCoordinates[i][2]);
//    }
//    fprintf(Lu_pBoundaryVTK,"POLYGONS numberOfelements numberofelements*(nodes+1) '\n");
//    for(i=0; i<numElements_; i++)
//    {
//        fprintf(Lu_pBoundaryVTK, "3 %d %d %d\n", elementsConnectivity[i][0],elementsConnectivity[i][1], elementsConnectivity[i][2]);
//    }

    fclose(Lu_pBoundary);
    fclose(Lu_vBoundary);
//    fclose(Lu_pBoundaryVTK);
}

// compute and write the sound pressure, velocity and energy intensities at nodes of the evaluation mesh (i. e. the internal points)
void NC_WriteResultsEvaluationGrid
(
	ofstream& NCout,
	Matrix<Complex>& cVpotele,  // velocity potential at collocnodes
	Matrix<Complex>& cVeloele   // particle velocity at collocnodes
)
{

  /* 11.12.2020
     kreiza: changed output of velocity to contain the velocitycomponent
     for each direction in 3D, at the evalgrid there is not obvious normal
     vector, thus the particle velocity is a vector of dimension 3
  */
	Complex zfacsourel, zquelinten, zprefree;
    FILE *lu_pEvalGrid, *lu_vEvalGrid;

    // open files
    sprintf(Dfname_3d, "be.out/be.%d/pEvalGrid",  currentFrequency_ + 1);
    lu_pEvalGrid = fopen(Dfname_3d, "w");
    if(!lu_pEvalGrid) NC_Error_Exit_1(NCout,  "Can not open the file ", Dfname_3d);
    sprintf(Dfname_3d, "be.out/be.%d/vEvalGrid",  currentFrequency_ + 1);
    lu_vEvalGrid = fopen(Dfname_3d, "w");
    if(!lu_vEvalGrid) NC_Error_Exit_1(NCout,  "Can not open the file ", Dfname_3d);
    
    fprintf(lu_pEvalGrid, "%s\n", versionNumber_.c_str());
    fprintf(lu_vEvalGrid, "%s\n", versionNumber_.c_str());

	// local variables and arrays
	int inp, i, j, ndip, idip;
	double enresu, re1;
	Complex zprip, z1;
	Vector<double> crdip(NDIM), veip(NDIM);
	Vector<Complex> zprint(numNodesOfEvaluationMesh_);
	Vector<int> nuinnode(numNodesOfEvaluationMesh_);
	Matrix<Complex> zveint(numNodesOfEvaluationMesh_, NDIM + 1);
	Matrix<double> denint(numNodes_, NDIM);


	// compute the vector NUINNODE
	j = 0;
	for(i=0; i<numNodes_; i++) if(isNodeMeshOrEvalNode[i] == -1) nuinnode[j++] = i;

	// initialize the result arrays
	for(i=0; i<numNodesOfEvaluationMesh_; i++)
	{
		zprint[i].set(0.0, 0.0);
		for(j=0; j<NDIM; j++) zveint(i, j).set(0.0, 0.0);
	}

	// computing potentials and velocities at internal points

	
	if(numInternalPointsClusters_ > 0) // SLFMBEM or MLFMBEM
	{
		// contributions of all elements to the results at the internal point
		NC_ContributionFMM(NCout, cVpotele, cVeloele, nuinnode, zprint,
			zveint);

		// contribution of the incident waves
		if(numIncidentPlaneWaves_ || numPointSources_) for(inp=0; inp<numNodesOfEvaluationMesh_; inp++)
		{
			// nodal number and coordinates of the current point
			ndip = nuinnode[inp];
			for(i=0; i<NDIM; i++) crdip[i] = nodesCoordinates[ndip][i];

			//z1 = NC_ContributionIncidentWaves(NCout, crdip, Zvelopoi_3d, 1);
			//zprint[inp] += z1;
			for(j=0; j<NDIM; j++) zveint(inp, j) += Zvelopoi_3d[j];
		}
	}
	else // TBEM
	{
		for(inp=0; inp<numNodesOfEvaluationMesh_; inp++)
		{
			// nodal number and coordinates of the current point
			ndip = nuinnode[inp];
			for(i=0; i<NDIM; i++) crdip[i] = nodesCoordinates[ndip][i];

			// contributions of all elements to the results at the internal point
			NC_ContributionTBEM(NCout, cVpotele, cVeloele, inp, ndip,
				crdip, zprint, zveint);

			// contribution of the incident waves
			if(numIncidentPlaneWaves_ || numPointSources_) 
			{
				//z1 = NC_ContributionIncidentWaves(NCout, crdip, Zvelopoi_3d, 1);
				//zprint[inp] += z1;
				for(j=0; j<NDIM; j++) zveint(inp, j) += Zvelopoi_3d[j];
			}
		} // end of loop INP
	} // end of ELSE

	/* 
	   zprint contains the pressure at the evalnode, zveint the components of the particle
           velocity
	*/
	
	// group number of the evaluation mesh
	int nugrinme = -1;
	for(i=0; i<numElements_; i++) if(listElementProperty[i] == 2)
	{
		nugrinme = listElementsElementGroup[i];
		break;
	}

	int ngrp = 0;
    if(numNodesOfEvaluationMesh_ > 0) ngrp = 1;
    fprintf(lu_vEvalGrid, "%5d\n", ngrp);

	if(numNodesOfEvaluationMesh_ > 0)
	{
        fprintf(lu_vEvalGrid, "%5d %5d\n", nugrinme, numNodesOfEvaluationMesh_);
	}

	// loop over nodes of the evaluation mesh
	for(inp=0; inp<numNodesOfEvaluationMesh_; inp++)
	  {
	    // external nodal number of the current node
	    ndip = extNumbersOfNodes[nuinnode[inp]];
	    
	    // components of velocity
	    for(j=0; j<NDIM; j++)
	      {
		NC_Magnitude2dBdeg(
				   Vmagdbph_v_3d, zveint(inp, j), 1
				   );
		veip[j] = Vmagdbph_v_3d[1];
	      }
	    
		// norm of the velocity
	    /*  changed by kreiza
	    z1.set(0.0, 0.0);
	    for(j=0; j<NDIM; j++) z1 += zveint(inp, j)*zveint(inp, j);
	    z1 = BLzsqrt(z1);
	    */
	    /* in dB */
	    /*  commented by kreiza dec 2020
	      NC_Magnitude2dBdeg(Vmagdbph_v_3d, z1, 1); 
	    
	      fprintf(lu_vEvalGrid, "%5d % E % E\n", ndip, z1.re(), z1.im());
	    */
	    fprintf(lu_vEvalGrid,"%5d %E %E %E %E %E %E\n", ndip,
		    zveint(inp,0).re(), zveint(inp,0).im(),
		    zveint(inp,1).re(), zveint(inp,1).im(),
		    zveint(inp,2).re(), zveint(inp,2).im());
	    
	} // end of loop INP
    
	fprintf(lu_pEvalGrid, "%5d\n", ngrp);
	if(numNodesOfEvaluationMesh_ > 0)
	  {
	    fprintf(lu_pEvalGrid, "%5d %5d\n", nugrinme, numNodesOfEvaluationMesh_);
	  }
	
	// loop over nodes of the evaluation mesh
	for(inp=0; inp<numNodesOfEvaluationMesh_; inp++)
	{
		// internal and external nodal number of the current node
		idip = nuinnode[inp];
		ndip = extNumbersOfNodes[idip];

		// sound pressure at the current node
		zprip = zprint[inp];
		zprip.mul_i(rpfact_);
        
		NC_Magnitude2dBdeg(Vmagdbph_pp_3d, zprip, 0);
		
		fprintf(lu_pEvalGrid, "%5d % E % E\n", ndip, zprip.re(), zprip.im());
	} // end of loop INP

	// loop over nodes of the evaluation mesh
	for(inp=0; inp<numNodesOfEvaluationMesh_; inp++)
	{
		// external nodal number of the current node
		ndip = extNumbersOfNodes[nuinnode[inp]];

		// compute the energy intensity at the current node
		z1 = zprint[inp];
		z1.mul_i(rpfact_);
		enresu = 0.0;
		for(j=0; j<NDIM; j++)
		{
			re1 = z1.re()*zveint(inp, j).re() + z1.im()*zveint(inp, j).im();
			denint(nuinnode[inp], j) = re1;
			enresu += re1*re1;
		}
		enresu = sqrt(enresu);
	} // end of loop INP

	// close the files
    fclose(lu_pEvalGrid);
    fclose(lu_vEvalGrid);
}

// compute the contributions of the incident waves to the velocity potential and the particle velocity at a given point
Complex NC_ContributionIncidentWaves
(
	ofstream& NCout,
	Vector<double>&  sourpoi, // I: coordinates of the point
	Vector<Complex>& zvelinw, // O: contributions to velocity
	const int& ifcvelo        // = 0: only compute potential
// = 1: compute potential and velocity
)
{
  	int i, j, kref;
    double re1, reffac, wavruim = harmonicTimeFactor_*waveNumbers_;
    Complex zrc(0.0, 0.0), zg, zini, zv;
    Vector<double> vect1(NDIM), vect2(NDIM);
    /*
    if(ifcvelo) for(i=0; i<NDIM; i++) zvelinw[i].set(0.0, 0.0);
    
    // the incident plane wave must be reflected with respect to the symmetric planes
    if(numIncidentPlaneWaves_ > 0)
    {
        Vector<double> sourposym(sourpoi, NDIM);
        if(numSymmetricPlanes_ > 0)
        {
            for(i=0; i<NDIM; i++) if(listSymmetryPlanes[i] != 0) sourposym[i] -=listSymmetryPlanesCoordinates[i];
        }
        
        // loop over incident plane waves
        for(i=0; i<numIncidentPlaneWaves_; i++)
        {
            // direction of the cuurent incident plane wave
            for(j=0; j<NDIM; j++) vect1[j] = dirinw[i][j];
            
            // reflection of the incident plane wave with respect to the symmetric planes
            zini = zinwav[i];
            reffac = 1.0;
            for(kref=0; kref<numReflectionsOfElements_; kref++)
            {
                if(kref > 0) reffac = NC_ReflectPointVector(vect1, kref, reffac, 0);
                
                re1 = vect1*sourposym*wavruim;
                zg.set(cos(re1)*reffac, sin(re1)*reffac);
                zg.mul_c(zini);
                zrc += zg;
                
                if(ifcvelo)
                {
                    zv = zg;
                    zv.mul_i(wavruim);
                    for(j=0; j<NDIM; j++) zvelinw[j] += zv*vect1[j];
                }
            }
        } // end of loop I
    }
    
    if(numPointSources_ > 0)
    {
        double aglen, re2;
        
        // loop over point sources
        for(i=0; i<numPointSources_; i++)
        {
            // coordinates of the point source
            for(j=0; j<NDIM; j++) vect1[j] = coorps[i][j];
            
            // reflection of the point source with respect to the symmetric planes
            zini = zposoc[i];
            reffac = 1.0;
            for(kref=0; kref<numReflectionsOfElements_; kref++)
            {
                if(kref > 0) reffac = NC_ReflectPointVector(vect1, kref, reffac, 1);
                
                for(j=0; j<NDIM; j++) vect2[j] = sourpoi[j] - vect1[j];
                Vcnorm_dim3_(vect2, aglen);
                
                re1 = wavruim*aglen;
                re2 = reffac/PI4/aglen;
                zg.set(cos(re1)*re2, sin(re1)*re2);
                zg.mul_c(zini);
                zrc += zg;
                
                if(ifcvelo) 
                {
                    zv.set(-1.0/aglen, wavruim);
                    zv *= zg;
                    for(j=0; j<NDIM; j++) zvelinw[j] += zv*vect2[j];
                }
            }
        } // end of loop I
    }
    */
    return zrc;
}

// compute contributions of the elements to the sound pressure, velocity and energy intensities at a given internal point by using the traditional BEM
void NC_ContributionTBEM
(
	ofstream& NCout,
	Matrix<Complex>& cVpotele,
	Matrix<Complex>& cVeloele,
	int& inp,                       // I: local number of the internal point
	int& ndip,                      // I: global number of the internal point
	Vector<double>& crdip,          // I: coordinate of the internal point
	Vector<Complex>& zprint,		// O: pressure at the internal points
	Matrix<Complex>& zveint			// O: velocity at the internal points
)
{
	// local variables and arrays
	int ie, i, j, inode, kref, nsbe, jsb, ireffac;
	double areli;
	Complex zprip;
	Vector<double> csip(NNODPE), etap(NNODPE);
	Vector<Complex> zveip(NDIM);

	// arrays and variables of elements
	int ivi_d, ifadmii;
	Complex zpotd, zpoto;
	Vector<int> idofaddre(2);
	Vector<double> centeri(NDIM), norveci(NDIM*(NNODPE + 1));
	Vector<Complex> zadmi(2);
	Vector<Complex> zvi_d(NNODPE);
	Matrix<double> crdeli(NNODPE, NDIM);

	// arrays of subelements
	Vector<double> xisbp(MSBE), etsbp(MSBE), fasbp(MSBE);
	Vector<int> ngsbp(MSBE);


	// loop over the elements
	for(ie=0; ie<numElements_; ie++)
	{
		if(listElementProperty[ie] == 2) continue;

		// data of the current element
		NC_ComputeElementData(ie, 0, inode, ivi_d,
			ifadmii, idofaddre, centeri, norveci, zadmi,
			zvi_d, crdeli);
		areli = areael[ie];

        zpoto = cVpotele(ie, 0);
        for(i=0; i<inode; i++) zvi_d[i] = cVeloele(ie, i);

		ivi_d = 0;
		for(i=0; i<inode; i++) 
		{
			if(zvi_d[i].norm() > EPSY)
			{
				ivi_d = 1;
				break;
			}
		}

		// loop over reflections of the current element
		ireffac = 1;
		for(kref=0; kref<numReflectionsOfElements_; kref++)
		{
			// reflection of the element
			if(kref > 0) NC_ReflectElement(kref, inode, 0, ireffac, centeri,
				norveci, zvi_d, crdeli);
			if(ireffac == -1) zpotd = zpoto*(-1.0); else zpotd = zpoto;
		
			// local coordinates of the nodes
			if(inode == NETYP3)	for(i=0; i<inode; i++)
			{
				csip[i] = CSI6[i];
				etap[i] = ETA6[i];
			}
			else for(i=0; i<inode; i++)
			{
				csip[i] = CSI8[i];
				etap[i] = ETA8[i];
			}

			// generate the subelements
			nsbe = NC_GenerateSubelements(NCout, ndip, crdip, ie, crdeli, areli,
				inode, csip, etap, inode, xisbp,
				etsbp, fasbp, ngsbp);

			// initialize the countributions of the current element
			zprip.set(0.0, 0.0);
			for(i=0; i<NDIM; i++) zveip[i].set(0.0, 0.0);

			// loop over the subelements
			for(jsb=0; jsb<nsbe; jsb++)
			{
				zprip += NC_IntegrationTBEM(NCout, crdip, xisbp[jsb], etsbp[jsb], ngsbp[jsb],
					fasbp[jsb], inode, crdeli, zpotd, ivi_d,
					zvi_d, zveip, inode, ie);
			}

			// add ZPRIP and ZVEIP to the global arrays
			if(Tao_ != 1.0)
			{
				zprip.mul_r(Tao_);
				for(j=0; j<NDIM; j++) zveip[j].mul_r(Tao_);
			}
			zprint[inp] += zprip;
			for(j=0; j<NDIM; j++) zveint(inp, j) += zveip[j]; 

		} // end of loop KREF
	} // end of loop IE
}

// compute contributions of the boundary elements to the sound pressures and velocities at the internal points by using the FMBEM
void NC_ContributionFMM
(
	ofstream& NCout,
	Matrix<Complex>& cVpotele,
	Matrix<Complex>& cVeloele,
	Vector<int>& nuinnode,
	Vector<Complex>& zprint,		// O: pressure at the internal points
	Matrix<Complex>& zveint			// O: velocity at the internal points
)
{
	// local variables and arrays
	int ie, i, j, inode, nsbe, jsb, ipc, ic, isum, k, ig, ipp, inp, ndip;
	double areli, di0, dik, scprd, wavruim = waveNumbers_*harmonicTimeFactor_;
	Complex zprip, zmun, zi, zk, *ztvcwork = nullptr;

	switch(methodFMM_)
	{
	case 1: // SLFMBEM
		ztvcwork = zT_vc;
		break;
	case 3: // DMLFMBEM

		ztvcwork = clulevarry[0].zwkT;

		// for MLFMBEM assign the coarsest level to be the level used in this function
		ClustArray = clulevarry[0].ClastArLv;

		numOriginalClusters_ = clulevarry[0].nClustOLv;
		numOriginalReflectedClusters_ = clulevarry[0].nClustSLv;

		numExpansionTerms_ = clulevarry[0].nExpaTermLv;
		numIntegrationPointsUnitSphere_ = clulevarry[0].nPoinSpheLv;
		maxClusterRadiusBE_ = clulevarry[0].RadiMaxLv;
		avgClusterRadiusBE_ = clulevarry[0].RadiAveLv;
		minClusterRadiusBE_ = clulevarry[0].RadiMinLv;

		weisphe = clulevarry[0].weisphe;
		uvcsphe = clulevarry[0].uvcsphe;
		break;
	}

	Vector<double> csip(NNODPE), etap(NNODPE), dwk0(NDIM), Pm(numExpansionTerms_),
		crdip(NDIM);
	Vector<Complex> zveip(NDIM), zH0m(numExpansionTerms_);
	Vector<bool> iffarclus(numOriginalReflectedClusters_);

	// arrays and variables of elements
	int ivi_d, ifadmii;
	Complex zpotd;
	Vector<int> idofaddre(2);
	Vector<double> centeri(NDIM), norveci(NDIM*(NNODPE + 1));
	Vector<Complex> zadmi(2);
	Vector<Complex> zvi_d(NNODPE), zdt_vct(numIntegrationPointsUnitSphere_);
	Matrix<double> crdeli(NNODPE, NDIM);

	// arrays of subelements
	Vector<double> xisbp(MSBE), etsbp(MSBE), fasbp(MSBE);
	Vector<int> ngsbp(MSBE);

	// cluster parameters
	int nrf_cl, rfa_cl;
	bool mir_cl;
	Vector<bool> rfdi_cl(NDIM);

    // loop over clusters of internal points
	for(ipc=0; ipc<numInternalPointsClusters_; ipc++) {
	  // compute the array indicating if a cluster is located in the far field
	  for(i=0; i<numOriginalReflectedClusters_; i++)
	    {
	      di0 = Dispoi_dim3_(ClustArray[i].CoorCent, ipcluarry[ipc].CoorCent);
	      dik = ClustArray[i].RadiClus + ipcluarry[ipc].RadiClus;
	      if(di0 > farFieldClusterFactor_*dik && di0 > minClusterDistance_) 
		iffarclus[i] = true; else iffarclus[i] = false;
	    }

	  // loop over clusters
	  isum = 0;
	  for(ic=0; ic<numOriginalReflectedClusters_; ic++)
	    {
	      if(iffarclus[ic]) // for a far field cluster using the FM-method
		{
		  // difference between the two cluster centers
		  for(j=0; j<NDIM; j++) dwk0[j] = ipcluarry[ipc].CoorCent[j] - 
					  ClustArray[ic].CoorCent[j];
		  Vcnorm_dim3_(dwk0, di0);
		  dik = di0*waveNumbers_;
		  
		  // compute the Hankel functions
		  NC_SphericalHankel(NCout, zH0m, numExpansionTerms_, dik);
		  
		  // compute the vector D*T
		  for(j=0; j<numIntegrationPointsUnitSphere_; j++)
		    {
		      scprd = Scprod_dim3_(dwk0, uvcsphe[j]);
		      
		      // compute the Legendre plynomes
		      NC_LegendrePolynomes(NCout, Pm, numExpansionTerms_, scprd);
		      
		      // compute the function Mu_n
		      zmun.set(0.0, 0.0);
		      zi.set(1.0, 0.0);
		      for(k=0; k<numExpansionTerms_; k++)
			{
			  zmun += zi*zH0m[k]*Pm[k]*(double)(2*k + 1);
			  zi.mul_i(harmonicTimeFactor_); // _ruimag!
			}
		      zdt_vct[j] = zmun*ztvcwork[isum + j]*weisphe[j];
		    } // end of loop J
		  
		      // loop over points of the current cluster of internal points
		  for(ipp=0; ipp<ipcluarry[ipc].NumOfIps; ipp++)
		    {
		      // initialize the countributions of the current cluster to potential and velocity of the current internal point
		      zprip.set(0.0, 0.0);
		      for(i=0; i<NDIM; i++) zveip[i].set(0.0, 0.0);
		      
		      // local number and nodal number of the current internal point
		      inp = ipcluarry[ipc].NumsOfIps[ipp];
		      ndip = nuinnode[inp];
		      
		      // coordinate difference of the current point and the cluster center
		      for(i=0; i<NDIM; i++) crdip[i] = nodesCoordinates[ndip][i] - 
					      ipcluarry[ipc].CoorCent[i];
		      
		      // compute the contributions of the current cluster
		      for(j=0; j<numIntegrationPointsUnitSphere_; j++) {
			scprd = Scprod_dim3_(crdip, uvcsphe[j])*wavruim;
			zi.set(cos(scprd), sin(scprd));
			zi *= zdt_vct[j];
			zprip += zi;
			for(k=0; k<NDIM; k++)
			  {
			    dik = wavruim*uvcsphe[j][k];
			    zk.set(-zi.im()*dik, zi.re()*dik);
			    zveip[k] += zk;
			  }
		      } // end of loop J
		      
		      // add ZPRIP and ZVEIP to the global arrays
		      zprint[inp] += zprip;
		      for(j=0; j<NDIM; j++) zveint(inp, j) += zveip[j];  
		    } // end of loop IPP
		}
	      else // for a near field cluster using the conventional method
		{
		  // cluster parameters
		  nrf_cl = ClustArray[ic].nuref;
		  rfa_cl = ClustArray[ic].rffac;
		  if(nrf_cl)
		    {
		      mir_cl = ClustArray[ic].ifmirro;
		      for(i=0; i<NDIM; i++) rfdi_cl[i] = ClustArray[ic].ifrfdi[i];
		    }
		  
		  // loop over members of the current cluster of internal points
		  for(ipp=0; ipp<ipcluarry[ipc].NumOfIps; ipp++)
		    {
		      // local number and nodal number of the current internal point
		      inp = ipcluarry[ipc].NumsOfIps[ipp];
		      ndip = nuinnode[inp];
		      
		      // coordinates of the current point
		      for(i=0; i<NDIM; i++) crdip[i] = nodesCoordinates[ndip][i];
		      
		      // loop over elements in the current cluster
		      for(ig=0; ig<ClustArray[ic].NumOfEl; ig++)
			{
			  ie = ClustArray[ic].NumsOfEl[ig]; // element number
			  
			  // data of the current element
			  NC_ComputeElementData(ie, 0, inode, ivi_d,
						ifadmii, idofaddre, centeri, norveci, zadmi,
						zvi_d, crdeli);
			  areli = areael[ie];
			  
			  zpotd = cVpotele(ie, 0);
			  for(k=0; k<inode; k++) zvi_d[k] = cVeloele(ie, k);
			  
			  ivi_d = 0;
			  for(k=0; k<inode; k++) 
			    {
			      if(zvi_d[k].norm() > EPSY)
				{
				  ivi_d = 1;
				  break;
				}
			    }
			  
			  // reflection of the element
			  if(nrf_cl) NC_ReflectElementFMM(rfa_cl, mir_cl, rfdi_cl, inode, 0,
							  centeri, norveci, zvi_d, crdeli);
			  if(nrf_cl && rfa_cl == -1) zpotd.nega();
			  
			  // local coordinates of the nodes
			  if(inode == NETYP3)	for(k=0; k<inode; k++)
						  {
						    csip[k] = CSI6[k];
						    etap[k] = ETA6[k];
						  }
			  else for(k=0; k<inode; k++)
				 {
				   csip[k] = CSI8[k];
				   etap[k] = ETA8[k];
				 }
			  
			  // generate the subelements
			  nsbe = NC_GenerateSubelements(NCout, ndip, crdip, ie, crdeli, areli,
							inode, csip, etap, inode, xisbp,
							etsbp, fasbp, ngsbp);
			  
			  // initialize the countributions of the current element
			  zprip.set(0.0, 0.0);
			  for(i=0; i<NDIM; i++) zveip[i].set(0.0, 0.0);
			  
			  // loop over the subelements
			  for(jsb=0; jsb<nsbe; jsb++)
			    zprip += NC_IntegrationTBEM(NCout, crdip, xisbp[jsb], etsbp[jsb],
							ngsbp[jsb], fasbp[jsb], inode, crdeli,
							zpotd, ivi_d, zvi_d, zveip,
							inode, ie);
			  
			  // add ZPRIP and ZVEIP to the global arrays
			  if(Tao_ != 1.0)
			    {
			      zprip.mul_r(Tao_);
			      for(j=0; j<NDIM; j++) zveip[j].mul_r(Tao_);
			    }
			  zprint[inp] += zprip;
			  for(j=0; j<NDIM; j++) zveint(inp, j) += zveip[j];
			      
			} // end of loop IG
		    } // end of loop IPP
		} // end of ELSE
	      isum += numIntegrationPointsUnitSphere_;
	    } // end of loop IC
	} // end of loop IPC
}

// perform the numerical integration over a subelemnt in order to compute its contributions to the potential and velocity at an internal point, return: contribution to the potential
Complex NC_IntegrationTBEM
(
	ofstream& NCout,
	Vector<double>& crdip,		//I: coordinates of the source point
	const double& xice,			//I: first local coord. of the ori. point of the subel.
	const double& etce,			//I: second local coord. of the ori. point of the subel.
	const int& ngse,			//I: order of the Gaussean integration
	const double& fase,			//I: factor of dimension of the current subelement
	const int& nvertj,			//I: number of vertex nodes of the current element
	Matrix<double> crdele,		//I: coordinates of nodes of the current element
	const Complex zpotee,		//I: potential(or diff. of potentials) at the element
	const int& ivi_d,			//I: flag for nonzero velocitis at the nodes of the el.
	Vector<Complex>& zveloe,	//I: velocity(or diff. of velocities) at nodes of the el.
	Vector<Complex>& zveine,	//O: contribution of the element to the velocity
	const int& inode,			//I: number of nodes of the current element
	const int& numel			//I: internal number of the current el.
)
{
    Complex zprine;
    int iforie, ngpoi, igp, j, i;
    double fase2, xis, ets, xio, eto, weih2, vjaco2, pleng, rk, re1, drdnx, vjacwe,
    wavruim = waveNumbers_*harmonicTimeFactor_;
    double csigau[36], etagau[36], weigau[36];
    Complex zvgp, zgrfu, z1, zgikr, zfac1, zfac2, zdgrdn, zdgrdx[NDIM],
    zdgdnx[NDIM];
    Vector<double> vect1p(NDIM);
    Vector<double> shfunx(NNODPE), shnds(NNODPE), shndt(NNODPE),
    elnorv(NDIM), crdpoi(NDIM), dxds(NDIM), dxdt(NDIM);
    
    // assign values to some variables
    fase2 = fase*fase;
    if(fabs(fase - 1.0) < EPSY) iforie = 1; // the element is not subdivided
    else iforie = 0; // the element is subdivided
    
    // generate the Gaussean points
    ngpoi = NC_ComputeGausseanPoints(
                                     csigau, etagau, weigau, nvertj, ngse
                                     );
    
    // loop over Gaussean points
    zprine.set(0.0, 0.0);
    for(igp=0; igp<ngpoi; igp++)
    {
        // local coordinates and weight of the current point in the original element
        xis = csigau[igp];
        ets = etagau[igp];
        if(iforie)
        {
            xio = xis;
            eto = ets;
            weih2 = weigau[igp];
        }
        else
        {
            xio = xice + xis*fase;
            eto = etce + ets*fase;
            weih2 = weigau[igp]*fase2;
        }
        
        // compute the parameters at the current Gaussean point
        vjaco2 = NC_ComputeParameters(shfunx, shnds, shndt, dxds, dxdt,
                                      crdpoi, elnorv, zvgp, crdele, ivi_d,
                                      zveloe, inode, xio, eto);
        
        // compute the Green function and its derivatives
        for(j=0; j<NDIM; j++) vect1p[j] = crdpoi[j] - crdip[j];
        Vcnorm_dim3_(vect1p, pleng);
        rk = pleng*wavruim;
        re1 = pleng*PI4;
        zgrfu.set(cos(rk)/re1, sin(rk)/re1);
        
        z1.set(-1.0/pleng, wavruim);
        zgikr = zgrfu*z1;
        
        drdnx = Scprod_dim3_(vect1p, elnorv);
        z1.set(drdnx*(waveNumbers_*waveNumbers_ - 3.0/pleng/pleng),
               wavruim*3.0*drdnx/pleng);
        zfac1 = zgrfu*z1;
        
        zfac2 = zgikr/(-pleng);
        
        // the first order derivatives of the Green function
        zdgrdn = zgikr*drdnx;
        for(j=0; j<NDIM; j++) zdgrdx[j] = zgikr*(-vect1p[j]);
        
        // the second derivatives of the Green function with respect to the normal direction and the coordinates of the source point
        for(j=0; j<NDIM; j++) zdgdnx[j] = zfac1*vect1p[j] + zfac2*elnorv[j];
        
        // product of the Jacobian and the weight factor
        vjacwe = vjaco2*weih2;
        
        // compute the contributions of the current Gaussean point to the potential and velocity at the internal point
        z1 = zpotee*vjacwe;
        zprine += zdgrdn*z1;
        for(j=0; j<NDIM; j++) zveine[j] += zdgdnx[j]*z1;
        if(ivi_d) for(i=0; i<inode; i++)
        {
            z1 = zveloe[i]*(vjacwe*shfunx[i]);
            zprine -= zgrfu*z1;
            for(j=0; j<NDIM; j++) zveine[j] -= zdgrdx[j]*z1;
        }
        
    } // end of loop IGP
    
    return zprine;
}

// compute the magnitude, dB value and phase of the sound pressure or the particle velocity
void NC_Magnitude2dBdeg
(
	Vector<double>& vmagdbph,  // O: magnitude, db and phase
	const Complex& zvalue,     // I: value of preesure or velocity
	const int& ifvelo          // I: = 0, pressure; != 0, velocity
)
{
    vmagdbph[0] = zvalue.norm();
    vmagdbph[2] = atan2(zvalue.im(), zvalue.re())*180/PI;
    
    switch(ifvelo)
    {
        case 0:  // for pressure
            if(vmagdbph[0] >= PREREF)
            {
                vmagdbph[1] = 20.0*(log10(vmagdbph[0]) - log10(PREREF));
            }
            else
            {
                vmagdbph[1] = 0.0;
            }
            break;
        case 1:
        default: // for velocity
            if(vmagdbph[0] >= VELREF)
            {
                vmagdbph[1] = 20.0*(log10(vmagdbph[0]) - log10(VELREF));
            }
            else
            {
                vmagdbph[1] = 0.0;
            }
            break;
    }
}
