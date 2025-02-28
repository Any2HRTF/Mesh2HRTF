/*===========================================================================*\
 *                                                                            *
 *  File name:      NC_Addresses.cpp                                          *
 *  Description:    computation of addresses of unknown DOFs of each element  *
 *                  in the global equation system                             *
 *  Author:         H. Ziegelwanger, W. Kreuzer and Z.S. Chen                 *
 *                                                                            *
 \*===========================================================================*/



//================================= I N C L U D E S ==================================
// local includes                                                                   //
#include"NC_ConstantsVariables.h"                                                   //
#include"NC_TypeDefinition.h"                                                       //
#include"NC_Arrays.h"                                                               //
//                                                                                  //
// system includes                                                                  //
#include<iostream>                                                                  //
#include<iomanip>                                                                   //
#include<string>                                                                    //
#include<fstream>  // file class functions                                          //
#include<math.h>                                                                    //
#include<cstddef>                                                                   //
using namespace std;                                                                //
//====================================================================================



// local functions
int NC_GenerateClustersSLFMM(ofstream& NCout);
void NC_AllocateNmtx(ofstream&);
void NC_ComputeClusterArraysSLFMM(ofstream&, Vector<int>&, Vector<int>&,int&);
int NC_GenerateClustersEvalGrid(ofstream&);
void NC_GenerateClustersFMM(ofstream&, int&, Vector<int>&,Vector<int>&, int&, Vector<int>&,Vector<int>&, const int&, Vector<int>&);
void NC_ClusterReflectionsSLFMM(ofstream&);
int NC_EstimateNumLevelsMLFMM(ofstream&);
void NC_GenerateClusterTreeMLFMM(ofstream&);
int NC_GenerateClustersAtLevelMLFMM(ofstream&, const int&);
void NC_GenerateClusterArrayAtLevelMLFMM(ofstream&, Vector<int>&, Vector<int>&,int&, const int&, Vector<int>&);
void NC_ClusterReflectionsAtLevelMLFMM(ofstream&, const int&);
void NC_AllocateDmtxMLFMM(ofstream&);
void NC_AllocateTmtxMLFMM(ofstream&);
void NC_AllocateSmtxMLFMM(ofstream&);
void NC_CheckClusters(ofstream&,	const int&,	const int&,	ElCluster*);
void NC_DeleteArrays(const int&, const int&, const int&);
void NC_CancelSmallClusters(Vector<int>&, Vector<int>&,  int&,  int&);



double ClusEdgL1;
int n_Pair_NeaF[10], // number of cluster pairs located in the frequncy near field
	n_Pair_FarD[10]; // number of cluster pairs located in the far field



// the main program for addressing computations
int NC_DeclareArrays
(
	ofstream& NCout,
	double* Freqs,
	bool estimate_ram
)
{
	double dLamdBeL = waveLength_/diameterBE_, d_ratio[10],
		thresholrat = 0.75,  // threshold for ratio n_Pair_NeaF[i]/n_Pair_FarD[i]
	                        // (i = 0, 1, ..., numClusterLevels_-1) _PARAMETER_
		thresholBEM = 80.0, // threshold of (wave length / diameter of BE mesh) for TBEM

		thresholFMM1 = 1.0; // threshold of (wave length / diameter of BE mesh) for SLFMBEM
	int Max_Dofs_TBEM = 20000, // maximum number of unknowns for traditional BEM
		Max_Dofs_FMM1 = 50000; // maximum number of unknowns for single level FMM
	int idifcomtyp = 1, i, npair_ne_sum = 1;
	static int imultipori, nlevmlfmori;

	// initialize the parameters to be computed
	numClusterLevels_ = 1;

	// determine the computation type for the current frequency
	if(methodBEM_ == 0) // for a TRBEM job
	{
		methodFMM_ = 0; // TRBEM is used
	}
	else // for a FMBEM job
	{
		if(dLamdBeL > thresholBEM)
		{
			methodFMM_ = 0; // TRBEM is used
		}
		else
		{
			if(methodBEM_ == 1) // for a SLFMBEM job
			{
				methodFMM_ = 1; // SLFMBEM is used
			}
			else if(methodBEM_ > 1) // for a MLFMBEM job
			{
				if(dLamdBeL > thresholFMM1)
				{
					methodFMM_ = 1; // SLFMBEM is used
				}
				else
				{
					numClusterLevels_ = NC_EstimateNumLevelsMLFMM(NCout);
					if(numClusterLevels_ == 1)
					{
						methodFMM_ = 1;
					} else {
						if(methodBEM_ == 4) {methodFMM_ = 3;}
					}
					if(numClusterLevels_ > 10) NC_Error_Exit_1(NCout, "numClusterLevels_ must <= 10!",
						"numClusterLevels_ = ", numClusterLevels_);
				}
			}
		}
	}

	// initialize
	if(currentFrequency_ == istart_) for(i=0; i<numClusterLevels_; i++) n_Pair_NeaF[i] = 1;

Labsettradibem:
	nlevtop_ = numClusterLevels_ - 1;

	if(currentFrequency_ > istart_) {
		// see if the computation type and number of levels is changed
		if(methodFMM_ == imultipori && numClusterLevels_ == nlevmlfmori) idifcomtyp = 0;

		// compute sum of cluster pairs which is located in the frequency near field
		npair_ne_sum = 0;
		for(i=0; i<nlevmlfmori; i++) npair_ne_sum += n_Pair_NeaF[i];

		// if the computation type is not changed and number cluster pairs in the
		// frequncy near field is zero RETURN
		if(idifcomtyp == 0 && npair_ne_sum == 0) return(0);

		// delete the old arrays generated in the following loop
		NC_DeleteArrays(imultipori, nlevmlfmori, 0);
	}

	// generate clusters or cluster tree
	for(i=0; i<numClusterLevels_; i++) n_Pair_FarD[i] = n_Pair_NeaF[i] = 0;
	numOriginalReflectedClusters_ = numOriginalClusters_ = 0;
	if(methodFMM_ == 1) // SLFMBEM
	{
		// generate clusters for SLFMBEM
		numOriginalClusters_ = NC_GenerateClustersSLFMM(NCout);

		if(numOriginalClusters_ > 0) {
			NC_CheckClusters(NCout, 0, numOriginalClusters_, clustarry);

			// reflection of the clusters
			if(numReflectionsOfElements_ > 1) NC_ClusterReflectionsSLFMM(NCout);
		} else {
			NC_Error_Exit_1(NCout, "Number of clusters must > 0!", "Number of clusters = ",
				numOriginalClusters_);
		}
	}
	else if(methodFMM_ >= 2) // MLFMBEM
	{
		// generate the cluster tree
		NC_GenerateClusterTreeMLFMM(NCout);

		numOriginalClusters_ = clulevarry[nlevtop_].nClustOLv;
		numOriginalReflectedClusters_ = clulevarry[nlevtop_].nClustSLv;
	}

	// set the cluster working array
	switch(methodFMM_)
	{
	case 1: // SLFMBEM
		ClustArray = clustarry;
		break;
	case 3: // DMLFMBEM
		ClustArray = clulevarry[nlevtop_].ClastArLv;
		break;
	}

    // compute the array JELIST (addresses of unknowns of each element in the global system)
    numRowsOfCoefficientMatrix_ = 0;  // number of rows of the equation system
    for(i=0; i<numElements_; i++)
    {
        if(listElementProperty[i] == 2) continue;
        jelist[i][0] = numRowsOfCoefficientMatrix_++;
    }

    // number of components of the coefficient matrix for the traditional BEM
    numComponentsOfCoefficientMatrix_ = numRowsOfCoefficientMatrix_*numRowsOfCoefficientMatrix_;
    
    // compute the near field coefficient matrix and
    // generate the clusters of internal points
    numInternalPointsClusters_ = 0;
    if(methodFMM_) // FMBEM is used
      {
	// compute the number of nonzeros of the near field coefficient matrix
	// kreiza: This is expensive an not really needed for the estimation
	//         of the entries. Trouble is that we have to estimate or
	//         determine the number of nonzero-entries later explicitely
	if( !estimate_ram )
	  NC_AllocateNmtx(NCout);
	
	// generate the clusters of the internal points
	if(numNodesOfEvaluationMesh_ > 0) numInternalPointsClusters_ = NC_GenerateClustersEvalGrid(NCout);
      }
    if ( !estimate_ram ) {
      // generate the arrays
      switch(methodFMM_)
	{
	case 1:
	  dmtxlev = new D_mtx_lev[1];
	  break;
	case 3:
	  dmtxlev = new D_mtx_lev[numClusterLevels_];
	  tmtxlev = new T_mtx_lev[numClusterLevels_];
	  smtxlev = new S_mtx_lev[numClusterLevels_];
	  break;
	}
      
      // create the right hand side vector
      zrhs = new Complex[numRowsOfCoefficientMatrix_];
    }
    // set original value of the computation type
    imultipori = methodFMM_;
    nlevmlfmori = numClusterLevels_;
    
    // compute the ratio of cluster paires in the frequency near field to the total number
    // of cluster pairs in the geometrical far field
    for(i=0; i<numClusterLevels_; i++) {
      if(n_Pair_NeaF[i] > 0) {
	d_ratio[i] = (double)n_Pair_NeaF[i]/(double)n_Pair_FarD[i];
      } else {
	d_ratio[i] = 0.0;
	//	NC_Error_Exit_1(NCout, "Frquency is too low, the single level FMBEM can not be used! Also, the number of elements for the traditional BEM is limited to 20000 Elements (~ 6.4Gb). Stopping NumCalc",
	//			"Frequency = ", frequency_);
      }
    }
    
    // write informations about canceled interacting cluster paires
    if(methodFMM_ && (n_Pair_NeaF[0] || n_Pair_NeaF[numClusterLevels_ - 1])) {
      cout << "\nInformations about the canceled interacting cluster pairs:" << endl;
      cout << "   Level   N. canceled pairs   N. interacting pairs   Cancel ratio (%)"
	   << endl;
      for(i=0; i<numClusterLevels_; i++) {
	cout << "   " << setw(3) << i << "     " << setw(10) << n_Pair_NeaF[i] <<
	  "           " << setw(10) << n_Pair_FarD[i] <<
	  "               " << setw(5) << d_ratio[i]*100.0 << endl;
      }
      NCout << "\nInformations about the canceled interacting cluster pairs:" << endl;
      NCout << "   Level   N. canceled pairs   N. interacting pairs   Cancel ratio (%)"
	    << endl;
      for(i=0; i<numClusterLevels_; i++) {
	NCout << "   " << setw(3) << i << "     " << setw(10) << n_Pair_NeaF[i] <<
	  "           " << setw(10) << n_Pair_FarD[i] <<
	  "               " << setw(5) << d_ratio[i]*100.0 << endl;
      }
    }
    
    if(methodFMM_ == 1) {
      if(d_ratio[0] > thresholrat) {
	methodFMM_ = 0;
	numClusterLevels_ = 1;
	idifcomtyp = 1;
	if(numRowsOfCoefficientMatrix_ < Max_Dofs_TBEM) {
	  goto Labsettradibem;
	} else {
	  NC_Error_Exit_1(NCout, "Frquency is too low, the single level FMBEM can not be used! Also, the number of elements for the traditioal BEM is limited to 20000 Elements (~ 6.4Gb). Stopping NumCalc",
			  "Frequency = ", frequency_);
	}
      }
    } else if(methodFMM_ > 1) {
      for(i=0; i<numClusterLevels_; i++)
	{
	  if(d_ratio[i] > thresholrat) {
	    methodFMM_ = 1;
	    numClusterLevels_ = 1;
	    idifcomtyp = 1;
	    if(numRowsOfCoefficientMatrix_ < Max_Dofs_FMM1) {
	      goto Labsettradibem;
	    } else {
	      NC_Error_Exit_1(NCout, "Frquency is too low, the multilevel FMBEM can not be used!",
			      "Frequency = ", frequency_);
	    }
	  }
	}
    }
    
    return(1);
}

// generate the element clusters for SLFMBEM
int NC_GenerateClustersSLFMM
(
	ofstream& NCout
)
{
	int i, j;

	int num_clus = 0; // number of clusters

	// number of elements of each cluster and numbers of these elements
	Vector<int> nel_clus(numElementsOfBoundaryMesh_), nuel_clus(numElementsOfBoundaryMesh_), nfath_clus(0);

	// number of BE groups
	int nbegrp = 0;
	for(i=0; i<numElementGroups_; i++) if(propertyOfGroup[i] == 0) nbegrp++;

	// numbers of BE groups
	Vector<int> nubegrp(nbegrp);
	j = 0;
	for(i=0; i<numElementGroups_; i++) if(propertyOfGroup[i] == 0) nubegrp[j++] = i;

	// maximal number of elements of an element group
	int melpergrp = 0;
	for(i=0; i<nbegrp; i++)
		if(numberOfElementsInGroup[nubegrp[i]] > melpergrp) melpergrp = numberOfElementsInGroup[nubegrp[i]];

	// numbers of elements of each BE group
	Vector<int> nuelbegrp(melpergrp);

	// average number of elements in a cluster
	double clus_nels = sqrt((double)(numElementsOfBoundaryMesh_)); //_TEST_!

	// compute the cluster edge length
	if(ClusEdgL0_ <= 0.0) // cluster edge length is not input by the user
	{
		ClusEdgL1 = sqrt(clus_nels*averageElementArea_);
	}
	else // it is input by the user
	{
		ClusEdgL1 = ClusEdgL0_;
		if(ClusEdgL1 < sqrt(4.0*averageElementArea_))
		{
			NC_Error_Exit_2(NCout, "The read cluster diameter is too small!",
				"Inputed value of cluster diameter = ", ClusEdgL1,
				"The smallst allowable value       = ", sqrt(4.0*averageElementArea_));
		}
		NCout << "\nClusEdgL (Inputed value of cluster diameter) = " << ClusEdgL1 << endl;
		cout << "ClusEdgL1 (Inputed value of cluster diameter) = " << ClusEdgL1 << endl;
		cout << endl;
	}

	// generate clusters
	NC_GenerateClustersFMM(NCout, num_clus, nel_clus, nuel_clus, nbegrp,
				 nubegrp, nuelbegrp, 0, nfath_clus);

	// compute the arrays for storing data of the clusters for SLFMBEM
	NC_ComputeClusterArraysSLFMM(NCout, nel_clus, nuel_clus, num_clus);

	return(num_clus);
}

// compute the arrays for storing data of the clusters for SLFMBEM
void NC_ComputeClusterArraysSLFMM
(
	ofstream& NCout,
	Vector<int>& nel_clus,		//I: number of elements in each cluster
	Vector<int>& nuel_clus,		//I: numbers of elements in the clusters
	int& num_clus				//I: number of clusters
)
{
    int i, j, k, k1, l_sum, l1, m, ie0, ind0;
    int nel_max = 0, n_nod;
    double wk0;
    double dareadis, dispo;
    Vector<double> crcent(NDIM);

    // maximum number of elements in a cluster
    for(i=0; i<num_clus; i++) if(nel_clus[i] > nel_max) nel_max = nel_clus[i];

    // vector to store nodal numbers of each cluster
    Vector<int> nu_nod_cl(nel_max*NNODPE);

    // create the cluster array
    numOriginalReflectedClusters_ = num_clus*numReflectionsOfElements_;
    clustarry = new ElCluster[numOriginalReflectedClusters_];

    // loop over clusters
    l_sum = 0;
    for(i=0; i<num_clus; i++)
    {
        // number of elements of the cluster
        clustarry[i].NumOfDOFs = clustarry[i].NumOfEl = nel_clus[i];

        // numbers of the elements of the cluster
        clustarry[i].NumsOfEl = new int[nel_clus[i]];
        for(j=0; j<nel_clus[i]; j++) clustarry[i].NumsOfEl[j] = nuel_clus[l_sum++];

        // number of the element group to which the cluster belongs
        ie0 = clustarry[i].NumsOfEl[0];
        clustarry[i].NuElGr = listElementsElementGroup[ie0];

        // surface or middle face elements ( = 0: surface els; = 1: middle face els)
        clustarry[i].listElementPropertyEl = listElementProperty[ie0];

        // if all elements of the cluster are of the same number of nodes
        clustarry[i].IfMonoEl = true;
        ind0 = listNumberNodesPerElement[ie0];
        for(j=1; j<nel_clus[i]; j++)
        {
            if(listNumberNodesPerElement[clustarry[i].NumsOfEl[j]] != ind0)
            {
                clustarry[i].IfMonoEl = false;
                break;
            }
        }

        // if admittance boundary condition are prescribed
        if(ibval[ie0] == 2 || ibval[ie0] == 4) clustarry[i].IfAdmiBc = true;
        else clustarry[i].IfAdmiBc = false;

        // number of unknown DOFs of the cluster
        clustarry[i].NDOFsPeEl = 1;

        // compute number and numbers of the nodes of the current cluster
        n_nod = 0;
        for(j=0; j<nel_clus[i]; j++)
        {
            k1 = clustarry[i].NumsOfEl[j];
            for(k=0; k<listNumberNodesPerElement[k1]; k++)
            {
                l1 = elementsConnectivity[k1][k];
                for(m=0; m<n_nod; m++) if(nu_nod_cl[m] == l1) goto LbclustAray0;
                nu_nod_cl[n_nod++] = l1;
            LbclustAray0: continue;
            }
        }

        // compute coordinates of the center of the cluster
        crcent = 0.0;
        for(j=0; j<n_nod; j++)
        {
            for(k=0; k<NDIM; k++) crcent[k] += nodesCoordinates[nu_nod_cl[j]][k];
        }
        for(k=0; k<NDIM; k++) clustarry[i].CoorCent[k] = crcent[k]/(double)n_nod;

        // compute radius of the cluster
        dareadis = 0.0;
        for(j=0; j<n_nod; j++)
        {
            dispo = Dispoi_dim3_(clustarry[i].CoorCent, nodesCoordinates[nu_nod_cl[j]]);
            if(dispo > dareadis) dareadis = dispo;
        }
        clustarry[i].RadiClus = dareadis;
    } // end of loop I

    // computer the maximum and minimum cluster radius
    avgClusterRadiusBE_ = maxClusterRadiusBE_ = 0;
    minClusterRadiusBE_ = 1.0e40;
    for(i=0; i<num_clus; i++)
    {
        if(clustarry[i].RadiClus > maxClusterRadiusBE_) maxClusterRadiusBE_ = clustarry[i].RadiClus;
        if(clustarry[i].RadiClus < minClusterRadiusBE_) minClusterRadiusBE_ = clustarry[i].RadiClus;
        avgClusterRadiusBE_ += clustarry[i].RadiClus;
    }
    avgClusterRadiusBE_ /= (double)num_clus;

    // computer the near and far field clusters
    if(numReflectionsOfElements_ == 1)
    {
        Matrix<bool> ifarclus(num_clus, num_clus, false);
        for(i=0; i<num_clus; i++) for(j=0; j<i; j++)
        {
            dispo = Dispoi_dim3_(clustarry[i].CoorCent, clustarry[j].CoorCent);
            wk0 = clustarry[i].RadiClus + clustarry[j].RadiClus;

            if(dispo > farFieldClusterFactor_*wk0) {
                if(dispo > minClusterDistance_) {
                    ifarclus(j, i) = ifarclus(i, j) = true;
                } else {
                    n_Pair_NeaF[0]++;
                }
                n_Pair_FarD[0]++;
            }
        }

        int nnea, nfar;
        for(i=0; i<num_clus; i++)
        {
            clustarry[i].NumNeaClus = clustarry[i].NumFarClus = 0;
            for(j=0; j<num_clus; j++)
            {
                if(ifarclus(i, j)) {clustarry[i].NumFarClus++;}
                else {clustarry[i].NumNeaClus++;}
            }

            clustarry[i].NumsNeaClus = new int[clustarry[i].NumNeaClus];
            clustarry[i].NumsFarClus = new int[clustarry[i].NumFarClus];
            nnea = nfar = 0;
            for(j=0; j<num_clus; j++)
            {
                if(ifarclus(i, j)) clustarry[i].NumsFarClus[nfar++] = j;
                else clustarry[i].NumsNeaClus[nnea++] = j;
            }
        }
    } // end of NREFL == 1

    // create reflection informations
    for(i=0; i<num_clus; i++)
    {
        clustarry[i].OriClust = i;
        clustarry[i].nuref = 0;
        clustarry[i].rffac = 1;
        clustarry[i].ifmirro = false;
        for(j=0; j<NDIM; j++) clustarry[i].ifrfdi[j] = false;
    }
}

// reflections of the element clusters with respect to the symmetric planes (for SLFMBEM)
void NC_ClusterReflectionsSLFMM
(
	ofstream& NCout
 )
{
    int i, j, i0 = 0, icl, nucl, kref, nrdir = 0, nel;
    double dispo, wk0;
    int nrefi[NDIM];
    int nrefdi[] = {0, 0, 1, 0, 2, 0, 1, 0};

    // generate the reflected clusters
    nucl = numOriginalClusters_;  // number of the generated cluster (updated)
    for(icl=0; icl<numOriginalClusters_; icl++)
    {
        for(i=0; i<NDIM; i++) nrefi[i] = 0;
        for(kref=1; kref<numReflectionsOfElements_; kref++)
        {
            // direction of the current reflection
            switch(numReflectionsOfElements_)
            {
                case 8:  // there are 3 symmetric planes
                    nrdir = nrefdi[kref];
                    nrefi[nrdir]++;
                    break;
                case 4:  // there are 2 symmetric planes
                    for(i=0; i<NDIM; i++) if(listSymmetryPlanes[i] == 0)
                    {
                        i0 = i + 1;
                        break;
                    }
                    nrdir = i0 + nrefdi[kref];
                    if(nrdir >= NDIM) nrdir -= NDIM;
                    nrefi[nrdir]++;
                    break;
                case 2:  // there is 1 symmetric plane
                    for(i=0; i<NDIM; i++) if(listSymmetryPlanes[i] != 0)
                    {
                        nrdir = i;
                        break;
                    }
                    nrefi[nrdir]++;
            } // end of SWITCH

            // parameters of the reflected cluster
            clustarry[nucl].NuElGr = clustarry[icl].NuElGr;
            clustarry[nucl].listElementPropertyEl = clustarry[icl].listElementPropertyEl;
            clustarry[nucl].IfMonoEl = clustarry[icl].IfMonoEl;
            clustarry[nucl].IfAdmiBc = clustarry[icl].IfAdmiBc;
            clustarry[nucl].NumOfDOFs = clustarry[icl].NumOfDOFs;
            clustarry[nucl].NDOFsPeEl = clustarry[icl].NDOFsPeEl;

            nel = clustarry[nucl].NumOfEl = clustarry[icl].NumOfEl;
            clustarry[nucl].NumsOfEl = new int[nel];
            for(i=0; i<nel; i++) clustarry[nucl].NumsOfEl[i] = clustarry[icl].NumsOfEl[i];

            clustarry[nucl].RadiClus = clustarry[icl].RadiClus;

            for(i=0; i<NDIM; i++)
            {
                if(nrefi[i]%2) clustarry[nucl].ifrfdi[i] = true;
                else clustarry[nucl].ifrfdi[i] = false;
                if(clustarry[nucl].ifrfdi[i])
                    clustarry[nucl].CoorCent[i] = 2.0*listSymmetryPlanesCoordinates[i] - clustarry[icl].CoorCent[i];
                else
                    clustarry[nucl].CoorCent[i] = clustarry[icl].CoorCent[i];
            }

            clustarry[nucl].OriClust = icl;
            clustarry[nucl].nuref = kref;
            if(kref%2) clustarry[nucl].ifmirro = true;
            else  clustarry[nucl].ifmirro = false;

            clustarry[nucl].rffac = 1;
            for(i=0; i<NDIM; i++) if(clustarry[nucl].ifrfdi[i] && listSymmetryPlanes[i] < 0)
                clustarry[nucl].rffac *= -1;

            nucl++;
        } // end of KREF
    } // end of loop ICL

    // compute the far and near field clusters of each cluster
    Matrix<bool> ifarclus(numOriginalReflectedClusters_, numOriginalReflectedClusters_, false);
    for(i=0; i<numOriginalReflectedClusters_; i++) for(j=0; j<i; j++)
    {
        dispo = Dispoi_dim3_(clustarry[i].CoorCent, clustarry[j].CoorCent);
        wk0 = clustarry[i].RadiClus + clustarry[j].RadiClus;

        if(dispo > farFieldClusterFactor_*wk0) {
            if(dispo > minClusterDistance_) {
                ifarclus(j, i) = ifarclus(i, j) = true;
            } else {
                n_Pair_NeaF[0]++;
            }
            n_Pair_FarD[0]++;
        }
    }

    int nnea, nfar;
    for(i=0; i<numOriginalReflectedClusters_; i++)
    {
        clustarry[i].NumNeaClus = clustarry[i].NumFarClus = 0;
        for(j=0; j<numOriginalReflectedClusters_; j++)
        {
            if(ifarclus(i, j)) clustarry[i].NumFarClus++;
            else clustarry[i].NumNeaClus++;
        }

        clustarry[i].NumsNeaClus = new int[clustarry[i].NumNeaClus];
        clustarry[i].NumsFarClus = new int[clustarry[i].NumFarClus];
        nnea = nfar = 0;
        for(j=0; j<numOriginalReflectedClusters_; j++)
        {
            if(ifarclus(i, j)) clustarry[i].NumsFarClus[nfar++] = j;
            else clustarry[i].NumsNeaClus[nnea++] = j;
        }
    }
}

// generate clusters for the SLFMBEM and MLFMBEM
void NC_GenerateClustersFMM
(
	ofstream& NCout,
	int& num_clus,				//O: number of clusters
	Vector<int>& nel_clus,		//O: number of elements in each cluster
	Vector<int>& nuel_clus,		//O: numbers of elements in each cluster
	int& nbegr_fcl,				//I: number of boundary element groups/father clusters
	Vector<int>& nubegrp,		//I: numbers of boundary element groups
	Vector<int>& nuelbegrp,		//w: numbers of elements of each BE group
	const int& nu_lev,			//I: number of the current level
	Vector<int>& nfath_clus		//O: number of the father clusters
)
{
  int i, j, k, l, ibg, nelgri, nodgri, ieli, ndij, inodi, nclus, l1, iclus, i1;
  int idexel = 0;
  Vector<int> nundgri(numNodesOfBoundaryMesh_), ndiv_xyz(NDIM);
  double d0, d1, dCluEdgLv = ClusEdgL1/pow(2.0, nu_lev);
  Vector<double> xyz_max(NDIM), xyz_min(NDIM), dif_xyz(NDIM);
  
  // compute (1) number of clusters; (2) number of elements of each cluster
  //         (3) numbers of elements of each cluster
  for(ibg=0; ibg<nbegr_fcl; ibg++) {
    
    // number of elements and numbers of elements of the current element group
    nelgri = 0;
    nuelbegrp = 0;
    if(nu_lev == 0) {
      for(j=0; j<numElements_; j++)
	if(indexOfElementGroup[nubegrp[ibg]] == listElementsElementGroup[j])
	  nuelbegrp[nelgri++] = j;
    }
        else
        {
            nelgri = clulevarry[nu_lev - 1].ClastArLv[ibg].NumOfEl;
            for(j=0; j<nelgri; j++) nuelbegrp[j] =
                clulevarry[nu_lev - 1].ClastArLv[ibg].NumsOfEl[j];
        }

        // compute number of nodes of the current element group and store numbers of these nodes
        nodgri = 0;
        for(i=0; i<nelgri; i++)
        {
            ieli = nuelbegrp[i];
            inodi = listNumberNodesPerElement[ieli];
            if(inodi > NETYP4) inodi /= 2;
            for(j=0; j<inodi; j++)
            {
                ndij = elementsConnectivity[ieli][j];
                for(k=0; k<nodgri; k++) if(nundgri[k] == ndij) goto Lab0clusMetho1;
                nundgri[nodgri++] = ndij;
            Lab0clusMetho1: continue;
            }
        }
	// nundgri contains now the nodenumbers of the nodes in the current element group


        // compute the maximum and minimum coordinates of the nodes of the current group
        xyz_max = -1.0e20;
        xyz_min = 1.0e20;
        for(i=0; i<nodgri; i++)
        {
            ndij = nundgri[i];
            for(j=0; j<NDIM; j++)
            {
                if(nodesCoordinates[ndij][j] > xyz_max[j]) xyz_max[j] = nodesCoordinates[ndij][j];
                if(nodesCoordinates[ndij][j] < xyz_min[j]) xyz_min[j] = nodesCoordinates[ndij][j];
            }
        }

        for(j=0; j<NDIM; j++)
        {
            xyz_max[j] += delta_;
            xyz_min[j] -= delta_;
        }
        for(j=0; j<NDIM; j++) dif_xyz[j] = xyz_max[j] - xyz_min[j];

        // compute number of divisions in each coordinate direction
        if(nu_lev == 0)
        {
            for(j=0; j<NDIM; j++)
            {
                d0 = dif_xyz[j]/ClusEdgL1;
                i = (int)(d0);
                if(i == 0) i = 1; else if(d0 - (double)(i) > 0.5) i++;
                ndiv_xyz[j] = i;
                dif_xyz[j] /= (double)(i);
            }
        }
        else
        {
            for(j=0; j<NDIM; j++)
            {
                d0 = dif_xyz[j]/dCluEdgLv;
                i = (int)(d0);
                if(i == 0) i = 1;
                else
                {
                    if(d0 - (double)(i) > 0.5) i++;
                    if(i > 2) i = 2;
                }
                ndiv_xyz[j] = i;
                dif_xyz[j] /= (double)(i);
            }
        }

        // number of clusters in the current element group
        nclus = ndiv_xyz[0]*ndiv_xyz[1]*ndiv_xyz[2];

        // number of elements in each cluster and numbers of these elements
        Vector<int> ne_clus(nclus, 0), nue_clus(nelgri);

        // compute numbers of the interval in the three coordinate directions for each element of the current group, in which it is located
        Matrix<int> nuindxyz(nelgri, NDIM, -1);
        for(j=0; j<NDIM; j++)
        {
            d1 = xyz_min[j];
            for(k=0; k<ndiv_xyz[j]; k++)
            {
	      d0 = d1; // changed by kreiza
	      d1 += dif_xyz[j];
	      //                d0 = d1 - dif_xyz[j];
                for(i=0; i<nelgri; i++)
                {
                    ieli = nuelbegrp[i];
                    if(nuindxyz(i, j) == -1 && centel[ieli][j] > d0 - epsilon_ &&
                       centel[ieli][j] < d1 + epsilon_) nuindxyz(i, j) = k;
                }
            }
        }

        // check the results
        for(i=0; i<nelgri; i++) for(j=0; j<NDIM; j++) {
            if(nuindxyz(i, j) == -1) {
                NC_Error_Exit_2(NCout, "The interval number of an element is not found!",
                                "number of the element = ", extNumbersOfElements[i],
                                "number of the direction = ", j);
            }
        }

        // compute (1) number of elements in each cluster
        //         (2) numbers of these elements
        iclus = 0;
        i1 = 0;
        for(i=0; i<ndiv_xyz[0]; i++) for(j=0; j<ndiv_xyz[1]; j++)
            for(k=0; k<ndiv_xyz[2]; k++)
            {
                l1 = 0;
                for(l=0; l<nelgri; l++)
                    if(nuindxyz(l, 0) == i && nuindxyz(l, 1) == j && nuindxyz(l, 2) == k)
                    {
                        l1++;
                        nue_clus[i1++] = nuelbegrp[l];
                    }
                ne_clus[iclus++] = l1;
            }
        // if there are very small clusters, combine them with the nearest
        // bigger one
        NC_CancelSmallClusters(ne_clus, nue_clus, nclus, ibg);



        // store the results
        i1 = 0;
        for(i=0; i<nclus; i++)
        {
            if(ne_clus[i] == 0) continue;

            if(nu_lev > 0) nfath_clus[num_clus] = ibg;
            nel_clus[num_clus++] = ne_clus[i];
            for(j=0; j<ne_clus[i]; j++) nuel_clus[idexel++] = nue_clus[i1++];
        }

    } // end of loop IBG
}


// estimate number of levels for MLFMBEM
int NC_EstimateNumLevelsMLFMM
(
	ofstream& NCout
)
{
	int n_levs;
	double ele_leng = sqrt(averageElementArea_), d_levs,
		d_min = 25.0; // _ACHTUNG_ _PARAMETER_
	double nclusroot, // estimated number of clusters at the root level
		clu_l_min;    // estimated diameter of clusters at the finest level

	// without interpolation/filtering procedures: unsymmetric clustering
    nclusroot = sqrt((double)numElementsOfBoundaryMesh_)*0.9;
    if(nclusroot < d_min) nclusroot = d_min;
    clu_l_min = sqrt(d_min)*ele_leng;

	// compute the cluster length on level 0
	if(ClusEdgL0_ <= 0.0) // the cluster length on level 0 is not input by the user
	{
		double clus_nels = (double)numElementsOfBoundaryMesh_/nclusroot;
		ClusEdgL1 = sqrt(clus_nels)*ele_leng;
	}
	else // cluster length of the level 0 is input by the user
	{
		ClusEdgL1 = ClusEdgL0_;
		if(ClusEdgL1 < clu_l_min)
		{
			NC_Error_Exit_2(NCout,
				"Inputed value of cluster diameter for level 0 is too small!",
				"Inputed value of cluster diameter for level 0 = ", ClusEdgL1,
				"Smallst allowable value of cluster diameter for level 0 = ", clu_l_min);
		}
	}

	d_levs = log10(ClusEdgL1/clu_l_min)/log10(2.0);
	n_levs = (int)d_levs;
	if(d_levs - (double)n_levs > 0.5) n_levs++;
	n_levs++;
	if(n_levs <= 0) n_levs = 1;

	return(n_levs);
}

// generate the cluster tree for MLFMBEM
void NC_GenerateClusterTreeMLFMM
(
	ofstream& NCout
)
{
	int i, ncl;

	// generate the array of levels of the tree
	clulevarry = new ClusterLev[numClusterLevels_];

	// compute the clusters at each level
	for(i=0; i<numClusterLevels_; i++)
	{
		ncl = NC_GenerateClustersAtLevelMLFMM(NCout, i);
		NC_CheckClusters(NCout, i, clulevarry[i].nClustOLv, clulevarry[i].ClastArLv);
	}
}

// generate clusters at given level
int NC_GenerateClustersAtLevelMLFMM
(
	ofstream& NCout,
	const int& nu_lev	// number of the level
)
{
	int i, j;
	int nbegrp_facl, melpergrp=0;
	int num_clus = 0; // number of clusters

	// number of elements of each cluster and numbers of these elements
	Vector<int> nel_clus(numElementsOfBoundaryMesh_), nuel_clus(numElementsOfBoundaryMesh_), nfath_clus(numElementsOfBoundaryMesh_);

	if(nu_lev == 0) // the TRUNK level
	{
		// number of BE groups
		int nbegrp_facl = 0;
		for(i=0; i<numElementGroups_; i++)
		{
			if(propertyOfGroup[i] == 0) nbegrp_facl++;
		}

		// numbers of BE groups
		Vector<int> nubegrp(nbegrp_facl);
		j = 0;
		for(i=0; i<numElementGroups_; i++)
		{
			if(propertyOfGroup[i] == 0) nubegrp[j++] = i;
		}

		// maximal number of elements of an element group
		for(i=0; i<nbegrp_facl; i++)
			if(numberOfElementsInGroup[nubegrp[i]] > melpergrp) melpergrp = numberOfElementsInGroup[nubegrp[i]];
		Vector<int> nuelbegrp(melpergrp);

		// generate clusters
		NC_GenerateClustersFMM(NCout, num_clus, nel_clus, nuel_clus, nbegrp_facl,
					 nubegrp, nuelbegrp, nu_lev, nfath_clus);
	}
	else // the BRANCH and LEAF levels
	{
		// number of original clusters at the father level
		nbegrp_facl = clulevarry[nu_lev - 1].nClustOLv;

		// does not be used
		Vector<int> nubegrp(0);

		// maximal number of elements of a father cluster
		for(i=0; i<nbegrp_facl; i++)
		{
			j = clulevarry[nu_lev - 1].ClastArLv[i].NumOfEl;
			if(j > melpergrp) melpergrp = j;
		}
		Vector<int> nuelbegrp(melpergrp);

		// generate clusters
		NC_GenerateClustersFMM(NCout, num_clus, nel_clus, nuel_clus, nbegrp_facl,
					 nubegrp, nuelbegrp, nu_lev, nfath_clus);
	} // end of BRANCH and LEAF level

	// generate the cluster structure at a given level
	NC_GenerateClusterArrayAtLevelMLFMM(NCout, nel_clus, nuel_clus, num_clus, nu_lev, nfath_clus);

	// reflection of the clusters with respect to the symmetry planes
	if(numReflectionsOfElements_ > 1) NC_ClusterReflectionsAtLevelMLFMM(NCout, nu_lev);

	return(num_clus);
}

// generate the cluster structure at a given level
void NC_GenerateClusterArrayAtLevelMLFMM
(
	ofstream& NCout,
	Vector<int>& nel_clus,		//I: number of elements in each cluster
	Vector<int>& nuel_clus,		//I: numbers of elements in the clusters
	int& num_clus,				//I: number of clusters
	const int& nu_lev,			//I: number of the level
	Vector<int>& nfath_clus		//I: number of the father clusters
)
{
	int i, j, k, k1, l_sum, l1, m, ie0, ind0;
	int nel_max = 0, n_nod;
	double wk0;
	double dareadis, dispo;
	Vector<double> crcent(NDIM);

	// maximum number of elements in a cluster
	for(i=0; i<num_clus; i++) if(nel_clus[i] > nel_max) nel_max = nel_clus[i];

	// vector to store nodal numbers of each cluster
	Vector<int> nu_nod_cl(nel_max*4);

	// number of original clusters and number of all clusters at the level
	clulevarry[nu_lev].nClustOLv = num_clus;
	clulevarry[nu_lev].nClustSLv = num_clus*numReflectionsOfElements_;

	// create the cluster structure array at the level
	clulevarry[nu_lev].ClastArLv = new ElCluster[num_clus*numReflectionsOfElements_];

	// loop over clusters
	l_sum = 0;

	for(i=0; i<num_clus; i++)
	{
		// number of elements of the cluster
		clulevarry[nu_lev].ClastArLv[i].NumOfDOFs =
			clulevarry[nu_lev].ClastArLv[i].NumOfEl = nel_clus[i];

		// numbers of the elements of the cluster
		clulevarry[nu_lev].ClastArLv[i].NumsOfEl = new int[nel_clus[i]];
		for(j=0; j<nel_clus[i]; j++) {
		  // for debugging
		  clulevarry[nu_lev].ClastArLv[i].NumsOfEl[j] = nuel_clus[l_sum++];
		}

		// number of the element group to which the cluster belongs
		ie0 = clulevarry[nu_lev].ClastArLv[i].NumsOfEl[0];
		clulevarry[nu_lev].ClastArLv[i].NuElGr = listElementsElementGroup[ie0];

		// surface or middle face elements ( = 0: surface els; = 1: middle face els)
		clulevarry[nu_lev].ClastArLv[i].listElementPropertyEl = listElementProperty[ie0];

		// if all elements of the cluster are of the same number of nodes
		clulevarry[nu_lev].ClastArLv[i].IfMonoEl = true;
		ind0 = listNumberNodesPerElement[ie0];
		for(j=1; j<nel_clus[i]; j++)
		{
			if(listNumberNodesPerElement[clulevarry[nu_lev].ClastArLv[i].NumsOfEl[j]] != ind0)
			{
				clulevarry[nu_lev].ClastArLv[i].IfMonoEl = false;
				break;
			}
		}

		// if admittance boundary condition are prescribed
		if(ibval[ie0] == 2 || ibval[ie0] == 4)
			clulevarry[nu_lev].ClastArLv[i].IfAdmiBc = true;
		else clulevarry[nu_lev].ClastArLv[i].IfAdmiBc = false;

		// number of the unknown DOFs of the cluster
        clulevarry[nu_lev].ClastArLv[i].NDOFsPeEl = 1;

		// compute number and numbers of the nodes of the current cluster
		n_nod = 0;
		for(j=0; j<nel_clus[i]; j++)
		{
			k1 = clulevarry[nu_lev].ClastArLv[i].NumsOfEl[j];
			for(k=0; k<listNumberNodesPerElement[k1]; k++)
			{
				l1 = elementsConnectivity[k1][k];
				for(m=0; m<n_nod; m++) if(nu_nod_cl[m] == l1) goto LbclustAray2;
				nu_nod_cl[n_nod++] = l1;
LbclustAray2: continue;
			}
		}

		// compute the coordinates of the center of the cluster
		crcent = 0.0;
		for(j=0; j<n_nod; j++)
		{
			for(k=0; k<NDIM; k++) crcent[k] += nodesCoordinates[nu_nod_cl[j]][k];
		}
		for(k=0; k<NDIM; k++)
			clulevarry[nu_lev].ClastArLv[i].CoorCent[k] = crcent[k]/(double)n_nod;

		// compute the radius of the cluster
		dareadis = 0.0;
		for(j=0; j<n_nod; j++)
		{
			dispo = Dispoi_dim3_(clulevarry[nu_lev].ClastArLv[i].CoorCent,
				nodesCoordinates[nu_nod_cl[j]]);
			if(dispo > dareadis) dareadis = dispo;
		}
		clulevarry[nu_lev].ClastArLv[i].RadiClus = dareadis;

	} // end of loop I

	// compute the maximum and minimum cluster radius
	avgClusterRadiusBE_ = maxClusterRadiusBE_ = 0;
	minClusterRadiusBE_ = 1.0e40;
	for(i=0; i<num_clus; i++)
	{
		if(clulevarry[nu_lev].ClastArLv[i].RadiClus > maxClusterRadiusBE_)
			maxClusterRadiusBE_ = clulevarry[nu_lev].ClastArLv[i].RadiClus;
		if(clulevarry[nu_lev].ClastArLv[i].RadiClus < minClusterRadiusBE_)
			minClusterRadiusBE_ = clulevarry[nu_lev].ClastArLv[i].RadiClus;
		avgClusterRadiusBE_ += clulevarry[nu_lev].ClastArLv[i].RadiClus;
	}
	avgClusterRadiusBE_ /= (double)num_clus;
	clulevarry[nu_lev].RadiMaxLv = maxClusterRadiusBE_;
	clulevarry[nu_lev].RadiAveLv = avgClusterRadiusBE_;
	clulevarry[nu_lev].RadiMinLv = minClusterRadiusBE_;

	// define the near and far field clusters
	if(numReflectionsOfElements_ == 1)
	{
		Matrix<bool> ifarclus(num_clus, num_clus, false);
		for(i=0; i<num_clus; i++) for(j=0; j<i; j++)
		{
			dispo = Dispoi_dim3_(clulevarry[nu_lev].ClastArLv[i].CoorCent,
				clulevarry[nu_lev].ClastArLv[j].CoorCent);

			wk0 = clulevarry[nu_lev].ClastArLv[i].RadiClus +
				clulevarry[nu_lev].ClastArLv[j].RadiClus;

			if(dispo > farFieldClusterFactor_*wk0) {
				if(dispo > minClusterDistance_) {
					ifarclus(j, i) = ifarclus(i, j) = true;
				} else {
					n_Pair_NeaF[nu_lev]++;
				}
				n_Pair_FarD[nu_lev]++;
			}
		}

		int nnea, nfar;
		for(i=0; i<num_clus; i++)
		{
			clulevarry[nu_lev].ClastArLv[i].NumNeaClus =
				clulevarry[nu_lev].ClastArLv[i].NumFarClus = 0;
			for(j=0; j<num_clus; j++)
				if(ifarclus(i, j)) clulevarry[nu_lev].ClastArLv[i].NumFarClus++;
				else clulevarry[nu_lev].ClastArLv[i].NumNeaClus++;

			clulevarry[nu_lev].ClastArLv[i].NumsNeaClus =
				new int[clulevarry[nu_lev].ClastArLv[i].NumNeaClus];
			clulevarry[nu_lev].ClastArLv[i].NumsFarClus =
				new int[clulevarry[nu_lev].ClastArLv[i].NumFarClus];
			nnea = nfar = 0;
			for(j=0; j<num_clus; j++) if(ifarclus(i, j))
				clulevarry[nu_lev].ClastArLv[i].NumsFarClus[nfar++] = j;
			else clulevarry[nu_lev].ClastArLv[i].NumsNeaClus[nnea++] = j;
		}
	} // end of NREFL == 1

	// create reflection informations
	for(i=0; i<num_clus; i++)
	{
		clulevarry[nu_lev].ClastArLv[i].OriClust = i;
		clulevarry[nu_lev].ClastArLv[i].nuref = 0;
		clulevarry[nu_lev].ClastArLv[i].rffac = 1;
		clulevarry[nu_lev].ClastArLv[i].ifmirro = false;
		for(j=0; j<NDIM; j++) clulevarry[nu_lev].ClastArLv[i].ifrfdi[j] = false;
	}

	// register the father-son relation
	if(nu_lev > 0)
	{
		int nfacl, nl_f = nu_lev - 1;

		for(i=0; i<clulevarry[nl_f].nClustOLv; i++)
			clulevarry[nl_f].ClastArLv[i].n_Son = 0;

		for(i=0; i<num_clus; i++)
		{
			nfacl = nfath_clus[i];
			clulevarry[nu_lev].ClastArLv[i].nuFather = nfacl;
			clulevarry[nl_f].ClastArLv[nfacl].nuSon
				[clulevarry[nl_f].ClastArLv[nfacl].n_Son++] = i;
		}
	}

	// level number
	for(i=0; i<num_clus; i++) clulevarry[nu_lev].ClastArLv[i].nuLev = nu_lev;
}

// generate the clusters of the internal points (nodes of the evaluation mesh)
int NC_GenerateClustersEvalGrid
(
	ofstream& NCout
 )
{
    int i, j, k, l, ind, iclus, nclus, l1, i1;
    double d0, d1, radi = 0.0;
    Vector<double> xyz_max(NDIM, -1.0e20), xyz_min(NDIM, 1.0e20), dif_xyz(NDIM);
    Vector<int> nuinnode(numNodesOfEvaluationMesh_), ndiv_xyz(NDIM);
    Matrix<int> nuindxyz(numNodesOfEvaluationMesh_, NDIM, -1);

    // compute the vector of the internal node number of each internal point
    j = 0;
    for(i=0; i<numNodes_; i++) if(isNodeMeshOrEvalNode[i] == -1) nuinnode[j++] = i;

    // compute the maximum and minimum coordinates of the internal points
    for(i=0; i<numNodesOfEvaluationMesh_; i++)
    {
        ind = nuinnode[i];
        for(j=0; j<NDIM; j++)
        {
            if(nodesCoordinates[ind][j] > xyz_max[j]) xyz_max[j] = nodesCoordinates[ind][j];
            if(nodesCoordinates[ind][j] < xyz_min[j]) xyz_min[j] = nodesCoordinates[ind][j];
        }
    }
    for(j=0; j<NDIM; j++)
    {
        xyz_max[j] += delta_;
        xyz_min[j] -= delta_;
    }
    for(j=0; j<NDIM; j++) dif_xyz[j] = xyz_max[j] - xyz_min[j];

    // compute number of divisions in each coordinate direction
    for(j=0; j<NDIM; j++)
    {
        d0 = dif_xyz[j]/(ClusEdgL1);
        i = (int)(d0);
        if(i == 0) i = 1;
        else if(d0 - (double)(i) > 0.5) i++;
        ndiv_xyz[j] = i;
        dif_xyz[j] /= (double)(i);
    }


    nclus = ndiv_xyz[0]*ndiv_xyz[1]*ndiv_xyz[2];
    if( (nclus > numNodesOfEvaluationMesh_) ) {
      /* kreiza: 17.08.22 
	 this will happen if the BEM mesh and the Eval mesh 
	 differ in size very much. Let's try a new approach, every 
	 evaluation node is its own clusters. This is not really ellegant, 
	 but it  saves trouble in the postprocessing */
      //cout << "Warning: This feature is very experimental, and has not been tested yet. There is a problem, that the evaluation grid is much bigger then the BEM grid, are you sure that you have the correct units for both of them?\n";
      iclus = numNodesOfEvaluationMesh_;
      if(iclus > 0) ipcluarry = new IpCluster[iclus];
      for (i = 0; i < iclus; i++) {
	ipcluarry[i].NumOfIps = 1;
	ipcluarry[i].NumsOfIps = new int[1];
	ipcluarry[i].NumsOfIps[0] = i;
	ind = nuinnode[i];
	for (j = 0; j < NDIM; j++) {
	  ipcluarry[i].CoorCent[j] = nodesCoordinates[ind][j];
	}
	ipcluarry[i].RadiClus = 0.0;
      }
    }
    else {
      // really cluster the evalgrid
      // compute the numbers of the intervals in the three coordinate directions for each internal point, in which it is located
      for(j=0; j<NDIM; j++)
	{
	  d0 = xyz_min[j];
	  for(k=0; k<ndiv_xyz[j]; k++)
	    {
	      d1 = d0 + dif_xyz[j];
	      for(i=0; i<numNodesOfEvaluationMesh_; i++)
		{
		  ind = nuinnode[i];
		  if(nuindxyz(i, j) == -1 && nodesCoordinates[ind][j] >= d0 && nodesCoordinates[ind][j] < d1)
                    nuindxyz(i, j) = k;
		}
	      d0 = d1;
	    }
	}
      
      // check the results of the interval nubers
      for(i=0; i<numNodesOfEvaluationMesh_; i++) for(j=0; j<NDIM; j++){
	  if(nuindxyz(i, j) == -1) {
            NC_Error_Exit_2(NCout, "The interval number of a internal point is not found!",
                            "number of the internal point = ", extNumbersOfNodes[nuinnode[i]],
                            "number of the direction = ", j);
	  }
	}
      
      // compute (1) number of the internal points in each cluster (nip_clus)
      //         (2) the sequential number of each point in the clusters (nuip_clus)
      
      Vector<int> nip_clus(nclus, 0), nuip_clus(numNodesOfEvaluationMesh_);
      
      iclus = 0;
      i1 = 0;
      for(i=0; i<ndiv_xyz[0]; i++) for(j=0; j<ndiv_xyz[1]; j++) for(k=0; k<ndiv_xyz[2]; k++)
								  {
								    l1 = 0;
								    for(l=0; l<numNodesOfEvaluationMesh_; l++)
								      if(nuindxyz(l, 0) == i && nuindxyz(l, 1) == j && nuindxyz(l, 2) == k)
									{
									  l1++;
									  nuip_clus[i1++] = l;
									}
								    nip_clus[iclus] = l1;
								    
								    iclus++;
								  }
      
      // compute the number of clusters
      iclus = 0;
      for(j=0; j<nclus; j++) if(nip_clus[j] > 0) iclus++;
      
      // generate the cluster array of the internal points
      if(iclus > 0) ipcluarry = new IpCluster[iclus];
      
      l1 = l = 0;
      minClusterRadiusRM_ = 1.0e15;
      avgClusterRadiusRM_ = maxClusterRadiusRM_ = 0.0;
      for(i=0; i<nclus; i++)
	{
	  // if there are not points in the cluster
	  if(nip_clus[i] == 0) continue;
	  
	  // number of innternal points in the cluster
	  ipcluarry[l1].NumOfIps = nip_clus[i];
	  
	  // the sequential numbers of internal points in the cluster
	  ipcluarry[l1].NumsOfIps = new int[nip_clus[i]];
	  for(j=0; j<nip_clus[i]; j++) ipcluarry[l1].NumsOfIps[j] = nuip_clus[l++];
	  
	  // coordinates of the center of the cluster
	  for(j=0; j<NDIM; j++) xyz_max[j] = 0.0;
	  for(j=0; j<ipcluarry[l1].NumOfIps; j++)
	    {
	      ind = nuinnode[ipcluarry[l1].NumsOfIps[j]];
	      for(k=0; k<NDIM; k++) xyz_max[k] += nodesCoordinates[ind][k];
	    }
	  for(j=0; j<NDIM; j++) ipcluarry[l1].CoorCent[j] =
				  xyz_max[j]/(double)(ipcluarry[l1].NumOfIps);
	  
	  // radius of the clusters
	  radi = 0.0;
	  for(j=0; j<ipcluarry[l1].NumOfIps; j++)
	    {
	      ind = nuinnode[ipcluarry[l1].NumsOfIps[j]];
	      d0 = 0.0;
	      for(k=0; k<NDIM; k++)
		{
		  d1 = nodesCoordinates[ind][k] - ipcluarry[l1].CoorCent[k];
		  d0 += d1*d1;
		}
	      d0 = sqrt(d0);
	      if(d0 > radi) radi = d0;
	    }
	  ipcluarry[l1].RadiClus = radi;
	  
	  // compute the maximum, minimum and average radius of the clusters
	  if(radi > maxClusterRadiusRM_) maxClusterRadiusRM_ = radi;
	  if(radi < minClusterRadiusRM_) minClusterRadiusRM_ = radi;
	  avgClusterRadiusRM_ += radi;
	  l1++;
	}
      avgClusterRadiusRM_ /= (double)l1;
      
      if(methodFMM_ >= 2) {radi = clulevarry[0].RadiMaxLv;}
      else if(methodFMM_ == 1) {radi = maxClusterRadiusBE_;}
      if(maxClusterRadiusRM_ > radi)
	{
	  NC_Error_Warning_2(NCout, "Maximum radius of evaluation clusters >  that of BE clusters!",
			     "Maximum radius of evaluation clusters = ", maxClusterRadiusRM_,
			     "Maximum radius of BE clusters = ", radi);
	}
    } // that should be the if statement about the regular clusters
    return(iclus);
}

// reflection of the element clusters of a given level with respect to the symmetry planes
void NC_ClusterReflectionsAtLevelMLFMM
(
	ofstream& NCout,
	const int& nu_lev			//I: number of the level
)
{
	int i, j, i0 = 0, icl, nucl, kref, nrdir = 0, nel;
	double dispo, wk0;
	int nrefi[NDIM];
	int nrefdi[] = {0, 0, 1, 0, 2, 0, 1, 0};

	// generate the reflected clusters
	nucl = clulevarry[nu_lev].nClustOLv; // number of the generated cluster (updated)
	for(icl=0; icl<clulevarry[nu_lev].nClustOLv; icl++)
	{
		for(i=0; i<NDIM; i++) nrefi[i] = 0;
		for(kref=1; kref<numReflectionsOfElements_; kref++)
		{
			// direction of the current reflection
			switch(numReflectionsOfElements_)
			{
			case 8:  // there are 3 symmetric planes
				nrdir = nrefdi[kref];
				nrefi[nrdir]++;
				break;
			case 4:  // there are 2 symmetric planes
				for(i=0; i<NDIM; i++) if(listSymmetryPlanes[i] == 0)
				{
					i0 = i + 1;
					break;
				}
				nrdir = i0 + nrefdi[kref];
				if(nrdir >= NDIM) nrdir -= NDIM;
				nrefi[nrdir]++;
				break;
			case 2:  // there is 1 symmetric plane
				for(i=0; i<NDIM; i++) if(listSymmetryPlanes[i] != 0)
				{
					nrdir = i;
					break;
				}
				nrefi[nrdir]++;
			} // end of SWITCH

			// parameters of the reflected cluster
			clulevarry[nu_lev].ClastArLv[nucl].NuElGr = clulevarry[nu_lev].ClastArLv[icl].NuElGr;
			clulevarry[nu_lev].ClastArLv[nucl].listElementPropertyEl = clulevarry[nu_lev].ClastArLv[icl].listElementPropertyEl;
			clulevarry[nu_lev].ClastArLv[nucl].IfMonoEl = clulevarry[nu_lev].ClastArLv[icl].IfMonoEl;
			clulevarry[nu_lev].ClastArLv[nucl].IfAdmiBc = clulevarry[nu_lev].ClastArLv[icl].IfAdmiBc;
			clulevarry[nu_lev].ClastArLv[nucl].NumOfDOFs = clulevarry[nu_lev].ClastArLv[icl].NumOfDOFs;
			clulevarry[nu_lev].ClastArLv[nucl].NDOFsPeEl = clulevarry[nu_lev].ClastArLv[icl].NDOFsPeEl;

			nel = clulevarry[nu_lev].ClastArLv[nucl].NumOfEl =
				clulevarry[nu_lev].ClastArLv[icl].NumOfEl;
			clulevarry[nu_lev].ClastArLv[nucl].NumsOfEl = new int[nel];
			for(i=0; i<nel; i++) clulevarry[nu_lev].ClastArLv[nucl].NumsOfEl[i] =
				clulevarry[nu_lev].ClastArLv[icl].NumsOfEl[i];

			clulevarry[nu_lev].ClastArLv[nucl].RadiClus = clulevarry[nu_lev].ClastArLv[icl].RadiClus;

			for(i=0; i<NDIM; i++)
			{
				if(nrefi[i]%2) clulevarry[nu_lev].ClastArLv[nucl].ifrfdi[i] = true;
				else clulevarry[nu_lev].ClastArLv[nucl].ifrfdi[i] = false;
				if(clulevarry[nu_lev].ClastArLv[nucl].ifrfdi[i])
					clulevarry[nu_lev].ClastArLv[nucl].CoorCent[i] = 2.0*listSymmetryPlanesCoordinates[i] -
						clulevarry[nu_lev].ClastArLv[icl].CoorCent[i];
				else
					clulevarry[nu_lev].ClastArLv[nucl].CoorCent[i] =
						clulevarry[nu_lev].ClastArLv[icl].CoorCent[i];
			}

			clulevarry[nu_lev].ClastArLv[nucl].OriClust = icl;
			clulevarry[nu_lev].ClastArLv[nucl].nuref = kref;
			if(kref%2) clulevarry[nu_lev].ClastArLv[nucl].ifmirro = true;
			else clulevarry[nu_lev].ClastArLv[nucl].ifmirro = false;

			clulevarry[nu_lev].ClastArLv[nucl].rffac = 1;
			for(i=0; i<NDIM; i++) {
				if(clulevarry[nu_lev].ClastArLv[nucl].ifrfdi[i] && listSymmetryPlanes[i] < 0)
					clulevarry[nu_lev].ClastArLv[nucl].rffac *= -1;
			}

			nucl++;
		} // end of KREF
	} // end of loop ICL

	// compute the far and near field clusters of each cluster
	int nclusLev = clulevarry[nu_lev].nClustSLv;
	Matrix<bool> ifarclus(nclusLev, nclusLev, false);

	for(i=0; i<nclusLev; i++) for(j=0; j<i; j++)
	{
		dispo = Dispoi_dim3_(clulevarry[nu_lev].ClastArLv[i].CoorCent,
			clulevarry[nu_lev].ClastArLv[j].CoorCent);
		wk0 = clulevarry[nu_lev].ClastArLv[i].RadiClus +
			clulevarry[nu_lev].ClastArLv[j].RadiClus;

		if(dispo > farFieldClusterFactor_*wk0) {
			if(dispo > minClusterDistance_) {
				ifarclus(j, i) = ifarclus(i, j) = true;
			} else {
				n_Pair_NeaF[nu_lev]++;
			}
			n_Pair_FarD[nu_lev]++;
		}
	}

	int nnea, nfar;
	for(i=0; i<nclusLev; i++)
	{
		clulevarry[nu_lev].ClastArLv[i].NumNeaClus =
			clulevarry[nu_lev].ClastArLv[i].NumFarClus = 0;
		for(j=0; j<nclusLev; j++)
			if(ifarclus(i, j)) clulevarry[nu_lev].ClastArLv[i].NumFarClus++;
			else clulevarry[nu_lev].ClastArLv[i].NumNeaClus++;

		clulevarry[nu_lev].ClastArLv[i].NumsNeaClus = new int[clulevarry[nu_lev].ClastArLv[i].NumNeaClus];
		clulevarry[nu_lev].ClastArLv[i].NumsFarClus = new int[clulevarry[nu_lev].ClastArLv[i].NumFarClus];
		nnea = nfar = 0;
		for(j=0; j<nclusLev; j++)
			if(ifarclus(i, j)) clulevarry[nu_lev].ClastArLv[i].NumsFarClus[nfar++] = j;
			else clulevarry[nu_lev].ClastArLv[i].NumsNeaClus[nnea++] = j;
	}


	// register the father-son relation
	if(nu_lev > 0)
	{
		int nfacl, nrefres = numReflectionsOfElements_ - 1, nclofa = clulevarry[nu_lev - 1].nClustOLv,
			nl_f = nu_lev - 1;

		nucl = clulevarry[nu_lev].nClustOLv;
		for(icl=0; icl<clulevarry[nu_lev].nClustOLv; icl++)
		{
			nfacl = nclofa + nrefres*clulevarry[nu_lev].ClastArLv[icl].nuFather;
			for(kref=1; kref<numReflectionsOfElements_; kref++)
			{
				clulevarry[nu_lev].ClastArLv[nucl++].nuFather = nfacl++;
			}
		}

		for(i=clulevarry[nl_f].nClustOLv; i<clulevarry[nl_f].nClustSLv; i++)
			clulevarry[nl_f].ClastArLv[i].n_Son = 0;

		for(i=clulevarry[nu_lev].nClustOLv; i<clulevarry[nu_lev].nClustSLv; i++)
		{
			nfacl = clulevarry[nu_lev].ClastArLv[i].nuFather;
			clulevarry[nl_f].ClastArLv[nfacl].nuSon
				[clulevarry[nl_f].ClastArLv[nfacl].n_Son++] = i;
		}
	}

	// level number
	for(i=clulevarry[nu_lev].nClustOLv; i<clulevarry[nu_lev].nClustSLv; i++)
		clulevarry[nu_lev].ClastArLv[i].nuLev = nu_lev;
}

// compute the auxiliary arrays used for storing the sparse far field matrices for the SLFMBEM (called by NC_ControlProgram)
void NC_AllocateSDTmtxsSLFMM
(
	ofstream& NCout
 )
{
    int i, j, k;
    int nnonzerdmtx, nnonzertmtx, nnonzersmtx;
    int inoze, nnoze, k1, idfs;
    int nnze, irowd, nrowdmtx;
    int ij, kcol, itr;

    // generate the T-matrix
    nnonzertmtx = 0; // number of the nonzeros
    for(i=0; i<numOriginalReflectedClusters_; i++) nnonzertmtx += clustarry[i].NumOfDOFs;
    nnonzertmtx *= numIntegrationPointsUnitSphere_;

    // compute the address arrays (for storing the column numbers of the nonzeros
    // and the global address of the first nonzero in every row)
    jcoltmtx = new int[nnonzertmtx];
    irowtmtx = new int[numIntegrationPointsUnitSphere_*numOriginalReflectedClusters_ + 1];

    inoze = nnoze = 0;
    for(i=0; i<numOriginalReflectedClusters_; i++)
    {
        idfs = clustarry[i].NDOFsPeEl;
        for(j=0; j<numIntegrationPointsUnitSphere_; j++)
        {
            irowtmtx[inoze++] = nnoze;
            for(k=0; k<clustarry[i].NumOfEl; k++)
            {
                k1 = clustarry[i].NumsOfEl[k];
                jcoltmtx[nnoze++] = jelist[k1][0];
                if(idfs == 2) jcoltmtx[nnoze++] = jelist[k1][1];
            }
        }
    }
    irowtmtx[inoze] = nnoze;
    // end of generate the T-matrix

    // generate the D-matrix
    maxRowNumberD_ = nrowdmtx = numOriginalClusters_*numIntegrationPointsUnitSphere_; // number of rows
    nnonzerdmtx = 0; // number of nonzeros
    for(i=0; i<numOriginalClusters_; i++)
    {
        nnonzerdmtx += clustarry[i].NumFarClus;
    }
    nnonzerdmtx *= numIntegrationPointsUnitSphere_;

    // compute the address arrays
    dmtxlev[0].jcolDmxLv = new int[nnonzerdmtx];
    dmtxlev[0].irowDmxLv = new int[nrowdmtx + 1];

    irowd = nnze = 0;
    for(i=0; i<numOriginalClusters_; i++)
    {
        for(k=0; k<numIntegrationPointsUnitSphere_; k++)
        {
            dmtxlev[0].irowDmxLv[irowd++] =  nnze;
            for(j=0; j<clustarry[i].NumFarClus; j++)
                dmtxlev[0].jcolDmxLv[nnze++] = clustarry[i].NumsFarClus[j]*numIntegrationPointsUnitSphere_ + k;
        }
    }
    dmtxlev[0].irowDmxLv[irowd] =  nnze;

    // compute the number of rows and the number of nonzero entries of the matrx D
    dmtxlev[0].nRowsD = numOriginalClusters_*numIntegrationPointsUnitSphere_;
    dmtxlev[0].nEntriesD = dmtxlev[0].irowDmxLv[numOriginalClusters_*numIntegrationPointsUnitSphere_];
    // end of generate the D-matrix

    // generate the S-matrix
    // compute the array for storing cluster number of each row
    Vector<int> nuclusrow(numRowsOfCoefficientMatrix_, -1);
    for(i=0; i<numOriginalClusters_; i++)
    {
        for(j=0; j<clustarry[i].NumOfEl; j++)
        {
            ij = clustarry[i].NumsOfEl[j];
            nuclusrow[jelist[ij][0]] = i;
            if(clustarry[i].NDOFsPeEl == 2) nuclusrow[jelist[ij][1]] = i;
        }
    }

    // compute number of nonzeros
    nnonzersmtx = 0;
    for(i=0; i<numRowsOfCoefficientMatrix_; i++) if(nuclusrow[i] >= 0) nnonzersmtx += numIntegrationPointsUnitSphere_;

    // compute the address arrays
    jcolsmtx = new int[nnonzersmtx];
    irowsmtx = new int[numRowsOfCoefficientMatrix_ + 1];

    // compute the address arrays
    itr = 0;
    for(i=0; i<numRowsOfCoefficientMatrix_; i++)
    {
        irowsmtx[i] = itr;
        if(nuclusrow[i] >= 0)
        {
            kcol = numIntegrationPointsUnitSphere_*nuclusrow[i];
            for(j=0; j<numIntegrationPointsUnitSphere_; j++) jcolsmtx[itr++] = kcol++;
        }
    }
    irowsmtx[numRowsOfCoefficientMatrix_] = itr;
    // end of generate the S-matrix

    // write the results into the output file
    int sum_nonzero = irownea[numRowsOfCoefficientMatrix_] + irowsmtx[numRowsOfCoefficientMatrix_] +
    dmtxlev[0].nEntriesD + irowtmtx[numOriginalReflectedClusters_*numIntegrationPointsUnitSphere_];

    NCout << endl;
    NCout << "Nonzeros of the near field matrix ..................... = " <<
    irownea[numRowsOfCoefficientMatrix_] << endl;
    NCout << "Nonzeros of the T-matrix .............................. = " <<
    irowtmtx[numOriginalReflectedClusters_*numIntegrationPointsUnitSphere_] << endl;
    NCout << "Nonzeros of the D-matrix .............................. = " <<
    dmtxlev[0].nEntriesD << endl;
    NCout << "Nonzeros of the S-matrix .............................. = " <<
    irowsmtx[numRowsOfCoefficientMatrix_] << endl;
    NCout << "Sum of nonzeros of the FMBEM matrices ................. = " <<
    sum_nonzero << endl;
}

// create the far field matrices for ML-FMBEM (called by NC_ControlProgram)
void NC_AllocateSDTmtxsMLFMM
(
	ofstream& NCout
)
{
	// generate the T-matrix
	NC_AllocateTmtxMLFMM(NCout);

	// generate the D-matrices
	NC_AllocateDmtxMLFMM(NCout);

	// generate the S-matrix
	NC_AllocateSmtxMLFMM(NCout);

	// output the results
    int i, sum_d = 0, sum_t = 0, sum_s = 0;

    for(i=0; i<numClusterLevels_; i++)
    {
        sum_t += tmtxlev[i].nEntriesT;
        sum_d += dmtxlev[i].nEntriesD;
        sum_s += smtxlev[i].nEntriesS;
    }

    NCout << endl;

    NCout << "Nonzeros of the near field matrix ..................... = " <<
        numNonZeroEntries_ << endl;

    NCout << "Nonzeros of the T-matrix .............................. = " <<
        sum_t << endl;
    for(i=0; i<numClusterLevels_; i++)
    {
        NCout << "    Nonzeros of T_lev_" << i <<
            " ............................... = " << tmtxlev[i].nEntriesT << endl;
    }

    NCout << "Nonzeros of the D-matrices ............................ = " <<
        sum_d << endl;
    for(i=0; i<numClusterLevels_; i++)
    {
        NCout << "    Nonzeros of D_lev_" << i <<
            " ............................... = " << dmtxlev[i].nEntriesD << endl;
    }

    NCout << "Nonzeros of the S-matrix .............................. = " <<
        sum_s << endl;
    for(i=0; i<numClusterLevels_; i++)
    {
        NCout << "    Nonzeros of S_lev_" << i <<
            " ............................... = " << smtxlev[i].nEntriesS << endl;
    }

    NCout << "Sum of nonzeros of the FMBEM matrices ................. = " <<
        numNonZeroEntries_ + sum_t + sum_d + sum_s << endl;
}

// generate the T-matrices at all levels
void NC_AllocateTmtxMLFMM
(
	ofstream& NCout
)
{
	int nlv, inoze, nnoze, i, j, k, k1, idfs;

	// loop over the levels
	for(nlv=0; nlv<numClusterLevels_; nlv++)
	{
		// compute number of rows of the T-matrix and the T-vector
		j = 0;
		for(i=0; i<clulevarry[nlv].nClustSLv; i++)
		{
			j += clulevarry[nlv].ClastArLv[i].NumOfDOFs;

		}
		tmtxlev[nlv].nEntriesT = j*clulevarry[nlv].nPoinSpheLv;
		tmtxlev[nlv].nRowsT = clulevarry[nlv].nClustSLv*clulevarry[nlv].nPoinSpheLv;

		tmtxlev[nlv].jcolTmxLv = new int[tmtxlev[nlv].nEntriesT];
		tmtxlev[nlv].irowTmxLv = new int[tmtxlev[nlv].nRowsT + 1];

		inoze = nnoze = 0;
		for(i=0; i<clulevarry[nlv].nClustSLv; i++)
		{
			idfs = clulevarry[nlv].ClastArLv[i].NDOFsPeEl;
			for(j=0; j<clulevarry[nlv].nPoinSpheLv; j++)
			{
				tmtxlev[nlv].irowTmxLv[inoze++] = nnoze;
				for(k=0; k<clulevarry[nlv].ClastArLv[i].NumOfEl; k++)
				{
					k1 = clulevarry[nlv].ClastArLv[i].NumsOfEl[k];
					tmtxlev[nlv].jcolTmxLv[nnoze++] = jelist[k1][0];
					if(idfs == 2) tmtxlev[nlv].jcolTmxLv[nnoze++] = jelist[k1][1];
				}
			}
		}
		tmtxlev[nlv].irowTmxLv[inoze] = nnoze;
	} // end of loop NLV
}

// generate the D-matrices at all levels
void NC_AllocateDmtxMLFMM
(
	ofstream& NCout
 )
{
    int i, j, k, nnze, irowd, nlv, nrosdmtx, nnonzerdmtx;

    Vector<int> nclufar(numOriginalClusters_);
    Matrix<int> nuclfar(numOriginalClusters_, numOriginalReflectedClusters_);

    // loop over levels
    maxRowNumberD_ = 0;
    for(nlv=0; nlv<numClusterLevels_; nlv++)
    {
        // compute the relevant far clusters for the current level
        if(nlv == 0) // the TRUNK level
        {
            for(i=0; i<clulevarry[nlv].nClustOLv; i++)
            {
                nclufar[i] = clulevarry[nlv].ClastArLv[i].NumFarClus;
                for(j=0; j<nclufar[i]; j++)
                    nuclfar(i,j) = clulevarry[nlv].ClastArLv[i].NumsFarClus[j];
            }
        }
        else // the BRANCH and LEAF levels
        {
            int nlf = nlv-1, ifa, jfa, jself;
            int ncluof = clulevarry[nlf].nClustOLv,
            nclusf = clulevarry[nlf].nClustSLv;
            Matrix<bool> ifnear(ncluof, nclusf, false);

            for(i=0; i<ncluof; i++)
            {
                for(j=0; j<clulevarry[nlf].ClastArLv[i].NumNeaClus; j++)
                    ifnear(i, clulevarry[nlf].ClastArLv[i].NumsNeaClus[j]) = true;
            }

            for(i=0; i<clulevarry[nlv].nClustOLv; i++)
            {
                ifa = clulevarry[nlv].ClastArLv[i].nuFather;
                nclufar[i] = 0;
                for(j=0; j<clulevarry[nlv].ClastArLv[i].NumFarClus; j++)
                {
                    jself = clulevarry[nlv].ClastArLv[i].NumsFarClus[j];
                    jfa = clulevarry[nlv].ClastArLv[jself].nuFather;

                    if(ifnear(ifa, jfa))
                    {
                        nuclfar(i, nclufar[i]++) = jself;
                    }
                }
            }
        } // end of ELSE

        for(i=0; i<clulevarry[nlv].nClustOLv; i++)
        {
            clulevarry[nlv].ClastArLv[i].NumFanClus = nclufar[i];

            clulevarry[nlv].ClastArLv[i].NumsFanClus = new int[nclufar[i]];
            for(j=0; j<nclufar[i]; j++)
                clulevarry[nlv].ClastArLv[i].NumsFanClus[j] = nuclfar(i, j);
        }

        nrosdmtx = clulevarry[nlv].nClustOLv*clulevarry[nlv].nPoinSpheLv;
        nnonzerdmtx = 0;
        for(i=0; i<clulevarry[nlv].nClustOLv; i++) nnonzerdmtx += nclufar[i];
        nnonzerdmtx *= clulevarry[nlv].nPoinSpheLv;

        dmtxlev[nlv].jcolDmxLv = new int[nnonzerdmtx];
        dmtxlev[nlv].irowDmxLv = new int[nrosdmtx + 1];

        irowd = nnze = 0;
        for(i=0; i<clulevarry[nlv].nClustOLv; i++)
        {
            for(k=0; k<clulevarry[nlv].nPoinSpheLv; k++)
            {
                dmtxlev[nlv].irowDmxLv[irowd++] =  nnze;
                for(j=0; j<nclufar[i]; j++)
                    dmtxlev[nlv].jcolDmxLv[nnze++] =
                    nuclfar(i,j)*clulevarry[nlv].nPoinSpheLv + k;
            }
        }

        dmtxlev[nlv].irowDmxLv[irowd] =  nnze;
        if(maxRowNumberD_ < nrosdmtx) maxRowNumberD_ = nrosdmtx;
        dmtxlev[nlv].nRowsD = nrosdmtx;
        dmtxlev[nlv].nEntriesD = nnonzerdmtx;
    } // end of loop NLV
}

// generate the S-matrices at all levels
void NC_AllocateSmtxMLFMM
(
	ofstream& NCout
)
{
	int nlv, i, j, ij, kcol, itr;
	Vector<int> nuclusrow(numRowsOfCoefficientMatrix_);

	// loop over the levels
	for(nlv=0; nlv<numClusterLevels_; nlv++)
	{
		// initialize
		for(i=0; i<numRowsOfCoefficientMatrix_; i++) nuclusrow[i] = -1;

		// compute the cluster number of each row
		for(i=0; i<clulevarry[nlv].nClustOLv; i++)
		{
			for(j=0; j<clulevarry[nlv].ClastArLv[i].NumOfEl; j++)
			{
				ij = clulevarry[nlv].ClastArLv[i].NumsOfEl[j];
				nuclusrow[jelist[ij][0]] = i;
				if(clulevarry[nlv].ClastArLv[i].NDOFsPeEl == 2)
					nuclusrow[jelist[ij][1]] = i;
			}
		}

		// compute number of nonzeros of the S-matrix
		smtxlev[nlv].nEntriesS = 0;
		for(i=0; i<numRowsOfCoefficientMatrix_; i++)
		{
			if(nuclusrow[i] >= 0) smtxlev[nlv].nEntriesS += clulevarry[nlv].nPoinSpheLv;
		}

		// create the address arrays of the S-matrix
		smtxlev[nlv].jcolSmxLv = new int[smtxlev[nlv].nEntriesS];
		smtxlev[nlv].irowSmxLv = new int[numRowsOfCoefficientMatrix_ + 1];

		// compute the address arrays for the S-matrix
		itr = 0;
		for(i=0; i<numRowsOfCoefficientMatrix_; i++)
		{
			smtxlev[nlv].irowSmxLv[i] = itr;
			if(nuclusrow[i] >= 0)
			{
				kcol = clulevarry[nlv].nPoinSpheLv*nuclusrow[i];
				for(j=0; j<clulevarry[nlv].nPoinSpheLv; j++)
					smtxlev[nlv].jcolSmxLv[itr++] = kcol++;
			}
		}
		smtxlev[nlv].irowSmxLv[numRowsOfCoefficientMatrix_] = itr;
	} // end of loop NLV
}

// create the arrays for storing the near field matrix of SLFMBEM or MLFMBEM
void NC_AllocateNmtx
(
	ofstream& NCout
 )
{
    int i, j, k, neli, nelij, nel0, nori_clu;
    Vector<bool> bneacluo(numOriginalClusters_);

    // compute number of the nonzero entries of the near field matrix
    numNonZeroEntries_ = 0;
    for(i=0; i<numOriginalClusters_; i++)
    {
        neli = ClustArray[i].NumOfDOFs;
        bneacluo = false;
        for(j=0; j<ClustArray[i].NumNeaClus; j++)
        {
            nori_clu = ClustArray[ClustArray[i].NumsNeaClus[j]].OriClust;
            if(!bneacluo[nori_clu]) bneacluo[nori_clu] = true;
        }
        nelij = 0;
        for(j=0; j<numOriginalClusters_; j++) if(bneacluo[j]) nelij += ClustArray[j].NumOfDOFs;
        numNonZeroEntries_ += neli*nelij;
    }


    jcolnea = new int[numNonZeroEntries_];
    irownea = new int[numRowsOfCoefficientMatrix_ + 1];
    int ij, ine, ij0, ij1, ijk, idfs;

    // compute the bool variable indicating if the T-vector should be computed
    boolComputeTVector_ = false;
    for(i=0; i<numOriginalReflectedClusters_; i++)
    {
        nel0 = ClustArray[i].NumsOfEl[0];
	// changed by kreiza, include also the PRES condition
	//        if(ClustArray[i].listElementPropertyEl == 0 && (ibval[nel0] == 0 || ibval[nel0] == 2))
	if( ClustArray[i].listElementPropertyEl == 0 && (ibval[nel0] < 3 ) )
        {
            boolComputeTVector_ = true;
            break;
        }
    }

    // compute for the near field matrix:
    //		(1) array of the column numbers of each nonzero
    //      (2) array of the adresses of the first nonzero in each row
    Matrix<bool> nuclus(numRowsOfCoefficientMatrix_, numOriginalClusters_, false);

    for(i=0; i<numOriginalClusters_; i++)
    {
        ine = ClustArray[i].NumNeaClus;
        idfs = ClustArray[i].NDOFsPeEl;
        for(j=0; j<ClustArray[i].NumOfEl; j++)
        {
            ij = ClustArray[i].NumsOfEl[j];
            ij0 = jelist[ij][0];
            for(k=0; k<ine; k++)
            {
                nori_clu = ClustArray[ClustArray[i].NumsNeaClus[k]].OriClust;
                if(!nuclus(ij0, nori_clu)) nuclus(ij0, nori_clu) =  true;
            }
            if(idfs == 2)
            {
                ij1 = jelist[ij][1];
                for(k=0; k<ine; k++)
                {
                    nori_clu = ClustArray[ClustArray[i].NumsNeaClus[k]].OriClust;
                    if(!nuclus(ij1, nori_clu)) nuclus(ij1, nori_clu) =  true;
                }
            }
        }
    }

    ijk = 0;
    for(i=0; i<numRowsOfCoefficientMatrix_; i++)
    {
        irownea[i] = ijk;
        for(j=0; j<numOriginalClusters_; j++)
        {
            if(!nuclus(i, j)) continue;
            for(k=0; k<ClustArray[j].NumOfEl; k++)
            {
                ij0 = ClustArray[j].NumsOfEl[k];
                jcolnea[ijk++] = jelist[ij0][0];
                if(ClustArray[j].NDOFsPeEl == 2) jcolnea[ijk++] = jelist[ij0][1];
            }
        }
    }
    irownea[i] = ijk;

    if(ijk != numNonZeroEntries_)
        NC_Error_Exit_2(NCout, "ijk must equal to numNonZeroEntries_!", "ijk = ", ijk,
                        "numNonZeroEntries_ = ", numNonZeroEntries_);
}

// check the generated clusters, see if each boundary elements is included by the clusters just once
void NC_CheckClusters
(
	ofstream& NCout,
	const int& n_levl,			//I: level number of the clusters
	const int& n_clustero,		//I: number of original clusters
	ElCluster *clust_arry		//I: array of clusters
)
{
	int i, j, nfalsregi = 0;
	Vector<int> ntimeselregi(numElements_, 0);

	// loop over clusters
	for(i=0; i<n_clustero; i++)
	{
		for(j=0; j<clust_arry[i].NumOfEl; j++)
		{
			ntimeselregi[clust_arry[i].NumsOfEl[j]]++;
		}

	}

	// loop over elements
	for(i=0; i<numElements_; i++)
	{
		if(listElementProperty[i] == 2) continue;
		if(ntimeselregi[i] == 0)
		{
			nfalsregi++;
			NCout << "\nThe element " << extNumbersOfElements[i] << " is not included in any cluster"
				<< endl;
			cout << "\nThe element " << extNumbersOfElements[i] << " is not included in any cluster"
				<< endl;
		}
		else if(ntimeselregi[i] > 1)
		{
			nfalsregi++;
			NCout << "\nThe element " << extNumbersOfElements[i] << " is included in " << ntimeselregi[i]
				<< " clusters" << endl;
			cout << "\nThe element " << extNumbersOfElements[i] << " is included in " << ntimeselregi[i]
				<< " clusters" << endl;
		}
	}

	// if any element is false included by clusters, give a warning or stop the job
	if(nfalsregi > 0)
	{
		NC_Error_Exit_2(NCout,
			"Elements are not or more than one times included by clusters",
			"level number of the clusters = ", n_levl,
			"number of such elements      = ", nfalsregi);
	}
}

// delete arrays generated by addressing computations
void NC_DeleteArrays(const int& imultpol, const int& nlevmlfm, const int& ifdelNFC)
{
	int i, j;

	// delete the D-, T- and S-matrices
	switch(imultpol)
	{
	case 1:
	case 2:
		delete [] dmtxlev;
		break;
	case 3:
		delete [] dmtxlev;
		delete [] tmtxlev;
		delete [] smtxlev;
		break;
	}

	if(imultpol && numInternalPointsClusters_ > 0)
	{
		for(i=0; i<numInternalPointsClusters_; i++) delete [] ipcluarry[i].NumsOfIps;
		delete [] ipcluarry;
	}

	switch(imultpol)
	{
	case 1: // SLFMBEM
		for(i=0; i<numOriginalReflectedClusters_; i++)
		{
			delete [] clustarry[i].NumsFarClus;
			delete [] clustarry[i].NumsNeaClus;
			delete [] clustarry[i].NumsOfEl;
		}
		delete [] clustarry;
		break;
	case 3: // DMLFMBEM
	  for(i=0; i<nlevmlfm; i++)
	    {
	      for(j=0; j<clulevarry[i].nClustSLv; j++)
		{
		  // kreiza 2024: there is the theoretical chance that
		  //      no Far Field Clusters exist, check for that
		  if( clulevarry[i].ClastArLv[j].NumsFarClus != NULL )
		    delete [] clulevarry[i].ClastArLv[j].NumsFarClus;
		  delete [] clulevarry[i].ClastArLv[j].NumsNeaClus;
		  delete [] clulevarry[i].ClastArLv[j].NumsOfEl;
		}
	      
	      if(ifdelNFC)
		{
		  for(j=0; j<clulevarry[i].nClustOLv; j++)
		    {
		      delete [] clulevarry[i].ClastArLv[j].NumsFanClus;
		    }
		}
	      
	      delete [] clulevarry[i].ClastArLv;
	      
	    }
	  delete [] clulevarry;
	  break;
	}

	if(imultpol) // FMBEM
	{
		delete [] jcolnea;
		delete [] irownea;
	}

	delete [] zrhs;
}


//--------------------------------------------------------------------------
// << A11 >> seek very small clusters, attach them to the nearst ones and cancel them
//
// original Code by Chen, modified by kreiza 2018
//
//  Notes and warnings:
//    *  thresfac1 is hardcoded with 3.0, the factor was set by Chen and
//       probably tested to some degree, so in general this factor can be
//       trusted
//--------------------------------------------------------------------------
void NC_CancelSmallClusters
(
 Vector<int>& ne_clus,	// I,O: number of elements in each cluster
 Vector<int>& nue_clus,	// I,O: numbers of elements in each cluster
 int& nclus,				// I  : number of clusters
 int& ibg                // I  : number of element group
 )
{
  int i, j, k, i1, iel, ine, nel_max = 0, n_nod, l1, m;
  double dwk, dmin, threshfac1 = 3.0, dispo, rdmax, rdmin, rdaver;
  bool iffewel = false;

  // maximum number of elements in a cluster
  for(i=0; i<nclus; i++) if(ne_clus[i] > nel_max) nel_max = ne_clus[i];

  // vector to store nodal numbers of each cluster
  Vector<int> nu_nod_cl(nel_max*NNODPE);

  // vector of radius of each cluster
  Vector<double> radi_cl(nclus);

  // create the array of coordinates of the center of each cluster
  double **crdcenters = new double*[nclus];
  //	crdcenters = new double*[nclus];
  if(!crdcenters) {
    fprintf(stderr,"Could not allocate crdcenters\n");
    exit(-1);
  }
  for(i=0; i<nclus; i++) {
    crdcenters[i] = new double[NDIM];
    if(!crdcenters[i]) {
      fprintf(stderr,"Could not allocate crdcenters[i]\n");
      exit(-1);
    }
  }
  for(i=0; i<nclus; i++) for(j=0; j<NDIM; j++) crdcenters[i][j] = 0.0;

  // loop over clusters
  i1 = 0;
  for(i=0; i<nclus; i++)
    {
      if(ne_clus[i] == 0) continue;

      // compute number and numbers of the nodes of the current cluster
      n_nod = 0;
      for(j=0; j<ne_clus[i]; j++)
	{
	  iel = nue_clus[i1++];
	  for(k=0; k<listNumberNodesPerElement[iel]; k++)
	    {
	      l1 = elementsConnectivity[iel][k];
	      for(m=0; m<n_nod; m++) if(nu_nod_cl[m] == l1) goto LbclustAray1;
	      nu_nod_cl[n_nod++] = l1;
	    LbclustAray1: continue;
	    }
	}

      // compute the coordinates of the center of the cluster
      for(j=0; j<n_nod; j++)
	{
	  for(k=0; k<NDIM; k++) crdcenters[i][k] += nodesCoordinates[nu_nod_cl[j]][k];
	}
      for(k=0; k<NDIM; k++) crdcenters[i][k] /= (double)n_nod;

      // compute the radius of the cluster
      radi_cl[i] = 0.0;
      for(j=0; j<n_nod; j++)
	{
	  dispo = Dispoi_dim3_(crdcenters[i], nodesCoordinates[nu_nod_cl[j]]);
	  if(dispo > radi_cl[i]) radi_cl[i] = dispo;
	}
    } // end of loop I

  // compute maximum, average and minimum radius of the clusters
  rdmax = rdaver = i1 = 0;
  rdmin = 1.0e20;
  for(i=0; i<nclus; i++)
    {
      if(ne_clus[i] == 0) continue;
      rdaver += radi_cl[i];
      i1++;
      if(radi_cl[i] < rdmin) rdmin = radi_cl[i];
      if(radi_cl[i] > rdmax) rdmax = radi_cl[i];
    }
  rdaver /= (double)i1;

  // see if there are small clusters
  if(rdaver > threshfac1*rdmin) iffewel = true;

  if(!iffewel) {
    for(i=0; i<nclus; i++) delete [] crdcenters[i];
    delete [] crdcenters;
    return;
  }

  // search the nearest cluster for each small cluster
  Vector<int> NearClus(nclus, -1);
  for(i=0; i<nclus; i++)
    {
      if(ne_clus[i] == 0 || threshfac1*radi_cl[i] >= rdaver) continue;
      dmin = 1.0e20;
      for(j=0; j<nclus; j++)
	{
	  if(i == j || ne_clus[j] == 0 || threshfac1*radi_cl[j] < rdaver) continue;
	  dwk = Dispoi_dim3_(crdcenters[i], crdcenters[j]);
	  if(dwk < dmin)
	    {
	      dmin = dwk;
	      NearClus[i] = j;
	    }
	}
    }

  // store the element numbers in a matrix
  int jdimtx = nel_max;
  for(i=0; i<nclus; i++)
    {
      if(NearClus[i] >= 0)
	{

	  /* kreiza 01.04.2020
	     the old version had the problem, that two clusters can
	     have the same neightbour cluster, and then jdimtx may be
	     too small, thus we decided to waste a bit more memory space and
	     add the dimension of the small cluster to jdimtx
	     old version:
	     j = ne_clus[i] + ne_clus[NearClus[i]];
	     if(j > jdimtx) jdimtx = j;
	  */
	  jdimtx = jdimtx + ne_clus[i];
	}
    }
  Matrix<int> NuElClust(nclus, jdimtx, -1);
  i1 = 0;
  for(i=0; i<nclus; i++)
    {
      if(ne_clus[i] == 0) continue;
      for(j=0; j<ne_clus[i]; j++) NuElClust(i, j) = nue_clus[i1++];
    }

  // attach the small clusters to the nearst clusters
  for(i=0; i<nclus; i++)
    {
      if(NearClus[i] == -1) continue;

      ine = NearClus[i];
      i1 = 0;
      for(j=0; j<jdimtx; j++)
	{
	  if(NuElClust(ine, j) == -1) break;
	  i1++;
	}
      for(j=0; j<ne_clus[i]; j++) NuElClust(ine, i1 + j) = NuElClust(i, j);
      ne_clus[ine] += ne_clus[i];
      ne_clus[i] = 0;
    }

  i1 = 0;
  for(i=0; i<nclus; i++)
    {
      if(ne_clus[i] == 0) continue;
      for(j=0; j<ne_clus[i]; j++) nue_clus[i1++] = NuElClust(i, j);
    }

  // delete the array of coordinates of the center of each cluster
  for(i=0; i<nclus; i++) delete [] crdcenters[i];
  delete [] crdcenters;
}

