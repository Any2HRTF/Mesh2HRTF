/*===========================================================================*\
 *                                                                            *
 *  File name:      NC_Main.cpp                                               *
 *  Description:    main and control programs of NumCalc (Acoustic BEM for    *
 *  the calculation of head-related transfer functions)                       *
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
#include<ctime>                                                                     //
#include<new>                                                                       //
#include<cstddef>                                                                   //
#include<fstream>  // file class functions                                          //
#include<time.h>                                                                    //
#include<stdlib.h>                                                                  //
#ifdef isWindows                                                                    //
#include<direct.h> // for mkdir WINDOWS                                             //
#else                                                                               //
#include<sys/stat.h> // for mkdir UNIX                                              //
#include<unistd.h>                                                                  //
#endif                                                                              //
using namespace std;                                                                //
typedef unsigned int uint;                                                          //
//====================================================================================



void NC_ControlProgram(ofstream&,int,bool,bool);
void NC_FrequencyInformations(ostream&, ofstream&);



extern void NC_Read(ofstream&,FILE *,char*,string[],double*,bool);
extern void NC_ReadBasicParametersA(ofstream&,FILE* inputFile_,char*, string[]);
extern int NC_DeclareArrays(ofstream&,double*,bool);
extern void NC_AllocateSDTmtxsSLFMM(ofstream&);
void NC_AllocateSDTmtxsMLFMM(ofstream&);
extern void NC_DeleteArrays(const int&,const int&,const int&);
extern void NC_SetupEquationSystem(ofstream&);
extern void NC_UpdateFreqCurves(ofstream&,double*);
extern void NC_PostProcessing(ofstream&);

// 17.11.22  Let's include LAPACK
#ifdef USE_LAPACK
extern "C" {
#include<lapacke_config.h>
#include<lapacke_utils.h>
#include<lapacke.h>
}
#endif


// variables declared in NC_ConstantsVariables.h
ofstream NCout;
FILE *inputFile_, *tmpFile_;
//char jobTitle_[SIZE_LINE], versionNumber_[30];
string jobTitle_, versionNumber_;
int currentFrequency_, numFrequencies_, istart_;
double firstFrequencyStep_, frequencyStepIncrement_;
int numElementGroups_, numElements_, numNodes_, minExpansionTermsFMM_, liorcn_, numSymmetricPlanes_, numIncidentPlaneWaves_, numPointSources_,
isInternalProblem_;
double speedOfSound_, densityOfMedium_, harmonicTimeFactor_ = 1.0, waveNumbers_, waveLength_, averageElementArea_, frequency_,
omega1_, rpfact_, Tao_, delta_, epsilon_, ClusEdgL0_;
int numNodesOfBoundaryMesh_, numNodesOfEvaluationMesh_, numElementsOfBoundaryMesh_, numElementsOfEvaluationMesh_, numRowsOfCoefficientMatrix_, numComponentsOfCoefficientMatrix_, numReflectionsOfElements_, numCurvesFrequency_;
int methodBEM_ = 0, numNonZeroEntries_ = 0, numClusterLevels_, nlevtop_,
maxRowNumberD_, numOriginalClusters_, numOriginalReflectedClusters_, numInternalPointsClusters_, methodPreconditioner_, scanningDegreeLU_, methodSolver_;
/*
  set the max number of iterations
*/
int niter_max_ = 250;
double farFieldClusterFactor_ = 2.0/sqrt(3.0), minClusterDistance_;
/*
kreiza 11.2021
fabian had some problems with a verrrrrrry big problem with respect to
stability. I found a distance factor in darve of 2/sqrt(3) that guarantees
absolute convergence of the fmm expansion, let's take this now
*/
//double farFieldClusterFactor_ = sqrt(5.0)/2.0, minClusterDistance_;
double maxClusterRadiusBE_, avgClusterRadiusBE_, minClusterRadiusBE_, maxClusterRadiusRM_, avgClusterRadiusRM_, minClusterRadiusRM_;
int numExpansionTerms_;
int numIntegrationPointsUnitSphere_, numIntegrationPointsThetaDirection_, numIntegrationPointsPhiDirection_;
int methodFMM_;
bool boolComputeTVector_;
int iMirrSour_ = 0, iSpringMass_ = 0;
double diameterBE_;
int *extNumbersOfNodes, *isNodeMeshOrEvalNode, *extNumbersOfElements;
double **nodesCoordinates;
int **elementsConnectivity;
int *listNumberNodesPerElement, *listElementProperty, *listElementsElementGroup;
int *listSymmetryPlanes;
double *listSymmetryPlanesCoordinates;
int *ibval;
Complex **zbval0, **zbvao0, *zbval1, *zbval2, *zbvao1, *zbvao2;
double **elenor, **centel, *areael;
int **jelist;
int *indexOfElementGroup, *numberOfElementsInGroup, *propertyOfGroup;
int *numinw, *numposo;
double **dirinw, **coorps;
int **inwacurv, **posocurv;
Complex *zinwav, *zposoc, *zinwao, *zposoo;
int **ibvcurv, *nucurv, *npcurv;
double **frqcurv, **faccurv;
ElCluster *clustarry;
IpCluster *ipcluarry;
int *jcolnea, *irownea;
int *jcoltmtx, *irowtmtx;
Complex *ztmtx;
Complex *ztvct;
D_mtx_lev *dmtxlev;
T_mtx_lev *tmtxlev;
S_mtx_lev *smtxlev;
int *jcolsmtx, *irowsmtx;
Complex *zsmtx;
double **uvcsphe;
double *weisphe;
Complex *zcoefl;
Complex *zrhs;
ClusterLev *clulevarry;
ElCluster *ClustArray;



// the main program for NumCalc
int main(int argc, char **argv)
{

  int i;          //* used for command line parameters
  int iend = 0;   //* last freq step
  char filename[200];  //* filename for the outputfile
  bool estimate_ram = false; //* parameter for ram estimation
  bool check_normals = false; //* checkflag for normalvector direction
  istart_ = 0;    //* first freq step


  /* check the caommand line */
  i = 1;
  while (i < argc) {
    if(!strcmp(argv[i],"-h")) {
      printf("-istart   int : start index of iteration\n");
      printf("-iend     int : end index of iteration\n");
      printf("-nitermax int : max number of CGS iterations\n");
      printf("-estimate_ram : estimation the RAM consumption of ML-FMM-BEM and write estimate to Memory.txt. Estimate is obtained from the number of non-zeros in the FMM matrices.\n");
      printf("-check_normals : check if all normals point to the same domain\n");
      printf("-h            : this message\n");
      exit(0);
    }
    else if(!strcmp(argv[i],"-istart")) {
      i++;
      if(i<argc) {
	istart_=atoi(argv[i]);
      }
    }
    else if(!strcmp(argv[i],"-iend")) {
      i++;
      if(i<argc)
	iend=atoi(argv[i]);
    }
    else if(!strcmp(argv[i],"-nitermax")) {
      i++;
      if( i < argc )
	niter_max_ = atoi( argv[i] );
    }
    else if(!strcmp(argv[i],"-estimate_ram")) {
      estimate_ram = true;
    }
    else if(!strcmp(argv[i],"-check_normals")) {
      check_normals = true;
    }
    else {
      cerr << "\nNumCalc was called with an unknown parameter or flag. Use NumCalc -h for help.\n";
      exit(-1);
    }
    i++;
  }


  // check if istart is okay
  if( istart_ < 0 ) {
    cerr << "I need to start with a step bigger 0\n";
    exit(-1);
  }

  fprintf(stdout, "istart= %d iend = %d\n", istart_, iend);

	// compute the start time
  time_t lot = time(NULL);
  tm *d_t = localtime(&lot);

	// Start informations
	cout << "\n---------- NumCalc started: " << d_t->tm_mday << "/" << d_t->tm_mon + 1 <<
		"/" << d_t->tm_year + 1900 << " " << d_t->tm_hour << ":" << d_t->tm_min <<
		":" << d_t->tm_sec << " ----------\n" << endl;

#ifdef isWindows

#else
    char Name[150];
    memset(Name, 0, 150);
    gethostname(Name, 150);
    cout << "Running on:         " << Name << "\n" << endl;
#endif


    // open the general output file
    if(istart_ > 0) {
      if( iend >= istart_ )
	sprintf(filename,"NC%d-%d.out",istart_,iend);
      else if (iend > 0) {
	cerr << "istart is bigger than iend. Stopping now!\n";
	exit(-1);
      }
      else
	sprintf(filename,"NCfrom%d.out",istart_);
    }
    else {
      if(iend > 0)
	sprintf(filename,"NCuntil%d.out",iend);
      else
	sprintf(filename,"NC.out");
    }

	ofstream NCout(filename);
	if(!NCout) NC_Error_Exit_0(NCout, "can not open the output stream NCout!");

	// start time
    NCout << "\nStart time: " << d_t->tm_mday << "/" << d_t->tm_mon + 1 << "/" <<
        d_t->tm_year + 1900 << " " << d_t->tm_hour << ":" << d_t->tm_min << ":" <<
        d_t->tm_sec << endl;


      // delete and create the output directories if all freqsteps
      // need to be calculated, otherwise keep the old directory
      // and the data in it
#ifdef isWindows
      int ifmkd;
      if( (istart_ == 0 && iend == 0) && (!estimate_ram)) { // no steps are given, remove the old directory
	if(system("rmdir /s /q be.out")==-1) cout << "\nCannot create directory be.out";
	//	if(system("rmdir /s /q fe.out")==-1) cout << "\nCannot create directory fe.out";
      }
      ifmkd = _mkdir("be.out"); // WINDOWS
      //   ifmkd = _mkdir("fe.out");
#else
      if(istart_ == 0 && iend == 0 && !estimate_ram) {
	if(system("rm -f -r be.out")==-1) cout << "\nCannot create directory be.out";
	//	if(system("rm -f -r fe.out")==-1) cout << "\nCannot create directory fe.out";
      }
      mkdir("be.out", S_IRUSR | S_IWUSR | S_IXUSR | S_IRGRP | S_IWGRP | S_IXGRP | S_IROTH | S_IXOTH);  // UNIX
      //  mkdir("fe.out", S_IRUSR | S_IWUSR | S_IXUSR | S_IRGRP | S_IWGRP | S_IXGRP | S_IROTH | S_IXOTH);
#endif

	// call the control program
      NC_ControlProgram(NCout,iend,estimate_ram,check_normals);

	// compute the end time
	lot = time(NULL);
	d_t = localtime(&lot);

	// end informations
    NCout << "\nEnd time: " << d_t->tm_mday << "/" << d_t->tm_mon + 1 << "/" <<
        d_t->tm_year + 1900 << " " << d_t->tm_hour << ":" << d_t->tm_min << ":" <<
        d_t->tm_sec << endl;

	cout << "\n---------- NumCalc ended: " << d_t->tm_mday << "/" << d_t->tm_mon + 1 <<
		"/" << d_t->tm_year + 1900 << " " << d_t->tm_hour << ":" << d_t->tm_min <<
		":" << d_t->tm_sec << " ----------" << endl;

	return(0);
}

// control program
void NC_ControlProgram(ofstream& NCout,int iend, bool estimate_ram, bool check_normals)
{
  double *Freqs = nullptr;

  int i, j, l, nptsouc, ifdiff, terml, npthet, npphi, npsph;
  double rw, rd_max;
  double H_Argum_Min = 0.0; // minimum allowable argument of the spherical Hankel function

  // terms of the input line
  string chterms[NTRM];
  // input line
  char* chinpline = new char[SIZE];

  // real time, U. Breymann p 629
  time_t ltim[7], lti_eqa = 0, lti_sol = 0, lti_pos = 0, lti_est, lti_sst, lti_pst;


  ofstream NCestim;

  if(estimate_ram) {
    NCestim.open("Memory.txt");
    if(!NCestim)
      if(!NCout) NC_Error_Exit_0(NCout, "can not open the output stream NCestim!");
  }

  time(&ltim[0]);





  // open the input file
  inputFile_ = fopen("NC.inp", "r");
  if(inputFile_) // If input file exists
    {
      // write the basic parameters
      NCout << "\n\n>> G E N E R A L   P A R A M E T E R S <<\n" << endl;

      cout << "3D analysis" << endl;
      NCout << "\n3D analysis" << endl;

      // read the basic parameters (part 1)
      NC_ReadBasicParametersA(NCout, inputFile_, chinpline, chterms);

      // create the frequency vector
      Freqs = new double[numFrequencies_];

      // read the input file
      NC_Read(NCout, inputFile_, chinpline, chterms, Freqs, check_normals);
    }
  else // If input file does not exist
    {
      NC_Error_Exit_0(NCout, "No input file found!");
    }

  // delete the input line
  delete [] chinpline;

  // absolute value of number of the point sources (when numPointSources_ = -1, a "test point source" is generated and used for testing the accuracy of the program)
  nptsouc = abs(numPointSources_);

  time(&ltim[1]);

  // Value of the smallst allowable argument of the spherical Hankel function (the bigst allowable norm of the spherical Hankel function = 3.0e13):
  if(minExpansionTermsFMM_ == 5) {H_Argum_Min = 0.0051;}
  else if(minExpansionTermsFMM_ == 6) {H_Argum_Min = 0.0177;}
  else if(minExpansionTermsFMM_ == 7) {H_Argum_Min = 0.0445;}
  else if(minExpansionTermsFMM_ == 8) {H_Argum_Min = 0.0905;}
  else if(minExpansionTermsFMM_ == 9) {H_Argum_Min = 0.1597;}
  else if(minExpansionTermsFMM_ == 10) {H_Argum_Min = 0.2547;}

  time(&ltim[2]);

  // loop over frequencies

  if(iend > 0)  {
    if(iend > numFrequencies_)
      iend = numFrequencies_;
  }
  if(iend == 0)
    iend = numFrequencies_;
  if(istart_ > 0)
    istart_--;



  for(currentFrequency_ = istart_; currentFrequency_ < iend; currentFrequency_++) {
    time(&ltim[3]);

    // compute some parameters
    frequency_ = Freqs[currentFrequency_];					// frequency
    omega1_ = 2.0*PI*frequency_;				// angular frequency
    rpfact_ = densityOfMedium_*omega1_*harmonicTimeFactor_;		// = sound pressure / velocity potential
    waveNumbers_ = omega1_/speedOfSound_;				// wave number
    waveLength_ = speedOfSound_/frequency_;				// wave length
    minClusterDistance_ = H_Argum_Min/waveNumbers_;        // minimum distance between interacting clusters

    // address computations

    // write the step informations
    cout << "\nStep " << currentFrequency_ + 1 << ", Frequency = " << frequency_ << " Hz" << endl;

    NCout << "\n\n\n>> S T E P   N U M B E R   A N D   F R E Q U E N C Y <<\n"
	  << endl;
    NCout << "Step " << currentFrequency_ + 1 << ", Frequency = " << frequency_ << " Hz"
	  << endl;

    // address computation
    ifdiff = NC_DeclareArrays(NCout, Freqs, estimate_ram);

    // write the analysis type and cluster informations
    if(ifdiff) NC_FrequencyInformations(cout, NCout);

    // Asign some arrays on the finest level to the corresponding actual arrays
    if(methodFMM_ >= 2 && currentFrequency_ > istart_)
      {
	ClustArray = clulevarry[nlevtop_].ClastArLv;

	numOriginalClusters_ = clulevarry[nlevtop_].nClustOLv;
	numOriginalReflectedClusters_ = clulevarry[nlevtop_].nClustSLv;

	numExpansionTerms_ = clulevarry[nlevtop_].nExpaTermLv;
	numIntegrationPointsUnitSphere_ = clulevarry[nlevtop_].nPoinSpheLv;
	maxClusterRadiusBE_ = clulevarry[nlevtop_].RadiMaxLv;
	avgClusterRadiusBE_ = clulevarry[nlevtop_].RadiAveLv;
	minClusterRadiusBE_ = clulevarry[nlevtop_].RadiMinLv;

	weisphe = clulevarry[nlevtop_].weisphe;
	uvcsphe = clulevarry[nlevtop_].uvcsphe;
      }

    // read the velocity boundary conditions from the output files of the structural analysis

    // If curves are used for prescribing the dependancies of the admittance boundary conditions and the incident waves on the frequencies, adate their walues according to the corresponding curves
    NC_UpdateFreqCurves(NCout, Freqs);

    // Compute parameters for FMBEM (number of the expansion terms, Gauss points at the unit sphere);
    switch(methodFMM_) // = 0: TBEM; 1: SLFMBEM; 3: DMLFMBEM
      {
      case 1: // SLFMBEM
	// order of expension
	rw = 2.0*maxClusterRadiusBE_*waveNumbers_ + 1.8*log10(2.0*maxClusterRadiusBE_*waveNumbers_ + PI);
	numExpansionTerms_ = (int)(rw);
	if(rw - (double)numExpansionTerms_ >= 0.5) numExpansionTerms_++;
	if(numExpansionTerms_ < minExpansionTermsFMM_) numExpansionTerms_ = minExpansionTermsFMM_;

	// number of the Gauss points in the Theta direction
	numIntegrationPointsThetaDirection_ = numExpansionTerms_;

	// number of the Gauss points must not excess the allowable maximum
	if(numIntegrationPointsThetaDirection_ > N_GAUORDER)
	  {
	    NC_Error_Exit_2(NCout,
			    "Too many integral points in the theta-direction!",
			    "Number of integration points = ", numIntegrationPointsThetaDirection_,
			    "Maximum allowable number of integration points = ", N_GAUORDER);
	  }

	// number of points in the Phi direction
	numIntegrationPointsPhiDirection_ = 2*numIntegrationPointsThetaDirection_;

	// number of points on the unit sphere
	numIntegrationPointsUnitSphere_ = numIntegrationPointsThetaDirection_*numIntegrationPointsPhiDirection_;
	break;
      case 3: // DMLFMBEM
	cout << endl;
	NCout << endl;
	for(l=0; l<numClusterLevels_; l++)
	  {
	    // order of expension
	    rd_max = clulevarry[l].RadiMaxLv;
	    rw = 2.0*rd_max*waveNumbers_ + 1.8*log10(2.0*rd_max*waveNumbers_ + PI);
	    terml = (int)(rw);
	    if(rw - (double)terml >= 0.5) terml++;
	    if(terml < minExpansionTermsFMM_) terml = minExpansionTermsFMM_;
	    npthet = clulevarry[l].nExpaTermLv = terml;

	    // For the second interpolation scheme of IMLFMBEM, numbers of Gaussean points in the theta-direction on two successive levels must be even at least on one of them
	    clulevarry[l].nPoinThetLv = npthet;

	    // number of the Gauss points must not excess the allowable maximum
	    if(npthet > N_GAUORDER)
	      {
		NC_Error_Exit_2(NCout,
				"Too many integral points in the theta-direction!",
				"Number of integration points = ", npthet,
				"Maximum allowable number of integration points = ", N_GAUORDER);
	      }

	    // number of points in the Phi direction
	    npphi = 2*npthet;
	    clulevarry[l].nPoinPhiLv = npphi;

	    // number of points on the unit sphere
	    clulevarry[l].nPoinSpheLv = npthet*npphi;
	  } // end of loop L

	// Asign some arrays on the finest level to the corresponding actual arrays
	numExpansionTerms_ = clulevarry[nlevtop_].nExpaTermLv;
	numIntegrationPointsThetaDirection_ = clulevarry[nlevtop_].nPoinThetLv;
	numIntegrationPointsPhiDirection_ = clulevarry[nlevtop_].nPoinPhiLv;
	numIntegrationPointsUnitSphere_ = clulevarry[nlevtop_].nPoinSpheLv;
	break;
      default:
	break;
      } // end of SWITCH

    // write infomations about FFM expansion and integration points on the unit sphere
    if(methodFMM_) {
      cout << "\nInformations about the FMM expansion:" << endl;
      cout <<
	"   Level    N. expa. terms    N. p. theta    N. p. phi    N. p. sphere"
	   << endl;
      NCout << "\nInformations about the FMM expansion:" << endl;
      NCout <<
	"   Level    N. expa. terms    N. p. theta    N. p. phi    N. p. sphere"
	    << endl;
      for(l=0; l<numClusterLevels_; l++) {
	if(methodFMM_ == 1) {
	  terml = numExpansionTerms_;
	  npthet = numIntegrationPointsThetaDirection_;
	  npphi = numIntegrationPointsPhiDirection_;
	  npsph = numIntegrationPointsUnitSphere_;
	} else {
	  terml = clulevarry[l].nExpaTermLv;
	  npthet = clulevarry[l].nPoinThetLv;
	  npphi = clulevarry[l].nPoinPhiLv;
	  npsph = clulevarry[l].nPoinSpheLv;
	}
	cout << "   " << setw(3) << l << "       " << setw(5) << terml <<
	  "           " << setw(5) << npthet <<
	  "           " << setw(5) << npphi <<
	  "         " << setw(5) << npsph << endl;
	NCout << "   " << setw(3) << l << "       " << setw(5) << terml <<
	  "           " << setw(5) << npthet <<
	  "           " << setw(5) << npphi <<
	  "         " << setw(5) << npsph << endl;
      }
    }

    if( estimate_ram ) {
      // do not reserve too much space. Just try to figure out
      // how much RAM is needed
      uint nint = 0;
      uint ncmplx = 0;
      int i,j,nnonzers;
      Vector<bool> bneacluo(numOriginalClusters_);
      Vector<int> nclufar(numOriginalClusters_);
      int neli, nelij, nori_clu, nlv;


      ncmplx = numRowsOfCoefficientMatrix_; // for zrhs
      switch( methodFMM_ ) {
      case 0: // TBEM
	cout << "TBEM: \n";
	cout << "Number of entries right hand side: " << ncmplx << "\n";
	cout << "Number of Matrix entries: " << ncmplx * ncmplx << "\n";
	ncmplx += numRowsOfCoefficientMatrix_ * numRowsOfCoefficientMatrix_; // for the matrix
	break;
      case 1: //SLFMM
	/* **********************************************
	   The unit sphere
	   the quadrature nodes on the unit sphere,
	   3 coordinates + 1 weight = 4, but dble instead of complex /2 = 2
	   ********************************************** */
	ncmplx += numIntegrationPointsUnitSphere_ * 2;
	/* ************************************************
	  the Nearfield matrix
	  ************************************************* */
	nnonzers = 0;
	for(i=0; i<numOriginalClusters_; i++) {
	  neli = ClustArray[i].NumOfDOFs;
	  bneacluo = false;
	  for(j=0; j<ClustArray[i].NumNeaClus; j++)
	    {
	      nori_clu = ClustArray[ClustArray[i].NumsNeaClus[j]].OriClust;
	      if(!bneacluo[nori_clu]) bneacluo[nori_clu] = true;
	    }
	  nelij = 0;
	  for(j=0; j<numOriginalClusters_; j++)
	    if(bneacluo[j])
	      nelij += ClustArray[j].NumOfDOFs;
	  nnonzers += neli*nelij;
	}

	cout << "SLFMM:\n";
	cout << "Nearfield: " << nnonzers << "\n";
	nint += (uint)nnonzers + (uint)numRowsOfCoefficientMatrix_ + 1;
	ncmplx += (uint)nnonzers;

	/* ******************************************************
	  the T matrix,
	  since we don't want to think about b.c. too much
	  we'll just do the calculations, in general we should end up with
	  constant * number_of_elements * number_nodes_sphere
	  ****************************************************** */

	nnonzers = 0;
	for(i=0; i<numOriginalReflectedClusters_; i++)
	  nnonzers += clustarry[i].NumOfDOFs;
	nnonzers *= numIntegrationPointsUnitSphere_;

	ncmplx += (uint)nnonzers;
	nint += (uint)nnonzers + (uint)numIntegrationPointsUnitSphere_*numOriginalReflectedClusters_ + 1;

	cout << "T Matrix: " << nnonzers << "\n";

	/* *********************************************
	   The D Matrix
	   ********************************************* */
	nnonzers = 0;
	for(i=0; i<numOriginalClusters_; i++)
	  {
	    nnonzers += clustarry[i].NumFarClus;
	  }
	nnonzers *= numIntegrationPointsUnitSphere_;
	nint += (uint)nnonzers + (uint)numOriginalClusters_*numIntegrationPointsUnitSphere_;
	ncmplx += (uint)nnonzers;

	cout << "D Matrix: " << nnonzers << "\n";
	/* **************************************************
	   The S Matrix
	   some estimate shortcut: all rows are assumed to have entries
	   ************************************************** */
	nnonzers = numRowsOfCoefficientMatrix_ * numIntegrationPointsUnitSphere_;
	nint += (uint)nnonzers + (uint)numRowsOfCoefficientMatrix_ + 1;
	ncmplx += (uint)nnonzers;

	cout << "S Matrix: " << nnonzers << "\n";
	break;
      case 3:  //MLFMM no interpolation
	/* **********************************************
	   The unit sphere
	   the quadrature nodes on the unit sphere,
	   3 coordinates + 1 weight = 4, but dble instead of complex /2 = 2
	   ********************************************** */
	ncmplx += numIntegrationPointsUnitSphere_ * 2;
	/* ************************************************
	   the Nearfield matrix
	   ************************************************* */
	nnonzers = 0;
	for(i=0; i<numOriginalClusters_; i++) {
	  neli = ClustArray[i].NumOfDOFs;
	  bneacluo = false;
	  for(j=0; j<ClustArray[i].NumNeaClus; j++)
	    {
	      nori_clu = ClustArray[ClustArray[i].NumsNeaClus[j]].OriClust;
	      if(!bneacluo[nori_clu]) bneacluo[nori_clu] = true;
	    }
	  nelij = 0;
	  for(j=0; j<numOriginalClusters_; j++)
	    if(bneacluo[j])
	      nelij += ClustArray[j].NumOfDOFs;
	  nnonzers += neli*nelij;
	}

	nint += (uint)nnonzers + (uint)numRowsOfCoefficientMatrix_ + 1;
	ncmplx += (uint)nnonzers;

	cout << "MLFMM: \n";
	cout << "Nearfield " << nnonzers << "\n";

	/* *************************************************
           T Matrix
	   ************************************************* */
	for(nlv=0; nlv<numClusterLevels_; nlv++) {
	  // compute number of rows of the T-matrix and the T-vector
	  nnonzers = 0;
	  for(i=0; i<clulevarry[nlv].nClustSLv; i++)
	    nnonzers += clulevarry[nlv].ClastArLv[i].NumOfDOFs;
	  nnonzers *= clulevarry[nlv].nPoinSpheLv;

	  cout << "Tmatrix level " << nlv << ": " << nnonzers << "\n";

	  nint += (uint)nnonzers + (uint)clulevarry[nlv].nClustSLv*clulevarry[nlv].nPoinSpheLv + 1;
	  ncmplx += (uint)nnonzers;
	}

	/* **************************************************
	   The D Matrices
	   ************************************************** */
	for(nlv=0; nlv<numClusterLevels_; nlv++) {
	  // compute the relevant far clusters for the current level
	  if(nlv == 0) {// the TRUNK level
	    for(i=0; i<clulevarry[nlv].nClustOLv; i++) {
	      nclufar[i] = clulevarry[nlv].ClastArLv[i].NumFarClus;
            }
	  }
	  else {// the BRANCH and LEAF levels
            int nlf = nlv-1, ifa, jfa, jself;
            int ncluof = clulevarry[nlf].nClustOLv,
	      nclusf = clulevarry[nlf].nClustSLv;
            Matrix<bool> ifnear(ncluof, nclusf, false);
	    for(i=0; i<ncluof; i++) {
	      for(j=0; j<clulevarry[nlf].ClastArLv[i].NumNeaClus; j++)
		ifnear(i, clulevarry[nlf].ClastArLv[i].NumsNeaClus[j]) = true;
	    }
	    for(i=0; i<clulevarry[nlv].nClustOLv; i++) {
	      ifa = clulevarry[nlv].ClastArLv[i].nuFather;
	      nclufar[i] = 0;
	      for(j=0; j<clulevarry[nlv].ClastArLv[i].NumFarClus; j++)
                {
		  jself = clulevarry[nlv].ClastArLv[i].NumsFarClus[j];
		  jfa = clulevarry[nlv].ClastArLv[jself].nuFather;
		  if(ifnear(ifa, jfa))
		    nclufar[i]++;
                }
            }
	  } // end of ELSE
	  nnonzers = 0;
	  for(i=0; i<clulevarry[nlv].nClustOLv; i++)
	    nnonzers += nclufar[i];
	  nnonzers *= clulevarry[nlv].nPoinSpheLv;

	  nint += (uint)nnonzers + (uint)clulevarry[nlv].nClustOLv*clulevarry[nlv].nPoinSpheLv + 1;
	  ncmplx += (uint)nnonzers;

	  cout << "D Matrix level " << nlv << ": " << nnonzers << "\n";
	} // loop over levels
	/* ************************************************
	   The S Matrix
	   ************************************************ */
	for(nlv=0; nlv<numClusterLevels_; nlv++) {
	  nnonzers = numRowsOfCoefficientMatrix_ * clulevarry[nlv].nPoinSpheLv;

	  nint += (uint)nnonzers + (uint)numRowsOfCoefficientMatrix_ + 1;
	  ncmplx += (uint)nnonzers;

	  cout << "S Matrix level " << nlv << ": " << nnonzers << "\n";
	}
	break;
      } // the FMM switch
      cout << "Total Nr. of Integers: " << nint << "\n";
      cout << "Total Nr. of Cmplxdble: " << ncmplx << "\n";
      /*      if( nint < 0 || ncmplx < 0 ) {
	// Oh boy, you really run into troubles now, you broke the integer
	// boundary
	cout << "Warning: This estimation is really just a rough lower limit, because of too many entries.\n";
	cout << "RAM Estimation: " << 2.2 * (double)sizeof(int) + 2.2 * (double)sizeof(Complex) << " GByte\n";

	NCestim << currentFrequency_ << " " << frequency_ << " ";
	NCestim << (double)(2.2 * sizeof(int) + 2.2 * sizeof(Complex)) << "\n";
	}*/

	cout << "RAM Estimation: " << (double)(nint * sizeof(int) + ncmplx * sizeof(Complex))/1000.0/1000.0/1000.0 << " GByte\n";
	NCestim << currentFrequency_ + 1 << " " << frequency_ << " ";
	NCestim << (double)(nint * sizeof(int) + ncmplx * sizeof(Complex))/1000.0/1000.0/1000.0 << "\n";

      continue; // end the freq loop
    }



    /* *******************************************************
     *
     * now we again are in the part without the ram estimation
     *
     ******************************************************* */

	// initialize the right hand side vector
    for(i=0; i<numRowsOfCoefficientMatrix_; i++)
      {
	zrhs[i].set(0.0, 0.0);
      }

        // generate and initialize the coefficient matrices
    switch(methodFMM_)
      {
      case 0: // TBEM
	// create and initialize the coefficient matrices
	zcoefl = new Complex[numComponentsOfCoefficientMatrix_];
	for(i=0; i<numComponentsOfCoefficientMatrix_; i++) zcoefl[i].set(0.0, 0.0);
	break;
      case 1: // SLFMBEM
	// compute the auxiliary arrays used for storing the sparse far field matrices
	NC_AllocateSDTmtxsSLFMM(NCout);

	// coordinates and weights of the integration points on the unit sphere
	uvcsphe = new double*[numIntegrationPointsUnitSphere_];
	for(i=0; i<numIntegrationPointsUnitSphere_; i++)
	  {
	    uvcsphe[i] = new double[NDIM];
	  }
	weisphe = new double[numIntegrationPointsUnitSphere_];

	// generate the "T-vector" ( = [T] * {x})
	if(boolComputeTVector_) ztvct = new Complex[numIntegrationPointsUnitSphere_*numOriginalReflectedClusters_];

	// generate and initialize the near field matrix
	zcoefl = new Complex[irownea[numRowsOfCoefficientMatrix_]];
	for(i=0; i<irownea[numRowsOfCoefficientMatrix_]; i++) zcoefl[i].set(0.0, 0.0);

	// generate the "T-matrix" ([T])
	ztmtx = new Complex[irowtmtx[numOriginalReflectedClusters_*numIntegrationPointsUnitSphere_]];

	// generate the "D-matrix" ([D])
	dmtxlev[0].zDmxLv = new Complex[dmtxlev[0].nEntriesD];

	// generate he "S-matrix" ([S])
	zsmtx = new Complex[irowsmtx[numRowsOfCoefficientMatrix_]];

	break;
      case 3: // DMLFMBEM
	// compute the auxiliary arrays used for storing the sparse far field matrices
	NC_AllocateSDTmtxsMLFMM(NCout);

	// create and initialize the near field matrix
	zcoefl = new Complex[irownea[numRowsOfCoefficientMatrix_]];
	for(i=0; i<irownea[numRowsOfCoefficientMatrix_]; i++) zcoefl[i].set(0.0, 0.0);

	// loop over levels
	for(j=0; j<numClusterLevels_; j++)
	  {
	    // number of the integration points on the unit sphere
	    int nthej = clulevarry[j].nPoinThetLv;
	    int npsh = clulevarry[j].nPoinSpheLv;

	    // create the working arrays for the T- and S-matrices
	    clulevarry[j].zwkT = new Complex[clulevarry[j].nClustSLv*npsh];
	    clulevarry[j].zwkS = new Complex[clulevarry[j].nClustOLv*npsh];

	    // create the arrays of coordinates and weights of the integral points on the unit sphere
	    clulevarry[j].uvcsphe = new double*[npsh];
	    for(i=0; i<npsh; i++) clulevarry[j].uvcsphe[i] = new double[NDIM];
	    clulevarry[j].weisphe = new double[npsh];

	    // create the arrays of the coorninates and weights of the Gauss points in the theta-direction
	    clulevarry[j].CrdGauLv = new double[nthej];
	    clulevarry[j].WeiGauLv = new double[nthej];

	    // create the T-matrix and the T-vector, the D- and S-matrices
	    tmtxlev[j].zTmxLv = new Complex[tmtxlev[j].nEntriesT];
	    if(boolComputeTVector_) tmtxlev[j].zTvcLv = new Complex[clulevarry[j].nClustSLv*npsh];
	    dmtxlev[j].zDmxLv = new Complex[dmtxlev[j].nEntriesD];
	    smtxlev[j].zSmxLv = new Complex[smtxlev[j].nEntriesS];
	  }

	break;
      } // end of SWITCH

    // set up the equation system
    NC_SetupEquationSystem(NCout);

    time(&ltim[4]);
    lti_est = ltim[4] - ltim[3];
    lti_eqa += lti_est;

    // solve the equation system
    switch(methodSolver_)
      {
      case 0: // CGS method
	NC_IterativeSolverCGS(NCout);
	break;
      case 4: // direct method, usable only to the TBEM factorize the coefficient matrix
#ifdef USE_LAPACK
	int info = 0;
	int* ipiv;
	cout << "Using LAPACK\n";
	ipiv = new int[numRowsOfCoefficientMatrix_];

	info = LAPACKE_zgetrf(LAPACK_ROW_MAJOR,numRowsOfCoefficientMatrix_, numRowsOfCoefficientMatrix_, (lapack_complex_double*)zcoefl, numRowsOfCoefficientMatrix_, ipiv);
	if(info != 0) {
	  cerr << "Problem with factorization of the stiffness matrix.\n";
	  cerr << "Info = " << info << "\n";
	  exit(-1);
	}
	info = LAPACKE_zgetrs(LAPACK_ROW_MAJOR, 'N', numRowsOfCoefficientMatrix_, 1, (lapack_complex_double*)zcoefl, numRowsOfCoefficientMatrix_, ipiv,  (lapack_complex_double*)zrhs, 1);
	if(info != 0) {
	  cerr << "Problem with the solution of the system (LAPACKE).\n";
	  cerr << "Info = "<< info << "\n";
	  exit(-1);
	}
	delete[] ipiv;
#else
	
	Tfactor_usy(zcoefl, numRowsOfCoefficientMatrix_);

	// forward substitution and backward eliminations
	Tfbelim(zcoefl, zrhs, numRowsOfCoefficientMatrix_);
#endif
	break;
      }

    // destroy the coefficient matrices
    delete [] zcoefl;

    switch(methodFMM_)
      {
      case 1:
	delete [] jcoltmtx;
	delete [] irowtmtx;

	delete [] jcolsmtx;
	delete [] irowsmtx;

	delete [] zsmtx;
	delete [] ztmtx;
	if(boolComputeTVector_) delete [] ztvct;
	break;
      case 3:
	if(boolComputeTVector_)
	  {
	    for(j=0; j<numClusterLevels_; j++) delete [] tmtxlev[j].zTvcLv;
	  }
	for(j=0; j<numClusterLevels_; j++)
	  {
	    delete [] tmtxlev[j].jcolTmxLv;
	    delete [] tmtxlev[j].irowTmxLv;
	    delete [] tmtxlev[j].zTmxLv;

	    delete [] smtxlev[j].jcolSmxLv;
	    delete [] smtxlev[j].irowSmxLv;
	    delete [] smtxlev[j].zSmxLv;
	  }
	break;
      }

    time(&ltim[5]);
    lti_sst = ltim[5] - ltim[4];
    lti_sol += lti_sst;
    
    // post process: compute and output the results
    NC_PostProcessing(NCout);

    switch(methodFMM_)
      {
      case 1: // SLFMBEM
	for(i=0; i<numIntegrationPointsUnitSphere_; i++) delete [] uvcsphe[i];
	delete [] uvcsphe;
	delete [] weisphe;

	delete [] dmtxlev[0].jcolDmxLv;
	delete [] dmtxlev[0].irowDmxLv;
	delete [] dmtxlev[0].zDmxLv;
	break;
      case 3: // DMLFMBEM
	for(j=0; j<numClusterLevels_; j++)
	  {
	    delete [] clulevarry[j].CrdGauLv;
	    delete [] clulevarry[j].WeiGauLv;

	    for(i=0; i<clulevarry[j].nPoinSpheLv; i++) delete [] clulevarry[j].uvcsphe[i];
	    delete [] clulevarry[j].uvcsphe;
	    delete [] clulevarry[j].weisphe;

	    delete [] clulevarry[j].zwkT;
	    delete [] clulevarry[j].zwkS;


	    delete [] dmtxlev[j].jcolDmxLv;
	    delete [] dmtxlev[j].irowDmxLv;
	    delete [] dmtxlev[j].zDmxLv;
	  }
	break;
      }

    time(&ltim[6]);
    lti_pst = ltim[6] - ltim[5];
    lti_pos += lti_pst;

    // time statistic for the current frequency
    NCout << "\nTime Statistic For The Current Step (In Second)\n" << endl;
    NCout << "Assembling the equation system         : " << lti_est << endl;
    NCout << "Solving the equation system            : " << lti_sst << endl;
    NCout << "Post processing                        : " << lti_pst << endl;
    NCout << "Total                                  : " <<
      lti_est + lti_sst + lti_pst << endl;

    cout << "Assembling the equation system         : " << lti_est << endl;
    cout << "Solving the equation system            : " << lti_sst << endl;
    cout << "Post processing                        : " << lti_pst << endl;
    cout << "Total                                  : " <<
      lti_est + lti_sst + lti_pst << endl;
    cout << endl;

    // delete arrays generated by 3-d address computations
    if(currentFrequency_ == iend - 1)  NC_DeleteArrays(methodFMM_, numClusterLevels_, 1);

  } // end of loop I_FREQ (Loop over Frequencies)

  if( !estimate_ram ) {
    NCout << "\n\n\n>> T I M E   S T A T I S T I C   (In   Second) <<\n" << endl;
    NCout << "Input                                  : " << ltim[1] - ltim[0] << endl;
    NCout << "Address computation                    : " << ltim[2] - ltim[1] << endl;
    NCout << "Assembling the equation system         : " << lti_eqa << endl;
    NCout << "Solving the equation system            : " << lti_sol << endl;
    NCout << "Post processing                        : " << lti_pos << endl;
    NCout << "The whole job                          : " << ltim[6] - ltim[0] << endl;

    cout << "\nTime Statistic for the job (in second):" << endl;
    cout << "Input ............................ " << ltim[1] - ltim[0] << endl;
    cout << "Address computation .............. " << ltim[2] - ltim[1] << endl;
    cout << "Assembling the equation system ... " << lti_eqa << endl;
    cout << "Solving the equation system ...... " << lti_sol << endl;
    cout << "Post processing .................. " << lti_pos << endl;
    cout << "The whole job .................... " << ltim[6] - ltim[0] << endl;
  }

  // delete
  delete [] extNumbersOfNodes;
  delete [] isNodeMeshOrEvalNode;
  delete [] extNumbersOfElements;

  for(i=0; i<numNodes_; i++) delete [] nodesCoordinates[i];
  delete [] nodesCoordinates;

  for(i=0; i<numElements_; i++) delete [] elementsConnectivity[i];
  delete [] elementsConnectivity;

  delete [] listNumberNodesPerElement;
  delete [] listElementProperty;
  delete [] listElementsElementGroup;

  delete [] listSymmetryPlanes;
  delete [] listSymmetryPlanesCoordinates;

  delete [] ibval;
  for(i=0; i<numElements_; i++) delete [] zbval0[i];
  delete [] zbval0;
  for(i=0; i<numElements_; i++) delete [] zbvao0[i];
  delete [] zbvao0;
  delete [] zbval1;
  delete [] zbval2;
  delete [] zbvao1;
  delete [] zbvao2;

  for(i=0; i<numElements_; i++) delete [] elenor[i];
  delete [] elenor;
  for(i=0; i<numElements_; i++) delete [] centel[i];
  delete [] centel;
  delete [] areael;
  for(i=0; i<numElements_; i++) delete [] jelist[i];
  delete [] jelist;

  delete [] indexOfElementGroup;
  delete [] numberOfElementsInGroup;
  delete [] propertyOfGroup;

  delete [] numinw;
  delete [] numposo;
  for(i=0; i<numIncidentPlaneWaves_; i++) delete [] dirinw[i];
  delete [] dirinw;
  for(i=0; i<nptsouc; i++) delete [] coorps[i];
  delete [] coorps;
  for(i=0; i<numIncidentPlaneWaves_; i++) delete [] inwacurv[i];
  delete [] inwacurv;
  for(i=0; i<nptsouc; i++) delete [] posocurv[i];
  delete [] posocurv;
  delete [] zinwav;
  delete [] zposoc;
  delete [] zinwao;
  delete [] zposoo;

  for(i=0; i<numElements_; i++) delete [] ibvcurv[i];
  delete [] ibvcurv;
  delete [] nucurv;
  delete [] npcurv;
  for(i=0; i<numCurvesFrequency_; i++) delete [] frqcurv[i];
  delete [] frqcurv;
  for(i=0; i<numCurvesFrequency_; i++) delete [] faccurv[i];
  delete [] faccurv;

  delete [] Freqs;
}

// output informations
void NC_FrequencyInformations
(
	ostream& xout,
	ofstream& yout
)
{
	int i, n;
	double r_min, r_max;

	switch(methodFMM_)
	{
	case 0:
		xout << "\nTraditional BEM" << endl;
		break;
	case 1:
		xout << "\nSingle level fast multipole BEM" << endl;
		break;
	case 3:
		xout << "\nDirect multilevel fast multipole BEM" << endl;
		break;
	}

	// write the parameters of the equation system
	if(currentFrequency_ == istart_)
	{
		xout << "\nNumber of equations = " << numRowsOfCoefficientMatrix_  << endl;
		if(methodFMM_ == 0) {
			xout << "Number of entries of the coefficient matrix = " << numComponentsOfCoefficientMatrix_ << endl;
		}
	}

	// write infomations of clusters
	if(methodFMM_) {
		xout << "\nInformations about clusters:" << endl;
		xout << "   Level    Num. of clusters    Maximum radius   Minimum radius" << endl;
		for(i=0; i<numClusterLevels_; i++) {
			if(methodFMM_ == 1) {
				n = numOriginalClusters_;
				r_max = maxClusterRadiusBE_;
				r_min = minClusterRadiusBE_;
			} else {
				n = clulevarry[i].nClustOLv;
				r_max = clulevarry[i].RadiMaxLv;
				r_min = clulevarry[i].RadiMinLv;
			}
			xout << "   " << setw(3) << i << "       " << setw(6) << n
				<< "                " << setw(6) << r_max
				<< "         " << setw(6) << r_min << endl;
		}
		if(numInternalPointsClusters_) {
			xout << " resu. nets" << "  " << setw(6) << numInternalPointsClusters_ <<
				"                " << setw(6) << maxClusterRadiusRM_  <<
				"         " << setw(6) << minClusterRadiusRM_ << endl;
		}
	}

	switch(methodFMM_)
	{
	case 0:
		yout << "\nTraditional BEM" << endl;
		break;
	case 1:
		yout << "\nSingle level fast multipole BEM" << endl;
		break;
	case 3:
		yout << "\nDirect multilevel fast multipole BEM" << endl;
		break;
	}

	// write the parameters of the equation system
	if(currentFrequency_ == istart_)
	{
		yout << "\nNumber of equations = " << numRowsOfCoefficientMatrix_  << endl;
		if(methodFMM_ == 0) {
			yout << "Number of entries of the coefficient matrix = " << numComponentsOfCoefficientMatrix_ << endl;
		}
	}

	// write infomations of clusters
	if(methodFMM_) {
		yout << "\nInformations about clusters:" << endl;
		yout << "   Level    Num. of clusters    Maximum radius   Minimum radius" << endl;
		for(i=0; i<numClusterLevels_; i++) {
			if(methodFMM_ == 1) {
				n = numOriginalClusters_;
				r_max = maxClusterRadiusBE_;
				r_min = minClusterRadiusBE_;
			} else {
				n = clulevarry[i].nClustOLv;
				r_max = clulevarry[i].RadiMaxLv;
				r_min = clulevarry[i].RadiMinLv;
			}
			yout << "   " << setw(3) << i << "       " << setw(6) << n
				<< "                " << setw(6) << r_max
				<< "         " << setw(6) << r_min << endl;
		}
		if(numInternalPointsClusters_) {
			yout << " resu. nets" << "  " << setw(6) << numInternalPointsClusters_ <<
				"                " << setw(6) << maxClusterRadiusRM_  <<
				"         " << setw(6) << minClusterRadiusRM_ << endl;
		}
	}
}
