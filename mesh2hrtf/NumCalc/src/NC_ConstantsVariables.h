/*===========================================================================*\
 *                                                                            *
 *  File name:      NC_ConstantsVariables.h                                   *
 *  Description:    Definitins of constants and variables                     *
 *  Author:         H. Ziegelwanger, W. Kreuzer and Z.S. Chen                 *
 *                                                                            *
 \*===========================================================================*/



#if !(defined(__APPLE__) || defined(__MACH__) || defined(unix) || defined(__unix__) ||defined(__unix))
#define _CRT_SECURE_NO_WARNINGS
#define _SCL_SECURE_NO_WARNINGS
#endif

#ifndef NC_ConstantsVariables_h
#define NC_ConstantsVariables_h



//================================= I N C L U D E S ==================================
// local includes                                                                   //
#include "NC_TypeDefinition.h"                                                      //
//                                                                                  //
// system includes                                                                  //
#include<stdio.h>                                                                   //
#include<stdlib.h>                                                                  //
#include<string.h>                                                                  //
//====================================================================================



/*
	defines for windows system
*/
#if !(defined(__APPLE__) || defined(__MACH__) || defined(unix) || defined(__unix__) ||defined(__unix))
#define isWindows
#endif

/*
	defines for input
*/
#define WHITE "\t \n\r"
#define DIGITS "0123456789.+-"
#define SIZE 300
#define NTRM 30

/*
	constants used in the program
*/
#define NDIM 3
#define NETYP3 3
#define NETYP4 4
#define SIZE_LINE 200
#define NNODPE 4
#define MSBE 220
#define PI 3.1415926535897932385
#define PI4 12.566370614359172944
#define PI2 6.283185307179586477
#define EPSY 1.0e-14

/*
	global variables
*/
extern ofstream NCout;  /* output stream */

extern FILE *inputFile_, *tmpFile_;   /* input and output files */

//extern char jobTitle_[SIZE_LINE];   /* title of the job */
extern string jobTitle_;
extern string versionNumber_;     /* version number */

extern int currentFrequency_,   /* number of the current frequency */
		   numFrequencies_;     /* number of frequencies */

extern int istart_;            /* index of the first frequency step */

extern double firstFrequencyStep_,      /* first frequncy step */
			  frequencyStepIncrement_;  /* frequency step increment */

extern int numElementGroups_,       /* number of element groups */
		   numElements_,            /* number of elements */
		   numNodes_,               /* number of nodes */
		   minExpansionTermsFMM_,   /* minimum number of terms in the FMM expansion */
	       numSymmetricPlanes_,     /* number of symmetric planes */
	       numIncidentPlaneWaves_,  /* number of the incident plane waves */
	       numPointSources_,        /* number of the point sources */
           isInternalProblem_;		/* if the job is internal problem */

extern double speedOfSound_,        /* sound speed in the medium */
  densityOfMedium_,     /* density of the medium */
  harmonicTimeFactor_,              /* factor for the unit of the imaginary numbers 
				       ( =  1.0  when G = exp(i k r)/(4 PI r)
				       = -1.0  when G = exp(-i k r)/(4 PI r) ) */
  waveNumbers_,         /* wave numbers */
  waveLength_,          /* wave length */
  averageElementArea_,  /* the average area of boundary elements */
  frequency_,           /* frequency */
  omega1_,              /* 2*PI*frequency */
  rpfact_,              /* sound pressure = rpfact * velocity potential */
  Tao_,                 /* = 1.0 for external problem; = -1 for internal problem */
  delta_,               /* a small distance, dependent to mesh */
  epsilon_,             /* a small distance (= 1.0e-3 mm) */
  ClusEdgL0_;           /* input value of average edge length of the clusters 
			   (of the coarsest level for MLFMBEM), if <= 0, a value
			   for average edge length of the clusters will be computed 
                                     by the program*/

extern int numNodesOfBoundaryMesh_,                 /* number of nodes of the boundary element mesh */
  numNodesOfEvaluationMesh_,               /* number of nodes of the evaluation mesh */
  numElementsOfBoundaryMesh_,              /* number of elements of the BE mesh */
  numElementsOfEvaluationMesh_,            /* number of elements of the evaluation mesh */
  numRowsOfCoefficientMatrix_,             /* number of rows of the coefficient matrix */
  numComponentsOfCoefficientMatrix_,       /* number of components of the coefficient matrix */
  numReflectionsOfElements_,               /* number of reflections of the elements with respect to the 
					      symmetric planes */
  numCurvesFrequency_;                     /* number of curves defining the dependance on freuqancy */

extern int methodBEM_,                      // input value for mehod to be used
// = 0: Traditional BEM
// = 1: Single level FMBEM (SLFMBEM)
// = 4: Direct MLFMBEM without interpolation/filtering procedures
//      (by storing all coef. matrices on all levels)
  numNonZeroEntries_,              // number of nonzero entries of the near field matrix
  numClusterLevels_,               // number of cluster levels
  nlevtop_,                        // = numClusterLevels_ - 1;
  maxRowNumberD_,                  // maximum row number of the D-matrix in each level
  numOriginalClusters_,            // number of original clusters of boundary elements
  numOriginalReflectedClusters_,   // sum of original clusters and reflected clusters
  numInternalPointsClusters_,      // number of clusters of nodes of the internal points
  methodPreconditioner_,           // flag for the preconditioners
// = -2: incomplete LU-decomposition preconditioner,
//      second order scanning degree
// = -1: incomplete LU-decomposition preconditioner,
//      first order scanning degree
// = 0: incomplete LU-decomposition preconditioner, 
//      zero order scanning degree
// = 1: scaling preconditioner
// = 2: no preconditioner
  scanningDegreeLU_,               // scanning degree for LU-decomposition preconditioning
// = 0: zero order scanning degree (fewer remainig components)
// = 1: first order scanning degree (more remainig components)
// = 2: second order scanning degree (most remaining components)
  methodSolver_,                   // = 0: CGS method (A. Meister p 168)
// = 1: QMRCGSTAB method (A. Meister p 186)
// = 2: stabilized BiCGSTAB method (Andreas Meister P 173)
// = 3: BiCGSTAB method (A. Meister P 172)
// = 4: Gaussean elimination method (direct method)
  niter_max_;
extern double farFieldClusterFactor_,   // distance factor used to define far field clusters
			  minClusterDistance_,      // minimum distance between two interacting clusters 
                                        // (in order to guarantee that the spherical Hankel function in (18) 
                                        // Chen et. al. is not too big)
			  maxClusterRadiusBE_,      // maximum radius of clusters of boundary elements
			  avgClusterRadiusBE_,      // average radius of clusters of boundary elements
			  minClusterRadiusBE_,      // minimum radius of clusters of boundary elements
			  maxClusterRadiusRM_,      // maximum radius of clusters of the result mesh
			  avgClusterRadiusRM_,      // average radius of clusters of the result mesh
			  minClusterRadiusRM_;      // minimum radius of clusters of the result mesh

extern int numExpansionTerms_;      // number of terms of the expansion expresions of the Green
                                    // functions

extern int numIntegrationPointsUnitSphere_,     // number of integration points on the unit sphere surface
           numIntegrationPointsThetaDirection_, // number of integral points in the thetha direction
           numIntegrationPointsPhiDirection_;   // number of integral points in the phi direction

extern int methodFMM_;    // = 0: traditional BEM (TBEM)
                          // = 1: single level FMBEM (SLFMBEM)
						  // = 3: multi level FMBEM (MLFMBEM)

extern bool	boolComputeTVector_;   // = true:  the T-vector for postprocessing should be computed 
                                   // = false: it is not necessary

extern double diameterBE_; // diagonal length of the box closing the BE-mesh

extern string fileMeshBoundarySurface_,
              fileMeshEvaluationGrid_,
              fileWavesPlanar_,
              fileWavesSpherical_;


#endif /* Don't add any thing after this line! */
