/*===========================================================================*\
 *                                                                            *
 *  File name:      NC_Arrays.h                                               *
 *  Description:    Definitins of vectors and matrices                        *
 *  Author:         H. Ziegelwanger, W. Kreuzer and Z.S. Chen                 *
 *                                                                            *
 \*===========================================================================*/



#ifndef NC_Arrays_h
#define NC_Arrays_h 



//================================= I N C L U D E S ==================================
// local includes                                                                   //
#include"NC_TypeDefinition.h"                                                       //
//                                                                                  //
// system includes                                                                  //
//====================================================================================



// define of global arrays
extern int *extNumbersOfNodes,       // external numbers of the nodes
           *isNodeMeshOrEvalNode,     // array indicating if a node is a node of the be mesh (= 1)
		                  // or of the evaluation mesh ( = -1)
		   *extNumbersOfElements;       // external numbers of the elements
extern double **nodesCoordinates;    // coordinates of the nodes
extern int **elementsConnectivity;      // connectivities of the elements
extern int *listNumberNodesPerElement,       // number of nodes of each element
           *listElementProperty,       // flag indicating the property of each element 
		                  //     = 0: surface el.
						  //	 = 2: evaluation el.
		   *listElementsElementGroup;       // number of the element group to which an element belongs
extern int *listSymmetryPlanes;       // array for symmetric or reflecting planes
                          // = 0: no symmetric or reflecting plane normal to the corresponding
                          //      coordinate axis
                          // = 1: symmetric plane
                          // = 2: antisymmetric plane   (converted to -1 in the program)
                          // = 3: hard reflecting plane (converted to 2 in the program)
                          // = 4: soft reflecting plane (converted to -2 in the program)
						  // = 5: admmitance/impedance of the ground surface (the infinite 
                          //      plane at z = 0.0) is prescribed 
                          //      (converted to 3 in the program)
extern double *listSymmetryPlanesCoordinates;    // coordinate of the intersection point of the symmetric (or
                          // reflecting) plane with the corresponding axis
extern int *ibval;        // array indicating the kind of the boundary condition for each
                          // element
                          // (= 0: velocity;
                          //  = 1: pressure;
                          //  = 2: velocity and surface admittance;
                          //  = 3: transfer admittance
                          //  = 4: transfer admittance and surface admittance)
extern Complex **zbval0,  // boundary conditions for the faces under the absorbing materials
			   **zbvao0,  // input values for zbval0 (frequency dependancies are not yet 
                          // considered)
			   *zbval1,   // boundary conditions at the positive sides
			   *zbval2,   // boundary conditions at the negative sides
			   *zbvao1,   // input values for zbval1 (frequency dependancies are not yet
						  // considered)
			   *zbvao2;   // input values for zbval2

extern double **elenor,   // normal vectors to each element at the centre and nodes
			  **centel,   // coordinates of the centres of the elements
			  *areael;    // area or length of each element

extern int    **jelist;   // addresses of the unknowns of each element in the global
                          // equation system

extern int    *indexOfElementGroup,  // number of each element group
			  *numberOfElementsInGroup,  // number of elements of each element group
			  *propertyOfGroup;  // property of each element group
						  // = 0: surface elements
                          // = 2: evaluation elements

extern int    *numinw,    // external numbers of the incident plane waves
			  *numposo;   // external numbers of the point sources (monopol)
extern double **dirinw,   // direction vectors of the incident plane waves
			  **coorps;   // coordinates of the point sources
extern int    **inwacurv, // numbers of curves to prescribe the dependancy of the incident
                          // plane waves on the frequency
			  **posocurv; // numbers of curves to prescribe the dependancy of the point
                          // soures on the frequency
extern Complex *zinwav,   // intensities of the incident plane waves (dependent on the
						  // frequncies)
			   *zposoc,   // intensities of the point sources (dependent on the frequncies)
			   *zinwao,   // input value of zinwav (frequency dependence not considered)
			   *zposoo;   // input value of zposoc (frequency dependence not considered)

extern int    **ibvcurv,  // numbers of the curves prescribibing the dependency of the 
                          // boundary conditions on the frequencies for each element
			  *nucurv,    // numbers of the cueves prescribibing the dependency of the
                          // boundary conditions on the frequencies
			  *npcurv;    // number of points of each such curve
extern double **frqcurv,  // frequencies corresponding to each point of each curve
              **faccurv;  // multiplicators corresponding to each point of each curve

extern ElCluster *clustarry; // array of clusters of boundary elements
extern IpCluster *ipcluarry; // array of clusters of internal points
extern int *jcolnea,      // column number of each nonzero of the near field matrix
           *irownea;      // number of the first nonzero of each row of the near field 
                          // matrix
extern int *jcoltmtx,     // column number of each nonzero of the T-matrix
           *irowtmtx;     // number of the first nonzero of each row of the T-matrix
extern Complex *ztmtx;    // array for nonzeros of the T-matrix
extern Complex *ztvct;    // T-vector to store the far field transform of the contributions
                          // of velocity boundary conditions to the right hand side
extern D_mtx_lev *dmtxlev; // D-matrix at each level
extern T_mtx_lev *tmtxlev; // D-matrix at each level
extern S_mtx_lev *smtxlev; // D-matrix at each level
extern int *jcolsmtx,     // column number of each nonzero of the S-matrix
           *irowsmtx;     // number of the first nonzero of each row of the S-matrix
extern Complex *zsmtx;    // array for nonzeros of the S-matrix
extern double **uvcsphe;  // untit vectors at the integarl points of the unit sphere
extern double *weisphe;   // weight of each integarl points of the unit sphere

extern Complex *zcoefl;   // coefficient matrix of BEM or near field matrix of FMBEM
extern Complex *zrhs;     // right hand side vector

extern ClusterLev *clulevarry; // array of cluster levels of a cluster tree
extern ElCluster *ClustArray;  // working array of clusters

// local coordinates of nodes of a 8/6-noded element
const double CSI8[8] = {1.0, -1.0, -1.0, 1.0, 0.0, -1.0, 0.0, 1.0},
			 ETA8[8] = {1.0, 1.0, -1.0, -1.0, 1.0, 0.0, -1.0, 0.0},
			 CSI6[6] = {1.0, 0.0, 0.0, 0.5, 0.0, 0.5},
			 ETA6[6] = {0.0, 1.0, 0.0, 0.5, 0.5, 0.0};


#endif /* Don't add any thing after this line! */