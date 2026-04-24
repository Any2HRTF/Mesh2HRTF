#include "NC_ConstantsVariables.h"
#include "NC_TypeDefinition.h"
#include "NC_MLFMM.h"
#include <sys/resource.h>
#include <cmath>

using namespace std;

/*
  Collection of functions for implementing the (ML)FMM

  (Semi)global matrices:
  zMmat: Allocated in cluster2clustermat, contrains all factors for
         the interaction between interaction clusters
	 zMmat[level][[clusterinlevel][interactionlist*nodesphere]

  zFmat: Contains the local expansions on the leaf level and the up pass values
         zFmat[ leafclusters ][pointssphere][numRowsOfCoefficientMatrix_]

  zFvec: zFvec[ leafclusters ][pointsphere]

  zSmat: zSmat[clusters at leaf][elems in cluster][points sphere]
  
*/
extern Complex NC_IncidentWaveRHS( ofstream&);
extern void NC_SingularIntegration(ofstream&,Vector<Complex>&,const int&,const int&,Vector<Complex>&,Matrix<double>&);
extern void NC_RegularIntegration(ofstream&,Vector<Complex>&,const int&,const int&,const int&,const int&,Vector<Complex>&,Vector<double>&,Matrix<double>&);
Complex NC_IncidentWaveRHS(ofstream&);


extern bool adapt_fmmlength_;
void Cleanup_MLFMM(bool deleteF) {
  /* frees the memory allocated by the 4 matrices */
  /* Global: zMmat, zFmat, zSmat, zNear */
  if( zMmat != NULL ) {
    for (int n = 0; n < numClusterLevels_; n++) {
      for (int i = 0; i < clulevarry[n].nClustOLv; i++)  {
	delete[] zMmat[n][i];
      }
      delete[] zMmat[n];
    }
    delete[] zMmat;
    zMmat = NULL;
  }

  if( zSmat != NULL ) {
    for (int i = 0; i < clulevarry[nlevtop_].nClustOLv; i++) 
      delete[] zSmat[i];
    delete[] zSmat;
    //zNear.~zSparsetype();
    zSmat = NULL;
  }
  
  if( deleteF && zFmat != NULL) {
    for( int i = 0; i < clulevarry[nlevtop_].nClustOLv; i++) {
      for ( int s = 0; s < clulevarry[nlevtop_].nPoinSpheLv; s++) {
	delete [] zFmat[i][s];
      }
      delete [] zFmat[i];
    }
    delete[] zFmat;
    zFmat = NULL;
  }

  if( zNear.zdata != NULL ) {
    zNear.delete_arrays();
  }
}


void LocalExpansionMat(int maxlevel, bool allocateFMM) {
  /*
    Calculation Routine for the local FMM expansion matrices at root level
    sum_{Gamma in C} int_Gamma e^{ik (y - z)^T s}d\Gamma x_i
    at leaf level, also calculate the right hand side contributions
    
    Variables
	        nclusters: Number of clusters at the level
	        Number quad nodes sphere: Number of quadrature nodes on the sphere
              Depends on the expansion length L
				    N = (L+1)^2
	maxlevel: max fmm level
	  
    The integral over Gamma can be calculated by discretizing
    y = v0 + e1 * xi + e2 * eta with xi in [0,1] and \eta in [0,xi]
    where v0 is the first vertex of the element, and e1,e2 are two edges
    (either triangle or quadrilateral)
    
    thus we have for the G part
    kappa_Gamma e^{ik (v0-z)*s) int_0^1 e^{ik e1*s xi} \int_0^{xi} e^{ik e2*s \eta}d\eta d\xi
    where kappa is the scaling from unit triangle/quadrilateral to Gamma
    this integrals can be calculated analytically
    
    using s1 = s*e1 and s2 = s*e2 the result(s) of the integral are

    s1 neq 0, s2 neq 0, s1 + s2 neq 0
    I =   (s1e^{iks1}(1 - e^{iks2}) - (1 - e^{iks1})s2) /
                                 k^2 ( s1s2^2 + s1^2 s2 )
    s1 = 0, s2 neq 0
    I = (iks2 + 1 - e^{iks2}) / (k^2 s2^2)
 
    s1 neq 0, s2 = 0
    I = ( (1 - iks1) e^{iks1} - 1 )/ (k^2 s1^2)
    
    s1 = 0, s2 = 0
    I = 1/2

    s1 + s2 = 0
    (iks1 + 1 - e^{iks1}) / ( k^2 s1^2 )

    
    Global variables:
       elementsConnectivity: int**, elementsConnectivity[elem][node] is the nodenumber of the vertices
       listNumberNodesPerElement: int* , listNumberNodesPerElement[elem] is the number of vertices of elem
       clulevarry:  Cluster info, input
          Variable (class) containing all the informaiton for the clustering
	  at the current level
	  contains:
	     Clust: pointer to ElCluster, and info about expansion lenght, points
	     on the sphere ....
	     use variables 
	     nPoinSpheLv: int: number of quadrature nodes on the spere
	     nExpaTermLv: int: expansion length on current level
	     uvcsphe: **double: coordinates of the quadrature nodes on the
	                        sphere
	     Elcluster: struct containing info about the cluster itself, e.g.
	         CoordCent: Coordinates of the cluster center
		 NumNeaClus: number of near field clusters
		 NumFarClus: interaction list
       zFmat: out: array containing complex double matrices
          dimension: number of clusters times number of quadrature nodes times number of BEM nodes
       zFvec: out: array of vectors zFmat * bcval, this has to be done
               for the rhs only once, and not every itertation,
	       thus it can be done here
              size: nclusters times Number quad nodes sphere times elements in the cluster

    Global: zbval0[iel][i]  pres or velo condition for i-th vertex for iel
	    zbval1[iel]	 admittance bc
	    ibval[iel]

   Notes: Looks more complicated than it is, but we have to consider the
          different boundary conditions and the burton miller method
	  thus, different combinations of the local expansions with or
	  without the derivative with respect to y

	  The matrix is frequency dependent, thus it needs to be calculated for
	  every frequency. However, if the clustering does not change between
	  frequencies, there is no need to allocate everything again.
	  maxlevel should be nlevtop_ = numClusterLevels_ - 1

	  It should not matter if we look at r = ||y - x|| or r = ||x - y||
	  because of the integral with respect to the unit sphere
	  There might be some trouble if the nodes on the sphere are not
	  symmetric

	  lets take the (y - x) = y - z1 + z1 - z2 + z2 - x version, just
	  to be on the same page with the nearfield calculation, however
	  do not forget that chen changes the sign on the whole BIE
   Written: Kreiza, Nov. 2024
  */

  int i,j,s,n,nn,iv,iv1,iv2;
  double v0[NDIM],e1[NDIM],e2[NDIM];
  double z0[NDIM];
  double s0,s1,s2,kappa,rval;
  Complex zs0,zs1,zs2,zval,bcval,admival;
  int Ibvj03 = 0;  
  bool Ifadmij = false;
  bool nonzerobc;
  int Gamma_j,nvert;
  int tentries = 0;  
  
  if( allocateFMM ) {
    /* ****************************************
                         zFmat
      ********************************* */

    zFmat = new Complex**[clulevarry[maxlevel].nClustOLv];
    if( zFmat == NULL ) {
      cerr << "Sorry could not allocate zFmat\n";
      exit(-1);
    }
    for ( i = 0; i < clulevarry[maxlevel].nClustOLv; i++ ) {
      zFmat[i] = new Complex*[ clulevarry[maxlevel].nPoinSpheLv ];
      
      if( zFmat[i] == NULL ) {
	cerr << "Sorry could not allocate zFmat\n";
	exit(-1);
      }
      for ( s = 0; s < clulevarry[maxlevel].nPoinSpheLv; s++ ) {
	zFmat[i][s] = new Complex[ clulevarry[maxlevel].ClustArLv[i].NumOfEl ];
	tentries += clulevarry[maxlevel].ClustArLv[i].NumOfEl;	
	if( zFmat[i][s] == NULL ) {
	  cerr << "Sorry could not allocate zFmat\n";
	  exit(-1);
	}	
      }
     
    }

    /* *******************************************************
                          zFvec
       ******************************************************** */
    // lets waste a little memory for efficiency
    zFvec.nonzeroblocks = 0;
    //    zFvec.clusterindx = new int[clulevarry[maxlevel].nClustOLv];
    for ( i = 0; i < clulevarry[maxlevel].nClustOLv; i++ ) {
      if( clulevarry[maxlevel].ClustArLv[i].IfNonZeroBc )
	zFvec.nonzeroblocks++;
    }
    if( zFvec.nonzeroblocks > 0 ) {
      zFvec.clusterindx = new int[zFvec.nonzeroblocks];
      if( zFvec.clusterindx == NULL) {
	cerr << "Sorry could not allocate zFvec\n";
	exit(-1);
      }
      zFvec.zdata = new Complex[zFvec.nonzeroblocks * clulevarry[maxlevel].nPoinSpheLv];
      if( zFvec.zdata == NULL) {
	cerr << "Sorry could not allocate zFvec\n";
	exit(-1);
      }
    }
  }

  cout << "Matrix T:\n";
  cout << "Complex: " << tentries << "\n";
  int nsphere = clulevarry[maxlevel].nPoinSpheLv;
  int zFvecblock = 0;  
  // loop over all leaf clusters
  for (n = 0; n < clulevarry[maxlevel].nClustOLv; n++) { 
    // coordinates of the cluster midpoint
    for (j = 0; j < NDIM; j++) 
      z0[j] = clulevarry[maxlevel].ClustArLv[n].CoorCent[j];  
    
    // loop over all elements in the cluster
    for (i = 0; i < clulevarry[maxlevel].ClustArLv[n].NumOfEl; i++) {
      // look at each single element in current cluster
      Gamma_j = clulevarry[maxlevel].ClustArLv[n].NumsOfEl[i];
      // set the flags for the bc
      // if zbval0 is > 0 at this element, thats a contribution to the rhs
      switch( ibval[Gamma_j] ) {
      case 0:         // velocity prescribed
	Ibvj03 = 0;  
	Ifadmij = false;
	break;
      case 1:         // pressure prescribed
	Ibvj03 = 1;
	Ifadmij = false;
	break;
      case 2:         // velocity and surface admittance prescribed
      case 5:
	Ibvj03 = 0;  
	Ifadmij = true;
	break;
      default:
	// should never happen, but just in case
	Ifadmij = false;
	Ibvj03 = 0;
      }


      // this could be made more efficient, if it is checked at the input
      // if there is a bc != 0 at all, if not, we would not need zFvec
      nonzerobc = false;
      // if all are false, we have a sound hard surface
      // zbval0 is defined for all vertices
      if( zbval0[Gamma_j][0].norm() > EPSY ) {
	bcval = zbval0[Gamma_j][0];
	nonzerobc = true;
      }
      if( Ifadmij ) {
	admival = zbval1[Gamma_j];
      }
	
      // get vertices of each element and two edges
      nvert = listNumberNodesPerElement[Gamma_j];
      // scaling from unit element to global element
      kappa = areael[Gamma_j];
      if( nvert == 3)
	kappa *= 2.0;
      
      if( nvert > 4 ) {
	cout << "Sorry not implemented yet \n";
	exit(-1);
      }
      // for the analytic solution of int_\Gamma_j e^{ik (z1-y)s}
      // as function of s 
      // we need one vertex and two edges of the triangle
      // chen uses r = |y - x| in the greens function, but here
      // we go to the more regular r = |x-y| formulation
      iv = elementsConnectivity[Gamma_j][0];
      iv1 = elementsConnectivity[Gamma_j][1];
      iv2 = elementsConnectivity[Gamma_j][2];
      // sign changed for debugging purposes  
      for( nn = 0; nn < NDIM; nn++) {
	v0[nn] = -(nodesCoordinates[iv][nn] - z0[nn]);  // vertex 0 - z0
	e1[nn] = -(nodesCoordinates[iv1][nn] - nodesCoordinates[iv][nn]);  // edge 1
	e2[nn] = -(nodesCoordinates[iv2][nn] - nodesCoordinates[iv1][nn]); // edge 2
      }  
            
    // loop over the quadrature nodes of the sphere
      for ( s = 0; s < nsphere; s++ ) {
	s0 = 0.0;
	s1 = 0.0;
	s2 = 0.0;

	for( nn = 0; nn < NDIM; nn++) {
	  s0 += clulevarry[maxlevel].uvcsphe[s][nn] * v0[nn];
	  s1 += clulevarry[maxlevel].uvcsphe[s][nn] * e1[nn];
	  s2 += clulevarry[maxlevel].uvcsphe[s][nn] * e2[nn];
	}

	s0 *= waveNumbers_;
	s1 *= waveNumbers_;
	s2 *= waveNumbers_;

	// look at 4 different cases
	zs0.set( -sin(s0), cos(s0) ); // -i * exp(ik (y - z)^T s)
	zFmat[n][s][i] = zs0 * kappa * waveNumbers_ / PI4;
	if( fabs(s1) < EPSY ) {
	  if( fabs(s2) < EPSY ) {
	    zFmat[n][s][i] *= -0.5;
	  }
	  else {
	    zs1.set( cos(s2), sin(s2) );
	    zs2.set(1.0, s2);
	    zFmat[n][s][i] *= -(zs2 - zs1)  ;
	    zFmat[n][s][i].div_r(s2 * s2);
	  }
	}
	else { //s1 != 0
	  if( fabs(s2) < EPSY ) {
	    zs1.set( cos( s1 ), sin(s1) );
	    zs2.set(1.0,-s1); // 1 - is1
	    zFmat[n][s][i] *= -(zs2 * zs1 - 1.0);
	    zFmat[n][s][i].div_r(s1*s1);
	  }
	  else if( fabs(s1 + s2) < EPSY ) {
	    zs1.set(cos(s1),sin(s1));
	    zs2.set(1.0,s1);
	    zFmat[n][s][i] *= -(zs2 - zs1);
	    zFmat[n][s][i].div_r(s1 * s1);
	  }
	  else {
	    // the exponential faktors exp(i(z_0 - y)\cdot s)
	    // minus because we need e^{ik( y - z ) }
	    zs1.set( cos(s1), sin(s1) );
	    zs2.set( cos(s1+s2), sin(s1+s2) );
	    // overwrite s0 to s2 with the denominator parts,
	    // wavenumber already included
	    // now that we have int_\Gamma_j e^{ik(z0 -y ) s}, for every sphere
	    // quadrate node, lets apply it
	    // depending on the boundary condition on Gamma_j we have different
	    // combinations of G,H,H', and E, plus for a nonzero bc value the
	    // integral is a contribution for the rhs not for the system matrix
	    
	    // no matter what bc we will always have at least one exp part in
	    // the system matrix either with or without derivative with respect
	    // to y:
	    // Dirichlet: Matrix: G + zBta3 * H', rhs: -p0/2 + Hp0 + zBta3 E p0
	    // Neumann: Matrix: H + zBta3 * E,  rhs: -zBta3 v0/2 - (G - zbta3 H')v0 
	    if( nvert == 3 ) {
	      zFmat[n][s][i] *= ( zs2 * s1  - (zs1 - 1.0) * s2 - zs1 * s1 );
	      zFmat[n][s][i].div_r(s1 * s2 * s2 + s1 * s1 * s2);
	    }
	    else if( nvert == 4) {
	      cerr << "Sorry not implemented yet\n";
	      exit(-1);
	    }
	    else {
	      cerr << "Sorry too many vertices/nodes in the element\n";
	      exit(-1);
	    }
	  } // if s2 == 0
	} // if s1 == 0
	// check the bc and adapt accordingly
	if( Ibvj03 == 1 ) {
	  // pressure boundary condition
	  // system consists of G and H', no derivative necessary
	  // rhs possibly on H and E

	  // matrix is (e_x + bta e'_x) M_L e_y, so no need to change anything
	  // for the matrix, but the rhs needs some changing
	  if( nonzerobc ) {
	    // H + bta E = (e_x + beta e'_x) M_L e'_y
	    // get the normal vector times the quadnode
	    rval = 0.0;
	    for (j = 0; j< NDIM; j++)
	      rval -= elenor[Gamma_j][j] * clulevarry[maxlevel].uvcsphe[s][j];
	    rval *= waveNumbers_;
	    bcval.mul_i(rval);
	    zFvec.zdata[zFvecblock * nsphere + s] +=  zFmat[n][s][i] * bcval;
	    if( s == 0 )
	      zFvec.clusterindx[zFvecblock] = n;
	    if( s == nsphere - 1)
	      zFvecblock++;
	  }
	}
	else {
	  // velocty bc
	  // we need H and E for the system and
	  //    G and H' for the rhs
	  // if we also have and admittance we need all of them
	  rval = 0.0;
	  for (j = 0; j< NDIM; j++)
	    rval += elenor[Gamma_j][j] * clulevarry[maxlevel].uvcsphe[s][j];
	  rval *= waveNumbers_; 
	  zval.set(0.0,rval);  // ikn
	  if( nonzerobc ) {
	    // admittance is not necessary here
	    // (G + \beta H')v_0, the beta factor comes later
	    zFvec.zdata[zFvecblock * nsphere + s] += zFmat[n][s][i] * bcval;
	    if( s == 0 )
	      zFvec.clusterindx[zFvecblock] = n;
	    if( s == nsphere - 1 )
		zFvecblock++;
	  }
	  if( Ifadmij )
	    // (e_x + beta e'_x) M_L ( alpha e_y +  e'_y )
	    // the plus for rval is because one negative sign from the
	    // derivative, the second from the BIE
	    zFmat[n][s][i] *= (zval  + admival);
	  else
	    zFmat[n][s][i] *= zval;
	}  // if bc condition
      } // loop over quad nodes
    } // loop elemes
  } // loop over all clusters in the level
}

#if 0

void apply_localExpansion(Complex* zF, Vector<Complex>& velopot, Vector<Complex>& partvelo) {
  /* calculates the local multipole expansions for the evaluation grid
     zF  inout  Complex, number of clusters times number of nodes on sphere
     velopot in   velocity potential on each boundary element
     partvelo in  particle velocity on each boundary element

     Global:
     zFmat: in, zFmat[cluster][number of elems in cluster][number nodes sphere]
         local expansion matrix on leaf level
     zFvec: 
  */
  int nsphere = clulevarry[nlevtop_].nPoinSpheLv;
  int Gamma_j;
  for( int n = 0; n < clulevarry[nlevtop_].nClustOLv; n++) {
    // local expansion in the leaf clusters
    for ( int s = 0; s < nsphere; s++) {
      zF[n*nsphere + s].set(0.0,0.0);
      for ( int j = 0; j < clulevarry[nlevtop_].ClustArLv[n].NumOfEl; j++) {
	Gamma_j  = clulevarry[nlevtop_].ClustArLv[n].NumsOfEl[j];
	// check the bc condition
	if( ibval[Gamma_j] == 0 || ibval[Gamma_j] == 2)
	  ztmp = velopot[Gamma_j];
	else
	  ztmp = partvelo[Gamma_j];
	
	zF[n*nsphere + s] += zFmat[n][s][j] * ztmp;
      }
    }
  }
  // add the velocity conditions PRES and VELO
  // ADMI has been already handled above
  for (int n = 0; n < zFvec.nonzeroblocks; n++) {
    C_i = zFvec.clusterindx[n];
    for( int s = 0; s < nsphere; s++) {
      zF[n*nsphere + s] += zFvec.zdata[n*nsphere + s];
    }
  }
}
#endif


void apply_localExpansion(Complex*** zF, Vector<Complex>& ztmp) {
  // applies the local FMM expansions at leaf level
  /* zFmat: in, zFmat[cluster][number of elems in cluster][number nodes sphere]
                local expansion matrix on leaf level
     zF:    out,   zFmat * x, where x is the approximation of the solution
                   in the iterative solver
		   
     ztmp:  in: [number of elements in mesh] contains the
                       "solution" vector x for all collocnodes

     Global
     zFmat: in, zFmat[cluster][number of elems in cluster][number nodes sphere]
                local expansion matrix on leaf level
  */
  int nsphere = clulevarry[nlevtop_].nPoinSpheLv;
  int n,s,j, Gamma_j;
  int maxGamma, minGamma;
  /*  maxGamma = -1;
  minGamma = 10000;
  */
  for( n = 0; n < clulevarry[nlevtop_].nClustOLv; n++) {
    // local expansion in the leaf clusters
    for ( s = 0; s < nsphere; s++) {
      zF[nlevtop_][n][s].set(0.0,0.0);
      for ( j = 0; j < clulevarry[nlevtop_].ClustArLv[n].NumOfEl; j++) {
	Gamma_j  = clulevarry[nlevtop_].ClustArLv[n].NumsOfEl[j];
	/*	if( Gamma_j > maxGamma )
	  maxGamma = Gamma_j;
	if( Gamma_j < minGamma )
	minGamma = Gamma_j;*/
	zF[nlevtop_][n][s] += zFmat[n][s][j] * ztmp[Gamma_j];
      }
    }
  }
}


void Get_Interpolation_Matrices(double** Ylev, int nlevels) {
  /* Calculates the interpolation matrix for the up(down) pass used between each level
     On each child level the interpolation matrix can be constructed by
     T_{ij} = sum_\ell (2l + 1)/(4pi) P_l( y_i x_j ),
     where the x_j are the quadrature nodes on the "children" sphere
     and y_i are the  nodes not the parent sphere
     currently the quadrture weights need to be drawn to the input vector 

     YLev: output: pointer to a Matrix: each matrix contains the
                   T_ij mentioned above
     
     nlevels: input: number of levels


     Warning: As arrays in c++ are in general from 0 to something, the
     Interpolation matrix from level n to level n - 1, will be at
     Ylev[n-1]

     Note: The multiplication with the quadrature weights on the sphere is
           done later for efficiency reasons
     It is assumed that the gsl is compiled with this code
     Written by kreiza: 19.11.24
     
     Global
     clastarry: Clustertree for all levels


  */
  int i,j,n,nn,lev,l;
  int currentL,parentL,Ncurrent,Nparent;
  double** QNodescurrent;
  double** QNodesparent;
  double Angle;

  // loop over multipole levels
  for ( lev = nlevels; lev > 0; lev-- ) {

    //if( lev == 1 )
    if( clulevarry[lev-1].nClustOLv == 1 ) //|| clulevarry[lev-1].isNearClust)
	continue;
    
    // expansion length of level and parent
    currentL = clulevarry[lev].nExpaTermLv;
    parentL = clulevarry[lev - 1].nExpaTermLv;
    Ncurrent = clulevarry[lev].nPoinSpheLv;
    Nparent = clulevarry[lev-1].nPoinSpheLv;
    
    QNodescurrent = clulevarry[lev].uvcsphe;
    QNodesparent = clulevarry[lev-1].uvcsphe;

    /*
    for( i = 0; i < Ncurrent; i++ ) {
      cout << clulevarry[0].uvcsphe[i][0] << " " << clulevarry[0].uvcsphe[i][1] << " " << clulevarry[0].uvcsphe[i][2] << "\n";
    }
    cout << "\n";
    for ( i = 0; i < Ncurrent; i++ )
      cout << clulevarry[0].weisphe[i] << "\n";

    exit(0);
    */
    
    if( currentL == parentL ) {
      // do nothing
      //      continue;
    }

    // loop over sphere points
    for( j = 0; j < Ncurrent; j++ )  { // nodes on the sphere 
      for ( i = 0; i < Nparent; i++ ) {
	Angle = 0.0;
	for ( n = 0; n < NDIM; n++ )
	  Angle += QNodescurrent[j][n] * QNodesparent[i][n];
	for ( l = 0; l < currentL; l++ ) {
	  if( l == 0 )
	    //	    Ylev[lev-1][i * Ncurrent + j] = 0.25/PI;
	    Ylev[lev-1][i * Ncurrent + j] = 1.0;
	  else
	    //	    Ylev[lev-1][i * Ncurrent + j] += (2.0 * double(l) + 1.0) / (4.0 * PI) * legendre(l,Angle);
	    Ylev[lev-1][i * Ncurrent + j] += (2.0 * double(l) + 1.0) * legendre(l,Angle);
	}
      }
    }
  }
}
 
void UpPass(Complex*** zF, const double* dY, int currentlevel) {
/* interpolate zFchild for the nodes necessary for the parent
   (in general different expansion order, thus different nodes on the sphere)
   and shift the interpolated zFchild to the midpoint of the parent

   Lets assume we alread have the interpolation matrix Y
   However, we still need to multiply Fold with the quadrature weights over
   the sphere. Also it is assumed that it has been checked if an UpwardPass
   is even necessary

   Loop over all clusters in the current level
      get the Fold (local expansion) over all nodes on the sphere
      interpolate Fold to Fnew
      scale with exp(1.0i * k * (zparent - zchild) * s)

   Variables:
      zFchild: Complex matrix: 
               dimension: number of clusters times
	                  number of nodes on the sphere
               contains the local expansion inside each cluster as function
	       of the quadrature nodes
      zFparent: twin to zFchild on the new level
      dY:    double matrix containing the transformation matrix.
             Use column wise storage
      currentlevel: number of the current level
      
   GLOBAL: ClusterLev* clulevarry: cluster tree containing
                                   the clusters at all levels

   Notes: Since dY contains real numbers, Fold is split in real and imaginary
          part and the multiplication is done seperately
	  Also, zFChild is not needed anymore after this procedure, thus it
	  can be manipulated here
*/
  int i,j,child,nclus;  // number of clusters in the level
  int msphere,nsphere; // number of quad nodes for parent and child
  double z1z2[3]; // difference of cluster midpoints
  double ddummy;
  Complex alpha;

  Complex zexpfact, zdummy; // e^{ik (z2-z1)^T s}
  int parent;
  double* weisphe;
  
  nclus = clulevarry[currentlevel].nClustOLv;
  msphere = clulevarry[currentlevel-1].nPoinSpheLv;
  nsphere = clulevarry[currentlevel].nPoinSpheLv;
  weisphe = clulevarry[currentlevel].weisphe;
  double ddummymat[2*nsphere]; // real + imag part for Fold
  double ddummymat2[2*msphere]; // read + imag part for Fnew

  // set father entries to zero
  // for the far future, think about allocating the space
  // for zF in one go, and use memset
  
  //if( clulevarry[currentlevel-1].nClustOLv == 1 || clulevarry[currentlevel-1].isNearClust)
  //  return;
  
  for (i = 0; i < clulevarry[currentlevel-1].nClustOLv; i++)
    for ( j = 0; j < msphere; j++ ) 
      zF[currentlevel-1][i][j].set(0.0,0.0);
    
  for (child = 0; child < nclus; child++) {
    parent = clulevarry[currentlevel].ClustArLv[child].nuFather;
    for( j = 0; j < NDIM; j++ ){
      z1z2[j] = (clulevarry[currentlevel-1].ClustArLv[parent].CoorCent[j] - 
		  clulevarry[currentlevel].ClustArLv[child].CoorCent[j]);
    
    }
    
    if(msphere == nsphere ) { // no interpolation necessary just shift the
                           // cluster midpoint
      for (i = 0; i < msphere; i++) {
	ddummy = 0.0;
	for (j = 0; j < NDIM; j++) {
	  ddummy += z1z2[j] * clulevarry[currentlevel-1].uvcsphe[i][j];
	}
	zexpfact.set(cos( waveNumbers_ *  ddummy),sin( waveNumbers_*ddummy ));
	zF[currentlevel-1][parent][i] += zexpfact * zF[currentlevel][child][i];
      }
    }
    else {
  // now we have different expansion length in parent and child
  
  // Do the interpolation
  // 1) multiply with the quadrature weights
      for( i = 0; i < nsphere; i++) {
	zF[currentlevel][child][i] *= weisphe[i];
      }
      
      // 2) multiplay with the interpolation matrix
      ddummy = 0.0;
#ifdef USE_LAPACK
      // since yMat is double, split zF into real and imaginary
      for ( i = 0; i < nsphere; i++ ) {
	ddummymat[2*i] = zF[currentlevel][child][i].re();
	ddummymat[2*i+1] = zF[currentlevel][child][i].im();
      }
      //      for (i = 0; i < nsphere; i++ )
      //	ddummymat[nsphere + i] = zF[currentlevel][child][i].im();
      
      //cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,msphere,2,nsphere,1.0,dY,msphere,ddummymat,nsphere,0.0,ddummymat2,msphere);
      cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, msphere, 2, nsphere, 1.0, dY, nsphere, ddummymat, 2, 0.0, ddummymat2, 2);
      /*
      for (i = 0; i < msphere; i++) {
	cout << ddummymat2[2*i] << " " << ddummymat2[2*i+1] << "\n";
      }
      exit(0);
      */
      for (i = 0; i < msphere; i++) {
	ddummy = 0.0;
	for (j = 0; j < NDIM; j++) {
	  ddummy += z1z2[j] * clulevarry[currentlevel-1].uvcsphe[i][j];
	}
	zexpfact.set(cos( waveNumbers_ *  ddummy),sin( waveNumbers_*ddummy ));
	//zdummy.set(ddummymat2[i],ddummymat2[i+msphere]);
	zdummy.set(ddummymat2[2*i], ddummymat2[2*i+1]);
	zexpfact.mul_c(zdummy);
	zF[currentlevel-1][parent][i] += zexpfact;
      }
#else
      for (i = 0; i < msphere; i++) {
	ddummy = 0.0;
	for(j = 0; j < NDIM; j++) {
	  ddummy += z1z2[j] * clulevarry[currentlevel-1].uvcsphe[i][j];
	}
	zexpfact.set(cos( waveNumbers_ *  ddummy),sin( waveNumbers_*ddummy ));
	for( j = 0; j < nsphere; j++)
	  zF[currentlevel-1][parent][i] += zexpfact * ( zF[currentlevel][child][j] * dY[i*nsphere + j] );
      }
#endif
    } // else msphere != nsphere
  } // loop over all clusters
}

void UpPasslocal(Complex* zFchild,Complex** zFparent, const double* dY, int currentlevel, int iclus) {
  int i,j;  // number of clusters in the level
  int msphere,nsphere; // number of quad nodes for parent and child
  double z1z2[3]; // difference of cluster midpoints
  double ddummy;
  Complex alpha;
 // read + imag part for Fnew
  Complex zexpfact, zdummy; // e^{ik (z2-z1)^T s}
  int parent;
  
  msphere = clulevarry[currentlevel-1].nPoinSpheLv;
  nsphere = clulevarry[currentlevel].nPoinSpheLv;

  double ddummymat[2*nsphere]; // real + imag part for Fold
  double ddummymat2[2*msphere]; 

  
  parent = clulevarry[currentlevel].ClustArLv[iclus].nuFather;
  if(msphere == nsphere) { // no interpolation necessary
    for( j = 0; j < NDIM; j++ ){
      z1z2[i] = clulevarry[currentlevel-1].ClustArLv[parent].CoorCent[j] - 
	clulevarry[currentlevel].ClustArLv[iclus].CoorCent[j];
    }
    for (i = 0; i < msphere; i++) {
      for (j = 0; j < NDIM; j++) {
	ddummy += z1z2[j] * clulevarry[currentlevel-1].uvcsphe[i][j];
      }
      zexpfact.set(cos( waveNumbers_ *  ddummy),sin( waveNumbers_*ddummy ));
      zFparent[parent][i] += zexpfact * zFchild[i];
    }
    return;
  }
  // change in expansion length, do an interpolationx
  for( i = 0; i < nsphere; i++) {
    zFchild[i] *= weisphe[i];
    ddummy = 0.0;
  }

#ifdef USE_LAPACK
  // split real and imag part of  zFchild in two columns
  // then multiply both with the real interpolation matrix dY
  for ( i = 0; i < nsphere; i++ ) {
    ddummymat[2*i] = zFchild[i].re();
    ddummymat[2*i+1] = zFchild[i].im();
  }
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, msphere, 2, nsphere, 1.0, dY, nsphere, ddummymat, 2, 0.0, ddummymat2, 2);
  for (i = 0; i < msphere; i++) {
    for (j = 0; j < NDIM; j++) {
      ddummy += z1z2[j] * clulevarry[currentlevel-1].uvcsphe[i][j];
    }
    zexpfact.set(cos( waveNumbers_ *  ddummy),sin( waveNumbers_*ddummy ));
    zdummy.set(ddummymat2[2*i],ddummymat2[2*i+1]);
    zexpfact.mul_c(zdummy);
    zFparent[parent][i] += zexpfact;
  }
#else
  for (i = 0; i < msphere; i++) {
    // dummy = (z1 - z2)^T s
    ddummy = 0.0;
    for(j = 0; j < NDIM; j++) {
      ddummy += z1z2[j] * clulevarry[currentlevel-1].uvcsphe[i][j];
    }
    zexpfact.set(cos( waveNumbers_ *  ddummy),sin( waveNumbers_*ddummy ));
    for( j = 0; j < nsphere; j++)
      zFparent[parent][i] += zexpfact * ( zFchild[j] * dY[i*nsphere + j] ); 
  }
#endif
}


void DownPass(Complex** zGparent,Complex** zGchild, const double* dY, int currentlevel) {
/* shift the midpoint of the cluster farfield expansion G and filters the 
   nodes of the sphere

   Lets assume we alread have the interpolation matrix Y from the uppass
   
   However, we still need to multiply Fold with the quadrature weights over
   the sphere. Also it is assumed that it has been checked if an UpwardPass
   is even necessary

   Loop over all clusters in the current level
      get the Fold (local expansion) over all nodes on the sphere
      interpolate Fold to Fnew
      scale with exp(1.0i * k * (zparent - zchild) * s)

   Variables:
      zFchild: Complex matrix: 
               dimension: number of clusters times
	                  number of nodes on the sphere
               contains the local expansion inside each cluster as function
	       of the quadrature nodes
      zFparent: twin to zFchild on the new level
      zY:    complex matrix containing the transformation matrix.
             Use column wise storage
      currentlevel: number of the current level
   GLOBAL: ClusterLev* clulevarry: cluster tree containing
                                   the clusters at all levels
*/
  int nclus;  // number of clusters in the level
  int msphere,nsphere; // number of quad nodes for parent and child
  double z1z2[3]; // difference of cluster midpoints
  double ddummy;
  Complex alpha;
 
  int nChild,childcluster;
  
  nclus = clulevarry[currentlevel].nClustOLv;
  msphere = clulevarry[currentlevel].nPoinSpheLv;
  nsphere = clulevarry[currentlevel+1].nPoinSpheLv;

  double ddummymat[2*msphere]; // real + imag part for Fold
  double ddummymat2[2*nsphere]; // read + imag part for Fnew
  Complex zexpfact; // e^{ik (z2-z1)^T s}
  Complex zdummy[msphere];
  Complex ztest;

  if( msphere != nsphere) { 
    // interpolation necessary 
    for( int n = 0;  n < nclus; n++) {
      for( int s = 0; s < msphere; s++)
	zGparent[n][s].mul_r( clulevarry[currentlevel].weisphe[s] );
    }
  }
  
  for (int n = 0;  n < nclus; n++) { // parent cluster
    nChild = clulevarry[currentlevel].ClustArLv[n].n_Son;
    for (int ichild = 0; ichild < nChild; ichild++) {
      childcluster = clulevarry[currentlevel].ClustArLv[n].nuSon[ichild];
      // difference between cluster midpoints
      for(int j = 0; j < NDIM; j++ ){
	z1z2[j] = clulevarry[currentlevel+1].ClustArLv[childcluster].CoorCent[j] -  clulevarry[currentlevel].ClustArLv[n].CoorCent[j];
      }
      // shift between different cluster midpoints (+ demodulate)
      for(int i = 0; i < msphere; i++ ) {
	ddummy = 0.0;
	for( int j = 0; j < NDIM; j++) 
	  ddummy += clulevarry[currentlevel].uvcsphe[i][j] * z1z2[j];
	zexpfact.set( cos(waveNumbers_ * ddummy), sin(waveNumbers_ * ddummy) );
	zdummy[i] = zexpfact * zGparent[n][i];
#ifdef USE_LAPACK
	// split in real and imag part, the interpolation is for doubles
	if( msphere != nsphere   ) {
	  ddummymat[2*i] = zdummy[i].re();
	  ddummymat[2*i + 1] = zdummy[i].im();
	}
#endif
      }
      // interpolate/filter from m_points to n_points
      // there is some weird error when using lapack, which i have not
      // figured out yet
#ifdef USE_LAPACK
      if( msphere != nsphere  ) {
	cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,nsphere,2,msphere,1.0,dY,nsphere,ddummymat,2,0.0,ddummymat2,2);

	for( int i = 0; i < nsphere; i++ ) {
	  ztest.set(ddummymat2[2*i],ddummymat2[2*i + 1]);
	  zGchild[childcluster][i] += ztest;
	}
      }  
      else
	for( int i = 0; i < nsphere; i++ )
	  zGchild[childcluster][i] += zdummy[i];
      
#else
      if( msphere != nsphere ) {  
	for (int i = 0; i < nsphere; i++)
	  for( int j = 0; j < msphere; j++)
	    zGchild[childcluster][i] += zdummy[j] * dY[j * nsphere + i];
      }
      else
	for (int i = 0; i < nsphere; i++)
	  zGchild[childcluster][i] += zdummy[i];
#endif
    } // loop over children
  } // loop over parent cluster
}


#if 0
void Cluster2ClusterEval(Complex** zF, Complex** zG) {
  /* calculates the cluster 2 cluster contributions in the
     postprocessing
     zF in: local expansions
     zG out: far field cluster values

     Global
     clulevarry: BE cluster tree
     ipcluarry[ipc] evalcluster array, ipc in [0,ninpclus_)
     ninpclus_ number of evalclusters

     No explicit calulcation of the zMmat is necessary
     
  */
  L = clulevarry[0].nExpaTermLv;
  nsphere = clulevarry[0].nPoinSpheLv;
  double jn[L+1],yn[L+1];
  Pl = new double*[nsphere]; // we could do the allocation before the loop
                               // with the max number of sphere points
  if( Pl == NULL) {
    cerr << "Cannot allocate the Legendrepolynomials\n";
    exit(-1);
  }
  for (int i = 0; i < nsphere; i++) {
    Pl[i] = new double[L+1];
    if( Pl[i] == NULL) {
      cerr << "Cannot allocate the Legendrepolynomials\n";
      exit(-1);
    }
  }
  // loop over all clusters in level 
  for ( int n = 0; n < clulevarry[0].nClustOLv; n++ ) {
    // cluster center
    for( int nn = 0; nn < NDIM; nn++) 
      z2[nn] = clulevarry[level].ClustArLv[n].CoorCent[nn];
    for ( int i = 0; i < ninpclus_; i++) {
      r = 0.0;
      for(int  nn = 0; nn < NDIM; nn++) {
	z1z2[nn] = (ipcluarry[i].CoorCent[nn] - z2[nn]);
	r += z1z2[nn] * z1z2[nn];
	// chen switches the sign in his definition, adapt accordingly
	// for debug purposes negativ sign
	//z1z2[nn] = -z1z2[nn];
      }
      r = sqrt(r);
#ifdef USE_GSL
      gsl_sf_bessel_jl_array(L,r * waveNumbers_, jn);
      gsl_sf_bessel_yl_array(L,r * waveNumbers_, yn);
#else
      jn = sph_bessel(L, r * waveNumbers_);
      yn = sph_neumann(L, r * waveNumbers_);
#endif
      // P_l( (z1 - z2) * s )
      for(int s = 0; s < nsphere; s++) {
	v = 0.0;
	for( nn = 0; nn < NDIM; nn++ )
	  v += z1z2[nn] * clulevarry[0].uvcsphe[s][nn];
	v = v / r;
#ifdef USE_GSL
	gsl_sf_legendre_Pl_array(L, v, Pl[s]);
#else
	for(int l = 0; l < L + 1; l++) 
	  Pl[s][l] = legendre(L,v);
#endif
      }
      
      for( s = 0; s < nsphere; s++ )
	zM[s].set(0.0,0.0);
      for (l = 0; l < L; l++) {
	switch( l % 4 ) {
	case 0:
	  zfact.set(1.0,0.0);
	  break;
	case 1:
	  zfact.set(0.0,1.0);
	  break;
	case 2:
	  zfact.set(-1.0,0.0);
	  break;
	case 3:
	  zfact.set(0.0,-1.0);
	  break;
	}
	hn.set(jn[l],yn[l]);
	zfact *= (Complex)(2 * l + 1) * hn;
	for( int s = 0; s < nsphere; s++ ) {
	  zM[s] += zfact * Pl[s][l];

	}
      }
      for( int s = 0; s < nsphere; s++ ) {
	zG[i][s] += zM[s] * zF[n][s];
      }
    } // loop for i 
  } // loop n
  for (s = 0; s < nsphere; s++) {
    delete[] Pl[s];
    Pl[s] = NULL;
  }
  delete[] Pl;
}
#endif

void Cluster2Clustermat(int maxlev, bool allocateFMM) {
  /* computers the cluster to cluster interaction matrix 
     for all levels
     zMmat[level][clusterinlevel][interactionlist*nodesphere]

     The way zMmat is constructed is using
     loop over all levels
       loop cluster in level
           loop interaction list for each cluster
	      loop nodes sphere

     here in general, maxlev should be numClusterLevels_
     global variables
     clulevarry  clustertree

     Global: zMmat;
   */

  double r,v;
  double rwfact_ = 1.8; 
  // we currently use one order more then chen, check if this yields
  // much difference
  double** Pl; //[nsphere* (L+1)];
  int i,j,l,level,L,nsphere;
  int Clustj,C_i;
  Complex zfact,hn;
  double z2[3];
  double z1z2[3];
  int s,nn;
  struct rusage usage;
  getrusage(RUSAGE_SELF, &usage);
  int explength;
  // currently wasting much RAM
  cout << "zMat alloc begin: " << usage.ru_maxrss/1000 << "Mb\n";
  int dentries = 0;
  if( allocateFMM ) {
    zMmat = new Complex**[maxlev];
    if(zMmat == NULL) {
      cerr << "Sorry could not allocate zMmat\n";
      exit(-1);
    }
    for( level = 0; level < maxlev; level++ ) {
      // for the root level NumsFanClus point to NumsFarClus
      // thus far field clusters are used
      nsphere = clulevarry[level].nPoinSpheLv;
      zMmat[level] = new Complex*[clulevarry[level].nClustOLv];
      if( zMmat[level] == NULL) {
	cerr << "Sorry, could not allocate zMmat\n";
	exit(-1);
      }
      for ( Clustj = 0; Clustj < clulevarry[level].nClustOLv; Clustj++ ) {
	// clusters at current level
	// interaction list at all other levels
	zMmat[level][Clustj] = new Complex[ clulevarry[level].ClustArLv[Clustj].NumFanClus * nsphere];
	dentries += clulevarry[level].ClustArLv[Clustj].NumFanClus * nsphere;
	if( zMmat[level][Clustj] == NULL ) {
	      cerr << "Sorry, could not allocate zMmat[i][n]\n";
	      exit(-1);
	}
      }
    }  //level loop

  } // ifallocate loop
  // part that has to be done individually for each wavenumber
  cout << "Matrix D:\n";
  cout << "Complex: " << dentries << "\n";

  getrusage(RUSAGE_SELF, &usage);
  cout << "zMat: " << usage.ru_maxrss/1000 << "Mb\n";
  for( level = 0; level < maxlev; level++ ) {
    L = clulevarry[level].nExpaTermLv;
    nsphere = clulevarry[level].nPoinSpheLv;
    double jn[L+1],yn[L+1];
    Pl = new double*[nsphere]; // we could do the allocation before the loop
                               // with the max number of sphere points
    if( Pl == NULL) {
      cerr << "Cannot allocate the Legendrepolynomials\n";
      exit(-1);
    }
    for (i = 0; i < nsphere; i++) {
      Pl[i] = new double[L+1];
      if( Pl[i] == NULL) {
	cerr << "Cannot allocate the Legendrepolynomials\n";
	exit(-1);
      }
    }
    // loop over all clusters in level 
    for ( Clustj = 0; Clustj < clulevarry[level].nClustOLv; Clustj++ ) {
      // cluster center
      for( nn = 0; nn < NDIM; nn++) 
	z2[nn] = clulevarry[level].ClustArLv[Clustj].CoorCent[nn];
      for ( j = 0; j < clulevarry[level].ClustArLv[Clustj].NumFanClus; j++) {
	// interaction list for cluster_n

	
	C_i = clulevarry[level].ClustArLv[Clustj].NumsFanClus[j];

	if( adapt_fmmlength_ ) {
	  double d = clulevarry[level].ClustArLv[C_i].RadiClus + clulevarry[level].ClustArLv[Clustj].RadiClus;
	  double rw = d*waveNumbers_ + rwfact_*log10(d*waveNumbers_ + PI);
	  explength = (int)(rw);
	  if(explength - (double)rw >= 0.5) explength++;
	  if(explength < minExpansionTermsFMM_) {
	    explength = minExpansionTermsFMM_;
	  }
	}
	else
	  explength = L;
       

	
	r = 0.0;
	for( nn = 0; nn < NDIM; nn++) {
	  // interaction cluster - current cluster
	  z1z2[nn] = (clulevarry[level].ClustArLv[C_i].CoorCent[nn] - z2[nn]);
	  r += z1z2[nn] * z1z2[nn];
	  // chen switches the sign in his definition, adapt accordingly
	  // for debug purposes negativ sign
	  //z1z2[nn] = -z1z2[nn];
	}
	r = sqrt(r);
#ifdef USE_GSL
	gsl_sf_bessel_jl_array(explength,r * waveNumbers_, jn);
	gsl_sf_bessel_yl_array(explength,r * waveNumbers_, yn);
#else
	for ( int l = 0; l  < explength + 1; l++) {
	  jn[l] = sph_bessel(l, r * waveNumbers_);
	  yn[l] = sph_neumann(l, r * waveNumbers_);
	}
#endif
	// P_l( (z1 - z2) * s )
	for(s = 0; s < nsphere; s++) {
	  v = 0.0;
	  for( nn = 0; nn < NDIM; nn++ )
	    v += z1z2[nn] * clulevarry[level].uvcsphe[s][nn];
	  v = v / r;
#ifdef USE_GSL
	  gsl_sf_legendre_Pl_array(L, v, Pl[s]);
#else
	  for ( int l = 0; l < L + 1; l++ ) 
	    Pl[s][l] = legendre(l, v);
#endif
	}
	
	for( s = 0; s < nsphere; s++ )
	  zMmat[level][Clustj][s + nsphere * j].set(0.0,0.0);
	for (l = 0; l < explength; l++) {
	  switch( l % 4 ) {
	  case 0:
	    zfact.set(1.0,0.0);
	    break;
	  case 1:
	    zfact.set(0.0,1.0);
	    break;
	  case 2:
	    zfact.set(-1.0,0.0);
	    break;
	  case 3:
	    zfact.set(0.0,-1.0);
	    break;
	  }
	  hn.set(jn[l],yn[l]);
	  zfact *= (Complex)(2 * l + 1) * hn;
	  for( s = 0; s < nsphere; s++ ) {
	    zMmat[level][Clustj][s + nsphere * j] += zfact * Pl[s][l];
	  }
	}
      } // loop interaction list
    } // loop for clusters at level
    for (i = 0; i < nsphere; i++) {
      delete[] Pl[i];
      Pl[i] = NULL;
    }
    delete[] Pl;
  } // loop level
}   
  
void cluster2clusterlv(Complex*** zF, Complex*** zG, int level) {
  // does the cluster to cluster interaction at a given level
  /*  
      zF: in: zF[level][cluster][nodes sphere], local expansion for each cluster
              we define that for each level, although the old zF is not needed
	      after the uppass, but the allocated space changes from level
	      to level
      zG: out: zG[level][cluster2][nodes on sphere], the transformed zF
      level:  in: current level
      
      Note: we use zG with the actual clusternumber Gamma_i of the
            interaction cluster

      Global: zMmat: in: zMmat[level][clusterinlevel][interactionlist * nodessphere] translation matrix    
  */
  int i,C_i,s,n;
  int nsphere = clulevarry[level].nPoinSpheLv;

  for (n = 0; n < clulevarry[level].nClustOLv; n++)
    for( s = 0; s < nsphere; s++)
      zG[level][n][s].set(0.0,0.0);
  
  for( n = 0; n < clulevarry[level].nClustOLv; n++) {
    
    // cluster2cluster
    for (i = 0; i < clulevarry[level].ClustArLv[n].NumFanClus; i++) {
      C_i = clulevarry[level].ClustArLv[n].NumsFanClus[i];
      for( s = 0; s < nsphere; s++)
	zG[level][C_i][s] += zMmat[level][n][s + nsphere * i] * zF[level][n][s];
    }
    
  }
}

void cluster2clusterVec() {
  /* does the cluster to cluster interaction for the right hand side
     contribution, compared to the regular version, we assume a sparse
     local expansion

     Global:
     zMmat: in: zMmat[level][clusterinlevel][interactionlist * nodessphere]
                 translation matrix
      dYmat: in: dYmat[n-1] interpolation matrix from level n to n-1
  */
  int nsphere = clulevarry[nlevtop_].nPoinSpheLv;
  int n, i, C_i, C_j, s, l;
  
  //Complex*** zG;

  /*
  zG = new Complex**[numClusterLevels_];
  if( zG == NULL ) {
    cerr << "Sorry could not allocate zG\n";
    exit(-1);
  }
  
  for (i = 0; i < numClusterLevels_; i++) {
    zG[i] = new Complex*[ clulevarry[i].nClustOLv ];
    if( zG[i] == NULL) {
      cerr << "Sorry, Could not allocate zG\n";
      exit(-1);
    }
    
    for( C_j = 0; i < clulevarry[nlevtop_].nClustOLv; C_j++) {
      zG[i][C_j] = new Complex[ clulevarry[i].nPoinSpheLv ];
      if( zG[i][C_j] == NULL ) {
	cerr << "Sorry could not allocate zG\n";
	exit(-1);
      }
    }
    
  }
  */
  /* ********************************************************
  **                     leaf level                        **
  ******************************************************** */
  for( n = 0; n < zFvec.nonzeroblocks; n++) {
    C_j = zFvec.clusterindx[n];
    for( i = 0; i < clulevarry[nlevtop_].ClustArLv[C_j].NumFanClus; i++) {
      C_i = clulevarry[l].ClustArLv[C_j].NumsFanClus[i];
      for( s = 0; s < nsphere; s++)
	zG[nlevtop_][C_i][s] += zMmat[nlevtop_][C_j][s + nsphere * C_i] * zFvec.zdata[n * nsphere + s];
    }
  }
  /* **********************************************************
  **                      Up and Down                        **
  ********************************************************** */
  // allocate zF
  
  // Complex*** zF;
  int iclus;
  /*
  zF = new Complex**[nlevtop_]; // leaf level not needed
  if( zF == NULL ) {
    cerr << "Sorry could not allocate zF\n";
    exit(-1);
  }
  for( l = 0; l < nlevtop_; l++) {
    zF[l] = new Complex*[clulevarry[l].nClustOLv];
    if( zF[l] == NULL ) {
      cerr << "Sorry could not allocate zF\n";
      exit(-1);
    }
    for (i = 0; i < clulevarry[l].nClustOLv; i++) {
      zF[l][i] = new Complex[ clulevarry[l].nPoinSpheLv ];
      if( zF[l][i] == NULL ) {
	cerr << "Sorry could not allocate zF\n";
	exit(-1);
      }
    }
  }
  */
  for( l = nlevtop_; l > 0; l-- ) {
    // upwardpass for F + FMM interaction
    //if( l == 1 ) {
    if ( clulevarry[l-1].nClustOLv == 1 ) //|| clulevarry[l-1].isNearClust )
	continue;
      //}
    if( l == nlevtop_ ) // leaf
      for (i = 0; i < zFvec.nonzeroblocks; i++) {
	iclus = zFvec.clusterindx[i];
	UpPasslocal(&zFvec.zdata[i * nsphere], zF[l-1], dYmat[l-1], l, iclus);
      }
    else {
      //if( l == 1 )
      //	if( clulevarry[l-1].isNearClust || clulevarry[l-1].nClustOLv == 1 )
      if( clulevarry[l-1].nClustOLv == 1 )
	  continue;
      UpPass( zF,  dYmat[l-1], l);
    }
    
    // zF[l] is not needed anymore
    /*
    for( C_j = 0; i < clulevarry[l].nClustOLv; C_j++) {
      delete[] zF[l][C_j];
    }
    delete[] zF[l];
    */
    // cluster2Cluster
    //nclus = clulevarry[l-1].nClustOLv;
    cluster2clusterlv(zF,  zG, l-1);
  }
  /*
  for( i = 0; i < clulevarry[0].nClustOLv; i++) 
    delete[] zF[0][i];
  delete[] zF[0];
  delete zF;
  zF = NULL;
  */
  // downpass
  for(l = 0; l < nlevtop_; l++) {
    // child , parent, Y, level
    //    if( l == 0 )
    if( clulevarry[i].nClustOLv == 1 ) // || clulevarry[i].isNearClust )
	continue;
    DownPass( zG[l], zG[l+1], dYmat[l], l);
    // delete zG[level]
    /*
    for (i = 0; i < clulevarry[l].nClustOLv; i++) 
      delete[] zG[l][i];
	  delete[] zG[l];
    */
  }
  // local expansion
  Expand2local(zG[nlevtop_], zrhs);
  /*
  for (i = 0; i < clulevarry[nlevtop_].nClustOLv; i++) 
    delete[] zG[nlevtop_][i];
  delete[] zG[nlevtop_];
  delete[] zG;
  zG = NULL;
  */
}
  
void Cluster2Cluster(Complex* zF,Complex* zG,int nsphere, int L, double* z1z2, double** uvcsphere) {
  /* computers the cluster to cluster interaction
     zF  in: complex vector length number of quad nodes on sphere, contains the local expansion
     zG  out: complex vector length number of quad nodes on sphere.cluster to cluster interaction result for each node on the sphere
     nsphere: in: number of nodes on the sphere
     L: in: expansion length
     z1z2: double[3]  difference z1 - z2
   */

  double r,v;
  double jn[L+1],yn[L+1];
  double** Pl; //[nsphere* (L+1)];
  int i,j,l,n;
  Complex zfact,hn;
  
  Pl = new double*[nsphere];
  if( Pl == NULL ) {
    cerr << "Cannot allocate the Legendrepolynomials\n";
    exit(-1);
  }
  for (i = 0; i < nsphere; i++) {
    Pl[i] = new double[L+1];
    if( Pl[i] == NULL ) {
      cerr << "Cannot allocate the Legendrepolynomials\n";
      exit(-1);
    }
  }

  r = 0.0;
  for( i = 0; i < NDIM; i++ ) 
    r += z1z2[i] * z1z2[i];
  r = sqrt(r);
#ifdef USE_GSL
  gsl_sf_bessel_jl_array(L,r * waveNumbers_, jn);
  gsl_sf_bessel_yl_array(L,r * waveNumbers_, yn);
#else
  for ( int l = 0; l < L + 1; l++) {
    jn[l] = sph_bessel(l, r * waveNumbers_);
    yn[l] = sph_neumann(l, r * waveNumbers_);
  }
#endif
  for(n = 0; n < nsphere; n++) {
    v = 0.0;
    for( j = 0; j < NDIM; j++ )
      v += z1z2[j] * uvcsphere[i][j];
    v = waveNumbers_ * v / r;
#ifdef USE_GSL
    gsl_sf_legendre_Pl_array(L, v, Pl[n]);
#else
    for( l = 0; l < L + 1; l++ )
      Pl[n][l] = legendre(L,v);
#endif
  }
  
  for (l = 0; l < L + 1; l++) {
    switch( l % 4 ) {
    case 0:
      zfact.set(1.0,0.0);
      break;
    case 1:
      zfact.set(0.0,1.0);
      break;
    case 2:
      zfact.set(-1.0,0.0);
      break;
    case 3:
      zfact.set(0.0,-1.0);
      break;
    }
    hn.set(jn[l],yn[l]);
    zfact *= (Complex)(2 * l + 1) * hn;
    for( i = 0; i < nsphere; i++ ) {
      zG[i] += zfact * Pl[i][l] * zF[i];
    }
  }
  for (i = 0; i < nsphere; i++) {
    delete[] Pl[i];
    Pl[i] = NULL;
  }
  delete[] Pl;
}   

void get_nearfield(int level, double* scalefact, Complex& zBta3, bool allocateFMM) {
  /* set up the Near field matrix for the leaf clusters

    
     level: in:  leaf level index
     scalefact: out: vector(numRowsOfCoefficientMatrix_)
                     factor used in the preconditioning for scaling an
                     individual row of system matrix and rhs 		 
     zBta3: in: complex, Burton Miller factor
     allocateFMM: in: bool, should the matrices be also allocated ?
     
     A major "problem" here is that depending on the boundary conditions
     the formulation needs to be adapted, thus we need to consider the integral
     over the boundary integral operators, the right hand side, and the
     u/2 and v/2 factors

     Note: In the original version, Chen used the system
          -u/2 + Hu - Gv = u_inc
	  
     Global:
     zrhs: in/out Complex array of length numRowsOfCoefficientMatrix_, contains the entries for
           the right hand side
     ibval: in: integer array of length numRowsOfCoefficientMatrix_, type of the bc for the element
                ibval = 0, velocity bc
		ibval = 1, press bc
		ivbal = 2,5  velobc + admi
     zbval0: in: complex array of length numRowsOfCoefficientMatrix_, value of the boundary condition
     Sourpoi3: out: center of the element Gamma_i
     Norvci3: out:  normal vector for element Gamma_i

     zNear: out: Complex array containing the entries of the
                 near field matrix, the entries are stored in
		 sparse rowwise format
     zrhs: out: nearfield contributions of boundary conditions to the right
                hand side
     
     
  */

  //  Matrix<double> crdelj(NNODPE, NDIM);
  double* rownorm2;
  rownorm2 = new double[numRowsOfCoefficientMatrix_];
  //  double rownorm2[numRowsOfCoefficientMatrix_];
  int* nentryinrow;
  nentryinrow = new int[numRowsOfCoefficientMatrix_];
  
  
  int totallength,blockstart,rowstart;
  int Ibvi03, Ibvj03, Ifcrh3;
  // for compatibility reasons we have to use a vector here instead of just
  // one complex value
  Vector<Complex> Zbvi03(4);
  Vector<Complex> Zbvj03(4);
  Complex Admia3; // admittance bc val
  bool Ifadmii, Ifadmij;
  int i, i1, n,  Gamma_i; // target stuff -> rows
  int j, j1, j2, nGamma_j, Gamma_j; // field stuff -> columns
  int Ci; // target cluster
  int bcval_i, bcval_j;  // orig def by chen as int
  Matrix<double> Crdej3(NNODPE, NDIM);
  Vector<double> center_j(3);
  Vector<Complex> zrsintel(6);
  int localstart,rowlength;
  for( i = 0; i < numRowsOfCoefficientMatrix_; i++) {
    nentryinrow[i] = 0;
    rownorm2[i] = 0.0;
  }
  // get the number of total nonzero entries in the near field matrix
  // loop over all cluster
  for ( n = 0; n < clulevarry[level].nClustOLv; n++) {
    nGamma_j = clulevarry[level].ClustArLv[n].NumOfEl; 
      // loop over all near field clusters
    for ( i = 0; i < clulevarry[level].ClustArLv[n].NumNeaClus; i++ ) {
      Ci = clulevarry[level].ClustArLv[n].NumsNeaClus[i]; // this would be the target clusters
      // loop over the targetcluster elements
      for (i1 = 0; i1 < clulevarry[level].ClustArLv[Ci].NumOfEl; i1++) {
	Gamma_i = clulevarry[level].ClustArLv[Ci].NumsOfEl[i1];
	nentryinrow[Gamma_i] += nGamma_j;
      }
    }
  }
  if( zNear.zdata != NULL )
    allocateFMM = false;
  
  totallength = 0;
  if( allocateFMM || zNear.startlist == NULL)  
    zNear.startlist = new int[numRowsOfCoefficientMatrix_+1];
  
  zNear.startlist[0] = 0;
  for (i = 0; i < numRowsOfCoefficientMatrix_; i++) {
    // now really ordered with respect to Gamma_i
    zNear.startlist[i+1] = zNear.startlist[i] + nentryinrow[i];  
    totallength += nentryinrow[i];
  }
  
  // now we are ready to initialize zNear
  if( allocateFMM || zNear.zdata == NULL) {
    if( zNear.zdata != NULL ) {
      cout << "there is a problem with zNear in BAMLFMM\n";
      exit(-1);
    }
    zNear.zdata = new Complex[totallength];
    if( zNear.zdata == NULL ) {
      cerr << "Sorry could not allocate zNear\n";
      exit(-1);
    }
    if( zNear.indxlist != NULL ) {
      cout << "there is a problem with zNear in BAMLFMM\n";
      exit(-1);
    }
    zNear.indxlist =new int[totallength];
    if( zNear.indxlist == NULL ) {
      cerr << "Sorry could not allocate zNear\n";
      exit(-1);
    }
    zNear.rowwise = true;

  }

  cout << "Nearfield matrix\n";
  cout << "Integer: " << numRowsOfCoefficientMatrix_+1 + totallength << "\n";
  cout << "Complex: " << totallength << "\n";
  
  for( i = 0; i < numRowsOfCoefficientMatrix_; i++)
    nentryinrow[i] = 0;


  // do the parts independent of Gamma_j
  for( i = 0; i < numRowsOfCoefficientMatrix_; i++ ) {
    // this would be the Gamma_i
    for( j = 0; j < NDIM; j++ ) {
      // global variables containing midpoint and normal vector
      // for Gamma_i
      Sourpoi3[j] = centel[i][j];
      Norvci3[j] = elenor[i][j];
    }
    switch(ibval[i]) {
    case 0:
      Ibvi03 = 0;
      Ifadmii = false;
      break;
    case 1:
      Ibvi03 = 1;
      Ifadmii = false;
      break;
    case 2:
      Ibvi03 = 0;
      Ifadmii = true;
      break;
    default:
      cerr << "Sorry wrong boundary condition\n";
      exit(-1);
    }
    for (j = 0; j < NNODPE; j++) 
      Zbvi03[j] = zbval0[i][j];

    if( Zbvi03[0].norm() > EPSY )
      bcval_i = 1;
    else
      bcval_i = 0;
    if(numIncidentPlaneWaves_ + numPointSources_ > 0) {
      // uses global Sourpoi3, right hand side, uinc
      zrhs[i] += NC_IncidentWaveRHS(NCout);
	//BAinwpos(NCout, 0);
    }
    if( Ifadmii )
      Admia3 = zbval1[i];
    // not really necessary for surface elements and constant
    // collocation, but sooner or later necessary for other cases
    // Dofpos_i = jelist[Gamma_i][0];
    
    // right hand side parts, don't forget chens uses -1/2u + ...
    if( Ibvi03 ) { // PRES 
      if( bcval_i ) 
	zrhs[i] += Zbvi03[0] * 0.5;
    }
    else { // VELO condition
      if( bcval_i )
	zrhs[i] += Zbvi03[0] * 0.5 * zBta3;
    }
  }

  // get the data for the near field matrix
  // loop over all clusters
  for ( n = 0; n < clulevarry[level].nClustOLv; n++) {
    // loop over all near field clusters
    for ( i = 0; i < clulevarry[level].ClustArLv[n].NumNeaClus; i++ ) {
      Ci = clulevarry[level].ClustArLv[n].NumsNeaClus[i];
      // near field element loop
      for (i1 = 0; i1 < clulevarry[level].ClustArLv[Ci].NumOfEl; i1++) {
	Gamma_i = clulevarry[level].ClustArLv[Ci].NumsOfEl[i1];

	for( j = 0; j < NDIM; j++ ) {
	  // global variables containing midpoint and normal vector
	  // for Gamma_i
	  Sourpoi3[j] = centel[Gamma_i][j];
	  Norvci3[j] = elenor[Gamma_i][j];
	}
	// nentryinrow will be updated on the fly for each entry
	rowstart = zNear.startlist[Gamma_i] + nentryinrow[Gamma_i];
	// get the boundary conditions
	
        // leaf cluster element loop
	for (j = 0; j < clulevarry[level].ClustArLv[n].NumOfEl; j++) {
	  // Number of field element
	  Gamma_j = clulevarry[level].ClustArLv[n].NumsOfEl[j];
	  // vertices of field element
	  for ( j1 = 0; j1 < listNumberNodesPerElement[Gamma_j]; j1++) {
	    for ( j2 = 0; j2 < NDIM; j2++ ) {
	      Crdej3(j1,j2) = nodesCoordinates[ elementsConnectivity[Gamma_j][j1] ][j2];
	    }
	  }
	  
	  // get the boundary conditions and the normalvec (kappa) and
	  // the distance to the collocnode
	  /*BLeldat(Gamma_j, 0, Inoj3, Thiej3, Ibvj_03,
		  Ifadmj3, Jdofaddr3, Centej3, Norvcj3, Admj3,
		  Zbvj03, Crdej3);
	  */
	  switch(ibval[Gamma_j]) {
	  case 0:
	    Ibvj03 = 0;
	    Ifadmij = false;
	    break;
	  case 1:
	    Ibvj03 = 1;
	    Ifadmij = false;
	    break;
	  case 2:
	    Ibvj03 = 0;
	    Ifadmij = true;
	    break;
	  default:
	    cerr << "Sorry wrong boundary condition\n";
	    exit(-1);
	  }

	  for( j2 = 0; j2 < NNODPE; j2++) 
	    Zbvj03[j2] = zbval0[Gamma_j][j2];
	  if( Zbvj03[0].norm() > EPSY )
	    bcval_j = 1;
	  else
	    bcval_j = 0;
	  if( Ifadmij )
	    Admia3 = zbval1[Gamma_j];
	  
	  /* Calculate the integrals over Gamma_j
	  do the sgl integration
	  zrsintel[0..5] has the results of the sgl integration
	  [0] -> G
	  [1] -> H
	  [2] -> H'
	  [3] -> E
	  [4] -> rhs
	  it may seem weird, but chen messes up the signs, and he
	  uses r = | y - x | */
	  if ( Gamma_i == Gamma_j ) 
	    NC_SingularIntegration(NCout, zrsintel, listNumberNodesPerElement[Gamma_i], bcval_j, Zbvj03, Crdej3);
	  else {
	    for (int n1 = 0; n1 < NDIM; n1++)
	      center_j[n1] = centel[Gamma_j][n1];
	    NC_RegularIntegration(NCout, zrsintel, Gamma_i, Gamma_j, listNumberNodesPerElement[Gamma_j],bcval_j, Zbvj03, center_j, Crdej3);
	  }
	  if( Ibvj03 ) // PRES CONDITION
	    zNear.zdata[rowstart + j] = -zrsintel[0] - zBta3 * zrsintel[2];
	  else {
	    zNear.zdata[rowstart + j] = zrsintel[1] + zBta3 * zrsintel[3];
	    if( Ifadmij ) {
	      zNear.zdata[rowstart + j] += zrsintel[0] * Admia3 + zrsintel[2] * Admia3 * zBta3;
	    }
	  }
	  
	  if( bcval_j ) // there are some contributions to the rhs
	    zrhs[Gamma_i] += zrsintel[4];
	  //cout << "Warning sign for the rhs may be wrong!\n";
	  
	  // Free terms u/2 and beta v/2 and (1 + beta alpha)/2
	  // depending on the boundary condition 
	  if( Gamma_i == Gamma_j ) {
	    if( Ibvj03 ) // PRES CONDITION
	      zNear.zdata[rowstart + j] -= zBta3 * 0.5;
	    else // VELO CONDITION
	      if( Ifadmij ) {
		zNear.zdata[rowstart + j] += (Admia3 * zBta3 - 1.0) * 0.5;
	      }
	      else
		zNear.zdata[rowstart + j] -= 0.5;
	  }
	  // nobody said that the entries per row need to be ordered
	  zNear.indxlist[rowstart + j] = Gamma_j;
	  rownorm2[Gamma_i] += zNear.zdata[rowstart + j].qnorm();
	} // loop Gamma_j
	nentryinrow[Gamma_i] += clulevarry[level].ClustArLv[n].NumOfEl;
      } // loop Gamma_i
    } // loop Near cluster
  } //loop cluster
  
  for (i = 0; i < numRowsOfCoefficientMatrix_; i++) {
    scalefact[i] = sqrt( (double)nentryinrow[i] / rownorm2[i]);
  }
  int maxlength = 0;
  for (i = 0; i < numRowsOfCoefficientMatrix_; i++) {
    if( zNear.startlist[i+1] - zNear.startlist[i] > maxlength)
      maxlength = zNear.startlist[i+1] - zNear.startlist[i];
  }
  int indxlist[maxlength];
  Complex dummy[maxlength];
  int idummy[maxlength];
  // sort the rows of zNear
  
  for (i = 0; i < numRowsOfCoefficientMatrix_; i++) {
    localstart = zNear.startlist[i];
    rowlength = zNear.startlist[i+1] - zNear.startlist[i];
    // maybe thing about 
    sortArr( &zNear.indxlist[localstart], rowlength, indxlist );
    for (j = 0; j < rowlength; j++) {
      dummy[j] = zNear.zdata[ localstart + indxlist[j] ];
      //idummy[j] = zNear.indxlist[ localstart + indxlist[j] ];
    }
    for (j = 0; j < rowlength; j++) {
      zNear.zdata[ localstart + j ] = dummy[j];
      //zNear.indxlist[ localstart + j ] = idummy[ j ];
    }

    
  }
  delete [] rownorm2;
  delete [] nentryinrow;
  
}

void setup_preconditioning(double* scalefact) {
/* sets up the lower and upper part of the incomplete LU decomposition of the
   near field matrix. L is stored in rowwise format, U in columnwise format
   The near field matrix is scaled with the Number of elements per row/ norm row
   If the entries are above a certain tolerance the non zero pattern will be
   set to 1 at this position

   scalefact: in, scaling factor for each row in the near field matrix
    set up the nonzero structure for L and U
        L contains the lower non-zero part of zcoefl
   U contains the upper non-zero part of zcoefl
   copy the right parts of the scaled matrix to L and U
   L = 1                    U =    1 2 4 7
       2 3                           3 5 8
       4 5 6                           6 9
       7 8 9 10                          10

   may seem a bit weird because of the order the entries are generated
   but it makes sense when computing them, because there is an additional
   vector product, and in the back subsitutions at the end of the calculations

    getting from row wise to columnwise:
	       go through rowstart vector, get the first entries per line
	       if jcol({rowstart[i]}) == 0
	          add i to the position list
		  rowstart[i] = rowstart[i+1]
	       end
	       this means only go to the next entry of rowstart[i] == 0
	       if the row is finished enter -1 in the rowstart array
	       repeat procedure and compare with 1
	       repead procedure and compare with 2

    we use two loops, one for counting the number of entries in  L and U
    and a second to copy the data from the nearfield matrix

    Global
       zNear: inout, sparse nearfield matrix, gets scaled by scalefact
       zL: out, sparse matrix containing the L part of the incompl. LU decomp
       zU: out,

  */

  double threshfac; // threshold for incomplete LU 
                    // non-zero pattern

  int n,i,i0,i1,j,j0,clusterstart,k;
  int Gamma_i, Gamma_j,Ci;
  int rowcounter = 0;
  int counter = 0;
  int rowstart = 0;
  int rowentry = 0;
  int roworder[numRowsOfCoefficientMatrix_];
  int entriesU, entriesL; // number of nonzero entries in L and U
  zSparsetype zU0;
  int ncols[numRowsOfCoefficientMatrix_-1];
  int m,currentmL, currentmU;

  if(methodFMM_ == 1) { // SLFMBEM
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
  }
  else {
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

  // Go through the rows of N and count the number of entries
  // so we can allocate the sparese matrix representation
  //
  // an alternative would be to use vectors and push_back()
  entriesU = 0;
  entriesL = 0;
  // loop over all rows = Gamma_i
  for (i = 0; i < numRowsOfCoefficientMatrix_; i++) {
    zrhs[i].mul_r(scalefact[i]);
    // loop over entries in row_i
    for (j = zNear.startlist[i]; j < zNear.startlist[i+1]; j++) {
      zNear.zdata[j].mul_r(scalefact[i]);
      if( i == zNear.indxlist[j] || zNear.zdata[j].norm() > threshfac ) {
	// scale the near field matrix for preconditioning
	if( zNear.indxlist[j] > i ) {
	  entriesU++;
	}
	else {
	  entriesL++;
	}
      }
    }
  }
  
  // initialise U and L
  
  zL.zdata = new Complex[entriesL];
  if( zL.zdata == NULL ) {
    cerr << "Sorry, could not allocate zL\n";
    exit(-1);
  }
  zL.startlist = new int[numRowsOfCoefficientMatrix_+1];
  zL.indxlist = new int[entriesL];
  if( zL.indxlist == NULL ) {
    cerr << "Sorry, could not allocate zL\n";
    exit(-1);
  }
  zU.zdata = new Complex[entriesU];
  if( zU.zdata == NULL ) {
    cerr << "Sorry, could not allocate zU\n";
    exit(-1);
  }
  zU.indxlist = new int[entriesU];
  zU.startlist = new int[numRowsOfCoefficientMatrix_ + 1];
  if( zU.indxlist == NULL ) {
    cerr << "Sorry, could not allocate zU\n";
    exit(-1);
  }
  zU0.zdata = new Complex[entriesU];
  if( zU0.zdata == NULL ) {
    cerr << "Sorry, could not allocate zU\n";
    exit(-1);
  }

  zU0.startlist = new int[numRowsOfCoefficientMatrix_ + 1];
  zU0.indxlist = new int[entriesU];
  if( zU0.indxlist == NULL ) {
    cerr << "Sorry, could not allocate zU\n";
    exit(-1);
  }
  
  // for U the diagonal is not counted
  // we first save U as U0 in rowwise format
  entriesL = 0;
  entriesU = 0;
  zL.startlist[0] = 0;
  zU0.startlist[0] = 0;
  
  for( i = 0; i < numRowsOfCoefficientMatrix_; i++ ) { // rows of N
    for( j = zNear.startlist[i]; j < zNear.startlist[i+1]; j++) {
      // i-th row, j-indices should be already ordered
      j0 = zNear.indxlist[j];  // columnindex N_{i,j0}, \Gamma_j
      if( i == j0 || zNear.zdata[j].norm() > threshfac ) {
	if( j0 > i ) {
	  // U part rowwise
	  zU0.zdata[entriesU] = zNear.zdata[j];
	  zU0.indxlist[entriesU] = j0;
	  entriesU++;
	}
	else {
	  // L part rowwise
	  zL.zdata[entriesL] = zNear.zdata[j];
	  zL.indxlist[entriesL] = j0;
	  entriesL++;
	}
      }
    }
    zU0.startlist[i+1] = entriesU;
    zL.startlist[i+1] = entriesL;

  }
  zL.rowwise = true;

  // U is a bit trickier,
  // 'cause we have to switch from rowwise to columnswise
  int colpos[numRowsOfCoefficientMatrix_-1];
  for (i = 0; i < numRowsOfCoefficientMatrix_-1; i++) {
    // first nonzero column in every row, remember: upper triangle
    // diagonal is not included
    colpos[i] =  i+1;
    // number of columns in each row
    ncols[i] = zU0.startlist[i+1] - zU0.startlist[i];
  }

  int colindx_in_row[numRowsOfCoefficientMatrix_];
  for ( i = 0; i < numRowsOfCoefficientMatrix_; i++ ) {
    if( zU0.startlist[i+1] - zU0.startlist[i] > 0 )
      colindx_in_row[i] = zU0.startlist[i];
    else
      colindx_in_row[i] = -1;
  }
  counter = 0;

  zU.startlist[0] = 0;
  for (j = 1; j < numRowsOfCoefficientMatrix_; j++) { // columns of U, diagonal is not part of U
    zU.startlist[j] = counter;
    for ( i = 0; i < j; i++ ) { // go through rows of U0
      if( colindx_in_row[i] == -1 )
	continue;
      if( zU0.indxlist[ colindx_in_row[i] ] == j ) { // the j-th col in row i is nonzero
	zU.indxlist[counter] = i;
	zU.zdata[counter++] = zU0.zdata[ colindx_in_row[i] ];
	colindx_in_row[i]++;
	if( colindx_in_row[i] == zU0.startlist[i+1] )
	    colindx_in_row[i] = -1;
      }
    }
  }
  
  
  zU.startlist[ numRowsOfCoefficientMatrix_ ] = counter;
  zU.rowwise = false;
  //zU0.delete(); destruct should be called anyway
  /* for debugging
  cout << "indx = [";
  for( i = 0; i < entriesU; i++ )
    cout << zU.indxlist[i] << "\n";
  
  cout << "];\n";
  cout << "startindx = [";
  for( i = 0; i < numRowsOfCoefficientMatrix_; i++ )
    cout << zU.startlist[i] << "\n";
  cout << "];\n";
  cout << "zvals = [";
  for( i = 0; i < entriesU; i++) {
    cout << zU.zdata[i].re() << " " << zU.zdata[i].im() << "\n";
  }
  cout << "]\n";
  exit(-1);
  */
  // now that everything is copied, calculate the factors for the incomplete
  // LU. see e.g. Meister or Saad for the algorithm
  for( i = 0; i < numRowsOfCoefficientMatrix_; i++) {
    i0 = 0;
    for( k = i; k < numRowsOfCoefficientMatrix_; k++ ) {
      // get the nonzero pattern for the k-th row
      if( !IsNonZero_LU(k, i, i0, zL) ) // finds the right i0
	continue;
      currentmL = 0;
      currentmU = 0;
      for (m = 0; m < i; m++) {
	if( IsNonZero_LU(k,m,currentmL,zL) ) {
	  if( IsNonZero_LU(m,i,currentmU,zU) ) {
	    // L_{ki} = A_{ki} - \sum L_{km} U_{mi}
	    zL.zdata[ zL.startlist[k] + i0 ] -= zL.zdata[ zL.startlist[k] + currentmL ] * zU.zdata[ zU.startlist[i] + currentmU ];
	  }
	}
      }
    }
    // the U part
    i0 = 0;
    for ( k = i+1; k < numRowsOfCoefficientMatrix_; k++ ) {
      if( !IsNonZero_LU(i,k,i0,zU) )
	continue;
      currentmL = 0;
      currentmU = 0;
      for ( m = 0; m < i; m++ ) {
	if( IsNonZero_LU(i,m,currentmL,zL) ) {
	  if( IsNonZero_LU(m,k,currentmU,zU) ) {
	    // U_{ik} = 1/L_{ii} (A_[ik} - sum L_{im}U_{mk}
	    zU.zdata[ zU.startlist[k] + i0 ] -= zL.zdata[ zL.startlist[i] + currentmL ] * zU.zdata[ zU.startlist[k] + currentmU];
	  }
	}
      }
      zU.zdata[ zU.startlist[k] + i0 ] /= zL.zdata[ zL.startlist[i+1] - 1];
    }
  }
}
  

bool IsNonZero_LU(int i, int j, int& kpos, zSparsetype& A) {
  /* checks if A_{ij} != 0 and returns the correct index in the sparse format in
     kpos, however this assumes that the indxlist is ordered per row/column
     variables
     i: input integer: rowindx
     j: input integer: colidnx
     kpos: out integer: position in the sparese matrix
     A: input : sparse Matrix
  */
  int k;
  int startindx;
  int endindx;
  
  if( A.rowwise ) {
    // we have L
    startindx = A.startlist[i];
    endindx = A.startlist[i+1];
    k = 0;
    if ( startindx + k >= endindx )
      return false;
    do {
      if( A.indxlist[startindx + k] == j ) {
	kpos = k;
	return true;
      }
      k++;
    } while( startindx + k < endindx && A.indxlist[startindx + k] <= j );
    return false;
  }
  else {
    // we have U, colwise format
    startindx = A.startlist[j];
    endindx = A.startlist[j+1];
    k = 0;
    if( startindx + k > endindx - 1 )
      return false;
    do {
      if( A.indxlist[startindx + k] == i ) {
	kpos = k;
	return true;
      }
      k++;
    } while( (startindx + k < endindx) && A.indxlist[startindx + k] <= i );
    return false;
  }
}

bool IsNonZero_LU(int i, int j, zSparsetype& A) {
  /* checks if A_{ij} != 0 and returns the correct index in the sparse format in
     kpos, however this assumes that the indxlist is ordered per row/column
  */
  int k;
  int start;
  int end;
  
  if( A.rowwise ) {
    // we have L
    start = A.startlist[i];
    end = A.startlist[i+1];
    k = 0;
    if ( start + k >= end )
      return false;
    do {
      if( A.indxlist[start + k] == j ) {
	return true;
      }
      k++;
    } while( A.indxlist[start + k] <= j && start + k < end);
    return false;
  }
  else {
    // we have U
    start = A.startlist[j];
    end = A.startlist[j+1];
    k = 0;
    if( start + k > end - 1 )
      return false;
    do {
      if( A.indxlist[start + k] == i ) {
	return true;
      }
      k++;
    } while(A.indxlist[start + k] <= i && start + k < end);
    return false;
  }
}



void sortArr(int x[], int n, int indx[]) {

  int i;
  vector< pair<int, int> > vp;
  for (i = 0; i < n; i++) {
    vp.push_back( make_pair(x[i],i) );
  }

 std:sort(vp.begin(),vp.end());   // std:sort
  for( i = 0; i < vp.size(); i++ ) {
    indx[i] = vp[i].second;
    x[i] = vp[i].first;
  }
}



void Cluster2Local(Complex** zGmat, double* scalefact, Complex* zx, Complex& zBta3) {
  /* local cluster to element expansion on the leaf level,
  
     zGmat: in: matrix containing the the value of the MLFMM at each
            quadrature node of the sphere
     scalefact: scaling factor from the preconditioning
     zx: output FMM value at each collocnode
     zBta: in: Burton Miller factor
     Global: clulevarray: clustertree for all levels
             nlevtop_: leaf level
             centel: midpoint of each element = collocnode

     Note: this routine has to be called at every iteration step, thus it may
           save memory, but takes longer than setting up the matrix beforehand
  */
  int nsphere;
  int n,i,j,Gamma_i,s;
  double z1[3],x[3]; // coordinates of the cluster center, coord of elem center
 
  double nvec[3]; // normal vector of the element
  int Ibvi03, Ifadmii; // boundary conditions of the target element
  double sn, sx; // node on the sphere times normal vector to element
                       // node on the sphere times x
  Complex sfact,expfact;
  

  //dont forget the multiplication with zbta3
  

  nsphere = clulevarry[nlevtop_].nPoinSpheLv;
  

  for (n = 0; n < clulevarry[nlevtop_].nClustOLv; n++) {
    for( s = 0; s < nsphere; s++) {
      zGmat[n][s] *= clulevarry[nlevtop_].weisphe[s];
    }
    for( j = 0; j < NDIM; j++) {
      z1[j] = clulevarry[nlevtop_].ClustArLv[n].CoorCent[j];
    }
    for (i = 0; i < clulevarry[nlevtop_].ClustArLv[n].NumOfEl; i++) {
      Gamma_i = clulevarry[nlevtop_].ClustArLv[n].NumsOfEl[i];
      for( j = 0; j < NDIM; j++) {
	x[j] = centel[Gamma_i][j] - z1[j];
	nvec[j] = elenor[Gamma_i][j];
      }
      // check the boundary condition at the collocnode
      switch( ibval[Gamma_i] ) {
	case 0:         // velocity prescribed
	  Ibvi03 = 0;  
	  Ifadmii = 0;
	  break;
	case 1:         // pressure prescribed
	  Ibvi03 = 1;
	  Ifadmii = 0;
	  break;
	case 2:         // velocity and surface admittance prescribed
	case 5:
	  Ibvi03 = 0;  
	  Ifadmii = 1;
	  break;
      }
      zx[Gamma_i].set(0.0,0.0);
      for( s = 0; s < nsphere; s++) {
	sn = 0.0;
	sx = 0.0;
	for (j = 0; j < NDIM; j++) {
	  sn += clulevarry[nlevtop_].uvcsphe[s][j] * nvec[j];
	  sx += clulevarry[nlevtop_].uvcsphe[s][j] * x[j];
	}
	// the factor (1 + zBta3 * ik * sn) is used for all bcs
	sfact.set(0.0,waveNumbers_ * sn);  // ik * sn
	sfact = sfact * zBta3;
	sfact += 1.0;
	expfact.set( cos( waveNumbers_ * sx ), sin( waveNumbers_ * sx ) );
	sfact = sfact * expfact;
	if( Ifadmii )
	  sfact = sfact * zbval1[Gamma_i];
	zx[Gamma_i] += sfact * zGmat[n][s];
      }
    }
  }
  // scaling for the preconditioning
  for (i = 0; i < numRowsOfCoefficientMatrix_; i++)
    zx[i] *= scalefact[i];
}

void Cluster2LocalMtx(Complex &zBta3,bool allocateFMM) {
  /* local cluster to element expansion on the leaf level, generates
       a matrix that can be used for all iteration steps, does not
       do the multiplication with the FMM field, *nor* the scaling for
     the preconditiong and the multiplication with the weights on the sphere
  

     zBta3: BurtonMillerfactor
     
     Note: We do not use Gamma_i explicitely for zSmat, thus a multiplication
           with zSmat will result in a vector that has to be reordered
	   explicitely to match the ordering of the elements
	   
     Global: clulevarray: clustertree for all levels
             nlevtop_: leaf level
	     centel: midpoint of each element = collocnode
             zSmat: out: matrix containing the the value of the MLFMM at each
	               quadrature node of the sphere
	               zSmat[cluster][elem][sphere] 

     Note: Do not forget, the quadrature nodes on the sphere and the scaling
           for the incomplete LU factorization are *not* included in the matrix

	   zSmat: out: matrix containing the the value of the MLFMM at each
	   quadrature node of the sphere
	   zSmat[cluster][elem][sphere] 
  */
  int nsphere; // number of quad nodes on the sphere
  int Gamma_i,s;
  double z1[NDIM]; // coordinates of cluster center
  double nvec[NDIM]; // normal vector for a single element
  double x[NDIM]; //  elem midpoint - clustermidpoint
  int Ibvi03, Ifadmii; // boundary conditions of the target element
  double sn,sx; // point sphere times normal vector, or difference vector
  Complex expfact;
  int sentries = 0;
  if( allocateFMM ) {
    
    nsphere = clulevarry[nlevtop_].nPoinSpheLv;
    zSmat = new Complex*[ clulevarry[nlevtop_].nClustOLv ];
    if( zSmat == NULL ) {
      cerr << "Sorry could not allocate Smat\n";
      exit(-1);
    }
  
    for (int i = 0; i < clulevarry[nlevtop_].nClustOLv; i++) {
      zSmat[i] = new Complex[  clulevarry[nlevtop_].ClustArLv[i].NumOfEl * nsphere];
      sentries += clulevarry[nlevtop_].ClustArLv[i].NumOfEl * nsphere;
      if( zSmat[i] == NULL ) {
	cerr << "Sorry could not allocate Smat\n";
	exit(-1);
      }
      
    }

  }

  cout << "Matrix S: \n";
  cout << "Complex: "	 << sentries << "\n";
  for (int n = 0; n < clulevarry[nlevtop_].nClustOLv; n++) {
    for(int j = 0; j < NDIM; j++) {
      z1[j] = clulevarry[nlevtop_].ClustArLv[n].CoorCent[j];
    }
    for (int i = 0; i < clulevarry[nlevtop_].ClustArLv[n].NumOfEl; i++) {
      Gamma_i = clulevarry[nlevtop_].ClustArLv[n].NumsOfEl[i];
      for(int j = 0; j < NDIM; j++) {
	x[j] = z1[j] - centel[Gamma_i][j];
	nvec[j] = elenor[Gamma_i][j];
	// for debugging
	x[j] = -x[j];
      }
      // check the boundary condition at the collocnode
      switch( ibval[Gamma_i] ) {
	case 0:         // velocity prescribed
	  Ibvi03 = 0;  
	  Ifadmii = 0;
	  break;
	case 1:         // pressure prescribed
	  Ibvi03 = 1;
	  Ifadmii = 0;
	  break;
	case 2:         // velocity and surface admittance prescribed
	case 5:
	  Ibvi03 = 0;  
	  Ifadmii = 1;
	  break;
      }

      for(int s = 0; s < nsphere; s++) {
	sn = 0.0;
	sx = 0.0;
	for (int j = 0; j < NDIM; j++) {
	  sn += clulevarry[nlevtop_].uvcsphe[s][j] * nvec[j];
	  sx += clulevarry[nlevtop_].uvcsphe[s][j] * x[j];
	}
	// the factor (1 + zBta3 * ik * sn) is used for all bcs
	expfact.set( cos( sx * waveNumbers_), sin( sx * waveNumbers_) );
	zSmat[n][i*nsphere + s] = expfact;
	expfact.set(0.0,1.0);
	expfact = expfact * zBta3 * waveNumbers_ * sn + 1.0;
	zSmat[n][i*nsphere + s] *= expfact;
	if( Ifadmii ) {
	  cerr << "Sorry not implemented yet\n";
	  exit(-1);
	}
	// multiplication is done on the zG level
	//zSmat[n][i * nsphere + s].mul_r(clulevarry[nlevtop_].weisphe[s]);
      }
    }
  }
}

#if 0
void Expand2Eval(Complex** zG, Complex* zy, int* glob2local) {
  /* see the description below, this is the routine for the evalnodes
     Global:
     
     ipcluarry[ipc] evalcluster array, ipc in [0,ninpclus_)
     ninpclus_ number of evalclusters
  */
  
  for(int n = 0; n < ninpclus_; n++) {
    for (int s = 0; s < nsphere; s++)
      G[n][s] *= clulevarry[0].weisphe[s];
    
    for( int nn = 0; nn < NDIM; nn++) {
      z0(nn) = ipcluarry[n].CoorCent[nn];
    }
    
    for( int j = 0; j < ipcluarry[n].NumOfIps; i++ ) {
      eglobal = ipcluarry[n].NumsOfIps[j]; // global number of the node
      elocal = glob2local[eglobal];
      zy[elocal].set(0.0,0.0);
      for( int s = 0; s < nsphere; s++) {
	dist = 0.0;
	for(int nni = 0; n < NDIM; n++)
	  dist += (nodesCoordinates[eglobal][nn] - z0) * clulevarry[0].uvcsphe[s];
	dist *= waveNumbers_;
	zfact.set( cos(dist), sin(dist) );
	zy[elocal] += zfact * zG[n][s];
      }
    }
  }
}
#endif

void Expand2local(Complex** zGmat, Complex* zy) {
   /* as the multiplication involves 3 parts, it is done here, just not to
    forget the 3 parts:
    1) zGmat has to be scaled with the quadrature weights
    2) multiplication with zSmat
    3) the result has to be scaled with scalefact for the incomplete LU
    zSmat[cluster][elem][nodes_sphere]  in: complex matrix
   zGmat[cluster][nodes_sphere]
    scalefact[elems]
    zy
    Lets assume that zy is already initialized, thus += instead of = and
    beta = 1.0

    Global zSmat[Cluster][
           zNearScalefact[number of elements]
	   zrhs right hand side
   */
  int nelinclus; // number of elements in a single cluster
  //  int n,s,i,j;
  int nsphere; // number of quad nodes on the sphere
  int Gamma_i;
  Complex alpha,beta;
  Complex zdummy[numRowsOfCoefficientMatrix_]; // max possible length
  alpha.set(1.0,0.0);
  beta.set(0.0,0.0);

  nsphere = clulevarry[nlevtop_].nPoinSpheLv;
  for (int n = 0; n < clulevarry[nlevtop_].nClustOLv; n++) {
    for(int s = 0; s < nsphere ; s++) 
      zGmat[n][s] *= clulevarry[nlevtop_].weisphe[s];
    
    nelinclus = clulevarry[nlevtop_].ClustArLv[n].NumOfEl;
    
#ifdef USE_LAPACK
    cblas_zgemv(CblasRowMajor, CblasNoTrans, nelinclus, nsphere, &alpha, zSmat[n],  nsphere, zGmat[n], 1, &beta, zdummy, 1);
    
#else
    for (int i = 0; i < nelinclus; i++) {
      zdummy[i].set(0.0,0.0);
      //Gamma_i = clulevarry[nlevtop_].ClustArLv[n].NumsOfEl[i];
      for (int s = 0; s < nsphere; s++) {
	zdummy[i] += zSmat[n][i * nsphere + s] * zGmat[n][s];
      }
    }
#endif

    for (int j = 0; j < nelinclus; j++) {
      Gamma_i = clulevarry[nlevtop_].ClustArLv[n].NumsOfEl[j];
      zy[Gamma_i] += zdummy[j] * zNearscalefact[Gamma_i];
    }
    
    // this needs to be done later
    /* for ( i = 0; i < numRowsOfCoefficientMatrix_; i++) {
      zy[i] *= scalefact[i];
    }
    */
  }
}



void allocate_zFG() {
  /* Does what it says */
     
  int i,C_j;
  zF = new Complex**[numClusterLevels_];
  if( zF == NULL ) {
    cerr << "Sorry could not allocate zF\n";
    exit(-1);
  }
  zG = new Complex**[numClusterLevels_];
  if( zG == NULL ) {
    cerr << "Sorry could not allocate zG\n";
    exit(-1);
  }
  // loop over all levels of the tree
  for (i = 0; i < numClusterLevels_; i++) {
    zF[i] = new Complex*[clulevarry[i].nClustOLv];
    if( zF[i] == NULL ) {
      cerr << "Sorry could not allocate zF\n";
      exit(-1);
    }
    zG[i] = new Complex*[ clulevarry[i].nClustOLv ];
    if( zG[i] == NULL ) {
      cerr << "Sorry, Could not allocate zG\n";
      exit(-1);
    }
    
    for( C_j = 0; C_j < clulevarry[i].nClustOLv; C_j++) {
      zF[i][C_j] = new Complex[ clulevarry[i].nPoinSpheLv ];
      if( zF[i][C_j] == NULL ) {
	cerr << "Sorry could not allocate zF\n";
	exit(-1);
      }
      zG[i][C_j] = new Complex[ clulevarry[i].nPoinSpheLv ];
      if( zG[i][C_j] == NULL ) {
	cerr << "Sorry could not allocate zG\n";
	exit(-1);
      }
      for( int s = 0; s < clulevarry[i].nPoinSpheLv; s++) {
	zF[i][C_j][s].set(0.0,0.0);
	zG[i][C_j][s].set(0.0,0.0);
      }
    }
  }
}

void get_interactionlist() {
  // sets the interaction list on the clusterlevel thus
  // clulevarry[level].NumFanClus and
  // clulevarry[level].NumsFanClus

  // global clulevarry inout: clustertree
  
  int i,j,k,n,parent,nchild,child;
  int ninter, nearcl;
  bool nearclust;
  //  get the number of interaction cluster
  for ( n = 0; n < numClusterLevels_; n++ ) {
    if( n == 0 ) {
      /* in therory a good idea however all breaks loose if you properby
	 want to free all the memory
      // root has not interaction list, thus we just copy the data
      // from the far field clusters,
      // currently this leads to a memory leak, because there is some problem
      // with the fan clusters on the root level
      for( i = 0; i < clulevarry[n].nClustOLv; i++ ) {
	clulevarry[n].ClustArLv[i].NumFanClus = clulevarry[n].ClustArLv[i].NumFarClus;
	clulevarry[n].ClustArLv[i].NumsFanClus = clulevarry[n].ClustArLv[i].NumsFarClus;
      }
      */
      for( i = 0; i < clulevarry[n].nClustOLv; i++ ) {
	clulevarry[n].ClustArLv[i].NumFanClus = clulevarry[n].ClustArLv[i].NumFarClus;
	clulevarry[n].ClustArLv[i].NumsFanClus = new int[clulevarry[n].ClustArLv[i].NumFanClus];
	for (j = 0; j < clulevarry[n].ClustArLv[i].NumFanClus; j++ )
	  clulevarry[n].ClustArLv[i].NumsFanClus[j] = clulevarry[n].ClustArLv[i].NumsFarClus[j];
      }
    }
    else {
      // all other levels
      for( i = 0; i < clulevarry[n].nClustOLv; i++ ) {
	ninter = 0;
	parent = clulevarry[n].ClustArLv[i].nuFather;
	for( j = 0; j < clulevarry[n-1].ClustArLv[parent].NumNeaClus; j++ ) {
	  nearcl = clulevarry[n-1].ClustArLv[parent].NumsNeaClus[j];
	  //	  if( nearcl == parent )
	  //  continue;
	  ninter += clulevarry[n-1].ClustArLv[  nearcl ].n_Son;
	}
	// allocate worst case
	// there may be some elements in the interaction list number
	// that are nearfield elements
	// these should not be too many, alternatively one could use
	// std:vector<int> that allows to add entries at the end of
	// a vector
	//clulevarry[n].ClustArLv[i].NumFanClus = ninter;
	clulevarry[n].ClustArLv[i].NumsFanClus = new int[ninter];
	if( clulevarry[n].ClustArLv[i].NumsFanClus == NULL) {
	  cerr << "Sorry could not allocate interaction list\n";
	  exit(-1);
	}
      }
    }
  }
    // get their numbers
  for ( n = 1; n < numClusterLevels_; n++ ) {
    for( i = 0; i < clulevarry[n].nClustOLv; i++ ) {
      ninter = 0;
      parent = clulevarry[n].ClustArLv[i].nuFather;
      for( j = 0; j < clulevarry[n-1].ClustArLv[parent].NumNeaClus; j++ ) {
	// nearfield cluster of the parent
	nearcl = clulevarry[n-1].ClustArLv[parent].NumsNeaClus[j];
	// if( nearcl == parent )
	//  continue;
	for (nchild = 0; nchild < clulevarry[n-1].ClustArLv[nearcl].n_Son; nchild++) {
	  nearclust = false;
	  child = clulevarry[n-1].ClustArLv[nearcl].nuSon[nchild];
	  for (k = 0; k < clulevarry[n].ClustArLv[child].NumNeaClus; k++) {
	    if( clulevarry[n].ClustArLv[child].NumsNeaClus[k] == i ) {
	      nearclust = true;
	      break;
	    }
	  }
	  if( !nearclust ) {
	    clulevarry[n].ClustArLv[i].NumsFanClus[ninter] = child;
	    ninter++;
	  }
	}
      }
      clulevarry[n].ClustArLv[i].NumFanClus = ninter;
    }
  }
}
/*
void modify_rhs(zSparseVec& zFvec, Complex* zrhs) {
  // modifies the rhs with the MLFMM contributions cause by nonzero boundary
  // conditions
  // zFvec: in: sparse vector: for each cluster with elements with non zero bc
  //            an expansion is done for all points on the sphere
  // zrhs: inout: right hand side of the system

  int i, C_i, j, i1;
  for( i = 0; i < zFvec.nonzeroblocks; i++ ) {
    C_i = zFvec.clusterindx[i];
    for( j = 0; j < clulevarry[nlevtop_].ClustArLv[C_i].NumOfElems; j++) {
      Gamma_j = clulevarry[nlevtop_].ClustArLv[C_i].NumsOfElems[j];
      for ( i1 = 0; i1 < clulevarry[nlevtop_].ClustArLv[C_i].NumFanClus; i1++) {
	zG[nlevtop_][Gamma_i][s] += 
      }
    }
  }
}
*/
