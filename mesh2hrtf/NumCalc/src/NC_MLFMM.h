#ifndef BAMLFMM_h
#define BAMLFMM_h
#include<gsl/gsl_sf_bessel.h>
#include<gsl/gsl_sf_legendre.h>
#include<vector>
#include <bits/stdc++.h>
//#include<x86_64-linux-gnu/cblas64_mangling.h>
#include "NC_TypeDefinition.h"
#include "NC_Arrays.h"
#ifdef USE_GSL
#include<gsl/gsl_sf_bessel.h>
#include<gsl/gsl_sf_legendre.h>
#endif
#include<vector>
#include <bits/stdc++.h>
#ifdef USE_LAPACK
#include<x86_64-linux-gnu/cblas64.h>
#endif
extern Vector<double> Sourpoi3;
extern Vector<double> Norvci3;
extern void BAsinguII(ofstream&, Vector<Complex>&, const int&,
	       const int&, Vector<Complex>&, Matrix<double>&,
	       const int&, const int&); 
extern void BAreguII(ofstream&, Vector<Complex>&, const int&,
	    const int&, const int&, const int&,
	    Vector<Complex>&, Matrix<double>&,
	    const int&, const int&);

void LocalExpansionMat(int, bool);
void apply_localExpansion(Complex***, Vector<Complex>&);
void Get_Interpolation_Matrices(double**, int);
void UpPass(Complex***, const double*, int);
void UpPasslocal(Complex*, Complex**,const double*, int , int);
void DownPass(Complex**, Complex**, const double*, int);
void Cluster2Cluster(Complex*, Complex*, int, int, double*, int, double**);
void Cluster2Clustermat(int,bool);
void get_nearfield(int, double*, Complex&, bool );
void setup_preconditioning(double*);
bool IsNonZero_LU(int, int, int&, zSparsetype&);
bool IsNonZero_LU(int, int, zSparsetype&);
void sortArr(int*, int, int* );
void Cluster2Local(Complex**, double*, Complex*);
void Cluster2LocalMtx(Complex&,bool);
void Expand2local( Complex**, Complex*, bool useprecond = true);
void allocate_zFG();
void delete_zFG();
void Cleanup_MLFMM(bool);
void cluster2clusterVec();
void cluster2clusterlv(Complex***, Complex***, int);
void get_interactionlist();
void  addbc2rhs(zSparseVec*);
#endif
