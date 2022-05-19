/*===========================================================================*\
 *                                                                            *
 *  File name:      NC_TypeDefinition.h                                       *
 *  Description:    Definition of types                                       *
 *  Author:         H. Ziegelwanger, W. Kreuzer and Z.S. Chen                 *
 *                                                                            *
 \*===========================================================================*/



#ifndef NC_TypeDefinition_h
#define NC_TypeDefinition_h



//================================= I N C L U D E S ==================================
// local includes                                                                   //
#include"NC_Macros.h"                                                               //
//                                                                                  //
// system includes                                                                  //
#include<cmath>                                                                     //
#include<iostream>                                                                  //
#include<assert.h>                                                                  //
using namespace std;                                                                //
//====================================================================================



// definition of complex class
class Complex {

protected:
	double re_, im_;

public:

	/* constructors and destructor */
	Complex() 
	{
		re_ = im_ = 0.0;
	}
	Complex(double rel, double img) 
	{
		re_ = rel;
		im_ = img;
	}
    Complex(const Complex &c)
    {
        re_ = c.re();
        im_ = c.im();
    }

    ~Complex(){}

    /* public functions */
    double re() const 
    {
        return(re_);
    }
    double im() const 
    {
        return(im_);
    }
    void set(double x, double y) 
    {
        re_ = x;
        im_ = y;
    }
    void set(const Complex z) 
    {
        re_ = z.re();
        im_ = z.im();
    }
    void mul_c(const Complex& x, const Complex& y)
    {
        re_ = x.re()*y.re() - x.im()*y.im();
        im_ = x.im()*y.re() + x.re()*y.im();
    }
    void mul_c(const Complex& x)
    {
        double tmp = re_;
        re_ = tmp*x.re() - im_*x.im();
        im_ = tmp*x.im() + im_*x.re();
    }
    void mul_r(const double& x)
    {
        re_ *= x;
        im_ *= x;
    }
    void mul_i(const double& y)
    {
        double tmp = re_;
        re_ = -im_*y;
        im_ = tmp*y;
    }
    void div_c(const Complex& x, const Complex& y)
    {
        double dy = y.re()*y.re() + y.im()*y.im();
        re_ = (x.re()*y.re() + x.im()*y.im())/dy;
        im_ = (x.im()*y.re() - x.re()*y.im())/dy;
    }
    void div_c(const Complex& x)
    {
        double tmp = re_, dx = x.re()*x.re() + x.im()*x.im();
        re_ = (tmp*x.re() + im_*x.im())/dx;
        im_ = (im_*x.re() - tmp*x.im())/dx;
    }
    void div_r(const double& x)
    {
        re_ /= x;
        im_ /= x;
    }
    void div_i(const double& y)
    {
        double tmp = re_;
        re_ = im_/y;
        im_ = -tmp/y;
    }
    void nega()
    {
        re_ = -re_;
        im_ = -im_;
    }
    void conju()
    {
        im_ = -im_;
    }
    Complex conjugate() const 
    {
        Complex z(re_, -im_);
        return(z);
    }
    double norm() const 
    {
        return sqrt(re_*re_ + im_*im_);
    }
    double qnorm() const 
    {
        return (re_*re_ + im_*im_);
    }
    double argu() const 
    {
        return atan2(im_, re_);
    }
    Complex zsqrt() const
    {
        double r, theta;
        r = sqrt(re_*re_ + im_*im_);
        theta = atan2(im_, re_)*0.5;
        r = sqrt(r);

        Complex z(r*cos(theta), r*sin(theta));
        return(z);
    }
    Complex zexp() const
    {
        double r = exp(re_);
        Complex z(r*cos(im_), r*sin(im_));
        return(z);
    }
    Complex zlog() const
    {
        Complex z(0.5*log(re_*re_ + im_*im_), atan2(im_, re_));
        return(z);
    }
    Complex zinv() const
    {
        double r = re_*re_ + im_*im_;
        Complex z(re_/r, -im_/r);
        return(z);
    }
    Complex znega() const
    {
        Complex z(-re_, -im_);
        return(z);
    }

    /* operators */
    Complex operator+(const Complex& x) const
    {
        Complex z(re_ + x.re(), im_ + x.im());
        return(z);
    }
    Complex operator+(const double& x) const
    {
        Complex z(re_ + x, im_);
        return(z);
    }
    Complex operator-(const Complex& x) const
    {
        Complex z(re_ - x.re(), im_ - x.im());
        return(z);
    }
    Complex operator-(const double& x) const
    {
        Complex z(re_ - x, im_);
        return(z);
    }
    Complex operator*(const Complex& x) const
    {
        Complex z(re_*x.re() - im_*x.im(), re_*x.im() + im_*x.re());
        return(z);
    }
    Complex operator*(const double& x) const
    {
        Complex z(re_*x, im_*x);
        return(z);
    }
    Complex operator/(const Complex& x) const
    {
        double tmp = x.re()*x.re() + x.im()*x.im();
        Complex z((re_*x.re() + im_*x.im())/tmp, (im_*x.re() - re_*x.im())/tmp);
        return(z); 
    }
    Complex operator/(const double& x) const
    {
        Complex z(re_/x, im_/x);
        return(z);
    }
    Complex& operator=(const int& i)
    {
        re_ = (double)i;
        im_ = 0.0; 
        return *this;
    }
    Complex& operator=(const double& d)
    {
        re_ = d;
        im_ = 0.0; 
        return *this;
    }
    Complex& operator=(const Complex& x)
    {
        re_ = x.re();
        im_ = x.im(); 
        return *this;
    }
    void operator+=(const Complex& x)
    {
        re_ += x.re(); 
        im_ += x.im(); 
    }
    void operator+=(const double& x)
    {
        re_ += x;  
    }
    void operator-=(const Complex& x)
    {
        re_ -= x.re(); 
        im_ -= x.im(); 
    }
    void operator-=(const double& x)
    {
        re_ -= x;  
    }
    void operator*=(const Complex& x)
    {
        double tmpr = re_;
        re_ = tmpr*x.re() - im_*x.im(); 
        im_ = tmpr*x.im() + im_*x.re(); 
    }
    void operator*=(const double& x)
    {
        re_ *= x; 
        im_ *= x; 
    }
    void operator/=(const Complex& x)
    {
        double tmpr = re_, tmp = x.re()*x.re() + x.im()*x.im();
        re_ = (tmpr*x.re() + im_*x.im())/tmp; 
        im_ = (im_*x.re() - tmpr*x.im())/tmp; 
    }
    void operator/=(const double& x)
    {
        re_ /= x; 
        im_ /= x; 
    }
};

// definition of the vector template
template <class T>
class Vector {

protected:

    int dim_;
    T *p_;

    // protected functions
    void Init(int n)
    {
        p_ = n > 0 ? new T[n] : NULL;
        dim_ = p_ != NULL ? n : 0;
    }
    void Reset()
    {
        if(p_ != NULL) delete [] p_;
        dim_ = 0;
        p_ = NULL;
    }

public:
    
    // constructors and destructor
    Vector() {
        Init(0);
    }
    Vector(int n)
    {
        Init(n);
    }
    Vector(int n, const T d)
    {
        Init(n);
        for(int i=0; i<dim_;) p_[i++] = d;
    }
    Vector(const Vector &B, int n) {
        Init(n);
        for(int i=0; i<B.size() && i<dim_; i++) p_[i] = B[i];
    }
    Vector(const Vector &B) {
        Init(B.size());
        for(int i=0; i<B.size(); i++) p_[i] = B[i];
    }
    ~Vector(){
        Reset();
    }
    
    // public functions
    int size() const
    {
        return(dim_);
    }
    int null() const
    {
        return dim_ == 0;
    }
    void sum(const Vector &A, const Vector &B) 
    {
        int n = min_(A.size(), B.size());
        for(int i=0; i<dim_ && i<n; i++) p_[i] = A[i] + B[i];
    }
    void sub(const Vector &A, const Vector &B) 
    {
        int n = min_(A.size(), B.size());
        for(int i=0; i<dim_ && i<n; i++) p_[i] = A[i] - B[i];
    }
    void mul_3(const Vector &A, const Vector &B) 
    {
        if(A.size() < 3 || B.size() < 3 || dim_ < 3) 
        {
            cout << "mul_3: dimention of the vectors must not less than 3!"
                << endl;
            exit(1);
        }
        p_[0] = A[1]*B[2] - B[1]*A[2];
        p_[1] = A[2]*B[0] - B[2]*A[0];
        p_[2] = A[0]*B[1] - B[0]*A[1];
    }

    // operators
    T&  operator[] (int i) const 
    {
        return p_[i];
    }
    Vector& operator=(const Vector &B)
    {
        for(int i=0; i<dim_ && i<B.size(); i++) p_[i] = B[i];
        return *this;
    }
    Vector& operator=(T a)
    {
        for(int i=0; i<dim_; i++) p_[i] = a; 
        return *this;
    }

    Vector operator+(const Vector &B);

    Vector operator-(const Vector &B);

    T operator*(const Vector &B);
    T operator*(const T a);
};

template <class T>
Vector<T> Vector<T>::operator+(const Vector<T> &B) 
{
    Vector E(dim_);
    for(int i=0; i<dim_ && i<B.size(); i++) E.p_[i] = p_[i] + B[i];
    return E;
}


template <class T>
Vector<T> Vector<T>::operator-(const Vector<T> &B) 
{
    Vector E(dim_);
    for(int i=0; i<dim_ && i<B.size(); i++) E.p_[i] = p_[i] - B[i];
    return E;
}

template <class T>
T Vector<T>::operator*(const Vector<T> &B) 
{
    int n = 0;
    T d;
    d = n;
    for(int i=0; i<dim_ && i<B.size(); i++) d += p_[i]*B[i];
    return(d);
}


template <class T>
T Vector<T>::operator*(const T a) 
{
    int n = 0;
    T d;
    d = n;
    for(int i=0; i<dim_ ; i++) d += p_[i]*a;
    return(d);
}

// definition of the matrix template
template <class T>
class Matrix {

protected:
    Vector<T> v_;
    int dim_[2]; // number of rows and columns

public:
    // constructors and destructor
    Matrix() : v_() 
    {
        dim_[0] = 0;
        dim_[1] = 0;
    }
    Matrix(int nr, int nc) : v_(nr*nc) 
    {
        dim_[0] = nr;
        dim_[1] = nc;
    }
    Matrix(int nr, int nc, T s) : v_(nr*nc, s) 
    {
        dim_[0] = nr;
        dim_[1] = nc;
    }
    Matrix(Vector<T> &B, int nr, int nc) : v_(B, nr*nc)
    {
        dim_[0] = nr;
        dim_[1] = nc;
    }
    Matrix(Matrix &B, int nr, int nc) : v_(B.v_, nr*nc)
    {
        dim_[0] = nr;
        dim_[1] = nc;
    }
    ~Matrix()
    {
        dim_[0] = dim_[1] = 0;
    }

    // public functions
    int size() const
    {
        return(dim_[0]*dim_[1]);
    }
    int null() const
    {
        return(dim_[0] == 0 || dim_[1] == 0);
    }
    void sum(const Matrix &A, const Matrix &B)
    {
        int n = min_(A.size(), B.size()), m = dim_[0]*dim_[1];
        for(int i=0; i<m && i<n; i++) this->v_[i] = A[i] + B[i];
    }
    void sub(const Matrix &A, const Matrix &B) {
        int n = min_(A.size(), B.size()), m = dim_[0]*dim_[1];
        for(int i=0; i<m && i<n; i++) this->v_[i] = A[i] - B[i];
    }

    void mul(const Matrix &, const Matrix &);

    void mul_trans(const Matrix &, const Matrix &);

    // operators
    T& operator[](int i) const
    {
        return (v_[i]);
    }
    T& operator()(int i, int j) const
    {
        return (v_[i*dim_[1] + j]);  // the elements arranged in row oder form
    }
    Matrix& operator=(const Matrix &B) {  
        for(int i=0; i<dim_[0]*dim_[1] && i<B.size(); i++) v_[i] = B[i];
        return(*this);
    }
    Matrix& operator=(T a) {
        for(int i=0; i<dim_[0]*dim_[1]; i++) v_[i] = a; 
        return(*this);
    }
};

// A*B
template <class T>
void Matrix<T>::mul(const Matrix<T> &A, const Matrix<T> &B) {
    assert(dim_[0] <= A.dim_[0]);
    assert(dim_[1] <= B.dim_[1]);
    int n = min_(A.dim_[1], B.dim_[0]), i, j ,k, l = 0, m;
    T dt;
    for(i=0; i<dim_[0]; i++) {
        for(j=0; j<dim_[1]; j++) {
            m = 0;
            dt = m;
            for(k=0; k<n; k++) dt += A(i, k)*B(k, j);
            v_[l] = dt;
            l++;
        }
    }
}

// A*B(T)
template <class T>
void Matrix<T>::mul_trans(const Matrix<T> &A, const Matrix<T> &B) {
    assert(dim_[0] <= A.dim_[0]);
    assert(dim_[1] <= B.dim_[0]);
    int n = min_(A.dim_[1], B.dim_[1]), i, j ,k, l = 0, m;
    T dt;
    for(i=0; i<dim_[0]; i++) {
        for(j=0; j<dim_[0]; j++) {
            m = 0;
            dt = m;
            for(k=0; k<n; k++) dt += A(i, k)*B(j, k);
            v_[l] = dt;
            l++;
        }
    }
}

// definition of the template for performing the factorization of an unsymmetric matrix stored as a vector
template <class T>
void Tfactor_usy(T *Coef, const int& n_row_col)
{
    int i, j, k, idia, i1, i0, j0;
    T C_dia, C_ji;
   
    i0 = 0;  // address of the first compunent of the row i
    for(i=0; i<n_row_col - 1; i++) // loop over rows
    {
        idia = i0 + i;     // address of the diagonal compunent of the row i
        C_dia = Coef[idia];
        i1 = i + 1;
        j0 = i1*n_row_col; // address of the first component of the row j
        for(j=i1; j<n_row_col; j++) // loop over rows
        {
            Coef[j0 + i] /= C_dia;
            C_ji = Coef[j0 + i];
            for(k=i1; k<n_row_col; k++) // loop over columns
            {
                Coef[j0 + k] -= Coef[i0 + k]*C_ji;
            }
            j0 += n_row_col;
        }
        for(k=i1; k<n_row_col; k++) Coef[i0 + k] /= C_dia;
        i0 += n_row_col;
    }   /* end of loop I */
}

// definition of the template for performing the forward elimination and backward substitution
template <class T>
void Tfbelim(T *Coef, T *rhsv, const int& n_row_col)
{
    int i, j, i0;

    // forward elimination
    i0 = n_row_col;
    for(i=1; i<n_row_col; i++) // loop over rows
    {
        for(j=0; j<i; j++)     // loop over columns
        {
            rhsv[i] -= Coef[i0 + j]*rhsv[j];
        }
        i0 += n_row_col;
    }

    for(i=0; i<n_row_col; i++) rhsv[i] /= Coef[i*n_row_col + i];

    // backward substitution
    for(i=n_row_col - 1; i>0; i--) // loop over columns
    { 
        for(j=0; j<i; j++) rhsv[j] -= rhsv[i]*Coef[j*n_row_col + i];
    }
}

// structure of an element cluster
struct ElCluster
{
    int NuElGr;            // number of the element group to which the cluster belongs
    int NumOfEl;           // number of elements in the cluster
    int listElementPropertyEl;          // = 0: surface elements
                           // = 1: middle face elements
                           // = 2: evaluation elements
    bool IfMonoEl;         // = true:  all elements of the cluster are of the same number of nodes 
                           // = false: the cluster contains elements with different number of nodes
    bool IfAdmiBc;         // = true: admittance boundary condition are prescribed
    int NumOfDOFs;         // number of the unknown DOFs of the cluster
    int NDOFsPeEl;         // number of the unknown DOFs per element
    int *NumsOfEl;         // numbers of the elements in the cluster
    double CoorCent[3];	   // coordinates of the center of the cluster
    double RadiClus;       // radius of the cluster
    int NumNeaClus;        // number of near clusters
    int NumFarClus;        // number of far clusters
    int *NumsNeaClus;      // numbers of the near clusters
    int *NumsFarClus;      // numbers of the far clusters
    int OriClust;          // number of the corresponding original cluster
    int nuref;             // number of the reflection (= one number in [0, 1, ..., 7])
    int rffac;             // factor for velocity boundary condition (= 1 or -1)
    bool ifmirro;          // = false: element conectivities are ident to original ones
                           // = true:  these are mirror image of the original ones
    bool ifrfdi[3];        // ifrfdi[i] = treu: elements must be reflected in the
                           //                   ith-direction
                           //           = false: do not reflected in this direction
    int nuFather;          // number of the father cluster
    int n_Son;             // number of the son clusters
    int nuSon[8];          // numbers of the son clusters
    int nuLev;             // number of the level to which the cluster belongs
    int NumFanClus;        // number of clusters, that is located in the far field of
                           // current cluster but the fathers are located in the near
                           // field
    int *NumsFanClus;      // numbers of the above clusters
};

// structure of a cluster of internal points
struct IpCluster
{
    int NumOfIps;          // number of internal points in the cluster
    int *NumsOfIps;        // local numbers of internal points in the cluster
    double CoorCent[3];	   // coordinates of the center of the cluster
    double RadiClus;       // radius of the cluster
};

// structure of a level of a cluster tree
struct ClusterLev
{
    // all the variables and arrays are defined for a given level in the tree
    int nClustSLv;         // number of all clusters 
    int nClustOLv;         // number of the original clusters
    double RadiMaxLv;      // maximum radius of clusters
    double RadiAveLv;      // average radius of clusters
    double RadiMinLv;      // minimum radius of clusters
    int nExpaTermLv;       // number of terms in the expresions of the Green functions
    int nPoinSpheLv;       // number of integration points on the unit sphere surface
    int nPoinThetLv;       // number of integration points in the thetha direction
    int nPoinPhiLv;        // number of integration points in the phi direction
    int Pow_532Lv[3];      // working array: nPoinPhiLv expressed in product of pwoers of
                           // 5, 3 and 2 is stored in this array
    int iplacsum;          // = Pow_532Lv[0] + Pow_532Lv[1] + Pow_532Lv[2]
    int **iplacpova;       // working array used by FFT 
    int **iregister;       // working array used by FFT
    double *CrdGauLv;      // coordinates of Gaussean points in the thetha-direction
    double *WeiGauLv;      // weights of Gaussean points in the thetha-direction
    double *XvThetLv;      // X-vector in the thetha-direction for spline interpolaton
    double *XvPhiLv;       // X-vector in the phi-direction for spline interpolaton
    double *CmThetLv;      // Coefficient matrix of the m-vector in the thetha-direction
                           // for spline interpolaton
    double *CmPhiLv;	   // Coefficient matrix of the m-vector in the phi-direction
                           // for spline interpolaton
    double *LmThetLv;      // Lagrange matrix for the level (stored as a vector)
    Complex *ZmThetLv;     // the m-vector in the thetha-direction for spline interpolaton
    Complex *ZmPhiLv;      // the m-vector in the phi-direction for spline interpolaton
    Complex *PbarnmLv;     // P_bar^{m}_{n} at all Gaussean points in the theta-direction
    Complex *PbarfaLv;     // P_bar^{m}_{n} at all Gaussean points of the father level
    Complex *zwkT;		   // working array for T-vectors
    Complex *zwkS;		   // working array for S-vectors
    double **uvcsphe;      // coordinates of the integration points on the unit sphere
    double *weisphe;       // weights of these integration points
    Complex **zIntpMtx;    // interpolation matrix from the current level to the father
                           // level
    Complex **zFiltMtx;    // filter matrix from the current level to the son level
    ElCluster *ClastArLv;  // array of clusters
};

// structure of the D-matrices
struct D_mtx_lev
{
    int nRowsD;			// number of rows of the matrx D
    int nEntriesD;		// number of nonzero entries of the matrx D
    int *jcolDmxLv;		// column number of each nonzero of the D-matrix
    int *irowDmxLv;		// number of the first nonzero of each row of the D-matrix
    Complex *zDmxLv;	// array for nonzeros of the D-matrix
};

// structure of the T-matrices
struct T_mtx_lev
{
    int nRowsT;			// number of rows of the matrx T
    int nEntriesT;		// number of nonzero entries of the matrx T
    int *jcolTmxLv;		// column number of each nonzero of the T-matrix
    int *irowTmxLv;		// number of the first nonzero of each row of the T-matrix
    Complex *zTmxLv;	// array for nonzeros of the T-matrix
    Complex *zTvcLv;    // the T-vector
};

// structure of the T-matrices
struct S_mtx_lev
{
    int nRowsS;			// number of rows of the matrx S
    int nEntriesS;		// number of nonzero entries of the matrx S
    int *jcolSmxLv;		// column number of each nonzero of the S-matrix
    int *irowSmxLv;		// number of the first nonzero of each row of the S-matrix
    Complex *zSmxLv;	// array for nonzeros of the S-matrix
};

// structure of the input line for boundary conditions
struct inputLineBoundaryCondition
{
    int		nLow,		// range of elements
            nHigh;
    //jw:ansi
    string	sKeyword;	// IMPE|ADMI|VELO|PRES|TRAD|TRIM
    double	dReal,		// real and imaginary boundary conditions
            dImag,	
            dNegReal,	// for IMPE and ADMI, negative conditions necessary
            dNegImag;
    bool	bRef;		// are curves referenced (default == no)
    bool    bNega;      // if the keyword NEGA appears
    int		nRealRef,	// reference to curve number (possible for IMPE|ADMI|TRIM|TRAD)
            nImagRef,
            nNegRealRef,
            nNegImagRef;
};



extern void NC_ReflectElementFMM(int&, bool&, Vector<bool>&,
                       const int&, const int&, Vector<double>&, 
                       Vector<double>&, Vector<Complex>&, Matrix<double>&);
extern void NC_ComputeGlobalCoordinates(Vector<double>&, Vector<double>&, Matrix<double>&,
                     const double&, const double&, const int&);
extern double NC_ComputeUnitVectorOnElement(Vector<double>&, Vector<double>&, Vector<double>&,
                       Vector<double>&, Vector<double>&, Matrix<double>&,
                       const double&, const double&, const int&);
extern double NC_ReflectPointVector(Vector<double>&, const int&, double&,
                      const int&);
extern void NC_ComputeElementData(const int&, const int&,
                    int&, int&, int&,
                    Vector<int>&, Vector<double>&, Vector<double>&,
                    Vector<Complex>&, Vector<Complex>&, Matrix<double>&); 
extern void NC_ReflectElement(const int&, const int&, const int&,
                    int&, Vector<double>&, Vector<double>&,
                    Vector<Complex>&, Matrix<double>&);
extern int NC_ComputeGausseanPoints(double*, double*, double*,
                    const int&, const int&);
extern double NC_ComputeParameters(Vector<double>&, Vector<double>&, Vector<double>&,
                      Vector<double>&, Vector<double>&, Vector<double>&,
                      Vector<double>&, Complex&, Matrix<double>&,
                      const int&, Vector<Complex>&, const int&,
                      const double&, const double&);
extern int NC_GenerateSubelements(ofstream&, const int&, Vector<double>&, const int&,
                    Matrix<double>&, const double&, const int&,
                    Vector<double>&, Vector<double>&, const int&,
                    Vector<double>&, Vector<double>&, Vector<double>&,
                    Vector<int>&);
void NC_ComputeUnitVectorOnUnitSphere(ofstream&, const int&, const int&, double*, double**);
void NC_SphericalHankel(ofstream&, Vector<Complex>&, const int&,const double&);
void NC_LegendrePolynomes(ofstream&, Vector<double>&, const int&,const double&);
extern void NC_Error_Exit_0(ofstream&, const string&);
extern void NC_Error_Warning_0(ofstream&, const string&);
extern void NC_Error_Exit_0(const string&);
extern void NC_Error_Warning_0(const string&);
void NC_Error_Exit_1(ofstream&, const string&, const string&,const int&);
void NC_Error_Exit_1(ofstream&, const string&, const string&,const double&);
void NC_Error_Exit_1(ofstream&, const string&, const string&);
void NC_Error_Warning_1(ofstream&, const string&, const string&,const int&);
void NC_Error_Warning_1(ofstream&, const string&, const string&,const double&);
void NC_Error_Exit_1(const string&, const string&, const int&);
void NC_Error_Exit_1(const string&, const string&, const double&);
void NC_Error_Exit_1(const string&, const string&);
void NC_Error_Warning_1(const string&, const string&, const int&);
void NC_Error_Warning_1(const string&, const string&, const double&);
void NC_Error_Exit_2(ofstream&, const string&, const string&,const int&, const string&, const int&);
void NC_Error_Exit_2(ofstream&, const string&, const string&,const double&, const string&, const double&);
void NC_Error_Warning_2(ofstream&, const string&, const string&,const int&, const string&, const int&);
void NC_Error_Warning_2(ofstream&, const string&, const string&,const double&, const string&, const double&);
void NC_Error_Warning_2(ofstream&, const string&, const string&,const int&, const string&, const double&);
void NC_Error_Exit_2(const string&, const string&, const int&,const string&, const int&);
void NC_Error_Exit_2(const string&, const string&, const double&,const string&, const double&);
void NC_Error_Warning_2(const string&, const string&, const int&,const string&, const int&);
void NC_Error_Warning_2(const string&, const string&, const double&,const string&, const double&);
extern ostream& operator<<(ostream &s, const Complex &u);
extern double BLdispoi(Vector<double>&, Vector<double>&);
extern Complex BLzsqrt(const Complex&);
void BLGauCooWe(ofstream&, const int&, Vector<double>&,Vector<double>&);
void NC_IterativeSolverCGS(ofstream&);
Vector<double> BLAverageVc_3(Vector<double>&, Vector<double>&, Vector<double>&);
double interpo_line(double&, double freq_poi[], double fac1_poi[], int&);



#endif
