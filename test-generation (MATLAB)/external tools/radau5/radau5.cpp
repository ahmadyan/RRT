/*********************************************************************
 * Demo.cpp
 *
 * This file shows the basics of setting up a mex file to work with
 * Matlab.  This example shows how to use 2D matricies.  This may
 * 
 * Keep in mind:
 * <> Use 0-based indexing as always in C or C++
 * <> Indexing is column-based as in Matlab (not row-based as in C)
 * <> Use linear indexing.  [x*dimy+y] instead of [x][y]
 *
 * For more information, see my site: www.shawnlankton.com
 * by: Shawn Lankton
 *
 ********************************************************************/
#include <matrix.h>
#include <mex.h>   
#include <iostream>
using namespace std;
#include "StiffIntegratorT.h"
#include "NonStiffIntegratorT.h"

/* Definitions to keep compatibility with earlier versions of ML */
#ifndef MWSIZE_MAX
typedef int mwSize;
typedef int mwIndex;
typedef int mwSignedIndex;

#if (defined(_LP64) || defined(_WIN64)) && !defined(MX_COMPAT_32)
/* Currently 2^48 based on hardware limitations */
# define MWSIZE_MAX    281474976710655UL
# define MWINDEX_MAX   281474976710655UL
# define MWSINDEX_MAX  281474976710655L
# define MWSINDEX_MIN -281474976710655L
#else
# define MWSIZE_MAX    2147483647UL
# define MWINDEX_MAX   2147483647UL
# define MWSINDEX_MAX  2147483647L
# define MWSINDEX_MIN -2147483647L
#endif
#define MWSIZE_MIN    0UL
#define MWINDEX_MIN   0UL
#endif
class mstream : public std::streambuf {
public:
protected:
  virtual std::streamsize xsputn(const char *s, std::streamsize n); 
  virtual int overflow(int c = EOF);
}; 

std::streamsize 
mstream::xsputn(const char *s, std::streamsize n) 
{
  mexPrintf("%.*s",n,s);
  return n;
}

int 
mstream::overflow(int c) 
{
    if (c != EOF) {
      mexPrintf("%.1s",&c);
    }
    return 1;
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
mstream matlab;
std::streambuf *outbuf = std::cout.rdbuf(&matlab); 
//declare variables
    mxArray *a_in_m, *b_in_m, *c_out_m, *d_out_m;
    const mwSize *dims;
    double *a, *b, *c, *d;
    int dimx, dimy, numdims;
    int i,j;

//associate inputs
    a_in_m = mxDuplicateArray(prhs[0]);
    b_in_m = mxDuplicateArray(prhs[1]);

//figure out dimensions
    dims = mxGetDimensions(prhs[0]);
    numdims = mxGetNumberOfDimensions(prhs[0]);
    dimy = (int)dims[0]; dimx = (int)dims[1];

//associate outputs
    c_out_m = plhs[0] = mxCreateDoubleMatrix(dimy,dimx,mxREAL);
    d_out_m = plhs[1] = mxCreateDoubleMatrix(dimy,dimx,mxREAL);

//associate pointers
    a = mxGetPr(a_in_m);
    b = mxGetPr(b_in_m);
    c = mxGetPr(c_out_m);
    d = mxGetPr(d_out_m);

//do something
    for(i=0;i<dimx;i++)
    {
        for(j=0;j<dimy;j++)
        {
            mexPrintf("element[%d][%d] = %f\n",j,i,a[i*dimy+j]);
            c[i*dimy+j] = a[i*dimy+j]+5; //adds 5 to every element in a
            d[i*dimy+j] = b[i*dimy+j]*b[i*dimy+j]; //squares b
        }
    }

    // dimension of problem
    const int n(8);
    // initial values for y
    double y[n] = {0.0, 6.0, 3.0, 3.0, 6.0, 3.0, 3.0, 0.0};
    // initial value for x
    double xbeg(0.0);
    // final value for x
    const double xend(0.03);
    // interval of x for printing output
    double dx(0.0025);
    // rtoler and atoler are scalars
    int itoler(0);
    // relative tolerance
    double *rtoler = new double(1.0e-5);
    // absolute tolerance
    double *atoler = new double(1.0e-11);
    // use SolutionOutput routine
    const int iout(1);
    // analytical Jacobian function provided
    const int ijac(1);
    // number of non-zero rows below main diagonal of Jacobian
    int mljac(1);
    // number of non-zero rows above main diagonal of Jacobian
    int mujac(2);
    // Mass matrix routine is provided
    const int imas(1);
    // number of non-zero rows below main diagonal of mass matrix
    int mlmas(1);
    // number of non-zero rows above main diagonal of mass matrix
    int mumas(1);
    
// Use default values (see header files) for these parameters:
    double hinit(0.0);
    double hmax(0.0);
    int nmax(0);
    double uround(0.0), safe(0.0), facl(0.0), facr(0.0);
    int nit(0);
    bool startn(false);
    int nind1(0), nind2(0), nind3(0), npred(0), m1(0), m2(0);
    bool hess(false);
    double fnewt(0.0), quot1(0.0), quot2(0.0), thet(0.0);

    StiffIntegratorT stiffT(n, y, xbeg, xend, dx, itoler, rtoler, atoler,
        iout, hinit, hmax, nmax, uround, safe, facl, facr, ijac, mljac,
        mujac, imas, mlmas, mumas, nit, startn, nind1, nind2, nind3, npred,
        m1, m2, hess, fnewt, quot1, quot2, thet);

    matlab << "\n\n*******Problem integrated with RADAU5*******\n\n";
    
    stiffT.Integrate();
    mexPrintf("result = %f\n",stiffT.NumFunction());
    // print statistics
    mexPrintf("result = %f\n",stiffT.NumFunction());
    mexPrintf("result = %f\n", stiffT.NumJacobian() );
    mexPrintf("result = %f\n", stiffT.NumStep() );
    mexPrintf("result = %f\n", stiffT.NumAccept() );
    mexPrintf("result = %f\n", stiffT.NumReject() );
    mexPrintf("result = %f\n", stiffT.NumDecomp() );
    mexPrintf("result = %f\n", stiffT.NumSol() ); 

    
    delete rtoler;
    delete atoler;

    std::cout.rdbuf(outbuf); 
    return;
}



void Function(double x, double *y, double *f)
{
    double ue, ub, uf, alpha, beta, r0, r1, r2, r3, r4, r5, r6, r7, r8, r9;
    double w, uet, fac1, fac2;

    ue = 0.1;
    ub = 6.0;
    uf = 0.026;
    alpha = 0.99;
    beta = 1.0e-6;
    r0 = 1000.0;
    r1 = r2 = r3 = r4 = r5 = r6 = r7 = r8 = r9 = 9000.0;
    w = 2.0*3.141592654*100.0;
    uet = ue*sin(w*x);
    fac1 = beta*(exp((y[3] - y[2])/uf) - 1.0);
    fac2 = beta*(exp((y[6] - y[5])/uf) - 1.0);

    f[0] = y[0]/r9;
    f[1] = (y[1] - ub)/r8 + alpha*fac1;
    f[2] = y[2]/r7 - fac1;
    f[3] = y[3]/r5 + (y[3] - ub)/r6 + (1.0 - alpha)*fac1;
    f[4] = (y[4] - ub)/r4 + alpha*fac2;
    f[5] = y[5]/r3 - fac2;
    f[6] = y[6]/r1 + (y[6] - ub)/r2 + (1.0 - alpha)*fac2;
    f[7] = (y[7] - uet)/r0;

    return;

} // Function

void Jacobian(double x, double *y, double **J)
{

    double uf, alpha, beta, r0, r1, r2, r3, r4, r5, r6, r7, r8, r9;
    double fac14, fac27;

    uf = 0.026;
    alpha = 0.99;
    beta = 1.0e-6;
    r0 = 1000.0;
    r1 = r2 = r3 = r4 = r5 = r6 = r7 = r8 = r9 = 9000.0;
    fac14 = beta*exp((y[3] - y[2])/uf)/uf;
    fac27 = beta*exp((y[6] - y[5])/uf)/uf;

    for (int j = 0; j < 8; j++) {
        J[0][j] = 0.0;
        J[1][j] = 0.0;
        J[3][j] = 0.0;
    }
    J[2][0] = 1.0/r9;
    J[2][1] = 1.0/r8;
    J[1][2] = -alpha*fac14;
    J[0][3] = alpha*fac14;
    J[2][2] = 1.0/r7 + fac14;
    J[1][3] = -fac14;
    J[2][3] = 1.0/r5 + 1.0/r6 + (1.0 - alpha)*fac14;
    J[3][2] = -(1.0 - alpha)*fac14;
    J[2][4] = 1.0/r4;
    J[1][5] = -alpha*fac27;
    J[0][6] = alpha*fac27;
    J[2][5] = 1.0/r3 + fac27;
    J[1][6] = -fac27;
    J[2][6] = 1.0/r1 + 1.0/r2 + (1.0 - alpha)*fac27;
    J[3][5] = -(1.0 - alpha)*fac27;
    J[2][7] = 1.0/r0;

    return;

} // Jacobian

void Mass(double **M)
{
    const double c1 = 1.0e-6;
    const double c2 = 2.0e-6;
    const double c3 = 3.0e-6;
    const double c4 = 4.0e-6;
    const double c5 = 5.0e-6;

    for (int i = 0; i < 8; i++) {
        M[0][i] = 0.0;
        M[2][i] = 0.0;
    }

    M[1][0] = -c5;
    M[0][1] = c5;
    M[2][0] = c5;
    M[1][1] = -c5;
    M[1][2] = -c4;
    M[1][3] = -c3;
    M[0][4] = c3;
    M[2][3] = c3;
    M[1][4] = -c3;
    M[1][5] = -c2;
    M[1][6] = -c1;
    M[0][7] = c1;
    M[2][6] = c1;
    M[1][7] = -c1;

    return;

} // Mass
