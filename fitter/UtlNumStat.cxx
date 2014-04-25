//____________________________________________________________________________
/*

Copyright (c) 2010-2014, The VALOR group. All rights reserved.

ChangeLog:

*/
//____________________________________________________________________________

//#define __show_cholesky_decomposition_mesg_debug__
//#define __show_cholesky_gen_corr_param_variations_mesg_debug__
//#define __show_xml_mesg_info__
//#define __show_xml_mesg_debug__
//#define __show_arr_idx_mesg_debug__

#include <cassert>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>

#include <TROOT.h>
#include <TMath.h>
#include <TStyle.h>
#include <TColor.h>
#include <TGraph.h>
#include <TH1.h>
#include <TH2.h>
#include <TAxis.h>
#include <TVectorD.h>
#include <TMatrixD.h>
#include <TRandom3.h>
#include <TXMLEngine.h>
#include <TFile.h>
#include <TSystem.h>

//~ #include "Conventions/Constants.h"
//~ #include "Conventions/PhysHypothesisId.h"
//~ #include "Conventions/Opt.h"
//~ #include "Utilities/UtlNumStat.h"
//~ #include "Utilities/Rnd.h"

#include "ParameterValues.h"
#include "Conventions.h"
#include "Opt.h"
#include "UtlNumStat.h"
#include "Rnd.h"

#ifdef __GSL_ENABLED__ 
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>
#endif

using std::cout;
using std::cerr;
using std::endl;
using std::setw;
using std::setfill;
using std::fixed;
using std::scientific;
using std::setprecision;
using std::ofstream;
using std::ostringstream;

using namespace t2kana;
using namespace t2kana::conv;

//____________________________________________________________________________
double t2kana::util::DeltaChisquare(
  ConfLevel_t cl, int nfitparams) 
{
  switch(cl) {
      case (kCLUndefined) : 
      {
	return -99999;  
        break;
      }
      case (kCL68) :
      {
        if      (nfitparams == 1)  { return  1.00;  }
        else if (nfitparams == 2)  { return  2.30;  }
        else if (nfitparams == 3)  { return  3.53;  }
        else                       { return -99999; }
        break;
      }
      case (kCL90) :
      {
        if      (nfitparams == 1)  { return  2.71;  }
        else if (nfitparams == 2)  { return  4.61;  }
        else if (nfitparams == 3)  { return  6.25;  }
        else                       { return -99999; }
        break;
      }
      case (kCL95) :
      {
        if      (nfitparams == 1)  { return  3.84;  }
        else if (nfitparams == 2)  { return  5.99;  }
        else if (nfitparams == 3)  { return  7.82;  }
        else                       { return -99999; }
        break;
      }
      case (kCL99) :
      {
        if      (nfitparams == 1)  { return  6.63;  }
        else if (nfitparams == 2)  { return  9.21;  }
        else if (nfitparams == 3)  { return 11.34;  }
        else                       { return -99999; }
        break;
      }
      default:
        return -99999;     
    }
  return -99999;
} 
//____________________________________________________________________________
double t2kana::util::ChisqGauss(
  double ndata, double nmc)
{
  double chisq = 0.;
  if(nmc>0) {
    chisq = TMath::Power(ndata-nmc,2.)/nmc;
  }
  return chisq;
}
//____________________________________________________________________________
double t2kana::util::ChisqPoisson(
  double ndata, double nmc)
{
  if(nmc <= 0.)
    nmc = kSetNegEventsTo;
  if(ndata <= 0.)
    return 2.*(nmc-ndata);
  return 2. * (nmc-ndata) + 2. * ndata * TMath::Log(ndata/nmc);
}
//____________________________________________________________________________
double t2kana::util::ArrayChisqGauss(
  const VALORArray * data, const VALORArray * mc, VALORArray * output)
{
  int size = data->GetSize();
  assert(size == mc->GetSize());
  
  const double * data_arr = data->GetArray();
  const double * mc_arr   = mc  ->GetArray();
  double * out_arr = NULL;
  if (output) {
    assert(size==output->GetSize());
    out_arr = output->GetArray();
  }
  
  double total_chisq = 0., chisq=0;
  double nmc,ndata;
  for (int ii = 0 ; ii < size ; ii++) {
    nmc   = mc_arr  [ii];
    ndata = data_arr[ii];
    if(nmc>0) {
      chisq = TMath::Power(ndata-nmc,2.)/nmc;
    }
    else {
      chisq=0;
    }
    if (output) out_arr[ii] = chisq;
    total_chisq+=chisq;
  }
  return total_chisq;
}
//____________________________________________________________________________
double t2kana::util::ArrayChisqPoisson(
  const VALORArray *  data, const VALORArray *  mc, VALORArray * output)
{
  // Warning - empty bins (0 MC, 0 data) have a small but non-zero chisq (kSetNegEventsTo).
  int size = data->GetSize();
  assert(size == mc->GetSize());
  
  const double * data_arr = data->GetArray();
  const double * mc_arr   = mc  ->GetArray();
  double * out_arr = NULL;
  if (output) {
    assert(size==output->GetSize());
    out_arr = output->GetArray();
  }
  
  double total_chisq = 0., chisq=0;
  double nmc,ndata;
  for (int ii = 0 ; ii < size ; ii++) {
    nmc   = mc_arr  [ii];
    ndata = data_arr[ii];
    chisq = 0.;
    if(nmc <= 0.)
      nmc = kSetNegEventsTo;
    if(ndata <= 0.)
      chisq = 2.*(nmc-ndata);
    else
      chisq = 2. * (nmc-ndata) + 2. * ndata * TMath::Log(ndata/nmc);
    
    if (output) out_arr[ii] = chisq;
    total_chisq+=chisq;
  }
  return total_chisq;
}
//____________________________________________________________________________
TMatrixD t2kana::util::CholeskyDecomposition(TMatrixD& cov_matrix)
{
// Perform a Cholesky decomposition of the input covariance matrix and
// return the lower triangular matrix
//
   int ncols = cov_matrix.GetNcols();
   int nrows = cov_matrix.GetNrows();

   assert(ncols==nrows);

   int n = nrows;

   TMatrixD L(n, n);

   for (int i = 0; i < n; ++i) {

     // calculate the diagonal term first
     L(i,i) = cov_matrix(i,i);
     for (int k = 0; k < i; ++k) {
        double tmp = L(k,i);
        L(i,i) -= tmp*tmp;
     }//k
     if(L(i,i) <= 0) {
       cerr << "Decomposed covariance matrix not positive-definite" << endl;
       exit(1);
     }
     L(i,i) = TMath::Sqrt(L(i,i));
     // then the off-diagonal terms
     for (int j = i+1; j < n; ++j) {
        L(i,j) = cov_matrix(i,j);
        for (int k = 0; k < i; ++k) {
           L(i,j) -= L(k,i)*L(k,j);
        }
        L(i,j) /= L(i,i);
     }//j
  }//i

  // create the transpose of L
  TMatrixD LT(TMatrixD::kTransposed,L);

#ifdef __show_cholesky_decomposition_mesg_debug__
  cout << "Cholesky matrix L^{T}" << endl;
  for(int ip=0; ip<n; ip++) {
    for(int jp=0; jp<n; jp++) {
      cout << setw(15) << setfill(' ') << setprecision(6) << LT[ip][jp] << " | ";
    }//jp
    cout << endl;
  }//ip
#endif

  return LT;
}
//____________________________________________________________________________
TVectorD t2kana::util::CholeskyGenerateCorrelatedParams (
    TMatrixD& cholesky_triangular, TVectorD& mean_params) 
{
// Generate a vector of correlated params

  int ncols = cholesky_triangular.GetNcols();
  int nrows = cholesky_triangular.GetNrows();
  int npars = mean_params.GetNrows();

  if(ncols != nrows) {
    cerr << "util::CholeskyGenerateCorrelatedParams >> Number of columns " << ncols << " & rows " << nrows << " disagree. Exiting..." << endl;
    exit(1);
  }
  if(npars != nrows) {
    cerr << "util::CholeskyGenerateCorrelatedParams >> Number of parameters " << npars << " & rows " << nrows << " disagree. Exiting..." << endl;
    exit(1);
  }

  int n = nrows;

  // create a vector of unit Gaussian variables
  // and multiply by Lt to introduce the appropriate correlations
  TVectorD g(n);
  for (int k = 0; k < n; ++k) {
    g(k) = Rnd::Instance()->RndOther().Gaus();
  }
  g *= cholesky_triangular;

#ifdef __show_cholesky_gen_corr_param_variations_mesg_debug__
  ofstream foutp1("cholesky_uncorrelated_params.txt",fstream::app);
  for (int i = 0; i < n; ++i) {
    foutp1 << g[i] << "\t";
  }
  foutp1 << endl;
  foutp1.close();
#endif

  // add the mean value offsets and store the results
  TVectorD correlated_params(n);
  for (int i = 0; i < n; ++i) {
     double v = mean_params[i];
     v += g(i);
     correlated_params[i] = v;
  }

#ifdef __show_cholesky_gen_corr_param_variations_mesg_debug__
  ofstream foutp2("cholesky_correlated_params.txt",fstream::app);
  for (int i = 0; i < n; ++i) {
    foutp2 << correlated_params[i] << "\t";
  }
  foutp2 << endl;
  foutp2.close();
#endif

  return correlated_params;
}
//____________________________________________________________________________
TVectorD t2kana::util::CholeskyGenerateCorrelatedParams (
 TMatrixD& cholesky_triangular, TVectorD& mean_params, TVectorD& g_uncorrelated) 
{
// Generate a vector of correlated params

  int ncols = cholesky_triangular.GetNcols();
  int nrows = cholesky_triangular.GetNrows();
  int npars = mean_params.GetNrows();
  int nunco = g_uncorrelated.GetNrows();

  if(ncols != nrows) {
    cerr << "util::CholeskyGenerateCorrelatedParams >> Number of columns " << ncols << " & rows " << nrows << " disagree. Exiting..." << endl;
    exit(1);
  }
  if(npars != nrows) {
    cerr << "util::CholeskyGenerateCorrelatedParams >> Number of parameters " << npars << " & rows " << nrows << " disagree. Exiting..." << endl;
    exit(1);
  }
  if(nunco != nrows) {
    cerr << "util::CholeskyGenerateCorrelatedParams >> Number of uncorrelated parameters " << nunco << " & rows " << nrows << " disagree. Exiting..." << endl;
    exit(1);
  }

  int n = nrows;

  // create a vector of unit Gaussian variables
  // and multiply by Lt to introduce the appropriate correlations
  g_uncorrelated *= cholesky_triangular;

#ifdef __show_cholesky_gen_corr_param_variations_mesg_debug__
  ofstream foutp1("cholesky_uncorrelated_params.txt",fstream::app);
  for (int i = 0; i < n; ++i) {
    foutp1 << g_uncorrelated[i] << "\t";
  }
  foutp1 << endl;
  foutp1.close();
#endif

  // add the mean value offsets and store the results
  TVectorD correlated_params(n);
  for (int i = 0; i < n; ++i) {
     double v = mean_params[i];
     v += g_uncorrelated(i);
     correlated_params[i] = v;
  }

#ifdef __show_cholesky_gen_corr_param_variations_mesg_debug__
  ofstream foutp2("cholesky_correlated_params.txt",fstream::app);
  for (int i = 0; i < n; ++i) {
    foutp2 << correlated_params[i] << "\t";
  }
  foutp2 << endl;
  foutp2.close();
#endif

  return correlated_params;
}
//____________________________________________________________________________
TVectorD t2kana::util::CholeskyGenerateCorrelatedParamVariations (
    TMatrixD& cholesky_triangular) 
{
// Generate a vector of correlated params

  int ncols = cholesky_triangular.GetNcols();
  int nrows = cholesky_triangular.GetNrows();

  assert(ncols==nrows);

  int n = nrows;

  // create a vector of unit Gaussian variables
  // and multiply by Lt to introduce the appropriate correlations
  TVectorD g(n);
  for (int k = 0; k < n; ++k) {
    g(k) = Rnd::Instance()->RndOther().Gaus();
  }
  g *= cholesky_triangular;

#ifdef __show_cholesky_gen_corr_param_variations_mesg_debug__
  cout << "Util >> Tweak vector size = " << g.GetNrows() << endl;
#endif

  return g;
}
//____________________________________________________________________________
TVectorD t2kana::util::CholeskyCalculateCorrelatedParamVariations (
    TMatrixD& cholesky_triangular, TVectorD & g_uncorrelated) 
{
// Generate a vector of correlated params

  int ncols = cholesky_triangular.GetNcols();
  int nrows = cholesky_triangular.GetNrows();
  int npars = g_uncorrelated.GetNrows();

  assert(ncols==nrows);
  assert(npars==nrows);

  // create a vector of unit Gaussian variables
  // and multiply by Lt to introduce the appropriate correlations
  TVectorD g(g_uncorrelated);
  g *= cholesky_triangular;

#ifdef __show_cholesky_gen_corr_param_variations_mesg_debug__
  cout << "Util >> Tweak vector size = " << g.GetNrows() << endl;
#endif

  return g;
}
//____________________________________________________________________________
double t2kana::util::RandomizeMixing12(void)
{
  double mix = 0;

  NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
  switch(omp) {
   case(kNuMixSinsq2Theta):
   {
     double mean  = conv::kGlobFitSinsq2Theta12; 
     double sigma = conv::kGlobFitSinsq2Theta12 * conv::kGlobFitErrSinsq2Theta12;
     double sinsq2theta = Rnd::Instance()->RndToyMCosc().Gaus(mean,sigma);
     mix = sinsq2theta;
     break;
   }
   case(kNuMixSinsqTheta):
   {
     double mean  = conv::kGlobFitSinsqTheta12; 
     double sigma = conv::kGlobFitSinsqTheta12 * conv::kGlobFitErrSinsqTheta12;
     double sinsqtheta = Rnd::Instance()->RndToyMCosc().Gaus(mean,sigma);
     mix = sinsqtheta;
     break;
   }
   default:
     exit(1);
  }

  // return if a physical value was generated
  if(mix >=0. && mix <= 1.) {
     return mix;
  }

  // or try again
  return util::RandomizeMixing12();
}
//____________________________________________________________________________
double t2kana::util::RandomizeMixing23(void)
{

  double mix = 0;

  NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
  switch(omp) {
   case(kNuMixSinsq2Theta):
   {
     double mean  = conv::kGlobFitSinsq2Theta23; 
     double sigma = conv::kGlobFitSinsq2Theta23 * conv::kGlobFitErrSinsq2Theta23;
     double sinsq2theta = Rnd::Instance()->RndToyMCosc().Gaus(mean,sigma);
     mix = sinsq2theta;
     break;
   }
   case(kNuMixSinsqTheta):
   {
     double mean  = conv::kGlobFitSinsqTheta23; 
     double sigma = conv::kGlobFitSinsqTheta23 * conv::kGlobFitErrSinsqTheta23;
     double sinsqtheta = Rnd::Instance()->RndToyMCosc().Gaus(mean,sigma);
     mix = sinsqtheta;
     break;
   }
   default:
     exit(1);
  }

  // return if a physical value was generated
  if(mix >=0. && mix <= 1.) {
     return mix;
  }

  // or try again
  return util::RandomizeMixing23();
}
//____________________________________________________________________________
double t2kana::util::RandomizeMixing13(void)
{

  double mix = 0;

  NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
  switch(omp) {
  case(kNuMixSinsq2Theta):
    {
      double mean  = conv::kGlobFitSinsq2Theta13;
      double sigma = conv::kGlobFitSinsq2Theta13 * conv::kGlobFitErrSinsq2Theta13;
      double sinsq2theta = Rnd::Instance()->RndToyMCosc().Gaus(mean,sigma);
      mix = sinsq2theta;
      break;
    }
  case(kNuMixSinsqTheta):
    {
      double mean  = conv::kGlobFitSinsqTheta13;
      double sigma = conv::kGlobFitSinsqTheta13 * conv::kGlobFitErrSinsqTheta13;
      double sinsqtheta = Rnd::Instance()->RndToyMCosc().Gaus(mean,sigma);
      mix = sinsqtheta;
      break;
    }
  default:
    exit(1);
  }

  // return if a physical value was generated
  if(mix >=0. && mix <= 1.) {
    return mix;
  }

  // or try again 
  return util::RandomizeMixing13();
}
//____________________________________________________________________________
double t2kana::util::RandomizeCPV(void)
{
  NuOscDeltaCP_t cpv = Opt::Instance()->gOscDeltaCP;
  switch(cpv) {

    case (kNuDeltaCPPhase) :
    {
       double min = -conv::kpi;
       double max =  conv::kpi;
       return Rnd::Instance()->RndToyMCosc().Uniform(min,max);
    }
    break;

    case (kNuDeltaCPSin) :
    {
       double min = -1.;
       double max =  1.;
       return Rnd::Instance()->RndToyMCosc().Uniform(min,max);
    }
    break;

    default:
      return -99999;
  }
  return -99999;
}
//____________________________________________________________________________
double t2kana::util::RandomizeDm12sq(void)
{
  double mean  = conv::kGlobFitDm12sq; 
  double sigma = conv::kGlobFitDm12sq * conv::kGlobFitErrDm12sq;
  double dm    = Rnd::Instance()->RndToyMCosc().Gaus(mean,sigma);
  // return if a physical value was generated
  if(dm >= 0.) {
    return dm;
  }
  //or try again
  return util::RandomizeDm12sq();
}
//____________________________________________________________________________
double t2kana::util::RandomizeDm23sq(void)
{
  double mean  = 0.0;
  NuOscDmsqScheme_t odmsq = Opt::Instance()->gOscDmsqScheme;
  switch (odmsq) {
    
  case(kNuDmsqFogliLisi):
    mean  =  kGlobFitDmsqFL;
    break;
  case(kNuDmsq12andAbs23):
    mean  = kGlobFitDm23sq; 
    break;
  case(kNuDmsq12andAbs31):
    if(Opt::Instance()->gMassHierarchyAssumed == kNmhNormal)
      mean  = kGlobFitDm23sq + kGlobFitDm12sq; 
    else if(Opt::Instance()->gMassHierarchyAssumed == kNmhInverted)
      mean  = kGlobFitDm23sq - kGlobFitDm12sq; 
    break;
  default:
    mean  = kGlobFitDm23sq;  
  }       
  double sigma = mean * kGlobFitErrDm23sq;
  double dm    = Rnd::Instance()->RndToyMCosc().Gaus(mean,sigma);
  // return if a physical value was generated
  if(dm >= 0.) {
    return dm;
  }
  //or try again
  return util::RandomizeDm23sq();
}
//____________________________________________________________________________
int t2kana::util::Sign(double tweaking_dial)
{
  if(fabs(tweaking_dial) < kASmallNum) return 0;
  if (tweaking_dial < 0) return -1;
  return +1;
}
//____________________________________________________________________________
int t2kana::util::DiagonaliseJacobi(
     std::complex<double> A[3][3],  std::complex<double> Q[3][3], double w[3])
// ----------------------------------------------------------------------------
// Calculates the eigenvalues and normalized eigenvectors of a hermitian 3x3
// matrix A using the Jacobi algorithm.
// The upper triangular part of A is destroyed during the calculation,
// the diagonal elements are read but not destroyed, and the lower
// triangular elements are not referenced at all.
// ----------------------------------------------------------------------------
// Parameters:
//   A: The hermitian input matrix
//   Q: Storage buffer for eigenvectors
//   w: Storage buffer for eigenvalues
// ----------------------------------------------------------------------------
// Return value:
//   0: Success
//  -1: Error (no convergence)
// ----------------------------------------------------------------------------
// Adapted from zheevj3 by J.Kopp for C++ by Steve Dennis, 2013 under LGPL.
// ----------------------------------------------------------------------------
// Numerical diagonalization of 3x3 matrcies
// Copyright (C) 2006  Joachim Kopp
// ----------------------------------------------------------------------------
{
  const int n = 3;
  double sd, so;                  // Sums of diagonal resp. off-diagonal elements
  std::complex<double> s, t;      // sin(phi), tan(phi) and temporary storage
  double c;                       // cos(phi)
  double g, h, z;                 // Temporary storage
  double thresh;
  
  // Initialize Q to the identitity matrix
  for (int i=0; i < n; i++)
  {
    Q[i][i] = 1.0;
    for (int j=0; j < i; j++)
      Q[i][j] = Q[j][i] = 0.0;
  }

  // Initialize w to diag(A)
  for (int i=0; i < n; i++)
    w[i] = real(A[i][i]);

  // Calculate tr(A) squared  
  sd = 0.0;
  for (int i=0; i < n; i++)
    sd += fabs(w[i]);
  sd = sd*sd;
 
  // Main iteration loop
  for (int nIter=0; nIter < 50; nIter++)
  {
    // Test for convergence 
    so = 0.0;
    for (int p=0; p < n; p++)
      for (int q=p+1; q < n; q++)
        so += fabs(real(A[p][q])) + fabs(imag(A[p][q]));
    if (so == 0.0)
      return 0;

    if (nIter < 4)
      thresh = 0.2 * so / (n*n);
    else
      thresh = 0.0;

    // Do sweep
    for (int p=0; p < n; p++)
      for (int q=p+1; q < n; q++)
      {
        g = 100.0 * (fabs(real(A[p][q])) + fabs(imag(A[p][q])));
        if (nIter > 4  &&  fabs(w[p]) + g == fabs(w[p])
                       &&  fabs(w[q]) + g == fabs(w[q]))
        {
          A[p][q] = 0.0;
        }
        else if (fabs(real(A[p][q])) + fabs(imag(A[p][q])) > thresh)
        {
          // Calculate Jacobi transformation
          h = w[q] - w[p];
          if (fabs(h) + g == fabs(h))
            t = A[p][q] / h;
          else
          {
            if (h < 0.0)
              t = -2.0 * A[p][q] / (sqrt(h*h + 4.0*norm(A[p][q])) - h);
            else if (h == 0.0)
              t = A[p][q] * (1.0 / abs(A[p][q]));  // A[p][q]/fabs(A[p][q]) could cause overflows
            else
              t = 2.0 * A[p][q] / (sqrt(h*h + 4.0*norm(A[p][q])) + h);
          }
          c = 1.0/sqrt(1.0 + norm(t));
          s = t * c;
          z = real(t * conj(A[p][q]));

          // Apply Jacobi transformation
          A[p][q] = 0.0;
          w[p] -= z;
          w[q] += z;
          for (int r=0; r < p; r++)
          {
            t = A[r][p];
            A[r][p] = c*t - conj(s)*A[r][q];
            A[r][q] = s*t + c*A[r][q];
          }
          for (int r=p+1; r < q; r++)
          {
            t = A[p][r];
            A[p][r] = c*t - s*conj(A[r][q]);
            A[r][q] = s*conj(t) + c*A[r][q];
          }
          for (int r=q+1; r < n; r++)
          {
            t = A[p][r];
            A[p][r] = c*t - s*A[q][r];
            A[q][r] = conj(s)*t + c*A[q][r];
          }

          // Update eigenvectors      
          for (int r=0; r < n; r++)
          {
            t = Q[r][p];
            Q[r][p] = c*t - conj(s)*Q[r][q];
            Q[r][q] = s*t + c*Q[r][q];
          }
        }
      }
  }

  return -1;
}
//____________________________________________________________________________
int t2kana::util::DiagonaliseGSL(
#ifdef __GSL_ENABLED__ 
int ndim, double * matrix_in, double * eigenvectors_out, double * eigenvalues_out
#else
int, double *, double *, double *
#endif
)
{
#ifdef __GSL_ENABLED__ 
  // Initialise GSL memory structures
  gsl_matrix_complex * matr = gsl_matrix_complex_alloc(ndim,ndim);
  gsl_matrix_complex * eigenvectors = gsl_matrix_complex_alloc(ndim,ndim);
  gsl_vector * eigenvalues  = gsl_vector_alloc(ndim);
  gsl_eigen_hermv_workspace * workspace = gsl_eigen_hermv_alloc(ndim);
  
  // Put our matrix in the GSL matrix container
  memcpy( (double*)gsl_matrix_complex_ptr(matr,0,0), matrix_in, ndim*ndim*2*sizeof(double) );
  
  // Diagonalise
  int status=gsl_eigen_hermv(matr, eigenvalues, eigenvectors, workspace);
  
  // Copy the results to the output matrices
  memcpy( eigenvectors_out, (const double*)gsl_matrix_complex_const_ptr(eigenvectors,0,0) , ndim*ndim*2*sizeof(double));
  memcpy( eigenvalues_out , (const double*)gsl_vector_const_ptr        (eigenvalues ,0)   , ndim*sizeof(double));
  
  // Free GSL memory structures
  gsl_matrix_complex_free(matr);
  gsl_matrix_complex_free(eigenvectors);
  gsl_vector_free(eigenvalues);
  gsl_eigen_hermv_free(workspace);
  
  return (status != GSL_SUCCESS);
#else
  std::cerr<<"GSL Disabled"<<std::endl;
  exit(-1);
#endif
}
//____________________________________________________________________________
