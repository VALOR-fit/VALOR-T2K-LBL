//____________________________________________________________________________
/*!

\class     t2kana::utils

\brief     Numerical and statistics utilities

\author    Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
           University of Liverpool & STFC Rutherford Appleton Laboratory
          
\created   Oct 20, 2010

\copyright (c) 2010-2014, The VALOR group.

*/
//____________________________________________________________________________

#ifndef _t2kana_numstat_utils_h_
#define _t2kana_numstat_utils_h_

#include <string>
#include <vector>
#include <complex>

#include <TXMLEngine.h>
#include <TFile.h>
#include <TTree.h>
#include <TVectorD.h>
#include <TMatrixD.h>
#include <TH1D.h>
#include <TSpline.h>

class TGraph;
class TH1;
class TH1F;
class TH1D;
class TH2D;
class TAxis;

#include "Conventions.h"
#include "VALORArray.h"

using std::string;
using std::vector;

namespace t2kana {
namespace util {

  //
  // Statistics
  //

  double DeltaChisquare   (conv::ConfLevel_t cl, int nfitparams);
  double ChisqGauss       (double ndata, double nmc);
  double ChisqPoisson     (double ndata, double nmc);
  
  double ArrayChisqGauss  (const VALORArray *  data, const VALORArray *  mc, VALORArray * output = 0);
  double ArrayChisqPoisson(const VALORArray *  data, const VALORArray *  mc, VALORArray * output = 0);

  //
  // Cholesky decomposition
  //

  TMatrixD CholeskyDecomposition                     (TMatrixD& cov);                 ///< Cholesky decomposition. Return lower triangular matrix.
  TVectorD CholeskyGenerateCorrelatedParams          (TMatrixD& Lch, TVectorD& mean); ///< Generate a vector of correlated params
  TVectorD CholeskyGenerateCorrelatedParams          (TMatrixD& Lch, TVectorD& mean, TVectorD& g_uncorrelated); ///< Generate a vector of correlated params
  TVectorD CholeskyGenerateCorrelatedParamVariations (TMatrixD& Lch);                 ///< Generate a vector of correlated param variations
  TVectorD CholeskyCalculateCorrelatedParamVariations(TMatrixD& Lch, TVectorD& g);    ///<

  //
  // Functions randomizing oscillation params (for use when params are marginalized)
  //

  double RandomizeMixing23 (void); ///< pull from gaussian with mean & sigma from global fit (mixing param is either sinsq2theta or sinsqtheta depending on config)
  double RandomizeMixing12 (void); ///< pull from gaussian with mean & sigma from global fit ( -//- )
  double RandomizeMixing13 (void); ///< uniform from 0 up to ~ small_factor x CHOOZ limit    ( -//- )
  double RandomizeCPV      (void); ///< uniform from -pi to +pi (or -1 to +1)
  double RandomizeDm12sq   (void); ///< pull from gaussian with mean & sigma from global fit
  double RandomizeDm23sq   (void); ///< pull from gaussian with mean & sigma from global fit
  
  //
  // 3x3 matrix diagonalisation
  //

  int DiagonaliseJacobi(std::complex<double> A[3][3],  std::complex<double> Q[3][3], double w[3]);

  // Full diagonalisation using GSL
  // matrix_in is a pointer to a block representing a n*n complex matrix, equivalent to a double[n][n][2]
  // eigenvectors_out is pointer to (assigned) memory in same format
  // eigenvalues_out is pointer to double[n]
  int DiagonaliseGSL(int ndim, double * matrix_in, double * eigenvectors_out, double * eigenvalues_out);

  //
  // Etc
  //

  int Sign(double x);
  
} // util
} // t2kana
  
#endif


