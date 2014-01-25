//____________________________________________________________________________
/*

\copyright (c) 2010-2013, VaLOR group.

30-12-2013  Nick
Added samples and modes for reverse horn current running

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

#include "Utilities.h"
#include "Opt.h"
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
  if(nmc <= 0. || ndata <= 0.) 
    {
      if(2.*(nmc-ndata)<0){
	return -2.*(nmc-ndata);
      }
      return 2.*(nmc-ndata);
    }
  if ((2. * (nmc-ndata) + 2. * ndata * TMath::Log(ndata/nmc))<0)
    {
      if(2.*(nmc-ndata)<0){
	return -2.*(nmc-ndata);
      }
      return 2.*(nmc-ndata);
    }
  
  return (2. * (nmc-ndata) + 2. * ndata * TMath::Log(ndata/nmc));
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
double t2kana::util::RandomizeDeltaCP(void)
{
  double min = -conv::kpi;
  double max =  conv::kpi;
  return Rnd::Instance()->RndToyMCosc().Uniform(min,max);
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
  return util::RandomizeDm12sq();
}
//____________________________________________________________________________
void t2kana::util::SetDefaultSKSpectrumStyle (
   t2kana::conv::SuperkSample_t s, t2kana::conv::SuperkMode_t m, TH1D * hist)
{
  if(!hist) return;

  if(s == kSKsUndefined) return;
  if(m == kSKmUndefined) return;

  int sidx = util::ArrIdx::Instance()->SuperkSampleIdx(s);
  int midx = util::ArrIdx::Instance()->SuperkModeIdx(s,m);

  if(sidx == -1 || midx == -1) return;

  int fill_style  = conv::kSKHistoAttributes[sidx][midx][conv::kGrAttFillStyle];
  int fill_colour = conv::kSKHistoAttributes[sidx][midx][conv::kGrAttFillColour];
  int line_style  = conv::kSKHistoAttributes[sidx][midx][conv::kGrAttLineStyle];
  int line_width  = conv::kSKHistoAttributes[sidx][midx][conv::kGrAttLineWidth];
  int line_colour = conv::kSKHistoAttributes[sidx][midx][conv::kGrAttLineColour];

  if(fill_style != -1) {
    hist->SetFillStyle(fill_style);
  }
  if(fill_colour != -1) {
    hist->SetFillColor(fill_colour);
  }
  if(line_style != -1) {
    hist->SetLineStyle(line_style);
  }
  if(line_width != -1) {
    hist->SetLineWidth(line_width);
  }
  if(line_colour != -1) {
    hist->SetLineColor(line_colour);
  }
}
//____________________________________________________________________________
string t2kana::util::AsString (
  t2kana::conv::PhysHypothesis_t hyp, 
  double mix12, double mix13, double mix23, 
  double deltacp, double dm12sq, double dm23sq)
{
// The code builds strings which are typically used in building names of
// files containing t2kana results
//

  NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
  NuOscOctant_t   oct = Opt::Instance()->gOscOctantChoice;
  string mixparam = "?";
  switch(omp) {
      case(kNuMixSinsq2Theta):
	 switch(oct) {
	   case(kNuOctantOne):
	     mixparam = "sinsq2theta_octantone";
	     break;
	   case(kNuOctantTwo):
	     mixparam = "sinsq2theta_octanttwo";
	     break;
	   default:
	     return "";
	 }
         break;
      case(kNuMixSinsqTheta):
         mixparam = "sinsqtheta";
         break;
      default:
         return "";
  }

  ostringstream out;

  //
  // 3-flavour oscillations
  // Build the string as mix[12_nnn,13_nnn,23_nnn],dcp_nnn,dmsq[12_nnn,23_nnn]eV2
  // where mix is either sinsqtheta or sinsq2theta and nnn is the parameter
  // value with the specified numerical accuracy.
  // Params with value of -999 (randomized during toy-MC generation) will not be 
  // included in the output string
  //
  if( hyp == t2kana::conv::kHypOsc3F )
  {
    double rndmz = -900;
    bool show_mix12   = (mix12   > rndmz);
    bool show_mix13   = (mix13   > rndmz);
    bool show_mix23   = (mix23   > rndmz);
    bool show_mix     = show_mix12 || show_mix13 || show_mix23;
    bool show_deltacp = (deltacp > rndmz);
    bool show_dm12sq  = (dm12sq  > rndmz);
    bool show_dm23sq  = (dm23sq  > rndmz);
    bool show_dmsq    = show_dm12sq || show_dm23sq;

    if(show_mix) 
    {
       out << mixparam << "[";
       if(show_mix12) {
          out << "12_" << Form("%.3f",mix12);
       }
       if(show_mix13) {
          if(show_mix12) { 
             out << ","; 
          }
          out << "13_" << Form("%.3f",mix13);
       }
       if(show_mix23) {
          if(show_mix13 || (!show_mix13 && show_mix12)) {
             out << ",";
          }
          out << "23_" << Form("%.3f",mix23);
       }
       out << "]";
    }

    if(show_mix && (show_deltacp || show_dmsq)) 
    {
      out << ",";
    } 

    if(show_deltacp)
    {
      out << "dcp_" << Form("%.3f",deltacp);
    }

    if(show_deltacp && show_dmsq)
    {
      out << ",";
    } 

    if(show_dmsq) 
    {
      out << "dmsq_" << util::AsString(Opt::Instance()->gMassHierarchyAssumed) << "[";
       if(show_dm12sq) {
          out << "12_" << Form("%.8f",dm12sq);
       }
       if(show_dm23sq) {
          if(show_dm12sq) { 
             out << ","; 
          }
          out << "23_" << Form("%.6f",dm23sq);
       }
       out << "]eV2";
    }

  }//3f

  //
  // 2-flavour oscillations, numu disappearance
  // Build the string as mix23_nnn,dm23sq_nnneV2
  // where mix is either sinsqtheta or sinsq2theta and nnn is the parameter
  // value with the specified numerical accuracy.
  //
  if( hyp == t2kana::conv::kHypNuMuDisapp2F)
  {
    out << mixparam << "23_" << Form("%.3f",mix23) << ",dm23sq_" << Form("%.5f",dm23sq) << "eV2";
  } 

  return out.str();
}
//____________________________________________________________________________
string t2kana::util::AsString(t2kana::conv::OscAnaType_t a)
{
  switch(a) {
   case (kAnaCombined) :
     return "Combined nue appearance + numu disappearance analysis";
     break;
   case (kAnaNuMuDisappOnly) :
     return "Numu disappearance analysis only";
     break;
   case (kAnaNuEAppOnly) :
     return "Nue appearance analysis only";
     break;
   default :
     return "Undefined"; 
     break;
  }
  return "";
}
//____________________________________________________________________________
string t2kana::util::AsString(t2kana::conv::NuMassHierarchy_t m)
{
  switch(m) {
   case (kNmhNormal) :
     return "Normal";
     break;
   case (kNmhInverted) :
     return "Inverted";
     break;
   default :
     return "Undefined"; 
     break;
  }
  return "";
}
//____________________________________________________________________________
string t2kana::util::AsString(t2kana::conv::NuOscMixParam_t p)
{
  switch(p) {   
   case (kNuMixSinsq2Theta) :
     return "sin^2(2*theta_{ij})";
     break;
   case (kNuMixSinsqTheta) :
     return "sin^2(theta_{ij})";
     break;
   default :
     return "Undefined"; 
     break;
  }
  return "";
}
//____________________________________________________________________________
string t2kana::util::AsString(t2kana::conv::NuOscOctant_t o)
{
  switch(o) {   
   case (kNuOctantOne) :
     return "sin^2(theta_{23}) <= 0.5";
     break;
   case (kNuOctantTwo) :
     return "sin^2(theta_{23}) >= 0.5";
     break;
   default :
     return "Undefined"; 
     break;
  }
  return "";
}
//____________________________________________________________________________                                                                                                                            
string t2kana::util::AsString(t2kana::conv::NuOscDmsqScheme_t dmsqs)
{
  switch(dmsqs) {
  case(kNuDmsq12andAbs23) :
    return "dm^2 = Dm^2_{21}, Dm^2 = |Dm^2_{32}|";
    break;
  case(kNuDmsq12andAbs31) :
    return "dm^2 = Dm^2_{21}, Dm^2 = |Dm^2_{31}|";
    break;
  case(kNuDmsqFogliLisi ) :
    return "dm^2 = Dm^2_{21}, Dm^2 = |m^2_{3} - (m^2_{1} + m^2_{1})/2|";
    break;
  default :
    return "Undefined";
    break;
  }
}
//____________________________________________________________________________                                                                                                                            
string t2kana::util::AsString(t2kana::conv::NuOscDeltaCP_t deltaCP)
{
  switch(deltaCP) {
  case(kNuDeltaCPPhase) :
    return "d_{CP}";
    break;
  case(kNuDeltaCPSin) :
    return "sin(d_{CP})";
    break;
  default :
    return "Undefined";
    break;
  }
}
//____________________________________________________________________________
string t2kana::util::AsString(t2kana::conv::PhysHypothesis_t h)
{
  switch(h) {
    case (kHypDecay) : 
       return "<Neutrino Decay>"; 
       break;
    case (kHypDecoh) : 
       return "<Neutrino Decoherence>"; 
       break;
    case (kHypNuMuDisapp2F) : 
       return "<2-Flavour Neutrino Oscillations: nu_mu disappearance>"; 
       break;
    case (kHypOsc3F) : 
       return "<3-flavour Neutrino Oscillations>"; 
       break;
    default:                
      return "undefined"; 
      break;
  }
  return "";
}
//____________________________________________________________________________
string t2kana::util::AsString(PhysParamId_t p)
{
  switch(p) {

    case ( kPhypMixing23 ) : 
    {
        NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
        switch(omp) {
           case(kNuMixSinsq2Theta) : return "sin^2(2*theta_23)"; break;
           case(kNuMixSinsqTheta ) : return "sin^2(theta_23)";   break;
           default                 : return "--";                break;
         }
    }
    break;
    case ( kPhypMixing13 ) : 
    {
        NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
        switch(omp) {
           case(kNuMixSinsq2Theta) : return "sin^2(2*theta_13)";  break;
           case(kNuMixSinsqTheta ) : return "sin^2(theta_13)";    break;
           default                 : return "--";                 break;
         }
    }
    break;
    case ( kPhypMixing12 ) :      
    {
        NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
        switch(omp) {
           case(kNuMixSinsq2Theta) : return "sin^2(2*theta_12)";  break;
           case(kNuMixSinsqTheta ) : return "sin^2(theta_12)";    break;
           default                 : return "--";                 break;
         }
    }
    break;

    case ( kPhypDeltaCP               ) : return "delta_cp";                   break;
    case ( kPhypDm23sq                ) :
      {
	NuOscDmsqScheme_t dmscheme = Opt::Instance()->gOscDmsqScheme;
	switch(dmscheme) {
	  case(kNuDmsq12andAbs23) : return "|Dm^2_32|";                  break;
	  case(kNuDmsq12andAbs31) : return "|Dm^2_31|";                  break;
	  case(kNuDmsqFogliLisi ) : return "|Dm^2|";                     break;
  	  default                 : return "--";                         break;
	}
      } 
    break;

    case ( kPhypDm12sq                ) : return "Dm^2_12";                    break;

    case ( kPhypAlphaDecay            ) : return "alpha_decay";                break;
    case ( kPhypSinsq2ThetaDecay      ) : return "sin^2(2*theta_decay)";       break;
    case ( kPhypMusqDecoh             ) : return "mu^2_decoh";                 break;
    case ( kPhypSinsq2ThetaDecoh      ) : return "sin^2(2*theta_decoh)";       break;

    /*
     * 2011-analysis nuisance parameters
     */

    //***************************** FORWARD HORN CURRENT **********************************

    //
    // SK detector
    //
    case ( kNuisSuperkEnergyScale     ) : return "SK energy scale";                 break;    
#ifndef SKFSITOGETHER
#if FSI_BINS != 3
    case ( kNuisSuperkEff1RmuCCQEoth  ) : return "SK efficiency; 1Rmu; CCQE/Other"; break;
    case ( kNuisSuperkEff1RmuCCQErc0u ) : return "SK efficiency; 1Rmu; CCQE/RC0";   break;
    case ( kNuisSuperkEff1RmuCCQErc1u ) : return "SK efficiency; 1Rmu; CCQE/RC1";   break;
    case ( kNuisSuperkEff1RmuCCQErc2u ) : return "SK efficiency; 1Rmu; CCQE/RC2";   break;
    case ( kNuisSuperkEff1RmuCCnonQE  ) : return "SK efficiency; 1Rmu; CCnonQE";    break;
    case ( kNuisSuperkEff1RmuNC       ) : return "SK efficiency; 1Rmu; NC";         break;
    case ( kNuisSuperkEff1RmuNueCC    ) : return "SK efficiency; 1Rmu; nue CC";     break;
#endif
#if FSI_BINS != 2
    case ( kNuisSuperkEff1ReNueSigCC0    ) : return "SK efficiency; 1Re; osc. nue CC E0";     break;
    case ( kNuisSuperkEff1ReNueSigCC1    ) : return "SK efficiency; 1Re; osc. nue CC E1";     break;
    case ( kNuisSuperkEff1ReNueSigCC2    ) : return "SK efficiency; 1Re; osc. nue CC E2";     break;
    case ( kNuisSuperkEff1ReNumuCC0    )   : return "SK efficiency; 1Re; numu CC E0";     break;
    case ( kNuisSuperkEff1ReNumuCC1    )   : return "SK efficiency; 1Re; numu CC E1";     break;
    case ( kNuisSuperkEff1ReNumuCC2    )   : return "SK efficiency; 1Re; numu CC E2";     break;
    case ( kNuisSuperkEff1ReNueCC0     )   : return "SK efficiency; 1Re; nue CC E0";     break;
    case ( kNuisSuperkEff1ReNueCC1     )   : return "SK efficiency; 1Re; nue CC E1";     break;
    case ( kNuisSuperkEff1ReNueCC2     )   : return "SK efficiency; 1Re; nue CC E2";     break;
    case ( kNuisSuperkEff1ReNC0        )   : return "SK efficiency; 1Re; NC E0";     break;
    case ( kNuisSuperkEff1ReNC1        )   : return "SK efficiency; 1Re; NC E1";     break;
    case ( kNuisSuperkEff1ReNC2        )   : return "SK efficiency; 1Re; NC E2";     break;
#endif
#endif // SuperK Detector FSI separated
    //
    // BANFF flux and cross-section params
    //

    case ( kNuisBANFF00 ) :
    case ( kNuisBANFF01 ) :
    case ( kNuisBANFF02 ) :
    case ( kNuisBANFF03 ) :
    case ( kNuisBANFF04 ) :
    case ( kNuisBANFF05 ) :
    case ( kNuisBANFF06 ) :
    case ( kNuisBANFF07 ) :
    case ( kNuisBANFF08 ) :
    case ( kNuisBANFF09 ) :
    case ( kNuisBANFF10 ) :
    case ( kNuisBANFF11 ) :
    case ( kNuisBANFF12 ) :
    case ( kNuisBANFF13 ) :
    case ( kNuisBANFF14 ) :
    case ( kNuisBANFF15 ) :
#if (BANFF_BINS != 2) && (BANFF_BINS != 4)
    case ( kNuisBANFF16 ) :
    case ( kNuisBANFF17 ) :
    case ( kNuisBANFF18 ) :
    case ( kNuisBANFF19 ) :
    case ( kNuisBANFF20 ) :
    case ( kNuisBANFF21 ) :
#endif
#if BANFF_BINS == 1
    case ( kNuisBANFF22 ) :
    case ( kNuisBANFF23 ) :
    case ( kNuisBANFF24 ) :
#endif
    {
      int ip = (int)p - (int)kNuisBANFF00;
      return Form("BANFF %2d; SK %s flux; %.1f - %.1f GeV",
                   ip, util::Pdg2String(kNuPdgCodeOfFluxParamInBANFFfit[ip]).c_str(), 
                   kEtrueRangeOfFluxParamInBANFFfit[ip][0],
                   kEtrueRangeOfFluxParamInBANFFfit[ip][1]);     
      break;
    }
#ifndef ONEOVERMAQE
    case ( kNuisBANFFMAQE    ) : return "BANFF; Ma QE";                           break;
#else
    case ( kNuisBANFFMAQE    ) : return "BANFF; 1 / (Ma QE)^2";                   break;
#endif
    case ( kNuisBANFFMARES   ) : return "BANFF; Ma RES";                          break;
    case ( kNuisBANFFCCQEE0  ) : return "BANFF; Norm; CC QE; 0.0 - 1.5 GeV";      break;
#if BANFF_BINS != 3
    case ( kNuisBANFFCCQEE1  ) : return "BANFF; Norm; CC QE; 1.5 - 3.5 GeV";      break;
    case ( kNuisBANFFCCQEE2  ) : return "BANFF; Norm; CC QE; 3.5 - 30.0 GeV";     break;
#endif
    case ( kNuisBANFFCC1piE0 ) : return "BANFF; Norm; CC 1pi; 0.0 - 2.5 GeV";     break;
#if BANFF_BINS != 3
    case ( kNuisBANFFCC1piE1 ) : return "BANFF; Norm; CC 1pi; 2.5 - 30.0 GeV";    break;
#endif
#if (BANFF_BINS != 2) && (BANFF_BINS != 4)
    case ( kNuisBANFFNC1pi0  ) : return "BANFF; Norm; NC 1pi0";   break;
#endif
#if BANFF_BINS == 4
    case ( kNuisBANFFMECNorm ) : return "BANFF; Norm; MEC";                       break;
    case ( kNuisBANFFPDDDoubleCount ) : return "BANFF; PDD double counting";      break;
#endif

    //
    // XSec params not correlated with flux
    //

    case ( kNuisXSecFermiMom   ) : return "Fermi Momentum";                  break;
    case ( kNuisXSecSpecFunc   ) : return "Spectral Func";                   break; 
    case ( kNuisXSecWshape     ) : return "W-shape";                         break;
    case ( kNuisXSecShapeCCoth ) : return "Shape; CC Oth";                   break;
    case ( kNuisXSecPiLessDelta) : return "Pi-less delta decay";             break;
#ifndef __NUEONLY__
    case ( kNuisXSecBindE      ) : return "Binding energy";                  break;
#endif
    case ( kNuisXSecNormCCcoh  ) : return "Norm; CC Coh";                    break;
#ifndef __NUEONLY__
    case ( kNuisXSecNormNC1piC ) : return "Norm; NC 1piC";                   break;
#endif
#ifndef __NUMUONLY__
    case ( kNuisXSecNormNCcoh  ) : return "Norm; NC Coh";                    break;
#endif
    case ( kNuisXSecNormNCoth  ) : return "Norm; NC Oth";                    break;
    case ( kNuisXSecNormNuEToNuMu  ) : return "Norm; NuE To NuMu";           break;
#ifndef __NUEONLY__
    case ( kNuisXSecNormAntiNuToNu ) : return "Norm; AntiNu To Nu";          break;
#endif
    //
    // FSI
    //
#ifndef SKFSITOGETHER
#if FSI_BINS != 3
    case ( kNuisFSI000 ) :
    case ( kNuisFSI001 ) :
    case ( kNuisFSI002 ) :
    case ( kNuisFSI003 ) :
    case ( kNuisFSI004 ) :
    case ( kNuisFSI005 ) :
    {
      int ip = (int)p - (int)kNuisFSI000;
      string nutypes = "";
      switch(p) {
        case (kNuisFSI000) :
        case (kNuisFSI001) :
        case (kNuisFSI002) :
	  nutypes = "numu/numubar CCQE (1Rmu)"; break;
        case (kNuisFSI003) :
	  nutypes = "numu/numubar CCnQE (1Rmu)"; break;
        case (kNuisFSI004) :
	  nutypes = "nue/nuebar/signue CC (1Rmu)"; break;
        case (kNuisFSI005) :
	  nutypes = "all NC (1Rmu)"; break;
        default : nutypes = "--";
      }
      return Form("FSI/SI %2d; Ereco range %.1f - %.1f GeV; %s",ip,
		  kFSIParamsErecoCut[ip][0],kFSIParamsErecoCut[ip][1],nutypes.c_str());     
      break;
    }//FSI for 1Rmu
#endif
#if FSI_BINS != 2
    case ( kNuisFSI006 ) :
    case ( kNuisFSI007 ) :
    case ( kNuisFSI008 ) :
    case ( kNuisFSI009 ) :
    case ( kNuisFSI010 ) :
    case ( kNuisFSI011 ) :
    case ( kNuisFSI012 ) :
    case ( kNuisFSI013 ) :
    case ( kNuisFSI014 ) :
    case ( kNuisFSI015 ) :
    case ( kNuisFSI016 ) :
    case ( kNuisFSI017 ) :
      {   
#if FSI_BINS != 3
      int ip = (int)p - (int)kNuisFSI000;
#else
      int ip = (int)p - (int)kNuisFSI006;
#endif
      string nutypes = "";
	switch(p) {
        case (kNuisFSI006) :
        case (kNuisFSI007) :
        case (kNuisFSI008) :
	  nutypes = "oscillated nue CC (1Re)"; break;
        case (kNuisFSI009) :
        case (kNuisFSI010) :
        case (kNuisFSI011) :
	  nutypes = "numu/numubar CC (1Re)"; break;
        case (kNuisFSI012) :
        case (kNuisFSI013) :
        case (kNuisFSI014) :
	  nutypes = "nue/nuebar CC (1Re)"; break;
        case (kNuisFSI015) :
        case (kNuisFSI016) :
        case (kNuisFSI017) :
	  nutypes = "all NC (1Re)"; break;
        default : nutypes = "--";
	}
      return Form("FSI/SI %2d; Ereco range %.1f - %.1f GeV; %s",ip,
		  kFSIParamsErecoCut[ip][0],kFSIParamsErecoCut[ip][1],nutypes.c_str());     
      break;
      }
#endif
#endif // SKDet / FSI separated
    //
    // SKDET+FSI
    //
#ifdef SKFSITOGETHER
#ifndef SKCORRELATIONS
#if FSI_BINS != 3
    case ( kNuisSKDetFSI000 ) :
    case ( kNuisSKDetFSI001 ) :
    case ( kNuisSKDetFSI002 ) :
    case ( kNuisSKDetFSI003 ) :
    case ( kNuisSKDetFSI004 ) :
    case ( kNuisSKDetFSI005 ) :
    {
      int ip = (int)p - (int)kNuisSKDetFSI000;
      string nutypes = "";
      switch(p) {
        case (kNuisSKDetFSI000) :
        case (kNuisSKDetFSI001) :
        case (kNuisSKDetFSI002) :
	  nutypes = "numu/numubar CCQE (1Rmu)"; break;
        case (kNuisSKDetFSI003) :
	  nutypes = "numu/numubar CCnQE (1Rmu)"; break;
        case (kNuisSKDetFSI004) :
	  nutypes = "nue/nuebar/signue CC (1Rmu)"; break;
        case (kNuisSKDetFSI005) :
	  nutypes = "all NC (1Rmu)"; break;
        default : nutypes = "--";
      }
      return Form("SKDet + FSI/SI %2d; Ereco range %.2f - %.2f GeV; %s",ip,
		  kSKDetFSIParamsErecoCut[ip][0],kSKDetFSIParamsErecoCut[ip][1],nutypes.c_str());     
      break;
    }//FSI for 1Rmu
#endif
#if FSI_BINS != 2
    case ( kNuisSKDetFSI006 ) :
    case ( kNuisSKDetFSI007 ) :
    case ( kNuisSKDetFSI008 ) :
    case ( kNuisSKDetFSI009 ) :
    case ( kNuisSKDetFSI010 ) :
    case ( kNuisSKDetFSI011 ) :
    case ( kNuisSKDetFSI012 ) :
    case ( kNuisSKDetFSI013 ) :
    case ( kNuisSKDetFSI014 ) :
    case ( kNuisSKDetFSI015 ) :
    case ( kNuisSKDetFSI016 ) :
    case ( kNuisSKDetFSI017 ) :
      {   
#if FSI_BINS == 3
      int ip = (int)p - (int)kNuisSKDetFSI006;
#elif FSI_BINS ==1
      int ip = (int)p - (int)kNuisSKDetFSI000;
#endif
      string nutypes = "";
	switch(p) {
        case (kNuisSKDetFSI006) :
        case (kNuisSKDetFSI007) :
        case (kNuisSKDetFSI008) :
	  nutypes = "oscillated nue CC (1Re)"; break;
        case (kNuisSKDetFSI009) :
        case (kNuisSKDetFSI010) :
        case (kNuisSKDetFSI011) :
	  nutypes = "numu/numubar CC (1Re)"; break;
        case (kNuisSKDetFSI012) :
        case (kNuisSKDetFSI013) :
        case (kNuisSKDetFSI014) :
	  nutypes = "nue/nuebar CC (1Re)"; break;
        case (kNuisSKDetFSI015) :
        case (kNuisSKDetFSI016) :
        case (kNuisSKDetFSI017) :
	  nutypes = "all NC (1Re)"; break;
        default : nutypes = "--";
	}
      return Form("SKDet + FSI/SI %2d; Ereco range %.2f - %.2f GeV; %s",ip,
		  kSKDetFSIParamsErecoCut[ip][0],kSKDetFSIParamsErecoCut[ip][1],nutypes.c_str());     
      break;
      }
#endif
#else
    case ( kNuisSKDetFSI000 ) :
    case ( kNuisSKDetFSI001 ) :
    case ( kNuisSKDetFSI002 ) :
    case ( kNuisSKDetFSI003 ) :
    case ( kNuisSKDetFSI004 ) :
    case ( kNuisSKDetFSI005 ) :
    case ( kNuisSKDetFSI006 ) :
    case ( kNuisSKDetFSI007 ) :
    case ( kNuisSKDetFSI008 ) :
    case ( kNuisSKDetFSI009 ) :
    case ( kNuisSKDetFSI010 ) :
    case ( kNuisSKDetFSI011 ) :
    case ( kNuisSKDetFSI012 ) :
    case ( kNuisSKDetFSI013 ) :
    case ( kNuisSKDetFSI014 ) :
    case ( kNuisSKDetFSI015 ) :
    case ( kNuisSKDetFSI016 ) :
    case ( kNuisSKDetFSI017 ) :
    case ( kNuisSKDetFSI018 ) :
    case ( kNuisSKDetFSI019 ) :
    {
      int ip = (int)p - (int)kNuisSKDetFSI000;
      string nutypes = "";
      switch(p) {
        case (kNuisSKDetFSI000) :
        case (kNuisSKDetFSI001) :
        case (kNuisSKDetFSI002) :
        case (kNuisSKDetFSI003) :
        case (kNuisSKDetFSI004) :
	  nutypes = "signal nue CC (1Rmu)"; break;
        case (kNuisSKDetFSI005) :
        case (kNuisSKDetFSI006) :
        case (kNuisSKDetFSI007) :
        case (kNuisSKDetFSI008) :
        case (kNuisSKDetFSI009) :
	  nutypes = "numu/numubar CC (1Rmu)"; break;
        case (kNuisSKDetFSI010) :
        case (kNuisSKDetFSI011) :
        case (kNuisSKDetFSI012) :
        case (kNuisSKDetFSI013) :
        case (kNuisSKDetFSI014) :
	  nutypes = "nue/nuebar  CC (1Rmu)"; break;
        case (kNuisSKDetFSI015) :
        case (kNuisSKDetFSI016) :
        case (kNuisSKDetFSI017) :
        case (kNuisSKDetFSI018) :
        case (kNuisSKDetFSI019) :
	  nutypes = "all NC (1Rmu)"; break;
        default : nutypes = "--";
      }
      return Form("SKDet + FSI/SI %2d; Ereco range %.2f - %.2f GeV; %s",ip,
		  kSKDetFSIParamsErecoCut[ip][0],kSKDetFSIParamsErecoCut[ip][1],nutypes.c_str());     
      break;
    }//FSI for 1Rmu
    case ( kNuisSKDetFSI020 ) :
    case ( kNuisSKDetFSI021 ) :
    case ( kNuisSKDetFSI022 ) :
    case ( kNuisSKDetFSI023 ) :
    case ( kNuisSKDetFSI024 ) :
    case ( kNuisSKDetFSI025 ) :
    case ( kNuisSKDetFSI026 ) :
    case ( kNuisSKDetFSI027 ) :
    case ( kNuisSKDetFSI028 ) :
    case ( kNuisSKDetFSI029 ) :
    case ( kNuisSKDetFSI030 ) :
    case ( kNuisSKDetFSI031 ) :
      {   
      int ip = (int)p - (int)kNuisSKDetFSI000;
      string nutypes = "";
	switch(p) {
        case (kNuisSKDetFSI020) :
        case (kNuisSKDetFSI021) :
        case (kNuisSKDetFSI022) :
	  nutypes = "oscillated nue CC (1Re)"; break;
        case (kNuisSKDetFSI023) :
        case (kNuisSKDetFSI024) :
        case (kNuisSKDetFSI025) :
	  nutypes = "numu/numubar CC (1Re)"; break;
        case (kNuisSKDetFSI026) :
        case (kNuisSKDetFSI027) :
        case (kNuisSKDetFSI028) :
	  nutypes = "nue/nuebar CC (1Re)"; break;
        case (kNuisSKDetFSI029) :
        case (kNuisSKDetFSI030) :
        case (kNuisSKDetFSI031) :
	  nutypes = "all NC (1Re)"; break;
        default : nutypes = "--";
	}
      return Form("SKDet + FSI/SI %2d; Ereco range %.2f - %.2f GeV; %s",ip,
		  kSKDetFSIParamsErecoCut[ip][0],kSKDetFSIParamsErecoCut[ip][1],nutypes.c_str());     
      break;
      }
#endif
#endif //SKDETFSI  

  case ( kNuisRHCToFHC ) : return "RHC to FHC normalisation"; break;

      //************************** REVERSE HORN CURRENT *****************************

    case ( kRHCNuisSuperkEnergyScale     ) : return "SK energy scale; RHC";                 break;
  
      //
      // BANFF flux and cross-section params
      //

  case ( kRHCNuisBANFF00 ) :
  case ( kRHCNuisBANFF01 ) :
  case ( kRHCNuisBANFF02 ) :
  case ( kRHCNuisBANFF03 ) :
  case ( kRHCNuisBANFF04 ) :
  case ( kRHCNuisBANFF05 ) :
  case ( kRHCNuisBANFF06 ) :
  case ( kRHCNuisBANFF07 ) :
  case ( kRHCNuisBANFF08 ) :
  case ( kRHCNuisBANFF09 ) :
  case ( kRHCNuisBANFF10 ) :
  case ( kRHCNuisBANFF11 ) :
  case ( kRHCNuisBANFF12 ) :
  case ( kRHCNuisBANFF13 ) :
  case ( kRHCNuisBANFF14 ) :
  case ( kRHCNuisBANFF15 ) :
#if (BANFF_BINS != 2) && (BANFF_BINS != 4)
  case ( kRHCNuisBANFF16 ) :
  case ( kRHCNuisBANFF17 ) :
  case ( kRHCNuisBANFF18 ) :
  case ( kRHCNuisBANFF19 ) :
  case ( kRHCNuisBANFF20 ) :
  case ( kRHCNuisBANFF21 ) :
#endif
#if BANFF_BINS == 1
  case ( kRHCNuisBANFF22 ) :
  case ( kRHCNuisBANFF23 ) :
  case ( kRHCNuisBANFF24 ) :
#endif
    {
      int ip = (int)p - (int)kRHCNuisBANFF00;
      return Form("BANFF %2d; SK %s flux; %.1f - %.1f GeV; RHC",
		  ip, util::Pdg2String(kNuPdgCodeOfRHCFluxParamInBANFFfit[ip]).c_str(), 
		  kEtrueRangeOfRHCFluxParamInBANFFfit[ip][0],
		  kEtrueRangeOfRHCFluxParamInBANFFfit[ip][1]);     
      break;
    }
#ifndef ONEOVERMAQE
  case ( kRHCNuisBANFFMAQE    ) : return "BANFF; Ma QE; RHC";                           break;
#else
  case ( kRHCNuisBANFFMAQE    ) : return "BANFF; 1 / (Ma QE)^2; RHC";                   break;
#endif
  case ( kRHCNuisBANFFMARES   ) : return "BANFF; Ma RES; RHC";                          break;
  case ( kRHCNuisBANFFCCQEE0  ) : return "BANFF; Norm; CC QE; 0.0 - 1.5 GeV; RHC";      break;
#if BANFF_BINS != 3
  case ( kRHCNuisBANFFCCQEE1  ) : return "BANFF; Norm; CC QE; 1.5 - 3.5 GeV; RHC";      break;
  case ( kRHCNuisBANFFCCQEE2  ) : return "BANFF; Norm; CC QE; 3.5 - 30.0 GeV; RHC";     break;
#endif
  case ( kRHCNuisBANFFCC1piE0 ) : return "BANFF; Norm; CC 1pi; 0.0 - 2.5 GeV; RHC";     break;
#if BANFF_BINS != 3
  case ( kRHCNuisBANFFCC1piE1 ) : return "BANFF; Norm; CC 1pi; 2.5 - 30.0 GeV; RHC";    break;
#endif
#if (BANFF_BINS != 2) && (BANFF_BINS != 4)
  case ( kRHCNuisBANFFNC1pi0  ) : return "BANFF; Norm; NC 1pi0; RHC";   break;
#endif
#if BANFF_BINS == 4
  case ( kRHCNuisBANFFMECNorm ) : return "BANFF; Norm; MEC; RHC";                       break;
  case ( kRHCNuisBANFFPDDDoubleCount ) : return "BANFF; PDD double counting; RHC";      break;
#endif

    //
    // XSec params not correlated with flux
    //

  case ( kRHCNuisXSecFermiMom   ) : return "Fermi Momentum; RHC";                  break;
  case ( kRHCNuisXSecSpecFunc   ) : return "Spectral Func; RHC";                   break; 
  case ( kRHCNuisXSecWshape     ) : return "W-shape; RHC";                         break;
  case ( kRHCNuisXSecShapeCCoth ) : return "Shape; CC Oth; RHC";                   break;
  case ( kRHCNuisXSecPiLessDelta) : return "Pi-less delta decay; RHC";             break;
#ifndef __NUEONLY__
  case ( kRHCNuisXSecBindE      ) : return "Binding energy; RHC";                  break;
#endif
  case ( kRHCNuisXSecNormCCcoh  ) : return "Norm; CC Coh; RHC";                    break;
#ifndef __NUEONLY__
  case ( kRHCNuisXSecNormNC1piC ) : return "Norm; NC 1piC; RHC";                   break;
#endif
#ifndef __NUMUONLY__
  case ( kRHCNuisXSecNormNCcoh  ) : return "Norm; NC Coh; RHC";                    break;
#endif
  case ( kRHCNuisXSecNormNCoth  ) : return "Norm; NC Oth; RHC";                    break;
  case ( kRHCNuisXSecNormNuEToNuMu  ) : return "Norm; NuE To NuMu; RHC";           break;
#ifndef __NUEONLY__
  case ( kRHCNuisXSecNormAntiNuToNu ) : return "Norm; AntiNu To Nu; RHC";          break;
#endif

    //
    // SKDET+FSI
    //
#ifdef SKFSITOGETHER
#ifndef SKCORRELATIONS
#if FSI_BINS != 3
  case ( kRHCNuisSKDetFSI000 ) :
  case ( kRHCNuisSKDetFSI001 ) :
  case ( kRHCNuisSKDetFSI002 ) :
  case ( kRHCNuisSKDetFSI003 ) :
  case ( kRHCNuisSKDetFSI004 ) :
  case ( kRHCNuisSKDetFSI005 ) :
    {
      int ip = (int)p - (int)kRHCNuisSKDetFSI000;
      string nutypes = "";
      switch(p) {
      case (kRHCNuisSKDetFSI000) :
      case (kRHCNuisSKDetFSI001) :
      case (kRHCNuisSKDetFSI002) :
	nutypes = "numu/numubar CCQE (1Rmu)"; break;
      case (kRHCNuisSKDetFSI003) :
	nutypes = "numu/numubar CCnQE (1Rmu)"; break;
      case (kRHCNuisSKDetFSI004) :
	nutypes = "nue/nuebar/signue CC (1Rmu)"; break;
      case (kRHCNuisSKDetFSI005) :
	nutypes = "all NC (1Rmu)"; break;
      default : nutypes = "--";
      }
      return Form("SKDet + FSI/SI %2d; Ereco range %.2f - %.2f GeV; %s; RHC",ip,
                  kRHCSKDetFSIParamsErecoCut[ip][0],kRHCSKDetFSIParamsErecoCut[ip][1],nutypes.c_str());     
      break;
    }//FSI for 1Rmu
#endif
#if FSI_BINS != 2
  case ( kRHCNuisSKDetFSI006 ) :
  case ( kRHCNuisSKDetFSI007 ) :
  case ( kRHCNuisSKDetFSI008 ) :
  case ( kRHCNuisSKDetFSI009 ) :
  case ( kRHCNuisSKDetFSI010 ) :
  case ( kRHCNuisSKDetFSI011 ) :
  case ( kRHCNuisSKDetFSI012 ) :
  case ( kRHCNuisSKDetFSI013 ) :
  case ( kRHCNuisSKDetFSI014 ) :
  case ( kRHCNuisSKDetFSI015 ) :
  case ( kRHCNuisSKDetFSI016 ) :
  case ( kRHCNuisSKDetFSI017 ) :
    {
#if FSI_BINS == 3
      int ip = (int)p - (int)kRHCNuisSKDetFSI006;
#elif FSI_BINS ==1
      int ip = (int)p - (int)kRHCNuisSKDetFSI000;
#endif
      string nutypes = "";
      switch(p) {
      case (kRHCNuisSKDetFSI006) :
      case (kRHCNuisSKDetFSI007) :
      case (kRHCNuisSKDetFSI008) :
	nutypes = "oscillated nue CC (1Re)"; break;
      case (kRHCNuisSKDetFSI009) :
      case (kRHCNuisSKDetFSI010) :
      case (kRHCNuisSKDetFSI011) :
	nutypes = "numu/numubar CC (1Re)"; break;
      case (kRHCNuisSKDetFSI012) :
      case (kRHCNuisSKDetFSI013) :
      case (kRHCNuisSKDetFSI014) :
	nutypes = "nue/nuebar CC (1Re)"; break;
      case (kRHCNuisSKDetFSI015) :
      case (kRHCNuisSKDetFSI016) :
      case (kRHCNuisSKDetFSI017) :
	nutypes = "all NC (1Re)"; break;
      default : nutypes = "--";
      }
      return Form("SKDet + FSI/SI %2d; Ereco range %.2f - %.2f GeV; %s; RHC",ip,
                  kRHCSKDetFSIParamsErecoCut[ip][0],kRHCSKDetFSIParamsErecoCut[ip][1],nutypes.c_str());     
      break;
    }
#endif
#else
  case ( kRHCNuisSKDetFSI000 ) :
  case ( kRHCNuisSKDetFSI001 ) :
  case ( kRHCNuisSKDetFSI002 ) :
  case ( kRHCNuisSKDetFSI003 ) :
  case ( kRHCNuisSKDetFSI004 ) :
  case ( kRHCNuisSKDetFSI005 ) :
  case ( kRHCNuisSKDetFSI006 ) :
  case ( kRHCNuisSKDetFSI007 ) :
  case ( kRHCNuisSKDetFSI008 ) :
  case ( kRHCNuisSKDetFSI009 ) :
  case ( kRHCNuisSKDetFSI010 ) :
  case ( kRHCNuisSKDetFSI011 ) :
  case ( kRHCNuisSKDetFSI012 ) :
  case ( kRHCNuisSKDetFSI013 ) :
  case ( kRHCNuisSKDetFSI014 ) :
  case ( kRHCNuisSKDetFSI015 ) :
  case ( kRHCNuisSKDetFSI016 ) :
  case ( kRHCNuisSKDetFSI017 ) :
  case ( kRHCNuisSKDetFSI018 ) :
  case ( kRHCNuisSKDetFSI019 ) :
    {
      {
	int ip = (int)p - (int)kRHCNuisSKDetFSI000;
	string nutypes = "";
	switch(p) {
        case (kRHCNuisSKDetFSI000) :
        case (kRHCNuisSKDetFSI001) :
        case (kRHCNuisSKDetFSI002) :
        case (kRHCNuisSKDetFSI003) :
        case (kRHCNuisSKDetFSI004) :
          nutypes = "signal nue CC (1Rmu)"; break;
        case (kRHCNuisSKDetFSI005) :
        case (kRHCNuisSKDetFSI006) :
        case (kRHCNuisSKDetFSI007) :
        case (kRHCNuisSKDetFSI008) :
        case (kRHCNuisSKDetFSI009) :
          nutypes = "numu/numubar CC (1Rmu)"; break;
        case (kRHCNuisSKDetFSI010) :
        case (kRHCNuisSKDetFSI011) :
        case (kRHCNuisSKDetFSI012) :
        case (kRHCNuisSKDetFSI013) :
        case (kRHCNuisSKDetFSI014) :
          nutypes = "nue/nuebar  CC (1Rmu)"; break;
        case (kRHCNuisSKDetFSI015) :
        case (kRHCNuisSKDetFSI016) :
        case (kRHCNuisSKDetFSI017) :
        case (kRHCNuisSKDetFSI018) :
        case (kRHCNuisSKDetFSI019) :
          nutypes = "all NC (1Rmu)"; break;
        default : nutypes = "--";
	}
	return Form("SKDet + FSI/SI %2d; Ereco range %.2f - %.2f GeV; %s; RHC",ip,
		    kRHCSKDetFSIParamsErecoCut[ip][0],kRHCSKDetFSIParamsErecoCut[ip][1],nutypes.c_str());     
	break;
      }//FSI for 1Rmu
    case ( kRHCNuisSKDetFSI020 ) :
    case ( kRHCNuisSKDetFSI021 ) :
    case ( kRHCNuisSKDetFSI022 ) :
    case ( kRHCNuisSKDetFSI023 ) :
    case ( kRHCNuisSKDetFSI024 ) :
    case ( kRHCNuisSKDetFSI025 ) :
    case ( kRHCNuisSKDetFSI026 ) :
    case ( kRHCNuisSKDetFSI027 ) :
    case ( kRHCNuisSKDetFSI028 ) :
    case ( kRHCNuisSKDetFSI029 ) :
    case ( kRHCNuisSKDetFSI030 ) :
    case ( kRHCNuisSKDetFSI031 ) :
      { 
	int ip = (int)p - (int)kRHCNuisSKDetFSI000;
	string nutypes = "";
        switch(p) {
        case (kRHCNuisSKDetFSI020) :
        case (kRHCNuisSKDetFSI021) :
        case (kRHCNuisSKDetFSI022) :
          nutypes = "oscillated nue CC (1Re)"; break;
        case (kRHCNuisSKDetFSI023) :
        case (kRHCNuisSKDetFSI024) :
        case (kRHCNuisSKDetFSI025) :
          nutypes = "numu/numubar CC (1Re)"; break;
        case (kRHCNuisSKDetFSI026) :
        case (kRHCNuisSKDetFSI027) :
        case (kRHCNuisSKDetFSI028) :
          nutypes = "nue/nuebar CC (1Re)"; break;
        case (kRHCNuisSKDetFSI029) :
        case (kRHCNuisSKDetFSI030) :
        case (kRHCNuisSKDetFSI031) :
          nutypes = "all NC (1Re)"; break;
        default : nutypes = "--";
        }
	return Form("SKDet + FSI/SI %2d; Ereco range %.2f - %.2f GeV; %s; RHC",ip,
		    kRHCSKDetFSIParamsErecoCut[ip][0],kRHCSKDetFSIParamsErecoCut[ip][1],nutypes.c_str());     
	break;
      }
#endif
#endif //SKDETFSI

    default : 
      return "--";
  }
  return "--";
}
//____________________________________________________________________________
string t2kana::util::AsStringTex(PhysParamId_t p)
{
  switch(p) {

    case ( kPhypMixing23 ) : 
    {
        NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
        switch(omp) {
           case(kNuMixSinsq2Theta) : return "$\\sin^2(2\\theta_{23})$"; break;
           case(kNuMixSinsqTheta ) : return "$\\sin^2(\\theta_{23})$";   break;
           default                 : return "--";                break;
         }
    }
    break;
    case ( kPhypMixing13 ) : 
    {
        NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
        switch(omp) {
           case(kNuMixSinsq2Theta) : return "$\\sin^2(2\\theta_{13})$";  break;
           case(kNuMixSinsqTheta ) : return "$\\sin^2(\\theta_{13})$";    break;
           default                 : return "--";                 break;
         }
    }
    break;
    case ( kPhypMixing12 ) :      
    {
        NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
        switch(omp) {
           case(kNuMixSinsq2Theta) : return "$\\sin^2(2\\theta_{12})$";  break;
           case(kNuMixSinsqTheta ) : return "$\\sin^2(\\theta_{12})$";    break;
           default                 : return "--";                 break;
         }
    }
    break;

    case ( kPhypDeltaCP               ) : return "$\\delta_{CP}$";                   break;
    case ( kPhypDm12sq                ) : return "$|\\Delta m^2_{21}|$";                  break;
    case ( kPhypDm23sq                ) :
      {
	NuOscDmsqScheme_t dmscheme = Opt::Instance()->gOscDmsqScheme;
	switch(dmscheme) {
	  case(kNuDmsq12andAbs23) : return "|\\Delta m^2_{32}|";                  break;
	  case(kNuDmsq12andAbs31) : return "|\\Delta m^2_{31}|";                  break;
	  case(kNuDmsqFogliLisi ) : return "|\\Delta m^2|";                       break;
  	  default                 : return "--";                                  break;
	}
      }
    break;

    case ( kPhypAlphaDecay            ) : return "$\\alpha_{decay}$";                break;
    case ( kPhypSinsq2ThetaDecay      ) : return "$\\sin^2(2\\theta_{decay})$";       break;
    case ( kPhypMusqDecoh             ) : return "$\\mu^2_{decoh}$";                 break;
    case ( kPhypSinsq2ThetaDecoh      ) : return "$\\sin^2(2\\theta_{decoh})$";       break;

    /*
     * 2011-analysis nuisance parameters
     */

    //
    // SK detector
    //
    case ( kNuisSuperkEnergyScale     ) : return "SK energy scale";                 break;    
#ifndef SKFSITOGETHER
#if FSI_BINS != 3
    case ( kNuisSuperkEff1RmuCCQEoth  ) : return "SK efficiency, 1Rmu, CCQE//Other"; break;
    case ( kNuisSuperkEff1RmuCCQErc0u ) : return "SK efficiency, 1Rmu, CCQE//RC0";   break;
    case ( kNuisSuperkEff1RmuCCQErc1u ) : return "SK efficiency, 1Rmu, CCQE//RC1";   break;
    case ( kNuisSuperkEff1RmuCCQErc2u ) : return "SK efficiency, 1Rmu, CCQE//RC2";   break;
    case ( kNuisSuperkEff1RmuCCnonQE  ) : return "SK efficiency, 1Rmu, CCnonQE";    break;
    case ( kNuisSuperkEff1RmuNC       ) : return "SK efficiency, 1Rmu, NC";         break;
    case ( kNuisSuperkEff1RmuNueCC    ) : return "SK efficiency, 1Rmu, nue CC";     break;
#endif
#if FSI_BINS != 2
    case ( kNuisSuperkEff1ReNueSigCC0 ) : return "SK efficiency, 1Re, Osc. $\\nu_{e}$ CC E0"; break;
    case ( kNuisSuperkEff1ReNueSigCC1 ) : return "SK efficiency, 1Re, Osc. $\\nu_{e}$ CC E1"; break;
    case ( kNuisSuperkEff1ReNueSigCC2 ) : return "SK efficiency, 1Re, Osc. $\\nu_{e}$ CC E2"; break;
    case ( kNuisSuperkEff1ReNumuCC0 )   : return "SK efficiency, 1Re, Osc. $\\nu_{\\mu}$ CC E0"; break;
    case ( kNuisSuperkEff1ReNumuCC1 )   : return "SK efficiency, 1Re, Osc. $\\nu_{\\mu}$ CC E1"; break;
    case ( kNuisSuperkEff1ReNumuCC2 )   : return "SK efficiency, 1Re, Osc. $\\nu_{\\mu}$ CC E2"; break;
    case ( kNuisSuperkEff1ReNueCC0 ) : return "SK efficiency, 1Re, $\\nu_{e}$ CC E0"; break;
    case ( kNuisSuperkEff1ReNueCC1 ) : return "SK efficiency, 1Re, $\\nu_{e}$ CC E1"; break;
    case ( kNuisSuperkEff1ReNueCC2 ) : return "SK efficiency, 1Re, $\\nu_{e}$ CC E2"; break;
    case ( kNuisSuperkEff1ReNC0 )    : return "SK efficiency, 1Re, NC E0"; break;
    case ( kNuisSuperkEff1ReNC1 )    : return "SK efficiency, 1Re, NC E1"; break;
    case ( kNuisSuperkEff1ReNC2 )    : return "SK efficiency, 1Re, NC E2"; break;
#endif
#endif
    //
    // BANFF flux and cross-section params
    //
    case ( kNuisBANFF00 ) :
    case ( kNuisBANFF01 ) :
    case ( kNuisBANFF02 ) :
    case ( kNuisBANFF03 ) :
    case ( kNuisBANFF04 ) :
    case ( kNuisBANFF05 ) :
    case ( kNuisBANFF06 ) :
    case ( kNuisBANFF07 ) :
    case ( kNuisBANFF08 ) :
    case ( kNuisBANFF09 ) :
    case ( kNuisBANFF10 ) :
    case ( kNuisBANFF11 ) :
    case ( kNuisBANFF12 ) :
    case ( kNuisBANFF13 ) :
    case ( kNuisBANFF14 ) :
    case ( kNuisBANFF15 ) :
#if (BANFF_BINS != 2) && (BANFF_BINS != 4)
    case ( kNuisBANFF16 ) :
    case ( kNuisBANFF17 ) :
    case ( kNuisBANFF18 ) :
    case ( kNuisBANFF19 ) :
    case ( kNuisBANFF20 ) :
    case ( kNuisBANFF21 ) :
#endif
#if BANFF_BINS == 1
    case ( kNuisBANFF22 ) :
    case ( kNuisBANFF23 ) :
    case ( kNuisBANFF24 ) :
#endif
    {
      int ip = (int)p - (int)kNuisBANFF00;
      return Form("BANFF %2d, SK %s flux, %.1f - %.1f GeV",
                   ip, util::Pdg2String(kNuPdgCodeOfFluxParamInBANFFfit[ip]).c_str(), 
                   kEtrueRangeOfFluxParamInBANFFfit[ip][0],
                   kEtrueRangeOfFluxParamInBANFFfit[ip][1]);     
      break;
    }
#ifndef ONEOVERMAQE
    case ( kNuisBANFFMAQE    ) : return "BANFF, Ma QE";                           break;
#else
    case ( kNuisBANFFMAQE    ) : return "BANFF, 1 / (Ma QE)$^2$";                 break;
#endif
    case ( kNuisBANFFMARES   ) : return "BANFF, Ma RES";                          break;
    case ( kNuisBANFFCCQEE0  ) : return "BANFF; Norm; CC QE; 0.0 - 1.5 GeV";      break;
#if BANFF_BINS != 3
    case ( kNuisBANFFCCQEE1  ) : return "BANFF; Norm; CC QE; 1.5 - 3.5 GeV";      break;
    case ( kNuisBANFFCCQEE2  ) : return "BANFF; Norm; CC QE; 3.5 - 30.0 GeV";     break;
#endif
    case ( kNuisBANFFCC1piE0 ) : return "BANFF; Norm; CC 1pi; 0.0 - 2.5 GeV";     break;
#if BANFF_BINS != 3
    case ( kNuisBANFFCC1piE1 ) : return "BANFF; Norm; CC 1pi; 2.5 - 30.0 GeV";    break;
#endif
#if (BANFF_BINS != 2) && (BANFF_BINS != 4)
    case ( kNuisBANFFNC1pi0  ) : return "BANFF; Norm; NC 1pi0";    break;
#endif
#if BANFF_BINS == 4
    case ( kNuisBANFFMECNorm ) : return "BANFF; Norm; MEC";                       break;
    case ( kNuisBANFFPDDDoubleCount ) : return "BANFF; PDD double counting";      break;
#endif

    //
    // XSec params not correlated with flux
    //

    case ( kNuisXSecFermiMom   ) : return "Fermi Momentum";                  break;
    case ( kNuisXSecSpecFunc   ) : return "Spectral Func";                   break; 
    case ( kNuisXSecWshape     ) : return "W-shape";                         break;
    case ( kNuisXSecShapeCCoth ) : return "Shape; CC Oth";                   break;
    case ( kNuisXSecPiLessDelta) : return "Pi-less delta decay";             break;
#ifndef __NUEONLY__
    case ( kNuisXSecBindE      ) : return "Binding energy";                  break;
#endif
    case ( kNuisXSecNormCCcoh  ) : return "Norm, CC Coh";                    break;
#ifndef __NUEONLY__
    case ( kNuisXSecNormNC1piC ) : return "Norm, NC 1piC";                   break;
#endif
#ifndef __NUMUONLY__
    case ( kNuisXSecNormNCcoh  ) : return "Norm, NC Coh";                    break;
#endif
    case ( kNuisXSecNormNCoth  ) : return "Norm, NC Oth";                    break;
    case ( kNuisXSecNormNuEToNuMu  ) : return "Norm, $\\nu_e$ To $\\nu_\\mu$";  break;
#ifndef __NUEONLY__
    case ( kNuisXSecNormAntiNuToNu ) : return "Norm, $\\bar{\\nu}$ To $\\nu$";  break;
#endif
#ifndef SKFSITOGETHER
#if FSI_BINS != 3
    case ( kNuisFSI000 ) :
    case ( kNuisFSI001 ) :
    case ( kNuisFSI002 ) :
    case ( kNuisFSI003 ) :
    case ( kNuisFSI004 ) :
    case ( kNuisFSI005 ) :
    {
      string nutypes = "";
      int ip = (int)p - (int)kNuisFSI000;
      switch(p) {
        case (kNuisFSI000) :
        case (kNuisFSI001) :
        case (kNuisFSI002) :
	  nutypes = "$\\nu_\\mu$/$\\bar{\\nu_\\mu}$ CCQE ($1R_{\\mu}$)"; break;
        case (kNuisFSI003) :
	  nutypes = "$\\nu_\\mu$/$\\bar{\\nu_\\mu}$ CCnQE ($1R_{\\mu}$)"; break;
        case (kNuisFSI004) :
	  nutypes = "$\\nu_e$/$\\bar{\\nu_e}$/sig$\\nu_e$ CC ($1R_{\\mu}$)"; break;
        case (kNuisFSI005) :
	  nutypes = "all NC ($1R_{\\mu}$)"; break;
        default : nutypes = "--";
      }
      return Form("FSI/SI %2d; Ereco range %.1f - %.1f GeV; %s",ip,
		  kFSIParamsErecoCut[ip][0],kFSIParamsErecoCut[ip][1],nutypes.c_str());     
      break;
    }//FSI for 1Rmu
#endif
#if FSI_BINS != 2
    case ( kNuisFSI006 ) :
    case ( kNuisFSI007 ) :
    case ( kNuisFSI008 ) :
    case ( kNuisFSI009 ) :
    case ( kNuisFSI010 ) :
    case ( kNuisFSI011 ) :
    case ( kNuisFSI012 ) :
    case ( kNuisFSI013 ) :
    case ( kNuisFSI014 ) :
    case ( kNuisFSI015 ) :
    case ( kNuisFSI016 ) :
    case ( kNuisFSI017 ) :

    {
      string nutypes = "";
#if FSI_BINS != 3
      int ip = (int)p - (int)kNuisFSI000;
#else
      int ip = (int)p - (int)kNuisFSI006;
#endif
      switch(p) {
        case (kNuisFSI006) :
        case (kNuisFSI007) :
        case (kNuisFSI008) :
	  nutypes = "Osc. $\\nu_e$ CC ($1R_{e}$)"; break;
        case (kNuisFSI009) :
        case (kNuisFSI010) :
        case (kNuisFSI011) :
	  nutypes = "$\\nu_\\mu$/$\\bar{\\nu_\\mu}$ CC ($1R_{e}$)"; break;
        case (kNuisFSI012) :
        case (kNuisFSI013) :
        case (kNuisFSI014) :
	  nutypes = "$\\nu_e$/$\\bar{\\nu_e}$ CC ($1R_{e}$)"; break;
        case (kNuisFSI015) :
        case (kNuisFSI016) :
        case (kNuisFSI017) :
	  nutypes = "all NC ($1R_{e}$)"; break;
        default : nutypes = "--";
      }
      return Form("FSI/SI %2d; Ereco range %.1f - %.1f GeV; %s",ip,
		  kFSIParamsErecoCut[ip][0],kFSIParamsErecoCut[ip][1],nutypes.c_str());     
      break;
    }//FSI for 1Re
#endif
#endif
#ifdef SKFSITOGETHER
#ifndef SKCORRELATIONS
#if FSI_BINS != 3
    case ( kNuisSKDetFSI000 ) :
    case ( kNuisSKDetFSI001 ) :
    case ( kNuisSKDetFSI002 ) :
    case ( kNuisSKDetFSI003 ) :
    case ( kNuisSKDetFSI004 ) :
    case ( kNuisSKDetFSI005 ) :
    {
      string nutypes = "";
      int ip = (int)p - (int)kNuisSKDetFSI000;
      switch(p) {
        case (kNuisSKDetFSI000) :
        case (kNuisSKDetFSI001) :
        case (kNuisSKDetFSI002) :
	  nutypes = "$\\nu_\\mu$/$\\bar{\\nu_\\mu}$ CCQE ($1R_{\\mu}$)"; break;
        case (kNuisSKDetFSI003) :
	  nutypes = "$\\nu_\\mu$/$\\bar{\\nu_\\mu}$ CCnQE ($1R_{\\mu}$)"; break;
        case (kNuisSKDetFSI004) :
	  nutypes = "$\\nu_e$/$\\bar{\\nu_e}$/sig$\\nu_e$ CC ($1R_{\\mu}$)"; break;
        case (kNuisSKDetFSI005) :
	  nutypes = "all NC ($1R_{\\mu}$)"; break;
        default : nutypes = "--";
      }
      return Form("SKDet + FSI/SI %2d; Ereco range %.2f - %.2f GeV; %s",ip,
		  kSKDetFSIParamsErecoCut[ip][0],kSKDetFSIParamsErecoCut[ip][1],nutypes.c_str());     
      break;
    }//FSI for 1Rmu
#endif
#if FSI_BINS != 2
    case ( kNuisSKDetFSI006 ) :
    case ( kNuisSKDetFSI007 ) :
    case ( kNuisSKDetFSI008 ) :
    case ( kNuisSKDetFSI009 ) :
    case ( kNuisSKDetFSI010 ) :
    case ( kNuisSKDetFSI011 ) :
    case ( kNuisSKDetFSI012 ) :
    case ( kNuisSKDetFSI013 ) :
    case ( kNuisSKDetFSI014 ) :
    case ( kNuisSKDetFSI015 ) :
    case ( kNuisSKDetFSI016 ) :
    case ( kNuisSKDetFSI017 ) :

    {
      string nutypes = "";
#if FSI_BINS != 3
      int ip = (int)p - (int)kNuisSKDetFSI000;
#else
      int ip = (int)p - (int)kNuisSKDetFSI006;
#endif
      switch(p) {
        case (kNuisSKDetFSI006) :
        case (kNuisSKDetFSI007) :
        case (kNuisSKDetFSI008) :
	  nutypes = "Osc. $\\nu_e$ CC ($1R_{e}$)"; break;
        case (kNuisSKDetFSI009) :
        case (kNuisSKDetFSI010) :
        case (kNuisSKDetFSI011) :
	  nutypes = "$\\nu_\\mu$/$\\bar{\\nu_\\mu}$ CC ($1R_{e}$)"; break;
        case (kNuisSKDetFSI012) :
        case (kNuisSKDetFSI013) :
        case (kNuisSKDetFSI014) :
	  nutypes = "$\\nu_e$/$\\bar{\\nu_e}$ CC ($1R_{e}$)"; break;
        case (kNuisSKDetFSI015) :
        case (kNuisSKDetFSI016) :
        case (kNuisSKDetFSI017) :
	  nutypes = "all NC ($1R_{e}$)"; break;
        default : nutypes = "--";
      }
      return Form("SKDet + FSI/SI %2d; Ereco range %.2f - %.2f GeV; %s",ip,
		  kSKDetFSIParamsErecoCut[ip][0],kSKDetFSIParamsErecoCut[ip][1],nutypes.c_str());     
      break;
    }//FSI for 1Re
#endif
#else
    case ( kNuisSKDetFSI000 ) :
    case ( kNuisSKDetFSI001 ) :
    case ( kNuisSKDetFSI002 ) :
    case ( kNuisSKDetFSI003 ) :
    case ( kNuisSKDetFSI004 ) :
    case ( kNuisSKDetFSI005 ) :
    case ( kNuisSKDetFSI006 ) :
    case ( kNuisSKDetFSI007 ) :
    case ( kNuisSKDetFSI008 ) :
    case ( kNuisSKDetFSI009 ) :
    case ( kNuisSKDetFSI010 ) :
    case ( kNuisSKDetFSI011 ) :
    case ( kNuisSKDetFSI012 ) :
    case ( kNuisSKDetFSI013 ) :
    case ( kNuisSKDetFSI014 ) :
    case ( kNuisSKDetFSI015 ) :
    case ( kNuisSKDetFSI016 ) :
    case ( kNuisSKDetFSI017 ) :
    case ( kNuisSKDetFSI018 ) :
    case ( kNuisSKDetFSI019 ) :
    {
      int ip = (int)p - (int)kNuisSKDetFSI000;
      string nutypes = "";
      switch(p) {
        case (kNuisSKDetFSI000) :
        case (kNuisSKDetFSI001) :
        case (kNuisSKDetFSI002) :
        case (kNuisSKDetFSI003) :
        case (kNuisSKDetFSI004) :
	  nutypes = "Osc. $\\nu_e$ CC ($1R_{\\mu}$)"; break;
        case (kNuisSKDetFSI005) :
        case (kNuisSKDetFSI006) :
        case (kNuisSKDetFSI007) :
        case (kNuisSKDetFSI008) :
        case (kNuisSKDetFSI009) :
	  nutypes = "$\\nu_\\mu$/$\\bar{\\nu_\\mu}$ CC ($1R_{\\mu}$)"; break;
        case (kNuisSKDetFSI010) :
        case (kNuisSKDetFSI011) :
        case (kNuisSKDetFSI012) :
        case (kNuisSKDetFSI013) :
        case (kNuisSKDetFSI014) :
	  nutypes = "$\\nu_e$/$\\bar{\\nu_e}$ CC ($1R_{\\mu}$)"; break;
        case (kNuisSKDetFSI015) :
        case (kNuisSKDetFSI016) :
        case (kNuisSKDetFSI017) :
        case (kNuisSKDetFSI018) :
        case (kNuisSKDetFSI019) :
	  nutypes = "all NC ($1R_{\\mu}$)"; break;
        default : nutypes = "--";
      }
      return Form("SKDet + FSI/SI %2d; Ereco range %.2f - %.2f GeV; %s",ip,
		  kSKDetFSIParamsErecoCut[ip][0],kSKDetFSIParamsErecoCut[ip][1],nutypes.c_str());     
      break;
    }//FSI for 1Rmu
    case ( kNuisSKDetFSI020 ) :
    case ( kNuisSKDetFSI021 ) :
    case ( kNuisSKDetFSI022 ) :
    case ( kNuisSKDetFSI023 ) :
    case ( kNuisSKDetFSI024 ) :
    case ( kNuisSKDetFSI025 ) :
    case ( kNuisSKDetFSI026 ) :
    case ( kNuisSKDetFSI027 ) :
    case ( kNuisSKDetFSI028 ) :
    case ( kNuisSKDetFSI029 ) :
    case ( kNuisSKDetFSI030 ) :
    case ( kNuisSKDetFSI031 ) :
      {   
      int ip = (int)p - (int)kNuisSKDetFSI000;
      string nutypes = "";
	switch(p) {
        case (kNuisSKDetFSI020) :
        case (kNuisSKDetFSI021) :
        case (kNuisSKDetFSI022) :
	  nutypes = "Osc. $\\nu_e$ CC ($1R_{e}$)"; break;
        case (kNuisSKDetFSI023) :
        case (kNuisSKDetFSI024) :
        case (kNuisSKDetFSI025) :
	  nutypes = "$\\nu_\\mu$/$\\bar{\\nu_\\mu}$ CC ($1R_{e}$)"; break;
        case (kNuisSKDetFSI026) :
        case (kNuisSKDetFSI027) :
        case (kNuisSKDetFSI028) :
	  nutypes = "$\\nu_e$/$\\bar{\\nu_e}$ CC ($1R_{e}$)"; break;
        case (kNuisSKDetFSI029) :
        case (kNuisSKDetFSI030) :
        case (kNuisSKDetFSI031) :
	  nutypes = "all NC ($1R_{e}$)"; break;
        default : nutypes = "--";
	}
      return Form("SKDet + FSI/SI %2d; Ereco range %.2f - %.2f GeV; %s",ip,
		  kSKDetFSIParamsErecoCut[ip][0],kSKDetFSIParamsErecoCut[ip][1],nutypes.c_str());     
      break;
      }
#endif
#endif

    case ( kNuisRHCToFHC ) : return "RHC to FHC normalisation"; break;

    default : 
      return "--";
  }
  return "--";
}
//____________________________________________________________________________
string t2kana::util::AsStringPresentTex(PhysParamId_t p)
{
  
  switch(p) {

    case ( kPhypMixing23 ) : 
    {
        NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
        switch(omp) {
           case(kNuMixSinsq2Theta) : return "$\\sin^2(2\\theta_{23})$"; break;
           case(kNuMixSinsqTheta ) : return "$\\sin^2(\\theta_{23})$";   break;
           default                 : return "--";                break;
         }
    }
    break;
    case ( kPhypMixing13 ) : 
    {
        NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
        switch(omp) {
           case(kNuMixSinsq2Theta) : return "$\\sin^2(2\\theta_{13})$";  break;
           case(kNuMixSinsqTheta ) : return "$\\sin^2(\\theta_{13})$";    break;
           default                 : return "--";                 break;
         }
    }
    break;
    case ( kPhypMixing12 ) :      
    {
        NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
        switch(omp) {
           case(kNuMixSinsq2Theta) : return "$\\sin^2(2\\theta_{12})$";  break;
           case(kNuMixSinsqTheta ) : return "$\\sin^2(\\theta_{12})$";    break;
           default                 : return "--";                 break;
         }
    }
    break;

    case ( kPhypDeltaCP               ) : return "$\\delta_{CP}$";                   break;
    case ( kPhypDm12sq                ) : return "$|\\Delta m^2_{21}|$";                  break;
    case ( kPhypDm23sq                ) :
      {
	NuOscDmsqScheme_t dmscheme = Opt::Instance()->gOscDmsqScheme;
	switch(dmscheme) {
	  case(kNuDmsq12andAbs23) : return "|\\Delta m^2_{32}|";                  break;
	  case(kNuDmsq12andAbs31) : return "|\\Delta m^2_{31}|";                  break;
	  case(kNuDmsqFogliLisi ) : return "|\\Delta m^2|";                       break;
  	  default                 : return "--";                                  break;
	}
      }
    break;

    case ( kPhypAlphaDecay            ) : return "$\\alpha_{decay}$";                break;
    case ( kPhypSinsq2ThetaDecay      ) : return "$\\sin^2(2\\theta_{decay})$";       break;
    case ( kPhypMusqDecoh             ) : return "$\\mu^2_{decoh}$";                 break;
    case ( kPhypSinsq2ThetaDecoh      ) : return "$\\sin^2(2\\theta_{decoh})$";       break;

    //
    // 2011-analysis nuisance parameters
    //
  
      //
    // SK detector
    //
#ifndef SKFSITOGETHER
#if FSI_BINS != 3
    case ( kNuisSuperkEnergyScale     ) : return "SK energy scale";                 break;    
    case ( kNuisSuperkEff1RmuCCQEoth  ) : return "SK efficiency, 1Rmu, CCQE//Other"; break;
    case ( kNuisSuperkEff1RmuCCQErc0u ) : return "SK efficiency, 1Rmu, CCQE//RC0";   break;
    case ( kNuisSuperkEff1RmuCCQErc1u ) : return "SK efficiency, 1Rmu, CCQE//RC1";   break;
    case ( kNuisSuperkEff1RmuCCQErc2u ) : return "SK efficiency, 1Rmu, CCQE//RC2";   break;
    case ( kNuisSuperkEff1RmuCCnonQE  ) : return "SK efficiency, 1Rmu, CCnonQE";    break;
    case ( kNuisSuperkEff1RmuNC       ) : return "SK efficiency, 1Rmu, NC";         break;
    case ( kNuisSuperkEff1RmuNueCC    ) : return "SK efficiency, 1Rmu, nue CC";     break;
#endif
#if FSI_BINS != 2
    case ( kNuisSuperkEff1ReNueSigCC0 ) : return "SK efficiency, 1Re, Osc. $\\nu_{e}$ CC E0"; break;
    case ( kNuisSuperkEff1ReNueSigCC1 ) : return "SK efficiency, 1Re, Osc. $\\nu_{e}$ CC E1"; break;
    case ( kNuisSuperkEff1ReNueSigCC2 ) : return "SK efficiency, 1Re, Osc. $\\nu_{e}$ CC E2"; break;
    case ( kNuisSuperkEff1ReNumuCC0 )   : return "SK efficiency, 1Re, Osc. $\\nu_{\\mu}$ CC E0"; break;
    case ( kNuisSuperkEff1ReNumuCC1 )   : return "SK efficiency, 1Re, Osc. $\\nu_{\\mu}$ CC E1"; break;
    case ( kNuisSuperkEff1ReNumuCC2 )   : return "SK efficiency, 1Re, Osc. $\\nu_{\\mu}$ CC E2"; break;
    case ( kNuisSuperkEff1ReNueCC0 ) : return "SK efficiency, 1Re, $\\nu_{e}$ CC E0"; break;
    case ( kNuisSuperkEff1ReNueCC1 ) : return "SK efficiency, 1Re, $\\nu_{e}$ CC E1"; break;
    case ( kNuisSuperkEff1ReNueCC2 ) : return "SK efficiency, 1Re, $\\nu_{e}$ CC E2"; break;
    case ( kNuisSuperkEff1ReNC0 )    : return "SK efficiency, 1Re, NC E0"; break;
    case ( kNuisSuperkEff1ReNC1 )    : return "SK efficiency, 1Re, NC E1"; break;
    case ( kNuisSuperkEff1ReNC2 )    : return "SK efficiency, 1Re, NC E2"; break;
#endif
#endif

    //
    // BANFF flux and cross-section params
    //
    case ( kNuisBANFF00 ) :
    case ( kNuisBANFF01 ) :
    case ( kNuisBANFF02 ) :
    case ( kNuisBANFF03 ) :
    case ( kNuisBANFF04 ) :
    case ( kNuisBANFF05 ) :
    case ( kNuisBANFF06 ) :
    case ( kNuisBANFF07 ) :
    case ( kNuisBANFF08 ) :
    case ( kNuisBANFF09 ) :
    case ( kNuisBANFF10 ) :
    case ( kNuisBANFF11 ) :
    case ( kNuisBANFF12 ) :
    case ( kNuisBANFF13 ) :
    case ( kNuisBANFF14 ) :
    case ( kNuisBANFF15 ) :
#if (BANFF_BINS != 2) && (BANFF_BINS != 4)
    case ( kNuisBANFF16 ) :
    case ( kNuisBANFF17 ) :
    case ( kNuisBANFF18 ) :
    case ( kNuisBANFF19 ) :
    case ( kNuisBANFF20 ) :
    case ( kNuisBANFF21 ) :
#endif
#if BANFF_BINS == 1
    case ( kNuisBANFF22 ) :
    case ( kNuisBANFF23 ) :
    case ( kNuisBANFF24 ) :
#endif
    {
      int ip = (int)p - (int)kNuisBANFF00;
      return Form("%s flux %.1f-%.1f GeV",
		  util::Pdg2StringTex(kNuPdgCodeOfFluxParamInBANFFfit[ip]).c_str(), 
		  kEtrueRangeOfFluxParamInBANFFfit[ip][0],
		  kEtrueRangeOfFluxParamInBANFFfit[ip][1]);     
      break;
    }
#ifndef ONEOVERMAQE
    case ( kNuisBANFFMAQE    ) : return "M$_{A}$QE";                           break;
#else
    case ( kNuisBANFFMAQE    ) : return "1/(M$_{A}$QE)$^2$";                 break;
#endif
    case ( kNuisBANFFMARES   ) : return "M$_{A}$RES";                          break;
    case ( kNuisBANFFCCQEE0  ) : return "CCQE 0.0-1.5 GeV";      break;
#if BANFF_BINS != 3
    case ( kNuisBANFFCCQEE1  ) : return "CCQE 1.5-3.5 GeV";      break;
    case ( kNuisBANFFCCQEE2  ) : return "CCQE 3.5-30 GeV";     break;
#endif
    case ( kNuisBANFFCC1piE0 ) : return "CC1$\\pi$ 0.0-2.5 GeV";     break;
#if BANFF_BINS != 3
    case ( kNuisBANFFCC1piE1 ) : return "CC1$\\pi$ 2.5-30 GeV";    break;
#endif
#if (BANFF_BINS != 2) && (BANFF_BINS != 4)
    case ( kNuisBANFFNC1pi0  ) : return "BANFF; Norm; NC 1pi0";    break;
#endif
#if BANFF_BINS == 4
    case ( kNuisBANFFMECNorm ) : return "BANFF; Norm; MEC";                       break;
    case ( kNuisBANFFPDDDoubleCount ) : return "BANFF; PDD double counting";      break;
#endif

    //
    // XSec params not correlated with flux
    //

    case ( kNuisXSecFermiMom   ) : return "Fermi Momentum";                  break;
    case ( kNuisXSecSpecFunc   ) : return "Spectral Func";                   break; 
    case ( kNuisXSecWshape     ) : return "W-shape";                         break;
    case ( kNuisXSecShapeCCoth ) : return "CC Other shape";                   break;
    case ( kNuisXSecPiLessDelta) : return "$\\pi$-less $\\Delta$-decay";             break;
#ifndef __NUEONLY__
    case ( kNuisXSecBindE      ) : return "Binding energy";                  break;
#endif
    case ( kNuisXSecNormCCcoh  ) : return "CC coherent norm";                    break;
#ifndef __NUEONLY__
    case ( kNuisXSecNormNC1piC ) : return "NC 1$\\pi^\\pm$ norm";                   break;
#endif
#ifndef __NUMUONLY__
    case ( kNuisXSecNormNCcoh  ) : return "NC coherent norm";                    break;
#endif
    case ( kNuisXSecNormNCoth  ) : return "NC other norm";                    break;
    case ( kNuisXSecNormNuEToNuMu  ) : return "CC $\\nu_e$/$\\nu_\\mu$ norm";  break;
#ifndef __NUEONLY__
    case ( kNuisXSecNormAntiNuToNu ) : return "CC $\\bar{\\nu}$/$\\nu$ norm";  break;
#endif
#ifndef SKFSITOGETHER
#if FSI_BINS != 3
    case ( kNuisFSI000 ) :
    case ( kNuisFSI001 ) :
    case ( kNuisFSI002 ) :
    case ( kNuisFSI003 ) :
    case ( kNuisFSI004 ) :
    case ( kNuisFSI005 ) :
    {
      string nutypes = "";
      int ip = (int)p - (int)kNuisFSI000;
      switch(p) {
        case (kNuisFSI000) :
        case (kNuisFSI001) :
        case (kNuisFSI002) :
	  nutypes = "$\\nu_\\mu$/$\\bar{\\nu_\\mu}$ CCQE ($1R_{\\mu}$)"; break;
        case (kNuisFSI003) :
	  nutypes = "$\\nu_\\mu$/$\\bar{\\nu_\\mu}$ CCnQE ($1R_{\\mu}$)"; break;
        case (kNuisFSI004) :
	  nutypes = "$\\nu_e$/$\\bar{\\nu_e}$/sig$\\nu_e$ CC ($1R_{\\mu}$)"; break;
        case (kNuisFSI005) :
	  nutypes = "all NC ($1R_{\\mu}$)"; break;
        default : nutypes = "--";
      }
      return Form("FSI/SI %2d; Ereco range %.1f - %.1f GeV; %s",ip,
		  kFSIParamsErecoCut[ip][0],kFSIParamsErecoCut[ip][1],nutypes.c_str());     
      break;
    }//FSI for 1Rmu
#endif
#if FSI_BINS != 2
    case ( kNuisFSI006 ) :
    case ( kNuisFSI007 ) :
    case ( kNuisFSI008 ) :
    case ( kNuisFSI009 ) :
    case ( kNuisFSI010 ) :
    case ( kNuisFSI011 ) :
    case ( kNuisFSI012 ) :
    case ( kNuisFSI013 ) :
    case ( kNuisFSI014 ) :
    case ( kNuisFSI015 ) :
    case ( kNuisFSI016 ) :
    case ( kNuisFSI017 ) :

    {
      string nutypes = "";
#if FSI_BINS != 3
      int ip = (int)p - (int)kNuisFSI000;
#else
      int ip = (int)p - (int)kNuisFSI006;
#endif
      switch(p) {
        case (kNuisFSI006) :
        case (kNuisFSI007) :
        case (kNuisFSI008) :
	  nutypes = "Osc. $\\nu_e$ CC ($1R_{e}$)"; break;
        case (kNuisFSI009) :
        case (kNuisFSI010) :
        case (kNuisFSI011) :
	  nutypes = "$\\nu_\\mu$/$\\bar{\\nu_\\mu}$ CC ($1R_{e}$)"; break;
        case (kNuisFSI012) :
        case (kNuisFSI013) :
        case (kNuisFSI014) :
	  nutypes = "$\\nu_e$/$\\bar{\\nu_e}$ CC ($1R_{e}$)"; break;
        case (kNuisFSI015) :
        case (kNuisFSI016) :
        case (kNuisFSI017) :
	  nutypes = "all NC ($1R_{e}$)"; break;
        default : nutypes = "--";
      }
      return Form("FSI/SI %2d; Ereco range %.1f - %.1f GeV; %s",ip,
		  kFSIParamsErecoCut[ip][0],kFSIParamsErecoCut[ip][1],nutypes.c_str());     
      break;
    }//FSI for 1Re
#endif
#endif
#ifdef SKFSITOGETHER
#ifndef SKCORRELATIONS
#if FSI_BINS != 3
    case ( kNuisSKDetFSI000 ) :
    case ( kNuisSKDetFSI001 ) :
    case ( kNuisSKDetFSI002 ) :
    case ( kNuisSKDetFSI003 ) :
    case ( kNuisSKDetFSI004 ) :
    case ( kNuisSKDetFSI005 ) :
    {
      string nutypes = "";
      switch(p) {
        case (kNuisSKDetFSI000) :
	  nutypes = "SK+FSI $\\nu_\\mu$/$\\bar{\\nu_\\mu}$CCQE $<$0.4 GeV"; break;
        case (kNuisSKDetFSI001) :
	  nutypes = "SK+FSI $\\nu_\\mu$/$\\bar{\\nu_\\mu}$CCQE 0.4-1.1 GeV"; break;
        case (kNuisSKDetFSI002) :
	  nutypes = "SK+FSI $\\nu_\\mu$/$\\bar{\\nu_\\mu}$CCQE $>$1.1 GeV"; break;
        case (kNuisSKDetFSI003) :
	  nutypes = "SK+FSI $\\nu_\\mu$/$\\bar{\\nu_\\mu}$ CCnQE"; break;
        case (kNuisSKDetFSI004) :
	  nutypes = "SK+FSI $\\nu_e$/$\\bar{\\nu_e}$/sig$\\nu_e$ CC"; break;
        case (kNuisSKDetFSI005) :
	  nutypes = "SK+FSI NC"; break;
        default : nutypes = "--";
      }
      return Form("%s",nutypes.c_str());     
      break;
    }//FSI for 1Rmu
#endif
#if FSI_BINS != 2
    case ( kNuisSKDetFSI006 ) :
    case ( kNuisSKDetFSI007 ) :
    case ( kNuisSKDetFSI008 ) :
    case ( kNuisSKDetFSI009 ) :
    case ( kNuisSKDetFSI010 ) :
    case ( kNuisSKDetFSI011 ) :
    case ( kNuisSKDetFSI012 ) :
    case ( kNuisSKDetFSI013 ) :
    case ( kNuisSKDetFSI014 ) :
    case ( kNuisSKDetFSI015 ) :
    case ( kNuisSKDetFSI016 ) :
    case ( kNuisSKDetFSI017 ) :

    {
      string nutypes = "";
#if FSI_BINS != 3
      int ip = (int)p - (int)kNuisSKDetFSI000;
#else
      int ip = (int)p - (int)kNuisSKDetFSI006;
#endif
      switch(p) {
        case (kNuisSKDetFSI006) :
        case (kNuisSKDetFSI007) :
        case (kNuisSKDetFSI008) :
	  nutypes = "Osc. $\\nu_e$ CC ($1R_{e}$)"; break;
        case (kNuisSKDetFSI009) :
        case (kNuisSKDetFSI010) :
        case (kNuisSKDetFSI011) :
	  nutypes = "$\\nu_\\mu$/$\\bar{\\nu_\\mu}$ CC ($1R_{e}$)"; break;
        case (kNuisSKDetFSI012) :
        case (kNuisSKDetFSI013) :
        case (kNuisSKDetFSI014) :
	  nutypes = "$\\nu_e$/$\\bar{\\nu_e}$ CC ($1R_{e}$)"; break;
        case (kNuisSKDetFSI015) :
        case (kNuisSKDetFSI016) :
        case (kNuisSKDetFSI017) :
	  nutypes = "all NC ($1R_{e}$)"; break;
        default : nutypes = "--";
      }
      return Form("SKDet + FSI/SI %2d; Ereco range %.1f - %.1f GeV; %s",ip,
		  kSKDetFSIParamsErecoCut[ip][0],kSKDetFSIParamsErecoCut[ip][1],nutypes.c_str());     
      break;
    }//FSI for 1Re
#endif
#else
    case ( kNuisSKDetFSI000 ) :
    case ( kNuisSKDetFSI001 ) :
    case ( kNuisSKDetFSI002 ) :
    case ( kNuisSKDetFSI003 ) :
    case ( kNuisSKDetFSI004 ) :
    case ( kNuisSKDetFSI005 ) :
    case ( kNuisSKDetFSI006 ) :
    case ( kNuisSKDetFSI007 ) :
    case ( kNuisSKDetFSI008 ) :
    case ( kNuisSKDetFSI009 ) :
    case ( kNuisSKDetFSI010 ) :
    case ( kNuisSKDetFSI011 ) :
    case ( kNuisSKDetFSI012 ) :
    case ( kNuisSKDetFSI013 ) :
    case ( kNuisSKDetFSI014 ) :
    case ( kNuisSKDetFSI015 ) :
    case ( kNuisSKDetFSI016 ) :
    case ( kNuisSKDetFSI017 ) :
    case ( kNuisSKDetFSI018 ) :
    case ( kNuisSKDetFSI019 ) :
    {
      int ip = (int)p - (int)kNuisSKDetFSI000;
      string nutypes = "";
      switch(p) {
        case (kNuisSKDetFSI000) :
        case (kNuisSKDetFSI001) :
        case (kNuisSKDetFSI002) :
        case (kNuisSKDetFSI003) :
        case (kNuisSKDetFSI004) :
	  nutypes = "Osc. $\\nu_e$ CC ($1R_{\\mu}$)"; break;
        case (kNuisSKDetFSI005) :
        case (kNuisSKDetFSI006) :
        case (kNuisSKDetFSI007) :
        case (kNuisSKDetFSI008) :
        case (kNuisSKDetFSI009) :
	  nutypes = "$\\nu_\\mu$/$\\bar{\\nu_\\mu}$ CC ($1R_{\\mu}$)"; break;
        case (kNuisSKDetFSI010) :
        case (kNuisSKDetFSI011) :
        case (kNuisSKDetFSI012) :
        case (kNuisSKDetFSI013) :
        case (kNuisSKDetFSI014) :
	  nutypes = "$\\nu_e$/$\\bar{\\nu_e}$ CC ($1R_{\\mu}$)"; break;
        case (kNuisSKDetFSI015) :
        case (kNuisSKDetFSI016) :
        case (kNuisSKDetFSI017) :
        case (kNuisSKDetFSI018) :
        case (kNuisSKDetFSI019) :
	  nutypes = "all NC ($1R_{\\mu}$)"; break;
        default : nutypes = "--";
      }
      return Form("SKDet + FSI/SI %2d; Ereco range %.2f - %.2f GeV; %s",ip,
		  kSKDetFSIParamsErecoCut[ip][0],kSKDetFSIParamsErecoCut[ip][1],nutypes.c_str());     
      break;
    }//FSI for 1Rmu
    case ( kNuisSKDetFSI020 ) :
    case ( kNuisSKDetFSI021 ) :
    case ( kNuisSKDetFSI022 ) :
    case ( kNuisSKDetFSI023 ) :
    case ( kNuisSKDetFSI024 ) :
    case ( kNuisSKDetFSI025 ) :
    case ( kNuisSKDetFSI026 ) :
    case ( kNuisSKDetFSI027 ) :
    case ( kNuisSKDetFSI028 ) :
    case ( kNuisSKDetFSI029 ) :
    case ( kNuisSKDetFSI030 ) :
    case ( kNuisSKDetFSI031 ) :
      {   
      int ip = (int)p - (int)kNuisSKDetFSI000;
      string nutypes = "";
	switch(p) {
        case (kNuisSKDetFSI020) :
        case (kNuisSKDetFSI021) :
        case (kNuisSKDetFSI022) :
	  nutypes = "Osc. $\\nu_e$ CC ($1R_{e}$)"; break;
        case (kNuisSKDetFSI023) :
        case (kNuisSKDetFSI024) :
        case (kNuisSKDetFSI025) :
	  nutypes = "$\\nu_\\mu$/$\\bar{\\nu_\\mu}$ CC ($1R_{e}$)"; break;
        case (kNuisSKDetFSI026) :
        case (kNuisSKDetFSI027) :
        case (kNuisSKDetFSI028) :
	  nutypes = "$\\nu_e$/$\\bar{\\nu_e}$ CC ($1R_{e}$)"; break;
        case (kNuisSKDetFSI029) :
        case (kNuisSKDetFSI030) :
        case (kNuisSKDetFSI031) :
	  nutypes = "all NC ($1R_{e}$)"; break;
        default : nutypes = "--";
	}
      return Form("SKDet + FSI/SI %2d; Ereco range %.2f - %.2f GeV; %s",ip,
		  kSKDetFSIParamsErecoCut[ip][0],kSKDetFSIParamsErecoCut[ip][1],nutypes.c_str());     
      break;
      }
#endif
#endif

    case ( kNuisRHCToFHC ) : return "RHC to FHC normalisation"; break;

    default : 
      return "--";
  }
  
  return "--";
}
//____________________________________________________________________________                                                                                                    
string t2kana::util::AsStringPaperTex(PhysParamId_t p)
{
  switch(p) {
    
  case ( kPhypMixing23 ) :
    {
      NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
      switch(omp) {
      case(kNuMixSinsq2Theta) : return "$\\sin^2(2\\theta_{23})$"; break;
      case(kNuMixSinsqTheta ) : return "$\\sin^2(\\theta_{23})$";   break;
      default                 : return "--";                break;
      }
    }
    break;
  case ( kPhypMixing13 ) :
    {
      NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
      switch(omp) {
      case(kNuMixSinsq2Theta) : return "$\\sin^2(2\\theta_{13})$";  break;
      case(kNuMixSinsqTheta ) : return "$\\sin^2(\\theta_{13})$";    break;
      default                 : return "--";                 break;
      }
    }
    break;
  case ( kPhypMixing12 ) :
    {
      NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
      switch(omp) {
      case(kNuMixSinsq2Theta) : return "$\\sin^2(2\\theta_{12})$";  break;
      case(kNuMixSinsqTheta ) : return "$\\sin^2(\\theta_{12})$";    break;
      default                 : return "--";                 break;
      }
    }
    break;
    
  case ( kPhypDeltaCP               ) : return "$\\delta_{CP}$";                   break;
  case ( kPhypDm23sq                ) : 
    {
      NuOscDmsqScheme_t omp = Opt::Instance()->gOscDmsqScheme;
      switch(omp) {
        case(kNuDmsqFogliLisi) :	return "$|\\Deltam^2| = |m^2_{3} - (m^2_{1} + m^2_{1})/2|$";                  break;
        case(kNuDmsq12andAbs23):	return "$|\\Deltam^2_{32}|$";                  break;
        case(kNuDmsq12andAbs31):	return "$|\\Deltam^2_{31}|$";                  break;
        default                :        return "--";                                   break;
      }
    }
    break;

  case ( kPhypDm12sq                ) : 
    {
      NuOscDmsqScheme_t omp = Opt::Instance()->gOscDmsqScheme;
      switch(omp) {
      case(kNuDmsqFogliLisi):
	return "$\\deltam^2 = \\Delta m^2_{21}$";      break;
      case(kNuDmsq12andAbs23):
      case(kNuDmsq12andAbs31):
	return "$\\Delta m^2_{21}$";                   break;
      default:
	return "--";                                   break;
      }
    }
    
  case ( kPhypAlphaDecay            ) : return "$\\alpha_{decay}$";                break;
  case ( kPhypSinsq2ThetaDecay      ) : return "$\\sin^2(2\\theta_{decay})$";       break;
  case ( kPhypMusqDecoh             ) : return "$\\mu^2_{decoh}$";                 break;
  case ( kPhypSinsq2ThetaDecoh      ) : return "$\\sin^2(2\\theta_{decoh})$";       break;
  
    /*                                                     
     * 2011-analysis nuisance parameters                                                                                                                                            
     */

    //
    // FHC systematics
    //

    // 
    // SK detector  
    //                                                                                                                                                                                 
  case ( kNuisSuperkEnergyScale     ) : return "$f_{E; r}^{SK}$";                 break;

    //                                                          
    // BANFF flux and cross-section params                                                                                                                                       
    //                                                                                                                                                                               
  case ( kNuisBANFF00 ) :
  case ( kNuisBANFF01 ) :
  case ( kNuisBANFF02 ) :
  case ( kNuisBANFF03 ) :
  case ( kNuisBANFF04 ) :
  case ( kNuisBANFF05 ) :
  case ( kNuisBANFF06 ) :
  case ( kNuisBANFF07 ) :
  case ( kNuisBANFF08 ) :
  case ( kNuisBANFF09 ) :
  case ( kNuisBANFF10 ) :
  case ( kNuisBANFF11 ) :
  case ( kNuisBANFF12 ) :
  case ( kNuisBANFF13 ) :
  case ( kNuisBANFF14 ) :
  case ( kNuisBANFF15 ) :
#if (BANFF_BINS != 2) && (BANFF_BINS != 4)
  case ( kNuisBANFF16 ) :
  case ( kNuisBANFF17 ) :
  case ( kNuisBANFF18 ) :
  case ( kNuisBANFF19 ) :
  case ( kNuisBANFF20 ) :
  case ( kNuisBANFF21 ) :
#endif
#if BANFF_BINS == 1
  case ( kNuisBANFF22 ) :
  case ( kNuisBANFF23 ) :
  case ( kNuisBANFF24 ) :
#endif
  case ( kNuisBANFFMAQE ) :
  case ( kNuisBANFFMARES ) :
  case ( kNuisBANFFCCQEE0 ) :
#if BANFF_BINS != 3
  case ( kNuisBANFFCCQEE1 ) :
  case ( kNuisBANFFCCQEE2 ) :
#endif
  case ( kNuisBANFFCC1piE0 ) :
#if BANFF_BINS != 3
  case ( kNuisBANFFCC1piE1 ) :
#endif
#if (BANFF_BINS != 2) && (BANFF_BINS != 4)
  case ( kNuisBANFFNC1pi0 ) :
#endif
#if BANFF_BINS == 4
    case ( kNuisBANFFMECNorm ) :
    case ( kNuisBANFFPDDDoubleCount ) :
#endif
    {
      int ip = (int)p - (int)kNuisBANFF00;
      return Form("$f_{%d; t, r}^{banff}$", ip);
      break;
    }
    
    //
    // XSec params not correlated with flux 
    // 

  case ( kNuisXSecFermiMom   ) : return "$f_{pF; t, r}$";                  break;
  case ( kNuisXSecSpecFunc   ) : return "$f_{SF; t, r}$";                  break;
  case ( kNuisXSecWshape     ) : return "$f_{Wshape; t, r}$";              break;
  case ( kNuisXSecShapeCCoth ) : return "$f_{CCothShape; t, r}$";          break;
  case ( kNuisXSecPiLessDelta) : return "$f_{\\pi-less\\Delta; t, r}$";    break;
#ifndef __NUEONLY__
  case ( kNuisXSecBindE      ) : return "$f_{bindE; t, r}$";               break;
#endif
  case ( kNuisXSecNormCCcoh  ) : return "$f_{CCcoh; t}$";                  break;
#ifndef __NUMUONLY__
  case ( kNuisXSecNormNCcoh  ) : return "$f_{NCcoh; t}$";                  break;
#endif
#ifndef __NUEONLY__
  case ( kNuisXSecNormNC1piC ) : return "$f_{NC1\\pi^\\pm; t}$";           break;
#endif
  case ( kNuisXSecNormNCoth  ) : return "$f_{NCoth; t}$";                  break;
  case ( kNuisXSecNormNuEToNuMu  ) : return "$f_{CC\\nu_e; t}$";           break;
#ifndef __NUEONLY__
  case ( kNuisXSecNormAntiNuToNu ) : return "$f_{CC\\bar{\\nu}; t}$";      break;
#endif

#if FSI_BINS != 3                
  case ( kNuisSKDetFSI000 ) :
  case ( kNuisSKDetFSI001 ) :
  case ( kNuisSKDetFSI002 ) :
  case ( kNuisSKDetFSI003 ) :
  case ( kNuisSKDetFSI004 ) :
  case ( kNuisSKDetFSI005 ) :
#endif
#if FSI_BINS != 2
  case ( kNuisSKDetFSI006 ) :
  case ( kNuisSKDetFSI007 ) :
  case ( kNuisSKDetFSI008 ) :
  case ( kNuisSKDetFSI009 ) :
  case ( kNuisSKDetFSI010 ) :
  case ( kNuisSKDetFSI011 ) :
  case ( kNuisSKDetFSI012 ) :
  case ( kNuisSKDetFSI013 ) :
  case ( kNuisSKDetFSI014 ) :
  case ( kNuisSKDetFSI015 ) :
  case ( kNuisSKDetFSI016 ) :
  case ( kNuisSKDetFSI017 ) :
#endif
#ifdef SKCORRELATIONS
  case ( kNuisSKDetFSI018 ) :
  case ( kNuisSKDetFSI019 ) :
  case ( kNuisSKDetFSI020 ) :
  case ( kNuisSKDetFSI021 ) :
  case ( kNuisSKDetFSI022 ) :
  case ( kNuisSKDetFSI023 ) :
  case ( kNuisSKDetFSI024 ) :
  case ( kNuisSKDetFSI025 ) :
  case ( kNuisSKDetFSI026 ) :
  case ( kNuisSKDetFSI027 ) :
  case ( kNuisSKDetFSI028 ) :
  case ( kNuisSKDetFSI029 ) :
  case ( kNuisSKDetFSI030 ) :
  case ( kNuisSKDetFSI031 ) :
#endif

    {
#if FSI_BINS != 3
      int ip = (int)p - (int)kNuisSKDetFSI000;
#else
      int ip = (int)p - (int)kNuisSKDetFSI006;
#endif
      return Form("$f_{%d; t, r}^{SK+FSI}$", ip);
      break;
    }//FSI

  case ( kNuisRHCToFHC ) : return "$f^{RHC/FHC}$"; break;


    //
    // RHC systematics
    //

    // 
    // SK detector  
    //                                                                                                                                                                                 
  case ( kRHCNuisSuperkEnergyScale     ) : return "$f_{E; r}^{SK}RHC$";                 break;

    //                                                          
    // BANFF flux and cross-section params                                                                                                                                       
    //                                                                                                                                                                               
  case ( kRHCNuisBANFF00 ) :
  case ( kRHCNuisBANFF01 ) :
  case ( kRHCNuisBANFF02 ) :
  case ( kRHCNuisBANFF03 ) :
  case ( kRHCNuisBANFF04 ) :
  case ( kRHCNuisBANFF05 ) :
  case ( kRHCNuisBANFF06 ) :
  case ( kRHCNuisBANFF07 ) :
  case ( kRHCNuisBANFF08 ) :
  case ( kRHCNuisBANFF09 ) :
  case ( kRHCNuisBANFF10 ) :
  case ( kRHCNuisBANFF11 ) :
  case ( kRHCNuisBANFF12 ) :
  case ( kRHCNuisBANFF13 ) :
  case ( kRHCNuisBANFF14 ) :
  case ( kRHCNuisBANFF15 ) :
#if (BANFF_BINS != 2) && (BANFF_BINS != 4)
  case ( kRHCNuisBANFF16 ) :
  case ( kRHCNuisBANFF17 ) :
  case ( kRHCNuisBANFF18 ) :
  case ( kRHCNuisBANFF19 ) :
  case ( kRHCNuisBANFF20 ) :
  case ( kRHCNuisBANFF21 ) :
#endif
#if BANFF_BINS == 1
  case ( kRHCNuisBANFF22 ) :
  case ( kRHCNuisBANFF23 ) :
  case ( kRHCNuisBANFF24 ) :
#endif
  case ( kRHCNuisBANFFMAQE ) :
  case ( kRHCNuisBANFFMARES ) :
  case ( kRHCNuisBANFFCCQEE0 ) :
#if BANFF_BINS != 3
  case ( kRHCNuisBANFFCCQEE1 ) :
  case ( kRHCNuisBANFFCCQEE2 ) :
#endif
  case ( kRHCNuisBANFFCC1piE0 ) :
#if BANFF_BINS != 3
  case ( kRHCNuisBANFFCC1piE1 ) :
#endif
#if (BANFF_BINS != 2) && (BANFF_BINS != 4)
  case ( kRHCNuisBANFFNC1pi0 ) :
#endif
#if BANFF_BINS == 4
    case ( kRHCNuisBANFFMECNorm ) :
    case ( kRHCNuisBANFFPDDDoubleCount ) :
#endif
    {
      int ip = (int)p - (int)kRHCNuisBANFF00;
      return Form("$f_{%d; t, r}^{banff}RHC$", ip);
      break;
    }
    
    //
    // XSec params not correlated with flux 
    // 

  case ( kRHCNuisXSecFermiMom   ) : return "$f_{pF; t, r}RHC$";                  break;
  case ( kRHCNuisXSecSpecFunc   ) : return "$f_{SF; t, r}RHC$";                  break;
  case ( kRHCNuisXSecWshape     ) : return "$f_{Wshape; t, r}RHC$";              break;
  case ( kRHCNuisXSecShapeCCoth ) : return "$f_{CCothShape; t, r}RHC$";          break;
  case ( kRHCNuisXSecPiLessDelta) : return "$f_{\\pi-less\\Delta; t, r}RHC$";    break;
#ifndef __NUEONLY__
  case ( kRHCNuisXSecBindE      ) : return "$f_{bindE; t, r}RHC$";               break;
#endif
  case ( kRHCNuisXSecNormCCcoh  ) : return "$f_{CCcoh; t}RHC$";                  break;
#ifndef __NUMUONLY__
  case ( kRHCNuisXSecNormNCcoh  ) : return "$f_{NCcoh; t}RHC$";                  break;
#endif
#ifndef __NUEONLY__
  case ( kRHCNuisXSecNormNC1piC ) : return "$f_{NC1\\pi^\\pm; t}RHC$";           break;
#endif
  case ( kRHCNuisXSecNormNCoth  ) : return "$f_{NCoth; t}RHC$";                  break;
  case ( kRHCNuisXSecNormNuEToNuMu  ) : return "$f_{CC\\nu_e; t}RHC$";           break;
#ifndef __NUEONLY__
  case ( kRHCNuisXSecNormAntiNuToNu ) : return "$f_{CC\\bar{\\nu}; t}RHC$";      break;
#endif

#if FSI_BINS != 3                
  case ( kRHCNuisSKDetFSI000 ) :
  case ( kRHCNuisSKDetFSI001 ) :
  case ( kRHCNuisSKDetFSI002 ) :
  case ( kRHCNuisSKDetFSI003 ) :
  case ( kRHCNuisSKDetFSI004 ) :
  case ( kRHCNuisSKDetFSI005 ) :
#endif
#if FSI_BINS != 2
  case ( kRHCNuisSKDetFSI006 ) :
  case ( kRHCNuisSKDetFSI007 ) :
  case ( kRHCNuisSKDetFSI008 ) :
  case ( kRHCNuisSKDetFSI009 ) :
  case ( kRHCNuisSKDetFSI010 ) :
  case ( kRHCNuisSKDetFSI011 ) :
  case ( kRHCNuisSKDetFSI012 ) :
  case ( kRHCNuisSKDetFSI013 ) :
  case ( kRHCNuisSKDetFSI014 ) :
  case ( kRHCNuisSKDetFSI015 ) :
  case ( kRHCNuisSKDetFSI016 ) :
  case ( kRHCNuisSKDetFSI017 ) :
#endif
#ifdef SKCORRELATIONS
  case ( kRHCNuisSKDetFSI018 ) :
  case ( kRHCNuisSKDetFSI019 ) :
  case ( kRHCNuisSKDetFSI020 ) :
  case ( kRHCNuisSKDetFSI021 ) :
  case ( kRHCNuisSKDetFSI022 ) :
  case ( kRHCNuisSKDetFSI023 ) :
  case ( kRHCNuisSKDetFSI024 ) :
  case ( kRHCNuisSKDetFSI025 ) :
  case ( kRHCNuisSKDetFSI026 ) :
  case ( kRHCNuisSKDetFSI027 ) :
  case ( kRHCNuisSKDetFSI028 ) :
  case ( kRHCNuisSKDetFSI029 ) :
  case ( kRHCNuisSKDetFSI030 ) :
  case ( kRHCNuisSKDetFSI031 ) :
#endif

    {
#if FSI_BINS != 3
      int ip = (int)p - (int)kRHCNuisSKDetFSI000;
#else
      int ip = (int)p - (int)kRHCNuisSKDetFSI006;
#endif
      return Form("$f_{%d; t, r}^{SK+FSI}RHC$", ip);
      break;
    }//FSI

  default :
    return "--";
  }
  return "--";
}
//____________________________________________________________________________                                                                                                    
string t2kana::util::AsStringPaperTexNoDollar(PhysParamId_t p)
{
    switch(p) {
    
  case ( kPhypMixing23 ) :
    {
      NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
      switch(omp) {
      case(kNuMixSinsq2Theta) : return "\\sin^2(2\\theta_{23})"; break;
      case(kNuMixSinsqTheta ) : return "\\sin^2(\\theta_{23})";   break;
      default                 : return "--";                break;
      }
    }
    break;
  case ( kPhypMixing13 ) :
    {
      NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
      switch(omp) {
      case(kNuMixSinsq2Theta) : return "\\sin^2(2\\theta_{13})";  break;
      case(kNuMixSinsqTheta ) : return "\\sin^2(\\theta_{13})";    break;
      default                 : return "--";                 break;
      }
    }
    break;
  case ( kPhypMixing12 ) :
    {
      NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
      switch(omp) {
      case(kNuMixSinsq2Theta) : return "\\sin^2(2\\theta_{12})";  break;
      case(kNuMixSinsqTheta ) : return "\\sin^2(\\theta_{12})";    break;
      default                 : return "--";                 break;
      }
    }
    break;
    
  case ( kPhypDeltaCP               ) : return "\\delta_{CP}";                   break;
  case ( kPhypDm23sq                ) : 
    {
      NuOscDmsqScheme_t omp = Opt::Instance()->gOscDmsqScheme;
      switch(omp) {
        case(kNuDmsqFogliLisi) :	return "|\\Deltam^2| = |m^2_{3} - (m^2_{1} + m^2_{1})/2|";                  break;
        case(kNuDmsq12andAbs23):	return "|\\Deltam^2_{32}|";                  break;
        case(kNuDmsq12andAbs31):	return "|\\Deltam^2_{31}|";                  break;
        default                :        return "--";                                   break;
      }
    }
    break;

  case ( kPhypDm12sq                ) : 
    {
      NuOscDmsqScheme_t omp = Opt::Instance()->gOscDmsqScheme;
      switch(omp) {
      case(kNuDmsqFogliLisi):
	return "\\deltam^2 = \\Delta m^2_{21}";      break;
      case(kNuDmsq12andAbs23):
      case(kNuDmsq12andAbs31):
	return "\\Delta m^2_{21}";                   break;
      default:
	return "--";                                   break;
      }
    }
    
  case ( kPhypAlphaDecay            ) : return "\\alpha_{decay}";                break;
  case ( kPhypSinsq2ThetaDecay      ) : return "\\sin^2(2\\theta_{decay})";       break;
  case ( kPhypMusqDecoh             ) : return "\\mu^2_{decoh}";                 break;
  case ( kPhypSinsq2ThetaDecoh      ) : return "\\sin^2(2\\theta_{decoh})";       break;
  
    //                                                    
    // 2011-analysis nuisance parameters                                                                                                                                            
    //

    // 
    // SK detector  
    //                                                                                                                                                                                 
#ifndef SKFSITOGETHER
#if FSI_BINS != 3  
  case ( kNuisSuperkEnergyScale     ) : return "f_{E; r}^{SK}  ";                 break;
  case ( kNuisSuperkEff1RmuCCQEoth  ) : return "f_{1R_{\\mu}: CCQEoth; r}^{SK}";             break;
  case ( kNuisSuperkEff1RmuCCQErc0u ) : return "f_{1R_{\\mu}: CCQERC0; r}^{SK}";             break;
  case ( kNuisSuperkEff1RmuCCQErc1u ) : return "f_{1R_{\\mu}: CCQERC1; r}^{SK}";             break;
  case ( kNuisSuperkEff1RmuCCQErc2u ) : return "f_{1R_{\\mu}: CCQERC2; r}^{SK}";             break;
  case ( kNuisSuperkEff1RmuCCnonQE  ) : return "f_{1R_{\\mu}: CCnQE; r}^{SK}";               break;
  case ( kNuisSuperkEff1RmuNC       ) : return "f_{1R_{\\mu}: NC; r}^{SK}";                  break;
  case ( kNuisSuperkEff1RmuNueCC    ) : return "f_{1R_{\\mu}: \\nu_e; r}^{SK}";              break;
#endif
#if FSI_BINS != 2
  case ( kNuisSuperkEff1ReNueSigCC0  )   : return "f_{1R_{e}: osc \\nu_e CC E0; r}^{SK}";     break;
  case ( kNuisSuperkEff1ReNueSigCC1  )   : return "f_{1R_{e}: osc \\nu_e CC E1; r}^{SK}";     break;
  case ( kNuisSuperkEff1ReNueSigCC2  )   : return "f_{1R_{e}: osc \\nu_e CC E2; r}^{SK}";     break;
  case ( kNuisSuperkEff1ReNumuCC0    )   : return "f_{1R_{e}: \\nu_{\\mu} CC E0; r}^{SK}";     break;
  case ( kNuisSuperkEff1ReNumuCC1    )   : return "f_{1R_{e}: \\nu_{\\mu} CC E1; r}^{SK}";   break;
  case ( kNuisSuperkEff1ReNumuCC2    )   : return "f_{1R_{e}: \\nu_{\\mu} CC E2; r}^{SK}";   break;
  case ( kNuisSuperkEff1ReNueCC0     )   : return "f_{1R_{e}: \\nu_{e} CC E0; r}^{SK}";   break;
  case ( kNuisSuperkEff1ReNueCC1     )   : return "f_{1R_{e}: \\nu_{e} CC E1; r}^{SK}";   break;
  case ( kNuisSuperkEff1ReNueCC2     )   : return "f_{1R_{e}: \\nu_{e} CC E2; r}^{SK}";   break;
  case ( kNuisSuperkEff1ReNC0        )   : return "f_{1R_{e}: NC E0; r}^{SK}";   break;
  case ( kNuisSuperkEff1ReNC1        )   : return "f_{1R_{e}: NC E1; r}^{SK}";   break;
  case ( kNuisSuperkEff1ReNC2        )   : return "f_{1R_{e}: NC E2; r}^{SK}";   break;
#endif
#endif
    //                                                          
    // BANFF flux and cross-section params                                                                                                                                       
    //                                                                                                                                                                               
  case ( kNuisBANFF00 ) :
  case ( kNuisBANFF01 ) :
  case ( kNuisBANFF02 ) :
  case ( kNuisBANFF03 ) :
  case ( kNuisBANFF04 ) :
  case ( kNuisBANFF05 ) :
  case ( kNuisBANFF06 ) :
  case ( kNuisBANFF07 ) :
  case ( kNuisBANFF08 ) :
  case ( kNuisBANFF09 ) :
  case ( kNuisBANFF10 ) :
  case ( kNuisBANFF11 ) :
  case ( kNuisBANFF12 ) :
  case ( kNuisBANFF13 ) :
  case ( kNuisBANFF14 ) :
  case ( kNuisBANFF15 ) :
#if (BANFF_BINS != 2) && (BANFF_BINS != 4)
  case ( kNuisBANFF16 ) :
  case ( kNuisBANFF17 ) :
  case ( kNuisBANFF18 ) :
  case ( kNuisBANFF19 ) :
  case ( kNuisBANFF20 ) :
  case ( kNuisBANFF21 ) :
#endif
#if BANFF_BINS == 1
  case ( kNuisBANFF22 ) :
  case ( kNuisBANFF23 ) :
  case ( kNuisBANFF24 ) :
#endif
  case ( kNuisBANFFMAQE ) :
  case ( kNuisBANFFMARES ) :
  case ( kNuisBANFFCCQEE0 ) :
#if BANFF_BINS != 3
  case ( kNuisBANFFCCQEE1 ) :
  case ( kNuisBANFFCCQEE2 ) :
#endif
  case ( kNuisBANFFCC1piE0 ) :
#if BANFF_BINS != 3
  case ( kNuisBANFFCC1piE1 ) :
#endif
#if (BANFF_BINS != 2) && (BANFF_BINS != 4)
  case ( kNuisBANFFNC1pi0 ) :
#endif
#if BANFF_BINS == 4
    case ( kNuisBANFFMECNorm ) :
    case ( kNuisBANFFPDDDoubleCount ) :
#endif
    {
      int ip = (int)p - (int)kNuisBANFF00;
      return Form("f_{%d; t, r}^{banff}", ip);
      break;
    }
    
    //
    // XSec params not correlated with flux 
    // 

  case ( kNuisXSecFermiMom   ) : return "f_{pF; t, r}";                  break;
  case ( kNuisXSecSpecFunc   ) : return "f_{SF; t, r}";                  break;
  case ( kNuisXSecWshape     ) : return "f_{Wshape; t, r}";              break;
  case ( kNuisXSecShapeCCoth ) : return "f_{CCothShape; t, r}";          break;
  case ( kNuisXSecPiLessDelta) : return "f_{\\pi-less\\Delta; t, r}";    break;
#ifndef __NUEONLY__
  case ( kNuisXSecBindE      ) : return "f_{bindE; t, r}";               break;
#endif
  case ( kNuisXSecNormCCcoh  ) : return "f_{CCcoh; t}";                  break;
#ifndef __NUMUONLY__
  case ( kNuisXSecNormNCcoh  ) : return "f_{NCcoh; t}";                  break;
#endif
#ifndef __NUEONLY__
  case ( kNuisXSecNormNC1piC ) : return "f_{NC1\\pi^\\pm; t}";           break;
#endif
  case ( kNuisXSecNormNCoth  ) : return "f_{NCoth; t}";                  break;
  case ( kNuisXSecNormNuEToNuMu  ) : return "f_{CC\\nu_e; t}";           break;
#ifndef __NUEONLY__
  case ( kNuisXSecNormAntiNuToNu ) : return "f_{CC\\bar{\\nu}; t}";      break;
#endif
#ifndef SKFSITOGETHER
#if FSI_BINS != 3
  case ( kNuisFSI000 ) :
  case ( kNuisFSI001 ) :
  case ( kNuisFSI002 ) :
  case ( kNuisFSI003 ) :
  case ( kNuisFSI004 ) :
  case ( kNuisFSI005 ) :
#endif
#if FSI_BINS != 2
  case ( kNuisFSI006 ) :
  case ( kNuisFSI007 ) :
  case ( kNuisFSI008 ) :
  case ( kNuisFSI009 ) :
  case ( kNuisFSI010 ) :
  case ( kNuisFSI011 ) :
  case ( kNuisFSI012 ) :
  case ( kNuisFSI013 ) :
  case ( kNuisFSI014 ) :
  case ( kNuisFSI015 ) :
  case ( kNuisFSI016 ) :
  case ( kNuisFSI017 ) :
#endif
#endif
    {
#if FSI_BINS != 3
      int ip = (int)p - (int)kNuisFSI000;
#else
      int ip = (int)p - (int)kNuisFSI006;
#endif
      return Form("f_{%d; r}^{FSI}", ip);
      break;
    }//FSI                
#ifdef SKFSITOGETHER
#if FSI_BINS != 3
  case ( kNuisSKDetFSI000 ) :
  case ( kNuisSKDetFSI001 ) :
  case ( kNuisSKDetFSI002 ) :
  case ( kNuisSKDetFSI003 ) :
  case ( kNuisSKDetFSI004 ) :
  case ( kNuisSKDetFSI005 ) :
#endif
#if FSI_BINS != 2
  case ( kNuisSKDetFSI006 ) :
  case ( kNuisSKDetFSI007 ) :
  case ( kNuisSKDetFSI008 ) :
  case ( kNuisSKDetFSI009 ) :
  case ( kNuisSKDetFSI010 ) :
  case ( kNuisSKDetFSI011 ) :
  case ( kNuisSKDetFSI012 ) :
  case ( kNuisSKDetFSI013 ) :
  case ( kNuisSKDetFSI014 ) :
  case ( kNuisSKDetFSI015 ) :
  case ( kNuisSKDetFSI016 ) :
  case ( kNuisSKDetFSI017 ) :
#endif
#ifdef SKCORRELATIONS
    case ( kNuisSKDetFSI018 ) :
    case ( kNuisSKDetFSI019 ) :
    case ( kNuisSKDetFSI020 ) :
    case ( kNuisSKDetFSI021 ) :
    case ( kNuisSKDetFSI022 ) :
    case ( kNuisSKDetFSI023 ) :
    case ( kNuisSKDetFSI024 ) :
    case ( kNuisSKDetFSI025 ) :
    case ( kNuisSKDetFSI026 ) :
    case ( kNuisSKDetFSI027 ) :
    case ( kNuisSKDetFSI028 ) :
    case ( kNuisSKDetFSI029 ) :
    case ( kNuisSKDetFSI030 ) :
    case ( kNuisSKDetFSI031 ) :
#endif
#endif
    {
#if FSI_BINS != 3
      int ip = (int)p - (int)kNuisSKDetFSI000;
#else
      int ip = (int)p - (int)kNuisSKDetFSI006;
#endif
      return Form("f_{%d; r}^{SK+FSI}", ip);
      break;
    }//FSI

  case ( kNuisRHCToFHC ) : return "f^{RHC/FHC}"; break;

  default :
    return "--";
  }
  return "--";
}
//____________________________________________________________________________                                                                                                    
string t2kana::util::AsStringPaperRooTex(PhysParamId_t p)
{
  switch(p) {

  case ( kPhypMixing23 ) :
    {
      NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
      switch(omp) {
      case(kNuMixSinsq2Theta) : return "sin^{2}(2#theta_{23})"; break;
      case(kNuMixSinsqTheta ) : return "sin^{2}(#theta_{23})";   break;
      default                 : return "--";                break;
      }
    }
    break;
  case ( kPhypMixing13 ) :
    {
      NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
      switch(omp) {
      case(kNuMixSinsq2Theta) : return "sin^{2}(2#theta_{13})";  break;
      case(kNuMixSinsqTheta ) : return "sin^{2}(#theta_{13})";    break;
      default                 : return "--";                 break;
      }
    }
    break;
  case ( kPhypMixing12 ) :
    {
      NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
      switch(omp) {
      case(kNuMixSinsq2Theta) : return "sin^{2}(2#theta_{12})";  break;
      case(kNuMixSinsqTheta ) : return "sin^{2}(#theta_{12})";    break;
      default                 : return "--";                 break;
      }
    }
    break;

  case ( kPhypDeltaCP               ) : return "#delta_{CP}";                   break;
  case ( kPhypDm23sq                ) : 
    {
      NuOscDmsqScheme_t omp = Opt::Instance()->gOscDmsqScheme;
      switch(omp) {
        case(kNuDmsqFogliLisi ):	return "|#Deltam^{2}| = |m^{2}_{3} - (m^{2}_{1} + m^{2}_{1})/2|";                  break;
        case(kNuDmsq12andAbs23):	return "|#Deltam^{2}_{32}|";                  break;
        case(kNuDmsq12andAbs31):	return "|#Deltam^{2}_{31}|";                  break;
        default                :        return "--";                                  break;
      }
    }
    break;

  case ( kPhypDm12sq                ) : 
    {
      NuOscDmsqScheme_t omp = Opt::Instance()->gOscDmsqScheme;
      switch(omp) {
      case(kNuDmsqFogliLisi):
	return "#deltam^{2} = #Deltam^{2}_{21}";     break;
      case(kNuDmsq12andAbs23):
      case(kNuDmsq12andAbs31):
	return "#Delta m^{2}_{21}";                  break;
      default                 :
	return "--";                                 break;
      }
    }

  case ( kPhypAlphaDecay            ) : return "#alpha_{decay}";                break;
  case ( kPhypSinsq2ThetaDecay      ) : return "sin^{2}(2#theta_{decay})";       break;
  case ( kPhypMusqDecoh             ) : return "#mu^{2}_{decoh}";                 break;
  case ( kPhypSinsq2ThetaDecoh      ) : return "sin^{2}(2#theta_{decoh})";       break;

    /*                                   
     * 2011-analysis nuisance parameters
     */

    // 
    // SK detector 
    //                                                                                                                                                                                 
  case ( kNuisSuperkEnergyScale     ) : return "f_{E; r}^{SK}  ";                 break;
#ifndef SKFSITOGETHER
#if FSI_BINS != 3
  case ( kNuisSuperkEff1RmuCCQEoth  ) : return "f_{1R_{#mu}: CCQEoth; r}^{SK}";             break;
  case ( kNuisSuperkEff1RmuCCQErc0u ) : return "f_{1R_{#mu}: CCQERC0; r}^{SK}";             break;
  case ( kNuisSuperkEff1RmuCCQErc1u ) : return "f_{1R_{#mu}: CCQERC1; r}^{SK}";             break;
  case ( kNuisSuperkEff1RmuCCQErc2u ) : return "f_{1R_{#mu}: CCQERC2; r}^{SK}";             break;
  case ( kNuisSuperkEff1RmuCCnonQE  ) : return "f_{1R_{#mu}: CCnQE; r}^{SK}";               break;
  case ( kNuisSuperkEff1RmuNC       ) : return "f_{1R_{#mu}: NC; r}^{SK}";                  break;
  case ( kNuisSuperkEff1RmuNueCC    ) : return "f_{1R_{#mu}: #nu_e; r}^{SK}";              break;
#endif
#if FSI_BINS != 2
  case ( kNuisSuperkEff1ReNueSigCC0  )   : return "f_{1R_{e}: osc #nu_e CC E0; r}^{SK}";     break;
  case ( kNuisSuperkEff1ReNueSigCC1  )   : return "f_{1R_{e}: osc #nu_e CC E1; r}^{SK}";     break;
  case ( kNuisSuperkEff1ReNueSigCC2  )   : return "f_{1R_{e}: osc #nu_e CC E2; r}^{SK}";     break;
  case ( kNuisSuperkEff1ReNumuCC0    )   : return "f_{1R_{e}: #nu_{#mu} CC E0; r}^{SK}";     break;
  case ( kNuisSuperkEff1ReNumuCC1    )   : return "f_{1R_{e}: #nu_{#mu} CC E1; r}^{SK}";   break;
  case ( kNuisSuperkEff1ReNumuCC2    )   : return "f_{1R_{e}: #nu_{#mu} CC E2; r}^{SK}";   break;
  case ( kNuisSuperkEff1ReNueCC0     )   : return "f_{1R_{e}: #nu_{e} CC E0; r}^{SK}";   break;
  case ( kNuisSuperkEff1ReNueCC1     )   : return "f_{1R_{e}: #nu_{e} CC E1; r}^{SK}";   break;
  case ( kNuisSuperkEff1ReNueCC2     )   : return "f_{1R_{e}: #nu_{e} CC E2; r}^{SK}";   break;
  case ( kNuisSuperkEff1ReNC0        )   : return "f_{1R_{e}: NC E0; r}^{SK}";   break;
  case ( kNuisSuperkEff1ReNC1        )   : return "f_{1R_{e}: NC E1; r}^{SK}";   break;
  case ( kNuisSuperkEff1ReNC2        )   : return "f_{1R_{e}: NC E2; r}^{SK}";   break;
#endif
#endif
    //
    // BANFF flux and cross-section params 
    //                                                                                                                                                                                 
  case ( kNuisBANFF00 ) :
  case ( kNuisBANFF01 ) :
  case ( kNuisBANFF02 ) :
  case ( kNuisBANFF03 ) :
  case ( kNuisBANFF04 ) :
  case ( kNuisBANFF05 ) :
  case ( kNuisBANFF06 ) :
  case ( kNuisBANFF07 ) :
  case ( kNuisBANFF08 ) :
  case ( kNuisBANFF09 ) :
  case ( kNuisBANFF10 ) :
  case ( kNuisBANFF11 ) :
  case ( kNuisBANFF12 ) :
  case ( kNuisBANFF13 ) :
  case ( kNuisBANFF14 ) :
  case ( kNuisBANFF15 ) :
#if (BANFF_BINS != 2) && (BANFF_BINS != 4)
  case ( kNuisBANFF16 ) :
  case ( kNuisBANFF17 ) :
  case ( kNuisBANFF18 ) :
  case ( kNuisBANFF19 ) :
  case ( kNuisBANFF20 ) :
  case ( kNuisBANFF21 ) :
#endif
#if BANFF_BINS == 1
  case ( kNuisBANFF22 ) :
  case ( kNuisBANFF23 ) :
  case ( kNuisBANFF24 ) :
#endif
  case ( kNuisBANFFMAQE ) :
  case ( kNuisBANFFMARES ) :
  case ( kNuisBANFFCCQEE0 ) :
#if BANFF_BINS != 3
  case ( kNuisBANFFCCQEE1 ) :
  case ( kNuisBANFFCCQEE2 ) :
#endif
  case ( kNuisBANFFCC1piE0 ) :
#if BANFF_BINS != 3
  case ( kNuisBANFFCC1piE1 ) :
#endif
#if (BANFF_BINS != 2) && (BANFF_BINS != 4)
  case ( kNuisBANFFNC1pi0 ) :
#endif
#if BANFF_BINS == 4
    case ( kNuisBANFFMECNorm ) :
    case ( kNuisBANFFPDDDoubleCount ) :
#endif
    {
      int ip = (int)p - (int)kNuisBANFF00;
      return Form("f_{%d; t, r}^{banff}", ip);
      break;
    }

    // 
    // XSec params not correlated with flux 
    //  

  case ( kNuisXSecFermiMom   ) : return "f_{pF; t, r}";                  break;
  case ( kNuisXSecSpecFunc   ) : return "f_{SF; t, r}";                  break;
  case ( kNuisXSecWshape     ) : return "f_{Wshape; t, r}";              break;
  case ( kNuisXSecShapeCCoth ) : return "f_{CCothShape; t, r}";          break;
  case ( kNuisXSecPiLessDelta) : return "f_{#pi-less#Delta; t, r}";    break;
#ifndef __NUEONLY__
  case ( kNuisXSecBindE      ) : return "f_{bindE; t, r}";               break;
#endif
  case ( kNuisXSecNormCCcoh  ) : return "f_{CCcoh; t}";                  break;
#ifndef __NUMUONLY__
  case ( kNuisXSecNormNCcoh  ) : return "f_{NCcoh; t}";                    break;
#endif
#ifndef __NUEONLY__
  case ( kNuisXSecNormNC1piC ) : return "f_{NC1#pi#pm; t}";           break;
#endif
  case ( kNuisXSecNormNCoth  ) : return "f_{NCoth; t}";                  break;
  case ( kNuisXSecNormNuEToNuMu  ) : return "f_{CC#nue; t}";           break;
#ifndef __NUEONLY__
  case ( kNuisXSecNormAntiNuToNu ) : return "f_{CC#bar{#nu}; t}";      break;
#endif
#ifndef SKFSITOGETHER
#if FSI_BINS != 3
  case ( kNuisFSI000 ) :
  case ( kNuisFSI001 ) :
  case ( kNuisFSI002 ) :
  case ( kNuisFSI003 ) :
  case ( kNuisFSI004 ) :
  case ( kNuisFSI005 ) :
#endif
#if FSI_BINS != 2
  case ( kNuisFSI006 ) :
  case ( kNuisFSI007 ) :
  case ( kNuisFSI008 ) :
  case ( kNuisFSI009 ) :
  case ( kNuisFSI010 ) :
  case ( kNuisFSI011 ) :
  case ( kNuisFSI012 ) :
  case ( kNuisFSI013 ) :
  case ( kNuisFSI014 ) :
  case ( kNuisFSI015 ) :
  case ( kNuisFSI016 ) :
  case ( kNuisFSI017 ) :
#endif
#endif
    {
#if FSI_BINS != 3
      int ip = (int)p - (int)kNuisFSI000;
#else
      int ip = (int)p - (int)kNuisFSI006;
#endif
      return Form("f_{%d; t, r}^{FSI}", ip);
      break;
    }//FSI      
#ifdef SKFSITOGETHER
#if FSI_BINS != 3
  case ( kNuisSKDetFSI000 ) :
  case ( kNuisSKDetFSI001 ) :
  case ( kNuisSKDetFSI002 ) :
  case ( kNuisSKDetFSI003 ) :
  case ( kNuisSKDetFSI004 ) :
  case ( kNuisSKDetFSI005 ) :
#endif
#if FSI_BINS != 2
  case ( kNuisSKDetFSI006 ) :
  case ( kNuisSKDetFSI007 ) :
  case ( kNuisSKDetFSI008 ) :
  case ( kNuisSKDetFSI009 ) :
  case ( kNuisSKDetFSI010 ) :
  case ( kNuisSKDetFSI011 ) :
  case ( kNuisSKDetFSI012 ) :
  case ( kNuisSKDetFSI013 ) :
  case ( kNuisSKDetFSI014 ) :
  case ( kNuisSKDetFSI015 ) :
  case ( kNuisSKDetFSI016 ) :
  case ( kNuisSKDetFSI017 ) :
#endif
#ifdef SKCORRELATIONS
  case ( kNuisSKDetFSI018 ) :
  case ( kNuisSKDetFSI019 ) :
  case ( kNuisSKDetFSI020 ) :
  case ( kNuisSKDetFSI021 ) :
  case ( kNuisSKDetFSI022 ) :
  case ( kNuisSKDetFSI023 ) :
  case ( kNuisSKDetFSI024 ) :
  case ( kNuisSKDetFSI025 ) :
  case ( kNuisSKDetFSI026 ) :
  case ( kNuisSKDetFSI027 ) :
  case ( kNuisSKDetFSI028 ) :
  case ( kNuisSKDetFSI029 ) :
  case ( kNuisSKDetFSI030 ) :
  case ( kNuisSKDetFSI031 ) :
#endif
#endif
    {
#if FSI_BINS != 3
      int ip = (int)p - (int)kNuisSKDetFSI000;
#else
      int ip = (int)p - (int)kNuisSKDetFSI006;
#endif
      return Form("f_{%d; t, r}^{SK+FSI}", ip);
      break;
    }//FSI

  case ( kNuisRHCToFHC ) : return "f^{RHC/FHC}"; break;

  default :
    return "--";
  }
  return "--";
}
//____________________________________________________________________________
string t2kana::util::AsStringCompact(PhysParamId_t p)
{
  switch(p) {

    case ( kPhypMixing23 ) : 
    {
        NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
        switch(omp) {
           case(kNuMixSinsq2Theta) : return "sinsq2theta23";  break;
           case(kNuMixSinsqTheta ) : return "sinsqtheta23";   break;
           default                 : return "--";             break;
         }
    }
    break;
    case ( kPhypMixing13 ) : 
    {
        NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
        switch(omp) {
           case(kNuMixSinsq2Theta) : return "sinsq2theta13";  break;
           case(kNuMixSinsqTheta ) : return "sinsqtheta13";   break;
           default                 : return "--";             break;
         }
    }
    break;
    case ( kPhypMixing12 ) :      
    {
        NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
        switch(omp) {
           case(kNuMixSinsq2Theta) : return "sinsq2theta12";  break;
           case(kNuMixSinsqTheta ) : return "sinsqtheta12";   break;
           default                 : return "--";             break;
         }
    }
    break;

    case ( kPhypDeltaCP               ) : return "delta_cp";                    break;
    case ( kPhypDm23sq                ) :
      {
	NuOscDmsqScheme_t dmscheme = Opt::Instance()->gOscDmsqScheme;
	switch(dmscheme) {
	  case(kNuDmsq12andAbs23) : return "|dm32sq|";                  break;
	  case(kNuDmsq12andAbs31) : return "|dm31sq|";                  break;
	  case(kNuDmsqFogliLisi ) : return "|Dmsq|";                    break;
   	  default                 : return "--";                        break;  
	}
      }
    break;
    case ( kPhypDm12sq                ) :
      {
	NuOscDmsqScheme_t dmscheme = Opt::Instance()->gOscDmsqScheme;
	switch(dmscheme) {
	  case(kNuDmsq12andAbs23) :
	  case(kNuDmsq12andAbs31) :
	    return "dm12sq";                  break;
	  case(kNuDmsqFogliLisi ) :
	    return "dmsq";                    break;
   	  default                 :
	    return "--";                        break;  
	}
      }
    break;

    case ( kPhypAlphaDecay            ) : return "alpha_decay";                 break;
    case ( kPhypSinsq2ThetaDecay      ) : return "sinsq2theta_decay";           break;
    case ( kPhypMusqDecoh             ) : return "musq_decoh";                  break;
    case ( kPhypSinsq2ThetaDecoh      ) : return "sinsq2theta_decoh";           break;

    /*
     * 2011-analysis nuisance parameters
     */

    //**************************** FORWARD HORN CURRENT *******************************

    case ( kNuisSuperkEnergyScale     ) : return "f_sk_e_scale";                   break;    
#ifndef SKFSITOGETHER
#if FSI_BINS != 3
    case ( kNuisSuperkEff1RmuCCQEoth  ) : return "f_sk_eff_1Rmu_ccqe_oth";         break;
    case ( kNuisSuperkEff1RmuCCQErc0u ) : return "f_sk_eff_1Rmu_ccqe_rc0";         break;
    case ( kNuisSuperkEff1RmuCCQErc1u ) : return "f_sk_eff_1Rmu_ccqe_rc1";         break;
    case ( kNuisSuperkEff1RmuCCQErc2u ) : return "f_sk_eff_1Rmu_ccqe_rc2";         break;
    case ( kNuisSuperkEff1RmuCCnonQE  ) : return "f_sk_eff_1Rmu_ccnonqe";          break;
    case ( kNuisSuperkEff1RmuNC       ) : return "f_sk_eff_1Rmu_nc";               break;
    case ( kNuisSuperkEff1RmuNueCC    ) : return "f_sk_eff_1Rmu_nuecc";            break;
#endif
#if FSI_BINS != 2
    case ( kNuisSuperkEff1ReNueSigCC0 ) : return "sk_1Re_eff_signue_0";            break;
    case ( kNuisSuperkEff1ReNueSigCC1 ) : return "sk_1Re_eff_signue_1";            break;
    case ( kNuisSuperkEff1ReNueSigCC2 ) : return "sk_1Re_eff_signue_2";            break;
    case ( kNuisSuperkEff1ReNumuCC0   ) : return "sk_1Re_eff_numu_0";              break;
    case ( kNuisSuperkEff1ReNumuCC1   ) : return "sk_1Re_eff_numu_1";              break;
    case ( kNuisSuperkEff1ReNumuCC2   ) : return "sk_1Re_eff_numu_2";              break;
    case ( kNuisSuperkEff1ReNueCC0    ) : return "sk_1Re_eff_nue_0";               break;
    case ( kNuisSuperkEff1ReNueCC1    ) : return "sk_1Re_eff_nue_1";               break;
    case ( kNuisSuperkEff1ReNueCC2    ) : return "sk_1Re_eff_nue_2";               break;
    case ( kNuisSuperkEff1ReNC0       ) : return "sk_1Re_eff_nc_0";                break;
    case ( kNuisSuperkEff1ReNC1       ) : return "sk_1Re_eff_nc_1";                break;
    case ( kNuisSuperkEff1ReNC2       ) : return "sk_1Re_eff_nc_2";                break;
#endif
#endif
    case ( kNuisBANFF00               ) : return "f_banff00";      break;
    case ( kNuisBANFF01               ) : return "f_banff01";      break;
    case ( kNuisBANFF02               ) : return "f_banff02";      break;
    case ( kNuisBANFF03               ) : return "f_banff03";      break;
    case ( kNuisBANFF04               ) : return "f_banff04";      break;
    case ( kNuisBANFF05               ) : return "f_banff05";      break;
    case ( kNuisBANFF06               ) : return "f_banff06";      break;
    case ( kNuisBANFF07               ) : return "f_banff07";      break;
    case ( kNuisBANFF08               ) : return "f_banff08";      break;
    case ( kNuisBANFF09               ) : return "f_banff09";      break;
    case ( kNuisBANFF10               ) : return "f_banff10";      break;
    case ( kNuisBANFF11               ) : return "f_banff11";      break;
    case ( kNuisBANFF12               ) : return "f_banff12";      break;
    case ( kNuisBANFF13               ) : return "f_banff13";      break;
    case ( kNuisBANFF14               ) : return "f_banff14";      break;
    case ( kNuisBANFF15               ) : return "f_banff15";      break;
#if (BANFF_BINS != 2) && (BANFF_BINS != 4)
    case ( kNuisBANFF16               ) : return "f_banff16";      break;
    case ( kNuisBANFF17               ) : return "f_banff17";      break;
    case ( kNuisBANFF18               ) : return "f_banff18";      break;
    case ( kNuisBANFF19               ) : return "f_banff19";      break;
    case ( kNuisBANFF20               ) : return "f_banff20";      break;
    case ( kNuisBANFF21               ) : return "f_banff21";      break;
#endif
#if BANFF_BINS == 1
    case ( kNuisBANFF22               ) : return "f_banff22";      break;
    case ( kNuisBANFF23               ) : return "f_banff23";      break;
    case ( kNuisBANFF24               ) : return "f_banff24";      break;
#endif

#ifndef ONEOVERMAQE
    case ( kNuisBANFFMAQE             ) : return "f_banffmaqe";    break;
#else
    case ( kNuisBANFFMAQE             ) : return "f_banff1overmaqe2";  break;
#endif
    case ( kNuisBANFFMARES            ) : return "f_banffmares";   break;
    case ( kNuisBANFFCCQEE0           ) : return "f_banffccqee0";  break;
#if BANFF_BINS != 3
    case ( kNuisBANFFCCQEE1           ) : return "f_banffccqee1";  break;
    case ( kNuisBANFFCCQEE2           ) : return "f_banffccqee2";  break;
#endif
    case ( kNuisBANFFCC1piE0          ) : return "f_banffcc1pie0"; break;
#if BANFF_BINS != 3
    case ( kNuisBANFFCC1piE1          ) : return "f_banffcc1pie1"; break;
#endif
#if (BANFF_BINS != 2) && (BANFF_BINS != 4)
    case ( kNuisBANFFNC1pi0           ) : return "f_banffnc1pi0";  break;
#endif
#if BANFF_BINS == 4
    case ( kNuisBANFFMECNorm          ) : return "f_banffmec";     break;
    case ( kNuisBANFFPDDDoubleCount   ) : return "f_banffpdd";     break;
#endif
    case ( kNuisXSecFermiMom          ) : return "f_pfermi";       break;
    case ( kNuisXSecSpecFunc          ) : return "f_specfunc";     break; 
    case ( kNuisXSecWshape            ) : return "f_wshape";       break;
    case ( kNuisXSecShapeCCoth        ) : return "f_shapeCCoth";   break;
    case ( kNuisXSecPiLessDelta       ) : return "f_pilessdelta";  break;
#ifndef __NUEONLY__
    case ( kNuisXSecBindE             ) : return "f_bindE";        break;
#endif
    case ( kNuisXSecNormCCcoh         ) : return "f_norm_cccoh";   break;
#ifndef __NUEONLY__
    case ( kNuisXSecNormNC1piC        ) : return "f_norm_nc1piC";  break;
#endif
#ifndef __NUMUONLY__
    case ( kNuisXSecNormNCcoh         ) : return "f_norm_nccoh";   break;
#endif
    case ( kNuisXSecNormNCoth         ) : return "f_norm_ncoth";   break;
    case ( kNuisXSecNormNuEToNuMu     ) : return "f_norm_nuetonumu";  break;
#ifndef __NUEONLY__
    case ( kNuisXSecNormAntiNuToNu    ) : return "f_norm_antinutonu"; break;
#endif
#ifndef SKFSITOGETHER
#if FSI_BINS != 3
    case ( kNuisFSI000 ) : return "f_fsi000"; break;
    case ( kNuisFSI001 ) : return "f_fsi001"; break;
    case ( kNuisFSI002 ) : return "f_fsi002"; break;
    case ( kNuisFSI003 ) : return "f_fsi003"; break;
    case ( kNuisFSI004 ) : return "f_fsi004"; break;
    case ( kNuisFSI005 ) : return "f_fsi005"; break;
#endif
#if FSI_BINS != 2
    case ( kNuisFSI006 ) : return "f_fsi1Re000"; break;
    case ( kNuisFSI007 ) : return "f_fsi1Re001"; break;
    case ( kNuisFSI008 ) : return "f_fsi1Re002"; break;
    case ( kNuisFSI009 ) : return "f_fsi1Re003"; break;
    case ( kNuisFSI010 ) : return "f_fsi1Re004"; break;
    case ( kNuisFSI011 ) : return "f_fsi1Re005"; break;
    case ( kNuisFSI012 ) : return "f_fsi1Re006"; break;
    case ( kNuisFSI013 ) : return "f_fsi1Re007"; break;
    case ( kNuisFSI014 ) : return "f_fsi1Re008"; break;
    case ( kNuisFSI015 ) : return "f_fsi1Re009"; break;
    case ( kNuisFSI016 ) : return "f_fsi1Re010"; break;
    case ( kNuisFSI017 ) : return "f_fsi1Re011"; break;
#endif
#endif
#ifdef SKFSITOGETHER
#if FSI_BINS != 3
    case ( kNuisSKDetFSI000 ) : return "f_skdetfsi000"; break;
    case ( kNuisSKDetFSI001 ) : return "f_skdetfsi001"; break;
    case ( kNuisSKDetFSI002 ) : return "f_skdetfsi002"; break;
    case ( kNuisSKDetFSI003 ) : return "f_skdetfsi003"; break;
    case ( kNuisSKDetFSI004 ) : return "f_skdetfsi004"; break;
    case ( kNuisSKDetFSI005 ) : return "f_skdetfsi005"; break;
#endif
#if FSI_BINS != 2
    case ( kNuisSKDetFSI006 ) : return "f_skdetfsi006"; break;
    case ( kNuisSKDetFSI007 ) : return "f_skdetfsi007"; break;
    case ( kNuisSKDetFSI008 ) : return "f_skdetfsi008"; break;
    case ( kNuisSKDetFSI009 ) : return "f_skdetfsi009"; break;
    case ( kNuisSKDetFSI010 ) : return "f_skdetfsi010"; break;
    case ( kNuisSKDetFSI011 ) : return "f_skdetfsi011"; break;
    case ( kNuisSKDetFSI012 ) : return "f_skdetfsi012"; break;
    case ( kNuisSKDetFSI013 ) : return "f_skdetfsi013"; break;
    case ( kNuisSKDetFSI014 ) : return "f_skdetfsi014"; break;
    case ( kNuisSKDetFSI015 ) : return "f_skdetfsi015"; break;
    case ( kNuisSKDetFSI016 ) : return "f_skdetfsi016"; break;
    case ( kNuisSKDetFSI017 ) : return "f_skdetfsi017"; break;
#endif
#ifdef SKCORRELATIONS
  case ( kNuisSKDetFSI018 ) : return "f_skdetfsi018"; break;
  case ( kNuisSKDetFSI019 ) : return "f_skdetfsi019"; break;
  case ( kNuisSKDetFSI020 ) : return "f_skdetfsi020"; break;
  case ( kNuisSKDetFSI021 ) : return "f_skdetfsi021"; break;
  case ( kNuisSKDetFSI022 ) : return "f_skdetfsi022"; break;
  case ( kNuisSKDetFSI023 ) : return "f_skdetfsi023"; break;
  case ( kNuisSKDetFSI024 ) : return "f_skdetfsi024"; break;
  case ( kNuisSKDetFSI025 ) : return "f_skdetfsi025"; break;
  case ( kNuisSKDetFSI026 ) : return "f_skdetfsi026"; break;
  case ( kNuisSKDetFSI027 ) : return "f_skdetfsi027"; break;
  case ( kNuisSKDetFSI028 ) : return "f_skdetfsi028"; break;
  case ( kNuisSKDetFSI029 ) : return "f_skdetfsi029"; break;
  case ( kNuisSKDetFSI030 ) : return "f_skdetfsi030"; break;
  case ( kNuisSKDetFSI031 ) : return "f_skdetfsi031"; break;
#endif
#endif

  case ( kNuisRHCToFHC ) : return "f_rhcfhcnorm"; break;

  //***************************** REVERSE HORN CURRENT *************************

  case ( kRHCNuisSuperkEnergyScale     ) : return "f_sk_e_scale_rhc";                   break;

  case ( kRHCNuisBANFF00               ) : return "f_banff00_rhc";      break;
  case ( kRHCNuisBANFF01               ) : return "f_banff01_rhc";      break;
  case ( kRHCNuisBANFF02               ) : return "f_banff02_rhc";      break;
  case ( kRHCNuisBANFF03               ) : return "f_banff03_rhc";      break;
  case ( kRHCNuisBANFF04               ) : return "f_banff04_rhc";      break;
  case ( kRHCNuisBANFF05               ) : return "f_banff05_rhc";      break;
  case ( kRHCNuisBANFF06               ) : return "f_banff06_rhc";      break;
  case ( kRHCNuisBANFF07               ) : return "f_banff07_rhc";      break;
  case ( kRHCNuisBANFF08               ) : return "f_banff08_rhc";      break;
  case ( kRHCNuisBANFF09               ) : return "f_banff09_rhc";      break;
  case ( kRHCNuisBANFF10               ) : return "f_banff10_rhc";      break;
  case ( kRHCNuisBANFF11               ) : return "f_banff11_rhc";      break;
  case ( kRHCNuisBANFF12               ) : return "f_banff12_rhc";      break;
  case ( kRHCNuisBANFF13               ) : return "f_banff13_rhc";      break;
  case ( kRHCNuisBANFF14               ) : return "f_banff14_rhc";      break;
  case ( kRHCNuisBANFF15               ) : return "f_banff15_rhc";      break;
#if (BANFF_BINS != 2) && (BANFF_BINS != 4)
  case ( kRHCNuisBANFF16               ) : return "f_banff16_rhc";      break;
  case ( kRHCNuisBANFF17               ) : return "f_banff17_rhc";      break;
  case ( kRHCNuisBANFF18               ) : return "f_banff18_rhc";      break;
  case ( kRHCNuisBANFF19               ) : return "f_banff19_rhc";      break;
  case ( kRHCNuisBANFF20               ) : return "f_banff20_rhc";      break;
  case ( kRHCNuisBANFF21               ) : return "f_banff21_rhc";      break;
#endif
#if BANFF_BINS == 1
  case ( kRHCNuisBANFF22               ) : return "f_banff22_rhc";      break;
  case ( kRHCNuisBANFF23               ) : return "f_banff23_rhc";      break;
  case ( kRHCNuisBANFF24               ) : return "f_banff24_rhc";      break;
#endif

#ifndef ONEOVERMAQE
  case ( kRHCNuisBANFFMAQE             ) : return "f_banffmaqe_rhc";    break;
#else
  case ( kRHCNuisBANFFMAQE             ) : return "f_banff1overmaqe2_rhc";  break;
#endif
  case ( kRHCNuisBANFFMARES            ) : return "f_banffmares_rhc";   break;
  case ( kRHCNuisBANFFCCQEE0           ) : return "f_banffccqee0_rhc";  break;
#if BANFF_BINS != 3
  case ( kRHCNuisBANFFCCQEE1           ) : return "f_banffccqee1_rhc";  break;
  case ( kRHCNuisBANFFCCQEE2           ) : return "f_banffccqee2_rhc";  break;
#endif
  case ( kRHCNuisBANFFCC1piE0          ) : return "f_banffcc1pie0_rhc"; break;
#if BANFF_BINS != 3
  case ( kRHCNuisBANFFCC1piE1          ) : return "f_banffcc1pie1_rhc"; break;
#endif
#if (BANFF_BINS != 2) && (BANFF_BINS != 4)
  case ( kRHCNuisBANFFNC1pi0           ) : return "f_banffnc1pi0_rhc";  break;
#endif
#if BANFF_BINS == 4
  case ( kRHCNuisBANFFMECNorm          ) : return "f_banffmec_rhc";     break;
  case ( kRHCNuisBANFFPDDDoubleCount   ) : return "f_banffpdd_rhc";     break;
#endif
  case ( kRHCNuisXSecFermiMom          ) : return "f_pfermi_rhc";       break;
  case ( kRHCNuisXSecSpecFunc          ) : return "f_specfunc_rhc";     break;
  case ( kRHCNuisXSecWshape            ) : return "f_wshape_rhc";       break;
  case ( kRHCNuisXSecShapeCCoth        ) : return "f_shapeCCoth_rhc";   break;
  case ( kRHCNuisXSecPiLessDelta       ) : return "f_pilessdelta_rhc";  break;
#ifndef __NUEONLY__
  case ( kRHCNuisXSecBindE             ) : return "f_bindE_rhc";        break;
#endif
  case ( kRHCNuisXSecNormCCcoh         ) : return "f_norm_cccoh_rhc";   break;
#ifndef __NUEONLY__
  case ( kRHCNuisXSecNormNC1piC        ) : return "f_norm_nc1piC_rhc";  break;
#endif
#ifndef __NUMUONLY__
  case ( kRHCNuisXSecNormNCcoh         ) : return "f_norm_nccoh_rhc";   break;
#endif
  case ( kRHCNuisXSecNormNCoth         ) : return "f_norm_ncoth_rhc";   break;
  case ( kRHCNuisXSecNormNuEToNuMu     ) : return "f_norm_nuetonumu_rhc";  break;
#ifndef __NUEONLY__
  case ( kRHCNuisXSecNormAntiNuToNu    ) : return "f_norm_antinutonu_rhc"; break;
#endif
#ifdef SKFSITOGETHER
#if FSI_BINS != 3
  case ( kRHCNuisSKDetFSI000 ) : return "f_skdetfsi000_rhc"; break;
  case ( kRHCNuisSKDetFSI001 ) : return "f_skdetfsi001_rhc"; break;
  case ( kRHCNuisSKDetFSI002 ) : return "f_skdetfsi002_rhc"; break;
  case ( kRHCNuisSKDetFSI003 ) : return "f_skdetfsi003_rhc"; break;
  case ( kRHCNuisSKDetFSI004 ) : return "f_skdetfsi004_rhc"; break;
  case ( kRHCNuisSKDetFSI005 ) : return "f_skdetfsi005_rhc"; break;
#endif
#if FSI_BINS != 2
  case ( kRHCNuisSKDetFSI006 ) : return "f_skdetfsi006_rhc"; break;
  case ( kRHCNuisSKDetFSI007 ) : return "f_skdetfsi007_rhc"; break;
  case ( kRHCNuisSKDetFSI008 ) : return "f_skdetfsi008_rhc"; break;
  case ( kRHCNuisSKDetFSI009 ) : return "f_skdetfsi009_rhc"; break;
  case ( kRHCNuisSKDetFSI010 ) : return "f_skdetfsi010_rhc"; break;
  case ( kRHCNuisSKDetFSI011 ) : return "f_skdetfsi011_rhc"; break;
  case ( kRHCNuisSKDetFSI012 ) : return "f_skdetfsi012_rhc"; break;
  case ( kRHCNuisSKDetFSI013 ) : return "f_skdetfsi013_rhc"; break;
  case ( kRHCNuisSKDetFSI014 ) : return "f_skdetfsi014_rhc"; break;
  case ( kRHCNuisSKDetFSI015 ) : return "f_skdetfsi015_rhc"; break;
  case ( kRHCNuisSKDetFSI016 ) : return "f_skdetfsi016_rhc"; break;
  case ( kRHCNuisSKDetFSI017 ) : return "f_skdetfsi017_rhc"; break;
#endif
#ifdef SKCORRELATIONS
  case ( kRHCNuisSKDetFSI018 ) : return "f_skdetfsi018_rhc"; break;
  case ( kRHCNuisSKDetFSI019 ) : return "f_skdetfsi019_rhc"; break;
  case ( kRHCNuisSKDetFSI020 ) : return "f_skdetfsi020_rhc"; break;
  case ( kRHCNuisSKDetFSI021 ) : return "f_skdetfsi021_rhc"; break;
  case ( kRHCNuisSKDetFSI022 ) : return "f_skdetfsi022_rhc"; break;
  case ( kRHCNuisSKDetFSI023 ) : return "f_skdetfsi023_rhc"; break;
  case ( kRHCNuisSKDetFSI024 ) : return "f_skdetfsi024_rhc"; break;
  case ( kRHCNuisSKDetFSI025 ) : return "f_skdetfsi025_rhc"; break;
  case ( kRHCNuisSKDetFSI026 ) : return "f_skdetfsi026_rhc"; break;
  case ( kRHCNuisSKDetFSI027 ) : return "f_skdetfsi027_rhc"; break;
  case ( kRHCNuisSKDetFSI028 ) : return "f_skdetfsi028_rhc"; break;
  case ( kRHCNuisSKDetFSI029 ) : return "f_skdetfsi029_rhc"; break;
  case ( kRHCNuisSKDetFSI030 ) : return "f_skdetfsi030_rhc"; break;
  case ( kRHCNuisSKDetFSI031 ) : return "f_skdetfsi031_rhc"; break;
#endif
#endif

    default : 
      return "--";
  }
  return "--";
}
//____________________________________________________________________________
string t2kana::util::AsStringCompact_splinehack(PhysParamId_t p)
{
  switch(p) {

  //***************************** REVERSE HORN CURRENT *************************

  case ( kRHCNuisSuperkEnergyScale     ) : return "f_sk_e_scale";                   break;

  case ( kRHCNuisBANFF00               ) : return "f_banff00";      break;
  case ( kRHCNuisBANFF01               ) : return "f_banff01";      break;
  case ( kRHCNuisBANFF02               ) : return "f_banff02";      break;
  case ( kRHCNuisBANFF03               ) : return "f_banff03";      break;
  case ( kRHCNuisBANFF04               ) : return "f_banff04";      break;
  case ( kRHCNuisBANFF05               ) : return "f_banff05";      break;
  case ( kRHCNuisBANFF06               ) : return "f_banff06";      break;
  case ( kRHCNuisBANFF07               ) : return "f_banff07";      break;
  case ( kRHCNuisBANFF08               ) : return "f_banff08";      break;
  case ( kRHCNuisBANFF09               ) : return "f_banff09";      break;
  case ( kRHCNuisBANFF10               ) : return "f_banff10";      break;
  case ( kRHCNuisBANFF11               ) : return "f_banff11";      break;
  case ( kRHCNuisBANFF12               ) : return "f_banff12";      break;
  case ( kRHCNuisBANFF13               ) : return "f_banff13";      break;
  case ( kRHCNuisBANFF14               ) : return "f_banff14";      break;
  case ( kRHCNuisBANFF15               ) : return "f_banff15";      break;
#if (BANFF_BINS != 2) && (BANFF_BINS != 4)
  case ( kRHCNuisBANFF16               ) : return "f_banff16";      break;
  case ( kRHCNuisBANFF17               ) : return "f_banff17";      break;
  case ( kRHCNuisBANFF18               ) : return "f_banff18";      break;
  case ( kRHCNuisBANFF19               ) : return "f_banff19";      break;
  case ( kRHCNuisBANFF20               ) : return "f_banff20";      break;
  case ( kRHCNuisBANFF21               ) : return "f_banff21";      break;
#endif
#if BANFF_BINS == 1
  case ( kRHCNuisBANFF22               ) : return "f_banff22";      break;
  case ( kRHCNuisBANFF23               ) : return "f_banff23";      break;
  case ( kRHCNuisBANFF24               ) : return "f_banff24";      break;
#endif

#ifndef ONEOVERMAQE
  case ( kRHCNuisBANFFMAQE             ) : return "f_banffmaqe";    break;
#else
  case ( kRHCNuisBANFFMAQE             ) : return "f_banff1overmaqe2";  break;
#endif
  case ( kRHCNuisBANFFMARES            ) : return "f_banffmares";   break;
  case ( kRHCNuisBANFFCCQEE0           ) : return "f_banffccqee0";  break;
#if BANFF_BINS != 3
  case ( kRHCNuisBANFFCCQEE1           ) : return "f_banffccqee1";  break;
  case ( kRHCNuisBANFFCCQEE2           ) : return "f_banffccqee2";  break;
#endif
  case ( kRHCNuisBANFFCC1piE0          ) : return "f_banffcc1pie0"; break;
#if BANFF_BINS != 3
  case ( kRHCNuisBANFFCC1piE1          ) : return "f_banffcc1pie1"; break;
#endif
#if (BANFF_BINS != 2) && (BANFF_BINS != 4)
  case ( kRHCNuisBANFFNC1pi0           ) : return "f_banffnc1pi0";  break;
#endif
#if BANFF_BINS == 4
  case ( kRHCNuisBANFFMECNorm          ) : return "f_banffmec";     break;
  case ( kRHCNuisBANFFPDDDoubleCount   ) : return "f_banffpdd";     break;
#endif
  case ( kRHCNuisXSecFermiMom          ) : return "f_pfermi";       break;
  case ( kRHCNuisXSecSpecFunc          ) : return "f_specfunc";     break;
  case ( kRHCNuisXSecWshape            ) : return "f_wshape";       break;
  case ( kRHCNuisXSecShapeCCoth        ) : return "f_shapeCCoth";   break;
  case ( kRHCNuisXSecPiLessDelta       ) : return "f_pilessdelta";  break;
#ifndef __NUEONLY__
  case ( kRHCNuisXSecBindE             ) : return "f_bindE";        break;
#endif
  case ( kRHCNuisXSecNormCCcoh         ) : return "f_norm_cccoh";   break;
#ifndef __NUEONLY__
  case ( kRHCNuisXSecNormNC1piC        ) : return "f_norm_nc1piC";  break;
#endif
#ifndef __NUMUONLY__
  case ( kRHCNuisXSecNormNCcoh         ) : return "f_norm_nccoh";   break;
#endif
  case ( kRHCNuisXSecNormNCoth         ) : return "f_norm_ncoth";   break;
  case ( kRHCNuisXSecNormNuEToNuMu     ) : return "f_norm_nuetonumu";  break;
#ifndef __NUEONLY__
  case ( kRHCNuisXSecNormAntiNuToNu    ) : return "f_norm_antinutonu"; break;
#endif
#ifdef SKFSITOGETHER
#if FSI_BINS != 3
  case ( kRHCNuisSKDetFSI000 ) : return "f_skdetfsi000"; break;
  case ( kRHCNuisSKDetFSI001 ) : return "f_skdetfsi001"; break;
  case ( kRHCNuisSKDetFSI002 ) : return "f_skdetfsi002"; break;
  case ( kRHCNuisSKDetFSI003 ) : return "f_skdetfsi003"; break;
  case ( kRHCNuisSKDetFSI004 ) : return "f_skdetfsi004"; break;
  case ( kRHCNuisSKDetFSI005 ) : return "f_skdetfsi005"; break;
#endif
#if FSI_BINS != 2
  case ( kRHCNuisSKDetFSI006 ) : return "f_skdetfsi006"; break;
  case ( kRHCNuisSKDetFSI007 ) : return "f_skdetfsi007"; break;
  case ( kRHCNuisSKDetFSI008 ) : return "f_skdetfsi008"; break;
  case ( kRHCNuisSKDetFSI009 ) : return "f_skdetfsi009"; break;
  case ( kRHCNuisSKDetFSI010 ) : return "f_skdetfsi010"; break;
  case ( kRHCNuisSKDetFSI011 ) : return "f_skdetfsi011"; break;
  case ( kRHCNuisSKDetFSI012 ) : return "f_skdetfsi012"; break;
  case ( kRHCNuisSKDetFSI013 ) : return "f_skdetfsi013"; break;
  case ( kRHCNuisSKDetFSI014 ) : return "f_skdetfsi014"; break;
  case ( kRHCNuisSKDetFSI015 ) : return "f_skdetfsi015"; break;
  case ( kRHCNuisSKDetFSI016 ) : return "f_skdetfsi016"; break;
  case ( kRHCNuisSKDetFSI017 ) : return "f_skdetfsi017"; break;
#endif
#ifdef SKCORRELATIONS
  case ( kRHCNuisSKDetFSI018 ) : return "f_skdetfsi018"; break;
  case ( kRHCNuisSKDetFSI019 ) : return "f_skdetfsi019"; break;
  case ( kRHCNuisSKDetFSI020 ) : return "f_skdetfsi020"; break;
  case ( kRHCNuisSKDetFSI021 ) : return "f_skdetfsi021"; break;
  case ( kRHCNuisSKDetFSI022 ) : return "f_skdetfsi022"; break;
  case ( kRHCNuisSKDetFSI023 ) : return "f_skdetfsi023"; break;
  case ( kRHCNuisSKDetFSI024 ) : return "f_skdetfsi024"; break;
  case ( kRHCNuisSKDetFSI025 ) : return "f_skdetfsi025"; break;
  case ( kRHCNuisSKDetFSI026 ) : return "f_skdetfsi026"; break;
  case ( kRHCNuisSKDetFSI027 ) : return "f_skdetfsi027"; break;
  case ( kRHCNuisSKDetFSI028 ) : return "f_skdetfsi028"; break;
  case ( kRHCNuisSKDetFSI029 ) : return "f_skdetfsi029"; break;
  case ( kRHCNuisSKDetFSI030 ) : return "f_skdetfsi030"; break;
  case ( kRHCNuisSKDetFSI031 ) : return "f_skdetfsi031"; break;
#endif
#endif

    default : 
      return AsStringCompact(p);
  }
  return AsStringCompact(p);
}
//____________________________________________________________________________
string t2kana::util::AsStringCompact_KDI(PhysParamId_t p)
{
  //Uses the KDI conventions (where different)

  switch(p) {

    case ( kNuisBANFFMAQE             ) : return "maqe";     break;
    case ( kNuisBANFFMARES            ) : return "mares";    break;
    case ( kNuisXSecFermiMom          ) : return "pf";       break;
    case ( kNuisXSecSpecFunc          ) : return "sf";       break; 
    case ( kNuisXSecWshape            ) : return "mdel";     break;
    case ( kNuisXSecShapeCCoth        ) : return "ccothshp"; break;
    case ( kNuisXSecPiLessDelta       ) : return "pdd";      break;
    case ( kNuisXSecBindE             ) : return "eb";       break;
#if BANFF_BINS == 4
    case ( kNuisBANFFPDDDoubleCount   ) : return "pddonly";  break;
#endif

    default : 
      return AsStringCompact(p);
  }
  return "--";
}
//____________________________________________________________________________
string t2kana::util::AsString(t2kana::conv::SuperkSample_t s)
{
  switch(s) {
     case ( kSKsFCFV1RingMuLike ) : return "FCFV 1-ring mu-like"; break;
     case ( kSKsFCFV1RingELike  ) : return "FCFV 1-ring e-like";  break;
     case ( kSKsRHCFCFV1RingMuLike ) : return "FCFV 1-ring mu-like (RHC running)"; break;
     case ( kSKsRHCFCFV1RingELike  ) : return "FCFV 1-ring e-like (RHC running)";  break;

     default                      : return "--";
  }
  return "";
}
//____________________________________________________________________________
string t2kana::util::AsString(t2kana::conv::SuperkMode_t m)
{
  switch(m) {
     case ( kSKmNuMuCCQE        ) : return "numu CC QE";                        break;
     case ( kSKmMEC             ) : return "numu CC MEC";                       break;
     case ( kSKmNuMuCC1pi       ) : return "numu CC 1pi";                       break;
     case ( kSKmNuMuCCcoh       ) : return "numu CC coh";                       break;
     case ( kSKmNuMuCCoth       ) : return "numu CC oth";                       break;
     case ( kSKmNuMuNC1piC      ) : return "numu NC 1piC";                      break;
     case ( kSKmNuMuNC1pi0      ) : return "numu/nutau NC 1pi0";                break;
     case ( kSKmNuMuNCcoh       ) : return "numu/nutau NC coh";                 break;
     case ( kSKmNuMuNCoth_e     ) : return "numu/nutau NC oth (nue app)";       break;
     case ( kSKmNuMuNCoth_mu    ) : return "numu/nutau NC oth (numu disapp)";   break;
     case ( kSKmOscNuMuCCQE     ) : return "oscillated numu CC QE";             break;
     case ( kSKmOscNuMuCC1pi    ) : return "oscillated numu CC 1pi";            break;  
     case ( kSKmOscNuMuCCcoh    ) : return "oscillated numu CC coh";            break;
     case ( kSKmOscNuMuCCoth    ) : return "oscillated numu CC oth";            break;
     case ( kSKmNuMuBarCCQE     ) : return "numubar CC QE";                     break;
     case ( kSKmNuMuBarCC1pi    ) : return "numubar CC 1pi";                    break;
     case ( kSKmNuMuBarCCcoh    ) : return "numubar CC coh";                    break;
     case ( kSKmNuMuBarCCoth    ) : return "numubar CC oth";                    break;
     case ( kSKmNuMuBarNC1piC   ) : return "numubar/nutaubar NC 1piC";          break;
     case ( kSKmNuMuBarNC1pi0   ) : return "numubar/nutaubar NC 1pi0";          break;
     case ( kSKmNuMuBarNCcoh    ) : return "numubar/nutaubar NC coh";           break;
     case ( kSKmNuMuBarNCoth_e  ) : return "numubar/nutaubar NC oth (nue app)"; break;
     case ( kSKmNuMuBarNCoth_mu ) : return "numubar/nutaubar NC oth (numu disapp)"; break;
     case ( kSKmOscNuMuBarCCQE  ) : return "oscillated numubar CC QE";          break;
     case ( kSKmOscNuMuBarCC1pi ) : return "oscillated numubar CC 1pi";         break;  
     case ( kSKmOscNuMuBarCCcoh ) : return "oscillated numubar CC coh";         break;
     case ( kSKmOscNuMuBarCCoth ) : return "oscillated numubar CC oth";         break;
     case ( kSKmNuECCQE         ) : return "nue CC QE";                         break;
     case ( kSKmNuECC1pi        ) : return "nue CC 1pi";                        break;
     case ( kSKmNuECCcoh        ) : return "nue CC coh";                        break;
     case ( kSKmNuECCoth        ) : return "nue CC oth";                        break;
     case ( kSKmNuENC1piC       ) : return "nue NC 1piC";                       break;
     case ( kSKmNuENC1pi0       ) : return "nue NC 1pi0";                       break;
     case ( kSKmNuENCcoh        ) : return "nue NC coh";                        break;
     case ( kSKmNuENCoth_e      ) : return "nue NC oth (nue app)";              break;
     case ( kSKmNuENCoth_mu     ) : return "nue NC oth (numu disapp)";          break;
     case ( kSKmOscNuECCQE      ) : return "oscillated nue CC QE";              break;
     case ( kSKmOscNuECC1pi     ) : return "oscillated nue CC 1pi";             break;
     case ( kSKmOscNuECCcoh     ) : return "oscillated nue CC coh";             break;
     case ( kSKmOscNuECCoth     ) : return "oscillated nue CC oth";             break;
     case ( kSKmNuEBarCCQE      ) : return "nuebar CC QE";                      break;
     case ( kSKmNuEBarCC1pi     ) : return "nuebar CC 1pi";                     break;
     case ( kSKmNuEBarCCcoh     ) : return "nuebar CC coh";                     break;
     case ( kSKmNuEBarCCoth     ) : return "nuebar CC oth";                     break;
     case ( kSKmNuEBarNC1piC    ) : return "nuebar NC 1piC";                    break;
     case ( kSKmNuEBarNC1pi0    ) : return "nuebar NC 1pi0";                    break;
     case ( kSKmNuEBarNCcoh     ) : return "nuebar NC coh";                     break;
     case ( kSKmNuEBarNCoth_e   ) : return "nuebar NC oth (nue app)";           break;
     case ( kSKmNuEBarNCoth_mu  ) : return "nuebar NC oth (numu disapp)";       break;
     case ( kSKmOscNuEBarCCQE   ) : return "oscillated nuebar CC QE";           break;
     case ( kSKmOscNuEBarCC1pi  ) : return "oscillated nuebar CC 1pi";          break;  
     case ( kSKmOscNuEBarCCcoh  ) : return "oscillated nuebar CC coh";          break;
     case ( kSKmOscNuEBarCCoth  ) : return "oscillated nuebar CC oth ";         break;

     default                      : return "--";
  }
  return "";
}
//____________________________________________________________________________
string t2kana::util::AsString(conv::SKMCTemplateType_t t)
{
  switch(t) {
    case ( kSKMCtNom ) : 
      return "nominal";          
      break;
    case ( kSKMCtTun ) : 
      return "tuned"; 
      break;
    case ( kSKMCtSys ) : 
      return "tuned + systematic variation";    
      break;
    case ( kSKMCtFin ) : 
      return "tuned + systematic variation + physics hypothesis";          
      break;
    default : 
      return "--";               
      break;
  }
  return "";
}
//____________________________________________________________________________
string t2kana::util::AsString(t2kana::conv::FluxType_t f)
{
  switch(f) {
    case ( kFlxNominal        ) : return "nominal";          break;
    case ( kFlxExtTuning      ) : return "tuned (external)"; break;
    case ( kFlxFinalTuning    ) : return "tuned (final)";    break;
    case ( kFlxCurrent        ) : return "current";          break;
    case ( kRHCFlxNominal     ) : return "nominal; RHC";          break;
    case ( kRHCFlxExtTuning   ) : return "tuned (external); RHC"; break;
    case ( kRHCFlxFinalTuning ) : return "tuned (final); RHC";    break;
    case ( kRHCFlxCurrent     ) : return "current; RHC";          break;
    default                     : return "--";               break;
  }
  return "";
}
//____________________________________________________________________________
string t2kana::util::AsString(t2kana::conv::OscChannelId_t osch)
{
  switch(osch) {
     case ( kOscNuE2NuE           ) : return "nue -> nue";            break;
     case ( kOscNuE2NuMu          ) : return "nue -> numu";           break;
     case ( kOscNuE2NuTau         ) : return "nue -> nutau";          break;
     case ( kOscNuMu2NuE          ) : return "numu -> nue";           break;  
     case ( kOscNuMu2NuMu         ) : return "numu -> numu";          break;
     case ( kOscNuMu2NuTau        ) : return "numu -> nutau";         break;
     case ( kOscNuTau2NuE         ) : return "nutau -> nue";          break;
     case ( kOscNuTau2NuMu        ) : return "nutau -> numu";         break;
     case ( kOscNuTau2NuTau       ) : return "nutau -> nutau";        break;
     case ( kOscNuEBar2NuEBar     ) : return "nuebar -> nuebar";      break;
     case ( kOscNuEBar2NuMuBar    ) : return "nuebar -> numubar";     break;
     case ( kOscNuEBar2NuTauBar   ) : return "nuebar -> nutaubar";    break;
     case ( kOscNuMuBar2NuEBar    ) : return "numubar -> nuebar";     break;
     case ( kOscNuMuBar2NuMuBar   ) : return "numubar -> numubar";    break;
     case ( kOscNuMuBar2NuTauBar  ) : return "numubar -> nutaubar";   break;
     case ( kOscNuTauBar2NuEBar   ) : return "nutaubar -> nuebar";    break;
     case ( kOscNuTauBar2NuMuBar  ) : return "nutaubar -> numubar";   break;
     case ( kOscNuTauBar2NuTauBar ) : return "nutaubar -> nutaubar";  break;
     default                        : return "--";
  }
  return "--";
}
//____________________________________________________________________________
string t2kana::util::AsStringTex(t2kana::conv::SuperkSample_t s)
{
  switch(s) {
     case ( kSKsFCFV1RingMuLike ) : return "1 $\\mu$-like ring"; break;
     case ( kSKsFCFV1RingELike  ) : return "1 e-like ring";      break;
     case ( kSKsRHCFCFV1RingMuLike ) : return "1 $\\mu$-like ring (RHC running)"; break;
     case ( kSKsRHCFCFV1RingELike  ) : return "1 e-like ring (RHC running)";      break;

     default                      : return "--";
  }
  return "";
}
//____________________________________________________________________________
string t2kana::util::AsStringTex(t2kana::conv::SuperkMode_t m)
{
  switch(m) {

     case ( kSKmNuMuCCQE        ) : return "$\\nu_\\mu$ CCQE";                                               break;
     case ( kSKmMEC             ) : return "$\\nu_\\mu$ CC MEC";                                             break;
     case ( kSKmNuMuCC1pi       ) : return "$\\nu_\\mu$ CC 1$\\pi$";                                         break;
     case ( kSKmNuMuCCcoh       ) : return "$\\nu_\\mu$ CC coherent";                                        break;
     case ( kSKmNuMuCCoth       ) : return "$\\nu_\\mu$ CC other";                                           break;
     case ( kSKmNuMuNC1piC      ) : return "$\\nu_\\mu$/$\\nu_\\tau$ NC 1$\\pi^{+/-}$";                      break;
     case ( kSKmNuMuNC1pi0      ) : return "$\\nu_\\mu$/$\\nu_\\tau$ NC 1$\\pi^0$";                          break;
     case ( kSKmNuMuNCcoh       ) : return "$\\nu_\\mu$/$\\nu_\\tau$ NC coherent";                           break;
     case ( kSKmNuMuNCoth_e     ) : return "$\\nu_\\mu$/$\\nu_\\tau$ NC other ($\\nu_e$ app)";               break;
     case ( kSKmNuMuNCoth_mu    ) : return "$\\nu_\\mu$/$\\nu_\\tau$ NC other ($\\nu_\\mu$ disap)";          break;
     case ( kSKmOscNuMuCCQE     ) : return "Osc. $\\nu_\\mu$ CCQE";                                          break;
     case ( kSKmOscNuMuCC1pi    ) : return "Osc. $\\nu_\\mu$ CC 1$\\pi$";                                    break;
     case ( kSKmOscNuMuCCcoh    ) : return "Osc. $\\nu_\\mu$ CC coherent";                                   break; 
     case ( kSKmOscNuMuCCoth    ) : return "Osc. $\\nu_\\mu$ CC other";                                      break;
     case ( kSKmNuMuBarCCQE     ) : return "$\\bar{\\nu}_\\mu$ CCQE";                                        break;
     case ( kSKmNuMuBarCC1pi    ) : return "$\\bar{\\nu}_\\mu$ CC 1$\\pi$";                                  break;
     case ( kSKmNuMuBarCCcoh    ) : return "$\\bar{\\nu}_\\mu$ CC coherent";                                 break;
     case ( kSKmNuMuBarCCoth    ) : return "$\\bar{\\nu}_\\mu$ CC other";                                    break;
     case ( kSKmNuMuBarNC1piC   ) : return "$\\bar{\\nu}_\\mu$/$\\bar{\\nu}_\\tau$ NC 1$\\pi^{+/-}$";        break;
     case ( kSKmNuMuBarNC1pi0   ) : return "$\\bar{\\nu}_\\mu$/$\\bar{\\nu}_\\tau$ NC 1$\\pi^0$";            break;
     case ( kSKmNuMuBarNCcoh    ) : return "$\\bar{\\nu}_\\mu$/$\\bar{\\nu}_\\tau$ NC coherent";             break;
     case ( kSKmNuMuBarNCoth_e  ) : return "$\\bar{\\nu}_\\mu$/$\\bar{\\nu}_\\tau$ NC other ($\\nu_e$ app)"; break;
     case ( kSKmNuMuBarNCoth_mu ) : return "$\\bar{\\nu}_\\mu$/$\\bar{\\nu}_\\tau$ NC other ($\\nu_\\mu$ disap)"; break;
     case ( kSKmOscNuMuBarCCQE  ) : return "Osc. $\\bar{\\nu}_\\mu$ CCQE";                                   break;
     case ( kSKmOscNuMuBarCC1pi ) : return "Osc. $\\bar{\\nu}_\\mu$ CC 1$\\pi$";                             break;
     case ( kSKmOscNuMuBarCCcoh ) : return "Osc. $\\bar{\\nu}_\\mu$ CC coherent";                            break; 
     case ( kSKmOscNuMuBarCCoth ) : return "Osc. $\\bar{\\nu}_\\mu$ CC other";                               break;
     case ( kSKmNuECCQE         ) : return "$\\nu_e$ CCQE";                                                  break;
     case ( kSKmNuECC1pi        ) : return "$\\nu_e$ CC 1$\\pi$";                                            break;
     case ( kSKmNuECCcoh        ) : return "$\\nu_e$ CC coherent";                                           break;
     case ( kSKmNuECCoth        ) : return "$\\nu_e$ CC other";                                              break;
     case ( kSKmNuENC1piC       ) : return "$\\nu_e$ NC 1$\\pi^{+/-}$";                                      break;
     case ( kSKmNuENC1pi0       ) : return "$\\nu_e$ NC 1$\\pi^0$";                                          break;
     case ( kSKmNuENCcoh        ) : return "$\\nu_e$ NC coherent";                                           break;
     case ( kSKmNuENCoth_e      ) : return "$\\nu_e$ NC other ($\\nu_e$ app)";                               break;
     case ( kSKmNuENCoth_mu     ) : return "$\\nu_e$ NC other ($\\nu_\\mu$ disap)";                          break;
     case ( kSKmOscNuECCQE      ) : return "Osc. $\\nu_e$ CCQE";                                             break;
     case ( kSKmOscNuECC1pi     ) : return "Osc. $\\nu_e$ CC 1$\\pi$";                                       break;
     case ( kSKmOscNuECCcoh     ) : return "Osc. $\\nu_e$ CC coherent";                                      break;
     case ( kSKmOscNuECCoth     ) : return "Osc. $\\nu_e$ CC other";                                         break;
     case ( kSKmNuEBarCCQE      ) : return "$\\bar{\\nu}_e$ CCQE";                                           break;
     case ( kSKmNuEBarCC1pi     ) : return "$\\bar{\\nu}_e$ CC 1$\\pi$";                                     break;
     case ( kSKmNuEBarCCcoh     ) : return "$\\bar{\\nu}_e$ CC coherent";                                    break;
     case ( kSKmNuEBarCCoth     ) : return "$\\bar{\\nu}_e$ CC other";                                       break;
     case ( kSKmNuEBarNC1piC    ) : return "$\\bar{\\nu}_e$ NC 1$\\pi^{+/-}$";                               break;
     case ( kSKmNuEBarNC1pi0    ) : return "$\\bar{\\nu}_e$ NC 1$\\pi^0$";                                   break;
     case ( kSKmNuEBarNCcoh     ) : return "$\\bar{\\nu}_e$ NC coherent";                                    break;
     case ( kSKmNuEBarNCoth_e   ) : return "$\\bar{\\nu}_e$ NC other ($\\nu_e$ app)";                        break;
     case ( kSKmNuEBarNCoth_mu  ) : return "$\\bar{\\nu}_e$ NC other ($\\nu_\\mu$ disap)";                   break;
     case ( kSKmOscNuEBarCCQE   ) : return "Osc. $\\bar{\\nu}_e$ CCQE";                                      break;
     case ( kSKmOscNuEBarCC1pi  ) : return "Osc. $\\bar{\\nu}_e$ CC 1$\\pi$";                                break;
     case ( kSKmOscNuEBarCCcoh  ) : return "Osc. $\\bar{\\nu}_e$ CC coherent";                               break; 
     case ( kSKmOscNuEBarCCoth  ) : return "Osc. $\\bar{\\nu}_e$ CC other";                                  break;

    default                      : return "--";
  }
  return "";
}
//____________________________________________________________________________
string t2kana::util::AsStringRooTex(t2kana::conv::SuperkSample_t s)
{
// As AsStringTex(SuperkSample_t) but taking into account ROOT TeX
// idiosyncracies (this method can be used for building TLegends etc).
//
  switch(s) {
     case ( kSKsFCFV1RingMuLike ) : return "1 #mu-like ring"; break;
     case ( kSKsFCFV1RingELike  ) : return "1 e-like ring";   break;
     case ( kSKsRHCFCFV1RingMuLike ) : return "1 #mu-like ring (RHC running)"; break;
     case ( kSKsRHCFCFV1RingELike  ) : return "1 e-like ring (RHC running)";   break;

     default                      : return "";
  }
  return "";
}
//____________________________________________________________________________
string t2kana::util::AsStringRooTex(t2kana::conv::SuperkMode_t m)
{
  switch(m) {

     case ( kSKmNuMuCCQE        ) : return "#nu_{#mu} CCQE";                                break;
     case ( kSKmMEC             ) : return "#nu_{#mu} CC MEC";                              break;
     case ( kSKmNuMuCC1pi       ) : return "#nu_{#mu} CC 1#pi";                             break;
     case ( kSKmNuMuCCcoh       ) : return "#nu_{#mu} CC coherent";                         break;
     case ( kSKmNuMuCCoth       ) : return "#nu_{#mu} CC other";                            break;
     case ( kSKmNuMuNC1piC      ) : return "#nu_{#mu}/#nu_{#tau} NC 1#pi^{+/-}";            break;
     case ( kSKmNuMuNC1pi0      ) : return "#nu_{#mu}/#nu_{#tau} NC 1#pi^{0}";              break;
     case ( kSKmNuMuNCcoh       ) : return "#nu_{#mu}/#nu_{#tau} NC coherent";              break;
     case ( kSKmNuMuNCoth_e     ) : return "#nu_{#mu}/#nu_{#tau} NC other";                 break;
     case ( kSKmNuMuNCoth_mu    ) : return "#nu_{#mu}/#nu_{#tau} NC other";                 break;
     case ( kSKmOscNuMuCCQE     ) : return "Osc. #nu_{#mu} CCQE";                           break;
     case ( kSKmOscNuMuCC1pi    ) : return "Osc. #nu_{#mu} CC 1#pi";                        break;
     case ( kSKmOscNuMuCCcoh    ) : return "Osc. #nu_{#mu} CC coherent";                    break;
     case ( kSKmOscNuMuCCoth    ) : return "Osc. #nu_{#mu} CC other";                       break;
     case ( kSKmNuMuBarCCQE     ) : return "#bar{#nu_{#mu}} CCQE";                          break;
     case ( kSKmNuMuBarCC1pi    ) : return "#bar{#nu_{#mu}} CC 1#pi";                       break;
     case ( kSKmNuMuBarCCcoh    ) : return "#bar{#nu_{#mu}} CC coherent";                   break;
     case ( kSKmNuMuBarCCoth    ) : return "#bar{#nu_{#mu}} CC other";                      break;
     case ( kSKmNuMuBarNC1piC   ) : return "#bar{#nu}_{#mu}/#bar{#nu}_{#tau} NC 1#pi^{+/-}";break;
     case ( kSKmNuMuBarNC1pi0   ) : return "#bar{#nu}_{#mu}/#bar{#nu}_{#tau} NC 1#pi^{0}";  break;
     case ( kSKmNuMuBarNCcoh    ) : return "#bar{#nu}_{#mu}/#bar{#nu}_{#tau} NC coherent";  break;
     case ( kSKmNuMuBarNCoth_e  ) : return "#bar{#nu}_{#mu}/#bar{#nu}_{#tau} NC other";     break;
     case ( kSKmNuMuBarNCoth_mu ) : return "#bar{#nu}_{#mu}/#bar{#nu}_{#tau} NC other";     break;
     case ( kSKmOscNuMuBarCCQE  ) : return "Osc. #bar{#nu}_{#mu} CCQE";                     break;
     case ( kSKmOscNuMuBarCC1pi ) : return "Osc. #bar{#nu}_{#mu} CC 1#pi";                  break;
     case ( kSKmOscNuMuBarCCcoh ) : return "Osc. #bar{#nu}_{#mu} CC coherent";              break;
     case ( kSKmOscNuMuBarCCoth ) : return "Osc. #bar{#nu_{#mu}} CC other";                 break;
     case ( kSKmNuECCQE         ) : return "#nu_{e} CCQE";                                  break;
     case ( kSKmNuECC1pi        ) : return "#nu_{e} CC 1#pi";                               break;
     case ( kSKmNuECCcoh        ) : return "#nu_{e} CC coherent";                           break;
     case ( kSKmNuECCoth        ) : return "#nu_{e} CC other";                              break;
     case ( kSKmNuENC1piC       ) : return "#nu_{e} NC 1#pi^{+/-}";                         break;
     case ( kSKmNuENC1pi0       ) : return "#nu_{e} NC 1#pi^{0}";                           break;
     case ( kSKmNuENCcoh        ) : return "#nu_{e} NC coherent";                           break;
     case ( kSKmNuENCoth_e      ) : return "#nu_{e} NC other";                              break;
     case ( kSKmNuENCoth_mu     ) : return "#nu_{e} NC other";                              break;
     case ( kSKmOscNuECCQE      ) : return "Osc. #nu_{e} CCQE";                             break;
     case ( kSKmOscNuECC1pi     ) : return "Osc. #nu_{e} CC 1#pi";                          break;
     case ( kSKmOscNuECCcoh     ) : return "Osc. #nu_{e} CC coherent";                      break;
     case ( kSKmOscNuECCoth     ) : return "Osc. #nu_{e} CC other";                         break;
     case ( kSKmNuEBarCCQE      ) : return "#bar{#nu}_{e} CCQE";                            break;
     case ( kSKmNuEBarCC1pi     ) : return "#bar{#nu}_{e} CC 1#pi";                         break;
     case ( kSKmNuEBarCCcoh     ) : return "#bar{#nu}_{e} CC coherent";                     break;
     case ( kSKmNuEBarCCoth     ) : return "#bar{#nu_{e}} CC other";                        break;
     case ( kSKmNuEBarNC1piC    ) : return "#bar{#nu}_{e} NC 1#pi^{+/-}";                   break;
     case ( kSKmNuEBarNC1pi0    ) : return "#bar{#nu}_{e} NC 1#pi^{0}";                     break;
     case ( kSKmNuEBarNCcoh     ) : return "#bar{#nu}_{e} NC coherent";                     break;
     case ( kSKmNuEBarNCoth_e   ) : return "#bar{#nu}_{e} NC other";                        break;
     case ( kSKmNuEBarNCoth_mu  ) : return "#bar{#nu}_{e} NC other";                        break;
     case ( kSKmOscNuEBarCCQE   ) : return "Osc. #bar{#nu}_{e} CCQE";                       break;
     case ( kSKmOscNuEBarCC1pi  ) : return "Osc. #bar{#nu}_{e} CC 1#pi";                    break;
     case ( kSKmOscNuEBarCCcoh  ) : return "Osc. #bar{#nu}_{e} CC coherent";                break;
     case ( kSKmOscNuEBarCCoth  ) : return "Osc. #bar{#nu}_{e} CC other";                   break;

     default                      : return "";
  }
  return "";
}
//____________________________________________________________________________
string t2kana::util::AsStringRooTex (
  t2kana::conv::SuperkSample_t s, t2kana::conv::SuperkMode_t m, 
  string append_begin, string append_end)
{
  ostringstream out;

  if(append_begin.size() > 0) {
    out << append_begin;
  }

  out << t2kana::util::AsStringRooTex(s);

  if(m != kSKmUndefined) {
    out << " (";
    out << t2kana::util::AsStringRooTex(m);
    out << ")";
  }

  if(append_end.size() > 0) {
    out << append_end;
  }

  return out.str();
}
//____________________________________________________________________________
string t2kana::util::AsStringCompact (
  t2kana::conv::SuperkSample_t s, t2kana::conv::SuperkMode_t m, 
  string append_begin, string append_end)
{
// Compact string used for naming output histograms stored ROOT files
//
  ostringstream out;

  if(append_begin.size() > 0) {
    out << append_begin;
  }

  out << t2kana::util::AsStringCompact(s);

  if(m != kSKmUndefined) {
    out << "_";
    out << t2kana::util::AsStringCompact(m);
  }

  if(append_end.size() > 0) {
    out << append_end;
  }

  return out.str();
}
//____________________________________________________________________________
string t2kana::util::AsStringCompact (t2kana::conv::SuperkSample_t s)
{
// Compact string used for naming output histograms stored ROOT files
//
  switch(s) {
     case ( kSKsFCFV1RingMuLike    ) : return "1Rmu";    break;
     case ( kSKsFCFV1RingELike     ) : return "1Re";     break;
     case ( kSKsRHCFCFV1RingMuLike ) : return "RHC1Rmu"; break;
     case ( kSKsRHCFCFV1RingELike  ) : return "RHC1Re";  break;

     default                         : return "";
  }
  return "";
}
//____________________________________________________________________________
string t2kana::util::AsStringCompact_splinehack (t2kana::conv::SuperkSample_t s)
{
// Compact string used for naming output histograms stored ROOT files
//
  switch(s) {
     case ( kSKsRHCFCFV1RingMuLike ) : return "1Rmu"; break;
     case ( kSKsRHCFCFV1RingELike  ) : return "1Re";  break;

     default                         : return AsStringCompact(s);
  }
  return AsStringCompact(s);
}
//____________________________________________________________________________
string t2kana::util::AsStringCompact (t2kana::conv::SuperkMode_t m)
{
// Compact string used for naming output histograms stored ROOT files
//
  switch(m) {
     case ( kSKmNuMuCCQE        ) : return "numuccqe";          break;
     case ( kSKmMEC             ) : return "mec";               break;
     case ( kSKmNuMuCC1pi       ) : return "numucc1pi";         break;
     case ( kSKmNuMuCCcoh       ) : return "numucccoh";         break;
     case ( kSKmNuMuCCoth       ) : return "numuccoth";         break;
     case ( kSKmNuMuNC1piC      ) : return "numunc1piC";        break;
     case ( kSKmNuMuNC1pi0      ) : return "numunc1pi0";        break;
     case ( kSKmNuMuNCcoh       ) : return "numunccoh";         break;
     case ( kSKmNuMuNCoth_e     ) : return "numuncoth";         break;
     case ( kSKmNuMuNCoth_mu    ) : return "numuncoth";         break;
     case ( kSKmOscNuMuCCQE     ) : return "oscnumuccqe";       break;
     case ( kSKmOscNuMuCC1pi    ) : return "oscnumucc1pi";      break;
     case ( kSKmOscNuMuCCcoh    ) : return "oscnumucccoh";      break;
     case ( kSKmOscNuMuCCoth    ) : return "oscnumuccoth";      break;
     case ( kSKmNuMuBarCCQE     ) : return "numubarccqe";       break;
     case ( kSKmNuMuBarCC1pi    ) : return "numubarcc1pi";      break;
     case ( kSKmNuMuBarCCcoh    ) : return "numubarcccoh";      break;
     case ( kSKmNuMuBarCCoth    ) : return "numubarccoth";      break;
     case ( kSKmNuMuBarNC1piC   ) : return "numubarnc1piC";     break;
     case ( kSKmNuMuBarNC1pi0   ) : return "numubarnc1pi0";     break;
     case ( kSKmNuMuBarNCcoh    ) : return "numubarnccoh";      break;
     case ( kSKmNuMuBarNCoth_e  ) : return "numubarncoth";      break;
     case ( kSKmNuMuBarNCoth_mu ) : return "numubarncoth";      break;
     case ( kSKmOscNuMuBarCCQE  ) : return "oscnumubarccqe";    break;
     case ( kSKmOscNuMuBarCC1pi ) : return "oscnumubarcc1pi";   break;
     case ( kSKmOscNuMuBarCCcoh ) : return "oscnumubarcccoh";   break;
     case ( kSKmOscNuMuBarCCoth ) : return "oscnumubarccoth";   break;
     case ( kSKmNuECCQE         ) : return "nueccqe";           break;
     case ( kSKmNuECC1pi        ) : return "nuecc1pi";          break;
     case ( kSKmNuECCcoh        ) : return "nuecccoh";          break;
     case ( kSKmNuECCoth        ) : return "nueccoth";          break;
     case ( kSKmNuENC1piC       ) : return "nuenc1piC";         break;
     case ( kSKmNuENC1pi0       ) : return "nuenc1pi0";         break;
     case ( kSKmNuENCcoh        ) : return "nuenccoh";          break;
     case ( kSKmNuENCoth_e      ) : return "nuencoth";          break;
     case ( kSKmNuENCoth_mu     ) : return "nuencoth";          break;
     case ( kSKmOscNuECCQE      ) : return "oscnueccqe";        break;
     case ( kSKmOscNuECC1pi     ) : return "oscnuecc1pi";       break;
     case ( kSKmOscNuECCcoh     ) : return "oscnuecccoh";       break;
     case ( kSKmOscNuECCoth     ) : return "oscnueccoth";       break;
     case ( kSKmNuEBarCCQE      ) : return "nuebarccqe";        break;
     case ( kSKmNuEBarCC1pi     ) : return "nuebarcc1pi";       break;
     case ( kSKmNuEBarCCcoh     ) : return "nuebarcccoh";       break;
     case ( kSKmNuEBarCCoth     ) : return "nuebarccoth";       break;
     case ( kSKmNuEBarNC1piC    ) : return "nuebarnc1piC";      break;
     case ( kSKmNuEBarNC1pi0    ) : return "nuebarnc1pi0";      break;
     case ( kSKmNuEBarNCcoh     ) : return "nuebarnccoh";       break;
     case ( kSKmNuEBarNCoth_e   ) : return "nuebarncoth";       break;
     case ( kSKmNuEBarNCoth_mu  ) : return "nuebarncoth";       break;
     case ( kSKmOscNuEBarCCQE   ) : return "oscnuebarccqe";     break;
     case ( kSKmOscNuEBarCC1pi  ) : return "oscnuebarcc1pi";    break;
     case ( kSKmOscNuEBarCCcoh  ) : return "oscnuebarcccoh";    break;
     case ( kSKmOscNuEBarCCoth  ) : return "oscnuebarccoth";    break;

     default                      : return "";
  }
  return "";
}
//____________________________________________________________________________
string t2kana::util::AsStringCompact_splinehack (t2kana::conv::SuperkMode_t m)
{
// TEMPORARY HACK using the old naming scheme for osc nue
//
  switch(m) {
     case ( kSKmOscNuECCQE      ) : return "signueccqe";        break;
     case ( kSKmOscNuECC1pi     ) : return "signuecc1pi";       break;
     case ( kSKmOscNuECCcoh     ) : return "signuecccoh";       break;
     case ( kSKmOscNuECCoth     ) : return "signueccoth";       break;

     default                      : return AsStringCompact(m);
  }
  return AsStringCompact(m);
}
//____________________________________________________________________________
string t2kana::util::AsStringCompact_KDI (t2kana::conv::SuperkMode_t m)
{
// Compact string used for naming output histograms stored ROOT files
// using the Japanese analysis naming convention
//
  switch(m) {
     case ( kSKmNuMuCCQE        ) : return "numu_ccqe";       break;
     case ( kSKmMEC             ) : return "mec";             break;
     case ( kSKmNuMuCC1pi       ) : return "numu_ccpi";       break;
     case ( kSKmNuMuCCcoh       ) : return "numu_cccoh";      break;
     case ( kSKmNuMuCCoth       ) : return "numu_ccoth";      break;
     case ( kSKmNuMuNC1piC      ) : return "numu_ncpi";       break;
     case ( kSKmNuMuNC1pi0      ) : return "numu_nc1pi0";     break;
     case ( kSKmNuMuNCcoh       ) : return "numu_nccoh";      break;
     case ( kSKmNuMuNCoth_e     ) : return "numu_ncoth";      break;
     case ( kSKmNuMuNCoth_mu    ) : return "numu_ncoth";      break;
     case ( kSKmOscNuMuCCQE     ) : return "do_not_know_kdi_name_numu_ccqe";       break;
     case ( kSKmOscNuMuCC1pi    ) : return "do_not_know_kdi_name_numu_ccpi";       break;
     case ( kSKmOscNuMuCCcoh    ) : return "do_not_know_kdi_name_numu_cccoh";      break;
     case ( kSKmOscNuMuCCoth    ) : return "do_not_know_kdi_name_numu_ccoth";      break;
     case ( kSKmNuMuBarCCQE     ) : return "numb_ccqe";       break;
     case ( kSKmNuMuBarCC1pi    ) : return "numb_ccpi";       break;
     case ( kSKmNuMuBarCCcoh    ) : return "numb_cccoh";      break;
     case ( kSKmNuMuBarCCoth    ) : return "numb_ccoth";      break;
     case ( kSKmNuMuBarNC1piC   ) : return "numb_ncpi";       break;
     case ( kSKmNuMuBarNC1pi0   ) : return "numb_nc1pi0";     break;
     case ( kSKmNuMuBarNCcoh    ) : return "numb_nccoh";      break;
     case ( kSKmNuMuBarNCoth_e  ) : return "numb_ncoth";      break;
     case ( kSKmNuMuBarNCoth_mu ) : return "numb_ncoth";      break;
     case ( kSKmNuECCQE         ) : return "nue_ccqe";        break;
     case ( kSKmNuECC1pi        ) : return "nue_ccpi";        break;
     case ( kSKmNuECCcoh        ) : return "nue_cccoh";       break;
     case ( kSKmNuECCoth        ) : return "nue_ccoth";       break;
     case ( kSKmNuENC1piC       ) : return "nue_ncpi";        break;
     case ( kSKmNuENC1pi0       ) : return "nue_nc1pi0";      break;
     case ( kSKmNuENCcoh        ) : return "nue_nccoh";       break;
     case ( kSKmNuENCoth_e      ) : return "nue_ncoth";       break;
     case ( kSKmNuENCoth_mu     ) : return "nue_ncoth";       break;
     case ( kSKmNuEBarCCQE      ) : return "nueb_ccqe";       break;
     case ( kSKmNuEBarCC1pi     ) : return "nueb_ccpi";       break;
     case ( kSKmNuEBarCCcoh     ) : return "nueb_cccoh";      break;
     case ( kSKmNuEBarCCoth     ) : return "nueb_ccoth";      break;
     case ( kSKmNuEBarNC1piC    ) : return "nueb_ncpi";       break;
     case ( kSKmNuEBarNC1pi0    ) : return "nueb_nc1pi0";     break;
     case ( kSKmNuEBarNCcoh     ) : return "nueb_nccoh";      break;
     case ( kSKmNuEBarNCoth_e   ) : return "nueb_ncoth";      break;
     case ( kSKmNuEBarNCoth_mu  ) : return "nueb_ncoth";      break;
     case ( kSKmOscNuECCQE      ) : return "signue_ccqe";     break;
     case ( kSKmOscNuECC1pi     ) : return "signue_ccpi";     break;
     case ( kSKmOscNuECCcoh     ) : return "signue_cccoh";    break;
     case ( kSKmOscNuECCoth     ) : return "signue_ccoth";    break;
     default                      : return "";
  }
  return "";
}
//____________________________________________________________________________
string t2kana::util::Pdg2String(int pdg_code)
{
  switch(pdg_code) {
     case ( kPdgNuE      ) : return "nue";      break; 
     case ( kPdgNuMu     ) : return "numu";     break;
     case ( kPdgNuTau    ) : return "nutau";    break;
     case ( kPdgNuEBar   ) : return "nuebar";   break;   
     case ( kPdgNuMuBar  ) : return "numubar";  break;   
     case ( kPdgNuTauBar ) : return "nutaubar"; break;   
     default               : return "";   
  }
  return "";
}
//____________________________________________________________________________
string t2kana::util::Pdg2StringTex(int pdg_code)
{
  switch(pdg_code) {
     case ( kPdgNuE      ) : return "$\\nu_e$";      break; 
     case ( kPdgNuMu     ) : return "$\\nu_\\mu$";     break;
     case ( kPdgNuTau    ) : return "$\\nu_\\tau$";    break;
     case ( kPdgNuEBar   ) : return "$\\bar{\nu_e}$";   break;   
     case ( kPdgNuMuBar  ) : return "$\\bar{\\nu_\\mu}$";  break;   
     case ( kPdgNuTauBar ) : return "$\\bar{\\nu_\\tau}$"; break;   
     default               : return "";   
  }
  return "";
}
//____________________________________________________________________________
t2kana::conv::SuperkMode_t t2kana::util::SuperkModeFromString(string sh)
{
  // WARNING: "NC oth" returns the disappearance mode if there is no disapp/app information given
  if(sh == "numu CC QE"  ) return kSKmNuMuCCQE;
  if(sh == "numu CC MEC" ) return kSKmMEC;
  if(sh == "numu CC 1pi" ) return kSKmNuMuCC1pi;
  if(sh == "numu CC coh" ) return kSKmNuMuCCcoh;
  if(sh == "numu CC oth" ) return kSKmNuMuCCoth;
  if(sh == "numu NC 1piC") return kSKmNuMuNC1piC;
  if(sh == "numu/nutau NC 1pi0") return kSKmNuMuNC1pi0;
  if(sh == "numu/nutau NC coh" ) return kSKmNuMuNCcoh;
  if(sh == "numu/nutau NC oth" ) return kSKmNuMuNCoth_mu;
  if(sh == "numu/nutau NC oth (numu disapp)") return kSKmNuMuNCoth_mu;
  if(sh == "numu/nutau NC oth (nue app)"    ) return kSKmNuMuNCoth_e;
  if(sh == "oscillated numu CC QE" ) return kSKmOscNuMuCCQE;
  if(sh == "oscillated numu CC 1pi") return kSKmOscNuMuCC1pi;
  if(sh == "oscillated numu CC coh") return kSKmOscNuMuCCcoh;
  if(sh == "oscillated numu CC oth") return kSKmOscNuMuCCoth;
  if(sh == "numubar CC QE" ) return kSKmNuMuBarCCQE;
  if(sh == "numubar CC 1pi") return kSKmNuMuBarCC1pi;
  if(sh == "numubar CC coh") return kSKmNuMuBarCCcoh;
  if(sh == "numubar CC oth") return kSKmNuMuBarCCoth;
  if(sh == "numubar/nutaubar NC 1piC") return kSKmNuMuBarNC1piC;
  if(sh == "numubar/nutaubar NC 1pi0") return kSKmNuMuBarNC1pi0;
  if(sh == "numubar/nutaubar NC coh" ) return kSKmNuMuBarNCcoh;
  if(sh == "numubar/nutaubar NC oth" ) return kSKmNuMuBarNCoth_mu;
  if(sh == "numubar/nutaubar NC oth (numu disapp)") return kSKmNuMuBarNCoth_mu;
  if(sh == "numubar/nutaubar NC oth (nue app)"    ) return kSKmNuMuBarNCoth_e;
  if(sh == "oscillated numubar CC QE" ) return kSKmOscNuMuBarCCQE;
  if(sh == "oscillated numubar CC 1pi") return kSKmOscNuMuBarCC1pi;
  if(sh == "oscillated numubar CC coh") return kSKmOscNuMuBarCCcoh;
  if(sh == "oscillated numubar CC oth") return kSKmOscNuMuBarCCoth;
  if(sh == "nue CC QE"  ) return kSKmNuECCQE;
  if(sh == "nue CC 1pi" ) return kSKmNuECC1pi;
  if(sh == "nue CC coh" ) return kSKmNuECCcoh;
  if(sh == "nue CC oth" ) return kSKmNuECCoth;
  if(sh == "nue NC 1piC") return kSKmNuENC1piC;
  if(sh == "nue NC 1pi0") return kSKmNuENC1pi0;
  if(sh == "nue NC coh" ) return kSKmNuENCcoh;
  if(sh == "nue NC oth" ) return kSKmNuENCoth_mu;
  if(sh == "nue NC oth (nue app)"    ) return kSKmNuENCoth_e;
  if(sh == "nue NC oth (numu disapp)") return kSKmNuENCoth_mu;
  if(sh == "oscillated nue CC QE" ) return kSKmOscNuECCQE;
  if(sh == "oscillated nue CC 1pi") return kSKmOscNuECC1pi;
  if(sh == "oscillated nue CC coh") return kSKmOscNuECCcoh;
  if(sh == "oscillated nue CC oth") return kSKmOscNuECCoth;
  if(sh == "nuebar CC QE"  ) return kSKmNuEBarCCQE;
  if(sh == "nuebar CC 1pi" ) return kSKmNuEBarCC1pi;
  if(sh == "nuebar CC coh" ) return kSKmNuEBarCCcoh;
  if(sh == "nuebar CC oth" ) return kSKmNuEBarCCoth;
  if(sh == "nuebar NC 1piC") return kSKmNuEBarNC1piC;
  if(sh == "nuebar NC 1pi0") return kSKmNuEBarNC1pi0;
  if(sh == "nuebar NC coh" ) return kSKmNuEBarNCcoh;
  if(sh == "nuebar NC oth" ) return kSKmNuEBarNCoth_mu;
  if(sh == "nuebar NC oth (numu disapp)") return kSKmNuEBarNCoth_mu;
  if(sh == "nuebar NC oth (nue app)"    ) return kSKmNuEBarNCoth_e;
  if(sh == "oscillated nuebar CC QE" ) return kSKmOscNuEBarCCQE;
  if(sh == "oscillated nuebar CC 1pi") return kSKmOscNuEBarCC1pi;
  if(sh == "oscillated nuebar CC coh") return kSKmOscNuEBarCCcoh;
  if(sh == "oscillated nuebar CC oth") return kSKmOscNuEBarCCoth;

  return kSKmUndefined;
}
//____________________________________________________________________________
t2kana::conv::PhysHypothesis_t t2kana::util::PhysHypFromString(string sh)
{
  if(sh == "decay"         ) return kHypDecay;
  if(sh == "decoh"         ) return kHypDecoh;
  if(sh == "2f_numu_disap" ) return kHypNuMuDisapp2F;
  if(sh == "3f_osc"        ) return kHypOsc3F;

  return kHypUndefined;
}
//____________________________________________________________________________
t2kana::conv::SuperkSample_t t2kana::util::SKSampleFromString(string sh)
{
  if(sh == "1Rmu") return kSKsFCFV1RingMuLike;
  if(sh == "1Re" ) return kSKsFCFV1RingELike;
  if(sh == "RHC1Rmu") return kSKsRHCFCFV1RingMuLike;
  if(sh == "RHC1Re" ) return kSKsRHCFCFV1RingELike;

  return kSKsUndefined;
}
//____________________________________________________________________________
t2kana::conv::OscAnaType_t  
   t2kana::util::OscAnaFromString(string sh) 
{
  if(sh == "combined"   ) return kAnaCombined;
  if(sh == "numu_disap" ) return kAnaNuMuDisappOnly;
  if(sh == "nue_app"    ) return kAnaNuEAppOnly;
  
  return kAnaUndefined;
}
//____________________________________________________________________________
t2kana::conv::NuOscDmsqScheme_t  
  t2kana::util::DmsqSchemeFromString(string sh) 
{
  if(sh == "abs23"    ) return kNuDmsq12andAbs23;
  if(sh == "abs31"    ) return kNuDmsq12andAbs31;
  if(sh == "fogli"    ) return kNuDmsqFogliLisi;
  
  return kNuDmsqUndefined;
}
//____________________________________________________________________________
t2kana::conv::NuMassHierarchy_t  
  t2kana::util::HierarchyFromString(string sh) 
{
  if(sh == "normal"   ) return kNmhNormal;
  if(sh == "inverted" ) return kNmhInverted;
  
  return kNmhUndefined;
}
//____________________________________________________________________________
t2kana::conv::NuOscOctant_t   
  t2kana::util::OctantFromString(string sh) 
{
  if(sh == "one") return kNuOctantOne;
  if(sh == "two") return kNuOctantTwo;
  
  return kNuOctantUndefined;
}
//____________________________________________________________________________
t2kana::conv::NuOscMixParam_t   
  t2kana::util::MixingTypeFromString(string sh) 
{
  if(sh == "sinsq2theta") return kNuMixSinsq2Theta;
  if(sh == "sinsqtheta" ) return kNuMixSinsqTheta;
  
  return kNuMixUndefined;
}
//____________________________________________________________________________
t2kana::conv::PhysParamId_t
  t2kana::util::PhysParamIdFromString(string sh) 
{

  //************************* FORWARD HORN CURRENT ***************************

  if(sh == "mixing23") return kPhypMixing23;
  if(sh == "mixing13") return kPhypMixing13;
  if(sh == "mixing12") return kPhypMixing12;
  if(sh == "deltacp") return kPhypDeltaCP;
  if(sh == "dm23sq") return kPhypDm23sq;
  if(sh == "dm31sq") return kPhypDm23sq;
  if(sh == "Dmsq")   return kPhypDm23sq;
  if(sh == "|dm23sq|") return kPhypDm23sq;
  if(sh == "|dm31sq|") return kPhypDm23sq;
  if(sh == "|Dmsq|")   return kPhypDm23sq;
  if(sh == "dm12sq") return kPhypDm12sq;
  if(sh == "alpha_decay") return kPhypAlphaDecay;
  if(sh == "sinsq2theta_decay") return kPhypSinsq2ThetaDecay;
  if(sh == "musq_decoh") return kPhypMusqDecoh;
  if(sh == "sinsq2theta_decoh") return kPhypSinsq2ThetaDecoh;

  if(sh == "sk_energy_scale") return kNuisSuperkEnergyScale;

  if(sh == "banff_00") return kNuisBANFF00;
  if(sh == "banff_01") return kNuisBANFF01;
  if(sh == "banff_02") return kNuisBANFF02;
  if(sh == "banff_03") return kNuisBANFF03;
  if(sh == "banff_04") return kNuisBANFF04;
  if(sh == "banff_05") return kNuisBANFF05;
  if(sh == "banff_06") return kNuisBANFF06;
  if(sh == "banff_07") return kNuisBANFF07;
  if(sh == "banff_08") return kNuisBANFF08;
  if(sh == "banff_09") return kNuisBANFF09;
  if(sh == "banff_10") return kNuisBANFF10;
  if(sh == "banff_11") return kNuisBANFF11;
  if(sh == "banff_12") return kNuisBANFF12;
  if(sh == "banff_13") return kNuisBANFF13;
  if(sh == "banff_14") return kNuisBANFF14;
  if(sh == "banff_15") return kNuisBANFF15;
#if (BANFF_BINS != 2) && (BANFF_BINS != 4)
  if(sh == "banff_16") return kNuisBANFF16;
  if(sh == "banff_17") return kNuisBANFF17;
  if(sh == "banff_18") return kNuisBANFF18;
  if(sh == "banff_19") return kNuisBANFF19;
  if(sh == "banff_20") return kNuisBANFF20;
  if(sh == "banff_21") return kNuisBANFF21;
#endif
#if BANFF_BINS == 1
  if(sh == "banff_22") return kNuisBANFF22;
  if(sh == "banff_23") return kNuisBANFF23;
  if(sh == "banff_24") return kNuisBANFF24;
#endif

  if(sh == "banff_maqe")    return kNuisBANFFMAQE; 
  if(sh == "banff_mares")   return kNuisBANFFMARES;
  if(sh == "banff_ccqeE0")  return kNuisBANFFCCQEE0;
#if BANFF_BINS != 3
  if(sh == "banff_ccqeE1")  return kNuisBANFFCCQEE1;
  if(sh == "banff_ccqeE2")  return kNuisBANFFCCQEE2;
#endif
  if(sh == "banff_cc1piE0") return kNuisBANFFCC1piE0;
#if BANFF_BINS != 3
  if(sh == "banff_cc1piE1") return kNuisBANFFCC1piE1;
#endif
#if (BANFF_BINS != 2) && (BANFF_BINS != 4)
  if(sh == "banff_nc1pi0") return kNuisBANFFNC1pi0;
#endif
#if BANFF_BINS == 4
  if(sh == "banff_mecnorm")  return kNuisBANFFMECNorm;
  if(sh == "banff_pdd")      return kNuisBANFFPDDDoubleCount;
#endif

#if BANFF_BINS == 1
  if(sh == "banff_25")    return kNuisBANFFMAQE; 
  if(sh == "banff_26")   return kNuisBANFFMARES;
  if(sh == "banff_27")  return kNuisBANFFCCQEE0;
  if(sh == "banff_28")  return kNuisBANFFCCQEE1;
  if(sh == "banff_29")  return kNuisBANFFCCQEE2;
  if(sh == "banff_30") return kNuisBANFFCC1piE0;
  if(sh == "banff_31") return kNuisBANFFCC1piE1;
  if(sh == "banff_32") return kNuisBANFFNC1pi0;
#elif BANFF_BINS != 3
  if(sh == "banff_16")    return kNuisBANFFMAQE; 
  if(sh == "banff_17")   return kNuisBANFFMARES;
  if(sh == "banff_18")  return kNuisBANFFCCQEE0;
  if(sh == "banff_19")  return kNuisBANFFCCQEE1;
  if(sh == "banff_20")  return kNuisBANFFCCQEE2;
  if(sh == "banff_21") return kNuisBANFFCC1piE0;
  if(sh == "banff_22") return kNuisBANFFCC1piE1;
#elif BANFF_BINS == 3
  if(sh == "banff_22")    return kNuisBANFFMAQE;
  if(sh == "banff_23")   return kNuisBANFFMARES;
  if(sh == "banff_24")  return kNuisBANFFCCQEE0;
  if(sh == "banff_25") return kNuisBANFFCC1piE0;
  if(sh == "banff_26") return kNuisBANFFNC1pi0;
#endif
#if BANFF_BINS == 4
  if(sh == "banff_23")  return kNuisBANFFMECNorm;
  if(sh == "banff_24")  return kNuisBANFFPDDDoubleCount;
#endif

  if(sh == "xsec_pf")           return kNuisXSecFermiMom;
  if(sh == "xsec_sf")           return kNuisXSecSpecFunc;
  if(sh == "xsec_wshape")       return kNuisXSecWshape;
  if(sh == "xsec_shape_ccoth")  return kNuisXSecShapeCCoth;
  if(sh == "xsec_pilessdelta")  return kNuisXSecPiLessDelta;
#ifndef __NUEONLY__
  if(sh == "xsec_bindE")        return kNuisXSecBindE;
#endif
  if(sh == "xsec_norm_cccoh")   return kNuisXSecNormCCcoh;
#ifndef __NUEONLY__
  if(sh == "xsec_norm_nc1piC")  return kNuisXSecNormNC1piC;
#endif
#ifndef __NUMUONLY__
  if(sh == "xsec_norm_nccoh")   return kNuisXSecNormNCcoh;
#endif
  if(sh == "xsec_norm_ncoth")   return kNuisXSecNormNCoth;
  if(sh == "xsec_norm_nuetonumu") return kNuisXSecNormNuEToNuMu;
#ifndef __NUEONLY__
  if(sh == "xsec_norm_antinutonu") return kNuisXSecNormAntiNuToNu;
#endif

#if FSI_BINS != 3
  if(sh == "skdetfsi_000") return kNuisSKDetFSI000;
  if(sh == "skdetfsi_001") return kNuisSKDetFSI001;
  if(sh == "skdetfsi_002") return kNuisSKDetFSI002;
  if(sh == "skdetfsi_003") return kNuisSKDetFSI003;
  if(sh == "skdetfsi_004") return kNuisSKDetFSI004;
  if(sh == "skdetfsi_005") return kNuisSKDetFSI005;
#endif
#if FSI_BINS != 2
  if(sh == "skdetfsi_006") return kNuisSKDetFSI006;
  if(sh == "skdetfsi_007") return kNuisSKDetFSI007;
  if(sh == "skdetfsi_008") return kNuisSKDetFSI008;
  if(sh == "skdetfsi_009") return kNuisSKDetFSI009;
  if(sh == "skdetfsi_010") return kNuisSKDetFSI010;
  if(sh == "skdetfsi_011") return kNuisSKDetFSI011;
  if(sh == "skdetfsi_012") return kNuisSKDetFSI012;
  if(sh == "skdetfsi_013") return kNuisSKDetFSI013;
  if(sh == "skdetfsi_014") return kNuisSKDetFSI014;
  if(sh == "skdetfsi_015") return kNuisSKDetFSI015;
  if(sh == "skdetfsi_016") return kNuisSKDetFSI016;
  if(sh == "skdetfsi_017") return kNuisSKDetFSI017;
#endif
#ifdef SKCORRELATIONS
  if(sh == "skdetfsi_018") return kNuisSKDetFSI018;
  if(sh == "skdetfsi_019") return kNuisSKDetFSI019;
  if(sh == "skdetfsi_020") return kNuisSKDetFSI020;
  if(sh == "skdetfsi_021") return kNuisSKDetFSI021;
  if(sh == "skdetfsi_022") return kNuisSKDetFSI022;
  if(sh == "skdetfsi_023") return kNuisSKDetFSI023;
  if(sh == "skdetfsi_024") return kNuisSKDetFSI024;
  if(sh == "skdetfsi_025") return kNuisSKDetFSI025;
  if(sh == "skdetfsi_026") return kNuisSKDetFSI026;
  if(sh == "skdetfsi_027") return kNuisSKDetFSI027;
  if(sh == "skdetfsi_028") return kNuisSKDetFSI028;
  if(sh == "skdetfsi_029") return kNuisSKDetFSI029;
  if(sh == "skdetfsi_030") return kNuisSKDetFSI030;
  if(sh == "skdetfsi_031") return kNuisSKDetFSI031;
#endif

  if(sh == "rhcfhcnorm") return kNuisRHCToFHC;

  //************************ REVERSE HORN CURRENT *****************************

  if(sh == "sk_energy_scale_rhc") return kRHCNuisSuperkEnergyScale;

  if(sh == "banff_00_rhc") return kRHCNuisBANFF00;
  if(sh == "banff_01_rhc") return kRHCNuisBANFF01;
  if(sh == "banff_02_rhc") return kRHCNuisBANFF02;
  if(sh == "banff_03_rhc") return kRHCNuisBANFF03;
  if(sh == "banff_04_rhc") return kRHCNuisBANFF04;
  if(sh == "banff_05_rhc") return kRHCNuisBANFF05;
  if(sh == "banff_06_rhc") return kRHCNuisBANFF06;
  if(sh == "banff_07_rhc") return kRHCNuisBANFF07;
  if(sh == "banff_08_rhc") return kRHCNuisBANFF08;
  if(sh == "banff_09_rhc") return kRHCNuisBANFF09;
  if(sh == "banff_10_rhc") return kRHCNuisBANFF10;
  if(sh == "banff_11_rhc") return kRHCNuisBANFF11;
  if(sh == "banff_12_rhc") return kRHCNuisBANFF12;
  if(sh == "banff_13_rhc") return kRHCNuisBANFF13;
  if(sh == "banff_14_rhc") return kRHCNuisBANFF14;
  if(sh == "banff_15_rhc") return kRHCNuisBANFF15;
#if (BANFF_BINS != 2) && (BANFF_BINS != 4)
  if(sh == "banff_16_rhc") return kRHCNuisBANFF16;
  if(sh == "banff_17_rhc") return kRHCNuisBANFF17;
  if(sh == "banff_18_rhc") return kRHCNuisBANFF18;
  if(sh == "banff_19_rhc") return kRHCNuisBANFF19;
  if(sh == "banff_20_rhc") return kRHCNuisBANFF20;
  if(sh == "banff_21_rhc") return kRHCNuisBANFF21;
#endif
#if BANFF_BINS == 1
  if(sh == "banff_22_rhc") return kRHCNuisBANFF22;
  if(sh == "banff_23_rhc") return kRHCNuisBANFF23;
  if(sh == "banff_24_rhc") return kRHCNuisBANFF24;
#endif

  if(sh == "banff_maqe_rhc")    return kRHCNuisBANFFMAQE;
  if(sh == "banff_mares_rhc")   return kRHCNuisBANFFMARES;
  if(sh == "banff_ccqeE0_rhc")  return kRHCNuisBANFFCCQEE0;
#if BANFF_BINS != 3
  if(sh == "banff_ccqeE1_rhc")  return kRHCNuisBANFFCCQEE1;
  if(sh == "banff_ccqeE2_rhc")  return kRHCNuisBANFFCCQEE2;
#endif
  if(sh == "banff_cc1piE0_rhc") return kRHCNuisBANFFCC1piE0;
#if BANFF_BINS != 3
  if(sh == "banff_cc1piE1_rhc") return kRHCNuisBANFFCC1piE1;
#endif
#if (BANFF_BINS != 2) && (BANFF_BINS != 4)
  if(sh == "banff_nc1pi0_rhc") return kRHCNuisBANFFNC1pi0;
#endif
#if BANFF_BINS == 4
  if(sh == "banff_mecnorm_rhc")  return kRHCNuisBANFFMECNorm;
  if(sh == "banff_pdd_rhc")      return kRHCNuisBANFFPDDDoubleCount;
#endif

#if BANFF_BINS == 1
  if(sh == "banff_25_rhc")    return kRHCNuisBANFFMAQE;
  if(sh == "banff_26_rhc")   return kRHCNuisBANFFMARES;
  if(sh == "banff_27_rhc")  return kRHCNuisBANFFCCQEE0;
  if(sh == "banff_28_rhc")  return kRHCNuisBANFFCCQEE1;
  if(sh == "banff_29_rhc")  return kRHCNuisBANFFCCQEE2;
  if(sh == "banff_30_rhc") return kRHCNuisBANFFCC1piE0;
  if(sh == "banff_31_rhc") return kRHCNuisBANFFCC1piE1;
  if(sh == "banff_32_rhc") return kRHCNuisBANFFNC1pi0;
#elif BANFF_BINS != 3
  if(sh == "banff_16_rhc")    return kRHCNuisBANFFMAQE;
  if(sh == "banff_17_rhc")   return kRHCNuisBANFFMARES;
  if(sh == "banff_18_rhc")  return kRHCNuisBANFFCCQEE0;
  if(sh == "banff_19_rhc")  return kRHCNuisBANFFCCQEE1;
  if(sh == "banff_20_rhc")  return kRHCNuisBANFFCCQEE2;
  if(sh == "banff_21_rhc") return kRHCNuisBANFFCC1piE0;
  if(sh == "banff_22_rhc") return kRHCNuisBANFFCC1piE1;
#elif BANFF_BINS == 3
  if(sh == "banff_22_rhc")    return kRHCNuisBANFFMAQE;
  if(sh == "banff_23_rhc")   return kRHCNuisBANFFMARES;
  if(sh == "banff_24_rhc")  return kRHCNuisBANFFCCQEE0;
  if(sh == "banff_25_rhc") return kRHCNuisBANFFCC1piE0;
  if(sh == "banff_26_rhc") return kRHCNuisBANFFNC1pi0;
#endif
#if BANFF_BINS == 4
  if(sh == "banff_23_rhc")  return kRHCNuisBANFFMECNorm;
  if(sh == "banff_24_rhc")  return kRHCNuisBANFFPDDDoubleCount;
#endif

  if(sh == "xsec_pf_rhc")           return kRHCNuisXSecFermiMom;
  if(sh == "xsec_sf_rhc")           return kRHCNuisXSecSpecFunc;
  if(sh == "xsec_wshape_rhc")       return kRHCNuisXSecWshape;
  if(sh == "xsec_shape_ccoth_rhc")  return kRHCNuisXSecShapeCCoth;
  if(sh == "xsec_pilessdelta_rhc")  return kRHCNuisXSecPiLessDelta;
#ifndef __NUEONLY__
  if(sh == "xsec_bindE_rhc")        return kRHCNuisXSecBindE;
#endif
  if(sh == "xsec_norm_cccoh_rhc")   return kRHCNuisXSecNormCCcoh;
#ifndef __NUEONLY__
  if(sh == "xsec_norm_nc1piC_rhc")  return kRHCNuisXSecNormNC1piC;
#endif
#ifndef __NUMUONLY__
  if(sh == "xsec_norm_nccoh_rhc")   return kRHCNuisXSecNormNCcoh;
#endif
  if(sh == "xsec_norm_ncoth_rhc")   return kRHCNuisXSecNormNCoth;
  if(sh == "xsec_norm_nuetonumu_rhc") return kRHCNuisXSecNormNuEToNuMu;
#ifndef __NUEONLY__
  if(sh == "xsec_norm_antinutonu_rhc") return kRHCNuisXSecNormAntiNuToNu;
#endif

#if FSI_BINS != 3
  if(sh == "skdetfsi_000_rhc") return kRHCNuisSKDetFSI000;
  if(sh == "skdetfsi_001_rhc") return kRHCNuisSKDetFSI001;
  if(sh == "skdetfsi_002_rhc") return kRHCNuisSKDetFSI002;
  if(sh == "skdetfsi_003_rhc") return kRHCNuisSKDetFSI003;
  if(sh == "skdetfsi_004_rhc") return kRHCNuisSKDetFSI004;
  if(sh == "skdetfsi_005_rhc") return kRHCNuisSKDetFSI005;
#endif
#if FSI_BINS != 2
  if(sh == "skdetfsi_006_rhc") return kRHCNuisSKDetFSI006;
  if(sh == "skdetfsi_007_rhc") return kRHCNuisSKDetFSI007;
  if(sh == "skdetfsi_008_rhc") return kRHCNuisSKDetFSI008;
  if(sh == "skdetfsi_009_rhc") return kRHCNuisSKDetFSI009;
  if(sh == "skdetfsi_010_rhc") return kRHCNuisSKDetFSI010;
  if(sh == "skdetfsi_011_rhc") return kRHCNuisSKDetFSI011;
  if(sh == "skdetfsi_012_rhc") return kRHCNuisSKDetFSI012;
  if(sh == "skdetfsi_013_rhc") return kRHCNuisSKDetFSI013;
  if(sh == "skdetfsi_014_rhc") return kRHCNuisSKDetFSI014;
  if(sh == "skdetfsi_015_rhc") return kRHCNuisSKDetFSI015;
  if(sh == "skdetfsi_016_rhc") return kRHCNuisSKDetFSI016;
  if(sh == "skdetfsi_017_rhc") return kRHCNuisSKDetFSI017;
#endif
#ifdef SKCORRELATIONS
  if(sh == "skdetfsi_018_rhc") return kRHCNuisSKDetFSI018;
  if(sh == "skdetfsi_019_rhc") return kRHCNuisSKDetFSI019;
  if(sh == "skdetfsi_020_rhc") return kRHCNuisSKDetFSI020;
  if(sh == "skdetfsi_021_rhc") return kRHCNuisSKDetFSI021;
  if(sh == "skdetfsi_022_rhc") return kRHCNuisSKDetFSI022;
  if(sh == "skdetfsi_023_rhc") return kRHCNuisSKDetFSI023;
  if(sh == "skdetfsi_024_rhc") return kRHCNuisSKDetFSI024;
  if(sh == "skdetfsi_025_rhc") return kRHCNuisSKDetFSI025;
  if(sh == "skdetfsi_026_rhc") return kRHCNuisSKDetFSI026;
  if(sh == "skdetfsi_027_rhc") return kRHCNuisSKDetFSI027;
  if(sh == "skdetfsi_028_rhc") return kRHCNuisSKDetFSI028;
  if(sh == "skdetfsi_029_rhc") return kRHCNuisSKDetFSI029;
  if(sh == "skdetfsi_030_rhc") return kRHCNuisSKDetFSI030;
  if(sh == "skdetfsi_031_rhc") return kRHCNuisSKDetFSI031;
#endif

  cout << sh << " not yet supported in util::PhysParamIdFromString()" << endl;
  return kPhypUndefined;
}
//____________________________________________________________________________
t2kana::conv::PhysParamId_t t2kana::util::PhysParamIdFromLongString(string sh)
{

  //************************ FORWARD HORN CURRENT ****************************

  if(sh == "sin^2(2*theta_23)"   ) return kPhypMixing23;
  if(sh == "sin^2(theta_23)"     ) return kPhypMixing23;
  if(sh == "sin^2(2*theta_13)"   ) return kPhypMixing13;
  if(sh == "sin^2(theta_13)"     ) return kPhypMixing13;
  if(sh == "sin^2(2*theta_12)"   ) return kPhypMixing12;
  if(sh == "sin^2(theta_12)"     ) return kPhypMixing12;
  if(sh == "delta_cp"            ) return kPhypDeltaCP;
  if(sh == "|Dm^2_32|"           ) return kPhypDm23sq;
  if(sh == "|Dm^2_31|"           ) return kPhypDm23sq;
  if(sh == "|Dm^2|"              ) return kPhypDm23sq;
  if(sh == "Dm^2_12"             ) return kPhypDm12sq;
  if(sh == "alpha_decay"         ) return kPhypAlphaDecay;
  if(sh == "sin^2(2*theta_decay)") return kPhypSinsq2ThetaDecay;
  if(sh == "mu^2_decoh"          ) return kPhypMusqDecoh;
  if(sh == "sin^2(2*theta_decoh)") return kPhypSinsq2ThetaDecoh;

  if(sh == "SK energy scale"                ) return kNuisSuperkEnergyScale;
#ifndef SKFSITOGETHER
#if FSI_BINS != 3
  if(sh == "SK efficiency; 1Rmu; CCQE/Other") return kNuisSuperkEff1RmuCCQEoth;
  if(sh == "SK efficiency; 1Rmu; CCQE/RC0"  ) return kNuisSuperkEff1RmuCCQErc0u;
  if(sh == "SK efficiency; 1Rmu; CCQE/RC1"  ) return kNuisSuperkEff1RmuCCQErc1u;
  if(sh == "SK efficiency; 1Rmu; CCQE/RC2"  ) return kNuisSuperkEff1RmuCCQErc2u;
  if(sh == "SK efficiency; 1Rmu; CCnonQE"   ) return kNuisSuperkEff1RmuCCnonQE;
  if(sh == "SK efficiency; 1Rmu; NC"        ) return kNuisSuperkEff1RmuNC;
  if(sh == "SK efficiency; 1Rmu; nue CC"    ) return kNuisSuperkEff1RmuNueCC;
#endif
#if FSI_BINS != 2
  if(sh == "SK efficiency; 1Re; osc. nue CC E0") return kNuisSuperkEff1ReNueSigCC0;
  if(sh == "SK efficiency; 1Re; osc. nue CC E1") return kNuisSuperkEff1ReNueSigCC1;
  if(sh == "SK efficiency; 1Re; osc. nue CC E2") return kNuisSuperkEff1ReNueSigCC2;
  if(sh == "SK efficiency; 1Re; numu CC E0") return kNuisSuperkEff1ReNumuCC0;
  if(sh == "SK efficiency; 1Re; numu CC E1") return kNuisSuperkEff1ReNumuCC1;
  if(sh == "SK efficiency; 1Re; numu CC E2") return kNuisSuperkEff1ReNumuCC2;
  if(sh == "SK efficiency; 1Re; nue CC E0") return kNuisSuperkEff1ReNueCC0;
  if(sh == "SK efficiency; 1Re; nue CC E1") return kNuisSuperkEff1ReNueCC1;
  if(sh == "SK efficiency; 1Re; nue CC E2") return kNuisSuperkEff1ReNueCC2;
  if(sh == "SK efficiency; 1Re; NC E0") return kNuisSuperkEff1ReNC0;
  if(sh == "SK efficiency; 1Re; NC E1") return kNuisSuperkEff1ReNC1;
  if(sh == "SK efficiency; 1Re; NC E2") return kNuisSuperkEff1ReNC2;
#endif
#endif
  for(int ip = 0; ip < kNumOfFluxParamsInBANFFfit; ip++) {
    string comp = Form("BANFF %2d; SK %s flux; %.1f - %.1f GeV",
		ip, util::Pdg2String(kNuPdgCodeOfFluxParamInBANFFfit[ip]).c_str(), 
		kEtrueRangeOfFluxParamInBANFFfit[ip][0],
		kEtrueRangeOfFluxParamInBANFFfit[ip][1]);
    if(sh == comp) return (PhysParamId_t)((int)kNuisBANFF00+ip);
  }

#ifndef ONEOVERMAQE
  if(sh == "BANFF; Ma QE"                       ) return kNuisBANFFMAQE;
#else
  if(sh == "BANFF; 1 / (Ma QE)^2"               ) return kNuisBANFFMAQE;
#endif
  if(sh == "BANFF; Ma RES"                      ) return kNuisBANFFMARES;
  if(sh == "BANFF; Norm; CC QE; 0.0 - 1.5 GeV"  ) return kNuisBANFFCCQEE0;
#if BANFF_BINS != 3
  if(sh == "BANFF; Norm; CC QE; 1.5 - 3.5 GeV"  ) return kNuisBANFFCCQEE1;
  if(sh == "BANFF; Norm; CC QE; 3.5 - 30.0 GeV" ) return kNuisBANFFCCQEE2;
#endif
  if(sh == "BANFF; Norm; CC 1pi; 0.0 - 2.5 GeV" ) return kNuisBANFFCC1piE0;
#if BANFF_BINS != 3
  if(sh == "BANFF; Norm; CC 1pi; 2.5 - 30.0 GeV") return kNuisBANFFCC1piE1;
#endif
#if (BANFF_BINS != 2) && (BANFF_BINS != 4)
  if(sh == "BANFF; Norm; NC 1pi0; 0.0 - 30.0 GeV") return kNuisBANFFNC1pi0;
#endif
#if BANFF_BINS == 4
  if(sh == "BANFF; Norm; MEC")                     return kNuisBANFFMECNorm;
  if(sh == "BANFF; PDD double counting")           return kNuisBANFFPDDDoubleCount;
#endif


    //
    // XSec params not correlated with flux
    //

 if(sh == "Fermi Momentum"     ) return kNuisXSecFermiMom;
 if(sh == "Spectral Func"      ) return kNuisXSecSpecFunc;
 if(sh == "W-shape"            ) return kNuisXSecWshape;
 if(sh == "Shape; CC Oth"      ) return kNuisXSecShapeCCoth;
 if(sh == "Pi-less delta decay") return kNuisXSecPiLessDelta;
#ifndef __NUEONLY__
 if(sh == "Binding energy"     ) return kNuisXSecBindE;
#endif
 if(sh == "Norm; CC Coh"      ) return kNuisXSecNormCCcoh;
#ifndef __NUEONLY__
 if(sh == "Norm; NC 1piC"     ) return kNuisXSecNormNC1piC;
#endif
#ifndef __NUMUONLY__
 if(sh == "Norm; NC Coh"      ) return kNuisXSecNormNCcoh;
#endif
 if(sh == "Norm; NC Oth"      ) return kNuisXSecNormNCoth;
 if(sh == "Norm; NuE To NuMu" ) return kNuisXSecNormNuEToNuMu;
#ifndef __NUEONLY__
 if(sh == "Norm; AntiNu To Nu") return kNuisXSecNormAntiNuToNu;
#endif
#ifndef SKFSITOGETHER
 for(int ip = 0; ip < kNumOfFSIParams; ip++) {
   string nutypes = "";
   if(ip < 3)
     nutypes = "numu/numubar CCQE (1Rmu)";
   else if(ip < 4)
     nutypes = "numu/numubar CCnQE (1Rmu)";
   else if(ip < 5)
     nutypes = "nue/nuebar/signue CC (1Rmu)";
   else if(ip < 6)
     nutypes = "all NC (1Rmu)";
   else if(ip < 9)
     nutypes = "Osc. nue CC (1Re)";
   else if(ip < 12)
     nutypes = "numu/numubar CC (1Re)";
   else if(ip < 15)
     nutypes = "nue/nuebar CC (1Re)";
   else if(ip < 18)
     nutypes = "all NC (1Re)";
   string comp = Form("FSI/SI %2d; Ereco range %.1f - %.1f GeV; %s",ip,
		      kFSIParamsErecoCut[ip][0],kFSIParamsErecoCut[ip][1],nutypes.c_str());     
#if FSI_BINS != 3
   if(sh == comp) return (PhysParamId_t)((int)kNuisFSI000+ip);
#else
   if(sh == comp) return (PhysParamId_t)((int)kNuisFSI006+ip);
#endif
 }
#endif
#ifdef SKFSITOGETHER
#ifndef SKCORRELATIONS
 for(int ip = 0; ip < kNumOfSKDetFSIParams; ip++) {
   string nutypes = "";
   if(ip < 3)
     nutypes = "numu/numubar CCQE (1Rmu)";
   else if(ip < 4)
     nutypes = "numu/numubar CCnQE (1Rmu)";
   else if(ip < 5)
     nutypes = "nue/nuebar/signue CC (1Rmu)";
   else if(ip < 6)
     nutypes = "all NC (1Rmu)";
   else if(ip < 9)
     nutypes = "Osc. nue CC (1Re)";
   else if(ip < 12)
     nutypes = "numu/numubar CC (1Re)";
   else if(ip < 15)
     nutypes = "nue/nuebar CC (1Re)";
   else if(ip < 18)
     nutypes = "all NC (1Re)";
   string comp = Form("SKDet + FSI/SI %2d; Ereco range %.2f - %.2f GeV; %s",ip,
		      kSKDetFSIParamsErecoCut[ip][0],kSKDetFSIParamsErecoCut[ip][1],nutypes.c_str());     
#if FSI_BINS != 3
   if(sh == comp) return (PhysParamId_t)((int)kNuisSKDetFSI000+ip);
#else
   if(sh == comp) return (PhysParamId_t)((int)kNuisSKDetFSI006+ip);
#endif
 }
#else
 for(int ip = 0; ip < kNumOfSKDetFSIParams; ip++) {
   string nutypes = "";
   if(ip < 5)
     nutypes = "Osc. nue CC (1Rmu)";
   else if(ip < 10)
     nutypes = "numu/numubar CC (1Rmu)";
   else if(ip < 15)
     nutypes = "nue/nuebar CC (1Rmu)";
   else if(ip < 20)
     nutypes = "all NC (1Rmu)";
   else if(ip < 23)
     nutypes = "Osc. nue CC (1Re)";
   else if(ip < 26)
     nutypes = "numu/numubar CC (1Re)";
   else if(ip < 29)
     nutypes = "nue/nuebar CC (1Re)";
   else if(ip < 32)
     nutypes = "all NC (1Re)";
   string comp = Form("SKDet + FSI/SI %2d; Ereco range %.2f - %.2f GeV; %s",ip,
		      kSKDetFSIParamsErecoCut[ip][0],kSKDetFSIParamsErecoCut[ip][1],nutypes.c_str());     
   if(sh == comp) return (PhysParamId_t)((int)kNuisSKDetFSI000+ip);
 }
#endif
#endif

 if(sh == "RHC to FHC normalisation") return kNuisRHCToFHC;

 //************************* REVERSE HORN CURRENT ***************************

 if(sh == "SK energy scale; RHC"                ) return kRHCNuisSuperkEnergyScale;

 for(int ip = 0; ip < kNumOfRHCFluxParamsInBANFFfit; ip++) {
   string comp = Form("BANFF %2d; SK %s flux; %.1f - %.1f GeV; RHC",
		      ip, util::Pdg2String(kNuPdgCodeOfRHCFluxParamInBANFFfit[ip]).c_str(),
		      kEtrueRangeOfRHCFluxParamInBANFFfit[ip][0],
		      kEtrueRangeOfRHCFluxParamInBANFFfit[ip][1]);
   if(sh == comp) return (PhysParamId_t)((int)kRHCNuisBANFF00+ip);
 }

#ifndef ONEOVERMAQE
 if(sh == "BANFF; Ma QE; RHC"                       ) return kRHCNuisBANFFMAQE;
#else
 if(sh == "BANFF; 1 / (Ma QE)^2; RHC"               ) return kRHCNuisBANFFMAQE;
#endif
 if(sh == "BANFF; Ma RES; RHC"                      ) return kRHCNuisBANFFMARES;
 if(sh == "BANFF; Norm; CC QE; 0.0 - 1.5 GeV; RHC"  ) return kRHCNuisBANFFCCQEE0;
#if BANFF_BINS != 3
 if(sh == "BANFF; Norm; CC QE; 1.5 - 3.5 GeV; RHC"  ) return kRHCNuisBANFFCCQEE1;
 if(sh == "BANFF; Norm; CC QE; 3.5 - 30.0 GeV; RHC" ) return kRHCNuisBANFFCCQEE2;
#endif
 if(sh == "BANFF; Norm; CC 1pi; 0.0 - 2.5 GeV; RHC" ) return kRHCNuisBANFFCC1piE0;
#if BANFF_BINS != 3
 if(sh == "BANFF; Norm; CC 1pi; 2.5 - 30.0 GeV; RHC") return kRHCNuisBANFFCC1piE1;
#endif
#if (BANFF_BINS != 2) && (BANFF_BINS != 4)
 if(sh == "BANFF; Norm; NC 1pi0; 0.0 - 30.0 GeV; RHC") return kRHCNuisBANFFNC1pi0;
#endif
#if BANFF_BINS == 4
 if(sh == "BANFF; Norm; MEC; RHC")                     return kRHCNuisBANFFMECNorm;
 if(sh == "BANFF; PDD double counting; RHC")           return kRHCNuisBANFFPDDDoubleCount;
#endif

 //
 // XSec params not correlated with flux
 //

 if(sh == "Fermi Momentum; RHC"     ) return kRHCNuisXSecFermiMom;
 if(sh == "Spectral Func; RHC"      ) return kRHCNuisXSecSpecFunc;
 if(sh == "W-shape; RHC"            ) return kRHCNuisXSecWshape;
 if(sh == "Shape; CC Oth; RHC"      ) return kRHCNuisXSecShapeCCoth;
 if(sh == "Pi-less delta decay; RHC") return kRHCNuisXSecPiLessDelta;
#ifndef __NUEONLY__
 if(sh == "Binding energy; RHC"     ) return kRHCNuisXSecBindE;
#endif
 if(sh == "Norm; CC Coh; RHC"      ) return kRHCNuisXSecNormCCcoh;
#ifndef __NUEONLY__
 if(sh == "Norm; NC 1piC; RHC"     ) return kRHCNuisXSecNormNC1piC;
#endif
#ifndef __NUMUONLY__
 if(sh == "Norm; NC Coh; RHC"      ) return kRHCNuisXSecNormNCcoh;
#endif
 if(sh == "Norm; NC Oth; RHC"      ) return kRHCNuisXSecNormNCoth;
 if(sh == "Norm; NuE To NuMu; RHC" ) return kRHCNuisXSecNormNuEToNuMu;
#ifndef __NUEONLY__
 if(sh == "Norm; AntiNu To Nu; RHC") return kRHCNuisXSecNormAntiNuToNu;
#endif
#ifdef SKFSITOGETHER
#ifndef SKCORRELATIONS
 for(int ip = 0; ip < kNumOfRHCSKDetFSIParams; ip++) {
   string nutypes = "";
   if(ip < 3)
     nutypes = "numu/numubar CCQE (1Rmu)";
   else if(ip < 4)
     nutypes = "numu/numubar CCnQE (1Rmu)";
   else if(ip < 5)
     nutypes = "nue/nuebar/signue CC (1Rmu)";
   else if(ip < 6)
     nutypes = "all NC (1Rmu)";
   else if(ip < 9)
     nutypes = "Osc. nue CC (1Re)";
   else if(ip < 12)
     nutypes = "numu/numubar CC (1Re)";
   else if(ip < 15)
     nutypes = "nue/nuebar CC (1Re)";
   else if(ip < 18)
     nutypes = "all NC (1Re)";
   string comp = Form("SKDet + FSI/SI %2d; Ereco range %.2f - %.2f GeV; %s; RHC",ip,
                      kRHCSKDetFSIParamsErecoCut[ip][0],kRHCSKDetFSIParamsErecoCut[ip][1],nutypes.c_str());
#if FSI_BINS != 3
   if(sh == comp) return (PhysParamId_t)((int)kRHCNuisSKDetFSI000+ip);
#else
   if(sh == comp) return (PhysParamId_t)((int)kRHCNuisSKDetFSI006+ip);
#endif
 }
#else
 for(int ip = 0; ip < kNumOfRHCSKDetFSIParams; ip++) {
   string nutypes = "";
   if(ip < 5)
     nutypes = "Osc. nue CC (1Rmu)";
   else if(ip < 10)
     nutypes = "numu/numubar CC (1Rmu)";
   else if(ip < 15)
     nutypes = "nue/nuebar CC (1Rmu)";
   else if(ip < 20)
     nutypes = "all NC (1Rmu)";
   else if(ip < 23)
     nutypes = "Osc. nue CC (1Re)";
   else if(ip < 26)
     nutypes = "numu/numubar CC (1Re)";
   else if(ip < 29)
     nutypes = "nue/nuebar CC (1Re)";
   else if(ip < 32)
     nutypes = "all NC (1Re)";
   string comp = Form("SKDet + FSI/SI %2d; Ereco range %.2f - %.2f GeV; %s; RHC",ip,
                      kRHCSKDetFSIParamsErecoCut[ip][0],kRHCSKDetFSIParamsErecoCut[ip][1],nutypes.c_str());
   if(sh == comp) return (PhysParamId_t)((int)kRHCNuisSKDetFSI000+ip);
 }
#endif
#endif

  cout << sh << " not yet supported in util::PhysParamIdFromLongString()" << endl;
  return kPhypUndefined;
}
//____________________________________________________________________________
t2kana::conv::SKMCTemplateType_t   
  t2kana::util::SKMCTemplateTypeFromString(string sh) 
{
  if(sh == "ori")        return kSKMCtOri;
  if(sh == "nom")        return kSKMCtNom;
  if(sh == "tun")        return kSKMCtTun;
  if(sh == "sys")        return kSKMCtSys;
  if(sh == "fin")        return kSKMCtFin;

  return kSKMCtUndefined;
}
//____________________________________________________________________________
t2kana::conv::OscChannelId_t
  t2kana::util::GetOscChannel(int nu_start, int nu_end)
{
  if ( nu_start == kPdgNuE      && nu_end == kPdgNuE      ) return kOscNuE2NuE;
  if ( nu_start == kPdgNuE      && nu_end == kPdgNuMu     ) return kOscNuE2NuMu;
  if ( nu_start == kPdgNuE      && nu_end == kPdgNuTau    ) return kOscNuE2NuTau;
  if ( nu_start == kPdgNuMu     && nu_end == kPdgNuE      ) return kOscNuMu2NuE;
  if ( nu_start == kPdgNuMu     && nu_end == kPdgNuMu     ) return kOscNuMu2NuMu;
  if ( nu_start == kPdgNuMu     && nu_end == kPdgNuTau    ) return kOscNuMu2NuTau;
  if ( nu_start == kPdgNuTau    && nu_end == kPdgNuE      ) return kOscNuTau2NuE;
  if ( nu_start == kPdgNuTau    && nu_end == kPdgNuMu     ) return kOscNuTau2NuMu;
  if ( nu_start == kPdgNuTau    && nu_end == kPdgNuTau    ) return kOscNuTau2NuTau;
     
  if ( nu_start == kPdgNuEBar   && nu_end == kPdgNuEBar   ) return kOscNuEBar2NuEBar;
  if ( nu_start == kPdgNuEBar   && nu_end == kPdgNuMuBar  ) return kOscNuEBar2NuMuBar;
  if ( nu_start == kPdgNuEBar   && nu_end == kPdgNuTauBar ) return kOscNuEBar2NuTauBar;
  if ( nu_start == kPdgNuMuBar  && nu_end == kPdgNuEBar   ) return kOscNuMuBar2NuEBar;
  if ( nu_start == kPdgNuMuBar  && nu_end == kPdgNuMuBar  ) return kOscNuMuBar2NuMuBar;
  if ( nu_start == kPdgNuMuBar  && nu_end == kPdgNuTauBar ) return kOscNuMuBar2NuTauBar;
  if ( nu_start == kPdgNuTauBar && nu_end == kPdgNuEBar   ) return kOscNuTauBar2NuEBar;
  if ( nu_start == kPdgNuTauBar && nu_end == kPdgNuMuBar  ) return kOscNuTauBar2NuMuBar;
  if ( nu_start == kPdgNuTauBar && nu_end == kPdgNuTauBar ) return kOscNuTauBar2NuTauBar;
   
  return kOscUndefined;
}
//____________________________________________________________________________
int t2kana::util::SuperkMode2FluxNuPdg(t2kana::conv::SuperkMode_t m)
{
  if ( util::IsNuMu      (m) ) return kPdgNuMu;
  if ( util::IsOscNuMu   (m) ) return kPdgNuE; /// nue -> numu
  if ( util::IsNuMuBar   (m) ) return kPdgNuMuBar;
  if ( util::IsOscNuMuBar(m) ) return kPdgNuEBar; /// nue -> numu
  if ( util::IsNuE       (m) ) return kPdgNuE;
  if ( util::IsNuEBar    (m) ) return kPdgNuEBar;
  if ( util::IsOscNuEBar (m) ) return kPdgNuMuBar; /// nue -> numu
  if ( util::IsOscNuE    (m) ) return kPdgNuMu; /// numu -> nue
   
  // NOTE: Beam NuEBars not handled

  return 0;
}
//____________________________________________________________________________
int t2kana::util::SuperkMode2SKMCPdg(t2kana::conv::SuperkMode_t m)
{
  if ( util::IsNuMu      (m) ) return kPdgNuMu;
  if ( util::IsOscNuMu   (m) ) return kPdgNuMu;
  if ( util::IsNuMuBar   (m) ) return kPdgNuMuBar;
  if ( util::IsOscNuMuBar(m) ) return kPdgNuMuBar;
  if ( util::IsNuE       (m) ) return kPdgNuE;
  if ( util::IsNuEBar    (m) ) return kPdgNuEBar;
  if ( util::IsOscNuEBar (m) ) return kPdgNuEBar;
  if ( util::IsOscNuE    (m) ) return kPdgNuE;
   
  // NOTE: Beam NuEBars not handled

  return 0;
}
//____________________________________________________________________________
bool t2kana::util::IsNuMu(t2kana::conv::SuperkMode_t  m)
{

  if(m >= kSKmNuMuCCQE && m <= kSKmNuMuNCoth_e)
    return true;
  /*
  switch(m) {
     case ( kSKmNuMuCCQE        ) :
     case ( kSKmMEC             ) : 
     case ( kSKmNuMuCC1pi       ) :
     case ( kSKmNuMuCCcoh       ) :
     case ( kSKmNuMuCCoth    ) :
       //     case ( kSKmNuMuCCoth     ) :
     case ( kSKmNuMuNC1piC      ) :
     case ( kSKmNuMuNC1pi0      ) :
     case ( kSKmNuMuNCcoh       ) :
     case ( kSKmNuMuNCoth_e     ) :
     case ( kSKmNuMuNCoth_mu    ) :
       return true;
       break;
     default :                
       return false;
       break;
  }
  */
  return false;
}
//____________________________________________________________________________
bool t2kana::util::IsOscNuMu(t2kana::conv::SuperkMode_t  m)
{
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
  if(m >= kSKmOscNuMuCCQE && m <= kSKmOscNuMuCCoth)
    return true;
#endif
  /*
  switch(m) {
     case ( kSKmOscNuMuCCQE     ) : 
     case ( kSKmOscNuMuCC1pi    ) :
     case ( kSKmOscNuMuCCcoh    ) :
       //     case ( kSKmOscNuMuCCoth ) :
     case ( kSKmOscNuMuCCoth  ) :
       return true;
       break;
     default :                
       return false;
       break;
  }
  */
  return false;
}
//____________________________________________________________________________
bool t2kana::util::IsOscNuMuBar(t2kana::conv::SuperkMode_t  m)
{
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
  if(m >= kSKmOscNuMuBarCCQE && m <= kSKmOscNuMuBarCCoth)
    return true;
#endif
  /*
  switch(m) {
     case ( kSKmOscNuMuBarCCQE     ) : 
     case ( kSKmOscNuMuBarCC1pi    ) :
     case ( kSKmOscNuMuBarCCcoh    ) :
       //     case ( kSKmOscNuMuCCoth ) :
     case ( kSKmOscNuMuBarCCoth  ) :
       return true;
       break;
     default :                
       return false;
       break;
       }*/
  return false;
}
//____________________________________________________________________________
bool t2kana::util::IsOscNuEBar(t2kana::conv::SuperkMode_t  m)
{
#ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
  if(m >= kSKmOscNuEBarCCQE && m <= kSKmOscNuEBarCCoth)
    return true;
#endif
  /*switch(m) {
     case ( kSKmOscNuEBarCCQE     ) : 
     case ( kSKmOscNuEBarCC1pi    ) :
     case ( kSKmOscNuEBarCCcoh    ) :
       //     case ( kSKmOscNuMuCCoth ) :
     case ( kSKmOscNuEBarCCoth  ) :
       return true;
       break;
     default :                
       return false;
       break;
       }*/
  return false;
}
//____________________________________________________________________________
bool t2kana::util::IsNuMuBar(t2kana::conv::SuperkMode_t m)
{
  if(m >= kSKmNuMuBarCCQE && m <= kSKmNuMuBarNCoth_e)
    return true;
  /*
  switch(m) {
     case (kSKmNuMuBarCCQE    ):
     case (kSKmNuMuBarCC1pi   ):
     case (kSKmNuMuBarCCcoh   ):
       //     case (kSKmNuMuBarCCoth):
     case (kSKmNuMuBarCCoth ):
     case (kSKmNuMuBarNC1piC  ):
     case (kSKmNuMuBarNC1pi0  ):
     case (kSKmNuMuBarNCcoh   ):
     case (kSKmNuMuBarNCoth_e ):
     case (kSKmNuMuBarNCoth_mu):
       return true;
       break;
     default :                
       return false;
       break;
  }
  */
  return false;
}
//____________________________________________________________________________
bool t2kana::util::IsNuE(t2kana::conv::SuperkMode_t m)
{
  if(m >= kSKmNuECCQE && m <= kSKmNuENCoth_e)
    return true;
  /*
  switch(m) {
     case (kSKmNuECCQE    ):
     case (kSKmNuECC1pi   ):
     case (kSKmNuECCcoh   ):
       //     case (kSKmNuECCoth):
     case (kSKmNuECCoth ):
     case (kSKmNuENC1piC  ):
     case (kSKmNuENC1pi0  ):
     case (kSKmNuENCcoh   ):
     case (kSKmNuENCoth_e ):
     case (kSKmNuENCoth_mu):
       return true;
       break;
     default :                
       return false;
       break;
  }
  */
  return false;
}
//____________________________________________________________________________
bool t2kana::util::IsNuEBar(t2kana::conv::SuperkMode_t m)
{
  if(m >= kSKmNuEBarCCQE && m <= kSKmNuEBarNCoth_e)
    return true;
  /*
  switch(m) {
     case (kSKmNuEBarCCQE    ):
     case (kSKmNuEBarCC1pi   ):
     case (kSKmNuEBarCCcoh   ):
       //     case (kSKmNuEBarCCoth):
     case (kSKmNuEBarCCoth ):
     case (kSKmNuEBarNC1piC  ):
     case (kSKmNuEBarNC1pi0  ):
     case (kSKmNuEBarNCcoh   ):
     case (kSKmNuEBarNCoth_e ):
     case (kSKmNuEBarNCoth_mu):
       return true;
       break;
     default :                
       return false;
       break;
  }
  */
  return false;
}
//____________________________________________________________________________
bool t2kana::util::IsOscNuE(t2kana::conv::SuperkMode_t m)
{
  if(m >= kSKmOscNuECCQE && m <= kSKmOscNuECCoth)
    return true;
  /*
  switch(m) {
     case ( kSKmSigNuECCQE     ) :
     case ( kSKmSigNuECC1pi    ) :
     case ( kSKmSigNuECCcoh    ) :
       //     case ( kSKmSigNuECCoth  ) :
     case ( kSKmSigNuECCoth ) :
       return true;
       break;
     default :                
       return false;
       break;
  }
  */
  return false;
}
//____________________________________________________________________________
bool t2kana::util::IsCC(t2kana::conv::SuperkMode_t m)
{
  return !(IsNC(m));
}
//____________________________________________________________________________
bool t2kana::util::IsNC(t2kana::conv::SuperkMode_t m)
{
  switch(m) {
     case ( kSKmNuMuNC1piC      ) : 
     case ( kSKmNuMuNC1pi0      ) : 
     case ( kSKmNuMuNCcoh       ) : 
     case ( kSKmNuMuNCoth_e     ) : 
     case ( kSKmNuMuNCoth_mu    ) : 
     case ( kSKmNuMuBarNC1piC   ) : 
     case ( kSKmNuMuBarNC1pi0   ) : 
     case ( kSKmNuMuBarNCcoh    ) : 
     case ( kSKmNuMuBarNCoth_e  ) : 
     case ( kSKmNuMuBarNCoth_mu ) : 
     case ( kSKmNuENC1piC       ) : 
     case ( kSKmNuENC1pi0       ) : 
     case ( kSKmNuENCcoh        ) : 
     case ( kSKmNuENCoth_e      ) : 
     case ( kSKmNuENCoth_mu     ) : 
     case ( kSKmNuEBarNC1piC    ) : 
     case ( kSKmNuEBarNC1pi0    ) : 
     case ( kSKmNuEBarNCcoh     ) : 
     case ( kSKmNuEBarNCoth_e   ) : 
     case ( kSKmNuEBarNCoth_mu  ) : 
	return true;
        break;
     default: 
	return false;
	break;
  }
  return false;
}
//____________________________________________________________________________
bool t2kana::util::IsMEC(t2kana::conv::SuperkMode_t m)
{
  switch(m) {
     case ( kSKmMEC      ) :
	return true;
	break;
     default: 
	return false;
	break;
  }
  return false;
}
//____________________________________________________________________________
bool t2kana::util::IsBANFF(PhysParamId_t p)
{
  if(p >= kNuisBANFF00
#if ((BANFF_BINS != 2)&&(BANFF_BINS != 4))
     && p <= kNuisBANFFNC1pi0
#elif BANFF_BINS == 4
     && p <= kNuisBANFFPDDDoubleCount
#else
     && p <= kNuisBANFFCC1piE1
#endif
     )
    return true;
  /*
  switch(p) {
    case ( kNuisBANFF00               ) : return true;  break;
    case ( kNuisBANFF01               ) : return true;  break;
    case ( kNuisBANFF02               ) : return true;  break;
    case ( kNuisBANFF03               ) : return true;  break;
    case ( kNuisBANFF04               ) : return true;  break;
    case ( kNuisBANFF05               ) : return true;  break;
    case ( kNuisBANFF06               ) : return true;  break;
    case ( kNuisBANFF07               ) : return true;  break;
    case ( kNuisBANFF08               ) : return true;  break;
    case ( kNuisBANFF09               ) : return true;  break;
    case ( kNuisBANFF10               ) : return true;  break;
    case ( kNuisBANFF11               ) : return true;  break;
    case ( kNuisBANFF12               ) : return true;  break;
    case ( kNuisBANFF13               ) : return true;  break;
    case ( kNuisBANFF14               ) : return true;  break;
    case ( kNuisBANFF15               ) : return true;  break;
#if (BANFF_BINS != 2) && (BANFF_BINS != 4)
    case ( kNuisBANFF16               ) : return true;  break;
    case ( kNuisBANFF17               ) : return true;  break;
    case ( kNuisBANFF18               ) : return true;  break;
    case ( kNuisBANFF19               ) : return true;  break;
    case ( kNuisBANFF20               ) : return true;  break;
    case ( kNuisBANFF21               ) : return true;  break;
#endif
#if BANFF_BINS == 1
    case ( kNuisBANFF22               ) : return true;  break;
    case ( kNuisBANFF23               ) : return true;  break;
    case ( kNuisBANFF24               ) : return true;  break;
#endif
    case ( kNuisBANFFMAQE             ) : return true;  break;
    case ( kNuisBANFFMARES            ) : return true;  break;
    case ( kNuisBANFFCCQEE0           ) : return true;  break;
#if BANFF_BINS != 3
    case ( kNuisBANFFCCQEE1           ) : return true;  break;
    case ( kNuisBANFFCCQEE2           ) : return true;  break;
#endif
    case ( kNuisBANFFCC1piE0          ) : return true;  break;
#if BANFF_BINS != 3
    case ( kNuisBANFFCC1piE1          ) : return true;  break;
#endif
#if (BANFF_BINS != 2) && (BANFF_BINS != 4)
    case ( kNuisBANFFNC1pi0           ) : return true;  break;
#endif
#if BANFF_BINS == 4
    case ( kNuisBANFFMECNorm          ) : return true;  break;
    case ( kNuisBANFFPDDDoubleCount   ) : return true;  break;
#endif

  default : 
    return false;
  }
  */
  return false;
  
}
//____________________________________________________________________________                                                                                
bool t2kana::util::IsXSec(PhysParamId_t p)
{
  if(p >= kNuisXSecFermiMom && p <= kNuisXSecNormAntiNuToNu)
    return true;
  /*
  switch(p) {
    case ( kNuisXSecFermiMom          ) : return true;  break;
    case ( kNuisXSecSpecFunc          ) : return true;  break; 
    case ( kNuisXSecWshape            ) : return true;  break;
    case ( kNuisXSecShapeCCoth        ) : return true;  break;
    case ( kNuisXSecPiLessDelta       ) : return true;  break;
#ifndef __NUEONLY__
    case ( kNuisXSecBindE             ) : return true;  break;
#endif
    case ( kNuisXSecNormCCcoh         ) : return true;  break;
#ifndef __NUEONLY__
    case ( kNuisXSecNormNC1piC        ) : return true;  break;
#endif
#ifndef __NUMUONLY__
    case ( kNuisXSecNormNCcoh         ) : return true;  break;
#endif
    case ( kNuisXSecNormNCoth         ) : return true;  break;
    case ( kNuisXSecNormNuEToNuMu     ) : return true;  break;
#ifndef __NUEONLY__
    case ( kNuisXSecNormAntiNuToNu    ) : return true;  break;
#endif
  default : 
    return false;
  }
  */
  return false;
}
//____________________________________________________________________________                                                                                
bool t2kana::util::IsUsingSplines(PhysParamId_t p)
{
  if(p >= kNuisXSecFermiMom && p <= kNuisXSecBindE)
    return true;
  if(p == kNuisBANFFMAQE || p == kNuisBANFFMAQE)
    return true;
#if BANFF_BINS == 4
  if(p == kNuisBANFFMECNorm || p == kNuisBANFFPDDDoubleCount)
    return true;
#endif
  /*
  switch(p) {
  case ( kNuisXSecFermiMom          ) : return true;  break;
  case ( kNuisXSecSpecFunc          ) : return true;  break;
  case ( kNuisXSecWshape            ) : return true;  break;
  case ( kNuisXSecShapeCCoth        ) : return true;  break;
  case ( kNuisXSecPiLessDelta       ) : return true;  break;
  case ( kNuisXSecBindE             ) : return true;  break;
  case ( kNuisBANFFMAQE             ) : return true;  break;
  case ( kNuisBANFFMARES            ) : return true;  break;
#if BANFF_BINS == 4
  case ( kNuisBANFFMECNorm          ) : return true;  break;
  case ( kNuisBANFFPDDDoubleCount   ) : return true;  break;
#endif

  default :
    return false;
  }
  */
  return false;
}
//____________________________________________________________________________                                                                                
bool t2kana::util::IsEff1Rmu(PhysParamId_t p)
{
  switch(p) {
#ifndef SKFSITOGETHER
#if FSI_BINS != 3
  case ( kNuisSuperkEnergyScale     ) : return true;  break;
  case ( kNuisSuperkEff1RmuCCQEoth  ) : return true;  break;
  case ( kNuisSuperkEff1RmuCCQErc0u ) : return true;  break;
  case ( kNuisSuperkEff1RmuCCQErc1u ) : return true;  break;
  case ( kNuisSuperkEff1RmuCCQErc2u ) : return true;  break;
  case ( kNuisSuperkEff1RmuCCnonQE  ) : return true;  break;
  case ( kNuisSuperkEff1RmuNC       ) : return true;  break;
  case ( kNuisSuperkEff1RmuNueCC    ) : return true;  break;
#endif
#endif
  default :
    return false;
  }
  return false; 
}
//____________________________________________________________________________                                                                                  
bool t2kana::util::IsEff1Re(PhysParamId_t p)
{
  switch(p) {
#if FSI_BINS != 2
#ifndef SKFSITOGETHER
  case ( kNuisSuperkEnergyScale     ) : return true;  break;
  case ( kNuisSuperkEff1ReNueSigCC0 ) : return true;  break;
  case ( kNuisSuperkEff1ReNueSigCC1 ) : return true;  break;
  case ( kNuisSuperkEff1ReNueSigCC2 ) : return true;  break;
  case ( kNuisSuperkEff1ReNumuCC0   ) : return true;  break;
  case ( kNuisSuperkEff1ReNumuCC1   ) : return true;  break;
  case ( kNuisSuperkEff1ReNumuCC2   ) : return true;  break;
  case ( kNuisSuperkEff1ReNueCC0    ) : return true;  break;
  case ( kNuisSuperkEff1ReNueCC1    ) : return true;  break;
  case ( kNuisSuperkEff1ReNueCC2    ) : return true;  break;
  case ( kNuisSuperkEff1ReNC0       ) : return true;  break;
  case ( kNuisSuperkEff1ReNC1       ) : return true;  break;
  case ( kNuisSuperkEff1ReNC2       ) : return true;  break;
#endif
#endif
  default :
    return false;
  }
  return false;
}
//____________________________________________________________________________                                                                                 
bool t2kana::util::IsEff(PhysParamId_t p)
{
  switch(p) {
  case ( kNuisSuperkEnergyScale     ) : return true;  break;
#ifndef SKFSITOGETHER
#if FSI_BINS != 3
  case ( kNuisSuperkEff1RmuCCQEoth  ) : return true;  break;
  case ( kNuisSuperkEff1RmuCCQErc0u ) : return true;  break;
  case ( kNuisSuperkEff1RmuCCQErc1u ) : return true;  break;
  case ( kNuisSuperkEff1RmuCCQErc2u ) : return true;  break;
  case ( kNuisSuperkEff1RmuCCnonQE  ) : return true;  break;
  case ( kNuisSuperkEff1RmuNC       ) : return true;  break;
  case ( kNuisSuperkEff1RmuNueCC    ) : return true;  break;
#endif
#if FSI_BINS != 2
  case ( kNuisSuperkEff1ReNueSigCC0 ) : return true;  break;
  case ( kNuisSuperkEff1ReNueSigCC1 ) : return true;  break;
  case ( kNuisSuperkEff1ReNueSigCC2 ) : return true;  break;
  case ( kNuisSuperkEff1ReNumuCC0   ) : return true;  break;
  case ( kNuisSuperkEff1ReNumuCC1   ) : return true;  break;
  case ( kNuisSuperkEff1ReNumuCC2   ) : return true;  break;
  case ( kNuisSuperkEff1ReNueCC0    ) : return true;  break;
  case ( kNuisSuperkEff1ReNueCC1    ) : return true;  break;
  case ( kNuisSuperkEff1ReNueCC2    ) : return true;  break;
  case ( kNuisSuperkEff1ReNC0       ) : return true;  break;
  case ( kNuisSuperkEff1ReNC1       ) : return true;  break;
  case ( kNuisSuperkEff1ReNC2       ) : return true;  break;
#endif
#endif
  default :
    return false;
  }
  return false;
}
//____________________________________________________________________________                                                                                   
bool t2kana::util::IsFSI(PhysParamId_t p)
{
  switch(p) {
#ifndef SKFSITOGETHER
#if FSI_BINS != 3
    case ( kNuisFSI000                ) : return true;  break;
    case ( kNuisFSI001                ) : return true;  break;
    case ( kNuisFSI002                ) : return true;  break;
    case ( kNuisFSI003                ) : return true;  break;
    case ( kNuisFSI004                ) : return true;  break;
    case ( kNuisFSI005                ) : return true;  break;
#endif
#if FSI_BINS != 2
    case ( kNuisFSI006             ) : return true;  break;
    case ( kNuisFSI007             ) : return true;  break;
    case ( kNuisFSI008             ) : return true;  break;
    case ( kNuisFSI009             ) : return true;  break;
    case ( kNuisFSI010             ) : return true;  break;
    case ( kNuisFSI011             ) : return true;  break;
    case ( kNuisFSI012             ) : return true;  break;
    case ( kNuisFSI013             ) : return true;  break;
    case ( kNuisFSI014             ) : return true;  break;
    case ( kNuisFSI015             ) : return true;  break;
    case ( kNuisFSI016             ) : return true;  break;
    case ( kNuisFSI017             ) : return true;  break;
#endif
#endif
  default :
    return false;
  }
  return false;
}
//____________________________________________________________________________                                                                                   
bool t2kana::util::IsSKDetFSI(PhysParamId_t p)
{
#ifdef SKCORRELATIONS
  if(p >= kNuisSKDetFSI000 && p <= kNuisSKDetFSI031)
    return true;
#else
#if FSI_BINS != 3
  if(p >= kNuisSKDetFSI000 && p <= kNuisSKDetFSI005)
    return true;
#endif
#if FSI_BINS != 2
  if(p >= kNuisSKDetFSI006 && p <= kNuisSKDetFSI017)
    return true;
#endif
#endif
  if(p == kNuisSuperkEnergyScale)
    return true;

  return false;
}
//____________________________________________________________________________
bool t2kana::util::IsSystematic(PhysParamId_t p)
{
  if(p <= kPhypSinsq2ThetaDecoh)
    return false;
  return true;
}
//____________________________________________________________________________
bool t2kana::util::IsFHCSystematic(PhysParamId_t p)
{
  return !IsRHCSystematic(p);
}
//____________________________________________________________________________
bool t2kana::util::IsRHCSystematic(PhysParamId_t p)
{
  if(!IsSystematic(p))
    return false;
  if(p < kNuisRHCToFHC)
    return false;
  return true;
}
//____________________________________________________________________________
bool t2kana::util::IsMixingAngle(PhysParamId_t p)
{
  switch(p) {
    case ( kPhypMixing23              ) : return true; break;
    case ( kPhypMixing13              ) : return true; break;
    case ( kPhypMixing12              ) : return true; break;
    case ( kPhypSinsq2ThetaDecay      ) : return true; break;
    case ( kPhypSinsq2ThetaDecoh      ) : return true; break;

    default : 
      return false;
  }
  return false;
}

//____________________________________________________________________________
bool t2kana::util::IsNeutrino(int pdg_code)
{
  if ( pdg_code == kPdgNuE   ||
       pdg_code == kPdgNuMu  ||
       pdg_code == kPdgNuTau ) 
  {
    return true;
  }
  return false;
}
//___________________________________________________________________________
bool t2kana::util::IsAntiNeutrino(int pdg_code)
{
  if ( pdg_code == kPdgNuEBar   ||
       pdg_code == kPdgNuMuBar  ||
       pdg_code == kPdgNuTauBar ) 
  {
    return true;
  }
  return false;
}
//___________________________________________________________________________
bool t2kana::util::IsFHCSample(conv::SuperkSample_t s)
{
  if(s == kSKsRHCFCFV1RingMuLike || s == kSKsRHCFCFV1RingELike)
    return false;
  return true;
}
//___________________________________________________________________________
bool t2kana::util::IsRHCSample(conv::SuperkSample_t s)
{
  if(s == kSKsRHCFCFV1RingMuLike || s == kSKsRHCFCFV1RingELike)
    return true;
  return false;
}
//___________________________________________________________________________
bool t2kana::util::Is1RMuSample(conv::SuperkSample_t s)
{
  if(s == kSKsRHCFCFV1RingMuLike || s == kSKsFCFV1RingMuLike)
    return true;
  return false;
}
//___________________________________________________________________________
bool t2kana::util::Is1RESample(conv::SuperkSample_t s)
{
  if(s == kSKsRHCFCFV1RingELike || s == kSKsFCFV1RingELike)
    return true;
  return false;
}
//___________________________________________________________________________
int t2kana::util::NumErecoBins(t2kana::conv::SuperkSample_t s)
{
  switch(s) {
     case ( kSKsFCFV1RingMuLike ) : return kNErecoBins1Rmu; break;
     case ( kSKsFCFV1RingELike  ) : return kNErecoBins1Re;  break;
     case ( kSKsRHCFCFV1RingMuLike ) : return kRHCNErecoBins1Rmu; break;
     case ( kSKsRHCFCFV1RingELike  ) : return kRHCNErecoBins1Re;  break;

     default                      : return 0;
  }
  return 0;
}
//___________________________________________________________________________
const double * t2kana::util::ErecoBinEdges(t2kana::conv::SuperkSample_t s)
{
  switch(s) {
     case ( kSKsFCFV1RingMuLike ) : return kErecoBinEdges1Rmu; break;
     case ( kSKsFCFV1RingELike  ) : return kErecoBinEdges1Re;  break;
     case ( kSKsRHCFCFV1RingMuLike ) : return kRHCErecoBinEdges1Rmu; break;
     case ( kSKsRHCFCFV1RingELike  ) : return kRHCErecoBinEdges1Re;  break;

     default                      : return 0;
  }
  return 0;
}
//___________________________________________________________________________
double t2kana::util::ErecoBinEdges(t2kana::conv::SuperkSample_t s, int i)
{
  switch(s) {
     case ( kSKsFCFV1RingMuLike ) : return kErecoBinEdges1Rmu[i]; break;
     case ( kSKsFCFV1RingELike  ) : return kErecoBinEdges1Re[i];  break;
     case ( kSKsRHCFCFV1RingMuLike ) : return kRHCErecoBinEdges1Rmu[i]; break;
     case ( kSKsRHCFCFV1RingELike  ) : return kRHCErecoBinEdges1Re[i];  break;

     default                      : return 0;
  }
  return 0;
}
//___________________________________________________________________________
int t2kana::util::NumXsecErecoBins(t2kana::conv::SuperkSample_t s)
{
  switch(s) {
     case ( kSKsFCFV1RingMuLike ) : return kXsNErecoBins1Rmu; break;
     case ( kSKsFCFV1RingELike  ) : return kXsNErecoBins1Re;  break;
     case ( kSKsRHCFCFV1RingMuLike ) : return kRHCXsNErecoBins1Rmu; break;
     case ( kSKsRHCFCFV1RingELike  ) : return kRHCXsNErecoBins1Re;  break;

     default                      : return 0;
  }
  return 0;
}
//___________________________________________________________________________
const double * t2kana::util::XsecErecoBinEdges(t2kana::conv::SuperkSample_t s)
{
  switch(s) {
     case ( kSKsFCFV1RingMuLike ) : return kXsErecoBinEdges1Rmu; break;
     case ( kSKsFCFV1RingELike  ) : return kXsErecoBinEdges1Re;  break;
     case ( kSKsRHCFCFV1RingMuLike ) : return kRHCXsErecoBinEdges1Rmu; break;
     case ( kSKsRHCFCFV1RingELike  ) : return kRHCXsErecoBinEdges1Re;  break;

     default                      : return 0;
  }
  return 0;
}
//___________________________________________________________________________
double t2kana::util::XsecErecoBinEdges(t2kana::conv::SuperkSample_t s, int i)
{
  switch(s) {
     case ( kSKsFCFV1RingMuLike ) : return kXsErecoBinEdges1Rmu[i]; break;
     case ( kSKsFCFV1RingELike  ) : return kXsErecoBinEdges1Re[i];  break;
     case ( kSKsRHCFCFV1RingMuLike ) : return kRHCXsErecoBinEdges1Rmu[i]; break;
     case ( kSKsRHCFCFV1RingELike  ) : return kRHCXsErecoBinEdges1Re[i];  break;

     default                      : return 0;
  }
  return 0;
}
//___________________________________________________________________________
int t2kana::util::NumEtrueBins(t2kana::conv::SuperkSample_t s)
{
  switch(s) {
     case ( kSKsFCFV1RingMuLike ) : return kNEtrueBins1Rmu; break;
     case ( kSKsFCFV1RingELike  ) : return kNEtrueBins1Re;  break;
     case ( kSKsRHCFCFV1RingMuLike ) : return kRHCNEtrueBins1Rmu; break;
     case ( kSKsRHCFCFV1RingELike  ) : return kRHCNEtrueBins1Re;  break;

     default                      : return 0;
  }
  return 0;
}
//___________________________________________________________________________
const double * t2kana::util::EtrueBinEdges(t2kana::conv::SuperkSample_t s)
{
  switch(s) {
     case ( kSKsFCFV1RingMuLike ) : return kEtrueBinEdges1Rmu; break;
     case ( kSKsFCFV1RingELike  ) : return kEtrueBinEdges1Re;  break;
     case ( kSKsRHCFCFV1RingMuLike ) : return kRHCEtrueBinEdges1Rmu; break;
     case ( kSKsRHCFCFV1RingELike  ) : return kRHCEtrueBinEdges1Re;  break;

     default                      : return 0;
  }
  return 0;
}
//___________________________________________________________________________
double  t2kana::util::EtrueBinEdges(t2kana::conv::SuperkSample_t s, int i)
{
  switch(s) {
     case ( kSKsFCFV1RingMuLike ) : return kEtrueBinEdges1Rmu[i]; break;
     case ( kSKsFCFV1RingELike  ) : return kEtrueBinEdges1Re[i];  break;
     case ( kSKsRHCFCFV1RingMuLike ) : return kRHCEtrueBinEdges1Rmu[i]; break;
     case ( kSKsRHCFCFV1RingELike  ) : return kRHCEtrueBinEdges1Re[i];  break;

     default                      : return 0;
  }
  return 0;
}
//___________________________________________________________________________
int t2kana::util::NumVariableXsEtrueBins(t2kana::conv::SuperkSample_t s,
                                         t2kana::conv::SuperkMode_t m)
{
  int is = util::ArrIdx::Instance()->SuperkSampleIdx(s);
  if(is == -1)
    return 0;
  int im = util::ArrIdx::Instance()->SuperkModeIdx(s,m);
  if(im == -1)
    return 0;

  //need to start counter at zero. Want nbins. Looping over bin edges.
  int nbins = -1;

  for(int it = 0; it <= t2kana::conv::kXsNEtrueBins; it++) {
    if(t2kana::conv::kXsEtrueBinEdges[is][im][it] < 0)
      break;
    nbins++;
    }

  return nbins;
}
//___________________________________________________________________________
const double * t2kana::util::VariableXsEtrueBinEdges(t2kana::conv::SuperkSample_t s,
                                                     t2kana::conv::SuperkMode_t m)
{
  int is = util::ArrIdx::Instance()->SuperkSampleIdx(s);
  if(is == -1)
    return 0;
  int im = util::ArrIdx::Instance()->SuperkModeIdx(s,m);
  if(im == -1)
    return 0;

  return t2kana::conv::kXsEtrueBinEdges[is][im];
}
//___________________________________________________________________________
int t2kana::util::NumVariableMCEtrueBins(t2kana::conv::SuperkSample_t s,
                                         t2kana::conv::SuperkMode_t m)
{
  int is = util::ArrIdx::Instance()->SuperkSampleIdx(s);
  if(is == -1)
    return 0;
  int im = util::ArrIdx::Instance()->SuperkModeIdx(s,m);
  if(im == -1)
    return 0;

  //need to start counter at zero. Want nbins. Looping over bin edges.
  int nbins = -1;

  for(int it = 0; it <= t2kana::conv::kMCNEtrueBins; it++) {
    if(t2kana::conv::kMCEtrueBinEdges[is][im][it] < 0)
      break;
    nbins++;
    }

  return nbins;
}
//___________________________________________________________________________
const double * t2kana::util::VariableMCEtrueBinEdges(t2kana::conv::SuperkSample_t s,
                                                     t2kana::conv::SuperkMode_t m)
{
  int is = util::ArrIdx::Instance()->SuperkSampleIdx(s);
  if(is == -1)
    return 0;
  int im = util::ArrIdx::Instance()->SuperkModeIdx(s,m);
  if(im == -1)
    return 0;

  return t2kana::conv::kMCEtrueBinEdges[is][im];
}
//___________________________________________________________________________
void t2kana::util::SetDefaultStyle(bool black_n_white)
{
  gROOT->SetStyle("Plain");

  gStyle -> SetPadTickX (1);
  gStyle -> SetPadTickY (1);

  //
  // Turn off all borders
  //
  gStyle -> SetCanvasBorderMode (0);
  gStyle -> SetFrameBorderMode  (0);
  gStyle -> SetPadBorderMode    (0);
  gStyle -> SetDrawBorder       (0);
  gStyle -> SetCanvasBorderSize (0);
  gStyle -> SetFrameBorderSize  (0);
  gStyle -> SetPadBorderSize    (0);
  gStyle -> SetTitleBorderSize  (0);

  //
  // Set the size of the default canvas
  //
  gStyle -> SetCanvasDefH (600);
  gStyle -> SetCanvasDefW (730);
  gStyle -> SetCanvasDefX  (10);
  gStyle -> SetCanvasDefY  (10);

  //
  // Set marker style
  //
  gStyle -> SetMarkerStyle (20);
  gStyle -> SetMarkerSize   (1);

  //
  // Set line widths
  //
  gStyle -> SetFrameLineWidth (1);
  gStyle -> SetFuncWidth      (2);
  gStyle -> SetHistLineWidth  (3);
  gStyle -> SetFuncColor      (2);
  gStyle -> SetFuncWidth      (3);

  //
  // Set margins 
  //
  gStyle -> SetPadTopMargin    (0.10);
  gStyle -> SetPadBottomMargin (0.20);
  gStyle -> SetPadLeftMargin   (0.15);
  gStyle -> SetPadRightMargin  (0.03);

  //
  // Set tick marks and turn off grids
  //
  gStyle -> SetNdivisions (505,"xyz");

  //
  // Adjust size and placement of axis labels
  //
  gStyle -> SetLabelSize   (0.050,  "xyz");
  gStyle -> SetLabelOffset (0.005,  "x"  );
  gStyle -> SetLabelOffset (0.005,  "y"  );
  gStyle -> SetLabelOffset (0.005,  "z"  );
  gStyle -> SetTitleSize   (0.060,  "xyz");
  gStyle -> SetTitleOffset (1.200,  "xz" );
  gStyle -> SetTitleOffset (1.000,  "y"  );

  //
  // Adjust fonts
  //
  gStyle -> SetLabelFont(42, "xyz");
  gStyle -> SetTitleFont(42, "xyz");
  gStyle -> SetStatFont (42);

  // Set Data/Stat/... and other options
  //
  gStyle -> SetOptDate          (0);
  gStyle -> SetOptFile          (0);
  gStyle -> SetOptStat          (0);
  gStyle -> SetStatFormat       ("6.2f");
  gStyle -> SetFitFormat        ("8.4f");
  gStyle -> SetOptFit           (1);
  gStyle -> SetStatH            (0.20);
  gStyle -> SetStatStyle        (0);
  gStyle -> SetStatW            (0.30);
  gStyle -> SetStatX            (0.845);
  gStyle -> SetStatY            (0.845);
  gStyle -> SetOptTitle         (0);
  gStyle -> SetTitleX           (0.15);
  gStyle -> SetTitleW           (0.75);
  gStyle -> SetTitleY           (0.90);
  gStyle -> SetPalette          (1);
  gStyle -> SetLegendBorderSize (0);

  //
  // Set paper size for life in the US or EU
  //
  gStyle -> SetPaperSize (TStyle::kA4);       //<-- tartes aux fraises
//gStyle -> SetPaperSize (TStyle::kUSLetter); //<-- donuts

  //
  // In B&W (papers)
  //
  if(black_n_white){
    const int ncol = 7;

    double red   [ncol];
    double green [ncol];
    double blue  [ncol];
    double stops [ncol];

    double dcol = -1/double(ncol);
    double gray = 1;
    for (int j = 0; j < ncol; j++) {
      // Define color with RGB equal to : gray, gray, gray
      stops[j] = double(j)/double(ncol-1);
      red  [j] = gray;
      blue [j] = gray;
      green[j] = gray;
      
      gray += dcol;
    }
    UInt_t totcol=50;
    TColor::CreateGradientColorTable(ncol,stops,red,green,blue,totcol); 

    gStyle -> SetFuncWidth     (1);
    gStyle -> SetHistLineWidth (1);
    gStyle -> SetFuncColor     (1);
    gStyle -> SetFuncWidth     (1);
  }//bw

  gROOT->ForceStyle();
}
//___________________________________________________________________________
void t2kana::util::Format(
   TGraph* gr, int lcol, int lsty, int lwid, int mcol, int msty, double msiz)
{
  if(!gr) return;

  if (lcol >= 0) gr -> SetLineColor   (lcol);
  if (lsty >= 0) gr -> SetLineStyle   (lsty);
  if (lwid >= 0) gr -> SetLineWidth   (lwid);


  if (mcol >= 0) gr -> SetMarkerColor (mcol);
  if (msty >= 0) gr -> SetMarkerStyle (msty);
  if (msiz >= 0) gr -> SetMarkerSize  (msiz);
}
//___________________________________________________________________________
void t2kana::util::Format(
  TH1* hst, int lcol, int lsty, int lwid, int mcol, int msty, double msiz)
{
  if(!hst) return;

  if (lcol >= 0) hst -> SetLineColor   (lcol);
  if (lsty >= 0) hst -> SetLineStyle   (lsty);
  if (lwid >= 0) hst -> SetLineWidth   (lwid);

  if (mcol >= 0) hst -> SetMarkerColor (mcol);
  if (msty >= 0) hst -> SetMarkerStyle (msty);
  if (msiz >= 0) hst -> SetMarkerSize  (msiz);
}
//___________________________________________________________________________
string t2kana::util::TrimSpaces(string input)
{
  if( input.find_first_not_of(" \n") != 0)
    input.erase( 0,  input.find_first_not_of(" \n")  );

  if( input.find_last_not_of(" \n") != input.length() )
    input.erase( input.find_last_not_of(" \n")+1, input.length() );

  return input;
}
//____________________________________________________________________________
vector<string> t2kana::util::Tokens (string input, string delimiter)
{
// split a string of 'delimiter' separated values and return each string
// part as a vector<string> element.

  vector<string> tokens;

  while(input.find_first_of(delimiter) < input.length()) {

    tokens.push_back( input.substr(0, input.find_first_of(delimiter)) );
    input = input.substr(input.find_first_of(delimiter)+1, input.length());
  }
  tokens.push_back(input);
  return tokens;
}
//___________________________________________________________________________
void t2kana::util::PrintCompact(
        const TVectorD & v, int w, int prec, Option_t * option, ostream & stream)
{
  TString opt = option;

  bool printing_banff_params = opt.Contains("BANFF fit params");
  bool printing_fsi_params   = opt.Contains("SKDetFSI params");
  bool printing_selected_banff_params = opt.Contains("BANFF selected fit params");

  int np = v.GetNoElements();
  stream << fixed;

  for(int ip = 0; ip < np; ip++) {
     stream << "> param = " << setw(3) << ip;
     if(printing_banff_params) {
         stream << " [ ";
         PhysParamId_t p = PhysParamId_t ( (int)kNuisBANFF00 + ip );
         stream << setw(42) << util::AsString(p);
         stream << " ]";
     }
     if(printing_fsi_params) {
         stream << " [ ";
#if FSI_BINS != 3
         PhysParamId_t p = PhysParamId_t ( (int)kNuisSKDetFSI000 + ip );
#else
	 PhysParamId_t p = PhysParamId_t ( (int)kNuisSKDetFSI006 + ip );
#endif
         stream << setw(30) << util::AsString(p);
         stream << " ]";
     }
     if(printing_selected_banff_params) {
         stream << " [ ";
         PhysParamId_t p = kPhypUndefined;
	 if(ip == 0)
	   p = kNuisBANFFMAQE;
	 else if(ip == 1)
	   p = kNuisBANFFCCQEE0;
	 else if(ip == 2)
	   p = kNuisBANFFCC1piE0;
         stream << setw(30) << util::AsString(p);
         stream << " ]";
     }
     stream << " - value = " << setw(w) << setprecision(prec) << v[ip] << endl;
  }
  stream << scientific;
}
//___________________________________________________________________________
void t2kana::util::PrintCompact(
   const TVectorD & current, const TVectorD & nominal, 
   int w, int prec, Option_t * option, ostream & stream)
{
  TString opt = option;

  bool printing_banff_params = opt.Contains("BANFF fit params");
  bool printing_fsi_params   = opt.Contains("FSI params");

  int np = current.GetNoElements();
  stream << fixed;

  for(int ip = 0; ip < np; ip++) {
     stream << "> param = " << setw(3) << ip;
     if(printing_banff_params) {
         stream << " [ ";
         PhysParamId_t p = PhysParamId_t ( (int)kNuisBANFF00 + ip );
         stream << setw(42) << util::AsString(p);
         stream << " ]";
     }
     if(printing_fsi_params) {
         stream << " [ ";
#if FSI_BINS != 3
         PhysParamId_t p = PhysParamId_t ( (int)kNuisFSI000 + ip );
#else
	 PhysParamId_t p = PhysParamId_t ( (int)kNuisFSI006 + ip );
#endif
         stream << setw(30) << util::AsString(p);
         stream << " ]";
     }
     stream << " - current = " 
          << setw(w) << setprecision(prec) << current[ip];
     stream << " - nominal = " 
          << setw(w) << setprecision(prec) << nominal[ip];
     if(nominal[ip] != 0) {
       stream << " - current/nominal = " 
            << setw(w) << setprecision(prec) << current[ip]/nominal[ip];
     }
     stream << endl;
  }
  stream << scientific;
}
//___________________________________________________________________________
void t2kana::util::PrintCompact(
  const TMatrixD & m, int w, int prec, ostream & stream)
{
  int nc = m.GetNcols();
  int nr = m.GetNrows();
  stream << fixed;
  for(int ir=0; ir<nr; ir++) {
    if(ir==0) {
       stream << setw(w) << setfill(' ') << "r\\c:" << " | "; 
       for(int ic=0; ic<nc; ic++) {
         stream << setw(w) << setfill(' ') << setprecision(0) << ic << " | ";
       }
       stream << endl;
    }
    stream << setw(w) << setfill(' ') << setprecision(0) << ir << " | "; 
    for(int ic=0; ic<nc; ic++) {
         stream << setw(w) << setfill(' ') << setprecision(prec) 
              << m[ir][ic] << " | ";
    }//rows
    stream << endl;
  }//columns
  stream << scientific;
}
//___________________________________________________________________________
void t2kana::util::PrintCompactTex(
				   const TMatrixD & m, int w, int prec, ostream & stream, TString opt)
{
  int nc = m.GetNcols();
  int nr = m.GetNrows();
  stream << fixed;
  string columntype = "c";
  if(opt.Contains("banff"))
    columntype = "p{18pt}";
  stream << "{" << endl;
  if(opt.Contains("skdetfsi"))
    stream << "\\small" << endl;
  else if(opt.Contains("banff"))
    stream << "\\tiny" << endl
	   << "\\begin{sideways}" << endl
	   << "\\begin{minipage}{23cm}" << endl
	   << "{" << endl
	   << "\\setlength{\\arraycolsep}{1pt}" << endl;
  stream << "\\begin{equation}" << endl
	 << "C^{";
  if(opt.Contains("skdetfsi"))
    stream << "SK+FSI";
  else if(opt.Contains("banff"))
    stream << "banff";
  stream << "} =";
  if(opt.Contains("banff"))
    stream << "\\end{equation}" << endl
	   << "\\begin{equation*}" << endl;  
  else
    stream << " \\\\" << endl;
  stream << "\\begin{blockarray}{c ";
  for(int ic = 0; ic < nc; ic++)
    stream << columntype;
  stream << "}" << endl;

  for(int ir=0; ir<nr; ir++) {
    if(ir==0) {
      if(!opt.Contains("skdetfsi") && !opt.Contains("banff"))
	stream << setw(w) << setfill(' ') << "r\\c:"; 
      else
	stream << setw(w) << setfill(' ');
      for(int ic=0; ic<nc; ic++) {
	if(opt.Contains("skdetfsi"))
	  stream << " & " << setw(w) << setfill(' ') << setprecision(0) << "{\\color{blue} " << util::AsStringPaperTexNoDollar(kListOfSKDetFSIParams[ic]) << "}";
	else if(opt.Contains("banff"))
	  stream << " & " << setw(w) << setfill(' ') << setprecision(0) << "\\begin{sideways} {\\color{blue} " << util::AsStringPaperTex((PhysParamId_t)(ic + (int)kNuisBANFF00)) << "} \\end{sideways}";
	else
	  stream << " & " << setw(w) << setfill(' ') << setprecision(0) << "{\\color{blue} " << ic << "}";
       }
       stream << " \\\\ " << endl;
       stream << "\\begin{block}{c(";
       for(int ic = 0; ic < nc; ic++)
	 stream << columntype;
       if(opt.Contains("banff"))
	 stream << "@{\\hspace*{1pt}})}" << endl;
       else
	 stream << "@{\\hspace*{2pt}})}" << endl;
    }
    if(opt.Contains("skdetfsi"))
      stream << setw(w) << setfill(' ') << setprecision(0) << "{\\color{blue} " << util::AsStringPaperTexNoDollar(kListOfSKDetFSIParams[ir]) << "}  \\vspace{0.2cm}";
    else if(opt.Contains("banff"))
      stream << setw(w) << setfill(' ') << setprecision(0) << "{\\color{blue} " << util::AsStringPaperTexNoDollar((PhysParamId_t)(ir + (int)kNuisBANFF00)) << "}  \\vspace{0.2cm}";
    else
      stream << setw(w) << setfill(' ') << setprecision(0) << "{\\color{blue} " << ir << "}  \\vspace{0.2cm}";
    for(int ic=0; ic<nc; ic++) {
      if(ic == ir)
	stream << " & " << setw(w) << setfill(' ') << setprecision(prec) 
	       << "{\\color{red} " << m[ir][ic] << "}";
      else
	stream << " & " << setw(w) << setfill(' ') << setprecision(prec) 
	       << m[ir][ic];
    }//rows
    stream << " \\\\ " << endl;
  }//columns
  stream << "\\end{block}" << endl;
  if(opt.Contains("skdetfsi"))
    stream << "\\label{eq:Cskfsi}" << endl;
  else if(opt.Contains("banff"))
    stream << "\\label{eq:Cbanff}" << endl;
  stream << "\\end{blockarray}" << endl;
  if(opt.Contains("banff"))
    stream << "\\end{equation*}" << endl;
  else
    stream << "\\end{equation}" << endl;
  if(opt.Contains("banff"))
    stream << "}" << endl
	   << "\\end{minipage}" << endl
	   << "\\end{sideways}" << endl;
  stream << "\\normalsize" << endl
	 << "}" << endl;
  stream << scientific;
}
//___________________________________________________________________________
int t2kana::util::Sign(double tweaking_dial)
{
  if(fabs(tweaking_dial) < kASmallNum) return 0;
  if (tweaking_dial < 0) return -1;
  return +1;
}
//___________________________________________________________________________
XMLNodePointer_t t2kana::util::FindNode(
  TXMLEngine * xmleng, XMLNodePointer_t topnode, string path)
{
// Find & return the XMLNodePointer_t at the specified path
//
  vector<string> node_path = Tokens(path,"/");

  unsigned int ndepth = node_path.size();
  unsigned int idepth = 0;

  XMLNodePointer_t curr_node = xmleng->GetChild(topnode);
  while(curr_node) {
    string curr_node_name = xmleng->GetNodeName(curr_node);
    if(node_path[idepth] == curr_node_name) {
        idepth++;
        if(idepth == ndepth) {
          return curr_node;
        } else {
          curr_node = xmleng->GetChild(curr_node);
        }
    } else {
       curr_node = xmleng->GetNext(curr_node);
    }
  }
  return 0;
}
//___________________________________________________________________________
bool t2kana::util::NodeExists(
  TXMLEngine * xmleng, XMLNodePointer_t topnode, string path)
{
// Find & return the XMLNodePointer_t at the specified path
//
  vector<string> node_path = Tokens(path,"/");

  unsigned int ndepth = node_path.size();
  unsigned int idepth = 0;

  XMLNodePointer_t curr_node = xmleng->GetChild(topnode);
  while(curr_node) {
    string curr_node_name = xmleng->GetNodeName(curr_node);
    if(node_path[idepth] == curr_node_name) {
        idepth++;
        if(idepth == ndepth) {
          return true;
        } else {
          curr_node = xmleng->GetChild(curr_node);
        }
    } else {
       curr_node = xmleng->GetNext(curr_node);
    }
  }
  return false;
}
//___________________________________________________________________________
double t2kana::util::XML2Dbl(
  TXMLEngine * xmleng, XMLNodePointer_t topnode, string path)
{
// Get the specified node content as a floating-point number
//

#ifdef __show_xml_mesg_info__
  cout << "utils/XML2Dbl >> Reading info at XML node path: " << path << endl;
#endif

  XMLNodePointer_t found_node = FindNode(xmleng,topnode,path);
  if(!found_node) {
    cerr << "utils/XML2Dbl >> No node was found at input XML node path: "
         << path << endl;
    exit(1);
  }
  string content = xmleng->GetNodeContent(found_node);
  double value = atof(content.c_str());
  return value;
}
//___________________________________________________________________________
vector<double> t2kana::util::XML2DblArray(
  TXMLEngine * xmleng, XMLNodePointer_t topnode, string path)
{
// Get the specified node content as a vector floating-point number
// The node content should be a comma-separed list of numbers

#ifdef __show_xml_mesg_info__
  cout << "utils/XML2Dbl >> Reading info at XML node path: " << path << endl;
#endif

  XMLNodePointer_t found_node = FindNode(xmleng,topnode,path);
  if(!found_node) {
    cerr << "utils/XML2Dbl >> No node was found at input XML node path: "
         << path << endl;
    exit(1);
  }
  string content = TrimSpaces(xmleng->GetNodeContent(found_node));

  vector<double> array;
  vector<string> str_tokens = Tokens(content, ",");
  vector<string>::const_iterator str_tokens_iter = str_tokens.begin();
  for( ; str_tokens_iter != str_tokens.end(); ++str_tokens_iter) {
     string token = *str_tokens_iter;
     array.push_back( atof(token.c_str()) );
  }
  return array;
}
//___________________________________________________________________________
string t2kana::util::XML2String(
  TXMLEngine * xmleng, XMLNodePointer_t topnode, string path)
{
// Get the specified node content as a floating-point number
//

#ifdef __show_xml_mesg_info__
  cout << "utils/XML2Str >> Reading info at XML node path: " << path << endl;
#endif

  XMLNodePointer_t found_node = FindNode(xmleng,topnode,path);
  if(!found_node) {
    cerr << "utils/XML2Str >> No node was found at input XML node path: "
         << path << endl;
    exit(1);
  }
  string content = xmleng->GetNodeContent(found_node);
  return TrimSpaces(content);
}
//___________________________________________________________________________
string t2kana::util::XML2RootFileName(
   TXMLEngine * xmleng, XMLNodePointer_t topnode, string path)
{
  string filename = XML2String(xmleng,topnode,path+"/filename");
#ifdef __show_xml_mesg_info__
  cout << "utils/XML2RootFileName >> filename : " << filename << endl;
#endif
  return filename;
}
//___________________________________________________________________________
string t2kana::util::XML2RootObjName (
   TXMLEngine * xmleng, XMLNodePointer_t topnode, string path)
{
  string objname = XML2String(xmleng,topnode,path+"/objname");
#ifdef __show_xml_mesg_info__
  cout << "utils/XML2RootObjName >> objname  : " << objname  << endl;
#endif
  return objname;
}
//___________________________________________________________________________
string t2kana::util::XML2RootObjType (
   TXMLEngine * xmleng, XMLNodePointer_t topnode, string path)
{
  string objtype = XML2String(xmleng,topnode,path+"/objtype");
#ifdef __show_xml_mesg_info__
  cout << "utils/XML2RootObjType >> objtype  : " << objtype  << endl;
#endif
  return objtype;
}
//___________________________________________________________________________
TFile * t2kana::util::XML2TFile(
  TXMLEngine * xmleng, XMLNodePointer_t topnode, string path)
{
#ifdef __show_xml_mesg_info__
  cout << "utils/XML2TFile >> Reading info at XML node path: " << path << endl;
#endif

  string filename = XML2RootFileName(xmleng,topnode,path);

  TFile * file = new TFile(
    Form("%s/%s", gSystem->Getenv("T2KOSCANA"), filename.c_str()), "read");
  if(!file) {
    cerr << "utils/XML2TFile >> No file: " << filename << " found" << endl;
    exit(1);
  }
  if(file->IsZombie()) {
    cerr << "utils/XML2TFile >> File is a zombie: " << filename << endl;
    exit(1);
  }

  return file;
}
//___________________________________________________________________________
TH1F * t2kana::util::XML2TH1F(
  TXMLEngine * xmleng, XMLNodePointer_t topnode, string path)
{
// Get a TH1F instance from the info in the input XML config file
//
#ifdef __show_xml_mesg_info__
  cout << "utils/XML2TH1F >> Reading info at XML node path: " << path << endl;
#endif

  string filename = XML2RootFileName(xmleng,topnode,path);
  string objname  = XML2RootObjName (xmleng,topnode,path);

  TFile file(
    Form("%s/%s", gSystem->Getenv("T2KOSCANA"), filename.c_str()), "read");
  TH1F * h = (TH1F*)file.Get(objname.c_str());
  if(!h) {
    cerr << "utils/XML2TH1F >> No " << objname << " in " << filename << endl;
    file.Close();
    //return 0;
    exit(1);
  }

  TH1F * hcopy = new TH1F(*h);
  hcopy->SetDirectory(0);

  file.Close();

  return hcopy;
}
//___________________________________________________________________________
TH1D * t2kana::util::XML2TH1D(
  TXMLEngine * xmleng, XMLNodePointer_t topnode, string path)
{
// Get a TH1D instance from the info in the input XML config file
//
#ifdef __show_xml_mesg_info__
  cout << "utils/XML2TH1D >> Reading info at XML node path: " << path << endl;
#endif

  string filename = XML2RootFileName(xmleng,topnode,path);
  string objname  = XML2RootObjName (xmleng,topnode,path);

  TFile file(
    Form("%s/%s", gSystem->Getenv("T2KOSCANA"), filename.c_str()), "read");
  TH1D * h = (TH1D*)file.Get(objname.c_str());
  if(!h) {
    cerr << "utils/XML2TH1D >> No " << objname << " in " << filename << endl;
    file.Close();
    //return 0;
    exit(1);
  }

  TH1D * hcopy = new TH1D(*h);
  hcopy->SetDirectory(0);

  file.Close();

  return hcopy;
}
//___________________________________________________________________________
TH2D * t2kana::util::XML2TH2D(
  TXMLEngine * xmleng, XMLNodePointer_t topnode, string path)
{
// Get a TH2D instance from the info in the input XML config file
//
#ifdef __show_xml_mesg_info__
  cout << "utils/XML2TH2D >> Reading info at XML node path: " << path << endl;
#endif

  string filename = XML2RootFileName(xmleng,topnode,path);
  string objname  = XML2RootObjName (xmleng,topnode,path);

  TFile file(
     Form("%s/%s", gSystem->Getenv("T2KOSCANA"), filename.c_str()), "read");    
  TH2D * h2 = (TH2D*)file.Get(objname.c_str());
  if(!h2) {
    cerr << "utils/XML2TH2D >> No " << objname << " in " << filename << endl;
    file.Close();
    //return 0;
    exit(1);
  }

  TH2D * h2copy = new TH2D(*h2);
  h2copy->SetDirectory(0);

  file.Close();

  return h2copy;
}
//___________________________________________________________________________
TTree * t2kana::util::XML2TTree(
   TXMLEngine * xmleng, XMLNodePointer_t topnode, string path)
{
// Get a TH2D instance from the info in the input XML config file
//
#ifdef __show_xml_mesg_info__
  cout << "utils/XML2TTree >> Reading info at XML node path: " << path << endl;
#endif

  string filename = XML2RootFileName(xmleng,topnode,path);
  string objname  = XML2RootObjName (xmleng,topnode,path);

/////  TFile file(filename.c_str(), "read");
/////  TTree * tree = (TTree*)file.Get(objname.c_str());

  TFile * file = new TFile(
     Form("%s/%s", gSystem->Getenv("T2KOSCANA"), filename.c_str()), "read");
  TTree * tree = (TTree*)file->Get(objname.c_str());
  if(!tree) {
    cerr << "utils/XML2TTree >> No " << objname << " in " << filename << endl;
//    file.Close();
    //return 0;
    exit(1);
  }

  return tree;
//  TTree * treecopy = tree->CloneTree();
//  file.Close();
//  return treecopy;
}
//___________________________________________________________________________
TVectorD * t2kana::util::XML2TVectorD(
   TXMLEngine * xmleng, XMLNodePointer_t topnode, string path)
{
// Get a TVectorD instance from the info in the input XML config file
//
#ifdef __show_xml_mesg_info__
  cout << "utils/XML2TVectorD >> "
       << "Reading info at XML node path: " << path << endl;
#endif

  string filename = XML2RootFileName(xmleng,topnode,path);
  string objname  = XML2RootObjName (xmleng,topnode,path);

  TFile file(
    Form("%s/%s", gSystem->Getenv("T2KOSCANA"), filename.c_str()), "read");
  TVectorD * vec = (TVectorD*)file.Get(objname.c_str());
  if(!vec) {
    cerr << "utils/XML2TVectorD >> No " 
         << objname << " in " << filename << endl;
    file.Close();
    //return 0;
    exit(1);
  }

  TVectorD * vec_copy = new TVectorD(*vec);

  file.Close();

  return vec_copy;
}
//___________________________________________________________________________
TMatrixDSym * t2kana::util::XML2TMatrixDSym(
   TXMLEngine * xmleng, XMLNodePointer_t topnode, string path)
{
// Get a TMatrixDSym instance from the info in the input XML config file
//
  string filename = XML2RootFileName(xmleng,topnode,path);
  string objname  = XML2RootObjName (xmleng,topnode,path);

  TFile file(
    Form("%s/%s", gSystem->Getenv("T2KOSCANA"), filename.c_str()), "read");
  TMatrixDSym * matrix = (TMatrixDSym*)file.Get(objname.c_str());
  if(!matrix) {
    cerr << "utils/XML2TMatrixDSym >> No " 
         << objname << " in " << filename << endl;
    file.Close();
    //return 0;
    exit(1);
  }

  TMatrixDSym * matrix_copy = new TMatrixDSym(*matrix);

  file.Close();

  return matrix_copy;
}
//___________________________________________________________________________
TMatrixD * t2kana::util::XML2TMatrixD(
    TXMLEngine * xmleng, XMLNodePointer_t topnode, string path)
{
// Get a TMatrixD instance from the info in the input XML config file
//
  string objtype  = XML2RootObjType (xmleng,topnode,path);
  if(objtype == "TH2D")
    return XML2TH2D2TMatrixD(xmleng,topnode,path);

  string filename = XML2RootFileName(xmleng,topnode,path);
  string objname  = XML2RootObjName (xmleng,topnode,path);

  TFile * file = new TFile( 
    Form("%s/%s", gSystem->Getenv("T2KOSCANA"), filename.c_str()), "read");
  TMatrixD * matrix = (TMatrixD*) file->Get(objname.c_str());
  if(!matrix) {
    cerr << "utils/XML2TMatrixD >> No " 
         << objname << " in " << filename << endl;
    file->Close();
    //return 0;
    exit(1);
  }

  TMatrixD * matrix_copy = new TMatrixD(*matrix);
  file->Close();

  return matrix_copy;
}
//___________________________________________________________________________
TMatrixD * t2kana::util::XML2TH2D2TMatrixD(
  TXMLEngine * xmleng, XMLNodePointer_t topnode, string path)
{
// Get a TH2D instance from the info in the input XML config file
//
#ifdef __show_xml_mesg_info__
  cout << "utils/XML2TH2D2TMatrixD >> Reading info at XML node path: " << path << endl;
#endif

  string filename = XML2RootFileName(xmleng,topnode,path);
  string objname  = XML2RootObjName (xmleng,topnode,path);

  TFile file(
     Form("%s/%s", gSystem->Getenv("T2KOSCANA"), filename.c_str()), "read");    
  TH2D * h2 = (TH2D*)file.Get(objname.c_str());
  if(!h2) {
    cerr << "utils/XML2TH2D2TMatrixD >> No " << objname << " in " << filename << endl;
    file.Close();
    //return 0;
    exit(1);
  }

  TH2D * h2copy = new TH2D(*h2);
  h2copy->SetDirectory(0);

  file.Close();

  const int nrows = h2copy->GetNbinsX();
  const int ncols = h2copy->GetNbinsY();

  TMatrixD matrix(nrows,ncols);
  for(int i = 1; i <= nrows; i++) {
    for(int j = 1; j <= ncols; j++) {
      matrix(i-1,j-1)=h2copy->GetBinContent(i,j);
#ifdef __show_xml_mesg_debug__
      cout << i << "\t" << j << "\t" << h2copy->GetBinContent(i,j) << endl;
#endif
    }
  }

  TMatrixD * matrix_copy = new TMatrixD(matrix);
#ifdef __show_xml_mesg_debug__
  matrix_copy->Print();
#endif

  return matrix_copy;
}
//___________________________________________________________________________
TSpline3 * t2kana::util::XML2TSpline3(
  TXMLEngine * xmleng, XMLNodePointer_t topnode, string path)
{
  string filename = XML2RootFileName(xmleng,topnode,path);
  string objname  = XML2RootObjName (xmleng,topnode,path);

  TFile * file = new TFile( 
    Form("%s/%s", gSystem->Getenv("T2KOSCANA"), filename.c_str()), "read");
  TSpline3 * spl = (TSpline3*) file->Get(objname.c_str());
  if(!spl) {
    cerr << "utils/XML2TSpline3 >> No " 
         << objname << " in " << filename << endl;
    file->Close();
    //return 0;
    exit(1);
  }

  TSpline3 * spl_copy = new TSpline3(*spl);
  file->Close();

  return spl_copy;
}
//___________________________________________________________________________
TAxis *  t2kana::util::XML2BinInfo(
  TXMLEngine * xmleng, XMLNodePointer_t topnode, string path)
{
//
//
#ifdef __show_xml_mesg_info__
  cout << "utils/XML2BinInfo >> Reading info at XML node path: " << path << endl;
#endif

  XMLNodePointer_t found_node = FindNode(xmleng,topnode,path);
  if(!found_node) {
    cerr << "utils/XML2BinInfo >> No node was found at input XML node path: "
         << path << endl;
    //return 0;
    exit(1);
  }

  string binsz = xmleng->GetAttr(found_node,"bin_size");
  if(binsz == "fixed") {
    int    nbins = 0;
    double min   = 0;
    double max   = 0;
    XMLNodePointer_t child_node = xmleng->GetChild(found_node);
    while(child_node) {
      string child_node_name = xmleng->GetNodeName(child_node);
      string content = xmleng->GetNodeContent(child_node);
      if(child_node_name == "nbins") { nbins = atoi(TrimSpaces(content).c_str()); }
      if(child_node_name == "min"  ) { min   = atof(TrimSpaces(content).c_str()); }
      if(child_node_name == "max"  ) { max   = atof(TrimSpaces(content).c_str()); }
      child_node = xmleng->GetNext(child_node);    
    }
#ifdef __show_xml_mesg_info__
    cout << "utils/XML2BinInfo >> nbins = " << nbins 
         << ", range = (" << min << ", " << max << ")" << endl;
#endif
    assert(nbins>0 && min<max);

    TAxis * axis = new TAxis(nbins,min,max);
    return axis;
  } 
  else   
  if(binsz == "variable") {
    cerr << "utils/XML2BinInfo >> Not handling variable size bins yet" << endl;
    exit(1);
  }
  else {
    cerr << "utils/XML2BinInfo >> Invalid bin_size attribute" << endl;
    exit(1);
  }

  return 0;
}
//___________________________________________________________________________

// A few helper functions for use in UseEventToFillInTemplate
// that make swapping mode definitions a little easier

bool t2kana::util::IsCCQE(int neut_reaction_mode) 
{
  return neut_reaction_mode ==  1;
}
//___________________________________________________________________________
bool t2kana::util::IsCC1pi(int neut_reaction_mode) 
{
  bool is_cc1pi = (neut_reaction_mode == 11 || 
		   neut_reaction_mode == 12 || 
		   neut_reaction_mode == 13);
  return is_cc1pi;
}
//___________________________________________________________________________
bool t2kana::util::IsCCcoh(int neut_reaction_mode) 
{
  return neut_reaction_mode ==  16;
}
//___________________________________________________________________________
bool t2kana::util::IsCCoth(int neut_reaction_mode) 
//this is a catch all for everything CC that is not QE, 1pi, or coh
// needs to be updated if another subcategory is added
{
  bool is_ccoth = (neut_reaction_mode == 17 ||
		     neut_reaction_mode == 21 ||
		     neut_reaction_mode == 22 || 
		     neut_reaction_mode == 23 || 
		     neut_reaction_mode == 26);
  return is_ccoth;
}
//___________________________________________________________________________
//bool t2kana::util::IsCCoth(int neut_reaction_mode) 
//{
//return IsCCoth(neut_reaction_mode);
//}
//___________________________________________________________________________
//bool t2kana::util::IsCCoth(int neut_reaction_mode) 
//{
//return IsCCoth(neut_reaction_mode);
//}
//___________________________________________________________________________
bool t2kana::util::IsNC1pi0(int neut_reaction_mode) 
{
  bool is_cc1pi0 = (neut_reaction_mode == 31 ||
		    neut_reaction_mode == 32);
  return is_cc1pi0;
}
//___________________________________________________________________________
bool t2kana::util::IsNC1piC(int neut_reaction_mode) 
{
  bool is_cc1pi0 = (neut_reaction_mode == 33 ||
		    neut_reaction_mode == 34);
  return is_cc1pi0;
}
//___________________________________________________________________________
bool t2kana::util::IsNCcoh(int neut_reaction_mode) 
{
  return neut_reaction_mode == 36;
}
//___________________________________________________________________________
bool t2kana::util::IsNCoth(int neut_reaction_mode) 
//this is a catch all for everything NC that is not 1pi0, 1piC or coh
// needs to be updated if another subcategory is added
{
  bool is_ncoth = (  neut_reaction_mode == 38 ||
		     neut_reaction_mode == 39 ||
		     neut_reaction_mode == 41 ||
		     neut_reaction_mode == 42 ||
		     neut_reaction_mode == 43 ||
		     neut_reaction_mode == 44 ||
		     neut_reaction_mode == 45 ||
		     neut_reaction_mode == 46 ||
		     neut_reaction_mode == 51 ||
		     neut_reaction_mode == 52);
  return is_ncoth;
}
//___________________________________________________________________________
bool t2kana::util::IsNCoth_mu(int neut_reaction_mode) 
{
  bool is_ncoth_mu = ( IsNCoth (neut_reaction_mode) ||
		       IsNCcoh (neut_reaction_mode) ||
		       IsNC1pi0(neut_reaction_mode) );
  return is_ncoth_mu;
}
//___________________________________________________________________________
bool t2kana::util::IsNCoth_e(int neut_reaction_mode) 
{
  bool is_ncoth_e = ( IsNCoth (neut_reaction_mode) ||
		      IsNC1piC(neut_reaction_mode) );
  return is_ncoth_e;
}
//___________________________________________________________________________
bool t2kana::util::UseEventToFillInTemplate(
  SuperkMode_t m, int nu_pdg_code, int neut_reaction_mode)
{

  bool use_it = false;

  //NuMu templates
  if(m == kSKmNuMuCCQE) 
  {
     bool is_numu = (nu_pdg_code == kPdgNuMu);
     use_it = is_numu && IsCCQE(neut_reaction_mode);
  } 

  else
  if(m == kSKmNuMuCC1pi) 
  {
     bool is_numu  = (nu_pdg_code == kPdgNuMu);
     use_it = is_numu && IsCC1pi(neut_reaction_mode);
  } 

  else
  if(m == kSKmNuMuCCcoh)
  {
    bool is_numu  = (nu_pdg_code == kPdgNuMu);
    use_it = is_numu && IsCCcoh(neut_reaction_mode);
  }

  else
  if(m == kSKmNuMuCCoth) 
  {
     bool is_numu     = (nu_pdg_code == kPdgNuMu);
     use_it = is_numu && IsCCoth(neut_reaction_mode);
  } 

  else
  if(m == kSKmNuMuCCoth)
  {
    bool is_numu    = (nu_pdg_code == kPdgNuMu);
    use_it = is_numu && IsCCoth(neut_reaction_mode);
  }

  else
  if(m == kSKmNuMuNC1piC) 
  {
     bool is_numu   = (nu_pdg_code == kPdgNuMu);
     use_it = is_numu && IsNC1piC(neut_reaction_mode);
  } 

  else
  if(m == kSKmNuMuNC1pi0) 
  {
     bool is_numu   = (nu_pdg_code == kPdgNuMu);
     use_it = is_numu && IsNC1pi0(neut_reaction_mode);
  } 

  else
  if(m == kSKmNuMuNCcoh) 
  {
     bool is_numu  = (nu_pdg_code == kPdgNuMu);
     use_it = is_numu && IsNCcoh(neut_reaction_mode);
  } 

  else
  if(m == kSKmNuMuNCoth_e) 
  {
     bool is_numu    = (nu_pdg_code == kPdgNuMu);
     use_it = is_numu && IsNCoth_e(neut_reaction_mode);
  } 

  else
  if(m == kSKmNuMuNCoth_mu) 
  {
     bool is_numu    = (nu_pdg_code == kPdgNuMu);
     use_it = is_numu && IsNCoth_mu(neut_reaction_mode);
  } 

  //NuMuBar templates
  else
  if(m == kSKmNuMuBarCCQE) 
  {
     bool is_numubar = (nu_pdg_code == kPdgNuMuBar);
     use_it = is_numubar && IsCCQE(-neut_reaction_mode);
  } 

  else
  if(m == kSKmNuMuBarCC1pi) 
  {
     bool is_numubar = (nu_pdg_code == kPdgNuMuBar);
     use_it = is_numubar && IsCC1pi(-neut_reaction_mode);
  } 

  else
  if(m == kSKmNuMuBarCCcoh) 
  {
     bool is_numubar = (nu_pdg_code == kPdgNuMuBar);
     use_it = is_numubar && IsCCcoh(-neut_reaction_mode);
  } 

  else
  if(m == kSKmNuMuBarCCoth) 
  {
     bool is_numubar  = (nu_pdg_code == kPdgNuMuBar);
     use_it = is_numubar && IsCCoth(-neut_reaction_mode);
  } 

  else
  if(m == kSKmNuMuBarCCoth) 
  {
     bool is_numubar = (nu_pdg_code == kPdgNuMuBar);
     use_it = is_numubar && IsCCoth(-neut_reaction_mode);
  } 

  else
  if(m == kSKmNuMuBarNC1piC) 
  {
     bool is_numubar  = (nu_pdg_code == kPdgNuMuBar);
     use_it = is_numubar && IsNC1piC(-neut_reaction_mode);
  } 

  else
  if(m == kSKmNuMuBarNC1pi0) 
  {
     bool is_numubar  = (nu_pdg_code == kPdgNuMuBar);
     use_it = is_numubar && IsNC1pi0(-neut_reaction_mode);
  } 

  else
  if(m == kSKmNuMuBarNCcoh) 
  {
     bool is_numubar  = (nu_pdg_code == kPdgNuMuBar);
     use_it = is_numubar && IsNCcoh(-neut_reaction_mode);
  } 

  else
  if(m == kSKmNuMuBarNCoth_e) 
  {
     bool is_numubar  = (nu_pdg_code == kPdgNuMuBar);
     use_it = is_numubar && IsNCoth_e(-neut_reaction_mode);
  } 

  else
  if(m == kSKmNuMuBarNCoth_mu) 
  {
     bool is_numubar  = (nu_pdg_code == kPdgNuMuBar);
     use_it = is_numubar && IsNCoth_mu(-neut_reaction_mode);
  } 

  //NuE (signal and beam) templates
  else
  if(m == kSKmNuECCQE || m == kSKmOscNuECCQE) 
  {
     bool is_nue  = (nu_pdg_code == kPdgNuE); 
     use_it = is_nue && IsCCQE(neut_reaction_mode);
  } 

  else
  if(m == kSKmNuECC1pi || m == kSKmOscNuECC1pi) 
  {
     bool is_nue   = (nu_pdg_code == kPdgNuE); 
     use_it = is_nue && IsCC1pi(neut_reaction_mode);
  }

  else
  if(m == kSKmNuECCcoh || m == kSKmOscNuECCcoh) 
  {
     bool is_nue   = (nu_pdg_code == kPdgNuE); 
     use_it = is_nue && IsCCcoh(neut_reaction_mode);
  } 

  else
  if(m == kSKmNuECCoth || m == kSKmOscNuECCoth) 
  {
     bool is_nue      = (nu_pdg_code == kPdgNuE); 
     use_it = is_nue && IsCCoth(neut_reaction_mode);
  }

  else
  if(m == kSKmNuECCoth || m == kSKmOscNuECCoth) 
  {
     bool is_nue     = (nu_pdg_code == kPdgNuE); 
     use_it = is_nue && IsCCoth(neut_reaction_mode);
  }

  //this is a little different. nuesig has no CCcoh in 1Re sample, whereas CCcoh is included
  // for the other species
  //  else
  //if(m == kSKmOscNuECCoth) 
  //{
  // bool is_nue     = (nu_pdg_code == kPdgNuE); 
  // use_it = is_nue && (IsCCoth(neut_reaction_mode) || IsCCcoh(neut_reaction_mode));
  //}

  else
  if(m == kSKmNuENC1piC) 
  {
     bool is_nue  = (nu_pdg_code == kPdgNuE); 
     use_it = is_nue && IsNC1piC(neut_reaction_mode);
  }

  else
  if(m == kSKmNuENC1pi0) 
  {
     bool is_nue    = (nu_pdg_code == kPdgNuE); 
     use_it = is_nue && IsNC1pi0(neut_reaction_mode);
  }

  else
  if(m == kSKmNuENCcoh) 
  {
     bool is_nue   = (nu_pdg_code == kPdgNuE); 
     use_it = is_nue && IsNCcoh(neut_reaction_mode);
  }

  else
  if(m == kSKmNuENCoth_e) 
  {
     bool is_nue     = (nu_pdg_code == kPdgNuE); 
     use_it = is_nue && IsNCoth_e(neut_reaction_mode);
  }

  else
  if(m == kSKmNuENCoth_mu) 
  {
     bool is_nue     = (nu_pdg_code == kPdgNuE); 
     use_it = is_nue && IsNCoth_mu(neut_reaction_mode);
  }

  //NuEBar (beam) templates
  else
  if(m == kSKmNuEBarCCQE) 
  {
     bool is_nue  = (nu_pdg_code == kPdgNuEBar); 
     use_it = is_nue && IsCCQE(-neut_reaction_mode);
  } 

  else
  if(m == kSKmNuEBarCC1pi) 
  {
     bool is_nue   = (nu_pdg_code == kPdgNuEBar); 
     use_it = is_nue && IsCC1pi(-neut_reaction_mode);
  }

  else
  if(m == kSKmNuEBarCCcoh) 
  {
     bool is_nue   = (nu_pdg_code == kPdgNuEBar); 
     use_it = is_nue && IsCCcoh(-neut_reaction_mode);
  } 


  else
  if(m == kSKmNuEBarCCoth) 
  {
     bool is_nue      = (nu_pdg_code == kPdgNuEBar); 
     use_it = is_nue && IsCCoth(-neut_reaction_mode);
  }

  else
  if(m == kSKmNuEBarCCoth) 
  {
     bool is_nue     = (nu_pdg_code == kPdgNuEBar); 
     use_it = is_nue && IsCCoth(-neut_reaction_mode);
  }

  else
  if(m == kSKmNuEBarNC1piC) 
  {
     bool is_nue  = (nu_pdg_code == kPdgNuEBar); 
     use_it = is_nue && IsNC1piC(-neut_reaction_mode);
  }

  else
  if(m == kSKmNuEBarNC1pi0) 
  {
     bool is_nue    = (nu_pdg_code == kPdgNuEBar); 
     use_it = is_nue && IsNC1pi0(-neut_reaction_mode);
  }

  else
  if(m == kSKmNuEBarNCcoh) 
  {
     bool is_nue   = (nu_pdg_code == kPdgNuEBar); 
     use_it = is_nue && IsNCcoh(-neut_reaction_mode);
  }

  else
  if(m == kSKmNuEBarNCoth_e) 
  {
     bool is_nue     = (nu_pdg_code == kPdgNuEBar); 
     use_it = is_nue && IsNCoth_e(-neut_reaction_mode);
  }

  else
  if(m == kSKmNuEBarNCoth_mu) 
  {
     bool is_nue     = (nu_pdg_code == kPdgNuEBar); 
     use_it = is_nue && IsNCoth_mu(-neut_reaction_mode);
  }

  return use_it;
}
//___________________________________________________________________________
double t2kana::util::NuMuSurvProbOsc2Flavour(
  double mixing23, double dm23sq, double Et, bool unphysical_osc_prob = false)
{
  if(Et <= 0.) return 0.;

  double sinsq2theta23 = 0.;

  NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
  switch(omp) {
   case(kNuMixSinsq2Theta):
       // mixing23 parameter is sinsq2theta23
       sinsq2theta23 = mixing23;
       break;
   case(kNuMixSinsqTheta):
       // mixing23 parameter is sinsqtheta23
       // convert to sinsq2theta23
       sinsq2theta23 = 4*mixing23*(1-mixing23);
       break;
   default:
       return 0.;
  }

  double osc   = TMath::Sin(1.267*dm23sq*kBaseline/Et);
  double posc  = sinsq2theta23 * TMath::Power(osc,2);
  double psurv = 1-posc;

  if(!unphysical_osc_prob) {
    psurv = TMath::Min(1.,psurv);
    psurv = TMath::Max(0.,psurv);
  }

  return psurv;
}
//____________________________________________________________________________
double t2kana::util::NuMuSurvProbOsc2Flavour(
  double mixing23, double dm23sq, double Et, t2kana::conv::SuperkMode_t m,
  bool unphysical_osc_prob = false)
{
// Note: Using only the numu disappearance analysis modes
//
  if(m == kSKmNuMuNC1piC      || 
     m == kSKmNuMuNC1pi0      ||
     m == kSKmNuMuNCcoh       ||
     m == kSKmNuMuNCoth_mu    ||
     m == kSKmNuMuNCoth_e     ||
     m == kSKmNuMuBarNC1piC   ||
     m == kSKmNuMuBarNC1pi0   ||
     m == kSKmNuMuBarNCcoh    ||
     m == kSKmNuMuBarNCoth_mu ||
     m == kSKmNuMuBarNCoth_e  ||
     m == kSKmNuECCQE         ||
     m == kSKmNuECC1pi        ||
     m == kSKmNuECCcoh        ||
     m == kSKmNuECCoth     ||
     //     m == kSKmNuECCoth      ||
     m == kSKmNuENC1piC       ||
     m == kSKmNuENC1pi0       ||
     m == kSKmNuENCcoh        ||
     m == kSKmNuENCoth_mu     ||
     m == kSKmNuENCoth_e      ||
     m == kSKmNuEBarCCQE      ||
     m == kSKmNuEBarCC1pi     ||
     m == kSKmNuEBarCCcoh     ||
     //     m == kSKmNuEBarCCoth  ||
     m == kSKmNuEBarCCoth   ||
     m == kSKmNuEBarNC1piC    ||
     m == kSKmNuEBarNC1pi0    ||
     m == kSKmNuEBarNCcoh     ||
     m == kSKmNuEBarNCoth_mu  ||
     m == kSKmNuEBarNCoth_e   )
  {
    return 1.;
  }
  else
  if(m == kSKmNuMuCCQE        || 
     m == kSKmNuMuCC1pi       || 
     m == kSKmNuMuCCcoh       || 
     //     m == kSKmNuMuCCoth    ||
     m == kSKmNuMuCCoth     ||
     m == kSKmNuMuBarCCQE     ||
     m == kSKmNuMuBarCC1pi    ||
     m == kSKmNuMuBarCCcoh    ||
     //     m == kSKmNuMuBarCCoth || 
     m == kSKmNuMuBarCCoth  ||
     m == kSKmMEC             )
  {
    return NuMuSurvProbOsc2Flavour(mixing23, dm23sq, Et, unphysical_osc_prob);
  }
  else 
  if(
     m == kSKmOscNuMuCCQE     ||
     m == kSKmOscNuMuCC1pi    ||
     m == kSKmOscNuMuCCcoh    ||
     //     m == kSKmOscNuMuCCoth ||
     m == kSKmOscNuMuCCoth  ||
     m == kSKmOscNuMuBarCCQE     ||
     m == kSKmOscNuMuBarCC1pi    ||
     m == kSKmOscNuMuBarCCcoh    ||
     //     m == kSKmOscNuMuCCoth ||
     m == kSKmOscNuMuBarCCoth  ||
     m == kSKmOscNuEBarCCQE     ||
     m == kSKmOscNuEBarCC1pi    ||
     m == kSKmOscNuEBarCCcoh    ||
     //     m == kSKmOscNuMuCCoth ||
     m == kSKmOscNuEBarCCoth 
)
  {
    return 0;
  }
  else {
    cerr << "NuMuSurvProbOsc2D > Can't handle SK mode = " << m << endl;
    exit(1);
  }
}
//____________________________________________________________________________
double t2kana::util::NuMuSurvProbDecay (
   double sinsq2theta, double alpha, double Et)
{
// Note: In oscillation formula (Dm^2 * L / 4E) in natural units becomes 
// 1.267 * Dm^2 * L / E where Dm^2 in eV2, L in km, E in GeV. 
// Apply the same conversion in the decoherence survival prob formula (Eq. 13) 
// given in Barger et al., astro-ph/9810121v1.
//
// a L / 2E (natural units) --> 2 * 1.267 a (eV^2) L(km) /E (GeV)
//

  double cossq2theta = TMath::Sqrt( TMath::Max(0.,1.-sinsq2theta) );
  double cos2theta   = TMath::Sqrt(cossq2theta);
  double sinsqtheta  = (1.-cos2theta)/2.;
  double cossqtheta  = (1.+cos2theta)/2.;
  double expaL       = TMath::Exp(-2*1.267*alpha*kBaseline/Et);

  double psurv = TMath::Power(sinsqtheta + cossqtheta*expaL, 2.);
  return psurv;
}
//____________________________________________________________________________
double t2kana::util::NuMuSurvProbDecay (
  double sinsq2theta, double alpha, double Et, 
  t2kana::conv::SuperkMode_t m)
{
// Note: Using only the numu disappearance analysis modes
//
  if(m == kSKmNuECCQE        ||
     m == kSKmNuECC1pi       ||
     m == kSKmNuECCcoh       ||
     m == kSKmNuECCoth    ||
     m == kSKmNuENC1piC      ||
     m == kSKmNuENCoth_mu    ||
     m == kSKmNuEBarCCQE     ||
     m == kSKmNuEBarCC1pi    ||
     m == kSKmNuEBarCCcoh    ||
     m == kSKmNuEBarCCoth ||
     m == kSKmNuEBarNC1piC   ||
     m == kSKmNuEBarNCoth_mu )
  {
    return 1.;
  }
  else
  if(m == kSKmNuMuCCQE        || 
     m == kSKmNuMuCC1pi       || 
     m == kSKmNuMuCCcoh       || 
     m == kSKmNuMuCCoth    ||
     m == kSKmNuMuNC1piC      ||
     m == kSKmNuMuNCoth_mu    ||
     m == kSKmNuMuBarCCQE     ||
     m == kSKmNuMuBarCC1pi    ||
     m == kSKmNuMuBarCCcoh    ||
     m == kSKmNuMuBarCCoth ||
     m == kSKmNuMuBarNC1piC   ||
     m == kSKmNuMuBarNCoth_mu ||
     m == kSKmMEC             )
  {
    return NuMuSurvProbDecay(sinsq2theta, alpha, Et);
  }
  else 
  if(
m == kSKmOscNuMuCCQE     ||
     m == kSKmOscNuMuCC1pi    ||
     m == kSKmOscNuMuCCcoh    ||
//     m == kSKmOscNuMuCCoth ||                                                                                                                                 
     m == kSKmOscNuMuCCoth  ||
     m == kSKmOscNuMuBarCCQE     ||
     m == kSKmOscNuMuBarCC1pi    ||
     m == kSKmOscNuMuBarCCcoh    ||
//     m == kSKmOscNuMuCCoth ||                                                                                                                                 
     m == kSKmOscNuMuBarCCoth  ||
     m == kSKmOscNuEBarCCQE     ||
     m == kSKmOscNuEBarCC1pi    ||
     m == kSKmOscNuEBarCCcoh    ||
//     m == kSKmOscNuMuCCoth ||                                                                                                                                 
     m == kSKmOscNuEBarCCoth
     )
  {
    return 0;
  }
  else {
    cerr << "NuMuSurvProbDecay > Can't handle SK mode = " << m << endl;
    exit(1);
  }
}
//____________________________________________________________________________
double t2kana::util::NuMuSurvProbDecoh (
   double sinsq2theta, double musq, double Et)
{
//
// Note: In oscillation formula (Dm^2 * L / 4E) in natural units becomes 
// 1.267 * Dm^2 * L / E where Dm^2 in eV2, L in km, E in GeV. 
// Apply the same conversion in the decoherence survival prob formula (Eq. 5) 
// given in Fogli et al., PRD67,093006(2003).
//
// \mu^2 L / 2 E (natural units) --> 2 * 1.267 \mu^2 (eV^2) L(km) /E (GeV)
//

  double expmusqLE = TMath::Exp(-2*1.267*musq*kBaseline/Et);
  double psurv = 1. - 0.5*sinsq2theta * ( 1 - expmusqLE );
  return psurv;
}
//____________________________________________________________________________
double t2kana::util::NuMuSurvProbDecoh (
  double sinsq2theta, double musq, double Et, 
  t2kana::conv::SuperkMode_t m)
{
// Note: Using only the numu disappearance analysis modes
//
  if(m == kSKmNuECCQE        ||
     m == kSKmNuECC1pi       ||
     m == kSKmNuECCcoh       ||
     m == kSKmNuECCoth    ||
     m == kSKmNuENC1piC      ||
     m == kSKmNuENCoth_mu    ||
     m == kSKmNuEBarCCQE     ||
     m == kSKmNuEBarCC1pi    ||
     m == kSKmNuEBarCCcoh    ||
     m == kSKmNuEBarCCoth ||
     m == kSKmNuEBarNC1piC   ||
     m == kSKmNuEBarNCoth_mu )
  {
    return 1.;
  }
  else
  if(m == kSKmNuMuCCQE        || 
     m == kSKmNuMuCC1pi       || 
     m == kSKmNuMuCCcoh       || 
     m == kSKmNuMuCCoth    ||
     m == kSKmNuMuNC1piC      ||
     m == kSKmNuMuNCoth_mu    ||
     m == kSKmNuMuBarCCQE     ||
     m == kSKmNuMuBarCC1pi    ||
     m == kSKmNuMuBarCCcoh    ||
     m == kSKmNuMuBarCCoth ||
     m == kSKmNuMuBarNC1piC   ||
     m == kSKmNuMuBarNCoth_mu ||
     m == kSKmMEC             )
   {
    return NuMuSurvProbDecoh(sinsq2theta, musq, Et);
  }
  else 
  if(
     m == kSKmOscNuMuCCQE     ||
     m == kSKmOscNuMuCC1pi    ||
     m == kSKmOscNuMuCCcoh    ||
//     m == kSKmOscNuMuCCoth ||                                                                                                                                 
     m == kSKmOscNuMuCCoth  ||
     m == kSKmOscNuMuBarCCQE     ||
     m == kSKmOscNuMuBarCC1pi    ||
     m == kSKmOscNuMuBarCCcoh    ||
//     m == kSKmOscNuMuCCoth ||                                                                                                                                 
     m == kSKmOscNuMuBarCCoth  ||
     m == kSKmOscNuEBarCCQE     ||
     m == kSKmOscNuEBarCC1pi    ||
     m == kSKmOscNuEBarCCcoh    ||
//     m == kSKmOscNuMuCCoth ||                                                                                                                                 
     m == kSKmOscNuEBarCCoth
     )
  {
    return 0;
  }
  else {
    cerr << "NuMuSurvProbDecoh > Can't handle SK mode = " << m << endl;
    exit(1);
  }
}
//____________________________________________________________________________
t2kana::util::ArrIdx * t2kana::util::ArrIdx::fInstance = 0;
//____________________________________________________________________________
t2kana::util::ArrIdx::ArrIdx()
{
  this->Init();
  fInstance = 0;
}
//____________________________________________________________________________
t2kana::util::ArrIdx * t2kana::util::ArrIdx::Instance()
{
  if(fInstance == 0) {
    static t2kana::util::ArrIdx::Cleaner cleaner;
    cleaner.DummyMethodAndSilentCompiler();
    fInstance = new t2kana::util::ArrIdx;
  }
  return fInstance;
}
//____________________________________________________________________________
void t2kana::util::ArrIdx::Init(void)
{
  //
  // sample/mode arrays
  //
  for(int i = 0; i < 2*kns; i++) {
     fSuperkSampleIdx [i] = -1;
     for(int j = 0; j < 2*knm; j++) {
         fSuperkModeIdx [i][j] = -1;
     }
  }
  for(int is = 0; is < kns; is++) {
     SuperkSample_t s = kSKSamples[is];
     if(s == kSKsUndefined) continue;
     fSuperkSampleIdx [s] = is;
#ifdef __show_arr_idx_mesg_debug__
     cout << "Util/ArrIdx >> Sample: " 
          << util::AsString(s)
          << " is at slot: " << fSuperkSampleIdx [s] << endl;
#endif
     for(int im = 0; im < knm; im++) {
         SuperkMode_t m = kSKReactionModes[is][im];
         if(m == kSKmUndefined) continue;
         fSuperkModeIdx [s][m] = im;
#ifdef __show_arr_idx_mesg_debug__
         cout << "Util/ArrIdx >>   Mode: " 
              << util::AsString(s) << "/" << util::AsString(m)
              << " is at slot: " << fSuperkModeIdx [s][m] << endl;
#endif
     }
  } 
  // 
  // phys param arrays
  // 
  //XSEC
  for(int i = 0; i < knppmax; i++) {
    fXSecSystIdx[i] = -1;
  }
  for(int ip = 0; ip < kNumOfXSecSyst; ip++) {
    PhysParamId_t p = kListOfXSecSyst[ip];
    if(p == kPhypUndefined) continue;
    fXSecSystIdx[p] = ip;
#ifdef __show_arr_idx_mesg_debug__
         cout << "Util/ArrIdx >> FHC Systematic: " 
              << util::AsString(p) 
              << " is at slot: " << fXSecSystIdx[p] << endl;
#endif
  }
  for(int i = 0; i < knppmax; i++) {
    fRHCXSecSystIdx[i] = -1;
  }
  for(int ip = 0; ip < kNumOfRHCXSecSyst; ip++) {
    PhysParamId_t p = kListOfRHCXSecSyst[ip];
    if(p == kPhypUndefined) continue;
    fRHCXSecSystIdx[p] = ip;
#ifdef __show_arr_idx_mesg_debug__
         cout << "Util/ArrIdx >> RHC Systematic: " 
              << util::AsString(p) 
              << " is at slot: " << fRHCXSecSystIdx[p] << endl;
#endif
  }

  //BANFF
  for(int i = 0; i < knppmax; i++) {
    fBANFFSystIdx[i] = -1;
  }
  for(int ip = 0; ip < kNumOfParamsInBANFFfit; ip++) {
    PhysParamId_t p = PhysParamId_t ((int) kNuisBANFF00 + ip);
    if(p == kPhypUndefined) continue;
    fBANFFSystIdx[p] = ip;
#ifdef __show_arr_idx_mesg_debug__
     cout << "Util/ArrIdx >> BANFF: " 
          << util::AsString(p)
          << " is at slot: " << fBANFFSystIdx [p] << endl;
#endif
  }
  for(int i = 0; i < knppmax; i++) {
    fRHCBANFFSystIdx[i] = -1;
  }
  for(int ip = 0; ip < kNumOfRHCParamsInBANFFfit; ip++) {
    PhysParamId_t p = PhysParamId_t ((int) kRHCNuisBANFF00 + ip);
    if(p == kPhypUndefined) continue;
    fRHCBANFFSystIdx[p] = ip;
#ifdef __show_arr_idx_mesg_debug__
    cout << "Util/ArrIdx >> RHC BANFF: "
	 << util::AsString(p)
	 << " is at slot: " << fRHCBANFFSystIdx [p] << endl;
#endif
  }

#ifndef SKFSITOGETHER
  //FSI
  for(int i = 0; i < knppmax; i++) {
    fFSISystIdx[i] = -1;
  }
  for(int ip = 0; ip < kNumOfFSIParams; ip++) {
#if FSI_BINS != 3
    PhysParamId_t p = PhysParamId_t ((int) kNuisFSI000 + ip);
#else
    PhysParamId_t p = PhysParamId_t ((int) kNuisFSI006 + ip);
#endif
    if(p == kPhypUndefined) continue;
    fFSISystIdx[p] = ip;
#ifdef __show_arr_idx_mesg_debug__
     cout << "Util/ArrIdx >> FSI: " 
          << util::AsString(p)
          << " is at slot: " << fFSISystIdx [p] << endl;
#endif
  }
#endif

  //SKDET + FSI
#ifdef SKFSITOGETHER
  for(int i = 0; i < knppmax; i++) {
    fSKDetFSISystIdx[i] = -1;
  }
  for(int ip = 0; ip < kNumOfSKDetFSIParams; ip++) {
#if FSI_BINS != 3
    PhysParamId_t p = PhysParamId_t ((int) kNuisSKDetFSI000 + ip);
#else
    PhysParamId_t p = PhysParamId_t ((int) kNuisSKDetFSI006 + ip);
#endif
    if(p == kPhypUndefined) continue;
    fSKDetFSISystIdx[p] = ip;
#ifdef __show_arr_idx_mesg_debug__
     cout << "Util/ArrIdx >> SKDet+FSI: " 
          << util::AsString(p)
          << " is at slot: " << fSKDetFSISystIdx [p] << endl;
#endif
  }
#endif
#ifdef SKFSITOGETHER
  for(int i = 0; i < knppmax; i++) {
    fRHCSKDetFSISystIdx[i] = -1;
  }
  for(int ip = 0; ip < kNumOfRHCSKDetFSIParams; ip++) {
#if FSI_BINS != 3
    PhysParamId_t p = PhysParamId_t ((int) kRHCNuisSKDetFSI000 + ip);
#else
    PhysParamId_t p = PhysParamId_t ((int) kRHCNuisSKDetFSI006 + ip);
#endif
    if(p == kPhypUndefined) continue;
    fRHCSKDetFSISystIdx[p] = ip;
#ifdef __show_arr_idx_mesg_debug__
    cout << "Util/ArrIdx >> RHC SKDet+FSI: "
	 << util::AsString(p)
	 << " is at slot: " << fRHCSKDetFSISystIdx [p] << endl;
#endif
  }
#endif

#ifndef SKFSITOGETHER
#if FSI_BINS != 2
  //1Re SK eff
  for(int i = 0; i < knppmax; i++) {
    fSKDet1ReSystIdx[i] = -1;
  }
  for(int ip = 0; ip < kNumOfSKDetParams1Re; ip++) {
    PhysParamId_t p = PhysParamId_t ((int) kNuisSuperkEff1ReNueSigCC0 + ip);
    if(p == kPhypUndefined) continue;
    fSKDet1ReSystIdx[p] = ip;
#ifdef __show_arr_idx_mesg_debug__
     cout << "Util/ArrIdx >> 1Re SK Det: " 
          << util::AsString(p)
          << " is at slot: " << fSKDet1ReSystIdx [p] << endl;
#endif
  }
#endif
#endif

  //PhysParam
  for(int i = 0; i < knppmax; i++) {
    fPhysParamIdx[i] = -1;
  }
  for(int ip = 0; ip < knppmax; ip++) {
    PhysParamId_t p = PhysParamId_t(ip);
    if(p == kPhypUndefined) continue;
    fPhysParamIdx[p] = ip;
#ifdef __show_arr_idx_mesg_debug__
     cout << "Util/ArrIdx >> PhysParam: " 
          << util::AsString(p)
          << " is at slot: " << fPhysParamIdx [p] << endl;
#endif
  }
}
//____________________________________________________________________________
int t2kana::util::DiagonaliseJacobi(std::complex<double> A[3][3],  std::complex<double> Q[3][3], double w[3])
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
int t2kana::util::DiagonaliseGSL(int ndim, double * matrix_in, double * eigenvectors_out, double * eigenvalues_out)
{
#ifdef __GSL_ENABLED__ 
  gsl_matrix_complex * matr = gsl_matrix_complex_alloc(ndim,ndim);
  memcpy( (double*)gsl_matrix_complex_const_ptr(matr,0,0), matrix_in, ndim*ndim*2*sizeof(double) );
  
  gsl_matrix_complex * eigenvectors = gsl_matrix_complex_alloc(ndim,ndim);
  gsl_vector * eigenvalues  = gsl_vector_alloc(ndim);
  gsl_eigen_hermv_workspace * workspace = gsl_eigen_hermv_alloc(ndim);
  
  int status=gsl_eigen_hermv(matr, eigenvalues, eigenvectors, workspace);
  
  memcpy( eigenvectors_out, (double*)gsl_matrix_complex_const_ptr(eigenvectors,0,0) , ndim*ndim*2*sizeof(double));
  memcpy( eigenvalues_out , (double*)gsl_vector_const_ptr        (eigenvalues ,0)   , ndim*sizeof(double));
                           
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
