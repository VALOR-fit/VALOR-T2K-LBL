//____________________________________________________________________________
/*

 >> ChangeLog:

 For changes up to version 2010_v3.0.0 (last implementation of the 2010-era analysis)
 see code in the corresponding SVN branch.

 Changes since 2011_v0.9.0:

 @
   
*/
//____________________________________________________________________________

//#define __show_cholesky_decomposition_mesg_debug__
//#define __show_cholesky_gen_corr_param_variations_mesg_debug__
//#define __show_xml_mesg_info__
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

TRandom3 gR(0);

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
  if(nmc <= 0. || ndata <= 0.) return 2.*(nmc-ndata);
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
    g(k) = gR.Gaus();
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
    g(k) = gR.Gaus();
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
  TRandom3 rand(0);

  double mix = 0;

  NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
  switch(omp) {
   case(kNuMixSinsq2Theta):
   {
     double mean  = conv::kGlobFitSinsq2Theta12; 
     double sigma = conv::kGlobFitSinsq2Theta12 * conv::kGlobFitErrSinsq2Theta12;
     double sinsq2theta = rand.Gaus(mean,sigma);
     mix = sinsq2theta;
     break;
   }
   case(kNuMixSinsqTheta):
   {
     double mean  = conv::kGlobFitSinsqTheta12; 
     double sigma = conv::kGlobFitSinsqTheta12 * conv::kGlobFitErrSinsqTheta12;
     double sinsqtheta = rand.Gaus(mean,sigma);
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
  TRandom3 rand(0);

  double mix = 0;

  NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
  switch(omp) {
   case(kNuMixSinsq2Theta):
   {
     double mean  = conv::kGlobFitSinsq2Theta23; 
     double sigma = conv::kGlobFitSinsq2Theta23 * conv::kGlobFitErrSinsq2Theta23;
     double sinsq2theta = rand.Gaus(mean,sigma);
     mix = sinsq2theta;
     break;
   }
   case(kNuMixSinsqTheta):
   {
     double mean  = conv::kGlobFitSinsqTheta23; 
     double sigma = conv::kGlobFitSinsqTheta23 * conv::kGlobFitErrSinsqTheta23;
     double sinsqtheta = rand.Gaus(mean,sigma);
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
  TRandom3 rand(0);

  double mix = 0;

  NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
  switch(omp) {
  case(kNuMixSinsq2Theta):
    {
      double mean  = conv::kGlobFitSinsq2Theta13;
      double sigma = conv::kGlobFitSinsq2Theta13 * conv::kGlobFitErrSinsq2Theta13;
      double sinsq2theta = rand.Gaus(mean,sigma);
      mix = sinsq2theta;
      break;
    }
  case(kNuMixSinsqTheta):
    {
      double mean  = conv::kGlobFitSinsqTheta13;
      double sigma = conv::kGlobFitSinsqTheta13 * conv::kGlobFitErrSinsqTheta13;
      double sinsqtheta = rand.Gaus(mean,sigma);
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
  TRandom3 rand(0);
  double min = -conv::kpi;
  double max =  conv::kpi;
  return rand.Uniform(min,max);
}
//____________________________________________________________________________
double t2kana::util::RandomizeDm12sq(void)
{
  TRandom3 rand(0);
  double mean  = conv::kGlobFitDm12sq; 
  double sigma = conv::kGlobFitDm12sq * conv::kGlobFitErrDm12sq;
  return rand.Gaus(mean,sigma);
}
//____________________________________________________________________________
double t2kana::util::RandomizeDm23sq(void)
{
  TRandom3 rand(0);
  double mean  = conv::kGlobFitDm23sq; 
  double sigma = conv::kGlobFitDm23sq * conv::kGlobFitErrDm23sq;
  return rand.Gaus(mean,sigma);
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
  string mixparam = "?";
  switch(omp) {
      case(kNuMixSinsq2Theta):
         mixparam = "sinsq2theta";
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
       out << "dmsq[";
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
    case ( kPhypDm23sq                ) : return "|Dm^2_32|";                  break;
    case ( kPhypDm12sq                ) : return "Dm^2_12";                    break;

    case ( kPhypAlphaDecay            ) : return "alpha_decay";                break;
    case ( kPhypSinsq2ThetaDecay      ) : return "sin^2(2*theta_decay)";       break;
    case ( kPhypMusqDecoh             ) : return "mu^2_decoh";                 break;
    case ( kPhypSinsq2ThetaDecoh      ) : return "sin^2(2*theta_decoh)";       break;

    /*
     * 2011-analysis nuisance parameters
     */

    //
    // SK detector
    //
    case ( kNuisSuperkEnergyScale     ) : return "SK energy scale";                 break;    
    case ( kNuisSuperkEff1RmuCCQEoth  ) : return "SK efficiency; 1Rmu; CCQE/Other"; break;
    case ( kNuisSuperkEff1RmuCCQErc0u ) : return "SK efficiency; 1Rmu; CCQE/RC0";   break;
    case ( kNuisSuperkEff1RmuCCQErc1u ) : return "SK efficiency; 1Rmu; CCQE/RC1";   break;
    case ( kNuisSuperkEff1RmuCCQErc2u ) : return "SK efficiency; 1Rmu; CCQE/RC2";   break;
    case ( kNuisSuperkEff1RmuCCnonQE  ) : return "SK efficiency; 1Rmu; CCnonQE";    break;
    case ( kNuisSuperkEff1RmuNC       ) : return "SK efficiency; 1Rmu; NC";         break;
    case ( kNuisSuperkEff1RmuNueCC    ) : return "SK efficiency; 1Rmu; nue CC";     break;

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
#if BANFF_BINS <= 2
    case ( kNuisBANFF16 ) :
    case ( kNuisBANFF17 ) :
    case ( kNuisBANFF18 ) :
    case ( kNuisBANFF19 ) : 
    case ( kNuisBANFF20 ) : 
    case ( kNuisBANFF21 ) : 
    case ( kNuisBANFF22 ) :
    case ( kNuisBANFF23 ) :
    case ( kNuisBANFF24 ) :
    case ( kNuisBANFF25 ) :
    case ( kNuisBANFF26 ) :
    case ( kNuisBANFF27 ) :
    case ( kNuisBANFF28 ) :
    case ( kNuisBANFF29 ) :
    case ( kNuisBANFF30 ) :
#if BANFF_BINS == -1
    case ( kNuisBANFF31 ) :
    case ( kNuisBANFF32 ) :
    case ( kNuisBANFF33 ) :
    case ( kNuisBANFF34 ) :
    case ( kNuisBANFF35 ) :
    case ( kNuisBANFF36 ) :
    case ( kNuisBANFF37 ) :
    case ( kNuisBANFF38 ) :
    case ( kNuisBANFF39 ) :
    case ( kNuisBANFF40 ) :
    case ( kNuisBANFF41 ) :
    case ( kNuisBANFF42 ) :
    case ( kNuisBANFF43 ) :
    case ( kNuisBANFF44 ) :
    case ( kNuisBANFF45 ) :
    case ( kNuisBANFF46 ) :
    case ( kNuisBANFF47 ) :
    case ( kNuisBANFF48 ) :
#endif
#endif
    {
      int ip = (int)p - (int)kNuisBANFF00;
      return Form("BANFF %2d; SK %s flux; %.1f - %.1f GeV",
                   ip, util::Pdg2String(kNuPdgCodeOfFluxParamInBANFFfit[ip]).c_str(), 
                   kEtrueRangeOfFluxParamInBANFFfit[ip][0],
                   kEtrueRangeOfFluxParamInBANFFfit[ip][1]);     
      break;
    }
#if BANFF_BINS == 4
    case ( kNuisBANFFMAQE    ) : return "BANFF; 1 / (Ma QE)^2";                   break;
#else
    case ( kNuisBANFFMAQE    ) : return "BANFF; Ma QE";                           break;
#endif
    case ( kNuisBANFFMARES   ) : return "BANFF; Ma RES";                          break;
    case ( kNuisBANFFCCQEE0  ) : return "BANFF; Norm; CC QE; 0.0 - 1.5 GeV";      break;
#if BANFF_BINS >= 2
    case ( kNuisBANFFCCQEE1  ) : return "BANFF; Norm; CC QE; 1.5 - 3.5 GeV";      break;
    case ( kNuisBANFFCCQEE2  ) : return "BANFF; Norm; CC QE; 3.5 - 30.0 GeV";     break;
#endif
    case ( kNuisBANFFCC1piE0 ) : return "BANFF; Norm; CC 1pi; 0.0 - 2.5 GeV";     break;
#if BANFF_BINS >= 2
    case ( kNuisBANFFCC1piE1 ) : return "BANFF; Norm; CC 1pi; 2.5 - 30.0 GeV";    break;
#endif
#if BANFF_BINS <= 2
    case ( kNuisBANFFNC1pi0  ) : return "BANFF; Norm; NC 1pi0";                   break;
#endif
#if BANFF_BINS == 5
    case ( kNuisBANFFMEC     ) : return "BANFF; MEC";                             break;
#endif

    //
    // XSec params not correlated with flux
    //

    case ( kNuisXSecFermiMom   ) : return "Fermi Momentum";                  break;
    case ( kNuisXSecSpecFunc   ) : return "Spectral Func";                   break; 
    case ( kNuisXSecWshape     ) : return "W-shape";                         break;
    case ( kNuisXSecShapeCCoth ) : return "Shape; CC Oth";                   break;
    case ( kNuisXSecPiLessDelta) : return "Pi-less delta decay";             break;
    case ( kNuisXSecENuShape   ) : return "Shape; CC1pi ENu";                break;
    case ( kNuisXSecBindE      ) : return "Binding energy";                  break;

#if BANFF_BINS <= 1
    case ( kNuisXSecNormCCQEE1 ) : return "Norm; CC QE Norm; 1.5 - 3.5 GeV"; break; 
    case ( kNuisXSecNormCCQEE2 ) : return "Norm; CC QE Norm; > 3.5 GeV";     break; 
    case ( kNuisXSecNormCC1piE1) : return "Norm; CC 1pi; > 2.5 GeV";         break; 
#endif
    case ( kNuisXSecNormCCcoh  ) : return "Norm; CC Coh";                    break;
    case ( kNuisXSecNormNC1piC ) : return "Norm; NC 1piC";                   break;
    case ( kNuisXSecNormNCcoh  ) : return "Norm; NC Coh";                    break;
    case ( kNuisXSecNormNCoth  ) : return "Norm; NC Oth";                    break;
    case ( kNuisXSecNormNuEToNuMu  ) : return "Norm; NuE To NuMu";           break;
    case ( kNuisXSecNormAntiNuToNu ) : return "Norm; AntiNu To Nu";          break;

    case ( kNuisFSI000 ) :
    case ( kNuisFSI001 ) :
    case ( kNuisFSI002 ) :
    case ( kNuisFSI003 ) :
    case ( kNuisFSI004 ) :
    case ( kNuisFSI005 ) :
#if FSI_BINS >= 2
    case ( kNuisFSI006 ) :
    case ( kNuisFSI007 ) :
    case ( kNuisFSI008 ) :
    case ( kNuisFSI009 ) :
    case ( kNuisFSI010 ) :
    case ( kNuisFSI011 ) :
    case ( kNuisFSI012 ) :
    case ( kNuisFSI013 ) :
    case ( kNuisFSI014 ) :
#if FSI_BINS >= 3
    case ( kNuisFSI015 ) :
    case ( kNuisFSI016 ) :
    case ( kNuisFSI017 ) :
    case ( kNuisFSI018 ) :
    case ( kNuisFSI019 ) :
    case ( kNuisFSI020 ) :
    case ( kNuisFSI021 ) :
    case ( kNuisFSI022 ) :
    case ( kNuisFSI023 ) :
    case ( kNuisFSI024 ) :
    case ( kNuisFSI025 ) :
    case ( kNuisFSI026 ) :
    case ( kNuisFSI027 ) :
    case ( kNuisFSI028 ) :
    case ( kNuisFSI029 ) :
    case ( kNuisFSI030 ) :
    case ( kNuisFSI031 ) :
    case ( kNuisFSI032 ) :
    case ( kNuisFSI033 ) :
    case ( kNuisFSI034 ) :
    case ( kNuisFSI035 ) :
    case ( kNuisFSI036 ) :
    case ( kNuisFSI037 ) :
    case ( kNuisFSI038 ) :
    case ( kNuisFSI039 ) :
    case ( kNuisFSI040 ) :
    case ( kNuisFSI041 ) :
    case ( kNuisFSI042 ) :
    case ( kNuisFSI043 ) :
    case ( kNuisFSI044 ) :
    case ( kNuisFSI045 ) :
    case ( kNuisFSI046 ) :
    case ( kNuisFSI047 ) :
    case ( kNuisFSI048 ) :
#endif
#endif
    {
      int ip = (int)p - (int)kNuisFSI000;
      string nutypes = "";
#if FSI_BINS == 2
      switch(p) {
        case (kNuisFSI000) :
        case (kNuisFSI001) :
        case (kNuisFSI002) :
	  nutypes = "numu"; break;
        case (kNuisFSI003) :
        case (kNuisFSI004) :
        case (kNuisFSI005) :
	  nutypes = "numubar"; break;
        case (kNuisFSI006) :
        case (kNuisFSI007) :
        case (kNuisFSI008) :
	  nutypes = "nue"; break;
        case (kNuisFSI009) :
        case (kNuisFSI010) :
        case (kNuisFSI011) :
	  nutypes = "nuebar"; break;
        case (kNuisFSI012) :
        case (kNuisFSI013) :
        case (kNuisFSI014) :
	  nutypes = "oscillated nue"; break;
        default : nutypes = "--";
      }
#elif FSI_BINS == 1
      switch(p) {
        case (kNuisFSI000) :
        case (kNuisFSI001) :
        case (kNuisFSI002) :
	  nutypes = "numu/numubar CCQE"; break;
        case (kNuisFSI003) :
	  nutypes = "numu/numubar CCnQE"; break;
        case (kNuisFSI004) :
	  nutypes = "nue/nuebar/signue CC"; break;
        case (kNuisFSI005) :
	  nutypes = "all NC"; break;
        default : nutypes = "--";
      }
#elif FSI_BINS == 3
      if(ip == (int)kNuisFSI000)
	nutypes = "numu/numubar CCQE";
      else if(ip <= (int)kNuisFSI016)
	nutypes = "numu/numubar CCnQE";
      else if(ip <= (int)kNuisFSI032)
	nutypes = "nue/nuebar/signue CC";
      else if(ip <= (int)kNuisFSI048)
	nutypes = "all NC";
#endif
      return Form("FSI/SI %2d; Ereco range %.1f - %.1f GeV; %s",ip,
		  kFSIParamsErecoCut[ip][0],kFSIParamsErecoCut[ip][1],nutypes.c_str());     
      break;
    }//FSI

    case ( kNuisUnBANFF00 ) :
    case ( kNuisUnBANFF01 ) :
    case ( kNuisUnBANFF02 ) :
    case ( kNuisUnBANFF03 ) :
    case ( kNuisUnBANFF04 ) :
    case ( kNuisUnBANFF05 ) :
    case ( kNuisUnBANFF06 ) :
    case ( kNuisUnBANFF07 ) :
    case ( kNuisUnBANFF08 ) :
    case ( kNuisUnBANFF09 ) :
    case ( kNuisUnBANFF10 ) :
    case ( kNuisUnBANFF11 ) :
    case ( kNuisUnBANFF12 ) :
    case ( kNuisUnBANFF13 ) :
    case ( kNuisUnBANFF14 ) :
    case ( kNuisUnBANFF15 ) :
    case ( kNuisUnBANFF16 ) :
    case ( kNuisUnBANFF17 ) :
    case ( kNuisUnBANFF18 ) :
    case ( kNuisUnBANFF19 ) : 
    case ( kNuisUnBANFF20 ) : 
    case ( kNuisUnBANFF21 ) : 
    case ( kNuisUnBANFF22 ) :
#if BANFF_BINS <= 2
    case ( kNuisUnBANFF23 ) :
    case ( kNuisUnBANFF24 ) :
    case ( kNuisUnBANFF25 ) :
    case ( kNuisUnBANFF26 ) :
    case ( kNuisUnBANFF27 ) :
    case ( kNuisUnBANFF28 ) :
    case ( kNuisUnBANFF29 ) :
    case ( kNuisUnBANFF30 ) :
    case ( kNuisUnBANFF31 ) :
    case ( kNuisUnBANFF32 ) :
    case ( kNuisUnBANFF33 ) :
    case ( kNuisUnBANFF34 ) :
    case ( kNuisUnBANFF35 ) :
#if BANFF_BINS != 1
    case ( kNuisUnBANFF36 ) :
    case ( kNuisUnBANFF37 ) :
    case ( kNuisUnBANFF38 ) :
#if BANFF_BINS == -1
    case ( kNuisUnBANFF39 ) :
    case ( kNuisUnBANFF40 ) :
    case ( kNuisUnBANFF41 ) :
    case ( kNuisUnBANFF42 ) :
    case ( kNuisUnBANFF43 ) :
    case ( kNuisUnBANFF44 ) :
    case ( kNuisUnBANFF45 ) :
    case ( kNuisUnBANFF46 ) :
    case ( kNuisUnBANFF47 ) :
    case ( kNuisUnBANFF48 ) :
    case ( kNuisUnBANFF49 ) :
    case ( kNuisUnBANFF50 ) :
    case ( kNuisUnBANFF51 ) :
    case ( kNuisUnBANFF52 ) :
    case ( kNuisUnBANFF53 ) :
#endif
#endif
#endif
    {
      int ip = (int)p - (int)kNuisUnBANFF00;
      return Form("Uncorrelated BANFF %2d", ip);
      break;
    }

    case ( kNuisEDBANFF00 ) :
    case ( kNuisEDBANFF01 ) :
    case ( kNuisEDBANFF02 ) :
    case ( kNuisEDBANFF03 ) :
    case ( kNuisEDBANFF04 ) :
    case ( kNuisEDBANFF05 ) :
    case ( kNuisEDBANFF06 ) :
    case ( kNuisEDBANFF07 ) :
    case ( kNuisEDBANFF08 ) :
    case ( kNuisEDBANFF09 ) :
    case ( kNuisEDBANFF10 ) :
    case ( kNuisEDBANFF11 ) :
    case ( kNuisEDBANFF12 ) :
    case ( kNuisEDBANFF13 ) :
    case ( kNuisEDBANFF14 ) :
    case ( kNuisEDBANFF15 ) :
    case ( kNuisEDBANFF16 ) :
    case ( kNuisEDBANFF17 ) :
    case ( kNuisEDBANFF18 ) :
    case ( kNuisEDBANFF19 ) : 
    case ( kNuisEDBANFF20 ) : 
    case ( kNuisEDBANFF21 ) : 
    case ( kNuisEDBANFF22 ) :
#if BANFF_BINS <= 2
    case ( kNuisEDBANFF23 ) :
    case ( kNuisEDBANFF24 ) :
    case ( kNuisEDBANFF25 ) :
    case ( kNuisEDBANFF26 ) :
    case ( kNuisEDBANFF27 ) :
    case ( kNuisEDBANFF28 ) :
    case ( kNuisEDBANFF29 ) :
    case ( kNuisEDBANFF30 ) :
    case ( kNuisEDBANFF31 ) :
    case ( kNuisEDBANFF32 ) :
    case ( kNuisEDBANFF33 ) :
    case ( kNuisEDBANFF34 ) :
    case ( kNuisEDBANFF35 ) :
#if BANFF_BINS != 1
    case ( kNuisEDBANFF36 ) :
    case ( kNuisEDBANFF37 ) :
    case ( kNuisEDBANFF38 ) :
#if BANFF_BINS == -1
    case ( kNuisEDBANFF39 ) :
    case ( kNuisEDBANFF40 ) :
    case ( kNuisEDBANFF41 ) :
    case ( kNuisEDBANFF42 ) :
    case ( kNuisEDBANFF43 ) :
    case ( kNuisEDBANFF44 ) :
    case ( kNuisEDBANFF45 ) :
    case ( kNuisEDBANFF46 ) :
    case ( kNuisEDBANFF47 ) :
    case ( kNuisEDBANFF48 ) :
    case ( kNuisEDBANFF49 ) :
    case ( kNuisEDBANFF50 ) :
    case ( kNuisEDBANFF51 ) :
    case ( kNuisEDBANFF52 ) :
    case ( kNuisEDBANFF53 ) :
#endif
#endif
#endif
    {
      int ip = (int)p - (int)kNuisEDBANFF00;
      return Form("EigenDecomposition BANFF %2d", ip);
      break;
    }

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

    case ( kPhypDeltaCP               ) : return "$\\delta_{cp}$";                   break;
    case ( kPhypDm23sq                ) : return "$|\\Deltam^2_{32}|$";                  break;
    case ( kPhypDm12sq                ) : return "$\\Deltam^2_{12}$";                    break;

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
    case ( kNuisSuperkEff1RmuCCQEoth  ) : return "SK efficiency, 1Rmu, CCQE//Other"; break;
    case ( kNuisSuperkEff1RmuCCQErc0u ) : return "SK efficiency, 1Rmu, CCQE//RC0";   break;
    case ( kNuisSuperkEff1RmuCCQErc1u ) : return "SK efficiency, 1Rmu, CCQE//RC1";   break;
    case ( kNuisSuperkEff1RmuCCQErc2u ) : return "SK efficiency, 1Rmu, CCQE//RC2";   break;
    case ( kNuisSuperkEff1RmuCCnonQE  ) : return "SK efficiency, 1Rmu, CCnonQE";    break;
    case ( kNuisSuperkEff1RmuNC       ) : return "SK efficiency, 1Rmu, NC";         break;
    case ( kNuisSuperkEff1RmuNueCC    ) : return "SK efficiency, 1Rmu, nue CC";     break;

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
#if BANFF_BINS <= 2
    case ( kNuisBANFF16 ) :
    case ( kNuisBANFF17 ) :
    case ( kNuisBANFF18 ) :
    case ( kNuisBANFF19 ) : 
    case ( kNuisBANFF20 ) : 
    case ( kNuisBANFF21 ) : 
    case ( kNuisBANFF22 ) :
    case ( kNuisBANFF23 ) :
    case ( kNuisBANFF24 ) :
    case ( kNuisBANFF25 ) :
    case ( kNuisBANFF26 ) :
    case ( kNuisBANFF27 ) :
    case ( kNuisBANFF28 ) :
    case ( kNuisBANFF29 ) :
    case ( kNuisBANFF30 ) :
#if BANFF_BINS == -1
    case ( kNuisBANFF31 ) :
    case ( kNuisBANFF32 ) :
    case ( kNuisBANFF33 ) :
    case ( kNuisBANFF34 ) :
    case ( kNuisBANFF35 ) :
    case ( kNuisBANFF36 ) :
    case ( kNuisBANFF37 ) :
    case ( kNuisBANFF38 ) :
    case ( kNuisBANFF39 ) :
    case ( kNuisBANFF40 ) :
    case ( kNuisBANFF41 ) :
    case ( kNuisBANFF42 ) :
    case ( kNuisBANFF43 ) :
    case ( kNuisBANFF44 ) :
    case ( kNuisBANFF45 ) :
    case ( kNuisBANFF46 ) :
    case ( kNuisBANFF47 ) :
    case ( kNuisBANFF48 ) :
#endif
#endif
    {
      int ip = (int)p - (int)kNuisBANFF00;
      return Form("BANFF %2d, SK %s flux, %.1f - %.1f GeV",
                   ip, util::Pdg2String(kNuPdgCodeOfFluxParamInBANFFfit[ip]).c_str(), 
                   kEtrueRangeOfFluxParamInBANFFfit[ip][0],
                   kEtrueRangeOfFluxParamInBANFFfit[ip][1]);     
      break;
    }
#if BANFF_BINS == 4
    case ( kNuisBANFFMAQE    ) : return "BANFF, 1 / (Ma QE)$^2$";                 break;
#else
    case ( kNuisBANFFMAQE    ) : return "BANFF, Ma QE";                           break;
#endif
    case ( kNuisBANFFMARES   ) : return "BANFF, Ma RES";                          break;
    case ( kNuisBANFFCCQEE0  ) : return "BANFF; Norm; CC QE; 0.0 - 1.5 GeV";      break;
#if BANFF_BINS >= 2
    case ( kNuisBANFFCCQEE1  ) : return "BANFF; Norm; CC QE; 1.5 - 3.5 GeV";      break;
    case ( kNuisBANFFCCQEE2  ) : return "BANFF; Norm; CC QE; 3.5 - 30.0 GeV";     break;
#endif
    case ( kNuisBANFFCC1piE0 ) : return "BANFF; Norm; CC 1pi; 0.0 - 2.5 GeV";     break;
#if BANFF_BINS >= 2
    case ( kNuisBANFFCC1piE1 ) : return "BANFF; Norm; CC 1pi; 2.5 - 30.0 GeV";    break;
#endif
#if BANFF_BINS <= 2
    case ( kNuisBANFFNC1pi0  ) : return "BANFF, Norm, NC 1pi0";                   break;
#endif
#if BANFF_BINS == 5
    case ( kNuisBANFFMEC     ) : return "BANFF; MEC";                             break;
#endif

    //
    // XSec params not correlated with flux
    //

    case ( kNuisXSecFermiMom   ) : return "Fermi Momentum";                  break;
    case ( kNuisXSecSpecFunc   ) : return "Spectral Func";                   break; 
    case ( kNuisXSecWshape     ) : return "W-shape";                         break;
    case ( kNuisXSecShapeCCoth ) : return "Shape; CC Oth";                   break;
    case ( kNuisXSecPiLessDelta) : return "Pi-less delta decay";             break;
    case ( kNuisXSecENuShape   ) : return "Shape; CC1pi ENu";                break;
    case ( kNuisXSecBindE      ) : return "Binding energy";                  break;

#if BANFF_BINS <= 1
    case ( kNuisXSecNormCCQEE1 ) : return "Norm, CC QE Norm, 1.5 - 3.5 GeV"; break; 
    case ( kNuisXSecNormCCQEE2 ) : return "Norm, CC QE Norm, $>$ 3.5 GeV";   break; 
    case ( kNuisXSecNormCC1piE1) : return "Norm, CC 1pi, $>$ 2.5 GeV";       break; 
#endif
    case ( kNuisXSecNormCCcoh  ) : return "Norm, CC Coh";                    break;
    case ( kNuisXSecNormNCcoh  ) : return "Norm, NC Coh";                    break;
    case ( kNuisXSecNormNC1piC ) : return "Norm, NC 1piC";                   break;
    case ( kNuisXSecNormNCoth  ) : return "Norm, NC Oth";                    break;
    case ( kNuisXSecNormNuEToNuMu  ) : return "Norm, $\\nu_e$ To $\\nu_\\mu$";  break;
    case ( kNuisXSecNormAntiNuToNu ) : return "Norm, $\\bar{\\nu}$ To $\\nu$";  break;

    case ( kNuisFSI000 ) :
    case ( kNuisFSI001 ) :
    case ( kNuisFSI002 ) :
    case ( kNuisFSI003 ) :
    case ( kNuisFSI004 ) :
    case ( kNuisFSI005 ) :
#if FSI_BINS >= 2
    case ( kNuisFSI006 ) :
    case ( kNuisFSI007 ) :
    case ( kNuisFSI008 ) :
    case ( kNuisFSI009 ) :
    case ( kNuisFSI010 ) :
    case ( kNuisFSI011 ) :
    case ( kNuisFSI012 ) :
    case ( kNuisFSI013 ) :
    case ( kNuisFSI014 ) :
#if FSI_BINS >= 3
    case ( kNuisFSI015 ) :
    case ( kNuisFSI016 ) :
    case ( kNuisFSI017 ) :
    case ( kNuisFSI018 ) :
    case ( kNuisFSI019 ) :
    case ( kNuisFSI020 ) :
    case ( kNuisFSI021 ) :
    case ( kNuisFSI022 ) :
    case ( kNuisFSI023 ) :
    case ( kNuisFSI024 ) :
    case ( kNuisFSI025 ) :
    case ( kNuisFSI026 ) :
    case ( kNuisFSI027 ) :
    case ( kNuisFSI028 ) :
    case ( kNuisFSI029 ) :
    case ( kNuisFSI030 ) :
    case ( kNuisFSI031 ) :
    case ( kNuisFSI032 ) :
    case ( kNuisFSI033 ) :
    case ( kNuisFSI034 ) :
    case ( kNuisFSI035 ) :
    case ( kNuisFSI036 ) :
    case ( kNuisFSI037 ) :
    case ( kNuisFSI038 ) :
    case ( kNuisFSI039 ) :
    case ( kNuisFSI040 ) :
    case ( kNuisFSI041 ) :
    case ( kNuisFSI042 ) :
    case ( kNuisFSI043 ) :
    case ( kNuisFSI044 ) :
    case ( kNuisFSI045 ) :
    case ( kNuisFSI046 ) :
    case ( kNuisFSI047 ) :
    case ( kNuisFSI048 ) :
#endif
#endif
    {
      string nutypes = "";
      int ip = (int)p - (int)kNuisFSI000;
#if FSI_BINS == 2
      switch(p) {
        case (kNuisFSI000) :
        case (kNuisFSI001) :
        case (kNuisFSI002) :
	  nutypes = "$\\nu_\\mu$"; break;
        case (kNuisFSI003) :
        case (kNuisFSI004) :
        case (kNuisFSI005) :
	  nutypes = "$\\bar{\\nu_\\mu}$"; break;
        case (kNuisFSI006) :
        case (kNuisFSI007) :
        case (kNuisFSI008) :
	  nutypes = "$\\nu_e$"; break;
        case (kNuisFSI009) :
        case (kNuisFSI010) :
        case (kNuisFSI011) :
	  nutypes = "$\\bar{\\nu_e}$"; break;
        case (kNuisFSI012) :
        case (kNuisFSI013) :
        case (kNuisFSI014) :
	  nutypes = "Oscillated $\\nu_e$"; break;
        default : nutypes = "--";
      }
#elif FSI_BINS == 1
      switch(p) {
        case (kNuisFSI000) :
        case (kNuisFSI001) :
        case (kNuisFSI002) :
	  nutypes = "$\\nu_\\mu$/$\\bar{\\nu_\\mu}$ CCQE"; break;
        case (kNuisFSI003) :
	  nutypes = "$\\nu_\\mu$/$\\bar{\\nu_\\mu}$ CCnQE"; break;
        case (kNuisFSI004) :
	  nutypes = "$\\nu_e$/$\\bar{\\nu_e}$/sig$\\nu_e$ CC"; break;
        case (kNuisFSI005) :
	  nutypes = "all NC"; break;
        default : nutypes = "--";
      }
#elif FSI_BINS == 3
      if(ip == (int)kNuisFSI000)
	nutypes = "$\\nu_\\mu$/$\\bar{\\nu_\\mu}$ CCQE";
      else if(ip <= (int)kNuisFSI016)
	nutypes = "$\\nu_\\mu$/$\\bar{\\nu_\\mu}$ CCnQE";
      else if(ip <= (int)kNuisFSI032)
	nutypes = "$\\nu_e$/$\\bar{\\nu_e}$/sig$\\nu_e$ CC";
      else if(ip <= (int)kNuisFSI048)
	nutypes = "all NC";
#endif
      return Form("FSI/SI %2d; Ereco range %.1f - %.1f GeV; %s",ip,
		  kFSIParamsErecoCut[ip][0],kFSIParamsErecoCut[ip][1],nutypes.c_str());     
      break;
    }//FSI

    case ( kNuisUnBANFF00 ) :
    case ( kNuisUnBANFF01 ) :
    case ( kNuisUnBANFF02 ) :
    case ( kNuisUnBANFF03 ) :
    case ( kNuisUnBANFF04 ) :
    case ( kNuisUnBANFF05 ) :
    case ( kNuisUnBANFF06 ) :
    case ( kNuisUnBANFF07 ) :
    case ( kNuisUnBANFF08 ) :
    case ( kNuisUnBANFF09 ) :
    case ( kNuisUnBANFF10 ) :
    case ( kNuisUnBANFF11 ) :
    case ( kNuisUnBANFF12 ) :
    case ( kNuisUnBANFF13 ) :
    case ( kNuisUnBANFF14 ) :
    case ( kNuisUnBANFF15 ) :
    case ( kNuisUnBANFF16 ) :
    case ( kNuisUnBANFF17 ) :
    case ( kNuisUnBANFF18 ) :
    case ( kNuisUnBANFF19 ) : 
    case ( kNuisUnBANFF20 ) : 
    case ( kNuisUnBANFF21 ) : 
    case ( kNuisUnBANFF22 ) :
#if BANFF_BINS <= 2
    case ( kNuisUnBANFF23 ) :
    case ( kNuisUnBANFF24 ) :
    case ( kNuisUnBANFF25 ) :
    case ( kNuisUnBANFF26 ) :
    case ( kNuisUnBANFF27 ) :
    case ( kNuisUnBANFF28 ) :
    case ( kNuisUnBANFF29 ) :
    case ( kNuisUnBANFF30 ) :
    case ( kNuisUnBANFF31 ) :
    case ( kNuisUnBANFF32 ) :
    case ( kNuisUnBANFF33 ) :
    case ( kNuisUnBANFF34 ) :
    case ( kNuisUnBANFF35 ) :
#if BANFF_BINS != 1
    case ( kNuisUnBANFF36 ) :
    case ( kNuisUnBANFF37 ) :
    case ( kNuisUnBANFF38 ) :
#if BANFF_BINS == -1
    case ( kNuisUnBANFF39 ) :
    case ( kNuisUnBANFF40 ) :
    case ( kNuisUnBANFF41 ) :
    case ( kNuisUnBANFF42 ) :
    case ( kNuisUnBANFF43 ) :
    case ( kNuisUnBANFF44 ) :
    case ( kNuisUnBANFF45 ) :
    case ( kNuisUnBANFF46 ) :
    case ( kNuisUnBANFF47 ) :
    case ( kNuisUnBANFF48 ) :
    case ( kNuisUnBANFF49 ) :
    case ( kNuisUnBANFF50 ) :
    case ( kNuisUnBANFF51 ) :
    case ( kNuisUnBANFF52 ) :
    case ( kNuisUnBANFF53 ) :
#endif
#endif
#endif
    {
      int ip = (int)p - (int)kNuisUnBANFF00;
      return Form("Uncorrelated BANFF %2d", ip);
      break;
    }

    case ( kNuisEDBANFF00 ) :
    case ( kNuisEDBANFF01 ) :
    case ( kNuisEDBANFF02 ) :
    case ( kNuisEDBANFF03 ) :
    case ( kNuisEDBANFF04 ) :
    case ( kNuisEDBANFF05 ) :
    case ( kNuisEDBANFF06 ) :
    case ( kNuisEDBANFF07 ) :
    case ( kNuisEDBANFF08 ) :
    case ( kNuisEDBANFF09 ) :
    case ( kNuisEDBANFF10 ) :
    case ( kNuisEDBANFF11 ) :
    case ( kNuisEDBANFF12 ) :
    case ( kNuisEDBANFF13 ) :
    case ( kNuisEDBANFF14 ) :
    case ( kNuisEDBANFF15 ) :
    case ( kNuisEDBANFF16 ) :
    case ( kNuisEDBANFF17 ) :
    case ( kNuisEDBANFF18 ) :
    case ( kNuisEDBANFF19 ) : 
    case ( kNuisEDBANFF20 ) : 
    case ( kNuisEDBANFF21 ) : 
    case ( kNuisEDBANFF22 ) :
#if BANFF_BINS <= 2
    case ( kNuisEDBANFF23 ) :
    case ( kNuisEDBANFF24 ) :
    case ( kNuisEDBANFF25 ) :
    case ( kNuisEDBANFF26 ) :
    case ( kNuisEDBANFF27 ) :
    case ( kNuisEDBANFF28 ) :
    case ( kNuisEDBANFF29 ) :
    case ( kNuisEDBANFF30 ) :
    case ( kNuisEDBANFF31 ) :
    case ( kNuisEDBANFF32 ) :
    case ( kNuisEDBANFF33 ) :
    case ( kNuisEDBANFF34 ) :
    case ( kNuisEDBANFF35 ) :
#if BANFF_BINS != 1
    case ( kNuisEDBANFF36 ) :
    case ( kNuisEDBANFF37 ) :
    case ( kNuisEDBANFF38 ) :
#if BANFF_BINS == -1
    case ( kNuisEDBANFF39 ) :
    case ( kNuisEDBANFF40 ) :
    case ( kNuisEDBANFF41 ) :
    case ( kNuisEDBANFF42 ) :
    case ( kNuisEDBANFF43 ) :
    case ( kNuisEDBANFF44 ) :
    case ( kNuisEDBANFF45 ) :
    case ( kNuisEDBANFF46 ) :
    case ( kNuisEDBANFF47 ) :
    case ( kNuisEDBANFF48 ) :
    case ( kNuisEDBANFF49 ) :
    case ( kNuisEDBANFF50 ) :
    case ( kNuisEDBANFF51 ) :
    case ( kNuisEDBANFF52 ) :
    case ( kNuisEDBANFF53 ) :
#endif
#endif
#endif
    {
      int ip = (int)p - (int)kNuisEDBANFF00;
      return Form("EigenDecomposition BANFF %2d", ip);
      break;
    }

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

    case ( kPhypDeltaCP               ) : return "$\\delta_{cp}$";                   break;
    case ( kPhypDm23sq                ) : return "$|\\Deltam^2_{32}|$";                  break;
    case ( kPhypDm12sq                ) : return "$\\Deltam^2_{12}$";                    break;

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
    case ( kNuisSuperkEnergyScale     ) : return "$f_{E; r}^{SK}$  ";                 break;    
    case ( kNuisSuperkEff1RmuCCQEoth  ) : return "$f_{CCQEoth; r}^{SK}$";             break;
    case ( kNuisSuperkEff1RmuCCQErc0u ) : return "$f_{CCQERC0; r}^{SK}$";             break;
    case ( kNuisSuperkEff1RmuCCQErc1u ) : return "$f_{CCQERC1; r}^{SK}$";             break;
    case ( kNuisSuperkEff1RmuCCQErc2u ) : return "$f_{CCQERC2; r}^{SK}$";             break;
    case ( kNuisSuperkEff1RmuCCnonQE  ) : return "$f_{CCnQE; r}^{SK}$";               break;
    case ( kNuisSuperkEff1RmuNC       ) : return "$f_{NC; r}^{SK}$";                  break;
    case ( kNuisSuperkEff1RmuNueCC    ) : return "$f_{\\nu_e; r}^{SK}$";              break;

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
#if BANFF_BINS <= 2
    case ( kNuisBANFF16 ) :
    case ( kNuisBANFF17 ) :
    case ( kNuisBANFF18 ) :
    case ( kNuisBANFF19 ) : 
    case ( kNuisBANFF20 ) : 
    case ( kNuisBANFF21 ) : 
    case ( kNuisBANFF22 ) :
    case ( kNuisBANFF23 ) :
    case ( kNuisBANFF24 ) :
    case ( kNuisBANFF25 ) :
    case ( kNuisBANFF26 ) :
    case ( kNuisBANFF27 ) :
    case ( kNuisBANFF28 ) :
    case ( kNuisBANFF29 ) :
    case ( kNuisBANFF30 ) :
#if BANFF_BINS == -1
    case ( kNuisBANFF31 ) :
    case ( kNuisBANFF32 ) :
    case ( kNuisBANFF33 ) :
    case ( kNuisBANFF34 ) :
    case ( kNuisBANFF35 ) :
    case ( kNuisBANFF36 ) :
    case ( kNuisBANFF37 ) :
    case ( kNuisBANFF38 ) :
    case ( kNuisBANFF39 ) :
    case ( kNuisBANFF40 ) :
    case ( kNuisBANFF41 ) :
    case ( kNuisBANFF42 ) :
    case ( kNuisBANFF43 ) :
    case ( kNuisBANFF44 ) :
    case ( kNuisBANFF45 ) :
    case ( kNuisBANFF46 ) :
    case ( kNuisBANFF47 ) :
    case ( kNuisBANFF48 ) :
#endif
#endif
#if BANFF_BINS == 4
    case ( kNuisBANFFMAQE    ) :
#else
    case ( kNuisBANFFMAQE    ) :
#endif
    case ( kNuisBANFFMARES   ) :
    case ( kNuisBANFFCCQEE0  ) :
#if BANFF_BINS >= 2
    case ( kNuisBANFFCCQEE1  ) :
    case ( kNuisBANFFCCQEE2  ) :
#endif
    case ( kNuisBANFFCC1piE0 ) :
#if BANFF_BINS >= 2
    case ( kNuisBANFFCC1piE1 ) :
#endif
#if BANFF_BINS <= 2
    case ( kNuisBANFFNC1pi0  ) :
#endif
#if BANFF_BINS == 5
    case ( kNuisBANFFMEC     ) : 
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
    case ( kNuisXSecENuShape   ) : return "Shape; CC1pi ENu";                break;
    case ( kNuisXSecBindE      ) : return "$f_{bindE; t, r}$";               break;

#if BANFF_BINS <= 1
    case ( kNuisXSecNormCCQEE1 ) : return "Norm, CC QE Norm, 1.5 - 3.5 GeV"; break; 
    case ( kNuisXSecNormCCQEE2 ) : return "Norm, CC QE Norm, $>$ 3.5 GeV";   break; 
    case ( kNuisXSecNormCC1piE1) : return "Norm, CC 1pi, $>$ 2.5 GeV";       break; 
#endif
    case ( kNuisXSecNormCCcoh  ) : return "$f_{CCcoh; t}$";                  break;
    case ( kNuisXSecNormNCcoh  ) : return "Norm, NC Coh";                    break;
    case ( kNuisXSecNormNC1piC ) : return "$f_{NC1\\pi^\\pm; t}$";           break;
    case ( kNuisXSecNormNCoth  ) : return "$f_{NCoth; t}$";                  break;
    case ( kNuisXSecNormNuEToNuMu  ) : return "$f_{CC\\nu_e; t}$";           break;
    case ( kNuisXSecNormAntiNuToNu ) : return "$f_{CC\\bar{\\nu}; t}$";      break;

    case ( kNuisFSI000 ) :
    case ( kNuisFSI001 ) :
    case ( kNuisFSI002 ) :
    case ( kNuisFSI003 ) :
    case ( kNuisFSI004 ) :
    case ( kNuisFSI005 ) :
#if FSI_BINS >= 2
    case ( kNuisFSI006 ) :
    case ( kNuisFSI007 ) :
    case ( kNuisFSI008 ) :
    case ( kNuisFSI009 ) :
    case ( kNuisFSI010 ) :
    case ( kNuisFSI011 ) :
    case ( kNuisFSI012 ) :
    case ( kNuisFSI013 ) :
    case ( kNuisFSI014 ) :
#if FSI_BINS >= 3
    case ( kNuisFSI015 ) :
    case ( kNuisFSI016 ) :
    case ( kNuisFSI017 ) :
    case ( kNuisFSI018 ) :
    case ( kNuisFSI019 ) :
    case ( kNuisFSI020 ) :
    case ( kNuisFSI021 ) :
    case ( kNuisFSI022 ) :
    case ( kNuisFSI023 ) :
    case ( kNuisFSI024 ) :
    case ( kNuisFSI025 ) :
    case ( kNuisFSI026 ) :
    case ( kNuisFSI027 ) :
    case ( kNuisFSI028 ) :
    case ( kNuisFSI029 ) :
    case ( kNuisFSI030 ) :
    case ( kNuisFSI031 ) :
    case ( kNuisFSI032 ) :
    case ( kNuisFSI033 ) :
    case ( kNuisFSI034 ) :
    case ( kNuisFSI035 ) :
    case ( kNuisFSI036 ) :
    case ( kNuisFSI037 ) :
    case ( kNuisFSI038 ) :
    case ( kNuisFSI039 ) :
    case ( kNuisFSI040 ) :
    case ( kNuisFSI041 ) :
    case ( kNuisFSI042 ) :
    case ( kNuisFSI043 ) :
    case ( kNuisFSI044 ) :
    case ( kNuisFSI045 ) :
    case ( kNuisFSI046 ) :
    case ( kNuisFSI047 ) :
    case ( kNuisFSI048 ) :
#endif
#endif
    {
      int ip = (int)p - (int)kNuisFSI000;
      return Form("$f_{%d; t, r}^{FSI}$", ip);
      break;
    }//FSI

    case ( kNuisUnBANFF00 ) :
    case ( kNuisUnBANFF01 ) :
    case ( kNuisUnBANFF02 ) :
    case ( kNuisUnBANFF03 ) :
    case ( kNuisUnBANFF04 ) :
    case ( kNuisUnBANFF05 ) :
    case ( kNuisUnBANFF06 ) :
    case ( kNuisUnBANFF07 ) :
    case ( kNuisUnBANFF08 ) :
    case ( kNuisUnBANFF09 ) :
    case ( kNuisUnBANFF10 ) :
    case ( kNuisUnBANFF11 ) :
    case ( kNuisUnBANFF12 ) :
    case ( kNuisUnBANFF13 ) :
    case ( kNuisUnBANFF14 ) :
    case ( kNuisUnBANFF15 ) :
    case ( kNuisUnBANFF16 ) :
    case ( kNuisUnBANFF17 ) :
    case ( kNuisUnBANFF18 ) :
    case ( kNuisUnBANFF19 ) : 
    case ( kNuisUnBANFF20 ) : 
    case ( kNuisUnBANFF21 ) : 
    case ( kNuisUnBANFF22 ) :
#if BANFF_BINS <= 2
    case ( kNuisUnBANFF23 ) :
    case ( kNuisUnBANFF24 ) :
    case ( kNuisUnBANFF25 ) :
    case ( kNuisUnBANFF26 ) :
    case ( kNuisUnBANFF27 ) :
    case ( kNuisUnBANFF28 ) :
    case ( kNuisUnBANFF29 ) :
    case ( kNuisUnBANFF30 ) :
    case ( kNuisUnBANFF31 ) :
    case ( kNuisUnBANFF32 ) :
    case ( kNuisUnBANFF33 ) :
    case ( kNuisUnBANFF34 ) :
    case ( kNuisUnBANFF35 ) :
#if BANFF_BINS != 1
    case ( kNuisUnBANFF36 ) :
    case ( kNuisUnBANFF37 ) :
    case ( kNuisUnBANFF38 ) :
#if BANFF_BINS == -1
    case ( kNuisUnBANFF39 ) :
    case ( kNuisUnBANFF40 ) :
    case ( kNuisUnBANFF41 ) :
    case ( kNuisUnBANFF42 ) :
    case ( kNuisUnBANFF43 ) :
    case ( kNuisUnBANFF44 ) :
    case ( kNuisUnBANFF45 ) :
    case ( kNuisUnBANFF46 ) :
    case ( kNuisUnBANFF47 ) :
    case ( kNuisUnBANFF48 ) :
    case ( kNuisUnBANFF49 ) :
    case ( kNuisUnBANFF50 ) :
    case ( kNuisUnBANFF51 ) :
    case ( kNuisUnBANFF52 ) :
    case ( kNuisUnBANFF53 ) :
#endif
#endif
#endif
    {
      int ip = (int)p - (int)kNuisUnBANFF00;
      return Form("Uncorrelated BANFF %2d", ip);
      break;
    }

    case ( kNuisEDBANFF00 ) :
    case ( kNuisEDBANFF01 ) :
    case ( kNuisEDBANFF02 ) :
    case ( kNuisEDBANFF03 ) :
    case ( kNuisEDBANFF04 ) :
    case ( kNuisEDBANFF05 ) :
    case ( kNuisEDBANFF06 ) :
    case ( kNuisEDBANFF07 ) :
    case ( kNuisEDBANFF08 ) :
    case ( kNuisEDBANFF09 ) :
    case ( kNuisEDBANFF10 ) :
    case ( kNuisEDBANFF11 ) :
    case ( kNuisEDBANFF12 ) :
    case ( kNuisEDBANFF13 ) :
    case ( kNuisEDBANFF14 ) :
    case ( kNuisEDBANFF15 ) :
    case ( kNuisEDBANFF16 ) :
    case ( kNuisEDBANFF17 ) :
    case ( kNuisEDBANFF18 ) :
    case ( kNuisEDBANFF19 ) : 
    case ( kNuisEDBANFF20 ) : 
    case ( kNuisEDBANFF21 ) : 
    case ( kNuisEDBANFF22 ) :
#if BANFF_BINS <= 2
    case ( kNuisEDBANFF23 ) :
    case ( kNuisEDBANFF24 ) :
    case ( kNuisEDBANFF25 ) :
    case ( kNuisEDBANFF26 ) :
    case ( kNuisEDBANFF27 ) :
    case ( kNuisEDBANFF28 ) :
    case ( kNuisEDBANFF29 ) :
    case ( kNuisEDBANFF30 ) :
    case ( kNuisEDBANFF31 ) :
    case ( kNuisEDBANFF32 ) :
    case ( kNuisEDBANFF33 ) :
    case ( kNuisEDBANFF34 ) :
    case ( kNuisEDBANFF35 ) :
#if BANFF_BINS != 1
    case ( kNuisEDBANFF36 ) :
    case ( kNuisEDBANFF37 ) :
    case ( kNuisEDBANFF38 ) :
#if BANFF_BINS == -1
    case ( kNuisEDBANFF39 ) :
    case ( kNuisEDBANFF40 ) :
    case ( kNuisEDBANFF41 ) :
    case ( kNuisEDBANFF42 ) :
    case ( kNuisEDBANFF43 ) :
    case ( kNuisEDBANFF44 ) :
    case ( kNuisEDBANFF45 ) :
    case ( kNuisEDBANFF46 ) :
    case ( kNuisEDBANFF47 ) :
    case ( kNuisEDBANFF48 ) :
    case ( kNuisEDBANFF49 ) :
    case ( kNuisEDBANFF50 ) :
    case ( kNuisEDBANFF51 ) :
    case ( kNuisEDBANFF52 ) :
    case ( kNuisEDBANFF53 ) :
#endif
#endif
#endif
    {
      int ip = (int)p - (int)kNuisEDBANFF00;
      return Form("EigenDecomposition BANFF %2d", ip);
      break;
    }

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
           case(kNuMixSinsq2Theta) : return "sin{^2}(2#theta_{23})"; break;
           case(kNuMixSinsqTheta ) : return "sin{^2}(#theta_{23})";   break;
           default                 : return "--";                break;
         }
    }
    break;
    case ( kPhypMixing13 ) : 
    {
        NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
        switch(omp) {
           case(kNuMixSinsq2Theta) : return "sin{^2}(2#theta_{13})";  break;
           case(kNuMixSinsqTheta ) : return "sin{^2}(#theta_{13})";    break;
           default                 : return "--";                 break;
         }
    }
    break;
    case ( kPhypMixing12 ) :      
    {
        NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
        switch(omp) {
           case(kNuMixSinsq2Theta) : return "sin{^2}(2#theta_{12})";  break;
           case(kNuMixSinsqTheta ) : return "sin{^2}(#theta_{12})";    break;
           default                 : return "--";                 break;
         }
    }
    break;

    case ( kPhypDeltaCP               ) : return "#delta_{cp}";                   break;
    case ( kPhypDm23sq                ) : return "|#Deltam{^2}_{32}|";                  break;
    case ( kPhypDm12sq                ) : return "#Deltam{^2}_{12}";                    break;

    case ( kPhypAlphaDecay            ) : return "#alpha_{decay}";                break;
    case ( kPhypSinsq2ThetaDecay      ) : return "sin{^2}(2#theta_{decay})";       break;
    case ( kPhypMusqDecoh             ) : return "#mu{^2}_{decoh}";                 break;
    case ( kPhypSinsq2ThetaDecoh      ) : return "sin{^2}(2#theta_{decoh})";       break;

    /*
     * 2011-analysis nuisance parameters
     */

    //
    // SK detector
    //
    case ( kNuisSuperkEnergyScale     ) : return "f_{E; r}^{SK}  ";                 break;    
    case ( kNuisSuperkEff1RmuCCQEoth  ) : return "f_{CCQEoth; r}^{SK}";             break;
    case ( kNuisSuperkEff1RmuCCQErc0u ) : return "f_{CCQERC0; r}^{SK}";             break;
    case ( kNuisSuperkEff1RmuCCQErc1u ) : return "f_{CCQERC1; r}^{SK}";             break;
    case ( kNuisSuperkEff1RmuCCQErc2u ) : return "f_{CCQERC2; r}^{SK}";             break;
    case ( kNuisSuperkEff1RmuCCnonQE  ) : return "f_{CCnQE; r}^{SK}";               break;
    case ( kNuisSuperkEff1RmuNC       ) : return "f_{NC; r}^{SK}";                  break;
    case ( kNuisSuperkEff1RmuNueCC    ) : return "f_{#nu_e; r}^{SK}";              break;

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
#if BANFF_BINS <= 2
    case ( kNuisBANFF16 ) :
    case ( kNuisBANFF17 ) :
    case ( kNuisBANFF18 ) :
    case ( kNuisBANFF19 ) : 
    case ( kNuisBANFF20 ) : 
    case ( kNuisBANFF21 ) : 
    case ( kNuisBANFF22 ) :
    case ( kNuisBANFF23 ) :
    case ( kNuisBANFF24 ) :
    case ( kNuisBANFF25 ) :
    case ( kNuisBANFF26 ) :
    case ( kNuisBANFF27 ) :
    case ( kNuisBANFF28 ) :
    case ( kNuisBANFF29 ) :
    case ( kNuisBANFF30 ) :
#if BANFF_BINS == -1
    case ( kNuisBANFF31 ) :
    case ( kNuisBANFF32 ) :
    case ( kNuisBANFF33 ) :
    case ( kNuisBANFF34 ) :
    case ( kNuisBANFF35 ) :
    case ( kNuisBANFF36 ) :
    case ( kNuisBANFF37 ) :
    case ( kNuisBANFF38 ) :
    case ( kNuisBANFF39 ) :
    case ( kNuisBANFF40 ) :
    case ( kNuisBANFF41 ) :
    case ( kNuisBANFF42 ) :
    case ( kNuisBANFF43 ) :
    case ( kNuisBANFF44 ) :
    case ( kNuisBANFF45 ) :
    case ( kNuisBANFF46 ) :
    case ( kNuisBANFF47 ) :
    case ( kNuisBANFF48 ) :
#endif
#endif
#if BANFF_BINS == 4
    case ( kNuisBANFFMAQE    ) :
#else
    case ( kNuisBANFFMAQE    ) :
#endif
    case ( kNuisBANFFMARES   ) :
    case ( kNuisBANFFCCQEE0  ) :
#if BANFF_BINS >= 2
    case ( kNuisBANFFCCQEE1  ) :
    case ( kNuisBANFFCCQEE2  ) :
#endif
    case ( kNuisBANFFCC1piE0 ) :
#if BANFF_BINS >= 2
    case ( kNuisBANFFCC1piE1 ) :
#endif
#if BANFF_BINS <= 2
    case ( kNuisBANFFNC1pi0  ) :
#endif
#if BANFF_BINS == 5
    case ( kNuisBANFFMEC     ) : 
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
    case ( kNuisXSecENuShape   ) : return "Shape; CC1pi ENu";                break;
    case ( kNuisXSecBindE      ) : return "f_{bindE; t, r}";               break;

#if BANFF_BINS <= 1
    case ( kNuisXSecNormCCQEE1 ) : return "Norm, CC QE Norm, 1.5 - 3.5 GeV"; break; 
    case ( kNuisXSecNormCCQEE2 ) : return "Norm, CC QE Norm, > 3.5 GeV";   break; 
    case ( kNuisXSecNormCC1piE1) : return "Norm, CC 1pi, > 2.5 GeV";       break; 
#endif
    case ( kNuisXSecNormCCcoh  ) : return "f_{CCcoh; t}";                  break;
    case ( kNuisXSecNormNCcoh  ) : return "Norm, NC Coh";                    break;
    case ( kNuisXSecNormNC1piC ) : return "f_{NC1#pi#pm; t}";           break;
    case ( kNuisXSecNormNCoth  ) : return "f_{NCoth; t}";                  break;
    case ( kNuisXSecNormNuEToNuMu  ) : return "f_{CC#nue; t}";           break;
    case ( kNuisXSecNormAntiNuToNu ) : return "f_{CC#bar{#nu}; t}";      break;

    case ( kNuisFSI000 ) :
    case ( kNuisFSI001 ) :
    case ( kNuisFSI002 ) :
    case ( kNuisFSI003 ) :
    case ( kNuisFSI004 ) :
    case ( kNuisFSI005 ) :
#if FSI_BINS >= 2
    case ( kNuisFSI006 ) :
    case ( kNuisFSI007 ) :
    case ( kNuisFSI008 ) :
    case ( kNuisFSI009 ) :
    case ( kNuisFSI010 ) :
    case ( kNuisFSI011 ) :
    case ( kNuisFSI012 ) :
    case ( kNuisFSI013 ) :
    case ( kNuisFSI014 ) :
#if FSI_BINS >= 3
    case ( kNuisFSI015 ) :
    case ( kNuisFSI016 ) :
    case ( kNuisFSI017 ) :
    case ( kNuisFSI018 ) :
    case ( kNuisFSI019 ) :
    case ( kNuisFSI020 ) :
    case ( kNuisFSI021 ) :
    case ( kNuisFSI022 ) :
    case ( kNuisFSI023 ) :
    case ( kNuisFSI024 ) :
    case ( kNuisFSI025 ) :
    case ( kNuisFSI026 ) :
    case ( kNuisFSI027 ) :
    case ( kNuisFSI028 ) :
    case ( kNuisFSI029 ) :
    case ( kNuisFSI030 ) :
    case ( kNuisFSI031 ) :
    case ( kNuisFSI032 ) :
    case ( kNuisFSI033 ) :
    case ( kNuisFSI034 ) :
    case ( kNuisFSI035 ) :
    case ( kNuisFSI036 ) :
    case ( kNuisFSI037 ) :
    case ( kNuisFSI038 ) :
    case ( kNuisFSI039 ) :
    case ( kNuisFSI040 ) :
    case ( kNuisFSI041 ) :
    case ( kNuisFSI042 ) :
    case ( kNuisFSI043 ) :
    case ( kNuisFSI044 ) :
    case ( kNuisFSI045 ) :
    case ( kNuisFSI046 ) :
    case ( kNuisFSI047 ) :
    case ( kNuisFSI048 ) :
#endif
#endif
    {
      int ip = (int)p - (int)kNuisFSI000;
      return Form("f_{%d; t, r}^{FSI}", ip);
      break;
    }//FSI

    case ( kNuisUnBANFF00 ) :
    case ( kNuisUnBANFF01 ) :
    case ( kNuisUnBANFF02 ) :
    case ( kNuisUnBANFF03 ) :
    case ( kNuisUnBANFF04 ) :
    case ( kNuisUnBANFF05 ) :
    case ( kNuisUnBANFF06 ) :
    case ( kNuisUnBANFF07 ) :
    case ( kNuisUnBANFF08 ) :
    case ( kNuisUnBANFF09 ) :
    case ( kNuisUnBANFF10 ) :
    case ( kNuisUnBANFF11 ) :
    case ( kNuisUnBANFF12 ) :
    case ( kNuisUnBANFF13 ) :
    case ( kNuisUnBANFF14 ) :
    case ( kNuisUnBANFF15 ) :
    case ( kNuisUnBANFF16 ) :
    case ( kNuisUnBANFF17 ) :
    case ( kNuisUnBANFF18 ) :
    case ( kNuisUnBANFF19 ) : 
    case ( kNuisUnBANFF20 ) : 
    case ( kNuisUnBANFF21 ) : 
    case ( kNuisUnBANFF22 ) :
#if BANFF_BINS <= 2
    case ( kNuisUnBANFF23 ) :
    case ( kNuisUnBANFF24 ) :
    case ( kNuisUnBANFF25 ) :
    case ( kNuisUnBANFF26 ) :
    case ( kNuisUnBANFF27 ) :
    case ( kNuisUnBANFF28 ) :
    case ( kNuisUnBANFF29 ) :
    case ( kNuisUnBANFF30 ) :
    case ( kNuisUnBANFF31 ) :
    case ( kNuisUnBANFF32 ) :
    case ( kNuisUnBANFF33 ) :
    case ( kNuisUnBANFF34 ) :
    case ( kNuisUnBANFF35 ) :
#if BANFF_BINS != 1
    case ( kNuisUnBANFF36 ) :
    case ( kNuisUnBANFF37 ) :
    case ( kNuisUnBANFF38 ) :
#if BANFF_BINS == -1
    case ( kNuisUnBANFF39 ) :
    case ( kNuisUnBANFF40 ) :
    case ( kNuisUnBANFF41 ) :
    case ( kNuisUnBANFF42 ) :
    case ( kNuisUnBANFF43 ) :
    case ( kNuisUnBANFF44 ) :
    case ( kNuisUnBANFF45 ) :
    case ( kNuisUnBANFF46 ) :
    case ( kNuisUnBANFF47 ) :
    case ( kNuisUnBANFF48 ) :
    case ( kNuisUnBANFF49 ) :
    case ( kNuisUnBANFF50 ) :
    case ( kNuisUnBANFF51 ) :
    case ( kNuisUnBANFF52 ) :
    case ( kNuisUnBANFF53 ) :
#endif
#endif
#endif
    {
      int ip = (int)p - (int)kNuisUnBANFF00;
      return Form("Uncorrelated BANFF %2d", ip);
      break;
    }

    case ( kNuisEDBANFF00 ) :
    case ( kNuisEDBANFF01 ) :
    case ( kNuisEDBANFF02 ) :
    case ( kNuisEDBANFF03 ) :
    case ( kNuisEDBANFF04 ) :
    case ( kNuisEDBANFF05 ) :
    case ( kNuisEDBANFF06 ) :
    case ( kNuisEDBANFF07 ) :
    case ( kNuisEDBANFF08 ) :
    case ( kNuisEDBANFF09 ) :
    case ( kNuisEDBANFF10 ) :
    case ( kNuisEDBANFF11 ) :
    case ( kNuisEDBANFF12 ) :
    case ( kNuisEDBANFF13 ) :
    case ( kNuisEDBANFF14 ) :
    case ( kNuisEDBANFF15 ) :
    case ( kNuisEDBANFF16 ) :
    case ( kNuisEDBANFF17 ) :
    case ( kNuisEDBANFF18 ) :
    case ( kNuisEDBANFF19 ) : 
    case ( kNuisEDBANFF20 ) : 
    case ( kNuisEDBANFF21 ) : 
    case ( kNuisEDBANFF22 ) :
#if BANFF_BINS <= 2
    case ( kNuisEDBANFF23 ) :
    case ( kNuisEDBANFF24 ) :
    case ( kNuisEDBANFF25 ) :
    case ( kNuisEDBANFF26 ) :
    case ( kNuisEDBANFF27 ) :
    case ( kNuisEDBANFF28 ) :
    case ( kNuisEDBANFF29 ) :
    case ( kNuisEDBANFF30 ) :
    case ( kNuisEDBANFF31 ) :
    case ( kNuisEDBANFF32 ) :
    case ( kNuisEDBANFF33 ) :
    case ( kNuisEDBANFF34 ) :
    case ( kNuisEDBANFF35 ) :
#if BANFF_BINS != 1
    case ( kNuisEDBANFF36 ) :
    case ( kNuisEDBANFF37 ) :
    case ( kNuisEDBANFF38 ) :
#if BANFF_BINS == -1
    case ( kNuisEDBANFF39 ) :
    case ( kNuisEDBANFF40 ) :
    case ( kNuisEDBANFF41 ) :
    case ( kNuisEDBANFF42 ) :
    case ( kNuisEDBANFF43 ) :
    case ( kNuisEDBANFF44 ) :
    case ( kNuisEDBANFF45 ) :
    case ( kNuisEDBANFF46 ) :
    case ( kNuisEDBANFF47 ) :
    case ( kNuisEDBANFF48 ) :
    case ( kNuisEDBANFF49 ) :
    case ( kNuisEDBANFF50 ) :
    case ( kNuisEDBANFF51 ) :
    case ( kNuisEDBANFF52 ) :
    case ( kNuisEDBANFF53 ) :
#endif
#endif
#endif
    {
      int ip = (int)p - (int)kNuisEDBANFF00;
      return Form("EigenDecomposition BANFF %2d", ip);
      break;
    }

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
    case ( kPhypDm23sq                ) : return "|dm32sq|";                    break;
    case ( kPhypDm12sq                ) : return "dm12sq";                      break;

    case ( kPhypAlphaDecay            ) : return "alpha_decay";                 break;
    case ( kPhypSinsq2ThetaDecay      ) : return "sinsq2theta_decay";           break;
    case ( kPhypMusqDecoh             ) : return "musq_decoh";                  break;
    case ( kPhypSinsq2ThetaDecoh      ) : return "sinsq2theta_decoh";           break;

    /*
     * 2011-analysis nuisance parameters
     */

    case ( kNuisSuperkEnergyScale     ) : return "f_sk_e_scale";                   break;    
    case ( kNuisSuperkEff1RmuCCQEoth  ) : return "f_sk_eff_1Rmu_ccqe_oth";         break;
    case ( kNuisSuperkEff1RmuCCQErc0u ) : return "f_sk_eff_1Rmu_ccqe_rc0";         break;
    case ( kNuisSuperkEff1RmuCCQErc1u ) : return "f_sk_eff_1Rmu_ccqe_rc1";         break;
    case ( kNuisSuperkEff1RmuCCQErc2u ) : return "f_sk_eff_1Rmu_ccqe_rc2";         break;
    case ( kNuisSuperkEff1RmuCCnonQE  ) : return "f_sk_eff_1Rmu_ccnonqe";          break;
    case ( kNuisSuperkEff1RmuNC       ) : return "f_sk_eff_1Rmu_nc";               break;
    case ( kNuisSuperkEff1RmuNueCC    ) : return "f_sk_eff_1Rmu_nuecc";            break;

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
#if BANFF_BINS <= 2
    case ( kNuisBANFF16               ) : return "f_banff16";      break;
    case ( kNuisBANFF17               ) : return "f_banff17";      break;
    case ( kNuisBANFF18               ) : return "f_banff18";      break;
    case ( kNuisBANFF19               ) : return "f_banff19";      break;
    case ( kNuisBANFF20               ) : return "f_banff20";      break;
    case ( kNuisBANFF21               ) : return "f_banff21";      break;
    case ( kNuisBANFF22               ) : return "f_banff22";      break;
    case ( kNuisBANFF23               ) : return "f_banff23";      break;
    case ( kNuisBANFF24               ) : return "f_banff24";      break;
    case ( kNuisBANFF25               ) : return "f_banff25";      break;
    case ( kNuisBANFF26               ) : return "f_banff26";      break;
    case ( kNuisBANFF27               ) : return "f_banff27";      break;
    case ( kNuisBANFF28               ) : return "f_banff28";      break;
    case ( kNuisBANFF29               ) : return "f_banff29";      break;
    case ( kNuisBANFF30               ) : return "f_banff30";      break;
#if BANFF_BINS == -1 
    case ( kNuisBANFF31               ) : return "f_banff31";      break;
    case ( kNuisBANFF32               ) : return "f_banff32";      break;
    case ( kNuisBANFF33               ) : return "f_banff33";      break;
    case ( kNuisBANFF34               ) : return "f_banff34";      break;
    case ( kNuisBANFF35               ) : return "f_banff35";      break;
    case ( kNuisBANFF36               ) : return "f_banff36";      break;
    case ( kNuisBANFF37               ) : return "f_banff37";      break;
    case ( kNuisBANFF38               ) : return "f_banff38";      break;
    case ( kNuisBANFF39               ) : return "f_banff39";      break;
    case ( kNuisBANFF40               ) : return "f_banff40";      break;
    case ( kNuisBANFF41               ) : return "f_banff41";      break;
    case ( kNuisBANFF42               ) : return "f_banff42";      break;
    case ( kNuisBANFF43               ) : return "f_banff43";      break;
    case ( kNuisBANFF44               ) : return "f_banff44";      break;
    case ( kNuisBANFF45               ) : return "f_banff45";      break;
    case ( kNuisBANFF46               ) : return "f_banff46";      break;
    case ( kNuisBANFF47               ) : return "f_banff47";      break;
    case ( kNuisBANFF48               ) : return "f_banff48";      break;
#endif
#endif

#if BANFF_BINS == 4
    case ( kNuisBANFFMAQE             ) : return "f_banff1overmaqe2";  break;
#else
    case ( kNuisBANFFMAQE             ) : return "f_banffmaqe";    break;
#endif
    case ( kNuisBANFFMARES            ) : return "f_banffmares";   break;
    case ( kNuisBANFFCCQEE0           ) : return "f_banffccqee0";  break;
#if BANFF_BINS >= 2
    case ( kNuisBANFFCCQEE1           ) : return "f_banffccqee1";  break;
    case ( kNuisBANFFCCQEE2           ) : return "f_banffccqee2";  break;
#endif
    case ( kNuisBANFFCC1piE0          ) : return "f_banffcc1pie0"; break;
#if BANFF_BINS >= 2
    case ( kNuisBANFFCC1piE1          ) : return "f_banffcc1pie1"; break;
#endif
#if BANFF_BINS <= 2
    case ( kNuisBANFFNC1pi0           ) : return "f_banffnc1pi0";  break;
#endif
#if BANFF_BINS == 5
    case ( kNuisBANFFMEC              ) : return "f_banffmec";     break;
#endif

    case ( kNuisXSecFermiMom          ) : return "f_pfermi";       break;
    case ( kNuisXSecSpecFunc          ) : return "f_specfunc";     break; 
    case ( kNuisXSecWshape            ) : return "f_wshape";       break;
    case ( kNuisXSecShapeCCoth        ) : return "f_shapeCCoth";   break;
    case ( kNuisXSecPiLessDelta       ) : return "f_pilessdelta";  break;
    case ( kNuisXSecENuShape          ) : return "f_enushape";     break;
    case ( kNuisXSecBindE             ) : return "f_bindE";        break;

#if BANFF_BINS <= 1
    case ( kNuisXSecNormCCQEE1        ) : return "f_norm_ccqe_1";  break; 
    case ( kNuisXSecNormCCQEE2        ) : return "f_norm_ccqe_2";  break; 
    case ( kNuisXSecNormCC1piE1       ) : return "f_norm_cc1pi_1"; break;
#endif
    case ( kNuisXSecNormCCcoh         ) : return "f_norm_cccoh";   break;
    case ( kNuisXSecNormNC1piC        ) : return "f_norm_nc1piC";  break;
    case ( kNuisXSecNormNCcoh         ) : return "f_norm_nccoh";   break;
    case ( kNuisXSecNormNCoth         ) : return "f_norm_ncoth";   break;
    case ( kNuisXSecNormNuEToNuMu     ) : return "f_norm_nuetonumu";  break;
    case ( kNuisXSecNormAntiNuToNu    ) : return "f_norm_antinutonu"; break;

    case ( kNuisFSI000 ) : return "f_fsi000"; break;
    case ( kNuisFSI001 ) : return "f_fsi001"; break;
    case ( kNuisFSI002 ) : return "f_fsi002"; break;
    case ( kNuisFSI003 ) : return "f_fsi003"; break;
    case ( kNuisFSI004 ) : return "f_fsi004"; break;
    case ( kNuisFSI005 ) : return "f_fsi005"; break;
#if FSI_BINS >= 2
    case ( kNuisFSI006 ) : return "f_fsi006"; break;
    case ( kNuisFSI007 ) : return "f_fsi007"; break;
    case ( kNuisFSI008 ) : return "f_fsi008"; break;
    case ( kNuisFSI009 ) : return "f_fsi009"; break;
    case ( kNuisFSI010 ) : return "f_fsi010"; break;
    case ( kNuisFSI011 ) : return "f_fsi011"; break;
    case ( kNuisFSI012 ) : return "f_fsi012"; break;
    case ( kNuisFSI013 ) : return "f_fsi013"; break;
    case ( kNuisFSI014 ) : return "f_fsi014"; break;
#if FSI_BINS >= 3
    case ( kNuisFSI015 ) : return "f_fsi015"; break;
    case ( kNuisFSI016 ) : return "f_fsi016"; break;
    case ( kNuisFSI017 ) : return "f_fsi017"; break;
    case ( kNuisFSI018 ) : return "f_fsi018"; break;
    case ( kNuisFSI019 ) : return "f_fsi019"; break;
    case ( kNuisFSI020 ) : return "f_fsi020"; break;
    case ( kNuisFSI021 ) : return "f_fsi021"; break;
    case ( kNuisFSI022 ) : return "f_fsi022"; break;
    case ( kNuisFSI023 ) : return "f_fsi023"; break;
    case ( kNuisFSI024 ) : return "f_fsi024"; break;
    case ( kNuisFSI025 ) : return "f_fsi025"; break;
    case ( kNuisFSI026 ) : return "f_fsi026"; break;
    case ( kNuisFSI027 ) : return "f_fsi027"; break;
    case ( kNuisFSI028 ) : return "f_fsi028"; break;
    case ( kNuisFSI029 ) : return "f_fsi029"; break;
    case ( kNuisFSI030 ) : return "f_fsi030"; break;
    case ( kNuisFSI031 ) : return "f_fsi031"; break;
    case ( kNuisFSI032 ) : return "f_fsi032"; break;
    case ( kNuisFSI033 ) : return "f_fsi033"; break;
    case ( kNuisFSI034 ) : return "f_fsi034"; break;
    case ( kNuisFSI035 ) : return "f_fsi035"; break;
    case ( kNuisFSI036 ) : return "f_fsi036"; break;
    case ( kNuisFSI037 ) : return "f_fsi037"; break;
    case ( kNuisFSI038 ) : return "f_fsi038"; break;
    case ( kNuisFSI039 ) : return "f_fsi039"; break;
    case ( kNuisFSI040 ) : return "f_fsi040"; break;
    case ( kNuisFSI041 ) : return "f_fsi041"; break;
    case ( kNuisFSI042 ) : return "f_fsi042"; break;
    case ( kNuisFSI043 ) : return "f_fsi043"; break;
    case ( kNuisFSI044 ) : return "f_fsi044"; break;
    case ( kNuisFSI045 ) : return "f_fsi045"; break;
    case ( kNuisFSI046 ) : return "f_fsi046"; break;
    case ( kNuisFSI047 ) : return "f_fsi047"; break;
    case ( kNuisFSI048 ) : return "f_fsi048"; break;
#endif
#endif

    case ( kNuisUnBANFF00 ) : return "f_unbanff00"; break;
    case ( kNuisUnBANFF01 ) : return "f_unbanff01"; break;
    case ( kNuisUnBANFF02 ) : return "f_unbanff02"; break;
    case ( kNuisUnBANFF03 ) : return "f_unbanff03"; break;
    case ( kNuisUnBANFF04 ) : return "f_unbanff04"; break;
    case ( kNuisUnBANFF05 ) : return "f_unbanff05"; break;
    case ( kNuisUnBANFF06 ) : return "f_unbanff06"; break;
    case ( kNuisUnBANFF07 ) : return "f_unbanff07"; break;
    case ( kNuisUnBANFF08 ) : return "f_unbanff08"; break;
    case ( kNuisUnBANFF09 ) : return "f_unbanff09"; break;
    case ( kNuisUnBANFF10 ) : return "f_unbanff10"; break;
    case ( kNuisUnBANFF11 ) : return "f_unbanff11"; break;
    case ( kNuisUnBANFF12 ) : return "f_unbanff12"; break;
    case ( kNuisUnBANFF13 ) : return "f_unbanff13"; break;
    case ( kNuisUnBANFF14 ) : return "f_unbanff14"; break;
    case ( kNuisUnBANFF15 ) : return "f_unbanff15"; break;
    case ( kNuisUnBANFF16 ) : return "f_unbanff16"; break;
    case ( kNuisUnBANFF17 ) : return "f_unbanff17"; break;
    case ( kNuisUnBANFF18 ) : return "f_unbanff18"; break;
    case ( kNuisUnBANFF19 ) : return "f_unbanff19"; break; 
    case ( kNuisUnBANFF20 ) : return "f_unbanff20"; break; 
    case ( kNuisUnBANFF21 ) : return "f_unbanff21"; break; 
    case ( kNuisUnBANFF22 ) : return "f_unbanff22"; break;
#if BANFF_BINS <= 2
    case ( kNuisUnBANFF23 ) : return "f_unbanff23"; break;
    case ( kNuisUnBANFF24 ) : return "f_unbanff24"; break;
    case ( kNuisUnBANFF25 ) : return "f_unbanff25"; break;
    case ( kNuisUnBANFF26 ) : return "f_unbanff26"; break;
    case ( kNuisUnBANFF27 ) : return "f_unbanff27"; break;
    case ( kNuisUnBANFF28 ) : return "f_unbanff28"; break;
    case ( kNuisUnBANFF29 ) : return "f_unbanff29"; break;
    case ( kNuisUnBANFF30 ) : return "f_unbanff30"; break;
    case ( kNuisUnBANFF31 ) : return "f_unbanff31"; break;
    case ( kNuisUnBANFF32 ) : return "f_unbanff32"; break;
    case ( kNuisUnBANFF33 ) : return "f_unbanff33"; break;
    case ( kNuisUnBANFF34 ) : return "f_unbanff34"; break;
    case ( kNuisUnBANFF35 ) : return "f_unbanff35"; break;
#if BANFF_BINS != 1
    case ( kNuisUnBANFF36 ) : return "f_unbanff36"; break;
    case ( kNuisUnBANFF37 ) : return "f_unbanff37"; break;
    case ( kNuisUnBANFF38 ) : return "f_unbanff38"; break;
#if BANFF_BINS == -1
    case ( kNuisUnBANFF39 ) : return "f_unbanff39"; break;
    case ( kNuisUnBANFF40 ) : return "f_unbanff40"; break;
    case ( kNuisUnBANFF41 ) : return "f_unbanff41"; break;
    case ( kNuisUnBANFF42 ) : return "f_unbanff42"; break;
    case ( kNuisUnBANFF43 ) : return "f_unbanff43"; break;
    case ( kNuisUnBANFF44 ) : return "f_unbanff44"; break;
    case ( kNuisUnBANFF45 ) : return "f_unbanff45"; break;
    case ( kNuisUnBANFF46 ) : return "f_unbanff46"; break;
    case ( kNuisUnBANFF47 ) : return "f_unbanff47"; break;
    case ( kNuisUnBANFF48 ) : return "f_unbanff48"; break;
    case ( kNuisUnBANFF49 ) : return "f_unbanff49"; break;
    case ( kNuisUnBANFF50 ) : return "f_unbanff50"; break;
    case ( kNuisUnBANFF51 ) : return "f_unbanff51"; break;
    case ( kNuisUnBANFF52 ) : return "f_unbanff52"; break;
    case ( kNuisUnBANFF53 ) : return "f_unbanff53"; break;
#endif
#endif
#endif

    case ( kNuisEDBANFF00 ) : return "f_edbanff00"; break;
    case ( kNuisEDBANFF01 ) : return "f_edbanff01"; break;
    case ( kNuisEDBANFF02 ) : return "f_edbanff02"; break;
    case ( kNuisEDBANFF03 ) : return "f_edbanff03"; break;
    case ( kNuisEDBANFF04 ) : return "f_edbanff04"; break;
    case ( kNuisEDBANFF05 ) : return "f_edbanff05"; break;
    case ( kNuisEDBANFF06 ) : return "f_edbanff06"; break;
    case ( kNuisEDBANFF07 ) : return "f_edbanff07"; break;
    case ( kNuisEDBANFF08 ) : return "f_edbanff08"; break;
    case ( kNuisEDBANFF09 ) : return "f_edbanff09"; break;
    case ( kNuisEDBANFF10 ) : return "f_edbanff10"; break;
    case ( kNuisEDBANFF11 ) : return "f_edbanff11"; break;
    case ( kNuisEDBANFF12 ) : return "f_edbanff12"; break;
    case ( kNuisEDBANFF13 ) : return "f_edbanff13"; break;
    case ( kNuisEDBANFF14 ) : return "f_edbanff14"; break;
    case ( kNuisEDBANFF15 ) : return "f_edbanff15"; break;
    case ( kNuisEDBANFF16 ) : return "f_edbanff16"; break;
    case ( kNuisEDBANFF17 ) : return "f_edbanff17"; break;
    case ( kNuisEDBANFF18 ) : return "f_edbanff18"; break;
    case ( kNuisEDBANFF19 ) : return "f_edbanff19"; break; 
    case ( kNuisEDBANFF20 ) : return "f_edbanff20"; break; 
    case ( kNuisEDBANFF21 ) : return "f_edbanff21"; break; 
    case ( kNuisEDBANFF22 ) : return "f_edbanff22"; break;
#if BANFF_BINS <= 2
    case ( kNuisEDBANFF23 ) : return "f_edbanff23"; break;
    case ( kNuisEDBANFF24 ) : return "f_edbanff24"; break;
    case ( kNuisEDBANFF25 ) : return "f_edbanff25"; break;
    case ( kNuisEDBANFF26 ) : return "f_edbanff26"; break;
    case ( kNuisEDBANFF27 ) : return "f_edbanff27"; break;
    case ( kNuisEDBANFF28 ) : return "f_edbanff28"; break;
    case ( kNuisEDBANFF29 ) : return "f_edbanff29"; break;
    case ( kNuisEDBANFF30 ) : return "f_edbanff30"; break;
    case ( kNuisEDBANFF31 ) : return "f_edbanff31"; break;
    case ( kNuisEDBANFF32 ) : return "f_edbanff32"; break;
    case ( kNuisEDBANFF33 ) : return "f_edbanff33"; break;
    case ( kNuisEDBANFF34 ) : return "f_edbanff34"; break;
    case ( kNuisEDBANFF35 ) : return "f_edbanff35"; break;
#if BANFF_BINS != 1
    case ( kNuisEDBANFF36 ) : return "f_edbanff36"; break;
    case ( kNuisEDBANFF37 ) : return "f_edbanff37"; break;
    case ( kNuisEDBANFF38 ) : return "f_edbanff38"; break;
#if BANFF_BINS == -1
    case ( kNuisEDBANFF39 ) : return "f_edbanff39"; break;
    case ( kNuisEDBANFF40 ) : return "f_edbanff40"; break;
    case ( kNuisEDBANFF41 ) : return "f_edbanff41"; break;
    case ( kNuisEDBANFF42 ) : return "f_edbanff42"; break;
    case ( kNuisEDBANFF43 ) : return "f_edbanff43"; break;
    case ( kNuisEDBANFF44 ) : return "f_edbanff44"; break;
    case ( kNuisEDBANFF45 ) : return "f_edbanff45"; break;
    case ( kNuisEDBANFF46 ) : return "f_edbanff46"; break;
    case ( kNuisEDBANFF47 ) : return "f_edbanff47"; break;
    case ( kNuisEDBANFF48 ) : return "f_edbanff48"; break;
    case ( kNuisEDBANFF49 ) : return "f_edbanff49"; break;
    case ( kNuisEDBANFF50 ) : return "f_edbanff50"; break;
    case ( kNuisEDBANFF51 ) : return "f_edbanff51"; break;
    case ( kNuisEDBANFF52 ) : return "f_edbanff52"; break;
    case ( kNuisEDBANFF53 ) : return "f_edbanff53"; break;
#endif
#endif
#endif

    default : 
      return "--";
  }
  return "--";
}
//____________________________________________________________________________
string t2kana::util::AsStringCompact_Japan(PhysParamId_t p)
{
  switch(p) {

    case ( kNuisBANFFMAQE             ) : return "maqe";     break;
    case ( kNuisBANFFMARES            ) : return "mares";    break;
    case ( kNuisXSecFermiMom          ) : return "pf";       break;
    case ( kNuisXSecSpecFunc          ) : return "sf";       break; 
    case ( kNuisXSecWshape            ) : return "mdel";     break;
    case ( kNuisXSecShapeCCoth        ) : return "ccothshp"; break;
    case ( kNuisXSecPiLessDelta       ) : return "pdd";      break;
    case ( kNuisXSecENuShape          ) : return "pishp";    break;

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
     default                      : return "--";
  }
  return "";
}
//____________________________________________________________________________
string t2kana::util::AsString(t2kana::conv::SuperkMode_t m)
{
  switch(m) {
     case ( kSKmNuMuCCQE        ) : return "numu CC QE";                        break;
     case ( kSKmNuMuCC1pi       ) : return "numu CC 1pi";                       break;
     case ( kSKmNuMuCCcoh       ) : return "numu CC coh";                       break;
     case ( kSKmNuMuCCoth_mu    ) : return "numu CC oth (numu disap)";          break;
     case ( kSKmNuMuCCoth_e     ) : return "numu CC oth (nue app)";             break;
     case ( kSKmNuMuNC1piC      ) : return "numu/nutau NC 1piC";                break;
     case ( kSKmNuMuNC1pi0      ) : return "numu/nutau NC 1pi0";                break;
     case ( kSKmNuMuNCcoh       ) : return "numu/nutau NC coh";                 break;
     case ( kSKmNuMuNCoth_e     ) : return "numu/nutau NC oth (nue app)";       break;
     case ( kSKmNuMuNCoth_mu    ) : return "numu/nutau NC oth (numu app)";      break;
     case ( kSKmNuMuBarCCQE     ) : return "numubar CC QE";                     break;
     case ( kSKmNuMuBarCC1pi    ) : return "numubar CC 1pi";                    break;
     case ( kSKmNuMuBarCCcoh    ) : return "numubar CC coh";                    break;
     case ( kSKmNuMuBarCCoth_mu ) : return "numubar CC oth (numu disap)";       break;
     case ( kSKmNuMuBarCCoth_e  ) : return "numubar CC oth (nue app)";          break;
     case ( kSKmNuMuBarNC1piC   ) : return "numubar/nutaubar NC 1piC";          break;
     case ( kSKmNuMuBarNC1pi0   ) : return "numubar/nutaubar NC 1pi0";          break;
     case ( kSKmNuMuBarNCcoh    ) : return "numubar/nutaubar NC coh";           break;
     case ( kSKmNuMuBarNCoth_e  ) : return "numubar/nutaubar NC oth (nue app)"; break;
     case ( kSKmNuMuBarNCoth_mu ) : return "numubar/nutaubar NC oth (numu app)";break;
     case ( kSKmNuECCQE         ) : return "nue CC QE";                         break;
     case ( kSKmNuECC1pi        ) : return "nue CC 1pi";                        break;
     case ( kSKmNuECCcoh        ) : return "nue CC coh";                        break;
     case ( kSKmNuECCoth_mu     ) : return "nue CC oth (numu disap)";           break;
     case ( kSKmNuECCoth_e      ) : return "nue CC oth (nue app)";              break;
     case ( kSKmNuENC1piC       ) : return "nue NC 1piC";                       break;
     case ( kSKmNuENC1pi0       ) : return "nue NC 1pi0";                       break;
     case ( kSKmNuENCcoh        ) : return "nue NC coh";                        break;
     case ( kSKmNuENCoth_e      ) : return "nue NC oth (nue app)";              break;
     case ( kSKmNuENCoth_mu     ) : return "nue NC oth (numu app)";             break;
     case ( kSKmNuEBarCCQE      ) : return "nuebar CC QE";                      break;
     case ( kSKmNuEBarCC1pi     ) : return "nuebar CC 1pi";                     break;
     case ( kSKmNuEBarCCcoh     ) : return "nuebar CC coh";                     break;
     case ( kSKmNuEBarCCoth_mu  ) : return "nuebar CC oth (numu disap)";        break;
     case ( kSKmNuEBarCCoth_e   ) : return "nuebar CC oth (nue app)";           break;
     case ( kSKmNuEBarNC1piC    ) : return "nuebar NC 1piC";                    break;
     case ( kSKmNuEBarNC1pi0    ) : return "nuebar NC 1pi0";                    break;
     case ( kSKmNuEBarNCcoh     ) : return "nuebar NC coh";                     break;
     case ( kSKmNuEBarNCoth_e   ) : return "nuebar NC oth (nue app)";           break;
     case ( kSKmNuEBarNCoth_mu  ) : return "nuebar NC oth (numu app)";          break;
     case ( kSKmSigNuECCQE      ) : return "oscillation nue CC QE";                  break;
     case ( kSKmSigNuECC1pi     ) : return "oscillation nue CC 1pi";                 break;
     case ( kSKmSigNuECCcoh     ) : return "oscillation nue CC coh";                 break;
     case ( kSKmSigNuECCoth_mu  ) : return "oscillation nue CC oth (numu disap)";    break;
     case ( kSKmSigNuECCoth_e   ) : return "oscillation nue CC oth (nue app)";       break;
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
    case ( kFlxNominal     ) : return "nominal";          break;
    case ( kFlxExtTuning   ) : return "tuned (external)"; break;
    case ( kFlxFinalTuning ) : return "tuned (final)";    break;
    case ( kFlxCurrent     ) : return "current";          break;
    default                  : return "--";               break;
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
     default                      : return "--";
  }
  return "";
}
//____________________________________________________________________________
string t2kana::util::AsStringTex(t2kana::conv::SuperkMode_t m)
{
  switch(m) {
     case ( kSKmNuMuCCQE        ) : return "$\\nu_\\mu$ CCQE";                                               break;
     case ( kSKmNuMuCC1pi       ) : return "$\\nu_\\mu$ CC 1$\\pi$";                                         break;
     case ( kSKmNuMuCCcoh       ) : return "$\\nu_\\mu$ CC coherent";                                        break;
     case ( kSKmNuMuCCoth_mu    ) : return "$\\nu_\\mu$ CC other ($\\nu_\\mu$ disap)";                       break;
     case ( kSKmNuMuCCoth_e     ) : return "$\\nu_\\mu$ CC other ($\\nu_e$ app)";                            break;
     case ( kSKmNuMuNC1piC      ) : return "$\\nu_\\mu$/$\\nu_\\tau$ NC 1$\\pi^{+/-}$";                      break;
     case ( kSKmNuMuNC1pi0      ) : return "$\\nu_\\mu$/$\\nu_\\tau$ NC 1$\\pi^0$";                          break;
     case ( kSKmNuMuNCcoh       ) : return "$\\nu_\\mu$/$\\nu_\\tau$ NC coherent";                           break;
     case ( kSKmNuMuNCoth_e     ) : return "$\\nu_\\mu$/$\\nu_\\tau$ NC other ($\\nu_e$ app)";               break;
     case ( kSKmNuMuNCoth_mu    ) : return "$\\nu_\\mu$/$\\nu_\\tau$ NC other ($\\nu_\\mu$ app)";            break;
     case ( kSKmNuMuBarCCQE     ) : return "$\\bar{\\nu_\\mu}$ CCQE";                                        break;
     case ( kSKmNuMuBarCC1pi    ) : return "$\\bar{\\nu_\\mu}$ CC 1$\\pi$";                                  break;
     case ( kSKmNuMuBarCCcoh    ) : return "$\\bar{\\nu_\\mu}$ CC coherent";                                 break;
     case ( kSKmNuMuBarCCoth_mu ) : return "$\\bar{\\nu_\\mu}$ CC other ($\\nu_\\mu$ disap)";                break;
     case ( kSKmNuMuBarCCoth_e  ) : return "$\\bar{\\nu_\\mu}$ CC other ($\\nu_e$ app)";                     break;
     case ( kSKmNuMuBarNC1piC   ) : return "$\\bar{\\nu_\\mu}$/$\\bar{\\nu_\\tau}$ NC 1$\\pi^{+/-}$";        break;
     case ( kSKmNuMuBarNC1pi0   ) : return "$\\bar{\\nu_\\mu}$/$\\bar{\\nu_\\tau}$ NC 1$\\pi^0$";            break;
     case ( kSKmNuMuBarNCcoh    ) : return "$\\bar{\\nu_\\mu}$/$\\bar{\\nu_\\tau}$ NC coherent";             break;
     case ( kSKmNuMuBarNCoth_e  ) : return "$\\bar{\\nu_\\mu}$/$\\bar{\\nu_\\tau}$ NC other ($\\nu_e$ app)"; break;
     case ( kSKmNuMuBarNCoth_mu ) : return "$\\bar{\\nu_\\mu}$/$\\bar{\\nu_\\tau}$ NC other ($\\nu_\\mu$ app)";break;
     case ( kSKmNuECCQE         ) : return "$\\nu_e$ CCQE";                                                  break;
     case ( kSKmNuECC1pi        ) : return "$\\nu_e$ CC 1$\\pi$";                                            break;
     case ( kSKmNuECCcoh        ) : return "$\\nu_e$ CC coherent";                                           break;
     case ( kSKmNuECCoth_mu     ) : return "$\\nu_e$ CC other ($\\nu_\\mu$ disap)";                          break;
     case ( kSKmNuECCoth_e      ) : return "$\\nu_e$ CC other ($\\nu_e$ app)";                               break;
     case ( kSKmNuENC1piC       ) : return "$\\nu_e$ NC 1$\\pi^{+/-}$";                                      break;
     case ( kSKmNuENC1pi0       ) : return "$\\nu_e$ NC 1$\\pi^0$";                                          break;
     case ( kSKmNuENCcoh        ) : return "$\\nu_e$ NC coherent";                                           break;
     case ( kSKmNuENCoth_e      ) : return "$\\nu_e$ NC other ($\\nu_e$ app)";                               break;
     case ( kSKmNuENCoth_mu     ) : return "$\\nu_e$ NC other ($\\nu_\\mu$ app)";                            break;
     case ( kSKmNuEBarCCQE      ) : return "$\\bar{\\nu_e}$ CCQE";                                           break;
     case ( kSKmNuEBarCC1pi     ) : return "$\\bar{\\nu_e}$ CC 1$\\pi$";                                     break;
     case ( kSKmNuEBarCCcoh     ) : return "$\\bar{\\nu_e}$ CC coherent";                                    break;
     case ( kSKmNuEBarCCoth_mu  ) : return "$\\bar{\\nu_e}$ CC other ($\\nu_\\mu$ disap)";                   break;
     case ( kSKmNuEBarCCoth_e   ) : return "$\\bar{\\nu_e}$ CC other ($\\nu_e$ app)";                        break;
     case ( kSKmNuEBarNC1piC    ) : return "$\\bar{\\nu_e}$ NC 1$\\pi^{+/-}$";                               break;
     case ( kSKmNuEBarNC1pi0    ) : return "$\\bar{\\nu_e}$ NC 1$\\pi^0$";                                   break;
     case ( kSKmNuEBarNCcoh     ) : return "$\\bar{\\nu_e}$ NC coherent";                                    break;
     case ( kSKmNuEBarNCoth_e   ) : return "$\\bar{\\nu_e}$ NC other ($\\nu_e$ app)";                        break;
     case ( kSKmNuEBarNCoth_mu  ) : return "$\\bar{\\nu_e}$ NC other ($\\nu_\\mu$ app)";                     break;
     case ( kSKmSigNuECCQE      ) : return "oscillation $\\nu_e$ CCQE";                                           break;
     case ( kSKmSigNuECC1pi     ) : return "oscillation $\\nu_e$ CC 1$\\pi$";                                     break;
     case ( kSKmSigNuECCcoh     ) : return "oscillation $\\nu_e$ CC coherent";                                    break;
     case ( kSKmSigNuECCoth_mu  ) : return "oscillation $\\nu_e$ CC other ($\\nu_\\mu$ app)";                     break;
     case ( kSKmSigNuECCoth_e   ) : return "oscillation $\\nu_e$ CC other ($\\nu_e$ app)";                        break;
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
     default                      : return "";
  }
  return "";
}
//____________________________________________________________________________
string t2kana::util::AsStringRooTex(t2kana::conv::SuperkMode_t m)
{
  switch(m) {
     case ( kSKmNuMuCCQE        ) : return "#nu_{#mu} CCQE";                                break;
     case ( kSKmNuMuCC1pi       ) : return "#nu_{#mu} CC 1#pi";                             break;
     case ( kSKmNuMuCCcoh       ) : return "#nu_{#mu} CC coherent";                         break;
     case ( kSKmNuMuCCoth_mu    ) : return "#nu_{#mu} CC other";                            break;
     case ( kSKmNuMuCCoth_e     ) : return "#nu_{#mu} CC other";                            break;
     case ( kSKmNuMuNC1piC      ) : return "#nu_{#mu}/#nu_{#tau} NC 1#pi^{+/-}";            break;
     case ( kSKmNuMuNC1pi0      ) : return "#nu_{#mu}/#nu_{#tau} NC 1#pi^{0}";              break;
     case ( kSKmNuMuNCcoh       ) : return "#nu_{#mu}/#nu_{#tau} NC coherent";              break;
     case ( kSKmNuMuNCoth_e     ) : return "#nu_{#mu}/#nu_{#tau} NC other";                 break;
     case ( kSKmNuMuNCoth_mu    ) : return "#nu_{#mu}/#nu_{#tau} NC other";                 break;
     case ( kSKmNuMuBarCCQE     ) : return "#bar{#nu_{#mu}} CCQE";                          break;
     case ( kSKmNuMuBarCC1pi    ) : return "#bar{#nu_{#mu}} CC 1#pi";                       break;
     case ( kSKmNuMuBarCCcoh    ) : return "#bar{#nu_{#mu}} CC coherent";                   break;
     case ( kSKmNuMuBarCCoth_mu ) : return "#bar{#nu_{#mu}} CC other";                      break;
     case ( kSKmNuMuBarCCoth_e  ) : return "#bar{#nu_{#mu}} CC other";                      break;
     case ( kSKmNuMuBarNC1piC   ) : return "#bar{#nu_{#mu}}/#bar{#nu_{#tau}} NC 1#pi^{+/-}";break;
     case ( kSKmNuMuBarNC1pi0   ) : return "#bar{#nu_{#mu}}/#bar{#nu_{#tau}} NC 1#pi^{0}";  break;
     case ( kSKmNuMuBarNCcoh    ) : return "#bar{#nu_{#mu}}/#bar{#nu_{#tau}} NC coherent";  break;
     case ( kSKmNuMuBarNCoth_e  ) : return "#bar{#nu_{#mu}}/#bar{#nu_{#tau}} NC other";     break;
     case ( kSKmNuMuBarNCoth_mu ) : return "#bar{#nu_{#mu}}/#bar{#nu_{#tau}} NC other";     break;
     case ( kSKmNuECCQE         ) : return "#nu_{e} CCQE";                                  break;
     case ( kSKmNuECC1pi        ) : return "#nu_{e} CC 1#pi";                               break;
     case ( kSKmNuECCcoh        ) : return "#nu_{e} CC coherent";                           break;
     case ( kSKmNuECCoth_mu     ) : return "#nu_{e} CC other";                              break;
     case ( kSKmNuECCoth_e      ) : return "#nu_{e} CC other";                              break;
     case ( kSKmNuENC1piC       ) : return "#nu_{e} NC 1#pi^{+/-}";                         break;
     case ( kSKmNuENC1pi0       ) : return "#nu_{e} NC 1#pi^{0}";                           break;
     case ( kSKmNuENCcoh        ) : return "#nu_{e} NC coherent";                           break;
     case ( kSKmNuENCoth_e      ) : return "#nu_{e} NC other";                              break;
     case ( kSKmNuENCoth_mu     ) : return "#nu_{e} NC other";                              break;
     case ( kSKmNuEBarCCQE      ) : return "#bar{#nu_{e}} CCQE";                            break;
     case ( kSKmNuEBarCC1pi     ) : return "#bar{#nu_{e}} CC 1#pi";                         break;
     case ( kSKmNuEBarCCcoh     ) : return "#bar{#nu_{e}} CC coherent";                     break;
     case ( kSKmNuEBarCCoth_mu  ) : return "#bar{#nu_{e}} CC other";                        break;
     case ( kSKmNuEBarCCoth_e   ) : return "#bar{#nu_{e}} CC other";                        break;
     case ( kSKmNuEBarNC1piC    ) : return "#bar{#nu_{e}} NC 1#pi^{+/-}";                   break;
     case ( kSKmNuEBarNC1pi0    ) : return "#bar{#nu_{e}} NC 1#pi^{0}";                     break;
     case ( kSKmNuEBarNCcoh     ) : return "#bar{#nu_{e}} NC coherent";                     break;
     case ( kSKmNuEBarNCoth_e   ) : return "#bar{#nu_{e}} NC other";                        break;
     case ( kSKmNuEBarNCoth_mu  ) : return "#bar{#nu_{e}} NC other";                        break;
     case ( kSKmSigNuECCQE      ) : return "oscillation #nu_{e} CCQE";                           break;
     case ( kSKmSigNuECC1pi     ) : return "oscillation #nu_{e} CC 1#pi";                        break;
     case ( kSKmSigNuECCcoh     ) : return "oscillation #nu_{e} CC coherent";                    break;
     case ( kSKmSigNuECCoth_mu  ) : return "oscillation #nu_{e} CC other";                       break;
     case ( kSKmSigNuECCoth_e   ) : return "oscillation #nu_{e} CC other";                       break;
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
     case ( kSKsFCFV1RingMuLike ) : return "1Rmu"; break;
     case ( kSKsFCFV1RingELike  ) : return "1Re";  break;
     default                      : return "";
  }
  return "";
}
//____________________________________________________________________________
string t2kana::util::AsStringCompact (t2kana::conv::SuperkMode_t m)
{
// Compact string used for naming output histograms stored ROOT files
//
  switch(m) {
     case ( kSKmNuMuCCQE        ) : return "numuccqe";       break;
     case ( kSKmNuMuCC1pi       ) : return "numucc1pi";      break;
     case ( kSKmNuMuCCcoh       ) : return "numucccoh";      break;
     case ( kSKmNuMuCCoth_mu    ) : return "numuccoth";      break;
     case ( kSKmNuMuCCoth_e     ) : return "numuccoth";      break;
     case ( kSKmNuMuNC1piC      ) : return "numunc1piC";     break;
     case ( kSKmNuMuNC1pi0      ) : return "numunc1pi0";     break;
     case ( kSKmNuMuNCcoh       ) : return "numunccoh";      break;
     case ( kSKmNuMuNCoth_e     ) : return "numuncoth";      break;
     case ( kSKmNuMuNCoth_mu    ) : return "numuncoth";      break;
     case ( kSKmNuMuBarCCQE     ) : return "numubarccqe";    break;
     case ( kSKmNuMuBarCC1pi    ) : return "numubarcc1pi";   break;
     case ( kSKmNuMuBarCCcoh    ) : return "numubarcccoh";   break;
     case ( kSKmNuMuBarCCoth_mu ) : return "numubarccoth";   break;
     case ( kSKmNuMuBarCCoth_e  ) : return "numubarccoth";   break;
     case ( kSKmNuMuBarNC1piC   ) : return "numubarnc1piC";  break;
     case ( kSKmNuMuBarNC1pi0   ) : return "numubarnc1pi0";  break;
     case ( kSKmNuMuBarNCcoh    ) : return "numubarnccoh";   break;
     case ( kSKmNuMuBarNCoth_e  ) : return "numubarncoth";   break;
     case ( kSKmNuMuBarNCoth_mu ) : return "numubarncoth";   break;
     case ( kSKmNuECCQE         ) : return "nueccqe";        break;
     case ( kSKmNuECC1pi        ) : return "nuecc1pi";       break;
     case ( kSKmNuECCcoh        ) : return "nuecccoh";       break;
     case ( kSKmNuECCoth_mu     ) : return "nueccoth";       break;
     case ( kSKmNuECCoth_e      ) : return "nueccoth";       break;
     case ( kSKmNuENC1piC       ) : return "nuenc1piC";      break;
     case ( kSKmNuENC1pi0       ) : return "nuenc1pi0";      break;
     case ( kSKmNuENCcoh        ) : return "nuenccoh";       break;
     case ( kSKmNuENCoth_e      ) : return "nuencoth";       break;
     case ( kSKmNuENCoth_mu     ) : return "nuencoth";       break;
     case ( kSKmNuEBarCCQE      ) : return "nuebarccqe";     break;
     case ( kSKmNuEBarCC1pi     ) : return "nuebarcc1pi";    break;
     case ( kSKmNuEBarCCcoh     ) : return "nuebarcccoh";    break;
     case ( kSKmNuEBarCCoth_mu  ) : return "nuebarccoth";    break;
     case ( kSKmNuEBarCCoth_e   ) : return "nuebarccoth";    break;
     case ( kSKmNuEBarNC1piC    ) : return "nuebarnc1piC";   break;
     case ( kSKmNuEBarNC1pi0    ) : return "nuebarnc1pi0";   break;
     case ( kSKmNuEBarNCcoh     ) : return "nuebarnccoh";    break;
     case ( kSKmNuEBarNCoth_e   ) : return "nuebarncoth";    break;
     case ( kSKmNuEBarNCoth_mu  ) : return "nuebarncoth";    break;
     case ( kSKmSigNuECCQE      ) : return "signueccqe";     break;
     case ( kSKmSigNuECC1pi     ) : return "signuecc1pi";    break;
     case ( kSKmSigNuECCcoh     ) : return "signuecccoh";    break;
     case ( kSKmSigNuECCoth_mu  ) : return "signueccoth";    break;
     case ( kSKmSigNuECCoth_e   ) : return "signueccoth";    break;
     default                      : return "";
  }
  return "";
}
//____________________________________________________________________________
string t2kana::util::AsStringCompact_Japan (t2kana::conv::SuperkMode_t m)
{
// Compact string used for naming output histograms stored ROOT files
// using the Japanese analysis naming convention
//
  switch(m) {
     case ( kSKmNuMuCCQE        ) : return "numu_ccqe";       break;
     case ( kSKmNuMuCC1pi       ) : return "numu_ccpi";       break;
     case ( kSKmNuMuCCcoh       ) : return "numu_cccoh";      break;
     case ( kSKmNuMuCCoth_mu    ) : return "numu_ccoth";      break;
     case ( kSKmNuMuCCoth_e     ) : return "numu_ccoth";      break;
     case ( kSKmNuMuNC1piC      ) : return "numu_ncpi";       break;
     case ( kSKmNuMuNC1pi0      ) : return "numu_nc1pi0";     break;
     case ( kSKmNuMuNCcoh       ) : return "numu_nccoh";      break;
     case ( kSKmNuMuNCoth_e     ) : return "numu_ncoth";      break;
     case ( kSKmNuMuNCoth_mu    ) : return "numu_ncoth";      break;
     case ( kSKmNuMuBarCCQE     ) : return "numb_ccqe";    break;
     case ( kSKmNuMuBarCC1pi    ) : return "numb_ccpi";    break;
     case ( kSKmNuMuBarCCcoh    ) : return "numb_cccoh";   break;
     case ( kSKmNuMuBarCCoth_mu ) : return "numb_ccoth";   break;
     case ( kSKmNuMuBarCCoth_e  ) : return "numb_ccoth";   break;
     case ( kSKmNuMuBarNC1piC   ) : return "numb_ncpi";    break;
     case ( kSKmNuMuBarNC1pi0   ) : return "numb_nc1pi0";  break;
     case ( kSKmNuMuBarNCcoh    ) : return "numb_nccoh";   break;
     case ( kSKmNuMuBarNCoth_e  ) : return "numb_ncoth";   break;
     case ( kSKmNuMuBarNCoth_mu ) : return "numb_ncoth";   break;
     case ( kSKmNuECCQE         ) : return "nue_ccqe";        break;
     case ( kSKmNuECC1pi        ) : return "nue_ccpi";        break;
     case ( kSKmNuECCcoh        ) : return "nue_cccoh";       break;
     case ( kSKmNuECCoth_mu     ) : return "nue_ccoth";       break;
     case ( kSKmNuECCoth_e      ) : return "nue_ccoth";       break;
     case ( kSKmNuENC1piC       ) : return "nue_ncpi";        break;
     case ( kSKmNuENC1pi0       ) : return "nue_nc1pi0";      break;
     case ( kSKmNuENCcoh        ) : return "nue_nccoh";       break;
     case ( kSKmNuENCoth_e      ) : return "nue_ncoth";       break;
     case ( kSKmNuENCoth_mu     ) : return "nue_ncoth";       break;
     case ( kSKmNuEBarCCQE      ) : return "nueb_ccqe";     break;
     case ( kSKmNuEBarCC1pi     ) : return "nueb_ccpi";     break;
     case ( kSKmNuEBarCCcoh     ) : return "nueb_cccoh";    break;
     case ( kSKmNuEBarCCoth_mu  ) : return "nueb_ccoth";    break;
     case ( kSKmNuEBarCCoth_e   ) : return "nueb_ccoth";    break;
     case ( kSKmNuEBarNC1piC    ) : return "nueb_ncpi";     break;
     case ( kSKmNuEBarNC1pi0    ) : return "nueb_nc1pi0";   break;
     case ( kSKmNuEBarNCcoh     ) : return "nueb_nccoh";    break;
     case ( kSKmNuEBarNCoth_e   ) : return "nueb_ncoth";    break;
     case ( kSKmNuEBarNCoth_mu  ) : return "nueb_ncoth";    break;
     case ( kSKmSigNuECCQE      ) : return "signue_ccqe";     break;
     case ( kSKmSigNuECC1pi     ) : return "signue_ccpi";     break;
     case ( kSKmSigNuECCcoh     ) : return "signue_cccoh";    break;
     case ( kSKmSigNuECCoth_mu  ) : return "signue_ccoth";    break;
     case ( kSKmSigNuECCoth_e   ) : return "signue_ccoth";    break;
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
t2kana::conv::NuMassHierarchy_t  
  t2kana::util::HierarchyFromString(string sh) 
{
  if(sh == "normal"   ) return kNmhNormal;
  if(sh == "inverted" ) return kNmhInverted;
  
  return kNmhUndefined;
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
  if(sh == "mixing23") return kPhypMixing23;
  if(sh == "mixing13") return kPhypMixing13;
  if(sh == "mixing12") return kPhypMixing12;
  if(sh == "deltacp") return kPhypDeltaCP;
  if(sh == "dm23sq") return kPhypDm23sq;
  if(sh == "dm12sq") return kPhypDm12sq;
  if(sh == "alpha_decay") return kPhypAlphaDecay;
  if(sh == "sinsq2theta_decay") return kPhypSinsq2ThetaDecay;
  if(sh == "musq_decoh") return kPhypMusqDecoh;
  if(sh == "sinsq2theta_decoh") return kPhypSinsq2ThetaDecoh;

  if(sh == "sk_energy_scale") return kNuisSuperkEnergyScale;
  if(sh == "sk_1Rmu_ccqe_other_eff") return kNuisSuperkEff1RmuCCQEoth;
  if(sh == "sk_1Rmu_ccqe_ring_counting_eff_0") return kNuisSuperkEff1RmuCCQErc0u;
  if(sh == "sk_1Rmu_ccqe_ring_counting_eff_1") return kNuisSuperkEff1RmuCCQErc1u;
  if(sh == "sk_1Rmu_ccqe_ring_counting_eff_2") return kNuisSuperkEff1RmuCCQErc2u;
  if(sh == "sk_1Rmu_ccnonqe_eff") return kNuisSuperkEff1RmuCCnonQE;
  if(sh == "sk_1Rmu_nc_eff") return kNuisSuperkEff1RmuNC;
  if(sh == "sk_1Rmu_nuecc_eff") return kNuisSuperkEff1RmuNueCC;

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
#if BANFF_BINS <= 2
  if(sh == "banff_16") return kNuisBANFF16;
  if(sh == "banff_17") return kNuisBANFF17;
  if(sh == "banff_18") return kNuisBANFF18;
  if(sh == "banff_19") return kNuisBANFF19;
  if(sh == "banff_20") return kNuisBANFF20;
  if(sh == "banff_21") return kNuisBANFF21;
  if(sh == "banff_22") return kNuisBANFF22;
  if(sh == "banff_23") return kNuisBANFF23;
  if(sh == "banff_24") return kNuisBANFF24;
  if(sh == "banff_25") return kNuisBANFF25;
  if(sh == "banff_26") return kNuisBANFF26;
  if(sh == "banff_27") return kNuisBANFF27;
  if(sh == "banff_28") return kNuisBANFF28;
  if(sh == "banff_29") return kNuisBANFF29;
  if(sh == "banff_30") return kNuisBANFF30;
#if BANFF_BINS == -1
  if(sh == "banff_31") return kNuisBANFF31;
  if(sh == "banff_32") return kNuisBANFF32;
  if(sh == "banff_33") return kNuisBANFF33;
  if(sh == "banff_34") return kNuisBANFF34;
  if(sh == "banff_35") return kNuisBANFF35;
  if(sh == "banff_36") return kNuisBANFF36;
  if(sh == "banff_37") return kNuisBANFF37;
  if(sh == "banff_38") return kNuisBANFF38;
  if(sh == "banff_39") return kNuisBANFF39;
  if(sh == "banff_40") return kNuisBANFF40;
  if(sh == "banff_41") return kNuisBANFF41;
  if(sh == "banff_42") return kNuisBANFF42;
  if(sh == "banff_43") return kNuisBANFF43;
  if(sh == "banff_44") return kNuisBANFF44;
  if(sh == "banff_45") return kNuisBANFF45;
  if(sh == "banff_46") return kNuisBANFF46;
  if(sh == "banff_47") return kNuisBANFF47;
  if(sh == "banff_48") return kNuisBANFF48;
#endif
#endif
#if BANFF_BINS >= 3
  if(sh == "banff_16")  return kNuisBANFFMAQE;
  if(sh == "banff_17")  return kNuisBANFFMARES;
  if(sh == "banff_18")  return kNuisBANFFCCQEE0;
  if(sh == "banff_19")  return kNuisBANFFCCQEE1;
  if(sh == "banff_20")  return kNuisBANFFCCQEE2;
  if(sh == "banff_21")  return kNuisBANFFCC1piE0;
  if(sh == "banff_22")  return kNuisBANFFCC1piE1;
#endif
  if(sh == "banff_maqe")    return kNuisBANFFMAQE;
  if(sh == "banff_mares")   return kNuisBANFFMARES;
  if(sh == "banff_ccqeE0")  return kNuisBANFFCCQEE0;
#if BANFF_BINS >= 2
  if(sh == "banff_ccqeE1")  return kNuisBANFFCCQEE1;
  if(sh == "banff_ccqeE2")  return kNuisBANFFCCQEE2;
#endif
  if(sh == "banff_cc1piE0") return kNuisBANFFCC1piE0;
#if BANFF_BINS >= 2
  if(sh == "banff_cc1piE1") return kNuisBANFFCC1piE1;
#endif
#if BANFF_BINS <= 2
  if(sh == "banff_nc1pi0")  return kNuisBANFFNC1pi0;
#endif
#if BANFF_BINS == 5
  if(sh == "banff_mec")     return kNuisBANFFMEC;
#endif

  if(sh == "xsec_pf")           return kNuisXSecFermiMom;
  if(sh == "xsec_sf")           return kNuisXSecSpecFunc;
  if(sh == "xsec_wshape")       return kNuisXSecWshape;
  if(sh == "xsec_shape_ccoth")  return kNuisXSecShapeCCoth;
  if(sh == "xsec_pilessdelta")  return kNuisXSecPiLessDelta;
  if(sh == "xsec_enushape")     return kNuisXSecENuShape;
  if(sh == "xsec_bindE")        return kNuisXSecBindE;
#if BANFF_BINS <= 1
  if(sh == "xsec_norm_ccqeE1")  return kNuisXSecNormCCQEE1;
  if(sh == "xsec_norm_ccqeE2")  return kNuisXSecNormCCQEE2;
  if(sh == "xsec_norm_cc1piE1") return kNuisXSecNormCC1piE1;
#endif
  if(sh == "xsec_norm_cccoh")   return kNuisXSecNormCCcoh;
  if(sh == "xsec_norm_nc1piC")  return kNuisXSecNormNC1piC;
  if(sh == "xsec_norm_nccoh")   return kNuisXSecNormNCcoh;
  if(sh == "xsec_norm_ncoth")   return kNuisXSecNormNCoth;
  if(sh == "xsec_norm_nuetonumu") return kNuisXSecNormNuEToNuMu;
  if(sh == "xsec_norm_antinutonu") return kNuisXSecNormAntiNuToNu;

  if(sh == "fsi_000") return kNuisFSI000;
  if(sh == "fsi_001") return kNuisFSI001;
  if(sh == "fsi_002") return kNuisFSI002;
  if(sh == "fsi_003") return kNuisFSI003;
  if(sh == "fsi_004") return kNuisFSI004;
  if(sh == "fsi_005") return kNuisFSI005;
#if FSI_BINS >= 2 
  if(sh == "fsi_006") return kNuisFSI006;
  if(sh == "fsi_007") return kNuisFSI007;
  if(sh == "fsi_008") return kNuisFSI008;
  if(sh == "fsi_009") return kNuisFSI009;
  if(sh == "fsi_010") return kNuisFSI010;
  if(sh == "fsi_011") return kNuisFSI011;
  if(sh == "fsi_012") return kNuisFSI012;
  if(sh == "fsi_013") return kNuisFSI013;
  if(sh == "fsi_014") return kNuisFSI014;
#if FSI_BINS >= 3
  if(sh == "fsi_015") return kNuisFSI015;
  if(sh == "fsi_016") return kNuisFSI016;
  if(sh == "fsi_017") return kNuisFSI017;
  if(sh == "fsi_018") return kNuisFSI018;
  if(sh == "fsi_019") return kNuisFSI019;
  if(sh == "fsi_020") return kNuisFSI020;
  if(sh == "fsi_021") return kNuisFSI021;
  if(sh == "fsi_022") return kNuisFSI022;
  if(sh == "fsi_023") return kNuisFSI023;
  if(sh == "fsi_024") return kNuisFSI024;
  if(sh == "fsi_025") return kNuisFSI025;
  if(sh == "fsi_026") return kNuisFSI026;
  if(sh == "fsi_027") return kNuisFSI027;
  if(sh == "fsi_028") return kNuisFSI028;
  if(sh == "fsi_029") return kNuisFSI029;
  if(sh == "fsi_030") return kNuisFSI030;
  if(sh == "fsi_031") return kNuisFSI031;
  if(sh == "fsi_032") return kNuisFSI032;
  if(sh == "fsi_033") return kNuisFSI033;
  if(sh == "fsi_034") return kNuisFSI034;
  if(sh == "fsi_035") return kNuisFSI035;
  if(sh == "fsi_036") return kNuisFSI036;
  if(sh == "fsi_037") return kNuisFSI037;
  if(sh == "fsi_038") return kNuisFSI038;
  if(sh == "fsi_039") return kNuisFSI039;
  if(sh == "fsi_040") return kNuisFSI040;
  if(sh == "fsi_041") return kNuisFSI041;
  if(sh == "fsi_042") return kNuisFSI042;
  if(sh == "fsi_043") return kNuisFSI043;
  if(sh == "fsi_044") return kNuisFSI044;
  if(sh == "fsi_045") return kNuisFSI045;
  if(sh == "fsi_046") return kNuisFSI046;
  if(sh == "fsi_047") return kNuisFSI047;
  if(sh == "fsi_048") return kNuisFSI048;
#endif
#endif

  if(sh == "unbanff_00") return kNuisUnBANFF00;
  if(sh == "unbanff_01") return kNuisUnBANFF01;
  if(sh == "unbanff_02") return kNuisUnBANFF02;
  if(sh == "unbanff_03") return kNuisUnBANFF03;
  if(sh == "unbanff_04") return kNuisUnBANFF04;
  if(sh == "unbanff_05") return kNuisUnBANFF05;
  if(sh == "unbanff_06") return kNuisUnBANFF06;
  if(sh == "unbanff_07") return kNuisUnBANFF07;
  if(sh == "unbanff_08") return kNuisUnBANFF08;
  if(sh == "unbanff_09") return kNuisUnBANFF09;
  if(sh == "unbanff_10") return kNuisUnBANFF10;
  if(sh == "unbanff_11") return kNuisUnBANFF11;
  if(sh == "unbanff_12") return kNuisUnBANFF12;
  if(sh == "unbanff_13") return kNuisUnBANFF13;
  if(sh == "unbanff_14") return kNuisUnBANFF14;
  if(sh == "unbanff_15") return kNuisUnBANFF15;
  if(sh == "unbanff_16") return kNuisUnBANFF16;
  if(sh == "unbanff_17") return kNuisUnBANFF17;
  if(sh == "unbanff_18") return kNuisUnBANFF18;
  if(sh == "unbanff_19") return kNuisUnBANFF19;
  if(sh == "unbanff_20") return kNuisUnBANFF20;
  if(sh == "unbanff_21") return kNuisUnBANFF21;
  if(sh == "unbanff_22") return kNuisUnBANFF22;
#if BANFF_BINS <= 2
  if(sh == "unbanff_23") return kNuisUnBANFF23;
  if(sh == "unbanff_24") return kNuisUnBANFF24;
  if(sh == "unbanff_25") return kNuisUnBANFF25;
  if(sh == "unbanff_26") return kNuisUnBANFF26;
  if(sh == "unbanff_27") return kNuisUnBANFF27;
  if(sh == "unbanff_28") return kNuisUnBANFF28;
  if(sh == "unbanff_29") return kNuisUnBANFF29;
  if(sh == "unbanff_30") return kNuisUnBANFF30;
  if(sh == "unbanff_31") return kNuisUnBANFF31;
  if(sh == "unbanff_32") return kNuisUnBANFF32;
  if(sh == "unbanff_33") return kNuisUnBANFF33;
  if(sh == "unbanff_34") return kNuisUnBANFF34;
  if(sh == "unbanff_35") return kNuisUnBANFF35;
#if BANFF_BINS != 1
  if(sh == "unbanff_36") return kNuisUnBANFF36;
  if(sh == "unbanff_37") return kNuisUnBANFF37;
  if(sh == "unbanff_38") return kNuisUnBANFF38;
#if BANFF_BINS == -1
  if(sh == "unbanff_39") return kNuisUnBANFF39;
  if(sh == "unbanff_40") return kNuisUnBANFF40;
  if(sh == "unbanff_41") return kNuisUnBANFF41;
  if(sh == "unbanff_42") return kNuisUnBANFF42;
  if(sh == "unbanff_43") return kNuisUnBANFF43;
  if(sh == "unbanff_44") return kNuisUnBANFF44;
  if(sh == "unbanff_45") return kNuisUnBANFF45;
  if(sh == "unbanff_46") return kNuisUnBANFF46;
  if(sh == "unbanff_47") return kNuisUnBANFF47;
  if(sh == "unbanff_48") return kNuisUnBANFF48;
  if(sh == "unbanff_49") return kNuisUnBANFF49;
  if(sh == "unbanff_50") return kNuisUnBANFF50;
  if(sh == "unbanff_51") return kNuisUnBANFF51;
  if(sh == "unbanff_52") return kNuisUnBANFF52;
  if(sh == "unbanff_53") return kNuisUnBANFF53;
#endif
#endif
#endif

  if(sh == "edbanff_00") return kNuisEDBANFF00;
  if(sh == "edbanff_01") return kNuisEDBANFF01;
  if(sh == "edbanff_02") return kNuisEDBANFF02;
  if(sh == "edbanff_03") return kNuisEDBANFF03;
  if(sh == "edbanff_04") return kNuisEDBANFF04;
  if(sh == "edbanff_05") return kNuisEDBANFF05;
  if(sh == "edbanff_06") return kNuisEDBANFF06;
  if(sh == "edbanff_07") return kNuisEDBANFF07;
  if(sh == "edbanff_08") return kNuisEDBANFF08;
  if(sh == "edbanff_09") return kNuisEDBANFF09;
  if(sh == "edbanff_10") return kNuisEDBANFF10;
  if(sh == "edbanff_11") return kNuisEDBANFF11;
  if(sh == "edbanff_12") return kNuisEDBANFF12;
  if(sh == "edbanff_13") return kNuisEDBANFF13;
  if(sh == "edbanff_14") return kNuisEDBANFF14;
  if(sh == "edbanff_15") return kNuisEDBANFF15;
  if(sh == "edbanff_16") return kNuisEDBANFF16;
  if(sh == "edbanff_17") return kNuisEDBANFF17;
  if(sh == "edbanff_18") return kNuisEDBANFF18;
  if(sh == "edbanff_19") return kNuisEDBANFF19;
  if(sh == "edbanff_20") return kNuisEDBANFF20;
  if(sh == "edbanff_21") return kNuisEDBANFF21;
  if(sh == "edbanff_22") return kNuisEDBANFF22;
#if BANFF_BINS <= 2
  if(sh == "edbanff_23") return kNuisEDBANFF23;
  if(sh == "edbanff_24") return kNuisEDBANFF24;
  if(sh == "edbanff_25") return kNuisEDBANFF25;
  if(sh == "edbanff_26") return kNuisEDBANFF26;
  if(sh == "edbanff_27") return kNuisEDBANFF27;
  if(sh == "edbanff_28") return kNuisEDBANFF28;
  if(sh == "edbanff_29") return kNuisEDBANFF29;
  if(sh == "edbanff_30") return kNuisEDBANFF30;
  if(sh == "edbanff_31") return kNuisEDBANFF31;
  if(sh == "edbanff_32") return kNuisEDBANFF32;
  if(sh == "edbanff_33") return kNuisEDBANFF33;
  if(sh == "edbanff_34") return kNuisEDBANFF34;
  if(sh == "edbanff_35") return kNuisEDBANFF35;
#if BANFF_BINS != 1
  if(sh == "edbanff_36") return kNuisEDBANFF36;
  if(sh == "edbanff_37") return kNuisEDBANFF37;
  if(sh == "edbanff_38") return kNuisEDBANFF38;
#if BANFF_BINS == -1
  if(sh == "edbanff_39") return kNuisEDBANFF39;
  if(sh == "edbanff_40") return kNuisEDBANFF40;
  if(sh == "edbanff_41") return kNuisEDBANFF41;
  if(sh == "edbanff_42") return kNuisEDBANFF42;
  if(sh == "edbanff_43") return kNuisEDBANFF43;
  if(sh == "edbanff_44") return kNuisEDBANFF44;
  if(sh == "edbanff_45") return kNuisEDBANFF45;
  if(sh == "edbanff_46") return kNuisEDBANFF46;
  if(sh == "edbanff_47") return kNuisEDBANFF47;
  if(sh == "edbanff_48") return kNuisEDBANFF48;
  if(sh == "edbanff_49") return kNuisEDBANFF49;
  if(sh == "edbanff_50") return kNuisEDBANFF50;
  if(sh == "edbanff_51") return kNuisEDBANFF51;
  if(sh == "edbanff_52") return kNuisEDBANFF52;
  if(sh == "edbanff_53") return kNuisEDBANFF53;
#endif
#endif
#endif

  cout << sh << " not yet supported in util::PhysParamIdFromString()" << endl;
  return kPhypUndefined;
}
//____________________________________________________________________________
t2kana::conv::PhysParamId_t t2kana::util::PhysParamIdFromLongString(string sh)
{
  if(sh == "sin^2(2*theta_23)"   ) return kPhypMixing23;
  if(sh == "sin^2(theta_23)"     ) return kPhypMixing23;
  if(sh == "sin^2(2*theta_13)"   ) return kPhypMixing13;
  if(sh == "sin^2(theta_13)"     ) return kPhypMixing13;
  if(sh == "sin^2(2*theta_12)"   ) return kPhypMixing12;
  if(sh == "sin^2(theta_12)"     ) return kPhypMixing12;
  if(sh == "delta_cp"            ) return kPhypDeltaCP;
  if(sh == "|Dm^2_32|"           ) return kPhypDm23sq;
  if(sh == "Dm^2_12"             ) return kPhypDm12sq;
  if(sh == "alpha_decay"         ) return kPhypAlphaDecay;
  if(sh == "sin^2(2*theta_decay)") return kPhypSinsq2ThetaDecay;
  if(sh == "mu^2_decoh"          ) return kPhypMusqDecoh;
  if(sh == "sin^2(2*theta_decoh)") return kPhypSinsq2ThetaDecoh;

  if(sh == "SK energy scale"                ) return kNuisSuperkEnergyScale;
  if(sh == "SK efficiency; 1Rmu; CCQE/Other") return kNuisSuperkEff1RmuCCQEoth;
  if(sh == "SK efficiency; 1Rmu; CCQE/RC0"  ) return kNuisSuperkEff1RmuCCQErc0u;
  if(sh == "SK efficiency; 1Rmu; CCQE/RC1"  ) return kNuisSuperkEff1RmuCCQErc1u;
  if(sh == "SK efficiency; 1Rmu; CCQE/RC2"  ) return kNuisSuperkEff1RmuCCQErc2u;
  if(sh == "SK efficiency; 1Rmu; CCnonQE"   ) return kNuisSuperkEff1RmuCCnonQE;
  if(sh == "SK efficiency; 1Rmu; NC"        ) return kNuisSuperkEff1RmuNC;
  if(sh == "SK efficiency; 1Rmu; nue CC"    ) return kNuisSuperkEff1RmuNueCC;

  for(int ip = 0; ip < kNumOfFluxParamsInBANFFfit; ip++) {
    string comp = Form("BANFF %2d; SK %s flux; %.1f - %.1f GeV",
		ip, util::Pdg2String(kNuPdgCodeOfFluxParamInBANFFfit[ip]).c_str(), 
		kEtrueRangeOfFluxParamInBANFFfit[ip][0],
		kEtrueRangeOfFluxParamInBANFFfit[ip][1]);
    if(sh == comp) return (PhysParamId_t)((int)kNuisBANFF00+ip);
  }

#if BANFF_BINS == 4
  if(sh == "BANFF; 1 / (Ma QE)^2"               ) return kNuisBANFFMAQE;
#else
  if(sh == "BANFF; Ma QE"                       ) return kNuisBANFFMAQE;
#endif
  if(sh == "BANFF; Ma RES"                      ) return kNuisBANFFMARES;
  if(sh == "BANFF; Norm; CC QE; 0.0 - 1.5 GeV"  ) return kNuisBANFFCCQEE0;
#if BANFF_BINS >= 2
  if(sh == "BANFF; Norm; CC QE; 1.5 - 3.5 GeV"  ) return kNuisBANFFCCQEE1;
  if(sh == "BANFF; Norm; CC QE; 3.5 - 30.0 GeV" ) return kNuisBANFFCCQEE2;
#endif
  if(sh == "BANFF; Norm; CC 1pi; 0.0 - 2.5 GeV" ) return kNuisBANFFCC1piE0;
#if BANFF_BINS >= 2
  if(sh == "BANFF; Norm; CC 1pi; 2.5 - 30.0 GeV") return kNuisBANFFCC1piE1;
#endif
#if BANFF_BINS <= 2
  if(sh == "BANFF; Norm; NC 1pi0"               ) return kNuisBANFFNC1pi0;
#endif
#if BANFF_BINS == 5
  if(sh == "BANFF; MEC"                         ) return kNuisBANFFMEC;
#endif

    //
    // XSec params not correlated with flux
    //

 if(sh == "Fermi Momentum"     ) return kNuisXSecFermiMom;
 if(sh == "Spectral Func"      ) return kNuisXSecSpecFunc;
 if(sh == "W-shape"            ) return kNuisXSecWshape;
 if(sh == "Shape; CC Oth"      ) return kNuisXSecShapeCCoth;
 if(sh == "Pi-less delta decay") return kNuisXSecPiLessDelta;
 if(sh == "Shape; CC1pi ENu"   ) return kNuisXSecENuShape;
 if(sh == "Binding energy"     ) return kNuisXSecBindE;

#if BANFF_BINS <= 1
 if(sh == "Norm; CC QE Norm; 1.5 - 3.5 GeV") return kNuisXSecNormCCQEE1;
 if(sh == "Norm; CC QE Norm; > 3.5 GeV"    ) return kNuisXSecNormCCQEE2;
 if(sh == "Norm; CC 1pi; > 2.5 GeV"        ) return kNuisXSecNormCC1piE1;
#endif
 if(sh == "Norm; CC Coh"      ) return kNuisXSecNormCCcoh;
 if(sh == "Norm; NC 1piC"     ) return kNuisXSecNormNC1piC;
 if(sh == "Norm; NC Coh"      ) return kNuisXSecNormNCcoh;
 if(sh == "Norm; NC Oth"      ) return kNuisXSecNormNCoth;
 if(sh == "Norm; NuE To NuMu" ) return kNuisXSecNormNuEToNuMu;
 if(sh == "Norm; AntiNu To Nu") return kNuisXSecNormAntiNuToNu;

 for(int ip = 0; ip < kNumOfFSIParams; ip++) {
   string nutypes = "";
#if FSI_BINS == 2
   if(ip < 3)
	  nutypes = "numu";
	else if(ip < 6)
	  nutypes = "numubar";
	else if(ip < 9)
	  nutypes = "nue";
	else if(ip < 12)
	  nutypes = "nuebar";
	else if(ip < 15)
	  nutypes = "oscillated nue";
#elif FSI_BINS == 1
	if(ip < 3)
	  nutypes = "numu/numubar CCQE";
	else if(ip < 4)
	  nutypes = "numu/numubar CCnQE";
	else if(ip < 5)
	  nutypes = "nue/nuebar/signue CC";
        else if(ip < 6)
	  nutypes = "all NC";
#elif FSI_BINS == 3
	if(ip < 1)
	  nutypes = "numu/numubar CCQE";
	else if(ip < 17)
	  nutypes = "numu/numubar CCnQE";
	else if(ip < 33)
	  nutypes = "nue/nuebar/signue CC";
	else if(ip < 49)
	  nutypes = "all NC";
#endif
	string comp = Form("FSI/SI %2d; Ereco range %.1f - %.1f GeV; %s",ip,
			   kFSIParamsErecoCut[ip][0],kFSIParamsErecoCut[ip][1],nutypes.c_str());     
	if(sh == comp) return (PhysParamId_t)((int)kNuisFSI000+ip);
 }

    if(sh == "Uncorrelated BANFF 00") return kNuisUnBANFF00;
    if(sh == "Uncorrelated BANFF 01") return kNuisUnBANFF01;
    if(sh == "Uncorrelated BANFF 02") return kNuisUnBANFF02;
    if(sh == "Uncorrelated BANFF 03") return kNuisUnBANFF03;
    if(sh == "Uncorrelated BANFF 04") return kNuisUnBANFF04;
    if(sh == "Uncorrelated BANFF 05") return kNuisUnBANFF05;
    if(sh == "Uncorrelated BANFF 06") return kNuisUnBANFF06;
    if(sh == "Uncorrelated BANFF 07") return kNuisUnBANFF07;
    if(sh == "Uncorrelated BANFF 08") return kNuisUnBANFF08;
    if(sh == "Uncorrelated BANFF 09") return kNuisUnBANFF09;
    if(sh == "Uncorrelated BANFF 10") return kNuisUnBANFF10;
    if(sh == "Uncorrelated BANFF 11") return kNuisUnBANFF11;
    if(sh == "Uncorrelated BANFF 12") return kNuisUnBANFF12;
    if(sh == "Uncorrelated BANFF 13") return kNuisUnBANFF13;
    if(sh == "Uncorrelated BANFF 14") return kNuisUnBANFF14;
    if(sh == "Uncorrelated BANFF 15") return kNuisUnBANFF15;
    if(sh == "Uncorrelated BANFF 16") return kNuisUnBANFF16;
    if(sh == "Uncorrelated BANFF 17") return kNuisUnBANFF17;
    if(sh == "Uncorrelated BANFF 18") return kNuisUnBANFF18;
    if(sh == "Uncorrelated BANFF 19") return kNuisUnBANFF19; 
    if(sh == "Uncorrelated BANFF 20") return kNuisUnBANFF20; 
    if(sh == "Uncorrelated BANFF 21") return kNuisUnBANFF21; 
    if(sh == "Uncorrelated BANFF 22") return kNuisUnBANFF22;
#if BANFF_BINS <= 2
    if(sh == "Uncorrelated BANFF 23") return kNuisUnBANFF23;
    if(sh == "Uncorrelated BANFF 24") return kNuisUnBANFF24;
    if(sh == "Uncorrelated BANFF 25") return kNuisUnBANFF25;
    if(sh == "Uncorrelated BANFF 26") return kNuisUnBANFF26;
    if(sh == "Uncorrelated BANFF 27") return kNuisUnBANFF27;
    if(sh == "Uncorrelated BANFF 28") return kNuisUnBANFF28;
    if(sh == "Uncorrelated BANFF 29") return kNuisUnBANFF29;
    if(sh == "Uncorrelated BANFF 30") return kNuisUnBANFF30;
    if(sh == "Uncorrelated BANFF 31") return kNuisUnBANFF31;
    if(sh == "Uncorrelated BANFF 32") return kNuisUnBANFF32;
    if(sh == "Uncorrelated BANFF 33") return kNuisUnBANFF33;
    if(sh == "Uncorrelated BANFF 34") return kNuisUnBANFF34;
    if(sh == "Uncorrelated BANFF 35") return kNuisUnBANFF35;
#if BANFF_BINS != 1
    if(sh == "Uncorrelated BANFF 36") return kNuisUnBANFF36;
    if(sh == "Uncorrelated BANFF 37") return kNuisUnBANFF37;
    if(sh == "Uncorrelated BANFF 38") return kNuisUnBANFF38;
#if BANFF_BINS == -1
    if(sh == "Uncorrelated BANFF 39") return kNuisUnBANFF39;
    if(sh == "Uncorrelated BANFF 40") return kNuisUnBANFF40;
    if(sh == "Uncorrelated BANFF 41") return kNuisUnBANFF41;
    if(sh == "Uncorrelated BANFF 42") return kNuisUnBANFF42;
    if(sh == "Uncorrelated BANFF 43") return kNuisUnBANFF43;
    if(sh == "Uncorrelated BANFF 44") return kNuisUnBANFF44;
    if(sh == "Uncorrelated BANFF 45") return kNuisUnBANFF45;
    if(sh == "Uncorrelated BANFF 46") return kNuisUnBANFF46;
    if(sh == "Uncorrelated BANFF 47") return kNuisUnBANFF47;
    if(sh == "Uncorrelated BANFF 48") return kNuisUnBANFF48;
    if(sh == "Uncorrelated BANFF 49") return kNuisUnBANFF49;
    if(sh == "Uncorrelated BANFF 50") return kNuisUnBANFF50;
    if(sh == "Uncorrelated BANFF 51") return kNuisUnBANFF51;
    if(sh == "Uncorrelated BANFF 52") return kNuisUnBANFF52;
    if(sh == "Uncorrelated BANFF 53") return kNuisUnBANFF53;
#endif
#endif
#endif

    if(sh == "EigenDecomposition BANFF 00") return kNuisEDBANFF00;
    if(sh == "EigenDecomposition BANFF 01") return kNuisEDBANFF01;
    if(sh == "EigenDecomposition BANFF 02") return kNuisEDBANFF02;
    if(sh == "EigenDecomposition BANFF 03") return kNuisEDBANFF03;
    if(sh == "EigenDecomposition BANFF 04") return kNuisEDBANFF04;
    if(sh == "EigenDecomposition BANFF 05") return kNuisEDBANFF05;
    if(sh == "EigenDecomposition BANFF 06") return kNuisEDBANFF06;
    if(sh == "EigenDecomposition BANFF 07") return kNuisEDBANFF07;
    if(sh == "EigenDecomposition BANFF 08") return kNuisEDBANFF08;
    if(sh == "EigenDecomposition BANFF 09") return kNuisEDBANFF09;
    if(sh == "EigenDecomposition BANFF 10") return kNuisEDBANFF10;
    if(sh == "EigenDecomposition BANFF 11") return kNuisEDBANFF11;
    if(sh == "EigenDecomposition BANFF 12") return kNuisEDBANFF12;
    if(sh == "EigenDecomposition BANFF 13") return kNuisEDBANFF13;
    if(sh == "EigenDecomposition BANFF 14") return kNuisEDBANFF14;
    if(sh == "EigenDecomposition BANFF 15") return kNuisEDBANFF15;
    if(sh == "EigenDecomposition BANFF 16") return kNuisEDBANFF16;
    if(sh == "EigenDecomposition BANFF 17") return kNuisEDBANFF17;
    if(sh == "EigenDecomposition BANFF 18") return kNuisEDBANFF18;
    if(sh == "EigenDecomposition BANFF 19") return kNuisEDBANFF19; 
    if(sh == "EigenDecomposition BANFF 20") return kNuisEDBANFF20; 
    if(sh == "EigenDecomposition BANFF 21") return kNuisEDBANFF21; 
    if(sh == "EigenDecomposition BANFF 22") return kNuisEDBANFF22;
#if BANFF_BINS <= 2
    if(sh == "EigenDecomposition BANFF 23") return kNuisEDBANFF23;
    if(sh == "EigenDecomposition BANFF 24") return kNuisEDBANFF24;
    if(sh == "EigenDecomposition BANFF 25") return kNuisEDBANFF25;
    if(sh == "EigenDecomposition BANFF 26") return kNuisEDBANFF26;
    if(sh == "EigenDecomposition BANFF 27") return kNuisEDBANFF27;
    if(sh == "EigenDecomposition BANFF 28") return kNuisEDBANFF28;
    if(sh == "EigenDecomposition BANFF 29") return kNuisEDBANFF29;
    if(sh == "EigenDecomposition BANFF 30") return kNuisEDBANFF30;
    if(sh == "EigenDecomposition BANFF 31") return kNuisEDBANFF31;
    if(sh == "EigenDecomposition BANFF 32") return kNuisEDBANFF32;
    if(sh == "EigenDecomposition BANFF 33") return kNuisEDBANFF33;
    if(sh == "EigenDecomposition BANFF 34") return kNuisEDBANFF34;
    if(sh == "EigenDecomposition BANFF 35") return kNuisEDBANFF35;
#if BANFF_BINS != 1
    if(sh == "EigenDecomposition BANFF 36") return kNuisEDBANFF36;
    if(sh == "EigenDecomposition BANFF 37") return kNuisEDBANFF37;
    if(sh == "EigenDecomposition BANFF 38") return kNuisEDBANFF38;
#if BANFF_BINS == -1
    if(sh == "EigenDecomposition BANFF 39") return kNuisEDBANFF39;
    if(sh == "EigenDecomposition BANFF 40") return kNuisEDBANFF40;
    if(sh == "EigenDecomposition BANFF 41") return kNuisEDBANFF41;
    if(sh == "EigenDecomposition BANFF 42") return kNuisEDBANFF42;
    if(sh == "EigenDecomposition BANFF 43") return kNuisEDBANFF43;
    if(sh == "EigenDecomposition BANFF 44") return kNuisEDBANFF44;
    if(sh == "EigenDecomposition BANFF 45") return kNuisEDBANFF45;
    if(sh == "EigenDecomposition BANFF 46") return kNuisEDBANFF46;
    if(sh == "EigenDecomposition BANFF 47") return kNuisEDBANFF47;
    if(sh == "EigenDecomposition BANFF 48") return kNuisEDBANFF48;
    if(sh == "EigenDecomposition BANFF 49") return kNuisEDBANFF49;
    if(sh == "EigenDecomposition BANFF 50") return kNuisEDBANFF50;
    if(sh == "EigenDecomposition BANFF 51") return kNuisEDBANFF51;
    if(sh == "EigenDecomposition BANFF 52") return kNuisEDBANFF52;
    if(sh == "EigenDecomposition BANFF 53") return kNuisEDBANFF53;
#endif
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
  if ( util::IsNuMuBar   (m) ) return kPdgNuMuBar;
  if ( util::IsNuE       (m) ) return kPdgNuE;
  if ( util::IsNuEBar    (m) ) return kPdgNuEBar;
  if ( util::IsSignalNuE (m) ) return kPdgNuMu; /// numu -> nue
   
  // NOTE: Beam NuEBars not handled

  return 0;
}
//____________________________________________________________________________
int t2kana::util::SuperkMode2SKMCPdg(t2kana::conv::SuperkMode_t m)
{
  if ( util::IsNuMu      (m) ) return kPdgNuMu;
  if ( util::IsNuMuBar   (m) ) return kPdgNuMuBar;
  if ( util::IsNuE       (m) ) return kPdgNuE;
  if ( util::IsNuEBar    (m) ) return kPdgNuEBar;
  if ( util::IsSignalNuE (m) ) return kPdgNuE;
   
  // NOTE: Beam NuEBars not handled

  return 0;
}
//____________________________________________________________________________
bool t2kana::util::IsNuMu(t2kana::conv::SuperkMode_t  m)
{
  switch(m) {
     case (kSKmNuMuCCQE    ):
     case (kSKmNuMuCC1pi   ):
     case (kSKmNuMuCCcoh   ):
     case (kSKmNuMuCCoth_mu):
     case (kSKmNuMuCCoth_e ):
     case (kSKmNuMuNC1piC  ):
     case (kSKmNuMuNC1pi0  ):
     case (kSKmNuMuNCcoh   ):
     case (kSKmNuMuNCoth_e ):
     case (kSKmNuMuNCoth_mu):
       return true;
       break;
     default :                
       return false;
       break;
  }
  return false;
}
//____________________________________________________________________________
bool t2kana::util::IsNuMuBar(t2kana::conv::SuperkMode_t m)
{
  switch(m) {
     case (kSKmNuMuBarCCQE    ):
     case (kSKmNuMuBarCC1pi   ):
     case (kSKmNuMuBarCCcoh   ):
     case (kSKmNuMuBarCCoth_mu):
     case (kSKmNuMuBarCCoth_e ):
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
  return false;
}
//____________________________________________________________________________
bool t2kana::util::IsNuE(t2kana::conv::SuperkMode_t m)
{
  switch(m) {
     case (kSKmNuECCQE    ):
     case (kSKmNuECC1pi   ):
     case (kSKmNuECCcoh   ):
     case (kSKmNuECCoth_mu):
     case (kSKmNuECCoth_e ):
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
  return false;
}
//____________________________________________________________________________
bool t2kana::util::IsNuEBar(t2kana::conv::SuperkMode_t m)
{
  switch(m) {
     case (kSKmNuEBarCCQE    ):
     case (kSKmNuEBarCC1pi   ):
     case (kSKmNuEBarCCcoh   ):
     case (kSKmNuEBarCCoth_mu):
     case (kSKmNuEBarCCoth_e ):
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
  return false;
}
//____________________________________________________________________________
bool t2kana::util::IsSignalNuE(t2kana::conv::SuperkMode_t m)
{
  switch(m) {
     case ( kSKmSigNuECCQE     ) :
     case ( kSKmSigNuECC1pi    ) :
     case ( kSKmSigNuECCcoh    ) :
     case ( kSKmSigNuECCoth_e  ) :
     case ( kSKmSigNuECCoth_mu ) :
       return true;
       break;
     default :                
       return false;
       break;
  }
  return false;
}
//____________________________________________________________________________
bool t2kana::util::IsCC(t2kana::conv::SuperkMode_t m)
{
  switch(m) {
     case ( kSKmNuMuCCQE        ) : 
     case ( kSKmNuMuCC1pi       ) : 
     case ( kSKmNuMuCCcoh       ) : 
     case ( kSKmNuMuCCoth_mu    ) : 
     case ( kSKmNuMuCCoth_e     ) : 
     case ( kSKmNuMuBarCCQE     ) : 
     case ( kSKmNuMuBarCC1pi    ) : 
     case ( kSKmNuMuBarCCcoh    ) : 
     case ( kSKmNuMuBarCCoth_mu ) : 
     case ( kSKmNuMuBarCCoth_e  ) : 
     case ( kSKmNuECCQE         ) : 
     case ( kSKmNuECC1pi        ) : 
     case ( kSKmNuECCcoh        ) : 
     case ( kSKmNuECCoth_mu     ) : 
     case ( kSKmNuECCoth_e      ) : 
     case ( kSKmNuEBarCCQE      ) : 
     case ( kSKmNuEBarCC1pi     ) : 
     case ( kSKmNuEBarCCcoh     ) : 
     case ( kSKmNuEBarCCoth_mu  ) : 
     case ( kSKmNuEBarCCoth_e   ) : 
     case ( kSKmSigNuECCQE      ) : 
     case ( kSKmSigNuECC1pi     ) : 
     case ( kSKmSigNuECCcoh     ) : 
     case ( kSKmSigNuECCoth_mu  ) : 
     case ( kSKmSigNuECCoth_e   ) : 
	return true;
	break;
     default: 
	return false;
	break;
  }
  return false;
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
bool t2kana::util::IsSystematic(PhysParamId_t p)
{
  switch(p) {
    case ( kPhypMixing23              ) : return false; break;
    case ( kPhypMixing13              ) : return false; break;
    case ( kPhypMixing12              ) : return false; break;
    case ( kPhypDeltaCP               ) : return false; break;
    case ( kPhypDm23sq                ) : return false; break;
    case ( kPhypDm12sq                ) : return false; break;
    case ( kPhypAlphaDecay            ) : return false; break;
    case ( kPhypSinsq2ThetaDecay      ) : return false; break;
    case ( kPhypMusqDecoh             ) : return false; break;
    case ( kPhypSinsq2ThetaDecoh      ) : return false; break;
    case ( kNuisSuperkEnergyScale     ) : return true;  break;
    case ( kNuisSuperkEff1RmuCCQEoth  ) : return true;  break;
    case ( kNuisSuperkEff1RmuCCQErc0u ) : return true;  break;
    case ( kNuisSuperkEff1RmuCCQErc1u ) : return true;  break;
    case ( kNuisSuperkEff1RmuCCQErc2u ) : return true;  break;
    case ( kNuisSuperkEff1RmuCCnonQE  ) : return true;  break;
    case ( kNuisSuperkEff1RmuNC       ) : return true;  break;
    case ( kNuisSuperkEff1RmuNueCC    ) : return true;  break;
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
#if BANFF_BINS <= 2
    case ( kNuisBANFF16               ) : return true;  break;
    case ( kNuisBANFF17               ) : return true;  break;
    case ( kNuisBANFF18               ) : return true;  break;
    case ( kNuisBANFF19               ) : return true;  break;
    case ( kNuisBANFF20               ) : return true;  break;
    case ( kNuisBANFF21               ) : return true;  break;
    case ( kNuisBANFF22               ) : return true;  break;
    case ( kNuisBANFF23               ) : return true;  break;
    case ( kNuisBANFF24               ) : return true;  break;
    case ( kNuisBANFF25               ) : return true;  break;
    case ( kNuisBANFF26               ) : return true;  break;
    case ( kNuisBANFF27               ) : return true;  break;
    case ( kNuisBANFF28               ) : return true;  break;
    case ( kNuisBANFF29               ) : return true;  break;
    case ( kNuisBANFF30               ) : return true;  break;
#if BANFF_BINS == -1
    case ( kNuisBANFF31               ) : return true;  break;
    case ( kNuisBANFF32               ) : return true;  break;
    case ( kNuisBANFF33               ) : return true;  break;
    case ( kNuisBANFF34               ) : return true;  break;
    case ( kNuisBANFF35               ) : return true;  break;
    case ( kNuisBANFF36               ) : return true;  break;
    case ( kNuisBANFF37               ) : return true;  break;
    case ( kNuisBANFF38               ) : return true;  break;
    case ( kNuisBANFF39               ) : return true;  break;
    case ( kNuisBANFF40               ) : return true;  break;
    case ( kNuisBANFF41               ) : return true;  break;
    case ( kNuisBANFF42               ) : return true;  break;
    case ( kNuisBANFF43               ) : return true;  break;
    case ( kNuisBANFF44               ) : return true;  break;
    case ( kNuisBANFF45               ) : return true;  break;
    case ( kNuisBANFF46               ) : return true;  break;
    case ( kNuisBANFF47               ) : return true;  break;
    case ( kNuisBANFF48               ) : return true;  break;
#endif
#endif
    case ( kNuisBANFFMAQE             ) : return true;  break;
    case ( kNuisBANFFMARES            ) : return true;  break;
    case ( kNuisBANFFCCQEE0           ) : return true;  break;
#if BANFF_BINS >= 2
    case ( kNuisBANFFCCQEE1           ) : return true;  break;
    case ( kNuisBANFFCCQEE2           ) : return true;  break;
#endif
    case ( kNuisBANFFCC1piE0          ) : return true;  break;
#if BANFF_BINS >= 2
    case ( kNuisBANFFCC1piE1          ) : return true;  break;
#endif
#if BANFF_BINS <= 2
    case ( kNuisBANFFNC1pi0           ) : return true;  break;
#endif
#if BANFF_BINS == 5
    case ( kNuisBANFFMEC              ) : return true;  break;
#endif

    case ( kNuisXSecFermiMom          ) : return true;  break;
    case ( kNuisXSecSpecFunc          ) : return true;  break; 
    case ( kNuisXSecWshape            ) : return true;  break;
    case ( kNuisXSecShapeCCoth        ) : return true;  break;
    case ( kNuisXSecPiLessDelta       ) : return true;  break;
    case ( kNuisXSecENuShape          ) : return true;  break;
    case ( kNuisXSecBindE             ) : return true;  break;

#if BANFF_BINS <= 1
    case ( kNuisXSecNormCCQEE1        ) : return true;  break; 
    case ( kNuisXSecNormCCQEE2        ) : return true;  break; 
    case ( kNuisXSecNormCC1piE1       ) : return true;  break; 
#endif
    case ( kNuisXSecNormCCcoh         ) : return true;  break;
    case ( kNuisXSecNormNC1piC        ) : return true;  break;
    case ( kNuisXSecNormNCcoh         ) : return true;  break;
    case ( kNuisXSecNormNCoth         ) : return true;  break;
    case ( kNuisXSecNormNuEToNuMu     ) : return true;  break;
    case ( kNuisXSecNormAntiNuToNu    ) : return true;  break;

    case ( kNuisFSI000                ) : return true;  break;
    case ( kNuisFSI001                ) : return true;  break;
    case ( kNuisFSI002                ) : return true;  break;
    case ( kNuisFSI003                ) : return true;  break;
    case ( kNuisFSI004                ) : return true;  break;
    case ( kNuisFSI005                ) : return true;  break;
#if FSI_BINS >= 2
    case ( kNuisFSI006                ) : return true;  break;
    case ( kNuisFSI007                ) : return true;  break;
    case ( kNuisFSI008                ) : return true;  break;
    case ( kNuisFSI009                ) : return true;  break;
    case ( kNuisFSI010                ) : return true;  break;
    case ( kNuisFSI011                ) : return true;  break;
    case ( kNuisFSI012                ) : return true;  break;
    case ( kNuisFSI013                ) : return true;  break;
    case ( kNuisFSI014                ) : return true;  break;
#if FSI_BINS >= 3
    case ( kNuisFSI015                ) : return true;  break;
    case ( kNuisFSI016                ) : return true;  break;
    case ( kNuisFSI017                ) : return true;  break;
    case ( kNuisFSI018                ) : return true;  break;
    case ( kNuisFSI019                ) : return true;  break;
    case ( kNuisFSI020                ) : return true;  break;
    case ( kNuisFSI021                ) : return true;  break;
    case ( kNuisFSI022                ) : return true;  break;
    case ( kNuisFSI023                ) : return true;  break;
    case ( kNuisFSI024                ) : return true;  break;
    case ( kNuisFSI025                ) : return true;  break;
    case ( kNuisFSI026                ) : return true;  break;
    case ( kNuisFSI027                ) : return true;  break;
    case ( kNuisFSI028                ) : return true;  break;
    case ( kNuisFSI029                ) : return true;  break;
    case ( kNuisFSI030                ) : return true;  break;
    case ( kNuisFSI031                ) : return true;  break;
    case ( kNuisFSI032                ) : return true;  break;
    case ( kNuisFSI033                ) : return true;  break;
    case ( kNuisFSI034                ) : return true;  break;
    case ( kNuisFSI035                ) : return true;  break;
    case ( kNuisFSI036                ) : return true;  break;
    case ( kNuisFSI037                ) : return true;  break;
    case ( kNuisFSI038                ) : return true;  break;
    case ( kNuisFSI039                ) : return true;  break;
    case ( kNuisFSI040                ) : return true;  break;
    case ( kNuisFSI041                ) : return true;  break;
    case ( kNuisFSI042                ) : return true;  break;
    case ( kNuisFSI043                ) : return true;  break;
    case ( kNuisFSI044                ) : return true;  break;
    case ( kNuisFSI045                ) : return true;  break;
    case ( kNuisFSI046                ) : return true;  break;
    case ( kNuisFSI047                ) : return true;  break;
    case ( kNuisFSI048                ) : return true;  break;
#endif
#endif
    case ( kNuisUnBANFF00             ) : return true; break;
    case ( kNuisUnBANFF01             ) : return true; break;
    case ( kNuisUnBANFF02             ) : return true; break;
    case ( kNuisUnBANFF03             ) : return true; break;
    case ( kNuisUnBANFF04             ) : return true; break;
    case ( kNuisUnBANFF05             ) : return true; break;
    case ( kNuisUnBANFF06             ) : return true; break;
    case ( kNuisUnBANFF07             ) : return true; break;
    case ( kNuisUnBANFF08             ) : return true; break;
    case ( kNuisUnBANFF09             ) : return true; break;
    case ( kNuisUnBANFF10             ) : return true; break;
    case ( kNuisUnBANFF11             ) : return true; break;
    case ( kNuisUnBANFF12             ) : return true; break;
    case ( kNuisUnBANFF13             ) : return true; break;
    case ( kNuisUnBANFF14             ) : return true; break;
    case ( kNuisUnBANFF15             ) : return true; break;
    case ( kNuisUnBANFF16             ) : return true; break;
    case ( kNuisUnBANFF17             ) : return true; break;
    case ( kNuisUnBANFF18             ) : return true; break;
    case ( kNuisUnBANFF19             ) : return true; break;
    case ( kNuisUnBANFF20             ) : return true; break;
    case ( kNuisUnBANFF21             ) : return true; break;
    case ( kNuisUnBANFF22             ) : return true; break;
#if BANFF_BINS <= 2
    case ( kNuisUnBANFF23             ) : return true; break;
    case ( kNuisUnBANFF24             ) : return true; break;
    case ( kNuisUnBANFF25             ) : return true; break;
    case ( kNuisUnBANFF26             ) : return true; break;
    case ( kNuisUnBANFF27             ) : return true; break;
    case ( kNuisUnBANFF28             ) : return true; break;
    case ( kNuisUnBANFF29             ) : return true; break;
    case ( kNuisUnBANFF30             ) : return true; break;
    case ( kNuisUnBANFF31             ) : return true; break;
    case ( kNuisUnBANFF32             ) : return true; break;
    case ( kNuisUnBANFF33             ) : return true; break;
    case ( kNuisUnBANFF34             ) : return true; break;
    case ( kNuisUnBANFF35             ) : return true; break;
#if BANFF_BINS != 1
    case ( kNuisUnBANFF36             ) : return true; break;
    case ( kNuisUnBANFF37             ) : return true; break;
    case ( kNuisUnBANFF38             ) : return true; break;
#if BANFF_BINS == -1
    case ( kNuisUnBANFF39             ) : return true; break;
    case ( kNuisUnBANFF40             ) : return true; break;
    case ( kNuisUnBANFF41             ) : return true; break;
    case ( kNuisUnBANFF42             ) : return true; break;
    case ( kNuisUnBANFF43             ) : return true; break;
    case ( kNuisUnBANFF44             ) : return true; break;
    case ( kNuisUnBANFF45             ) : return true; break;
    case ( kNuisUnBANFF46             ) : return true; break;
    case ( kNuisUnBANFF47             ) : return true; break;
    case ( kNuisUnBANFF48             ) : return true; break;
    case ( kNuisUnBANFF49             ) : return true; break;
    case ( kNuisUnBANFF50             ) : return true; break;
    case ( kNuisUnBANFF51             ) : return true; break;
    case ( kNuisUnBANFF52             ) : return true; break;
    case ( kNuisUnBANFF53             ) : return true; break;
#endif
#endif
#endif
    case ( kNuisEDBANFF00             ) : return true; break;
    case ( kNuisEDBANFF01             ) : return true; break;
    case ( kNuisEDBANFF02             ) : return true; break;
    case ( kNuisEDBANFF03             ) : return true; break;
    case ( kNuisEDBANFF04             ) : return true; break;
    case ( kNuisEDBANFF05             ) : return true; break;
    case ( kNuisEDBANFF06             ) : return true; break;
    case ( kNuisEDBANFF07             ) : return true; break;
    case ( kNuisEDBANFF08             ) : return true; break;
    case ( kNuisEDBANFF09             ) : return true; break;
    case ( kNuisEDBANFF10             ) : return true; break;
    case ( kNuisEDBANFF11             ) : return true; break;
    case ( kNuisEDBANFF12             ) : return true; break;
    case ( kNuisEDBANFF13             ) : return true; break;
    case ( kNuisEDBANFF14             ) : return true; break;
    case ( kNuisEDBANFF15             ) : return true; break;
    case ( kNuisEDBANFF16             ) : return true; break;
    case ( kNuisEDBANFF17             ) : return true; break;
    case ( kNuisEDBANFF18             ) : return true; break;
    case ( kNuisEDBANFF19             ) : return true; break;
    case ( kNuisEDBANFF20             ) : return true; break;
    case ( kNuisEDBANFF21             ) : return true; break;
    case ( kNuisEDBANFF22             ) : return true; break;
#if BANFF_BINS <= 2
    case ( kNuisEDBANFF23             ) : return true; break;
    case ( kNuisEDBANFF24             ) : return true; break;
    case ( kNuisEDBANFF25             ) : return true; break;
    case ( kNuisEDBANFF26             ) : return true; break;
    case ( kNuisEDBANFF27             ) : return true; break;
    case ( kNuisEDBANFF28             ) : return true; break;
    case ( kNuisEDBANFF29             ) : return true; break;
    case ( kNuisEDBANFF30             ) : return true; break;
    case ( kNuisEDBANFF31             ) : return true; break;
    case ( kNuisEDBANFF32             ) : return true; break;
    case ( kNuisEDBANFF33             ) : return true; break;
    case ( kNuisEDBANFF34             ) : return true; break;
    case ( kNuisEDBANFF35             ) : return true; break;
#if BANFF_BINS != 1
    case ( kNuisEDBANFF36             ) : return true; break;
    case ( kNuisEDBANFF37             ) : return true; break;
    case ( kNuisEDBANFF38             ) : return true; break;
#if BANFF_BINS == -1
    case ( kNuisEDBANFF39             ) : return true; break;
    case ( kNuisEDBANFF40             ) : return true; break;
    case ( kNuisEDBANFF41             ) : return true; break;
    case ( kNuisEDBANFF42             ) : return true; break;
    case ( kNuisEDBANFF43             ) : return true; break;
    case ( kNuisEDBANFF44             ) : return true; break;
    case ( kNuisEDBANFF45             ) : return true; break;
    case ( kNuisEDBANFF46             ) : return true; break;
    case ( kNuisEDBANFF47             ) : return true; break;
    case ( kNuisEDBANFF48             ) : return true; break;
    case ( kNuisEDBANFF49             ) : return true; break;
    case ( kNuisEDBANFF50             ) : return true; break;
    case ( kNuisEDBANFF51             ) : return true; break;
    case ( kNuisEDBANFF52             ) : return true; break;
    case ( kNuisEDBANFF53             ) : return true; break;
#endif
#endif
#endif

    default : 
      return false;
  }
  return false;
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
const int t2kana::util::NumErecoBins(t2kana::conv::SuperkSample_t s)
{
  switch(s) {
     case ( kSKsFCFV1RingMuLike ) : return kNErecoBins1Rmu; break;
     case ( kSKsFCFV1RingELike  ) : return kNErecoBins1Re;  break;
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
     default                      : return 0;
  }
  return 0;
}
//___________________________________________________________________________
const int t2kana::util::NumVariableXsEtrueBins(t2kana::conv::SuperkSample_t s,
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
const int t2kana::util::NumVariableMCEtrueBins(t2kana::conv::SuperkSample_t s,
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
  bool printing_fsi_params   = opt.Contains("FSI params");
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
         PhysParamId_t p = PhysParamId_t ( (int)kNuisFSI000 + ip );
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
         PhysParamId_t p = PhysParamId_t ( (int)kNuisFSI000 + ip );
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
int t2kana::util::Sign(double tweaking_dial)
{
  if(TMath::Abs(tweaking_dial) < kASmallNum) return 0;
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

/*
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
}*/
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
bool t2kana::util::IsCCoth_mu(int neut_reaction_mode) 
{
  return IsCCoth(neut_reaction_mode);
}
//___________________________________________________________________________
bool t2kana::util::IsCCoth_e(int neut_reaction_mode) 
{
  return IsCCoth(neut_reaction_mode);
}
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
  if(m == kSKmNuMuCCoth_mu) 
  {
     bool is_numu     = (nu_pdg_code == kPdgNuMu);
     use_it = is_numu && IsCCoth_mu(neut_reaction_mode);
  } 

  else
  if(m == kSKmNuMuCCoth_e)
  {
    bool is_numu    = (nu_pdg_code == kPdgNuMu);
    use_it = is_numu && IsCCoth_e(neut_reaction_mode);
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
  if(m == kSKmNuMuBarCCoth_mu) 
  {
     bool is_numubar  = (nu_pdg_code == kPdgNuMuBar);
     use_it = is_numubar && IsCCoth_mu(-neut_reaction_mode);
  } 

  else
  if(m == kSKmNuMuBarCCoth_e) 
  {
     bool is_numubar = (nu_pdg_code == kPdgNuMuBar);
     use_it = is_numubar && IsCCoth_e(-neut_reaction_mode);
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
  if(m == kSKmNuECCQE || m == kSKmSigNuECCQE) 
  {
     bool is_nue  = (nu_pdg_code == kPdgNuE); 
     use_it = is_nue && IsCCQE(neut_reaction_mode);
  } 

  else
  if(m == kSKmNuECC1pi || m == kSKmSigNuECC1pi) 
  {
     bool is_nue   = (nu_pdg_code == kPdgNuE); 
     use_it = is_nue && IsCC1pi(neut_reaction_mode);
  }

  else
  if(m == kSKmNuECCcoh || m == kSKmSigNuECCcoh) 
  {
     bool is_nue   = (nu_pdg_code == kPdgNuE); 
     use_it = is_nue && IsCCcoh(neut_reaction_mode);
  } 

  else
  if(m == kSKmNuECCoth_mu || m == kSKmSigNuECCoth_mu) 
  {
     bool is_nue      = (nu_pdg_code == kPdgNuE); 
     use_it = is_nue && IsCCoth_mu(neut_reaction_mode);
  }

  else
  if(m == kSKmNuECCoth_e) 
  {
     bool is_nue     = (nu_pdg_code == kPdgNuE); 
     use_it = is_nue && IsCCoth_e(neut_reaction_mode);
  }

  //this is a little different. nuesig has no CCcoh in 1Re sample, whereas CCcoh is included
  // for the other species
  else
  if(m == kSKmSigNuECCoth_e) 
  {
     bool is_nue     = (nu_pdg_code == kPdgNuE); 
     use_it = is_nue && (IsCCoth_e(neut_reaction_mode) || IsCCcoh(neut_reaction_mode));
  }

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
  if(m == kSKmNuEBarCCoth_mu) 
  {
     bool is_nue      = (nu_pdg_code == kPdgNuEBar); 
     use_it = is_nue && IsCCoth_mu(-neut_reaction_mode);
  }

  else
  if(m == kSKmNuEBarCCoth_e) 
  {
     bool is_nue     = (nu_pdg_code == kPdgNuEBar); 
     use_it = is_nue && IsCCoth_e(-neut_reaction_mode);
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
  double mixing23, double dm23sq, double Et)
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

  psurv = TMath::Min(1.,psurv);
  psurv = TMath::Max(0.,psurv);

  return psurv;
}
//____________________________________________________________________________
double t2kana::util::NuMuSurvProbOsc2Flavour(
  double mixing23, double dm23sq, double Et, t2kana::conv::SuperkMode_t m)
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
     m == kSKmNuECCoth_mu     ||
     m == kSKmNuENC1piC       ||
     m == kSKmNuENC1pi0       ||
     m == kSKmNuENCcoh        ||
     m == kSKmNuENCoth_mu     ||
     m == kSKmNuENCoth_e      ||
     m == kSKmNuEBarCCQE      ||
     m == kSKmNuEBarCC1pi     ||
     m == kSKmNuEBarCCcoh     ||
     m == kSKmNuEBarCCoth_mu  ||
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
     m == kSKmNuMuCCoth_mu    ||
     m == kSKmNuMuCCoth_e     ||
     m == kSKmNuMuBarCCQE     ||
     m == kSKmNuMuBarCC1pi    ||
     m == kSKmNuMuBarCCcoh    ||
     m == kSKmNuMuBarCCoth_mu || 
     m == kSKmNuMuBarCCoth_e) 
  {
    return NuMuSurvProbOsc2Flavour(mixing23, dm23sq, Et);
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
     m == kSKmNuECCoth_mu    ||
     m == kSKmNuENC1piC      ||
     m == kSKmNuENCoth_mu    ||
     m == kSKmNuEBarCCQE     ||
     m == kSKmNuEBarCC1pi    ||
     m == kSKmNuEBarCCcoh    ||
     m == kSKmNuEBarCCoth_mu ||
     m == kSKmNuEBarNC1piC   ||
     m == kSKmNuEBarNCoth_mu )
  {
    return 1.;
  }
  else
  if(m == kSKmNuMuCCQE        || 
     m == kSKmNuMuCC1pi       || 
     m == kSKmNuMuCCcoh       || 
     m == kSKmNuMuCCoth_mu    ||
     m == kSKmNuMuNC1piC      ||
     m == kSKmNuMuNCoth_mu    ||
     m == kSKmNuMuBarCCQE     ||
     m == kSKmNuMuBarCC1pi    ||
     m == kSKmNuMuBarCCcoh    ||
     m == kSKmNuMuBarCCoth_mu ||
     m == kSKmNuMuBarNC1piC   ||
     m == kSKmNuMuBarNCoth_mu ) 
  {
    return NuMuSurvProbDecay(sinsq2theta, alpha, Et);
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
     m == kSKmNuECCoth_mu    ||
     m == kSKmNuENC1piC      ||
     m == kSKmNuENCoth_mu    ||
     m == kSKmNuEBarCCQE     ||
     m == kSKmNuEBarCC1pi    ||
     m == kSKmNuEBarCCcoh    ||
     m == kSKmNuEBarCCoth_mu ||
     m == kSKmNuEBarNC1piC   ||
     m == kSKmNuEBarNCoth_mu )
  {
    return 1.;
  }
  else
  if(m == kSKmNuMuCCQE        || 
     m == kSKmNuMuCC1pi       || 
     m == kSKmNuMuCCcoh       || 
     m == kSKmNuMuCCoth_mu    ||
     m == kSKmNuMuNC1piC      ||
     m == kSKmNuMuNCoth_mu    ||
     m == kSKmNuMuBarCCQE     ||
     m == kSKmNuMuBarCC1pi    ||
     m == kSKmNuMuBarCCcoh    ||
     m == kSKmNuMuBarCCoth_mu ||
     m == kSKmNuMuBarNC1piC   ||
     m == kSKmNuMuBarNCoth_mu ) 
  {
    return NuMuSurvProbDecoh(sinsq2theta, musq, Et);
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
  //UnBANFF
  for(int i = 0; i < knppmax; i++) {
    fUnBANFFSystIdx[i] = -1;
  }
  for(int ip = 0; ip < kNumOfParamsInBANFFfit; ip++) {
    PhysParamId_t p = PhysParamId_t ((int) kNuisUnBANFF00 + ip);
    if(p == kPhypUndefined) continue;
    fUnBANFFSystIdx[p] = ip;
#ifdef __show_arr_idx_mesg_debug__
     cout << "Util/ArrIdx >> UnBANFF: " 
          << util::AsString(p)
          << " is at slot: " << fUnBANFFSystIdx [p] << endl;
#endif
  }
  //EDBANFF
  for(int i = 0; i < knppmax; i++) {
    fEDBANFFSystIdx[i] = -1;
  }
  for(int ip = 0; ip < kNumOfParamsInBANFFfit; ip++) {
    PhysParamId_t p = PhysParamId_t ((int) kNuisEDBANFF00 + ip);
    if(p == kPhypUndefined) continue;
    fEDBANFFSystIdx[p] = ip;
#ifdef __show_arr_idx_mesg_debug__
     cout << "Util/ArrIdx >> EDBANFF: " 
          << util::AsString(p)
          << " is at slot: " << fEDBANFFSystIdx [p] << endl;
#endif
  }
  //FSI
  for(int i = 0; i < knppmax; i++) {
    fFSISystIdx[i] = -1;
  }
  for(int ip = 0; ip < kNumOfFSIParams; ip++) {
    PhysParamId_t p = PhysParamId_t ((int) kNuisFSI000 + ip);
    if(p == kPhypUndefined) continue;
    fFSISystIdx[p] = ip;
#ifdef __show_arr_idx_mesg_debug__
     cout << "Util/ArrIdx >> FSI: " 
          << util::AsString(p)
          << " is at slot: " << fFSISystIdx [p] << endl;
#endif
  }
  //PhysParam
  for(int i = 0; i < knppmax; i++) {
    fPhysParamIdx[i] = -1;
  }
  for(int ip = 0; ip < kNumOfFSIParams; ip++) {
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

