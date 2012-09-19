//____________________________________________________________________________
/*

 >> ChangeLog:

 For changes up to version 2010_v3.0.0 (last implementation of the 2010-era analysis)
 see code in the corresponding SVN branch.

 Changes since 2011_v0.9.0:

 @
   
*/
//____________________________________________________________________________

#define __show_cholesky_decomposition_mesg_debug__
#define __show_cholesky_gen_corr_param_variations_mesg_debug__
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

TRandom3 gR(89289281);

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
    {
      int ip = (int)p - (int)kNuisBANFF00;
      return Form("BANFF %2d; SK %s flux; %.1f - %.1f GeV",
                   ip, util::Pdg2String(kNuPdgCodeOfFluxParamInBANFFfit[ip]).c_str(), 
                   kEtrueRangeOfFluxParamInBANFFfit[ip][0],
                   kEtrueRangeOfFluxParamInBANFFfit[ip][1]);     
      break;
    }
    case ( kNuisBANFFMAQE    ) : return "BANFF; Ma QE";                           break;
    case ( kNuisBANFFMARES   ) : return "BANFF; Ma RES";                          break;
    case ( kNuisBANFFCCQEE0  ) : return "BANFF; Norm; CC QE; 0.0 - 1.5 GeV";      break;
#if BANFF_BINS == 2
    case ( kNuisBANFFCCQEE1  ) : return "BANFF; Norm; CC QE; 1.5 - 3.5 GeV";      break;
    case ( kNuisBANFFCCQEE2  ) : return "BANFF; Norm; CC QE; 3.5 - 30.0 GeV";     break;
#endif
    case ( kNuisBANFFCC1piE0 ) : return "BANFF; Norm; CC 1pi; 0.0 - 2.5 GeV";     break;
#if BANFF_BINS == 2
    case ( kNuisBANFFCC1piE1 ) : return "BANFF; Norm; CC 1pi; 2.5 - 30.0 GeV";    break;
#endif
    case ( kNuisBANFFNC1pi0  ) : return "BANFF; Norm; NC 1pi0";                   break;

    //
    // XSec params not correlated with flux
    //

    case ( kNuisXSecFermiMom   ) : return "Fermi Momentum";                  break;
    case ( kNuisXSecSpecFunc   ) : return "Spectral Func";                   break; 
    case ( kNuisXSecWshape     ) : return "W-shape";                         break;
    case ( kNuisXSecShapeCCoth ) : return "Shape; CC Oth";                   break;
    case ( kNuisXSecPiLessDelta) : return "Pi-less delta decay";             break;
    case ( kNuisXSecENuShape   ) : return "Shape; CC1pi ENu";                break;

#if BANFF_BINS != 2
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
    case ( kNuisFSI049 ) :
    case ( kNuisFSI050 ) :
    case ( kNuisFSI051 ) :
    case ( kNuisFSI052 ) :
    case ( kNuisFSI053 ) :
    case ( kNuisFSI054 ) :
    case ( kNuisFSI055 ) :
    case ( kNuisFSI056 ) :
    case ( kNuisFSI057 ) :
    case ( kNuisFSI058 ) :
    case ( kNuisFSI059 ) :
    case ( kNuisFSI060 ) :
    case ( kNuisFSI061 ) :
    case ( kNuisFSI062 ) :
    case ( kNuisFSI063 ) :
    case ( kNuisFSI064 ) :
    case ( kNuisFSI065 ) :
    case ( kNuisFSI066 ) :
    case ( kNuisFSI067 ) :
    case ( kNuisFSI068 ) :
    case ( kNuisFSI069 ) :
    case ( kNuisFSI070 ) :
    case ( kNuisFSI071 ) :
    case ( kNuisFSI072 ) :
    case ( kNuisFSI073 ) :
    case ( kNuisFSI074 ) :
    case ( kNuisFSI075 ) :
    case ( kNuisFSI076 ) :
    case ( kNuisFSI077 ) :
    case ( kNuisFSI078 ) :
#if FSI_BINS >= 4
    case ( kNuisFSI079 ) :
    case ( kNuisFSI080 ) :
    case ( kNuisFSI081 ) :
    case ( kNuisFSI082 ) :
    case ( kNuisFSI083 ) :
    case ( kNuisFSI084 ) :
    case ( kNuisFSI085 ) :
    case ( kNuisFSI086 ) :
    case ( kNuisFSI087 ) :
    case ( kNuisFSI088 ) :
    case ( kNuisFSI089 ) :
    case ( kNuisFSI090 ) :
    case ( kNuisFSI091 ) :
    case ( kNuisFSI092 ) :
    case ( kNuisFSI093 ) :
    case ( kNuisFSI094 ) :
    case ( kNuisFSI095 ) :
    case ( kNuisFSI096 ) :
    case ( kNuisFSI097 ) :
    case ( kNuisFSI098 ) :
    case ( kNuisFSI099 ) :
    case ( kNuisFSI100 ) :
    case ( kNuisFSI101 ) :
    case ( kNuisFSI102 ) :
    case ( kNuisFSI103 ) :
    case ( kNuisFSI104 ) :
    case ( kNuisFSI105 ) :
    case ( kNuisFSI106 ) :
    case ( kNuisFSI107 ) :
    case ( kNuisFSI108 ) :
    case ( kNuisFSI109 ) :
    case ( kNuisFSI110 ) :
    case ( kNuisFSI111 ) :
    case ( kNuisFSI112 ) :
    case ( kNuisFSI113 ) :
    case ( kNuisFSI114 ) :
    case ( kNuisFSI115 ) :
    case ( kNuisFSI116 ) :
    case ( kNuisFSI117 ) :
    case ( kNuisFSI118 ) :
    case ( kNuisFSI119 ) :
    case ( kNuisFSI120 ) :
    case ( kNuisFSI121 ) :
    case ( kNuisFSI122 ) :
    case ( kNuisFSI123 ) :
    case ( kNuisFSI124 ) :
    case ( kNuisFSI125 ) :
    case ( kNuisFSI126 ) :
    case ( kNuisFSI127 ) :
    case ( kNuisFSI128 ) :
    case ( kNuisFSI129 ) :
    case ( kNuisFSI130 ) :
    case ( kNuisFSI131 ) :
    case ( kNuisFSI132 ) :
    case ( kNuisFSI133 ) :
    case ( kNuisFSI134 ) :
    case ( kNuisFSI135 ) :
    case ( kNuisFSI136 ) :
    case ( kNuisFSI137 ) :
    case ( kNuisFSI138 ) :
    case ( kNuisFSI139 ) :
    case ( kNuisFSI140 ) :
    case ( kNuisFSI141 ) :
    case ( kNuisFSI142 ) :
    case ( kNuisFSI143 ) :
    case ( kNuisFSI144 ) :
    case ( kNuisFSI145 ) :
    case ( kNuisFSI146 ) :
    case ( kNuisFSI147 ) :
    case ( kNuisFSI148 ) :
    case ( kNuisFSI149 ) :
    case ( kNuisFSI150 ) :
    case ( kNuisFSI151 ) :
    case ( kNuisFSI152 ) :
    case ( kNuisFSI153 ) :
    case ( kNuisFSI154 ) :
    case ( kNuisFSI155 ) :
    case ( kNuisFSI156 ) :
    case ( kNuisFSI157 ) :
    case ( kNuisFSI158 ) :
    case ( kNuisFSI159 ) :
    case ( kNuisFSI160 ) :
    case ( kNuisFSI161 ) :
    case ( kNuisFSI162 ) :
    case ( kNuisFSI163 ) :
    case ( kNuisFSI164 ) :
    case ( kNuisFSI165 ) :
    case ( kNuisFSI166 ) :
    case ( kNuisFSI167 ) :
    case ( kNuisFSI168 ) :
    case ( kNuisFSI169 ) :
    case ( kNuisFSI170 ) :
    case ( kNuisFSI171 ) :
    case ( kNuisFSI172 ) :
    case ( kNuisFSI173 ) :
    case ( kNuisFSI174 ) :
    case ( kNuisFSI175 ) :
    case ( kNuisFSI176 ) :
    case ( kNuisFSI177 ) :
    case ( kNuisFSI178 ) :
    case ( kNuisFSI179 ) :
    case ( kNuisFSI180 ) :
    case ( kNuisFSI181 ) :
    case ( kNuisFSI182 ) :
    case ( kNuisFSI183 ) :
    case ( kNuisFSI184 ) :
    case ( kNuisFSI185 ) :
    case ( kNuisFSI186 ) :
    case ( kNuisFSI187 ) :
    case ( kNuisFSI188 ) :
    case ( kNuisFSI189 ) :
    case ( kNuisFSI190 ) :
    case ( kNuisFSI191 ) :
    case ( kNuisFSI192 ) :
    case ( kNuisFSI193 ) :
    case ( kNuisFSI194 ) :
    case ( kNuisFSI195 ) :
    case ( kNuisFSI196 ) :
    case ( kNuisFSI197 ) :
    case ( kNuisFSI198 ) :
    case ( kNuisFSI199 ) :
    case ( kNuisFSI200 ) :
    case ( kNuisFSI201 ) :
    case ( kNuisFSI202 ) :
    case ( kNuisFSI203 ) :
    case ( kNuisFSI204 ) :
    case ( kNuisFSI205 ) :
    case ( kNuisFSI206 ) :
    case ( kNuisFSI207 ) :
    case ( kNuisFSI208 ) :
    case ( kNuisFSI209 ) :
    case ( kNuisFSI210 ) :
    case ( kNuisFSI211 ) :
    case ( kNuisFSI212 ) :
    case ( kNuisFSI213 ) :
    case ( kNuisFSI214 ) :
    case ( kNuisFSI215 ) :
    case ( kNuisFSI216 ) :
    case ( kNuisFSI217 ) :
    case ( kNuisFSI218 ) :
    case ( kNuisFSI219 ) :
    case ( kNuisFSI220 ) :
    case ( kNuisFSI221 ) :
    case ( kNuisFSI222 ) :
    case ( kNuisFSI223 ) :
    case ( kNuisFSI224 ) :
    case ( kNuisFSI225 ) :
    case ( kNuisFSI226 ) :
    case ( kNuisFSI227 ) :
    case ( kNuisFSI228 ) :
    case ( kNuisFSI229 ) :
    case ( kNuisFSI230 ) :
    case ( kNuisFSI231 ) :
    case ( kNuisFSI232 ) :
    case ( kNuisFSI233 ) :
    case ( kNuisFSI234 ) :
    case ( kNuisFSI235 ) :
    case ( kNuisFSI236 ) :
    case ( kNuisFSI237 ) :
    case ( kNuisFSI238 ) :
    case ( kNuisFSI239 ) :
    case ( kNuisFSI240 ) :
    case ( kNuisFSI241 ) :
    case ( kNuisFSI242 ) :
    case ( kNuisFSI243 ) :
    case ( kNuisFSI244 ) :
    case ( kNuisFSI245 ) :
    case ( kNuisFSI246 ) :
    case ( kNuisFSI247 ) :
    case ( kNuisFSI248 ) :
    case ( kNuisFSI249 ) :
    case ( kNuisFSI250 ) :
    case ( kNuisFSI251 ) :
    case ( kNuisFSI252 ) :
    case ( kNuisFSI253 ) :
    case ( kNuisFSI254 ) :
    case ( kNuisFSI255 ) :
    case ( kNuisFSI256 ) :
    case ( kNuisFSI257 ) :
    case ( kNuisFSI258 ) :
    case ( kNuisFSI259 ) :
    case ( kNuisFSI260 ) :
    case ( kNuisFSI261 ) :
    case ( kNuisFSI262 ) :
    case ( kNuisFSI263 ) :
    case ( kNuisFSI264 ) :
    case ( kNuisFSI265 ) :
    case ( kNuisFSI266 ) :
    case ( kNuisFSI267 ) :
    case ( kNuisFSI268 ) :
    case ( kNuisFSI269 ) :
    case ( kNuisFSI270 ) :
    case ( kNuisFSI271 ) :
    case ( kNuisFSI272 ) :
    case ( kNuisFSI273 ) :
    case ( kNuisFSI274 ) :
    case ( kNuisFSI275 ) :
    case ( kNuisFSI276 ) :
    case ( kNuisFSI277 ) :
    case ( kNuisFSI278 ) :
    case ( kNuisFSI279 ) :
    case ( kNuisFSI280 ) :
    case ( kNuisFSI281 ) :
    case ( kNuisFSI282 ) :
    case ( kNuisFSI283 ) :
    case ( kNuisFSI284 ) :
    case ( kNuisFSI285 ) :
    case ( kNuisFSI286 ) :
    case ( kNuisFSI287 ) :
    case ( kNuisFSI288 ) :
    case ( kNuisFSI289 ) :
    case ( kNuisFSI290 ) :
    case ( kNuisFSI291 ) :
    case ( kNuisFSI292 ) :
    case ( kNuisFSI293 ) :
    case ( kNuisFSI294 ) :
    case ( kNuisFSI295 ) :
    case ( kNuisFSI296 ) :
    case ( kNuisFSI297 ) :
    case ( kNuisFSI298 ) :
    case ( kNuisFSI299 ) :
    case ( kNuisFSI300 ) :
    case ( kNuisFSI301 ) :
    case ( kNuisFSI302 ) :
    case ( kNuisFSI303 ) :
    case ( kNuisFSI304 ) :
    case ( kNuisFSI305 ) :
    case ( kNuisFSI306 ) :
    case ( kNuisFSI307 ) :
    case ( kNuisFSI308 ) :
    case ( kNuisFSI309 ) :
    case ( kNuisFSI310 ) :
    case ( kNuisFSI311 ) :
    case ( kNuisFSI312 ) :
    case ( kNuisFSI313 ) :
    case ( kNuisFSI314 ) :
    case ( kNuisFSI315 ) :
    case ( kNuisFSI316 ) :
    case ( kNuisFSI317 ) :
    case ( kNuisFSI318 ) :
    case ( kNuisFSI319 ) :
    case ( kNuisFSI320 ) :
    case ( kNuisFSI321 ) :
    case ( kNuisFSI322 ) :
    case ( kNuisFSI323 ) :
    case ( kNuisFSI324 ) :
    case ( kNuisFSI325 ) :
    case ( kNuisFSI326 ) :
    case ( kNuisFSI327 ) :
    case ( kNuisFSI328 ) :
    case ( kNuisFSI329 ) :
    case ( kNuisFSI330 ) :
    case ( kNuisFSI331 ) :
    case ( kNuisFSI332 ) :
    case ( kNuisFSI333 ) :
    case ( kNuisFSI334 ) :
    case ( kNuisFSI335 ) :
    case ( kNuisFSI336 ) :
    case ( kNuisFSI337 ) :
    case ( kNuisFSI338 ) :
    case ( kNuisFSI339 ) :
    case ( kNuisFSI340 ) :
    case ( kNuisFSI341 ) :
    case ( kNuisFSI342 ) :
    case ( kNuisFSI343 ) :
    case ( kNuisFSI344 ) :
    case ( kNuisFSI345 ) :
    case ( kNuisFSI346 ) :
    case ( kNuisFSI347 ) :
    case ( kNuisFSI348 ) :
    case ( kNuisFSI349 ) :
    case ( kNuisFSI350 ) :
    case ( kNuisFSI351 ) :
    case ( kNuisFSI352 ) :
    case ( kNuisFSI353 ) :
    case ( kNuisFSI354 ) :
    case ( kNuisFSI355 ) :
    case ( kNuisFSI356 ) :
    case ( kNuisFSI357 ) :
    case ( kNuisFSI358 ) :
    case ( kNuisFSI359 ) :
    case ( kNuisFSI360 ) :
    case ( kNuisFSI361 ) :
    case ( kNuisFSI362 ) :
    case ( kNuisFSI363 ) :
    case ( kNuisFSI364 ) :
    case ( kNuisFSI365 ) :
    case ( kNuisFSI366 ) :
    case ( kNuisFSI367 ) :
    case ( kNuisFSI368 ) :
    case ( kNuisFSI369 ) :
    case ( kNuisFSI370 ) :
    case ( kNuisFSI371 ) :
    case ( kNuisFSI372 ) :
    case ( kNuisFSI373 ) :
    case ( kNuisFSI374 ) :
    case ( kNuisFSI375 ) :
    case ( kNuisFSI376 ) :
    case ( kNuisFSI377 ) :
    case ( kNuisFSI378 ) :
    case ( kNuisFSI379 ) :
    case ( kNuisFSI380 ) :
    case ( kNuisFSI381 ) :
    case ( kNuisFSI382 ) :
    case ( kNuisFSI383 ) :
    case ( kNuisFSI384 ) :
    case ( kNuisFSI385 ) :
    case ( kNuisFSI386 ) :
    case ( kNuisFSI387 ) :
    case ( kNuisFSI388 ) :
    case ( kNuisFSI389 ) :
    case ( kNuisFSI390 ) :
    case ( kNuisFSI391 ) :
    case ( kNuisFSI392 ) :
    case ( kNuisFSI393 ) :
    case ( kNuisFSI394 ) :
    case ( kNuisFSI395 ) :
    case ( kNuisFSI396 ) :
    case ( kNuisFSI397 ) :
    case ( kNuisFSI398 ) :
    case ( kNuisFSI399 ) :
    case ( kNuisFSI400 ) :
    case ( kNuisFSI401 ) :
    case ( kNuisFSI402 ) :
    case ( kNuisFSI403 ) :
    case ( kNuisFSI404 ) :
    case ( kNuisFSI405 ) :
    case ( kNuisFSI406 ) :
    case ( kNuisFSI407 ) :
    case ( kNuisFSI408 ) :
    case ( kNuisFSI409 ) :
    case ( kNuisFSI410 ) :
    case ( kNuisFSI411 ) :
    case ( kNuisFSI412 ) :
    case ( kNuisFSI413 ) :
    case ( kNuisFSI414 ) :
    case ( kNuisFSI415 ) :
    case ( kNuisFSI416 ) :
    case ( kNuisFSI417 ) :
    case ( kNuisFSI418 ) :
    case ( kNuisFSI419 ) :
    case ( kNuisFSI420 ) :
    case ( kNuisFSI421 ) :
    case ( kNuisFSI422 ) :
    case ( kNuisFSI423 ) :
    case ( kNuisFSI424 ) :
    case ( kNuisFSI425 ) :
    case ( kNuisFSI426 ) :
    case ( kNuisFSI427 ) :
    case ( kNuisFSI428 ) :
    case ( kNuisFSI429 ) :
    case ( kNuisFSI430 ) :
    case ( kNuisFSI431 ) :
    case ( kNuisFSI432 ) :
    case ( kNuisFSI433 ) :
    case ( kNuisFSI434 ) :
    case ( kNuisFSI435 ) :
    case ( kNuisFSI436 ) :
    case ( kNuisFSI437 ) :
    case ( kNuisFSI438 ) :
    case ( kNuisFSI439 ) :
    case ( kNuisFSI440 ) :
    case ( kNuisFSI441 ) :
    case ( kNuisFSI442 ) :
    case ( kNuisFSI443 ) :
    case ( kNuisFSI444 ) :
    case ( kNuisFSI445 ) :
    case ( kNuisFSI446 ) :
    case ( kNuisFSI447 ) :
    case ( kNuisFSI448 ) :
    case ( kNuisFSI449 ) :
    case ( kNuisFSI450 ) :
    case ( kNuisFSI451 ) :
    case ( kNuisFSI452 ) :
    case ( kNuisFSI453 ) :
    case ( kNuisFSI454 ) :
    case ( kNuisFSI455 ) :
    case ( kNuisFSI456 ) :
    case ( kNuisFSI457 ) :
    case ( kNuisFSI458 ) :
    case ( kNuisFSI459 ) :
    case ( kNuisFSI460 ) :
    case ( kNuisFSI461 ) :
    case ( kNuisFSI462 ) :
    case ( kNuisFSI463 ) :
    case ( kNuisFSI464 ) :
    case ( kNuisFSI465 ) :
    case ( kNuisFSI466 ) :
    case ( kNuisFSI467 ) :
    case ( kNuisFSI468 ) :
    case ( kNuisFSI469 ) :
    case ( kNuisFSI470 ) :
    case ( kNuisFSI471 ) :
    case ( kNuisFSI472 ) :
    case ( kNuisFSI473 ) :
    case ( kNuisFSI474 ) :
    case ( kNuisFSI475 ) :
    case ( kNuisFSI476 ) :
    case ( kNuisFSI477 ) :
    case ( kNuisFSI478 ) :
    case ( kNuisFSI479 ) :
#endif
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
      if(ip <= (int)kNuisFSI020)
	nutypes = "numu/numubar CCQE";
      else if(ip <= (int)kNuisFSI041)
	nutypes = "numu/numubar CCnQE";
      else if(ip <= (int)kNuisFSI062)
	nutypes = "nue/nuebar/signue CC";
      else if(ip <= (int)kNuisFSI078)
	nutypes = "all NC";
#elif FSI_BINS == 4
      if(ip <= (int)kNuisFSI095)
	nutypes = "numu";
      else if(ip <= (int)kNuisFSI191)
	nutypes = "numubar";
      else if(ip <= (int)kNuisFSI287)
	nutypes = "nue";
      else if(ip <= (int)kNuisFSI384)
	nutypes = "nuebar";
      else if(ip <= (int)kNuisFSI479)
	nutypes = "oscillated nue";
#endif
      return Form("FSI/SI %2d; Ereco range %.1f - %.1f GeV; %s",ip,
		  kFSIParamsErecoCut[ip][0],kFSIParamsErecoCut[ip][1],nutypes.c_str());     
      break;
    }//FSI
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
    {
      int ip = (int)p - (int)kNuisBANFF00;
      return Form("BANFF %2d, SK %s flux, %.1f - %.1f GeV",
                   ip, util::Pdg2String(kNuPdgCodeOfFluxParamInBANFFfit[ip]).c_str(), 
                   kEtrueRangeOfFluxParamInBANFFfit[ip][0],
                   kEtrueRangeOfFluxParamInBANFFfit[ip][1]);     
      break;
    }
    case ( kNuisBANFFMAQE    ) : return "BANFF, Ma QE";                           break;
    case ( kNuisBANFFMARES   ) : return "BANFF, Ma RES";                          break;
    case ( kNuisBANFFCCQEE0  ) : return "BANFF; Norm; CC QE; 0.0 - 1.5 GeV";      break;
#if BANFF_BINS == 2
    case ( kNuisBANFFCCQEE1  ) : return "BANFF; Norm; CC QE; 1.5 - 3.5 GeV";      break;
    case ( kNuisBANFFCCQEE2  ) : return "BANFF; Norm; CC QE; 3.5 - 30.0 GeV";     break;
#endif
    case ( kNuisBANFFCC1piE0 ) : return "BANFF; Norm; CC 1pi; 0.0 - 2.5 GeV";     break;
#if BANFF_BINS == 2
    case ( kNuisBANFFCC1piE1 ) : return "BANFF; Norm; CC 1pi; 2.5 - 30.0 GeV";    break;
#endif
    case ( kNuisBANFFNC1pi0  ) : return "BANFF, Norm, NC 1pi0";                   break;

    //
    // XSec params not correlated with flux
    //

    case ( kNuisXSecFermiMom   ) : return "Fermi Momentum";                  break;
    case ( kNuisXSecSpecFunc   ) : return "Spectral Func";                   break; 
    case ( kNuisXSecWshape     ) : return "W-shape";                         break;
    case ( kNuisXSecShapeCCoth ) : return "Shape; CC Oth";                   break;
    case ( kNuisXSecPiLessDelta) : return "Pi-less delta decay";             break;
    case ( kNuisXSecENuShape   ) : return "Shape; CC1pi ENu";                break;

#if BANFF_BINS != 2
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
    case ( kNuisFSI049 ) :
    case ( kNuisFSI050 ) :
    case ( kNuisFSI051 ) :
    case ( kNuisFSI052 ) :
    case ( kNuisFSI053 ) :
    case ( kNuisFSI054 ) :
    case ( kNuisFSI055 ) :
    case ( kNuisFSI056 ) :
    case ( kNuisFSI057 ) :
    case ( kNuisFSI058 ) :
    case ( kNuisFSI059 ) :
    case ( kNuisFSI060 ) :
    case ( kNuisFSI061 ) :
    case ( kNuisFSI062 ) :
    case ( kNuisFSI063 ) :
    case ( kNuisFSI064 ) :
    case ( kNuisFSI065 ) :
    case ( kNuisFSI066 ) :
    case ( kNuisFSI067 ) :
    case ( kNuisFSI068 ) :
    case ( kNuisFSI069 ) :
    case ( kNuisFSI070 ) :
    case ( kNuisFSI071 ) :
    case ( kNuisFSI072 ) :
    case ( kNuisFSI073 ) :
    case ( kNuisFSI074 ) :
    case ( kNuisFSI075 ) :
    case ( kNuisFSI076 ) :
    case ( kNuisFSI077 ) :
    case ( kNuisFSI078 ) :
#if FSI_BINS >= 4
    case ( kNuisFSI079 ) :
    case ( kNuisFSI080 ) :
    case ( kNuisFSI081 ) :
    case ( kNuisFSI082 ) :
    case ( kNuisFSI083 ) :
    case ( kNuisFSI084 ) :
    case ( kNuisFSI085 ) :
    case ( kNuisFSI086 ) :
    case ( kNuisFSI087 ) :
    case ( kNuisFSI088 ) :
    case ( kNuisFSI089 ) :
    case ( kNuisFSI090 ) :
    case ( kNuisFSI091 ) :
    case ( kNuisFSI092 ) :
    case ( kNuisFSI093 ) :
    case ( kNuisFSI094 ) :
    case ( kNuisFSI095 ) :
    case ( kNuisFSI096 ) :
    case ( kNuisFSI097 ) :
    case ( kNuisFSI098 ) :
    case ( kNuisFSI099 ) :
    case ( kNuisFSI100 ) :
    case ( kNuisFSI101 ) :
    case ( kNuisFSI102 ) :
    case ( kNuisFSI103 ) :
    case ( kNuisFSI104 ) :
    case ( kNuisFSI105 ) :
    case ( kNuisFSI106 ) :
    case ( kNuisFSI107 ) :
    case ( kNuisFSI108 ) :
    case ( kNuisFSI109 ) :
    case ( kNuisFSI110 ) :
    case ( kNuisFSI111 ) :
    case ( kNuisFSI112 ) :
    case ( kNuisFSI113 ) :
    case ( kNuisFSI114 ) :
    case ( kNuisFSI115 ) :
    case ( kNuisFSI116 ) :
    case ( kNuisFSI117 ) :
    case ( kNuisFSI118 ) :
    case ( kNuisFSI119 ) :
    case ( kNuisFSI120 ) :
    case ( kNuisFSI121 ) :
    case ( kNuisFSI122 ) :
    case ( kNuisFSI123 ) :
    case ( kNuisFSI124 ) :
    case ( kNuisFSI125 ) :
    case ( kNuisFSI126 ) :
    case ( kNuisFSI127 ) :
    case ( kNuisFSI128 ) :
    case ( kNuisFSI129 ) :
    case ( kNuisFSI130 ) :
    case ( kNuisFSI131 ) :
    case ( kNuisFSI132 ) :
    case ( kNuisFSI133 ) :
    case ( kNuisFSI134 ) :
    case ( kNuisFSI135 ) :
    case ( kNuisFSI136 ) :
    case ( kNuisFSI137 ) :
    case ( kNuisFSI138 ) :
    case ( kNuisFSI139 ) :
    case ( kNuisFSI140 ) :
    case ( kNuisFSI141 ) :
    case ( kNuisFSI142 ) :
    case ( kNuisFSI143 ) :
    case ( kNuisFSI144 ) :
    case ( kNuisFSI145 ) :
    case ( kNuisFSI146 ) :
    case ( kNuisFSI147 ) :
    case ( kNuisFSI148 ) :
    case ( kNuisFSI149 ) :
    case ( kNuisFSI150 ) :
    case ( kNuisFSI151 ) :
    case ( kNuisFSI152 ) :
    case ( kNuisFSI153 ) :
    case ( kNuisFSI154 ) :
    case ( kNuisFSI155 ) :
    case ( kNuisFSI156 ) :
    case ( kNuisFSI157 ) :
    case ( kNuisFSI158 ) :
    case ( kNuisFSI159 ) :
    case ( kNuisFSI160 ) :
    case ( kNuisFSI161 ) :
    case ( kNuisFSI162 ) :
    case ( kNuisFSI163 ) :
    case ( kNuisFSI164 ) :
    case ( kNuisFSI165 ) :
    case ( kNuisFSI166 ) :
    case ( kNuisFSI167 ) :
    case ( kNuisFSI168 ) :
    case ( kNuisFSI169 ) :
    case ( kNuisFSI170 ) :
    case ( kNuisFSI171 ) :
    case ( kNuisFSI172 ) :
    case ( kNuisFSI173 ) :
    case ( kNuisFSI174 ) :
    case ( kNuisFSI175 ) :
    case ( kNuisFSI176 ) :
    case ( kNuisFSI177 ) :
    case ( kNuisFSI178 ) :
    case ( kNuisFSI179 ) :
    case ( kNuisFSI180 ) :
    case ( kNuisFSI181 ) :
    case ( kNuisFSI182 ) :
    case ( kNuisFSI183 ) :
    case ( kNuisFSI184 ) :
    case ( kNuisFSI185 ) :
    case ( kNuisFSI186 ) :
    case ( kNuisFSI187 ) :
    case ( kNuisFSI188 ) :
    case ( kNuisFSI189 ) :
    case ( kNuisFSI190 ) :
    case ( kNuisFSI191 ) :
    case ( kNuisFSI192 ) :
    case ( kNuisFSI193 ) :
    case ( kNuisFSI194 ) :
    case ( kNuisFSI195 ) :
    case ( kNuisFSI196 ) :
    case ( kNuisFSI197 ) :
    case ( kNuisFSI198 ) :
    case ( kNuisFSI199 ) :
    case ( kNuisFSI200 ) :
    case ( kNuisFSI201 ) :
    case ( kNuisFSI202 ) :
    case ( kNuisFSI203 ) :
    case ( kNuisFSI204 ) :
    case ( kNuisFSI205 ) :
    case ( kNuisFSI206 ) :
    case ( kNuisFSI207 ) :
    case ( kNuisFSI208 ) :
    case ( kNuisFSI209 ) :
    case ( kNuisFSI210 ) :
    case ( kNuisFSI211 ) :
    case ( kNuisFSI212 ) :
    case ( kNuisFSI213 ) :
    case ( kNuisFSI214 ) :
    case ( kNuisFSI215 ) :
    case ( kNuisFSI216 ) :
    case ( kNuisFSI217 ) :
    case ( kNuisFSI218 ) :
    case ( kNuisFSI219 ) :
    case ( kNuisFSI220 ) :
    case ( kNuisFSI221 ) :
    case ( kNuisFSI222 ) :
    case ( kNuisFSI223 ) :
    case ( kNuisFSI224 ) :
    case ( kNuisFSI225 ) :
    case ( kNuisFSI226 ) :
    case ( kNuisFSI227 ) :
    case ( kNuisFSI228 ) :
    case ( kNuisFSI229 ) :
    case ( kNuisFSI230 ) :
    case ( kNuisFSI231 ) :
    case ( kNuisFSI232 ) :
    case ( kNuisFSI233 ) :
    case ( kNuisFSI234 ) :
    case ( kNuisFSI235 ) :
    case ( kNuisFSI236 ) :
    case ( kNuisFSI237 ) :
    case ( kNuisFSI238 ) :
    case ( kNuisFSI239 ) :
    case ( kNuisFSI240 ) :
    case ( kNuisFSI241 ) :
    case ( kNuisFSI242 ) :
    case ( kNuisFSI243 ) :
    case ( kNuisFSI244 ) :
    case ( kNuisFSI245 ) :
    case ( kNuisFSI246 ) :
    case ( kNuisFSI247 ) :
    case ( kNuisFSI248 ) :
    case ( kNuisFSI249 ) :
    case ( kNuisFSI250 ) :
    case ( kNuisFSI251 ) :
    case ( kNuisFSI252 ) :
    case ( kNuisFSI253 ) :
    case ( kNuisFSI254 ) :
    case ( kNuisFSI255 ) :
    case ( kNuisFSI256 ) :
    case ( kNuisFSI257 ) :
    case ( kNuisFSI258 ) :
    case ( kNuisFSI259 ) :
    case ( kNuisFSI260 ) :
    case ( kNuisFSI261 ) :
    case ( kNuisFSI262 ) :
    case ( kNuisFSI263 ) :
    case ( kNuisFSI264 ) :
    case ( kNuisFSI265 ) :
    case ( kNuisFSI266 ) :
    case ( kNuisFSI267 ) :
    case ( kNuisFSI268 ) :
    case ( kNuisFSI269 ) :
    case ( kNuisFSI270 ) :
    case ( kNuisFSI271 ) :
    case ( kNuisFSI272 ) :
    case ( kNuisFSI273 ) :
    case ( kNuisFSI274 ) :
    case ( kNuisFSI275 ) :
    case ( kNuisFSI276 ) :
    case ( kNuisFSI277 ) :
    case ( kNuisFSI278 ) :
    case ( kNuisFSI279 ) :
    case ( kNuisFSI280 ) :
    case ( kNuisFSI281 ) :
    case ( kNuisFSI282 ) :
    case ( kNuisFSI283 ) :
    case ( kNuisFSI284 ) :
    case ( kNuisFSI285 ) :
    case ( kNuisFSI286 ) :
    case ( kNuisFSI287 ) :
    case ( kNuisFSI288 ) :
    case ( kNuisFSI289 ) :
    case ( kNuisFSI290 ) :
    case ( kNuisFSI291 ) :
    case ( kNuisFSI292 ) :
    case ( kNuisFSI293 ) :
    case ( kNuisFSI294 ) :
    case ( kNuisFSI295 ) :
    case ( kNuisFSI296 ) :
    case ( kNuisFSI297 ) :
    case ( kNuisFSI298 ) :
    case ( kNuisFSI299 ) :
    case ( kNuisFSI300 ) :
    case ( kNuisFSI301 ) :
    case ( kNuisFSI302 ) :
    case ( kNuisFSI303 ) :
    case ( kNuisFSI304 ) :
    case ( kNuisFSI305 ) :
    case ( kNuisFSI306 ) :
    case ( kNuisFSI307 ) :
    case ( kNuisFSI308 ) :
    case ( kNuisFSI309 ) :
    case ( kNuisFSI310 ) :
    case ( kNuisFSI311 ) :
    case ( kNuisFSI312 ) :
    case ( kNuisFSI313 ) :
    case ( kNuisFSI314 ) :
    case ( kNuisFSI315 ) :
    case ( kNuisFSI316 ) :
    case ( kNuisFSI317 ) :
    case ( kNuisFSI318 ) :
    case ( kNuisFSI319 ) :
    case ( kNuisFSI320 ) :
    case ( kNuisFSI321 ) :
    case ( kNuisFSI322 ) :
    case ( kNuisFSI323 ) :
    case ( kNuisFSI324 ) :
    case ( kNuisFSI325 ) :
    case ( kNuisFSI326 ) :
    case ( kNuisFSI327 ) :
    case ( kNuisFSI328 ) :
    case ( kNuisFSI329 ) :
    case ( kNuisFSI330 ) :
    case ( kNuisFSI331 ) :
    case ( kNuisFSI332 ) :
    case ( kNuisFSI333 ) :
    case ( kNuisFSI334 ) :
    case ( kNuisFSI335 ) :
    case ( kNuisFSI336 ) :
    case ( kNuisFSI337 ) :
    case ( kNuisFSI338 ) :
    case ( kNuisFSI339 ) :
    case ( kNuisFSI340 ) :
    case ( kNuisFSI341 ) :
    case ( kNuisFSI342 ) :
    case ( kNuisFSI343 ) :
    case ( kNuisFSI344 ) :
    case ( kNuisFSI345 ) :
    case ( kNuisFSI346 ) :
    case ( kNuisFSI347 ) :
    case ( kNuisFSI348 ) :
    case ( kNuisFSI349 ) :
    case ( kNuisFSI350 ) :
    case ( kNuisFSI351 ) :
    case ( kNuisFSI352 ) :
    case ( kNuisFSI353 ) :
    case ( kNuisFSI354 ) :
    case ( kNuisFSI355 ) :
    case ( kNuisFSI356 ) :
    case ( kNuisFSI357 ) :
    case ( kNuisFSI358 ) :
    case ( kNuisFSI359 ) :
    case ( kNuisFSI360 ) :
    case ( kNuisFSI361 ) :
    case ( kNuisFSI362 ) :
    case ( kNuisFSI363 ) :
    case ( kNuisFSI364 ) :
    case ( kNuisFSI365 ) :
    case ( kNuisFSI366 ) :
    case ( kNuisFSI367 ) :
    case ( kNuisFSI368 ) :
    case ( kNuisFSI369 ) :
    case ( kNuisFSI370 ) :
    case ( kNuisFSI371 ) :
    case ( kNuisFSI372 ) :
    case ( kNuisFSI373 ) :
    case ( kNuisFSI374 ) :
    case ( kNuisFSI375 ) :
    case ( kNuisFSI376 ) :
    case ( kNuisFSI377 ) :
    case ( kNuisFSI378 ) :
    case ( kNuisFSI379 ) :
    case ( kNuisFSI380 ) :
    case ( kNuisFSI381 ) :
    case ( kNuisFSI382 ) :
    case ( kNuisFSI383 ) :
    case ( kNuisFSI384 ) :
    case ( kNuisFSI385 ) :
    case ( kNuisFSI386 ) :
    case ( kNuisFSI387 ) :
    case ( kNuisFSI388 ) :
    case ( kNuisFSI389 ) :
    case ( kNuisFSI390 ) :
    case ( kNuisFSI391 ) :
    case ( kNuisFSI392 ) :
    case ( kNuisFSI393 ) :
    case ( kNuisFSI394 ) :
    case ( kNuisFSI395 ) :
    case ( kNuisFSI396 ) :
    case ( kNuisFSI397 ) :
    case ( kNuisFSI398 ) :
    case ( kNuisFSI399 ) :
    case ( kNuisFSI400 ) :
    case ( kNuisFSI401 ) :
    case ( kNuisFSI402 ) :
    case ( kNuisFSI403 ) :
    case ( kNuisFSI404 ) :
    case ( kNuisFSI405 ) :
    case ( kNuisFSI406 ) :
    case ( kNuisFSI407 ) :
    case ( kNuisFSI408 ) :
    case ( kNuisFSI409 ) :
    case ( kNuisFSI410 ) :
    case ( kNuisFSI411 ) :
    case ( kNuisFSI412 ) :
    case ( kNuisFSI413 ) :
    case ( kNuisFSI414 ) :
    case ( kNuisFSI415 ) :
    case ( kNuisFSI416 ) :
    case ( kNuisFSI417 ) :
    case ( kNuisFSI418 ) :
    case ( kNuisFSI419 ) :
    case ( kNuisFSI420 ) :
    case ( kNuisFSI421 ) :
    case ( kNuisFSI422 ) :
    case ( kNuisFSI423 ) :
    case ( kNuisFSI424 ) :
    case ( kNuisFSI425 ) :
    case ( kNuisFSI426 ) :
    case ( kNuisFSI427 ) :
    case ( kNuisFSI428 ) :
    case ( kNuisFSI429 ) :
    case ( kNuisFSI430 ) :
    case ( kNuisFSI431 ) :
    case ( kNuisFSI432 ) :
    case ( kNuisFSI433 ) :
    case ( kNuisFSI434 ) :
    case ( kNuisFSI435 ) :
    case ( kNuisFSI436 ) :
    case ( kNuisFSI437 ) :
    case ( kNuisFSI438 ) :
    case ( kNuisFSI439 ) :
    case ( kNuisFSI440 ) :
    case ( kNuisFSI441 ) :
    case ( kNuisFSI442 ) :
    case ( kNuisFSI443 ) :
    case ( kNuisFSI444 ) :
    case ( kNuisFSI445 ) :
    case ( kNuisFSI446 ) :
    case ( kNuisFSI447 ) :
    case ( kNuisFSI448 ) :
    case ( kNuisFSI449 ) :
    case ( kNuisFSI450 ) :
    case ( kNuisFSI451 ) :
    case ( kNuisFSI452 ) :
    case ( kNuisFSI453 ) :
    case ( kNuisFSI454 ) :
    case ( kNuisFSI455 ) :
    case ( kNuisFSI456 ) :
    case ( kNuisFSI457 ) :
    case ( kNuisFSI458 ) :
    case ( kNuisFSI459 ) :
    case ( kNuisFSI460 ) :
    case ( kNuisFSI461 ) :
    case ( kNuisFSI462 ) :
    case ( kNuisFSI463 ) :
    case ( kNuisFSI464 ) :
    case ( kNuisFSI465 ) :
    case ( kNuisFSI466 ) :
    case ( kNuisFSI467 ) :
    case ( kNuisFSI468 ) :
    case ( kNuisFSI469 ) :
    case ( kNuisFSI470 ) :
    case ( kNuisFSI471 ) :
    case ( kNuisFSI472 ) :
    case ( kNuisFSI473 ) :
    case ( kNuisFSI474 ) :
    case ( kNuisFSI475 ) :
    case ( kNuisFSI476 ) :
    case ( kNuisFSI477 ) :
    case ( kNuisFSI478 ) :
    case ( kNuisFSI479 ) :
#endif
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
      if(ip <= (int)kNuisFSI020)
	nutypes = "$\\nu_\\mu$/$\\bar{\\nu_\\mu}$ CCQE";
      else if(ip <= (int)kNuisFSI041)
	nutypes = "$\\nu_\\mu$/$\\bar{\\nu_\\mu}$ CCnQE";
      else if(ip <= (int)kNuisFSI062)
	nutypes = "$\\nu_e$/$\\bar{\\nu_e}$/sig$\\nu_e$ CC";
      else if(ip <= (int)kNuisFSI078)
	nutypes = "all NC";
#elif FSI_BINS == 4
      if(ip <= (int)kNuisFSI095)
	nutypes = "$\\nu_\\mu$";
      else if(ip <= (int)kNuisFSI191)
	nutypes = "$\\bar{\\nu_\\mu}$";
      else if(ip <= (int)kNuisFSI287)
	nutypes = "$\\nu_e$";
      else if(ip <= (int)kNuisFSI384)
	nutypes = "$\\bar{\\nu_e}$";
      else if(ip <= (int)kNuisFSI479)
	nutypes = "Oscillated $\\nu_e$";
#endif
      return Form("FSI/SI %2d; Ereco range %.1f - %.1f GeV; %s",ip,
		  kFSIParamsErecoCut[ip][0],kFSIParamsErecoCut[ip][1],nutypes.c_str());     
      break;
    }//FSI
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

    case ( kNuisBANFFMAQE             ) : return "f_banffmaqe";    break;
    case ( kNuisBANFFMARES            ) : return "f_banffmares";   break;
    case ( kNuisBANFFCCQEE0           ) : return "f_banffccqee0";  break;
#if BANFF_BINS == 2
    case ( kNuisBANFFCCQEE1           ) : return "f_banffccqee1";  break;
    case ( kNuisBANFFCCQEE2           ) : return "f_banffccqee2";  break;
#endif
    case ( kNuisBANFFCC1piE0          ) : return "f_banffcc1pie0"; break;
#if BANFF_BINS == 2
    case ( kNuisBANFFCC1piE1          ) : return "f_banffcc1pie1"; break;
#endif
    case ( kNuisBANFFNC1pi0           ) : return "f_banffnc1pi0";  break;

    case ( kNuisXSecFermiMom          ) : return "f_pfermi";       break;
    case ( kNuisXSecSpecFunc          ) : return "f_specfunc";     break; 
    case ( kNuisXSecWshape            ) : return "f_wshape";       break;
    case ( kNuisXSecShapeCCoth        ) : return "f_shapeCCoth";   break;
    case ( kNuisXSecPiLessDelta       ) : return "f_pilessdelta";  break;
    case ( kNuisXSecENuShape          ) : return "f_enushape";     break;

#if BANFF_BINS != 2
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
    case ( kNuisFSI049 ) : return "f_fsi049"; break;
    case ( kNuisFSI050 ) : return "f_fsi050"; break;
    case ( kNuisFSI051 ) : return "f_fsi051"; break;
    case ( kNuisFSI052 ) : return "f_fsi052"; break;
    case ( kNuisFSI053 ) : return "f_fsi053"; break;
    case ( kNuisFSI054 ) : return "f_fsi054"; break;
    case ( kNuisFSI055 ) : return "f_fsi055"; break;
    case ( kNuisFSI056 ) : return "f_fsi056"; break;
    case ( kNuisFSI057 ) : return "f_fsi057"; break;
    case ( kNuisFSI058 ) : return "f_fsi058"; break;
    case ( kNuisFSI059 ) : return "f_fsi059"; break;
    case ( kNuisFSI060 ) : return "f_fsi060"; break;
    case ( kNuisFSI061 ) : return "f_fsi061"; break;
    case ( kNuisFSI062 ) : return "f_fsi062"; break;
    case ( kNuisFSI063 ) : return "f_fsi063"; break;
    case ( kNuisFSI064 ) : return "f_fsi064"; break;
    case ( kNuisFSI065 ) : return "f_fsi065"; break;
    case ( kNuisFSI066 ) : return "f_fsi066"; break;
    case ( kNuisFSI067 ) : return "f_fsi067"; break;
    case ( kNuisFSI068 ) : return "f_fsi068"; break;
    case ( kNuisFSI069 ) : return "f_fsi069"; break;
    case ( kNuisFSI070 ) : return "f_fsi070"; break;
    case ( kNuisFSI071 ) : return "f_fsi071"; break;
    case ( kNuisFSI072 ) : return "f_fsi072"; break;
    case ( kNuisFSI073 ) : return "f_fsi073"; break;
    case ( kNuisFSI074 ) : return "f_fsi074"; break;
    case ( kNuisFSI075 ) : return "f_fsi075"; break;
    case ( kNuisFSI076 ) : return "f_fsi076"; break;
    case ( kNuisFSI077 ) : return "f_fsi077"; break;
    case ( kNuisFSI078 ) : return "f_fsi078"; break;
#if FSI_BINS >= 4
    case ( kNuisFSI079 ) : return "f_fsi079"; break;
    case ( kNuisFSI080 ) : return "f_fsi080"; break;
    case ( kNuisFSI081 ) : return "f_fsi081"; break;
    case ( kNuisFSI082 ) : return "f_fsi082"; break;
    case ( kNuisFSI083 ) : return "f_fsi083"; break;
    case ( kNuisFSI084 ) : return "f_fsi084"; break;
    case ( kNuisFSI085 ) : return "f_fsi085"; break;
    case ( kNuisFSI086 ) : return "f_fsi086"; break;
    case ( kNuisFSI087 ) : return "f_fsi087"; break;
    case ( kNuisFSI088 ) : return "f_fsi088"; break;
    case ( kNuisFSI089 ) : return "f_fsi089"; break;
    case ( kNuisFSI090 ) : return "f_fsi090"; break;
    case ( kNuisFSI091 ) : return "f_fsi091"; break;
    case ( kNuisFSI092 ) : return "f_fsi092"; break;
    case ( kNuisFSI093 ) : return "f_fsi093"; break;
    case ( kNuisFSI094 ) : return "f_fsi094"; break;
    case ( kNuisFSI095 ) : return "f_fsi095"; break;
    case ( kNuisFSI096 ) : return "f_fsi096"; break;
    case ( kNuisFSI097 ) : return "f_fsi097"; break;
    case ( kNuisFSI098 ) : return "f_fsi098"; break;
    case ( kNuisFSI099 ) : return "f_fsi099"; break;
    case ( kNuisFSI100 ) : return "f_fsi100"; break;
    case ( kNuisFSI101 ) : return "f_fsi101"; break;
    case ( kNuisFSI102 ) : return "f_fsi102"; break;
    case ( kNuisFSI103 ) : return "f_fsi103"; break;
    case ( kNuisFSI104 ) : return "f_fsi104"; break;
    case ( kNuisFSI105 ) : return "f_fsi105"; break;
    case ( kNuisFSI106 ) : return "f_fsi106"; break;
    case ( kNuisFSI107 ) : return "f_fsi107"; break;
    case ( kNuisFSI108 ) : return "f_fsi108"; break;
    case ( kNuisFSI109 ) : return "f_fsi109"; break;
    case ( kNuisFSI110 ) : return "f_fsi110"; break;
    case ( kNuisFSI111 ) : return "f_fsi111"; break;
    case ( kNuisFSI112 ) : return "f_fsi112"; break;
    case ( kNuisFSI113 ) : return "f_fsi113"; break;
    case ( kNuisFSI114 ) : return "f_fsi114"; break;
    case ( kNuisFSI115 ) : return "f_fsi115"; break;
    case ( kNuisFSI116 ) : return "f_fsi116"; break;
    case ( kNuisFSI117 ) : return "f_fsi117"; break;
    case ( kNuisFSI118 ) : return "f_fsi118"; break;
    case ( kNuisFSI119 ) : return "f_fsi119"; break;
    case ( kNuisFSI120 ) : return "f_fsi120"; break;
    case ( kNuisFSI121 ) : return "f_fsi121"; break;
    case ( kNuisFSI122 ) : return "f_fsi122"; break;
    case ( kNuisFSI123 ) : return "f_fsi123"; break;
    case ( kNuisFSI124 ) : return "f_fsi124"; break;
    case ( kNuisFSI125 ) : return "f_fsi125"; break;
    case ( kNuisFSI126 ) : return "f_fsi126"; break;
    case ( kNuisFSI127 ) : return "f_fsi127"; break;
    case ( kNuisFSI128 ) : return "f_fsi128"; break;
    case ( kNuisFSI129 ) : return "f_fsi129"; break;
    case ( kNuisFSI130 ) : return "f_fsi130"; break;
    case ( kNuisFSI131 ) : return "f_fsi131"; break;
    case ( kNuisFSI132 ) : return "f_fsi132"; break;
    case ( kNuisFSI133 ) : return "f_fsi133"; break;
    case ( kNuisFSI134 ) : return "f_fsi134"; break;
    case ( kNuisFSI135 ) : return "f_fsi135"; break;
    case ( kNuisFSI136 ) : return "f_fsi136"; break;
    case ( kNuisFSI137 ) : return "f_fsi137"; break;
    case ( kNuisFSI138 ) : return "f_fsi138"; break;
    case ( kNuisFSI139 ) : return "f_fsi139"; break;
    case ( kNuisFSI140 ) : return "f_fsi140"; break;
    case ( kNuisFSI141 ) : return "f_fsi141"; break;
    case ( kNuisFSI142 ) : return "f_fsi142"; break;
    case ( kNuisFSI143 ) : return "f_fsi143"; break;
    case ( kNuisFSI144 ) : return "f_fsi144"; break;
    case ( kNuisFSI145 ) : return "f_fsi145"; break;
    case ( kNuisFSI146 ) : return "f_fsi146"; break;
    case ( kNuisFSI147 ) : return "f_fsi147"; break;
    case ( kNuisFSI148 ) : return "f_fsi148"; break;
    case ( kNuisFSI149 ) : return "f_fsi149"; break;
    case ( kNuisFSI150 ) : return "f_fsi150"; break;
    case ( kNuisFSI151 ) : return "f_fsi151"; break;
    case ( kNuisFSI152 ) : return "f_fsi152"; break;
    case ( kNuisFSI153 ) : return "f_fsi153"; break;
    case ( kNuisFSI154 ) : return "f_fsi154"; break;
    case ( kNuisFSI155 ) : return "f_fsi155"; break;
    case ( kNuisFSI156 ) : return "f_fsi156"; break;
    case ( kNuisFSI157 ) : return "f_fsi157"; break;
    case ( kNuisFSI158 ) : return "f_fsi158"; break;
    case ( kNuisFSI159 ) : return "f_fsi159"; break;
    case ( kNuisFSI160 ) : return "f_fsi160"; break;
    case ( kNuisFSI161 ) : return "f_fsi161"; break;
    case ( kNuisFSI162 ) : return "f_fsi162"; break;
    case ( kNuisFSI163 ) : return "f_fsi163"; break;
    case ( kNuisFSI164 ) : return "f_fsi164"; break;
    case ( kNuisFSI165 ) : return "f_fsi165"; break;
    case ( kNuisFSI166 ) : return "f_fsi166"; break;
    case ( kNuisFSI167 ) : return "f_fsi167"; break;
    case ( kNuisFSI168 ) : return "f_fsi168"; break;
    case ( kNuisFSI169 ) : return "f_fsi169"; break;
    case ( kNuisFSI170 ) : return "f_fsi170"; break;
    case ( kNuisFSI171 ) : return "f_fsi171"; break;
    case ( kNuisFSI172 ) : return "f_fsi172"; break;
    case ( kNuisFSI173 ) : return "f_fsi173"; break;
    case ( kNuisFSI174 ) : return "f_fsi174"; break;
    case ( kNuisFSI175 ) : return "f_fsi175"; break;
    case ( kNuisFSI176 ) : return "f_fsi176"; break;
    case ( kNuisFSI177 ) : return "f_fsi177"; break;
    case ( kNuisFSI178 ) : return "f_fsi178"; break;
    case ( kNuisFSI179 ) : return "f_fsi179"; break;
    case ( kNuisFSI180 ) : return "f_fsi180"; break;
    case ( kNuisFSI181 ) : return "f_fsi181"; break;
    case ( kNuisFSI182 ) : return "f_fsi182"; break;
    case ( kNuisFSI183 ) : return "f_fsi183"; break;
    case ( kNuisFSI184 ) : return "f_fsi184"; break;
    case ( kNuisFSI185 ) : return "f_fsi185"; break;
    case ( kNuisFSI186 ) : return "f_fsi186"; break;
    case ( kNuisFSI187 ) : return "f_fsi187"; break;
    case ( kNuisFSI188 ) : return "f_fsi188"; break;
    case ( kNuisFSI189 ) : return "f_fsi189"; break;
    case ( kNuisFSI190 ) : return "f_fsi190"; break;
    case ( kNuisFSI191 ) : return "f_fsi191"; break;
    case ( kNuisFSI192 ) : return "f_fsi192"; break;
    case ( kNuisFSI193 ) : return "f_fsi193"; break;
    case ( kNuisFSI194 ) : return "f_fsi194"; break;
    case ( kNuisFSI195 ) : return "f_fsi195"; break;
    case ( kNuisFSI196 ) : return "f_fsi196"; break;
    case ( kNuisFSI197 ) : return "f_fsi197"; break;
    case ( kNuisFSI198 ) : return "f_fsi198"; break;
    case ( kNuisFSI199 ) : return "f_fsi199"; break;
    case ( kNuisFSI200 ) : return "f_fsi200"; break;
    case ( kNuisFSI201 ) : return "f_fsi201"; break;
    case ( kNuisFSI202 ) : return "f_fsi202"; break;
    case ( kNuisFSI203 ) : return "f_fsi203"; break;
    case ( kNuisFSI204 ) : return "f_fsi204"; break;
    case ( kNuisFSI205 ) : return "f_fsi205"; break;
    case ( kNuisFSI206 ) : return "f_fsi206"; break;
    case ( kNuisFSI207 ) : return "f_fsi207"; break;
    case ( kNuisFSI208 ) : return "f_fsi208"; break;
    case ( kNuisFSI209 ) : return "f_fsi209"; break;
    case ( kNuisFSI210 ) : return "f_fsi210"; break;
    case ( kNuisFSI211 ) : return "f_fsi211"; break;
    case ( kNuisFSI212 ) : return "f_fsi212"; break;
    case ( kNuisFSI213 ) : return "f_fsi213"; break;
    case ( kNuisFSI214 ) : return "f_fsi214"; break;
    case ( kNuisFSI215 ) : return "f_fsi215"; break;
    case ( kNuisFSI216 ) : return "f_fsi216"; break;
    case ( kNuisFSI217 ) : return "f_fsi217"; break;
    case ( kNuisFSI218 ) : return "f_fsi218"; break;
    case ( kNuisFSI219 ) : return "f_fsi219"; break;
    case ( kNuisFSI220 ) : return "f_fsi220"; break;
    case ( kNuisFSI221 ) : return "f_fsi221"; break;
    case ( kNuisFSI222 ) : return "f_fsi222"; break;
    case ( kNuisFSI223 ) : return "f_fsi223"; break;
    case ( kNuisFSI224 ) : return "f_fsi224"; break;
    case ( kNuisFSI225 ) : return "f_fsi225"; break;
    case ( kNuisFSI226 ) : return "f_fsi226"; break;
    case ( kNuisFSI227 ) : return "f_fsi227"; break;
    case ( kNuisFSI228 ) : return "f_fsi228"; break;
    case ( kNuisFSI229 ) : return "f_fsi229"; break;
    case ( kNuisFSI230 ) : return "f_fsi230"; break;
    case ( kNuisFSI231 ) : return "f_fsi231"; break;
    case ( kNuisFSI232 ) : return "f_fsi232"; break;
    case ( kNuisFSI233 ) : return "f_fsi233"; break;
    case ( kNuisFSI234 ) : return "f_fsi234"; break;
    case ( kNuisFSI235 ) : return "f_fsi235"; break;
    case ( kNuisFSI236 ) : return "f_fsi236"; break;
    case ( kNuisFSI237 ) : return "f_fsi237"; break;
    case ( kNuisFSI238 ) : return "f_fsi238"; break;
    case ( kNuisFSI239 ) : return "f_fsi239"; break;
    case ( kNuisFSI240 ) : return "f_fsi240"; break;
    case ( kNuisFSI241 ) : return "f_fsi241"; break;
    case ( kNuisFSI242 ) : return "f_fsi242"; break;
    case ( kNuisFSI243 ) : return "f_fsi243"; break;
    case ( kNuisFSI244 ) : return "f_fsi244"; break;
    case ( kNuisFSI245 ) : return "f_fsi245"; break;
    case ( kNuisFSI246 ) : return "f_fsi246"; break;
    case ( kNuisFSI247 ) : return "f_fsi247"; break;
    case ( kNuisFSI248 ) : return "f_fsi248"; break;
    case ( kNuisFSI249 ) : return "f_fsi249"; break;
    case ( kNuisFSI250 ) : return "f_fsi250"; break;
    case ( kNuisFSI251 ) : return "f_fsi251"; break;
    case ( kNuisFSI252 ) : return "f_fsi252"; break;
    case ( kNuisFSI253 ) : return "f_fsi253"; break;
    case ( kNuisFSI254 ) : return "f_fsi254"; break;
    case ( kNuisFSI255 ) : return "f_fsi255"; break;
    case ( kNuisFSI256 ) : return "f_fsi256"; break;
    case ( kNuisFSI257 ) : return "f_fsi257"; break;
    case ( kNuisFSI258 ) : return "f_fsi258"; break;
    case ( kNuisFSI259 ) : return "f_fsi259"; break;
    case ( kNuisFSI260 ) : return "f_fsi260"; break;
    case ( kNuisFSI261 ) : return "f_fsi261"; break;
    case ( kNuisFSI262 ) : return "f_fsi262"; break;
    case ( kNuisFSI263 ) : return "f_fsi263"; break;
    case ( kNuisFSI264 ) : return "f_fsi264"; break;
    case ( kNuisFSI265 ) : return "f_fsi265"; break;
    case ( kNuisFSI266 ) : return "f_fsi266"; break;
    case ( kNuisFSI267 ) : return "f_fsi267"; break;
    case ( kNuisFSI268 ) : return "f_fsi268"; break;
    case ( kNuisFSI269 ) : return "f_fsi269"; break;
    case ( kNuisFSI270 ) : return "f_fsi270"; break;
    case ( kNuisFSI271 ) : return "f_fsi271"; break;
    case ( kNuisFSI272 ) : return "f_fsi272"; break;
    case ( kNuisFSI273 ) : return "f_fsi273"; break;
    case ( kNuisFSI274 ) : return "f_fsi274"; break;
    case ( kNuisFSI275 ) : return "f_fsi275"; break;
    case ( kNuisFSI276 ) : return "f_fsi276"; break;
    case ( kNuisFSI277 ) : return "f_fsi277"; break;
    case ( kNuisFSI278 ) : return "f_fsi278"; break;
    case ( kNuisFSI279 ) : return "f_fsi279"; break;
    case ( kNuisFSI280 ) : return "f_fsi280"; break;
    case ( kNuisFSI281 ) : return "f_fsi281"; break;
    case ( kNuisFSI282 ) : return "f_fsi282"; break;
    case ( kNuisFSI283 ) : return "f_fsi283"; break;
    case ( kNuisFSI284 ) : return "f_fsi284"; break;
    case ( kNuisFSI285 ) : return "f_fsi285"; break;
    case ( kNuisFSI286 ) : return "f_fsi286"; break;
    case ( kNuisFSI287 ) : return "f_fsi287"; break;
    case ( kNuisFSI288 ) : return "f_fsi288"; break;
    case ( kNuisFSI289 ) : return "f_fsi289"; break;
    case ( kNuisFSI290 ) : return "f_fsi290"; break;
    case ( kNuisFSI291 ) : return "f_fsi291"; break;
    case ( kNuisFSI292 ) : return "f_fsi292"; break;
    case ( kNuisFSI293 ) : return "f_fsi293"; break;
    case ( kNuisFSI294 ) : return "f_fsi294"; break;
    case ( kNuisFSI295 ) : return "f_fsi295"; break;
    case ( kNuisFSI296 ) : return "f_fsi296"; break;
    case ( kNuisFSI297 ) : return "f_fsi297"; break;
    case ( kNuisFSI298 ) : return "f_fsi298"; break;
    case ( kNuisFSI299 ) : return "f_fsi299"; break;
    case ( kNuisFSI300 ) : return "f_fsi300"; break;
    case ( kNuisFSI301 ) : return "f_fsi301"; break;
    case ( kNuisFSI302 ) : return "f_fsi302"; break;
    case ( kNuisFSI303 ) : return "f_fsi303"; break;
    case ( kNuisFSI304 ) : return "f_fsi304"; break;
    case ( kNuisFSI305 ) : return "f_fsi305"; break;
    case ( kNuisFSI306 ) : return "f_fsi306"; break;
    case ( kNuisFSI307 ) : return "f_fsi307"; break;
    case ( kNuisFSI308 ) : return "f_fsi308"; break;
    case ( kNuisFSI309 ) : return "f_fsi309"; break;
    case ( kNuisFSI310 ) : return "f_fsi310"; break;
    case ( kNuisFSI311 ) : return "f_fsi311"; break;
    case ( kNuisFSI312 ) : return "f_fsi312"; break;
    case ( kNuisFSI313 ) : return "f_fsi313"; break;
    case ( kNuisFSI314 ) : return "f_fsi314"; break;
    case ( kNuisFSI315 ) : return "f_fsi315"; break;
    case ( kNuisFSI316 ) : return "f_fsi316"; break;
    case ( kNuisFSI317 ) : return "f_fsi317"; break;
    case ( kNuisFSI318 ) : return "f_fsi318"; break;
    case ( kNuisFSI319 ) : return "f_fsi319"; break;
    case ( kNuisFSI320 ) : return "f_fsi320"; break;
    case ( kNuisFSI321 ) : return "f_fsi321"; break;
    case ( kNuisFSI322 ) : return "f_fsi322"; break;
    case ( kNuisFSI323 ) : return "f_fsi323"; break;
    case ( kNuisFSI324 ) : return "f_fsi324"; break;
    case ( kNuisFSI325 ) : return "f_fsi325"; break;
    case ( kNuisFSI326 ) : return "f_fsi326"; break;
    case ( kNuisFSI327 ) : return "f_fsi327"; break;
    case ( kNuisFSI328 ) : return "f_fsi328"; break;
    case ( kNuisFSI329 ) : return "f_fsi329"; break;
    case ( kNuisFSI330 ) : return "f_fsi330"; break;
    case ( kNuisFSI331 ) : return "f_fsi331"; break;
    case ( kNuisFSI332 ) : return "f_fsi332"; break;
    case ( kNuisFSI333 ) : return "f_fsi333"; break;
    case ( kNuisFSI334 ) : return "f_fsi334"; break;
    case ( kNuisFSI335 ) : return "f_fsi335"; break;
    case ( kNuisFSI336 ) : return "f_fsi336"; break;
    case ( kNuisFSI337 ) : return "f_fsi337"; break;
    case ( kNuisFSI338 ) : return "f_fsi338"; break;
    case ( kNuisFSI339 ) : return "f_fsi339"; break;
    case ( kNuisFSI340 ) : return "f_fsi340"; break;
    case ( kNuisFSI341 ) : return "f_fsi341"; break;
    case ( kNuisFSI342 ) : return "f_fsi342"; break;
    case ( kNuisFSI343 ) : return "f_fsi343"; break;
    case ( kNuisFSI344 ) : return "f_fsi344"; break;
    case ( kNuisFSI345 ) : return "f_fsi345"; break;
    case ( kNuisFSI346 ) : return "f_fsi346"; break;
    case ( kNuisFSI347 ) : return "f_fsi347"; break;
    case ( kNuisFSI348 ) : return "f_fsi348"; break;
    case ( kNuisFSI349 ) : return "f_fsi349"; break;
    case ( kNuisFSI350 ) : return "f_fsi350"; break;
    case ( kNuisFSI351 ) : return "f_fsi351"; break;
    case ( kNuisFSI352 ) : return "f_fsi352"; break;
    case ( kNuisFSI353 ) : return "f_fsi353"; break;
    case ( kNuisFSI354 ) : return "f_fsi354"; break;
    case ( kNuisFSI355 ) : return "f_fsi355"; break;
    case ( kNuisFSI356 ) : return "f_fsi356"; break;
    case ( kNuisFSI357 ) : return "f_fsi357"; break;
    case ( kNuisFSI358 ) : return "f_fsi358"; break;
    case ( kNuisFSI359 ) : return "f_fsi359"; break;
    case ( kNuisFSI360 ) : return "f_fsi360"; break;
    case ( kNuisFSI361 ) : return "f_fsi361"; break;
    case ( kNuisFSI362 ) : return "f_fsi362"; break;
    case ( kNuisFSI363 ) : return "f_fsi363"; break;
    case ( kNuisFSI364 ) : return "f_fsi364"; break;
    case ( kNuisFSI365 ) : return "f_fsi365"; break;
    case ( kNuisFSI366 ) : return "f_fsi366"; break;
    case ( kNuisFSI367 ) : return "f_fsi367"; break;
    case ( kNuisFSI368 ) : return "f_fsi368"; break;
    case ( kNuisFSI369 ) : return "f_fsi369"; break;
    case ( kNuisFSI370 ) : return "f_fsi370"; break;
    case ( kNuisFSI371 ) : return "f_fsi371"; break;
    case ( kNuisFSI372 ) : return "f_fsi372"; break;
    case ( kNuisFSI373 ) : return "f_fsi373"; break;
    case ( kNuisFSI374 ) : return "f_fsi374"; break;
    case ( kNuisFSI375 ) : return "f_fsi375"; break;
    case ( kNuisFSI376 ) : return "f_fsi376"; break;
    case ( kNuisFSI377 ) : return "f_fsi377"; break;
    case ( kNuisFSI378 ) : return "f_fsi378"; break;
    case ( kNuisFSI379 ) : return "f_fsi379"; break;
    case ( kNuisFSI380 ) : return "f_fsi380"; break;
    case ( kNuisFSI381 ) : return "f_fsi381"; break;
    case ( kNuisFSI382 ) : return "f_fsi382"; break;
    case ( kNuisFSI383 ) : return "f_fsi383"; break;
    case ( kNuisFSI384 ) : return "f_fsi384"; break;
    case ( kNuisFSI385 ) : return "f_fsi385"; break;
    case ( kNuisFSI386 ) : return "f_fsi386"; break;
    case ( kNuisFSI387 ) : return "f_fsi387"; break;
    case ( kNuisFSI388 ) : return "f_fsi388"; break;
    case ( kNuisFSI389 ) : return "f_fsi389"; break;
    case ( kNuisFSI390 ) : return "f_fsi390"; break;
    case ( kNuisFSI391 ) : return "f_fsi391"; break;
    case ( kNuisFSI392 ) : return "f_fsi392"; break;
    case ( kNuisFSI393 ) : return "f_fsi393"; break;
    case ( kNuisFSI394 ) : return "f_fsi394"; break;
    case ( kNuisFSI395 ) : return "f_fsi395"; break;
    case ( kNuisFSI396 ) : return "f_fsi396"; break;
    case ( kNuisFSI397 ) : return "f_fsi397"; break;
    case ( kNuisFSI398 ) : return "f_fsi398"; break;
    case ( kNuisFSI399 ) : return "f_fsi399"; break;
    case ( kNuisFSI400 ) : return "f_fsi400"; break;
    case ( kNuisFSI401 ) : return "f_fsi401"; break;
    case ( kNuisFSI402 ) : return "f_fsi402"; break;
    case ( kNuisFSI403 ) : return "f_fsi403"; break;
    case ( kNuisFSI404 ) : return "f_fsi404"; break;
    case ( kNuisFSI405 ) : return "f_fsi405"; break;
    case ( kNuisFSI406 ) : return "f_fsi406"; break;
    case ( kNuisFSI407 ) : return "f_fsi407"; break;
    case ( kNuisFSI408 ) : return "f_fsi408"; break;
    case ( kNuisFSI409 ) : return "f_fsi409"; break;
    case ( kNuisFSI410 ) : return "f_fsi410"; break;
    case ( kNuisFSI411 ) : return "f_fsi411"; break;
    case ( kNuisFSI412 ) : return "f_fsi412"; break;
    case ( kNuisFSI413 ) : return "f_fsi413"; break;
    case ( kNuisFSI414 ) : return "f_fsi414"; break;
    case ( kNuisFSI415 ) : return "f_fsi415"; break;
    case ( kNuisFSI416 ) : return "f_fsi416"; break;
    case ( kNuisFSI417 ) : return "f_fsi417"; break;
    case ( kNuisFSI418 ) : return "f_fsi418"; break;
    case ( kNuisFSI419 ) : return "f_fsi419"; break;
    case ( kNuisFSI420 ) : return "f_fsi420"; break;
    case ( kNuisFSI421 ) : return "f_fsi421"; break;
    case ( kNuisFSI422 ) : return "f_fsi422"; break;
    case ( kNuisFSI423 ) : return "f_fsi423"; break;
    case ( kNuisFSI424 ) : return "f_fsi424"; break;
    case ( kNuisFSI425 ) : return "f_fsi425"; break;
    case ( kNuisFSI426 ) : return "f_fsi426"; break;
    case ( kNuisFSI427 ) : return "f_fsi427"; break;
    case ( kNuisFSI428 ) : return "f_fsi428"; break;
    case ( kNuisFSI429 ) : return "f_fsi429"; break;
    case ( kNuisFSI430 ) : return "f_fsi430"; break;
    case ( kNuisFSI431 ) : return "f_fsi431"; break;
    case ( kNuisFSI432 ) : return "f_fsi432"; break;
    case ( kNuisFSI433 ) : return "f_fsi433"; break;
    case ( kNuisFSI434 ) : return "f_fsi434"; break;
    case ( kNuisFSI435 ) : return "f_fsi435"; break;
    case ( kNuisFSI436 ) : return "f_fsi436"; break;
    case ( kNuisFSI437 ) : return "f_fsi437"; break;
    case ( kNuisFSI438 ) : return "f_fsi438"; break;
    case ( kNuisFSI439 ) : return "f_fsi439"; break;
    case ( kNuisFSI440 ) : return "f_fsi440"; break;
    case ( kNuisFSI441 ) : return "f_fsi441"; break;
    case ( kNuisFSI442 ) : return "f_fsi442"; break;
    case ( kNuisFSI443 ) : return "f_fsi443"; break;
    case ( kNuisFSI444 ) : return "f_fsi444"; break;
    case ( kNuisFSI445 ) : return "f_fsi445"; break;
    case ( kNuisFSI446 ) : return "f_fsi446"; break;
    case ( kNuisFSI447 ) : return "f_fsi447"; break;
    case ( kNuisFSI448 ) : return "f_fsi448"; break;
    case ( kNuisFSI449 ) : return "f_fsi449"; break;
    case ( kNuisFSI450 ) : return "f_fsi450"; break;
    case ( kNuisFSI451 ) : return "f_fsi451"; break;
    case ( kNuisFSI452 ) : return "f_fsi452"; break;
    case ( kNuisFSI453 ) : return "f_fsi453"; break;
    case ( kNuisFSI454 ) : return "f_fsi454"; break;
    case ( kNuisFSI455 ) : return "f_fsi455"; break;
    case ( kNuisFSI456 ) : return "f_fsi456"; break;
    case ( kNuisFSI457 ) : return "f_fsi457"; break;
    case ( kNuisFSI458 ) : return "f_fsi458"; break;
    case ( kNuisFSI459 ) : return "f_fsi459"; break;
    case ( kNuisFSI460 ) : return "f_fsi460"; break;
    case ( kNuisFSI461 ) : return "f_fsi461"; break;
    case ( kNuisFSI462 ) : return "f_fsi462"; break;
    case ( kNuisFSI463 ) : return "f_fsi463"; break;
    case ( kNuisFSI464 ) : return "f_fsi464"; break;
    case ( kNuisFSI465 ) : return "f_fsi465"; break;
    case ( kNuisFSI466 ) : return "f_fsi466"; break;
    case ( kNuisFSI467 ) : return "f_fsi467"; break;
    case ( kNuisFSI468 ) : return "f_fsi468"; break;
    case ( kNuisFSI469 ) : return "f_fsi469"; break;
    case ( kNuisFSI470 ) : return "f_fsi470"; break;
    case ( kNuisFSI471 ) : return "f_fsi471"; break;
    case ( kNuisFSI472 ) : return "f_fsi472"; break;
    case ( kNuisFSI473 ) : return "f_fsi473"; break;
    case ( kNuisFSI474 ) : return "f_fsi474"; break;
    case ( kNuisFSI475 ) : return "f_fsi475"; break;
    case ( kNuisFSI476 ) : return "f_fsi476"; break;
    case ( kNuisFSI477 ) : return "f_fsi477"; break;
    case ( kNuisFSI478 ) : return "f_fsi478"; break;
    case ( kNuisFSI479 ) : return "f_fsi479"; break;
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
     case ( kSKmSigNuECCQE      ) : return "appearance nue CC QE";                  break;
     case ( kSKmSigNuECC1pi     ) : return "appearance nue CC 1pi";                 break;
     case ( kSKmSigNuECCcoh     ) : return "appearance nue CC coh";                 break;
     case ( kSKmSigNuECCoth_mu  ) : return "appearance nue CC oth (numu disap)";    break;
     case ( kSKmSigNuECCoth_e   ) : return "appearance nue CC oth (nue app)";       break;
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
     case ( kSKmSigNuECCQE      ) : return "appearance $\\nu_e$ CCQE";                                           break;
     case ( kSKmSigNuECC1pi     ) : return "appearance $\\nu_e$ CC 1$\\pi$";                                     break;
     case ( kSKmSigNuECCcoh     ) : return "appearance $\\nu_e$ CC coherent";                                    break;
     case ( kSKmSigNuECCoth_mu  ) : return "appearance $\\nu_e$ CC other ($\\nu_\\mu$ app)";                     break;
     case ( kSKmSigNuECCoth_e   ) : return "appearance $\\nu_e$ CC other ($\\nu_e$ app)";                        break;
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
     case ( kSKmSigNuECCQE      ) : return "appearance #nu_{e} CCQE";                           break;
     case ( kSKmSigNuECC1pi     ) : return "appearance #nu_{e} CC 1#pi";                        break;
     case ( kSKmSigNuECCcoh     ) : return "appearance #nu_{e} CC coherent";                    break;
     case ( kSKmSigNuECCoth_mu  ) : return "appearance #nu_{e} CC other";                       break;
     case ( kSKmSigNuECCoth_e   ) : return "appearance #nu_{e} CC other";                       break;
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
  if(sh == "banff_maqe")    return kNuisBANFFMAQE;
  if(sh == "banff_mares")   return kNuisBANFFMARES;
  if(sh == "banff_ccqeE0")  return kNuisBANFFCCQEE0;
#if BANFF_BINS == 2
  if(sh == "banff_ccqeE1")  return kNuisBANFFCCQEE1;
  if(sh == "banff_ccqeE2")  return kNuisBANFFCCQEE2;
#endif
  if(sh == "banff_cc1piE0") return kNuisBANFFCC1piE0;
#if BANFF_BINS == 2
  if(sh == "banff_cc1piE1") return kNuisBANFFCC1piE1;
#endif
  if(sh == "banff_nc1pi0")  return kNuisBANFFNC1pi0;

  if(sh == "xsec_pf")           return kNuisXSecFermiMom;
  if(sh == "xsec_sf")           return kNuisXSecSpecFunc;
  if(sh == "xsec_wshape")       return kNuisXSecWshape;
  if(sh == "xsec_shape_ccoth")  return kNuisXSecShapeCCoth;
  if(sh == "xsec_pilessdelta")  return kNuisXSecPiLessDelta;
  if(sh == "xsec_enushape")     return kNuisXSecENuShape;
#if BANFF_BINS != 2
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
  if(sh == "fsi_049") return kNuisFSI049;
  if(sh == "fsi_050") return kNuisFSI050;
  if(sh == "fsi_051") return kNuisFSI051;
  if(sh == "fsi_052") return kNuisFSI052;
  if(sh == "fsi_053") return kNuisFSI053;
  if(sh == "fsi_054") return kNuisFSI054;
  if(sh == "fsi_055") return kNuisFSI055;
  if(sh == "fsi_056") return kNuisFSI056;
  if(sh == "fsi_057") return kNuisFSI057;
  if(sh == "fsi_058") return kNuisFSI058;
  if(sh == "fsi_059") return kNuisFSI059;
  if(sh == "fsi_060") return kNuisFSI060;
  if(sh == "fsi_061") return kNuisFSI061;
  if(sh == "fsi_062") return kNuisFSI062;
  if(sh == "fsi_063") return kNuisFSI063;
  if(sh == "fsi_064") return kNuisFSI064;
  if(sh == "fsi_065") return kNuisFSI065;
  if(sh == "fsi_066") return kNuisFSI066;
  if(sh == "fsi_067") return kNuisFSI067;
  if(sh == "fsi_068") return kNuisFSI068;
  if(sh == "fsi_069") return kNuisFSI069;
  if(sh == "fsi_070") return kNuisFSI070;
  if(sh == "fsi_071") return kNuisFSI071;
  if(sh == "fsi_072") return kNuisFSI072;
  if(sh == "fsi_073") return kNuisFSI073;
  if(sh == "fsi_074") return kNuisFSI074;
  if(sh == "fsi_075") return kNuisFSI075;
  if(sh == "fsi_076") return kNuisFSI076;
  if(sh == "fsi_077") return kNuisFSI077;
  if(sh == "fsi_078") return kNuisFSI078;
#if FSI_BINS >= 4
  if(sh == "fsi_079") return kNuisFSI079;
  if(sh == "fsi_080") return kNuisFSI080;
  if(sh == "fsi_081") return kNuisFSI081;
  if(sh == "fsi_082") return kNuisFSI082;
  if(sh == "fsi_083") return kNuisFSI083;
  if(sh == "fsi_084") return kNuisFSI084;
  if(sh == "fsi_085") return kNuisFSI085;
  if(sh == "fsi_086") return kNuisFSI086;
  if(sh == "fsi_087") return kNuisFSI087;
  if(sh == "fsi_088") return kNuisFSI088;
  if(sh == "fsi_089") return kNuisFSI089;
  if(sh == "fsi_090") return kNuisFSI090;
  if(sh == "fsi_091") return kNuisFSI091;
  if(sh == "fsi_092") return kNuisFSI092;
  if(sh == "fsi_093") return kNuisFSI093;
  if(sh == "fsi_094") return kNuisFSI094;
  if(sh == "fsi_095") return kNuisFSI095;
  if(sh == "fsi_096") return kNuisFSI096;
  if(sh == "fsi_097") return kNuisFSI097;
  if(sh == "fsi_098") return kNuisFSI098;
  if(sh == "fsi_099") return kNuisFSI099;
  if(sh == "fsi_100") return kNuisFSI100;
  if(sh == "fsi_101") return kNuisFSI101;
  if(sh == "fsi_102") return kNuisFSI102;
  if(sh == "fsi_103") return kNuisFSI103;
  if(sh == "fsi_104") return kNuisFSI104;
  if(sh == "fsi_105") return kNuisFSI105;
  if(sh == "fsi_106") return kNuisFSI106;
  if(sh == "fsi_107") return kNuisFSI107;
  if(sh == "fsi_108") return kNuisFSI108;
  if(sh == "fsi_109") return kNuisFSI109;
  if(sh == "fsi_110") return kNuisFSI110;
  if(sh == "fsi_111") return kNuisFSI111;
  if(sh == "fsi_112") return kNuisFSI112;
  if(sh == "fsi_113") return kNuisFSI113;
  if(sh == "fsi_114") return kNuisFSI114;
  if(sh == "fsi_115") return kNuisFSI115;
  if(sh == "fsi_116") return kNuisFSI116;
  if(sh == "fsi_117") return kNuisFSI117;
  if(sh == "fsi_118") return kNuisFSI118;
  if(sh == "fsi_119") return kNuisFSI119;
  if(sh == "fsi_120") return kNuisFSI120;
  if(sh == "fsi_121") return kNuisFSI121;
  if(sh == "fsi_122") return kNuisFSI122;
  if(sh == "fsi_123") return kNuisFSI123;
  if(sh == "fsi_124") return kNuisFSI124;
  if(sh == "fsi_125") return kNuisFSI125;
  if(sh == "fsi_126") return kNuisFSI126;
  if(sh == "fsi_127") return kNuisFSI127;
  if(sh == "fsi_128") return kNuisFSI128;
  if(sh == "fsi_129") return kNuisFSI129;
  if(sh == "fsi_130") return kNuisFSI130;
  if(sh == "fsi_131") return kNuisFSI131;
  if(sh == "fsi_132") return kNuisFSI132;
  if(sh == "fsi_133") return kNuisFSI133;
  if(sh == "fsi_134") return kNuisFSI134;
  if(sh == "fsi_135") return kNuisFSI135;
  if(sh == "fsi_136") return kNuisFSI136;
  if(sh == "fsi_137") return kNuisFSI137;
  if(sh == "fsi_138") return kNuisFSI138;
  if(sh == "fsi_139") return kNuisFSI139;
  if(sh == "fsi_140") return kNuisFSI140;
  if(sh == "fsi_141") return kNuisFSI141;
  if(sh == "fsi_142") return kNuisFSI142;
  if(sh == "fsi_143") return kNuisFSI143;
  if(sh == "fsi_144") return kNuisFSI144;
  if(sh == "fsi_145") return kNuisFSI145;
  if(sh == "fsi_146") return kNuisFSI146;
  if(sh == "fsi_147") return kNuisFSI147;
  if(sh == "fsi_148") return kNuisFSI148;
  if(sh == "fsi_149") return kNuisFSI149;
  if(sh == "fsi_150") return kNuisFSI150;
  if(sh == "fsi_151") return kNuisFSI151;
  if(sh == "fsi_152") return kNuisFSI152;
  if(sh == "fsi_153") return kNuisFSI153;
  if(sh == "fsi_154") return kNuisFSI154;
  if(sh == "fsi_155") return kNuisFSI155;
  if(sh == "fsi_156") return kNuisFSI156;
  if(sh == "fsi_157") return kNuisFSI157;
  if(sh == "fsi_158") return kNuisFSI158;
  if(sh == "fsi_159") return kNuisFSI159;
  if(sh == "fsi_160") return kNuisFSI160;
  if(sh == "fsi_161") return kNuisFSI161;
  if(sh == "fsi_162") return kNuisFSI162;
  if(sh == "fsi_163") return kNuisFSI163;
  if(sh == "fsi_164") return kNuisFSI164;
  if(sh == "fsi_165") return kNuisFSI165;
  if(sh == "fsi_166") return kNuisFSI166;
  if(sh == "fsi_167") return kNuisFSI167;
  if(sh == "fsi_168") return kNuisFSI168;
  if(sh == "fsi_169") return kNuisFSI169;
  if(sh == "fsi_170") return kNuisFSI170;
  if(sh == "fsi_171") return kNuisFSI171;
  if(sh == "fsi_172") return kNuisFSI172;
  if(sh == "fsi_173") return kNuisFSI173;
  if(sh == "fsi_174") return kNuisFSI174;
  if(sh == "fsi_175") return kNuisFSI175;
  if(sh == "fsi_176") return kNuisFSI176;
  if(sh == "fsi_177") return kNuisFSI177;
  if(sh == "fsi_178") return kNuisFSI178;
  if(sh == "fsi_179") return kNuisFSI179;
  if(sh == "fsi_180") return kNuisFSI180;
  if(sh == "fsi_181") return kNuisFSI181;
  if(sh == "fsi_182") return kNuisFSI182;
  if(sh == "fsi_183") return kNuisFSI183;
  if(sh == "fsi_184") return kNuisFSI184;
  if(sh == "fsi_185") return kNuisFSI185;
  if(sh == "fsi_186") return kNuisFSI186;
  if(sh == "fsi_187") return kNuisFSI187;
  if(sh == "fsi_188") return kNuisFSI188;
  if(sh == "fsi_189") return kNuisFSI189;
  if(sh == "fsi_190") return kNuisFSI190;
  if(sh == "fsi_191") return kNuisFSI191;
  if(sh == "fsi_192") return kNuisFSI192;
  if(sh == "fsi_193") return kNuisFSI193;
  if(sh == "fsi_194") return kNuisFSI194;
  if(sh == "fsi_195") return kNuisFSI195;
  if(sh == "fsi_196") return kNuisFSI196;
  if(sh == "fsi_197") return kNuisFSI197;
  if(sh == "fsi_198") return kNuisFSI198;
  if(sh == "fsi_199") return kNuisFSI199;
  if(sh == "fsi_200") return kNuisFSI200;
  if(sh == "fsi_201") return kNuisFSI201;
  if(sh == "fsi_202") return kNuisFSI202;
  if(sh == "fsi_203") return kNuisFSI203;
  if(sh == "fsi_204") return kNuisFSI204;
  if(sh == "fsi_205") return kNuisFSI205;
  if(sh == "fsi_206") return kNuisFSI206;
  if(sh == "fsi_207") return kNuisFSI207;
  if(sh == "fsi_208") return kNuisFSI208;
  if(sh == "fsi_209") return kNuisFSI209;
  if(sh == "fsi_210") return kNuisFSI210;
  if(sh == "fsi_211") return kNuisFSI211;
  if(sh == "fsi_212") return kNuisFSI212;
  if(sh == "fsi_213") return kNuisFSI213;
  if(sh == "fsi_214") return kNuisFSI214;
  if(sh == "fsi_215") return kNuisFSI215;
  if(sh == "fsi_216") return kNuisFSI216;
  if(sh == "fsi_217") return kNuisFSI217;
  if(sh == "fsi_218") return kNuisFSI218;
  if(sh == "fsi_219") return kNuisFSI219;
  if(sh == "fsi_220") return kNuisFSI220;
  if(sh == "fsi_221") return kNuisFSI221;
  if(sh == "fsi_222") return kNuisFSI222;
  if(sh == "fsi_223") return kNuisFSI223;
  if(sh == "fsi_224") return kNuisFSI224;
  if(sh == "fsi_225") return kNuisFSI225;
  if(sh == "fsi_226") return kNuisFSI226;
  if(sh == "fsi_227") return kNuisFSI227;
  if(sh == "fsi_228") return kNuisFSI228;
  if(sh == "fsi_229") return kNuisFSI229;
  if(sh == "fsi_230") return kNuisFSI230;
  if(sh == "fsi_231") return kNuisFSI231;
  if(sh == "fsi_232") return kNuisFSI232;
  if(sh == "fsi_233") return kNuisFSI233;
  if(sh == "fsi_234") return kNuisFSI234;
  if(sh == "fsi_235") return kNuisFSI235;
  if(sh == "fsi_236") return kNuisFSI236;
  if(sh == "fsi_237") return kNuisFSI237;
  if(sh == "fsi_238") return kNuisFSI238;
  if(sh == "fsi_239") return kNuisFSI239;
  if(sh == "fsi_240") return kNuisFSI240;
  if(sh == "fsi_241") return kNuisFSI241;
  if(sh == "fsi_242") return kNuisFSI242;
  if(sh == "fsi_243") return kNuisFSI243;
  if(sh == "fsi_244") return kNuisFSI244;
  if(sh == "fsi_245") return kNuisFSI245;
  if(sh == "fsi_246") return kNuisFSI246;
  if(sh == "fsi_247") return kNuisFSI247;
  if(sh == "fsi_248") return kNuisFSI248;
  if(sh == "fsi_249") return kNuisFSI249;
  if(sh == "fsi_250") return kNuisFSI250;
  if(sh == "fsi_251") return kNuisFSI251;
  if(sh == "fsi_252") return kNuisFSI252;
  if(sh == "fsi_253") return kNuisFSI253;
  if(sh == "fsi_254") return kNuisFSI254;
  if(sh == "fsi_255") return kNuisFSI255;
  if(sh == "fsi_256") return kNuisFSI256;
  if(sh == "fsi_257") return kNuisFSI257;
  if(sh == "fsi_258") return kNuisFSI258;
  if(sh == "fsi_259") return kNuisFSI259;
  if(sh == "fsi_260") return kNuisFSI260;
  if(sh == "fsi_261") return kNuisFSI261;
  if(sh == "fsi_262") return kNuisFSI262;
  if(sh == "fsi_263") return kNuisFSI263;
  if(sh == "fsi_264") return kNuisFSI264;
  if(sh == "fsi_265") return kNuisFSI265;
  if(sh == "fsi_266") return kNuisFSI266;
  if(sh == "fsi_267") return kNuisFSI267;
  if(sh == "fsi_268") return kNuisFSI268;
  if(sh == "fsi_269") return kNuisFSI269;
  if(sh == "fsi_270") return kNuisFSI270;
  if(sh == "fsi_271") return kNuisFSI271;
  if(sh == "fsi_272") return kNuisFSI272;
  if(sh == "fsi_273") return kNuisFSI273;
  if(sh == "fsi_274") return kNuisFSI274;
  if(sh == "fsi_275") return kNuisFSI275;
  if(sh == "fsi_276") return kNuisFSI276;
  if(sh == "fsi_277") return kNuisFSI277;
  if(sh == "fsi_278") return kNuisFSI278;
  if(sh == "fsi_279") return kNuisFSI279;
  if(sh == "fsi_280") return kNuisFSI280;
  if(sh == "fsi_281") return kNuisFSI281;
  if(sh == "fsi_282") return kNuisFSI282;
  if(sh == "fsi_283") return kNuisFSI283;
  if(sh == "fsi_284") return kNuisFSI284;
  if(sh == "fsi_285") return kNuisFSI285;
  if(sh == "fsi_286") return kNuisFSI286;
  if(sh == "fsi_287") return kNuisFSI287;
  if(sh == "fsi_288") return kNuisFSI288;
  if(sh == "fsi_289") return kNuisFSI289;
  if(sh == "fsi_290") return kNuisFSI290;
  if(sh == "fsi_291") return kNuisFSI291;
  if(sh == "fsi_292") return kNuisFSI292;
  if(sh == "fsi_293") return kNuisFSI293;
  if(sh == "fsi_294") return kNuisFSI294;
  if(sh == "fsi_295") return kNuisFSI295;
  if(sh == "fsi_296") return kNuisFSI296;
  if(sh == "fsi_297") return kNuisFSI297;
  if(sh == "fsi_298") return kNuisFSI298;
  if(sh == "fsi_299") return kNuisFSI299;
  if(sh == "fsi_300") return kNuisFSI300;
  if(sh == "fsi_301") return kNuisFSI301;
  if(sh == "fsi_302") return kNuisFSI302;
  if(sh == "fsi_303") return kNuisFSI303;
  if(sh == "fsi_304") return kNuisFSI304;
  if(sh == "fsi_305") return kNuisFSI305;
  if(sh == "fsi_306") return kNuisFSI306;
  if(sh == "fsi_307") return kNuisFSI307;
  if(sh == "fsi_308") return kNuisFSI308;
  if(sh == "fsi_309") return kNuisFSI309;
  if(sh == "fsi_310") return kNuisFSI310;
  if(sh == "fsi_311") return kNuisFSI311;
  if(sh == "fsi_312") return kNuisFSI312;
  if(sh == "fsi_313") return kNuisFSI313;
  if(sh == "fsi_314") return kNuisFSI314;
  if(sh == "fsi_315") return kNuisFSI315;
  if(sh == "fsi_316") return kNuisFSI316;
  if(sh == "fsi_317") return kNuisFSI317;
  if(sh == "fsi_318") return kNuisFSI318;
  if(sh == "fsi_319") return kNuisFSI319;
  if(sh == "fsi_320") return kNuisFSI320;
  if(sh == "fsi_321") return kNuisFSI321;
  if(sh == "fsi_322") return kNuisFSI322;
  if(sh == "fsi_323") return kNuisFSI323;
  if(sh == "fsi_324") return kNuisFSI324;
  if(sh == "fsi_325") return kNuisFSI325;
  if(sh == "fsi_326") return kNuisFSI326;
  if(sh == "fsi_327") return kNuisFSI327;
  if(sh == "fsi_328") return kNuisFSI328;
  if(sh == "fsi_329") return kNuisFSI329;
  if(sh == "fsi_330") return kNuisFSI330;
  if(sh == "fsi_331") return kNuisFSI331;
  if(sh == "fsi_332") return kNuisFSI332;
  if(sh == "fsi_333") return kNuisFSI333;
  if(sh == "fsi_334") return kNuisFSI334;
  if(sh == "fsi_335") return kNuisFSI335;
  if(sh == "fsi_336") return kNuisFSI336;
  if(sh == "fsi_337") return kNuisFSI337;
  if(sh == "fsi_338") return kNuisFSI338;
  if(sh == "fsi_339") return kNuisFSI339;
  if(sh == "fsi_340") return kNuisFSI340;
  if(sh == "fsi_341") return kNuisFSI341;
  if(sh == "fsi_342") return kNuisFSI342;
  if(sh == "fsi_343") return kNuisFSI343;
  if(sh == "fsi_344") return kNuisFSI344;
  if(sh == "fsi_345") return kNuisFSI345;
  if(sh == "fsi_346") return kNuisFSI346;
  if(sh == "fsi_347") return kNuisFSI347;
  if(sh == "fsi_348") return kNuisFSI348;
  if(sh == "fsi_349") return kNuisFSI349;
  if(sh == "fsi_350") return kNuisFSI350;
  if(sh == "fsi_351") return kNuisFSI351;
  if(sh == "fsi_352") return kNuisFSI352;
  if(sh == "fsi_353") return kNuisFSI353;
  if(sh == "fsi_354") return kNuisFSI354;
  if(sh == "fsi_355") return kNuisFSI355;
  if(sh == "fsi_356") return kNuisFSI356;
  if(sh == "fsi_357") return kNuisFSI357;
  if(sh == "fsi_358") return kNuisFSI358;
  if(sh == "fsi_359") return kNuisFSI359;
  if(sh == "fsi_360") return kNuisFSI360;
  if(sh == "fsi_361") return kNuisFSI361;
  if(sh == "fsi_362") return kNuisFSI362;
  if(sh == "fsi_363") return kNuisFSI363;
  if(sh == "fsi_364") return kNuisFSI364;
  if(sh == "fsi_365") return kNuisFSI365;
  if(sh == "fsi_366") return kNuisFSI366;
  if(sh == "fsi_367") return kNuisFSI367;
  if(sh == "fsi_368") return kNuisFSI368;
  if(sh == "fsi_369") return kNuisFSI369;
  if(sh == "fsi_370") return kNuisFSI370;
  if(sh == "fsi_371") return kNuisFSI371;
  if(sh == "fsi_372") return kNuisFSI372;
  if(sh == "fsi_373") return kNuisFSI373;
  if(sh == "fsi_374") return kNuisFSI374;
  if(sh == "fsi_375") return kNuisFSI375;
  if(sh == "fsi_376") return kNuisFSI376;
  if(sh == "fsi_377") return kNuisFSI377;
  if(sh == "fsi_378") return kNuisFSI378;
  if(sh == "fsi_379") return kNuisFSI379;
  if(sh == "fsi_380") return kNuisFSI380;
  if(sh == "fsi_381") return kNuisFSI381;
  if(sh == "fsi_382") return kNuisFSI382;
  if(sh == "fsi_383") return kNuisFSI383;
  if(sh == "fsi_384") return kNuisFSI384;
  if(sh == "fsi_385") return kNuisFSI385;
  if(sh == "fsi_386") return kNuisFSI386;
  if(sh == "fsi_387") return kNuisFSI387;
  if(sh == "fsi_388") return kNuisFSI388;
  if(sh == "fsi_389") return kNuisFSI389;
  if(sh == "fsi_390") return kNuisFSI390;
  if(sh == "fsi_391") return kNuisFSI391;
  if(sh == "fsi_392") return kNuisFSI392;
  if(sh == "fsi_393") return kNuisFSI393;
  if(sh == "fsi_394") return kNuisFSI394;
  if(sh == "fsi_395") return kNuisFSI395;
  if(sh == "fsi_396") return kNuisFSI396;
  if(sh == "fsi_397") return kNuisFSI397;
  if(sh == "fsi_398") return kNuisFSI398;
  if(sh == "fsi_399") return kNuisFSI399;
  if(sh == "fsi_400") return kNuisFSI400;
  if(sh == "fsi_401") return kNuisFSI401;
  if(sh == "fsi_402") return kNuisFSI402;
  if(sh == "fsi_403") return kNuisFSI403;
  if(sh == "fsi_404") return kNuisFSI404;
  if(sh == "fsi_405") return kNuisFSI405;
  if(sh == "fsi_406") return kNuisFSI406;
  if(sh == "fsi_407") return kNuisFSI407;
  if(sh == "fsi_408") return kNuisFSI408;
  if(sh == "fsi_409") return kNuisFSI409;
  if(sh == "fsi_410") return kNuisFSI410;
  if(sh == "fsi_411") return kNuisFSI411;
  if(sh == "fsi_412") return kNuisFSI412;
  if(sh == "fsi_413") return kNuisFSI413;
  if(sh == "fsi_414") return kNuisFSI414;
  if(sh == "fsi_415") return kNuisFSI415;
  if(sh == "fsi_416") return kNuisFSI416;
  if(sh == "fsi_417") return kNuisFSI417;
  if(sh == "fsi_418") return kNuisFSI418;
  if(sh == "fsi_419") return kNuisFSI419;
  if(sh == "fsi_420") return kNuisFSI420;
  if(sh == "fsi_421") return kNuisFSI421;
  if(sh == "fsi_422") return kNuisFSI422;
  if(sh == "fsi_423") return kNuisFSI423;
  if(sh == "fsi_424") return kNuisFSI424;
  if(sh == "fsi_425") return kNuisFSI425;
  if(sh == "fsi_426") return kNuisFSI426;
  if(sh == "fsi_427") return kNuisFSI427;
  if(sh == "fsi_428") return kNuisFSI428;
  if(sh == "fsi_429") return kNuisFSI429;
  if(sh == "fsi_430") return kNuisFSI430;
  if(sh == "fsi_431") return kNuisFSI431;
  if(sh == "fsi_432") return kNuisFSI432;
  if(sh == "fsi_433") return kNuisFSI433;
  if(sh == "fsi_434") return kNuisFSI434;
  if(sh == "fsi_435") return kNuisFSI435;
  if(sh == "fsi_436") return kNuisFSI436;
  if(sh == "fsi_437") return kNuisFSI437;
  if(sh == "fsi_438") return kNuisFSI438;
  if(sh == "fsi_439") return kNuisFSI439;
  if(sh == "fsi_440") return kNuisFSI440;
  if(sh == "fsi_441") return kNuisFSI441;
  if(sh == "fsi_442") return kNuisFSI442;
  if(sh == "fsi_443") return kNuisFSI443;
  if(sh == "fsi_444") return kNuisFSI444;
  if(sh == "fsi_445") return kNuisFSI445;
  if(sh == "fsi_446") return kNuisFSI446;
  if(sh == "fsi_447") return kNuisFSI447;
  if(sh == "fsi_448") return kNuisFSI448;
  if(sh == "fsi_449") return kNuisFSI449;
  if(sh == "fsi_450") return kNuisFSI450;
  if(sh == "fsi_451") return kNuisFSI451;
  if(sh == "fsi_452") return kNuisFSI452;
  if(sh == "fsi_453") return kNuisFSI453;
  if(sh == "fsi_454") return kNuisFSI454;
  if(sh == "fsi_455") return kNuisFSI455;
  if(sh == "fsi_456") return kNuisFSI456;
  if(sh == "fsi_457") return kNuisFSI457;
  if(sh == "fsi_458") return kNuisFSI458;
  if(sh == "fsi_459") return kNuisFSI459;
  if(sh == "fsi_460") return kNuisFSI460;
  if(sh == "fsi_461") return kNuisFSI461;
  if(sh == "fsi_462") return kNuisFSI462;
  if(sh == "fsi_463") return kNuisFSI463;
  if(sh == "fsi_464") return kNuisFSI464;
  if(sh == "fsi_465") return kNuisFSI465;
  if(sh == "fsi_466") return kNuisFSI466;
  if(sh == "fsi_467") return kNuisFSI467;
  if(sh == "fsi_468") return kNuisFSI468;
  if(sh == "fsi_469") return kNuisFSI469;
  if(sh == "fsi_470") return kNuisFSI470;
  if(sh == "fsi_471") return kNuisFSI471;
  if(sh == "fsi_472") return kNuisFSI472;
  if(sh == "fsi_473") return kNuisFSI473;
  if(sh == "fsi_474") return kNuisFSI474;
  if(sh == "fsi_475") return kNuisFSI475;
  if(sh == "fsi_476") return kNuisFSI476;
  if(sh == "fsi_477") return kNuisFSI477;
  if(sh == "fsi_478") return kNuisFSI478;
  if(sh == "fsi_479") return kNuisFSI479;
#endif
#endif
#endif

  cout << sh << " not yet supported in util::PhysParamIdFromString()" << endl;
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
    case ( kNuisBANFFMAQE             ) : return true;  break;
    case ( kNuisBANFFMARES            ) : return true;  break;
    case ( kNuisBANFFCCQEE0           ) : return true;  break;
#if BANFF_BINS == 2
    case ( kNuisBANFFCCQEE1           ) : return true;  break;
    case ( kNuisBANFFCCQEE2           ) : return true;  break;
#endif
    case ( kNuisBANFFCC1piE0          ) : return true;  break;
#if BANFF_BINS == 2
    case ( kNuisBANFFCC1piE1          ) : return true;  break;
#endif
    case ( kNuisBANFFNC1pi0           ) : return true;  break;

    case ( kNuisXSecFermiMom          ) : return true;  break;
    case ( kNuisXSecSpecFunc          ) : return true;  break; 
    case ( kNuisXSecWshape            ) : return true;  break;
    case ( kNuisXSecShapeCCoth        ) : return true;  break;
    case ( kNuisXSecPiLessDelta       ) : return true;  break;
    case ( kNuisXSecENuShape          ) : return true;  break;

#if BANFF_BINS != 2
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
    case ( kNuisFSI049                ) : return true;  break;
    case ( kNuisFSI050                ) : return true;  break;
    case ( kNuisFSI051                ) : return true;  break;
    case ( kNuisFSI052                ) : return true;  break;
    case ( kNuisFSI053                ) : return true;  break;
    case ( kNuisFSI054                ) : return true;  break;
    case ( kNuisFSI055                ) : return true;  break;
    case ( kNuisFSI056                ) : return true;  break;
    case ( kNuisFSI057                ) : return true;  break;
    case ( kNuisFSI058                ) : return true;  break;
    case ( kNuisFSI059                ) : return true;  break;
    case ( kNuisFSI060                ) : return true;  break;
    case ( kNuisFSI061                ) : return true;  break;
    case ( kNuisFSI062                ) : return true;  break;
    case ( kNuisFSI063                ) : return true;  break;
    case ( kNuisFSI064                ) : return true;  break;
    case ( kNuisFSI065                ) : return true;  break;
    case ( kNuisFSI066                ) : return true;  break;
    case ( kNuisFSI067                ) : return true;  break;
    case ( kNuisFSI068                ) : return true;  break;
    case ( kNuisFSI069                ) : return true;  break;
    case ( kNuisFSI070                ) : return true;  break;
    case ( kNuisFSI071                ) : return true;  break;
    case ( kNuisFSI072                ) : return true;  break;
    case ( kNuisFSI073                ) : return true;  break;
    case ( kNuisFSI074                ) : return true;  break;
    case ( kNuisFSI075                ) : return true;  break;
    case ( kNuisFSI076                ) : return true;  break;
    case ( kNuisFSI077                ) : return true;  break;
    case ( kNuisFSI078                ) : return true;  break;
#if FSI_BINS >= 4
    case ( kNuisFSI079                ) : return true;  break;
    case ( kNuisFSI080                ) : return true;  break;
    case ( kNuisFSI081                ) : return true;  break;
    case ( kNuisFSI082                ) : return true;  break;
    case ( kNuisFSI083                ) : return true;  break;
    case ( kNuisFSI084                ) : return true;  break;
    case ( kNuisFSI085                ) : return true;  break;
    case ( kNuisFSI086                ) : return true;  break;
    case ( kNuisFSI087                ) : return true;  break;
    case ( kNuisFSI088                ) : return true;  break;
    case ( kNuisFSI089                ) : return true;  break;
    case ( kNuisFSI090                ) : return true;  break;
    case ( kNuisFSI091                ) : return true;  break;
    case ( kNuisFSI092                ) : return true;  break;
    case ( kNuisFSI093                ) : return true;  break;
    case ( kNuisFSI094                ) : return true;  break;
    case ( kNuisFSI095                ) : return true;  break;
    case ( kNuisFSI096                ) : return true;  break;
    case ( kNuisFSI097                ) : return true;  break;
    case ( kNuisFSI098                ) : return true;  break;
    case ( kNuisFSI099                ) : return true;  break;
    case ( kNuisFSI100                ) : return true;  break;
    case ( kNuisFSI101                ) : return true;  break;
    case ( kNuisFSI102                ) : return true;  break;
    case ( kNuisFSI103                ) : return true;  break;
    case ( kNuisFSI104                ) : return true;  break;
    case ( kNuisFSI105                ) : return true;  break;
    case ( kNuisFSI106                ) : return true;  break;
    case ( kNuisFSI107                ) : return true;  break;
    case ( kNuisFSI108                ) : return true;  break;
    case ( kNuisFSI109                ) : return true;  break;
    case ( kNuisFSI110                ) : return true;  break;
    case ( kNuisFSI111                ) : return true;  break;
    case ( kNuisFSI112                ) : return true;  break;
    case ( kNuisFSI113                ) : return true;  break;
    case ( kNuisFSI114                ) : return true;  break;
    case ( kNuisFSI115                ) : return true;  break;
    case ( kNuisFSI116                ) : return true;  break;
    case ( kNuisFSI117                ) : return true;  break;
    case ( kNuisFSI118                ) : return true;  break;
    case ( kNuisFSI119                ) : return true;  break;
    case ( kNuisFSI120                ) : return true;  break;
    case ( kNuisFSI121                ) : return true;  break;
    case ( kNuisFSI122                ) : return true;  break;
    case ( kNuisFSI123                ) : return true;  break;
    case ( kNuisFSI124                ) : return true;  break;
    case ( kNuisFSI125                ) : return true;  break;
    case ( kNuisFSI126                ) : return true;  break;
    case ( kNuisFSI127                ) : return true;  break;
    case ( kNuisFSI128                ) : return true;  break;
    case ( kNuisFSI129                ) : return true;  break;
    case ( kNuisFSI130                ) : return true;  break;
    case ( kNuisFSI131                ) : return true;  break;
    case ( kNuisFSI132                ) : return true;  break;
    case ( kNuisFSI133                ) : return true;  break;
    case ( kNuisFSI134                ) : return true;  break;
    case ( kNuisFSI135                ) : return true;  break;
    case ( kNuisFSI136                ) : return true;  break;
    case ( kNuisFSI137                ) : return true;  break;
    case ( kNuisFSI138                ) : return true;  break;
    case ( kNuisFSI139                ) : return true;  break;
    case ( kNuisFSI140                ) : return true;  break;
    case ( kNuisFSI141                ) : return true;  break;
    case ( kNuisFSI142                ) : return true;  break;
    case ( kNuisFSI143                ) : return true;  break;
    case ( kNuisFSI144                ) : return true;  break;
    case ( kNuisFSI145                ) : return true;  break;
    case ( kNuisFSI146                ) : return true;  break;
    case ( kNuisFSI147                ) : return true;  break;
    case ( kNuisFSI148                ) : return true;  break;
    case ( kNuisFSI149                ) : return true;  break;
    case ( kNuisFSI150                ) : return true;  break;
    case ( kNuisFSI151                ) : return true;  break;
    case ( kNuisFSI152                ) : return true;  break;
    case ( kNuisFSI153                ) : return true;  break;
    case ( kNuisFSI154                ) : return true;  break;
    case ( kNuisFSI155                ) : return true;  break;
    case ( kNuisFSI156                ) : return true;  break;
    case ( kNuisFSI157                ) : return true;  break;
    case ( kNuisFSI158                ) : return true;  break;
    case ( kNuisFSI159                ) : return true;  break;
    case ( kNuisFSI160                ) : return true;  break;
    case ( kNuisFSI161                ) : return true;  break;
    case ( kNuisFSI162                ) : return true;  break;
    case ( kNuisFSI163                ) : return true;  break;
    case ( kNuisFSI164                ) : return true;  break;
    case ( kNuisFSI165                ) : return true;  break;
    case ( kNuisFSI166                ) : return true;  break;
    case ( kNuisFSI167                ) : return true;  break;
    case ( kNuisFSI168                ) : return true;  break;
    case ( kNuisFSI169                ) : return true;  break;
    case ( kNuisFSI170                ) : return true;  break;
    case ( kNuisFSI171                ) : return true;  break;
    case ( kNuisFSI172                ) : return true;  break;
    case ( kNuisFSI173                ) : return true;  break;
    case ( kNuisFSI174                ) : return true;  break;
    case ( kNuisFSI175                ) : return true;  break;
    case ( kNuisFSI176                ) : return true;  break;
    case ( kNuisFSI177                ) : return true;  break;
    case ( kNuisFSI178                ) : return true;  break;
    case ( kNuisFSI179                ) : return true;  break;
    case ( kNuisFSI180                ) : return true;  break;
    case ( kNuisFSI181                ) : return true;  break;
    case ( kNuisFSI182                ) : return true;  break;
    case ( kNuisFSI183                ) : return true;  break;
    case ( kNuisFSI184                ) : return true;  break;
    case ( kNuisFSI185                ) : return true;  break;
    case ( kNuisFSI186                ) : return true;  break;
    case ( kNuisFSI187                ) : return true;  break;
    case ( kNuisFSI188                ) : return true;  break;
    case ( kNuisFSI189                ) : return true;  break;
    case ( kNuisFSI190                ) : return true;  break;
    case ( kNuisFSI191                ) : return true;  break;
    case ( kNuisFSI192                ) : return true;  break;
    case ( kNuisFSI193                ) : return true;  break;
    case ( kNuisFSI194                ) : return true;  break;
    case ( kNuisFSI195                ) : return true;  break;
    case ( kNuisFSI196                ) : return true;  break;
    case ( kNuisFSI197                ) : return true;  break;
    case ( kNuisFSI198                ) : return true;  break;
    case ( kNuisFSI199                ) : return true;  break;
    case ( kNuisFSI200                ) : return true;  break;
    case ( kNuisFSI201                ) : return true;  break;
    case ( kNuisFSI202                ) : return true;  break;
    case ( kNuisFSI203                ) : return true;  break;
    case ( kNuisFSI204                ) : return true;  break;
    case ( kNuisFSI205                ) : return true;  break;
    case ( kNuisFSI206                ) : return true;  break;
    case ( kNuisFSI207                ) : return true;  break;
    case ( kNuisFSI208                ) : return true;  break;
    case ( kNuisFSI209                ) : return true;  break;
    case ( kNuisFSI210                ) : return true;  break;
    case ( kNuisFSI211                ) : return true;  break;
    case ( kNuisFSI212                ) : return true;  break;
    case ( kNuisFSI213                ) : return true;  break;
    case ( kNuisFSI214                ) : return true;  break;
    case ( kNuisFSI215                ) : return true;  break;
    case ( kNuisFSI216                ) : return true;  break;
    case ( kNuisFSI217                ) : return true;  break;
    case ( kNuisFSI218                ) : return true;  break;
    case ( kNuisFSI219                ) : return true;  break;
    case ( kNuisFSI220                ) : return true;  break;
    case ( kNuisFSI221                ) : return true;  break;
    case ( kNuisFSI222                ) : return true;  break;
    case ( kNuisFSI223                ) : return true;  break;
    case ( kNuisFSI224                ) : return true;  break;
    case ( kNuisFSI225                ) : return true;  break;
    case ( kNuisFSI226                ) : return true;  break;
    case ( kNuisFSI227                ) : return true;  break;
    case ( kNuisFSI228                ) : return true;  break;
    case ( kNuisFSI229                ) : return true;  break;
    case ( kNuisFSI230                ) : return true;  break;
    case ( kNuisFSI231                ) : return true;  break;
    case ( kNuisFSI232                ) : return true;  break;
    case ( kNuisFSI233                ) : return true;  break;
    case ( kNuisFSI234                ) : return true;  break;
    case ( kNuisFSI235                ) : return true;  break;
    case ( kNuisFSI236                ) : return true;  break;
    case ( kNuisFSI237                ) : return true;  break;
    case ( kNuisFSI238                ) : return true;  break;
    case ( kNuisFSI239                ) : return true;  break;
    case ( kNuisFSI240                ) : return true;  break;
    case ( kNuisFSI241                ) : return true;  break;
    case ( kNuisFSI242                ) : return true;  break;
    case ( kNuisFSI243                ) : return true;  break;
    case ( kNuisFSI244                ) : return true;  break;
    case ( kNuisFSI245                ) : return true;  break;
    case ( kNuisFSI246                ) : return true;  break;
    case ( kNuisFSI247                ) : return true;  break;
    case ( kNuisFSI248                ) : return true;  break;
    case ( kNuisFSI249                ) : return true;  break;
    case ( kNuisFSI250                ) : return true;  break;
    case ( kNuisFSI251                ) : return true;  break;
    case ( kNuisFSI252                ) : return true;  break;
    case ( kNuisFSI253                ) : return true;  break;
    case ( kNuisFSI254                ) : return true;  break;
    case ( kNuisFSI255                ) : return true;  break;
    case ( kNuisFSI256                ) : return true;  break;
    case ( kNuisFSI257                ) : return true;  break;
    case ( kNuisFSI258                ) : return true;  break;
    case ( kNuisFSI259                ) : return true;  break;
    case ( kNuisFSI260                ) : return true;  break;
    case ( kNuisFSI261                ) : return true;  break;
    case ( kNuisFSI262                ) : return true;  break;
    case ( kNuisFSI263                ) : return true;  break;
    case ( kNuisFSI264                ) : return true;  break;
    case ( kNuisFSI265                ) : return true;  break;
    case ( kNuisFSI266                ) : return true;  break;
    case ( kNuisFSI267                ) : return true;  break;
    case ( kNuisFSI268                ) : return true;  break;
    case ( kNuisFSI269                ) : return true;  break;
    case ( kNuisFSI270                ) : return true;  break;
    case ( kNuisFSI271                ) : return true;  break;
    case ( kNuisFSI272                ) : return true;  break;
    case ( kNuisFSI273                ) : return true;  break;
    case ( kNuisFSI274                ) : return true;  break;
    case ( kNuisFSI275                ) : return true;  break;
    case ( kNuisFSI276                ) : return true;  break;
    case ( kNuisFSI277                ) : return true;  break;
    case ( kNuisFSI278                ) : return true;  break;
    case ( kNuisFSI279                ) : return true;  break;
    case ( kNuisFSI280                ) : return true;  break;
    case ( kNuisFSI281                ) : return true;  break;
    case ( kNuisFSI282                ) : return true;  break;
    case ( kNuisFSI283                ) : return true;  break;
    case ( kNuisFSI284                ) : return true;  break;
    case ( kNuisFSI285                ) : return true;  break;
    case ( kNuisFSI286                ) : return true;  break;
    case ( kNuisFSI287                ) : return true;  break;
    case ( kNuisFSI288                ) : return true;  break;
    case ( kNuisFSI289                ) : return true;  break;
    case ( kNuisFSI290                ) : return true;  break;
    case ( kNuisFSI291                ) : return true;  break;
    case ( kNuisFSI292                ) : return true;  break;
    case ( kNuisFSI293                ) : return true;  break;
    case ( kNuisFSI294                ) : return true;  break;
    case ( kNuisFSI295                ) : return true;  break;
    case ( kNuisFSI296                ) : return true;  break;
    case ( kNuisFSI297                ) : return true;  break;
    case ( kNuisFSI298                ) : return true;  break;
    case ( kNuisFSI299                ) : return true;  break;
    case ( kNuisFSI300                ) : return true;  break;
    case ( kNuisFSI301                ) : return true;  break;
    case ( kNuisFSI302                ) : return true;  break;
    case ( kNuisFSI303                ) : return true;  break;
    case ( kNuisFSI304                ) : return true;  break;
    case ( kNuisFSI305                ) : return true;  break;
    case ( kNuisFSI306                ) : return true;  break;
    case ( kNuisFSI307                ) : return true;  break;
    case ( kNuisFSI308                ) : return true;  break;
    case ( kNuisFSI309                ) : return true;  break;
    case ( kNuisFSI310                ) : return true;  break;
    case ( kNuisFSI311                ) : return true;  break;
    case ( kNuisFSI312                ) : return true;  break;
    case ( kNuisFSI313                ) : return true;  break;
    case ( kNuisFSI314                ) : return true;  break;
    case ( kNuisFSI315                ) : return true;  break;
    case ( kNuisFSI316                ) : return true;  break;
    case ( kNuisFSI317                ) : return true;  break;
    case ( kNuisFSI318                ) : return true;  break;
    case ( kNuisFSI319                ) : return true;  break;
    case ( kNuisFSI320                ) : return true;  break;
    case ( kNuisFSI321                ) : return true;  break;
    case ( kNuisFSI322                ) : return true;  break;
    case ( kNuisFSI323                ) : return true;  break;
    case ( kNuisFSI324                ) : return true;  break;
    case ( kNuisFSI325                ) : return true;  break;
    case ( kNuisFSI326                ) : return true;  break;
    case ( kNuisFSI327                ) : return true;  break;
    case ( kNuisFSI328                ) : return true;  break;
    case ( kNuisFSI329                ) : return true;  break;
    case ( kNuisFSI330                ) : return true;  break;
    case ( kNuisFSI331                ) : return true;  break;
    case ( kNuisFSI332                ) : return true;  break;
    case ( kNuisFSI333                ) : return true;  break;
    case ( kNuisFSI334                ) : return true;  break;
    case ( kNuisFSI335                ) : return true;  break;
    case ( kNuisFSI336                ) : return true;  break;
    case ( kNuisFSI337                ) : return true;  break;
    case ( kNuisFSI338                ) : return true;  break;
    case ( kNuisFSI339                ) : return true;  break;
    case ( kNuisFSI340                ) : return true;  break;
    case ( kNuisFSI341                ) : return true;  break;
    case ( kNuisFSI342                ) : return true;  break;
    case ( kNuisFSI343                ) : return true;  break;
    case ( kNuisFSI344                ) : return true;  break;
    case ( kNuisFSI345                ) : return true;  break;
    case ( kNuisFSI346                ) : return true;  break;
    case ( kNuisFSI347                ) : return true;  break;
    case ( kNuisFSI348                ) : return true;  break;
    case ( kNuisFSI349                ) : return true;  break;
    case ( kNuisFSI350                ) : return true;  break;
    case ( kNuisFSI351                ) : return true;  break;
    case ( kNuisFSI352                ) : return true;  break;
    case ( kNuisFSI353                ) : return true;  break;
    case ( kNuisFSI354                ) : return true;  break;
    case ( kNuisFSI355                ) : return true;  break;
    case ( kNuisFSI356                ) : return true;  break;
    case ( kNuisFSI357                ) : return true;  break;
    case ( kNuisFSI358                ) : return true;  break;
    case ( kNuisFSI359                ) : return true;  break;
    case ( kNuisFSI360                ) : return true;  break;
    case ( kNuisFSI361                ) : return true;  break;
    case ( kNuisFSI362                ) : return true;  break;
    case ( kNuisFSI363                ) : return true;  break;
    case ( kNuisFSI364                ) : return true;  break;
    case ( kNuisFSI365                ) : return true;  break;
    case ( kNuisFSI366                ) : return true;  break;
    case ( kNuisFSI367                ) : return true;  break;
    case ( kNuisFSI368                ) : return true;  break;
    case ( kNuisFSI369                ) : return true;  break;
    case ( kNuisFSI370                ) : return true;  break;
    case ( kNuisFSI371                ) : return true;  break;
    case ( kNuisFSI372                ) : return true;  break;
    case ( kNuisFSI373                ) : return true;  break;
    case ( kNuisFSI374                ) : return true;  break;
    case ( kNuisFSI375                ) : return true;  break;
    case ( kNuisFSI376                ) : return true;  break;
    case ( kNuisFSI377                ) : return true;  break;
    case ( kNuisFSI378                ) : return true;  break;
    case ( kNuisFSI379                ) : return true;  break;
    case ( kNuisFSI380                ) : return true;  break;
    case ( kNuisFSI381                ) : return true;  break;
    case ( kNuisFSI382                ) : return true;  break;
    case ( kNuisFSI383                ) : return true;  break;
    case ( kNuisFSI384                ) : return true;  break;
    case ( kNuisFSI385                ) : return true;  break;
    case ( kNuisFSI386                ) : return true;  break;
    case ( kNuisFSI387                ) : return true;  break;
    case ( kNuisFSI388                ) : return true;  break;
    case ( kNuisFSI389                ) : return true;  break;
    case ( kNuisFSI390                ) : return true;  break;
    case ( kNuisFSI391                ) : return true;  break;
    case ( kNuisFSI392                ) : return true;  break;
    case ( kNuisFSI393                ) : return true;  break;
    case ( kNuisFSI394                ) : return true;  break;
    case ( kNuisFSI395                ) : return true;  break;
    case ( kNuisFSI396                ) : return true;  break;
    case ( kNuisFSI397                ) : return true;  break;
    case ( kNuisFSI398                ) : return true;  break;
    case ( kNuisFSI399                ) : return true;  break;
    case ( kNuisFSI400                ) : return true;  break;
    case ( kNuisFSI401                ) : return true;  break;
    case ( kNuisFSI402                ) : return true;  break;
    case ( kNuisFSI403                ) : return true;  break;
    case ( kNuisFSI404                ) : return true;  break;
    case ( kNuisFSI405                ) : return true;  break;
    case ( kNuisFSI406                ) : return true;  break;
    case ( kNuisFSI407                ) : return true;  break;
    case ( kNuisFSI408                ) : return true;  break;
    case ( kNuisFSI409                ) : return true;  break;
    case ( kNuisFSI410                ) : return true;  break;
    case ( kNuisFSI411                ) : return true;  break;
    case ( kNuisFSI412                ) : return true;  break;
    case ( kNuisFSI413                ) : return true;  break;
    case ( kNuisFSI414                ) : return true;  break;
    case ( kNuisFSI415                ) : return true;  break;
    case ( kNuisFSI416                ) : return true;  break;
    case ( kNuisFSI417                ) : return true;  break;
    case ( kNuisFSI418                ) : return true;  break;
    case ( kNuisFSI419                ) : return true;  break;
    case ( kNuisFSI420                ) : return true;  break;
    case ( kNuisFSI421                ) : return true;  break;
    case ( kNuisFSI422                ) : return true;  break;
    case ( kNuisFSI423                ) : return true;  break;
    case ( kNuisFSI424                ) : return true;  break;
    case ( kNuisFSI425                ) : return true;  break;
    case ( kNuisFSI426                ) : return true;  break;
    case ( kNuisFSI427                ) : return true;  break;
    case ( kNuisFSI428                ) : return true;  break;
    case ( kNuisFSI429                ) : return true;  break;
    case ( kNuisFSI430                ) : return true;  break;
    case ( kNuisFSI431                ) : return true;  break;
    case ( kNuisFSI432                ) : return true;  break;
    case ( kNuisFSI433                ) : return true;  break;
    case ( kNuisFSI434                ) : return true;  break;
    case ( kNuisFSI435                ) : return true;  break;
    case ( kNuisFSI436                ) : return true;  break;
    case ( kNuisFSI437                ) : return true;  break;
    case ( kNuisFSI438                ) : return true;  break;
    case ( kNuisFSI439                ) : return true;  break;
    case ( kNuisFSI440                ) : return true;  break;
    case ( kNuisFSI441                ) : return true;  break;
    case ( kNuisFSI442                ) : return true;  break;
    case ( kNuisFSI443                ) : return true;  break;
    case ( kNuisFSI444                ) : return true;  break;
    case ( kNuisFSI445                ) : return true;  break;
    case ( kNuisFSI446                ) : return true;  break;
    case ( kNuisFSI447                ) : return true;  break;
    case ( kNuisFSI448                ) : return true;  break;
    case ( kNuisFSI449                ) : return true;  break;
    case ( kNuisFSI450                ) : return true;  break;
    case ( kNuisFSI451                ) : return true;  break;
    case ( kNuisFSI452                ) : return true;  break;
    case ( kNuisFSI453                ) : return true;  break;
    case ( kNuisFSI454                ) : return true;  break;
    case ( kNuisFSI455                ) : return true;  break;
    case ( kNuisFSI456                ) : return true;  break;
    case ( kNuisFSI457                ) : return true;  break;
    case ( kNuisFSI458                ) : return true;  break;
    case ( kNuisFSI459                ) : return true;  break;
    case ( kNuisFSI460                ) : return true;  break;
    case ( kNuisFSI461                ) : return true;  break;
    case ( kNuisFSI462                ) : return true;  break;
    case ( kNuisFSI463                ) : return true;  break;
    case ( kNuisFSI464                ) : return true;  break;
    case ( kNuisFSI465                ) : return true;  break;
    case ( kNuisFSI466                ) : return true;  break;
    case ( kNuisFSI467                ) : return true;  break;
    case ( kNuisFSI468                ) : return true;  break;
    case ( kNuisFSI469                ) : return true;  break;
    case ( kNuisFSI470                ) : return true;  break;
    case ( kNuisFSI471                ) : return true;  break;
    case ( kNuisFSI472                ) : return true;  break;
    case ( kNuisFSI473                ) : return true;  break;
    case ( kNuisFSI474                ) : return true;  break;
    case ( kNuisFSI475                ) : return true;  break;
    case ( kNuisFSI476                ) : return true;  break;
    case ( kNuisFSI477                ) : return true;  break;
    case ( kNuisFSI478                ) : return true;  break;
    case ( kNuisFSI479                ) : return true;  break;
#endif
#endif
#endif

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
   const TVectorD & v, int w, int prec, Option_t * option)
{
  TString opt = option;

  bool printing_banff_params = opt.Contains("BANFF fit params");
  bool printing_fsi_params   = opt.Contains("FSI params");

  int np = v.GetNoElements();
  cout << fixed;

  for(int ip = 0; ip < np; ip++) {
     cout << "> param = " << setw(3) << ip;
     if(printing_banff_params) {
         cout << " [ ";
         PhysParamId_t p = PhysParamId_t ( (int)kNuisBANFF00 + ip );
         cout << setw(42) << util::AsString(p);
         cout << " ]";
     }
     if(printing_fsi_params) {
         cout << " [ ";
         PhysParamId_t p = PhysParamId_t ( (int)kNuisFSI000 + ip );
         cout << setw(30) << util::AsString(p);
         cout << " ]";
     }
     cout << " - value = " << setw(w) << setprecision(prec) << v[ip] << endl;
  }
  cout << scientific;
}
//___________________________________________________________________________
void t2kana::util::PrintCompact(
   const TVectorD & current, const TVectorD & nominal, 
   int w, int prec, Option_t * option)
{
  TString opt = option;

  bool printing_banff_params = opt.Contains("BANFF fit params");
  bool printing_fsi_params   = opt.Contains("FSI params");

  int np = current.GetNoElements();
  cout << fixed;

  for(int ip = 0; ip < np; ip++) {
     cout << "> param = " << setw(3) << ip;
     if(printing_banff_params) {
         cout << " [ ";
         PhysParamId_t p = PhysParamId_t ( (int)kNuisBANFF00 + ip );
         cout << setw(42) << util::AsString(p);
         cout << " ]";
     }
     if(printing_fsi_params) {
         cout << " [ ";
         PhysParamId_t p = PhysParamId_t ( (int)kNuisFSI000 + ip );
         cout << setw(30) << util::AsString(p);
         cout << " ]";
     }
     cout << " - current = " 
          << setw(w) << setprecision(prec) << current[ip];
     cout << " - nominal = " 
          << setw(w) << setprecision(prec) << nominal[ip];
     if(nominal[ip] != 0) {
       cout << " - current/nominal = " 
            << setw(w) << setprecision(prec) << current[ip]/nominal[ip];
     }
     cout << endl;
  }
  cout << scientific;
}
//___________________________________________________________________________
void t2kana::util::PrintCompact(
  const TMatrixD & m, int w, int prec)
{
  int nc = m.GetNcols();
  int nr = m.GetNrows();
  cout << fixed;
  for(int ir=0; ir<nr; ir++) {
    if(ir==0) {
       cout << setw(w) << setfill(' ') << "r\\c:" << " | "; 
       for(int ic=0; ic<nc; ic++) {
         cout << setw(w) << setfill(' ') << setprecision(0) << ic << " | ";
       }
       cout << endl;
    }
    cout << setw(w) << setfill(' ') << setprecision(0) << ir << " | "; 
    for(int ic=0; ic<nc; ic++) {
         cout << setw(w) << setfill(' ') << setprecision(prec) 
              << m[ir][ic] << " | ";
    }//rows
    cout << endl;
  }//columns
  cout << scientific;
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
     m == kSKmNuMuNCoth_mu    ||
     m == kSKmNuMuBarNC1piC   ||
     m == kSKmNuMuBarNCoth_mu ||
     m == kSKmNuECCQE         ||
     m == kSKmNuECC1pi        ||
     m == kSKmNuECCcoh        ||
     m == kSKmNuECCoth_mu     ||
     m == kSKmNuENC1piC       ||
     m == kSKmNuENCoth_mu     ||
     m == kSKmNuEBarCCQE      ||
     m == kSKmNuEBarCC1pi     ||
     m == kSKmNuEBarCCcoh     ||
     m == kSKmNuEBarCCoth_mu  ||
     m == kSKmNuEBarNC1piC    ||
     m == kSKmNuEBarNCoth_mu  )
  {
    return 1.;
  }
  else
  if(m == kSKmNuMuCCQE        || 
     m == kSKmNuMuCC1pi       || 
     m == kSKmNuMuCCcoh       || 
     m == kSKmNuMuCCoth_mu    ||
     m == kSKmNuMuBarCCQE     ||
     m == kSKmNuMuBarCC1pi    ||
     m == kSKmNuMuBarCCcoh    ||
     m == kSKmNuMuBarCCoth_mu) 
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
  for(int i = 0; i < knppmax; i++) {
    fXSecSystIdx[i] = -1;
  }
  for(int ip = 0; ip < kNumOfXSecSyst; ip++) {
    PhysParamId_t p = kListOfXSecSyst[ip];
    if(p == kPhypUndefined) continue;
    fXSecSystIdx[p] = ip;
  }
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
}
//____________________________________________________________________________

