//____________________________________________________________________________
/*
 Changes since v1.0.0:

 @ 06/06/11
 Added UseEventToFillInTemplate(SuperkMode_t m, int pdg, int reaction_mode)
 to simplify SuperkMC class.
 
*/
//____________________________________________________________________________

//#define __debug_level_1__
//#define __debug_level_2__
//#define __debug_level_3__

#include <cassert>
#include <cstdlib>
#include <fstream>
#include <iostream>
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

using std::cout;
using std::cerr;
using std::endl;
using std::setw;
using std::setfill;
using std::setprecision;
using std::ofstream;

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
//  if(ndata>0) {
//    chisq = TMath::Power(ndata-nmc,2.)/ndata;
//  }
  if(nmc>0) {
    chisq = TMath::Power(ndata-nmc,2.)/nmc;
  }
  return chisq;
}
//____________________________________________________________________________
double t2kana::util::ChisqPoisson(
  double ndata, double nmc)
{
  //assert(nmc >= 0. && ndata >= 0.);

  if(nmc <= 0. || ndata <= 0.) return 2.*(nmc-ndata);
  return (2. * (nmc-ndata) + 2. * ndata * TMath::Log(ndata/nmc));
}
//____________________________________________________________________________
double t2kana::util::NuMuSurvProbOsc2D(
  double sinsq2theta, double dmsq, double Et)
{
  if(Et <= 0.) return 0.;

  double osc   = TMath::Sin(1.267*dmsq*kBaseline/Et);
  double posc  = sinsq2theta * TMath::Power(osc,2);
  double psurv = 1-posc;

  psurv = TMath::Min(1.,psurv);
  psurv = TMath::Max(0.,psurv);

  return psurv;
}
//____________________________________________________________________________
double t2kana::util::NuMuSurvProbOsc2D(
  double sinsq2theta, double dmsq, double Et, 
  t2kana::conv::SuperkMode_t m)
{
  if(m == kSKmNC         || 
     m == kSKmNuECCQE    ||
     m == kSKmNuECC1pi   ||
     m == kSKmNuECCoth   ||
     m == kSKmNuENC      ||
     m == kSKmNuMuBarNC)
  {
    return 1.;
  }
  else
  if(m == kSKmNuMuCCQE     || 
     m == kSKmNuMuCC1pi    || 
     m == kSKmNuMuCCoth    ||
     m == kSKmNuMuBarCCQE  ||
     m == kSKmNuMuBarCC1pi ||
     m == kSKmNuMuBarCCoth) 
  {
    return NuMuSurvProbOsc2D(sinsq2theta, dmsq, Et);
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
//
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
  if(m == kSKmNuECCQE  ||
     m == kSKmNuECC1pi ||
     m == kSKmNuECCoth ||
     m == kSKmNuENC)
  {
    return 1.;
  }
  else
  if(m == kSKmNuMuCCQE       || 
     m == kSKmNuMuCC1pi      || 
     m == kSKmNuMuCCoth      ||
     m == kSKmNC             ||
     m == kSKmNuMuBarCCQE    ||
     m == kSKmNuMuBarCC1pi   ||
     m == kSKmNuMuBarCCoth   ||
     m == kSKmNuMuBarNC) 
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
  if(m == kSKmNuECCQE  ||
     m == kSKmNuECC1pi ||
     m == kSKmNuECCoth ||
     m == kSKmNuENC)
  {
    return 1.;
  }
  else
  if(m == kSKmNuMuCCQE       || 
     m == kSKmNuMuCC1pi      || 
     m == kSKmNuMuCCoth      ||
     m == kSKmNC             ||
     m == kSKmNuMuBarCCQE    ||
     m == kSKmNuMuBarCC1pi   ||
     m == kSKmNuMuBarCCoth   ||
     m == kSKmNuMuBarNC) 
  {
    return NuMuSurvProbDecoh(sinsq2theta, musq, Et);
  }
  else {
    cerr << "NuMuSurvProbDecoh > Can't handle SK mode = " << m << endl;
    exit(1);
  }
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
     if(L(i,i) < 0) {
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

#ifdef __debug_level_1___
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

  assert(ncols==nrows);
  assert(npars==nrows);

  int n = nrows;

  // create a vector of unit Gaussian variables
  // and multiply by Lt to introduce the appropriate correlations
  TVectorD g(n);
  for (int k = 0; k < n; ++k) {
    g(k) = gR.Gaus();
  }
  g *= cholesky_triangular;

#ifdef __debug_level_3__
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

#ifdef __debug_level_3__
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

#ifdef __debug_level_3__
  cout << "Util >> Tweak vector size = " << g.GetNrows() << endl;
#endif

  return g;
}
//____________________________________________________________________________
string t2kana::util::AsString(PhysParamId_t p)
{
  switch(p) {
    case ( kPhypSinsq2Theta23  ) : return "sinsq2theta23"; break;
    case ( kPhypDm23sq         ) : return "dm23sq";        break;
    case ( kPhypNorm           ) : return "norm";          break;
    case ( kPhypAlphaDecay     ) : return "alpha_decay";   break;
    case ( kPhypMusqDecoherence) : return "musq_decoh";    break;
    default                      : return "--";
  }
  return "--";
}
//____________________________________________________________________________
string t2kana::util::AsString(
    t2kana::conv::SuperkSample_t s)
{
  switch(s) {
     case ( kSKsFCFV1RingMuLike ) : return "FCFV1RingMuLike"; break;
     case ( kSKsFCFVmRingMuLike ) : return "FCFVmRingMuLike"; break;
     default                      : return "--";
  }
  return "--";
}
//____________________________________________________________________________
string t2kana::util::AsString(
     t2kana::conv::SuperkMode_t m)
{
  switch(m) {
     case ( kSKmNuMuCCQE      ) : return "numu CCQE";            break;
     case ( kSKmNuMuCC1pi     ) : return "numu CC1pi";           break;
     case ( kSKmNuMuCCoth     ) : return "numu CCoth";           break;
     case ( kSKmNC            ) : return "numu/nutau NC";        break;
     case ( kSKmNuMuBarCCQE   ) : return "numubar CCQE";         break;
     case ( kSKmNuMuBarCC1pi  ) : return "numubar CC1pi";        break;
     case ( kSKmNuMuBarCCoth  ) : return "numubar CCoth";        break;
     case ( kSKmNuMuBarNC     ) : return "numubar/nutaubar NC";  break;
     case ( kSKmNuECCQE       ) : return "nue CCQE";             break;
     case ( kSKmNuECC1pi      ) : return "nue CC1pi";            break;
     case ( kSKmNuECCoth      ) : return "nue CCoth";            break;
     case ( kSKmNuENC         ) : return "nue NC";               break;
     default                    : return "--";
  }
  return "--";
}
//____________________________________________________________________________
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
int t2kana::util::FluxFix2VarBinNumber(int fixBinNo)
{
  // Given flux cov matrix is 19x19, made for a variable binning scheme:
  // (0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 1.2, 1.5, 2.0, 
  //  2.5, 3.0, 3.5, 4.0, 5.0, 7.0, 10.0) GeV
  // We use fixed-size 0.05 GeV binning.
  // Bin numbers (from the fixed binning scheme to the variable one) are
  // mapped as follows:
  //
  // 0.0(2)0.1; 0.1(2)0.2; .....; 07(2)08; =>16 50Mev bins
  // 0.8(4)1.0; 1.0(4)1.2; = > 8 50MeV bins
  // 1.2(6)1.5; = > 6 50MeV bins
  // 1.5(10)2.0; ....;3.5(10)4.0 => 50 50MeV bins
  // 4.0(20)5.0 =>20 50MeV bins
  // 5(40)7 =>40 50MeV bins
  // 7(40)10 =>60 50MeV bins

/*      
  std::vector<int> binMapping;
            
  //lets zero element to be 0
  binMapping.push_back(0);
  for( int i=1  ;  i<=16; i++){binMapping.push_back((i+1)/2);}
  for( int i=17 ; i<=20; i++){ binMapping.push_back( 9);}
  for( int i=21 ; i<=24; i++){ binMapping.push_back( 10);}
  for( int i=25 ; i<=30; i++){ binMapping.push_back( 11);}
  for( int i=31 ; i<=40; i++){ binMapping.push_back( 12);}
  for( int i=41 ; i<=50; i++){ binMapping.push_back( 13);}
  for( int i=51 ; i<=60; i++){ binMapping.push_back( 14);}
  for( int i=61 ; i<=70; i++){ binMapping.push_back( 15);}
  for( int i=71 ; i<=80; i++){ binMapping.push_back( 16);}
  for( int i=81 ; i<=100; i++){ binMapping.push_back( 17);}
  for( int i=101;i<=140; i++){ binMapping.push_back( 18);}
  for( int i=141;i<=200; i++){ binMapping.push_back( 19);}
      
  return binMapping[fixBinNo];
*/

  if (fixBinNo >=   1 && fixBinNo <=  16) return (fixBinNo+1)/2;
  if (fixBinNo >=  17 && fixBinNo <=  20) return 9;
  if (fixBinNo >=  21 && fixBinNo <=  24) return 10;
  if (fixBinNo >=  25 && fixBinNo <=  30) return 11;
  if (fixBinNo >=  31 && fixBinNo <=  40) return 12;
  if (fixBinNo >=  41 && fixBinNo <=  50) return 13;
  if (fixBinNo >=  51 && fixBinNo <=  60) return 14;
  if (fixBinNo >=  61 && fixBinNo <=  70) return 15;
  if (fixBinNo >=  71 && fixBinNo <=  80) return 16;
  if (fixBinNo >=  81 && fixBinNo <= 100) return 17;
  if (fixBinNo >= 101 && fixBinNo <= 140) return 18;
  if (fixBinNo >= 141 && fixBinNo <= 200) return 19;

  return -1;
}
//___________________________________________________________________________
double t2kana::util::FluxFractionsFixedBinRelToVarBin (
  int fixBinNo, TH1D* fixBinFluxHisto, TH1D* varBinFluxHisto)
{
// ratio 50 MeV bin content relative to the larger bin content (in the 
// var size binning scheme) that contains it.
//
      
  if(fixBinFluxHisto->GetBinCenter(fixBinNo) < 10.0 /* GeV */)
  {
     double fixBinFlux = fixBinFluxHisto->GetBinContent(fixBinNo);
     int ibin = FluxFix2VarBinNumber(fixBinNo);
     double varBinFlux = varBinFluxHisto->GetBinContent(ibin);
     double ratio = (varBinFlux>0) ? fixBinFlux/varBinFlux : 1;
     return ratio;
   }

  cout << "Util/FluxFractionsFixedBinRelToVarBin >> "
       << "Note! Input bin center above 10 GeV" << endl;

   return 1;
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
TH1F * t2kana::util::XML2TH1F(
  TXMLEngine * xmleng, XMLNodePointer_t top_node, string path)
{
// Get a TH1F instance from the info in the input XML config file
//
#ifdef __debug_level_2__
  cout << "utils/XML2TH1F >> Reading info at XML node path: " << path << endl;
#endif

  XMLNodePointer_t found_node = FindNode(xmleng,top_node,path);
  if(!found_node) {
    cerr << "utils/XML2TH1F >> No node was found at input XML node path: "
         << path << endl;
    //return 0;
    exit(1);
  }

  string filename = "";
  string objname  = "";

  XMLNodePointer_t child_node = xmleng->GetChild(found_node);
  while(child_node) {
    string child_node_name = xmleng->GetNodeName(child_node);
    string content = xmleng->GetNodeContent(child_node);
#ifdef __debug_level_2__
    cout << "utils/XML2TH1F >> node = " << child_node_name << endl;
    cout << "utils/XML2TH1F >> \t content = " << content << endl;
#endif
    if(child_node_name == "filename") { filename = TrimSpaces(content); }
    if(child_node_name == "objname" ) { objname  = TrimSpaces(content); }
    child_node = xmleng->GetNext(child_node);    
  }

#ifdef __debug_level_2__
  cout << "utils/XML2TH1F >> filename : " << filename << endl;
  cout << "utils/XML2TH1F >> objname  : " << objname  << endl;
#endif

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
  TXMLEngine * xmleng, XMLNodePointer_t top_node, string path)
{
// Get a TH1D instance from the info in the input XML config file
//
#ifdef __debug_level_2__
  cout << "utils/XML2TH1D >> Reading info at XML node path: " << path << endl;
#endif

  XMLNodePointer_t found_node = FindNode(xmleng,top_node,path);
  if(!found_node) {
    cerr << "utils/XML2TH1D >> No node was found at input XML node path: "
         << path << endl;
    //return 0;
    exit(1);
  }

  string filename = "";
  string objname  = "";

  XMLNodePointer_t child_node = xmleng->GetChild(found_node);
  while(child_node) {
    string child_node_name = xmleng->GetNodeName(child_node);
    string content = xmleng->GetNodeContent(child_node);
#ifdef __debug_level_2__
    cout << "utils/XML2TH1D >> node = " << child_node_name << endl;
    cout << "utils/XML2TH1D >> \t content = " << content << endl;
#endif
    if(child_node_name == "filename") { filename = TrimSpaces(content); }
    if(child_node_name == "objname" ) { objname  = TrimSpaces(content); }
    child_node = xmleng->GetNext(child_node);    
  }

#ifdef __debug_level_2__
  cout << "utils/XML2TH1D >> filename : " << filename << endl;
  cout << "utils/XML2TH1D >> objname  : " << objname  << endl;
#endif

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
  TXMLEngine * xmleng, XMLNodePointer_t top_node, string path)
{
// Get a TH2D instance from the info in the input XML config file
//
#ifdef __debug_level_2__
  cout << "utils/XML2TH2D >> Reading info at XML node path: " << path << endl;
#endif

  XMLNodePointer_t found_node = FindNode(xmleng,top_node,path);
  if(!found_node) {
    cerr << "utils/XML2TH2D >> No node was found at input XML node path: "
         << path << endl;
    //return 0;
    exit(1);
  }

  string filename = "";
  string objname  = "";

  XMLNodePointer_t child_node = xmleng->GetChild(found_node);
  while(child_node) {
    string child_node_name = xmleng->GetNodeName(child_node);
    string content = xmleng->GetNodeContent(child_node);
#ifdef __debug_level_2__
    cout << "utils/XML2TH2D >> node = " << child_node_name << endl;
    cout << "utils/XML2TH2D >> \t content = " << content << endl;
#endif
    if(child_node_name == "filename") { filename = TrimSpaces(content); }
    if(child_node_name == "objname" ) { objname  = TrimSpaces(content); }
    child_node = xmleng->GetNext(child_node);    
  }

#ifdef __debug_level_2__
  cout << "utils/XML2TH1D >> filename : " << filename << endl;
  cout << "utils/XML2TH1D >> objname  : " << objname  << endl;
#endif

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
   TXMLEngine * xmleng, XMLNodePointer_t top_node, string path)
{
// Get a TH2D instance from the info in the input XML config file
//
#ifdef __debug_level_2__
  cout << "utils/XML2TTree >> Reading info at XML node path: " << path << endl;
#endif

  XMLNodePointer_t found_node = FindNode(xmleng,top_node,path);
  if(!found_node) {
    cerr << "utils/XML2TTree >> No node was found at input XML node path: "
         << path << endl;
    //return 0;
    exit(1);
  }

  string filename = "";
  string objname  = "";

  XMLNodePointer_t child_node = xmleng->GetChild(found_node);
  while(child_node) {
    string child_node_name = xmleng->GetNodeName(child_node);
    string content = xmleng->GetNodeContent(child_node);
#ifdef __debug_level_2__
    cout << "utils/XML2TTree >> node = " << child_node_name << endl;
    cout << "utils/XML2TTree >> \t content = " << content << endl;
#endif
    if(child_node_name == "filename") { filename = TrimSpaces(content); }
    if(child_node_name == "objname" ) { objname  = TrimSpaces(content); }
    child_node = xmleng->GetNext(child_node);    
  }

#ifdef __debug_level_2__
  cout << "utils/XML2TTree >> filename : " << filename << endl;
  cout << "utils/XML2TTree >> objname  : " << objname  << endl;
#endif

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
TAxis *  t2kana::util::XML2BinInfo(
  TXMLEngine * xmleng, XMLNodePointer_t top_node, string path)
{
//
//

#ifdef __debug_level_2__
  cout << "utils/XML2BinInfo >> Reading info at XML node path: " << path << endl;
#endif

  XMLNodePointer_t found_node = FindNode(xmleng,top_node,path);
  if(!found_node) {
    cerr << "utils/XML2TH2D >> No node was found at input XML node path: "
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
#ifdef __debug_level_2__
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
double t2kana::util::XML2Dbl(
  TXMLEngine * xmleng, XMLNodePointer_t top_node, string path)
{
// Get the specified node content as a floating-point number
//

#ifdef __debug_level_2__
  cout << "utils/XML2Dbl >> Reading info at XML node path: " << path << endl;
#endif

  XMLNodePointer_t found_node = FindNode(xmleng,top_node,path);
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
  TXMLEngine * xmleng, XMLNodePointer_t top_node, string path)
{
// Get the specified node content as a vector floating-point number
// The node content should be a comma-separed list of numbers

#ifdef __debug_level_2__
  cout << "utils/XML2Dbl >> Reading info at XML node path: " << path << endl;
#endif

  XMLNodePointer_t found_node = FindNode(xmleng,top_node,path);
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
  TXMLEngine * xmleng, XMLNodePointer_t top_node, string path)
{
// Get the specified node content as a floating-point number
//

#ifdef __debug_level_2__
  cout << "utils/XML2Str >> Reading info at XML node path: " << path << endl;
#endif

  XMLNodePointer_t found_node = FindNode(xmleng,top_node,path);
  if(!found_node) {
    cerr << "utils/XML2Str >> No node was found at input XML node path: "
         << path << endl;
    exit(1);
  }
  string content = xmleng->GetNodeContent(found_node);
  return TrimSpaces(content);
}
//___________________________________________________________________________
TMatrixDSym * t2kana::util::XML2TMatrixDSym(
   TXMLEngine * xmleng, XMLNodePointer_t top_node, string path)
{
// Get a TMatrixDSym instance from the info in the input XML config file
//

#ifdef __debug_level_2__
  cout << "utils/XML2TMatrixDSym >> "
       << "Reading info at XML node path: " << path << endl;
#endif

  XMLNodePointer_t found_node = FindNode(xmleng,top_node,path);
  if(!found_node) {
    cerr << "utils/XML2TMatrixDSym >> "
         << "No node was found at input XML node path: " << path << endl;
    //return 0;
    exit(1);
  }

  string filename = "";
  string objname  = "";

  XMLNodePointer_t child_node = xmleng->GetChild(found_node);
  while(child_node) {
    string child_node_name = xmleng->GetNodeName(child_node);
    string content = xmleng->GetNodeContent(child_node);
    if(child_node_name == "filename") { filename = TrimSpaces(content); }
    if(child_node_name == "objname" ) { objname  = TrimSpaces(content); }
    child_node = xmleng->GetNext(child_node);    
  }

#ifdef __debug_level_2__
  cout << "utils/XML2TMatrixDSym >> filename : " << filename << endl;
  cout << "utils/XML2TMatrixDSym >> objname  : " << objname  << endl;
#endif

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
    TXMLEngine * xmleng, XMLNodePointer_t top_node, string path)
{
// Get a TMatrixD instance from the info in the input XML config file
//

#ifdef __debug_level_2__
  cout << "utils/XML2TMatrixD >> "
       << "Reading info at XML node path: " << path << endl;
#endif

  XMLNodePointer_t found_node = FindNode(xmleng,top_node,path);
  if(!found_node) {
    cerr << "utils/XML2TMatrixD >> "
         << "No node was found at input XML node path: " << path << endl;
    //return 0;
    exit(1);
  }

  string filename = "";
  string objname  = "";

  XMLNodePointer_t child_node = xmleng->GetChild(found_node);
  while(child_node) {
    string child_node_name = xmleng->GetNodeName(child_node);
    string content = xmleng->GetNodeContent(child_node);
    if(child_node_name == "filename") { filename = TrimSpaces(content); }
    if(child_node_name == "objname" ) { objname  = TrimSpaces(content); }
    child_node = xmleng->GetNext(child_node);    
  }

#ifdef __debug_level_2__
  cout << "utils/XML2TMatrixD >> filename : " << filename << endl;
  cout << "utils/XML2TMatrixD >> objname  : " << objname  << endl;
#endif

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
bool t2kana::util::UseEventToFillInTemplate(
  SuperkMode_t m, int nu_pdg_code, int neut_reaction_mode)
{
  bool use_it = false;

  if(m == kSKmNuMuCCQE) 
  {
     bool is_numu = (nu_pdg_code == kPdgNuMu);
     bool is_ccqe = (neut_reaction_mode == 1);
     use_it = is_numu && is_ccqe;
  } 

  else
  if(m == kSKmNuMuCC1pi) 
  {
     bool is_numu  = (nu_pdg_code == kPdgNuMu);
     bool is_cc1pi = (neut_reaction_mode == 11 || 
                      neut_reaction_mode == 12 || 
                      neut_reaction_mode == 13);
     use_it = is_numu && is_cc1pi;
  } 

  else
  if(m == kSKmNuMuCCoth) 
  {
     bool is_numu  = (nu_pdg_code == kPdgNuMu);
     bool is_ccoth = (neut_reaction_mode == 16 || 
                      neut_reaction_mode == 17 || 
                      neut_reaction_mode == 21 ||
                      neut_reaction_mode == 22 || 
                      neut_reaction_mode == 23 || 
                      neut_reaction_mode == 26);
     use_it = is_numu && is_ccoth;
  } 

  else
  if(m == kSKmNC) 
  {
     bool is_numu  = (nu_pdg_code == kPdgNuMu);
     bool is_nc    = (neut_reaction_mode == 31 ||
                      neut_reaction_mode == 32 ||
                      neut_reaction_mode == 36 ||
                      neut_reaction_mode == 33 ||
                      neut_reaction_mode == 34 ||
                      neut_reaction_mode == 38 ||
                      neut_reaction_mode == 39 ||
                      neut_reaction_mode == 41 ||
                      neut_reaction_mode == 42 ||
                      neut_reaction_mode == 43 ||
                      neut_reaction_mode == 46 ||
                      neut_reaction_mode == 51 ||
                      neut_reaction_mode == 52);
     use_it = is_numu && is_nc;
  } 

  else
  if(m == kSKmNuMuBarCCQE) 
  {
     bool is_numubar = (nu_pdg_code == kPdgNuMuBar);
     bool is_ccqe    = (neut_reaction_mode == -1);
     use_it = is_numubar && is_ccqe;
  } 

  else
  if(m == kSKmNuMuBarCC1pi) 
  {
     bool is_numubar = (nu_pdg_code == kPdgNuMuBar);
     bool is_cc1pi   = (neut_reaction_mode == -11 || 
                        neut_reaction_mode == -12 || 
                        neut_reaction_mode == -13);
     use_it = is_numubar && is_cc1pi;
  } 

  else
  if(m == kSKmNuMuBarCCoth) 
  {
     bool is_numubar = (nu_pdg_code == kPdgNuMuBar);
     bool is_ccqe    = (neut_reaction_mode == -1);
     bool is_cc1pi   = (neut_reaction_mode == -11 || 
                        neut_reaction_mode == -12 || 
                        neut_reaction_mode == -13);
     bool is_nc      = (neut_reaction_mode == -31 ||
                        neut_reaction_mode == -32 ||
                        neut_reaction_mode == -36 ||
                        neut_reaction_mode == -33 ||
                        neut_reaction_mode == -34 ||
                        neut_reaction_mode == -38 ||
                        neut_reaction_mode == -39 ||
                        neut_reaction_mode == -41 ||
                        neut_reaction_mode == -42 ||
                        neut_reaction_mode == -43 ||
                        neut_reaction_mode == -44 ||
                        neut_reaction_mode == -45 ||
                        neut_reaction_mode == -46 ||
                        neut_reaction_mode == -51 ||
                        neut_reaction_mode == -52);
     use_it = is_numubar && !(is_nc||is_ccqe||is_cc1pi);
  } 

  else
  if(m == kSKmNuMuBarNC) 
  {
     bool is_numubar  = (nu_pdg_code == kPdgNuMuBar);
     bool is_nc       = (neut_reaction_mode == -31 ||
                         neut_reaction_mode == -32 ||
                         neut_reaction_mode == -36 ||
                         neut_reaction_mode == -33 ||
                         neut_reaction_mode == -34 ||
                         neut_reaction_mode == -38 ||
                         neut_reaction_mode == -39 ||
                         neut_reaction_mode == -41 ||
                         neut_reaction_mode == -42 ||
                         neut_reaction_mode == -43 ||
                         neut_reaction_mode == -44 ||
                         neut_reaction_mode == -45 ||
                         neut_reaction_mode == -46 ||
                         neut_reaction_mode == -51 ||
                         neut_reaction_mode == -52);
     use_it = is_numubar && is_nc;
  } 

  else
  if(m == kSKmNuECCQE) 
  {
     bool is_nue  = (nu_pdg_code == kPdgNuE || 
                     nu_pdg_code == kPdgNuEBar);
     bool is_ccqe = (neut_reaction_mode == 1 || 
                     neut_reaction_mode == -1);
     use_it = is_nue && is_ccqe;
  } 

  else
  if(m == kSKmNuECC1pi) 
  {
     bool is_nue   = (nu_pdg_code == kPdgNuE || 
                      nu_pdg_code == kPdgNuEBar);
     bool is_cc1pi = (neut_reaction_mode ==  11 || 
                      neut_reaction_mode ==  12 || 
                      neut_reaction_mode ==  13 ||
                      neut_reaction_mode == -11 || 
                      neut_reaction_mode == -12 || 
                      neut_reaction_mode == -13);
     use_it = is_nue && is_cc1pi;
  }

  else
  if(m == kSKmNuECCoth) 
  {
     bool is_nue    = (nu_pdg_code == kPdgNuE || 
                       nu_pdg_code == kPdgNuEBar);
     bool is_ccqe   = (neut_reaction_mode == 1   || 
                       neut_reaction_mode == -1);
     bool is_cc1pi  = (neut_reaction_mode ==  11 || 
                       neut_reaction_mode ==  12 || 
                       neut_reaction_mode ==  13 ||
                       neut_reaction_mode == -11 || 
                       neut_reaction_mode == -12 || 
                       neut_reaction_mode == -13);
     bool is_nc     = (neut_reaction_mode ==  31 ||
                       neut_reaction_mode ==  32 ||
                       neut_reaction_mode ==  36 ||
                       neut_reaction_mode ==  33 ||
                       neut_reaction_mode ==  34 ||
                       neut_reaction_mode ==  38 ||
                       neut_reaction_mode ==  39 ||
                       neut_reaction_mode ==  41 ||
                       neut_reaction_mode ==  42 ||
                       neut_reaction_mode ==  43 ||
                       neut_reaction_mode ==  46 ||
                       neut_reaction_mode ==  51 ||
                       neut_reaction_mode ==  52 ||
                       neut_reaction_mode == -31 ||
                       neut_reaction_mode == -32 ||
                       neut_reaction_mode == -36 ||
                       neut_reaction_mode == -33 ||
                       neut_reaction_mode == -34 ||
                       neut_reaction_mode == -38 ||
                       neut_reaction_mode == -39 ||
                       neut_reaction_mode == -41 ||
                       neut_reaction_mode == -42 ||
                       neut_reaction_mode == -43 ||
                       neut_reaction_mode == -44 ||
                       neut_reaction_mode == -45 ||
                       neut_reaction_mode == -46 ||
                       neut_reaction_mode == -51 ||
                       neut_reaction_mode == -52);
     use_it = is_nue && !(is_ccqe||is_cc1pi||is_nc);
  }

  else
  if(m == kSKmNuENC) 
  {
     bool is_nue  = (nu_pdg_code == kPdgNuE || 
                     nu_pdg_code == kPdgNuEBar);
     bool is_nc   = (neut_reaction_mode ==  31 ||
                     neut_reaction_mode ==  32 ||
                     neut_reaction_mode ==  36 ||
                     neut_reaction_mode ==  33 ||
                     neut_reaction_mode ==  34 ||
                     neut_reaction_mode ==  38 ||
                     neut_reaction_mode ==  39 ||
                     neut_reaction_mode ==  41 ||
                     neut_reaction_mode ==  42 ||
                     neut_reaction_mode ==  43 ||
                     neut_reaction_mode ==  46 ||
                     neut_reaction_mode ==  51 ||
                     neut_reaction_mode ==  52 ||
                     neut_reaction_mode == -31 ||
                     neut_reaction_mode == -32 ||
                     neut_reaction_mode == -36 ||
                     neut_reaction_mode == -33 ||
                     neut_reaction_mode == -34 ||
                     neut_reaction_mode == -38 ||
                     neut_reaction_mode == -39 ||
                     neut_reaction_mode == -41 ||
                     neut_reaction_mode == -42 ||
                     neut_reaction_mode == -43 ||
                     neut_reaction_mode == -44 ||
                     neut_reaction_mode == -45 ||
                     neut_reaction_mode == -46 ||
                     neut_reaction_mode == -51 ||
                     neut_reaction_mode == -52);
     use_it = is_nue && is_nc;
  }
 
  return use_it;
}
//___________________________________________________________________________

