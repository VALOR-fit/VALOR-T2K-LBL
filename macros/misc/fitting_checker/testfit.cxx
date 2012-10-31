//#define __print_generated_toy_expt__
//#define __print_fit_func_debug_mesgs__
//#define __print_chisq_debug__mesgs__
//#define SYST_FLUCT_MODE 0 //no fluctuations
#define SYST_FLUCT_MODE 1 //uncorrelated fluctuations
//#define SYST_FLUCT_MODE 2 //correlated fluctuations
//#define DONT_FIT_ALL 0 // fit all systematics
#define DONT_FIT_ALL 1 // fit all but 1 systematic
//#define DONT_FIT_ALL 2 // fit all but 2 systematics
#define FIX_NOT_FIT // fix the systematics not fit to their bestfit. otherwise, tweak according to correlations
#define CORRELATED_CHISQ // use the correlation matrix to define the chisq. otherwise, use uncorrelated chisq
//#define FLUCT_NOT_FIT // don't fluctuate the parameters that aren't fit

#include <iostream>
#include <cassert>

#include "TFile.h"
#include "TH1.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TRandom3.h"
#include "TVirtualFitter.h"
#include "TMath.h"
#include "TMinuit.h"
#include "TNtuple.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TStyle.h"

const int    kntoymc = 10000;
const int    knbins  = 200;
const double kxmin   = 0.;
const double kxmax   = 10.;
const double ka1err  = 0.02;
const double ka2err  = 0.06;
const double ka3err  = 0.10;
const double ka4err  = 0.01;
const double ka12cor = -0.015;
const double ka13cor = -0.002;
const double ka14cor = -0.006;
const double ka23cor = 0.008;
const double ka24cor = 0.001;
const double ka34cor = 0.01;
const double klevel  = 20;
const double kslope  = 6;
const double kbump   = 4;
const double true_a1   = 1.0;
const double true_a2   = 1.0;
const double true_a3   = 1.0;
const double true_a4   = 1.0;
const double true_mean = 4.0;
const double true_dip  = 0.8;
TMatrixD cov(4,4);
TMatrixD covInv(4,4);

using namespace std;

double pdf          (double x, double mean, double dip, double a1, double a2, double a3, double a4);
TH1D * pdf          (double mean, double dip, double a1, double a2, double a3, double a4);
TH1D * mock_data    (double mean, double dip, double a1, double a2, double a3, double a4);
int    fit          (TH1D * data, double & mean, double & dip, double & a1, double & a2, double & a3, double & a4, double & mean_err, double & dip_err, double & a1_err, double & a2_err, double & a3_err, double & a4_err);
void   fitfunc      (Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t iflag);
double chisq        (TH1D * data, double mean, double dip, double a1, double a2, double a3, double a4, bool penalty_only);
void   init_fitter  (void);
void   set_fit_data (TH1D * data);
TVectorD CholeskyGenerateCorrelatedParams (TMatrixD& cholesky_triangular, TVectorD& mean_params);
TMatrixD CholeskyDecomposition(TMatrixD& cov_matrix);
int    main         (void);

//.............................................................................
double pdf(double x, double mean, double dip, double a1, double a2, double a3, double a4) 
{
// function is a constant level with a gaussian dip on top of linear background
// and a gaussian bump towards the bottom of the range
//
   double xreco = a1*x;
   double xrmmsq = TMath::Power(xreco-mean,2);
   double xrsq   = TMath::Power(xreco,2);
   double fdip   = a2 * klevel * (1 - dip * TMath::Exp(-1*xrmmsq));
   double flin   = x * (a3*kslope / kxmax);
   double fbump  = a4 * kbump * TMath::Exp(-1*xrsq);
   double f      = fdip+flin+fbump;
   return f;
}
//.............................................................................
TH1D * pdf(double mean, double dip, double a1, double a2, double a3, double a4)
{
  TH1D * h = new TH1D("hpdf","pdf",knbins,kxmin,kxmax);
  h->SetDirectory(0);

  for(int ibin = 1; ibin <= knbins;  ibin++) {
	double x = h->GetBinCenter(ibin);
        double f = pdf(x,mean,dip,a1,a2,a3,a4);
        h->SetBinContent(ibin,f);
  }
  return h;
}
//.............................................................................
void fitfunc(Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t iflag)
{
  TVirtualFitter * fitter = TVirtualFitter::GetFitter();

  TH1D * data = (TH1D *) fitter->GetObjectFit();
  if(!data) {
	cout << " *** no data! ***" << endl;
	exit(1);
  }

  // get fit params 
  double mean = par[0];
  double dip  = par[1];
  double a4   = par[2];
  double a3   = par[3];

#if DONT_FIT_ALL == 2
#ifdef FIX_NOT_FIT
  double a2 = true_a2;
#else
  double a2 = true_a2 + (((a3 - true_a3) / cov[2][2]) * cov[1][2]) + (((a4 - true_a4) / cov[3][3]) * cov[1][3]);
#endif
#elif DONT_FIX_ALL <= 1
  double a2 = par[4];
#endif

#if DONT_FIT_ALL >= 1
#ifdef FIX_NOT_FIT
  double a1 = true_a1;
#else
#if DONT_FIT_ALL == 1
  double a1 = true_a1 + (((a2 - true_a2) / cov[1][1]) * cov[0][1]) + (((a3 - true_a3) / cov[2][2]) * cov[0][2]) + (((a4 - true_a4) / cov[3][3]) * cov[0][3]);
#elif DONT_FIT_ALL == 2
  double a1 = true_a1 + (((a3 - true_a3) / cov[2][2]) * cov[0][2]) + (((a4 - true_a4) / cov[3][3]) * cov[0][3]);
#endif
#endif
#elif DONT_FIX_ALL == 0
  double a1 = par[5];
#endif

#ifdef __print_fit_func_debug_mesgs__
  cout << "Minimization function / Trying parameters: mean = " 
       << mean << ", dip = " << dip 
       << "a1 = " << a1 << ", a2 = " << a2 << ", a3 = " << a3 ", a4 = " << a4 << endl;
#endif

  f = chisq(data, mean, dip, a1, a2, a3, a4, false);
}
//.............................................................................
double chisq(TH1D * data, double mean, double dip, double a1, double a2, double a3, double a4, bool penalty_only)
{
  double chisq = 0;

  if(!penalty_only) {
    int nbins = data->GetNbinsX();
    for(int ibin = 1; ibin <= nbins; ibin++) {       
       double x    = data->GetBinCenter(ibin);
       double nobs = data->GetBinContent(ibin);
       double nexp = pdf(x,mean,dip,a1,a2,a3,a4);
       double delta_chisq = 0;  
       if(nexp <= 0. || nobs <= 0.) {
          delta_chisq = 2.*(nexp-nobs);
       } else {
          delta_chisq = (2. * (nexp-nobs) + 2. * nobs * TMath::Log(nobs/nexp));
       }
       chisq += delta_chisq;     
    }//ibin
  }

  // penalty term
#ifdef CORRELATED_CHISQ
  //Vector with deviations of current values from  best fit values
  const int np = 4;
  TMatrixD dev  (1,np); //vector of one column in matrix form for multiplication
  TMatrixD devT (np,1); //vector of one row in matrix form for multiplication
  dev [0][0] = a1-true_a1;
  devT[0][0] = a1-true_a1;
  dev [0][1] = a2-true_a2;
  devT[1][0] = a2-true_a2;
  dev [0][2] = a3-true_a3;
  devT[2][0] = a3-true_a3;
  dev [0][3] = a4-true_a4;
  devT[3][0] = a4-true_a4;
  TMatrixD chisqpenalty =  dev * covInv * devT;
  double chisq_penalty = chisqpenalty[0][0]; //The final matrix is 1x1, just a number
  /*
  dev.Print();
  covInv.Print();
  devT.Print();
  chisqpenalty.Print();
  */
#else
  double chisq_penalty = 0;
#if DONT_FIT_ALL == 0
  chisq_penalty += TMath::Power(true_a1-a1,2.) / TMath::Power(ka1err,2.);
#endif
#if DONT_FIT_ALL <= 1
  chisq_penalty += TMath::Power(true_a2-a2,2.) / TMath::Power(ka2err,2.);
#endif
  chisq_penalty += TMath::Power(true_a3-a3,2.) / TMath::Power(ka3err,2.);
  chisq_penalty += TMath::Power(true_a4-a4,2.) / TMath::Power(ka4err,2.);
#endif

  chisq += chisq_penalty;

#ifdef __print_chisq_debug__mesgs__
  cout << "Minimization function / chisq = " << chisq << " (penalty term = " << chisq_penalty << ")" << endl;
#endif

  return chisq;
}
//.............................................................................
TH1D * mock_data(double mean, double dip, double a1, double a2, double a3, double a4)
{
  TRandom3 rand(0);

  cout << "True parameters: " 
       << "mean = " << mean << "dip = " << dip 
       << ", a1 = " << a1 << ", a2 = " << a2 << ", a3 = " << a3 << ", a4 = " << a4 << endl;

  TH1D * h = new TH1D("h","toy MC expt",knbins,kxmin,kxmax);
  h->SetDirectory(0);

  for(int ibin = 1; ibin <= knbins;  ibin++) {
	double x = h->GetBinCenter(ibin);
        double f = pdf(x,mean,dip,a1,a2,a3,a4);
        double n = rand.Poisson(f);
        h->SetBinContent(ibin,n);
#ifdef __print_generated_toy_expt__
	cout << "     x = " << x << " (bin = " << ibin << ") : pdf = " << f 
             << ", this toy-MC expt = " << n << endl;
#endif
  }
  return h;
}
//.............................................................................
int fit(TH1D * data, 
  double & mean, double & dip, double & a1, double & a2, double & a3, double & a4, double & mean_err, double & dip_err, double & a1_err, double & a2_err, double & a3_err, double & a4_err)
{
  //cout << "Getting fitter" << endl;
  TVirtualFitter * fitter = TVirtualFitter::GetFitter();
  assert(fitter);

  // MINUIT's minimization step
  //cout << "Running MIGRAD..." << endl;
  Double_t arglist[5];
  arglist[0] = 5000; // num of func calls
  arglist[1] = 0.01; // tolerance
  fitter->ExecuteCommand("MIGRAD",arglist, 2);
 
  //Use Hesse to make the errors better
  fitter->ExecuteCommand("HESSE",arglist, 1);

  //cout << "Fit results: " << endl;
  Double_t amin = 0;
  fitter->PrintResults(3,amin);
   
  // Get error matrix status
  double ha(0.0);     // dummy vars
  double edm, errdef; // dummy vars
  int nvpar, nparx;   // dummy vars
  int status_code = fitter->GetStats(ha, edm, errdef, nvpar, nparx);
  
  mean = fitter->GetParameter(0);
  dip  = fitter->GetParameter(1);
  a4   = fitter->GetParameter(2);
  a3   = fitter->GetParameter(3);
#if DONT_FIT_ALL <= 1
  a2   = fitter->GetParameter(4);
#else
  a2   = true_a2;
#endif
#if DONT_FIT_ALL == 0
  a1   = fitter->GetParameter(5);
#else
  a1   = true_a1;
#endif

  mean_err = fitter->GetParError(0);
  dip_err  = fitter->GetParError(1);
  a4_err   = fitter->GetParError(2);
  a3_err   = fitter->GetParError(3);
#if DONT_FIT_ALL <= 1
  a2_err   = fitter->GetParError(4);
#endif
#if DONT_FIT_ALL == 0
  a1_err   = fitter->GetParError(5);
#endif

  cout << "Done fitting! Error matrix status = " << status_code << endl;

  return status_code;
}
//.............................................................................
void init_fitter(void)
{
  // get fitter
  //cout << "Setting fitter (MINUIT)" << endl;
  TVirtualFitter::SetDefaultFitter("Minuit");
  TVirtualFitter * fitter = TVirtualFitter::Fitter(0,6);
  if(!fitter) return;
  
  // set print level
  //cout << "Setting print-level" << endl;
  Double_t arglist[5];
  arglist[0] = 0;
  fitter->ExecuteCommand("SET PRINT",arglist,1);
    
  // set params
  //cout << "Setting fit params" << endl;
  fitter->Clear();
  fitter->SetParameter(0, "mean", 5.00, 0.01, 0.00, 10.00 );
  fitter->SetParameter(1, "dip",  0.80, 0.01, 0.30,  1.00 );
  fitter->SetParameter(2, "a4",   true_a1, 0.01, 0.00,  2.00 );
  fitter->SetParameter(3, "a3",   true_a2, 0.01, 0.00,  2.00 );
#if DONT_FIT_ALL <= 1
  fitter->SetParameter(4, "a2",   true_a3, 0.01, 0.00,  2.00 );
#endif
#if DONT_FIT_ALL == 0
  fitter->SetParameter(5, "a1",   true_a4, 0.01, 0.00,  2.00 );
#endif

  // set minimization function
  //cout << "Setting functions" << endl;
  fitter->SetFCN(fitfunc);
}
//.............................................................................
void set_fit_data (TH1D * data)
{
  TVirtualFitter * fitter = TVirtualFitter::GetFitter();
  assert(fitter);
  fitter->SetObjectFit(data);
}
//.............................................................................
TMatrixD CholeskyDecomposition(TMatrixD& cov_matrix)
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
//.............................................................................
TVectorD CholeskyGenerateCorrelatedParams (
    TMatrixD& cholesky_triangular, TVectorD& mean_params) 
{
// Generate a vector of correlated params

  int ncols = cholesky_triangular.GetNcols();
  int nrows = cholesky_triangular.GetNrows();
  int npars = mean_params.GetNrows();

  TRandom3 rand(0);

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
    g(k) = rand.Gaus();
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
//.............................................................................
int main(void)
{

#ifndef SYST_FLUCT_MODE
  return 1;
#endif
#ifndef DONT_FIT_ALL
  return 1;
#endif

  cout << "Start test..." << endl;

  cov.Zero();
  cov[0][0] = pow(ka1err,2);
  cov[1][1] = pow(ka2err,2);
  cov[2][2] = pow(ka3err,2);
  cov[3][3] = pow(ka4err,2);
  cov[0][1] = ka1err * ka2err * ka12cor;
  cov[1][0] = ka1err * ka2err * ka12cor;
  cov[0][2] = ka1err * ka3err * ka13cor;
  cov[2][0] = ka1err * ka3err * ka13cor;
  cov[0][3] = ka1err * ka4err * ka14cor;
  cov[3][0] = ka1err * ka4err * ka14cor;
  cov[1][2] = ka1err * ka2err * ka23cor;
  cov[2][1] = ka1err * ka2err * ka23cor;
  cov[1][3] = ka1err * ka2err * ka24cor;
  cov[3][1] = ka1err * ka2err * ka24cor;
  cov[2][3] = ka1err * ka2err * ka23cor;
  cov[3][2] = ka1err * ka2err * ka23cor;
  cov.Print();

  double det;
  covInv = cov;
  covInv.Invert(&det);
  covInv.Print();

  TMatrixD unit = covInv * cov;
  unit.Print();

  TMatrixD covL = CholeskyDecomposition(cov);
  assert(covL.GetNrows() == 4);
  covL.Print();

#if SYST_FLUCT_MODE == 2
  TVectorD bf(4);
  bf[0] = true_a1;
  bf[1] = true_a2;
  bf[2] = true_a3;
  bf[3] = true_a4;
#endif

  TFile f("./out.root","recreate");
  TNtuple * bestfit_params = new TNtuple("bestfit_params","","status:m:d:a1:a2:a3:a4:chisq:chisqp:merr:derr:a1err:a2err:a3err:a4err");
  TH1D * pull_mean = new TH1D("pull_mean","pull_mean",200,-5,5);
  TH1D * pull_dip  = new TH1D("pull_dip", "pull_dip", 200,-5,5);
  TH1D * pull_a1   = new TH1D("pull_a1",  "pull_a1",  200,-5,5);
  TH1D * pull_a2   = new TH1D("pull_a2",  "pull_a2",  200,-5,5);
  TH1D * pull_a3   = new TH1D("pull_a3",  "pull_a3",  200,-5,5);
  TH1D * pull_a4   = new TH1D("pull_a4",  "pull_a4",  200,-5,5);
  TH1D * pullalt_a1   = new TH1D("pullalt_a1",  "pullalt_a1",  200,-5,5);
  TH1D * pullalt_a2   = new TH1D("pullalt_a2",  "pullalt_a2",  200,-5,5);
  TH1D * pullalt_a3   = new TH1D("pullalt_a3",  "pullalt_a3",  200,-5,5);
  TH1D * pullalt_a4   = new TH1D("pullalt_a4",  "pullalt_a4",  200,-5,5);
  TH1D * throw_a1 = new TH1D("throw_a1", "throw_a1", 1000,-5,5);
  TH1D * throw_a2 = new TH1D("throw_a2", "throw_a2", 1000,-5,5);
  TH1D * throw_a3 = new TH1D("throw_a3", "throw_a3", 1000,-5,5);
  TH1D * throw_a4 = new TH1D("throw_a4", "throw_a4", 1000,-5,5);

  cout << "Initialize fitter..." << endl;
  init_fitter();

  for(int i=0; i < kntoymc; i++) {  

     cout << "* Generating toy-MC experiment: " << i << endl;
     double input_mean = true_mean;
     double input_dip  = true_dip;
#if SYST_FLUCT_MODE == 0
     double input_a1 = true_a1;
     double input_a2 = true_a2;
     double input_a3 = true_a3;
     double input_a4 = true_a4;
#elif SYST_FLUCT_MODE == 1
     TRandom3 rand(0);
     double input_a1 = rand.Gaus(true_a1,ka1err);
     double input_a2 = rand.Gaus(true_a2,ka2err);
     double input_a3 = rand.Gaus(true_a3,ka3err);
     double input_a4 = rand.Gaus(true_a4,ka4err);
#elif SYST_FLUCT_MODE == 2
     TVectorD twp = CholeskyGenerateCorrelatedParams(covL, bf);
     double input_a1 = twp[0];
     double input_a2 = twp[1];
     double input_a3 = twp[2];
     double input_a4 = twp[3];
#endif

     //reset not tweaked parameters to true
#ifndef FLUCT_NOT_FIT
#if DONT_FIT_ALL >= 1
     input_a1 = true_a1;
#if DONT_FIT_ALL >= 2
     input_a2 = true_a2;
#endif
#endif
#endif

     TH1D * data = mock_data(input_mean, input_dip, input_a1, input_a2, input_a3, input_a4);
     set_fit_data(data);

     cout << "Fitting..." << endl;
     double bestfit_mean = -1;
     double bestfit_dip  = -1;
     double bestfit_a1   = -1;
     double bestfit_a2   = -1;
     double bestfit_a3   = -1;
     double bestfit_a4   = -1;
     double error_mean = -1;
     double error_dip  = -1;
     double error_a1   = -1;
     double error_a2   = -1;
     double error_a3   = -1;
     double error_a4   = -1;
     int status_code = fit(data, bestfit_mean, bestfit_dip, bestfit_a1, bestfit_a2, bestfit_a3, bestfit_a4, error_mean, error_dip, error_a1, error_a2, error_a3, error_a4);

     cout << "Best-fit parameters: "
          << "mean = "  << bestfit_mean 
          << ", dip = " << bestfit_dip
          << ", a1 = "  << bestfit_a1 
          << ", a2 = "  << bestfit_a2 
          << ", a3 = "  << bestfit_a3 
          << ", a4 = "  << bestfit_a4 << endl
	  << "Best-fit errors:     "
          << "mean = "  << error_mean 
          << ", dip = " << error_dip
          << ", a1 = "  << error_a1 
          << ", a2 = "  << error_a2 
          << ", a3 = "  << error_a3 
          << ", a4 = "  << error_a4 << endl;


     TH1D * bestfit = pdf(bestfit_mean, bestfit_dip, bestfit_a1, bestfit_a2, bestfit_a3, bestfit_a4);

     //chisq
     double bestfit_chisq  = chisq(data, bestfit_mean, bestfit_dip, bestfit_a1, bestfit_a2, bestfit_a3, bestfit_a4, false);
     //only chisq penalty term
     double bestfit_chisqp = chisq(data, bestfit_mean, bestfit_dip, bestfit_a1, bestfit_a2, bestfit_a3, bestfit_a4, true);

     bestfit_params->Fill(status_code, bestfit_mean, bestfit_dip, bestfit_a1, bestfit_a2, bestfit_a3, bestfit_a4, bestfit_chisq, bestfit_chisqp, error_mean, error_dip, error_a1, error_a2, error_a3, error_a4);

     data   ->Write(Form("data_%d",   i));
     bestfit->Write(Form("bestfit_%d",i));

     throw_a1->Fill(input_a1);
     throw_a2->Fill(input_a2);
     throw_a3->Fill(input_a3);
     throw_a4->Fill(input_a4);

     pull_mean->Fill((bestfit_mean - input_mean) / (error_mean));
     pull_dip-> Fill((bestfit_dip  - input_dip)  / (error_dip) );
     /*
     if(error_a1 < ka1err) {
       pull_a1   ->  Fill((bestfit_a1   - input_a1)   / sqrt(- error_a1 + ka1err));
       pullalt_a1->  Fill((bestfit_a1   - input_a1)   / sqrt(error_a1));
     }
     if(error_a2 < ka2err) {
       pull_a2   ->  Fill((bestfit_a2   - input_a2)   / sqrt(- error_a2 + ka2err));
       pullalt_a2->  Fill((bestfit_a2   - input_a2)   / sqrt(error_a2));
     }
     if(error_a3 < ka3err) {
       pull_a3   ->  Fill((bestfit_a3   - input_a3)   / sqrt(- error_a3 + ka3err));
       pullalt_a3->  Fill((bestfit_a3   - input_a3)   / sqrt(error_a3));
     }
     if(error_a4 < ka4err) {
       pull_a4   ->  Fill((bestfit_a4   - input_a4)   / sqrt(- error_a4 + ka4err));
       pullalt_a4->  Fill((bestfit_a4   - input_a4)   / sqrt(error_a4));
     }
     */
     if(error_a1 < ka1err) {
       pull_a1   ->  Fill((bestfit_a1   - input_a1)   / sqrt(- pow(error_a1,2) + pow(ka1err,2)));
       pullalt_a1->  Fill((bestfit_a1   - input_a1)   / (error_a1));
     }
     if(error_a2 < ka2err) {
       pull_a2   ->  Fill((bestfit_a2   - input_a2)   / sqrt(- pow(error_a2,2) + pow(ka2err,2)));
       pullalt_a2->  Fill((bestfit_a2   - input_a2)   / (error_a2));
     }
     if(error_a3 < ka3err) {
       pull_a3   ->  Fill((bestfit_a3   - input_a3)   / sqrt(- pow(error_a3,2) + pow(ka3err,2)));
       pullalt_a3->  Fill((bestfit_a3   - input_a3)   / (error_a3));
     }
     if(error_a4 < ka4err) {
       pull_a4   ->  Fill((bestfit_a4   - input_a4)   / sqrt(- pow(error_a4,2) + pow(ka4err,2)));
       pullalt_a4->  Fill((bestfit_a4   - input_a4)   / (error_a4));
     }
  }//loop over toys

  gStyle->SetOptFit(1);

  TCanvas * c_mean = new TCanvas(); pull_mean->Draw(); pull_mean->Fit("gaus"); c_mean->SaveAs("mean.png");
  TCanvas * c_dip  = new TCanvas(); pull_dip ->Draw(); pull_dip ->Fit("gaus"); c_dip ->SaveAs("dip.png");
  //TCanvas * c_a1   = new TCanvas(); pull_a1  ->Draw(); pull_a1  ->Fit("gaus"); c_a1  ->SaveAs("a1.png");
  //TCanvas * c_a2   = new TCanvas(); pull_a2  ->Draw(); pull_a2  ->Fit("gaus"); c_a2  ->SaveAs("a2.png");
  //TCanvas * c_a3   = new TCanvas(); pull_a3  ->Draw(); pull_a3  ->Fit("gaus"); c_a3  ->SaveAs("a3.png");
  //TCanvas * c_a4   = new TCanvas(); pull_a4  ->Draw(); pull_a4  ->Fit("gaus"); c_a4  ->SaveAs("a4.png");
  TCanvas * calt_a1   = new TCanvas(); pullalt_a1  ->Draw(); pullalt_a1  ->Fit("gaus"); calt_a1  ->SaveAs("a1alt.png");
  TCanvas * calt_a2   = new TCanvas(); pullalt_a2  ->Draw(); pullalt_a2  ->Fit("gaus"); calt_a2  ->SaveAs("a2alt.png");
  TCanvas * calt_a3   = new TCanvas(); pullalt_a3  ->Draw(); pullalt_a3  ->Fit("gaus"); calt_a3  ->SaveAs("a3alt.png");
  TCanvas * calt_a4   = new TCanvas(); pullalt_a4  ->Draw(); pullalt_a4  ->Fit("gaus"); calt_a4  ->SaveAs("a4alt.png");
  TCanvas * cthrow_a1   = new TCanvas(); throw_a1  ->Draw(); throw_a1  ->Fit("gaus"); cthrow_a1  ->SaveAs("throwa1.png");
  TCanvas * cthrow_a2   = new TCanvas(); throw_a2  ->Draw(); throw_a2  ->Fit("gaus"); cthrow_a2  ->SaveAs("throwa2.png");
  TCanvas * cthrow_a3   = new TCanvas(); throw_a3  ->Draw(); throw_a3  ->Fit("gaus"); cthrow_a3  ->SaveAs("throwa3.png");
  TCanvas * cthrow_a4   = new TCanvas(); throw_a4  ->Draw(); throw_a4  ->Fit("gaus"); cthrow_a4  ->SaveAs("throwa4.png");

  TH1D * meansummary = new TH1D("meansum","meansum",6,0,6); 
  TH1D * sigmasummary = new TH1D("sigmasum","sigmasum",6,0,6);
  TF1 *fit = 0;
  int ifit = 1;
  fit = pull_mean->GetFunction("gaus");
  meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
  sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2)); ifit++;
  fit = pull_dip->GetFunction("gaus");
  meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
  sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2)); ifit++;
#if DONT_FIT_ALL == 0
  fit = pullalt_a1->GetFunction("gaus");
  meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
  sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2)); ifit++;
#else
  meansummary->SetBinContent(ifit, 0);
  sigmasummary->SetBinContent(ifit, 1); ifit++;
#endif
#if DONT_FIT_ALL <= 1
  fit = pullalt_a2->GetFunction("gaus");
  meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
  sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2)); ifit++;
#else
  meansummary->SetBinContent(ifit, 0);
  sigmasummary->SetBinContent(ifit, 1); ifit++;
#endif
  fit = pullalt_a3->GetFunction("gaus");
  meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
  sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2)); ifit++;
  fit = pullalt_a4->GetFunction("gaus");
  meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
  sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2)); ifit++;
  TCanvas * csum = new TCanvas(); csum->Divide(1,2); csum->cd(1); meansummary->Draw(); csum->cd(2); sigmasummary->Draw(); csum->SaveAs("1.png");

  bestfit_params->Write();
  pull_mean->Write();
  pull_dip ->Write();
  pull_a1  ->Write();
  pull_a2  ->Write();
  pull_a3  ->Write();
  pull_a4  ->Write();
  pullalt_a1  ->Write();
  pullalt_a2  ->Write();
  pullalt_a3  ->Write();
  pullalt_a4  ->Write();
  throw_a1  ->Write();
  throw_a2  ->Write();
  throw_a3  ->Write();
  throw_a4  ->Write();
  meansummary->Write();
  sigmasummary->Write();
  f.Close();

  cout << "Done" << endl;

  return 0;
}
//.............................................................................


