//____________________________________________________________________________ 
/*

 >> ChangeLog:

 For changes up to version 2010_v3.0.0 (last implementation of the 2010-era analysis)
 see code in the corresponding SVN branch.

 Changes since 2011_v0.9.0:

 @

*/ 
//____________________________________________________________________________

#define __show_fit_init_info__ 
#define __show_fit_init_debug__ 
#define __show_fit_result__ 
#define __show_fit_steps__ 
#define __show_fit_cpu_time__

#include <cassert> 
#include <cstdlib> 
#include <iostream> 
#include <sstream>

#include <TMath.h> 
#include <TGraph.h> 
#include <TH2D.h> 
#include <THStack.h> 
#include <TMarker.h> 
#include <TLegend.h> 
#include <TString.h> 
#include <TGaxis.h> 
#include <TCanvas.h> 
#include <TFile.h>
#include <TSystem.h>
#include <TStopwatch.h>

#include <Minuit2/MnUserParameters.h>
#include <Minuit2/MnUserCovariance.h>
#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MnApplication.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/MnHesse.h>
#include <Minuit2/MnMinos.h>
#include <Minuit2/MnPrint.h>

#include "FitterKernelMinuit2.h" 
#include "SuperkObservation.h" 
#include "SuperkExpectation.h" 
#include "FitInputs.h" 
#include "FitParamList.h" 
#include "Utilities.h"
#include "Conventions.h"

using std::cerr; 
using std::cout; 
using std::endl; 
using std::ostringstream;

using namespace t2kana; 
using namespace t2kana::conv;
   
//____________________________________________________________________________ 
FitterKernelMinuit2::FitterKernelMinuit2() :
FitterABC()
{
  this->Init();
}
//____________________________________________________________________________ 
FitterKernelMinuit2::FitterKernelMinuit2(FitInputs * fit_inputs, Option_t * option) :
FitterABC()
{
  this->Init(fit_inputs, option);
}
//____________________________________________________________________________ 
FitterKernelMinuit2::~FitterKernelMinuit2() 
{
  if(fResults)    { delete fResults;    }
  if(fUserParams) { delete fUserParams; }
}
//____________________________________________________________________________ 
void FitterKernelMinuit2::Init(void) 
{
  fFitInputs      = 0;
  fOption         = "";
  fNFreeFitParams = 0;
  fNFitBins       = 0;
  fResults        = 0;
  fUserParams     = 0;
}
//____________________________________________________________________________ 
void FitterKernelMinuit2::Init(FitInputs * fit_inputs, Option_t * option) 
{
#ifdef __show_fit_init_info__
  cout << "FitterKernelMinuit2/Init >> Fitter initialization..." << endl;
#endif

  this->Init();

  fResults = new TObjArray();

  fFitInputs = fit_inputs;
  if(!fFitInputs) {
     cerr << "FitterKernelMinuit2/Init >> ** Null fit inputs!" << endl;
     exit(1);
  }

  FitParamList & fit_params = fFitInputs->FitParams();
  int param_list_size = fit_params.Size();

#ifdef __show_fit_init_info__
  cout << "FitterKernelMinuit2/Init >> Fit parameter list size = " 
       << param_list_size << endl;
#endif
  
  // set params
  ROOT::Minuit2::MnUserParameters * fUserParams = new ROOT::Minuit2::MnUserParameters();
  fit_params.SetTo(fUserParams);

#ifdef __show_fit_init_info__
  cout << "FitterKernelMinuit2/Init >> Constructing covariance matrix " << endl;
#endif


#ifdef __show_fit_init_debug__
  cout << "FitterKernelMinuit2/Init >> Covariance matrix:" << endl;
  //covariance_matrix.Print();
#endif

  // calculate degrees of freedom

  fNFitBins       = kNErecoBins1Rmu;
  fNFreeFitParams = fit_params.GetNumberFreeParameters();

#ifdef __show_fit_init_info__
  cout << "FitterKernelMinuit2/Init >> NDoF = " << this->NDOF()
       << " (Number of bins = " << this->NBinsInFit()
       << ", Number of free parameters = " << this->NFitParams() << ")"
       << endl;
#endif
}
//____________________________________________________________________________ 
bool FitterKernelMinuit2::DoTheFit(void) 
{
#ifdef __show_fit_cpu_time__
  TStopwatch w;
  w.Start();
#endif

  // Run the fitter once - we have the option here to check the status code
  // and if the fitter failed, re-seed the fitter and try again (not needed
  // so far).
  bool ok = this->RunOnce();

#ifdef __show_fit_cpu_time__
  w.Stop();
  cout << "FitterKernelMinuit2/DoTheFit >> Fit performed in " 
       << w.CpuTime() << " sec." << endl;
#endif

  return ok;
}
//____________________________________________________________________________ 
bool FitterKernelMinuit2::RunOnce(void) 
{
  // a few arrays to store parameters
  const int nparamsmax = knppmax;
  double bestfit_arr   [nparamsmax] = {0.};
  double err_migrad_cov[nparamsmax] = {0.};
  double err_migrad_err[nparamsmax] = {0.};
  double err_hesse_cov [nparamsmax] = {0.};
#ifdef __get_minos_errors__
  double err_minosplus [nparamsmax] = {0.};
  double err_minosminus[nparamsmax] = {0.};
#endif

  // Set fitter strategy
  // In general, low values of level mean fewer function calls and high values mean more reliable minimization.
  // Currently allowed values are 0 (low), 1 (default), and 2 (high).
  ROOT::Minuit2::MnStrategy strategy(1);

  // Get the fitter & pass it the inputs
#ifdef __show_fit_steps__
  cout << "FitterKernelMinuit2/RunOnce >> Getting the fitter & giving it inputs (data & pdf)..." << endl;
#endif

  ROOT::Minuit2::MnMigrad migrad(*fFitInputs, *fUserParams);
  
  // Run MINUIT's minimization step
#ifdef __show_fit_steps__
  cout << "FitterKernelMinuit2/RunOnce >> Running MIGRAD..." << endl;
#endif

  // max number of function calls
  // if zero, max_calls = 200 + 100 * npar + 5 * npar**2
  int max_calls = 0;

  // value used for terminating iteration procedure.
  // For example, MIGRAD will stop iterating when edm (expected distance from minimum) will be:
  // edm < tolerance * 10**-3
  double tolerance = 0.1;

  //Minimise!
  ROOT::Minuit2::FunctionMinimum bestfit = migrad(max_calls, tolerance);

  // Done! Copy current --> best-fit values
  fFitInputs->FitParams().CopyCurrToBestFitValues();
  std::vector<double> bestfit_vec;
  bestfit_vec = bestfit.UserParameters().Params();
  for(unsigned int ip = 0; ip < bestfit_vec.size(); ip++) {
    bestfit_arr[ip] = bestfit_vec[ip];
    cout << "Bestfit for par#" << ip << " = " << bestfit_arr[ip] << endl;
  }

#ifdef __show_fit_steps__
  cout << "FitterKernelMinuit2/RunOnce >> Getting fit success status..." << endl;
#endif
  //Check if fit succeeded
  if(!bestfit.IsValid()){
    std::cout << "Fit Failed" << std::endl;
    return false;
  }
  else {
    std::cout << "Fit succeeded" << std::endl;
    fStatusCode = 3;
  }

  // Get errors
#ifdef __show_fit_steps__
  cout << "FitterKernelMinuit2/RunOnce >> Getting errors..." << endl;
#endif

  // migrad covariance
  ROOT::Minuit2::MnUserCovariance cov_migrad(bestfit.UserCovariance());
  for(unsigned int ip = 0; ip < bestfit_vec.size(); ip++) {
    err_migrad_cov[ip] = cov_migrad(ip,ip);
    cout << "Migrad covariance for par#" << ip << " = " << bestfit_arr[ip] << " +/- " << err_migrad_cov[ip] << endl;
  }

  // migrad error
  std::vector<double> error_vec;
  error_vec =  bestfit.UserParameters().Errors(); 
  //fFitInputs->FitParams().SetMigradErrors(error_vec);
  for(unsigned int ip = 0; ip < bestfit_vec.size(); ip++) {
    err_migrad_err[ip] = error_vec[ip];
    cout << "Migrad error for par#" << ip << " = " << bestfit_arr[ip] << " +/- " << err_migrad_err[ip] << endl;
  }

  //Get Hesse errors
#ifdef __show_fit_steps__
  cout << "FitterKernelMinuit2/RunOnce >> Getting Hesse errors..." << endl;
#endif
  // argument is strategy
  ROOT::Minuit2::MnHesse hesse(2);
  hesse(*fFitInputs, bestfit);

  // hesse covariance
  ROOT::Minuit2::MnUserCovariance cov_hesse(bestfit.UserCovariance());
  for(unsigned int ip = 0; ip < bestfit_vec.size(); ip++) {
    err_hesse_cov[ip] = cov_hesse(ip,ip);
    cout << "Migrad covariance for par#" << ip << " = " << bestfit_arr[ip] << " +/- " << err_hesse_cov[ip] << endl;
  }

#ifdef __get_minos_errors__
#ifdef __show_fit_steps__
  cout << "FitterKernelMinuit2/RunOnce >> Getting Minos errors..." << endl;
#endif

  //Get Minos Errors
  ROOT::Minuit2::MnMinos minos(*migrad,bestfit);
  for(int ip = 0; ip < fFitInputs->FitParams().Size(); ip++) {
      double up = minos_nd280.Upper(ip);
      double lo = minos_nd280.Lower(ip);
      err_minosminus[ip] = pow(lo,2);
      err_minosplus [ip] = pow(up,2);
      cout << "FitterKernelMinuit2/RunOnce >> par#" << ip << " = " <<  << " +" << up << " -" << lo << " (minos)" << endl;
  }
#endif

  //Get the chi^2
  cout << "fcn = " << bestfit.Fval() << endl;
  cout << "chi2 = " << fFitInputs->Chisquare() << endl;

  // Print fit result
#ifdef __show_fit_result__
#ifdef __show_fit_steps__
  cout << "FitterKernelMinuit2/RunOnce >> Printing fit results..." << endl;
#endif
  cout << bestfit << endl;
#endif

#ifdef __show_fit_steps__
  cout << "FitterKernelMinuit2/RunOnce >> Done" << endl;
#endif

  return (fStatusCode == 3);
}
//____________________________________________________________________________ 
