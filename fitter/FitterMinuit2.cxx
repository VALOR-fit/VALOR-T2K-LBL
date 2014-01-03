//____________________________________________________________________________ 
/*

 \copyright (c) 2010-2013, VaLOR group.

*/ 
//____________________________________________________________________________

#define __show_fit_init_info__ 
#define __show_fit_cov_matrix_debug__
#define __show_fit_result__ 
#define __show_fit_steps__ 
#define __show_fit_steps_debug__
#define __show_fit_cpu_time__
#define __get_hesse_errors__
#define __get_minos_errors__
//#define __give_fitter_covariance_matrix__

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

#include "FitterMinuit2.h" 
#include "SuperkObservation.h" 
#include "SuperkExpectation.h" 
#include "FitInputs.h" 
#include "FitParamList.h" 
#include "Utilities.h"
#include "Conventions.h"
#include "BANFF.h"
#include "FSI.h"
#include "SuperkDetector.h"

using std::cerr; 
using std::cout; 
using std::endl; 
using std::ostringstream;

using namespace t2kana; 
using namespace t2kana::conv;
   
//____________________________________________________________________________ 
FitterMinuit2::FitterMinuit2() :
FitterABC()
{
  this->Init();
}
//____________________________________________________________________________ 
FitterMinuit2::FitterMinuit2(FitInputs * fit_inputs, Option_t * option) :
FitterABC()
{
  this->Init(fit_inputs, option);
}
//____________________________________________________________________________ 
FitterMinuit2::~FitterMinuit2() 
{
  if(fCovMtrx)    { delete fCovMtrx;    }
}
//____________________________________________________________________________ 
void FitterMinuit2::Init(void) 
{
  fFitInputs      = 0;
  fOption         = "";
  fNFreeFitParams = 0;
  fNFitBins       = 0;
  fStatusCode     = -999;
  fFitTries       = 0;
}
//____________________________________________________________________________ 
void FitterMinuit2::Init(FitInputs * fit_inputs, Option_t * /*option*/) 
{
#ifdef __show_fit_init_info__
  cout << "FitterMinuit2/Init >> Fitter initialization..." << endl;
#endif

  this->Init();

  fFitInputs = fit_inputs;
  if(!fFitInputs) {
     cerr << "FitterMinuit2/Init >> ** Null fit inputs!" << endl;
     exit(1);
  }

  FitParamList & fit_params = fFitInputs->FitParams();
  int param_list_size = fit_params.Size();

#ifdef __show_fit_init_info__
  cout << "FitterMinuit2/Init >> Fit parameter list size = " 
       << param_list_size << endl;
#endif
  
  // set params
  fit_params.SetTo(fUserParams);
#ifdef __show_fit_init_info__
  cout << fUserParams << endl;
#endif

  // calculate degrees of freedom
  fNFitBins = 0;
  for(int is = 0; is < kns; is++) {
    SuperkSample_t s = kSKSamples[is];
    if(s == kSKsFCFV1RingMuLike) fNFitBins += kNErecoBins1Rmu;
    if(s == kSKsFCFV1RingELike ) fNFitBins += kNErecoBins1Re;
    if(s == kSKsRHCFCFV1RingMuLike) fNFitBins += kRHCNErecoBins1Rmu;
    if(s == kSKsRHCFCFV1RingELike ) fNFitBins += kRHCNErecoBins1Re;
  }
  fNFreeFitParams = fit_params.GetNumberFreeParameters();

#ifdef __show_fit_init_info__
  cout << "FitterMinuit2/Init >> NDoF = " << this->NDOF()
       << " (Number of bins = " << this->NBinsInFit()
       << ", Number of free parameters = " << this->NFitParams() << ")"
       << endl;
#endif
}
//____________________________________________________________________________ 
void FitterMinuit2::ConstructCovMatrix(void)
{
#ifdef __show_fit_cov_matrix_debug__
  cout << "FitterMinuit2/ConstructCovMatrix >> Constructing covariance matrix" << endl;
#endif

  //Get number of free parameters
  int nfree = fNFreeFitParams;
  TMatrixD total_covariance(nfree, nfree);
  total_covariance.Zero();

  //Loop over parameters, find the PhysParamId, and get the corresponding error (if free)
  int ipf = 0;
  for(int ip = 0; ip < fFitInputs->FitParams().Size(); ip++) {
    string name = fUserParams.Name((unsigned int)ip);
    PhysParamId_t p = util::PhysParamIdFromLongString(name);
    if(fFitInputs->FitParams().IsFixed(p))
      continue;
    if(!util::IsSystematic(p))
      total_covariance(ipf, ipf) = fFitInputs->FitParams().StepSize(p);
    //now find the appropriate error for the parameter
    //BANFF
    else if(util::ArrIdx::Instance()->BANFFSystIdx(p) >= 0) {
      if(util::ArrIdx::Instance()->BANFFSystIdx(p) == 0) {
	const TMatrixD & sub_cov = fFitInputs->Expectation().GetBANFF()->CovMatrix();
	total_covariance.SetSub(ipf, ipf, sub_cov);
      }
    }
    //SKDet
#if FSI_BINS != 3
    else if(p == kNuisSuperkEff1RmuCCQErc0u) {
      const TMatrixD & sub_cov = fFitInputs->Expectation().GetSuperkDetector()->Eff1RmuCCQERCCovarRelErr();
      total_covariance.SetSub(ipf, ipf, sub_cov);
    }
#endif
    else if(p == kNuisSuperkEnergyScale)
      total_covariance(ipf, ipf) = 1;
#if FSI_BINS !=3
    else if(p == kNuisSuperkEff1RmuCCQEoth)
      total_covariance(ipf, ipf) = 1;
    else if(p == kNuisSuperkEff1RmuCCnonQE)
      total_covariance(ipf, ipf) = 1;
    else if(p == kNuisSuperkEff1RmuNC)
      total_covariance(ipf, ipf) = 1;
    else if(p == kNuisSuperkEff1RmuNueCC)
      total_covariance(ipf, ipf) = 1;
#endif
    //XSec
    else if(p == kNuisXSecFermiMom)
      total_covariance(ipf, ipf) = 1;
    else if(p == kNuisXSecSpecFunc)
      total_covariance(ipf, ipf) = 1;
    else if(p == kNuisXSecWshape)
      total_covariance(ipf, ipf) = 1;
    else if(p == kNuisXSecShapeCCoth)
      total_covariance(ipf, ipf) = 1;
    else if(p == kNuisXSecPiLessDelta)
      total_covariance(ipf, ipf) = 1;
    else if(p == kNuisXSecBindE)
      total_covariance(ipf, ipf) = 1;
    else if(p == kNuisXSecNormCCcoh)
      total_covariance(ipf, ipf) = 1;
    else if(p == kNuisXSecNormNC1piC)
      total_covariance(ipf, ipf) = 1;
    else if(p == kNuisXSecNormNCoth)
      total_covariance(ipf, ipf) = 1;
    else if(p == kNuisXSecNormNuEToNuMu)
      total_covariance(ipf, ipf) = 1;
    else if(p == kNuisXSecNormAntiNuToNu)
      total_covariance(ipf, ipf) = 1;
    //FSI
    else if(util::ArrIdx::Instance()->FSISystIdx(p) >= 0) {
      if(util::ArrIdx::Instance()->FSISystIdx(p) == 0) {
	const TMatrixD & sub_cov = fFitInputs->Expectation().GetFSI()->CovMatrix();
	total_covariance.SetSub(ipf, ipf, sub_cov);
      }
    }
    else
      continue;
    //increment position in the covariance matrix we're constructing
    ipf++;
  }//end loop over FitParamList

  if(fCovMtrx)
    delete fCovMtrx;
  fCovMtrx = new TMatrixD(total_covariance);

#ifdef __show_fit_cov_matrix_debug__
  cout << "FitterMinuit2/ConstructCovMatrix >> Covariance matrix:" << endl;
  total_covariance.Print();
#endif

}
//____________________________________________________________________________ 
int FitterMinuit2::FitStatusCode(Option_t * /*option*/)
{
  return fStatusCode;
}
//____________________________________________________________________________
bool FitterMinuit2::DoTheFit(void) 
{
#ifdef __show_fit_cpu_time__
  TStopwatch w;
  w.Start();
#endif
  fFitTime = 0;

  // Run the fitter once - we have the option here to check the status code
  // and if the fitter failed, re-seed the fitter and try again (not needed
  // so far).
  bool ok = this->RunOnce();

#ifdef __show_fit_cpu_time__
  w.Stop();
  cout << "FitterMinuit2/DoTheFit >> Fit performed in " 
       << w.CpuTime() << " sec." << endl;
  fFitTime = w.CpuTime();
#endif

  return ok;
}
//____________________________________________________________________________ 
bool FitterMinuit2::RunOnce(void) 
{
  // a few arrays to store parameters
  unsigned int nfreeparams = (unsigned int)fNFreeFitParams;
  unsigned int nparams     = (unsigned int)fFitInputs->FitParams().Size();
  vector<double> bestfit_free;
  vector<double> err_migrad_cov;
  vector<double> err_migrad_err;
#ifdef __get_hesse_errors__
  vector<double> err_hesse_cov;
  vector<double> err_hesse_err;
#endif
#ifdef __get_minos_errors__
  vector<double> err_minosplus;
  vector<double> err_minosminus;
#endif

  // Set fitter strategy
  // In general, low values of level mean fewer function calls and high values mean more reliable minimization.
  // Currently allowed values are 0 (low), 1 (default), and 2 (high).
  ROOT::Minuit2::MnStrategy strategy(1);

#ifdef __give_fitter_covariance_matrix__
  // We have a TMatrixD covariance matrix, need to translate it into
  // the Minuit2::MnUserCovariance
  this->ConstructCovMatrix();
  ROOT::Minuit2::MnUserCovariance input_cov(nfreeparams);
  TMatrixD covmtrx(*fCovMtrx);
  for(unsigned int i = 0; i < nfreeparams; i++) {
    for(unsigned int j = 0; j < nfreeparams; j++) {
      input_cov(i,j) = covmtrx(i,j);
    }
  }
#ifdef __show_fit_cov_matrix_debug__
  cout << "FitterMinuit2/RunOnce >> Input covariance matrix:" << endl;
  cout << input_cov << endl;
#endif
#else
  fCovMtrx = new TMatrixD(); // need to do this to make the destructor work. confusing...
#endif

  // Get the fitter & pass it the inputs
#ifdef __show_fit_steps__
  cout << "FitterMinuit2/RunOnce >> Getting the fitter & giving it inputs (data & pdf)..." << endl;
#endif

#ifdef __give_fitter_covariance_matrix__
  ROOT::Minuit2::MnMigrad migrad(*fFitInputs, fUserParams, input_cov);
#else
  ROOT::Minuit2::MnMigrad migrad(*fFitInputs, fUserParams);
#endif

  // Run MINUIT's minimization step
#ifdef __show_fit_steps__
  cout << "FitterMinuit2/RunOnce >> Running MIGRAD..." << endl;
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

#ifdef __show_fit_steps__
  cout << "FitterMinuit2/RunOnce >> MIGRAD minimisation complete" << endl;
#endif

#ifdef __show_fit_steps__
  cout << "FitterMinuit2/RunOnce >> Getting fit success status..." << endl;
#endif

  // Done! Copy current --> best-fit values
  fFitInputs->FitParams().CopyCurrToBestFitValues();

  //Check if fit succeeded
  GetStatusCode(bestfit);
  if(fStatusCode == 9) {
    cerr << "FitterMinuit2/RunOnce >> Fit failed. Returning..." << endl;
    return false;
  }

  //get the best fit vector, for use in this function
  std::vector<ROOT::Minuit2::MinuitParameter> bestfit_params = bestfit.UserParameters().Parameters();
  for(unsigned int ip = 0; ip < bestfit_params.size(); ip++) {
    if(!bestfit_params[ip].IsFixed()) {
      bestfit_free.push_back(bestfit_params[ip].Value());
    }
#ifdef __show_fit_steps_debug__
    cout << "FitterMinuit2/RunOnce >> Bestfit for par#" << ip << " = " << bestfit_params[ip].Value() << endl;
#endif
  }

  // Get errors
#ifdef __show_fit_steps__
  cout << endl << "FitterMinuit2/RunOnce >> Getting Migrad errors..." << endl;
#endif

  // migrad covariance
  if(bestfit.HasCovariance()) {
    ROOT::Minuit2::MnUserCovariance cov_migrad(bestfit.UserCovariance());
#ifdef __show_fit_steps_debug__
    cout << cov_migrad << endl;
#endif
    for(unsigned int ip = 0; ip < nfreeparams; ip++) {
      double tmp = sqrt(cov_migrad(ip,ip));
      err_migrad_cov.push_back(tmp);
#ifdef __show_fit_steps_debug__
      cout << "FitterMinuit2/RunOnce >> Migrad covariance for par#" << ip << " = " << bestfit_free[ip] << " +/- " << err_migrad_cov[ip] << endl;
#endif
    }
  }

  // migrad error
  err_migrad_err =  bestfit.UserParameters().Errors(); 
  for(unsigned int ip = 0; ip < nparams; ip++) {
#ifdef __show_fit_steps_debug__
    cout << "FitterMinuit2/RunOnce >> Migrad error for par#" << ip << " = " << bestfit_free[ip] << " +/- " << err_migrad_err[ip] << endl;
#endif
  }
  fFitInputs->FitParams().SetMigradErrors(err_migrad_err);

#ifdef __get_hesse_errors__
  //Get Hesse errors
#ifdef __show_fit_steps__
  cout << endl << "FitterMinuit2/RunOnce >> Getting Hesse errors..." << endl;
#endif
  // argument is strategy
  ROOT::Minuit2::MnHesse hesse(2);
  // perform hesse error calculation
  hesse(*fFitInputs, bestfit);
  //Check if fit succeeded
  GetStatusCode(bestfit);

  // hesse covariance
  if(bestfit.HasCovariance()) {
    ROOT::Minuit2::MnUserCovariance cov_hesse(bestfit.UserCovariance());
#ifdef __show_fit_steps_debug__
    cout << cov_hesse << endl;
#endif
    for(unsigned int ip = 0; ip < nfreeparams && ip < cov_hesse.Nrow(); ip++) {
      err_hesse_cov.push_back(sqrt(cov_hesse(ip,ip)));
#ifdef __show_fit_steps_debug__
      cout << "FitterMinuit2/RunOnce >> Hesse covariance for par#" << ip << " = " << bestfit_free[ip] << " +/- " << err_hesse_cov[ip] << endl;
#endif
    }
  }

  // hesse error
  err_hesse_err =  bestfit.UserParameters().Errors(); 
  fFitInputs->FitParams().SetHesseErrors(err_hesse_err);
  for(unsigned int ip = 0; ip < nparams; ip++) {
#ifdef __show_fit_steps_debug__
    cout << "FitterMinuit2/RunOnce >> Hesse error for par#" << ip << " = " << bestfit_free[ip] << " +/- " << err_hesse_err[ip] << endl;
#endif
  }
  GetStatusCode(bestfit);
#endif

#ifdef __get_minos_errors__
#ifdef __show_fit_steps__
  cout << endl << "FitterMinuit2/RunOnce >> Getting Minos errors..." << endl;
#endif

  //Get Minos Errors
  ROOT::Minuit2::MnMinos minos(*fFitInputs,bestfit);
  for(unsigned int ip = 0; ip < nparams; ip++) {
    if(bestfit_params[ip].IsFixed())
      continue;
    double up = minos.Upper(ip);
    double lo = minos.Lower(ip);
    err_minosminus.push_back(lo);
    err_minosplus. push_back(up);
#ifdef __show_fit_steps_debug__
    cout << "FitterMinuit2/RunOnce >> par#" << ip << " = " << bestfit_params[ip].Value() << " +" << up << " -" << lo << " (minos)" << endl;
#endif
  }
  fFitInputs->FitParams().SetMinosErrors(err_minosplus, err_minosminus);
#endif

  fFitInputs->FitParams().SetFitTime ( fFitTime );
  fFitInputs->FitParams().SetBestFitChisq ( fFitInputs->BestFitChisquare());


  // Print fit result
#ifdef __show_fit_result__
#ifdef __show_fit_steps__
  cout << endl << "FitterMinuit2/RunOnce >> Printing fit results..." << endl;
#endif
  cout << "FitterMinuit2/RunOnce >> fcn =  " << bestfit.Fval() << endl;
  cout << "FitterMinuit2/RunOnce >> chi2 = " << fFitInputs->Chisquare() << endl;
  cout << bestfit << endl;
  cout << "FitterMinuit2/RunOnce >> fit status code = " << fStatusCode << endl;
#endif

#ifdef __show_fit_steps__
  cout << "FitterMinuit2/RunOnce >> Done" << endl;
#endif

  return (fStatusCode == 0);
}
//____________________________________________________________________________ 
void FitterMinuit2::GetStatusCode(ROOT::Minuit2::FunctionMinimum bestfit) 
{
  if(!bestfit.IsValid()) {
#ifdef __show_fit_steps__
    cout << "FitterMinuit2/GetStatusCode >> Minuit failed." << endl;
#endif
    fStatusCode = 9;
  }
  else if (bestfit.HasReachedCallLimit()) {
    #ifdef __show_fit_steps__
    cout << "FitterMinuit2/GetStatusCode >> Minuit reached the specified call limit (maxcalls)." << endl;
#endif
    fStatusCode = 1;
  }
  else if (bestfit.IsAboveMaxEdm()) {
#ifdef __show_fit_steps__
    cout << "FitterMinuit2/GetStatusCode >> Function value is above the specified estimated distance to the minimum (edm)." << endl;
#endif
    fStatusCode = 2;
  }
  else if (!bestfit.HasPosDefCovar()) {
#ifdef __show_fit_steps__
    cout << "FitterMinuit2/GetStatusCode >> Covariance is not positive definite." << endl;
#endif
    fStatusCode = 3;
  }
  else if (!bestfit.HasMadePosDefCovar()) {
#ifdef __show_fit_steps__
    cout << "FitterMinuit2/GetStatusCode >> Covariance could not be made positive definite." << endl;
#endif
    fStatusCode = 4;
  }
  else if (!bestfit.HasAccurateCovar()) {
#ifdef __show_fit_steps__
    cout << "FitterMinuit2/GetStatusCode >> Covariance is not accurate." << endl;
#endif
    fStatusCode = 5;
  }
  else if (!bestfit.HasValidCovariance()) {
#ifdef __show_fit_steps__
    cout << "FitterMinuit2/GetStatusCode >> Covariance is not valid." << endl;
#endif
    fStatusCode = 6;
  }
  else if (bestfit.HesseFailed()) {
#ifdef __show_fit_steps__
    cout << "FitterMinuit2/GetStatusCode >> HESSE failed." << endl;
#endif
    fStatusCode = 7;
  }
  else if (!bestfit.HasValidParameters()) {
#ifdef __show_fit_steps__
    cout << "FitterMinuit2/GetStatusCode >> Parameters are not valid." << endl;
#endif
    fStatusCode = 8;
  }
  else {
#ifdef __show_fit_steps__
    cout << "FitterMinuit2/GetStatusCode >> Fit valid" << endl;
#endif
    fStatusCode = 0;
  }
}
//____________________________________________________________________________ 
