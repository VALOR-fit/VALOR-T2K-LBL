//____________________________________________________________________________ 
/*

 \copyright (c) 2010-2013, VaLOR group.

*/ 
//____________________________________________________________________________

//#define __show_fit_init_info__ 
#define __show_fit_result__ 
#define __show_fit_steps__ 
#define __show_fit_cpu_time__

#include <cassert> 
#include <cstdlib> 
#include <iostream> 
#include <sstream>

#include <TMath.h> 
#include <TRandom3.h> 
#include <TVirtualFitter.h> 
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

#include "FitterMinuit.h" 
#include "SuperkObservation.h" 
#include "SuperkExpectation.h" 
#include "FitInputs.h" 
#include "FitParamList.h" 
#include "Utilities.h"
#include "Conventions.h"
#include "Opt.h"
#include "Rnd.h"

using std::cerr; 
using std::cout; 
using std::endl; 
using std::ostringstream;

using namespace t2kana; 
using namespace t2kana::conv;
   
//____________________________________________________________________________ 
FitterMinuit::FitterMinuit() :
FitterABC()
{
  this->Init();
}
//____________________________________________________________________________ 
FitterMinuit::FitterMinuit(FitInputs * fit_inputs, Option_t * option) :
FitterABC()
{
  this->Init(fit_inputs, option);
}
//____________________________________________________________________________ 
FitterMinuit::~FitterMinuit() 
{

}
//____________________________________________________________________________ 
void FitterMinuit::Init(void) 
{
  fFitInputs           = 0;
  fOption              = "";
  fNFreeFitParams      = 0;
  fNFitBins            = 0;
  fCalcHesseErrors     = false;
  fCalcMinosErrors     = false;
//fSimplex             = false;
  fTol                 = 0.1;
  fStr                 = 1;
  fMigradMinStatusCode = -999;
  fMigradErrStatusCode = -999;      
  fHesseErrStatusCode  = -999; 
  fMinosErrStatusCode  = -999; 
  fFitTime             = 0;
}
//____________________________________________________________________________ 
int FitterMinuit::FitStatusCode (Option_t * option) 
{ 
  if(!option) return 9999;

  TString opt(option);

  if ( opt.Contains("migrad_min") ) return fMigradMinStatusCode;
  if ( opt.Contains("migrad_err") ) return fMigradErrStatusCode;
  if ( opt.Contains("hesse_err")  ) return fHesseErrStatusCode;
  if ( opt.Contains("minos_err")  ) return fMinosErrStatusCode;
  
  return -999;
}
//____________________________________________________________________________ 
void FitterMinuit::Init(FitInputs * fit_inputs, Option_t * option) 
{
#ifdef __show_fit_init_info__
  cout << "FitterMinuit/Init >> Fitter initialization..." << endl;
#endif

  this->Init();

  fFitInputs = fit_inputs;
  if(!fFitInputs) {
     cerr << "FitterMinuit/Init >> ** Null fit inputs!" << endl;
     exit(1);
  }

  FitParamList & fit_params = fFitInputs->FitParams();
  int param_list_size = fit_params.Size();

#ifdef __show_fit_init_info__
  cout << "FitterMinuit/Init >> Fit parameter list size = " 
       << param_list_size << endl;
#endif

  fOption = option;

  // calc HESSE or MINOS errors after fit
  fCalcHesseErrors = fOption.Contains("hesse");
  fCalcMinosErrors = fOption.Contains("minos");

  if(fOption.Contains("minuit")) {
    // get fitter
    TVirtualFitter::SetDefaultFitter("Minuit");
    TVirtualFitter * fitter = TVirtualFitter::Fitter(0,param_list_size);
    
    // set print level
    Double_t arglist[5];
    arglist[0] = -1;
    fitter->ExecuteCommand("SET PRINT",arglist,1);
    
    // set params
    fit_params.SetTo(fitter);
    
    // set minimization function
    fitter->SetFCN(FitFunc);            
  }

/*
  if(fOption.Contains("simplex")) {
    TVirtualFitter::SetDefaultFitter("Simplex");
    TVirtualFitter * fitter = TVirtualFitter::Fitter(0,param_list_size);
    
    // set print level              
    Double_t arglist[5];
    arglist[0] = -1;
    fitter->ExecuteCommand("SET PRINT",arglist,1);
    
    // set params                
    fit_params.SetTo(fitter);
    
    // set minimization function     
    fitter->SetFCN(FitFunc);
    fSimplex = true;
  }
*/

  // calculate degrees of freedom  
  fNFitBins = 0;  
  for(int is = 0; is < kns; is++) {
    SuperkSample_t s = kSKSamples[is];
    if(s == kSKsFCFV1RingMuLike) fNFitBins += kNErecoBins1Rmu;
    if(s == kSKsFCFV1RingELike ) fNFitBins += kNErecoBins1Re;
  }  
  fNFreeFitParams = fit_params.GetNumberFreeParameters();
  
#ifdef __show_fit_init_info__
  cout << "FitterMinuit/Init >> NDoF = " << this->NDOF()
       << " (Number of bins = " << this->NBinsInFit()
       << ", Number of free parameters = " << this->NFitParams() << ")"
       << endl;
#endif
  
}
//____________________________________________________________________________ 
bool FitterMinuit::DoTheFit(void) 
{
  // start clock
  TStopwatch w;
  w.Start();

  fFitTime = 0;

  // reset fit time, bestfit chisq and status codes in
  // persistable FitParamList obj
  fFitInputs->FitParams().SetFitTime(0);
  fFitInputs->FitParams().SetBestFitChisq(0);
  fFitInputs->FitParams().ClearStatusCodes();

  // Run the fitter once - we have the option here to check the status code
  // and if the fitter failed, re-seed the fitter and try again
  bool ok = this->RunOnce(false);

  // 4: migrad not converged
  int ntries = 0;
  int maxntries = 10;
  OscAnaType_t ana = Opt::Instance()->gOscAnaTypeChoice;
  if(ana != kAnaNuMuDisappOnly)
    maxntries = 25;
  while(fMigradMinStatusCode == 4 && ntries < maxntries) {
    cout << "Original fit failed; reseeding & trying again. Attempt: " << ntries << endl;
    NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
    switch(omp) {
    case(conv::kNuMixSinsq2Theta):
      if(ana == kAnaNuMuDisappOnly || ana == kAnaCombined)
	fFitInputs->FitParams().SetInitValue(conv::kPhypMixing23, Rnd::Instance()->RndFitter().Uniform(0.80,0.98));
      if(ana == kAnaNuEAppOnly     || ana == kAnaCombined)
	fFitInputs->FitParams().SetInitValue(conv::kPhypMixing13, Rnd::Instance()->RndFitter().Uniform(0.05,0.15));
      break;
    case(conv::kNuMixSinsqTheta):
      if(ana == kAnaNuMuDisappOnly || ana == kAnaCombined)
	if (fFitInputs->FitParams().InitValue(conv::kPhypMixing23)< 0.515)
	  fFitInputs->FitParams().SetInitValue(conv::kPhypMixing23, Rnd::Instance()->RndFitter().Uniform(0.35,0.5));
	else
	  fFitInputs->FitParams().SetInitValue(conv::kPhypMixing23, Rnd::Instance()->RndFitter().Uniform(0.5,0.65));
      if(ana == kAnaNuEAppOnly || ana == kAnaCombined)
	fFitInputs->FitParams().SetInitValue(conv::kPhypMixing13, Rnd::Instance()->RndFitter().Uniform(0.015,0.035));
      break;
    case(conv::kNuMixUndefined):
      cerr << "Mixing parameter type undefinded. Exiting..." << endl;
      exit(1);
      break;
    }
    fFitInputs->FitParams().SetInitValue(conv::kPhypDm23sq, Rnd::Instance()->RndFitter().Uniform(0.0020,0.0030));
    this->RunOnce(false);
    ntries++;
  }

  if(fMigradMinStatusCode == 4) {
    cout << "Fits still failed; trying SIMPLEX" << endl;
    this->RunOnce(true);
  }  

  // stop clock
  w.Stop();
  fFitTime = w.CpuTime();
#ifdef __show_fit_cpu_time__
  cout << "FitterMinuit/DoTheFit >> Fit performed in " 
       << fFitTime << " sec." << endl;
#endif

  // copy fit time and bestfit chisq and status codes 
  // to persistable FitParamList obj

  fFitInputs->FitParams().SetFitTime ( fFitTime );
  fFitInputs->FitParams().SetBestFitChisq ( fFitInputs->BestFitChisquare());

  fFitInputs->FitParams().SetStatusCode("migrad_min", fMigradMinStatusCode);
  fFitInputs->FitParams().SetStatusCode("migrad_err", fMigradErrStatusCode);

  if(fCalcHesseErrors) {
    fFitInputs->FitParams().SetStatusCode("hesse_err", fHesseErrStatusCode);
  }
  if(fCalcMinosErrors) {
    fFitInputs->FitParams().SetStatusCode("minos_err", fMinosErrStatusCode);
  }

  return ok;
}
//____________________________________________________________________________ 
bool FitterMinuit::RunOnce(bool use_simplex) 
{
  fHesseErrStatusCode = -1;

  // Get the fitter
#ifdef __show_fit_steps__
  cout << "FitterMinuit/RunOnce >> Getting the fitter..." << endl;
#endif
  TVirtualFitter * fitter = TVirtualFitter::GetFitter();

  Double_t arg[5];
  arg[0] = fStr;
  fitter->ExecuteCommand("SET STRategy",arg,1);    

  // Pass fit inputs (data, pdf) to MINUIT
#ifdef __show_fit_steps__
  cout << "FitterMinuit/RunOnce >> Passing fit inputs (data & pdf)..." << endl;
#endif
  gMinuit->SetObjectFit(fFitInputs);

  // Run MINUIT's minimization step
#ifdef __show_fit_steps__
  cout << "FitterMinuit/RunOnce >> Running MIGRAD..." << endl;
#endif
  Double_t arglist[100];
  arglist[0] = 50000; // num of func calls
  arglist[1] = fTol;  // tolerance
  if(!use_simplex)
    fMigradMinStatusCode = fitter->ExecuteCommand("MIGRAD",arglist, 2);
  //SIMPLEX
  else {
    arglist[1] = 0.1; // higher tolerance
    fMigradMinStatusCode = fitter->ExecuteCommand("SIMPLEX",arglist, 2);
  }
  
  //return on not converged
  if(fMigradMinStatusCode == 4 && !use_simplex)
    return false;
  else if(fMigradMinStatusCode == 4 && use_simplex) {
    //copy current to best fit, to stop horrible crashes
    fFitInputs->FitParams().CopyCurrToBestFitValues();
    return false;
  }

  // Print fit result
#ifdef __show_fit_result__
#ifdef __show_fit_steps__
  cout << "FitterMinuit/RunOnce >> Printing MIGRAD fit results..." << endl;
#endif
  Double_t aminn = 0;
  fitter->PrintResults(3,aminn);
#endif

  // Found best fit! Copy current --> best-fit values
  fFitInputs->FitParams().CopyCurrToBestFitValues();

  // Get the error matrix status
  double ha(0.0);     // dummy vars
  double edm, errdef; // dummy vars
  int nvpar, nparx;   // dummy vars
  int errMatrixStatus = fitter->GetStats(ha, edm, errdef, nvpar, nparx);
  fMigradErrStatusCode = errMatrixStatus;

  // Calculate HESSE or MINOS errors
  if(fCalcHesseErrors) {
    fHesseErrStatusCode = fitter->ExecuteCommand("HESSE",arglist, 1);
  }
  else if(fCalcMinosErrors) {
    fMinosErrStatusCode = fitter->ExecuteCommand("MINOS",arglist, 1);
  }

  // Get error matrix status
#ifdef __show_fit_steps__
  cout << "FitterMinuit/RunOnce >> Getting error matrix status..." << endl;
#endif
  errMatrixStatus = fitter->GetStats(ha, edm, errdef, nvpar, nparx);
  fMigradErrStatusCode = errMatrixStatus;

  // Get errors
#ifdef __show_fit_steps__
  cout << "FitterMinuit/RunOnce >> Getting errors..." << endl;
#endif
  const int nparamsmax = knppmax;
  double errplus [nparamsmax] = {0.};
  double errminus[nparamsmax] = {0.};
  double errpara [nparamsmax] = {0.};
  for(int ip = 0; ip < nparamsmax; ip++) {
    double glob_corr; // dummy vars
    fitter->GetErrors(ip,errplus[ip],errminus[ip],errpara[ip],glob_corr);
  }
  fFitInputs->FitParams().SetErrors(errplus, errminus, errpara);

  // Print fit result
  if(fCalcHesseErrors) {
#ifdef __show_fit_result__
#ifdef __show_fit_steps__
    cout << "FitterMinuit/RunOnce >> Printing HESSE fit results..." << endl;
#endif
    Double_t amin = 0;
    fitter->PrintResults(3,amin);
  //fitter->PrintResults(1,amin); // prints limits
  }
  cout << "FitterMinuit/RunOnce >> Final error matrix status = " 
       << fMigradErrStatusCode << endl; 
#endif

#ifdef __show_fit_steps__
  cout << "FitterMinuit/RunOnce >> Done" << endl;
#endif

  for(int i = 0; i <= 47; i++)
    cout << i << "\t" << fitter->GetParameter(i) << endl;

  return (fMigradErrStatusCode == 3);
}
//____________________________________________________________________________ 
