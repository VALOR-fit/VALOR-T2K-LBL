//____________________________________________________________________________ 
/*

 \copyright (c) 2010-2013, VaLOR group.

*/ 
//____________________________________________________________________________

#define __show_fit_init_info__ 
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
  fCalcEScaleErrors    = false;
//fSimplex             = false;
  fTol                 = 0.01;
  fStr                 = 1;
  fMigradMinStatusCode = -999;
  fMigradErrStatusCode = -999;      
  fHesseErrStatusCode  = -999; 
  fMinosErrStatusCode  = -999; 
  fFitTime             = 0;
  fFitTries            = 0;
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
  fCalcHesseErrors  = fOption.Contains("hesse");
  fCalcMinosErrors  = fOption.Contains("minos");
  fCalcEScaleErrors = (fCalcMinosErrors) ? false : fOption.Contains("escale"); //don't need to calculate minos errors twice for escale

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
    if(s == kSKsRHCFCFV1RingMuLike) fNFitBins += kRHCNErecoBins1Rmu;
    if(s == kSKsRHCFCFV1RingELike ) fNFitBins += kRHCNErecoBins1Re;
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
  fFitTries = -1;

  // reset fit time, bestfit chisq and status codes in
  // persistable FitParamList obj
  fFitInputs->FitParams().SetFitTime(0);
  fFitInputs->FitParams().SetBestFitChisq(0);
  fFitInputs->FitParams().ClearStatusCodes();

  // Run the fitter once - we have the option here to check the status code
  // and if the fitter failed, re-seed the fitter and try again
  bool ok = this->RunOnce(false);

  // 4: migrad not converged
  int maxntries = 10;
  fFitTries = 1;
  OscAnaType_t ana = Opt::Instance()->gOscAnaTypeChoice;
  if(ana != kAnaNuMuDisappOnly)
    maxntries = 25;
  while(fMigradMinStatusCode == 4 && fFitTries < maxntries) {
    cout << "Original fit failed; reseeding & trying again. Attempt: " << fFitTries << endl;
    //randomise oscillation parameters
    NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
    switch(omp) {
    case(conv::kNuMixSinsq2Theta):
      if((ana == kAnaNuMuDisappOnly || ana == kAnaCombined) && fFitInputs->FitParams().Exists(conv::kPhypMixing23) && !fFitInputs->FitParams().IsFixed(conv::kPhypMixing23))
	fFitInputs->FitParams().SetInitValue(conv::kPhypMixing23, Rnd::Instance()->RndFitter().Uniform(0.80,0.98));
      if((ana == kAnaNuEAppOnly     || ana == kAnaCombined) && fFitInputs->FitParams().Exists(conv::kPhypMixing13) && !fFitInputs->FitParams().IsFixed(conv::kPhypMixing13))
	fFitInputs->FitParams().SetInitValue(conv::kPhypMixing13, Rnd::Instance()->RndFitter().Uniform(0.05,0.15));
      break;
    case(conv::kNuMixSinsqTheta):
      if((ana == kAnaNuMuDisappOnly || ana == kAnaCombined) && fFitInputs->FitParams().Exists(conv::kPhypMixing23) && !fFitInputs->FitParams().IsFixed(conv::kPhypMixing23)) {
	if (fFitInputs->FitParams().InitValue(conv::kPhypMixing23)< 0.515)
	  fFitInputs->FitParams().SetInitValue(conv::kPhypMixing23, Rnd::Instance()->RndFitter().Uniform(0.35,0.50));
	else
	  fFitInputs->FitParams().SetInitValue(conv::kPhypMixing23, Rnd::Instance()->RndFitter().Uniform(0.53,0.68));
      }
      if((ana == kAnaNuEAppOnly || ana == kAnaCombined) && fFitInputs->FitParams().Exists(conv::kPhypMixing13) && !fFitInputs->FitParams().IsFixed(conv::kPhypMixing13))
	fFitInputs->FitParams().SetInitValue(conv::kPhypMixing13, Rnd::Instance()->RndFitter().Uniform(0.015,0.035));
      break;
    case(conv::kNuMixUndefined):
      cerr << "Mixing parameter type undefinded. Exiting..." << endl;
      exit(1);
      break;
    }
    if(fFitInputs->FitParams().Exists(conv::kPhypDm23sq) && !fFitInputs->FitParams().IsFixed(conv::kPhypDm23sq))
      fFitInputs->FitParams().SetInitValue(conv::kPhypDm23sq, Rnd::Instance()->RndFitter().Uniform(0.0020,0.0030));
    if(fFitInputs->FitParams().Exists(conv::kPhypDm12sq) && !fFitInputs->FitParams().IsFixed(conv::kPhypDm12sq))
      fFitInputs->FitParams().SetInitValue(conv::kPhypDm12sq, util::RandomizeDm12sq());
    if(fFitInputs->FitParams().Exists(conv::kPhypMixing12) && !fFitInputs->FitParams().IsFixed(conv::kPhypMixing12))
      fFitInputs->FitParams().SetInitValue(conv::kPhypMixing12, util::RandomizeMixing12());
    if(fFitInputs->FitParams().Exists(conv::kPhypDeltaCP) && !fFitInputs->FitParams().IsFixed(conv::kPhypDeltaCP))
      fFitInputs->FitParams().SetInitValue(conv::kPhypDeltaCP, util::RandomizeDeltaCP());
    //randomise systematics
    //first check if there is a floating systematic in the fit (only need to look at one, since we're only ever fitting either 0 or all systematics)
    if(fFitInputs->FitParams().Exists(conv::kNuisXSecNormNuEToNuMu) && !fFitInputs->FitParams().IsFixed(conv::kNuisXSecNormNuEToNuMu)) {
      //tweak everything
      fFitInputs->Expectation().TweakAllSystematics();
      //get the tweaks
      PhysParamList tweaks = fFitInputs->Expectation().Tweaks();
      //initialise the FitParamList with these random tweaks
      fFitInputs->FitParams().FromPhysParamList(tweaks,"initial",true);
      fFitInputs->FitParams().FromPhysParamList(tweaks,"current",true);
    }//systematics exist in fit
    //cout << fFitInputs->FitParams() << endl;
    this->RunOnce(false);
    fFitTries++;
  }

  if(fMigradMinStatusCode == 4) {
    cout << "Fits still failed; trying SIMPLEX" << endl;
    this->RunOnce(true);
    fFitTries = 99;
  }  

  // stop clock
  w.Stop();
  fFitTime = w.CpuTime();
#ifdef __show_fit_cpu_time__
  cout << "FitterMinuit/DoTheFit >> Fit performed in " 
       << fFitTime << " sec." << endl;
  cout << "FitterMinuit/DoTheFit >> Fit performed in " 
       << fFitTries << " attempts" << endl;
#endif

  // copy fit time and bestfit chisq and status codes 
  // to persistable FitParamList obj

  fFitInputs->FitParams().SetFitTime ( fFitTime  );
  fFitInputs->FitParams().SetFitTries( fFitTries );
  fFitInputs->FitParams().SetBestFitChisq ( fFitInputs->BestFitChisquare());

  fFitInputs->FitParams().SetStatusCode("migrad_min", fMigradMinStatusCode);
  fFitInputs->FitParams().SetStatusCode("migrad_err", fMigradErrStatusCode);

  if(fCalcHesseErrors) {
    fFitInputs->FitParams().SetStatusCode("hesse_err", fHesseErrStatusCode);
  }
  if(fCalcMinosErrors) {
    fFitInputs->FitParams().SetStatusCode("minos_err", fMinosErrStatusCode);
  }
  else if(fCalcEScaleErrors) {
    fFitInputs->FitParams().SetStatusCode("minos_err", fMinosErrStatusCode);
  }

  return ok;
}
//____________________________________________________________________________ 
bool FitterMinuit::RunOnce(bool use_simplex) 
{

  fHesseErrStatusCode = -1;
  double fcntemp;

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
  if(!use_simplex) {
    /*
    //in disappearance mode, when fitting for deltacp, do two fits - fit everything but deltacp, then fit only deltacp
    OscAnaType_t ana = Opt::Instance()->gOscAnaTypeChoice;
    if(ana == kAnaNuMuDisappOnly && fFitInputs->FitParams().GetParListName() == "nuosc3f_numudisapp_nuis4+nuisance") {
#ifdef __show_fit_steps__
      cout << "FitterMinuit/RunOnce >> Performing two step fit. Fit everything but dcp, fit only dcp..." << endl;
#endif
      //fix dcp
      fFitInputs->FitParams().FixParameter(conv::kPhypDeltaCP);
      double dcp_slot = (double)(fFitInputs->FitParams().Id(conv::kPhypDeltaCP) + 1); //FitParamList Id runs from 0. MINUIT Id runs from 1.
      fitter->ExecuteCommand("FIX",&dcp_slot,1);
      //fit
      fMigradMinStatusCode = fitter->ExecuteCommand("MIGRAD",arglist, 2);
      //return if the not-dcp fit failed
      if(fMigradMinStatusCode == 4) {
	fFitInputs->FitParams().FreeParameter(conv::kPhypDeltaCP);
	fitter->ExecuteCommand("RELEASE",&dcp_slot,1);
	cout << "FitterMinuit/RunOnce >> Fit did not converge when fixing dcp" << endl;
	return false;
      }
      //print
#ifdef __show_fit_result__
#ifdef __show_fit_steps__
      cout << "FitterMinuit/RunOnce >> Printing MIGRAD intermediate fit results..." << endl;
#endif
      fitter->ExecuteCommand("SHOW FCNVALUE",&fcntemp,0); //print the bestfit chisq
      Double_t aminn = 0;
      fitter->PrintResults(3,aminn);
#endif
      //fix everything else
      double fixed_slot[100] = {0};
      int nfix = 0;
      for(int ip = 0; ip < fFitInputs->FitParams().Size(); ip++) {
	PhysParamId_t p = fFitInputs->FitParams().WhichParam(ip);
	if(!fFitInputs->FitParams().IsFixed(p)) {
	  fFitInputs->FitParams().FixParameter(ip);
 	  fixed_slot[nfix] = ip+1; //FitParamList Id runs from 0. MINUIT Id runs from 1.
	  nfix++;
	}
      }//ip
      //free dcp
      fFitInputs->FitParams().FreeParameter(conv::kPhypDeltaCP);
      fitter->ExecuteCommand("RELEASE",&dcp_slot,1);
      //fit
      fMigradMinStatusCode = fitter->ExecuteCommand("MIGRAD",arglist, 2);
      //free everything that I just fixed
      for(int ip = 0; ip <= nfix; ip++) {
	fFitInputs->FitParams().FreeParameter(ip);
      }
      fitter->ExecuteCommand("RELEASE",fixed_slot,nfix);
    }//disapp & free dcp
    else
    */
      fMigradMinStatusCode = fitter->ExecuteCommand("MIGRAD",arglist, 2);
  }//!_use_simples
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

  // Found best fit! Copy current --> best-fit values
  for(int ip = 0; ip < fFitInputs->FitParams().Size(); ip++) {
    double value, error;
    gMinuit->GetParameter(ip, value, error);
    PhysParamId_t p = fFitInputs->FitParams().WhichParam(ip);
    fFitInputs->FitParams().Update(p, value);
  }
  fFitInputs->FitParams().CopyCurrToBestFitValues();

  // Print fit result
#ifdef __show_fit_result__
#ifdef __show_fit_steps__
  cout << "FitterMinuit/RunOnce >> Printing MIGRAD fit results..." << endl;
#endif
  fitter->ExecuteCommand("SHOW FCNVALUE",&fcntemp,0);
  Double_t aminn = 0;
  fitter->PrintResults(3,aminn);
#endif

  // Get the error matrix status
  double ha(0.0);     // dummy vars
  double edm, errdef; // dummy vars
  int nvpar, nparx;   // dummy vars
  int errMatrixStatus = fitter->GetStats(ha, edm, errdef, nvpar, nparx);
  fMigradErrStatusCode = errMatrixStatus;

  // Calculate HESSE or MINOS errors
  if(fCalcHesseErrors) {
#ifdef __show_fit_steps__
    cout << "FitterMinuit/RunOnce >> Running HESSE..." << endl;
#endif
    fHesseErrStatusCode = fitter->ExecuteCommand("HESSE",arglist, 1);
  }
  //else 
  if(fCalcMinosErrors) {
#ifdef __show_fit_steps__
    cout << "FitterMinuit/RunOnce >> Running MINOS..." << endl;
#endif
    fMinosErrStatusCode = fitter->ExecuteCommand("MINOS",arglist, 1);
  }
  //only run this if minos has not already been run
  else if(fCalcEScaleErrors) {
    int escale_slot = fFitInputs->FitParams().Id(conv::kNuisSuperkEnergyScale) + 1; //FitParamList Id runs from 0. MINUIT Id runs from 1.
#ifdef __show_fit_steps__
    cout << "FitterMinuit/RunOnce >> Running MINOS for energy scale (MINUIT parameter " << escale_slot << ") only..." << endl;
#endif
    arglist[1] = escale_slot;    //min parno
    arglist[2] = escale_slot;    //max parno
    fMinosErrStatusCode = fitter->ExecuteCommand("MINOS",arglist, 3);
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
    if(fCalcMinosErrors)
      fitter->GetErrors(ip,errplus[ip],errminus[ip],errpara[ip],glob_corr);
    else
      errpara[ip] = fitter->GetParError(ip);
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

  return true;//(fMigradErrStatusCode == 3);
}
//____________________________________________________________________________ 
