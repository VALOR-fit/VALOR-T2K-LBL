//____________________________________________________________________________ 
/*

 >> ChangeLog:

 For changes up to version 2010_v3.0.0 (last implementation of the 2010-era analysis)
 see code in the corresponding SVN branch.

 Changes since 2011_v0.9.0:

 @

*/ 
//____________________________________________________________________________

#include <cassert> 
#include <cstdlib> 
#include <iostream> 
#include <sstream>

#include <TMath.h> 
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

#include "Superk1ReNFitter.h" 
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
Superk1ReNFitter::Superk1ReNFitter() :
FitterKernel() 
{

}
//____________________________________________________________________________ 
Superk1ReNFitter::Superk1ReNFitter(FitInputs * fi, Option_t * opt) :
FitterKernel(fi,opt)
{

}
//____________________________________________________________________________ 
Superk1ReNFitter::~Superk1ReNFitter() 
{

}
//____________________________________________________________________________ 
void Superk1ReNFitter::Analyze(Option_t * option) 
{
  util::SetDefaultStyle();

  TString opt = option;
  cout << "Superk1ReNFitter/Analyze >> Option: " << opt.Data() << endl;

  //
  // Print best-fit results used for generating the plots
  //
  fFitInputs->FitParams().ShowOnlyBestFitValuesInNextPrint();
  cout << "Superk1ReNFitter/Analyze >> "
       << "Using the following best-fit results: " << fFitInputs->FitParams();

  //
  // Declare (potential) outpts
  //
  TH1D * bestfit_tot            = 0;
  TH1D * bestfit_numuccqe       = 0;
  TH1D * bestfit_numucc1pi      = 0;
  TH1D * bestfit_numuccoth      = 0;
  TH1D * bestfit_numucccoh      = 0;
  TH1D * bestfit_numunc1pi0     = 0;
  TH1D * bestfit_numunccoh      = 0;
  TH1D * bestfit_numuncoth      = 0;
  TH1D * bestfit_numubarccqe    = 0;
  TH1D * bestfit_numubarcc1pi   = 0;
  TH1D * bestfit_numubarcccoh   = 0;
  TH1D * bestfit_numubarccoth   = 0;
  TH1D * bestfit_numubarnc1pi0  = 0;
  TH1D * bestfit_numubarnccoh   = 0;
  TH1D * bestfit_numubarncoth   = 0;
  TH1D * bestfit_nueccqe        = 0;
  TH1D * bestfit_nuecc1pi       = 0;
  TH1D * bestfit_nuecccoh       = 0;
  TH1D * bestfit_nueccoth       = 0;
  TH1D * bestfit_nuenc1pi0      = 0;
  TH1D * bestfit_nuenccoh       = 0;
  TH1D * bestfit_nuencoth       = 0;
  TH1D * bestfit_nuebarccqe     = 0;
  TH1D * bestfit_nuebarcc1pi    = 0;
  TH1D * bestfit_nuebarcccoh    = 0;
  TH1D * bestfit_nuebarccoth    = 0;
  TH1D * bestfit_nuebarnc1pi0   = 0;
  TH1D * bestfit_nuebarnccoh    = 0;
  TH1D * bestfit_nuebarncoth    = 0;
  TH1D * bestfit_signueccqe     = 0;
  TH1D * bestfit_signuecc1pi    = 0;
  TH1D * bestfit_signueccoth    = 0;
  TH2D * chisq_surface          = 0;
  TCanvas * c1 = 0;
  TCanvas * c2 = 0;
  TCanvas * c3 = 0;
  TCanvas * c4 = 0;

  //
  // Write the best-fit spectra (total & and its components)
  //
  if(opt.Contains("BestFitSpectra")) {

     // `best-fit' --> `current' values
     fFitInputs->FitParams().CopyBestFitToCurrValues();

     // generate expectation
     fFitInputs->Expectation().Generate(fFitInputs->FitParams());
     conv::SuperkSample_t s = conv::kSKsFCFV1RingELike; 
     int idx = util::ArrIdx::Instance()->SuperkSampleIdx(s);
     if(idx != -1) {

     // copy & format Ereco histograms
     bestfit_tot            = new TH1D(*fFitInputs->Expectation().CurrPredictedErecoSpectrum(s));
     bestfit_numuccqe       = new TH1D(*fFitInputs->Expectation().CurrPredictedErecoSpectrum(s, kSKmNuMuCCQE));
     bestfit_numucc1pi      = new TH1D(*fFitInputs->Expectation().CurrPredictedErecoSpectrum(s, kSKmNuMuCC1pi));
     bestfit_numucccoh      = new TH1D(*fFitInputs->Expectation().CurrPredictedErecoSpectrum(s, kSKmNuMuCCcoh));
     bestfit_numuccoth      = new TH1D(*fFitInputs->Expectation().CurrPredictedErecoSpectrum(s, kSKmNuMuCCoth_e));
     bestfit_numunc1pi0     = new TH1D(*fFitInputs->Expectation().CurrPredictedErecoSpectrum(s, kSKmNuMuNC1pi0));
     bestfit_numunccoh      = new TH1D(*fFitInputs->Expectation().CurrPredictedErecoSpectrum(s, kSKmNuMuNCcoh));
     bestfit_numuncoth      = new TH1D(*fFitInputs->Expectation().CurrPredictedErecoSpectrum(s, kSKmNuMuNCoth_e));
     bestfit_numubarccqe    = new TH1D(*fFitInputs->Expectation().CurrPredictedErecoSpectrum(s, kSKmNuMuBarCCQE));
     bestfit_numubarcc1pi   = new TH1D(*fFitInputs->Expectation().CurrPredictedErecoSpectrum(s, kSKmNuMuBarCC1pi));
     bestfit_numubarcccoh   = new TH1D(*fFitInputs->Expectation().CurrPredictedErecoSpectrum(s, kSKmNuMuBarCCcoh));
     bestfit_numubarccoth   = new TH1D(*fFitInputs->Expectation().CurrPredictedErecoSpectrum(s, kSKmNuMuBarCCoth_e));
     bestfit_numubarnc1pi0  = new TH1D(*fFitInputs->Expectation().CurrPredictedErecoSpectrum(s, kSKmNuMuBarNC1pi0));
     bestfit_numubarnccoh   = new TH1D(*fFitInputs->Expectation().CurrPredictedErecoSpectrum(s, kSKmNuMuBarNCcoh));
     bestfit_numubarncoth   = new TH1D(*fFitInputs->Expectation().CurrPredictedErecoSpectrum(s, kSKmNuMuBarNCoth_e));
     bestfit_nueccqe        = new TH1D(*fFitInputs->Expectation().CurrPredictedErecoSpectrum(s, kSKmNuECCQE));
     bestfit_nuecc1pi       = new TH1D(*fFitInputs->Expectation().CurrPredictedErecoSpectrum(s, kSKmNuECC1pi));
     bestfit_nuecccoh       = new TH1D(*fFitInputs->Expectation().CurrPredictedErecoSpectrum(s, kSKmNuECCcoh));
     bestfit_nueccoth       = new TH1D(*fFitInputs->Expectation().CurrPredictedErecoSpectrum(s, kSKmNuECCoth_e));
     bestfit_nuenc1pi0      = new TH1D(*fFitInputs->Expectation().CurrPredictedErecoSpectrum(s, kSKmNuENC1pi0));
     bestfit_nuenccoh       = new TH1D(*fFitInputs->Expectation().CurrPredictedErecoSpectrum(s, kSKmNuENCcoh));
     bestfit_nuencoth       = new TH1D(*fFitInputs->Expectation().CurrPredictedErecoSpectrum(s, kSKmNuENCoth_e));
     bestfit_nuebarccqe     = new TH1D(*fFitInputs->Expectation().CurrPredictedErecoSpectrum(s, kSKmNuEBarCCQE));
     bestfit_nuebarcc1pi    = new TH1D(*fFitInputs->Expectation().CurrPredictedErecoSpectrum(s, kSKmNuEBarCC1pi));
     bestfit_nuebarcccoh    = new TH1D(*fFitInputs->Expectation().CurrPredictedErecoSpectrum(s, kSKmNuEBarCCcoh));
     bestfit_nuebarccoth    = new TH1D(*fFitInputs->Expectation().CurrPredictedErecoSpectrum(s, kSKmNuEBarCCoth_e));
     bestfit_nuebarnc1pi0   = new TH1D(*fFitInputs->Expectation().CurrPredictedErecoSpectrum(s, kSKmNuEBarNC1pi0));
     bestfit_nuebarnccoh    = new TH1D(*fFitInputs->Expectation().CurrPredictedErecoSpectrum(s, kSKmNuEBarNCcoh));
     bestfit_nuebarncoth    = new TH1D(*fFitInputs->Expectation().CurrPredictedErecoSpectrum(s, kSKmNuEBarNCoth_e));
     bestfit_signueccqe     = new TH1D(*fFitInputs->Expectation().CurrPredictedErecoSpectrum(s, kSKmSigNuECCQE));
     bestfit_signuecc1pi    = new TH1D(*fFitInputs->Expectation().CurrPredictedErecoSpectrum(s, kSKmSigNuECC1pi));
     bestfit_signueccoth    = new TH1D(*fFitInputs->Expectation().CurrPredictedErecoSpectrum(s, kSKmSigNuECCoth));


     bestfit_tot            -> SetDirectory ( 0        );
     bestfit_numuccqe       -> SetDirectory ( 0        );
     bestfit_numucc1pi      -> SetDirectory ( 0        );
     bestfit_numucccoh      -> SetDirectory ( 0        );
     bestfit_numuccoth      -> SetDirectory ( 0        );
     bestfit_numunc1pi0     -> SetDirectory ( 0        );
     bestfit_numunccoh      -> SetDirectory ( 0        );
     bestfit_numuncoth      -> SetDirectory ( 0        );
     bestfit_numubarccqe    -> SetDirectory ( 0        );
     bestfit_numubarcc1pi   -> SetDirectory ( 0        );
     bestfit_numubarcccoh   -> SetDirectory ( 0        );
     bestfit_numubarccoth   -> SetDirectory ( 0        );
     bestfit_numubarnc1pi0  -> SetDirectory ( 0        );
     bestfit_numubarnccoh   -> SetDirectory ( 0        );
     bestfit_numubarncoth   -> SetDirectory ( 0        );
     bestfit_nueccqe        -> SetDirectory ( 0        );
     bestfit_nuecc1pi       -> SetDirectory ( 0        );
     bestfit_nuecccoh       -> SetDirectory ( 0        );
     bestfit_nueccoth       -> SetDirectory ( 0        );
     bestfit_nuenc1pi0      -> SetDirectory ( 0        );
     bestfit_nuenccoh       -> SetDirectory ( 0        );
     bestfit_nuencoth       -> SetDirectory ( 0        );
     bestfit_nuebarccqe     -> SetDirectory ( 0        );
     bestfit_nuebarcc1pi    -> SetDirectory ( 0        );
     bestfit_nuebarcccoh    -> SetDirectory ( 0        );
     bestfit_nuebarccoth    -> SetDirectory ( 0        );
     bestfit_nuebarnc1pi0   -> SetDirectory ( 0        );
     bestfit_nuebarnccoh    -> SetDirectory ( 0        );
     bestfit_nuebarncoth    -> SetDirectory ( 0        );
     bestfit_signueccqe     -> SetDirectory ( 0        );
     bestfit_signuecc1pi    -> SetDirectory ( 0        );
     bestfit_signueccoth    -> SetDirectory ( 0        );

     bestfit_numucc1pi      -> SetFillStyle ( 3544     );
     bestfit_numucccoh      -> SetFillStyle ( 3244     );
     bestfit_numuccoth      -> SetFillStyle ( 3244     );
     bestfit_numunc1pi0     -> SetFillStyle ( 3006     );
     bestfit_numunccoh      -> SetFillStyle ( 3006     );
     bestfit_numuncoth      -> SetFillStyle ( 3006     );
     bestfit_numubarccqe    -> SetFillStyle ( 3001     );
     bestfit_numubarcc1pi   -> SetFillStyle ( 3345     );
     bestfit_numubarcccoh   -> SetFillStyle ( 3354     );
     bestfit_numubarccoth   -> SetFillStyle ( 3354     );
     bestfit_numubarnc1pi0  -> SetFillStyle ( 3006     );
     bestfit_numubarnccoh   -> SetFillStyle ( 3006     );
     bestfit_numubarncoth   -> SetFillStyle ( 3006     );
     bestfit_nueccqe        -> SetFillStyle ( 3001     );
     bestfit_nuecc1pi       -> SetFillStyle ( 3345     );
     bestfit_nuecccoh       -> SetFillStyle ( 3354     );
     bestfit_nueccoth       -> SetFillStyle ( 3354     );
     bestfit_nuenc1pi0      -> SetFillStyle ( 3006     );
     bestfit_nuenccoh       -> SetFillStyle ( 3006     );
     bestfit_nuencoth       -> SetFillStyle ( 3006     );
     bestfit_signueccqe     -> SetFillStyle ( 3001     );
     bestfit_signuecc1pi    -> SetFillStyle ( 3345     );
     bestfit_signueccoth    -> SetFillStyle ( 3354     );

     bestfit_numuccqe       -> SetLineColor ( kBlack   );

     bestfit_numucc1pi      -> SetFillColor ( kBlack   );
     bestfit_numucccoh      -> SetFillColor ( kGreen   );
     bestfit_numuccoth      -> SetFillColor ( kBlue    );
     bestfit_numunc1pi0     -> SetFillColor ( kRed     );
     bestfit_numunccoh      -> SetFillColor ( kRed-1   );
     bestfit_numuncoth      -> SetFillColor ( kRed-2   );
     bestfit_numubarccqe    -> SetFillColor ( kAzure-3 );
     bestfit_numubarcc1pi   -> SetFillColor ( kAzure-4 );
     bestfit_numubarcccoh   -> SetFillColor ( kAzure-6 );
     bestfit_numubarccoth   -> SetFillColor ( 6        );
     bestfit_numubarnc1pi0  -> SetFillColor ( kAzure-7 );
     bestfit_numubarnccoh   -> SetFillColor ( kAzure-9 );
     bestfit_numubarncoth   -> SetFillColor ( 6        );
     bestfit_nueccqe        -> SetFillColor ( kGreen+1 );
     bestfit_nuecc1pi       -> SetFillColor ( kGreen+2 );
     bestfit_nuecccoh       -> SetFillColor ( kGreen+2 );
     bestfit_nueccoth       -> SetFillColor ( kBlack   );
     bestfit_nuenc1pi0      -> SetFillColor ( kGreen+5 );
     bestfit_nuenccoh       -> SetFillColor ( kGreen+7 );
     bestfit_nuencoth       -> SetFillColor ( kBlack   );
     bestfit_signueccqe     -> SetFillColor ( kGreen+8 );
     bestfit_signuecc1pi    -> SetFillColor ( kGreen+10);
     bestfit_signueccoth    -> SetFillColor ( kRed     );

     bestfit_tot            -> SetLineWidth ( 2 );
     bestfit_numuccqe       -> SetLineWidth ( 1 );
     bestfit_numucc1pi      -> SetLineWidth ( 1 );
     bestfit_numucccoh      -> SetLineWidth ( 1 );
     bestfit_numuccoth      -> SetLineWidth ( 1 );
     bestfit_numunc1pi0     -> SetLineWidth ( 1 );
     bestfit_numunccoh      -> SetLineWidth ( 1 );
     bestfit_numuncoth      -> SetLineWidth ( 1 );
     bestfit_numubarccqe    -> SetLineWidth ( 1 );
     bestfit_numubarcc1pi   -> SetLineWidth ( 1 );
     bestfit_numubarcccoh   -> SetLineWidth ( 1 );
     bestfit_numubarccoth   -> SetLineWidth ( 1 );
     bestfit_numubarnc1pi0  -> SetLineWidth ( 1 );
     bestfit_numubarnccoh   -> SetLineWidth ( 1 );
     bestfit_numubarncoth   -> SetLineWidth ( 1 );
     bestfit_nueccqe        -> SetLineWidth ( 1 );
     bestfit_nuecc1pi       -> SetLineWidth ( 1 );
     bestfit_nuecccoh       -> SetLineWidth ( 1 );
     bestfit_nueccoth       -> SetLineWidth ( 1 );
     bestfit_nuenc1pi0      -> SetLineWidth ( 1 );
     bestfit_nuenccoh       -> SetLineWidth ( 1 );
     bestfit_nuencoth       -> SetLineWidth ( 1 );
     bestfit_signueccqe     -> SetLineWidth ( 1 );
     bestfit_signuecc1pi    -> SetLineWidth ( 1 );
     bestfit_signueccoth    -> SetLineWidth ( 1 );

     bestfit_tot            -> SetName("bestfit_1Re_tot");
     bestfit_numuccqe       -> SetName("bestfit_1Re_numuccqe");
     bestfit_numucc1pi      -> SetName("bestfit_1Re_numucc1pi");
     bestfit_numucccoh      -> SetName("bestfit_1Re_numucccoh");
     bestfit_numuccoth      -> SetName("bestfit_1Re_numuccoth");
     bestfit_numunc1pi0     -> SetName("bestfit_1Re_numunc1pi0");
     bestfit_numunccoh      -> SetName("bestfit_1Re_numunccoh");
     bestfit_numuncoth      -> SetName("bestfit_1Re_numuncoth");
     bestfit_numubarccqe    -> SetName("bestfit_1Re_numubarcc");
     bestfit_numubarcc1pi   -> SetName("bestfit_1Re_numubarcc1pi");
     bestfit_numubarcccoh   -> SetName("bestfit_1Re_numubarcccoh");
     bestfit_numubarccoth   -> SetName("bestfit_1Re_numubarccoth");
     bestfit_numubarnc1pi0  -> SetName("bestfit_1Re_numubarnc1pi0");
     bestfit_numubarnccoh   -> SetName("bestfit_1Re_numubarnccoh");
     bestfit_numubarncoth   -> SetName("bestfit_1Re_numubarncoth");
     bestfit_nueccqe        -> SetName("bestfit_1Re_nueccqe");
     bestfit_nuecc1pi       -> SetName("bestfit_1Re_nuecc1pi");
     bestfit_nuecccoh       -> SetName("bestfit_1Re_nuecccoh");
     bestfit_nueccoth       -> SetName("bestfit_1Re_nueccoth");
     bestfit_nuenc1pi0      -> SetName("bestfit_1Re_nuenc1pi0");
     bestfit_nuenccoh       -> SetName("bestfit_1Re_nuenccoh");
     bestfit_nuencoth       -> SetName("bestfit_1Re_nuencoth");
     bestfit_signueccqe     -> SetName("bestfit_1Re_signueccqe");
     bestfit_signuecc1pi    -> SetName("bestfit_1Re_signuecc1pi");
     bestfit_signueccoth    -> SetName("bestfit_1Re_signueccoth");

     fResults -> Add (bestfit_tot);
     fResults -> Add (bestfit_numuccqe);
     fResults -> Add (bestfit_numucc1pi);
     fResults -> Add (bestfit_numucccoh);
     fResults -> Add (bestfit_numuccoth);
     fResults -> Add (bestfit_numunc1pi0);
     fResults -> Add (bestfit_numunccoh);
     fResults -> Add (bestfit_numuncoth);
     fResults -> Add (bestfit_numubarccqe);
     fResults -> Add (bestfit_numubarcc1pi);
     fResults -> Add (bestfit_numubarcccoh);
     fResults -> Add (bestfit_numubarccoth);
     fResults -> Add (bestfit_numubarnc1pi0);
     fResults -> Add (bestfit_numubarnccoh);
     fResults -> Add (bestfit_numubarncoth);
     fResults -> Add (bestfit_nueccqe);
     fResults -> Add (bestfit_nuecc1pi);
     fResults -> Add (bestfit_nuecccoh);
     fResults -> Add (bestfit_nueccoth);
     fResults -> Add (bestfit_nuenc1pi0);
     fResults -> Add (bestfit_nuenccoh);
     fResults -> Add (bestfit_nuencoth);
     fResults -> Add (bestfit_signueccqe);
     fResults -> Add (bestfit_signuecc1pi);
     fResults -> Add (bestfit_signueccoth);
  }

  }
  if(opt.Contains("ChisqSurface")) {
    cout << "ChisqSurface" << endl;

    PhysHypothesis_t phyp = fFitInputs->FitParams().PhysicsHypothesis();
     if(phyp == kHypOsc3F) {   
       chisq_surface = this->BuildChisqSurfaceNuEApp2ParamFit();
     }else 
       {
	 cerr << "Superk1ReNFitter/Analyze >> Can't handle hypothesis: " 
	      << util::AsString(phyp) << endl;
	 exit(1);
       }
     
     chisq_surface -> SetName("chisq_surface_1Re");
     fResults -> Add(chisq_surface);
  }

  //
  // Save canvases with some standard plots
  //
  if(opt.Contains("Canvases")) {
    //    if(opt.Contains("BestFitSpectra") && opt.Contains("Observation")) {
    if(opt.Contains("BestFitSpectra")) {
    cout << "Canvases" << endl;

      THStack * contrib_stack = new THStack("stack","spectrum contributions");
      contrib_stack->Add(bestfit_signuecc1pi);
      contrib_stack->Add(bestfit_signueccoth);
      contrib_stack->Add(bestfit_signueccqe);
      contrib_stack->Add(bestfit_nuenc1pi0);
      contrib_stack->Add(bestfit_nuenccoh);
      contrib_stack->Add(bestfit_nuencoth);
      contrib_stack->Add(bestfit_nuecc1pi);
      contrib_stack->Add(bestfit_nueccoth);
      contrib_stack->Add(bestfit_nuecccoh);
      contrib_stack->Add(bestfit_nueccqe);
      contrib_stack->Add(bestfit_numubarnc1pi0);
      contrib_stack->Add(bestfit_numubarnccoh);
      contrib_stack->Add(bestfit_numubarncoth);
      contrib_stack->Add(bestfit_numubarcc1pi);
      contrib_stack->Add(bestfit_numubarcccoh);
      contrib_stack->Add(bestfit_numubarccoth);
      contrib_stack->Add(bestfit_numubarccqe);
      contrib_stack->Add(bestfit_numunc1pi0);
      contrib_stack->Add(bestfit_numunccoh);
      contrib_stack->Add(bestfit_numuncoth);
      contrib_stack->Add(bestfit_numuccoth);
      contrib_stack->Add(bestfit_numucccoh);
      contrib_stack->Add(bestfit_numucc1pi);
      contrib_stack->Add(bestfit_numuccqe);

      //added another stack with just signal
      THStack * contrib_stack_sig = new THStack("stack_sig","spectrum signal contributions");
      contrib_stack_sig->Add(bestfit_signuecc1pi);
      contrib_stack_sig->Add(bestfit_signueccoth);
      contrib_stack_sig->Add(bestfit_signueccqe);

      //and last one with background
      THStack * contrib_stack_back = new THStack("stack_back","spectrum background contributions");
      contrib_stack_back->Add(bestfit_nuenc1pi0);
      contrib_stack_back->Add(bestfit_nuenccoh);
      contrib_stack_back->Add(bestfit_nuencoth);
      contrib_stack_back->Add(bestfit_nuecc1pi);
      contrib_stack_back->Add(bestfit_nueccoth);
      contrib_stack_back->Add(bestfit_nuecccoh);
      contrib_stack_back->Add(bestfit_nueccqe);
      contrib_stack_back->Add(bestfit_numubarnc1pi0);
      contrib_stack_back->Add(bestfit_numubarnccoh);
      contrib_stack_back->Add(bestfit_numubarncoth);
      contrib_stack_back->Add(bestfit_numubarcc1pi);
      contrib_stack_back->Add(bestfit_numubarcccoh);
      contrib_stack_back->Add(bestfit_numubarccoth);
      contrib_stack_back->Add(bestfit_numubarccqe);
      contrib_stack_back->Add(bestfit_numunc1pi0);
      contrib_stack_back->Add(bestfit_numunccoh);
      contrib_stack_back->Add(bestfit_numuncoth);
      contrib_stack_back->Add(bestfit_numuccoth);
      contrib_stack_back->Add(bestfit_numucccoh);
      contrib_stack_back->Add(bestfit_numucc1pi);
      contrib_stack_back->Add(bestfit_numuccqe);

	
      c1 = new TCanvas("c1","Best-fit spectrum (total)",20,20,500,500);
      c1->SetBorderMode(0);
      c1->SetFillColor (0);
      bestfit_tot -> Draw("perr");
      c1->Update();

      c2 = new TCanvas("c2","Best-fit spectrum (total & components)",20,20,500,500);
      c2->SetBorderMode(0);
      c2->SetFillColor (0);
      bestfit_tot   -> Draw("perr");
      contrib_stack -> Draw("same");
      TLegend * spectrum_legend = new TLegend(0.6,0.5,0.9,0.9);
      spectrum_legend->SetFillColor (0);
      spectrum_legend->SetBorderSize(0);
      spectrum_legend->AddEntry(bestfit_tot,             "best fit",                                                    "L" );
      spectrum_legend->AddEntry(bestfit_numuccqe,        "best fit (#nu_{#mu} CCQE)",                                   "L" );
      spectrum_legend->AddEntry(bestfit_numucc1pi,       "best fit (#nu_{#mu} CC1#pi)",                                 "LF");
      spectrum_legend->AddEntry(bestfit_numucccoh,       "best fit (#nu_{#mu} CC coherence)",                           "LF");
      spectrum_legend->AddEntry(bestfit_numuccoth,       "best fit (#nu_{#mu} CC other)",                               "LF");
      spectrum_legend->AddEntry(bestfit_numunc1pi0,      "best fit (#nu_{#mu}/#nu_{#tau} NC 1pi0)",                     "LF");
      spectrum_legend->AddEntry(bestfit_numunccoh,       "best fit (#nu_{#mu}/#nu_{#tau} NC coherence)",                "LF");
      spectrum_legend->AddEntry(bestfit_numuncoth,       "best fit (#nu_{#mu}/#nu_{#tau} NC other)",                    "LF");
      spectrum_legend->AddEntry(bestfit_numubarccqe,     "best fit (#bar{#nu_{#mu}} CCQE)",                             "LF");
      spectrum_legend->AddEntry(bestfit_numubarcc1pi,    "best fit (#bar{#nu_{#mu}} CC1#pi)",                           "LF");
      spectrum_legend->AddEntry(bestfit_numubarcccoh,    "best fit (#bar{#nu_{#mu}} CC coherence)",                     "LF");
      spectrum_legend->AddEntry(bestfit_numubarccoth,    "best fit (#bar{#nu_{#mu}} CC other)",                         "LF");
      spectrum_legend->AddEntry(bestfit_numubarnc1pi0,   "best fit (#bar{#nu_{#mu}}/#bar{#nu_{#tau}} NC 1pi0)",         "LF");
      spectrum_legend->AddEntry(bestfit_numubarnccoh,    "best fit (#bar{#nu_{#mu}}/#bar{#nu_{#tau}} NC coherence)",    "LF");
      spectrum_legend->AddEntry(bestfit_numubarncoth,    "best fit (#bar{#nu_{#mu}}/#bar{#nu_{#tau}} NC other)",        "LF");
      spectrum_legend->AddEntry(bestfit_nueccqe,         "best fit (#nu_{e} CCQE)",                                     "LF");
      spectrum_legend->AddEntry(bestfit_nuecc1pi,        "best fit (#nu_{e} CC1#pi)",                                   "LF");
      spectrum_legend->AddEntry(bestfit_nuecccoh,        "best fit (#nu_{e} CC coherence)",                             "LF");
      spectrum_legend->AddEntry(bestfit_nueccoth,        "best fit (#nu_{e} CC other)",                                 "LF");
      spectrum_legend->AddEntry(bestfit_nuenc1pi0,       "best fit (#nu_{e} NC 1pi0)",                                  "LF");
      spectrum_legend->AddEntry(bestfit_nuenccoh,        "best fit (#nu_{e} NC coherence)",                             "LF");
      spectrum_legend->AddEntry(bestfit_nuencoth,        "best fit (#nu_{e} NC other)",                                 "LF");
      spectrum_legend->AddEntry(bestfit_signueccqe,      "best fit (signal-#nu_{e} CCQE)",                              "LF");
      spectrum_legend->AddEntry(bestfit_signuecc1pi,     "best fit (signal-#nu_{e} CC1#pi)",                            "LF");
      spectrum_legend->AddEntry(bestfit_signueccoth,     "best fit (signal-#nu_{e} CC other)",                          "LF");
      spectrum_legend->Draw();
      c2->Update();

      //just signal
      c3 = new TCanvas("c3","Best-fit spectrum (total & signal components)",20,20,500,500);
      c3->SetBorderMode(0);
      c3->SetFillColor (0);
      bestfit_tot       -> Draw("perr");
      contrib_stack_sig -> Draw("same");
      TLegend * spectrum_legend_sig = new TLegend(0.6,0.5,0.9,0.9);
      spectrum_legend_sig->SetFillColor (0);
      spectrum_legend_sig->SetBorderSize(0);
      spectrum_legend_sig->AddEntry(bestfit_tot,             "best fit",                                                    "L" );
      spectrum_legend_sig->AddEntry(bestfit_signueccqe,      "best fit (signal-#nu_{e} CCQE)",                              "LF");
      spectrum_legend_sig->AddEntry(bestfit_signuecc1pi,     "best fit (signal-#nu_{e} CC1#pi)",                            "LF");
      spectrum_legend_sig->AddEntry(bestfit_signueccoth,     "best fit (signal-#nu_{e} CC other)",                          "LF");
      spectrum_legend_sig->Draw();
      c3->Update();

      //background
      c4 = new TCanvas("c4","Best-fit spectrum (total & background components)",20,20,500,500);
      c4->SetBorderMode(0);
      c4->SetFillColor (0);
      bestfit_tot   -> Draw("perr");
      contrib_stack_back -> Draw("same");
      TLegend * spectrum_legend_back = new TLegend(0.6,0.5,0.9,0.9);
      spectrum_legend_back->SetFillColor (0);
      spectrum_legend_back->SetBorderSize(0);
      spectrum_legend_back->AddEntry(bestfit_tot,             "best fit",                                                    "L" );
      spectrum_legend_back->AddEntry(bestfit_numuccqe,        "best fit (#nu_{#mu} CCQE)",                                   "L" );
      spectrum_legend_back->AddEntry(bestfit_numucc1pi,       "best fit (#nu_{#mu} CC1#pi)",                                 "LF");
      spectrum_legend_back->AddEntry(bestfit_numucccoh,       "best fit (#nu_{#mu} CC coherence)",                           "LF");
      spectrum_legend_back->AddEntry(bestfit_numuccoth,       "best fit (#nu_{#mu} CC other)",                               "LF");
      spectrum_legend_back->AddEntry(bestfit_numunc1pi0,      "best fit (#nu_{#mu}/#nu_{#tau} NC 1pi0)",                     "LF");
      spectrum_legend_back->AddEntry(bestfit_numunccoh,       "best fit (#nu_{#mu}/#nu_{#tau} NC coherence)",                "LF");
      spectrum_legend_back->AddEntry(bestfit_numuncoth,       "best fit (#nu_{#mu}/#nu_{#tau} NC other)",                    "LF");
      spectrum_legend_back->AddEntry(bestfit_numubarccqe,     "best fit (#bar{#nu_{#mu}} CCQE)",                             "LF");
      spectrum_legend_back->AddEntry(bestfit_numubarcc1pi,    "best fit (#bar{#nu_{#mu}} CC1#pi)",                           "LF");
      spectrum_legend_back->AddEntry(bestfit_numubarcccoh,    "best fit (#bar{#nu_{#mu}} CC coherence)",                     "LF");
      spectrum_legend_back->AddEntry(bestfit_numubarccoth,    "best fit (#bar{#nu_{#mu}} CC other)",                         "LF");
      spectrum_legend_back->AddEntry(bestfit_numubarnc1pi0,   "best fit (#bar{#nu_{#mu}}/#bar{#nu_{#tau}} NC 1pi0)",         "LF");
      spectrum_legend_back->AddEntry(bestfit_numubarnccoh,    "best fit (#bar{#nu_{#mu}}/#bar{#nu_{#tau}} NC coherence)",    "LF");
      spectrum_legend_back->AddEntry(bestfit_numubarncoth,    "best fit (#bar{#nu_{#mu}}/#bar{#nu_{#tau}} NC other)",        "LF");
      spectrum_legend_back->AddEntry(bestfit_nueccqe,         "best fit (#nu_{e} CCQE)",                                     "LF");
      spectrum_legend_back->AddEntry(bestfit_nuecc1pi,        "best fit (#nu_{e} CC1#pi)",                                   "LF");
      spectrum_legend_back->AddEntry(bestfit_nuecccoh,        "best fit (#nu_{e} CC coherence)",                             "LF");
      spectrum_legend_back->AddEntry(bestfit_nueccoth,        "best fit (#nu_{e} CC other)",                                 "LF");
      spectrum_legend_back->AddEntry(bestfit_nuenc1pi0,       "best fit (#nu_{e} NC 1pi0)",                                  "LF");
      spectrum_legend_back->AddEntry(bestfit_nuenccoh,        "best fit (#nu_{e} NC coherence)",                             "LF");
      spectrum_legend_back->AddEntry(bestfit_nuencoth,        "best fit (#nu_{e} NC other)",                                 "LF");
      spectrum_legend_back->Draw();
      c4->Update();
	
      c1->SetName("canvas_1Re_1");
      c2->SetName("canvas_1Re_2");
      c3->SetName("canvas_1Re_3");
      c4->SetName("canvas_1Re_4");

      fResults->Add(c1);
      fResults->Add(c2);
      fResults->Add(c3);
      fResults->Add(c4);

   } // BestFitSpectra & Observation
  } // Canvases

}

//____________________________________________________________________________ 
void Superk1ReNFitter::Save(const char * filename) 
{
  Save(filename, "recreate");
}

//____________________________________________________________________________ 
void Superk1ReNFitter::Save(const char * filename, const char* option_file) 
{
  util::SetDefaultStyle();

  cout << "Superk1ReNFitter/Save >> Writing fit results in: " 
       << filename << endl;
  TFile * file_out = new TFile(filename, option_file);

  TH1D * bestfit_tot            = (TH1D*) fResults->FindObject("bestfit_1Re_tot");
  TH1D * bestfit_numuccqe       = (TH1D*) fResults->FindObject("bestfit_1Re_numuccqe");
  TH1D * bestfit_numucc1pi      = (TH1D*) fResults->FindObject("bestfit_1Re_numucc1pi");
  TH1D * bestfit_numucccoh      = (TH1D*) fResults->FindObject("bestfit_1Re_numucccoh");
  TH1D * bestfit_numuccoth      = (TH1D*) fResults->FindObject("bestfit_1Re_numuccoth");
  TH1D * bestfit_numunc1pi0     = (TH1D*) fResults->FindObject("bestfit_1Re_numunc1pi0");
  TH1D * bestfit_numunccoh      = (TH1D*) fResults->FindObject("bestfit_1Re_numunccoh");
  TH1D * bestfit_numuncoth      = (TH1D*) fResults->FindObject("bestfit_1Re_numuncoth");
  TH1D * bestfit_numubarccqe    = (TH1D*) fResults->FindObject("bestfit_1Re_numubarccqe");
  TH1D * bestfit_numubarcc1pi   = (TH1D*) fResults->FindObject("bestfit_1Re_numubarcc1pi");
  TH1D * bestfit_numubarcccoh   = (TH1D*) fResults->FindObject("bestfit_1Re_numubarcccoh");
  TH1D * bestfit_numubarccoth   = (TH1D*) fResults->FindObject("bestfit_1Re_numubarccoth");
  TH1D * bestfit_numubarnc1pi0  = (TH1D*) fResults->FindObject("bestfit_1Re_numubarnc1pi0");
  TH1D * bestfit_numubarnccoh   = (TH1D*) fResults->FindObject("bestfit_1Re_numubarnccoh");
  TH1D * bestfit_numubarncoth   = (TH1D*) fResults->FindObject("bestfit_1Re_numubarncoth");
  TH1D * bestfit_nueccqe        = (TH1D*) fResults->FindObject("bestfit_1Re_nueccqe");
  TH1D * bestfit_nuecc1pi       = (TH1D*) fResults->FindObject("bestfit_1Re_nuecc1pi");
  TH1D * bestfit_nuecccoh       = (TH1D*) fResults->FindObject("bestfit_1Re_nuecccoh");
  TH1D * bestfit_nueccoth       = (TH1D*) fResults->FindObject("bestfit_1Re_nueccoth");
  TH1D * bestfit_nuenc1pi0      = (TH1D*) fResults->FindObject("bestfit_1Re_nuenc1pi0");
  TH1D * bestfit_nuenccoh       = (TH1D*) fResults->FindObject("bestfit_1Re_nuenccoh");
  TH1D * bestfit_nuencoth       = (TH1D*) fResults->FindObject("bestfit_1Re_nuencoth");
  TH1D * bestfit_signueccqe     = (TH1D*) fResults->FindObject("bestfit_1Re_signueccqe");
  TH1D * bestfit_signuecc1pi    = (TH1D*) fResults->FindObject("bestfit_1Re_signuecc1pi");
  TH1D * bestfit_signueccoth    = (TH1D*) fResults->FindObject("bestfit_1Re_signueccoth");
  TH2D * chisq_surface          = (TH2D*) fResults->FindObject("chisq_surface_1Re");
  TCanvas * c1 = (TCanvas*) fResults->FindObject("canvas_1Re_1");
  TCanvas * c2 = (TCanvas*) fResults->FindObject("canvas_1Re_2");
  TCanvas * c3 = (TCanvas*) fResults->FindObject("canvas_1Re_3");
  TCanvas * c4 = (TCanvas*) fResults->FindObject("canvas_1Re_4");

  if ( bestfit_tot          ) { bestfit_tot            -> Write(); }
  if ( bestfit_numuccqe     ) { bestfit_numuccqe       -> Write(); }
  if ( bestfit_numucc1pi    ) { bestfit_numucc1pi      -> Write(); }
  if ( bestfit_numucccoh    ) { bestfit_numucccoh      -> Write(); }
  if ( bestfit_numuccoth    ) { bestfit_numuccoth      -> Write(); }
  if ( bestfit_numunc1pi0   ) { bestfit_numunc1pi0     -> Write(); }
  if ( bestfit_numunccoh    ) { bestfit_numunccoh      -> Write(); }
  if ( bestfit_numuncoth    ) { bestfit_numuncoth      -> Write(); }
  if ( bestfit_numubarccqe  ) { bestfit_numubarccqe    -> Write(); }
  if ( bestfit_numubarcc1pi ) { bestfit_numubarcc1pi   -> Write(); }
  if ( bestfit_numubarcccoh ) { bestfit_numubarcccoh   -> Write(); }
  if ( bestfit_numubarccoth ) { bestfit_numubarccoth   -> Write(); }
  if ( bestfit_numubarnc1pi0) { bestfit_numubarnc1pi0  -> Write(); }
  if ( bestfit_numubarnccoh ) { bestfit_numubarnccoh   -> Write(); }
  if ( bestfit_numubarncoth ) { bestfit_numubarncoth   -> Write(); }
  if ( bestfit_nueccqe      ) { bestfit_nueccqe        -> Write(); }
  if ( bestfit_nuecc1pi     ) { bestfit_nuecc1pi       -> Write(); }
  if ( bestfit_nuecccoh     ) { bestfit_nuecccoh       -> Write(); }
  if ( bestfit_nueccoth     ) { bestfit_nueccoth       -> Write(); }
  if ( bestfit_nuenc1pi0    ) { bestfit_nuenc1pi0      -> Write(); }
  if ( bestfit_nuenccoh     ) { bestfit_nuenccoh       -> Write(); }
  if ( bestfit_nuencoth     ) { bestfit_nuencoth       -> Write(); }
  if ( bestfit_signueccqe   ) { bestfit_signueccqe     -> Write(); }
  if ( bestfit_signuecc1pi  ) { bestfit_signuecc1pi    -> Write(); }
  if ( bestfit_signueccoth  ) { bestfit_signueccoth    -> Write(); }

  if ( chisq_surface        ) { chisq_surface          -> Write(); }
  if ( c1                   ) { c1                     -> Write(); }
  if ( c2                   ) { c2                     -> Write(); }
  if ( c3                   ) { c3                     -> Write(); }
  if ( c4                   ) { c4                     -> Write(); }

  file_out->Close();
  delete file_out;
}
//____________________________________________________________________________ 
TH2D * Superk1ReNFitter::BuildChisqSurfaceNuEApp2ParamFit(void) 
{
   FitParamList & fp = fFitInputs->FitParams();
   //lorena -check parameters
   double mixing13_min   = fp.MinValue(conv::kPhypMixing13);
   double mixing13_max   = fp.MaxValue(conv::kPhypMixing13);
   double mixing13_step  = fp.StepSize(conv::kPhypMixing13);
   double deltacp_min         = fp.MinValue(conv::kPhypDeltaCP);
   double deltacp_max         = fp.MaxValue(conv::kPhypDeltaCP);
   double deltacp_step        = fp.StepSize(conv::kPhypDeltaCP);

   // make sure that min and max are at the bin centre
   int mixing13_nbins  = 1 + (int) ((mixing13_max - mixing13_min) / mixing13_step);
   int deltacp_nbins        = 1 + (int) ((deltacp_max -        deltacp_min )       / deltacp_step );

   TH2D * chisq_surface = new TH2D(
				   "chisq_surface","#chi^{2} surface",
				   mixing13_nbins, mixing13_min - 0.5*mixing13_step, mixing13_max + 0.5*mixing13_step,
				   deltacp_nbins, deltacp_min - 0.5*deltacp_step, deltacp_max + 0.5*deltacp_step);
   chisq_surface->SetDirectory(0);
   chisq_surface->GetXaxis()->SetRangeUser(0.,1.);
   
   for(int i = 1; i <= chisq_surface->GetNbinsX(); i++) {
     for(int j = 1; j <= chisq_surface->GetNbinsY(); j++) {
       double mixing13  = chisq_surface->GetXaxis()->GetBinCenter(i);
       double deltacp        = chisq_surface->GetYaxis()->GetBinCenter(j);
       fp.Update( conv::kPhypMixing13, mixing13 );
       fp.Update( conv::kPhypDeltaCP,      deltacp        );
       double chisq = fFitInputs->Chisquare();
       //       cout << "mixing13 = " << mixing13 << "   and deltacp =" << deltacp << "  with chisq=" << chisq << endl;
       chisq_surface->SetBinContent(i,j,chisq);
     }//j
   }//i
   
   chisq_surface->GetXaxis()->SetTitle("sin^{2}(2#theta_{13})");
   chisq_surface->GetYaxis()->SetTitle("#delta CP");
   chisq_surface->GetYaxis()->SetTitleOffset(1);
   chisq_surface->GetYaxis()->SetTitleSize(0.05);
   chisq_surface->GetXaxis()->SetTitleSize(0.05);

   return chisq_surface;            
}


