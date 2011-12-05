//____________________________________________________________________________ 
/*
 Changes since v1.0.0:

 @ 16/03/2011 - CA
   Removed the numu10a namespace. Renamed: OscillationFitter->NuMuDisappFitter.
   The FitInputs classs is no longer nested in the fitter class.
 @ 08/04/2011 - CA
   Update MINOS contour to the 2011 one.
 @ 05/07/2011 - CA
   PhysParamList was renamed to FitParamList so that distinction between the two
   SuperkExpectation interfaces is clear. Adapt the fitter.
   Use FitterABC as base class.
 @ 28/07/2011 - CA
   Added method to get MINUIT err matrix status code.
 @ 28/09/2011 - CA
   Actual fitting code factorized out of the existing fitters and put in
   FitterKernel class.  
 @ 06/10/2011 - CA
   Functionality of Write() was spread to Analyze(), GetResults(), Save() to
   better serve the different needs of various toy-MC fitting studies.
 @ 05/12/2011 - CA
   Add capability to use either sin^2(2*theta_ij) or sin^2(theta_ij) as the
   parameter of choice for describing neutrino mixing.
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

#include "Superk1RmuErecoFitter.h" 
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
Superk1RmuErecoFitter::Superk1RmuErecoFitter() :
FitterKernel() 
{

}
//____________________________________________________________________________ 
Superk1RmuErecoFitter::Superk1RmuErecoFitter(FitInputs * fi, Option_t * opt) :
FitterKernel(fi,opt)
{

}
//____________________________________________________________________________ 
Superk1RmuErecoFitter::~Superk1RmuErecoFitter() 
{

}
//____________________________________________________________________________ 
void Superk1RmuErecoFitter::Analyze(Option_t * option) 
{
  util::SetDefaultStyle();

  TString opt = option;
  cout << "Superk1RmuErecoFitter/Analyze >> Option: " << opt.Data() << endl;

  //
  // Print best-fit results used for generating the plots
  //
  fFitInputs->FitParams().ShowOnlyBestFitValuesInNextPrint();
  cout << "Superk1RmuErecoFitter/Analyze >> "
       << "Using the following best-fit results: " << fFitInputs->FitParams();

  //
  // Declare (potential) outpts
  //
  TH1D * data                   = 0;
  TH1D * bestfit_tot            = 0;
  TH1D * bestfit_numuccqe       = 0;
  TH1D * bestfit_numucc1pi      = 0;
  TH1D * bestfit_numuccoth      = 0;
  TH1D * bestfit_numunc         = 0;
  TH1D * bestfit_numubarccqe    = 0;
  TH1D * bestfit_numubarcc1pi   = 0;
  TH1D * bestfit_numubarccoth   = 0;
  TH1D * bestfit_numubarnc      = 0;
  TH1D * bestfit_nueccqe        = 0;
  TH1D * bestfit_nuecc1pi       = 0;
  TH1D * bestfit_nueccoth       = 0;
  TH1D * bestfit_nuenc          = 0;
  TH1D * bestfit_signueccqe     = 0;
  TH1D * bestfit_signuecc1pi    = 0;
  TH1D * bestfit_signueccoth    = 0;
  TH2D * chisq_surface          = 0;
  TCanvas * c1 = 0;
  TCanvas * c2 = 0;

  //
  // Write the best-fit spectra (total & and its components)
  //
  if(opt.Contains("BestFitSpectra")) {

     // `best-fit' --> `current' values
     fFitInputs->FitParams().CopyBestFitToCurrValues();

     // generate expectation
     fFitInputs->Expectation().Generate(fFitInputs->FitParams());

     // copy & format Ereco histograms
     bestfit_tot            = new TH1D(*fFitInputs->Expectation().CurrPredictedErecoSpectrum(kSKsFCFV1RingMuLike));
     bestfit_numuccqe       = new TH1D(*fFitInputs->Expectation().CurrPredictedErecoSpectrum(kSKsFCFV1RingMuLike, kSKmNuMuCCQE));
     bestfit_numucc1pi      = new TH1D(*fFitInputs->Expectation().CurrPredictedErecoSpectrum(kSKsFCFV1RingMuLike, kSKmNuMuCC1pi));
     bestfit_numuccoth      = new TH1D(*fFitInputs->Expectation().CurrPredictedErecoSpectrum(kSKsFCFV1RingMuLike, kSKmNuMuCCoth_mu));
     bestfit_numunc         = new TH1D(*fFitInputs->Expectation().CurrPredictedErecoSpectrum(kSKsFCFV1RingMuLike, kSKmNuMuNC));
     bestfit_numubarccqe    = new TH1D(*fFitInputs->Expectation().CurrPredictedErecoSpectrum(kSKsFCFV1RingMuLike, kSKmNuMuBarCCQE));
     bestfit_numubarcc1pi   = new TH1D(*fFitInputs->Expectation().CurrPredictedErecoSpectrum(kSKsFCFV1RingMuLike, kSKmNuMuBarCC1pi));
     bestfit_numubarccoth   = new TH1D(*fFitInputs->Expectation().CurrPredictedErecoSpectrum(kSKsFCFV1RingMuLike, kSKmNuMuBarCCoth_mu));
     bestfit_numubarnc      = new TH1D(*fFitInputs->Expectation().CurrPredictedErecoSpectrum(kSKsFCFV1RingMuLike, kSKmNuMuBarNC));
     bestfit_nueccqe        = new TH1D(*fFitInputs->Expectation().CurrPredictedErecoSpectrum(kSKsFCFV1RingMuLike, kSKmNuECCQE));
     bestfit_nuecc1pi       = new TH1D(*fFitInputs->Expectation().CurrPredictedErecoSpectrum(kSKsFCFV1RingMuLike, kSKmNuECC1pi));
     bestfit_nueccoth       = new TH1D(*fFitInputs->Expectation().CurrPredictedErecoSpectrum(kSKsFCFV1RingMuLike, kSKmNuECCoth_mu));
     bestfit_nuenc          = new TH1D(*fFitInputs->Expectation().CurrPredictedErecoSpectrum(kSKsFCFV1RingMuLike, kSKmNuENC));
     bestfit_signueccqe     = new TH1D(*fFitInputs->Expectation().CurrPredictedErecoSpectrum(kSKsFCFV1RingMuLike, kSKmSigNuECCQE));
     bestfit_signuecc1pi    = new TH1D(*fFitInputs->Expectation().CurrPredictedErecoSpectrum(kSKsFCFV1RingMuLike, kSKmSigNuECC1pi));
     bestfit_signueccoth    = new TH1D(*fFitInputs->Expectation().CurrPredictedErecoSpectrum(kSKsFCFV1RingMuLike, kSKmSigNuECCoth));

     bestfit_tot            -> SetDirectory ( 0        );
     bestfit_numuccqe       -> SetDirectory ( 0        );
     bestfit_numucc1pi      -> SetDirectory ( 0        );
     bestfit_numuccoth      -> SetDirectory ( 0        );
     bestfit_numunc         -> SetDirectory ( 0        );
     bestfit_numubarccqe    -> SetDirectory ( 0        );
     bestfit_numubarcc1pi   -> SetDirectory ( 0        );
     bestfit_numubarccoth   -> SetDirectory ( 0        );
     bestfit_numubarnc      -> SetDirectory ( 0        );
     bestfit_nueccqe        -> SetDirectory ( 0        );
     bestfit_nuecc1pi       -> SetDirectory ( 0        );
     bestfit_nueccoth       -> SetDirectory ( 0        );
     bestfit_nuenc          -> SetDirectory ( 0        );
     bestfit_signueccqe     -> SetDirectory ( 0        );
     bestfit_signuecc1pi    -> SetDirectory ( 0        );
     bestfit_signueccoth    -> SetDirectory ( 0        );

     bestfit_numucc1pi      -> SetFillStyle ( 3544     );
     bestfit_numuccoth      -> SetFillStyle ( 3244     );
     bestfit_numunc         -> SetFillStyle ( 3006     );
     bestfit_numubarccqe    -> SetFillStyle ( 3001     );
     bestfit_numubarcc1pi   -> SetFillStyle ( 3345     );
     bestfit_numubarccoth   -> SetFillStyle ( 3354     );
     bestfit_numubarnc      -> SetFillStyle ( 3006     );
     bestfit_nueccqe        -> SetFillStyle ( 3001     );
     bestfit_nuecc1pi       -> SetFillStyle ( 3345     );
     bestfit_nueccoth       -> SetFillStyle ( 3354     );
     bestfit_nuenc          -> SetFillStyle ( 3006     );
     bestfit_signueccqe     -> SetFillStyle ( 3001     );
     bestfit_signuecc1pi    -> SetFillStyle ( 3345     );
     bestfit_signueccoth    -> SetFillStyle ( 3354     );

     bestfit_numuccqe       -> SetLineColor ( kBlack   );

     bestfit_numucc1pi      -> SetFillColor ( kBlack   );
     bestfit_numuccoth      -> SetFillColor ( kBlack   );
     bestfit_numunc         -> SetFillColor ( kRed     );
     bestfit_numubarccqe    -> SetFillColor ( kAzure-3 );
     bestfit_numubarcc1pi   -> SetFillColor ( kAzure-3 );
     bestfit_numubarccoth   -> SetFillColor ( kAzure-3 );
     bestfit_numubarnc      -> SetFillColor ( kAzure-3 );
     bestfit_nueccqe        -> SetFillColor ( kGreen+1 );
     bestfit_nuecc1pi       -> SetFillColor ( kGreen+1 );
     bestfit_nueccoth       -> SetFillColor ( kGreen+1 );
     bestfit_nuenc          -> SetFillColor ( kGreen+1 );
     bestfit_signueccqe     -> SetFillColor ( kGreen-3 );
     bestfit_signuecc1pi    -> SetFillColor ( kGreen-3 );
     bestfit_signueccoth    -> SetFillColor ( kGreen-3 );

     bestfit_tot            -> SetLineWidth ( 2 );
     bestfit_numuccqe       -> SetLineWidth ( 1 );
     bestfit_numucc1pi      -> SetLineWidth ( 1 );
     bestfit_numuccoth      -> SetLineWidth ( 1 );
     bestfit_numunc         -> SetLineWidth ( 1 );
     bestfit_numubarccqe    -> SetLineWidth ( 1 );
     bestfit_numubarcc1pi   -> SetLineWidth ( 1 );
     bestfit_numubarccoth   -> SetLineWidth ( 1 );
     bestfit_numubarnc      -> SetLineWidth ( 1 );
     bestfit_nueccqe        -> SetLineWidth ( 1 );
     bestfit_nuecc1pi       -> SetLineWidth ( 1 );
     bestfit_nueccoth       -> SetLineWidth ( 1 );
     bestfit_nuenc          -> SetLineWidth ( 1 );
     bestfit_signueccqe     -> SetLineWidth ( 1 );
     bestfit_signuecc1pi    -> SetLineWidth ( 1 );
     bestfit_signueccoth    -> SetLineWidth ( 1 );

     bestfit_tot            -> SetName("bestfit_tot");
     bestfit_numuccqe       -> SetName("bestfit_numuccqe");
     bestfit_numucc1pi      -> SetName("bestfit_numucc1pi");
     bestfit_numuccoth      -> SetName("bestfit_numuccoth");
     bestfit_numunc         -> SetName("bestfit_numunc");
     bestfit_numubarccqe    -> SetName("bestfit_numubarcc");
     bestfit_numubarcc1pi   -> SetName("bestfit_numubarcc1pi");
     bestfit_numubarccoth   -> SetName("bestfit_numubarccoth");
     bestfit_numubarnc      -> SetName("bestfit_numubarnc");
     bestfit_nueccqe        -> SetName("bestfit_nueccqe");
     bestfit_nuecc1pi       -> SetName("bestfit_nuecc1pi");
     bestfit_nueccoth       -> SetName("bestfit_nueccoth");
     bestfit_nuenc          -> SetName("bestfit_nuenc");
     bestfit_signueccqe     -> SetName("bestfit_signueccqe");
     bestfit_signuecc1pi    -> SetName("bestfit_signuecc1pi");
     bestfit_signueccoth    -> SetName("bestfit_signueccoth");

     fResults -> Add (bestfit_tot);
     fResults -> Add (bestfit_numuccqe);
     fResults -> Add (bestfit_numucc1pi);
     fResults -> Add (bestfit_numuccoth);
     fResults -> Add (bestfit_numunc);
     fResults -> Add (bestfit_numubarccqe);
     fResults -> Add (bestfit_numubarcc1pi);
     fResults -> Add (bestfit_numubarccoth);
     fResults -> Add (bestfit_numubarnc);
     fResults -> Add (bestfit_nueccqe);
     fResults -> Add (bestfit_nuecc1pi);
     fResults -> Add (bestfit_nueccoth);
     fResults -> Add (bestfit_nuenc);
     fResults -> Add (bestfit_signueccqe);
     fResults -> Add (bestfit_signuecc1pi);
     fResults -> Add (bestfit_signueccoth);
  }

  //
  // Write the fit data
  //
  if(opt.Contains("Observation")) {
     data = new TH1D(*fFitInputs->Observation().ErecoSpectrum(kSKsFCFV1RingMuLike));
     data -> SetDirectory(0);
     data -> SetMarkerStyle(20);
     data -> SetMarkerSize(1);
     data -> GetXaxis()->SetTitle("E_{reco} (GeV)");
     data -> GetYaxis()->SetTitle("Number of 1Rmu events");
     data -> SetName("data");
     fResults -> Add(data);
  }

  //
  // Write out the chisq surface
  //

  if(opt.Contains("ChisqSurface")) {

     PhysHypothesis_t phyp = fFitInputs->FitParams().PhysicsHypothesis();

     if(phyp == kHypNuMuDisapp2F || phyp == kHypOsc3F) {
        chisq_surface = this->BuildChisqSurfaceNuMuDisap2ParamFit();
     }
     else
     if(phyp == kHypDecay) {   
        chisq_surface = this->BuildChisqSurfaceNuMuDecay2ParamFit();
     }
     else
     if(phyp == kHypDecoh) {   
        chisq_surface = this->BuildChisqSurfaceNuMuDecoh2ParamFit();
     } 
     else 
     {
        cerr << "Superk1RmuErecoFitter/Analyze >> Can't handle hypothesis: " 
             << util::AsString(phyp) << endl;
        exit(1);
     }

     chisq_surface -> SetName("chisq_surface");
     fResults -> Add(chisq_surface);
  }

  //
  // Save canvases with some standard plots
  //
  if(opt.Contains("Canvases")) {
    if(opt.Contains("BestFitSpectra") && opt.Contains("Observation")) {

        THStack * contrib_stack = new THStack("stack","spectrum contributions");
        contrib_stack->Add(bestfit_signuecc1pi);
        contrib_stack->Add(bestfit_signueccoth);
        contrib_stack->Add(bestfit_signueccqe);
        contrib_stack->Add(bestfit_nuenc);
        contrib_stack->Add(bestfit_nuecc1pi);
        contrib_stack->Add(bestfit_nueccoth);
        contrib_stack->Add(bestfit_nueccqe);
        contrib_stack->Add(bestfit_numubarnc);
        contrib_stack->Add(bestfit_numubarcc1pi);
        contrib_stack->Add(bestfit_numubarccoth);
        contrib_stack->Add(bestfit_numubarccqe);
        contrib_stack->Add(bestfit_numunc);
        contrib_stack->Add(bestfit_numuccoth);
        contrib_stack->Add(bestfit_numucc1pi);
        contrib_stack->Add(bestfit_numuccqe);

        c1 = new TCanvas("c1","Data and best-fit spectrum (total)",20,20,500,500);
        c1->SetBorderMode(0);
        c1->SetFillColor (0);
        data        -> Draw("perr");
        bestfit_tot -> Draw("same");
        c1->Update();

        c2 = new TCanvas("c2","Data and best-fit spectrum (total & components)",20,20,500,500);
        c2->SetBorderMode(0);
        c2->SetFillColor (0);
        data          -> Draw("perr");
        bestfit_tot   -> Draw("same");
        contrib_stack -> Draw("same");
        TLegend * spectrum_legend = new TLegend(0.6,0.5,0.9,0.9);
        spectrum_legend->SetFillColor (0);
        spectrum_legend->SetBorderSize(0);
        spectrum_legend->AddEntry(data,                    "data",                                           "P" );
        spectrum_legend->AddEntry(bestfit_tot,             "best fit",                                       "L" );
        spectrum_legend->AddEntry(bestfit_numuccqe,        "best fit (#nu_{#mu} CCQE)",                      "L" );
        spectrum_legend->AddEntry(bestfit_numucc1pi,       "best fit (#nu_{#mu} CC1#pi)",                    "LF");
        spectrum_legend->AddEntry(bestfit_numuccoth,       "best fit (#nu_{#mu} CC other)",                  "LF");
        spectrum_legend->AddEntry(bestfit_numunc,          "best fit (#nu_{#mu}/#nu_{#tau} NC)",             "LF");
        spectrum_legend->AddEntry(bestfit_numubarccqe,     "best fit (#bar{#nu_{#mu}} CCQE)",                "LF");
        spectrum_legend->AddEntry(bestfit_numubarcc1pi,    "best fit (#bar{#nu_{#mu}} CC1#pi)",              "LF");
        spectrum_legend->AddEntry(bestfit_numubarccoth,    "best fit (#bar{#nu_{#mu}} CC other)",            "LF");
        spectrum_legend->AddEntry(bestfit_numubarnc,       "best fit (#bar{#nu_{#mu}}/#bar{#nu_{#tau}} NC)", "LF");
        spectrum_legend->AddEntry(bestfit_nueccqe,         "best fit (#nu_{e} CCQE)",                        "LF");
        spectrum_legend->AddEntry(bestfit_nuecc1pi,        "best fit (#nu_{e} CC1#pi)",                      "LF");
        spectrum_legend->AddEntry(bestfit_nueccoth,        "best fit (#nu_{e} CC other)",                    "LF");
        spectrum_legend->AddEntry(bestfit_nuenc,           "best fit (#nu_{e} NC)",                          "LF");
        spectrum_legend->AddEntry(bestfit_signueccqe,      "best fit (signal-#nu_{e} CCQE)",                 "LF");
        spectrum_legend->AddEntry(bestfit_signuecc1pi,     "best fit (signal-#nu_{e} CC1#pi)",               "LF");
        spectrum_legend->AddEntry(bestfit_signueccoth,     "best fit (signal-#nu_{e} CC other)",             "LF");
        spectrum_legend->Draw();
        c2->Update();

        c1->SetName("c1");
        c2->SetName("c2");

        fResults->Add(c1);
        fResults->Add(c2);

   } // BestFitSpectra & Observation
  } // Canvases

}
//____________________________________________________________________________ 
void Superk1RmuErecoFitter::Save(const char * filename) 
{
  util::SetDefaultStyle();

  cout << "Superk1RmuErecoFitter/Save >> Writing fit results in: " 
       << filename << endl;
  TFile * file_out = new TFile(filename, "recreate");

  TH1D * data                   = (TH1D*) fResults->FindObject("data");
  TH1D * bestfit_tot            = (TH1D*) fResults->FindObject("bestfit_tot");
  TH1D * bestfit_numuccqe       = (TH1D*) fResults->FindObject("bestfit_numuccqe");
  TH1D * bestfit_numucc1pi      = (TH1D*) fResults->FindObject("bestfit_numucc1pi");
  TH1D * bestfit_numuccoth      = (TH1D*) fResults->FindObject("bestfit_numuccoth");
  TH1D * bestfit_numunc         = (TH1D*) fResults->FindObject("bestfit_numunc");
  TH1D * bestfit_numubarccqe    = (TH1D*) fResults->FindObject("bestfit_numubarccqe");
  TH1D * bestfit_numubarcc1pi   = (TH1D*) fResults->FindObject("bestfit_numubarcc1pi");
  TH1D * bestfit_numubarccoth   = (TH1D*) fResults->FindObject("bestfit_numubarccoth");
  TH1D * bestfit_numubarnc      = (TH1D*) fResults->FindObject("bestfit_numubarnc");
  TH1D * bestfit_nueccqe        = (TH1D*) fResults->FindObject("bestfit_nueccqe");
  TH1D * bestfit_nuecc1pi       = (TH1D*) fResults->FindObject("bestfit_nuecc1pi");
  TH1D * bestfit_nueccoth       = (TH1D*) fResults->FindObject("bestfit_nueccoth");
  TH1D * bestfit_nuenc          = (TH1D*) fResults->FindObject("bestfit_nuenc");
  TH1D * bestfit_signueccqe     = (TH1D*) fResults->FindObject("bestfit_signueccqe");
  TH1D * bestfit_signuecc1pi    = (TH1D*) fResults->FindObject("bestfit_signuecc1pi");
  TH1D * bestfit_signueccoth    = (TH1D*) fResults->FindObject("bestfit_signueccoth");
  TH2D * chisq_surface          = (TH2D*) fResults->FindObject("chisq_surface");
  TCanvas * c1 = (TCanvas*) fResults->FindObject("c1");
  TCanvas * c2 = (TCanvas*) fResults->FindObject("c2");

  if ( bestfit_tot          ) { bestfit_tot            -> Write(); }
  if ( bestfit_numuccqe     ) { bestfit_numuccqe       -> Write(); }
  if ( bestfit_numucc1pi    ) { bestfit_numucc1pi      -> Write(); }
  if ( bestfit_numuccoth    ) { bestfit_numuccoth      -> Write(); }
  if ( bestfit_numunc       ) { bestfit_numunc         -> Write(); }
  if ( bestfit_numubarccqe  ) { bestfit_numubarccqe    -> Write(); }
  if ( bestfit_numubarcc1pi ) { bestfit_numubarcc1pi   -> Write(); }
  if ( bestfit_numubarccoth ) { bestfit_numubarccoth   -> Write(); }
  if ( bestfit_numubarnc    ) { bestfit_numubarnc      -> Write(); }
  if ( bestfit_nueccqe      ) { bestfit_nueccqe        -> Write(); }
  if ( bestfit_nuecc1pi     ) { bestfit_nuecc1pi       -> Write(); }
  if ( bestfit_nueccoth     ) { bestfit_nueccoth       -> Write(); }
  if ( bestfit_nuenc        ) { bestfit_nuenc          -> Write(); }
  if ( bestfit_signueccqe   ) { bestfit_signueccqe     -> Write(); }
  if ( bestfit_signuecc1pi  ) { bestfit_signuecc1pi    -> Write(); }
  if ( bestfit_signueccoth  ) { bestfit_signueccoth    -> Write(); }
  if ( data                 ) { data                   -> Write(); }
  if ( chisq_surface        ) { chisq_surface          -> Write(); }
  if ( c1                   ) { c1                     -> Write(); }
  if ( c2                   ) { c2                     -> Write(); }

  file_out->Close();
  delete file_out;
}
//____________________________________________________________________________ 
TH2D * Superk1RmuErecoFitter::BuildChisqSurfaceNuMuDisap2ParamFit(void) 
{
   FitParamList & fp = fFitInputs->FitParams();

   double mixing23_min  = fp.MinValue(conv::kPhypMixing23);
   double mixing23_max  = fp.MaxValue(conv::kPhypMixing23);
   double mixing23_step = fp.StepSize(conv::kPhypMixing23);
   double dm23sq_min    = fp.MinValue(conv::kPhypDm23sq);
   double dm23sq_max    = fp.MaxValue(conv::kPhypDm23sq);
   double dm23sq_step   = fp.StepSize(conv::kPhypDm23sq);

   // make sure that min and max are at the bin centre
   int mixing23_nbins = 1 + (int) ((mixing23_max - mixing23_min) / mixing23_step);
   int dm23sq_nbins   = 1 + (int) ((dm23sq_max -   dm23sq_min )  / dm23sq_step );

   TH2D * chisq_surface = new TH2D(
      "chisq_surface","#chi^{2} surface",
         mixing23_nbins, mixing23_min - 0.5*mixing23_step, mixing23_max + 0.5*mixing23_step,
         dm23sq_nbins, dm23sq_min - 0.5*dm23sq_step, dm23sq_max + 0.5*dm23sq_step);
   chisq_surface->SetDirectory(0);
   chisq_surface->GetXaxis()->SetRangeUser(0.,1.);

   for(int i = 1; i <= chisq_surface->GetNbinsX(); i++) {
     for(int j = 1; j <= chisq_surface->GetNbinsY(); j++) {
         double mixing23 = chisq_surface->GetXaxis()->GetBinCenter(i);
         double dm23sq   = chisq_surface->GetYaxis()->GetBinCenter(j);
         fp.Update( conv::kPhypMixing23, mixing23 );
         fp.Update( conv::kPhypDm23sq,   dm23sq   );
         double chisq = fFitInputs->Chisquare();
         chisq_surface->SetBinContent(i,j,chisq);
     }//j
   }//i

   chisq_surface->GetXaxis()->SetTitle("sin^{2}(2#theta_{23})");
   chisq_surface->GetYaxis()->SetTitle("|#Deltam_{32}^{2}| (eV^{2}/c^{4})");
   chisq_surface->GetYaxis()->SetTitleOffset(1);
   chisq_surface->GetYaxis()->SetTitleSize(0.05);
   chisq_surface->GetXaxis()->SetTitleSize(0.05);

   return chisq_surface;            
}
//____________________________________________________________________________ 
TH2D * Superk1RmuErecoFitter::BuildChisqSurfaceNuMuDecay2ParamFit(void) 
{
   FitParamList & fp = fFitInputs->FitParams();

   double sinsq2theta_min  = fp.MinValue(conv::kPhypSinsq2ThetaDecay);
   double sinsq2theta_max  = fp.MaxValue(conv::kPhypSinsq2ThetaDecay);
   double sinsq2theta_step = fp.StepSize(conv::kPhypSinsq2ThetaDecay);
   double alpha_min        = fp.MinValue(conv::kPhypAlphaDecay);
   double alpha_max        = fp.MaxValue(conv::kPhypAlphaDecay);
   double alpha_step       = fp.StepSize(conv::kPhypAlphaDecay);

   // make sure that min and max are at the bin centre
   int sinsq2theta_nbins = 1 + (int) ((sinsq2theta_max - sinsq2theta_min) / sinsq2theta_step);
   int alpha_nbins       = 1 + (int) ((alpha_max -       alpha_min )      / alpha_step );

   TH2D * chisq_surface = new TH2D(
      "chisq_surface","#chi^{2} surface",
         sinsq2theta_nbins, sinsq2theta_min - 0.5*sinsq2theta_step, sinsq2theta_max + 0.5*sinsq2theta_step,
         alpha_nbins, alpha_min - 0.5*alpha_step, alpha_max + 0.5*alpha_step);
   chisq_surface->SetDirectory(0);
   chisq_surface->GetXaxis()->SetRangeUser(0.,1.);

   for(int i = 1; i <= chisq_surface->GetNbinsX(); i++) {
     for(int j = 1; j <= chisq_surface->GetNbinsY(); j++) {
         double sinsq2theta = chisq_surface->GetXaxis()->GetBinCenter(i);
         double alpha       = chisq_surface->GetYaxis()->GetBinCenter(j);
         fp.Update( conv::kPhypSinsq2ThetaDecay, sinsq2theta );
         fp.Update( conv::kPhypAlphaDecay,       alpha       );
         double chisq = fFitInputs->Chisquare();
         chisq_surface->SetBinContent(i,j,chisq);
     }//j
   }//i

   chisq_surface->GetXaxis()->SetTitle("sin^{2}(2#theta)");
   chisq_surface->GetYaxis()->SetTitle("#alpha (eV^{2})");
   chisq_surface->GetYaxis()->SetTitleOffset(1);
   chisq_surface->GetYaxis()->SetTitleSize(0.05);
   chisq_surface->GetXaxis()->SetTitleSize(0.05);

   return chisq_surface;                                                             
}
//____________________________________________________________________________ 
TH2D * Superk1RmuErecoFitter::BuildChisqSurfaceNuMuDecoh2ParamFit(void) 
{
   FitParamList & fp = fFitInputs->FitParams();

   double sinsq2theta_min  = fp.MinValue(conv::kPhypSinsq2ThetaDecoh);
   double sinsq2theta_max  = fp.MaxValue(conv::kPhypSinsq2ThetaDecoh);
   double sinsq2theta_step = fp.StepSize(conv::kPhypSinsq2ThetaDecoh);
   double musq_min         = fp.MinValue(conv::kPhypMusqDecoh);
   double musq_max         = fp.MaxValue(conv::kPhypMusqDecoh);
   double musq_step        = fp.StepSize(conv::kPhypMusqDecoh);

   // make sure that min and max are at the bin centre
   int sinsq2theta_nbins = 1 + (int) ((sinsq2theta_max - sinsq2theta_min) / sinsq2theta_step);
   int musq_nbins        = 1 + (int) ((musq_max -        musq_min )       / musq_step );

   TH2D * chisq_surface = new TH2D(
      "chisq_surface","#chi^{2} surface",
         sinsq2theta_nbins, sinsq2theta_min - 0.5*sinsq2theta_step, sinsq2theta_max + 0.5*sinsq2theta_step,
         musq_nbins, musq_min - 0.5*musq_step, musq_max + 0.5*musq_step);
   chisq_surface->SetDirectory(0);
   chisq_surface->GetXaxis()->SetRangeUser(0.,1.);

   for(int i = 1; i <= chisq_surface->GetNbinsX(); i++) {
     for(int j = 1; j <= chisq_surface->GetNbinsY(); j++) {
         double sinsq2theta = chisq_surface->GetXaxis()->GetBinCenter(i);
         double musq        = chisq_surface->GetYaxis()->GetBinCenter(j);
         fp.Update( conv::kPhypSinsq2ThetaDecoh,  sinsq2theta );
         fp.Update( conv::kPhypMusqDecoh,         musq        );
         double chisq = fFitInputs->Chisquare();
         chisq_surface->SetBinContent(i,j,chisq);
     }//j
   }//i

   chisq_surface->GetXaxis()->SetTitle("sin^{2}(2#theta)");
   chisq_surface->GetYaxis()->SetTitle("#mu^{2} (eV^{2})");
   chisq_surface->GetYaxis()->SetTitleOffset(1);
   chisq_surface->GetYaxis()->SetTitleSize(0.05);
   chisq_surface->GetXaxis()->SetTitleSize(0.05);

   return chisq_surface;                                                             
}
//____________________________________________________________________________ 



