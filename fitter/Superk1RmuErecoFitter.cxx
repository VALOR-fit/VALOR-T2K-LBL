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
void Superk1RmuErecoFitter::Write(Option_t * option, const char * outtag) 
{
  util::SetDefaultStyle();

  TString opt = option;
  cout << "Superk1RmuErecoFitter/Write >> Input option: " 
       << opt.Data() << endl;

  TFile * file_out = 0;
  if(outtag) {
     ostringstream outroot;
     outroot << outtag << ".root";
     cout << "Superk1RmuErecoFitter/Write >> Writing fit results in: " 
          << outroot.str() << endl; 
     file_out = new TFile(outroot.str().c_str(), "recreate");
  }

  //
  // Print best-fit results used for generating the plots
  //
  fFitInputs->FitParams().ShowOnlyBestFitValuesInNextPrint();
  cout << "Superk1RmuErecoFitter/Write >> Using the following best-fit results: " 
       << fFitInputs->FitParams();

  //
  // Declare (potential) outpts
  //
  TH1D * data                   = 0;
  TH1D * bestfit_tot            = 0;
  TH1D * bestfit_numuccqe       = 0;
  TH1D * bestfit_numucc1pi      = 0;
  TH1D * bestfit_numuccoth      = 0;
  TH1D * bestfit_nc             = 0;
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
     bestfit_numuccoth      = new TH1D(*fFitInputs->Expectation().CurrPredictedErecoSpectrum(kSKsFCFV1RingMuLike, kSKmNuMuCCoth));
     bestfit_nc             = new TH1D(*fFitInputs->Expectation().CurrPredictedErecoSpectrum(kSKsFCFV1RingMuLike, kSKmNC));
     bestfit_numubarccqe    = new TH1D(*fFitInputs->Expectation().CurrPredictedErecoSpectrum(kSKsFCFV1RingMuLike, kSKmNuMuBarCCQE));
     bestfit_numubarcc1pi   = new TH1D(*fFitInputs->Expectation().CurrPredictedErecoSpectrum(kSKsFCFV1RingMuLike, kSKmNuMuBarCC1pi));
     bestfit_numubarccoth   = new TH1D(*fFitInputs->Expectation().CurrPredictedErecoSpectrum(kSKsFCFV1RingMuLike, kSKmNuMuBarCCoth));
     bestfit_numubarnc      = new TH1D(*fFitInputs->Expectation().CurrPredictedErecoSpectrum(kSKsFCFV1RingMuLike, kSKmNuMuBarNC));
     bestfit_nueccqe        = new TH1D(*fFitInputs->Expectation().CurrPredictedErecoSpectrum(kSKsFCFV1RingMuLike, kSKmNuECCQE));
     bestfit_nuecc1pi       = new TH1D(*fFitInputs->Expectation().CurrPredictedErecoSpectrum(kSKsFCFV1RingMuLike, kSKmNuECC1pi));
     bestfit_nueccoth       = new TH1D(*fFitInputs->Expectation().CurrPredictedErecoSpectrum(kSKsFCFV1RingMuLike, kSKmNuECCoth));
     bestfit_nuenc          = new TH1D(*fFitInputs->Expectation().CurrPredictedErecoSpectrum(kSKsFCFV1RingMuLike, kSKmNuENC));
     bestfit_signueccqe     = new TH1D(*fFitInputs->Expectation().CurrPredictedErecoSpectrum(kSKsFCFV1RingMuLike, kSKmSigNuECCQE));
     bestfit_signuecc1pi    = new TH1D(*fFitInputs->Expectation().CurrPredictedErecoSpectrum(kSKsFCFV1RingMuLike, kSKmSigNuECC1pi));
     bestfit_signueccoth    = new TH1D(*fFitInputs->Expectation().CurrPredictedErecoSpectrum(kSKsFCFV1RingMuLike, kSKmSigNuECCoth));

     bestfit_tot            -> SetDirectory ( 0        );
     bestfit_numuccqe       -> SetDirectory ( 0        );
     bestfit_numucc1pi      -> SetDirectory ( 0        );
     bestfit_numuccoth      -> SetDirectory ( 0        );
     bestfit_nc             -> SetDirectory ( 0        );
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
     bestfit_nc             -> SetFillStyle ( 3006     );
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
     bestfit_nc             -> SetFillColor ( kRed     );
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
     bestfit_nc             -> SetLineWidth ( 1 );
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

     if(file_out) {
       bestfit_tot            -> Write("bestfit_tot");
       bestfit_numuccqe       -> Write("bestfit_numuccqe");
       bestfit_numucc1pi      -> Write("bestfit_numucc1pi");
       bestfit_numuccoth      -> Write("bestfit_numuccoth");
       bestfit_nc             -> Write("bestfit_numunc");
       bestfit_numubarccqe    -> Write("bestfit_numubarcc");
       bestfit_numubarcc1pi   -> Write("bestfit_numubarcc1pi");
       bestfit_numubarccoth   -> Write("bestfit_numubarccoth");
       bestfit_numubarnc      -> Write("bestfit_numubarnc");
       bestfit_nueccqe        -> Write("bestfit_nueccqe");
       bestfit_nuecc1pi       -> Write("bestfit_nuecc1pi");
       bestfit_nueccoth       -> Write("bestfit_nueccoth");
       bestfit_nuenc          -> Write("bestfit_nuenc");
       bestfit_signueccqe     -> Write("bestfit_signueccqe");
       bestfit_signuecc1pi    -> Write("bestfit_signuecc1pi");
       bestfit_signueccoth    -> Write("bestfit_signueccoth");
     }
  }

  //
  // Write the fit data
  //
  if(opt.Contains("Observation")) {

     // get & format data
     data = new TH1D(*fFitInputs->Observation().ErecoSpectrum(kSKsFCFV1RingMuLike));
     data -> SetDirectory(0);
     data -> SetMarkerStyle(20);
     data -> SetMarkerSize(1);
     data -> GetXaxis()->SetTitle("E_{reco} (GeV)");
     data -> GetYaxis()->SetTitle("Number of 1Rmu events");
     if(file_out) {
        data -> Write("data");
     }
  }

  //
  // Write out the chisq surface
  //

  if(opt.Contains("ChisqSurface")) {

     PhysHypothesis_t phyp = fFitInputs->FitParams().PhysicsHypothesis();

     if(phyp == kHypNuMuDisapp2F || phyp == kHypOsc3F) 
     {
        chisq_surface = this->BuildChisqSurfaceNuMuDisap2ParamFit();
     }
     else
     if(phyp == kHypDecay) 
     {   
        chisq_surface = this->BuildChisqSurfaceNuMuDecay2ParamFit();
     }
     else
     if(phyp == kHypDecoh) {   
        chisq_surface = this->BuildChisqSurfaceNuMuDecoh2ParamFit();
     } 
     else 
     {
        cerr << "Superk1RmuErecoFitter/Write >> Can not handle hypothesis: " 
             << util::AsString(phyp) << endl;
        exit(1);
     }

     if(file_out) {
       chisq_surface -> Write("chisq_surface");
     }
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
        contrib_stack->Add(bestfit_nc);
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
        spectrum_legend->AddEntry(bestfit_nc,              "best fit (#nu_{#mu}/#nu_{#tau} NC)",             "LF");
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

        if(file_out) {
          c1->Write("c1");
          c2->Write("c2");
        }

   } // BestFitSpectra & Observation
  } // Canvases

  if(file_out) {
    file_out->Close();
    delete file_out;
  }
}
//____________________________________________________________________________ 
TH2D * Superk1RmuErecoFitter::BuildChisqSurfaceNuMuDisap2ParamFit(void) 
{
   FitParamList & fp = fFitInputs->FitParams();

   double sinsq2theta23_min  = fp.MinValue(conv::kPhypSinsq2Theta23);
   double sinsq2theta23_max  = fp.MaxValue(conv::kPhypSinsq2Theta23);
   double sinsq2theta23_step = fp.StepSize(conv::kPhypSinsq2Theta23);
   double dm23sq_min         = fp.MinValue(conv::kPhypDm23sq);
   double dm23sq_max         = fp.MaxValue(conv::kPhypDm23sq);
   double dm23sq_step        = fp.StepSize(conv::kPhypDm23sq);

   // make sure that min and max are at the bin centre
   int sinsq2theta23_nbins = 1 + (int) ((sinsq2theta23_max - sinsq2theta23_min) / sinsq2theta23_step);
   int dm23sq_nbins        = 1 + (int) ((dm23sq_max -        dm23sq_min )       / dm23sq_step );

   TH2D * chisq_surface = new TH2D(
      "chisq_surface","#chi^{2} surface",
         sinsq2theta23_nbins, sinsq2theta23_min - 0.5*sinsq2theta23_step, sinsq2theta23_max + 0.5*sinsq2theta23_step,
         dm23sq_nbins, dm23sq_min - 0.5*dm23sq_step, dm23sq_max + 0.5*dm23sq_step);
   chisq_surface->SetDirectory(0);
   chisq_surface->GetXaxis()->SetRangeUser(0.,1.);

   for(int i = 1; i <= chisq_surface->GetNbinsX(); i++) {
     for(int j = 1; j <= chisq_surface->GetNbinsY(); j++) {
         double sinsq2theta23 = chisq_surface->GetXaxis()->GetBinCenter(i);
         double dm23sq        = chisq_surface->GetYaxis()->GetBinCenter(j);
         fp.Update( conv::kPhypSinsq2Theta23, sinsq2theta23 );
         fp.Update( conv::kPhypDm23sq,        dm23sq        );
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



/*

//
// SOME OLD CODE - REMOVE
//


     int bin_sinsq_bestfit = fChisqSurface->GetXaxis()->FindBin(sinsq2theta23_bestfit);
     int bin_dmsq_bestfit = fChisqSurface->GetYaxis()->FindBin(dm23sq_bestfit);
     double chisq_min = fChisqSurface->GetBinContent(bin_sinsq_bestfit,bin_dmsq_bestfit);
     TH2D * chisq_surf_shift = new TH2D(*fChisqSurface); 
     chisq_surf_shift->SetDirectory(0); 
     for(int i = 1; i <= fChisqSurface->GetNbinsX(); i++) {
          for(int j = 1; j <= fChisqSurface->GetNbinsY(); j++) {
             double chisq = fChisqSurface->GetBinContent(i,j); 
             chisq_surf_shift->SetBinContent(i,j, TMath::Max(0., chisq-chisq_min));
          }//j
     }//i

     TGaxis::SetMaxDigits(2);
  
     chisq_surf_shift->GetXaxis()->SetTitle("sin^{2}(2#theta_{23})");
     chisq_surf_shift->GetYaxis()->SetTitle("|#Deltam_{32}^{2}| (eV^{2}/c^{4})");
     chisq_surf_shift->GetYaxis()->SetTitleOffset(1);
     chisq_surf_shift->GetYaxis()->SetTitleSize(0.05);
     chisq_surf_shift->GetXaxis()->SetTitleSize(0.05);

     const int nlev = 2;
     double levels[nlev] = {
         util::DeltaChisquare(conv::kCL68, fNFreeFitParams),
         util::DeltaChisquare(conv::kCL90, fNFreeFitParams)
     };

     chisq_surf_shift->SetContour(nlev,levels); 
     chisq_surf_shift->Draw("cont2"); 

     if(file_out) {
       fChisqSurface    -> Write("chisq_surface");
       chisq_surf_shift -> Write("chisq_minus_chisqmin_surface");
     }

     TMarker * best_fit_point =
          new TMarker(sinsq2theta23_bestfit, dm23sq_bestfit, 29);
     best_fit_point->SetMarkerSize(1.4);
     best_fit_point->Draw();
     if(file_out) {
        best_fit_point->Write("point_best_fit");
     }
     
     if(opt.Contains("Minos90CL")) {
        string contour_file = gSystem->Getenv("T2KOSCANA") + string("/inputs/external/contours/numu_disappearance.root");
        TFile f(contour_file.c_str(),"read");
        TGraph * gr = (TGraph*)f.Get("minos2011_90cl");
        TGraph * minos_90cl = new TGraph(*gr);
        minos_90cl->SetLineColor(kRed);
        minos_90cl->SetLineStyle(kSolid);
        minos_90cl->SetLineWidth(1);
        minos_90cl->Draw("l");
	f.Close();
     }

   }

     //
     // plot chisq = f (sinsq2theta23)
     //

     c->cd(3);
     TGraph * gr1 = 0;
     if(fOption.Contains("minuit")) {
        gMinuit->Command("SCAn 1");
        gr1 = (TGraph *) gMinuit->GetPlot();
         for(int i=0; i < gr1->GetN(); i++) {
             double x,y;
             gr1->GetPoint(i,x,y);
             gr1->SetPoint(i,x,y/fNdof);
         }
     }//minuit?
     else
     if(fOption.Contains("custom")) {
        int n = chisq_surf_shift->GetXaxis()->GetNbins();
        double * lchisq = new double[n];
        double * lsinsq = new double[n];
        for(int i=1; i<=n; i++) {
	    lchisq[i-1] = (fChisqSurface->GetBinContent(i,bin_dmsq_bestfit))/fNdof;
	    lsinsq[i-1] = fChisqSurface->GetXaxis()->GetBinCenter(i);
        }
        gr1 = new TGraph(n,lsinsq,lchisq);
     }//custom?
     xmin = ( gr1->GetX() )[TMath::LocMin(gr1->GetN(),gr1->GetX())];
     xmax = ( gr1->GetX() )[TMath::LocMax(gr1->GetN(),gr1->GetX())];
     ymin = ( gr1->GetY() )[TMath::LocMin(gr1->GetN(),gr1->GetY())];
     ymax = ( gr1->GetY() )[TMath::LocMax(gr1->GetN(),gr1->GetY())];
     TH1F* hframe1 = c->GetPad(3)->DrawFrame(xmin,0.8*ymin,xmax,1.2*ymax);
     hframe1->GetXaxis()->SetTitle("sin^{2}(2#theta_{23})");
     hframe1->GetYaxis()->SetTitle("#chi^{2}/ndof");
     gr1->Draw("l");

     //
     // plot chisq = f (Dm23sq)
     //
     c->cd(4);
     TGraph * gr2 = 0;
     if(fOption.Contains("minuit")) {
        gMinuit->Command("SCAn 2");
        gr2 = (TGraph *) gMinuit->GetPlot();
        for(int i=0; i < gr2->GetN(); i++) {
            double x,y;
            gr2->GetPoint(i,x,y);
            gr2->SetPoint(i,x,y/fNdof);
        }
     }//minuit?
     if(fOption.Contains("custom")) {
        int n = fChisqSurface->GetYaxis()->GetNbins();
        double * lchisq = new double[n];
        double * ldmsq = new double[n];
        for(int i=1; i<=n; i++) {
	    lchisq[i-1] = (fChisqSurface->GetBinContent(bin_sinsq_bestfit,i))/fNdof;
	    ldmsq [i-1] = fChisqSurface->GetYaxis()->GetBinCenter(i);
	
        }
        gr2 = new TGraph(n,ldmsq,lchisq);
     }//custom?
     xmin = ( gr2->GetX() )[TMath::LocMin(gr2->GetN(),gr2->GetX())];
     xmax = ( gr2->GetX() )[TMath::LocMax(gr2->GetN(),gr2->GetX())];
     ymin = ( gr2->GetY() )[TMath::LocMin(gr2->GetN(),gr2->GetY())];
     ymax = ( gr2->GetY() )[TMath::LocMax(gr2->GetN(),gr2->GetY())];
     TH1F* hframe2 = c->GetPad(4)->DrawFrame(xmin,0.8*ymin,xmax,1.2*ymax);
     hframe2->GetXaxis()->SetTitle("|#Deltam_{32}^{2}| (eV^{2}/c^{4})");
     hframe2->GetYaxis()->SetTitle("#chi^{2}/ndof");
     gr2->Draw("l");

     if(file_out) {
        gr1->Write("chisq_vs_sinsq2theta23");
        gr2->Write("chisq_vs_dm23sq");
     }
  }

  if(opt.Contains("ChisqSurface")) {
     if(fOption.Contains("minuit")) {
       this->BuildChisqSurface();
     } //minuit?
     fChisqSurface->Draw("colz");
  }

  if(opt.Contains("MarginalizedChisq")) {
     TH1D * hdm23sq = 
        this->MarginalizeChisqSurface("marginalize_sinsq2theta23+subtract_chisqmin");
     TH1D * hsinsq2theta23 = 
        this->MarginalizeChisqSurface("marginalize_dm23sq+subtract_chisqmin");
     hdm23sq ->GetXaxis()->SetTitle("|#Deltam_{32}^{2}| (eV^{2}/c^{4})");
     hsinsq2theta23 ->GetXaxis()->SetTitle("sin^{2}(2#theta_{23})");
     c->Divide(1,2);
     c->cd(1);
     hdm23sq->Draw("c");
     c->cd(2);
     hsinsq2theta23->Draw("c");
     if(file_out) {
        hdm23sq        -> Write("marginalized_chisq_dm23sq");
        hsinsq2theta23 -> Write("marginalized_chisq_sinsq2theta23");
     }
  }

  // Save canvas
  c->Update();
  if(outtag) {
    ostringstream outgif;
    ostringstream outroot;
    outgif << outtag << "_canvas.gif";
    outroot << outtag << "_canvas.root";
    c->SaveAs( outgif.str().c_str() );
    c->SaveAs( outroot.str().c_str() );
  }
*/

