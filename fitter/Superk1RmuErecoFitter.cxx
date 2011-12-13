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
 @ 13/12/2011 - CA
   Removed hardcoded list of modes for the 1Rmu sample. Considered modes are
   looked up from Conventions.h
*/ 
//____________________________________________________________________________

//#define __show_bestfit_values_used_for_postfit_analysis__

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
#ifdef __show_bestfit_values_used_for_postfit_analysis__
  fFitInputs->FitParams().ShowOnlyBestFitValuesInNextPrint();
  cout << "Superk1RmuErecoFitter/Analyze >> "
       << "Using the following best-fit results: " << fFitInputs->FitParams();
#endif

  SuperkSample_t s = kSKsFCFV1RingMuLike;
  int isample = util::FindSampleSlot(s);

  //
  // Declare (potential) outputs
  //
  TH1D * data = 0;
  TH1D * bestfit_tot = 0;
  TH1D * bestfit_mode[knm] = {0};
  TH2D * chisq_surface = 0;
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

     // copy & format Ereco spectra
     bestfit_tot = new TH1D(*fFitInputs->Expectation().CurrPredictedErecoSpectrum(s));
     bestfit_tot->SetDirectory (0);
     for(int imode = 0; imode < knm; imode++) {
        SuperkMode_t m = kSKReactionModes[isample][imode];
        if(m == kSKmUndefined) continue;
        bestfit_mode[imode] = new TH1D(*fFitInputs->Expectation().CurrPredictedErecoSpectrum(s, m));
        bestfit_mode[imode]->SetDirectory (0);
     }//imode

     // add spectra to the results array
     bestfit_tot->SetName(util::AsStringCompact(kSKsFCFV1RingMuLike,kSKmUndefined,"bestfit_").c_str());
     fResults->Add(bestfit_tot);
     for(int imode = 0; imode < knm; imode++) {
        SuperkMode_t m = kSKReactionModes[isample][imode];
        if(m == kSKmUndefined) continue;
        bestfit_mode[imode]->SetName(util::AsStringCompact(kSKsFCFV1RingMuLike,m,"bestfit_").c_str());
        fResults->Add(bestfit_mode[imode]);
     }//imode
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
     data -> SetName("data_1Rmu");
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

     chisq_surface -> SetName("chisq_surface_1Rmu");
     fResults -> Add(chisq_surface);
  }

  //
  // Save canvases with some standard plots
  //
  if(opt.Contains("Canvases")) {
    if(opt.Contains("BestFitSpectra") && opt.Contains("Observation")) {
        THStack * contrib_stack = new THStack("stack","spectrum contributions");
        for(int imode = 0; imode < knm; imode++) {
           SuperkMode_t m = kSKReactionModesOrderInStack[isample][imode]; 
           if(m == kSKmUndefined) continue;
           int jmode = util::FindModeSlot(s,m);
           contrib_stack->Add(bestfit_mode[jmode]);
        }
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
        spectrum_legend->AddEntry(data, "data", "P" );
        spectrum_legend->AddEntry(bestfit_tot,  util::AsStringRooTex(s,kSKmUndefined,"best fit: ").c_str(), "L" );
        for(int imode = 0; imode < knm; imode++) {
           SuperkMode_t m = kSKReactionModes[isample][imode];
           if(m == kSKmUndefined) continue;
           spectrum_legend->AddEntry(bestfit_mode[imode], util::AsStringRooTex(s,m,"best fit: ").c_str(), "LF" );
        }
        spectrum_legend->Draw();

        c1->SetName("canvas_1Rmu_1");
        c2->SetName("canvas_1Rmu_2");

        fResults->Add(c1);
        fResults->Add(c2);

   } // BestFitSpectra & Observation
  } // Canvases

}
//____________________________________________________________________________ 
void Superk1RmuErecoFitter::Save(const char * filename) 
{
  util::SetDefaultStyle();

  SuperkSample_t s = kSKsFCFV1RingMuLike;
  int isample = util::FindSampleSlot(s);

  cout << "Superk1RmuErecoFitter/Save >> Writing fit results in: " 
       << filename << endl;
  TFile * file_out = new TFile(filename, "recreate");

  TH1D * data = (TH1D*) fResults->FindObject("data_1Rmu");
  if(data) { 
    data->Write(); 
  }

  TH1D * bestfit_tot = (TH1D*) fResults->FindObject(util::AsStringCompact(kSKsFCFV1RingMuLike,kSKmUndefined,"bestfit_").c_str());
  if(bestfit_tot) {
   bestfit_tot->Write();
  }

  TH1D * bestfit_mode[knm];
  for(int imode = 0; imode < knm; imode++) {
     SuperkMode_t m = kSKReactionModes[isample][imode];
     if(m == kSKmUndefined) continue;
     bestfit_mode[imode] = (TH1D*) fResults->FindObject(util::AsStringCompact(kSKsFCFV1RingMuLike,m,"bestfit_").c_str());
     if(bestfit_mode[imode]) {
       bestfit_mode[imode]->Write();
     }
  }//imode

  TH2D * chisurf = (TH2D*) fResults->FindObject("chisq_surface_1Rmu");
  if(chisurf) {
    chisurf->Write();
  }

  TCanvas * c1 = (TCanvas*) fResults->FindObject("canvas_1Rmu_1");
  TCanvas * c2 = (TCanvas*) fResults->FindObject("canvas_1Rmu_2");
  if ( c1 ) { c1 -> Write(); }
  if ( c2 ) { c2 -> Write(); }

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



