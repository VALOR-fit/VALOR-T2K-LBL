//____________________________________________________________________________
/*

*/
//____________________________________________________________________________

#define _debug_

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <sstream>

#include <TMath.h>
#include <TVirtualFitter.h>
#include <TGraph.h>
#include <TH2D.h>
#include <TMarker.h>
#include <TLegend.h>
#include <TString.h>
#include <TGAxis.h>
#include <TCanvas.h>

#include "year1_numu_osc_analysis/OscillationFitter.h"
#include "year1_numu_osc_analysis/SuperkObservation.h"
#include "year1_numu_osc_analysis/SuperkExpectation.h"
#include "year1_numu_osc_analysis/PhysParamList.h"
#include "year1_numu_osc_analysis/Utilities.h"

using std::cerr;
using std::cout;
using std::endl;
using std::ostringstream;

using namespace t2kana;
using namespace t2kana::yr1numu;
   
//____________________________________________________________________________
OscillationFitter::OscillationFitter()
{
  this->Init();
}
//____________________________________________________________________________
OscillationFitter::OscillationFitter(
   FitInputs * fit_inputs, Option_t * option)
{
  this->Init();
  this->SetInputs(fit_inputs);
  this->PrepareFitter(option);
//  this->SetObservation(skobs);
//  this->SetExpectation(skexp);
//  this->SetFitParams(params);
}
//____________________________________________________________________________
OscillationFitter::~OscillationFitter()
{

}
//____________________________________________________________________________
void OscillationFitter::Init(void)
{
//  fObs              = 0;
//  fExp              = 0;
//  fFitParams        = 0;

  fFitInputs        = 0;

  fKnowActualParams = false;
  fSinsq2Theta23    = -1;
  fDm23sq           = -1;

  fNdof             = 0;
  fNFreeFitParams   = 0;
  fNonZeroBins      = 0;
}
//____________________________________________________________________________
void OscillationFitter::ActualOscParams (double sinsq2theta23, double dm23sq)
{
  fKnowActualParams = true;
  fSinsq2Theta23    = sinsq2theta23;
  fDm23sq           = dm23sq;
}
//____________________________________________________________________________
void OscillationFitter::SetInputs(FitInputs * fit_inputs)
{
  fFitInputs = fit_inputs;
}
//____________________________________________________________________________
/*
void OscillationFitter::SetObservation(SuperkObservation * skobs)
{
  cout << "OscillationFitter> Setting SK observation..." << endl;
  assert(skobs);

  fObs = skobs;

  fNonZeroBins = 0;

  assert( fData->ErecoCCQElike() );

  int nbins = fData->ErecoCCQElike()->GetNbinsX();
  for(int i = 1; i <= nbins; i++) {
    double y = fData->ErecoCCQElike()->GetBinContent(i);
	cout << i << " --> " << y << endl;
    if(y>0) {
      fNonZeroBins++;
    }
  }
}
//____________________________________________________________________________
void OscillationFitter::SetExpectation(SuperkExpectation * skexp)
{
  cout << "OscillationFitter> Setting SK expectation..." << endl;
  assert(skexp);

  fExp = skexp;
}
*/
//____________________________________________________________________________
void OscillationFitter::PrepareFitter(Option_t * option)
{
  cout << "OscillationFitter> Fitter initialization..." << endl;

  if(!fFitInputs) {
     cout << "OscillationFitter> Set the fit inputs first!" << endl;
     return;
  }

  TString opt = option;

  if(opt.Contains("minuit")) {
    // get fitter
    TVirtualFitter::SetDefaultFitter("Minuit");
  } else {
    cerr << "** We have not tested any fitter other than MINUIT." << endl;
    cerr << "** Modify the code and proceed at your own risk." << endl;
    exit(1);
  }

  PhysParamList & fit_params = fFitInputs->FitParams();
  TVirtualFitter * fitter = TVirtualFitter::Fitter(0,fit_params.Size());

  // set params
  fit_params.SetTo(fitter);

  // set minimization function
  fitter->SetFCN(FitFunc1);

  // set print level
  Double_t arglist[100];
  arglist[0] = 0;
  fitter->ExecuteCommand("SET PRINT",arglist,2);
}
//____________________________________________________________________________
void OscillationFitter::DoTheFit(void)
{
  cout << "OscillationFitter> Running the fitter..." << endl;

  TVirtualFitter * fitter = TVirtualFitter::GetFitter();

  fNdof           = fNonZeroBins - fitter->GetNumberFreeParameters();
  fNFreeFitParams = fitter->GetNumberFreeParameters();

  cout << "OscillationFitter> ndof = " << fNdof
       << ", nfreeparams = " << fNFreeFitParams << endl;

  // pass fitted data & expectation to MINUIT
  //OscillationFitter::FitInputs fit_inputs(fObs,fExp,fFitParams);
  gMinuit->SetObjectFit(fFitInputs);

  // MINUIT's minimization step
  cout << "OscillationFitter> Running MIGRAD..." << endl;
  Double_t arglist[100];
  arglist[0] = 5000;  // num of func calls
  arglist[1] = 0.01;  // tolerance
  fitter->ExecuteCommand("MIGRAD",arglist, 2);

  Double_t amin = 0;
  fitter->PrintResults(3,amin);

  cout << "OscillationFitter> Done!" << endl;
}
//____________________________________________________________________________
void OscillationFitter::Draw(
  TVirtualPad * c, Option_t * option, const char * outtag)
{
/*
  style::SetDefaultStyle();

  if(!c) return;
  c->cd();
  c->SetBorderMode(0);
  c->SetFillColor (0);

  TString opt = option;

  //
  // draw fitter inputs
  //

  if(opt.Contains("fit:inputs")) { 
     c->Divide(2,2);
     fMC->Draw(c->GetPad(1),  "Etrue:def");
     if(fKnowActualParams) {
       double xmax = fMC->ErecoCCQElike()->GetXaxis()->GetXmax();
       double ymax = 1.2 * fMC->ErecoCCQElike()->GetBinContent(fMC->ErecoCCQElike()->GetMaximumBin());
       TGaxis *axis = new TGaxis(xmax, 0, xmax, ymax, 0, 1, 50510, "-L");
       axis->SetLabelColor(kRed);
       axis->SetLabelOffset(0.030);
       axis->SetLabelSize(0.025);
       axis->Draw();
       // scale probability (Pmax := 1 --> scale) for drawing purposes
       // only (so that the oscillation probability can be superimposed
       // on the spectrum plot - both vertical axis scales will be shown)
       double scale = ymax;
       const int n = 2000;
       double e[n], p[n];
       double de = (sksel::kEnu_max-sksel::kEnu_min)/(n-1);
       for(int i=0; i<n; i++) {
          e[i] = sksel::kEnu_min + i*de;
          p[i] = utils::osc::NuMuSurvProb2D(fSinsq2Theta23, fDm23sq, e[i]);
          p[i] *= scale; // for plotting on canvas with two superimposed y-axes
       }
       TGraph * grosc = new TGraph(n,e,p);
       grosc->SetLineColor(kRed);
       grosc->Draw("lsame");
     }
     fMC->Draw(c->GetPad(2),  "Ereco:def");
     fMC->Draw(c->GetPad(3),  "Ereco2Etrue");

  }//option = fit:inputs

  //
  // draw fitter results
  //

  if(opt.Contains("fit:result")) { 

     TVirtualFitter * fitter = TVirtualFitter::GetFitter();

     // best fit results
     double sinsq2theta23_bestfit = fitter->GetParameter(0);
     double dm23sq_bestfit        = fitter->GetParameter(1);
     double norm_bestfit          = fitter->GetParameter(2);

     double xmin, ymin, xmax, ymax;

     c->Divide(2,2);

     //
     // plot unoscillated + oscillated spectra (sig + bkg) and best fit
     //

     c->cd(1);

     // get oscillated MC
     fFitParams->Update( kPhypSinsq2Theta23, sinsq2theta23_bestfit );
     fFitParams->Update( kPhypDm23sq,        dm23sq_bestfit        );
     fFitParams->Update( kPhypNorm,          norm_bestfit          );
     fMC->Update(fFitParams);

     fData -> ErecoCCQElike() -> SetMarkerStyle(20);
     fData -> ErecoCCQElike() -> SetMarkerSize(0.8);

     TH1D * bestfit_tot   = utils::misc::Copy( fMC -> CurrErecoCCQElike() );
     TH1D * bestfit_ccqe  = utils::misc::Copy( fMC -> CurrErecoSigCCQE()  ); 
     TH1D * bestfit_bkgcc = utils::misc::Copy( fMC -> CurrErecoBkgCC()    ); 
     TH1D * bestfit_bkgnc = utils::misc::Copy( fMC -> CurrErecoBkgNC()    ); 

     bestfit_bkgcc -> SetFillStyle ( 3006    );
     bestfit_bkgnc -> SetFillStyle ( 3006    );
     bestfit_bkgcc -> SetFillColor ( kBlue   );
     bestfit_bkgnc -> SetFillColor ( kRed    );
     bestfit_ccqe  -> SetLineColor ( kBlack  );
     bestfit_ccqe  -> SetLineWidth ( 1       );
     bestfit_tot   -> SetLineWidth ( 2       );
   
     fData->ErecoCCQElike()->Draw("PERR");
     bestfit_tot   -> Draw("SAME");
     bestfit_ccqe  -> Draw("SAME");
     bestfit_bkgcc -> Draw("SAME");
     bestfit_bkgnc -> Draw("SAME");

     TLegend * spectrum_legend = new TLegend(0.6,0.5,0.9,0.9);
     spectrum_legend->SetFillColor (0);
     spectrum_legend->SetBorderSize(0);
     spectrum_legend->AddEntry(fData,         "mock data",         "P" );
     spectrum_legend->AddEntry(bestfit_tot,   "best fit",          "L" );
     spectrum_legend->AddEntry(bestfit_ccqe,  "best fit (ccqe)",   "L" );
     spectrum_legend->AddEntry(bestfit_bkgcc, "best fit (bkg cc)", "LF");
     spectrum_legend->AddEntry(bestfit_bkgnc, "best fit (bkg nc)", "LF");
     spectrum_legend->Draw();

     fData->ErecoCCQElike()->GetXaxis()->SetTitle("E_{reco} (GeV)");
     fData->ErecoCCQElike()->GetYaxis()->SetTitle("number of events");

     //
     // plot contour
     //

     c->cd(2);
     
     // chisquare surface
  
     int nDm23sqBins       = utils::misc::NBins(osc23::kDm23sqMin,        osc23::kDm23sqMax,        osc23::kDm23sqStep       );
     int nSinsqTheta23Bins = utils::misc::NBins(osc23::kSinsq2Theta23Min, osc23::kSinsq2Theta23Max, osc23::kSinsq2Theta23Step);
  
     TH2D * delta_chisq_surf = new TH2D("delta_chisq_surf","#delta{#chi^{2}} surface",
        nSinsqTheta23Bins, osc23::kSinsq2Theta23Min, osc23::kSinsq2Theta23Max,
        nDm23sqBins,       osc23::kDm23sqMin,        osc23::kDm23sqMax);
  
     double chisq_min = 999999;
     for(int i = 1; i <= delta_chisq_surf->GetNbinsX(); i++) {
        for(int j = 1; j <= delta_chisq_surf->GetNbinsY(); j++) {

           fFitParams->Update( kPhypSinsq2Theta23, delta_chisq_surf->GetXaxis()->GetBinCenter(i) );
           fFitParams->Update( kPhypDm23sq,        delta_chisq_surf->GetYaxis()->GetBinCenter(j) );
           fFitParams->Update( kPhypNorm,          norm_bestfit );
           fMC->Update(fFitParams);
           double chisq = fMC->Chisquare(*fData);
           chisq_min = TMath::Min(chisq_min, chisq);
       
           delta_chisq_surf->Fill( 
               delta_chisq_surf->GetXaxis()->GetBinCenter(i),
               delta_chisq_surf->GetYaxis()->GetBinCenter(j),
               chisq);
        }//j
     }//i

     for(int i = 1; i <= delta_chisq_surf->GetNbinsX(); i++) {
        for(int j = 1; j <= delta_chisq_surf->GetNbinsY(); j++) {
            double chisq = delta_chisq_surf->GetBinContent(i,j);
            delta_chisq_surf->SetBinContent(i,j, TMath::Max(0., chisq-chisq_min));
        }//j
     }//i

     TGaxis::SetMaxDigits(2);
  
     delta_chisq_surf->GetXaxis()->SetTitle("sin^{2}(2#theta_{23})");
     delta_chisq_surf->GetYaxis()->SetTitle("|#Deltam_{32}^{2}| (eV^{2}/c^{4})");
     delta_chisq_surf->GetYaxis()->SetTitleOffset(1);
     delta_chisq_surf->GetYaxis()->SetTitleSize(0.05);
     delta_chisq_surf->GetXaxis()->SetTitleSize(0.05);

     const int nlev = 2;
     double levels[nlev] = { 
         T2KConfLevel::DeltaChisquare(kCL90, fNFreeFitParams),
         T2KConfLevel::DeltaChisquare(kCL99, fNFreeFitParams)
     };

     delta_chisq_surf->SetContour(nlev,levels);
     delta_chisq_surf->Draw("cont2");

     TMarker * best_fit_point =
          new TMarker(sinsq2theta23_bestfit, dm23sq_bestfit, 29);
     best_fit_point->SetMarkerSize(1.4);
     best_fit_point->Draw();

     if(fKnowActualParams) {
        TMarker * true_point =
            new TMarker(fSinsq2Theta23, fDm23sq, 20);
        true_point->SetMarkerSize(1.2);
        true_point->SetMarkerColor(kRed);
        true_point->Draw();
      }
     
      c->GetPad(2)->Update();

      //
      // plot chisq = f (sinsq2theta23)
      //

      c->cd(3);
      gMinuit->Command("SCAn 1");
      TGraph * gr1 = (TGraph *) gMinuit->GetPlot();
      for(int i=0; i < gr1->GetN(); i++) {
          double x,y;
          gr1->GetPoint(i,x,y);
          gr1->SetPoint(i,x,y/fNdof);
      }
      xmin  = ( gr1->GetX() )[TMath::LocMin(gr1->GetN(),gr1->GetX())];
      xmax  = ( gr1->GetX() )[TMath::LocMax(gr1->GetN(),gr1->GetX())];
      ymin  = ( gr1->GetY() )[TMath::LocMin(gr1->GetN(),gr1->GetY())];
      ymax  = ( gr1->GetY() )[TMath::LocMax(gr1->GetN(),gr1->GetY())];
      TH1F* hframe1 = c->GetPad(3)->DrawFrame(xmin,0.8*ymin,xmax,1.2*ymax);
      hframe1->GetXaxis()->SetTitle("sin^{2}(2#theta_{23})");
      hframe1->GetYaxis()->SetTitle("#chi^{2}");
      gr1->Draw("l");

      //
      // plot chisq = f (Dm23sq)
      //

      c->cd(4);
      gMinuit->Command("SCAn 2");
      TGraph * gr2 = (TGraph *) gMinuit->GetPlot();
      for(int i=0; i < gr2->GetN(); i++) {
          double x,y;
          gr2->GetPoint(i,x,y); 
          gr2->SetPoint(i,x,y/fNdof);
      }
      xmin  = ( gr2->GetX() )[TMath::LocMin(gr2->GetN(),gr2->GetX())];
      xmax  = ( gr2->GetX() )[TMath::LocMax(gr2->GetN(),gr2->GetX())];
      ymin  = ( gr2->GetY() )[TMath::LocMin(gr2->GetN(),gr2->GetY())];
      ymax  = ( gr2->GetY() )[TMath::LocMax(gr2->GetN(),gr2->GetY())];
      TH1F* hframe2 = c->GetPad(4)->DrawFrame(xmin,0.8*ymin,xmax,1.2*ymax);
      hframe2->GetXaxis()->SetTitle("|#Deltam_{32}^{2}| (eV^{2}/c^{4})");
      hframe2->GetYaxis()->SetTitle("#chi^{2}");
      gr2->Draw("l");
  }//option = fit:result

  // Save canvas
  c->Update();
  if(outtag) {
    ostringstream outgif;
    ostringstream outroot;
    outgif  << outtag << ".gif";
    outroot << outtag << ".root";
    c->SaveAs( outgif.str().c_str()  );
    c->SaveAs( outroot.str().c_str() );
  }
*/
}
//____________________________________________________________________________
// FitInputs utility class 
//
double OscillationFitter::FitInputs::Chisquare(Option_t * option)
{
  double chisq = 0;
  TString opt = option;

  const TAxis & ereco_binning =
          this->Observation().ErecoBinning();

  int nbins = ereco_binning.GetNbins(); 
  for(int ibin = 1; ibin <= nbins; ibin++) {
     double Ereco = ereco_binning.GetBinCenter(ibin);
     double delta = 0.;
     double nobs = this->Observation().NumOfEvents();
     double nexp = this->Expectation().Generate(Ereco, this->FitParams());
     if      (opt.Contains("poisson")) { delta = util::ChisqPoisson(nobs,nexp); }
     else if (opt.Contains("gauss")  ) { delta = util::ChisqGauss  (nobs,nexp); }
     chisq += delta;
  }//ibin  

  return chisq;
}
//____________________________________________________________________________
// MINUIT fit functions
//
void t2kana::yr1numu::FitFunc1 (
  Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t /*iflag*/)
{
  OscillationFitter::FitInputs * inputs = 
      (OscillationFitter::FitInputs *) gMinuit->GetObjectFit();

//SuperkObservation & skobs = inputs->Observation(); // data
//SuperkExpectation & skexp = inputs->Expectation(); // model
  PhysParamList &     fpar  = inputs->FitParams();   // params

  // synchronize PhysParamList with the paramerer array used by MINUIT        
  fpar.Sync(par); 

  // judge the agreement between the data and the current prediction
  f = inputs->Chisquare("poisson");
}
//____________________________________________________________________________

