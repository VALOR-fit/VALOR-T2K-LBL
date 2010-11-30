//____________________________________________________________________________
/*

*/
//____________________________________________________________________________

//#define __debug_level_1__
//#define __debug_level_2__
//#define __debug_level_3__

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
}
//____________________________________________________________________________
OscillationFitter::~OscillationFitter()
{

}
//____________________________________________________________________________
void OscillationFitter::Init(void)
{
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
void OscillationFitter::PrepareFitter(Option_t * option)
{
#ifdef __debug_level_1__
  cout << "OscillationFitter >> Fitter initialization..." << endl;
#endif

  if(!fFitInputs) {
     cerr << "OscillationFitter >> ** Null fit inputs!" << endl;
     return;
  }

  TString opt = option;

  if(opt.Contains("minuit")) {
    // get fitter
    TVirtualFitter::SetDefaultFitter("Minuit");
  } else {
    cerr << "OscillationFitter >> ** We have not tested any fitter other than MINUIT." << endl;
    cerr << "OscillationFitter >> ** Modify the code and proceed at your own risk." << endl;
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

  // calculate degrees of freedom

  fNonZeroBins = 0;
  const TAxis & ereco_binning = fFitInputs->Observation().ErecoBinning();
  for(int bin = 1; bin <= ereco_binning.GetNbins(); bin++) {
     double Etrue = ereco_binning.GetBinCenter(bin);
     double nev = fFitInputs->Observation().NumOfEvents(Etrue);
     if(nev>0) {
       fNonZeroBins++;
     }
  }
  fNFreeFitParams = fitter->GetNumberFreeParameters();
  fNdof           = fNonZeroBins - fNFreeFitParams;

#ifdef __debug_level_1__
  cout << "OscillationFitter >> ndof = " << fNdof
       << ", nfreeparams = " << fNFreeFitParams << endl;
#endif
}
//____________________________________________________________________________
void OscillationFitter::DoTheFit(void)
{
#ifdef __debug_level_1__
  cout << "OscillationFitter/DoTheFit >> Running the fitter..." << endl;
#endif

  TVirtualFitter * fitter = TVirtualFitter::GetFitter();

  // pass fitted data & expectation to MINUIT
  //OscillationFitter::FitInputs fit_inputs(fObs,fExp,fFitParams);
  gMinuit->SetObjectFit(fFitInputs);

  // MINUIT's minimization step
#ifdef __debug_level_1__
  cout << "OscillationFitter/DoTheFit >> Running MIGRAD..." << endl;
#endif
  Double_t arglist[100];
  arglist[0] = 5000;  // num of func calls
  arglist[1] = 0.01;  // tolerance
  fitter->ExecuteCommand("MIGRAD",arglist, 2);

  Double_t amin = 0;
  fitter->PrintResults(3,amin);

#ifdef __debug_level_1__
  cout << "OscillationFitter/DoTheFit >> Done!" << endl;
#endif
}
//____________________________________________________________________________
void OscillationFitter::Draw(
  TVirtualPad * c, Option_t * option, const char * outtag)
{
  util::SetDefaultStyle();

  if(!c) return;
  c->cd();
  c->SetBorderMode(0);
  c->SetFillColor (0);

  TString opt = option;

#ifdef __debug_level_1__
  cout << "OscillationFitter/Draw >> Drawing with option = " << opt.Data() << endl;
#endif

  //
  // draw fitter results
  //

  if(opt.Contains("Result")) { 

     TVirtualFitter * fitter = TVirtualFitter::GetFitter();

     // best fit results
     double sinsq2theta23_bestfit = fitter->GetParameter(0);
     double dm23sq_bestfit        = fitter->GetParameter(1);

     cout << "OscillationFitter/Draw >> Best fit results: " 
          << "sinsq2theta23 = " << sinsq2theta23_bestfit 
          << ", Dmsq23 = " << dm23sq_bestfit << endl;


     double xmin, ymin, xmax, ymax;

     c->Divide(2,2);

     //
     // plot unoscillated + oscillated spectra (sig + bkg) and best fit
     //

     c->cd(1);

     // get oscillated MC
     fFitInputs->FitParams().Update( conv::kPhypSinsq2Theta23, sinsq2theta23_bestfit );
     fFitInputs->FitParams().Update( conv::kPhypDm23sq,        dm23sq_bestfit        );

     TH1D * data = fFitInputs->Observation().AsTH1D();

     TH1D * bestfit_tot       = fFitInputs->Expectation().AsTH1D(fFitInputs->FitParams());
     TH1D * bestfit_numuccqe  = fFitInputs->Expectation().AsTH1D(fFitInputs->FitParams(), conv::kSKmNuMuCCQE);
     TH1D * bestfit_numucc1pi = fFitInputs->Expectation().AsTH1D(fFitInputs->FitParams(), conv::kSKmNuMuCC1pi);
     TH1D * bestfit_numuccoth = fFitInputs->Expectation().AsTH1D(fFitInputs->FitParams(), conv::kSKmNuMuCCoth);
     TH1D * bestfit_nc        = fFitInputs->Expectation().AsTH1D(fFitInputs->FitParams(), conv::kSKmNC);
     TH1D * bestfit_numubar   = fFitInputs->Expectation().AsTH1D(fFitInputs->FitParams(), conv::kSKmNuMuBar);
     TH1D * bestfit_nue       = fFitInputs->Expectation().AsTH1D(fFitInputs->FitParams(), conv::kSKmNuE);

     data -> SetMarkerStyle(20);
     data -> SetMarkerSize(1);

     bestfit_numucc1pi -> SetFillStyle ( 3354     );
     bestfit_numuccoth -> SetFillStyle ( 3345     );
     bestfit_nc        -> SetFillStyle ( 3006     );
     bestfit_numubar   -> SetFillStyle ( 3001     );
     bestfit_nue       -> SetFillStyle ( 3001     );
     bestfit_numucc1pi -> SetFillColor ( kBlue    );
     bestfit_numuccoth -> SetFillColor ( kMagenta );
     bestfit_nc        -> SetFillColor ( kRed     );
     bestfit_numubar   -> SetFillColor ( kYellow  );
     bestfit_nue       -> SetFillColor ( kGreen   );
     bestfit_numuccqe  -> SetLineColor ( kBlack   );
     bestfit_numucc1pi -> SetLineWidth ( 1        );
     bestfit_numuccoth -> SetLineWidth ( 1        );
     bestfit_nc        -> SetLineWidth ( 1        );
     bestfit_numuccqe  -> SetLineWidth ( 1        );
     bestfit_numubar   -> SetLineWidth ( 1        );
     bestfit_nue       -> SetLineWidth ( 1        );
     bestfit_tot       -> SetLineWidth ( 2        );
     bestfit_nc        -> SetLineStyle ( kDashed  );

     data -> GetXaxis()->SetTitle("E_{reco} (GeV)");
     data -> GetYaxis()->SetTitle("Number of events");
     data -> Draw("PERR");
     //data -> Draw("P");
     bestfit_tot   -> Draw("SAME");

     THStack * stack = new THStack("stack","spectrum contributions");
     stack->Add(bestfit_numuccqe);
     stack->Add(bestfit_numucc1pi);
     stack->Add(bestfit_numuccoth);
     stack->Add(bestfit_nc);
     stack->Add(bestfit_numubar);
     stack->Add(bestfit_nue);
     stack->Draw("same");

//     bestfit_ccqe  -> Draw("SAME");
//     bestfit_cc1pi -> Draw("SAME");
//     bestfit_ccoth -> Draw("SAME");
//     bestfit_nc    -> Draw("SAME");

     TLegend * spectrum_legend = new TLegend(0.6,0.5,0.9,0.9);
     spectrum_legend->SetFillColor (0);
     spectrum_legend->SetBorderSize(0);
     spectrum_legend->AddEntry(data,              "data",                               "P" );
     spectrum_legend->AddEntry(bestfit_tot,       "best fit",                           "L" );
     spectrum_legend->AddEntry(bestfit_numuccqe,  "best fit (#nu_{#mu} CCQE)",          "L" );
     spectrum_legend->AddEntry(bestfit_numucc1pi, "best fit (#nu_{#mu} CC1#pi)",        "LF");
     spectrum_legend->AddEntry(bestfit_numuccoth, "best fit (#nu_{#mu} CC other)",      "LF");
     spectrum_legend->AddEntry(bestfit_nc,        "best fit (#nu_{#mu}/#nu_{#tau} NC)", "LF");
     spectrum_legend->AddEntry(bestfit_numubar,   "best fit (#bar{#nu_{#mu}})",         "LF");
     spectrum_legend->AddEntry(bestfit_nue,       "best fit (#nu_{e})", "LF");
     spectrum_legend->Draw();

     c->Update();

     //
     // plot contour
     //

     c->cd(2);
     
     // chisquare surface

     double sinsq2theta23_min  = fFitInputs->FitParams().MinValue(conv::kPhypSinsq2Theta23);
     double sinsq2theta23_max  = fFitInputs->FitParams().MaxValue(conv::kPhypSinsq2Theta23);
     double sinsq2theta23_step = fFitInputs->FitParams().StepSize(conv::kPhypSinsq2Theta23);
     double dm23sq_min         = fFitInputs->FitParams().MinValue(conv::kPhypDm23sq);
     double dm23sq_max         = fFitInputs->FitParams().MaxValue(conv::kPhypDm23sq);
     double dm23sq_step        = fFitInputs->FitParams().StepSize(conv::kPhypDm23sq);

//     int sinsq2theta23_nbins   = (int) ((sinsq2theta23_max - sinsq2theta23_min) / sinsq2theta23_step); 
//     int dm23sq_nbins          = (int) ((dm23sq_max        - dm23sq_min       ) / dm23sq_step       );
     int sinsq2theta23_nbins   = 40;
     int dm23sq_nbins          = 40;
  
     TH2D * delta_chisq_surf = new TH2D(
        "delta_chisq_surf","#delta{#chi^{2}} surface",
        sinsq2theta23_nbins, sinsq2theta23_min, sinsq2theta23_max, dm23sq_nbins, dm23sq_min, dm23sq_max);
  
     double chisq_min = 999999;
     for(int i = 1; i <= delta_chisq_surf->GetNbinsX(); i++) {
        for(int j = 1; j <= delta_chisq_surf->GetNbinsY(); j++) {

           fFitInputs->FitParams().Update( conv::kPhypSinsq2Theta23, delta_chisq_surf->GetXaxis()->GetBinCenter(i) );
           fFitInputs->FitParams().Update( conv::kPhypDm23sq,        delta_chisq_surf->GetYaxis()->GetBinCenter(j) );
           double chisq = fFitInputs->Chisquare();
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

     const int nlev = 3;
     double levels[nlev] = { 
         util::DeltaChisquare(conv::kCL68, fNFreeFitParams),
         util::DeltaChisquare(conv::kCL90, fNFreeFitParams),
         util::DeltaChisquare(conv::kCL99, fNFreeFitParams)
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
     double nobs = this->Observation().NumOfEvents(Ereco);
     double nexp = this->Expectation().Generate(Ereco, this->FitParams());
     if      (opt.Contains("poisson")) { delta = util::ChisqPoisson(nobs,nexp); }
     else if (opt.Contains("gauss")  ) { delta = util::ChisqGauss  (nobs,nexp); }
     chisq += delta;

#ifdef __debug_level_3__
     cout << "OscillationFitter> chisq (bin = " << ibin << " / Ereco = " 
          << Ereco << " GeV) = " << delta << endl;
#endif

  }//ibin  

#ifdef __debug_level_2__
  cout << "OscillationFitter> chisq = " << chisq << endl;
#endif

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

