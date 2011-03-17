//____________________________________________________________________________ 
/*
 Changes since v1.0.0:

 @ 16/03/2011 - CA
   Removed the numu10a namespace. Renamed: OscillationFitter->NuMuDisappFitter.
   The FitInputs classs is no longer nested in the fitter class.

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
#include <TGaxis.h> 
#include <TCanvas.h> 
#include <TFile.h>

#include "NuMuDisappFitter.h" 
#include "SuperkObservation.h" 
#include "SuperkExpectation.h" 
#include "FitInputs.h" 
#include "PhysParamList.h" 
#include "Utilities.h"
#include "Conventions.h"

using std::cerr; 
using std::cout; 
using std::endl; 
using std::ostringstream;

using namespace t2kana; 
using namespace t2kana::conv;
   
//____________________________________________________________________________ 
NuMuDisappFitter::NuMuDisappFitter() 
{
  this->Init();
}
//____________________________________________________________________________ 
NuMuDisappFitter::NuMuDisappFitter(
  FitInputs * fit_inputs, Option_t * option) 
{
  this->Init();
  this->SetInputs(fit_inputs);
  this->PrepareFitter(option);
}
//____________________________________________________________________________ 
NuMuDisappFitter::~NuMuDisappFitter() 
{

}
//____________________________________________________________________________ 
void NuMuDisappFitter::Init(void) 
{
  fFitInputs = 0;

  fKnowActualParams = false;
  fSinsq2Theta23 = -1;
  fDm23sq = -1;

  fNdof           = 0;
  fNFreeFitParams = 0;
  fNFitBins       = 0;

  fChisqSurface = 0;
}
//____________________________________________________________________________ 
void  NuMuDisappFitter::ActualOscParams 
(double sinsq2theta23, double dm23sq) {
  fKnowActualParams = true;
  fSinsq2Theta23 = sinsq2theta23;
  fDm23sq = dm23sq;
}
//____________________________________________________________________________ 
void NuMuDisappFitter::SetInputs(FitInputs * fit_inputs) 
{
  fFitInputs = fit_inputs;
}
//____________________________________________________________________________ 
void  NuMuDisappFitter::PrepareFitter(Option_t * option) 
{
#ifdef __debug_level_1__
  cout << "NuMuDisappFitter >> Fitter initialization..." << endl;
#endif

  if(!fFitInputs) {
     cerr << "NuMuDisappFitter >> ** Null fit inputs!" << endl;
     return;
  }

  fOption = option;

  PhysParamList & fit_params = fFitInputs->FitParams();

  if(fOption.Contains("minuit")) {
    // get fitter
    TVirtualFitter::SetDefaultFitter("Minuit");
    TVirtualFitter * fitter = TVirtualFitter::Fitter(0,fit_params.Size());

   // set params
   fit_params.SetTo(fitter);

   // set minimization function
   fitter->SetFCN(FitFunc);

   // set print level
   Double_t arglist[100];
   arglist[0] = 0;
   fitter->ExecuteCommand("SET PRINT",arglist,2);
  }

  // calculate degrees of freedom

  fNFitBins = kNvErecoBins;
  fNFreeFitParams = fit_params.Size(); //fitter->GetNumberFreeParameters();
  fNdof = fNFitBins - fNFreeFitParams;

  cout << "NuMuDisappFitter >> ndof = " << fNdof
       << ", nfreeparams = " << fNFreeFitParams << endl;
}
//____________________________________________________________________________ 
bool NuMuDisappFitter::DoTheFit(void) 
{
  if(fOption.Contains("minuit")) {
    return this->DoTheFitMinuit();
  }
  else
  if(fOption.Contains("custom")) {
    return this->DoTheFitCustom();
  }
  return false;
}
//____________________________________________________________________________ 
bool NuMuDisappFitter::DoTheFitMinuit(void) 
{
#ifdef __debug_level_1__
  cout << "NuMuDisappFitter/DoTheFit >> Running the fitter..." << endl;
#endif

  TVirtualFitter * fitter = TVirtualFitter::GetFitter();

  // pass fitted data & expectation to MINUIT
  gMinuit->SetObjectFit(fFitInputs);

  // MINUIT's minimization step
#ifdef __debug_level_1__
  cout << "NuMuDisappFitter/DoTheFit >> Running MIGRAD..." << endl;
#endif
  Double_t arglist[100];
  arglist[0] = 5000; // num of func calls
  arglist[1] = 0.01; // tolerance
  fitter->ExecuteCommand("MIGRAD",arglist, 2);

  Double_t amin = 0;
  fitter->PrintResults(3,amin);

  // Dummy variables - need to feed them to the function
  // used for getting error matrix status
  double ha(0.0);
  double edm, errdef;
  int nvpar, nparx;
  int errMatrixStatus = fitter->GetStats(ha, edm, errdef, nvpar, nparx);

#ifdef __debug_level_2__
  cout << "Final error matrix status = " << errMatrixStatus << endl;
#endif

#ifdef __debug_level_1__
  cout << "NuMuDisappFitter/DoTheFit >> Done!" << endl;
#endif

   return (errMatrixStatus == 3);
}
//____________________________________________________________________________ 
bool NuMuDisappFitter::DoTheFitCustom(void) 
{
#ifdef __debug_level_1__
  cout << "NuMuDisappFitter/DoTheFit >> Running the fitter..." << endl;
#endif

  double sinsq2theta23_bestfit = 0;
  double dm23sq_bestfit = 0;

  this->BuildChisqSurface();
  
  double chisq_min = 999999;
  for(int i = 1; i <= fChisqSurface->GetNbinsX(); i++) {
    for(int j = 1; j <= fChisqSurface->GetNbinsY(); j++) {
         double chisq = fChisqSurface->GetBinContent(i,j);
         if(chisq < chisq_min) {
	    chisq_min = chisq;
            sinsq2theta23_bestfit = fChisqSurface->GetXaxis()->GetBinCenter(i);
            dm23sq_bestfit = fChisqSurface->GetYaxis()->GetBinCenter(j);
         }       
     }//j
   }//i

  fFitInputs->FitParams().Update( conv::kPhypSinsq2Theta23, sinsq2theta23_bestfit );
  fFitInputs->FitParams().Update( conv::kPhypDm23sq, dm23sq_bestfit );

  cout << "NuMuDisappFitter >> Best fit params: " << endl;
  cout << "NuMuDisappFitter >> Dmsq = " << dm23sq_bestfit << endl;
  cout << "NuMuDisappFitter >> sinsq2theta = " << sinsq2theta23_bestfit << endl;
  cout << "NuMuDisappFitter >> chisqmin/ndof = " << chisq_min/fNdof << endl;

  return true;
}
//____________________________________________________________________________ 
void  NuMuDisappFitter::BuildChisqSurface(void) 
{
   if(fChisqSurface) delete fChisqSurface;

   double sinsq2theta23_save = fFitInputs->FitParams().CurrValue( conv::kPhypSinsq2Theta23 );
   double dm23sq_save = fFitInputs->FitParams().CurrValue( conv::kPhypDm23sq );

   double sinsq2theta23_min = fFitInputs->FitParams().MinValue(conv::kPhypSinsq2Theta23);
   double sinsq2theta23_max = fFitInputs->FitParams().MaxValue(conv::kPhypSinsq2Theta23);
   double sinsq2theta23_step = fFitInputs->FitParams().StepSize(conv::kPhypSinsq2Theta23);
   double dm23sq_min = fFitInputs->FitParams().MinValue(conv::kPhypDm23sq);
   double dm23sq_max = fFitInputs->FitParams().MaxValue(conv::kPhypDm23sq);
   double dm23sq_step = fFitInputs->FitParams().StepSize(conv::kPhypDm23sq);

   // make sure that min and max are at the bin centre
   int sinsq2theta23_nbins = 1 + (int) ((sinsq2theta23_max - sinsq2theta23_min) / sinsq2theta23_step);
   int dm23sq_nbins = 1 + (int) ((dm23sq_max - dm23sq_min ) / dm23sq_step );

   fChisqSurface = new TH2D(
      "fChisqSurface","#delta{#chi^{2}} surface",
         sinsq2theta23_nbins, sinsq2theta23_min - 0.5*sinsq2theta23_step, sinsq2theta23_max + 0.5*sinsq2theta23_step,
         dm23sq_nbins, dm23sq_min - 0.5*dm23sq_step, dm23sq_max + 0.5*dm23sq_step);
   fChisqSurface->SetDirectory(0);
   fChisqSurface->GetXaxis()->SetRangeUser(0.,1.);

   for(int i = 1; i <= fChisqSurface->GetNbinsX(); i++) {
     for(int j = 1; j <= fChisqSurface->GetNbinsY(); j++) {
         fFitInputs->FitParams().Update( conv::kPhypSinsq2Theta23, fChisqSurface->GetXaxis()->GetBinCenter(i) );
         fFitInputs->FitParams().Update( conv::kPhypDm23sq, fChisqSurface->GetYaxis()->GetBinCenter(j) );
         double chisq = fFitInputs->Chisquare();
         fChisqSurface->SetBinContent(i,j,chisq);
     }//j
   }//i

   fChisqSurface->GetXaxis()->SetTitle("sin^{2}(2#theta_{23})");
   fChisqSurface->GetYaxis()->SetTitle("|#Deltam_{32}^{2}| (eV^{2}/c^{4})");
   fChisqSurface->GetYaxis()->SetTitleOffset(1);
   fChisqSurface->GetYaxis()->SetTitleSize(0.05);
   fChisqSurface->GetXaxis()->SetTitleSize(0.05);

  // restore inputs osc params
  fFitInputs->FitParams().Update( conv::kPhypSinsq2Theta23, sinsq2theta23_save );
  fFitInputs->FitParams().Update( conv::kPhypDm23sq, dm23sq_save );
}
//____________________________________________________________________________ 
TH1D* NuMuDisappFitter::MarginalizeChisqSurface(Option_t * option) 
{
  TString opt = option;

  TH1D * marginalized_chisq = 0;

  int    sinsq2theta23_nbins = fChisqSurface->GetXaxis()->GetNbins();
  double sinsq2theta23_min   = fChisqSurface->GetXaxis()->GetXmin();
  double sinsq2theta23_max   = fChisqSurface->GetXaxis()->GetXmax();
  int    dm23sq_nbins        = fChisqSurface->GetYaxis()->GetNbins();
  double dm23sq_min          = fChisqSurface->GetYaxis()->GetXmin();
  double dm23sq_max          = fChisqSurface->GetYaxis()->GetXmax();

  //
  // pick minimum chisq along marginalised parameter
  //
  if(opt.Contains("marginalize_dm23sq")) 
  {
     //
     // integrate out Dm23sq
     //
     marginalized_chisq = new TH1D("","", sinsq2theta23_nbins, sinsq2theta23_min, sinsq2theta23_max);
     for(int ibin = 1; ibin <= sinsq2theta23_nbins; ibin++) {
        double chisq = 999999999;
	for(int jbin = 1; jbin <= dm23sq_nbins; jbin++) {
	    chisq = TMath::Min(chisq, fChisqSurface->GetBinContent(ibin,jbin));
        }//jbin
	marginalized_chisq->SetBinContent(ibin, chisq);
     }//ibin
  } 
  else
  if(opt.Contains("marginalize_sinsq2theta23")) 
  {
     //
     // integrate out sinsq2theta23
     //
     marginalized_chisq = new TH1D("","", dm23sq_nbins, dm23sq_min, dm23sq_max);
     for(int jbin = 1; jbin <= dm23sq_nbins; jbin++) {
        double chisq = 999999999;
	for(int ibin = 1; ibin <= sinsq2theta23_nbins; ibin++) {
	    chisq = TMath::Min(chisq, fChisqSurface->GetBinContent(ibin,jbin));
        }//ibin
	marginalized_chisq->SetBinContent(jbin, chisq);
     }//jbin
  } 

/*
  //
  // takes average chisq over the marginalized parameter
  //
  if(opt.Contains("marginalize_dm23sq")) 
  {
     //
     // integrate out Dm23sq
     //
     marginalized_chisq = new TH1D("","", sinsq2theta23_nbins, sinsq2theta23_min, sinsq2theta23_max);
     for(int ibin = 1; ibin <= sinsq2theta23_nbins; ibin++) {
        double chisq = 0;
	for(int jbin = 1; jbin <= dm23sq_nbins; jbin++) {
	    chisq += fChisqSurface->GetBinContent(ibin,jbin);
        }//jbin
        chisq /= dm23sq_nbins;
	marginalized_chisq->SetBinContent(ibin, chisq);
     }//ibin
  } 
  else
  if(opt.Contains("marginalize_sinsq2theta23")) 
  {
     //
     // integrate out sinsq2theta23
     //
     marginalized_chisq = new TH1D("","", dm23sq_nbins, dm23sq_min, dm23sq_max);
     for(int jbin = 1; jbin <= dm23sq_nbins; jbin++) {
        double chisq = 0;
	for(int ibin = 1; ibin <= sinsq2theta23_nbins; ibin++) {
	    chisq += fChisqSurface->GetBinContent(ibin,jbin);
        }//ibin
        chisq /= sinsq2theta23_nbins;
	marginalized_chisq->SetBinContent(jbin, chisq);
     }//jbin
  } 
*/

  if(opt.Contains("subtract_chisqmin")) {
    double chisqmin = marginalized_chisq->GetMinimum();
     for(int ibin = 1; ibin <= marginalized_chisq->GetNbinsX(); ibin++) {
       double chisq = marginalized_chisq->GetBinContent(ibin);
       chisq -= chisqmin;
       chisq = TMath::Max(0., chisq);
       marginalized_chisq->SetBinContent(ibin, chisq);
     }	
  }

  return marginalized_chisq;
}
//____________________________________________________________________________ 
void NuMuDisappFitter::Draw(
  TVirtualPad * c, Option_t * option, const char * outtag) 
{
  util::SetDefaultStyle();

  if(!c) return;
  c->cd();
  c->SetBorderMode(0);
  c->SetFillColor (0);

  TString opt = option;


  TFile * file_out = 0;
  if(outtag) {
     ostringstream outroot;
     outroot << outtag << ".root";
     file_out = new TFile(outroot.str().c_str(), "recreate");
  }


#ifdef __debug_level_1__
  cout << "NuMuDisappFitter/Draw >> Drawing with option = " << opt.Data() << endl;
#endif

  //
  // draw fitter results
  //

  if(opt.Contains("StandardOscFitPlots")) {

     // best fit results
     double sinsq2theta23_bestfit = fFitInputs->FitParams().CurrValue( conv::kPhypSinsq2Theta23 );
     double dm23sq_bestfit = fFitInputs->FitParams().CurrValue( conv::kPhypDm23sq );

     cout << "NuMuDisappFitter/Draw >> Best fit results: "
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
     fFitInputs->FitParams().Update( conv::kPhypDm23sq, dm23sq_bestfit               );
     fFitInputs->Expectation().Generate(fFitInputs->FitParams());

     TH1D * bestfit_tot            = new TH1D(*fFitInputs->Expectation().CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike));
     TH1D * bestfit_numuccqe       = new TH1D(*fFitInputs->Expectation().CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuMuCCQE));
     TH1D * bestfit_numucc1pi      = new TH1D(*fFitInputs->Expectation().CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuMuCC1pi));
     TH1D * bestfit_numuccoth      = new TH1D(*fFitInputs->Expectation().CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuMuCCoth));
     TH1D * bestfit_nc             = new TH1D(*fFitInputs->Expectation().CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNC));
     TH1D * bestfit_numubarccqe    = new TH1D(*fFitInputs->Expectation().CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuMuBarCCQE));
     TH1D * bestfit_numubarcc1pi   = new TH1D(*fFitInputs->Expectation().CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuMuBarCC1pi));
     TH1D * bestfit_numubarccoth   = new TH1D(*fFitInputs->Expectation().CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuMuBarCCoth));
     TH1D * bestfit_numubarnc      = new TH1D(*fFitInputs->Expectation().CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuMuBarNC));
     TH1D * bestfit_nueccqe        = new TH1D(*fFitInputs->Expectation().CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuECCQE));
     TH1D * bestfit_nuecc1pi       = new TH1D(*fFitInputs->Expectation().CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuECC1pi));
     TH1D * bestfit_nueccoth       = new TH1D(*fFitInputs->Expectation().CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuECCoth));
     TH1D * bestfit_nuenc          = new TH1D(*fFitInputs->Expectation().CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuENC));

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
     bestfit_numucc1pi      -> SetFillColor ( kBlack   );
     bestfit_numuccoth      -> SetFillColor ( kBlack   );
     bestfit_nc             -> SetFillColor ( kRed     );
     bestfit_numubarccqe    -> SetFillColor ( kAzure-3 );
     bestfit_numubarcc1pi   -> SetFillColor ( kAzure-3 );
     bestfit_numubarccoth   -> SetFillColor ( kAzure-3 );
     bestfit_numubarnc      -> SetFillColor ( kAzure-3 );
     bestfit_nueccqe        -> SetFillColor ( kGreen   );
     bestfit_nuecc1pi       -> SetFillColor ( kGreen   );
     bestfit_nueccoth       -> SetFillColor ( kGreen   );
     bestfit_nuenc          -> SetFillColor ( kGreen   );
     bestfit_numuccqe       -> SetLineColor ( kBlack   );
     bestfit_numucc1pi      -> SetLineWidth ( 1 );
     bestfit_numuccoth      -> SetLineWidth ( 1 );
     bestfit_nc             -> SetLineWidth ( 1 );
     bestfit_numuccqe       -> SetLineWidth ( 1 );
     bestfit_numubarccqe    -> SetLineWidth ( 1 );
     bestfit_numubarcc1pi   -> SetLineWidth ( 1 );
     bestfit_numubarccoth   -> SetLineWidth ( 1 );
     bestfit_numubarnc      -> SetLineWidth ( 1 );
     bestfit_nueccqe        -> SetLineWidth ( 1 );
     bestfit_nuecc1pi       -> SetLineWidth ( 1 );
     bestfit_nueccoth       -> SetLineWidth ( 1 );
     bestfit_nuenc          -> SetLineWidth ( 1 );
     bestfit_tot            -> SetLineWidth ( 2 );

     // get data
     TH1D * data = new TH1D(*fFitInputs->Observation().ErecoSpectrum(conv::kSKsFCFV1RingMuLike));
     data -> SetDirectory(0);
     data -> SetMarkerStyle(20);
     data -> SetMarkerSize(1);

     data -> GetXaxis()->SetTitle("E_{reco} (GeV)");
     data -> GetYaxis()->SetTitle("Number of events");
     //data -> Draw("P");
     data -> Draw("PERR");
     bestfit_tot -> Draw("same");

     if(opt.Contains("ShowContributions")) {
        THStack * stack = new THStack("stack","spectrum contributions");
        stack->Add(bestfit_nuenc);
        stack->Add(bestfit_nuecc1pi);
        stack->Add(bestfit_nueccoth);
        stack->Add(bestfit_nueccqe);
        stack->Add(bestfit_numubarnc);
        stack->Add(bestfit_numubarcc1pi);
        stack->Add(bestfit_numubarccoth);
        stack->Add(bestfit_numubarccqe);
        stack->Add(bestfit_nc);
        stack->Add(bestfit_numuccoth);
        stack->Add(bestfit_numucc1pi);
        stack->Add(bestfit_numuccqe);
        stack->Draw("same");
        data->Draw("PERRSAME");
        TLegend * spectrum_legend = new TLegend(0.6,0.5,0.9,0.9);
        spectrum_legend->SetFillColor (0);
        spectrum_legend->SetBorderSize(0);
        spectrum_legend->AddEntry(data,                    "data", "P" );
        spectrum_legend->AddEntry(bestfit_tot,             "best fit", "L" );
        spectrum_legend->AddEntry(bestfit_numuccqe,        "best fit (#nu_{#mu} CCQE)", "L" );
        spectrum_legend->AddEntry(bestfit_numucc1pi,       "best fit (#nu_{#mu} CC1#pi)", "LF");
        spectrum_legend->AddEntry(bestfit_numuccoth,       "best fit (#nu_{#mu} CC other)", "LF");
        spectrum_legend->AddEntry(bestfit_nc,              "best fit (#nu_{#mu}/#nu_{#tau} NC)", "LF");
        spectrum_legend->AddEntry(bestfit_numubarccqe,     "best fit (#bar{#nu_{#mu}} CCQE)", "LF");
        spectrum_legend->AddEntry(bestfit_numubarcc1pi,    "best fit (#bar{#nu_{#mu}} CC1#pi)", "LF");
        spectrum_legend->AddEntry(bestfit_numubarccoth,    "best fit (#bar{#nu_{#mu}} CC other)", "LF");
        spectrum_legend->AddEntry(bestfit_numubarnc,       "best fit (#bar{#nu_{#mu}}/#bar{#nu_{#tau}} NC)", "LF");
        spectrum_legend->AddEntry(bestfit_nueccqe,         "best fit (#nu_{e}) CCQE", "LF");
        spectrum_legend->AddEntry(bestfit_nuecc1pi,        "best fit (#nu_{e}) CC1#pi", "LF");
        spectrum_legend->AddEntry(bestfit_nueccoth,        "best fit (#nu_{e}) CC other", "LF");
        spectrum_legend->AddEntry(bestfit_nuenc,           "best fit (#nu_{e}) NC", "LF");
        spectrum_legend->Draw();
     }//ShowContributions?

     c->Update();

     if(file_out) {
       data                   -> Write("data");
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
     }

     //
     // plot contour
     //

     c->cd(2);

     // don't trust MINUIT, just get the chisq surface by yourself
     if(fOption.Contains("minuit")) {
       this->BuildChisqSurface();
     } //minuit?

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

     if(fKnowActualParams) {
        TMarker * true_point =
            new TMarker(fSinsq2Theta23, fDm23sq, 20);
        true_point->SetMarkerSize(1);
        true_point->SetMarkerColor(kRed);
        true_point->Draw();
        if(file_out) {
           true_point->Write("point_true");
        }
     }
     
     if(opt.Contains("Minos90CL")) {
         TFile f("./inputs/test/minos_90cl_contour.root","read");
         TGraph * gr = (TGraph*)f.Get("Graph");
	 TGraph * minos_90cl = new TGraph(*gr);
         minos_90cl->SetLineColor(kRed);
         minos_90cl->SetLineStyle(kSolid);
         minos_90cl->SetLineWidth(1);
         minos_90cl->Draw("l");
	 f.Close();
     }

     c->GetPad(2)->Update();

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

  if(file_out) {
    file_out->Close();
    delete file_out;
  }
}
//____________________________________________________________________________ 
