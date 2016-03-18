
//example of macro illustrating how to superimpose two histograms
#include "TCanvas.h"
#include "TStyle.h"
#include "TH1.h"
#include "TGaxis.h"
#include "TRandom.h"


//
// It takes expected histogram from test_sk_ereco_spectrum_nominal_mainmode_nuebarapp.C output
//

using namespace std;
using namespace valor;
using namespace valor::conv;

void plot_bestfit_on_top_of_data_and_prediction
(
 const char * sk_sample            = "1Rmu",                         /* options: "1Rmu", "1Re", "MultiRmu" */
 
 //const char * expfile = "spectra/template_ErecoVsTheta_2D_RHC1Rmu_RHC1Rmu.root",
 //const char * expfile = "spectra/template_ErecoVsTheta_2D_notosc_RHC1Rmu.root",

 //const char * expfile = "spectra/template_ErecoVsTheta_2D_1Rmu_1Rmu.root",
 //const char * expfile = "spectra/template_ErecoVsTheta_2D_notosc_1Rmu.root",

 //const char * expfile = "spectra/template_ErecoVsTheta_2D_RHC1Re_RHC1Re.root",
 //const char * expfile = "spectra/template_ErecoVsTheta_2D_notosc_RHC1Re.root",

 //const char * expfile = "spectra/template_ErecoVsTheta_2D_1Re_1Re.root",
 //const char * expfile = "spectra/template_ErecoVsTheta_2D_notosc_1Re.root",

 const char * datafile = "data_run1to6/SKdata_Run16.root",
 //const char * datafile = "../outputs/test/SKdata_Run16.root",

 const char * bestfitfile = "bestfit/spectrum_bestfit_fitmarg_joint_dcp_react_nh.root",

 const bool draw_preliminary = true
 
 )
{ 
  // Get SK sample
  
  conv::SuperkSample_t s = util::SKSampleFromString(sk_sample); // requested sample
  if(s == conv::kSKsUndefined) {
    cerr << "** Undefined sample!" << endl;
    return;
  }
  int is = util::ArrIdx::Instance()->SuperkSampleIdx(s); // index of specified sample

  
  // Get Best-Fit
  
  TFile *fbestfit = new TFile(bestfitfile,"READ");
  TH1D *hBestFit;
  if(s==kSKsFCFV1RingELike)          hBestFit = (TH1D*) fbestfit->Get("h_bestfit_marg_1Re");
  else if(s==kSKsRHCFCFV1RingELike)  hBestFit = (TH1D*) fbestfit->Get("h_bestfit_marg_RHC1Re");
  else if(s==kSKsFCFV1RingMuLike)    hBestFit = (TH1D*) fbestfit->Get("h_bestfit_marg_1Rmu");
  else if(s==kSKsRHCFCFV1RingMuLike) hBestFit = (TH1D*) fbestfit->Get("h_bestfit_marg_RHC1Rmu");
  else{
    cerr << endl;
    cerr << "Sample " << util::AsStringCompact(s) << " is not defined" << endl;
    cerr << endl;
    exit(1);
  }


  
  // Get Expectation
  
  TFile *fexp = new TFile(expfile,"READ");
  TH1D *hErecoExp = (TH1D*) fexp->Get("hErecoTot");
  TH1D *hQ2Exp = (TH1D*) fexp->Get("hQ2Tot");
  
  
  // Get Data
    
  TString datafile2(Form("%s/inputs/%s", gSystem->Getenv("T2KOSCANA"), datafile));
  cout << endl;
  cout << "Reading data from: " << datafile2 << endl; 
  cout << endl;

  SuperkObservation * skobs_fakedata = new SuperkObservation;
  skobs_fakedata->Read(datafile2.Data(),0,false);

  TH2F *hErecoVsQ2Data = (TH2F*)skobs_fakedata->TwoDSpectrum(s);
  TH1F *hQ2Data = (TH1F*)hErecoVsQ2Data->ProjectionX()->Clone("hQ2Data");
  TH1F *hErecoData = (TH1F*)hErecoVsQ2Data->ProjectionY()->Clone("hErecoData");
  
  
  // Print # of events
  
  cout << endl;
  cout << "Tot # of events (expected) = " << hErecoExp->Integral() << endl;
  cout << "Tot # of events (best-fit) = " << hBestFit->Integral() << endl;
  cout << "Tot # of events (observed) = " << hErecoVsQ2Data->Integral() << endl;
  cout << endl;

  

  // Draw 
  
  hErecoVsQ2Data->SetMarkerColor(kGray);
  //hErecoVsQ2Data->SetMarkerStyle(1);
  hErecoVsQ2Data->SetMarkerSize(3);

  TCanvas *Ereco_ExpVsData = new TCanvas("Ereco_ExpVsData","Ereco_ExpVsData");
  hErecoExp->GetXaxis()->SetTitle("Reconstructed Neutrino Energy (GeV)");
  hErecoExp->GetYaxis()->SetTitle("Number of events");  
  hErecoExp->SetLineColor(kBlue);
  hErecoExp->SetLineWidth(2);
  hErecoExp->DrawCopy("histo");  
  hBestFit->SetLineColor(kRed);
  hBestFit->SetLineWidth(2);
  hBestFit->Draw("same");
  hErecoData->DrawCopy("pe same");
  
  TLegend *leg_Ereco = new TLegend(0.6, 0.4, 0.9, 0.7);
  leg_Ereco->SetFillColor(0);
  leg_Ereco->SetBorderSize(0);
  leg_Ereco->AddEntry(hErecoExp,"Expected","L");
  leg_Ereco->AddEntry(hErecoData,"Data","P");
  leg_Ereco->AddEntry(hBestFit,"Best-fit","L");
  leg_Ereco->Draw();


  TLatex *latPreliminary = new TLatex(hBestFit->GetXaxis()->GetBinLowEdge(5),hBestFit->GetMaximum(),"PRELIMINARY");
  if(draw_preliminary) latPreliminary->Draw("same");

  return;
}
