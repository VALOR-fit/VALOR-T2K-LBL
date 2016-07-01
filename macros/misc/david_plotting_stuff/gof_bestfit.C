
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

void gof_bestfit
(
 const int NfreePar = 3, // sinsqth23 Vs dcp Vs MH
 //const int NfreePar = 2, // sinsqth23 Vs dcp (no MH)

 const char * bestfitfile = "../TECH_NOTE_14_3_16/BESTFIT/RUN17/spectrum_bestfit_fitmarg_joint_sinsqth23Vsdcp_react_nh_Run17.root",
 //const char * bestfitfile = "../TECH_NOTE_14_3_16/BESTFIT/RUN16/NEW_VARRAY_TH1D/spectrum_bestfit_fitmarg_joint_sinsqth23Vsdcp_react_nh.root",
 
 const char * datafile = "data_run1to7/SKdata_Run17.root"   // RUN 1-7  
 //const char * datafile = "data_run1to6/SKdata_Run16.root" // RUN 1-6 
 )
{ 

 const int Nnewbins_1Rmu    = 5;
 const double newbinedges_1Rmu[Nnewbins_1Rmu+1]       = {0., 0.4, 0.7, 1.0, 2.0, 30.0};
 const int Nnewbins_RHC1Rmu = 5;
 const double newbinedges_RHC1Rmu[Nnewbins_RHC1Rmu+1] = {0., 0.4, 0.7, 1.0, 2.0, 30.0};
 //const int Nnewbins_RHC1Rmu = 1;
 //const double newbinedges_RHC1Rmu[Nnewbins_RHC1Rmu+1] = {0., 30.0};
 const int Nnewbins_1Re     = 1;
 const double newbinedges_1Re[Nnewbins_1Re+1]         = {0., 30.0};
 const int Nnewbins_RHC1Re  = 1;
 const double newbinedges_RHC1Re[Nnewbins_RHC1Re+1]   = {0., 30.0};
   
  // Get Best-Fit

  TFile *fbestfit = new TFile(bestfitfile,"READ");
  TH1D *hBestFit_1Rmu    = (TH1D*) fbestfit->Get("h_bestfit_marg_1Rmu");
  TH1D *hBestFit_RHC1Rmu = (TH1D*) fbestfit->Get("h_bestfit_marg_RHC1Rmu");
  TH1D *hBestFit_1Re     = (TH1D*) fbestfit->Get("h_bestfit_marg_1Re");
  TH1D *hBestFit_RHC1Re  = (TH1D*) fbestfit->Get("h_bestfit_marg_RHC1Re");

  // New bin histograms (best-fit)
  hBestFit_1Rmu   ->Rebin(Nnewbins_1Rmu,"New_bf_1Rmu",newbinedges_1Rmu);
  hBestFit_RHC1Rmu->Rebin(Nnewbins_RHC1Rmu,"New_bf_RHC1Rmu",newbinedges_RHC1Rmu);
  hBestFit_1Re    ->Rebin(Nnewbins_1Re,"New_bf_1Re",newbinedges_1Re);
  hBestFit_RHC1Re ->Rebin(Nnewbins_RHC1Re,"New_bf_RHC1Re",newbinedges_RHC1Re);

  TH1D *hBestFit_New_1Rmu    = (TH1D*) New_bf_1Rmu->Clone("hBestFit_New_1Rmu");
  TH1D *hBestFit_New_1Re     = (TH1D*) New_bf_1Re->Clone("hBestFit_New_1Re");
  TH1D *hBestFit_New_RHC1Rmu = (TH1D*) New_bf_RHC1Rmu->Clone("hBestFit_New_RHC1Rmu");
  TH1D *hBestFit_New_RHC1Re  = (TH1D*) New_bf_RHC1Re->Clone("hBestFit_New_RHC1Re");
  
 
  // Get Data
    
  TString datafile2(Form("%s/inputs/%s", gSystem->Getenv("T2KOSCANA"), datafile));
  cout << endl;
  cout << "Reading data from: " << datafile2 << endl; 
  cout << endl;

  SuperkObservation * skobs_fakedata = new SuperkObservation;
  skobs_fakedata->Read(datafile2.Data(),0,false);
  
  TH2D *hErecoVsQ2Data_1Rmu = (TH2D*)skobs_fakedata->TwoDSpectrum(kSKsFCFV1RingMuLike);
  TH1D *hErecoData_1Rmu     = (TH1D*)hErecoVsQ2Data_1Rmu->ProjectionY()->Clone("hErecoData_1Rmu");

  TH2D *hErecoVsQ2Data_1Re = (TH2D*)skobs_fakedata->TwoDSpectrum(kSKsFCFV1RingELike);
  TH1D *hErecoData_1Re     = (TH1D*)hErecoVsQ2Data_1Re->ProjectionY()->Clone("hErecoData_1Re");

  TH2D *hErecoVsQ2Data_RHC1Rmu = (TH2D*)skobs_fakedata->TwoDSpectrum(kSKsRHCFCFV1RingMuLike);
  TH1D *hErecoData_RHC1Rmu     = (TH1D*)hErecoVsQ2Data_RHC1Rmu->ProjectionY()->Clone("hErecoData_RHC1Rmu");

  TH2D *hErecoVsQ2Data_RHC1Re = (TH2D*)skobs_fakedata->TwoDSpectrum(kSKsRHCFCFV1RingELike);
  TH1D *hErecoData_RHC1Re     = (TH1D*)hErecoVsQ2Data_RHC1Re->ProjectionY()->Clone("hErecoData_RHC1Re");

  // New bin histograms (observed)
  
  hErecoData_1Rmu   ->Rebin(Nnewbins_1Rmu,"New_data_1Rmu",newbinedges_1Rmu);
  hErecoData_1Re    ->Rebin(Nnewbins_1Re,"New_data_1Re",newbinedges_1Re);
  hErecoData_RHC1Rmu->Rebin(Nnewbins_RHC1Rmu,"New_data_RHC1Rmu",newbinedges_RHC1Rmu);
  hErecoData_RHC1Re ->Rebin(Nnewbins_RHC1Re,"New_data_RHC1Re",newbinedges_RHC1Re);
 
  TH1D *hErecoData_New_1Rmu    = (TH1D*) New_data_1Rmu->Clone("hErecoData_New_1Rmu");
  TH1D *hErecoData_New_1Re     = (TH1D*) New_data_1Re->Clone("hErecoData_New_1Re");
  TH1D *hErecoData_New_RHC1Rmu = (TH1D*) New_data_RHC1Rmu->Clone("hErecoData_New_RHC1Rmu");
  TH1D *hErecoData_New_RHC1Re  = (TH1D*) New_data_RHC1Re->Clone("hErecoData_New_RHC1Re");


  // Calculate gof

  double chisq_tot = 0.;
  double gof_tot = 0.;

  for(int ibin=0;ibin<hErecoData_New_1Rmu->GetNbinsX();ibin++){
    double obs = hErecoData_New_1Rmu->GetBinContent(ibin+1);
    double exp = hBestFit_New_1Rmu->GetBinContent(ibin+1);
    chisq_tot += valor::util::ChisqPoisson(obs,exp);
  }
  for(int ibin=0;ibin<hErecoData_New_1Re->GetNbinsX();ibin++){
    double obs = hErecoData_New_1Re->GetBinContent(ibin+1);
    double exp = hBestFit_New_1Re->GetBinContent(ibin+1);
    chisq_tot += valor::util::ChisqPoisson(obs,exp);
  }
 for(int ibin=0;ibin<hErecoData_New_RHC1Rmu->GetNbinsX();ibin++){
    double obs = hErecoData_New_RHC1Rmu->GetBinContent(ibin+1);
    double exp = hBestFit_New_RHC1Rmu->GetBinContent(ibin+1);
    chisq_tot += valor::util::ChisqPoisson(obs,exp);
  }
 for(int ibin=0;ibin<hErecoData_New_RHC1Re->GetNbinsX();ibin++){
    double obs = hErecoData_New_RHC1Re->GetBinContent(ibin+1);
    double exp = hBestFit_New_RHC1Re->GetBinContent(ibin+1);
    chisq_tot += valor::util::ChisqPoisson(obs,exp);
  }

 double ndf = (Nnewbins_1Rmu + Nnewbins_1Re + Nnewbins_RHC1Rmu + Nnewbins_RHC1Re) - NfreePar;


  // Print # of events
  
  cout << endl;
  cout << "///////////////// 1Rmu ////////////////////" << endl;
  cout << endl;
  cout << "- Original binning:" << endl;
  cout << "  Tot # of events (observed) = " << hErecoData_1Rmu->Integral() << endl;
  cout << "  Tot # of events (best-fit) = " << hBestFit_1Rmu->Integral() << endl;
  // cout << "  Best-fit:" << endl;
  // for(int ibin=0;ibin<hBestFit_1Rmu->GetNbinsX();ibin++){
  //   cout << hBestFit_1Rmu->GetBinContent(ibin+1) << ", ";
  // }
  // cout << endl;
  cout << "- New binning:" << endl;
  cout << "   Tot # of events (observed) = " << hErecoData_New_1Rmu->Integral() << endl;
  cout << "   Tot # of events (best-fit) = " << hBestFit_New_1Rmu->Integral() << endl;
  // cout << "   Best-fit:" << endl;
  // for(int ibin=0;ibin<hBestFit_New_1Rmu->GetNbinsX();ibin++){
  //   cout << hBestFit_New_1Rmu->GetBinContent(ibin+1) << ", ";
  // }
  // cout << endl;
  cout << endl;
  cout << "///////////////// 1Re ////////////////////" << endl;
  cout << endl;
  cout << "- Original binning:" << endl;
  cout << "  Tot # of events (observed) = " << hErecoData_1Re->Integral() << endl;
  cout << "  Tot # of events (best-fit) = " << hBestFit_1Re->Integral() << endl;
  cout << "- New binning:" << endl;
  cout << "   Tot # of events (observed) = " << hErecoData_New_1Re->Integral() << endl;
  cout << "   Tot # of events (best-fit) = " << hBestFit_New_1Re->Integral() << endl;
  cout << endl;
  cout << "///////////////// RHC1Rmu ////////////////////" << endl;
  cout << endl;
  cout << "- Original binning:" << endl;
  cout << "  Tot # of events (observed) = " << hErecoData_RHC1Rmu->Integral() << endl;
  cout << "  Tot # of events (best-fit) = " << hBestFit_RHC1Rmu->Integral() << endl;
  cout << "- New binning:" << endl;
  cout << "   Tot # of events (observed) = " << hErecoData_New_RHC1Rmu->Integral() << endl;
  cout << "   Tot # of events (best-fit) = " << hBestFit_New_RHC1Rmu->Integral() << endl;
  cout << endl;
  cout << "///////////////// RHC1Re ////////////////////" << endl;
  cout << endl;
  cout << "- Original binning:" << endl;
  cout << "  Tot # of events (observed) = " << hErecoData_RHC1Re->Integral() << endl;
  cout << "  Tot # of events (best-fit) = " << hBestFit_RHC1Re->Integral() << endl;
  cout << endl;
  cout << "- New binning:" << endl;
  cout << "   Tot # of events (observed) = " << hErecoData_New_RHC1Re->Integral() << endl;
  cout << "   Tot # of events (best-fit) = " << hBestFit_New_RHC1Re->Integral() << endl;
  cout << endl;
  cout << "///////////////// Goodness of fit result ////////////////////" << endl;
  cout << endl;
  cout << "chisq / ndf = " << chisq_tot << " / " << ndf << " = " << chisq_tot / ndf  << endl;
  cout << "prob = " << TMath::Prob(chisq_tot,ndf) << endl;
  cout << endl;

  

  // Draw histograms

  // 1Rmu

  TCanvas *cEreco_BestFitVsData_1Rmu = new TCanvas("cEreco_BestFitVsData_1Rmu","cEreco_BestFitVsData_1Rmu");
  hErecoData_1Rmu->DrawCopy("pe");
  hBestFit_1Rmu->SetLineColor(2);
  hBestFit_1Rmu->SetLineWidth(2);
  hBestFit_1Rmu->DrawCopy("same");
  hErecoData_1Rmu->SetLineColor(1);
  hErecoData_1Rmu->SetLineWidth(2);
  hErecoData_1Rmu->DrawCopy("pe same");
  TLegend *leg_Ereco_1Rmu = new TLegend(0.7, 0.5, 0.9, 0.7);
  leg_Ereco_1Rmu->SetFillColor(0);
  leg_Ereco_1Rmu->SetBorderSize(0);
  leg_Ereco_1Rmu->AddEntry(hErecoData_1Rmu,"Data","P");
  leg_Ereco_1Rmu->AddEntry(hBestFit_1Rmu,"Exp","L");
  leg_Ereco_1Rmu->Draw();

  TCanvas *cEreco_BestFitVsData_New_1Rmu = new TCanvas("cEreco_BestFitVsData_New_1Rmu","cEreco_BestFitVsData_New_1Rmu");
  hErecoData_New_1Rmu->DrawCopy("pe");
  hBestFit_New_1Rmu->SetLineColor(2);
  hBestFit_New_1Rmu->SetLineWidth(2);
  hBestFit_New_1Rmu->DrawCopy("same");
  hErecoData_New_1Rmu->SetLineColor(1);
  hErecoData_New_1Rmu->SetLineWidth(2);
  hErecoData_New_1Rmu->DrawCopy("pe same");
  TLegend *leg_Ereco_New_1Rmu = new TLegend(0.7, 0.5, 0.9, 0.7);
  leg_Ereco_New_1Rmu->SetFillColor(0);
  leg_Ereco_New_1Rmu->SetBorderSize(0);
  leg_Ereco_New_1Rmu->AddEntry(hErecoData_New_1Rmu,"Data","P");
  leg_Ereco_New_1Rmu->AddEntry(hBestFit_New_1Rmu,"Exp","L");
  leg_Ereco_New_1Rmu->Draw();

  // 1Re

  TCanvas *cEreco_BestFitVsData_1Re = new TCanvas("cEreco_BestFitVsData_1Re","cEreco_BestFitVsData_1Re");
  hErecoData_1Re->DrawCopy("pe");
  hBestFit_1Re->SetLineColor(2);
  hBestFit_1Re->SetLineWidth(2);
  hBestFit_1Re->DrawCopy("same");
  hErecoData_1Re->SetLineColor(1);
  hErecoData_1Re->SetLineWidth(2);
  hErecoData_1Re->DrawCopy("pe same");
  TLegend *leg_Ereco_1Re = new TLegend(0.7, 0.5, 0.9, 0.7);
  leg_Ereco_1Re->SetFillColor(0);
  leg_Ereco_1Re->SetBorderSize(0);
  leg_Ereco_1Re->AddEntry(hErecoData_1Re,"Data","P");
  leg_Ereco_1Re->AddEntry(hBestFit_1Re,"Exp","L");
  leg_Ereco_1Re->Draw();

  TCanvas *cEreco_BestFitVsData_New_1Re = new TCanvas("cEreco_BestFitVsData_New_1Re","cEreco_BestFitVsData_New_1Re");
  hErecoData_New_1Re->DrawCopy("pe");
  hBestFit_New_1Re->SetLineColor(2);
  hBestFit_New_1Re->SetLineWidth(2);
  hBestFit_New_1Re->DrawCopy("same");
  hErecoData_New_1Re->SetLineColor(1);
  hErecoData_New_1Re->SetLineWidth(2);
  hErecoData_New_1Re->DrawCopy("pe same");
  TLegend *leg_Ereco_New_1Re = new TLegend(0.7, 0.5, 0.9, 0.7);
  leg_Ereco_New_1Re->SetFillColor(0);
  leg_Ereco_New_1Re->SetBorderSize(0);
  leg_Ereco_New_1Re->AddEntry(hErecoData_New_1Re,"Data","P");
  leg_Ereco_New_1Re->AddEntry(hBestFit_New_1Re,"Exp","L");
  leg_Ereco_New_1Re->Draw();

  // RHC1Rmu

  TCanvas *cEreco_BestFitVsData_RHC1Rmu = new TCanvas("cEreco_BestFitVsData_RHC1Rmu","cEreco_BestFitVsData_RHC1Rmu");
  hErecoData_RHC1Rmu->DrawCopy("pe");
  hBestFit_RHC1Rmu->SetLineColor(2);
  hBestFit_RHC1Rmu->SetLineWidth(2);
  hBestFit_RHC1Rmu->DrawCopy("same");
  hErecoData_RHC1Rmu->SetLineColor(1);
  hErecoData_RHC1Rmu->SetLineWidth(2);
  hErecoData_RHC1Rmu->DrawCopy("pe same");
  TLegend *leg_Ereco_RHC1Rmu = new TLegend(0.7, 0.5, 0.9, 0.7);
  leg_Ereco_RHC1Rmu->SetFillColor(0);
  leg_Ereco_RHC1Rmu->SetBorderSize(0);
  leg_Ereco_RHC1Rmu->AddEntry(hErecoData_RHC1Rmu,"Data","P");
  leg_Ereco_RHC1Rmu->AddEntry(hBestFit_RHC1Rmu,"Exp","L");
  leg_Ereco_RHC1Rmu->Draw();

  TCanvas *cEreco_BestFitVsData_New_RHC1Rmu = new TCanvas("cEreco_BestFitVsData_New_RHC1Rmu","cEreco_BestFitVsData_New_RHC1Rmu");
  hErecoData_New_RHC1Rmu->DrawCopy("pe");
  hBestFit_New_RHC1Rmu->SetLineColor(2);
  hBestFit_New_RHC1Rmu->SetLineWidth(2);
  hBestFit_New_RHC1Rmu->DrawCopy("same");
  hErecoData_New_RHC1Rmu->SetLineColor(1);
  hErecoData_New_RHC1Rmu->SetLineWidth(2);
  hErecoData_New_RHC1Rmu->DrawCopy("pe same");
  TLegend *leg_Ereco_New_RHC1Rmu = new TLegend(0.7, 0.5, 0.9, 0.7);
  leg_Ereco_New_RHC1Rmu->SetFillColor(0);
  leg_Ereco_New_RHC1Rmu->SetBorderSize(0);
  leg_Ereco_New_RHC1Rmu->AddEntry(hErecoData_New_RHC1Rmu,"Data","P");
  leg_Ereco_New_RHC1Rmu->AddEntry(hBestFit_New_RHC1Rmu,"Exp","L");
  leg_Ereco_New_RHC1Rmu->Draw();


  // RHC1Re

  TCanvas *cEreco_BestFitVsData_RHC1Re = new TCanvas("cEreco_BestFitVsData_RHC1Re","cEreco_BestFitVsData_RHC1Re");
  hErecoData_RHC1Re->DrawCopy("pe");
  hBestFit_RHC1Re->SetLineColor(2);
  hBestFit_RHC1Re->SetLineWidth(2);
  hBestFit_RHC1Re->DrawCopy("same");
  hErecoData_RHC1Re->SetLineColor(1);
  hErecoData_RHC1Re->SetLineWidth(2);
  hErecoData_RHC1Re->DrawCopy("pe same");
  TLegend *leg_Ereco_RHC1Re = new TLegend(0.7, 0.5, 0.9, 0.7);
  leg_Ereco_RHC1Re->SetFillColor(0);
  leg_Ereco_RHC1Re->SetBorderSize(0);
  leg_Ereco_RHC1Re->AddEntry(hErecoData_RHC1Re,"Data","P");
  leg_Ereco_RHC1Re->AddEntry(hBestFit_RHC1Re,"Exp","L");
  leg_Ereco_RHC1Re->Draw();

  TCanvas *cEreco_BestFitVsData_New_RHC1Re = new TCanvas("cEreco_BestFitVsData_New_RHC1Re","cEreco_BestFitVsData_New_RHC1Re");
  hErecoData_New_RHC1Re->DrawCopy("pe");
  hBestFit_New_RHC1Re->SetLineColor(2);
  hBestFit_New_RHC1Re->SetLineWidth(2);
  hBestFit_New_RHC1Re->DrawCopy("same");
  hErecoData_New_RHC1Re->SetLineColor(1);
  hErecoData_New_RHC1Re->SetLineWidth(2);
  hErecoData_New_RHC1Re->DrawCopy("pe same");
  TLegend *leg_Ereco_New_RHC1Re = new TLegend(0.7, 0.5, 0.9, 0.7);
  leg_Ereco_New_RHC1Re->SetFillColor(0);
  leg_Ereco_New_RHC1Re->SetBorderSize(0);
  leg_Ereco_New_RHC1Re->AddEntry(hErecoData_New_RHC1Re,"Data","P");
  leg_Ereco_New_RHC1Re->AddEntry(hBestFit_New_RHC1Re,"Exp","L");
  leg_Ereco_New_RHC1Re->Draw();

    
  return;
}
