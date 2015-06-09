#include <iomanip.h>
#include <TXMLEngine.h>
#include <string>

void compare_systematic_error_propagation_all_2015
(
 const char * infile_1 = "$T2KOSCANA/outputs/test/fouterr_prefit.root",
 const char * infile_2 = "$T2KOSCANA/outputs/test/fouterr_postfit.root",
 bool write_to_eps = false
)
{
  using namespace std;
  using namespace valor;
  using namespace valor::conv;
  
  
  int canv_width  = 1280;
  int canv_height = 720;
  
  double title_size = 0.08;
  double label_size = 0.07;
  double title_offset = 0.8;
 
  TFile *file_1 = new TFile(infile_1,"READ");
  TFile *file_2 = new TFile(infile_2,"READ");
  
  TH1D *hVarFrac_p1s_file1 = (TH1D*) file_1->Get("hVarFrac_p1s");
  TH1D *hVarFrac_m1s_file1 = (TH1D*) file_1->Get("hVarFrac_m1s");
  TH1D *hVarFrac_p3s_file1 = (TH1D*) file_1->Get("hVarFrac_p3s");
  TH1D *hVarFrac_m3s_file1 = (TH1D*) file_1->Get("hVarFrac_m3s");

  TH1D *hVarTot_p1s_file1 = (TH1D*) file_1->Get("hVarTot_p1s");
  TH1D *hVarTot_m1s_file1 = (TH1D*) file_1->Get("hVarTot_m1s");
  TH1D *hVarTot_p3s_file1 = (TH1D*) file_1->Get("hVarTot_p3s");
  TH1D *hVarTot_m3s_file1 = (TH1D*) file_1->Get("hVarTot_m3s");

  TH1D *hVarFrac_p1s_file2 = (TH1D*) file_2->Get("hVarFrac_p1s");
  TH1D *hVarFrac_m1s_file2 = (TH1D*) file_2->Get("hVarFrac_m1s");
  TH1D *hVarFrac_p3s_file2 = (TH1D*) file_2->Get("hVarFrac_p3s");
  TH1D *hVarFrac_m3s_file2 = (TH1D*) file_2->Get("hVarFrac_m3s");

  TH1D *hVarTot_p1s_file2 = (TH1D*) file_2->Get("hVarTot_p1s");
  TH1D *hVarTot_m1s_file2 = (TH1D*) file_2->Get("hVarTot_m1s");
  TH1D *hVarTot_p3s_file2 = (TH1D*) file_2->Get("hVarTot_p3s");
  TH1D *hVarTot_m3s_file2 = (TH1D*) file_2->Get("hVarTot_m3s");

  // All parameters

  const int Npt = hVarFrac_p1s_file1->GetNbinsX();

  TGraphAsymmErrors *gVarFrac_file1 = new TGraphAsymmErrors(Npt);
  for(int ipt=1;ipt<=Npt;ipt++){
    double VarP1 = hVarFrac_p1s_file1->GetBinContent(ipt) * 100;
    double VarM1 = hVarFrac_m1s_file1->GetBinContent(ipt) * 100; // must be positive
    gVarFrac_file1->SetPoint(ipt-1,ipt-1,100.);
    gVarFrac_file1->SetPointError(ipt-1,0.,0.,VarM1,VarP1);
  }
  TGraphAsymmErrors *gVarFrac_file2 = new TGraphAsymmErrors(Npt);
  for(int ipt=1;ipt<=Npt;ipt++){
    double VarP1 = hVarFrac_p1s_file2->GetBinContent(ipt) * 100;
    double VarM1 = hVarFrac_m1s_file2->GetBinContent(ipt) * 100; // must be positive
    gVarFrac_file2->SetPoint(ipt-1,ipt-1,100.);
    gVarFrac_file2->SetPointError(ipt-1,0.,0.,VarM1,VarP1);
  }
  
  TCanvas *cVarFrac_all = new TCanvas("cVarFrac_all","cVarFrac_all",canv_width,canv_height);
  gVarFrac_file1->GetXaxis()->SetTitle("# Parameter");
  gVarFrac_file1->GetYaxis()->SetTitle("Scale of event rate [%]");
  gVarFrac_file1->GetXaxis()->SetTitleSize(title_size);
  gVarFrac_file1->GetYaxis()->SetTitleSize(title_size);
  gVarFrac_file1->GetXaxis()->SetLabelSize(label_size);
  gVarFrac_file1->GetYaxis()->SetLabelSize(label_size);
  gVarFrac_file1->GetXaxis()->SetTitleOffset(title_offset);
  gVarFrac_file1->GetYaxis()->SetTitleOffset(title_offset);
  gVarFrac_file1->SetLineColor(kRed);
  gVarFrac_file2->SetLineColor(kBlue);
  gVarFrac_file1->Draw("ape");
  gVarFrac_file2->Draw("pe same");
  
  cVarFrac_all->Update();
  if (write_to_eps) {
    cVarFrac_all->Print("error_fraction_all.eps");
  }
    
  // Flux parameters

  const int NptFlux = 25;

  int gpt=0;
  TGraphAsymmErrors *gVarFrac_Flux_file1 = new TGraphAsymmErrors(NptFlux);
  for(int ipt=1;ipt<=NptFlux;ipt++){
    double VarP1 = hVarFrac_p1s_file1->GetBinContent(ipt) * 100;
    double VarM1 = hVarFrac_m1s_file1->GetBinContent(ipt) * 100; // must be positive
    gVarFrac_Flux_file1->SetPoint(gpt,gpt+1,100.);
    gVarFrac_Flux_file1->SetPointError(gpt,0.2,0.2,VarM1,VarP1);
    gpt++;
  }
  gpt=0;
  TGraphAsymmErrors *gVarFrac_Flux_file2 = new TGraphAsymmErrors(NptFlux);
  for(int ipt=1;ipt<=NptFlux;ipt++){
    double VarP1 = hVarFrac_p1s_file2->GetBinContent(ipt) * 100;
    double VarM1 = hVarFrac_m1s_file2->GetBinContent(ipt) * 100; // must be positive
    gVarFrac_Flux_file2->SetPoint(gpt,gpt+1,100.);
    gVarFrac_Flux_file2->SetPointError(gpt,0.2,0.2,VarM1,VarP1);
    gpt++;
  }

  TCanvas * cVarFrac_flux = new TCanvas("cVarFrac_flux","cVarFrac_flux",canv_width,canv_height);
  gVarFrac_Flux_file1->GetXaxis()->SetTitle("# Flux Parameter");
  gVarFrac_Flux_file1->GetYaxis()->SetTitle("Scale of event rate [%]");
  gVarFrac_Flux_file1->GetXaxis()->SetTitleSize(title_size);
  gVarFrac_Flux_file1->GetYaxis()->SetTitleSize(title_size);
  gVarFrac_Flux_file1->GetXaxis()->SetLabelSize(label_size);
  gVarFrac_Flux_file1->GetYaxis()->SetLabelSize(label_size);
  gVarFrac_Flux_file1->GetXaxis()->SetTitleOffset(title_offset);
  gVarFrac_Flux_file1->GetYaxis()->SetTitleOffset(title_offset);
  gVarFrac_Flux_file1->SetLineColor(kBlue);
  gVarFrac_Flux_file1->SetFillColor(kBlue);
  gVarFrac_Flux_file1->SetLineWidth(3);
  gVarFrac_Flux_file2->SetFillColor(kRed);
  gVarFrac_Flux_file2->SetLineColor(kRed);
  gVarFrac_Flux_file2->SetLineWidth(3);
  gVarFrac_Flux_file1->Draw("ap2");
  gVarFrac_Flux_file2->Draw("p2");
  //gVarFrac_Flux_file1->Draw("ape");
  //gVarFrac_Flux_file2->Draw("pe same");
  
  cVarFrac_flux->Update();
  if (write_to_eps) {
    cVarFrac_flux->Print("error_fraction_flux.eps");
  }


  // XSec parameters

  const int NptXSec = 12;

  gpt=0;
  TGraphAsymmErrors *gVarFrac_XSec_file1 = new TGraphAsymmErrors(NptXSec);
  for(int ipt=NptFlux+1;ipt<=(NptFlux+NptXSec);ipt++){
    double VarP1 = hVarFrac_p1s_file1->GetBinContent(ipt) * 100;
    double VarM1 = hVarFrac_m1s_file1->GetBinContent(ipt) * 100; // must be positive
    gVarFrac_XSec_file1->SetPoint(gpt,gpt+1,100.);
    gVarFrac_XSec_file1->SetPointError(gpt,0.3,0.3,VarM1,VarP1);
    gpt++;
  }
  gpt=0;
  TGraphAsymmErrors *gVarFrac_XSec_file2 = new TGraphAsymmErrors(NptXSec);
  for(int ipt=NptFlux+1;ipt<=(NptFlux+NptXSec);ipt++){
    double VarP1 = hVarFrac_p1s_file2->GetBinContent(ipt) * 100;
    double VarM1 = hVarFrac_m1s_file2->GetBinContent(ipt) * 100; // must be positive
    gVarFrac_XSec_file2->SetPoint(gpt,gpt+1,100.);
    gVarFrac_XSec_file2->SetPointError(gpt,0.3,0.3,VarM1,VarP1);
    gpt++;
  }

  TCanvas *cVarFrac_xsec = new TCanvas("cVarFrac_xsec","cVarFrac_xsec",canv_width,canv_height);
  gVarFrac_XSec_file1->GetXaxis()->SetTitle("# XSec Parameter");
  gVarFrac_XSec_file1->GetYaxis()->SetTitle("Scale of event rate [%]");
  gVarFrac_XSec_file1->GetXaxis()->SetTitleSize(title_size);
  gVarFrac_XSec_file1->GetYaxis()->SetTitleSize(title_size);
  gVarFrac_XSec_file1->GetXaxis()->SetLabelSize(label_size);
  gVarFrac_XSec_file1->GetYaxis()->SetLabelSize(label_size);
  gVarFrac_XSec_file1->GetXaxis()->SetTitleOffset(title_offset);
  gVarFrac_XSec_file1->GetYaxis()->SetTitleOffset(title_offset);
  gVarFrac_XSec_file1->SetLineColor(kBlue);
  gVarFrac_XSec_file1->SetFillColor(kBlue);
  gVarFrac_XSec_file1->SetLineWidth(3);
  gVarFrac_XSec_file2->SetFillColor(kRed);
  gVarFrac_XSec_file2->SetLineColor(kRed);
  gVarFrac_XSec_file2->SetLineWidth(3);
  gVarFrac_XSec_file1->Draw("ap2");
  gVarFrac_XSec_file2->Draw("p2");
  //gVarFrac_XSec_file1->Draw("ape");
  //gVarFrac_XSec_file2->Draw("pe same");
  
  cVarFrac_xsec->Update();
  if (write_to_eps) {
    cVarFrac_xsec->Print("error_fraction_xsec.eps");
  }


  // Det parameters

  const int NptDet = 7;

  gpt=0;
  TGraphAsymmErrors *gVarFrac_Det_file1 = new TGraphAsymmErrors(NptDet);
  for(int ipt=NptFlux+NptXSec+1;ipt<=(NptFlux+NptXSec+NptDet);ipt++){
    double VarP1 = hVarFrac_p1s_file1->GetBinContent(ipt) * 100;
    double VarM1 = hVarFrac_m1s_file1->GetBinContent(ipt) * 100; // must be positive
    gVarFrac_Det_file1->SetPoint(gpt,gpt+1,100.);
    gVarFrac_Det_file1->SetPointError(gpt,0.3,0.3,VarM1,VarP1);
    gpt++;
  }
  gpt=0;
  TGraphAsymmErrors *gVarFrac_Det_file2 = new TGraphAsymmErrors(NptDet);
  for(int ipt=NptFlux+NptXSec+1;ipt<=(NptFlux+NptXSec+NptDet);ipt++){
    double VarP1 = hVarFrac_p1s_file2->GetBinContent(ipt) * 100;
    double VarM1 = hVarFrac_m1s_file2->GetBinContent(ipt) * 100; // must be positive
    gVarFrac_Det_file2->SetPoint(gpt,gpt+1,100.);
    gVarFrac_Det_file2->SetPointError(gpt,0.3,0.3,VarM1,VarP1);
    gpt++;
  }

  TCanvas *cVarFrac_det = new TCanvas("cVarFrac_det","cVarFrac_det",canv_width,canv_height);
  gVarFrac_Det_file1->GetXaxis()->SetTitle("# Detector Parameter");
  gVarFrac_Det_file1->GetYaxis()->SetTitle("Scale of event rate [%]");
  gVarFrac_Det_file1->GetXaxis()->SetTitleSize(title_size);
  gVarFrac_Det_file1->GetYaxis()->SetTitleSize(title_size);
  gVarFrac_Det_file1->GetXaxis()->SetLabelSize(label_size);
  gVarFrac_Det_file1->GetYaxis()->SetLabelSize(label_size);
  gVarFrac_Det_file1->GetXaxis()->SetTitleOffset(title_offset);
  gVarFrac_Det_file1->GetYaxis()->SetTitleOffset(title_offset);
  gVarFrac_Det_file1->SetLineColor(kBlue);
  gVarFrac_Det_file1->SetFillColor(kBlue);
  gVarFrac_Det_file1->SetLineWidth(3);
  gVarFrac_Det_file2->SetFillColor(kRed);
  gVarFrac_Det_file2->SetLineColor(kRed);
  gVarFrac_Det_file2->SetLineWidth(3);
  gVarFrac_Det_file1->Draw("ap2");
  gVarFrac_Det_file2->Draw("p2");
  //gVarFrac_Det_file1->Draw("ape");
  //gVarFrac_Det_file2->Draw("pe same");
  
  cVarFrac_det->Update();
  if (write_to_eps) {
    cVarFrac_det->Print("error_fraction_skdet.eps");
  }



  return;
}
