
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

void plot_ereco_on_top_of_prediction
(
 const bool draw_preliminary = true,
 const int preliminary_color = kBlack,
 
 const char * sk_sample            = "RHC1Rmu",                         /* options: "1Rmu", "1Re", "RHC1Rmu", "RHC1Re", "MultiRmu" */
 
 //const char * expfile = "../TECH_NOTE_14_3_16/SPECTRA/RUN17/files/ErecoVsTheta_Run17_Asimov_A_1Rmu.root",
 //const char * expfile = "../TECH_NOTE_14_3_16/SPECTRA/RUN17/files/ErecoVsTheta_Run17_Asimov_A_1Re.root",
 const char * expfile = "../TECH_NOTE_14_3_16/SPECTRA/RUN17/files/ErecoVsTheta_Run17_Asimov_A_RHC1Rmu.root",
 //const char * expfile = "../TECH_NOTE_14_3_16/SPECTRA/RUN17/files/ErecoVsTheta_Run17_Asimov_A_RHC1Re.root",
 
 //const char * rootfile = "data_run1to6/SKdata_Run16.root",
 const char * rootfile = "data_run1to7/SKdata_Run17.root",

 const bool drawdata = true,

 float myscale = 1.,

 // Check Sig/Bkg ratio for this point and compare Ereco-Theta Vs Ereco only 
 //const double pt_ereco = 0.347578,
 //const double pt_theta = 70.9962,
 const double pt_ereco = 0.36448,
 const double pt_theta = 113.233,


 )
{  


  // Get Expectation
  
  TFile *fexp = new TFile(expfile,"READ");
  TH1D *hErecoExp = (TH1D*) fexp->Get("hErecoTot");
  TH1D *hQ2Exp = (TH1D*) fexp->Get("hQ2Tot");

  TH1D *hErecoMode_NueBeam = (TH1D*) fexp->Get("hErecoMode_NueBeam");
  TH1D *hErecoMode_NumuBeam = (TH1D*) fexp->Get("hErecoMode_NumuBeam");
  TH1D *hErecoMode_NumuBarBeam = (TH1D*) fexp->Get("hErecoMode_NumuBarBeam");
  TH1D *hErecoMode_NumuBarToNueBar = (TH1D*) fexp->Get("hErecoMode_NumuBarToNueBar");
  TH1D *hErecoMode_NumuToNue = (TH1D*) fexp->Get("hErecoMode_NumuToNue");
  TH1D *hErecoMode_NC = (TH1D*) fexp->Get("hErecoMode_NC"); 
  TH1D *hEreco_AllBkg = (TH1D*) hErecoMode_NC->Clone("hEreco_AllBkg");
  hEreco_AllBkg->Add(hErecoMode_NumuBarBeam);
  hEreco_AllBkg->Add(hErecoMode_NumuBeam);
  hEreco_AllBkg->Add(hErecoMode_NueBeam);

  TH1D *hQ2Mode_NueBeam = (TH1D*) fexp->Get("hQ2Mode_NueBeam");
  TH1D *hQ2Mode_NumuBeam = (TH1D*) fexp->Get("hQ2Mode_NumuBeam");
  TH1D *hQ2Mode_NumuBarBeam = (TH1D*) fexp->Get("hQ2Mode_NumuBarBeam");
  TH1D *hQ2Mode_NumuBarToNueBar = (TH1D*) fexp->Get("hQ2Mode_NumuBarToNueBar");
  TH1D *hQ2Mode_NumuToNue = (TH1D*) fexp->Get("hQ2Mode_NumuToNue");
  TH1D *hQ2Mode_NC = (TH1D*) fexp->Get("hQ2Mode_NC"); 
  TH1D *hQ2_AllBkg = (TH1D*) hQ2Mode_NC->Clone("hQ2_AllBkg");
  hQ2_AllBkg->Add(hQ2Mode_NumuBarBeam);
  hQ2_AllBkg->Add(hQ2Mode_NumuBeam);
  hQ2_AllBkg->Add(hQ2Mode_NueBeam);
 
  TH2D *hErecoVsQ2Mode_NueBeam = (TH2D*) fexp->Get("hErecoVsQ2Mode_NueBeam");
  TH2D *hErecoVsQ2Mode_NumuBeam = (TH2D*) fexp->Get("hErecoVsQ2Mode_NumuBeam");
  TH2D *hErecoVsQ2Mode_NumuBarBeam = (TH2D*) fexp->Get("hErecoVsQ2Mode_NumuBarBeam");
  TH2D *hErecoVsQ2Mode_NumuBarToNueBar = (TH2D*) fexp->Get("hErecoVsQ2Mode_NumuBarToNueBar");
  TH2D *hErecoVsQ2Mode_NumuToNue = (TH2D*) fexp->Get("hErecoVsQ2Mode_NumuToNue");
  TH2D *hErecoVsQ2Mode_NC = (TH2D*) fexp->Get("hErecoVsQ2Mode_NC"); 
  TH2D *hErecoVsQ2_AllBkg = (TH2D*) hErecoVsQ2Mode_NC->Clone("hErecoVsQ2_AllBkg");
  hErecoVsQ2_AllBkg->Add(hErecoVsQ2Mode_NumuBarBeam);
  hErecoVsQ2_AllBkg->Add(hErecoVsQ2Mode_NumuBeam);
  hErecoVsQ2_AllBkg->Add(hErecoVsQ2Mode_NueBeam);
  TH2D *hErecoVsQ2_All = (TH2D*) hErecoVsQ2Mode_NC->Clone("hErecoVsQ2_All");
  hErecoVsQ2_All->Add(hErecoVsQ2Mode_NumuBarBeam);
  hErecoVsQ2_All->Add(hErecoVsQ2Mode_NumuBeam);
  hErecoVsQ2_All->Add(hErecoVsQ2Mode_NueBeam);
  hErecoVsQ2_All->Add(hErecoVsQ2Mode_NumuToNue);
  hErecoVsQ2_All->Add(hErecoVsQ2Mode_NumuBarToNueBar);

  // Normalized histograms

  TH1D *hEreco_AllBkg_norm = (TH1D*) hEreco_AllBkg->Clone("hEreco_AllBkg_norm");
  hEreco_AllBkg_norm->Scale(1./hEreco_AllBkg_norm->Integral());
  TH1D *hErecoMode_NumuToNue_norm = (TH1D*) hErecoMode_NumuToNue->Clone("hErecoMode_NumuToNue_norm");
  hErecoMode_NumuToNue_norm->Scale(1./hErecoMode_NumuToNue->Integral());
  TH1D *hErecoMode_NumuBarToNueBar_norm = (TH1D*) hErecoMode_NumuBarToNueBar->Clone("hErecoMode_NumuBarToNueBar_norm");
  hErecoMode_NumuBarToNueBar_norm->Scale(1./hErecoMode_NumuBarToNueBar->Integral());

  TH1D *hQ2_AllBkg_norm = (TH1D*) hQ2_AllBkg->Clone("hQ2_AllBkg_norm");
  hQ2_AllBkg_norm->Scale(1./hQ2_AllBkg_norm->Integral());
  TH1D *hQ2Mode_NumuToNue_norm = (TH1D*) hQ2Mode_NumuToNue->Clone("hQ2Mode_NumuToNue_norm");
  hQ2Mode_NumuToNue_norm->Scale(1./hQ2Mode_NumuToNue->Integral());
  TH1D *hQ2Mode_NumuBarToNueBar_norm = (TH1D*) hQ2Mode_NumuBarToNueBar->Clone("hQ2Mode_NumuBarToNueBar_norm");
  hQ2Mode_NumuBarToNueBar_norm->Scale(1./hQ2Mode_NumuBarToNueBar->Integral());


  
  // Get Data
    
  TString rootfile2(Form("%s/inputs/%s", gSystem->Getenv("T2KOSCANA"), rootfile));
  cout << endl;
  cout << "Reading data from: " << rootfile2 << endl; 
  cout << endl;

  SuperkObservation * skobs_fakedata = new SuperkObservation;
  skobs_fakedata->Read(rootfile2.Data(),0,false);

  conv::SuperkSample_t s = util::SKSampleFromString(sk_sample); // requested sample
  if(s == conv::kSKsUndefined) {
    cerr << "** Undefined sample!" << endl;
    return;
  }
  int is = util::ArrIdx::Instance()->SuperkSampleIdx(s); // index of specified sample
  TH2D *hErecoVsQ2Data = (TH2D*)skobs_fakedata->TwoDSpectrum(s);

  TH1D *hQ2Data = (TH1D*)hErecoVsQ2Data->ProjectionX()->Clone("hQ2Data");
  TH1D *hErecoData = (TH1D*)hErecoVsQ2Data->ProjectionY()->Clone("hErecoData");
  
  
  // Print # of events
  
  cout << endl;
  cout << "Tot # of events = " << hErecoVsQ2_All->Integral() << endl;
  cout << "Tot # of events 2 = " 
       << 
    hErecoVsQ2Mode_NumuBarToNueBar->Integral() + 
    hErecoVsQ2Mode_NumuToNue->Integral() + 
    hErecoVsQ2_AllBkg->Integral()
       << endl; 
  cout << endl;
  cout << "Check expected # of events of each component:" << endl;
  cout << " - NumuBar -> NueBar = " << hErecoVsQ2Mode_NumuBarToNueBar->Integral() << endl; 
  cout << " - Numu -> Nue = " << hErecoVsQ2Mode_NumuToNue->Integral() << endl; 
  cout << " - NC = " << hErecoVsQ2Mode_NC->Integral() << endl;
  cout << " - NumuBar beam = " << hErecoVsQ2Mode_NumuBarBeam->Integral() << endl;
  cout << " - Numu beam = " << hErecoVsQ2Mode_NumuBeam->Integral() << endl;
  cout << " - Nue beam = " << hErecoVsQ2Mode_NueBeam->Integral() << endl;
  cout << endl;
  cout << "# of observed events = " << hErecoVsQ2Data->Integral() << endl;
  cout << endl;

  int data_erecoVstheta_bin = hErecoVsQ2_All->FindBin(pt_theta,pt_ereco); 
  double Sig_erecoVstheta = hErecoVsQ2Mode_NumuBarToNueBar->GetBinContent(data_erecoVstheta_bin);
  double Bkg_erecoVstheta = hErecoVsQ2_All->GetBinContent(data_erecoVstheta_bin) - Sig_erecoVstheta;
  cout << endl;
  cout << "Check RHC1Re template components of Ereco = " << pt_ereco << " and Theta = " << pt_theta << endl;  
  cout << " - Ereco Vs Theta:" << endl;
  cout << "   NumuBar -> NueBar = " << hErecoVsQ2Mode_NumuBarToNueBar->GetBinContent(data_erecoVstheta_bin) << endl; 
  cout << "   Numu -> Nue = "       << hErecoVsQ2Mode_NumuToNue->GetBinContent(data_erecoVstheta_bin) << endl; 
  cout << "   NC = "                << hErecoVsQ2Mode_NC->GetBinContent(data_erecoVstheta_bin) << endl;
  cout << "   NumuBar beam = "      << hErecoVsQ2Mode_NumuBarBeam->GetBinContent(data_erecoVstheta_bin) << endl;
  cout << "   Numu beam = "         << hErecoVsQ2Mode_NumuBeam->GetBinContent(data_erecoVstheta_bin) << endl;
  cout << "   Nue beam = "          << hErecoVsQ2Mode_NueBeam->GetBinContent(data_erecoVstheta_bin) << endl;
  cout << "   Signal = " << Sig_erecoVstheta << endl;
  cout << "   Background = " << Bkg_erecoVstheta << endl;
   cout << "   Signal / Ratio: " << Sig_erecoVstheta / Bkg_erecoVstheta << endl;
  cout << endl;
  int data_ereco_bin = hErecoVsQ2_All->ProjectionY()->FindBin(pt_ereco); 
  double Sig_ereco = hErecoVsQ2Mode_NumuBarToNueBar->ProjectionY()->GetBinContent(data_ereco_bin);
  double Bkg_ereco = hErecoVsQ2_All->ProjectionY()->GetBinContent(data_ereco_bin) - Sig_ereco;
  cout << " - Ereco only:" << endl;
  cout << "   NumuBar -> NueBar = " << hErecoVsQ2Mode_NumuBarToNueBar->ProjectionY()->GetBinContent(data_ereco_bin) << endl; 
  cout << "   Numu -> Nue = "       << hErecoVsQ2Mode_NumuToNue->ProjectionY()->GetBinContent(data_ereco_bin) << endl; 
  cout << "   NC = "                << hErecoVsQ2Mode_NC->ProjectionY()->GetBinContent(data_ereco_bin) << endl;
  cout << "   NumuBar beam = "      << hErecoVsQ2Mode_NumuBarBeam->ProjectionY()->GetBinContent(data_ereco_bin) << endl;
  cout << "   Numu beam = "         << hErecoVsQ2Mode_NumuBeam->ProjectionY()->GetBinContent(data_ereco_bin) << endl;
  cout << "   Nue beam = "          << hErecoVsQ2Mode_NueBeam->ProjectionY()->GetBinContent(data_ereco_bin) << endl;
  cout << "   Signal = " << Sig_ereco << endl;
  cout << "   Background = " << Bkg_ereco << endl;
  cout << "   Signal / Ratio: " << Sig_ereco / Bkg_ereco << endl;
  cout << endl;
  

  // Draw 
  
  hErecoVsQ2Data->SetMarkerColor(kGray);
  //hErecoVsQ2Data->SetMarkerStyle(1);
  hErecoVsQ2Data->SetMarkerSize(3);

  TCanvas *Ereco_ExpVsData = new TCanvas("Ereco_ExpVsData","Ereco_ExpVsData");
  hErecoData->GetXaxis()->SetTitle("Reconstructed Neutrino Energy (GeV)");
  hErecoData->GetYaxis()->SetTitle("Number of events");  
  hErecoExp->SetLineColor(kBlue);
  hErecoExp->SetLineWidth(2);
  hErecoExp->DrawCopy("histo");  
  if(drawdata) hErecoData->DrawCopy("pe same");
  TLegend *leg_Ereco = new TLegend(0.6, 0.4, 0.9, 0.7);
  leg_Ereco->SetFillColor(0);
  leg_Ereco->SetBorderSize(0);
  leg_Ereco->AddEntry(hErecoExp,"Expected","L");
  leg_Ereco->AddEntry(hErecoData,"Data","P");
  leg_Ereco->Draw();
  double posX = hErecoData->GetXaxis()->GetBinCenter(2);
  double posY = hErecoData->GetMaximum()*0.7;
  TLatex *latPreliminary_1d = new TLatex(posX,posY,"T2K Run1-7b preliminary");
  latPreliminary_1d->SetTextSize(0.05);
  latPreliminary_1d->SetTextColor(preliminary_color);
  if(draw_preliminary) latPreliminary_1d->Draw("same");

  TCanvas *Q2_ExpVsData = new TCanvas("Q2_ExpVsData","Q2_ExpVsData");
  hQ2Data->GetXaxis()->SetTitle("#theta (degrees)");
  hQ2Data->GetYaxis()->SetTitle("Number of events");  
  hQ2Exp->SetLineColor(kBlue);
  hQ2Exp->SetLineWidth(2);
  hQ2Data->DrawCopy("pe");
  hQ2Exp->DrawCopy("histo same");  
  TLegend *leg_Q2 = new TLegend(0.6, 0.4, 0.9, 0.7);
  leg_Q2->SetFillColor(0);
  leg_Q2->SetBorderSize(0);
  leg_Q2->AddEntry(hQ2Exp,"Expected","L");
  leg_Q2->AddEntry(hQ2Data,"Data","P");
  leg_Q2->Draw();

  TCanvas *ErecoVsQ2_ExpVsData = new TCanvas("ErecoVsQ2_ExpVsData","ErecoVsQ2_ExpVsData");
  hErecoVsQ2Data->GetYaxis()->SetTitle("Reconstructed Neutrino Energy (GeV)");
  hErecoVsQ2Data->GetXaxis()->SetTitle("#theta (degrees)");
  hErecoVsQ2Data->GetZaxis()->SetTitle("Number of events");  
  hErecoVsQ2_All->SetLineColor(kBlue);
  hErecoVsQ2_All->SetLineWidth(2);
  hErecoVsQ2_All->DrawCopy("colz");  
  if(drawdata) hErecoVsQ2Data->DrawCopy("p same");
  posX = hErecoVsQ2Data->GetXaxis()->GetBinCenter(1);
  posY = hErecoVsQ2Data->GetYaxis()->GetBinCenter(2);
  TLatex *latPreliminary_2d = new TLatex(posX,posY,"T2K Run1-7b preliminary");
  latPreliminary_2d->SetTextSize(0.05);
  latPreliminary_2d->SetTextColor(preliminary_color);
  if(draw_preliminary) latPreliminary_2d->Draw("same");
  
  TCanvas *ErecoVsQ2_NueBeam_ExpVsData = new TCanvas("ErecoVsQ2_NueBeam_ExpVsData","ErecoVsQ2_NueBeam_ExpVsData");
  hErecoVsQ2Data->GetYaxis()->SetTitle("Reconstructed Neutrino Energy (GeV)");
  hErecoVsQ2Data->GetXaxis()->SetTitle("#theta (degrees)");
  hErecoVsQ2Data->GetZaxis()->SetTitle("Number of events");  
  hErecoVsQ2Mode_NueBeam->SetLineColor(kBlue);
  hErecoVsQ2Mode_NueBeam->SetLineWidth(2);
  hErecoVsQ2Mode_NueBeam->DrawCopy("colz");  
  hErecoVsQ2Data->DrawCopy("p same");

  TCanvas *ErecoVsQ2_NumuBeam_ExpVsData = new TCanvas("ErecoVsQ2_NumuBeam_ExpVsData","ErecoVsQ2_NumuBeam_ExpVsData");
  hErecoVsQ2Data->GetYaxis()->SetTitle("Reconstructed Neutrino Energy (GeV)");
  hErecoVsQ2Data->GetXaxis()->SetTitle("#theta (degrees)");
  hErecoVsQ2Data->GetZaxis()->SetTitle("Number of events");  
  hErecoVsQ2Mode_NumuBeam->SetLineColor(kBlue);
  hErecoVsQ2Mode_NumuBeam->SetLineWidth(2);
  hErecoVsQ2Mode_NumuBeam->DrawCopy("colz");  
  hErecoVsQ2Data->DrawCopy("p same");

  TCanvas *ErecoVsQ2_NumuBarBeam_ExpVsData = new TCanvas("ErecoVsQ2_NumuBarBeam_ExpVsData","ErecoVsQ2_NumuBarBeam_ExpVsData");
  hErecoVsQ2Data->GetYaxis()->SetTitle("Reconstructed Neutrino Energy (GeV)");
  hErecoVsQ2Data->GetXaxis()->SetTitle("#theta (degrees)");
  hErecoVsQ2Data->GetZaxis()->SetTitle("Number of events");  
  hErecoVsQ2Mode_NumuBarBeam->SetLineColor(kBlue);
  hErecoVsQ2Mode_NumuBarBeam->SetLineWidth(2);
  hErecoVsQ2Mode_NumuBarBeam->DrawCopy("colz");  
  hErecoVsQ2Data->DrawCopy("p same");

  TCanvas *ErecoVsQ2_NumuBarToNueBar_ExpVsData = new TCanvas("ErecoVsQ2_NumuBarToNueBar_ExpVsData","ErecoVsQ2_NumuBarToNueBar_ExpVsData");
  hErecoVsQ2Data->GetYaxis()->SetTitle("Reconstructed Neutrino Energy (GeV)");
  hErecoVsQ2Data->GetXaxis()->SetTitle("#theta (degrees)");
  hErecoVsQ2Data->GetZaxis()->SetTitle("Number of events");  
  hErecoVsQ2Mode_NumuBarToNueBar->SetLineColor(kBlue);
  hErecoVsQ2Mode_NumuBarToNueBar->SetLineWidth(2);
  hErecoVsQ2Mode_NumuBarToNueBar->DrawCopy("colz");  
  hErecoVsQ2Data->DrawCopy("p same");

  TCanvas *ErecoVsQ2_NumuToNue_ExpVsData = new TCanvas("ErecoVsQ2_NumuToNue_ExpVsData","ErecoVsQ2_NumuToNue_ExpVsData");
  hErecoVsQ2Data->GetYaxis()->SetTitle("Reconstructed Neutrino Energy (GeV)");
  hErecoVsQ2Data->GetXaxis()->SetTitle("#theta (degrees)");
  hErecoVsQ2Data->GetZaxis()->SetTitle("Number of events");  
  hErecoVsQ2Mode_NumuToNue->SetLineColor(kBlue);
  hErecoVsQ2Mode_NumuToNue->SetLineWidth(2);
  hErecoVsQ2Mode_NumuToNue->DrawCopy("colz");  
  hErecoVsQ2Data->DrawCopy("p same");

  TCanvas *ErecoVsQ2_NC_ExpVsData = new TCanvas("ErecoVsQ2_NC_ExpVsData","ErecoVsQ2_NC_ExpVsData");
  hErecoVsQ2Data->GetYaxis()->SetTitle("Reconstructed Neutrino Energy (GeV)");
  hErecoVsQ2Data->GetXaxis()->SetTitle("#theta (degrees)");
  hErecoVsQ2Data->GetZaxis()->SetTitle("Number of events");  
  hErecoVsQ2Mode_NC->SetLineColor(kBlue);
  hErecoVsQ2Mode_NC->SetLineWidth(2);
  hErecoVsQ2Mode_NC->DrawCopy("colz");  
  hErecoVsQ2Data->DrawCopy("p same");
  
  TCanvas *ErecoVsQ2_Components_ExpVsData = new TCanvas("ErecoVsQ2_Components_ExpVsData","ErecoVsQ2_Components_ExpVsData");
  hErecoVsQ2Data->SetMarkerSize(1.5);
  ErecoVsQ2_Components_ExpVsData->Divide(2,3);
  ErecoVsQ2_Components_ExpVsData->cd(1);
  hErecoVsQ2Mode_NumuBarToNueBar->Draw("colz");
  hErecoVsQ2Data->DrawCopy("p same");
  ErecoVsQ2_Components_ExpVsData->cd(2);
  hErecoVsQ2Mode_NumuToNue->Draw("colz");
  hErecoVsQ2Data->DrawCopy("p same");
  ErecoVsQ2_Components_ExpVsData->cd(3);
  hErecoVsQ2Mode_NC->Draw("colz");
  hErecoVsQ2Data->DrawCopy("p same");
  ErecoVsQ2_Components_ExpVsData->cd(4);
  hErecoVsQ2Mode_NueBeam->Draw("colz");
  hErecoVsQ2Data->DrawCopy("p same");
  ErecoVsQ2_Components_ExpVsData->cd(5);
  hErecoVsQ2Mode_NumuBeam->Draw("colz");
  hErecoVsQ2Data->DrawCopy("p same");
  ErecoVsQ2_Components_ExpVsData->cd(6);
  hErecoVsQ2Mode_NumuBarBeam->Draw("colz");
  hErecoVsQ2Data->DrawCopy("p same");

  









  
  TCanvas *cEreco_ExpVsData_Component = new TCanvas("cEreco_ExpVsData_Component","cEreco_ExpVsData_Component");
  //hErecoMode_NC->SetFillColor(0);
  //hErecoMode_NC->SetLineColor(kRed);
  hEreco_AllBkg->SetFillColor(0);
  hEreco_AllBkg->SetLineColor(kBlue);
  hErecoMode_NumuToNue->SetFillColor(0);
  hErecoMode_NumuToNue->SetLineColor(kGray+3);
  hErecoMode_NumuBarToNueBar->SetFillColor(0);
  hErecoMode_NumuBarToNueBar->SetLineColor(kGreen);
  //scale hint1 to the pad coordinates
  Float_t rightmax = 1.1*hErecoData->GetMaximum();
  //Float_t scale = gPad->GetUymax()/rightmax * myscale;
  Float_t scale = myscale;
  hErecoMode_NumuBarToNueBar->Scale(scale);
  hErecoMode_NumuBarToNueBar->Draw("");
  hEreco_AllBkg->Scale(scale);
  hEreco_AllBkg->DrawCopy("same");
  //hErecoMode_NC->Scale(scale);
  //hErecoMode_NC->DrawCopy("same");
  hErecoMode_NumuToNue->Scale(scale);
  hErecoMode_NumuToNue->DrawCopy("same");
  if(drawdata) hErecoData->DrawCopy("pe same");
 
  TLegend *leg_component = new TLegend(0.6, 0.4, 0.9, 0.7);
  leg_component->SetFillColor(0);
  leg_component->SetBorderSize(0);
  leg_component->AddEntry(hErecoMode_NumuToNue,"#nu_{#mu} #rightarrow #nu_{e}","L");
  leg_component->AddEntry(hErecoMode_NumuBarToNueBar,"#bar{#nu}_{#mu} #rightarrow #bar{#nu}_{e}","L");
  leg_component->AddEntry(hEreco_AllBkg,"Bkg","L");
  if(drawdata) leg_component->AddEntry(hErecoData,"Data","P");
  leg_component->Draw();

  char name[200];
  sprintf(name,"MC #times %i",scale);
  TLatex *tex_comp = new TLatex(0.7,0.7,name);
  tex_comp->SetTextSize(0.07);
  tex_comp->Draw();

  return;


}
