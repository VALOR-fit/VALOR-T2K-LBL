
#include "Build.h"
#include "AnalysisOptions.h"
#include <vector>

using namespace std;
using namespace valor;
using namespace valor::conv;

void DrawGraph_Coverage(){
  
  const double CRITVAL_1SIGMA_1D = 1.;
  const double CRITVAL_90CL_1D =  2.71;
  const double CRITVAL_95CL_1D = 3.84;
  const double CRITVAL_2SIGMA_1D = 4.;

  TLine *lcritval_90CL = new TLine(-TMath::Pi()/2.,CRITVAL_90CL_1D,+TMath::Pi()/2.,CRITVAL_90CL_1D);
  lcritval_90CL->SetLineWidth(2);
  lcritval_90CL->SetLineStyle(2);
  lcritval_90CL->SetLineColor(kBlack);
  TLine *lcritval_2sigma = new TLine(-TMath::Pi()/2.,CRITVAL_2SIGMA_1D,+TMath::Pi()/2.,CRITVAL_2SIGMA_1D);
  lcritval_2sigma->SetLineWidth(2);
  lcritval_2sigma->SetLineStyle(2);
  lcritval_2sigma->SetLineColor(kRed);

  const double CLfrac_1sigma = 0.6827;  // 1sigma% CL
  const double CLfrac_90CL   = 0.9;  // 90% CL
  const double CLfrac_95CL   = 0.95;  // 95% CL
  const double CLfrac_2sigma = 0.9545;  // 2sigma% CL

  TLine *lfrac_90CL = new TLine(-TMath::Pi()/2.,CLfrac_90CL,+TMath::Pi()/2.,CLfrac_90CL);
  lfrac_90CL->SetLineWidth(2);
  lfrac_90CL->SetLineStyle(2);
  lfrac_90CL->SetLineColor(kBlack);
  TLine *lfrac_2sigma = new TLine(-TMath::Pi()/2.,CLfrac_2sigma,+TMath::Pi()/2.,CLfrac_2sigma);
  lfrac_2sigma->SetLineWidth(2);
  lfrac_2sigma->SetLineStyle(2);
  lfrac_2sigma->SetLineColor(kRed);


  
  // Critical values Normal Hierarchy

  TGraph *gcritval_90CL_nh = new TGraph();
  gcritval_90CL_nh->SetPoint(0,-TMath::Pi()/2.,2.02406);
  gcritval_90CL_nh->SetPoint(1,0.,2.45842);
  gcritval_90CL_nh->SetPoint(2,+TMath::Pi()/2.,2.2354701);
  gcritval_90CL_nh->GetYaxis()->SetRangeUser(0,5);
  gcritval_90CL_nh->GetYaxis()->SetTitle("#Delta #chi^{2}_{crit}");
  gcritval_90CL_nh->GetXaxis()->SetTitle("#delta_{CP}"); 
  gcritval_90CL_nh->SetLineColor(kBlack);
  gcritval_90CL_nh->SetLineWidth(2);

  TGraph *gcritval_2sigma_nh = new TGraph();
  gcritval_2sigma_nh->SetPoint(0,-TMath::Pi()/2.,3.1667);
  gcritval_2sigma_nh->SetPoint(1,0.,3.18826);
  gcritval_2sigma_nh->SetPoint(2,+TMath::Pi()/2.,3.0152642);
  gcritval_2sigma_nh->GetYaxis()->SetRangeUser(0,5);
  gcritval_2sigma_nh->GetYaxis()->SetTitle("#Delta #chi^{2}_{crit}");
  gcritval_2sigma_nh->GetXaxis()->SetTitle("#delta_{CP}"); 
  gcritval_2sigma_nh->SetLineColor(kRed);
  gcritval_2sigma_nh->SetLineWidth(2);

  TCanvas *ccritval_NH = new TCanvas("ccritval_NH","ccritval_NH");
  gcritval_2sigma_nh->Draw("apl");
  gcritval_90CL_nh->Draw("pl same");
  lcritval_2sigma->Draw();
  lcritval_90CL->Draw();
  TLegend *legcritval_NH = new TLegend(0.6,0.45,0.8,0.75);
  legcritval_NH->SetFillColor(0);
  legcritval_NH->SetLineColor(0);
  legcritval_NH->SetTextFont(132);
  legcritval_NH->SetTextSize(0.07);
  legcritval_NH->AddEntry(lcritval_90CL,"Constant","l");
  legcritval_NH->AddEntry(gcritval_90CL_nh,"Feldman-Cousins","l");
  legcritval_NH->AddEntry(gcritval_90CL_nh,"90% CL","l");
  legcritval_NH->AddEntry(gcritval_2sigma_nh,"2#sigma CL","l");
  legcritval_NH->Draw();





  // Critical values Inverted Hierarchy

  TGraph *gcritval_90CL_ih = new TGraph();
  gcritval_90CL_ih->SetPoint(0,-TMath::Pi()/2.,2.18394);
  gcritval_90CL_ih->SetPoint(1,0.,2.4278416);
  gcritval_90CL_ih->SetPoint(2,+TMath::Pi()/2.,2.5150831);
  gcritval_90CL_ih->GetYaxis()->SetRangeUser(0,5);
  gcritval_90CL_ih->GetYaxis()->SetTitle("#Delta #chi^{2}_{crit}");
  gcritval_90CL_ih->GetXaxis()->SetTitle("#delta_{CP}"); 
  gcritval_90CL_ih->SetLineColor(kBlack);
  gcritval_90CL_ih->SetLineWidth(2);

  // TGraph *gcritval_2sigma_ih = new TGraph();
  // gcritval_2sigma_ih->SetPoint(0,-TMath::Pi()/2.,);
  // gcritval_2sigma_ih->SetPoint(1,0.,);
  // gcritval_2sigma_ih->SetPoint(2,+TMath::Pi()/2.,);
  // gcritval_2sigma_ih->GetYaxis()->SetRangeUser(0,5);
  // gcritval_2sigma_ih->GetYaxis()->SetTitle("#Delta #chi^{2}_{crit}");
  // gcritval_2sigma_ih->GetXaxis()->SetTitle("#delta_{CP}"); 
  // gcritval_2sigma_ih->SetLineColor(kRed);
  // gcritval_2sigma_ih->SetLineWidth(2);

  TCanvas *ccritval_IH = new TCanvas("ccritval_IH","ccritval_IH");
  //gcritval_2sigma_ih->Draw("apl");
  gcritval_90CL_ih->Draw("apl");
  //lcritval_2sigma->Draw();
  lcritval_90CL->Draw();
  TLegend *legcritval_IH = new TLegend(0.6,0.45,0.8,0.75);
  legcritval_IH->SetFillColor(0);
  legcritval_IH->SetLineColor(0);
  legcritval_IH->SetTextFont(132);
  legcritval_IH->SetTextSize(0.07);
  legcritval_IH->AddEntry(lcritval_90CL,"Constant","l");
  legcritval_IH->AddEntry(gcritval_90CL_ih,"Feldman-Cousins","l");
  legcritval_IH->AddEntry(gcritval_90CL_ih,"90% CL","l");
  //legcritval_IH->AddEntry(gcritval_2sigma_ih,"2#sigma CL","l");
  legcritval_IH->Draw();



  


  // Coverage Normal Hierarchy
  
  TGraph *gcov_90CL_nh = new TGraph();
  gcov_90CL_nh->SetPoint(0,-TMath::Pi()/2.,0.9368);
  gcov_90CL_nh->SetPoint(1,0.,0.925);
  gcov_90CL_nh->SetPoint(2,+TMath::Pi()/2.,0.9422);
  gcov_90CL_nh->GetYaxis()->SetRangeUser(0.5,1);
  gcov_90CL_nh->GetYaxis()->SetTitle("Coverage");
  gcov_90CL_nh->GetXaxis()->SetTitle("#delta_{CP}"); 
  gcov_90CL_nh->SetLineColor(kBlack);
  gcov_90CL_nh->SetLineWidth(2);
 
  TGraph *gcov_2sigma_nh = new TGraph();
  gcov_2sigma_nh->SetPoint(0,-TMath::Pi()/2.,0.9732);
  gcov_2sigma_nh->SetPoint(1,0.,0.9808);
  gcov_2sigma_nh->SetPoint(2,+TMath::Pi()/2.,0.9786);
  gcov_2sigma_nh->GetYaxis()->SetRangeUser(0.5,1);
  gcov_2sigma_nh->GetYaxis()->SetTitle("Coverage");
  gcov_2sigma_nh->GetXaxis()->SetTitle("#delta_{CP}"); 
  gcov_2sigma_nh->SetLineColor(kRed);
  gcov_2sigma_nh->SetLineWidth(2);

  TCanvas *ccov_NH = new TCanvas("ccov_NH","ccov_NH");
  gcov_2sigma_nh->Draw("apl");
  gcov_90CL_nh->Draw("pl same");
  lfrac_90CL->Draw();
  lfrac_2sigma->Draw();
  TLegend *legcov_NH = new TLegend(0.6,0.45,0.8,0.75);
  legcov_NH->SetFillColor(0);
  legcov_NH->SetLineColor(0);
  legcov_NH->SetTextFont(132);
  legcov_NH->SetTextSize(0.07);
  legcov_NH->AddEntry(gcov_90CL_nh,"90% CL","l");
  legcov_NH->AddEntry(gcov_2sigma_nh,"2#sigma CL","l");
  legcov_NH->AddEntry(lfrac_90CL,"Constant","l");
  legcov_NH->Draw();



  // Coverage Inverted Hierarchy
  
  TGraph *gcov_90CL_ih = new TGraph();
  gcov_90CL_ih->SetPoint(0,-TMath::Pi()/2.,0.9405);
  gcov_90CL_ih->SetPoint(1,0.,0.925);
  gcov_90CL_ih->SetPoint(2,+TMath::Pi()/2.,0.9185);
  gcov_90CL_ih->GetYaxis()->SetRangeUser(0.5,1);
  gcov_90CL_ih->GetYaxis()->SetTitle("Coverage");
  gcov_90CL_ih->GetXaxis()->SetTitle("#delta_{CP}"); 
  gcov_90CL_ih->SetLineColor(kBlack);
  gcov_90CL_ih->SetLineWidth(2);
 
  // TGraph *gcov_2sigma_ih = new TGraph();
  // gcov_2sigma_ih->SetPoint(0,-TMath::Pi()/2.,0.9732);
  // gcov_2sigma_ih->SetPoint(1,0.,0.9808);
  // gcov_2sigma_ih->SetPoint(2,+TMath::Pi()/2.,0.9786);
  // gcov_2sigma_ih->GetYaxis()->SetRangeUser(0.5,1);
  // gcov_2sigma_ih->GetYaxis()->SetTitle("Coverage");
  // gcov_2sigma_ih->GetXaxis()->SetTitle("#delta_{CP}"); 
  // gcov_2sigma_ih->SetLineColor(kRed);
  // gcov_2sigma_ih->SetLineWidth(2);

  TCanvas *ccov_IH = new TCanvas("ccov_IH","ccov_IH");
  //gcov_2sigma_ih->Draw("apl");
  gcov_90CL_ih->Draw("apl");
  lfrac_90CL->Draw();
  //lfrac_2sigma->Draw();
  TLegend *legcov_IH = new TLegend(0.6,0.45,0.8,0.75);
  legcov_IH->SetFillColor(0);
  legcov_IH->SetLineColor(0);
  legcov_IH->SetTextFont(132);
  legcov_IH->SetTextSize(0.07);
  legcov_IH->AddEntry(gcov_90CL_ih,"90% CL","l");
  //legcov_IH->AddEntry(gcov_2sigma_ih,"2#sigma CL","l");
  legcov_IH->AddEntry(lfrac_90CL,"Constant","l");
  legcov_IH->Draw();

  return;
}
