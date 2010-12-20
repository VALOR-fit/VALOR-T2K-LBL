//
// ------------------------------------------------------------------------------------------------
// Plot contributions from various event modes in the predicted SuperK spectrum.
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

#include "Conventions.h"

void test_surv_prob()
{
  gSystem->Load("libGpad");  
  gSystem->Load("libPhysics");  
  gSystem->Load("libMatrix");   
  gSystem->Load("libHist");     
  gSystem->Load("libGraf");     
  gSystem->Load("libTree");     
  gSystem->Load("libRIO");   
  gSystem->Load("libXMLIO");
  gSystem->Load("libMinuit");   

  gSystem->Load("$T2KANA/lib/libT2KAnaYr1NuMuOscAnalysis.so");

  using namespace std;
  using namespace t2kana;
  using namespace t2kana::yr1numu;

  const int n = 1000;
  double emin = 0.05;
  double emax = 4.00;
  double de   = (emax-emin)/(n-1);

  double surv_prob_disap[n];
  double surv_prob_decay[n];
  double surv_prob_decoh[n];
  double e[n];
  double l_e[n];

  // oscillation & decay params
  double input_sinsq2theta  = 1.00;
  double input_dmsq         = 2.4E-3;
  double input_alpha_decay  = .001;
  double input_musq_decoh   = 0.003;

  for(int i=0; i<n; i++) {
    e[i] = emin + i*de;
    l_e[i] = TMath::Log10(295.3/e[i]);
    surv_prob_disap[i] = util::NuMuSurvProbOsc2D(input_sinsq2theta, input_dmsq,        e[i]);
    surv_prob_decay[i] = util::NuMuSurvProbDecay(input_sinsq2theta, input_alpha_decay, e[i]);
    surv_prob_decoh[i] = util::NuMuSurvProbDecoh(input_sinsq2theta, input_musq_decoh,  e[i]);
  }

  TGraph * gr_surv_prob_disap_e  = new TGraph(n,e,  surv_prob_disap);
  TGraph * gr_surv_prob_decay_e  = new TGraph(n,e,  surv_prob_decay);
  TGraph * gr_surv_prob_decoh_e  = new TGraph(n,e,  surv_prob_decoh);
  TGraph * gr_surv_prob_disap_le = new TGraph(n,l_e,surv_prob_disap);
  TGraph * gr_surv_prob_decay_le = new TGraph(n,l_e,surv_prob_decay);
  TGraph * gr_surv_prob_decoh_le = new TGraph(n,l_e,surv_prob_decoh);

  gStyle->SetOptStat(0);

  TCanvas * c1 = new TCanvas("c1","",20,20,500,500);
  TLegend * legend1 = new TLegend(0.3, 0.7, 0.85, 0.9);
  legend1->SetFillColor(0);
  legend1->SetBorderSize(0);
  gr_surv_prob_decay_e->SetLineStyle(kDotted);
  gr_surv_prob_decoh_e->SetLineStyle(kDashed);
  gr_surv_prob_disap_e->Draw("al");
  gr_surv_prob_decay_e->Draw("l");
  gr_surv_prob_decoh_e->Draw("l");
  legend1->AddEntry(gr_surv_prob_disap_e, Form("Oscillations: #Deltam^{2}=%.5f eV^{2}, sin^{2}2#theta=%.2f",input_dmsq,input_sinsq2theta), "L");
  legend1->AddEntry(gr_surv_prob_decay_e, Form("Decay: #alpha=%.5f eV^{2}, sin^{2}2#theta=%.2f",input_alpha_decay,input_sinsq2theta), "L");
  legend1->AddEntry(gr_surv_prob_decoh_e, Form("Decoherence: #mu^{2}=%.5f eV^{2}, sin^{2}2#theta=%.2f",input_musq_decoh,input_sinsq2theta), "L");
  legend1->Draw();
  c1->Update();

  TCanvas * c2 = new TCanvas("c2","",120,120,600,600);
  TLegend * legend2 = new TLegend(0.3, 0.7, 0.85, 0.9);
  legend2->SetFillColor(0);
  legend2->SetBorderSize(0);
  gr_surv_prob_decay_le->SetLineStyle(kDotted);
  gr_surv_prob_decoh_le->SetLineStyle(kDashed);
  gr_surv_prob_disap_le->Draw("al");
  gr_surv_prob_decay_le->Draw("l");
  gr_surv_prob_decoh_le->Draw("l");
  legend2->AddEntry(gr_surv_prob_disap_le, Form("Oscillations: #Deltam^{2}=%.5f eV^{2}, sin^{2}2#theta=%.2f",input_dmsq,input_sinsq2theta), "L");
  legend2->AddEntry(gr_surv_prob_decay_le, Form("Decay: #alpha=%.5f eV^{2}, sin^{2}2#theta=%.2f",input_alpha_decay,input_sinsq2theta), "L");
  legend2->AddEntry(gr_surv_prob_decoh_le, Form("Decoherence: #mu^{2}=%.5f eV^{2}, sin^{2}2#theta=%.2f",input_musq_decoh,input_sinsq2theta), "L");
  legend2->Draw();
  c2->Update();
}
