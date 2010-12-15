//
// ------------------------------------------------------------------------------------------------
// Plot predicted SuperK reconstructed neutrino energy spectrum for various tweaks of
// nd280 efficiency curves
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

#include "Conventions.h"

void test_superk_spectrum_prediction_twk_nd280eff()
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

  // oscillation params
  double input_sinsq2theta23 = 1.0;
  double input_dm23sq        = 2.4E-3;

  bool correlated = true;
  
  // n tweaks
  const int n = 1000;

  // superk spectrum predictor
  SuperkExpectation skexp;
  skexp.Read("inputs/analysis_configuration.xml");

  // phys params
  PhysParamList * fparm = PhysParamList::NumuDisapFitNoSyst();
  fparm->Update(t2kana::yr1numu::conv::kPhypSinsq2Theta23, input_sinsq2theta23);
  fparm->Update(t2kana::yr1numu::conv::kPhypDm23sq,        input_dm23sq);

  TRandom3 rand(0);

  //
  // generate predictions
  //

  TH1D * hEreco[n];

  for(int i = 0; i < n; i++) {
     cout << "Generating SuperK spectrum number: " << i << endl;
 
     // tweak spectrum generation model
     if(i>0) {
        // use i=0 as the default (untweaked)
        if(!correlated) {
          skexp.GetND280Inputs()->TweakEfficiency( conv::kNDmNuMuCCQE,     rand.Gaus(0,1) );
          skexp.GetND280Inputs()->TweakEfficiency( conv::kNDmNuMuCCRES,    rand.Gaus(0,1) );
          skexp.GetND280Inputs()->TweakEfficiency( conv::kNDmNuMuCCCOH,    rand.Gaus(0,1) );
          skexp.GetND280Inputs()->TweakEfficiency( conv::kNDmNuMuCCOth,    rand.Gaus(0,1) );
          skexp.GetND280Inputs()->TweakEfficiency( conv::kNDmNuMuNCEL,     rand.Gaus(0,1) );
          skexp.GetND280Inputs()->TweakEfficiency( conv::kNDmNuMuNCRES,    rand.Gaus(0,1) );
          skexp.GetND280Inputs()->TweakEfficiency( conv::kNDmNuMuNCCOH,    rand.Gaus(0,1) );
          skexp.GetND280Inputs()->TweakEfficiency( conv::kNDmNuMuNCOth,    rand.Gaus(0,1) );
          skexp.GetND280Inputs()->TweakEfficiency( conv::kNDmNuMuBarCCQE,  rand.Gaus(0,1) );
          skexp.GetND280Inputs()->TweakEfficiency( conv::kNDmNuMuBarCCRES, rand.Gaus(0,1) );
          skexp.GetND280Inputs()->TweakEfficiency( conv::kNDmNuMuBarCCCOH, rand.Gaus(0,1) );
          skexp.GetND280Inputs()->TweakEfficiency( conv::kNDmNuMuBarCCOth, rand.Gaus(0,1) );
          skexp.GetND280Inputs()->TweakEfficiency( conv::kNDmNuMuBarNCEL,  rand.Gaus(0,1) );
          skexp.GetND280Inputs()->TweakEfficiency( conv::kNDmNuMuBarNCRES, rand.Gaus(0,1) );
          skexp.GetND280Inputs()->TweakEfficiency( conv::kNDmNuMuBarNCCOH, rand.Gaus(0,1) );
          skexp.GetND280Inputs()->TweakEfficiency( conv::kNDmNuMuBarNCOth, rand.Gaus(0,1) );
          skexp.GetND280Inputs()->TweakEfficiency( conv::kNDmNuECCQE,      rand.Gaus(0,1) );
          skexp.GetND280Inputs()->TweakEfficiency( conv::kNDmNuECCRES,     rand.Gaus(0,1) );
          skexp.GetND280Inputs()->TweakEfficiency( conv::kNDmNuECCCOH,     rand.Gaus(0,1) );
          skexp.GetND280Inputs()->TweakEfficiency( conv::kNDmNuECCOth,     rand.Gaus(0,1) );
          skexp.GetND280Inputs()->TweakEfficiency( conv::kNDmNuENCEL,      rand.Gaus(0,1) );
          skexp.GetND280Inputs()->TweakEfficiency( conv::kNDmNuENCRES,     rand.Gaus(0,1) );
          skexp.GetND280Inputs()->TweakEfficiency( conv::kNDmNuENCCOH,     rand.Gaus(0,1) );
          skexp.GetND280Inputs()->TweakEfficiency( conv::kNDmNuENCOth,     rand.Gaus(0,1) );
       } else {
          double twkdial = (i==1) ? +1 : -1;
          skexp.GetND280Inputs()->TweakEfficiency(twkdial);
       }
     }

     // generate predictions  
     skexp.Generate(*fparm);
  
     // store spectra
     hEreco[i] = new TH1D(*skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike));
     hEreco[i]->SetDirectory(0);
  }//i

  // calculate err envelopes
  TH1D * alo = new TH1D("alo",  "spectrum err envelope: low  boundary",            conv::kNvErecoBins,conv::kNvErecoBinEdges);
  TH1D * ahi = new TH1D("ahi",  "spectrum err envelope: high boundary",            conv::kNvErecoBins,conv::kNvErecoBinEdges);
  TH1D * flo = new TH1D("flo",  "spectrum fractional err envelope: low  boundary", conv::kNvErecoBins,conv::kNvErecoBinEdges);
  TH1D * fhi = new TH1D("fhi",  "spectrum fractional err envelope: high boundary", conv::kNvErecoBins,conv::kNvErecoBinEdges);
  for(int ibin=1; ibin <= conv::kNvErecoBins; ibin++) {
    double nlo =  999999;
    double nhi = -999999;
    double def = hEreco[0]->GetBinContent(ibin);
    for(int i = 0; i < n; i++) {    
       double ncurr = hEreco[i]->GetBinContent(ibin);
       nlo = TMath::Min(nlo, ncurr);
       nhi = TMath::Max(nhi, ncurr);
    }
    alo->SetBinContent(ibin,nlo);
    ahi->SetBinContent(ibin,nhi);
    flo->SetBinContent(ibin, ((def>0) ? (nlo-def)/def : 0) );
    fhi->SetBinContent(ibin, ((def>0) ? (nhi-def)/def : 0) );
  }

  //
  // draw
  //
  TCanvas * c = new TCanvas("c","",20,20,500,500);
  gStyle->SetOptStat(0);
  c->Divide(1,2);
  c->cd(1);
  hEreco[0]->Draw();
  hEreco[0]->GetXaxis()->SetTitle("E_{reco} (GeV)");
  hEreco[0]->GetYaxis()->SetTitle("Number of events");
  alo->SetLineStyle(kDashed);
  ahi->SetLineStyle(kDashed);
  alo->Draw("same");
  ahi->Draw("same");
  c->cd(2);
  TH1F * hframe = (TH1F*) c->GetPad(2)->DrawFrame(0,-0.15,conv::kNvErecoBinEdges[conv::kNvErecoBins],0.15);
  hframe->GetXaxis()->SetTitle("E_{reco} (GeV)");
  hframe->Draw();
  flo->Draw("same");
  fhi->Draw("same");
}
