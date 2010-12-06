//
// ------------------------------------------------------------------------------------------------
// Plot predicted SuperK reconstructed neutrino energy spectrum for various sets of SuperK
// nuisance parameter variations that affect the SuperK efficiency.
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

#include "Conventions.h"

void test_superk_spectrum_prediction_twk_skeff_many()
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

  // histogram booking
  const int N = 300;
  TH1D * hEreco   [N];
  for(int i = 0; i < N; i++) {
     hEreco[i] = new TH1D("","",conv::kNvErecoBins,conv::kNvErecoBinEdges);
     hEreco[i]->SetDirectory(0);
  }
  TH1D * hErecoErrPos = new TH1D("hErecoErrPos","",conv::kNvErecoBins,conv::kNvErecoBinEdges);
  hErecoErrPos->SetDirectory(0);
  TH1D * hErecoErrNeg = new TH1D("hErecoErrNeg","",conv::kNvErecoBins,conv::kNvErecoBinEdges);
  hErecoErrNeg->SetDirectory(0);

  // superk spectrum predictor
  SuperkExpectation skexp;
  skexp.Read("inputs/analysis_configuration.xml");

  // phys params
  PhysParamList * fparm = PhysParamList::NumuDisapFitNoSyst();
  fparm->Update(t2kana::yr1numu::conv::kPhypSinsq2Theta23, 0);
  fparm->Update(t2kana::yr1numu::conv::kPhypDm23sq,        0);

  TRandom3 rand(0);

  // generate predictions
  double maxNexp = -1;
  for(int i = 0; i < N; i++) {
     double ccqe_eff_twk    = rand.Gaus(0.,1.);
     double ccnonqe_eff_twk = rand.Gaus(0.,1.);
     double nc_eff_twk      = rand.Gaus(0.,1.);
     skexp.GetSuperkInputs()->TweakEfficiency(
        t2kana::yr1numu::conv::kSKsFCFV1RingMuLike, t2kana::yr1numu::conv::kSKmNuMuCCQE, ccqe_eff_twk);
     skexp.GetSuperkInputs()->TweakEfficiency(
        t2kana::yr1numu::conv::kSKsFCFV1RingMuLike, t2kana::yr1numu::conv::kSKmNuMuCC1pi, ccnonqe_eff_twk);
     skexp.GetSuperkInputs()->TweakEfficiency(
        t2kana::yr1numu::conv::kSKsFCFV1RingMuLike, t2kana::yr1numu::conv::kSKmNC, nc_eff_twk);
     cout << " ** i = " << i << " : "
             << "ccqe eff twk = " << ccqe_eff_twk 
             << ", ccnonqe_eff_twk = " << ccnonqe_eff_twk 
             << ", nc_eff_twk = " << nc_eff_twk << endl;
     for(int ibin=1; ibin <= conv::kNvErecoBins; ibin++) {
        double Ereco = hEreco[i]->GetBinCenter(ibin);
        double Nexp  = skexp.Generate(Ereco, *fparm);
        maxNexp = TMath::Max(maxNexp,Nexp);
//      cout << "NumExpSK (Ereco = " << Ereco << ") = " << Nexp << endl;
        hEreco[i]->SetBinContent(ibin,Nexp);
     }//ibin
  }//i
  // get fractional err
  double global_maxerrpos = -1;
  double global_maxerrneg = +1;
  for(int ibin=1; ibin <= conv::kNvErecoBins; ibin++) {
    double max = -1.;
    double min = +1.;
    double ndeflt = hEreco[0]->GetBinContent(ibin);
    for(int i = 1; i < N; i++) {
      double ncurr  = hEreco[i]->GetBinContent(ibin);
      max = TMath::Max(max, ncurr);
      min = TMath::Min(min, ncurr);      
    }//i
    if(ndeflt>0) {
       double maxerrpos = (max-ndeflt)/ndeflt;
       double maxerrneg = (min-ndeflt)/ndeflt;
       global_maxerrpos = TMath::Max(maxerrpos, global_maxerrpos);
       global_maxerrneg = TMath::Min(maxerrneg, global_maxerrneg);
       hErecoErrPos->SetBinContent(ibin, maxerrpos);
       hErecoErrNeg->SetBinContent(ibin, maxerrneg);
    }
  }//ibin

  // draw
  TCanvas * c = new TCanvas("c","",20,20,900,500);
  c->Divide(1,2);
  gStyle->SetOptStat(0);
  c->cd(1);
  TH1F * hframe1 = (TH1F*) c->DrawFrame(
     0., 0., hEreco[0]->GetXaxis()->GetXmax(), 1.2*maxNexp);
  hframe1->GetXaxis()->SetTitle("E_{reco} (GeV)");                                
  hframe1->GetYaxis()->SetTitle("Number of events");
  hframe1->Draw();
  for(int i = 0; i < N; i++) {
    hEreco[i]->SetLineWidth( ((i==0) ? 2      : 1      ) );
    hEreco[i]->SetLineStyle( ((i==0) ? kSolid : kDotted) );
    hEreco[i]->Draw("same");
  }//i
  c->cd(2);
  TH1F * hframe2 = (TH1F*) c->DrawFrame(
     0., 1.2*global_maxerrneg, hEreco[0]->GetXaxis()->GetXmax(), 1.2*global_maxerrpos);
  hframe2->GetXaxis()->SetTitle("E_{reco} (GeV)");                                
  hframe2->GetYaxis()->SetTitle("Fractional variation");
  hframe2->Draw();
  hErecoErrPos->Draw("same");
  hErecoErrNeg->Draw("same");
  c->Update();
}
