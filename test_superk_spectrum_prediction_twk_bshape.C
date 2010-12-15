//
// ------------------------------------------------------------------------------------------------
// Plot predicted SuperK reconstructed neutrino energy spectrum for various flux shape tweaks.
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

#include "Conventions.h"

void test_superk_spectrum_prediction_twk_bshape()
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

  const int nexp = 500; // number of 'flux tweaks'

  // oscillation params
  double input_sinsq2theta23 = 0.;
  double input_dm23sq        = 0.;

  // histogram booking
  TH1D * hEreco [nexp];
  TH1D * hErecoMax = 0;
  TH1D * hErecoMin = 0;

  // superk spectrum predictor
  SuperkExpectation skexp;
  skexp.Read("inputs/analysis_configuration.xml");

  // phys params
  PhysParamList * fparm = PhysParamList::NumuDisapFitNoSyst();
  fparm->Update(t2kana::yr1numu::conv::kPhypSinsq2Theta23, input_sinsq2theta23);
  fparm->Update(t2kana::yr1numu::conv::kPhypDm23sq,        input_dm23sq);

  for(int iexp = 0; iexp < nexp; iexp++) {

     // tweak spectrum generation model
     cout << "*** twk.nu: " << iexp << endl;
     if(iexp>0) {
       // use 0^th entry to store the untweaked SK reco energy spectrum
       skexp.GetBeamInputs()->TweakShape();       
     }

     // generate predictions
     skexp.Generate(*fparm);

     // store spectrum
     hEreco[iexp] = new TH1D(*skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike));
     hEreco[iexp]->SetDirectory(0);
  }

  // find min/max and calculate fractional deviation from nominal
  TH1D * hErecoMax = new TH1D(*hEreco[0]);
  TH1D * hErecoMin = new TH1D(*hEreco[0]);
  hErecoMax->SetDirectory(0);
  hErecoMin->SetDirectory(0);
  for(int iexp = 1; iexp < nexp; iexp++) {
    for(int ibin=1; ibin <= conv::kNvErecoBins; ibin++) {
       double N = hEreco[iexp]->GetBinContent(ibin);
       hErecoMax->SetBinContent(ibin, TMath::Max(N, hErecoMax->GetBinContent(ibin)));
       hErecoMin->SetBinContent(ibin, TMath::Min(N, hErecoMin->GetBinContent(ibin)));
     }//ibin
  }//iexp
  for(int ibin=1; ibin <= conv::kNvErecoBins; ibin++) {
     double ndef = hEreco[0]->GetBinContent(ibin);
     double nmax = hErecoMax->GetBinContent(ibin);
     double nmin = hErecoMin->GetBinContent(ibin);
     double frac_max = (ndef > 0) ? (nmax-ndef)/ndef : 0;
     double frac_min = (ndef > 0) ? (nmin-ndef)/ndef : 0;
     hErecoMax->SetBinContent(ibin, frac_max);
     hErecoMin->SetBinContent(ibin, frac_min);
  }
  hErecoMax->GetYaxis()->SetRangeUser(-1,1);
  hErecoMin->GetYaxis()->SetRangeUser(-1,1);

  //
  // draw
  //
  TCanvas * c = new TCanvas("c","",20,20,500,500);
  gStyle->SetOptStat(0);
  c->Divide(1,2);
  c->cd(1);
  for(int iexp = 0; iexp < nexp; iexp++) {
    if(iexp==0) {
      hEreco[iexp]->GetXaxis()->SetTitle("E_{reco} (GeV)");
      hEreco[iexp]->GetYaxis()->SetTitle("Number of events");
      hEreco[iexp]->SetLineWidth(2);
      hEreco[iexp]->Draw();
    } else {
      hEreco[iexp]->Draw("same");
   }
  }//iexp
  c->cd(2);
  hErecoMax->Draw();
  hErecoMin->Draw("same");
  c->Update();
}
