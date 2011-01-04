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
  double input_sinsq2theta23 = 1.;
  double input_dm23sq        = 2.4E-3;

  // histogram booking
  TH1D * hEreco [nexp];

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
  for(int ibin=1; ibin <= conv::kNvErecoBins; ibin++) {
    double nom   = hEreco[0]->GetBinContent(ibin);
    double sumsq = 0;
    for(int iexp = 1; iexp < nexp; iexp++) {
       double N = hEreco[iexp]->GetBinContent(ibin);
       sumsq += TMath::Power(N-nom,2);	
    }//iexp
    double sig = TMath::Sqrt(sumsq/(nexp-1));
    double min = nom-sig;
    double max = nom+sig; 
    if(nom>0){
      hErecoMax -> SetBinContent (ibin,  (max-nom)/nom);
      hErecoMin -> SetBinContent (ibin,  (min-nom)/nom);
    }
  }//ibin
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
