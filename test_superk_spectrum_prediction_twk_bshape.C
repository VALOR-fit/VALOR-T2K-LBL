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

  // histogram booking
  TH1D * hEreco [nexp];
  for(int iexp = 0; iexp < nexp; iexp++) {
     hEreco[iexp] = new TH1D("","",conv::kNErecoBins,conv::kErecoMin,conv::kErecoMax);
     hEreco[iexp]->SetDirectory(0);
  }

  // superk spectrum predictor
  SuperkExpectation skexp;
  skexp.Read("inputs/analysis_configuration.xml");

  // phys params
  PhysParamList * fparm = PhysParamList::NumuDisapFitNoSyst();
  fparm->Update(t2kana::yr1numu::conv::kPhypSinsq2Theta23, 0.0);
  fparm->Update(t2kana::yr1numu::conv::kPhypDm23sq,        0.0);

  // generate predictions
  for(int iexp = 0; iexp < nexp; iexp++) {
      cout << "*** twk.nu: " << iexp << endl;
     if(iexp>0) {
       // use 0^th entry to store the untweaked SK reco energy spectrum
       skexp.GetBeamInputs()->TweakShape();       
     }
     for(int ibin=1; ibin <= conv::kNErecoBins; ibin++) {
        double Ereco = hEreco[iexp]->GetBinCenter(ibin);
        double Nexp  = skexp.Generate(Ereco, *fparm);
        //cout << "NumExpSK (Ereco = " << Ereco << ") = " << Nexp << endl;
        hEreco[iexp]->SetBinContent(ibin,Nexp);
     }//ibin
  }//iexp

  // draw
  TCanvas * c = new TCanvas("c","",20,20,500,500);
  gStyle->SetOptStat(0);
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
  c->Update();
}
