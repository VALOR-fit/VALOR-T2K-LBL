//
// ------------------------------------------------------------------------------------------------
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

#include <iomanip.h>

#include "Conventions.h"

void test_nd280_nevents_twk_eff_many()
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

  // superk spectrum predictor
  SuperkExpectation skexp;
  skexp.Read("inputs/analysis_configuration.xml");

  TRandom3 rand(0);

  const int n = 10000;

  TH1D * hst_flux_correction_factor = new TH1D("hst_flux_correction_factor","",70,0.8,1.5);

  for(int i = 0; i < n; i++) {
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

     double nobscclike_nd280  = skexp.GetND280Inputs()->NObsCCLike();
     double nexpncbkg_nd280   = skexp.GetND280Inputs()->NExpNCBkg();
     double nexpcc_nd280      = skexp.GetND280Inputs()->NExpCC();

     flux_correction_factor = (nobscclike_nd280 - nexpncbkg_nd280) / nexpcc_nd280;

     hst_flux_correction_factor->Fill(flux_correction_factor);

     cout << " >> " << setw(8) << setfill(' ') << i
          << ", nobscclike_nd280 = " 
          << setw(8) << setfill(' ') << setprecision(6) << nobscclike_nd280
          << ", nexpncbkg_nd280 = " 
          << setw(8) << setfill(' ') << setprecision(6) << nexpncbkg_nd280
          << ", nexpcc_nd280 = " 
          << setw(8) << setfill(' ') << setprecision(6) << nexpcc_nd280 
          << " ==> flux_correction_factor = " 
          << setw(8) << setfill(' ') << setprecision(6) << flux_correction_factor
          << endl;
  }

  TCanvas * c = new TCanvas("c","",20,20,500,500);
  hst_flux_correction_factor->Draw();
  c->Update();
}
