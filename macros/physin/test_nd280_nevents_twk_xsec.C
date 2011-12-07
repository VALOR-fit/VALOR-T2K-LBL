//
// ------------------------------------------------------------------------------------------------
// Tset change in number of ND280 events as we tweak cross section systematic params
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

#include <iomanip.h>

//#include "Conventions.h"

void test_nd280_nevents_twk_xsec
(
   const char * systematic  = "cc1pi_to_ccqe",                          /* xsec systematic to tweak. allowed values are:*/
								       //cc1pi_to_ccqe ccoth_to_ccqe nc_to_ccqe nue_to_numu
   const char * config_file = "analysis_configuration.xml"  /* analysis configuration file */
)
{
  using namespace std;
  using namespace t2kana;

  const char * config = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);
  cout << "Reading analysis configuration from: " << config << endl;

  // superk spectrum predictor
  SuperkExpectation skexp;
  skexp.Read(config);

  TRandom3 rand(0);

  const int n = 1000;

  TH1D * hst_flux_correction_factor = new TH1D("hst_flux_correction_factor","",70,0.8,1.5);

  for(int i = 0; i < n; i++) {

     double twkdial = rand.Gaus(0,1);
     if(systematic == "cc1pi_to_ccqe")
       skexp.GetNuXSec()->TweakCC1piToCCQE(twkdial);
     else if(systematic == "ccoth_to_ccqe")
       skexp.GetNuXSec()->TweakCCOthToCCQE(twkdial);
     else if(systematic == "nc_to_ccqe")
       skexp.GetNuXSec()->TweakNCToCCQE   (twkdial);
     else if(systematic == "nue_to_numu")
       skexp.GetNuXSec()->TweakNuEToNuMu  (twkdial);
     else {
       cout << "\nUnknown systematic parameter " << systematic 
	    << ". Valid values are: cc1pi_to_ccqe ccoth_to_ccqe nc_to_ccqe nue_to_numu\n" << endl;
       return;
     }

     double nobscclike_nd280  = skexp.GetND280Measurement()->NObsCCLike();
     double nexpncbkg_nd280   = skexp.GetND280Measurement()->NExpNCBkg();
     double nexpcc_nd280      = skexp.GetND280Measurement()->NExpCC();

     flux_correction_factor = (nobscclike_nd280 - nexpncbkg_nd280) / nexpcc_nd280;

     hst_flux_correction_factor->Fill(flux_correction_factor);

     cout << " >> " << setw(8) << setfill(' ') << i 
          << ", twkdial = " 
          << setw(8) << setfill(' ') << setprecision(6) << twkdial
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
