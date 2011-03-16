//
// ------------------------------------------------------------------------------------------------
// Tset change in number of ND280 events as we tweak the ND280 efficiency
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

#include <iomanip.h>

#include "Conventions.h"

void test_nd280_nevents_twk_eff()
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

  gSystem->Load("libT2KNuOsc");

  using namespace std;
  using namespace t2kana;

  const int ntwkdials = 5;
  double twkdial[ntwkdials] = { -1.0, -0.5, 0.0, +0.5, +1.0 };

  // histogram drawing options
  int hwidth [ntwkdials] = { 1,       1,       2,       1,       1      };
  int hcolor [ntwkdials] = { kBlue,   kBlue,   kBlack,  kRed,    kRed   };
  int hstyle [ntwkdials] = { kDotted, kDashed, kSolid,  kDashed, kDotted};

  double flux_correction_factor[ntwkdials];

  // superk spectrum predictor
  SuperkExpectation skexp;
  skexp.Read("inputs/analysis_configuration.xml");

  for(int idial = 0; idial < ntwkdials; idial++) {

     skexp.GetND280Measurement()->TweakEfficiency(twkdial[idial]);

     double nobscclike_nd280  = skexp.GetND280Measurement()->NObsCCLike();
     double nexpncbkg_nd280   = skexp.GetND280Measurement()->NExpNCBkg();
     double nexpcc_nd280      = skexp.GetND280Measurement()->NExpCC();

     flux_correction_factor[idial] = nobscclike_nd280 / (nexpncbkg_nd280 + nexpcc_nd280);

     cout << " > twk_dial = " 
          << setw(6) << setfill(' ') << setprecision(3) << twkdial[idial] 
          << ", nobscclike_nd280 = " 
          << setw(8) << setfill(' ') << setprecision(6) << nobscclike_nd280
          << ", nexpncbkg_nd280 = " 
          << setw(8) << setfill(' ') << setprecision(6) << nexpncbkg_nd280
          << ", nexpcc_nd280 = " 
          << setw(8) << setfill(' ') << setprecision(6) << nexpcc_nd280 
          << " ==> flux_correction_factor = " 
          << setw(8) << setfill(' ') << setprecision(6) << flux_correction_factor[idial] << endl;
  }

  TCanvas * c = new TCanvas("c","",20,20,500,500);
  TGraph * gr = new TGraph(ntwkdials, twkdial, flux_correction_factor);
  gr->SetMarkerStyle(20);
  gr->SetMarkerSize(1);
  gr->Draw("alp");
  c->Update();
}
