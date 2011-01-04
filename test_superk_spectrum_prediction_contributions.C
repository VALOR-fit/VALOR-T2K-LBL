//
// ------------------------------------------------------------------------------------------------
// Plot contributions from various event modes in the predicted SuperK spectrum.
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

#include <iomanip.h>
#include "Conventions.h"

void test_superk_spectrum_prediction_contributions()
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
  double input_sinsq2theta23 = 0.00;
  double input_dm23sq        = 2.4E-3;

  const unsigned int nm = 10;
  conv::SuperkMode_t modes [nm] = {
     conv::kSKmNuMuCCQE, 
     conv::kSKmNuMuCC1pi,
     conv::kSKmNuMuCCoth, 
     conv::kSKmNC, 
     conv::kSKmNuMuBarCCQE, 
     conv::kSKmNuMuBarCCnonQE, 
     conv::kSKmNuMuBarNC, 
     conv::kSKmNuECCQE,
     conv::kSKmNuECCnonQE,
     conv::kSKmNuENC
  };

  // histogram booking
  TH1D * hErecoTot = new TH1D("","",conv::kNvErecoBins,conv::kNvErecoBinEdges);
  TH1D * hErecoMode[nm]; //< contribution from each mode  

  // superk spectrum predictor
  SuperkExpectation skexp;
  skexp.Read("inputs/analysis_configuration.xml");

  // fit params
  PhysParamList * fparm = PhysParamList::NumuDisapFitNoSyst();
  fparm->Update(t2kana::yr1numu::conv::kPhypSinsq2Theta23, input_sinsq2theta23);
  fparm->Update(t2kana::yr1numu::conv::kPhypDm23sq,        input_dm23sq);

  // generate prediction
  skexp.Generate(*fparm);

  // get energy spectra
  for(int im = 0; im < nm; im++) {
    hErecoMode[im] = new TH1D(*skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike, modes[im]));
  }
  hErecoTot = new TH1D(*skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike));

  //
  // draw
  //

  gStyle->SetOptStat(0);

  hErecoTot                            -> SetLineWidth ( 2        );
  hErecoMode[conv::kSKmNuMuCC1pi]      -> SetFillStyle ( 3544     );
  hErecoMode[conv::kSKmNuMuCCoth]      -> SetFillStyle ( 3244     );
  hErecoMode[conv::kSKmNC]             -> SetFillStyle ( 3006     );
  hErecoMode[conv::kSKmNuMuBarCCQE]    -> SetFillStyle ( 3001     );
  hErecoMode[conv::kSKmNuMuBarCCnonQE] -> SetFillStyle ( 3345     );
  hErecoMode[conv::kSKmNuMuBarNC]      -> SetFillStyle ( 3006     );
  hErecoMode[conv::kSKmNuECCQE]        -> SetFillStyle ( 3001     );
  hErecoMode[conv::kSKmNuECCnonQE]     -> SetFillStyle ( 3345     );
  hErecoMode[conv::kSKmNuENC]          -> SetFillStyle ( 3006     );
  hErecoMode[conv::kSKmNuMuCC1pi]      -> SetFillColor ( kBlack   );
  hErecoMode[conv::kSKmNuMuCCoth]      -> SetFillColor ( kBlack   );
  hErecoMode[conv::kSKmNC]             -> SetFillColor ( kRed     );
  hErecoMode[conv::kSKmNuMuBarCCQE]    -> SetFillColor ( kAzure-3 );
  hErecoMode[conv::kSKmNuMuBarCCnonQE] -> SetFillColor ( kAzure-3 );
  hErecoMode[conv::kSKmNuMuBarNC]      -> SetFillColor ( kAzure-3 );
  hErecoMode[conv::kSKmNuECCQE]        -> SetFillColor ( kGreen+1 );
  hErecoMode[conv::kSKmNuECCnonQE]     -> SetFillColor ( kGreen+1 );
  hErecoMode[conv::kSKmNuENC]          -> SetFillColor ( kGreen+1 );
  hErecoMode[conv::kSKmNuMuCCQE]       -> SetLineColor ( kBlack   );
  hErecoMode[conv::kSKmNuMuCC1pi]      -> SetLineWidth ( 1        );
  hErecoMode[conv::kSKmNuMuCCoth]      -> SetLineWidth ( 1        );
  hErecoMode[conv::kSKmNC]             -> SetLineWidth ( 1        );
  hErecoMode[conv::kSKmNuMuCCQE]       -> SetLineWidth ( 1        );
  hErecoMode[conv::kSKmNuMuBarCCQE]    -> SetLineWidth ( 1        );
  hErecoMode[conv::kSKmNuMuBarCCnonQE] -> SetLineWidth ( 1        );
  hErecoMode[conv::kSKmNuMuBarNC]      -> SetLineWidth ( 1        );
  hErecoMode[conv::kSKmNuECCQE]        -> SetLineWidth ( 1        );
  hErecoMode[conv::kSKmNuECCnonQE]     -> SetLineWidth ( 1        );
  hErecoMode[conv::kSKmNuENC]          -> SetLineWidth ( 1        );

  TCanvas * c = new TCanvas("c","",20,20,500,500);
  TLegend * legend = new TLegend(0.3, 0.7, 0.85, 0.9);
  legend->SetFillColor(0);
  legend->SetBorderSize(0);
  hErecoTot->Draw();
  hErecoTot->GetXaxis()->SetTitle("E_{reco} (GeV)");
  hErecoTot->GetYaxis()->SetTitle("Number of events");
  THStack * stack = new THStack("stack","spectrum contributions");
  stack->Add(hErecoMode[conv::kSKmNuENC]);
  stack->Add(hErecoMode[conv::kSKmNuECCnonQE]);
  stack->Add(hErecoMode[conv::kSKmNuECCQE]);
  stack->Add(hErecoMode[conv::kSKmNuMuBarNC]);
  stack->Add(hErecoMode[conv::kSKmNuMuBarCCnonQE]);
  stack->Add(hErecoMode[conv::kSKmNuMuBarCCQE]);
  stack->Add(hErecoMode[conv::kSKmNC]);
  stack->Add(hErecoMode[conv::kSKmNuMuCCoth]);
  stack->Add(hErecoMode[conv::kSKmNuMuCC1pi]);
  stack->Add(hErecoMode[conv::kSKmNuMuCCQE]);
  stack->Draw("same");
  legend->AddEntry(hErecoTot,                            "total",                                "L" );
  legend->AddEntry(hErecoMode[conv::kSKmNuMuCCQE],       "#nu_{#mu} CCQE",                       "L" );
  legend->AddEntry(hErecoMode[conv::kSKmNuMuCC1pi],      "#nu_{#mu} CC1#pi",                     "LF");
  legend->AddEntry(hErecoMode[conv::kSKmNuMuCCoth],      "#nu_{#mu} CC other",                   "LF");
  legend->AddEntry(hErecoMode[conv::kSKmNC],             "#nu_{#mu}(#nu_{#tau}) NC",             "LF");
  legend->AddEntry(hErecoMode[conv::kSKmNuMuBarCCQE],    "#bar{#nu_{#mu}} CCQE",                 "LF");
  legend->AddEntry(hErecoMode[conv::kSKmNuMuBarCCnonQE], "#bar{#nu_{#mu}} CCnonQE",              "LF");
  legend->AddEntry(hErecoMode[conv::kSKmNuMuBarNC],      "#bar{#nu_{#mu}}(#bar{#nu_{#tau}}) NC", "LF");
  legend->AddEntry(hErecoMode[conv::kSKmNuECCQE],        "#nu_{e} CCQE",                         "LF");
  legend->AddEntry(hErecoMode[conv::kSKmNuECCnonQE],     "#nu_{e} CCnonQE",                      "LF");
  legend->AddEntry(hErecoMode[conv::kSKmNuENC],          "#nu_{e} NC",                           "LF");
  legend->Draw();
  c->Update();
  double E0 = 0.50; // GeV
  double E1 = 0.75; // GeV
  double E2 = 1.00; // GeV
  double E3 = 2.00; // GeV
  double ntot = hErecoTot->Integral();
  cout << "Total number of events: " 
       << setw(10) << setfill(' ') << setprecision(6) << ntot << endl;
  cout << "Total number of events (< " << setw(8) << setfill(' ') << setprecision(4) << E0 << " GeV): " 
       << setw(10) << setfill(' ') << setprecision(6) << hErecoTot->Integral(1,hErecoTot->FindBin(E0)) << endl;
  cout << "Total number of events (< " << setw(8) << setfill(' ') << setprecision(4) << E1 << " GeV): " 
       << setw(10) << setfill(' ') << setprecision(6) << hErecoTot->Integral(1,hErecoTot->FindBin(E1)) << endl;
  cout << "Total number of events (< " << setw(8) << setfill(' ') << setprecision(4) << E2 << " GeV): " 
       << setw(10) << setfill(' ') << setprecision(6) << hErecoTot->Integral(1,hErecoTot->FindBin(E2)) << endl;
  cout << "Total number of events (< " << setw(8) << setfill(' ') << setprecision(4) << E3 << " GeV): " 
       << setw(10) << setfill(' ') << setprecision(6) << hErecoTot->Integral(1,hErecoTot->FindBin(E3)) << endl;
  for(int im = 0; im < nm; im++) {
     double nmode = hErecoMode[im]->Integral();
     cout << "Contrib from mode: " 
          << setw(20) << setfill(' ') << util::AsString(modes[im]) << " : " 
          << setw(15) << setfill(' ') << setprecision(6) << nmode << " events ("
          << setw(10) << setfill(' ') << setprecision(6) << 100*nmode/ntot << " %)" << endl;
  }

  TFile f(Form("sk_ereco_spectrum_components-sinsq2theta23_%.2f-dm23sq_%.5feV2.root",input_sinsq2theta23,input_dm23sq),"recreate");
  hErecoTot                            -> Write("total");
  hErecoMode[conv::kSKmNuMuCCQE]       -> Write("numu_ccqe");
  hErecoMode[conv::kSKmNuMuCC1pi]      -> Write("numu_cc1pi"); 
  hErecoMode[conv::kSKmNuMuCCoth]      -> Write("numu_ccoth");
  hErecoMode[conv::kSKmNC]             -> Write("numu_and_nutau_nc");
  hErecoMode[conv::kSKmNuMuBarCCQE]    -> Write("numubar_ccqe");
  hErecoMode[conv::kSKmNuMuBarCCnonQE] -> Write("numubar_ccnonqe");
  hErecoMode[conv::kSKmNuMuBarNC]      -> Write("numubar_and_nutaubar_nc");
  hErecoMode[conv::kSKmNuECCQE]        -> Write("nue_ccqe");
  hErecoMode[conv::kSKmNuECCnonQE]     -> Write("nue_ccnonqe");
  hErecoMode[conv::kSKmNuENC]          -> Write("nue_nc");
  f.Close();
}
