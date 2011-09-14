//
// ------------------------------------------------------------------------------------------------
// Plot contributions from various event modes in the predicted SuperK spectrum.
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

#include <iomanip.h>
#include "Conventions.h"

void test_superk_spectrum_prediction_contributions
(
  const char * hypothesis              = "3f_osc",                       /* options: "2f_numu_disap", "3f_osc" */
  double       input_sinsq2theta12     = 0.86,                           /* */
  double       input_sinsq2theta13     = 0.10,                           /* */
  double       input_sinsq2theta23     = 1.00,                           /* */
  double       input_deltacp           = 0.00,                           /* */
  double       input_dm12sq            = 7.58E-5,                        /* */
  double       input_dm23sq            = 2.32E-3,                        /* */
  const char * config_file             = "analysis_configuration.xml",   /* */
  bool         print_spectrum          = false                           /* */
)
{
  using namespace std;
  using namespace t2kana;

  // define list of contributing modes
  const unsigned int nm = 15;
  conv::SuperkMode_t modes [nm] = {
     conv::kSKmNuMuCCQE, 
     conv::kSKmNuMuCC1pi,
     conv::kSKmNuMuCCoth, 
     conv::kSKmNC, 
     conv::kSKmNuMuBarCCQE, 
     conv::kSKmNuMuBarCC1pi, 
     conv::kSKmNuMuBarCCoth, 
     conv::kSKmNuMuBarNC, 
     conv::kSKmNuECCQE,
     conv::kSKmNuECC1pi,
     conv::kSKmNuECCoth,
     conv::kSKmNuENC,
     conv::kSKmSigNuECCQE,
     conv::kSKmSigNuECC1pi,
     conv::kSKmSigNuECCoth
  };

  // histogram booking
  TH1D * hErecoTot = new TH1D("","",conv::kNvErecoBins,conv::kNvErecoBinEdges);
  TH1D * hErecoMode[nm]; //< contribution from each mode  

  // set physics hypothesis
  PhysHypothesis ph(hypothesis);
  if(ph.Hypothesis() == conv::kHypNuMuDisapp2F) {
    ph.SetCurrValue(conv::kPhypSinsq2Theta23,  input_sinsq2theta23);
    ph.SetCurrValue(conv::kPhypDm23sq,         input_dm23sq);
  } else
  if(ph.Hypothesis() == conv::kHypOsc3F) {
    ph.SetCurrValue(conv::kPhypSinsq2Theta12,  input_sinsq2theta12);
    ph.SetCurrValue(conv::kPhypSinsq2Theta13,  input_sinsq2theta13);
    ph.SetCurrValue(conv::kPhypSinsq2Theta23,  input_sinsq2theta23);
    ph.SetCurrValue(conv::kPhypDeltaCP,        input_deltacp);
    ph.SetCurrValue(conv::kPhypDm12sq,         input_dm12sq);
    ph.SetCurrValue(conv::kPhypDm23sq,         input_dm23sq);
  }
  cout << "Physics hypothesis set: " << endl;
  cout << ph;

  // create & config superk spectrum predictor
  const char * config = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);
  cout << "Reading analysis configuration from: " << config << endl;
  SuperkExpectation skexp;
  skexp.Read(config);

  // generate prediction
  skexp.Generate(ph);

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
  hErecoMode[conv::kSKmNuMuBarCC1pi]   -> SetFillStyle ( 3345     );
  hErecoMode[conv::kSKmNuMuBarCCoth]   -> SetFillStyle ( 3354     );
  hErecoMode[conv::kSKmNuMuBarNC]      -> SetFillStyle ( 3006     );
  hErecoMode[conv::kSKmNuECCQE]        -> SetFillStyle ( 3001     );
  hErecoMode[conv::kSKmNuECC1pi]       -> SetFillStyle ( 3345     );
  hErecoMode[conv::kSKmNuECCoth]       -> SetFillStyle ( 3354     );
  hErecoMode[conv::kSKmNuENC]          -> SetFillStyle ( 3006     );
  hErecoMode[conv::kSKmSigNuECCQE]     -> SetFillStyle ( 3001     );
  hErecoMode[conv::kSKmSigNuECC1pi]    -> SetFillStyle ( 3345     );
  hErecoMode[conv::kSKmSigNuECCoth]    -> SetFillStyle ( 3354     );
  hErecoMode[conv::kSKmNuMuCCQE]       -> SetFillColor ( kWhite   );
  hErecoMode[conv::kSKmNuMuCC1pi]      -> SetFillColor ( kBlack   );
  hErecoMode[conv::kSKmNuMuCCoth]      -> SetFillColor ( kBlack   );
  hErecoMode[conv::kSKmNC]             -> SetFillColor ( kRed     );
  hErecoMode[conv::kSKmNuMuBarCCQE]    -> SetFillColor ( kAzure-3 );
  hErecoMode[conv::kSKmNuMuBarCC1pi]   -> SetFillColor ( kAzure-3 );
  hErecoMode[conv::kSKmNuMuBarCCoth]   -> SetFillColor ( kAzure-3 );
  hErecoMode[conv::kSKmNuMuBarNC]      -> SetFillColor ( kAzure-3 );
  hErecoMode[conv::kSKmNuECCQE]        -> SetFillColor ( kGreen+1 );
  hErecoMode[conv::kSKmNuECC1pi]       -> SetFillColor ( kGreen+1 );
  hErecoMode[conv::kSKmNuECCoth]       -> SetFillColor ( kGreen+1 );
  hErecoMode[conv::kSKmNuENC]          -> SetFillColor ( kGreen+1 );
  hErecoMode[conv::kSKmSigNuECCQE]     -> SetFillColor ( kGreen+3 );
  hErecoMode[conv::kSKmSigNuECC1pi]    -> SetFillColor ( kGreen+3 );
  hErecoMode[conv::kSKmSigNuECCoth]    -> SetFillColor ( kGreen+3 );
  hErecoMode[conv::kSKmNuMuCCQE]       -> SetLineColor ( kBlack   );
  hErecoMode[conv::kSKmNuMuCC1pi]      -> SetLineWidth ( 1        );
  hErecoMode[conv::kSKmNuMuCCoth]      -> SetLineWidth ( 1        );
  hErecoMode[conv::kSKmNC]             -> SetLineWidth ( 1        );
  hErecoMode[conv::kSKmNuMuCCQE]       -> SetLineWidth ( 1        );
  hErecoMode[conv::kSKmNuMuBarCCQE]    -> SetLineWidth ( 1        );
  hErecoMode[conv::kSKmNuMuBarCC1pi]   -> SetLineWidth ( 1        );
  hErecoMode[conv::kSKmNuMuBarCCoth]   -> SetLineWidth ( 1        );
  hErecoMode[conv::kSKmNuMuBarNC]      -> SetLineWidth ( 1        );
  hErecoMode[conv::kSKmNuECCQE]        -> SetLineWidth ( 1        );
  hErecoMode[conv::kSKmNuECC1pi]       -> SetLineWidth ( 1        );
  hErecoMode[conv::kSKmNuECCoth]       -> SetLineWidth ( 1        );
  hErecoMode[conv::kSKmNuENC]          -> SetLineWidth ( 1        );
  hErecoMode[conv::kSKmSigNuECCQE]     -> SetLineWidth ( 1        );
  hErecoMode[conv::kSKmSigNuECC1pi]    -> SetLineWidth ( 1        );
  hErecoMode[conv::kSKmSigNuECCoth]    -> SetLineWidth ( 1        );

  TCanvas * c = new TCanvas("c","",20,20,500,500);
  TLegend * legend = new TLegend(0.3, 0.7, 0.85, 0.9);
  legend->SetFillColor(0);
  legend->SetBorderSize(0);
  hErecoTot->Draw();
  hErecoTot->GetXaxis()->SetTitle("E_{reco} (GeV)");
  hErecoTot->GetYaxis()->SetTitle("Number of events");
  THStack * stack = new THStack("stack","spectrum contributions");
  stack->Add(hErecoMode[conv::kSKmSigNuECC1pi]);
  stack->Add(hErecoMode[conv::kSKmSigNuECCoth]);
  stack->Add(hErecoMode[conv::kSKmSigNuECCQE]);
  stack->Add(hErecoMode[conv::kSKmNuENC]);
  stack->Add(hErecoMode[conv::kSKmNuECC1pi]);
  stack->Add(hErecoMode[conv::kSKmNuECCoth]);
  stack->Add(hErecoMode[conv::kSKmNuECCQE]);
  stack->Add(hErecoMode[conv::kSKmNuMuBarNC]);
  stack->Add(hErecoMode[conv::kSKmNuMuBarCC1pi]);
  stack->Add(hErecoMode[conv::kSKmNuMuBarCCoth]);
  stack->Add(hErecoMode[conv::kSKmNuMuBarCCQE]);
  stack->Add(hErecoMode[conv::kSKmNC]);
  stack->Add(hErecoMode[conv::kSKmNuMuCCoth]);
  stack->Add(hErecoMode[conv::kSKmNuMuCC1pi]);
  stack->Add(hErecoMode[conv::kSKmNuMuCCQE]);
  stack->Draw("same");
  legend->AddEntry(hErecoTot,                            "total",                                "L" );
  legend->AddEntry(hErecoMode[conv::kSKmNuMuCCQE],       "#nu_{#mu} CCQE",                       "LF");
  legend->AddEntry(hErecoMode[conv::kSKmNuMuCC1pi],      "#nu_{#mu} CC1#pi",                     "LF");
  legend->AddEntry(hErecoMode[conv::kSKmNuMuCCoth],      "#nu_{#mu} CC other",                   "LF");
  legend->AddEntry(hErecoMode[conv::kSKmNC],             "#nu_{#mu}(#nu_{#tau}) NC",             "LF");
  legend->AddEntry(hErecoMode[conv::kSKmNuMuBarCCQE],    "#bar{#nu_{#mu}} CCQE",                 "LF");
  legend->AddEntry(hErecoMode[conv::kSKmNuMuBarCC1pi],   "#bar{#nu_{#mu}} CC1#pi",               "LF");
  legend->AddEntry(hErecoMode[conv::kSKmNuMuBarCCoth],   "#bar{#nu_{#mu}} CC other",             "LF");
  legend->AddEntry(hErecoMode[conv::kSKmNuMuBarNC],      "#bar{#nu_{#mu}}(#bar{#nu_{#tau}}) NC", "LF");
  legend->AddEntry(hErecoMode[conv::kSKmNuECCQE],        "#nu_{e} CCQE",                         "LF");
  legend->AddEntry(hErecoMode[conv::kSKmNuECC1pi],       "#nu_{e} CC1#pi",                       "LF");
  legend->AddEntry(hErecoMode[conv::kSKmNuECCoth],       "#nu_{e} CC other",                     "LF");
  legend->AddEntry(hErecoMode[conv::kSKmNuENC],          "#nu_{e} NC",                           "LF");
  legend->AddEntry(hErecoMode[conv::kSKmSigNuECCQE],     "#nu_{e} (sig) CCQE",                   "LF");
  legend->AddEntry(hErecoMode[conv::kSKmSigNuECC1pi],    "#nu_{e} (sig) CC1#pi",                 "LF");
  legend->AddEntry(hErecoMode[conv::kSKmSigNuECCoth],    "#nu_{e} (sig) CC other",               "LF");
  legend->Draw();
  c->Update();
  cout << fixed;
  double ntot = hErecoTot->Integral();
  cout << "Total number of events: " 
       << setw(10) << setfill(' ') << setprecision(6) << ntot << endl;
  const int nslices=6;
  double slice_max[nslices] = { 0.50, 0.75, 1.00, 1.25, 1.50, 2.00 };
  for(int islice=0; islice<nslices; islice++) {
    double esmax = slice_max[islice];
    cout << "Total number of events (< " << setw(8) << setfill(' ') << setprecision(4) << esmax << " GeV): " 
         << setw(10) << setfill(' ') << setprecision(6) << hErecoTot->Integral(1,hErecoTot->FindBin(esmax)) << endl;
  }
  for(int im = 0; im < nm; im++) {
     double nmode = hErecoMode[im]->Integral();
     cout << "Contrib from mode: " 
          << setw(20) << setfill(' ') << util::AsString(modes[im]) << " : " 
          << setw(15) << setfill(' ') << setprecision(6) << nmode << " events ("
          << setw(10) << setfill(' ') << setprecision(6) << 100*nmode/ntot << " %)" << endl;
  }

  //
  // save histograms in a ROOT file
  //
  string oscparams = util::AsString(hypothesis,
     input_sinsq2theta12,input_sinsq2theta13,input_sinsq2theta23,input_deltacp,input_dm12sq,input_dm23sq);
  TFile f(Form("sk_ereco_spectrum-%s.root",oscparams.c_str()),"recreate");
  hErecoTot                            -> Write("total");
  hErecoMode[conv::kSKmNuMuCCQE]       -> Write("numu_ccqe");
  hErecoMode[conv::kSKmNuMuCC1pi]      -> Write("numu_cc1pi"); 
  hErecoMode[conv::kSKmNuMuCCoth]      -> Write("numu_ccoth");
  hErecoMode[conv::kSKmNC]             -> Write("numu_and_nutau_nc");
  hErecoMode[conv::kSKmNuMuBarCCQE]    -> Write("numubar_ccqe");
  hErecoMode[conv::kSKmNuMuBarCC1pi]   -> Write("numubar_cc1pi");
  hErecoMode[conv::kSKmNuMuBarCCoth]   -> Write("numubar_ccoth");
  hErecoMode[conv::kSKmNuMuBarNC]      -> Write("numubar_and_nutaubar_nc");
  hErecoMode[conv::kSKmNuECCQE]        -> Write("nue_ccqe");
  hErecoMode[conv::kSKmNuECC1pi]       -> Write("nue_cc1pi");
  hErecoMode[conv::kSKmNuECCoth]       -> Write("nue_ccoth");
  hErecoMode[conv::kSKmNuENC]          -> Write("nue_nc");
  hErecoMode[conv::kSKmSigNuECCQE]     -> Write("nuesig_ccqe");
  hErecoMode[conv::kSKmSigNuECC1pi]    -> Write("nuesig_cc1pi");
  hErecoMode[conv::kSKmSigNuECCoth]    -> Write("nuesig_ccoth");
  f.Close();

  //
  // print-out numbers requested during the analysis review
  //

  if(print_spectrum) {  
    int nprec=5;
    cout << fixed;
    for(int ibin = 1; ibin <= hErecoTot->GetXaxis()->GetNbins(); ibin++) {
      cout << setw(6) << setfill(' ') << ibin
           << " : "
           << setw(6) << setfill(' ') << setprecision(3) << hErecoTot->GetXaxis()->GetBinLowEdge(ibin)
           << " - "
           << setw(6) << setfill(' ') << setprecision(3) << hErecoTot->GetXaxis()->GetBinUpEdge(ibin)
           << " : "
           << setw(6) << setfill(' ') << setprecision(nprec) << hErecoTot->GetBinContent(ibin)
           << endl;
    }//ibin
  }

}
