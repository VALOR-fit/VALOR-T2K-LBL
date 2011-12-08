//
// ------------------------------------------------------------------------------------------------
// Plot contributions from various event modes in the predicted SuperK spectrum.
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

#include <iomanip.h>

using namespace std;
using namespace t2kana;

void test_sk_ereco_spectrum_nominal_nue
(
  conv::SuperkSample_t s               = conv::kSKsFCFV1RingELike,      /* */
  const char * hypothesis              = "3f_osc",                       /* options: "2f_numu_disap", "3f_osc" */
  double       input_mixing12          = 0.8704,                           /* */
  double       input_mixing13          = 0.10,                           /* */
  double       input_mixing23          = 1.00,                           /* */
  double       input_deltacp           = 0.00,                           /* */
  double       input_dm12sq            = 7.6E-5,                        /* */
  double       input_dm23sq            = 2.4E-3,                        /* */
  const char * config_file             = "analysis_configuration.xml",   /* */
  bool         print_spectrum          = false                           /* */
)
{

  int idx = util::FindSampleSlot(s);
  if(idx == -1){
    cerr << "Sample 1Re not defined" << endl;
    exit(1); 
  }
  // define list of contributing modes
  const unsigned int nm = 24;
  conv::SuperkMode_t modes [nm] = {
    conv::kSKmNuMuCCQE,
    conv::kSKmNuMuCC1pi,
    conv::kSKmNuMuCCcoh,
    conv::kSKmNuMuCCoth_e,
    conv::kSKmNuMuNC1pi0,
    conv::kSKmNuMuNCcoh,
    conv::kSKmNuMuNCoth_e,
    conv::kSKmNuMuBarCCQE,
    conv::kSKmNuMuBarCC1pi,
    conv::kSKmNuMuBarCCcoh,
    conv::kSKmNuMuBarCCoth_e,
    conv::kSKmNuMuBarNC1pi0,
    conv::kSKmNuMuBarNCcoh,
    conv::kSKmNuMuBarNCoth_e,
    conv::kSKmNuECCQE,
    conv::kSKmNuECC1pi,
    conv::kSKmNuECCcoh,
    conv::kSKmNuECCoth_e,
    conv::kSKmNuENC1pi0,
    conv::kSKmNuENCcoh,
    conv::kSKmNuENCoth_e,
    conv::kSKmSigNuECCQE,
    conv::kSKmSigNuECC1pi,
    conv::kSKmSigNuECCoth
  };

  // histogram booking
  TH1D * hErecoTot = new TH1D("","",conv::kNErecoBins1Re,conv::kErecoBinEdges1Re);
  TH1D * hErecoMode[conv::knm]; //< contribution from each mode  

  // set physics hypothesis
  PhysHypothesis ph(hypothesis);
  if(ph.Hypothesis() == conv::kHypNuMuDisapp2F) {
    ph.SetCurrValue(conv::kPhypMixing23,       input_mixing23);
    ph.SetCurrValue(conv::kPhypDm23sq,         input_dm23sq);
  } else
  if(ph.Hypothesis() == conv::kHypOsc3F) {
    ph.SetCurrValue(conv::kPhypMixing12,       input_mixing12);
    ph.SetCurrValue(conv::kPhypMixing13,       input_mixing13);
    ph.SetCurrValue(conv::kPhypMixing23,       input_mixing23);
    ph.SetCurrValue(conv::kPhypDeltaCP,        input_deltacp);
    ph.SetCurrValue(conv::kPhypDm12sq,         input_dm12sq);
    ph.SetCurrValue(conv::kPhypDm23sq,         input_dm23sq);
  }

  // create & config superk spectrum predictor
  const char * config = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);
  cout << "Reading analysis configuration from: " << config << endl;

  // Global options
  Opt::Instance()->gOscMixingParamChoice = conv::kNuMixSinsq2Theta;
  Opt::Instance()->gMassHierarchyAssumed = conv::kNmhNormal;

  SuperkExpectation skexp;
  skexp.Read(config);


  // Print options
  cout << "\n[Job inputs]" << endl;
  cout << "- Reading analysis configuration from: " << config << endl;
  cout << "- SuperK sample: " << util::AsString(s) << endl;
  cout << "- Physics hypothesis: " << ph << endl;
  cout << *Opt::Instance() << endl;


  // generate prediction
  skexp.Generate(ph);

  // get energy spectra
  for(int im = 0; im < nm; im++) {
    conv::SuperkMode_t m = modes[im];
    int midx = util::FindModeSlot(s,m);
    if(midx == -1) continue;
    hErecoMode[m] = new TH1D(*skexp.CurrPredictedErecoSpectrum(s, m));
  }
  hErecoTot = new TH1D(*skexp.CurrPredictedErecoSpectrum(s));

  //
  // draw
  //

  gStyle->SetOptStat(0);
  //  gStyle->SetFillColor(kWhite);


  hErecoTot                             -> SetLineWidth ( 2        );
  hErecoMode[conv::kSKmNuMuCC1pi]       -> SetFillStyle ( 3544     );
  hErecoMode[conv::kSKmNuMuCCcoh]    -> SetFillStyle ( 3244     );
  hErecoMode[conv::kSKmNuMuCCoth_e]    -> SetFillStyle ( 3244     );
  hErecoMode[conv::kSKmNuMuNC1pi0]          -> SetFillStyle ( 3006     );
  hErecoMode[conv::kSKmNuMuNCcoh]          -> SetFillStyle ( 3006     );
  hErecoMode[conv::kSKmNuMuNCoth_e]          -> SetFillStyle ( 3006     );
  hErecoMode[conv::kSKmNuMuBarCCQE]     -> SetFillStyle ( 3001     );
  hErecoMode[conv::kSKmNuMuBarCC1pi]    -> SetFillStyle ( 3345     );
  hErecoMode[conv::kSKmNuMuBarCCcoh] -> SetFillStyle ( 3354     );
  hErecoMode[conv::kSKmNuMuBarCCoth_e] -> SetFillStyle ( 3354     );
  hErecoMode[conv::kSKmNuMuBarNC1pi0]       -> SetFillStyle ( 3006     );
  hErecoMode[conv::kSKmNuMuBarNCcoh]       -> SetFillStyle ( 3006     );
  hErecoMode[conv::kSKmNuMuBarNCoth_e]       -> SetFillStyle ( 3006     );
  hErecoMode[conv::kSKmNuECCQE]         -> SetFillStyle ( 3001     );
  hErecoMode[conv::kSKmNuECC1pi]        -> SetFillStyle ( 3345     );
  hErecoMode[conv::kSKmNuECCcoh]     -> SetFillStyle ( 3354     );
  hErecoMode[conv::kSKmNuECCoth_e]     -> SetFillStyle ( 3354     );
  hErecoMode[conv::kSKmNuENC1pi0]           -> SetFillStyle ( 3006     );
  hErecoMode[conv::kSKmNuENCcoh]           -> SetFillStyle ( 3006     );
  hErecoMode[conv::kSKmNuENCoth_e]           -> SetFillStyle ( 3006     );
  hErecoMode[conv::kSKmSigNuECCQE]      -> SetFillStyle ( 3001     );
  hErecoMode[conv::kSKmSigNuECC1pi]     -> SetFillStyle ( 3345     );
  hErecoMode[conv::kSKmSigNuECCoth]     -> SetFillStyle ( 3354     );

  hErecoMode[conv::kSKmNuMuCCQE]        -> SetFillColor ( kWhite   );
  hErecoMode[conv::kSKmNuMuCC1pi]       -> SetFillColor ( kBlack   );
  hErecoMode[conv::kSKmNuMuCCcoh]    -> SetFillColor ( kGreen   );
  hErecoMode[conv::kSKmNuMuCCoth_e]    -> SetFillColor ( kBlue   );
  hErecoMode[conv::kSKmNuMuNC1pi0]          -> SetFillColor ( kRed     );
  hErecoMode[conv::kSKmNuMuNCcoh]          -> SetFillColor ( kRed -1     );  
  hErecoMode[conv::kSKmNuMuNCoth_e]          -> SetFillColor ( kRed -2    );
  hErecoMode[conv::kSKmNuMuBarCCQE]     -> SetFillColor ( kAzure-3 );
  hErecoMode[conv::kSKmNuMuBarCC1pi]    -> SetFillColor ( kAzure-4 );
  hErecoMode[conv::kSKmNuMuBarCCcoh] -> SetFillColor ( kAzure-6 );
  hErecoMode[conv::kSKmNuMuBarCCoth_e] -> SetFillColor ( 6 );
  hErecoMode[conv::kSKmNuMuBarNC1pi0]       -> SetFillColor ( kAzure-7 );
  hErecoMode[conv::kSKmNuMuBarNCcoh]       -> SetFillColor ( kAzure-9 );
  hErecoMode[conv::kSKmNuMuBarNCoth_e]       -> SetFillColor ( 6 );
  hErecoMode[conv::kSKmNuECCQE]         -> SetFillColor ( kGreen+1 );
  hErecoMode[conv::kSKmNuECC1pi]        -> SetFillColor ( kGreen+2 );
  hErecoMode[conv::kSKmNuECCcoh]     -> SetFillColor ( kGreen+2 );
  hErecoMode[conv::kSKmNuECCoth_e]     -> SetFillColor ( kBlack );
  hErecoMode[conv::kSKmNuENC1pi0]           -> SetFillColor ( kGreen+5 );
  hErecoMode[conv::kSKmNuENCcoh]           -> SetFillColor ( kGreen+7 );
  hErecoMode[conv::kSKmNuENCoth_e]           -> SetFillColor ( kBlack );
  hErecoMode[conv::kSKmSigNuECCQE]      -> SetFillColor ( kGreen+8 );
  hErecoMode[conv::kSKmSigNuECC1pi]     -> SetFillColor ( kGreen+10 );
  hErecoMode[conv::kSKmSigNuECCoth]     -> SetFillColor ( kRed );
  hErecoMode[conv::kSKmNuMuCCQE]        -> SetLineColor ( kBlack   );

  hErecoMode[conv::kSKmNuMuCC1pi]       -> SetLineWidth ( 1        );
  hErecoMode[conv::kSKmNuMuCCcoh]    -> SetLineWidth ( 1        );
  hErecoMode[conv::kSKmNuMuCCoth_e]    -> SetLineWidth ( 1        );
  hErecoMode[conv::kSKmNuMuNC1pi0]          -> SetLineWidth ( 1        );
  hErecoMode[conv::kSKmNuMuNCcoh]          -> SetLineWidth ( 1        );
  hErecoMode[conv::kSKmNuMuNCoth_e]          -> SetLineWidth ( 1        );
  hErecoMode[conv::kSKmNuMuCCQE]        -> SetLineWidth ( 1        );
  hErecoMode[conv::kSKmNuMuBarCCQE]     -> SetLineWidth ( 1        );
  hErecoMode[conv::kSKmNuMuBarCC1pi]    -> SetLineWidth ( 1        );
  hErecoMode[conv::kSKmNuMuBarCCcoh] -> SetLineWidth ( 1        );
  hErecoMode[conv::kSKmNuMuBarCCoth_e] -> SetLineWidth ( 1        );
  hErecoMode[conv::kSKmNuMuBarNC1pi0]       -> SetLineWidth ( 1        );
  hErecoMode[conv::kSKmNuMuBarNCcoh]       -> SetLineWidth ( 1        );
  hErecoMode[conv::kSKmNuMuBarNCoth_e]       -> SetLineWidth ( 1        );
  hErecoMode[conv::kSKmNuECCQE]         -> SetLineWidth ( 1        );
  hErecoMode[conv::kSKmNuECC1pi]        -> SetLineWidth ( 1        );
  hErecoMode[conv::kSKmNuECCcoh]     -> SetLineWidth ( 1        );
  hErecoMode[conv::kSKmNuECCoth_e]     -> SetLineWidth ( 1        );
  hErecoMode[conv::kSKmNuENC1pi0]           -> SetLineWidth ( 1        );
  hErecoMode[conv::kSKmNuENCcoh]           -> SetLineWidth ( 1        );
  hErecoMode[conv::kSKmNuENCoth_e]           -> SetLineWidth ( 1        );
  hErecoMode[conv::kSKmSigNuECCQE]      -> SetLineWidth ( 1        );
  hErecoMode[conv::kSKmSigNuECC1pi]     -> SetLineWidth ( 1        );
  hErecoMode[conv::kSKmSigNuECCoth]     -> SetLineWidth ( 1        );

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
  stack->Add(hErecoMode[conv::kSKmNuENC1pi0]);
  stack->Add(hErecoMode[conv::kSKmNuENCcoh]);
  stack->Add(hErecoMode[conv::kSKmNuENCoth_e]);
  stack->Add(hErecoMode[conv::kSKmNuECC1pi]);
  stack->Add(hErecoMode[conv::kSKmNuECCcoh]);
  stack->Add(hErecoMode[conv::kSKmNuECCoth_e]);
  stack->Add(hErecoMode[conv::kSKmNuECCQE]);
  stack->Add(hErecoMode[conv::kSKmNuMuBarNC1pi0]);
  stack->Add(hErecoMode[conv::kSKmNuMuBarNCcoh]);
  stack->Add(hErecoMode[conv::kSKmNuMuBarNCoth_e]);
  stack->Add(hErecoMode[conv::kSKmNuMuBarCC1pi]);
  stack->Add(hErecoMode[conv::kSKmNuMuBarCCcoh]);
  stack->Add(hErecoMode[conv::kSKmNuMuBarCCoth_e]);
  stack->Add(hErecoMode[conv::kSKmNuMuBarCCQE]);
  stack->Add(hErecoMode[conv::kSKmNuMuNC1pi0]);
  stack->Add(hErecoMode[conv::kSKmNuMuNCcoh]);
  stack->Add(hErecoMode[conv::kSKmNuMuNCoth_e]);
  stack->Add(hErecoMode[conv::kSKmNuMuCCcoh]);
  stack->Add(hErecoMode[conv::kSKmNuMuCCoth_e]);
  stack->Add(hErecoMode[conv::kSKmNuMuCC1pi]);
  stack->Add(hErecoMode[conv::kSKmNuMuCCQE]);
  stack->Draw("same");
  //lorena
  //hErecoMode[conv::kSKmSigNuECC1pi]->Draw("same");

  legend->AddEntry(hErecoTot,                              "total",                                "L" );
  legend->AddEntry(hErecoMode[conv::kSKmNuMuCCQE],         "#nu_{#mu} CCQE",                       "LF");
  legend->AddEntry(hErecoMode[conv::kSKmNuMuCC1pi],        "#nu_{#mu} CC1#pi",                     "LF");
  legend->AddEntry(hErecoMode[conv::kSKmNuMuCCcoh],        "#nu_{#mu} CC coherence",                   "LF");
  legend->AddEntry(hErecoMode[conv::kSKmNuMuCCoth_e],      "#nu_{#mu} CC other",                   "LF");
  legend->AddEntry(hErecoMode[conv::kSKmNuMuNC1pi0],       "#nu_{#mu}(#nu_{#tau}) NC1pi0",             "LF");
  legend->AddEntry(hErecoMode[conv::kSKmNuMuNCcoh],        "#nu_{#mu}(#nu_{#tau}) NC coherence",             "LF");
  legend->AddEntry(hErecoMode[conv::kSKmNuMuNCoth_e],      "#nu_{#mu}(#nu_{#tau}) NC other",             "LF");
  legend->AddEntry(hErecoMode[conv::kSKmNuMuBarCCQE],      "#bar{#nu_{#mu}} CCQE",                 "LF");
  legend->AddEntry(hErecoMode[conv::kSKmNuMuBarCC1pi],     "#bar{#nu_{#mu}} CC1#pi",               "LF");
  legend->AddEntry(hErecoMode[conv::kSKmNuMuBarCCcoh],     "#bar{#nu_{#mu}} CC coherence",             "LF");
  legend->AddEntry(hErecoMode[conv::kSKmNuMuBarCCoth_e],   "#bar{#nu_{#mu}} CC other",             "LF");
  legend->AddEntry(hErecoMode[conv::kSKmNuMuBarNC1pi0],    "#bar{#nu_{#mu}}(#bar{#nu_{#tau}}) NC1pi0", "LF");
  legend->AddEntry(hErecoMode[conv::kSKmNuMuBarNCcoh],     "#bar{#nu_{#mu}}(#bar{#nu_{#tau}}) NC coherence", "LF");
  legend->AddEntry(hErecoMode[conv::kSKmNuMuBarNCoth_e],   "#bar{#nu_{#mu}}(#bar{#nu_{#tau}}) NC other", "LF");
  legend->AddEntry(hErecoMode[conv::kSKmNuECCQE],          "#nu_{e} CCQE",                         "LF");
  legend->AddEntry(hErecoMode[conv::kSKmNuECC1pi],         "#nu_{e} CC1#pi",                       "LF");
  legend->AddEntry(hErecoMode[conv::kSKmNuECCcoh],         "#nu_{e} CC coherence",                     "LF");
  legend->AddEntry(hErecoMode[conv::kSKmNuECCoth_e],       "#nu_{e} CC other",                     "LF");
  legend->AddEntry(hErecoMode[conv::kSKmNuENC1pi0],        "#nu_{e} NC1pi0",                           "LF");
  legend->AddEntry(hErecoMode[conv::kSKmNuENCcoh],         "#nu_{e} NC coherence",                           "LF");
  legend->AddEntry(hErecoMode[conv::kSKmNuENCoth_e],       "#nu_{e} NC other",                           "LF");
  legend->AddEntry(hErecoMode[conv::kSKmSigNuECCQE],       "#nu_{e} (sig) CCQE",                   "LF");
  legend->AddEntry(hErecoMode[conv::kSKmSigNuECC1pi],      "#nu_{e} (sig) CC1#pi",                 "LF");
  legend->AddEntry(hErecoMode[conv::kSKmSigNuECCoth],      "#nu_{e} (sig) CC other",               "LF");
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
    conv::SuperkMode_t m = modes[im];
     double nmode = hErecoMode[m]->Integral();
     cout << "Contrib from mode: " 
          << setw(20) << setfill(' ') << util::AsString(m) << " : " 
          << setw(15) << setfill(' ') << setprecision(6) << nmode << " events ("
          << setw(10) << setfill(' ') << setprecision(6) << 100*nmode/ntot << " %)" << endl;
  }

  //
  // save histograms in a ROOT file
  //
  string oscparams = util::AsString(ph.Hypothesis(),
     input_mixing12,input_mixing13,input_mixing23,input_deltacp,input_dm12sq,input_dm23sq);
  TFile f(Form("sk_ereco_spectrum-%s.root",oscparams.c_str()),"recreate");
  hErecoTot                             -> Write("total");
  hErecoMode[conv::kSKmNuMuCCQE]        -> Write("numu_ccqe");
  hErecoMode[conv::kSKmNuMuCC1pi]       -> Write("numu_cc1pi"); 
  hErecoMode[conv::kSKmNuMuCCcoh]       -> Write("numu_cccoh");
  hErecoMode[conv::kSKmNuMuCCoth_e]     -> Write("numu_ccoth");
  hErecoMode[conv::kSKmNuMuNC1pi0]      -> Write("numu_and_nutau_nc_1pi0");
  hErecoMode[conv::kSKmNuMuNCcoh]       -> Write("numu_and_nutau_nc_coh");
  hErecoMode[conv::kSKmNuMuNCoth_e]     -> Write("numu_and_nutau_nc_oth");
  hErecoMode[conv::kSKmNuMuBarCCQE]     -> Write("numubar_ccqe");
  hErecoMode[conv::kSKmNuMuBarCC1pi]    -> Write("numubar_cc1pi");
  hErecoMode[conv::kSKmNuMuBarCCcoh]    -> Write("numubar_cccoh");
  hErecoMode[conv::kSKmNuMuBarCCoth_e]  -> Write("numubar_ccoth");
  hErecoMode[conv::kSKmNuMuBarNC1pi0]   -> Write("numubar_and_nutaubar_nc_1pi0");
  hErecoMode[conv::kSKmNuMuBarNCcoh]    -> Write("numubar_and_nutaubar_nc_coh");
  hErecoMode[conv::kSKmNuMuBarNCoth_e]  -> Write("numubar_and_nutaubar_nc_oth");
  hErecoMode[conv::kSKmNuECCQE]         -> Write("nue_ccqe");
  hErecoMode[conv::kSKmNuECC1pi]        -> Write("nue_cc1pi");
  hErecoMode[conv::kSKmNuECCcoh]        -> Write("nue_cccoh");
  hErecoMode[conv::kSKmNuECCoth_e]      -> Write("nue_ccoth");
  hErecoMode[conv::kSKmNuENC1pi0]       -> Write("nue_nc_1pi0");
  hErecoMode[conv::kSKmNuENCcoh]        -> Write("nue_nc_coh");
  hErecoMode[conv::kSKmNuENCoth_e]      -> Write("nue_nc_oth");
  hErecoMode[conv::kSKmSigNuECCQE]      -> Write("nuesig_ccqe");
  hErecoMode[conv::kSKmSigNuECC1pi]     -> Write("nuesig_cc1pi");
  hErecoMode[conv::kSKmSigNuECCoth]     -> Write("nuesig_ccoth");
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
