//
// ------------------------------------------------------------------------------------------------
// Plot predicted SuperK reconstructed neutrino energy spectrum for various sets of SuperK
// energy scale variations.
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

void test_sk_ereco_spectrum_twk_sk_escale
(
  const char * hypothesis         = "3f_osc",                       /* options: "2f_numu_disap", "3f_osc" */
  double       input_mixing12     = -99,                            /* */
  double       input_mixing13     = -99,                            /* */
  double       input_mixing23     = -99,                            /* */
  double       input_deltacp      = -99,                            /* */
  double       input_dm12sq       = -99,                            /* */
  double       input_dm23sq       = -99,                            /* */
  const char * config_file        = "analysis_configuration.xml"    /* */
)
{
  using namespace std;
  using namespace t2kana;

  conv::SuperkSample_t s = conv::kSKsFCFV1RingMuLike;

  // Analysis configuration
  const char * config  = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);
    
  // Global options
  Opt::Instance()->gOscMixingParamChoice = conv::kNuMixSinsq2Theta;
  Opt::Instance()->gMassHierarchyAssumed = conv::kNmhNormal;

  // Set physics hypothesis
  if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsq2Theta) {
    if(input_mixing13 < -1) input_mixing13 = conv::kGlobFitSinsq2Theta13;
    if(input_mixing23 < -1) input_mixing23 = conv::kGlobFitSinsq2Theta23;
    if(input_mixing12 < -1) input_mixing12 = conv::kGlobFitSinsq2Theta12;
  }
  else {
    if(input_mixing13 < -1) input_mixing13 = conv::kGlobFitSinsqTheta13;
    if(input_mixing23 < -1) input_mixing23 = conv::kGlobFitSinsqTheta23;
    if(input_mixing12 < -1) input_mixing12 = conv::kGlobFitSinsqTheta12;
  }
  if(input_dm23sq < -1)  input_dm23sq = conv::kGlobFitDm23sq;
  if(input_dm12sq < -1)  input_dm12sq = conv::kGlobFitDm12sq;
  if(input_deltacp < -90) input_deltacp = conv::kGlobFitDeltaCP;
  PhysHypothesis ph(hypothesis);
  if(ph.Hypothesis() == conv::kHypNuMuDisapp2F) {
    ph.SetCurrValue(conv::kPhypMixing23,  input_mixing23);
    ph.SetCurrValue(conv::kPhypDm23sq,    input_dm23sq);
  } else
  if(ph.Hypothesis() == conv::kHypOsc3F) {
    ph.SetCurrValue(conv::kPhypMixing12,  input_mixing12);
    ph.SetCurrValue(conv::kPhypMixing13,  input_mixing13);
    ph.SetCurrValue(conv::kPhypMixing23,  input_mixing23);
    ph.SetCurrValue(conv::kPhypDeltaCP,   input_deltacp);
    ph.SetCurrValue(conv::kPhypDm12sq,    input_dm12sq);
    ph.SetCurrValue(conv::kPhypDm23sq,    input_dm23sq);
  }

  // Print options
  cout << "\n[Job inputs]" << endl;
  cout << "- Reading analysis configuration from: " << config << endl;
  cout << "- Physics hypothesis: \n" << ph << endl;
  cout << *Opt::Instance() << endl;

  // Tweaks  
  const int ntwkdials = 3;
  double twkdial[ntwkdials] = { -1.0, 0.0, +1.0 };

  int deflt = ntwkdials/2; //untweaked case

  // Histogram drawing options
  int hwidth [ntwkdials] = { 1,       1,       2      };
  int hcolor [ntwkdials] = { kBlue,   kBlack,  kRed   };
  int hstyle [ntwkdials] = { kDashed, kSolid, kDashed };
  
  // Create & config superk spectrum predictor
  SuperkExpectation skexp;
  skexp.Read(config);
  
  // Calc spectra
  TH1D * hEreco[ntwkdials];
  for(int idial = 0; idial < ntwkdials; idial++) {
     // tweak spectrum generation model
     skexp.GetSuperkDetector()->TweakEnergyScale( twkdial[idial] );
     // generate predictions
     skexp.Generate(ph);
     // store spectrum
     hEreco[idial] = new TH1D(*skexp.CurrPredictedErecoSpectrum(s));
     hEreco[idial]->SetDirectory(0);
  }//idial
  // get fractional err
  double maxerrpos = -1.;
  double maxerrneg = +1.;
  TH1D * hErecoErr[ntwkdials];
  for(int idial = 0; idial < ntwkdials; idial++) {
    hErecoErr[idial] = new TH1D("","",conv::kNErecoBins1Rmu,conv::kErecoBinEdges1Rmu);
    hErecoErr[idial]->SetDirectory(0);
    for(int ibin=1; ibin <= conv::kNErecoBins1Rmu; ibin++) {
      double ncurr  = hEreco[idial]->GetBinContent(ibin);
      double ndeflt = hEreco[deflt]->GetBinContent(ibin);
      if(ndeflt>0) {
        double err = (ncurr-ndeflt)/ndeflt;
        maxerrpos = TMath::Max(maxerrpos, err);
        maxerrneg = TMath::Min(maxerrneg, err);          
        hErecoErr[idial]->SetBinContent(ibin, (ncurr-ndeflt)/ndeflt);
      }
    }//ibin
  }//idial

  //
  // Draw
  //
  for(int idial = 0; idial < ntwkdials; idial++) {
     hEreco   [idial]->SetLineWidth(hwidth[idial]);
     hEreco   [idial]->SetLineColor(hcolor[idial]);
     hEreco   [idial]->SetLineStyle(hstyle[idial]);
     hErecoErr[idial]->SetLineWidth(hwidth[idial]);
     hErecoErr[idial]->SetLineColor(hcolor[idial]);
     hErecoErr[idial]->SetLineStyle(hstyle[idial]);
  }
  gStyle->SetOptStat(0);
  TCanvas * c = new TCanvas("c","",20,20,500,500);
  c->Divide(1,2);
  c->cd(1);
  TLegend * legend = new TLegend(0.3, 0.7, 0.85, 0.9);
  legend->SetFillColor(0);
  legend->SetBorderSize(0);
  legend->SetHeader("SuperK energy scale");
  for(int idial = 0; idial < ntwkdials; idial++) {
    const char * opt   = (idial==0) ? "" : "same";
    const char * label = Form("Tweaking dial = %.2f", twkdial[idial]);
    if(idial==0) {
      hEreco[idial]->GetXaxis()->SetTitle("E_{reco} (GeV)");
      hEreco[idial]->GetYaxis()->SetTitle("Number of events");
    }
    hEreco[idial]->Draw(opt);
    legend->AddEntry(hEreco[idial], label, "L");
    cout << " *** " << label << " : NumEvents = " << hEreco[idial]->Integral() << endl;
  }//idial
  legend->Draw();
  c->cd(2);
  TH1F * hframe = (TH1F*) c->GetPad(2)->DrawFrame(
     0., 1.2*maxerrneg, hErecoErr[deflt]->GetXaxis()->GetXmax(), 1.2*maxerrpos);
  hframe->GetXaxis()->SetTitle("E_{reco} (GeV)");
  hframe->GetYaxis()->SetTitle("Fractional variation");
  hframe->Draw();
  for(int idial = 0; idial < ntwkdials; idial++) {
    if(idial==deflt) continue;
    hErecoErr[idial]->Draw("same");
  }
  c->Update();

  // write-out in file
  string oscparams = util::AsString(ph.Hypothesis(),
     input_mixing12,input_mixing13,input_mixing23,input_deltacp,input_dm12sq,input_dm23sq);
  TFile f(Form("sk_ereco_spectrum-twk_sk_escale-%s.root",oscparams.c_str()),"recreate");
  for(int idial = 0; idial < ntwkdials; idial++) {
     const char * dial = "";
     if( TMath::Abs(twkdial[idial]-1.0) < conv::kASmallNum ) dial="plus_one_sigma";
     if( TMath::Abs(twkdial[idial]-0.5) < conv::kASmallNum ) dial="plus_half_sigma";
     if( TMath::Abs(twkdial[idial]    ) < conv::kASmallNum ) dial="nominal";
     if( TMath::Abs(twkdial[idial]+1.0) < conv::kASmallNum ) dial="minus_one_sigma";
     if( TMath::Abs(twkdial[idial]+0.5) < conv::kASmallNum ) dial="minus_half_sigma";
     const char * title = Form("sk_ereco_spectrum_twk_sk_energy_scale_%s",dial);
     hEreco[idial]->Write(title);
  }
  f.Close();

}
