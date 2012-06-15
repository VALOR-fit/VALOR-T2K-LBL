//
// ------------------------------------------------------------------------------------------------
// Generates histograms of nominal SK spectrum with/without BANFF tuning
// with either Numu disappearance, three flavour oscillatin or no oscillation
//
// Adapted from test_nsk.C by Steve Dennis
// ------------------------------------------------------------------------------------------------

#include <iomanip>
using namespace std;
using namespace t2kana;
using namespace t2kana::conv;
  
void nsk_check
(
  bool withbanff           = true,
  const string hyp         = "3f_osc",
  const string sample_str  = "1Rmu",
  double input_mixing13    = -99,
  double input_mixing23    = -99,
  double input_mixing12    = -99,
  double input_dm23sq      = -99,
  double input_dm12sq      = -99,
  double input_deltacp     = -99,
  const char * config_file = "analysis_configuration.xml"
)
{
  conv::SuperkSample_t s = util::SKSampleFromString(sample_str);
  int is = util::ArrIdx::Instance()->SuperkSampleIdx(s);

  const char * config = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);
  cout << "Reading analysis configuration from: " << config << endl;

  // superk spectrum predictor
  SuperkExpectation skexp;
  skexp.Read(config);

  if (!withbanff) {
    skexp.GetBANFF()->RestoreNominalMC();
  }
  
  // generate predictions
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
  PhysHypothesis_t phyp = util::PhysHypFromString(hyp);
  if (phyp == kHypNuMuDisapp2F || phyp == kHypOsc3F) {
    // fit params
    PhysHypothesis ph(phyp);

    // generate predictions
    ph.SetCurrValue(conv::kPhypMixing23,  input_mixing23);
    ph.SetCurrValue(conv::kPhypDm23sq,    input_dm23sq);
    if(ph.Hypothesis() == conv::kHypOsc3F) {
      ph.SetCurrValue(conv::kPhypMixing12,  input_mixing12);
      ph.SetCurrValue(conv::kPhypMixing13,  input_mixing13);
      ph.SetCurrValue(conv::kPhypDeltaCP,   input_deltacp);
      ph.SetCurrValue(conv::kPhypDm12sq,    input_dm12sq);
    }

    skexp.Generate(ph);
  }
  else { /* Adapted from SuperkExpectation */
    // update SK MC templates
    skexp.GetSuperkMC()->UpdateTemplates();

    // update Ereco spectrum prediction
    skexp.GetSuperkMC()->UpdateErecoSpectra();

    // copy histograms
    for(int is = 0; is < kns; is++) {
      SuperkSample_t s = kSKSamples[is];
      if(s == kSKsUndefined) continue;
      const int nbins = util::NumErecoBins(s);
      for(int ibin = 1; ibin <= nbins; ibin++) {
         double ns = skexp.GetSuperkMC()->CurrErecoSpectrum(s)->GetBinContent(ibin);
         skexp.CurrPredictedErecoSpectrum(s)->SetBinContent(ibin, ns);
         for(int im = 0; im < knm; im++) {
             SuperkMode_t m = kSKReactionModes[is][im];
             if(m == kSKmUndefined) continue;
             double nsm = skexp.GetSuperkMC()->CurrErecoSpectrum(s,m)->GetBinContent(ibin);
             skexp.CurrPredictedErecoSpectrum(s,m)->SetBinContent(ibin, nsm);       
         }//m	
      }//ibin
    }//s
  }

  ostringstream outfile_name;
  
  string phyp_short;
  if (phyp == kHypOsc3F) {
    phyp_short = "3A";
  }
  else if (phyp == kHypNuMuDisapp2F) {
    phyp_short = "2A";
  }
  else {
    phyp_short = "N";
  }
  
  outfile_name << "N-" << phyp_short << "-" << (withbanff ? "W" : "N");
  TFile outfile( (outfile_name.str() + ".root").c_str() , "RECREATE");
  
  //write the ereco spectra to file
  skexp.CurrPredictedErecoSpectrum(s)->Write("Total");
  for(int im = 0; im < conv::knm; im++) {
    SuperkMode_t m = kSKReactionModes[is][im];
    if(m == kSKmUndefined) continue;
    skexp.CurrPredictedErecoSpectrum(s, m)->Write(util::AsStringCompact(m).c_str());
  }

  TH1D * tot_hist = skexp.CurrPredictedErecoSpectrum(s);
  
  TCanvas * c = new TCanvas("c","",20,20,500,500);
  int highn = tot_hist->GetNbinsX();
  double emax  = tot_hist->GetBinCenter(highn)+ tot_hist->GetBinWidth(highn)/2.;
  
  TH1F * hframe = (TH1F*)c->DrawFrame(0,0, emax, tot_hist->GetMaximum());
  hframe->GetXaxis()->SetTitle("E_{reco} (GeV)");
  hframe->GetYaxis()->SetTitle("Number of events");
  c->SetLogy();
  
  //draw the total & individual mode spectra
  tot_hist->Draw();
  TH1D * hErecoMode[conv::knm]; //< contribution from each mode 
  for(int im = 0; im < conv::knm; im++) {
    conv::SuperkMode_t m = conv::kSKReactionModes[is][im];
    if(m == conv::kSKmUndefined) continue;
    cout << "- Retrieving spectrum prediction for mode: " << util::AsString(m) << endl;
    TH1D * h = skexp.CurrPredictedErecoSpectrum(s, m);  
    assert(h);
    hErecoMode[im] = new TH1D(*h);
  }
  THStack * stack = new THStack("stack","spectrum contributions");
  for(int imstack = 0; imstack < conv::knm; imstack++) {
    conv::SuperkMode_t m = conv::kSKReactionModesOrderInStack[is][imstack];	
    if(m == conv::kSKmUndefined) continue;
    int im = util::ArrIdx::Instance()->SuperkModeIdx(s,m);
    stack->Add(hErecoMode[im]);
  }
  stack->Draw("same");
  TLegend * legend = new TLegend(0.3, 0.4, 0.9, 0.9);
  legend->AddEntry(tot_hist, util::AsStringRooTex(s,conv::kSKmUndefined,""," total").c_str(), "L" );
  for(int im = 0; im < conv::knm; im++) {
    conv::SuperkMode_t m = conv::kSKReactionModes[is][im];
    if(m == conv::kSKmUndefined) continue;
    legend->AddEntry(hErecoMode[im], util::AsStringRooTex(s,m).c_str(), "LF");
  }
  legend->Draw();
  c->Update();
  c->Write("All");
  stack->Write("Stack");
  
  ofstream returnstr( (outfile_name.str() + ".txt").c_str() );
  returnstr << "EReco_min\tEreco_max\tn_events\n";
  
  //~ returnstr.precision(5);
  returnstr.setf(ios::fixed,ios::floatfield);

  for (int ii = 1 ; ii <= tot_hist->GetNbinsX() ; ii++) {
    
    double bin_width = tot_hist->GetBinWidth(ii);
    double bin_center = tot_hist->GetBinCenter(ii);
    double bin_min = bin_center - bin_width/2.;
    double bin_max = bin_center + bin_width/2.;;
    
    double n_events = tot_hist->GetBinContent(ii);
    
    returnstr << setprecision(2) << bin_min << "         \t" << bin_max
              << "         \t"   << setprecision(5) << n_events << endl;
              
    //~ returnstr.unsetf(ios::floatfield);
  }
  
  
}
