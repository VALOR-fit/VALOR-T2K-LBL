//
// ------------------------------------------------------------------------------------------------
// Plot contributions from various event modes in the predicted SuperK 1Re or 1Rmu spectrum.
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

#include <iomanip.h>

void test_sk_ereco_spectrum_nominal
(
  const char * sk_sample            = "1Rmu",                         /* options: "1Rmu", "1Re" */
  const char * hypothesis           = "3f_osc",                       /* options: "2f_numu_disap", "3f_osc" */
  double       input_mixing12       = -99,                            /* */
  double       input_mixing13       = -99,                            /* */
  double       input_mixing23       = -99,                            /* */
  double       input_deltacp        = -99,                            /* */
  double       input_dm12sq         = -99,                            /* */
  double       input_dm23sq         = -99,                            /* */
  bool         disable_banff_tuning = false,                          /* */
  const char * config_file          = "analysis_configuration.xml",   /* */
  bool         print_spectrum       = false,                          /* */
  const char * mixparamchoice       = "sinsqtheta",
  const char * hierarchy            = "normal"
)
{
  using namespace std;
  using namespace t2kana;

  // Global options
  Opt::Instance()->gMassHierarchyAssumed = util::HierarchyFromString(hierarchy);
  Opt::Instance()->gOscMixingParamChoice = util::MixingTypeFromString(mixparamchoice);
  if(Opt::Instance()->Undefined()) {
    cerr << "Global options not defined correctly. Exiting..." << endl;
    cerr << *Opt::Instance() << endl;
    exit(1);
  }

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
    ph.SetCurrValue(conv::kPhypCPV    ,   input_deltacp);
    ph.SetCurrValue(conv::kPhypDm12sq,    input_dm12sq);
    ph.SetCurrValue(conv::kPhypDm23sq,    input_dm23sq);
  }

  string oscparams = util::AsString(ph.Hypothesis(),
     input_mixing12,input_mixing13,input_mixing23,input_deltacp,input_dm12sq,input_dm23sq);

  conv::SuperkSample_t s = util::SKSampleFromString(sk_sample); // requested sample
  if(s == conv::kSKsUndefined) {
    cerr << "** Undefined sample!" << endl;
    return;
  }
  int is = util::ArrIdx::Instance()->SuperkSampleIdx(s); // index of specified sample

  // Analysis configuration
  const char * config  = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);

  // Print options
  cout << "\n[Job inputs]" << endl;
  cout << "- Reading analysis configuration from: " << config << endl;
  cout << "- Plotting spectra for sample: " << util::AsString(s) << endl;
  cout << "- Disabling BANFF tuning: " << ((disable_banff_tuning) ? "YES" : "NO") << endl;
  cout << "- Physics hypothesis: " << ph << endl;
  cout << *Opt::Instance() << endl;

  // Histogram booking
  const int nbins = util::NumErecoBins(s);
  const double * ereco_bin_edges = util::ErecoBinEdges(s);
  TH1D * hErecoTot = new TH1D("","",nbins,ereco_bin_edges);
  TH1D * hErecoMode[conv::knm]; //< contribution from each mode  
  const int nbins_coarse = util::NumVariableXsEtrueBins(s,conv::kSKmNuMuCCQE);
  const double * ereco_bin_edges_coarse = util::VariableXsEtrueBinEdges(s,conv::kSKmNuMuCCQE);
  TH1D * hErecoTot_coarse;
  TH1D * hErecoMode_coarse[conv::knm]; //< contribution from each mode

  // Create & config superk spectrum predictor
  SuperkExpectation skexp;
  skexp.Read(config);

  // Disable BANFF tuning if requested
  if(disable_banff_tuning) {
    cout << "- Performing special systematic parameter variation which cancels the effect of BANFF tuning..." << endl;
    skexp.GetBANFF()->RestoreNominalMC();
  }

  // Generate prediction
  cout << "- Generating predictions..." << endl;
  skexp.GetSuperkMC()->SetXMEC(1);
  skexp.Generate(ph);

  // Get energy spectra
  for(int im = 0; im < conv::knm; im++) {
    conv::SuperkMode_t m = conv::kSKReactionModes[is][im];
    if(m == conv::kSKmUndefined) continue;
    cout << "- Retrieving spectrum prediction for mode: " << util::AsString(m) << endl;
    TH1D * h = skexp.CurrPredictedErecoSpectrum(s, m);  
    assert(h);
    hErecoMode[im] = new TH1D(*h);
    hErecoMode_coarse[im] = (TH1D)hErecoMode[im]->Rebin(nbins_coarse,Form("%s_coarse",hErecoMode[im]->GetName()),ereco_bin_edges_coarse);
  }
  hErecoTot = new TH1D(*skexp.CurrPredictedErecoSpectrum(s));
  hErecoTot_coarse = (TH1D)hErecoTot->Rebin(nbins_coarse,Form("%s_coarse",hErecoTot->GetName()),ereco_bin_edges_coarse);

  //
  // Draw
  //

  gStyle->SetOptStat(0);
  TCanvas * c = new TCanvas("c","",20,20,500,500);
  TLegend * legend = new TLegend(0.6, 0.4, 0.9, 0.9);
  legend->SetFillColor(0);
  legend->SetBorderSize(0);
  hErecoTot->GetXaxis()->SetRangeUser(0,7-conv::kOneMeV);
  hErecoTot->Draw();
  hErecoTot->GetXaxis()->SetTitle("E_{reco} (GeV)");
  hErecoTot->GetYaxis()->SetTitle("Number of events");
  THStack * stack = new THStack("stack","spectrum contributions");
  // Get spectrum of each component mode as defined in kSKReactionModesOrderInStack array
  // so as to enhance visibility in the stack of smaller components
  for(int imstack = 0; imstack < conv::knm; imstack++) {
    conv::SuperkMode_t m = conv::kSKReactionModesOrderInStack[is][imstack];	
    if(m == conv::kSKmUndefined) continue;
    int im = util::ArrIdx::Instance()->SuperkModeIdx(s,m);
    stack->Add(hErecoMode[im]);
  }
  stack->Draw("same");
  legend->AddEntry(hErecoTot, util::AsStringRooTex(s,conv::kSKmUndefined,""," total").c_str(), "L" );
  for(int im = 0; im < conv::knm; im++) {
    conv::SuperkMode_t m = conv::kSKReactionModes[is][im];
    if(m == conv::kSKmUndefined) continue;
    legend->AddEntry(hErecoMode[im], util::AsStringRooTex(s,m).c_str(), "LF");
  }
  legend->Draw();
  c->Update();
  c->SaveAs(Form("sk_ereco_spectrum-%s_%f_%f.pdf",input_mixing23, input_dm23sq, (disable_banff_tuning ? "_nobanff" : "")));
  hErecoTot->GetXaxis()->UnZoom();

  //draw coarse binned histograms
  TCanvas * c_coarse = new TCanvas("c_coarse","",20,20,500,500);
  TLegend * legend_coarse = new TLegend(0.6, 0.4, 0.9, 0.9);
  legend_coarse->SetName("legend_coarse");
  legend_coarse->SetFillColor(0);
  legend_coarse->SetBorderSize(0);
  hErecoTot_coarse->GetXaxis()->SetRangeUser(0,7-conv::kASmallNum);
  hErecoTot_coarse->Draw();
  hErecoTot_coarse->GetXaxis()->SetTitle("E_{reco} (GeV)");
  hErecoTot_coarse->GetYaxis()->SetTitle("Number of events");
  THStack * stack_coarse = new THStack("stack_coarse","spectrum contributions");
  // Get spectrum of each component mode as defined in kSKReactionModesOrderInStack array
  // so as to enhance visibility in the stack of smaller components
  for(int imstack = 0; imstack < conv::knm; imstack++) {
    conv::SuperkMode_t m = conv::kSKReactionModesOrderInStack[is][imstack];	
    if(m == conv::kSKmUndefined) continue;
    int im = util::ArrIdx::Instance()->SuperkModeIdx(s,m);
    stack_coarse->Add(hErecoMode_coarse[im]);
  }
  stack_coarse->Draw("same");
  legend_coarse->AddEntry(hErecoTot, util::AsStringRooTex(s,conv::kSKmUndefined,""," total").c_str(), "L" );
  for(int im = 0; im < conv::knm; im++) {
    conv::SuperkMode_t m = conv::kSKReactionModes[is][im];
    if(m == conv::kSKmUndefined) continue;
    legend_coarse->AddEntry(hErecoMode[im], util::AsStringRooTex(s,m).c_str(), "LF");
  }
  legend_coarse->Draw();
  c_coarse->Update();
  c_coarse->SaveAs(Form("sk_ereco_spectrum-%s_%f_%f.pdf",input_mixing23, input_dm23sq, (disable_banff_tuning ? "_nobanff" : "")));


  //can save this to file &/or print it
  stringstream streamer;

  streamer << fixed;
  double ntot = hErecoTot->Integral();
  double ntot_coarse = hErecoTot->Integral();
  if(TMath::Abs(ntot - ntot_coarse) > conv::kASmallNum) {
    cerr << "Integrals are different between differently binned histograms "
	 << "50MeV: " << ntot
	 << "BANFF binning: " << ntot_coarse << endl;
    exit(1);
  }

  streamer << "Total number of events: " 
       << setw(10) << setfill(' ') << setprecision(6) << ntot << endl;
  const int nslices=6;
  double slice_max[nslices] = { 0.50, 0.75, 1.00, 1.25, 1.50, 2.00 };
  for(int islice=0; islice<nslices; islice++) {
    double esmax = slice_max[islice];
    streamer << "Total number of events (< " << setw(8) << setfill(' ') << setprecision(4) << esmax << " GeV): " 
         << setw(10) << setfill(' ') << setprecision(6) << hErecoTot->Integral(1,hErecoTot->FindBin(esmax)) << endl;
  }
  for(int im = 0; im < conv::knm; im++) {
     conv::SuperkMode_t m = conv::kSKReactionModes[is][im];
     if(m == conv::kSKmUndefined) continue;
     double nmode = hErecoMode[im]->Integral();
     streamer << "Contrib from mode: " 
          << setw(40) << setfill(' ') << util::AsString(m) << " : " 
          << setw(15) << setfill(' ') << setprecision(6) << nmode << " events ("
          << setw(10) << setfill(' ') << setprecision(6) << 100*nmode/ntot << " %)" << endl;
  }

  //combine modes
  double nccqe = 0, nccnqe = 0, nnuecc = 0, nnc = 0, nmec = 0;
  double nnc1pic = 0, n1pi;
  for(int im = 0; im < conv::knm; im++) {
    conv::SuperkMode_t m = conv::kSKReactionModes[is][im];
    if(m == conv::kSKmUndefined) continue;
    double nmode = hErecoMode[im]->Integral();
    if(util::IsNC(m))
      nnc += nmode;
    else if(util::IsMEC(m))
      nmec += nmode;
    else if(util::IsNuE(m) || util::IsNuEBar(m) || util::IsOscNuE(m))
      nnuecc += nmode;
    else if(m == conv::kSKmNuMuCCQE || m == conv::kSKmNuMuBarCCQE
	    || m == conv::kSKmOscNuMuCCQE)// || m == conv::kSKmOscNuMuBarCCQE)
      nccqe += nmode;
    else
      nccnqe += nmode;
    if(m == conv::kSKmNuMuNC1piC
       || m == conv::kSKmNuMuBarNC1piC
       || m == conv::kSKmNuENC1piC
       || m == conv::kSKmNuEBarNC1piC
       )
      nnc1pic += nmode;
    if(m == conv::kSKmNuMuNCoth_mu
       || m == conv::kSKmNuMuBarNCoth_mu
       || m == conv::kSKmNuENCoth_mu
       || m == conv::kSKmNuEBarNCoth_mu
       || m == conv::kSKmNuMuNC1piC
       || m == conv::kSKmNuMuBarNC1piC
       || m == conv::kSKmNuENC1piC
       || m == conv::kSKmNuEBarNC1piC
       || m == conv::kSKmNuMuCC1pi
       || m == conv::kSKmNuMuBarCC1pi
       || m == conv::kSKmNuECC1pi
       || m == conv::kSKmNuEBarCC1pi
       || m == conv::kSKmOscNuMuCC1pi
       || m == conv::kSKmOscNuECC1pi
       )
      n1pi += nmode;
  }
  streamer << endl;
  streamer << "Contrib from modes: " 
	   << setw(40) << setfill(' ') << "CCQE" << " : " 
	   << setw(15) << setfill(' ') << setprecision(6) << nccqe << " events ("
	   << setw(10) << setfill(' ') << setprecision(6) << 100*nccqe/ntot << " %)" << endl;
  streamer << "Contrib from modes: " 
	   << setw(40) << setfill(' ') << "CCnQE" << " : " 
	   << setw(15) << setfill(' ') << setprecision(6) << nccnqe << " events ("
	   << setw(10) << setfill(' ') << setprecision(6) << 100*nccnqe/ntot << " %)" << endl;
  streamer << "Contrib from modes: " 
	   << setw(40) << setfill(' ') << "NC" << " : " 
	   << setw(15) << setfill(' ') << setprecision(6) << nnc << " events ("
	   << setw(10) << setfill(' ') << setprecision(6) << 100*nnc/ntot << " %)" << endl;
  streamer << "Contrib from modes: " 
	   << setw(40) << setfill(' ') << "$\\nu_e$CC" << " : " 
	   << setw(15) << setfill(' ') << setprecision(6) << nnuecc << " events ("
	   << setw(10) << setfill(' ') << setprecision(6) << 100*nnuecc/ntot << " %)" << endl;
  streamer << "Contrib from modes: " 
	   << setw(40) << setfill(' ') << "MEC" << " : " 
	   << setw(15) << setfill(' ') << setprecision(6) << nmec << " events ("
	   << setw(10) << setfill(' ') << setprecision(6) << 100*nmec/ntot << " %)" << endl;

  streamer << endl;
  streamer << "Contrib from modes: " 
	   << setw(40) << setfill(' ') << "Non-MEC" << " : " 
	   << setw(15) << setfill(' ') << setprecision(6) << nnuecc+nnc+nccnqe+nccqe << " events ("
	   << setw(10) << setfill(' ') << setprecision(6) << 100*(nnuecc+nnc+nccnqe+nccqe)/ntot << " %)" << endl;
  streamer << "Contrib from modes: " 
	   << setw(40) << setfill(' ') << "Non-MEC or NC" << " : " 
	   << setw(15) << setfill(' ') << setprecision(6) << nnuecc+nccnqe+nccqe << " events ("
	   << setw(10) << setfill(' ') << setprecision(6) << 100*(nnuecc+nccnqe+nccqe)/ntot << " %)" << endl;
  streamer << nnuecc+nnc+nccnqe+nccqe+nmec << endl;

  streamer << endl;
  streamer << "Contrib from modes: " 
	   << setw(40) << setfill(' ') << "NC1piC" << " : " 
	   << setw(15) << setfill(' ') << setprecision(6) << nnc1pic << " events ("
	   << setw(10) << setfill(' ') << setprecision(6) << 100*nnc1pic/ntot << " %)" << endl;
  streamer << "Contrib from modes: " 
	   << setw(40) << setfill(' ') << "1pi" << " : " 
	   << setw(15) << setfill(' ') << setprecision(6) << n1pi << " events ("
	   << setw(10) << setfill(' ') << setprecision(6) << 100*n1pi/ntot << " %)" << endl;

  // get tweaks
  PhysParamList nuisance_params = skexp.Tweaks();

  //
  // Save histograms in a ROOT file
  //

  TFile f(Form("sk_ereco_spectrum-%s%s-%s.root",oscparams.c_str(),(disable_banff_tuning ? "_nobanff" : ""),sk_sample) ,"recreate");
  nuisance_params.Write("nuisance_params");
  hErecoTot -> Write(util::AsStringCompact(s,conv::kSKmUndefined,"sk_","total").c_str());
  hErecoTot -> Write("erec_sk"); //KDI
  hErecoTot_coarse->Write(Form("%s_coarse",util::AsStringCompact(s,conv::kSKmUndefined,"sk_","total").c_str()));
  for(int im = 0; im < conv::knm; im++) {
    conv::SuperkMode_t m = conv::kSKReactionModes[is][im];
    if(m == conv::kSKmUndefined) continue;
    hErecoMode[im] -> Write(util::AsStringCompact(s,m,"sk_","").c_str());
    hErecoMode_coarse[im] -> Write(Form("sk_%s_coarse",util::AsStringCompact(s,m).c_str()));
  }
  stack->Write();
  legend->Write();
  stack_coarse->Write();
  legend_coarse->Write();
  f.Close();

  //
  // Print-out numbers requested during the analysis review
  //

  if(print_spectrum) {  
    int nprec=5;
    streamer << fixed;

    streamer << endl;
    for(int ibin = 1; ibin <= hErecoTot_coarse->GetXaxis()->GetNbins(); ibin++) {
      streamer << setw(6) << setfill(' ') << ibin
	       << " : "
	       << setw(6) << setfill(' ') << setprecision(3) << hErecoTot_coarse->GetXaxis()->GetBinLowEdge(ibin)
	       << " - "
	       << setw(6) << setfill(' ') << setprecision(3) << hErecoTot_coarse->GetXaxis()->GetBinUpEdge(ibin)
	       << " : "
	       << setw(6) << setfill(' ') << setprecision(nprec) << hErecoTot_coarse->GetBinContent(ibin)
	       << endl;
    }//ibin

    streamer << endl;
    for(int ibin = 1; ibin <= hErecoTot->GetXaxis()->GetNbins(); ibin++) {
      streamer << setw(6) << setfill(' ') << ibin
	       << " : "
	       << setw(6) << setfill(' ') << setprecision(3) << hErecoTot->GetXaxis()->GetBinLowEdge(ibin)
	       << " - "
	       << setw(6) << setfill(' ') << setprecision(3) << hErecoTot->GetXaxis()->GetBinUpEdge(ibin)
	       << " : "
	       << setw(6) << setfill(' ') << setprecision(nprec) << hErecoTot->GetBinContent(ibin)
	       << endl;
    }//ibin
  }

  //actually print
  cout << streamer.str() << endl;
  
  // get size of streamer
  long size;
  streamer.seekg(0,ifstream::end);
  size=streamer.tellg();
  streamer.seekg(0);

  //write to file
  ofstream outfile (Form("sk_ereco_spectrum-%s%s-%s.txt",oscparams.c_str(),(disable_banff_tuning ? "_nobanff" : ""),sk_sample));
  outfile.write (streamer.str().c_str(),size);
  outfile.close();

}
