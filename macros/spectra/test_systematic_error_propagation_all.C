//
// ------------------------------------------------------------------------------------------------
// Test systematic errors. Looks at all known error's so far.
//


#include <iomanip.h>
#include <TXMLEngine.h>

void test_systematic_error_propagation_all
(
  const char * sksample         = "1Rmu",                         /* options: "1Rmu", "1Re" */
  double       input_mixing12   = -99,                            /* */
  double       input_mixing13   = -99,                            /* */
  double       input_mixing23   = -99,                            /* */
  double       input_deltacp    = -99,                            /* */
  double       input_dm12sq     = -99,                            /* */
  double       input_dm23sq     = -99,                            /* */
  const char * config_file      = "analysis_configuration.xml"    /* */
)
{
  using namespace std;
  using namespace t2kana;
  using namespace t2kana::util;

  // Analysis configuration
  const char * config  = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);

  // Global options
  Opt::Instance()->gOscMixingParamChoice = conv::kNuMixSinsqTheta;
  Opt::Instance()->gMassHierarchyAssumed = conv::kNmhNormal;

  // Get requested sample
  conv::SuperkSample_t s = util::SKSampleFromString(sksample);
  if(s == conv::kSKsUndefined)
    exit(1);

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
  char* hypothesis("3f_osc");
  PhysHypothesis ph(hypothesis);
  ph.SetCurrValue(conv::kPhypMixing12,  input_mixing12);
  ph.SetCurrValue(conv::kPhypMixing13,  input_mixing13);
  ph.SetCurrValue(conv::kPhypMixing23,  input_mixing23);
  ph.SetCurrValue(conv::kPhypDeltaCP,   input_deltacp);
  ph.SetCurrValue(conv::kPhypDm12sq,    input_dm12sq);
  ph.SetCurrValue(conv::kPhypDm23sq,    input_dm23sq);

  //output file name
  string oscparams = util::AsString(ph.Hypothesis(),
	  input_mixing12,input_mixing13,input_mixing23,input_deltacp,input_dm12sq,input_dm23sq);
  string output_file (Form("systematic_tables-%s-%s.tex",oscparams.c_str(),sksample));

  // write document header to outfile
  ofstream outfile (output_file.c_str());
  string f;
  f += "\\documentclass[10pt,a4paper]{article}\n";
  f += ("\\usepackage{amsmath}\n");
  f += ("\\usepackage{graphicx}\n");
  f += ("\\usepackage{dcolumn}\n");
  f += ("\\usepackage[cm]{fullpage}\n");
  f += ("\\begin{document}\n");
  outfile.write (f.c_str(),f.size());
  outfile.close();

  // Print options
  cout << "\n[Job inputs]" << endl;
  cout << "- Reading analysis configuration from: " << config << endl;
  cout << "- SuperK sample: " << util::AsString(s) << endl;
  cout << "- Physics hypothesis: " << ph << endl;
  cout << *Opt::Instance() << endl;

  // Precison for printed-out numerical results
  int nprecisionnev = 4; // <-- for number of events
  int nprecisionpct = 4; // <-- for % changes

  // Create & config superk spectrum predictor
  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read(config);

  // Get the EReco bin corresponding to the edges of SK CCQE efficiency systematic
  Long64_t ErecoBin1NuMuCCQE    = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuMuCCQE    )->FindBin(0.4-conv::kASmallNum);
  Long64_t ErecoBin2NuMuCCQE    = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuMuCCQE    )->FindBin(1.1-conv::kASmallNum);
  Long64_t ErecoBin1NuMuBarCCQE = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuMuBarCCQE )->FindBin(0.4-conv::kASmallNum);
  Long64_t ErecoBin2NuMuBarCCQE = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuMuBarCCQE )->FindBin(1.1-conv::kASmallNum);
  Long64_t ErecoBin1NuECCQE     = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuECCQE )->FindBin(0.35-conv::kASmallNum);
  Long64_t ErecoBin2NuECCQE     = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuECCQE )->FindBin(0.80-conv::kASmallNum);
  
  // Get the ETrue bin edges for XSec normalisation systematics
  Long64_t EtrueBin1CCQENorm = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuCCQE, conv::kSKMCtFin)->GetYaxis()->FindBin(1.5-conv::kASmallNum);
  Long64_t EtrueBin2CCQENorm = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuCCQE, conv::kSKMCtFin)->GetYaxis()->FindBin(3.5-conv::kASmallNum);
  Long64_t EtrueBin1CCResNorm = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuCCQE, conv::kSKMCtFin)->GetYaxis()->FindBin(2.5-conv::kASmallNum);

  // Get the ETrue bin edges for NuMuCCQE
  int NEtrueBinTempl = conv::kNEtrueBins1Rmu;
  int NEtrueBinNuMu = conv::kXsNEtrueBins;
  Long64_t EtrueBinNuMu[conv::kXsNEtrueBins+1];
  for(int itb = 0; itb < conv::kXsNEtrueBins+1; itb++) {
    EtrueBinNuMu[itb] = -1;
    EtrueBinNuMu[itb] = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuCCQE, conv::kSKMCtFin)->GetYaxis()->FindBin(conv::kXsEtrueBinEdges[util::ArrIdx::Instance()->SuperkSampleIdx(s)][util::ArrIdx::Instance()->SuperkModeIdx(s,conv::kSKmNuMuCCQE)][itb]-conv::kASmallNum);
    if (EtrueBinNuMu[itb] > NEtrueBinTempl) EtrueBinNuMu[itb] = NEtrueBinTempl;
    if (EtrueBinNuMu[itb] < 1 && itb != 0) {NEtrueBinNuMu = itb-1; break;}
  }

  // Get the ETrue bin edges for NuMuBarCCQE
  int NEtrueBinNuMuBar = conv::kXsNEtrueBins;
  Long64_t EtrueBinNuMuBar[conv::kXsNEtrueBins+1];
  for(int itb = 0; itb < conv::kXsNEtrueBins+1; itb++) {
    EtrueBinNuMuBar[itb] = -1;
    EtrueBinNuMuBar[itb] = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuBarCCQE, conv::kSKMCtFin)->GetYaxis()->FindBin(conv::kXsEtrueBinEdges[util::ArrIdx::Instance()->SuperkSampleIdx(s)][util::ArrIdx::Instance()->SuperkModeIdx(s,conv::kSKmNuMuBarCCQE)][itb]-conv::kASmallNum);
    if (EtrueBinNuMuBar[itb] > NEtrueBinTempl) EtrueBinNuMuBar[itb] = NEtrueBinTempl;
    if (EtrueBinNuMuBar[itb] < 1 && itb != 0) {NEtrueBinNuMuBar = itb-1; break;}
  }

  // Get the ETrue bin edges for NuECCQE
  int NEtrueBinNuE = conv::kXsNEtrueBins;
  Long64_t EtrueBinNuE[conv::kXsNEtrueBins+1];
  for(int itb = 0; itb < conv::kXsNEtrueBins+1; itb++) {
    EtrueBinNuE[itb] = -1;
    EtrueBinNuE[itb] = skexp->GetSuperkMC()->Template(s, conv::kSKmNuECCQE, conv::kSKMCtFin)->GetYaxis()->FindBin(conv::kXsEtrueBinEdges[util::ArrIdx::Instance()->SuperkSampleIdx(s)][util::ArrIdx::Instance()->SuperkModeIdx(s,conv::kSKmNuECCQE)][itb]-conv::kASmallNum);
    if (EtrueBinNuE[itb] > NEtrueBinTempl) EtrueBinNuE[itb] = NEtrueBinTempl;
    if (EtrueBinNuE[itb] < 1 && itb != 0) {NEtrueBinNuE = itb-1; break;}
  }

  // Get the ETrue bin edges for NuEBarCCQE
  int NEtrueBinNuEBar = conv::kXsNEtrueBins;
  Long64_t EtrueBinNuEBar[conv::kXsNEtrueBins+1];
  for(int itb = 0; itb < conv::kXsNEtrueBins+1; itb++) {
    EtrueBinNuEBar[itb] = -1;
    EtrueBinNuEBar[itb] = skexp->GetSuperkMC()->Template(s, conv::kSKmNuEBarCCQE, conv::kSKMCtFin)->GetYaxis()->FindBin(conv::kXsEtrueBinEdges[util::ArrIdx::Instance()->SuperkSampleIdx(s)][util::ArrIdx::Instance()->SuperkModeIdx(s,conv::kSKmNuEBarCCQE)][itb]-conv::kASmallNum);
    if (EtrueBinNuEBar[itb] > NEtrueBinTempl) EtrueBinNuEBar[itb] = NEtrueBinTempl;
    if (EtrueBinNuEBar[itb] < 1 && itb != 0) {NEtrueBinNuEBar = itb-1; break;}
  }

  // Get the ETrue bin edges for OscNuMuCCQE
  int NEtrueBinOscNuMu = conv::kXsNEtrueBins;
  Long64_t EtrueBinOscNuMu[conv::kXsNEtrueBins+1];
  for(int itb = 0; itb < conv::kXsNEtrueBins+1; itb++) {
    EtrueBinOscNuMu[itb] = -1;
    EtrueBinOscNuMu[itb] = skexp->GetSuperkMC()->Template(s, conv::kSKmOscNuMuCCQE, conv::kSKMCtFin)->GetYaxis()->FindBin(conv::kXsEtrueBinEdges[util::ArrIdx::Instance()->SuperkSampleIdx(s)][util::ArrIdx::Instance()->SuperkModeIdx(s,conv::kSKmNuMuCCQE)][itb]-conv::kASmallNum);
    if (EtrueBinOscNuMu[itb] > NEtrueBinTempl) EtrueBinOscNuMu[itb] = NEtrueBinTempl;
    if (EtrueBinOscNuMu[itb] < 1 && itb != 0) {NEtrueBinOscNuMu = itb-1; break;}
  }

  // Get the ETrue bin edges for OscNuMuBarCCQE
  int NEtrueBinOscNuMuBar = conv::kXsNEtrueBins;
  Long64_t EtrueBinOscNuMuBar[conv::kXsNEtrueBins+1];
  for(int itb = 0; itb < conv::kXsNEtrueBins+1; itb++) {
    EtrueBinOscNuMuBar[itb] = -1;
    EtrueBinOscNuMuBar[itb] = skexp->GetSuperkMC()->Template(s, conv::kSKmOscNuMuBarCCQE, conv::kSKMCtFin)->GetYaxis()->FindBin(conv::kXsEtrueBinEdges[util::ArrIdx::Instance()->SuperkSampleIdx(s)][util::ArrIdx::Instance()->SuperkModeIdx(s,conv::kSKmNuMuCCQE)][itb]-conv::kASmallNum);
    if (EtrueBinOscNuMuBar[itb] > NEtrueBinTempl) EtrueBinOscNuMuBar[itb] = NEtrueBinTempl;
    if (EtrueBinOscNuMuBar[itb] < 1 && itb != 0) {NEtrueBinOscNuMuBar = itb-1; break;}
  }

  // Get the ETrue bin edges for OscNuECCQE
  int NEtrueBinOscNuE = conv::kXsNEtrueBins;
  Long64_t EtrueBinOscNuE[conv::kXsNEtrueBins+1];
  for(int itb = 0; itb < conv::kXsNEtrueBins+1; itb++) {
    EtrueBinOscNuE[itb] = -1;
    EtrueBinOscNuE[itb] = skexp->GetSuperkMC()->Template(s, conv::kSKmOscNuECCQE, conv::kSKMCtFin)->GetYaxis()->FindBin(conv::kXsEtrueBinEdges[util::ArrIdx::Instance()->SuperkSampleIdx(s)][util::ArrIdx::Instance()->SuperkModeIdx(s,conv::kSKmNuMuCCQE)][itb]-conv::kASmallNum);
    if (EtrueBinOscNuE[itb] > NEtrueBinTempl) EtrueBinOscNuE[itb] = NEtrueBinTempl;
    if (EtrueBinOscNuE[itb] < 1 && itb != 0) {NEtrueBinOscNuE = itb-1; break;}
  }

  // Get the ETrue bin edges for OscNuEBarCCQE
  int NEtrueBinOscNuEBar = conv::kXsNEtrueBins;
  Long64_t EtrueBinOscNuEBar[conv::kXsNEtrueBins+1];
  for(int itb = 0; itb < conv::kXsNEtrueBins+1; itb++) {
    EtrueBinOscNuEBar[itb] = -1;
    EtrueBinOscNuEBar[itb] = skexp->GetSuperkMC()->Template(s, conv::kSKmOscNuEBarCCQE, conv::kSKMCtFin)->GetYaxis()->FindBin(conv::kXsEtrueBinEdges[util::ArrIdx::Instance()->SuperkSampleIdx(s)][util::ArrIdx::Instance()->SuperkModeIdx(s,conv::kSKmNuMuCCQE)][itb]-conv::kASmallNum);
    if (EtrueBinOscNuEBar[itb] > NEtrueBinTempl) EtrueBinOscNuEBar[itb] = NEtrueBinTempl;
    if (EtrueBinOscNuEBar[itb] < 1 && itb != 0) {NEtrueBinOscNuEBar = itb-1; break;}
  }

  vector<string> systematics;

  ////////////
  //FHC SYSTEMATICS
  ////////////

  //systematics.push_back("sk_1Rmu_ccqe_ring_counting_eff"); //not supported in this macro
  systematics.push_back("sk_energy_scale");

  //systematics.push_back("skdetfsi"); //not supported in this macro
  for(int ib = 0; ib < conv::kNumOfSKDetFSIParams; ib++) {
    systematics.push_back(Form("skdetfsi_%03d",ib));
  }

  //systematics.push_back("banff"); //not supported in this macro
  for(int ib = 0; ib < conv::kNumOfFluxParamsInBANFFfit; ib++) {
    systematics.push_back(Form("banff_%02d",ib));
  }

  systematics.push_back("banff_maqe");
  systematics.push_back("banff_mares");
  systematics.push_back("banff_ccqeE0");
  systematics.push_back("banff_ccqeE1");
  systematics.push_back("banff_ccqeE2");
  systematics.push_back("banff_cc1piE0");
  systematics.push_back("banff_cc1piE1");
  systematics.push_back("banff_nc1pi0");

  systematics.push_back("xsec_pf");
  systematics.push_back("xsec_sf");
  systematics.push_back("xsec_wshape");
  systematics.push_back("xsec_shape_ccoth");
  systematics.push_back("xsec_pilessdelta");
  systematics.push_back("xsec_bindE");

  systematics.push_back("xsec_norm_cccoh");
  systematics.push_back("xsec_norm_nc1piC");
  systematics.push_back("xsec_norm_nccoh");
  systematics.push_back("xsec_norm_ncoth");
  systematics.push_back("xsec_norm_nuetonumu");
  systematics.push_back("xsec_norm_antinutonu");

  ////////////
  //RHC SYSTEMATICS
  ////////////

  //systematics.push_back("sk_1Rmu_ccqe_ring_counting_eff"); //not supported in this macro
  systematics.push_back("sk_energy_scale_rhc");

  //systematics.push_back("skdetfsi_rhc"); //not supported in this macro
  for(int ib = 0; ib < conv::kNumOfRHCSKDetFSIParams; ib++) {
    systematics.push_back(Form("skdetfsi_%03d_rhc",ib));
  }

  //systematics.push_back("banff_rhc"); //not supported in this macro
  for(int ib = 0; ib < conv::kNumOfFluxParamsInBANFFfit; ib++) {
    systematics.push_back(Form("banff_%02d_rhc",ib));
  }

  systematics.push_back("banff_maqe_rhc");
  systematics.push_back("banff_mares_rhc");
  systematics.push_back("banff_ccqeE0_rhc");
  systematics.push_back("banff_ccqeE1_rhc");
  systematics.push_back("banff_ccqeE2_rhc");
  systematics.push_back("banff_cc1piE0_rhc");
  systematics.push_back("banff_cc1piE1_rhc");
  systematics.push_back("banff_nc1pi0_rhc");
  //banff unbinned norms
  //systematics.push_back("banff_mecnorm");
  //systematics.push_back("banff_pdd");

  systematics.push_back("xsec_pf_rhc");
  systematics.push_back("xsec_sf_rhc");
  systematics.push_back("xsec_wshape_rhc");
  systematics.push_back("xsec_shape_ccoth_rhc");
  systematics.push_back("xsec_pilessdelta_rhc");
  systematics.push_back("xsec_bindE_rhc");

  systematics.push_back("xsec_norm_cccoh_rhc");
  systematics.push_back("xsec_norm_nc1piC_rhc");
  systematics.push_back("xsec_norm_nccoh_rhc");
  systematics.push_back("xsec_norm_ncoth_rhc");
  systematics.push_back("xsec_norm_nuetonumu_rhc");
  systematics.push_back("xsec_norm_antinutonu_rhc");

  systematics.push_back("rhcfhcnorm");

  ////////////
  //OSCILLATION SYSTEMATICS
  ////////////
  systematics.push_back("mixing13");
  systematics.push_back("mixing12");
  systematics.push_back("dm12sq");
  systematics.push_back("deltacp");

  //TESTING
  //systematics.clear();
  //for(int ib = 0; ib < 11; ib++) { systematics.push_back(Form("banff_%02d",ib)); }
  //for(int ib = 0; ib < conv::kNumOfRHCSKDetFSIParams; ib++) { systematics.push_back(Form("skdetfsi_%03d_rhc",ib)); }
  //for(int ib = 0; ib < conv::kNumOfSKDetFSIParams; ib++) { systematics.push_back(Form("skdetfsi_%03d",ib)); }




  bool use_tex = true;
  stringstream streamer_summary;
  streamer_summary << setprecision(2);

  for(int isys = 0; isys < systematics.size(); isys++) {
    string systematic = systematics[isys];

    bool show_sk_ccqe_spectral_info_ereco = false;
    bool show_sk_spectral_info_etrue_flux_binning = false;
    bool show_sk_spectral_info_etrue_flux_binning_mu = false;
    bool show_sk_spectral_info_etrue_norm_binning = false;

    //sk systematics
    if(systematic.find("sk_energy_scale") == 0)
      show_sk_ccqe_spectral_info_ereco = true;
    //banff binned norms
    else if(systematic.find("banff_ccqeE") == 0 || systematic.find("banff_cc1piE") == 0)
      show_sk_spectral_info_etrue_norm_binning = true;
    //banff flux
    else if(systematic.find("banff_") == 0) {
      if(systematic.find("banff_ma") == string::npos && systematic.find("banff_nc") == string::npos) {
	//we know its a banff flux parameter. now check whether it effects numu(bar) or nue(bar) fluxes
	conv::PhysParamId_t p_flux_check = util::PhysParamIdFromString(systematic);
	cout << p_flux_check << "\t" << util::AsString(p_flux_check) << endl;
	int flux_type;
	if(util::IsFHCSystematic(p_flux_check)) {
	  int ip_flux_check = util::ArrIdx::Instance()->BANFFSystIdx(p_flux_check);
	  flux_type = conv::kNuPdgCodeOfFluxParamInBANFFfit[ip_flux_check];
	}
	else {
	  int ip_flux_check = util::ArrIdx::Instance()->RHCBANFFSystIdx(p_flux_check);
	  flux_type = conv::kNuPdgCodeOfRHCFluxParamInBANFFfit[ip_flux_check];
	}
	bool effects_numu_fluxes = (flux_type == conv::kPdgNuMu) || (flux_type == conv::kPdgNuMuBar);
	if(effects_numu_fluxes)
	  show_sk_spectral_info_etrue_flux_binning_mu = true;
	else
	  show_sk_spectral_info_etrue_flux_binning = true;
      }//banff flux param
    }
    //skdet
    else if(systematic.find("skdetfsi_") == 0) {
      show_sk_ccqe_spectral_info_ereco = true;
    }

    //
    // nominal MC (pre BANFF)
    //

    //remove BANFF tuning
    skexp->GetBANFF()->RestoreNominalMC();
    skexp->GetBANFF()->RestoreRHCNominalMC();
  
    skexp->Generate(ph);

    // total 
    double nsk_n = skexp->CurrPredictedErecoSpectrum(s)->Integral();

    // individual modes & cc/nc species totals
    double nskmodes_n[conv::knm];
    double nsknumucc_n = 0;
    double nsknumunc_n = 0;
    double nsknumubarcc_n = 0;
    double nsknumubarnc_n = 0;
    double nsknuecc_n = 0;
    double nsknuenc_n = 0;
    double nsknuebarcc_n = 0;
    double nsknuebarnc_n = 0;
    double nsknuesig_n = 0;
    double nskoscnumu_n = 0;

    double nskoscnumubar_n = 0;
    double nsknuebarsig_n =0;	

    for (int im = 0; im < conv::knm; im++) {
      conv::SuperkMode_t m = conv::kSKReactionModes[s][im];
      if(m == conv::kSKmUndefined) continue;
      if(util::ArrIdx::Instance()->SuperkModeIdx(s,m) < 0) continue;
      // use s, m (not is,im) as indices, eg fMCErecoSpectrum[s][m]
      nskmodes_n[m] = skexp->CurrPredictedErecoSpectrum(s,m)->Integral();
      if(util::IsCC(m)) {
	if(util::IsNuMu(m))
	  nsknumucc_n += nskmodes_n[m];
	else if(util::IsNuMuBar(m))
	  nsknumubarcc_n += nskmodes_n[m];
	else if(util::IsNuE(m))
	  nsknuecc_n += nskmodes_n[m];
	else if(util::IsNuEBar(m))
	  nsknuebarcc_n += nskmodes_n[m];
	else if(util::IsOscNuE(m))
	  nsknuesig_n += nskmodes_n[m];
	else if(util::IsOscNuMu(m))
	  nskoscnumu_n += nskmodes_n[m];
	else if(util::IsOscNuMuBar(m))
	  nskoscnumubar_n += nskmodes_n[m];
	else if(util::IsOscNuEBar(m))
	  nsknuebarsig_n += nskmodes_n[m];
	else {
	  cerr << "CC mode " << AsString(m) << " is not numu, numubar, nue, nuebar, nuesig, oscnumu!" << endl;
	  exit(1);
	}
      }//CC?
      else if(util::IsNC(m)) {
	if(util::IsNuMu(m))
	  nsknumunc_n += nskmodes_n[m];
	else if(util::IsNuMuBar(m))
	  nsknumubarnc_n += nskmodes_n[m];
	else if(util::IsNuE(m))
	  nsknuenc_n += nskmodes_n[m];
	else if(util::IsNuEBar(m))
	  nsknuebarnc_n += nskmodes_n[m];
	else {
	  cerr << "NC mode " << m << AsString(m) << " is not numu, numubar, nue, nuebar!" << endl;
	  exit(1);
	}
      }//NC?
      else {
	cerr << "mode " << m << AsString(m) << " is not NC or CC!" << endl;
	exit(1);
      }//not CC or NC
    } //loop over modes
 
    // spectral numu ccqe (ereco)
    double nsknumuccqe_na        = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuMuCCQE    )->Integral(1,                      ErecoBin1NuMuCCQE);
    double nsknumuccqe_nb        = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuMuCCQE    )->Integral(ErecoBin1NuMuCCQE+1,    ErecoBin2NuMuCCQE);
    double nsknumuccqe_nc        = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuMuCCQE    )->Integral(ErecoBin2NuMuCCQE+1,    conv::kNErecoBins1Rmu);
    double nsknumubarccqe_na     = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuMuBarCCQE )->Integral(1,                      ErecoBin1NuMuBarCCQE);
    double nsknumubarccqe_nb     = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuMuBarCCQE )->Integral(ErecoBin1NuMuBarCCQE+1, ErecoBin2NuMuBarCCQE);
    double nsknumubarccqe_nc     = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuMuBarCCQE )->Integral(ErecoBin2NuMuBarCCQE+1, conv::kNErecoBins1Rmu);
    double nskoscnue_ereco_na = 0, nskoscnue_ereco_nb = 0, nskoscnue_ereco_nc = 0;
    double nsknue_ereco_na = 0, nsknue_ereco_nb = 0, nsknue_ereco_nc = 0;
    double nsknumu_ereco_na = 0, nsknumu_ereco_nb = 0, nsknumu_ereco_nc = 0;
    double nsknc_ereco_na = 0, nsknc_ereco_nb = 0, nsknc_ereco_nc = 0;
    for (int im = 0; im < conv::knm; im++) {
      conv::SuperkMode_t m = conv::kSKReactionModes[s][im];
      if(m == conv::kSKmUndefined) continue;
      if(util::IsCC(m)) {
	if(util::IsNuMu(m) || util::IsNuMuBar(m) || util::IsOscNuMu(m) || util::IsOscNuMuBar(m)) {
	  nsknumu_ereco_na[itb] += skexp->CurrPredictedErecoSpectrum(s, m)->Integral(1,                  ErecoBin1NuECCQE);
	  nsknumu_ereco_nb[itb] += skexp->CurrPredictedErecoSpectrum(s, m)->Integral(ErecoBin1NuECCQE+1, ErecoBin2NuECCQE);
	  nsknumu_ereco_nc[itb] += skexp->CurrPredictedErecoSpectrum(s, m)->Integral(ErecoBin2NuECCQE+1, conv::kNErecoBins1Rmu);
	}
	else if(util::IsNuE(m) || util::IsNuEBar(m)) {
	  nsknue_ereco_na[itb] += skexp->CurrPredictedErecoSpectrum(s, m)->Integral(1,                  ErecoBin1NuECCQE);
	  nsknue_ereco_nb[itb] += skexp->CurrPredictedErecoSpectrum(s, m)->Integral(ErecoBin1NuECCQE+1, ErecoBin2NuECCQE);
	  nsknue_ereco_nc[itb] += skexp->CurrPredictedErecoSpectrum(s, m)->Integral(ErecoBin2NuECCQE+1, conv::kNErecoBins1Rmu);
	}
	else if(util::IsOscNuE(m) || util::IsOscNuEBar(m)) {
	  nskoscnue_ereco_na[itb] += skexp->CurrPredictedErecoSpectrum(s, m)->Integral(1,                  ErecoBin1NuECCQE);
	  nskoscnue_ereco_nb[itb] += skexp->CurrPredictedErecoSpectrum(s, m)->Integral(ErecoBin1NuECCQE+1, ErecoBin2NuECCQE);
	  nskoscnue_ereco_nc[itb] += skexp->CurrPredictedErecoSpectrum(s, m)->Integral(ErecoBin2NuECCQE+1, conv::kNErecoBins1Rmu);
	}
      }//cc
      else {
	  nsknc_ereco_na[itb] += skexp->CurrPredictedErecoSpectrum(s, m)->Integral(1,                  ErecoBin1NuECCQE);
	  nsknc_ereco_nb[itb] += skexp->CurrPredictedErecoSpectrum(s, m)->Integral(ErecoBin1NuECCQE+1, ErecoBin2NuECCQE);
	  nsknc_ereco_nc[itb] += skexp->CurrPredictedErecoSpectrum(s, m)->Integral(ErecoBin2NuECCQE+1, conv::kNErecoBins1Rmu);
      }//nc
    }//im

    //spectral etrue (norm)
    double nsknumuccqe_etrue_na = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCQENorm);
    double nsknumuccqe_etrue_nb = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCQENorm+1, EtrueBin2CCQENorm);
    double nsknumuccqe_etrue_nc = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin2CCQENorm+1, NEtrueBinTempl);
    double nsknumuccres_etrue_na = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuCC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCResNorm);
    double nsknumuccres_etrue_nb = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuCC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCResNorm+1, NEtrueBinTempl);
    double nsknumubarccqe_etrue_na = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuBarCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCQENorm);
    double nsknumubarccqe_etrue_nb = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuBarCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCQENorm+1, EtrueBin2CCQENorm);
    double nsknumubarccqe_etrue_nc = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuBarCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin2CCQENorm+1, NEtrueBinTempl);
    double nsknumubarccres_etrue_na = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuBarCC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCResNorm);
    double nsknumubarccres_etrue_nb = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuBarCC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCResNorm+1, NEtrueBinTempl);
    double nsknueccqe_etrue_na = skexp->GetSuperkMC()->Template(s, conv::kSKmNuECCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCQENorm);
    double nsknueccqe_etrue_nb = skexp->GetSuperkMC()->Template(s, conv::kSKmNuECCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCQENorm+1, EtrueBin2CCQENorm);
    double nsknueccqe_etrue_nc = skexp->GetSuperkMC()->Template(s, conv::kSKmNuECCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin2CCQENorm+1, NEtrueBinTempl);
    double nsknueccres_etrue_na = skexp->GetSuperkMC()->Template(s, conv::kSKmNuECC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCResNorm);
    double nsknueccres_etrue_nb = skexp->GetSuperkMC()->Template(s, conv::kSKmNuECC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCResNorm+1, NEtrueBinTempl);
    double nsknuebarccqe_etrue_na = skexp->GetSuperkMC()->Template(s, conv::kSKmNuEBarCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCQENorm);
    double nsknuebarccqe_etrue_nb = skexp->GetSuperkMC()->Template(s, conv::kSKmNuEBarCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCQENorm+1, EtrueBin2CCQENorm);
    double nsknuebarccqe_etrue_nc = skexp->GetSuperkMC()->Template(s, conv::kSKmNuEBarCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin2CCQENorm+1, NEtrueBinTempl);
    double nsknuebarccres_etrue_na = skexp->GetSuperkMC()->Template(s, conv::kSKmNuEBarCC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCResNorm);
    double nsknuebarccres_etrue_nb = skexp->GetSuperkMC()->Template(s, conv::kSKmNuEBarCC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCResNorm+1, NEtrueBinTempl);
    double nsksignueccqe_etrue_na = skexp->GetSuperkMC()->Template(s, conv::kSKmOscNuECCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCQENorm);
    double nsksignueccqe_etrue_nb = skexp->GetSuperkMC()->Template(s, conv::kSKmOscNuECCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCQENorm+1, EtrueBin2CCQENorm);
    double nsksignueccqe_etrue_nc = skexp->GetSuperkMC()->Template(s, conv::kSKmOscNuECCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin2CCQENorm+1, NEtrueBinTempl);
    double nsksignueccres_etrue_na = skexp->GetSuperkMC()->Template(s, conv::kSKmOscNuECC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCResNorm);
    double nsksignueccres_etrue_nb = skexp->GetSuperkMC()->Template(s, conv::kSKmOscNuECC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCResNorm+1, NEtrueBinTempl);
    double nskoscnumuccqe_etrue_na = skexp->GetSuperkMC()->Template(s, conv::kSKmOscNuMuCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCQENorm);
    double nskoscnumuccqe_etrue_nb = skexp->GetSuperkMC()->Template(s, conv::kSKmOscNuMuCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCQENorm+1, EtrueBin2CCQENorm);
    double nskoscnumuccqe_etrue_nc = skexp->GetSuperkMC()->Template(s, conv::kSKmOscNuMuCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin2CCQENorm+1, NEtrueBinTempl);
    double nskoscnumuccres_etrue_na = skexp->GetSuperkMC()->Template(s, conv::kSKmOscNuMuCC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCResNorm);
    double nskoscnumuccres_etrue_nb = skexp->GetSuperkMC()->Template(s, conv::kSKmOscNuMuCC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCResNorm+1, NEtrueBinTempl);

    double nsksignuebarccqe_etrue_na = skexp->GetSuperkMC()->Template(s, conv::kSKmOscNuEBarCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCQENorm);
    double nsksignuebarccqe_etrue_nb = skexp->GetSuperkMC()->Template(s, conv::kSKmOscNuEBarCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCQENorm+1, EtrueBin2CCQENorm);
    double nsksignuebarccqe_etrue_nc = skexp->GetSuperkMC()->Template(s, conv::kSKmOscNuEBarCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin2CCQENorm+1, NEtrueBinTempl);
    double nsksignuebarccres_etrue_na = skexp->GetSuperkMC()->Template(s, conv::kSKmOscNuEBarCC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCResNorm);
    double nsksignuebarccres_etrue_nb = skexp->GetSuperkMC()->Template(s, conv::kSKmOscNuEBarCC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCResNorm+1, NEtrueBinTempl);
    double nskoscnumubarccqe_etrue_na = skexp->GetSuperkMC()->Template(s, conv::kSKmOscNuMuBarCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCQENorm);
    double nskoscnumubarccqe_etrue_nb = skexp->GetSuperkMC()->Template(s, conv::kSKmOscNuMuBarCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCQENorm+1, EtrueBin2CCQENorm);
    double nskoscnumubarccqe_etrue_nc = skexp->GetSuperkMC()->Template(s, conv::kSKmOscNuMuBarCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin2CCQENorm+1, NEtrueBinTempl);
    double nskoscnumubarccres_etrue_na = skexp->GetSuperkMC()->Template(s, conv::kSKmOscNuMuBarCC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCResNorm);
    double nskoscnumubarccres_etrue_nb = skexp->GetSuperkMC()->Template(s, conv::kSKmOscNuMuBarCC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCResNorm+1, NEtrueBinTempl);

    // spectral etrue (flux binning)
    double nsknumu_etrue_n[conv::kXsNEtrueBins], nsknumubar_etrue_n[conv::kXsNEtrueBins], nsknue_etrue_n[conv::kXsNEtrueBins], nsknuebar_etrue_n[conv::kXsNEtrueBins];
    double nsksignue_etrue_n[conv::kXsNEtrueBins], nsksignuebar_etrue_n[conv::kXsNEtrueBins], nskoscnumu_etrue_n[conv::kXsNEtrueBins], nskoscnumubar_etrue_n[conv::kXsNEtrueBins];
    for(int itb = 0; itb < NEtrueBinNuMu; itb++) {
      nsknumu_etrue_n[itb] = 0.;
      for (int im = 0; im < conv::knm; im++) {
	conv::SuperkMode_t m = conv::kSKReactionModes[s][im];
	if(m == conv::kSKmUndefined) continue;
	if(util::IsNuMu(m)) {
	  nsknumu_etrue_n[itb] += skexp->GetSuperkMC()->Template(s, m, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBinNuMu[itb]+1, EtrueBinNuMu[itb+1]);
	}
      }
      //cout << "numu etrue bin edges are from " << skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuCCQE, conv::kSKMCtFin)->GetYaxis()->GetBinLowEdge(EtrueBinNuMu[itb]+1) << " to " << skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuCCQE, conv::kSKMCtFin)->GetYaxis()->GetBinUpEdge(EtrueBinNuMu[itb+1]) << endl;;
    }
    for(int itb = 0; itb < NEtrueBinNuMuBar; itb++) {
      nsknumubar_etrue_n[itb] = 0.;
      for (int im = 0; im < conv::knm; im++) {
	conv::SuperkMode_t m = conv::kSKReactionModes[s][im];
	if(m == conv::kSKmUndefined) continue;
	if(util::IsNuMuBar(m)) {
	  nsknumubar_etrue_n[itb] += skexp->GetSuperkMC()->Template(s, m, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBinNuMuBar[itb]+1, EtrueBinNuMuBar[itb+1]);
	}
      }
    }
    for(int itb = 0; itb < NEtrueBinNuE; itb++) {
      nsknue_etrue_n[itb] = 0.;
      for (int im = 0; im < conv::knm; im++) {
	conv::SuperkMode_t m = conv::kSKReactionModes[s][im];
	if(m == conv::kSKmUndefined) continue;
	if(util::IsNuE(m)) {
	  nsknue_etrue_n[itb] += skexp->GetSuperkMC()->Template(s, m, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBinNuE[itb]+1, EtrueBinNuE[itb+1]);
	}
      }
    }
    for(int itb = 0; itb < NEtrueBinNuEBar; itb++) {
      nsknuebar_etrue_n[itb] = 0.;
      for (int im = 0; im < conv::knm; im++) {
	conv::SuperkMode_t m = conv::kSKReactionModes[s][im];
	if(m == conv::kSKmUndefined) continue;
	if(util::IsNuEBar(m)) {
	  nsknuebar_etrue_n[itb] += skexp->GetSuperkMC()->Template(s, m, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBinNuEBar[itb]+1, EtrueBinNuEBar[itb+1]);
	}
      }
    }
    for(int itb = 0; itb < NEtrueBinOscNuMu; itb++) {
      nskoscnumu_etrue_n[itb] = 0.;
      for (int im = 0; im < conv::knm; im++) {
	conv::SuperkMode_t m = conv::kSKReactionModes[s][im];
	if(m == conv::kSKmUndefined) continue;
	if(util::IsOscNuMu(m)) {
	  nskoscnumu_etrue_n[itb] += skexp->GetSuperkMC()->Template(s, m, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBinOscNuMu[itb]+1, EtrueBinOscNuMu[itb+1]);
	}
      }
    }
    for(int itb = 0; itb < NEtrueBinOscNuMuBar; itb++) {
      nskoscnumubar_etrue_n[itb] = 0.;
      for (int im = 0; im < conv::knm; im++) {
	conv::SuperkMode_t m = conv::kSKReactionModes[s][im];
	if(m == conv::kSKmUndefined) continue;
	if(util::IsOscNuMuBar(m)) {
	  nskoscnumubar_etrue_n[itb] += skexp->GetSuperkMC()->Template(s, m, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBinOscNuMuBar[itb]+1, EtrueBinOscNuMuBar[itb+1]);
	}
      }
    }
    for(int itb = 0; itb < NEtrueBinOscNuE; itb++) {
      nsksignue_etrue_n[itb] = 0.;
      for (int im = 0; im < conv::knm; im++) {
	conv::SuperkMode_t m = conv::kSKReactionModes[s][im];
	if(m == conv::kSKmUndefined) continue;
	if(util::IsOscNuE(m)) {
	  nsksignue_etrue_n[itb] += skexp->GetSuperkMC()->Template(s, m, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBinOscNuE[itb]+1, EtrueBinOscNuE[itb+1]);
	}
      }
    }
    for(int itb = 0; itb < NEtrueBinOscNuEBar; itb++) {
      nsksignuebar_etrue_n[itb] = 0.;
      for (int im = 0; im < conv::knm; im++) {
	conv::SuperkMode_t m = conv::kSKReactionModes[s][im];
	if(m == conv::kSKmUndefined) continue;
	if(util::IsOscNuEBar(m)) {
	  nsksignuebar_etrue_n[itb] += skexp->GetSuperkMC()->Template(s, m, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBinOscNuEBar[itb]+1, EtrueBinOscNuEBar[itb+1]);
	}
      }
    }

    //
    // nominal BANFF tuned MC
    //

    //restore BANFF tuning
    TVectorD twkv(skexp->GetBANFF()->BestFitValues());
    skexp->GetBANFF()->Tweak(true,&twkv);
    TVectorD twkv_rhc(skexp->GetBANFF()->RHCBestFitValues());
    skexp->GetBANFF()->TweakRHC(true,&twkv_rhc);

    skexp->Generate(ph);

    // total 
    double nsk_0 = skexp->CurrPredictedErecoSpectrum(s)->Integral();

    // individual modes & cc/nc species totals
    double nskmodes_0[conv::knm];
    double nsknumucc_0 = 0;
    double nsknumunc_0 = 0;
    double nsknumubarcc_0 = 0;
    double nsknumubarnc_0 = 0;
    double nsknuecc_0 = 0;
    double nsknuenc_0 = 0;
    double nsknuebarcc_0 = 0;
    double nsknuebarnc_0 = 0;
    double nsknuesig_0 = 0;
    double nskoscnumu_0 = 0;

    double nskoscnumubar_0 = 0;
    double nsknuebarsig_0 = 0;

    for (int im = 0; im < conv::knm; im++) {
      conv::SuperkMode_t m = conv::kSKReactionModes[s][im];
      if(m == conv::kSKmUndefined) continue;
      if(util::ArrIdx::Instance()->SuperkModeIdx(s,m) < 0) continue;
      // use s, m (not is,im) as indices, eg fMCErecoSpectrum[s][m]
      nskmodes_0[m] = skexp->CurrPredictedErecoSpectrum(s,m)->Integral();

      if(util::IsCC(m)) {
	if(util::IsNuMu(m))
	  nsknumucc_0 += nskmodes_0[m];
	else if(util::IsNuMuBar(m))
	  nsknumubarcc_0 += nskmodes_0[m];
	else if(util::IsNuE(m))
	  nsknuecc_0 += nskmodes_0[m];
	else if(util::IsNuEBar(m))
	  nsknuebarcc_0 += nskmodes_0[m];
	else if(util::IsOscNuE(m))
	  nsknuesig_0 += nskmodes_0[m];
	else if(util::IsOscNuMu(m))
	  nskoscnumu_0 += nskmodes_0[m];
	else if(util::IsOscNuMuBar(m))
	  nskoscnumubar_0 += nskmodes_0[m];
	else if(util::IsOscNuEBar(m))
	  nsknuebarsig_0 += nskmodes_0[m];
	else {
	  cerr << "CC mode " << AsString(m) << " is not numu, numubar, nue, nuebar, nuesig, oscnumu!2" << endl;
	  exit(1);
	}
      }//CC?
      else if(util::IsNC(m)) {
	if(util::IsNuMu(m))
	  nsknumunc_0 += nskmodes_0[m];
	else if(util::IsNuMuBar(m))
	  nsknumubarnc_0 += nskmodes_0[m];
	else if(util::IsNuE(m))
	  nsknuenc_0 += nskmodes_0[m];
	else if(util::IsNuEBar(m))
	  nsknuebarnc_0 += nskmodes_0[m];
	else {
	  cerr << "NC mode " << AsString(m) << " is not numu, numubar, nue, nuebar!" << endl;
	  exit(1);
	}
      }//NC?
      else {
	cerr << "NC mode " << AsString(m) << " is not numu, numubar, nue, nuebar!" << endl;
	exit(1);
      }//not CC or NC
    } //loop over modes

    // spectral numu ccqe (ereco)
    double nsknumuccqe_0a        = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuMuCCQE    )->Integral(1,                  ErecoBin1NuMuCCQE);
    double nsknumuccqe_0b        = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuMuCCQE    )->Integral(ErecoBin1NuMuCCQE+1,    ErecoBin2NuMuCCQE);
    double nsknumuccqe_0c        = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuMuCCQE    )->Integral(ErecoBin2NuMuCCQE+1,    conv::kNErecoBins1Rmu);
    double nsknumubarccqe_0a     = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuMuBarCCQE )->Integral(1,                  ErecoBin1NuMuBarCCQE);
    double nsknumubarccqe_0b     = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuMuBarCCQE )->Integral(ErecoBin1NuMuBarCCQE+1, ErecoBin2NuMuBarCCQE);
    double nsknumubarccqe_0c     = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuMuBarCCQE )->Integral(ErecoBin2NuMuBarCCQE+1, conv::kNErecoBins1Rmu);
    double nskoscnue_ereco_0a = 0, nskoscnue_ereco_0b = 0, nskoscnue_ereco_0c = 0;
    double nsknue_ereco_0a = 0, nsknue_ereco_0b = 0, nsknue_ereco_0c = 0;
    double nsknumu_ereco_0a = 0, nsknumu_ereco_0b = 0, nsknumu_ereco_0c = 0;
    double nsknc_ereco_0a = 0, nsknc_ereco_0b = 0, nsknc_ereco_0c = 0;
    for (int im = 0; im < conv::knm; im++) {
      conv::SuperkMode_t m = conv::kSKReactionModes[s][im];
      if(m == conv::kSKmUndefined) continue;
      if(util::IsCC(m)) {
	if(util::IsNuMu(m) || util::IsNuMuBar(m) || util::IsOscNuMu(m) || util::IsOscNuMuBar(m)) {
	  nsknumu_ereco_0a[itb] += skexp->CurrPredictedErecoSpectrum(s, m)->Integral(1,                  ErecoBin1NuECCQE);
	  nsknumu_ereco_0b[itb] += skexp->CurrPredictedErecoSpectrum(s, m)->Integral(ErecoBin1NuECCQE+1, ErecoBin2NuECCQE);
	  nsknumu_ereco_0c[itb] += skexp->CurrPredictedErecoSpectrum(s, m)->Integral(ErecoBin2NuECCQE+1, conv::kNErecoBins1Rmu);
	}
	else if(util::IsNuE(m) || util::IsNuEBar(m)) {
	  nsknue_ereco_0a[itb] += skexp->CurrPredictedErecoSpectrum(s, m)->Integral(1,                  ErecoBin1NuECCQE);
	  nsknue_ereco_0b[itb] += skexp->CurrPredictedErecoSpectrum(s, m)->Integral(ErecoBin1NuECCQE+1, ErecoBin2NuECCQE);
	  nsknue_ereco_0c[itb] += skexp->CurrPredictedErecoSpectrum(s, m)->Integral(ErecoBin2NuECCQE+1, conv::kNErecoBins1Rmu);
	}
	else if(util::IsOscNuE(m) || util::IsOscNuEBar(m)) {
	  nskoscnue_ereco_0a[itb] += skexp->CurrPredictedErecoSpectrum(s, m)->Integral(1,                  ErecoBin1NuECCQE);
	  nskoscnue_ereco_0b[itb] += skexp->CurrPredictedErecoSpectrum(s, m)->Integral(ErecoBin1NuECCQE+1, ErecoBin2NuECCQE);
	  nskoscnue_ereco_0c[itb] += skexp->CurrPredictedErecoSpectrum(s, m)->Integral(ErecoBin2NuECCQE+1, conv::kNErecoBins1Rmu);
	}
      }//cc
      else {
	  nsknc_ereco_0a[itb] += skexp->CurrPredictedErecoSpectrum(s, m)->Integral(1,                  ErecoBin1NuECCQE);
	  nsknc_ereco_0b[itb] += skexp->CurrPredictedErecoSpectrum(s, m)->Integral(ErecoBin1NuECCQE+1, ErecoBin2NuECCQE);
	  nsknc_ereco_0c[itb] += skexp->CurrPredictedErecoSpectrum(s, m)->Integral(ErecoBin2NuECCQE+1, conv::kNErecoBins1Rmu);
      }//nc
    }//im

    //spectral etrue (norm)
    double nsknumuccqe_etrue_0a = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCQENorm);
    double nsknumuccqe_etrue_0b = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCQENorm+1, EtrueBin2CCQENorm);
    double nsknumuccqe_etrue_0c = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin2CCQENorm+1, NEtrueBinTempl);
    double nsknumuccres_etrue_0a = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuCC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCResNorm);
    double nsknumuccres_etrue_0b = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuCC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCResNorm+1, NEtrueBinTempl);
    double nsknumubarccqe_etrue_0a = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuBarCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCQENorm);
    double nsknumubarccqe_etrue_0b = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuBarCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCQENorm+1, EtrueBin2CCQENorm);
    double nsknumubarccqe_etrue_0c = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuBarCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin2CCQENorm+1, NEtrueBinTempl);
    double nsknumubarccres_etrue_0a = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuBarCC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCResNorm);
    double nsknumubarccres_etrue_0b = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuBarCC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCResNorm+1, NEtrueBinTempl);
    double nsknueccqe_etrue_0a = skexp->GetSuperkMC()->Template(s, conv::kSKmNuECCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCQENorm);
    double nsknueccqe_etrue_0b = skexp->GetSuperkMC()->Template(s, conv::kSKmNuECCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCQENorm+1, EtrueBin2CCQENorm);
    double nsknueccqe_etrue_0c = skexp->GetSuperkMC()->Template(s, conv::kSKmNuECCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin2CCQENorm+1, NEtrueBinTempl);
    double nsknueccres_etrue_0a = skexp->GetSuperkMC()->Template(s, conv::kSKmNuECC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCResNorm);
    double nsknueccres_etrue_0b = skexp->GetSuperkMC()->Template(s, conv::kSKmNuECC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCResNorm+1, NEtrueBinTempl);
    double nsknuebarccqe_etrue_0a = skexp->GetSuperkMC()->Template(s, conv::kSKmNuEBarCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCQENorm);
    double nsknuebarccqe_etrue_0b = skexp->GetSuperkMC()->Template(s, conv::kSKmNuEBarCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCQENorm+1, EtrueBin2CCQENorm);
    double nsknuebarccqe_etrue_0c = skexp->GetSuperkMC()->Template(s, conv::kSKmNuEBarCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin2CCQENorm+1, NEtrueBinTempl);
    double nsknuebarccres_etrue_0a = skexp->GetSuperkMC()->Template(s, conv::kSKmNuEBarCC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCResNorm);
    double nsknuebarccres_etrue_0b = skexp->GetSuperkMC()->Template(s, conv::kSKmNuEBarCC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCResNorm+1, NEtrueBinTempl);
    double nsksignueccqe_etrue_0a = skexp->GetSuperkMC()->Template(s, conv::kSKmOscNuECCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCQENorm);
    double nsksignueccqe_etrue_0b = skexp->GetSuperkMC()->Template(s, conv::kSKmOscNuECCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCQENorm+1, EtrueBin2CCQENorm);
    double nsksignueccqe_etrue_0c = skexp->GetSuperkMC()->Template(s, conv::kSKmOscNuECCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin2CCQENorm+1, NEtrueBinTempl);
    double nsksignueccres_etrue_0a = skexp->GetSuperkMC()->Template(s, conv::kSKmOscNuECC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCResNorm);
    double nsksignueccres_etrue_0b = skexp->GetSuperkMC()->Template(s, conv::kSKmOscNuECC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCResNorm+1, NEtrueBinTempl);
    double nskoscnumuccqe_etrue_0a = skexp->GetSuperkMC()->Template(s, conv::kSKmOscNuMuCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCQENorm);
    double nskoscnumuccqe_etrue_0b = skexp->GetSuperkMC()->Template(s, conv::kSKmOscNuMuCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCQENorm+1, EtrueBin2CCQENorm);
    double nskoscnumuccqe_etrue_0c = skexp->GetSuperkMC()->Template(s, conv::kSKmOscNuMuCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin2CCQENorm+1, NEtrueBinTempl);
    double nskoscnumuccres_etrue_0a = skexp->GetSuperkMC()->Template(s, conv::kSKmOscNuMuCC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCResNorm);
    double nskoscnumuccres_etrue_0b = skexp->GetSuperkMC()->Template(s, conv::kSKmOscNuMuCC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCResNorm+1, NEtrueBinTempl);

    double nskoscnumubarccqe_etrue_0a = skexp->GetSuperkMC()->Template(s, conv::kSKmOscNuMuBarCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCQENorm);
    double nskoscnumubarccqe_etrue_0b = skexp->GetSuperkMC()->Template(s, conv::kSKmOscNuMuBarCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCQENorm+1, EtrueBin2CCQENorm);
    double nskoscnumubarccqe_etrue_0c = skexp->GetSuperkMC()->Template(s, conv::kSKmOscNuMuBarCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin2CCQENorm+1, NEtrueBinTempl);
    double nskoscnumubarccres_etrue_0a = skexp->GetSuperkMC()->Template(s, conv::kSKmOscNuMuBarCC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCResNorm);
    double nskoscnumubarccres_etrue_0b = skexp->GetSuperkMC()->Template(s, conv::kSKmOscNuMuBarCC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCResNorm+1, NEtrueBinTempl);
    
    double nsksignuebarccqe_etrue_0a = skexp->GetSuperkMC()->Template(s, conv::kSKmOscNuEBarCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCQENorm);
    double nsksignuebarccqe_etrue_0b = skexp->GetSuperkMC()->Template(s, conv::kSKmOscNuEBarCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCQENorm+1, EtrueBin2CCQENorm);
    double nsksignuebarccqe_etrue_0c = skexp->GetSuperkMC()->Template(s, conv::kSKmOscNuEBarCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin2CCQENorm+1, NEtrueBinTempl);
    double nsksignuebarccres_etrue_0a = skexp->GetSuperkMC()->Template(s, conv::kSKmOscNuEBarCC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCResNorm);
    double nsksignuebarccres_etrue_0b = skexp->GetSuperkMC()->Template(s, conv::kSKmOscNuEBarCC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCResNorm+1, NEtrueBinTempl);


  
    // spectral etrue (flux binning)
    double nsknumu_etrue_0[conv::kXsNEtrueBins], nsknumubar_etrue_0[conv::kXsNEtrueBins], nsknue_etrue_0[conv::kXsNEtrueBins], nsknuebar_etrue_0[conv::kXsNEtrueBins];
    double nsksignue_etrue_0[conv::kXsNEtrueBins], nsksignuebar_etrue_0[conv::kXsNEtrueBins], nskoscnumu_etrue_0[conv::kXsNEtrueBins], nskoscnumubar_etrue_0[conv::kXsNEtrueBins];
    for(int itb = 0; itb < NEtrueBinNuMu; itb++) {
      nsknumu_etrue_0[itb] = 0.;
      for (int im = 0; im < conv::knm; im++) {
	conv::SuperkMode_t m = conv::kSKReactionModes[s][im];
	if(m == conv::kSKmUndefined) continue;
	if(util::IsNuMu(m)) {
	  nsknumu_etrue_0[itb] += skexp->GetSuperkMC()->Template(s, m, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBinNuMu[itb]+1, EtrueBinNuMu[itb+1]);
	}
      }
    }
    for(int itb = 0; itb < NEtrueBinNuMuBar; itb++) {
      nsknumubar_etrue_0[itb] = 0.;
      for (int im = 0; im < conv::knm; im++) {
	conv::SuperkMode_t m = conv::kSKReactionModes[s][im];
	if(m == conv::kSKmUndefined) continue;
	if(util::IsNuMuBar(m)) {
	  nsknumubar_etrue_0[itb] += skexp->GetSuperkMC()->Template(s, m, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBinNuMuBar[itb]+1, EtrueBinNuMuBar[itb+1]);
	}
      }
    }
    for(int itb = 0; itb < NEtrueBinNuE; itb++) {
      nsknue_etrue_0[itb] = 0.;
      for (int im = 0; im < conv::knm; im++) {
	conv::SuperkMode_t m = conv::kSKReactionModes[s][im];
	if(m == conv::kSKmUndefined) continue;
	if(util::IsNuE(m)) {
	  nsknue_etrue_0[itb] += skexp->GetSuperkMC()->Template(s, m, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBinNuE[itb]+1, EtrueBinNuE[itb+1]);
	}
      }
    }
    for(int itb = 0; itb < NEtrueBinNuEBar; itb++) {
      nsknuebar_etrue_0[itb] = 0.;
      for (int im = 0; im < conv::knm; im++) {
	conv::SuperkMode_t m = conv::kSKReactionModes[s][im];
	if(m == conv::kSKmUndefined) continue;
	if(util::IsNuEBar(m)) {
	  nsknuebar_etrue_0[itb] += skexp->GetSuperkMC()->Template(s, m, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBinNuEBar[itb]+1, EtrueBinNuEBar[itb+1]);
	}
      }
    }
    for(int itb = 0; itb < NEtrueBinOscNuMu; itb++) {
      nskoscnumu_etrue_0[itb] = 0.;
      for (int im = 0; im < conv::knm; im++) {
	conv::SuperkMode_t m = conv::kSKReactionModes[s][im];
	if(m == conv::kSKmUndefined) continue;
	if(util::IsOscNuMu(m)) {
	  nskoscnumu_etrue_0[itb] += skexp->GetSuperkMC()->Template(s, m, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBinOscNuMu[itb]+1, EtrueBinOscNuMu[itb+1]);
	}
      }
    }
    for(int itb = 0; itb < NEtrueBinOscNuMuBar; itb++) {
      nskoscnumubar_etrue_0[itb] = 0.;
      for (int im = 0; im < conv::knm; im++) {
	conv::SuperkMode_t m = conv::kSKReactionModes[s][im];
	if(m == conv::kSKmUndefined) continue;
	if(util::IsOscNuMuBar(m)) {
	  nskoscnumubar_etrue_0[itb] += skexp->GetSuperkMC()->Template(s, m, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBinOscNuMuBar[itb]+1, EtrueBinOscNuMuBar[itb+1]);
	}
      }
    }
    for(int itb = 0; itb < NEtrueBinOscNuE; itb++) {
      nsksignue_etrue_0[itb] = 0.;
      for (int im = 0; im < conv::knm; im++) {
	conv::SuperkMode_t m = conv::kSKReactionModes[s][im];
	if(m == conv::kSKmUndefined) continue;
	if(util::IsOscNuE(m)) {
	  nsksignue_etrue_0[itb] += skexp->GetSuperkMC()->Template(s, m, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBinOscNuE[itb]+1, EtrueBinOscNuE[itb+1]);
	}
      }
    }
    for(int itb = 0; itb < NEtrueBinOscNuEBar; itb++) {
      nsksignuebar_etrue_0[itb] = 0.;
      for (int im = 0; im < conv::knm; im++) {
	conv::SuperkMode_t m = conv::kSKReactionModes[s][im];
	if(m == conv::kSKmUndefined) continue;
	if(util::IsOscNuEBar(m)) {
	  nsksignuebar_etrue_0[itb] += skexp->GetSuperkMC()->Template(s, m, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBinOscNuEBar[itb]+1, EtrueBinOscNuEBar[itb+1]);
	}
      }
    }

    //
    // +1 sigma (relative to untuned BANFF)
    //
    cout << "\n\nApplying +1 sigma tweak" << endl;
    if(systematic.find("banff_") == 0) {
      conv::PhysParamId_t p = util::PhysParamIdFromString(systematic);
      if(systematic.find("_rhc") == string::npos) {
	int ip = util::ArrIdx::Instance()->BANFFSystIdx(p);
	// curr/nominal = 1 + (nsigma * x)
	skexp->TweakSystematic(systematic, skexp->GetBANFF()->BestFitValue(p) + (1 * sqrt(skexp->GetBANFF()->CovMatrix()(ip,ip))));
	skexp->GetBANFF()->PrintSystematicVariations();
      }
      else {
	int ip = util::ArrIdx::Instance()->RHCBANFFSystIdx(p);
	// curr/nominal = 1 + (nsigma * x)
	skexp->TweakSystematic(systematic, skexp->GetBANFF()->RHCBestFitValue(p) + (1 * sqrt(skexp->GetBANFF()->RHCCovMatrix()(ip,ip))));
	skexp->GetBANFF()->PrintRHCSystematicVariations();
      }
    }
    else if(systematic.find("skdetfsi_") == 0) {
      conv::PhysParamId_t p = util::PhysParamIdFromString(systematic);
      if(systematic.find("_rhc") == string::npos) {
	int ip = util::ArrIdx::Instance()->SKDetFSISystIdx(p);
	skexp->TweakSystematic(systematic, 1 + (1 * sqrt(skexp->GetSuperkDetectorFSI()->CovMatrix()(ip,ip))));
	skexp->GetSuperkDetectorFSI()->PrintSystematicVariations();
      }
      else {
	int ip = util::ArrIdx::Instance()->RHCSKDetFSISystIdx(p);
	skexp->TweakSystematic(systematic, 1 + (1 * sqrt(skexp->GetSuperkDetectorFSI()->RHCCovMatrix()(ip,ip))));
	skexp->GetSuperkDetectorFSI()->PrintRHCSystematicVariations();
      }
    }
    else if(systematic.find("xsec_") == 0 || systematic.find("rhcfhcnorm") == 0) {
      skexp->TweakSystematic(systematic,+1.);
      if(systematic.find("_rhc") == string::npos && systematic.find("rhcfhcnorm") != 0)
	skexp->GetNuXSec()->PrintSystematicVariations();
      else
	skexp->GetNuXSec()->PrintRHCSystematicVariations();
    }
    else if(systematic.find("mixing13") == 0) {
      if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsq2Theta)
	ph.SetCurrValue(conv::kPhypMixing13,  input_mixing13 * (1 + conv::kGlobFitErrSinsq2Theta13));
      else
	ph.SetCurrValue(conv::kPhypMixing13,  input_mixing13 * (1 + conv::kGlobFitErrSinsqTheta13));
      cout << ph << endl;
    }
    else if(systematic.find("mixing12") == 0) {
      if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsq2Theta)
	ph.SetCurrValue(conv::kPhypMixing12,  input_mixing12 * (1 + conv::kGlobFitErrSinsq2Theta12));
      else
	ph.SetCurrValue(conv::kPhypMixing12,  input_mixing12 * (1 + conv::kGlobFitErrSinsqTheta12));
      cout << ph << endl;
    }
    else if(systematic.find("dm12sq") == 0) {
      ph.SetCurrValue(conv::kPhypDm12sq,    input_dm12sq * (1 + conv::kGlobFitDm12sq));
      cout << ph << endl;
    }
    else if(systematic.find("deltacp") == 0) {
      ph.SetCurrValue(conv::kPhypDeltaCP,   input_deltacp + conv::kpi/2);
      cout << ph << endl;
    }
    else
      skexp->TweakSystematic(systematic,+1.);
    skexp->Generate(ph);

    // total 
    double nsk_p = skexp->CurrPredictedErecoSpectrum(s)->Integral();

    // individual modes & cc/nc species totals
    double nskmodes_p[conv::knm];
    double nsknumucc_p = 0;
    double nsknumunc_p = 0;
    double nsknumubarcc_p = 0;
    double nsknumubarnc_p = 0;
    double nsknuecc_p = 0;
    double nsknuenc_p = 0;
    double nsknuebarcc_p = 0;
    double nsknuebarnc_p = 0;
    double nsknuesig_p = 0;
    double nskoscnumu_p = 0;

    double nsknuebarsig_p = 0;
    double nskoscnumubar_p = 0;

    for (int im = 0; im < conv::knm; im++) {
      conv::SuperkMode_t m = conv::kSKReactionModes[s][im];
      if(m == conv::kSKmUndefined) continue;
      if(util::ArrIdx::Instance()->SuperkModeIdx(s,m) < 0) continue;
      // use s, m (not is,im) as indices, eg fMCErecoSpectrum[s][m]
      nskmodes_p[m] = skexp->CurrPredictedErecoSpectrum(s,m)->Integral();

      if(util::IsCC(m)) {
	if(util::IsNuMu(m))
	  nsknumucc_p += nskmodes_p[m];
	else if(util::IsNuMuBar(m))
	  nsknumubarcc_p += nskmodes_p[m];
	else if(util::IsNuE(m))
	  nsknuecc_p += nskmodes_p[m];
	else if(util::IsNuEBar(m))
	  nsknuebarcc_p += nskmodes_p[m];
	else if(util::IsOscNuE(m))
	  nsknuesig_p += nskmodes_p[m];
	else if(util::IsOscNuMu(m))
	  nskoscnumu_p += nskmodes_p[m];
        else if(util::IsOscNuEBar(m))
	  nsknuebarsig_p += nskmodes_p[m];
	else if(util::IsOscNuMuBar(m))
	  nskoscnumubar_p += nskmodes_p[m];
	else {
	  cerr << "CC mode " << AsString(m) << " is not numu, numubar, nue, nuebar, nuesig, oscnumu!" << endl;
	  exit(1);
	}
      }//CC?
      else if(util::IsNC(m)) {
	if(util::IsNuMu(m))
	  nsknumunc_p += nskmodes_p[m];
	else if(util::IsNuMuBar(m))
	  nsknumubarnc_p += nskmodes_p[m];
	else if(util::IsNuE(m))
	  nsknuenc_p += nskmodes_p[m];
	else if(util::IsNuEBar(m))
	  nsknuebarnc_p += nskmodes_p[m];
	else {
	  cerr << "NC mode " << AsString(m) << " is not numu, numubar, nue, nuebar!" << endl;
	  exit(1);
	}
      }//NC?
      else {
	cerr << "NC mode " << AsString(m) << " is not numu, numubar, nue, nuebar!" << endl;
	exit(1);
      }//not CC or NC
    } //loop over modes

    // spectral numu ccqe (ereco)
    double nsknumuccqe_pa        = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuMuCCQE    )->Integral(1,                  ErecoBin1NuMuCCQE);
    double nsknumuccqe_pb        = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuMuCCQE    )->Integral(ErecoBin1NuMuCCQE+1,    ErecoBin2NuMuCCQE);
    double nsknumuccqe_pc        = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuMuCCQE    )->Integral(ErecoBin2NuMuCCQE+1,    conv::kNErecoBins1Rmu);
    double nsknumubarccqe_pa     = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuMuBarCCQE )->Integral(1,                  ErecoBin1NuMuBarCCQE);
    double nsknumubarccqe_pb     = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuMuBarCCQE )->Integral(ErecoBin1NuMuBarCCQE+1, ErecoBin2NuMuBarCCQE);
    double nsknumubarccqe_pc     = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuMuBarCCQE )->Integral(ErecoBin2NuMuBarCCQE+1, conv::kNErecoBins1Rmu);
    double nskoscnue_ereco_pa = 0, nskoscnue_ereco_pb = 0, nskoscnue_ereco_pc = 0;
    double nsknue_ereco_pa = 0, nsknue_ereco_pb = 0, nsknue_ereco_pc = 0;
    double nsknumu_ereco_pa = 0, nsknumu_ereco_pb = 0, nsknumu_ereco_pc = 0;
    double nsknc_ereco_pa = 0, nsknc_ereco_pb = 0, nsknc_ereco_pc = 0;
    for (int im = 0; im < conv::knm; im++) {
      conv::SuperkMode_t m = conv::kSKReactionModes[s][im];
      if(m == conv::kSKmUndefined) continue;
      if(util::IsCC(m)) {
	if(util::IsNuMu(m) || util::IsNuMuBar(m) || util::IsOscNuMu(m) || util::IsOscNuMuBar(m)) {
	  nsknumu_ereco_pa[itb] += skexp->CurrPredictedErecoSpectrum(s, m)->Integral(1,                  ErecoBin1NuECCQE);
	  nsknumu_ereco_pb[itb] += skexp->CurrPredictedErecoSpectrum(s, m)->Integral(ErecoBin1NuECCQE+1, ErecoBin2NuECCQE);
	  nsknumu_ereco_pc[itb] += skexp->CurrPredictedErecoSpectrum(s, m)->Integral(ErecoBin2NuECCQE+1, conv::kNErecoBins1Rmu);
	}
	else if(util::IsNuE(m) || util::IsNuEBar(m)) {
	  nsknue_ereco_pa[itb] += skexp->CurrPredictedErecoSpectrum(s, m)->Integral(1,                  ErecoBin1NuECCQE);
	  nsknue_ereco_pb[itb] += skexp->CurrPredictedErecoSpectrum(s, m)->Integral(ErecoBin1NuECCQE+1, ErecoBin2NuECCQE);
	  nsknue_ereco_pc[itb] += skexp->CurrPredictedErecoSpectrum(s, m)->Integral(ErecoBin2NuECCQE+1, conv::kNErecoBins1Rmu);
	}
	else if(util::IsOscNuE(m) || util::IsOscNuEBar(m)) {
	  nskoscnue_ereco_pa[itb] += skexp->CurrPredictedErecoSpectrum(s, m)->Integral(1,                  ErecoBin1NuECCQE);
	  nskoscnue_ereco_pb[itb] += skexp->CurrPredictedErecoSpectrum(s, m)->Integral(ErecoBin1NuECCQE+1, ErecoBin2NuECCQE);
	  nskoscnue_ereco_pc[itb] += skexp->CurrPredictedErecoSpectrum(s, m)->Integral(ErecoBin2NuECCQE+1, conv::kNErecoBins1Rmu);
	}
      }//cc
      else {
	  nsknc_ereco_pa[itb] += skexp->CurrPredictedErecoSpectrum(s, m)->Integral(1,                  ErecoBin1NuECCQE);
	  nsknc_ereco_pb[itb] += skexp->CurrPredictedErecoSpectrum(s, m)->Integral(ErecoBin1NuECCQE+1, ErecoBin2NuECCQE);
	  nsknc_ereco_pc[itb] += skexp->CurrPredictedErecoSpectrum(s, m)->Integral(ErecoBin2NuECCQE+1, conv::kNErecoBins1Rmu);
      }//nc
    }//im

    //spectral etrue (norm)
    double nsknumuccqe_etrue_pa = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCQENorm);
    double nsknumuccqe_etrue_pb = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCQENorm+1, EtrueBin2CCQENorm);
    double nsknumuccqe_etrue_pc = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin2CCQENorm+1, NEtrueBinTempl);
    double nsknumuccres_etrue_pa = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuCC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCResNorm);
    double nsknumuccres_etrue_pb = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuCC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCResNorm+1, NEtrueBinTempl);
    double nsknumubarccqe_etrue_pa = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuBarCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCQENorm);
    double nsknumubarccqe_etrue_pb = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuBarCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCQENorm+1, EtrueBin2CCQENorm);
    double nsknumubarccqe_etrue_pc = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuBarCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin2CCQENorm+1, NEtrueBinTempl);
    double nsknumubarccres_etrue_pa = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuBarCC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCResNorm);
    double nsknumubarccres_etrue_pb = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuBarCC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCResNorm+1, NEtrueBinTempl);
    double nsknueccqe_etrue_pa = skexp->GetSuperkMC()->Template(s, conv::kSKmNuECCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCQENorm);
    double nsknueccqe_etrue_pb = skexp->GetSuperkMC()->Template(s, conv::kSKmNuECCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCQENorm+1, EtrueBin2CCQENorm);
    double nsknueccqe_etrue_pc = skexp->GetSuperkMC()->Template(s, conv::kSKmNuECCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin2CCQENorm+1, NEtrueBinTempl);
    double nsknueccres_etrue_pa = skexp->GetSuperkMC()->Template(s, conv::kSKmNuECC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCResNorm);
    double nsknueccres_etrue_pb = skexp->GetSuperkMC()->Template(s, conv::kSKmNuECC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCResNorm+1, NEtrueBinTempl);
    double nsknuebarccqe_etrue_pa = skexp->GetSuperkMC()->Template(s, conv::kSKmNuEBarCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCQENorm);
    double nsknuebarccqe_etrue_pb = skexp->GetSuperkMC()->Template(s, conv::kSKmNuEBarCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCQENorm+1, EtrueBin2CCQENorm);
    double nsknuebarccqe_etrue_pc = skexp->GetSuperkMC()->Template(s, conv::kSKmNuEBarCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin2CCQENorm+1, NEtrueBinTempl);
    double nsknuebarccres_etrue_pa = skexp->GetSuperkMC()->Template(s, conv::kSKmNuEBarCC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCResNorm);
    double nsknuebarccres_etrue_pb = skexp->GetSuperkMC()->Template(s, conv::kSKmNuEBarCC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCResNorm+1, NEtrueBinTempl);
    double nsksignueccqe_etrue_pa = skexp->GetSuperkMC()->Template(s, conv::kSKmOscNuECCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCQENorm);
    double nsksignueccqe_etrue_pb = skexp->GetSuperkMC()->Template(s, conv::kSKmOscNuECCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCQENorm+1, EtrueBin2CCQENorm);
    double nsksignueccqe_etrue_pc = skexp->GetSuperkMC()->Template(s, conv::kSKmOscNuECCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin2CCQENorm+1, NEtrueBinTempl);
    double nsksignueccres_etrue_pa = skexp->GetSuperkMC()->Template(s, conv::kSKmOscNuECC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCResNorm);
    double nsksignueccres_etrue_pb = skexp->GetSuperkMC()->Template(s, conv::kSKmOscNuECC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCResNorm+1, NEtrueBinTempl);
    double nskoscnumuccqe_etrue_pa = skexp->GetSuperkMC()->Template(s, conv::kSKmOscNuMuCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCQENorm);
    double nskoscnumuccqe_etrue_pb = skexp->GetSuperkMC()->Template(s, conv::kSKmOscNuMuCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCQENorm+1, EtrueBin2CCQENorm);
    double nskoscnumuccqe_etrue_pc = skexp->GetSuperkMC()->Template(s, conv::kSKmOscNuMuCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin2CCQENorm+1, NEtrueBinTempl);
    double nskoscnumuccres_etrue_pa = skexp->GetSuperkMC()->Template(s, conv::kSKmOscNuMuCC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCResNorm);
    double nskoscnumuccres_etrue_pb = skexp->GetSuperkMC()->Template(s, conv::kSKmOscNuMuCC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCResNorm+1, NEtrueBinTempl);
  

    double nskoscnumubarccqe_etrue_pa = skexp->GetSuperkMC()->Template(s, conv::kSKmOscNuMuBarCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCQENorm);
    double nskoscnumubarccqe_etrue_pb = skexp->GetSuperkMC()->Template(s, conv::kSKmOscNuMuBarCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCQENorm+1, EtrueBin2CCQENorm);
    double nskoscnumubarccqe_etrue_pc = skexp->GetSuperkMC()->Template(s, conv::kSKmOscNuMuBarCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin2CCQENorm+1, NEtrueBinTempl);
    double nskoscnumubarccres_etrue_pa = skexp->GetSuperkMC()->Template(s, conv::kSKmOscNuMuBarCC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCResNorm);
    double nskoscnumubarccres_etrue_pb = skexp->GetSuperkMC()->Template(s, conv::kSKmOscNuMuBarCC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCResNorm+1, NEtrueBinTempl);
    double nsksignuebarccqe_etrue_pa = skexp->GetSuperkMC()->Template(s, conv::kSKmOscNuEBarCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCQENorm);
    double nsksignuebarccqe_etrue_pb = skexp->GetSuperkMC()->Template(s, conv::kSKmOscNuEBarCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCQENorm+1, EtrueBin2CCQENorm);
    double nsksignuebarccqe_etrue_pc = skexp->GetSuperkMC()->Template(s, conv::kSKmOscNuEBarCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin2CCQENorm+1, NEtrueBinTempl);
    double nsksignuebarccres_etrue_pa = skexp->GetSuperkMC()->Template(s, conv::kSKmOscNuEBarCC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCResNorm);
    double nsksignuebarccres_etrue_pb = skexp->GetSuperkMC()->Template(s, conv::kSKmOscNuEBarCC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCResNorm+1, NEtrueBinTempl);

    // spectral etrue (flux binning)
    double nsknumu_etrue_p[conv::kXsNEtrueBins], nsknumubar_etrue_p[conv::kXsNEtrueBins], nsknue_etrue_p[conv::kXsNEtrueBins], nsknuebar_etrue_p[conv::kXsNEtrueBins];
    double nsksignue_etrue_p[conv::kXsNEtrueBins], nsksignuebar_etrue_p[conv::kXsNEtrueBins], nskoscnumu_etrue_p[conv::kXsNEtrueBins], nskoscnumubar_etrue_p[conv::kXsNEtrueBins];
    for(int itb = 0; itb < NEtrueBinNuMu; itb++) {
      nsknumu_etrue_p[itb] = 0.;
      for (int im = 0; im < conv::knm; im++) {
	conv::SuperkMode_t m = conv::kSKReactionModes[s][im];
	if(m == conv::kSKmUndefined) continue;
	if(util::IsNuMu(m)) {
	  nsknumu_etrue_p[itb] += skexp->GetSuperkMC()->Template(s, m, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBinNuMu[itb]+1, EtrueBinNuMu[itb+1]);
	}
      }
    }
    for(int itb = 0; itb < NEtrueBinNuMuBar; itb++) {
      nsknumubar_etrue_p[itb] = 0.;
      for (int im = 0; im < conv::knm; im++) {
	conv::SuperkMode_t m = conv::kSKReactionModes[s][im];
	if(m == conv::kSKmUndefined) continue;
	if(util::IsNuMuBar(m)) {
	  nsknumubar_etrue_p[itb] += skexp->GetSuperkMC()->Template(s, m, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBinNuMuBar[itb]+1, EtrueBinNuMuBar[itb+1]);
	}
      }
    }
    for(int itb = 0; itb < NEtrueBinNuE; itb++) {
      nsknue_etrue_p[itb] = 0.;
      for (int im = 0; im < conv::knm; im++) {
	conv::SuperkMode_t m = conv::kSKReactionModes[s][im];
	if(m == conv::kSKmUndefined) continue;
	if(util::IsNuE(m)) {
	  nsknue_etrue_p[itb] += skexp->GetSuperkMC()->Template(s, m, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBinNuE[itb]+1, EtrueBinNuE[itb+1]);
	}
      }
    }
    for(int itb = 0; itb < NEtrueBinNuEBar; itb++) {
      nsknuebar_etrue_p[itb] = 0.;
      for (int im = 0; im < conv::knm; im++) {
	conv::SuperkMode_t m = conv::kSKReactionModes[s][im];
	if(m == conv::kSKmUndefined) continue;
	if(util::IsNuEBar(m)) {
	  nsknuebar_etrue_p[itb] += skexp->GetSuperkMC()->Template(s, m, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBinNuEBar[itb]+1, EtrueBinNuEBar[itb+1]);
	}
      }
    }
    for(int itb = 0; itb < NEtrueBinOscNuMu; itb++) {
      nskoscnumu_etrue_p[itb] = 0.;
      for (int im = 0; im < conv::knm; im++) {
	conv::SuperkMode_t m = conv::kSKReactionModes[s][im];
	if(m == conv::kSKmUndefined) continue;
	if(util::IsOscNuMu(m)) {
	  nskoscnumu_etrue_p[itb] += skexp->GetSuperkMC()->Template(s, m, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBinOscNuMu[itb]+1, EtrueBinOscNuMu[itb+1]);
	}
      }
    }
    for(int itb = 0; itb < NEtrueBinOscNuMuBar; itb++) {
      nskoscnumubar_etrue_p[itb] = 0.;
      for (int im = 0; im < conv::knm; im++) {
	conv::SuperkMode_t m = conv::kSKReactionModes[s][im];
	if(m == conv::kSKmUndefined) continue;
	if(util::IsOscNuMuBar(m)) {
	  nskoscnumubar_etrue_p[itb] += skexp->GetSuperkMC()->Template(s, m, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBinOscNuMuBar[itb]+1, EtrueBinOscNuMuBar[itb+1]);
	}
      }
    }
    for(int itb = 0; itb < NEtrueBinOscNuE; itb++) {
      nsksignue_etrue_p[itb] = 0.;
      for (int im = 0; im < conv::knm; im++) {
	conv::SuperkMode_t m = conv::kSKReactionModes[s][im];
	if(m == conv::kSKmUndefined) continue;
	if(util::IsOscNuE(m)) {
	  nsksignue_etrue_p[itb] += skexp->GetSuperkMC()->Template(s, m, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBinOscNuE[itb]+1, EtrueBinOscNuE[itb+1]);
	}
      }
    }
    for(int itb = 0; itb < NEtrueBinOscNuEBar; itb++) {
      nsksignuebar_etrue_p[itb] = 0.;
      for (int im = 0; im < conv::knm; im++) {
	conv::SuperkMode_t m = conv::kSKReactionModes[s][im];
	if(m == conv::kSKmUndefined) continue;
	if(util::IsOscNuEBar(m)) {
	  nsksignuebar_etrue_p[itb] += skexp->GetSuperkMC()->Template(s, m, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBinOscNuEBar[itb]+1, EtrueBinOscNuEBar[itb+1]);
	}
      }
    }

    //
    // -1 sigma (relative to untuned BANFF)
    //
    cout << "\n\nApplying -1 sigma tweak" << endl;
    if(systematic.find("banff_") == 0) {
      conv::PhysParamId_t p = util::PhysParamIdFromString(systematic);
      if(systematic.find("_rhc") == string::npos) {
	int ip = util::ArrIdx::Instance()->BANFFSystIdx(p);
	// curr/nominal = 1 + (nsigma * x)
	skexp->TweakSystematic(systematic, skexp->GetBANFF()->BestFitValue(p) + (-1 * sqrt(skexp->GetBANFF()->CovMatrix()(ip,ip))));
	skexp->GetBANFF()->PrintSystematicVariations();
      }
      else {
	int ip = util::ArrIdx::Instance()->RHCBANFFSystIdx(p);
	// curr/nominal = 1 + (nsigma * x)
	skexp->TweakSystematic(systematic, skexp->GetBANFF()->RHCBestFitValue(p) + (-1 * sqrt(skexp->GetBANFF()->RHCCovMatrix()(ip,ip))));
	skexp->GetBANFF()->PrintRHCSystematicVariations();
      }
    }
    else if(systematic.find("skdetfsi_") == 0) {
      conv::PhysParamId_t p = util::PhysParamIdFromString(systematic);
      if(systematic.find("_rhc") == string::npos) {
	int ip = util::ArrIdx::Instance()->SKDetFSISystIdx(p);
	skexp->TweakSystematic(systematic, 1 + (-1 * sqrt(skexp->GetSuperkDetectorFSI()->CovMatrix()(ip,ip))));
	skexp->GetSuperkDetectorFSI()->PrintSystematicVariations();
      }
      else {
	int ip = util::ArrIdx::Instance()->RHCSKDetFSISystIdx(p);
	skexp->TweakSystematic(systematic, 1 + (-1 * sqrt(skexp->GetSuperkDetectorFSI()->RHCCovMatrix()(ip,ip))));
	skexp->GetSuperkDetectorFSI()->PrintRHCSystematicVariations();
      }
    }
    else if(systematic.find("xsec_") == 0 || systematic.find("rhcfhcnorm") == 0) {
      skexp->TweakSystematic(systematic,-1.);
      if(systematic.find("_rhc") == string::npos && systematic.find("rhcfhcnorm") != 0)
	skexp->GetNuXSec()->PrintSystematicVariations();
      else
	skexp->GetNuXSec()->PrintRHCSystematicVariations();
    }
    else if(systematic.find("mixing13") == 0) {
      if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsq2Theta)
	ph.SetCurrValue(conv::kPhypMixing13,  input_mixing13 * (1 - conv::kGlobFitErrSinsq2Theta13));
      else
	ph.SetCurrValue(conv::kPhypMixing13,  input_mixing13 * (1 - conv::kGlobFitErrSinsqTheta13));
      cout << ph << endl;
    }
    else if(systematic.find("mixing12") == 0) {
      if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsq2Theta)
	ph.SetCurrValue(conv::kPhypMixing12,  input_mixing12 * (1 - conv::kGlobFitErrSinsq2Theta12));
      else
	ph.SetCurrValue(conv::kPhypMixing12,  input_mixing12 * (1 - conv::kGlobFitErrSinsqTheta12));
      cout << ph << endl;
    }
    else if(systematic.find("dm12sq") == 0) {
      ph.SetCurrValue(conv::kPhypDm12sq,    input_dm12sq * (1 - conv::kGlobFitDm12sq));
      cout << ph << endl;
    }
    else if(systematic.find("deltacp") == 0) {
      ph.SetCurrValue(conv::kPhypDeltaCP,   input_deltacp - conv::kpi/2);
      cout << ph << endl;
    }
    else
      skexp->TweakSystematic(systematic,-1.);
    skexp->Generate(ph);

    // total 
    double nsk_m = skexp->CurrPredictedErecoSpectrum(s)->Integral();

    // individual modes & cc/nc species totals
    double nskmodes_m[conv::knm];
    double nsknumucc_m = 0;
    double nsknumunc_m = 0;
    double nsknumubarcc_m = 0;
    double nsknumubarnc_m = 0;
    double nsknuecc_m = 0;
    double nsknuenc_m = 0;
    double nsknuebarcc_m = 0;
    double nsknuebarnc_m = 0;
    double nsknuesig_m = 0;
    double nskoscnumu_m = 0;

    double nsknuebarsig_m = 0;
    double nskoscnumubar_m = 0;

    for (int im = 0; im < conv::knm; im++) {
      conv::SuperkMode_t m = conv::kSKReactionModes[s][im];
      if(m == conv::kSKmUndefined) continue;
      if(util::ArrIdx::Instance()->SuperkModeIdx(s,m) < 0) continue;
      // use s, m (not is,im) as indices, eg fMCErecoSpectrum[s][m]
      nskmodes_m[m] = skexp->CurrPredictedErecoSpectrum(s,m)->Integral();

      if(util::IsCC(m)) {
	if(util::IsNuMu(m))
	  nsknumucc_m += nskmodes_m[m];
	else if(util::IsNuMuBar(m))
	  nsknumubarcc_m += nskmodes_m[m];
	else if(util::IsNuE(m))
	  nsknuecc_m += nskmodes_m[m];
	else if(util::IsNuEBar(m))
	  nsknuebarcc_m += nskmodes_m[m];
	else if(util::IsOscNuE(m))
	  nsknuesig_m += nskmodes_m[m];
	else if(util::IsOscNuMu(m))
	  nskoscnumu_m += nskmodes_m[m];
	else if(util::IsOscNuEBar(m))
	  nsknuebarsig_m += nskmodes_m[m];
	else if(util::IsOscNuMuBar(m))
	  nskoscnumubar_m += nskmodes_m[m];
	else {
	  cerr << "CC mode " << AsString(m) << " is not numu, numubar, nue, nuebar, nuesig, oscnumu!3" << endl;
	  exit(1);
	}
      }//CC?
      else if(util::IsNC(m)) {
	if(util::IsNuMu(m))
	  nsknumunc_m += nskmodes_m[m];
	else if(util::IsNuMuBar(m))
	  nsknumubarnc_m += nskmodes_m[m];
	else if(util::IsNuE(m))
	  nsknuenc_m += nskmodes_m[m];
	else if(util::IsNuEBar(m))
	  nsknuebarnc_m += nskmodes_m[m];
	else {
	  cerr << "NC mode " << AsString(m) << " is not numu, numubar, nue, nuebar!" << endl;
	  exit(1);
	}
      }//NC?
      else {
	cerr << "NC mode " << AsString(m) << " is not numu, numubar, nue, nuebar!" << endl;
	exit(1);
      }//not CC or NC
    } //loop over modes

    // spectral numu ccqe (ereco)
    double nsknumuccqe_ma        = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuMuCCQE    )->Integral(1,                  ErecoBin1NuMuCCQE);
    double nsknumuccqe_mb        = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuMuCCQE    )->Integral(ErecoBin1NuMuCCQE+1,    ErecoBin2NuMuCCQE);
    double nsknumuccqe_mc        = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuMuCCQE    )->Integral(ErecoBin2NuMuCCQE+1,    conv::kNErecoBins1Rmu);
    double nsknumubarccqe_ma     = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuMuBarCCQE )->Integral(1,                  ErecoBin1NuMuBarCCQE);
    double nsknumubarccqe_mb     = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuMuBarCCQE )->Integral(ErecoBin1NuMuBarCCQE+1, ErecoBin2NuMuBarCCQE);
    double nsknumubarccqe_mc     = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuMuBarCCQE )->Integral(ErecoBin2NuMuBarCCQE+1, conv::kNErecoBins1Rmu);
    double nskoscnue_ereco_ma = 0, nskoscnue_ereco_mb = 0, nskoscnue_ereco_mc = 0;
    double nsknue_ereco_ma = 0, nsknue_ereco_mb = 0, nsknue_ereco_mc = 0;
    double nsknumu_ereco_ma = 0, nsknumu_ereco_mb = 0, nsknumu_ereco_mc = 0;
    double nsknc_ereco_ma = 0, nsknc_ereco_mb = 0, nsknc_ereco_mc = 0;
    for (int im = 0; im < conv::knm; im++) {
      conv::SuperkMode_t m = conv::kSKReactionModes[s][im];
      if(m == conv::kSKmUndefined) continue;
      if(util::IsCC(m)) {
	if(util::IsNuMu(m) || util::IsNuMuBar(m) || util::IsOscNuMu(m) || util::IsOscNuMuBar(m)) {
	  nsknumu_ereco_ma[itb] += skexp->CurrPredictedErecoSpectrum(s, m)->Integral(1,                  ErecoBin1NuECCQE);
	  nsknumu_ereco_mb[itb] += skexp->CurrPredictedErecoSpectrum(s, m)->Integral(ErecoBin1NuECCQE+1, ErecoBin2NuECCQE);
	  nsknumu_ereco_mc[itb] += skexp->CurrPredictedErecoSpectrum(s, m)->Integral(ErecoBin2NuECCQE+1, conv::kNErecoBins1Rmu);
	}
	else if(util::IsNuE(m) || util::IsNuEBar(m)) {
	  nsknue_ereco_ma[itb] += skexp->CurrPredictedErecoSpectrum(s, m)->Integral(1,                  ErecoBin1NuECCQE);
	  nsknue_ereco_mb[itb] += skexp->CurrPredictedErecoSpectrum(s, m)->Integral(ErecoBin1NuECCQE+1, ErecoBin2NuECCQE);
	  nsknue_ereco_mc[itb] += skexp->CurrPredictedErecoSpectrum(s, m)->Integral(ErecoBin2NuECCQE+1, conv::kNErecoBins1Rmu);
	}
	else if(util::IsOscNuE(m) || util::IsOscNuEBar(m)) {
	  nskoscnue_ereco_ma[itb] += skexp->CurrPredictedErecoSpectrum(s, m)->Integral(1,                  ErecoBin1NuECCQE);
	  nskoscnue_ereco_mb[itb] += skexp->CurrPredictedErecoSpectrum(s, m)->Integral(ErecoBin1NuECCQE+1, ErecoBin2NuECCQE);
	  nskoscnue_ereco_mc[itb] += skexp->CurrPredictedErecoSpectrum(s, m)->Integral(ErecoBin2NuECCQE+1, conv::kNErecoBins1Rmu);
	}
      }//cc
      else {
	  nsknc_ereco_ma[itb] += skexp->CurrPredictedErecoSpectrum(s, m)->Integral(1,                  ErecoBin1NuECCQE);
	  nsknc_ereco_mb[itb] += skexp->CurrPredictedErecoSpectrum(s, m)->Integral(ErecoBin1NuECCQE+1, ErecoBin2NuECCQE);
	  nsknc_ereco_mc[itb] += skexp->CurrPredictedErecoSpectrum(s, m)->Integral(ErecoBin2NuECCQE+1, conv::kNErecoBins1Rmu);
      }//nc
    }//im

    //spectral etrue (norm)
    double nsknumuccqe_etrue_ma = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCQENorm);
    double nsknumuccqe_etrue_mb = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCQENorm+1, EtrueBin2CCQENorm);
    double nsknumuccqe_etrue_mc = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin2CCQENorm+1, NEtrueBinTempl);
    double nsknumuccres_etrue_ma = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuCC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCResNorm);
    double nsknumuccres_etrue_mb = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuCC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCResNorm+1, NEtrueBinTempl);
    double nsknumubarccqe_etrue_ma = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuBarCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCQENorm);
    double nsknumubarccqe_etrue_mb = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuBarCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCQENorm+1, EtrueBin2CCQENorm);
    double nsknumubarccqe_etrue_mc = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuBarCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin2CCQENorm+1, NEtrueBinTempl);
    double nsknumubarccres_etrue_ma = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuBarCC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCResNorm);
    double nsknumubarccres_etrue_mb = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuBarCC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCResNorm+1, NEtrueBinTempl);
    double nsknueccqe_etrue_ma = skexp->GetSuperkMC()->Template(s, conv::kSKmNuECCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCQENorm);
    double nsknueccqe_etrue_mb = skexp->GetSuperkMC()->Template(s, conv::kSKmNuECCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCQENorm+1, EtrueBin2CCQENorm);
    double nsknueccqe_etrue_mc = skexp->GetSuperkMC()->Template(s, conv::kSKmNuECCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin2CCQENorm+1, NEtrueBinTempl);
    double nsknueccres_etrue_ma = skexp->GetSuperkMC()->Template(s, conv::kSKmNuECC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCResNorm);
    double nsknueccres_etrue_mb = skexp->GetSuperkMC()->Template(s, conv::kSKmNuECC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCResNorm+1, NEtrueBinTempl);
    double nsknuebarccqe_etrue_ma = skexp->GetSuperkMC()->Template(s, conv::kSKmNuEBarCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCQENorm);
    double nsknuebarccqe_etrue_mb = skexp->GetSuperkMC()->Template(s, conv::kSKmNuEBarCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCQENorm+1, EtrueBin2CCQENorm);
    double nsknuebarccqe_etrue_mc = skexp->GetSuperkMC()->Template(s, conv::kSKmNuEBarCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin2CCQENorm+1, NEtrueBinTempl);
    double nsknuebarccres_etrue_ma = skexp->GetSuperkMC()->Template(s, conv::kSKmNuEBarCC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCResNorm);
    double nsknuebarccres_etrue_mb = skexp->GetSuperkMC()->Template(s, conv::kSKmNuEBarCC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCResNorm+1, NEtrueBinTempl);
    double nsksignueccqe_etrue_ma = skexp->GetSuperkMC()->Template(s, conv::kSKmOscNuECCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCQENorm);
    double nsksignueccqe_etrue_mb = skexp->GetSuperkMC()->Template(s, conv::kSKmOscNuECCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCQENorm+1, EtrueBin2CCQENorm);
    double nsksignueccqe_etrue_mc = skexp->GetSuperkMC()->Template(s, conv::kSKmOscNuECCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin2CCQENorm+1, NEtrueBinTempl);
    double nsksignueccres_etrue_ma = skexp->GetSuperkMC()->Template(s, conv::kSKmOscNuECC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCResNorm);
    double nsksignueccres_etrue_mb = skexp->GetSuperkMC()->Template(s, conv::kSKmOscNuECC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCResNorm+1, NEtrueBinTempl);
    double nskoscnumuccqe_etrue_ma = skexp->GetSuperkMC()->Template(s, conv::kSKmOscNuMuCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCQENorm);
    double nskoscnumuccqe_etrue_mb = skexp->GetSuperkMC()->Template(s, conv::kSKmOscNuMuCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCQENorm+1, EtrueBin2CCQENorm);
    double nskoscnumuccqe_etrue_mc = skexp->GetSuperkMC()->Template(s, conv::kSKmOscNuMuCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin2CCQENorm+1, NEtrueBinTempl);
    double nskoscnumuccres_etrue_ma = skexp->GetSuperkMC()->Template(s, conv::kSKmOscNuMuCC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCResNorm);
    double nskoscnumuccres_etrue_mb = skexp->GetSuperkMC()->Template(s, conv::kSKmOscNuMuCC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCResNorm+1, NEtrueBinTempl);

    double nsksignuebarccqe_etrue_ma = skexp->GetSuperkMC()->Template(s, conv::kSKmOscNuEBarCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCQENorm);
    double nsksignuebarccqe_etrue_mb = skexp->GetSuperkMC()->Template(s, conv::kSKmOscNuEBarCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCQENorm+1, EtrueBin2CCQENorm);
    double nsksignuebarccqe_etrue_mc = skexp->GetSuperkMC()->Template(s, conv::kSKmOscNuEBarCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin2CCQENorm+1, NEtrueBinTempl);
    double nsksignuebarccres_etrue_ma = skexp->GetSuperkMC()->Template(s, conv::kSKmOscNuEBarCC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCResNorm);
    double nsksignuebarccres_etrue_mb = skexp->GetSuperkMC()->Template(s, conv::kSKmOscNuEBarCC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCResNorm+1, NEtrueBinTempl);
    double nskoscnumubarccqe_etrue_ma = skexp->GetSuperkMC()->Template(s, conv::kSKmOscNuMuBarCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCQENorm);
    double nskoscnumubarccqe_etrue_mb = skexp->GetSuperkMC()->Template(s, conv::kSKmOscNuMuBarCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCQENorm+1, EtrueBin2CCQENorm);
    double nskoscnumubarccqe_etrue_mc = skexp->GetSuperkMC()->Template(s, conv::kSKmOscNuMuBarCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin2CCQENorm+1, NEtrueBinTempl);
    double nskoscnumubarccres_etrue_ma = skexp->GetSuperkMC()->Template(s, conv::kSKmOscNuMuBarCC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCResNorm);
    double nskoscnumubarccres_etrue_mb = skexp->GetSuperkMC()->Template(s, conv::kSKmOscNuMuBarCC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCResNorm+1, NEtrueBinTempl);
  
    // spectral etrue (flux binning)
    double nsknumu_etrue_m[conv::kXsNEtrueBins], nsknumubar_etrue_m[conv::kXsNEtrueBins], nsknue_etrue_m[conv::kXsNEtrueBins], nsknuebar_etrue_m[conv::kXsNEtrueBins];
    double nsksignue_etrue_m[conv::kXsNEtrueBins], nsksignuebar_etrue_m[conv::kXsNEtrueBins], nskoscnumu_etrue_m[conv::kXsNEtrueBins], nskoscnumubar_etrue_m[conv::kXsNEtrueBins];
    for(int itb = 0; itb < NEtrueBinNuMu; itb++) {
      nsknumu_etrue_m[itb] = 0.;
      for (int im = 0; im < conv::knm; im++) {
	conv::SuperkMode_t m = conv::kSKReactionModes[s][im];
	if(m == conv::kSKmUndefined) continue;
	if(util::IsNuMu(m)) {
	  nsknumu_etrue_m[itb] += skexp->GetSuperkMC()->Template(s, m, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBinNuMu[itb]+1, EtrueBinNuMu[itb+1]);
	}
      }
    }
    for(int itb = 0; itb < NEtrueBinNuMuBar; itb++) {
      nsknumubar_etrue_m[itb] = 0.;
      for (int im = 0; im < conv::knm; im++) {
	conv::SuperkMode_t m = conv::kSKReactionModes[s][im];
	if(m == conv::kSKmUndefined) continue;
	if(util::IsNuMuBar(m)) {
	  nsknumubar_etrue_m[itb] += skexp->GetSuperkMC()->Template(s, m, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBinNuMuBar[itb]+1, EtrueBinNuMuBar[itb+1]);
	}
      }
    }
    for(int itb = 0; itb < NEtrueBinNuE; itb++) {
      nsknue_etrue_m[itb] = 0.;
      for (int im = 0; im < conv::knm; im++) {
	conv::SuperkMode_t m = conv::kSKReactionModes[s][im];
	if(m == conv::kSKmUndefined) continue;
	if(util::IsNuE(m)) {
	  nsknue_etrue_m[itb] += skexp->GetSuperkMC()->Template(s, m, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBinNuE[itb]+1, EtrueBinNuE[itb+1]);
	}
      }
    }
    for(int itb = 0; itb < NEtrueBinNuEBar; itb++) {
      nsknuebar_etrue_m[itb] = 0.;
      for (int im = 0; im < conv::knm; im++) {
	conv::SuperkMode_t m = conv::kSKReactionModes[s][im];
	if(m == conv::kSKmUndefined) continue;
	if(util::IsNuEBar(m)) {
	  nsknuebar_etrue_m[itb] += skexp->GetSuperkMC()->Template(s, m, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBinNuEBar[itb]+1, EtrueBinNuEBar[itb+1]);
	}
      }
    }
    for(int itb = 0; itb < NEtrueBinOscNuMu; itb++) {
      nskoscnumu_etrue_m[itb] = 0.;
      for (int im = 0; im < conv::knm; im++) {
	conv::SuperkMode_t m = conv::kSKReactionModes[s][im];
	if(m == conv::kSKmUndefined) continue;
	if(util::IsOscNuMu(m)) {
	  nskoscnumu_etrue_m[itb] += skexp->GetSuperkMC()->Template(s, m, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBinOscNuMu[itb]+1, EtrueBinOscNuMu[itb+1]);
	}
      }
    }
    for(int itb = 0; itb < NEtrueBinOscNuMuBar; itb++) {
      nskoscnumubar_etrue_m[itb] = 0.;
      for (int im = 0; im < conv::knm; im++) {
	conv::SuperkMode_t m = conv::kSKReactionModes[s][im];
	if(m == conv::kSKmUndefined) continue;
	if(util::IsOscNuMuBar(m)) {
	  nskoscnumubar_etrue_m[itb] += skexp->GetSuperkMC()->Template(s, m, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBinOscNuMuBar[itb]+1, EtrueBinOscNuMuBar[itb+1]);
	}
      }
    }
    for(int itb = 0; itb < NEtrueBinOscNuE; itb++) {
      nsksignue_etrue_m[itb] = 0.;
      for (int im = 0; im < conv::knm; im++) {
	conv::SuperkMode_t m = conv::kSKReactionModes[s][im];
	if(m == conv::kSKmUndefined) continue;
	if(util::IsOscNuE(m)) {
	  nsksignue_etrue_m[itb] += skexp->GetSuperkMC()->Template(s, m, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBinOscNuE[itb]+1, EtrueBinOscNuE[itb+1]);
	}
      }
    }
    for(int itb = 0; itb < NEtrueBinOscNuEBar; itb++) {
      nsksignuebar_etrue_m[itb] = 0.;
      for (int im = 0; im < conv::knm; im++) {
	conv::SuperkMode_t m = conv::kSKReactionModes[s][im];
	if(m == conv::kSKmUndefined) continue;
	if(util::IsOscNuEBar(m)) {
	  nsksignuebar_etrue_m[itb] += skexp->GetSuperkMC()->Template(s, m, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBinOscNuEBar[itb]+1, EtrueBinOscNuEBar[itb+1]);
	}
      }
    }

    //
    // calc fractional changes
    //

    // totals
    double frac_nsk_p                  = (nsk_0                  > 0) ? 100.*(nsk_p                 - nsk_0                ) / nsk_0                 : 0.;
    double frac_nsknumucc_p            = (nsknumucc_0            > 0) ? 100.*(nsknumucc_p           - nsknumucc_0          ) / nsknumucc_0           : 0.;
    double frac_nsknumunc_p            = (nsknumunc_0            > 0) ? 100.*(nsknumunc_p           - nsknumunc_0          ) / nsknumunc_0           : 0.;
    double frac_nsknumubarcc_p         = (nsknumubarcc_0         > 0) ? 100.*(nsknumubarcc_p        - nsknumubarcc_0       ) / nsknumubarcc_0        : 0.;
    double frac_nsknumubarnc_p         = (nsknumubarnc_0         > 0) ? 100.*(nsknumubarnc_p        - nsknumubarnc_0       ) / nsknumubarnc_0        : 0.;
    double frac_nsknuecc_p             = (nsknuecc_0             > 0) ? 100.*(nsknuecc_p            - nsknuecc_0           ) / nsknuecc_0            : 0.;
    double frac_nsknuenc_p             = (nsknuenc_0             > 0) ? 100.*(nsknuenc_p            - nsknuenc_0           ) / nsknuenc_0            : 0.;
    double frac_nsknuebarcc_p          = (nsknuebarcc_0          > 0) ? 100.*(nsknuebarcc_p         - nsknuebarcc_0        ) / nsknuebarcc_0         : 0.;
    double frac_nsknuebarnc_p          = (nsknuebarnc_0          > 0) ? 100.*(nsknuebarnc_p         - nsknuebarnc_0        ) / nsknuebarnc_0         : 0.;
    double frac_nsknuesig_p            = (nsknuesig_0            > 0) ? 100.*(nsknuesig_p           - nsknuesig_0          ) / nsknuesig_0           : 0.;
    double frac_nskoscnumu_p            = (nskoscnumu_0            > 0) ? 100.*(nskoscnumu_p           - nskoscnumu_0          ) / nskoscnumu_0           : 0.;

    double frac_nsknuebarsig_p            = (nsknuebarsig_0            > 0) ? 100.*(nsknuebarsig_p           - nsknuebarsig_0          ) / nsknuebarsig_0           : 0.;
    
    double frac_nskoscnumubar_p            = (nskoscnumubar_0            > 0) ? 100.*(nskoscnumubar_p           - nskoscnumubar_0          ) / nskoscnumubar_0           : 0.;





    // individual modes
    double frac_nskmodes_p[conv::knm];
    for (int im = 0; im < conv::knm; im++) {
      conv::SuperkMode_t m = conv::kSKReactionModes[s][im];
      if(m == conv::kSKmUndefined) continue;
      //   use s, m (not is,im) as indices, eg fMCErecoSpectrum[s][m]
      frac_nskmodes_p[m] = (nskmodes_0[m] > 0) ?  100.*(nskmodes_p[m] - nskmodes_0[m]) / nskmodes_0[m]  :  0.;
    }
    // spectral ccqe
    double frac_nsknumuccqe_pa         = (nsknumuccqe_0a         > 0) ? 100.*(nsknumuccqe_pa        - nsknumuccqe_0a       ) / nsknumuccqe_0a        : 0.;
    double frac_nsknumuccqe_pb         = (nsknumuccqe_0b         > 0) ? 100.*(nsknumuccqe_pb        - nsknumuccqe_0b       ) / nsknumuccqe_0b        : 0.;
    double frac_nsknumuccqe_pc         = (nsknumuccqe_0c         > 0) ? 100.*(nsknumuccqe_pc        - nsknumuccqe_0c       ) / nsknumuccqe_0c        : 0.;
    double frac_nsknumubarccqe_pa      = (nsknumubarccqe_0a      > 0) ? 100.*(nsknumubarccqe_pa     - nsknumubarccqe_0a    ) / nsknumubarccqe_0a     : 0.;
    double frac_nsknumubarccqe_pb      = (nsknumubarccqe_0b      > 0) ? 100.*(nsknumubarccqe_pb     - nsknumubarccqe_0b    ) / nsknumubarccqe_0b     : 0.;
    double frac_nsknumubarccqe_pc      = (nsknumubarccqe_0c      > 0) ? 100.*(nsknumubarccqe_pc     - nsknumubarccqe_0c    ) / nsknumubarccqe_0c     : 0.;
    double frac_nsknumu_ereco_pa         = (nsknumu_ereco_0a         > 0) ? 100.*(nsknumu_ereco_pa        - nsknumu_ereco_0a       ) / nsknumu_ereco_0a        : 0.;
    double frac_nsknumu_ereco_pb         = (nsknumu_ereco_0b         > 0) ? 100.*(nsknumu_ereco_pb        - nsknumu_ereco_0b       ) / nsknumu_ereco_0b        : 0.;
    double frac_nsknumu_ereco_pc         = (nsknumu_ereco_0c         > 0) ? 100.*(nsknumu_ereco_pc        - nsknumu_ereco_0c       ) / nsknumu_ereco_0c        : 0.;
    double frac_nsknue_ereco_pa         = (nsknue_ereco_0a         > 0) ? 100.*(nsknue_ereco_pa        - nsknue_ereco_0a       ) / nsknue_ereco_0a        : 0.;
    double frac_nsknue_ereco_pb         = (nsknue_ereco_0b         > 0) ? 100.*(nsknue_ereco_pb        - nsknue_ereco_0b       ) / nsknue_ereco_0b        : 0.;
    double frac_nsknue_ereco_pc         = (nsknue_ereco_0c         > 0) ? 100.*(nsknue_ereco_pc        - nsknue_ereco_0c       ) / nsknue_ereco_0c        : 0.;
    double frac_nskoscnue_ereco_pa         = (nskoscnue_ereco_0a         > 0) ? 100.*(nskoscnue_ereco_pa        - nskoscnue_ereco_0a       ) / nskoscnue_ereco_0a        : 0.;
    double frac_nskoscnue_ereco_pb         = (nskoscnue_ereco_0b         > 0) ? 100.*(nskoscnue_ereco_pb        - nskoscnue_ereco_0b       ) / nskoscnue_ereco_0b        : 0.;
    double frac_nskoscnue_ereco_pc         = (nskoscnue_ereco_0c         > 0) ? 100.*(nskoscnue_ereco_pc        - nskoscnue_ereco_0c       ) / nskoscnue_ereco_0c        : 0.;
    double frac_nsknc_ereco_pa         = (nsknc_ereco_0a         > 0) ? 100.*(nsknc_ereco_pa        - nsknc_ereco_0a       ) / nsknc_ereco_0a        : 0.;
    double frac_nsknc_ereco_pb         = (nsknc_ereco_0b         > 0) ? 100.*(nsknc_ereco_pb        - nsknc_ereco_0b       ) / nsknc_ereco_0b        : 0.;
    double frac_nsknc_ereco_pc         = (nsknc_ereco_0c         > 0) ? 100.*(nsknc_ereco_pc        - nsknc_ereco_0c       ) / nsknc_ereco_0c        : 0.;

    //spectral etrue (norm)
    double frac_nsknumuccqe_etrue_pa = (nsknumuccqe_etrue_0a > 0) ? 100.*(nsknumuccqe_etrue_pa - nsknumuccqe_etrue_0a) / nsknumuccqe_etrue_0a : 0.;
    double frac_nsknumuccqe_etrue_pb = (nsknumuccqe_etrue_0b > 0) ? 100.*(nsknumuccqe_etrue_pb - nsknumuccqe_etrue_0b) / nsknumuccqe_etrue_0b : 0.;
    double frac_nsknumuccqe_etrue_pc = (nsknumuccqe_etrue_0c > 0) ? 100.*(nsknumuccqe_etrue_pc - nsknumuccqe_etrue_0c) / nsknumuccqe_etrue_0c : 0.;
    double frac_nsknumuccres_etrue_pa = (nsknumuccres_etrue_0a > 0) ? 100.*(nsknumuccres_etrue_pa - nsknumuccres_etrue_0a) / nsknumuccres_etrue_0a : 0.;
    double frac_nsknumuccres_etrue_pb = (nsknumuccres_etrue_0b > 0) ? 100.*(nsknumuccres_etrue_pb - nsknumuccres_etrue_0b) / nsknumuccres_etrue_0b : 0.;
    double frac_nsknumubarccqe_etrue_pa = (nsknumubarccqe_etrue_0a > 0) ? 100.*(nsknumubarccqe_etrue_pa - nsknumubarccqe_etrue_0a) / nsknumubarccqe_etrue_0a : 0.;
    double frac_nsknumubarccqe_etrue_pb = (nsknumubarccqe_etrue_0b > 0) ? 100.*(nsknumubarccqe_etrue_pb - nsknumubarccqe_etrue_0b) / nsknumubarccqe_etrue_0b : 0.;
    double frac_nsknumubarccqe_etrue_pc = (nsknumubarccqe_etrue_0c > 0) ? 100.*(nsknumubarccqe_etrue_pc - nsknumubarccqe_etrue_0c) / nsknumubarccqe_etrue_0c : 0.;
    double frac_nsknumubarccres_etrue_pa = (nsknumubarccres_etrue_0a > 0) ? 100.*(nsknumubarccres_etrue_pa - nsknumubarccres_etrue_0a) / nsknumubarccres_etrue_0a : 0.;
    double frac_nsknumubarccres_etrue_pb = (nsknumubarccres_etrue_0b > 0) ? 100.*(nsknumubarccres_etrue_pb - nsknumubarccres_etrue_0b) / nsknumubarccres_etrue_0b : 0.;
    double frac_nsknueccqe_etrue_pa = (nsknueccqe_etrue_0a > 0) ? 100.*(nsknueccqe_etrue_pa - nsknueccqe_etrue_0a) / nsknueccqe_etrue_0a : 0.;
    double frac_nsknueccqe_etrue_pb = (nsknueccqe_etrue_0b > 0) ? 100.*(nsknueccqe_etrue_pb - nsknueccqe_etrue_0b) / nsknueccqe_etrue_0b : 0.;
    double frac_nsknueccqe_etrue_pc = (nsknueccqe_etrue_0c > 0) ? 100.*(nsknueccqe_etrue_pc - nsknueccqe_etrue_0c) / nsknueccqe_etrue_0c : 0.;
    double frac_nsknueccres_etrue_pa = (nsknueccres_etrue_0a > 0) ? 100.*(nsknueccres_etrue_pa - nsknueccres_etrue_0a) / nsknueccres_etrue_0a : 0.;
    double frac_nsknueccres_etrue_pb = (nsknueccres_etrue_0b > 0) ? 100.*(nsknueccres_etrue_pb - nsknueccres_etrue_0b) / nsknueccres_etrue_0b : 0.;
    double frac_nsknuebarccqe_etrue_pa = (nsknuebarccqe_etrue_0a > 0) ? 100.*(nsknuebarccqe_etrue_pa - nsknuebarccqe_etrue_0a) / nsknuebarccqe_etrue_0a : 0.;
    double frac_nsknuebarccqe_etrue_pb = (nsknuebarccqe_etrue_0b > 0) ? 100.*(nsknuebarccqe_etrue_pb - nsknuebarccqe_etrue_0b) / nsknuebarccqe_etrue_0b : 0.;
    double frac_nsknuebarccqe_etrue_pc = (nsknuebarccqe_etrue_0c > 0) ? 100.*(nsknuebarccqe_etrue_pc - nsknuebarccqe_etrue_0c) / nsknuebarccqe_etrue_0c : 0.;
    double frac_nsknuebarccres_etrue_pa = (nsknuebarccres_etrue_0a > 0) ? 100.*(nsknuebarccres_etrue_pa - nsknuebarccres_etrue_0a) / nsknuebarccres_etrue_0a : 0.;
    double frac_nsknuebarccres_etrue_pb = (nsknuebarccres_etrue_0b > 0) ? 100.*(nsknuebarccres_etrue_pb - nsknuebarccres_etrue_0b) / nsknuebarccres_etrue_0b : 0.;
    double frac_nsksignueccqe_etrue_pa = (nsksignueccqe_etrue_0a > 0) ? 100.*(nsksignueccqe_etrue_pa - nsksignueccqe_etrue_0a) / nsksignueccqe_etrue_0a : 0.;
    double frac_nsksignueccqe_etrue_pb = (nsksignueccqe_etrue_0b > 0) ? 100.*(nsksignueccqe_etrue_pb - nsksignueccqe_etrue_0b) / nsksignueccqe_etrue_0b : 0.;
    double frac_nsksignueccqe_etrue_pc = (nsksignueccqe_etrue_0c > 0) ? 100.*(nsksignueccqe_etrue_pc - nsksignueccqe_etrue_0c) / nsksignueccqe_etrue_0c : 0.;
    double frac_nsksignueccres_etrue_pa = (nsksignueccres_etrue_0a > 0) ? 100.*(nsksignueccres_etrue_pa - nsksignueccres_etrue_0a) / nsksignueccres_etrue_0a : 0.;
    double frac_nsksignueccres_etrue_pb = (nsksignueccres_etrue_0b > 0) ? 100.*(nsksignueccres_etrue_pb - nsksignueccres_etrue_0b) / nsksignueccres_etrue_0b : 0.;
    double frac_nskoscnumuccqe_etrue_pa = (nskoscnumuccqe_etrue_0a > 0) ? 100.*(nskoscnumuccqe_etrue_pa - nskoscnumuccqe_etrue_0a) / nskoscnumuccqe_etrue_0a : 0.;
    double frac_nskoscnumuccqe_etrue_pb = (nskoscnumuccqe_etrue_0b > 0) ? 100.*(nskoscnumuccqe_etrue_pb - nskoscnumuccqe_etrue_0b) / nskoscnumuccqe_etrue_0b : 0.;
    double frac_nskoscnumuccqe_etrue_pc = (nskoscnumuccqe_etrue_0c > 0) ? 100.*(nskoscnumuccqe_etrue_pc - nskoscnumuccqe_etrue_0c) / nskoscnumuccqe_etrue_0c : 0.;
    double frac_nskoscnumuccres_etrue_pa = (nskoscnumuccres_etrue_0a > 0) ? 100.*(nskoscnumuccres_etrue_pa - nskoscnumuccres_etrue_0a) / nskoscnumuccres_etrue_0a : 0.;
    double frac_nskoscnumuccres_etrue_pb = (nskoscnumuccres_etrue_0b > 0) ? 100.*(nskoscnumuccres_etrue_pb - nskoscnumuccres_etrue_0b) / nskoscnumuccres_etrue_0b : 0.;

    double frac_nsksignuebarccqe_etrue_pa = (nsksignuebarccqe_etrue_0a > 0) ? 100.*(nsksignuebarccqe_etrue_pa - nsksignuebarccqe_etrue_0a) / nsksignuebarccqe_etrue_0a : 0.;
    double frac_nsksignuebarccqe_etrue_pb = (nsksignuebarccqe_etrue_0b > 0) ? 100.*(nsksignuebarccqe_etrue_pb - nsksignuebarccqe_etrue_0b) / nsksignuebarccqe_etrue_0b : 0.;
    double frac_nsksignuebarccqe_etrue_pc = (nsksignuebarccqe_etrue_0c > 0) ? 100.*(nsksignuebarccqe_etrue_pc - nsksignuebarccqe_etrue_0c) / nsksignuebarccqe_etrue_0c : 0.;
    double frac_nsksignuebarccres_etrue_pa = (nsksignuebarccres_etrue_0a > 0) ? 100.*(nsksignuebarccres_etrue_pa - nsksignuebarccres_etrue_0a) / nsksignuebarccres_etrue_0a : 0.;
    double frac_nsksignuebarccres_etrue_pb = (nsksignuebarccres_etrue_0b > 0) ? 100.*(nsksignuebarccres_etrue_pb - nsksignuebarccres_etrue_0b) / nsksignuebarccres_etrue_0b : 0.;
    double frac_nskoscnumubarccqe_etrue_pa = (nskoscnumubarccqe_etrue_0a > 0) ? 100.*(nskoscnumubarccqe_etrue_pa - nskoscnumubarccqe_etrue_0a) / nskoscnumubarccqe_etrue_0a : 0.;
    double frac_nskoscnumubarccqe_etrue_pb = (nskoscnumubarccqe_etrue_0b > 0) ? 100.*(nskoscnumubarccqe_etrue_pb - nskoscnumubarccqe_etrue_0b) / nskoscnumubarccqe_etrue_0b : 0.;
    double frac_nskoscnumubarccqe_etrue_pc = (nskoscnumubarccqe_etrue_0c > 0) ? 100.*(nskoscnumubarccqe_etrue_pc - nskoscnumubarccqe_etrue_0c) / nskoscnumubarccqe_etrue_0c : 0.;
    double frac_nskoscnumubarccres_etrue_pa = (nskoscnumubarccres_etrue_0a > 0) ? 100.*(nskoscnumubarccres_etrue_pa - nskoscnumubarccres_etrue_0a) / nskoscnumubarccres_etrue_0a : 0.;
    double frac_nskoscnumubarccres_etrue_pb = (nskoscnumubarccres_etrue_0b > 0) ? 100.*(nskoscnumubarccres_etrue_pb - nskoscnumubarccres_etrue_0b) / nskoscnumubarccres_etrue_0b : 0.;
  
    // spectral etrue (flux binning)
    double frac_nsknumu_etrue_p[conv::kXsNEtrueBins], frac_nsknumubar_etrue_p[conv::kXsNEtrueBins], frac_nsknue_etrue_p[conv::kXsNEtrueBins], frac_nsknuebar_etrue_p[conv::kXsNEtrueBins];
    double frac_nsksignue_etrue_p[conv::kXsNEtrueBins], frac_nsksignuebar_etrue_p[conv::kXsNEtrueBins], frac_nskoscnumu_etrue_p[conv::kXsNEtrueBins], frac_nskoscnumubar_etrue_p[conv::kXsNEtrueBins];
    for(int itb = 0; itb < NEtrueBinNuMu; itb++) {
      frac_nsknumu_etrue_p[itb]        = (nsknumu_etrue_0[itb]    > 0) ? 100.*(nsknumu_etrue_p[itb]   - nsknumu_etrue_0[itb]  ) / nsknumu_etrue_0[itb]   : 0.;
    }
    for(int itb = 0; itb < NEtrueBinNuMuBar; itb++) {
      frac_nsknumubar_etrue_p[itb]        = (nsknumubar_etrue_0[itb]    > 0) ? 100.*(nsknumubar_etrue_p[itb]   - nsknumubar_etrue_0[itb]  ) / nsknumubar_etrue_0[itb]   : 0.;
    }
    for(int itb = 0; itb < NEtrueBinNuE; itb++) {
      frac_nsknue_etrue_p[itb]        = (nsknue_etrue_0[itb]    > 0) ? 100.*(nsknue_etrue_p[itb]   - nsknue_etrue_0[itb]  ) / nsknue_etrue_0[itb]   : 0.;
    }
    for(int itb = 0; itb < NEtrueBinNuEBar; itb++) {
      frac_nsknuebar_etrue_p[itb]        = (nsknuebar_etrue_0[itb]    > 0) ? 100.*(nsknuebar_etrue_p[itb]   - nsknuebar_etrue_0[itb]  ) / nsknuebar_etrue_0[itb]   : 0.;
    }
    for(int itb = 0; itb < NEtrueBinOscNuMu; itb++) {
      frac_nskoscnumu_etrue_p[itb]        = (nskoscnumu_etrue_0[itb]    > 0) ? 100.*(nskoscnumu_etrue_p[itb]   - nskoscnumu_etrue_0[itb]  ) / nskoscnumu_etrue_0[itb]   : 0.;
    }
    for(int itb = 0; itb < NEtrueBinOscNuMuBar; itb++) {
      frac_nskoscnumubar_etrue_p[itb]        = (nskoscnumubar_etrue_0[itb]    > 0) ? 100.*(nskoscnumubar_etrue_p[itb]   - nskoscnumubar_etrue_0[itb]  ) / nskoscnumubar_etrue_0[itb]   : 0.;
    }
    for(int itb = 0; itb < NEtrueBinOscNuE; itb++) {
      frac_nsksignue_etrue_p[itb]        = (nsksignue_etrue_0[itb]    > 0) ? 100.*(nsksignue_etrue_p[itb]   - nsksignue_etrue_0[itb]  ) / nsksignue_etrue_0[itb]   : 0.;
    }
    for(int itb = 0; itb < NEtrueBinOscNuEBar; itb++) {
      frac_nsksignuebar_etrue_p[itb]        = (nsksignuebar_etrue_0[itb]    > 0) ? 100.*(nsksignuebar_etrue_p[itb]   - nsksignuebar_etrue_0[itb]  ) / nsksignuebar_etrue_0[itb]   : 0.;
    }



    // totals
    double frac_nsk_m                  = (nsk_0                  > 0) ? 100.*(nsk_m                 - nsk_0                ) / nsk_0                 : 0.;
    double frac_nsknumucc_m            = (nsknumucc_0            > 0) ? 100.*(nsknumucc_m           - nsknumucc_0          ) / nsknumucc_0           : 0.;
    double frac_nsknumunc_m            = (nsknumunc_0            > 0) ? 100.*(nsknumunc_m           - nsknumunc_0          ) / nsknumunc_0           : 0.;
    double frac_nsknumubarcc_m         = (nsknumubarcc_0         > 0) ? 100.*(nsknumubarcc_m        - nsknumubarcc_0       ) / nsknumubarcc_0        : 0.;
    double frac_nsknumubarnc_m         = (nsknumubarnc_0         > 0) ? 100.*(nsknumubarnc_m        - nsknumubarnc_0       ) / nsknumubarnc_0        : 0.;
    double frac_nsknuecc_m             = (nsknuecc_0             > 0) ? 100.*(nsknuecc_m            - nsknuecc_0           ) / nsknuecc_0            : 0.;
    double frac_nsknuenc_m             = (nsknuenc_0             > 0) ? 100.*(nsknuenc_m            - nsknuenc_0           ) / nsknuenc_0            : 0.;
    double frac_nsknuebarcc_m          = (nsknuebarcc_0          > 0) ? 100.*(nsknuebarcc_m         - nsknuebarcc_0        ) / nsknuebarcc_0         : 0.;
    double frac_nsknuebarnc_m          = (nsknuebarnc_0          > 0) ? 100.*(nsknuebarnc_m         - nsknuebarnc_0        ) / nsknuebarnc_0         : 0.;
    double frac_nsknuesig_m            = (nsknuesig_0            > 0) ? 100.*(nsknuesig_m           - nsknuesig_0          ) / nsknuesig_0           : 0.;
    double frac_nskoscnumu_m           = (nskoscnumu_0           > 0) ? 100.*(nskoscnumu_m           - nskoscnumu_0          ) / nskoscnumu_0           : 0.;

    double frac_nsknuebarsig_m            = (nsknuebarsig_0            > 0) ? 100.*(nsknuebarsig_m           - nsknuebarsig_0          ) / nsknuebarsig_0           : 0.;
    double frac_nskoscnumubar_m           = (nskoscnumubar_0           > 0) ? 100.*(nskoscnumubar_m           - nskoscnumubar_0          ) / nskoscnumubar_0           : 0.;

    // individual modes
    double frac_nskmodes_m[conv::knm];
    for (int im = 0; im < conv::knm; im++) {
      conv::SuperkMode_t m = conv::kSKReactionModes[s][im];
      if(m == conv::kSKmUndefined) continue;
      //   use s, m (not is,im) as indices, eg fMCErecoSpectrum[s][m]
      frac_nskmodes_m[m] = (nskmodes_0[m] > 0) ?  100.*(nskmodes_m[m] - nskmodes_0[m]) / nskmodes_0[m]  :  0.;
    }

    // spectral ccqe
    double frac_nsknumuccqe_ma         = (nsknumuccqe_0a         > 0) ? 100.*(nsknumuccqe_ma        - nsknumuccqe_0a       ) / nsknumuccqe_0a        : 0.;
    double frac_nsknumuccqe_mb         = (nsknumuccqe_0b         > 0) ? 100.*(nsknumuccqe_mb        - nsknumuccqe_0b       ) / nsknumuccqe_0b        : 0.;
    double frac_nsknumuccqe_mc         = (nsknumuccqe_0c         > 0) ? 100.*(nsknumuccqe_mc        - nsknumuccqe_0c       ) / nsknumuccqe_0c        : 0.;
    double frac_nsknumubarccqe_ma      = (nsknumubarccqe_0a      > 0) ? 100.*(nsknumubarccqe_ma     - nsknumubarccqe_0a    ) / nsknumubarccqe_0a     : 0.;
    double frac_nsknumubarccqe_mb      = (nsknumubarccqe_0b      > 0) ? 100.*(nsknumubarccqe_mb     - nsknumubarccqe_0b    ) / nsknumubarccqe_0b     : 0.;
    double frac_nsknumubarccqe_mc      = (nsknumubarccqe_0c      > 0) ? 100.*(nsknumubarccqe_mc     - nsknumubarccqe_0c    ) / nsknumubarccqe_0c     : 0.;
    double frac_nsknumu_ereco_ma         = (nsknumu_ereco_0a         > 0) ? 100.*(nsknumu_ereco_ma        - nsknumu_ereco_0a       ) / nsknumu_ereco_0a        : 0.;
    double frac_nsknumu_ereco_mb         = (nsknumu_ereco_0b         > 0) ? 100.*(nsknumu_ereco_mb        - nsknumu_ereco_0b       ) / nsknumu_ereco_0b        : 0.;
    double frac_nsknumu_ereco_mc         = (nsknumu_ereco_0c         > 0) ? 100.*(nsknumu_ereco_mc        - nsknumu_ereco_0c       ) / nsknumu_ereco_0c        : 0.;
    double frac_nsknue_ereco_ma         = (nsknue_ereco_0a         > 0) ? 100.*(nsknue_ereco_ma        - nsknue_ereco_0a       ) / nsknue_ereco_0a        : 0.;
    double frac_nsknue_ereco_mb         = (nsknue_ereco_0b         > 0) ? 100.*(nsknue_ereco_mb        - nsknue_ereco_0b       ) / nsknue_ereco_0b        : 0.;
    double frac_nsknue_ereco_mc         = (nsknue_ereco_0c         > 0) ? 100.*(nsknue_ereco_mc        - nsknue_ereco_0c       ) / nsknue_ereco_0c        : 0.;
    double frac_nskoscnue_ereco_ma         = (nskoscnue_ereco_0a         > 0) ? 100.*(nskoscnue_ereco_ma        - nskoscnue_ereco_0a       ) / nskoscnue_ereco_0a        : 0.;
    double frac_nskoscnue_ereco_mb         = (nskoscnue_ereco_0b         > 0) ? 100.*(nskoscnue_ereco_mb        - nskoscnue_ereco_0b       ) / nskoscnue_ereco_0b        : 0.;
    double frac_nskoscnue_ereco_mc         = (nskoscnue_ereco_0c         > 0) ? 100.*(nskoscnue_ereco_mc        - nskoscnue_ereco_0c       ) / nskoscnue_ereco_0c        : 0.;
    double frac_nsknc_ereco_ma         = (nsknc_ereco_0a         > 0) ? 100.*(nsknc_ereco_ma        - nsknc_ereco_0a       ) / nsknc_ereco_0a        : 0.;
    double frac_nsknc_ereco_mb         = (nsknc_ereco_0b         > 0) ? 100.*(nsknc_ereco_mb        - nsknc_ereco_0b       ) / nsknc_ereco_0b        : 0.;
    double frac_nsknc_ereco_mc         = (nsknc_ereco_0c         > 0) ? 100.*(nsknc_ereco_mc        - nsknc_ereco_0c       ) / nsknc_ereco_0c        : 0.;

    //spectral etrue (norm)
    double frac_nsknumuccqe_etrue_ma = (nsknumuccqe_etrue_0a > 0) ? 100.*(nsknumuccqe_etrue_ma - nsknumuccqe_etrue_0a) / nsknumuccqe_etrue_0a : 0.;
    double frac_nsknumuccqe_etrue_mb = (nsknumuccqe_etrue_0b > 0) ? 100.*(nsknumuccqe_etrue_mb - nsknumuccqe_etrue_0b) / nsknumuccqe_etrue_0b : 0.;
    double frac_nsknumuccqe_etrue_mc = (nsknumuccqe_etrue_0c > 0) ? 100.*(nsknumuccqe_etrue_mc - nsknumuccqe_etrue_0c) / nsknumuccqe_etrue_0c : 0.;
    double frac_nsknumuccres_etrue_ma = (nsknumuccres_etrue_0a > 0) ? 100.*(nsknumuccres_etrue_ma - nsknumuccres_etrue_0a) / nsknumuccres_etrue_0a : 0.;
    double frac_nsknumuccres_etrue_mb = (nsknumuccres_etrue_0b > 0) ? 100.*(nsknumuccres_etrue_mb - nsknumuccres_etrue_0b) / nsknumuccres_etrue_0b : 0.;
    double frac_nsknumubarccqe_etrue_ma = (nsknumubarccqe_etrue_0a > 0) ? 100.*(nsknumubarccqe_etrue_ma - nsknumubarccqe_etrue_0a) / nsknumubarccqe_etrue_0a : 0.;
    double frac_nsknumubarccqe_etrue_mb = (nsknumubarccqe_etrue_0b > 0) ? 100.*(nsknumubarccqe_etrue_mb - nsknumubarccqe_etrue_0b) / nsknumubarccqe_etrue_0b : 0.;
    double frac_nsknumubarccqe_etrue_mc = (nsknumubarccqe_etrue_0c > 0) ? 100.*(nsknumubarccqe_etrue_mc - nsknumubarccqe_etrue_0c) / nsknumubarccqe_etrue_0c : 0.;
    double frac_nsknumubarccres_etrue_ma = (nsknumubarccres_etrue_0a > 0) ? 100.*(nsknumubarccres_etrue_ma - nsknumubarccres_etrue_0a) / nsknumubarccres_etrue_0a : 0.;
    double frac_nsknumubarccres_etrue_mb = (nsknumubarccres_etrue_0b > 0) ? 100.*(nsknumubarccres_etrue_mb - nsknumubarccres_etrue_0b) / nsknumubarccres_etrue_0b : 0.;
    double frac_nsknueccqe_etrue_ma = (nsknueccqe_etrue_0a > 0) ? 100.*(nsknueccqe_etrue_ma - nsknueccqe_etrue_0a) / nsknueccqe_etrue_0a : 0.;
    double frac_nsknueccqe_etrue_mb = (nsknueccqe_etrue_0b > 0) ? 100.*(nsknueccqe_etrue_mb - nsknueccqe_etrue_0b) / nsknueccqe_etrue_0b : 0.;
    double frac_nsknueccqe_etrue_mc = (nsknueccqe_etrue_0c > 0) ? 100.*(nsknueccqe_etrue_mc - nsknueccqe_etrue_0c) / nsknueccqe_etrue_0c : 0.;
    double frac_nsknueccres_etrue_ma = (nsknueccres_etrue_0a > 0) ? 100.*(nsknueccres_etrue_ma - nsknueccres_etrue_0a) / nsknueccres_etrue_0a : 0.;
    double frac_nsknueccres_etrue_mb = (nsknueccres_etrue_0b > 0) ? 100.*(nsknueccres_etrue_mb - nsknueccres_etrue_0b) / nsknueccres_etrue_0b : 0.;
    double frac_nsknuebarccqe_etrue_ma = (nsknuebarccqe_etrue_0a > 0) ? 100.*(nsknuebarccqe_etrue_ma - nsknuebarccqe_etrue_0a) / nsknuebarccqe_etrue_0a : 0.;
    double frac_nsknuebarccqe_etrue_mb = (nsknuebarccqe_etrue_0b > 0) ? 100.*(nsknuebarccqe_etrue_mb - nsknuebarccqe_etrue_0b) / nsknuebarccqe_etrue_0b : 0.;
    double frac_nsknuebarccqe_etrue_mc = (nsknuebarccqe_etrue_0c > 0) ? 100.*(nsknuebarccqe_etrue_mc - nsknuebarccqe_etrue_0c) / nsknuebarccqe_etrue_0c : 0.;
    double frac_nsknuebarccres_etrue_ma = (nsknuebarccres_etrue_0a > 0) ? 100.*(nsknuebarccres_etrue_ma - nsknuebarccres_etrue_0a) / nsknuebarccres_etrue_0a : 0.;
    double frac_nsknuebarccres_etrue_mb = (nsknuebarccres_etrue_0b > 0) ? 100.*(nsknuebarccres_etrue_mb - nsknuebarccres_etrue_0b) / nsknuebarccres_etrue_0b : 0.;
    double frac_nsksignueccqe_etrue_ma = (nsksignueccqe_etrue_0a > 0) ? 100.*(nsksignueccqe_etrue_ma - nsksignueccqe_etrue_0a) / nsksignueccqe_etrue_0a : 0.;
    double frac_nsksignueccqe_etrue_mb = (nsksignueccqe_etrue_0b > 0) ? 100.*(nsksignueccqe_etrue_mb - nsksignueccqe_etrue_0b) / nsksignueccqe_etrue_0b : 0.;
    double frac_nsksignueccqe_etrue_mc = (nsksignueccqe_etrue_0c > 0) ? 100.*(nsksignueccqe_etrue_mc - nsksignueccqe_etrue_0c) / nsksignueccqe_etrue_0c : 0.;
    double frac_nsksignueccres_etrue_ma = (nsksignueccres_etrue_0a > 0) ? 100.*(nsksignueccres_etrue_ma - nsksignueccres_etrue_0a) / nsksignueccres_etrue_0a : 0.;
    double frac_nsksignueccres_etrue_mb = (nsksignueccres_etrue_0b > 0) ? 100.*(nsksignueccres_etrue_mb - nsksignueccres_etrue_0b) / nsksignueccres_etrue_0b : 0.;
    double frac_nskoscnumuccqe_etrue_ma = (nskoscnumuccqe_etrue_0a > 0) ? 100.*(nskoscnumuccqe_etrue_ma - nskoscnumuccqe_etrue_0a) / nskoscnumuccqe_etrue_0a : 0.;
    double frac_nskoscnumuccqe_etrue_mb = (nskoscnumuccqe_etrue_0b > 0) ? 100.*(nskoscnumuccqe_etrue_mb - nskoscnumuccqe_etrue_0b) / nskoscnumuccqe_etrue_0b : 0.;
    double frac_nskoscnumuccqe_etrue_mc = (nskoscnumuccqe_etrue_0c > 0) ? 100.*(nskoscnumuccqe_etrue_mc - nskoscnumuccqe_etrue_0c) / nskoscnumuccqe_etrue_0c : 0.;
    double frac_nskoscnumuccres_etrue_ma = (nskoscnumuccres_etrue_0a > 0) ? 100.*(nskoscnumuccres_etrue_ma - nskoscnumuccres_etrue_0a) / nskoscnumuccres_etrue_0a : 0.;
    double frac_nskoscnumuccres_etrue_mb = (nskoscnumuccres_etrue_0b > 0) ? 100.*(nskoscnumuccres_etrue_mb - nskoscnumuccres_etrue_0b) / nskoscnumuccres_etrue_0b : 0.;

    double frac_nsksignuebarccqe_etrue_ma = (nsksignuebarccqe_etrue_0a > 0) ? 100.*(nsksignuebarccqe_etrue_ma - nsksignuebarccqe_etrue_0a) / nsksignuebarccqe_etrue_0a : 0.;
    double frac_nsksignuebarccqe_etrue_mb = (nsksignuebarccqe_etrue_0b > 0) ? 100.*(nsksignuebarccqe_etrue_mb - nsksignuebarccqe_etrue_0b) / nsksignuebarccqe_etrue_0b : 0.;
    double frac_nsksignuebarccqe_etrue_mc = (nsksignuebarccqe_etrue_0c > 0) ? 100.*(nsksignuebarccqe_etrue_mc - nsksignuebarccqe_etrue_0c) / nsksignuebarccqe_etrue_0c : 0.;
    double frac_nsksignuebarccres_etrue_ma = (nsksignuebarccres_etrue_0a > 0) ? 100.*(nsksignuebarccres_etrue_ma - nsksignuebarccres_etrue_0a) / nsksignuebarccres_etrue_0a : 0.;
    double frac_nsksignuebarccres_etrue_mb = (nsksignuebarccres_etrue_0b > 0) ? 100.*(nsksignuebarccres_etrue_mb - nsksignuebarccres_etrue_0b) / nsksignuebarccres_etrue_0b : 0.;
    double frac_nskoscnumubarccqe_etrue_ma = (nskoscnumubarccqe_etrue_0a > 0) ? 100.*(nskoscnumubarccqe_etrue_ma - nskoscnumubarccqe_etrue_0a) / nskoscnumubarccqe_etrue_0a : 0.;
    double frac_nskoscnumubarccqe_etrue_mb = (nskoscnumubarccqe_etrue_0b > 0) ? 100.*(nskoscnumubarccqe_etrue_mb - nskoscnumubarccqe_etrue_0b) / nskoscnumubarccqe_etrue_0b : 0.;
    double frac_nskoscnumubarccqe_etrue_mc = (nskoscnumubarccqe_etrue_0c > 0) ? 100.*(nskoscnumubarccqe_etrue_mc - nskoscnumubarccqe_etrue_0c) / nskoscnumubarccqe_etrue_0c : 0.;
    double frac_nskoscnumubarccres_etrue_ma = (nskoscnumubarccres_etrue_0a > 0) ? 100.*(nskoscnumubarccres_etrue_ma - nskoscnumubarccres_etrue_0a) / nskoscnumubarccres_etrue_0a : 0.;
    double frac_nskoscnumubarccres_etrue_mb = (nskoscnumubarccres_etrue_0b > 0) ? 100.*(nskoscnumubarccres_etrue_mb - nskoscnumubarccres_etrue_0b) / nskoscnumubarccres_etrue_0b : 0.;



  
    // spectral etrue (flux binning)
    double frac_nsknumu_etrue_m[conv::kXsNEtrueBins], frac_nsknumubar_etrue_m[conv::kXsNEtrueBins], frac_nsknue_etrue_m[conv::kXsNEtrueBins], frac_nsknuebar_etrue_m[conv::kXsNEtrueBins], frac_nsksignue_etrue_m[conv::kXsNEtrueBins], frac_nsksignuebar_etrue_m[conv::kXsNEtrueBins], frac_nskoscnumu_etrue_m[conv::kXsNEtrueBins], frac_nskoscnumubar_etrue_m[conv::kXsNEtrueBins];
    for(int itb = 0; itb < NEtrueBinNuMu; itb++) {
      frac_nsknumu_etrue_m[itb]        = (nsknumu_etrue_0[itb]    > 0) ? 100.*(nsknumu_etrue_m[itb]   - nsknumu_etrue_0[itb]  ) / nsknumu_etrue_0[itb]   : 0.;
    }
    for(int itb = 0; itb < NEtrueBinNuMuBar; itb++) {
      frac_nsknumubar_etrue_m[itb]        = (nsknumubar_etrue_0[itb]    > 0) ? 100.*(nsknumubar_etrue_m[itb]   - nsknumubar_etrue_0[itb]  ) / nsknumubar_etrue_0[itb]   : 0.;
    }
    for(int itb = 0; itb < NEtrueBinNuE; itb++) {
      frac_nsknue_etrue_m[itb]        = (nsknue_etrue_0[itb]    > 0) ? 100.*(nsknue_etrue_m[itb]   - nsknue_etrue_0[itb]  ) / nsknue_etrue_0[itb]   : 0.;
    }
    for(int itb = 0; itb < NEtrueBinNuEBar; itb++) {
      frac_nsknuebar_etrue_m[itb]        = (nsknuebar_etrue_0[itb]    > 0) ? 100.*(nsknuebar_etrue_m[itb]   - nsknuebar_etrue_0[itb]  ) / nsknuebar_etrue_0[itb]   : 0.;
    }
    for(int itb = 0; itb < NEtrueBinOscNuMu; itb++) {
      frac_nskoscnumu_etrue_m[itb]        = (nskoscnumu_etrue_0[itb]    > 0) ? 100.*(nskoscnumu_etrue_m[itb]   - nskoscnumu_etrue_0[itb]  ) / nskoscnumu_etrue_0[itb]   : 0.;
    }
    for(int itb = 0; itb < NEtrueBinOscNuMuBar; itb++) {
      frac_nskoscnumubar_etrue_m[itb]        = (nskoscnumubar_etrue_0[itb]    > 0) ? 100.*(nskoscnumubar_etrue_m[itb]   - nskoscnumubar_etrue_0[itb]  ) / nskoscnumubar_etrue_0[itb]   : 0.;
    }
    for(int itb = 0; itb < NEtrueBinOscNuE; itb++) {
      frac_nsksignue_etrue_m[itb]        = (nsksignue_etrue_0[itb]    > 0) ? 100.*(nsksignue_etrue_m[itb]   - nsksignue_etrue_0[itb]  ) / nsksignue_etrue_0[itb]   : 0.;
    }
    for(int itb = 0; itb < NEtrueBinOscNuEBar; itb++) {
      frac_nsksignuebar_etrue_m[itb]        = (nsksignuebar_etrue_0[itb]    > 0) ? 100.*(nsksignuebar_etrue_m[itb]   - nsksignuebar_etrue_0[itb]  ) / nsksignuebar_etrue_0[itb]   : 0.;
    }


    //
    // Print-out
    //
    
    string pctsymbol  = (use_tex) ? " \\% "   : "% ";
    string endcolumn  = (use_tex) ? " & "     : "" ;
    string endrow     = (use_tex) ? " \\\\ "  : "" ;

    int ncol1   = 80;
    int ncolnev = 15;
    int ncolfrc =  8;

    cout << "\n\nSystematic: " << systematic << endl;
    cout << endl;
    cout << fixed;

    stringstream streamer;

    if(use_tex) {
      streamer << "\\begin{table}"   << endl;
      streamer << "\\tiny" << endl;
      //streamer << "\\centering"     << endl;
      streamer << "\\begin{tabular}{| l | r | r | r  r | r  r |}" << endl;
      streamer << "\\hline \\hline" << endl;
      streamer << "{\\bf Component} &  "
	       << "{\\bf Nominal}   &  "
	       << "{\\bf BANFF tuned} &  "
	       << "\\multicolumn {2}{c|} { {\\bf BANFF tuned +1$\\sigma$} } &  "
	       << "\\multicolumn {2}{c|} { {\\bf BANFF tuned -1$\\sigma$} } \\\\"
	       << endl;
      streamer << "\\hline" << endl;
    }
  
    string temp = "$N_{SK}$ (" + util::AsStringTex(s) + ")";
    streamer << setw(ncol1)  << setfill(' ') << temp
	     << endcolumn
	     << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsk_n 
	     << endcolumn
	     << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsk_0 
	     << endcolumn
	     << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsk_p
	     << endcolumn
	     << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsk_p << pctsymbol 
	     << endcolumn
	     << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsk_m 
	     << endcolumn
	     << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsk_m << pctsymbol 
	     << endrow
	     << endl;

    if(use_tex) {
      streamer << "\\hline" << endl;
    } else {
      streamer << setw(150)  << setfill('-') << " " << endl;
    }

    for (int im = 0; im < conv::knm; im++) {
      conv::SuperkMode_t m = conv::kSKReactionModes[s][im];
      if(m == conv::kSKmUndefined) continue;
      //get the title of the row (or part of it for E binned rows)
      string component = "$N_{SK}$ (";
      component += util::AsStringTex(s);
      component += "): ";
      component +=  util::AsStringTex(m);
      //extra hline above ccqe modes
      if(use_tex)
	if(m == conv::kSKmNuMuCCQE || m == conv::kSKmNuMuBarCCQE || m == conv::kSKmNuECCQE || m == conv::kSKmNuEBarCCQE || m == conv::kSKmOscNuECCQE)
	  streamer << "\\hline" << endl;
      ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // spectral etrue ccqe and ccres (all species)
      ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //spectral ccqe etrue (numu)
      if ( m == conv::kSKmNuMuCCQE && show_sk_spectral_info_etrue_norm_binning) {
	temp = component + " ($E_{t} <$ 1.5 GeV)";
	streamer << setw(ncol1)  << setfill(' ') << temp
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumuccqe_etrue_na
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumuccqe_etrue_0a
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumuccqe_etrue_pa
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknumuccqe_etrue_pa << pctsymbol 
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumuccqe_etrue_ma 
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknumuccqe_etrue_ma << pctsymbol 
		 << endrow
		 << endl;
	temp = component + " (1.5 $\\leq E_{t} <$ 3.5 GeV)";
	streamer << setw(ncol1)  << setfill(' ') << temp
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumuccqe_etrue_nb
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumuccqe_etrue_0b
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumuccqe_etrue_pb
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknumuccqe_etrue_pb << pctsymbol 
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumuccqe_etrue_mb 
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknumuccqe_etrue_mb << pctsymbol 
		 << endrow
		 << endl;
	temp = component + " ($E_{t} \\geq$ 3.5 GeV)";
	streamer << setw(ncol1)  << setfill(' ') << temp
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumuccqe_etrue_nc
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumuccqe_etrue_0c
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumuccqe_etrue_pc
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknumuccqe_etrue_pc << pctsymbol 
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumuccqe_etrue_mc 
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknumuccqe_etrue_mc << pctsymbol 
		 << endrow
		 << endl;
	if(use_tex) streamer << "\\hline" << endl;
      }
      //spectral ccqe etrue (numubar)
      if ( m == conv::kSKmNuMuBarCCQE && show_sk_spectral_info_etrue_norm_binning) {
	temp = component + " ($E_{t} <$ 1.5 GeV)";
	streamer << setw(ncol1)  << setfill(' ') << temp
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumubarccqe_etrue_na
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumubarccqe_etrue_0a
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumubarccqe_etrue_pa
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknumubarccqe_etrue_pa << pctsymbol 
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumubarccqe_etrue_ma 
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknumubarccqe_etrue_ma << pctsymbol 
		 << endrow
		 << endl;
	temp = component + " (1.5 $\\leq E_{t} <$ 3.5 GeV)";
	streamer << setw(ncol1)  << setfill(' ') << temp
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumubarccqe_etrue_nb
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumubarccqe_etrue_0b
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumubarccqe_etrue_pb
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknumubarccqe_etrue_pb << pctsymbol 
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumubarccqe_etrue_mb 
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknumubarccqe_etrue_mb << pctsymbol 
		 << endrow
		 << endl;
	temp = component + " ($E_{t} \\geq$ 3.5 GeV)";
	streamer << setw(ncol1)  << setfill(' ') << temp
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumubarccqe_etrue_nc
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumubarccqe_etrue_0c
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumubarccqe_etrue_pc
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknumubarccqe_etrue_pc << pctsymbol 
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumubarccqe_etrue_mc 
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknumubarccqe_etrue_mc << pctsymbol 
		 << endrow
		 << endl;
	if(use_tex) streamer << "\\hline" << endl;
      }
      //spectral ccqe etrue (nue)
      if ( m == conv::kSKmNuECCQE && show_sk_spectral_info_etrue_norm_binning) {
	temp = component + " ($E_{t} <$ 1.5 GeV)";
	streamer << setw(ncol1)  << setfill(' ') << temp
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknueccqe_etrue_na
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknueccqe_etrue_0a
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknueccqe_etrue_pa
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknueccqe_etrue_pa << pctsymbol 
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknueccqe_etrue_ma 
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknueccqe_etrue_ma << pctsymbol 
		 << endrow
		 << endl;
	temp = component + " (1.5 $\\leq E_{t} <$ 3.5 GeV)";
	streamer << setw(ncol1)  << setfill(' ') << temp
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknueccqe_etrue_nb
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknueccqe_etrue_0b
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknueccqe_etrue_pb
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknueccqe_etrue_pb << pctsymbol 
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknueccqe_etrue_mb 
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknueccqe_etrue_mb << pctsymbol 
		 << endrow
		 << endl;
	temp = component + " ($E_{t} \\geq$ 3.5 GeV)";
	streamer << setw(ncol1)  << setfill(' ') << temp
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknueccqe_etrue_nc
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknueccqe_etrue_0c
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknueccqe_etrue_pc
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknueccqe_etrue_pc << pctsymbol 
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknueccqe_etrue_mc 
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknueccqe_etrue_mc << pctsymbol 
		 << endrow
		 << endl;
	if(use_tex) streamer << "\\hline" << endl;
      }
      //spectral ccqe etrue (nuebar)
      if ( m == conv::kSKmNuEBarCCQE && show_sk_spectral_info_etrue_norm_binning) {
	temp = component + " ($E_{t} <$ 1.5 GeV)";
	streamer << setw(ncol1)  << setfill(' ') << temp
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknuebarccqe_etrue_na
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknuebarccqe_etrue_0a
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknuebarccqe_etrue_pa
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknuebarccqe_etrue_pa << pctsymbol 
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknuebarccqe_etrue_ma 
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknuebarccqe_etrue_ma << pctsymbol 
		 << endrow
		 << endl;
	temp = component + " (1.5 $\\leq E_{t} <$ 3.5 GeV)";
	streamer << setw(ncol1)  << setfill(' ') << temp
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknuebarccqe_etrue_nb
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknuebarccqe_etrue_0b
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknuebarccqe_etrue_pb
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknuebarccqe_etrue_pb << pctsymbol 
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknuebarccqe_etrue_mb 
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknuebarccqe_etrue_mb << pctsymbol 
		 << endrow
		 << endl;
	temp = component + " ($E_{t} \\geq$ 3.5 GeV)";
	streamer << setw(ncol1)  << setfill(' ') << temp
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknuebarccqe_etrue_nc
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknuebarccqe_etrue_0c
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknuebarccqe_etrue_pc
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknuebarccqe_etrue_pc << pctsymbol 
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknuebarccqe_etrue_mc 
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknuebarccqe_etrue_mc << pctsymbol 
		 << endrow
		 << endl;
	if(use_tex) streamer << "\\hline" << endl;
      }
      //spectral ccqe etrue (sig nue)
      if ( m == conv::kSKmOscNuECCQE && show_sk_spectral_info_etrue_norm_binning) {
	temp = component + " ($E_{t} <$ 1.5 GeV)";
	streamer << setw(ncol1)  << setfill(' ') << temp
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsksignueccqe_etrue_na
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsksignueccqe_etrue_0a
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsksignueccqe_etrue_pa
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsksignueccqe_etrue_pa << pctsymbol 
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsksignueccqe_etrue_ma 
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsksignueccqe_etrue_ma << pctsymbol 
		 << endrow
		 << endl;
	temp = component + " (1.5 $\\leq E_{t} <$ 3.5 GeV)";
	streamer << setw(ncol1)  << setfill(' ') << temp
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsksignueccqe_etrue_nb
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsksignueccqe_etrue_0b
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsksignueccqe_etrue_pb
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsksignueccqe_etrue_pb << pctsymbol 
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsksignueccqe_etrue_mb 
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsksignueccqe_etrue_mb << pctsymbol 
		 << endrow
		 << endl;
	temp = component + " ($E_{t} \\geq$ 3.5 GeV)";
	streamer << setw(ncol1)  << setfill(' ') << temp
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsksignueccqe_etrue_nc
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsksignueccqe_etrue_0c
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsksignueccqe_etrue_pc
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsksignueccqe_etrue_pc << pctsymbol 
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsksignueccqe_etrue_mc 
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsksignueccqe_etrue_mc << pctsymbol 
		 << endrow
		 << endl;
	if(use_tex) streamer << "\\hline" << endl;
      }

    //spectral ccqe etrue (sig nuebar)
      if ( m == conv::kSKmOscNuEBarCCQE && show_sk_spectral_info_etrue_norm_binning) {
	temp = component + " ($E_{t} <$ 1.5 GeV)";
	streamer << setw(ncol1)  << setfill(' ') << temp
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsksignuebarccqe_etrue_na
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsksignuebarccqe_etrue_0a
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsksignuebarccqe_etrue_pa
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsksignuebarccqe_etrue_pa << pctsymbol 
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsksignuebarccqe_etrue_ma 
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsksignuebarccqe_etrue_ma << pctsymbol 
		 << endrow
		 << endl;
	temp = component + " (1.5 $\\leq E_{t} <$ 3.5 GeV)";
	streamer << setw(ncol1)  << setfill(' ') << temp
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsksignuebarccqe_etrue_nb
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsksignuebarccqe_etrue_0b
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsksignuebarccqe_etrue_pb
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsksignuebarccqe_etrue_pb << pctsymbol 
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsksignuebarccqe_etrue_mb 
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsksignuebarccqe_etrue_mb << pctsymbol 
		 << endrow
		 << endl;
	temp = component + " ($E_{t} \\geq$ 3.5 GeV)";
	streamer << setw(ncol1)  << setfill(' ') << temp
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsksignuebarccqe_etrue_nc
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsksignuebarccqe_etrue_0c
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsksignuebarccqe_etrue_pc
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsksignuebarccqe_etrue_pc << pctsymbol 
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsksignuebarccqe_etrue_mc 
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsksignuebarccqe_etrue_mc << pctsymbol 
		 << endrow
		 << endl;
	if(use_tex) streamer << "\\hline" << endl;
      }      

	//spectral ccqe etrue (oscnumubar)
      if ( m == conv::kSKmOscNuMuBarCCQE && show_sk_spectral_info_etrue_norm_binning) {
	temp = component + " ($E_{t} <$ 1.5 GeV)";
	streamer << setw(ncol1)  << setfill(' ') << temp
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nskoscnumubarccqe_etrue_na
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nskoscnumubarccqe_etrue_0a
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nskoscnumubarccqe_etrue_pa
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nskoscnumubarccqe_etrue_pa << pctsymbol 
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nskoscnumubarccqe_etrue_ma 
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nskoscnumubarccqe_etrue_ma << pctsymbol 
		 << endrow
		 << endl;
	temp = component + " (1.5 $\\leq E_{t} <$ 3.5 GeV)";
	streamer << setw(ncol1)  << setfill(' ') << temp
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nskoscnumubarccqe_etrue_nb
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nskoscnumubarccqe_etrue_0b
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nskoscnumubarccqe_etrue_pb
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nskoscnumubarccqe_etrue_pb << pctsymbol 
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nskoscnumubarccqe_etrue_mb 
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nskoscnumubarccqe_etrue_mb << pctsymbol 
		 << endrow
		 << endl;
	temp = component + " ($E_{t} \\geq$ 3.5 GeV)";
	streamer << setw(ncol1)  << setfill(' ') << temp
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nskoscnumubarccqe_etrue_nc
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nskoscnumubarccqe_etrue_0c
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nskoscnumubarccqe_etrue_pc
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nskoscnumubarccqe_etrue_pc << pctsymbol 
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nskoscnumubarccqe_etrue_mc 
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nskoscnumubarccqe_etrue_mc << pctsymbol 
		 << endrow
		 << endl;
	if(use_tex) streamer << "\\hline" << endl;
      }

      //spectral ccqe etrue (oscnumu)
      if ( m == conv::kSKmOscNuMuCCQE && show_sk_spectral_info_etrue_norm_binning) {
	temp = component + " ($E_{t} <$ 1.5 GeV)";
	streamer << setw(ncol1)  << setfill(' ') << temp
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nskoscnumuccqe_etrue_na
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nskoscnumuccqe_etrue_0a
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nskoscnumuccqe_etrue_pa
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nskoscnumuccqe_etrue_pa << pctsymbol 
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nskoscnumuccqe_etrue_ma 
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nskoscnumuccqe_etrue_ma << pctsymbol 
		 << endrow
		 << endl;
	temp = component + " (1.5 $\\leq E_{t} <$ 3.5 GeV)";
	streamer << setw(ncol1)  << setfill(' ') << temp
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nskoscnumuccqe_etrue_nb
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nskoscnumuccqe_etrue_0b
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nskoscnumuccqe_etrue_pb
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nskoscnumuccqe_etrue_pb << pctsymbol 
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nskoscnumuccqe_etrue_mb 
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nskoscnumuccqe_etrue_mb << pctsymbol 
		 << endrow
		 << endl;
	temp = component + " ($E_{t} \\geq$ 3.5 GeV)";
	streamer << setw(ncol1)  << setfill(' ') << temp
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nskoscnumuccqe_etrue_nc
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nskoscnumuccqe_etrue_0c
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nskoscnumuccqe_etrue_pc
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nskoscnumuccqe_etrue_pc << pctsymbol 
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nskoscnumuccqe_etrue_mc 
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nskoscnumuccqe_etrue_mc << pctsymbol 
		 << endrow
		 << endl;
	if(use_tex) streamer << "\\hline" << endl;
      }
      //spectral ccres etrue (numu)
      if ( m == conv::kSKmNuMuCC1pi && show_sk_spectral_info_etrue_norm_binning) {
	if(use_tex) streamer << "\\hline" << endl;
	temp = component + " ($E_{t} <$ 2.5 GeV)";
	streamer << setw(ncol1)  << setfill(' ') << temp
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumuccres_etrue_na
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumuccres_etrue_0a
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumuccres_etrue_pa
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknumuccres_etrue_pa << pctsymbol 
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumuccres_etrue_ma 
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknumuccres_etrue_ma << pctsymbol 
		 << endrow
		 << endl;
	temp = component + " ($E_{t} \\geq$ 2.5 GeV)";
	streamer << setw(ncol1)  << setfill(' ') << temp
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumuccres_etrue_nb
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumuccres_etrue_0b
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumuccres_etrue_pb
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknumuccres_etrue_pb << pctsymbol 
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumuccres_etrue_mb 
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknumuccres_etrue_mb << pctsymbol 
		 << endrow
		 << endl;
	if(use_tex) streamer << "\\hline" << endl;
      }
      //spectral ccres etrue (numubar)
      if ( m == conv::kSKmNuMuBarCC1pi && show_sk_spectral_info_etrue_norm_binning) {
	if(use_tex) streamer << "\\hline" << endl;
	temp = component + " ($E_{t} <$ 2.5 GeV)";
	streamer << setw(ncol1)  << setfill(' ') << temp
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumubarccres_etrue_na
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumubarccres_etrue_0a
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumubarccres_etrue_pa
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknumubarccres_etrue_pa << pctsymbol 
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumubarccres_etrue_ma 
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknumubarccres_etrue_ma << pctsymbol 
		 << endrow
		 << endl;
	temp = component + " ($E_{t} \\geq$ 2.5 GeV)";
	streamer << setw(ncol1)  << setfill(' ') << temp
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumubarccres_etrue_nb
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumubarccres_etrue_0b
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumubarccres_etrue_pb
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknumubarccres_etrue_pb << pctsymbol 
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumubarccres_etrue_mb 
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknumubarccres_etrue_mb << pctsymbol 
		 << endrow
		 << endl;
	if(use_tex) streamer << "\\hline" << endl;
      }
      //spectral ccres etrue (nue)
      if ( m == conv::kSKmNuECC1pi && show_sk_spectral_info_etrue_norm_binning) {
	if(use_tex) streamer << "\\hline" << endl;
	temp = component + " ($E_{t} <$ 2.5 GeV)";
	streamer << setw(ncol1)  << setfill(' ') << temp
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknueccres_etrue_na
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknueccres_etrue_0a
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknueccres_etrue_pa
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknueccres_etrue_pa << pctsymbol 
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknueccres_etrue_ma 
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknueccres_etrue_ma << pctsymbol 
		 << endrow
		 << endl;
	temp = component + " ($E_{t} \\geq$ 2.5 GeV)";
	streamer << setw(ncol1)  << setfill(' ') << temp
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknueccres_etrue_nb
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknueccres_etrue_0b
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknueccres_etrue_pb
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknueccres_etrue_pb << pctsymbol 
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknueccres_etrue_mb 
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknueccres_etrue_mb << pctsymbol 
		 << endrow
		 << endl;
	temp = component + " ($E_{t} \\geq$ 2.5 GeV)";
	if(use_tex) streamer << "\\hline" << endl;
      }
      //spectral ccres etrue (nuebar)
      if ( m == conv::kSKmNuEBarCC1pi && show_sk_spectral_info_etrue_norm_binning) {
	if(use_tex) streamer << "\\hline" << endl;
	temp = component + " ($E_{t} <$ 2.5 GeV)";
	streamer << setw(ncol1)  << setfill(' ') << temp
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknuebarccres_etrue_na
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknuebarccres_etrue_0a
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknuebarccres_etrue_pa
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknuebarccres_etrue_pa << pctsymbol 
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknuebarccres_etrue_ma 
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknuebarccres_etrue_ma << pctsymbol 
		 << endrow
		 << endl;
	temp = component + " ($E_{t} \\geq$ 2.5 GeV)";
	streamer << setw(ncol1)  << setfill(' ') << temp
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknuebarccres_etrue_nb
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknuebarccres_etrue_0b
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknuebarccres_etrue_pb
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknuebarccres_etrue_pb << pctsymbol 
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknuebarccres_etrue_mb 
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknuebarccres_etrue_mb << pctsymbol 
		 << endrow
		 << endl;
	if(use_tex) streamer << "\\hline" << endl;
      }
      //spectral ccres etrue (sig nue)
      if ( m == conv::kSKmOscNuECC1pi && show_sk_spectral_info_etrue_norm_binning) {
	if(use_tex) streamer << "\\hline" << endl;
	temp = component + " ($E_{t} <$ 2.5 GeV)";
	streamer << setw(ncol1)  << setfill(' ') << temp
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsksignueccres_etrue_na
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsksignueccres_etrue_0a
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsksignueccres_etrue_pa
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsksignueccres_etrue_pa << pctsymbol 
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsksignueccres_etrue_ma 
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsksignueccres_etrue_ma << pctsymbol 
		 << endrow
		 << endl;
	temp = component + " ($E_{t} \\geq$ 2.5 GeV)";
	streamer << setw(ncol1)  << setfill(' ') << temp
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsksignueccres_etrue_nb
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsksignueccres_etrue_0b
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsksignueccres_etrue_pb
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsksignueccres_etrue_pb << pctsymbol 
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsksignueccres_etrue_mb 
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsksignueccres_etrue_mb << pctsymbol 
		 << endrow
		 << endl;
	if(use_tex) streamer << "\\hline" << endl;
      }

//spectral ccres etrue (sig nuebar)
      if ( m == conv::kSKmOscNuEBarCC1pi && show_sk_spectral_info_etrue_norm_binning) {
	if(use_tex) streamer << "\\hline" << endl;
	temp = component + " ($E_{t} <$ 2.5 GeV)";
	streamer << setw(ncol1)  << setfill(' ') << temp
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsksignuebarccres_etrue_na
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsksignuebarccres_etrue_0a
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsksignuebarccres_etrue_pa
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsksignuebarccres_etrue_pa << pctsymbol 
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsksignuebarccres_etrue_ma 
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsksignuebarccres_etrue_ma << pctsymbol 
		 << endrow
		 << endl;
	temp = component + " ($E_{t} \\geq$ 2.5 GeV)";
	streamer << setw(ncol1)  << setfill(' ') << temp
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsksignuebarccres_etrue_nb
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsksignuebarccres_etrue_0b
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsksignuebarccres_etrue_pb
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsksignuebarccres_etrue_pb << pctsymbol 
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsksignuebarccres_etrue_mb 
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsksignuebarccres_etrue_mb << pctsymbol 
		 << endrow
		 << endl;
	if(use_tex) streamer << "\\hline" << endl;
      }
      //spectral ccres etrue (oscnumu)
      if ( m == conv::kSKmOscNuMuCC1pi && show_sk_spectral_info_etrue_norm_binning) {
	if(use_tex) streamer << "\\hline" << endl;
	temp = component + " ($E_{t} <$ 2.5 GeV)";
	streamer << setw(ncol1)  << setfill(' ') << temp
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nskoscnumuccres_etrue_na
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nskoscnumuccres_etrue_0a
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nskoscnumuccres_etrue_pa
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nskoscnumuccres_etrue_pa << pctsymbol 
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nskoscnumuccres_etrue_ma 
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nskoscnumuccres_etrue_ma << pctsymbol 
		 << endrow
		 << endl;
	temp = component + " ($E_{t} \\geq$ 2.5 GeV)";
	streamer << setw(ncol1)  << setfill(' ') << temp
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nskoscnumuccres_etrue_nb
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nskoscnumuccres_etrue_0b
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nskoscnumuccres_etrue_pb
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nskoscnumuccres_etrue_pb << pctsymbol 
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nskoscnumuccres_etrue_mb 
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nskoscnumuccres_etrue_mb << pctsymbol 
		 << endrow
		 << endl;
	if(use_tex) streamer << "\\hline" << endl;
      }
      //spectral ccres etrue (oscnumubar)
      if ( m == conv::kSKmOscNuMuBarCC1pi && show_sk_spectral_info_etrue_norm_binning) {
	if(use_tex) streamer << "\\hline" << endl;
	temp = component + " ($E_{t} <$ 2.5 GeV)";
	streamer << setw(ncol1)  << setfill(' ') << temp
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nskoscnumubarccres_etrue_na
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nskoscnumubarccres_etrue_0a
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nskoscnumubarccres_etrue_pa
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nskoscnumubarccres_etrue_pa << pctsymbol 
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nskoscnumubarccres_etrue_ma 
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nskoscnumubarccres_etrue_ma << pctsymbol 
		 << endrow
		 << endl;
	temp = component + " ($E_{t} \\geq$ 2.5 GeV)";
	streamer << setw(ncol1)  << setfill(' ') << temp
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nskoscnumubarccres_etrue_nb
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nskoscnumubarccres_etrue_0b
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nskoscnumubarccres_etrue_pb
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nskoscnumubarccres_etrue_pb << pctsymbol 
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nskoscnumubarccres_etrue_mb 
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nskoscnumubarccres_etrue_mb << pctsymbol 
		 << endrow
		 << endl;
	if(use_tex) streamer << "\\hline" << endl;
      }
      ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // spectral ereco ccqe (numu and numubar only)
      ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //spectral ccqe ereco (numu)
      if ( m == conv::kSKmNuMuCCQE && show_sk_ccqe_spectral_info_ereco) {
	temp = component + " ($E_{r} <$ 0.4 GeV)";
	streamer << setw(ncol1)  << setfill(' ') << temp
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumuccqe_na
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumuccqe_0a
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumuccqe_pa
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknumuccqe_pa << pctsymbol 
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumuccqe_ma 
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknumuccqe_ma << pctsymbol 
		 << endrow
		 << endl;
	temp = component + " (0.4 $\\leq E_{r} <$ 1.1 GeV)";
	streamer << setw(ncol1)  << setfill(' ') << temp
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumuccqe_nb
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumuccqe_0b
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumuccqe_pb
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknumuccqe_pb << pctsymbol 
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumuccqe_mb 
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknumuccqe_mb << pctsymbol 
		 << endrow
		 << endl;
	temp = component + " ($E_{r} \\geq$ 1.1 GeV)";
	streamer << setw(ncol1)  << setfill(' ') << temp
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumuccqe_nc
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumuccqe_0c
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumuccqe_pc
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknumuccqe_pc << pctsymbol 
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumuccqe_mc 
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknumuccqe_mc << pctsymbol 
		 << endrow
		 << endl;
	if(use_tex) streamer << "\\hline" << endl;
      }
      //ccqe spectral ereco (numubar)
      if ( m == conv::kSKmNuMuBarCCQE && show_sk_ccqe_spectral_info_ereco) {
	temp = component + " ($E_{r} <$ 0.4 GeV)";
	streamer << setw(ncol1)  << setfill(' ') << temp
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumubarccqe_na
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumubarccqe_0a
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumubarccqe_pa
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknumubarccqe_pa << pctsymbol 
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumubarccqe_ma 
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknumubarccqe_ma << pctsymbol 
		 << endrow
		 << endl;
	temp = component + " (0.4 $\\leq E_{r} <$ 1.1 GeV)";
	streamer << setw(ncol1)  << setfill(' ') << temp
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumubarccqe_nb
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumubarccqe_0b
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumubarccqe_pb
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknumubarccqe_pb << pctsymbol 
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumubarccqe_mb 
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknumubarccqe_mb << pctsymbol 
		 << endrow
		 << endl;
	temp = component + " ($E_{r} \\geq$ 1.1 GeV)";
	streamer << setw(ncol1)  << setfill(' ') << temp
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumubarccqe_nc
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumubarccqe_0c
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumubarccqe_pc
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknumubarccqe_pc << pctsymbol 
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumubarccqe_mc 
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknumubarccqe_mc << pctsymbol 
		 << endrow
		 << endl;
	if(use_tex) streamer << "\\hline" << endl;
      }
      ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // spectral ereco ccqe (nue binning)
      ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //spectral ccqe ereco (numu+numubar+oscnumu+oscnumubar)
      if ( m == conv::kSKmOscNuMuBarCCoth && show_sk_ccqe_spectral_info_ereco) {
	temp = "$N_{SK}$ (";
	temp += util::AsStringTex(s);
	temp += "): $\\nu_\\mu+\\bar{\\nu_\\mu}$";
	temp += " ($E_{r} <$ 0.35 GeV)";
	streamer << setw(ncol1)  << setfill(' ') << temp
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumu_ereco_na
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumu_ereco_0a
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumu_ereco_pa
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknumu_ereco_pa << pctsymbol 
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumu_ereco_ma 
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknumu_ereco_ma << pctsymbol 
		 << endrow
		 << endl;
	temp = "$N_{SK}$ (";
	temp += util::AsStringTex(s);
	temp += "): $\\nu_\\mu+\\bar{\\nu_\\mu}$";
	temp += " (0.35 $\\leq E_{r} <$ 0.8 GeV)";
	streamer << setw(ncol1)  << setfill(' ') << temp
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumu_ereco_nb
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumu_ereco_0b
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumu_ereco_pb
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknumu_ereco_pb << pctsymbol 
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumu_ereco_mb 
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknumu_ereco_mb << pctsymbol 
		 << endrow
		 << endl;
	temp = "$N_{SK}$ (";
	temp += util::AsStringTex(s);
	temp += "): $\\nu_\\mu+\\bar{\\nu_\\mu}$";
	temp += " ($E_{r} \\geq$ 0.8 GeV)";
	streamer << setw(ncol1)  << setfill(' ') << temp
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumu_ereco_nc
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumu_ereco_0c
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumu_ereco_pc
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknumu_ereco_pc << pctsymbol 
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumu_ereco_mc 
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknumu_ereco_mc << pctsymbol 
		 << endrow
		 << endl;
	if(use_tex) streamer << "\\hline" << endl;
      }
      //spectral ccqe ereco (nue+nuebar)
      if ( m == conv::kSKmNuEBarCCoth && show_sk_ccqe_spectral_info_ereco) {
	temp = "$N_{SK}$ (";
	temp += util::AsStringTex(s);
	temp += "): $\\nu_e+\\bar{\\nu_e}$";
	temp += " ($E_{r} <$ 0.35 GeV)";
	streamer << setw(ncol1)  << setfill(' ') << temp
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknue_ereco_na
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknue_ereco_0a
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknue_ereco_pa
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknue_ereco_pa << pctsymbol 
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknue_ereco_ma 
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknue_ereco_ma << pctsymbol 
		 << endrow
		 << endl;
	temp = "$N_{SK}$ (";
	temp += util::AsStringTex(s);
	temp += "): $\\nu_e+\\bar{\\nu_e}$";
	temp += " (0.35 $\\leq E_{r} <$ 0.8 GeV)";
	streamer << setw(ncol1)  << setfill(' ') << temp
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknue_ereco_nb
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknue_ereco_0b
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknue_ereco_pb
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknue_ereco_pb << pctsymbol 
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknue_ereco_mb 
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknue_ereco_mb << pctsymbol 
		 << endrow
		 << endl;
	temp = "$N_{SK}$ (";
	temp += util::AsStringTex(s);
	temp += "): $\\nu_e+\\bar{\\nu_e}$";
	temp += " ($E_{r} \\geq$ 0.8 GeV)";
	streamer << setw(ncol1)  << setfill(' ') << temp
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknue_ereco_nc
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknue_ereco_0c
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknue_ereco_pc
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknue_ereco_pc << pctsymbol 
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknue_ereco_mc 
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknue_ereco_mc << pctsymbol 
		 << endrow
		 << endl;
	if(use_tex) streamer << "\\hline" << endl;
      }
      //spectral ccqe ereco (oscnue+oscnuebar)
      if ( m == conv::kSKmOscNuEBarCCoth && show_sk_ccqe_spectral_info_ereco) {
	temp = "$N_{SK}$ (";
	temp += util::AsStringTex(s);
	temp += "): Oscillated $\\nu_e+\\bar{\\nu_e}$";
	temp += " ($E_{r} <$ 0.35 GeV)";
	streamer << setw(ncol1)  << setfill(' ') << temp
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nskoscnue_ereco_na
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nskoscnue_ereco_0a
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nskoscnue_ereco_pa
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nskoscnue_ereco_pa << pctsymbol 
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nskoscnue_ereco_ma 
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nskoscnue_ereco_ma << pctsymbol 
		 << endrow
		 << endl;
	temp = "$N_{SK}$ (";
	temp += util::AsStringTex(s);
	temp += "): Oscillated $\\nu_e+\\bar{\\nu_e}$";
	temp += " (0.35 $\\leq E_{r} <$ 0.8 GeV)";
	streamer << setw(ncol1)  << setfill(' ') << temp
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nskoscnue_ereco_nb
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nskoscnue_ereco_0b
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nskoscnue_ereco_pb
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nskoscnue_ereco_pb << pctsymbol 
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nskoscnue_ereco_mb 
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nskoscnue_ereco_mb << pctsymbol 
		 << endrow
		 << endl;
	temp = "$N_{SK}$ (";
	temp += util::AsStringTex(s);
	temp += "): Oscillated $\\nu_e+\\bar{\\nu_e}$";
	temp += " ($E_{r} \\geq$ 0.8 GeV)";
	streamer << setw(ncol1)  << setfill(' ') << temp
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nskoscnue_ereco_nc
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nskoscnue_ereco_0c
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nskoscnue_ereco_pc
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nskoscnue_ereco_pc << pctsymbol 
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nskoscnue_ereco_mc 
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nskoscnue_ereco_mc << pctsymbol 
		 << endrow
		 << endl;
	if(use_tex) streamer << "\\hline" << endl;
      }
      //spectral ccqe ereco (nc)
      if ((m == conv::kSKmNuENCoth_e || m == conv::kSKmNuENCoth_mu) && show_sk_ccqe_spectral_info_ereco) {
	temp = "$N_{SK}$ (";
	temp += util::AsStringTex(s);
	temp += "): NC";
	temp += " ($E_{r} <$ 0.35 GeV)";
	streamer << setw(ncol1)  << setfill(' ') << temp
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknc_ereco_na
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknc_ereco_0a
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknc_ereco_pa
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknc_ereco_pa << pctsymbol 
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknc_ereco_ma 
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknc_ereco_ma << pctsymbol 
		 << endrow
		 << endl;
	temp = "$N_{SK}$ (";
	temp += util::AsStringTex(s);
	temp += "): NC";
	temp += " (0.35 $\\leq E_{r} <$ 0.8 GeV)";
	streamer << setw(ncol1)  << setfill(' ') << temp
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknc_ereco_nb
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknc_ereco_0b
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknc_ereco_pb
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknc_ereco_pb << pctsymbol 
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknc_ereco_mb 
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknc_ereco_mb << pctsymbol 
		 << endrow
		 << endl;
	temp = "$N_{SK}$ (";
	temp += util::AsStringTex(s);
	temp += "): NC";
	temp += " ($E_{r} \\geq$ 0.8 GeV)";
	streamer << setw(ncol1)  << setfill(' ') << temp
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknc_ereco_nc
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknc_ereco_0c
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknc_ereco_pc
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknc_ereco_pc << pctsymbol 
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknc_ereco_mc 
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknc_ereco_mc << pctsymbol 
		 << endrow
		 << endl;
	if(use_tex) streamer << "\\hline" << endl;
      }
      ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // regular mode fill
      ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //a normal individual mode table fill
      streamer << setw(ncol1)  << setfill(' ') << component
	       << endcolumn
	       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nskmodes_n[m]
	       << endcolumn
	       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nskmodes_0[m]
	       << endcolumn
	       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nskmodes_p[m]
	       << endcolumn
	       << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nskmodes_p[m] << pctsymbol 
	       << endcolumn
	       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nskmodes_m[m] 
	       << endcolumn
	       << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nskmodes_m[m] << pctsymbol 
	       << endrow
	       << endl;
      ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // total of indvidual species for a given mode
      ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //cc total (numu)
      if ((util::Is1RMuSample(s) && m == conv::kSKmNuMuCCoth) || (util::Is1RESample(s) && m == conv::kSKmNuMuCCoth)) {
	temp = "$N_{SK}$ (" + util::AsStringTex(s) + "): $\\nu_{\\mu}$ CC all";
	streamer << setw(ncol1)  << setfill(' ') << temp
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumucc_n
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumucc_0
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumucc_p
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknumucc_p << pctsymbol 
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumucc_m 
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknumucc_m << pctsymbol 
		 << endrow
		 << endl;
	if(use_tex) streamer << "\\hline" << endl;
      }
      //nc total (numu)
      if ((util::Is1RMuSample(s) && m == conv::kSKmNuMuNCoth_mu) || (util::Is1RESample(s) && m == conv::kSKmNuMuNCoth_e)) {
	temp = "$N_{SK}$ (" + util::AsStringTex(s) + "): $\\nu_{\\mu}$ NC all";
	streamer << setw(ncol1)  << setfill(' ') << temp
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumunc_n
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumunc_0
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumunc_p
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknumunc_p << pctsymbol 
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumunc_m 
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknumunc_m << pctsymbol 
		 << endrow
		 << endl;
	if(use_tex) streamer << "\\hline" << endl;
      }
      //cc total (osc numu)
      if ((util::Is1RMuSample(s) && m == conv::kSKmOscNuMuCCoth) || (util::Is1RESample(s) && m == conv::kSKmOscNuMuCCoth)) {
	temp = "$N_{SK}$ (" + util::AsStringTex(s) + "): oscillation $\\nu_\\mu$ CC all";
	streamer << setw(ncol1)  << setfill(' ') << temp
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nskoscnumu_n
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nskoscnumu_0
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nskoscnumu_p
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nskoscnumu_p << pctsymbol 
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nskoscnumu_m 
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nskoscnumu_m << pctsymbol 
		 << endrow
		 << endl;
	if(use_tex) streamer << "\\hline" << endl;
      }
      //cc total (osc numubar)
      if ((util::Is1RMuSample(s) && m == conv::kSKmOscNuMuBarCCoth) || (util::Is1RESample(s) && m == conv::kSKmOscNuMuBarCCoth)) {
	temp = "$N_{SK}$ (" + util::AsStringTex(s) + "): oscillation $\\bar{\\nu_{\\mu}}$ CC all";
	streamer << setw(ncol1)  << setfill(' ') << temp
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nskoscnumubar_n
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nskoscnumubar_0
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nskoscnumubar_p
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nskoscnumubar_p << pctsymbol 
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nskoscnumubar_m 
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nskoscnumubar_m << pctsymbol 
		 << endrow
		 << endl;
	if(use_tex) streamer << "\\hline" << endl;
      }
      //cc total (numubar)
      if ((util::Is1RMuSample(s) && m == conv::kSKmNuMuBarCCoth) || (util::Is1RESample(s) && m == conv::kSKmNuMuBarCCoth)) {
	temp = "$N_{SK}$ (" + util::AsStringTex(s) + "): $\\bar{\\nu_{\\mu}}$ CC all";
	streamer << setw(ncol1)  << setfill(' ') << temp
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumubarcc_n
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumubarcc_0
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumubarcc_p
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknumubarcc_p << pctsymbol 
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumubarcc_m 
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknumubarcc_m << pctsymbol 
		 << endrow
		 << endl;
	if(use_tex) streamer << "\\hline" << endl;
      }
      //nc total (numubar)
      if ((util::Is1RMuSample(s) && m == conv::kSKmNuMuBarNCoth_mu) || (util::Is1RESample(s) && m == conv::kSKmNuMuBarNCoth_e)) {
	temp = "$N_{SK}$ (" + util::AsStringTex(s) + "): $\\bar{\\nu_{\\mu}}$ NC all";
	streamer << setw(ncol1)  << setfill(' ') << temp
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumubarnc_n
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumubarnc_0
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumubarnc_p
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknumubarnc_p << pctsymbol 
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumubarnc_m 
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknumubarnc_m << pctsymbol 
		 << endrow
		 << endl;
	if(use_tex) streamer << "\\hline" << endl;
      }
      //cc total (nue)
      if ((util::Is1RMuSample(s) && m == conv::kSKmNuECCoth) || (util::Is1RESample(s) && m == conv::kSKmNuECCoth)) {
	temp = "$N_{SK}$ (" + util::AsStringTex(s) + "): $\\nu_{e}$ CC all";
	streamer << setw(ncol1)  << setfill(' ') << temp
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknuecc_n
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknuecc_0
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknuecc_p
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknuecc_p << pctsymbol 
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknuecc_m 
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknuecc_m << pctsymbol 
		 << endrow
		 << endl;
	if(use_tex) streamer << "\\hline" << endl;
      }
      //nc total (nue)
      if ((util::Is1RMuSample(s) && m == conv::kSKmNuENCoth_mu) || (util::Is1RESample(s) && m == conv::kSKmNuENCoth_e)) {
	temp = "$N_{SK}$ (" + util::AsStringTex(s) + "): $\\nu_{e}$ NC all";
	streamer << setw(ncol1)  << setfill(' ') << temp
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknuenc_n
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknuenc_0
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknuenc_p
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknuenc_p << pctsymbol 
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknuenc_m 
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknuenc_m << pctsymbol 
		 << endrow
		 << endl;
	if(use_tex) streamer << "\\hline" << endl;
      }
      //cc total (nuebar)
      if ((util::Is1RMuSample(s) && m == conv::kSKmNuEBarCCoth) || (util::Is1RESample(s) && m == conv::kSKmNuEBarCCoth)) {
	temp = "$N_{SK}$ (" + util::AsStringTex(s) + "): $\\bar{\\nu_{e}}$ CC all";
	streamer << setw(ncol1)  << setfill(' ') << temp
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknuebarcc_n
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknuebarcc_0
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknuebarcc_p
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknuebarcc_p << pctsymbol 
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknuebarcc_m 
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknuebarcc_m << pctsymbol 
		 << endrow
		 << endl;
	if(use_tex) streamer << "\\hline" << endl;
      }
      //nc total (nuebar)
      if ((util::Is1RMuSample(s) && m == conv::kSKmNuEBarNCoth_mu) || (util::Is1RESample(s) && m == conv::kSKmNuEBarNCoth_e)) {
	temp = "$N_{SK}$ (" + util::AsStringTex(s) + "): $\\bar{\\nu_{e}}$ NC all";
	streamer << setw(ncol1)  << setfill(' ') << temp
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknuebarnc_n
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknuebarnc_0
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknuebarnc_p
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknuebarnc_p << pctsymbol 
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknuebarnc_m 
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknuebarnc_m << pctsymbol 
		 << endrow
		 << endl;
	if(use_tex) streamer << "\\hline" << endl;
      }
      //cc total (sig nue)
      if ((util::Is1RMuSample(s) && m == conv::kSKmOscNuECCoth) || (util::Is1RESample(s) && m == conv::kSKmOscNuECCoth)) {
	temp = "$N_{SK}$ (" + util::AsStringTex(s) + "): oscillation $\\nu_e$ CC all";
	streamer << setw(ncol1)  << setfill(' ') << temp
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknuesig_n
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknuesig_0
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknuesig_p
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknuesig_p << pctsymbol 
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknuesig_m 
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknuesig_m << pctsymbol 
		 << endrow
		 << endl;
	if(use_tex) streamer << "\\hline" << endl;
      }

      //cc total (sig nuebar)
      if ((util::Is1RMuSample(s) && m == conv::kSKmOscNuEBarCCoth) || (util::Is1RESample(s) && m == conv::kSKmOscNuEBarCCoth)) {
	temp = "$N_{SK}$ (" + util::AsStringTex(s) + "): oscillation $\\bar{\\nu_{e}}$ CC all";
	streamer << setw(ncol1)  << setfill(' ') << temp
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknuebarsig_n
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknuebarsig_0
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknuebarsig_p
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknuebarsig_p << pctsymbol 
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknuebarsig_m 
		 << endcolumn
		 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknuebarsig_m << pctsymbol 
		 << endrow
		 << endl;
	if(use_tex) streamer << "\\hline" << endl;
      }
      

      ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // spectral total for each species etrue
      ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //spectral total etrue (numu)
      if (((util::Is1RMuSample(s) && m == conv::kSKmNuMuNCoth_mu) || (util::Is1RESample(s) && m == conv::kSKmNuMuNCoth_e)) && show_sk_spectral_info_etrue_flux_binning_mu) {
	char bufl[100], bufu[100];
	for(int itb = 0; itb < NEtrueBinNuMu; itb++) {
	  snprintf(bufl, sizeof(bufl), "%g",  skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuCCQE, conv::kSKMCtFin)->GetYaxis()->GetBinLowEdge(EtrueBinNuMu[itb]+1));
	  snprintf(bufu, sizeof(bufu), "%g",  skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuCCQE, conv::kSKMCtFin)->GetYaxis()->GetBinUpEdge(EtrueBinNuMu[itb+1]));
	  temp = "$N_{SK}$ (" + util::AsStringTex(s) + "): $\\nu_{\\mu}$ all (" + bufl + " $\\leq E_{t} <$ " + bufu + "GeV)";
	  streamer << setw(ncol1)  << setfill(' ') << temp
		   << endcolumn
		   << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumu_etrue_n[itb]
		   << endcolumn
		   << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumu_etrue_0[itb]
		   << endcolumn
		   << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumu_etrue_p[itb]
		   << endcolumn
		   << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknumu_etrue_p[itb] << pctsymbol 
		   << endcolumn
		   << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumu_etrue_m[itb]
		   << endcolumn
		   << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknumu_etrue_m[itb] << pctsymbol 
		   << endrow
		   << endl;
	}//loop over etrue bins
	if(use_tex) streamer << "\\hline" << endl;
      }
      //spectral total etrue (numu)
      if (((util::Is1RMuSample(s) && m == conv::kSKmNuMuBarNCoth_mu) || (util::Is1RESample(s) && m == conv::kSKmNuMuBarNCoth_e)) && show_sk_spectral_info_etrue_flux_binning_mu) {
	char bufl[100], bufu[100];
	for(int itb = 0; itb < NEtrueBinNuMuBar; itb++) {
	  snprintf(bufl, sizeof(bufl), "%g",  skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuBarCCQE, conv::kSKMCtFin)->GetYaxis()->GetBinLowEdge(EtrueBinNuMuBar[itb]+1));
	  snprintf(bufu, sizeof(bufu), "%g",  skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuBarCCQE, conv::kSKMCtFin)->GetYaxis()->GetBinUpEdge(EtrueBinNuMuBar[itb+1]));
	  temp = "$N_{SK}$ (" + util::AsStringTex(s) + "): $\\bar{\\nu_{\\mu}}$ all (" + bufl + " $\\leq E_{t} <$ " + bufu + "GeV)";
	  streamer << setw(ncol1)  << setfill(' ') << temp
		   << endcolumn
		   << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumubar_etrue_n[itb]
		   << endcolumn
		   << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumubar_etrue_0[itb]
		   << endcolumn
		   << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumubar_etrue_p[itb]
		   << endcolumn
		   << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknumubar_etrue_p[itb] << pctsymbol 
		   << endcolumn
		   << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumubar_etrue_m[itb]
		   << endcolumn
		   << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknumubar_etrue_m[itb] << pctsymbol 
		   << endrow
		   << endl;
	}//loop over etrue bins
	if(use_tex) streamer << "\\hline" << endl;
      }
      //spectral total etrue (nue)
      if (((util::Is1RMuSample(s) && m == conv::kSKmNuENCoth_mu) || (util::Is1RESample(s) && m == conv::kSKmNuENCoth_e)) && show_sk_spectral_info_etrue_flux_binning) {
	char bufl[100], bufu[100];
	for(int itb = 0; itb < NEtrueBinNuE; itb++) {
	  snprintf(bufl, sizeof(bufl), "%g",  skexp->GetSuperkMC()->Template(s, conv::kSKmNuECCQE, conv::kSKMCtFin)->GetYaxis()->GetBinLowEdge(EtrueBinNuE[itb]+1));
	  snprintf(bufu, sizeof(bufu), "%g",  skexp->GetSuperkMC()->Template(s, conv::kSKmNuECCQE, conv::kSKMCtFin)->GetYaxis()->GetBinUpEdge(EtrueBinNuE[itb+1]));
	  temp = "$N_{SK}$ (" + util::AsStringTex(s) + "): $\\nu_{e}$ all (" + bufl + " $\\leq E_{t} <$ " + bufu + "GeV)";
	  streamer << setw(ncol1)  << setfill(' ') << temp
		   << endcolumn
		   << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknue_etrue_n[itb]
		   << endcolumn
		   << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknue_etrue_0[itb]
		   << endcolumn
		   << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknue_etrue_p[itb]
		   << endcolumn
		   << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknue_etrue_p[itb] << pctsymbol 
		   << endcolumn
		   << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknue_etrue_m[itb]
		   << endcolumn
		   << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknue_etrue_m[itb] << pctsymbol 
		   << endrow
		   << endl;
	}//loop over etrue bins
	if(use_tex) streamer << "\\hline" << endl;
      }
      //spectral total etrue (nuebar)
      if (((util::Is1RMuSample(s) && m == conv::kSKmNuEBarNCoth_mu) || (util::Is1RESample(s) && m == conv::kSKmNuEBarNCoth_e)) && show_sk_spectral_info_etrue_flux_binning) {
	char bufl[100], bufu[100];
	for(int itb = 0; itb < NEtrueBinNuEBar; itb++) {
	  snprintf(bufl, sizeof(bufl), "%g",  skexp->GetSuperkMC()->Template(s, conv::kSKmNuEBarCCQE, conv::kSKMCtFin)->GetYaxis()->GetBinLowEdge(EtrueBinNuEBar[itb]+1));
	  snprintf(bufu, sizeof(bufu), "%g",  skexp->GetSuperkMC()->Template(s, conv::kSKmNuEBarCCQE, conv::kSKMCtFin)->GetYaxis()->GetBinUpEdge(EtrueBinNuEBar[itb+1]));
	  temp = "$N_{SK}$ (" + util::AsStringTex(s) + "): $\\bar{\\nu_{e}}$ all (" + bufl + " $\\leq E_{t} <$ " + bufu + "GeV)";
	  streamer << setw(ncol1)  << setfill(' ') << temp
		   << endcolumn
		   << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknuebar_etrue_n[itb]
		   << endcolumn
		   << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknuebar_etrue_0[itb]
		   << endcolumn
		   << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknuebar_etrue_p[itb]
		   << endcolumn
		   << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknuebar_etrue_p[itb] << pctsymbol 
		   << endcolumn
		   << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknuebar_etrue_m[itb]
		   << endcolumn
		   << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknuebar_etrue_m[itb] << pctsymbol 
		   << endrow
		   << endl;
	}//loop over etrue bins
	if(use_tex) streamer << "\\hline" << endl;
      }
      //spectral total etrue (oscnumu)
      if (((util::Is1RMuSample(s) && m == conv::kSKmOscNuMuCCoth) || (util::Is1RESample(s) && m == conv::kSKmOscNuMuCCoth)) && show_sk_spectral_info_etrue_flux_binning) {
	char bufl[100], bufu[100];
	for(int itb = 0; itb < NEtrueBinOscNuMu; itb++) {
	  snprintf(bufl, sizeof(bufl), "%g",  skexp->GetSuperkMC()->Template(s, conv::kSKmOscNuMuCCQE, conv::kSKMCtFin)->GetYaxis()->GetBinLowEdge(EtrueBinOscNuMu[itb]+1));
	  snprintf(bufu, sizeof(bufu), "%g",  skexp->GetSuperkMC()->Template(s, conv::kSKmOscNuMuCCQE, conv::kSKMCtFin)->GetYaxis()->GetBinUpEdge(EtrueBinOscNuMu[itb+1]));
	  temp = "$N_{SK}$ (" + util::AsStringTex(s) + "): oscillation $\\nu_{\\mu}$ all (" + bufl + " $\\leq E_{t} <$ " + bufu + "GeV)";
	  streamer << setw(ncol1)  << setfill(' ') << temp
		   << endcolumn
		   << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nskoscnumu_etrue_n[itb]
		   << endcolumn
		   << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nskoscnumu_etrue_0[itb]
		   << endcolumn
		   << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nskoscnumu_etrue_p[itb]
		   << endcolumn
		   << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nskoscnumu_etrue_p[itb] << pctsymbol 
		   << endcolumn
		   << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nskoscnumu_etrue_m[itb]
		   << endcolumn
		   << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nskoscnumu_etrue_m[itb] << pctsymbol 
		   << endrow
		   << endl;
	}//loop over etrue bins
	if(use_tex) streamer << "\\hline" << endl;
      }
      //spectral total etrue (oscnumubar)
      if (((util::Is1RMuSample(s) && m == conv::kSKmOscNuMuBarCCoth) || (util::Is1RESample(s) && m == conv::kSKmOscNuMuBarCCoth)) && show_sk_spectral_info_etrue_flux_binning) {
	char bufl[100], bufu[100];
	for(int itb = 0; itb < NEtrueBinOscNuMuBar; itb++) {
	  snprintf(bufl, sizeof(bufl), "%g",  skexp->GetSuperkMC()->Template(s, conv::kSKmOscNuMuBarCCQE, conv::kSKMCtFin)->GetYaxis()->GetBinLowEdge(EtrueBinOscNuMuBar[itb]+1));
	  snprintf(bufu, sizeof(bufu), "%g",  skexp->GetSuperkMC()->Template(s, conv::kSKmOscNuMuBarCCQE, conv::kSKMCtFin)->GetYaxis()->GetBinUpEdge(EtrueBinOscNuMuBar[itb+1]));
	  temp = "$N_{SK}$ (" + util::AsStringTex(s) + "): oscillation $\\bar{\\nu_{\\mu}}$ all (" + bufl + " $\\leq E_{t} <$ " + bufu + "GeV)";
	  streamer << setw(ncol1)  << setfill(' ') << temp
		   << endcolumn
		   << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nskoscnumubar_etrue_n[itb]
		   << endcolumn
		   << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nskoscnumubar_etrue_0[itb]
		   << endcolumn
		   << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nskoscnumubar_etrue_p[itb]
		   << endcolumn
		   << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nskoscnumubar_etrue_p[itb] << pctsymbol 
		   << endcolumn
		   << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nskoscnumubar_etrue_m[itb]
		   << endcolumn
		   << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nskoscnumubar_etrue_m[itb] << pctsymbol 
		   << endrow
		   << endl;
	}//loop over etrue bins
	if(use_tex) streamer << "\\hline" << endl;
      }
      //spectral total etrue (signue)
      if (((util::Is1RMuSample(s) && m == conv::kSKmOscNuECCoth) || (util::Is1RESample(s) && m == conv::kSKmOscNuECCoth)) && show_sk_spectral_info_etrue_flux_binning_mu) {
	char bufl[100], bufu[100];
	for(int itb = 0; itb < NEtrueBinOscNuE; itb++) {
	  snprintf(bufl, sizeof(bufl), "%g",  skexp->GetSuperkMC()->Template(s, conv::kSKmOscNuECCQE, conv::kSKMCtFin)->GetYaxis()->GetBinLowEdge(EtrueBinOscNuE[itb]+1));
	  snprintf(bufu, sizeof(bufu), "%g",  skexp->GetSuperkMC()->Template(s, conv::kSKmOscNuECCQE, conv::kSKMCtFin)->GetYaxis()->GetBinUpEdge(EtrueBinOscNuE[itb+1]));
	  temp = "$N_{SK}$ (" + util::AsStringTex(s) + "): oscillation $\\nu_{e}$ all (" + bufl + " $\\leq E_{t} <$ " + bufu + "GeV)";
	  streamer << setw(ncol1)  << setfill(' ') << temp
		   << endcolumn
		   << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsksignue_etrue_n[itb]
		   << endcolumn
		   << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsksignue_etrue_0[itb]
		   << endcolumn
		   << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsksignue_etrue_p[itb]
		   << endcolumn
		   << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsksignue_etrue_p[itb] << pctsymbol 
		   << endcolumn
		   << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsksignue_etrue_m[itb]
		   << endcolumn
		   << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsksignue_etrue_m[itb] << pctsymbol 
		   << endrow
		   << endl;
	}//loop over etrue bins
	if(use_tex) streamer << "\\hline" << endl;
      }
      //spectral total etrue (signuebar)
      if (((util::Is1RMuSample(s) && m == conv::kSKmOscNuEBarCCoth) || (util::Is1RESample(s) && m == conv::kSKmOscNuEBarCCoth)) && show_sk_spectral_info_etrue_flux_binning_mu) {
	char bufl[100], bufu[100];
	for(int itb = 0; itb < NEtrueBinOscNuEBar; itb++) {
	  snprintf(bufl, sizeof(bufl), "%g",  skexp->GetSuperkMC()->Template(s, conv::kSKmOscNuEBarCCQE, conv::kSKMCtFin)->GetYaxis()->GetBinLowEdge(EtrueBinOscNuEBar[itb]+1));
	  snprintf(bufu, sizeof(bufu), "%g",  skexp->GetSuperkMC()->Template(s, conv::kSKmOscNuEBarCCQE, conv::kSKMCtFin)->GetYaxis()->GetBinUpEdge(EtrueBinOscNuEBar[itb+1]));
	  temp = "$N_{SK}$ (" + util::AsStringTex(s) + "): oscillation $\\bar{\\nu_{e}}$ all (" + bufl + " $\\leq E_{t} <$ " + bufu + "GeV)";
	  streamer << setw(ncol1)  << setfill(' ') << temp
		   << endcolumn
		   << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsksignuebar_etrue_n[itb]
		   << endcolumn
		   << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsksignuebar_etrue_0[itb]
		   << endcolumn
		   << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsksignuebar_etrue_p[itb]
		   << endcolumn
		   << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsksignuebar_etrue_p[itb] << pctsymbol 
		   << endcolumn
		   << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsksignuebar_etrue_m[itb]
		   << endcolumn
		   << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsksignuebar_etrue_m[itb] << pctsymbol 
		   << endrow
		   << endl;
	}//loop over etrue bins
	if(use_tex) streamer << "\\hline" << endl;
      }//signue etrue (flux binning)
      if (use_tex && util::Is1RMuSample(s))
	switch(m) {
	case conv::kSKmNuMuNCoth_mu:
	case conv::kSKmNuMuBarNCoth_mu:
	case conv::kSKmNuENCoth_mu:
	case conv::kSKmNuEBarNCoth_mu:
	  streamer << "\\hline" << endl;
	  break;
	}
      else if (use_tex && util::Is1RESample(s))
	switch(m) {
	case conv::kSKmNuMuNCoth_e:
	case conv::kSKmNuMuBarNCoth_e:
	case conv::kSKmNuENCoth_e:
	case conv::kSKmNuEBarNCoth_e:
	  streamer << "\\hline" << endl;
	  break;
	}
    }//end loop over sk modes
  
    if(!use_tex) {
      streamer << setw(150)  << setfill('-') << " " << endl;
    }
    else
      streamer << "\\hline" << endl;

    if(use_tex) {
      streamer << "\\end{tabular}" << endl;
      string hypothesisstr(hypothesis);
      string sys = util::AsStringPaperTex(util::PhysParamIdFromString(systematic));
      size_t found=hypothesisstr.find_first_of("_");
      while(found!=string::npos) {
	hypothesisstr[found]='-';
	found=hypothesisstr.find_first_of("_",found+1);
      }
      streamer << "\\caption{ Effect of " << sys << " for ";
      if(input_mixing12 < conv::kASmallNum && input_mixing13 < conv::kASmallNum 
	 && input_mixing23 < conv::kASmallNum && input_dm12sq < conv::kASmallNum
	 && input_dm23sq < conv::kASmallNum && input_dm12sq < conv::kASmallNum) {
	streamer << "no oscillations (all mixing parameters set to 0)";
      }
      else {
	streamer << "for 3 flavour oscillations,"
		 << " using oscillation parameters: " 
		 << " $\\sin^2\\theta_{12} = $"  << input_mixing12
		 << "; $\\sin^2\\theta_{13} = $" << input_mixing13
		 << "; $\\sin^2\\theta_{23} = $" << input_mixing23
	  //<< " $\\sin^22\\theta_{12} = $"  << input_mixing12
	  //<< "; $\\sin^22\\theta_{13} = $" << input_mixing13
	  //	 << "; $\\sin^22\\theta_{23} = $" << input_mixing23
		 << "; $\\delta_{CP} = $        " << input_deltacp
		 << "; $\\Delta m^2_{21} = $    " << input_dm12sq
		 << "; $\\Delta m^2_{32} = $    " << input_dm23sq;
      }
      streamer << ". }" << endl;
      string shortsamplestr = util::Is1RMuSample(s) ? "1RmuT_" : "1ReT_";
      string spectralstrreco = show_sk_ccqe_spectral_info_ereco ? "_spectralreco" : "";
      string spectralstrtrue = show_sk_ccqe_spectral_info_ereco ? "_spectraltrue" : "";
      streamer << "\\label{tab:num_" << systematic << "_" << shortsamplestr << hypothesis 
	       << "_t13:" << input_mixing13 << spectralstrreco << spectralstrtrue << "}" << endl;
      streamer << "\\end{table}" << endl;
      streamer << "\\clearpage" << endl;
    }

    cout << streamer.str() << endl;

    // write the table?
    if(output_file != "") {
      long size;

      // get size of file
      streamer.seekg(0,ifstream::end);
      size=streamer.tellg();
      streamer.seekg(0);

      ofstream outfile (output_file.c_str(),ios_base::app);

      // write to outfile
      outfile.write (streamer.str().c_str(),size);
  
      outfile.close();
    }

    streamer_summary << setw(ncol1)  << setfill(' ')
		     << util::AsStringPaperTex(util::PhysParamIdFromString(systematic))
		     << endcolumn
		     << setw(ncolnev)  << setfill(' ') << setprecision(2)
		     << TMath::Max(TMath::Abs(frac_nsk_p), TMath::Abs(frac_nsk_m))
		     << endl;

    cout << "Done sytematic! " << systematic << endl;

    //reset the systematic
    skexp->ResetSystematicVariations();
    ph.SetCurrValue(conv::kPhypMixing12,  input_mixing12);
    ph.SetCurrValue(conv::kPhypMixing13,  input_mixing13);
    ph.SetCurrValue(conv::kPhypMixing23,  input_mixing23);
    ph.SetCurrValue(conv::kPhypDeltaCP,   input_deltacp);
    ph.SetCurrValue(conv::kPhypDm12sq,    input_dm12sq);
    ph.SetCurrValue(conv::kPhypDm23sq,    input_dm23sq);
  }//isys

  // don't forget to write document footer to outfile
  ofstream outfile (output_file.c_str(),ios_base::trunc&&ios_base::out);
  string f;
  f += ("\\end{document}\n");
  outfile.write (f.c_str(),f.size());
  outfile.close();

  cout << streamer_summary.str() << endl;
  ofstream outfilesummary(Form("%s.summary",output_file.c_str()),ios_base::out);
  long size;
  streamer_summary.seekg(0,ifstream::end);
  size=streamer_summary.tellg();
  streamer_summary.seekg(0);
  outfilesummary.write(streamer_summary.str().c_str(),size);
  outfilesummary.close();
}
