//
// ------------------------------------------------------------------------------------------------
// Test systematic errors.
//
// The input (systematic) string can be:
//	- sk_energy_scale 
//	- sk_1Rmu_ccqe_other_eff
//	- sk_1Rmu_ccqe_ring_counting_eff : 
//             This energy-dependent correlated error can not be used in this test macro. Instead,
//             you can tweak specific components of the efficiency tweaking vector (neglecting,
//             correlation which, in the fitter, are handled using the Cholesky method). Use:
//	       - sk_1Rmu_ccqe_ring_counting_eff_0 
//	       - sk_1Rmu_ccqe_ring_counting_eff_1
//	       - sk_1Rmu_ccqe_ring_counting_eff_2
//	- sk_1Rmu_ccnonqe_eff
//	- sk_1Rmu_nc_eff
//	- sk_1Rmu_nuecc_eff
//      - banff
//             This correlated error can not be used in this test macro. Instead,
//             you can tweak specific components of the BANFF tweaking vector (neglecting,
//             correlation which, in the fitter, are handled using the Cholesky method). Use:
//             - banff_00
//             - ...
//             - banff_25
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//
// Major revisions by TD. 25/11/11.
//  Now takes SK sample as input
//  Can write output tables to tex file
//  Restructuring of SK modes part of code means its easier to change in the future:
//   - if order or number of kSKReactionModes[s][m] changes, just need to modify lines ~120-160
// Major revisions by TD. 15/03/12.
//  Takes advantage of util classes for determining which 'total' to put a mode in -> code greatly simplified
//  Can run with or without BANFF best fit values
//  Now prints Etrue bins (for numu CCQE) to see effect of BANFF flux bins

#include <iomanip.h>

void test_systematic_error_propagation
(
  const char * systematic,                                                  /* options: see above */
  const char * sksample                   = "1Rmu",                         /* options: "1Rmu", "1Re" */
  const char * hypothesis                 = "3f_osc",                       /* options: "2f_numu_disap", "3f_osc" */
  double       input_mixing12             = 0.86,                           /* */
  double       input_mixing13             = 0.10,                           /* */
  double       input_mixing23             = 1.00,                           /* */
  double       input_deltacp              = 0.00,                           /* */
  double       input_dm12sq               = 7.58E-5,                        /* */
  double       input_dm23sq               = 2.32E-3,                        /* */
  bool         use_banff_bestfit          = true,                           /* if false, BANFF::RestoreNominalMC() is called */
  bool         use_tex                    = false,                          /* print-out results as a TeX table (for inclusion in tech. note)*/
  bool         show_sk_ccqe_spectral_info_ereco = false,                    /* */
  bool         show_sk_ccqe_spectral_info_etrue = false,                    /* */
  const char * output_file                = "",                             /* */
  const char * config_file                = "analysis_configuration.xml"    /* */
)
{
  using namespace std;
  using namespace t2kana;

  // Analysis configuration
  const char * config  = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);
    
  // Global options
  Opt::Instance()->gOscMixingParamChoice = conv::kNuMixSinsq2Theta;
  Opt::Instance()->gMassHierarchyAssumed = conv::kNmhNormal;

  // Get requested sample
  conv::SuperkSample_t s = util::SKSampleFromString(sksample);
  assert(s != conv::kSKsUndefined);

  // Set physics hypothesis
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
  cout << "- SuperK sample: " << util::AsString(s) << endl;
  cout << "- Varied systematic parameter: " << systematic << endl;
  cout << "- Physics hypothesis: " << ph << endl;
  cout << *Opt::Instance() << endl;

  // Precison for printed-out numerical results
  int nprecisionnev = 4; // <-- for number of events
  int nprecisionpct = 2; // <-- for % changes

  // Create & config superk spectrum predictor
  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read(config);
  
  // Use BANFF best fit?
  if(!use_banff_bestfit)
    skexp->GetBANFF()->RestoreNominalMC();

  // Get the EReco bin corresponding to the edges of SK CCQE efficiency systematic
  Long64_t ErecoBin1NuMu    = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuMuCCQE    )->FindBin(0.4-conv::kASmallNum);
  Long64_t ErecoBin2NuMu    = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuMuCCQE    )->FindBin(1.1-conv::kASmallNum);
  Long64_t ErecoBin1NuMuBar = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuMuBarCCQE )->FindBin(0.4-conv::kASmallNum);
  Long64_t ErecoBin2NuMuBar = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuMuBarCCQE )->FindBin(1.1-conv::kASmallNum);

  // Get the ETrue bin edges for NuMuCCQE
  int NEtrueBinNuMu = conv::kXsNEtrueBins;
  Long64_t EtrueBinNuMu[conv::kXsNEtrueBins+1];
  for(int itb = 0; itb < conv::kXsNEtrueBins+1; itb++) {
    EtrueBinNuMu[itb] = -1;
    EtrueBinNuMu[itb] = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuCCQE, conv::kSKMCtFin)->GetYaxis()->FindBin(conv::kXsEtrueBinEdges[util::FindSampleSlot(s)][util::FindModeSlot(s,conv::kSKmNuMuCCQE)][itb]-conv::kASmallNum);
    if (EtrueBinNuMu[itb] > conv::kNEtrueBins) EtrueBinNuMu[itb] = conv::kNEtrueBins;
    if (EtrueBinNuMu[itb] < 1 && itb != 0) {NEtrueBinNuMu = itb-1; break;}
  }

  // Is it a BANFF systematic? (if it is, 1 sigma errors are not definied)
  string sys(systematic);
  bool sys_is_banff = sys.find("banff") >= 0;

  //
  // nominal
  //

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

  for (int im = 0; im < conv::knm; im++) {
    conv::SuperkMode_t m = conv::kSKReactionModes[s][im];
    if(m == conv::kSKmUndefined) continue;
    if(util::FindModeSlot(s,m) < 0) continue;
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
      else if(util::IsSignalNuE(m))
	nsknuesig_0 += nskmodes_0[m];
      else {
	cerr << "CC mode " << AsString(m) << " is not numu, numubar, nue, nuebar, nuesig!" << endl;
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
  double nsknumuccqe_0a        = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuMuCCQE    )->Integral(1,                  ErecoBin1NuMu);
  double nsknumuccqe_0b        = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuMuCCQE    )->Integral(ErecoBin1NuMu+1,    ErecoBin2NuMu);
  double nsknumuccqe_0c        = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuMuCCQE    )->Integral(ErecoBin2NuMu+1,    conv::kNErecoBins1Rmu);
  double nsknumubarccqe_0a     = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuMuBarCCQE )->Integral(1,                  ErecoBin1NuMuBar);
  double nsknumubarccqe_0b     = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuMuBarCCQE )->Integral(ErecoBin1NuMuBar+1, ErecoBin2NuMuBar);
  double nsknumubarccqe_0c     = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuMuBarCCQE )->Integral(ErecoBin2NuMuBar+1, conv::kNErecoBins1Rmu);

  // spectral numu ccqe (etrue)
  double nsknumuccqe_0t[conv::kXsNEtrueBins];
  for(int itb = 0; itb < NEtrueBinNuMu; itb++) {
    nsknumuccqe_0t[itb] = 0.;
    nsknumuccqe_0t[itb] = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBinNuMu[itb]+1, EtrueBinNuMu[itb+1]);
    cout << "bin edges are from " << skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuCCQE, conv::kSKMCtFin)->GetYaxis()->GetBinLowEdge(EtrueBinNuMu[itb]+1) << " to " << skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuCCQE, conv::kSKMCtFin)->GetYaxis()->GetBinUpEdge(EtrueBinNuMu[itb+1]) << endl;;
  }

  //
  // +1 sigma
  //
  
  if(sys_is_banff)
    skexp->TweakSystematic(systematic,+2.);
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

  for (int im = 0; im < conv::knm; im++) {
    conv::SuperkMode_t m = conv::kSKReactionModes[s][im];
    if(m == conv::kSKmUndefined) continue;
    if(util::FindModeSlot(s,m) < 0) continue;
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
      else if(util::IsSignalNuE(m))
	nsknuesig_p += nskmodes_p[m];
      else {
	cerr << "CC mode " << AsString(m) << " is not numu, numubar, nue, nuebar, nuesig!" << endl;
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
  double nsknumuccqe_pa        = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuMuCCQE    )->Integral(1,                  ErecoBin1NuMu);
  double nsknumuccqe_pb        = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuMuCCQE    )->Integral(ErecoBin1NuMu+1,    ErecoBin2NuMu);
  double nsknumuccqe_pc        = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuMuCCQE    )->Integral(ErecoBin2NuMu+1,    conv::kNErecoBins1Rmu);
  double nsknumubarccqe_pa     = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuMuBarCCQE )->Integral(1,                  ErecoBin1NuMuBar);
  double nsknumubarccqe_pb     = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuMuBarCCQE )->Integral(ErecoBin1NuMuBar+1, ErecoBin2NuMuBar);
  double nsknumubarccqe_pc     = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuMuBarCCQE )->Integral(ErecoBin2NuMuBar+1, conv::kNErecoBins1Rmu);

  // spectral numu ccqe (etrue)
  double nsknumuccqe_pt[conv::kXsNEtrueBins];
  for(int itb = 0; itb < NEtrueBinNuMu; itb++) {
    nsknumuccqe_pt[itb] = 0.;
    nsknumuccqe_pt[itb] = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBinNuMu[itb]+1, EtrueBinNuMu[itb+1]);
  }

  //
  // -1 sigma
  //

  if(sys_is_banff)
    skexp->TweakSystematic(systematic,+0.5);
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

  for (int im = 0; im < conv::knm; im++) {
    conv::SuperkMode_t m = conv::kSKReactionModes[s][im];
    if(m == conv::kSKmUndefined) continue;
    if(util::FindModeSlot(s,m) < 0) continue;
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
      else if(util::IsSignalNuE(m))
	nsknuesig_m += nskmodes_m[m];
      else {
	cerr << "CC mode " << AsString(m) << " is not numu, numubar, nue, nuebar, nuesig!" << endl;
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
  double nsknumuccqe_ma        = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuMuCCQE    )->Integral(1,                  ErecoBin1NuMu);
  double nsknumuccqe_mb        = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuMuCCQE    )->Integral(ErecoBin1NuMu+1,    ErecoBin2NuMu);
  double nsknumuccqe_mc        = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuMuCCQE    )->Integral(ErecoBin2NuMu+1,    conv::kNErecoBins1Rmu);
  double nsknumubarccqe_ma     = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuMuBarCCQE )->Integral(1,                  ErecoBin1NuMuBar);
  double nsknumubarccqe_mb     = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuMuBarCCQE )->Integral(ErecoBin1NuMuBar+1, ErecoBin2NuMuBar);
  double nsknumubarccqe_mc     = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuMuBarCCQE )->Integral(ErecoBin2NuMuBar+1, conv::kNErecoBins1Rmu);

  // spectral numu ccqe (etrue)
  double nsknumuccqe_mt[conv::kXsNEtrueBins];
  for(int itb = 0; itb < NEtrueBinNuMu; itb++) {
    nsknumuccqe_mt[itb] = 0.;
    nsknumuccqe_mt[itb] = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBinNuMu[itb]+1, EtrueBinNuMu[itb+1]);
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
  //spectral ccqe (etrue)
  double frac_nsknumuccqe_pt[conv::kXsNEtrueBins];
   for(int itb = 0; itb < NEtrueBinNuMu; itb++) {
     frac_nsknumuccqe_pt[itb]        = (nsknumuccqe_0t[itb]    > 0) ? 100.*(nsknumuccqe_pt[itb]   - nsknumuccqe_0t[itb]  ) / nsknumuccqe_0t[itb]   : 0.;
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
  //spectral ccqe (etrue)
  double frac_nsknumuccqe_mt[conv::kXsNEtrueBins];
   for(int itb = 0; itb < NEtrueBinNuMu; itb++) {
     frac_nsknumuccqe_mt[itb]        = (nsknumuccqe_0t[itb]    > 0) ? 100.*(nsknumuccqe_mt[itb]   - nsknumuccqe_0t[itb]  ) / nsknumuccqe_0t[itb]   : 0.;
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
     streamer << "\\footnotesize" << endl;
     streamer << "\\centering"     << endl;
     streamer << "\\begin{tabular}{| l | r | r  r | r  r |}" << endl;
     streamer << "\\hline" << endl;
     streamer << "{\\bf Component} &  "
          << "{\\bf Nominal}   &  "
          << "\\multicolumn {2}{c|} { {\\bf +1$\\sigma$} } &  "
          << "\\multicolumn {2}{c|} { {\\bf -1$\\sigma$} } \\\\"
          << endl;
     streamer << "\\hline" << endl;
  }
 
  string temp = "$N_{SK}$ (" + util::AsStringTex(s) + ")";
  streamer << setw(ncol1)  << setfill(' ') << temp
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
    //a normal individual mode table fill
    string component = "$N_{SK}$ (";
    component += util::AsStringTex(s);
    component += "): ";
    component +=  util::AsStringTex(m);
    streamer << setw(ncol1)  << setfill(' ') << component
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
    //spectral ccqe (ereco)
    if ( m == conv::kSKmNuMuCCQE && show_sk_ccqe_spectral_info_ereco) {
      temp = component + " ($E_{r} <$ 0.4 GeV)";
      streamer << setw(ncol1)  << setfill(' ') << temp
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
    }//numu ccqe spectral (ereco)
    if ( m == conv::kSKmNuMuCCQE && show_sk_ccqe_spectral_info_etrue) {
      char bufl[100], bufu[100];
      for(int itb = 0; itb < NEtrueBinNuMu; itb++) {
	snprintf(bufl, sizeof(bufl), "%g",  skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuCCQE, conv::kSKMCtFin)->GetYaxis()->GetBinLowEdge(EtrueBinNuMu[itb]+1));
	snprintf(bufu, sizeof(bufu), "%g",  skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuCCQE, conv::kSKMCtFin)->GetYaxis()->GetBinUpEdge(EtrueBinNuMu[itb+1]));
	temp = component + " (" + bufl + " $\\leq E_{t} <$ " + bufu + "GeV)";
	streamer << setw(ncol1)  << setfill(' ') << temp
	   << endcolumn
           << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumuccqe_0t[itb]
           << endcolumn
           << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumuccqe_pt[itb]
           << endcolumn
           << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknumuccqe_pt[itb] << pctsymbol 
           << endcolumn
           << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumuccqe_mt[itb]
           << endcolumn
           << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknumuccqe_mt[itb] << pctsymbol 
           << endrow
           << endl;
      }//loop over etrue bins
    }//numu ccqe spectral (etrue)
    if ( m == conv::kSKmNuMuBarCCQE && show_sk_ccqe_spectral_info_ereco) {
      temp = component + " ($E_{r} <$ 0.4 GeV)";
      streamer << setw(ncol1)  << setfill(' ') << temp
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
    }//numbar ccqe spectral (ereco)
    //total flavour-nc/cc
    if ((s == conv::kSKsFCFV1RingMuLike && m == conv::kSKmNuMuCCoth_mu) || (s == conv::kSKsFCFV1RingMuLike && m == conv::kSKmNuMuCCoth_e)) {
      temp = "$N_{SK}$ (" + util::AsStringTex(s) + "): $\\nu_{\\mu}$ CC all";
      streamer << setw(ncol1)  << setfill(' ') << temp
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
    }
    if (s == conv::kSKsFCFV1RingMuLike && m == conv::kSKmNuMuNCoth_e)) {
      temp = "$N_{SK}$ (" + util::AsStringTex(s) + "): $\\nu_{\\mu}$ NC all";
      streamer << setw(ncol1)  << setfill(' ') << temp
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
    }
    if ((s == conv::kSKsFCFV1RingMuLike && m == conv::kSKmNuMuBarCCoth_mu) || (s == conv::kSKsFCFV1RingMuLike && m == conv::kSKmNuMuBarCCoth_e)) {
      temp = "$N_{SK}$ (" + util::AsStringTex(s) + "): $\\bar{\\nu_{\\mu}}$ CC all";
      streamer << setw(ncol1)  << setfill(' ') << temp
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
    }
    if ((s == conv::kSKsFCFV1RingMuLike && m == conv::kSKmNuMuBarNCoth_e)) {
      temp = "$N_{SK}$ (" + util::AsStringTex(s) + "): $\\bar{\\nu_{\\mu}}$ NC all";
      streamer << setw(ncol1)  << setfill(' ') << temp
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
    }
    if ((s == conv::kSKsFCFV1RingMuLike && m == conv::kSKmNuECCoth_mu) || (s == conv::kSKsFCFV1RingMuLike && m == conv::kSKmNuECCoth_e)) {
      temp = "$N_{SK}$ (" + util::AsStringTex(s) + "): $\\nu_{e}$ CC all";
      streamer << setw(ncol1)  << setfill(' ') << temp
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
    }
    if ((s == conv::kSKsFCFV1RingMuLike && m == conv::kSKmNuENCoth_e)) {
      temp = "$N_{SK}$ (" + util::AsStringTex(s) + "): $\\nu_{e}$ NC all";
      streamer << setw(ncol1)  << setfill(' ') << temp
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
    }
    if ((s == conv::kSKsFCFV1RingMuLike && m == conv::kSKmNuEBarCCoth_mu) || (s == conv::kSKsFCFV1RingMuLike && m == conv::kSKmNuEBarCCoth_e)) {
      temp = "$N_{SK}$ (" + util::AsStringTex(s) + "): $\\bar{\\nu_{e}}$ CC all";
      streamer << setw(ncol1)  << setfill(' ') << temp
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
    }
    if ((s == conv::kSKsFCFV1RingMuLike && m == conv::kSKmNuEBarNCoth_e)) {
      temp = "$N_{SK}$ (" + util::AsStringTex(s) + "): $\\bar{\\nu_{e}}$ NC all";
      streamer << setw(ncol1)  << setfill(' ') << temp
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
    }
    if ((s == conv::kSKsFCFV1RingMuLike && m == conv::kSKmSigNuECCoth) || (s == conv::kSKsFCFV1RingMuLike && m == conv::kSKmSigNuECCoth)) {
      temp = "$N_{SK}$ (" + util::AsStringTex(s) + "): $\\nu_{e}$ (signal) CC all";
      streamer << setw(ncol1)  << setfill(' ') << temp
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
    }
  }//end loop over sk modes

  if(use_tex) {
     streamer << "\\hline" << endl;
  } else {
     streamer << setw(150)  << setfill('-') << " " << endl;
  }


  if(use_tex) {
     streamer << "\\hline" << endl;
     streamer << "\\end{tabular}" << endl;
     string systematicstr(systematic), hypothesisstr(hypothesis);
     size_t found;
     found=systematicstr.find_first_of("_");
     while(found!=string::npos) {
       systematicstr[found]='-';
       found=systematicstr.find_first_of("_",found+1);
     }
     found=hypothesisstr.find_first_of("_");
     while(found!=string::npos) {
       hypothesisstr[found]='-';
       found=hypothesisstr.find_first_of("_",found+1);
     }
     string banffstr = use_banff_bestfit ? "Using BANFF bestfit values." : "Using nominal MC (BANFF values set to 1)";
     streamer << "\\caption{ Effect of " << systematicstr << " for sample " << util::AsStringTex(s) 
	      << " with oscillation hypothesis " <<  hypothesisstr
	      << " using oscillation parameters: " 
	      << " $\\sin^22\\theta_{12} = $"  << input_mixing12
	      << "; $\\sin^22\\theta_{13} = $" << input_mixing13
	      << "; $\\sin^22\\theta_{23} = $" << input_mixing23
	      << "; $\\delta_{CP} = $        " << input_deltacp
	      << "; $\\Delta m^2_{12} = $    " << input_dm12sq
	      << "; $\\Delta m^2_{23} = $    " << input_dm23sq
	      << ". " << banffstr
	      << " }" << endl;
     string shortsamplestr = s==conv::kSKsFCFV1RingMuLike ? "1RmuT_" : "1ReT_";
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

    ofstream outfile (output_file,ios_base::app);

    // write to outfile
    outfile.write (streamer.str().c_str(),size);
  
    outfile.close();
  }

  cout << "Done!" << endl;

}
