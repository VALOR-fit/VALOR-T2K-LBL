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
//      - expected
//             This is a special option. Skip the main macro & print out a table with the 
//             expected 1 sigma percentage change for each systematic
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
// TD 23/03/12
//  Added spectral information pertaining to normalisation binning (CCQE and CC1pi)
//  Added "expected" option for systematics

#include <iomanip.h>
#include <TXMLEngine.h>

void print_expected_errors(bool, const char*, const char*);

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
  bool         use_tex                    = false,                          /* print-out results as a TeX table (for inclusion in tech. note)*/
  bool         show_sk_ccqe_spectral_info_ereco = false,                    /* */
  bool         show_sk_spectral_info_etrue_flux_binning = false,            /* */
  bool         show_sk_spectral_info_etrue_norm_binning = false,            /* */
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

  //string for convieniance
  string sys(systematic);

  //special mode
  if(sys == "expected") {
    print_expected_errors(use_tex, output_file, config);
    return;
  }

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

  // Get the EReco bin corresponding to the edges of SK CCQE efficiency systematic
  Long64_t ErecoBin1NuMuCCQE    = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuMuCCQE    )->FindBin(0.4-conv::kASmallNum);
  Long64_t ErecoBin2NuMuCCQE    = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuMuCCQE    )->FindBin(1.1-conv::kASmallNum);
  Long64_t ErecoBin1NuMuBarCCQE = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuMuBarCCQE )->FindBin(0.4-conv::kASmallNum);
  Long64_t ErecoBin2NuMuBarCCQE = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuMuBarCCQE )->FindBin(1.1-conv::kASmallNum);

  // Get the ETrue bin edges for XSec normalisation systematics
  Long64_t EtrueBin1CCQENorm = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuCCQE, conv::kSKMCtFin)->GetYaxis()->FindBin(conv::kXsEtrueBinEdges[util::ArrIdx::Instance()->SuperkSampleIdx(s)][util::ArrIdx::Instance()->SuperkModeIdx(s,conv::kSKmNuMuCCQE)][1.5-conv::kASmallNum]-conv::kASmallNum);
  Long64_t EtrueBin2CCQENorm = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuCCQE, conv::kSKMCtFin)->GetYaxis()->FindBin(conv::kXsEtrueBinEdges[util::ArrIdx::Instance()->SuperkSampleIdx(s)][util::ArrIdx::Instance()->SuperkModeIdx(s,conv::kSKmNuMuCCQE)][3.5-conv::kASmallNum]-conv::kASmallNum);
  Long64_t EtrueBin1CCResNorm = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuCCQE, conv::kSKMCtFin)->GetYaxis()->FindBin(conv::kXsEtrueBinEdges[util::ArrIdx::Instance()->SuperkSampleIdx(s)][util::ArrIdx::Instance()->SuperkModeIdx(s,conv::kSKmNuMuCCQE)][2.5-conv::kASmallNum]-conv::kASmallNum);

  // Get the ETrue bin edges for NuMuCCQE
  int NEtrueBinNuMu = conv::kXsNEtrueBins;
  Long64_t EtrueBinNuMu[conv::kXsNEtrueBins+1];
  for(int itb = 0; itb < conv::kXsNEtrueBins+1; itb++) {
    EtrueBinNuMu[itb] = -1;
    EtrueBinNuMu[itb] = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuCCQE, conv::kSKMCtFin)->GetYaxis()->FindBin(conv::kXsEtrueBinEdges[util::ArrIdx::Instance()->SuperkSampleIdx(s)][util::ArrIdx::Instance()->SuperkModeIdx(s,conv::kSKmNuMuCCQE)][itb]-conv::kASmallNum);
    if (EtrueBinNuMu[itb] > conv::kNEtrueBins) EtrueBinNuMu[itb] = conv::kNEtrueBins;
    if (EtrueBinNuMu[itb] < 1 && itb != 0) {NEtrueBinNuMu = itb-1; break;}
  }

  // Get the ETrue bin edges for NuMuBarCCQE
  int NEtrueBinNuMuBar = conv::kXsNEtrueBins;
  Long64_t EtrueBinNuMuBar[conv::kXsNEtrueBins+1];
  for(int itb = 0; itb < conv::kXsNEtrueBins+1; itb++) {
    EtrueBinNuMuBar[itb] = -1;
    EtrueBinNuMuBar[itb] = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuBarCCQE, conv::kSKMCtFin)->GetYaxis()->FindBin(conv::kXsEtrueBinEdges[util::ArrIdx::Instance()->SuperkSampleIdx(s)][util::ArrIdx::Instance()->SuperkModeIdx(s,conv::kSKmNuMuBarCCQE)][itb]-conv::kASmallNum);
    if (EtrueBinNuMuBar[itb] > conv::kNEtrueBins) EtrueBinNuMuBar[itb] = conv::kNEtrueBins;
    if (EtrueBinNuMuBar[itb] < 1 && itb != 0) {NEtrueBinNuMuBar = itb-1; break;}
  }

  // Get the ETrue bin edges for NuECCQE
  int NEtrueBinNuE = conv::kXsNEtrueBins;
  Long64_t EtrueBinNuE[conv::kXsNEtrueBins+1];
  for(int itb = 0; itb < conv::kXsNEtrueBins+1; itb++) {
    EtrueBinNuE[itb] = -1;
    EtrueBinNuE[itb] = skexp->GetSuperkMC()->Template(s, conv::kSKmNuECCQE, conv::kSKMCtFin)->GetYaxis()->FindBin(conv::kXsEtrueBinEdges[util::ArrIdx::Instance()->SuperkSampleIdx(s)][util::ArrIdx::Instance()->SuperkModeIdx(s,conv::kSKmNuECCQE)][itb]-conv::kASmallNum);
    if (EtrueBinNuE[itb] > conv::kNEtrueBins) EtrueBinNuE[itb] = conv::kNEtrueBins;
    if (EtrueBinNuE[itb] < 1 && itb != 0) {NEtrueBinNuE = itb-1; break;}
  }

  // Get the ETrue bin edges for NuEBarCCQE
  int NEtrueBinNuEBar = conv::kXsNEtrueBins;
  Long64_t EtrueBinNuEBar[conv::kXsNEtrueBins+1];
  for(int itb = 0; itb < conv::kXsNEtrueBins+1; itb++) {
    EtrueBinNuEBar[itb] = -1;
    EtrueBinNuEBar[itb] = skexp->GetSuperkMC()->Template(s, conv::kSKmNuEBarCCQE, conv::kSKMCtFin)->GetYaxis()->FindBin(conv::kXsEtrueBinEdges[util::ArrIdx::Instance()->SuperkSampleIdx(s)][util::ArrIdx::Instance()->SuperkModeIdx(s,conv::kSKmNuEBarCCQE)][itb]-conv::kASmallNum);
    if (EtrueBinNuEBar[itb] > conv::kNEtrueBins) EtrueBinNuEBar[itb] = conv::kNEtrueBins;
    if (EtrueBinNuEBar[itb] < 1 && itb != 0) {NEtrueBinNuEBar = itb-1; break;}
  }

  // Get the ETrue bin edges for SigNuECCQE
  int NEtrueBinSigNuE = conv::kXsNEtrueBins;
  Long64_t EtrueBinSigNuE[conv::kXsNEtrueBins+1];
  for(int itb = 0; itb < conv::kXsNEtrueBins+1; itb++) {
    EtrueBinSigNuE[itb] = -1;
    EtrueBinSigNuE[itb] = skexp->GetSuperkMC()->Template(s, conv::kSKmSigNuECCQE, conv::kSKMCtFin)->GetYaxis()->FindBin(conv::kXsEtrueBinEdges[util::ArrIdx::Instance()->SuperkSampleIdx(s)][util::ArrIdx::Instance()->SuperkModeIdx(s,conv::kSKmNuMuCCQE)][itb]-conv::kASmallNum);
    if (EtrueBinSigNuE[itb] > conv::kNEtrueBins) EtrueBinSigNuE[itb] = conv::kNEtrueBins;
    if (EtrueBinSigNuE[itb] < 1 && itb != 0) {NEtrueBinSigNuE = itb-1; break;}
  }

  //
  // nominal MC (pre BANFF)
  //

  //remove BANFF tuning
  skexp->GetBANFF()->RestoreNominalMC();

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
      else if(util::IsSignalNuE(m))
	nsknuesig_n += nskmodes_n[m];
      else {
	cerr << "CC mode " << AsString(m) << " is not numu, numubar, nue, nuebar, nuesig!" << endl;
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
  double nsknumuccqe_na        = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuMuCCQE    )->Integral(1,                  ErecoBin1NuMuCCQE);
  double nsknumuccqe_nb        = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuMuCCQE    )->Integral(ErecoBin1NuMuCCQE+1,    ErecoBin2NuMuCCQE);
  double nsknumuccqe_nc        = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuMuCCQE    )->Integral(ErecoBin2NuMuCCQE+1,    conv::kNErecoBins1Rmu);
  double nsknumubarccqe_na     = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuMuBarCCQE )->Integral(1,                  ErecoBin1NuMuBarCCQE);
  double nsknumubarccqe_nb     = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuMuBarCCQE )->Integral(ErecoBin1NuMuBarCCQE+1, ErecoBin2NuMuBarCCQE);
  double nsknumubarccqe_nc     = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuMuBarCCQE )->Integral(ErecoBin2NuMuBarCCQE+1, conv::kNErecoBins1Rmu);

  //spectral etrue (norm)
  double nsknumuccqe_etrue_na = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCQENorm);
  double nsknumuccqe_etrue_nb = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCQENorm+1, EtrueBin2CCQENorm);
  double nsknumuccqe_etrue_nc = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin2CCQENorm+1, conv::kNEtrueBins);
  double nsknumuccres_etrue_na = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuCC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCResNorm);
  double nsknumuccres_etrue_nb = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuCC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCResNorm+1, conv::kNEtrueBins);
  double nsknumubarccqe_etrue_na = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuBarCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCQENorm);
  double nsknumubarccqe_etrue_nb = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuBarCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCQENorm+1, EtrueBin2CCQENorm);
  double nsknumubarccqe_etrue_nc = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuBarCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin2CCQENorm+1, conv::kNEtrueBins);
  double nsknumubarccres_etrue_na = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuBarCC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCResNorm);
  double nsknumubarccres_etrue_nb = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuBarCC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCResNorm+1, conv::kNEtrueBins);
  double nsknueccqe_etrue_na = skexp->GetSuperkMC()->Template(s, conv::kSKmNuECCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCQENorm);
  double nsknueccqe_etrue_nb = skexp->GetSuperkMC()->Template(s, conv::kSKmNuECCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCQENorm+1, EtrueBin2CCQENorm);
  double nsknueccqe_etrue_nc = skexp->GetSuperkMC()->Template(s, conv::kSKmNuECCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin2CCQENorm+1, conv::kNEtrueBins);
  double nsknueccres_etrue_na = skexp->GetSuperkMC()->Template(s, conv::kSKmNuECC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCResNorm);
  double nsknueccres_etrue_nb = skexp->GetSuperkMC()->Template(s, conv::kSKmNuECC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCResNorm+1, conv::kNEtrueBins);
  double nsknuebarccqe_etrue_na = skexp->GetSuperkMC()->Template(s, conv::kSKmNuEBarCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCQENorm);
  double nsknuebarccqe_etrue_nb = skexp->GetSuperkMC()->Template(s, conv::kSKmNuEBarCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCQENorm+1, EtrueBin2CCQENorm);
  double nsknuebarccqe_etrue_nc = skexp->GetSuperkMC()->Template(s, conv::kSKmNuEBarCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin2CCQENorm+1, conv::kNEtrueBins);
  double nsknuebarccres_etrue_na = skexp->GetSuperkMC()->Template(s, conv::kSKmNuEBarCC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCResNorm);
  double nsknuebarccres_etrue_nb = skexp->GetSuperkMC()->Template(s, conv::kSKmNuEBarCC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCResNorm+1, conv::kNEtrueBins);
  double nsksignueccqe_etrue_na = skexp->GetSuperkMC()->Template(s, conv::kSKmSigNuECCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCQENorm);
  double nsksignueccqe_etrue_nb = skexp->GetSuperkMC()->Template(s, conv::kSKmSigNuECCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCQENorm+1, EtrueBin2CCQENorm);
  double nsksignueccqe_etrue_nc = skexp->GetSuperkMC()->Template(s, conv::kSKmSigNuECCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin2CCQENorm+1, conv::kNEtrueBins);
  double nsksignueccres_etrue_na = skexp->GetSuperkMC()->Template(s, conv::kSKmSigNuECC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCResNorm);
  double nsksignueccres_etrue_nb = skexp->GetSuperkMC()->Template(s, conv::kSKmSigNuECC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCResNorm+1, conv::kNEtrueBins);
  
  // spectral etrue (flux binning)
  double nsknumu_etrue_n[conv::kXsNEtrueBins], nsknumubar_etrue_n[conv::kXsNEtrueBins], nsknue_etrue_n[conv::kXsNEtrueBins], nsknuebar_etrue_n[conv::kXsNEtrueBins], nsksignue_etrue_n[conv::kXsNEtrueBins];
  for(int itb = 0; itb < NEtrueBinNuMu; itb++) {
    nsknumu_etrue_n[itb] = 0.;
    nsknumu_etrue_n[itb] = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBinNuMu[itb]+1, EtrueBinNuMu[itb+1]);
    cout << "numu etrue bin edges are from " << skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuCCQE, conv::kSKMCtFin)->GetYaxis()->GetBinLowEdge(EtrueBinNuMu[itb]+1) << " to " << skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuCCQE, conv::kSKMCtFin)->GetYaxis()->GetBinUpEdge(EtrueBinNuMu[itb+1]) << endl;;
  }
  for(int itb = 0; itb < NEtrueBinNuMuBar; itb++) {
    nsknumubar_etrue_n[itb] = 0.;
    nsknumubar_etrue_n[itb] = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuBarCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBinNuMuBar[itb]+1, EtrueBinNuMuBar[itb+1]);
  }
  for(int itb = 0; itb < NEtrueBinNuE; itb++) {
    nsknue_etrue_n[itb] = 0.;
    nsknue_etrue_n[itb] = skexp->GetSuperkMC()->Template(s, conv::kSKmNuECCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBinNuE[itb]+1, EtrueBinNuE[itb+1]);
  }
  for(int itb = 0; itb < NEtrueBinNuEBar; itb++) {
    nsknuebar_etrue_n[itb] = 0.;
    nsknuebar_etrue_n[itb] = skexp->GetSuperkMC()->Template(s, conv::kSKmNuEBarCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBinNuEBar[itb]+1, EtrueBinNuEBar[itb+1]);
  }
  for(int itb = 0; itb < NEtrueBinSigNuE; itb++) {
    nsksignue_etrue_n[itb] = 0.;
    nsksignue_etrue_n[itb] = skexp->GetSuperkMC()->Template(s, conv::kSKmSigNuECCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBinSigNuE[itb]+1, EtrueBinSigNuE[itb+1]);
  }

  //
  // nominal BANFF tuned MC
  //

  //restore BANFF tuning
  TVectorD twkv(skexp->GetBANFF()->BestFitValues());
  skexp->GetBANFF()->Tweak(true,&twkv);

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
  double nsknumuccqe_0a        = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuMuCCQE    )->Integral(1,                  ErecoBin1NuMuCCQE);
  double nsknumuccqe_0b        = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuMuCCQE    )->Integral(ErecoBin1NuMuCCQE+1,    ErecoBin2NuMuCCQE);
  double nsknumuccqe_0c        = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuMuCCQE    )->Integral(ErecoBin2NuMuCCQE+1,    conv::kNErecoBins1Rmu);
  double nsknumubarccqe_0a     = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuMuBarCCQE )->Integral(1,                  ErecoBin1NuMuBarCCQE);
  double nsknumubarccqe_0b     = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuMuBarCCQE )->Integral(ErecoBin1NuMuBarCCQE+1, ErecoBin2NuMuBarCCQE);
  double nsknumubarccqe_0c     = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuMuBarCCQE )->Integral(ErecoBin2NuMuBarCCQE+1, conv::kNErecoBins1Rmu);

  //spectral etrue (norm)
  double nsknumuccqe_etrue_0a = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCQENorm);
  double nsknumuccqe_etrue_0b = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCQENorm+1, EtrueBin2CCQENorm);
  double nsknumuccqe_etrue_0c = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin2CCQENorm+1, conv::kNEtrueBins);
  double nsknumuccres_etrue_0a = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuCC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCResNorm);
  double nsknumuccres_etrue_0b = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuCC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCResNorm+1, conv::kNEtrueBins);
  double nsknumubarccqe_etrue_0a = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuBarCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCQENorm);
  double nsknumubarccqe_etrue_0b = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuBarCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCQENorm+1, EtrueBin2CCQENorm);
  double nsknumubarccqe_etrue_0c = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuBarCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin2CCQENorm+1, conv::kNEtrueBins);
  double nsknumubarccres_etrue_0a = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuBarCC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCResNorm);
  double nsknumubarccres_etrue_0b = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuBarCC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCResNorm+1, conv::kNEtrueBins);
  double nsknueccqe_etrue_0a = skexp->GetSuperkMC()->Template(s, conv::kSKmNuECCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCQENorm);
  double nsknueccqe_etrue_0b = skexp->GetSuperkMC()->Template(s, conv::kSKmNuECCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCQENorm+1, EtrueBin2CCQENorm);
  double nsknueccqe_etrue_0c = skexp->GetSuperkMC()->Template(s, conv::kSKmNuECCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin2CCQENorm+1, conv::kNEtrueBins);
  double nsknueccres_etrue_0a = skexp->GetSuperkMC()->Template(s, conv::kSKmNuECC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCResNorm);
  double nsknueccres_etrue_0b = skexp->GetSuperkMC()->Template(s, conv::kSKmNuECC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCResNorm+1, conv::kNEtrueBins);
  double nsknuebarccqe_etrue_0a = skexp->GetSuperkMC()->Template(s, conv::kSKmNuEBarCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCQENorm);
  double nsknuebarccqe_etrue_0b = skexp->GetSuperkMC()->Template(s, conv::kSKmNuEBarCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCQENorm+1, EtrueBin2CCQENorm);
  double nsknuebarccqe_etrue_0c = skexp->GetSuperkMC()->Template(s, conv::kSKmNuEBarCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin2CCQENorm+1, conv::kNEtrueBins);
  double nsknuebarccres_etrue_0a = skexp->GetSuperkMC()->Template(s, conv::kSKmNuEBarCC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCResNorm);
  double nsknuebarccres_etrue_0b = skexp->GetSuperkMC()->Template(s, conv::kSKmNuEBarCC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCResNorm+1, conv::kNEtrueBins);
  double nsksignueccqe_etrue_0a = skexp->GetSuperkMC()->Template(s, conv::kSKmSigNuECCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCQENorm);
  double nsksignueccqe_etrue_0b = skexp->GetSuperkMC()->Template(s, conv::kSKmSigNuECCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCQENorm+1, EtrueBin2CCQENorm);
  double nsksignueccqe_etrue_0c = skexp->GetSuperkMC()->Template(s, conv::kSKmSigNuECCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin2CCQENorm+1, conv::kNEtrueBins);
  double nsksignueccres_etrue_0a = skexp->GetSuperkMC()->Template(s, conv::kSKmSigNuECC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCResNorm);
  double nsksignueccres_etrue_0b = skexp->GetSuperkMC()->Template(s, conv::kSKmSigNuECC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCResNorm+1, conv::kNEtrueBins);
  
  // spectral etrue (flux binning)
  double nsknumu_etrue_0[conv::kXsNEtrueBins], nsknumubar_etrue_0[conv::kXsNEtrueBins], nsknue_etrue_0[conv::kXsNEtrueBins], nsknuebar_etrue_0[conv::kXsNEtrueBins], nsksignue_etrue_0[conv::kXsNEtrueBins];
  for(int itb = 0; itb < NEtrueBinNuMu; itb++) {
    nsknumu_etrue_0[itb] = 0.;
    nsknumu_etrue_0[itb] = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBinNuMu[itb]+1, EtrueBinNuMu[itb+1]);
  }
  for(int itb = 0; itb < NEtrueBinNuMuBar; itb++) {
    nsknumubar_etrue_0[itb] = 0.;
    nsknumubar_etrue_0[itb] = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuBarCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBinNuMuBar[itb]+1, EtrueBinNuMuBar[itb+1]);
  }
  for(int itb = 0; itb < NEtrueBinNuE; itb++) {
    nsknue_etrue_0[itb] = 0.;
    nsknue_etrue_0[itb] = skexp->GetSuperkMC()->Template(s, conv::kSKmNuECCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBinNuE[itb]+1, EtrueBinNuE[itb+1]);
  }
  for(int itb = 0; itb < NEtrueBinNuEBar; itb++) {
    nsknuebar_etrue_0[itb] = 0.;
    nsknuebar_etrue_0[itb] = skexp->GetSuperkMC()->Template(s, conv::kSKmNuEBarCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBinNuEBar[itb]+1, EtrueBinNuEBar[itb+1]);
  }
  for(int itb = 0; itb < NEtrueBinSigNuE; itb++) {
    nsksignue_etrue_0[itb] = 0.;
    nsksignue_etrue_0[itb] = skexp->GetSuperkMC()->Template(s, conv::kSKmSigNuECCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBinSigNuE[itb]+1, EtrueBinSigNuE[itb+1]);
  }

  //
  // +1 sigma (relative to untuned BANFF)
  //
  cout << "\n\nApplying +1 sigma tweak" << endl;
  if(sys.find("banff_") == 0) {
    size_t found;
    found = sys.find_first_of("012");
    string bnum = sys.substr(found,2);
    int bindex = atoi(bnum.c_str());
    cout << bindex << endl;
    skexp->TweakSystematic(systematic, skexp->GetBANFF()->BestFitValue(util::PhysParamIdFromString(systematic)) + skexp->GetBANFF()->CovMatrix()(bindex,bindex));
    skexp->GetBANFF()->PrintSystematicVariations();
  }
  else if(sys == "sk_1Rmu_ccqe_ring_counting_eff_0")
    skexp->TweakSystematic(systematic,skexp->GetSuperkDetector()->Eff1RmuCCQERCCovarRelErr()(0,0));
  else if(sys == "sk_1Rmu_ccqe_ring_counting_eff_1")
    skexp->TweakSystematic(systematic,skexp->GetSuperkDetector()->Eff1RmuCCQERCCovarRelErr()(1,1));
  else if(sys == "sk_1Rmu_ccqe_ring_counting_eff_2")
    skexp->TweakSystematic(systematic,skexp->GetSuperkDetector()->Eff1RmuCCQERCCovarRelErr()(2,2));
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
  double nsknumuccqe_pa        = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuMuCCQE    )->Integral(1,                  ErecoBin1NuMuCCQE);
  double nsknumuccqe_pb        = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuMuCCQE    )->Integral(ErecoBin1NuMuCCQE+1,    ErecoBin2NuMuCCQE);
  double nsknumuccqe_pc        = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuMuCCQE    )->Integral(ErecoBin2NuMuCCQE+1,    conv::kNErecoBins1Rmu);
  double nsknumubarccqe_pa     = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuMuBarCCQE )->Integral(1,                  ErecoBin1NuMuBarCCQE);
  double nsknumubarccqe_pb     = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuMuBarCCQE )->Integral(ErecoBin1NuMuBarCCQE+1, ErecoBin2NuMuBarCCQE);
  double nsknumubarccqe_pc     = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuMuBarCCQE )->Integral(ErecoBin2NuMuBarCCQE+1, conv::kNErecoBins1Rmu);

  //spectral etrue (norm)
  double nsknumuccqe_etrue_pa = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCQENorm);
  double nsknumuccqe_etrue_pb = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCQENorm+1, EtrueBin2CCQENorm);
  double nsknumuccqe_etrue_pc = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin2CCQENorm+1, conv::kNEtrueBins);
  double nsknumuccres_etrue_pa = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuCC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCResNorm);
  double nsknumuccres_etrue_pb = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuCC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCResNorm+1, conv::kNEtrueBins);
  double nsknumubarccqe_etrue_pa = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuBarCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCQENorm);
  double nsknumubarccqe_etrue_pb = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuBarCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCQENorm+1, EtrueBin2CCQENorm);
  double nsknumubarccqe_etrue_pc = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuBarCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin2CCQENorm+1, conv::kNEtrueBins);
  double nsknumubarccres_etrue_pa = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuBarCC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCResNorm);
  double nsknumubarccres_etrue_pb = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuBarCC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCResNorm+1, conv::kNEtrueBins);
  double nsknueccqe_etrue_pa = skexp->GetSuperkMC()->Template(s, conv::kSKmNuECCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCQENorm);
  double nsknueccqe_etrue_pb = skexp->GetSuperkMC()->Template(s, conv::kSKmNuECCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCQENorm+1, EtrueBin2CCQENorm);
  double nsknueccqe_etrue_pc = skexp->GetSuperkMC()->Template(s, conv::kSKmNuECCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin2CCQENorm+1, conv::kNEtrueBins);
  double nsknueccres_etrue_pa = skexp->GetSuperkMC()->Template(s, conv::kSKmNuECC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCResNorm);
  double nsknueccres_etrue_pb = skexp->GetSuperkMC()->Template(s, conv::kSKmNuECC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCResNorm+1, conv::kNEtrueBins);
  double nsknuebarccqe_etrue_pa = skexp->GetSuperkMC()->Template(s, conv::kSKmNuEBarCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCQENorm);
  double nsknuebarccqe_etrue_pb = skexp->GetSuperkMC()->Template(s, conv::kSKmNuEBarCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCQENorm+1, EtrueBin2CCQENorm);
  double nsknuebarccqe_etrue_pc = skexp->GetSuperkMC()->Template(s, conv::kSKmNuEBarCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin2CCQENorm+1, conv::kNEtrueBins);
  double nsknuebarccres_etrue_pa = skexp->GetSuperkMC()->Template(s, conv::kSKmNuEBarCC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCResNorm);
  double nsknuebarccres_etrue_pb = skexp->GetSuperkMC()->Template(s, conv::kSKmNuEBarCC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCResNorm+1, conv::kNEtrueBins);
  double nsksignueccqe_etrue_pa = skexp->GetSuperkMC()->Template(s, conv::kSKmSigNuECCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCQENorm);
  double nsksignueccqe_etrue_pb = skexp->GetSuperkMC()->Template(s, conv::kSKmSigNuECCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCQENorm+1, EtrueBin2CCQENorm);
  double nsksignueccqe_etrue_pc = skexp->GetSuperkMC()->Template(s, conv::kSKmSigNuECCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin2CCQENorm+1, conv::kNEtrueBins);
  double nsksignueccres_etrue_pa = skexp->GetSuperkMC()->Template(s, conv::kSKmSigNuECC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCResNorm);
  double nsksignueccres_etrue_pb = skexp->GetSuperkMC()->Template(s, conv::kSKmSigNuECC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCResNorm+1, conv::kNEtrueBins);
  
  // spectral etrue (flux binning)
  double nsknumu_etrue_p[conv::kXsNEtrueBins], nsknumubar_etrue_p[conv::kXsNEtrueBins], nsknue_etrue_p[conv::kXsNEtrueBins], nsknuebar_etrue_p[conv::kXsNEtrueBins], nsksignue_etrue_p[conv::kXsNEtrueBins];
  for(int itb = 0; itb < NEtrueBinNuMu; itb++) {
    nsknumu_etrue_p[itb] = 0.;
    nsknumu_etrue_p[itb] = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBinNuMu[itb]+1, EtrueBinNuMu[itb+1]);
  }
  for(int itb = 0; itb < NEtrueBinNuMuBar; itb++) {
    nsknumubar_etrue_p[itb] = 0.;
    nsknumubar_etrue_p[itb] = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuBarCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBinNuMuBar[itb]+1, EtrueBinNuMuBar[itb+1]);
  }
  for(int itb = 0; itb < NEtrueBinNuE; itb++) {
    nsknue_etrue_p[itb] = 0.;
    nsknue_etrue_p[itb] = skexp->GetSuperkMC()->Template(s, conv::kSKmNuECCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBinNuE[itb]+1, EtrueBinNuE[itb+1]);
  }
  for(int itb = 0; itb < NEtrueBinNuEBar; itb++) {
    nsknuebar_etrue_p[itb] = 0.;
    nsknuebar_etrue_p[itb] = skexp->GetSuperkMC()->Template(s, conv::kSKmNuEBarCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBinNuEBar[itb]+1, EtrueBinNuEBar[itb+1]);
  }
  for(int itb = 0; itb < NEtrueBinSigNuE; itb++) {
    nsksignue_etrue_p[itb] = 0.;
    nsksignue_etrue_p[itb] = skexp->GetSuperkMC()->Template(s, conv::kSKmSigNuECCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBinSigNuE[itb]+1, EtrueBinSigNuE[itb+1]);
  }

  //
  // -1 sigma (relative to untuned BANFF)
  //
  cout << "\n\nApplying -1 sigma tweak" << endl;
  if(sys.find("banff_") == 0) {
    size_t found;
    found = sys.find_first_of("012");
    string bnum = sys.substr(found,2);
    int bindex = atoi(bnum.c_str());
    cout << bindex << endl;
    skexp->TweakSystematic(systematic, skexp->GetBANFF()->BestFitValue(util::PhysParamIdFromString(systematic)) - skexp->GetBANFF()->CovMatrix()(bindex,bindex));
    skexp->GetBANFF()->PrintSystematicVariations();
  }
  else if(sys == "sk_1Rmu_ccqe_ring_counting_eff_0")
    skexp->TweakSystematic(systematic,-skexp->GetSuperkDetector()->Eff1RmuCCQERCCovarRelErr()(0,0));
  else if(sys == "sk_1Rmu_ccqe_ring_counting_eff_1")
    skexp->TweakSystematic(systematic,-skexp->GetSuperkDetector()->Eff1RmuCCQERCCovarRelErr()(1,1));
  else if(sys == "sk_1Rmu_ccqe_ring_counting_eff_2")
    skexp->TweakSystematic(systematic,-skexp->GetSuperkDetector()->Eff1RmuCCQERCCovarRelErr()(2,2));
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
  double nsknumuccqe_ma        = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuMuCCQE    )->Integral(1,                  ErecoBin1NuMuCCQE);
  double nsknumuccqe_mb        = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuMuCCQE    )->Integral(ErecoBin1NuMuCCQE+1,    ErecoBin2NuMuCCQE);
  double nsknumuccqe_mc        = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuMuCCQE    )->Integral(ErecoBin2NuMuCCQE+1,    conv::kNErecoBins1Rmu);
  double nsknumubarccqe_ma     = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuMuBarCCQE )->Integral(1,                  ErecoBin1NuMuBarCCQE);
  double nsknumubarccqe_mb     = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuMuBarCCQE )->Integral(ErecoBin1NuMuBarCCQE+1, ErecoBin2NuMuBarCCQE);
  double nsknumubarccqe_mc     = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuMuBarCCQE )->Integral(ErecoBin2NuMuBarCCQE+1, conv::kNErecoBins1Rmu);

  //spectral etrue (norm)
  double nsknumuccqe_etrue_ma = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCQENorm);
  double nsknumuccqe_etrue_mb = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCQENorm+1, EtrueBin2CCQENorm);
  double nsknumuccqe_etrue_mc = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin2CCQENorm+1, conv::kNEtrueBins);
  double nsknumuccres_etrue_ma = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuCC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCResNorm);
  double nsknumuccres_etrue_mb = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuCC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCResNorm+1, conv::kNEtrueBins);
  double nsknumubarccqe_etrue_ma = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuBarCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCQENorm);
  double nsknumubarccqe_etrue_mb = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuBarCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCQENorm+1, EtrueBin2CCQENorm);
  double nsknumubarccqe_etrue_mc = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuBarCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin2CCQENorm+1, conv::kNEtrueBins);
  double nsknumubarccres_etrue_ma = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuBarCC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCResNorm);
  double nsknumubarccres_etrue_mb = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuBarCC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCResNorm+1, conv::kNEtrueBins);
  double nsknueccqe_etrue_ma = skexp->GetSuperkMC()->Template(s, conv::kSKmNuECCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCQENorm);
  double nsknueccqe_etrue_mb = skexp->GetSuperkMC()->Template(s, conv::kSKmNuECCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCQENorm+1, EtrueBin2CCQENorm);
  double nsknueccqe_etrue_mc = skexp->GetSuperkMC()->Template(s, conv::kSKmNuECCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin2CCQENorm+1, conv::kNEtrueBins);
  double nsknueccres_etrue_ma = skexp->GetSuperkMC()->Template(s, conv::kSKmNuECC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCResNorm);
  double nsknueccres_etrue_mb = skexp->GetSuperkMC()->Template(s, conv::kSKmNuECC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCResNorm+1, conv::kNEtrueBins);
  double nsknuebarccqe_etrue_ma = skexp->GetSuperkMC()->Template(s, conv::kSKmNuEBarCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCQENorm);
  double nsknuebarccqe_etrue_mb = skexp->GetSuperkMC()->Template(s, conv::kSKmNuEBarCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCQENorm+1, EtrueBin2CCQENorm);
  double nsknuebarccqe_etrue_mc = skexp->GetSuperkMC()->Template(s, conv::kSKmNuEBarCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin2CCQENorm+1, conv::kNEtrueBins);
  double nsknuebarccres_etrue_ma = skexp->GetSuperkMC()->Template(s, conv::kSKmNuEBarCC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCResNorm);
  double nsknuebarccres_etrue_mb = skexp->GetSuperkMC()->Template(s, conv::kSKmNuEBarCC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCResNorm+1, conv::kNEtrueBins);
  double nsksignueccqe_etrue_ma = skexp->GetSuperkMC()->Template(s, conv::kSKmSigNuECCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCQENorm);
  double nsksignueccqe_etrue_mb = skexp->GetSuperkMC()->Template(s, conv::kSKmSigNuECCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCQENorm+1, EtrueBin2CCQENorm);
  double nsksignueccqe_etrue_mc = skexp->GetSuperkMC()->Template(s, conv::kSKmSigNuECCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin2CCQENorm+1, conv::kNEtrueBins);
  double nsksignueccres_etrue_ma = skexp->GetSuperkMC()->Template(s, conv::kSKmSigNuECC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, 1, EtrueBin1CCResNorm);
  double nsksignueccres_etrue_mb = skexp->GetSuperkMC()->Template(s, conv::kSKmSigNuECC1pi, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBin1CCResNorm+1, conv::kNEtrueBins);
  
  // spectral etrue (flux binning)
  double nsknumu_etrue_m[conv::kXsNEtrueBins], nsknumubar_etrue_m[conv::kXsNEtrueBins], nsknue_etrue_m[conv::kXsNEtrueBins], nsknuebar_etrue_m[conv::kXsNEtrueBins], nsksignue_etrue_m[conv::kXsNEtrueBins];
  for(int itb = 0; itb < NEtrueBinNuMu; itb++) {
    nsknumu_etrue_m[itb] = 0.;
    nsknumu_etrue_m[itb] = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBinNuMu[itb]+1, EtrueBinNuMu[itb+1]);
  }
  for(int itb = 0; itb < NEtrueBinNuMuBar; itb++) {
    nsknumubar_etrue_m[itb] = 0.;
    nsknumubar_etrue_m[itb] = skexp->GetSuperkMC()->Template(s, conv::kSKmNuMuBarCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBinNuMuBar[itb]+1, EtrueBinNuMuBar[itb+1]);
  }
  for(int itb = 0; itb < NEtrueBinNuE; itb++) {
    nsknue_etrue_m[itb] = 0.;
    nsknue_etrue_m[itb] = skexp->GetSuperkMC()->Template(s, conv::kSKmNuECCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBinNuE[itb]+1, EtrueBinNuE[itb+1]);
  }
  for(int itb = 0; itb < NEtrueBinNuEBar; itb++) {
    nsknuebar_etrue_m[itb] = 0.;
    nsknuebar_etrue_m[itb] = skexp->GetSuperkMC()->Template(s, conv::kSKmNuEBarCCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBinNuEBar[itb]+1, EtrueBinNuEBar[itb+1]);
  }
  for(int itb = 0; itb < NEtrueBinSigNuE; itb++) {
    nsksignue_etrue_m[itb] = 0.;
    nsksignue_etrue_m[itb] = skexp->GetSuperkMC()->Template(s, conv::kSKmSigNuECCQE, conv::kSKMCtFin)->Integral(1, conv::kNErecoBins1Rmu, EtrueBinSigNuE[itb]+1, EtrueBinSigNuE[itb+1]);
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
  
  // spectral etrue (flux binning)
  double frac_nsknumu_etrue_p[conv::kXsNEtrueBins], frac_nsknumubar_etrue_p[conv::kXsNEtrueBins], frac_nsknue_etrue_p[conv::kXsNEtrueBins], frac_nsknuebar_etrue_p[conv::kXsNEtrueBins], frac_nsksignue_etrue_p[conv::kXsNEtrueBins];
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
  for(int itb = 0; itb < NEtrueBinSigNuE; itb++) {
     frac_nsksignue_etrue_p[itb]        = (nsksignue_etrue_0[itb]    > 0) ? 100.*(nsksignue_etrue_p[itb]   - nsksignue_etrue_0[itb]  ) / nsksignue_etrue_0[itb]   : 0.;
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
  
  // spectral etrue (flux binning)
  double frac_nsknumu_etrue_m[conv::kXsNEtrueBins], frac_nsknumubar_etrue_m[conv::kXsNEtrueBins], frac_nsknue_etrue_m[conv::kXsNEtrueBins], frac_nsknuebar_etrue_m[conv::kXsNEtrueBins], frac_nsksignue_etrue_m[conv::kXsNEtrueBins];
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
  for(int itb = 0; itb < NEtrueBinSigNuE; itb++) {
     frac_nsksignue_etrue_m[itb]        = (nsksignue_etrue_0[itb]    > 0) ? 100.*(nsksignue_etrue_m[itb]   - nsksignue_etrue_0[itb]  ) / nsksignue_etrue_0[itb]   : 0.;
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
     streamer << "\\scriptsize" << endl;
     //streamer << "\\centering"     << endl;
     streamer << "\\begin{tabular}{| l | r | r | r  r | r  r |}" << endl;
     streamer << "\\hline" << endl;
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
     streamer << "\\hline \\hline" << endl;
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
    //spectral ccqe (etrue)
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
      temp = component + " (1.5 $\\leq E_{t} <$ 2.5 GeV)";
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
      temp = component + " ($E_{t} \\geq$ 2.5 GeV)";
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
      temp = component + " (1.5 $\\leq E_{t} <$ 2.5 GeV)";
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
      temp = component + " ($E_{t} \\geq$ 2.5 GeV)";
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
      temp = component + " (1.5 $\\leq E_{t} <$ 2.5 GeV)";
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
      temp = component + " ($E_{t} \\geq$ 2.5 GeV)";
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
      temp = component + " (1.5 $\\leq E_{t} <$ 2.5 GeV)";
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
      temp = component + " ($E_{t} \\geq$ 2.5 GeV)";
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
    if ( m == conv::kSKmSigNuECCQE && show_sk_spectral_info_etrue_norm_binning) {
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
      temp = component + " (1.5 $\\leq E_{t} <$ 2.5 GeV)";
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
      temp = component + " ($E_{t} \\geq$ 2.5 GeV)";
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
    //spectral ccres (etrue)
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
    if ( m == conv::kSKmSigNuECC1pi && show_sk_spectral_info_etrue_norm_binning) {
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

    //spectral ccqe (ereco)
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
    }//numu ccqe spectral (ereco)
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
    }//numbar ccqe spectral (ereco)
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
    if (((s == conv::kSKsFCFV1RingMuLike && m == conv::kSKmNuMuNC) || (s == conv::kSKsFCFV1RingELike && m == conv::kSKmNuMuNCoth_e)) && show_sk_spectral_info_etrue_flux_binning) {
      if(use_tex) streamer << "\\hline" << endl;
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
    }//numu etrue (flux binning)
    if (((s == conv::kSKsFCFV1RingMuLike && m == conv::kSKmNuMuBarNC) || (s == conv::kSKsFCFV1RingELike && m == conv::kSKmNuMuBarNCoth_e)) && show_sk_spectral_info_etrue_flux_binning) {
      if(use_tex) streamer << "\\hline" << endl;
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
    }//numubar etrue (flux binning)
    if (((s == conv::kSKsFCFV1RingMuLike && m == conv::kSKmNuENC) || (s == conv::kSKsFCFV1RingELike && m == conv::kSKmNuENCoth_e)) && show_sk_spectral_info_etrue_flux_binning) {
      if(use_tex) streamer << "\\hline" << endl;
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
    }//nue etrue (flux binning)
    if (((s == conv::kSKsFCFV1RingMuLike && m == conv::kSKmNuEBarNC) || (s == conv::kSKsFCFV1RingELike && m == conv::kSKmNuEBarNCoth_e)) && show_sk_spectral_info_etrue_flux_binning) {
      if(use_tex) streamer << "\\hline" << endl;
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
    }//nuebar etrue (flux binning)
    //total flavour-nc/cc
    if ((s == conv::kSKsFCFV1RingMuLike && m == conv::kSKmNuMuCCoth_mu) || (s == conv::kSKsFCFV1RingELike && m == conv::kSKmNuMuCCoth_e)) {
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
    }
    if ((s == conv::kSKsFCFV1RingELike && m == conv::kSKmNuMuNCoth_e)) {
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
    }
    if ((s == conv::kSKsFCFV1RingMuLike && m == conv::kSKmNuMuBarCCoth_mu) || (s == conv::kSKsFCFV1RingELike && m == conv::kSKmNuMuBarCCoth_e)) {
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
    }
    if ((s == conv::kSKsFCFV1RingELike && m == conv::kSKmNuMuBarNCoth_e)) {
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
    }
    if ((s == conv::kSKsFCFV1RingMuLike && m == conv::kSKmNuECCoth_mu) || (s == conv::kSKsFCFV1RingELike && m == conv::kSKmNuECCoth_e)) {
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
    }
    if ((s == conv::kSKsFCFV1RingELike && m == conv::kSKmNuENCoth_e)) {
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
    }
    if ((s == conv::kSKsFCFV1RingMuLike && m == conv::kSKmNuEBarCCoth_mu) || (s == conv::kSKsFCFV1RingELike && m == conv::kSKmNuEBarCCoth_e)) {
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
    }
    if ((s == conv::kSKsFCFV1RingELike && m == conv::kSKmNuEBarNCoth_e)) {
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
    }
    if ((s == conv::kSKsFCFV1RingMuLike && m == conv::kSKmSigNuECCoth) || (s == conv::kSKsFCFV1RingELike && m == conv::kSKmSigNuECCoth)) {
      temp = "$N_{SK}$ (" + util::AsStringTex(s) + "): $\\nu_{e}$ (signal) CC all";
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
    }
    if (((s == conv::kSKsFCFV1RingMuLike && m == conv::kSKmSigNuECCoth) || (s == conv::kSKsFCFV1RingELike && m == conv::kSKmSigNuECCoth)) && show_sk_spectral_info_etrue_flux_binning) {
      if(use_tex) streamer << "\\hline" << endl;
      char bufl[100], bufu[100];
      for(int itb = 0; itb < NEtrueBinSigNuE; itb++) {
	snprintf(bufl, sizeof(bufl), "%g",  skexp->GetSuperkMC()->Template(s, conv::kSKmSigNuECCQE, conv::kSKMCtFin)->GetYaxis()->GetBinLowEdge(EtrueBinSigNuE[itb]+1));
	snprintf(bufu, sizeof(bufu), "%g",  skexp->GetSuperkMC()->Template(s, conv::kSKmSigNuECCQE, conv::kSKMCtFin)->GetYaxis()->GetBinUpEdge(EtrueBinSigNuE[itb+1]));
	temp = "$N_{SK}$ (" + util::AsStringTex(s) + "): $\\nu_{e}$ (signal) all (" + bufl + " $\\leq E_{t} <$ " + bufu + "GeV)";
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
    }//signue etrue (flux binning)
    if (use_tex && s == conv::kSKsFCFV1RingMuLike)
      switch(m) {
      case conv::kSKmNuMuNC:
      case conv::kSKmNuMuBarNC:
      case conv::kSKmNuENC:
      case conv::kSKmNuEBarNC:
	streamer << "\\hline" << endl;
      }
    else if (use_tex && s == conv::kSKsFCFV1RingELike)
      switch(m) {
      case conv::kSKmNuMuNCoth_e:
      case conv::kSKmNuMuBarNCoth_e:
      case conv::kSKmNuENCoth_e:
      case conv::kSKmNuEBarNCoth_e:
	streamer << "\\hline" << endl;
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
     if(sys.find("banff") == 0)
       sys = util::AsString(util::PhysParamIdFromString(systematic));
     size_t found;
     found=sys.find_first_of("_");
     while(found!=string::npos) {
       sys[found]='-';
       found=sys.find_first_of("_",found+1);
     }
     found=hypothesisstr.find_first_of("_");
     while(found!=string::npos) {
       hypothesisstr[found]='-';
       found=hypothesisstr.find_first_of("_",found+1);
     }
     streamer << "\\caption{ Effect of " << sys << " for sample " << util::AsStringTex(s) 
	      << " with oscillation hypothesis " <<  hypothesisstr
	      << " using oscillation parameters: " 
	      << " $\\sin^22\\theta_{12} = $"  << input_mixing12
	      << "; $\\sin^22\\theta_{13} = $" << input_mixing13
	      << "; $\\sin^22\\theta_{23} = $" << input_mixing23
	      << "; $\\delta_{CP} = $        " << input_deltacp
	      << "; $\\Delta m^2_{12} = $    " << input_dm12sq
	      << "; $\\Delta m^2_{23} = $    " << input_dm23sq
	      << ". }" << endl;
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

print_expected_errors(bool use_tex, const char * output_file, const char * config) {
  using namespace std;
  using namespace t2kana;
  
  // Parse the input XML configuration file
  TXMLEngine* xml = new TXMLEngine;   
  XMLDocPointer_t xmldoc = xml->ParseFile(config);
  if (xmldoc==0) {
    cerr << "Exp >> " << "Can not parse input: " << config << endl;
    delete xml;
    exit(1);
  }
  XMLNodePointer_t root_node = xml->DocGetRootElement(xmldoc);

  // Get the expected systematic variations
  double SKDetEnergyScale1RmuRelErr = 
    util::XML2Dbl(xml,root_node,
     "detector_systematics/superk/energy_scale_relative_error/one_mu_ring");
  double SKDetEff1RmuCCQEOthRelErr =
    util::XML2Dbl(xml,root_node,
      "detector_systematics/superk/efficiency_relative_error/one_mu_ring/ccqe/other_error");
  TMatrixD *SKDetEff1RmuCCQERCCovarRelErr = 
    util::XML2TMatrixD(xml,root_node,
      "detector_systematics/superk/efficiency_relative_error/one_mu_ring/ccqe/rc_error");
  double SKDetEff1RmuCCnonQERelErr   =
    util::XML2Dbl(xml,root_node,
      "detector_systematics/superk/efficiency_relative_error/one_mu_ring/ccnonqe");
  double SKDetEff1RmuNCRelErr        =
    util::XML2Dbl(xml,root_node,
      "detector_systematics/superk/efficiency_relative_error/one_mu_ring/nc");
  double SKDetEff1RmuNuECCRelErr =
    util::XML2Dbl(xml,root_node, 
      "detector_systematics/superk/efficiency_relative_error/one_mu_ring/nuecc");

  TMatrixD *BANFFCovMtrx = 
    util::XML2TMatrixD(xml,root_node,
      "banff/covariance_matrix");

  // Release memory
  xml->FreeDoc(xmldoc);
  delete xml;

  //print out
  string pctsymbol  = (use_tex) ? " \\% "   : "% ";
  string endcolumn  = (use_tex) ? " & "     : "" ;
  string endrow     = (use_tex) ? " \\\\ "  : "" ;
  int nprecisionnev = 4; // <-- for number of events
  int nprecisionpct = 2; // <-- for % changes
  int ncol1   = 80;
  int ncolnev = 15;
  int ncolfrc =  8;

  stringstream streamer;

  if(use_tex) {
     streamer << "\\begin{table}"   << endl;
     streamer << "\\scriptsize" << endl;
     //streamer << "\\centering"     << endl;
     streamer << "\\begin{tabular}{| l | r |}" << endl;
     streamer << "\\hline" << endl;
     streamer << "{\\bf Systematic} &  "
	  << "{\\bf Expected change}   "
          << endl;
     streamer << "\\hline" << endl;
  }
 
  string systematic;

  for (int i = 0; i < 3; i++) {
    systematic = Form("sk_1Rmu_ccqe_ring_counting_eff_%d",i);
    streamer << setw(ncol1)  << setfill(' ') << systematic << endcolumn
	     << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev)
	     << SKDetEff1RmuCCQERCCovarRelErr(i,i) * 100 << pctsymbol
	     << endrow << endl;
  }
  systematic = "sk_1Rmu_ccqe_other_eff";
  streamer << setw(ncol1)  << setfill(' ') << systematic << endcolumn
	   << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev)
	   << SKDetEff1RmuCCQEOthRelErr * 100 << pctsymbol
	   << endrow << endl;
  systematic = "sk_1Rmu_ccnonqe_eff";
  streamer << setw(ncol1)  << setfill(' ') << systematic << endcolumn
	   << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev)
	   << SKDetEff1RmuCCnonQERelErr * 100 << pctsymbol
	   << endrow << endl;
  systematic = "sk_1Rmu_nc_eff";
  streamer << setw(ncol1)  << setfill(' ') << systematic << endcolumn
	   << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev)
	   << SKDetEff1RmuNCRelErr * 100 << pctsymbol
	   << endrow << endl;
  systematic = "sk_1Rmu_nuecc_eff";
  streamer << setw(ncol1)  << setfill(' ') << systematic << endcolumn
	   << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev)
	   << SKDetEff1RmuNuECCRelErr * 100 << pctsymbol
	   << endrow << endl;
  systematic = "sk_energy_scale";
  streamer << setw(ncol1)  << setfill(' ') << systematic << endcolumn
	   << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev)
	   << SKDetEnergyScale1RmuRelErr * 100 << pctsymbol
	   << endrow << endl;
  for (int i = 0; i < 26; i++) {
    systematic = Form("banff_%d",i);
    streamer << setw(ncol1)  << setfill(' ') << systematic << endcolumn
	     << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev)
	     << BANFFCovMtrx(i,i) * 100 << pctsymbol
	     << endrow << endl;
  }

  if(!use_tex) {
     streamer << setw(150)  << setfill('-') << " " << endl;
  }

  if(use_tex) {
     streamer << "\\hline" << endl;
     streamer << "\\end{tabular}" << endl;
     streamer << "\\label{tab:expected}" << endl;
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

}
