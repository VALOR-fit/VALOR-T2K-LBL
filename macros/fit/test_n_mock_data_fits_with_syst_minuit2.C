//
// ------------------------------------------------------------------------------------------------
// Run numu disappearance fits to N toy experiments and store best fit values.
// The toy experiments are generated on the fly.
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//
/*
  #include <string>
  #include <sstream>
  #include <cstdlib>
  #include <TSystem.h>
  #include <TH1D.h>
  #include <TRandom3.h>
  #include <TObject.h>
  #include <TNtupleD.h>

  #include "Conventions.h"
  #include "SuperkExpectation.h"
  #include "FitParamList.h"
  #include "Utilities.h"
  #include "SuperkObservation.h"
  #include "FitInputs.h"
  #include "FitterABC.h"
  #include "Opt.h"
  #include "Superk1RmuErecoFitterMinuit2.h"
  #include "Superk1RmuErecoFitter.h"
*/

void test_n_mock_data_fits_with_syst_minuit2
(
 const char * tag                     = "toyexptfitswithsyst",          /* tag for output file */
 int          syst_mode               = 6,                              /* which systematics to tweak? 0 none. 1 all. 2 BANFF. 3 Unconstrained XSec. 4 SK detector. 5 FSI. 6 Selected. */
 int          run                     = 0,                              /* run number */
 int          nexperiments            = 1,                              /* number of toy MC experiment fits */
 double       input_mixing12          = -99,                            /* 12-mixing, sinq2theta or sinsqtheta depending on user choice, if -999 then randomize */
 double       input_mixing13          = -99,                            /* 13-mixing, sinq2theta or sinsqtheta depending on user choice, if -999 then randomize */
 double       input_mixing23          = -99,                            /* 23-mixing, sinq2theta or sinsqtheta depending on user choice */
 double       input_deltacp           = -99,                            /* deltacp, if -999 then randomize */
 double       input_dm12sq            = -99,                            /* Dm^2_{12}, if -999 then randomize */
 double       input_dm23sq            = -99,                            /* |Dm^2_{23}| */
 const char * config_file             = "analysis_configuration.xml",   /* analysis configuration file */
 bool         stat_fluctuations       = true,                           /* include statistical fluctuation in toy MC expts? */
 string       fit_param_list          = "nuosc3f_numudisapp_fix4+nuisanceselection", /* fit parameter list for the given physics hypothesis (see run_fitter.C) */
 const char * fit_option              = "minuit",                       /* fitter option (see run_fitter.C) */
 bool         store_chisq_surf        = false,                          /* by default storing only best-fit values, enable this to store the chisq surfaces   */
 bool         store_withsyst_chisq_surf = false,                        /* by default storing only best-fit values, enable this to store the chisq with systematics surfaces   */
 bool         store_toy_expts         = false,                          /* by default storing only best-fit values, enable this to store the toy MC expts     */
 bool         store_bestfit_spectrum  = false,                          /* by default storing only best-fit values, enable this to store the best-fit spectra */
 bool         use_minuit2             = true                            /* use minuit2, or minuit */
 )
{
  using namespace std;
  using namespace t2kana;
  using namespace t2kana::conv;

  // Global options
  Opt::Instance()->gOscMixingParamChoice = kNuMixSinsq2Theta;
  Opt::Instance()->gMassHierarchyAssumed = kNmhNormal;

  //Set to defaults
  if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsq2Theta) {
    if(input_mixing13 < -1 && input_mixing13 > -900) input_mixing13 = conv::kGlobFitSinsq2Theta13;
    if(input_mixing23 < -1 && input_mixing23 > -900) input_mixing23 = conv::kGlobFitSinsq2Theta23;
    if(input_mixing12 < -1 && input_mixing12 > -900) input_mixing12 = conv::kGlobFitSinsq2Theta12;
  }
  else {
    if(input_mixing13 < -1 && input_mixing13 > -900) input_mixing13 = conv::kGlobFitSinsqTheta13;
    if(input_mixing23 < -1 && input_mixing23 > -900) input_mixing23 = conv::kGlobFitSinsqTheta23;
    if(input_mixing12 < -1 && input_mixing12 > -900) input_mixing12 = conv::kGlobFitSinsqTheta12;
  }
  if(input_dm23sq < -1 && input_dm23sq > -900)  input_dm23sq = conv::kGlobFitDm23sq;
  if(input_dm12sq < -1 && input_dm12sq > -900)  input_dm12sq = conv::kGlobFitDm12sq;
  if(input_deltacp < -90 && input_deltacp > -900) input_deltacp = conv::kGlobFitDeltaCP;
  // Check whether some osc. params need to be randomized
  bool randomize_mixing12 = (TMath::Abs(input_mixing12 + 999) < kASmallNum);
  bool randomize_mixing13 = (TMath::Abs(input_mixing13 + 999) < kASmallNum);
  bool randomize_deltacp  = (TMath::Abs(input_deltacp  + 999) < kASmallNum);
  bool randomize_dm12sq   = (TMath::Abs(input_dm12sq   + 999) < kASmallNum);

  // Get syst mode string
  string syst_modestr = "";
  switch(syst_mode) {
  case 0:
    break;
  case 1:
    syst_modestr = "All";
    break;
  case 2:
    syst_modestr = "BANFF";
    break;
  case 3:
    syst_modestr = "Unconstrained XSec";
    break;
  case 4:
    syst_modestr = "SK detector";
    break;
  case 5:
    syst_modestr = "FSI";
    break;
  case 6:
    syst_modestr = "Selected";
    break;
  default:
    cerr << "Invalid systematic mode " << syst_mode << ". Exiting" << endl;
    exit(1);
    break;
  }

  // Print settings
  cout << "\n[Job inputs]" << endl;
  cout << "- Job: " << tag << "/" << run << endl;
  cout << "- Toy-MC generation" << endl;
  cout << "   - Number of toy experiments: " << nexperiments << endl; 
  cout << "   - Randomize systematics : "      << ((syst_mode > 0)     ?      "YES " : "NO") << syst_modestr << endl;
  cout << "   - Statistical fluctuations: "    << ((stat_fluctuations) ?      "YES " : "NO") << endl;
  cout << "   - True oscillation parameters: " << endl;
  cout << "       > " << util::AsString(kPhypMixing12) << " = " << ((randomize_mixing12) ? "* randomized *" : input_mixing12) << endl;
  cout << "       > " << util::AsString(kPhypMixing13) << " = " << ((randomize_mixing13) ? "* randomized *" : input_mixing13) << endl;
  cout << "       > " << util::AsString(kPhypMixing23) << " = " << input_mixing23 << endl;
  cout << "       > " << util::AsString(kPhypDeltaCP)  << " = " << ((randomize_deltacp) ? "* randomized *" : input_deltacp) << endl;
  cout << "       > " << util::AsString(kPhypDm12sq)   << " = " << ((randomize_dm12sq)  ? "* randomized *" : input_dm12sq ) << " eV^2/c^4" << endl;
  cout << "       > " << util::AsString(kPhypDm23sq)   << " = " << input_dm23sq << " eV^2/c^4" << endl;
  cout << "- Fit" << endl;
  cout << "   - Fit Parameter List: " << fit_param_list << endl;
  cout << "   - Fitter option: " << fit_option << endl;
  cout << "   - Minuit" << ((use_minuit2) ? "2" : "") << endl;
  cout << "- Outputs" << endl;
  cout << "     - Storing chisq surfaces: "           << ((store_chisq_surf) ?          "YES " : "NO") << endl;
  cout << "     - Storing chisq with syst surfaces: " << ((store_withsyst_chisq_surf) ? "YES " : "NO") << endl;
  cout << "     - Storing toy MC experiments:  "      << ((store_toy_expts) ?           "YES " : "NO") << endl;
  cout << "     - Storing best-fit spectra: "         << ((store_bestfit_spectrum) ?    "YES " : "NO") << endl;
  cout << endl;
  cout << *Opt::Instance() << endl;

  // Fit parameter list
  FitParamList * fparm = FitParamList::NamedParameterList(fit_param_list.c_str());
  int fit_syst_mode = 0;
  if(fit_param_list.find("nuisanceselection") != string::npos)
    fit_syst_mode = 1;
  else if(fit_param_list.find("nuisance") != string::npos)
    fit_syst_mode = 2;
  cout << " Using fit systematic mode " << fit_syst_mode << endl; 
  PhysHypothesis_t hyp = fparm->PhysicsHypothesis();
  cout << "- The following physics hypothesis is used in the fit: " << util::AsString(hyp) << endl;

  // Analysis configuration
  const char * config = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);
  cout << "- Reading analysis configuration from: " << config << endl;

  // Get superk spectrum predictor & configure
  SuperkExpectation * skexp = new SuperkExpectation;
  cout << config_file << endl;
  cout << config << endl;
  skexp->Read(config);

  // Create file to store results
  string oscstr = util::AsString(hyp,input_mixing12,input_mixing13,input_mixing23,input_deltacp,input_dm12sq,input_dm23sq);
  string syststr = syst_modestr;
  if(syststr == "")
    syststr = "none";
  string statstr = "none";
  if(stat_fluctuations)
    statstr = "yes";
  TFile f(Form("%s-run_%d-systfluc_%s-statfluc_%s-fitparam_%s-%s.root", tag, run, syststr.c_str(), statstr.c_str(), fit_param_list.c_str(), oscstr.c_str()), "recreate");

  // Create ntuple to store oscillation fit results
  TTree * nfitnt = new TTree("nfitnt","nfitnt");
  
  ///////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////
  // ***** variables for the branches *******
  ///////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////
  //fit info
  int iexp;
  int fitstatus;
  double bestfit_chisq;
  double true_chisq;
  double nsk_nom, nsk_bf;
  //inputs
  /* already initialised
     double input_mixing12;
     double input_mixing13;
     double input_mixing23;
     double input_deltacp;
     double input_dm12sq;
     double input_dm23sq;
  */
  double bestfit_mixing12, errminosplus_mixing12, errminosminus_mixing12, errmigrad_mixing12, errhesse_mixing12;
  double bestfit_mixing13, errminosplus_mixing13, errminosminus_mixing13, errmigrad_mixing13, errhesse_mixing13;
  double bestfit_mixing23, errminosplus_mixing23, errminosminus_mixing23, errmigrad_mixing23, errhesse_mixing23;
  double bestfit_deltacp, errminosplus_deltacp, errminosminus_deltacp, errmigrad_deltacp, errhesse_deltacp;
  double bestfit_dm12sq, errminosplus_dm12sq, errminosminus_dm12sq, errmigrad_dm12sq, errhesse_dm12sq;
  double bestfit_dm23sq, errminosplus_dm23sq, errminosminus_dm23sq, errmigrad_dm23sq, errhesse_dm23sq;
  double input_banff_00, bestfit_banff_00, errminosplus_banff_00, errminosminus_banff_00, errmigrad_banff_00, errhesse_banff_00;
  double input_banff_01, bestfit_banff_01, errminosplus_banff_01, errminosminus_banff_01, errmigrad_banff_01, errhesse_banff_01;
  double input_banff_02, bestfit_banff_02, errminosplus_banff_02, errminosminus_banff_02, errmigrad_banff_02, errhesse_banff_02;
  double input_banff_03, bestfit_banff_03, errminosplus_banff_03, errminosminus_banff_03, errmigrad_banff_03, errhesse_banff_03;
  double input_banff_04, bestfit_banff_04, errminosplus_banff_04, errminosminus_banff_04, errmigrad_banff_04, errhesse_banff_04;
  double input_banff_05, bestfit_banff_05, errminosplus_banff_05, errminosminus_banff_05, errmigrad_banff_05, errhesse_banff_05;
  double input_banff_06, bestfit_banff_06, errminosplus_banff_06, errminosminus_banff_06, errmigrad_banff_06, errhesse_banff_06;
  double input_banff_07, bestfit_banff_07, errminosplus_banff_07, errminosminus_banff_07, errmigrad_banff_07, errhesse_banff_07;
  double input_banff_08, bestfit_banff_08, errminosplus_banff_08, errminosminus_banff_08, errmigrad_banff_08, errhesse_banff_08;
  double input_banff_09, bestfit_banff_09, errminosplus_banff_09, errminosminus_banff_09, errmigrad_banff_09, errhesse_banff_09;
  double input_banff_10, bestfit_banff_10, errminosplus_banff_10, errminosminus_banff_10, errmigrad_banff_10, errhesse_banff_10;
  double input_banff_11, bestfit_banff_11, errminosplus_banff_11, errminosminus_banff_11, errmigrad_banff_11, errhesse_banff_11;
  double input_banff_12, bestfit_banff_12, errminosplus_banff_12, errminosminus_banff_12, errmigrad_banff_12, errhesse_banff_12;
  double input_banff_13, bestfit_banff_13, errminosplus_banff_13, errminosminus_banff_13, errmigrad_banff_13, errhesse_banff_13;
  double input_banff_14, bestfit_banff_14, errminosplus_banff_14, errminosminus_banff_14, errmigrad_banff_14, errhesse_banff_14;
  double input_banff_15, bestfit_banff_15, errminosplus_banff_15, errminosminus_banff_15, errmigrad_banff_15, errhesse_banff_15;
  double input_banff_maqe, bestfit_banff_maqe, errminosplus_banff_maqe, errminosminus_banff_maqe, errmigrad_banff_maqe, errhesse_banff_maqe;
  double input_banff_mares, bestfit_banff_mares, errminosplus_banff_mares, errminosminus_banff_mares, errmigrad_banff_mares, errhesse_banff_mares;
  double input_banff_ccqeE0, bestfit_banff_ccqeE0, errminosplus_banff_ccqeE0, errminosminus_banff_ccqeE0, errmigrad_banff_ccqeE0, errhesse_banff_ccqeE0;
  double input_banff_ccqeE1, bestfit_banff_ccqeE1, errminosplus_banff_ccqeE1, errminosminus_banff_ccqeE1, errmigrad_banff_ccqeE1, errhesse_banff_ccqeE1;
  double input_banff_ccqeE2, bestfit_banff_ccqeE2, errminosplus_banff_ccqeE2, errminosminus_banff_ccqeE2, errmigrad_banff_ccqeE2, errhesse_banff_ccqeE2;
  double input_banff_cc1piE0, bestfit_banff_cc1piE0, errminosplus_banff_cc1piE0, errminosminus_banff_cc1piE0, errmigrad_banff_cc1piE0, errhesse_banff_cc1piE0;
  double input_banff_cc1piE1, bestfit_banff_cc1piE1, errminosplus_banff_cc1piE1, errminosminus_banff_cc1piE1, errmigrad_banff_cc1piE1, errhesse_banff_cc1piE1;
  double input_unbanff_00, bestfit_unbanff_00, errminosplus_unbanff_00, errminosminus_unbanff_00, errmigrad_unbanff_00, errhesse_unbanff_00;
  double input_unbanff_01, bestfit_unbanff_01, errminosplus_unbanff_01, errminosminus_unbanff_01, errmigrad_unbanff_01, errhesse_unbanff_01;
  double input_unbanff_02, bestfit_unbanff_02, errminosplus_unbanff_02, errminosminus_unbanff_02, errmigrad_unbanff_02, errhesse_unbanff_02;
  double input_unbanff_03, bestfit_unbanff_03, errminosplus_unbanff_03, errminosminus_unbanff_03, errmigrad_unbanff_03, errhesse_unbanff_03;
  double input_unbanff_04, bestfit_unbanff_04, errminosplus_unbanff_04, errminosminus_unbanff_04, errmigrad_unbanff_04, errhesse_unbanff_04;
  double input_unbanff_05, bestfit_unbanff_05, errminosplus_unbanff_05, errminosminus_unbanff_05, errmigrad_unbanff_05, errhesse_unbanff_05;
  double input_unbanff_06, bestfit_unbanff_06, errminosplus_unbanff_06, errminosminus_unbanff_06, errmigrad_unbanff_06, errhesse_unbanff_06;
  double input_unbanff_07, bestfit_unbanff_07, errminosplus_unbanff_07, errminosminus_unbanff_07, errmigrad_unbanff_07, errhesse_unbanff_07;
  double input_unbanff_08, bestfit_unbanff_08, errminosplus_unbanff_08, errminosminus_unbanff_08, errmigrad_unbanff_08, errhesse_unbanff_08;
  double input_unbanff_09, bestfit_unbanff_09, errminosplus_unbanff_09, errminosminus_unbanff_09, errmigrad_unbanff_09, errhesse_unbanff_09;
  double input_unbanff_10, bestfit_unbanff_10, errminosplus_unbanff_10, errminosminus_unbanff_10, errmigrad_unbanff_10, errhesse_unbanff_10;
  double input_unbanff_11, bestfit_unbanff_11, errminosplus_unbanff_11, errminosminus_unbanff_11, errmigrad_unbanff_11, errhesse_unbanff_11;
  double input_unbanff_12, bestfit_unbanff_12, errminosplus_unbanff_12, errminosminus_unbanff_12, errmigrad_unbanff_12, errhesse_unbanff_12;
  double input_unbanff_13, bestfit_unbanff_13, errminosplus_unbanff_13, errminosminus_unbanff_13, errmigrad_unbanff_13, errhesse_unbanff_13;
  double input_unbanff_14, bestfit_unbanff_14, errminosplus_unbanff_14, errminosminus_unbanff_14, errmigrad_unbanff_14, errhesse_unbanff_14;
  double input_unbanff_15, bestfit_unbanff_15, errminosplus_unbanff_15, errminosminus_unbanff_15, errmigrad_unbanff_15, errhesse_unbanff_15;
  double input_unbanff_16, bestfit_unbanff_16, errminosplus_unbanff_16, errminosminus_unbanff_16, errmigrad_unbanff_16, errhesse_unbanff_16;
  double input_unbanff_17, bestfit_unbanff_17, errminosplus_unbanff_17, errminosminus_unbanff_17, errmigrad_unbanff_17, errhesse_unbanff_17;
  double input_unbanff_18, bestfit_unbanff_18, errminosplus_unbanff_18, errminosminus_unbanff_18, errmigrad_unbanff_18, errhesse_unbanff_18;
  double input_unbanff_19, bestfit_unbanff_19, errminosplus_unbanff_19, errminosminus_unbanff_19, errmigrad_unbanff_19, errhesse_unbanff_19;
  double input_unbanff_20, bestfit_unbanff_20, errminosplus_unbanff_20, errminosminus_unbanff_20, errmigrad_unbanff_20, errhesse_unbanff_20;
  double input_unbanff_21, bestfit_unbanff_21, errminosplus_unbanff_21, errminosminus_unbanff_21, errmigrad_unbanff_21, errhesse_unbanff_21;
  double input_unbanff_22, bestfit_unbanff_22, errminosplus_unbanff_22, errminosminus_unbanff_22, errmigrad_unbanff_22, errhesse_unbanff_22;
  double input_xsec_pf, bestfit_xsec_pf, errminosplus_xsec_pf, errminosminus_xsec_pf, errmigrad_xsec_pf, errhesse_xsec_pf;
  double input_xsec_sf, bestfit_xsec_sf, errminosplus_xsec_sf, errminosminus_xsec_sf, errmigrad_xsec_sf, errhesse_xsec_sf;
  double input_xsec_binde, bestfit_xsec_binde, errminosplus_xsec_binde, errminosminus_xsec_binde, errmigrad_xsec_binde, errhesse_xsec_binde;
  double input_xsec_wshape, bestfit_xsec_wshape, errminosplus_xsec_wshape, errminosminus_xsec_wshape, errmigrad_xsec_wshape, errhesse_xsec_wshape;
  double input_xsec_ccothshape, bestfit_xsec_ccothshape, errminosplus_xsec_ccothshape, errminosminus_xsec_ccothshape, errmigrad_xsec_ccothshape, errhesse_xsec_ccothshape;
  double input_xsec_cc1pishape, bestfit_xsec_cc1pishape, errminosplus_xsec_cc1pishape, errminosminus_xsec_cc1pishape, errmigrad_xsec_cc1pishape, errhesse_xsec_cc1pishape;
  double input_xsec_pilessdelta, bestfit_xsec_pilessdelta, errminosplus_xsec_pilessdelta, errminosminus_xsec_pilessdelta, errmigrad_xsec_pilessdelta, errhesse_xsec_pilessdelta;
  double input_xsec_normcccoh, bestfit_xsec_normcccoh, errminosplus_xsec_normcccoh, errminosminus_xsec_normcccoh, errmigrad_xsec_normcccoh, errhesse_xsec_normcccoh;
  double input_xsec_normnc1piC, bestfit_xsec_normnc1piC, errminosplus_xsec_normnc1piC, errminosminus_xsec_normnc1piC, errmigrad_xsec_normnc1piC, errhesse_xsec_normnc1piC;
  double input_xsec_normnccoh, bestfit_xsec_normnccoh, errminosplus_xsec_normnccoh, errminosminus_xsec_normnccoh, errmigrad_xsec_normnccoh, errhesse_xsec_normnccoh;
  double input_xsec_normncoth, bestfit_xsec_normncoth, errminosplus_xsec_normncoth, errminosminus_xsec_normncoth, errmigrad_xsec_normncoth, errhesse_xsec_normncoth;
  double input_xsec_normnuetonumu, bestfit_xsec_normnuetonumu, errminosplus_xsec_normnuetonumu, errminosminus_xsec_normnuetonumu, errmigrad_xsec_normnuetonumu, errhesse_xsec_normnuetonumu;
  double input_xsec_normantinutonu, bestfit_xsec_normantinutonu, errminosplus_xsec_normantinutonu, errminosminus_xsec_normantinutonu, errmigrad_xsec_normantinutonu, errhesse_xsec_normantinutonu;
  double input_skdet_energyscale, bestfit_skdet_energyscale, errminosplus_skdet_energyscale, errminosminus_skdet_energyscale, errmigrad_skdet_energyscale, errhesse_skdet_energyscale;
  double input_skdet_ccqerc0, bestfit_skdet_ccqerc0, errminosplus_skdet_ccqerc0, errminosminus_skdet_ccqerc0, errmigrad_skdet_ccqerc0, errhesse_skdet_ccqerc0;
  double input_skdet_ccqerc1, bestfit_skdet_ccqerc1, errminosplus_skdet_ccqerc1, errminosminus_skdet_ccqerc1, errmigrad_skdet_ccqerc1, errhesse_skdet_ccqerc1;
  double input_skdet_ccqerc2, bestfit_skdet_ccqerc2, errminosplus_skdet_ccqerc2, errminosminus_skdet_ccqerc2, errmigrad_skdet_ccqerc2, errhesse_skdet_ccqerc2;
  double input_skdet_ccqeoth, bestfit_skdet_ccqeoth, errminosplus_skdet_ccqeoth, errminosminus_skdet_ccqeoth, errmigrad_skdet_ccqeoth, errhesse_skdet_ccqeoth;
  double input_skdet_ccnqe, bestfit_skdet_ccnqe, errminosplus_skdet_ccnqe, errminosminus_skdet_ccnqe, errmigrad_skdet_ccnqe, errhesse_skdet_ccnqe;
  double input_skdet_nc, bestfit_skdet_nc, errminosplus_skdet_nc, errminosminus_skdet_nc, errmigrad_skdet_nc, errhesse_skdet_nc;
  double input_skdet_nuecc, bestfit_skdet_nuecc, errminosplus_skdet_nuecc, errminosminus_skdet_nuecc, errmigrad_skdet_nuecc, errhesse_skdet_nuecc;
  double input_fsi_000, bestfit_fsi_000, errminosplus_fsi_000, errminosminus_fsi_000, errmigrad_fsi_000, errhesse_fsi_000;
  double input_fsi_001, bestfit_fsi_001, errminosplus_fsi_001, errminosminus_fsi_001, errmigrad_fsi_001, errhesse_fsi_001;
  double input_fsi_002, bestfit_fsi_002, errminosplus_fsi_002, errminosminus_fsi_002, errmigrad_fsi_002, errhesse_fsi_002;
  double input_fsi_003, bestfit_fsi_003, errminosplus_fsi_003, errminosminus_fsi_003, errmigrad_fsi_003, errhesse_fsi_003;
  double input_fsi_004, bestfit_fsi_004, errminosplus_fsi_004, errminosminus_fsi_004, errmigrad_fsi_004, errhesse_fsi_004;
  double input_fsi_005, bestfit_fsi_005, errminosplus_fsi_005, errminosminus_fsi_005, errmigrad_fsi_005, errhesse_fsi_005;


  ///////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////
  // ***** setup the branches *******
  ///////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////
  //fit info
  nfitnt->Branch("iexp",&iexp,"iexp/I");
  nfitnt->Branch("fitstatus",&fitstatus,"fitstatus/I");
  nfitnt->Branch("true_chisq",&true_chisq,"true_chisq/D");
  nfitnt->Branch("bestfit_chisq",&bestfit_chisq,"bestfit_chisq/D");
  nfitnt->Branch("nsk_nom",&nsk_nom,"nsk_nom/D");
  nfitnt->Branch("nsk_bf",&nsk_bf,"nsk_bf/D");
  //inputs
  nfitnt->Branch("input_mixing12",&input_mixing12,"input_mixing12/D");
  nfitnt->Branch("input_mixing13",&input_mixing13,"input_mixing13/D");
  nfitnt->Branch("input_mixing23",&input_mixing23,"input_mixing23/D");
  nfitnt->Branch("input_deltacp",&input_deltacp,"input_deltacp/D");
  nfitnt->Branch("input_dm12sq",&input_dm12sq,"input_dm12sq/D");
  nfitnt->Branch("input_dm23sq",&input_dm23sq,"input_dm23sq/D");
  nfitnt->Branch("input_xsec_sf",&input_xsec_sf,"input_xsec_sf/D");
  nfitnt->Branch("input_xsec_binde",&input_xsec_binde,"input_xsec_binde/D");
  nfitnt->Branch("input_xsec_pilessdelta",&input_xsec_pilessdelta,"input_xsec_pilessdelta/D");
  nfitnt->Branch("input_skdet_ccnqe",&input_skdet_ccnqe,"input_skdet_ccnqe/D");
  nfitnt->Branch("input_skdet_nc",&input_skdet_nc,"input_skdet_nc/D");
  nfitnt->Branch("input_skdet_energyscale",&input_skdet_energyscale,"input_skdet_energyscale/D");
  nfitnt->Branch("input_banff_00",&input_banff_00,"input_banff_00/D");
  nfitnt->Branch("input_banff_01",&input_banff_01,"input_banff_01/D");
  nfitnt->Branch("input_banff_02",&input_banff_02,"input_banff_02/D");
  nfitnt->Branch("input_banff_03",&input_banff_03,"input_banff_03/D");
  nfitnt->Branch("input_banff_04",&input_banff_04,"input_banff_04/D");
  nfitnt->Branch("input_banff_05",&input_banff_05,"input_banff_05/D");
  nfitnt->Branch("input_banff_06",&input_banff_06,"input_banff_06/D");
  nfitnt->Branch("input_banff_07",&input_banff_07,"input_banff_07/D");
  nfitnt->Branch("input_banff_08",&input_banff_08,"input_banff_08/D");
  nfitnt->Branch("input_banff_09",&input_banff_09,"input_banff_09/D");
  nfitnt->Branch("input_banff_10",&input_banff_10,"input_banff_10/D");
  nfitnt->Branch("input_banff_11",&input_banff_11,"input_banff_11/D");
  nfitnt->Branch("input_banff_12",&input_banff_12,"input_banff_12/D");
  nfitnt->Branch("input_banff_13",&input_banff_13,"input_banff_13/D");
  nfitnt->Branch("input_banff_14",&input_banff_14,"input_banff_14/D");
  nfitnt->Branch("input_banff_15",&input_banff_15,"input_banff_15/D");
  nfitnt->Branch("input_banff_maqe",&input_banff_maqe,"input_banff_maqe/D");
  nfitnt->Branch("input_banff_mares",&input_banff_mares,"input_banff_mares/D");
  nfitnt->Branch("input_banff_ccqeE0",&input_banff_ccqeE0,"input_banff_ccqeE0/D");
  nfitnt->Branch("input_banff_ccqeE1",&input_banff_ccqeE1,"input_banff_ccqeE1/D");
  nfitnt->Branch("input_banff_ccqeE2",&input_banff_ccqeE2,"input_banff_ccqeE2/D");
  nfitnt->Branch("input_banff_cc1piE0",&input_banff_cc1piE0,"input_banff_cc1piE0/D");
  nfitnt->Branch("input_banff_cc1piE1",&input_banff_cc1piE1,"input_banff_cc1piE1/D");
  nfitnt->Branch("input_unbanff_00",&input_unbanff_00,"input_unbanff_00/D");
  nfitnt->Branch("input_unbanff_01",&input_unbanff_01,"input_unbanff_01/D");
  nfitnt->Branch("input_unbanff_02",&input_unbanff_02,"input_unbanff_02/D");
  nfitnt->Branch("input_unbanff_03",&input_unbanff_03,"input_unbanff_03/D");
  nfitnt->Branch("input_unbanff_04",&input_unbanff_04,"input_unbanff_04/D");
  nfitnt->Branch("input_unbanff_05",&input_unbanff_05,"input_unbanff_05/D");
  nfitnt->Branch("input_unbanff_06",&input_unbanff_06,"input_unbanff_06/D");
  nfitnt->Branch("input_unbanff_07",&input_unbanff_07,"input_unbanff_07/D");
  nfitnt->Branch("input_unbanff_08",&input_unbanff_08,"input_unbanff_08/D");
  nfitnt->Branch("input_unbanff_09",&input_unbanff_09,"input_unbanff_09/D");
  nfitnt->Branch("input_unbanff_10",&input_unbanff_10,"input_unbanff_10/D");
  nfitnt->Branch("input_unbanff_11",&input_unbanff_11,"input_unbanff_11/D");
  nfitnt->Branch("input_unbanff_12",&input_unbanff_12,"input_unbanff_12/D");
  nfitnt->Branch("input_unbanff_13",&input_unbanff_13,"input_unbanff_13/D");
  nfitnt->Branch("input_unbanff_14",&input_unbanff_14,"input_unbanff_14/D");
  nfitnt->Branch("input_unbanff_15",&input_unbanff_15,"input_unbanff_15/D");
  nfitnt->Branch("input_unbanff_16",&input_unbanff_16,"input_unbanff_16/D");
  nfitnt->Branch("input_unbanff_17",&input_unbanff_17,"input_unbanff_17/D");
  nfitnt->Branch("input_unbanff_18",&input_unbanff_18,"input_unbanff_18/D");
  nfitnt->Branch("input_unbanff_19",&input_unbanff_19,"input_unbanff_19/D");
  nfitnt->Branch("input_unbanff_20",&input_unbanff_20,"input_unbanff_20/D");
  nfitnt->Branch("input_unbanff_21",&input_unbanff_21,"input_unbanff_21/D");
  nfitnt->Branch("input_unbanff_22",&input_unbanff_22,"input_unbanff_22/D");
  nfitnt->Branch("input_xsec_pf",&input_xsec_pf,"input_xsec_pf/D");
  nfitnt->Branch("input_xsec_wshape",&input_xsec_wshape,"input_xsec_wshape/D");
  nfitnt->Branch("input_xsec_ccothshape",&input_xsec_ccothshape,"input_xsec_ccothshape/D");
  nfitnt->Branch("input_xsec_cc1pishape",&input_xsec_cc1pishape,"input_xsec_cc1pishape/D");
  nfitnt->Branch("input_xsec_normcccoh",&input_xsec_normcccoh,"input_xsec_normcccoh/D");
  nfitnt->Branch("input_xsec_normnc1piC",&input_xsec_normnc1piC,"input_xsec_normnc1piC/D");
  nfitnt->Branch("input_xsec_normnccoh",&input_xsec_normnccoh,"input_xsec_normnccoh/D");
  nfitnt->Branch("input_xsec_normncoth",&input_xsec_normncoth,"input_xsec_normncoth/D");
  nfitnt->Branch("input_xsec_normnuetonumu",&input_xsec_normnuetonumu,"input_xsec_normnuetonumu/D");
  nfitnt->Branch("input_xsec_normantinutonu",&input_xsec_normantinutonu,"input_xsec_normantinutonu/D");
  nfitnt->Branch("input_skdet_ccqerc0",&input_skdet_ccqerc0,"input_skdet_ccqerc0/D");
  nfitnt->Branch("input_skdet_ccqerc1",&input_skdet_ccqerc1,"input_skdet_ccqerc1/D");
  nfitnt->Branch("input_skdet_ccqerc2",&input_skdet_ccqerc2,"input_skdet_ccqerc2/D");
  nfitnt->Branch("input_skdet_ccqeoth",&input_skdet_ccqeoth,"input_skdet_ccqeoth/D");
  nfitnt->Branch("input_skdet_nuecc",&input_skdet_nuecc,"input_skdet_nuecc/D");
  nfitnt->Branch("input_fsi_000",&input_fsi_000,"input_fsi_000/D");
  nfitnt->Branch("input_fsi_001",&input_fsi_001,"input_fsi_001/D");
  nfitnt->Branch("input_fsi_002",&input_fsi_002,"input_fsi_002/D");
  nfitnt->Branch("input_fsi_003",&input_fsi_003,"input_fsi_003/D");
  nfitnt->Branch("input_fsi_004",&input_fsi_004,"input_fsi_004/D");
  nfitnt->Branch("input_fsi_005",&input_fsi_005,"input_fsi_005/D");
  //bestfits
  nfitnt->Branch("bestfit_mixing12",&bestfit_mixing12,"bestfit_mixing12/D");
  nfitnt->Branch("bestfit_mixing13",&bestfit_mixing13,"bestfit_mixing13/D");
  nfitnt->Branch("bestfit_mixing23",&bestfit_mixing23,"bestfit_mixing23/D");
  nfitnt->Branch("bestfit_deltacp",&bestfit_deltacp,"bestfit_deltacp/D");
  nfitnt->Branch("bestfit_dm12sq",&bestfit_dm12sq,"bestfit_dm12sq/D");
  nfitnt->Branch("bestfit_dm23sq",&bestfit_dm23sq,"bestfit_dm23sq/D");
  if(fit_syst_mode >= 1) {
    nfitnt->Branch("bestfit_xsec_sf",&bestfit_xsec_sf,"bestfit_xsec_sf/D");
    nfitnt->Branch("bestfit_xsec_binde",&bestfit_xsec_binde,"bestfit_xsec_binde/D");
    nfitnt->Branch("bestfit_xsec_pilessdelta",&bestfit_xsec_pilessdelta,"bestfit_xsec_pilessdelta/D");
    nfitnt->Branch("bestfit_skdet_ccnqe",&bestfit_skdet_ccnqe,"bestfit_skdet_ccnqe/D");
    nfitnt->Branch("bestfit_skdet_nc",&bestfit_skdet_nc,"bestfit_skdet_nc/D");
    nfitnt->Branch("bestfit_skdet_energyscale",&bestfit_skdet_energyscale,"bestfit_skdet_energyscale/D");
    nfitnt->Branch("bestfit_banff_00",&bestfit_banff_00,"bestfit_banff_00/D");
    nfitnt->Branch("bestfit_banff_01",&bestfit_banff_01,"bestfit_banff_01/D");
    nfitnt->Branch("bestfit_banff_02",&bestfit_banff_02,"bestfit_banff_02/D");
    nfitnt->Branch("bestfit_banff_03",&bestfit_banff_03,"bestfit_banff_03/D");
    nfitnt->Branch("bestfit_banff_04",&bestfit_banff_04,"bestfit_banff_04/D");
    nfitnt->Branch("bestfit_banff_05",&bestfit_banff_05,"bestfit_banff_05/D");
    nfitnt->Branch("bestfit_banff_06",&bestfit_banff_06,"bestfit_banff_06/D");
    nfitnt->Branch("bestfit_banff_07",&bestfit_banff_07,"bestfit_banff_07/D");
    nfitnt->Branch("bestfit_banff_08",&bestfit_banff_08,"bestfit_banff_08/D");
    nfitnt->Branch("bestfit_banff_09",&bestfit_banff_09,"bestfit_banff_09/D");
    nfitnt->Branch("bestfit_banff_10",&bestfit_banff_10,"bestfit_banff_10/D");
    nfitnt->Branch("bestfit_banff_11",&bestfit_banff_11,"bestfit_banff_11/D");
    nfitnt->Branch("bestfit_banff_12",&bestfit_banff_12,"bestfit_banff_12/D");
    nfitnt->Branch("bestfit_banff_13",&bestfit_banff_13,"bestfit_banff_13/D");
    nfitnt->Branch("bestfit_banff_14",&bestfit_banff_14,"bestfit_banff_14/D");
    nfitnt->Branch("bestfit_banff_15",&bestfit_banff_15,"bestfit_banff_15/D");
    nfitnt->Branch("bestfit_banff_maqe",&bestfit_banff_maqe,"bestfit_banff_maqe/D");
    nfitnt->Branch("bestfit_banff_mares",&bestfit_banff_mares,"bestfit_banff_mares/D");
    nfitnt->Branch("bestfit_banff_ccqeE0",&bestfit_banff_ccqeE0,"bestfit_banff_ccqeE0/D");
    nfitnt->Branch("bestfit_banff_ccqeE1",&bestfit_banff_ccqeE1,"bestfit_banff_ccqeE1/D");
    nfitnt->Branch("bestfit_banff_ccqeE2",&bestfit_banff_ccqeE2,"bestfit_banff_ccqeE2/D");
    nfitnt->Branch("bestfit_banff_cc1piE0",&bestfit_banff_cc1piE0,"bestfit_banff_cc1piE0/D");
    nfitnt->Branch("bestfit_banff_cc1piE1",&bestfit_banff_cc1piE1,"bestfit_banff_cc1piE1/D");
  }
  if(fit_syst_mode >= 1 && fit_param_list.find("unbanff") != string::npos) {
    nfitnt->Branch("bestfit_unbanff_00",&bestfit_unbanff_00,"bestfit_unbanff_00/D");
    nfitnt->Branch("bestfit_unbanff_01",&bestfit_unbanff_01,"bestfit_unbanff_01/D");
    nfitnt->Branch("bestfit_unbanff_02",&bestfit_unbanff_02,"bestfit_unbanff_02/D");
    nfitnt->Branch("bestfit_unbanff_03",&bestfit_unbanff_03,"bestfit_unbanff_03/D");
    nfitnt->Branch("bestfit_unbanff_04",&bestfit_unbanff_04,"bestfit_unbanff_04/D");
    nfitnt->Branch("bestfit_unbanff_05",&bestfit_unbanff_05,"bestfit_unbanff_05/D");
    nfitnt->Branch("bestfit_unbanff_06",&bestfit_unbanff_06,"bestfit_unbanff_06/D");
    nfitnt->Branch("bestfit_unbanff_07",&bestfit_unbanff_07,"bestfit_unbanff_07/D");
    nfitnt->Branch("bestfit_unbanff_08",&bestfit_unbanff_08,"bestfit_unbanff_08/D");
    nfitnt->Branch("bestfit_unbanff_09",&bestfit_unbanff_09,"bestfit_unbanff_09/D");
    nfitnt->Branch("bestfit_unbanff_10",&bestfit_unbanff_10,"bestfit_unbanff_10/D");
    nfitnt->Branch("bestfit_unbanff_11",&bestfit_unbanff_11,"bestfit_unbanff_11/D");
    nfitnt->Branch("bestfit_unbanff_12",&bestfit_unbanff_12,"bestfit_unbanff_12/D");
    nfitnt->Branch("bestfit_unbanff_13",&bestfit_unbanff_13,"bestfit_unbanff_13/D");
    nfitnt->Branch("bestfit_unbanff_14",&bestfit_unbanff_14,"bestfit_unbanff_14/D");
    nfitnt->Branch("bestfit_unbanff_15",&bestfit_unbanff_15,"bestfit_unbanff_15/D");
    nfitnt->Branch("bestfit_unbanff_16",&bestfit_unbanff_16,"bestfit_unbanff_16/D");
    nfitnt->Branch("bestfit_unbanff_17",&bestfit_unbanff_17,"bestfit_unbanff_17/D");
    nfitnt->Branch("bestfit_unbanff_18",&bestfit_unbanff_18,"bestfit_unbanff_18/D");
    nfitnt->Branch("bestfit_unbanff_19",&bestfit_unbanff_19,"bestfit_unbanff_19/D");
    nfitnt->Branch("bestfit_unbanff_20",&bestfit_unbanff_20,"bestfit_unbanff_20/D");
    nfitnt->Branch("bestfit_unbanff_21",&bestfit_unbanff_21,"bestfit_unbanff_21/D");
    nfitnt->Branch("bestfit_unbanff_22",&bestfit_unbanff_22,"bestfit_unbanff_22/D");
  }
  if(fit_syst_mode == 2) {
    nfitnt->Branch("bestfit_xsec_pf",&bestfit_xsec_pf,"bestfit_xsec_pf/D");
    nfitnt->Branch("bestfit_xsec_wshape",&bestfit_xsec_wshape,"bestfit_xsec_wshape/D");
    nfitnt->Branch("bestfit_xsec_ccothshape",&bestfit_xsec_ccothshape,"bestfit_xsec_ccothshape/D");
    nfitnt->Branch("bestfit_xsec_cc1pishape",&bestfit_xsec_cc1pishape,"bestfit_xsec_cc1pishape/D");
    nfitnt->Branch("bestfit_xsec_normcccoh",&bestfit_xsec_normcccoh,"bestfit_xsec_normcccoh/D");
    nfitnt->Branch("bestfit_xsec_normnc1piC",&bestfit_xsec_normnc1piC,"bestfit_xsec_normnc1piC/D");
    nfitnt->Branch("bestfit_xsec_normnccoh",&bestfit_xsec_normnccoh,"bestfit_xsec_normnccoh/D");
    nfitnt->Branch("bestfit_xsec_normncoth",&bestfit_xsec_normncoth,"bestfit_xsec_normncoth/D");
    nfitnt->Branch("bestfit_xsec_normnuetonumu",&bestfit_xsec_normnuetonumu,"bestfit_xsec_normnuetonumu/D");
    nfitnt->Branch("bestfit_xsec_normantinutonu",&bestfit_xsec_normantinutonu,"bestfit_xsec_normantinutonu/D");
    nfitnt->Branch("bestfit_skdet_ccqerc0",&bestfit_skdet_ccqerc0,"bestfit_skdet_ccqerc0/D");
    nfitnt->Branch("bestfit_skdet_ccqerc1",&bestfit_skdet_ccqerc1,"bestfit_skdet_ccqerc1/D");
    nfitnt->Branch("bestfit_skdet_ccqerc2",&bestfit_skdet_ccqerc2,"bestfit_skdet_ccqerc2/D");
    nfitnt->Branch("bestfit_skdet_ccqeoth",&bestfit_skdet_ccqeoth,"bestfit_skdet_ccqeoth/D");
    nfitnt->Branch("bestfit_skdet_nuecc",&bestfit_skdet_nuecc,"bestfit_skdet_nuecc/D");
    nfitnt->Branch("bestfit_fsi_000",&bestfit_fsi_000,"bestfit_fsi_000/D");
    nfitnt->Branch("bestfit_fsi_001",&bestfit_fsi_001,"bestfit_fsi_001/D");
    nfitnt->Branch("bestfit_fsi_002",&bestfit_fsi_002,"bestfit_fsi_002/D");
    nfitnt->Branch("bestfit_fsi_003",&bestfit_fsi_003,"bestfit_fsi_003/D");
    nfitnt->Branch("bestfit_fsi_004",&bestfit_fsi_004,"bestfit_fsi_004/D");
    nfitnt->Branch("bestfit_fsi_005",&bestfit_fsi_005,"bestfit_fsi_005/D");
  }
  //errminospluss
  nfitnt->Branch("errminosplus_mixing12",&errminosplus_mixing12,"errminosplus_mixing12/D");
  nfitnt->Branch("errminosplus_mixing13",&errminosplus_mixing13,"errminosplus_mixing13/D");
  nfitnt->Branch("errminosplus_mixing23",&errminosplus_mixing23,"errminosplus_mixing23/D");
  nfitnt->Branch("errminosplus_deltacp",&errminosplus_deltacp,"errminosplus_deltacp/D");
  nfitnt->Branch("errminosplus_dm12sq",&errminosplus_dm12sq,"errminosplus_dm12sq/D");
  nfitnt->Branch("errminosplus_dm23sq",&errminosplus_dm23sq,"errminosplus_dm23sq/D");
  if(fit_syst_mode >= 1) {
    nfitnt->Branch("errminosplus_xsec_sf",&errminosplus_xsec_sf,"errminosplus_xsec_sf/D");
    nfitnt->Branch("errminosplus_xsec_binde",&errminosplus_xsec_binde,"errminosplus_xsec_binde/D");
    nfitnt->Branch("errminosplus_xsec_pilessdelta",&errminosplus_xsec_pilessdelta,"errminosplus_xsec_pilessdelta/D");
    nfitnt->Branch("errminosplus_skdet_ccnqe",&errminosplus_skdet_ccnqe,"errminosplus_skdet_ccnqe/D");
    nfitnt->Branch("errminosplus_skdet_nc",&errminosplus_skdet_nc,"errminosplus_skdet_nc/D");
    nfitnt->Branch("errminosplus_skdet_energyscale",&errminosplus_skdet_energyscale,"errminosplus_skdet_energyscale/D");
    nfitnt->Branch("errminosplus_banff_00",&errminosplus_banff_00,"errminosplus_banff_00/D");
    nfitnt->Branch("errminosplus_banff_01",&errminosplus_banff_01,"errminosplus_banff_01/D");
    nfitnt->Branch("errminosplus_banff_02",&errminosplus_banff_02,"errminosplus_banff_02/D");
    nfitnt->Branch("errminosplus_banff_03",&errminosplus_banff_03,"errminosplus_banff_03/D");
    nfitnt->Branch("errminosplus_banff_04",&errminosplus_banff_04,"errminosplus_banff_04/D");
    nfitnt->Branch("errminosplus_banff_05",&errminosplus_banff_05,"errminosplus_banff_05/D");
    nfitnt->Branch("errminosplus_banff_06",&errminosplus_banff_06,"errminosplus_banff_06/D");
    nfitnt->Branch("errminosplus_banff_07",&errminosplus_banff_07,"errminosplus_banff_07/D");
    nfitnt->Branch("errminosplus_banff_08",&errminosplus_banff_08,"errminosplus_banff_08/D");
    nfitnt->Branch("errminosplus_banff_09",&errminosplus_banff_09,"errminosplus_banff_09/D");
    nfitnt->Branch("errminosplus_banff_10",&errminosplus_banff_10,"errminosplus_banff_10/D");
    nfitnt->Branch("errminosplus_banff_11",&errminosplus_banff_11,"errminosplus_banff_11/D");
    nfitnt->Branch("errminosplus_banff_12",&errminosplus_banff_12,"errminosplus_banff_12/D");
    nfitnt->Branch("errminosplus_banff_13",&errminosplus_banff_13,"errminosplus_banff_13/D");
    nfitnt->Branch("errminosplus_banff_14",&errminosplus_banff_14,"errminosplus_banff_14/D");
    nfitnt->Branch("errminosplus_banff_15",&errminosplus_banff_15,"errminosplus_banff_15/D");
    nfitnt->Branch("errminosplus_banff_maqe",&errminosplus_banff_maqe,"errminosplus_banff_maqe/D");
    nfitnt->Branch("errminosplus_banff_mares",&errminosplus_banff_mares,"errminosplus_banff_mares/D");
    nfitnt->Branch("errminosplus_banff_ccqeE0",&errminosplus_banff_ccqeE0,"errminosplus_banff_ccqeE0/D");
    nfitnt->Branch("errminosplus_banff_ccqeE1",&errminosplus_banff_ccqeE1,"errminosplus_banff_ccqeE1/D");
    nfitnt->Branch("errminosplus_banff_ccqeE2",&errminosplus_banff_ccqeE2,"errminosplus_banff_ccqeE2/D");
    nfitnt->Branch("errminosplus_banff_cc1piE0",&errminosplus_banff_cc1piE0,"errminosplus_banff_cc1piE0/D");
    nfitnt->Branch("errminosplus_banff_cc1piE1",&errminosplus_banff_cc1piE1,"errminosplus_banff_cc1piE1/D");
  }
  if(fit_syst_mode >= 1 && fit_param_list.find("unbanff") != string::npos) {
    nfitnt->Branch("errminosplus_unbanff_00",&errminosplus_unbanff_00,"errminosplus_unbanff_00/D");
    nfitnt->Branch("errminosplus_unbanff_01",&errminosplus_unbanff_01,"errminosplus_unbanff_01/D");
    nfitnt->Branch("errminosplus_unbanff_02",&errminosplus_unbanff_02,"errminosplus_unbanff_02/D");
    nfitnt->Branch("errminosplus_unbanff_03",&errminosplus_unbanff_03,"errminosplus_unbanff_03/D");
    nfitnt->Branch("errminosplus_unbanff_04",&errminosplus_unbanff_04,"errminosplus_unbanff_04/D");
    nfitnt->Branch("errminosplus_unbanff_05",&errminosplus_unbanff_05,"errminosplus_unbanff_05/D");
    nfitnt->Branch("errminosplus_unbanff_06",&errminosplus_unbanff_06,"errminosplus_unbanff_06/D");
    nfitnt->Branch("errminosplus_unbanff_07",&errminosplus_unbanff_07,"errminosplus_unbanff_07/D");
    nfitnt->Branch("errminosplus_unbanff_08",&errminosplus_unbanff_08,"errminosplus_unbanff_08/D");
    nfitnt->Branch("errminosplus_unbanff_09",&errminosplus_unbanff_09,"errminosplus_unbanff_09/D");
    nfitnt->Branch("errminosplus_unbanff_10",&errminosplus_unbanff_10,"errminosplus_unbanff_10/D");
    nfitnt->Branch("errminosplus_unbanff_11",&errminosplus_unbanff_11,"errminosplus_unbanff_11/D");
    nfitnt->Branch("errminosplus_unbanff_12",&errminosplus_unbanff_12,"errminosplus_unbanff_12/D");
    nfitnt->Branch("errminosplus_unbanff_13",&errminosplus_unbanff_13,"errminosplus_unbanff_13/D");
    nfitnt->Branch("errminosplus_unbanff_14",&errminosplus_unbanff_14,"errminosplus_unbanff_14/D");
    nfitnt->Branch("errminosplus_unbanff_15",&errminosplus_unbanff_15,"errminosplus_unbanff_15/D");
    nfitnt->Branch("errminosplus_unbanff_16",&errminosplus_unbanff_16,"errminosplus_unbanff_16/D");
    nfitnt->Branch("errminosplus_unbanff_17",&errminosplus_unbanff_17,"errminosplus_unbanff_17/D");
    nfitnt->Branch("errminosplus_unbanff_18",&errminosplus_unbanff_18,"errminosplus_unbanff_18/D");
    nfitnt->Branch("errminosplus_unbanff_19",&errminosplus_unbanff_19,"errminosplus_unbanff_19/D");
    nfitnt->Branch("errminosplus_unbanff_20",&errminosplus_unbanff_20,"errminosplus_unbanff_20/D");
    nfitnt->Branch("errminosplus_unbanff_21",&errminosplus_unbanff_21,"errminosplus_unbanff_21/D");
    nfitnt->Branch("errminosplus_unbanff_22",&errminosplus_unbanff_22,"errminosplus_unbanff_22/D");
  }
  if(fit_syst_mode == 2) {
    nfitnt->Branch("errminosplus_xsec_pf",&errminosplus_xsec_pf,"errminosplus_xsec_pf/D");
    nfitnt->Branch("errminosplus_xsec_wshape",&errminosplus_xsec_wshape,"errminosplus_xsec_wshape/D");
    nfitnt->Branch("errminosplus_xsec_ccothshape",&errminosplus_xsec_ccothshape,"errminosplus_xsec_ccothshape/D");
    nfitnt->Branch("errminosplus_xsec_cc1pishape",&errminosplus_xsec_cc1pishape,"errminosplus_xsec_cc1pishape/D");
    nfitnt->Branch("errminosplus_xsec_normcccoh",&errminosplus_xsec_normcccoh,"errminosplus_xsec_normcccoh/D");
    nfitnt->Branch("errminosplus_xsec_normnc1piC",&errminosplus_xsec_normnc1piC,"errminosplus_xsec_normnc1piC/D");
    nfitnt->Branch("errminosplus_xsec_normnccoh",&errminosplus_xsec_normnccoh,"errminosplus_xsec_normnccoh/D");
    nfitnt->Branch("errminosplus_xsec_normncoth",&errminosplus_xsec_normncoth,"errminosplus_xsec_normncoth/D");
    nfitnt->Branch("errminosplus_xsec_normnuetonumu",&errminosplus_xsec_normnuetonumu,"errminosplus_xsec_normnuetonumu/D");
    nfitnt->Branch("errminosplus_xsec_normantinutonu",&errminosplus_xsec_normantinutonu,"errminosplus_xsec_normantinutonu/D");
    nfitnt->Branch("errminosplus_skdet_ccqerc0",&errminosplus_skdet_ccqerc0,"errminosplus_skdet_ccqerc0/D");
    nfitnt->Branch("errminosplus_skdet_ccqerc1",&errminosplus_skdet_ccqerc1,"errminosplus_skdet_ccqerc1/D");
    nfitnt->Branch("errminosplus_skdet_ccqerc2",&errminosplus_skdet_ccqerc2,"errminosplus_skdet_ccqerc2/D");
    nfitnt->Branch("errminosplus_skdet_ccqeoth",&errminosplus_skdet_ccqeoth,"errminosplus_skdet_ccqeoth/D");
    nfitnt->Branch("errminosplus_skdet_nuecc",&errminosplus_skdet_nuecc,"errminosplus_skdet_nuecc/D");
    nfitnt->Branch("errminosplus_fsi_000",&errminosplus_fsi_000,"errminosplus_fsi_000/D");
    nfitnt->Branch("errminosplus_fsi_001",&errminosplus_fsi_001,"errminosplus_fsi_001/D");
    nfitnt->Branch("errminosplus_fsi_002",&errminosplus_fsi_002,"errminosplus_fsi_002/D");
    nfitnt->Branch("errminosplus_fsi_003",&errminosplus_fsi_003,"errminosplus_fsi_003/D");
    nfitnt->Branch("errminosplus_fsi_004",&errminosplus_fsi_004,"errminosplus_fsi_004/D");
    nfitnt->Branch("errminosplus_fsi_005",&errminosplus_fsi_005,"errminosplus_fsi_005/D");
  }
  //errms
  nfitnt->Branch("errminosminus_mixing12",&errminosminus_mixing12,"errminosminus_mixing12/D");
  nfitnt->Branch("errminosminus_mixing13",&errminosminus_mixing13,"errminosminus_mixing13/D");
  nfitnt->Branch("errminosminus_mixing23",&errminosminus_mixing23,"errminosminus_mixing23/D");
  nfitnt->Branch("errminosminus_deltacp",&errminosminus_deltacp,"errminosminus_deltacp/D");
  nfitnt->Branch("errminosminus_dm12sq",&errminosminus_dm12sq,"errminosminus_dm12sq/D");
  nfitnt->Branch("errminosminus_dm23sq",&errminosminus_dm23sq,"errminosminus_dm23sq/D");
  if(fit_syst_mode >= 1) {
    nfitnt->Branch("errminosminus_xsec_sf",&errminosminus_xsec_sf,"errminosminus_xsec_sf/D");
    nfitnt->Branch("errminosminus_xsec_binde",&errminosminus_xsec_binde,"errminosminus_xsec_binde/D");
    nfitnt->Branch("errminosminus_xsec_pilessdelta",&errminosminus_xsec_pilessdelta,"errminosminus_xsec_pilessdelta/D");
    nfitnt->Branch("errminosminus_skdet_ccnqe",&errminosminus_skdet_ccnqe,"errminosminus_skdet_ccnqe/D");
    nfitnt->Branch("errminosminus_skdet_nc",&errminosminus_skdet_nc,"errminosminus_skdet_nc/D");
    nfitnt->Branch("errminosminus_skdet_energyscale",&errminosminus_skdet_energyscale,"errminosminus_skdet_energyscale/D");
    nfitnt->Branch("errminosminus_banff_00",&errminosminus_banff_00,"errminosminus_banff_00/D");
    nfitnt->Branch("errminosminus_banff_01",&errminosminus_banff_01,"errminosminus_banff_01/D");
    nfitnt->Branch("errminosminus_banff_02",&errminosminus_banff_02,"errminosminus_banff_02/D");
    nfitnt->Branch("errminosminus_banff_03",&errminosminus_banff_03,"errminosminus_banff_03/D");
    nfitnt->Branch("errminosminus_banff_04",&errminosminus_banff_04,"errminosminus_banff_04/D");
    nfitnt->Branch("errminosminus_banff_05",&errminosminus_banff_05,"errminosminus_banff_05/D");
    nfitnt->Branch("errminosminus_banff_06",&errminosminus_banff_06,"errminosminus_banff_06/D");
    nfitnt->Branch("errminosminus_banff_07",&errminosminus_banff_07,"errminosminus_banff_07/D");
    nfitnt->Branch("errminosminus_banff_08",&errminosminus_banff_08,"errminosminus_banff_08/D");
    nfitnt->Branch("errminosminus_banff_09",&errminosminus_banff_09,"errminosminus_banff_09/D");
    nfitnt->Branch("errminosminus_banff_10",&errminosminus_banff_10,"errminosminus_banff_10/D");
    nfitnt->Branch("errminosminus_banff_11",&errminosminus_banff_11,"errminosminus_banff_11/D");
    nfitnt->Branch("errminosminus_banff_12",&errminosminus_banff_12,"errminosminus_banff_12/D");
    nfitnt->Branch("errminosminus_banff_13",&errminosminus_banff_13,"errminosminus_banff_13/D");
    nfitnt->Branch("errminosminus_banff_14",&errminosminus_banff_14,"errminosminus_banff_14/D");
    nfitnt->Branch("errminosminus_banff_15",&errminosminus_banff_15,"errminosminus_banff_15/D");
    nfitnt->Branch("errminosminus_banff_maqe",&errminosminus_banff_maqe,"errminosminus_banff_maqe/D");
    nfitnt->Branch("errminosminus_banff_mares",&errminosminus_banff_mares,"errminosminus_banff_mares/D");
    nfitnt->Branch("errminosminus_banff_ccqeE0",&errminosminus_banff_ccqeE0,"errminosminus_banff_ccqeE0/D");
    nfitnt->Branch("errminosminus_banff_ccqeE1",&errminosminus_banff_ccqeE1,"errminosminus_banff_ccqeE1/D");
    nfitnt->Branch("errminosminus_banff_ccqeE2",&errminosminus_banff_ccqeE2,"errminosminus_banff_ccqeE2/D");
    nfitnt->Branch("errminosminus_banff_cc1piE0",&errminosminus_banff_cc1piE0,"errminosminus_banff_cc1piE0/D");
    nfitnt->Branch("errminosminus_banff_cc1piE1",&errminosminus_banff_cc1piE1,"errminosminus_banff_cc1piE1/D");
  }
  if(fit_syst_mode >= 1 && fit_param_list.find("unbanff") != string::npos) {
    nfitnt->Branch("errminosminus_unbanff_00",&errminosminus_unbanff_00,"errminosminus_unbanff_00/D");
    nfitnt->Branch("errminosminus_unbanff_01",&errminosminus_unbanff_01,"errminosminus_unbanff_01/D");
    nfitnt->Branch("errminosminus_unbanff_02",&errminosminus_unbanff_02,"errminosminus_unbanff_02/D");
    nfitnt->Branch("errminosminus_unbanff_03",&errminosminus_unbanff_03,"errminosminus_unbanff_03/D");
    nfitnt->Branch("errminosminus_unbanff_04",&errminosminus_unbanff_04,"errminosminus_unbanff_04/D");
    nfitnt->Branch("errminosminus_unbanff_05",&errminosminus_unbanff_05,"errminosminus_unbanff_05/D");
    nfitnt->Branch("errminosminus_unbanff_06",&errminosminus_unbanff_06,"errminosminus_unbanff_06/D");
    nfitnt->Branch("errminosminus_unbanff_07",&errminosminus_unbanff_07,"errminosminus_unbanff_07/D");
    nfitnt->Branch("errminosminus_unbanff_08",&errminosminus_unbanff_08,"errminosminus_unbanff_08/D");
    nfitnt->Branch("errminosminus_unbanff_09",&errminosminus_unbanff_09,"errminosminus_unbanff_09/D");
    nfitnt->Branch("errminosminus_unbanff_10",&errminosminus_unbanff_10,"errminosminus_unbanff_10/D");
    nfitnt->Branch("errminosminus_unbanff_11",&errminosminus_unbanff_11,"errminosminus_unbanff_11/D");
    nfitnt->Branch("errminosminus_unbanff_12",&errminosminus_unbanff_12,"errminosminus_unbanff_12/D");
    nfitnt->Branch("errminosminus_unbanff_13",&errminosminus_unbanff_13,"errminosminus_unbanff_13/D");
    nfitnt->Branch("errminosminus_unbanff_14",&errminosminus_unbanff_14,"errminosminus_unbanff_14/D");
    nfitnt->Branch("errminosminus_unbanff_15",&errminosminus_unbanff_15,"errminosminus_unbanff_15/D");
    nfitnt->Branch("errminosminus_unbanff_16",&errminosminus_unbanff_16,"errminosminus_unbanff_16/D");
    nfitnt->Branch("errminosminus_unbanff_17",&errminosminus_unbanff_17,"errminosminus_unbanff_17/D");
    nfitnt->Branch("errminosminus_unbanff_18",&errminosminus_unbanff_18,"errminosminus_unbanff_18/D");
    nfitnt->Branch("errminosminus_unbanff_19",&errminosminus_unbanff_19,"errminosminus_unbanff_19/D");
    nfitnt->Branch("errminosminus_unbanff_20",&errminosminus_unbanff_20,"errminosminus_unbanff_20/D");
    nfitnt->Branch("errminosminus_unbanff_21",&errminosminus_unbanff_21,"errminosminus_unbanff_21/D");
    nfitnt->Branch("errminosminus_unbanff_22",&errminosminus_unbanff_22,"errminosminus_unbanff_22/D");
  }
  if(fit_syst_mode == 2) {
    nfitnt->Branch("errminosminus_xsec_pf",&errminosminus_xsec_pf,"errminosminus_xsec_pf/D");
    nfitnt->Branch("errminosminus_xsec_wshape",&errminosminus_xsec_wshape,"errminosminus_xsec_wshape/D");
    nfitnt->Branch("errminosminus_xsec_ccothshape",&errminosminus_xsec_ccothshape,"errminosminus_xsec_ccothshape/D");
    nfitnt->Branch("errminosminus_xsec_cc1pishape",&errminosminus_xsec_cc1pishape,"errminosminus_xsec_cc1pishape/D");
    nfitnt->Branch("errminosminus_xsec_normcccoh",&errminosminus_xsec_normcccoh,"errminosminus_xsec_normcccoh/D");
    nfitnt->Branch("errminosminus_xsec_normnc1piC",&errminosminus_xsec_normnc1piC,"errminosminus_xsec_normnc1piC/D");
    nfitnt->Branch("errminosminus_xsec_normnccoh",&errminosminus_xsec_normnccoh,"errminosminus_xsec_normnccoh/D");
    nfitnt->Branch("errminosminus_xsec_normncoth",&errminosminus_xsec_normncoth,"errminosminus_xsec_normncoth/D");
    nfitnt->Branch("errminosminus_xsec_normnuetonumu",&errminosminus_xsec_normnuetonumu,"errminosminus_xsec_normnuetonumu/D");
    nfitnt->Branch("errminosminus_xsec_normantinutonu",&errminosminus_xsec_normantinutonu,"errminosminus_xsec_normantinutonu/D");
    nfitnt->Branch("errminosminus_skdet_ccqerc0",&errminosminus_skdet_ccqerc0,"errminosminus_skdet_ccqerc0/D");
    nfitnt->Branch("errminosminus_skdet_ccqerc1",&errminosminus_skdet_ccqerc1,"errminosminus_skdet_ccqerc1/D");
    nfitnt->Branch("errminosminus_skdet_ccqerc2",&errminosminus_skdet_ccqerc2,"errminosminus_skdet_ccqerc2/D");
    nfitnt->Branch("errminosminus_skdet_ccqeoth",&errminosminus_skdet_ccqeoth,"errminosminus_skdet_ccqeoth/D");
    nfitnt->Branch("errminosminus_skdet_nuecc",&errminosminus_skdet_nuecc,"errminosminus_skdet_nuecc/D");
    nfitnt->Branch("errminosminus_fsi_000",&errminosminus_fsi_000,"errminosminus_fsi_000/D");
    nfitnt->Branch("errminosminus_fsi_001",&errminosminus_fsi_001,"errminosminus_fsi_001/D");
    nfitnt->Branch("errminosminus_fsi_002",&errminosminus_fsi_002,"errminosminus_fsi_002/D");
    nfitnt->Branch("errminosminus_fsi_003",&errminosminus_fsi_003,"errminosminus_fsi_003/D");
    nfitnt->Branch("errminosminus_fsi_004",&errminosminus_fsi_004,"errminosminus_fsi_004/D");
    nfitnt->Branch("errminosminus_fsi_005",&errminosminus_fsi_005,"errminosminus_fsi_005/D");
  }
  //errmigrads
  nfitnt->Branch("errmigrad_mixing12",&errmigrad_mixing12,"errmigrad_mixing12/D");
  nfitnt->Branch("errmigrad_mixing13",&errmigrad_mixing13,"errmigrad_mixing13/D");
  nfitnt->Branch("errmigrad_mixing23",&errmigrad_mixing23,"errmigrad_mixing23/D");
  nfitnt->Branch("errmigrad_deltacp",&errmigrad_deltacp,"errmigrad_deltacp/D");
  nfitnt->Branch("errmigrad_dm12sq",&errmigrad_dm12sq,"errmigrad_dm12sq/D");
  nfitnt->Branch("errmigrad_dm23sq",&errmigrad_dm23sq,"errmigrad_dm23sq/D");
  if(fit_syst_mode >= 1) {
    nfitnt->Branch("errmigrad_xsec_sf",&errmigrad_xsec_sf,"errmigrad_xsec_sf/D");
    nfitnt->Branch("errmigrad_xsec_binde",&errmigrad_xsec_binde,"errmigrad_xsec_binde/D");
    nfitnt->Branch("errmigrad_xsec_pilessdelta",&errmigrad_xsec_pilessdelta,"errmigrad_xsec_pilessdelta/D");
    nfitnt->Branch("errmigrad_skdet_ccnqe",&errmigrad_skdet_ccnqe,"errmigrad_skdet_ccnqe/D");
    nfitnt->Branch("errmigrad_skdet_nc",&errmigrad_skdet_nc,"errmigrad_skdet_nc/D");
    nfitnt->Branch("errmigrad_skdet_energyscale",&errmigrad_skdet_energyscale,"errmigrad_skdet_energyscale/D");
    nfitnt->Branch("errmigrad_banff_00",&errmigrad_banff_00,"errmigrad_banff_00/D");
    nfitnt->Branch("errmigrad_banff_01",&errmigrad_banff_01,"errmigrad_banff_01/D");
    nfitnt->Branch("errmigrad_banff_02",&errmigrad_banff_02,"errmigrad_banff_02/D");
    nfitnt->Branch("errmigrad_banff_03",&errmigrad_banff_03,"errmigrad_banff_03/D");
    nfitnt->Branch("errmigrad_banff_04",&errmigrad_banff_04,"errmigrad_banff_04/D");
    nfitnt->Branch("errmigrad_banff_05",&errmigrad_banff_05,"errmigrad_banff_05/D");
    nfitnt->Branch("errmigrad_banff_06",&errmigrad_banff_06,"errmigrad_banff_06/D");
    nfitnt->Branch("errmigrad_banff_07",&errmigrad_banff_07,"errmigrad_banff_07/D");
    nfitnt->Branch("errmigrad_banff_08",&errmigrad_banff_08,"errmigrad_banff_08/D");
    nfitnt->Branch("errmigrad_banff_09",&errmigrad_banff_09,"errmigrad_banff_09/D");
    nfitnt->Branch("errmigrad_banff_10",&errmigrad_banff_10,"errmigrad_banff_10/D");
    nfitnt->Branch("errmigrad_banff_11",&errmigrad_banff_11,"errmigrad_banff_11/D");
    nfitnt->Branch("errmigrad_banff_12",&errmigrad_banff_12,"errmigrad_banff_12/D");
    nfitnt->Branch("errmigrad_banff_13",&errmigrad_banff_13,"errmigrad_banff_13/D");
    nfitnt->Branch("errmigrad_banff_14",&errmigrad_banff_14,"errmigrad_banff_14/D");
    nfitnt->Branch("errmigrad_banff_15",&errmigrad_banff_15,"errmigrad_banff_15/D");
    nfitnt->Branch("errmigrad_banff_maqe",&errmigrad_banff_maqe,"errmigrad_banff_maqe/D");
    nfitnt->Branch("errmigrad_banff_mares",&errmigrad_banff_mares,"errmigrad_banff_mares/D");
    nfitnt->Branch("errmigrad_banff_ccqeE0",&errmigrad_banff_ccqeE0,"errmigrad_banff_ccqeE0/D");
    nfitnt->Branch("errmigrad_banff_ccqeE1",&errmigrad_banff_ccqeE1,"errmigrad_banff_ccqeE1/D");
    nfitnt->Branch("errmigrad_banff_ccqeE2",&errmigrad_banff_ccqeE2,"errmigrad_banff_ccqeE2/D");
    nfitnt->Branch("errmigrad_banff_cc1piE0",&errmigrad_banff_cc1piE0,"errmigrad_banff_cc1piE0/D");
    nfitnt->Branch("errmigrad_banff_cc1piE1",&errmigrad_banff_cc1piE1,"errmigrad_banff_cc1piE1/D");
  }
  if(fit_syst_mode >= 1 && fit_param_list.find("unbanff") != string::npos) {
    nfitnt->Branch("errmigrad_unbanff_00",&errmigrad_unbanff_00,"errmigrad_unbanff_00/D");
    nfitnt->Branch("errmigrad_unbanff_01",&errmigrad_unbanff_01,"errmigrad_unbanff_01/D");
    nfitnt->Branch("errmigrad_unbanff_02",&errmigrad_unbanff_02,"errmigrad_unbanff_02/D");
    nfitnt->Branch("errmigrad_unbanff_03",&errmigrad_unbanff_03,"errmigrad_unbanff_03/D");
    nfitnt->Branch("errmigrad_unbanff_04",&errmigrad_unbanff_04,"errmigrad_unbanff_04/D");
    nfitnt->Branch("errmigrad_unbanff_05",&errmigrad_unbanff_05,"errmigrad_unbanff_05/D");
    nfitnt->Branch("errmigrad_unbanff_06",&errmigrad_unbanff_06,"errmigrad_unbanff_06/D");
    nfitnt->Branch("errmigrad_unbanff_07",&errmigrad_unbanff_07,"errmigrad_unbanff_07/D");
    nfitnt->Branch("errmigrad_unbanff_08",&errmigrad_unbanff_08,"errmigrad_unbanff_08/D");
    nfitnt->Branch("errmigrad_unbanff_09",&errmigrad_unbanff_09,"errmigrad_unbanff_09/D");
    nfitnt->Branch("errmigrad_unbanff_10",&errmigrad_unbanff_10,"errmigrad_unbanff_10/D");
    nfitnt->Branch("errmigrad_unbanff_11",&errmigrad_unbanff_11,"errmigrad_unbanff_11/D");
    nfitnt->Branch("errmigrad_unbanff_12",&errmigrad_unbanff_12,"errmigrad_unbanff_12/D");
    nfitnt->Branch("errmigrad_unbanff_13",&errmigrad_unbanff_13,"errmigrad_unbanff_13/D");
    nfitnt->Branch("errmigrad_unbanff_14",&errmigrad_unbanff_14,"errmigrad_unbanff_14/D");
    nfitnt->Branch("errmigrad_unbanff_15",&errmigrad_unbanff_15,"errmigrad_unbanff_15/D");
    nfitnt->Branch("errmigrad_unbanff_16",&errmigrad_unbanff_16,"errmigrad_unbanff_16/D");
    nfitnt->Branch("errmigrad_unbanff_17",&errmigrad_unbanff_17,"errmigrad_unbanff_17/D");
    nfitnt->Branch("errmigrad_unbanff_18",&errmigrad_unbanff_18,"errmigrad_unbanff_18/D");
    nfitnt->Branch("errmigrad_unbanff_19",&errmigrad_unbanff_19,"errmigrad_unbanff_19/D");
    nfitnt->Branch("errmigrad_unbanff_20",&errmigrad_unbanff_20,"errmigrad_unbanff_20/D");
    nfitnt->Branch("errmigrad_unbanff_21",&errmigrad_unbanff_21,"errmigrad_unbanff_21/D");
    nfitnt->Branch("errmigrad_unbanff_22",&errmigrad_unbanff_22,"errmigrad_unbanff_22/D");
  }
  if(fit_syst_mode == 2) {
    nfitnt->Branch("errmigrad_xsec_pf",&errmigrad_xsec_pf,"errmigrad_xsec_pf/D");
    nfitnt->Branch("errmigrad_xsec_wshape",&errmigrad_xsec_wshape,"errmigrad_xsec_wshape/D");
    nfitnt->Branch("errmigrad_xsec_ccothshape",&errmigrad_xsec_ccothshape,"errmigrad_xsec_ccothshape/D");
    nfitnt->Branch("errmigrad_xsec_cc1pishape",&errmigrad_xsec_cc1pishape,"errmigrad_xsec_cc1pishape/D");
    nfitnt->Branch("errmigrad_xsec_normcccoh",&errmigrad_xsec_normcccoh,"errmigrad_xsec_normcccoh/D");
    nfitnt->Branch("errmigrad_xsec_normnc1piC",&errmigrad_xsec_normnc1piC,"errmigrad_xsec_normnc1piC/D");
    nfitnt->Branch("errmigrad_xsec_normnccoh",&errmigrad_xsec_normnccoh,"errmigrad_xsec_normnccoh/D");
    nfitnt->Branch("errmigrad_xsec_normncoth",&errmigrad_xsec_normncoth,"errmigrad_xsec_normncoth/D");
    nfitnt->Branch("errmigrad_xsec_normnuetonumu",&errmigrad_xsec_normnuetonumu,"errmigrad_xsec_normnuetonumu/D");
    nfitnt->Branch("errmigrad_xsec_normantinutonu",&errmigrad_xsec_normantinutonu,"errmigrad_xsec_normantinutonu/D");
    nfitnt->Branch("errmigrad_skdet_ccqerc0",&errmigrad_skdet_ccqerc0,"errmigrad_skdet_ccqerc0/D");
    nfitnt->Branch("errmigrad_skdet_ccqerc1",&errmigrad_skdet_ccqerc1,"errmigrad_skdet_ccqerc1/D");
    nfitnt->Branch("errmigrad_skdet_ccqerc2",&errmigrad_skdet_ccqerc2,"errmigrad_skdet_ccqerc2/D");
    nfitnt->Branch("errmigrad_skdet_ccqeoth",&errmigrad_skdet_ccqeoth,"errmigrad_skdet_ccqeoth/D");
    nfitnt->Branch("errmigrad_skdet_nuecc",&errmigrad_skdet_nuecc,"errmigrad_skdet_nuecc/D");
    nfitnt->Branch("errmigrad_fsi_000",&errmigrad_fsi_000,"errmigrad_fsi_000/D");
    nfitnt->Branch("errmigrad_fsi_001",&errmigrad_fsi_001,"errmigrad_fsi_001/D");
    nfitnt->Branch("errmigrad_fsi_002",&errmigrad_fsi_002,"errmigrad_fsi_002/D");
    nfitnt->Branch("errmigrad_fsi_003",&errmigrad_fsi_003,"errmigrad_fsi_003/D");
    nfitnt->Branch("errmigrad_fsi_004",&errmigrad_fsi_004,"errmigrad_fsi_004/D");
    nfitnt->Branch("errmigrad_fsi_005",&errmigrad_fsi_005,"errmigrad_fsi_005/D");
  }
  //errhesses
  nfitnt->Branch("errhesse_mixing12",&errhesse_mixing12,"errhesse_mixing12/D");
  nfitnt->Branch("errhesse_mixing13",&errhesse_mixing13,"errhesse_mixing13/D");
  nfitnt->Branch("errhesse_mixing23",&errhesse_mixing23,"errhesse_mixing23/D");
  nfitnt->Branch("errhesse_deltacp",&errhesse_deltacp,"errhesse_deltacp/D");
  nfitnt->Branch("errhesse_dm12sq",&errhesse_dm12sq,"errhesse_dm12sq/D");
  nfitnt->Branch("errhesse_dm23sq",&errhesse_dm23sq,"errhesse_dm23sq/D");
  if(fit_syst_mode >= 1) {
    nfitnt->Branch("errhesse_xsec_sf",&errhesse_xsec_sf,"errhesse_xsec_sf/D");
    nfitnt->Branch("errhesse_xsec_binde",&errhesse_xsec_binde,"errhesse_xsec_binde/D");
    nfitnt->Branch("errhesse_xsec_pilessdelta",&errhesse_xsec_pilessdelta,"errhesse_xsec_pilessdelta/D");
    nfitnt->Branch("errhesse_skdet_ccnqe",&errhesse_skdet_ccnqe,"errhesse_skdet_ccnqe/D");
    nfitnt->Branch("errhesse_skdet_nc",&errhesse_skdet_nc,"errhesse_skdet_nc/D");
    nfitnt->Branch("errhesse_skdet_energyscale",&errhesse_skdet_energyscale,"errhesse_skdet_energyscale/D");
    nfitnt->Branch("errhesse_banff_00",&errhesse_banff_00,"errhesse_banff_00/D");
    nfitnt->Branch("errhesse_banff_01",&errhesse_banff_01,"errhesse_banff_01/D");
    nfitnt->Branch("errhesse_banff_02",&errhesse_banff_02,"errhesse_banff_02/D");
    nfitnt->Branch("errhesse_banff_03",&errhesse_banff_03,"errhesse_banff_03/D");
    nfitnt->Branch("errhesse_banff_04",&errhesse_banff_04,"errhesse_banff_04/D");
    nfitnt->Branch("errhesse_banff_05",&errhesse_banff_05,"errhesse_banff_05/D");
    nfitnt->Branch("errhesse_banff_06",&errhesse_banff_06,"errhesse_banff_06/D");
    nfitnt->Branch("errhesse_banff_07",&errhesse_banff_07,"errhesse_banff_07/D");
    nfitnt->Branch("errhesse_banff_08",&errhesse_banff_08,"errhesse_banff_08/D");
    nfitnt->Branch("errhesse_banff_09",&errhesse_banff_09,"errhesse_banff_09/D");
    nfitnt->Branch("errhesse_banff_10",&errhesse_banff_10,"errhesse_banff_10/D");
    nfitnt->Branch("errhesse_banff_11",&errhesse_banff_11,"errhesse_banff_11/D");
    nfitnt->Branch("errhesse_banff_12",&errhesse_banff_12,"errhesse_banff_12/D");
    nfitnt->Branch("errhesse_banff_13",&errhesse_banff_13,"errhesse_banff_13/D");
    nfitnt->Branch("errhesse_banff_14",&errhesse_banff_14,"errhesse_banff_14/D");
    nfitnt->Branch("errhesse_banff_15",&errhesse_banff_15,"errhesse_banff_15/D");
    nfitnt->Branch("errhesse_banff_maqe",&errhesse_banff_maqe,"errhesse_banff_maqe/D");
    nfitnt->Branch("errhesse_banff_mares",&errhesse_banff_mares,"errhesse_banff_mares/D");
    nfitnt->Branch("errhesse_banff_ccqeE0",&errhesse_banff_ccqeE0,"errhesse_banff_ccqeE0/D");
    nfitnt->Branch("errhesse_banff_ccqeE1",&errhesse_banff_ccqeE1,"errhesse_banff_ccqeE1/D");
    nfitnt->Branch("errhesse_banff_ccqeE2",&errhesse_banff_ccqeE2,"errhesse_banff_ccqeE2/D");
    nfitnt->Branch("errhesse_banff_cc1piE0",&errhesse_banff_cc1piE0,"errhesse_banff_cc1piE0/D");
    nfitnt->Branch("errhesse_banff_cc1piE1",&errhesse_banff_cc1piE1,"errhesse_banff_cc1piE1/D");
  }
  if(fit_syst_mode >= 1 && fit_param_list.find("unbanff") != string::npos) {
    nfitnt->Branch("errhesse_unbanff_00",&errhesse_unbanff_00,"errhesse_unbanff_00/D");
    nfitnt->Branch("errhesse_unbanff_01",&errhesse_unbanff_01,"errhesse_unbanff_01/D");
    nfitnt->Branch("errhesse_unbanff_02",&errhesse_unbanff_02,"errhesse_unbanff_02/D");
    nfitnt->Branch("errhesse_unbanff_03",&errhesse_unbanff_03,"errhesse_unbanff_03/D");
    nfitnt->Branch("errhesse_unbanff_04",&errhesse_unbanff_04,"errhesse_unbanff_04/D");
    nfitnt->Branch("errhesse_unbanff_05",&errhesse_unbanff_05,"errhesse_unbanff_05/D");
    nfitnt->Branch("errhesse_unbanff_06",&errhesse_unbanff_06,"errhesse_unbanff_06/D");
    nfitnt->Branch("errhesse_unbanff_07",&errhesse_unbanff_07,"errhesse_unbanff_07/D");
    nfitnt->Branch("errhesse_unbanff_08",&errhesse_unbanff_08,"errhesse_unbanff_08/D");
    nfitnt->Branch("errhesse_unbanff_09",&errhesse_unbanff_09,"errhesse_unbanff_09/D");
    nfitnt->Branch("errhesse_unbanff_10",&errhesse_unbanff_10,"errhesse_unbanff_10/D");
    nfitnt->Branch("errhesse_unbanff_11",&errhesse_unbanff_11,"errhesse_unbanff_11/D");
    nfitnt->Branch("errhesse_unbanff_12",&errhesse_unbanff_12,"errhesse_unbanff_12/D");
    nfitnt->Branch("errhesse_unbanff_13",&errhesse_unbanff_13,"errhesse_unbanff_13/D");
    nfitnt->Branch("errhesse_unbanff_14",&errhesse_unbanff_14,"errhesse_unbanff_14/D");
    nfitnt->Branch("errhesse_unbanff_15",&errhesse_unbanff_15,"errhesse_unbanff_15/D");
    nfitnt->Branch("errhesse_unbanff_16",&errhesse_unbanff_16,"errhesse_unbanff_16/D");
    nfitnt->Branch("errhesse_unbanff_17",&errhesse_unbanff_17,"errhesse_unbanff_17/D");
    nfitnt->Branch("errhesse_unbanff_18",&errhesse_unbanff_18,"errhesse_unbanff_18/D");
    nfitnt->Branch("errhesse_unbanff_19",&errhesse_unbanff_19,"errhesse_unbanff_19/D");
    nfitnt->Branch("errhesse_unbanff_20",&errhesse_unbanff_20,"errhesse_unbanff_20/D");
    nfitnt->Branch("errhesse_unbanff_21",&errhesse_unbanff_21,"errhesse_unbanff_21/D");
    nfitnt->Branch("errhesse_unbanff_22",&errhesse_unbanff_22,"errhesse_unbanff_22/D");
  }
  if(fit_syst_mode == 2) {
    nfitnt->Branch("errhesse_xsec_pf",&errhesse_xsec_pf,"errhesse_xsec_pf/D");
    nfitnt->Branch("errhesse_xsec_wshape",&errhesse_xsec_wshape,"errhesse_xsec_wshape/D");
    nfitnt->Branch("errhesse_xsec_ccothshape",&errhesse_xsec_ccothshape,"errhesse_xsec_ccothshape/D");
    nfitnt->Branch("errhesse_xsec_cc1pishape",&errhesse_xsec_cc1pishape,"errhesse_xsec_cc1pishape/D");
    nfitnt->Branch("errhesse_xsec_normcccoh",&errhesse_xsec_normcccoh,"errhesse_xsec_normcccoh/D");
    nfitnt->Branch("errhesse_xsec_normnc1piC",&errhesse_xsec_normnc1piC,"errhesse_xsec_normnc1piC/D");
    nfitnt->Branch("errhesse_xsec_normnccoh",&errhesse_xsec_normnccoh,"errhesse_xsec_normnccoh/D");
    nfitnt->Branch("errhesse_xsec_normncoth",&errhesse_xsec_normncoth,"errhesse_xsec_normncoth/D");
    nfitnt->Branch("errhesse_xsec_normnuetonumu",&errhesse_xsec_normnuetonumu,"errhesse_xsec_normnuetonumu/D");
    nfitnt->Branch("errhesse_xsec_normantinutonu",&errhesse_xsec_normantinutonu,"errhesse_xsec_normantinutonu/D");
    nfitnt->Branch("errhesse_skdet_ccqerc0",&errhesse_skdet_ccqerc0,"errhesse_skdet_ccqerc0/D");
    nfitnt->Branch("errhesse_skdet_ccqerc1",&errhesse_skdet_ccqerc1,"errhesse_skdet_ccqerc1/D");
    nfitnt->Branch("errhesse_skdet_ccqerc2",&errhesse_skdet_ccqerc2,"errhesse_skdet_ccqerc2/D");
    nfitnt->Branch("errhesse_skdet_ccqeoth",&errhesse_skdet_ccqeoth,"errhesse_skdet_ccqeoth/D");
    nfitnt->Branch("errhesse_skdet_nuecc",&errhesse_skdet_nuecc,"errhesse_skdet_nuecc/D");
    nfitnt->Branch("errhesse_fsi_000",&errhesse_fsi_000,"errhesse_fsi_000/D");
    nfitnt->Branch("errhesse_fsi_001",&errhesse_fsi_001,"errhesse_fsi_001/D");
    nfitnt->Branch("errhesse_fsi_002",&errhesse_fsi_002,"errhesse_fsi_002/D");
    nfitnt->Branch("errhesse_fsi_003",&errhesse_fsi_003,"errhesse_fsi_003/D");
    nfitnt->Branch("errhesse_fsi_004",&errhesse_fsi_004,"errhesse_fsi_004/D");
    nfitnt->Branch("errhesse_fsi_005",&errhesse_fsi_005,"errhesse_fsi_005/D");
  }


  ///////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////
  // Loop over toy-MC experiments
  ///////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////
  for(int iexperiment = 0; iexperiment < nexperiments; iexperiment++) {

    cout << "\n\n *** Pseudo-experiment = " << iexperiment << endl;
    iexp = iexperiment;

    // randomize specified osc params for current toy MC expt
    if (randomize_mixing12) { input_mixing12 = util::RandomizeMixing12(); }
    if (randomize_mixing13) { input_mixing13 = util::RandomizeMixing13(); }
    if (randomize_deltacp ) { input_deltacp  = util::RandomizeDeltaCP();  }
    if (randomize_dm12sq  ) { input_dm12sq   = util::RandomizeDm12sq();   }

    // set oscillation parameter inputs as necessary (used to generate toy-MC expt)
    if(hyp == kHypOsc3F) {
      fparm->Update(kPhypMixing12, input_mixing12);
      fparm->Update(kPhypMixing13, input_mixing13);
      fparm->Update(kPhypMixing23, input_mixing23);
      fparm->Update(kPhypDeltaCP,  input_deltacp);
      fparm->Update(kPhypDm12sq,   input_dm12sq);
      fparm->Update(kPhypDm23sq,   input_dm23sq);
    }
    else
      if(hyp == kHypNuMuDisapp2F) {
	fparm->Update(kPhypMixing23, input_mixing23);
	fparm->Update(kPhypDm23sq,   input_dm23sq);
      }

    //
    // get superk spectrum observation (toy experiment)
    //
    SuperkObservation * skobs = new SuperkObservation;
    bool applied_negevents_condition = true;
    while(applied_negevents_condition) {
      if(syst_mode > 0) {
	// tweak all sources of systematic err (gaussian err)
	cout << "- Randomizing systematic params in pdf..." << endl;
	switch(syst_mode) {
	case 1:
	  skexp->TweakAllSystematics();
	  break;
	case 2:
	  skexp->TweakAllBANFFSystematics();
	  break;
	case 3:
	  skexp->TweakAllUncorrelatedXSecSystematics();
	  break;
	case 4:
	  skexp->TweakAllSKDetectorSystematics();
	  break;
	case 5:
	  skexp->TweakAllFSISystematics();
	  break;
	case 6:
	  TRandom3 rand(0);
	  skexp->GetNuXSec()->TweakSpecFunc      (rand.Uniform(0.,1.));
	  skexp->GetNuXSec()->TweakBindE         (rand.Gaus(0.,1.));
	  skexp->GetNuXSec()->TweakPiLessDelta   (rand.Gaus(0.,1.));
	  skexp->GetSuperkDetector()->TweakEfficiency1RmuCCnonQE    (rand.Gaus(0.,1.));
	  skexp->GetSuperkDetector()->TweakEfficiency1RmuNC         (rand.Gaus(0.,1.));
	  //skexp->GetSuperkDetector()->TweakEnergyScale              (rand.Gaus(0.,1.));
	  TVectorD twkv(skexp->GetBANFF()->BestFitValues());
	  TVectorD cov (skexp->GetBANFF()->CovMatrix());
	  PhysParamId_t p_maqe = conv::kNuisBANFFMAQE;
	  int ip_maqe    = util::ArrIdx::Instance()->BANFFSystIdx(kNuisBANFFMAQE);
	  twkv[ip_maqe]  = rand.Gaus(twkv[ip_maqe], cov[ip_maqe]);
	  double sigma_maqe = (twkv[ip_maqe] - skexp->GetBANFF()->BestFitValue(p_maqe)) / cov[ip_maqe][ip_maqe];
	  PhysParamId_t p_ccqe0 = kNuisBANFFCCQEE0;
	  int ip_ccqe0       = util::ArrIdx::Instance()->BANFFSystIdx(kNuisBANFFCCQEE0);
	  twkv[ip_ccqe0] = rand.Gaus(twkv[ip_ccqe0], cov[ip_ccqe0]);
	  int ip_ccqe0 = util::ArrIdx::Instance()->BANFFSystIdx(p_ccqe0);
	  double sigma_ccqe0 = (twkv[ip_ccqe0] - skexp->GetBANFF()->BestFitValue(p_ccqe0)) / cov[ip_ccqe0][ip_ccqe0];
	  // now, tweak the other parameters according to the correlations between them and the fitted parameters
	  for(int ip = 0; ip < nbanff; ip++) {
	    int p_i = ip + (int)conv::kNuisBANFF00;
	    conv::PhysParamId_t p = (conv::PhysParamId_t)p_i;
	    if(p == kNuisBANFFMAQE || p == kNuisBANFFCCQEE0)
	      continue;
	    twkv[ip] = skexp->GetBANFF()->BestFitValue(p);
	    twkv[ip] += sigma_maqe  * cov[ip_maqe] [ip];
	    twkv[ip] += sigma_ccqe0 * cov[ip_ccqe0][ip];
	  }
	  skexp->GetBANFF()->Tweak(true,&twkv);
	  break;
	default:
	  cerr << "Invalid systematic mode " << syst_mode << ". Exiting" << endl;
	  exit(1);
	  break;
	}
	skexp->PrintSystematicVariations();
	// save the input systematics to the tree
	input_banff_00 = skexp->GetBANFF()->CurrentValue(conv::kNuisBANFF00);
	input_banff_01 = skexp->GetBANFF()->CurrentValue(conv::kNuisBANFF01);
	input_banff_02 = skexp->GetBANFF()->CurrentValue(conv::kNuisBANFF02);
	input_banff_03 = skexp->GetBANFF()->CurrentValue(conv::kNuisBANFF03);
	input_banff_04 = skexp->GetBANFF()->CurrentValue(conv::kNuisBANFF04);
	input_banff_05 = skexp->GetBANFF()->CurrentValue(conv::kNuisBANFF05);
	input_banff_06 = skexp->GetBANFF()->CurrentValue(conv::kNuisBANFF06);
	input_banff_07 = skexp->GetBANFF()->CurrentValue(conv::kNuisBANFF07);
	input_banff_08 = skexp->GetBANFF()->CurrentValue(conv::kNuisBANFF08);
	input_banff_09 = skexp->GetBANFF()->CurrentValue(conv::kNuisBANFF09);
	input_banff_10 = skexp->GetBANFF()->CurrentValue(conv::kNuisBANFF10);
	input_banff_11 = skexp->GetBANFF()->CurrentValue(conv::kNuisBANFF11);
	input_banff_12 = skexp->GetBANFF()->CurrentValue(conv::kNuisBANFF12);
	input_banff_13 = skexp->GetBANFF()->CurrentValue(conv::kNuisBANFF13);
	input_banff_14 = skexp->GetBANFF()->CurrentValue(conv::kNuisBANFF14);
	input_banff_15 = skexp->GetBANFF()->CurrentValue(conv::kNuisBANFF15);
	input_banff_maqe = skexp->GetBANFF()->CurrentValue(conv::kNuisBANFFMAQE);
	input_banff_mares = skexp->GetBANFF()->CurrentValue(conv::kNuisBANFFMARES);
	input_banff_ccqeE0 = skexp->GetBANFF()->CurrentValue(conv::kNuisBANFFCCQEE0);
	input_banff_ccqeE1 = skexp->GetBANFF()->CurrentValue(conv::kNuisBANFFCCQEE1);
	input_banff_ccqeE2 = skexp->GetBANFF()->CurrentValue(conv::kNuisBANFFCCQEE2);
	input_banff_cc1piE0 = skexp->GetBANFF()->CurrentValue(conv::kNuisBANFFCC1piE0);
	input_banff_cc1piE1 = skexp->GetBANFF()->CurrentValue(conv::kNuisBANFFCC1piE1);
	input_unbanff_00 = skexp->GetBANFF()->CurrentUncorrelatedValue(conv::kNuisUnBANFF00);
	input_unbanff_01 = skexp->GetBANFF()->CurrentUncorrelatedValue(conv::kNuisUnBANFF01);
	input_unbanff_02 = skexp->GetBANFF()->CurrentUncorrelatedValue(conv::kNuisUnBANFF02);
	input_unbanff_03 = skexp->GetBANFF()->CurrentUncorrelatedValue(conv::kNuisUnBANFF03);
	input_unbanff_04 = skexp->GetBANFF()->CurrentUncorrelatedValue(conv::kNuisUnBANFF04);
	input_unbanff_05 = skexp->GetBANFF()->CurrentUncorrelatedValue(conv::kNuisUnBANFF05);
	input_unbanff_06 = skexp->GetBANFF()->CurrentUncorrelatedValue(conv::kNuisUnBANFF06);
	input_unbanff_07 = skexp->GetBANFF()->CurrentUncorrelatedValue(conv::kNuisUnBANFF07);
	input_unbanff_08 = skexp->GetBANFF()->CurrentUncorrelatedValue(conv::kNuisUnBANFF08);
	input_unbanff_09 = skexp->GetBANFF()->CurrentUncorrelatedValue(conv::kNuisUnBANFF09);
	input_unbanff_10 = skexp->GetBANFF()->CurrentUncorrelatedValue(conv::kNuisUnBANFF10);
	input_unbanff_11 = skexp->GetBANFF()->CurrentUncorrelatedValue(conv::kNuisUnBANFF11);
	input_unbanff_12 = skexp->GetBANFF()->CurrentUncorrelatedValue(conv::kNuisUnBANFF12);
	input_unbanff_13 = skexp->GetBANFF()->CurrentUncorrelatedValue(conv::kNuisUnBANFF13);
	input_unbanff_14 = skexp->GetBANFF()->CurrentUncorrelatedValue(conv::kNuisUnBANFF14);
	input_unbanff_15 = skexp->GetBANFF()->CurrentUncorrelatedValue(conv::kNuisUnBANFF15);
	input_unbanff_16 = skexp->GetBANFF()->CurrentUncorrelatedValue(conv::kNuisUnBANFF16);
	input_unbanff_17 = skexp->GetBANFF()->CurrentUncorrelatedValue(conv::kNuisUnBANFF17);
	input_unbanff_18 = skexp->GetBANFF()->CurrentUncorrelatedValue(conv::kNuisUnBANFF18);
	input_unbanff_19 = skexp->GetBANFF()->CurrentUncorrelatedValue(conv::kNuisUnBANFF19);
	input_unbanff_20 = skexp->GetBANFF()->CurrentUncorrelatedValue(conv::kNuisUnBANFF20);
	input_unbanff_21 = skexp->GetBANFF()->CurrentUncorrelatedValue(conv::kNuisUnBANFF21);
	input_unbanff_22 = skexp->GetBANFF()->CurrentUncorrelatedValue(conv::kNuisUnBANFF22);
	input_xsec_pf = skexp->GetNuXSec()->CurrentValue(conv::kNuisXSecFermiMom);
	input_xsec_sf = skexp->GetNuXSec()->CurrentValue(conv::kNuisXSecSpecFunc);
	input_xsec_binde = skexp->GetNuXSec()->CurrentValue(conv::kNuisXSecBindE);
	input_xsec_wshape = skexp->GetNuXSec()->CurrentValue(conv::kNuisXSecWshape);
	input_xsec_ccothshape = skexp->GetNuXSec()->CurrentValue(conv::kNuisXSecShapeCCoth);
	input_xsec_cc1pishape = skexp->GetNuXSec()->CurrentValue(conv::kNuisXSecENuShape);
	input_xsec_pilessdelta = skexp->GetNuXSec()->CurrentValue(conv::kNuisXSecPiLessDelta);
	input_xsec_normcccoh = skexp->GetNuXSec()->CurrentValue(conv::kNuisXSecNormCCcoh);
	input_xsec_normnc1piC = skexp->GetNuXSec()->CurrentValue(conv::kNuisXSecNormNC1piC);
	input_xsec_normnccoh = skexp->GetNuXSec()->CurrentValue(conv::kNuisXSecNormNCcoh);
	input_xsec_normncoth = skexp->GetNuXSec()->CurrentValue(conv::kNuisXSecNormNCoth);
	input_xsec_normnuetonumu = skexp->GetNuXSec()->CurrentValue(conv::kNuisXSecNormNuEToNuMu);
	input_xsec_normantinutonu = skexp->GetNuXSec()->CurrentValue(conv::kNuisXSecNormAntiNuToNu);
	input_skdet_energyscale = skexp->GetSuperkDetector()->CurrentValue(conv::kNuisSuperkEnergyScale);
	input_skdet_ccqerc0 = skexp->GetSuperkDetector()->CurrentValue(conv::kNuisSuperkEff1RmuCCQErc0u);
	input_skdet_ccqerc1 = skexp->GetSuperkDetector()->CurrentValue(conv::kNuisSuperkEff1RmuCCQErc1u);
	input_skdet_ccqerc2 = skexp->GetSuperkDetector()->CurrentValue(conv::kNuisSuperkEff1RmuCCQErc2u);
	input_skdet_ccqeoth = skexp->GetSuperkDetector()->CurrentValue(conv::kNuisSuperkEff1RmuCCQEoth);
	input_skdet_ccnqe = skexp->GetSuperkDetector()->CurrentValue(conv::kNuisSuperkEff1RmuCCnonQE);
	input_skdet_nc = skexp->GetSuperkDetector()->CurrentValue(conv::kNuisSuperkEff1RmuNC);
	input_skdet_nuecc = skexp->GetSuperkDetector()->CurrentValue(conv::kNuisSuperkEff1RmuNueCC);
	input_fsi_000 = skexp->GetFSI()->CurrentValue(conv::kNuisFSI000);
	input_fsi_001 = skexp->GetFSI()->CurrentValue(conv::kNuisFSI001);
	input_fsi_002 = skexp->GetFSI()->CurrentValue(conv::kNuisFSI002);
	input_fsi_003 = skexp->GetFSI()->CurrentValue(conv::kNuisFSI003);
	input_fsi_004 = skexp->GetFSI()->CurrentValue(conv::kNuisFSI004);
	input_fsi_005 = skexp->GetFSI()->CurrentValue(conv::kNuisFSI005);
	// generate mock data using the current PDF (incl. systematic variations) 
	// + statistical fluctuations if requested (Poisson fluctuation in each bin)
	cout << "- Generating toy-MC..." << endl;
	applied_negevents_condition = skobs->CreateMockData(*skexp, PhysHypothesis(*fparm), stat_fluctuations);
	// reset all systematics so that the PDF used in the fitter is the nominal one
	cout << "- Reseting all systematic params in pdf..." << endl;
	skexp->ResetSystematicVariations();
	//skexp->PrintSystematicVariations();
      } //systmode >0
      else {       
	// generate mock data using the current PDF (without systematic variations) 
	// with statistical fluctuations if requested (Poisson fluctuation in each bin)
	skexp->ResetSystematicVariations(); //make sure there are no systematic variations
	cout << "- Generating toy-MC..." << endl;
	applied_negevents_condition = skobs->CreateMockData(*skexp, PhysHypothesis(*fparm), stat_fluctuations); 
      }
    }//applied_negevents_condition

    // Get nominal number of events
    nsk_nom = skexp->CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike)->Integral();

    //Seed the FitParamList with true values
    if(fit_param_list == "nuosc3f_numudisapp_fix4+nuisanceselection") {
      fparm->SetInitValue(kPhypMixing12, input_mixing12);
      fparm->SetInitValue(kPhypMixing13, input_mixing13);
      if(input_mixing23 > 0.99)
	fparm->SetInitValue(kPhypMixing23, 0.95);
      else
	fparm->SetInitValue(kPhypMixing23, input_mixing23);
      fparm->SetInitValue(kPhypDeltaCP,  input_deltacp);
      fparm->SetInitValue(kPhypDm12sq,   input_dm12sq);
      fparm->SetInitValue(kPhypDm23sq,   input_dm23sq);
      fparm->SetInitValue(conv::kNuisBANFFMAQE,            input_banff_maqe);
      fparm->SetInitValue(conv::kNuisBANFFCCQEE0,          input_banff_ccqeE0);
      fparm->SetInitValue(conv::kNuisXSecSpecFunc,         input_xsec_sf);
      fparm->SetInitValue(conv::kNuisXSecBindE,            input_xsec_binde);
      fparm->SetInitValue(conv::kNuisXSecPiLessDelta,      input_xsec_pilessdelta);
      //fparm->SetInitValue(conv::kNuisSuperkEnergyScale,    input_skdet_energyscale);
      fparm->SetInitValue(conv::kNuisSuperkEff1RmuCCnonQE, input_skdet_ccnqe);
      fparm->SetInitValue(conv::kNuisSuperkEff1RmuNC,      input_skdet_nc);
    }

    // Get fitter
    FitInputs * fit_inputs = new FitInputs(skobs,skexp,fparm);
    FitterABC * fitter = 0;
    if(use_minuit2)
      fitter = new Superk1RmuErecoFitterMinuit2(fit_inputs,fit_option);
    else
      fitter = new Superk1RmuErecoFitter(fit_inputs,fit_option);

    // Chisq for true (input) params
    true_chisq = fit_inputs->Chisquare();

    // Print-out true params
    //fparm->ShowOnlyCurrValuesInNextPrint();
    //cout << "- True oscillation params (below) yield chisq = " << true_chisq
    //	   << *fparm;

    // Do the fit
    fparm->ResetCurrToInitValues(); // reset osc params
    fitter->DoTheFit();
    fitstatus = fitter->FitStatusCode();
    bestfit_chisq = fit_inputs->BestFitChisquare();

    // Get bestfit number of events
    nsk_bf = skexp->CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike)->Integral();

    // Print-out best-fit params
    fparm->ShowOnlyBestFitValuesInNextPrint();
    cout << "- Best-fit oscillation params (below) yield chisq = " << bestfit_chisq 
	 << *fparm;

    ////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////
    // Get & store best-fit params
    ////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////
    // Init best-fit values
    bestfit_mixing12 = -1;
    bestfit_mixing13 = -1;
    bestfit_mixing23 = -1;
    bestfit_deltacp  = -1;
    bestfit_dm12sq   = -1;
    bestfit_dm23sq   = -1;
    bestfit_unbanff_00 = -1;
    bestfit_unbanff_01 = -1;
    bestfit_unbanff_02 = -1;
    bestfit_unbanff_03 = -1;
    bestfit_unbanff_04 = -1;
    bestfit_unbanff_05 = -1;
    bestfit_unbanff_06 = -1;
    bestfit_unbanff_07 = -1;
    bestfit_unbanff_08 = -1;
    bestfit_unbanff_09 = -1;
    bestfit_unbanff_10 = -1;
    bestfit_unbanff_11 = -1;
    bestfit_unbanff_12 = -1;
    bestfit_unbanff_13 = -1;
    bestfit_unbanff_14 = -1;
    bestfit_unbanff_15 = -1;
    bestfit_unbanff_16 = -1;
    bestfit_unbanff_17 = -1;
    bestfit_unbanff_18 = -1;
    bestfit_unbanff_19 = -1;
    bestfit_unbanff_20 = -1;
    bestfit_unbanff_21 = -1;
    bestfit_unbanff_22 = -1;
    bestfit_banff_00 = -1;
    bestfit_banff_01 = -1;
    bestfit_banff_02 = -1;
    bestfit_banff_03 = -1;
    bestfit_banff_04 = -1;
    bestfit_banff_05 = -1;
    bestfit_banff_06 = -1;
    bestfit_banff_07 = -1;
    bestfit_banff_08 = -1;
    bestfit_banff_09 = -1;
    bestfit_banff_10 = -1;
    bestfit_banff_11 = -1;
    bestfit_banff_12 = -1;
    bestfit_banff_13 = -1;
    bestfit_banff_14 = -1;
    bestfit_banff_15 = -1;
    bestfit_banff_maqe = -1;
    bestfit_banff_mares = -1;
    bestfit_banff_ccqeE0 = -1;
    bestfit_banff_ccqeE1 = -1;
    bestfit_banff_ccqeE2 = -1;
    bestfit_banff_cc1piE0 = -1;
    bestfit_banff_cc1piE1 = -1;
    bestfit_xsec_pf = -1;
    bestfit_xsec_sf = -1;
    bestfit_xsec_wshape = -1;
    bestfit_xsec_ccothshape = -1;
    bestfit_xsec_cc1pishape = -1;
    bestfit_xsec_pilessdelta = -1;
    bestfit_xsec_normcccoh = -1;
    bestfit_xsec_normnc1piC = -1;
    bestfit_xsec_normnccoh = -1;
    bestfit_xsec_normncoth = -1;
    bestfit_xsec_normnuetonumu = -1;
    bestfit_xsec_normantinutonu = -1;
    bestfit_skdet_energyscale = -1;
    bestfit_skdet_ccqerc0 = -1;
    bestfit_skdet_ccqerc1 = -1;
    bestfit_skdet_ccqerc2 = -1;
    bestfit_skdet_ccqeoth = -1;
    bestfit_skdet_ccnqe = -1;
    bestfit_skdet_nc = -1;
    bestfit_skdet_nuecc = -1;
    bestfit_fsi_000 = -1;
    bestfit_fsi_001 = -1;
    bestfit_fsi_002 = -1;
    bestfit_fsi_003 = -1;
    bestfit_fsi_004 = -1;
    bestfit_fsi_005 = -1;
    // Init plus errors
    errminosplus_mixing12 = -1;
    errminosplus_mixing13 = -1;
    errminosplus_mixing23 = -1;
    errminosplus_deltacp  = -1;
    errminosplus_dm12sq   = -1;
    errminosplus_dm23sq   = -1;
    errminosplus_unbanff_00 = -1;
    errminosplus_unbanff_01 = -1;
    errminosplus_unbanff_02 = -1;
    errminosplus_unbanff_03 = -1;
    errminosplus_unbanff_04 = -1;
    errminosplus_unbanff_05 = -1;
    errminosplus_unbanff_06 = -1;
    errminosplus_unbanff_07 = -1;
    errminosplus_unbanff_08 = -1;
    errminosplus_unbanff_09 = -1;
    errminosplus_unbanff_10 = -1;
    errminosplus_unbanff_11 = -1;
    errminosplus_unbanff_12 = -1;
    errminosplus_unbanff_13 = -1;
    errminosplus_unbanff_14 = -1;
    errminosplus_unbanff_15 = -1;
    errminosplus_unbanff_16 = -1;
    errminosplus_unbanff_17 = -1;
    errminosplus_unbanff_18 = -1;
    errminosplus_unbanff_19 = -1;
    errminosplus_unbanff_20 = -1;
    errminosplus_unbanff_21 = -1;
    errminosplus_unbanff_22 = -1;
    errminosplus_banff_00 = -1;
    errminosplus_banff_01 = -1;
    errminosplus_banff_02 = -1;
    errminosplus_banff_03 = -1;
    errminosplus_banff_04 = -1;
    errminosplus_banff_05 = -1;
    errminosplus_banff_06 = -1;
    errminosplus_banff_07 = -1;
    errminosplus_banff_08 = -1;
    errminosplus_banff_09 = -1;
    errminosplus_banff_10 = -1;
    errminosplus_banff_11 = -1;
    errminosplus_banff_12 = -1;
    errminosplus_banff_13 = -1;
    errminosplus_banff_14 = -1;
    errminosplus_banff_15 = -1;
    errminosplus_banff_maqe = -1;
    errminosplus_banff_mares = -1;
    errminosplus_banff_ccqeE0 = -1;
    errminosplus_banff_ccqeE1 = -1;
    errminosplus_banff_ccqeE2 = -1;
    errminosplus_banff_cc1piE0 = -1;
    errminosplus_banff_cc1piE1 = -1;
    errminosplus_xsec_pf = -1;
    errminosplus_xsec_sf = -1;
    errminosplus_xsec_wshape = -1;
    errminosplus_xsec_ccothshape = -1;
    errminosplus_xsec_cc1pishape = -1;
    errminosplus_xsec_pilessdelta = -1;
    errminosplus_xsec_normcccoh = -1;
    errminosplus_xsec_normnc1piC = -1;
    errminosplus_xsec_normnccoh = -1;
    errminosplus_xsec_normncoth = -1;
    errminosplus_xsec_normnuetonumu = -1;
    errminosplus_xsec_normantinutonu = -1;
    errminosplus_skdet_energyscale = -1;
    errminosplus_skdet_ccqerc0 = -1;
    errminosplus_skdet_ccqerc1 = -1;
    errminosplus_skdet_ccqerc2 = -1;
    errminosplus_skdet_ccqeoth = -1;
    errminosplus_skdet_ccnqe = -1;
    errminosplus_skdet_nc = -1;
    errminosplus_skdet_nuecc = -1;
    errminosplus_fsi_000 = -1;
    errminosplus_fsi_001 = -1;
    errminosplus_fsi_002 = -1;
    errminosplus_fsi_003 = -1;
    errminosplus_fsi_004 = -1;
    errminosplus_fsi_005 = -1;
    // Init minus errors
    errminosminus_mixing12 = -1;
    errminosminus_mixing13 = -1;
    errminosminus_mixing23 = -1;
    errminosminus_deltacp  = -1;
    errminosminus_dm12sq   = -1;
    errminosminus_dm23sq   = -1;
    errminosminus_unbanff_00 = -1;
    errminosminus_unbanff_01 = -1;
    errminosminus_unbanff_02 = -1;
    errminosminus_unbanff_03 = -1;
    errminosminus_unbanff_04 = -1;
    errminosminus_unbanff_05 = -1;
    errminosminus_unbanff_06 = -1;
    errminosminus_unbanff_07 = -1;
    errminosminus_unbanff_08 = -1;
    errminosminus_unbanff_09 = -1;
    errminosminus_unbanff_10 = -1;
    errminosminus_unbanff_11 = -1;
    errminosminus_unbanff_12 = -1;
    errminosminus_unbanff_13 = -1;
    errminosminus_unbanff_14 = -1;
    errminosminus_unbanff_15 = -1;
    errminosminus_unbanff_16 = -1;
    errminosminus_unbanff_17 = -1;
    errminosminus_unbanff_18 = -1;
    errminosminus_unbanff_19 = -1;
    errminosminus_unbanff_20 = -1;
    errminosminus_unbanff_21 = -1;
    errminosminus_unbanff_22 = -1;
    errminosminus_banff_00 = -1;
    errminosminus_banff_01 = -1;
    errminosminus_banff_02 = -1;
    errminosminus_banff_03 = -1;
    errminosminus_banff_04 = -1;
    errminosminus_banff_05 = -1;
    errminosminus_banff_06 = -1;
    errminosminus_banff_07 = -1;
    errminosminus_banff_08 = -1;
    errminosminus_banff_09 = -1;
    errminosminus_banff_10 = -1;
    errminosminus_banff_11 = -1;
    errminosminus_banff_12 = -1;
    errminosminus_banff_13 = -1;
    errminosminus_banff_14 = -1;
    errminosminus_banff_15 = -1;
    errminosminus_banff_maqe = -1;
    errminosminus_banff_mares = -1;
    errminosminus_banff_ccqeE0 = -1;
    errminosminus_banff_ccqeE1 = -1;
    errminosminus_banff_ccqeE2 = -1;
    errminosminus_banff_cc1piE0 = -1;
    errminosminus_banff_cc1piE1 = -1;
    errminosminus_xsec_pf = -1;
    errminosminus_xsec_sf = -1;
    errminosminus_xsec_wshape = -1;
    errminosminus_xsec_ccothshape = -1;
    errminosminus_xsec_cc1pishape = -1;
    errminosminus_xsec_pilessdelta = -1;
    errminosminus_xsec_normcccoh = -1;
    errminosminus_xsec_normnc1piC = -1;
    errminosminus_xsec_normnccoh = -1;
    errminosminus_xsec_normncoth = -1;
    errminosminus_xsec_normnuetonumu = -1;
    errminosminus_xsec_normantinutonu = -1;
    errminosminus_skdet_energyscale = -1;
    errminosminus_skdet_ccqerc0 = -1;
    errminosminus_skdet_ccqerc1 = -1;
    errminosminus_skdet_ccqerc2 = -1;
    errminosminus_skdet_ccqeoth = -1;
    errminosminus_skdet_ccnqe = -1;
    errminosminus_skdet_nc = -1;
    errminosminus_skdet_nuecc = -1;
    errminosminus_fsi_000 = -1;
    errminosminus_fsi_001 = -1;
    errminosminus_fsi_002 = -1;
    errminosminus_fsi_003 = -1;
    errminosminus_fsi_004 = -1;
    errminosminus_fsi_005 = -1;
    // Init migrad errors
    errmigrad_mixing12 = -1;
    errmigrad_mixing13 = -1;
    errmigrad_mixing23 = -1;
    errmigrad_deltacp  = -1;
    errmigrad_dm12sq   = -1;
    errmigrad_dm23sq   = -1;
    errmigrad_unbanff_00 = -1;
    errmigrad_unbanff_01 = -1;
    errmigrad_unbanff_02 = -1;
    errmigrad_unbanff_03 = -1;
    errmigrad_unbanff_04 = -1;
    errmigrad_unbanff_05 = -1;
    errmigrad_unbanff_06 = -1;
    errmigrad_unbanff_07 = -1;
    errmigrad_unbanff_08 = -1;
    errmigrad_unbanff_09 = -1;
    errmigrad_unbanff_10 = -1;
    errmigrad_unbanff_11 = -1;
    errmigrad_unbanff_12 = -1;
    errmigrad_unbanff_13 = -1;
    errmigrad_unbanff_14 = -1;
    errmigrad_unbanff_15 = -1;
    errmigrad_unbanff_16 = -1;
    errmigrad_unbanff_17 = -1;
    errmigrad_unbanff_18 = -1;
    errmigrad_unbanff_19 = -1;
    errmigrad_unbanff_20 = -1;
    errmigrad_unbanff_21 = -1;
    errmigrad_unbanff_22 = -1;
    errmigrad_banff_00 = -1;
    errmigrad_banff_01 = -1;
    errmigrad_banff_02 = -1;
    errmigrad_banff_03 = -1;
    errmigrad_banff_04 = -1;
    errmigrad_banff_05 = -1;
    errmigrad_banff_06 = -1;
    errmigrad_banff_07 = -1;
    errmigrad_banff_08 = -1;
    errmigrad_banff_09 = -1;
    errmigrad_banff_10 = -1;
    errmigrad_banff_11 = -1;
    errmigrad_banff_12 = -1;
    errmigrad_banff_13 = -1;
    errmigrad_banff_14 = -1;
    errmigrad_banff_15 = -1;
    errmigrad_banff_maqe = -1;
    errmigrad_banff_mares = -1;
    errmigrad_banff_ccqeE0 = -1;
    errmigrad_banff_ccqeE1 = -1;
    errmigrad_banff_ccqeE2 = -1;
    errmigrad_banff_cc1piE0 = -1;
    errmigrad_banff_cc1piE1 = -1;
    errmigrad_xsec_pf = -1;
    errmigrad_xsec_sf = -1;
    errmigrad_xsec_wshape = -1;
    errmigrad_xsec_ccothshape = -1;
    errmigrad_xsec_cc1pishape = -1;
    errmigrad_xsec_pilessdelta = -1;
    errmigrad_xsec_normcccoh = -1;
    errmigrad_xsec_normnc1piC = -1;
    errmigrad_xsec_normnccoh = -1;
    errmigrad_xsec_normncoth = -1;
    errmigrad_xsec_normnuetonumu = -1;
    errmigrad_xsec_normantinutonu = -1;
    errmigrad_skdet_energyscale = -1;
    errmigrad_skdet_ccqerc0 = -1;
    errmigrad_skdet_ccqerc1 = -1;
    errmigrad_skdet_ccqerc2 = -1;
    errmigrad_skdet_ccqeoth = -1;
    errmigrad_skdet_ccnqe = -1;
    errmigrad_skdet_nc = -1;
    errmigrad_skdet_nuecc = -1;
    errmigrad_fsi_000 = -1;
    errmigrad_fsi_001 = -1;
    errmigrad_fsi_002 = -1;
    errmigrad_fsi_003 = -1;
    errmigrad_fsi_004 = -1;
    errmigrad_fsi_005 = -1;
    // Init hesse errors
    errhesse_mixing12 = -1;
    errhesse_mixing13 = -1;
    errhesse_mixing23 = -1;
    errhesse_deltacp  = -1;
    errhesse_dm12sq   = -1;
    errhesse_dm23sq   = -1;
    errhesse_unbanff_00 = -1;
    errhesse_unbanff_01 = -1;
    errhesse_unbanff_02 = -1;
    errhesse_unbanff_03 = -1;
    errhesse_unbanff_04 = -1;
    errhesse_unbanff_05 = -1;
    errhesse_unbanff_06 = -1;
    errhesse_unbanff_07 = -1;
    errhesse_unbanff_08 = -1;
    errhesse_unbanff_09 = -1;
    errhesse_unbanff_10 = -1;
    errhesse_unbanff_11 = -1;
    errhesse_unbanff_12 = -1;
    errhesse_unbanff_13 = -1;
    errhesse_unbanff_14 = -1;
    errhesse_unbanff_15 = -1;
    errhesse_unbanff_16 = -1;
    errhesse_unbanff_17 = -1;
    errhesse_unbanff_18 = -1;
    errhesse_unbanff_19 = -1;
    errhesse_unbanff_20 = -1;
    errhesse_unbanff_21 = -1;
    errhesse_unbanff_22 = -1;
    errhesse_banff_00 = -1;
    errhesse_banff_01 = -1;
    errhesse_banff_02 = -1;
    errhesse_banff_03 = -1;
    errhesse_banff_04 = -1;
    errhesse_banff_05 = -1;
    errhesse_banff_06 = -1;
    errhesse_banff_07 = -1;
    errhesse_banff_08 = -1;
    errhesse_banff_09 = -1;
    errhesse_banff_10 = -1;
    errhesse_banff_11 = -1;
    errhesse_banff_12 = -1;
    errhesse_banff_13 = -1;
    errhesse_banff_14 = -1;
    errhesse_banff_15 = -1;
    errhesse_banff_maqe = -1;
    errhesse_banff_mares = -1;
    errhesse_banff_ccqeE0 = -1;
    errhesse_banff_ccqeE1 = -1;
    errhesse_banff_ccqeE2 = -1;
    errhesse_banff_cc1piE0 = -1;
    errhesse_banff_cc1piE1 = -1;
    errhesse_xsec_pf = -1;
    errhesse_xsec_sf = -1;
    errhesse_xsec_wshape = -1;
    errhesse_xsec_ccothshape = -1;
    errhesse_xsec_cc1pishape = -1;
    errhesse_xsec_pilessdelta = -1;
    errhesse_xsec_normcccoh = -1;
    errhesse_xsec_normnc1piC = -1;
    errhesse_xsec_normnccoh = -1;
    errhesse_xsec_normncoth = -1;
    errhesse_xsec_normnuetonumu = -1;
    errhesse_xsec_normantinutonu = -1;
    errhesse_skdet_energyscale = -1;
    errhesse_skdet_ccqerc0 = -1;
    errhesse_skdet_ccqerc1 = -1;
    errhesse_skdet_ccqerc2 = -1;
    errhesse_skdet_ccqeoth = -1;
    errhesse_skdet_ccnqe = -1;
    errhesse_skdet_nc = -1;
    errhesse_skdet_nuecc = -1;
    errhesse_fsi_000 = -1;
    errhesse_fsi_001 = -1;
    errhesse_fsi_002 = -1;
    errhesse_fsi_003 = -1;
    errhesse_fsi_004 = -1;
    errhesse_fsi_005 = -1;
    //get the bestfit values
    if(hyp == kHypOsc3F) {
      bestfit_mixing12 = fparm->BestFitValue(kPhypMixing12);
      bestfit_mixing13 = fparm->BestFitValue(kPhypMixing13);
      bestfit_mixing23 = fparm->BestFitValue(kPhypMixing23);
      bestfit_deltacp  = fparm->BestFitValue(kPhypDeltaCP);
      bestfit_dm12sq   = fparm->BestFitValue(kPhypDm12sq);
      bestfit_dm23sq   = fparm->BestFitValue(kPhypDm23sq);
    } 
    else if(hyp == kHypNuMuDisapp2F) {
      bestfit_mixing23 = fparm->BestFitValue(kPhypMixing23);
      bestfit_dm23sq   = fparm->BestFitValue(kPhypDm23sq);
    }
    if(fparm->Exists(conv::kNuisUnBANFF00)) bestfit_unbanff_00 = fparm->BestFitValue(conv::kNuisUnBANFF00);
    if(fparm->Exists(conv::kNuisUnBANFF01)) bestfit_unbanff_01 = fparm->BestFitValue(conv::kNuisUnBANFF01);
    if(fparm->Exists(conv::kNuisUnBANFF02)) bestfit_unbanff_02 = fparm->BestFitValue(conv::kNuisUnBANFF02);
    if(fparm->Exists(conv::kNuisUnBANFF03)) bestfit_unbanff_03 = fparm->BestFitValue(conv::kNuisUnBANFF03);
    if(fparm->Exists(conv::kNuisUnBANFF04)) bestfit_unbanff_04 = fparm->BestFitValue(conv::kNuisUnBANFF04);
    if(fparm->Exists(conv::kNuisUnBANFF05)) bestfit_unbanff_05 = fparm->BestFitValue(conv::kNuisUnBANFF05);
    if(fparm->Exists(conv::kNuisUnBANFF06)) bestfit_unbanff_06 = fparm->BestFitValue(conv::kNuisUnBANFF06);
    if(fparm->Exists(conv::kNuisUnBANFF07)) bestfit_unbanff_07 = fparm->BestFitValue(conv::kNuisUnBANFF07);
    if(fparm->Exists(conv::kNuisUnBANFF08)) bestfit_unbanff_08 = fparm->BestFitValue(conv::kNuisUnBANFF08);
    if(fparm->Exists(conv::kNuisUnBANFF09)) bestfit_unbanff_09 = fparm->BestFitValue(conv::kNuisUnBANFF09);
    if(fparm->Exists(conv::kNuisUnBANFF10)) bestfit_unbanff_10 = fparm->BestFitValue(conv::kNuisUnBANFF10);
    if(fparm->Exists(conv::kNuisUnBANFF11)) bestfit_unbanff_11 = fparm->BestFitValue(conv::kNuisUnBANFF11);
    if(fparm->Exists(conv::kNuisUnBANFF12)) bestfit_unbanff_12 = fparm->BestFitValue(conv::kNuisUnBANFF12);
    if(fparm->Exists(conv::kNuisUnBANFF13)) bestfit_unbanff_13 = fparm->BestFitValue(conv::kNuisUnBANFF13);
    if(fparm->Exists(conv::kNuisUnBANFF14)) bestfit_unbanff_14 = fparm->BestFitValue(conv::kNuisUnBANFF14);
    if(fparm->Exists(conv::kNuisUnBANFF15)) bestfit_unbanff_15 = fparm->BestFitValue(conv::kNuisUnBANFF15);
    if(fparm->Exists(conv::kNuisUnBANFF16)) bestfit_unbanff_16 = fparm->BestFitValue(conv::kNuisUnBANFF16);
    if(fparm->Exists(conv::kNuisUnBANFF17)) bestfit_unbanff_17 = fparm->BestFitValue(conv::kNuisUnBANFF17);
    if(fparm->Exists(conv::kNuisUnBANFF18)) bestfit_unbanff_18 = fparm->BestFitValue(conv::kNuisUnBANFF18);
    if(fparm->Exists(conv::kNuisUnBANFF19)) bestfit_unbanff_19 = fparm->BestFitValue(conv::kNuisUnBANFF19);
    if(fparm->Exists(conv::kNuisUnBANFF20)) bestfit_unbanff_20 = fparm->BestFitValue(conv::kNuisUnBANFF20);
    if(fparm->Exists(conv::kNuisUnBANFF21)) bestfit_unbanff_21 = fparm->BestFitValue(conv::kNuisUnBANFF21);
    if(fparm->Exists(conv::kNuisUnBANFF22)) bestfit_unbanff_22 = fparm->BestFitValue(conv::kNuisUnBANFF22);
    if(fparm->Exists(conv::kNuisBANFF00)) bestfit_banff_00 = fparm->BestFitValue(conv::kNuisBANFF00);
    if(fparm->Exists(conv::kNuisBANFF01)) bestfit_banff_01 = fparm->BestFitValue(conv::kNuisBANFF01);
    if(fparm->Exists(conv::kNuisBANFF02)) bestfit_banff_02 = fparm->BestFitValue(conv::kNuisBANFF02);
    if(fparm->Exists(conv::kNuisBANFF03)) bestfit_banff_03 = fparm->BestFitValue(conv::kNuisBANFF03);
    if(fparm->Exists(conv::kNuisBANFF04)) bestfit_banff_04 = fparm->BestFitValue(conv::kNuisBANFF04);
    if(fparm->Exists(conv::kNuisBANFF05)) bestfit_banff_05 = fparm->BestFitValue(conv::kNuisBANFF05);
    if(fparm->Exists(conv::kNuisBANFF06)) bestfit_banff_06 = fparm->BestFitValue(conv::kNuisBANFF06);
    if(fparm->Exists(conv::kNuisBANFF07)) bestfit_banff_07 = fparm->BestFitValue(conv::kNuisBANFF07);
    if(fparm->Exists(conv::kNuisBANFF08)) bestfit_banff_08 = fparm->BestFitValue(conv::kNuisBANFF08);
    if(fparm->Exists(conv::kNuisBANFF09)) bestfit_banff_09 = fparm->BestFitValue(conv::kNuisBANFF09);
    if(fparm->Exists(conv::kNuisBANFF10)) bestfit_banff_10 = fparm->BestFitValue(conv::kNuisBANFF10);
    if(fparm->Exists(conv::kNuisBANFF11)) bestfit_banff_11 = fparm->BestFitValue(conv::kNuisBANFF11);
    if(fparm->Exists(conv::kNuisBANFF12)) bestfit_banff_12 = fparm->BestFitValue(conv::kNuisBANFF12);
    if(fparm->Exists(conv::kNuisBANFF13)) bestfit_banff_13 = fparm->BestFitValue(conv::kNuisBANFF13);
    if(fparm->Exists(conv::kNuisBANFF14)) bestfit_banff_14 = fparm->BestFitValue(conv::kNuisBANFF14);
    if(fparm->Exists(conv::kNuisBANFF15)) bestfit_banff_15 = fparm->BestFitValue(conv::kNuisBANFF15);
    if(fparm->Exists(conv::kNuisBANFFMAQE)) bestfit_banff_maqe = fparm->BestFitValue(conv::kNuisBANFFMAQE);
    if(fparm->Exists(conv::kNuisBANFFMARES)) bestfit_banff_mares = fparm->BestFitValue(conv::kNuisBANFFMARES);
    if(fparm->Exists(conv::kNuisBANFFCCQEE0)) bestfit_banff_ccqeE0 = fparm->BestFitValue(conv::kNuisBANFFCCQEE0);
    if(fparm->Exists(conv::kNuisBANFFCCQEE1)) bestfit_banff_ccqeE1 = fparm->BestFitValue(conv::kNuisBANFFCCQEE1);
    if(fparm->Exists(conv::kNuisBANFFCCQEE2)) bestfit_banff_ccqeE2 = fparm->BestFitValue(conv::kNuisBANFFCCQEE2);
    if(fparm->Exists(conv::kNuisBANFFCC1piE0)) bestfit_banff_cc1piE0 = fparm->BestFitValue(conv::kNuisBANFFCC1piE0);
    if(fparm->Exists(conv::kNuisBANFFCC1piE1)) bestfit_banff_cc1piE1 = fparm->BestFitValue(conv::kNuisBANFFCC1piE1);
    if(fparm->Exists(conv::kNuisXSecFermiMom)) bestfit_xsec_pf = fparm->BestFitValue(conv::kNuisXSecFermiMom);
    if(fparm->Exists(conv::kNuisXSecSpecFunc)) bestfit_xsec_sf = fparm->BestFitValue(conv::kNuisXSecSpecFunc);
    if(fparm->Exists(conv::kNuisXSecBindE)) bestfit_xsec_binde = fparm->BestFitValue(conv::kNuisXSecBindE);
    if(fparm->Exists(conv::kNuisXSecWshape)) bestfit_xsec_wshape = fparm->BestFitValue(conv::kNuisXSecWshape);
    if(fparm->Exists(conv::kNuisXSecShapeCCoth)) bestfit_xsec_ccothshape = fparm->BestFitValue(conv::kNuisXSecShapeCCoth);
    if(fparm->Exists(conv::kNuisXSecPiLessDelta)) bestfit_xsec_pilessdelta = fparm->BestFitValue(conv::kNuisXSecPiLessDelta);
    if(fparm->Exists(conv::kNuisXSecENuShape)) bestfit_xsec_enushape = fparm->BestFitValue(conv::kNuisXSecENuShape);
    if(fparm->Exists(conv::kNuisXSecNormCCcoh)) bestfit_xsec_normcccoh = fparm->BestFitValue(conv::kNuisXSecNormCCcoh);
    if(fparm->Exists(conv::kNuisXSecNormNC1piC)) bestfit_xsec_normnc1piC = fparm->BestFitValue(conv::kNuisXSecNormNC1piC);
    if(fparm->Exists(conv::kNuisXSecNormNCcoh)) bestfit_xsec_normnccoh = fparm->BestFitValue(conv::kNuisXSecNormNCcoh);
    if(fparm->Exists(conv::kNuisXSecNormNCoth)) bestfit_xsec_normncoth = fparm->BestFitValue(conv::kNuisXSecNormNCoth);
    if(fparm->Exists(conv::kNuisXSecNormNuEToNuMu)) bestfit_xsec_normnuetonumu = fparm->BestFitValue(conv::kNuisXSecNormNuEToNuMu);
    if(fparm->Exists(conv::kNuisXSecNormAntiNuToNu)) bestfit_xsec_normantinutonu = fparm->BestFitValue(conv::kNuisXSecNormAntiNuToNu);
    if(fparm->Exists(conv::kNuisSuperkEnergyScale)) bestfit_skdet_energyscale = fparm->BestFitValue(conv::kNuisSuperkEnergyScale);
    if(fparm->Exists(conv::kNuisSuperkEff1RmuCCQErc0u)) bestfit_skdet_ccqerc0 = fparm->BestFitValue(conv::kNuisSuperkEff1RmuCCQErc0u);
    if(fparm->Exists(conv::kNuisSuperkEff1RmuCCQErc1u)) bestfit_skdet_ccqerc1 = fparm->BestFitValue(conv::kNuisSuperkEff1RmuCCQErc1u);
    if(fparm->Exists(conv::kNuisSuperkEff1RmuCCQErc2u)) bestfit_skdet_ccqerc2 = fparm->BestFitValue(conv::kNuisSuperkEff1RmuCCQErc2u);
    if(fparm->Exists(conv::kNuisSuperkEff1RmuCCQEoth)) bestfit_skdet_ccqeoth = fparm->BestFitValue(conv::kNuisSuperkEff1RmuCCQEoth);
    if(fparm->Exists(conv::kNuisSuperkEff1RmuCCnonQE)) bestfit_skdet_ccnqe = fparm->BestFitValue(conv::kNuisSuperkEff1RmuCCnonQE);
    if(fparm->Exists(conv::kNuisSuperkEff1RmuNC)) bestfit_skdet_nc = fparm->BestFitValue(conv::kNuisSuperkEff1RmuNC);
    if(fparm->Exists(conv::kNuisSuperkEff1RmuNueCC)) bestfit_skdet_nuecc = fparm->BestFitValue(conv::kNuisSuperkEff1RmuNueCC);
    if(fparm->Exists(conv::kNuisFSI000)) bestfit_fsi_000 = fparm->BestFitValue(conv::kNuisFSI000);
    if(fparm->Exists(conv::kNuisFSI001)) bestfit_fsi_001 = fparm->BestFitValue(conv::kNuisFSI001);
    if(fparm->Exists(conv::kNuisFSI002)) bestfit_fsi_002 = fparm->BestFitValue(conv::kNuisFSI002);
    if(fparm->Exists(conv::kNuisFSI003)) bestfit_fsi_003 = fparm->BestFitValue(conv::kNuisFSI003);
    if(fparm->Exists(conv::kNuisFSI004)) bestfit_fsi_004 = fparm->BestFitValue(conv::kNuisFSI004);
    if(fparm->Exists(conv::kNuisFSI005)) bestfit_fsi_005 = fparm->BestFitValue(conv::kNuisFSI005);
    //MINOS errors
    if(fparm->MinosErrorsSet()) {
      //get the plus errors
      if(hyp == kHypOsc3F) {
	errminosplus_mixing12 = fparm->PlusError(kPhypMixing12);
	errminosplus_mixing13 = fparm->PlusError(kPhypMixing13);
	errminosplus_mixing23 = fparm->PlusError(kPhypMixing23);
	errminosplus_deltacp  = fparm->PlusError(kPhypDeltaCP);
	errminosplus_dm12sq   = fparm->PlusError(kPhypDm12sq);
	errminosplus_dm23sq   = fparm->PlusError(kPhypDm23sq);
      } else
	if(hyp == kHypNuMuDisapp2F) {
	  errminosplus_mixing23 = fparm->PlusError(kPhypMixing23);
	  errminosplus_dm23sq   = fparm->PlusError(kPhypDm23sq);
	}
      if(fparm->Exists(conv::kNuisUnBANFF00)) errminosplus_unbanff_00 = fparm->PlusError(conv::kNuisUnBANFF00);
      if(fparm->Exists(conv::kNuisUnBANFF01)) errminosplus_unbanff_01 = fparm->PlusError(conv::kNuisUnBANFF01);
      if(fparm->Exists(conv::kNuisUnBANFF02)) errminosplus_unbanff_02 = fparm->PlusError(conv::kNuisUnBANFF02);
      if(fparm->Exists(conv::kNuisUnBANFF03)) errminosplus_unbanff_03 = fparm->PlusError(conv::kNuisUnBANFF03);
      if(fparm->Exists(conv::kNuisUnBANFF04)) errminosplus_unbanff_04 = fparm->PlusError(conv::kNuisUnBANFF04);
      if(fparm->Exists(conv::kNuisUnBANFF05)) errminosplus_unbanff_05 = fparm->PlusError(conv::kNuisUnBANFF05);
      if(fparm->Exists(conv::kNuisUnBANFF06)) errminosplus_unbanff_06 = fparm->PlusError(conv::kNuisUnBANFF06);
      if(fparm->Exists(conv::kNuisUnBANFF07)) errminosplus_unbanff_07 = fparm->PlusError(conv::kNuisUnBANFF07);
      if(fparm->Exists(conv::kNuisUnBANFF08)) errminosplus_unbanff_08 = fparm->PlusError(conv::kNuisUnBANFF08);
      if(fparm->Exists(conv::kNuisUnBANFF09)) errminosplus_unbanff_09 = fparm->PlusError(conv::kNuisUnBANFF09);
      if(fparm->Exists(conv::kNuisUnBANFF10)) errminosplus_unbanff_10 = fparm->PlusError(conv::kNuisUnBANFF10);
      if(fparm->Exists(conv::kNuisUnBANFF11)) errminosplus_unbanff_11 = fparm->PlusError(conv::kNuisUnBANFF11);
      if(fparm->Exists(conv::kNuisUnBANFF12)) errminosplus_unbanff_12 = fparm->PlusError(conv::kNuisUnBANFF12);
      if(fparm->Exists(conv::kNuisUnBANFF13)) errminosplus_unbanff_13 = fparm->PlusError(conv::kNuisUnBANFF13);
      if(fparm->Exists(conv::kNuisUnBANFF14)) errminosplus_unbanff_14 = fparm->PlusError(conv::kNuisUnBANFF14);
      if(fparm->Exists(conv::kNuisUnBANFF15)) errminosplus_unbanff_15 = fparm->PlusError(conv::kNuisUnBANFF15);
      if(fparm->Exists(conv::kNuisUnBANFF16)) errminosplus_unbanff_16 = fparm->PlusError(conv::kNuisUnBANFF16);
      if(fparm->Exists(conv::kNuisUnBANFF17)) errminosplus_unbanff_17 = fparm->PlusError(conv::kNuisUnBANFF17);
      if(fparm->Exists(conv::kNuisUnBANFF18)) errminosplus_unbanff_18 = fparm->PlusError(conv::kNuisUnBANFF18);
      if(fparm->Exists(conv::kNuisUnBANFF19)) errminosplus_unbanff_19 = fparm->PlusError(conv::kNuisUnBANFF19);
      if(fparm->Exists(conv::kNuisUnBANFF20)) errminosplus_unbanff_20 = fparm->PlusError(conv::kNuisUnBANFF20);
      if(fparm->Exists(conv::kNuisUnBANFF21)) errminosplus_unbanff_21 = fparm->PlusError(conv::kNuisUnBANFF21);
      if(fparm->Exists(conv::kNuisUnBANFF22)) errminosplus_unbanff_22 = fparm->PlusError(conv::kNuisUnBANFF22);
      if(fparm->Exists(conv::kNuisBANFF00)) errminosplus_banff_00 = fparm->PlusError(conv::kNuisBANFF00);
      if(fparm->Exists(conv::kNuisBANFF01)) errminosplus_banff_01 = fparm->PlusError(conv::kNuisBANFF01);
      if(fparm->Exists(conv::kNuisBANFF02)) errminosplus_banff_02 = fparm->PlusError(conv::kNuisBANFF02);
      if(fparm->Exists(conv::kNuisBANFF03)) errminosplus_banff_03 = fparm->PlusError(conv::kNuisBANFF03);
      if(fparm->Exists(conv::kNuisBANFF04)) errminosplus_banff_04 = fparm->PlusError(conv::kNuisBANFF04);
      if(fparm->Exists(conv::kNuisBANFF05)) errminosplus_banff_05 = fparm->PlusError(conv::kNuisBANFF05);
      if(fparm->Exists(conv::kNuisBANFF06)) errminosplus_banff_06 = fparm->PlusError(conv::kNuisBANFF06);
      if(fparm->Exists(conv::kNuisBANFF07)) errminosplus_banff_07 = fparm->PlusError(conv::kNuisBANFF07);
      if(fparm->Exists(conv::kNuisBANFF08)) errminosplus_banff_08 = fparm->PlusError(conv::kNuisBANFF08);
      if(fparm->Exists(conv::kNuisBANFF09)) errminosplus_banff_09 = fparm->PlusError(conv::kNuisBANFF09);
      if(fparm->Exists(conv::kNuisBANFF10)) errminosplus_banff_10 = fparm->PlusError(conv::kNuisBANFF10);
      if(fparm->Exists(conv::kNuisBANFF11)) errminosplus_banff_11 = fparm->PlusError(conv::kNuisBANFF11);
      if(fparm->Exists(conv::kNuisBANFF12)) errminosplus_banff_12 = fparm->PlusError(conv::kNuisBANFF12);
      if(fparm->Exists(conv::kNuisBANFF13)) errminosplus_banff_13 = fparm->PlusError(conv::kNuisBANFF13);
      if(fparm->Exists(conv::kNuisBANFF14)) errminosplus_banff_14 = fparm->PlusError(conv::kNuisBANFF14);
      if(fparm->Exists(conv::kNuisBANFF15)) errminosplus_banff_15 = fparm->PlusError(conv::kNuisBANFF15);
      if(fparm->Exists(conv::kNuisBANFFMAQE)) errminosplus_banff_maqe = fparm->PlusError(conv::kNuisBANFFMAQE);
      if(fparm->Exists(conv::kNuisBANFFMARES)) errminosplus_banff_mares = fparm->PlusError(conv::kNuisBANFFMARES);
      if(fparm->Exists(conv::kNuisBANFFCCQEE0)) errminosplus_banff_ccqeE0 = fparm->PlusError(conv::kNuisBANFFCCQEE0);
      if(fparm->Exists(conv::kNuisBANFFCCQEE1)) errminosplus_banff_ccqeE1 = fparm->PlusError(conv::kNuisBANFFCCQEE1);
      if(fparm->Exists(conv::kNuisBANFFCCQEE2)) errminosplus_banff_ccqeE2 = fparm->PlusError(conv::kNuisBANFFCCQEE2);
      if(fparm->Exists(conv::kNuisBANFFCC1piE0)) errminosplus_banff_cc1piE0 = fparm->PlusError(conv::kNuisBANFFCC1piE0);
      if(fparm->Exists(conv::kNuisBANFFCC1piE1)) errminosplus_banff_cc1piE1 = fparm->PlusError(conv::kNuisBANFFCC1piE1);
      if(fparm->Exists(conv::kNuisXSecFermiMom)) errminosplus_xsec_pf = fparm->PlusError(conv::kNuisXSecFermiMom);
      if(fparm->Exists(conv::kNuisXSecSpecFunc)) errminosplus_xsec_sf = fparm->PlusError(conv::kNuisXSecSpecFunc);
      if(fparm->Exists(conv::kNuisXSecBindE)) errminosplus_xsec_binde = fparm->PlusError(conv::kNuisXSecBindE);
      if(fparm->Exists(conv::kNuisXSecWshape)) errminosplus_xsec_wshape = fparm->PlusError(conv::kNuisXSecWshape);
      if(fparm->Exists(conv::kNuisXSecShapeCCoth)) errminosplus_xsec_ccothshape = fparm->PlusError(conv::kNuisXSecShapeCCoth);
      if(fparm->Exists(conv::kNuisXSecPiLessDelta)) errminosplus_xsec_pilessdelta = fparm->PlusError(conv::kNuisXSecPiLessDelta);
      if(fparm->Exists(conv::kNuisXSecENuShape)) errminosplus_xsec_enushape = fparm->PlusError(conv::kNuisXSecENuShape);
      if(fparm->Exists(conv::kNuisXSecNormCCcoh)) errminosplus_xsec_normcccoh = fparm->PlusError(conv::kNuisXSecNormCCcoh);
      if(fparm->Exists(conv::kNuisXSecNormNC1piC)) errminosplus_xsec_normnc1piC = fparm->PlusError(conv::kNuisXSecNormNC1piC);
      if(fparm->Exists(conv::kNuisXSecNormNCcoh)) errminosplus_xsec_normnccoh = fparm->PlusError(conv::kNuisXSecNormNCcoh);
      if(fparm->Exists(conv::kNuisXSecNormNCoth)) errminosplus_xsec_normncoth = fparm->PlusError(conv::kNuisXSecNormNCoth);
      if(fparm->Exists(conv::kNuisXSecNormNuEToNuMu)) errminosplus_xsec_normnuetonumu = fparm->PlusError(conv::kNuisXSecNormNuEToNuMu);
      if(fparm->Exists(conv::kNuisXSecNormAntiNuToNu)) errminosplus_xsec_normantinutonu = fparm->PlusError(conv::kNuisXSecNormAntiNuToNu);
      if(fparm->Exists(conv::kNuisSuperkEnergyScale)) errminosplus_skdet_energyscale = fparm->PlusError(conv::kNuisSuperkEnergyScale);
      if(fparm->Exists(conv::kNuisSuperkEff1RmuCCQErc0u)) errminosplus_skdet_ccqerc0 = fparm->PlusError(conv::kNuisSuperkEff1RmuCCQErc0u);
      if(fparm->Exists(conv::kNuisSuperkEff1RmuCCQErc1u)) errminosplus_skdet_ccqerc1 = fparm->PlusError(conv::kNuisSuperkEff1RmuCCQErc1u);
      if(fparm->Exists(conv::kNuisSuperkEff1RmuCCQErc2u)) errminosplus_skdet_ccqerc2 = fparm->PlusError(conv::kNuisSuperkEff1RmuCCQErc2u);
      if(fparm->Exists(conv::kNuisSuperkEff1RmuCCQEoth)) errminosplus_skdet_ccqeoth = fparm->PlusError(conv::kNuisSuperkEff1RmuCCQEoth);
      if(fparm->Exists(conv::kNuisSuperkEff1RmuCCnonQE)) errminosplus_skdet_ccnqe = fparm->PlusError(conv::kNuisSuperkEff1RmuCCnonQE);
      if(fparm->Exists(conv::kNuisSuperkEff1RmuNC)) errminosplus_skdet_nc = fparm->PlusError(conv::kNuisSuperkEff1RmuNC);
      if(fparm->Exists(conv::kNuisSuperkEff1RmuNueCC)) errminosplus_skdet_nuecc = fparm->PlusError(conv::kNuisSuperkEff1RmuNueCC);
      if(fparm->Exists(conv::kNuisFSI000)) errminosplus_fsi_000 = fparm->PlusError(conv::kNuisFSI000);
      if(fparm->Exists(conv::kNuisFSI001)) errminosplus_fsi_001 = fparm->PlusError(conv::kNuisFSI001);
      if(fparm->Exists(conv::kNuisFSI002)) errminosplus_fsi_002 = fparm->PlusError(conv::kNuisFSI002);
      if(fparm->Exists(conv::kNuisFSI003)) errminosplus_fsi_003 = fparm->PlusError(conv::kNuisFSI003);
      if(fparm->Exists(conv::kNuisFSI004)) errminosplus_fsi_004 = fparm->PlusError(conv::kNuisFSI004);
      if(fparm->Exists(conv::kNuisFSI005)) errminosplus_fsi_005 = fparm->PlusError(conv::kNuisFSI005);
      //get the minus errors
      if(hyp == kHypOsc3F) {
	errminosminus_mixing12 = fparm->MinusError(kPhypMixing12);
	errminosminus_mixing13 = fparm->MinusError(kPhypMixing13);
	errminosminus_mixing23 = fparm->MinusError(kPhypMixing23);
	errminosminus_deltacp  = fparm->MinusError(kPhypDeltaCP);
	errminosminus_dm12sq   = fparm->MinusError(kPhypDm12sq);
	errminosminus_dm23sq   = fparm->MinusError(kPhypDm23sq);
      }
      else if(hyp == kHypNuMuDisapp2F) {
	errminosminus_mixing23 = fparm->MinusError(kPhypMixing23);
	errminosminus_dm23sq   = fparm->MinusError(kPhypDm23sq);
      }
      if(fparm->Exists(conv::kNuisUnBANFF00)) errminosminus_unbanff_00 = fparm->MinusError(conv::kNuisUnBANFF00);
      if(fparm->Exists(conv::kNuisUnBANFF01)) errminosminus_unbanff_01 = fparm->MinusError(conv::kNuisUnBANFF01);
      if(fparm->Exists(conv::kNuisUnBANFF02)) errminosminus_unbanff_02 = fparm->MinusError(conv::kNuisUnBANFF02);
      if(fparm->Exists(conv::kNuisUnBANFF03)) errminosminus_unbanff_03 = fparm->MinusError(conv::kNuisUnBANFF03);
      if(fparm->Exists(conv::kNuisUnBANFF04)) errminosminus_unbanff_04 = fparm->MinusError(conv::kNuisUnBANFF04);
      if(fparm->Exists(conv::kNuisUnBANFF05)) errminosminus_unbanff_05 = fparm->MinusError(conv::kNuisUnBANFF05);
      if(fparm->Exists(conv::kNuisUnBANFF06)) errminosminus_unbanff_06 = fparm->MinusError(conv::kNuisUnBANFF06);
      if(fparm->Exists(conv::kNuisUnBANFF07)) errminosminus_unbanff_07 = fparm->MinusError(conv::kNuisUnBANFF07);
      if(fparm->Exists(conv::kNuisUnBANFF08)) errminosminus_unbanff_08 = fparm->MinusError(conv::kNuisUnBANFF08);
      if(fparm->Exists(conv::kNuisUnBANFF09)) errminosminus_unbanff_09 = fparm->MinusError(conv::kNuisUnBANFF09);
      if(fparm->Exists(conv::kNuisUnBANFF10)) errminosminus_unbanff_10 = fparm->MinusError(conv::kNuisUnBANFF10);
      if(fparm->Exists(conv::kNuisUnBANFF11)) errminosminus_unbanff_11 = fparm->MinusError(conv::kNuisUnBANFF11);
      if(fparm->Exists(conv::kNuisUnBANFF12)) errminosminus_unbanff_12 = fparm->MinusError(conv::kNuisUnBANFF12);
      if(fparm->Exists(conv::kNuisUnBANFF13)) errminosminus_unbanff_13 = fparm->MinusError(conv::kNuisUnBANFF13);
      if(fparm->Exists(conv::kNuisUnBANFF14)) errminosminus_unbanff_14 = fparm->MinusError(conv::kNuisUnBANFF14);
      if(fparm->Exists(conv::kNuisUnBANFF15)) errminosminus_unbanff_15 = fparm->MinusError(conv::kNuisUnBANFF15);
      if(fparm->Exists(conv::kNuisUnBANFF16)) errminosminus_unbanff_16 = fparm->MinusError(conv::kNuisUnBANFF16);
      if(fparm->Exists(conv::kNuisUnBANFF17)) errminosminus_unbanff_17 = fparm->MinusError(conv::kNuisUnBANFF17);
      if(fparm->Exists(conv::kNuisUnBANFF18)) errminosminus_unbanff_18 = fparm->MinusError(conv::kNuisUnBANFF18);
      if(fparm->Exists(conv::kNuisUnBANFF19)) errminosminus_unbanff_19 = fparm->MinusError(conv::kNuisUnBANFF19);
      if(fparm->Exists(conv::kNuisUnBANFF20)) errminosminus_unbanff_20 = fparm->MinusError(conv::kNuisUnBANFF20);
      if(fparm->Exists(conv::kNuisUnBANFF21)) errminosminus_unbanff_21 = fparm->MinusError(conv::kNuisUnBANFF21);
      if(fparm->Exists(conv::kNuisUnBANFF22)) errminosminus_unbanff_22 = fparm->MinusError(conv::kNuisUnBANFF22);
      if(fparm->Exists(conv::kNuisBANFF00)) errminosminus_banff_00 = fparm->MinusError(conv::kNuisBANFF00);
      if(fparm->Exists(conv::kNuisBANFF01)) errminosminus_banff_01 = fparm->MinusError(conv::kNuisBANFF01);
      if(fparm->Exists(conv::kNuisBANFF02)) errminosminus_banff_02 = fparm->MinusError(conv::kNuisBANFF02);
      if(fparm->Exists(conv::kNuisBANFF03)) errminosminus_banff_03 = fparm->MinusError(conv::kNuisBANFF03);
      if(fparm->Exists(conv::kNuisBANFF04)) errminosminus_banff_04 = fparm->MinusError(conv::kNuisBANFF04);
      if(fparm->Exists(conv::kNuisBANFF05)) errminosminus_banff_05 = fparm->MinusError(conv::kNuisBANFF05);
      if(fparm->Exists(conv::kNuisBANFF06)) errminosminus_banff_06 = fparm->MinusError(conv::kNuisBANFF06);
      if(fparm->Exists(conv::kNuisBANFF07)) errminosminus_banff_07 = fparm->MinusError(conv::kNuisBANFF07);
      if(fparm->Exists(conv::kNuisBANFF08)) errminosminus_banff_08 = fparm->MinusError(conv::kNuisBANFF08);
      if(fparm->Exists(conv::kNuisBANFF09)) errminosminus_banff_09 = fparm->MinusError(conv::kNuisBANFF09);
      if(fparm->Exists(conv::kNuisBANFF10)) errminosminus_banff_10 = fparm->MinusError(conv::kNuisBANFF10);
      if(fparm->Exists(conv::kNuisBANFF11)) errminosminus_banff_11 = fparm->MinusError(conv::kNuisBANFF11);
      if(fparm->Exists(conv::kNuisBANFF12)) errminosminus_banff_12 = fparm->MinusError(conv::kNuisBANFF12);
      if(fparm->Exists(conv::kNuisBANFF13)) errminosminus_banff_13 = fparm->MinusError(conv::kNuisBANFF13);
      if(fparm->Exists(conv::kNuisBANFF14)) errminosminus_banff_14 = fparm->MinusError(conv::kNuisBANFF14);
      if(fparm->Exists(conv::kNuisBANFF15)) errminosminus_banff_15 = fparm->MinusError(conv::kNuisBANFF15);
      if(fparm->Exists(conv::kNuisBANFFMAQE)) errminosminus_banff_maqe = fparm->MinusError(conv::kNuisBANFFMAQE);
      if(fparm->Exists(conv::kNuisBANFFMARES)) errminosminus_banff_mares = fparm->MinusError(conv::kNuisBANFFMARES);
      if(fparm->Exists(conv::kNuisBANFFCCQEE0)) errminosminus_banff_ccqeE0 = fparm->MinusError(conv::kNuisBANFFCCQEE0);
      if(fparm->Exists(conv::kNuisBANFFCCQEE1)) errminosminus_banff_ccqeE1 = fparm->MinusError(conv::kNuisBANFFCCQEE1);
      if(fparm->Exists(conv::kNuisBANFFCCQEE2)) errminosminus_banff_ccqeE2 = fparm->MinusError(conv::kNuisBANFFCCQEE2);
      if(fparm->Exists(conv::kNuisBANFFCC1piE0)) errminosminus_banff_cc1piE0 = fparm->MinusError(conv::kNuisBANFFCC1piE0);
      if(fparm->Exists(conv::kNuisBANFFCC1piE1)) errminosminus_banff_cc1piE1 = fparm->MinusError(conv::kNuisBANFFCC1piE1);
      if(fparm->Exists(conv::kNuisXSecFermiMom)) errminosminus_xsec_pf = fparm->MinusError(conv::kNuisXSecFermiMom);
      if(fparm->Exists(conv::kNuisXSecSpecFunc)) errminosminus_xsec_sf = fparm->MinusError(conv::kNuisXSecSpecFunc);
      if(fparm->Exists(conv::kNuisXSecBindE)) errminosminus_xsec_binde = fparm->MinusError(conv::kNuisXSecBindE);
      if(fparm->Exists(conv::kNuisXSecWshape)) errminosminus_xsec_wshape = fparm->MinusError(conv::kNuisXSecWshape);
      if(fparm->Exists(conv::kNuisXSecShapeCCoth)) errminosminus_xsec_ccothshape = fparm->MinusError(conv::kNuisXSecShapeCCoth);
      if(fparm->Exists(conv::kNuisXSecPiLessDelta)) errminosminus_xsec_pilessdelta = fparm->MinusError(conv::kNuisXSecPiLessDelta);
      if(fparm->Exists(conv::kNuisXSecENuShape)) errminosminus_xsec_enushape = fparm->MinusError(conv::kNuisXSecENuShape);
      if(fparm->Exists(conv::kNuisXSecNormCCcoh)) errminosminus_xsec_normcccoh = fparm->MinusError(conv::kNuisXSecNormCCcoh);
      if(fparm->Exists(conv::kNuisXSecNormNC1piC)) errminosminus_xsec_normnc1piC = fparm->MinusError(conv::kNuisXSecNormNC1piC);
      if(fparm->Exists(conv::kNuisXSecNormNCcoh)) errminosminus_xsec_normnccoh = fparm->MinusError(conv::kNuisXSecNormNCcoh);
      if(fparm->Exists(conv::kNuisXSecNormNCoth)) errminosminus_xsec_normncoth = fparm->MinusError(conv::kNuisXSecNormNCoth);
      if(fparm->Exists(conv::kNuisXSecNormNuEToNuMu)) errminosminus_xsec_normnuetonumu = fparm->MinusError(conv::kNuisXSecNormNuEToNuMu);
      if(fparm->Exists(conv::kNuisXSecNormAntiNuToNu)) errminosminus_xsec_normantinutonu = fparm->MinusError(conv::kNuisXSecNormAntiNuToNu);
      if(fparm->Exists(conv::kNuisSuperkEnergyScale)) errminosminus_skdet_energyscale = fparm->MinusError(conv::kNuisSuperkEnergyScale);
      if(fparm->Exists(conv::kNuisSuperkEff1RmuCCQErc0u)) errminosminus_skdet_ccqerc0 = fparm->MinusError(conv::kNuisSuperkEff1RmuCCQErc0u);
      if(fparm->Exists(conv::kNuisSuperkEff1RmuCCQErc1u)) errminosminus_skdet_ccqerc1 = fparm->MinusError(conv::kNuisSuperkEff1RmuCCQErc1u);
      if(fparm->Exists(conv::kNuisSuperkEff1RmuCCQErc2u)) errminosminus_skdet_ccqerc2 = fparm->MinusError(conv::kNuisSuperkEff1RmuCCQErc2u);
      if(fparm->Exists(conv::kNuisSuperkEff1RmuCCQEoth)) errminosminus_skdet_ccqeoth = fparm->MinusError(conv::kNuisSuperkEff1RmuCCQEoth);
      if(fparm->Exists(conv::kNuisSuperkEff1RmuCCnonQE)) errminosminus_skdet_ccnqe = fparm->MinusError(conv::kNuisSuperkEff1RmuCCnonQE);
      if(fparm->Exists(conv::kNuisSuperkEff1RmuNC)) errminosminus_skdet_nc = fparm->MinusError(conv::kNuisSuperkEff1RmuNC);
      if(fparm->Exists(conv::kNuisSuperkEff1RmuNueCC)) errminosminus_skdet_nuecc = fparm->MinusError(conv::kNuisSuperkEff1RmuNueCC);
      if(fparm->Exists(conv::kNuisFSI000)) errminosminus_fsi_000 = fparm->MinusError(conv::kNuisFSI000);
      if(fparm->Exists(conv::kNuisFSI001)) errminosminus_fsi_001 = fparm->MinusError(conv::kNuisFSI001);
      if(fparm->Exists(conv::kNuisFSI002)) errminosminus_fsi_002 = fparm->MinusError(conv::kNuisFSI002);
      if(fparm->Exists(conv::kNuisFSI003)) errminosminus_fsi_003 = fparm->MinusError(conv::kNuisFSI003);
      if(fparm->Exists(conv::kNuisFSI004)) errminosminus_fsi_004 = fparm->MinusError(conv::kNuisFSI004);
      if(fparm->Exists(conv::kNuisFSI005)) errminosminus_fsi_005 = fparm->MinusError(conv::kNuisFSI005);
    }// minos errors set?
    //get the migrad errors
    if(hyp == kHypOsc3F) {
      errmigrad_mixing12 = fparm->MigradError(kPhypMixing12);
      errmigrad_mixing13 = fparm->MigradError(kPhypMixing13);
      errmigrad_mixing23 = fparm->MigradError(kPhypMixing23);
      errmigrad_deltacp  = fparm->MigradError(kPhypDeltaCP);
      errmigrad_dm12sq   = fparm->MigradError(kPhypDm12sq);
      errmigrad_dm23sq   = fparm->MigradError(kPhypDm23sq);
    }
    else if(hyp == kHypNuMuDisapp2F) {
      errmigrad_mixing23 = fparm->MigradError(kPhypMixing23);
      errmigrad_dm23sq   = fparm->MigradError(kPhypDm23sq);
    }
    if(fparm->Exists(conv::kNuisUnBANFF00)) errmigrad_unbanff_00 = fparm->MigradError(conv::kNuisUnBANFF00);
    if(fparm->Exists(conv::kNuisUnBANFF01)) errmigrad_unbanff_01 = fparm->MigradError(conv::kNuisUnBANFF01);
    if(fparm->Exists(conv::kNuisUnBANFF02)) errmigrad_unbanff_02 = fparm->MigradError(conv::kNuisUnBANFF02);
    if(fparm->Exists(conv::kNuisUnBANFF03)) errmigrad_unbanff_03 = fparm->MigradError(conv::kNuisUnBANFF03);
    if(fparm->Exists(conv::kNuisUnBANFF04)) errmigrad_unbanff_04 = fparm->MigradError(conv::kNuisUnBANFF04);
    if(fparm->Exists(conv::kNuisUnBANFF05)) errmigrad_unbanff_05 = fparm->MigradError(conv::kNuisUnBANFF05);
    if(fparm->Exists(conv::kNuisUnBANFF06)) errmigrad_unbanff_06 = fparm->MigradError(conv::kNuisUnBANFF06);
    if(fparm->Exists(conv::kNuisUnBANFF07)) errmigrad_unbanff_07 = fparm->MigradError(conv::kNuisUnBANFF07);
    if(fparm->Exists(conv::kNuisUnBANFF08)) errmigrad_unbanff_08 = fparm->MigradError(conv::kNuisUnBANFF08);
    if(fparm->Exists(conv::kNuisUnBANFF09)) errmigrad_unbanff_09 = fparm->MigradError(conv::kNuisUnBANFF09);
    if(fparm->Exists(conv::kNuisUnBANFF10)) errmigrad_unbanff_10 = fparm->MigradError(conv::kNuisUnBANFF10);
    if(fparm->Exists(conv::kNuisUnBANFF11)) errmigrad_unbanff_11 = fparm->MigradError(conv::kNuisUnBANFF11);
    if(fparm->Exists(conv::kNuisUnBANFF12)) errmigrad_unbanff_12 = fparm->MigradError(conv::kNuisUnBANFF12);
    if(fparm->Exists(conv::kNuisUnBANFF13)) errmigrad_unbanff_13 = fparm->MigradError(conv::kNuisUnBANFF13);
    if(fparm->Exists(conv::kNuisUnBANFF14)) errmigrad_unbanff_14 = fparm->MigradError(conv::kNuisUnBANFF14);
    if(fparm->Exists(conv::kNuisUnBANFF15)) errmigrad_unbanff_15 = fparm->MigradError(conv::kNuisUnBANFF15);
    if(fparm->Exists(conv::kNuisUnBANFF16)) errmigrad_unbanff_16 = fparm->MigradError(conv::kNuisUnBANFF16);
    if(fparm->Exists(conv::kNuisUnBANFF17)) errmigrad_unbanff_17 = fparm->MigradError(conv::kNuisUnBANFF17);
    if(fparm->Exists(conv::kNuisUnBANFF18)) errmigrad_unbanff_18 = fparm->MigradError(conv::kNuisUnBANFF18);
    if(fparm->Exists(conv::kNuisUnBANFF19)) errmigrad_unbanff_19 = fparm->MigradError(conv::kNuisUnBANFF19);
    if(fparm->Exists(conv::kNuisUnBANFF20)) errmigrad_unbanff_20 = fparm->MigradError(conv::kNuisUnBANFF20);
    if(fparm->Exists(conv::kNuisUnBANFF21)) errmigrad_unbanff_21 = fparm->MigradError(conv::kNuisUnBANFF21);
    if(fparm->Exists(conv::kNuisUnBANFF22)) errmigrad_unbanff_22 = fparm->MigradError(conv::kNuisUnBANFF22);
    if(fparm->Exists(conv::kNuisBANFF00)) errmigrad_banff_00 = fparm->MigradError(conv::kNuisBANFF00);
    if(fparm->Exists(conv::kNuisBANFF01)) errmigrad_banff_01 = fparm->MigradError(conv::kNuisBANFF01);
    if(fparm->Exists(conv::kNuisBANFF02)) errmigrad_banff_02 = fparm->MigradError(conv::kNuisBANFF02);
    if(fparm->Exists(conv::kNuisBANFF03)) errmigrad_banff_03 = fparm->MigradError(conv::kNuisBANFF03);
    if(fparm->Exists(conv::kNuisBANFF04)) errmigrad_banff_04 = fparm->MigradError(conv::kNuisBANFF04);
    if(fparm->Exists(conv::kNuisBANFF05)) errmigrad_banff_05 = fparm->MigradError(conv::kNuisBANFF05);
    if(fparm->Exists(conv::kNuisBANFF06)) errmigrad_banff_06 = fparm->MigradError(conv::kNuisBANFF06);
    if(fparm->Exists(conv::kNuisBANFF07)) errmigrad_banff_07 = fparm->MigradError(conv::kNuisBANFF07);
    if(fparm->Exists(conv::kNuisBANFF08)) errmigrad_banff_08 = fparm->MigradError(conv::kNuisBANFF08);
    if(fparm->Exists(conv::kNuisBANFF09)) errmigrad_banff_09 = fparm->MigradError(conv::kNuisBANFF09);
    if(fparm->Exists(conv::kNuisBANFF10)) errmigrad_banff_10 = fparm->MigradError(conv::kNuisBANFF10);
    if(fparm->Exists(conv::kNuisBANFF11)) errmigrad_banff_11 = fparm->MigradError(conv::kNuisBANFF11);
    if(fparm->Exists(conv::kNuisBANFF12)) errmigrad_banff_12 = fparm->MigradError(conv::kNuisBANFF12);
    if(fparm->Exists(conv::kNuisBANFF13)) errmigrad_banff_13 = fparm->MigradError(conv::kNuisBANFF13);
    if(fparm->Exists(conv::kNuisBANFF14)) errmigrad_banff_14 = fparm->MigradError(conv::kNuisBANFF14);
    if(fparm->Exists(conv::kNuisBANFF15)) errmigrad_banff_15 = fparm->MigradError(conv::kNuisBANFF15);
    if(fparm->Exists(conv::kNuisBANFFMAQE)) errmigrad_banff_maqe = fparm->MigradError(conv::kNuisBANFFMAQE);
    if(fparm->Exists(conv::kNuisBANFFMARES)) errmigrad_banff_mares = fparm->MigradError(conv::kNuisBANFFMARES);
    if(fparm->Exists(conv::kNuisBANFFCCQEE0)) errmigrad_banff_ccqeE0 = fparm->MigradError(conv::kNuisBANFFCCQEE0);
    if(fparm->Exists(conv::kNuisBANFFCCQEE1)) errmigrad_banff_ccqeE1 = fparm->MigradError(conv::kNuisBANFFCCQEE1);
    if(fparm->Exists(conv::kNuisBANFFCCQEE2)) errmigrad_banff_ccqeE2 = fparm->MigradError(conv::kNuisBANFFCCQEE2);
    if(fparm->Exists(conv::kNuisBANFFCC1piE0)) errmigrad_banff_cc1piE0 = fparm->MigradError(conv::kNuisBANFFCC1piE0);
    if(fparm->Exists(conv::kNuisBANFFCC1piE1)) errmigrad_banff_cc1piE1 = fparm->MigradError(conv::kNuisBANFFCC1piE1);
    if(fparm->Exists(conv::kNuisXSecFermiMom)) errmigrad_xsec_pf = fparm->MigradError(conv::kNuisXSecFermiMom);
    if(fparm->Exists(conv::kNuisXSecSpecFunc)) errmigrad_xsec_sf = fparm->MigradError(conv::kNuisXSecSpecFunc);
    if(fparm->Exists(conv::kNuisXSecBindE)) errmigrad_xsec_binde = fparm->MigradError(conv::kNuisXSecBindE);
    if(fparm->Exists(conv::kNuisXSecWshape)) errmigrad_xsec_wshape = fparm->MigradError(conv::kNuisXSecWshape);
    if(fparm->Exists(conv::kNuisXSecShapeCCoth)) errmigrad_xsec_ccothshape = fparm->MigradError(conv::kNuisXSecShapeCCoth);
    if(fparm->Exists(conv::kNuisXSecPiLessDelta)) errmigrad_xsec_pilessdelta = fparm->MigradError(conv::kNuisXSecPiLessDelta);
    if(fparm->Exists(conv::kNuisXSecENuShape)) errmigrad_xsec_enushape = fparm->MigradError(conv::kNuisXSecENuShape);
    if(fparm->Exists(conv::kNuisXSecNormCCcoh)) errmigrad_xsec_normcccoh = fparm->MigradError(conv::kNuisXSecNormCCcoh);
    if(fparm->Exists(conv::kNuisXSecNormNC1piC)) errmigrad_xsec_normnc1piC = fparm->MigradError(conv::kNuisXSecNormNC1piC);
    if(fparm->Exists(conv::kNuisXSecNormNCcoh)) errmigrad_xsec_normnccoh = fparm->MigradError(conv::kNuisXSecNormNCcoh);
    if(fparm->Exists(conv::kNuisXSecNormNCoth)) errmigrad_xsec_normncoth = fparm->MigradError(conv::kNuisXSecNormNCoth);
    if(fparm->Exists(conv::kNuisXSecNormNuEToNuMu)) errmigrad_xsec_normnuetonumu = fparm->MigradError(conv::kNuisXSecNormNuEToNuMu);
    if(fparm->Exists(conv::kNuisXSecNormAntiNuToNu)) errmigrad_xsec_normantinutonu = fparm->MigradError(conv::kNuisXSecNormAntiNuToNu);
    if(fparm->Exists(conv::kNuisSuperkEnergyScale)) errmigrad_skdet_energyscale = fparm->MigradError(conv::kNuisSuperkEnergyScale);
    if(fparm->Exists(conv::kNuisSuperkEff1RmuCCQErc0u)) errmigrad_skdet_ccqerc0 = fparm->MigradError(conv::kNuisSuperkEff1RmuCCQErc0u);
    if(fparm->Exists(conv::kNuisSuperkEff1RmuCCQErc1u)) errmigrad_skdet_ccqerc1 = fparm->MigradError(conv::kNuisSuperkEff1RmuCCQErc1u);
    if(fparm->Exists(conv::kNuisSuperkEff1RmuCCQErc2u)) errmigrad_skdet_ccqerc2 = fparm->MigradError(conv::kNuisSuperkEff1RmuCCQErc2u);
    if(fparm->Exists(conv::kNuisSuperkEff1RmuCCQEoth)) errmigrad_skdet_ccqeoth = fparm->MigradError(conv::kNuisSuperkEff1RmuCCQEoth);
    if(fparm->Exists(conv::kNuisSuperkEff1RmuCCnonQE)) errmigrad_skdet_ccnqe = fparm->MigradError(conv::kNuisSuperkEff1RmuCCnonQE);
    if(fparm->Exists(conv::kNuisSuperkEff1RmuNC)) errmigrad_skdet_nc = fparm->MigradError(conv::kNuisSuperkEff1RmuNC);
    if(fparm->Exists(conv::kNuisSuperkEff1RmuNueCC)) errmigrad_skdet_nuecc = fparm->MigradError(conv::kNuisSuperkEff1RmuNueCC);
    if(fparm->Exists(conv::kNuisFSI000)) errmigrad_fsi_000 = fparm->MigradError(conv::kNuisFSI000);
    if(fparm->Exists(conv::kNuisFSI001)) errmigrad_fsi_001 = fparm->MigradError(conv::kNuisFSI001);
    if(fparm->Exists(conv::kNuisFSI002)) errmigrad_fsi_002 = fparm->MigradError(conv::kNuisFSI002);
    if(fparm->Exists(conv::kNuisFSI003)) errmigrad_fsi_003 = fparm->MigradError(conv::kNuisFSI003);
    if(fparm->Exists(conv::kNuisFSI004)) errmigrad_fsi_004 = fparm->MigradError(conv::kNuisFSI004);
    if(fparm->Exists(conv::kNuisFSI005)) errmigrad_fsi_005 = fparm->MigradError(conv::kNuisFSI005);

    //get the hesse errors
    if(fparm->HesseErrorsSet()) {
      if(hyp == kHypOsc3F) {
	errhesse_mixing12 = fparm->HesseError(kPhypMixing12);
	errhesse_mixing13 = fparm->HesseError(kPhypMixing13);
	errhesse_mixing23 = fparm->HesseError(kPhypMixing23);
	errhesse_deltacp  = fparm->HesseError(kPhypDeltaCP);
	errhesse_dm12sq   = fparm->HesseError(kPhypDm12sq);
	errhesse_dm23sq   = fparm->HesseError(kPhypDm23sq);
      }
      else if(hyp == kHypNuMuDisapp2F) {
	errhesse_mixing23 = fparm->HesseError(kPhypMixing23);
	errhesse_dm23sq   = fparm->HesseError(kPhypDm23sq);
      }
      if(fparm->Exists(conv::kNuisUnBANFF00)) errhesse_unbanff_00 = fparm->HesseError(conv::kNuisUnBANFF00);
      if(fparm->Exists(conv::kNuisUnBANFF01)) errhesse_unbanff_01 = fparm->HesseError(conv::kNuisUnBANFF01);
      if(fparm->Exists(conv::kNuisUnBANFF02)) errhesse_unbanff_02 = fparm->HesseError(conv::kNuisUnBANFF02);
      if(fparm->Exists(conv::kNuisUnBANFF03)) errhesse_unbanff_03 = fparm->HesseError(conv::kNuisUnBANFF03);
      if(fparm->Exists(conv::kNuisUnBANFF04)) errhesse_unbanff_04 = fparm->HesseError(conv::kNuisUnBANFF04);
      if(fparm->Exists(conv::kNuisUnBANFF05)) errhesse_unbanff_05 = fparm->HesseError(conv::kNuisUnBANFF05);
      if(fparm->Exists(conv::kNuisUnBANFF06)) errhesse_unbanff_06 = fparm->HesseError(conv::kNuisUnBANFF06);
      if(fparm->Exists(conv::kNuisUnBANFF07)) errhesse_unbanff_07 = fparm->HesseError(conv::kNuisUnBANFF07);
      if(fparm->Exists(conv::kNuisUnBANFF08)) errhesse_unbanff_08 = fparm->HesseError(conv::kNuisUnBANFF08);
      if(fparm->Exists(conv::kNuisUnBANFF09)) errhesse_unbanff_09 = fparm->HesseError(conv::kNuisUnBANFF09);
      if(fparm->Exists(conv::kNuisUnBANFF10)) errhesse_unbanff_10 = fparm->HesseError(conv::kNuisUnBANFF10);
      if(fparm->Exists(conv::kNuisUnBANFF11)) errhesse_unbanff_11 = fparm->HesseError(conv::kNuisUnBANFF11);
      if(fparm->Exists(conv::kNuisUnBANFF12)) errhesse_unbanff_12 = fparm->HesseError(conv::kNuisUnBANFF12);
      if(fparm->Exists(conv::kNuisUnBANFF13)) errhesse_unbanff_13 = fparm->HesseError(conv::kNuisUnBANFF13);
      if(fparm->Exists(conv::kNuisUnBANFF14)) errhesse_unbanff_14 = fparm->HesseError(conv::kNuisUnBANFF14);
      if(fparm->Exists(conv::kNuisUnBANFF15)) errhesse_unbanff_15 = fparm->HesseError(conv::kNuisUnBANFF15);
      if(fparm->Exists(conv::kNuisUnBANFF16)) errhesse_unbanff_16 = fparm->HesseError(conv::kNuisUnBANFF16);
      if(fparm->Exists(conv::kNuisUnBANFF17)) errhesse_unbanff_17 = fparm->HesseError(conv::kNuisUnBANFF17);
      if(fparm->Exists(conv::kNuisUnBANFF18)) errhesse_unbanff_18 = fparm->HesseError(conv::kNuisUnBANFF18);
      if(fparm->Exists(conv::kNuisUnBANFF19)) errhesse_unbanff_19 = fparm->HesseError(conv::kNuisUnBANFF19);
      if(fparm->Exists(conv::kNuisUnBANFF20)) errhesse_unbanff_20 = fparm->HesseError(conv::kNuisUnBANFF20);
      if(fparm->Exists(conv::kNuisUnBANFF21)) errhesse_unbanff_21 = fparm->HesseError(conv::kNuisUnBANFF21);
      if(fparm->Exists(conv::kNuisUnBANFF22)) errhesse_unbanff_22 = fparm->HesseError(conv::kNuisUnBANFF22);
      if(fparm->Exists(conv::kNuisBANFF00)) errhesse_banff_00 = fparm->HesseError(conv::kNuisBANFF00);
      if(fparm->Exists(conv::kNuisBANFF01)) errhesse_banff_01 = fparm->HesseError(conv::kNuisBANFF01);
      if(fparm->Exists(conv::kNuisBANFF02)) errhesse_banff_02 = fparm->HesseError(conv::kNuisBANFF02);
      if(fparm->Exists(conv::kNuisBANFF03)) errhesse_banff_03 = fparm->HesseError(conv::kNuisBANFF03);
      if(fparm->Exists(conv::kNuisBANFF04)) errhesse_banff_04 = fparm->HesseError(conv::kNuisBANFF04);
      if(fparm->Exists(conv::kNuisBANFF05)) errhesse_banff_05 = fparm->HesseError(conv::kNuisBANFF05);
      if(fparm->Exists(conv::kNuisBANFF06)) errhesse_banff_06 = fparm->HesseError(conv::kNuisBANFF06);
      if(fparm->Exists(conv::kNuisBANFF07)) errhesse_banff_07 = fparm->HesseError(conv::kNuisBANFF07);
      if(fparm->Exists(conv::kNuisBANFF08)) errhesse_banff_08 = fparm->HesseError(conv::kNuisBANFF08);
      if(fparm->Exists(conv::kNuisBANFF09)) errhesse_banff_09 = fparm->HesseError(conv::kNuisBANFF09);
      if(fparm->Exists(conv::kNuisBANFF10)) errhesse_banff_10 = fparm->HesseError(conv::kNuisBANFF10);
      if(fparm->Exists(conv::kNuisBANFF11)) errhesse_banff_11 = fparm->HesseError(conv::kNuisBANFF11);
      if(fparm->Exists(conv::kNuisBANFF12)) errhesse_banff_12 = fparm->HesseError(conv::kNuisBANFF12);
      if(fparm->Exists(conv::kNuisBANFF13)) errhesse_banff_13 = fparm->HesseError(conv::kNuisBANFF13);
      if(fparm->Exists(conv::kNuisBANFF14)) errhesse_banff_14 = fparm->HesseError(conv::kNuisBANFF14);
      if(fparm->Exists(conv::kNuisBANFF15)) errhesse_banff_15 = fparm->HesseError(conv::kNuisBANFF15);
      if(fparm->Exists(conv::kNuisBANFFMAQE)) errhesse_banff_maqe = fparm->HesseError(conv::kNuisBANFFMAQE);
      if(fparm->Exists(conv::kNuisBANFFMARES)) errhesse_banff_mares = fparm->HesseError(conv::kNuisBANFFMARES);
      if(fparm->Exists(conv::kNuisBANFFCCQEE0)) errhesse_banff_ccqeE0 = fparm->HesseError(conv::kNuisBANFFCCQEE0);
      if(fparm->Exists(conv::kNuisBANFFCCQEE1)) errhesse_banff_ccqeE1 = fparm->HesseError(conv::kNuisBANFFCCQEE1);
      if(fparm->Exists(conv::kNuisBANFFCCQEE2)) errhesse_banff_ccqeE2 = fparm->HesseError(conv::kNuisBANFFCCQEE2);
      if(fparm->Exists(conv::kNuisBANFFCC1piE0)) errhesse_banff_cc1piE0 = fparm->HesseError(conv::kNuisBANFFCC1piE0);
      if(fparm->Exists(conv::kNuisBANFFCC1piE1)) errhesse_banff_cc1piE1 = fparm->HesseError(conv::kNuisBANFFCC1piE1);
      if(fparm->Exists(conv::kNuisXSecFermiMom)) errhesse_xsec_pf = fparm->HesseError(conv::kNuisXSecFermiMom);
      if(fparm->Exists(conv::kNuisXSecSpecFunc)) errhesse_xsec_sf = fparm->HesseError(conv::kNuisXSecSpecFunc);
      if(fparm->Exists(conv::kNuisXSecBindE)) errhesse_xsec_binde = fparm->HesseError(conv::kNuisXSecBindE);
      if(fparm->Exists(conv::kNuisXSecWshape)) errhesse_xsec_wshape = fparm->HesseError(conv::kNuisXSecWshape);
      if(fparm->Exists(conv::kNuisXSecShapeCCoth)) errhesse_xsec_ccothshape = fparm->HesseError(conv::kNuisXSecShapeCCoth);
      if(fparm->Exists(conv::kNuisXSecPiLessDelta)) errhesse_xsec_pilessdelta = fparm->HesseError(conv::kNuisXSecPiLessDelta);
      if(fparm->Exists(conv::kNuisXSecENuShape)) errhesse_xsec_enushape = fparm->HesseError(conv::kNuisXSecENuShape);
      if(fparm->Exists(conv::kNuisXSecNormCCcoh)) errhesse_xsec_normcccoh = fparm->HesseError(conv::kNuisXSecNormCCcoh);
      if(fparm->Exists(conv::kNuisXSecNormNC1piC)) errhesse_xsec_normnc1piC = fparm->HesseError(conv::kNuisXSecNormNC1piC);
      if(fparm->Exists(conv::kNuisXSecNormNCcoh)) errhesse_xsec_normnccoh = fparm->HesseError(conv::kNuisXSecNormNCcoh);
      if(fparm->Exists(conv::kNuisXSecNormNCoth)) errhesse_xsec_normncoth = fparm->HesseError(conv::kNuisXSecNormNCoth);
      if(fparm->Exists(conv::kNuisXSecNormNuEToNuMu)) errhesse_xsec_normnuetonumu = fparm->HesseError(conv::kNuisXSecNormNuEToNuMu);
      if(fparm->Exists(conv::kNuisXSecNormAntiNuToNu)) errhesse_xsec_normantinutonu = fparm->HesseError(conv::kNuisXSecNormAntiNuToNu);
      if(fparm->Exists(conv::kNuisSuperkEnergyScale)) errhesse_skdet_energyscale = fparm->HesseError(conv::kNuisSuperkEnergyScale);
      if(fparm->Exists(conv::kNuisSuperkEff1RmuCCQErc0u)) errhesse_skdet_ccqerc0 = fparm->HesseError(conv::kNuisSuperkEff1RmuCCQErc0u);
      if(fparm->Exists(conv::kNuisSuperkEff1RmuCCQErc1u)) errhesse_skdet_ccqerc1 = fparm->HesseError(conv::kNuisSuperkEff1RmuCCQErc1u);
      if(fparm->Exists(conv::kNuisSuperkEff1RmuCCQErc2u)) errhesse_skdet_ccqerc2 = fparm->HesseError(conv::kNuisSuperkEff1RmuCCQErc2u);
      if(fparm->Exists(conv::kNuisSuperkEff1RmuCCQEoth)) errhesse_skdet_ccqeoth = fparm->HesseError(conv::kNuisSuperkEff1RmuCCQEoth);
      if(fparm->Exists(conv::kNuisSuperkEff1RmuCCnonQE)) errhesse_skdet_ccnqe = fparm->HesseError(conv::kNuisSuperkEff1RmuCCnonQE);
      if(fparm->Exists(conv::kNuisSuperkEff1RmuNC)) errhesse_skdet_nc = fparm->HesseError(conv::kNuisSuperkEff1RmuNC);
      if(fparm->Exists(conv::kNuisSuperkEff1RmuNueCC)) errhesse_skdet_nuecc = fparm->HesseError(conv::kNuisSuperkEff1RmuNueCC);
      if(fparm->Exists(conv::kNuisFSI000)) errhesse_fsi_000 = fparm->HesseError(conv::kNuisFSI000);
      if(fparm->Exists(conv::kNuisFSI001)) errhesse_fsi_001 = fparm->HesseError(conv::kNuisFSI001);
      if(fparm->Exists(conv::kNuisFSI002)) errhesse_fsi_002 = fparm->HesseError(conv::kNuisFSI002);
      if(fparm->Exists(conv::kNuisFSI003)) errhesse_fsi_003 = fparm->HesseError(conv::kNuisFSI003);
      if(fparm->Exists(conv::kNuisFSI004)) errhesse_fsi_004 = fparm->HesseError(conv::kNuisFSI004);
      if(fparm->Exists(conv::kNuisFSI005)) errhesse_fsi_005 = fparm->HesseError(conv::kNuisFSI005);
    }//hesse errors set?

    ///////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////
    // Fill tree
    ///////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////
    nfitnt->Fill();

    string output_option = "";
    if ( store_toy_expts          ) { output_option += "Observation";     }
    if ( store_bestfit_spectrum   ) { output_option += "+BestFitSpectra"; }
    if ( store_chisq_surf         ) { output_option += "+ChisqSurface";   }
    fitter->Analyze(output_option.c_str());

    if ( store_withsyst_chisq_surf) {

      //set up binning
      double mixing23_min  = kMinChisqSurfaceSinsq2Theta23;
      double mixing23_max  = kMaxChisqSurfaceSinsq2Theta23;
      double mixing23_step = fparm->StepSize(conv::kPhypMixing23);
      double dm23sq_min    = kMinChisqSurfaceDm23sq;
      double dm23sq_max    = kMaxChisqSurfaceDm23sq;
      double dm23sq_step   = fparm->StepSize(conv::kPhypDm23sq);

      // make sure that min and max are at the bin centre
      int mixing23_nbins = 1 + (int) ((mixing23_max - mixing23_min) / mixing23_step);
      int dm23sq_nbins   = 1 + (int) ((dm23sq_max -   dm23sq_min )  / dm23sq_step );

      //FIRST
      // Use the 7 best-fit systematics, the N-7 nominal systematics, loop over oscillation params and calculate chisq (this is the chisq surface
      // without systematics). The chisq at the minimum is the chisq_min. Calculate the deltachisq surface.

      TH2D * chisq_surface = new TH2D("chisq_surface","#chi^{2} surface",
		mixing23_nbins, mixing23_min - 0.5*mixing23_step, mixing23_max + 0.5*mixing23_step,
		dm23sq_nbins, dm23sq_min - 0.5*dm23sq_step, dm23sq_max + 0.5*dm23sq_step);
      chisq_surface->SetDirectory(0);
      chisq_surface->GetXaxis()->SetRangeUser(0.,1.);
      chisq_surface->GetXaxis()->SetTitle("sin^{2}(2#theta_{23})");
      chisq_surface->GetYaxis()->SetTitle("|#Deltam_{32}^{2}| (eV^{2}/c^{4})");
      chisq_surface->GetYaxis()->SetTitleOffset(1);
      chisq_surface->GetYaxis()->SetTitleSize(0.05);
      chisq_surface->GetXaxis()->SetTitleSize(0.05);

      TH2D * dchisq_surface = new TH2D("dchisq_surface","#Delta#chi^{2} surface",
		mixing23_nbins, mixing23_min - 0.5*mixing23_step, mixing23_max + 0.5*mixing23_step,
		dm23sq_nbins, dm23sq_min - 0.5*dm23sq_step, dm23sq_max + 0.5*dm23sq_step);
      dchisq_surface->SetDirectory(0);
      dchisq_surface->GetXaxis()->SetRangeUser(0.,1.);
      dchisq_surface->GetXaxis()->SetTitle("sin^{2}(2#theta_{23})");
      dchisq_surface->GetYaxis()->SetTitle("|#Deltam_{32}^{2}| (eV^{2}/c^{4})");
      dchisq_surface->GetYaxis()->SetTitleOffset(1);
      dchisq_surface->GetYaxis()->SetTitleSize(0.05);
      dchisq_surface->GetXaxis()->SetTitleSize(0.05);

      //make sure that the systematics are set at their bestfit value
      fparm->CopyBestFitToCurrValues();

      //Get the chisq surface
      double nfits_done = 0;
      double nfits_all  = chisq_surface->GetNbinsX() * chisq_surface->GetNbinsY();
      for(int i = 1; i <= chisq_surface->GetNbinsX(); i++) {
	for(int j = 1; j <= chisq_surface->GetNbinsY(); j++) {
	  double mixing23 = chisq_surface->GetXaxis()->GetBinCenter(i);
	  double dm23sq   = chisq_surface->GetYaxis()->GetBinCenter(j);
	  cout << "> Making chisq surface, fixing Dm23sq := " << dm23sq << " eV^2 and mixing23 := " << mixing23 << endl;
	  fparm->Update( conv::kPhypMixing23, mixing23 );
	  fparm->Update( conv::kPhypDm23sq,   dm23sq   );
	  double chisq = fit_inputs->Chisquare();
	  chisq_surface->SetBinContent(i,j,chisq);
	  nfits_done += 1;
	}
	cout << "> ..... " << 100*nfits_done/nfits_all << "% done" << endl;
      }

      //Get the min(chisq)
      double min_chisq = chisq_surface->GetBinContent(chisq_surface->GetMinimumBin());
      
      //Fill delta(chisq) surface
      for(int i = 1; i <= chisq_surface->GetNbinsX(); i++) {
	for(int j = 1; j <= chisq_surface->GetNbinsY(); j++) {
	  double chi = chisq_surface->GetBinContent(i,j);
	  dchisq_surface->SetBinContent(i,j, chi - min_chisq);
	}
      }

      //save
      chisq_surface  -> SetName(Form("chisq_surface_1Rmu_%d",iexperiment));
      dchisq_surface -> SetName(Form("dchisq_surface_1Rmu_%d",iexperiment));
      f.cd();
      chisq_surface ->Write();
      dchisq_surface->Write();

      //NEXT
      //Loop over oscillation param bins, fix oscillation params and fit for the 7 systematics (*note* that was N before).
      //This will give lower chisq than before so you can calculate by how much you need to correct the deltachisq surface
      /* Turn off for now.
      TH2D * chisq_surface_syst = new TH2D("chisq_surface_syst","#chi^{2} surface",
		mixing23_nbins, mixing23_min - 0.5*mixing23_step, mixing23_max + 0.5*mixing23_step,
		dm23sq_nbins, dm23sq_min - 0.5*dm23sq_step, dm23sq_max + 0.5*dm23sq_step);
      chisq_surface_syst->SetDirectory(0);
      chisq_surface_syst->GetXaxis()->SetRangeUser(0.,1.);
      chisq_surface_syst->GetXaxis()->SetTitle("sin^{2}(2#theta_{23})");
      chisq_surface_syst->GetYaxis()->SetTitle("|#Deltam_{32}^{2}| (eV^{2}/c^{4})");
      chisq_surface_syst->GetYaxis()->SetTitleOffset(1);
      chisq_surface_syst->GetYaxis()->SetTitleSize(0.05);
      chisq_surface_syst->GetXaxis()->SetTitleSize(0.05);

      TH2D * dchisq_surface_syst = new TH2D("dchisq_surface_syst","#Delta#chi^{2} surface",
		mixing23_nbins, mixing23_min - 0.5*mixing23_step, mixing23_max + 0.5*mixing23_step,
		dm23sq_nbins, dm23sq_min - 0.5*dm23sq_step, dm23sq_max + 0.5*dm23sq_step);
      dchisq_surface_syst->SetDirectory(0);
      dchisq_surface_syst->GetXaxis()->SetRangeUser(0.,1.);
      dchisq_surface_syst->GetXaxis()->SetTitle("sin^{2}(2#theta_{23})");
      dchisq_surface_syst->GetYaxis()->SetTitle("|#Deltam_{32}^{2}| (eV^{2}/c^{4})");
      dchisq_surface_syst->GetYaxis()->SetTitleOffset(1);
      dchisq_surface_syst->GetYaxis()->SetTitleSize(0.05);
      dchisq_surface_syst->GetXaxis()->SetTitleSize(0.05);

      // Get the chisq surface with systematics (fit for systematics at every point)
      nfits_done = 0;
      for(int i = 1; i <= chisq_surface_syst->GetNbinsX(); i++) {
	for(int j = 1; j <= chisq_surface_syst->GetNbinsY(); j++) {
	  double mixing23 = chisq_surface_syst->GetXaxis()->GetBinCenter(i);
	  double dm23sq   = chisq_surface_syst->GetYaxis()->GetBinCenter(j);
	  cout << "> Fit for systematics, fixing Dm23sq := " << dm23sq << " eV^2 and mixing23 := " << mixing23 << endl;
	  //          cout << *fparm;
	  //          cout << "Resetting current to the initial values" << endl;
	  fparm->ResetCurrToInitValues();
	  //          cout << *fparm;
	  //          cout << "Free-ing 23-params" << endl;
	  fparm->FreeParameter(conv::kPhypMixing23);
	  fparm->FreeParameter(conv::kPhypDm23sq);
	  //          cout << *fparm;
	  //          cout << "Setting new initial value for 23-params" << endl;
	  fparm->SetInitValue( conv::kPhypMixing23, mixing23 );
	  fparm->SetInitValue( conv::kPhypDm23sq,   dm23sq   );
	  //          cout << *fparm;
	  //          cout << "Fixing 23-params" << endl;
	  fparm->FixParameter(conv::kPhypMixing23);
	  fparm->FixParameter(conv::kPhypDm23sq);
	  //          cout << *fparm;
	  //for some reason, need to make a new fitter when changing the FitParamList
	  FitInputs * fit_inputs2 = new FitInputs(skobs,skexp,fparm);
	  FitterABC * fitter2 = new Superk1RmuErecoFitterMinuit2(fit_inputs,fit_option);
	  fitter2->DoTheFit();
	  double chisq = fit_inputs2->Chisquare();
	  chisq_surface_syst->SetBinContent(i,j,chisq);
	  dchisq_surface_syst->SetBinContent(i,j, chisq - min_chisq);
	  //if(fparm->CurrValue(conv::kPhypMixing23) != mixing23 || fparm->CurrValue(conv::kPhypDm23sq) != dm23sq) {}
	  delete fitter2;
	  delete fit_inputs2;
          nfits_done += 1;
          cout << "> ..... " << 100*nfits_done/nfits_all << "% done" << endl;
	}//j
      }//i

      //Check the min(chisq)
      double min_chisq_syst = chisq_surface_syst->GetBinContent(chisq_surface_syst->GetMinimumBin());
      cout << "Minima of the 2 chisq surface should agree: " << min_chisq << "\t" << min_chisq_syst << endl;

      //save
      chisq_surface_syst  -> SetName(Form("chisq_surface_syst_1Rmu_%d",iexperiment));
      dchisq_surface_syst -> SetName(Form("dchisq_surface_syst_1Rmu_%d",iexperiment));
      f.cd();
      chisq_surface_syst ->Write();
      dchisq_surface_syst->Write();

      //FINALLY
      //correct the dchisq surface, using the dchisq surface with syst
      TH2D * dchisq_surface_corr = new TH2D("dchisq_surface_corr","#Delta#chi^{2} surface",
		mixing23_nbins, mixing23_min - 0.5*mixing23_step, mixing23_max + 0.5*mixing23_step,
		dm23sq_nbins, dm23sq_min - 0.5*dm23sq_step, dm23sq_max + 0.5*dm23sq_step);
      dchisq_surface_corr->SetDirectory(0);
      dchisq_surface_corr->GetXaxis()->SetRangeUser(0.,1.);
      dchisq_surface_corr->GetXaxis()->SetTitle("sin^{2}(2#theta_{23})");
      dchisq_surface_corr->GetYaxis()->SetTitle("|#Deltam_{32}^{2}| (eV^{2}/c^{4})");
      dchisq_surface_corr->GetYaxis()->SetTitleOffset(1);
      dchisq_surface_corr->GetYaxis()->SetTitleSize(0.05);
      dchisq_surface_corr->GetXaxis()->SetTitleSize(0.05);

      // what to actually do?

      //save
      dchisq_surface_corr -> SetName(Form("dchisq_surface_corr_1Rmu_%d",iexperiment));
      f.cd();
      dchisq_surface_corr->Write();
      */
    }    

    // Store toy-MC expt (if asked)
    if(store_toy_expts) {
      TH1D * h = (TH1D*) fitter->GetResults().FindObject("data_1Rmu");
      TH1D * toyexpspectrum = new TH1D(*h);
      f.cd();
      toyexpspectrum->SetName(Form("toyexp_ereco_spectrum_1Rmu_%d",iexperiment));
      toyexpspectrum->Write();
    }

    // Store best-fit spectrum (if asked)
    if(store_bestfit_spectrum) {
      TH1D * h = (TH1D*) fitter->GetResults().FindObject("bestfit_1Rmu");
      TH1D * bestfitspectrum = new TH1D(*h);
      bestfitspectrum->SetName(Form("bestfit_ereco_spectrum_1Rmu_%d",iexperiment));
      f.cd();
      bestfitspectrum->Write();
    }
    // Store chisq surface (if asked)
    if(store_chisq_surf) {
      TH2D * h2 = (TH2D*) fitter->GetResults().FindObject("chisq_surface_1Rmu");
      if(h2) {
	TH2D * chisqsurf = new TH2D(*h2); 
	chisqsurf->SetName(Form("chisqsurf_%d",iexperiment)); 
	f.cd();
	chisqsurf->Write();
      }
      TH2D * h2a = (TH2D*) fitter->GetResults().FindObject("dchisq_surface_1Rmu");
      if(h2a) {
	TH2D * dchisqsurf = new TH2D(*h2a); 
	dchisqsurf->SetName(Form("dchisqsurf_%d",iexperiment)); 
	f.cd();
	dchisqsurf->Write();
      }
    }
    /* doesn't currently work
    // Store chisq syst surface (if asked)
    if(store_withsyst_chisq_surf) {
      TH2D * h2 = (TH2D*) fitter->GetResults().FindObject("chisq_withsyst_surface_1Rmu");
      TH2D * chisqsystsurf = new TH2D(*h2); 
      chisqsystsurf->SetName(Form("chisqsystsurf_%d",iexperiment)); 
      f.cd();
      chisqsystsurf->Write();
    }
    */

    delete skobs;
    delete fit_inputs;
    delete fitter;
  }//iexperiment

  // Store ntuple with best-fit results
  f.cd();
  nfitnt->Write("nfitnt");

  f.Close();
}
