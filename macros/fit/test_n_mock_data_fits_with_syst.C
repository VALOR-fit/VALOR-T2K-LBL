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
  #include "Superk1RmuErecoFitter.h"
*/

void test_n_mock_data_fits_with_syst
(
 const char * tag                     = "toyexptfitswithsyst",          /* tag for output file */
 int          syst_mode               = 0,                              /* which systematics to tweak? 0 none. 1 all. 2 BANFF. 3 Unconstrained XSec. 4 SK detector. 5 FSI. */
 int          run                     = 0,                              /* run number */
 int          nexperiments            = 1,                            /* number of toy MC experiment fits */
 double       input_mixing12          = 0.86,                           /* 12-mixing, sinq2theta or sinsqtheta depending on user choice, if -999 then randomize */
 double       input_mixing13          = 0.10,                           /* 13-mixing, sinq2theta or sinsqtheta depending on user choice, if -999 then randomize */
 double       input_mixing23          = 1.00,                           /* 23-mixing, sinq2theta or sinsqtheta depending on user choice */
 double       input_deltacp           = 0.00,                           /* deltacp, if -999 then randomize */
 double       input_dm12sq            = 7.58E-5,                        /* Dm^2_{12}, if -999 then randomize */
 double       input_dm23sq            = 2.32E-3,                        /* |Dm^2_{23}| */
 const char * config_file             = "analysis_configuration.xml",   /* analysis configuration file */
 bool         stat_fluctuations       = true,                           /* include statistical fluctuation in toy MC expts? */
 string       fit_param_list          = "nuosc3f_numudisapp_fix4+nuisanceselection", /* fit parameter list for the given physics hypothesis (see run_fitter.C) */
 const char * fit_option              = "minuit",                       /* fitter option (see run_fitter.C) */
 bool         store_chisq_surf        = true,                          /* by default storing only best-fit values, enable this to store the chisq surfaces   */
 bool         store_withsyst_chisq_surf = true,                        /* by default storing only best-fit values, enable this to store the chisq with systematics surfaces   */
 bool         store_toy_expts         = true,                          /* by default storing only best-fit values, enable this to store the toy MC expts     */
 bool         store_bestfit_spectrum  = true                           /* by default storing only best-fit values, enable this to store the best-fit spectra */
 )
{
  using namespace std;
  using namespace t2kana;
  using namespace t2kana::conv;

  // Global options
  Opt::Instance()->gOscMixingParamChoice = kNuMixSinsq2Theta;
  Opt::Instance()->gMassHierarchyAssumed = kNmhNormal;

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
  
  // ***** variables for the branches *******
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
  double bestfit_mixing12, errp_mixing12, errm_mixing12, errpar_mixing12;
  double bestfit_mixing13, errp_mixing13, errm_mixing13, errpar_mixing13;
  double bestfit_mixing23, errp_mixing23, errm_mixing23, errpar_mixing23;
  double bestfit_deltacp, errp_deltacp, errm_deltacp, errpar_deltacp;
  double bestfit_dm12sq, errp_dm12sq, errm_dm12sq, errpar_dm12sq;
  double bestfit_dm23sq, errp_dm23sq, errm_dm23sq, errpar_dm23sq;
  double input_banff_00, bestfit_banff_00, errp_banff_00, errm_banff_00, errpar_banff_00;
  double input_banff_01, bestfit_banff_01, errp_banff_01, errm_banff_01, errpar_banff_01;
  double input_banff_02, bestfit_banff_02, errp_banff_02, errm_banff_02, errpar_banff_02;
  double input_banff_03, bestfit_banff_03, errp_banff_03, errm_banff_03, errpar_banff_03;
  double input_banff_04, bestfit_banff_04, errp_banff_04, errm_banff_04, errpar_banff_04;
  double input_banff_05, bestfit_banff_05, errp_banff_05, errm_banff_05, errpar_banff_05;
  double input_banff_06, bestfit_banff_06, errp_banff_06, errm_banff_06, errpar_banff_06;
  double input_banff_07, bestfit_banff_07, errp_banff_07, errm_banff_07, errpar_banff_07;
  double input_banff_08, bestfit_banff_08, errp_banff_08, errm_banff_08, errpar_banff_08;
  double input_banff_09, bestfit_banff_09, errp_banff_09, errm_banff_09, errpar_banff_09;
  double input_banff_10, bestfit_banff_10, errp_banff_10, errm_banff_10, errpar_banff_10;
  double input_banff_11, bestfit_banff_11, errp_banff_11, errm_banff_11, errpar_banff_11;
  double input_banff_12, bestfit_banff_12, errp_banff_12, errm_banff_12, errpar_banff_12;
  double input_banff_13, bestfit_banff_13, errp_banff_13, errm_banff_13, errpar_banff_13;
  double input_banff_14, bestfit_banff_14, errp_banff_14, errm_banff_14, errpar_banff_14;
  double input_banff_15, bestfit_banff_15, errp_banff_15, errm_banff_15, errpar_banff_15;
  double input_banff_16, bestfit_banff_16, errp_banff_16, errm_banff_16, errpar_banff_16;
  double input_banff_17, bestfit_banff_17, errp_banff_17, errm_banff_17, errpar_banff_17;
  double input_banff_18, bestfit_banff_18, errp_banff_18, errm_banff_18, errpar_banff_18;
  double input_banff_19, bestfit_banff_19, errp_banff_19, errm_banff_19, errpar_banff_19;
  double input_banff_20, bestfit_banff_20, errp_banff_20, errm_banff_20, errpar_banff_20;
  double input_banff_21, bestfit_banff_21, errp_banff_21, errm_banff_21, errpar_banff_21;
  double input_banff_22, bestfit_banff_22, errp_banff_22, errm_banff_22, errpar_banff_22;
  double input_banff_23, bestfit_banff_23, errp_banff_23, errm_banff_23, errpar_banff_23;
  double input_banff_24, bestfit_banff_24, errp_banff_24, errm_banff_24, errpar_banff_24;
  double input_banff_25, bestfit_banff_25, errp_banff_25, errm_banff_25, errpar_banff_25;
  double input_banff_26, bestfit_banff_26, errp_banff_26, errm_banff_26, errpar_banff_26;
  double input_banff_27, bestfit_banff_27, errp_banff_27, errm_banff_27, errpar_banff_27;
  double input_banff_28, bestfit_banff_28, errp_banff_28, errm_banff_28, errpar_banff_28;
  double input_banff_29, bestfit_banff_29, errp_banff_29, errm_banff_29, errpar_banff_29;
  double input_banff_30, bestfit_banff_30, errp_banff_30, errm_banff_30, errpar_banff_30;
  double input_banff_maqe, bestfit_banff_maqe, errp_banff_maqe, errm_banff_maqe, errpar_banff_maqe;
  double input_banff_mares, bestfit_banff_mares, errp_banff_mares, errm_banff_mares, errpar_banff_mares;
  double input_banff_ccqeE0, bestfit_banff_ccqeE0, errp_banff_ccqeE0, errm_banff_ccqeE0, errpar_banff_ccqeE0;
  double input_banff_ccqeE1, bestfit_banff_ccqeE1, errp_banff_ccqeE1, errm_banff_ccqeE1, errpar_banff_ccqeE1;
  double input_banff_ccqeE2, bestfit_banff_ccqeE2, errp_banff_ccqeE2, errm_banff_ccqeE2, errpar_banff_ccqeE2;
  double input_banff_cc1piE0, bestfit_banff_cc1piE0, errp_banff_cc1piE0, errm_banff_cc1piE0, errpar_banff_cc1piE0;
  double input_banff_cc1piE1, bestfit_banff_cc1piE1, errp_banff_cc1piE1, errm_banff_cc1piE1, errpar_banff_cc1piE1;
  double input_banff_nc1pi0, bestfit_banff_nc1pi0, errp_banff_nc1pi0, errm_banff_nc1pi0, errpar_banff_nc1pi0;
  double input_xsec_pf, bestfit_xsec_pf, errp_xsec_pf, errm_xsec_pf, errpar_xsec_pf;
  double input_xsec_sf, bestfit_xsec_sf, errp_xsec_sf, errm_xsec_sf, errpar_xsec_sf;
  double input_xsec_wshape, bestfit_xsec_wshape, errp_xsec_wshape, errm_xsec_wshape, errpar_xsec_wshape;
  double input_xsec_ccothshape, bestfit_xsec_ccothshape, errp_xsec_ccothshape, errm_xsec_ccothshape, errpar_xsec_ccothshape;
  double input_xsec_cc1pishape, bestfit_xsec_cc1pishape, errp_xsec_cc1pishape, errm_xsec_cc1pishape, errpar_xsec_cc1pishape;
  double input_xsec_pilessdelta, bestfit_xsec_pilessdelta, errp_xsec_pilessdelta, errm_xsec_pilessdelta, errpar_xsec_pilessdelta;
  double input_xsec_normcccoh, bestfit_xsec_normcccoh, errp_xsec_normcccoh, errm_xsec_normcccoh, errpar_xsec_normcccoh;
  double input_xsec_normnc1piC, bestfit_xsec_normnc1piC, errp_xsec_normnc1piC, errm_xsec_normnc1piC, errpar_xsec_normnc1piC;
  double input_xsec_normnccoh, bestfit_xsec_normnccoh, errp_xsec_normnccoh, errm_xsec_normnccoh, errpar_xsec_normnccoh;
  double input_xsec_normncoth, bestfit_xsec_normncoth, errp_xsec_normncoth, errm_xsec_normncoth, errpar_xsec_normncoth;
  double input_xsec_normnuetonumu, bestfit_xsec_normnuetonumu, errp_xsec_normnuetonumu, errm_xsec_normnuetonumu, errpar_xsec_normnuetonumu;
  double input_xsec_normantinutonu, bestfit_xsec_normantinutonu, errp_xsec_normantinutonu, errm_xsec_normantinutonu, errpar_xsec_normantinutonu;
  double input_skdet_energyscale, bestfit_skdet_energyscale, errp_skdet_energyscale, errm_skdet_energyscale, errpar_skdet_energyscale;
  double input_skdet_ccqerc0, bestfit_skdet_ccqerc0, errp_skdet_ccqerc0, errm_skdet_ccqerc0, errpar_skdet_ccqerc0;
  double input_skdet_ccqerc1, bestfit_skdet_ccqerc1, errp_skdet_ccqerc1, errm_skdet_ccqerc1, errpar_skdet_ccqerc1;
  double input_skdet_ccqerc2, bestfit_skdet_ccqerc2, errp_skdet_ccqerc2, errm_skdet_ccqerc2, errpar_skdet_ccqerc2;
  double input_skdet_ccqeoth, bestfit_skdet_ccqeoth, errp_skdet_ccqeoth, errm_skdet_ccqeoth, errpar_skdet_ccqeoth;
  double input_skdet_ccnqe, bestfit_skdet_ccnqe, errp_skdet_ccnqe, errm_skdet_ccnqe, errpar_skdet_ccnqe;
  double input_skdet_nc, bestfit_skdet_nc, errp_skdet_nc, errm_skdet_nc, errpar_skdet_nc;
  double input_skdet_nuecc, bestfit_skdet_nuecc, errp_skdet_nuecc, errm_skdet_nuecc, errpar_skdet_nuecc;
  double input_fsi_000, bestfit_fsi_000, errp_fsi_000, errm_fsi_000, errpar_fsi_000;
  double input_fsi_001, bestfit_fsi_001, errp_fsi_001, errm_fsi_001, errpar_fsi_001;
  double input_fsi_002, bestfit_fsi_002, errp_fsi_002, errm_fsi_002, errpar_fsi_002;
  double input_fsi_003, bestfit_fsi_003, errp_fsi_003, errm_fsi_003, errpar_fsi_003;
  double input_fsi_004, bestfit_fsi_004, errp_fsi_004, errm_fsi_004, errpar_fsi_004;
  double input_fsi_005, bestfit_fsi_005, errp_fsi_005, errm_fsi_005, errpar_fsi_005;
  double input_fsi_006, bestfit_fsi_006, errp_fsi_006, errm_fsi_006, errpar_fsi_006;
  double input_fsi_007, bestfit_fsi_007, errp_fsi_007, errm_fsi_007, errpar_fsi_007;
  double input_fsi_008, bestfit_fsi_008, errp_fsi_008, errm_fsi_008, errpar_fsi_008;
  double input_fsi_009, bestfit_fsi_009, errp_fsi_009, errm_fsi_009, errpar_fsi_009;
  double input_fsi_010, bestfit_fsi_010, errp_fsi_010, errm_fsi_010, errpar_fsi_010;
  double input_fsi_011, bestfit_fsi_011, errp_fsi_011, errm_fsi_011, errpar_fsi_011;
  double input_fsi_012, bestfit_fsi_012, errp_fsi_012, errm_fsi_012, errpar_fsi_012;
  double input_fsi_013, bestfit_fsi_013, errp_fsi_013, errm_fsi_013, errpar_fsi_013;
  double input_fsi_014, bestfit_fsi_014, errp_fsi_014, errm_fsi_014, errpar_fsi_014;

  // ***** setup the branches *******
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
  if(fit_syst_mode == 1) {
    nfitnt->Branch("input_banff_maqe",&input_banff_maqe,"input_banff_maqe/D");
    nfitnt->Branch("input_banff_ccqeE0",&input_banff_ccqeE0,"input_banff_ccqeE0/D");
    nfitnt->Branch("input_banff_cc1piE0",&input_banff_cc1piE0,"input_banff_cc1piE0/D");
    nfitnt->Branch("input_xsec_sf",&input_xsec_sf,"input_xsec_sf/D");
    nfitnt->Branch("input_xsec_pilessdelta",&input_xsec_pilessdelta,"input_xsec_pilessdelta/D");
    nfitnt->Branch("input_skdet_ccnqe",&input_skdet_ccnqe,"input_skdet_ccnqe/D");
    nfitnt->Branch("input_skdet_nc",&input_skdet_nc,"input_skdet_nc/D");
  }
  else if(fit_syst_mode == 2) {
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
    nfitnt->Branch("input_banff_16",&input_banff_16,"input_banff_16/D");
    nfitnt->Branch("input_banff_17",&input_banff_17,"input_banff_17/D");
    nfitnt->Branch("input_banff_18",&input_banff_18,"input_banff_18/D");
    nfitnt->Branch("input_banff_19",&input_banff_19,"input_banff_19/D");
    nfitnt->Branch("input_banff_20",&input_banff_20,"input_banff_20/D");
    nfitnt->Branch("input_banff_21",&input_banff_21,"input_banff_21/D");
    nfitnt->Branch("input_banff_22",&input_banff_22,"input_banff_22/D");
    nfitnt->Branch("input_banff_23",&input_banff_23,"input_banff_23/D");
    nfitnt->Branch("input_banff_24",&input_banff_24,"input_banff_24/D");
    nfitnt->Branch("input_banff_25",&input_banff_25,"input_banff_25/D");
    nfitnt->Branch("input_banff_26",&input_banff_26,"input_banff_26/D");
    nfitnt->Branch("input_banff_27",&input_banff_27,"input_banff_27/D");
    nfitnt->Branch("input_banff_28",&input_banff_28,"input_banff_28/D");
    nfitnt->Branch("input_banff_29",&input_banff_29,"input_banff_29/D");
    nfitnt->Branch("input_banff_30",&input_banff_30,"input_banff_30/D");
    nfitnt->Branch("input_banff_maqe",&input_banff_maqe,"input_banff_maqe/D");
    nfitnt->Branch("input_banff_mares",&input_banff_mares,"input_banff_mares/D");
    nfitnt->Branch("input_banff_ccqeE0",&input_banff_ccqeE0,"input_banff_ccqeE0/D");
    nfitnt->Branch("input_banff_ccqeE1",&input_banff_ccqeE1,"input_banff_ccqeE1/D");
    nfitnt->Branch("input_banff_ccqeE2",&input_banff_ccqeE2,"input_banff_ccqeE2/D");
    nfitnt->Branch("input_banff_cc1piE0",&input_banff_cc1piE0,"input_banff_cc1piE0/D");
    nfitnt->Branch("input_banff_cc1piE1",&input_banff_cc1piE1,"input_banff_cc1piE1/D");
    nfitnt->Branch("input_banff_nc1pi0",&input_banff_nc1pi0,"input_banff_nc1pi0/D");
    nfitnt->Branch("input_xsec_pf",&input_xsec_pf,"input_xsec_pf/D");
    nfitnt->Branch("input_xsec_sf",&input_xsec_sf,"input_xsec_sf/D");
    nfitnt->Branch("input_xsec_wshape",&input_xsec_wshape,"input_xsec_wshape/D");
    nfitnt->Branch("input_xsec_ccothshape",&input_xsec_ccothshape,"input_xsec_ccothshape/D");
    nfitnt->Branch("input_xsec_cc1pishape",&input_xsec_cc1pishape,"input_xsec_cc1pishape/D");
    nfitnt->Branch("input_xsec_pilessdelta",&input_xsec_pilessdelta,"input_xsec_pilessdelta/D");
    nfitnt->Branch("input_xsec_normcccoh",&input_xsec_normcccoh,"input_xsec_normcccoh/D");
    nfitnt->Branch("input_xsec_normnc1piC",&input_xsec_normnc1piC,"input_xsec_normnc1piC/D");
    nfitnt->Branch("input_xsec_normnccoh",&input_xsec_normnccoh,"input_xsec_normnccoh/D");
    nfitnt->Branch("input_xsec_normncoth",&input_xsec_normncoth,"input_xsec_normncoth/D");
    nfitnt->Branch("input_xsec_normnuetonumu",&input_xsec_normnuetonumu,"input_xsec_normnuetonumu/D");
    nfitnt->Branch("input_xsec_normantinutonu",&input_xsec_normantinutonu,"input_xsec_normantinutonu/D");
    nfitnt->Branch("input_skdet_energyscale",&input_skdet_energyscale,"input_skdet_energyscale/D");
    nfitnt->Branch("input_skdet_ccqerc0",&input_skdet_ccqerc0,"input_skdet_ccqerc0/D");
    nfitnt->Branch("input_skdet_ccqerc1",&input_skdet_ccqerc1,"input_skdet_ccqerc1/D");
    nfitnt->Branch("input_skdet_ccqerc2",&input_skdet_ccqerc2,"input_skdet_ccqerc2/D");
    nfitnt->Branch("input_skdet_ccqeoth",&input_skdet_ccqeoth,"input_skdet_ccqeoth/D");
    nfitnt->Branch("input_skdet_ccnqe",&input_skdet_ccnqe,"input_skdet_ccnqe/D");
    nfitnt->Branch("input_skdet_nc",&input_skdet_nc,"input_skdet_nc/D");
    nfitnt->Branch("input_skdet_nuecc",&input_skdet_nuecc,"input_skdet_nuecc/D");
    nfitnt->Branch("input_fsi_000",&input_fsi_000,"input_fsi_000/D");
    nfitnt->Branch("input_fsi_001",&input_fsi_001,"input_fsi_001/D");
    nfitnt->Branch("input_fsi_002",&input_fsi_002,"input_fsi_002/D");
    nfitnt->Branch("input_fsi_003",&input_fsi_003,"input_fsi_003/D");
    nfitnt->Branch("input_fsi_004",&input_fsi_004,"input_fsi_004/D");
    nfitnt->Branch("input_fsi_005",&input_fsi_005,"input_fsi_005/D");
    if(conv::kNumOfFSIParams == 15) {
      nfitnt->Branch("input_fsi_006",&input_fsi_006,"input_fsi_006/D");
      nfitnt->Branch("input_fsi_007",&input_fsi_007,"input_fsi_007/D");
      nfitnt->Branch("input_fsi_008",&input_fsi_008,"input_fsi_008/D");
      nfitnt->Branch("input_fsi_009",&input_fsi_009,"input_fsi_009/D");
      nfitnt->Branch("input_fsi_010",&input_fsi_010,"input_fsi_010/D");
      nfitnt->Branch("input_fsi_011",&input_fsi_011,"input_fsi_011/D");
      nfitnt->Branch("input_fsi_012",&input_fsi_012,"input_fsi_012/D");
      nfitnt->Branch("input_fsi_013",&input_fsi_013,"input_fsi_013/D");
      nfitnt->Branch("input_fsi_014",&input_fsi_014,"input_fsi_014/D");
    }
  }
  //bestfits
  nfitnt->Branch("bestfit_mixing12",&bestfit_mixing12,"bestfit_mixing12/D");
  nfitnt->Branch("bestfit_mixing13",&bestfit_mixing13,"bestfit_mixing13/D");
  nfitnt->Branch("bestfit_mixing23",&bestfit_mixing23,"bestfit_mixing23/D");
  nfitnt->Branch("bestfit_deltacp",&bestfit_deltacp,"bestfit_deltacp/D");
  nfitnt->Branch("bestfit_dm12sq",&bestfit_dm12sq,"bestfit_dm12sq/D");
  nfitnt->Branch("bestfit_dm23sq",&bestfit_dm23sq,"bestfit_dm23sq/D");
  if(fit_syst_mode == 1) {
    nfitnt->Branch("bestfit_banff_maqe",&bestfit_banff_maqe,"bestfit_banff_maqe/D");
    nfitnt->Branch("bestfit_banff_ccqeE0",&bestfit_banff_ccqeE0,"bestfit_banff_ccqeE0/D");
    nfitnt->Branch("bestfit_banff_cc1piE0",&bestfit_banff_cc1piE0,"bestfit_banff_cc1piE0/D");
    nfitnt->Branch("bestfit_xsec_sf",&bestfit_xsec_sf,"bestfit_xsec_sf/D");
    nfitnt->Branch("bestfit_xsec_pilessdelta",&bestfit_xsec_pilessdelta,"bestfit_xsec_pilessdelta/D");
    nfitnt->Branch("bestfit_skdet_ccnqe",&bestfit_skdet_ccnqe,"bestfit_skdet_ccnqe/D");
    nfitnt->Branch("bestfit_skdet_nc",&bestfit_skdet_nc,"bestfit_skdet_nc/D");
  }
  else if(fit_syst_mode == 2) {
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
    nfitnt->Branch("bestfit_banff_16",&bestfit_banff_16,"bestfit_banff_16/D");
    nfitnt->Branch("bestfit_banff_17",&bestfit_banff_17,"bestfit_banff_17/D");
    nfitnt->Branch("bestfit_banff_18",&bestfit_banff_18,"bestfit_banff_18/D");
    nfitnt->Branch("bestfit_banff_19",&bestfit_banff_19,"bestfit_banff_19/D");
    nfitnt->Branch("bestfit_banff_20",&bestfit_banff_20,"bestfit_banff_20/D");
    nfitnt->Branch("bestfit_banff_21",&bestfit_banff_21,"bestfit_banff_21/D");
    nfitnt->Branch("bestfit_banff_22",&bestfit_banff_22,"bestfit_banff_22/D");
    nfitnt->Branch("bestfit_banff_23",&bestfit_banff_23,"bestfit_banff_23/D");
    nfitnt->Branch("bestfit_banff_24",&bestfit_banff_24,"bestfit_banff_24/D");
    nfitnt->Branch("bestfit_banff_25",&bestfit_banff_25,"bestfit_banff_25/D");
    nfitnt->Branch("bestfit_banff_26",&bestfit_banff_26,"bestfit_banff_26/D");
    nfitnt->Branch("bestfit_banff_27",&bestfit_banff_27,"bestfit_banff_27/D");
    nfitnt->Branch("bestfit_banff_28",&bestfit_banff_28,"bestfit_banff_28/D");
    nfitnt->Branch("bestfit_banff_29",&bestfit_banff_29,"bestfit_banff_29/D");
    nfitnt->Branch("bestfit_banff_30",&bestfit_banff_30,"bestfit_banff_30/D");
    nfitnt->Branch("bestfit_banff_maqe",&bestfit_banff_maqe,"bestfit_banff_maqe/D");
    nfitnt->Branch("bestfit_banff_mares",&bestfit_banff_mares,"bestfit_banff_mares/D");
    nfitnt->Branch("bestfit_banff_ccqeE0",&bestfit_banff_ccqeE0,"bestfit_banff_ccqeE0/D");
    nfitnt->Branch("bestfit_banff_ccqeE1",&bestfit_banff_ccqeE1,"bestfit_banff_ccqeE1/D");
    nfitnt->Branch("bestfit_banff_ccqeE2",&bestfit_banff_ccqeE2,"bestfit_banff_ccqeE2/D");
    nfitnt->Branch("bestfit_banff_cc1piE0",&bestfit_banff_cc1piE0,"bestfit_banff_cc1piE0/D");
    nfitnt->Branch("bestfit_banff_cc1piE1",&bestfit_banff_cc1piE1,"bestfit_banff_cc1piE1/D");
    nfitnt->Branch("bestfit_banff_nc1pi0",&bestfit_banff_nc1pi0,"bestfit_banff_nc1pi0/D");
    nfitnt->Branch("bestfit_xsec_pf",&bestfit_xsec_pf,"bestfit_xsec_pf/D");
    nfitnt->Branch("bestfit_xsec_sf",&bestfit_xsec_sf,"bestfit_xsec_sf/D");
    nfitnt->Branch("bestfit_xsec_wshape",&bestfit_xsec_wshape,"bestfit_xsec_wshape/D");
    nfitnt->Branch("bestfit_xsec_ccothshape",&bestfit_xsec_ccothshape,"bestfit_xsec_ccothshape/D");
    nfitnt->Branch("bestfit_xsec_cc1pishape",&bestfit_xsec_cc1pishape,"bestfit_xsec_cc1pishape/D");
    nfitnt->Branch("bestfit_xsec_pilessdelta",&bestfit_xsec_pilessdelta,"bestfit_xsec_pilessdelta/D");
    nfitnt->Branch("bestfit_xsec_normcccoh",&bestfit_xsec_normcccoh,"bestfit_xsec_normcccoh/D");
    nfitnt->Branch("bestfit_xsec_normnc1piC",&bestfit_xsec_normnc1piC,"bestfit_xsec_normnc1piC/D");
    nfitnt->Branch("bestfit_xsec_normnccoh",&bestfit_xsec_normnccoh,"bestfit_xsec_normnccoh/D");
    nfitnt->Branch("bestfit_xsec_normncoth",&bestfit_xsec_normncoth,"bestfit_xsec_normncoth/D");
    nfitnt->Branch("bestfit_xsec_normnuetonumu",&bestfit_xsec_normnuetonumu,"bestfit_xsec_normnuetonumu/D");
    nfitnt->Branch("bestfit_xsec_normantinutonu",&bestfit_xsec_normantinutonu,"bestfit_xsec_normantinutonu/D");
    nfitnt->Branch("bestfit_skdet_energyscale",&bestfit_skdet_energyscale,"bestfit_skdet_energyscale/D");
    nfitnt->Branch("bestfit_skdet_ccqerc0",&bestfit_skdet_ccqerc0,"bestfit_skdet_ccqerc0/D");
    nfitnt->Branch("bestfit_skdet_ccqerc1",&bestfit_skdet_ccqerc1,"bestfit_skdet_ccqerc1/D");
    nfitnt->Branch("bestfit_skdet_ccqerc2",&bestfit_skdet_ccqerc2,"bestfit_skdet_ccqerc2/D");
    nfitnt->Branch("bestfit_skdet_ccqeoth",&bestfit_skdet_ccqeoth,"bestfit_skdet_ccqeoth/D");
    nfitnt->Branch("bestfit_skdet_ccnqe",&bestfit_skdet_ccnqe,"bestfit_skdet_ccnqe/D");
    nfitnt->Branch("bestfit_skdet_nc",&bestfit_skdet_nc,"bestfit_skdet_nc/D");
    nfitnt->Branch("bestfit_skdet_nuecc",&bestfit_skdet_nuecc,"bestfit_skdet_nuecc/D");
    nfitnt->Branch("bestfit_fsi_000",&bestfit_fsi_000,"bestfit_fsi_000/D");
    nfitnt->Branch("bestfit_fsi_001",&bestfit_fsi_001,"bestfit_fsi_001/D");
    nfitnt->Branch("bestfit_fsi_002",&bestfit_fsi_002,"bestfit_fsi_002/D");
    nfitnt->Branch("bestfit_fsi_003",&bestfit_fsi_003,"bestfit_fsi_003/D");
    nfitnt->Branch("bestfit_fsi_004",&bestfit_fsi_004,"bestfit_fsi_004/D");
    nfitnt->Branch("bestfit_fsi_005",&bestfit_fsi_005,"bestfit_fsi_005/D");
    if(conv::kNumOfFSIParams == 15) {
      nfitnt->Branch("bestfit_fsi_006",&bestfit_fsi_006,"bestfit_fsi_006/D");
      nfitnt->Branch("bestfit_fsi_007",&bestfit_fsi_007,"bestfit_fsi_007/D");
      nfitnt->Branch("bestfit_fsi_008",&bestfit_fsi_008,"bestfit_fsi_008/D");
      nfitnt->Branch("bestfit_fsi_009",&bestfit_fsi_009,"bestfit_fsi_009/D");
      nfitnt->Branch("bestfit_fsi_010",&bestfit_fsi_010,"bestfit_fsi_010/D");
      nfitnt->Branch("bestfit_fsi_011",&bestfit_fsi_011,"bestfit_fsi_011/D");
      nfitnt->Branch("bestfit_fsi_012",&bestfit_fsi_012,"bestfit_fsi_012/D");
      nfitnt->Branch("bestfit_fsi_013",&bestfit_fsi_013,"bestfit_fsi_013/D");
      nfitnt->Branch("bestfit_fsi_014",&bestfit_fsi_014,"bestfit_fsi_014/D");
    }
  }
  //plus errors
  nfitnt->Branch("errp_mixing12",&errp_mixing12,"errp_mixing12/D");
  nfitnt->Branch("errp_mixing13",&errp_mixing13,"errp_mixing13/D");
  nfitnt->Branch("errp_mixing23",&errp_mixing23,"errp_mixing23/D");
  nfitnt->Branch("errp_deltacp",&errp_deltacp,"errp_deltacp/D");
  nfitnt->Branch("errp_dm12sq",&errp_dm12sq,"errp_dm12sq/D");
  nfitnt->Branch("errp_dm23sq",&errp_dm23sq,"errp_dm23sq/D");
  if(fit_syst_mode == 1) {
    nfitnt->Branch("errp_banff_maqe",&errp_banff_maqe,"errp_banff_maqe/D");
    nfitnt->Branch("errp_banff_ccqeE0",&errp_banff_ccqeE0,"errp_banff_ccqeE0/D");
    nfitnt->Branch("errp_banff_cc1piE0",&errp_banff_cc1piE0,"errp_banff_cc1piE0/D");
    nfitnt->Branch("errp_xsec_sf",&errp_xsec_sf,"errp_xsec_sf/D");
    nfitnt->Branch("errp_xsec_pilessdelta",&errp_xsec_pilessdelta,"errp_xsec_pilessdelta/D");
    nfitnt->Branch("errp_skdet_ccnqe",&errp_skdet_ccnqe,"errp_skdet_ccnqe/D");
    nfitnt->Branch("errp_skdet_nc",&errp_skdet_nc,"errp_skdet_nc/D");
  }
  else if(fit_syst_mode == 2) {
    nfitnt->Branch("errp_banff_00",&errp_banff_00,"errp_banff_00/D");
    nfitnt->Branch("errp_banff_01",&errp_banff_01,"errp_banff_01/D");
    nfitnt->Branch("errp_banff_02",&errp_banff_02,"errp_banff_02/D");
    nfitnt->Branch("errp_banff_03",&errp_banff_03,"errp_banff_03/D");
    nfitnt->Branch("errp_banff_04",&errp_banff_04,"errp_banff_04/D");
    nfitnt->Branch("errp_banff_05",&errp_banff_05,"errp_banff_05/D");
    nfitnt->Branch("errp_banff_06",&errp_banff_06,"errp_banff_06/D");
    nfitnt->Branch("errp_banff_07",&errp_banff_07,"errp_banff_07/D");
    nfitnt->Branch("errp_banff_08",&errp_banff_08,"errp_banff_08/D");
    nfitnt->Branch("errp_banff_09",&errp_banff_09,"errp_banff_09/D");
    nfitnt->Branch("errp_banff_10",&errp_banff_10,"errp_banff_10/D");
    nfitnt->Branch("errp_banff_11",&errp_banff_11,"errp_banff_11/D");
    nfitnt->Branch("errp_banff_12",&errp_banff_12,"errp_banff_12/D");
    nfitnt->Branch("errp_banff_13",&errp_banff_13,"errp_banff_13/D");
    nfitnt->Branch("errp_banff_14",&errp_banff_14,"errp_banff_14/D");
    nfitnt->Branch("errp_banff_15",&errp_banff_15,"errp_banff_15/D");
    nfitnt->Branch("errp_banff_16",&errp_banff_16,"errp_banff_16/D");
    nfitnt->Branch("errp_banff_17",&errp_banff_17,"errp_banff_17/D");
    nfitnt->Branch("errp_banff_18",&errp_banff_18,"errp_banff_18/D");
    nfitnt->Branch("errp_banff_19",&errp_banff_19,"errp_banff_19/D");
    nfitnt->Branch("errp_banff_20",&errp_banff_20,"errp_banff_20/D");
    nfitnt->Branch("errp_banff_21",&errp_banff_21,"errp_banff_21/D");
    nfitnt->Branch("errp_banff_22",&errp_banff_22,"errp_banff_22/D");
    nfitnt->Branch("errp_banff_23",&errp_banff_23,"errp_banff_23/D");
    nfitnt->Branch("errp_banff_24",&errp_banff_24,"errp_banff_24/D");
    nfitnt->Branch("errp_banff_25",&errp_banff_25,"errp_banff_25/D");
    nfitnt->Branch("errp_banff_26",&errp_banff_26,"errp_banff_26/D");
    nfitnt->Branch("errp_banff_27",&errp_banff_27,"errp_banff_27/D");
    nfitnt->Branch("errp_banff_28",&errp_banff_28,"errp_banff_28/D");
    nfitnt->Branch("errp_banff_29",&errp_banff_29,"errp_banff_29/D");
    nfitnt->Branch("errp_banff_30",&errp_banff_30,"errp_banff_30/D");
    nfitnt->Branch("errp_banff_maqe",&errp_banff_maqe,"errp_banff_maqe/D");
    nfitnt->Branch("errp_banff_mares",&errp_banff_mares,"errp_banff_mares/D");
    nfitnt->Branch("errp_banff_ccqeE0",&errp_banff_ccqeE0,"errp_banff_ccqeE0/D");
    nfitnt->Branch("errp_banff_ccqeE1",&errp_banff_ccqeE1,"errp_banff_ccqeE1/D");
    nfitnt->Branch("errp_banff_ccqeE2",&errp_banff_ccqeE2,"errp_banff_ccqeE2/D");
    nfitnt->Branch("errp_banff_cc1piE0",&errp_banff_cc1piE0,"errp_banff_cc1piE0/D");
    nfitnt->Branch("errp_banff_cc1piE1",&errp_banff_cc1piE1,"errp_banff_cc1piE1/D");
    nfitnt->Branch("errp_banff_nc1pi0",&errp_banff_nc1pi0,"errp_banff_nc1pi0/D");
    nfitnt->Branch("errp_xsec_pf",&errp_xsec_pf,"errp_xsec_pf/D");
    nfitnt->Branch("errp_xsec_sf",&errp_xsec_sf,"errp_xsec_sf/D");
    nfitnt->Branch("errp_xsec_wshape",&errp_xsec_wshape,"errp_xsec_wshape/D");
    nfitnt->Branch("errp_xsec_ccothshape",&errp_xsec_ccothshape,"errp_xsec_ccothshape/D");
    nfitnt->Branch("errp_xsec_cc1pishape",&errp_xsec_cc1pishape,"errp_xsec_cc1pishape/D");
    nfitnt->Branch("errp_xsec_pilessdelta",&errp_xsec_pilessdelta,"errp_xsec_pilessdelta/D");
    nfitnt->Branch("errp_xsec_normcccoh",&errp_xsec_normcccoh,"errp_xsec_normcccoh/D");
    nfitnt->Branch("errp_xsec_normnc1piC",&errp_xsec_normnc1piC,"errp_xsec_normnc1piC/D");
    nfitnt->Branch("errp_xsec_normnccoh",&errp_xsec_normnccoh,"errp_xsec_normnccoh/D");
    nfitnt->Branch("errp_xsec_normncoth",&errp_xsec_normncoth,"errp_xsec_normncoth/D");
    nfitnt->Branch("errp_xsec_normnuetonumu",&errp_xsec_normnuetonumu,"errp_xsec_normnuetonumu/D");
    nfitnt->Branch("errp_xsec_normantinutonu",&errp_xsec_normantinutonu,"errp_xsec_normantinutonu/D");
    nfitnt->Branch("errp_skdet_energyscale",&errp_skdet_energyscale,"errp_skdet_energyscale/D");
    nfitnt->Branch("errp_skdet_ccqerc0",&errp_skdet_ccqerc0,"errp_skdet_ccqerc0/D");
    nfitnt->Branch("errp_skdet_ccqerc1",&errp_skdet_ccqerc1,"errp_skdet_ccqerc1/D");
    nfitnt->Branch("errp_skdet_ccqerc2",&errp_skdet_ccqerc2,"errp_skdet_ccqerc2/D");
    nfitnt->Branch("errp_skdet_ccqeoth",&errp_skdet_ccqeoth,"errp_skdet_ccqeoth/D");
    nfitnt->Branch("errp_skdet_ccnqe",&errp_skdet_ccnqe,"errp_skdet_ccnqe/D");
    nfitnt->Branch("errp_skdet_nc",&errp_skdet_nc,"errp_skdet_nc/D");
    nfitnt->Branch("errp_skdet_nuecc",&errp_skdet_nuecc,"errp_skdet_nuecc/D");
    nfitnt->Branch("errp_fsi_000",&errp_fsi_000,"errp_fsi_000/D");
    nfitnt->Branch("errp_fsi_001",&errp_fsi_001,"errp_fsi_001/D");
    nfitnt->Branch("errp_fsi_002",&errp_fsi_002,"errp_fsi_002/D");
    nfitnt->Branch("errp_fsi_003",&errp_fsi_003,"errp_fsi_003/D");
    nfitnt->Branch("errp_fsi_004",&errp_fsi_004,"errp_fsi_004/D");
    nfitnt->Branch("errp_fsi_005",&errp_fsi_005,"errp_fsi_005/D");
    if(conv::kNumOfFSIParams == 15) {
      nfitnt->Branch("errp_fsi_006",&errp_fsi_006,"errp_fsi_006/D");
      nfitnt->Branch("errp_fsi_007",&errp_fsi_007,"errp_fsi_007/D");
      nfitnt->Branch("errp_fsi_008",&errp_fsi_008,"errp_fsi_008/D");
      nfitnt->Branch("errp_fsi_009",&errp_fsi_009,"errp_fsi_009/D");
      nfitnt->Branch("errp_fsi_010",&errp_fsi_010,"errp_fsi_010/D");
      nfitnt->Branch("errp_fsi_011",&errp_fsi_011,"errp_fsi_011/D");
      nfitnt->Branch("errp_fsi_012",&errp_fsi_012,"errp_fsi_012/D");
      nfitnt->Branch("errp_fsi_013",&errp_fsi_013,"errp_fsi_013/D");
      nfitnt->Branch("errp_fsi_014",&errp_fsi_014,"errp_fsi_014/D");
    }
  }
  //minus errors
  nfitnt->Branch("errm_mixing12",&errm_mixing12,"errm_mixing12/D");
  nfitnt->Branch("errm_mixing13",&errm_mixing13,"errm_mixing13/D");
  nfitnt->Branch("errm_mixing23",&errm_mixing23,"errm_mixing23/D");
  nfitnt->Branch("errm_deltacp",&errm_deltacp,"errm_deltacp/D");
  nfitnt->Branch("errm_dm12sq",&errm_dm12sq,"errm_dm12sq/D");
  nfitnt->Branch("errm_dm23sq",&errm_dm23sq,"errm_dm23sq/D");
  if(fit_syst_mode == 1) {
    nfitnt->Branch("errm_banff_maqe",&errm_banff_maqe,"errm_banff_maqe/D");
    nfitnt->Branch("errm_banff_ccqeE0",&errm_banff_ccqeE0,"errm_banff_ccqeE0/D");
    nfitnt->Branch("errm_banff_cc1piE0",&errm_banff_cc1piE0,"errm_banff_cc1piE0/D");
    nfitnt->Branch("errm_xsec_sf",&errm_xsec_sf,"errm_xsec_sf/D");
    nfitnt->Branch("errm_xsec_pilessdelta",&errm_xsec_pilessdelta,"errm_xsec_pilessdelta/D");
    nfitnt->Branch("errm_skdet_ccnqe",&errm_skdet_ccnqe,"errm_skdet_ccnqe/D");
    nfitnt->Branch("errm_skdet_nc",&errm_skdet_nc,"errm_skdet_nc/D");
  }
  else if(fit_syst_mode == 2) {
    nfitnt->Branch("errm_banff_00",&errm_banff_00,"errm_banff_00/D");
    nfitnt->Branch("errm_banff_01",&errm_banff_01,"errm_banff_01/D");
    nfitnt->Branch("errm_banff_02",&errm_banff_02,"errm_banff_02/D");
    nfitnt->Branch("errm_banff_03",&errm_banff_03,"errm_banff_03/D");
    nfitnt->Branch("errm_banff_04",&errm_banff_04,"errm_banff_04/D");
    nfitnt->Branch("errm_banff_05",&errm_banff_05,"errm_banff_05/D");
    nfitnt->Branch("errm_banff_06",&errm_banff_06,"errm_banff_06/D");
    nfitnt->Branch("errm_banff_07",&errm_banff_07,"errm_banff_07/D");
    nfitnt->Branch("errm_banff_08",&errm_banff_08,"errm_banff_08/D");
    nfitnt->Branch("errm_banff_09",&errm_banff_09,"errm_banff_09/D");
    nfitnt->Branch("errm_banff_10",&errm_banff_10,"errm_banff_10/D");
    nfitnt->Branch("errm_banff_11",&errm_banff_11,"errm_banff_11/D");
    nfitnt->Branch("errm_banff_12",&errm_banff_12,"errm_banff_12/D");
    nfitnt->Branch("errm_banff_13",&errm_banff_13,"errm_banff_13/D");
    nfitnt->Branch("errm_banff_14",&errm_banff_14,"errm_banff_14/D");
    nfitnt->Branch("errm_banff_15",&errm_banff_15,"errm_banff_15/D");
    nfitnt->Branch("errm_banff_16",&errm_banff_16,"errm_banff_16/D");
    nfitnt->Branch("errm_banff_17",&errm_banff_17,"errm_banff_17/D");
    nfitnt->Branch("errm_banff_18",&errm_banff_18,"errm_banff_18/D");
    nfitnt->Branch("errm_banff_19",&errm_banff_19,"errm_banff_19/D");
    nfitnt->Branch("errm_banff_20",&errm_banff_20,"errm_banff_20/D");
    nfitnt->Branch("errm_banff_21",&errm_banff_21,"errm_banff_21/D");
    nfitnt->Branch("errm_banff_22",&errm_banff_22,"errm_banff_22/D");
    nfitnt->Branch("errm_banff_23",&errm_banff_23,"errm_banff_23/D");
    nfitnt->Branch("errm_banff_24",&errm_banff_24,"errm_banff_24/D");
    nfitnt->Branch("errm_banff_25",&errm_banff_25,"errm_banff_25/D");
    nfitnt->Branch("errm_banff_26",&errm_banff_26,"errm_banff_26/D");
    nfitnt->Branch("errm_banff_27",&errm_banff_27,"errm_banff_27/D");
    nfitnt->Branch("errm_banff_28",&errm_banff_28,"errm_banff_28/D");
    nfitnt->Branch("errm_banff_29",&errm_banff_29,"errm_banff_29/D");
    nfitnt->Branch("errm_banff_30",&errm_banff_30,"errm_banff_30/D");
    nfitnt->Branch("errm_banff_maqe",&errm_banff_maqe,"errm_banff_maqe/D");
    nfitnt->Branch("errm_banff_mares",&errm_banff_mares,"errm_banff_mares/D");
    nfitnt->Branch("errm_banff_ccqeE0",&errm_banff_ccqeE0,"errm_banff_ccqeE0/D");
    nfitnt->Branch("errm_banff_ccqeE1",&errm_banff_ccqeE1,"errm_banff_ccqeE1/D");
    nfitnt->Branch("errm_banff_ccqeE2",&errm_banff_ccqeE2,"errm_banff_ccqeE2/D");
    nfitnt->Branch("errm_banff_cc1piE0",&errm_banff_cc1piE0,"errm_banff_cc1piE0/D");
    nfitnt->Branch("errm_banff_cc1piE1",&errm_banff_cc1piE1,"errm_banff_cc1piE1/D");
    nfitnt->Branch("errm_banff_nc1pi0",&errm_banff_nc1pi0,"errm_banff_nc1pi0/D");
    nfitnt->Branch("errm_xsec_pf",&errm_xsec_pf,"errm_xsec_pf/D");
    nfitnt->Branch("errm_xsec_sf",&errm_xsec_sf,"errm_xsec_sf/D");
    nfitnt->Branch("errm_xsec_wshape",&errm_xsec_wshape,"errm_xsec_wshape/D");
    nfitnt->Branch("errm_xsec_ccothshape",&errm_xsec_ccothshape,"errm_xsec_ccothshape/D");
    nfitnt->Branch("errm_xsec_cc1pishape",&errm_xsec_cc1pishape,"errm_xsec_cc1pishape/D");
    nfitnt->Branch("errm_xsec_pilessdelta",&errm_xsec_pilessdelta,"errm_xsec_pilessdelta/D");
    nfitnt->Branch("errm_xsec_normcccoh",&errm_xsec_normcccoh,"errm_xsec_normcccoh/D");
    nfitnt->Branch("errm_xsec_normnc1piC",&errm_xsec_normnc1piC,"errm_xsec_normnc1piC/D");
    nfitnt->Branch("errm_xsec_normnccoh",&errm_xsec_normnccoh,"errm_xsec_normnccoh/D");
    nfitnt->Branch("errm_xsec_normncoth",&errm_xsec_normncoth,"errm_xsec_normncoth/D");
    nfitnt->Branch("errm_xsec_normnuetonumu",&errm_xsec_normnuetonumu,"errm_xsec_normnuetonumu/D");
    nfitnt->Branch("errm_xsec_normantinutonu",&errm_xsec_normantinutonu,"errm_xsec_normantinutonu/D");
    nfitnt->Branch("errm_skdet_energyscale",&errm_skdet_energyscale,"errm_skdet_energyscale/D");
    nfitnt->Branch("errm_skdet_ccqerc0",&errm_skdet_ccqerc0,"errm_skdet_ccqerc0/D");
    nfitnt->Branch("errm_skdet_ccqerc1",&errm_skdet_ccqerc1,"errm_skdet_ccqerc1/D");
    nfitnt->Branch("errm_skdet_ccqerc2",&errm_skdet_ccqerc2,"errm_skdet_ccqerc2/D");
    nfitnt->Branch("errm_skdet_ccqeoth",&errm_skdet_ccqeoth,"errm_skdet_ccqeoth/D");
    nfitnt->Branch("errm_skdet_ccnqe",&errm_skdet_ccnqe,"errm_skdet_ccnqe/D");
    nfitnt->Branch("errm_skdet_nc",&errm_skdet_nc,"errm_skdet_nc/D");
    nfitnt->Branch("errm_skdet_nuecc",&errm_skdet_nuecc,"errm_skdet_nuecc/D");
    nfitnt->Branch("errm_fsi_000",&errm_fsi_000,"errm_fsi_000/D");
    nfitnt->Branch("errm_fsi_001",&errm_fsi_001,"errm_fsi_001/D");
    nfitnt->Branch("errm_fsi_002",&errm_fsi_002,"errm_fsi_002/D");
    nfitnt->Branch("errm_fsi_003",&errm_fsi_003,"errm_fsi_003/D");
    nfitnt->Branch("errm_fsi_004",&errm_fsi_004,"errm_fsi_004/D");
    nfitnt->Branch("errm_fsi_005",&errm_fsi_005,"errm_fsi_005/D");
    if(conv::kNumOfFSIParams == 15) {
      nfitnt->Branch("errm_fsi_006",&errm_fsi_006,"errm_fsi_006/D");
      nfitnt->Branch("errm_fsi_007",&errm_fsi_007,"errm_fsi_007/D");
      nfitnt->Branch("errm_fsi_008",&errm_fsi_008,"errm_fsi_008/D");
      nfitnt->Branch("errm_fsi_009",&errm_fsi_009,"errm_fsi_009/D");
      nfitnt->Branch("errm_fsi_010",&errm_fsi_010,"errm_fsi_010/D");
      nfitnt->Branch("errm_fsi_011",&errm_fsi_011,"errm_fsi_011/D");
      nfitnt->Branch("errm_fsi_012",&errm_fsi_012,"errm_fsi_012/D");
      nfitnt->Branch("errm_fsi_013",&errm_fsi_013,"errm_fsi_013/D");
      nfitnt->Branch("errm_fsi_014",&errm_fsi_014,"errm_fsi_014/D");
    }
  }
  //parabolic errors
  nfitnt->Branch("errpar_mixing12",&errpar_mixing12,"errpar_mixing12/D");
  nfitnt->Branch("errpar_mixing13",&errpar_mixing13,"errpar_mixing13/D");
  nfitnt->Branch("errpar_mixing23",&errpar_mixing23,"errpar_mixing23/D");
  nfitnt->Branch("errpar_deltacp",&errpar_deltacp,"errpar_deltacp/D");
  nfitnt->Branch("errpar_dm12sq",&errpar_dm12sq,"errpar_dm12sq/D");
  nfitnt->Branch("errpar_dm23sq",&errpar_dm23sq,"errpar_dm23sq/D");
  if(fit_syst_mode == 1) {
    nfitnt->Branch("errpar_banff_maqe",&errpar_banff_maqe,"errpar_banff_maqe/D");
    nfitnt->Branch("errpar_banff_ccqeE0",&errpar_banff_ccqeE0,"errpar_banff_ccqeE0/D");
    nfitnt->Branch("errpar_banff_cc1piE0",&errpar_banff_cc1piE0,"errpar_banff_cc1piE0/D");
    nfitnt->Branch("errpar_xsec_sf",&errpar_xsec_sf,"errpar_xsec_sf/D");
    nfitnt->Branch("errpar_xsec_pilessdelta",&errpar_xsec_pilessdelta,"errpar_xsec_pilessdelta/D");
    nfitnt->Branch("errpar_skdet_ccnqe",&errpar_skdet_ccnqe,"errpar_skdet_ccnqe/D");
    nfitnt->Branch("errpar_skdet_nc",&errpar_skdet_nc,"errpar_skdet_nc/D");
  }
  else if(fit_syst_mode == 2) {
    nfitnt->Branch("errpar_banff_00",&errpar_banff_00,"errpar_banff_00/D");
    nfitnt->Branch("errpar_banff_01",&errpar_banff_01,"errpar_banff_01/D");
    nfitnt->Branch("errpar_banff_02",&errpar_banff_02,"errpar_banff_02/D");
    nfitnt->Branch("errpar_banff_03",&errpar_banff_03,"errpar_banff_03/D");
    nfitnt->Branch("errpar_banff_04",&errpar_banff_04,"errpar_banff_04/D");
    nfitnt->Branch("errpar_banff_05",&errpar_banff_05,"errpar_banff_05/D");
    nfitnt->Branch("errpar_banff_06",&errpar_banff_06,"errpar_banff_06/D");
    nfitnt->Branch("errpar_banff_07",&errpar_banff_07,"errpar_banff_07/D");
    nfitnt->Branch("errpar_banff_08",&errpar_banff_08,"errpar_banff_08/D");
    nfitnt->Branch("errpar_banff_09",&errpar_banff_09,"errpar_banff_09/D");
    nfitnt->Branch("errpar_banff_10",&errpar_banff_10,"errpar_banff_10/D");
    nfitnt->Branch("errpar_banff_11",&errpar_banff_11,"errpar_banff_11/D");
    nfitnt->Branch("errpar_banff_12",&errpar_banff_12,"errpar_banff_12/D");
    nfitnt->Branch("errpar_banff_13",&errpar_banff_13,"errpar_banff_13/D");
    nfitnt->Branch("errpar_banff_14",&errpar_banff_14,"errpar_banff_14/D");
    nfitnt->Branch("errpar_banff_15",&errpar_banff_15,"errpar_banff_15/D");
    nfitnt->Branch("errpar_banff_16",&errpar_banff_16,"errpar_banff_16/D");
    nfitnt->Branch("errpar_banff_17",&errpar_banff_17,"errpar_banff_17/D");
    nfitnt->Branch("errpar_banff_18",&errpar_banff_18,"errpar_banff_18/D");
    nfitnt->Branch("errpar_banff_19",&errpar_banff_19,"errpar_banff_19/D");
    nfitnt->Branch("errpar_banff_20",&errpar_banff_20,"errpar_banff_20/D");
    nfitnt->Branch("errpar_banff_21",&errpar_banff_21,"errpar_banff_21/D");
    nfitnt->Branch("errpar_banff_22",&errpar_banff_22,"errpar_banff_22/D");
    nfitnt->Branch("errpar_banff_23",&errpar_banff_23,"errpar_banff_23/D");
    nfitnt->Branch("errpar_banff_24",&errpar_banff_24,"errpar_banff_24/D");
    nfitnt->Branch("errpar_banff_25",&errpar_banff_25,"errpar_banff_25/D");
    nfitnt->Branch("errpar_banff_26",&errpar_banff_26,"errpar_banff_26/D");
    nfitnt->Branch("errpar_banff_27",&errpar_banff_27,"errpar_banff_27/D");
    nfitnt->Branch("errpar_banff_28",&errpar_banff_28,"errpar_banff_28/D");
    nfitnt->Branch("errpar_banff_29",&errpar_banff_29,"errpar_banff_29/D");
    nfitnt->Branch("errpar_banff_30",&errpar_banff_30,"errpar_banff_30/D");
    nfitnt->Branch("errpar_banff_maqe",&errpar_banff_maqe,"errpar_banff_maqe/D");
    nfitnt->Branch("errpar_banff_mares",&errpar_banff_mares,"errpar_banff_mares/D");
    nfitnt->Branch("errpar_banff_ccqeE0",&errpar_banff_ccqeE0,"errpar_banff_ccqeE0/D");
    nfitnt->Branch("errpar_banff_ccqeE1",&errpar_banff_ccqeE1,"errpar_banff_ccqeE1/D");
    nfitnt->Branch("errpar_banff_ccqeE2",&errpar_banff_ccqeE2,"errpar_banff_ccqeE2/D");
    nfitnt->Branch("errpar_banff_cc1piE0",&errpar_banff_cc1piE0,"errpar_banff_cc1piE0/D");
    nfitnt->Branch("errpar_banff_cc1piE1",&errpar_banff_cc1piE1,"errpar_banff_cc1piE1/D");
    nfitnt->Branch("errpar_banff_nc1pi0",&errpar_banff_nc1pi0,"errpar_banff_nc1pi0/D");
    nfitnt->Branch("errpar_xsec_pf",&errpar_xsec_pf,"errpar_xsec_pf/D");
    nfitnt->Branch("errpar_xsec_sf",&errpar_xsec_sf,"errpar_xsec_sf/D");
    nfitnt->Branch("errpar_xsec_wshape",&errpar_xsec_wshape,"errpar_xsec_wshape/D");
    nfitnt->Branch("errpar_xsec_ccothshape",&errpar_xsec_ccothshape,"errpar_xsec_ccothshape/D");
    nfitnt->Branch("errpar_xsec_cc1pishape",&errpar_xsec_cc1pishape,"errpar_xsec_cc1pishape/D");
    nfitnt->Branch("errpar_xsec_pilessdelta",&errpar_xsec_pilessdelta,"errpar_xsec_pilessdelta/D");
    nfitnt->Branch("errpar_xsec_normcccoh",&errpar_xsec_normcccoh,"errpar_xsec_normcccoh/D");
    nfitnt->Branch("errpar_xsec_normnc1piC",&errpar_xsec_normnc1piC,"errpar_xsec_normnc1piC/D");
    nfitnt->Branch("errpar_xsec_normnccoh",&errpar_xsec_normnccoh,"errpar_xsec_normnccoh/D");
    nfitnt->Branch("errpar_xsec_normncoth",&errpar_xsec_normncoth,"errpar_xsec_normncoth/D");
    nfitnt->Branch("errpar_xsec_normnuetonumu",&errpar_xsec_normnuetonumu,"errpar_xsec_normnuetonumu/D");
    nfitnt->Branch("errpar_xsec_normantinutonu",&errpar_xsec_normantinutonu,"errpar_xsec_normantinutonu/D");
    nfitnt->Branch("errpar_skdet_energyscale",&errpar_skdet_energyscale,"errpar_skdet_energyscale/D");
    nfitnt->Branch("errpar_skdet_ccqerc0",&errpar_skdet_ccqerc0,"errpar_skdet_ccqerc0/D");
    nfitnt->Branch("errpar_skdet_ccqerc1",&errpar_skdet_ccqerc1,"errpar_skdet_ccqerc1/D");
    nfitnt->Branch("errpar_skdet_ccqerc2",&errpar_skdet_ccqerc2,"errpar_skdet_ccqerc2/D");
    nfitnt->Branch("errpar_skdet_ccqeoth",&errpar_skdet_ccqeoth,"errpar_skdet_ccqeoth/D");
    nfitnt->Branch("errpar_skdet_ccnqe",&errpar_skdet_ccnqe,"errpar_skdet_ccnqe/D");
    nfitnt->Branch("errpar_skdet_nc",&errpar_skdet_nc,"errpar_skdet_nc/D");
    nfitnt->Branch("errpar_skdet_nuecc",&errpar_skdet_nuecc,"errpar_skdet_nuecc/D");
    nfitnt->Branch("errpar_fsi_000",&errpar_fsi_000,"errpar_fsi_000/D");
    nfitnt->Branch("errpar_fsi_001",&errpar_fsi_001,"errpar_fsi_001/D");
    nfitnt->Branch("errpar_fsi_002",&errpar_fsi_002,"errpar_fsi_002/D");
    nfitnt->Branch("errpar_fsi_003",&errpar_fsi_003,"errpar_fsi_003/D");
    nfitnt->Branch("errpar_fsi_004",&errpar_fsi_004,"errpar_fsi_004/D");
    nfitnt->Branch("errpar_fsi_005",&errpar_fsi_005,"errpar_fsi_005/D");
    if(conv::kNumOfFSIParams == 15) {
      nfitnt->Branch("errpar_fsi_006",&errpar_fsi_006,"errpar_fsi_006/D");
      nfitnt->Branch("errpar_fsi_007",&errpar_fsi_007,"errpar_fsi_007/D");
      nfitnt->Branch("errpar_fsi_008",&errpar_fsi_008,"errpar_fsi_008/D");
      nfitnt->Branch("errpar_fsi_009",&errpar_fsi_009,"errpar_fsi_009/D");
      nfitnt->Branch("errpar_fsi_010",&errpar_fsi_010,"errpar_fsi_010/D");
      nfitnt->Branch("errpar_fsi_011",&errpar_fsi_011,"errpar_fsi_011/D");
      nfitnt->Branch("errpar_fsi_012",&errpar_fsi_012,"errpar_fsi_012/D");
      nfitnt->Branch("errpar_fsi_013",&errpar_fsi_013,"errpar_fsi_013/D");
      nfitnt->Branch("errpar_fsi_014",&errpar_fsi_014,"errpar_fsi_014/D");
    }
  }

  // Loop over toy-MC experiments
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
      default:
	cerr << "Invalid systematic mode " << syst_mode << ". Exiting" << endl;
	exit(1);
	break;
      }
      skexp->PrintSystematicVariations();
      // save the input systematics to the tree
      if(fit_syst_mode == 1) {
	input_banff_maqe = skexp->GetBANFF()->CurrentValue(conv::kNuisBANFFMAQE);
	input_banff_ccqeE0 = skexp->GetBANFF()->CurrentValue(conv::kNuisBANFFCCQEE0);
	input_banff_cc1piE0 = skexp->GetBANFF()->CurrentValue(conv::kNuisBANFFCC1piE0);
	input_xsec_sf = skexp->GetNuXSec()->CurrentValue(conv::kNuisXSecSpecFunc);
	input_xsec_pilessdelta = skexp->GetNuXSec()->CurrentValue(conv::kNuisXSecENuShape);
	input_skdet_ccnqe = skexp->GetSuperkDetector()->CurrentValue(conv::kNuisSuperkEff1RmuCCnonQE);
	input_skdet_nc = skexp->GetSuperkDetector()->CurrentValue(conv::kNuisSuperkEff1RmuNC);
      }
      else if(fit_syst_mode == 2) {
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
	input_banff_16 = skexp->GetBANFF()->CurrentValue(conv::kNuisBANFF16);
	input_banff_17 = skexp->GetBANFF()->CurrentValue(conv::kNuisBANFF17);
	input_banff_18 = skexp->GetBANFF()->CurrentValue(conv::kNuisBANFF18);
	input_banff_19 = skexp->GetBANFF()->CurrentValue(conv::kNuisBANFF19);
	input_banff_20 = skexp->GetBANFF()->CurrentValue(conv::kNuisBANFF20);
	input_banff_21 = skexp->GetBANFF()->CurrentValue(conv::kNuisBANFF21);
	input_banff_22 = skexp->GetBANFF()->CurrentValue(conv::kNuisBANFF22);
	input_banff_23 = skexp->GetBANFF()->CurrentValue(conv::kNuisBANFF23);
	input_banff_24 = skexp->GetBANFF()->CurrentValue(conv::kNuisBANFF24);
	input_banff_25 = skexp->GetBANFF()->CurrentValue(conv::kNuisBANFF25);
	input_banff_26 = skexp->GetBANFF()->CurrentValue(conv::kNuisBANFF26);
	input_banff_27 = skexp->GetBANFF()->CurrentValue(conv::kNuisBANFF27);
	input_banff_28 = skexp->GetBANFF()->CurrentValue(conv::kNuisBANFF28);
	input_banff_29 = skexp->GetBANFF()->CurrentValue(conv::kNuisBANFF29);
	input_banff_30 = skexp->GetBANFF()->CurrentValue(conv::kNuisBANFF30);
	input_banff_maqe = skexp->GetBANFF()->CurrentValue(conv::kNuisBANFFMAQE);
	input_banff_mares = skexp->GetBANFF()->CurrentValue(conv::kNuisBANFFMARES);
	input_banff_ccqeE0 = skexp->GetBANFF()->CurrentValue(conv::kNuisBANFFCCQEE0);
	input_banff_ccqeE1 = skexp->GetBANFF()->CurrentValue(conv::kNuisBANFFCCQEE1);
	input_banff_ccqeE2 = skexp->GetBANFF()->CurrentValue(conv::kNuisBANFFCCQEE2);
	input_banff_cc1piE0 = skexp->GetBANFF()->CurrentValue(conv::kNuisBANFFCC1piE0);
	input_banff_cc1piE1 = skexp->GetBANFF()->CurrentValue(conv::kNuisBANFFCC1piE1);
	input_banff_nc1pi0 = skexp->GetBANFF()->CurrentValue(conv::kNuisBANFFNC1pi0);
	input_xsec_pf = skexp->GetNuXSec()->CurrentValue(conv::kNuisXSecFermiMom);
	input_xsec_sf = skexp->GetNuXSec()->CurrentValue(conv::kNuisXSecSpecFunc);
	input_xsec_wshape = skexp->GetNuXSec()->CurrentValue(conv::kNuisXSecWshape);
	input_xsec_ccothshape = skexp->GetNuXSec()->CurrentValue(conv::kNuisXSecShapeCCoth);
	input_xsec_cc1pishape = skexp->GetNuXSec()->CurrentValue(conv::kNuisXSecPiLessDelta);
	input_xsec_pilessdelta = skexp->GetNuXSec()->CurrentValue(conv::kNuisXSecENuShape);
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
	if(conv::kNumOfFSIParams == 15) {
	  input_fsi_006 = skexp->GetFSI()->CurrentValue(conv::kNuisFSI006);
	  input_fsi_007 = skexp->GetFSI()->CurrentValue(conv::kNuisFSI007);
	  input_fsi_008 = skexp->GetFSI()->CurrentValue(conv::kNuisFSI008);
	  input_fsi_009 = skexp->GetFSI()->CurrentValue(conv::kNuisFSI009);
	  input_fsi_010 = skexp->GetFSI()->CurrentValue(conv::kNuisFSI010);
	  input_fsi_011 = skexp->GetFSI()->CurrentValue(conv::kNuisFSI011);
	  input_fsi_012 = skexp->GetFSI()->CurrentValue(conv::kNuisFSI012);
	  input_fsi_013 = skexp->GetFSI()->CurrentValue(conv::kNuisFSI013);
	  input_fsi_014 = skexp->GetFSI()->CurrentValue(conv::kNuisFSI014);
	}
      }
      //update the FitParamList with the fluctuated nuisance parameters
      if(fit_syst_mode == 1) {
	if(fparm->Exists(conv::kNuisBANFFMAQE)) fparm->Update(conv::kNuisBANFFMAQE, input_banff_maqe);
	if(fparm->Exists(conv::kNuisBANFFCCQEE0)) fparm->Update(conv::kNuisBANFFCCQEE0, input_banff_ccqeE0);
	if(fparm->Exists(conv::kNuisBANFFCC1piE0)) fparm->Update(conv::kNuisBANFFCC1piE0, input_banff_cc1piE0);
	if(fparm->Exists(conv::kNuisXSecSpecFunc)) fparm->Update(conv::kNuisXSecSpecFunc, input_xsec_sf);
	if(fparm->Exists(conv::kNuisXSecPiLessDelta)) fparm->Update(conv::kNuisXSecPiLessDelta, input_xsec_cc1pishape);
	if(fparm->Exists(conv::kNuisSuperkEff1RmuCCnonQE)) fparm->Update(conv::kNuisSuperkEff1RmuCCnonQE, input_skdet_ccnqe);
	if(fparm->Exists(conv::kNuisSuperkEff1RmuNC)) fparm->Update(conv::kNuisSuperkEff1RmuNC, input_skdet_nc);
      }
      else if(fit_syst_mode == 2) {
	if(fparm->Exists(conv::kNuisBANFF00)) fparm->Update(conv::kNuisBANFF00, input_banff_00);
	if(fparm->Exists(conv::kNuisBANFF01)) fparm->Update(conv::kNuisBANFF01, input_banff_01);
	if(fparm->Exists(conv::kNuisBANFF02)) fparm->Update(conv::kNuisBANFF02, input_banff_02);
	if(fparm->Exists(conv::kNuisBANFF03)) fparm->Update(conv::kNuisBANFF03, input_banff_03);
	if(fparm->Exists(conv::kNuisBANFF04)) fparm->Update(conv::kNuisBANFF04, input_banff_04);
	if(fparm->Exists(conv::kNuisBANFF05)) fparm->Update(conv::kNuisBANFF05, input_banff_05);
	if(fparm->Exists(conv::kNuisBANFF06)) fparm->Update(conv::kNuisBANFF06, input_banff_06);
	if(fparm->Exists(conv::kNuisBANFF07)) fparm->Update(conv::kNuisBANFF07, input_banff_07);
	if(fparm->Exists(conv::kNuisBANFF08)) fparm->Update(conv::kNuisBANFF08, input_banff_08);
	if(fparm->Exists(conv::kNuisBANFF09)) fparm->Update(conv::kNuisBANFF09, input_banff_09);
	if(fparm->Exists(conv::kNuisBANFF10)) fparm->Update(conv::kNuisBANFF10, input_banff_10);
	if(fparm->Exists(conv::kNuisBANFF11)) fparm->Update(conv::kNuisBANFF11, input_banff_11);
	if(fparm->Exists(conv::kNuisBANFF12)) fparm->Update(conv::kNuisBANFF12, input_banff_12);
	if(fparm->Exists(conv::kNuisBANFF13)) fparm->Update(conv::kNuisBANFF13, input_banff_13);
	if(fparm->Exists(conv::kNuisBANFF14)) fparm->Update(conv::kNuisBANFF14, input_banff_14);
	if(fparm->Exists(conv::kNuisBANFF15)) fparm->Update(conv::kNuisBANFF15, input_banff_15);
	if(fparm->Exists(conv::kNuisBANFF16)) fparm->Update(conv::kNuisBANFF16, input_banff_16);
	if(fparm->Exists(conv::kNuisBANFF17)) fparm->Update(conv::kNuisBANFF17, input_banff_17);
	if(fparm->Exists(conv::kNuisBANFF18)) fparm->Update(conv::kNuisBANFF18, input_banff_18);
	if(fparm->Exists(conv::kNuisBANFF19)) fparm->Update(conv::kNuisBANFF19, input_banff_19);
	if(fparm->Exists(conv::kNuisBANFF20)) fparm->Update(conv::kNuisBANFF20, input_banff_20);
	if(fparm->Exists(conv::kNuisBANFF21)) fparm->Update(conv::kNuisBANFF21, input_banff_21);
	if(fparm->Exists(conv::kNuisBANFF22)) fparm->Update(conv::kNuisBANFF22, input_banff_22);
	if(fparm->Exists(conv::kNuisBANFF23)) fparm->Update(conv::kNuisBANFF23, input_banff_23);
	if(fparm->Exists(conv::kNuisBANFF24)) fparm->Update(conv::kNuisBANFF24, input_banff_24);
	if(fparm->Exists(conv::kNuisBANFF25)) fparm->Update(conv::kNuisBANFF25, input_banff_25);
	if(fparm->Exists(conv::kNuisBANFF26)) fparm->Update(conv::kNuisBANFF26, input_banff_26);
	if(fparm->Exists(conv::kNuisBANFF27)) fparm->Update(conv::kNuisBANFF27, input_banff_27);
	if(fparm->Exists(conv::kNuisBANFF28)) fparm->Update(conv::kNuisBANFF28, input_banff_28);
	if(fparm->Exists(conv::kNuisBANFF29)) fparm->Update(conv::kNuisBANFF29, input_banff_29);
	if(fparm->Exists(conv::kNuisBANFF30)) fparm->Update(conv::kNuisBANFF30, input_banff_30);
	if(fparm->Exists(conv::kNuisBANFFMAQE)) fparm->Update(conv::kNuisBANFFMAQE, input_banff_maqe);
	if(fparm->Exists(conv::kNuisBANFFMARES)) fparm->Update(conv::kNuisBANFFMARES, input_banff_mares);
	if(fparm->Exists(conv::kNuisBANFFCCQEE0)) fparm->Update(conv::kNuisBANFFCCQEE0, input_banff_ccqeE0);
	if(fparm->Exists(conv::kNuisBANFFCCQEE1)) fparm->Update(conv::kNuisBANFFCCQEE1, input_banff_ccqeE1);
	if(fparm->Exists(conv::kNuisBANFFCCQEE2)) fparm->Update(conv::kNuisBANFFCCQEE2, input_banff_ccqeE2);
	if(fparm->Exists(conv::kNuisBANFFCC1piE0)) fparm->Update(conv::kNuisBANFFCC1piE0, input_banff_cc1piE0);
	if(fparm->Exists(conv::kNuisBANFFCC1piE1)) fparm->Update(conv::kNuisBANFFCC1piE1, input_banff_cc1piE1);
	if(fparm->Exists(conv::kNuisBANFFNC1pi0)) fparm->Update(conv::kNuisBANFFNC1pi0, input_banff_nc1pi0);
	if(fparm->Exists(conv::kNuisXSecFermiMom)) fparm->Update(conv::kNuisXSecFermiMom, input_xsec_pf);
	if(fparm->Exists(conv::kNuisXSecSpecFunc)) fparm->Update(conv::kNuisXSecSpecFunc, input_xsec_sf);
	if(fparm->Exists(conv::kNuisXSecWshape)) fparm->Update(conv::kNuisXSecWshape, input_xsec_wshape);
	if(fparm->Exists(conv::kNuisXSecShapeCCoth)) fparm->Update(conv::kNuisXSecShapeCCoth, input_xsec_ccothshape);
	if(fparm->Exists(conv::kNuisXSecPiLessDelta)) fparm->Update(conv::kNuisXSecPiLessDelta, input_xsec_cc1pishape);
	if(fparm->Exists(conv::kNuisXSecENuShape)) fparm->Update(conv::kNuisXSecENuShape, input_xsec_pilessdelta);
	if(fparm->Exists(conv::kNuisXSecNormCCcoh)) fparm->Update(conv::kNuisXSecNormCCcoh, input_xsec_normcccoh);
	if(fparm->Exists(conv::kNuisXSecNormNC1piC)) fparm->Update(conv::kNuisXSecNormNC1piC, input_xsec_normnc1piC);
	if(fparm->Exists(conv::kNuisXSecNormNCcoh)) fparm->Update(conv::kNuisXSecNormNCcoh, input_xsec_normnccoh);
	if(fparm->Exists(conv::kNuisXSecNormNCoth)) fparm->Update(conv::kNuisXSecNormNCoth, input_xsec_normncoth);
	if(fparm->Exists(conv::kNuisXSecNormNuEToNuMu)) fparm->Update(conv::kNuisXSecNormNuEToNuMu, input_xsec_normnuetonumu);
	if(fparm->Exists(conv::kNuisXSecNormAntiNuToNu)) fparm->Update(conv::kNuisXSecNormAntiNuToNu, input_xsec_normantinutonu);
	if(fparm->Exists(conv::kNuisSuperkEnergyScale)) fparm->Update(conv::kNuisSuperkEnergyScale, input_skdet_energyscale);
	if(fparm->Exists(conv::kNuisSuperkEff1RmuCCQErc0u)) fparm->Update(conv::kNuisSuperkEff1RmuCCQErc0u, input_skdet_ccqerc0);
	if(fparm->Exists(conv::kNuisSuperkEff1RmuCCQErc1u)) fparm->Update(conv::kNuisSuperkEff1RmuCCQErc1u, input_skdet_ccqerc1);
	if(fparm->Exists(conv::kNuisSuperkEff1RmuCCQErc2u)) fparm->Update(conv::kNuisSuperkEff1RmuCCQErc2u, input_skdet_ccqerc2);
	if(fparm->Exists(conv::kNuisSuperkEff1RmuCCQEoth)) fparm->Update(conv::kNuisSuperkEff1RmuCCQEoth, input_skdet_ccqeoth);
	if(fparm->Exists(conv::kNuisSuperkEff1RmuCCnonQE)) fparm->Update(conv::kNuisSuperkEff1RmuCCnonQE, input_skdet_ccnqe);
	if(fparm->Exists(conv::kNuisSuperkEff1RmuNC)) fparm->Update(conv::kNuisSuperkEff1RmuNC, input_skdet_nc);
	if(fparm->Exists(conv::kNuisSuperkEff1RmuNueCC)) fparm->Update(conv::kNuisSuperkEff1RmuNueCC, input_skdet_nuecc);
	if(fparm->Exists(conv::kNuisFSI000)) fparm->Update(conv::kNuisFSI000, input_fsi_000);
	if(fparm->Exists(conv::kNuisFSI001)) fparm->Update(conv::kNuisFSI001, input_fsi_001);
	if(fparm->Exists(conv::kNuisFSI002)) fparm->Update(conv::kNuisFSI002, input_fsi_002);
	if(fparm->Exists(conv::kNuisFSI003)) fparm->Update(conv::kNuisFSI003, input_fsi_003);
	if(fparm->Exists(conv::kNuisFSI004)) fparm->Update(conv::kNuisFSI004, input_fsi_004);
	if(fparm->Exists(conv::kNuisFSI005)) fparm->Update(conv::kNuisFSI005, input_fsi_005);
	if(conv::kNumOfFSIParams == 15) {
	  if(fparm->Exists(conv::kNuisFSI006)) fparm->Update(conv::kNuisFSI006, input_fsi_006);
	  if(fparm->Exists(conv::kNuisFSI007)) fparm->Update(conv::kNuisFSI007, input_fsi_007);
	  if(fparm->Exists(conv::kNuisFSI008)) fparm->Update(conv::kNuisFSI008, input_fsi_008);
	  if(fparm->Exists(conv::kNuisFSI009)) fparm->Update(conv::kNuisFSI009, input_fsi_009);
	  if(fparm->Exists(conv::kNuisFSI010)) fparm->Update(conv::kNuisFSI010, input_fsi_010);
	  if(fparm->Exists(conv::kNuisFSI011)) fparm->Update(conv::kNuisFSI011, input_fsi_011);
	  if(fparm->Exists(conv::kNuisFSI012)) fparm->Update(conv::kNuisFSI012, input_fsi_012);
	  if(fparm->Exists(conv::kNuisFSI013)) fparm->Update(conv::kNuisFSI013, input_fsi_013);
	  if(fparm->Exists(conv::kNuisFSI014)) fparm->Update(conv::kNuisFSI014, input_fsi_014);
	}
      }
      // generate mock data using the current PDF (incl. systematic variations) 
      // + statistical fluctuations if requested (Poisson fluctuation in each bin)
      cout << "- Generating toy-MC..." << endl;
      skobs->CreateMockData(*skexp, *fparm, stat_fluctuations);
      // reset all systematics so that the PDF used in the fitter is the nominal one
      cout << "- Reseting all systematic params in pdf..." << endl;
      skexp->ResetSystematicVariations();
      //skexp->PrintSystematicVariations();
    } else {       
      // generate mock data using the current PDF (without systematic variations) 
      // with statistical fluctuations if requested (Poisson fluctuation in each bin)
      cout << "- Generating toy-MC..." << endl;
      skobs->CreateMockData(*skexp, *fparm, stat_fluctuations); 
    }

    // Get nominal number of events
    nsk_nom = skexp->CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike)->Integral();

    // Get fitter
    FitInputs * fit_inputs = new FitInputs(skobs,skexp,fparm);
    FitterABC * fitter = new Superk1RmuErecoFitter(fit_inputs,fit_option);

    // Chisq for true (input) params
    true_chisq = fit_inputs->Chisquare();

    // Print-out true params
    fparm->ShowOnlyCurrValuesInNextPrint();
    cout << "- True oscillation params (below) yield chisq = " << true_chisq
	 << *fparm;

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

    // Get & store best-fit params
    // Init best-fit values
    bestfit_mixing12 = -1;
    bestfit_mixing13 = -1;
    bestfit_mixing23 = -1;
    bestfit_deltacp  = -1;
    bestfit_dm12sq   = -1;
    bestfit_dm23sq   = -1;
    if(fit_syst_mode == 1) {
      bestfit_banff_maqe = -1;
      bestfit_banff_ccqeE0 = -1;
      bestfit_banff_cc1piE0 = -1;
      bestfit_xsec_sf = -1;
      bestfit_xsec_pilessdelta = -1;
      bestfit_skdet_ccnqe = -1;
      bestfit_skdet_nc = -1;
    }
    else if(fit_syst_mode == 2) {
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
      bestfit_banff_16 = -1;
      bestfit_banff_17 = -1;
      bestfit_banff_18 = -1;
      bestfit_banff_19 = -1;
      bestfit_banff_20 = -1;
      bestfit_banff_21 = -1;
      bestfit_banff_22 = -1;
      bestfit_banff_23 = -1;
      bestfit_banff_24 = -1;
      bestfit_banff_25 = -1;
      bestfit_banff_26 = -1;
      bestfit_banff_27 = -1;
      bestfit_banff_28 = -1;
      bestfit_banff_29 = -1;
      bestfit_banff_30 = -1;
      bestfit_banff_maqe = -1;
      bestfit_banff_mares = -1;
      bestfit_banff_ccqeE0 = -1;
      bestfit_banff_ccqeE1 = -1;
      bestfit_banff_ccqeE2 = -1;
      bestfit_banff_cc1piE0 = -1;
      bestfit_banff_cc1piE1 = -1;
      bestfit_banff_nc1pi0 = -1;
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
      if(conv::kNumOfFSIParams == 15) {
	bestfit_fsi_006 = -1;
	bestfit_fsi_007 = -1;
	bestfit_fsi_008 = -1;
	bestfit_fsi_009 = -1;
	bestfit_fsi_010 = -1;
	bestfit_fsi_011 = -1;
	bestfit_fsi_012 = -1;
	bestfit_fsi_013 = -1;
	bestfit_fsi_014 = -1;
      }
    }
    // Init plus errors
    errp_mixing12 = -1;
    errp_mixing13 = -1;
    errp_mixing23 = -1;
    errp_deltacp  = -1;
    errp_dm12sq   = -1;
    errp_dm23sq   = -1;
    if(fit_syst_mode == 1) {
      errp_banff_maqe = -1;
      errp_banff_ccqeE0 = -1;
      errp_banff_cc1piE0 = -1;
      errp_xsec_sf = -1;
      errp_xsec_pilessdelta = -1;
      errp_skdet_ccnqe = -1;
      errp_skdet_nc = -1;
    }
    else if(fit_syst_mode == 2) {
      errp_banff_00 = -1;
      errp_banff_01 = -1;
      errp_banff_02 = -1;
      errp_banff_03 = -1;
      errp_banff_04 = -1;
      errp_banff_05 = -1;
      errp_banff_06 = -1;
      errp_banff_07 = -1;
      errp_banff_08 = -1;
      errp_banff_09 = -1;
      errp_banff_10 = -1;
      errp_banff_11 = -1;
      errp_banff_12 = -1;
      errp_banff_13 = -1;
      errp_banff_14 = -1;
      errp_banff_15 = -1;
      errp_banff_16 = -1;
      errp_banff_17 = -1;
      errp_banff_18 = -1;
      errp_banff_19 = -1;
      errp_banff_20 = -1;
      errp_banff_21 = -1;
      errp_banff_22 = -1;
      errp_banff_23 = -1;
      errp_banff_24 = -1;
      errp_banff_25 = -1;
      errp_banff_26 = -1;
      errp_banff_27 = -1;
      errp_banff_28 = -1;
      errp_banff_29 = -1;
      errp_banff_30 = -1;
      errp_banff_maqe = -1;
      errp_banff_mares = -1;
      errp_banff_ccqeE0 = -1;
      errp_banff_ccqeE1 = -1;
      errp_banff_ccqeE2 = -1;
      errp_banff_cc1piE0 = -1;
      errp_banff_cc1piE1 = -1;
      errp_banff_nc1pi0 = -1;
      errp_xsec_pf = -1;
      errp_xsec_sf = -1;
      errp_xsec_wshape = -1;
      errp_xsec_ccothshape = -1;
      errp_xsec_cc1pishape = -1;
      errp_xsec_pilessdelta = -1;
      errp_xsec_normcccoh = -1;
      errp_xsec_normnc1piC = -1;
      errp_xsec_normnccoh = -1;
      errp_xsec_normncoth = -1;
      errp_xsec_normnuetonumu = -1;
      errp_xsec_normantinutonu = -1;
      errp_skdet_energyscale = -1;
      errp_skdet_ccqerc0 = -1;
      errp_skdet_ccqerc1 = -1;
      errp_skdet_ccqerc2 = -1;
      errp_skdet_ccqeoth = -1;
      errp_skdet_ccnqe = -1;
      errp_skdet_nc = -1;
      errp_skdet_nuecc = -1;
      errp_fsi_000 = -1;
      errp_fsi_001 = -1;
      errp_fsi_002 = -1;
      errp_fsi_003 = -1;
      errp_fsi_004 = -1;
      errp_fsi_005 = -1;
      if(conv::kNumOfFSIParams == 15) {
	errp_fsi_006 = -1;
	errp_fsi_007 = -1;
	errp_fsi_008 = -1;
	errp_fsi_009 = -1;
	errp_fsi_010 = -1;
	errp_fsi_011 = -1;
	errp_fsi_012 = -1;
	errp_fsi_013 = -1;
	errp_fsi_014 = -1;
      }
    }
    // Init minus errors
    errm_mixing12 = -1;
    errm_mixing13 = -1;
    errm_mixing23 = -1;
    errm_deltacp  = -1;
    errm_dm12sq   = -1;
    errm_dm23sq   = -1;
    if(fit_syst_mode == 1) {
      errm_banff_maqe = -1;
      errm_banff_ccqeE0 = -1;
      errm_banff_cc1piE0 = -1;
      errm_xsec_sf = -1;
      errm_xsec_pilessdelta = -1;
      errm_skdet_ccnqe = -1;
      errm_skdet_nc = -1;
    }
    else if(fit_syst_mode == 2) {
      errm_banff_00 = -1;
      errm_banff_01 = -1;
      errm_banff_02 = -1;
      errm_banff_03 = -1;
      errm_banff_04 = -1;
      errm_banff_05 = -1;
      errm_banff_06 = -1;
      errm_banff_07 = -1;
      errm_banff_08 = -1;
      errm_banff_09 = -1;
      errm_banff_10 = -1;
      errm_banff_11 = -1;
      errm_banff_12 = -1;
      errm_banff_13 = -1;
      errm_banff_14 = -1;
      errm_banff_15 = -1;
      errm_banff_16 = -1;
      errm_banff_17 = -1;
      errm_banff_18 = -1;
      errm_banff_19 = -1;
      errm_banff_20 = -1;
      errm_banff_21 = -1;
      errm_banff_22 = -1;
      errm_banff_23 = -1;
      errm_banff_24 = -1;
      errm_banff_25 = -1;
      errm_banff_26 = -1;
      errm_banff_27 = -1;
      errm_banff_28 = -1;
      errm_banff_29 = -1;
      errm_banff_30 = -1;
      errm_banff_maqe = -1;
      errm_banff_mares = -1;
      errm_banff_ccqeE0 = -1;
      errm_banff_ccqeE1 = -1;
      errm_banff_ccqeE2 = -1;
      errm_banff_cc1piE0 = -1;
      errm_banff_cc1piE1 = -1;
      errm_banff_nc1pi0 = -1;
      errm_xsec_pf = -1;
      errm_xsec_sf = -1;
      errm_xsec_wshape = -1;
      errm_xsec_ccothshape = -1;
      errm_xsec_cc1pishape = -1;
      errm_xsec_pilessdelta = -1;
      errm_xsec_normcccoh = -1;
      errm_xsec_normnc1piC = -1;
      errm_xsec_normnccoh = -1;
      errm_xsec_normncoth = -1;
      errm_xsec_normnuetonumu = -1;
      errm_xsec_normantinutonu = -1;
      errm_skdet_energyscale = -1;
      errm_skdet_ccqerc0 = -1;
      errm_skdet_ccqerc1 = -1;
      errm_skdet_ccqerc2 = -1;
      errm_skdet_ccqeoth = -1;
      errm_skdet_ccnqe = -1;
      errm_skdet_nc = -1;
      errm_skdet_nuecc = -1;
      errm_fsi_000 = -1;
      errm_fsi_001 = -1;
      errm_fsi_002 = -1;
      errm_fsi_003 = -1;
      errm_fsi_004 = -1;
      errm_fsi_005 = -1;
      if(conv::kNumOfFSIParams == 15) {
	errm_fsi_006 = -1;
	errm_fsi_007 = -1;
	errm_fsi_008 = -1;
	errm_fsi_009 = -1;
	errm_fsi_010 = -1;
	errm_fsi_011 = -1;
	errm_fsi_012 = -1;
	errm_fsi_013 = -1;
	errm_fsi_014 = -1;
      }
    }
    // Init parabolic errors
    errpar_mixing12 = -1;
    errpar_mixing13 = -1;
    errpar_mixing23 = -1;
    errpar_deltacp  = -1;
    errpar_dm12sq   = -1;
    errpar_dm23sq   = -1;
    if(fit_syst_mode == 1) {
      errpar_banff_maqe = -1;
      errpar_banff_ccqeE0 = -1;
      errpar_banff_cc1piE0 = -1;
      errpar_xsec_sf = -1;
      errpar_xsec_pilessdelta = -1;
      errpar_skdet_ccnqe = -1;
      errpar_skdet_nc = -1;
    }
    else if(fit_syst_mode == 2) {
      errpar_banff_00 = -1;
      errpar_banff_01 = -1;
      errpar_banff_02 = -1;
      errpar_banff_03 = -1;
      errpar_banff_04 = -1;
      errpar_banff_05 = -1;
      errpar_banff_06 = -1;
      errpar_banff_07 = -1;
      errpar_banff_08 = -1;
      errpar_banff_09 = -1;
      errpar_banff_10 = -1;
      errpar_banff_11 = -1;
      errpar_banff_12 = -1;
      errpar_banff_13 = -1;
      errpar_banff_14 = -1;
      errpar_banff_15 = -1;
      errpar_banff_16 = -1;
      errpar_banff_17 = -1;
      errpar_banff_18 = -1;
      errpar_banff_19 = -1;
      errpar_banff_20 = -1;
      errpar_banff_21 = -1;
      errpar_banff_22 = -1;
      errpar_banff_23 = -1;
      errpar_banff_24 = -1;
      errpar_banff_25 = -1;
      errpar_banff_26 = -1;
      errpar_banff_27 = -1;
      errpar_banff_28 = -1;
      errpar_banff_29 = -1;
      errpar_banff_30 = -1;
      errpar_banff_maqe = -1;
      errpar_banff_mares = -1;
      errpar_banff_ccqeE0 = -1;
      errpar_banff_ccqeE1 = -1;
      errpar_banff_ccqeE2 = -1;
      errpar_banff_cc1piE0 = -1;
      errpar_banff_cc1piE1 = -1;
      errpar_banff_nc1pi0 = -1;
      errpar_xsec_pf = -1;
      errpar_xsec_sf = -1;
      errpar_xsec_wshape = -1;
      errpar_xsec_ccothshape = -1;
      errpar_xsec_cc1pishape = -1;
      errpar_xsec_pilessdelta = -1;
      errpar_xsec_normcccoh = -1;
      errpar_xsec_normnc1piC = -1;
      errpar_xsec_normnccoh = -1;
      errpar_xsec_normncoth = -1;
      errpar_xsec_normnuetonumu = -1;
      errpar_xsec_normantinutonu = -1;
      errpar_skdet_energyscale = -1;
      errpar_skdet_ccqerc0 = -1;
      errpar_skdet_ccqerc1 = -1;
      errpar_skdet_ccqerc2 = -1;
      errpar_skdet_ccqeoth = -1;
      errpar_skdet_ccnqe = -1;
      errpar_skdet_nc = -1;
      errpar_skdet_nuecc = -1;
      errpar_fsi_000 = -1;
      errpar_fsi_001 = -1;
      errpar_fsi_002 = -1;
      errpar_fsi_003 = -1;
      errpar_fsi_004 = -1;
      errpar_fsi_005 = -1;
      if(conv::kNumOfFSIParams == 15) {
	errpar_fsi_006 = -1;
	errpar_fsi_007 = -1;
	errpar_fsi_008 = -1;
	errpar_fsi_009 = -1;
	errpar_fsi_010 = -1;
	errpar_fsi_011 = -1;
	errpar_fsi_012 = -1;
	errpar_fsi_013 = -1;
	errpar_fsi_014 = -1;
      }
    }
    //get the bestfit values
    if(hyp == kHypOsc3F) {
      bestfit_mixing12 = fparm->BestFitValue(kPhypMixing12);
      bestfit_mixing13 = fparm->BestFitValue(kPhypMixing13);
      bestfit_mixing23 = fparm->BestFitValue(kPhypMixing23);
      bestfit_deltacp  = fparm->BestFitValue(kPhypDeltaCP);
      bestfit_dm12sq   = fparm->BestFitValue(kPhypDm12sq);
      bestfit_dm23sq   = fparm->BestFitValue(kPhypDm23sq);
    } else
      if(hyp == kHypNuMuDisapp2F) {
        bestfit_mixing23 = fparm->BestFitValue(kPhypMixing23);
        bestfit_dm23sq   = fparm->BestFitValue(kPhypDm23sq);
      }
    if(fit_syst_mode == 1) {
      if(fparm->Exists(conv::kNuisBANFFMAQE)) bestfit_banff_maqe = fparm->BestFitValue(conv::kNuisBANFFMAQE);
      if(fparm->Exists(conv::kNuisBANFFCCQEE0)) bestfit_banff_ccqeE0 = fparm->BestFitValue(conv::kNuisBANFFCCQEE0);
      if(fparm->Exists(conv::kNuisBANFFCC1piE0)) bestfit_banff_cc1piE0 = fparm->BestFitValue(conv::kNuisBANFFCC1piE0);
      if(fparm->Exists(conv::kNuisXSecSpecFunc)) bestfit_xsec_sf = fparm->BestFitValue(conv::kNuisXSecSpecFunc);
      if(fparm->Exists(conv::kNuisXSecPiLessDelta)) bestfit_xsec_cc1pishape = fparm->BestFitValue(conv::kNuisXSecPiLessDelta);
      if(fparm->Exists(conv::kNuisSuperkEff1RmuCCnonQE)) bestfit_skdet_ccnqe = fparm->BestFitValue(conv::kNuisSuperkEff1RmuCCnonQE);
      if(fparm->Exists(conv::kNuisSuperkEff1RmuNC)) bestfit_skdet_nc = fparm->BestFitValue(conv::kNuisSuperkEff1RmuNC);
    }
    else if(fit_syst_mode == 2) {
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
      if(fparm->Exists(conv::kNuisBANFF16)) bestfit_banff_16 = fparm->BestFitValue(conv::kNuisBANFF16);
      if(fparm->Exists(conv::kNuisBANFF17)) bestfit_banff_17 = fparm->BestFitValue(conv::kNuisBANFF17);
      if(fparm->Exists(conv::kNuisBANFF18)) bestfit_banff_18 = fparm->BestFitValue(conv::kNuisBANFF18);
      if(fparm->Exists(conv::kNuisBANFF19)) bestfit_banff_19 = fparm->BestFitValue(conv::kNuisBANFF19);
      if(fparm->Exists(conv::kNuisBANFF20)) bestfit_banff_20 = fparm->BestFitValue(conv::kNuisBANFF20);
      if(fparm->Exists(conv::kNuisBANFF21)) bestfit_banff_21 = fparm->BestFitValue(conv::kNuisBANFF21);
      if(fparm->Exists(conv::kNuisBANFF22)) bestfit_banff_22 = fparm->BestFitValue(conv::kNuisBANFF22);
      if(fparm->Exists(conv::kNuisBANFF23)) bestfit_banff_23 = fparm->BestFitValue(conv::kNuisBANFF23);
      if(fparm->Exists(conv::kNuisBANFF24)) bestfit_banff_24 = fparm->BestFitValue(conv::kNuisBANFF24);
      if(fparm->Exists(conv::kNuisBANFF25)) bestfit_banff_25 = fparm->BestFitValue(conv::kNuisBANFF25);
      if(fparm->Exists(conv::kNuisBANFF26)) bestfit_banff_26 = fparm->BestFitValue(conv::kNuisBANFF26);
      if(fparm->Exists(conv::kNuisBANFF27)) bestfit_banff_27 = fparm->BestFitValue(conv::kNuisBANFF27);
      if(fparm->Exists(conv::kNuisBANFF28)) bestfit_banff_28 = fparm->BestFitValue(conv::kNuisBANFF28);
      if(fparm->Exists(conv::kNuisBANFF29)) bestfit_banff_29 = fparm->BestFitValue(conv::kNuisBANFF29);
      if(fparm->Exists(conv::kNuisBANFF30)) bestfit_banff_30 = fparm->BestFitValue(conv::kNuisBANFF30);
      if(fparm->Exists(conv::kNuisBANFFMAQE)) bestfit_banff_maqe = fparm->BestFitValue(conv::kNuisBANFFMAQE);
      if(fparm->Exists(conv::kNuisBANFFMARES)) bestfit_banff_mares = fparm->BestFitValue(conv::kNuisBANFFMARES);
      if(fparm->Exists(conv::kNuisBANFFCCQEE0)) bestfit_banff_ccqeE0 = fparm->BestFitValue(conv::kNuisBANFFCCQEE0);
      if(fparm->Exists(conv::kNuisBANFFCCQEE1)) bestfit_banff_ccqeE1 = fparm->BestFitValue(conv::kNuisBANFFCCQEE1);
      if(fparm->Exists(conv::kNuisBANFFCCQEE2)) bestfit_banff_ccqeE2 = fparm->BestFitValue(conv::kNuisBANFFCCQEE2);
      if(fparm->Exists(conv::kNuisBANFFCC1piE0)) bestfit_banff_cc1piE0 = fparm->BestFitValue(conv::kNuisBANFFCC1piE0);
      if(fparm->Exists(conv::kNuisBANFFCC1piE1)) bestfit_banff_cc1piE1 = fparm->BestFitValue(conv::kNuisBANFFCC1piE1);
      if(fparm->Exists(conv::kNuisBANFFNC1pi0)) bestfit_banff_nc1pi0 = fparm->BestFitValue(conv::kNuisBANFFNC1pi0);
      if(fparm->Exists(conv::kNuisXSecFermiMom)) bestfit_xsec_pf = fparm->BestFitValue(conv::kNuisXSecFermiMom);
      if(fparm->Exists(conv::kNuisXSecSpecFunc)) bestfit_xsec_sf = fparm->BestFitValue(conv::kNuisXSecSpecFunc);
      if(fparm->Exists(conv::kNuisXSecWshape)) bestfit_xsec_wshape = fparm->BestFitValue(conv::kNuisXSecWshape);
      if(fparm->Exists(conv::kNuisXSecShapeCCoth)) bestfit_xsec_ccothshape = fparm->BestFitValue(conv::kNuisXSecShapeCCoth);
      if(fparm->Exists(conv::kNuisXSecPiLessDelta)) bestfit_xsec_cc1pishape = fparm->BestFitValue(conv::kNuisXSecPiLessDelta);
      if(fparm->Exists(conv::kNuisXSecENuShape)) bestfit_xsec_pilessdelta = fparm->BestFitValue(conv::kNuisXSecENuShape);
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
      if(conv::kNumOfFSIParams == 15) {
	if(fparm->Exists(conv::kNuisFSI006)) bestfit_fsi_006 = fparm->BestFitValue(conv::kNuisFSI006);
	if(fparm->Exists(conv::kNuisFSI007)) bestfit_fsi_007 = fparm->BestFitValue(conv::kNuisFSI007);
	if(fparm->Exists(conv::kNuisFSI008)) bestfit_fsi_008 = fparm->BestFitValue(conv::kNuisFSI008);
	if(fparm->Exists(conv::kNuisFSI009)) bestfit_fsi_009 = fparm->BestFitValue(conv::kNuisFSI009);
	if(fparm->Exists(conv::kNuisFSI010)) bestfit_fsi_010 = fparm->BestFitValue(conv::kNuisFSI010);
	if(fparm->Exists(conv::kNuisFSI011)) bestfit_fsi_011 = fparm->BestFitValue(conv::kNuisFSI011);
	if(fparm->Exists(conv::kNuisFSI012)) bestfit_fsi_012 = fparm->BestFitValue(conv::kNuisFSI012);
	if(fparm->Exists(conv::kNuisFSI013)) bestfit_fsi_013 = fparm->BestFitValue(conv::kNuisFSI013);
	if(fparm->Exists(conv::kNuisFSI014)) bestfit_fsi_014 = fparm->BestFitValue(conv::kNuisFSI014);
      }
    }
    //get the plus errors
    if(hyp == kHypOsc3F) {
      errp_mixing12 = fparm->PlusError(kPhypMixing12);
      errp_mixing13 = fparm->PlusError(kPhypMixing13);
      errp_mixing23 = fparm->PlusError(kPhypMixing23);
      errp_deltacp  = fparm->PlusError(kPhypDeltaCP);
      errp_dm12sq   = fparm->PlusError(kPhypDm12sq);
      errp_dm23sq   = fparm->PlusError(kPhypDm23sq);
    } else
      if(hyp == kHypNuMuDisapp2F) {
        errp_mixing23 = fparm->PlusError(kPhypMixing23);
        errp_dm23sq   = fparm->PlusError(kPhypDm23sq);
      }
    if(fit_syst_mode == 1) {
      if(fparm->Exists(conv::kNuisBANFFMAQE)) errp_banff_maqe = fparm->PlusError(conv::kNuisBANFFMAQE);
      if(fparm->Exists(conv::kNuisBANFFCCQEE0)) errp_banff_ccqeE0 = fparm->PlusError(conv::kNuisBANFFCCQEE0);
      if(fparm->Exists(conv::kNuisBANFFCC1piE0)) errp_banff_cc1piE0 = fparm->PlusError(conv::kNuisBANFFCC1piE0);
      if(fparm->Exists(conv::kNuisXSecSpecFunc)) errp_xsec_sf = fparm->PlusError(conv::kNuisXSecSpecFunc);
      if(fparm->Exists(conv::kNuisXSecPiLessDelta)) errp_xsec_cc1pishape = fparm->PlusError(conv::kNuisXSecPiLessDelta);
      if(fparm->Exists(conv::kNuisSuperkEff1RmuCCnonQE)) errp_skdet_ccnqe = fparm->PlusError(conv::kNuisSuperkEff1RmuCCnonQE);
      if(fparm->Exists(conv::kNuisSuperkEff1RmuNC)) errp_skdet_nc = fparm->PlusError(conv::kNuisSuperkEff1RmuNC);
    }
    else if(fit_syst_mode == 2) {
      if(fparm->Exists(conv::kNuisBANFF00)) errp_banff_00 = fparm->PlusError(conv::kNuisBANFF00);
      if(fparm->Exists(conv::kNuisBANFF01)) errp_banff_01 = fparm->PlusError(conv::kNuisBANFF01);
      if(fparm->Exists(conv::kNuisBANFF02)) errp_banff_02 = fparm->PlusError(conv::kNuisBANFF02);
      if(fparm->Exists(conv::kNuisBANFF03)) errp_banff_03 = fparm->PlusError(conv::kNuisBANFF03);
      if(fparm->Exists(conv::kNuisBANFF04)) errp_banff_04 = fparm->PlusError(conv::kNuisBANFF04);
      if(fparm->Exists(conv::kNuisBANFF05)) errp_banff_05 = fparm->PlusError(conv::kNuisBANFF05);
      if(fparm->Exists(conv::kNuisBANFF06)) errp_banff_06 = fparm->PlusError(conv::kNuisBANFF06);
      if(fparm->Exists(conv::kNuisBANFF07)) errp_banff_07 = fparm->PlusError(conv::kNuisBANFF07);
      if(fparm->Exists(conv::kNuisBANFF08)) errp_banff_08 = fparm->PlusError(conv::kNuisBANFF08);
      if(fparm->Exists(conv::kNuisBANFF09)) errp_banff_09 = fparm->PlusError(conv::kNuisBANFF09);
      if(fparm->Exists(conv::kNuisBANFF10)) errp_banff_10 = fparm->PlusError(conv::kNuisBANFF10);
      if(fparm->Exists(conv::kNuisBANFF11)) errp_banff_11 = fparm->PlusError(conv::kNuisBANFF11);
      if(fparm->Exists(conv::kNuisBANFF12)) errp_banff_12 = fparm->PlusError(conv::kNuisBANFF12);
      if(fparm->Exists(conv::kNuisBANFF13)) errp_banff_13 = fparm->PlusError(conv::kNuisBANFF13);
      if(fparm->Exists(conv::kNuisBANFF14)) errp_banff_14 = fparm->PlusError(conv::kNuisBANFF14);
      if(fparm->Exists(conv::kNuisBANFF15)) errp_banff_15 = fparm->PlusError(conv::kNuisBANFF15);
      if(fparm->Exists(conv::kNuisBANFF16)) errp_banff_16 = fparm->PlusError(conv::kNuisBANFF16);
      if(fparm->Exists(conv::kNuisBANFF17)) errp_banff_17 = fparm->PlusError(conv::kNuisBANFF17);
      if(fparm->Exists(conv::kNuisBANFF18)) errp_banff_18 = fparm->PlusError(conv::kNuisBANFF18);
      if(fparm->Exists(conv::kNuisBANFF19)) errp_banff_19 = fparm->PlusError(conv::kNuisBANFF19);
      if(fparm->Exists(conv::kNuisBANFF20)) errp_banff_20 = fparm->PlusError(conv::kNuisBANFF20);
      if(fparm->Exists(conv::kNuisBANFF21)) errp_banff_21 = fparm->PlusError(conv::kNuisBANFF21);
      if(fparm->Exists(conv::kNuisBANFF22)) errp_banff_22 = fparm->PlusError(conv::kNuisBANFF22);
      if(fparm->Exists(conv::kNuisBANFF23)) errp_banff_23 = fparm->PlusError(conv::kNuisBANFF23);
      if(fparm->Exists(conv::kNuisBANFF24)) errp_banff_24 = fparm->PlusError(conv::kNuisBANFF24);
      if(fparm->Exists(conv::kNuisBANFF25)) errp_banff_25 = fparm->PlusError(conv::kNuisBANFF25);
      if(fparm->Exists(conv::kNuisBANFF26)) errp_banff_26 = fparm->PlusError(conv::kNuisBANFF26);
      if(fparm->Exists(conv::kNuisBANFF27)) errp_banff_27 = fparm->PlusError(conv::kNuisBANFF27);
      if(fparm->Exists(conv::kNuisBANFF28)) errp_banff_28 = fparm->PlusError(conv::kNuisBANFF28);
      if(fparm->Exists(conv::kNuisBANFF29)) errp_banff_29 = fparm->PlusError(conv::kNuisBANFF29);
      if(fparm->Exists(conv::kNuisBANFF30)) errp_banff_30 = fparm->PlusError(conv::kNuisBANFF30);
      if(fparm->Exists(conv::kNuisBANFFMAQE)) errp_banff_maqe = fparm->PlusError(conv::kNuisBANFFMAQE);
      if(fparm->Exists(conv::kNuisBANFFMARES)) errp_banff_mares = fparm->PlusError(conv::kNuisBANFFMARES);
      if(fparm->Exists(conv::kNuisBANFFCCQEE0)) errp_banff_ccqeE0 = fparm->PlusError(conv::kNuisBANFFCCQEE0);
      if(fparm->Exists(conv::kNuisBANFFCCQEE1)) errp_banff_ccqeE1 = fparm->PlusError(conv::kNuisBANFFCCQEE1);
      if(fparm->Exists(conv::kNuisBANFFCCQEE2)) errp_banff_ccqeE2 = fparm->PlusError(conv::kNuisBANFFCCQEE2);
      if(fparm->Exists(conv::kNuisBANFFCC1piE0)) errp_banff_cc1piE0 = fparm->PlusError(conv::kNuisBANFFCC1piE0);
      if(fparm->Exists(conv::kNuisBANFFCC1piE1)) errp_banff_cc1piE1 = fparm->PlusError(conv::kNuisBANFFCC1piE1);
      if(fparm->Exists(conv::kNuisBANFFNC1pi0)) errp_banff_nc1pi0 = fparm->PlusError(conv::kNuisBANFFNC1pi0);
      if(fparm->Exists(conv::kNuisXSecFermiMom)) errp_xsec_pf = fparm->PlusError(conv::kNuisXSecFermiMom);
      if(fparm->Exists(conv::kNuisXSecSpecFunc)) errp_xsec_sf = fparm->PlusError(conv::kNuisXSecSpecFunc);
      if(fparm->Exists(conv::kNuisXSecWshape)) errp_xsec_wshape = fparm->PlusError(conv::kNuisXSecWshape);
      if(fparm->Exists(conv::kNuisXSecShapeCCoth)) errp_xsec_ccothshape = fparm->PlusError(conv::kNuisXSecShapeCCoth);
      if(fparm->Exists(conv::kNuisXSecPiLessDelta)) errp_xsec_cc1pishape = fparm->PlusError(conv::kNuisXSecPiLessDelta);
      if(fparm->Exists(conv::kNuisXSecENuShape)) errp_xsec_pilessdelta = fparm->PlusError(conv::kNuisXSecENuShape);
      if(fparm->Exists(conv::kNuisXSecNormCCcoh)) errp_xsec_normcccoh = fparm->PlusError(conv::kNuisXSecNormCCcoh);
      if(fparm->Exists(conv::kNuisXSecNormNC1piC)) errp_xsec_normnc1piC = fparm->PlusError(conv::kNuisXSecNormNC1piC);
      if(fparm->Exists(conv::kNuisXSecNormNCcoh)) errp_xsec_normnccoh = fparm->PlusError(conv::kNuisXSecNormNCcoh);
      if(fparm->Exists(conv::kNuisXSecNormNCoth)) errp_xsec_normncoth = fparm->PlusError(conv::kNuisXSecNormNCoth);
      if(fparm->Exists(conv::kNuisXSecNormNuEToNuMu)) errp_xsec_normnuetonumu = fparm->PlusError(conv::kNuisXSecNormNuEToNuMu);
      if(fparm->Exists(conv::kNuisXSecNormAntiNuToNu)) errp_xsec_normantinutonu = fparm->PlusError(conv::kNuisXSecNormAntiNuToNu);
      if(fparm->Exists(conv::kNuisSuperkEnergyScale)) errp_skdet_energyscale = fparm->PlusError(conv::kNuisSuperkEnergyScale);
      if(fparm->Exists(conv::kNuisSuperkEff1RmuCCQErc0u)) errp_skdet_ccqerc0 = fparm->PlusError(conv::kNuisSuperkEff1RmuCCQErc0u);
      if(fparm->Exists(conv::kNuisSuperkEff1RmuCCQErc1u)) errp_skdet_ccqerc1 = fparm->PlusError(conv::kNuisSuperkEff1RmuCCQErc1u);
      if(fparm->Exists(conv::kNuisSuperkEff1RmuCCQErc2u)) errp_skdet_ccqerc2 = fparm->PlusError(conv::kNuisSuperkEff1RmuCCQErc2u);
      if(fparm->Exists(conv::kNuisSuperkEff1RmuCCQEoth)) errp_skdet_ccqeoth = fparm->PlusError(conv::kNuisSuperkEff1RmuCCQEoth);
      if(fparm->Exists(conv::kNuisSuperkEff1RmuCCnonQE)) errp_skdet_ccnqe = fparm->PlusError(conv::kNuisSuperkEff1RmuCCnonQE);
      if(fparm->Exists(conv::kNuisSuperkEff1RmuNC)) errp_skdet_nc = fparm->PlusError(conv::kNuisSuperkEff1RmuNC);
      if(fparm->Exists(conv::kNuisSuperkEff1RmuNueCC)) errp_skdet_nuecc = fparm->PlusError(conv::kNuisSuperkEff1RmuNueCC);
      if(fparm->Exists(conv::kNuisFSI000)) errp_fsi_000 = fparm->PlusError(conv::kNuisFSI000);
      if(fparm->Exists(conv::kNuisFSI001)) errp_fsi_001 = fparm->PlusError(conv::kNuisFSI001);
      if(fparm->Exists(conv::kNuisFSI002)) errp_fsi_002 = fparm->PlusError(conv::kNuisFSI002);
      if(fparm->Exists(conv::kNuisFSI003)) errp_fsi_003 = fparm->PlusError(conv::kNuisFSI003);
      if(fparm->Exists(conv::kNuisFSI004)) errp_fsi_004 = fparm->PlusError(conv::kNuisFSI004);
      if(fparm->Exists(conv::kNuisFSI005)) errp_fsi_005 = fparm->PlusError(conv::kNuisFSI005);
      if(conv::kNumOfFSIParams == 15) {
	if(fparm->Exists(conv::kNuisFSI006)) errp_fsi_006 = fparm->PlusError(conv::kNuisFSI006);
	if(fparm->Exists(conv::kNuisFSI007)) errp_fsi_007 = fparm->PlusError(conv::kNuisFSI007);
	if(fparm->Exists(conv::kNuisFSI008)) errp_fsi_008 = fparm->PlusError(conv::kNuisFSI008);
	if(fparm->Exists(conv::kNuisFSI009)) errp_fsi_009 = fparm->PlusError(conv::kNuisFSI009);
	if(fparm->Exists(conv::kNuisFSI010)) errp_fsi_010 = fparm->PlusError(conv::kNuisFSI010);
	if(fparm->Exists(conv::kNuisFSI011)) errp_fsi_011 = fparm->PlusError(conv::kNuisFSI011);
	if(fparm->Exists(conv::kNuisFSI012)) errp_fsi_012 = fparm->PlusError(conv::kNuisFSI012);
	if(fparm->Exists(conv::kNuisFSI013)) errp_fsi_013 = fparm->PlusError(conv::kNuisFSI013);
	if(fparm->Exists(conv::kNuisFSI014)) errp_fsi_014 = fparm->PlusError(conv::kNuisFSI014);
      }
    }
    //get the minus errors
    if(hyp == kHypOsc3F) {
      errm_mixing12 = fparm->MinusError(kPhypMixing12);
      errm_mixing13 = fparm->MinusError(kPhypMixing13);
      errm_mixing23 = fparm->MinusError(kPhypMixing23);
      errm_deltacp  = fparm->MinusError(kPhypDeltaCP);
      errm_dm12sq   = fparm->MinusError(kPhypDm12sq);
      errm_dm23sq   = fparm->MinusError(kPhypDm23sq);
    } else
      if(hyp == kHypNuMuDisapp2F) {
        errm_mixing23 = fparm->MinusError(kPhypMixing23);
        errm_dm23sq   = fparm->MinusError(kPhypDm23sq);
      }
    if(fit_syst_mode == 1) {
      if(fparm->Exists(conv::kNuisBANFFMAQE)) errm_banff_maqe = fparm->MinusError(conv::kNuisBANFFMAQE);
      if(fparm->Exists(conv::kNuisBANFFCCQEE0)) errm_banff_ccqeE0 = fparm->MinusError(conv::kNuisBANFFCCQEE0);
      if(fparm->Exists(conv::kNuisBANFFCC1piE0)) errm_banff_cc1piE0 = fparm->MinusError(conv::kNuisBANFFCC1piE0);
      if(fparm->Exists(conv::kNuisXSecSpecFunc)) errm_xsec_sf = fparm->MinusError(conv::kNuisXSecSpecFunc);
      if(fparm->Exists(conv::kNuisXSecPiLessDelta)) errm_xsec_cc1pishape = fparm->MinusError(conv::kNuisXSecPiLessDelta);
      if(fparm->Exists(conv::kNuisSuperkEff1RmuCCnonQE)) errm_skdet_ccnqe = fparm->MinusError(conv::kNuisSuperkEff1RmuCCnonQE);
      if(fparm->Exists(conv::kNuisSuperkEff1RmuNC)) errm_skdet_nc = fparm->MinusError(conv::kNuisSuperkEff1RmuNC);
    }
    else if(fit_syst_mode == 2) {
      if(fparm->Exists(conv::kNuisBANFF00)) errm_banff_00 = fparm->MinusError(conv::kNuisBANFF00);
      if(fparm->Exists(conv::kNuisBANFF01)) errm_banff_01 = fparm->MinusError(conv::kNuisBANFF01);
      if(fparm->Exists(conv::kNuisBANFF02)) errm_banff_02 = fparm->MinusError(conv::kNuisBANFF02);
      if(fparm->Exists(conv::kNuisBANFF03)) errm_banff_03 = fparm->MinusError(conv::kNuisBANFF03);
      if(fparm->Exists(conv::kNuisBANFF04)) errm_banff_04 = fparm->MinusError(conv::kNuisBANFF04);
      if(fparm->Exists(conv::kNuisBANFF05)) errm_banff_05 = fparm->MinusError(conv::kNuisBANFF05);
      if(fparm->Exists(conv::kNuisBANFF06)) errm_banff_06 = fparm->MinusError(conv::kNuisBANFF06);
      if(fparm->Exists(conv::kNuisBANFF07)) errm_banff_07 = fparm->MinusError(conv::kNuisBANFF07);
      if(fparm->Exists(conv::kNuisBANFF08)) errm_banff_08 = fparm->MinusError(conv::kNuisBANFF08);
      if(fparm->Exists(conv::kNuisBANFF09)) errm_banff_09 = fparm->MinusError(conv::kNuisBANFF09);
      if(fparm->Exists(conv::kNuisBANFF10)) errm_banff_10 = fparm->MinusError(conv::kNuisBANFF10);
      if(fparm->Exists(conv::kNuisBANFF11)) errm_banff_11 = fparm->MinusError(conv::kNuisBANFF11);
      if(fparm->Exists(conv::kNuisBANFF12)) errm_banff_12 = fparm->MinusError(conv::kNuisBANFF12);
      if(fparm->Exists(conv::kNuisBANFF13)) errm_banff_13 = fparm->MinusError(conv::kNuisBANFF13);
      if(fparm->Exists(conv::kNuisBANFF14)) errm_banff_14 = fparm->MinusError(conv::kNuisBANFF14);
      if(fparm->Exists(conv::kNuisBANFF15)) errm_banff_15 = fparm->MinusError(conv::kNuisBANFF15);
      if(fparm->Exists(conv::kNuisBANFF16)) errm_banff_16 = fparm->MinusError(conv::kNuisBANFF16);
      if(fparm->Exists(conv::kNuisBANFF17)) errm_banff_17 = fparm->MinusError(conv::kNuisBANFF17);
      if(fparm->Exists(conv::kNuisBANFF18)) errm_banff_18 = fparm->MinusError(conv::kNuisBANFF18);
      if(fparm->Exists(conv::kNuisBANFF19)) errm_banff_19 = fparm->MinusError(conv::kNuisBANFF19);
      if(fparm->Exists(conv::kNuisBANFF20)) errm_banff_20 = fparm->MinusError(conv::kNuisBANFF20);
      if(fparm->Exists(conv::kNuisBANFF21)) errm_banff_21 = fparm->MinusError(conv::kNuisBANFF21);
      if(fparm->Exists(conv::kNuisBANFF22)) errm_banff_22 = fparm->MinusError(conv::kNuisBANFF22);
      if(fparm->Exists(conv::kNuisBANFF23)) errm_banff_23 = fparm->MinusError(conv::kNuisBANFF23);
      if(fparm->Exists(conv::kNuisBANFF24)) errm_banff_24 = fparm->MinusError(conv::kNuisBANFF24);
      if(fparm->Exists(conv::kNuisBANFF25)) errm_banff_25 = fparm->MinusError(conv::kNuisBANFF25);
      if(fparm->Exists(conv::kNuisBANFF26)) errm_banff_26 = fparm->MinusError(conv::kNuisBANFF26);
      if(fparm->Exists(conv::kNuisBANFF27)) errm_banff_27 = fparm->MinusError(conv::kNuisBANFF27);
      if(fparm->Exists(conv::kNuisBANFF28)) errm_banff_28 = fparm->MinusError(conv::kNuisBANFF28);
      if(fparm->Exists(conv::kNuisBANFF29)) errm_banff_29 = fparm->MinusError(conv::kNuisBANFF29);
      if(fparm->Exists(conv::kNuisBANFF30)) errm_banff_30 = fparm->MinusError(conv::kNuisBANFF30);
      if(fparm->Exists(conv::kNuisBANFFMAQE)) errm_banff_maqe = fparm->MinusError(conv::kNuisBANFFMAQE);
      if(fparm->Exists(conv::kNuisBANFFMARES)) errm_banff_mares = fparm->MinusError(conv::kNuisBANFFMARES);
      if(fparm->Exists(conv::kNuisBANFFCCQEE0)) errm_banff_ccqeE0 = fparm->MinusError(conv::kNuisBANFFCCQEE0);
      if(fparm->Exists(conv::kNuisBANFFCCQEE1)) errm_banff_ccqeE1 = fparm->MinusError(conv::kNuisBANFFCCQEE1);
      if(fparm->Exists(conv::kNuisBANFFCCQEE2)) errm_banff_ccqeE2 = fparm->MinusError(conv::kNuisBANFFCCQEE2);
      if(fparm->Exists(conv::kNuisBANFFCC1piE0)) errm_banff_cc1piE0 = fparm->MinusError(conv::kNuisBANFFCC1piE0);
      if(fparm->Exists(conv::kNuisBANFFCC1piE1)) errm_banff_cc1piE1 = fparm->MinusError(conv::kNuisBANFFCC1piE1);
      if(fparm->Exists(conv::kNuisBANFFNC1pi0)) errm_banff_nc1pi0 = fparm->MinusError(conv::kNuisBANFFNC1pi0);
      if(fparm->Exists(conv::kNuisXSecFermiMom)) errm_xsec_pf = fparm->MinusError(conv::kNuisXSecFermiMom);
      if(fparm->Exists(conv::kNuisXSecSpecFunc)) errm_xsec_sf = fparm->MinusError(conv::kNuisXSecSpecFunc);
      if(fparm->Exists(conv::kNuisXSecWshape)) errm_xsec_wshape = fparm->MinusError(conv::kNuisXSecWshape);
      if(fparm->Exists(conv::kNuisXSecShapeCCoth)) errm_xsec_ccothshape = fparm->MinusError(conv::kNuisXSecShapeCCoth);
      if(fparm->Exists(conv::kNuisXSecPiLessDelta)) errm_xsec_cc1pishape = fparm->MinusError(conv::kNuisXSecPiLessDelta);
      if(fparm->Exists(conv::kNuisXSecENuShape)) errm_xsec_pilessdelta = fparm->MinusError(conv::kNuisXSecENuShape);
      if(fparm->Exists(conv::kNuisXSecNormCCcoh)) errm_xsec_normcccoh = fparm->MinusError(conv::kNuisXSecNormCCcoh);
      if(fparm->Exists(conv::kNuisXSecNormNC1piC)) errm_xsec_normnc1piC = fparm->MinusError(conv::kNuisXSecNormNC1piC);
      if(fparm->Exists(conv::kNuisXSecNormNCcoh)) errm_xsec_normnccoh = fparm->MinusError(conv::kNuisXSecNormNCcoh);
      if(fparm->Exists(conv::kNuisXSecNormNCoth)) errm_xsec_normncoth = fparm->MinusError(conv::kNuisXSecNormNCoth);
      if(fparm->Exists(conv::kNuisXSecNormNuEToNuMu)) errm_xsec_normnuetonumu = fparm->MinusError(conv::kNuisXSecNormNuEToNuMu);
      if(fparm->Exists(conv::kNuisXSecNormAntiNuToNu)) errm_xsec_normantinutonu = fparm->MinusError(conv::kNuisXSecNormAntiNuToNu);
      if(fparm->Exists(conv::kNuisSuperkEnergyScale)) errm_skdet_energyscale = fparm->MinusError(conv::kNuisSuperkEnergyScale);
      if(fparm->Exists(conv::kNuisSuperkEff1RmuCCQErc0u)) errm_skdet_ccqerc0 = fparm->MinusError(conv::kNuisSuperkEff1RmuCCQErc0u);
      if(fparm->Exists(conv::kNuisSuperkEff1RmuCCQErc1u)) errm_skdet_ccqerc1 = fparm->MinusError(conv::kNuisSuperkEff1RmuCCQErc1u);
      if(fparm->Exists(conv::kNuisSuperkEff1RmuCCQErc2u)) errm_skdet_ccqerc2 = fparm->MinusError(conv::kNuisSuperkEff1RmuCCQErc2u);
      if(fparm->Exists(conv::kNuisSuperkEff1RmuCCQEoth)) errm_skdet_ccqeoth = fparm->MinusError(conv::kNuisSuperkEff1RmuCCQEoth);
      if(fparm->Exists(conv::kNuisSuperkEff1RmuCCnonQE)) errm_skdet_ccnqe = fparm->MinusError(conv::kNuisSuperkEff1RmuCCnonQE);
      if(fparm->Exists(conv::kNuisSuperkEff1RmuNC)) errm_skdet_nc = fparm->MinusError(conv::kNuisSuperkEff1RmuNC);
      if(fparm->Exists(conv::kNuisSuperkEff1RmuNueCC)) errm_skdet_nuecc = fparm->MinusError(conv::kNuisSuperkEff1RmuNueCC);
      if(fparm->Exists(conv::kNuisFSI000)) errm_fsi_000 = fparm->MinusError(conv::kNuisFSI000);
      if(fparm->Exists(conv::kNuisFSI001)) errm_fsi_001 = fparm->MinusError(conv::kNuisFSI001);
      if(fparm->Exists(conv::kNuisFSI002)) errm_fsi_002 = fparm->MinusError(conv::kNuisFSI002);
      if(fparm->Exists(conv::kNuisFSI003)) errm_fsi_003 = fparm->MinusError(conv::kNuisFSI003);
      if(fparm->Exists(conv::kNuisFSI004)) errm_fsi_004 = fparm->MinusError(conv::kNuisFSI004);
      if(fparm->Exists(conv::kNuisFSI005)) errm_fsi_005 = fparm->MinusError(conv::kNuisFSI005);
      if(conv::kNumOfFSIParams == 15) {
	if(fparm->Exists(conv::kNuisFSI006)) errm_fsi_006 = fparm->MinusError(conv::kNuisFSI006);
	if(fparm->Exists(conv::kNuisFSI007)) errm_fsi_007 = fparm->MinusError(conv::kNuisFSI007);
	if(fparm->Exists(conv::kNuisFSI008)) errm_fsi_008 = fparm->MinusError(conv::kNuisFSI008);
	if(fparm->Exists(conv::kNuisFSI009)) errm_fsi_009 = fparm->MinusError(conv::kNuisFSI009);
	if(fparm->Exists(conv::kNuisFSI010)) errm_fsi_010 = fparm->MinusError(conv::kNuisFSI010);
	if(fparm->Exists(conv::kNuisFSI011)) errm_fsi_011 = fparm->MinusError(conv::kNuisFSI011);
	if(fparm->Exists(conv::kNuisFSI012)) errm_fsi_012 = fparm->MinusError(conv::kNuisFSI012);
	if(fparm->Exists(conv::kNuisFSI013)) errm_fsi_013 = fparm->MinusError(conv::kNuisFSI013);
	if(fparm->Exists(conv::kNuisFSI014)) errm_fsi_014 = fparm->MinusError(conv::kNuisFSI014);
      }
    }
    //get the parabolic errors
    if(hyp == kHypOsc3F) {
      errpar_mixing12 = fparm->ParabolicError(kPhypMixing12);
      errpar_mixing13 = fparm->ParabolicError(kPhypMixing13);
      errpar_mixing23 = fparm->ParabolicError(kPhypMixing23);
      errpar_deltacp  = fparm->ParabolicError(kPhypDeltaCP);
      errpar_dm12sq   = fparm->ParabolicError(kPhypDm12sq);
      errpar_dm23sq   = fparm->ParabolicError(kPhypDm23sq);
    } else
      if(hyp == kHypNuMuDisapp2F) {
        errpar_mixing23 = fparm->ParabolicError(kPhypMixing23);
        errpar_dm23sq   = fparm->ParabolicError(kPhypDm23sq);
      }
    if(fit_syst_mode == 1) {
      if(fparm->Exists(conv::kNuisBANFFMAQE)) errpar_banff_maqe = fparm->ParabolicError(conv::kNuisBANFFMAQE);
      if(fparm->Exists(conv::kNuisBANFFCCQEE0)) errpar_banff_ccqeE0 = fparm->ParabolicError(conv::kNuisBANFFCCQEE0);
      if(fparm->Exists(conv::kNuisBANFFCC1piE0)) errpar_banff_cc1piE0 = fparm->ParabolicError(conv::kNuisBANFFCC1piE0);
      if(fparm->Exists(conv::kNuisXSecSpecFunc)) errpar_xsec_sf = fparm->ParabolicError(conv::kNuisXSecSpecFunc);
      if(fparm->Exists(conv::kNuisXSecPiLessDelta)) errpar_xsec_cc1pishape = fparm->ParabolicError(conv::kNuisXSecPiLessDelta);
      if(fparm->Exists(conv::kNuisSuperkEff1RmuCCnonQE)) errpar_skdet_ccnqe = fparm->ParabolicError(conv::kNuisSuperkEff1RmuCCnonQE);
      if(fparm->Exists(conv::kNuisSuperkEff1RmuNC)) errpar_skdet_nc = fparm->ParabolicError(conv::kNuisSuperkEff1RmuNC);
    }
    else if(fit_syst_mode == 2) {
      if(fparm->Exists(conv::kNuisBANFF00)) errpar_banff_00 = fparm->ParabolicError(conv::kNuisBANFF00);
      if(fparm->Exists(conv::kNuisBANFF01)) errpar_banff_01 = fparm->ParabolicError(conv::kNuisBANFF01);
      if(fparm->Exists(conv::kNuisBANFF02)) errpar_banff_02 = fparm->ParabolicError(conv::kNuisBANFF02);
      if(fparm->Exists(conv::kNuisBANFF03)) errpar_banff_03 = fparm->ParabolicError(conv::kNuisBANFF03);
      if(fparm->Exists(conv::kNuisBANFF04)) errpar_banff_04 = fparm->ParabolicError(conv::kNuisBANFF04);
      if(fparm->Exists(conv::kNuisBANFF05)) errpar_banff_05 = fparm->ParabolicError(conv::kNuisBANFF05);
      if(fparm->Exists(conv::kNuisBANFF06)) errpar_banff_06 = fparm->ParabolicError(conv::kNuisBANFF06);
      if(fparm->Exists(conv::kNuisBANFF07)) errpar_banff_07 = fparm->ParabolicError(conv::kNuisBANFF07);
      if(fparm->Exists(conv::kNuisBANFF08)) errpar_banff_08 = fparm->ParabolicError(conv::kNuisBANFF08);
      if(fparm->Exists(conv::kNuisBANFF09)) errpar_banff_09 = fparm->ParabolicError(conv::kNuisBANFF09);
      if(fparm->Exists(conv::kNuisBANFF10)) errpar_banff_10 = fparm->ParabolicError(conv::kNuisBANFF10);
      if(fparm->Exists(conv::kNuisBANFF11)) errpar_banff_11 = fparm->ParabolicError(conv::kNuisBANFF11);
      if(fparm->Exists(conv::kNuisBANFF12)) errpar_banff_12 = fparm->ParabolicError(conv::kNuisBANFF12);
      if(fparm->Exists(conv::kNuisBANFF13)) errpar_banff_13 = fparm->ParabolicError(conv::kNuisBANFF13);
      if(fparm->Exists(conv::kNuisBANFF14)) errpar_banff_14 = fparm->ParabolicError(conv::kNuisBANFF14);
      if(fparm->Exists(conv::kNuisBANFF15)) errpar_banff_15 = fparm->ParabolicError(conv::kNuisBANFF15);
      if(fparm->Exists(conv::kNuisBANFF16)) errpar_banff_16 = fparm->ParabolicError(conv::kNuisBANFF16);
      if(fparm->Exists(conv::kNuisBANFF17)) errpar_banff_17 = fparm->ParabolicError(conv::kNuisBANFF17);
      if(fparm->Exists(conv::kNuisBANFF18)) errpar_banff_18 = fparm->ParabolicError(conv::kNuisBANFF18);
      if(fparm->Exists(conv::kNuisBANFF19)) errpar_banff_19 = fparm->ParabolicError(conv::kNuisBANFF19);
      if(fparm->Exists(conv::kNuisBANFF20)) errpar_banff_20 = fparm->ParabolicError(conv::kNuisBANFF20);
      if(fparm->Exists(conv::kNuisBANFF21)) errpar_banff_21 = fparm->ParabolicError(conv::kNuisBANFF21);
      if(fparm->Exists(conv::kNuisBANFF22)) errpar_banff_22 = fparm->ParabolicError(conv::kNuisBANFF22);
      if(fparm->Exists(conv::kNuisBANFF23)) errpar_banff_23 = fparm->ParabolicError(conv::kNuisBANFF23);
      if(fparm->Exists(conv::kNuisBANFF24)) errpar_banff_24 = fparm->ParabolicError(conv::kNuisBANFF24);
      if(fparm->Exists(conv::kNuisBANFF25)) errpar_banff_25 = fparm->ParabolicError(conv::kNuisBANFF25);
      if(fparm->Exists(conv::kNuisBANFF26)) errpar_banff_26 = fparm->ParabolicError(conv::kNuisBANFF26);
      if(fparm->Exists(conv::kNuisBANFF27)) errpar_banff_27 = fparm->ParabolicError(conv::kNuisBANFF27);
      if(fparm->Exists(conv::kNuisBANFF28)) errpar_banff_28 = fparm->ParabolicError(conv::kNuisBANFF28);
      if(fparm->Exists(conv::kNuisBANFF29)) errpar_banff_29 = fparm->ParabolicError(conv::kNuisBANFF29);
      if(fparm->Exists(conv::kNuisBANFF30)) errpar_banff_30 = fparm->ParabolicError(conv::kNuisBANFF30);
      if(fparm->Exists(conv::kNuisBANFFMAQE)) errpar_banff_maqe = fparm->ParabolicError(conv::kNuisBANFFMAQE);
      if(fparm->Exists(conv::kNuisBANFFMARES)) errpar_banff_mares = fparm->ParabolicError(conv::kNuisBANFFMARES);
      if(fparm->Exists(conv::kNuisBANFFCCQEE0)) errpar_banff_ccqeE0 = fparm->ParabolicError(conv::kNuisBANFFCCQEE0);
      if(fparm->Exists(conv::kNuisBANFFCCQEE1)) errpar_banff_ccqeE1 = fparm->ParabolicError(conv::kNuisBANFFCCQEE1);
      if(fparm->Exists(conv::kNuisBANFFCCQEE2)) errpar_banff_ccqeE2 = fparm->ParabolicError(conv::kNuisBANFFCCQEE2);
      if(fparm->Exists(conv::kNuisBANFFCC1piE0)) errpar_banff_cc1piE0 = fparm->ParabolicError(conv::kNuisBANFFCC1piE0);
      if(fparm->Exists(conv::kNuisBANFFCC1piE1)) errpar_banff_cc1piE1 = fparm->ParabolicError(conv::kNuisBANFFCC1piE1);
      if(fparm->Exists(conv::kNuisBANFFNC1pi0)) errpar_banff_nc1pi0 = fparm->ParabolicError(conv::kNuisBANFFNC1pi0);
      if(fparm->Exists(conv::kNuisXSecFermiMom)) errpar_xsec_pf = fparm->ParabolicError(conv::kNuisXSecFermiMom);
      if(fparm->Exists(conv::kNuisXSecSpecFunc)) errpar_xsec_sf = fparm->ParabolicError(conv::kNuisXSecSpecFunc);
      if(fparm->Exists(conv::kNuisXSecWshape)) errpar_xsec_wshape = fparm->ParabolicError(conv::kNuisXSecWshape);
      if(fparm->Exists(conv::kNuisXSecShapeCCoth)) errpar_xsec_ccothshape = fparm->ParabolicError(conv::kNuisXSecShapeCCoth);
      if(fparm->Exists(conv::kNuisXSecPiLessDelta)) errpar_xsec_cc1pishape = fparm->ParabolicError(conv::kNuisXSecPiLessDelta);
      if(fparm->Exists(conv::kNuisXSecENuShape)) errpar_xsec_pilessdelta = fparm->ParabolicError(conv::kNuisXSecENuShape);
      if(fparm->Exists(conv::kNuisXSecNormCCcoh)) errpar_xsec_normcccoh = fparm->ParabolicError(conv::kNuisXSecNormCCcoh);
      if(fparm->Exists(conv::kNuisXSecNormNC1piC)) errpar_xsec_normnc1piC = fparm->ParabolicError(conv::kNuisXSecNormNC1piC);
      if(fparm->Exists(conv::kNuisXSecNormNCcoh)) errpar_xsec_normnccoh = fparm->ParabolicError(conv::kNuisXSecNormNCcoh);
      if(fparm->Exists(conv::kNuisXSecNormNCoth)) errpar_xsec_normncoth = fparm->ParabolicError(conv::kNuisXSecNormNCoth);
      if(fparm->Exists(conv::kNuisXSecNormNuEToNuMu)) errpar_xsec_normnuetonumu = fparm->ParabolicError(conv::kNuisXSecNormNuEToNuMu);
      if(fparm->Exists(conv::kNuisXSecNormAntiNuToNu)) errpar_xsec_normantinutonu = fparm->ParabolicError(conv::kNuisXSecNormAntiNuToNu);
      if(fparm->Exists(conv::kNuisSuperkEnergyScale)) errpar_skdet_energyscale = fparm->ParabolicError(conv::kNuisSuperkEnergyScale);
      if(fparm->Exists(conv::kNuisSuperkEff1RmuCCQErc0u)) errpar_skdet_ccqerc0 = fparm->ParabolicError(conv::kNuisSuperkEff1RmuCCQErc0u);
      if(fparm->Exists(conv::kNuisSuperkEff1RmuCCQErc1u)) errpar_skdet_ccqerc1 = fparm->ParabolicError(conv::kNuisSuperkEff1RmuCCQErc1u);
      if(fparm->Exists(conv::kNuisSuperkEff1RmuCCQErc2u)) errpar_skdet_ccqerc2 = fparm->ParabolicError(conv::kNuisSuperkEff1RmuCCQErc2u);
      if(fparm->Exists(conv::kNuisSuperkEff1RmuCCQEoth)) errpar_skdet_ccqeoth = fparm->ParabolicError(conv::kNuisSuperkEff1RmuCCQEoth);
      if(fparm->Exists(conv::kNuisSuperkEff1RmuCCnonQE)) errpar_skdet_ccnqe = fparm->ParabolicError(conv::kNuisSuperkEff1RmuCCnonQE);
      if(fparm->Exists(conv::kNuisSuperkEff1RmuNC)) errpar_skdet_nc = fparm->ParabolicError(conv::kNuisSuperkEff1RmuNC);
      if(fparm->Exists(conv::kNuisSuperkEff1RmuNueCC)) errpar_skdet_nuecc = fparm->ParabolicError(conv::kNuisSuperkEff1RmuNueCC);
      if(fparm->Exists(conv::kNuisFSI000)) errpar_fsi_000 = fparm->ParabolicError(conv::kNuisFSI000);
      if(fparm->Exists(conv::kNuisFSI001)) errpar_fsi_001 = fparm->ParabolicError(conv::kNuisFSI001);
      if(fparm->Exists(conv::kNuisFSI002)) errpar_fsi_002 = fparm->ParabolicError(conv::kNuisFSI002);
      if(fparm->Exists(conv::kNuisFSI003)) errpar_fsi_003 = fparm->ParabolicError(conv::kNuisFSI003);
      if(fparm->Exists(conv::kNuisFSI004)) errpar_fsi_004 = fparm->ParabolicError(conv::kNuisFSI004);
      if(fparm->Exists(conv::kNuisFSI005)) errpar_fsi_005 = fparm->ParabolicError(conv::kNuisFSI005);
      if(conv::kNumOfFSIParams == 15) {
	if(fparm->Exists(conv::kNuisFSI006)) errpar_fsi_006 = fparm->ParabolicError(conv::kNuisFSI006);
	if(fparm->Exists(conv::kNuisFSI007)) errpar_fsi_007 = fparm->ParabolicError(conv::kNuisFSI007);
	if(fparm->Exists(conv::kNuisFSI008)) errpar_fsi_008 = fparm->ParabolicError(conv::kNuisFSI008);
	if(fparm->Exists(conv::kNuisFSI009)) errpar_fsi_009 = fparm->ParabolicError(conv::kNuisFSI009);
	if(fparm->Exists(conv::kNuisFSI010)) errpar_fsi_010 = fparm->ParabolicError(conv::kNuisFSI010);
	if(fparm->Exists(conv::kNuisFSI011)) errpar_fsi_011 = fparm->ParabolicError(conv::kNuisFSI011);
	if(fparm->Exists(conv::kNuisFSI012)) errpar_fsi_012 = fparm->ParabolicError(conv::kNuisFSI012);
	if(fparm->Exists(conv::kNuisFSI013)) errpar_fsi_013 = fparm->ParabolicError(conv::kNuisFSI013);
	if(fparm->Exists(conv::kNuisFSI014)) errpar_fsi_014 = fparm->ParabolicError(conv::kNuisFSI014);
      }
    }

    nfitnt->Fill();

    string output_option = "";
    if ( store_toy_expts          ) { output_option += "Observation";     }
    if ( store_bestfit_spectrum   ) { output_option += "+BestFitSpectra"; }
    if ( store_chisq_surf         ) { output_option += "+ChisqSurface";   }
    fitter->Analyze(output_option.c_str());

    if ( store_withsyst_chisq_surf) {
      double mixing23_min  = 0.8;//fparm.MinValue(conv::kPhypMixing23);
      double mixing23_max  = fparm.MaxValue(conv::kPhypMixing23);
      double mixing23_step = fparm.StepSize(conv::kPhypMixing23);
      double dm23sq_min    = 1.5E-3;//fparm.MinValue(conv::kPhypDm23sq);
      double dm23sq_max    = 3.5E-3;//fparm.MaxValue(conv::kPhypDm23sq);
      double dm23sq_step   = fparm.StepSize(conv::kPhypDm23sq);

      // make sure that min and max are at the bin centre
      int mixing23_nbins = 1 + (int) ((mixing23_max - mixing23_min) / mixing23_step);
      int dm23sq_nbins   = 1 + (int) ((dm23sq_max -   dm23sq_min )  / dm23sq_step );

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

      // Do the fit
      for(int i = 1; i <= chisq_surface->GetNbinsX(); i++) {
	for(int j = 1; j <= chisq_surface->GetNbinsY(); j++) {
	  fparm->ResetCurrToInitValues();
	  double mixing23 = chisq_surface->GetXaxis()->GetBinCenter(i);
	  double dm23sq   = chisq_surface->GetYaxis()->GetBinCenter(j);
	  fparm->Update( conv::kPhypMixing23, mixing23 );
	  fparm->Update( conv::kPhypDm23sq,   dm23sq   );
	  fparm->FixParameter(conv::kPhypMixing23);
	  fparm->FixParameter(conv::kPhypDm23sq);
	  //for some reason, need to make a new fitter when changing the FitParamList
	  FitInputs * fit_inputs2 = new FitInputs(skobs,skexp,fparm);
	  FitterABC * fitter2 = new Superk1RmuErecoFitter(fit_inputs,fit_option);
	  fitter2->DoTheFit();
	  double chisq = fit_inputs2->Chisquare();
	  chisq_surface->SetBinContent(i,j,chisq);
	  delete fitter2;
	  delete fit_inputs2;
	}//j
      }//i
      chisq_surface -> SetName(Form("chisq_withsyst_surface_1Rmu_%d",iexperiment));
      f.cd();
      chisq_surface->Write();
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
      TH2D * chisqsurf = new TH2D(*h2); 
      chisqsurf->SetName(Form("chisqsurf_%d",iexperiment)); 
      f.cd();
      chisqsurf->Write();
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
