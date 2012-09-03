//
// ------------------------------------------------------------------------------------------------
// Run numu disappearance fits to N toy experiments and store best fit values.
// The toy experiments are generated on the fly.
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

/*
TChain csyst("nfitnt");
csyst.Add("syst-run_*");
TCanvas *cansyst = new TCanvas();
cansyst->SetTitle("syst");
cansyst->Divide(2,2);
cansyst->cd(1);
csyst.Draw("bestfit_dm23sq");
cansyst->cd(2);
csyst.Draw("bestfit_dm23sq","fitstatus==3");
TH1F *hsyst_dm23sq = (TH1F*)gPad->GetPrimitive("htemp");
double dm23_mean_syst = hsyst_dm23sq->GetMean();
double dm23_err_syst = hsyst_dm23sq->GetMeanError();
cansyst->cd(3);
csyst.Draw("bestfit_mixing23");
cansyst->cd(4);
csyst.Draw("bestfit_mixing23","fitstatus==3");
TH1F *hsyst_mixing23 = (TH1F*)gPad->GetPrimitive("htemp");
double mixing23_mean_syst = hsyst_mixing23->GetMean();
double mixing23_err_syst = hsyst_mixing23->GetMeanError();

TChain cnosyst("nfitnt");
cnosyst.Add("nosyst-run_*");
TCanvas *cannosyst = new TCanvas();
cannosyst->SetTitle("nosyst");
cannosyst->Divide(2,2);
cannosyst->cd(1);
cnosyst.Draw("bestfit_dm23sq");
cannosyst->cd(2);
cnosyst.Draw("bestfit_dm23sq","fitstatus==3");
TH1F *hnosyst_dm23sq = (TH1F*)gPad->GetPrimitive("htemp");
double dm23_mean_nosyst = hnosyst_dm23sq->GetMean();
double dm23_err_nosyst = hnosyst_dm23sq->GetMeanError();
cannosyst->cd(3);
cnosyst.Draw("bestfit_mixing23");
cannosyst->cd(4);
cnosyst.Draw("bestfit_mixing23","fitstatus==3");
TH1F *hnosyst_mixing23 = (TH1F*)gPad->GetPrimitive("htemp");
double mixing23_mean_nosyst = hnosyst_mixing23->GetMean();
double mixing23_err_nosyst = hnosyst_mixing23->GetMeanError();

//csyst->Scan("*","fitstatus!=3");
//cnosyst->Scan("*","fitstatus!=3");

cout << dm23_err_syst << endl;

TCanvas *pull = new TCanvas();
csyst.Draw("(bestfit_dm23sq-input_dm23sq)/sqrt(pow(9.6E-5,2)-pow(6.72199e-06,2))","fitstatus==3");
TLine l1(htemp->GetMean(),0,htemp->GetMean(),1000);
l1->Draw();
*/


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
  int          run                     = 0,                              /* run number */
  int          nexperiments            = 100,                            /* number of toy MC experiment fits */
  double       input_mixing12          = 0.86,                           /* 12-mixing, sinq2theta or sinsqtheta depending on user choice, if -999 then randomize */
  double       input_mixing13          = 0.10,                           /* 13-mixing, sinq2theta or sinsqtheta depending on user choice, if -999 then randomize */
  double       input_mixing23          = 1.00,                           /* 23-mixing, sinq2theta or sinsqtheta depending on user choice */
  double       input_deltacp           = 0.00,                           /* deltacp, if -999 then randomize */
  double       input_dm12sq            = 7.58E-5,                        /* Dm^2_{12}, if -999 then randomize */
  double       input_dm23sq            = 2.32E-3,                        /* |Dm^2_{23}| */
  const char * config_file             = "analysis_configuration.xml",   /* analysis configuration file */
  bool         stat_fluctuations       = true,                           /* include statistical fluctuation in toy MC expts? */
  bool         tweak_systematics       = true,                           /* generate toy MC with randomized systematics? */
  const char * fit_param_list          = "nuosc3f_numudisapp_fix4+nuisance", /* fit parameter list for the given physics hypothesis (see run_fitter.C) */
  const char * fit_option              = "minuit",                       /* fitter option (see run_fitter.C) */
  bool         store_chisq_surf        = false,                          /* by default storing only best-fit values, enable this to store the chisq surfaces   */
  bool         store_toy_expts         = false,                          /* by default storing only best-fit values, enable this to store the toy MC expts     */
  bool         store_bestfit_spectrum  = false                           /* by default storing only best-fit values, enable this to store the best-fit spectra */
)
{
  using namespace std;
  using namespace t2kana;
  using namespace t2kana::conv;

  // Analysis configuration
  const char * config = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);

  // Global options
  Opt::Instance()->gOscMixingParamChoice = kNuMixSinsq2Theta;
  Opt::Instance()->gMassHierarchyAssumed = kNmhNormal;

  // Check whether some osc. params need to be randomized
  bool randomize_mixing12 = (TMath::Abs(input_mixing12 + 999) < kASmallNum);
  bool randomize_mixing13 = (TMath::Abs(input_mixing13 + 999) < kASmallNum);
  bool randomize_deltacp  = (TMath::Abs(input_deltacp  + 999) < kASmallNum);
  bool randomize_dm12sq   = (TMath::Abs(input_dm12sq   + 999) < kASmallNum);

  // Print settings
  cout << "\n[Job inputs]" << endl;
  cout << "- Job: " << tag << "/" << run << endl;
  cout << "- Reading analysis configuration from: " << config << endl;
  cout << "- Toy-MC generation" << endl;
  cout << "   - Number of toy experiments: " << nexperiments << endl; 
  cout << "   - Randomize systematics : "      << ((tweak_systematics) ?      "YES " : "NO") << endl;
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
  cout << "     - Storing chisq surfaces: "      << ((store_chisq_surf) ?       "YES " : "NO") << endl;
  cout << "     - Storing toy MC experiments:  " << ((store_toy_expts) ?        "YES " : "NO") << endl;
  cout << "     - Storing best-fit spectra: "    << ((store_bestfit_spectrum) ? "YES " : "NO") << endl;
  cout << endl;
  cout << *Opt::Instance() << endl;

  // Fit parameter list
  FitParamList * fparm = FitParamList::NamedParameterList(fit_param_list);
  PhysHypothesis_t hyp = fparm->PhysicsHypothesis();
  cout << "- The following physics hypothesis is used in the fit: " << util::AsString(hyp) << endl;

  // Get superk spectrum predictor & configure
  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read(config);

  // Create file to store results
  string oscstr = util::AsString(hyp,input_mixing12,input_mixing13,input_mixing23,input_deltacp,input_dm12sq,input_dm23sq);
  TFile f(Form("%s-run_%d-%s.root", tag, run, oscstr.c_str()), "recreate");

  // Create ntuple to store oscillation fit results
  TTree * nfitnt = new TTree("nfitnt","nfitnt");
  
  // ***** variables for the branches *******
  //fit info
  int iexp;
  int fitstatus;
  double bestfit_chisq;
  double true_chisq;
  //inputs
  /* already initialised
  double input_mixing12;
  double input_mixing13;
  double input_mixing23;
  double input_deltacp;
  double input_dm12sq;
  double input_dm23sq;
  */
  double input_banff_00;
  double input_banff_01;
  double input_banff_02;
  double input_banff_03;
  double input_banff_04;
  double input_banff_05;
  double input_banff_06;
  double input_banff_07;
  double input_banff_08;
  double input_banff_09;
  double input_banff_10;
  double input_banff_11;
  double input_banff_12;
  double input_banff_13;
  double input_banff_14;
  double input_banff_15;
  double input_banff_16;
  double input_banff_17;
  double input_banff_18;
  double input_banff_19;
  double input_banff_20;
  double input_banff_21;
  double input_banff_22;
  double input_banff_23;
  double input_banff_24;
  double input_banff_25;
  double input_banff_26;
  double input_banff_27;
  double input_banff_28;
  double input_banff_29;
  double input_banff_30;
  double input_banff_maqe;
  double input_banff_mares;
  double input_banff_ccqeE0;
  double input_banff_ccqeE1;
  double input_banff_ccqeE2;
  double input_banff_cc1piE0;
  double input_banff_cc1piE1;
  double input_banff_nc1pi0;
  double input_xsec_pf;
  double input_xsec_sf;
  double input_xsec_wshape;
  double input_xsec_ccothshape;
  double input_xsec_cc1pishape;
  double input_xsec_pilessdelta;
  double input_xsec_normcccoh;
  double input_xsec_normnc1piC;
  double input_xsec_normnccoh;
  double input_xsec_normncoth;
  double input_xsec_normnuetonumu;
  double input_skdet_energyscale;
  double input_skdet_ccqerc0;
  double input_skdet_ccqerc1;
  double input_skdet_ccqerc2;
  double input_skdet_ccqeoth;
  double input_skdet_ccnqe;
  double input_skdet_nc;
  double input_skdet_nuecc;
  double input_fsi_00;
  double input_fsi_01;
  double input_fsi_02;
  double input_fsi_03;
  double input_fsi_04;
  double input_fsi_05;
  double input_fsi_06;
  double input_fsi_07;
  double input_fsi_08;
  double input_fsi_09;
  double input_fsi_10;
  double input_fsi_11;
  double input_fsi_12;
  double input_fsi_13;
  double input_fsi_14;
  //bestfits
  double bestfit_mixing12;
  double bestfit_mixing13;
  double bestfit_mixing23;
  double bestfit_deltacp;
  double bestfit_dm12sq;
  double bestfit_dm23sq;
  double bestfit_banff_00;
  double bestfit_banff_01;
  double bestfit_banff_02;
  double bestfit_banff_03;
  double bestfit_banff_04;
  double bestfit_banff_05;
  double bestfit_banff_06;
  double bestfit_banff_07;
  double bestfit_banff_08;
  double bestfit_banff_09;
  double bestfit_banff_10;
  double bestfit_banff_11;
  double bestfit_banff_12;
  double bestfit_banff_13;
  double bestfit_banff_14;
  double bestfit_banff_15;
  double bestfit_banff_16;
  double bestfit_banff_17;
  double bestfit_banff_18;
  double bestfit_banff_19;
  double bestfit_banff_20;
  double bestfit_banff_21;
  double bestfit_banff_22;
  double bestfit_banff_23;
  double bestfit_banff_24;
  double bestfit_banff_25;
  double bestfit_banff_26;
  double bestfit_banff_27;
  double bestfit_banff_28;
  double bestfit_banff_29;
  double bestfit_banff_30;
  double bestfit_banff_maqe;
  double bestfit_banff_mares;
  double bestfit_banff_ccqeE0;
  double bestfit_banff_ccqeE1;
  double bestfit_banff_ccqeE2;
  double bestfit_banff_cc1piE0;
  double bestfit_banff_cc1piE1;
  double bestfit_banff_nc1pi0;
  double bestfit_xsec_pf;
  double bestfit_xsec_sf;
  double bestfit_xsec_wshape;
  double bestfit_xsec_ccothshape;
  double bestfit_xsec_cc1pishape;
  double bestfit_xsec_pilessdelta;
  double bestfit_xsec_normcccoh;
  double bestfit_xsec_normnc1piC;
  double bestfit_xsec_normnccoh;
  double bestfit_xsec_normncoth;
  double bestfit_xsec_normnuetonumu;
  double bestfit_skdet_energyscale;
  double bestfit_skdet_ccqerc0;
  double bestfit_skdet_ccqerc1;
  double bestfit_skdet_ccqerc2;
  double bestfit_skdet_ccqeoth;
  double bestfit_skdet_ccnqe;
  double bestfit_skdet_nc;
  double bestfit_skdet_nuecc;
  double bestfit_fsi_00;
  double bestfit_fsi_01;
  double bestfit_fsi_02;
  double bestfit_fsi_03;
  double bestfit_fsi_04;
  double bestfit_fsi_05;
  double bestfit_fsi_06;
  double bestfit_fsi_07;
  double bestfit_fsi_08;
  double bestfit_fsi_09;
  double bestfit_fsi_10;
  double bestfit_fsi_11;
  double bestfit_fsi_12;
  double bestfit_fsi_13;
  double bestfit_fsi_14;

  // ***** setup the branches *******
  //fit info
  nfitnt->Branch("iexp",&iexp,"iexp/I");
  nfitnt->Branch("fitstatus",&fitstatus,"fitstatus/I");
  nfitnt->Branch("true_chisq",&true_chisq,"true_chisq/D");
  nfitnt->Branch("bestfit_chisq",&bestfit_chisq,"bestfit_chisq/D");
  //inputs
  nfitnt->Branch("input_mixing12",&input_mixing12,"input_mixing12/D");
  nfitnt->Branch("input_mixing13",&input_mixing13,"input_mixing13/D");
  nfitnt->Branch("input_mixing23",&input_mixing23,"input_mixing23/D");
  nfitnt->Branch("input_deltacp",&input_deltacp,"input_deltacp/D");
  nfitnt->Branch("input_dm12sq",&input_dm12sq,"input_dm12sq/D");
  nfitnt->Branch("input_dm23sq",&input_dm23sq,"input_dm23sq/D");
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
  nfitnt->Branch("input_skdet_energyscale",&input_skdet_energyscale,"input_skdet_energyscale/D");
  nfitnt->Branch("input_skdet_ccqerc0",&input_skdet_ccqerc0,"input_skdet_ccqerc0/D");
  nfitnt->Branch("input_skdet_ccqerc1",&input_skdet_ccqerc1,"input_skdet_ccqerc1/D");
  nfitnt->Branch("input_skdet_ccqerc2",&input_skdet_ccqerc2,"input_skdet_ccqerc2/D");
  nfitnt->Branch("input_skdet_ccqeoth",&input_skdet_ccqeoth,"input_skdet_ccqeoth/D");
  nfitnt->Branch("input_skdet_ccnqe",&input_skdet_ccnqe,"input_skdet_ccnqe/D");
  nfitnt->Branch("input_skdet_nc",&input_skdet_nc,"input_skdet_nc/D");
  nfitnt->Branch("input_skdet_nuecc",&input_skdet_nuecc,"input_skdet_nuecc/D");
  nfitnt->Branch("input_fsi_00",&input_fsi_00,"input_fsi_00/D");
  nfitnt->Branch("input_fsi_01",&input_fsi_01,"input_fsi_01/D");
  nfitnt->Branch("input_fsi_02",&input_fsi_02,"input_fsi_02/D");
  nfitnt->Branch("input_fsi_03",&input_fsi_03,"input_fsi_03/D");
  nfitnt->Branch("input_fsi_04",&input_fsi_04,"input_fsi_04/D");
  nfitnt->Branch("input_fsi_05",&input_fsi_05,"input_fsi_05/D");
  nfitnt->Branch("input_fsi_06",&input_fsi_06,"input_fsi_06/D");
  nfitnt->Branch("input_fsi_07",&input_fsi_07,"input_fsi_07/D");
  nfitnt->Branch("input_fsi_08",&input_fsi_08,"input_fsi_08/D");
  nfitnt->Branch("input_fsi_09",&input_fsi_09,"input_fsi_09/D");
  nfitnt->Branch("input_fsi_10",&input_fsi_10,"input_fsi_10/D");
  nfitnt->Branch("input_fsi_11",&input_fsi_11,"input_fsi_11/D");
  nfitnt->Branch("input_fsi_12",&input_fsi_12,"input_fsi_12/D");
  nfitnt->Branch("input_fsi_13",&input_fsi_13,"input_fsi_13/D");
  nfitnt->Branch("input_fsi_14",&input_fsi_14,"input_fsi_14/D");
  //bestfits
  nfitnt->Branch("bestfit_mixing12",&bestfit_mixing12,"bestfit_mixing12/D");
  nfitnt->Branch("bestfit_mixing13",&bestfit_mixing13,"bestfit_mixing13/D");
  nfitnt->Branch("bestfit_mixing23",&bestfit_mixing23,"bestfit_mixing23/D");
  nfitnt->Branch("bestfit_deltacp",&bestfit_deltacp,"bestfit_deltacp/D");
  nfitnt->Branch("bestfit_dm12sq",&bestfit_dm12sq,"bestfit_dm12sq/D");
  nfitnt->Branch("bestfit_dm23sq",&bestfit_dm23sq,"bestfit_dm23sq/D");
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
  nfitnt->Branch("bestfit_skdet_energyscale",&bestfit_skdet_energyscale,"bestfit_skdet_energyscale/D");
  nfitnt->Branch("bestfit_skdet_ccqerc0",&bestfit_skdet_ccqerc0,"bestfit_skdet_ccqerc0/D");
  nfitnt->Branch("bestfit_skdet_ccqerc1",&bestfit_skdet_ccqerc1,"bestfit_skdet_ccqerc1/D");
  nfitnt->Branch("bestfit_skdet_ccqerc2",&bestfit_skdet_ccqerc2,"bestfit_skdet_ccqerc2/D");
  nfitnt->Branch("bestfit_skdet_ccqeoth",&bestfit_skdet_ccqeoth,"bestfit_skdet_ccqeoth/D");
  nfitnt->Branch("bestfit_skdet_ccnqe",&bestfit_skdet_ccnqe,"bestfit_skdet_ccnqe/D");
  nfitnt->Branch("bestfit_skdet_nc",&bestfit_skdet_nc,"bestfit_skdet_nc/D");
  nfitnt->Branch("bestfit_skdet_nuecc",&bestfit_skdet_nuecc,"bestfit_skdet_nuecc/D");
  nfitnt->Branch("bestfit_fsi_00",&bestfit_fsi_00,"bestfit_fsi_00/D");
  nfitnt->Branch("bestfit_fsi_01",&bestfit_fsi_01,"bestfit_fsi_01/D");
  nfitnt->Branch("bestfit_fsi_02",&bestfit_fsi_02,"bestfit_fsi_02/D");
  nfitnt->Branch("bestfit_fsi_03",&bestfit_fsi_03,"bestfit_fsi_03/D");
  nfitnt->Branch("bestfit_fsi_04",&bestfit_fsi_04,"bestfit_fsi_04/D");
  nfitnt->Branch("bestfit_fsi_05",&bestfit_fsi_05,"bestfit_fsi_05/D");
  nfitnt->Branch("bestfit_fsi_06",&bestfit_fsi_06,"bestfit_fsi_06/D");
  nfitnt->Branch("bestfit_fsi_07",&bestfit_fsi_07,"bestfit_fsi_07/D");
  nfitnt->Branch("bestfit_fsi_08",&bestfit_fsi_08,"bestfit_fsi_08/D");
  nfitnt->Branch("bestfit_fsi_09",&bestfit_fsi_09,"bestfit_fsi_09/D");
  nfitnt->Branch("bestfit_fsi_10",&bestfit_fsi_10,"bestfit_fsi_10/D");
  nfitnt->Branch("bestfit_fsi_11",&bestfit_fsi_11,"bestfit_fsi_11/D");
  nfitnt->Branch("bestfit_fsi_12",&bestfit_fsi_12,"bestfit_fsi_12/D");
  nfitnt->Branch("bestfit_fsi_13",&bestfit_fsi_13,"bestfit_fsi_13/D");
  nfitnt->Branch("bestfit_fsi_14",&bestfit_fsi_14,"bestfit_fsi_14/D");

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
     if(tweak_systematics) {
        // tweak all sources of systematic err (gaussian err)
        cout << "- Randomizing all systematic params in pdf..." << endl;
        skexp->TweakAllSystematics();     
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
        input_skdet_energyscale = skexp->GetSuperkDetector()->CurrentValue(conv::kNuisSuperkEnergyScale);
        input_skdet_ccqerc0 = skexp->GetSuperkDetector()->CurrentValue(conv::kNuisSuperkEff1RmuCCQErc0u);
        input_skdet_ccqerc1 = skexp->GetSuperkDetector()->CurrentValue(conv::kNuisSuperkEff1RmuCCQErc1u);
        input_skdet_ccqerc2 = skexp->GetSuperkDetector()->CurrentValue(conv::kNuisSuperkEff1RmuCCQErc2u);
        input_skdet_ccqeoth = skexp->GetSuperkDetector()->CurrentValue(conv::kNuisSuperkEff1RmuCCQEoth);
        input_skdet_ccnqe = skexp->GetSuperkDetector()->CurrentValue(conv::kNuisSuperkEff1RmuCCnonQE);
        input_skdet_nc = skexp->GetSuperkDetector()->CurrentValue(conv::kNuisSuperkEff1RmuNC);
        input_skdet_nuecc = skexp->GetSuperkDetector()->CurrentValue(conv::kNuisSuperkEff1RmuNueCC);
	input_fsi_00 = skexp->GetFSI()->CurrentValue(conv::kNuisFSI00);
	input_fsi_01 = skexp->GetFSI()->CurrentValue(conv::kNuisFSI01);
	input_fsi_02 = skexp->GetFSI()->CurrentValue(conv::kNuisFSI02);
	input_fsi_03 = skexp->GetFSI()->CurrentValue(conv::kNuisFSI03);
	input_fsi_04 = skexp->GetFSI()->CurrentValue(conv::kNuisFSI04);
	input_fsi_05 = skexp->GetFSI()->CurrentValue(conv::kNuisFSI05);
	input_fsi_06 = skexp->GetFSI()->CurrentValue(conv::kNuisFSI06);
	input_fsi_07 = skexp->GetFSI()->CurrentValue(conv::kNuisFSI07);
	input_fsi_08 = skexp->GetFSI()->CurrentValue(conv::kNuisFSI08);
	input_fsi_09 = skexp->GetFSI()->CurrentValue(conv::kNuisFSI09);
	input_fsi_10 = skexp->GetFSI()->CurrentValue(conv::kNuisFSI10);
	input_fsi_11 = skexp->GetFSI()->CurrentValue(conv::kNuisFSI11);
	input_fsi_12 = skexp->GetFSI()->CurrentValue(conv::kNuisFSI12);
	input_fsi_13 = skexp->GetFSI()->CurrentValue(conv::kNuisFSI13);
	input_fsi_14 = skexp->GetFSI()->CurrentValue(conv::kNuisFSI14);
	//update the FitParamList with the fluctuated nuisance parameters
	fparm->Update(conv::kNuisBANFF00, input_banff_00);
	fparm->Update(conv::kNuisBANFF01, input_banff_01);
	fparm->Update(conv::kNuisBANFF02, input_banff_02);
	fparm->Update(conv::kNuisBANFF03, input_banff_03);
	fparm->Update(conv::kNuisBANFF04, input_banff_04);
	fparm->Update(conv::kNuisBANFF05, input_banff_05);
	fparm->Update(conv::kNuisBANFF06, input_banff_06);
	fparm->Update(conv::kNuisBANFF07, input_banff_07);
	fparm->Update(conv::kNuisBANFF08, input_banff_08);
	fparm->Update(conv::kNuisBANFF09, input_banff_09);
	fparm->Update(conv::kNuisBANFF10, input_banff_10);
	fparm->Update(conv::kNuisBANFF11, input_banff_11);
	fparm->Update(conv::kNuisBANFF12, input_banff_12);
	fparm->Update(conv::kNuisBANFF13, input_banff_13);
	fparm->Update(conv::kNuisBANFF14, input_banff_14);
	fparm->Update(conv::kNuisBANFF15, input_banff_15);
	fparm->Update(conv::kNuisBANFF16, input_banff_16);
	fparm->Update(conv::kNuisBANFF17, input_banff_17);
	fparm->Update(conv::kNuisBANFF18, input_banff_18);
	fparm->Update(conv::kNuisBANFF19, input_banff_19);
	fparm->Update(conv::kNuisBANFF20, input_banff_20);
	fparm->Update(conv::kNuisBANFF21, input_banff_21);
	fparm->Update(conv::kNuisBANFF22, input_banff_22);
	fparm->Update(conv::kNuisBANFF23, input_banff_23);
	fparm->Update(conv::kNuisBANFF24, input_banff_24);
	fparm->Update(conv::kNuisBANFF25, input_banff_25);
	fparm->Update(conv::kNuisBANFF26, input_banff_26);
	fparm->Update(conv::kNuisBANFF27, input_banff_27);
	fparm->Update(conv::kNuisBANFF28, input_banff_28);
	fparm->Update(conv::kNuisBANFF29, input_banff_29);
	fparm->Update(conv::kNuisBANFF30, input_banff_30);
	fparm->Update(conv::kNuisBANFFMAQE, input_banff_maqe);
	fparm->Update(conv::kNuisBANFFMARES, input_banff_mares);
	fparm->Update(conv::kNuisBANFFCCQEE0, input_banff_ccqeE0);
	fparm->Update(conv::kNuisBANFFCCQEE1, input_banff_ccqeE1);
	fparm->Update(conv::kNuisBANFFCCQEE2, input_banff_ccqeE2);
	fparm->Update(conv::kNuisBANFFCC1piE0, input_banff_cc1piE0);
	fparm->Update(conv::kNuisBANFFCC1piE1, input_banff_cc1piE1);
	fparm->Update(conv::kNuisBANFFNC1pi0, input_banff_nc1pi0);
	fparm->Update(conv::kNuisXSecFermiMom, input_xsec_pf);
	fparm->Update(conv::kNuisXSecSpecFunc, input_xsec_sf);
	fparm->Update(conv::kNuisXSecWshape, input_xsec_wshape);
	fparm->Update(conv::kNuisXSecShapeCCoth, input_xsec_ccothshape);
	fparm->Update(conv::kNuisXSecPiLessDelta, input_xsec_cc1pishape);
	fparm->Update(conv::kNuisXSecENuShape, input_xsec_pilessdelta);
	fparm->Update(conv::kNuisXSecNormCCcoh, input_xsec_normcccoh);
	fparm->Update(conv::kNuisXSecNormNC1piC, input_xsec_normnc1piC);
	fparm->Update(conv::kNuisXSecNormNCcoh, input_xsec_normnccoh);
	fparm->Update(conv::kNuisXSecNormNCoth, input_xsec_normncoth);
	cout << "ha" << endl;
	fparm->Update(conv::kNuisXSecNormNuEToNuMu, input_xsec_normnuetonumu);
        fparm->Update(conv::kNuisSuperkEnergyScale, input_skdet_energyscale);
        fparm->Update(conv::kNuisSuperkEff1RmuCCQErc0u, input_skdet_ccqerc0);
        fparm->Update(conv::kNuisSuperkEff1RmuCCQErc1u, input_skdet_ccqerc1);
        fparm->Update(conv::kNuisSuperkEff1RmuCCQErc2u, input_skdet_ccqerc2);
        fparm->Update(conv::kNuisSuperkEff1RmuCCQEoth, input_skdet_ccqeoth);
        fparm->Update(conv::kNuisSuperkEff1RmuCCnonQE, input_skdet_ccnqe);
        fparm->Update(conv::kNuisSuperkEff1RmuNC, input_skdet_nc);
        fparm->Update(conv::kNuisSuperkEff1RmuNueCC, input_skdet_nuecc);
	fparm->Update(conv::kNuisFSI00, input_fsi_00);
	fparm->Update(conv::kNuisFSI01, input_fsi_01);
	fparm->Update(conv::kNuisFSI02, input_fsi_02);
	fparm->Update(conv::kNuisFSI03, input_fsi_03);
	fparm->Update(conv::kNuisFSI04, input_fsi_04);
	fparm->Update(conv::kNuisFSI05, input_fsi_05);
	fparm->Update(conv::kNuisFSI06, input_fsi_06);
	fparm->Update(conv::kNuisFSI07, input_fsi_07);
	fparm->Update(conv::kNuisFSI08, input_fsi_08);
	fparm->Update(conv::kNuisFSI09, input_fsi_09);
	fparm->Update(conv::kNuisFSI10, input_fsi_10);
	fparm->Update(conv::kNuisFSI11, input_fsi_11);
	fparm->Update(conv::kNuisFSI12, input_fsi_12);
	fparm->Update(conv::kNuisFSI13, input_fsi_13);
	fparm->Update(conv::kNuisFSI14, input_fsi_14);
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
     bestfit_skdet_energyscale = -1;
     bestfit_skdet_ccqerc0 = -1;
     bestfit_skdet_ccqerc1 = -1;
     bestfit_skdet_ccqerc2 = -1;
     bestfit_skdet_ccqeoth = -1;
     bestfit_skdet_ccnqe = -1;
     bestfit_skdet_nc = -1;
     bestfit_skdet_nuecc = -1;
     bestfit_fsi_00 = -1;
     bestfit_fsi_01 = -1;
     bestfit_fsi_02 = -1;
     bestfit_fsi_03 = -1;
     bestfit_fsi_04 = -1;
     bestfit_fsi_05 = -1;
     bestfit_fsi_06 = -1;
     bestfit_fsi_07 = -1;
     bestfit_fsi_08 = -1;
     bestfit_fsi_09 = -1;
     bestfit_fsi_10 = -1;
     bestfit_fsi_11 = -1;
     bestfit_fsi_12 = -1;
     bestfit_fsi_13 = -1;
     bestfit_fsi_14 = -1;
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
     bestfit_banff_00 = fparm->BestFitValue(conv::kNuisBANFF00);
     bestfit_banff_01 = fparm->BestFitValue(conv::kNuisBANFF01);
     bestfit_banff_02 = fparm->BestFitValue(conv::kNuisBANFF02);
     bestfit_banff_03 = fparm->BestFitValue(conv::kNuisBANFF03);
     bestfit_banff_04 = fparm->BestFitValue(conv::kNuisBANFF04);
     bestfit_banff_05 = fparm->BestFitValue(conv::kNuisBANFF05);
     bestfit_banff_06 = fparm->BestFitValue(conv::kNuisBANFF06);
     bestfit_banff_07 = fparm->BestFitValue(conv::kNuisBANFF07);
     bestfit_banff_08 = fparm->BestFitValue(conv::kNuisBANFF08);
     bestfit_banff_09 = fparm->BestFitValue(conv::kNuisBANFF09);
     bestfit_banff_10 = fparm->BestFitValue(conv::kNuisBANFF10);
     bestfit_banff_11 = fparm->BestFitValue(conv::kNuisBANFF11);
     bestfit_banff_12 = fparm->BestFitValue(conv::kNuisBANFF12);
     bestfit_banff_13 = fparm->BestFitValue(conv::kNuisBANFF13);
     bestfit_banff_14 = fparm->BestFitValue(conv::kNuisBANFF14);
     bestfit_banff_15 = fparm->BestFitValue(conv::kNuisBANFF15);
     bestfit_banff_16 = fparm->BestFitValue(conv::kNuisBANFF16);
     bestfit_banff_17 = fparm->BestFitValue(conv::kNuisBANFF17);
     bestfit_banff_18 = fparm->BestFitValue(conv::kNuisBANFF18);
     bestfit_banff_19 = fparm->BestFitValue(conv::kNuisBANFF19);
     bestfit_banff_20 = fparm->BestFitValue(conv::kNuisBANFF20);
     bestfit_banff_21 = fparm->BestFitValue(conv::kNuisBANFF21);
     bestfit_banff_22 = fparm->BestFitValue(conv::kNuisBANFF22);
     bestfit_banff_23 = fparm->BestFitValue(conv::kNuisBANFF23);
     bestfit_banff_24 = fparm->BestFitValue(conv::kNuisBANFF24);
     bestfit_banff_25 = fparm->BestFitValue(conv::kNuisBANFF25);
     bestfit_banff_26 = fparm->BestFitValue(conv::kNuisBANFF26);
     bestfit_banff_27 = fparm->BestFitValue(conv::kNuisBANFF27);
     bestfit_banff_28 = fparm->BestFitValue(conv::kNuisBANFF28);
     bestfit_banff_29 = fparm->BestFitValue(conv::kNuisBANFF29);
     bestfit_banff_30 = fparm->BestFitValue(conv::kNuisBANFF30);
     bestfit_banff_maqe = fparm->BestFitValue(conv::kNuisBANFFMAQE);
     bestfit_banff_mares = fparm->BestFitValue(conv::kNuisBANFFMARES);
     bestfit_banff_ccqeE0 = fparm->BestFitValue(conv::kNuisBANFFCCQEE0);
     bestfit_banff_ccqeE1 = fparm->BestFitValue(conv::kNuisBANFFCCQEE1);
     bestfit_banff_ccqeE2 = fparm->BestFitValue(conv::kNuisBANFFCCQEE2);
     bestfit_banff_cc1piE0 = fparm->BestFitValue(conv::kNuisBANFFCC1piE0);
     bestfit_banff_cc1piE1 = fparm->BestFitValue(conv::kNuisBANFFCC1piE1);
     bestfit_banff_nc1pi0 = fparm->BestFitValue(conv::kNuisBANFFNC1pi0);
     bestfit_xsec_pf = fparm->BestFitValue(conv::kNuisXSecFermiMom);
     bestfit_xsec_sf = fparm->BestFitValue(conv::kNuisXSecSpecFunc);
     bestfit_xsec_wshape = fparm->BestFitValue(conv::kNuisXSecWshape);
     bestfit_xsec_ccothshape = fparm->BestFitValue(conv::kNuisXSecShapeCCoth);
     bestfit_xsec_cc1pishape = fparm->BestFitValue(conv::kNuisXSecPiLessDelta);
     bestfit_xsec_pilessdelta = fparm->BestFitValue(conv::kNuisXSecENuShape);
     bestfit_xsec_normcccoh = fparm->BestFitValue(conv::kNuisXSecNormCCcoh);
     bestfit_xsec_normnc1piC = fparm->BestFitValue(conv::kNuisXSecNormNC1piC);
     bestfit_xsec_normnccoh = fparm->BestFitValue(conv::kNuisXSecNormNCcoh);
     bestfit_xsec_normncoth = fparm->BestFitValue(conv::kNuisXSecNormNCoth);
     bestfit_xsec_normnuetonumu = fparm->BestFitValue(conv::kNuisXSecNormNuEToNuMu);
     bestfit_skdet_energyscale = fparm->BestFitValue(conv::kNuisSuperkEnergyScale);
     bestfit_skdet_ccqerc0 = fparm->BestFitValue(conv::kNuisSuperkEff1RmuCCQErc0u);
     bestfit_skdet_ccqerc1 = fparm->BestFitValue(conv::kNuisSuperkEff1RmuCCQErc1u);
     bestfit_skdet_ccqerc2 = fparm->BestFitValue(conv::kNuisSuperkEff1RmuCCQErc2u);
     bestfit_skdet_ccqeoth = fparm->BestFitValue(conv::kNuisSuperkEff1RmuCCQEoth);
     bestfit_skdet_ccnqe = fparm->BestFitValue(conv::kNuisSuperkEff1RmuCCnonQE);
     bestfit_skdet_nc = fparm->BestFitValue(conv::kNuisSuperkEff1RmuNC);
     bestfit_skdet_nuecc = fparm->BestFitValue(conv::kNuisSuperkEff1RmuNueCC);
     bestfit_fsi_00 = fparm->BestFitValue(conv::kNuisFSI00);
     bestfit_fsi_01 = fparm->BestFitValue(conv::kNuisFSI01);
     bestfit_fsi_02 = fparm->BestFitValue(conv::kNuisFSI02);
     bestfit_fsi_03 = fparm->BestFitValue(conv::kNuisFSI03);
     bestfit_fsi_04 = fparm->BestFitValue(conv::kNuisFSI04);
     bestfit_fsi_05 = fparm->BestFitValue(conv::kNuisFSI05);
     bestfit_fsi_06 = fparm->BestFitValue(conv::kNuisFSI06);
     bestfit_fsi_07 = fparm->BestFitValue(conv::kNuisFSI07);
     bestfit_fsi_08 = fparm->BestFitValue(conv::kNuisFSI08);
     bestfit_fsi_09 = fparm->BestFitValue(conv::kNuisFSI09);
     bestfit_fsi_10 = fparm->BestFitValue(conv::kNuisFSI10);
     bestfit_fsi_11 = fparm->BestFitValue(conv::kNuisFSI11);
     bestfit_fsi_12 = fparm->BestFitValue(conv::kNuisFSI12);
     bestfit_fsi_13 = fparm->BestFitValue(conv::kNuisFSI13);
     bestfit_fsi_14 = fparm->BestFitValue(conv::kNuisFSI14);

     nfitnt->Fill();

     string output_option = "";
     if ( store_toy_expts        ) { output_option += "Observation";     }
     if ( store_bestfit_spectrum ) { output_option += "+BestFitSpectra"; }
     if ( store_chisq_surf       ) { output_option += "+ChisqSurface";   }
     fitter->Analyze(output_option.c_str());

     // Store toy-MC expt (if asked)
     if(store_toy_expts) {
        TH1D * h = (TH1D*) fitter->GetResults().FindObject("data");
        TH1D * toyexpspectrum = new TH1D(*h);
        f.cd();
        toyexpspectrum->SetName(Form("toyexp_ereco_spectrum_1Rmu_%d",iexperiment));
        toyexpspectrum->Write();
     }

     // Store best-fit spectrum (if asked)
     if(store_bestfit_spectrum) {
        TH1D * h = (TH1D*) fitter->GetResults().FindObject("bestfit_tot");
        TH1D * bestfitspectrum = new TH1D(*h);
        bestfitspectrum->SetName(Form("bestfit_ereco_spectrum_1Rmu_%d",iexperiment));
        f.cd();
        bestfitspectrum->Write();
     }
     // Store chisq surface (if asked)
     if(store_chisq_surf) {
        TH2D * h2 = (TH2D*) fitter->GetResults().FindObject("chisq_surface");
        TH2D * chisqsurf = new TH2D(*h2); 
	chisqsurf->SetName(Form("chisqsurf_%d",iexperiment)); 
	f.cd();
	chisqsurf->Write();
     }

     delete skobs;
     delete fit_inputs;
     delete fitter;
  }

  // Store ntuple with best-fit results
  f.cd();
  nfitnt->Write("nfitnt");

  f.Close();
}
