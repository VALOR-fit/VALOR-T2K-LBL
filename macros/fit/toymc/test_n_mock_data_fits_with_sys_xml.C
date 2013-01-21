//
// ------------------------------------------------------------------------------------------------
// Run numu disappearance fits to N toy experiments and store best fit values.
// The toy experiments are read in from an xml file
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

void test_n_mock_data_fits_with_sys_xml
(
 const char * xmlfile,                                                  /* input XML file with toy MC experiments */
 const char * tag                     = "toyexptfitswithsyst",          /* tag for output file */
 int          run                     = 0,                              /* run number */
 int          first_experiment        = 0,                              /* */
 int          last_experiment         = 0,                              /* */
 double       input_mixing12          = -99,                            /* 12-mixing, sinq2theta or sinsqtheta depending on user choice, if -999 then randomize */
 double       input_mixing13          = -99,                            /* 13-mixing, sinq2theta or sinsqtheta depending on user choice, if -999 then randomize */
 double       input_mixing23          = -99,                            /* 23-mixing, sinq2theta or sinsqtheta depending on user choice */
 double       input_deltacp           = -99,                            /* deltacp, if -999 then randomize */
 double       input_dm12sq            = -99,                            /* Dm^2_{12}, if -999 then randomize */
 double       input_dm23sq            = -99,                            /* |Dm^2_{23}| */
 const char * config_file             = "analysis_configuration.xml",   /* analysis configuration file */
 string       fit_param_list          = "nuosc3f_numudisapp_fix4+nuisanceselection", /* fit parameter list for the given physics hypothesis (see run_fitter.C) */
 const char * fit_option              = "minuit"                        /* fitter option (see run_fitter.C) */
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

  // Print settings
  cout << "\n[Job inputs]" << endl;
  cout << "- Job: " << tag << "/" << run << endl;
  cout << "- Toy-MC generation" << endl;
  cout << "   - True oscillation parameters: " << endl;
  cout << "       > " << util::AsString(kPhypMixing12) << " = " << input_mixing12 << endl;
  cout << "       > " << util::AsString(kPhypMixing13) << " = " << input_mixing13 << endl;
  cout << "       > " << util::AsString(kPhypMixing23) << " = " << input_mixing23 << endl;
  cout << "       > " << util::AsString(kPhypDeltaCP)  << " = " << input_deltacp  << endl;
  cout << "       > " << util::AsString(kPhypDm12sq)   << " = " << input_dm12sq   << " eV^2/c^4" << endl;
  cout << "       > " << util::AsString(kPhypDm23sq)   << " = " << input_dm23sq   << " eV^2/c^4" << endl;
  cout << "- Fit" << endl;
  cout << "   - Fit Parameter List: " << fit_param_list << endl;
  cout << "   - Fitter option: " << fit_option << endl;
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
  TFile f(Form("%s-run_%d-fitparam_%s-%s.root", tag, run, fit_param_list.c_str(), oscstr.c_str()), "recreate");

  // Create ntuple to store oscillation fit results
  TTree * nfitnt = new TTree("nfitnt","nfitnt");
  
  ///////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////
  // ***** variables for the branches *******
  ///////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////
  //fit info
  int iexp;
  int fitstatus, migradfitstatus, hessefitstatus;
  double fittime;
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

  ///////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////
  // ***** setup the branches *******
  ///////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////
  //fit info
  nfitnt->Branch("iexp",&iexp,"iexp/I");
  nfitnt->Branch("fitstatus",&fitstatus,"fitstatus/I");
  nfitnt->Branch("migradfitstatus",&migradfitstatus,"migradfitstatus/I");
  nfitnt->Branch("hessefitstatus",&hessefitstatus,"hessefitstatus/I");
  nfitnt->Branch("fittime",&fittime,"fittime/D");
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
  //bestfits
  nfitnt->Branch("bestfit_mixing12",&bestfit_mixing12,"bestfit_mixing12/D");
  nfitnt->Branch("bestfit_mixing13",&bestfit_mixing13,"bestfit_mixing13/D");
  nfitnt->Branch("bestfit_mixing23",&bestfit_mixing23,"bestfit_mixing23/D");
  nfitnt->Branch("bestfit_deltacp",&bestfit_deltacp,"bestfit_deltacp/D");
  nfitnt->Branch("bestfit_dm12sq",&bestfit_dm12sq,"bestfit_dm12sq/D");
  nfitnt->Branch("bestfit_dm23sq",&bestfit_dm23sq,"bestfit_dm23sq/D");
  //errminospluss
  nfitnt->Branch("errminosplus_mixing12",&errminosplus_mixing12,"errminosplus_mixing12/D");
  nfitnt->Branch("errminosplus_mixing13",&errminosplus_mixing13,"errminosplus_mixing13/D");
  nfitnt->Branch("errminosplus_mixing23",&errminosplus_mixing23,"errminosplus_mixing23/D");
  nfitnt->Branch("errminosplus_deltacp",&errminosplus_deltacp,"errminosplus_deltacp/D");
  nfitnt->Branch("errminosplus_dm12sq",&errminosplus_dm12sq,"errminosplus_dm12sq/D");
  nfitnt->Branch("errminosplus_dm23sq",&errminosplus_dm23sq,"errminosplus_dm23sq/D");
  //errms
  nfitnt->Branch("errminosminus_mixing12",&errminosminus_mixing12,"errminosminus_mixing12/D");
  nfitnt->Branch("errminosminus_mixing13",&errminosminus_mixing13,"errminosminus_mixing13/D");
  nfitnt->Branch("errminosminus_mixing23",&errminosminus_mixing23,"errminosminus_mixing23/D");
  nfitnt->Branch("errminosminus_deltacp",&errminosminus_deltacp,"errminosminus_deltacp/D");
  nfitnt->Branch("errminosminus_dm12sq",&errminosminus_dm12sq,"errminosminus_dm12sq/D");
  nfitnt->Branch("errminosminus_dm23sq",&errminosminus_dm23sq,"errminosminus_dm23sq/D");
  //errmigrads
  nfitnt->Branch("errmigrad_mixing12",&errmigrad_mixing12,"errmigrad_mixing12/D");
  nfitnt->Branch("errmigrad_mixing13",&errmigrad_mixing13,"errmigrad_mixing13/D");
  nfitnt->Branch("errmigrad_mixing23",&errmigrad_mixing23,"errmigrad_mixing23/D");
  nfitnt->Branch("errmigrad_deltacp",&errmigrad_deltacp,"errmigrad_deltacp/D");
  nfitnt->Branch("errmigrad_dm12sq",&errmigrad_dm12sq,"errmigrad_dm12sq/D");
  nfitnt->Branch("errmigrad_dm23sq",&errmigrad_dm23sq,"errmigrad_dm23sq/D");
  //errhesses
  nfitnt->Branch("errhesse_mixing12",&errhesse_mixing12,"errhesse_mixing12/D");
  nfitnt->Branch("errhesse_mixing13",&errhesse_mixing13,"errhesse_mixing13/D");
  nfitnt->Branch("errhesse_mixing23",&errhesse_mixing23,"errhesse_mixing23/D");
  nfitnt->Branch("errhesse_deltacp",&errhesse_deltacp,"errhesse_deltacp/D");
  nfitnt->Branch("errhesse_dm12sq",&errhesse_dm12sq,"errhesse_dm12sq/D");
  nfitnt->Branch("errhesse_dm23sq",&errhesse_dm23sq,"errhesse_dm23sq/D");


  ///////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////
  // Loop over toy-MC experiments
  ///////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////
  for(int iexperiment = first_experiment; iexperiment <= last_experiment; iexperiment++) {

    cout << "\n\n *** Pseudo-experiment = " << iexperiment << endl;
     
    // Re-seed fitter (reset current values to specified init values)
    fparm->ResetCurrToInitValues();

    // Get superk spectrum observation (toy experiment)
    SuperkObservation * skobs = new SuperkObservation;
    skobs->Read(xmlfile,iexperiment);

    // Get fitter
    FitInputs * fit_inputs = new FitInputs(skobs,skexp,fparm);
    FitterABC * fitter = new Superk1RmuErecoFitter(fit_inputs,fit_option);

    // Get nominal number of events
    nsk_nom = 0;
    for(int ir = 0; ir < kNErecoBins1Rmu; ir++)
      nsk_nom += skobs->ErecoSpectrumArray(conv::kSKsFCFV1RingMuLike)[ir];

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
    }

    // Chisq for true (input) params
    true_chisq = fit_inputs->Chisquare();

    // Print-out true params
    //fparm->ShowOnlyCurrValuesInNextPrint();
    //cout << "- True oscillation params (below) yield chisq = " << true_chisq
    //	   << *fparm;

    // Do the fit
    fparm->ResetCurrToInitValues(); // reset osc params
    fitter->DoTheFit();
    fitstatus       = fitter->FitStatusCode();
    migradfitstatus = fitter->FitMigradStatusCode();
    hessefitstatus  = fitter->FitHesseStatusCode();
    fittime = fitter->FitTime();
    bestfit_chisq = fit_inputs->BestFitChisquare();

    // Get bestfit number of events
    nsk_bf = 0;
    for(int ir = 0; ir < kNErecoBins1Rmu; ir++)
      nsk_bf += skexp->CurrPredictedErecoSpectrumArray(conv::kSKsFCFV1RingMuLike)[ir];

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
    // Init plus errors
    errminosplus_mixing12 = -1;
    errminosplus_mixing13 = -1;
    errminosplus_mixing23 = -1;
    errminosplus_deltacp  = -1;
    errminosplus_dm12sq   = -1;
    errminosplus_dm23sq   = -1;
    // Init minus errors
    errminosminus_mixing12 = -1;
    errminosminus_mixing13 = -1;
    errminosminus_mixing23 = -1;
    errminosminus_deltacp  = -1;
    errminosminus_dm12sq   = -1;
    errminosminus_dm23sq   = -1;
    // Init migrad errors
    errmigrad_mixing12 = -1;
    errmigrad_mixing13 = -1;
    errmigrad_mixing23 = -1;
    errmigrad_deltacp  = -1;
    errmigrad_dm12sq   = -1;
    errmigrad_dm23sq   = -1;
    // Init hesse errors
    errhesse_mixing12 = -1;
    errhesse_mixing13 = -1;
    errhesse_mixing23 = -1;
    errhesse_deltacp  = -1;
    errhesse_dm12sq   = -1;
    errhesse_dm23sq   = -1;
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
    }//hesse errors set?

    ///////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////
    // Fill tree
    ///////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////
    nfitnt->Fill();

    delete skobs;
    delete fit_inputs;
    delete fitter;
  }//iexperiment

  // Store ntuple with best-fit results
  f.cd();
  nfitnt->Write("nfitnt");

  f.Close();
}
