//
// ------------------------------------------------------------------------------------------------
// Run numu disappearance systematic study for given oscillation parameters.
//	
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
// TD: Added 3 flavour oscillation support
//

//#include "Conventions.h"

void bestfit_systematic_study_fitting_job
(
   string       systematic,                                          /* tweaked systematic param */
   const char * fit_param_list      = "nuosc3f_numudisapp_fix4",     /* */
   int          irun                = 0,                             /* run # */
   double       input_mixing12      = -99,                           /* */
   double       input_mixing13      = -99,                           /* */
   double       input_mixing23      = -99,                           /* */
   double       input_deltacp       = -99,                           /* */
   double       input_dm12sq        = -99,                           /* */
   double       input_dm23sq        = -99,                           /* */
   int          nexperiments        = 500,                           /* number of toy MC expt fits in this run */
   const char * config_file         = "analysis_configuration.xml"   /* analysis configuration */
)
{
  using namespace std;
  using namespace valor;

  TRandom3 rand(0);

  //get appropriate input osc parameters
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

  // Analysis configuration
  const char * config = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);
  cout << "Reading analysis configuration from: " << config << endl;

  // Global options
  //Opt::Instance()->gOscMixingParamChoice = conv::kNuMixSinsq2Theta;
  //Opt::Instance()->gMassHierarchyAssumed = conv::kNmhNormal;
  //cout << *Opt::Instance() << endl;

  // Create ntuple to store results
  TNtupleD * systnt = new TNtupleD("sysnt","", 
      "irun:iexp:input_mixing12:input_mixing13:input_mixing23:input_dm12sq:input_dm23sq:mixing12:mixing13:mixing23:deltacp:dm12sq:dm23sq:chisq:status");

  // Get superk spectrum predictor & configure
  SuperkExpectation skexp;
  skexp.Read(config);

  // Fit params
  FitParamList * fparm = util::NamedParameterList(fit_param_list);
  conv::PhysHypothesis_t hyp = fparm->PhysicsHypothesis();

  // Initialise bestfits for the oscillation parameters
  double bestfit_mixing12 = -10, bestfit_mixing13 = -10, bestfit_mixing12 = 0;
  double bestfit_deltacp = -10, bestfit_dm12sq = -10, bestfit_dm23sq = 0;

  for(int iexperiment = 0; iexperiment < nexperiments; iexperiment++) {

     cout << " ******* experiment = " << iexperiment << endl;

     // reset oscillation params following last fit
     if(hyp == conv::kHypOsc3F) {
       fparm->Update(conv::kPhypMixing12, input_mixing12);
       fparm->Update(conv::kPhypMixing13, input_mixing13);
       fparm->Update(conv::kPhypCPV    ,  input_deltacp);
       fparm->Update(conv::kPhypDm12sq,   input_dm12sq);
     }
     fparm->Update(conv::kPhypMixing23, input_mixing23);
     fparm->Update(conv::kPhypDm23sq,   input_dm23sq);

     // tweak systematic in spectrum predictor
     skexp.TweakSystematic(systematic, rand.Gaus(0.,1.));

     // generate mock data using tweaked spectrum predictor - don't include statistical fluctuations
     SuperkObservation skobs;
     skobs.CreateMockData(skexp, *fparm, false);       

     // reset systematic variations so that the nominal spectrum prediction is used for the fit
     skexp.ResetSystematicVariations();

     // get fitter
     FitInputs fit_inputs(&skobs, &skexp, fparm);
     Superk1RmuErecoFitter fitter(&fit_inputs,"minuit");

     // do the fit
     bool ok = fitter.DoTheFit();

     // store params
     if(hyp == conv::kHypOsc3F) {
       bestfit_mixing12 = fparm->CurrValue(conv::kPhypMixing12);
       bestfit_mixing13 = fparm->CurrValue(conv::kPhypMixing13);
       bestfit_deltacp  = fparm->CurrValue(conv::kPhypCPV    );
       bestfit_dm12sq   = fparm->CurrValue(conv::kPhypDm12sq);
     }
     double bestfit_mixing23 = fparm->CurrValue(conv::kPhypMixing23);
     double bestfit_dm23sq   = fparm->CurrValue(conv::kPhypDm23sq);
     double bestfit_chisq    = fit_inputs.Chisquare();
     double fit_status       = (ok) ? +1 : -1;

     systnt->Fill(irun, 
                  iexperiment,
		  input_mixing12,
		  input_mixing13,
                  input_mixing23,
		  input_deltacp,
		  input_dm12sq, 
                  input_dm23sq,
		  bestfit_mixing12,
		  bestfit_mixing13,
                  bestfit_mixing23,
		  // bestfit_deltacp, //can only call this function with max 15 arguments
		  bestfit_dm12sq, 
                  bestfit_dm23sq,
                  bestfit_chisq, 
                  fit_status);
  }


  string oscparams = util::AsString(hyp,
	  input_mixing12,input_mixing13,input_mixing23,input_deltacp,input_dm12sq,input_dm23sq);
  TFile f(Form("syst_%s-%s_%d.root",
          systematic.c_str(),oscparams.c_str(),irun), "recreate");
  systnt->Write("sysnt");
  f.Close();
}
