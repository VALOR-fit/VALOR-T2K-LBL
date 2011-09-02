///////////////////////////////////
// A macro for testing purposes only.
// "makeFitInputs" gets a FitInputs by either creating mock data or reading in an xml file
// The main function then uses this to create a chisq surface by go over each point in the 23-grid
//  and getting the chisq. Will use minuit to find the minima if required (if there are any free parameters left)
// The chisq surface building code will be added to NuMuDisappFitter when validated
//
// Tom Dealtry
///////////////////////////////////

#include "Conventions.h"
using namespace t2kana;

FitInputs & makeFitInputs(
  const char * fit_param_list          = "nuosc3f_numudisapp_marginalize",/* fit parameter list for the given physics hypothesis (see run_fitter.C) */
  bool         calculatebestfit        = true,                            /* if true, calculate best fit */
  bool         usexml                  = false,                           /* if true, read xml file, otherwise make toy expts*/
  double       input_sinsq2theta12     = 0.86,                            /* if -1 then randomize */
  double       input_sinsq2theta13     = 0.10,                            /* if -1 then randomize */
  double       input_sinsq2theta23     = 1.00,                            /* */
  double       input_deltacp           = 0.00,                            /* if -1 then randomize */
  double       input_dm12sq            = 7.58E-5,                         /* if -1 then randomize */
  double       input_dm23sq            = 2.32E-3,                         /* */
  int          dataset_id              = 0,                               /* number in <data_set_%d> XML tag, see inputs/data.xml */
  const char * data_file               = "inputs/data.xml",               /* input SK data */
  const char * fitter_name             = "minuit",                        /* fitter name (see run_fitter.C) */
  bool         fluctuate_stat          = true,                            /* include statistical fluctuations? */
  bool         include_systematics     = false,                           /* include systematic variation in generated toy MC expt? */
  const char * config_file             = "analysis_configuration.xml"     /* analysis configuration file */
)
{
  FitParamList * fparam = FitParamList::NamedParameterList(fit_param_list);
  conv::PhysHypothesis_t hyp = fparam->PhysicsHypothesis();

  // get superk spectrum predictor & configure
  SuperkExpectation * skexp = new SuperkExpectation;
  const char * config = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);
  skexp->Read(config);

  SuperkObservation * skobs = new SuperkObservation;
  if(usexml) {
    skobs->Read(data_file, dataset_id);
  }//usexml?
  else {
    // set oscillation parameter inputs as necessary (used to generate toy-MC expt)
    if(hyp == conv::kHypOsc3F) {
      fparam->Update(t2kana::conv::kPhypSinsq2Theta12, input_sinsq2theta12);
      fparam->Update(t2kana::conv::kPhypSinsq2Theta13, input_sinsq2theta13);
      fparam->Update(t2kana::conv::kPhypSinsq2Theta23, input_sinsq2theta23);
      fparam->Update(t2kana::conv::kPhypDeltaCP,       input_deltacp);
      fparam->Update(t2kana::conv::kPhypDm12sq,        input_dm12sq);
      fparam->Update(t2kana::conv::kPhypDm23sq,        input_dm23sq);
    }//3flavour?
    else if(hyp == conv::kHypNuMuDisapp2F) {
      fparam->Update(t2kana::conv::kPhypSinsq2Theta23, input_sinsq2theta23);
      fparam->Update(t2kana::conv::kPhypDm23sq,        input_dm23sq);
    }//2flavour?
    else {
      cout << "Unknown physics hypothesis, exiting" << endl;
      return;
    }//unknown?

    //and generate the toy mc
    if(include_systematics) {
      // tweak all sources of systematic err (gaussian err)
      skexp->TweakAllSystematics();     
      // generate mock data using the current PDF (incl. systematic variations) 
      // + statistical fluctuations (Poisson fluctuation in each bin)
      skobs->CreateMockData(*skexp, *fparam, fluctuate_stat);
      // reset all systematics so that the PDF used in the fitter is the nominal one
      skexp->ResetSystematicVariations();
    }//include_systematics?
    else {       
      // generate mock data using the current PDF (without systematic variations) 
      // with statistical fluctuations (Poisson fluctuation in each bin)
      skobs->CreateMockData(*skexp, *fparam, fluctuate_stat); 
    }//include_systematics?
  }//usexml?

  // get fitter
  FitInputs * fit_inputs = new FitInputs(skobs,skexp,fparam);

  // get specified fit parameter list & build fit inputs
//   FitParamList * fparam = FitParamList::NamedParameterList(fit_param_list);
//   FitInputs * fit_inputs = new FitInputs(skobs,skexp,fparam);

  if(calculatebestfit) {
    // get the fitter
    NuMuDisappFitter * fitter = new NuMuDisappFitter(fit_inputs,fitter_name);
    
    // do the fit
    fitter->DoTheFit();
    int bestfit_fitstatus = fitter->FitStatusCode();

    // get the best-fit values
    double bestfit_sinsq2theta12 = -1;
    double bestfit_sinsq2theta13 = -1;
    double bestfit_sinsq2theta23 = -1;
    double bestfit_deltacp       = -1;
    double bestfit_dm12sq        = -1;
    double bestfit_dm23sq        = -1;
    double bestfit_chisq         = -1;
    conv::PhysHypothesis_t hyp = fparam->PhysicsHypothesis();
    if(hyp == conv::kHypOsc3F) {
      bestfit_sinsq2theta12 = fparam->CurrValue(t2kana::conv::kPhypSinsq2Theta12);
      bestfit_sinsq2theta13 = fparam->CurrValue(t2kana::conv::kPhypSinsq2Theta13);
      bestfit_sinsq2theta23 = fparam->CurrValue(t2kana::conv::kPhypSinsq2Theta23);
      bestfit_deltacp       = fparam->CurrValue(t2kana::conv::kPhypDeltaCP);
      bestfit_dm12sq        = fparam->CurrValue(t2kana::conv::kPhypDm12sq);
      bestfit_dm23sq        = fparam->CurrValue(t2kana::conv::kPhypDm23sq);
      cout << "Best-fit:"
	   << " sinsq2theta12 = "  << bestfit_sinsq2theta12
	   << ", sinsq2theta13 = " << bestfit_sinsq2theta13
	   << ", sinsq2theta23 = " << bestfit_sinsq2theta23
	   << ", delta_cp = " << bestfit_deltacp      
	   << ", Dm12^2 = "   << bestfit_dm12sq << " eV^2/c^4"
	   << ", |Dm23^2| = " << bestfit_dm23sq << " eV^2/c^4"
	   << " ==> chisq = " << bestfit_chisq << endl;
    }//3f?
    else if(hyp == conv::kHypNuMuDisapp2F) {
      bestfit_sinsq2theta23 = fparam->CurrValue(t2kana::conv::kPhypSinsq2Theta23);
      bestfit_dm23sq        = fparam->CurrValue(t2kana::conv::kPhypDm23sq);
      cout << "Best-fit:"
	   << " sinsq2theta23 = " << bestfit_sinsq2theta23
	   << ", |Dm23^2| = " << bestfit_dm23sq << " eV^2/c^4"
	   << " ==> chisq = " << bestfit_chisq << endl;
    }//2f?
  }//calculatebestfit?



  return *fit_inputs;
}


void build_chisq_surface(
 		  const char * tag                     = "chisq_surface",                    /* */
 		  double       sinsq2theta_min_bin     = 81,                               /*4 variables to define the size of the chisq surface you want to build*/
 		  double       sinsq2theta_max_bin     = 101,                              /* */
 		  double       dmsq_min_bin            = 11,                               /* */
 		  double       dmsq_max_bin            = 29,                               /* */
 		  const char * fit_param_list          = "nuosc3f_numudisapp_marginalize", /* fit parameter list for the given physics hypothesis (see run_fitter.C) */
		  bool         calculatebestfit        = true,
 		  const char * fitter_name             = "minuit",                         /* fitter name (see run_fitter.C) */
 		  const char * config_file             = "analysis_configuration.xml"      /* analysis configuration file */
		  ) 
{
  using namespace std;
  using namespace t2kana;
  using namespace t2kana::conv;
  
  // get the fit inputs and fitter
  FitInputs         * fit_inputs = 0;
  fit_inputs = makeFitInputs(fit_param_list);
  FitParamList      * fparam      = 0;
  fparam      = fit_inputs->FitParams();
  SuperkExpectation * skexp      = 0;
  skexp      = fit_inputs->Expectation();
  SuperkObservation * skobs      = 0;
  skobs      = fit_inputs->Observation();

  // create file to store results
  TFile f(Form("%s.root",tag), "recreate");

  //make a histogram to store the chisqsurface
  double sinsq2theta23_min = 0.;//fit_inputs->FitParams().MinValue(t2kana::conv::kPhypSinsq2Theta23);
  double sinsq2theta23_max = 1.;//fit_inputs->FitParams().MaxValue(t2kana::conv::kPhypSinsq2Theta23);
  double sinsq2theta23_step = 0.01;//fit_inputs->FitParams().StepSize(t2kana::conv::kPhypSinsq2Theta23);
  double dm23sq_min = 1e-3;//fit_inputs->FitParams().MinValue(t2kana::conv::kPhypDm23sq);
  double dm23sq_max = 6e-3;//fit_inputs->FitParams().MaxValue(t2kana::conv::kPhypDm23sq);
  double dm23sq_step = 1e-4;//fit_inputs->FitParams().StepSize(t2kana::conv::kPhypDm23sq);
  // make sure that min and max are at the bin centre
  int sinsq2theta23_nbins = 1 + (int) ((sinsq2theta23_max - sinsq2theta23_min) / sinsq2theta23_step);
  int dm23sq_nbins = 1 + (int) ((dm23sq_max - dm23sq_min ) / dm23sq_step );
  TH2D *ChisqSurface = new TH2D("chisqsurf","#chi^{2} surface",
	 sinsq2theta23_nbins, sinsq2theta23_min - 0.5*sinsq2theta23_step, sinsq2theta23_max + 0.5*sinsq2theta23_step,
	 dm23sq_nbins, dm23sq_min - 0.5*dm23sq_step, dm23sq_max + 0.5*dm23sq_step);
  ChisqSurface->SetDirectory(0);
  ChisqSurface->GetXaxis()->SetRangeUser(0.,1.);

  //fix the oscillation parameters of the grid points
  fparam->FixParameter( t2kana::conv::kPhypSinsq2Theta23 );
  fparam->FixParameter( t2kana::conv::kPhypDm23sq );

  bool needtofit = fparam->GetNumberFreeParameters() > 0;
  cout << fparam->GetNumberFreeParameters() << " free parameters to fit" << endl;

  for(int i = 1; i <= ChisqSurface->GetNbinsX(); i++) {
    bool in_range_theta = (i >= sinsq2theta_min_bin && i <= sinsq2theta_max_bin);
    if(!in_range_theta) continue;
    cout << "___________for loop for bin " << i << " in x" << endl;

    for(int j = 1; j <= ChisqSurface->GetNbinsY(); j++) {
      bool in_range_dm = (j >= dmsq_min_bin && j <= dmsq_max_bin);
      if(!in_range_dm) continue;

      int fitstatus = 0;
      fparam->Update( t2kana::conv::kPhypSinsq2Theta23, ChisqSurface->GetXaxis()->GetBinCenter(i) );
      fparam->Update( t2kana::conv::kPhypDm23sq, ChisqSurface->GetYaxis()->GetBinCenter(j) );
      if (needtofit) {
	fparam->SetInitValue( conv::kPhypSinsq2Theta23, ChisqSurface->GetXaxis()->GetBinCenter(i) );
	fparam->SetInitValue( conv::kPhypDm23sq, ChisqSurface->GetYaxis()->GetBinCenter(j) );
	fparam->SetInitValue( t2kana::conv::kPhypSinsq2Theta13, bestfit_sinsq2theta13);
	fparam->SetInitValue( t2kana::conv::kPhypSinsq2Theta12, bestfit_sinsq2theta12);
	fparam->SetInitValue( t2kana::conv::kPhypDeltaCP, bestfit_deltacp);
	fparam->SetInitValue( t2kana::conv::kPhypDm12sq, bestfit_dm23sq);

	FitInputs * surfacefit_inputs = new FitInputs(skobs,skexp,fparam);
	NuMuDisappFitter * surfacefitter = new NuMuDisappFitter(surfacefit_inputs,fitter_name);
	surfacefitter->DoTheFit();
	fitstatus = surfacefitter->FitStatusCode(), refitcount = 0;

	//refit if we need to
	while( fitstatus!=3 && refitcount<10 ) {
	  cout << "______________refitting____________" << endl;
	  fparam->SetInitValue( t2kana::conv::kPhypSinsq2Theta13, util::RandomizeSinsq2Theta13() );
	  fparam->SetInitValue( t2kana::conv::kPhypSinsq2Theta12, util::RandomizeSinsq2Theta12() );
	  fparam->SetInitValue( t2kana::conv::kPhypDeltaCP, util::RandomizeDeltaCP() );
	  fparam->SetInitValue( t2kana::conv::kPhypDm12sq, util::RandomizeDm12sq() );
	  delete surfacefit_inputs;
	  delete surfacefitter;
	  FitInputs * surfacefit_inputs = new FitInputs(skobs,skexp,fparam);
	  NuMuDisappFitter * surfacefitter = new NuMuDisappFitter(surfacefit_inputs,fitter_name);
	  surfacefitter->DoTheFit();
	  fitstatus = surfacefitter->FitStatusCode();
	  refitcount++;
	}//refitting loop
      }//needtofit?
      else
	fitstatus = 3;

      double chisq = fit_inputs->Chisquare();
      if(fitstatus!=3)
	cout << "Fit failure " << fitstatus << "for sin:dm = " << ChisqSurface->GetXaxis()->GetBinCenter(i) << ":" 
	     << ChisqSurface->GetYaxis()->GetBinCenter(i) << endl;
      ChisqSurface->SetBinContent(i,j,chisq);
      if(needtofit) {
	delete surfacefit_inputs;
	delete surfacefitter;
      }//needtofit?
    }//j
  }//i

  ChisqSurface->GetXaxis()->SetTitle("sin^{2}(2#theta_{23})");
  ChisqSurface->GetYaxis()->SetTitle("|#Deltam_{32}^{2}| (eV^{2}/c^{4})");
  ChisqSurface->GetYaxis()->SetTitleOffset(1);
  ChisqSurface->GetYaxis()->SetTitleSize(0.05);
  ChisqSurface->GetXaxis()->SetTitleSize(0.05);

  f.cd();
  ChisqSurface->Draw("colz");
  ChisqSurface->Write();
  f.Close();
}
