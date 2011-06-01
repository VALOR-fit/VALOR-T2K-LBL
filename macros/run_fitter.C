//
// ------------------------------------------------------------------------------------------------
// Run fitter on data
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

void run_fitter
(
  int          dataset_id  = 0,                            /* number in <data_set_%d> XML tag, see inputs/data.xml */
  const char * data_file   = "inputs/data.xml",            /* input SK data */
  const char * config_file = "analysis_configuration.xml", /* analysis configuration file */
  const char * hypothesis  = "disappearance",              /* ("disappearance", "decoherence", "decay")  */
  const char * fitter_name = "custom"                      /* fitter_name ("minuit", "custom" (2-D grid scan) */
)
{
  using namespace std;
  using namespace t2kana;

  cout << "Hypothesis : " << hypothesis << " / Fitter : " << fitter_name << endl;

  //
  // get superk spectrum predictor & configure
  //

  const char * config = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);
  cout << "Reading analysis configuration from: " << config << endl;

  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read(config);

  //
  // get superk spectrum observation
  //

  const char * data_file_fullpath = Form("%s/%s", gSystem->Getenv("T2KOSCANA"), data_file);
  cout << "Reading dataset: " << dataset_id << " from: " << data_file_fullpath << endl;

  SuperkObservation * skobs = new SuperkObservation;
  skobs->Read(data_file, dataset_id);

///  // fit params
//  PhysParamList * fparm = 0;

  //
  // disappearance fit
  //

  if(hypothesis == "disappearance") {
    PhysParamList * fparm = PhysParamList::NumuDisapFitNoSyst();
    FitInputs * fit_inputs = new FitInputs(skobs,skexp,fparm);
    NuMuDisappFitter * fitter = new NuMuDisappFitter(fit_inputs,fitter_name);

    // do the fit
    fitter->DoTheFit();

    // see inputs
    const char * outfile = Form("disap_fit_%d",dataset_id);
    fitter->Draw(new TCanvas(), "StandardOscFitPlots+Minos90CL",outfile);
    fitter->Draw(new TCanvas(), "ChisqSurface");
    fitter->Draw(new TCanvas(), "MarginalizedChisq");
  }

  //
  // decoherence fit
  //

  if(hypothesis == "decoherence") {
    PhysParamList * fparm = PhysParamList::NumuDecohFitNoSyst();
    FitInputs * fit_inputs = new FitInputs(skobs,skexp,fparm);
    NuMuDecohFitter * fitter = new NuMuDecohFitter(fit_inputs,fitter_name);

    // do the fit
    fitter->DoTheFit();

    // see inputs
    const char * outfile = Form("decoh_fit_%d",dataset_id);
    fitter->Draw(new TCanvas(), "StandardOscFitPlots",outfile);
    fitter->Draw(new TCanvas(), "ChisqSurface");
    fitter->Draw(new TCanvas(), "MarginalizedChisq");
  }

  //
  // decay fit
  //

  if(hypothesis == "decay") {
    PhysParamList * fparm = PhysParamList::NumuDecayFitNoSyst();
    FitInputs * fit_inputs = new FitInputs(skobs,skexp,fparm);
    NuMuDecayFitter * fitter = new NuMuDecayFitter(fit_inputs,fitter_name);

    // do the fit
    fitter->DoTheFit();

    // see inputs
    const char * outfile = Form("decay_fit_%d",dataset_id);
    fitter->Draw(new TCanvas(), "StandardOscFitPlots",outfile);
    fitter->Draw(new TCanvas(), "ChisqSurface");
    fitter->Draw(new TCanvas(), "MarginalizedChisq");
  }

}
