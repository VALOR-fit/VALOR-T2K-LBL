//
// ------------------------------------------------------------------------------------------------
// Generate and fit ensembles of toy-MC experiments.
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//
	
void run_fit
(
  const char * tag              = "toymc",                               /* tag for output file */
  int          run              = 0,                                     /* run # */
  int          ntoymc           = 1,                                     /* number of toy-MC experiments */
  const char * hypname          = "3f_osc",                              /* physics hypothesis (see UtlParam.cxx) */
  const char * param_values     = "randomize-all;prefit-errors",         /* decide how to set the values of params used for toy-MC generation */
  bool         stat_fluct       = true,                                  /* include statistical fluctuation in toy MC generation? */
  const char * fit_option       = "minuit",                              /* fitter option */
  bool         store_spectra    = true,                                  /* store toy-MC and best-fit spectra? (default: store only true and best-fit parameters) */
  const char * config_file      = "dune_ndfit-2015v1-FGTv1.xml",         /* analysis configuration file */
  const char * fit_param_opt    = ""                                     /* FitParamList options */
)

{
  using namespace std;
  using namespace valor;
  using namespace valor::conv;

  // Analysis configuration
  TString config = Form("%s/config/%s", gSystem->Getenv("VALOR_EXPT"), config_file);
  
  TString param_values_ts(param_values); // TString version for useful methods

  //
  // Print settings
  //

  cout << "\n[Job inputs]" << endl;
  cout << "- Job: " << tag << "/" << run << endl;
  cout << "- Reading analysis configuration from: " << config << endl;
  cout << "- Toy-MC generation" << endl;
  cout << "   - Number of toy experiments: " << ntoymc << endl;
  cout << "   - Statistical fluctuations: "  << ((stat_fluct) ? "YES " : "NO") << endl;
  cout << "   - Option for setting parameter values: " << param_values << endl;
  cout << "- Fit" << endl;
  cout << "   - Fit Parameter List: " << hypname << endl;
  cout << "   - Fitter option: " << fit_option << endl;
  cout << endl;
  cout << *Rnd::Instance() << endl;

  //
  // Create file and tree to store outputs
  //

  const char * outfilename = Form("%s_%d.root",tag,run);
  TFile * outfile = new TFile(outfilename,"recreate");

  PhysParamList * true_param_array     = new PhysParamList;
  FitParamList *  best_fit_param_array = new FitParamList;
  TObjArray *     best_fit_model_array = new TObjArray;
  TObjArray *     data_array           = new TObjArray;
  best_fit_model_array -> SetOwner(false);
  data_array           -> SetOwner(false);
  int migradmin, migraderr;

  TTree * fitnt = new TTree("fitnt","VALOR toy-MC experiment fit results");
  fitnt->Branch("true_params",     "valor::PhysParamList", &true_param_array,     1280000,0);
  fitnt->Branch("best_fit_params", "valor::FitParamList",  &best_fit_param_array, 1280000,0);
  fitnt->Branch("best_fit_model",  "TObjArray",            &best_fit_model_array, 1280000,0);
  fitnt->Branch("data",            "TObjArray",            &data_array,           1280000,0);
//  fitnt->Branch("migradmin",                               &migradmin);
//  fitnt->Branch("migraderr",                               &migraderr);

  // Get ND model & configure
  cout << "Getting and configuring fit model" << endl;
  FitModelABC * model = new SDModel;
  model->Read(config.Data());

  // Use ND model to generate toy-MC experiments & fit
  for(int itoymc = 0; itoymc < ntoymc; itoymc++) {
    cout << "\n\n *** Toy-MC experiment #" << itoymc << endl;

    // Reset arrays
    true_param_array     -> Clear();
    best_fit_param_array -> Clear();
    best_fit_model_array -> Clear();
    data_array           -> Clear();

    // Getting the fit parameter list that corresponds to the input option.
    // The second option determined the `current' values which are used for toy-MC generation below.
    FitParamList * fp = util::GetFitParameterList(hypname, fit_param_opt);
    
    
    // Randomise the systematics in the model, and copy them to the FPL
    if (param_values_ts.Contains("randomize")) {
      model->RandomizeSystematics(/*param_values*/);
    }
    fp->FromPhysParamList(model->Tweaks(),"current");
    
    // Generate a toy-MC experiment 
    cout << "Generating toy-MC data" << endl;
    FitDataABC * data = new SDData;
    data->Generate(model,*fp,stat_fluct);

    // The toy-MC was generated using the `current' values of the input parameter list.
    // Copy those values now (and they will be overwritten during the fit) to save 
    // them in the output fit file.
    true_param_array->Append(fp->AsPhysParamList("current"));
    
    fp->ResetCurrToInitValues();

    // Get fitter
    cout << "Preparing fitter" << endl;
    FitInputs * fit_inputs = new FitInputs(data,model,fp);
    FitterABC * fitter = util::GetFitter(fit_inputs,fit_option);

    // Do the fit
    cout << "Running fitter" << endl;
    fitter->DoTheFit();

    migradmin = fitter->FitStatusCode("migrad_min");
    migraderr = fitter->FitStatusCode("migrad_err");

    // Copy best-fit params
    best_fit_param_array->Copy(*fp);

    if(store_spectra) {
      // Copy best-fit to current values and re-generate the best-fit prediction
      // with special option (not used in the fit) to generate versions of the fit model
      // using a fine binning scheme (for plots).
      // Then, write-out the best-fit predictions
      fp->CopyBestFitToCurrValues();
      model->Generate(*fp,"+SB");
      model->Write(best_fit_model_array,"best_fit","pdfc_rf,pdft_rf,pdfts_rf,pdfcs_rf +SB");

      // Write-out the toy-MC
      data->Write(data_array,"toymc","all");
    }

    // Write to tree
    fitnt->Fill();

    delete data;
    delete fp;
    delete fit_inputs;
    delete fitter;

  }//itoymc

  outfile->cd();
  fitnt->Write();
  outfile->Close();
}
