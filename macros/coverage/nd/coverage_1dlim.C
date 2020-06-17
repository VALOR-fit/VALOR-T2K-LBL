//
// ------------------------------------------------------------------------------------------------
// Generate and fit ensembles of toy-MC experiments with the purpose of calculating the coverage
// (or the critical \delta\chisq_crit for 1-dimensional confidence regions
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

void coverage_1dlim
(
  const char * tag                   = "toymc",                               /* tag for output file */
  int          run                   = 0,                                     /* run # */
  int          ntoymc                = 1,                                     /* number of toy-MC experiments */
  const char * hypname               = "3f_osc",                              /* fit parameter list (see UtlParam.cxx) */
  const char * param_values          = "randomize-all;prefit-errors",         /* decide how to set the values of params used for toy-MC generation */
  const char * lim_param_name        = 0,                                     /* specify the param name (see Parameter.h) for which a confidence limit is computed */
  int          lim_param_npoints     = 11,                                    /* coverage will be checked for a number of points within a range */
  double       lim_param_value_min   = 0.9,                                   /* min value in above-mentioned range */
  double       lim_param_value_max   = 1.1,                                   /* max value in avove-mentioned range */
  bool         stat_fluct            = true,                                  /* include statistical fluctuation in toy MC generation? */
  const char * fit_option            = "minuit",                              /* fitter option */
  bool         store_spectra         = false,                                 /* store toy-MC and best-fit spectra? (default: store only true and best-fit parameters) */
  const char * config_file           = "dune_ndfit-2015v1-FGTv1.xml",         /* analysis configuration file */
  const char * fit_param_opt         = ""                                     /* FitParamList options */
)

{
  using namespace std;
  using namespace valor;
  using namespace valor::conv;

  // Analysis configuration
  TString config = Form("%s/config/%s", gSystem->Getenv("VALOR_EXPT"), config_file);

  // Read the config file to populate the ParameterPool
  CompoundFitElements * elements = new CompoundFitElements();
  elements->Read(config.Data());

  assert(lim_param_name);    

  PhysParamId_t sp = util::NameRegistry::Instance()->PhysParamFromString(lim_param_name);
  cout << util::AsString(sp) << endl;
  if(sp == kPhypUndefined) {
    cerr << "Undefined parameter: " << lim_param_name << endl;
    exit(1);
  }

  //
  // Print settings
  //

  cout << "\n[Job inputs]" << endl;
  cout << "- Job: " << tag << "/" << run << endl;
  cout << "- Reading analysis configuration from: " << config << endl;
  cout << "- Studying sensitivity to parameter: " << lim_param_name << endl;
  cout << "- Toy-MC generation" << endl;
  cout << "   - Number of toy experiments: " << ntoymc << endl;
  cout << "   - Statistical fluctuations: "  << ((stat_fluct) ? "YES " : "NO") << endl;
  cout << "   - Option for setting parameter values: " << param_values << endl;
  cout << "- Fit" << endl;
  cout << "   - Physics hypothesis: " << hypname << endl;
  cout << "   - Fit Parameter List option: " << fit_param_opt << endl;
  cout << "   - Fitter option: " << fit_option << endl;
  cout << "- Coverage check" << endl;  
  cout << "   - For parameter: " << lim_param_name << endl;  
  cout << "   - Using " << lim_param_npoints << " points in range [" << lim_param_value_min << ", " << lim_param_value_max << "]" << endl;
  cout << endl;
  cout << *Rnd::Instance() << endl;

  //
  // Create file and tree to store outputs
  //

  const char * outfilename = Form("%s_%d.root",tag,run);
  TFile * outfile = new TFile(outfilename,"recreate");

  double          fit_statistic_true   = 0;
  PhysParamList * true_param_array     = new PhysParamList;
  FitParamList *  best_fit_param_array = new FitParamList;
  TObjArray *     best_fit_model_array = new TObjArray;
  TObjArray *     data_array           = new TObjArray;
  best_fit_model_array -> SetOwner(false);
  data_array           -> SetOwner(false);

  TTree * fitnt = new TTree("fitnt","VALOR toy-MC experiment fit results");
  fitnt->Branch("fit_statistic_true", &fit_statistic_true);
  fitnt->Branch("true_params",     "valor::PhysParamList", &true_param_array,     1280000,0);
  fitnt->Branch("best_fit_params", "valor::FitParamList",  &best_fit_param_array, 1280000,0);
  fitnt->Branch("best_fit_model",  "TObjArray",            &best_fit_model_array, 1280000,0);
  fitnt->Branch("data",            "TObjArray",            &data_array,           1280000,0);

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
    // The secondd option determined the `current' values which are used for toy-MC generatio below.
    FitParamList * fp = util::GetFitParameterList(hypname, fit_param_opt);

    // Make sure the specified parameter exists in the specified list of parameters
    if(! fp->Exists(sp)) {
       cerr << "Parameter: " << lim_param_name << " does not exist in list" << endl;
       exit(1);
    }

    // Choose randomly one of the possible values of the specified parameter
    int ival = (int) (lim_param_npoints * Rnd::Instance()->RndOther().Rndm());
    double lim_param_value = lim_param_value_min + ival * (lim_param_value_max - lim_param_value_min)/(lim_param_npoints-1);

    fp->Update      (sp, lim_param_value); // matters when I call SDModel::Generate()
    fp->SetInitValue(sp, lim_param_value); // matters when the fitter starts calling SDModel::Generate()
    cout << "Toy-MC generated for a fixed value of parameter: " << lim_param_name << " (= " << lim_param_value << ")" << endl;
    //cout << *fp << endl;

    // Generate a toy-MC experiment 
    cout << "Generating toy-MC data" << endl;
    FitDataABC * data = new SDData;
    data->Generate(model,*fp,stat_fluct);

    // The toy-MC was generated using the `current' values of the input parameter list.
    // Copy those values now (and they will be overwritten during the fit) to save 
    // them in the output fit file.
    true_param_array->Append(fp->AsPhysParamList("current"));
    cout << "Toy-MC #" << itoymc << " was generated with true parameters: " << *true_param_array; 

    // Calculate the fit-statistic for the true parameters
    fit_statistic_true = data->Statistic(model, *fp, "likelihood-ratio");
    cout << "fit_statistic (true parameters) = " << fit_statistic_true << endl;

    // Get fitter
    cout << "Preparing fitter" << endl;
    FitInputs * fit_inputs = new FitInputs(data,model,fp);
    FitterABC * fitter = util::GetFitter(fit_inputs,fit_option);

    // Do the fit
    cout << "Running fitter" << endl;
    cout << "Fit parameter list before running the fitter: " << *fp;
    fitter->DoTheFit();

    cout << "Fit parameter list after running the fitter: " << *fp;

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

    delete fit_inputs;
    delete fitter;
    delete fp;
    delete data;

  }//itoymc

  delete model;

  outfile->cd();
  fitnt->Write();
  outfile->Close();
}
