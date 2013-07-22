//
// ------------------------------------------------------------------------------------------------
// Run fits to a single experiment read from an input XML file and store best-fit parameters and spectra.
// Seed the fit with the chisq surface best fit, and run the fit multiple times to get the best best-fit possible
//
// Tom Dealtry
// ------------------------------------------------------------------------------------------------
//

void run_fit_xml_seeded
(
  const char * chisqfile,                                                   /* input ROOT file with chisq surface */
  const int    nfits,                                                       /* number of attempts to fit */         
  const char * tag,                                                         /* */
  int          run,                                                         /* */
  const char * xmlfile,                                                     /* input XML file with toy MC experiments */
  int          experiment_id,                                               /* */
  const char * config_file         = "analysis_configuration.xml",          /* analysis configuration file */
  const char * fit_param_list      = "nuosc3f_numudisapp_nuis4+nuisance",   /* fit parameter list for the given physics hypothesis (see FitParamList.h) */
  const char * fit_option          = "minuit",                              /* fitter option */
  const char * dmsqscheme          = "abs23",                               /* abs23==kNuDmsq12andAbs23 fogli==kNuDmsqFogliLisi */
  const char * hierarchy           = "normal",                              /* normal or inverted */
  const char * mixparamchoice      = "sinsq2theta",                         /* sinsq2theta or sinsqtheta */
  const char * mixparamoctant      = "one",                                 /* one==kNuOctantOne two==kNuOctantTwo */
  bool         apply_2nu_boundary  = true                                   /* if using the 2nu flavour osc prob, apply the physical boundary? */
)
{
  using namespace std;
  using namespace t2kana;

  // Analysis configuration
  const char * config = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);
  TString xmlfile2 = Form("%s/inputs/%s", gSystem->Getenv("T2KOSCANA"), xmlfile);

  // Global options
  Opt::Instance()->gMassHierarchyAssumed = util::HierarchyFromString(hierarchy);
  Opt::Instance()->gOscDmsqScheme        = util::DmsqSchemeFromString(dmsqscheme);
  Opt::Instance()->gOscMixingParamChoice = util::MixingTypeFromString(mixparamchoice);
  Opt::Instance()->gOscOctantChoice      = util::OctantFromString(mixparamoctant);
  if(Opt::Instance()->Undefined()) {
    cerr << "Global options not defined correctly. Exiting..." << endl;
    cerr << *Opt::Instance() << endl;
    exit(1);
  }

  // Print settings
  cout << "\n[Job inputs]" << endl;
  cout << "- Job: " << tag << "/" << run << endl;
  cout << "- Reading analysis configuration from: " << config << endl;
  cout << "- Toy-MC" << endl;
  cout << "   - Reading experiments from: " << xmlfile2 << endl; 
  cout << "   - Using experiments " << experiment_id << endl;
  cout << "- Fit" << endl;
  cout << "   - Fit param list: " << fit_param_list << endl;
  cout << "   - Fitter option: " << fit_option << endl;
  cout << endl;
  cout << *Opt::Instance() << endl;

  // Fit param list
  FitParamList * fparm = FitParamList::NamedParameterList(fit_param_list);

  conv::PhysHypothesis_t hyp = fparm->PhysicsHypothesis();
  cout << "- The following physics hypothesis is used in the fit: " << util::AsString(hyp) << endl;

  // Get superk spectrum predictor & configure
  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read(config);

  // Set whether the osc prob can be > 1 in 2 flavour world
  if(!apply_2nu_boundary && hyp == conv::kHypNuMuDisapp2F) {
    skexp->GetSuperkMC()->SetUnphysicalOscProb(true);
    fparm->SetMaxValue( conv::kPhypMixing23, 2);
  }

  // Arrays for storing all pieces of fit data and the best-fit pdf (and its components) for each expt.
  TObjArray * dat_array = new TObjArray;
  TObjArray * pdf_array = new TObjArray;
  dat_array->SetOwner(false);
  pdf_array->SetOwner(false);

  // Create file to store results
  TFile f(Form("%s-run_%d-expt_%d.root", tag, run, experiment_id), "recreate");

  // Open the file to read the chisq map from it
  TFile inputfile(chisqfile,"READ");
  TH2D *hchisq = 0;
  inputfile.GetObject("chisq_surface",hchisq);
  if(!hchisq) {
    cerr << "Chisq surface not found. Exiting..." << endl;
    return;
  }
  TH2D * hdchisq = new TH2D(*hchisq);
  
  // Find the minimum, and convert into mixing23 & dm23sq. This will give the seeds for the osc parameters
  int globalbin = hchisq->GetMinimumBin();
  int xbin, ybin, zbin;
  hchisq->GetBinXYZ(globalbin, xbin, ybin, zbin);
  double input_mixing23 = hchisq->GetXaxis()->GetBinCenter(xbin);
  double input_dm23sq   = hchisq->GetYaxis()->GetBinCenter(ybin);
  //think about the centre of symmetry
  double input_mixing23_below = input_mixing23, input_mixing23_above = input_mixing23;
  if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsqTheta) {
    double centre = 1 / (2 * (1 - fparm->InitValue(conv::kPhypMixing13)));
    if     (input_mixing23 < centre)
      input_mixing23_above = centre + TMath::Abs(centre - input_mixing23);
    else if(input_mixing23 > centre)
      input_mixing23_below = centre - TMath::Abs(centre - input_mixing23);
  }
  cout << "Best fit on chisq surface is " << hchisq->GetMinimum() << " at (mixing23,dm32sq) = "
       << input_mixing23 << "," << input_dm23sq << endl;

  // Seed the oscillation parameters
  fparm->SetInitValue(conv::kPhypMixing23, input_mixing23);
  fparm->SetInitValue(conv::kPhypDm23sq,   input_dm23sq);

  // Get superk spectrum observation (toy experiment)
  SuperkObservation * skobs = new SuperkObservation;
  skobs->Read(xmlfile2.Data(),experiment_id);

  FitParamList * global_bestfit_fparm = 0;
  double         global_bestfit_chisq = 99999999;

  for(int ifit = 0; ifit < nfits; ifit++) {
    cout << "Fit number " << ifit << endl;

    //randomise the seed (after fitting once with the default seeds)
    if(ifit != 0) {
      //randomise systematic seeds
      skexp->TweakAllSystematics();
      fparm->FromPhysParamList(skexp->Tweaks(),"initial", true);

      //randomise nuisance osc param seeds
      fparm->SetInitValue(conv::kPhypMixing12, util::RandomizeMixing12());
      fparm->SetInitValue(conv::kPhypMixing13, util::RandomizeMixing13());
      fparm->SetInitValue(conv::kPhypDm12sq,   util::RandomizeDm12sq());
      fparm->SetInitValue(conv::kPhypDeltaCP,  util::RandomizeDeltaCP());

      // Set the osc params seed to be within 3 bins of the chisq surfaces' best fit
      if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsqTheta)
	input_mixing23 = (Rnd::Instance()->RndOther().Uniform()) > 0.5 ? input_mixing23_below : input_mixing23_above;
      fparm->SetInitValue(conv::kPhypMixing23, Rnd::Instance()->RndOther().Uniform(input_mixing23 - 0.015,  TMath::Min(input_mixing23 + 0.015, 1.0)));
      fparm->SetInitValue(conv::kPhypDm23sq,   Rnd::Instance()->RndOther().Uniform(input_dm23sq   - 0.00015, input_dm23sq - 0.00015));
    }//ifit!=0

    // Re-seed fitter (reset current values to specified init values)
    fparm->ResetCurrToInitValues();

    // Get fitter
    FitInputs * fit_inputs = new FitInputs(skobs,skexp,fparm);
    FitterABC * fitter = util::GetFitter(fit_inputs,fit_option);

    // Do the fit
    //skexp->PrintSystematicVariations();
    fitter->DoTheFit();
    int fitstatus = fitter->FitStatusCode();
    double bestfit_chisq = fit_inputs->BestFitChisquare();

    // Print-out best-fit params
    fparm->ShowOnlyBestFitValuesInNextPrint();
    cout << "- Best-fit oscillation params (below) yield chisq = " << bestfit_chisq 
	 << *fparm;

    //get info to save
    if(bestfit_chisq < global_bestfit_chisq) {
      // Reset
      dat_array->Clear();
      pdf_array->Clear();
      delete global_bestfit_fparm;
      //add new
      skobs->Write(dat_array, "obs");
      skexp->Write(pdf_array, "pdf");
      global_bestfit_fparm = new FitParamList(*fparm);
    }

    //cleanup
    delete fit_inputs;
    delete fitter;
  }
  delete skobs;

  // Now we have the best fit, make the dchisq surface
  for(int ix = 1; ix <= hdchisq->GetNbinsX(); ix++) {
    for(int iy = 1; iy <= hdchisq->GetNbinsY(); iy++) {
      hdchisq->SetBinContent(ix, iy, hdchisq->GetBinContent(ix, iy) - bestfit_chisq);
    }//iy
  }//ix

  // Store everything (SKobs spectrum, SKexp spectrum, bestfit FitParamList, chisq surface, dchisq surface)
  f.cd();
  dat_array->Write();
  pdf_array->Write();
  global_bestfit_fparm->Write("fparm");
  hchisq->Write("chisq");
  hdchisq->Write("dchisq");
  hdchisq->SetContour(1);
  hdchisq->SetContourLevel(0,2.30);
  hdchisq->Write("dchisq_cont68");
  hdchisq->SetContourLevel(0,4.61);
  hdchisq->Write("dchisq_cont90");
  f.Close();
}
