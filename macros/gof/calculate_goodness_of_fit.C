//
// ------------------------------------------------------------------------------------------------
// Run numu disappearance fits in the 5 bin mode, calculate the best fit oscillation parameters,
// then run N (~1000) toy experiments with those values (throwing stat&syst).
// Output is a plot of chisq of each of the N fits, with the bestfit of the xml data marked,
// and the probabilty of the data occuring.
//
// Tom Dealtry
// ------------------------------------------------------------------------------------------------
//

void calculate_goodness_of_fit
(
  const char * tag,                                                             /* */
  const char * datafit_file,                                                    /* root file containing output of run_fit_xml_seeded.C */
  const int    ntoys,                                                           /* number of toys to run */
  const char * xmlfile,                                                         /* input XML file with toy MC experiments */
  int          iexperiment,                                                     /* */
  const char * config_file             = "analysis_configuration.xml",          /* analysis configuration file */
  const char * fit_param_list          = "nuosc3f_numudisapp_fix4+nuisance",    /* fit parameter list for the given physics hypothesis (see FitParamList.h) */
  const char * fit_option              = "minuit",                              /* fitter option */
  bool         store_toy_expts         = true,                                  /* disable this if you don't want to store the toy MC expts     */
  bool         store_bestfit_spectrum  = true,                                  /* disable this if you don't want to store the best-fit spectra */
  const char * dmsqscheme              = "abs23",                               /* abs23==kNuDmsq12andAbs23 fogli==kNuDmsqFogliLisi */
  const char * hierarchy               = "normal",                              /* normal or inverted */
  const char * mixparamchoice          = "sinsqtheta",                          /* sinsq2theta or sinsqtheta */
  const char * mixparamoctant          = "one"                                  /* one==kNuOctantOne two==kNuOctantTwo */
  //bool         apply_2nu_boundary      = true                                   /* if using the 2nu flavour osc prob, apply the physical boundary? */
)
{
  using namespace std;
  using namespace t2kana;

  // Analysis configuration
  TString config   = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);
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
  cout << "- Job: " << tag << endl;
  cout << "- Reading analysis configuration from: " << config << endl;
  cout << "- Reading experiment: " << iexperiment << " from: " << xmlfile2 << endl; 
  cout << "- Reading data fit result from: " << datafit_file << endl; 
  cout << "- Fit" << endl;
  cout << "   - Fit param list: " << fit_param_list << endl;
  cout << "   - Fitter option: " << fit_option << endl;
  cout << endl;
  cout << *Opt::Instance() << endl;

  // Fit param list
  FitParamList * fparm  = FitParamList::NamedParameterList(fit_param_list); //only used to initialise branch
  FitParamList * fparm1 = FitParamList::NamedParameterList(fit_param_list);
  FitParamList * fparm2 = FitParamList::NamedParameterList(fit_param_list);

  conv::PhysHypothesis_t hyp = fparm1->PhysicsHypothesis();
  cout << "- The following physics hypothesis is used in the fit: " << util::AsString(hyp) << endl;

  // Get superk spectrum predictor & configure
  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read(config.Data());

  // Object for storing input oscillation and nuisance params in each toy-MC
  PhysParamList * iparm = new PhysParamList;

  // Arrays for storing all pieces of fit data and the best-fit pdf (and its components) for each expt.
  TObjArray * dat_array = new TObjArray; //only used to initialise branch
  TObjArray * pdf_array = new TObjArray; //only used to initialise branch
  dat_array->SetOwner(false);
  pdf_array->SetOwner(false);
  TObjArray * dat_array1 = new TObjArray;
  TObjArray * pdf_array1 = new TObjArray;
  dat_array1->SetOwner(false);
  pdf_array1->SetOwner(false);
  TObjArray * dat_array2 = new TObjArray;
  TObjArray * pdf_array2 = new TObjArray;
  dat_array2->SetOwner(false);
  pdf_array2->SetOwner(false);

  // Create file to store results
  TFile f(Form("gof-%s-expt_%d.root", tag, iexperiment), "recreate");

  // Create tree to store outputs
  TTree * fitnt = new TTree("fitnt","Results of VaLOR goodness of fit test");
  fitnt->Branch("best_fit_params",      "t2kana::FitParamList",  &fparm,      1280000,0);
  fitnt->Branch("best_fit_params_fit1", "t2kana::FitParamList",  &fparm1,     1280000,0);
  fitnt->Branch("best_fit_params_fit2", "t2kana::FitParamList",  &fparm2,     1280000,0);
  fitnt->Branch("input_params",         "t2kana::PhysParamList", &iparm,      1280000,0);
  fitnt->Branch("data",                 "TObjArray",             &dat_array,  1280000,0);
  fitnt->Branch("best_fit_pdf",         "TObjArray",             &pdf_array,  1280000,0);
  fitnt->Branch("data_fit1",            "TObjArray",             &dat_array1, 1280000,0);
  fitnt->Branch("best_fit_pdf_fit1",    "TObjArray",             &pdf_array1, 1280000,0);
  fitnt->Branch("data_fit2",            "TObjArray",             &dat_array2, 1280000,0);
  fitnt->Branch("best_fit_pdf_fit2",    "TObjArray",             &pdf_array2, 1280000,0);
  TBranch * br_fparm     = fitnt->GetBranch("best_fit_params");
  TBranch * br_obs_array = fitnt->GetBranch("data");
  TBranch * br_bf_array  = fitnt->GetBranch("best_fit_pdf");

  // Get the data-fit best-fit point
  TFile fin(datafit_file);
  FitParamList * fparm_input = 0;
  fin.GetObject("fparm",fparm_input);
  if(!fparm_input) {
    cerr << "Input FitParamList not found. Exiting..." << endl;
    exit();
  }
  double bestfit_mixing23 = fparm_input->BestFitValue(conv::kPhypMixing23);
  double bestfit_dm32sq   = fparm_input->BestFitValue(conv::kPhypDm23sq);
  //and get the corresponding 5 bin chisq
  fparm_input->CopyBestFitToCurrValues();
  SuperkObservation * skobs_input = new SuperkObservation;
  skobs_input->Read(xmlfile2.Data(),iexperiment);
  skexp->Generate(*fparm_input);
  FitInputs * fit_inputs_input = new FitInputs(skobs_input,skexp,fparm_input);
  fit_inputs_input->SetChisqDefinition(true);
  double bestfit_chisq    = fit_inputs_input->Chisquare();
  delete skobs_input;
  delete fit_inputs_input;

  //book histogram
  TH1D * hgof = new TH1D("gof", "gof", 100000, 0, 1000);
  hgof->GetXaxis()->SetTitle("#chi^{2}_{gof}");

  //counter for how many chisqs are greater than the data chisq
  int dataset_counter = 0;

  // Loop over toy-MC experiments, generated on the fly, at the data best-fit point
  for(int itoy = 0; itoy < ntoys; itoy++) {

     cout << "\n\n *** Fit = " << itoy << endl;

     // Reset output spectra
     dat_array ->Clear();
     pdf_array ->Clear();
     dat_array1->Clear();
     pdf_array1->Clear();
     dat_array2->Clear();
     pdf_array2->Clear();

     // Re-seed fitter (reset current values to specified init values)
     fparm1->ResetCurrToInitValues();

     // Get superk spectrum observation (toy experiment)
     SuperkObservation * skobs = new SuperkObservation;
     // tweak all sources of systematic error
     cout << "- Randomizing all systematic params in pdf..." << endl;
     skexp->TweakAllSystematics();     
     skexp->PrintSystematicVariations();
     // keep track of all nuisance params set so as to add in output tree
     PhysParamList true_nuisance_params = skexp->Tweaks();
     iparm->Append(true_nuisance_params);
     // generate mock data using the current PDF (incl. systematic variations) 
     // + statistical fluctuations if requested (Poisson fluctuation in each bin)
     cout << "- Generating toy-MC..." << endl;
     //note that fparm does not contain the tweaks, so use only the PhysHyp here
     bool stat_fluctuations = true;
     skobs->CreateMockData(*skexp, PhysHypothesis(*fparm1), stat_fluctuations);
     // reset all systematics so that the PDF used in the fitter is the nominal one
     cout << "- Reseting all systematic params in pdf..." << endl;
     skexp->ResetSystematicVariations();
     //skexp->PrintSystematicVariations();

     // Get fitter
     FitInputs * fit_inputs = new FitInputs(skobs,skexp,fparm1);
     fit_inputs->SetChisqDefinition(true); //5-bin chisq
     FitterABC * fitter = util::GetFitter(fit_inputs,fit_option);

     // Do the fit
     skexp->PrintSystematicVariations();
     fitter->DoTheFit();
     int fitstatus = fitter->FitStatusCode();
     double toy_chisq = fit_inputs->BestFitChisquare();

     // Save spectra (data and best-fit prediction)
     if(store_toy_expts) {
        skobs->Write(dat_array1, "obs");
     }
     if(store_bestfit_spectrum) {
        skexp->Write(pdf_array1, "pdf");
     }

     //Do the fit again
     if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsqTheta) {
       //get the correct intial values (best fit points from first fit)
       fparm2->CopyBestFitToInitValues(*fparm1);
       //set the mixing23 seed to be the mirror point of the initial best fit
       double sinsqtheta23_bestfit = fparm1->BestFitValue(conv::kPhypMixing23);
       double sinsqtheta13_bestfit = fparm1->BestFitValue(conv::kPhypMixing13);
       double centre = 1 / (2 * (1 - sinsqtheta13_bestfit));
       double reseed = sinsqtheta23_bestfit;
       if     (reseed < centre)
	 reseed = centre + TMath::Abs(centre - sinsqtheta23_bestfit);
       else if(reseed > centre)
	 reseed = centre - TMath::Abs(centre - sinsqtheta23_bestfit);
       fparm2->SetInitValue(conv::kPhypMixing23, reseed);
       // Get fitter
       FitInputs * fit_inputs2 = new FitInputs(skobs,skexp,fparm2);
       fit_inputs2->SetChisqDefinition(true); //5-bin chisq
       FitterABC * fitter2 = util::GetFitter(fit_inputs2,fit_option);
       fparm2->ResetCurrToInitValues(); // reset osc params
       fitter2->DoTheFit();
       toy_chisq2 = fit_inputs2->BestFitChisquare();
       // Save spectra (data and best-fit prediction)
       if(store_toy_expts) {
	 skobs->Write(dat_array2, "obs");
       }
       if(store_bestfit_spectrum) {
	 skexp->Write(pdf_array2, "pdf");
       }
     }

     // Print-out best-fit params, and set the "best of the two fits" branches correctly
     cout << "- Best-fit oscillation params (below) yield chisq = ";
     if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsqTheta) {
       if(toy_chisq2 < toy_chisq) {
	 cout << toy_chisq2 << " from the second seed fit" << endl;
	 if(toy_chisq2 > bestfit_chisq)
	   dataset_counter++;
	 hgof->Fill(toy_chisq2);
	 fparm2->ShowOnlyBestFitValuesInNextPrint();
	 cout << *fparm2;
	 br_fparm->SetAddress(&fparm2);
	 br_obs_array->SetAddress(&dat_array2);
	 br_bf_array->SetAddress(&pdf_array2);
       }
       else {
	 cout << toy_chisq << " from the first seed fit" << endl;
	 if(toy_chisq > bestfit_chisq)
	   dataset_counter++;
	 hgof->Fill(toy_chisq);
	 fparm1->ShowOnlyBestFitValuesInNextPrint();
	 cout << *fparm1;
	 br_fparm->SetAddress(&fparm1);
	 br_obs_array->SetAddress(&dat_array1);
	 br_bf_array->SetAddress(&pdf_array1);
       }
     }
     else {
       cout << toy_chisq << endl;
       if(toy_chisq > bestfit_chisq)
	 dataset_counter++;
       hgof->Fill(toy_chisq);
       fparm1->ShowOnlyBestFitValuesInNextPrint();
       cout << *fparm1;
       br_fparm->SetAddress(&fparm1);
       br_obs_array->SetAddress(&dat_array1);
       br_bf_array->SetAddress(&pdf_array1);
     }

     // Add all results in the tree
     fitnt->Fill();

     // Clean-up
     delete skobs;
     delete fit_inputs;
     delete fitter;
     if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsqTheta) {
       delete fit_inputs2;
       delete fitter2;
     }

  }//itoy

  //get the dataset probability
  double dataset_prob = double(dataset_counter) / double(ntoys);
  cout << "*** Data set probability is " << dataset_prob << endl;
  cout << "*** Data chisq is " << bestfit_chisq << endl;
  cout << fparm_input->GetBestFitChisq() << endl;

  // Text file, to save the important numbers to
  ofstream outfile (Form("gof-%s-expt_%d.txt", tag, iexperiment), ios_base::trunc);
  TString outstring(Form("%d\t%f\t%f\n", iexperiment, dataset_prob, bestfit_chisq));
  outfile.write (outstring.Data(), outstring.Sizeof());
  outfile.close();

  //store the important numbers
  TVectorD v(2);
  v[0] = dataset_prob;
  v[1] = bestfit_chisq;
  TLine l(bestfit_chisq, 0, bestfit_chisq, hgof->GetMaximum()*1.05);
  l.SetLineColor(kRed);
  TPaveText pt(15,hgof->GetMaximum()*0.85,22,hgof->GetMaximum()*1.00,"tr");
  pt.AddText(Form("#chi^{2}_{data} = %.2f", bestfit_chisq));
  pt.AddText(Form("p-value = %.2f", dataset_prob));
  pt.SetShadowColor(0);
  pt.SetFillColor(kWhite);
  //pt.SetFillStyle(0);

  // Store ntuple with best-fit results
  f.cd();
  fitnt->Write();
  hgof->Write();
  l.Write("line");
  v.Write("goodness_params");
  pt.Write("label");
  f.Close();

  //draw the histogram
  TCanvas * c = new TCanvas();
  hgof->Draw();
  hgof->Rebin(10);
  hgof->GetXaxis()->SetRangeUser(0,25);
  l.Draw();
  pt.Draw();
  c->SaveAs(Form("chisq_%d.pdf",iexperiment));
}
