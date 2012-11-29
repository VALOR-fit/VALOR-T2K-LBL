//
// ------------------------------------------------------------------------------------------------
// Use the minimum of a chisq map to seed the fitter, to find the true minimum of the data
// Output is the values and pulls of the fitted parameters, and a deltachisq map
//
// Tom Dealtry
// ------------------------------------------------------------------------------------------------
//

void fit_xml_seeded_from_chisq_map
(
  const char * xmlfile,                                            /* input XML file with toy MC experiments */
  const int    dataset_id,                                         /* which dataset to use*/
  const char * chisqfile,                                          /* root file where the chisq map is */
  const int    nfits,                                              /* number of fits to do (each with different seeds of the systematic parameters */
  const char * config_file         = "analysis_configuration.xml", /* analysis configuration file */
  const char * fit_param_list      = "nuosc3f_numudisapp_fix4+nuisanceselection",    /* fit parameter list for the given physics hypothesis (see run_fitter.C) */
  const char * fit_option          = "minuit"                      /* fitter option (see run_fitter.C) */
)
{
  using namespace std;
  using namespace t2kana;

  // Analysis configuration
  const char * config = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);

  // Global options
  Opt::Instance()->gOscMixingParamChoice = conv::kNuMixSinsq2Theta;
  Opt::Instance()->gMassHierarchyAssumed = conv::kNmhNormal;

  // Print settings
  cout << "\n[Job inputs]" << endl;
  cout << "- Reading analysis configuration from: " << config << endl;
  cout << "- Toy-MC" << endl;
  cout << "   - Reading experiment " << dataset_id << " from: " << xmlfile << endl;
  cout << "   - Reading chisq map from: " << chisqfile << endl;
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

  // Open the file to read the chisq map from it
  TFile inputfile(chisqfile,"READ");
  TH2D *hchisq = 0;
  inputfile.GetObject(Form("chisqwithsyst_%d",dataset_id),hchisq);
  TH2D & hdchisq = *hchisq;
  
  // Find the minimum, and convert into mixing23 & dm23sq. This will give the seeds for the osc parameters
  int globalbin = hchisq->GetMinimumBin();
  int xbin, ybin, zbin;
  hchisq->GetBinXYZ(globalbin, xbin, ybin, zbin);
  double input_mixing23 = hchisq->GetXaxis()->GetBinCenter(xbin);
  double input_dm23sq   = hchisq->GetYaxis()->GetBinCenter(ybin);

  // Seed the oscillation parameters
  fparm->SetInitValue(conv::kPhypMixing23, input_mixing23);
  fparm->SetInitValue(conv::kPhypDm23sq,   input_dm23sq);

  // random numbers
  TRandom3 rand(0);

  // Get the BANFF errors & bestfits
  double err_banff_maqe   = skexp->GetBANFF()->OneSigmaError(conv::kNuisBANFFMAQE);
  double err_banff_ccqeE0 = skexp->GetBANFF()->OneSigmaError(conv::kNuisBANFFCCQEE0);
  double bf_banff_maqe    = skexp->GetBANFF()->BestFitValue(conv::kNuisBANFFMAQE);
  double bf_banff_ccqeE0  = skexp->GetBANFF()->BestFitValue(conv::kNuisBANFFCCQEE0);

  // Variables to save bestfits & errors in
  double bestfit_chisq = 9999999;
  double bestfit_mixing23 = -1;
  double bestfit_dm23sq   = -1;
  double bestfit_sf       = -1;
  double bestfit_piless   = -1;
  double bestfit_binde    = -1;
  double bestfit_ccnqeeff = -1;
  double bestfit_nceff    = -1;
  //double bestfit_escale   = -1;
  double bestfit_maqe     = -1;
  double bestfit_ccqeE0   = -1;
  double err_mixing23 = -1;
  double err_dm23sq   = -1;
  double err_sf       = -1;
  double err_piless   = -1;
  double err_binde    = -1;
  double err_ccnqeeff = -1;
  double err_nceff    = -1;
  //double err_escale   = -1;
  double err_maqe     = -1;
  double err_ccqeE0   = -1;

  // Get superk spectrum observation (toy experiment)
  SuperkObservation * skobs = new SuperkObservation;
  skobs->Read(xmlfile,dataset_id);

  //Get spectra
  TH1D bestfit_spectrum;
  TH1D observation_spectrum(*(skobs->ErecoSpectrum(conv::kSKsFCFV1RingMuLike)));

  for(int ifit = 0; ifit < nfits; ifit++) {

    // Randomise seeds
    fparm->SetInitValue(conv::kNuisXSecSpecFunc,         rand.Uniform(0,1));
    fparm->SetInitValue(conv::kNuisXSecPiLessDelta,      rand.Gaus(0,1));
    fparm->SetInitValue(conv::kNuisXSecBindE,            rand.Gaus(0,1));
    fparm->SetInitValue(conv::kNuisSuperkEff1RmuCCnonQE, rand.Gaus(0,1));
    fparm->SetInitValue(conv::kNuisSuperkEff1RmuNC,      rand.Gaus(0,1));
    //fparm->SetInitValue(conv::kNuisSuperkEnergyScale,    rand.Gaus(0,1));
    fparm->SetInitValue(conv::kNuisBANFFMAQE,            rand.Gaus(bf_banff_maqe,   err_banff_maqe));
    fparm->SetInitValue(conv::kNuisBANFFCCQEE0,          rand.Gaus(bf_banff_ccqeE0, err_banff_ccqeE0));

    // Re-seed fitter (reset current values to specified init values)
    fparm->ResetCurrToInitValues();

    // Get fitter
    FitInputs * fit_inputs = new FitInputs(skobs,skexp,fparm);
    FitterABC * fitter = new Superk1RmuErecoFitter(fit_inputs,fit_option);

    // Do the fit
    fitter->DoTheFit();
    int fitstatus = fitter->FitStatusCode();

    // Get best-fit params (if its a better minimum)
    if(fit_inputs->BestFitChisquare() < bestfit_chisq) {
      bestfit_chisq = fit_inputs->BestFitChisquare();
      bestfit_mixing23 = fparm->BestFitValue(t2kana::conv::kPhypMixing23);
      bestfit_dm23sq   = fparm->BestFitValue(t2kana::conv::kPhypDm23sq);
      bestfit_sf       = fparm->BestFitValue(conv::kNuisXSecSpecFunc);
      bestfit_piless   = fparm->BestFitValue(conv::kNuisXSecPiLessDelta);
      bestfit_binde    = fparm->BestFitValue(conv::kNuisXSecBindE);
      bestfit_ccnqeeff = fparm->BestFitValue(conv::kNuisSuperkEff1RmuCCnonQE);
      bestfit_nceff    = fparm->BestFitValue(conv::kNuisSuperkEff1RmuNC);
      //bestfit_escale   = fparm->BestFitValue(conv::kNuisSuperkEnergyScale);
      bestfit_maqe     = fparm->BestFitValue(conv::kNuisBANFFMAQE);
      bestfit_ccqeE0   = fparm->BestFitValue(conv::kNuisBANFFCCQEE0);
      err_mixing23 = fparm->MigradError(t2kana::conv::kPhypMixing23);
      err_dm23sq   = fparm->MigradError(t2kana::conv::kPhypDm23sq);
      err_sf       = fparm->MigradError(conv::kNuisXSecSpecFunc);
      err_piless   = fparm->MigradError(conv::kNuisXSecPiLessDelta);
      err_binde    = fparm->MigradError(conv::kNuisXSecBindE);
      err_ccnqeeff = fparm->MigradError(conv::kNuisSuperkEff1RmuCCnonQE);
      err_nceff    = fparm->MigradError(conv::kNuisSuperkEff1RmuNC);
      //err_escale   = fparm->MigradError(conv::kNuisSuperkEnergyScale);
      err_maqe     = fparm->MigradError(conv::kNuisBANFFMAQE);
      err_ccqeE0   = fparm->MigradError(conv::kNuisBANFFCCQEE0);
      // Print-out best-fit params
      fparm->ShowOnlyBestFitValuesInNextPrint();
      cout << "- Best-fit oscillation params (below) yield chisq = " << bestfit_chisq 
	   << *fparm;
      //Get the bestfit spectrum
      bestfit_spectrum =  *(skexp->CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike));
    }
    delete fit_inputs;
    delete fitter;
  }//ifit
  delete skobs;

  //Output file
  //inputfile.Close();
  TFile outputfile(Form("fitresults_%d.root",dataset_id),"RECREATE");

  // Now we have the best fits, lets plot them
  for(int ix = 1; ix <= hdchisq.GetNbinsX(); ix++) {
    for(int iy = 1; iy <= hdchisq.GetNbinsY(); iy++) {
      hdchisq.SetBinContent(ix, iy, hdchisq.GetBinContent(ix, iy) - bestfit_chisq);
    }//iy
  }//ix

  gStyle->SetOptStat(0);

  TCanvas *c0 = new TCanvas();
  hdchisq.Draw("colz");
  c0->SaveAs(Form("dchisq_%d.pdf",dataset_id));
  double constchisq[2] = {2.30, 4.61};
  hdchisq.SetContour(2, constchisq);
  hdchisq.SetDirectory(0);
  hdchisq.Write("hdchisq");
  TCanvas *c1 = new TCanvas();
  hdchisq.SetTitle("");
  hdchisq.Draw("cont1");
  c1->SaveAs(Form("contour_%d.pdf",dataset_id));

  // Also plot the pulls of the systematics
  TH1D * hsystpull = new TH1D("hsystpull","hsystpull",7,0,7);
  hsystpull->GetXaxis()->SetBinLabel(1, "SF");
  hsystpull->GetXaxis()->SetBinLabel(2, "Piless-Delta");
  hsystpull->GetXaxis()->SetBinLabel(3, "Binding-E");
  hsystpull->GetXaxis()->SetBinLabel(4, "CCnQE-eff");
  hsystpull->GetXaxis()->SetBinLabel(5, "NC-eff");
  hsystpull->GetXaxis()->SetBinLabel(6, "MAQE");
  hsystpull->GetXaxis()->SetBinLabel(7, "CCQE-E0");
  //hsystpull->GetXaxis()->SetBinLabel(8,"E-scale");
  hsystpull->SetBinContent(1, bestfit_sf          / err_sf);
  hsystpull->SetBinContent(2, bestfit_piless      / err_piless);
  hsystpull->SetBinContent(3, bestfit_binde       / err_binde);
  hsystpull->SetBinContent(4, bestfit_ccnqeeff    / err_ccnqeeff);
  hsystpull->SetBinContent(5, bestfit_nceff       / err_nceff);
  hsystpull->SetBinContent(6, (bestfit_maqe   - bf_banff_maqe  )  / err_maqe);
  hsystpull->SetBinContent(7, (bestfit_ccqeE0 - bf_banff_ccqeE0)  / err_maqe);
  //hsystpull->SetBinContent(8 ,bestfit_escale      / err_escale);
  TCanvas *c2 = new TCanvas();
  hsystpull->Draw();
  c2->SaveAs(Form("pulls_%d.pdf",dataset_id));
  hsystpull->Write();
  
  //And plot the best fit spectra with the observation
  TCanvas *c3 = new TCanvas();
  bestfit_spectrum.GetXaxis()->SetRangeUser(0,9.99);
  observation_spectrum.GetXaxis()->SetRangeUser(0,9.99);
  observation_spectrum.Draw("E1");
  bestfit_spectrum.Draw("same");
  c3->SaveAs(Form("spectrum_%d.pdf",dataset_id));
  bestfit_spectrum.Write("bf_spectrum");
  observation_spectrum.Write("obs_spectrum");

  //And not forgetting the value & error of the mixing angles
  TVectorD v(4);
  v[0] = bestfit_mixing23;
  v[1] = err_mixing23;
  v[2] = bestfit_dm23sq;
  v[3] = err_dm23sq;
  cout << endl << "Bestfit of oscillation parameters yield chisq: " << bestfit_chisq << endl
       << "\tMixing23: " << bestfit_mixing23 << " +- " << err_mixing23 << endl
       << "\tDm32sq:   "   << bestfit_dm23sq << " +- "   << err_dm23sq << endl << endl;
  v.Write("oscparams");

  //inputfile.Close();
  outputfile.Close();
}
