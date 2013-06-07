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
  const char * tag,                                                /* */
  const int    run,                                                /* */
  const char * fitresults_file,                                    /* input file with results from fit_xml_seeded_from_chisq_map.C */
  const char * config_file         = "analysis_configuration.xml", /* analysis configuration file */
  const char * fit_param_list      = "nuosc3f_numudisapp_fix4+nuisanceselection",    /* fit parameter list for the given physics hypothesis (see run_fitter.C) */
  const char * fit_option          = "minuit",                     /* fitter option (see run_fitter.C) */
  const int    ntoys               = 10                          /* number of toys to create to put into the histogram */
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
  cout << "- Job: " << tag << "/" << run << endl;
  cout << "- Reading analysis configuration from: " << config << endl;
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

  // random numbers
  TRandom3 rand(0);

  // Get the BANFF errors & bestfits
  double err_banff_maqe   = skexp->GetBANFF()->OneSigmaError(conv::kNuisBANFFMAQE);
  double err_banff_ccqeE0 = skexp->GetBANFF()->OneSigmaError(conv::kNuisBANFFCCQEE0);
  double bf_banff_maqe    = skexp->GetBANFF()->BestFitValue(conv::kNuisBANFFMAQE);
  double bf_banff_ccqeE0  = skexp->GetBANFF()->BestFitValue(conv::kNuisBANFFCCQEE0);

  // Get the bestfits
  TFile fin(fitresults_file);
  if(fin.IsZombie()) return;
  TVectorD * bfvecp = 0;
  fin.GetObject("oscparams",bfvecp);
  const TVectorD & bfvec = *bfvecp;
  double bestfit_chisq    = bfvec[5];
  double bestfit_mixing23 = bfvec[0];
  double bestfit_dm23sq   = bfvec[2];
  fin.Close();

  //now we have the bestfit, we'll fit lots of toys (seeding them a few times each) and find the chisq of each
  fparm->SetInitValue(conv::kPhypMixing23, bestfit_mixing23);
  fparm->SetInitValue(conv::kPhypDm23sq,   bestfit_dm23sq);

  //book histogram
  TH1D * hchisq = new TH1D(Form("chisq_%d",run), Form("chisq_%d",run), 100000, 0, 1000);

  //counter for how many chisqs are greater than the data chisq
  int dataset_counter = 0;
  for(int itoy = 0; itoy < ntoys; itoy++) {

    //make some mock data, with statistical & all systematics fluctuations
    skexp->TweakAllSystematics();
    SuperkObservation * skobs = new SuperkObservation;
    skobs->CreateMockData(*skexp, PhysHypothesis(*fparm), true);

    //now do the fit, using the thrown values of the systematics as seeds
    fparm->SetInitValue(conv::kNuisXSecSpecFunc,         skexp->GetNuXSec()->CurrentValue(conv::kNuisXSecSpecFunc));
    fparm->SetInitValue(conv::kNuisXSecPiLessDelta,      skexp->GetNuXSec()->CurrentValue(conv::kNuisXSecPiLessDelta));
    fparm->SetInitValue(conv::kNuisXSecBindE,            skexp->GetNuXSec()->CurrentValue(conv::kNuisXSecBindE));
    fparm->SetInitValue(conv::kNuisSuperkEff1RmuCCnonQE, skexp->GetSuperkDetector()->CurrentValue(conv::kNuisSuperkEff1RmuCCnonQE));
    fparm->SetInitValue(conv::kNuisSuperkEff1RmuNC,      skexp->GetSuperkDetector()->CurrentValue(conv::kNuisSuperkEff1RmuNC));
    fparm->SetInitValue(conv::kNuisSuperkEnergyScale,    skexp->GetSuperkDetector()->CurrentValue(conv::kNuisSuperkEnergyScale));
    fparm->SetInitValue(conv::kNuisBANFFMAQE,            skexp->GetBANFF()->CurrentValue(conv::kNuisBANFFMAQE));
    fparm->SetInitValue(conv::kNuisBANFFCCQEE0,          skexp->GetBANFF()->CurrentValue(conv::kNuisBANFFCCQEE0));
     
    // Re-seed fitter (reset current values to specified init values)
    fparm->ResetCurrToInitValues();

    // Get fitter
    FitInputs * fit_inputs = new FitInputs(skobs,skexp,fparm);
    FitterABC * fitter = new Superk1RmuErecoFitter(fit_inputs,fit_option);
       
    //Make sure we're getting the correct chisq (i.e. the 5 bin one)
    fit_inputs->SetChisqDefinition(true);

    // Do the fit
    fitter->DoTheFit();
    int fitstatus = fitter->FitStatusCode();

    double toy_chisq = fit_inputs->BestFitChisquare();
       
    //cleanup
    delete fit_inputs;
    delete fitter;
    delete skobs;

    //fill the histogram & increment counter (if required)
    hchisq->Fill(toy_chisq);
    if(toy_chisq > bestfit_chisq)
      dataset_counter++;
  }//itoy

  //get the dataset probability
  double dataset_prob = double(dataset_counter) / double(ntoys);
  cout << "*** Data set probability is " << dataset_prob << endl;
  cout << "*** Data chisq is " << bestfit_chisq << endl;

  // File, to save things to
  TFile    f       (Form("goodness-of-fit-%s-run_%d.root", tag, run), "recreate");
  ofstream outfile (Form("goodness-of-fit-%s-run_%d.txt",  tag, run), ios_base::trunc);

  //save the histogram
  hchisq->Write();
  TCanvas *c = new TCanvas();
  hchisq->Draw();
  TLine l(bestfit_chisq, 0, bestfit_chisq, 200);
  l.SetLineColor(kRed);
  l.Draw();
  l.Write(Form("line_%d",run));
  c->SaveAs(Form("chisq_%d.pdf",run));
  delete hchisq;
  delete c;
  TVectorD v(2);
  v[0] = dataset_prob;
  v[1] = bestfit_chisq;
  v.Write("goodness_params");

  //save the important numbers
  TString outstring(Form("%d\t%f\t%f\n",run, dataset_prob, bestfit_chisq));
  outfile.write (outstring.Data(), outstring.Sizeof());

  //close
  outfile.close();
  f.Close();

}
