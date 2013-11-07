//
// ------------------------------------------------------------------------------------------------
// Run a job with some of the toy expts fits needed to calculate the \Delta\chi^{2} values
// for constructing a Feldman-Cousins contour 
//
// Costas Andreopoulos
//
// ChangeLog
//
// 16/8/2013 Nick
// Modified to do fits for 1D confidence intervals
// Oscillation parameter min and max values are read from a 2D delta chisq surface
// since make_feldman_cousins_coverage_1D.C uses numbers of jobs in each bin that 
// are proportional to exp(-chisq/2).
// 5 fits are made of each experiment:
// 1. 1+N fit where the 1 is sinsqtheta23 and |dm32sq| is fixed to the true input value
// 2. 1+N fit where the 1 is sinsqtheta23 with initial value at the mirror point in the other 
//    octant and |dm32sq| is fixed to the true input value
// 3. 1+N fit where the 1 is |dm32sq| and sinsqtheta23 is fixed to the true input value
// 4. 2+N fit to find best-fit chisq with sinsqtheta23 and |dm32sq| as free parameters
// 5. 2+N fit to find best-fit chisq with sinsqtheta23 and |dm32sq| as free parameters, with 
//    initial value of sinsqtheta23 being at the mirror point in the other octant 
// Added extra branches to output ROOT tree for the chisq and status of these fits
// 
// ------------------------------------------------------------------------------------------------
//

void feldman_cousins_fitting_job_1D
(
  const char * tag                     = "fc",                                /* tag for output file */
  int          run                     = 0,                                   /* run number */
  int          nexperiments            = 3,                                   /* number of toy MC experiment fits */
  double       dm23sq_input            = -999.9,                                 /* generate toy MC with this range of Dm^2_{23} */
  double       mixing23_input          = 0.375,                                 /* generate toy MC with this range of 23 mixing (expressed in terms of either sinsq2theta, sinsqtheta, ...) */
  const char * config_file             = "2013a_numu_run1+2+3+4_flux11bv3.2_neut13a_banff13v5_sk+fsi.xml",        /* analysis configuration file */
  bool         stat_fluctuations       = true,                                /* include statistical fluctuation in toy MC expts? */
  bool         tweak_systematics       = true,                                /* generate toy MC with randomized systematics? */
  bool         tweak_non23oscparams    = true,                                /* generate toy MC with randomized dmsq12,theta12,theta13,dcp / only in 3-flavour fmwk */
  const char * fit_param_list          = "nuosc3f_numudisapp_nuis4+nuisance",  /* fit parameter list for the given physics hypothesis (see FitParamList.h) */
  const char * fit_option              = "minuit",                            /* fitter option */
  bool         store_toy_expts         = false,                                /* disable this if you don't want to store the toy MC expts     */
  bool         store_bestfit_spectrum  = false,                                /* disable this if you don't want to store the best-fit spectra */
  const char * dmsqscheme              = "abs23",                             /* abs23==kNuDmsq12andAbs23 fogli==kNuDmsqFogliLisi */
  const char * hierarchy               = "normal",                            /* normal or inverted */
  const char * mixparamchoice          = "sinsqtheta",                        /* sinsq2theta or sinsqtheta */
  const char * mixparamoctant          = "one",                                /* one==kNuOctantOne two==kNuOctantTwo */
  const char * dchisqsurf_file         = "/home/ppd/tdealtry/RAL_repoistory/t2k/t2kana/nuosc/trunk/2013results/DATA/FULLRUN1-4FIXSPLINE/run1-4_nuis4_sintheta_normalfine/combination/final.root"
  //bool         apply_2nu_boundary    = true                                 /* if using the 2nu flavour osc prob, apply the physical boundary? */
)
{
  using namespace std;
  using namespace t2kana;
  using namespace t2kana::conv;

  // Analysis configuration
  TString config = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);

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

  TFile *dchisqfile = TFile::Open(dchisqsurf_file);
  TH2D *dchisq2D = (TH2D*)dchisqfile->Get("chisq_surface");

  // dm23sq, mixing23 range and step
  double mixing23_min   = dchisq2D->GetXaxis()->GetBinCenter(1);
  double mixing23_max   = dchisq2D->GetXaxis()->GetBinCenter(dchisq2D->GetNbinsX());
  double mixing23_step  = dchisq2D->GetXaxis()->GetBinCenter(2) - dchisq2D->GetXaxis()->GetBinCenter(1);
  double dm23sq_min     = dchisq2D->GetYaxis()->GetBinCenter(1);
  double dm23sq_max     = dchisq2D->GetYaxis()->GetBinCenter(dchisq2D->GetNbinsY());
  double dm23sq_step    = dchisq2D->GetYaxis()->GetBinCenter(2) - dchisq2D->GetYaxis()->GetBinCenter(1);
  // make sure that min and max are at the bin centre and not at the low/high edge of the first/last bin
  //add 0.001 to avoid rounding error when finding (int)
  int nbins_mixing23 = 1 + (int)(((mixing23_max - mixing23_min) / mixing23_step) + 0.001);
  int nbins_dm23sq   = 1 + (int)(((dm23sq_max   - dm23sq_min  ) / dm23sq_step  ) + 0.001);

  const int mixing23_nbins = nbins_mixing23;
  const int dm23sq_nbins = nbins_dm23sq;

  // param space
  TH2D * osc_param_space = new TH2D("osc_param_space", "",
       mixing23_nbins, 
       mixing23_min - 0.5*mixing23_step, 
       mixing23_max + 0.5*mixing23_step,
       dm23sq_nbins, 
       dm23sq_min - 0.5*dm23sq_step, 
       dm23sq_max + 0.5*dm23sq_step);

  int input_mixing23_bin, input_dm23sq_bin;
  double input_mixing23, input_dm23sq;

  if(mixing23_input > 0 && dm23sq_input < 0)
    {
    input_mixing23 = mixing23_input;
    input_mixing23_bin = dchisq2D->GetXaxis()->FindBin(input_mixing23);
    cout<<"Input value of mixing23 = "<<input_mixing23<<" in bin "<<input_mixing23_bin<<endl;

  double totalLikelihoodDm23sq = 0.0;

  for(int dm23sq_bin=1; dm23sq_bin<=dm23sq_nbins; dm23sq_bin++)
    totalLikelihoodDm23sq += exp(-0.5 * dchisq2D->GetBinContent(input_mixing23_bin, dm23sq_bin));

  double cumulNormLikelihoodDm23sq[dm23sq_nbins+1];

  cumulNormLikelihoodDm23sq[0] = 0.0;

  for(int dm23sq_bin=1; dm23sq_bin<=dm23sq_nbins; dm23sq_bin++)  
    {
    cumulNormLikelihoodDm23sq[dm23sq_bin] = cumulNormLikelihoodDm23sq[dm23sq_bin-1] 
                                           + (exp(-0.5 * dchisq2D->GetBinContent(input_mixing23_bin, dm23sq_bin)) / totalLikelihoodDm23sq);
    //prevent cumulative normalised likelihood being > 1.0
    cumulNormLikelihoodDm23sq[dm23sq_bin] = TMath::Min(1.0, cumulNormLikelihoodDm23sq[dm23sq_bin]);
    }

  cumulNormLikelihoodDm23sq[dm23sq_nbins] = 1.0;

  cout<<endl;
  cout<<"This table shows |dm32sq| bin number, value of |dm32sq| at centre of bin, likelihood, normalised likelihood and cumulative normalised likelihood from the 2D chisq surface for each bin."
      <<endl;
  cout<<endl;

  for(int dm23sq_bin=1; dm23sq_bin<=dm23sq_nbins; dm23sq_bin++)
    {
    cout<<dm23sq_bin<<"    "<<dchisq2D->GetYaxis()->GetBinCenter(dm23sq_bin)
	<<"    "<<exp(-0.5 * dchisq2D->GetBinContent(input_mixing23_bin, dm23sq_bin))
        <<"    "<<"    "<<(exp(-0.5 * dchisq2D->GetBinContent(input_mixing23_bin, dm23sq_bin)) / totalLikelihoodDm23sq)
	<<"    "<<cumulNormLikelihoodDm23sq[dm23sq_bin]<<endl;
    }
    }//end of if(mixing23_input > 0 && dm23sq_input < 0)
  else if(mixing23_input < 0 && dm23sq_input > 0)
    {
      input_dm23sq = dm23sq_input;
      input_dm23sq_bin = dchisq2D->GetYaxis()->FindBin(input_dm23sq);
      cout<<"Input value of dm23sq = "<<input_dm23sq<<" in bin "<<input_dm23sq_bin<<endl;

      double totalLikelihoodMixing23 = 0.0;

      for(int mixing23_bin=1; mixing23_bin<=mixing23_nbins; mixing23_bin++)
	totalLikelihoodMixing23 += exp(-0.5 * dchisq2D->GetBinContent(mixing23_bin, input_dm23sq_bin));

      double cumulNormLikelihoodMixing23[mixing23_nbins+1];

      cumulNormLikelihoodMixing23[0] = 0.0;

      for(int mixing23_bin=1; mixing23_bin<=mixing23_nbins; mixing23_bin++)
	{
    cumulNormLikelihoodMixing23[mixing23_bin] = cumulNormLikelihoodMixing23[mixing23_bin-1]
      + (exp(-0.5 * dchisq2D->GetBinContent(mixing23_bin, input_dm23sq_bin)) / totalLikelihoodMixing23);
    //prevent cumulative normalised likelihood being > 1.0
        cumulNormLikelihoodMixing23[mixing23_bin] = TMath::Min(1.0, cumulNormLikelihoodMixing23[mixing23_bin]);
	}

      cumulNormLikelihoodMixing23[mixing23_nbins] = 1.0;

      cout<<endl;
      cout<<"This table shows mixing23 bin number, value of mixing23 at centre of bin, likelihood, normalised likelihood and cumulative normalised likelihood from the 2D chisq surface for each bin."
          <<endl;
      cout<<endl;

      for(int mixing23_bin=1; mixing23_bin<=mixing23_nbins; mixing23_bin++)
	{
	  cout<<mixing23_bin<<"    "<<dchisq2D->GetXaxis()->GetBinCenter(mixing23_bin)
	      <<"    "<<exp(-0.5 * dchisq2D->GetBinContent(mixing23_bin, input_dm23sq_bin))
	      <<"    "<<"    "<<(exp(-0.5 * dchisq2D->GetBinContent(mixing23_bin, input_dm23sq_bin)) / totalLikelihoodMixing23)
	      <<"    "<<cumulNormLikelihoodMixing23[mixing23_bin]<<endl;
	}
    }//end of else if(mixing23_input < 0 && dm23sq_input > 0)
  else
    {
      cerr<<"You must choose (mixing23_input > 0 and dm23sq_input < 0) or (mixing23_input < 0 and dm23sq_input > 0)"<<endl;
      exit(1);
    }

  // Print settings
  cout << "\n[Job inputs]" << endl;
  cout << "- Job: " << tag << "/" << run << endl;
  cout << "- Reading analysis configuration from: " << config << endl;
  cout << "- Toy-MC generation" << endl;
  cout << "   - Number of toy experiments: " << nexperiments << endl; 
  cout << "   - Randomize systematics : "      << ((tweak_systematics) ?      "YES " : "NO") << endl;
  cout << "   - Statistical fluctuations: "    << ((stat_fluctuations) ?      "YES " : "NO") << endl;
  cout << "   - Randomizing non-'23'-sector oscillation params: " << ((tweak_non23oscparams) ? "YES" : "NO") << endl;
  cout << "- Fit" << endl;
  cout << "   - Fit Parameter List: " << fit_param_list << endl;
  cout << "   - Fitter option: " << fit_option << endl;
  cout << "- Oscillation parameter space binning" << endl;
  cout << "   - " << util::AsString(conv::kPhypMixing23) << " : " 
       << osc_param_space->GetXaxis()->GetNbins() << " bins in (" 
       << osc_param_space->GetXaxis()->GetXmin() << ", "
       << osc_param_space->GetXaxis()->GetXmax() << ")" << endl;
  cout << "   - " << util::AsString(conv::kPhypDm23sq) << " : " 
       << osc_param_space->GetYaxis()->GetNbins() << " bins in (" 
       << osc_param_space->GetYaxis()->GetXmin() << ", "
       << osc_param_space->GetYaxis()->GetXmax() << ") eV^2/c^4" << endl;
  cout << "- Input bin range " << endl;
  cout << "   - " << util::AsString(conv::kPhypMixing23) 
    //<< " : [" << min_mixing23_bin << ", " << max_mixing23_bin << "]"
       << " : [" << osc_param_space->GetXaxis()->GetBinCenter(1)
       << ", "   << osc_param_space->GetXaxis()->GetBinCenter(mixing23_nbins)
       << "]" << endl;
  cout << "   - " << util::AsString(conv::kPhypDm23sq) 
    // << " : [" << min_dm23sq_bin   << ", " << max_dm23sq_bin   << "]"
       << " : [" << osc_param_space->GetYaxis()->GetBinCenter(1)
       << ", "   << osc_param_space->GetYaxis()->GetBinCenter(dm23sq_nbins)
       << "]" << endl;
  cout << "- Input (" << util::AsString(conv::kPhypMixing23) << ", " << util::AsString(conv::kPhypDm23sq) << ") bin list " << endl; 
  cout << "   - none" << endl;
  cout << "Doing a 2D FC scanning dm32sq and mixing23" << endl;
  cout << "- Outputs" << endl;
  cout << "     - Storing toy MC experiments:  " << ((store_toy_expts) ?        "YES " : "NO") << endl;
  cout << "     - Storing best-fit spectra: "    << ((store_bestfit_spectrum) ? "YES " : "NO") << endl;
  cout << endl;
  cout << *Opt::Instance() << endl;
  cout << *Rnd::Instance() << endl;
  
  // Fit parameter list
  FitParamList * fparm      = FitParamList::NamedParameterList(fit_param_list); //only used to initialise branch
  FitParamList * fparm1     = FitParamList::NamedParameterList(fit_param_list); // for 1+N fit of mixing23 using the first seed
  FitParamList * fparm2     = FitParamList::NamedParameterList(fit_param_list); // for 1+N fit of mixing23 using the mirror of best fit seed
  FitParamList * fparm3     = FitParamList::NamedParameterList(fit_param_list); // for 1+N fit of |dm32sq| 
  FitParamList * fparm4     = FitParamList::NamedParameterList(fit_param_list); // for 2+N fit using the first seed of mixing23
  FitParamList * fparm5     = FitParamList::NamedParameterList(fit_param_list); // for 2+N fit using the mirror of best fit seed of mixing23
  PhysHypothesis_t hyp = fparm1->PhysicsHypothesis();
  cout << "- The following physics hypothesis is used in the fit: " << util::AsString(hyp) << endl;

  // Get superk spectrum predictor & configure
  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read(config.Data());

  // Object for storing input oscillation and nuisance params in each toy-MC
  PhysParamList * iparm = new PhysParamList;

  // Arrays for storing all pieces of fit data and the best-fit pdf (and its components) for each toy-MC.
  TObjArray * dat_array = new TObjArray;
  TObjArray * pdf_array = new TObjArray; //only used to initialise branch
  dat_array->SetOwner(false);
  pdf_array->SetOwner(false);
  TObjArray * pdf_array1 = new TObjArray;
  pdf_array1->SetOwner(false);
  TObjArray * pdf_array2 = new TObjArray;
  pdf_array2->SetOwner(false);
  TObjArray * pdf_array3 = new TObjArray;
  pdf_array3->SetOwner(false);
  TObjArray * pdf_array4 = new TObjArray;
  pdf_array4->SetOwner(false);
  TObjArray * pdf_array5 = new TObjArray;
  pdf_array5->SetOwner(false);

  // Create file to store results
  TFile f(Form("%s-run_%d.root",tag, run), "recreate"); 

  //numbers to store results in
  int iexperiment;
  int nsk;
  int fittries_mixing23_1N, fittries_mixing23_1N_alt, fittries_dm23sq_1N, fittries_2N, fittries_2N_alt;
  double fittime_mixing23_1N, fittime_mixing23_1N_alt, fittime_dm23sq_1N, fittime_2N, fittime_2N_alt;
  int migradmin_mixing23_1N, migradmin_mixing23_1N_alt, migradmin_dm23sq_1N, migradmin_2N, migradmin_2N_alt;
  int migraderr_mixing23_1N, migraderr_mixing23_1N_alt, migraderr_dm23sq_1N, migraderr_2N, migraderr_2N_alt;
  int migradmin1, migradmin2, migradmin3, migradmin4, migradmin5;
  int migraderr1, migraderr2, migraderr3, migraderr4, migraderr5; 
  double bestfit_chisq_seed1, bestfit_chisq_seed2;
  double sinsqtheta23_bestfit, sinsqtheta13_bestfit, centre, reseed;
  double bestfit_mixing23_1N, bestfit_mixing23_1N_alt, bestfit_dm23sq_1N, bestfit_chisq_mixing23_1N, bestfit_chisq_mixing23_1N_alt,  bestfit_chisq_dm23sq_1N;
  double bestfit_mixing23_2N, bestfit_mixing23_2N_alt, bestfit_dm23sq_2N, bestfit_dm23sq_2N_alt, bestfit_chisq_2N, bestfit_chisq_2N_alt;


  // Create tree to store outputs
  TTree * fitnt = new TTree("fitnt","Results of VaLOR data and/or toy-MC experiment fits");
  /*
  fitnt->Branch("best_fit_params",           "t2kana::FitParamList",  &fparm,          1280000,0);
  fitnt->Branch("best_fit_params_fit1",      "t2kana::FitParamList",  &fparm1,         1280000,0);
  fitnt->Branch("best_fit_params_fit2",      "t2kana::FitParamList",  &fparm2,         1280000,0);
  fitnt->Branch("best_fit_params_fit3",      "t2kana::FitParamList",  &fparm3,         1280000,0);
  fitnt->Branch("best_fit_params_fit4",      "t2kana::FitParamList",  &fparm4,         1280000,0);
  fitnt->Branch("best_fit_params_fit5",      "t2kana::FitParamList",  &fparm5,         1280000,0);
  fitnt->Branch("input_params",              "t2kana::PhysParamList", &iparm,          1280000,0);
  */
  if(store_toy_expts)
    {
    fitnt->Branch("data",                      "TObjArray",             &dat_array,      1280000,0);
    }
  if(store_bestfit_spectrum)
    {
    fitnt->Branch("best_fit_pdf",              "TObjArray",             &pdf_array,      1280000,0);
    fitnt->Branch("best_fit_pdf_fit1",         "TObjArray",             &pdf_array1,     1280000,0);
    fitnt->Branch("best_fit_pdf_fit2",         "TObjArray",             &pdf_array2,     1280000,0);
    fitnt->Branch("best_fit_pdf_fit3",         "TObjArray",             &pdf_array3,     1280000,0);
    fitnt->Branch("best_fit_pdf_fit4",         "TObjArray",             &pdf_array4,     1280000,0);
    fitnt->Branch("best_fit_pdf_fit5",         "TObjArray",             &pdf_array5,     1280000,0);
    }
  
  fitnt->Branch("iexp",                      &iexperiment);
  fitnt->Branch("irun",                      &run);
  fitnt->Branch("nsk",                       &nsk);
  fitnt->Branch("input_mixing23_bin",        &input_mixing23_bin);
  fitnt->Branch("input_dm23sq_bin",          &input_dm23sq_bin);
  fitnt->Branch("input_mixing23",            &input_mixing23);
  fitnt->Branch("input_dm23sq",              &input_dm23sq);
  fitnt->Branch("fittries_mixing23_1N",      &fittries_mixing23_1N);
  fitnt->Branch("fittries_mixing23_1N_alt",  &fittries_mixing23_1N_alt);
  fitnt->Branch("fittries_dm23sq_1N",        &fittries_dm23sq_1N);
  fitnt->Branch("fittries_2N",               &fittries_2N);
  fitnt->Branch("fittries_2N_alt",           &fittries_2N_alt);
  fitnt->Branch("fittime_mixing23_1N",       &fittime_mixing23_1N);
  fitnt->Branch("fittime_mixing23_1N_alt",   &fittime_mixing23_1N_alt);
  fitnt->Branch("fittime_dm23sq_1N",         &fittime_dm23sq_1N);
  fitnt->Branch("fittime_2N",                &fittime_2N);
  fitnt->Branch("fittime_2N_alt",            &fittime_2N_alt);
  fitnt->Branch("bestfit_mixing23_1N",       &bestfit_mixing23_1N);
  fitnt->Branch("bestfit_mixing23_1N_alt",   &bestfit_mixing23_1N_alt);
  fitnt->Branch("bestfit_dm23sq_1N",         &bestfit_dm23sq_1N);
  fitnt->Branch("bestfit_chisq_mixing23_1N", &bestfit_chisq_mixing23_1N);
  fitnt->Branch("bestfit_chisq_mixing23_1N_alt", &bestfit_chisq_mixing23_1N_alt);
  fitnt->Branch("bestfit_chisq_dm23sq_1N",   &bestfit_chisq_dm23sq_1N);
  fitnt->Branch("migradmin_mixing23_1N",     &migradmin_mixing23_1N);
  fitnt->Branch("migradmin_mixing23_1N_alt", &migradmin_mixing23_1N_alt);
  fitnt->Branch("migradmin_dm23sq_1N",       &migradmin_dm23sq_1N);
  fitnt->Branch("migraderr_mixing23_1N",     &migraderr_mixing23_1N);
  fitnt->Branch("migraderr_mixing23_1N_alt", &migraderr_mixing23_1N_alt);
  fitnt->Branch("migraderr_dm23sq_1N",       &migraderr_dm23sq_1N);
  fitnt->Branch("bestfit_mixing23_2N",       &bestfit_mixing23_2N);
  fitnt->Branch("bestfit_mixing23_2N_alt",   &bestfit_mixing23_2N_alt);
  fitnt->Branch("bestfit_dm23sq_2N",         &bestfit_dm23sq_2N);
  fitnt->Branch("bestfit_dm23sq_2N_alt",     &bestfit_dm23sq_2N_alt);
  fitnt->Branch("bestfit_chisq_2N",          &bestfit_chisq_2N);
  fitnt->Branch("bestfit_chisq_2N_alt",      &bestfit_chisq_2N_alt);
  fitnt->Branch("migradmin_2N",              &migradmin_2N);
  fitnt->Branch("migradmin_2N_alt",          &migradmin_2N_alt);
  fitnt->Branch("migraderr_2N",              &migraderr_2N);
  fitnt->Branch("migraderr_2N_alt",          &migraderr_2N_alt);


  //TBranch * br_fparm     = fitnt->GetBranch("best_fit_params");
  if(store_bestfit_spectrum)
    TBranch * br_bf_array  = fitnt->GetBranch("best_fit_pdf");

  double uniformRand;

  TRandom3 rand(0);

  // Loop over toy-MC experiments
  for(iexperiment = 0; iexperiment < nexperiments; iexperiment++) {

     cout << "\n\n *** Pseudo-experiment = " << iexperiment << endl;

     //     
     // Random selection of dm23sq bin if mixing23 is fixed or mixing23 bin if |dm32sq| is fixed
     //     
     
     if(mixing23_input > 0 && dm23sq_input < 0)
       {
     input_dm23sq_bin   = -1;

     uniformRand = rand.Uniform(0.0, 1.0);

     for(int dm23sq_bin=1; dm23sq_bin<=dm23sq_nbins; dm23sq_bin++)
       {
	 if(uniformRand >= cumulNormLikelihoodDm23sq[dm23sq_bin-1] && uniformRand <= cumulNormLikelihoodDm23sq[dm23sq_bin])
	   input_dm23sq_bin = dm23sq_bin;
       }

     cout<<"Uniform random number thrown between 0 and 1 = "<<uniformRand<<" which means the input |dm32sq| bin is "<<input_dm23sq_bin<<endl;
       }//end of if(mixing23_input > 0 && dm23sq_input < 0)
     else if(mixing23_input < 0 && dm23sq_input > 0)
       {
	 input_mixing23_bin   = -1;

	 uniformRand = rand.Uniform(0.0, 1.0);

	 for(int mixing23_bin=1; mixing23_bin<=mixing23_nbins; mixing23_bin++)
	   {
	     if(uniformRand >= cumulNormLikelihoodMixing23[mixing23_bin-1] && uniformRand <= cumulNormLikelihoodMixing23[mixing23_bin])
	       input_mixing23_bin = mixing23_bin;
	   }

	 cout<<"Uniform random number thrown between 0 and 1 = "<<uniformRand<<" which means the input mixing23 bin is "<<input_mixing23_bin<<endl;
       }//end of else if(mixing23_input < 0 && dm23sq_input > 0)
     else
       {
         cerr<<"You must choose (mixing23_input > 0 and dm23sq_input < 0) or (mixing23_input < 0 and dm23sq_input > 0)"<<endl;
         exit(1);
       }

     // Sanity check on bin range
     input_mixing23_bin = TMath::Max(1,              input_mixing23_bin);
     input_mixing23_bin = TMath::Min(mixing23_nbins, input_mixing23_bin);
     input_dm23sq_bin   = TMath::Max(1,              input_dm23sq_bin  );
     input_dm23sq_bin   = TMath::Min(dm23sq_nbins,   input_dm23sq_bin  );

     cout << "- Toy-MC > Selected 23-sector params bin: " 
          << util::AsString(conv::kPhypMixing23) << " bin # = " << input_mixing23_bin << ", "
          << util::AsString(conv::kPhypDm23sq)   << " bin # = " << input_dm23sq_bin
          << endl;

     // Get dm23sq or mixing23 at the centre of the randomly selected bin
     if(mixing23_input > 0 && dm23sq_input < 0)
       input_dm23sq   = dchisq2D->GetYaxis()->GetBinCenter(input_dm23sq_bin);
     else if(mixing23_input < 0 && dm23sq_input > 0)
       input_mixing23 = dchisq2D->GetXaxis()->GetBinCenter(input_mixing23_bin);
     else
       {
	 cerr<<"You must choose (mixing23_input > 0 and dm23sq_input < 0) or (mixing23_input < 0 and dm23sq_input > 0)"<<endl;
	 exit(1); 
       }

     input_mixing23 = TMath::Max(0.,input_mixing23); // protection from rnd-err
     input_mixing23 = TMath::Min(1.,input_mixing23); // protection from rnd-err

     cout << "- Toy-MC > True 23-sector params: " 
          << util::AsString(conv::kPhypMixing23) << " = " << input_mixing23 << ", "
          << util::AsString(conv::kPhypDm23sq)   << " = " << input_dm23sq   << " eV^2"
          << endl;

     //
     // if we are working in 3-flavour framework, randomize mixing12,mixing13,deltacp and dm12sq
     // (these parameters are marginalised in the 3-flavour numu disappearance analysis)
     //

     double input_mixing12 = -1;
     double input_mixing13 = -1;
     double input_deltacp  = -1;
     double input_dm12sq   = -1;
     if(hyp == conv::kHypOsc3F) {
       if(tweak_non23oscparams) {
         input_mixing12 = util::RandomizeMixing12(); 
         input_mixing13 = util::RandomizeMixing13(); 
         input_deltacp  = util::RandomizeDeltaCP();
         input_dm12sq   = util::RandomizeDm12sq(); 
       }
       else {
	 if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsq2Theta) {
	   input_mixing12 = conv::kGlobFitSinsq2Theta12;
	   input_mixing13 = conv::kGlobFitSinsq2Theta13;
	 }
	 else {
	   input_mixing12 = conv::kGlobFitSinsqTheta12;
	   input_mixing13 = conv::kGlobFitSinsqTheta13;
	 }
         input_deltacp  = conv::kGlobFitDeltaCP;
         input_dm12sq   = conv::kGlobFitDm12sq;
       }
       cout << "- Toy-MC > True (randomized) non-23-sector params: " 
	    << util::AsString(conv::kPhypMixing12) << " = " << input_mixing12 << ", "
	    << util::AsString(conv::kPhypDm12sq)   << " = " << input_dm12sq   << " eV^2, "
	    << util::AsString(conv::kPhypMixing12) << " = " << input_mixing13 << ", "
	    << util::AsString(conv::kPhypDeltaCP)  << " = " << input_deltacp  << " rad "
	    << endl;
     }

     //
     // Set the generated oscillation parameter values (used to generate toy-MC expt)
     //
     fparm1->Update(conv::kPhypMixing23, input_mixing23);
     fparm1->Update(conv::kPhypDm23sq,   input_dm23sq);
     if(hyp == conv::kHypOsc3F) {
       fparm1->Update(conv::kPhypMixing12, input_mixing12);
       fparm1->Update(conv::kPhypMixing13, input_mixing13);
       fparm1->Update(conv::kPhypDeltaCP,  input_deltacp);
       fparm1->Update(conv::kPhypDm12sq,   input_dm12sq);
     }

     // Reset all output spectra
     if(store_toy_expts)
       dat_array     ->Clear();
     if(store_bestfit_spectrum)
       {
       pdf_array     ->Clear();
       pdf_array1    ->Clear();
       pdf_array2    ->Clear();
       pdf_array3    ->Clear();
       pdf_array4    ->Clear();
       pdf_array5    ->Clear();
       }

     // Reset input true params and add true oscillation params
     iparm->Clear();
     PhysParamList true_oscillation_params = fparm1->AsPhysParamList("current");
     iparm->Append(true_oscillation_params);

     //
     // get toy MC with/without systematics
     //
     SuperkObservation * skobs = new SuperkObservation;
     bool applied_negevents_condition = true;
     while(applied_negevents_condition) {
       if(tweak_systematics) {
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
	 applied_negevents_condition = skobs->CreateMockData(*skexp, PhysHypothesis(*fparm1), stat_fluctuations);
	 // reset all systematics so that the PDF used in the fitter is the nominal one
	 cout << "- Reseting all systematic params in pdf..." << endl;
	 skexp->ResetSystematicVariations();
	 //skexp->PrintSystematicVariations();
       }
       else {
	 //make sure systematics are reset
	 skexp->ResetSystematicVariations();
	 // generate mock data using the current PDF (without systematic variations) 
	 // with statistical fluctuations if requested (Poisson fluctuation in each bin)
	 cout << "- Generating toy-MC..." << endl;
	 applied_negevents_condition = skobs->CreateMockData(*skexp, PhysHypothesis(*fparm1), stat_fluctuations); 
       }
     }//applied_negevents_condition

     nsk = skobs->NumOfEvents(conv::kSKsFCFV1RingMuLike);

     cout<<endl;
     cout<<endl;
     cout<<endl;
     cout<<endl;
     cout<<"NSK = "<<skobs->NumOfEvents(conv::kSKsFCFV1RingMuLike)<<endl;
     cout<<endl;
     cout<<endl;
     cout<<endl;
     cout<<endl;

     //
     // Get the 1+N fit results for mixing23 (with |dm32sq| fixed to true input value)
     //

     fparm1->SetInitValue(conv::kPhypMixing23, input_mixing23);
     fparm1->SetInitValue(conv::kPhypDm23sq,   input_dm23sq);
     fparm1->FixParameter(conv::kPhypDm23sq);

     // Get fitter 
     FitInputs * fit_inputs = new FitInputs(skobs,skexp,fparm1);
     FitterABC * fitter = util::GetFitter(fit_inputs,fit_option);
     fparm1->ResetCurrToInitValues(); // reset osc params 

     fitter->DoTheFit();
     fparm1->ResetCurrToInitValues();
     migradmin1 = fitter->FitStatusCode("migrad_min");
     migraderr1 = fitter->FitStatusCode("migrad_err");
     cout<<"Migrad min and err codes for fit 1:   "<<migradmin1<<"    "<<migraderr1<<endl;
     bestfit_chisq_seed1  = fit_inputs->BestFitChisquare();
     bestfit_chisq_seed2 = 0;

     // Save spectra (data and best-fit prediction)
     if(store_toy_expts) {
        skobs->Write(dat_array, "obs");
     }
     if(store_bestfit_spectrum) {
        skexp->Write(pdf_array1, "pdf");
     }

     //Do the fit again at the mirror point
     migradmin2 = -1;
     if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsqTheta) {
       //get the correct initial values (best fit points from first fit)
       fparm2->CopyBestFitToInitValues(*fparm1);
       //set the mixing23 seed to be the mirror point of the initial best fit
       sinsqtheta23_bestfit = fparm1->BestFitValue(conv::kPhypMixing23);
       sinsqtheta13_bestfit = fparm1->BestFitValue(conv::kPhypMixing13);
       centre = 1 / (2 * (1 - sinsqtheta13_bestfit));
       reseed = sinsqtheta23_bestfit;
       if     (reseed < centre)
	 reseed = centre + TMath::Abs(centre - sinsqtheta23_bestfit);
       else if(reseed > centre)
	 reseed = centre - TMath::Abs(centre - sinsqtheta23_bestfit);
       fparm2->SetInitValue(conv::kPhypMixing23, reseed);
       fparm2->SetInitValue(conv::kPhypDm23sq,   input_dm23sq);
       fparm2->FixParameter(conv::kPhypDm23sq);

       // Get fitter
       FitInputs * fit_inputs2 = new FitInputs(skobs,skexp,fparm2);
       FitterABC * fitter2 = util::GetFitter(fit_inputs2,fit_option);
       fparm2->ResetCurrToInitValues(); // reset osc params
       fitter2->DoTheFit();
       bestfit_chisq_seed2 = fit_inputs2->BestFitChisquare();
       migradmin2 = fitter2->FitStatusCode("migrad_min");
       migraderr2 = fitter2->FitStatusCode("migrad_err");
       cout<<"Migrad min and err codes for fit 2:   "<<migradmin2<<"    "<<migraderr2<<endl;
       // Save spectra (best-fit prediction)
       if(store_bestfit_spectrum) {
	 skexp->Write(pdf_array2, "pdf");
       }
     }

     // Print-out best-fit params, and set the "best of the two fits" branches correctly
     cout << "- Best-fit mixing23 (below) yields chisq = ";
     if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsqTheta) {
       if(bestfit_chisq_seed2 < bestfit_chisq_seed1) {
	 cout << bestfit_chisq_seed2 << " from the second seed fit" << endl;
         fittries_mixing23_1N = fparm2->GetFitTries();
         fittries_mixing23_1N_alt = fparm1->GetFitTries();
         fittime_mixing23_1N = fparm2->GetFitTime();
	 fittime_mixing23_1N_alt = fparm1->GetFitTime();
	 bestfit_chisq_mixing23_1N    = bestfit_chisq_seed2;
	 bestfit_chisq_mixing23_1N_alt    = bestfit_chisq_seed1;
	 bestfit_mixing23_1N = fparm2->BestFitValue(conv::kPhypMixing23);
	 bestfit_mixing23_1N_alt = fparm1->BestFitValue(conv::kPhypMixing23);
	 migradmin_mixing23_1N        = migradmin2;
	 migradmin_mixing23_1N_alt        = migradmin1;
         migraderr_mixing23_1N        = migraderr2;
         migraderr_mixing23_1N_alt        = migraderr1;
	 fparm2->ShowOnlyBestFitValuesInNextPrint();
	 cout << *fparm2;
	 //br_fparm->SetAddress(&fparm2);
	 if(store_bestfit_spectrum)
	   br_bf_array->SetAddress(&pdf_array2);
       }
       else {
	 cout << bestfit_chisq_seed1 << " from the first seed fit" << endl;
         fittries_mixing23_1N = fparm1->GetFitTries();
	 fittries_mixing23_1N_alt = fparm2->GetFitTries();
	 fittime_mixing23_1N = fparm1->GetFitTime();
         fittime_mixing23_1N_alt = fparm2->GetFitTime();
         bestfit_chisq_mixing23_1N    = bestfit_chisq_seed1;
         bestfit_chisq_mixing23_1N_alt    = bestfit_chisq_seed2;
	 bestfit_mixing23_1N = fparm1->BestFitValue(conv::kPhypMixing23);
	 bestfit_mixing23_1N_alt = fparm2->BestFitValue(conv::kPhypMixing23);
	 migradmin_mixing23_1N        = migradmin1;
         migradmin_mixing23_1N_alt        = migradmin2;
         migraderr_mixing23_1N        = migraderr1;
         migraderr_mixing23_1N_alt        = migraderr2;
	 fparm1->ShowOnlyBestFitValuesInNextPrint();
	 cout << *fparm1;
	 //br_fparm->SetAddress(&fparm1);
	 if(store_bestfit_spectrum)
	   br_bf_array->SetAddress(&pdf_array1);
       }
     }
     else {
       cout << bestfit_chisq_seed1 << endl;
       fittries_mixing23_1N = fparm1->GetFitTries();
       fittime_mixing23_1N = fparm1->GetFitTime();
       bestfit_chisq_mixing23_1N    = bestfit_chisq_seed1;
       bestfit_mixing23_1N = fparm1->BestFitValue(conv::kPhypMixing23);
       migradmin_mixing23_1N        = migradmin1;
       migraderr_mixing23_1N        = migraderr1;
       fparm1->ShowOnlyBestFitValuesInNextPrint();
       cout << *fparm1;
       //br_fparm->SetAddress(&fparm1);
       if(store_bestfit_spectrum)
         br_bf_array->SetAddress(&pdf_array1);
     }
     
     //
     // Get the 1+N fit results for |dm32sq| (with mixing23 fixed to true input value)
     //
     
     fparm3->SetInitValue(conv::kPhypMixing23, input_mixing23);
     fparm3->SetInitValue(conv::kPhypDm23sq,   input_dm23sq);
     fparm3->FixParameter(conv::kPhypMixing23);

     //Get fitter
     FitInputs * fit_inputs3 = new FitInputs(skobs,skexp,fparm3);
     FitterABC * fitter3 = util::GetFitter(fit_inputs3,fit_option);
     fparm3->ResetCurrToInitValues(); // reset osc params     
     fitter3->DoTheFit();
     fparm3->ResetCurrToInitValues();
     migradmin3 = fitter3->FitStatusCode("migrad_min");
     migraderr3 = fitter3->FitStatusCode("migrad_err");
     cout<<"Migrad min and err codes for fit 3:   "<<migradmin3<<"    "<<migraderr3<<endl;
     bestfit_chisq_dm23sq_1N  = fit_inputs3->BestFitChisquare();

     if(store_bestfit_spectrum) {
       skexp->Write(pdf_array3, "pdf");
     }

     // Print-out best-fit params
     cout << "- Best-fit dm23sq (below) yields chisq = ";
     cout << bestfit_chisq_dm23sq_1N<<endl;
     fittries_dm23sq_1N = fparm3->GetFitTries();
     fittime_dm23sq_1N = fparm3->GetFitTime();
     bestfit_dm23sq_1N = fparm3->BestFitValue(conv::kPhypDm23sq);
     migradmin_dm23sq_1N        = migradmin3;
     migraderr_dm23sq_1N        = migraderr3;
     fparm3->ShowOnlyBestFitValuesInNextPrint();
     cout << *fparm3;

     //
     // Get the 2+N fit results
     //

     fparm4->SetInitValue(conv::kPhypMixing23, input_mixing23);
     fparm4->SetInitValue(conv::kPhypDm23sq,   input_dm23sq);

     // Get fitter
     FitInputs * fit_inputs4 = new FitInputs(skobs,skexp,fparm4);
     FitterABC * fitter4 = util::GetFitter(fit_inputs4,fit_option);

     // do the fit using the true osc param as a seed
     fparm4->ResetCurrToInitValues(); // reset osc params
     fitter4->DoTheFit();
     fparm4->ResetCurrToInitValues();
     migradmin4 = fitter4->FitStatusCode("migrad_min");
     migraderr4 = fitter4->FitStatusCode("migrad_err");
     cout<<"Migrad min and err codes for fit 4:   "<<migradmin4<<"    "<<migraderr4<<endl;
     bestfit_chisq_seed1  = fit_inputs4->BestFitChisquare();
     bestfit_chisq_seed2 = 0;

     // Save spectra (data and best-fit prediction)
     if(store_toy_expts) {
       skobs->Write(dat_array4, "obs");
     }
     if(store_bestfit_spectrum) {
       skexp->Write(pdf_array4, "pdf");
     }

     //Do the fit again at the mirror point
     migradmin5 = -1;
     if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsqTheta) {
       //get the correct intial values (best fit points from first fit)
       fparm5->CopyBestFitToInitValues(*fparm4);
       //set the mixing23 seed to be the mirror point of the initial best fit
       sinsqtheta23_bestfit = fparm4->BestFitValue(conv::kPhypMixing23);
       sinsqtheta13_bestfit = fparm4->BestFitValue(conv::kPhypMixing13);
       centre = 1 / (2 * (1 - sinsqtheta13_bestfit));
       reseed = sinsqtheta23_bestfit;
       if     (reseed < centre)
         reseed = centre + TMath::Abs(centre - sinsqtheta23_bestfit);
       else if(reseed > centre)
         reseed = centre - TMath::Abs(centre - sinsqtheta23_bestfit);
       fparm5->SetInitValue(conv::kPhypMixing23, reseed);
       // Get fitter
       FitInputs * fit_inputs5 = new FitInputs(skobs,skexp,fparm5);
       FitterABC * fitter5 = util::GetFitter(fit_inputs5,fit_option);
       fparm5->ResetCurrToInitValues(); // reset osc params
       fitter5->DoTheFit();
       bestfit_chisq_seed2 = fit_inputs5->BestFitChisquare();
       migradmin5 = fitter5->FitStatusCode("migrad_min");
       migraderr5 = fitter5->FitStatusCode("migrad_err");
       cout<<"Migrad min and err codes for fit 5:   "<<migradmin5<<"    "<<migraderr5<<endl;
       // Save spectra (best-fit prediction)
       if(store_bestfit_spectrum) {
         skexp->Write(pdf_array5, "pdf");
       }
     }
        
       // Print-out best-fit params, and set the "best of the two fits" branches correctly
       cout << "- Best-fit oscillation params (below) yield chisq = ";
       if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsqTheta) {
	 if(bestfit_chisq_seed2 < bestfit_chisq_seed1) {
	   cout << bestfit_chisq_seed2 << " from the second seed fit" << endl;
	   fittries_2N = fparm5->GetFitTries();
	   fittries_2N_alt = fparm4->GetFitTries();
	   fittime_2N = fparm5->GetFitTime();
	   fittime_2N_alt = fparm4->GetFitTime();
	   bestfit_chisq_2N    = bestfit_chisq_seed2;
	   bestfit_chisq_2N_alt    = bestfit_chisq_seed1;
	   bestfit_mixing23_2N = fparm5->BestFitValue(conv::kPhypMixing23);
           bestfit_mixing23_2N_alt = fparm4->BestFitValue(conv::kPhypMixing23);
	   bestfit_dm23sq_2N   = fparm5->BestFitValue(conv::kPhypDm23sq);
	   bestfit_dm23sq_2N_alt   = fparm4->BestFitValue(conv::kPhypDm23sq);
	   migradmin_2N        = migradmin5;
	   migradmin_2N_alt        = migradmin4;
           migraderr_2N        = migraderr5;
           migraderr_2N_alt        = migraderr4;
	   fparm5->ShowOnlyBestFitValuesInNextPrint();
	   cout << *fparm5;
	   //br_fparm->SetAddress(&fparm5);
	   if(store_bestfit_spectrum)
	     br_bf_array->SetAddress(&pdf_array5);
	 }
	 else {
	   cout << bestfit_chisq_seed1 << " from the first seed fit" << endl;
	   fittries_2N = fparm4->GetFitTries();
           fittries_2N_alt = fparm5->GetFitTries();
           fittime_2N = fparm4->GetFitTime();
           fittime_2N_alt = fparm5->GetFitTime();
           bestfit_chisq_2N    = bestfit_chisq_seed1;
	   bestfit_chisq_2N_alt    = bestfit_chisq_seed2;
	   bestfit_mixing23_2N = fparm4->BestFitValue(conv::kPhypMixing23);
	   bestfit_mixing23_2N_alt = fparm5->BestFitValue(conv::kPhypMixing23);
	   bestfit_dm23sq_2N   = fparm4->BestFitValue(conv::kPhypDm23sq);
	   bestfit_dm23sq_2N_alt   = fparm5->BestFitValue(conv::kPhypDm23sq);
	   migradmin_2N        = migradmin4;
	   migradmin_2N_alt        = migradmin5;
	   migraderr_2N        = migraderr4;
           migraderr_2N_alt        = migraderr5;
	   fparm4->ShowOnlyBestFitValuesInNextPrint();
	   cout << *fparm4;
	   //br_fparm->SetAddress(&fparm4);
	   if(store_bestfit_spectrum)
	     br_bf_array->SetAddress(&pdf_array4);
	 }
       }
       else {
	 cout << bestfit_chisq_seed1 << endl;
         fittries_2N = fparm4->GetFitTries();
	 fittime_2N = fparm4->GetFitTime();
	 bestfit_chisq_2N    = bestfit_chisq_seed1;
	 bestfit_mixing23_2N = fparm4->BestFitValue(conv::kPhypMixing23);
	 bestfit_dm23sq_2N   = fparm4->BestFitValue(conv::kPhypDm23sq);
	 migradmin_2N        = migradmin4;
         migraderr_2N        = migraderr4;
	 fparm4->ShowOnlyBestFitValuesInNextPrint();
	 cout << *fparm4;
	 //br_fparm->SetAddress(&fparm4);
	 if(store_bestfit_spectrum)
	   br_bf_array->SetAddress(&pdf_array4);
       }
     
     cout<<endl;

     //  Print out result
     cout<<"Overall results"<<endl;
     cout << "- Chisq (2+N fit) = " << bestfit_chisq_2N << endl;
     cout << "- Chisq for mixing23 (Best-fit param in 1+N fit:" 
	  << util::AsString(conv::kPhypMixing23) << " = " << bestfit_mixing23_1N
	  << ") = " << bestfit_chisq_mixing23_1N << endl;
     cout << "- DeltaChisq for mixing23 = " << bestfit_chisq_mixing23_1N - bestfit_chisq_2N << endl;
     cout << "- Chisq for dm23sq (Best-fit param in 1+N fit:"
          << util::AsString(conv::kPhypDm23sq)   << " = " << bestfit_dm23sq_1N   << " eV^2"
          << ") = " << bestfit_chisq_dm23sq_1N << endl;
     cout << "- DeltaChisq for dm23sq = " << bestfit_chisq_dm23sq_1N - bestfit_chisq_2N << endl;
       

     // Add all results in the tree
     fitnt->Fill();

     // Clean-up
     delete skobs;
     delete fit_inputs;
     delete fitter;
     delete fit_inputs3;
     delete fitter3;
     delete fit_inputs4;
     delete fitter4;
     if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsqTheta) {
       delete fit_inputs2;
       delete fitter2;
       delete fit_inputs5;
       delete fitter5;
     }

     }//iexperiment

  // Store ntuple with best-fit results
  f.cd();
  fitnt->Write();
  f.Close();
}
