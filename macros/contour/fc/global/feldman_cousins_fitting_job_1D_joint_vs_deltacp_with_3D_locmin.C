//
// ------------------------------------------------------------------------------------------------
// Run a job with some of the toy expts fits needed to calculate the \Delta\chi^{2} values
// for constructing a Feldman-Cousins contour in 1D vs deltaCP for the joint fit 
// applying reactor constraint
//
// Lorena Escudero
//
// ------------------------------------------------------------------------------------------------
//

void feldman_cousins_fitting_job_1D_joint_vs_deltacp_with_3D_locmin
( 
 const char * tag                     = "fc_test_wo_syst_pointX_0",
 //const char * tag                     = "fc_prova",
 int run                              = 0,                                      /* run number */
  int          nexperiments            = 1,                                  /* number of toy MC experiment fits */
 double       input_deltacp           = -1.,                                 /* generate toy MC with this range of 23 mixing (expressed in terms of either sinsq2theta, sinsqtheta, ...) */
  const char * config_file             = "analysis_configuration_test.xml",        /* analysis configuration file */
  bool         stat_fluctuations       = true,                                /* include statistical fluctuation in toy MC expts? */
  bool         tweak_systematics       = 1,                                /* generate toy MC with randomized systematics? */
 //const char * fit_param_list          = "nuosc3f_12fix+nuisance",  /* fit parameter list for the given physics hypothesis (see FitParamList.h) */
 const char * fit_param_list          = "nuosc3f_12fix",  /* fit parameter list for the given physics hypothesis (see FitParamList.h) */
  const char * fit_option              = "minuit+hesse",                            /* fitter option */
  double       tolerance               = 0.01,
  bool         store_toy_expts         = false,                                /* disable this if you don't want to store the toy MC expts     */
  bool         store_bestfit_spectrum  = false,                                /* disable this if you don't want to store the best-fit spectra */
  const char * dmsqscheme              = "abs23",                             /* abs23==kNuDmsq12andAbs23 fogli==kNuDmsqFogliLisi */
  const char * hierarchy               = "normal",                            /* normal or inverted */
  const char * mixparamchoice          = "sinsqtheta",                        /* sinsq2theta or sinsqtheta */
  const char * mixparamoctant          = "one",                                /* one==kNuOctantOne two==kNuOctantTwo */
 //const char * dchisqsurf_file         = "data_theta2theta_fit_joint_nh_final_fix_reactor_abs23_chisq_vs_theta23_dmsq_theta13_extra_3D.root" 
 //const char * dchisqsurf_file         = "LikeGrid_mixing23_Nbins5_dmsq_Nbins5_mixing13_Nbins5.root"
 const char * dchisqsurf_file         = "../../inputs/fc/LikeGrid_Data_Run1to4_NH_mixing23_dmsq_mixing13_41x21x21.root"
)
{
  using namespace std;
  using namespace valor;
  using namespace valor::conv;

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
  
  Opt::Instance()->SetDefaultParameterValues();

  // Take input oscillation likelihood grid for toy extraction
  TFile *dchisqfile = TFile::Open(dchisqsurf_file);
  TH3D *dchisq3D = (TH3D*)dchisqfile->Get("hLike3D");
  TH2D *dchisq2D = (TH2D*)dchisqfile->Get("hLike2D");

  // dm23sq, mixing23, mixing13 range and step
  double mixing23_min   = dchisq3D->GetXaxis()->GetBinCenter(1);
  double mixing23_max   = dchisq3D->GetXaxis()->GetBinCenter(dchisq3D->GetNbinsX());
  double mixing23_step  = dchisq3D->GetXaxis()->GetBinCenter(2) - dchisq3D->GetXaxis()->GetBinCenter(1);
  double dmsq_min     = dchisq3D->GetYaxis()->GetBinCenter(1);

  double dmsq_max     = dchisq3D->GetYaxis()->GetBinCenter(dchisq3D->GetNbinsY());
  double dmsq_step    = dchisq3D->GetYaxis()->GetBinCenter(2) - dchisq3D->GetYaxis()->GetBinCenter(1);
  double mixing13_min   = dchisq3D->GetZaxis()->GetBinCenter(1);
  double mixing13_max   = dchisq3D->GetZaxis()->GetBinCenter(dchisq3D->GetNbinsZ());
  double mixing13_step  = dchisq3D->GetZaxis()->GetBinCenter(2) - dchisq3D->GetZaxis()->GetBinCenter(1);
  
  // make sure that min and max are at the bin centre and not at the low/high edge of the first/last bin
  //add 0.001 to avoid rounding error when finding (int)
  int nbins_mixing23 = 1 + (int)(((mixing23_max - mixing23_min) / mixing23_step) + 0.001);
  int nbins_mixing13 = 1 + (int)(((mixing13_max - mixing13_min) / mixing13_step) + 0.001);
  int nbins_dmsq   = 1 + (int)(((dmsq_max   - dmsq_min  ) / dmsq_step  ) + 0.001);

  cout << " mixing23_min = " << mixing23_min-mixing23_step/2 << " , mixing23_max = " << mixing23_max-mixing23_step/2 << " and mixing23_step = " << mixing23_step << " and nbins_mixing23 = " << nbins_mixing23 << endl;
  cout << " dmsq_min = " << dmsq_min-dmsq_step/2 << " , dmsq_max = " << dmsq_max-dmsq_step/2 << " and dmsq_step = " << dmsq_step << " and nbins_dmsq = " << nbins_dmsq << endl;
  cout << " mixing13_min = " << mixing13_min-mixing13_step/2 << " , mixing13_max = " << mixing13_max-mixing13_step/2 << " and mixing13_step = " << mixing13_step << " and nbins_mixing13 = " << nbins_mixing13 << endl;

  const int mixing23_nbins = nbins_mixing23;
  const int mixing13_nbins = nbins_mixing13;
  const int dmsq_nbins = nbins_dmsq;
 
  // param space
  TH2D * osc_param_space = new TH2D("osc_param_space", "",
       mixing23_nbins, 
       mixing23_min - 0.5*mixing23_step, 
       mixing23_max + 0.5*mixing23_step,
       dmsq_nbins, 
       dmsq_min - 0.5*dmsq_step, 
       dmsq_max + 0.5*dmsq_step);

 TH1D * theta13_param_space = new TH1D("theta13_param_space", "",
 				       mixing13_nbins, 
 				       mixing13_min - 0.5*mixing13_step, 
 				       mixing13_max + 0.5*mixing13_step);

 TH1D * theta23_param_space = new TH1D("theta23_param_space", "",
 				       mixing23_nbins, 
 				       mixing23_min - 0.5*mixing23_step, 
 				       mixing23_max + 0.5*mixing23_step);

 TH1D * dmsq_param_space = new TH1D("dmsq_param_space", "",
 				       dmsq_nbins, 
 				       dmsq_min - 0.5*dmsq_step, 
 				       dmsq_max + 0.5*dmsq_step);

  int input_mixing23_bin,input_mixing13_bin, input_dmsq_bin;
  double input_mixing23,input_mixing13,input_dmsq;

  double totalLikelihood = 0.0;
  for(int mixing13_bin=1; mixing13_bin<=mixing13_nbins; mixing13_bin++){
    for(int dmsq_bin=1; dmsq_bin<=dmsq_nbins; dmsq_bin++){
      for(int mixing23_bin=1; mixing23_bin<=mixing23_nbins; mixing23_bin++){
 	totalLikelihood += exp(-0.5 * dchisq3D->GetBinContent(mixing23_bin, dmsq_bin,mixing13_bin));	
      }
    }
  }
  
  double cumulNormLikelihood[dmsq_nbins*mixing23_nbins*mixing13_nbins+1];

  cumulNormLikelihood[0] = 0.0;

  for(int mixing13_bin=1; mixing13_bin<=mixing13_nbins; mixing13_bin++){
    for(int dmsq_bin=1; dmsq_bin<=dmsq_nbins; dmsq_bin++){
      for(int mixing23_bin=1; mixing23_bin<=mixing23_nbins; mixing23_bin++){
 	int index = (mixing13_bin-1)*mixing23_nbins*dmsq_nbins + (dmsq_bin-1)*mixing23_nbins + mixing23_bin;
 	//	cout << "index = " << index << endl;
 	cumulNormLikelihood[index] = cumulNormLikelihood[index-1] 
 	  + (exp(-0.5 * dchisq3D->GetBinContent(mixing23_bin, dmsq_bin, mixing13_bin)) / totalLikelihood);
 	//cout << " cumulNormLikelihood["<< index<<"] = " << cumulNormLikelihood[index] << ", with exp = " << (exp(-0.5 * dchisq3D->GetBinContent(mixing23_bin, dmsq_bin, mixing13_bin))) << endl;
 	//prevent cumulative normalised likelihood being > 1.0
 	cumulNormLikelihood[index] = TMath::Min(1.0, cumulNormLikelihood[index]);
      }
    }
  }
  
  cumulNormLikelihood[dmsq_nbins*mixing23_nbins*mixing13_nbins] = 1.0;
 
  // cout<<endl;
  // cout<<"This table shows likelihood, normalised likelihood and cumulative normalised likelihood for the 2D chisq surface for each bin" <<endl;
  // cout<<endl;
  // for(int mixing13_bin=1; mixing13_bin<=mixing13_nbins; mixing13_bin++){
  //   for(int dmsq_bin=1; dmsq_bin<=dmsq_nbins; dmsq_bin++){
  //     for(int mixing23_bin=1; mixing23_bin<=mixing23_nbins; mixing23_bin++){
  // 	int index = (mixing13_bin-1)*mixing23_nbins*dmsq_nbins + (dmsq_bin-1)*mixing23_nbins + mixing23_bin;
  // 	cout<<mixing23_bin<<"  ,  "<<dmsq_bin << "  ,  " << mixing13_bin
  // 	    <<"    "<<dchisq3D->GetBinContent(mixing23_bin, dmsq_bin, mixing13_bin) // -2 log likelihood
  // 	    <<"    "<<exp(-0.5 * dchisq3D->GetBinContent(mixing23_bin, dmsq_bin, mixing13_bin)) // likelihood
  // 	    <<"    "<<"    "<<(exp(-0.5 * dchisq3D->GetBinContent(mixing23_bin, dmsq_bin, mixing13_bin)) / totalLikelihood) // normalized likelihood
  // 	    <<"    "<<cumulNormLikelihood[index]<<endl; // cumulative likelihood
  //     }
  //   }
  // }

  // Print settings
  cout << "\n[Job inputs]" << endl;
  cout << "- Job: " << tag << endl;
  cout << "- Reading analysis configuration from: " << config << endl;
  cout << "- Toy-MC generation" << endl;
  cout << "   - Number of toy experiments: " << nexperiments << endl; 
  cout << "   - Randomize systematics : "      << ((tweak_systematics) ?      "YES " : "NO") << endl;
  cout << "   - Statistical fluctuations: "    << ((stat_fluctuations) ?      "YES " : "NO") << endl;
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
       << ", "   << osc_param_space->GetYaxis()->GetBinCenter(dmsq_nbins)
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
  // dcp
  FitParamList * fparm      = util::NamedParameterList(fit_param_list); //only used to initialise branch
  FitParamList * fparm1     = util::NamedParameterList(fit_param_list); // for 1+N fit of mixing23 using the first seed
  FitParamList * fparm2     = util::NamedParameterList(fit_param_list); // for 1+N fit of mixing23 using the mirror of best fit seed
  FitParamList * fparm3     = util::NamedParameterList(fit_param_list); // for 1+N fit of |dm32sq| 
  FitParamList * fparm4     = util::NamedParameterList(fit_param_list); // for 2+N fit using the first seed of mixing23
  FitParamList * fparm5     = util::NamedParameterList(fit_param_list); // for 2+N fit using the mirror of best fit seed of mixing23
  FitParamList * fparm6     = util::NamedParameterList(fit_param_list); // for 2+N fit using the mirror of best fit seed of mixing23
  // dcp + Pi
  FitParamList * fparm1_dcp     = util::NamedParameterList(fit_param_list); // for 1+N fit of mixing23 using the first seed
  FitParamList * fparm2_dcp     = util::NamedParameterList(fit_param_list); // for 1+N fit of mixing23 using the mirror of best fit seed
  FitParamList * fparm3_dcp     = util::NamedParameterList(fit_param_list); // for 1+N fit of |dm32sq| 
  FitParamList * fparm4_dcp     = util::NamedParameterList(fit_param_list); // for 2+N fit using the first seed of mixing23
  FitParamList * fparm5_dcp     = util::NamedParameterList(fit_param_list); // for 2+N fit using the mirror of best fit seed of mixing23
  FitParamList * fparm6_dcp     = util::NamedParameterList(fit_param_list); // for 2+N fit using the mirror of best fit seed of mixing23

  PhysHypothesis_t hyp = fparm1->PhysicsHypothesis();
  cout << "- The following physics hypothesis is used in the fit: " << util::AsString(hyp) << endl;

  // Get superk spectrum predictor & configure
  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read(config.Data());

  if(hierarchy == "normal")
    Opt::Instance()->gMassHierarchyAssumed = util::HierarchyFromString("inverted");
  else Opt::Instance()->gMassHierarchyAssumed = util::HierarchyFromString("normal");

  SuperkExpectation * skexp_mh = new SuperkExpectation;
  skexp_mh->Read(config.Data());

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
  TObjArray * pdf_array6 = new TObjArray;
  pdf_array6->SetOwner(false);

  TObjArray * dat_array_dcp = new TObjArray;
  TObjArray * pdf_array_dcp = new TObjArray; //only used to initialise branch
  dat_array_dcp->SetOwner(false);
  pdf_array_dcp->SetOwner(false);
  TObjArray * pdf_array1_dcp = new TObjArray;
  pdf_array1_dcp->SetOwner(false);
  TObjArray * pdf_array2_dcp = new TObjArray;
  pdf_array2_dcp->SetOwner(false);
  TObjArray * pdf_array3_dcp = new TObjArray;
  pdf_array3_dcp->SetOwner(false);
  TObjArray * pdf_array4_dcp = new TObjArray;
  pdf_array4_dcp->SetOwner(false);
  TObjArray * pdf_array5_dcp = new TObjArray;
  pdf_array5_dcp->SetOwner(false);
  TObjArray * pdf_array6_dcp = new TObjArray;
  pdf_array6_dcp->SetOwner(false);

  // Create file to store results
  TFile f(Form("%s_run_%d.root",tag,run), "recreate"); 

  //numbers to store results in
  int nsk_1Rmu,nsk_1Re;
  int fittries_true, fittries_true_alt, fittries_min, fittries_min_oct, fittries_min_mh,fittries_min_mh_oct,
    fittries_true_dcp, fittries_true_alt_dcp, fittries_min_dcp, fittries_min_oct_dcp, fittries_min_mh_dcp,fittries_min_mh_oct_dcp;
  double fittime_true, fittime_true_alt, fittime_min, fittime_min_oct, fittime_min_mh,fittime_min_mh_oct,
    fittime_true_dcp, fittime_true_alt_dcp, fittime_min_dcp, fittime_min_oct_dcp, fittime_min_mh_dcp,fittime_min_mh_oct_dcp;
  int migradmin_true, migradmin_true_alt, migradmin_min, migradmin_min_oct, migradmin_min_mh,migradmin_min_mh_oct,
    migradmin_true_dcp, migradmin_true_alt_dcp, migradmin_min_dcp, migradmin_min_oct_dcp, migradmin_min_mh_dcp,migradmin_min_mh_oct_dcp;
  int migraderr_true, migraderr_true_alt, migraderr_min, migraderr_min_oct, migraderr_min_mh,migraderr_min_mh_oct,
    migraderr_true_dcp, migraderr_true_alt_dcp, migraderr_min_dcp, migraderr_min_oct_dcp, migraderr_min_mh_dcp,migraderr_min_mh_oct_dcp;
  int migradmin1, migradmin2, migradmin3, migradmin4, migradmin5, migradmin6,
    migradmin1_dcp, migradmin2_dcp, migradmin3_dcp, migradmin4_dcp, migradmin5_dcp, migradmin6_dcp;
  int migraderr1, migraderr2, migraderr3, migraderr4, migraderr5, migraderr6,
    migraderr1_dcp, migraderr2_dcp, migraderr3_dcp, migraderr4_dcp, migraderr5_dcp, migraderr6_dcp;  //////
  double bestfit_chisq_seed1, bestfit_chisq_seed2, bestfit_chisq_seed3,bestfit_chisq_seed4,
    bestfit_chisq_seed1_dcp, bestfit_chisq_seed2_dcp, bestfit_chisq_seed3_dcp,bestfit_chisq_seed4_dcp;
  double sinsqtheta23_bestfit, sinsqtheta13_bestfit, centre, reseed;
  double bestfit_chisq_true, bestfit_chisq_true_alt, bestfit_chisq_min, bestfit_chisq_min_oct, bestfit_chisq_min_mh,bestfit_chisq_min_mh_oct,
    bestfit_chisq_true_dcp, bestfit_chisq_true_alt_dcp, bestfit_chisq_min_dcp, bestfit_chisq_min_oct_dcp, bestfit_chisq_min_mh_dcp,bestfit_chisq_min_mh_oct_dcp;
  double bestfit_mixing23_true, bestfit_mixing23_true_alt, bestfit_mixing23_min, bestfit_mixing23_min_oct, bestfit_mixing23_min_mh,bestfit_mixing23_min_mh_oct,
    bestfit_mixing23_true_dcp, bestfit_mixing23_true_alt_dcp, bestfit_mixing23_min_dcp, bestfit_mixing23_min_oct_dcp, bestfit_mixing23_min_mh_dcp,bestfit_mixing23_min_mh_oct_dcp;
  double bestfit_mixing13_true, bestfit_mixing13_true_alt, bestfit_mixing13_min, bestfit_mixing13_min_oct, bestfit_mixing13_min_mh,bestfit_mixing13_min_mh_oct,
    bestfit_mixing13_true_dcp, bestfit_mixing13_true_alt_dcp, bestfit_mixing13_min_dcp, bestfit_mixing13_min_oct_dcp, bestfit_mixing13_min_mh_dcp,bestfit_mixing13_min_mh_oct_dcp;
  double bestfit_dmsq_true, bestfit_dmsq_true_alt, bestfit_dmsq_min, bestfit_dmsq_min_oct, bestfit_dmsq_min_mh,bestfit_dmsq_min_mh_oct,
    bestfit_dmsq_true_dcp, bestfit_dmsq_true_alt_dcp, bestfit_dmsq_min_dcp, bestfit_dmsq_min_oct_dcp, bestfit_dmsq_min_mh_dcp,bestfit_dmsq_min_mh_oct_dcp;
  double bestfit_deltacp_true, bestfit_deltacp_true_alt, bestfit_deltacp_min, bestfit_deltacp_min_oct, bestfit_deltacp_min_mh,bestfit_deltacp_min_mh_oct,
    bestfit_deltacp_true_dcp, bestfit_deltacp_true_alt_dcp, bestfit_deltacp_min_dcp, bestfit_deltacp_min_oct_dcp, bestfit_deltacp_min_mh_dcp,bestfit_deltacp_min_mh_oct_dcp;

  // Create tree to store outputs
  TTree * fitnt = new TTree("fitnt","Results of VaLOR data and/or toy-MC experiment fits");
  /*
  fitnt->Branch("best_fit_params",           "t2kana::FitParamList",  &fparm,          1280000,0);
  fitnt->Branch("best_fit_params_fit1",      "t2kana::FitParamList",  &fparm1,         1280000,0);
  fitnt->Branch("best_fit_params_fit2",      "t2kana::FitParamList",  &fparm2,         1280000,0);
  fitnt->Branch("best_fit_params_fit3",      "t2kana::FitParamList",  &fparm3,         1280000,0);
  fitnt->Branch("best_fit_params_fit4",      "t2kana::FitParamList",  &fparm4,         1280000,0);
  fitnt->Branch("best_fit_params_fit5",      "t2kana::FitParamList",  &fparm5,         1280000,0);
  fitnt->Branch("best_fit_params_fit6",      "t2kana::FitParamList",  &fparm5,         1280000,0);
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
  
  fitnt->Branch("nsk_1Rmu",                  &nsk_1Rmu);
  fitnt->Branch("nsk_1Re",                   &nsk_1Re);
  fitnt->Branch("input_mixing23_bin",        &input_mixing23_bin);
  fitnt->Branch("input_mixing13_bin",        &input_mixing13_bin);
  fitnt->Branch("input_dmsq_bin",            &input_dmsq_bin);
  fitnt->Branch("input_mixing23",            &input_mixing23);
  fitnt->Branch("input_mixing13",            &input_mixing13);
  fitnt->Branch("input_dmsq",                &input_dmsq);
  fitnt->Branch("input_deltacp",             &input_deltacp);

  fitnt->Branch("fittries_true",             &fittries_true);
  fitnt->Branch("fittries_true_alt",         &fittries_true_alt);
  fitnt->Branch("fittries_min",              &fittries_min);
  fitnt->Branch("fittries_min_oct",          &fittries_min_oct);
  fitnt->Branch("fittries_min_mh",           &fittries_min_mh);
  fitnt->Branch("fittries_min_mh_oct",       &fittries_min_mh_oct);
  fitnt->Branch("fittime_true",              &fittime_true);
  fitnt->Branch("fittime_true_alt",          &fittime_true_alt);
  fitnt->Branch("fittime_min",               &fittime_min);
  fitnt->Branch("fittime_min_oct",           &fittime_min_oct);
  fitnt->Branch("fittime_min_mh",            &fittime_min_mh);
  fitnt->Branch("fittime_min_mh_oct",        &fittime_min_mh_oct);
  fitnt->Branch("migradmin_true",            &migradmin_true);
  fitnt->Branch("migradmin_true_alt",        &migradmin_true_alt);
  fitnt->Branch("migradmin_min",             &migradmin_min);
  fitnt->Branch("migradmin_min_oct",         &migradmin_min_oct);
  fitnt->Branch("migradmin_min_mh",          &migradmin_min_mh);
  fitnt->Branch("migradmin_min_mh_oct",      &migradmin_min_mh_oct);
  fitnt->Branch("migraderr_true",            &migraderr_true);
  fitnt->Branch("migraderr_true_alt",        &migraderr_true_alt);
  fitnt->Branch("migraderr_min",             &migraderr_min);
  fitnt->Branch("migraderr_min_oct",         &migraderr_min_oct);
  fitnt->Branch("migraderr_min_mh",          &migraderr_min_mh);
  fitnt->Branch("migraderr_min_mh_oct",      &migraderr_min_mh_oct);
  fitnt->Branch("bestfit_chisq_true",        &bestfit_chisq_true);
  fitnt->Branch("bestfit_chisq_true_alt",    &bestfit_chisq_true_alt);
  fitnt->Branch("bestfit_chisq_min",         &bestfit_chisq_min);
  fitnt->Branch("bestfit_chisq_min_oct",     &bestfit_chisq_min_oct);
  fitnt->Branch("bestfit_chisq_min_mh",      &bestfit_chisq_min_mh);
  fitnt->Branch("bestfit_chisq_min_mh_oct",  &bestfit_chisq_min_mh_oct);
  fitnt->Branch("bestfit_mixing23_true",     &bestfit_mixing23_true);
  fitnt->Branch("bestfit_mixing23_true_alt", &bestfit_mixing23_true_alt);
  fitnt->Branch("bestfit_mixing23_min",      &bestfit_mixing23_min);
  fitnt->Branch("bestfit_mixing23_min_oct",  &bestfit_mixing23_min_oct);
  fitnt->Branch("bestfit_mixing23_min_mh",   &bestfit_mixing23_min_mh);
  fitnt->Branch("bestfit_mixing23_min_mh_oct",  &bestfit_mixing23_min_mh_oct);
  fitnt->Branch("bestfit_mixing13_true",     &bestfit_mixing13_true);
  fitnt->Branch("bestfit_mixing13_true_alt", &bestfit_mixing13_true_alt);
  fitnt->Branch("bestfit_mixing13_min",      &bestfit_mixing13_min);
  fitnt->Branch("bestfit_mixing13_min_oct",  &bestfit_mixing13_min_oct);
  fitnt->Branch("bestfit_mixing13_min_mh",   &bestfit_mixing13_min_mh);
  fitnt->Branch("bestfit_mixing13_min_mh_oct",  &bestfit_mixing13_min_mh_oct);
  fitnt->Branch("bestfit_dmsq_true",         &bestfit_dmsq_true);
  fitnt->Branch("bestfit_dmsq_true_alt",     &bestfit_dmsq_true_alt);
  fitnt->Branch("bestfit_dmsq_min",          &bestfit_dmsq_min);
  fitnt->Branch("bestfit_dmsq_min_oct",      &bestfit_dmsq_min_oct);
  fitnt->Branch("bestfit_dmsq_min_mh",       &bestfit_dmsq_min_mh);
  fitnt->Branch("bestfit_dmsq_min_mh_oct",   &bestfit_dmsq_min_mh_oct);
  fitnt->Branch("bestfit_deltacp_true",      &bestfit_deltacp_true);
  fitnt->Branch("bestfit_deltacp_true_alt",  &bestfit_deltacp_true_alt);
  fitnt->Branch("bestfit_deltacp_min",       &bestfit_deltacp_min);
  fitnt->Branch("bestfit_deltacp_min_oct",   &bestfit_deltacp_min_oct);
  fitnt->Branch("bestfit_deltacp_min_mh",    &bestfit_deltacp_min_mh);
  fitnt->Branch("bestfit_deltacp_min_mh_oct",  &bestfit_deltacp_min_mh_oct);

  //TBranch * br_fparm     = fitnt->GetBranch("best_fit_params");
  if(store_bestfit_spectrum)
    TBranch * br_bf_array  = fitnt->GetBranch("best_fit_pdf");

  double uniformRand;

  TRandom3 rand(0);

  // Loop over toy-MC experiments
  for(int iexperiment = 0; iexperiment < nexperiments; iexperiment++) {

    cout << endl;
    cout << "//////////////////////////////////" << endl;
    cout << "//                              //" << endl;
    cout << "//   *** Pseudo-experiment = " << iexperiment << " //"<< endl;
    cout << "//                              //" << endl;
    cout << "//////////////////////////////////" << endl;
    cout << endl;






    // //     
    // // Random selection of (mixing23,dmsq,mixing13) bin 
    // //     
    
    // uniformRand = rand.Uniform(0.0, 1.0);
   
    // cout << "uniformRand = " << uniformRand << endl;
    
    //  for(int mixing13_bin=1; mixing13_bin<=mixing13_nbins; mixing13_bin++){
    //    for(int dmsq_bin=1; dmsq_bin<=dmsq_nbins; dmsq_bin++){
    // 	 for(int mixing23_bin=1; mixing23_bin<=mixing23_nbins; mixing23_bin++){
    // 	   int index = (mixing13_bin-1)*mixing23_nbins*dmsq_nbins + (dmsq_bin-1)*mixing23_nbins + mixing23_bin;
    // 	   if(uniformRand >= cumulNormLikelihood[index-1] && uniformRand <= cumulNormLikelihood[index]){
    // 	     input_dmsq_bin = dmsq_bin;
    // 	     input_mixing23_bin = mixing23_bin;
    // 	     input_mixing13_bin = mixing13_bin;
    // 	   }
    // 	 }
    //    }
    //  }

    //  cout << endl;
    //  cout << endl;
    //  cout<<"Uniform random number thrown between 0 and 1 = "<<uniformRand<<" which means the input (mixing23,dmsq,mixing13) bin is "<<input_mixing23_bin <<" , " << input_dmsq_bin<< " and " << input_mixing13_bin << endl;
    //  cout << "Extracted Likelihood = " << dchisq3D->GetBinContent(input_mixing23_bin, input_dmsq_bin, input_mixing13_bin) << endl;
    //  cout << "mixing23 = " << dchisq3D->GetXaxis()->GetBinCenter(input_mixing23_bin) << endl;
    //  cout << "dmsq = " << dchisq3D->GetYaxis()->GetBinCenter(input_dmsq_bin) << endl;
    //  cout << "mixing13 = " << dchisq3D->GetZaxis()->GetBinCenter(input_mixing13_bin) << endl;
    //  cout << endl;
    //  cout << endl;

     // // Sanity check on bin range
     // input_mixing23_bin = TMath::Max(1,              input_mixing23_bin);
     // input_mixing23_bin = TMath::Min(mixing23_nbins, input_mixing23_bin);
     // input_mixing13_bin = TMath::Max(1,              input_mixing13_bin);
     // input_mixing13_bin = TMath::Min(mixing13_nbins, input_mixing13_bin);
     // input_dmsq_bin     = TMath::Max(1,              input_dmsq_bin  );
     // input_dmsq_bin     = TMath::Min(dmsq_nbins,     input_dmsq_bin  );

     // cout << "- Toy-MC > Selected 23-sector and theta13 params bin: " 
     //      << util::AsString(conv::kPhypMixing23) << " bin # = " << input_mixing23_bin << ", "
     //      << util::AsString(conv::kPhypMixing13) << " bin # = " << input_mixing13_bin << ", "
     //      << util::AsString(conv::kPhypDm23sq)   << " bin # = " << input_dmsq_bin
     //      << endl;

     // //input_mixing13 = rand.Gaus(glob_theta13, sigma_theta13);

     // // Get dm23sq or mixing23 at the centre of the randomly selected bin
     // input_dmsq   = dchisq3D->GetYaxis()->GetBinCenter(input_dmsq_bin);
     // input_mixing23 = dchisq3D->GetXaxis()->GetBinCenter(input_mixing23_bin);
     // input_mixing13 = dchisq3D->GetZaxis()->GetBinCenter(input_mixing13_bin);

     // osc_param_space->Fill(input_mixing23,input_dmsq);
     // theta13_param_space->Fill(input_mixing13);
     // theta23_param_space->Fill(input_mixing23);
     // dmsq_param_space->Fill(input_dmsq);

     // //input_dmsq = input_dmsq*1.E-3;

     // input_mixing23 = TMath::Max(0.,input_mixing23); // protection from rnd-err
     // input_mixing23 = TMath::Min(1.,input_mixing23); // protection from rnd-err
     // input_mixing13 = TMath::Max(0.,input_mixing13); // protection from rnd-err
     // input_mixing13 = TMath::Min(1.,input_mixing13); // protection from rnd-err






















    input_mixing13          = conv::kGlobFitSinsqTheta13;                                /* 13-mixing, sinq2theta or sinsqtheta depending on user choice, if -999 then randomize */
    input_mixing23          = conv::kGlobFitSinsqTheta23;                                /* 23-mixing, sinq2theta or sinsqtheta depending on user choice */
    input_dmsq            = conv::kGlobFitDm23sq;                                        /* |Dm^2_{23}| */
 
    fparm1->Update(kPhypMixing13, input_mixing13);
    fparm1->Update(kPhypMixing23, input_mixing23);
    fparm1->Update(kPhypDm23sq,   input_dmsq);

    fparm2->Update(kPhypMixing13, input_mixing13);
    fparm2->Update(kPhypMixing23, input_mixing23);
    fparm2->Update(kPhypDm23sq,   input_dmsq);

    fparm3->Update(kPhypMixing13, input_mixing13);
    fparm3->Update(kPhypMixing23, input_mixing23);
    fparm3->Update(kPhypDm23sq,   input_dmsq);

    fparm4->Update(kPhypMixing13, input_mixing13);
    fparm4->Update(kPhypMixing23, input_mixing23);
    fparm4->Update(kPhypDm23sq,   input_dmsq);

    fparm5->Update(kPhypMixing13, input_mixing13);
    fparm5->Update(kPhypMixing23, input_mixing23);
    fparm5->Update(kPhypDm23sq,   input_dmsq);

    fparm6->Update(kPhypMixing13, input_mixing13);
    fparm6->Update(kPhypMixing23, input_mixing23);
    fparm6->Update(kPhypDm23sq,   input_dmsq);

    

    




    // Set deltacp min and max values
    //double dcp_max = conv::kpi+conv::kpi/2.;
    //double dcp_min = -conv::kpi+conv::kpi/2.;
    double dcp_max = +2*conv::kpi;
    double dcp_min = -2*conv::kpi;
    fparm1->SetMinValue(conv::kPhypCPV,dcp_min);
    fparm1->SetMaxValue(conv::kPhypCPV,dcp_max);
    fparm2->SetMinValue(conv::kPhypCPV,dcp_min);
    fparm2->SetMaxValue(conv::kPhypCPV,dcp_max);
    fparm3->SetMinValue(conv::kPhypCPV,dcp_min);
    fparm3->SetMaxValue(conv::kPhypCPV,dcp_max);
    fparm4->SetMinValue(conv::kPhypCPV,dcp_min);
    fparm4->SetMaxValue(conv::kPhypCPV,dcp_max);
    fparm5->SetMinValue(conv::kPhypCPV,dcp_min);
    fparm5->SetMaxValue(conv::kPhypCPV,dcp_max);
    fparm6->SetMinValue(conv::kPhypCPV,dcp_min);
    fparm6->SetMaxValue(conv::kPhypCPV,dcp_max);
    fparm3_dcp->SetMinValue(conv::kPhypCPV,dcp_min);
    fparm3_dcp->SetMaxValue(conv::kPhypCPV,dcp_max);
    fparm4_dcp->SetMinValue(conv::kPhypCPV,dcp_min);
    fparm4_dcp->SetMaxValue(conv::kPhypCPV,dcp_max);
    fparm5_dcp->SetMinValue(conv::kPhypCPV,dcp_min);
    fparm5_dcp->SetMaxValue(conv::kPhypCPV,dcp_max);
    fparm6_dcp->SetMinValue(conv::kPhypCPV,dcp_min);
    fparm6_dcp->SetMaxValue(conv::kPhypCPV,dcp_max);











     cout << "- Toy-MC > True 23-sector params: " 
          << util::AsString(conv::kPhypMixing23) << " = " << input_mixing23 << ", "
          << util::AsString(conv::kPhypMixing13) << " = " << input_mixing13 << ", "
          << util::AsString(conv::kPhypDm23sq)   << " = " << input_dmsq   << " eV^2"
          << endl
  	  << endl;
      
     //
     // Set the generated oscillation parameter values (used to generate toy-MC expt)
     //
     fparm1->Update(conv::kPhypMixing23, input_mixing23);
     fparm1->Update(conv::kPhypDm23sq,   input_dmsq);
     fparm1->Update(conv::kPhypMixing13, input_mixing13);
     fparm1->Update(conv::kPhypCPV,  input_deltacp);

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

     nsk_1Rmu = skobs->NumOfEvents(conv::kSKsFCFV1RingMuLike);
     nsk_1Re  = skobs->NumOfEvents(conv::kSKsFCFV1RingELike);



     ////////////////////////////// delta_cp //////////////////////////////////

     //
     // Get the 1+N fit results for deltacp (with deltacp fixed to true input value)
     //

     cout << endl;
     cout << "////////////////////////////////////////////////////////////////" << endl;
     cout << endl;
     cout << "deltacp: True" << endl;
     cout << endl;

     fparm1->SetInitValue(conv::kPhypMixing23, input_mixing23);
     fparm1->SetInitValue(conv::kPhypDm23sq,   input_dmsq);
     fparm1->SetInitValue(conv::kPhypMixing13, input_mixing13);
     fparm1->SetInitValue(conv::kPhypCPV,  input_deltacp);
     fparm1->FixParameter(conv::kPhypCPV);    
     fparm1->SetNuisance(conv::kPhypMixing13,true);//reactor!
     
     // stay in [-Pi;+Pi] region
     if(fparm1->InitValue(conv::kPhypCPV) < -conv::kpi)
       fparm1->SetInitValue(conv::kPhypCPV,fparm1->InitValue(conv::kPhypCPV)+2*conv::kpi);
     else if(fparm1->InitValue(conv::kPhypCPV) > conv::kpi)
       fparm1->SetInitValue(conv::kPhypCPV,fparm1->InitValue(conv::kPhypCPV)-2*conv::kpi);
 
     fparm1->ResetCurrToInitValues(); // reset osc params 
     cout << *fparm1 << endl;
     
     // Get fitter 
     FitInputs * fit_inputs = new FitInputs(skobs,skexp,fparm1);
     FitterABC * fitter = util::GetFitter(fit_inputs,fit_option);
     fitter->SetTol(tolerance);
     fitter->DoTheFit();
     fparm1->ResetCurrToInitValues();
     migradmin1 = fitter->FitStatusCode("migrad_min");
     migraderr1 = fitter->FitStatusCode("migrad_err");
     cout<<"Migrad min and err codes for fit 1:   "<<migradmin1<<"    "<<migraderr1<<endl;
     bestfit_chisq_seed1  = fparm1->GetBestFitStatistic();
     bestfit_chisq_seed2 = 0;

     // Save spectra (data and best-fit prediction)
     if(store_toy_expts) {
        skobs->Write(dat_array, "obs");
     }
     if(store_bestfit_spectrum) {
        skexp->Write(pdf_array1, "pdf");
     }

     
     //Do the fit again at the mirror point

     cout << endl;
     cout << "////////////////////////////////////////////////////////////////" << endl;
     cout << endl;
     cout << "deltacp: True - mirror" << endl;
     cout << endl;

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
       fparm2->SetInitValue(conv::kPhypDm23sq,   input_dmsq);
       fparm2->SetInitValue(conv::kPhypMixing13, input_mixing13);
       fparm2->SetInitValue(conv::kPhypCPV,  input_deltacp);
       fparm2->FixParameter(conv::kPhypCPV);
       fparm2->SetNuisance(conv::kPhypMixing13,true);//reactor!
       
       // stay in [-Pi;+Pi] region
       if(fparm2->InitValue(conv::kPhypCPV) < -conv::kpi)
	 fparm2->SetInitValue(conv::kPhypCPV,fparm2->InitValue(conv::kPhypCPV)+2*conv::kpi);
       else if(fparm2->InitValue(conv::kPhypCPV) > conv::kpi)
	 fparm2->SetInitValue(conv::kPhypCPV,fparm2->InitValue(conv::kPhypCPV)-2*conv::kpi);

       fparm2->ResetCurrToInitValues(); // reset osc params 
       cout << *fparm2 << endl;
     
       // Get fitter
       FitInputs * fit_inputs2 = new FitInputs(skobs,skexp,fparm2);
       FitterABC * fitter2 = util::GetFitter(fit_inputs2,fit_option);
       fitter2->SetTol(tolerance);
       fitter2->DoTheFit();
       bestfit_chisq_seed2 = fparm2->GetBestFitStatistic();
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
         fittries_true = fparm2->GetFitTries();
         fittries_true_alt = fparm1->GetFitTries();
         fittime_true = fparm2->GetFitTime();
     	 fittime_true_alt = fparm1->GetFitTime();
     	 bestfit_chisq_true    = bestfit_chisq_seed2;
     	 bestfit_chisq_true_alt    = bestfit_chisq_seed1;
     	 bestfit_mixing23_true = fparm2->BestFitValue(conv::kPhypMixing23);
     	 bestfit_mixing23_true_alt = fparm1->BestFitValue(conv::kPhypMixing23);
     	 bestfit_mixing13_true = fparm2->BestFitValue(conv::kPhypMixing13);
     	 bestfit_mixing13_true_alt = fparm1->BestFitValue(conv::kPhypMixing13);	 
     	 bestfit_dmsq_true = fparm2->BestFitValue(conv::kPhypDm23sq);
     	 bestfit_dmsq_true_alt = fparm1->BestFitValue(conv::kPhypDm23sq);
     	 bestfit_deltacp_true = fparm2->BestFitValue(conv::kPhypCPV);
     	 bestfit_deltacp_true_alt = fparm1->BestFitValue(conv::kPhypCPV);
     	 migradmin_true        = migradmin2;
     	 migradmin_true_alt        = migradmin1;
         migraderr_true        = migraderr2;
         migraderr_true_alt        = migraderr1;
     	 fparm2->ShowOnlyBestFitValuesInNextPrint();
     	 cout << *fparm2;
     	 //br_fparm->SetAddress(&fparm2);
     	 if(store_bestfit_spectrum)
     	   br_bf_array->SetAddress(&pdf_array2);
       }
       else {
     	 cout << bestfit_chisq_seed1 << " from the first seed fit" << endl;
         fittries_true = fparm1->GetFitTries();
         fittries_true_alt = fparm2->GetFitTries();
         fittime_true = fparm1->GetFitTime();
     	 fittime_true_alt = fparm2->GetFitTime();
     	 bestfit_chisq_true    = bestfit_chisq_seed1;
     	 bestfit_chisq_true_alt    = bestfit_chisq_seed2;
     	 bestfit_mixing23_true = fparm1->BestFitValue(conv::kPhypMixing23);
     	 bestfit_mixing23_true_alt = fparm2->BestFitValue(conv::kPhypMixing23);
     	 bestfit_mixing13_true = fparm1->BestFitValue(conv::kPhypMixing13);
     	 bestfit_mixing13_true_alt = fparm2->BestFitValue(conv::kPhypMixing13);	 
     	 bestfit_dmsq_true = fparm1->BestFitValue(conv::kPhypDm23sq);
     	 bestfit_dmsq_true_alt = fparm2->BestFitValue(conv::kPhypDm23sq);
     	 bestfit_deltacp_true = fparm1->BestFitValue(conv::kPhypCPV);
     	 bestfit_deltacp_true_alt = fparm2->BestFitValue(conv::kPhypCPV);
     	 migradmin_true        = migradmin1;
     	 migradmin_true_alt        = migradmin2;
         migraderr_true        = migraderr1;
         migraderr_true_alt        = migraderr2;
     	 fparm1->ShowOnlyBestFitValuesInNextPrint();
     	 cout << *fparm1;
     	 //br_fparm->SetAddress(&fparm1);
     	 if(store_bestfit_spectrum)
     	   br_bf_array->SetAddress(&pdf_array1);
       }
     }
     else {
       cout << bestfit_chisq_seed1 << endl;
       fittries_true = fparm1->GetFitTries();
       fittime_true = fparm1->GetFitTime();
       bestfit_chisq_true    = bestfit_chisq_seed1;
       bestfit_mixing23_true = fparm1->BestFitValue(conv::kPhypMixing23);
       bestfit_mixing13_true = fparm1->BestFitValue(conv::kPhypMixing13);
       bestfit_dmsq_true = fparm1->BestFitValue(conv::kPhypDm23sq);
       bestfit_deltacp_true = fparm1->BestFitValue(conv::kPhypCPV);
       migradmin_true        = migradmin1;
       migraderr_true        = migraderr1;
       fparm1->ShowOnlyBestFitValuesInNextPrint();
       cout << *fparm1;
       //br_fparm->SetAddress(&fparm1);
       if(store_bestfit_spectrum)
         br_bf_array->SetAddress(&pdf_array1);
     }
 
     //
     // Get the 2+N fit results
     //

     cout << endl;
     cout << "////////////////////////////////////////////////////////////////" << endl;
     cout << endl;
     cout << "deltacp: Min" << endl;
     cout << endl;

     fparm3->SetInitValue(conv::kPhypMixing23, input_mixing23);
     fparm3->SetInitValue(conv::kPhypDm23sq,   input_dmsq);
     fparm3->SetInitValue(conv::kPhypMixing13, input_mixing13);
     fparm3->SetInitValue(conv::kPhypCPV,  input_deltacp);
     fparm3->SetNuisance(conv::kPhypMixing13,true);//reactor!

     // stay in [-Pi;+Pi] region
     if(fparm3->InitValue(conv::kPhypCPV) < -conv::kpi)
       fparm3->SetInitValue(conv::kPhypCPV,fparm3->InitValue(conv::kPhypCPV)+2*conv::kpi);
     else if(fparm3->InitValue(conv::kPhypCPV) > conv::kpi)
       fparm3->SetInitValue(conv::kPhypCPV,fparm3->InitValue(conv::kPhypCPV)-2*conv::kpi);

     fparm3->ResetCurrToInitValues(); // reset osc params
     
     cout << *fparm3 << endl;

     // Get fitter
     FitInputs * fit_inputs3 = new FitInputs(skobs,skexp,fparm3);
     FitterABC * fitter3 = util::GetFitter(fit_inputs3,fit_option);

     // do the fit using the true osc param as a seed
     fitter3->SetTol(tolerance);
     fitter3->DoTheFit();
     fparm3->ResetCurrToInitValues();
     migradmin3 = fitter3->FitStatusCode("migrad_min");
     migraderr3 = fitter3->FitStatusCode("migrad_err");
     cout<<"Migrad min and err codes for fit 3:   "<<migradmin3<<"    "<<migraderr3<<endl;
     bestfit_chisq_seed1  = fparm3->GetBestFitStatistic();

     // Opposite hierarchy 

     cout << endl;
     cout << "////////////////////////////////////////////////////////////////" << endl;
     cout << endl;
     cout << "deltacp: Min - opposite hierarchy" << endl;
     cout << endl;

     if(hierarchy == "normal")
       Opt::Instance()->gMassHierarchyAssumed = util::HierarchyFromString("inverted");
     else Opt::Instance()->gMassHierarchyAssumed = util::HierarchyFromString("normal");
     
     fparm4->SetInitValue(conv::kPhypMixing23, input_mixing23);
     fparm4->SetInitValue(conv::kPhypDm23sq,   input_dmsq);
     fparm4->SetInitValue(conv::kPhypMixing13, input_mixing13);
     fparm4->SetInitValue(conv::kPhypCPV,  input_deltacp);
     fparm4->SetNuisance(conv::kPhypMixing13,true);//reactor!
     
     // stay in [-Pi;+Pi] region
     if(fparm4->InitValue(conv::kPhypCPV) < -conv::kpi)
       fparm4->SetInitValue(conv::kPhypCPV,fparm4->InitValue(conv::kPhypCPV)+2*conv::kpi);
     else if(fparm4->InitValue(conv::kPhypCPV) > conv::kpi)
       fparm4->SetInitValue(conv::kPhypCPV,fparm4->InitValue(conv::kPhypCPV)-2*conv::kpi);
     
     fparm4->ResetCurrToInitValues(); // reset osc params
     
     cout << *fparm4 << endl;
     
     // Get fitter
     FitInputs * fit_inputs4 = new FitInputs(skobs,skexp_mh,fparm4);
     FitterABC * fitter4 = util::GetFitter(fit_inputs4,fit_option);

     // do the fit using the true osc param as a seed
     fitter4->SetTol(tolerance);
     fitter4->DoTheFit();
     fparm4->ResetCurrToInitValues();
     migradmin4 = fitter4->FitStatusCode("migrad_min");
     migraderr4 = fitter4->FitStatusCode("migrad_err");
     cout<<"Migrad min and err codes for fit 4:   "<<migradmin4<<"    "<<migraderr4<<endl;
     bestfit_chisq_seed2  = fparm4->GetBestFitStatistic();

     bestfit_chisq_seed3 = 0;
     bestfit_chisq_seed4 = 0;

     // Save spectra (data and best-fit prediction)
     if(store_toy_expts) {
       skobs->Write(dat_array4, "obs");
     }
     if(store_bestfit_spectrum) {
       skexp->Write(pdf_array4, "pdf");
     }

     //Do the fit again at the mirror point

     cout << endl;
     cout << "////////////////////////////////////////////////////////////////" << endl;
     cout << endl;
     cout << "deltacp: Min - mirror" << endl;
     cout << endl;

     migradmin5 = -1;
     migradmin6 = -1;
     if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsqTheta) {
       //get the correct intial values (best fit points from first fit)
       //Opt::Instance()->gMassHierarchyAssumed = util::HierarchyFromString(hierarchy);
       fparm5->CopyBestFitToInitValues(*fparm3);
       //set the mixing23 seed to be the mirror point of the initial best fit
       sinsqtheta23_bestfit = fparm3->BestFitValue(conv::kPhypMixing23);
       sinsqtheta13_bestfit = fparm3->BestFitValue(conv::kPhypMixing13);
       centre = 1 / (2 * (1 - sinsqtheta13_bestfit));
       reseed = sinsqtheta23_bestfit;
       if     (reseed < centre)
         reseed = centre + TMath::Abs(centre - sinsqtheta23_bestfit);
       else if(reseed > centre)
         reseed = centre - TMath::Abs(centre - sinsqtheta23_bestfit);
       fparm5->SetInitValue(conv::kPhypMixing23, reseed);
       fparm5->SetNuisance(conv::kPhypMixing13,true);//reactor!

       // stay in [-Pi;+Pi] region
       if(fparm5->InitValue(conv::kPhypCPV) < -conv::kpi)
	 fparm5->SetInitValue(conv::kPhypCPV,fparm5->InitValue(conv::kPhypCPV)+2*conv::kpi);
       else if(fparm5->InitValue(conv::kPhypCPV) > conv::kpi)
	 fparm5->SetInitValue(conv::kPhypCPV,fparm5->InitValue(conv::kPhypCPV)-2*conv::kpi);

       fparm5->ResetCurrToInitValues(); // reset osc params
     
       cout << *fparm5 << endl;

       // Get fitter
       FitInputs * fit_inputs5 = new FitInputs(skobs,skexp,fparm5);
       FitterABC * fitter5 = util::GetFitter(fit_inputs5,fit_option);
       fitter5->SetTol(tolerance);
       fitter5->DoTheFit();
       bestfit_chisq_seed3 = fparm5->GetBestFitStatistic();
       migradmin5 = fitter5->FitStatusCode("migrad_min");
       migraderr5 = fitter5->FitStatusCode("migrad_err");
       cout<<"Migrad min and err codes for fit 5:   "<<migradmin5<<"    "<<migraderr5<<endl;
       // Save spectra (best-fit prediction)
       if(store_bestfit_spectrum) {
         skexp->Write(pdf_array5, "pdf");
       }

       
       // Opposite hierarchy

       cout << endl;
       cout << "////////////////////////////////////////////////////////////////" << endl;
       cout << endl;
       cout << "deltacp: Min - mirror point - opposite hierarchy" << endl;
       cout << endl;
       
       if(hierarchy == "normal")
     	 Opt::Instance()->gMassHierarchyAssumed = util::HierarchyFromString("inverted");
       else Opt::Instance()->gMassHierarchyAssumed = util::HierarchyFromString("normal");

       fparm6->CopyBestFitToInitValues(*fparm4);
       //set the mixing23 seed to be the mirror point of the initial best fit
       sinsqtheta23_bestfit = fparm4->BestFitValue(conv::kPhypMixing23);
       sinsqtheta13_bestfit = fparm4->BestFitValue(conv::kPhypMixing13);
       centre = 1 / (2 * (1 - sinsqtheta13_bestfit));
       reseed = sinsqtheta23_bestfit;
       if     (reseed < centre)
         reseed = centre + TMath::Abs(centre - sinsqtheta23_bestfit);
       else if(reseed > centre)
         reseed = centre - TMath::Abs(centre - sinsqtheta23_bestfit);
       fparm6->SetInitValue(conv::kPhypMixing23, reseed);
       fparm6->SetNuisance(conv::kPhypMixing13,true);//reactor!
       
      // stay in [-Pi;+Pi] region
      if(fparm6->InitValue(conv::kPhypCPV) < -conv::kpi)
	 fparm6->SetInitValue(conv::kPhypCPV,fparm6->InitValue(conv::kPhypCPV)+2*conv::kpi);
       else if(fparm6->InitValue(conv::kPhypCPV) > conv::kpi)
	 fparm6->SetInitValue(conv::kPhypCPV,fparm6->InitValue(conv::kPhypCPV)-2*conv::kpi);
       
       fparm6->ResetCurrToInitValues(); // reset osc params
       
       cout << *fparm6 << endl;

       // Get fitter
       FitInputs * fit_inputs6 = new FitInputs(skobs,skexp_mh,fparm6);
       FitterABC * fitter6 = util::GetFitter(fit_inputs6,fit_option);
       fitter6->SetTol(tolerance);
       fitter6->DoTheFit();
       bestfit_chisq_seed4 = fparm6->GetBestFitStatistic();
       migradmin6 = fitter6->FitStatusCode("migrad_min");
       migraderr6 = fitter6->FitStatusCode("migrad_err");
       cout<<"Migrad min and err codes for fit 6:   "<<migradmin6<<"    "<<migraderr6<<endl;
       // Save spectra (best-fit prediction)
       if(store_bestfit_spectrum) {
         skexp->Write(pdf_array6, "pdf");
       }
     }
      

     ////////////////////////////// delta_cp + Pi //////////////////////////////////
 
     //
     // Get the 2+N fit results
     //

     cout << endl;
     cout << "////////////////////////////////////////////////////////////////" << endl;
     cout << endl;
     cout << "deltacp+Pi: Min" << endl;
     cout << endl;

     fparm3_dcp->SetInitValue(conv::kPhypMixing23, input_mixing23);
     fparm3_dcp->SetInitValue(conv::kPhypDm23sq,   input_dmsq);
     fparm3_dcp->SetInitValue(conv::kPhypMixing13, input_mixing13);
     fparm3_dcp->SetNuisance(conv::kPhypMixing13,true);//reactor!

     // +Pi --> other dcp minimum
     if(input_deltacp<=0.) fparm3_dcp->SetInitValue(conv::kPhypCPV,  input_deltacp+conv::kpi);
     else                  fparm3_dcp->SetInitValue(conv::kPhypCPV,  input_deltacp-conv::kpi);

     // stay in [-Pi;+Pi] region
     if(fparm3_dcp->InitValue(conv::kPhypCPV) < -conv::kpi)
       fparm3_dcp->SetInitValue(conv::kPhypCPV,fparm3_dcp->InitValue(conv::kPhypCPV)+2*conv::kpi);
     else if(fparm3_dcp->InitValue(conv::kPhypCPV) > conv::kpi)
       fparm3_dcp->SetInitValue(conv::kPhypCPV,fparm3_dcp->InitValue(conv::kPhypCPV)-2*conv::kpi);

     fparm3_dcp->ResetCurrToInitValues(); // reset osc params

     cout << *fparm3_dcp << endl;

     // Get fitter
     FitInputs * fit_inputs3_dcp = new FitInputs(skobs,skexp,fparm3_dcp);
     FitterABC * fitter3_dcp = util::GetFitter(fit_inputs3_dcp,fit_option);

     // do the fit using the true osc param as a seed
     fitter3_dcp->SetTol(tolerance);
     fitter3_dcp->DoTheFit();
     fparm3_dcp->ResetCurrToInitValues();
     migradmin3_dcp = fitter3_dcp->FitStatusCode("migrad_min");
     migraderr3_dcp = fitter3_dcp->FitStatusCode("migrad_err");
     cout<<"Migrad min and err codes for fit 3 (deltacp+Pi):   "<<migradmin3_dcp<<"    "<<migraderr3_dcp<<endl;
     bestfit_chisq_seed1_dcp  = fparm3_dcp->GetBestFitStatistic();

 
     // Opposite hierarchy 

     cout << endl;
     cout << "////////////////////////////////////////////////////////////////" << endl;
     cout << endl;
     cout << "deltacp+Pi: Min - opposite hierarchy" << endl;
     cout << endl;

     if(hierarchy == "normal")
       Opt::Instance()->gMassHierarchyAssumed = util::HierarchyFromString("inverted");
     else Opt::Instance()->gMassHierarchyAssumed = util::HierarchyFromString("normal");
     
     fparm4_dcp->SetInitValue(conv::kPhypMixing23, input_mixing23);
     fparm4_dcp->SetInitValue(conv::kPhypDm23sq,   input_dmsq);
     fparm4_dcp->SetInitValue(conv::kPhypMixing13, input_mixing13);
     fparm4_dcp->SetNuisance(conv::kPhypMixing13,true);//reactor!
 
     // +Pi --> other dcp minimum
     if(input_deltacp<=0) fparm4_dcp->SetInitValue(conv::kPhypCPV,  input_deltacp+conv::kpi);
     else                 fparm4_dcp->SetInitValue(conv::kPhypCPV,  input_deltacp-conv::kpi);
    
     // stay in [-Pi;+Pi] region
     if(fparm4_dcp->InitValue(conv::kPhypCPV) < -conv::kpi)
       fparm4_dcp->SetInitValue(conv::kPhypCPV,fparm4_dcp->InitValue(conv::kPhypCPV)+2*conv::kpi);
     else if(fparm4_dcp->InitValue(conv::kPhypCPV) > conv::kpi)
       fparm4_dcp->SetInitValue(conv::kPhypCPV,fparm4_dcp->InitValue(conv::kPhypCPV)-2*conv::kpi);

     fparm4_dcp->ResetCurrToInitValues(); // reset osc params

     cout << *fparm4_dcp << endl;

     // Get fitter
     FitInputs * fit_inputs4_dcp = new FitInputs(skobs,skexp_mh,fparm4_dcp);
     FitterABC * fitter4_dcp = util::GetFitter(fit_inputs4_dcp,fit_option);

     // do the fit using the true osc param as a seed
     fitter4_dcp->SetTol(tolerance);
     fitter4_dcp->DoTheFit();
     fparm4_dcp->ResetCurrToInitValues();
     migradmin4_dcp = fitter4_dcp->FitStatusCode("migrad_min");
     migraderr4_dcp = fitter4_dcp->FitStatusCode("migrad_err");
     cout<<"Migrad min and err codes for fit 4 (deltacp+Pi):   "<<migradmin4_dcp<<"    "<<migraderr4_dcp<<endl;
     bestfit_chisq_seed2_dcp  = fparm4_dcp->GetBestFitStatistic();

     bestfit_chisq_seed3_dcp = 0;
     bestfit_chisq_seed4_dcp = 0;

     // Save spectra (data and best-fit prediction)
     if(store_toy_expts) {
       skobs->Write(dat_array4_dcp, "obs");
     }
     if(store_bestfit_spectrum) {
       skexp->Write(pdf_array4_dcp, "pdf");
     }

     //Do the fit again at the mirror point

     cout << endl;
     cout << "////////////////////////////////////////////////////////////////" << endl;
     cout << endl;
     cout << "deltacp+Pi: Min - mirror" << endl;
     cout << endl;

     migradmin5_dcp = -1;
     migradmin6_dcp = -1;
     if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsqTheta) {
       //get the correct intial values (best fit points from first fit)
       //Opt::Instance()->gMassHierarchyAssumed = util::HierarchyFromString(hierarchy);
       fparm5_dcp->CopyBestFitToInitValues(*fparm3_dcp);
       //set the mixing23 seed to be the mirror point of the initial best fit
       sinsqtheta23_bestfit = fparm3_dcp->BestFitValue(conv::kPhypMixing23);
       sinsqtheta13_bestfit = fparm3_dcp->BestFitValue(conv::kPhypMixing13);
       centre = 1 / (2 * (1 - sinsqtheta13_bestfit));
       reseed = sinsqtheta23_bestfit;
       if     (reseed < centre)
         reseed = centre + TMath::Abs(centre - sinsqtheta23_bestfit);
       else if(reseed > centre)
         reseed = centre - TMath::Abs(centre - sinsqtheta23_bestfit);
       fparm5_dcp->SetInitValue(conv::kPhypMixing23, reseed);
       fparm5_dcp->SetNuisance(conv::kPhypMixing13,true);//reactor!
       
       // stay in [-Pi;+Pi] region
       if(fparm5_dcp->InitValue(conv::kPhypCPV) < -conv::kpi)
	 fparm5_dcp->SetInitValue(conv::kPhypCPV,fparm5_dcp->InitValue(conv::kPhypCPV)+2*conv::kpi);
       else if(fparm5_dcp->InitValue(conv::kPhypCPV) > conv::kpi)
	 fparm5_dcp->SetInitValue(conv::kPhypCPV,fparm5_dcp->InitValue(conv::kPhypCPV)-2*conv::kpi);
       fparm5_dcp->ResetCurrToInitValues(); // reset osc params
       
       cout << *fparm5_dcp << endl;

       // Get fitter
       FitInputs * fit_inputs5_dcp = new FitInputs(skobs,skexp,fparm5_dcp);
       FitterABC * fitter5_dcp = util::GetFitter(fit_inputs5_dcp,fit_option);
       fitter5_dcp->SetTol(tolerance);
       fitter5_dcp->DoTheFit();
       bestfit_chisq_seed3_dcp = fparm5_dcp->GetBestFitStatistic();
       migradmin5_dcp = fitter5_dcp->FitStatusCode("migrad_min");
       migraderr5_dcp = fitter5_dcp->FitStatusCode("migrad_err");
       cout<<"Migrad min and err codes for fit 5:   "<<migradmin5_dcp<<"    "<<migraderr5_dcp<<endl;
       // Save spectra (best-fit prediction)
       if(store_bestfit_spectrum) {
         skexp->Write(pdf_array5_dcp, "pdf");
       }
       
       // Opposite hierarchy
       
       cout << endl;
       cout << "////////////////////////////////////////////////////////////////" << endl;
       cout << endl;
       cout << "deltacp+Pi: Min - mirror point - opposite hierarchy" << endl;
       cout << endl;
       
       if(hierarchy == "normal")
     	 Opt::Instance()->gMassHierarchyAssumed = util::HierarchyFromString("inverted");
       else Opt::Instance()->gMassHierarchyAssumed = util::HierarchyFromString("normal");

       fparm6_dcp->CopyBestFitToInitValues(*fparm4_dcp);
       //set the mixing23 seed to be the mirror point of the initial best fit
       sinsqtheta23_bestfit = fparm4_dcp->BestFitValue(conv::kPhypMixing23);
       sinsqtheta13_bestfit = fparm4_dcp->BestFitValue(conv::kPhypMixing13);
       centre = 1 / (2 * (1 - sinsqtheta13_bestfit));
       reseed = sinsqtheta23_bestfit;
       if     (reseed < centre)
         reseed = centre + TMath::Abs(centre - sinsqtheta23_bestfit);
       else if(reseed > centre)
         reseed = centre - TMath::Abs(centre - sinsqtheta23_bestfit);
       fparm6_dcp->SetInitValue(conv::kPhypMixing23, reseed);
       fparm6_dcp->SetNuisance(conv::kPhypMixing13,true);//reactor!
       
       // stay in [-Pi;+Pi] region
       if(fparm6_dcp->InitValue(conv::kPhypCPV) < -conv::kpi)
	 fparm6_dcp->SetInitValue(conv::kPhypCPV,fparm6_dcp->InitValue(conv::kPhypCPV)+2*conv::kpi);
       else if(fparm6_dcp->InitValue(conv::kPhypCPV) > conv::kpi)
	 fparm6_dcp->SetInitValue(conv::kPhypCPV,fparm6_dcp->InitValue(conv::kPhypCPV)-2*conv::kpi);
       
       fparm6_dcp->ResetCurrToInitValues(); // reset osc params
       
       cout << *fparm6_dcp << endl;

       // Get fitter
       FitInputs * fit_inputs6_dcp = new FitInputs(skobs,skexp_mh,fparm6_dcp);
       FitterABC * fitter6_dcp = util::GetFitter(fit_inputs6_dcp,fit_option);
       fitter6_dcp->SetTol(tolerance);
       fitter6_dcp->DoTheFit();
       bestfit_chisq_seed4_dcp = fparm6_dcp->GetBestFitStatistic();
       migradmin6_dcp = fitter6_dcp->FitStatusCode("migrad_min");
       migraderr6_dcp = fitter6_dcp->FitStatusCode("migrad_err");
       cout<<"Migrad min and err codes for fit 6:   "<<migradmin6_dcp<<"    "<<migraderr6_dcp<<endl;
       // Save spectra (best-fit prediction)
       if(store_bestfit_spectrum) {
         skexp->Write(pdf_array6_dcp, "pdf");
       }
     }
     
      
     //set the "best fit of the four fits" correctly
     double bestfit_chisq_total = 999;
     int index_chisq = 0;
     if (bestfit_chisq_seed1<bestfit_chisq_total)  {
       bestfit_chisq_total = bestfit_chisq_seed1;
       index_chisq = 1;
     }
     if (bestfit_chisq_seed2<bestfit_chisq_total)  {
       bestfit_chisq_total = bestfit_chisq_seed2;
       index_chisq = 2;
     }
     if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsqTheta) {
       if (bestfit_chisq_seed3<bestfit_chisq_total)  {
  	 bestfit_chisq_total = bestfit_chisq_seed3;
  	 index_chisq = 3;
       }
       if (bestfit_chisq_seed4<bestfit_chisq_total)  {
  	 bestfit_chisq_total = bestfit_chisq_seed4;
  	 index_chisq = 4;
       }
     }

     if (bestfit_chisq_seed1_dcp<bestfit_chisq_total)  {
       bestfit_chisq_total = bestfit_chisq_seed1_dcp;
       index_chisq = 5;
     }
     if (bestfit_chisq_seed2_dcp<bestfit_chisq_total)  {
       bestfit_chisq_total = bestfit_chisq_seed2_dcp;
       index_chisq = 6;
     }
     if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsqTheta) {
       if (bestfit_chisq_seed3_dcp<bestfit_chisq_total)  {
  	 bestfit_chisq_total = bestfit_chisq_seed3_dcp;
  	 index_chisq = 7;
       }
       if (bestfit_chisq_seed4_dcp<bestfit_chisq_total)  {
  	 bestfit_chisq_total = bestfit_chisq_seed4_dcp;
  	 index_chisq = 8;
       }
     }

     
     if(index_chisq == 1){
       cout << bestfit_chisq_seed1 << " from the first seed fit" << endl;
       fittries_min = fparm3->GetFitTries();
       fittries_min_oct = fparm5->GetFitTries();
       fittries_min_mh = fparm4->GetFitTries();
       fittries_min_mh_oct = fparm6->GetFitTries();
       fittime_min = fparm3->GetFitTime();
       fittime_min_oct = fparm5->GetFitTime();
       fittime_min_mh = fparm4->GetFitTime();
       fittime_min_mh_oct = fparm6->GetFitTime();
       bestfit_chisq_min    = bestfit_chisq_seed1;
       bestfit_chisq_min_oct    = bestfit_chisq_seed3;
       bestfit_chisq_min_mh    = bestfit_chisq_seed2;
       bestfit_chisq_min_mh_oct = bestfit_chisq_seed4;
       bestfit_mixing23_min = fparm3->BestFitValue(conv::kPhypMixing23);
       bestfit_mixing23_min_oct = fparm5->BestFitValue(conv::kPhypMixing23);
       bestfit_mixing23_min_mh = fparm4->BestFitValue(conv::kPhypMixing23);
       bestfit_mixing23_min_mh_oct = fparm6->BestFitValue(conv::kPhypMixing23);
       bestfit_mixing13_min = fparm3->BestFitValue(conv::kPhypMixing13);
       bestfit_mixing13_min_oct = fparm5->BestFitValue(conv::kPhypMixing13);
       bestfit_mixing13_min_mh = fparm4->BestFitValue(conv::kPhypMixing13);
       bestfit_mixing13_min_mh_oct = fparm6->BestFitValue(conv::kPhypMixing13);
       bestfit_dmsq_min = fparm3->BestFitValue(conv::kPhypDm23sq);
       bestfit_dmsq_min_oct = fparm5->BestFitValue(conv::kPhypDm23sq);
       bestfit_dmsq_min_mh = fparm4->BestFitValue(conv::kPhypDm23sq);
       bestfit_dmsq_min_mh_oct = fparm6->BestFitValue(conv::kPhypDm23sq);
       bestfit_deltacp_min = fparm3->BestFitValue(conv::kPhypCPV);
       bestfit_deltacp_min_oct = fparm5->BestFitValue(conv::kPhypCPV);
       bestfit_deltacp_min_mh = fparm4->BestFitValue(conv::kPhypCPV);
       bestfit_deltacp_min_mh_oct = fparm6->BestFitValue(conv::kPhypCPV);
       migradmin_min        = migradmin1;
       migradmin_min_oct    = migradmin3;
       migradmin_min_mh     = migradmin2;
       migradmin_min_mh_oct = migradmin4;
       migraderr_min        = migraderr1;
       migraderr_min_oct    = migraderr3;
       migraderr_min_mh     = migraderr2;
       migraderr_min_mh_oct = migraderr4;
       fparm1->ShowOnlyBestFitValuesInNextPrint();
       cout << *fparm3;
       //br_fparm->SetAddress(&fparm3);
       //if(store_bestfit_spectrum)
  	 //br_bf_array->SetAddress(&pdf_array1);
     }
     else if(index_chisq == 2){
     cout << bestfit_chisq_seed2 << " from the second seed fit" << endl;
       fittries_min = fparm4->GetFitTries();
       fittries_min_oct = fparm6->GetFitTries();
       fittries_min_mh = fparm3->GetFitTries();
       fittries_min_mh_oct = fparm5->GetFitTries();
       fittime_min = fparm4->GetFitTime();
       fittime_min_oct = fparm6->GetFitTime();
       fittime_min_mh = fparm3->GetFitTime();
       fittime_min_mh_oct = fparm5->GetFitTime();
       bestfit_chisq_min    = bestfit_chisq_seed2;
       bestfit_chisq_min_oct    = bestfit_chisq_seed4;
       bestfit_chisq_min_mh    = bestfit_chisq_seed1;
       bestfit_chisq_min_mh_oct = bestfit_chisq_seed3;
       bestfit_mixing23_min = fparm4->BestFitValue(conv::kPhypMixing23);
       bestfit_mixing23_min_oct = fparm6->BestFitValue(conv::kPhypMixing23);
       bestfit_mixing23_min_mh = fparm3->BestFitValue(conv::kPhypMixing23);
       bestfit_mixing23_min_mh_oct = fparm5->BestFitValue(conv::kPhypMixing23);
       bestfit_mixing13_min = fparm4->BestFitValue(conv::kPhypMixing13);
       bestfit_mixing13_min_oct = fparm6->BestFitValue(conv::kPhypMixing13);
       bestfit_mixing13_min_mh = fparm3->BestFitValue(conv::kPhypMixing13);
       bestfit_mixing13_min_mh_oct = fparm5->BestFitValue(conv::kPhypMixing13);
       bestfit_dmsq_min = fparm4->BestFitValue(conv::kPhypDm23sq);
       bestfit_dmsq_min_oct = fparm6->BestFitValue(conv::kPhypDm23sq);
       bestfit_dmsq_min_mh = fparm3->BestFitValue(conv::kPhypDm23sq);
       bestfit_dmsq_min_mh_oct = fparm5->BestFitValue(conv::kPhypDm23sq);
       bestfit_deltacp_min = fparm4->BestFitValue(conv::kPhypCPV);
       bestfit_deltacp_min_oct = fparm6->BestFitValue(conv::kPhypCPV);
       bestfit_deltacp_min_mh = fparm3->BestFitValue(conv::kPhypCPV);
       bestfit_deltacp_min_mh_oct = fparm5->BestFitValue(conv::kPhypCPV);
       migradmin_min        = migradmin2;
       migradmin_min_oct    = migradmin4;
       migradmin_min_mh     = migradmin1;
       migradmin_min_mh_oct = migradmin3;
       migraderr_min        = migraderr2;
       migraderr_min_oct    = migraderr4;
       migraderr_min_mh     = migraderr1;
       migraderr_min_mh_oct = migraderr3;
       fparm4->ShowOnlyBestFitValuesInNextPrint();
       cout << *fparm4;
       //br_fparm->SetAddress(&fparm4);
       //if(store_bestfit_spectrum)
       //br_bf_array->SetAddress(&pdf_array2);
     }
     else if(index_chisq == 3){
     cout << bestfit_chisq_seed3 << " from the third seed fit" << endl;
       fittries_min = fparm5->GetFitTries();
       fittries_min_oct = fparm3->GetFitTries();
       fittries_min_mh = fparm6->GetFitTries();
       fittries_min_mh_oct = fparm4->GetFitTries();
       fittime_min = fparm5->GetFitTime();
       fittime_min_oct = fparm3->GetFitTime();
       fittime_min_mh = fparm6->GetFitTime();
       fittime_min_mh_oct = fparm4->GetFitTime();
       bestfit_chisq_min    = bestfit_chisq_seed3;
       bestfit_chisq_min_oct    = bestfit_chisq_seed1;
       bestfit_chisq_min_mh    = bestfit_chisq_seed4;
       bestfit_chisq_min_mh_oct = bestfit_chisq_seed2;
       bestfit_mixing23_min = fparm5->BestFitValue(conv::kPhypMixing23);
       bestfit_mixing23_min_oct = fparm3->BestFitValue(conv::kPhypMixing23);
       bestfit_mixing23_min_mh = fparm6->BestFitValue(conv::kPhypMixing23);
       bestfit_mixing23_min_mh_oct = fparm4->BestFitValue(conv::kPhypMixing23);
       bestfit_mixing13_min = fparm5->BestFitValue(conv::kPhypMixing13);
       bestfit_mixing13_min_oct = fparm3->BestFitValue(conv::kPhypMixing13);
       bestfit_mixing13_min_mh = fparm6->BestFitValue(conv::kPhypMixing13);
       bestfit_mixing13_min_mh_oct = fparm4->BestFitValue(conv::kPhypMixing13);
       bestfit_dmsq_min = fparm5->BestFitValue(conv::kPhypDm23sq);
       bestfit_dmsq_min_oct = fparm3->BestFitValue(conv::kPhypDm23sq);
       bestfit_dmsq_min_mh = fparm6->BestFitValue(conv::kPhypDm23sq);
       bestfit_dmsq_min_mh_oct = fparm4->BestFitValue(conv::kPhypDm23sq);
       bestfit_deltacp_min = fparm5->BestFitValue(conv::kPhypCPV);
       bestfit_deltacp_min_oct = fparm3->BestFitValue(conv::kPhypCPV);
       bestfit_deltacp_min_mh = fparm6->BestFitValue(conv::kPhypCPV);
       bestfit_deltacp_min_mh_oct = fparm4->BestFitValue(conv::kPhypCPV);
       migradmin_min        = migradmin3;
       migradmin_min_oct    = migradmin1;
       migradmin_min_mh     = migradmin4;
       migradmin_min_mh_oct = migradmin2;
       migraderr_min        = migraderr3;
       migraderr_min_oct    = migraderr1;
       migraderr_min_mh     = migraderr4;
       migraderr_min_mh_oct = migraderr2;
       fparm5->ShowOnlyBestFitValuesInNextPrint();
       cout << *fparm5;
       //br_fparm->SetAddress(&fparm5);
       //if(store_bestfit_spectrum)
       //br_bf_array->SetAddress(&pdf_array3);
     }
     else if(index_chisq == 4){
     cout << bestfit_chisq_seed4 << " from the forth seed fit" << endl;
       fittries_min = fparm6->GetFitTries();
       fittries_min_oct = fparm4->GetFitTries();
       fittries_min_mh = fparm5->GetFitTries();
       fittries_min_mh_oct = fparm3->GetFitTries();
       fittime_min = fparm6->GetFitTime();
       fittime_min_oct = fparm4->GetFitTime();
       fittime_min_mh = fparm5->GetFitTime();
       fittime_min_mh_oct = fparm3->GetFitTime();
       bestfit_chisq_min    = bestfit_chisq_seed4;
       bestfit_chisq_min_oct    = bestfit_chisq_seed2;
       bestfit_chisq_min_mh    = bestfit_chisq_seed3;
       bestfit_chisq_min_mh_oct = bestfit_chisq_seed1;
       bestfit_mixing23_min = fparm6->BestFitValue(conv::kPhypMixing23);
       bestfit_mixing23_min_oct = fparm4->BestFitValue(conv::kPhypMixing23);
       bestfit_mixing23_min_mh = fparm5->BestFitValue(conv::kPhypMixing23);
       bestfit_mixing23_min_mh_oct = fparm3->BestFitValue(conv::kPhypMixing23);
       bestfit_mixing13_min = fparm6->BestFitValue(conv::kPhypMixing13);
       bestfit_mixing13_min_oct = fparm4->BestFitValue(conv::kPhypMixing13);
       bestfit_mixing13_min_mh = fparm5->BestFitValue(conv::kPhypMixing13);
       bestfit_mixing13_min_mh_oct = fparm3->BestFitValue(conv::kPhypMixing13);
       bestfit_dmsq_min = fparm6->BestFitValue(conv::kPhypDm23sq);
       bestfit_dmsq_min_oct = fparm4->BestFitValue(conv::kPhypDm23sq);
       bestfit_dmsq_min_mh = fparm5->BestFitValue(conv::kPhypDm23sq);
       bestfit_dmsq_min_mh_oct = fparm3->BestFitValue(conv::kPhypDm23sq);
       bestfit_deltacp_min = fparm6->BestFitValue(conv::kPhypCPV);
       bestfit_deltacp_min_oct = fparm4->BestFitValue(conv::kPhypCPV);
       bestfit_deltacp_min_mh = fparm5->BestFitValue(conv::kPhypCPV);
       bestfit_deltacp_min_mh_oct = fparm3->BestFitValue(conv::kPhypCPV);
       migradmin_min        = migradmin4;
       migradmin_min_oct    = migradmin2;
       migradmin_min_mh     = migradmin3;
       migradmin_min_mh_oct = migradmin1;
       migraderr_min        = migraderr4;
       migraderr_min_oct    = migraderr2;
       migraderr_min_mh     = migraderr3;
       migraderr_min_mh_oct = migraderr1;
       fparm6->ShowOnlyBestFitValuesInNextPrint();
       cout << *fparm6;
       //br_fparm->SetAddress(&fparm5);
       //if(store_bestfit_spectrum)
       //br_bf_array->SetAddress(&pdf_array4);
     }
     
     else if(index_chisq == 5){
       cout << bestfit_chisq_seed1_dcp << " from the fifth seed fit (dcp+pi)" << endl;
       fittries_min = fparm3_dcp->GetFitTries();
       fittries_min_oct = fparm5_dcp->GetFitTries();
       fittries_min_mh = fparm4_dcp->GetFitTries();
       fittries_min_mh_oct = fparm6_dcp->GetFitTries();
       fittime_min = fparm3_dcp->GetFitTime();
       fittime_min_oct = fparm5_dcp->GetFitTime();
       fittime_min_mh = fparm4_dcp->GetFitTime();
       fittime_min_mh_oct = fparm6_dcp->GetFitTime();
       bestfit_chisq_min    = bestfit_chisq_seed1_dcp;
       bestfit_chisq_min_oct    = bestfit_chisq_seed3_dcp;
       bestfit_chisq_min_mh    = bestfit_chisq_seed2_dcp;
       bestfit_chisq_min_mh_oct = bestfit_chisq_seed4_dcp;
       bestfit_mixing23_min = fparm3_dcp->BestFitValue(conv::kPhypMixing23);
       bestfit_mixing23_min_oct = fparm5_dcp->BestFitValue(conv::kPhypMixing23);
       bestfit_mixing23_min_mh = fparm4_dcp->BestFitValue(conv::kPhypMixing23);
       bestfit_mixing23_min_mh_oct = fparm6_dcp->BestFitValue(conv::kPhypMixing23);
       bestfit_mixing13_min = fparm3_dcp->BestFitValue(conv::kPhypMixing13);
       bestfit_mixing13_min_oct = fparm5_dcp->BestFitValue(conv::kPhypMixing13);
       bestfit_mixing13_min_mh = fparm4_dcp->BestFitValue(conv::kPhypMixing13);
       bestfit_mixing13_min_mh_oct = fparm6_dcp->BestFitValue(conv::kPhypMixing13);
       bestfit_dmsq_min = fparm3_dcp->BestFitValue(conv::kPhypDm23sq);
       bestfit_dmsq_min_oct = fparm5_dcp->BestFitValue(conv::kPhypDm23sq);
       bestfit_dmsq_min_mh = fparm4_dcp->BestFitValue(conv::kPhypDm23sq);
       bestfit_dmsq_min_mh_oct = fparm6_dcp->BestFitValue(conv::kPhypDm23sq);
       bestfit_deltacp_min = fparm3_dcp->BestFitValue(conv::kPhypCPV);
       bestfit_deltacp_min_oct = fparm5_dcp->BestFitValue(conv::kPhypCPV);
       bestfit_deltacp_min_mh = fparm4_dcp->BestFitValue(conv::kPhypCPV);
       bestfit_deltacp_min_mh_oct = fparm6_dcp->BestFitValue(conv::kPhypCPV);
       migradmin_min        = migradmin1_dcp;
       migradmin_min_oct    = migradmin3_dcp;
       migradmin_min_mh     = migradmin2_dcp;
       migradmin_min_mh_oct = migradmin4_dcp;
       migraderr_min        = migraderr1_dcp;
       migraderr_min_oct    = migraderr3_dcp;
       migraderr_min_mh     = migraderr2_dcp;
       migraderr_min_mh_oct = migraderr4_dcp;
       fparm3_dcp->ShowOnlyBestFitValuesInNextPrint();
       cout << *fparm3_dcp;
       //br_fparm->SetAddress(&fparm3_dcp);
       //if(store_bestfit_spectrum)
       //br_bf_array->SetAddress(&pdf_array1_dcp);
     }
     else if(index_chisq == 6){
     cout << bestfit_chisq_seed2_dcp << " from the sixth seed fit (dcp+pi)" << endl;
       fittries_min = fparm4_dcp->GetFitTries();
       fittries_min_oct = fparm6_dcp->GetFitTries();
       fittries_min_mh = fparm3_dcp->GetFitTries();
       fittries_min_mh_oct = fparm5_dcp->GetFitTries();
       fittime_min = fparm4_dcp->GetFitTime();
       fittime_min_oct = fparm6_dcp->GetFitTime();
       fittime_min_mh = fparm3_dcp->GetFitTime();
       fittime_min_mh_oct = fparm5_dcp->GetFitTime();
       bestfit_chisq_min    = bestfit_chisq_seed2_dcp;
       bestfit_chisq_min_oct    = bestfit_chisq_seed4_dcp;
       bestfit_chisq_min_mh    = bestfit_chisq_seed1_dcp;
       bestfit_chisq_min_mh_oct = bestfit_chisq_seed3_dcp;
       bestfit_mixing23_min = fparm4_dcp->BestFitValue(conv::kPhypMixing23);
       bestfit_mixing23_min_oct = fparm6_dcp->BestFitValue(conv::kPhypMixing23);
       bestfit_mixing23_min_mh = fparm3_dcp->BestFitValue(conv::kPhypMixing23);
       bestfit_mixing23_min_mh_oct = fparm5_dcp->BestFitValue(conv::kPhypMixing23);
       bestfit_mixing13_min = fparm4_dcp->BestFitValue(conv::kPhypMixing13);
       bestfit_mixing13_min_oct = fparm6_dcp->BestFitValue(conv::kPhypMixing13);
       bestfit_mixing13_min_mh = fparm3_dcp->BestFitValue(conv::kPhypMixing13);
       bestfit_mixing13_min_mh_oct = fparm5_dcp->BestFitValue(conv::kPhypMixing13);
       bestfit_dmsq_min = fparm4_dcp->BestFitValue(conv::kPhypDm23sq);
       bestfit_dmsq_min_oct = fparm6_dcp->BestFitValue(conv::kPhypDm23sq);
       bestfit_dmsq_min_mh = fparm3_dcp->BestFitValue(conv::kPhypDm23sq);
       bestfit_dmsq_min_mh_oct = fparm5_dcp->BestFitValue(conv::kPhypDm23sq);
       bestfit_deltacp_min = fparm4_dcp->BestFitValue(conv::kPhypCPV);
       bestfit_deltacp_min_oct = fparm6_dcp->BestFitValue(conv::kPhypCPV);
       bestfit_deltacp_min_mh = fparm3_dcp->BestFitValue(conv::kPhypCPV);
       bestfit_deltacp_min_mh_oct = fparm5_dcp->BestFitValue(conv::kPhypCPV);
       migradmin_min        = migradmin2_dcp;
       migradmin_min_oct    = migradmin4_dcp;
       migradmin_min_mh     = migradmin1_dcp;
       migradmin_min_mh_oct = migradmin3_dcp;
       migraderr_min        = migraderr2_dcp;
       migraderr_min_oct    = migraderr4_dcp;
       migraderr_min_mh     = migraderr1_dcp;
       migraderr_min_mh_oct = migraderr3_dcp;
       fparm4_dcp->ShowOnlyBestFitValuesInNextPrint();
       cout << *fparm4_dcp;
       //br_fparm->SetAddress(&fparm4_dcp);
       //if(store_bestfit_spectrum)
       //br_bf_array->SetAddress(&pdf_array2_dcp);
     }
     else if(index_chisq == 7){
     cout << bestfit_chisq_seed3_dcp << " from the seventh seed fit (dcp+pi)" << endl;
       fittries_min = fparm5_dcp->GetFitTries();
       fittries_min_oct = fparm3_dcp->GetFitTries();
       fittries_min_mh = fparm6_dcp->GetFitTries();
       fittries_min_mh_oct = fparm4_dcp->GetFitTries();
       fittime_min = fparm5_dcp->GetFitTime();
       fittime_min_oct = fparm3_dcp->GetFitTime();
       fittime_min_mh = fparm6_dcp->GetFitTime();
       fittime_min_mh_oct = fparm4_dcp->GetFitTime();
       bestfit_chisq_min    = bestfit_chisq_seed3_dcp;
       bestfit_chisq_min_oct    = bestfit_chisq_seed1_dcp;
       bestfit_chisq_min_mh    = bestfit_chisq_seed4_dcp;
       bestfit_chisq_min_mh_oct = bestfit_chisq_seed2_dcp;
       bestfit_mixing23_min = fparm5_dcp->BestFitValue(conv::kPhypMixing23);
       bestfit_mixing23_min_oct = fparm3_dcp->BestFitValue(conv::kPhypMixing23);
       bestfit_mixing23_min_mh = fparm6_dcp->BestFitValue(conv::kPhypMixing23);
       bestfit_mixing23_min_mh_oct = fparm4_dcp->BestFitValue(conv::kPhypMixing23);
       bestfit_mixing13_min = fparm5_dcp->BestFitValue(conv::kPhypMixing13);
       bestfit_mixing13_min_oct = fparm3_dcp->BestFitValue(conv::kPhypMixing13);
       bestfit_mixing13_min_mh = fparm6_dcp->BestFitValue(conv::kPhypMixing13);
       bestfit_mixing13_min_mh_oct = fparm4_dcp->BestFitValue(conv::kPhypMixing13);
       bestfit_dmsq_min = fparm5_dcp->BestFitValue(conv::kPhypDm23sq);
       bestfit_dmsq_min_oct = fparm3_dcp->BestFitValue(conv::kPhypDm23sq);
       bestfit_dmsq_min_mh = fparm6_dcp->BestFitValue(conv::kPhypDm23sq);
       bestfit_dmsq_min_mh_oct = fparm4_dcp->BestFitValue(conv::kPhypDm23sq);
       bestfit_deltacp_min = fparm5_dcp->BestFitValue(conv::kPhypCPV);
       bestfit_deltacp_min_oct = fparm3_dcp->BestFitValue(conv::kPhypCPV);
       bestfit_deltacp_min_mh = fparm6_dcp->BestFitValue(conv::kPhypCPV);
       bestfit_deltacp_min_mh_oct = fparm4_dcp->BestFitValue(conv::kPhypCPV);
       migradmin_min        = migradmin3_dcp;
       migradmin_min_oct    = migradmin1_dcp;
       migradmin_min_mh     = migradmin4_dcp;
       migradmin_min_mh_oct = migradmin2_dcp;
       migraderr_min        = migraderr3_dcp;
       migraderr_min_oct    = migraderr1_dcp;
       migraderr_min_mh     = migraderr4_dcp;
       migraderr_min_mh_oct = migraderr2_dcp;
       fparm5_dcp->ShowOnlyBestFitValuesInNextPrint();
       cout << *fparm5_dcp;
       //br_fparm->SetAddress(&fparm5_dcp);
       //if(store_bestfit_spectrum)
       //br_bf_array->SetAddress(&pdf_array3);
     }
     else if(index_chisq == 8){
     cout << bestfit_chisq_seed4_dcp << " from the eighth seed fit (dcp+pi)" << endl;
       fittries_min = fparm6_dcp->GetFitTries();
       fittries_min_oct = fparm4_dcp->GetFitTries();
       fittries_min_mh = fparm5_dcp->GetFitTries();
       fittries_min_mh_oct = fparm3_dcp->GetFitTries();
       fittime_min = fparm6_dcp->GetFitTime();
       fittime_min_oct = fparm4_dcp->GetFitTime();
       fittime_min_mh = fparm5_dcp->GetFitTime();
       fittime_min_mh_oct = fparm3_dcp->GetFitTime();
       bestfit_chisq_min    = bestfit_chisq_seed4_dcp;
       bestfit_chisq_min_oct    = bestfit_chisq_seed2_dcp;
       bestfit_chisq_min_mh    = bestfit_chisq_seed3_dcp;
       bestfit_chisq_min_mh_oct = bestfit_chisq_seed1_dcp;
       bestfit_mixing23_min = fparm6_dcp->BestFitValue(conv::kPhypMixing23);
       bestfit_mixing23_min_oct = fparm4_dcp->BestFitValue(conv::kPhypMixing23);
       bestfit_mixing23_min_mh = fparm5_dcp->BestFitValue(conv::kPhypMixing23);
       bestfit_mixing23_min_mh_oct = fparm3_dcp->BestFitValue(conv::kPhypMixing23);
       bestfit_mixing13_min = fparm6_dcp->BestFitValue(conv::kPhypMixing13);
       bestfit_mixing13_min_oct = fparm4_dcp->BestFitValue(conv::kPhypMixing13);
       bestfit_mixing13_min_mh = fparm5_dcp->BestFitValue(conv::kPhypMixing13);
       bestfit_mixing13_min_mh_oct = fparm3_dcp->BestFitValue(conv::kPhypMixing13);
       bestfit_dmsq_min = fparm6_dcp->BestFitValue(conv::kPhypDm23sq);
       bestfit_dmsq_min_oct = fparm4_dcp->BestFitValue(conv::kPhypDm23sq);
       bestfit_dmsq_min_mh = fparm5_dcp->BestFitValue(conv::kPhypDm23sq);
       bestfit_dmsq_min_mh_oct = fparm3_dcp->BestFitValue(conv::kPhypDm23sq);
       bestfit_deltacp_min = fparm6_dcp->BestFitValue(conv::kPhypCPV);
       bestfit_deltacp_min_oct = fparm4_dcp->BestFitValue(conv::kPhypCPV);
       bestfit_deltacp_min_mh = fparm5_dcp->BestFitValue(conv::kPhypCPV);
       bestfit_deltacp_min_mh_oct = fparm3_dcp->BestFitValue(conv::kPhypCPV);
       migradmin_min        = migradmin4_dcp;
       migradmin_min_oct    = migradmin2_dcp;
       migradmin_min_mh     = migradmin3_dcp;
       migradmin_min_mh_oct = migradmin1_dcp;
       migraderr_min        = migraderr4_dcp;
       migraderr_min_oct    = migraderr2_dcp;
       migraderr_min_mh     = migraderr3_dcp;
       migraderr_min_mh_oct = migraderr1_dcp;
       fparm6_dcp->ShowOnlyBestFitValuesInNextPrint();
       cout << *fparm6_dcp;
       //br_fparm->SetAddress(&fparm6_dcp);
       //if(store_bestfit_spectrum)
       //br_bf_array->SetAddress(&pdf_array4);
     }
    
     cout<<endl;

     cout << "True:                                " << endl 
	  << "    dcp_init 1 = " << fparm1->InitValue(conv::kPhypCPV) << endl
	  << "    dcp_bf 1 = " << fparm1->BestFitValue(conv::kPhypCPV) << endl
	  << "    bf chisq 1 = " << fparm1->GetBestFitStatistic() << endl
	  << endl;
     cout << "True - mirror:                                " << endl 
	  << "    dcp_init 2 = " << fparm2->InitValue(conv::kPhypCPV) << endl
	  << "    dcp_bf 2 = " << fparm2->BestFitValue(conv::kPhypCPV) << endl
	  << "    bf chisq 2 = " << fparm2->GetBestFitStatistic() << endl
	  << endl;
     cout << "Min:                                " << endl 
	  << "    dcp_init 1 = " << fparm3->InitValue(conv::kPhypCPV) << endl
	  << "    dcp_bf 1 = " << fparm3->BestFitValue(conv::kPhypCPV) << endl
	  << "    bf chisq 1 = " << bestfit_chisq_seed1 << endl
	  << endl;
     cout << "Min - opp hierarchy:                " << endl
	  << "    dcp_init 2 = " << fparm4->InitValue(conv::kPhypCPV) << endl
	  << "    dcp_bf 2 = " << fparm4->BestFitValue(conv::kPhypCPV) << endl
	  << "    bf chisq 2 = " << bestfit_chisq_seed2 << endl
	  << endl;
     cout << "Min - mirror:                       " << endl
	  << "    dcp_init 3 = " << fparm5->InitValue(conv::kPhypCPV) << endl
	  << "    dcp_bf 3 = " << fparm5->BestFitValue(conv::kPhypCPV) << endl
	  << "    bf chisq 3 = " << bestfit_chisq_seed3 << endl
	  << endl;
     cout << "Min - opp hierarchy - mirror:       " << endl
	  << "    dcp_init 4 = " << fparm6->InitValue(conv::kPhypCPV) << endl
	  << "    dcp_bf 4 = " << fparm6->BestFitValue(conv::kPhypCPV) << endl
	  << "    bf chisq 4 = " << bestfit_chisq_seed4 << endl
	  << endl;
     cout << "Min - +Pi:                          " << endl
	  << "    dcp_init 5 = " << fparm3_dcp->InitValue(conv::kPhypCPV) << endl
	  << "    dcp_bf 5 = " << fparm3_dcp->BestFitValue(conv::kPhypCPV) << endl
	  << "    bf chisq 1 (dcp+pi) = " << bestfit_chisq_seed1_dcp << endl
	  << endl;
     cout << "Min - +Pi - opp hierarchy:          " << endl
	  << "    dcp_init 6 = " << fparm4_dcp->InitValue(conv::kPhypCPV) << endl
	  << "    dcp_bf 6 = " << fparm4_dcp->BestFitValue(conv::kPhypCPV) << endl
	  << "    bf chisq 2 (dcp+pi) = " << bestfit_chisq_seed2_dcp << endl
	  << endl;
     cout << "Min - +Pi - mirror:                 " << endl
	  << "    dcp_init 7 = " << fparm5_dcp->InitValue(conv::kPhypCPV) << endl
	  << "    dcp_bf 7 = " << fparm5_dcp->BestFitValue(conv::kPhypCPV) << endl
	  << "    bf chisq 3 (dcp+pi) = " << bestfit_chisq_seed3_dcp << endl
	  << endl;
     cout << "Min - +Pi - opp hierarchy - mirror: " << endl
	  << "    dcp_init 6 = " << fparm6_dcp->InitValue(conv::kPhypCPV) << endl
	  << "    dcp_bf 8 = " << fparm6_dcp->BestFitValue(conv::kPhypCPV) << endl
	  << "    bf chisq 4 (dcp+pi) = " << bestfit_chisq_seed4_dcp << endl
	  << endl;
     cout << endl;
     
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
     delete fit_inputs3_dcp;
     delete fitter3_dcp;
     delete fit_inputs4_dcp;
     delete fitter4_dcp;
     if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsqTheta) {
       delete fit_inputs2;
       delete fitter2;
       delete fit_inputs5;
       delete fitter5;
       delete fit_inputs6;
       delete fitter6;
       delete fit_inputs5_dcp;
       delete fitter5_dcp;
       delete fit_inputs6_dcp;
       delete fitter6_dcp;
     }
     
     }//iexperiment

  // Store ntuple with best-fit results
  f.cd();
  fitnt->Write();
  osc_param_space->Write();
  theta13_param_space->Write();
  theta23_param_space->Write();
  dmsq_param_space->Write();
  f.Close();
}
