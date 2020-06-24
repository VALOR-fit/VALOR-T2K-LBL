void test_prior_throws_mixing12()
{
  using std::cout;
  using std::endl;
  using namespace valor;
  using namespace valor::conv;

  Env::Instance()->AutoSet("t2xk");

  const int nthrows = 10;

  // Global options
  Opt::Instance()->gMassHierarchyAssumed = kNmhNormal;
  Opt::Instance()->gOscMixingParamChoice = kNuMixSinsqTheta;
  Opt::Instance()->gOscOctantChoice      = kNuOctantTwo;
  Opt::Instance()->gOscDmsqScheme        = kNuDmsqAuto;
  Opt::Instance()->SetDefaultParameterValues();

  util::OscParams::Instance()->SetValues(-999, -999, -999, -conv::kpi, -999, -999, -99,
                                         -999, -999, -999, -conv::kpi, -999, -999,
                                         0.5, false,
                                         -999, -999, -999, -999, -999, -999,
                                         -999, -999, -999, -999, -999, -999
  );

  util::OscParams::Instance()->SetPrior(conv::kPhypMixing12, "gaus", conv::kPrior_Sinsq2th12_Used_Mean, conv::kPrior_Sinsq2th12_Used_RMS, "sinsq2th"); // Convert to sinsqth

  SuperkExpectation * skexp = new SuperkExpectation;
  const char * config = TString::Format("%s/config/%s", gSystem->Getenv("T2KOSCANA"), "analysis_configuration_oa2018_run1-9d_v1_postfit.xml").Data();
  skexp->Read(config);

  // Fit param list for fakedata sets
  FitParamList * fparm  = util::NamedParameterList("nuosc3f+nuisance");
  PhysHypothesis_t hyp = fparm->PhysicsHypothesis();

  // Update physics fitparamlist with chosen value of oscillation parameters
  bool set_init = false;
  set_init = true;
  util::OscParams::Instance()->UpdateFitParamList(fparm,set_init);
  set_init = false;
  util::OscParams::Instance()->UpdateFitParamList(fparm,set_init);

  Rnd::Instance()->SetSeeds(10, 20, 30, 40, 50);  // stat, syst, osc, fitter, other
  for(int ithrow = 0; ithrow < nthrows; ++ithrow) {
    double par1           = util::OscParams::Instance()->params_mixing12[0];
    double par2           = util::OscParams::Instance()->params_mixing12[1];
    std::string priormode = util::OscParams::Instance()->priormode_mixing12;
    std::string mixmode   = util::OscParams::Instance()->mixmode_mixing12;
    mixing12_nh = util::OscParams::Instance()->ExtractParam(conv::kPhypMixing12,priormode,par1,par2,mixmode);
    cout << "Throw #" << ithrow << ": sinsq(theta12) = " << mixing12_nh << endl;
  }
}
