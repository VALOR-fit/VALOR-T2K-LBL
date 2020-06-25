//
// ------------------------------------------------------------------------------------------------
// The macro studies how the value of the fit statistic changes as function a single fit 
// (flux or cross-section) parameter. The fit statistic is calcuated for fixed values of all
// other fit parameters (i.e. they are not marginalized). For a similar macro that marginalizes
// all other parameters see: fit_param_sensitivity.C
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//
	
void fit_param_scan
(
  const char * tag            = "scan",                                /* tag for output file */
  const char * pname          = 0,                                     /* parameter name */
  double       pmin           = 0.5,                                   /* min parameter value */
  double       pmax           = 1.5,                                   /* max parameter value */
  int          npv            = 21,                                    /* number of points between minimum and maximum */
  int          ntoymc         = 1,                                     /* number of toy-MC experiments */
  bool         stat_fluct     = false,                                 /* include statistical fluctuation in toy MC expts? */
  const char * param_values   = "nominal",                             /* decide how to set the values of params used for toy-MC generation (see UtlParam.cxx) */
  const char * param_list     = "standard",                            /* model parameter list (see UtlParam.cxx) */
  const char * config_file    = "dune_ndfit-2015v1-FGTv1.xml"    /* analysis configuration file */
)
{
  using namespace std;
  using namespace valor;
  using namespace valor::conv;

  // Analysis configuration
  TString config = Form("%s/config/%s", gSystem->Getenv("VALOR_EXPT"), config_file);

  // Print settings
  cout << "\n[Job inputs]" << endl;
  cout << "- Job: " << tag << endl;
  cout << "- Scanning parameter: " << pname << " in range [" << pmin << ", " << pmax << "] using " << npv << " points" << endl;
  cout << "- Reading analysis configuration from: " << config << endl;
  cout << "- Toy-MC generation" << endl;
  cout << "   - Number of toy experiments: " << ntoymc << endl;
  cout << "   - Statistical fluctuations: "    << ((stat_fluct) ? "YES " : "NO") << endl;
  cout << "   - Option for setting parameter values: " << param_values << endl;
  cout << "- Model" << endl;
  cout << "   - Parameter List: " << param_list << endl;
  cout << endl;
  cout << *Rnd::Instance() << endl;

  // Get ND model & configure
  cout << "Getting and configuring fit model" << endl;
  CompoundFitElements * elements = new CompoundFitElements();
  elements->Read(config.Data());

  PhysParamId_t p = util::NameRegistry::Instance()->PhysParamFromString(pname);
  int pidx = util::ArrIdx::Instance()->FitParamIdx(p);

  // Create histogram to plot \Delta{fit statistic} vs parameter
  // Arrange so that the first bin is centered at pmin and the last bin is centered at pmax
  
  //
  double pstep = (pmax - pmin) / (npv-1);
  const char * hndc_name    = Form("nd_constraint_%s", util::NameRegistry::Instance()->ParamAsString(p).c_str());
  const char * hndc_title   = Form("\Delta(lnL) as function of parameter %s in [%.2f,%.2f] range", util::NameRegistry::Instance()->ParamAsString(p).c_str(),pmin,pmax);
  const char * hprior_name  = Form("prior_%s", util::NameRegistry::Instance()->ParamAsString(p).c_str());
  const char * hprior_title = Form("\Delta(lnL) penalty term as function of parameter %s in [%.2f,%.2f] range", util::NameRegistry::Instance()->ParamAsString(p).c_str(),pmin,pmax);
  TH1D * hndc   = new TH1D(hndc_name,  hndc_title,  npv,pmin-0.5*pstep,pmax+0.5*pstep);
  TH1D * hprior = new TH1D(hprior_name,hprior_title,npv,pmin-0.5*pstep,pmax+0.5*pstep);

  for(int ipv = 0; ipv < npv; ipv++) {
    double pcurr = pmin + ipv * pstep;
    cout << "Parameter: " << util::NameRegistry::Instance()->ParamAsString(p) << " -> value = " << pcurr << endl;

    double avg_fit_statistic = 0;

    cout << "Averaging over " << ntoymc << " toy-MC experiments" << endl;
    for(int itoy = 0; itoy < ntoymc; itoy++) {
       cout << "Toy-MC #" << itoy << endl;

       // Get model parameters and set the parameter values according to the
       // input options (ie params can be randomized).
       PhysHypothesis ph = util::GetPhysHypothesis(param_list,param_values);
       //cout << ph;

       // Generate a toy-MC experiment
       FitDataABC * data = new SDData;
       data->Generate(model,ph,stat_fluct);

       // Set input param to current value and calculate the fit statistic
       cout << "Calculating fit-statistic" << endl;
       PhysHypothesis ph2 = util::GetPhysHypothesis(param_list,"nominal");
       ph2.SetCurrValue(p, pcurr);
       //cout << ph2;
       double fit_statistic = data->Statistic(model,ph2,"likelihood-ratio");

       avg_fit_statistic += fit_statistic;

       delete data;
    }//itoy

    avg_fit_statistic /= ntoymc;
    hndc   -> Fill(pcurr, avg_fit_statistic);
    
    double sig = kNDFitParamAPrioriFracError[pidx];
    double prior = TMath::Power((pcurr-1.)/sig, 2.); 
    hprior -> Fill(pcurr, prior);
  }//ipv

  double min = hndc->GetMinimum();
  for(int ibin=1; ibin <= hndc->GetXaxis()->GetNbins(); ibin++) {
    double val = hndc->GetBinContent(ibin);
    double valnew = val-min;
    hndc->SetBinContent(ibin,valnew);
  }

  hndc->GetXaxis()->SetTitle(util::NameRegistry::Instance()->AsStringRooTeX(p).c_str());
  hndc->GetYaxis()->SetTitle("-#Delta(lnL)");
  hndc->SetMarkerSize(1.3);
  hndc->SetMarkerStyle(20);
  hndc->Draw("LP");
  hprior->SetLineColor(kRed);
  hprior->Draw("LSAME");

  const char * outfilename =  Form("%s-%s.root", tag, util::NameRegistry::Instance()->ParamAsString(p).c_str());
  TFile outfile(outfilename, "recreate");
  hndc->Write();
  hprior->Write();
  outfile.Close();

  delete model;
}
