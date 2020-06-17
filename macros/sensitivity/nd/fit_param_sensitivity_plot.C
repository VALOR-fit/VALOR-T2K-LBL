//
// ------------------------------------------------------------------------------------------------
// Processes fit data produced by fit_param_sensitivity.C
// For each toy-MC experiment fit, it extracts the best-fit statistic as function of pv_fit (see
// comments in fit_param_sensitivity.C) and then averages over all toy-MC experiments.
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

#include <algorithm>
#include <multimap>
#include <vector>


TGraph * get_sensitivity     (const char * pname, const char * filename, double & ptrue, double & pmin, double & pmax);
TF1 *    get_sensitivity_fit (TGraph * sensitivity, double pmin, double pmax);
TGraph * get_prior           (const char * pname, double pmin, double pmax, double & err);

void fit_param_sensitivity_plot
(
  const char * pname,                        /* parameter name */
  bool         fit_sensitivity_curve,        /* plot averaged 'Delta\chi^2 directly, or show fit with quadratic function? */
  const char * filename1,                    /* run_fit.C output ROOT file with fit results (use wildcards if there are many) */
  const char * label1       = "Sensitivity", /*  */
  const char * filename2    = 0,             /*  */
  const char * label2       = 0,             /*  */
  bool         custom_range = false,         /* by default plots will be shown in range of fit data, but if sensitivity is fitted with pol2 it can be extrapolated - define custom range if needed */
  double       pmin_in      =  9999,
  double       pmax_in      = -9999,
  const char * config       = "analysis/dune_tdr.xml"
)
{
  using namespace valor;

  // Read the config file to populate the ParameterPool
  std::string config_path = util::CompleteConfigPath(config);
  CompoundFitElements * elements = new CompoundFitElements();
  elements->Read(config_path.c_str());

  double epsilon = 1E-5;

  double ptrue = -999;
  double pmin  = +999;
  double pmax  = -999;
  double err   = -999;

  //
  // get sensitivity from first file
  //
  double ptrue1 = -999;
  double pmin1  =  999;
  double pmax1  = -999;
  TGraph * sensitivity1 = get_sensitivity(pname,filename1,ptrue1,pmin1,pmax1);  
  sensitivity1->SetLineColor(kBlack);
  sensitivity1->SetLineStyle(kSolid);
  sensitivity1->SetMarkerStyle(20);
  sensitivity1->GetXaxis()->SetTitle(pname);
  sensitivity1->GetYaxis()->SetTitle("#delta(-2lnL)");
  TF1 * sensitivity1_fit_func = 0;
  if(fit_sensitivity_curve) {
    double pmin_fit = (custom_range) ? pmin_in : pmin1;
    double pmax_fit = (custom_range) ? pmax_in : pmax1;
    sensitivity1_fit_func = get_sensitivity_fit (sensitivity1,pmin_fit,pmax_fit);
    sensitivity1_fit_func->SetLineColor(kBlack);
    sensitivity1_fit_func->SetLineStyle(kSolid);
    sensitivity1_fit_func->GetXaxis()->SetTitle(pname);
    sensitivity1_fit_func->GetYaxis()->SetTitle("#delta(-2lnL)");
  }
  ptrue = ptrue1;
  pmin  = pmin1;
  pmax  = pmax1;

  //
  // get alternative sensitivity from second file (optional)
  //
  double ptrue2 = -999;
  double pmin2  =  999;
  double pmax2  = -999;
  TGraph * sensitivity2 = 0;
  TF1 * sensitivity2_fit_func = 0;
  if(filename2) {
     sensitivity2 = get_sensitivity(pname,filename2,ptrue2,pmin2,pmax2);  
     assert(TMath::Abs(ptrue1-ptrue2) < epsilon);
     sensitivity2->SetLineColor(kBlack);
     sensitivity2->SetLineStyle(kDashed);
     sensitivity2->SetMarkerStyle(20);
     sensitivity2->GetXaxis()->SetTitle(pname);
     sensitivity2->GetYaxis()->SetTitle("#delta(-2lnL)");
     if(fit_sensitivity_curve) {
       double pmin_fit = (custom_range) ? pmin_in : pmin2;
       double pmax_fit = (custom_range) ? pmax_in : pmax2;
       sensitivity2_fit_func = get_sensitivity_fit (sensitivity2,pmin_fit,pmax_fit);
       sensitivity2_fit_func->SetLineColor(kBlack);
       sensitivity2_fit_func->SetLineStyle(kDashed);
       sensitivity2_fit_func->GetXaxis()->SetTitle(pname);
       sensitivity2_fit_func->GetYaxis()->SetTitle("#delta(-2lnL)");
     }
  }
  if(custom_range) {
    pmin = pmin_in;
    pmax = pmax_in;
  } else {
    pmin = TMath::Min(pmin,pmin2);
    pmax = TMath::Max(pmax,pmax2);
  }

  //
  // get prior
  //
  TGraph * prior = get_prior(pname,pmin,pmax,err);
  prior->SetLineColor(kRed);
  prior->SetLineWidth(2);

  //
  // Draw
  //

  TCanvas * c = new TCanvas("c","",20,20,600,600);
  if(fit_sensitivity_curve) {
    sensitivity1_fit_func->Draw();
    if(sensitivity2_fit_func) sensitivity2_fit_func->Draw("same");
    prior->Draw("lsame");
    TLegend * leg = new TLegend(0.4, 0.6, 0.8, 0.9);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetHeader(Form("True %s = %.2f",pname,ptrue));
    leg->AddEntry(sensitivity1_fit_func,label1,"L");
    if(sensitivity2_fit_func) leg->AddEntry("sensitivity2_fit_func",label2,"L");
    leg->AddEntry(prior, Form("Prior constraint (%.1f\%\%)",err), "L");
    leg->Draw();
  } else {
    sensitivity1->Draw("alp");
    prior->Draw("lsame");
    TLegend * leg = new TLegend(0.4, 0.6, 0.8, 0.9);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetHeader(Form("True %s = %.2f",pname,ptrue));
    leg->AddEntry(sensitivity1,label1,"L");
    if(sensitivity2) leg->AddEntry("sensitivity2",label2,"L");
    leg->AddEntry(prior, Form("Prior constraint (%.1f\%\%)",err), "L");
    leg->Draw();
  }
  TLine *  l68cl = new TLine (pmin,1.00,pmax,1.00);
  TLine *  l90cl = new TLine (pmin,2.71,pmax,2.71);
  TLine *  l99cl = new TLine (pmin,6.63,pmax,6.63);
  TLatex * t68cl = new TLatex(1.01*pmin,1.01*1.00,"68\% C.L.");
  TLatex * t90cl = new TLatex(1.01*pmin,1.01*2.71,"90\% C.L.");
  TLatex * t99cl = new TLatex(1.01*pmin,1.01*6.63,"99\% C.L.");
  l68cl->SetLineStyle(kDashed);
  l90cl->SetLineStyle(kDashed);
  l99cl->SetLineStyle(kDashed);
  t68cl->SetTextSize(0.04);
  t90cl->SetTextSize(0.04);
  t99cl->SetTextSize(0.04);
  l68cl->Draw();
  l90cl->Draw();
  l99cl->Draw();
  t68cl->Draw();
  t90cl->Draw();
  t99cl->Draw();
  c->Update();

  cout << "Done!" << endl;
}

// ------------------------------------------------------------------------------------------------
TGraph * get_sensitivity( 
  const char * pname, const char * filename, double & p_true_out, double & p_min_out, double & p_max_out)
{
  conv::PhysParamId_t p = util::NameRegistry::Instance()->PhysParamFromString(pname);
  assert(p != conv::kPhypUndefined);
  
  cout << "Plotting sensitivity to parameter: " << pname << " using toy-MC fit results stored in: " << filename << endl;
  
  // Open input files
  TChain * fitnt = new TChain("fitnt");
  fitnt->Add(filename);

  // Get tree
  // Tree structure below is what is written out by run_fit.C
  // Similar (or adapt trivially) for other macros

  PhysParamList * true_param_array     = 0; // input true params
  FitParamList *  best_fit_param_array = 0; // best-fit params
//TObjArray *     best_fit_model_array = 0; // an array of various best-fit model predictions for the current toy-MC experiment 
//TObjArray *     data_array           = 0; // an array of all fitted data distributions for the current toy-MC experiment

  fitnt->SetBranchAddress("true_params",        &true_param_array);
  fitnt->SetBranchAddress("best_fit_params",    &best_fit_param_array);
//fitnt->SetBranchAddress("best_fit_model",     &best_fit_model_array);
//fitnt->SetBranchAddress("data",               &data_array);

  // Get entries
  Long_t ntoymc = fitnt->GetEntries();
  cout << "Found " << ntoymc << " toy-MC experiment fits in the input tree" << endl;

  // Use a multimap to store "param value (fixed in fit)" --> "best-fit statistic" 
  // mappings for all toy-MC experiments included in the input ROOT file.
  multimap<double, double> inpfitdata;

  bool   first       = true;
  double p_true_save = -999999999;
  double epsilon     = 1E-5;

  // Loop over entries (toy-MC experiments) in the tree
  for(Long_t itoymc = 0; itoymc < ntoymc; itoymc++) {
    fitnt->GetEntry(itoymc);  
    double p_true = true_param_array->Value(p); 
    if(first) {
      p_true_save = p_true;
      first = false;
    } else {
      if( TMath::Abs(p_true-p_true_save) > epsilon ) {
        cout << "All toy-MC experiments should have been generated with the same true valye of parameter: " << pname << endl;
        cout << "Current value = " << p_true << ", but a value of " << p_true_save << " was seen earlier" << endl;
        //exit(1);
      }
    }
    double p_bf     = best_fit_param_array->BestFitValue(p); // this param was fixed, it was the other N-1 that were fitted
    double chisq_bf = best_fit_param_array->GetBestFitStatistic();
    cout << "> Toy-MC: #" << itoymc << " : true param = " << p_true 
         << ", best-fit param = " << p_bf << ", best-fit \chi^2 = " << chisq_bf << endl;
    inpfitdata.insert(pair<double,double>(p_bf,chisq_bf));
  }//itoy

  //
  // Now process the input fit data
  //

  // Get a sorted vector of all unique paramerer values
  cout << "Getting unique param values (keys)" << endl;
  vector<double> p_vec;
  for(multimap<double, double>::iterator it_mmp = inpfitdata.begin(); it_mmp != inpfitdata.end(); ++it_mmp) {
     double p_value = it_mmp->first;
   //cout << "Parameter value: " << p_value << endl;
     if(count(p_vec.begin(), p_vec.end(), p_value) == 0) {
         p_vec.push_back(p_value);
         cout << "Adding new value (" << p_value << ")" << endl;
     }
  }
  sort(p_vec.begin(), p_vec.end());

  // For each true value, get all toy-MC results and average the best-fit statistic
  cout << "Averaging best-fit statistic for each true param value" << endl;
  map<double,double> chisq_avg_mp; // true param value -> average chi^2
  double chisq_avg_min = 9999999;
  for(vector<double>::iterator it_v = p_vec.begin(); it_v != p_vec.end(); ++it_v) {
    double p_value = *it_v;
    cout << "Averaging best-fit statistic for true value = " << p_value << endl;
    cout << "There are " << inpfitdata.count(p_value) << " toy-MC experiments with that true value" << endl;
    double chisq_sum = 0;
    double n = 0;   
    for(multimap<double, double>::iterator it_mmp2 = inpfitdata.begin(); it_mmp2 != inpfitdata.end(); ++it_mmp2) {
      double p_value_curr = it_mmp2->first;
      if(TMath::Abs(p_value-p_value_curr) > epsilon ) {
          continue;
      }
      double chisq = it_mmp2->second;
      cout << "  toy-MC #" << (int)n  << ", true value = " << p_value_curr << " --> chi^2 = " << chisq << endl;
      n = n+1;
      chisq_sum = chisq_sum + chisq;
    }
    double chisq_avg = chisq_sum/n;
    cout << "Average chisq (true param value = " << p_value << ") =" << chisq_avg << endl;
    chisq_avg_min = TMath::Min(chisq_avg_min, chisq_avg);
    chisq_avg_mp.insert(pair<double,double>(p_value,chisq_avg));
  }//unique keys

  // Subtract the minimum best-fit statistic
  map<double,double> delta_chisq_avg_mp;  // true param value -> \delta {average chi^2}
  for(map<double,double>::iterator it_mp = chisq_avg_mp.begin(); it_mp != chisq_avg_mp.end(); ++it_mp) {
    double p_value = it_mp->first;
    double chisq_avg = it_mp->second;
    double delta_chisq_avg = TMath::Max(0., chisq_avg-chisq_avg_min);
    delta_chisq_avg_mp.insert(pair<double,double>(p_value,delta_chisq_avg));
  }

  // Print
  for(map<double,double>::iterator it_mp = delta_chisq_avg_mp.begin(); it_mp != delta_chisq_avg_mp.end(); ++it_mp) {
    cout << "*** True value = " << it_mp->first << " ==> \delta{average \chi^2} = " << it_mp->second << endl;
  }

  // Format as TGraphs
  int pidx = util::ArrIdx::Instance()->FitParamIdx(p);
  double sig = conv::kNDFitParamAPrioriFracError[pidx];
  int np = delta_chisq_avg_mp.size();
  double * p_value_array = new double[np];
  double * dchisq_array = new double[np];
  double * dchisq_penalty_array = new double[np];
  int ip = 0;
  for(map<double,double>::iterator it_mp = delta_chisq_avg_mp.begin(); it_mp != delta_chisq_avg_mp.end(); ++it_mp) {
    p_value_array[ip] = it_mp->first;
    dchisq_array[ip] = it_mp->second;
    dchisq_penalty_array[ip] =  TMath::Power((p_value_array[ip]-1.)/sig, 2.);
    ip++;
  }

  TGraph * sensitivity = new TGraph(np,p_value_array,dchisq_array);

  p_true_out = p_true_save;
  p_min_out  = p_vec[0];
  p_max_out  = p_vec[p_vec.size()-1];
  
  delete fitnt;

  return sensitivity;
}
// ------------------------------------------------------------------------------------------------
TF1 * get_sensitivity_fit (TGraph * sensitivity, double pmin, double pmax)
{
  cout << "Fitting sensitivity with a 2nd order polynomial" << endl;
  sensitivity->Fit("pol2","","",pmin,pmax);
  TF1 * sensitivity_fit_func = sensitivity->GetFunction("pol2");
  double min = sensitivity_fit_func->GetMinimum(pmin,pmax);
  cout << "fit func min = " << min << endl;
  sensitivity_fit_func->SetParameter(0, sensitivity_fit_func->GetParameter(0)-min);
  return sensitivity_fit_func;
}
// ------------------------------------------------------------------------------------------------
TGraph * get_prior (const char * pname, double pmin, double pmax, double & err_out)
{
  conv::PhysParamId_t p = util::NameRegistry::Instance()->PhysParamFromString(pname);
  assert(p != conv::kPhypUndefined);
  
  cout << "Plotting prior for parameter: " << pname << " in range: [" << pmin << ", " << pmax << "]" << endl;

  int pidx = util::ArrIdx::Instance()->FitParamIdx(p);
  double sig = conv::kNDFitParamAPrioriFracError[pidx];

  const int np = 100;
  const double dp = (pmax-pmin)/(np-1);

  double pvalue_array [np];
  double prior_array  [np];

  for(int i=0; i<np; i++) {
     pvalue_array[i] = pmin + i*dp;
     prior_array[i]  = TMath::Power((pvalue_array[i]-1.)/sig, 2.);
  }

  TGraph * prior = new TGraph(np,pvalue_array,prior_array);

  err_out = 100*sig;

  return prior;
}
// ------------------------------------------------------------------------------------------------
