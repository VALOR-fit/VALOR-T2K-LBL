//
// ------------------------------------------------------------------------------------------------
// Calculate the critical value of the fit statistic (using the FC method) 
// for confidence limits involving a single parameter, 
// using the fit results in the input file (generated from coverage_1dlim.C jobs)
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

#include <algorithm>
#include <vector>

void calc_critical_fit_statistic
(
  const char * filename,
  const char * lim_param_name,  /* specify the param name (see UtlParam.cxx) */
  const char * config
)
{
  using namespace valor;

  // Read the config file to populate the ParameterPool
  std::string config_path = util::CompleteConfigPath(config);
  CompoundFitElements * elements = new CompoundFitElements();
  elements->Read(config_path.c_str());

  conv::PhysParamId_t p = util::NameRegistry::Instance()->PhysParamFromString(lim_param_name);
  assert(p != conv::kPhypUndefined);
  
  cout << "Calculating critical values of fit statistic for parameter: " << lim_param_name << endl;
  cout << "Using toy-MC fit results stored in: " << filename << endl;
  
  // Open input file
  TFile f(filename, "read");

  // Get tree
  // Tree structure below is what is written out by run_fit.C
  // Similar (or adapt trivially) for other macros

  double          fit_statistic_true   = 0;
  PhysParamList * true_param_array     = 0; // input true params
  FitParamList *  best_fit_param_array = 0; // best-fit params
//TObjArray *     best_fit_model_array = 0; // an array of various best-fit model predictions for the current toy-MC experiment 
//TObjArray *     data_array           = 0; // an array of all fitted data distributions for the current toy-MC experiment

  TTree * fitnt = (TTree*)f.Get("fitnt");
  fitnt->SetBranchAddress("fit_statistic_true", &fit_statistic_true);
  fitnt->SetBranchAddress("true_params",        &true_param_array);
  fitnt->SetBranchAddress("best_fit_params",    &best_fit_param_array);
//fitnt->SetBranchAddress("best_fit_model",     &best_fit_model_array);
//fitnt->SetBranchAddress("data",               &data_array);

  Long_t ntoymc = fitnt->GetEntries();

  // Get unique values of parameter p used for toy-MC generation
  cout << "Getting unique true param values" << endl;
  vector<double> p_true_vec;
  for(Long_t itoymc = 0; itoymc < ntoymc; itoymc++) {
    fitnt->GetEntry(itoymc);
    double p_true = true_param_array->Value(p); 
    if(count(p_true_vec.begin(), p_true_vec.end(), p_true) == 0) {
         p_true_vec.push_back(p_true);
       //cout << "Adding new true value (" << p_true << ")" << endl;
    }
  } //itoy
  // Sort and print-out
  sort(p_true_vec.begin(), p_true_vec.end());
  int npv = p_true_vec.size();
  cout << "There are " << npv << " true values for parameter: " << lim_param_name << endl;
  for(vector<double>::const_iterator it = p_true_vec.begin(); it != p_true_vec.end(); ++it) {
    cout << " > " << *it << endl;
  }
  
  // loop over toy-MC and calculate/write out the critical value of the fit statistic
  // as function of the parameter value
  const double epsilon = 1E-8;
  double * x   = new double[npv];
  double * y68 = new double[npv];
  double * y90 = new double[npv];
  double * y95 = new double[npv];
  int ipv=0;
  for(vector<double>::const_iterator it = p_true_vec.begin(); it != p_true_vec.end(); ++it) {
     double p_true = *it;
     cout << "Calculating critical fit statistic values for true parameter value: " << p_true << endl;
     vector<double> dfitstat_vec;
     for(Long_t itoymc = 0; itoymc < ntoymc; itoymc++) {
        fitnt->GetEntry(itoymc);
        double p_true_toy = true_param_array->Value(p); 
      //cout << "True parameter value for current toy-MC: " << p_true_toy << endl;
        if (TMath::Abs(p_true-p_true_toy) > epsilon) continue;
        double fit_statistic_best = best_fit_param_array->GetBestFitStatistic();
        double delta = fit_statistic_true - fit_statistic_best;
        cout << "toy-mc #" << itoymc << " fit statistic ==> "
             << "true = " << fit_statistic_true 
             << ", best fit = " << fit_statistic_best
             << " (\delta{fit statistic} = " << delta << ")" << endl;
        dfitstat_vec.push_back(delta);
     }//itoy
     // sort
     sort(dfitstat_vec.begin(), dfitstat_vec.end());
     double nfits = dfitstat_vec.size();
     int i68 = (int) (0.68*nfits);
     int i90 = (int) (0.90*nfits);
     int i95 = (int) (0.95*nfits);
     cout << "Indices (nfits = [0 - " << nfits << "]): 68%CL -> " << i68 << ", 90%CL -> " << i90 << ", 95%CL -> " << i95 << endl;
     x  [ipv] = p_true;
     y68[ipv] = dfitstat_vec[i68];
     y90[ipv] = dfitstat_vec[i90];
     y95[ipv] = dfitstat_vec[i95];
     cout << "Critical values (p = " << p_true << "): 68%CL -> " << y68[ipv] << ", 90%CL -> " << y90[ipv] << ", 95%CL -> " << y95[ipv] << endl;
     ipv++;
  }//it
  TGraph * dfitstat_crit_68cl = new TGraph(npv,x,y68);
  TGraph * dfitstat_crit_90cl = new TGraph(npv,x,y90);
  TGraph * dfitstat_crit_95cl = new TGraph(npv,x,y95);
  TFile fout("./crit.root","recreate");
  dfitstat_crit_68cl->Write("dfitstat_crit_68cl");
  dfitstat_crit_90cl->Write("dfitstat_crit_90cl");
  dfitstat_crit_95cl->Write("dfitstat_crit_95cl");
  fout.Close();
}
