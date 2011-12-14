//
// ------------------------------------------------------------------------------------------------
// Calculate chi^2 = f(13-mixing, 23-mixing, atmospheric Dm^2 mass splitting, delta_cp)
//
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

#include <iostream>

using namespace std;

void build_chisq_surface_4d
(
  int          dataset_id            = 0,                                                  /* number in <data_set_%d> XML tag, see inputs/data.xml */
  const char * data_file             = "inputs/data.xml",                                  /* input SK data */
  const char * config_file           = "analysis_configuration.xml",                       /* analysis configuration file */
  const char * fit_param_list        = "nuosc3f",                                          /* fit parameter list  */

  // Cuts to be applied in the oscillation parameter space considered in this job.
  // -1 means that no cut will be applied.
  int input_dcp_min_bin     = -1,
  int input_dcp_max_bin     = -1,
  int input_mix13_min_bin   = -1,
  int input_mix13_max_bin   = -1,
  int input_mix23_min_bin   = -1,
  int input_mix23_max_bin   = -1,
  int input_dm23sq_min_bin  = -1,
  int input_dm23sq_max_bin  = -1  
)
{
  using namespace std;
  using namespace t2kana;

  // Analysis configuration
  const char * config  = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);

  // Data file
  const char * data_file_fullpath = Form("%s/%s", gSystem->Getenv("T2KOSCANA"), data_file);

  // Global options
  Opt::Instance()->gOscMixingParamChoice = conv::kNuMixSinsqTheta;
  Opt::Instance()->gMassHierarchyAssumed = conv::kNmhNormal;

  // Parameter range and binning
  double dcp_min     = -1 * conv::kpi;
  double dcp_max     = conv::kpi;
  double dcp_step    = 2*conv::kpi/90.;
  double mix13_min   = 0.0;
  double mix13_max   = 0.5;
  double mix13_step  = 0.005;
  double mix23_min   = 0.3;
  double mix23_max   = 0.7;
  double mix23_step  = 0.005;
  double dm23sq_min  = 1.0E-3;
  double dm23sq_max  = 5.0E-3;
  double dm23sq_step = 1E-4;

  // Make sure that min and max are at the bin centre and not at the low/high edge of the first/last bin
  int dcp_nbins    = 1 + (int) ((dcp_max    - dcp_min   ) / dcp_step   );
  int mix13_nbins  = 1 + (int) ((mix13_max  - mix13_min ) / mix13_step );
  int mix23_nbins  = 1 + (int) ((mix23_max  - mix23_min ) / mix23_step );
  int dm23sq_nbins = 1 + (int) ((dm23sq_max - dm23sq_min) / dm23sq_step);

  // Apply user cuts in the bin range
  int dcp_min_bin     = (input_dcp_min_bin    == -1) ? 1            : input_dcp_min_bin    ;
  int dcp_max_bin     = (input_dcp_max_bin    == -1) ? dcp_nbins    : input_dcp_max_bin    ;
  int mix13_min_bin   = (input_mix13_min_bin  == -1) ? 1            : input_mix13_min_bin  ;
  int mix13_max_bin   = (input_mix13_max_bin  == -1) ? mix13_nbins  : input_mix13_max_bin  ;
  int mix23_min_bin   = (input_mix23_min_bin  == -1) ? 1            : input_mix23_min_bin  ;
  int mix23_max_bin   = (input_mix23_max_bin  == -1) ? mix23_nbins  : input_mix23_max_bin  ;
  int dm23sq_min_bin  = (input_dm23sq_min_bin == -1) ? 1            : input_dm23sq_min_bin ;
  int dm23sq_max_bin  = (input_dm23sq_max_bin == -1) ? dm23sq_nbins : input_dm23sq_max_bin ;

  // Print settings
  cout << "[Job inputs]" << endl;
  cout << " - Reading analysis configuration from: " << config << endl;
  cout << " - Reading dataset: " << dataset_id << " from: " << data_file_fullpath << endl;
  cout << " - Fit parameter list: " << fit_param_list << endl;
  cout << "- Oscillation parameter space binning" << endl;
  cout << "   - " << util::AsString(conv::kPhypMixing13) << " : " << mix13_nbins 
       << " bins in (" << mix13_min << ", " << mix13_max << ")" << endl;
  cout << "   - " << util::AsString(conv::kPhypMixing23) << " : " << mix23_nbins 
       << " bins in (" << mix23_min << ", " << mix23_max << ")" << endl;
  cout << "   - " << util::AsString(conv::kPhypDeltaCP) << " : " << dcp_nbins 
       << " bins in (" << dcp_min << ", " << dcp_max << ")" << endl;
  cout << "   - " << util::AsString(conv::kPhypDm23sq) << " : " << dm23sq_nbins 
       << " bins in (" << dm23sq_min << ", " << dm23sq_max << ") eV^2/c^4" << endl;
  cout << "- Oscillation parameter space considered for current job (incl. input cuts)" << endl;
  cout << "   - " << util::AsString(conv::kPhypMixing13) << " : " 
       << " considering bins [" << mix13_min_bin << ", " << mix13_max_bin << "]" << endl;
  cout << "   - " << util::AsString(conv::kPhypMixing23) << " : " 
       << " considering bins [" << mix23_min_bin << ", " << mix23_max_bin << "]" << endl;
  cout << "   - " << util::AsString(conv::kPhypDeltaCP) << " : " 
       << " considering bins [" << dcp_min_bin << ", " << dcp_max_bin << "]" << endl;
  cout << "   - " << util::AsString(conv::kPhypDm23sq) << " : " 
       << " considering bins [" << dm23sq_min_bin << ", " << dm23sq_max_bin << "]" << endl;
  cout << *Opt::Instance() << endl;

  // Output
  TNtupleD * ntchisq = new TNtupleD("ntchisq","","dcp:mix13:mix23:dm23sq:chisq:chisq1Rmu:chisq1Re");

  // Get superk spectrum predictor & configure
  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read(config);

  // Get superk spectrum observation
  SuperkObservation * skobs = new SuperkObservation;
  skobs->Read(data_file, dataset_id);

  // Set osc param list
  FitParamList * fparm = FitParamList::NamedParameterList(fit_param_list);
  fparm->Update(conv::kPhypMixing12, 0.313);    // fix solar osc params
  fparm->Update(conv::kPhypDm12sq,   7.58E-5);  // fix solar osc params

  for(int idcp = dcp_min_bin; idcp <= dcp_max_bin; idcp++) {
    double dcp = dcp_min + (idcp-1) * dcp_step;
    fparm->Update(conv::kPhypDeltaCP, dcp);

    for(int imix13 = mix13_min_bin; imix13 <= mix13_max_bin; imix13++) {
      double mix13 = mix13_min + (imix13-1) * mix13_step;
      fparm->Update(conv::kPhypMixing13, mix13);

      for(int imix23 = mix23_min_bin; imix23 <= mix23_max_bin; imix23++) {
          double mix23 = mix23_min + (imix23-1) * mix23_step; 
          fparm->Update(conv::kPhypMixing23, mix23);

          for(int idmsq = dm23sq_min_bin; idmsq <= dm23sq_max_bin; idmsq++) {
             double dm23sq = dm23sq_min + (idmsq-1) * dm23sq_step;
             fparm->Update(conv::kPhypDm23sq, dm23sq);

             // Generate prediction
             skexp->Generate(*fparm);

             // Get fit inputs object (chisq is calculated there)
             FitInputs fit_inputs(skobs,skexp,fparm);
             double chisq     = fit_inputs.Chisquare("poisson");
             double chisq1Rmu = fit_inputs.Chisquare("poisson", conv::kSKsFCFV1RingMuLike);
             double chisq1Re  = fit_inputs.Chisquare("poisson", conv::kSKsFCFV1RingELike);

             // Report
             fparm->ShowOnlyCurrValuesInNextPrint();
             cout << *fparm;
             cout << "chisq        = " << chisq      << endl;
             cout << "chisq [1Rmu] = " << chisq1Rmu  << endl;
             cout << "chisq [1Re ] = " << chisq1Re   << endl;

             // Save
             ntchisq->Fill(dcp,mix13,mix23,dm23sq,chisq,chisq1Rmu,chisq1Re);
          }
       }
    }
  }

  // Write-out
  const char * filename = 
    Form("chisq_dcp[%d,%d],mix13[%d,%d],mix23[%d,%d],dm23sq[%d,%d].root",
          dcp_min_bin,    dcp_max_bin,
          mix13_min_bin,  mix13_max_bin,
          mix23_min_bin,  mix23_max_bin,
          dm23sq_min_bin, dm23sq_max_bin
        );
  TFile f(filename,"recreate");
  ntchisq->Write();
  f.Close();

}
