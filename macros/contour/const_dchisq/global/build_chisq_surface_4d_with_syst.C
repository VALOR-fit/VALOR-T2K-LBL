//
// ------------------------------------------------------------------------------------------------
// Fit all nuisance parameters in the 3f analysis in two steps: 
//                first fit only the shape parameters 
//                fix the fitted parameters and release remaining nuisance parameters. 
//                Repeat the fit
// Calculate chi^2 = f(13-mixing, 23-mixing, atmospheric Dm^2 mass splitting, delta_cp)
//
//
// In average a fit per bin takes ~1h. It is better to run one job per bin!
//
// Jelena Ilic
// ------------------------------------------------------------------------------------------------
//

#include <iostream>

using namespace std;

void build_chisq_surface_4d_with_syst
(  
  int          dataset_id            = 0,                                              /* number in <data_set_%d> XML tag, see inputs/data.xml */
  const char * data_file             = "inputs/data.xml",                              /* input SK data */
  const char * config_file           = "analysis_configuration.xml",                   /* analysis configuration file */
  const char * fit_param_list        = "nuosc3f_12fix+nuisance",                       /* fit parameter list  */
  const char * fit_option            = "minuit",                                       /* fitter option (see run_fitter.C) */
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
      double mix13_max   = 0.12;
      double mix13_step  = 0.0012;
      double mix23_min   = 0.3;
      double mix23_max   = 0.7;
      double mix23_step  = 0.005;
      double dm23sq_min  = 1.0E-3;
      double dm23sq_max  = 5.0E-3;
      double dm23sq_step = 1E-4;

      // Make sure that min and max are at the bin centre and not at the low/high edge of the first/last bin
      int dcp_nbins    = 1 + TMath::Nint( (dcp_max    - dcp_min   ) / dcp_step    );
      int mix13_nbins  = 1 + TMath::Nint( (mix13_max  - mix13_min ) / mix13_step  );
      int mix23_nbins  = 1 + TMath::Nint( (mix23_max  - mix23_min ) / mix23_step  );
      int dm23sq_nbins = 1 + TMath::Nint( (dm23sq_max - dm23sq_min) / dm23sq_step );

      // Apply user cuts in the bin range
      int dcp_min_bin     = (input_dcp_min_bin    == -1) ? 1            : input_dcp_min_bin    ;
      int dcp_max_bin     = (input_dcp_max_bin    == -1) ? dcp_nbins    : input_dcp_max_bin    ;
      int mix13_min_bin   = (input_mix13_min_bin  == -1) ? 1            : input_mix13_min_bin  ;
      int mix13_max_bin   = (input_mix13_max_bin  == -1) ? mix13_nbins  : input_mix13_max_bin  ;
      int mix23_min_bin   = (input_mix23_min_bin  == -1) ? 1            : input_mix23_min_bin  ;
      int mix23_max_bin   = (input_mix23_max_bin  == -1) ? mix23_nbins  : input_mix23_max_bin  ;
      int dm23sq_min_bin  = (input_dm23sq_min_bin == -1) ? 1            : input_dm23sq_min_bin ;
      int dm23sq_max_bin  = (input_dm23sq_max_bin == -1) ? dm23sq_nbins : input_dm23sq_max_bin ;

      // Calculate total number of bins considered in this job
      int nbins = (dcp_max_bin    - dcp_min_bin    + 1) *
            (mix13_max_bin  - mix13_min_bin  + 1) *
            (mix23_max_bin  - mix23_min_bin  + 1) *
            (dm23sq_max_bin - dm23sq_min_bin + 1);

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
      cout << "- Current job will consider in total: " << nbins 
            << " oscillation parameter bins" << endl;
      cout << *Opt::Instance() << endl;

      // Output
      TNtupleD * ntchisq = new TNtupleD("ntchisq","","dcp:mix13:mix23:dm23sq:chisq:chisq1Rmu:chisq1Re,chisq_f:chisq1Rmu_f:chisq1Re_f");

      // Get superk spectrum predictor & configure
      SuperkExpectation * skexp = new SuperkExpectation;
      skexp->Read(config);

      // Data file
      const char * data_file_fullpath = Form("%s/%s", gSystem->Getenv("T2KOSCANA"), data_file);
      
      // Get superk spectrum observation
      SuperkObservation * skobs = new SuperkObservation;
      skobs->Read(data_file_fullpath, dataset_id);

      // Set osc param list
      FitParamList * fparm = FitParamList::NamedParameterList(fit_param_list);
      fparm->Update(conv::kPhypMixing12, 0.313);    // fix solar osc params
      fparm->Update(conv::kPhypDm12sq,   7.58E-5);  // fix solar osc params
      fparm->SetInitValue(conv::kPhypMixing12, 0.313);    // fix solar osc params
      fparm->SetInitValue(conv::kPhypDm12sq,   7.58E-5);  // fix solar osc params

      int ibin=0;

      for(int idcp = dcp_min_bin; idcp <= dcp_max_bin; idcp++) {
            double dcp = dcp_min + (idcp-1) * dcp_step;
            fparm->Update(conv::kPhypDeltaCP, dcp);
            fparm->SetInitValue(conv::kPhypDeltaCP, dcp);
            fparm->FixParameter(conv::kPhypDeltaCP);

            for(int imix13 = mix13_min_bin; imix13 <= mix13_max_bin; imix13++) {
                  double mix13 = mix13_min + (imix13-1) * mix13_step;
                  fparm->Update(conv::kPhypMixing13, mix13);
                  fparm->SetInitValue(conv::kPhypMixing13, mix13);
                  fparm->FixParameter(conv::kPhypMixing13);

                  for(int imix23 = mix23_min_bin; imix23 <= mix23_max_bin; imix23++) {
                        double mix23 = mix23_min + (imix23-1) * mix23_step; 
                        fparm->Update(conv::kPhypMixing23, mix23);
                        fparm->SetInitValue(conv::kPhypMixing23, mix23);
                        fparm->FixParameter(conv::kPhypMixing23);

                        for(int idmsq = dm23sq_min_bin; idmsq <= dm23sq_max_bin; idmsq++) {
                              double dm23sq = dm23sq_min + (idmsq-1) * dm23sq_step;
                              fparm->Update(conv::kPhypDm23sq, dm23sq);
                              fparm->SetInitValue(conv::kPhypDm23sq, dm23sq);
                              fparm->FixParameter(conv::kPhypDm23sq);

                              // Set the sample
                              fparm->SetTheSample("numu");
                              cout<<"sample ====="<<fparm->GetTheSample()<<endl; 



                              //fix some nuisance parameters 
                              fparm->FixParameter(conv::kNuisSuperkEff1RmuCCQEoth);
                              fparm->FixParameter(conv::kNuisSuperkEff1RmuCCQErc0u);
                              fparm->FixParameter(conv::kNuisSuperkEff1RmuCCQErc1u);
                              fparm->FixParameter(conv::kNuisSuperkEff1RmuCCQErc2u);
                              fparm->FixParameter(conv::kNuisSuperkEff1RmuCCnonQE);
                              fparm->FixParameter(conv::kNuisSuperkEff1RmuNC);
                              fparm->FixParameter(conv::kNuisSuperkEff1RmuNueCC);
                              fparm->FixParameter(conv::kNuisSuperkEnergyScale);
                              fparm->FixParameter(conv::kNuisNd280Eff);
                              fparm->FixParameter(conv::kNuisFluxNorm);
                              fparm->FixParameter(conv::kNuisXSecCCQE);
                              fparm->FixParameter(conv::kNuisXSecCC1piToCCQE);
                              fparm->FixParameter(conv::kNuisXSecCCothToCCQE);
                              fparm->FixParameter(conv::kNuisXSecNCToCCQE);
                              fparm->FixParameter(conv::kNuisXSecNuclEffect);
                              //additional, only important for nue
                              fparm->FixParameter(conv::kNuisSuperkEff1ReSig);
                              fparm->FixParameter(conv::kNuisSuperkEff1ReBkg);
                              fparm->FixParameter(conv::kNuisXSecCCcohToCCQE);
                              fparm->FixParameter(conv::kNuisXSecNC1pi0ToCCQE);
                              fparm->FixParameter(conv::kNuisXSecNCcohToCCQE);
                              fparm->FixParameter(conv::kNuisXSecNCothToCCQE);
                              fparm->FixParameter(conv::kNuisXSecNuEToNuMu);

                              //get the "true chisq" - before the fit
                              FitInputs *fit_inputs = new FitInputs(skobs,skexp,fparm);
                              double chisq     = fit_inputs->Chisquare(1);//"poisson");
                              double chisq1Rmu = fit_inputs->Chisquare(1, conv::kSKsFCFV1RingMuLike);
                              double chisq1Re  = fit_inputs->Chisquare(1, conv::kSKsFCFV1RingELike);

                              // Report
                              cout << "chisq        = " << chisq      << endl;
                              cout << "chisq [1Rmu] = " << chisq1Rmu  << endl;
                              cout << "chisq [1Re ] = " << chisq1Re   << endl;

                              //fitter
                              Superk1ReNFitter fitter(fit_inputs,fit_option);
                              //do the fit
                              fitter.DoTheFit();

                              //get the best fit results from the previous fit and fix the shape variables
                              fparm->SetInitValue(conv::kNuisFluxShapeNuMu0u ,fparm->CurrValue(conv::kNuisFluxShapeNuMu0u) ); fparm->FixParameter(conv::kNuisFluxShapeNuMu0u); 
                              fparm->SetInitValue(conv::kNuisFluxShapeNuMu1u ,fparm->CurrValue(conv::kNuisFluxShapeNuMu1u) ); fparm->FixParameter(conv::kNuisFluxShapeNuMu1u); 
                              fparm->SetInitValue(conv::kNuisFluxShapeNuMu2u ,fparm->CurrValue(conv::kNuisFluxShapeNuMu2u) ); fparm->FixParameter(conv::kNuisFluxShapeNuMu2u); 
                              fparm->SetInitValue(conv::kNuisFluxShapeNuMu3u ,fparm->CurrValue(conv::kNuisFluxShapeNuMu3u) ); fparm->FixParameter(conv::kNuisFluxShapeNuMu3u); 
                              fparm->SetInitValue(conv::kNuisFluxShapeNuMu4u ,fparm->CurrValue(conv::kNuisFluxShapeNuMu4u) ); fparm->FixParameter(conv::kNuisFluxShapeNuMu4u); 
                              fparm->SetInitValue(conv::kNuisFluxShapeNuMu5u ,fparm->CurrValue(conv::kNuisFluxShapeNuMu5u) ); fparm->FixParameter(conv::kNuisFluxShapeNuMu5u); 
                              fparm->SetInitValue(conv::kNuisFluxShapeNuMu6u ,fparm->CurrValue(conv::kNuisFluxShapeNuMu6u) ); fparm->FixParameter(conv::kNuisFluxShapeNuMu6u); 
                              fparm->SetInitValue(conv::kNuisFluxShapeNuMu7u ,fparm->CurrValue(conv::kNuisFluxShapeNuMu7u) ); fparm->FixParameter(conv::kNuisFluxShapeNuMu7u); 
                              fparm->SetInitValue(conv::kNuisFluxShapeNuMu8u ,fparm->CurrValue(conv::kNuisFluxShapeNuMu8u) ); fparm->FixParameter(conv::kNuisFluxShapeNuMu8u); 
                              fparm->SetInitValue(conv::kNuisFluxShapeNuMu9u ,fparm->CurrValue(conv::kNuisFluxShapeNuMu9u) ); fparm->FixParameter(conv::kNuisFluxShapeNuMu9u); 
                              fparm->SetInitValue(conv::kNuisFluxShapeNuMu10u,fparm->CurrValue(conv::kNuisFluxShapeNuMu10u)); fparm->FixParameter(conv::kNuisFluxShapeNuMu10u);
                              fparm->SetInitValue(conv::kNuisFluxShapeNuMu11u,fparm->CurrValue(conv::kNuisFluxShapeNuMu11u)); fparm->FixParameter(conv::kNuisFluxShapeNuMu11u);
                              fparm->SetInitValue(conv::kNuisFluxShapeNuMu12u,fparm->CurrValue(conv::kNuisFluxShapeNuMu12u)); fparm->FixParameter(conv::kNuisFluxShapeNuMu12u);
                              fparm->SetInitValue(conv::kNuisFluxShapeNuMu13u,fparm->CurrValue(conv::kNuisFluxShapeNuMu13u)); fparm->FixParameter(conv::kNuisFluxShapeNuMu13u);
                              fparm->SetInitValue(conv::kNuisFluxShapeNuMu14u,fparm->CurrValue(conv::kNuisFluxShapeNuMu14u)); fparm->FixParameter(conv::kNuisFluxShapeNuMu14u);
                              fparm->SetInitValue(conv::kNuisFluxShapeNuMu15u,fparm->CurrValue(conv::kNuisFluxShapeNuMu15u)); fparm->FixParameter(conv::kNuisFluxShapeNuMu15u);
                              fparm->SetInitValue(conv::kNuisFluxShapeNuMu16u,fparm->CurrValue(conv::kNuisFluxShapeNuMu16u)); fparm->FixParameter(conv::kNuisFluxShapeNuMu16u);
                              fparm->SetInitValue(conv::kNuisFluxShapeNuMu17u,fparm->CurrValue(conv::kNuisFluxShapeNuMu17u)); fparm->FixParameter(conv::kNuisFluxShapeNuMu17u);
                              fparm->SetInitValue(conv::kNuisFluxShapeNuMu18u,fparm->CurrValue(conv::kNuisFluxShapeNuMu18u)); fparm->FixParameter(conv::kNuisFluxShapeNuMu18u);

                              // release the other parameters and repeat the fit
                              fparm->FreeParameter(conv::kNuisSuperkEff1RmuCCQEoth);
                              fparm->FreeParameter(conv::kNuisSuperkEff1RmuCCQErc0u);
                              fparm->FreeParameter(conv::kNuisSuperkEff1RmuCCQErc1u);
                              fparm->FreeParameter(conv::kNuisSuperkEff1RmuCCQErc2u);
                              fparm->FreeParameter(conv::kNuisSuperkEff1RmuCCnonQE);
                              fparm->FreeParameter(conv::kNuisSuperkEff1RmuNC);
                              fparm->FreeParameter(conv::kNuisSuperkEff1RmuNueCC);
                              fparm->FreeParameter(conv::kNuisSuperkEnergyScale);
                              fparm->FreeParameter(conv::kNuisNd280Eff);
                              fparm->FreeParameter(conv::kNuisFluxNorm);
                              fparm->FreeParameter(conv::kNuisXSecCCQE);
                              fparm->FreeParameter(conv::kNuisXSecCC1piToCCQE);
                              fparm->FreeParameter(conv::kNuisXSecCCothToCCQE);
                              fparm->FreeParameter(conv::kNuisXSecNCToCCQE);
                              fparm->FreeParameter(conv::kNuisXSecNuclEffect);
                              fparm->FreeParameter(conv::kNuisSuperkEff1ReSig);
                              fparm->FreeParameter(conv::kNuisSuperkEff1ReBkg);
                              fparm->FreeParameter(conv::kNuisXSecCCcohToCCQE);
                              fparm->FreeParameter(conv::kNuisXSecNC1pi0ToCCQE);
                              fparm->FreeParameter(conv::kNuisXSecNCcohToCCQE);
                              fparm->FreeParameter(conv::kNuisXSecNCothToCCQE);
                              fparm->FreeParameter(conv::kNuisXSecNuEToNuMu);

                              //Set the sample
                              fparm->SetTheSample(" ");
                              cout<<"sample ====="<<fparm->GetTheSample()<<endl; 
                              //fitter

                              Superk1ReNFitter fitter2(fit_inputs,fit_option);
                              //do the fit
                              fitter2.DoTheFit();

                              double chisq_f     = fit_inputs->Chisquare(1);//"poisson");
                              double chisq1Rmu_f = fit_inputs->Chisquare(1, conv::kSKsFCFV1RingMuLike);
                              double chisq1Re_f  = fit_inputs->Chisquare(1, conv::kSKsFCFV1RingELike);

                              cout << "chisq_f        = " << chisq_f      << endl;
                              cout << "chisq_f [1Rmu] = " << chisq1Rmu_f  << endl;
                              cout << "chisq_f [1Re ] = " << chisq1Re_f   << endl;


                              // Save
                              ntchisq->Fill(dcp,mix13,mix23,dm23sq,chisq,chisq1Rmu,chisq1Re,chisq_f,chisq1Rmu_f,chisq1Re_f);

                              ibin++;
                              double frac_done = 100 * (double)ibin / (double)nbins;
                              cout << " **** " << frac_done << "% done" << endl;
                        }
                  }
            }
      }

      // Write-out
      const char * filename = 
    Form("chisq_dcp%d_mix13%d_mix23%d_,dm23sq%d.root",
                        idcp_bin,
                        imix13_bin,
                        imix23_bin,
                        idmsq_bin
                );
      TFile f(filename,"recreate");
      ntchisq->Write();
      f.Close();

}
