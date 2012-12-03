//
// -----------------------------------------------------------------------------------
// fit the nuisance parameters for the pearticular value of sinsq2theta23 and dm23sq and get the chiSq of the fit
// Used for the nonFC CL contours (syst included) construction 
// fit per bin takes around 1h => run separate jobs for each bin
// ..................................................................
//
// Jelena Ilic
// ---------------------------------------------------------------------------------------    
//

void 2d_nuisance_fit
(
 int          input_sinsq2theta23_bin = 1,                                /* */
 int          input_dm23sq_bin        = 2,                                /* */
 int          min_dm23sq_bin          = 11,                               /* */
 int          max_dm23sq_bin          = 25,                               /* */
 int          min_sinsq2theta23_bin   = 81,                               /* */
 int          max_sinsq2theta23_bin   = 101,                              /* */
 int          dataset_id              = 0,                                /* number in <data_set_%d> XML tag, see inputs/data.xml */
 const char * data_file               = "inputs/data.xml",                /* input SK data */
 const char * config_file             = "analysis_configuration.xml",     /* analysis configuration file */
 const char * fit_param_list          = "numudisapp2f+nuisance",          /* fit parameter list  */
 const char * fit_option              = "minuit",                         /* fitter option */
 const char * tag                     = "CL-fit"                           /* */ 
 )

{
      using namespace std;
      using namespace t2kana;

      const char * config  = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);
      const char * data_file_fullpath = Form("%s/%s", gSystem->Getenv("T2KOSCANA"), data_file);

      // dm23sq, sinsq2theta23 range and step
      double sinsq2theta23_min   = 0.0;
      double sinsq2theta23_max   = 1.0;
      double sinsq2theta23_step  = 0.01;
      double dm23sq_min          = 1E-3;
      double dm23sq_max          = 6E-3;
      double dm23sq_step         = 1E-4;

      // make sure that min and max are at the bin centre and not at the low/high edge of the first/last bin
      int sinsq2theta23_nbins = 1 + (int) ((sinsq2theta23_max - sinsq2theta23_min) / sinsq2theta23_step);
      int dm23sq_nbins        = 1 + (int) ((dm23sq_max        - dm23sq_min       ) / dm23sq_step       );
      cout<< "sinsq2theta23_nbins =" << sinsq2theta23_nbins<<endl;
      cout<< "dm23sq_nbins =" << dm23sq_nbins<<endl;
      // param space - use a histo to select dmSq and sinSq2theta
      TH2D * osc_param_space = new TH2D("osc_param_space", "",
                  sinsq2theta23_nbins, sinsq2theta23_min - 0.5*sinsq2theta23_step, sinsq2theta23_max + 0.5*sinsq2theta23_step,
                  dm23sq_nbins, dm23sq_min - 0.5*dm23sq_step, dm23sq_max + 0.5*dm23sq_step);


      cout << "=========================== "<< endl;
      cout << "- Oscillation parameter space binning" << endl;
      cout << "   - sinsq2theta23 : " << osc_param_space->GetXaxis()->GetNbins() << " bins in (" 
            << osc_param_space->GetXaxis()->GetXmin() << ", "
            << osc_param_space->GetXaxis()->GetXmax() << ")" << endl;
      cout << "   - dm23sq        : " << osc_param_space->GetYaxis()->GetNbins() << " bins in (" 
            << osc_param_space->GetYaxis()->GetXmin() << ", "
            << osc_param_space->GetYaxis()->GetXmax() << ") eV^2/c^4" << endl;
      cout << "- Input bin range " << endl;
      cout << "   - sinsq2theta23 : [" << min_sinsq2theta23_bin << ", " << max_sinsq2theta23_bin << "]" << endl;
      cout << "   - dm23sq        : [" << min_dm23sq_bin        << ", " << max_dm23sq_bin        << "]" << endl;
      cout << "  Fit parameter list: " << fit_param_list << endl;
      cout << "=========================== "<< endl;


      // get superk spectrum predictor & configure
      SuperkExpectation * skexp = new SuperkExpectation;
      skexp->Read(config);

      // get superk spectrum observation
      SuperkObservation * skobs = new SuperkObservation;
 
      //data file - full path
      const char * data = Form("%s/%s", gSystem->Getenv("T2KOSCANA"), data_file);
      skobs->Read(data, dataset_id);

      // ntuple to store results
      TNtupleD * clnt = new TNtupleD("clnt","", 
                  "input_sinsq2theta23_bin:input_dm23sq_bin:input_sinsq2theta23:input_dm23sq:chisq_stat:chisq_sys");

      //Some variables I'll be needed
      double chisq_stat= -99999.0;
      double chisq_sys= -99999.0;
      double input_sinsq2theta23 = -99.0;
      double input_dm23sq= -99.0;

      // get dm23sq and sinsq2theta23 at the centre of the given bin

      if (input_sinsq2theta23_bin < min_sinsq2theta23_bin || input_sinsq2theta23_bin > max_sinsq2theta23_bin){
            cerr<<"input_sinsq2theta23_bin Out Of Range! Resend the job with the correct input!"<<endl;
            gApplication->Terminate();
      }else{
            input_sinsq2theta23 = osc_param_space->GetXaxis()->GetBinCenter(input_sinsq2theta23_bin);
      }
      
      if ( (input_dm23sq_bin < min_dm23sq_bin) || (input_dm23sq_bin > max_dm23sq_bin) ){
            cerr<<"input_dm23sq_bin Out Of Range! Resend the job with the correct input!"<<endl;
            gApplication->Terminate();
      }else{
            input_dm23sq = osc_param_space->GetYaxis()->GetBinCenter(input_dm23sq_bin);
      }
      
      cout<<"input_sinsq2theta23_bin = "<<input_sinsq2theta23_bin<<endl;
      cout<<"input_sinsq2theta23 = "<<input_sinsq2theta23<<endl;
      cout<<"input_dm23sq_bin = "<<input_dm23sq_bin<<endl;
      cout<<"input_dm23sq = "<<input_dm23sq<<endl;
      
      //==============

      // get specified fit parameter list
      FitParamList * fparm = FitParamList::NamedParameterList(fit_param_list);
      assert(fparm);

      conv::PhysHypothesis_t hyp = fparm->PhysicsHypothesis();
      cout << "Physics hypothesis for current job: " << util::AsString(hyp) << endl;

      // Update the oscillation parameter values
      fparm->SetInitValue(t2kana::conv::kPhypSinsq2Theta23, input_sinsq2theta23);
      fparm->SetInitValue(t2kana::conv::kPhypDm23sq,        input_dm23sq);
      
      //and then fix them
      fparm->FixParameter( t2kana::conv::kPhypSinsq2Theta23 );
      fparm->FixParameter( t2kana::conv::kPhypDm23sq );
    

      //====== Calculate the chiSq ===================

      // create the fit inputs
      FitInputs * fit_inputs = new FitInputs(skobs,skexp,fparm);

      //Get the ChiSq
      chisq_stat = fit_inputs->Chisquare();
      cout<<"chisq_stat= " <<chisq_stat<<endl;

      //==== sys ====== 
      // get specified fitter
      FitterABC * fitter = new Superk1RmuErecoFitter(fit_inputs,fit_option);
      assert(fitter);

      //do the fit
      fitter->DoTheFit();

      // get best-chisq
      chisq_sys = fit_inputs->BestFitChisquare();
      cout<<"chisq_sys = "<<chisq_sys<<endl;

      //fill the ntuple
      clnt->Fill(input_sinsq2theta23_bin,input_dm23sq_bin,input_sinsq2theta23,input_dm23sq,chisq_stat,chisq_sys);
      delete fit_inputs;
      delete fparm;

      // write-out results
      TFile f(Form("%s_binS%d_binM%d.root",tag,input_sinsq2theta23_bin,input_dm23sq_bin), "recreate"); 
      clnt->Write("clnt");
      f.Close();
}
