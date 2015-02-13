//
// ------------------------------------------------------------------------------------------------
// Run fits to toy MC experiments generated on-the-fly and store the toy MC, the best-fit parameter
// values and the best-fit spectra, and the true input parameter values.
//
// Costas Andreopoulos
//
// Fit for nubar mode analysis: different oscillation parameters and probabilities are used to
// oscillate neutrinos and anti-neutrinos. Used in general for the RHC-only mode to search for 
// anti-nue appearance and anti-numu disappearance 
//
// Davide Sgalaberna
//
// ------------------------------------------------------------------------------------------------
//

#include <stdio.h>

void run_fit_xml_toys_joint_bar_MARGINALIZATION
(
 const int    ntoys_marg             = 1000,
 const int    nbins_tot_param1       = 200,

 const char * outfilename            = "ToyMarginalization_1k_x_200.bin",
 const char * infilename             = "ToyMarginalization_1k_x_200.bin",
 const bool readtoys                 = false,  // true->read binary, false-> make toys        

 const char * xmlfile = "fakedata/201411.fakedata/fakedata_201411_0p170e21_F.xml",       /* input XML file Run5-6 fakedata */ 
 bool using_xml_file = true,
 

 const char * tag                     = "toyexptfits",    /* tag for output file */
 int          run                     = 0,                /* run number */
 int          first_experiment        = 0,                       /* */
 int          last_experiment         = 0,                       /* */


  double       input_mixing12          = -99, //0.306,                                /* 12-mixing, sinq2theta or sinsqtheta depending on user choice, if -999 then randomize */
  double       input_mixing13          = -99, //0.02146, //0.0251,                                /* 13-mixing, sinq2theta or sinsqtheta depending on user choice, if -999 then randomize */
  double       input_mixing23          = -99, //0.50,                                /* 23-mixing, sinq2theta or sinsqtheta depending on user choice */
  double       input_deltacp           = -99, //0.,                                /* deltacp, if -999 then randomize x TMath::Pi()!!!!*/
  double       input_dm12sq            = -99, //7.58E-5,                             /* Dm^2_{12}, if -999 then randomize */
  double       input_dm23sq            = -99, //2.32E-3,                             /* |Dm^2_{23}| */

  double       input_mixing12_bar      = -99, //0.306,                                /* 12-mixing, sinq2theta or sinsqtheta depending on user choice, if -999 then randomize */
  double       input_mixing13_bar      = -99, //0.02146, //0.0251,                    /* 13-mixing, sinq2theta or sinsqtheta depending on user choice, if -999 then randomize */
 double       input_mixing23_bar       = -99, //-99, //0.7, //0.50,                                /* 23-mixing, sinq2theta or sinsqtheta depending on user choice */
  double       input_deltacp_bar       = -99, //0.,                                /* deltacp, if -999 then randomize x TMath::Pi()!!!!*/
  double       input_dm12sq_bar        = -99, //7.58E-5,                             /* Dm^2_{12}, if -999 then randomize */
 double       input_dm23sq_bar         = -99, //-99, //2.7E-3, //2.32E-3,                             /* |Dm^2_{23}| */


 const char * config_file      = "analysis_configuration_1.7E20_prefit.xml",   /* analysis configuration file */

  bool         stat_fluctuations       = false,                                /* include statistical fluctuation in toy MC expts? */
  int          tweak_systematics       = 1,                                /* generate toy MC with randomized systematics? */



 const char * fit_param_list    = "nuosc3fbar_numudisapp_fixnonbar_fix4+nuisance",

  bool         store_toy_expts         = true,                                /* disable this if you don't want to store the toy MC expts     */
  bool         store_bestfit_spectrum  = false,                               /* disable this if you don't want to store the best-fit spectra */
  const char * dmsqscheme              = "abs23",                             /* abs23==kNuDmsq12andAbs23 abs31==kNuDmsq12andAbs31 fogli==kNuDmsqFogliLisi */
  const char * hierarchy               = "inverted",                            /* normal or inverted */
  const char * mixparamchoice          = "sinsqtheta",                        /* sinsq2theta or sinsqtheta */
  const char * mixparamoctant          = "one",                               /* one==kNuOctantOne two==kNuOctantTwo */

 const char *optfixpar         = "mixing12, mixing23, dmsq, dmsq12, mixing13, deltacp"  // list of the params to be fixed in the fit   
)
{
  using namespace std;
  using namespace valor;
  using namespace valor::conv;

  TString opt(fit_param_list); 
  TString optfixparams(optfixpar);

  double bestfit_chisq,bestfit_chisq2;

  // Analysis configuration
  TString config = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);
  TString xmlfile2 = Form("%s/inputs/%s", gSystem->Getenv("T2KOSCANA"), xmlfile);

  FileStat_t dummy;

  // Check for config file, and XML data file.
  // If not found in inputs or config, see if it works as a normal path
  if (gSystem->GetPathInfo(config.Data(), dummy)) 
  {
    config = config_file;
    if (gSystem->GetPathInfo(config.Data(), dummy)) {
      cerr<<"Config file not found: "<<config_file<<endl;
      exit(-1);
    }
  }
  
  if (gSystem->GetPathInfo(xmlfile2.Data(), dummy)) 
  {
    xmlfile2 = xmlfile;
    if (gSystem->GetPathInfo(xmlfile2.Data(), dummy)) {
      cerr<<"Data input file not found: "<<xmlfile<<endl;
      exit(-1);
    }
  }
  

  int ifakedata = 0;
  SuperkObservation * skobs_fakedata = new SuperkObservation;
  skobs_fakedata->Read(xmlfile2.Data(),ifakedata,using_xml_file);
  //skobs_fakedata->ErecoSpectrum(kSKsRHCFCFV1RingMuLike)->Draw();
  //return;

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


  // Get superk spectrum predictor & configure
  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read(config.Data());

  // Fit parameter list
  FitParamList * fparm1  = util::NamedParameterList(fit_param_list); //only used to initialise branch
  PhysHypothesis_t hyp = fparm1->PhysicsHypothesis();
  cout << "- The following physics hypothesis is used in the fit: " << util::AsString(hyp) << endl;
 
 
  ////////////////////////
  //                    //
  //   Fix parameters   //
  //                    //
  ////////////////////////
  
  if(hyp == kHypOsc3F){	
    if(optfixparams.Contains("mixing23")){ 
      fparm1->FixParameter(kPhypMixing23);
    }
    if(optfixparams.Contains("mixing13")){ 
      fparm1->FixParameter(kPhypMixing13);
    }
    if(optfixparams.Contains("dmsq32")){ 
      fparm1->FixParameter(kPhypDm23sq);
    }
    if(optfixparams.Contains("dmsq12")){ 
      fparm1->FixParameter(kPhypDm12sq);
    }
    if(optfixparams.Contains("deltacp")){ 
      fparm1->FixParameter(kPhypCPV);
    }
  }
  else if(hyp == kHypOsc3FBar){	
    if(optfixparams.Contains("mixing12")){ 
      fparm1->FixParameter(kPhypMixing12Bar);
    }	
    if(optfixparams.Contains("mixing23")){ 
      fparm1->FixParameter(kPhypMixing23Bar);
    }
    if(optfixparams.Contains("mixing13")){ 
      fparm1->FixParameter(kPhypMixing13Bar);
    }
    if(optfixparams.Contains("dmsq32")){ 
      fparm1->FixParameter(kPhypDm23sqBar);
    }
    if(optfixparams.Contains("dmsq12")){ 
      fparm1->FixParameter(kPhypDm12sqBar);
    }
    if(optfixparams.Contains("deltacp")){ 
      fparm1->FixParameter(kPhypCPVBar);
    }
  }
   

  // Get toys for marginalization (statistical fluctuations must be turned off!!!)

  TH1F *hToyMarg[ntoys_marg][nbins_tot_param1];
  double vToyMarg[ntoys_marg][util::NumErecoBins(kSKsRHCFCFV1RingMuLike)][nbins_tot_param1];
  
  if(!readtoys){ // make toys
  
    cout << endl;
    cout << "Make the toys and store in a binary file" << endl;
    cout << endl;
  
    for(int iexperiment = 0; iexperiment < ntoys_marg; iexperiment++) {
    
      cout << endl;
      cout << "iexperiment " << iexperiment << endl;
      cout << "# params " << nbins_tot_param1 << endl;
      
      fparm1->ResetCurrToInitValues();
      
      SuperkObservation * skobs = new SuperkObservation;
      
      bool applied_negevents_condition = true;
      while(applied_negevents_condition) {
	if(tweak_systematics) {
	  
	  // tweak all sources of systematic error
	  if(tweak_systematics == 1) {
	    cout << "- Randomizing all systematic params in pdf..." << endl;
	    skexp->TweakAllSystematics(); 
	  }
	  else if(tweak_systematics == 2) {
	    cout << "- Randomizing BANFF systematic params in pdf..." << endl;
	    skexp->TweakAllBANFFSystematics();
	  }
	  else if(tweak_systematics == 3) {
	    cout << "- Randomizing uncorrelated-xsec systematic params in pdf..." << endl;
	    skexp->TweakAllUncorrelatedXSecSystematics();
	  }
	  else if(tweak_systematics == 4) {
	    cout << "- Randomizing SKDetFSI systematic params in pdf..." << endl;
	    skexp->TweakAllSKDetectorFSISystematics();
	  }
	  //skexp->PrintSystematicVariations();
	  //cout << "- Reseting all systematic params in pdf..." << endl;
	  //skexp->ResetSystematicVariations();
	}
	else { // no toys
	  skexp->ResetSystematicVariations();
	  cout << "- Generating toy-MC without systematic variation" << endl;
	}
	
	// loop over oscillation parameter points to scan and create a different distribution
	
	double param1_min = 0.;
	double param1_max = 1.;
	double param1_step;
	if(nbins_tot_param1==1) param1_step = 0.;
	else param1_step = (param1_max-param1_min)/(nbins_tot_param1-1);
	
	for(int ipar=0; ipar<nbins_tot_param1; ipar++){
	  
	  double currval = param1_min + ipar*param1_step;
	  //cout << "currval = " << currval << endl;
	  
	  EPhysParamId parvalue1 = conv::kPhypMixing23Bar;
	  //fparm1->SetInitValue(parvalue1,currval);
	  fparm1->Update(parvalue1,currval);
	  //fparm1->FixParameter(parvalue1);
	  //cout << fparm1->CurrValue(parvalue1) << endl;
	  
	  applied_negevents_condition = skobs->CreateMockData(*skexp, PhysHypothesis(*fparm1), false); // build marginalization distribution w/ syst + osc param
	  
	  for(int ib = 0; ib < util::NumErecoBins(kSKsRHCFCFV1RingMuLike); ib++){ 
	    vToyMarg[iexperiment][ib][ipar] = skobs->ErecoSpectrumArray(kSKsRHCFCFV1RingMuLike)[ib];
	  }
	  
	} // nbins_tot_param
      } // applied_negevents_condition
    } // iexperiment
    
    
    // Store the toys in a file
    int ok = 0;
    FILE* out = fopen(outfilename, "wb");  
    if (out != NULL)
      {
	ok = fwrite(vToyMarg, ntoys_marg*util::NumErecoBins(kSKsRHCFCFV1RingMuLike)*nbins_tot_param1 * sizeof (double), 1, out) == 1;
	fclose(out);
      }
 
  } // ! read toys --> make toys
  


  else{ // read toys

    // Read the toys from a binary file
    
    cout << endl;
    cout << "Read the toys from a binary file" << endl;
    cout << endl;

    double vToyMarg[ntoys_marg][util::NumErecoBins(kSKsRHCFCFV1RingMuLike)][nbins_tot_param1];
    int ok = 0;
    FILE* in = fopen(infilename,"rb");
    if (in != NULL)
      {
	ok = fread(vToyMarg, ntoys_marg*util::NumErecoBins(kSKsRHCFCFV1RingMuLike)*nbins_tot_param1 * sizeof (double), 1, in) == 1;
	fclose(in);
      }  
   
    //for(int ipar=0;ipar<nbins_tot_param1;ipar++){
    //for(int iexperiment = 0; iexperiment < ntoys_marg; iexperiment++) {
    //for(int ibin=0;ibin<util::NumErecoBins(kSKsRHCFCFV1RingMuLike);ibin++){
    //cout << vToyMarg[iexperiment][ibin][ipar] << " " << vToyMarg_in[iexperiment][ibin][ipar] << endl;
    //}
    //}
    //}
  }

 
  // Now all the toys are stored in an array. Compute the marginalized likelihood
  
  cout << endl;
  cout << "Now all the toys are stored in an array. Compute the likelihood and average over all the marginalization toys." << endl;
  cout << endl;

  TH1F *hchisq_marg = new TH1F("hchisq_marg","hchisq_marg",nbins_tot_param1,0,nbins_tot_param1);

  double like[nbins_tot_param1];
  double chisq_marg[nbins_tot_param1];

  for(int ipar=0;ipar<nbins_tot_param1;ipar++){
    
    like[ipar]=0;
    chisq_marg[ipar]=0.;

    for(int iexperiment = 0; iexperiment < ntoys_marg; iexperiment++) {

      // compute the likelihood for a single marginalization toy
      double chisq = 0.;
      for(int ibin=0;ibin<util::NumErecoBins(kSKsRHCFCFV1RingMuLike);ibin++){
	double nexp = vToyMarg[iexperiment][ibin][ipar];
	double nobs = skobs_fakedata->ErecoSpectrumArray(kSKsRHCFCFV1RingMuLike)[ibin];
	chisq += valor::util::ChisqPoisson(nobs,nexp); // penalty term is a simple constant --> not needed
      }
      
      like[ipar] += exp(-chisq/2.) / ntoys_marg;
      
      //cout << chisq << ", ";
      //cout << like[ipar] << " ";
      //cout << -TMath::Log(like[ipar]) << " ";
    }
    //cout << endl;
    
    chisq_marg[ipar] = -2*TMath::Log(like[ipar]);
    
    hchisq_marg->SetBinContent(ipar+1,chisq_marg[ipar]);
  }

  
  hchisq_marg->Draw();


    
  return;
}
