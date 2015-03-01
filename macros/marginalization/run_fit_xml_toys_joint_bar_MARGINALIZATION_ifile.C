//
// ------------------------------------------------------------------------------------------------
// 
// 
// Macro to test the feasibility of marginalization in the T2K-VALOR framework.
// Start from the macro for nubar mode analysis: different oscillation parameters and probabilities are used to
// oscillate neutrinos and anti-neutrinos. Used in general for the RHC-only mode to search for 
// anti-nue appearance and anti-numu disappearance 
//
// Davide Sgalaberna
//
// ------------------------------------------------------------------------------------------------
//





//
// STRATEGY TO FOLLOW FOR A SINGLE FIT
//
// 1) Store toys in binary files 
// 
// - Run the macro N times to produce N binary files, for a total of 10k toys 
// - Binary files differ between each other for the contained toys and/or 
//   for the osc parameter values contained 
// - Each binary file must be tagged consistently with the contained toys, i.e. number
// - Each binary file must be tagged consistently with the osc parameter values contained
// - Need a bash script to do this in parallel on cluster 
// - Marginalization toys must be always performed w/o statistical fluctuations 
// 
// 2) Read 
//
// - Run the macro to read the N binary files, for a total of 10k toys
// - The macro can be run for sub-groups of osc parameter values and compute the 
//   the marginalized likelihood for some vales of the osc parameters, stored in a histogram. 
//   Then merge the histogram results obtained by running the macro for all the sub-groups.
//   This is also the only possible approach for a 2D fit, that otherwise takes too long.
// - Probably 2D Feldman-Cousins is hard to do, since averaging the likelihood over 
//   10k x 200 x 200 toys for a single fit takes too much. 
//   Also splitting by sub-groups needs probably too many jobs and time. 
//   In this case the Minuit minimization should be faster.
//
// 3) Marginalization of the oscillation parameters
//
// - It should be trivial, since they are done in the same way as the systematic parameters
// - Need to provide a particular distribution (e.g. dcp T2K w/ reactors) if not simply 
//   treated as gaussian or flat. Once done, in a single toy the osc param values are fixed
//   to the extracted values. The output is the distribution of the same kind as when only
//   systematic parameters are extracted. 
// - Randomize all the nuisance parameters. If 1D fit the other free parameters which has to be 
//   marginalized (e.g. fit only mixing23 implies marginalization of dmsq32, flat) 
// - Add method MarginalizeParameter() and IsMarginalized to FitParamList
// - Add randomizing methods: RandomizeFlat(), RandomizeGaus() and RandomizePrior() 
//
//
// 5) 1D Feldman-Cousins
//
// - 1D Feldman-Cousins should take about the same time of a Minuit minimization (a bit more???), 
//   but any multiple scan due to local minima (e.g. along dm2) is not needed any longer.
//   e.g. Single fit w/ 1k x 200 toys takes about 10s and we need 10k x 200 toys.
// - Maybe a smart approach is to implement an algorithm that performs the scan of the osc parameter
//   only for the region around the minimum (a kind of Minuit minimization). 
//   Such method could avoid to skip many of the 200 iterations of the osc param grid scan. 
//   Otherwise we can start with a coarser scan and detect the region where most probably the minimum is. 
//   Then the full scan is performed only for a reduced number of points around the minimum.
//   To check the feasibility comparing the full scan with reduced one.
//   Reduced scan in 1D fit, like mixing23, mixing13 or dcp should be safe. 
//   More complicated for dmsq where there could be several minima. 
//   Such reduced scan in 2D fit needs to best carefully tested.
//
// 6) Todo
//
// - Move the code to the framework
//
// - Allow to either fix the oscillation parameters or fit, i.e. make the scan of the grid (DONE, mixing23 or dmsq32)
//
// - Marginalize oscillation parameters using a predefined prior probability distribution 
//   (e.g. dcp T2K w/ reactors) or as flat/gaussian
//
// - Think and implement approach to reduce the number iterations in the scan of the osc param grid
//
// - Use vectors instead of arrays in order to get number of elements when binary files are read 
//
// - Write script to make marginalization toys on cluster as well as perform a single fit 
//
// - Allow the script to compute the marginalized chisq for sub-ranges of osc param.   (DONE)
//   Instead toys are made for all osc param range together                            
//
// - Find the best-fit param value and chisq (DONE)
//


#include <stdio.h>

//#define param1_min 0.           // first value of param1 
//#define param1_max 1.           // last value of param1 
//#define nbins_tot_param1 21     // tot # grid values for param1

#define mixing23Bar_min 0.           // first value of param1 
#define mixing23Bar_max 1.           // last value of param1 
#define nbins_tot_mixing23Bar 21     // tot # grid values for param1

#define dmsq32Bar_min 2.0E-4           // first value of param1 
#define dmsq32Bar_max 2.0E-2           // last value of param1 
#define nbins_tot_dmsq32Bar 41     // tot # grid values for param1






using namespace std;
using namespace valor;
using namespace valor::conv;

double RandomizeGaus(EPhysParamId oscpar,double mean,double sigma){
  //double mean  = Opt::Instance()->GetGlobalBestFit(oscpar); 
  //double sigma = mean * Opt::Instance()->GetGlobalFitError(oscpar);
  double dm    = Rnd::Instance()->RndToyMCosc().Gaus(mean,sigma);
  // return if a physical value was generated
  if(dm >= 0.){
    return dm;
  }
  //or try again
  return RandomizeGaus();
}
double RandomizeFlat(EPhysParamId oscpar,double min,double max){
  return Rnd::Instance()->RndToyMCosc().Uniform(min,max);
}
double RandomizePrior(TH1F *hprior){  // Todo
}




void run_fit_xml_toys_joint_bar_MARGINALIZATION_ifile
(
 const int    nfiles_marg          = 1,     // number of files where marginalization toys are stored
 const int    ntoys_marg           = 10,    // total number of marginalization toys

 const int    idfile               = 0,    // ID of the margin toy file

 const char * param1               = "mixing23Bar",         // 1st param, mixing23, dmsq, mixing13, deltacp,  
 //const char * param1               = "dmsq32Bar",         // 1st param, mixing23, dmsq32, mixing13, deltacp,  

 int  min_bin_sub_param1          = 0,     // min bin for param1 (read only). Start from 0 
 int  max_bin_sub_param1          = 20,  // max bin for param1 (read only). Last nbins_tot_param1 - 1


 
 const char * outfilename            = "framework",
 const char * infilename             = "framework",
 //const char * outfilename            = "prova",
 //const char * infilename             = "prova",
 //const char * outfilename            = "ToyMarginalization_1k_x_200.bin",
 //const char * infilename             = "ToyMarginalization_1k_x_200.bin",
 
 const bool   readtoys               = true,                                        // true->read binary, false-> make toys        

 const char * xmlfile = "fakedata/201411.fakedata/fakedata_201411_0p170e21_A.xml",  // input XML file 
 bool using_xml_file  = true,
   
 const char * config_file      = "analysis_configuration_1.7E20_prefit.xml",   /* analysis configuration file */
 
 int          tweak_systematics       = 1,                                /* generate toy MC with randomized systematics? */
 
 const char * fit_param_list   = "nuosc3fbar_fixnonbar_12fix+nuisance",  
 //const char * fit_param_list   = "nuosc3fbar_fixnonbar_fix5",  
 //const char * fit_param_list   = "nuosc3fbar_fixnonbar_fix5+nuisance", 
 //const char * fit_param_list    = "nuosc3fbar_numudisapp_fixnonbar_fix4",
 //const char * fit_param_list   = "nuosc3fbar_numudisapp_fixnonbar_fix4+nuisance",                                                                                                                         

 const char * dmsqscheme              = "abs23",                             /* abs23==kNuDmsq12andAbs23 abs31==kNuDmsq12andAbs31 fogli==kNuDmsqFogliLisi */
 const char * hierarchy               = "inverted",                            /* normal or inverted */
 const char * mixparamchoice          = "sinsqtheta",                        /* sinsq2theta or sinsqtheta */
 const char * mixparamoctant          = "one",                               /* one==kNuOctantOne two==kNuOctantTwo */
 
)
{

  TString opt(fit_param_list); 
  TString nameparam1(param1);

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
  





  //
  // Randomize all the nuisance parameters and the 
  // other free parameter which has to be marginalized

  //input_mixing12_bar = util::RandomizeMixing12Bar()
  
  double input_dmsq32_bar = RandomizeFlat(kPhypDm23sq,dmsq32Bar_min,dmsq32Bar_max);
  double input_mixing23_bar = RandomizeFlat(kPhypMixing23,mixing23Bar_min,mixing23Bar_max);
  //cout << "input_dmsq32_bar = " << input_dmsq32_bar << endl;
  //cout << "input_mixing23_bar = " << input_mixing23_bar << endl;
  //return;


  





  if(hyp == kHypOsc3F){	
    if( !nameparam1.Contains("mixing12") ){ 
      fparm1->FixParameter(kPhypMixing12);
    }	
    if( !nameparam1.Contains("mixing23") ){ 
      fparm1->FixParameter(kPhypMixing23);
    }
    if( !nameparam1.Contains("mixing13") ){ 
      fparm1->FixParameter(kPhypMixing13);
    }
    if( !nameparam1.Contains("dmsq32") ){ 
      fparm1->FixParameter(kPhypDm23sq);
    }
    if( !nameparam1.Contains("dmsq12") ){ 
      fparm1->FixParameter(kPhypDm12sq);
    }
    if( !nameparam1.Contains("deltacp") ){ 
      fparm1->FixParameter(kPhypCPV);
    }
  }
  else if(hyp == kHypOsc3FBar){	
    if( !nameparam1.Contains("mixing12Bar") ){ 
      fparm1->FixParameter(kPhypMixing12Bar);
    }	
    if( !nameparam1.Contains("mixing23Bar") ){ 
      fparm1->FixParameter(kPhypMixing23Bar);
    }
    if( !nameparam1.Contains("mixing13Bar") ){ 
      fparm1->FixParameter(kPhypMixing13Bar);
    }
    if( !nameparam1.Contains("dmsq32Bar") ){ 
      fparm1->FixParameter(kPhypDm23sqBar);
    }
    if( !nameparam1.Contains("dmsq12Bar") ){ 
      fparm1->FixParameter(kPhypDm12sqBar);
    }
    if( !nameparam1.Contains("deltacpBar") ){ 
      fparm1->FixParameter(kPhypCPVBar);
    }
  }



   



  // generalize in order to initialize arrays/vectors w/o any problem

  //int nbins_tot_param1 = -999;
  const int nbins_tot_param1 = nbins_tot_mixing23Bar;
  //const int nbins_tot_param1 = nbins_tot_dmsq32Bar;

  ///////////////////////



  
  PhysParamId_t parvalue1;
  double param1_min = -999.;
  double param1_max = -999.;
  double step1 = -999.;

    
  if(hyp == kHypOsc3F){
  }
  else if(hyp == kHypOsc3FBar){
 
    if( !fparm1->IsFixed(kPhypMixing23Bar) ){
      parvalue1 = kPhypMixing23Bar;
      param1_min = mixing23Bar_min;
      param1_max = mixing23Bar_max;
      step1 = (param1_max-param1_min)/(nbins_tot_param1-1);
      //nbins_tot_param1 = nbins_tot_mixing23Bar;
      cout << endl;
      cout << "Minimize oscillation parameter " << param1 << endl;
      cout << "Range: " << param1_min << " - " << param1_max << endl;
      cout << "# grid points: " << nbins_tot_param1 << endl; 
      cout << "step: " << step1 << endl;
      cout << endl;
    }
    else if( !fparm1->IsFixed(kPhypDm23sqBar) ){
      parvalue1 = kPhypDm23sqBar;
      param1_min = dmsq32Bar_min;
      param1_max = dmsq32Bar_max;
      step1 = (param1_max-param1_min)/(nbins_tot_param1-1);
      //nbins_tot_param1 = nbins_tot_dmsq32Bar;
      cout << endl;
      cout << "Minimize oscillation parameter " << param1 << endl;
      cout << "Range: " << param1_min << " - " << param1_max << endl;
      cout << "# grid points: " << nbins_tot_param1 << endl; 
      cout << "step: " << step1 << endl;
      cout << endl;
    }
  }

  

 

  // Make sub-range of osc param1
    
  int nbins_param1 = (max_bin_sub_param1-min_bin_sub_param1)+1; // range from 0 to nbins_tot_param1-1
  
  if( nbins_param1 > nbins_tot_param1 ) {
    cout << "nbins_param1 = " << nbins_param1 << endl;
    cout << endl;
    cout << "nbins_tot_param1 = " << nbins_tot_param1 << endl;
    cout << "Partial # of grid bins > Total # of grid bins!!!" << endl;
    exit(1);
  }
    

  // Set histograms of likelihood grid
  int n_param1 = (param1_min - param1_min)/step1 + 1; // +1 for last point
  double min1 = param1_min;
  double max1 = param1_max+step1;
  TH1F *hchisq_marg = new TH1F("hchisq_marg","hchisq_marg",nbins_tot_param1,min1,max1); // the grid point value is the binlowedge 

  double like[nbins_tot_param1];       // convert to vector
  double chisq_marg[nbins_tot_param1]; // convert to vector
  
  TH1F *hToyMarg[ntoys_marg][nbins_tot_param1];

  //double vToyMarg[ntoys_marg][util::NumErecoBins(kSKsRHCFCFV1RingMuLike)][nbins_tot_param1];
  double vToyMarg[ntoys_marg][nbins_tot_param1][util::NumErecoBins(kSKsRHCFCFV1RingMuLike)];


  // Get toys for marginalization (statistical fluctuations must be turned off!!!)

  
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
	

	// Produce toys for each osc param 1 value

	for(int ipar=0; ipar<nbins_tot_param1; ipar++){ 
	  
	  double currval = param1_min + ipar*step1;
	  //cout << "currval = " << currval << endl;
	  
	  fparm1->Update(parvalue1,currval);
	  

	  applied_negevents_condition = skobs->CreateMockData(*skexp, PhysHypothesis(*fparm1), false); // build marginalization distribution w/ syst + osc param but no stat fluctuations
	  
	  for(int ib = 0; ib < util::NumErecoBins(kSKsRHCFCFV1RingMuLike); ib++){ 
	    //vToyMarg[iexperiment][ib][ipar] = skobs->ErecoSpectrumArray(kSKsRHCFCFV1RingMuLike)[ib];
	    vToyMarg[iexperiment][ipar][ib] = skobs->ErecoSpectrumArray(kSKsRHCFCFV1RingMuLike)[ib];
	  }
	  
	} // nbins_tot_param


      } // applied_negevents_condition
    } // iexperiment
    
    
    // Store the toys in a file
    int ok = 0;
    char outname[200];
    sprintf(outname,"%s_file%d_toy%d_x_oscpar%d.bin",infilename,idfile,ntoys_marg,nbins_tot_param1);
    FILE* out = fopen(outname,"wb");  
    if (out != NULL){
      ok = fwrite(vToyMarg, ntoys_marg*util::NumErecoBins(kSKsRHCFCFV1RingMuLike)*nbins_tot_param1 * sizeof (double), 1, out) == 1;
      fclose(out);
    }
  } // ! read toys --> make toys
  









  else{ // read toys

    // Read the toys from a binary file
    
    cout << endl;
    cout << "Read the toys from a binary file" << endl;
    cout << endl;

    int NtoysTot = 0;
    int NparTot = 0;
    
    for(int ipar=0;ipar<nbins_tot_param1;ipar++){
      like[ipar]=0;
      chisq_marg[ipar]=0.;
    }

    for(int ifile=0;ifile<nfiles_marg;ifile++){

      // Read the file
      
      //double vToyMarg[ntoys_marg][util::NumErecoBins(kSKsRHCFCFV1RingMuLike)][nbins_tot_param1];
      double vToyMarg[ntoys_marg][nbins_tot_param1][util::NumErecoBins(kSKsRHCFCFV1RingMuLike)];
      
      int ok = 0;
      char inname[200];
      sprintf(inname,"%s_file%d_toy%d_x_oscpar%d.bin",infilename,ifile,ntoys_marg,nbins_tot_param1);
      FILE* in = fopen(inname,"rb");
      if (in != NULL){
	ok = fread(vToyMarg, ntoys_marg*util::NumErecoBins(kSKsRHCFCFV1RingMuLike)*nbins_tot_param1 * sizeof (double), 1, in) == 1;
	fclose(in);
      }  
      
      // Now all the toys are stored in an array. Compute the marginalized likelihood
      
      cout << endl;
      cout << "Now all the toys of file " << inname << " are stored in an array. Compute the likelihood and average over all the marginalization toys." << endl;
      cout << endl;
      
      for(int ipar=min_bin_sub_param1;ipar<=max_bin_sub_param1;ipar++){
	for(int iexperiment = 0; iexperiment < ntoys_marg; iexperiment++) {
	
	  // compute the likelihood for a single marginalization toy
	  double chisq = 0.;
	  for(int ibin=0;ibin<util::NumErecoBins(kSKsRHCFCFV1RingMuLike);ibin++){
	    //double nexp = vToyMarg[iexperiment][ibin][ipar];
	    double nexp = vToyMarg[iexperiment][ipar][ibin];
	    double nobs = skobs_fakedata->ErecoSpectrumArray(kSKsRHCFCFV1RingMuLike)[ibin];
	    chisq += valor::util::ChisqPoisson(nobs,nexp); // penalty term is a simple constant --> not needed
	    
	  }
	  like[ipar] += exp(-chisq/2.) / ntoys_marg;
	  
	  NtoysTot++;
	} // iexperiment
	NparTot++;
      } // ipar
    } // ifile
    
    // Now the marginalized likelihood is computed using toys from all the binary files
    // Compute the chisq for each free osc parameter value
    
    for(int ipar=min_bin_sub_param1;ipar<=max_bin_sub_param1;ipar++){
      chisq_marg[ipar] = -2*log(like[ipar]);
      hchisq_marg->SetBinContent(ipar+1,chisq_marg[ipar]);
    }
    

    


    // Find the minimum and printout the best-fit osc parameter value
    
    int min_bin = hchisq_marg->GetMinimumBin();
    
    double param1_bf = hchisq_marg->GetBinLowEdge(min_bin);
    double chisq_bf = hchisq_marg->GetBinContent(min_bin);

    cout << endl;
    cout << "param1 best-fit = " << param1_bf << endl;
    cout << "chisq best-fit = " << chisq_bf << endl;
    cout << endl;

    hchisq_marg->Draw();
    
    


    // Write output file
    char foutname[200];
    sprintf(foutname,"fchisq_marg_min%d_max%d.root",min_bin_sub_param1,max_bin_sub_param1);
    TFile *fout = new TFile(foutname,"RECREATE");
    hchisq_marg->Write();
    fout->Close();


  } // else read toys
  
  
  return;
}





