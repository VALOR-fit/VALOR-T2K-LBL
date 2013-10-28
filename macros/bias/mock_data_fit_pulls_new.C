//
// ------------------------------------------------------------------------------------------------
// Read macros/fit/run_fit.C ROOT outputs with input and best-fit values for fits to toy-MC expts.
// Plot various residual and pull distributions for fit parameters.
// ------------------------------------------------------------------------------------------------
//

//#define FITPARMLIST 0 // nuosc3f_numudisapp_fix4+nuisance
//#define FITPARMLIST 1 // nuosc3f_numudisapp_nuis3fixdcp+nuisance
//#define FITPARMLIST 2 // nuosc3f_numudisapp_nuis4+nuisance
//#define FITPARMLIST -3 // only want to look at sk escale minos errors
#define FITPARMLIST -2 // mix23 & dm32sq only

//#define VERBOSE //print out the best fit and input of every single parameter in every single event
//#define VERYVERBOSE //print out parameter name too

void save_canvas(TCanvas * c, const char * name, const char * tag) {
  string name2(name); //don't ask me why I need to do this. I have no idea...
  c->SaveAs(Form("%s_%s.pdf",tag,name));
  c->SaveAs(Form("%s_%s.eps",tag,name));
  c->SaveAs(Form("%s_%s.C",  tag,name2.c_str()));
}

void mock_data_fit_pulls_new
(
  const char * filenames,                                              /* add these filenames to the tree */
  const char * tag,                                                    /* tag for the output files */
  int          err_mode             = 1,                               /* which errors to use for pulls (1: migrad, 2 hesse, 3 minos) */
  int          fitter               = 0,                               /* which fitter has been used (0: minuit, 1: minuit2) */
  const char * config_file          = "analysis_configuration.xml",    /* file to get the pre errors from */
  bool         do_gaus_fits         = false,
  bool         is_sinsqtheta        = true                             /* global best fit comes from sinsqtheta or sinsq2theta */
)
{
  using namespace std;
  using namespace std::string;
  using namespace t2kana;
  using namespace t2kana::conv;

  if(is_sinsqtheta)
    Opt::Instance()->gOscMixingParamChoice = conv::kNuMixSinsqTheta;
  else
    Opt::Instance()->gOscMixingParamChoice = conv::kNuMixSinsq2Theta;

  // First get the expected systematic variations
  if(do_gaus_fits)
    gStyle->SetOptFit(1);         // show fit results
  gStyle->SetOptStat(1100);     // show mean & RMS only
  //gStyle->SetOptStat(1101);     // show mean, RMS, name
  //gStyle->SetOptStat(1111);     // show mean, RMS, number of entries, name
  //gStyle->SetStatStyle(0);      // make stat boxes see through
  gStyle->SetStatBorderSize(0); //
  gStyle->SetTextFont(42);      // less bold
  gStyle->SetFitFormat("20.6g");//controls format of the fit parameters in the stat box

  //
  // Parameters in the fit
  //
#if FITPARMLIST == 0
  const int nparams = 43;
#elif FITPARMLIST == 1
  const int nparams = 46;
#elif FITPARMLIST == 2
  const int nparams = 47;
#elif FITPARMLIST == -3
  const int nparams = 1;
#elif FITPARMLIST == -2
  const int nparams = 2;
#endif

  const int nosc     = 2; // nu of osc params in list below (osc params are plotted on separate canvas)
  const int nsystcol = 3; // number of columns on canvases showing residual/pulls for systematics
  const int nsystrow = 4; // number of rows    on canvases showing residual/pulls for systematics

  PhysParamId_t params[nparams] = {
#if FITPARMLIST != -3
     kPhypMixing23,
     kPhypDm23sq,  
#endif
#if FITPARMLIST >= 1
     kPhypMixing13,  
     kPhypMixing12,   
     kPhypDm12sq,
#if FITPARMLIST == 2
     kPhypDeltaCP,
#endif //2
#endif //>=1
#if FITPARMLIST >= 0
     kNuisBANFF00,
     kNuisBANFF01,
     kNuisBANFF02,
     kNuisBANFF03,
     kNuisBANFF04,
     kNuisBANFF05,
     kNuisBANFF06,
     kNuisBANFF07,
     kNuisBANFF08,
     kNuisBANFF09,
     kNuisBANFF10,
     kNuisBANFF11,
     kNuisBANFF12,
     kNuisBANFF13,
     kNuisBANFF14,
     kNuisBANFF15,
     kNuisBANFFMAQE,         
     kNuisBANFFMARES,        
     kNuisBANFFCCQEE0,       
     kNuisBANFFCCQEE1,        
     kNuisBANFFCCQEE2,        
     kNuisBANFFCC1piE0,       
     kNuisBANFFCC1piE1,       
   //kNuisBANFFMECNorm,        
   //kNuisBANFFPDDDoubleCount, 
     kNuisXSecFermiMom,         
     kNuisXSecSpecFunc,         
     kNuisXSecWshape,           
     kNuisXSecShapeCCoth,       
     kNuisXSecPiLessDelta,      
     kNuisXSecBindE,            
     kNuisXSecNormCCcoh,         
     kNuisXSecNormNC1piC,        
     //kNuisXSecNormNCcoh,         
     kNuisXSecNormNCoth,         
     kNuisXSecNormNuEToNuMu,     
     kNuisXSecNormAntiNuToNu,    
#endif //>= 0
#if FITPARMLIST != -2
     kNuisSuperkEnergyScale,
#endif
#if FITPARMLIST >= 0
     kNuisSKDetFSI000,
     kNuisSKDetFSI001,
     kNuisSKDetFSI002,
     kNuisSKDetFSI003,
     kNuisSKDetFSI004,
     kNuisSKDetFSI005
#endif
  };

  //
  // Get 1-sigma (prefit) errors & nominal values for all fit params
  //

  // Parse the input XML configuration file
  const char * config  = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);
  TXMLEngine* xml = new TXMLEngine;   
  XMLDocPointer_t xmldoc = xml->ParseFile(config);
  if (xmldoc==0) {
    cerr << "Exp >> " << "Can not parse input: " << config << endl;
    delete xml;
    exit(1);
  }
  XMLNodePointer_t root_node = xml->DocGetRootElement(xmldoc);

  TMatrixD * banffcov_ptr = util::XML2TMatrixD(xml,root_node, "banff/covariance_matrix");
  const TMatrixD & banffcov = *banffcov_ptr;
  TVectorD * banffbf_ptr = util::XML2TVectorD(xml,root_node, "banff/best_fit_values");
  const TVectorD & banffbf = *banffbf_ptr;
  TMatrixD * skfsicov_ptr = util::XML2TMatrixD(xml,root_node, "skdetfsi/covariance_matrix");
  const TMatrixD & skfsicov = *skfsicov_ptr;

  double escale_err = 1;//util::XML2Dbl(xml,root_node,"detector_systematics/superk/energy_scale_relative_error/one_mu_ring");

  double mix23_err   = 0;
  double mix12_err   = 0;
  double mix13_err   = 0;
  double dm32sq_err  = kGlobFitDm23sq        * kGlobFitErrDm23sq;
  double dm21sq_err  = kGlobFitDm12sq        * kGlobFitErrDm12sq;
  double deltacp_err = conv::kpi             * conv::kpi;
  if(is_sinsqtheta) {
    mix23_err   = kGlobFitSinsqTheta23 * kGlobFitErrSinsqTheta23;
    mix12_err   = kGlobFitSinsqTheta12 * kGlobFitErrSinsqTheta12;
    mix13_err   = kGlobFitSinsqTheta13 * kGlobFitErrSinsqTheta13;
    
  }
  else {
    mix23_err   = kGlobFitSinsq2Theta23 * kGlobFitErrSinsq2Theta23;
    mix12_err   = kGlobFitSinsq2Theta12 * kGlobFitErrSinsq2Theta12;
    mix13_err   = kGlobFitSinsq2Theta13 * kGlobFitErrSinsq2Theta13;
  }

  double prefit_err_array[nparams] = 
  {
#if FITPARMLIST != -3
     mix23_err,                    // kPhypMixing23
     dm32sq_err,                   // kPhypDm23sq
#endif
#if FITPARMLIST >= 1
     mix13_err,                    // kPhypMixing13
     mix12_err,                    // kPhypMixing12
     dm21sq_err,                   // kPhypDm12sq
#if FITPARAMLIST == 2
     deltacp_err,                  // kPhypDeltaCP
#endif //==2
#endif //>=1
#if FITPARMLIST >= 0
     TMath::Sqrt(banffcov( 0, 0)), // kNuisBANFF00
     TMath::Sqrt(banffcov( 1, 1)), // kNuisBANFF01
     TMath::Sqrt(banffcov( 2, 2)), // kNuisBANFF02
     TMath::Sqrt(banffcov( 3, 3)), // kNuisBANFF03
     TMath::Sqrt(banffcov( 4, 4)), // kNuisBANFF04
     TMath::Sqrt(banffcov( 5, 5)), // kNuisBANFF05
     TMath::Sqrt(banffcov( 6, 6)), // kNuisBANFF06
     TMath::Sqrt(banffcov( 7, 7)), // kNuisBANFF07
     TMath::Sqrt(banffcov( 8, 8)), // kNuisBANFF08
     TMath::Sqrt(banffcov( 9, 9)), // kNuisBANFF09
     TMath::Sqrt(banffcov(10,10)), // kNuisBANFF10
     TMath::Sqrt(banffcov(11,11)), // kNuisBANFF11
     TMath::Sqrt(banffcov(12,12)), // kNuisBANFF12
     TMath::Sqrt(banffcov(13,13)), // kNuisBANFF13
     TMath::Sqrt(banffcov(14,14)), // kNuisBANFF14
     TMath::Sqrt(banffcov(15,15)), // kNuisBANFF15
     TMath::Sqrt(banffcov(16,16)), // kNuisBANFFMAQE
     TMath::Sqrt(banffcov(17,17)), // kNuisBANFFMARES        
     TMath::Sqrt(banffcov(18,18)), // kNuisBANFFCCQEE0       
     TMath::Sqrt(banffcov(19,19)), // kNuisBANFFCCQEE1        
     TMath::Sqrt(banffcov(20,20)), // kNuisBANFFCCQEE2        
     TMath::Sqrt(banffcov(21,21)), // kNuisBANFFCC1piE0       
     TMath::Sqrt(banffcov(22,22)), // kNuisBANFFCC1piE1       
                                   // kNuisBANFFMECNorm        
                                   // kNuisBANFFPDDDoubleCount 
     1,                            // kNuisXSecFermiMom         
     1,                            // kNuisXSecSpecFunc         
     1,                            // kNuisXSecWshape           
     1,                            // kNuisXSecShapeCCoth       
     1,                            // kNuisXSecPiLessDelta      
     1,                            // kNuisXSecBindE            
     1,                            // kNuisXSecNormCCcoh         
     1,                            // kNuisXSecNormNC1piC        
     //1,                            // kNuisXSecNormNCcoh         
     1,                            // kNuisXSecNormNCoth         
     1,                            // kNuisXSecNormNuEToNuMu     
     1,                            // kNuisXSecNormAntiNuToNu
#endif //>=0
#if FITPARMLIST != -2
     escale_err,                   // kNuisSuperkEnergyScale
#endif
#if FITPARMLIST >= 0
     TMath::Sqrt(skfsicov(0,0)),   // kNuisSKDetFSI000
     TMath::Sqrt(skfsicov(1,1)),   // kNuisSKDetFSI001
     TMath::Sqrt(skfsicov(2,2)),   // kNuisSKDetFSI002
     TMath::Sqrt(skfsicov(3,3)),   // kNuisSKDetFSI003
     TMath::Sqrt(skfsicov(4,4)),   // kNuisSKDetFSI004
     TMath::Sqrt(skfsicov(5,5))    // kNuisSKDetFSI005
#endif
  };


  double global_bestfit_mix23 = 0;
  double global_bestfit_mix13 = 0;
  double global_bestfit_mix12 = 0;
  if(is_sinsqtheta) {
    global_bestfit_mix23 = kGlobFitSinsqTheta23;
    global_bestfit_mix13 = kGlobFitSinsqTheta13;
    global_bestfit_mix12 = kGlobFitSinsqTheta12;
  }
  else {
    global_bestfit_mix23 = kGlobFitSinsq2Theta23;
    global_bestfit_mix13 = kGlobFitSinsq2Theta13;
    global_bestfit_mix12 = kGlobFitSinsq2Theta12;
  }

  double nominal_array[nparams] = 
  {
#if FITPARMLIST != -3
     global_bestfit_mix23,     // kPhypMixing23
     kGlobFitDm23sq,           // kPhypDm23sq
#endif
#if FITPARMLIST >= 1
     global_bestfit_mix13,     // kPhypMixing13
     global_bestfit_mix12,     // kPhypMixing12
     kGlobFitDm12sq,           // kPhypDm12sq
#if FITPARMLIST == 2
     kGlobFitDeltaCP,          // kPhypDeltaCP
#endif
#endif
#if FITPARMLIST >= 0
     TMath::Sqrt(banffbf[ 0]), // kNuisBANFF00
     TMath::Sqrt(banffbf[ 1]), // kNuisBANFF01
     TMath::Sqrt(banffbf[ 2]), // kNuisBANFF02
     TMath::Sqrt(banffbf[ 3]), // kNuisBANFF03
     TMath::Sqrt(banffbf[ 4]), // kNuisBANFF04
     TMath::Sqrt(banffbf[ 5]), // kNuisBANFF05
     TMath::Sqrt(banffbf[ 6]), // kNuisBANFF06
     TMath::Sqrt(banffbf[ 7]), // kNuisBANFF07
     TMath::Sqrt(banffbf[ 8]), // kNuisBANFF08
     TMath::Sqrt(banffbf[ 9]), // kNuisBANFF09
     TMath::Sqrt(banffbf[10]), // kNuisBANFF10
     TMath::Sqrt(banffbf[11]), // kNuisBANFF11
     TMath::Sqrt(banffbf[12]), // kNuisBANFF12
     TMath::Sqrt(banffbf[13]), // kNuisBANFF13
     TMath::Sqrt(banffbf[14]), // kNuisBANFF14
     TMath::Sqrt(banffbf[15]), // kNuisBANFF15
     TMath::Sqrt(banffbf[16]), // kNuisBANFFMAQE
     TMath::Sqrt(banffbf[17]), // kNuisBANFFMARES        
     TMath::Sqrt(banffbf[18]), // kNuisBANFFCCQEE0       
     TMath::Sqrt(banffbf[19]), // kNuisBANFFCCQEE1        
     TMath::Sqrt(banffbf[20]), // kNuisBANFFCCQEE2        
     TMath::Sqrt(banffbf[21]), // kNuisBANFFCC1piE0       
     TMath::Sqrt(banffbf[22]), // kNuisBANFFCC1piE1       
                               // kNuisBANFFMECNorm        
                               // kNuisBANFFPDDDoubleCount 
     0,                        // kNuisXSecFermiMom         
     0,                        // kNuisXSecSpecFunc         
     0,                        // kNuisXSecWshape           
     0,                        // kNuisXSecShapeCCoth       
     0,                        // kNuisXSecPiLessDelta      
     0,                        // kNuisXSecBindE            
     0,                        // kNuisXSecNormCCcoh         
     0,                        // kNuisXSecNormNC1piC        
     //0,                        // kNuisXSecNormNCcoh         
     0,                        // kNuisXSecNormNCoth         
     0,                        // kNuisXSecNormNuEToNuMu     
     0,                        // kNuisXSecNormAntiNuToNu
#endif //>= 0
#if FITPARMLIST != -2
     0,                        // kNuisSuperkEnergyScale
#endif
#if FITPARMLIST >= 0
     1,                        // kNuisSKDetFSI000
     1,                        // kNuisSKDetFSI001
     1,                        // kNuisSKDetFSI002
     1,                        // kNuisSKDetFSI003
     1,                        // kNuisSKDetFSI004
     1                         // kNuisSKDetFSI005
#endif
  };

  // Release memory
  xml->FreeDoc(xmldoc);
  delete xml;


  //
  // Book histograms
  //

  const int    nbins_pull = 100;
  const double min_pull   = -5;
  const double max_pull   =  5;
  const int    nbins_err  = 100;

  TH1I * h_migradminstat  = new TH1I("h_migradminstat", "migrad minimization status",  11, -1,   10);
  TH1I * h_migraderrstat  = new TH1I("h_migraderrstat", "migrad error status",         14, -1,   13);
  TH1I * h_hesseerrstat   = new TH1I("h_hesseerrstat",  "hesse error status",          14, -1,   13);
  TH1D * h_fittime        = new TH1D("h_fittime",       "fit time",                  1000,  0, 1000);
  TH1D * h_chisqbf        = new TH1D("h_chisqbf",       "best-fit chisq",             615, -5,  200);
  TH1D * h_nsk_nom        = new TH1D("h_nsk_nom",       "nsk(1Rmu) nom",              500,  0,  500);
  TH1D * h_nsk_bf         = new TH1D("h_nsk_bf",        "nsk(1Rmu) best-fit",         500,  0,  500);
  // -3 is for mixing23, dm32sq, spectral function
  TH1D * h_summary_nuis   = new TH1D("h_summary_nuis",  "summary of nuisance parameters", nparams-3, 0, nparams-3);

  TH1D * h_bestfit    [nparams] = { 0 }; // best-fit param distributions
  TH1D * h_residuals  [nparams] = { 0 }; // residuals
  TH1D * h_pull_1     [nparams] = { 0 }; // pulls - 1st definition (as in T2K-TN-141v6, p.320)
  TH1D * h_pull_2     [nparams] = { 0 }; //  ...    2nd     ...          ...
  TH1D * h_pull_3     [nparams] = { 0 }; //  ...    3rd     ...          ...
  TH1D * h_err_migrad [nparams] = { 0 }; // migrad err distribution
  TH1D * h_err_hesse  [nparams] = { 0 }; // hesse  err distribution
  TH1D * h_err_minos  [nparams] = { 0 }; // minos  err distribution
  TH1D * h_err_minos_plus   [nparams] = { 0 }; // minos  plus err distribution
  TH1D * h_err_minos_minus  [nparams] = { 0 }; // minos  minus err distribution
  TH2D * h_corr_inputdcp_residual[nparams] = {0}; //residual correlation with input dcp
  TH2D * h_corr_inputdcp_bestfit [nparams] = {0}; //best fit correlation with input dcp
  TH2D * h_corr_fitmix23_residual[nparams] = {0}; //residual correlation with fitted mix23
  TH2D * h_corr_fitmix23_bestfit [nparams] = {0}; //best fit correlation with fitted mix23
  for(int ip = 0; ip < nparams; ip++) 
  {
     PhysParamId_t p = params[ip];

     string name  = util::AsStringCompact(p); 

     // The name might contain unsuitable chars for ROOT obj name, such as `|'.
     // Fiter those chars.
     string name2 = name;
     string filtlist = "|";
     string::size_type pos = 0;
     while( (pos = name2.find_first_of(filtlist, pos)) != string::npos )
     {
       name2.erase(pos,1);
     }

     h_bestfit[ip] = new TH1D( 
         Form("h_bestfit_%s",name2.c_str()), 
         Form("Best-fit distribution for %s", util::AsStringCompact(p).c_str()), 
         nbins_pull, min_pull, max_pull);
     double residual_lo = -5*prefit_err_array[ip];
     double residual_hi =  5*prefit_err_array[ip];
     if(p == kPhypMixing23) {
       residual_lo = -0.5;
       residual_hi =  0.5;
     }
     else if(p == kPhypDeltaCP) {
       residual_lo = -10;
       residual_hi =  10;
     }
     h_residuals[ip] = new TH1D( 
         Form("h_residuals_%s",name2.c_str()), 
         Form("Residuals for %s", util::AsStringCompact(p).c_str()), 
         nbins_pull, residual_lo, residual_hi);
     h_corr_inputdcp_residual[ip] = new TH2D( 
         Form("h_corr_inputdcp_residual_%s",name2.c_str()), 
         Form("Correlation with input #delta_{CP} for %s residual", util::AsStringCompact(p).c_str()),
	 1000,-5,+5,
         nbins_pull, residual_lo, residual_hi);
     h_corr_inputdcp_bestfit[ip] = new TH2D( 
         Form("h_corr_inputdcp_bestfit_%s",name2.c_str()), 
         Form("Correlation with input #delta_{CP} for %s best fit", util::AsStringCompact(p).c_str()),
	 1000,-5,+5,
         nbins_pull, nominal_array[ip] - (2* prefit_err_array[ip]), nominal_array[ip] + (2 * prefit_err_array[ip]));
     h_corr_fitmix23_residual[ip] = new TH2D( 
         Form("h_corr_fitmix23_residual_%s",name2.c_str()), 
         Form("Correlation with fitted sin^2#theta_{23} for %s residual", util::AsStringCompact(p).c_str()),
	 1000,0,1,
         nbins_pull, residual_lo, residual_hi);
     h_corr_fitmix23_bestfit[ip] = new TH2D( 
         Form("h_corr_fitmix23_bestfit_%s",name2.c_str()), 
         Form("Correlation with fitted sin^2#theta_{23} for %s best fit", util::AsStringCompact(p).c_str()),
	 1000,0,1,
         nbins_pull, nominal_array[ip] - (2* prefit_err_array[ip]), nominal_array[ip] + (2 * prefit_err_array[ip]));
     h_pull_1[ip] = new TH1D( 
         Form("h_pull_1_%s",name2.c_str()), 
         Form("Pull_1 for %s", util::AsStringCompact(p).c_str()), 
         nbins_pull, min_pull, max_pull);
     h_pull_2[ip] = new TH1D( 
         Form("h_pull_2_%s",name2.c_str()), 
         Form("Pull_2 for %s", util::AsStringCompact(p).c_str()), 
         nbins_pull, min_pull, max_pull);
     h_pull_3[ip] = new TH1D( 
         Form("h_pull_3_%s",name2.c_str()), 
         Form("Pull_3 for %s", util::AsStringCompact(p).c_str()), 
         nbins_pull, min_pull, max_pull);
     h_err_migrad[ip] = new TH1D( 
         Form("h_err_migrad_%s",name2.c_str()), 
         Form("Migrad error for %s", util::AsStringCompact(p).c_str()), 
         nbins_err, 0., 2*prefit_err_array[ip]);
     h_err_hesse[ip] = new TH1D( 
         Form("h_err_hesse_%s",name2.c_str()), 
         Form("Hesse error for %s", util::AsStringCompact(p).c_str()), 
         nbins_err, 0., 2*prefit_err_array[ip]);
     h_err_minos[ip] = new TH1D( 
         Form("h_err_minos_%s",name2.c_str()), 
         Form("Minos error for %s", util::AsStringCompact(p).c_str()), 
         nbins_err, 0., 2*prefit_err_array[ip]);
     h_err_minos_plus[ip] = new TH1D( 
         Form("h_err_minos_plus_%s",name2.c_str()), 
         Form("Minos plus error for %s", util::AsStringCompact(p).c_str()), 
         2*nbins_err, -2*prefit_err_array[ip], 2*prefit_err_array[ip]);
     h_err_minos_minus[ip] = new TH1D( 
         Form("h_err_minos_minus_%s",name2.c_str()), 
         Form("Minos minus error for %s", util::AsStringCompact(p).c_str()), 
         2*nbins_err, -2*prefit_err_array[ip], 2*prefit_err_array[ip]);
     h_bestfit   [ip]->GetXaxis()->SetTitle(util::AsStringPaperRooTex(p).c_str());
     h_residuals [ip]->GetXaxis()->SetTitle(util::AsStringPaperRooTex(p).c_str());
     h_pull_1    [ip]->GetXaxis()->SetTitle(util::AsStringPaperRooTex(p).c_str());
     h_pull_2    [ip]->GetXaxis()->SetTitle(util::AsStringPaperRooTex(p).c_str());
     h_pull_3    [ip]->GetXaxis()->SetTitle(util::AsStringPaperRooTex(p).c_str());
     h_err_migrad[ip]->GetXaxis()->SetTitle(util::AsStringPaperRooTex(p).c_str());
     h_err_hesse [ip]->GetXaxis()->SetTitle(util::AsStringPaperRooTex(p).c_str());
     h_err_minos [ip]->GetXaxis()->SetTitle(util::AsStringPaperRooTex(p).c_str());
     h_err_minos_plus  [ip]->GetXaxis()->SetTitle(util::AsStringPaperRooTex(p).c_str());
     h_err_minos_minus [ip]->GetXaxis()->SetTitle(util::AsStringPaperRooTex(p).c_str());
  }//p

  //
  // Set up chain
  //

  TChain ch("fitnt");
  ch.Add(filenames);

  FitParamList *  fparm = 0;
  PhysParamList * iparm = 0;
  TObjArray * dat_array = 0;
  TObjArray * pdf_array = 0;

  ch.SetBranchAddress("best_fit_params", &fparm);
  ch.SetBranchAddress("input_params",    &iparm);
  ch.SetBranchAddress("data",            &dat_array);
  ch.SetBranchAddress("best_fit_pdf",    &pdf_array);

  const int nentries = ch.GetEntries();
  cout << nentries << " entries in TChain from files " << filenames << endl;
  if(nentries == 0) return;

  //
  // Loop over events & fill histograms
  //

  for(int i = 0; i < nentries; i++) 
  {
    if(i%100 == 0) { cout << "Getting entry... " << i << endl; }
    ch.GetEntry(i);

    // 
    // ignore failed fits!    
    //
    // MINUIT
    if(fitter == 0) {
      int migradminstat = fparm->GetStatusCode("migrad_min");
      int migraderrstat = fparm->GetStatusCode("migrad_err");
      int hesseerrstat  = fparm->GetStatusCode("hesse_err");
      h_migradminstat -> Fill(migradminstat);
      h_migraderrstat -> Fill(migraderrstat);
      h_hesseerrstat  -> Fill(hesseerrstat);
      if (migradminstat == 4) continue; // migrad didn't converge
      if (migraderrstat != 3 && err_mode != 3) continue; // error matrix inaccurate
    }
    // MINUIT2
    else 
    if(fitter == 1) {
      //int fitstat = fparm->GetStatusCode("??"); // currently not set by minuit2
      //if (fitstat == 9) continue; //
    }
  
    double fittime = fparm->GetFitTime();
    double chisqbf = fparm->GetBestFitChisq();
    h_fittime -> Fill(fittime);
    h_chisqbf -> Fill(chisqbf);

    TH1D * data1Rmu    = (TH1D*) dat_array->FindObject("obs_1Rmu");
    TH1D * bestfit1Rmu = (TH1D*) pdf_array->FindObject("pdf_1Rmu");
    double nsk_nom  = data1Rmu->Integral();
    double nsk_bf   = bestfit1Rmu->Integral();
    h_nsk_nom  -> Fill(nsk_nom);
    h_nsk_bf   -> Fill(nsk_bf);

    // loop over params in the fit
    for(int ip = 0; ip < nparams; ip++) 
    {
       PhysParamId_t p = params[ip];

       if(!fparm->Exists(p)) {
          if(i<100) {
            cerr << "**** Parameter " << util::AsStringCompact(p) << " doesn't exist in input parameter list" << endl;
          }
          continue;
       }

       double bestfit    = fparm->BestFitValue(p);
       double truevalue  = iparm->Value(p);

#ifdef VERBOSE
       cout << "#" << ip << " => best-fit = " << bestfit << ", true = " << truevalue;
#ifdef VERYVERBOSE
       cout << "\t" << util::AsString(p);
#endif
       cout << endl;
#endif

       double nominal    = nominal_array[ip];
       double err_prefit = prefit_err_array[ip];
       double err_migrad = fparm->MigradError(p);
       double err_hesse  = fparm->HesseError(p);
       double err_minos_plus  = fparm->PlusError(p);
       double err_minos_minus = fparm->MinusError(p);
       double err_minos  = 0;
       // err_minos  = TMath::Max(err_minos_plus, err_minos_plus);
       // err_minos  = 0.5*(err_minos_plus + err_minos_minus);
       if((bestfit - truevalue) < 0)
	 err_minos = err_minos_plus;
       else
	 err_minos = TMath::Abs(err_minos_minus);
       double err_fit    = 0;
       if      (err_mode == 1) err_fit = err_migrad;
       else if (err_mode == 2) err_fit = err_hesse;
       else if (err_mode == 3) err_fit = err_minos;
       else {
         cerr << "** Unknown err_mode = " << err_mode << endl;
         exit(1);
       }

       if(err_fit <= 0 && err_mode != 3) {
         cerr << "** err_fit[p = " << util::AsStringCompact(p) 
              << "; err_mode = " << err_mode << "] = " << err_fit << endl;
         exit(1);
       }

       double residuals  = (bestfit-truevalue);
       double pull_1     = (bestfit-truevalue)/err_fit;
       double pull_2     = (bestfit-truevalue)/TMath::Sqrt(err_prefit*err_prefit-err_fit*err_fit);
       double pull_3     = (bestfit-nominal)/err_fit;

       h_bestfit    [ip] -> Fill(bestfit);
       h_corr_inputdcp_residual[ip] -> Fill(iparm->Value(conv::kPhypDeltaCP),         residuals);
       h_corr_inputdcp_bestfit [ip] -> Fill(iparm->Value(conv::kPhypDeltaCP),         bestfit);
       h_corr_fitmix23_residual[ip] -> Fill(fparm->BestFitValue(conv::kPhypMixing23), residuals);
       h_corr_fitmix23_bestfit [ip] -> Fill(fparm->BestFitValue(conv::kPhypMixing23), bestfit);
       h_residuals  [ip] -> Fill(residuals);
       h_pull_1     [ip] -> Fill(pull_1);
       h_pull_2     [ip] -> Fill(pull_2);
       h_pull_3     [ip] -> Fill(pull_3);
       h_err_migrad [ip] -> Fill(err_migrad);
       h_err_hesse  [ip] -> Fill(err_hesse);
       h_err_minos  [ip] -> Fill(err_minos);
       h_err_minos_plus   [ip] -> Fill(err_minos_plus);
       h_err_minos_minus  [ip] -> Fill(err_minos_minus);


    }//p

  }//end loop over chain

  //
  // Summarise & fit
  //
  // loop over params in the fit
  TF1 *fit = 0;
  for(int ip = 0, ipp = 0; ip < nparams; ip++) {
    PhysParamId_t p = params[ip];
    if(p == conv::kPhypDm23sq || p == conv::kPhypMixing23 || p == conv::kNuisXSecSpecFunc)
      continue;
    TString label;
    if(ipp%2)
      label = Form("%s                         ",util::AsStringPaperRooTex(p).c_str());
    else
      label = Form("%s",util::AsStringPaperRooTex(p).c_str());
    h_summary_nuis->GetXaxis()->SetBinLabel(ipp+1, label.Data());
    if(!do_gaus_fits) {
      h_summary_nuis->SetBinContent(ipp+1, h_pull_1[ip]->GetMean());
      h_summary_nuis->SetBinError(  ipp+1, h_pull_1[ip]->GetRMS());
    }
    else {
      h_pull_1[ip]->Draw();
      cout << util::AsString(p) << endl;
      h_pull_1[ip]->Fit("gaus");
      fit = h_pull_1[ip]->GetFunction("gaus");
      if(fit) {
	h_summary_nuis->SetBinContent(ipp+1, fit->GetParameter(1));
	h_summary_nuis->SetBinError(  ipp+1, fit->GetParameter(2));
      }
    }//do_gaus_fits
    ipp++;
  }

  //
  // Draw
  //
  TCanvas * c = 0;

  // >> Plot etc
  c = new TCanvas(); h_migradminstat->Draw();
  save_canvas(c, "migrad_minimisation_status", tag);
  c = new TCanvas(); h_migraderrstat->Draw();
  save_canvas(c, "migrad_error_status", tag);
  c = new TCanvas(); h_hesseerrstat->Draw();
  save_canvas(c, "hesse_error_status", tag);
  c = new TCanvas(); h_fittime->Draw();
  save_canvas(c, "fit_time", tag);
  c = new TCanvas(); h_chisqbf->Draw();
  save_canvas(c, "bf_chisq", tag);
  c = new TCanvas();
  h_nsk_bf->SetLineColor(kRed);
  if(h_nsk_nom->GetMaximum() > h_nsk_bf->GetMaximum()) {
    h_nsk_nom->Draw();
    h_nsk_bf->Draw("same");
  }
  else {
    h_nsk_bf->Draw();
    h_nsk_nom->Draw("same");
  }
  save_canvas(c, "nsk", tag);

  // >> Plot residuals 
  // put all osc params in one canvas (1 columnn x nosc rows)
  // put all syst params in few canvases with (nsystcol columns x nsystrow rows)
  int ipc=0;//plot nu in curr canvas
  int ipcmax=0;//max nu of plot in curr canvas
  int icanvas=0;//number of the canvas
  //int iposc=0;//plot nu in osc canvas
  //TCanvas * cosc = new TCanvas();
  //cosc->Divide(1,nosc);
  for(int ip = 0; ip < nparams; ip++) 
  {
     PhysParamId_t p = params[ip];
     if(ipc==0 && ip==0) {
        c = new TCanvas();
        c->Divide(1,nosc);
        ipcmax = nosc;
     }
     if(ipc==0 && ip>=nosc) {
        c = new TCanvas();
        c->Divide(nsystcol,nsystrow);
        ipcmax = nsystcol*nsystrow;
     }
     c->cd(ipc+1);
     h_residuals[ip]->Draw();
     ipc++;
     if(ipc >= ipcmax || ip == nparams-1) { 
       ipc=0; c->Update();
       save_canvas(c,Form("residual_%d",icanvas),tag);
       icanvas++;
     }
  }//p
  //save_canvas(cosc,"residual_osc",tag);

  // >> Plot pulls (all pull 1,2 and 3 on same plot), same logic as before
  ipc=0;//plot nu in curr canvas
  ipcmax=0;//max nu of plot in curr canvas
  icanvas=0;//number of the canvas
  for(int ip = 0; ip < nparams; ip++) 
  {
     PhysParamId_t p = params[ip];
     if(ipc==0 && ip==0) {
        c = new TCanvas();
        c->Divide(1,nosc);
        ipcmax = nosc;
     }
     if(ipc==0 && ip>=nosc) {
        c = new TCanvas();
        c->Divide(nsystcol,nsystrow);
        ipcmax = nsystcol*nsystrow;
     }
     c->cd(ipc+1);
     h_pull_1[ip]->SetLineColor(kBlack);
     h_pull_2[ip]->SetLineColor(kRed);
     h_pull_3[ip]->SetLineColor(kGreen);
     h_pull_1[ip]->Draw();
     h_pull_2[ip]->Draw("same");
     h_pull_3[ip]->Draw("same");
     ipc++;
     if(ipc >= ipcmax || ip == nparams-1) { 
       ipc=0; c->Update();
       save_canvas(c,Form("pull_all_%d",icanvas),tag);
       icanvas++;
     }
  }//p
  // >> Plot pulls, same logic as before
  // >> Plot pulls (only pull 1), same logic as before
  ipc=0;//plot nu in curr canvas
  ipcmax=0;//max nu of plot in curr canvas
  icanvas=0;//number of the canvas
  for(int ip = 0; ip < nparams; ip++) 
  {
     PhysParamId_t p = params[ip];
     if(ipc==0 && ip==0) {
        c = new TCanvas();
        c->Divide(1,nosc);
        ipcmax = nosc;
     }
     if(ipc==0 && ip>=nosc) {
        c = new TCanvas();
        c->Divide(nsystcol,nsystrow);
        ipcmax = nsystcol*nsystrow;
     }
     c->cd(ipc+1);
     h_pull_1[ip]->SetLineColor(kBlack);
     h_pull_1[ip]->Draw();
     ipc++;
     if(ipc >= ipcmax || ip == nparams-1) { 
       ipc=0; c->Update();
       save_canvas(c,Form("pull_1_%d",icanvas),tag);
       icanvas++;
     }
  }//p
  // >> Plot pulls, same logic as before

#if FITPARMLIST == -3
  c = new TCanvas();
  h_pull_1[0]->SetLineColor(kBlack);
  h_pull_1[0]->Draw();
  save_canvas(c,"pull_1_escale",tag);
#endif

  //plot summary at the end, because need to turn stat box off & also increase margin
  c = new TCanvas(); c->SetBottomMargin(0.35);
  gStyle->SetOptStat(0);
  h_summary_nuis->GetXaxis()->LabelsOption("v");
  h_summary_nuis->Draw();
  save_canvas(c, "summary", tag);

  //
  // Store
  //
  TFile foutput(Form("./valor_fit_bias_study_%s.root",tag),"recreate");
  h_migradminstat       -> Write();
  h_migraderrstat -> Write();
  h_hesseerrstat  -> Write();
  h_fittime       -> Write();
  h_chisqbf       -> Write();
  h_nsk_nom       -> Write();
  h_nsk_bf        -> Write();
  h_summary_nuis  -> Write();
  for(int ip = 0; ip < nparams; ip++) 
    {
      h_bestfit    [ip] -> Write();
      h_corr_inputdcp_residual[ip] -> Write();
      h_corr_inputdcp_bestfit [ip] -> Write();
      h_corr_fitmix23_residual[ip] -> Write();
      h_corr_fitmix23_bestfit [ip] -> Write();
      h_residuals  [ip] -> Write();
      h_pull_1     [ip] -> Write();
      h_pull_2     [ip] -> Write();
      h_pull_3     [ip] -> Write();
      h_err_migrad [ip] -> Write();
      h_err_hesse  [ip] -> Write();
      h_err_minos  [ip] -> Write();
      h_err_minos_plus   [ip] -> Write();
      h_err_minos_minus  [ip] -> Write();
    }//ip
  foutput.Close();
}

