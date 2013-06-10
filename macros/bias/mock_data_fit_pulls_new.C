//
// ------------------------------------------------------------------------------------------------
// Read macros/fit/run_fit.C ROOT outputs with input and best-fit values for fits to toy-MC expts.
// Plot various residual and pull distributions for fit parameters.
// ------------------------------------------------------------------------------------------------
//

void mock_data_fit_pulls_new
(
  const char * filenames,                                              /* add these filenames to the tree */
  const char * tag,                                                    /* tag for the output files */
  int          err_mode             = 1,                               /* which errors to use for pulls (1: migrad, 2 hesse, 3 minos) */
  int          fitter               = 0,                               /* which fitter has been used (0: minuit, 1: minuit2) */
  const char * config_file          = "analysis_configuration.xml",    /* file to get the pre errors from */
  bool         do_gaus_fits         = false
)
{
  using namespace std;
  using namespace std::string;
  using namespace t2kana;
  using namespace t2kana::conv;

  // First get the expected systematic variations
//gStyle->SetOptFit(1);         // show fit results
  gStyle->SetOptStat(1100);     // show mean & RMS only
//gStyle->SetStatStyle(0);      // make stat boxes see through
  gStyle->SetStatBorderSize(0); //
  gStyle->SetTextFont(42);      // less bold

  //
  // Parameters in the fit
  //
  const int nparams = 44;

  PhysParamId_t params[nparams] = {
     kPhypMixing23,             
   //kPhypMixing13,  
   //kPhypMixing12,   
   //kPhypDeltaCP, 
     kPhypDm23sq,  
   //kPhypDm12sq, 
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
     kNuisXSecNormNCcoh,         
     kNuisXSecNormNCoth,         
     kNuisXSecNormNuEToNuMu,     
     kNuisXSecNormAntiNuToNu,    
     kNuisSuperkEnergyScale, 
     kNuisSKDetFSI000,
     kNuisSKDetFSI001,
     kNuisSKDetFSI002,
     kNuisSKDetFSI003,
     kNuisSKDetFSI004,
     kNuisSKDetFSI005
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

  double escale_err = util::XML2Dbl(xml,root_node,
    "detector_systematics/superk/energy_scale_relative_error/one_mu_ring");

  double mix23_err  = kGlobFitSinsq2Theta23 * kGlobFitErrSinsq2Theta23;
  double dm32sq_err = kGlobFitDm23sq        * kGlobFitErrDm23sq;

  double prefit_err_array[nparams] = 
  {
     mix23_err,                    // kPhypMixing23
                                   // kPhypMixing13
                                   // kPhypMixing12
                                   // kPhypDeltaCP
     dm32sq_err,                   // kPhypDm23sq
                                   // kPhypDm12sq
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
     1,                            // kNuisXSecNormNCcoh         
     1,                            // kNuisXSecNormNCoth         
     1,                            // kNuisXSecNormNuEToNuMu     
     1,                            // kNuisXSecNormAntiNuToNu
     escale_err,                   // kNuisSuperkEnergyScale
     TMath::Sqrt(skfsicov(0,0)),   // kNuisSKDetFSI000
     TMath::Sqrt(skfsicov(1,1)),   // kNuisSKDetFSI001
     TMath::Sqrt(skfsicov(2,2)),   // kNuisSKDetFSI002
     TMath::Sqrt(skfsicov(3,3)),   // kNuisSKDetFSI003
     TMath::Sqrt(skfsicov(4,4)),   // kNuisSKDetFSI004
     TMath::Sqrt(skfsicov(5,5))    // kNuisSKDetFSI005
  };


  double nominal_array[nparams] = 
  {
     kGlobFitSinsq2Theta23,    // kPhypMixing23
                               // kPhypMixing13
                               // kPhypMixing12
                               // kPhypDeltaCP
     kGlobFitDm23sq,           // kPhypDm23sq
                               // kPhypDm12sq
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
     0,                        // kNuisXSecNormNCcoh         
     0,                        // kNuisXSecNormNCoth         
     0,                        // kNuisXSecNormNuEToNuMu     
     0,                        // kNuisXSecNormAntiNuToNu
     0,                        // kNuisSuperkEnergyScale
     1,                        // kNuisSKDetFSI000
     1,                        // kNuisSKDetFSI001
     1,                        // kNuisSKDetFSI002
     1,                        // kNuisSKDetFSI003
     1,                        // kNuisSKDetFSI004
     1                         // kNuisSKDetFSI005
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

  TH1I * h_fitstat        = new TH1I("h_fitstat",       "fit status",           11, -1,   10);
  TH1I * h_migradfitstat  = new TH1I("h_migradfitstat", "migrad fit status",    14, -1,   13);
  TH1I * h_hessefitstat   = new TH1I("h_hessefitstat",  "hesse fit status",     14, -1,   13);
  TH1D * h_fittime        = new TH1D("h_fittime",       "fit time",           1000,  0, 1000);
  TH1D * h_nsk_nom        = new TH1D("h_nsk_nom",       "nsk(1Rmu) nom",       500,  0,  500);
  TH1D * h_nsk_bf         = new TH1D("h_nsk_bf",        "nsk(1Rmu) best-fit",  500,  0,  500);
  TH1D * h_dchisq         = new TH1D("h_dchisq",        "h_dchisq",            300, -5,   25);

  TH1D * h_bestfit    [nparams] = { 0 }; // best-fit param distributions
  TH1D * h_residuals  [nparams] = { 0 }; // residuals
  TH1D * h_pull_1     [nparams] = { 0 }; // pulls - 1st definition (as in T2K-TN-141v6, p.320)
  TH1D * h_pull_2     [nparams] = { 0 }; //  ...    2nd     ...          ...
  TH1D * h_pull_3     [nparams] = { 0 }; //  ...    3rd     ...          ...
  TH1D * h_err_migrad [nparams] = { 0 }; // migrad err distribution
  TH1D * h_err_hesse  [nparams] = { 0 }; // hesse  err distribution
  TH1D * h_err_minos  [nparams] = { 0 }; // minos  err distribution

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
         Form("Best-fit distribution for %s", util::AsString(p).c_str()), 
         nbins_pull, min_pull, max_pull);
     h_residuals[ip] = new TH1D( 
         Form("h_residuals_%s",name2.c_str()), 
         Form("Residuals for %s", util::AsString(p).c_str()), 
         nbins_pull, -5*prefit_err_array[ip], +5*prefit_err_array[ip]);
     h_pull_1[ip] = new TH1D( 
         Form("h_pull_1_%s",name2.c_str()), 
         Form("Pull_1 for %s", util::AsString(p).c_str()), 
         nbins_pull, min_pull, max_pull);
     h_pull_2[ip] = new TH1D( 
         Form("h_pull_2_%s",name2.c_str()), 
         Form("Pull_2 for %s", util::AsString(p).c_str()), 
         nbins_pull, min_pull, max_pull);
     h_pull_3[ip] = new TH1D( 
         Form("h_pull_3_%s",name2.c_str()), 
         Form("Pull_3 for %s", util::AsString(p).c_str()), 
         nbins_pull, min_pull, max_pull);
     h_err_migrad[ip] = new TH1D( 
         Form("h_err_migrad_%s",name2.c_str()), 
         Form("Migrad error for %s", util::AsString(p).c_str()), 
         nbins_err, 0., 2*prefit_err_array[ip]);
     h_err_hesse[ip] = new TH1D( 
         Form("h_err_hesse_%s",name2.c_str()), 
         Form("Hesse error for %s", util::AsString(p).c_str()), 
         nbins_err, 0., 2*prefit_err_array[ip]);
     h_err_minos[ip] = new TH1D( 
         Form("h_err_minos_%s",name2.c_str()), 
         Form("Minos error for %s", util::AsString(p).c_str()), 
         nbins_err, 0., 2*prefit_err_array[ip]);
  }//p


  //
  // Set up chain
  //

  TChain ch("fitnt");
  ch.Add(filenames);

  FitParamList *  fparm = 0;
  PhysParamList * iparm = 0;

  ch.SetBranchAddress("best_fit_params", &fparm);
  ch.SetBranchAddress("input_params",    &iparm);

  const int nentries = ch.GetEntries();
  cout << nentries << " entries in TChain" << endl;
  if(nentries == 0) return;

  //
  // Loop over events & fill histograms
  //

  for(int i = 0; i < nentries; i++) 
  {
    if(i%1 == 0) { cout << "Getting entry... " << i << endl; }
    ch.GetEntry(i);

    int    fitstat        = 0; // get from where?
    int    migradfitstat  = 0; // get from where?
    int    hessefitstat   = 0; // get from where?
    double fittime        = 0; // get from where?

    h_fitstat        -> Fill(fitstat);
    h_migradfitstat  -> Fill(migradfitstat);
    h_hessefitstat   -> Fill(hessefitstat);
    h_fittime        -> Fill(fittime);

    //
    // ignore failed fits!    
    //
/*
    // MINUIT
    if(fitter == 0) {
      if (migradfitstat == 4) continue; // migrad didn't converge
      if (fitstat       != 3) continue; // error matrix inaccurate
    }
    // MINUIT2
    else 
    if(fitter == 1) {
      if (fitstat == 9) continue; //
    }
*/
  
    double nsk_nom  = 0; // get from where?
    double nsk_bf   = 0; // get from where?
    double dchisq   = 0; // get from where?

    h_nsk_nom  -> Fill(nsk_nom);
    h_nsk_bf   -> Fill(nsk_bf);
    h_dchisq   -> Fill(dchisq);

    // loop over params in the fit
    for(int ip = 0; ip < nparams; ip++) 
    {
       PhysParamId_t p = params[ip];
  
       double bestfit    = fparm->BestFitValue(p);
       double true       = iparm->Value(p);

       cout << ip << "best: " << bestfit << ", true: " << true << endl;

       double nominal    = nominal_array[ip];
       double err_prefit = prefit_err_array[ip];
       double err_migrad = 0;
       double err_hesse  = 0;
       double err_minos  = 0;
       double err_fit    = 0;

       if      (err_mode == 1) err_fit = err_migrad;
       else if (err_mode == 2) err_fit = err_hesse;
       else if (err_mode == 3) err_fit = err_minos;
       else {
         cerr << "** Unknown err_mode = " << err_mode << endl;
         exit(1);
       }

       if(err_fit <= 0) {
         cerr << "** err_fit[p = " << util::AsStringCompact(p) 
              << "; err_mode = " << err_mode << "] = " << err_fit << endl;
//         exit(1);
         err_fit = 1;//HERE!
       }

       double residuals  = (bestfit-true);
       double pull_1     = (bestfit-true)/err_fit;
       double pull_2     = (bestfit-true)/TMath::Sqrt(err_prefit*err_prefit-err_fit*err_fit);
       double pull_3     = (bestfit-nominal)/err_fit;

       h_bestfit    [ip] -> Fill(bestfit);
       h_residuals  [ip] -> Fill(residuals);
       h_pull_1     [ip] -> Fill(pull_1);
       h_pull_2     [ip] -> Fill(pull_2);
       h_pull_3     [ip] -> Fill(pull_3);
       h_err_migrad [ip] -> Fill(err_migrad);
       h_err_hesse  [ip] -> Fill(err_hesse);
       h_err_minos  [ip] -> Fill(err_minos);
    }//p

  }//end loop over chain

  //
  // Draw
  //

  new TCanvas(); h_fitstat->Draw();
  new TCanvas(); h_migradfitstat->Draw();
  new TCanvas(); h_hessefitstat->Draw();
  new TCanvas(); h_fittime->Draw();
  new TCanvas(); h_nsk_nom->Draw();
  new TCanvas(); h_nsk_bf->Draw();
  new TCanvas(); h_dchisq->Draw();

  for(int ip = 0; ip < nparams; ip++) 
  {
     PhysParamId_t p = params[ip];
 
  }//p

  //
  // Store
  //
  TFile foutput("./valor_fit_bias_study.root","recreate");
  h_fitstat       -> Write();
  h_migradfitstat -> Write();
  h_hessefitstat  -> Write();
  h_fittime       -> Write();
  h_nsk_nom       -> Write();
  h_nsk_bf        -> Write();
  h_dchisq        -> Write();
  for(int ip = 0; ip < nparams; ip++) 
  {
     h_bestfit    [ip] -> Write();
     h_residuals  [ip] -> Write();
     h_pull_1     [ip] -> Write();
     h_pull_2     [ip] -> Write();
     h_pull_3     [ip] -> Write();
     h_err_migrad [ip] -> Write();
     h_err_hesse  [ip] -> Write();
     h_err_minos  [ip] -> Write();
  }
  foutput.Close();
}

