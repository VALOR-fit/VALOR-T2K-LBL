//
// ------------------------------------------------------------------------------------------------
// 
//
// ------------------------------------------------------------------------------------------------
//

void mock_data_fit_pulls
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
  using namespace t2kana;

  // First get the expected systematic variations
  //gStyle->SetOptFit(1); // show fit results
  gStyle->SetOptStat(1100); // show mean & RMS only
  //gStyle->SetStatStyle(0); //make stat boxes see through
  gStyle->SetStatBorderSize(0);
  gStyle->SetTextFont(42); //less bold

  //
  // Parameters in the fit
  //
  const int nparams = 46;
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
     kNuisBANFFMECNorm,        
     kNuisBANFFPDDDoubleCount, 
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
  // Get 1-sigma errors for all systematics
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

  //SK
  double preerr_skdet_energyscale = 1;//util::XML2Dbl(xml,root_node,
     //"detector_systematics/superk/energy_scale_relative_error/one_mu_ring");
  double preerr_skdet_ccqeoth = 1;//util::XML2Dbl(xml,root_node,
      //"detector_systematics/superk/efficiency_relative_error/one_mu_ring/ccqe/other_error");
  //TMatrixD *SKDetEff1RmuCCQERCCovarRelErr = util::XML2TMatrixD(xml,root_node,
  //  "detector_systematics/superk/efficiency_relative_error/one_mu_ring/ccqe/rc_error");
  double preerr_skdet_ccqerc0 = 1;//SKDetEff1RmuCCQERCCovarRelErr[0];
  double preerr_skdet_ccqerc1 = 1;//SKDetEff1RmuCCQERCCovarRelErr[1];
  double preerr_skdet_ccqerc2 = 1;//SKDetEff1RmuCCQERCCovarRelErr[2];
  double preerr_skdet_ccnqe   = 1;//util::XML2Dbl(xml,root_node,
      //"detector_systematics/superk/efficiency_relative_error/one_mu_ring/ccnonqe");
  double preerr_skdet_nc = 1;//util::XML2Dbl(xml,root_node,
      //"detector_systematics/superk/efficiency_relative_error/one_mu_ring/nc");
  double preerr_skdet_nuecc = 1;//util::XML2Dbl(xml,root_node, 
      //"detector_systematics/superk/efficiency_relative_error/one_mu_ring/nuecc");

  //BANFF
  TMatrixD * BANFFCovMtrxP = util::XML2TMatrixD(xml,root_node, "banff/covariance_matrix");
  const TMatrixD & BANFFCovMtrx = *BANFFCovMtrxP;
  double preerr_banff_00 = sqrt(BANFFCovMtrx(0,0));
  double preerr_banff_01 = sqrt(BANFFCovMtrx(1,1));
  double preerr_banff_02 = sqrt(BANFFCovMtrx(2,2));
  double preerr_banff_03 = sqrt(BANFFCovMtrx(3,3));
  double preerr_banff_04 = sqrt(BANFFCovMtrx(4,4));
  double preerr_banff_05 = sqrt(BANFFCovMtrx(5,5));
  double preerr_banff_06 = sqrt(BANFFCovMtrx(6,6));
  double preerr_banff_07 = sqrt(BANFFCovMtrx(7,7));
  double preerr_banff_08 = sqrt(BANFFCovMtrx(8,8));
  double preerr_banff_09 = sqrt(BANFFCovMtrx(9,9));
  double preerr_banff_10 = sqrt(BANFFCovMtrx(10,10));
  double preerr_banff_11 = sqrt(BANFFCovMtrx(11,11));
  double preerr_banff_12 = sqrt(BANFFCovMtrx(12,12));
  double preerr_banff_13 = sqrt(BANFFCovMtrx(13,13));
  double preerr_banff_14 = sqrt(BANFFCovMtrx(14,14));
  double preerr_banff_15 = sqrt(BANFFCovMtrx(15,15));
  double preerr_banff_maqe = sqrt(BANFFCovMtrx(16,16));
  double preerr_banff_mares = sqrt(BANFFCovMtrx(17,17));
  double preerr_banff_ccqeE0 = sqrt(BANFFCovMtrx(18,18));
  double preerr_banff_ccqeE1 = sqrt(BANFFCovMtrx(19,19));
  double preerr_banff_ccqeE2 = sqrt(BANFFCovMtrx(20,20));
  double preerr_banff_cc1piE0 = sqrt(BANFFCovMtrx(21,21));
  double preerr_banff_cc1piE1 = sqrt(BANFFCovMtrx(22,22));
  //bestfit
  TVectorD * BANFFBFP = util::XML2TVectorD(xml,root_node, "banff/best_fit_values");
  const TVectorD & BANFFBF = *BANFFBFP;
  double prebf_banff_00 = BANFFBF[0];
  double prebf_banff_01 = BANFFBF[1];
  double prebf_banff_02 = BANFFBF[2];
  double prebf_banff_03 = BANFFBF[3];
  double prebf_banff_04 = BANFFBF[4];
  double prebf_banff_05 = BANFFBF[5];
  double prebf_banff_06 = BANFFBF[6];
  double prebf_banff_07 = BANFFBF[7];
  double prebf_banff_08 = BANFFBF[8];
  double prebf_banff_09 = BANFFBF[9];
  double prebf_banff_10 = BANFFBF[10];
  double prebf_banff_11 = BANFFBF[11];
  double prebf_banff_12 = BANFFBF[12];
  double prebf_banff_13 = BANFFBF[13];
  double prebf_banff_14 = BANFFBF[14];
  double prebf_banff_15 = BANFFBF[15];
  double prebf_banff_maqe = BANFFBF[16];
  double prebf_banff_mares = BANFFBF[17];
  double prebf_banff_ccqeE0 = BANFFBF[18];
  double prebf_banff_ccqeE1 = BANFFBF[19];
  double prebf_banff_ccqeE2 = BANFFBF[20];
  double prebf_banff_cc1piE0 = BANFFBF[21];
  double prebf_banff_cc1piE1 = BANFFBF[22];

  //FSI
  TMatrixD *FSICovMtrxP = util::XML2TMatrixD(xml,root_node, "fsi/covariance_matrix");
  const TMatrixD & FSICovMtrx = *FSICovMtrxP;
  double preerr_fsi_000 = sqrt(FSICovMtrx(0,0));
  double preerr_fsi_001 = sqrt(FSICovMtrx(1,1));
  double preerr_fsi_002 = sqrt(FSICovMtrx(2,2));
  double preerr_fsi_003 = sqrt(FSICovMtrx(3,3));
  double preerr_fsi_004 = sqrt(FSICovMtrx(4,4));
  double preerr_fsi_005 = sqrt(FSICovMtrx(5,5));

  //XSec
  double preerr_xsec_normcccoh = 1;//util::XML2Dbl(xml,root_node,"xsec/normalisations/cccoh");
  double preerr_xsec_normnc1piC = 1;//util::XML2Dbl(xml,root_node,"xsec/normalisations/nc1piC");
  double preerr_xsec_normnccoh = 1;//util::XML2Dbl(xml,root_node,"xsec/normalisations/nccoh");
  double preerr_xsec_normncoth = 1;//util::XML2Dbl(xml,root_node,"xsec/normalisations/ncoth");
  double preerr_xsec_normnuetonumu = 1;//util::XML2Dbl(xml,root_node,"xsec/normalisations/numu_to_nue");
  double preerr_xsec_normantinutonu = 1;//util::XML2Dbl(xml,root_node,"");
  double preerr_xsec_pf = 1;//util::XML2Dbl(xml,root_node,"xsec/splines/pf");
  double preerr_xsec_sf = 1;//util::XML2Dbl(xml,root_node,"xsec/splines/sf");
  double preerr_xsec_binde = 1;//util::XML2Dbl(xml,root_node,"xsec/splines/bindE");
  double preerr_xsec_wshape = 1;//util::XML2Dbl(xml,root_node,"xsec/splines/wshape");
  double preerr_xsec_ccothshape = 1;//util::XML2Dbl(xml,root_node,"xsec/splines/ccothshape");
  double preerr_xsec_cc1pishape = 1;//util::XML2Dbl(xml,root_node,"xsec/splines/cc1pishape");
  double preerr_xsec_pilessdelta = 1;//util::XML2Dbl(xml,root_node,"xsec/splines/pilessdelta");

  // Release memory
  xml->FreeDoc(xmldoc);
  delete xml;

  //
  // Book histograms
  //

  const int    nbins_pull = 100;
  const double min_pull   = -5;
  const double max_pull   =  5;

  TH1I * h_fitstat        = new TH1I("h_fitstat",       "fit status",         11, -1,   10);
  TH1I * h_migradfitstat  = new TH1I("h_migradfitstat", "migrad fit status",  14, -1,   13);
  TH1I * h_hessefitstat   = new TH1I("h_hessefitstat",  "hesse fit status",   14, -1,   13);
  TH1D * h_fittime        = new TH1D("h_fittime",       "fit time",         1000,  0, 1000);
  TH1D * h_nsk_nom        = new TH1D("h_nsk_nom",       "nsk nom",           150,  0,  150);
  TH1D * h_nsk_bf         = new TH1D("h_nsk_bf",        "nsk best-fit",      150,  0,  150);
  TH1D * h_dchisq         = new TH1D("h_dchisq",        "h_dchisq",          300, -5,   25);

  TH1D * h_bestfit   [nparams] = { 0 }; //
  TH1D * h_pull_1    [nparams] = { 0 }; //
  TH1D * h_pull_2    [nparams] = { 0 }; //
  TH1D * h_pull_3    [nparams] = { 0 }; //
  TH1D * h_err_minos [nparams] = { 0 }; //
  TH1D * h_err_hesse [nparams] = { 0 }; //

  for(int ip = 0; ip < nparams; ip++) 
  {
     PhysParamId_t p params[ip];

     h_bestfit[ip] = new TH1D( 
         Form("h_bestfit_%s",util::AsStringCompact(p).c_str()), 
         Form("Best-fit distribution for %s", util::AsString(p).c_str()), 
         nbins_pull, min_pull, max_pull);
     h_pull_1[ip] = new TH1D( 
         Form("h_pull_1_%s",util::AsStringCompact(p).c_str()), 
         Form("Pull_1 for %s", util::AsString(p).c_str()), 
         nbins_pull, min_pull, max_pull);
     h_pull_2[ip] = new TH1D( 
         Form("h_pull_2_%s",util::AsStringCompact(p).c_str()), 
         Form("Pull_2 for %s", util::AsString(p).c_str()), 
         nbins_pull, min_pull, max_pull);
     h_pull_3[ip] = new TH1D( 
         Form("h_pull_3_%s",util::AsStringCompact(p).c_str()), 
         Form("Pull_3 for %s", util::AsString(p).c_str()), 
         nbins_pull, min_pull, max_pull);
     h_err_minos[ip] = new TH1D( 
         Form("h_err_minos_%s",util::AsStringCompact(p).c_str()), 
         Form("Minos error for %s", util::AsString(p).c_str()), 
         nbins_pull, min_pull, max_pull);
     h_err_hesse[ip] = new TH1D( 
         Form("h_err_hesse_%s",util::AsStringCompact(p).c_str()), 
         Form("Hesse error for %s", util::AsString(p).c_str()), 
         nbins_pull, min_pull, max_pull);
  }//p


  //
  // Set up chain
  //

  TChain ch("fitnt");
  ch.Add(filenames);

  FitParamList *  fparm = 0;
  PhysParamList * iparm = 0;

//  TTree * fitnt = (TTree*)f.Get("fitnt");
  ch.SetBranchAddress("best_fit_params", &fparm);
  ch.SetBranchAddress("input_params",    &iparm);

  cout << ch.GetEntries() << " entries in TChain" << endl;
  if(ch.GetEntries() == 0) return;

  // Loop over events & fill histograms
  for(int i = 0; i < ch.GetEntries(); i++) {

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

    // MINUIT
    if(fitter == 0) {
      if (migradfitstatus == 4) continue; // migrad didn't converge
      if (fitstatus       != 3) continue; // error matrix inaccurate
    }
    // MINUIT2
    else 
    if(fitter == 1) {
      if (fitstatus == 9) continue; //
    }
  
    double nsk_nom  = 0; // get from where?
    double nsk_bf   = 0; // get from where?
    double dchisq   = 0; // get from where?

    h_nsk_nom  -> Fill(nsk_nom);
    h_nsk_bf   -> Fill(nsk_bf);
    h_dchisq   -> Fill(dchisq);

    for(int ip = 0; ip < nparams; ip++) 
    {
       PhysParamId_t p params[ip];
  
       double bestfit    = 0;
       double true       = 0;
       double err_prefit = 0;
       double err_minos  = 0;
       double err_hesse  = 0;
       double pull_1     = 0;
       double pull_2     = 0;
       double pull_3     = 0;

       h_bestfit   [ip] -> Fill();
       h_pull_1    [ip] -> Fill();
       h_pull_2    [ip] -> Fill();
       h_pull_3    [ip] -> Fill();
       h_err_minos [ip] -> Fill();
       h_err_hesse [ip] -> Fill();
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
     PhysParamId_t p params[ip];
  
       double bestfit    = 0;
       double true       = 0;
       double err_prefit = 0;
       double err_minos  = 0;
       double err_hesse  = 0;
       double pull_1     = 0;
       double pull_2     = 0;
       double pull_3     = 0;

       h_bestfit   [ip] -> Fill();
       h_pull_1    [ip] -> Fill();
       h_pull_2    [ip] -> Fill();
       h_pull_3    [ip] -> Fill();
       h_err_minos [ip] -> Fill();
       h_err_hesse [ip] -> Fill();
    }//p




  TH1D * meansummary = new TH1D("mean","mean",100,0,100);
  TH1D * sigmasummary = new TH1D("sigma","sigma",100,0,100);
  int ifit = 1;
  TF1 *fit = 0;
  //gStyle->SetOptFit();
  

  TCanvas *coscpull = new TCanvas(); coscpull->Divide(1,2);
  coscpull->cd(1); hpull_mixing23->Draw(); hpull_mixing23->Fit("gaus");
  //fit = hpull_mixing23->GetFunction("gaus");
  if(fit) {
    meansummary->GetXaxis()->SetBinLabel(ifit,"mixing23"); sigmasummary->GetXaxis()->SetBinLabel(ifit,"mixing23");
    meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
    sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2));
    ifit++;
  }
  coscpull->cd(2); hpull_dm23sq->Draw(); hpull_dm23sq->Fit("gaus");
  //fit = hpull_dm23sq->GetFunction("gaus");
  if(fit) {
    meansummary->GetXaxis()->SetBinLabel(ifit,"dm32sq"); sigmasummary->GetXaxis()->SetBinLabel(ifit,"dm32sq");
    meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
    sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2));
    ifit++;
  }
  coscpull->SaveAs(Form("%s-%f_osc_pull.pdf", tag, input_mixing23));

  TCanvas *coscerr = new TCanvas(); coscerr->Divide(1,2);
  coscerr->cd(1); herr_mixing23->Draw();
  coscerr->cd(2); herr_dm23sq->Draw();

  if(fit_systematics_mode >= 1 && do_gaus_fits) {

    TCanvas *cbanff = new TCanvas(); cbanff->Divide(1,2);
    cbanff->cd(1); hpull_banff_maqe->Draw(); hpull_banff_maqe->Fit("gaus");
    fit = hpull_banff_maqe->GetFunction("gaus");
    if(fit) {
      meansummary->GetXaxis()->SetBinLabel(ifit,"maqe"); sigmasummary->GetXaxis()->SetBinLabel(ifit,"maqe");
      meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
      sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2));
      ifit++;
    }
    cbanff->cd(2); hpull_banff_ccqeE0->Draw(); hpull_banff_ccqeE0->Fit("gaus");
    fit = hpull_banff_ccqeE0->GetFunction("gaus");
    if(fit) {
      meansummary->GetXaxis()->SetBinLabel(ifit,"ccqeE0"); sigmasummary->GetXaxis()->SetBinLabel(ifit,"ccqeE0");
      meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
      sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2));
      ifit++;
    }
    TCanvas *cxsec = new TCanvas(); cxsec->Divide(1,2);
    cxsec->cd(1); hpull_xsec_sf->Draw();
    /* hpull_xsec_sf->Fit("gaus");
    if(fit) {
      fit = hpull_xsec_sf->GetFunction("gaus");
      meansummary->GetXaxis()->SetBinLabel(ifit,"SF"); sigmasummary->GetXaxis()->SetBinLabel(ifit,"SF");
      meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
      sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2));
      ifit++;
    }
    */
    cxsec->cd(2); hpull_xsec_pilessdelta->Draw(); hpull_xsec_pilessdelta->Fit("gaus");
    fit = hpull_xsec_pilessdelta->GetFunction("gaus");
    if(fit) {
      meansummary->GetXaxis()->SetBinLabel(ifit,"piless"); sigmasummary->GetXaxis()->SetBinLabel(ifit,"piless");
      meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
      sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2));
      ifit++;
    }
    TCanvas *cextra = new TCanvas(); cextra->Divide(1,2);
    cextra->cd(1); hpull_xsec_binde->Draw(); hpull_xsec_binde->Fit("gaus");
    fit = hpull_xsec_binde->GetFunction("gaus");
    if(fit) {
      meansummary->GetXaxis()->SetBinLabel(ifit,"bindE"); sigmasummary->GetXaxis()->SetBinLabel(ifit,"bindE");
      meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
      sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2));
      ifit++;
    }
    cextra->cd(2); hpull_skdet_energyscale->Draw(); hpull_skdet_energyscale->Fit("gaus");
      fit = hpull_skdet_energyscale->GetFunction("gaus");
    if(fit) {
      meansummary->GetXaxis()->SetBinLabel(ifit,"Escale"); sigmasummary->GetXaxis()->SetBinLabel(ifit,"Escale");
      meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
      sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2));
      ifit++;
    }
    TCanvas *cskdet = new TCanvas(); cskdet->Divide(1,2);
    cskdet->cd(1); hpull_skdet_ccnqe->Draw(); hpull_skdet_ccnqe->Fit("gaus");
    fit = hpull_skdet_ccnqe->GetFunction("gaus");
    if(fit) {
      meansummary->GetXaxis()->SetBinLabel(ifit,"CCNQEeff"); sigmasummary->GetXaxis()->SetBinLabel(ifit,"CCNQEeff");
      meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
      sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2));
      ifit++;
    }
    cskdet->cd(2); hpull_skdet_nc->Draw(); hpull_skdet_nc->Fit("gaus");
    fit = hpull_skdet_nc->GetFunction("gaus");
    if(fit) {
      meansummary->GetXaxis()->SetBinLabel(ifit,"NCeff"); sigmasummary->GetXaxis()->SetBinLabel(ifit,"NCeff");
      meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
      sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2));
      ifit++;
    }

    cbanff->SaveAs(Form("%s_%f_banff_pull.pdf", tag, input_mixing23));
    cxsec->SaveAs(Form("%s_%f_xsec_pull.pdf", tag, input_mixing23));
    cextra->SaveAs(Form("%s_%f_extra_pull.pdf", tag, input_mixing23));
    cskdet->SaveAs(Form("%s_%f_skdet_pull.pdf", tag, input_mixing23));

    //altpull
    ifit++;
    TCanvas *caltbanff = new TCanvas(); caltbanff->Divide(1,2);
    caltbanff->cd(1); haltpull_banff_maqe->Draw(); haltpull_banff_maqe->Fit("gaus");
    fit = haltpull_banff_maqe->GetFunction("gaus");
    if(fit) {
      meansummary->GetXaxis()->SetBinLabel(ifit,"maqe"); sigmasummary->GetXaxis()->SetBinLabel(ifit,"maqe");
      meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
      sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2));
      ifit++;
    }
    caltbanff->cd(2); haltpull_banff_ccqeE0->Draw(); haltpull_banff_ccqeE0->Fit("gaus");
    fit = haltpull_banff_ccqeE0->GetFunction("gaus");
    if(fit) {
      meansummary->GetXaxis()->SetBinLabel(ifit,"ccqeE0"); sigmasummary->GetXaxis()->SetBinLabel(ifit,"ccqeE0");
      meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
      sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2));
      ifit++;
    }
    TCanvas *caltxsec = new TCanvas(); caltxsec->Divide(1,2);
    caltxsec->cd(1); haltpull_xsec_sf->Draw();
    /* haltpull_xsec_sf->Fit("gaus");
    fit = haltpull_xsec_sf->GetFunction("gaus");
    if(fit) {
      meansummary->GetXaxis()->SetBinLabel(ifit,"SF"); sigmasummary->GetXaxis()->SetBinLabel(ifit,"SF");
      meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
      sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2));
      ifit++;
    }
    */
    caltxsec->cd(2); haltpull_xsec_pilessdelta->Draw(); haltpull_xsec_pilessdelta->Fit("gaus");
    if(fit) {
      fit = haltpull_xsec_pilessdelta->GetFunction("gaus");
      meansummary->GetXaxis()->SetBinLabel(ifit,"piless"); sigmasummary->GetXaxis()->SetBinLabel(ifit,"piless");
      meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
      sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2));
      ifit++;
    }
    TCanvas *caltextra = new TCanvas(); caltextra->Divide(1,2);
    caltextra->cd(1); haltpull_xsec_binde->Draw(); haltpull_xsec_binde->Fit("gaus");
    fit = haltpull_xsec_binde->GetFunction("gaus");
    if(fit) {
      meansummary->GetXaxis()->SetBinLabel(ifit,"bindE"); sigmasummary->GetXaxis()->SetBinLabel(ifit,"bindE");
      meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
      sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2));
      ifit++;
    }
    caltextra->cd(2); haltpull_skdet_energyscale->Draw(); haltpull_skdet_energyscale->Fit("gaus");
      fit = haltpull_skdet_energyscale->GetFunction("gaus");
    if(fit) {
      meansummary->GetXaxis()->SetBinLabel(ifit,"Escale"); sigmasummary->GetXaxis()->SetBinLabel(ifit,"Escale");
      meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
      sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2));
      ifit++;
    }
    TCanvas *caltskdet = new TCanvas(); caltskdet->Divide(1,2);
    caltskdet->cd(1); haltpull_skdet_ccnqe->Draw(); haltpull_skdet_ccnqe->Fit("gaus");
    fit = haltpull_skdet_ccnqe->GetFunction("gaus");
    if(fit) {
      meansummary->GetXaxis()->SetBinLabel(ifit,"CCNQEeff"); sigmasummary->GetXaxis()->SetBinLabel(ifit,"CCNQEeff");
      meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
      sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2));
      ifit++;
    }
    caltskdet->cd(2); haltpull_skdet_nc->Draw(); haltpull_skdet_nc->Fit("gaus");
    fit = haltpull_skdet_nc->GetFunction("gaus");
    if(fit) {
      meansummary->GetXaxis()->SetBinLabel(ifit,"NCeff"); sigmasummary->GetXaxis()->SetBinLabel(ifit,"NCeff");
      meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
      sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2));
      ifit++;
    }

    caltbanff->SaveAs(Form("%s_%f_banff_altpull.pdf", tag, input_mixing23));
    caltxsec->SaveAs(Form("%s_%f_xsec_altpull.pdf", tag, input_mixing23));
    caltextra->SaveAs(Form("%s_%f_extra_altpull.pdf", tag, input_mixing23));
    caltskdet->SaveAs(Form("%s_%f_skdet_altpull.pdf", tag, input_mixing23));

    //altpull2
    haltpull_banff_maqe->Draw(); haltpull_banff_maqe->Fit("gaus");
    haltpull_banff_ccqeE0->Draw(); haltpull_banff_ccqeE0->Fit("gaus");
    haltpull_xsec_sf->Draw(); //haltpull_xsec_sf->Fit("gaus");
    haltpull_xsec_pilessdelta->Draw(); haltpull_xsec_pilessdelta->Fit("gaus");
    haltpull_xsec_binde->Draw(); haltpull_xsec_binde->Fit("gaus");
    haltpull_skdet_energyscale->Draw(); haltpull_skdet_energyscale->Fit("gaus");
    haltpull_skdet_ccnqe->Draw(); haltpull_skdet_ccnqe->Fit("gaus");
    haltpull_skdet_nc->Draw(); haltpull_skdet_nc->Fit("gaus");
  }
  if(fit_systematics_mode == 2 && do_gaus_fits) {
    new TCanvas(); hpull_banff_00->Draw(); hpull_banff_00->Fit("gaus"); fit = hpull_banff_00->GetFunction("gaus");
    if(fit) {
      meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
      sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2)); ifit++;
      new TCanvas(); hpull_banff_01->Draw(); hpull_banff_01->Fit("gaus"); fit = hpull_banff_01->GetFunction("gaus");
    }
    if(fit) {
      meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
      sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2)); ifit++;
      new TCanvas(); hpull_banff_02->Draw(); hpull_banff_02->Fit("gaus"); fit = hpull_banff_02->GetFunction("gaus");
    }
    if(fit) {
      meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
      sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2)); ifit++;
      new TCanvas(); hpull_banff_03->Draw(); hpull_banff_03->Fit("gaus"); fit = hpull_banff_03->GetFunction("gaus");
    }
    if(fit) {
      meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
      sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2)); ifit++;
      new TCanvas(); hpull_banff_04->Draw(); hpull_banff_04->Fit("gaus"); fit = hpull_banff_04->GetFunction("gaus");
    }
    if(fit) {
      meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
      sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2)); ifit++;
      new TCanvas(); hpull_banff_05->Draw(); hpull_banff_05->Fit("gaus"); fit = hpull_banff_05->GetFunction("gaus");
    }
    if(fit) {
      meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
      sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2)); ifit++;
      new TCanvas(); hpull_banff_06->Draw(); hpull_banff_06->Fit("gaus"); fit = hpull_banff_06->GetFunction("gaus");
    }
    if(fit) {
      meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
      sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2)); ifit++;
      new TCanvas(); hpull_banff_07->Draw(); hpull_banff_07->Fit("gauss"); fit = hpull_banff_07->GetFunction("gaus");
    }
    if(fit) {
      meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
      sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2)); ifit++;
      new TCanvas(); hpull_banff_08->Draw(); hpull_banff_08->Fit("gaus"); fit = hpull_banff_08->GetFunction("gaus");
    }
    if(fit) {
      meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
      sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2)); ifit++;
      new TCanvas(); hpull_banff_09->Draw(); hpull_banff_09->Fit("gaus"); fit = hpull_banff_09->GetFunction("gaus");
    }
    if(fit) {
      meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
      sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2)); ifit++;
      new TCanvas(); hpull_banff_10->Draw(); hpull_banff_10->Fit("gaus"); fit = hpull_banff_10->GetFunction("gaus");
    }
    if(fit) {
      meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
      sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2)); ifit++;
      new TCanvas(); hpull_banff_11->Draw(); hpull_banff_11->Fit("gaus"); fit = hpull_banff_11->GetFunction("gaus");
    }
    if(fit) {
      meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
      sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2)); ifit++;
      new TCanvas(); hpull_banff_12->Draw(); hpull_banff_12->Fit("gaus"); fit = hpull_banff_12->GetFunction("gaus");
    }
    if(fit) {
      meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
      sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2)); ifit++;
      new TCanvas(); hpull_banff_13->Draw(); hpull_banff_13->Fit("gaus"); fit = hpull_banff_13->GetFunction("gaus");
    }
    if(fit) {
      meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
      sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2)); ifit++;
      new TCanvas(); hpull_banff_14->Draw(); hpull_banff_14->Fit("gaus"); fit = hpull_banff_14->GetFunction("gaus");
    }
    if(fit) {
      meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
      sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2)); ifit++;
      new TCanvas(); hpull_banff_15->Draw(); hpull_banff_15->Fit("gaus"); fit = hpull_banff_15->GetFunction("gaus");
    }
    if(fit) {
      meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
      sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2)); ifit++;
      new TCanvas(); hpull_banff_mares->Draw(); hpull_banff_mares->Fit("gaus"); fit = hpull_banff_mares->GetFunction("gaus");
    }
    if(fit) {
      meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
      sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2)); ifit++;
      new TCanvas(); hpull_banff_ccqeE1->Draw(); hpull_banff_ccqeE1->Fit("gaus"); fit = hpull_banff_ccqeE1->GetFunction("gaus");
    }
    if(fit) {
      meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
      sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2)); ifit++;
      new TCanvas(); hpull_banff_ccqeE2->Draw(); hpull_banff_ccqeE2->Fit("gaus"); fit = hpull_banff_ccqeE2->GetFunction("gaus");
    }
    if(fit) {
      meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
      sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2)); ifit++;
      new TCanvas(); hpull_banff_cc1piE0->Draw(); hpull_banff_cc1piE0->Fit("gaus"); fit = hpull_banff_cc1piE0->GetFunction("gaus");
    }
    if(fit) {
      meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
      sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2)); ifit++;
      new TCanvas(); hpull_banff_cc1piE1->Draw(); hpull_banff_cc1piE1->Fit("gaus"); fit = hpull_banff_cc1piE1->GetFunction("gaus");
    }
    if(fit) {
      meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
      sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2)); ifit++;
    }
    new TCanvas(); hpull_xsec_pf->Draw(); hpull_xsec_pf->Fit("gaus"); fit = hpull_xsec_pf->GetFunction("gaus");
    if(fit) {
      meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
      sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2)); ifit++;
    }
    new TCanvas(); hpull_xsec_wshape->Draw(); hpull_xsec_wshape->Fit("gaus"); fit = hpull_xsec_wshape->GetFunction("gaus");
    if(fit) {
      meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
      sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2)); ifit++;
    }
    new TCanvas(); hpull_xsec_ccothshape->Draw(); hpull_xsec_ccothshape->Fit("gaus"); fit = hpull_xsec_ccothshape->GetFunction("gaus");
    if(fit) {
      meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
      sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2)); ifit++;
    }
    //new TCanvas(); hpull_xsec_cc1pishape->Draw(); hpull_xsec_cc1pishape->Fit("gaus"); fit = hpull_xsec_cc1pishape->GetFunction("gaus");
    //if(fit) {
    //meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
    //sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2)); ifit++;
    //}
    new TCanvas(); hpull_xsec_normcccoh->Draw(); hpull_xsec_normcccoh->Fit("gaus"); fit = hpull_xsec_normcccoh->GetFunction("gaus");
    if(fit) {
      meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
      sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2)); ifit++;
    }
    new TCanvas(); hpull_xsec_normnc1piC->Draw(); hpull_xsec_normnc1piC->Fit("gaus"); fit = hpull_xsec_normnc1piC->GetFunction("gaus");
    if(fit) {
      meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
      sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2)); ifit++;
    }
    //new TCanvas(); hpull_xsec_normnccoh->Draw(); hpull_xsec_normnccoh->Fit("gaus"); fit = hpull_xsec_normnccoh->GetFunction("gaus");
    //if(fit) {
    //meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
    //sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2)); ifit++;
    //}
    new TCanvas(); hpull_xsec_normncoth->Draw(); hpull_xsec_normncoth->Fit("gaus"); fit = hpull_xsec_normncoth->GetFunction("gaus");
    if(fit) {
      meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
      sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2)); ifit++;
    }
    new TCanvas(); hpull_xsec_normnuetonumu->Draw(); hpull_xsec_normnuetonumu->Fit("gaus"); fit = hpull_xsec_normnuetonumu->GetFunction("gaus");
    if(fit) {
      meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
      sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2)); ifit++;
    }
    new TCanvas(); hpull_xsec_normantinutonu->Draw(); hpull_xsec_normantinutonu->Fit("gaus"); fit = hpull_xsec_normantinutonu->GetFunction("gaus");
    if(fit) {
      meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
      sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2)); ifit++;
    }
    new TCanvas(); hpull_skdet_ccqerc0->Draw(); hpull_skdet_ccqerc0->Fit("gaus"); fit = hpull_skdet_ccqerc0->GetFunction("gaus");
    if(fit) {
      meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
      sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2)); ifit++;
    }
    new TCanvas(); hpull_skdet_ccqerc1->Draw(); hpull_skdet_ccqerc1->Fit("gaus"); fit = hpull_skdet_ccqerc1->GetFunction("gaus");
    if(fit) {
      meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
      sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2)); ifit++;
    }
    new TCanvas(); hpull_skdet_ccqerc2->Draw(); hpull_skdet_ccqerc2->Fit("gaus"); fit = hpull_skdet_ccqerc2->GetFunction("gaus");
    if(fit) {
      meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
      sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2)); ifit++;
    }
    new TCanvas(); hpull_skdet_ccqeoth->Draw(); hpull_skdet_ccqeoth->Fit("gaus"); fit = hpull_skdet_ccqeoth->GetFunction("gaus");
    if(fit) {
      meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
      sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2)); ifit++;
    }
    new TCanvas(); hpull_skdet_nuecc->Draw(); hpull_skdet_nuecc->Fit("gaus"); fit = hpull_skdet_nuecc->GetFunction("gaus");
    if(fit) {
      meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
      sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2)); ifit++;
    }
    new TCanvas(); hpull_fsi_000->Draw(); hpull_fsi_000->Fit("gaus"); fit = hpull_fsi_000->GetFunction("gaus");
    if(fit) {
      meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
      sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2)); ifit++;
    }
    new TCanvas(); hpull_fsi_001->Draw(); hpull_fsi_001->Fit("gaus"); fit = hpull_fsi_001->GetFunction("gaus");
    if(fit) {
      meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
      sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2)); ifit++;
    }
    new TCanvas(); hpull_fsi_002->Draw(); hpull_fsi_002->Fit("gaus"); fit = hpull_fsi_002->GetFunction("gaus");
    if(fit) {
      meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
      sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2)); ifit++;
    }
    new TCanvas(); hpull_fsi_003->Draw(); hpull_fsi_003->Fit("gaus"); fit = hpull_fsi_003->GetFunction("gaus");
    if(fit) {
      meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
      sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2)); ifit++;
    }
    new TCanvas(); hpull_fsi_004->Draw(); hpull_fsi_004->Fit("gaus"); fit = hpull_fsi_004->GetFunction("gaus");
    if(fit) {
      meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
      sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2)); ifit++;
    }
    new TCanvas(); hpull_fsi_005->Draw(); hpull_fsi_005->Fit("gaus"); fit = hpull_fsi_005->GetFunction("gaus");
    if(fit) {
      meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
      sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2)); ifit++;
    }
  }

  if(fit_systematics_mode > 0) {
    TCanvas * csummary = new TCanvas(); csummary->Divide(1,2);
    //gStyle->SetOptStat(0);
    if(fit_systematics_mode == 1) {
      meansummary->GetXaxis()->SetRangeUser(0,20);
      sigmasummary->GetXaxis()->SetRangeUser(0,20);
    }
    csummary->cd(1); meansummary->Draw();
    csummary->cd(2); sigmasummary->Draw();
    csummary->SaveAs(Form("%s_%f_summary_pull.pdf", tag, input_mixing23));
    //gStyle->SetOptStat(1);
  }


  //some more plots (better versions, for the tech note)

  //for some reason, need to set the axis titles
  h_mixing23->GetXaxis()->SetTitle("sin^{2}(2#theta_{23})");
  h_dm23sq->GetXaxis()->SetTitle("|#Delta m^{2}_{32}|");
  hpull_mixing23->GetXaxis()->SetTitle("pull_{1}");
  hpull_dm23sq->GetXaxis()->SetTitle("pull_{1}");
  if(fit_systematics_mode == 1) {
    h_banff_maqe->GetXaxis()->SetTitle("#sigma");
    h_banff_ccqeE0->GetXaxis()->SetTitle("#sigma");
    h_xsec_sf->GetXaxis()->SetTitle("#sigma");
    h_xsec_pilessdelta->GetXaxis()->SetTitle("#sigma");
    h_xsec_binde->GetXaxis()->SetTitle("#sigma");
    h_skdet_energyscale->GetXaxis()->SetTitle("#sigma");
    h_skdet_ccnqe->GetXaxis()->SetTitle("#sigma");
    h_skdet_nc->GetXaxis()->SetTitle("#sigma");

    haltpull_banff_maqe->GetXaxis()->SetTitle("pull_{1}");
    haltpull_banff_ccqeE0->GetXaxis()->SetTitle("pull_{1}");
    haltpull_xsec_sf->GetXaxis()->SetTitle("pull_{1}");
    haltpull_xsec_pilessdelta->GetXaxis()->SetTitle("pull_{1}");
    haltpull_xsec_binde->GetXaxis()->SetTitle("pull_{1}");
    haltpull_skdet_energyscale->GetXaxis()->SetTitle("pull_{1}");
    haltpull_skdet_ccnqe->GetXaxis()->SetTitle("pull_{1}");
    haltpull_skdet_nc->GetXaxis()->SetTitle("pull_{1}");

    haltpull2_banff_maqe->GetXaxis()->SetTitle("pull_{3}");
    haltpull2_banff_ccqeE0->GetXaxis()->SetTitle("pull_{3}");
    haltpull2_xsec_sf->GetXaxis()->SetTitle("pull_{3}");
    haltpull2_xsec_pilessdelta->GetXaxis()->SetTitle("pull_{3}");
    haltpull2_xsec_binde->GetXaxis()->SetTitle("pull_{3}");
    haltpull2_skdet_energyscale->GetXaxis()->SetTitle("pull_{3}");
    haltpull2_skdet_ccnqe->GetXaxis()->SetTitle("pull_{3}");
    haltpull2_skdet_nc->GetXaxis()->SetTitle("pull_{3}");

    hpull_banff_maqe->GetXaxis()->SetTitle("pull_{2}");
    hpull_banff_ccqeE0->GetXaxis()->SetTitle("pull_{2}");
    hpull_xsec_sf->GetXaxis()->SetTitle("pull_{2}");
    hpull_xsec_pilessdelta->GetXaxis()->SetTitle("pull_{2}");
    hpull_xsec_binde->GetXaxis()->SetTitle("pull_{2}");
    hpull_skdet_energyscale->GetXaxis()->SetTitle("pull_{2}");
    hpull_skdet_ccnqe->GetXaxis()->SetTitle("pull_{2}");
    hpull_skdet_nc->GetXaxis()->SetTitle("pull_{2}");

    cout << "hahahahha" << endl;
  }//systmode 1
  else if(fit_systematics_mode == 2) {
    h_banff_00->GetXaxis()->SetTitle("#sigma");
    h_banff_01->GetXaxis()->SetTitle("#sigma");
    h_banff_02->GetXaxis()->SetTitle("#sigma");
    h_banff_03->GetXaxis()->SetTitle("#sigma");
    h_banff_04->GetXaxis()->SetTitle("#sigma");
    h_banff_05->GetXaxis()->SetTitle("#sigma");
    h_banff_06->GetXaxis()->SetTitle("#sigma");
    h_banff_07->GetXaxis()->SetTitle("#sigma");
    h_banff_08->GetXaxis()->SetTitle("#sigma");
    h_banff_09->GetXaxis()->SetTitle("#sigma");
    h_banff_10->GetXaxis()->SetTitle("#sigma");
    h_banff_11->GetXaxis()->SetTitle("#sigma");
    h_banff_12->GetXaxis()->SetTitle("#sigma");
    h_banff_13->GetXaxis()->SetTitle("#sigma");
    h_banff_14->GetXaxis()->SetTitle("#sigma");
    h_banff_15->GetXaxis()->SetTitle("#sigma");
    h_banff_mares->GetXaxis()->SetTitle("#sigma");
    h_banff_ccqeE1->GetXaxis()->SetTitle("#sigma");
    h_banff_ccqeE2->GetXaxis()->SetTitle("#sigma");
    h_banff_cc1piE0->GetXaxis()->SetTitle("#sigma");
    h_banff_cc1piE1->GetXaxis()->SetTitle("#sigma");
    h_xsec_pf->GetXaxis()->SetTitle("#sigma");
    h_xsec_wshape->GetXaxis()->SetTitle("#sigma");
    h_xsec_ccothshape->GetXaxis()->SetTitle("#sigma");
    h_xsec_normcccoh->GetXaxis()->SetTitle("#sigma");
    h_xsec_normnc1piC->GetXaxis()->SetTitle("#sigma");
    h_xsec_normncoth->GetXaxis()->SetTitle("#sigma");
    h_xsec_normnuetonumu->GetXaxis()->SetTitle("#sigma");
    h_xsec_normantinutonu->GetXaxis()->SetTitle("#sigma");
    h_skdet_ccqerc0->GetXaxis()->SetTitle("#sigma");
    h_skdet_ccqerc1->GetXaxis()->SetTitle("#sigma");
    h_skdet_ccqerc2->GetXaxis()->SetTitle("#sigma");
    h_skdet_ccqeoth->GetXaxis()->SetTitle("#sigma");
    h_skdet_nuecc->GetXaxis()->SetTitle("#sigma");
    h_fsi_000->GetXaxis()->SetTitle("#sigma");
    h_fsi_001->GetXaxis()->SetTitle("#sigma");
    h_fsi_002->GetXaxis()->SetTitle("#sigma");
    h_fsi_003->GetXaxis()->SetTitle("#sigma");
    h_fsi_004->GetXaxis()->SetTitle("#sigma");
    h_fsi_005->GetXaxis()->SetTitle("#sigma");

    hpull_banff_00->GetXaxis()->SetTitle("pull_{2}");
    hpull_banff_01->GetXaxis()->SetTitle("pull_{2}");
    hpull_banff_02->GetXaxis()->SetTitle("pull_{2}");
    hpull_banff_03->GetXaxis()->SetTitle("pull_{2}");
    hpull_banff_04->GetXaxis()->SetTitle("pull_{2}");
    hpull_banff_05->GetXaxis()->SetTitle("pull_{2}");
    hpull_banff_06->GetXaxis()->SetTitle("pull_{2}");
    hpull_banff_07->GetXaxis()->SetTitle("pull_{2}");
    hpull_banff_08->GetXaxis()->SetTitle("pull_{2}");
    hpull_banff_09->GetXaxis()->SetTitle("pull_{2}");
    hpull_banff_10->GetXaxis()->SetTitle("pull_{2}");
    hpull_banff_11->GetXaxis()->SetTitle("pull_{2}");
    hpull_banff_12->GetXaxis()->SetTitle("pull_{2}");
    hpull_banff_13->GetXaxis()->SetTitle("pull_{2}");
    hpull_banff_14->GetXaxis()->SetTitle("pull_{2}");
    hpull_banff_15->GetXaxis()->SetTitle("pull_{2}");
    hpull_banff_mares->GetXaxis()->SetTitle("pull_{2}");
    hpull_banff_ccqeE1->GetXaxis()->SetTitle("pull_{2}");
    hpull_banff_ccqeE2->GetXaxis()->SetTitle("pull_{2}");
    hpull_banff_cc1piE0->GetXaxis()->SetTitle("pull_{2}");
    hpull_banff_cc1piE1->GetXaxis()->SetTitle("pull_{2}");
    hpull_xsec_pf->GetXaxis()->SetTitle("pull_{2}");
    hpull_xsec_wshape->GetXaxis()->SetTitle("pull_{2}");
    hpull_xsec_ccothshape->GetXaxis()->SetTitle("pull_{2}");
    hpull_xsec_normcccoh->GetXaxis()->SetTitle("pull_{2}");
    hpull_xsec_normnc1piC->GetXaxis()->SetTitle("pull_{2}");
    hpull_xsec_normncoth->GetXaxis()->SetTitle("pull_{2}");
    hpull_xsec_normnuetonumu->GetXaxis()->SetTitle("pull_{2}");
    hpull_xsec_normantinutonu->GetXaxis()->SetTitle("pull_{2}");
    hpull_skdet_ccqerc0->GetXaxis()->SetTitle("pull_{2}");
    hpull_skdet_ccqerc1->GetXaxis()->SetTitle("pull_{2}");
    hpull_skdet_ccqerc2->GetXaxis()->SetTitle("pull_{2}");
    hpull_skdet_ccqeoth->GetXaxis()->SetTitle("pull_{2}");
    hpull_skdet_nuecc->GetXaxis()->SetTitle("pull_{2}");
    hpull_fsi_000->GetXaxis()->SetTitle("pull_{2}");
    hpull_fsi_001->GetXaxis()->SetTitle("pull_{2}");
    hpull_fsi_002->GetXaxis()->SetTitle("pull_{2}");
    hpull_fsi_003->GetXaxis()->SetTitle("pull_{2}");
    hpull_fsi_004->GetXaxis()->SetTitle("pull_{2}");
    hpull_fsi_005->GetXaxis()->SetTitle("pull_{2}");

    haltpull_banff_00->GetXaxis()->SetTitle("pull_{1}");
    haltpull_banff_01->GetXaxis()->SetTitle("pull_{1}");
    haltpull_banff_02->GetXaxis()->SetTitle("pull_{1}");
    haltpull_banff_03->GetXaxis()->SetTitle("pull_{1}");
    haltpull_banff_04->GetXaxis()->SetTitle("pull_{1}");
    haltpull_banff_05->GetXaxis()->SetTitle("pull_{1}");
    haltpull_banff_06->GetXaxis()->SetTitle("pull_{1}");
    haltpull_banff_07->GetXaxis()->SetTitle("pull_{1}");
    haltpull_banff_08->GetXaxis()->SetTitle("pull_{1}");
    haltpull_banff_09->GetXaxis()->SetTitle("pull_{1}");
    haltpull_banff_10->GetXaxis()->SetTitle("pull_{1}");
    haltpull_banff_11->GetXaxis()->SetTitle("pull_{1}");
    haltpull_banff_12->GetXaxis()->SetTitle("pull_{1}");
    haltpull_banff_13->GetXaxis()->SetTitle("pull_{1}");
    haltpull_banff_14->GetXaxis()->SetTitle("pull_{1}");
    haltpull_banff_15->GetXaxis()->SetTitle("pull_{1}");
    haltpull_banff_mares->GetXaxis()->SetTitle("pull_{1}");
    haltpull_banff_ccqeE1->GetXaxis()->SetTitle("pull_{1}");
    haltpull_banff_ccqeE2->GetXaxis()->SetTitle("pull_{1}");
    haltpull_banff_cc1piE0->GetXaxis()->SetTitle("pull_{1}");
    haltpull_banff_cc1piE1->GetXaxis()->SetTitle("pull_{1}");
    haltpull_xsec_pf->GetXaxis()->SetTitle("pull_{1}");
    haltpull_xsec_wshape->GetXaxis()->SetTitle("pull_{1}");
    haltpull_xsec_ccothshape->GetXaxis()->SetTitle("pull_{1}");
    haltpull_xsec_normcccoh->GetXaxis()->SetTitle("pull_{1}");
    haltpull_xsec_normnc1piC->GetXaxis()->SetTitle("pull_{1}");
    haltpull_xsec_normncoth->GetXaxis()->SetTitle("pull_{1}");
    haltpull_xsec_normnuetonumu->GetXaxis()->SetTitle("pull_{1}");
    haltpull_xsec_normantinutonu->GetXaxis()->SetTitle("pull_{1}");
    haltpull_skdet_ccqerc0->GetXaxis()->SetTitle("pull_{1}");
    haltpull_skdet_ccqerc1->GetXaxis()->SetTitle("pull_{1}");
    haltpull_skdet_ccqerc2->GetXaxis()->SetTitle("pull_{1}");
    haltpull_skdet_ccqeoth->GetXaxis()->SetTitle("pull_{1}");
    haltpull_skdet_nuecc->GetXaxis()->SetTitle("pull_{1}");
    haltpull_fsi_000->GetXaxis()->SetTitle("pull_{1}");
    haltpull_fsi_001->GetXaxis()->SetTitle("pull_{1}");
    haltpull_fsi_002->GetXaxis()->SetTitle("pull_{1}");
    haltpull_fsi_003->GetXaxis()->SetTitle("pull_{1}");
    haltpull_fsi_004->GetXaxis()->SetTitle("pull_{1}");
    haltpull_fsi_005->GetXaxis()->SetTitle("pull_{1}");

    haltpull2_banff_00->GetXaxis()->SetTitle("pull_{3}");
    haltpull2_banff_01->GetXaxis()->SetTitle("pull_{3}");
    haltpull2_banff_02->GetXaxis()->SetTitle("pull_{3}");
    haltpull2_banff_03->GetXaxis()->SetTitle("pull_{3}");
    haltpull2_banff_04->GetXaxis()->SetTitle("pull_{3}");
    haltpull2_banff_05->GetXaxis()->SetTitle("pull_{3}");
    haltpull2_banff_06->GetXaxis()->SetTitle("pull_{3}");
    haltpull2_banff_07->GetXaxis()->SetTitle("pull_{3}");
    haltpull2_banff_08->GetXaxis()->SetTitle("pull_{3}");
    haltpull2_banff_09->GetXaxis()->SetTitle("pull_{3}");
    haltpull2_banff_10->GetXaxis()->SetTitle("pull_{3}");
    haltpull2_banff_11->GetXaxis()->SetTitle("pull_{3}");
    haltpull2_banff_12->GetXaxis()->SetTitle("pull_{3}");
    haltpull2_banff_13->GetXaxis()->SetTitle("pull_{3}");
    haltpull2_banff_14->GetXaxis()->SetTitle("pull_{3}");
    haltpull2_banff_15->GetXaxis()->SetTitle("pull_{3}");
    haltpull2_banff_mares->GetXaxis()->SetTitle("pull_{3}");
    haltpull2_banff_ccqeE1->GetXaxis()->SetTitle("pull_{3}");
    haltpull2_banff_ccqeE2->GetXaxis()->SetTitle("pull_{3}");
    haltpull2_banff_cc1piE0->GetXaxis()->SetTitle("pull_{3}");
    haltpull2_banff_cc1piE1->GetXaxis()->SetTitle("pull_{3}");
    haltpull2_xsec_pf->GetXaxis()->SetTitle("pull_{3}");
    haltpull2_xsec_wshape->GetXaxis()->SetTitle("pull_{3}");
    haltpull2_xsec_ccothshape->GetXaxis()->SetTitle("pull_{3}");
    haltpull2_xsec_normcccoh->GetXaxis()->SetTitle("pull_{3}");
    haltpull2_xsec_normnc1piC->GetXaxis()->SetTitle("pull_{3}");
    haltpull2_xsec_normncoth->GetXaxis()->SetTitle("pull_{3}");
    haltpull2_xsec_normnuetonumu->GetXaxis()->SetTitle("pull_{3}");
    haltpull2_xsec_normantinutonu->GetXaxis()->SetTitle("pull_{3}");
    haltpull2_skdet_ccqerc0->GetXaxis()->SetTitle("pull_{3}");
    haltpull2_skdet_ccqerc1->GetXaxis()->SetTitle("pull_{3}");
    haltpull2_skdet_ccqerc2->GetXaxis()->SetTitle("pull_{3}");
    haltpull2_skdet_ccqeoth->GetXaxis()->SetTitle("pull_{3}");
    haltpull2_skdet_nuecc->GetXaxis()->SetTitle("pull_{3}");
    haltpull2_fsi_000->GetXaxis()->SetTitle("pull_{3}");
    haltpull2_fsi_001->GetXaxis()->SetTitle("pull_{3}");
    haltpull2_fsi_002->GetXaxis()->SetTitle("pull_{3}");
    haltpull2_fsi_003->GetXaxis()->SetTitle("pull_{3}");
    haltpull2_fsi_004->GetXaxis()->SetTitle("pull_{3}");
    haltpull2_fsi_005->GetXaxis()->SetTitle("pull_{3}");
  }//systmode 2

  TCanvas * c = new TCanvas();
  h_mixing23->Draw(); 
  c->SaveAs(Form("%s_%f_mixing23.eps", tag, input_mixing23));
  c->SaveAs(Form("%s_%f_mixing23.pdf", tag, input_mixing23));
  h_dm23sq->Draw();
  c->SaveAs(Form("%s_%f_dm23sq.eps", tag, input_mixing23));
  c->SaveAs(Form("%s_%f_dm23sq.pdf", tag, input_mixing23));

  hpull_mixing23->Draw();
  c->SaveAs(Form("%s_%f_mixing23_pull.eps", tag, input_mixing23));
  c->SaveAs(Form("%s_%f_mixing23_pull.pdf", tag, input_mixing23));
  hpull_dm23sq->Draw();
  c->SaveAs(Form("%s_%f_dm23sq_pull.eps", tag, input_mixing23));
  c->SaveAs(Form("%s_%f_dm23sq_pull.pdf", tag, input_mixing23));

  herr_mixing23->Draw();
  c->SaveAs(Form("%s_%f_mixing23_err.eps", tag, input_mixing23));
  c->SaveAs(Form("%s_%f_mixing23_err.pdf", tag, input_mixing23));
  herr_dm23sq->Draw();
  c->SaveAs(Form("%s_%f_dm23sq_err.eps", tag, input_mixing23));
  c->SaveAs(Form("%s_%f_dm23sq_err.pdf", tag, input_mixing23));

  if(fit_systematics_mode == 1) {
    TCanvas * csum = new TCanvas();
    sigmasummary->SetLineColor(kRed);
    meansummary->GetXaxis()->SetRangeUser(0,10);
    sigmasummary->GetXaxis()->SetRangeUser(0,10);
    meansummary->GetYaxis()->SetRangeUser(0,5);
    meansummary->Draw();
    sigmasummary->Draw("same");
    csum->SaveAs(Form("%s_%f_summary_altpull.pdf", tag, input_mixing23));
    meansummary->GetXaxis()->SetRangeUser(11,20);
    sigmasummary->GetXaxis()->SetRangeUser(11,20);
    meansummary->GetYaxis()->SetRangeUser(0,5);
    meansummary->Draw();
    sigmasummary->Draw("same");
    csum->SaveAs(Form("%s_%f_summary_pull.pdf", tag, input_mixing23));
  }
  else if(fit_systematics_mode == 2) {
    TCanvas * csum = new TCanvas();
    sigmasummary->SetLineColor(kRed);
    meansummary->GetXaxis()->SetRangeUser(0,50);
    sigmasummary->GetXaxis()->SetRangeUser(49,100);
    meansummary->GetYaxis()->SetRangeUser(0,5);
    meansummary->Draw();
    sigmasummary->Draw("same");
    csum->SaveAs(Form("%s_%f_summary_altpull.pdf", tag, input_mixing23));
    meansummary->GetXaxis()->SetRangeUser(11,20);
    sigmasummary->GetXaxis()->SetRangeUser(11,20);
    meansummary->GetYaxis()->SetRangeUser(0,5);
    meansummary->Draw();
    sigmasummary->Draw("same");
    csum->SaveAs(Form("%s_%f_summary_pull.pdf", tag, input_mixing23));
  }

  if(fit_systematics_mode == 1) {
    TCanvas * cimp = new TCanvas("","",500,700);
    cimp->Divide(2,4);

    cimp->cd(1); h_skdet_energyscale->Draw();
    cimp->cd(2); h_skdet_ccnqe->Draw();
    cimp->cd(3); h_skdet_nc->Draw();
    cimp->cd(4); h_banff_maqe->Draw();
    cimp->cd(5); h_banff_ccqeE0->Draw();
    cimp->cd(6); h_xsec_sf->Draw();
    cimp->cd(7); h_xsec_binde->Draw();
    cimp->cd(8); h_xsec_pilessdelta->Draw();
    cimp->SaveAs(Form("%s_%f_important.eps", tag, input_mixing23));
    cimp->SaveAs(Form("%s_%f_important.pdf", tag, input_mixing23));

    cimp->cd(1); hpull_skdet_energyscale->Draw();
    cimp->cd(2); hpull_skdet_ccnqe->Draw();
    cimp->cd(3); hpull_skdet_nc->Draw();
    cimp->cd(4); hpull_banff_maqe->Draw();
    cimp->cd(5); hpull_banff_ccqeE0->Draw();
    cimp->cd(6); hpull_xsec_sf->Draw();
    cimp->cd(7); hpull_xsec_binde->Draw();
    cimp->cd(8); hpull_xsec_pilessdelta->Draw();
    cimp->SaveAs(Form("%s_%f_important_pull.eps", tag, input_mixing23));
    cimp->SaveAs(Form("%s_%f_important_pull.pdf", tag, input_mixing23));

    cimp->cd(1); haltpull_skdet_energyscale->Draw();
    cimp->cd(2); haltpull_skdet_ccnqe->Draw();
    cimp->cd(3); haltpull_skdet_nc->Draw();
    cimp->cd(4); haltpull_banff_maqe->Draw();
    cimp->cd(5); haltpull_banff_ccqeE0->Draw();
    cimp->cd(6); haltpull_xsec_sf->Draw();
    cimp->cd(7); haltpull_xsec_binde->Draw();
    cimp->cd(8); haltpull_xsec_pilessdelta->Draw();
    cimp->SaveAs(Form("%s_%f_important_altpull.eps", tag, input_mixing23));
    cimp->SaveAs(Form("%s_%f_important_altpull.pdf", tag, input_mixing23));

    cimp->cd(1); haltpull2_skdet_energyscale->Draw();
    cimp->cd(2); haltpull2_skdet_ccnqe->Draw();
    cimp->cd(3); haltpull2_skdet_nc->Draw();
    cimp->cd(4); haltpull2_banff_maqe->Draw();
    cimp->cd(5); haltpull2_banff_ccqeE0->Draw();
    cimp->cd(6); haltpull2_xsec_sf->Draw();
    cimp->cd(7); haltpull2_xsec_binde->Draw();
    cimp->cd(8); haltpull2_xsec_pilessdelta->Draw();
    cimp->SaveAs(Form("%s_%f_important_altpull2.eps", tag, input_mixing23));
    cimp->SaveAs(Form("%s_%f_important_altpull2.pdf", tag, input_mixing23));

    cimp->cd(1); herr_skdet_energyscale->Draw();
    cimp->cd(2); herr_skdet_ccnqe->Draw();
    cimp->cd(3); herr_skdet_nc->Draw();
    cimp->cd(4); herr_banff_maqe->Draw();
    cimp->cd(5); herr_banff_ccqeE0->Draw();
    cimp->cd(6); herr_xsec_sf->Draw();
    cimp->cd(7); herr_xsec_binde->Draw();
    cimp->cd(8); herr_xsec_pilessdelta->Draw();
    cimp->SaveAs(Form("%s_%f_important_err.eps", tag, input_mixing23));
    cimp->SaveAs(Form("%s_%f_important_err.pdf", tag, input_mixing23));
  }
  else if(fit_systematics_mode == 2) {
    TCanvas * cban = new TCanvas();
    cban->Divide(4,6);
    TCanvas * coth = new TCanvas();
    coth->Divide(4,6);

    cban->cd(1); h_banff_00->Draw();
    cban->cd(2); h_banff_01->Draw();
    cban->cd(3); h_banff_02->Draw();
    cban->cd(4); h_banff_03->Draw();
    cban->cd(5); h_banff_04->Draw();
    cban->cd(6); h_banff_05->Draw();
    cban->cd(7); h_banff_06->Draw();
    cban->cd(8); h_banff_07->Draw();
    cban->cd(9); h_banff_08->Draw();
    cban->cd(10); h_banff_09->Draw();
    cban->cd(11); h_banff_10->Draw();
    cban->cd(12); h_banff_11->Draw();
    cban->cd(13); h_banff_12->Draw();
    cban->cd(14); h_banff_13->Draw();
    cban->cd(15); h_banff_14->Draw();
    cban->cd(16); h_banff_15->Draw();
    cban->cd(17); h_banff_maqe->Draw();
    cban->cd(18); h_banff_mares->Draw();
    cban->cd(19); h_banff_ccqeE0->Draw();
    cban->cd(20); h_banff_ccqeE1->Draw();
    cban->cd(21); h_banff_ccqeE2->Draw();
    cban->cd(22); h_banff_cc1piE0->Draw();
    cban->cd(23); h_banff_cc1piE1->Draw();
    cban->cd(24); h_skdet_energyscale->Draw();
    cban->SaveAs(Form("%s_%f_banff.eps", tag, input_mixing23));
    cban->SaveAs(Form("%s_%f_banff.pdf", tag, input_mixing23));

    coth->cd(1); h_xsec_pf->Draw();
    coth->cd(2); h_xsec_sf->Draw();
    coth->cd(3); h_xsec_wshape->Draw();
    coth->cd(4); h_xsec_ccothshape->Draw();
    coth->cd(5); h_xsec_pilessdelta->Draw();
    coth->cd(6); h_xsec_binde->Draw();
    coth->cd(7); h_xsec_normcccoh->Draw();
    coth->cd(8); h_xsec_normnc1piC->Draw();
    coth->cd(9); h_xsec_normncoth->Draw();
    coth->cd(10); h_xsec_normnuetonumu->Draw();
    coth->cd(11); h_xsec_normantinutonu->Draw();
    coth->cd(12); h_skdet_ccqerc0->Draw();
    coth->cd(13); h_skdet_ccqerc1->Draw();
    coth->cd(14); h_skdet_ccqerc2->Draw();
    coth->cd(15); h_skdet_ccqeoth->Draw();
    coth->cd(16); h_skdet_ccnqe->Draw();
    coth->cd(17); h_skdet_nc->Draw();
    coth->cd(18); h_skdet_nuecc->Draw();
    coth->cd(19); h_fsi_000->Draw();
    coth->cd(20); h_fsi_001->Draw();
    coth->cd(21); h_fsi_002->Draw();
    coth->cd(22); h_fsi_003->Draw();
    coth->cd(23); h_fsi_004->Draw();
    coth->cd(24); h_fsi_005->Draw();
    coth->SaveAs(Form("%s_%f_other.eps", tag, input_mixing23));
    coth->SaveAs(Form("%s_%f_other.pdf", tag, input_mixing23));


    cban->cd(1); hpull_banff_00->Draw();
    cban->cd(2); hpull_banff_01->Draw();
    cban->cd(3); hpull_banff_02->Draw();
    cban->cd(4); hpull_banff_03->Draw();
    cban->cd(5); hpull_banff_04->Draw();
    cban->cd(6); hpull_banff_05->Draw();
    cban->cd(7); hpull_banff_06->Draw();
    cban->cd(8); hpull_banff_07->Draw();
    cban->cd(9); hpull_banff_08->Draw();
    cban->cd(10); hpull_banff_09->Draw();
    cban->cd(11); hpull_banff_10->Draw();
    cban->cd(12); hpull_banff_11->Draw();
    cban->cd(13); hpull_banff_12->Draw();
    cban->cd(14); hpull_banff_13->Draw();
    cban->cd(15); hpull_banff_14->Draw();
    cban->cd(16); hpull_banff_15->Draw();
    cban->cd(17); hpull_banff_maqe->Draw();
    cban->cd(18); hpull_banff_mares->Draw();
    cban->cd(19); hpull_banff_ccqeE0->Draw();
    cban->cd(20); hpull_banff_ccqeE1->Draw();
    cban->cd(21); hpull_banff_ccqeE2->Draw();
    cban->cd(22); hpull_banff_cc1piE0->Draw();
    cban->cd(23); hpull_banff_cc1piE1->Draw();
    cban->cd(24); hpull_skdet_energyscale->Draw();
    cban->SaveAs(Form("%s_%f_banff_pull.eps", tag, input_mixing23));
    cban->SaveAs(Form("%s_%f_banff_pull.pdf", tag, input_mixing23));

    coth->cd(1); hpull_xsec_pf->Draw();
    coth->cd(2); hpull_xsec_sf->Draw();
    coth->cd(3); hpull_xsec_wshape->Draw();
    coth->cd(4); hpull_xsec_ccothshape->Draw();
    coth->cd(5); hpull_xsec_pilessdelta->Draw();
    coth->cd(6); hpull_xsec_binde->Draw();
    coth->cd(7); hpull_xsec_normcccoh->Draw();
    coth->cd(8); hpull_xsec_normnc1piC->Draw();
    coth->cd(9); hpull_xsec_normncoth->Draw();
    coth->cd(10); hpull_xsec_normnuetonumu->Draw();
    coth->cd(11); hpull_xsec_normantinutonu->Draw();
    coth->cd(12); hpull_skdet_ccqerc0->Draw();
    coth->cd(13); hpull_skdet_ccqerc1->Draw();
    coth->cd(14); hpull_skdet_ccqerc2->Draw();
    coth->cd(15); hpull_skdet_ccqeoth->Draw();
    coth->cd(16); hpull_skdet_ccnqe->Draw();
    coth->cd(17); hpull_skdet_nc->Draw();
    coth->cd(18); hpull_skdet_nuecc->Draw();
    coth->cd(19); hpull_fsi_000->Draw();
    coth->cd(20); hpull_fsi_001->Draw();
    coth->cd(21); hpull_fsi_002->Draw();
    coth->cd(22); hpull_fsi_003->Draw();
    coth->cd(23); hpull_fsi_004->Draw();
    coth->cd(24); hpull_fsi_005->Draw();
    coth->SaveAs(Form("%s_%f_other_pull.eps", tag, input_mixing23));
    coth->SaveAs(Form("%s_%f_other_pull.pdf", tag, input_mixing23));


    cban->cd(1); haltpull_banff_00->Draw();
    cban->cd(2); haltpull_banff_01->Draw();
    cban->cd(3); haltpull_banff_02->Draw();
    cban->cd(4); haltpull_banff_03->Draw();
    cban->cd(5); haltpull_banff_04->Draw();
    cban->cd(6); haltpull_banff_05->Draw();
    cban->cd(7); haltpull_banff_06->Draw();
    cban->cd(8); haltpull_banff_07->Draw();
    cban->cd(9); haltpull_banff_08->Draw();
    cban->cd(10); haltpull_banff_09->Draw();
    cban->cd(11); haltpull_banff_10->Draw();
    cban->cd(12); haltpull_banff_11->Draw();
    cban->cd(13); haltpull_banff_12->Draw();
    cban->cd(14); haltpull_banff_13->Draw();
    cban->cd(15); haltpull_banff_14->Draw();
    cban->cd(16); haltpull_banff_15->Draw();
    cban->cd(17); haltpull_banff_maqe->Draw();
    cban->cd(18); haltpull_banff_mares->Draw();
    cban->cd(19); haltpull_banff_ccqeE0->Draw();
    cban->cd(20); haltpull_banff_ccqeE1->Draw();
    cban->cd(21); haltpull_banff_ccqeE2->Draw();
    cban->cd(22); haltpull_banff_cc1piE0->Draw();
    cban->cd(23); haltpull_banff_cc1piE1->Draw();
    cban->cd(24); haltpull_skdet_energyscale->Draw();
    cban->SaveAs(Form("%s_%f_banff_altpull.eps", tag, input_mixing23));
    cban->SaveAs(Form("%s_%f_banff_altpull.pdf", tag, input_mixing23));

    coth->cd(1); haltpull_xsec_pf->Draw();
    coth->cd(2); haltpull_xsec_sf->Draw();
    coth->cd(3); haltpull_xsec_wshape->Draw();
    coth->cd(4); haltpull_xsec_ccothshape->Draw();
    coth->cd(5); haltpull_xsec_pilessdelta->Draw();
    coth->cd(6); haltpull_xsec_binde->Draw();
    coth->cd(7); haltpull_xsec_normcccoh->Draw();
    coth->cd(8); haltpull_xsec_normnc1piC->Draw();
    coth->cd(9); haltpull_xsec_normncoth->Draw();
    coth->cd(10); haltpull_xsec_normnuetonumu->Draw();
    coth->cd(11); haltpull_xsec_normantinutonu->Draw();
    coth->cd(12); haltpull_skdet_ccqerc0->Draw();
    coth->cd(13); haltpull_skdet_ccqerc1->Draw();
    coth->cd(14); haltpull_skdet_ccqerc2->Draw();
    coth->cd(15); haltpull_skdet_ccqeoth->Draw();
    coth->cd(16); haltpull_skdet_ccnqe->Draw();
    coth->cd(17); haltpull_skdet_nc->Draw();
    coth->cd(18); haltpull_skdet_nuecc->Draw();
    coth->cd(19); haltpull_fsi_000->Draw();
    coth->cd(20); haltpull_fsi_001->Draw();
    coth->cd(21); haltpull_fsi_002->Draw();
    coth->cd(22); haltpull_fsi_003->Draw();
    coth->cd(23); haltpull_fsi_004->Draw();
    coth->cd(24); haltpull_fsi_005->Draw();
    coth->SaveAs(Form("%s_%f_other_altpull.eps", tag, input_mixing23));
    coth->SaveAs(Form("%s_%f_other_altpull.pdf", tag, input_mixing23));

    cban->cd(1); haltpull2_banff_00->Draw();
    cban->cd(2); haltpull2_banff_01->Draw();
    cban->cd(3); haltpull2_banff_02->Draw();
    cban->cd(4); haltpull2_banff_03->Draw();
    cban->cd(5); haltpull2_banff_04->Draw();
    cban->cd(6); haltpull2_banff_05->Draw();
    cban->cd(7); haltpull2_banff_06->Draw();
    cban->cd(8); haltpull2_banff_07->Draw();
    cban->cd(9); haltpull2_banff_08->Draw();
    cban->cd(10); haltpull2_banff_09->Draw();
    cban->cd(11); haltpull2_banff_10->Draw();
    cban->cd(12); haltpull2_banff_11->Draw();
    cban->cd(13); haltpull2_banff_12->Draw();
    cban->cd(14); haltpull2_banff_13->Draw();
    cban->cd(15); haltpull2_banff_14->Draw();
    cban->cd(16); haltpull2_banff_15->Draw();
    cban->cd(17); haltpull2_banff_maqe->Draw();
    cban->cd(18); haltpull2_banff_mares->Draw();
    cban->cd(19); haltpull2_banff_ccqeE0->Draw();
    cban->cd(20); haltpull2_banff_ccqeE1->Draw();
    cban->cd(21); haltpull2_banff_ccqeE2->Draw();
    cban->cd(22); haltpull2_banff_cc1piE0->Draw();
    cban->cd(23); haltpull2_banff_cc1piE1->Draw();
    cban->cd(24); haltpull2_skdet_energyscale->Draw();
    cban->SaveAs(Form("%s_%f_banff_altpull2.eps", tag, input_mixing23));
    cban->SaveAs(Form("%s_%f_banff_altpull2.pdf", tag, input_mixing23));

    coth->cd(1); haltpull2_xsec_pf->Draw();
    coth->cd(2); haltpull2_xsec_sf->Draw();
    coth->cd(3); haltpull2_xsec_wshape->Draw();
    coth->cd(4); haltpull2_xsec_ccothshape->Draw();
    coth->cd(5); haltpull2_xsec_pilessdelta->Draw();
    coth->cd(6); haltpull2_xsec_binde->Draw();
    coth->cd(7); haltpull2_xsec_normcccoh->Draw();
    coth->cd(8); haltpull2_xsec_normnc1piC->Draw();
    coth->cd(9); haltpull2_xsec_normncoth->Draw();
    coth->cd(10); haltpull2_xsec_normnuetonumu->Draw();
    coth->cd(11); haltpull2_xsec_normantinutonu->Draw();
    coth->cd(12); haltpull2_skdet_ccqerc0->Draw();
    coth->cd(13); haltpull2_skdet_ccqerc1->Draw();
    coth->cd(14); haltpull2_skdet_ccqerc2->Draw();
    coth->cd(15); haltpull2_skdet_ccqeoth->Draw();
    coth->cd(16); haltpull2_skdet_ccnqe->Draw();
    coth->cd(17); haltpull2_skdet_nc->Draw();
    coth->cd(18); haltpull2_skdet_nuecc->Draw();
    coth->cd(19); haltpull2_fsi_000->Draw();
    coth->cd(20); haltpull2_fsi_001->Draw();
    coth->cd(21); haltpull2_fsi_002->Draw();
    coth->cd(22); haltpull2_fsi_003->Draw();
    coth->cd(23); haltpull2_fsi_004->Draw();
    coth->cd(24); haltpull2_fsi_005->Draw();
    coth->SaveAs(Form("%s_%f_other_altpull2.eps", tag, input_mixing23));
    coth->SaveAs(Form("%s_%f_other_altpull2.pdf", tag, input_mixing23));

    cban->cd(1); herr_banff_00->Draw();
    cban->cd(2); herr_banff_01->Draw();
    cban->cd(3); herr_banff_02->Draw();
    cban->cd(4); herr_banff_03->Draw();
    cban->cd(5); herr_banff_04->Draw();
    cban->cd(6); herr_banff_05->Draw();
    cban->cd(7); herr_banff_06->Draw();
    cban->cd(8); herr_banff_07->Draw();
    cban->cd(9); herr_banff_08->Draw();
    cban->cd(10); herr_banff_09->Draw();
    cban->cd(11); herr_banff_10->Draw();
    cban->cd(12); herr_banff_11->Draw();
    cban->cd(13); herr_banff_12->Draw();
    cban->cd(14); herr_banff_13->Draw();
    cban->cd(15); herr_banff_14->Draw();
    cban->cd(16); herr_banff_15->Draw();
    cban->cd(17); herr_banff_maqe->Draw();
    cban->cd(18); herr_banff_mares->Draw();
    cban->cd(19); herr_banff_ccqeE0->Draw();
    cban->cd(20); herr_banff_ccqeE1->Draw();
    cban->cd(21); herr_banff_ccqeE2->Draw();
    cban->cd(22); herr_banff_cc1piE0->Draw();
    cban->cd(23); herr_banff_cc1piE1->Draw();
    cban->cd(24); herr_skdet_energyscale->Draw();
    cban->SaveAs(Form("%s_%f_banff_err.eps", tag, input_mixing23));
    cban->SaveAs(Form("%s_%f_banff_err.pdf", tag, input_mixing23));

    coth->cd(1); herr_xsec_pf->Draw();
    coth->cd(2); herr_xsec_sf->Draw();
    coth->cd(3); herr_xsec_wshape->Draw();
    coth->cd(4); herr_xsec_ccothshape->Draw();
    coth->cd(5); herr_xsec_pilessdelta->Draw();
    coth->cd(6); herr_xsec_binde->Draw();
    coth->cd(7); herr_xsec_normcccoh->Draw();
    coth->cd(8); herr_xsec_normnc1piC->Draw();
    coth->cd(9); herr_xsec_normncoth->Draw();
    coth->cd(10); herr_xsec_normnuetonumu->Draw();
    coth->cd(11); herr_xsec_normantinutonu->Draw();
    coth->cd(12); herr_skdet_ccqerc0->Draw();
    coth->cd(13); herr_skdet_ccqerc1->Draw();
    coth->cd(14); herr_skdet_ccqerc2->Draw();
    coth->cd(15); herr_skdet_ccqeoth->Draw();
    coth->cd(16); herr_skdet_ccnqe->Draw();
    coth->cd(17); herr_skdet_nc->Draw();
    coth->cd(18); herr_skdet_nuecc->Draw();
    coth->cd(19); herr_fsi_000->Draw();
    coth->cd(20); herr_fsi_001->Draw();
    coth->cd(21); herr_fsi_002->Draw();
    coth->cd(22); herr_fsi_003->Draw();
    coth->cd(23); herr_fsi_004->Draw();
    coth->cd(24); herr_fsi_005->Draw();
    coth->SaveAs(Form("%s_%f_other_err.eps", tag, input_mixing23));
    coth->SaveAs(Form("%s_%f_other_err.pdf", tag, input_mixing23));

    gStyle->SetOptStat(0);
    TH1D * hnewsummary = new TH1D("","",48,0-0.5,48-0.5);
    //banff flux
    hnewsummary->SetBinContent(1,  haltpull_banff_00->GetMean());
    hnewsummary->SetBinError  (1,  haltpull_banff_00->GetRMS());
    hnewsummary->SetBinContent(2,  haltpull_banff_01->GetMean());
    hnewsummary->SetBinError  (2,  haltpull_banff_01->GetRMS());
    hnewsummary->SetBinContent(3,  haltpull_banff_02->GetMean());
    hnewsummary->SetBinError  (3,  haltpull_banff_02->GetRMS());
    hnewsummary->SetBinContent(4,  haltpull_banff_03->GetMean());
    hnewsummary->SetBinError  (4,  haltpull_banff_03->GetRMS());
    hnewsummary->SetBinContent(5,  haltpull_banff_04->GetMean());
    hnewsummary->SetBinError  (5,  haltpull_banff_04->GetRMS());
    hnewsummary->SetBinContent(6,  haltpull_banff_05->GetMean());
    hnewsummary->SetBinError  (6,  haltpull_banff_05->GetRMS());
    hnewsummary->SetBinContent(7,  haltpull_banff_06->GetMean());
    hnewsummary->SetBinError  (7,  haltpull_banff_06->GetRMS());
    hnewsummary->SetBinContent(8,  haltpull_banff_07->GetMean());
    hnewsummary->SetBinError  (8,  haltpull_banff_07->GetRMS());
    hnewsummary->SetBinContent(9,  haltpull_banff_08->GetMean());
    hnewsummary->SetBinError  (9,  haltpull_banff_08->GetRMS());
    hnewsummary->SetBinContent(10, haltpull_banff_09->GetMean());
    hnewsummary->SetBinError  (10, haltpull_banff_09->GetRMS());
    hnewsummary->SetBinContent(11, haltpull_banff_10->GetMean());
    hnewsummary->SetBinError  (11, haltpull_banff_10->GetRMS());
    hnewsummary->SetBinContent(12, haltpull_banff_11->GetMean());
    hnewsummary->SetBinError  (12, haltpull_banff_11->GetRMS());
    hnewsummary->SetBinContent(13, haltpull_banff_12->GetMean());
    hnewsummary->SetBinError  (13, haltpull_banff_12->GetRMS());
    hnewsummary->SetBinContent(14, haltpull_banff_13->GetMean());
    hnewsummary->SetBinError  (14, haltpull_banff_13->GetRMS());
    hnewsummary->SetBinContent(15, haltpull_banff_14->GetMean());
    hnewsummary->SetBinError  (15, haltpull_banff_14->GetRMS());
    hnewsummary->SetBinContent(16, haltpull_banff_15->GetMean());
    hnewsummary->SetBinError  (16, haltpull_banff_15->GetRMS());
    //banff xsec
    hnewsummary->SetBinContent(17, haltpull_banff_maqe->GetMean());
    hnewsummary->SetBinError  (17, haltpull_banff_maqe->GetRMS());
    hnewsummary->SetBinContent(18, haltpull_banff_mares->GetMean());
    hnewsummary->SetBinError  (18, haltpull_banff_mares->GetRMS());
    hnewsummary->SetBinContent(19, haltpull_banff_ccqeE0->GetMean());
    hnewsummary->SetBinError  (19, haltpull_banff_ccqeE0->GetRMS());
    hnewsummary->SetBinContent(20, haltpull_banff_ccqeE1->GetMean());
    hnewsummary->SetBinError  (20, haltpull_banff_ccqeE1->GetRMS());
    hnewsummary->SetBinContent(21, haltpull_banff_ccqeE2->GetMean());
    hnewsummary->SetBinError  (21, haltpull_banff_ccqeE2->GetRMS());
    hnewsummary->SetBinContent(22, haltpull_banff_cc1piE0->GetMean());
    hnewsummary->SetBinError  (22, haltpull_banff_cc1piE0->GetRMS());
    hnewsummary->SetBinContent(23, haltpull_banff_cc1piE1->GetMean());
    hnewsummary->SetBinError  (23, haltpull_banff_cc1piE1->GetRMS());
    //sk e scale
    hnewsummary->SetBinContent(24, haltpull_skdet_energyscale->GetMean());
    hnewsummary->SetBinError  (24, haltpull_skdet_energyscale->GetRMS());
    //xsec splines
    hnewsummary->SetBinContent(25, haltpull_xsec_pf->GetMean());
    hnewsummary->SetBinError  (25, haltpull_xsec_pf->GetRMS());
    hnewsummary->SetBinContent(26, haltpull_xsec_sf->GetMean());
    hnewsummary->SetBinError  (26, haltpull_xsec_sf->GetRMS());
    hnewsummary->SetBinContent(27, haltpull_xsec_wshape->GetMean());
    hnewsummary->SetBinError  (27, haltpull_xsec_wshape->GetRMS());
    hnewsummary->SetBinContent(28, haltpull_xsec_ccothshape->GetMean());
    hnewsummary->SetBinError  (28, haltpull_xsec_ccothshape->GetRMS());
    hnewsummary->SetBinContent(29, haltpull_xsec_pilessdelta->GetMean());
    hnewsummary->SetBinError  (29, haltpull_xsec_pilessdelta->GetRMS());
    hnewsummary->SetBinContent(30, haltpull_xsec_binde->GetMean());
    hnewsummary->SetBinError  (30, haltpull_xsec_binde->GetRMS());
    //xsec norm
    hnewsummary->SetBinContent(31, haltpull_xsec_normcccoh->GetMean());
    hnewsummary->SetBinError  (31, haltpull_xsec_normcccoh->GetRMS());
    hnewsummary->SetBinContent(32, haltpull_xsec_normnc1piC->GetMean());
    hnewsummary->SetBinError  (32, haltpull_xsec_normnc1piC->GetRMS());
    hnewsummary->SetBinContent(33, haltpull_xsec_normncoth->GetMean());
    hnewsummary->SetBinError  (33, haltpull_xsec_normncoth->GetRMS());
    hnewsummary->SetBinContent(34, haltpull_xsec_normnuetonumu->GetMean());
    hnewsummary->SetBinError  (34, haltpull_xsec_normnuetonumu->GetRMS());
    hnewsummary->SetBinContent(35, haltpull_xsec_normantinutonu->GetMean());
    hnewsummary->SetBinError  (35, haltpull_xsec_normantinutonu->GetRMS());
    //skdet
    hnewsummary->SetBinContent(36, haltpull_skdet_ccqerc0->GetMean());
    hnewsummary->SetBinError  (36, haltpull_skdet_ccqerc0->GetRMS());
    hnewsummary->SetBinContent(37, haltpull_skdet_ccqerc1->GetMean());
    hnewsummary->SetBinError  (37, haltpull_skdet_ccqerc1->GetRMS());
    hnewsummary->SetBinContent(38, haltpull_skdet_ccqerc2->GetMean());
    hnewsummary->SetBinError  (38, haltpull_skdet_ccqerc2->GetRMS());
    hnewsummary->SetBinContent(39, haltpull_skdet_ccqeoth->GetMean());
    hnewsummary->SetBinError  (39, haltpull_skdet_ccqeoth->GetRMS());
    hnewsummary->SetBinContent(40, haltpull_skdet_ccnqe->GetMean());
    hnewsummary->SetBinError  (40, haltpull_skdet_ccnqe->GetRMS());
    hnewsummary->SetBinContent(41, haltpull_skdet_nc->GetMean());
    hnewsummary->SetBinError  (41, haltpull_skdet_nc->GetRMS());
    hnewsummary->SetBinContent(42, haltpull_skdet_nuecc->GetMean());
    hnewsummary->SetBinError  (42, haltpull_skdet_nuecc->GetRMS());
    //fsi
    hnewsummary->SetBinContent(43, haltpull_fsi_000->GetMean());
    hnewsummary->SetBinError  (43, haltpull_fsi_000->GetRMS());
    hnewsummary->SetBinContent(44, haltpull_fsi_001->GetMean());
    hnewsummary->SetBinError  (44, haltpull_fsi_001->GetRMS());
    hnewsummary->SetBinContent(45, haltpull_fsi_002->GetMean());
    hnewsummary->SetBinError  (45, haltpull_fsi_002->GetRMS());
    hnewsummary->SetBinContent(46, haltpull_fsi_003->GetMean());
    hnewsummary->SetBinError  (46, haltpull_fsi_003->GetRMS());
    hnewsummary->SetBinContent(47, haltpull_fsi_004->GetMean());
    hnewsummary->SetBinError  (47, haltpull_fsi_004->GetRMS());
    hnewsummary->SetBinContent(48, haltpull_fsi_005->GetMean());
    hnewsummary->SetBinError  (48, haltpull_fsi_005->GetRMS());

    TCanvas * cnewsummary = new TCanvas();
    hnewsummary->Draw();
    cnewsummary->SaveAs(Form("%s_%f_newsummary.eps", tag, input_mixing23));
    cnewsummary->SaveAs(Form("%s_%f_newsummary.pdf", tag, input_mixing23));
  }//sysmode==2

}
