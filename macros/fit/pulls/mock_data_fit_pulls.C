#include <TXMLEngine.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TTree.h>
#include <TChain.h>
#include <TMatrixD.h>
#include <TVectorD.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TStyle.h>
#include <TSystem.h>

#include <Utilities.h>
#include <Conventions.h>

#include <iomanip.h>
#include <stdio.h>
#include <math.h>

void mock_data_fit_pulls(
     const char * filenames,                                              /* add these filenames to the tree */
     int          fit_systematics_mode = 1,                               /* what was fit for? 0: osc only. 1: osc+selected nuisance. 2: osc+all nuisance */
     int          err_mode             = 1,                               /* which errors to use for pulls. 1: migrad, 2 hesse, 3 minos*/
     int          fittermode           = 0,                               /* which fitter has been used? 0: minuit. 1: minuit2 */
     const char * config_file          = "analysis_configuration.xml"     /* file to get the pre errors from */
)
{
  using namespace std;
  using namespace t2kana;

  // First get the expected systematic variations

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

  //Set up chain
  TChain ch("nfitnt");
  ch.Add(filenames);

  // Declaration of leaf types
  //fit info
  int iexp;
  int fitstatus, migradfitstatus, hessefitstatus;
  double fittime;
  double bestfit_chisq, true_chisq;
  double nsk_nom, nsk_bf;
  //oscillation parameters
  double input_mixing12, bestfit_mixing12, errminosplus_mixing12, errminosminus_mixing12, errmigrad_mixing12, errhesse_mixing12, err_mixing12;
  double input_mixing13, bestfit_mixing13, errminosplus_mixing13, errminosminus_mixing13, errmigrad_mixing13, errhesse_mixing13, err_mixing13;
  double input_mixing23, bestfit_mixing23, errminosplus_mixing23, errminosminus_mixing23, errmigrad_mixing23, errhesse_mixing23, err_mixing23;
  double input_deltacp, bestfit_deltacp, errminosplus_deltacp, errminosminus_deltacp, errmigrad_deltacp, errhesse_deltacp, err_deltacp;
  double input_dm12sq, bestfit_dm12sq, errminosplus_dm12sq, errminosminus_dm12sq, errmigrad_dm12sq, errhesse_dm12sq, err_dm12sq;
  double input_dm23sq, bestfit_dm23sq, errminosplus_dm23sq, errminosminus_dm23sq, errmigrad_dm23sq, errhesse_dm23sq, err_dm23sq;
  double input_banff_00, bestfit_banff_00, errminosplus_banff_00, errminosminus_banff_00, errmigrad_banff_00, errhesse_banff_00, err_banff_00;
  double input_banff_01, bestfit_banff_01, errminosplus_banff_01, errminosminus_banff_01, errmigrad_banff_01, errhesse_banff_01, err_banff_01;
  double input_banff_02, bestfit_banff_02, errminosplus_banff_02, errminosminus_banff_02, errmigrad_banff_02, errhesse_banff_02, err_banff_02;
  double input_banff_03, bestfit_banff_03, errminosplus_banff_03, errminosminus_banff_03, errmigrad_banff_03, errhesse_banff_03, err_banff_03;
  double input_banff_04, bestfit_banff_04, errminosplus_banff_04, errminosminus_banff_04, errmigrad_banff_04, errhesse_banff_04, err_banff_04;
  double input_banff_05, bestfit_banff_05, errminosplus_banff_05, errminosminus_banff_05, errmigrad_banff_05, errhesse_banff_05, err_banff_05;
  double input_banff_06, bestfit_banff_06, errminosplus_banff_06, errminosminus_banff_06, errmigrad_banff_06, errhesse_banff_06, err_banff_06;
  double input_banff_07, bestfit_banff_07, errminosplus_banff_07, errminosminus_banff_07, errmigrad_banff_07, errhesse_banff_07, err_banff_07;
  double input_banff_08, bestfit_banff_08, errminosplus_banff_08, errminosminus_banff_08, errmigrad_banff_08, errhesse_banff_08, err_banff_08;
  double input_banff_09, bestfit_banff_09, errminosplus_banff_09, errminosminus_banff_09, errmigrad_banff_09, errhesse_banff_09, err_banff_09;
  double input_banff_10, bestfit_banff_10, errminosplus_banff_10, errminosminus_banff_10, errmigrad_banff_10, errhesse_banff_10, err_banff_10;
  double input_banff_11, bestfit_banff_11, errminosplus_banff_11, errminosminus_banff_11, errmigrad_banff_11, errhesse_banff_11, err_banff_11;
  double input_banff_12, bestfit_banff_12, errminosplus_banff_12, errminosminus_banff_12, errmigrad_banff_12, errhesse_banff_12, err_banff_12;
  double input_banff_13, bestfit_banff_13, errminosplus_banff_13, errminosminus_banff_13, errmigrad_banff_13, errhesse_banff_13, err_banff_13;
  double input_banff_14, bestfit_banff_14, errminosplus_banff_14, errminosminus_banff_14, errmigrad_banff_14, errhesse_banff_14, err_banff_14;
  double input_banff_15, bestfit_banff_15, errminosplus_banff_15, errminosminus_banff_15, errmigrad_banff_15, errhesse_banff_15, err_banff_15;
  double input_banff_maqe, bestfit_banff_maqe, errminosplus_banff_maqe, errminosminus_banff_maqe, errmigrad_banff_maqe, errhesse_banff_maqe, err_banff_maqe;
  double input_banff_mares, bestfit_banff_mares, errminosplus_banff_mares, errminosminus_banff_mares, errmigrad_banff_mares, errhesse_banff_mares, err_banff_mares;
  double input_banff_ccqeE0, bestfit_banff_ccqeE0, errminosplus_banff_ccqeE0, errminosminus_banff_ccqeE0, errmigrad_banff_ccqeE0, errhesse_banff_ccqeE0, err_banff_ccqeE0;
  double input_banff_ccqeE1, bestfit_banff_ccqeE1, errminosplus_banff_ccqeE1, errminosminus_banff_ccqeE1, errmigrad_banff_ccqeE1, errhesse_banff_ccqeE1, err_banff_ccqeE1;
  double input_banff_ccqeE2, bestfit_banff_ccqeE2, errminosplus_banff_ccqeE2, errminosminus_banff_ccqeE2, errmigrad_banff_ccqeE2, errhesse_banff_ccqeE2, err_banff_ccqeE2;
  double input_banff_cc1piE0, bestfit_banff_cc1piE0, errminosplus_banff_cc1piE0, errminosminus_banff_cc1piE0, errmigrad_banff_cc1piE0, errhesse_banff_cc1piE0, err_banff_cc1piE0;
  double input_banff_cc1piE1, bestfit_banff_cc1piE1, errminosplus_banff_cc1piE1, errminosminus_banff_cc1piE1, errmigrad_banff_cc1piE1, errhesse_banff_cc1piE1, err_banff_cc1piE1;
  double input_xsec_pf, bestfit_xsec_pf, errminosplus_xsec_pf, errminosminus_xsec_pf, errmigrad_xsec_pf, errhesse_xsec_pf, err_xsec_pf;
  double input_xsec_sf, bestfit_xsec_sf, errminosplus_xsec_sf, errminosminus_xsec_sf, errmigrad_xsec_sf, errhesse_xsec_sf, err_xsec_sf;
  double input_xsec_binde, bestfit_xsec_binde, errminosplus_xsec_binde, errminosminus_xsec_binde, errmigrad_xsec_binde, errhesse_xsec_binde, err_xsec_binde;
  double input_xsec_wshape, bestfit_xsec_wshape, errminosplus_xsec_wshape, errminosminus_xsec_wshape, errmigrad_xsec_wshape, errhesse_xsec_wshape, err_xsec_wshape;
  double input_xsec_ccothshape, bestfit_xsec_ccothshape, errminosplus_xsec_ccothshape, errminosminus_xsec_ccothshape, errmigrad_xsec_ccothshape, errhesse_xsec_ccothshape, err_xsec_ccothshape;
  double input_xsec_cc1pishape, bestfit_xsec_cc1pishape, errminosplus_xsec_cc1pishape, errminosminus_xsec_cc1pishape, errmigrad_xsec_cc1pishape, errhesse_xsec_cc1pishape, err_xsec_cc1pishape;
  double input_xsec_pilessdelta, bestfit_xsec_pilessdelta, errminosplus_xsec_pilessdelta, errminosminus_xsec_pilessdelta, errmigrad_xsec_pilessdelta, errhesse_xsec_pilessdelta, err_xsec_pilessdelta;
  double input_xsec_normcccoh, bestfit_xsec_normcccoh, errminosplus_xsec_normcccoh, errminosminus_xsec_normcccoh, errmigrad_xsec_normcccoh, errhesse_xsec_normcccoh, err_xsec_normcccoh;
  double input_xsec_normnc1piC, bestfit_xsec_normnc1piC, errminosplus_xsec_normnc1piC, errminosminus_xsec_normnc1piC, errmigrad_xsec_normnc1piC, errhesse_xsec_normnc1piC, err_xsec_normnc1piC;
  double input_xsec_normnccoh, bestfit_xsec_normnccoh, errminosplus_xsec_normnccoh, errminosminus_xsec_normnccoh, errmigrad_xsec_normnccoh, errhesse_xsec_normnccoh, err_xsec_normnccoh;
  double input_xsec_normncoth, bestfit_xsec_normncoth, errminosplus_xsec_normncoth, errminosminus_xsec_normncoth, errmigrad_xsec_normncoth, errhesse_xsec_normncoth, err_xsec_normncoth;
  double input_xsec_normnuetonumu, bestfit_xsec_normnuetonumu, errminosplus_xsec_normnuetonumu, errminosminus_xsec_normnuetonumu, errmigrad_xsec_normnuetonumu, errhesse_xsec_normnuetonumu, err_xsec_normnuetonumu;
  double input_skdet_energyscale, bestfit_skdet_energyscale, errminosplus_skdet_energyscale, errminosminus_skdet_energyscale, errmigrad_skdet_energyscale, errhesse_skdet_energyscale, err_skdet_energyscale;
  double input_skdet_ccqerc0, bestfit_skdet_ccqerc0, errminosplus_skdet_ccqerc0, errminosminus_skdet_ccqerc0, errmigrad_skdet_ccqerc0, errhesse_skdet_ccqerc0, err_skdet_ccqerc0;
  double input_skdet_ccqerc1, bestfit_skdet_ccqerc1, errminosplus_skdet_ccqerc1, errminosminus_skdet_ccqerc1, errmigrad_skdet_ccqerc1, errhesse_skdet_ccqerc1, err_skdet_ccqerc1;
  double input_skdet_ccqerc2, bestfit_skdet_ccqerc2, errminosplus_skdet_ccqerc2, errminosminus_skdet_ccqerc2, errmigrad_skdet_ccqerc2, errhesse_skdet_ccqerc2, err_skdet_ccqerc2;
  double input_skdet_ccqeoth, bestfit_skdet_ccqeoth, errminosplus_skdet_ccqeoth, errminosminus_skdet_ccqeoth, errmigrad_skdet_ccqeoth, errhesse_skdet_ccqeoth, err_skdet_ccqeoth;
  double input_skdet_ccnqe, bestfit_skdet_ccnqe, errminosplus_skdet_ccnqe, errminosminus_skdet_ccnqe, errmigrad_skdet_ccnqe, errhesse_skdet_ccnqe, err_skdet_ccnqe;
  double input_skdet_nc, bestfit_skdet_nc, errminosplus_skdet_nc, errminosminus_skdet_nc, errmigrad_skdet_nc, errhesse_skdet_nc, err_skdet_nc;
  double input_skdet_nuecc, bestfit_skdet_nuecc, errminosplus_skdet_nuecc, errminosminus_skdet_nuecc, errmigrad_skdet_nuecc, errhesse_skdet_nuecc, err_skdet_nuecc;
  double input_fsi_000, bestfit_fsi_000, errminosplus_fsi_000, errminosminus_fsi_000, errmigrad_fsi_000, errhesse_fsi_000, err_fsi_000;
  double input_fsi_001, bestfit_fsi_001, errminosplus_fsi_001, errminosminus_fsi_001, errmigrad_fsi_001, errhesse_fsi_001, err_fsi_001;
  double input_fsi_002, bestfit_fsi_002, errminosplus_fsi_002, errminosminus_fsi_002, errmigrad_fsi_002, errhesse_fsi_002, err_fsi_002;
  double input_fsi_003, bestfit_fsi_003, errminosplus_fsi_003, errminosminus_fsi_003, errmigrad_fsi_003, errhesse_fsi_003, err_fsi_003;
  double input_fsi_004, bestfit_fsi_004, errminosplus_fsi_004, errminosminus_fsi_004, errmigrad_fsi_004, errhesse_fsi_004, err_fsi_004;
  double input_fsi_005, bestfit_fsi_005, errminosplus_fsi_005, errminosminus_fsi_005, errmigrad_fsi_005, errhesse_fsi_005, err_fsi_005;

  // Set branch addresses and branch pointers
  //fit info
  ch.GetEntry(0);
  ch.SetBranchAddress("iexp",&iexp);
  ch.SetBranchAddress("fitstatus",&fitstatus);
  ch.SetBranchAddress("migradfitstatus",&migradfitstatus);
  ch.SetBranchAddress("hessefitstatus",&hessefitstatus);
  ch.SetBranchAddress("fittime",&fittime);
  ch.SetBranchAddress("bestfit_chisq",&bestfit_chisq);
  ch.SetBranchAddress("true_chisq",&true_chisq);
  ch.SetBranchAddress("nsk_nom",&nsk_nom);
  ch.SetBranchAddress("nsk_bf",&nsk_bf);
  //oscillation parameters
  ch.SetBranchAddress("input_mixing12",&input_mixing12);
  ch.SetBranchAddress("input_mixing13",&input_mixing13);
  ch.SetBranchAddress("input_mixing23",&input_mixing23);
  ch.SetBranchAddress("input_deltacp",&input_deltacp);
  ch.SetBranchAddress("input_dm12sq",&input_dm12sq);
  ch.SetBranchAddress("input_dm23sq",&input_dm23sq);
  ch.SetBranchAddress("bestfit_mixing12",&bestfit_mixing12);
  ch.SetBranchAddress("bestfit_mixing13",&bestfit_mixing13);
  ch.SetBranchAddress("bestfit_mixing23",&bestfit_mixing23);
  ch.SetBranchAddress("bestfit_deltacp",&bestfit_deltacp);
  ch.SetBranchAddress("bestfit_dm12sq",&bestfit_dm12sq);
  ch.SetBranchAddress("bestfit_dm23sq",&bestfit_dm23sq);
  ch.SetBranchAddress("errminosplus_mixing12",&errminosplus_mixing12);
  ch.SetBranchAddress("errminosplus_mixing13",&errminosplus_mixing13);
  ch.SetBranchAddress("errminosplus_mixing23",&errminosplus_mixing23);
  ch.SetBranchAddress("errminosplus_deltacp",&errminosplus_deltacp);
  ch.SetBranchAddress("errminosplus_dm12sq",&errminosplus_dm12sq);
  ch.SetBranchAddress("errminosplus_dm23sq",&errminosplus_dm23sq);
  ch.SetBranchAddress("errminosminus_mixing12",&errminosminus_mixing12);
  ch.SetBranchAddress("errminosminus_mixing13",&errminosminus_mixing13);
  ch.SetBranchAddress("errminosminus_mixing23",&errminosminus_mixing23);
  ch.SetBranchAddress("errminosminus_deltacp",&errminosminus_deltacp);
  ch.SetBranchAddress("errminosminus_dm12sq",&errminosminus_dm12sq);
  ch.SetBranchAddress("errminosminus_dm23sq",&errminosminus_dm23sq);
  ch.SetBranchAddress("errmigrad_mixing12",&errmigrad_mixing12);
  ch.SetBranchAddress("errmigrad_mixing13",&errmigrad_mixing13);
  ch.SetBranchAddress("errmigrad_mixing23",&errmigrad_mixing23);
  ch.SetBranchAddress("errmigrad_deltacp",&errmigrad_deltacp);
  ch.SetBranchAddress("errmigrad_dm12sq",&errmigrad_dm12sq);
  ch.SetBranchAddress("errmigrad_dm23sq",&errmigrad_dm23sq);
  ch.SetBranchAddress("errhesse_mixing12",&errhesse_mixing12);
  ch.SetBranchAddress("errhesse_mixing13",&errhesse_mixing13);
  ch.SetBranchAddress("errhesse_mixing23",&errhesse_mixing23);
  ch.SetBranchAddress("errhesse_deltacp",&errhesse_deltacp);
  ch.SetBranchAddress("errhesse_dm12sq",&errhesse_dm12sq);
  ch.SetBranchAddress("errhesse_dm23sq",&errhesse_dm23sq);
  if(fit_systematics_mode >= 1) {
    //inputs
    ch.SetBranchAddress("input_xsec_sf",&input_xsec_sf);
    ch.SetBranchAddress("input_xsec_binde",&input_xsec_binde);
    ch.SetBranchAddress("input_xsec_pilessdelta",&input_xsec_pilessdelta);
    ch.SetBranchAddress("input_skdet_ccnqe",&input_skdet_ccnqe);
    ch.SetBranchAddress("input_skdet_nc",&input_skdet_nc);
    ch.SetBranchAddress("input_skdet_energyscale",&input_skdet_energyscale);
    //bestfits
    ch.SetBranchAddress("bestfit_xsec_sf",&bestfit_xsec_sf);
    ch.SetBranchAddress("bestfit_xsec_binde",&bestfit_xsec_binde);
    ch.SetBranchAddress("bestfit_xsec_pilessdelta",&bestfit_xsec_pilessdelta);
    ch.SetBranchAddress("bestfit_skdet_ccnqe",&bestfit_skdet_ccnqe);
    ch.SetBranchAddress("bestfit_skdet_nc",&bestfit_skdet_nc);
    ch.SetBranchAddress("bestfit_skdet_energyscale",&bestfit_skdet_energyscale);
    //errminospluss
    ch.SetBranchAddress("errminosplus_xsec_sf",&errminosplus_xsec_sf);
    ch.SetBranchAddress("errminosplus_xsec_binde",&errminosplus_xsec_binde);
    ch.SetBranchAddress("errminosplus_xsec_pilessdelta",&errminosplus_xsec_pilessdelta);
    ch.SetBranchAddress("errminosplus_skdet_ccnqe",&errminosplus_skdet_ccnqe);
    ch.SetBranchAddress("errminosplus_skdet_nc",&errminosplus_skdet_nc);
    ch.SetBranchAddress("errminosplus_skdet_energyscale",&errminosplus_skdet_energyscale);
    //errminosminuss
    ch.SetBranchAddress("errminosminus_xsec_sf",&errminosminus_xsec_sf);
    ch.SetBranchAddress("errminosminus_xsec_binde",&errminosminus_xsec_binde);
    ch.SetBranchAddress("errminosminus_xsec_pilessdelta",&errminosminus_xsec_pilessdelta);
    ch.SetBranchAddress("errminosminus_skdet_ccnqe",&errminosminus_skdet_ccnqe);
    ch.SetBranchAddress("errminosminus_skdet_nc",&errminosminus_skdet_nc);
    ch.SetBranchAddress("errminosminus_skdet_energyscale",&errminosminus_skdet_energyscale);
    //errmigrads
    ch.SetBranchAddress("errmigrad_xsec_sf",&errmigrad_xsec_sf);
    ch.SetBranchAddress("errmigrad_xsec_binde",&errmigrad_xsec_binde);
    ch.SetBranchAddress("errmigrad_xsec_pilessdelta",&errmigrad_xsec_pilessdelta);
    ch.SetBranchAddress("errmigrad_skdet_ccnqe",&errmigrad_skdet_ccnqe);
    ch.SetBranchAddress("errmigrad_skdet_nc",&errmigrad_skdet_nc);
    ch.SetBranchAddress("errmigrad_skdet_energyscale",&errmigrad_skdet_energyscale);
    //errhesses
    ch.SetBranchAddress("errhesse_xsec_sf",&errhesse_xsec_sf);
    ch.SetBranchAddress("errhesse_xsec_binde",&errhesse_xsec_binde);
    ch.SetBranchAddress("errhesse_xsec_pilessdelta",&errhesse_xsec_pilessdelta);
    ch.SetBranchAddress("errhesse_skdet_ccnqe",&errhesse_skdet_ccnqe);
    ch.SetBranchAddress("errhesse_skdet_nc",&errhesse_skdet_nc);
    ch.SetBranchAddress("errhesse_skdet_energyscale",&errhesse_skdet_energyscale);

    //inputs
    ch.SetBranchAddress("input_banff_00",&input_banff_00);
    ch.SetBranchAddress("input_banff_01",&input_banff_01);
    ch.SetBranchAddress("input_banff_02",&input_banff_02);
    ch.SetBranchAddress("input_banff_03",&input_banff_03);
    ch.SetBranchAddress("input_banff_04",&input_banff_04);
    ch.SetBranchAddress("input_banff_05",&input_banff_05);
    ch.SetBranchAddress("input_banff_06",&input_banff_06);
    ch.SetBranchAddress("input_banff_07",&input_banff_07);
    ch.SetBranchAddress("input_banff_08",&input_banff_08);
    ch.SetBranchAddress("input_banff_09",&input_banff_09);
    ch.SetBranchAddress("input_banff_10",&input_banff_10);
    ch.SetBranchAddress("input_banff_11",&input_banff_11);
    ch.SetBranchAddress("input_banff_12",&input_banff_12);
    ch.SetBranchAddress("input_banff_13",&input_banff_13);
    ch.SetBranchAddress("input_banff_14",&input_banff_14);
    ch.SetBranchAddress("input_banff_15",&input_banff_15);
    ch.SetBranchAddress("input_banff_maqe",&input_banff_maqe);
    ch.SetBranchAddress("input_banff_mares",&input_banff_mares);
    ch.SetBranchAddress("input_banff_ccqeE0",&input_banff_ccqeE0);
    ch.SetBranchAddress("input_banff_ccqeE1",&input_banff_ccqeE1);
    ch.SetBranchAddress("input_banff_ccqeE2",&input_banff_ccqeE2);
    ch.SetBranchAddress("input_banff_cc1piE0",&input_banff_cc1piE0);
    ch.SetBranchAddress("input_banff_cc1piE1",&input_banff_cc1piE1);
    //bestfits
    ch.SetBranchAddress("bestfit_banff_00",&bestfit_banff_00);
    ch.SetBranchAddress("bestfit_banff_01",&bestfit_banff_01);
    ch.SetBranchAddress("bestfit_banff_02",&bestfit_banff_02);
    ch.SetBranchAddress("bestfit_banff_03",&bestfit_banff_03);
    ch.SetBranchAddress("bestfit_banff_04",&bestfit_banff_04);
    ch.SetBranchAddress("bestfit_banff_05",&bestfit_banff_05);
    ch.SetBranchAddress("bestfit_banff_06",&bestfit_banff_06);
    ch.SetBranchAddress("bestfit_banff_07",&bestfit_banff_07);
    ch.SetBranchAddress("bestfit_banff_08",&bestfit_banff_08);
    ch.SetBranchAddress("bestfit_banff_09",&bestfit_banff_09);
    ch.SetBranchAddress("bestfit_banff_10",&bestfit_banff_10);
    ch.SetBranchAddress("bestfit_banff_11",&bestfit_banff_11);
    ch.SetBranchAddress("bestfit_banff_12",&bestfit_banff_12);
    ch.SetBranchAddress("bestfit_banff_13",&bestfit_banff_13);
    ch.SetBranchAddress("bestfit_banff_14",&bestfit_banff_14);
    ch.SetBranchAddress("bestfit_banff_15",&bestfit_banff_15);
    ch.SetBranchAddress("bestfit_banff_maqe",&bestfit_banff_maqe);
    ch.SetBranchAddress("bestfit_banff_mares",&bestfit_banff_mares);
    ch.SetBranchAddress("bestfit_banff_ccqeE0",&bestfit_banff_ccqeE0);
    ch.SetBranchAddress("bestfit_banff_ccqeE1",&bestfit_banff_ccqeE1);
    ch.SetBranchAddress("bestfit_banff_ccqeE2",&bestfit_banff_ccqeE2);
    ch.SetBranchAddress("bestfit_banff_cc1piE0",&bestfit_banff_cc1piE0);
    ch.SetBranchAddress("bestfit_banff_cc1piE1",&bestfit_banff_cc1piE1);
    //errminospluss
    ch.SetBranchAddress("errminosplus_banff_00",&errminosplus_banff_00);
    ch.SetBranchAddress("errminosplus_banff_01",&errminosplus_banff_01);
    ch.SetBranchAddress("errminosplus_banff_02",&errminosplus_banff_02);
    ch.SetBranchAddress("errminosplus_banff_03",&errminosplus_banff_03);
    ch.SetBranchAddress("errminosplus_banff_04",&errminosplus_banff_04);
    ch.SetBranchAddress("errminosplus_banff_05",&errminosplus_banff_05);
    ch.SetBranchAddress("errminosplus_banff_06",&errminosplus_banff_06);
    ch.SetBranchAddress("errminosplus_banff_07",&errminosplus_banff_07);
    ch.SetBranchAddress("errminosplus_banff_08",&errminosplus_banff_08);
    ch.SetBranchAddress("errminosplus_banff_09",&errminosplus_banff_09);
    ch.SetBranchAddress("errminosplus_banff_10",&errminosplus_banff_10);
    ch.SetBranchAddress("errminosplus_banff_11",&errminosplus_banff_11);
    ch.SetBranchAddress("errminosplus_banff_12",&errminosplus_banff_12);
    ch.SetBranchAddress("errminosplus_banff_13",&errminosplus_banff_13);
    ch.SetBranchAddress("errminosplus_banff_14",&errminosplus_banff_14);
    ch.SetBranchAddress("errminosplus_banff_15",&errminosplus_banff_15);
    ch.SetBranchAddress("errminosplus_banff_maqe",&errminosplus_banff_maqe);
    ch.SetBranchAddress("errminosplus_banff_mares",&errminosplus_banff_mares);
    ch.SetBranchAddress("errminosplus_banff_ccqeE0",&errminosplus_banff_ccqeE0);
    ch.SetBranchAddress("errminosplus_banff_ccqeE1",&errminosplus_banff_ccqeE1);
    ch.SetBranchAddress("errminosplus_banff_ccqeE2",&errminosplus_banff_ccqeE2);
    ch.SetBranchAddress("errminosplus_banff_cc1piE0",&errminosplus_banff_cc1piE0);
    ch.SetBranchAddress("errminosplus_banff_cc1piE1",&errminosplus_banff_cc1piE1);
    //errminosminuss
    ch.SetBranchAddress("errminosminus_banff_00",&errminosminus_banff_00);
    ch.SetBranchAddress("errminosminus_banff_01",&errminosminus_banff_01);
    ch.SetBranchAddress("errminosminus_banff_02",&errminosminus_banff_02);
    ch.SetBranchAddress("errminosminus_banff_03",&errminosminus_banff_03);
    ch.SetBranchAddress("errminosminus_banff_04",&errminosminus_banff_04);
    ch.SetBranchAddress("errminosminus_banff_05",&errminosminus_banff_05);
    ch.SetBranchAddress("errminosminus_banff_06",&errminosminus_banff_06);
    ch.SetBranchAddress("errminosminus_banff_07",&errminosminus_banff_07);
    ch.SetBranchAddress("errminosminus_banff_08",&errminosminus_banff_08);
    ch.SetBranchAddress("errminosminus_banff_09",&errminosminus_banff_09);
    ch.SetBranchAddress("errminosminus_banff_10",&errminosminus_banff_10);
    ch.SetBranchAddress("errminosminus_banff_11",&errminosminus_banff_11);
    ch.SetBranchAddress("errminosminus_banff_12",&errminosminus_banff_12);
    ch.SetBranchAddress("errminosminus_banff_13",&errminosminus_banff_13);
    ch.SetBranchAddress("errminosminus_banff_14",&errminosminus_banff_14);
    ch.SetBranchAddress("errminosminus_banff_15",&errminosminus_banff_15);
    ch.SetBranchAddress("errminosminus_banff_maqe",&errminosminus_banff_maqe);
    ch.SetBranchAddress("errminosminus_banff_mares",&errminosminus_banff_mares);
    ch.SetBranchAddress("errminosminus_banff_ccqeE0",&errminosminus_banff_ccqeE0);
    ch.SetBranchAddress("errminosminus_banff_ccqeE1",&errminosminus_banff_ccqeE1);
    ch.SetBranchAddress("errminosminus_banff_ccqeE2",&errminosminus_banff_ccqeE2);
    ch.SetBranchAddress("errminosminus_banff_cc1piE0",&errminosminus_banff_cc1piE0);
    ch.SetBranchAddress("errminosminus_banff_cc1piE1",&errminosminus_banff_cc1piE1);
    //errmigrad
    ch.SetBranchAddress("errmigrad_banff_00",&errmigrad_banff_00);
    ch.SetBranchAddress("errmigrad_banff_01",&errmigrad_banff_01);
    ch.SetBranchAddress("errmigrad_banff_02",&errmigrad_banff_02);
    ch.SetBranchAddress("errmigrad_banff_03",&errmigrad_banff_03);
    ch.SetBranchAddress("errmigrad_banff_04",&errmigrad_banff_04);
    ch.SetBranchAddress("errmigrad_banff_05",&errmigrad_banff_05);
    ch.SetBranchAddress("errmigrad_banff_06",&errmigrad_banff_06);
    ch.SetBranchAddress("errmigrad_banff_07",&errmigrad_banff_07);
    ch.SetBranchAddress("errmigrad_banff_08",&errmigrad_banff_08);
    ch.SetBranchAddress("errmigrad_banff_09",&errmigrad_banff_09);
    ch.SetBranchAddress("errmigrad_banff_10",&errmigrad_banff_10);
    ch.SetBranchAddress("errmigrad_banff_11",&errmigrad_banff_11);
    ch.SetBranchAddress("errmigrad_banff_12",&errmigrad_banff_12);
    ch.SetBranchAddress("errmigrad_banff_13",&errmigrad_banff_13);
    ch.SetBranchAddress("errmigrad_banff_14",&errmigrad_banff_14);
    ch.SetBranchAddress("errmigrad_banff_15",&errmigrad_banff_15);
    ch.SetBranchAddress("errmigrad_banff_maqe",&errmigrad_banff_maqe);
    ch.SetBranchAddress("errmigrad_banff_mares",&errmigrad_banff_mares);
    ch.SetBranchAddress("errmigrad_banff_ccqeE0",&errmigrad_banff_ccqeE0);
    ch.SetBranchAddress("errmigrad_banff_ccqeE1",&errmigrad_banff_ccqeE1);
    ch.SetBranchAddress("errmigrad_banff_ccqeE2",&errmigrad_banff_ccqeE2);
    ch.SetBranchAddress("errmigrad_banff_cc1piE0",&errmigrad_banff_cc1piE0);
    ch.SetBranchAddress("errmigrad_banff_cc1piE1",&errmigrad_banff_cc1piE1);
    //errhesse
    ch.SetBranchAddress("errhesse_banff_00",&errhesse_banff_00);
    ch.SetBranchAddress("errhesse_banff_01",&errhesse_banff_01);
    ch.SetBranchAddress("errhesse_banff_02",&errhesse_banff_02);
    ch.SetBranchAddress("errhesse_banff_03",&errhesse_banff_03);
    ch.SetBranchAddress("errhesse_banff_04",&errhesse_banff_04);
    ch.SetBranchAddress("errhesse_banff_05",&errhesse_banff_05);
    ch.SetBranchAddress("errhesse_banff_06",&errhesse_banff_06);
    ch.SetBranchAddress("errhesse_banff_07",&errhesse_banff_07);
    ch.SetBranchAddress("errhesse_banff_08",&errhesse_banff_08);
    ch.SetBranchAddress("errhesse_banff_09",&errhesse_banff_09);
    ch.SetBranchAddress("errhesse_banff_10",&errhesse_banff_10);
    ch.SetBranchAddress("errhesse_banff_11",&errhesse_banff_11);
    ch.SetBranchAddress("errhesse_banff_12",&errhesse_banff_12);
    ch.SetBranchAddress("errhesse_banff_13",&errhesse_banff_13);
    ch.SetBranchAddress("errhesse_banff_14",&errhesse_banff_14);
    ch.SetBranchAddress("errhesse_banff_15",&errhesse_banff_15);
    ch.SetBranchAddress("errhesse_banff_maqe",&errhesse_banff_maqe);
    ch.SetBranchAddress("errhesse_banff_mares",&errhesse_banff_mares);
    ch.SetBranchAddress("errhesse_banff_ccqeE0",&errhesse_banff_ccqeE0);
    ch.SetBranchAddress("errhesse_banff_ccqeE1",&errhesse_banff_ccqeE1);
    ch.SetBranchAddress("errhesse_banff_ccqeE2",&errhesse_banff_ccqeE2);
    ch.SetBranchAddress("errhesse_banff_cc1piE0",&errhesse_banff_cc1piE0);
    ch.SetBranchAddress("errhesse_banff_cc1piE1",&errhesse_banff_cc1piE1);
  }
  if(fit_systematics_mode == 2) {
    //inputs
    ch.SetBranchAddress("input_xsec_pf",&input_xsec_pf);
    ch.SetBranchAddress("input_xsec_wshape",&input_xsec_wshape);
    ch.SetBranchAddress("input_xsec_ccothshape",&input_xsec_ccothshape);
    ch.SetBranchAddress("input_xsec_cc1pishape",&input_xsec_cc1pishape);
    ch.SetBranchAddress("input_xsec_normcccoh",&input_xsec_normcccoh);
    ch.SetBranchAddress("input_xsec_normnc1piC",&input_xsec_normnc1piC);
    ch.SetBranchAddress("input_xsec_normnccoh",&input_xsec_normnccoh);
    ch.SetBranchAddress("input_xsec_normncoth",&input_xsec_normncoth);
    ch.SetBranchAddress("input_xsec_normnuetonumu",&input_xsec_normnuetonumu);
    ch.SetBranchAddress("input_skdet_ccqerc0",&input_skdet_ccqerc0);
    ch.SetBranchAddress("input_skdet_ccqerc1",&input_skdet_ccqerc1);
    ch.SetBranchAddress("input_skdet_ccqerc2",&input_skdet_ccqerc2);
    ch.SetBranchAddress("input_skdet_ccqeoth",&input_skdet_ccqeoth);
    ch.SetBranchAddress("input_skdet_nuecc",&input_skdet_nuecc);
    ch.SetBranchAddress("input_fsi_000",&input_fsi_000);
    ch.SetBranchAddress("input_fsi_001",&input_fsi_001);
    ch.SetBranchAddress("input_fsi_002",&input_fsi_002);
    ch.SetBranchAddress("input_fsi_003",&input_fsi_003);
    ch.SetBranchAddress("input_fsi_004",&input_fsi_004);
    ch.SetBranchAddress("input_fsi_005",&input_fsi_005);
    //bestfits
    ch.SetBranchAddress("bestfit_xsec_pf",&bestfit_xsec_pf);
    ch.SetBranchAddress("bestfit_xsec_wshape",&bestfit_xsec_wshape);
    ch.SetBranchAddress("bestfit_xsec_ccothshape",&bestfit_xsec_ccothshape);
    ch.SetBranchAddress("bestfit_xsec_cc1pishape",&bestfit_xsec_cc1pishape);
    ch.SetBranchAddress("bestfit_xsec_normcccoh",&bestfit_xsec_normcccoh);
    ch.SetBranchAddress("bestfit_xsec_normnc1piC",&bestfit_xsec_normnc1piC);
    ch.SetBranchAddress("bestfit_xsec_normnccoh",&bestfit_xsec_normnccoh);
    ch.SetBranchAddress("bestfit_xsec_normncoth",&bestfit_xsec_normncoth);
    ch.SetBranchAddress("bestfit_xsec_normnuetonumu",&bestfit_xsec_normnuetonumu);
    ch.SetBranchAddress("bestfit_skdet_ccqerc0",&bestfit_skdet_ccqerc0);
    ch.SetBranchAddress("bestfit_skdet_ccqerc1",&bestfit_skdet_ccqerc1);
    ch.SetBranchAddress("bestfit_skdet_ccqerc2",&bestfit_skdet_ccqerc2);
    ch.SetBranchAddress("bestfit_skdet_ccqeoth",&bestfit_skdet_ccqeoth);
    ch.SetBranchAddress("bestfit_skdet_nuecc",&bestfit_skdet_nuecc);
    ch.SetBranchAddress("bestfit_fsi_000",&bestfit_fsi_000);
    ch.SetBranchAddress("bestfit_fsi_001",&bestfit_fsi_001);
    ch.SetBranchAddress("bestfit_fsi_002",&bestfit_fsi_002);
    ch.SetBranchAddress("bestfit_fsi_003",&bestfit_fsi_003);
    ch.SetBranchAddress("bestfit_fsi_004",&bestfit_fsi_004);
    ch.SetBranchAddress("bestfit_fsi_005",&bestfit_fsi_005);
    //errminospluss
    ch.SetBranchAddress("errminosplus_xsec_pf",&errminosplus_xsec_pf);
    ch.SetBranchAddress("errminosplus_xsec_wshape",&errminosplus_xsec_wshape);
    ch.SetBranchAddress("errminosplus_xsec_ccothshape",&errminosplus_xsec_ccothshape);
    ch.SetBranchAddress("errminosplus_xsec_cc1pishape",&errminosplus_xsec_cc1pishape);
    ch.SetBranchAddress("errminosplus_xsec_normcccoh",&errminosplus_xsec_normcccoh);
    ch.SetBranchAddress("errminosplus_xsec_normnc1piC",&errminosplus_xsec_normnc1piC);
    ch.SetBranchAddress("errminosplus_xsec_normnccoh",&errminosplus_xsec_normnccoh);
    ch.SetBranchAddress("errminosplus_xsec_normncoth",&errminosplus_xsec_normncoth);
    ch.SetBranchAddress("errminosplus_xsec_normnuetonumu",&errminosplus_xsec_normnuetonumu);
    ch.SetBranchAddress("errminosplus_skdet_ccqerc0",&errminosplus_skdet_ccqerc0);
    ch.SetBranchAddress("errminosplus_skdet_ccqerc1",&errminosplus_skdet_ccqerc1);
    ch.SetBranchAddress("errminosplus_skdet_ccqerc2",&errminosplus_skdet_ccqerc2);
    ch.SetBranchAddress("errminosplus_skdet_ccqeoth",&errminosplus_skdet_ccqeoth);
    ch.SetBranchAddress("errminosplus_skdet_nuecc",&errminosplus_skdet_nuecc);
    ch.SetBranchAddress("errminosplus_fsi_000",&errminosplus_fsi_000);
    ch.SetBranchAddress("errminosplus_fsi_001",&errminosplus_fsi_001);
    ch.SetBranchAddress("errminosplus_fsi_002",&errminosplus_fsi_002);
    ch.SetBranchAddress("errminosplus_fsi_003",&errminosplus_fsi_003);
    ch.SetBranchAddress("errminosplus_fsi_004",&errminosplus_fsi_004);
    ch.SetBranchAddress("errminosplus_fsi_005",&errminosplus_fsi_005);
    //errminosminuss
    ch.SetBranchAddress("errminosminus_xsec_pf",&errminosminus_xsec_pf);
    ch.SetBranchAddress("errminosminus_xsec_wshape",&errminosminus_xsec_wshape);
    ch.SetBranchAddress("errminosminus_xsec_ccothshape",&errminosminus_xsec_ccothshape);
    ch.SetBranchAddress("errminosminus_xsec_cc1pishape",&errminosminus_xsec_cc1pishape);
    ch.SetBranchAddress("errminosminus_xsec_normcccoh",&errminosminus_xsec_normcccoh);
    ch.SetBranchAddress("errminosminus_xsec_normnc1piC",&errminosminus_xsec_normnc1piC);
    ch.SetBranchAddress("errminosminus_xsec_normnccoh",&errminosminus_xsec_normnccoh);
    ch.SetBranchAddress("errminosminus_xsec_normncoth",&errminosminus_xsec_normncoth);
    ch.SetBranchAddress("errminosminus_xsec_normnuetonumu",&errminosminus_xsec_normnuetonumu);
    ch.SetBranchAddress("errminosminus_skdet_ccqerc0",&errminosminus_skdet_ccqerc0);
    ch.SetBranchAddress("errminosminus_skdet_ccqerc1",&errminosminus_skdet_ccqerc1);
    ch.SetBranchAddress("errminosminus_skdet_ccqerc2",&errminosminus_skdet_ccqerc2);
    ch.SetBranchAddress("errminosminus_skdet_ccqeoth",&errminosminus_skdet_ccqeoth);
    ch.SetBranchAddress("errminosminus_skdet_nuecc",&errminosminus_skdet_nuecc);
    ch.SetBranchAddress("errminosminus_fsi_000",&errminosminus_fsi_000);
    ch.SetBranchAddress("errminosminus_fsi_001",&errminosminus_fsi_001);
    ch.SetBranchAddress("errminosminus_fsi_002",&errminosminus_fsi_002);
    ch.SetBranchAddress("errminosminus_fsi_003",&errminosminus_fsi_003);
    ch.SetBranchAddress("errminosminus_fsi_004",&errminosminus_fsi_004);
    ch.SetBranchAddress("errminosminus_fsi_005",&errminosminus_fsi_005);
    //errmigrads
    ch.SetBranchAddress("errmigrad_xsec_pf",&errmigrad_xsec_pf);
    ch.SetBranchAddress("errmigrad_xsec_wshape",&errmigrad_xsec_wshape);
    ch.SetBranchAddress("errmigrad_xsec_ccothshape",&errmigrad_xsec_ccothshape);
    ch.SetBranchAddress("errmigrad_xsec_cc1pishape",&errmigrad_xsec_cc1pishape);
    ch.SetBranchAddress("errmigrad_xsec_normcccoh",&errmigrad_xsec_normcccoh);
    ch.SetBranchAddress("errmigrad_xsec_normnc1piC",&errmigrad_xsec_normnc1piC);
    ch.SetBranchAddress("errmigrad_xsec_normnccoh",&errmigrad_xsec_normnccoh);
    ch.SetBranchAddress("errmigrad_xsec_normncoth",&errmigrad_xsec_normncoth);
    ch.SetBranchAddress("errmigrad_xsec_normnuetonumu",&errmigrad_xsec_normnuetonumu);
    ch.SetBranchAddress("errmigrad_skdet_ccqerc0",&errmigrad_skdet_ccqerc0);
    ch.SetBranchAddress("errmigrad_skdet_ccqerc1",&errmigrad_skdet_ccqerc1);
    ch.SetBranchAddress("errmigrad_skdet_ccqerc2",&errmigrad_skdet_ccqerc2);
    ch.SetBranchAddress("errmigrad_skdet_ccqeoth",&errmigrad_skdet_ccqeoth);
    ch.SetBranchAddress("errmigrad_skdet_nuecc",&errmigrad_skdet_nuecc);
    ch.SetBranchAddress("errmigrad_fsi_000",&errmigrad_fsi_000);
    ch.SetBranchAddress("errmigrad_fsi_001",&errmigrad_fsi_001);
    ch.SetBranchAddress("errmigrad_fsi_002",&errmigrad_fsi_002);
    ch.SetBranchAddress("errmigrad_fsi_003",&errmigrad_fsi_003);
    ch.SetBranchAddress("errmigrad_fsi_004",&errmigrad_fsi_004);
    ch.SetBranchAddress("errmigrad_fsi_005",&errmigrad_fsi_005);
    //errhesses
    ch.SetBranchAddress("errhesse_xsec_pf",&errhesse_xsec_pf);
    ch.SetBranchAddress("errhesse_xsec_wshape",&errhesse_xsec_wshape);
    ch.SetBranchAddress("errhesse_xsec_ccothshape",&errhesse_xsec_ccothshape);
    ch.SetBranchAddress("errhesse_xsec_cc1pishape",&errhesse_xsec_cc1pishape);
    ch.SetBranchAddress("errhesse_xsec_normcccoh",&errhesse_xsec_normcccoh);
    ch.SetBranchAddress("errhesse_xsec_normnc1piC",&errhesse_xsec_normnc1piC);
    ch.SetBranchAddress("errhesse_xsec_normnccoh",&errhesse_xsec_normnccoh);
    ch.SetBranchAddress("errhesse_xsec_normncoth",&errhesse_xsec_normncoth);
    ch.SetBranchAddress("errhesse_xsec_normnuetonumu",&errhesse_xsec_normnuetonumu);
    ch.SetBranchAddress("errhesse_skdet_ccqerc0",&errhesse_skdet_ccqerc0);
    ch.SetBranchAddress("errhesse_skdet_ccqerc1",&errhesse_skdet_ccqerc1);
    ch.SetBranchAddress("errhesse_skdet_ccqerc2",&errhesse_skdet_ccqerc2);
    ch.SetBranchAddress("errhesse_skdet_ccqeoth",&errhesse_skdet_ccqeoth);
    ch.SetBranchAddress("errhesse_skdet_nuecc",&errhesse_skdet_nuecc);
    ch.SetBranchAddress("errhesse_fsi_000",&errhesse_fsi_000);
    ch.SetBranchAddress("errhesse_fsi_001",&errhesse_fsi_001);
    ch.SetBranchAddress("errhesse_fsi_002",&errhesse_fsi_002);
    ch.SetBranchAddress("errhesse_fsi_003",&errhesse_fsi_003);
    ch.SetBranchAddress("errhesse_fsi_004",&errhesse_fsi_004);
    ch.SetBranchAddress("errhesse_fsi_005",&errhesse_fsi_005);
  }

  //disable a few branches
  ch.SetBranchStatus("*",0);
  ch.SetBranchStatus("iexp",1);
  ch.SetBranchStatus("fitstatus",1);
  ch.SetBranchStatus("migradfitstatus",1);
  ch.SetBranchStatus("hessefitstatus",1);
  ch.SetBranchStatus("fittime",1);
  ch.SetBranchStatus("bestfit_chisq",1);
  ch.SetBranchStatus("true_chisq",1);
  ch.SetBranchStatus("nsk_nom",1);
  ch.SetBranchStatus("nsk_bf",1);
  //oscillation parameters
  ch.SetBranchStatus("input_mixing12",1);
  ch.SetBranchStatus("input_mixing13",1);
  ch.SetBranchStatus("input_mixing23",1);
  ch.SetBranchStatus("input_deltacp",1);
  ch.SetBranchStatus("input_dm12sq",1);
  ch.SetBranchStatus("input_dm23sq",1);
  ch.SetBranchStatus("bestfit_mixing12",1);
  ch.SetBranchStatus("bestfit_mixing13",1);
  ch.SetBranchStatus("bestfit_mixing23",1);
  ch.SetBranchStatus("bestfit_deltacp",1);
  ch.SetBranchStatus("bestfit_dm12sq",1);
  ch.SetBranchStatus("bestfit_dm23sq",1);
  ch.SetBranchStatus("errminosplus_mixing12",1);
  ch.SetBranchStatus("errminosplus_mixing13",1);
  ch.SetBranchStatus("errminosplus_mixing23",1);
  ch.SetBranchStatus("errminosplus_deltacp",1);
  ch.SetBranchStatus("errminosplus_dm12sq",1);
  ch.SetBranchStatus("errminosplus_dm23sq",1);
  ch.SetBranchStatus("errminosminus_mixing12",1);
  ch.SetBranchStatus("errminosminus_mixing13",1);
  ch.SetBranchStatus("errminosminus_mixing23",1);
  ch.SetBranchStatus("errminosminus_deltacp",1);
  ch.SetBranchStatus("errminosminus_dm12sq",1);
  ch.SetBranchStatus("errminosminus_dm23sq",1);
  ch.SetBranchStatus("errmigrad_mixing12",1);
  ch.SetBranchStatus("errmigrad_mixing13",1);
  ch.SetBranchStatus("errmigrad_mixing23",1);
  ch.SetBranchStatus("errmigrad_deltacp",1);
  ch.SetBranchStatus("errmigrad_dm12sq",1);
  ch.SetBranchStatus("errmigrad_dm23sq",1);
  ch.SetBranchStatus("errhesse_mixing12",1);
  ch.SetBranchStatus("errhesse_mixing13",1);
  ch.SetBranchStatus("errhesse_mixing23",1);
  ch.SetBranchStatus("errhesse_deltacp",1);
  ch.SetBranchStatus("errhesse_dm12sq",1);
  ch.SetBranchStatus("errhesse_dm23sq",1);

  if(fit_systematics_mode >= 1) {
    //inputs
    ch.SetBranchStatus("input_xsec_sf",1);
    ch.SetBranchStatus("input_xsec_binde",1);
    ch.SetBranchStatus("input_xsec_pilessdelta",1);
    ch.SetBranchStatus("input_skdet_ccnqe",1);
    ch.SetBranchStatus("input_skdet_nc",1);
    ch.SetBranchStatus("input_skdet_energyscale",1);
    //bestfits
    ch.SetBranchStatus("bestfit_xsec_sf",1);
    ch.SetBranchStatus("bestfit_xsec_binde",1);
    ch.SetBranchStatus("bestfit_xsec_pilessdelta",1);
    ch.SetBranchStatus("bestfit_skdet_ccnqe",1);
    ch.SetBranchStatus("bestfit_skdet_nc",1);
    ch.SetBranchStatus("bestfit_skdet_energyscale",1);
    //errminospluss
    ch.SetBranchStatus("errminosplus_xsec_sf",1);
    ch.SetBranchStatus("errminosplus_xsec_binde",1);
    ch.SetBranchStatus("errminosplus_xsec_pilessdelta",1);
    ch.SetBranchStatus("errminosplus_skdet_ccnqe",1);
    ch.SetBranchStatus("errminosplus_skdet_nc",1);
    ch.SetBranchStatus("errminosplus_skdet_energyscale",1);
    //errminosminuss
    ch.SetBranchStatus("errminosminus_xsec_sf",1);
    ch.SetBranchStatus("errminosminus_xsec_binde",1);
    ch.SetBranchStatus("errminosminus_xsec_pilessdelta",1);
    ch.SetBranchStatus("errminosminus_skdet_ccnqe",1);
    ch.SetBranchStatus("errminosminus_skdet_nc",1);
    ch.SetBranchStatus("errminosminus_skdet_energyscale",1);
    //errmigrads
    ch.SetBranchStatus("errmigrad_xsec_sf",1);
    ch.SetBranchStatus("errmigrad_xsec_binde",1);
    ch.SetBranchStatus("errmigrad_xsec_pilessdelta",1);
    ch.SetBranchStatus("errmigrad_skdet_ccnqe",1);
    ch.SetBranchStatus("errmigrad_skdet_nc",1);
    ch.SetBranchStatus("errmigrad_skdet_energyscale",1);
    //errhesses
    ch.SetBranchStatus("errhesse_xsec_sf",1);
    ch.SetBranchStatus("errhesse_xsec_binde",1);
    ch.SetBranchStatus("errhesse_xsec_pilessdelta",1);
    ch.SetBranchStatus("errhesse_skdet_ccnqe",1);
    ch.SetBranchStatus("errhesse_skdet_nc",1);
    ch.SetBranchStatus("errhesse_skdet_energyscale",1);

    //inputs
    ch.SetBranchStatus("input_banff_00",1);
    ch.SetBranchStatus("input_banff_01",1);
    ch.SetBranchStatus("input_banff_02",1);
    ch.SetBranchStatus("input_banff_03",1);
    ch.SetBranchStatus("input_banff_04",1);
    ch.SetBranchStatus("input_banff_05",1);
    ch.SetBranchStatus("input_banff_06",1);
    ch.SetBranchStatus("input_banff_07",1);
    ch.SetBranchStatus("input_banff_08",1);
    ch.SetBranchStatus("input_banff_09",1);
    ch.SetBranchStatus("input_banff_10",1);
    ch.SetBranchStatus("input_banff_11",1);
    ch.SetBranchStatus("input_banff_12",1);
    ch.SetBranchStatus("input_banff_13",1);
    ch.SetBranchStatus("input_banff_14",1);
    ch.SetBranchStatus("input_banff_15",1);
    ch.SetBranchStatus("input_banff_maqe",1);
    ch.SetBranchStatus("input_banff_mares",1);
    ch.SetBranchStatus("input_banff_ccqeE0",1);
    ch.SetBranchStatus("input_banff_ccqeE1",1);
    ch.SetBranchStatus("input_banff_ccqeE2",1);
    ch.SetBranchStatus("input_banff_cc1piE0",1);
    ch.SetBranchStatus("input_banff_cc1piE1",1);
    //bestfits
    ch.SetBranchStatus("bestfit_banff_00",1);
    ch.SetBranchStatus("bestfit_banff_01",1);
    ch.SetBranchStatus("bestfit_banff_02",1);
    ch.SetBranchStatus("bestfit_banff_03",1);
    ch.SetBranchStatus("bestfit_banff_04",1);
    ch.SetBranchStatus("bestfit_banff_05",1);
    ch.SetBranchStatus("bestfit_banff_06",1);
    ch.SetBranchStatus("bestfit_banff_07",1);
    ch.SetBranchStatus("bestfit_banff_08",1);
    ch.SetBranchStatus("bestfit_banff_09",1);
    ch.SetBranchStatus("bestfit_banff_10",1);
    ch.SetBranchStatus("bestfit_banff_11",1);
    ch.SetBranchStatus("bestfit_banff_12",1);
    ch.SetBranchStatus("bestfit_banff_13",1);
    ch.SetBranchStatus("bestfit_banff_14",1);
    ch.SetBranchStatus("bestfit_banff_15",1);
    ch.SetBranchStatus("bestfit_banff_maqe",1);
    ch.SetBranchStatus("bestfit_banff_mares",1);
    ch.SetBranchStatus("bestfit_banff_ccqeE0",1);
    ch.SetBranchStatus("bestfit_banff_ccqeE1",1);
    ch.SetBranchStatus("bestfit_banff_ccqeE2",1);
    ch.SetBranchStatus("bestfit_banff_cc1piE0",1);
    ch.SetBranchStatus("bestfit_banff_cc1piE1",1);
    //errminospluss
    ch.SetBranchStatus("errminosplus_banff_00",1);
    ch.SetBranchStatus("errminosplus_banff_01",1);
    ch.SetBranchStatus("errminosplus_banff_02",1);
    ch.SetBranchStatus("errminosplus_banff_03",1);
    ch.SetBranchStatus("errminosplus_banff_04",1);
    ch.SetBranchStatus("errminosplus_banff_05",1);
    ch.SetBranchStatus("errminosplus_banff_06",1);
    ch.SetBranchStatus("errminosplus_banff_07",1);
    ch.SetBranchStatus("errminosplus_banff_08",1);
    ch.SetBranchStatus("errminosplus_banff_09",1);
    ch.SetBranchStatus("errminosplus_banff_10",1);
    ch.SetBranchStatus("errminosplus_banff_11",1);
    ch.SetBranchStatus("errminosplus_banff_12",1);
    ch.SetBranchStatus("errminosplus_banff_13",1);
    ch.SetBranchStatus("errminosplus_banff_14",1);
    ch.SetBranchStatus("errminosplus_banff_15",1);
    ch.SetBranchStatus("errminosplus_banff_maqe",1);
    ch.SetBranchStatus("errminosplus_banff_mares",1);
    ch.SetBranchStatus("errminosplus_banff_ccqeE0",1);
    ch.SetBranchStatus("errminosplus_banff_ccqeE1",1);
    ch.SetBranchStatus("errminosplus_banff_ccqeE2",1);
    ch.SetBranchStatus("errminosplus_banff_cc1piE0",1);
    ch.SetBranchStatus("errminosplus_banff_cc1piE1",1);

    //errminosminuss
    ch.SetBranchStatus("errminosminus_banff_00",1);
    ch.SetBranchStatus("errminosminus_banff_01",1);
    ch.SetBranchStatus("errminosminus_banff_02",1);
    ch.SetBranchStatus("errminosminus_banff_03",1);
    ch.SetBranchStatus("errminosminus_banff_04",1);
    ch.SetBranchStatus("errminosminus_banff_05",1);
    ch.SetBranchStatus("errminosminus_banff_06",1);
    ch.SetBranchStatus("errminosminus_banff_07",1);
    ch.SetBranchStatus("errminosminus_banff_08",1);
    ch.SetBranchStatus("errminosminus_banff_09",1);
    ch.SetBranchStatus("errminosminus_banff_10",1);
    ch.SetBranchStatus("errminosminus_banff_11",1);
    ch.SetBranchStatus("errminosminus_banff_12",1);
    ch.SetBranchStatus("errminosminus_banff_13",1);
    ch.SetBranchStatus("errminosminus_banff_14",1);
    ch.SetBranchStatus("errminosminus_banff_15",1);
    ch.SetBranchStatus("errminosminus_banff_maqe",1);
    ch.SetBranchStatus("errminosminus_banff_mares",1);
    ch.SetBranchStatus("errminosminus_banff_ccqeE0",1);
    ch.SetBranchStatus("errminosminus_banff_ccqeE1",1);
    ch.SetBranchStatus("errminosminus_banff_ccqeE2",1);
    ch.SetBranchStatus("errminosminus_banff_cc1piE0",1);
    ch.SetBranchStatus("errminosminus_banff_cc1piE1",1);
    //errmigrad
    ch.SetBranchStatus("errmigrad_banff_00",1);
    ch.SetBranchStatus("errmigrad_banff_01",1);
    ch.SetBranchStatus("errmigrad_banff_02",1);
    ch.SetBranchStatus("errmigrad_banff_03",1);
    ch.SetBranchStatus("errmigrad_banff_04",1);
    ch.SetBranchStatus("errmigrad_banff_05",1);
    ch.SetBranchStatus("errmigrad_banff_06",1);
    ch.SetBranchStatus("errmigrad_banff_07",1);
    ch.SetBranchStatus("errmigrad_banff_08",1);
    ch.SetBranchStatus("errmigrad_banff_09",1);
    ch.SetBranchStatus("errmigrad_banff_10",1);
    ch.SetBranchStatus("errmigrad_banff_11",1);
    ch.SetBranchStatus("errmigrad_banff_12",1);
    ch.SetBranchStatus("errmigrad_banff_13",1);
    ch.SetBranchStatus("errmigrad_banff_14",1);
    ch.SetBranchStatus("errmigrad_banff_15",1);
    ch.SetBranchStatus("errmigrad_banff_maqe",1);
    ch.SetBranchStatus("errmigrad_banff_mares",1);
    ch.SetBranchStatus("errmigrad_banff_ccqeE0",1);
    ch.SetBranchStatus("errmigrad_banff_ccqeE1",1);
    ch.SetBranchStatus("errmigrad_banff_ccqeE2",1);
    ch.SetBranchStatus("errmigrad_banff_cc1piE0",1);
    ch.SetBranchStatus("errmigrad_banff_cc1piE1",1);
    //errhesse
    ch.SetBranchStatus("errhesse_banff_00",1);
    ch.SetBranchStatus("errhesse_banff_01",1);
    ch.SetBranchStatus("errhesse_banff_02",1);
    ch.SetBranchStatus("errhesse_banff_03",1);
    ch.SetBranchStatus("errhesse_banff_04",1);
    ch.SetBranchStatus("errhesse_banff_05",1);
    ch.SetBranchStatus("errhesse_banff_06",1);
    ch.SetBranchStatus("errhesse_banff_07",1);
    ch.SetBranchStatus("errhesse_banff_08",1);
    ch.SetBranchStatus("errhesse_banff_09",1);
    ch.SetBranchStatus("errhesse_banff_10",1);
    ch.SetBranchStatus("errhesse_banff_11",1);
    ch.SetBranchStatus("errhesse_banff_12",1);
    ch.SetBranchStatus("errhesse_banff_13",1);
    ch.SetBranchStatus("errhesse_banff_14",1);
    ch.SetBranchStatus("errhesse_banff_15",1);
    ch.SetBranchStatus("errhesse_banff_maqe",1);
    ch.SetBranchStatus("errhesse_banff_mares",1);
    ch.SetBranchStatus("errhesse_banff_ccqeE0",1);
    ch.SetBranchStatus("errhesse_banff_ccqeE1",1);
    ch.SetBranchStatus("errhesse_banff_ccqeE2",1);
    ch.SetBranchStatus("errhesse_banff_cc1piE0",1);
    ch.SetBranchStatus("errhesse_banff_cc1piE1",1);
  }
  if(fit_systematics_mode == 2) {
    //inputs
    ch.SetBranchStatus("input_xsec_pf",1);
    ch.SetBranchStatus("input_xsec_wshape",1);
    ch.SetBranchStatus("input_xsec_ccothshape",1);
    ch.SetBranchStatus("input_xsec_cc1pishape",1);
    ch.SetBranchStatus("input_xsec_normcccoh",1);
    ch.SetBranchStatus("input_xsec_normnc1piC",1);
    ch.SetBranchStatus("input_xsec_normnccoh",1);
    ch.SetBranchStatus("input_xsec_normncoth",1);
    ch.SetBranchStatus("input_xsec_normnuetonumu",1);
    ch.SetBranchStatus("input_skdet_ccqerc0",1);
    ch.SetBranchStatus("input_skdet_ccqerc1",1);
    ch.SetBranchStatus("input_skdet_ccqerc2",1);
    ch.SetBranchStatus("input_skdet_ccqeoth",1);
    ch.SetBranchStatus("input_skdet_nuecc",1);
    ch.SetBranchStatus("input_fsi_000",1);
    ch.SetBranchStatus("input_fsi_001",1);
    ch.SetBranchStatus("input_fsi_002",1);
    ch.SetBranchStatus("input_fsi_003",1);
    ch.SetBranchStatus("input_fsi_004",1);
    ch.SetBranchStatus("input_fsi_005",1);
    //bestfits
    ch.SetBranchStatus("bestfit_xsec_pf",1);
    ch.SetBranchStatus("bestfit_xsec_wshape",1);
    ch.SetBranchStatus("bestfit_xsec_ccothshape",1);
    ch.SetBranchStatus("bestfit_xsec_cc1pishape",1);
    ch.SetBranchStatus("bestfit_xsec_normcccoh",1);
    ch.SetBranchStatus("bestfit_xsec_normnc1piC",1);
    ch.SetBranchStatus("bestfit_xsec_normnccoh",1);
    ch.SetBranchStatus("bestfit_xsec_normncoth",1);
    ch.SetBranchStatus("bestfit_xsec_normnuetonumu",1);
    ch.SetBranchStatus("bestfit_skdet_ccqerc0",1);
    ch.SetBranchStatus("bestfit_skdet_ccqerc1",1);
    ch.SetBranchStatus("bestfit_skdet_ccqerc2",1);
    ch.SetBranchStatus("bestfit_skdet_ccqeoth",1);
    ch.SetBranchStatus("bestfit_skdet_nuecc",1);
    ch.SetBranchStatus("bestfit_fsi_000",1);
    ch.SetBranchStatus("bestfit_fsi_001",1);
    ch.SetBranchStatus("bestfit_fsi_002",1);
    ch.SetBranchStatus("bestfit_fsi_003",1);
    ch.SetBranchStatus("bestfit_fsi_004",1);
    ch.SetBranchStatus("bestfit_fsi_005",1);
    //errminospluss
    ch.SetBranchStatus("errminosplus_xsec_pf",1);
    ch.SetBranchStatus("errminosplus_xsec_wshape",1);
    ch.SetBranchStatus("errminosplus_xsec_ccothshape",1);
    ch.SetBranchStatus("errminosplus_xsec_cc1pishape",1);
    ch.SetBranchStatus("errminosplus_xsec_normcccoh",1);
    ch.SetBranchStatus("errminosplus_xsec_normnc1piC",1);
    ch.SetBranchStatus("errminosplus_xsec_normnccoh",1);
    ch.SetBranchStatus("errminosplus_xsec_normncoth",1);
    ch.SetBranchStatus("errminosplus_xsec_normnuetonumu",1);
    ch.SetBranchStatus("errminosplus_skdet_ccqerc0",1);
    ch.SetBranchStatus("errminosplus_skdet_ccqerc1",1);
    ch.SetBranchStatus("errminosplus_skdet_ccqerc2",1);
    ch.SetBranchStatus("errminosplus_skdet_ccqeoth",1);
    ch.SetBranchStatus("errminosplus_skdet_nuecc",1);
    ch.SetBranchStatus("errminosplus_fsi_000",1);
    ch.SetBranchStatus("errminosplus_fsi_001",1);
    ch.SetBranchStatus("errminosplus_fsi_002",1);
    ch.SetBranchStatus("errminosplus_fsi_003",1);
    ch.SetBranchStatus("errminosplus_fsi_004",1);
    ch.SetBranchStatus("errminosplus_fsi_005",1);
    //errminosminuss
    ch.SetBranchStatus("errminosminus_xsec_pf",1);
    ch.SetBranchStatus("errminosminus_xsec_wshape",1);
    ch.SetBranchStatus("errminosminus_xsec_ccothshape",1);
    ch.SetBranchStatus("errminosminus_xsec_cc1pishape",1);
    ch.SetBranchStatus("errminosminus_xsec_normcccoh",1);
    ch.SetBranchStatus("errminosminus_xsec_normnc1piC",1);
    ch.SetBranchStatus("errminosminus_xsec_normnccoh",1);
    ch.SetBranchStatus("errminosminus_xsec_normncoth",1);
    ch.SetBranchStatus("errminosminus_xsec_normnuetonumu",1);
    ch.SetBranchStatus("errminosminus_skdet_ccqerc0",1);
    ch.SetBranchStatus("errminosminus_skdet_ccqerc1",1);
    ch.SetBranchStatus("errminosminus_skdet_ccqerc2",1);
    ch.SetBranchStatus("errminosminus_skdet_ccqeoth",1);
    ch.SetBranchStatus("errminosminus_skdet_nuecc",1);
    ch.SetBranchStatus("errminosminus_fsi_000",1);
    ch.SetBranchStatus("errminosminus_fsi_001",1);
    ch.SetBranchStatus("errminosminus_fsi_002",1);
    ch.SetBranchStatus("errminosminus_fsi_003",1);
    ch.SetBranchStatus("errminosminus_fsi_004",1);
    ch.SetBranchStatus("errminosminus_fsi_005",1);
    //errmigrads
    ch.SetBranchStatus("errmigrad_xsec_pf",1);
    ch.SetBranchStatus("errmigrad_xsec_wshape",1);
    ch.SetBranchStatus("errmigrad_xsec_ccothshape",1);
    ch.SetBranchStatus("errmigrad_xsec_cc1pishape",1);
    ch.SetBranchStatus("errmigrad_xsec_normcccoh",1);
    ch.SetBranchStatus("errmigrad_xsec_normnc1piC",1);
    ch.SetBranchStatus("errmigrad_xsec_normnccoh",1);
    ch.SetBranchStatus("errmigrad_xsec_normncoth",1);
    ch.SetBranchStatus("errmigrad_xsec_normnuetonumu",1);
    ch.SetBranchStatus("errmigrad_skdet_ccqerc0",1);
    ch.SetBranchStatus("errmigrad_skdet_ccqerc1",1);
    ch.SetBranchStatus("errmigrad_skdet_ccqerc2",1);
    ch.SetBranchStatus("errmigrad_skdet_ccqeoth",1);
    ch.SetBranchStatus("errmigrad_skdet_nuecc",1);
    ch.SetBranchStatus("errmigrad_fsi_000",1);
    ch.SetBranchStatus("errmigrad_fsi_001",1);
    ch.SetBranchStatus("errmigrad_fsi_002",1);
    ch.SetBranchStatus("errmigrad_fsi_003",1);
    ch.SetBranchStatus("errmigrad_fsi_004",1);
    ch.SetBranchStatus("errmigrad_fsi_005",1);
    //errhesses
    ch.SetBranchStatus("errhesse_xsec_pf",1);
    ch.SetBranchStatus("errhesse_xsec_wshape",1);
    ch.SetBranchStatus("errhesse_xsec_ccothshape",1);
    ch.SetBranchStatus("errhesse_xsec_cc1pishape",1);
    ch.SetBranchStatus("errhesse_xsec_normcccoh",1);
    ch.SetBranchStatus("errhesse_xsec_normnc1piC",1);
    ch.SetBranchStatus("errhesse_xsec_normnccoh",1);
    ch.SetBranchStatus("errhesse_xsec_normncoth",1);
    ch.SetBranchStatus("errhesse_xsec_normnuetonumu",1);
    ch.SetBranchStatus("errhesse_skdet_ccqerc0",1);
    ch.SetBranchStatus("errhesse_skdet_ccqerc1",1);
    ch.SetBranchStatus("errhesse_skdet_ccqerc2",1);
    ch.SetBranchStatus("errhesse_skdet_ccqeoth",1);
    ch.SetBranchStatus("errhesse_skdet_nuecc",1);
    ch.SetBranchStatus("errhesse_fsi_000",1);
    ch.SetBranchStatus("errhesse_fsi_001",1);
    ch.SetBranchStatus("errhesse_fsi_002",1);
    ch.SetBranchStatus("errhesse_fsi_003",1);
    ch.SetBranchStatus("errhesse_fsi_004",1);
    ch.SetBranchStatus("errhesse_fsi_005",1);
  }
 

  // Book histograms
  TH1I * h_fitstat        = new TH1I("h_fitstat","h_fitstat",11,-1,10);
  TH1I * h_migradfitstat  = new TH1I("h_migradfitstat","h_migradfitstat",14,-1,13);
  TH1I * h_hessefitstat   = new TH1I("h_hessefitstat","h_hessefitstat",14,-1,13);
  TH1D * h_fittime  = new TH1D("h_fittime","h_fittime",1000,0,1000);
  TH1D * h_nsk_nom  = new TH1D("h_nsk_nom","h_nsk_nom",150,0,150);
  TH1D * h_nsk_bf   = new TH1D("h_nsk_bf","h_nsk_bf",150,0,150);
  TH1D * h_dchisq   = new TH1D("h_dchisq","h_dchisq",300,-5,25);
  TH1D * h_mixing23 = new TH1D("h_mixing23","h_mixing23",500,0.5,1.0);
  TH1D * h_dm23sq   = new TH1D("h_dm23sq","h_dm23sq",500,1.4E-3,3.4E-3);
  //pull histograms
  TH1D * hpull_mixing23 = new TH1D("hpull_mixing23","hpull_mixing23",250,-5,5);
  TH1D * hpull_dm23sq   = new TH1D("hpull_dm23sq","hpull_dm23sq",250,-5,5);
  //error histograms
  TH1D * herr_mixing23 = new TH1D("herr_mixing23","herr_mixing23",300,0,0.05);
  TH1D * herr_dm23sq   = new TH1D("herr_dm23sq","herr_dm23sq",100,0,1E-4);
  //distribution histograms
  TH1D * h_xsec_sf = new TH1D("h_xsec_sf","h_xsec_sf",200,-3,3);
  TH1D * h_xsec_binde = new TH1D("h_xsec_binde","h_xsec_binde",200,-3,3);
  TH1D * h_xsec_pilessdelta = new TH1D("h_xsec_pilessdelta","h_xsec_pilessdelta",200,-3,3);
  TH1D * h_skdet_ccnqe = new TH1D("h_skdet_ccnqe","h_skdet_ccnqe",200,-3,3);
  TH1D * h_skdet_nc = new TH1D("h_skdet_nc","h_skdet_nc",200,-3,3);
  TH1D * h_skdet_energyscale = new TH1D("h_skdet_energyscale","h_skdet_energyscale",200,-3,3);
  //pull histograms
  TH1D * hpull_xsec_sf = new TH1D("hpull_xsec_sf","hpull_xsec_sf",100,-5,5);
  TH1D * hpull_xsec_binde = new TH1D("hpull_xsec_binde","hpull_xsec_binde",100,-10,10);
  TH1D * hpull_xsec_pilessdelta = new TH1D("hpull_xsec_pilessdelta","hpull_xsec_pilessdelta",100,-5,5);
  TH1D * hpull_skdet_ccnqe = new TH1D("hpull_skdet_ccnqe","hpull_skdet_ccnqe",100,-5,5);
  TH1D * hpull_skdet_nc = new TH1D("hpull_skdet_nc","hpull_skdet_nc",100,-5,5);
  TH1D * hpull_skdet_energyscale = new TH1D("hpull_skdet_energyscale","hpull_skdet_energyscale",100,-5,5);
  //alterative pull histograms
  TH1D * haltpull_xsec_sf = new TH1D("haltpull_xsec_sf","haltpull_xsec_sf",100,-5,5);
  TH1D * haltpull_xsec_binde = new TH1D("haltpull_xsec_binde","haltpull_xsec_binde",100,-10,10);
  TH1D * haltpull_xsec_pilessdelta = new TH1D("haltpull_xsec_pilessdelta","haltpull_xsec_pilessdelta",100,-5,5);
  TH1D * haltpull_skdet_ccnqe = new TH1D("haltpull_skdet_ccnqe","haltpull_skdet_ccnqe",100,-5,5);
  TH1D * haltpull_skdet_nc = new TH1D("haltpull_skdet_nc","haltpull_skdet_nc",100,-5,5);
  TH1D * haltpull_skdet_energyscale = new TH1D("haltpull_skdet_energyscale","haltpull_skdet_energyscale",100,-5,5);
  //distribution histograms
  TH1D * h_banff_00 = new TH1D("h_banff_00","h_banff_00",200,-3,3);
  TH1D * h_banff_01 = new TH1D("h_banff_01","h_banff_01",200,-3,3);
  TH1D * h_banff_02 = new TH1D("h_banff_02","h_banff_02",200,-3,3);
  TH1D * h_banff_03 = new TH1D("h_banff_03","h_banff_03",200,-3,3);
  TH1D * h_banff_04 = new TH1D("h_banff_04","h_banff_04",200,-3,3);
  TH1D * h_banff_05 = new TH1D("h_banff_05","h_banff_05",200,-3,3);
  TH1D * h_banff_06 = new TH1D("h_banff_06","h_banff_06",200,-3,3);
  TH1D * h_banff_07 = new TH1D("h_banff_07","h_banff_07",200,-3,3);
  TH1D * h_banff_08 = new TH1D("h_banff_08","h_banff_08",200,-3,3);
  TH1D * h_banff_09 = new TH1D("h_banff_09","h_banff_09",200,-3,3);
  TH1D * h_banff_10 = new TH1D("h_banff_10","h_banff_10",200,-3,3);
  TH1D * h_banff_11 = new TH1D("h_banff_11","h_banff_11",200,-3,3);
  TH1D * h_banff_12 = new TH1D("h_banff_12","h_banff_12",200,-3,3);
  TH1D * h_banff_13 = new TH1D("h_banff_13","h_banff_13",200,-3,3);
  TH1D * h_banff_14 = new TH1D("h_banff_14","h_banff_14",200,-3,3);
  TH1D * h_banff_15 = new TH1D("h_banff_15","h_banff_15",200,-3,3);
  TH1D * h_banff_maqe = new TH1D("h_banff_maqe","h_banff_maqe",200,-3,3);
  TH1D * h_banff_mares = new TH1D("h_banff_mares","h_banff_mares",200,-3,3);
  TH1D * h_banff_ccqeE0 = new TH1D("h_banff_ccqeE0","h_banff_ccqeE0",200,-3,3);
  TH1D * h_banff_ccqeE1 = new TH1D("h_banff_ccqeE1","h_banff_ccqeE1",200,-3,3);
  TH1D * h_banff_ccqeE2 = new TH1D("h_banff_ccqeE2","h_banff_ccqeE2",200,-3,3);
  TH1D * h_banff_cc1piE0 = new TH1D("h_banff_cc1piE0","h_banff_cc1piE0",200,-3,3);
  TH1D * h_banff_cc1piE1 = new TH1D("h_banff_cc1piE1","h_banff_cc1piE1",200,-3,3);
  //pull histograms
  TH1D * hpull_banff_00 = new TH1D("hpull_banff_00","hpull_banff_00",100,-5,5);
  TH1D * hpull_banff_01 = new TH1D("hpull_banff_01","hpull_banff_01",100,-5,5);
  TH1D * hpull_banff_02 = new TH1D("hpull_banff_02","hpull_banff_02",100,-5,5);
  TH1D * hpull_banff_03 = new TH1D("hpull_banff_03","hpull_banff_03",100,-5,5);
  TH1D * hpull_banff_04 = new TH1D("hpull_banff_04","hpull_banff_04",100,-5,5);
  TH1D * hpull_banff_05 = new TH1D("hpull_banff_05","hpull_banff_05",100,-5,5);
  TH1D * hpull_banff_06 = new TH1D("hpull_banff_06","hpull_banff_06",100,-5,5);
  TH1D * hpull_banff_07 = new TH1D("hpull_banff_07","hpull_banff_07",100,-5,5);
  TH1D * hpull_banff_08 = new TH1D("hpull_banff_08","hpull_banff_08",100,-5,5);
  TH1D * hpull_banff_09 = new TH1D("hpull_banff_09","hpull_banff_09",100,-5,5);
  TH1D * hpull_banff_10 = new TH1D("hpull_banff_10","hpull_banff_10",100,-5,5);
  TH1D * hpull_banff_11 = new TH1D("hpull_banff_11","hpull_banff_11",100,-5,5);
  TH1D * hpull_banff_12 = new TH1D("hpull_banff_12","hpull_banff_12",100,-5,5);
  TH1D * hpull_banff_13 = new TH1D("hpull_banff_13","hpull_banff_13",100,-5,5);
  TH1D * hpull_banff_14 = new TH1D("hpull_banff_14","hpull_banff_14",100,-5,5);
  TH1D * hpull_banff_15 = new TH1D("hpull_banff_15","hpull_banff_15",100,-5,5);
  TH1D * hpull_banff_maqe = new TH1D("hpull_banff_maqe","hpull_banff_maqe",100,-5,5);
  TH1D * hpull_banff_mares = new TH1D("hpull_banff_mares","hpull_banff_mares",100,-5,5);
  TH1D * hpull_banff_ccqeE0 = new TH1D("hpull_banff_ccqeE0","hpull_banff_ccqeE0",100,-5,5);
  TH1D * hpull_banff_ccqeE1 = new TH1D("hpull_banff_ccqeE1","hpull_banff_ccqeE1",100,-5,5);
  TH1D * hpull_banff_ccqeE2 = new TH1D("hpull_banff_ccqeE2","hpull_banff_ccqeE2",100,-5,5);
  TH1D * hpull_banff_cc1piE0 = new TH1D("hpull_banff_cc1piE0","hpull_banff_cc1piE0",100,-5,5);
  TH1D * hpull_banff_cc1piE1 = new TH1D("hpull_banff_cc1piE1","hpull_banff_cc1piE1",100,-5,5);
  //alt pull
  TH1D * haltpull_banff_maqe = new TH1D("haltpull_banff_maqe","haltpull_banff_maqe",100,-5,5);
  TH1D * haltpull_banff_ccqeE0 = new TH1D("haltpull_banff_ccqeE0","haltpull_banff_ccqeE0",100,-5,5);
  //fit_systematics_mode == 2
  TH1D * h_xsec_pf = new TH1D("h_xsec_pf","h_xsec_pf",50,-5,5);
  TH1D * h_xsec_wshape = new TH1D("h_xsec_wshape","h_xsec_wshape",50,-5,5);
  TH1D * h_xsec_ccothshape = new TH1D("h_xsec_ccothshape","h_xsec_ccothshape",50,-5,5);
  TH1D * h_xsec_cc1pishape = new TH1D("h_xsec_cc1pishape","h_xsec_cc1pishape",50,-5,5);
  TH1D * h_xsec_normcccoh = new TH1D("h_xsec_normcccoh","h_xsec_normcccoh",50,-5,5);
  TH1D * h_xsec_normnc1piC = new TH1D("h_xsec_normnc1piC","h_xsec_normnc1piC",50,-5,5);
  TH1D * h_xsec_normnccoh = new TH1D("h_xsec_normnccoh","h_xsec_normnccoh",50,-5,5);
  TH1D * h_xsec_normncoth = new TH1D("h_xsec_normncoth","h_xsec_normncoth",50,-5,5);
  TH1D * h_xsec_normnuetonumu = new TH1D("h_xsec_normnuetonumu","h_xsec_normnuetonumu",50,-5,5);
  TH1D * h_skdet_ccqerc0 = new TH1D("h_skdet_ccqerc0","h_skdet_ccqerc0",50,-5,5);
  TH1D * h_skdet_ccqerc1 = new TH1D("h_skdet_ccqerc1","h_skdet_ccqerc1",50,-5,5);
  TH1D * h_skdet_ccqerc2 = new TH1D("h_skdet_ccqerc2","h_skdet_ccqerc2",50,-5,5);
  TH1D * h_skdet_ccqeoth = new TH1D("h_skdet_ccqeoth","h_skdet_ccqeoth",50,-5,5);
  TH1D * h_skdet_nuecc = new TH1D("h_skdet_nuecc","h_skdet_nuecc",50,-5,5);
  TH1D * h_fsi_000 = new TH1D("h_fsi_000","h_fsi_000",50,-5,5);
  TH1D * h_fsi_001 = new TH1D("h_fsi_001","h_fsi_001",50,-5,5);
  TH1D * h_fsi_002 = new TH1D("h_fsi_002","h_fsi_002",50,-5,5);
  TH1D * h_fsi_003 = new TH1D("h_fsi_003","h_fsi_003",50,-5,5);
  TH1D * h_fsi_004 = new TH1D("h_fsi_004","h_fsi_004",50,-5,5);
  TH1D * h_fsi_005 = new TH1D("h_fsi_005","h_fsi_005",50,-5,5); 
  //pull histograms
  TH1D * hpull_xsec_pf = new TH1D("hpull_xsec_pf","hpull_xsec_pf",50,-5,5);
  TH1D * hpull_xsec_wshape = new TH1D("hpull_xsec_wshape","hpull_xsec_wshape",50,-5,5);
  TH1D * hpull_xsec_ccothshape = new TH1D("hpull_xsec_ccothshape","hpull_xsec_ccothshape",50,-5,5);
  TH1D * hpull_xsec_cc1pishape = new TH1D("hpull_xsec_cc1pishape","hpull_xsec_cc1pishape",50,-5,5);
  TH1D * hpull_xsec_normcccoh = new TH1D("hpull_xsec_normcccoh","hpull_xsec_normcccoh",50,-5,5);
  TH1D * hpull_xsec_normnc1piC = new TH1D("hpull_xsec_normnc1piC","hpull_xsec_normnc1piC",50,-5,5);
  TH1D * hpull_xsec_normnccoh = new TH1D("hpull_xsec_normnccoh","hpull_xsec_normnccoh",50,-5,5);
  TH1D * hpull_xsec_normncoth = new TH1D("hpull_xsec_normncoth","hpull_xsec_normncoth",50,-5,5);
  TH1D * hpull_xsec_normnuetonumu = new TH1D("hpull_xsec_normnuetonumu","hpull_xsec_normnuetonumu",50,-5,5);
  TH1D * hpull_skdet_ccqerc0 = new TH1D("hpull_skdet_ccqerc0","hpull_skdet_ccqerc0",50,-5,5);
  TH1D * hpull_skdet_ccqerc1 = new TH1D("hpull_skdet_ccqerc1","hpull_skdet_ccqerc1",50,-5,5);
  TH1D * hpull_skdet_ccqerc2 = new TH1D("hpull_skdet_ccqerc2","hpull_skdet_ccqerc2",50,-5,5);
  TH1D * hpull_skdet_ccqeoth = new TH1D("hpull_skdet_ccqeoth","hpull_skdet_ccqeoth",50,-5,5);
  TH1D * hpull_skdet_nuecc = new TH1D("hpull_skdet_nuecc","hpull_skdet_nuecc",50,-5,5);
  TH1D * hpull_fsi_000 = new TH1D("hpull_fsi_000","hpull_fsi_000",50,-5,5);
  TH1D * hpull_fsi_001 = new TH1D("hpull_fsi_001","hpull_fsi_001",50,-5,5);
  TH1D * hpull_fsi_002 = new TH1D("hpull_fsi_002","hpull_fsi_002",50,-5,5);
  TH1D * hpull_fsi_003 = new TH1D("hpull_fsi_003","hpull_fsi_003",50,-5,5);
  TH1D * hpull_fsi_004 = new TH1D("hpull_fsi_004","hpull_fsi_004",50,-5,5);
  TH1D * hpull_fsi_005 = new TH1D("hpull_fsi_005","hpull_fsi_005",50,-5,5);

  cout << ch.GetEntries() << " entries in TChain" << endl;
  if(ch.GetEntries() == 0)
    return;

  ch.SetBranchStatus("*",1);
  // Loop over events & fill histograms
  for(int i = 0; i < ch.GetEntries(); i++) {
    ch.GetEntry(i);
    //ignore failed fits!
    h_fitstat      ->Fill(fitstatus);
    h_migradfitstat->Fill(migradfitstatus);
    h_hessefitstat ->Fill(hessefitstatus);
    h_fittime      ->Fill(fittime);
    //minuit
    if(fittermode == 0) {
      if(migradfitstatus == 4) continue; //migrad didn't converge
      if(fitstatus != 3) continue; //error matrix inaccurate
    }
    //minuit2
    //else if(fittermode == 1) {
    //if(fitstatus == 9)
    //continue; //
    //}
  
    /*
    if(bestfit_mixing23 > 0.997)
      continue; 
    if(TMath::Abs(bestfit_mixing23 - 1.) < 1E-3)
      continue; 
    */

    if(err_mode == 1) {
      err_mixing12 = errmigrad_mixing12;
      err_mixing13 = errmigrad_mixing13;
      err_mixing23 = errmigrad_mixing23;
      err_deltacp = errmigrad_deltacp;
      err_dm12sq = errmigrad_dm12sq;
      err_dm23sq = errmigrad_dm23sq;
      err_banff_00 = errmigrad_banff_00;
      err_banff_01 = errmigrad_banff_01;
      err_banff_02 = errmigrad_banff_02;
      err_banff_03 = errmigrad_banff_03;
      err_banff_04 = errmigrad_banff_04;
      err_banff_05 = errmigrad_banff_05;
      err_banff_06 = errmigrad_banff_06;
      err_banff_07 = errmigrad_banff_07;
      err_banff_08 = errmigrad_banff_08;
      err_banff_09 = errmigrad_banff_09;
      err_banff_10 = errmigrad_banff_10;
      err_banff_11 = errmigrad_banff_11;
      err_banff_12 = errmigrad_banff_12;
      err_banff_13 = errmigrad_banff_13;
      err_banff_14 = errmigrad_banff_14;
      err_banff_15 = errmigrad_banff_15;
      err_banff_maqe = errmigrad_banff_maqe;
      err_banff_mares = errmigrad_banff_mares;
      err_banff_ccqeE0 = errmigrad_banff_ccqeE0;
      err_banff_ccqeE1 = errmigrad_banff_ccqeE1;
      err_banff_ccqeE2 = errmigrad_banff_ccqeE2;
      err_banff_cc1piE0 = errmigrad_banff_cc1piE0;
      err_banff_cc1piE1 = errmigrad_banff_cc1piE1;
      err_xsec_pf = errmigrad_xsec_pf;
      err_xsec_sf = errmigrad_xsec_sf;
      err_xsec_binde = errmigrad_xsec_binde;
      err_xsec_wshape = errmigrad_xsec_wshape;
      err_xsec_ccothshape = errmigrad_xsec_ccothshape;
      err_xsec_cc1pishape = errmigrad_xsec_cc1pishape;
      err_xsec_pilessdelta = errmigrad_xsec_pilessdelta;
      err_xsec_normcccoh = errmigrad_xsec_normcccoh;
      err_xsec_normnc1piC = errmigrad_xsec_normnc1piC;
      err_xsec_normnccoh = errmigrad_xsec_normnccoh;
      err_xsec_normncoth = errmigrad_xsec_normncoth;
      err_xsec_normnuetonumu = errmigrad_xsec_normnuetonumu;
      err_skdet_energyscale = errmigrad_skdet_energyscale;
      err_skdet_ccqerc0 = errmigrad_skdet_ccqerc0;
      err_skdet_ccqerc1 = errmigrad_skdet_ccqerc1;
      err_skdet_ccqerc2 = errmigrad_skdet_ccqerc2;
      err_skdet_ccqeoth = errmigrad_skdet_ccqeoth;
      err_skdet_ccnqe = errmigrad_skdet_ccnqe;
      err_skdet_nc = errmigrad_skdet_nc;
      err_skdet_nuecc = errmigrad_skdet_nuecc;
      err_fsi_000 = errmigrad_fsi_000;
      err_fsi_001 = errmigrad_fsi_001;
      err_fsi_002 = errmigrad_fsi_002;
      err_fsi_003 = errmigrad_fsi_003;
      err_fsi_004 = errmigrad_fsi_004;
      err_fsi_005 = errmigrad_fsi_005;
    }//errmode == 1
    else if(err_mode == 2) {
      err_mixing12 = errhesse_mixing12;
      err_mixing13 = errhesse_mixing13;
      err_mixing23 = errhesse_mixing23;
      err_deltacp = errhesse_deltacp;
      err_dm12sq = errhesse_dm12sq;
      err_dm23sq = errhesse_dm23sq;
      err_banff_00 = errhesse_banff_00;
      err_banff_01 = errhesse_banff_01;
      err_banff_02 = errhesse_banff_02;
      err_banff_03 = errhesse_banff_03;
      err_banff_04 = errhesse_banff_04;
      err_banff_05 = errhesse_banff_05;
      err_banff_06 = errhesse_banff_06;
      err_banff_07 = errhesse_banff_07;
      err_banff_08 = errhesse_banff_08;
      err_banff_09 = errhesse_banff_09;
      err_banff_10 = errhesse_banff_10;
      err_banff_11 = errhesse_banff_11;
      err_banff_12 = errhesse_banff_12;
      err_banff_13 = errhesse_banff_13;
      err_banff_14 = errhesse_banff_14;
      err_banff_15 = errhesse_banff_15;
      err_banff_maqe = errhesse_banff_maqe;
      err_banff_mares = errhesse_banff_mares;
      err_banff_ccqeE0 = errhesse_banff_ccqeE0;
      err_banff_ccqeE1 = errhesse_banff_ccqeE1;
      err_banff_ccqeE2 = errhesse_banff_ccqeE2;
      err_banff_cc1piE0 = errhesse_banff_cc1piE0;
      err_banff_cc1piE1 = errhesse_banff_cc1piE1;
      err_xsec_pf = errhesse_xsec_pf;
      err_xsec_sf = errhesse_xsec_sf;
      err_xsec_binde = errhesse_xsec_binde;
      err_xsec_wshape = errhesse_xsec_wshape;
      err_xsec_ccothshape = errhesse_xsec_ccothshape;
      err_xsec_cc1pishape = errhesse_xsec_cc1pishape;
      err_xsec_pilessdelta = errhesse_xsec_pilessdelta;
      err_xsec_normcccoh = errhesse_xsec_normcccoh;
      err_xsec_normnc1piC = errhesse_xsec_normnc1piC;
      err_xsec_normnccoh = errhesse_xsec_normnccoh;
      err_xsec_normncoth = errhesse_xsec_normncoth;
      err_xsec_normnuetonumu = errhesse_xsec_normnuetonumu;
      err_skdet_energyscale = errhesse_skdet_energyscale;
      err_skdet_ccqerc0 = errhesse_skdet_ccqerc0;
      err_skdet_ccqerc1 = errhesse_skdet_ccqerc1;
      err_skdet_ccqerc2 = errhesse_skdet_ccqerc2;
      err_skdet_ccqeoth = errhesse_skdet_ccqeoth;
      err_skdet_ccnqe = errhesse_skdet_ccnqe;
      err_skdet_nc = errhesse_skdet_nc;
      err_skdet_nuecc = errhesse_skdet_nuecc;
      err_fsi_000 = errhesse_fsi_000;
      err_fsi_001 = errhesse_fsi_001;
      err_fsi_002 = errhesse_fsi_002;
      err_fsi_003 = errhesse_fsi_003;
      err_fsi_004 = errhesse_fsi_004;
      err_fsi_005 = errhesse_fsi_005;
    }//errmode == 2
    else if(err_mode == 3) {
      if(bestfit_mixing12 <= input_mixing12) err_mixing12 = errminosplus_mixing12; else err_mixing12 = -errminosminus_mixing12;
      if(bestfit_mixing13 <= input_mixing13) err_mixing13 = errminosplus_mixing13; else err_mixing13 = -errminosminus_mixing13;
      if(bestfit_mixing23 <= input_mixing23) err_mixing23 = errminosplus_mixing23; else err_mixing23 = -errminosminus_mixing23;
      if(bestfit_deltacp <= input_deltacp) err_deltacp = errminosplus_deltacp; else err_deltacp = -errminosminus_deltacp;
      if(bestfit_dm12sq <= input_dm12sq) err_dm12sq = errminosplus_dm12sq; else err_dm12sq = -errminosminus_dm12sq;
      if(bestfit_dm23sq <= input_dm23sq) err_dm23sq = errminosplus_dm23sq; else err_dm23sq = -errminosminus_dm23sq;
      if(bestfit_banff_00 <= input_banff_00) err_banff_00 = errminosplus_banff_00; else err_banff_00 = errminosminus_banff_00;
      if(bestfit_banff_01 <= input_banff_01) err_banff_01 = errminosplus_banff_01; else err_banff_01 = errminosminus_banff_01;
      if(bestfit_banff_02 <= input_banff_02) err_banff_02 = errminosplus_banff_02; else err_banff_02 = errminosminus_banff_02;
      if(bestfit_banff_03 <= input_banff_03) err_banff_03 = errminosplus_banff_03; else err_banff_03 = errminosminus_banff_03;
      if(bestfit_banff_04 <= input_banff_04) err_banff_04 = errminosplus_banff_04; else err_banff_04 = errminosminus_banff_04;
      if(bestfit_banff_05 <= input_banff_05) err_banff_05 = errminosplus_banff_05; else err_banff_05 = errminosminus_banff_05;
      if(bestfit_banff_06 <= input_banff_06) err_banff_06 = errminosplus_banff_06; else err_banff_06 = errminosminus_banff_06;
      if(bestfit_banff_07 <= input_banff_07) err_banff_07 = errminosplus_banff_07; else err_banff_07 = errminosminus_banff_07;
      if(bestfit_banff_08 <= input_banff_08) err_banff_08 = errminosplus_banff_08; else err_banff_08 = errminosminus_banff_08;
      if(bestfit_banff_09 <= input_banff_09) err_banff_09 = errminosplus_banff_09; else err_banff_09 = errminosminus_banff_09;
      if(bestfit_banff_10 <= input_banff_10) err_banff_10 = errminosplus_banff_10; else err_banff_10 = errminosminus_banff_10;
      if(bestfit_banff_11 <= input_banff_11) err_banff_11 = errminosplus_banff_11; else err_banff_11 = errminosminus_banff_11;
      if(bestfit_banff_12 <= input_banff_12) err_banff_12 = errminosplus_banff_12; else err_banff_12 = errminosminus_banff_12;
      if(bestfit_banff_13 <= input_banff_13) err_banff_13 = errminosplus_banff_13; else err_banff_13 = errminosminus_banff_13;
      if(bestfit_banff_14 <= input_banff_14) err_banff_14 = errminosplus_banff_14; else err_banff_14 = errminosminus_banff_14;
      if(bestfit_banff_15 <= input_banff_15) err_banff_15 = errminosplus_banff_15; else err_banff_15 = errminosminus_banff_15;
      if(bestfit_banff_maqe <= input_banff_maqe) err_banff_maqe = errminosplus_banff_maqe; else err_banff_maqe = errminosminus_banff_maqe;
      if(bestfit_banff_mares <= input_banff_mares) err_banff_mares = errminosplus_banff_mares; else err_banff_mares = errminosminus_banff_mares;
      if(bestfit_banff_ccqeE0 <= input_banff_ccqeE0) err_banff_ccqeE0 = errminosplus_banff_ccqeE0; else err_banff_ccqeE0 = errminosminus_banff_ccqeE0;
      if(bestfit_banff_ccqeE1 <= input_banff_ccqeE1) err_banff_ccqeE1 = errminosplus_banff_ccqeE1; else err_banff_ccqeE1 = errminosminus_banff_ccqeE1;
      if(bestfit_banff_ccqeE2 <= input_banff_ccqeE2) err_banff_ccqeE2 = errminosplus_banff_ccqeE2; else err_banff_ccqeE2 = errminosminus_banff_ccqeE2;
      if(bestfit_banff_cc1piE0 <= input_banff_cc1piE0) err_banff_cc1piE0 = errminosplus_banff_cc1piE0; else err_banff_cc1piE0 = errminosminus_banff_cc1piE0;
      if(bestfit_banff_cc1piE1 <= input_banff_cc1piE1) err_banff_cc1piE1 = errminosplus_banff_cc1piE1; else err_banff_cc1piE1 = errminosminus_banff_cc1piE1;
      if(bestfit_xsec_pf <= input_xsec_pf) err_xsec_pf = errminosplus_xsec_pf; else err_xsec_pf = errminosminus_xsec_pf;
      if(bestfit_xsec_sf <= input_xsec_sf) err_xsec_sf = errminosplus_xsec_sf; else err_xsec_sf = errminosminus_xsec_sf;
      if(bestfit_xsec_binde <= input_xsec_binde) err_xsec_binde = errminosplus_xsec_binde; else err_xsec_binde = errminosminus_xsec_binde;
      if(bestfit_xsec_wshape <= input_xsec_wshape) err_xsec_wshape = errminosplus_xsec_wshape; else err_xsec_wshape = errminosminus_xsec_wshape;
      if(bestfit_xsec_ccothshape <= input_xsec_ccothshape) err_xsec_ccothshape = errminosplus_xsec_ccothshape; else err_xsec_ccothshape = errminosminus_xsec_ccothshape;
      if(bestfit_xsec_cc1pishape <= input_xsec_cc1pishape) err_xsec_cc1pishape = errminosplus_xsec_cc1pishape; else err_xsec_cc1pishape = errminosminus_xsec_cc1pishape;
      if(bestfit_xsec_pilessdelta <= input_xsec_pilessdelta) err_xsec_pilessdelta = errminosplus_xsec_pilessdelta; else err_xsec_pilessdelta = errminosminus_xsec_pilessdelta;
      if(bestfit_xsec_normcccoh <= input_xsec_normcccoh) err_xsec_normcccoh = errminosplus_xsec_normcccoh; else err_xsec_normcccoh = errminosminus_xsec_normcccoh;
      if(bestfit_xsec_normnc1piC <= input_xsec_normnc1piC) err_xsec_normnc1piC = errminosplus_xsec_normnc1piC; else err_xsec_normnc1piC = errminosminus_xsec_normnc1piC;
      if(bestfit_xsec_normnccoh <= input_xsec_normnccoh) err_xsec_normnccoh = errminosplus_xsec_normnccoh; else err_xsec_normnccoh = errminosminus_xsec_normnccoh;
      if(bestfit_xsec_normncoth <= input_xsec_normncoth) err_xsec_normncoth = errminosplus_xsec_normncoth; else err_xsec_normncoth = errminosminus_xsec_normncoth;
      if(bestfit_xsec_normnuetonumu <= input_xsec_normnuetonumu) err_xsec_normnuetonumu = errminosplus_xsec_normnuetonumu; else err_xsec_normnuetonumu = errminosminus_xsec_normnuetonumu;
      if(bestfit_skdet_energyscale <= input_skdet_energyscale) err_skdet_energyscale = errminosplus_skdet_energyscale; else err_skdet_energyscale = errminosminus_skdet_energyscale;
      if(bestfit_skdet_ccqerc0 <= input_skdet_ccqerc0) err_skdet_ccqerc0 = errminosplus_skdet_ccqerc0; else err_skdet_ccqerc0 = errminosminus_skdet_ccqerc0;
      if(bestfit_skdet_ccqerc1 <= input_skdet_ccqerc1) err_skdet_ccqerc1 = errminosplus_skdet_ccqerc1; else err_skdet_ccqerc1 = errminosminus_skdet_ccqerc1;
      if(bestfit_skdet_ccqerc2 <= input_skdet_ccqerc2) err_skdet_ccqerc2 = errminosplus_skdet_ccqerc2; else err_skdet_ccqerc2 = errminosminus_skdet_ccqerc2;
      if(bestfit_skdet_ccqeoth <= input_skdet_ccqeoth) err_skdet_ccqeoth = errminosplus_skdet_ccqeoth; else err_skdet_ccqeoth = errminosminus_skdet_ccqeoth;
      if(bestfit_skdet_ccnqe <= input_skdet_ccnqe) err_skdet_ccnqe = errminosplus_skdet_ccnqe; else err_skdet_ccnqe = errminosminus_skdet_ccnqe;
      if(bestfit_skdet_nc <= input_skdet_nc) err_skdet_nc = errminosplus_skdet_nc; else err_skdet_nc = errminosminus_skdet_nc;
      if(bestfit_skdet_nuecc <= input_skdet_nuecc) err_skdet_nuecc = errminosplus_skdet_nuecc; else err_skdet_nuecc = errminosminus_skdet_nuecc;
      if(bestfit_fsi_000 <= input_fsi_000) err_fsi_000 = errminosplus_fsi_000; else err_fsi_000 = errminosminus_fsi_000;
      if(bestfit_fsi_001 <= input_fsi_001) err_fsi_001 = errminosplus_fsi_001; else err_fsi_001 = errminosminus_fsi_001;
      if(bestfit_fsi_002 <= input_fsi_002) err_fsi_002 = errminosplus_fsi_002; else err_fsi_002 = errminosminus_fsi_002;
      if(bestfit_fsi_003 <= input_fsi_003) err_fsi_003 = errminosplus_fsi_003; else err_fsi_003 = errminosminus_fsi_003;
      if(bestfit_fsi_004 <= input_fsi_004) err_fsi_004 = errminosplus_fsi_004; else err_fsi_004 = errminosminus_fsi_004;
      if(bestfit_fsi_005 <= input_fsi_005) err_fsi_005 = errminosplus_fsi_005; else err_fsi_005 = errminosminus_fsi_005;
    }//errmode == 3

    h_nsk_nom->Fill(nsk_nom);
    h_nsk_bf->Fill(nsk_bf);
    h_dchisq->Fill(true_chisq - bestfit_chisq);
    h_mixing23->Fill((bestfit_mixing23));
    h_dm23sq->Fill((bestfit_dm23sq));
    //pull histograms
    hpull_mixing23->Fill((bestfit_mixing23 - input_mixing23) / err_mixing23);
    hpull_dm23sq->Fill((bestfit_dm23sq - input_dm23sq) / err_dm23sq);
    //error histograms
    herr_mixing23->Fill(err_mixing23);
    herr_dm23sq->Fill(err_dm23sq);
    if(fit_systematics_mode >= 1) {
      h_xsec_sf->Fill((bestfit_xsec_sf - input_xsec_sf));
      h_xsec_binde->Fill((bestfit_xsec_binde - input_xsec_binde));
      h_xsec_pilessdelta->Fill((bestfit_xsec_pilessdelta - input_xsec_pilessdelta));
      h_skdet_ccnqe->Fill((bestfit_skdet_ccnqe - input_skdet_ccnqe));
      h_skdet_nc->Fill((bestfit_skdet_nc - input_skdet_nc));
      h_skdet_energyscale->Fill((bestfit_skdet_energyscale - input_skdet_energyscale));
      //pull histograms
      if(err_xsec_sf < preerr_xsec_sf) hpull_xsec_sf->Fill((bestfit_xsec_sf - input_xsec_sf) / sqrt(pow(preerr_xsec_sf,2) - pow(err_xsec_sf,2)));
      if(err_xsec_binde < preerr_xsec_binde) hpull_xsec_binde->Fill((bestfit_xsec_binde - input_xsec_binde) / sqrt(pow(preerr_xsec_binde,2) - pow(err_xsec_binde,2)));
      if(err_xsec_pilessdelta < preerr_xsec_pilessdelta) hpull_xsec_pilessdelta->Fill((bestfit_xsec_pilessdelta - input_xsec_pilessdelta) / sqrt(pow(preerr_xsec_pilessdelta,2) - pow(err_xsec_pilessdelta,2)));
      if(err_skdet_ccnqe < preerr_skdet_ccnqe) hpull_skdet_ccnqe->Fill((bestfit_skdet_ccnqe - input_skdet_ccnqe) / sqrt(pow(preerr_skdet_ccnqe,2) - pow(err_skdet_ccnqe,2)));
      if(err_skdet_nc < preerr_skdet_nc) hpull_skdet_nc->Fill((bestfit_skdet_nc - input_skdet_nc) / sqrt(pow(preerr_skdet_nc,2) - pow(err_skdet_nc,2)));
      if(err_skdet_energyscale < preerr_skdet_energyscale) hpull_skdet_energyscale->Fill((bestfit_skdet_energyscale - input_skdet_energyscale) / sqrt(pow(preerr_skdet_energyscale,2) - pow(err_skdet_energyscale,2)));
      //altpull
      haltpull_xsec_sf->Fill((bestfit_xsec_sf - input_xsec_sf) / err_xsec_sf);
      haltpull_xsec_binde->Fill((bestfit_xsec_binde - input_xsec_binde) / err_xsec_binde);
      haltpull_xsec_pilessdelta->Fill((bestfit_xsec_pilessdelta - input_xsec_pilessdelta) / err_xsec_pilessdelta);
      haltpull_skdet_ccnqe->Fill((bestfit_skdet_ccnqe - input_skdet_ccnqe) / err_skdet_ccnqe);
      haltpull_skdet_nc->Fill((bestfit_skdet_nc - input_skdet_nc) / err_skdet_nc);
      haltpull_skdet_energyscale->Fill((bestfit_skdet_energyscale - input_skdet_energyscale) / err_skdet_energyscale);
      //distribution histograms
      h_banff_00->Fill((bestfit_banff_00 - input_banff_00)/preerr_banff_00);
      h_banff_01->Fill((bestfit_banff_01 - input_banff_01)/preerr_banff_01);
      h_banff_02->Fill((bestfit_banff_02 - input_banff_02)/preerr_banff_02);
      h_banff_03->Fill((bestfit_banff_03 - input_banff_03)/preerr_banff_03);
      h_banff_04->Fill((bestfit_banff_04 - input_banff_04)/preerr_banff_04);
      h_banff_05->Fill((bestfit_banff_05 - input_banff_05)/preerr_banff_05);
      h_banff_06->Fill((bestfit_banff_06 - input_banff_06)/preerr_banff_06);
      h_banff_07->Fill((bestfit_banff_07 - input_banff_07)/preerr_banff_07);
      h_banff_08->Fill((bestfit_banff_08 - input_banff_08)/preerr_banff_08);
      h_banff_09->Fill((bestfit_banff_09 - input_banff_09)/preerr_banff_09);
      h_banff_10->Fill((bestfit_banff_10 - input_banff_10)/preerr_banff_10);
      h_banff_11->Fill((bestfit_banff_11 - input_banff_11)/preerr_banff_11);
      h_banff_12->Fill((bestfit_banff_12 - input_banff_12)/preerr_banff_12);
      h_banff_13->Fill((bestfit_banff_13 - input_banff_13)/preerr_banff_13);
      h_banff_14->Fill((bestfit_banff_14 - input_banff_14)/preerr_banff_14);
      h_banff_15->Fill((bestfit_banff_15 - input_banff_15)/preerr_banff_15);
      h_banff_maqe->Fill((bestfit_banff_maqe - input_banff_maqe)/preerr_banff_maqe);
      h_banff_mares->Fill((bestfit_banff_mares - input_banff_mares)/preerr_banff_mares);
      h_banff_ccqeE0->Fill((bestfit_banff_ccqeE0 - input_banff_ccqeE0)/preerr_banff_ccqeE0);
      h_banff_ccqeE1->Fill((bestfit_banff_ccqeE1 - input_banff_ccqeE1)/preerr_banff_ccqeE1);
      h_banff_ccqeE2->Fill((bestfit_banff_ccqeE2 - input_banff_ccqeE2)/preerr_banff_ccqeE2);
      h_banff_cc1piE0->Fill((bestfit_banff_cc1piE0 - input_banff_cc1piE0)/preerr_banff_cc1piE0);
      h_banff_cc1piE1->Fill((bestfit_banff_cc1piE1 - input_banff_cc1piE1)/preerr_banff_cc1piE1);
      //pull histograms
      if(err_banff_00 < preerr_banff_00) hpull_banff_00->Fill((bestfit_banff_00 - input_banff_00) / sqrt(pow(preerr_banff_00,2) - pow(err_banff_00,2)));
      if(err_banff_01 < preerr_banff_01) hpull_banff_01->Fill((bestfit_banff_01 - input_banff_01) / sqrt(pow(preerr_banff_01,2) - pow(err_banff_01,2)));
      if(err_banff_02 < preerr_banff_02) hpull_banff_02->Fill((bestfit_banff_02 - input_banff_02) / sqrt(pow(preerr_banff_02,2) - pow(err_banff_02,2)));
      if(err_banff_03 < preerr_banff_03) hpull_banff_03->Fill((bestfit_banff_03 - input_banff_03) / sqrt(pow(preerr_banff_03,2) - pow(err_banff_03,2)));
      if(err_banff_04 < preerr_banff_04) hpull_banff_04->Fill((bestfit_banff_04 - input_banff_04) / sqrt(pow(preerr_banff_04,2) - pow(err_banff_04,2)));
      if(err_banff_05 < preerr_banff_05) hpull_banff_05->Fill((bestfit_banff_05 - input_banff_05) / sqrt(pow(preerr_banff_05,2) - pow(err_banff_05,2)));
      if(err_banff_06 < preerr_banff_06) hpull_banff_06->Fill((bestfit_banff_06 - input_banff_06) / sqrt(pow(preerr_banff_06,2) - pow(err_banff_06,2)));
      if(err_banff_07 < preerr_banff_07) hpull_banff_07->Fill((bestfit_banff_07 - input_banff_07) / sqrt(pow(preerr_banff_07,2) - pow(err_banff_07,2)));
      if(err_banff_08 < preerr_banff_08) hpull_banff_08->Fill((bestfit_banff_08 - input_banff_08) / sqrt(pow(preerr_banff_08,2) - pow(err_banff_08,2)));
      if(err_banff_09 < preerr_banff_09) hpull_banff_09->Fill((bestfit_banff_09 - input_banff_09) / sqrt(pow(preerr_banff_09,2) - pow(err_banff_09,2)));
      if(err_banff_10 < preerr_banff_10) hpull_banff_10->Fill((bestfit_banff_10 - input_banff_10) / sqrt(pow(preerr_banff_10,2) - pow(err_banff_10,2)));
      if(err_banff_11 < preerr_banff_11) hpull_banff_11->Fill((bestfit_banff_11 - input_banff_11) / sqrt(pow(preerr_banff_11,2) - pow(err_banff_11,2)));
      if(err_banff_12 < preerr_banff_12) hpull_banff_12->Fill((bestfit_banff_12 - input_banff_12) / sqrt(pow(preerr_banff_12,2) - pow(err_banff_12,2)));
      if(err_banff_13 < preerr_banff_13) hpull_banff_13->Fill((bestfit_banff_13 - input_banff_13) / sqrt(pow(preerr_banff_13,2) - pow(err_banff_13,2)));
      if(err_banff_14 < preerr_banff_14) hpull_banff_14->Fill((bestfit_banff_14 - input_banff_14) / sqrt(pow(preerr_banff_14,2) - pow(err_banff_14,2)));
      if(err_banff_15 < preerr_banff_15) hpull_banff_15->Fill((bestfit_banff_15 - input_banff_15) / sqrt(pow(preerr_banff_15,2) - pow(err_banff_15,2)));
      if(err_banff_maqe < preerr_banff_maqe) hpull_banff_maqe->Fill((bestfit_banff_maqe - input_banff_maqe) / sqrt(pow(preerr_banff_maqe,2) - pow(err_banff_maqe,2)));
      if(err_banff_mares < preerr_banff_mares) hpull_banff_mares->Fill((bestfit_banff_mares - input_banff_mares) / sqrt(pow(preerr_banff_mares,2) - pow(err_banff_mares,2)));
      if(err_banff_ccqeE0 < preerr_banff_ccqeE0) hpull_banff_ccqeE0->Fill((bestfit_banff_ccqeE0 - input_banff_ccqeE0) / sqrt(pow(preerr_banff_ccqeE0,2) - pow(err_banff_ccqeE0,2)));
      if(err_banff_ccqeE1 < preerr_banff_ccqeE1) hpull_banff_ccqeE1->Fill((bestfit_banff_ccqeE1 - input_banff_ccqeE1) / sqrt(pow(preerr_banff_ccqeE1,2) - pow(err_banff_ccqeE1,2)));
      if(err_banff_ccqeE2 < preerr_banff_ccqeE2) hpull_banff_ccqeE2->Fill((bestfit_banff_ccqeE2 - input_banff_ccqeE2) / sqrt(pow(preerr_banff_ccqeE2,2) - pow(err_banff_ccqeE2,2)));
      if(err_banff_cc1piE0 < preerr_banff_cc1piE0) hpull_banff_cc1piE0->Fill((bestfit_banff_cc1piE0 - input_banff_cc1piE0) / sqrt(pow(preerr_banff_cc1piE0,2) - pow(err_banff_cc1piE0,2)));
      if(err_banff_cc1piE1 < preerr_banff_cc1piE1) hpull_banff_cc1piE1->Fill((bestfit_banff_cc1piE1 - input_banff_cc1piE1) / sqrt(pow(preerr_banff_cc1piE1,2) - pow(err_banff_cc1piE1,2)));
      //alt pull
      haltpull_banff_maqe->Fill((bestfit_banff_maqe - input_banff_maqe) / err_banff_maqe);
      haltpull_banff_ccqeE0->Fill((bestfit_banff_ccqeE0 - input_banff_ccqeE0) / err_banff_ccqeE0);
    }//fit_systematics_mode >= 1
    if(fit_systematics_mode == 2) {
      h_xsec_pf->Fill((bestfit_xsec_pf - input_xsec_pf));
      h_xsec_wshape->Fill((bestfit_xsec_wshape - input_xsec_wshape));
      h_xsec_ccothshape->Fill((bestfit_xsec_ccothshape - input_xsec_ccothshape));
      h_xsec_cc1pishape->Fill((bestfit_xsec_cc1pishape - input_xsec_cc1pishape));
      h_xsec_normcccoh->Fill((bestfit_xsec_normcccoh - input_xsec_normcccoh));
      h_xsec_normnc1piC->Fill((bestfit_xsec_normnc1piC - input_xsec_normnc1piC));
      h_xsec_normnccoh->Fill((bestfit_xsec_normnccoh - input_xsec_normnccoh));
      h_xsec_normncoth->Fill((bestfit_xsec_normncoth - input_xsec_normncoth));
      h_xsec_normnuetonumu->Fill((bestfit_xsec_normnuetonumu - input_xsec_normnuetonumu));
      h_skdet_ccqerc0->Fill((bestfit_skdet_ccqerc0 - input_skdet_ccqerc0));
      h_skdet_ccqerc1->Fill((bestfit_skdet_ccqerc1 - input_skdet_ccqerc1));
      h_skdet_ccqerc2->Fill((bestfit_skdet_ccqerc2 - input_skdet_ccqerc2));
      h_skdet_ccqeoth->Fill((bestfit_skdet_ccqeoth - input_skdet_ccqeoth));
      h_skdet_nuecc->Fill((bestfit_skdet_nuecc - input_skdet_nuecc));
      h_fsi_000->Fill((bestfit_fsi_000 - input_fsi_000));
      h_fsi_001->Fill((bestfit_fsi_001 - input_fsi_001));
      h_fsi_002->Fill((bestfit_fsi_002 - input_fsi_002));
      h_fsi_003->Fill((bestfit_fsi_003 - input_fsi_003));
      h_fsi_004->Fill((bestfit_fsi_004 - input_fsi_004));
      h_fsi_005->Fill((bestfit_fsi_005 - input_fsi_005));
      //pull histograms
      if(err_xsec_pf < preerr_xsec_pf) hpull_xsec_pf->Fill((bestfit_xsec_pf - input_xsec_pf) / sqrt(pow(preerr_xsec_pf,2) - pow(err_xsec_pf,2)));
      if(err_xsec_wshape < preerr_xsec_wshape) hpull_xsec_wshape->Fill((bestfit_xsec_wshape - input_xsec_wshape) / sqrt(pow(preerr_xsec_wshape,2) - pow(err_xsec_wshape,2)));
      if(err_xsec_ccothshape < preerr_xsec_ccothshape) hpull_xsec_ccothshape->Fill((bestfit_xsec_ccothshape - input_xsec_ccothshape) / sqrt(pow(preerr_xsec_ccothshape,2) - pow(err_xsec_ccothshape,2)));
      if(err_xsec_cc1pishape < preerr_xsec_cc1pishape) hpull_xsec_cc1pishape->Fill((bestfit_xsec_cc1pishape - input_xsec_cc1pishape) / sqrt(pow(preerr_xsec_cc1pishape,2) - pow(err_xsec_cc1pishape,2)));
      if(err_xsec_normcccoh < preerr_xsec_normcccoh) hpull_xsec_normcccoh->Fill((bestfit_xsec_normcccoh - input_xsec_normcccoh) / sqrt(pow(preerr_xsec_normcccoh,2) - pow(err_xsec_normcccoh,2)));
      if(err_xsec_normnc1piC < preerr_xsec_normnc1piC) hpull_xsec_normnc1piC->Fill((bestfit_xsec_normnc1piC - input_xsec_normnc1piC) / sqrt(pow(preerr_xsec_normnc1piC,2) - pow(err_xsec_normnc1piC,2)));
      if(err_xsec_normnccoh < preerr_xsec_normnccoh) hpull_xsec_normnccoh->Fill((bestfit_xsec_normnccoh - input_xsec_normnccoh) / sqrt(pow(preerr_xsec_normnccoh,2) - pow(err_xsec_normnccoh,2)));
      if(err_xsec_normncoth < preerr_xsec_normncoth) hpull_xsec_normncoth->Fill((bestfit_xsec_normncoth - input_xsec_normncoth) / sqrt(pow(preerr_xsec_normncoth,2) - pow(err_xsec_normncoth,2)));
      if(err_xsec_normnuetonumu < preerr_xsec_normnuetonumu) hpull_xsec_normnuetonumu->Fill((bestfit_xsec_normnuetonumu - input_xsec_normnuetonumu) / sqrt(pow(preerr_xsec_normnuetonumu,2) - pow(err_xsec_normnuetonumu,2)));
      if(err_skdet_ccqerc0 < preerr_skdet_ccqerc0) hpull_skdet_ccqerc0->Fill((bestfit_skdet_ccqerc0 - input_skdet_ccqerc0) / sqrt(pow(preerr_skdet_ccqerc0,2) - pow(err_skdet_ccqerc0,2)));
      if(err_skdet_ccqerc1 < preerr_skdet_ccqerc1) hpull_skdet_ccqerc1->Fill((bestfit_skdet_ccqerc1 - input_skdet_ccqerc1) / sqrt(pow(preerr_skdet_ccqerc1,2) - pow(err_skdet_ccqerc1,2)));
      if(err_skdet_ccqerc2 < preerr_skdet_ccqerc2) hpull_skdet_ccqerc2->Fill((bestfit_skdet_ccqerc2 - input_skdet_ccqerc2) / sqrt(pow(preerr_skdet_ccqerc2,2) - pow(err_skdet_ccqerc2,2)));
      if(err_skdet_ccqeoth < preerr_skdet_ccqeoth) hpull_skdet_ccqeoth->Fill((bestfit_skdet_ccqeoth - input_skdet_ccqeoth) / sqrt(pow(preerr_skdet_ccqeoth,2) - pow(err_skdet_ccqeoth,2)));
      if(err_skdet_nuecc < preerr_skdet_nuecc) hpull_skdet_nuecc->Fill((bestfit_skdet_nuecc - input_skdet_nuecc) / sqrt(pow(preerr_skdet_nuecc,2) - pow(err_skdet_nuecc,2)));
      if(err_fsi_000 < preerr_fsi_000) hpull_fsi_000->Fill((bestfit_fsi_000 - input_fsi_000) / sqrt(pow(preerr_fsi_000,2) - pow(err_fsi_000,2)));
      if(err_fsi_001 < preerr_fsi_001) hpull_fsi_001->Fill((bestfit_fsi_001 - input_fsi_001) / sqrt(pow(preerr_fsi_001,2) - pow(err_fsi_001,2)));
      if(err_fsi_002 < preerr_fsi_002) hpull_fsi_002->Fill((bestfit_fsi_002 - input_fsi_002) / sqrt(pow(preerr_fsi_002,2) - pow(err_fsi_002,2)));
      if(err_fsi_003 < preerr_fsi_003) hpull_fsi_003->Fill((bestfit_fsi_003 - input_fsi_003) / sqrt(pow(preerr_fsi_003,2) - pow(err_fsi_003,2)));
      if(err_fsi_004 < preerr_fsi_004) hpull_fsi_004->Fill((bestfit_fsi_004 - input_fsi_004) / sqrt(pow(preerr_fsi_004,2) - pow(err_fsi_004,2)));
      if(err_fsi_005 < preerr_fsi_005) hpull_fsi_005->Fill((bestfit_fsi_005 - input_fsi_005) / sqrt(pow(preerr_fsi_005,2) - pow(err_fsi_005,2)));
    }//fit_systematic_mode == 2
  }//end loop over chain

  // Draw
  new TCanvas(); h_fitstat->Draw();
  new TCanvas(); h_migradfitstat->Draw();
  new TCanvas(); h_hessefitstat->Draw();
  new TCanvas(); h_fittime->Draw();
  new TCanvas(); h_nsk_nom->Draw();
  new TCanvas(); h_nsk_bf->Draw();
  new TCanvas(); h_dchisq->Draw();
  TCanvas *cosc = new TCanvas(); cosc->Divide(1,2);
  h_mixing23->GetXaxis()->SetRangeUser(0.8,1.0);
  h_dm23sq->GetXaxis()->SetRangeUser(0.0018,0.0030);
  cosc->cd(1); h_mixing23->Draw();
  cosc->cd(2); h_dm23sq->Draw();
  cosc->SaveAs(Form("%f_osc.pdf", input_mixing23));

  if(fit_systematics_mode >= 1) {
    TCanvas *cbanff = new TCanvas(); cbanff->Divide(1,2);
    cbanff->cd(1); h_banff_maqe->Draw();
    cbanff->cd(2); h_banff_ccqeE0->Draw();
    TCanvas *cxsec = new TCanvas(); cxsec->Divide(1,2);
    cxsec->cd(1); h_xsec_sf->Draw();
    cxsec->cd(2); h_xsec_pilessdelta->Draw();
    TCanvas *cextra = new TCanvas(); cextra->Divide(1,2);
    cextra->cd(1); h_xsec_binde->Draw();
    cextra->cd(2); h_skdet_energyscale->Draw();
    TCanvas *cskdet = new TCanvas(); cskdet->Divide(1,2);
    cskdet->cd(1); h_skdet_ccnqe->Draw();
    cskdet->cd(2); h_skdet_nc->Draw();
    
    cbanff->SaveAs(Form("%f_banff.pdf", input_mixing23));
    cxsec->SaveAs(Form("%f_xsec.pdf", input_mixing23));
    cextra->SaveAs(Form("%f_extra.pdf", input_mixing23));
    cskdet->SaveAs(Form("%f_skdet.pdf", input_mixing23));
    /*
    new TCanvas(); h_banff_00->Draw();
    new TCanvas(); h_banff_01->Draw();
    new TCanvas(); h_banff_02->Draw();
    new TCanvas(); h_banff_03->Draw();
    new TCanvas(); h_banff_04->Draw();
    new TCanvas(); h_banff_05->Draw();
    new TCanvas(); h_banff_06->Draw();
    new TCanvas(); h_banff_07->Draw();
    new TCanvas(); h_banff_08->Draw();
    new TCanvas(); h_banff_09->Draw();
    new TCanvas(); h_banff_10->Draw();
    new TCanvas(); h_banff_11->Draw();
    new TCanvas(); h_banff_12->Draw();
    new TCanvas(); h_banff_13->Draw();
    new TCanvas(); h_banff_14->Draw();
    new TCanvas(); h_banff_15->Draw();
    new TCanvas(); h_banff_maqe->Draw();
    new TCanvas(); h_banff_mares->Draw();
    new TCanvas(); h_banff_ccqeE0->Draw();
    new TCanvas(); h_banff_ccqeE1->Draw();
    new TCanvas(); h_banff_ccqeE2->Draw();
    new TCanvas(); h_banff_cc1piE0->Draw();
    new TCanvas(); h_banff_cc1piE1->Draw();
    */
  }
  if(fit_systematics_mode == 2) {
    new TCanvas(); h_xsec_pf->Draw();
    new TCanvas(); h_xsec_sf->Draw();
    new TCanvas(); h_xsec_wshape->Draw();
    new TCanvas(); h_xsec_ccothshape->Draw();
    new TCanvas(); h_xsec_cc1pishape->Draw();
    new TCanvas(); h_xsec_pilessdelta->Draw();
    new TCanvas(); h_xsec_normcccoh->Draw();
    new TCanvas(); h_xsec_normnc1piC->Draw();
    new TCanvas(); h_xsec_normnccoh->Draw();
    new TCanvas(); h_xsec_normncoth->Draw();
    new TCanvas(); h_xsec_normnuetonumu->Draw();
    new TCanvas(); h_skdet_energyscale->Draw();
    new TCanvas(); h_skdet_ccqerc0->Draw();
    new TCanvas(); h_skdet_ccqerc1->Draw();
    new TCanvas(); h_skdet_ccqerc2->Draw();
    new TCanvas(); h_skdet_ccqeoth->Draw();
    new TCanvas(); h_skdet_ccnqe->Draw();
    new TCanvas(); h_skdet_nc->Draw();
    new TCanvas(); h_skdet_nuecc->Draw();
    new TCanvas(); h_fsi_000->Draw();
    new TCanvas(); h_fsi_001->Draw();
    new TCanvas(); h_fsi_002->Draw();
    new TCanvas(); h_fsi_003->Draw();
    new TCanvas(); h_fsi_004->Draw();
    new TCanvas(); h_fsi_005->Draw();
  }


  // pull histograms

  TH1D * meansummary = new TH1D("mean","mean",100,0,100);
  TH1D * sigmasummary = new TH1D("sigma","sigma",100,0,100);
  int ifit = 1;
  TF1 *fit = 0;
  gStyle->SetOptFit();
  

  TCanvas *coscpull = new TCanvas(); coscpull->Divide(1,2);
  coscpull->cd(1); hpull_mixing23->Draw(); hpull_mixing23->Fit("gaus");
  fit = hpull_mixing23->GetFunction("gaus");
  if(fit) {
    meansummary->GetXaxis()->SetBinLabel(ifit,"mixing23"); sigmasummary->GetXaxis()->SetBinLabel(ifit,"mixing23");
    meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
    sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2));
    ifit++;
  }
  coscpull->cd(2); hpull_dm23sq->Draw(); hpull_dm23sq->Fit("gaus");
  fit = hpull_dm23sq->GetFunction("gaus");
  if(fit) {
    meansummary->GetXaxis()->SetBinLabel(ifit,"dm32sq"); sigmasummary->GetXaxis()->SetBinLabel(ifit,"dm32sq");
    meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
    sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2));
    ifit++;
  }
  coscpull->SaveAs(Form("%f_osc_pull.pdf", input_mixing23));

  TCanvas *coscerr = new TCanvas(); coscerr->Divide(1,2);
  coscerr->cd(1); herr_mixing23->Draw();
  coscerr->cd(2); herr_dm23sq->Draw();

  if(fit_systematics_mode >= 1) {

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
    cxsec->cd(1); hpull_xsec_sf->Draw(); hpull_xsec_sf->Fit("gaus");
    if(fit) {
      fit = hpull_xsec_sf->GetFunction("gaus");
      meansummary->GetXaxis()->SetBinLabel(ifit,"SF"); sigmasummary->GetXaxis()->SetBinLabel(ifit,"SF");
      meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
      sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2));
      ifit++;
    }
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

    cbanff->SaveAs(Form("%f_banff_pull.pdf", input_mixing23));
    cxsec->SaveAs(Form("%f_xsec_pull.pdf", input_mixing23));
    cextra->SaveAs(Form("%f_extra_pull.pdf", input_mixing23));
    cskdet->SaveAs(Form("%f_skdet_pull.pdf", input_mixing23));

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
    caltxsec->cd(1); haltpull_xsec_sf->Draw(); haltpull_xsec_sf->Fit("gaus");
    fit = haltpull_xsec_sf->GetFunction("gaus");
    if(fit) {
      meansummary->GetXaxis()->SetBinLabel(ifit,"SF"); sigmasummary->GetXaxis()->SetBinLabel(ifit,"SF");
      meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
      sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2));
      ifit++;
    }
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

    caltbanff->SaveAs(Form("%f_banff_altpull.pdf", input_mixing23));
    caltxsec->SaveAs(Form("%f_xsec_altpull.pdf", input_mixing23));
    caltextra->SaveAs(Form("%f_extra_altpull.pdf", input_mixing23));
    caltskdet->SaveAs(Form("%f_skdet_altpull.pdf", input_mixing23));

    /*
    new TCanvas(); hpull_banff_00->Draw(); hpull_banff_00->Fit("gaus"); fit = hpull_banff_00->GetFunction("gaus");
    meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
    sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2)); ifit++;
    new TCanvas(); hpull_banff_01->Draw(); hpull_banff_01->Fit("gaus"); fit = hpull_banff_01->GetFunction("gaus");
    meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
    sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2)); ifit++;
    new TCanvas(); hpull_banff_02->Draw(); hpull_banff_02->Fit("gaus"); fit = hpull_banff_02->GetFunction("gaus");
    meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
    sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2)); ifit++;
    new TCanvas(); hpull_banff_03->Draw(); hpull_banff_03->Fit("gaus"); fit = hpull_banff_03->GetFunction("gaus");
    meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
    sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2)); ifit++;
    new TCanvas(); hpull_banff_04->Draw(); hpull_banff_04->Fit("gaus"); fit = hpull_banff_04->GetFunction("gaus");
    meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
    sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2)); ifit++;
    new TCanvas(); hpull_banff_05->Draw(); hpull_banff_05->Fit("gaus"); fit = hpull_banff_05->GetFunction("gaus");
    meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
    sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2)); ifit++;
    new TCanvas(); hpull_banff_06->Draw(); hpull_banff_06->Fit("gaus"); fit = hpull_banff_06->GetFunction("gaus");
    meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
    sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2)); ifit++;
    new TCanvas(); hpull_banff_07->Draw(); hpull_banff_07->Fit("gaus"); fit = hpull_banff_07->GetFunction("gaus");
    meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
    sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2)); ifit++;
    new TCanvas(); hpull_banff_08->Draw(); hpull_banff_08->Fit("gaus"); fit = hpull_banff_08->GetFunction("gaus");
    meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
    sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2)); ifit++;
    new TCanvas(); hpull_banff_09->Draw(); hpull_banff_09->Fit("gaus"); fit = hpull_banff_09->GetFunction("gaus");
    meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
    sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2)); ifit++;
    new TCanvas(); hpull_banff_10->Draw(); hpull_banff_10->Fit("gaus"); fit = hpull_banff_10->GetFunction("gaus");
    meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
    sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2)); ifit++;
    new TCanvas(); hpull_banff_11->Draw(); hpull_banff_11->Fit("gaus"); fit = hpull_banff_11->GetFunction("gaus");
    meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
    sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2)); ifit++;
    new TCanvas(); hpull_banff_12->Draw(); hpull_banff_12->Fit("gaus"); fit = hpull_banff_12->GetFunction("gaus");
    meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
    sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2)); ifit++;
    new TCanvas(); hpull_banff_13->Draw(); hpull_banff_13->Fit("gaus"); fit = hpull_banff_13->GetFunction("gaus");
    meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
    sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2)); ifit++;
    new TCanvas(); hpull_banff_14->Draw(); hpull_banff_14->Fit("gaus"); fit = hpull_banff_14->GetFunction("gaus");
    meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
    sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2)); ifit++;
    new TCanvas(); hpull_banff_15->Draw(); hpull_banff_15->Fit("gaus"); fit = hpull_banff_15->GetFunction("gaus");
    meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
    sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2)); ifit++;
    new TCanvas(); hpull_banff_maqe->Draw(); hpull_banff_maqe->Fit("gaus"); fit = hpull_banff_maqe->GetFunction("gaus");
    meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
    sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2)); ifit++;
    new TCanvas(); hpull_banff_mares->Draw(); hpull_banff_mares->Fit("gaus"); fit = hpull_banff_mares->GetFunction("gaus");
    meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
    sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2)); ifit++;
    new TCanvas(); hpull_banff_ccqeE0->Draw(); hpull_banff_ccqeE0->Fit("gaus"); fit = hpull_banff_ccqeE0->GetFunction("gaus");
    meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
    sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2)); ifit++;
    new TCanvas(); hpull_banff_ccqeE1->Draw(); hpull_banff_ccqeE1->Fit("gaus"); fit = hpull_banff_ccqeE1->GetFunction("gaus");
    meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
    sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2)); ifit++;
    new TCanvas(); hpull_banff_ccqeE2->Draw(); hpull_banff_ccqeE2->Fit("gaus"); fit = hpull_banff_ccqeE2->GetFunction("gaus");
    meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
    sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2)); ifit++;
    new TCanvas(); hpull_banff_cc1piE0->Draw(); hpull_banff_cc1piE0->Fit("gaus"); fit = hpull_banff_cc1piE0->GetFunction("gaus");
    meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
    sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2)); ifit++;
    new TCanvas(); hpull_banff_cc1piE1->Draw(); hpull_banff_cc1piE1->Fit("gaus"); fit = hpull_banff_cc1piE1->GetFunction("gaus");
    meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
    sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2)); ifit++;
    */
  }
  if(fit_systematics_mode == 2) {
    new TCanvas(); hpull_xsec_pf->Draw(); hpull_xsec_pf->Fit("gaus"); fit = hpull_xsec_pf->GetFunction("gaus");
    meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
    sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2)); ifit++;
    new TCanvas(); hpull_xsec_wshape->Draw(); hpull_xsec_wshape->Fit("gaus"); fit = hpull_xsec_wshape->GetFunction("gaus");
    meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
    sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2)); ifit++;
    new TCanvas(); hpull_xsec_ccothshape->Draw(); hpull_xsec_ccothshape->Fit("gaus"); fit = hpull_xsec_ccothshape->GetFunction("gaus");
    meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
    sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2)); ifit++;
    new TCanvas(); hpull_xsec_cc1pishape->Draw(); hpull_xsec_cc1pishape->Fit("gaus"); fit = hpull_xsec_cc1pishape->GetFunction("gaus");
    meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
    sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2)); ifit++;
    new TCanvas(); hpull_xsec_normcccoh->Draw(); hpull_xsec_normcccoh->Fit("gaus"); fit = hpull_xsec_normcccoh->GetFunction("gaus");
    meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
    sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2)); ifit++;
    new TCanvas(); hpull_xsec_normnc1piC->Draw(); hpull_xsec_normnc1piC->Fit("gaus"); fit = hpull_xsec_normnc1piC->GetFunction("gaus");
    meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
    sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2)); ifit++;
    new TCanvas(); hpull_xsec_normnccoh->Draw(); hpull_xsec_normnccoh->Fit("gaus"); fit = hpull_xsec_normnccoh->GetFunction("gaus");
    meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
    sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2)); ifit++;
    new TCanvas(); hpull_xsec_normncoth->Draw(); hpull_xsec_normncoth->Fit("gaus"); fit = hpull_xsec_normncoth->GetFunction("gaus");
    meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
    sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2)); ifit++;
    new TCanvas(); hpull_xsec_normnuetonumu->Draw(); hpull_xsec_normnuetonumu->Fit("gaus"); fit = hpull_xsec_normnuetonumu->GetFunction("gaus");
    meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
    sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2)); ifit++;
    new TCanvas(); hpull_skdet_ccqerc0->Draw(); hpull_skdet_ccqerc0->Fit("gaus"); fit = hpull_skdet_ccqerc0->GetFunction("gaus");
    meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
    sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2)); ifit++;
    new TCanvas(); hpull_skdet_ccqerc1->Draw(); hpull_skdet_ccqerc1->Fit("gaus"); fit = hpull_skdet_ccqerc1->GetFunction("gaus");
    meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
    sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2)); ifit++;
    new TCanvas(); hpull_skdet_ccqerc2->Draw(); hpull_skdet_ccqerc2->Fit("gaus"); fit = hpull_skdet_ccqerc2->GetFunction("gaus");
    meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
    sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2)); ifit++;
    new TCanvas(); hpull_skdet_ccqeoth->Draw(); hpull_skdet_ccqeoth->Fit("gaus"); fit = hpull_skdet_ccqeoth->GetFunction("gaus");
    meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
    sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2)); ifit++;
    new TCanvas(); hpull_skdet_nuecc->Draw(); hpull_skdet_nuecc->Fit("gaus"); fit = hpull_skdet_nuecc->GetFunction("gaus");
    meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
    sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2)); ifit++;
    new TCanvas(); hpull_fsi_000->Draw(); hpull_fsi_000->Fit("gaus"); fit = hpull_fsi_000->GetFunction("gaus");
    meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
    sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2)); ifit++;
    new TCanvas(); hpull_fsi_001->Draw(); hpull_fsi_001->Fit("gaus"); fit = hpull_fsi_001->GetFunction("gaus");
    meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
    sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2)); ifit++;
    new TCanvas(); hpull_fsi_002->Draw(); hpull_fsi_002->Fit("gaus"); fit = hpull_fsi_002->GetFunction("gaus");
    meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
    sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2)); ifit++;
    new TCanvas(); hpull_fsi_003->Draw(); hpull_fsi_003->Fit("gaus"); fit = hpull_fsi_003->GetFunction("gaus");
    meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
    sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2)); ifit++;
    new TCanvas(); hpull_fsi_004->Draw(); hpull_fsi_004->Fit("gaus"); fit = hpull_fsi_004->GetFunction("gaus");
    meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
    sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2)); ifit++;
    new TCanvas(); hpull_fsi_005->Draw(); hpull_fsi_005->Fit("gaus"); fit = hpull_fsi_005->GetFunction("gaus");
    meansummary->SetBinContent(ifit, fit->GetParameter(1)); meansummary->SetBinError(ifit, fit->GetParError(1));
    sigmasummary->SetBinContent(ifit, fit->GetParameter(2)); sigmasummary->SetBinError(ifit, fit->GetParError(2)); ifit++;
  }

  if(fit_systematics_mode > 0) {
    TCanvas * csummary = new TCanvas(); csummary->Divide(1,2);
    gStyle->SetOptStat(0);
    if(fit_systematics_mode == 1) {
      meansummary->GetXaxis()->SetRangeUser(0,20);
      sigmasummary->GetXaxis()->SetRangeUser(0,20);
    }
    csummary->cd(1); meansummary->Draw();
    csummary->cd(2); sigmasummary->Draw();
    csummary->SaveAs(Form("%f_summary_pull.pdf", input_mixing23));
    gStyle->SetOptStat(1);
  }
}
