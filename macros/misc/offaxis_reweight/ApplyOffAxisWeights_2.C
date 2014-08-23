#include <fstream.h>
#include <iomanip.h>

using namespace std;
//using namespace valor;

const char *outfilename_spectra_offaxis2p0 = "sk_tuned11bv3.2_11anom_run1-run4_fine_offaxis_2p0.root";
const char *outfilename_spectra_offaxis2p25 = "sk_tuned11bv3.2_11anom_run1-run4_fine_offaxis_2p25.root";

const char *outfilename_spectra_offaxis2p0_rhc = "sk_tuned11bv3.2_11anom_m250ka_fine_2p0.root";
const char *outfilename_spectra_offaxis2p25_rhc = "sk_tuned11bv3.2_11anom_m250ka_fine_2p25.root";

const char *infilename_weights = "flux_offaxis_weights.root";
/*
TFile**		flux_offaxis_weights.root	
 TFile*		flux_offaxis_weights.root	
  OBJ: TH1D	weights_oa_2p0_numode_numub	 : 0 at: 0x7fd2fb8ba0a0
  OBJ: TH1D	weights_oa_2p0_numode_numub	 : 0 at: 0x7fd2fb96f730
  KEY: TH1D	weights_oa_2p0_numode_numu;1	
  KEY: TH1D	weights_oa_2p25_numode_numu;1	
  KEY: TH1D	weights_oa_2p0_numode_numub;1	
  KEY: TH1D	weights_oa_2p25_numode_numub;1	
  KEY: TH1D	weights_oa_2p0_numode_nue;1	
  KEY: TH1D	weights_oa_2p25_numode_nue;1	
  KEY: TH1D	weights_oa_2p0_numode_nueb;1	
  KEY: TH1D	weights_oa_2p25_numode_nueb;1	
  KEY: TH1D	weights_oa_2p0_anumode_numu;1	
  KEY: TH1D	weights_oa_2p25_anumode_numu;1	
  KEY: TH1D	weights_oa_2p0_anumode_numub;1	
  KEY: TH1D	weights_oa_2p25_anumode_numub;1	
  KEY: TH1D	weights_oa_2p0_anumode_nue;1	
  KEY: TH1D	weights_oa_2p25_anumode_nue;1	
  KEY: TH1D	weights_oa_2p0_anumode_nueb;1	
  KEY: TH1D	weights_oa_2p25_anumode_nueb;1	
*/

const char *infilename_spectra = "/Users/davidesgalaberna/Desktop/ETH_phd/T2K/ANA/VALOR/valor/softw/t2xk/trunk/inputs/official/2013/flux/11b_v3.2_run1-4a/sk_tuned11bv3.2_11anom_run1-run4_fine.root";
/*
  TFile**		sk_tuned11bv3.2_11anom_run1-run4_fine.root	
  TFile*		sk_tuned11bv3.2_11anom_run1-run4_fine.root	
  KEY: TH1D	enu_sk_tuned11b_numu_ratio;1	11b tuning v3.2 flux (Run1+2+3+4) / 11a nominal flux --
  KEY: TH1D	enu_sk_tuned11b_numub_ratio;1	11b tuning v3.2 flux (Run1+2+3+4) / 11a nominal flux
  KEY: TH1D	enu_sk_tuned11b_nue_ratio;1	11b tuning v3.2 flux (Run1+2+3+4) / 11a nominal flux
  KEY: TH1D	enu_sk_tuned11b_nueb_ratio;1	11b tuning v3.2 flux (Run1+2+3+4) / 11a nominal flux
  KEY: TH1D	enu_sk_11b_real_numu_ratio;1	real p beam profile flux (Run1+2+3+4) / 11a nominal flux
  KEY: TH1D	enu_sk_11b_real_numub_ratio;1	real p beam profile flux (Run1+2+3+4) / 11a nominal flux
  KEY: TH1D	enu_sk_11b_real_nue_ratio;1	real p beam profile flux (Run1+2+3+4) / 11a nominal flux
  KEY: TH1D	enu_sk_11b_real_nueb_ratio;1	real p beam profile flux (Run1+2+3+4) / 11a nominal flux
  KEY: TH1D	enu_sk_11a_norm_numu;1	11a flux : nominal p beam profile (sigma=0.42 cm)
  KEY: TH1D	enu_sk_11a_norm_numub;1	11a flux : nominal p beam profile (sigma=0.42 cm)
  KEY: TH1D	enu_sk_11a_norm_nue;1	11a flux : nominal p beam profile (sigma=0.42 cm)
  KEY: TH1D	enu_sk_11a_norm_nueb;1	11a flux : nominal p beam profile (sigma=0.42 cm)
  KEY: TH1D	enu_sk_11b_real_numu;1	11b flux : real p beam profile (Run1+2+3+4)
  KEY: TH1D	enu_sk_11b_real_numub;1	11b flux : real p beam profile (Run1+2+3+4)
  KEY: TH1D	enu_sk_11b_real_nue;1	11b flux : real p beam profile (Run1+2+3+4)
  KEY: TH1D	enu_sk_11b_real_nueb;1	11b flux : real p beam profile (Run1+2+3+4)
  KEY: TH1D	enu_sk_tuned11b_numu;1	11b flux : 11b tuning ver.3.2 (Run1+2+3+4)
  KEY: TH1D	enu_sk_tuned11b_numub;1	11b flux : 11b tuning ver.3.2 (Run1+2+3+4)
  KEY: TH1D	enu_sk_tuned11b_nue;1	11b flux : 11b tuning ver.3.2 (Run1+2+3+4)
  KEY: TH1D	enu_sk_tuned11b_nueb;1	11b flux : 11b tuning ver.3.2 (Run1+2+3+4)
*/

const char *infilename_spectra_rhc = "/Users/davidesgalaberna/Desktop/ETH_phd/T2K/ANA/VALOR/valor/softw/t2xk/trunk/inputs/official/2013/flux/11b_v3.2_antinu/sk_tuned11bv3.2_11anom_m250ka_fine.root";

/*
  TFile**		inputs/official/2013/flux/11b_v3.2_antinu/sk_tuned11bv3.2_11anom_m250ka_fine.root	
  TFile*		inputs/official/2013/flux/11b_v3.2_antinu/sk_tuned11bv3.2_11anom_m250ka_fine.root	
  KEY: TH1D	enu_sk_tuned11b_numu_ratio;1	11b tuning v3.2 flux (-250kA) / 11a nominal flux
  KEY: TH1D	enu_sk_tuned11b_numub_ratio;1	11b tuning v3.2 flux (-250kA) / 11a nominal flux
  KEY: TH1D	enu_sk_tuned11b_nue_ratio;1	11b tuning v3.2 flux (-250kA) / 11a nominal flux
  KEY: TH1D	enu_sk_tuned11b_nueb_ratio;1	11b tuning v3.2 flux (-250kA) / 11a nominal flux
  KEY: TH1D	enu_sk_11b_real_numu_ratio;1	11b nominal flux (-250kA) / 11a nominal flux
  KEY: TH1D	enu_sk_11b_real_numub_ratio;1	11b nominal flux (-250kA) / 11a nominal flux
  KEY: TH1D	enu_sk_11b_real_nue_ratio;1	11b nominal flux (-250kA) / 11a nominal flux
  KEY: TH1D	enu_sk_11b_real_nueb_ratio;1	11b nominal flux (-250kA) / 11a nominal flux
  KEY: TH1D	enu_sk_11a_norm_numu;1	11a flux : nominal p beam profile (sigma=0.42 cm)
  KEY: TH1D	enu_sk_11a_norm_numub;1	11a flux : nominal p beam profile (sigma=0.42 cm)
  KEY: TH1D	enu_sk_11a_norm_nue;1	11a flux : nominal p beam profile (sigma=0.42 cm)
  KEY: TH1D	enu_sk_11a_norm_nueb;1	11a flux : nominal p beam profile (sigma=0.42 cm)
  KEY: TH1D	enu_sk_11b_real_numu;1	11b flux : 11b nominal flux (-250kA)
  KEY: TH1D	enu_sk_11b_real_numub;1	11b flux : 11b nominal flux (-250kA)
  KEY: TH1D	enu_sk_11b_real_nue;1	11b flux : 11b nominal flux (-250kA)
  KEY: TH1D	enu_sk_11b_real_nueb;1	11b flux : 11b nominal flux (-250kA)
  KEY: TH1D	enu_sk_tuned11b_numu;1	11b flux : 11b tuning ver.3.2 (-250kA)
  KEY: TH1D	enu_sk_tuned11b_numub;1	11b flux : 11b tuning ver.3.2 (-250kA)
  KEY: TH1D	enu_sk_tuned11b_nue;1	11b flux : 11b tuning ver.3.2 (-250kA)
  KEY: TH1D	enu_sk_tuned11b_nueb;1	11b flux : 11b tuning ver.3.2 (-250kA)
*/

void ApplyOffAxisWeights_2(){

  // Take off-axis weights (2.0 degree)
  TFile *infile_weights = new TFile(infilename_weights,"READ");

  TH1D *hWeight_numu_fhc_2p0 = (TH1D*) infile_weights->Get("weights_oa_2p0_numode_numu");
  TH1D *hWeight_numub_fhc_2p0 = (TH1D*) infile_weights->Get("weights_oa_2p0_numode_numub");
  TH1D *hWeight_nue_fhc_2p0 = (TH1D*) infile_weights->Get("weights_oa_2p0_numode_nue");
  TH1D *hWeight_nueb_fhc_2p0 = (TH1D*) infile_weights->Get("weights_oa_2p0_numode_nueb");

  TH1D *hWeight_numu_rhc_2p0 = (TH1D*) infile_weights->Get("weights_oa_2p0_anumode_numu");
  TH1D *hWeight_numub_rhc_2p0 = (TH1D*) infile_weights->Get("weights_oa_2p0_anumode_numub");
  TH1D *hWeight_nue_rhc_2p0 = (TH1D*) infile_weights->Get("weights_oa_2p0_anumode_nue");
  TH1D *hWeight_nueb_rhc_2p0 = (TH1D*) infile_weights->Get("weights_oa_2p0_anumode_nueb");

  TH1D *hWeight_numu_fhc_2p25 = (TH1D*) infile_weights->Get("weights_oa_2p25_numode_numu");
  TH1D *hWeight_numub_fhc_2p25 = (TH1D*) infile_weights->Get("weights_oa_2p25_numode_numub");
  TH1D *hWeight_nue_fhc_2p25 = (TH1D*) infile_weights->Get("weights_oa_2p25_numode_nue");
  TH1D *hWeight_nueb_fhc_2p25 = (TH1D*) infile_weights->Get("weights_oa_2p25_numode_nueb");

  TH1D *hWeight_numu_rhc_2p25 = (TH1D*) infile_weights->Get("weights_oa_2p25_anumode_numu");
  TH1D *hWeight_numub_rhc_2p25 = (TH1D*) infile_weights->Get("weights_oa_2p25_anumode_numub");
  TH1D *hWeight_nue_rhc_2p25 = (TH1D*) infile_weights->Get("weights_oa_2p25_anumode_nue");
  TH1D *hWeight_nueb_rhc_2p25 = (TH1D*) infile_weights->Get("weights_oa_2p25_anumode_nueb");

  // Take spectra histograms (fhc)
  TFile *infile_spectra = new TFile(infilename_spectra,"READ");  
  TH1D *hNumuSpec = (TH1D*) infile_spectra->Get("enu_sk_tuned11b_numu");
  TH1D *hNumubSpec = (TH1D*) infile_spectra->Get("enu_sk_tuned11b_numub");
  TH1D *hNueSpec = (TH1D*) infile_spectra->Get("enu_sk_tuned11b_nue");
  TH1D *hNuebSpec = (TH1D*) infile_spectra->Get("enu_sk_tuned11b_nueb");

  TH1D *hNumuRatio = (TH1D*) infile_spectra->Get("enu_sk_tuned11b_numu_ratio");
  TH1D *hNumubRatio = (TH1D*) infile_spectra->Get("enu_sk_tuned11b_numub_ratio");
  TH1D *hNueRatio = (TH1D*) infile_spectra->Get("enu_sk_tuned11b_nue_ratio");
  TH1D *hNuebRatio = (TH1D*) infile_spectra->Get("enu_sk_tuned11b_nueb_ratio");
  
  TH1D *hNumuSpec_w2p0 = (TH1D*) hNumuSpec->Clone();
  TH1D *hNumubSpec_w2p0 = (TH1D*) hNumubSpec->Clone();
  TH1D *hNueSpec_w2p0 = (TH1D*) hNueSpec->Clone();
  TH1D *hNuebSpec_w2p0 = (TH1D*) hNuebSpec->Clone();

  TH1D *hNumuRatio_w2p0 = (TH1D*) hNumuRatio->Clone();
  TH1D *hNumubRatio_w2p0 = (TH1D*) hNumubRatio->Clone();
  TH1D *hNueRatio_w2p0 = (TH1D*) hNueRatio->Clone();
  TH1D *hNuebRatio_w2p0 = (TH1D*) hNuebRatio->Clone();

  TH1D *hNumuSpec_w2p25 = (TH1D*) hNumuSpec->Clone();
  TH1D *hNumubSpec_w2p25 = (TH1D*) hNumubSpec->Clone();
  TH1D *hNueSpec_w2p25 = (TH1D*) hNueSpec->Clone();
  TH1D *hNuebSpec_w2p25 = (TH1D*) hNuebSpec->Clone();

  TH1D *hNumuRatio_w2p25 = (TH1D*) hNumuRatio->Clone();
  TH1D *hNumubRatio_w2p25 = (TH1D*) hNumubRatio->Clone();
  TH1D *hNueRatio_w2p25 = (TH1D*) hNueRatio->Clone();
  TH1D *hNuebRatio_w2p25 = (TH1D*) hNuebRatio->Clone();

  // Take spectra histograms (rhc)
  TFile *infile_spectra_rhc = new TFile(infilename_spectra_rhc,"READ");  
  TH1D *hNumuSpec_rhc = (TH1D*) infile_spectra_rhc->Get("enu_sk_tuned11b_numu");
  TH1D *hNumubSpec_rhc = (TH1D*) infile_spectra_rhc->Get("enu_sk_tuned11b_numub");
  TH1D *hNueSpec_rhc = (TH1D*) infile_spectra_rhc->Get("enu_sk_tuned11b_nue");
  TH1D *hNuebSpec_rhc = (TH1D*) infile_spectra_rhc->Get("enu_sk_tuned11b_nueb");

  TH1D *hNumuRatio_rhc = (TH1D*) infile_spectra_rhc->Get("enu_sk_tuned11b_numu_ratio");
  TH1D *hNumubRatio_rhc = (TH1D*) infile_spectra_rhc->Get("enu_sk_tuned11b_numub_ratio");
  TH1D *hNueRatio_rhc = (TH1D*) infile_spectra_rhc->Get("enu_sk_tuned11b_nue_ratio");
  TH1D *hNuebRatio_rhc = (TH1D*) infile_spectra_rhc->Get("enu_sk_tuned11b_nueb_ratio");

  TH1D *hNumuSpec_w2p0_rhc = (TH1D*) hNumuSpec_rhc->Clone();
  TH1D *hNumubSpec_w2p0_rhc = (TH1D*) hNumubSpec_rhc->Clone();
  TH1D *hNueSpecw2p0_rhc_ = (TH1D*) hNueSpec_rhc->Clone();
  TH1D *hNuebSpec_w2p0_rhc = (TH1D*) hNuebSpec_rhc->Clone();

  TH1D *hNumuRatio_w2p0_rhc = (TH1D*) hNumuRatio_rhc->Clone();
  TH1D *hNumubRatio_w2p0_rhc = (TH1D*) hNumubRatio_rhc->Clone();
  TH1D *hNueRatio_w2p0_rhc = (TH1D*) hNueRatio_rhc->Clone();
  TH1D *hNuebRatio_w2p0_rhc = (TH1D*) hNuebRatio_rhc->Clone();

  TH1D *hNumuSpec_w2p25_rhc = (TH1D*) hNumuSpec_rhc->Clone();
  TH1D *hNumubSpec_w2p25_rhc = (TH1D*) hNumubSpec_rhc->Clone();
  TH1D *hNueSpec_w2p25_rhc = (TH1D*) hNueSpec_rhc->Clone();
  TH1D *hNuebSpec_w2p25_rhc = (TH1D*) hNuebSpec_rhc->Clone();
  
  TH1D *hNumuRatio_w2p25_rhc = (TH1D*) hNumuRatio_rhc->Clone();
  TH1D *hNumubRatio_w2p25_rhc = (TH1D*) hNumubRatio_rhc->Clone();
  TH1D *hNueRatio_w2p25_rhc = (TH1D*) hNueRatio_rhc->Clone();
  TH1D *hNuebRatio_w2p25_rhc = (TH1D*) hNuebRatio_rhc->Clone();

  // Weight histograms

  hNumuSpec_w2p0 = (TH1D*) ApplyOffAxisWeights(hNumuSpec,hWeight_numu_fhc_2p0);
  hNumubSpec_w2p0 = (TH1D*) ApplyOffAxisWeights(hNumubSpec,hWeight_numub_fhc_2p0);
  hNueSpec_w2p0 = (TH1D*) ApplyOffAxisWeights(hNueSpec,hWeight_nue_fhc_2p0);
  hNuebSpec_w2p0 = (TH1D*) ApplyOffAxisWeights(hNuebSpec,hWeight_nueb_fhc_2p0);

  hNumuSpec_w2p0_rhc = (TH1D*) ApplyOffAxisWeights(hNumuSpec_rhc,hWeight_numu_rhc_2p0);
  hNumubSpec_w2p0_rhc = (TH1D*) ApplyOffAxisWeights(hNumubSpec_rhc,hWeight_numub_rhc_2p0);
  hNueSpec_w2p0_rhc = (TH1D*) ApplyOffAxisWeights(hNueSpec_rhc,hWeight_nue_rhc_2p0);
  hNuebSpec_w2p0_rhc = (TH1D*) ApplyOffAxisWeights(hNuebSpec_rhc,hWeight_nueb_rhc_2p0);

  hNumuSpec_w2p25 = (TH1D*) ApplyOffAxisWeights(hNumuSpec,hWeight_numu_fhc_2p25);
  hNumubSpec_w2p25 = (TH1D*) ApplyOffAxisWeights(hNumubSpec,hWeight_numub_fhc_2p25);
  hNueSpec_w2p25 = (TH1D*) ApplyOffAxisWeights(hNueSpec,hWeight_nue_fhc_2p25);
  hNuebSpec_w2p25 = (TH1D*) ApplyOffAxisWeights(hNuebSpec,hWeight_nueb_fhc_2p25);

  hNumuSpec_w2p25_rhc = (TH1D*) ApplyOffAxisWeights(hNumuSpec_rhc,hWeight_numu_rhc_2p25);
  hNumubSpec_w2p25_rhc = (TH1D*) ApplyOffAxisWeights(hNumubSpec_rhc,hWeight_numub_rhc_2p25);
  hNueSpec_w2p25_rhc = (TH1D*) ApplyOffAxisWeights(hNueSpec_rhc,hWeight_nue_rhc_2p25);
  hNuebSpec_w2p25_rhc = (TH1D*) ApplyOffAxisWeights(hNuebSpec_rhc,hWeight_nueb_rhc_2p25);

  // Get Ratio weights
    
  hNumuRatio_w2p0 = (TH1D*) ApplyOffAxisWeights(hNumuRatio,hWeight_numu_fhc_2p0);
  hNumubRatio_w2p0 = (TH1D*) ApplyOffAxisWeights(hNumubRatio,hWeight_numub_fhc_2p0);
  hNueRatio_w2p0 = (TH1D*) ApplyOffAxisWeights(hNueRatio,hWeight_nue_fhc_2p0);
  hNuebRatio_w2p0 = (TH1D*) ApplyOffAxisWeights(hNuebRatio,hWeight_nueb_fhc_2p0);

  hNumuRatio_w2p0_rhc = (TH1D*) ApplyOffAxisWeights(hNumuRatio_rhc,hWeight_numu_rhc_2p0);
  hNumubRatio_w2p0_rhc = (TH1D*) ApplyOffAxisWeights(hNumubRatio_rhc,hWeight_numub_rhc_2p0);
  hNueRatio_w2p0_rhc = (TH1D*) ApplyOffAxisWeights(hNueRatio_rhc,hWeight_nue_rhc_2p0);
  hNuebRatio_w2p0_rhc = (TH1D*) ApplyOffAxisWeights(hNuebRatio_rhc,hWeight_nueb_rhc_2p0);

  hNumuRatio_w2p25 = (TH1D*) ApplyOffAxisWeights(hNumuRatio,hWeight_numu_fhc_2p25);
  hNumubRatio_w2p25 = (TH1D*) ApplyOffAxisWeights(hNumubRatio,hWeight_numub_fhc_2p25);
  hNueRatio_w2p25 = (TH1D*) ApplyOffAxisWeights(hNueRatio,hWeight_nue_fhc_2p25);
  hNuebRatio_w2p25 = (TH1D*) ApplyOffAxisWeights(hNuebRatio,hWeight_nueb_fhc_2p25);

  hNumuRatio_w2p25_rhc = (TH1D*) ApplyOffAxisWeights(hNumuRatio_rhc,hWeight_numu_rhc_2p25);
  hNumubRatio_w2p25_rhc = (TH1D*) ApplyOffAxisWeights(hNumubRatio_rhc,hWeight_numub_rhc_2p25);
  hNueRatio_w2p25_rhc = (TH1D*) ApplyOffAxisWeights(hNueRatio_rhc,hWeight_nue_rhc_2p25);
  hNuebRatio_w2p25_rhc = (TH1D*) ApplyOffAxisWeights(hNuebRatio_rhc,hWeight_nueb_rhc_2p25);




  // Draw Spectra

  TLegend *leglike = new TLegend(.51,.51,.75,.75);
  leglike->SetTextSize(0.04);
  leglike->SetFillColor(0);
  hNumuSpec_w2p0->SetLineColor(kRed);
  hNumuSpec_w2p25->SetLineColor(kBlue);
  hNumuSpec->SetLineColor(kBlack);
  leglike->AddEntry(hNumuSpec_w2p0,"#theta_{OA} = 2.0 #circ","L");
  leglike->AddEntry(hNumuSpec_w2p25,"#theta_{OA} = 2.25 #circ","L");
  leglike->AddEntry(hNumuSpec,"#theta_{OA} = 2.25 #circ","L");

  // fhc
  TCanvas *cNumuSpec = new TCanvas("cNumuSpec","cNumuSpec");
  hNumuSpec_w2p0->SetLineColor(kRed);
  hNumuSpec_w2p25->SetLineColor(kBlue);
  hNumuSpec->SetLineColor(kBlack);
  hNumuSpec_w2p0->Draw("histo");
  hNumuSpec_w2p25->Draw("histo same");
  hNumuSpec->Draw("histo same");
  leglike->Draw("same");

  
  TCanvas *cNumubSpec = new TCanvas("cNumubSpec","cNumubSpec");
  hNumubSpec_w2p0->SetLineColor(kRed);
  hNumubSpec_w2p25->SetLineColor(kBlue);
  hNumubSpec->SetLineColor(kBlack);
  hNumubSpec_w2p0->Draw("histo");
  hNumubSpec_w2p25->Draw("histo same");
  hNumubSpec->Draw("histo same");
  leglike->Draw("same");

  TCanvas *cNueSpec = new TCanvas("cNueSpec","cNueSpec");
  hNueSpec_w2p0->SetLineColor(kRed);
  hNueSpec_w2p25->SetLineColor(kBlue);
  hNueSpec->SetLineColor(kBlack);
  hNueSpec_w2p0->Draw("histo");
  hNueSpec_w2p25->Draw("histo same");
  hNueSpec->Draw("histo same");
  leglike->Draw("same");

  TCanvas *cNuebSpec = new TCanvas("cNuebSpec","cNuebSpec");
  hNuebSpec_w2p0->SetLineColor(kRed);
  hNuebSpec_w2p25->SetLineColor(kBlue);
  hNuebSpec->SetLineColor(kBlack);
  hNuebSpec_w2p0->Draw("histo");
  hNuebSpec_w2p25->Draw("histo same");
  hNuebSpec->Draw("histo same");
  leglike->Draw("same");

  // rhc
  TCanvas *cNumuSpec_rhc = new TCanvas("cNumuSpec_rhc","cNumuSpec_rhc");
  hNumuSpec_w2p0_rhc->SetLineColor(kRed);
  hNumuSpec_w2p25_rhc->SetLineColor(kBlue);
  hNumuSpec_rhc->SetLineColor(kBlack);
  hNumuSpec_w2p0_rhc->Draw("histo");
  hNumuSpec_w2p25_rhc->Draw("histo same");
  hNumuSpec_rhc->Draw("histo same");
  leglike->Draw("same");

  TCanvas *cNumubSpec_rhc = new TCanvas("cNumubSpec_rhc","cNumubSpec_rhc");
  hNumubSpec_w2p0_rhc->SetLineColor(kRed);
  hNumubSpec_w2p25_rhc->SetLineColor(kBlue);
  hNumubSpec_rhc->SetLineColor(kBlack);
  hNumubSpec_w2p0_rhc->Draw("histo");
  hNumubSpec_w2p25_rhc->Draw("histo same");
  hNumubSpec_rhc->Draw("histo same");
  leglike->Draw("same");

  TCanvas *cNueSpec_rhc = new TCanvas("cNueSpec_rhc","cNueSpec_rhc");
  hNueSpec_w2p0_rhc->SetLineColor(kRed);
  hNueSpec_w2p25_rhc->SetLineColor(kBlue);
  hNueSpec_rhc->SetLineColor(kBlack);
  hNueSpec_w2p0_rhc->Draw("histo");
  hNueSpec_w2p25_rhc->Draw("histo same");
  hNueSpec_rhc->Draw("histo same");
  leglike->Draw("same");

  TCanvas *cNuebSpec_rhc = new TCanvas("cNuebSpec_rhc","cNuebSpec_rhc");
  hNuebSpec_w2p0_rhc->SetLineColor(kRed);
  hNuebSpec_w2p25_rhc->SetLineColor(kBlue);
  hNuebSpec_rhc->SetLineColor(kBlack);
  hNuebSpec_w2p0_rhc->Draw("histo");
  hNuebSpec_w2p25_rhc->Draw("histo same");
  hNuebSpec_rhc->Draw("histo same");
  leglike->Draw("same");



  // Draw Ratio

  // fhc
  TCanvas *cNumuRatio = new TCanvas("cNumuRatio","cNumuRatio");
  hNumuRatio_w2p0->SetLineColor(kRed);
  hNumuRatio_w2p25->SetLineColor(kBlue);
  hNumuRatio->SetLineColor(kBlack);
  hNumuRatio_w2p0->Draw("histo");
  hNumuRatio_w2p25->Draw("histo same");
  hNumuRatio->Draw("histo same");
  leglike->Draw("same");




  // // Store spectra in the output files

  // TFile *outfile_spectra_offaxis2p0 = new TFile(outfilename_spectra_offaxis2p0,"RECREATE");
  // hNumuSpec_w2p0->Write();
  // hNumubSpec_w2p0->Write();
  // hNueSpec_w2p0->Write();
  // hNuebSpec_w2p0->Write();
  // hNumuRatio_w2p0->Write(); 
  // hNumubRatio_w2p0->Write();
  // hNueRatio_w2p0->Write();
  // hNuebRatio_w2p0->Write();

  // TFile *outfile_spectra_offaxis2p0_rhc = new TFile(outfilename_spectra_offaxis2p0_rhc,"RECREATE");
  // hNumuSpec_w2p0_rhc->Write();
  // hNumubSpec_w2p0_rhc->Write();
  // hNueSpec_w2p0_rhc->Write();
  // hNuebSpec_w2p0_rhc->Write();
  // hNumuRatio_w2p0_rhc->Write(); 
  // hNumubRatio_w2p0_rhc->Write();
  // hNueRatio_w2p0_rhc->Write();
  // hNuebRatio_w2p0_rhc->Write();

  // TFile *outfile_spectra_offaxis2p25 = new TFile(outfilename_spectra_offaxis2p25,"RECREATE");
  // hNumuSpec_w2p25->Write();
  // hNumubSpec_w2p25->Write();
  // hNueSpec_w2p25->Write();
  // hNuebSpec_w2p25->Write();
  // hNumuRatio_w2p25->Write(); 
  // hNumubRatio_w2p25->Write();
  // hNueRatio_w2p25->Write();
  // hNuebRatio_w2p25->Write();

  // TFile *outfile_spectra_offaxis2p25_rhc = new TFile(outfilename_spectra_offaxis2p25_rhc,"RECREATE");
  // hNumuSpec_w2p25_rhc->Write();
  // hNumubSpec_w2p25_rhc->Write();
  // hNueSpec_w2p25_rhc->Write();
  // hNuebSpec_w2p25_rhc->Write();
  // hNumuRatio_w2p25_rhc->Write(); 
  // hNumubRatio_w2p25_rhc->Write();
  // hNueRatio_w2p25_rhc->Write();
  // hNuebRatio_w2p25_rhc->Write();

  return;
}



TH1D *ApplyOffAxisWeights(TH1D *hspectra,TH1D *hweights){
  TH1D *hspectra_w = (TH1D*) hspectra->Clone();

  for(int ispec=0;ispec<hspectra_w->GetNbinsX();ispec++){
   
    double BinCont_spec = hspectra_w->GetBinContent(ispec+1);
    double Cent_spec = hspectra_w->GetBinCenter(ispec+1);
    double LowEdge_spec = hspectra_w->GetBinLowEdge(ispec+1);
    double UpEdge_spec = hspectra_w->GetBinLowEdge(ispec+1)+hspectra_w->GetBinWidth(ispec+1);

    int Bin_w = hweights->FindBin(Cent_spec);
    double w = hweights->GetBinContent(Bin_w);
   
    hspectra_w->SetBinContent(ispec+1,BinCont_spec * w);
    //hspectra_w->SetBinContent(ispec+1,BinCont_spec);

    //cout << BinCont_spec << " " << w << " ";
    //cout << hspectra_w->GetBinContent(ispec+1) << endl;
   }
    
  return hspectra_w;
}
