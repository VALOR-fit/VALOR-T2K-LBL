//
// ------------------------------------------------------------------------------------------------
// Test systematic err cancellation.
//
// The input (systematic) string can be any of:
// 	- nd280_eff
//	- sk_energy_scale
//	- sk_ccqe_eff
//	- sk_ccnonqe_eff
//	- sk_nc_eff
//	- flux_norm
//	- flux_shape
//      - xsec_ccqe
//      - xsec_cc1pi_to_ccqe
//      - xsec_ccoth_to_ccqe
//      - xsec_nc_to_ccqe
//      - nuclear_effect
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

#include "Conventions.h"

void test_syst_err_cancellation(string systematic)
{
  gSystem->Load("libGpad");  
  gSystem->Load("libPhysics");  
  gSystem->Load("libMatrix");   
  gSystem->Load("libHist");     
  gSystem->Load("libGraf");     
  gSystem->Load("libTree");     
  gSystem->Load("libRIO");   
  gSystem->Load("libXMLIO");
  gSystem->Load("libMinuit");   

  gSystem->Load("$T2KANA/lib/libT2KAnaYr1NuMuOscAnalysis.so");

  using namespace std;
  using namespace t2kana;
  using namespace t2kana::yr1numu;

  const double input_sinsq2theta23 = 1.00;
  const double input_dm23sq        = 2.4E-3;

  TRandom3 rand(0);

  // fit params
  PhysParamList * fparm = PhysParamList::NumuDisapFitNoSyst();
  fparm->Update(t2kana::yr1numu::conv::kPhypSinsq2Theta23, input_sinsq2theta23);
  fparm->Update(t2kana::yr1numu::conv::kPhypDm23sq,        input_dm23sq);

  // get superk spectrum predictor & configure
  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read("inputs/analysis_configuration.xml");

  TH1D * effect_plus_1_sigma = new TH1D("effect_plus_1_sigma","",8,0,8);
  TH1D * effect_mins_1_sigma = new TH1D("effect_mins_1_sigma","",8,0,8);

  effect_plus_1_sigma->SetLineColor(kRed );
  effect_mins_1_sigma->SetLineColor(kBlue);

  effect_plus_1_sigma->GetXaxis()->SetBinLabel(1,"nd280 MC: CC");
  effect_plus_1_sigma->GetXaxis()->SetBinLabel(2,"nd280 MC: CC-like NC (bkg)");
  effect_plus_1_sigma->GetXaxis()->SetBinLabel(3,"nd280: CC obs/MC");
  effect_plus_1_sigma->GetXaxis()->SetBinLabel(4,"SK MC: #nu_{#mu} CCQE");
  effect_plus_1_sigma->GetXaxis()->SetBinLabel(5,"SK MC: #nu_{#mu} CC1pi");
  effect_plus_1_sigma->GetXaxis()->SetBinLabel(6,"SK MC: #nu_{#mu} CCOth");
  effect_plus_1_sigma->GetXaxis()->SetBinLabel(7,"SK MC: #nu_{#mu}, #nu_{#tau} NC");
  effect_plus_1_sigma->GetXaxis()->SetBinLabel(8,"Nsk");

  effect_mins_1_sigma->GetXaxis()->SetBinLabel(1,"nd280 MC: CC");
  effect_mins_1_sigma->GetXaxis()->SetBinLabel(2,"nd280 MC: CC-like NC (bkg)");
  effect_mins_1_sigma->GetXaxis()->SetBinLabel(3,"nd280: CC obs/MC");
  effect_mins_1_sigma->GetXaxis()->SetBinLabel(4,"SK MC: #nu_{#mu} CCQE");
  effect_mins_1_sigma->GetXaxis()->SetBinLabel(5,"SK MC: #nu_{#mu} CC1pi");
  effect_mins_1_sigma->GetXaxis()->SetBinLabel(6,"SK MC: #nu_{#mu} CCOth");
  effect_mins_1_sigma->GetXaxis()->SetBinLabel(7,"SK MC: #nu_{#mu}, #nu_{#tau} NC");
  effect_mins_1_sigma->GetXaxis()->SetBinLabel(8,"Nsk");

  //
  // nominal
  //

  double nd280cc_0       = skexp->GetND280Inputs()->NExpCC();
  double nd280ncbkg_0    = skexp->GetND280Inputs()->NExpNCBkg();
  double nd280corr_0     = skexp->GetND280Inputs()->CorrectionFactor();
  double nsknumuccqe_0   = skexp->AsTH1D(*fparm, conv::kSKsFCFV1RingMuLike, conv::kSKmNuMuCCQE,  false)->Integral();
  double nsknumucc1pi_0  = skexp->AsTH1D(*fparm, conv::kSKsFCFV1RingMuLike, conv::kSKmNuMuCC1pi, false)->Integral();
  double nsknumuccoth_0  = skexp->AsTH1D(*fparm, conv::kSKsFCFV1RingMuLike, conv::kSKmNuMuCCoth, false)->Integral();
  double nsknumunc_0     = skexp->AsTH1D(*fparm, conv::kSKsFCFV1RingMuLike, conv::kSKmNC,        false)->Integral();
  double nsk_0           = skexp->AsTH1D(*fparm, conv::kSKsFCFV1RingMuLike)->Integral();

  skexp->TweakSystematic(systematic,+1.);

  double nd280cc_p       = skexp->GetND280Inputs()->NExpCC();
  double nd280ncbkg_p    = skexp->GetND280Inputs()->NExpNCBkg();
  double nd280corr_p     = skexp->GetND280Inputs()->CorrectionFactor();
  double nsknumuccqe_p   = skexp->AsTH1D(*fparm, conv::kSKsFCFV1RingMuLike, conv::kSKmNuMuCCQE,  false)->Integral();
  double nsknumucc1pi_p  = skexp->AsTH1D(*fparm, conv::kSKsFCFV1RingMuLike, conv::kSKmNuMuCC1pi, false)->Integral();
  double nsknumuccoth_p  = skexp->AsTH1D(*fparm, conv::kSKsFCFV1RingMuLike, conv::kSKmNuMuCCoth, false)->Integral();
  double nsknumunc_p     = skexp->AsTH1D(*fparm, conv::kSKsFCFV1RingMuLike, conv::kSKmNC,        false)->Integral();
  double nsk_p           = skexp->AsTH1D(*fparm, conv::kSKsFCFV1RingMuLike)->Integral();

  skexp->TweakSystematic(systematic,-1.);

  double nd280cc_m       = skexp->GetND280Inputs()->NExpCC();
  double nd280ncbkg_m    = skexp->GetND280Inputs()->NExpNCBkg();
  double nd280corr_m     = skexp->GetND280Inputs()->CorrectionFactor();
  double nsknumuccqe_m   = skexp->AsTH1D(*fparm, conv::kSKsFCFV1RingMuLike, conv::kSKmNuMuCCQE,  false)->Integral();
  double nsknumucc1pi_m  = skexp->AsTH1D(*fparm, conv::kSKsFCFV1RingMuLike, conv::kSKmNuMuCC1pi, false)->Integral();
  double nsknumuccoth_m  = skexp->AsTH1D(*fparm, conv::kSKsFCFV1RingMuLike, conv::kSKmNuMuCCoth, false)->Integral();
  double nsknumunc_m     = skexp->AsTH1D(*fparm, conv::kSKsFCFV1RingMuLike, conv::kSKmNC,        false)->Integral();
  double nsk_m           = skexp->AsTH1D(*fparm, conv::kSKsFCFV1RingMuLike)->Integral();

  double frac_nd280cc_p       = (nd280cc_0      > 0) ? (nd280cc_p      - nd280cc_0     ) / nd280cc_0      : 0.;
  double frac_nd280ncbkg_p    = (nd280ncbkg_0   > 0) ? (nd280ncbkg_p   - nd280ncbkg_0  ) / nd280ncbkg_0   : 0.;
  double frac_nd280corr_p     = (nd280corr_0    > 0) ? (nd280corr_p    - nd280corr_0   ) / nd280corr_0    : 0.;
  double frac_nsknumuccqe_p   = (nsknumuccqe_0  > 0) ? (nsknumuccqe_p  - nsknumuccqe_0 ) / nsknumuccqe_0  : 0.;
  double frac_nsknumucc1pi_p  = (nsknumucc1pi_0 > 0) ? (nsknumucc1pi_p - nsknumucc1pi_0) / nsknumucc1pi_0 : 0.;
  double frac_nsknumuccoth_p  = (nsknumuccoth_0 > 0) ? (nsknumuccoth_p - nsknumuccoth_0) / nsknumuccoth_0 : 0.;
  double frac_nsknumunc_p     = (nsknumunc_0    > 0) ? (nsknumunc_p    - nsknumunc_0   ) / nsknumunc_0    : 0.;
  double frac_nsk_p           = (nsk_0          > 0) ? (nsk_p          - nsk_0         ) / nsk_0          : 0.;

  double frac_nd280cc_m       = (nd280cc_0      > 0) ? (nd280cc_m      - nd280cc_0     ) / nd280cc_0      : 0.;
  double frac_nd280ncbkg_m    = (nd280ncbkg_0   > 0) ? (nd280ncbkg_m   - nd280ncbkg_0  ) / nd280ncbkg_0   : 0.;
  double frac_nd280corr_m     = (nd280corr_0    > 0) ? (nd280corr_m    - nd280corr_0   ) / nd280corr_0    : 0.;
  double frac_nsknumuccqe_m   = (nsknumuccqe_0  > 0) ? (nsknumuccqe_m  - nsknumuccqe_0 ) / nsknumuccqe_0  : 0.;
  double frac_nsknumucc1pi_m  = (nsknumucc1pi_0 > 0) ? (nsknumucc1pi_m - nsknumucc1pi_0) / nsknumucc1pi_0 : 0.;
  double frac_nsknumuccoth_m  = (nsknumuccoth_0 > 0) ? (nsknumuccoth_m - nsknumuccoth_0) / nsknumuccoth_0 : 0.;
  double frac_nsknumunc_m     = (nsknumunc_0    > 0) ? (nsknumunc_m    - nsknumunc_0   ) / nsknumunc_0    : 0.;
  double frac_nsk_m           = (nsk_0          > 0) ? (nsk_m          - nsk_0         ) / nsk_0          : 0.;

  effect_plus_1_sigma -> SetBinContent(1, 100*frac_nd280cc_p       );
  effect_plus_1_sigma -> SetBinContent(2, 100*frac_nd280ncbkg_p    ); 
  effect_plus_1_sigma -> SetBinContent(3, 100*frac_nd280corr_p     ); 
  effect_plus_1_sigma -> SetBinContent(4, 100*frac_nsknumuccqe_p   ); 
  effect_plus_1_sigma -> SetBinContent(5, 100*frac_nsknumucc1pi_p  ); 
  effect_plus_1_sigma -> SetBinContent(6, 100*frac_nsknumuccoth_p  ); 
  effect_plus_1_sigma -> SetBinContent(7, 100*frac_nsknumunc_p     ); 
  effect_plus_1_sigma -> SetBinContent(8, 100*frac_nsk_p           ); 

  effect_mins_1_sigma -> SetBinContent(1, 100*frac_nd280cc_m       );
  effect_mins_1_sigma -> SetBinContent(2, 100*frac_nd280ncbkg_m    ); 
  effect_mins_1_sigma -> SetBinContent(3, 100*frac_nd280corr_m     ); 
  effect_mins_1_sigma -> SetBinContent(4, 100*frac_nsknumuccqe_m   ); 
  effect_mins_1_sigma -> SetBinContent(5, 100*frac_nsknumucc1pi_m  ); 
  effect_mins_1_sigma -> SetBinContent(6, 100*frac_nsknumuccoth_m  ); 
  effect_mins_1_sigma -> SetBinContent(7, 100*frac_nsknumunc_m     ); 
  effect_mins_1_sigma -> SetBinContent(8, 100*frac_nsk_m           ); 

  double max = 1.4 * TMath::Max( effect_plus_1_sigma->GetMaximum(), effect_mins_1_sigma->GetMaximum() );
  double min = 1.4 * TMath::Min( effect_plus_1_sigma->GetMinimum(), effect_mins_1_sigma->GetMinimum() );

  effect_plus_1_sigma->GetYaxis()->SetRangeUser(min,max);
  effect_mins_1_sigma->GetYaxis()->SetRangeUser(min,max);
  effect_plus_1_sigma->GetYaxis()->SetTitle("%");
  effect_mins_1_sigma->GetYaxis()->SetTitle("%");

  TCanvas * c = new TCanvas("c","",20,20,600,600);
  effect_plus_1_sigma->Draw();
  effect_mins_1_sigma->Draw("SAME");

  gStyle->SetOptStat(0);

  TLegend * legend = new TLegend(0.2,.6,0.4,.9);
  legend->SetFillColor(0);
  legend->SetBorderSize(0);
  legend->SetHeader(systematic.c_str());
  legend->AddEntry(effect_plus_1_sigma, "+1 #sigma", "L");
  legend->AddEntry(effect_mins_1_sigma, "-1 #sigma", "L");
  legend->Draw();

  c->Update();
}
