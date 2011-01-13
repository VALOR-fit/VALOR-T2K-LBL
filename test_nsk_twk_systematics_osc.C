//
// ------------------------------------------------------------------------------------------------
// Plot fractional error on number of events predicted at SuperK (NSK) vs oscillation params.
// The `flux shape' and 'sk energy scale' systematics (conserve NSK) are not included at all.
// In default mode, the error includes only cross section systematics (so the estimated error 
// can be compared against the NEUT/GENIE differences)
// 
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

#include "Conventions.h"

void test_nsk_twk_systematics_osc()
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

  gSystem->Load("libT2KDisapp2010a");

  using namespace std;
  using namespace t2kana;
  using namespace t2kana::numu10a;

  double sinsq2theta23_min  = 0.0;
  double sinsq2theta23_max  = 1.0;
  double sinsq2theta23_step = 0.01;
  double dm23sq_min         = 1E-3;
  double dm23sq_max         = 6E-3;
  double dm23sq_step        = 1E-4;

 // make sure that min and max are at the bin centre
 int sinsq2theta23_nbins = 1 + (int) ((sinsq2theta23_max - sinsq2theta23_min) / sinsq2theta23_step);
 int dm23sq_nbins        = 1 + (int) ((dm23sq_max        - dm23sq_min       ) / dm23sq_step       );

 TH2D * frac_nsk_err = 
   new TH2D("frac_nsk_err", "",
       sinsq2theta23_nbins, 
       sinsq2theta23_min - 0.5*sinsq2theta23_step, 
       sinsq2theta23_max + 0.5*sinsq2theta23_step,
       dm23sq_nbins, 
       dm23sq_min - 0.5*dm23sq_step, 
       dm23sq_max + 0.5*dm23sq_step);

  const int nsystematics = 11;
  string systematic[nsystematics] = {
        "sk_ccqe_eff",
        "sk_ccnonqe_eff",
        "sk_nc_eff",
        "sk_nuecc_eff",
        "nd280_eff",
        "flux_norm",
        "xsec_ccqe",
        "xsec_cc1pi_to_ccqe",
        "xsec_ccoth_to_ccqe",
        "xsec_nc_to_ccqe",
        "nuclear_effect"
  };
  bool skip[nsystematics] = {
        true,   /* sk_ccqe_eff        */
        true,   /* sk_ccnonqe_eff     */
        true,   /* sk_nc_eff          */
        true,   /* sk_nuecc_eff       */
        true,   /* nd280_eff          */
        true,   /* flux_norm          */
        false,  /* xsec_ccqe          */
        false,  /* xsec_cc1pi_to_ccqe */
        false,  /* xsec_ccoth_to_ccqe */
        false,  /* xsec_nc_to_ccqe    */
        false   /* nuclear_effect     */
  };


  // superk spectrum predictor
  SuperkExpectation skexp;
  skexp.Read("inputs/analysis_configuration.xml");

  // fit params
  PhysParamList * fparm = PhysParamList::NumuDisapFitNoSyst();

  // generate predictions
  for(int i = 1; i <= frac_nsk_err->GetNbinsX(); i++) {
    for(int j = 1; j <= frac_nsk_err->GetNbinsY(); j++) {
       double sinsq2theta23 = frac_nsk_err->GetXaxis()->GetBinCenter(i);
       double dm23sq        = frac_nsk_err->GetYaxis()->GetBinCenter(j);
       sinsq2theta23 = TMath::Max(0.,sinsq2theta23);
       fparm->Update(t2kana::numu10a::conv::kPhypSinsq2Theta23, sinsq2theta23);
       fparm->Update(t2kana::numu10a::conv::kPhypDm23sq,        dm23sq);

       cout << "* Dm23^2 = " << dm23sq 
            << " eV^2, sin^2(2theta) = " << sinsq2theta23 << endl;

       // generate predictions
       skexp.ResetSystematicVariations();
       skexp.Generate(*fparm);

       // integrate Ereco distribution & get nominal nsk
       double nsk = skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike)->Integral();

       cout << "\t nsk = " << nsk << endl;

       // get predictions with tweaked systematic params
       double nsk_min[nsystematics];
       double nsk_max[nsystematics];
       for(int isystematic=0; isystematic<nsystematics; isystematic++) {

            nsk_min[isystematic] = nsk;
            nsk_max[isystematic] = nsk;

            if(skip[isystematic]) {
              continue;
            }

            skexp.ResetSystematicVariations();
            skexp.TweakSystematic(systematic[isystematic], +1);
            skexp.Generate(*fparm);
            double nsk1 = skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike)->Integral();
            skexp.ResetSystematicVariations();
            skexp.TweakSystematic(systematic[isystematic], -1);
            skexp.Generate(*fparm);
            double nsk2 = skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike)->Integral();
            nsk_min[isystematic] = TMath::Min(nsk1,nsk2);
            nsk_max[isystematic] = TMath::Max(nsk1,nsk2);
	    cout << "\t " << systematic[isystematic] 
                 << " -> (" << nsk_min[isystematic] << ", " << nsk_max[isystematic] << ")" << endl;
       }//isystematic
       double dnsk_syst_pos = 0;
       double dnsk_syst_neg = 0;
       for(int isystematic=0; isystematic<nsystematics; isystematic++) {
            double pos = nsk_max[isystematic]-nsk;
            double neg = nsk_min[isystematic]-nsk;
            dnsk_syst_pos += (pos*pos);
            dnsk_syst_neg += (neg*neg);
       }
       dnsk_syst_pos = TMath::Sqrt(dnsk_syst_pos);
       dnsk_syst_neg = TMath::Sqrt(dnsk_syst_neg);
       double dnsk_syst = (dnsk_syst_pos+dnsk_syst_pos)/2.;

       double frac_err = 100*dnsk_syst/nsk;

       cout << "\t frac error = " << frac_err << endl;

       frac_nsk_err->SetBinContent(i,j,frac_err);
    }//j
  }//i

  gStyle->SetOptStat(0);

  frac_nsk_err->GetXaxis()->SetTitle("sin^{2}2#theta_{23}");
  frac_nsk_err->GetYaxis()->SetTitle("|#Deltam^{2}_{23}| (eV^{2}/c^{4})");

  frac_nsk_err->Draw("colz");

  TFile f("./frac_nsk_err_vs_osc_params.root","recreate");
  frac_nsk_err->Write("frac_nsk_err");
  f.Close();
}
