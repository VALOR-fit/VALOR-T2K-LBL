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

void test_nsk_twk_systematics_osc
(
  const char * config_file = "analysis_configuration.xml"
)
{
  using namespace std;
  using namespace t2kana;

  const char * config = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);
  cout << "Reading analysis configuration from: " << config << endl;

  const int nexperiments = 1000;

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

  const int nsystematics = 13;
  string systematic[nsystematics] = {
        "sk_ccqe_other_eff",
        "sk_ccqe_ring_counting_eff",
        "sk_ccnonqe_eff",
        "sk_nc_eff",
        "sk_nuecc_eff",
        "nd280_eff",
        "nd280_norm",
        "flux_norm",
        "xsec_ccqe",
        "xsec_cc1pi_to_ccqe",
        "xsec_ccoth_to_ccqe",
        "xsec_nc_to_ccqe",
        "nuclear_effect"
  };
  bool generate_random[nsystematics] = {
        false,  /* sk_ccqe_other_eff          */
        true,   /* sk_ccqe_ring_counting_eff  */
        false,  /* sk_ccnonqe_eff             */
        false,  /* sk_nc_eff                  */
        false,  /* sk_nuecc_eff               */
        false,  /* nd280_eff                  */
        false,  /* nd280_norm                 */
        false,  /* flux_norm                  */
        false,  /* xsec_ccqe                  */
        false,  /* xsec_cc1pi_to_ccqe         */
        false,  /* xsec_ccoth_to_ccqe         */
        false,  /* xsec_nc_to_ccqe            */
        false   /* nuclear_effect             */
  };
  bool skip[nsystematics] = {
        true,   /* sk_ccqe_other_eff         */
        true,   /* sk_ccqe_ring_counting_eff */
        true,   /* sk_ccnonqe_eff            */
        true,   /* sk_nc_eff                 */
        true,   /* sk_nuecc_eff              */
        true,   /* nd280_eff                 */
        true,   /* nd280_norm                */
        true,   /* flux_norm                 */
        false,  /* xsec_ccqe                 */
        false,  /* xsec_cc1pi_to_ccqe        */
        false,  /* xsec_ccoth_to_ccqe        */
        false,  /* xsec_nc_to_ccqe           */
        false   /* nuclear_effect            */
  };


  // superk spectrum predictor
  SuperkExpectation skexp;
  skexp.Read(config);

  // phys hypothesis
  PhysHypothesis ph(t2kana::conv::kHypNuMuDisapp2F);

  // generate predictions
  for(int i = 1; i <= frac_nsk_err->GetNbinsX(); i++) {
    for(int j = 1; j <= frac_nsk_err->GetNbinsY(); j++) {

       double sinsq2theta23 = frac_nsk_err->GetXaxis()->GetBinCenter(i);
       double dm23sq        = frac_nsk_err->GetYaxis()->GetBinCenter(j);
       sinsq2theta23 = TMath::Max(0.,sinsq2theta23);
       ph.SetCurrValue(t2kana::conv::kPhypSinsq2Theta23, sinsq2theta23);
       ph.SetCurrValue(t2kana::conv::kPhypDm23sq,        dm23sq);

       cout << "* Dm23^2 = " << dm23sq 
            << " eV^2, sin^2(2theta) = " << sinsq2theta23 << endl;

       // generate predictions
       skexp.ResetSystematicVariations();
       skexp.Generate(ph);

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

            if(generate_random[isystematic]) {
                // generate toy MC expts and calculate Nsk  
                double nskexp[nexperiments];
                for(int iexperiment = 0; iexperiment < nexperiments; iexperiment++) {
                    if(iexperiment%100==0) {
                        cout << "Current experiment: " << iexperiment << endl;
                    }
                    skexp.ResetSystematicVariations(); 
                    skexp.TweakSystematic(systematic[isystematic]);
                    skexp.Generate(ph);
                    double nsk = skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike)->Integral();
                    nskexp[iexperiment] = nsk;
                }//iexperiment
                // calculate standard deviation
                double nsk_sigma = 0;
                for(int iexperiment = 0; iexperiment < nexperiments; iexperiment++) {
                    nsk_sigma += TMath::Power(nskexp[iexperiment] - nsk, 2.);
                }
                nsk_sigma = TMath::Sqrt(nsk_sigma/nexperiments);
                // set Nsk range
                nsk_min[isystematic] = nsk - nsk_sigma;
                nsk_max[isystematic] = nsk + nsk_sigma;
            }
            else {
                skexp.ResetSystematicVariations();
                skexp.TweakSystematic(systematic[isystematic], +1);
                skexp.Generate(ph);
                double nsk1 = skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike)->Integral();
                skexp.ResetSystematicVariations();
                skexp.TweakSystematic(systematic[isystematic], -1);
                skexp.Generate(ph);
                double nsk2 = skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike)->Integral();
                nsk_min[isystematic] = TMath::Min(nsk1,nsk2);
                nsk_max[isystematic] = TMath::Max(nsk1,nsk2);
            }
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
