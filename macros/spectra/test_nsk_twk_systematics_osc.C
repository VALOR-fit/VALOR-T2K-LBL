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

//#include "Conventions.h"

void test_nsk_twk_systematics_osc
(
  const char * config_file = "analysis_configuration.xml"
)
{
  using namespace std;
  using namespace t2kana;

  conv::SuperkSample_t s = conv::kSKsFCFV1RingMuLike;

  const char * config = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);
  cout << "Reading analysis configuration from: " << config << endl;

  const int nexperiments = 1000;

  double mixing23_min  = 0.0;
  double mixing23_max  = 1.0;
  double mixing23_step = 0.01;
  double dm23sq_min    = 1E-3;
  double dm23sq_max    = 6E-3;
  double dm23sq_step   = 1E-4;

  // make sure that min and max are at the bin centre
  int mixing23_nbins = 1 + (int) ((mixing23_max - mixing23_min) / mixing23_step);
  int dm23sq_nbins        = 1 + (int) ((dm23sq_max        - dm23sq_min       ) / dm23sq_step       );

  TH2D * frac_nsk_err = 
    new TH2D("frac_nsk_err", "",
       mixing23_nbins, 
       mixing23_min - 0.5*mixing23_step, 
       mixing23_max + 0.5*mixing23_step,
       dm23sq_nbins, 
       dm23sq_min - 0.5*dm23sq_step, 
       dm23sq_max + 0.5*dm23sq_step);

  const int nsystematics = 21;
  string systematic[nsystematics] = {
    "sk_1Rmu_ccqe_ring_counting_eff", //0
    "sk_1Rmu_ccqe_other_eff",         //1
    "sk_1Rmu_ccnonqe_eff",            //2
    "sk_1Rmu_nc_eff",                 //3
    "sk_1Rmu_nuecc_eff",              //4
    "sk_energy_scale",           //5
    "banff",                     //6
    "xsec_pf",                   //7
    "xsec_sf",                   //8
    "xsec_wshape",               //9
    "xsec_shape_ccoth",          //10
    "xsec_pilessdelta",          //11
    "xsec_enushape",             //12
    "xsec_norm_ccqeE1",          //13
    "xsec_norm_ccqeE2",          //14
    "xsec_norm_cc1piE1",         //15
    "xsec_norm_cccoh",           //16
    "xsec_norm_nccoh",           //17
    "xsec_norm_ncoth",           //18
    "xsec_norm_nuetonumu",       //19
    "xsec_norm_antinutonu",      //19
    "fsi"                        //20
  };
  bool generate_random[nsystematics] = {
    true,  /* sk_1Rmu_ccqe_ring_counting_eff */
    false, /* sk_1Rmu_ccqe_other_eff         */
    false, /* sk_1Rmu_ccnonqe_eff            */
    false, /* sk_1Rmu_nc_eff                 */
    false, /* sk_1Rmu_nuecc_eff              */
    false, /* sk_energy_scale           */
    true,  /* banff                     */
    false, /* xsec_pf                   */
    false, /* xsec_sf                   */
    false, /* xsec_wshape               */
    false, /* xsec_shape_ccoth          */
    false, /* xsec_pilessdelta          */
    false, /* xsec_enushape             */
    false, /* xsec_norm_ccqeE1          */
    false, /* xsec_norm_ccqeE2          */
    false, /* xsec_norm_cc1piE1         */
    false, /* xsec_norm_cccoh           */
    false, /* xsec_norm_nccoh           */
    false, /* xsec_norm_ncoth           */
    false, /* xsec_norm_nuetonumu       */
    false, /* xsec_norm_antinutonu      */
    true   /* fsi                       */
  };
  bool skip[nsystematics] = {
    false, /* sk_1Rmu_ccqe_ring_counting_eff */
    false, /* sk_1Rmu_ccqe_other_eff         */
    false, /* sk_1Rmu_ccnonqe_eff            */
    false, /* sk_1Rmu_nc_eff                 */
    false, /* sk_1Rmu_nuecc_eff              */
    false, /* sk_energy_scale           */
    false, /* banff                     */
    false, /* xsec_pf                   */
    false, /* xsec_sf                   */
    false, /* xsec_wshape               */
    false, /* xsec_shape_ccoth          */
    false, /* xsec_pilessdelta          */
    false, /* xsec_enushape             */
    false, /* xsec_norm_ccqeE1          */
    false, /* xsec_norm_ccqeE2          */
    false, /* xsec_norm_cc1piE1         */
    false, /* xsec_norm_cccoh           */
    false, /* xsec_norm_nccoh           */
    false, /* xsec_norm_ncoth           */
    false, /* xsec_norm_nuetonumu       */
    false, /* xsec_norm_antinutonu      */
    false  /* fsi                       */
  };


  // superk spectrum predictor
  SuperkExpectation skexp;
  skexp.Read(config);

  // phys hypothesis
  PhysHypothesis ph(t2kana::conv::kHypNuMuDisapp2F);

  // generate predictions
  for(int i = 1; i <= frac_nsk_err->GetNbinsX(); i++) {
    for(int j = 1; j <= frac_nsk_err->GetNbinsY(); j++) {

       double mixing23 = frac_nsk_err->GetXaxis()->GetBinCenter(i);
       double dm23sq   = frac_nsk_err->GetYaxis()->GetBinCenter(j);
       mixing23 = TMath::Max(0.,mixing23);
       ph.SetCurrValue(t2kana::conv::kPhypMixing23, mixing23);
       ph.SetCurrValue(t2kana::conv::kPhypDm23sq,   dm23sq);

       cout << "* Dm23^2 = " << dm23sq 
            << " eV^2, sin^2(2theta) = " << mixing23 << endl;

       // generate predictions
       skexp.ResetSystematicVariations();
       skexp.Generate(ph);

       // integrate Ereco distribution & get nominal nsk
       double nsk = skexp.CurrPredictedErecoSpectrum(s)->Integral();

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
                    double nsk = skexp.CurrPredictedErecoSpectrum(s)->Integral();
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
                double nsk1 = skexp.CurrPredictedErecoSpectrum(s)->Integral();
                skexp.ResetSystematicVariations();
                skexp.TweakSystematic(systematic[isystematic], -1);
                skexp.Generate(ph);
                double nsk2 = skexp.CurrPredictedErecoSpectrum(s)->Integral();
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
