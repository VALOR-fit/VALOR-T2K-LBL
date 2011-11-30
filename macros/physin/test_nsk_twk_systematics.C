//
// ------------------------------------------------------------------------------------------------
// Summarize Nsk change as a result of tweaking systematic parameters.
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

//#include "Conventions.h"

void test_nsk_twk_systematics
(
  const string sample_str          = "1Rmu",
  const char * hyp                 = "3f_osc",
  const double input_sinsq2theta23 = 1.00,
  const double input_sinsq2theta13 = 0.1,
  const double input_sinsq2theta12 = 0.8704,
  const double input_dm23sq        = 2.4E-3,
  const double input_dm12sq        = 7.6E-5,
  const double input_deltacp       = 0,
  const char * config_file         = "analysis_configuration.xml"
)
{
  gROOT->Reset();

  using namespace std;
  using namespace t2kana;

  conv::SuperkSample_t s = util::SKSampleFromString(sample_str);

  const int nexperiments = 10;//5000;

  const int nsystematics = 22;
  string systematic[nsystematics] = {
    "sk_ccqe_ring_counting_eff", //0
    "sk_ccqe_other_eff",         //1
    "sk_ccnonqe_eff",            //2
    "sk_nc_eff",                 //3
    "sk_nuecc_eff",              //4
    "sk_1Resig_eff",             //5
    "sk_1Rebkg_eff",             //6
    "sk_energy_scale",           //7
    "nd280_eff",                 //8
    "nd280_norm",                //9
    "flux_norm",                 //10
    "flux_shape",                //11 
    "xsec_ccqe",                 //12
    "xsec_cc1pi_to_ccqe",        //13
    "xsec_cccoh_to_ccqe",        //14
    "xsec_ccoth_to_ccqe",        //15
    "xsec_nc_to_ccqe",           //16
    "xsec_nc1pi0_to_ccqe",       //17
    "xsec_nccoh_to_ccqe",        //18
    "xsec_ncoth_to_ccqe",        //19
    "xsec_nue_to_numu",          //20
    "nuclear_effect"             //21
  };
  bool generate_random[nsystematics] = {
	true,   /* sk_ccqe_ring_counting_eff  */
	false,  /* sk_ccqe_other_eff          */
	false,  /* sk_ccnonqe_eff             */
	false,  /* sk_nc_eff                  */
	false,  /* sk_nuecc_eff               */
	false,  /* sk_1Resig_eff              */
	false,  /* sk_1Rebkg_eff              */
	false,  /* sk_energy_scale            */
	false,  /* nd280_eff                  */
	false,  /* nd280_norm                 */
	false,  /* flux_norm                  */
	true,   /* flux_shape                 */
	false,  /* xsec_ccqe                  */
	false,  /* xsec_cc1pi_to_ccqe         */
	false,  /* xsec_cccoh_to_ccqe         */
	false,  /* xsec_ccoth_to_ccqe         */
	false,  /* xsec_nc_to_ccqe            */
	false,  /* xsec_nc1pi0_to_ccqe        */
	false,  /* xsec_nccoh_to_ccqe         */
	false,  /* xsec_ncoth_to_ccqe         */
	false,  /* xsec_nue_to_numu           */
	false   /* nuclear_effect             */
  };
  bool skip[nsystematics] = {
	false,  /* sk_ccqe_ring_counting_eff          */
	false,  /* sk_ccqe_other_eff                  */
	false,  /* sk_ccnonqe_eff                     */
	false,  /* sk_nc_eff                          */
	false,  /* sk_nuecc_eff                       */
	false,  /* sk_1Resig_eff                      */
	false,  /* sk_1Rebkg_eff                      */
	true,   /* sk_energy_scale (no effect on Nsk) */
	false,  /* nd280_eff                          */
	false,  /* nd280_norm                         */
	false,  /* flux_norm                          */
	true,   /* flux_shape (no effect on Nsk)      */
	false,  /* xsec_ccqe                          */
	false,  /* xsec_cc1pi_to_ccqe                 */
	false,  /* xsec_cccoh_to_ccqe                 */
	false,  /* xsec_ccoth_to_ccqe                 */
	false,  /* xsec_nc_to_ccqe                    */
	false,  /* xsec_nc1pi0_to_ccqe                */
	false,  /* xsec_nccoh_to_ccqe                 */
	false,  /* xsec_ncoth_to_ccqe                 */
	false,  /* xsec_nue_to_numu                   */
	false   /* nuclear_effect                     */
  };

  if(s == conv::kSKsFCFV1RingMuLike) {
    skip[5]  = true;
    skip[6]  = true;
    skip[14] = true;
    skip[17] = true;
    skip[18] = true;
    skip[19] = true;
    skip[20] = true;
  }
  else if (s == conv::kSKsFCFV1RingELike) {
    skip[0]  = true;
    skip[1]  = true;
    skip[2]  = true;
    skip[3]  = true;
    skip[4]  = true;
    skip[16] = true;
  }
  else if (s == conv::kSKsUndefined) {
    cout << "Sample undefined" << endl;
    return;
  }
  else {
    cout << "Sample not supported" << endl;
    return;
  }

  double nsk_min[nsystematics];
  double nsk_max[nsystematics];

  // phys hypothesis
  PhysHypothesis ph(hyp);
  ph.SetCurrValue(t2kana::conv::kPhypSinsq2Theta23, input_sinsq2theta23);
  ph.SetCurrValue(t2kana::conv::kPhypDm23sq,        input_dm23sq);
  if(hyp == conv::kHypOsc3F) {
    ph.SetCurrValue(t2kana::conv::kPhypSinsq2Theta13, input_sinsq2theta13);
    ph.SetCurrValue(t2kana::conv::kPhypSinsq2Theta12, input_sinsq2theta12);
    ph.SetCurrValue(t2kana::conv::kPhypDm12sq,        input_dm12sq);
    ph.SetCurrValue(t2kana::conv::kPhypDeltaCP,       input_deltacp);
  }

  // get & configure spectrum predictor
  const char * config = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);
  cout << "Reading analysis configuration from: " << config << endl;
  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read(config);
 
  // get default prediction
  skexp->Generate(ph);
  double nsk_deflt = skexp->CurrPredictedErecoSpectrum(s)->Integral();
  
  // get predictions with tweaked systematic params
  for(int isystematic=0; isystematic<nsystematics; isystematic++) {

	cout << "Current systematic: " << systematic[isystematic] << endl;

        nsk_min[isystematic] = nsk_deflt;
        nsk_max[isystematic] = nsk_deflt;

	if(skip[isystematic]) {
            cout << "... skipping " << endl;
	    continue;
        }

	if(generate_random[isystematic]) {
            nsk_min[isystematic] = 99999999.;
            nsk_max[isystematic] = 0.;
            // generate toy MC expts and calculate Nsk
            double nskexp[nexperiments];
            for(int iexperiment = 0; iexperiment < nexperiments; iexperiment++) {
                    if(iexperiment%100==0) {
			cout << "Current experiment: " << iexperiment << endl;
                    }
		    skexp->ResetSystematicVariations();
	     	    skexp->TweakSystematic(systematic[isystematic]);
	  	    skexp->Generate(ph);
  		    double nsk = skexp->CurrPredictedErecoSpectrum(s)->Integral();
                    nskexp[iexperiment] = nsk;
            }//iexperiment
            // calculate standard deviation
            double nsk_sigma = 0;
            for(int iexperiment = 0; iexperiment < nexperiments; iexperiment++) {
                nsk_sigma += TMath::Power(nskexp[iexperiment] - nsk_deflt, 2.);
            }
            nsk_sigma = TMath::Sqrt(nsk_sigma/nexperiments);
            // set Nsk range
	    nsk_min[isystematic] = nsk_deflt - nsk_sigma;
            nsk_max[isystematic] = nsk_deflt + nsk_sigma;

//	    nsk_min[isystematic] = TMath::Min(nsk, nsk_min[isystematic]);
//          nsk_max[isystematic] = TMath::Max(nsk, nsk_max[isystematic]);

	} else {
	    skexp->ResetSystematicVariations();
     	    skexp->TweakSystematic(systematic[isystematic], +1);
  	    skexp->Generate(ph);
  	    double nsk1 = skexp->CurrPredictedErecoSpectrum(s)->Integral();
	    skexp->ResetSystematicVariations();
     	    skexp->TweakSystematic(systematic[isystematic], -1);
  	    skexp->Generate(ph);
  	    double nsk2 = skexp->CurrPredictedErecoSpectrum(s)->Integral();
            nsk_min[isystematic] = TMath::Min(nsk1,nsk2);
            nsk_max[isystematic] = TMath::Max(nsk1,nsk2);
	}	
  }//isystematic

  cout << endl;
  cout << "Nsk(default) = " << nsk_deflt << endl;
  for(int isystematic=0; isystematic<nsystematics; isystematic++) {
        if(skip[isystematic] == true) continue;
	cout << " * " << systematic[isystematic] << " : Nsk = (" 
	     << nsk_min[isystematic] << ", " << nsk_max[isystematic] << ")" << endl;
  }
  cout << endl;
  cout << "Fractional change: " << endl;
  for(int isystematic=0; isystematic<nsystematics; isystematic++) {
        if(skip[isystematic] == true) continue;
	cout << " * " << systematic[isystematic] << " : dNsk/Nsk = (" 
             << 100*(nsk_max[isystematic]-nsk_deflt)/nsk_deflt << " %," 
             << 100*(nsk_min[isystematic]-nsk_deflt)/nsk_deflt << " %)" 
             << endl;
  }
  cout << endl;
  cout << "Total factional change: " << endl;
  double dnsk_syst_pos = 0;
  double dnsk_syst_neg = 0;
  for(int isystematic=0; isystematic<nsystematics; isystematic++) {
     double pos = nsk_max[isystematic]-nsk_deflt;
     double neg = nsk_min[isystematic]-nsk_deflt;
     dnsk_syst_pos += (pos*pos);
     dnsk_syst_neg += (neg*neg);
  }
  dnsk_syst_pos = TMath::Sqrt(dnsk_syst_pos);
  dnsk_syst_neg = TMath::Sqrt(dnsk_syst_neg);
  cout << "(+" << 100*dnsk_syst_pos/nsk_deflt << "%, -" << 100*dnsk_syst_neg/nsk_deflt << ")" << endl;

  double dnsk_stat = TMath::Sqrt(nsk_deflt);

  cout << "Nsk = " << nsk_deflt 
       << " +/- " << dnsk_stat << " (stat)"
       << " + " << dnsk_syst_pos 
       << " - " << dnsk_syst_neg << " (syst)" << endl;

  cout << "+dNsk(total) = " << TMath::Sqrt(dnsk_stat*dnsk_stat + dnsk_syst_pos*dnsk_syst_pos) << endl;
  cout << "-dNsk(total) = " << TMath::Sqrt(dnsk_stat*dnsk_stat + dnsk_syst_neg*dnsk_syst_neg) << endl;
}
