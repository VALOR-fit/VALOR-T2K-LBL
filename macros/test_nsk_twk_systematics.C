//
// ------------------------------------------------------------------------------------------------
// Summarize Nsk change as a result of tweaking systematic parameters.
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

#include "Conventions.h"

void test_nsk_twk_systematics()
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

  gSystem->Load("libT2KNuOsc");

  using namespace std;
  using namespace t2kana;

  const double input_sinsq2theta23 = 0.00;
  const double input_dm23sq        = 2.4E-3;

  const int nexperiments = 5000;

  const int nsystematics = 13;
  string systematic[nsystematics] = {
	"sk_ccqe_eff",
	"sk_ccnonqe_eff",
	"sk_nc_eff",
	"sk_nuecc_eff",
	"sk_energy_scale",
	"nd280_eff",
	"flux_norm",
	"flux_shape",
	"xsec_ccqe",
	"xsec_cc1pi_to_ccqe",
	"xsec_ccoth_to_ccqe",
	"xsec_nc_to_ccqe",
	"nuclear_effect"
  };
  bool generate_random[nsystematics] = {
	false,  /* sk_ccqe_eff        */
	false,  /* sk_ccnonqe_eff     */
	false,  /* sk_nc_eff          */
	false,  /* sk_nuecc_eff       */
	false,  /* sk_energy_scale    */
	false,  /* nd280_eff          */
	false,  /* flux_norm          */
	true,   /* flux_shape         */
	false,  /* xsec_ccqe          */
	false,  /* xsec_cc1pi_to_ccqe */
	false,  /* xsec_ccoth_to_ccqe */
	false,  /* xsec_nc_to_ccqe    */
	false   /* nuclear_effect     */
  };
  bool skip[nsystematics] = {
	false,  /* sk_ccqe_eff        */
	false,  /* sk_ccnonqe_eff     */
	false,  /* sk_nc_eff          */
	false,  /* sk_nuecc_eff       */
	true,   /* sk_energy_scale    (was shown to maintain Nsk) */
	false,  /* nd280_eff          */
	false,  /* flux_norm          */
	true,   /* flux_shape         (was shown to maintain Nsk)*/
	false,  /* xsec_ccqe          */
	false,  /* xsec_cc1pi_to_ccqe */
	false,  /* xsec_ccoth_to_ccqe */
	false,  /* xsec_nc_to_ccqe    */
	false   /* nuclear_effect     */
  };

  double nsk_min[nsystematics];
  double nsk_max[nsystematics];

  // osc params
  PhysParamList * fparm = PhysParamList::NumuDisapFitNoSyst();
  fparm->Update(t2kana::conv::kPhypSinsq2Theta23, input_sinsq2theta23);
  fparm->Update(t2kana::conv::kPhypDm23sq,        input_dm23sq);

  // get superk spectrum predictor & configure
  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read("inputs/analysis_configuration.xml");

  // get default prediction
  skexp->Generate(*fparm);
  double nsk_deflt = skexp->CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike)->Integral();
  
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
            for(int iexperiment = 0; iexperiment < nexperiments; iexperiment++) {
                    if(iexperiment%100==0) {
			cout << "Current experiment: " << iexperiment << endl;
                    }
		    skexp->ResetSystematicVariations();
	     	    skexp->TweakSystematic(systematic[isystematic]);
	  	    skexp->Generate(*fparm);
  		    double nsk = skexp->CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike)->Integral();
	            nsk_min[isystematic] = TMath::Min(nsk, nsk_min[isystematic]);
        	    nsk_max[isystematic] = TMath::Max(nsk, nsk_max[isystematic]);
            }//iexperiment
	} else {
	    skexp->ResetSystematicVariations();
     	    skexp->TweakSystematic(systematic[isystematic], +1);
  	    skexp->Generate(*fparm);
  	    double nsk1 = skexp->CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike)->Integral();
	    skexp->ResetSystematicVariations();
     	    skexp->TweakSystematic(systematic[isystematic], -1);
  	    skexp->Generate(*fparm);
  	    double nsk2 = skexp->CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike)->Integral();
            nsk_min[isystematic] = TMath::Min(nsk1,nsk2);
            nsk_max[isystematic] = TMath::Max(nsk1,nsk2);
	}	
  }//isystematic

  cout << endl;
  cout << "Nsk(default) = " << nsk_deflt << endl;
  for(int isystematic=0; isystematic<nsystematics; isystematic++) {
	cout << " * " << systematic[isystematic] << " : Nsk = (" 
	     << nsk_min[isystematic] << ", " << nsk_max[isystematic] << ")" << endl;
  }
  cout << endl;
  cout << "Fractional change: " << endl;
  for(int isystematic=0; isystematic<nsystematics; isystematic++) {
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
