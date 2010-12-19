//
// ------------------------------------------------------------------------------------------------
// Summarize Nsk change as a result of tweaking systematic parameters.
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

#include "Conventions.h"

void summarize_nsk_change()
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

  const int nexperiments = 100;

  const int nsystematics = 10;
  string systematic[nsystematics] = {
	"nd280_eff",
	"sk_ccqe_eff",
	"sk_ccnonqe_eff",
	"sk_nc_eff",
	"flux_norm",
	"xsec_ccqe",
	"xsec_cc1pi_to_ccqe",
	"xsec_ccoth_to_ccqe",
	"xsec_nc_to_ccqe",
	"nuclear_effect"
  };
  bool generate_random[nsystematics] = {
	false,
	false,
	false,
	false,
	false,
	false,
	false,
	false,
	false,
	false
  };

  double nsk_min[nsystematics];
  double nsk_max[nsystematics];

  // osc params
  PhysParamList * fparm = PhysParamList::NumuDisapFitNoSyst();
  fparm->Update(t2kana::yr1numu::conv::kPhypSinsq2Theta23, input_sinsq2theta23);
  fparm->Update(t2kana::yr1numu::conv::kPhypDm23sq,        input_dm23sq);

  // get superk spectrum predictor & configure
  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read("inputs/analysis_configuration.xml");

  // get default prediction
  skexp->Generate(*fparm);
  double nsk_deflt = skexp->CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike)->Integral();
  
  // get predictions with tweaked systematic params
  for(int isystematic=0; isystematic<nsystematics; isystematic++) {

	if(generate_random[isystematic]) {

            nsk_min[isystematic] = 99999999.;
            nsk_max[isystematic] = 0.;
            for(int iexperiment = 0; iexperiment < nexperiments; iexperiment++) {
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
             << 100*(nsk_min[isystematic]-nsk_deflt)/nsk_deflt << " %, " 
             << 100*(nsk_max[isystematic]-nsk_deflt)/nsk_deflt << " %) " << endl;
  }
  cout << endl;
  cout << "Total factional change: " << endl;
  double tot_pos = 0;
  double tot_neg = 0;
  for(int isystematic=0; isystematic<nsystematics; isystematic++) {
     double pos = nsk_max[isystematic]-nsk_deflt;
     double neg = nsk_min[isystematic]-nsk_deflt;
     tot_pos += (pos*pos);
     tot_neg += (neg*neg);
  }
  tot_pos = TMath::Sqrt(tot_pos);
  tot_neg = TMath::Sqrt(tot_neg);
  cout << "(-" << 100*tot_neg/nsk_deflt << "%, +" << 100*tot_pos/nsk_deflt << ")" << endl;
}
