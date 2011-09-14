//
// ------------------------------------------------------------------------------------------------
// Test systematic errors.
//
// The input (systematic) string can be:
//	- sk_energy_scale 
//	- sk_ccqe_other_eff
//	- sk_ccqe_ring_counting_eff : 
//             This energy-dependent correlated error can not be used in this test macro. Instead,
//             you can tweak specific components of the efficiency tweaking vector (neglecting,
//             correlation which, in the fitter, are handled using the Cholesky method). Use:
//	       - sk_ccqe_ring_counting_eff_0 
//	       - sk_ccqe_ring_counting_eff_1
//	       - sk_ccqe_ring_counting_eff_2
//	- sk_ccnonqe_eff
//	- sk_nc_eff
//	- sk_nuecc_eff
// 	- nd280_eff
// 	- nd280_norm
//	- flux_norm
//	- flux_shape : Can not be used in this test macro
//      - xsec_ccqe
//      - xsec_cc1pi_to_ccqe
//      - xsec_ccoth_to_ccqe
//      - xsec_nc_to_ccqe
//      - nuclear_effect
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

#include <iomanip.h>

#include "Conventions.h"

void test_systematic_error_propagation
(
  const char * systematic,                                                  /* options: see above */
  const char * hypothesis                 = "3f_osc",                       /* options: "2f_numu_disap", "3f_osc" */
  double       input_sinsq2theta12        = 0.86,                           /* */
  double       input_sinsq2theta13        = 0.10,                           /* */
  double       input_sinsq2theta23        = 1.00,                           /* */
  double       input_deltacp              = 0.00,                           /* */
  double       input_dm12sq               = 7.58E-5,                        /* */
  double       input_dm23sq               = 2.32E-3,                        /* */
  bool         use_tex                    = false,                          /* print-out results as a TeX table (for inclusion in tech. note)*/
  bool         show_sk_ccqe_spectral_info = false,                          /* */
  const char * config_file                = "analysis_configuration.xml"    /* */
)
{
  using namespace std;
  using namespace t2kana;

  // set physics hypothesis
  PhysHypothesis ph(hypothesis);
  if(ph.Hypothesis() == conv::kHypNuMuDisapp2F) {
    ph.SetCurrValue(conv::kPhypSinsq2Theta23,  input_sinsq2theta23);
    ph.SetCurrValue(conv::kPhypDm23sq,         input_dm23sq);
  } else
  if(ph.Hypothesis() == conv::kHypOsc3F) {
    ph.SetCurrValue(conv::kPhypSinsq2Theta12,  input_sinsq2theta12);
    ph.SetCurrValue(conv::kPhypSinsq2Theta13,  input_sinsq2theta13);
    ph.SetCurrValue(conv::kPhypSinsq2Theta23,  input_sinsq2theta23);
    ph.SetCurrValue(conv::kPhypDeltaCP,        input_deltacp);
    ph.SetCurrValue(conv::kPhypDm12sq,         input_dm12sq);
    ph.SetCurrValue(conv::kPhypDm23sq,         input_dm23sq);
  }
  cout << "Physics hypothesis set: " << endl;
  cout << ph;

  // create & config superk spectrum predictor
  const char * config = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);
  cout << "Reading analysis configuration from: " << config << endl;
  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read(config);


  double alpha = skexp->GetNuXSec()->OnePionInRES();
  cout << "alpha = " << alpha << endl;

  int ErecoBin1 =  9; // 0.4 GeV (incl.)
  int ErecoBin2 = 22; // 1.1 GeV (incl.)

  //
  // nominal
  //

  skexp->Generate(ph);

  double nsk_0                 = skexp->CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike)->Integral();

  double nsknumuccqe_0         = skexp->GetSuperkMC()->CurrErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuMuCCQE      )->Integral();
  double nsknumuccqe_0a        = skexp->GetSuperkMC()->CurrErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuMuCCQE      )->Integral(1,           ErecoBin1-1 );
  double nsknumuccqe_0b        = skexp->GetSuperkMC()->CurrErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuMuCCQE      )->Integral(ErecoBin1,   ErecoBin2   );
  double nsknumuccqe_0c        = skexp->GetSuperkMC()->CurrErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuMuCCQE      )->Integral(ErecoBin2+1, conv::kNvErecoBins);
  double nsknumucc1pi_0        = skexp->GetSuperkMC()->CurrErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuMuCC1pi     )->Integral();
  double nsknumuccoth_0        = skexp->GetSuperkMC()->CurrErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuMuCCoth     )->Integral();
  double nsknumucc_0           = nsknumuccqe_0 + nsknumucc1pi_0 + nsknumuccoth_0;
  double nsknumunc_0           = skexp->GetSuperkMC()->CurrErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNC            )->Integral();
  double nsknumubarccqe_0      = skexp->GetSuperkMC()->CurrErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuMuBarCCQE   )->Integral();
  double nsknumubarccqe_0a     = skexp->GetSuperkMC()->CurrErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuMuBarCCQE   )->Integral(1,           ErecoBin1-1 );
  double nsknumubarccqe_0b     = skexp->GetSuperkMC()->CurrErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuMuBarCCQE   )->Integral(ErecoBin1,   ErecoBin2   );
  double nsknumubarccqe_0c     = skexp->GetSuperkMC()->CurrErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuMuBarCCQE   )->Integral(ErecoBin2+1, conv::kNvErecoBins);
  double nsknumubarcc1pi_0     = skexp->GetSuperkMC()->CurrErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuMuBarCC1pi  )->Integral();
  double nsknumubarccoth_0     = skexp->GetSuperkMC()->CurrErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuMuBarCCoth  )->Integral();
  double nsknumubarcc_0        = nsknumubarccqe_0 + nsknumubarcc1pi_0 + nsknumubarccoth_0;
  double nsknumubarnc_0        = skexp->GetSuperkMC()->CurrErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuMuBarNC     )->Integral();
  double nsknueccqe_0          = skexp->GetSuperkMC()->CurrErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuECCQE       )->Integral();
  double nsknuecc1pi_0         = skexp->GetSuperkMC()->CurrErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuECC1pi      )->Integral();
  double nsknueccoth_0         = skexp->GetSuperkMC()->CurrErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuECCoth      )->Integral();
  double nsknuecc_0            = nsknueccqe_0 + nsknuecc1pi_0 + nsknueccoth_0;
  double nsknuenc_0            = skexp->GetSuperkMC()->CurrErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuENC         )->Integral();
  double nsksignueccqe_0       = skexp->GetSuperkMC()->CurrErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmSigNuECCQE    )->Integral();
  double nsksignuecc1pi_0      = skexp->GetSuperkMC()->CurrErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmSigNuECC1pi   )->Integral();
  double nsksignueccoth_0      = skexp->GetSuperkMC()->CurrErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmSigNuECCoth   )->Integral();
  double nsksignuecc_0         = nsksignueccqe_0 + nsksignuecc1pi_0 + nsksignueccoth_0;
  double nskall_0              = nsknumucc_0 + nsknumunc_0 + nsknumubarcc_0 + nsknumubarnc_0 + nsknuecc_0 + nsknuenc_0 + nsksignuecc_0;

  double nd280obs_0            = skexp->GetND280Measurement()->NObsCCLike();
  double nd280cc_0             = skexp->GetND280Measurement()->NExpCC();
  double nd280ncbkg_0          = skexp->GetND280Measurement()->NExpNCBkg();
  double nd280mccclike_0       = nd280cc_0 + nd280ncbkg_0;
  double nd280corr_0           = skexp->GetND280Measurement()->CorrectionFactor();

  double nd280numuccqe_0       =             skexp->GetND280Measurement()->NumOfEvents(conv::kNDmNuMuCCQE);
  double nd280numucc1pi_0      =     alpha * skexp->GetND280Measurement()->NumOfEvents(conv::kNDmNuMuCCRES);
  double nd280numuccoth_0      =  (1-alpha)* skexp->GetND280Measurement()->NumOfEvents(conv::kNDmNuMuCCRES) +   
                                             skexp->GetND280Measurement()->NumOfEvents(conv::kNDmNuMuCCCOH) +
                                             skexp->GetND280Measurement()->NumOfEvents(conv::kNDmNuMuCCOth);
  double nd280numunc_0         =             skexp->GetND280Measurement()->NumOfEvents(conv::kNDmNuMuNCRES) +
                                             skexp->GetND280Measurement()->NumOfEvents(conv::kNDmNuMuNCCOH) +
                                             skexp->GetND280Measurement()->NumOfEvents(conv::kNDmNuMuNCOth);
  double nd280numucc_0         = nd280numuccqe_0 + nd280numucc1pi_0 + nd280numuccoth_0; 
  double nd280numu_0           = nd280numucc_0 + nd280numunc_0;

  double nd280numubarccqe_0    =             skexp->GetND280Measurement()->NumOfEvents(conv::kNDmNuMuBarCCQE);
  double nd280numubarcc1pi_0   =     alpha * skexp->GetND280Measurement()->NumOfEvents(conv::kNDmNuMuBarCCRES);
  double nd280numubarccoth_0   =  (1-alpha)* skexp->GetND280Measurement()->NumOfEvents(conv::kNDmNuMuBarCCRES) +
                                             skexp->GetND280Measurement()->NumOfEvents(conv::kNDmNuMuBarCCCOH) +
                                             skexp->GetND280Measurement()->NumOfEvents(conv::kNDmNuMuBarCCOth);
  double nd280numubarnc_0      =             skexp->GetND280Measurement()->NumOfEvents(conv::kNDmNuMuBarNCRES) +
                                             skexp->GetND280Measurement()->NumOfEvents(conv::kNDmNuMuBarNCCOH) +
                                             skexp->GetND280Measurement()->NumOfEvents(conv::kNDmNuMuBarNCOth);
  double nd280numubarcc_0      = nd280numubarccqe_0 + nd280numubarcc1pi_0 + nd280numubarccoth_0;
  double nd280numubar_0        = nd280numubarcc_0 + nd280numubarnc_0;

  double nd280nueccqe_0        =            skexp->GetND280Measurement()->NumOfEvents(conv::kNDmNuECCQE);
  double nd280nuecc1pi_0       =    alpha * skexp->GetND280Measurement()->NumOfEvents(conv::kNDmNuECCRES);
  double nd280nueccoth_0       = (1-alpha)* skexp->GetND280Measurement()->NumOfEvents(conv::kNDmNuECCRES)  +
                                            skexp->GetND280Measurement()->NumOfEvents(conv::kNDmNuECCCOH)  +
                                            skexp->GetND280Measurement()->NumOfEvents(conv::kNDmNuECCOth);
  double nd280nuenc_0          =            skexp->GetND280Measurement()->NumOfEvents(conv::kNDmNuENCRES)  +
                                            skexp->GetND280Measurement()->NumOfEvents(conv::kNDmNuENCCOH)  +
                                            skexp->GetND280Measurement()->NumOfEvents(conv::kNDmNuENCOth);
  double nd280nuecc_0          = nd280nueccqe_0 + nd280nuecc1pi_0 + nd280nueccoth_0;
  double nd280nue_0            = nd280nuecc_0 + nd280nuenc_0;


  //
  // +1 sigma
  //

  skexp->TweakSystematic(systematic,+1.);

  skexp->Generate(ph);


  double nsk_p                 = skexp->CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike)->Integral();

  double nsknumuccqe_p         = skexp->GetSuperkMC()->CurrErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuMuCCQE      )->Integral();
  double nsknumuccqe_pa        = skexp->GetSuperkMC()->CurrErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuMuCCQE      )->Integral(1,           ErecoBin1-1 );
  double nsknumuccqe_pb        = skexp->GetSuperkMC()->CurrErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuMuCCQE      )->Integral(ErecoBin1,   ErecoBin2   );
  double nsknumuccqe_pc        = skexp->GetSuperkMC()->CurrErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuMuCCQE      )->Integral(ErecoBin2+1, conv::kNvErecoBins);
  double nsknumucc1pi_p        = skexp->GetSuperkMC()->CurrErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuMuCC1pi     )->Integral();
  double nsknumuccoth_p        = skexp->GetSuperkMC()->CurrErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuMuCCoth     )->Integral();
  double nsknumucc_p           = nsknumuccqe_p + nsknumucc1pi_p + nsknumuccoth_p;
  double nsknumunc_p           = skexp->GetSuperkMC()->CurrErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNC            )->Integral();
  double nsknumubarccqe_p      = skexp->GetSuperkMC()->CurrErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuMuBarCCQE   )->Integral();
  double nsknumubarccqe_pa     = skexp->GetSuperkMC()->CurrErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuMuBarCCQE   )->Integral(1,           ErecoBin1-1 );
  double nsknumubarccqe_pb     = skexp->GetSuperkMC()->CurrErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuMuBarCCQE   )->Integral(ErecoBin1,   ErecoBin2   );
  double nsknumubarccqe_pc     = skexp->GetSuperkMC()->CurrErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuMuBarCCQE   )->Integral(ErecoBin2+1, conv::kNvErecoBins);
  double nsknumubarcc1pi_p     = skexp->GetSuperkMC()->CurrErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuMuBarCC1pi  )->Integral();
  double nsknumubarccoth_p     = skexp->GetSuperkMC()->CurrErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuMuBarCCoth  )->Integral();
  double nsknumubarcc_p        = nsknumubarccqe_p + nsknumubarcc1pi_p + nsknumubarccoth_p;
  double nsknumubarnc_p        = skexp->GetSuperkMC()->CurrErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuMuBarNC     )->Integral();
  double nsknueccqe_p          = skexp->GetSuperkMC()->CurrErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuECCQE       )->Integral();
  double nsknuecc1pi_p         = skexp->GetSuperkMC()->CurrErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuECC1pi      )->Integral();
  double nsknueccoth_p         = skexp->GetSuperkMC()->CurrErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuECCoth      )->Integral();
  double nsknuecc_p            = nsknueccqe_p + nsknuecc1pi_p + nsknueccoth_p;
  double nsknuenc_p            = skexp->GetSuperkMC()->CurrErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuENC         )->Integral();
  double nsksignueccqe_p       = skexp->GetSuperkMC()->CurrErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmSigNuECCQE    )->Integral();
  double nsksignuecc1pi_p      = skexp->GetSuperkMC()->CurrErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmSigNuECC1pi   )->Integral();
  double nsksignueccoth_p      = skexp->GetSuperkMC()->CurrErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmSigNuECCoth   )->Integral();
  double nsksignuecc_p         = nsksignueccqe_p + nsksignuecc1pi_p + nsksignueccoth_p;
  double nskall_p              = nsknumucc_p + nsknumunc_p + nsknumubarcc_p + nsknumubarnc_p + nsknuecc_p + nsknuenc_p + nsksignuecc_p;

  double nd280obs_p            = skexp->GetND280Measurement()->NObsCCLike();
  double nd280cc_p             = skexp->GetND280Measurement()->NExpCC();
  double nd280ncbkg_p          = skexp->GetND280Measurement()->NExpNCBkg();
  double nd280mccclike_p       = nd280cc_p + nd280ncbkg_p;
  double nd280corr_p           = skexp->GetND280Measurement()->CorrectionFactor();

  double nd280numuccqe_p       =             skexp->GetND280Measurement()->NumOfEvents(conv::kNDmNuMuCCQE);
  double nd280numucc1pi_p      =     alpha * skexp->GetND280Measurement()->NumOfEvents(conv::kNDmNuMuCCRES);
  double nd280numuccoth_p      =  (1-alpha)* skexp->GetND280Measurement()->NumOfEvents(conv::kNDmNuMuCCRES) +   
                                             skexp->GetND280Measurement()->NumOfEvents(conv::kNDmNuMuCCCOH) +
                                             skexp->GetND280Measurement()->NumOfEvents(conv::kNDmNuMuCCOth);
  double nd280numunc_p         =             skexp->GetND280Measurement()->NumOfEvents(conv::kNDmNuMuNCRES) +
                                             skexp->GetND280Measurement()->NumOfEvents(conv::kNDmNuMuNCCOH) +
                                             skexp->GetND280Measurement()->NumOfEvents(conv::kNDmNuMuNCOth);
  double nd280numucc_p         = nd280numuccqe_p + nd280numucc1pi_p + nd280numuccoth_p; 
  double nd280numu_p           = nd280numucc_p + nd280numunc_p;

  double nd280numubarccqe_p    =             skexp->GetND280Measurement()->NumOfEvents(conv::kNDmNuMuBarCCQE);
  double nd280numubarcc1pi_p   =     alpha * skexp->GetND280Measurement()->NumOfEvents(conv::kNDmNuMuBarCCRES);
  double nd280numubarccoth_p   =  (1-alpha)* skexp->GetND280Measurement()->NumOfEvents(conv::kNDmNuMuBarCCRES) +
                                             skexp->GetND280Measurement()->NumOfEvents(conv::kNDmNuMuBarCCCOH) +
                                             skexp->GetND280Measurement()->NumOfEvents(conv::kNDmNuMuBarCCOth);
  double nd280numubarnc_p      =             skexp->GetND280Measurement()->NumOfEvents(conv::kNDmNuMuBarNCRES) +
                                             skexp->GetND280Measurement()->NumOfEvents(conv::kNDmNuMuBarNCCOH) +
                                             skexp->GetND280Measurement()->NumOfEvents(conv::kNDmNuMuBarNCOth);
  double nd280numubarcc_p      = nd280numubarccqe_p + nd280numubarcc1pi_p + nd280numubarccoth_p;
  double nd280numubar_p        = nd280numubarcc_p + nd280numubarnc_p;

  double nd280nueccqe_p        =            skexp->GetND280Measurement()->NumOfEvents(conv::kNDmNuECCQE);
  double nd280nuecc1pi_p       =    alpha * skexp->GetND280Measurement()->NumOfEvents(conv::kNDmNuECCRES);
  double nd280nueccoth_p       = (1-alpha)* skexp->GetND280Measurement()->NumOfEvents(conv::kNDmNuECCRES)  +
                                            skexp->GetND280Measurement()->NumOfEvents(conv::kNDmNuECCCOH)  +
                                            skexp->GetND280Measurement()->NumOfEvents(conv::kNDmNuECCOth);
  double nd280nuenc_p          =            skexp->GetND280Measurement()->NumOfEvents(conv::kNDmNuENCRES)  +
                                            skexp->GetND280Measurement()->NumOfEvents(conv::kNDmNuENCCOH)  +
                                            skexp->GetND280Measurement()->NumOfEvents(conv::kNDmNuENCOth);
  double nd280nuecc_p          = nd280nueccqe_p + nd280nuecc1pi_p + nd280nueccoth_p;
  double nd280nue_p            = nd280nuecc_p + nd280nuenc_p;

  //
  // -1 sigma
  //

  skexp->TweakSystematic(systematic,-1.);

  skexp->Generate(ph);

  double nsk_m                 = skexp->CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike)->Integral();

  double nsknumuccqe_m         = skexp->GetSuperkMC()->CurrErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuMuCCQE      )->Integral();
  double nsknumuccqe_ma        = skexp->GetSuperkMC()->CurrErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuMuCCQE      )->Integral(1,           ErecoBin1-1 );
  double nsknumuccqe_mb        = skexp->GetSuperkMC()->CurrErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuMuCCQE      )->Integral(ErecoBin1,   ErecoBin2   );
  double nsknumuccqe_mc        = skexp->GetSuperkMC()->CurrErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuMuCCQE      )->Integral(ErecoBin2+1, conv::kNvErecoBins);
  double nsknumucc1pi_m        = skexp->GetSuperkMC()->CurrErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuMuCC1pi     )->Integral();
  double nsknumuccoth_m        = skexp->GetSuperkMC()->CurrErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuMuCCoth     )->Integral();
  double nsknumucc_m           = nsknumuccqe_m + nsknumucc1pi_m + nsknumuccoth_m;
  double nsknumunc_m           = skexp->GetSuperkMC()->CurrErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNC            )->Integral();
  double nsknumubarccqe_m      = skexp->GetSuperkMC()->CurrErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuMuBarCCQE   )->Integral();
  double nsknumubarccqe_ma     = skexp->GetSuperkMC()->CurrErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuMuBarCCQE   )->Integral(1,           ErecoBin1-1 );
  double nsknumubarccqe_mb     = skexp->GetSuperkMC()->CurrErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuMuBarCCQE   )->Integral(ErecoBin1,   ErecoBin2   );
  double nsknumubarccqe_mc     = skexp->GetSuperkMC()->CurrErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuMuBarCCQE   )->Integral(ErecoBin2+1, conv::kNvErecoBins);
  double nsknumubarcc1pi_m     = skexp->GetSuperkMC()->CurrErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuMuBarCC1pi  )->Integral();
  double nsknumubarccoth_m     = skexp->GetSuperkMC()->CurrErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuMuBarCCoth  )->Integral();
  double nsknumubarcc_m        = nsknumubarccqe_m + nsknumubarcc1pi_m + nsknumubarccoth_m;
  double nsknumubarnc_m        = skexp->GetSuperkMC()->CurrErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuMuBarNC     )->Integral();
  double nsknueccqe_m          = skexp->GetSuperkMC()->CurrErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuECCQE       )->Integral();
  double nsknuecc1pi_m         = skexp->GetSuperkMC()->CurrErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuECC1pi      )->Integral();
  double nsknueccoth_m         = skexp->GetSuperkMC()->CurrErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuECCoth      )->Integral();
  double nsknuecc_m            = nsknueccqe_m + nsknuecc1pi_m + nsknueccoth_m;
  double nsknuenc_m            = skexp->GetSuperkMC()->CurrErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuENC         )->Integral();
  double nsksignueccqe_m       = skexp->GetSuperkMC()->CurrErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmSigNuECCQE    )->Integral();
  double nsksignuecc1pi_m      = skexp->GetSuperkMC()->CurrErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmSigNuECC1pi   )->Integral();
  double nsksignueccoth_m      = skexp->GetSuperkMC()->CurrErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmSigNuECCoth   )->Integral();
  double nsksignuecc_m         = nsksignueccqe_m + nsksignuecc1pi_m + nsksignueccoth_m;
  double nskall_m              = nsknumucc_m + nsknumunc_m + nsknumubarcc_m + nsknumubarnc_m + nsknuecc_m + nsknuenc_m + nsksignuecc_m;

  double nd280obs_m            = skexp->GetND280Measurement()->NObsCCLike();
  double nd280cc_m             = skexp->GetND280Measurement()->NExpCC();
  double nd280ncbkg_m          = skexp->GetND280Measurement()->NExpNCBkg();
  double nd280mccclike_m       = nd280cc_m + nd280ncbkg_m;
  double nd280corr_m           = skexp->GetND280Measurement()->CorrectionFactor();

  double nd280numuccqe_m       =             skexp->GetND280Measurement()->NumOfEvents(conv::kNDmNuMuCCQE);
  double nd280numucc1pi_m      =     alpha * skexp->GetND280Measurement()->NumOfEvents(conv::kNDmNuMuCCRES);
  double nd280numuccoth_m      =  (1-alpha)* skexp->GetND280Measurement()->NumOfEvents(conv::kNDmNuMuCCRES) +   
                                             skexp->GetND280Measurement()->NumOfEvents(conv::kNDmNuMuCCCOH) +
                                             skexp->GetND280Measurement()->NumOfEvents(conv::kNDmNuMuCCOth);
  double nd280numunc_m         =             skexp->GetND280Measurement()->NumOfEvents(conv::kNDmNuMuNCRES) +
                                             skexp->GetND280Measurement()->NumOfEvents(conv::kNDmNuMuNCCOH) +
                                             skexp->GetND280Measurement()->NumOfEvents(conv::kNDmNuMuNCOth);
  double nd280numucc_m         = nd280numuccqe_m + nd280numucc1pi_m + nd280numuccoth_m; 
  double nd280numu_m           = nd280numucc_m + nd280numunc_m;

  double nd280numubarccqe_m    =             skexp->GetND280Measurement()->NumOfEvents(conv::kNDmNuMuBarCCQE);
  double nd280numubarcc1pi_m   =     alpha * skexp->GetND280Measurement()->NumOfEvents(conv::kNDmNuMuBarCCRES);
  double nd280numubarccoth_m   =  (1-alpha)* skexp->GetND280Measurement()->NumOfEvents(conv::kNDmNuMuBarCCRES) +
                                             skexp->GetND280Measurement()->NumOfEvents(conv::kNDmNuMuBarCCCOH) +
                                             skexp->GetND280Measurement()->NumOfEvents(conv::kNDmNuMuBarCCOth);
  double nd280numubarnc_m      =             skexp->GetND280Measurement()->NumOfEvents(conv::kNDmNuMuBarNCRES) +
                                             skexp->GetND280Measurement()->NumOfEvents(conv::kNDmNuMuBarNCCOH) +
                                             skexp->GetND280Measurement()->NumOfEvents(conv::kNDmNuMuBarNCOth);
  double nd280numubarcc_m      = nd280numubarccqe_m + nd280numubarcc1pi_m + nd280numubarccoth_m;
  double nd280numubar_m        = nd280numubarcc_m + nd280numubarnc_m;

  double nd280nueccqe_m        =             skexp->GetND280Measurement()->NumOfEvents(conv::kNDmNuECCQE);
  double nd280nuecc1pi_m       =    alpha * skexp->GetND280Measurement()->NumOfEvents(conv::kNDmNuECCRES);
  double nd280nueccoth_m       = (1-alpha)* skexp->GetND280Measurement()->NumOfEvents(conv::kNDmNuECCRES)  +
                                            skexp->GetND280Measurement()->NumOfEvents(conv::kNDmNuECCCOH)  +
                                            skexp->GetND280Measurement()->NumOfEvents(conv::kNDmNuECCOth);
  double nd280nuenc_m          =            skexp->GetND280Measurement()->NumOfEvents(conv::kNDmNuENCRES)  +
                                            skexp->GetND280Measurement()->NumOfEvents(conv::kNDmNuENCCOH)  +
                                            skexp->GetND280Measurement()->NumOfEvents(conv::kNDmNuENCOth);
  double nd280nuecc_m          = nd280nueccqe_m + nd280nuecc1pi_m + nd280nueccoth_m;
  double nd280nue_m            = nd280nuecc_m + nd280nuenc_m;


  //
  // calc fractional changes
  //

  double frac_nsk_p                  = (nsk_0                  > 0) ? 100.*(nsk_p                 - nsk_0                ) / nsk_0                 : 0.;
  double frac_nskall_p               = (nskall_0               > 0) ? 100.*(nskall_p              - nskall_0             ) / nskall_0              : 0.;
  double frac_nsknumuccqe_p          = (nsknumuccqe_0          > 0) ? 100.*(nsknumuccqe_p         - nsknumuccqe_0        ) / nsknumuccqe_0         : 0.;
  double frac_nsknumuccqe_pa         = (nsknumuccqe_0a         > 0) ? 100.*(nsknumuccqe_pa        - nsknumuccqe_0a       ) / nsknumuccqe_0a        : 0.;
  double frac_nsknumuccqe_pb         = (nsknumuccqe_0b         > 0) ? 100.*(nsknumuccqe_pb        - nsknumuccqe_0b       ) / nsknumuccqe_0b        : 0.;
  double frac_nsknumuccqe_pc         = (nsknumuccqe_0c         > 0) ? 100.*(nsknumuccqe_pc        - nsknumuccqe_0c       ) / nsknumuccqe_0c        : 0.;
  double frac_nsknumucc1pi_p         = (nsknumucc1pi_0         > 0) ? 100.*(nsknumucc1pi_p        - nsknumucc1pi_0       ) / nsknumucc1pi_0        : 0.;
  double frac_nsknumuccoth_p         = (nsknumuccoth_0         > 0) ? 100.*(nsknumuccoth_p        - nsknumuccoth_0       ) / nsknumuccoth_0        : 0.;
  double frac_nsknumucc_p            = (nsknumucc_0            > 0) ? 100.*(nsknumucc_p           - nsknumucc_0          ) / nsknumucc_0           : 0.;
  double frac_nsknumunc_p            = (nsknumunc_0            > 0) ? 100.*(nsknumunc_p           - nsknumunc_0          ) / nsknumunc_0           : 0.;
  double frac_nsknumubarccqe_p       = (nsknumubarccqe_0       > 0) ? 100.*(nsknumubarccqe_p      - nsknumubarccqe_0     ) / nsknumubarccqe_0      : 0.;
  double frac_nsknumubarccqe_pa      = (nsknumubarccqe_0a      > 0) ? 100.*(nsknumubarccqe_pa     - nsknumubarccqe_0a    ) / nsknumubarccqe_0a     : 0.;
  double frac_nsknumubarccqe_pb      = (nsknumubarccqe_0b      > 0) ? 100.*(nsknumubarccqe_pb     - nsknumubarccqe_0b    ) / nsknumubarccqe_0b     : 0.;
  double frac_nsknumubarccqe_pc      = (nsknumubarccqe_0c      > 0) ? 100.*(nsknumubarccqe_pc     - nsknumubarccqe_0c    ) / nsknumubarccqe_0c     : 0.;
  double frac_nsknumubarcc1pi_p      = (nsknumubarcc1pi_0      > 0) ? 100.*(nsknumubarcc1pi_p     - nsknumubarcc1pi_0    ) / nsknumubarcc1pi_0     : 0.;
  double frac_nsknumubarccoth_p      = (nsknumubarccoth_0      > 0) ? 100.*(nsknumubarccoth_p     - nsknumubarccoth_0    ) / nsknumubarccoth_0     : 0.;
  double frac_nsknumubarcc_p         = (nsknumubarcc_0         > 0) ? 100.*(nsknumubarcc_p        - nsknumubarcc_0       ) / nsknumubarcc_0        : 0.;
  double frac_nsknumubarnc_p         = (nsknumubarnc_0         > 0) ? 100.*(nsknumubarnc_p        - nsknumubarnc_0       ) / nsknumubarnc_0        : 0.;
  double frac_nsknueccqe_p           = (nsknueccqe_0           > 0) ? 100.*(nsknueccqe_p          - nsknueccqe_0         ) / nsknueccqe_0          : 0.;
  double frac_nsknuecc1pi_p          = (nsknuecc1pi_0          > 0) ? 100.*(nsknuecc1pi_p         - nsknuecc1pi_0        ) / nsknuecc1pi_0         : 0.;
  double frac_nsknueccoth_p          = (nsknueccoth_0          > 0) ? 100.*(nsknueccoth_p         - nsknueccoth_0        ) / nsknueccoth_0         : 0.;
  double frac_nsknuecc_p             = (nsknuecc_0             > 0) ? 100.*(nsknuecc_p            - nsknuecc_0           ) / nsknuecc_0            : 0.;
  double frac_nsknuenc_p             = (nsknuenc_0             > 0) ? 100.*(nsknuenc_p            - nsknuenc_0           ) / nsknuenc_0            : 0.;
  double frac_nsksignueccqe_p        = (nsksignueccqe_0        > 0) ? 100.*(nsksignueccqe_p       - nsksignueccqe_0      ) / nsksignueccqe_0       : 0.;
  double frac_nsksignuecc1pi_p       = (nsksignuecc1pi_0       > 0) ? 100.*(nsksignuecc1pi_p      - nsksignuecc1pi_0     ) / nsksignuecc1pi_0      : 0.;
  double frac_nsksignueccoth_p       = (nsksignueccoth_0       > 0) ? 100.*(nsksignueccoth_p      - nsksignueccoth_0     ) / nsksignueccoth_0      : 0.;
  double frac_nsksignuecc_p          = (nsksignuecc_0          > 0) ? 100.*(nsksignuecc_p         - nsksignuecc_0        ) / nsksignuecc_0         : 0.;
  double frac_nd280obs_p             = (nd280obs_0             > 0) ? 100.*(nd280obs_p            - nd280obs_0           ) / nd280obs_0            : 0.;
  double frac_nd280mccclike_p        = (nd280mccclike_0        > 0) ? 100.*(nd280mccclike_p       - nd280mccclike_0      ) / nd280mccclike_0       : 0.;
  double frac_nd280cc_p              = (nd280cc_0              > 0) ? 100.*(nd280cc_p             - nd280cc_0            ) / nd280cc_0             : 0.;
  double frac_nd280ncbkg_p           = (nd280ncbkg_0           > 0) ? 100.*(nd280ncbkg_p          - nd280ncbkg_0         ) / nd280ncbkg_0          : 0.;
  double frac_nd280corr_p            = (nd280corr_0            > 0) ? 100.*(nd280corr_p           - nd280corr_0          ) / nd280corr_0           : 0.;
  double frac_nd280numuccqe_p        = (nd280numuccqe_0        > 0) ? 100.*(nd280numuccqe_p       - nd280numuccqe_0      ) / nd280numuccqe_0       : 0.;
  double frac_nd280numucc1pi_p       = (nd280numucc1pi_0       > 0) ? 100.*(nd280numucc1pi_p      - nd280numucc1pi_0     ) / nd280numucc1pi_0      : 0.;
  double frac_nd280numuccoth_p       = (nd280numuccoth_0       > 0) ? 100.*(nd280numuccoth_p      - nd280numuccoth_0     ) / nd280numuccoth_0      : 0.;
  double frac_nd280numucc_p          = (nd280numucc_0          > 0) ? 100.*(nd280numucc_p         - nd280numucc_0        ) / nd280numucc_0         : 0.;
  double frac_nd280numunc_p          = (nd280numunc_0          > 0) ? 100.*(nd280numunc_p         - nd280numunc_0        ) / nd280numunc_0         : 0.;
  double frac_nd280numubarccqe_p     = (nd280numubarccqe_0     > 0) ? 100.*(nd280numubarccqe_p    - nd280numubarccqe_0   ) / nd280numubarccqe_0    : 0.;
  double frac_nd280numubarcc1pi_p    = (nd280numubarcc1pi_0    > 0) ? 100.*(nd280numubarcc1pi_p   - nd280numubarcc1pi_0  ) / nd280numubarcc1pi_0   : 0.;
  double frac_nd280numubarccoth_p    = (nd280numubarccoth_0    > 0) ? 100.*(nd280numubarccoth_p   - nd280numubarccoth_0  ) / nd280numubarccoth_0   : 0.;
  double frac_nd280numubarcc_p       = (nd280numubarcc_0       > 0) ? 100.*(nd280numubarcc_p      - nd280numubarcc_0     ) / nd280numubarcc_0      : 0.;
  double frac_nd280numubarnc_p       = (nd280numubarnc_0       > 0) ? 100.*(nd280numubarnc_p      - nd280numubarnc_0     ) / nd280numubarnc_0      : 0.;
  double frac_nd280nueccqe_p         = (nd280nueccqe_0         > 0) ? 100.*(nd280nueccqe_p        - nd280nueccqe_0       ) / nd280nueccqe_0        : 0.;
  double frac_nd280nuecc1pi_p        = (nd280nuecc1pi_0        > 0) ? 100.*(nd280nuecc1pi_p       - nd280nuecc1pi_0      ) / nd280nuecc1pi_0       : 0.;
  double frac_nd280nueccoth_p        = (nd280nueccoth_0        > 0) ? 100.*(nd280nueccoth_p       - nd280nueccoth_0      ) / nd280nueccoth_0       : 0.;
  double frac_nd280nuecc_p           = (nd280nuecc_0           > 0) ? 100.*(nd280nuecc_p          - nd280nuecc_0         ) / nd280nuecc_0          : 0.;
  double frac_nd280nuenc_p           = (nd280nuenc_0           > 0) ? 100.*(nd280nuenc_p          - nd280nuenc_0         ) / nd280nuenc_0          : 0.;

  double frac_nsk_m                  = (nsk_0                  > 0) ? 100.*(nsk_m                 - nsk_0                ) / nsk_0                 : 0.;
  double frac_nskall_m               = (nskall_0               > 0) ? 100.*(nskall_m              - nskall_0             ) / nskall_0              : 0.;
  double frac_nsknumuccqe_m          = (nsknumuccqe_0          > 0) ? 100.*(nsknumuccqe_m         - nsknumuccqe_0        ) / nsknumuccqe_0         : 0.;
  double frac_nsknumuccqe_ma         = (nsknumuccqe_0a         > 0) ? 100.*(nsknumuccqe_ma        - nsknumuccqe_0a       ) / nsknumuccqe_0a        : 0.;
  double frac_nsknumuccqe_mb         = (nsknumuccqe_0b         > 0) ? 100.*(nsknumuccqe_mb        - nsknumuccqe_0b       ) / nsknumuccqe_0b        : 0.;
  double frac_nsknumuccqe_mc         = (nsknumuccqe_0c         > 0) ? 100.*(nsknumuccqe_mc        - nsknumuccqe_0c       ) / nsknumuccqe_0c        : 0.;
  double frac_nsknumucc1pi_m         = (nsknumucc1pi_0         > 0) ? 100.*(nsknumucc1pi_m        - nsknumucc1pi_0       ) / nsknumucc1pi_0        : 0.;
  double frac_nsknumuccoth_m         = (nsknumuccoth_0         > 0) ? 100.*(nsknumuccoth_m        - nsknumuccoth_0       ) / nsknumuccoth_0        : 0.;
  double frac_nsknumucc_m            = (nsknumucc_0            > 0) ? 100.*(nsknumucc_m           - nsknumucc_0          ) / nsknumucc_0           : 0.;
  double frac_nsknumunc_m            = (nsknumunc_0            > 0) ? 100.*(nsknumunc_m           - nsknumunc_0          ) / nsknumunc_0           : 0.;
  double frac_nsknumubarccqe_m       = (nsknumubarccqe_0       > 0) ? 100.*(nsknumubarccqe_m      - nsknumubarccqe_0     ) / nsknumubarccqe_0      : 0.;
  double frac_nsknumubarccqe_ma      = (nsknumubarccqe_0a      > 0) ? 100.*(nsknumubarccqe_ma     - nsknumubarccqe_0a    ) / nsknumubarccqe_0a     : 0.;
  double frac_nsknumubarccqe_mb      = (nsknumubarccqe_0b      > 0) ? 100.*(nsknumubarccqe_mb     - nsknumubarccqe_0b    ) / nsknumubarccqe_0b     : 0.;
  double frac_nsknumubarccqe_mc      = (nsknumubarccqe_0c      > 0) ? 100.*(nsknumubarccqe_mc     - nsknumubarccqe_0c    ) / nsknumubarccqe_0c     : 0.;
  double frac_nsknumubarcc1pi_m      = (nsknumubarcc1pi_0      > 0) ? 100.*(nsknumubarcc1pi_m     - nsknumubarcc1pi_0    ) / nsknumubarcc1pi_0     : 0.;
  double frac_nsknumubarccoth_m      = (nsknumubarccoth_0      > 0) ? 100.*(nsknumubarccoth_m     - nsknumubarccoth_0    ) / nsknumubarccoth_0     : 0.;
  double frac_nsknumubarcc_m         = (nsknumubarcc_0         > 0) ? 100.*(nsknumubarcc_m        - nsknumubarcc_0       ) / nsknumubarcc_0        : 0.;
  double frac_nsknumubarnc_m         = (nsknumubarnc_0         > 0) ? 100.*(nsknumubarnc_m        - nsknumubarnc_0       ) / nsknumubarnc_0        : 0.;
  double frac_nsknueccqe_m           = (nsknueccqe_0           > 0) ? 100.*(nsknueccqe_m          - nsknueccqe_0         ) / nsknueccqe_0          : 0.;
  double frac_nsknuecc1pi_m          = (nsknuecc1pi_0          > 0) ? 100.*(nsknuecc1pi_m         - nsknuecc1pi_0        ) / nsknuecc1pi_0         : 0.;
  double frac_nsknueccoth_m          = (nsknueccoth_0          > 0) ? 100.*(nsknueccoth_m         - nsknueccoth_0        ) / nsknueccoth_0         : 0.;
  double frac_nsknuecc_m             = (nsknuecc_0             > 0) ? 100.*(nsknuecc_m            - nsknuecc_0           ) / nsknuecc_0            : 0.;
  double frac_nsknuenc_m             = (nsknuenc_0             > 0) ? 100.*(nsknuenc_m            - nsknuenc_0           ) / nsknuenc_0            : 0.;
  double frac_nsksignueccqe_m        = (nsksignueccqe_0        > 0) ? 100.*(nsksignueccqe_m       - nsksignueccqe_0      ) / nsksignueccqe_0       : 0.;
  double frac_nsksignuecc1pi_m       = (nsksignuecc1pi_0       > 0) ? 100.*(nsksignuecc1pi_m      - nsksignuecc1pi_0     ) / nsksignuecc1pi_0      : 0.;
  double frac_nsksignueccoth_m       = (nsksignueccoth_0       > 0) ? 100.*(nsksignueccoth_m      - nsksignueccoth_0     ) / nsksignueccoth_0      : 0.;
  double frac_nsksignuecc_m          = (nsksignuecc_0          > 0) ? 100.*(nsksignuecc_m         - nsksignuecc_0        ) / nsksignuecc_0         : 0.;
  double frac_nd280obs_m             = (nd280obs_0             > 0) ? 100.*(nd280obs_m            - nd280obs_0           ) / nd280obs_0            : 0.;
  double frac_nd280mccclike_m        = (nd280mccclike_0        > 0) ? 100.*(nd280mccclike_m       - nd280mccclike_0      ) / nd280mccclike_0       : 0.;
  double frac_nd280cc_m              = (nd280cc_0              > 0) ? 100.*(nd280cc_m             - nd280cc_0            ) / nd280cc_0             : 0.;
  double frac_nd280ncbkg_m           = (nd280ncbkg_0           > 0) ? 100.*(nd280ncbkg_m          - nd280ncbkg_0         ) / nd280ncbkg_0          : 0.;
  double frac_nd280corr_m            = (nd280corr_0            > 0) ? 100.*(nd280corr_m           - nd280corr_0          ) / nd280corr_0           : 0.;
  double frac_nd280numuccqe_m        = (nd280numuccqe_0        > 0) ? 100.*(nd280numuccqe_m       - nd280numuccqe_0      ) / nd280numuccqe_0       : 0.;
  double frac_nd280numucc1pi_m       = (nd280numucc1pi_0       > 0) ? 100.*(nd280numucc1pi_m      - nd280numucc1pi_0     ) / nd280numucc1pi_0      : 0.;
  double frac_nd280numuccoth_m       = (nd280numuccoth_0       > 0) ? 100.*(nd280numuccoth_m      - nd280numuccoth_0     ) / nd280numuccoth_0      : 0.;
  double frac_nd280numucc_m          = (nd280numucc_0          > 0) ? 100.*(nd280numucc_m         - nd280numucc_0        ) / nd280numucc_0         : 0.;
  double frac_nd280numunc_m          = (nd280numunc_0          > 0) ? 100.*(nd280numunc_m         - nd280numunc_0        ) / nd280numunc_0         : 0.;
  double frac_nd280numubarccqe_m     = (nd280numubarccqe_0     > 0) ? 100.*(nd280numubarccqe_m    - nd280numubarccqe_0   ) / nd280numubarccqe_0    : 0.;
  double frac_nd280numubarcc1pi_m    = (nd280numubarcc1pi_0    > 0) ? 100.*(nd280numubarcc1pi_m   - nd280numubarcc1pi_0  ) / nd280numubarcc1pi_0   : 0.;
  double frac_nd280numubarccoth_m    = (nd280numubarccoth_0    > 0) ? 100.*(nd280numubarccoth_m   - nd280numubarccoth_0  ) / nd280numubarccoth_0   : 0.;
  double frac_nd280numubarcc_m       = (nd280numubarcc_0       > 0) ? 100.*(nd280numubarcc_m      - nd280numubarcc_0     ) / nd280numubarcc_0      : 0.;
  double frac_nd280numubarnc_m       = (nd280numubarnc_0       > 0) ? 100.*(nd280numubarnc_m      - nd280numubarnc_0     ) / nd280numubarnc_0      : 0.;
  double frac_nd280nueccqe_m         = (nd280nueccqe_0         > 0) ? 100.*(nd280nueccqe_m        - nd280nueccqe_0       ) / nd280nueccqe_0        : 0.;
  double frac_nd280nuecc1pi_m        = (nd280nuecc1pi_0        > 0) ? 100.*(nd280nuecc1pi_m       - nd280nuecc1pi_0      ) / nd280nuecc1pi_0       : 0.;
  double frac_nd280nueccoth_m        = (nd280nueccoth_0        > 0) ? 100.*(nd280nueccoth_m       - nd280nueccoth_0      ) / nd280nueccoth_0       : 0.;
  double frac_nd280nuecc_m           = (nd280nuecc_0           > 0) ? 100.*(nd280nuecc_m          - nd280nuecc_0         ) / nd280nuecc_0          : 0.;
  double frac_nd280nuenc_m           = (nd280nuenc_0           > 0) ? 100.*(nd280nuenc_m          - nd280nuenc_0         ) / nd280nuenc_0          : 0.;

  //
  // print-out
  //

  int    nprecisionnev = 4;
  int    nprecisionpct = 2;
  string pctsymbol  = (use_tex) ? " \\% "   : "% ";
  string endcolumn  = (use_tex) ? " & "     : "" ;
  string endrow     = (use_tex) ? " \\\\ "  : "" ;

  int ncol1   = 80;
  int ncolnev = 15;
  int ncolfrc =  8;

  cout << "\n\nSystematic: " << systematic << endl;
  cout << endl;
  cout << fixed;

  if(use_tex) {
     cout << "\\begin{table}"   << endl;
     cout << "\\centering"     << endl;
     cout << "\\begin{tabular}{| l | r | r  r | r  r |}" << endl;
     cout << "\\hline" << endl;
     cout << "{\\bf Component} &  "
          << "{\\bf Nominal}   &  "
          << "\\multicolumn {2}{c|} { {\\bf +1$\\sigma$} } &  "
          << "\\multicolumn {2}{c|} { {\\bf -1$\\sigma$} } \\\\"
          << endl;
     cout << "\\hline" << endl;
  }

  cout << setw(ncol1)  << setfill(' ') << "$N_{SK}$ (1 $\\mu$-like ring)" 
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsk_0 
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsk_p
       << endcolumn
       << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsk_p << pctsymbol 
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsk_m 
       << endcolumn
       << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsk_m << pctsymbol 
       << endrow
       << endl;

  if(use_tex) {
     cout << "\\hline" << endl;
  } else {
     cout << setw(150)  << setfill('-') << " " << endl;
  }

  cout << setw(ncol1)  << setfill(' ') << "SK MC (1 $\\mu$-like ring)"
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nskall_0
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nskall_p
       << endcolumn
       << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nskall_p << pctsymbol 
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nskall_m 
       << endcolumn
       << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nskall_m << pctsymbol 
       << endrow
       << endl;

  if(use_tex) {
     cout << "\\hline" << endl;
  } else {
     cout << setw(150)  << setfill('-') << " " << endl;
  }

  cout << setw(ncol1)  << setfill(' ') << "SK MC (1 $\\mu$-like ring): $\\nu_{\\mu}$ CCQE" 
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumuccqe_0
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumuccqe_p
       << endcolumn
       << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknumuccqe_p << pctsymbol 
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumuccqe_m 
       << endcolumn
       << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknumuccqe_m << pctsymbol 
       << endrow
       << endl;
  if(show_sk_ccqe_spectral_info) {
      cout << setw(ncol1)  << setfill(' ') << "SK MC (1 $\\mu$-like ring): $\\nu_{\\mu}$ CCQE  ($E_{r}$ < 0.4 GeV)"
           << endcolumn
           << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumuccqe_0a
           << endcolumn
           << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumuccqe_pa
           << endcolumn
           << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknumuccqe_pa << pctsymbol 
           << endcolumn
           << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumuccqe_ma 
           << endcolumn
           << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknumuccqe_ma << pctsymbol 
           << endrow
           << endl;
      cout << setw(ncol1)  << setfill(' ') << "SK MC (1 $\\mu$-like ring): $\\nu_{\\mu}$ CCQE  (0.4 <= $E_{r}$ <= 1.1 GeV)"
           << endcolumn
           << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumuccqe_0b
           << endcolumn
           << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumuccqe_pb
           << endcolumn
           << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknumuccqe_pb << pctsymbol 
           << endcolumn
           << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumuccqe_mb 
           << endcolumn
           << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknumuccqe_mb << pctsymbol 
           << endrow
           << endl;
      cout << setw(ncol1)  << setfill(' ') << "SK MC (1 $\\mu$-like ring): $\\nu_{\\mu}$ CCQE  ($E_{r}$ > 1.1 GeV)"
           << endcolumn
           << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumuccqe_0c
           << endcolumn
           << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumuccqe_pc
           << endcolumn
           << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknumuccqe_pc << pctsymbol 
           << endcolumn
           << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumuccqe_mc 
           << endcolumn
           << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknumuccqe_mc << pctsymbol 
           << endrow
           << endl;
  }
  cout << setw(ncol1)  << setfill(' ') << "SK MC (1 $\\mu$-like ring): $\\nu_{\\mu}$ CC1$\\pi$" 
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumucc1pi_0
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumucc1pi_p
       << endcolumn
       << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknumucc1pi_p << pctsymbol 
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumucc1pi_m 
       << endcolumn
       << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknumucc1pi_m << pctsymbol 
       << endrow
       << endl;
  cout << setw(ncol1)  << setfill(' ') << "SK MC (1 $\\mu$-like ring): $\\nu_{\\mu}$ CC other" 
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumuccoth_0
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumuccoth_p
       << endcolumn
       << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknumuccoth_p << pctsymbol 
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumuccoth_m 
       << endcolumn
       << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknumuccoth_m << pctsymbol 
       << endrow
       << endl;
  cout << setw(ncol1)  << setfill(' ') << "SK MC (1 $\\mu$-like ring): $\\nu_{\\mu}$ CC all" 
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumucc_0
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumucc_p
       << endcolumn
       << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknumucc_p << pctsymbol 
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumucc_m 
       << endcolumn
       << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknumucc_m << pctsymbol 
       << endrow
       << endl;
  cout << setw(ncol1)  << setfill(' ') << "SK MC (1 $\\mu$-like ring): $\\nu_{\\mu}$ NC all" 
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumunc_0
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumunc_p
       << endcolumn
       << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknumunc_p << pctsymbol 
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumunc_m 
       << endcolumn
       << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknumunc_m << pctsymbol 
       << endrow
       << endl;
  cout << setw(ncol1)  << setfill(' ') << "SK MC (1 $\\mu$-like ring): $\\bar{\\nu_{\\mu}}$ CCQE" 
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumubarccqe_0
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumubarccqe_p
       << endcolumn
       << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknumubarccqe_p << pctsymbol 
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumubarccqe_m 
       << endcolumn
       << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknumubarccqe_m << pctsymbol 
       << endrow
       << endl;
  if(show_sk_ccqe_spectral_info) {
      cout << setw(ncol1)  << setfill(' ') << "SK MC (1 $\\mu$-like ring): $\\bar{\\nu_{\\mu}}$ CCQE ($E_{r}$ < 0.4 GeV)" 
           << endcolumn
           << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumubarccqe_0a
           << endcolumn
           << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumubarccqe_pa
           << endcolumn
           << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknumubarccqe_pa << pctsymbol 
           << endcolumn
           << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumubarccqe_ma 
           << endcolumn
           << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknumubarccqe_ma << pctsymbol 
           << endrow
           << endl;
      cout << setw(ncol1)  << setfill(' ') << "SK MC (1 $\\mu$-like ring): $\\bar{\\nu_{\\mu}}$ CCQE (0.4 <= $E_{r}$ <= 1.1 GeV)" 
           << endcolumn
           << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumubarccqe_0b
           << endcolumn
           << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumubarccqe_pb
           << endcolumn
           << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknumubarccqe_pb << pctsymbol 
           << endcolumn
           << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumubarccqe_mb 
           << endcolumn
           << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknumubarccqe_mb << pctsymbol 
           << endrow
           << endl;
      cout << setw(ncol1)  << setfill(' ') << "SK MC (1 $\\mu$-like ring): $\\bar{\\nu_{\\mu}}$ CCQE ($E_{r}$ > 1.1 GeV)" 
           << endcolumn
           << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumubarccqe_0c
           << endcolumn
           << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumubarccqe_pc
           << endcolumn
           << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknumubarccqe_pc << pctsymbol 
           << endcolumn
           << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumubarccqe_mc 
           << endcolumn
           << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknumubarccqe_mc << pctsymbol 
           << endrow
           << endl;
  }
  cout << setw(ncol1)  << setfill(' ') << "SK MC (1 $\\mu$-like ring): $\\bar{\\nu_{\\mu}}$ CC1$\\pi$" 
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumubarcc1pi_0
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumubarcc1pi_p
       << endcolumn
       << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknumubarcc1pi_p << pctsymbol 
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumubarcc1pi_m 
       << endcolumn
       << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknumubarcc1pi_m << pctsymbol 
       << endrow
       << endl;
  cout << setw(ncol1)  << setfill(' ') << "SK MC (1 $\\mu$-like ring): $\\bar{\\nu_{\\mu}}$ CC other" 
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumubarccoth_0
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumubarccoth_p
       << endcolumn
       << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknumubarccoth_p << pctsymbol 
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumubarccoth_m 
       << endcolumn
       << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknumubarccoth_m << pctsymbol 
       << endrow
       << endl;
  cout << setw(ncol1)  << setfill(' ') << "SK MC (1 $\\mu$-like ring): $\\bar{\\nu_{\\mu}}$ CC all" 
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumubarcc_0
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumubarcc_p
       << endcolumn
       << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknumubarcc_p << pctsymbol 
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumubarcc_m 
       << endcolumn
       << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknumubarcc_m << pctsymbol 
       << endrow
       << endl;
  cout << setw(ncol1)  << setfill(' ') << "SK MC (1 $\\mu$-like ring): $\\bar{\\nu_{\\mu}}$ NC all"
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumubarnc_0
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumubarnc_p
       << endcolumn
       << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknumubarnc_p << pctsymbol 
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknumubarnc_m 
       << endcolumn
       << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknumubarnc_m << pctsymbol 
       << endrow
       << endl;
  cout << setw(ncol1)  << setfill(' ') << "SK MC (1 $\\mu$-like ring): $\\nu_{e}$ CCQE" 
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknueccqe_0
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknueccqe_p
       << endcolumn
       << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknueccqe_p << pctsymbol 
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknueccqe_m 
       << endcolumn
       << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknueccqe_m << pctsymbol 
       << endrow
       << endl;
  cout << setw(ncol1)  << setfill(' ') << "SK MC (1 $\\mu$-like ring): $\\nu_{e}$ CC1$\\pi$" 
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknuecc1pi_0
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknuecc1pi_p
       << endcolumn
       << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknuecc1pi_p << pctsymbol 
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknuecc1pi_m 
       << endcolumn
       << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknuecc1pi_m << pctsymbol 
       << endrow
       << endl;
  cout << setw(ncol1)  << setfill(' ') << "SK MC (1 $\\mu$-like ring): $\\nu_{e}$ CC other" 
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknueccoth_0
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknueccoth_p
       << endcolumn
       << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknueccoth_p << pctsymbol 
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknueccoth_m 
       << endcolumn
       << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknueccoth_m << pctsymbol 
       << endrow
       << endl;
  cout << setw(ncol1)  << setfill(' ') << "SK MC (1 $\\mu$-like ring): $\\nu_{e}$ CC all" 
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknuecc_0
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknuecc_p
       << endcolumn
       << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknuecc_p << pctsymbol 
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknuecc_m 
       << endcolumn
       << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknuecc_m << pctsymbol 
       << endrow
       << endl;
  cout << setw(ncol1)  << setfill(' ') << "SK MC (1 $\\mu$-like ring): $\\nu_{e}$ NC all" 
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknuenc_0
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknuenc_p
       << endcolumn
       << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknuenc_p << pctsymbol 
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknuenc_m 
       << endcolumn
       << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknuenc_m << pctsymbol 
       << endrow
       << endl;
  cout << setw(ncol1)  << setfill(' ') << "SK MC (1 $\\mu$-like ring): $\\nu_{e}$ (signal) CCQE" 
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsksignueccqe_0
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsksignueccqe_p
       << endcolumn
       << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsksignueccqe_p << pctsymbol 
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsksignueccqe_m 
       << endcolumn
       << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsksignueccqe_m << pctsymbol 
       << endrow
       << endl;
  cout << setw(ncol1)  << setfill(' ') << "SK MC (1 $\\mu$-like ring): $\\nu_{e}$ (signal) CC1$\\pi$" 
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsksignuecc1pi_0
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsksignuecc1pi_p
       << endcolumn
       << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsksignuecc1pi_p << pctsymbol 
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsksignuecc1pi_m 
       << endcolumn
       << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsksignuecc1pi_m << pctsymbol 
       << endrow
       << endl;
  cout << setw(ncol1)  << setfill(' ') << "SK MC (1 $\\mu$-like ring): $\\nu_{e}$ (signal) CC other" 
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsksignueccoth_0
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsksignueccoth_p
       << endcolumn
       << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsksignueccoth_p << pctsymbol 
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsksignueccoth_m 
       << endcolumn
       << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsksignueccoth_m << pctsymbol 
       << endrow
       << endl;
  cout << setw(ncol1)  << setfill(' ') << "SK MC (1 $\\mu$-like ring): $\\nu_{e}$ (signal) CC all" 
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsksignuecc_0
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsksignuecc_p
       << endcolumn
       << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsksignuecc_p << pctsymbol 
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsksignuecc_m 
       << endcolumn
       << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsksignuecc_m << pctsymbol 
       << endrow
       << endl;

  if(use_tex) {
     cout << "\\hline" << endl;
  } else {
     cout << setw(150)  << setfill('-') << " " << endl;
  }

  cout << setw(ncol1)  << setfill(' ') << "ND280 CC-like Observed/MC ratio" 
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nd280corr_0
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nd280corr_p
       << endcolumn
       << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nd280corr_p << pctsymbol 
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nd280corr_m 
       << endcolumn
       << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nd280corr_m << pctsymbol 
       << endrow
       << endl;
  cout << setw(ncol1)  << setfill(' ') << "ND280 Observed CC-like" 
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nd280obs_0
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nd280obs_p
       << endcolumn
       << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nd280obs_p << pctsymbol 
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nd280obs_m 
       << endcolumn
       << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nd280obs_m << pctsymbol 
       << endrow
       << endl;
  cout << setw(ncol1)  << setfill(' ') << "ND280 MC: Total CC-like" 
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nd280mccclike_0
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nd280mccclike_p
       << endcolumn
       << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nd280mccclike_p << pctsymbol 
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nd280mccclike_m 
       << endcolumn
       << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nd280mccclike_m << pctsymbol 
       << endrow
       << endl;
  cout << setw(ncol1)  << setfill(' ') << "ND280 MC: Total CC" 
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nd280cc_0
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nd280cc_p
       << endcolumn
       << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nd280cc_p << pctsymbol 
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nd280cc_m 
       << endcolumn
       << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nd280cc_m << pctsymbol 
       << endrow
       << endl;
  cout << setw(ncol1)  << setfill(' ') << "ND280 MC: Total NC background" 
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nd280ncbkg_0
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nd280ncbkg_p
       << endcolumn
       << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nd280ncbkg_p << pctsymbol 
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nd280ncbkg_m 
       << endcolumn
       << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nd280ncbkg_m << pctsymbol 
       << endrow
       << endl;

  if(use_tex) {
     cout << "\\hline" << endl;
  } else {
     cout << setw(150)  << setfill('-') << " " << endl;
  }

  cout << setw(ncol1)  << setfill(' ') << "ND280 MC: $\\nu_{\\mu}$ CCQE" 
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nd280numuccqe_0
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nd280numuccqe_p
       << endcolumn
       << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nd280numuccqe_p << pctsymbol 
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nd280numuccqe_m 
       << endcolumn
       << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nd280numuccqe_m << pctsymbol 
       << endrow
       << endl;
  cout << setw(ncol1)  << setfill(' ') << "ND280 MC: $\\nu_{\\mu}$ CC1$\\pi$" 
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nd280numucc1pi_0
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nd280numucc1pi_p
       << endcolumn
       << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nd280numucc1pi_p << pctsymbol 
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nd280numucc1pi_m 
       << endcolumn
       << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nd280numucc1pi_m << pctsymbol 
       << endrow
       << endl;
  cout << setw(ncol1)  << setfill(' ') << "ND280 MC: $\\nu_{\\mu}$ CC other" 
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nd280numuccoth_0
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nd280numuccoth_p
       << endcolumn
       << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nd280numuccoth_p << pctsymbol 
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nd280numuccoth_m 
       << endcolumn
       << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nd280numuccoth_m << pctsymbol 
       << endrow
       << endl;
  cout << setw(ncol1)  << setfill(' ') << "ND280 MC: $\\nu_{\\mu}$ CC all" 
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nd280numucc_0
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nd280numucc_p
       << endcolumn
       << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nd280numucc_p << pctsymbol 
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nd280numucc_m 
       << endcolumn
       << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nd280numucc_m << pctsymbol 
       << endrow
       << endl;
  cout << setw(ncol1)  << setfill(' ') << "ND280 MC: $\\nu_{\\mu}$ NC all" 
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nd280numunc_0
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nd280numunc_p
       << endcolumn
       << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nd280numunc_p << pctsymbol 
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nd280numunc_m 
       << endcolumn
       << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nd280numunc_m << pctsymbol 
       << endrow
       << endl;
  cout << setw(ncol1)  << setfill(' ') << "ND280 MC: $\\bar{\\nu_{\\mu}}$ CCQE" 
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nd280numubarccqe_0
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nd280numubarccqe_p
       << endcolumn
       << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nd280numubarccqe_p << pctsymbol 
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nd280numubarccqe_m 
       << endcolumn
       << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nd280numubarccqe_m << pctsymbol 
       << endrow
       << endl;
  cout << setw(ncol1)  << setfill(' ') << "ND280 MC: $\\bar{\\nu_{\\mu}}$ CC1$\\pi$" 
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nd280numubarcc1pi_0
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nd280numubarcc1pi_p
       << endcolumn
       << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nd280numubarcc1pi_p << pctsymbol 
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nd280numubarcc1pi_m 
       << endcolumn
       << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nd280numubarcc1pi_m << pctsymbol 
       << endrow
       << endl;
  cout << setw(ncol1)  << setfill(' ') << "ND280 MC: $\\bar{\\nu_{\\mu}}$ CC other" 
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nd280numubarccoth_0
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nd280numubarccoth_p
       << endcolumn
       << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nd280numubarccoth_p << pctsymbol 
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nd280numubarccoth_m 
       << endcolumn
       << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nd280numubarccoth_m << pctsymbol 
       << endrow
       << endl;
  cout << setw(ncol1)  << setfill(' ') << "ND280 MC: $\\bar{\\nu_{\\mu}}$ CC all" 
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nd280numubarcc_0
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nd280numubarcc_p
       << endcolumn
       << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nd280numubarcc_p << pctsymbol 
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nd280numubarcc_m 
       << endcolumn
       << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nd280numubarcc_m << pctsymbol 
       << endrow
       << endl;
  cout << setw(ncol1)  << setfill(' ') << "ND280 MC: $\\bar{\\nu_{\\mu}}$ NC all" 
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nd280numubarnc_0
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nd280numubarnc_p
       << endcolumn
       << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nd280numubarnc_p << pctsymbol 
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nd280numubarnc_m 
       << endcolumn
       << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nd280numubarnc_m << pctsymbol 
       << endrow
       << endl;
  cout << setw(ncol1)  << setfill(' ') << "ND280 MC: $\\nu_{e}$ CCQE" 
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nd280nueccqe_0
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nd280nueccqe_p
       << endcolumn
       << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nd280nueccqe_p << pctsymbol 
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nd280nueccqe_m 
       << endcolumn
       << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nd280nueccqe_m << pctsymbol 
       << endrow
       << endl;
  cout << setw(ncol1)  << setfill(' ') << "ND280 MC: $\\nu_{e}$ CC1$\\pi$" 
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nd280nuecc1pi_0
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nd280nuecc1pi_p
       << endcolumn
       << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nd280nuecc1pi_p << pctsymbol 
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nd280nuecc1pi_m 
       << endcolumn
       << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nd280nuecc1pi_m << pctsymbol 
       << endrow
       << endl;
  cout << setw(ncol1)  << setfill(' ') << "ND280 MC: $\\nu_{e}$ CC other" 
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nd280nueccoth_0
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nd280nueccoth_p
       << endcolumn
       << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nd280nueccoth_p << pctsymbol 
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nd280nueccoth_m 
       << endcolumn
       << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nd280nueccoth_m << pctsymbol 
       << endrow
       << endl;
  cout << setw(ncol1)  << setfill(' ') << "ND280 MC: $\\nu_{e}$ CC all" 
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nd280nuecc_0
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nd280nuecc_p
       << endcolumn
       << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nd280nuecc_p << pctsymbol 
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nd280nuecc_m 
       << endcolumn
       << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nd280nuecc_m << pctsymbol 
       << endrow
       << endl;
  cout << setw(ncol1)  << setfill(' ') << "ND280 MC: $\\nu_{e}$ NC all" 
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nd280nuenc_0
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nd280nuenc_p
       << endcolumn
       << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nd280nuenc_p << pctsymbol 
       << endcolumn
       << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nd280nuenc_m 
       << endcolumn
       << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nd280nuenc_m << pctsymbol 
       << endrow
       << endl;

  if(use_tex) {
     cout << "\\hline" << endl;
     cout << "\\end{tabular}" << endl;
     cout << "\\caption{ Effect of " << systematic << " }" << endl;
     cout << "\\label{tab:num_" << systematic << "}" << endl;
     cout << "\\end{table}" << endl;
  }


  cout << "Done!" << endl;
}
