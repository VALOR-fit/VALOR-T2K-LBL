//
// ------------------------------------------------------------------------------------------------
// Test systematic errors.
//
// The input (systematic) string can be:
//	- sk_energy_scale 
//	- sk_1Rmu_ccqe_other_eff
//	- sk_1Rmu_ccqe_ring_counting_eff : 
//             This energy-dependent correlated error can not be used in this test macro. Instead,
//             you can tweak specific components of the efficiency tweaking vector (neglecting,
//             correlation which, in the fitter, are handled using the Cholesky method). Use:
//	       - sk_1Rmu_ccqe_ring_counting_eff_0 
//	       - sk_1Rmu_ccqe_ring_counting_eff_1
//	       - sk_1Rmu_ccqe_ring_counting_eff_2
//	- sk_1Rmu_ccnonqe_eff
//	- sk_1Rmu_nc_eff
//	- sk_1Rmu_nuecc_eff
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//
// Major revisions by TD. 25/11/11.
//  Now takes SK sample as input
//  Can write output tables to tex file
//  Restructuring of SK modes part of code means its easier to change in the future:
//   - if order or number of kSKReactionModes[s][m] changes, just need to modify lines ~120-160
//

#include <iomanip.h>

void test_systematic_error_propagation
(
  const char * systematic,                                                  /* options: see above */
  const char * sksample                   = "1Rmu",                         /* options: "1Rmu", "1Re" */
  const char * hypothesis                 = "3f_osc",                       /* options: "2f_numu_disap", "3f_osc" */
  double       input_mixing12             = 0.86,                           /* */
  double       input_mixing13             = 0.10,                           /* */
  double       input_mixing23             = 1.00,                           /* */
  double       input_deltacp              = 0.00,                           /* */
  double       input_dm12sq               = 7.58E-5,                        /* */
  double       input_dm23sq               = 2.32E-3,                        /* */
  bool         use_tex                    = false,                          /* print-out results as a TeX table (for inclusion in tech. note)*/
  bool         show_sk_ccqe_spectral_info = false,                          /* */
  const char * output_file                = "",                             /* */
  const char * config_file                = "analysis_configuration.xml"    /* */
)
{
  //need to declare static variables for array size at the start of function due to CINT bug
  //this is the max number of modes that make up a certain category in the table (total nc/cc for each flavour)
  static const int nmodesallmax = 4;

  using namespace std;
  using namespace t2kana;

  // Analysis configuration
  const char * config  = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);
    
  // Global options
  Opt::Instance()->gOscMixingParamChoice = conv::kNuMixSinsq2Theta;
  Opt::Instance()->gMassHierarchyAssumed = conv::kNmhNormal;

  // Get requested sample
  conv::SuperkSample_t s = util::SKSampleFromString(sksample);
  assert(s != conv::kSKsUndefined);

  // Set physics hypothesis
  PhysHypothesis ph(hypothesis);
  if(ph.Hypothesis() == conv::kHypNuMuDisapp2F) {
    ph.SetCurrValue(conv::kPhypMixing23,  input_mixing23);
    ph.SetCurrValue(conv::kPhypDm23sq,    input_dm23sq);
  } else
  if(ph.Hypothesis() == conv::kHypOsc3F) {
    ph.SetCurrValue(conv::kPhypMixing12,  input_mixing12);
    ph.SetCurrValue(conv::kPhypMixing13,  input_mixing13);
    ph.SetCurrValue(conv::kPhypMixing23,  input_mixing23);
    ph.SetCurrValue(conv::kPhypDeltaCP,   input_deltacp);
    ph.SetCurrValue(conv::kPhypDm12sq,    input_dm12sq);
    ph.SetCurrValue(conv::kPhypDm23sq,    input_dm23sq);
  }

  // Print options
  cout << "\n[Job inputs]" << endl;
  cout << "- Reading analysis configuration from: " << config << endl;
  cout << "- SuperK sample: " << util::AsString(s) << endl;
  cout << "- Varied systematic parameter: " << systematic << endl;
  cout << "- Physics hypothesis: " << ph << endl;
  cout << *Opt::Instance() << endl;

  // Precison for printed-out numerical results
  int nprecisionnev = 4; // <-- for number of events
  int nprecisionpct = 2; // <-- for % changes

  // Create & config superk spectrum predictor
  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read(config);

  //
  int ErecoBin1 =  9; // 0.4 GeV (incl.)
  int ErecoBin2 = 22; // 1.1 GeV (incl.)

  // Define arrays with entries the modes that make up the total nc/cc-flavour category.
  // This makes the code easier to adapt if the list of modes considered for each sample changes.
  int insknumuccall   [nmodesallmax];
  int insknumuncall   [nmodesallmax];
  int insknumubarccall[nmodesallmax];
  int insknumubarncall[nmodesallmax];
  int insknueccall    [nmodesallmax];
  int insknuencall    [nmodesallmax];
  int insknuesigall   [nmodesallmax];
  // initialize arrays
  for (int i = 0; i < nmodesallmax; i++) {
       insknumuccall   [i] = -1;
       insknumuncall   [i] = -1;
       insknumubarccall[i] = -1;
       insknumubarncall[i] = -1;
       insknueccall    [i] = -1;
       insknuencall    [i] = -1;
       insknuesigall   [i] = -1;
  }
  // fill arrays
  switch(s) {
     case (conv::kSKsFCFV1RingMuLike): {
       insknumuccall   [0] = 0;
       insknumuccall   [1] = 1;
       insknumuccall   [2] = 2;
       insknumubarccall[0] = 4;
       insknumubarccall[1] = 5;
       insknumubarccall[2] = 6;
       insknueccall    [0] = 8;
       insknueccall    [1] = 9;
       insknueccall    [2] = 10;
       insknuesigall   [0] = 12;
       insknuesigall   [1] = 13;
       insknuesigall   [2] = 14;
       break;
     }
     case (conv::kSKsFCFV1RingELike): {
       insknumuccall   [0] = 0;
       insknumuccall   [1] = 1;
       insknumuccall   [2] = 2;
       insknumuccall   [3] = 3;
       insknumuncall   [0] = 4;
       insknumuncall   [1] = 5;
       insknumuncall   [2] = 6;
       insknumubarccall[0] = 7;
       insknumubarccall[1] = 8;
       insknumubarccall[2] = 9;
       insknumubarccall[3] = 10;
       insknumubarncall[0] = 11;
       insknumubarncall[1] = 12;
       insknumubarncall[2] = 13;
       insknueccall    [0] = 14;
       insknueccall    [1] = 15;
       insknueccall    [2] = 16;
       insknueccall    [3] = 17;
       insknuencall    [0] = 18;
       insknuencall    [1] = 19;
       insknuencall    [2] = 20;
       insknuesigall   [0] = 21;
       insknuesigall   [1] = 22;
       insknuesigall   [2] = 23;
       break;
     }
     default: 
       cout << "*** Sample " << util::AsString(s) << " not supported yet" << endl; 
       exit(1); 
  }
  //and get the last non zero entry in the array (for use to tell you when to print the total)
  int sknumuccall_printat=-1;
  int sknumuncall_printat=-1;
  int sknumubarccall_printat=-1;
  int sknumubarncall_printat=-1;
  int sknueccall_printat=-1;
  int sknuencall_printat=-1;
  int sknuesigall_printat=-1;
  for( int i = nmodesallmax-1; i>=0; i-- ) {
    if(insknumuccall[i]>=0 && sknumuccall_printat<0)
      sknumuccall_printat = insknumuccall[i];
    if(insknumuncall[i]>=0 && sknumuncall_printat<0)
      sknumuncall_printat = insknumuncall[i];
    if(insknumubarccall[i]>=0 && sknumubarccall_printat<0)
      sknumubarccall_printat = insknumubarccall[i];
    if(insknumubarncall[i]>=0 && sknumubarncall_printat<0)
      sknumubarncall_printat = insknumubarncall[i];
    if(insknueccall[i]>=0 && sknueccall_printat<0)
      sknueccall_printat = insknueccall[i];
    if(insknuencall[i]>=0 && sknuencall_printat<0)
      sknuencall_printat = insknuencall[i];
    if(insknuesigall[i]>=0 && sknuesigall_printat<0)
      sknuesigall_printat = insknuesigall[i];
  }

  //
  // nominal
  //

  skexp->Generate(ph);

  // total 
  double nsk_0 = skexp->CurrPredictedErecoSpectrum(s)->Integral();

  // individual modes 
  double nskmodes_0[conv::knm];
  for (int im = 0; im < conv::knm; im++) {
    conv::SuperkMode_t m = conv::kSKReactionModes[s][im];
    if(m == conv::kSKmUndefined) continue;
    // use s, m (not is,im) as indices, eg fMCErecoSpectrum[s][m]
    nskmodes_0[m] = skexp->CurrPredictedErecoSpectrum(s,m)->Integral();
  }

  // ccnumu
  double nsknumucc_0 = 0;
  for (int i = 0; i<nmodesallmax; i++) {
    if(insknumuccall[i] < 0) continue;
    conv::SuperkMode_t m = conv::kSKReactionModes[s][insknumuccall[i]];
    if(m == conv::kSKmUndefined) continue;
    nsknumucc_0 += nskmodes_0[m];
  }

  // ncnumu
  double nsknumunc_0 = 0;
  for (int i = 0; i<nmodesallmax; i++) {
    if(insknumuncall[i] < 0) continue;
    conv::SuperkMode_t m = conv::kSKReactionModes[s][insknumuncall[i]];
    if(m == conv::kSKmUndefined) continue;
    nsknumunc_0 += nskmodes_0[m];
  }

  // ccnumubar
  double nsknumubarcc_0 = 0;
  for (int i = 0; i<nmodesallmax; i++) {
    if(insknumubarccall[i] < 0) continue;
    conv::SuperkMode_t m = conv::kSKReactionModes[s][insknumubarccall[i]];
    if(m == conv::kSKmUndefined) continue;
    nsknumubarcc_0 += nskmodes_0[m];
  }

  // ncnumubar
  double nsknumubarnc_0 = 0;
  for (int i = 0; i<nmodesallmax; i++) {
    if(insknumubarncall[i] < 0) continue;
    conv::SuperkMode_t m = conv::kSKReactionModes[s][insknumubarncall[i]];
    if(m == conv::kSKmUndefined) continue;
    nsknumubarnc_0 += nskmodes_0[m];
  }

  // ccnue
  double nsknuecc_0 = 0;
  for (int i = 0; i<nmodesallmax; i++) {
    if(insknueccall[i] < 0) continue;
    conv::SuperkMode_t m = conv::kSKReactionModes[s][insknueccall[i]];
    if(m == conv::kSKmUndefined) continue;
    nsknuecc_0 += nskmodes_0[m];
  }

  // ncnue
  double nsknuenc_0 = 0;
  for (int i = 0; i<nmodesallmax; i++) {
    if(insknuencall[i] < 0) continue;
    conv::SuperkMode_t m = conv::kSKReactionModes[s][insknuencall[i]];
    if(m == conv::kSKmUndefined) continue;
    nsknuenc_0 += nskmodes_0[m];
  }

  // signue
  double nsknuesig_0 = 0;
  for (int i = 0; i<nmodesallmax; i++) {
    if(insknuesigall[i] < 0) continue;
    conv::SuperkMode_t m = conv::kSKReactionModes[s][insknuesigall[i]];
    if(m == conv::kSKmUndefined) continue;
    nsknuesig_0 += nskmodes_0[m];
  }

  // spectral numu ccqe
  double nsknumuccqe_0a        = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuMuCCQE       )->Integral(1,           ErecoBin1-1          );
  double nsknumuccqe_0b        = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuMuCCQE       )->Integral(ErecoBin1,   ErecoBin2            );
  double nsknumuccqe_0c        = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuMuCCQE       )->Integral(ErecoBin2+1, conv::kNErecoBins1Rmu);
  double nsknumubarccqe_0a     = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuMuBarCCQE    )->Integral(1,           ErecoBin1-1          );
  double nsknumubarccqe_0b     = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuMuBarCCQE    )->Integral(ErecoBin1,   ErecoBin2            );
  double nsknumubarccqe_0c     = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuMuBarCCQE    )->Integral(ErecoBin2+1, conv::kNErecoBins1Rmu);


  //
  // +1 sigma
  //

  skexp->TweakSystematic(systematic,+1.);
  skexp->Generate(ph);

  // total 
  double nsk_p = skexp->CurrPredictedErecoSpectrum(s)->Integral();

  // individual modes 
  double nskmodes_p[conv::knm];
  for (int im = 0; im < conv::knm; im++) {
    conv::SuperkMode_t m = conv::kSKReactionModes[s][im];
    if(m == conv::kSKmUndefined) continue;
    // use s, m (not is,im) as indices, eg fMCErecoSpectrum[s][m]
    nskmodes_p[m] = skexp->CurrPredictedErecoSpectrum(s, m)->Integral();
  }

  // ccnumu
  double nsknumucc_p = 0;
  for (int i = 0; i<nmodesallmax; i++) {
    if (insknumuccall[i] < 0) continue;
    conv::SuperkMode_t m = conv::kSKReactionModes[s][insknumuccall[i]];
    if(m == conv::kSKmUndefined) continue;
    nsknumucc_p += nskmodes_p[m];
  }

  // ncnumu
  double nsknumunc_p = 0;
  for (int i = 0; i<nmodesallmax; i++) {
    if (insknumuncall[i] < 0) continue;
    conv::SuperkMode_t m = conv::kSKReactionModes[s][insknumuncall[i]];
    if(m == conv::kSKmUndefined) continue;
    nsknumunc_p += nskmodes_p[m];
  }

  // ccnumubar
  double nsknumubarcc_p = 0;
  for (int i = 0; i<nmodesallmax; i++) {
    if (insknumubarccall[i] < 0) continue;
    conv::SuperkMode_t m = conv::kSKReactionModes[s][insknumubarccall[i]];
    if(m == conv::kSKmUndefined) continue;
    nsknumubarcc_p += nskmodes_p[m];
  }

  // ncnumubar
  double nsknumubarnc_p = 0;
  for (int i = 0; i<nmodesallmax; i++) {
    if (insknumubarncall[i] < 0) continue;
    conv::SuperkMode_t m = conv::kSKReactionModes[s][insknumubarncall[i]];
    if(m == conv::kSKmUndefined) continue;
    nsknumubarnc_p += nskmodes_p[m];
  }

  // ccnue
  double nsknuecc_p = 0;
  for (int i = 0; i<nmodesallmax; i++) {
    if (insknueccall[i] < 0) continue;
    conv::SuperkMode_t m = conv::kSKReactionModes[s][insknueccall[i]];
    if(m == conv::kSKmUndefined) continue;
    nsknuecc_p += nskmodes_p[m];
  }

  // ncnue
  double nsknuenc_p = 0;
  for (int i = 0; i<nmodesallmax; i++) {
    if (insknuencall[i] < 0) continue;
    conv::SuperkMode_t m = conv::kSKReactionModes[s][insknuencall[i]];
    if(m == conv::kSKmUndefined) continue;
    nsknuenc_p += nskmodes_p[m];
  }

  // signue
  double nsknuesig_p = 0;
  for (int i = 0; i<nmodesallmax; i++) {
    if (insknuesigall[i] < 0) continue;
    conv::SuperkMode_t m = conv::kSKReactionModes[s][insknuesigall[i]];
    if(m == conv::kSKmUndefined) continue;
    nsknuesig_p += nskmodes_p[m];
  }

  // spectral numu ccqe
  double nsknumuccqe_pa        = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuMuCCQE       )->Integral(1,           ErecoBin1-1          );
  double nsknumuccqe_pb        = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuMuCCQE       )->Integral(ErecoBin1,   ErecoBin2            );
  double nsknumuccqe_pc        = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuMuCCQE       )->Integral(ErecoBin2+1, conv::kNErecoBins1Rmu);
  double nsknumubarccqe_pa     = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuMuBarCCQE    )->Integral(1,           ErecoBin1-1          );
  double nsknumubarccqe_pb     = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuMuBarCCQE    )->Integral(ErecoBin1,   ErecoBin2            );
  double nsknumubarccqe_pc     = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuMuBarCCQE    )->Integral(ErecoBin2+1, conv::kNErecoBins1Rmu);

  //
  // -1 sigma
  //

  skexp->TweakSystematic(systematic,-1.);
  skexp->Generate(ph);

  // total 
  double nsk_m = skexp->CurrPredictedErecoSpectrum(s)->Integral();

  // individual modes 
  double nskmodes_m[conv::knm];
  for (int im = 0; im < conv::knm; im++) {
    conv::SuperkMode_t m = conv::kSKReactionModes[s][im];
    if(m == conv::kSKmUndefined) continue;
    //   use s, m (not is,im) as indices, eg fMCErecoSpectrum[s][m]
    nskmodes_m[m] = skexp->CurrPredictedErecoSpectrum(s, m)->Integral();
  }

  // ccnumu
  double nsknumucc_m = 0;
  for (int i = 0; i<nmodesallmax; i++) {
    if (insknumuccall[i] < 0) continue;
    conv::SuperkMode_t m = conv::kSKReactionModes[s][insknumuccall[i]];
    if(m == conv::kSKmUndefined) continue;
    nsknumucc_m += nskmodes_m[m];
  }

  // ncnumu
  double nsknumunc_m = 0;
  for (int i = 0; i<nmodesallmax; i++) {
    if (insknumuncall[i] < 0) continue;
    conv::SuperkMode_t m = conv::kSKReactionModes[s][insknumuncall[i]];
    if(m == conv::kSKmUndefined) continue;
    nsknumunc_m += nskmodes_m[m];
  }

  // ccnumubar
  double nsknumubarcc_m = 0;
  for (int i = 0; i<nmodesallmax; i++) {
    if (insknumubarccall[i] < 0) continue;
    conv::SuperkMode_t m = conv::kSKReactionModes[s][insknumubarccall[i]];
    if(m == conv::kSKmUndefined) continue;
    nsknumubarcc_m += nskmodes_m[m];
  }

  // ncnumubar
  double nsknumubarnc_m = 0;
  for (int i = 0; i<nmodesallmax; i++) {
    if (insknumubarncall[i] < 0) continue;
    conv::SuperkMode_t m = conv::kSKReactionModes[s][insknumubarncall[i]];
    if(m == conv::kSKmUndefined) continue;
    nsknumubarnc_m += nskmodes_m[m];
  }

  // ccnue
  double nsknuecc_m = 0;
  for (int i = 0; i<nmodesallmax; i++) {
    if (insknueccall[i] < 0) continue;
    conv::SuperkMode_t m = conv::kSKReactionModes[s][insknueccall[i]];
    if(m == conv::kSKmUndefined) continue;
    nsknuecc_m += nskmodes_m[m];
  }

  // ncnue
  double nsknuenc_m = 0;
  for (int i = 0; i<nmodesallmax; i++) {
    if (insknuencall[i] < 0) continue;
    conv::SuperkMode_t m = conv::kSKReactionModes[s][insknuencall[i]];
    if(m == conv::kSKmUndefined) continue;
    nsknuenc_m += nskmodes_m[m];
  }

  // signue
  double nsknuesig_m = 0;
  for (int i = 0; i<nmodesallmax; i++) {
    if (insknuesigall[i] < 0) continue;
    conv::SuperkMode_t m = conv::kSKReactionModes[s][insknuesigall[i]];
    if(m == conv::kSKmUndefined) continue;
    nsknuesig_m += nskmodes_m[m];
  }

  // spectral numu ccqe
  double nsknumuccqe_ma        = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuMuCCQE       )->Integral(1,           ErecoBin1-1          );
  double nsknumuccqe_mb        = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuMuCCQE       )->Integral(ErecoBin1,   ErecoBin2            );
  double nsknumuccqe_mc        = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuMuCCQE       )->Integral(ErecoBin2+1, conv::kNErecoBins1Rmu);
  double nsknumubarccqe_ma     = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuMuBarCCQE    )->Integral(1,           ErecoBin1-1          );
  double nsknumubarccqe_mb     = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuMuBarCCQE    )->Integral(ErecoBin1,   ErecoBin2            );
  double nsknumubarccqe_mc     = skexp->CurrPredictedErecoSpectrum(s, conv::kSKmNuMuBarCCQE    )->Integral(ErecoBin2+1, conv::kNErecoBins1Rmu);


  //
  // calc fractional changes
  //

  // totals
  double frac_nsk_p                  = (nsk_0                  > 0) ? 100.*(nsk_p                 - nsk_0                ) / nsk_0                 : 0.;
  double frac_nsknumucc_p            = (nsknumucc_0            > 0) ? 100.*(nsknumucc_p           - nsknumucc_0          ) / nsknumucc_0           : 0.;
  double frac_nsknumunc_p            = (nsknumunc_0            > 0) ? 100.*(nsknumunc_p           - nsknumunc_0          ) / nsknumunc_0           : 0.;
  double frac_nsknumubarcc_p         = (nsknumubarcc_0         > 0) ? 100.*(nsknumubarcc_p        - nsknumubarcc_0       ) / nsknumubarcc_0        : 0.;
  double frac_nsknumubarnc_p         = (nsknumubarnc_0         > 0) ? 100.*(nsknumubarnc_p        - nsknumubarnc_0       ) / nsknumubarnc_0        : 0.;
  double frac_nsknuecc_p             = (nsknuecc_0             > 0) ? 100.*(nsknuecc_p            - nsknuecc_0           ) / nsknuecc_0            : 0.;
  double frac_nsknuenc_p             = (nsknuenc_0             > 0) ? 100.*(nsknuenc_p            - nsknuenc_0           ) / nsknuenc_0            : 0.;
  double frac_nsknuesig_p            = (nsknuesig_0            > 0) ? 100.*(nsknuesig_p           - nsknuesig_0          ) / nsknuesig_0           : 0.;

  // individual modes
  double frac_nskmodes_p[conv::knm];
  for (int im = 0; im < conv::knm; im++) {
    conv::SuperkMode_t m = conv::kSKReactionModes[s][im];
    if(m == conv::kSKmUndefined) continue;
    //   use s, m (not is,im) as indices, eg fMCErecoSpectrum[s][m]
    frac_nskmodes_p[m] = (nskmodes_0[m] > 0) ?  100.*(nskmodes_p[m] - nskmodes_0[m]) / nskmodes_0[m]  :  0.;
  }
  // spectral ccqe
  double frac_nsknumuccqe_pa         = (nsknumuccqe_0a         > 0) ? 100.*(nsknumuccqe_pa        - nsknumuccqe_0a       ) / nsknumuccqe_0a        : 0.;
  double frac_nsknumuccqe_pb         = (nsknumuccqe_0b         > 0) ? 100.*(nsknumuccqe_pb        - nsknumuccqe_0b       ) / nsknumuccqe_0b        : 0.;
  double frac_nsknumuccqe_pc         = (nsknumuccqe_0c         > 0) ? 100.*(nsknumuccqe_pc        - nsknumuccqe_0c       ) / nsknumuccqe_0c        : 0.;
  double frac_nsknumubarccqe_pa      = (nsknumubarccqe_0a      > 0) ? 100.*(nsknumubarccqe_pa     - nsknumubarccqe_0a    ) / nsknumubarccqe_0a     : 0.;
  double frac_nsknumubarccqe_pb      = (nsknumubarccqe_0b      > 0) ? 100.*(nsknumubarccqe_pb     - nsknumubarccqe_0b    ) / nsknumubarccqe_0b     : 0.;
  double frac_nsknumubarccqe_pc      = (nsknumubarccqe_0c      > 0) ? 100.*(nsknumubarccqe_pc     - nsknumubarccqe_0c    ) / nsknumubarccqe_0c     : 0.;


  // totals
  double frac_nsk_m                  = (nsk_0                  > 0) ? 100.*(nsk_m                 - nsk_0                ) / nsk_0                 : 0.;
  double frac_nsknumucc_m            = (nsknumucc_0            > 0) ? 100.*(nsknumucc_m           - nsknumucc_0          ) / nsknumucc_0           : 0.;
  double frac_nsknumunc_m            = (nsknumunc_0            > 0) ? 100.*(nsknumunc_m           - nsknumunc_0          ) / nsknumunc_0           : 0.;
  double frac_nsknumubarcc_m         = (nsknumubarcc_0         > 0) ? 100.*(nsknumubarcc_m        - nsknumubarcc_0       ) / nsknumubarcc_0        : 0.;
  double frac_nsknumubarnc_m         = (nsknumubarnc_0         > 0) ? 100.*(nsknumubarnc_m        - nsknumubarnc_0       ) / nsknumubarnc_0        : 0.;
  double frac_nsknuecc_m             = (nsknuecc_0             > 0) ? 100.*(nsknuecc_m            - nsknuecc_0           ) / nsknuecc_0            : 0.;
  double frac_nsknuenc_m             = (nsknuenc_0             > 0) ? 100.*(nsknuenc_m            - nsknuenc_0           ) / nsknuenc_0            : 0.;
  double frac_nsknuesig_m            = (nsknuesig_0            > 0) ? 100.*(nsknuesig_m           - nsknuesig_0          ) / nsknuesig_0           : 0.;

  // individual modes
  double frac_nskmodes_m[conv::knm];
  for (int im = 0; im < conv::knm; im++) {
    conv::SuperkMode_t m = conv::kSKReactionModes[s][im];
    if(m == conv::kSKmUndefined) continue;
    //   use s, m (not is,im) as indices, eg fMCErecoSpectrum[s][m]
    frac_nskmodes_m[m] = (nskmodes_0[m] > 0) ?  100.*(nskmodes_m[m] - nskmodes_0[m]) / nskmodes_0[m]  :  0.;
  }

  // spectral ccqe
  double frac_nsknumuccqe_ma         = (nsknumuccqe_0a         > 0) ? 100.*(nsknumuccqe_ma        - nsknumuccqe_0a       ) / nsknumuccqe_0a        : 0.;
  double frac_nsknumuccqe_mb         = (nsknumuccqe_0b         > 0) ? 100.*(nsknumuccqe_mb        - nsknumuccqe_0b       ) / nsknumuccqe_0b        : 0.;
  double frac_nsknumuccqe_mc         = (nsknumuccqe_0c         > 0) ? 100.*(nsknumuccqe_mc        - nsknumuccqe_0c       ) / nsknumuccqe_0c        : 0.;
  double frac_nsknumubarccqe_ma      = (nsknumubarccqe_0a      > 0) ? 100.*(nsknumubarccqe_ma     - nsknumubarccqe_0a    ) / nsknumubarccqe_0a     : 0.;
  double frac_nsknumubarccqe_mb      = (nsknumubarccqe_0b      > 0) ? 100.*(nsknumubarccqe_mb     - nsknumubarccqe_0b    ) / nsknumubarccqe_0b     : 0.;
  double frac_nsknumubarccqe_mc      = (nsknumubarccqe_0c      > 0) ? 100.*(nsknumubarccqe_mc     - nsknumubarccqe_0c    ) / nsknumubarccqe_0c     : 0.;

  //
  // Print-out
  //

  string pctsymbol  = (use_tex) ? " \\% "   : "% ";
  string endcolumn  = (use_tex) ? " & "     : "" ;
  string endrow     = (use_tex) ? " \\\\ "  : "" ;

  int ncol1   = 80;
  int ncolnev = 15;
  int ncolfrc =  8;

  cout << "\n\nSystematic: " << systematic << endl;
  cout << endl;
  cout << fixed;

  stringstream streamer;

  if(use_tex) {
     streamer << "\\begin{table}"   << endl;
     streamer << "\\footnotesize" << endl;
     streamer << "\\centering"     << endl;
     streamer << "\\begin{tabular}{| l | r | r  r | r  r |}" << endl;
     streamer << "\\hline" << endl;
     streamer << "{\\bf Component} &  "
          << "{\\bf Nominal}   &  "
          << "\\multicolumn {2}{c|} { {\\bf +1$\\sigma$} } &  "
          << "\\multicolumn {2}{c|} { {\\bf -1$\\sigma$} } \\\\"
          << endl;
     streamer << "\\hline" << endl;
  }
 
  string temp = "$N_{SK}$ (" + util::AsStringTex(s) + ")";
  streamer << setw(ncol1)  << setfill(' ') << temp
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
     streamer << "\\hline" << endl;
  } else {
     streamer << setw(150)  << setfill('-') << " " << endl;
  }

  for (int im = 0; im < conv::knm; im++) {
    conv::SuperkMode_t m = conv::kSKReactionModes[s][im];
    if(m == conv::kSKmUndefined) continue;
    //a normal individual mode table fill
    string component = "$N_{SK}$ (";
    component += util::AsStringTex(s);
    component += "): ";
    component +=  util::AsStringTex(m);
    streamer << setw(ncol1)  << setfill(' ') << component
	 << endcolumn
	 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nskmodes_0[m]
	 << endcolumn
	 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nskmodes_p[m]
	 << endcolumn
	 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nskmodes_p[m] << pctsymbol 
	 << endcolumn
	 << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nskmodes_m[m] 
	 << endcolumn
	 << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nskmodes_m[m] << pctsymbol 
	 << endrow
	 << endl;
    //spectral ccqe
    if ( m == conv::kSKmNuMuCCQE && show_sk_ccqe_spectral_info) {
      temp = component + " ($E_{r} <$ 0.4 GeV)";
      streamer << setw(ncol1)  << setfill(' ') << temp
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
      temp = component + " (0.4 $\leq E_{r} \geq$ 1.1 GeV)";
      streamer << setw(ncol1)  << setfill(' ') << temp
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
      temp = component + " ($E_{r} >$ 1.1 GeV)";
      streamer << setw(ncol1)  << setfill(' ') << temp
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
    }//numu ccqe spectral
    if ( m == conv::kSKmNuMuBarCCQE && show_sk_ccqe_spectral_info) {
      temp = component + " ($E_{r} <$ 0.4 GeV)";
      streamer << setw(ncol1)  << setfill(' ') << temp
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
      temp = component + " (0.4 $\leq E_{r} \geq$ 1.1 GeV)";
      streamer << setw(ncol1)  << setfill(' ') << temp
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
      temp = component + " ($E_{r} >$ 1.1 GeV)";
      streamer << setw(ncol1)  << setfill(' ') << temp
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
    }//numbar ccqe spectral
    //total flavour-nc/cc
    if (sknumuccall_printat > 0 && m == conv::kSKReactionModes[s][sknumuccall_printat]) {
      temp = "$N_{SK}$ (" + util::AsStringTex(s) + "): $\\nu_{\\mu}$ CC all";
      streamer << setw(ncol1)  << setfill(' ') << temp
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
    }
    if (sknumuncall_printat > 0 && m == conv::kSKReactionModes[s][sknumuncall_printat]) {
      temp = "$N_{SK}$ (" + util::AsStringTex(s) + "): $\\nu_{\\mu}$ NC all";
      streamer << setw(ncol1)  << setfill(' ') << temp
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
    }
    if (sknumubarccall_printat > 0 && m == conv::kSKReactionModes[s][sknumubarccall_printat]) {
      temp = "$N_{SK}$ (" + util::AsStringTex(s) + "): $\\bar{\\nu_{\\mu}}$ CC all";
      streamer << setw(ncol1)  << setfill(' ') << temp
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
    }
    if (sknumubarncall_printat > 0 && m == conv::kSKReactionModes[s][sknumubarncall_printat]) {
      temp = "$N_{SK}$ (" + util::AsStringTex(s) + "): $\\bar{\\nu_{\\mu}}$ NC all";
      streamer << setw(ncol1)  << setfill(' ') << temp
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
    }
    if (sknueccall_printat > 0 && m == conv::kSKReactionModes[s][sknueccall_printat]) {
      temp = "$N_{SK}$ (" + util::AsStringTex(s) + "): $\\nu_{e}$ CC all";
      streamer << setw(ncol1)  << setfill(' ') << temp
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
    }
    if (sknuencall_printat > 0 && m == conv::kSKReactionModes[s][sknuencall_printat]) {
      temp = "$N_{SK}$ (" + util::AsStringTex(s) + "): $\\nu_{e}$ NC all";
      streamer << setw(ncol1)  << setfill(' ') << temp
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
    }
    if (sknuesigall_printat > 0 && m == conv::kSKReactionModes[s][sknuesigall_printat]) {
      temp = "$N_{SK}$ (" + util::AsStringTex(s) + "): $\\nu_{e}$ (signal) CC all";
      streamer << setw(ncol1)  << setfill(' ') << temp
	   << endcolumn
	   << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknuesig_0
	   << endcolumn
	   << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknuesig_p
	   << endcolumn
	   << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknuesig_p << pctsymbol 
	   << endcolumn
	   << setw(ncolnev)  << setfill(' ') << setprecision(nprecisionnev) << nsknuesig_m 
	   << endcolumn
	   << setw(ncolfrc)  << setfill(' ') << setprecision(nprecisionpct) << frac_nsknuesig_m << pctsymbol 
	   << endrow
	   << endl;
    }
  }//end loop over sk modes

  if(use_tex) {
     streamer << "\\hline" << endl;
  } else {
     streamer << setw(150)  << setfill('-') << " " << endl;
  }


  if(use_tex) {
     streamer << "\\hline" << endl;
     streamer << "\\end{tabular}" << endl;
     string systematicstr(systematic), hypothesisstr(hypothesis);
     size_t found;
     found=systematicstr.find_first_of("_");
     while(found!=string::npos) {
       systematicstr[found]='-';
       found=systematicstr.find_first_of("_",found+1);
     }
     found=hypothesisstr.find_first_of("_");
     while(found!=string::npos) {
       hypothesisstr[found]='-';
       found=hypothesisstr.find_first_of("_",found+1);
     }
     streamer << "\\caption{ Effect of " << systematicstr << " for sample " << util::AsStringTex(s) 
	      << " with oscillation hypothesis " <<  hypothesisstr
	      << " using oscillation parameters " 
	      << " $\\sin^22\\theta_{12} = $" << input_mixing12
	      << " $\\sin^22\\theta_{13} = $" << input_mixing13
	      << " $\\sin^22\\theta_{23} = $" << input_mixing23
	      << " $\\delta_{CP} = $        " << input_deltacp
	      << " $\\Delta m^2_{12} = $    " << input_dm12sq
	      << " $\\Delta m^2_{23} = $    " << input_dm23sq
	      << " }" << endl;
     string shortsamplestr = s==conv::kSKsFCFV1RingMuLike ? "1RmuT_" : "1ReT_";
     string spectralstr = show_sk_ccqe_spectral_info ? "_spectral" : "";
     streamer << "\\label{tab:num_" << systematic << "_" << shortsamplestr << hypothesis 
	      << "_t13:" << input_mixing13 << spectralstr << "}" << endl;
     streamer << "\\end{table}" << endl;
     streamer << "\\clearpage" << endl;
  }

  cout << streamer.str() << endl;

  // write the table?
  if(output_file != "") {
    long size;

    // get size of file
    streamer.seekg(0,ifstream::end);
    size=streamer.tellg();
    streamer.seekg(0);

    ofstream outfile (output_file,ios_base::app);

    // write to outfile
    outfile.write (streamer.str().c_str(),size);
  
    outfile.close();
  }

  cout << "Done!" << endl;

}
