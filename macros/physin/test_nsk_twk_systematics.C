//
// ------------------------------------------------------------------------------------------------
// Summarize Nsk change as a result of tweaking systematic parameters.
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

void test_nsk_twk_systematics
(
  const char * sksample           = "1Rmu",
  const char * hypothesis         = "3f_osc",
  double input_mixing12           = -99,
  double input_mixing13           = -99,
  double input_mixing23           = -99,
  double input_deltacp            = -99,
  double input_dm12sq             = -99,
  double input_dm23sq             = -99,
  const char * config_file        = "analysis_configuration.xml"
)
{
  gROOT->Reset();

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
  if (s != conv::kSKsFCFV1RingMuLike) {
    cout << "Sample not supported: " << util::AsString(s) << endl;
    return;
  }

  // Set physics hypothesis
  if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsq2Theta) {
    if(input_mixing13 < -1) input_mixing13 = conv::kGlobFitSinsq2Theta13;
    if(input_mixing23 < -1) input_mixing23 = conv::kGlobFitSinsq2Theta23;
    if(input_mixing12 < -1) input_mixing12 = conv::kGlobFitSinsq2Theta12;
  }
  else {
    if(input_mixing13 < -1) input_mixing13 = conv::kGlobFitSinsqTheta13;
    if(input_mixing23 < -1) input_mixing23 = conv::kGlobFitSinsqTheta23;
    if(input_mixing12 < -1) input_mixing12 = conv::kGlobFitSinsqTheta12;
  }
  if(input_dm23sq < -1)  input_dm23sq = conv::kGlobFitDm23sq;
  if(input_dm12sq < -1)  input_dm12sq = conv::kGlobFitDm12sq;
  if(input_deltacp < -90) input_deltacp = conv::kGlobFitDeltaCP;
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
  cout << "- Physics hypothesis: " << ph << endl;
  cout << *Opt::Instance() << endl;
  
  const int nexperiments = 5000;

  const int nsystematics = 21;
  string systematic[nsystematics] = {
    "sk_1Rmu_ccqe_ring_counting_eff", //0
    "sk_1Rmu_ccqe_other_eff",    //1
    "sk_1Rmu_ccnonqe_eff",       //2
    "sk_1Rmu_nc_eff",            //3
    "sk_1Rmu_nuecc_eff",         //4
    "sk_1Rmu_energy_scale",      //5
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
    "fsi"                        //20
  };
  bool generate_random[nsystematics] = {
    true,  /* sk_1Rmu_ccqe_ring_counting_eff */
    false, /* sk_1Rmu_ccqe_other_eff    */
    false, /* sk_1Rmu_ccnonqe_eff       */
    false, /* sk_1Rmu_nc_eff            */
    false, /* sk_1Rmu_nuecc_eff         */
    false, /* sk_1Rmu_energy_scale      */
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
    true   /* fsi                       */
  };
  bool skip[nsystematics] = {
    false, /* sk_1Rmu_ccqe_ring_counting_eff */
    false, /* sk_1Rmu_ccqe_other_eff    */
    false, /* sk_1Rmu_ccnonqe_eff       */
    false, /* sk_1Rmu_nc_eff            */
    false, /* sk_1Rmu_nuecc_eff         */
    false, /* sk_1Rmu_energy_scale      */
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
    true   /* fsi                       */
  };

  double nsknnd_min[nsystematics];
  double nsknnd_max[nsystematics];
  double nsk_min[nsystematics];
  double nsk_max[nsystematics];

  // get & configure spectrum predictor
  const char * config = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);
  cout << "Reading analysis configuration from: " << config << endl;
  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read(config);
 
  // get default prediction
  skexp->Generate(ph);
  double nsknnd_deflt = skexp->CurrPredictedErecoSpectrum(s)->Integral();
  double nsk_deflt = 0;
  for (int im = 0; im < conv::knm; im++) {
    conv::SuperkMode_t m = conv::kSKReactionModes[s][im];
    if(m == conv::kSKmUndefined) continue;
    nsk_deflt += skexp->GetSuperkMC()->CurrErecoSpectrum(s,m)->Integral();
  }
  
  // get predictions with tweaked systematic params
  for(int isystematic=0; isystematic<nsystematics; isystematic++) {

	cout << "Current systematic: " << systematic[isystematic] << endl;

        nsknnd_min[isystematic] = nsknnd_deflt;
        nsknnd_max[isystematic] = nsknnd_deflt;
        nsk_min[isystematic] = nsk_deflt;
        nsk_max[isystematic] = nsk_deflt;

	if(skip[isystematic]) {
            cout << "... skipping " << endl;
	    continue;
        }

	if(generate_random[isystematic]) {
            nsknnd_min[isystematic] = 99999999.;
            nsknnd_max[isystematic] = 0.;
            nsk_min[isystematic] = 99999999.;
            nsk_max[isystematic] = 0.;
            // generate toy MC expts and calculate Nsknnd
            double nsknndexp[nexperiments];
            double nskexp[nexperiments];
            for(int iexperiment = 0; iexperiment < nexperiments; iexperiment++) {
                    if(iexperiment%100==0) {
			cout << "Current experiment: " << iexperiment << endl;
                    }
		    skexp->ResetSystematicVariations();
	     	    skexp->TweakSystematic(systematic[isystematic]);
	  	    skexp->Generate(ph);
  		    double nsknnd = skexp->CurrPredictedErecoSpectrum(s)->Integral();
                    nsknndexp[iexperiment] = nsknnd;
		    double nsk = 0;
		    for (int im = 0; im < conv::knm; im++) {
		      conv::SuperkMode_t m = conv::kSKReactionModes[s][im];
		      if(m == conv::kSKmUndefined) continue;
		      nsk += skexp->GetSuperkMC()->CurrErecoSpectrum(s,m)->Integral();
		    }
                    nskexp[iexperiment] = nsk;
            }//iexperiment
            // calculate standard deviation
            double nsknnd_sigma = 0;
            double nsk_sigma = 0;
            for(int iexperiment = 0; iexperiment < nexperiments; iexperiment++) {
                nsknnd_sigma += TMath::Power(nsknndexp[iexperiment] - nsknnd_deflt, 2.);
                nsk_sigma += TMath::Power(nskexp[iexperiment] - nsk_deflt, 2.);
            }
            nsknnd_sigma = TMath::Sqrt(nsknnd_sigma/nexperiments);
            nsk_sigma = TMath::Sqrt(nsk_sigma/nexperiments);
            // set Nsknnd range
	    nsknnd_min[isystematic] = nsknnd_deflt - nsknnd_sigma;
            nsknnd_max[isystematic] = nsknnd_deflt + nsknnd_sigma;
	    nsk_min[isystematic] = nsk_deflt - nsk_sigma;
            nsk_max[isystematic] = nsk_deflt + nsk_sigma;
	} else {
	    skexp->ResetSystematicVariations();
     	    skexp->TweakSystematic(systematic[isystematic], +1);
  	    skexp->Generate(ph);
  	    double nsknnd1 = skexp->CurrPredictedErecoSpectrum(s)->Integral();
  	    double nsk1 = 0;
	    for (int im = 0; im < conv::knm; im++) {
	      conv::SuperkMode_t m = conv::kSKReactionModes[s][im];
	      if(m == conv::kSKmUndefined) continue;
	      nsk1 += skexp->GetSuperkMC()->CurrErecoSpectrum(s,m)->Integral();
	    }
	    skexp->ResetSystematicVariations();
     	    skexp->TweakSystematic(systematic[isystematic], -1);
  	    skexp->Generate(ph);
  	    double nsknnd2 = skexp->CurrPredictedErecoSpectrum(s)->Integral();
            nsknnd_min[isystematic] = TMath::Min(nsknnd1,nsknnd2);
            nsknnd_max[isystematic] = TMath::Max(nsknnd1,nsknnd2);
  	    double nsk2 = 0;
	    for (int im = 0; im < conv::knm; im++) {
	      conv::SuperkMode_t m = conv::kSKReactionModes[s][im];
	      if(m == conv::kSKmUndefined) continue;
	      nsk2 += skexp->GetSuperkMC()->CurrErecoSpectrum(s,m)->Integral();
	    }
            nsk_min[isystematic] = TMath::Min(nsk1,nsk2);
            nsk_max[isystematic] = TMath::Max(nsk1,nsk2);
	}	
  }//isystematic

  cout << endl;
  cout << "Nsk/Nnd(default) = " << nsknnd_deflt << endl;
  cout << "Nsk(default) = " << nsk_deflt << endl;
  for(int isystematic=0; isystematic<nsystematics; isystematic++) {
        if(skip[isystematic] == true) continue;
	cout << " * " << systematic[isystematic] << " :  \tNsk/Nnd = (" 
	     << nsknnd_min[isystematic] << ", " << nsknnd_max[isystematic] << ") \tNsk = ("
	     << nsk_min[isystematic] << ", " << nsk_max[isystematic] << ")" << endl;
  }
  cout << endl;
  cout << "Fractional change: " << endl;
  for(int isystematic=0; isystematic<nsystematics; isystematic++) {
        if(skip[isystematic] == true) continue;
	cout << " * " << systematic[isystematic] << " :  \tdNnd/Nsk = (" 
             << 100*(nsknnd_max[isystematic]-nsknnd_deflt)/nsknnd_deflt << " %," 
             << 100*(nsknnd_min[isystematic]-nsknnd_deflt)/nsknnd_deflt << " %) \tdNsk = (" 
             << 100*(nsk_max[isystematic]-nsk_deflt)/nsk_deflt << " %," 
             << 100*(nsk_min[isystematic]-nsk_deflt)/nsk_deflt << " %) \td" << endl;
  }
  cout << endl;
  cout << "Total factional change: " << endl;
  double dnsknnd_syst_pos = 0;
  double dnsknnd_syst_neg = 0;
  double dnsk_syst_pos = 0;
  double dnsk_syst_neg = 0;
  for(int isystematic=0; isystematic<nsystematics; isystematic++) {
     double pos = nsknnd_max[isystematic]-nsknnd_deflt;
     double neg = nsknnd_min[isystematic]-nsknnd_deflt;
     dnsknnd_syst_pos += (pos*pos);
     dnsknnd_syst_neg += (neg*neg);
     double pos = nsk_max[isystematic]-nsk_deflt;
     double neg = nsk_min[isystematic]-nsk_deflt;
     dnsk_syst_pos += (pos*pos);
     dnsk_syst_neg += (neg*neg);
  }
  dnsknnd_syst_pos = TMath::Sqrt(dnsknnd_syst_pos);
  dnsknnd_syst_neg = TMath::Sqrt(dnsknnd_syst_neg);
  cout << "dNsk/Nnd = (+" << 100*dnsknnd_syst_pos/nsknnd_deflt << "%, -" << 100*dnsknnd_syst_neg/nsknnd_deflt << ")" << endl;
  dnsk_syst_pos = TMath::Sqrt(dnsk_syst_pos);
  dnsk_syst_neg = TMath::Sqrt(dnsk_syst_neg);
  cout << "dNsk = (+" << 100*dnsk_syst_pos/nsk_deflt << "%, -" << 100*dnsk_syst_neg/nsk_deflt << ")" << endl;

  double dnsknnd_stat = TMath::Sqrt(nsknnd_deflt);
  double dnsk_stat = TMath::Sqrt(nsk_deflt);

  cout << "Nsk/Nnd = " << nsknnd_deflt 
       << " +/- " << dnsknnd_stat << " (stat)"
       << " + " << dnsknnd_syst_pos 
       << " - " << dnsknnd_syst_neg << " (syst)" << endl;
  cout << "Nsk = " << nsk_deflt 
       << " +/- " << dnsk_stat << " (stat)"
       << " + " << dnsk_syst_pos 
       << " - " << dnsk_syst_neg << " (syst)" << endl;

  cout << "+dNsk/Nnd(total) = " << TMath::Sqrt(dnsknnd_stat*dnsknnd_stat + dnsknnd_syst_pos*dnsknnd_syst_pos) << endl;
  cout << "-dNsk/Nnd(total) = " << TMath::Sqrt(dnsknnd_stat*dnsknnd_stat + dnsknnd_syst_neg*dnsknnd_syst_neg) << endl;
  cout << "+dNsk(total) = " << TMath::Sqrt(dnsk_stat*dnsk_stat + dnsk_syst_pos*dnsk_syst_pos) << endl;
  cout << "-dNsk(total) = " << TMath::Sqrt(dnsk_stat*dnsk_stat + dnsk_syst_neg*dnsk_syst_neg) << endl;
}
