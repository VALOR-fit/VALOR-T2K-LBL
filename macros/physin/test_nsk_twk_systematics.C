//
// ------------------------------------------------------------------------------------------------
// Summarize Nsk change as a result of tweaking systematic parameters.
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

void test_nsk_twk_systematics
(
  const char * sksample            = "1Rmu",
  const char * hypothesis          = "3f_osc",
  const double input_mixing23      = 1.00,
  const double input_mixing13      = 0.1,
  const double input_mixing12      = 0.8704,
  const double input_dm23sq        = 2.4E-3,
  const double input_dm12sq        = 7.6E-5,
  const double input_deltacp       = 0,
  const char * config_file         = "analysis_configuration.xml"
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
  cout << "- Physics hypothesis: " << ph << endl;
  cout << *Opt::Instance() << endl;
  
  const int nexperiments = 5000;

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
	false,  /* sk_energy_scale (no effect on Nsk) */
	false,  /* nd280_eff                          */
	false,  /* nd280_norm                         */
	false,  /* flux_norm                          */
	false,  /* flux_shape                         */
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
  }
  else if (s == conv::kSKsFCFV1RingELike) {
    skip[0]  = true;
    skip[1]  = true;
    skip[2]  = true;
    skip[3]  = true;
    skip[4]  = true;
    //skip[16] = true; //xsec_nc_to_ccqe is used for nd280 systematics
  }
  else if (s == conv::kSKsUndefined) {
    cout << "Sample undefined" << endl;
    return;
  }
  else {
    cout << "Sample not supported" << endl;
    return;
  }

  double nsknnd_min[nsystematics];
  double nsknnd_max[nsystematics];
  double nsk_min[nsystematics];
  double nsk_max[nsystematics];
  double nnd_min[nsystematics];
  double nnd_max[nsystematics];

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
  double nnd_deflt = skexp->GetND280Measurement()->NExpCC() + skexp->GetND280Measurement()->NExpNCBkg();
  
  // get predictions with tweaked systematic params
  for(int isystematic=0; isystematic<nsystematics; isystematic++) {

	cout << "Current systematic: " << systematic[isystematic] << endl;

        nsknnd_min[isystematic] = nsknnd_deflt;
        nsknnd_max[isystematic] = nsknnd_deflt;
        nsk_min[isystematic] = nsk_deflt;
        nsk_max[isystematic] = nsk_deflt;
        nnd_min[isystematic] = nnd_deflt;
        nnd_max[isystematic] = nnd_deflt;

	if(skip[isystematic]) {
            cout << "... skipping " << endl;
	    continue;
        }

	if(generate_random[isystematic]) {
            nsknnd_min[isystematic] = 99999999.;
            nsknnd_max[isystematic] = 0.;
            nsk_min[isystematic] = 99999999.;
            nsk_max[isystematic] = 0.;
            nnd_min[isystematic] = 99999999.;
            nnd_max[isystematic] = 0.;
            // generate toy MC expts and calculate Nsknnd
            double nsknndexp[nexperiments];
            double nskexp[nexperiments];
            double nndexp[nexperiments];
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
  		    double nnd = skexp->GetND280Measurement()->NExpCC() + skexp->GetND280Measurement()->NExpNCBkg();
                    nndexp[iexperiment] = nnd;
            }//iexperiment
            // calculate standard deviation
            double nsknnd_sigma = 0;
            double nsk_sigma = 0;
            double nnd_sigma = 0;
            for(int iexperiment = 0; iexperiment < nexperiments; iexperiment++) {
                nsknnd_sigma += TMath::Power(nsknndexp[iexperiment] - nsknnd_deflt, 2.);
                nsk_sigma += TMath::Power(nskexp[iexperiment] - nsk_deflt, 2.);
                nnd_sigma += TMath::Power(nndexp[iexperiment] - nnd_deflt, 2.);
            }
            nsknnd_sigma = TMath::Sqrt(nsknnd_sigma/nexperiments);
            nsk_sigma = TMath::Sqrt(nsk_sigma/nexperiments);
            nnd_sigma = TMath::Sqrt(nnd_sigma/nexperiments);
            // set Nsknnd range
	    nsknnd_min[isystematic] = nsknnd_deflt - nsknnd_sigma;
            nsknnd_max[isystematic] = nsknnd_deflt + nsknnd_sigma;
	    nsk_min[isystematic] = nsk_deflt - nsk_sigma;
            nsk_max[isystematic] = nsk_deflt + nsk_sigma;
	    nnd_min[isystematic] = nnd_deflt - nnd_sigma;
            nnd_max[isystematic] = nnd_deflt + nnd_sigma;
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
  	    double nnd1 = skexp->GetND280Measurement()->NExpCC() + skexp->GetND280Measurement()->NExpNCBkg();
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
  	    double nnd2 = skexp->GetND280Measurement()->NExpCC() + skexp->GetND280Measurement()->NExpNCBkg();
            nnd_min[isystematic] = TMath::Min(nnd1,nnd2);
            nnd_max[isystematic] = TMath::Max(nnd1,nnd2);
	}	
  }//isystematic

  cout << endl;
  cout << "Nsk/Nnd(default) = " << nsknnd_deflt << endl;
  cout << "Nsk(default) = " << nsk_deflt << endl;
  cout << "Nnd(default) = " << nnd_deflt << endl;
  for(int isystematic=0; isystematic<nsystematics; isystematic++) {
        if(skip[isystematic] == true) continue;
	cout << " * " << systematic[isystematic] << " :  \tNsk/Nnd = (" 
	     << nsknnd_min[isystematic] << ", " << nsknnd_max[isystematic] << ") \tNsk = ("
	     << nsk_min[isystematic] << ", " << nsk_max[isystematic] << ") \tNnd = ("
	     << nnd_min[isystematic] << ", " << nnd_max[isystematic] << ")" << endl;
  }
  cout << endl;
  cout << "Fractional change: " << endl;
  for(int isystematic=0; isystematic<nsystematics; isystematic++) {
        if(skip[isystematic] == true) continue;
	cout << " * " << systematic[isystematic] << " :  \tdNnd/Nsk = (" 
             << 100*(nsknnd_max[isystematic]-nsknnd_deflt)/nsknnd_deflt << " %," 
             << 100*(nsknnd_min[isystematic]-nsknnd_deflt)/nsknnd_deflt << " %) \tdNsk = (" 
             << 100*(nsk_max[isystematic]-nsk_deflt)/nsk_deflt << " %," 
             << 100*(nsk_min[isystematic]-nsk_deflt)/nsk_deflt << " %) \tdNnd = (" 
             << 100*(nnd_max[isystematic]-nnd_deflt)/nnd_deflt << " %," 
             << 100*(nnd_min[isystematic]-nnd_deflt)/nnd_deflt << " %)" 
             << endl;
  }
  cout << endl;
  cout << "Total factional change: " << endl;
  double dnsknnd_syst_pos = 0;
  double dnsknnd_syst_neg = 0;
  double dnsk_syst_pos = 0;
  double dnsk_syst_neg = 0;
  double dnnd_syst_pos = 0;
  double dnnd_syst_neg = 0;
  for(int isystematic=0; isystematic<nsystematics; isystematic++) {
     double pos = nsknnd_max[isystematic]-nsknnd_deflt;
     double neg = nsknnd_min[isystematic]-nsknnd_deflt;
     dnsknnd_syst_pos += (pos*pos);
     dnsknnd_syst_neg += (neg*neg);
     double pos = nsk_max[isystematic]-nsk_deflt;
     double neg = nsk_min[isystematic]-nsk_deflt;
     dnsk_syst_pos += (pos*pos);
     dnsk_syst_neg += (neg*neg);
     double pos = nnd_max[isystematic]-nnd_deflt;
     double neg = nnd_min[isystematic]-nnd_deflt;
     dnnd_syst_pos += (pos*pos);
     dnnd_syst_neg += (neg*neg);
  }
  dnsknnd_syst_pos = TMath::Sqrt(dnsknnd_syst_pos);
  dnsknnd_syst_neg = TMath::Sqrt(dnsknnd_syst_neg);
  cout << "dNsk/Nnd = (+" << 100*dnsknnd_syst_pos/nsknnd_deflt << "%, -" << 100*dnsknnd_syst_neg/nsknnd_deflt << ")" << endl;
  dnsk_syst_pos = TMath::Sqrt(dnsk_syst_pos);
  dnsk_syst_neg = TMath::Sqrt(dnsk_syst_neg);
  cout << "dNsk = (+" << 100*dnsk_syst_pos/nsk_deflt << "%, -" << 100*dnsk_syst_neg/nsk_deflt << ")" << endl;
  dnnd_syst_pos = TMath::Sqrt(dnnd_syst_pos);
  dnnd_syst_neg = TMath::Sqrt(dnnd_syst_neg);
  cout << "dNnd = (+" << 100*dnnd_syst_pos/nnd_deflt << "%, -" << 100*dnnd_syst_neg/nnd_deflt << ")" << endl;

  double dnsknnd_stat = TMath::Sqrt(nsknnd_deflt);
  double dnsk_stat = TMath::Sqrt(nsk_deflt);
  double dnnd_stat = TMath::Sqrt(nnd_deflt);

  cout << "Nsk/Nnd = " << nsknnd_deflt 
       << " +/- " << dnsknnd_stat << " (stat)"
       << " + " << dnsknnd_syst_pos 
       << " - " << dnsknnd_syst_neg << " (syst)" << endl;
  cout << "Nsk = " << nsk_deflt 
       << " +/- " << dnsk_stat << " (stat)"
       << " + " << dnsk_syst_pos 
       << " - " << dnsk_syst_neg << " (syst)" << endl;
  cout << "Nnd = " << nnd_deflt 
       << " +/- " << dnnd_stat << " (stat)"
       << " + " << dnnd_syst_pos 
       << " - " << dnnd_syst_neg << " (syst)" << endl;

  cout << "+dNsk/Nnd(total) = " << TMath::Sqrt(dnsknnd_stat*dnsknnd_stat + dnsknnd_syst_pos*dnsknnd_syst_pos) << endl;
  cout << "-dNsk/Nnd(total) = " << TMath::Sqrt(dnsknnd_stat*dnsknnd_stat + dnsknnd_syst_neg*dnsknnd_syst_neg) << endl;
  cout << "+dNsk(total) = " << TMath::Sqrt(dnsk_stat*dnsk_stat + dnsk_syst_pos*dnsk_syst_pos) << endl;
  cout << "-dNsk(total) = " << TMath::Sqrt(dnsk_stat*dnsk_stat + dnsk_syst_neg*dnsk_syst_neg) << endl;
  cout << "+dNnd(total) = " << TMath::Sqrt(dnnd_stat*dnnd_stat + dnnd_syst_pos*dnnd_syst_pos) << endl;
  cout << "-dNnd(total) = " << TMath::Sqrt(dnnd_stat*dnnd_stat + dnnd_syst_neg*dnnd_syst_neg) << endl;
}
