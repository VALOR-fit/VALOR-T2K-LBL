//
// ------------------------------------------------------------------------------------------------
// Print out the number of events predicted at SuperK (NSK)
// Also print out the number of events for certain categories. Using the same categories as
// table 16 in TN52 - which gives the results for run1+2 nue analysis B
//
// Tom Dealtry
// ------------------------------------------------------------------------------------------------
//

//#include "Conventions.h"

void test_nsk
(
 double input_mixing13    = -99,
 double input_mixing23    = -99,
 double input_mixing12    = -99,
 double input_dm23sq      = -99,
 double input_dm12sq      = -99,
 double input_deltacp     = -99,
 const string sample_str  = "1Rmu",
 const char * config_file = "analysis_configuration.xml"
)
{
  using namespace std;
  using namespace t2kana;

  conv::SuperkSample_t s = util::SKSampleFromString(sample_str);
  int is = util::ArrIdx::Instance()->SuperkSampleIdx(s);

  const char * config = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);
  cout << "Reading analysis configuration from: " << config << endl;

  // superk spectrum predictor
  SuperkExpectation skexp;
  skexp.Read(config);

  // fit params
  PhysHypothesis ph(t2kana::conv::kHypOsc3F);

  // generate predictions
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

  ph.SetCurrValue(t2kana::conv::kPhypMixing23, input_mixing23);
  ph.SetCurrValue(t2kana::conv::kPhypMixing13, input_mixing13);
  ph.SetCurrValue(t2kana::conv::kPhypMixing12, input_mixing12);
  ph.SetCurrValue(t2kana::conv::kPhypDm23sq,   input_dm23sq);
  ph.SetCurrValue(t2kana::conv::kPhypDm12sq,   input_dm12sq);
  ph.SetCurrValue(t2kana::conv::kPhypDeltaCP,  input_deltacp);

  // generate predictions
  skexp.Generate(ph);

  // integrate Ereco distribution
  double nsk = skexp.CurrPredictedErecoSpectrum(s)->Integral();
  double nskmc = 0, nskmucc = 0, nskmubarcc = 0, nskecc = 0, nskebarcc = 0, nsksigecc = 0, nsknc = 0;
  for(int im = 0; im < conv::knm; im++) {
    conv::SuperkMode_t m = conv::kSKReactionModes[is][im];
    if(m == conv::kSKmUndefined)
      continue;
    nskmc += skexp.GetSuperkMC()->Template(s, m, conv::kSKMCtFin)->Integral();

    //split into categories (all NC, CC numu, CC numubar, CC nue, CC nuebar, CC signue)
    if(util::IsNC(m))
      nsknc += skexp.GetSuperkMC()->Template(s, m, conv::kSKMCtFin)->Integral();
    else {
      if(util::IsNuMu(m))
	nskmucc += skexp.GetSuperkMC()->Template(s, m, conv::kSKMCtFin)->Integral();
      else if(util::IsNuMuBar(m))
	nskmubarcc += skexp.GetSuperkMC()->Template(s, m, conv::kSKMCtFin)->Integral();
      else if(util::IsNuE(m))
	nskecc += skexp.GetSuperkMC()->Template(s, m, conv::kSKMCtFin)->Integral();
      else if(util::IsNuEBar(m))
	nskebarcc += skexp.GetSuperkMC()->Template(s, m, conv::kSKMCtFin)->Integral();
      else if(util::IsSignalNuE(m))
	nsksigecc += skexp.GetSuperkMC()->Template(s, m, conv::kSKMCtFin)->Integral();
    }//CC?
  }

  //output in form of table 16, TN52
  cout << endl << " dm23sq is " << input_dm23sq << endl;
  cout << "     " << nskmucc + nskmubarcc + nskecc + nskebarcc + nsksigecc + nsknc << " beam total" << endl;
  cout << "     " << nskmucc + nskmubarcc << " numu+numubar CC" << endl;
  cout << "     " << nskecc + nskebarcc << " nue+nuebar beam CC" << endl;
  cout << "     " << nsknc << " all NC" << endl;
  cout << "     " << nsksigecc << " nue signal CC" << endl;

  //output a bit more broken down
  cout << endl << " dm23sq is " << input_dm23sq << endl;
  cout << "     " << nskmucc + nskmubarcc + nskecc + nskebarcc + nsksigecc + nsknc << " beam total" << endl;
  cout << "     " << nskmucc << " numu CC" << endl;
  cout << "     " << nskmubarcc << " numubar CC" << endl;
  cout << "     " << nskecc << " nue beam CC" << endl;
  cout << "     " << nskebarcc << " nuebar beam CC" << endl;
  cout << "     " << nsknc << " all NC" << endl;
  cout << "     " << nsksigecc << " nue signal CC" << endl;

  //the total
  cout << endl << "Number of predicted events at SK is: " << nsk << "     = " 
       << nskmucc + nskmubarcc + nskecc + nskebarcc + nsksigecc + nsknc << "     = " 
       << nskmc << endl;
}
