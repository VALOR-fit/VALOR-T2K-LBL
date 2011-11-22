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
  double sinsq2theta13 = 0,
  double sinsq2theta23 = 1,
  double sinsq2theta12 = 0.8704,
  double dm23sq = 0.0024,
  double dm12sq = 0.000076,
  double deltacp = 0,
  t2kana::conv::NuMassHierarchy_t hier = t2kana::conv::kNmhNormal, //or t2kana::conv::kNMhInverted
  const char * config_file = "analysis_configuration.xml"
)
{
  using namespace std;
  using namespace t2kana;

  const char * config = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);
  cout << "Reading analysis configuration from: " << config << endl;

  if(hier == t2kana::conv::kNmhInverted)
    dm23sq = -dm23sq;

  // superk spectrum predictor
  SuperkExpectation skexp;
  skexp.Read(config);

  // fit params
  PhysHypothesis ph(t2kana::conv::kHypOsc3F);

  // generate predictions
  ph.SetCurrValue(t2kana::conv::kPhypSinsq2Theta23, sinsq2theta23);
  ph.SetCurrValue(t2kana::conv::kPhypSinsq2Theta13, sinsq2theta13);
  ph.SetCurrValue(t2kana::conv::kPhypSinsq2Theta12, sinsq2theta12);
  ph.SetCurrValue(t2kana::conv::kPhypDm23sq,        dm23sq);
  ph.SetCurrValue(t2kana::conv::kPhypDm12sq,        dm12sq);
  ph.SetCurrValue(t2kana::conv::kPhypDeltaCP,       deltacp);

  // generate predictions
  skexp.Generate(ph);


  // integrate Ereco distribution
  double nsk = skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike)->Integral();
  double nskmucc = skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuMuCCQE)->Integral();
  nskmucc += skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuMuCC1pi)->Integral();
  nskmucc += skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuMuCCoth_mu)->Integral();
  nskmucc += skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuMuBarCCQE)->Integral();
  nskmucc += skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuMuBarCC1pi)->Integral();
  nskmucc += skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuMuBarCCoth_mu)->Integral();
  double nsknc = skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuMuNC)->Integral();
  nskmucc += skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuMuBarNC)->Integral();
  nskmucc += skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuENC)->Integral();
  double nskebeamcc = skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuECCQE)->Integral();
  nskebeamcc += skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuECC1pi)->Integral();
  nskebeamcc += skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmNuECCoth_mu)->Integral();
  double nskesigcc = skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmSigNuECCQE)->Integral();
  nskesigcc += skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmSigNuECC1pi)->Integral();
  nskesigcc += skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike, conv::kSKmSigNuECCoth)->Integral();


  //output in form of table 16, TN52
  cout << "Hierarchy is " << t2kana::util::AsString(hier) << " dm23sq is " << dm23sq << endl;
  cout << "     " << nskmucc+nskebeamcc+nsknc << " beam total" << endl;
  cout << "     " << nskmucc << " numu+numubar CC" << endl;
  cout << "     " << nskebeamcc << " nue beam CC" << endl;
  cout << "     " << nsknc << " all NC" << endl;
  cout << "     " << nskesigcc << " nue signal CC" << endl;
  cout << "Number of predicted events at SK is: " << nsk << "     = " << nskmucc+nsknc+nskebeamcc+nskesigcc << endl;

}
