//
// ------------------------------------------------------------------------------------------------
// Generate total beam error tables for nue appearance analysis
// To check with Josh's tables
// These tables may be used as input instead of calculating the error in each point
//
// Lorena
// ------------------------------------------------------------------------------------------------
//
//

#include <iomanip.h>

using namespace t2kana; //need to put this here, to allow a SuperkSample_t as input
void test_nsk_twk_beam_norm_nue
(
  const string sample_str                 = "1Re",                          /* options: "1Rmu", "1Re" */
  const char * hypothesis                 = "3f_osc",                       /* options: "2f_numu_disap", "3f_osc" */
  double       input_mixing12             = -99,                            /* */
  double       input_mixing13             = -99,                            /* */
  double       input_mixing23             = -99,                            /* */
  double       input_deltacp              = -99,                            /* */
  double       input_dm12sq               = -99,                            /* */
  double       input_dm23sq               = -99,                            /* */
  const char * config_file                = "analysis_configuration.xml"    /* */
)
{
  using namespace std;
  using namespace t2kana;

  conv::SuperkSample_t s = util::SKSampleFromString(sample_str);
  int sidx = util::ArrIdx::Instance()->SuperkSampleIdx(s);
  if(sidx == -1) {
    cerr << "Sample " << sample_str << " not defined. Check conv::kSKSamples[kns]" << endl;
    exit(1);
  }
  
  // set physics hypothesis
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
    ph.SetCurrValue(conv::kPhypDm23sq,         input_dm23sq);
  } else
  if(ph.Hypothesis() == conv::kHypOsc3F) {
    ph.SetCurrValue(conv::kPhypMixing12,       input_mixing12);
    ph.SetCurrValue(conv::kPhypMixing13,       input_mixing13);
    ph.SetCurrValue(conv::kPhypMixing23,       input_mixing23);
    ph.SetCurrValue(conv::kPhypDeltaCP,        input_deltacp);
    ph.SetCurrValue(conv::kPhypDm12sq,         input_dm12sq);
    ph.SetCurrValue(conv::kPhypDm23sq,         input_dm23sq);
  }

  // create & config superk spectrum predictor
  const char * config = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);
  cout << "Reading analysis configuration from: " << config << endl;

  // Global options
  Opt::Instance()->gOscMixingParamChoice = conv::kNuMixSinsq2Theta;
  Opt::Instance()->gMassHierarchyAssumed = conv::kNmhNormal;
  
  // Print options
  cout << "\n[Job inputs]" << endl;
  cout << "- Reading analysis configuration from: " << config << endl;
  cout << "- SuperK sample: " << util::AsString(s) << endl;
  cout << "- Physics hypothesis: " << ph << endl;
  cout << *Opt::Instance() << endl;

  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read(config);

  //FIRST - DELTACP VS MIXING13

  // define the grid
  // deltacp, mixing13 range and step
  double mixing13_min        = 0.0;
  double mixing13_max        = 0.4;
  double mixing13_step       = 0.004;
  double deltacp_min         = -1 * conv::kpi;
  double deltacp_max         = conv::kpi;
  double deltacp_step        = 2*conv::kpi/90.;
  // make sure that min and max are at the bin centre and not at the low/high edge of the first/last bin
  int mixing13_nbins      = 1 + (int) ((mixing13_max - mixing13_min) / mixing13_step);
  int deltacp_nbins       = 1 + (int) ((deltacp_max       - deltacp_min      ) / deltacp_step      );
  
  // book the histogram for the table
  TH2D *  mixing13_deltacp_table = new TH2D("mixing13_deltacp_table", "",
       mixing13_nbins, 
       mixing13_min - 0.5*mixing13_step, 
       mixing13_max + 0.5*mixing13_step,
       deltacp_nbins, 
       deltacp_min - 0.5*deltacp_step, 
       deltacp_max + 0.5*deltacp_step);

  //loop over points in grid

  for(int i = 1; i <= mixing13_nbins; i++ ){
    for (int j =1; j <= deltacp_nbins; j++) {
  //for(int i = 1; i <= 3; i++ ){
  //for (int j =1; j <= 3; j++) {
      
      input_mixing13 = mixing13_deltacp_table->GetXaxis()->GetBinCenter(i);
      input_deltacp       = mixing13_deltacp_table->GetYaxis()->GetBinCenter(j);
      
      ph.SetCurrValue(conv::kPhypMixing13,       input_mixing13);
      ph.SetCurrValue(conv::kPhypDeltaCP,        input_deltacp);
      

      //
      // nominal
      //
      skexp->Generate(ph);
      //sk total * nd280 ratio
      double nsk_0                 = skexp->CurrPredictedErecoSpectrum(s)->Integral();
      
      //
      // +1 sigma // it is the same with -1 sigma
      //
      skexp->TweakSystematic("flux_norm",+1.);
      skexp->Generate(ph);
      //sk total * nd280 ratio
      double nsk                 = skexp->CurrPredictedErecoSpectrum(s)->Integral();
      
      //fraction of change (percentile)
      double frac_nsk            = (nsk_0 > 0) ? 100.*(nsk - nsk_0) / nsk_0 : 0.;
      
      //Now store the fraction
      mixing13_deltacp_table->SetBinContent(i,j,frac_nsk);
      
      skexp->ResetSystematicVariations();
      
    }
  }//end of loop
  
  gStyle->SetPalette(1);
  mixing13_deltacp_table->Draw("colz");
  
  //now write the histogram
  TFile f("test_table.root", "recreate"); 
  mixing13_deltacp_table->Write("table");
  f.Close();
  
  //now do in log scale for X axis (?)

  //  mixing13_deltacp_table->GetXaxis()->SetLog();
  //c1->Update()

  cout << "Done!" << endl;

}
