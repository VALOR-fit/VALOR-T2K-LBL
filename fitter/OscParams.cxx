//____________________________________________________________________________
/*
  \copyright (c) 2010-2014, The VALOR group.

  >> ChangeLog:

*/
//____________________________________________________________________________

#include <cstdlib>

#include "OscParams.h"

#include "UtlForMacros.h"
#include "T2KPhysParamId.h"
#include "Conventions/Opt.h"
#include "Conventions/PhysHypothesisId.h"
#include "Conventions/Constants.h"
#include "Utilities/UtlNumStat.h"

using std::cout;
using std::cerr;
using std::endl;
using std::string;

using namespace valor;
using namespace valor::conv;
  
//____________________________________________________________________________
namespace valor {
  namespace util {
    ostream & operator << (ostream & stream, const OscParams & oscparams)
    {
      oscparams.Print(stream);
      return stream;
    }
  }
}
//____________________________________________________________________________
valor::util::OscParams * valor::util::OscParams::fInstance = 0;
//____________________________________________________________________________
valor::util::OscParams::OscParams()
{
  fInstance = 0;

  this->Init();
}
//____________________________________________________________________________
valor::util::OscParams::~OscParams()
{
  fInstance = 0;

  this->Cleanup();
}
//____________________________________________________________________________
valor::util::OscParams * valor::util::OscParams::Instance()
{
  if(fInstance == 0) {
    static valor::util::OscParams::Cleaner cleaner;
    cleaner.DummyMethodAndSilentCompiler();
    fInstance = new OscParams;
  }
  return fInstance;
}
//____________________________________________________________________________
void valor::util::OscParams::Init()
{
  hPrior_th13_dcp_nh = 0;
  hPrior_th13_dcp_ih = 0;
  hPrior_23_nh = 0;
  hPrior_23_ih = 0;
  hToy_th13_dcp_nh = 0;
  hToy_th13_dcp_ih = 0;
  hToy_23_nh = 0;
  hToy_23_ih = 0;
  hToy_MH = 0;
}
//____________________________________________________________________________
void valor::util::OscParams::Cleanup()
{
  if(hPrior_th13_dcp_nh) delete hPrior_th13_dcp_nh;
  if(hPrior_th13_dcp_ih) delete hPrior_th13_dcp_ih;
  if(hPrior_23_nh) delete hPrior_23_nh;
  if(hPrior_23_ih) delete hPrior_23_ih;
  if(hToy_th13_dcp_nh) delete hToy_th13_dcp_nh;
  if(hToy_th13_dcp_ih) delete hToy_th13_dcp_ih;
  if(hToy_23_nh) delete hToy_23_nh;
  if(hToy_23_ih) delete hToy_23_ih;
  if(hToy_MH) delete hToy_MH;
}
//____________________________________________________________________________
void valor::util::OscParams::SetValues(double input_mixing12_nh, double input_mixing13_nh, double input_mixing23_nh,
				       double input_dcp_nh,      double input_dm12sq_nh,   double input_dm23sq_nh,
				       double input_beta_bar,
				       double input_mixing12_ih, double input_mixing13_ih, double input_mixing23_ih,
				       double input_dcp_ih,      double input_dm12sq_ih,   double input_dm23sq_ih,
				       double input_mh_nh_prior, bool marginalize_MH,
				       string file_priors)
{  
  // Check whether some osc. params need to be randomized (i.e. if == -999)
  //NH
  randomize_mixing12_nh = (TMath::Abs(input_mixing12_nh + 999) < kASmallNum);
  randomize_mixing13_nh = (TMath::Abs(input_mixing13_nh + 999) < kASmallNum);
  randomize_mixing23_nh = (TMath::Abs(input_mixing23_nh + 999) < kASmallNum);
  randomize_dcp_nh      = (TMath::Abs(input_dcp_nh      + 999) < kASmallNum);
  randomize_dm12sq_nh   = (TMath::Abs(input_dm12sq_nh   + 999) < kASmallNum);
  randomize_dm23sq_nh   = (TMath::Abs(input_dm23sq_nh   + 999) < kASmallNum);
  //IH
  randomize_mixing12_ih = (TMath::Abs(input_mixing12_ih + 999) < kASmallNum);
  randomize_mixing13_ih = (TMath::Abs(input_mixing13_ih + 999) < kASmallNum);
  randomize_mixing23_ih = (TMath::Abs(input_mixing23_ih + 999) < kASmallNum);
  randomize_dcp_ih      = (TMath::Abs(input_dcp_ih      + 999) < kASmallNum);
  randomize_dm12sq_ih   = (TMath::Abs(input_dm12sq_ih   + 999) < kASmallNum);
  randomize_dm23sq_ih   = (TMath::Abs(input_dm23sq_ih   + 999) < kASmallNum);
  //Beta bar
  randomize_beta_bar    = (TMath::Abs(input_beta_bar    + 999) < kASmallNum);

  //Set to defaults (only if enough negative (standard value is -99)
  //NH
  if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsq2Theta) {
    if(input_mixing13_nh < -1 && input_mixing13_nh > -900) input_mixing13_nh = conv::kGlobFitSinsq2Theta13;
    if(input_mixing23_nh < -1 && input_mixing23_nh > -900) input_mixing23_nh = conv::kGlobFitSinsq2Theta23;
    if(input_mixing12_nh < -1 && input_mixing12_nh > -900) input_mixing12_nh = conv::kGlobFitSinsq2Theta12;
  }
  else {
    if(input_mixing13_nh < -1 && input_mixing13_nh > -900) input_mixing13_nh = conv::kGlobFitSinsqTheta13;
    if(input_mixing23_nh < -1 && input_mixing23_nh > -900) input_mixing23_nh = conv::kGlobFitSinsqTheta23;
    if(input_mixing12_nh < -1 && input_mixing12_nh > -900) input_mixing12_nh = conv::kGlobFitSinsqTheta12;
  }
  if(input_dm23sq_nh < -1  && input_dm23sq_nh > -900) input_dm23sq_nh = conv::kGlobFitDm23sq;
  if(input_dm12sq_nh < -1  && input_dm12sq_nh > -900) input_dm12sq_nh = conv::kGlobFitDm12sq;
  if(input_dcp_nh    < -90 && input_dcp_nh    > -900) input_dcp_nh    = conv::kGlobFitDeltaCP;
  //IH
  if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsq2Theta) {
    if(input_mixing13_ih < -1 && input_mixing13_ih > -900) input_mixing13_ih = conv::kGlobFitSinsq2Theta13;
    if(input_mixing23_ih < -1 && input_mixing23_ih > -900) input_mixing23_ih = conv::kGlobFitSinsq2Theta23;
    if(input_mixing12_ih < -1 && input_mixing12_ih > -900) input_mixing12_ih = conv::kGlobFitSinsq2Theta12;
  }
  else {
    if(input_mixing13_ih < -1 && input_mixing13_ih > -900) input_mixing13_ih = conv::kGlobFitSinsqTheta13;
    if(input_mixing23_ih < -1 && input_mixing23_ih > -900) input_mixing23_ih = conv::kGlobFitSinsqTheta23;
    if(input_mixing12_ih < -1 && input_mixing12_ih > -900) input_mixing12_ih = conv::kGlobFitSinsqTheta12;
  }
  if(input_dm23sq_ih < -1  && input_dm23sq_ih > -900) input_dm23sq_ih = conv::kGlobFitDm23sq;
  if(input_dm12sq_ih < -1  && input_dm12sq_ih > -900) input_dm12sq_ih = conv::kGlobFitDm12sq;
  if(input_dcp_ih    < -90 && input_dcp_ih    > -900) input_dcp_ih = conv::kGlobFitDeltaCP;
  //Beta bar
  if(input_beta_bar < -90 && input_beta_bar > -900) input_beta_bar = conv::kGlobFitBetaBar;
  //MH
  if(input_mh_nh_prior < -1 && input_mh_nh_prior > -900) input_mh_nh_prior = 0.684; //Mach3
  cerr << "Default value for input_mh_nh_prior is hard coded. This should be moved into the relevant Conventions file" << endl;
  
  //Set the input values in the OscParams
  //NH
  mixing12_nh = input_mixing12_nh;
  mixing13_nh = input_mixing13_nh;
  mixing23_nh = input_mixing23_nh;
  dcp_nh      = input_dcp_nh;
  dm12sq_nh   = input_dm12sq_nh;
  dm23sq_nh   = input_dm23sq_nh;
  //IH
  mixing12_ih = input_mixing12_ih;
  mixing13_ih = input_mixing13_ih;
  mixing23_ih = input_mixing23_ih;
  dcp_ih      = input_dcp_ih;
  dm12sq_ih   = input_dm12sq_ih;
  dm23sq_ih   = input_dm23sq_ih;
  //Other
  beta_bar       = input_beta_bar;
  marginalize_mh = marginalize_MH;
  mh_nh_prior    = input_mh_nh_prior;
  MH             = Opt::Instance()->gMassHierarchyAssumed;
  iMH            = (MH == kNmhNormal) ? 0 : 1;

  //Get the prior histograms
  TFile * fPrior = new TFile(file_priors.c_str(), "READ");
  if(!fPrior || !fPrior->IsOpen()){
    cerr << "Error: the file " << file_priors << " is not open!" << endl;
    exit(1);
  }
  hPrior_th13_dcp_nh = new TH2D(*GetTH2DFromTFile(fPrior, "hPrior_th13_dcp_nh"));
  hPrior_th13_dcp_ih = new TH2D(*GetTH2DFromTFile(fPrior, "hPrior_th13_dcp_ih"));
  hPrior_23_nh       = new TH2D(*GetTH2DFromTFile(fPrior, "hPrior_L23_nh"));
  hPrior_23_ih       = new TH2D(*GetTH2DFromTFile(fPrior, "hPrior_L23_ih"));

  //Setup the 'values used' histograms
  //MH
  hToy_MH = new TH1D("hToy_MH", "0-->NH, 1-->IH", 2, 0, 2);
  hToy_MH->GetXaxis()->SetBinLabel(1, "NH");
  hToy_MH->GetXaxis()->SetBinLabel(2, "IH");
  //nh
  if(hPrior_th13_dcp_nh) {
    hToy_th13_dcp_nh = new TH2D(*(hPrior_th13_dcp_nh));
    hToy_th13_dcp_nh->Reset();
    hToy_th13_dcp_nh->SetName("hToy_th13_dcp_nh");
  }
  if(hPrior_23_nh) {
    hToy_23_nh = new TH2D(*hPrior_23_nh);
    hToy_23_nh->Reset();
    hToy_23_nh->SetName("hToy_23_nh");
  }
  //ih
  if(hPrior_th13_dcp_ih) {
    hToy_th13_dcp_ih = new TH2D(*hPrior_th13_dcp_ih);
    hToy_th13_dcp_ih->Reset();
    hToy_th13_dcp_ih->SetName("hToy_th13_dcp_ih");
  }
  if(hPrior_23_ih) {
    hToy_23_ih = new TH2D(*hPrior_23_ih);
    hToy_23_ih->Reset();
    hToy_23_ih->SetName("hToy_23_ih");
  }
}
//____________________________________________________________________________
void valor::util::OscParams::Randomize()
{
  //NH
  if (randomize_mixing12_nh) { mixing12_nh = util::RandomizeMixing12(); }
  if (randomize_dm12sq_nh  ) { dm12sq_nh   = util::RandomizeDm12sq();   }
  //theta13 vs dcp
  if(hPrior_th13_dcp_nh && randomize_dcp_nh && randomize_mixing13_nh) {
    hPrior_th13_dcp_nh->GetRandom2(mixing13_nh, dcp_nh);
  }
  else {
    if (randomize_mixing13_nh) { mixing13_nh = util::RandomizeMixing13(); }
    if (randomize_dcp_nh     ) { dcp_nh      = util::RandomizeCPV();      }
  }
  //theta23 vs dm32sq
  if(hPrior_23_nh && randomize_dm23sq_nh && randomize_mixing23_nh) {
    hPrior_23_nh->GetRandom2(mixing23_nh, dm23sq_nh);
  }
  else {
    if (randomize_mixing23_nh) { mixing23_nh = util::RandomizeMixing23(); }
    if (randomize_dm23sq_nh  ) { dm23sq_nh   = util::RandomizeDm23sq();   }
  }

  //IH
  if (randomize_mixing12_ih) { mixing12_ih = util::RandomizeMixing12(); }
  if (randomize_dm12sq_ih  ) { dm12sq_ih   = util::RandomizeDm12sq();   }
  //theta13 vs dcp
  if(hPrior_th13_dcp_ih && randomize_dcp_ih && randomize_mixing13_ih) {
    hPrior_th13_dcp_ih->GetRandom2(mixing13_ih, dcp_ih);
  }
  else {
    if (randomize_mixing13_ih) { mixing13_ih = util::RandomizeMixing13(); }
    if (randomize_dcp_ih     ) { dcp_ih      = util::RandomizeCPV();      }
  }
  //theta23 vs dm32sq
  if(hPrior_23_ih && randomize_dm23sq_ih && randomize_mixing23_ih) {
    hPrior_23_ih->GetRandom2(mixing23_ih, dm23sq_ih);
  }
  else {
    if (randomize_mixing23_ih) { mixing23_ih = util::RandomizeMixing23(); }
    if (randomize_dm23sq_ih  ) { dm23sq_ih   = util::RandomizeDm23sq();   }
  }

  //other
  if (randomize_beta_bar   ) { beta_bar    = util::RandomizeBetaBar();  }
}
//____________________________________________________________________________
void valor::util::OscParams::UpdatePhysHypothesis(PhysHypothesis& ph, bool print)
{
  double input_mixing12_used, input_mixing13_used, input_mixing23_used;
  double input_dcp_used,      input_dm12sq_used,   input_dm23sq_used;
  double input_beta_bar_used;

  //Get the values to use (based on MH)
  if(iMH == 0) {
    //NH
    input_mixing12_used = mixing12_nh;
    input_mixing13_used = mixing13_nh;
    input_mixing23_used = mixing23_nh;
    input_dcp_used      = dcp_nh;
    input_dm12sq_used   = dm12sq_nh;
    input_dm23sq_used   = dm23sq_nh;
  }
  else if(iMH == 1) {
    //IH
    input_mixing12_used = mixing12_ih;
    input_mixing13_used = mixing13_ih;
    input_mixing23_used = mixing23_ih;
    input_dcp_used      = dcp_ih;
    input_dm12sq_used   = dm12sq_ih;
    input_dm23sq_used   = dm23sq_ih;
  }
  else {
    cerr << "SetOscParams has unknown value of iMH: " << iMH << endl;
    exit(1);
  }
  //Other
  input_beta_bar_used = beta_bar;
  
  //Update the PhysHypothesis
  if(ph.Hypothesis() == kHypNuMuDisapp2F) {
    ph.SetCurrValue(kPhypMixing23,  input_mixing23_used);
    ph.SetCurrValue(kPhypDm23sq,    input_dm23sq_used);
  }//kHypNuMuDisapp2F
  else if(ph.Hypothesis() == kHypOsc3F) {  
    ph.SetCurrValue(kPhypMixing12,  input_mixing12_used);
    ph.SetCurrValue(kPhypMixing13,  input_mixing13_used);
    ph.SetCurrValue(kPhypMixing23,  input_mixing23_used);
    ph.SetCurrValue(kPhypCPV    ,   input_dcp_used);
    ph.SetCurrValue(kPhypDm12sq,    input_dm12sq_used);
    ph.SetCurrValue(kPhypDm23sq,    input_dm23sq_used);
  }//kHypOsc3F
  else if(ph.Hypothesis() == kHypOsc3FBar) {
    ph.SetCurrValue(kPhypMixing12,  input_mixing12_used);
    ph.SetCurrValue(kPhypMixing13,  input_mixing13_used);
    ph.SetCurrValue(kPhypMixing23,  input_mixing23_used);
    ph.SetCurrValue(kPhypCPV    ,   input_dcp_used);
    ph.SetCurrValue(kPhypDm12sq,    input_dm12sq_used);
    ph.SetCurrValue(kPhypDm23sq,    input_dm23sq_used);
       
    ph.SetCurrValue(kPhypMixing12Bar,  input_mixing12_used);
    ph.SetCurrValue(kPhypMixing13Bar,  input_mixing13_used);
    ph.SetCurrValue(kPhypMixing23Bar,  input_mixing23_used);
    ph.SetCurrValue(kPhypCPVBar    ,   input_dcp_used);
    ph.SetCurrValue(kPhypDm12sqBar,    input_dm12sq_used);
    ph.SetCurrValue(kPhypDm23sqBar,    input_dm23sq_used);
  }//kHypOsc3FBar
  else if(ph.Hypothesis() == kHypOsc3FNormBar) {
    ph.SetCurrValue(kPhypMixing12,  input_mixing12_used);
    ph.SetCurrValue(kPhypMixing13,  input_mixing13_used);
    ph.SetCurrValue(kPhypMixing23,  input_mixing23_used);
    ph.SetCurrValue(kPhypCPV    ,   input_dcp_used);
    ph.SetCurrValue(kPhypDm12sq,    input_dm12sq_used);
    ph.SetCurrValue(kPhypDm23sq,    input_dm23sq_used);
       
    ph.SetCurrValue(kPhypBetaBar,    input_beta_bar_used);
  }//kHypOsc3FNormBar
  else if(ph.Hypothesis() == kHypOsc3FBarNormBar) {
    ph.SetCurrValue(kPhypMixing12,  input_mixing12_used);
    ph.SetCurrValue(kPhypMixing13,  input_mixing13_used);
    ph.SetCurrValue(kPhypMixing23,  input_mixing23_used);
    ph.SetCurrValue(kPhypCPV    ,   input_dcp_used);
    ph.SetCurrValue(kPhypDm12sq,    input_dm12sq_used);
    ph.SetCurrValue(kPhypDm23sq,    input_dm23sq_used);
       
    ph.SetCurrValue(kPhypMixing12Bar,  input_mixing12_used);
    ph.SetCurrValue(kPhypMixing13Bar,  input_mixing13_used);
    ph.SetCurrValue(kPhypMixing23Bar,  input_mixing23_used);
    ph.SetCurrValue(kPhypCPVBar    ,   input_dcp_used);
    ph.SetCurrValue(kPhypDm12sqBar,    input_dm12sq_used);
    ph.SetCurrValue(kPhypDm23sqBar,    input_dm23sq_used);
       
    ph.SetCurrValue(kPhypBetaBar,    input_beta_bar_used);
  }//kHypOsc3FBarNormBar
  else {
    cout << "SetOscParams() " << util::AsString(ph.Hypothesis()) << " not implemented" << endl;
    exit(1);
  }
  
  if(print)
    cout << "- Physics hypothesis set: " << ph << endl;

  //save the values used
  if(iMH == 0) {
    if(hToy_th13_dcp_nh)
      hToy_th13_dcp_nh->Fill(input_mixing13_used, input_dcp_used);
    if(hToy_23_nh)
      hToy_23_nh->Fill(input_mixing23_used, input_dm23sq_used);
  }//NH
  else if(iMH == 1) {
    if(hToy_th13_dcp_ih)
      hToy_th13_dcp_ih->Fill(input_mixing13_used, input_dcp_used);
    if(hToy_23_ih)
      hToy_23_ih->Fill(input_mixing23_used, input_dm23sq_used);
  }//IH

}
//____________________________________________________________________________
void valor::util::OscParams::UpdateFitParamList(FitParamList& fparm, bool set_init, bool print)
{
  double input_mixing12_used, input_mixing13_used, input_mixing23_used;
  double input_dcp_used,      input_dm12sq_used,   input_dm23sq_used;
  double input_beta_bar_used;

  //Get the values to use (based on MH)
  if(iMH == 0) {
    //NH
    input_mixing12_used = mixing12_nh;
    input_mixing13_used = mixing13_nh;
    input_mixing23_used = mixing23_nh;
    input_dcp_used      = dcp_nh;
    input_dm12sq_used   = dm12sq_nh;
    input_dm23sq_used   = dm23sq_nh;
  }
  else if(iMH == 1) {
    //IH
    input_mixing12_used = mixing12_ih;
    input_mixing13_used = mixing13_ih;
    input_mixing23_used = mixing23_ih;
    input_dcp_used      = dcp_ih;
    input_dm12sq_used   = dm12sq_ih;
    input_dm23sq_used   = dm23sq_ih;
  }
  else {
    cerr << "SetOscParams has unknown value of iMH: " << iMH << endl;
    exit(1);
  }
  //Other
  input_beta_bar_used = beta_bar;
  
  //Update the FitParamList
  if(fparm.PhysicsHypothesis() == kHypNuMuDisapp2F) {
    if(set_init) {
      fparm.SetInitValue(kPhypMixing23,  input_mixing23_used);
      fparm.SetInitValue(kPhypDm23sq,    input_dm23sq_used);
    }
    else {
      fparm.SetCurrentValue(kPhypMixing23,  input_mixing23_used);
      fparm.SetCurrentValue(kPhypDm23sq,    input_dm23sq_used);
    }
  }//kHypNuMuDisapp2F
  else if(fparm.PhysicsHypothesis() == kHypOsc3F) {  
    if(set_init) {
      fparm.SetInitValue(kPhypMixing12,  input_mixing12_used);
      fparm.SetInitValue(kPhypMixing13,  input_mixing13_used);
      fparm.SetInitValue(kPhypMixing23,  input_mixing23_used);
      fparm.SetInitValue(kPhypCPV    ,   input_dcp_used);
      fparm.SetInitValue(kPhypDm12sq,    input_dm12sq_used);
      fparm.SetInitValue(kPhypDm23sq,    input_dm23sq_used);
    }
    else {
      fparm.SetCurrentValue(kPhypMixing12,  input_mixing12_used);
      fparm.SetCurrentValue(kPhypMixing13,  input_mixing13_used);
      fparm.SetCurrentValue(kPhypMixing23,  input_mixing23_used);
      fparm.SetCurrentValue(kPhypCPV    ,   input_dcp_used);
      fparm.SetCurrentValue(kPhypDm12sq,    input_dm12sq_used);
      fparm.SetCurrentValue(kPhypDm23sq,    input_dm23sq_used);
    }
  }//kHypOsc3F
  else if(fparm.PhysicsHypothesis() == kHypOsc3FBar) {
    if(set_init) {
      fparm.SetInitValue(kPhypMixing12,  input_mixing12_used);
      fparm.SetInitValue(kPhypMixing13,  input_mixing13_used);
      fparm.SetInitValue(kPhypMixing23,  input_mixing23_used);
      fparm.SetInitValue(kPhypCPV    ,   input_dcp_used);
      fparm.SetInitValue(kPhypDm12sq,    input_dm12sq_used);
      fparm.SetInitValue(kPhypDm23sq,    input_dm23sq_used);
       
      fparm.SetInitValue(kPhypMixing12Bar,  input_mixing12_used);
      fparm.SetInitValue(kPhypMixing13Bar,  input_mixing13_used);
      fparm.SetInitValue(kPhypMixing23Bar,  input_mixing23_used);
      fparm.SetInitValue(kPhypCPVBar    ,   input_dcp_used);
      fparm.SetInitValue(kPhypDm12sqBar,    input_dm12sq_used);
      fparm.SetInitValue(kPhypDm23sqBar,    input_dm23sq_used);
    }
    else {
      fparm.SetCurrentValue(kPhypMixing12,  input_mixing12_used);
      fparm.SetCurrentValue(kPhypMixing13,  input_mixing13_used);
      fparm.SetCurrentValue(kPhypMixing23,  input_mixing23_used);
      fparm.SetCurrentValue(kPhypCPV    ,   input_dcp_used);
      fparm.SetCurrentValue(kPhypDm12sq,    input_dm12sq_used);
      fparm.SetCurrentValue(kPhypDm23sq,    input_dm23sq_used);
       
      fparm.SetCurrentValue(kPhypMixing12Bar,  input_mixing12_used);
      fparm.SetCurrentValue(kPhypMixing13Bar,  input_mixing13_used);
      fparm.SetCurrentValue(kPhypMixing23Bar,  input_mixing23_used);
      fparm.SetCurrentValue(kPhypCPVBar    ,   input_dcp_used);
      fparm.SetCurrentValue(kPhypDm12sqBar,    input_dm12sq_used);
      fparm.SetCurrentValue(kPhypDm23sqBar,    input_dm23sq_used);
    }
  }//kHypOsc3FBar
  else if(fparm.PhysicsHypothesis() == kHypOsc3FNormBar) {
    if(set_init) {
      fparm.SetInitValue(kPhypMixing12,  input_mixing12_used);
      fparm.SetInitValue(kPhypMixing13,  input_mixing13_used);
      fparm.SetInitValue(kPhypMixing23,  input_mixing23_used);
      fparm.SetInitValue(kPhypCPV    ,   input_dcp_used);
      fparm.SetInitValue(kPhypDm12sq,    input_dm12sq_used);
      fparm.SetInitValue(kPhypDm23sq,    input_dm23sq_used);
       
      fparm.SetInitValue(kPhypBetaBar,    input_beta_bar_used);
    }
    else {
      fparm.SetCurrentValue(kPhypMixing12,  input_mixing12_used);
      fparm.SetCurrentValue(kPhypMixing13,  input_mixing13_used);
      fparm.SetCurrentValue(kPhypMixing23,  input_mixing23_used);
      fparm.SetCurrentValue(kPhypCPV    ,   input_dcp_used);
      fparm.SetCurrentValue(kPhypDm12sq,    input_dm12sq_used);
      fparm.SetCurrentValue(kPhypDm23sq,    input_dm23sq_used);
       
      fparm.SetCurrentValue(kPhypBetaBar,    input_beta_bar_used);
    }
  }//kHypOsc3FNormBar
  else if(fparm.PhysicsHypothesis() == kHypOsc3FBarNormBar) {
    if(set_init) {
      fparm.SetInitValue(kPhypMixing12,  input_mixing12_used);
      fparm.SetInitValue(kPhypMixing13,  input_mixing13_used);
      fparm.SetInitValue(kPhypMixing23,  input_mixing23_used);
      fparm.SetInitValue(kPhypCPV    ,   input_dcp_used);
      fparm.SetInitValue(kPhypDm12sq,    input_dm12sq_used);
      fparm.SetInitValue(kPhypDm23sq,    input_dm23sq_used);
       
      fparm.SetInitValue(kPhypMixing12Bar,  input_mixing12_used);
      fparm.SetInitValue(kPhypMixing13Bar,  input_mixing13_used);
      fparm.SetInitValue(kPhypMixing23Bar,  input_mixing23_used);
      fparm.SetInitValue(kPhypCPVBar    ,   input_dcp_used);
      fparm.SetInitValue(kPhypDm12sqBar,    input_dm12sq_used);
      fparm.SetInitValue(kPhypDm23sqBar,    input_dm23sq_used);
       
      fparm.SetInitValue(kPhypBetaBar,    input_beta_bar_used);
    }
    else {
      fparm.SetCurrentValue(kPhypMixing12,  input_mixing12_used);
      fparm.SetCurrentValue(kPhypMixing13,  input_mixing13_used);
      fparm.SetCurrentValue(kPhypMixing23,  input_mixing23_used);
      fparm.SetCurrentValue(kPhypCPV    ,   input_dcp_used);
      fparm.SetCurrentValue(kPhypDm12sq,    input_dm12sq_used);
      fparm.SetCurrentValue(kPhypDm23sq,    input_dm23sq_used);
       
      fparm.SetCurrentValue(kPhypMixing12Bar,  input_mixing12_used);
      fparm.SetCurrentValue(kPhypMixing13Bar,  input_mixing13_used);
      fparm.SetCurrentValue(kPhypMixing23Bar,  input_mixing23_used);
      fparm.SetCurrentValue(kPhypCPVBar    ,   input_dcp_used);
      fparm.SetCurrentValue(kPhypDm12sqBar,    input_dm12sq_used);
      fparm.SetCurrentValue(kPhypDm23sqBar,    input_dm23sq_used);
       
      fparm.SetCurrentValue(kPhypBetaBar,    input_beta_bar_used);
    }
  }//kHypOsc3FBarNormBar
  else {
    cout << "SetOscParams() " << util::AsString(fparm.PhysicsHypothesis()) << " not implemented" << endl;
    exit(1);
  }
  
  if(print)
    cout << "- Oscillation parameters in physics param list updated: " << fparm << endl;

  //save the values used
  if(iMH == 0) {
    if(hToy_th13_dcp_nh)
      hToy_th13_dcp_nh->Fill(input_mixing13_used, input_dcp_used);
    if(hToy_23_nh)
      hToy_23_nh->Fill(input_mixing23_used, input_dm23sq_used);
  }//NH
  else if(iMH == 1) {
    if(hToy_th13_dcp_ih)
      hToy_th13_dcp_ih->Fill(input_mixing13_used, input_dcp_used);
    if(hToy_23_ih)
      hToy_23_ih->Fill(input_mixing23_used, input_dm23sq_used);
  }//IH

}
//____________________________________________________________________________
void valor::util::OscParams::Print(ostream & stream) const
{
  stream << "   - True oscillation parameters: " << endl
	 << "     - Normal ordering: " << endl
	 << "       > " << util::AsString(kPhypMixing12) << " = "; ((randomize_mixing12_nh) ? cout << "* randomized *" : cout << mixing12_nh) << endl
	 << "       > " << util::AsString(kPhypMixing13) << " = "; ((randomize_mixing13_nh) ? cout << "* randomized *" : cout << mixing13_nh) << endl
	 << "       > " << util::AsString(kPhypMixing23) << " = "; ((randomize_mixing23_nh) ? cout << "* randomized *" : cout << mixing23_nh) << endl
	 << "       > " << util::AsString(kPhypCPV    )  << " = "; ((randomize_dcp_nh)      ? cout << "* randomized *" : cout << dcp_nh) << endl
	 << "       > " << util::AsString(kPhypDm12sq)   << " = "; ((randomize_dm12sq_nh)   ? cout << "* randomized *" : cout << dm12sq_nh )  << " eV^2/c^4" << endl
	 << "       > " << util::AsString(kPhypDm23sq)   << " = "; ((randomize_dm23sq_nh)   ? cout << "* randomized *" : cout << dm23sq_nh )  << " eV^2/c^4" << endl
	 << "     - Inverted ordering: " << endl
	 << "       > " << util::AsString(kPhypMixing12) << " = "; ((randomize_mixing12_ih) ? cout << "* randomized *" : cout << mixing12_ih) << endl
	 << "       > " << util::AsString(kPhypMixing13) << " = "; ((randomize_mixing13_ih) ? cout << "* randomized *" : cout << mixing13_ih) << endl
	 << "       > " << util::AsString(kPhypMixing23) << " = "; ((randomize_mixing23_ih) ? cout << "* randomized *" : cout << mixing23_ih) << endl
	 << "       > " << util::AsString(kPhypCPV    )  << " = "; ((randomize_dcp_ih)      ? cout << "* randomized *" : cout << dcp_ih) << endl
	 << "       > " << util::AsString(kPhypDm12sq)   << " = "; ((randomize_dm12sq_ih)   ? cout << "* randomized *" : cout << dm12sq_ih )  << " eV^2/c^4" << endl
	 << "       > " << util::AsString(kPhypDm23sq)   << " = "; ((randomize_dm23sq_ih)   ? cout << "* randomized *" : cout << dm23sq_ih )  << " eV^2/c^4" << endl
	 << "     - Other: " << endl
	 << "       > " << util::AsString(kPhypBetaBar)  << " = "; ((randomize_beta_bar)    ? cout << "* randomized *" : cout << beta_bar) << endl
	 << "       > Mass hierarchy = "; ((marginalize_mh) ? cout << "Margninalized (NH prior = " << mh_nh_prior << ")" : cout << util::AsString(MH)) << endl;
  if(!priors_file.empty()) {
    if((randomize_mixing13_nh && randomize_dcp_nh)
       || (randomize_mixing13_ih && randomize_dcp_ih)
       || (randomize_mixing23_nh && randomize_dm23sq_nh)
       || (randomize_mixing23_ih && randomize_dm23sq_ih)
       ) {
      stream << "       > Using priors from file: " << priors_file << " for parameters:" << endl;
      if(randomize_mixing13_nh && randomize_dcp_nh)
	stream << "          - " << util::AsString(kPhypMixing13) << " vs " << util::AsString(kPhypCPV) << " (NH)" << endl;
      if(randomize_mixing13_ih && randomize_dcp_ih)
	stream << "          - " << util::AsString(kPhypMixing13) << " vs " << util::AsString(kPhypCPV) << " (IH)" << endl;
      if(randomize_mixing23_nh && randomize_dm23sq_nh)
	stream << "          - " << util::AsString(kPhypMixing23) << " vs " << util::AsString(kPhypDm23sq) << " (NH)" << endl;
      if(randomize_mixing23_ih && randomize_dm23sq_ih)
	stream << "          - " << util::AsString(kPhypMixing23) << " vs " << util::AsString(kPhypDm23sq) << " (IH)" << endl;
    }
  }
  stream << endl;
  
}
//____________________________________________________________________________
void valor::util::OscParams::Write(TFile * f)
{
  f->cd();

  //priors
  if(hPrior_th13_dcp_nh) {
    hPrior_th13_dcp_nh->SetDirectory(0);
    hPrior_th13_dcp_nh->Write();
  }
  if(hPrior_th13_dcp_ih) {
    hPrior_th13_dcp_ih->SetDirectory(0);
    hPrior_th13_dcp_ih->Write();
  }
  if(hPrior_23_nh) {
    hPrior_23_nh->SetDirectory(0);
    hPrior_23_nh->Write();
  }
  if(hPrior_23_ih) {
    hPrior_23_ih->SetDirectory(0);
    hPrior_23_ih->Write();
  }
  //toys used
  if(hToy_th13_dcp_nh) {
    hToy_th13_dcp_nh->SetDirectory(0);
    hToy_th13_dcp_nh->Write();
  }
  if(hToy_th13_dcp_ih) {
    hToy_th13_dcp_ih->SetDirectory(0);
    hToy_th13_dcp_ih->Write();
  }
  if(hToy_23_nh) {
    hToy_23_nh->SetDirectory(0);
    hToy_23_nh->Write();
  }
  if(hToy_23_ih) {
    hToy_23_ih->SetDirectory(0);
    hToy_23_ih->Write();
  }
  if(hToy_MH) {
    hToy_MH->SetDirectory(0);
    hToy_MH->Write();
  }
}
//____________________________________________________________________________
//____________________________________________________________________________
//____________________________________________________________________________
//____________________________________________________________________________
//____________________________________________________________________________
//____________________________________________________________________________
//____________________________________________________________________________
//____________________________________________________________________________
