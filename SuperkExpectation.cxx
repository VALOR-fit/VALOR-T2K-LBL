//________________________________________________________________________________________
/*

*/
//________________________________________________________________________________________

#include <TFile.h>
#include <TXMLEngine.h>
#include <TString.h>

#include <cassert>
#include <cstdlib>
#include <iostream>

#include "year1_numu_osc_analysis/SuperkExpectation.h"
#include "year1_numu_osc_analysis/Conventions.h"
#include "year1_numu_osc_analysis/Utilities.h"

using std::cout;
using std::cerr;
using std::endl;

using namespace t2kana;
using namespace t2kana::yr1numu;
using namespace t2kana::yr1numu::conv;

ClassImp(SuperkExpectation)

//________________________________________________________________________________________
SuperkExpectation::SuperkExpectation(void)
{
  this->Init();
}
//________________________________________________________________________________________
SuperkExpectation::~SuperkExpectation(void)
{
  this->CleanUp();
}
//________________________________________________________________________________________
double SuperkExpectation::Generate(double Ereco, const PhysParamList & fp)
{
// Calculate dN^{exp}_{SuperK}(E_{reco}) / dE_{reco} for the input E_{reco} and 
// the input physics parameter list.
//

  double expectation     = 0.;

  // input physics params (params included in the fit / systematics params)
  double dmsq23          = fp.CurrValue(kPhypDm23sq);
  double sinsq2theta23   = fp.CurrValue(kPhypSinsq2Theta23);

  // calc stuff with no energy dependency
  double Mfiducial_nd280 = fND280Inputs  -> FiducialMass();
  double Npot_nd280      = fND280Inputs  -> AnalyzedSamplePOT();
  double A_nd280         = 12;///////
  double Mfiducial_sk    = fSuperkInputs -> FiducialMass();
  double Npot_sk         = fSuperkInputs -> AnalyzedSamplePOT();
  double A_sk            = 16;///////
  double C0 = (Npot_sk    * Mfiducial_sk    * A_nd280) / 
              (Npot_nd280 * Mfiducial_nd280 * A_sk);

  SuperkSample_t sk_samples[kNSuperkSamples] =  { kSKsFCFV1RingMuLike, kSKsFCFVmRingMuLike };
  SuperkMode_t sk_modes[kNSuperkModes] = { kSKmNuMuCCQE, kSKmNuMuNonCCQE, kSKmNC };

  const TAxis & etrue_binning = fSuperkInputs->EtrueBinning();

  // loop over SK samples & modes/sample
  for(int is = 0; is < kNSuperkSamples; is++) {
    SuperkSample_t s = sk_samples[is];
    for(int im = 0; im < kNSuperkModes; im++) {
       SuperkMode_t m = sk_modes[im];

       double IntgEtrue = 0; // integral over E_{true}

       // E_{true} bin loop
       for(int itbin = 1; itbin <= etrue_binning.GetNbins(); itbin++) {
         double Etrue = etrue_binning.GetBinCenter(itbin);

         // calc stuff with E_{true} energy dependency only
         double fn_ratio          = fBeamInputs->FNRatio(Etrue);
         double xsec_totcc_nd280  = fGeneratorInputs -> SigmaTotCCND280(Etrue);
         double xsec_tot_sk       = fGeneratorInputs -> SigmaTotWater(Etrue);
         double evtfrac_m_sk      = fGeneratorInputs -> ModeFraction(Etrue,m);
         double xsec_m_sk         = xsec_tot_sk * evtfrac_m_sk;
         double prob_surv         = util::NuMuSurvProb2D(sinsq2theta23, dmsq23, Etrue, m);
         double R_s_m_sk          = fSuperkInputs -> EnergyResolutionMatrix(Etrue,Ereco,s,m);
         double eff_sel_s_m_sk    = fSuperkInputs -> SelectionEfficiency(Ereco,s,m);
         double eff_reco_m_sk     = fSuperkInputs -> ReconstructionEfficiency(Etrue,m);
         double nobscclike_nd280  = fND280Inputs -> NumOfObsCCLikeEvents ();
         double nexpnc_nd280      = fND280Inputs -> NumOfEstNCBkgEvents  ();
         double eff_nd280         = fND280Inputs  -> Efficiency(Etrue); 

         double Ct = fn_ratio * 
                     (xsec_m_sk / xsec_totcc_nd280) * 
                     prob_surv * 
                     (nobscclike_nd280 - nexpnc_nd280) *
                     R_s_m_sk * (eff_sel_s_m_sk * eff_reco_m_sk / eff_nd280);     

         IntgEtrue += Ct;
       }//E_{true}

       expectation += (C0*IntgEtrue);
    }//m
  }//s

  return expectation;
}
//________________________________________________________________________________________
void SuperkExpectation::Init(void)
{
  fND280Inputs     = new ND280Inputs;
  fSuperkInputs    = new SuperkInputs;
  fBeamInputs      = new BeamInputs; 
  fGeneratorInputs = new GeneratorInputs; 
}
//________________________________________________________________________________________
void SuperkExpectation::CleanUp(void)
{
  delete fND280Inputs;
  delete fSuperkInputs;
  delete fBeamInputs; 
  delete fGeneratorInputs;
}
//________________________________________________________________________________________
//........................................................................................
// nested class ND280Inputs
//........................................................................................
//________________________________________________________________________________________
SuperkExpectation::ND280Inputs::ND280Inputs()
{
  this->Init();
}
//________________________________________________________________________________________
SuperkExpectation::ND280Inputs::~ND280Inputs()
{
  this->CleanUp();
}
//________________________________________________________________________________________
void SuperkExpectation::ND280Inputs::Init(void)
{
  fPOT             = 0;
  fFidMass         = 0;
  fEtBinning       = 0;
  fEfficiency      = 0;
}
//________________________________________________________________________________________
void SuperkExpectation::ND280Inputs::CleanUp(void)
{
  if ( fEtBinning   ) delete fEtBinning;
  if ( fEfficiency  ) delete fEfficiency;
}
//________________________________________________________________________________________
bool SuperkExpectation::ND280Inputs::Read(string config_file)
{
  return true;
}
//________________________________________________________________________________________
double SuperkExpectation::ND280Inputs::NumOfObsCCLikeEvents(void) const
{
// observed number of CC-like events

  return fNumObsCCLike;
}
//________________________________________________________________________________________
double SuperkExpectation::ND280Inputs::NumOfEstNCBkgEvents(void) const
{
// estimated NC background

  double Nnc0 = fNumEstNCBkg; // nominal estimation 

  //
  // ...
  //

  return Nnc0;
}
//________________________________________________________________________________________
double SuperkExpectation::ND280Inputs::FiducialMass(void) const
{
  return fFidMass;
}
//________________________________________________________________________________________
double SuperkExpectation::ND280Inputs::AnalyzedSamplePOT(void) const
{
  return fPOT;
}
//________________________________________________________________________________________
double SuperkExpectation::ND280Inputs::Efficiency(double Etrue) const
{
  return 1; ///TODO
}
//________________________________________________________________________________________
//........................................................................................
// nested class SuperkInputs
//........................................................................................
//________________________________________________________________________________________
SuperkExpectation::SuperkInputs::SuperkInputs()
{
  this->Init();
}
//________________________________________________________________________________________
SuperkExpectation::SuperkInputs::~SuperkInputs()
{
  this->CleanUp();
}
//________________________________________________________________________________________
void SuperkExpectation::SuperkInputs::Init(void)
{
  fPOT       = 0;  
  fFidMass   = 0; 
  fErBinning = 0;
  fEtBinning = 0;

  for(int im=0; im<knm; im++) {
    fRecoEff[im] = 0;
  }
  for(int is=0; is<kns; is++) {
    for(int im=0; im<knm; im++) {
      fSelEff   [is][im] = 0; 
      fResolMtrx[is][im] = 0; 
    }
  }
}
//________________________________________________________________________________________
void SuperkExpectation::SuperkInputs::CleanUp(void)
{
  if (fErBinning ) delete fErBinning;
  if (fEtBinning ) delete fEtBinning;

  for(int im=0; im<knm; im++) {
    if(fRecoEff[im]) delete fRecoEff[im];
  }
  for(int is=0; is<kns; is++) {
    for(int im=0; im<knm; im++) {
      if ( fSelEff   [is][im] ) delete  fSelEff   [is][im];
      if ( fResolMtrx[is][im] ) delete  fResolMtrx[is][im];
    }
  }
}
//________________________________________________________________________________________
bool SuperkExpectation::SuperkInputs::Read(string config_file)
{
  return true;
}
//________________________________________________________________________________________
double SuperkExpectation::SuperkInputs::FiducialMass(void) const
{
  return 0;
}
//________________________________________________________________________________________
double SuperkExpectation::SuperkInputs::AnalyzedSamplePOT(void) const
{
  return 0;
}
//________________________________________________________________________________________
double SuperkExpectation::SuperkInputs::ReconstructionEfficiency(
    double Etrue, SuperkMode_t m) const
{
  return 0;
}
//________________________________________________________________________________________
double SuperkExpectation::SuperkInputs::SelectionEfficiency(
    double Ereco, SuperkSample_t s, SuperkMode_t m) const
{
  return 0;
}
//________________________________________________________________________________________
double SuperkExpectation::SuperkInputs::EnergyResolutionMatrix(
    double Etrue, double Ereco, SuperkSample_t s, SuperkMode_t m) const
{
  return 0;
}
//________________________________________________________________________________________
//........................................................................................
// nested class BeamInputs
//........................................................................................
//________________________________________________________________________________________
SuperkExpectation::BeamInputs::BeamInputs()
{
  this->Init();
}
//________________________________________________________________________________________
SuperkExpectation::BeamInputs::~BeamInputs()
{
  this->CleanUp();
}
//________________________________________________________________________________________
void SuperkExpectation::BeamInputs::Init(void)
{
  fSuperkFlux = 0; 
  fND280Flux  = 0; 
}
//________________________________________________________________________________________
void SuperkExpectation::BeamInputs::CleanUp(void)
{
  if (fSuperkFlux) delete fSuperkFlux;
  if (fND280Flux ) delete fND280Flux;
}
//________________________________________________________________________________________
bool SuperkExpectation::BeamInputs::Read(string config_file)
{
// Read beam analysis inputs from the input XML config file
//
  cout << "SKexp/Beam >> " << "Reading config file: " << config_file << endl;

  // Parse the input configuration file
  TXMLEngine* xml = new TXMLEngine;   
  XMLDocPointer_t xmldoc = xml->ParseFile(config_file.c_str());
  if (xmldoc==0) {
    cerr << "SKexp/Beam >> " << "Can not parse input: " << config_file << endl;
    delete xml;
    return false;  
  }
  XMLNodePointer_t root_node = xml->DocGetRootElement(xmldoc);

  // Get 
  fSuperkFlux = util::XML2TH1D(xml,root_node,"beam/superk_flux");
  if(!fSuperkFlux) {
    cerr << "SKexp/Beam >> No input SuperK flux histogram!" << endl;
    exit(1);
  }

  fND280Flux  = util::XML2TH1D(xml,root_node,"beam/nd280_flux");
  if(!fND280Flux) {
    cerr << "SKexp/Beam >> No input nd280 flux histogram!" << endl;
    exit(1);
  }

  // Release memory
  xml->FreeDoc(xmldoc);
  delete xml;

  cout << "SKexp/Beam >> Done reading beam inputs!" << endl;

  return true;
}
//________________________________________________________________________________________
void SuperkExpectation::BeamInputs::Draw(TVirtualPad * c, Option_t * option) const
{
  if(!c) return;

  c->cd();

  TString opt = option;
  util::SetDefaultStyle();

  if(opt.Contains("Nd280Flux")) { 
  } 
  else
  if(opt.Contains("SKFlux")) { 
    fSuperkFlux->Draw();
  } 

  c->Update();
}
//________________________________________________________________________________________
double SuperkExpectation::BeamInputs::FNRatio(double Etrue) const 
{ 
  double f = this->SuperkFlux(Etrue);
  double n = this->ND280Flux(Etrue);

  double R = (n>0 && f>=0) ? (f/n) : 0;
  return R;
}
//________________________________________________________________________________________
double SuperkExpectation::BeamInputs::SuperkFlux (double Etrue) const
{
  return 1;//TODO
}
//________________________________________________________________________________________
double SuperkExpectation::BeamInputs::ND280Flux(double Etrue) const
{
  return 1;//TODO
}
//________________________________________________________________________________________
//........................................................................................
// nested class GeneratorInputs
//........................................................................................
//________________________________________________________________________________________
SuperkExpectation::GeneratorInputs::GeneratorInputs()
{
  this->Init();
}
//________________________________________________________________________________________
SuperkExpectation::GeneratorInputs::~GeneratorInputs()
{
  this->CleanUp();
}
//________________________________________________________________________________________
void SuperkExpectation::GeneratorInputs::Init(void)
{
  for(int im=0; im<knm; im++) {
    fModeFraction[im] = 0;
  }
  fSigmaTotWater   = 0;
  fSigmaTotCCND280 = 0;
}
//________________________________________________________________________________________
void SuperkExpectation::GeneratorInputs::CleanUp(void)
{
  for(int im=0; im<knm; im++) {
    if ( fModeFraction[im] ) delete fModeFraction[im];
  }
  if ( fSigmaTotWater   ) delete fSigmaTotWater;
  if ( fSigmaTotCCND280 ) delete fSigmaTotCCND280;
}
//________________________________________________________________________________________
bool SuperkExpectation::GeneratorInputs::Read(string config_file)
{
  return true;  
}
//________________________________________________________________________________________
double SuperkExpectation::GeneratorInputs::ModeFraction(double Etrue, SuperkMode_t m) const
{
  return 0;
}
//________________________________________________________________________________________
double SuperkExpectation::GeneratorInputs::SigmaTotWater(double Etrue) const
{
  return 0;
}
//________________________________________________________________________________________
double SuperkExpectation::GeneratorInputs::SigmaTotCCND280(double Etrue) const
{
  return 0;
}
//________________________________________________________________________________________


