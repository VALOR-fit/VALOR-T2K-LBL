//________________________________________________________________________________________
/*

*/
//________________________________________________________________________________________

#include <cassert>

#include "year1_numu_osc_analysis/SuperkExpectation.h"
#include "year1_numu_osc_analysis/Conventions.h"

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
double SuperkExpectation::Generate(double Er, const PhysParamList & fp)
{
// Calculate dN^{exp}_{SuperK}(E_{reco}) / dE_{reco} for the input E_{reco} and 
// the input physics parameter list.
//

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

  const TAxis & ereco_binning = fSuperkInputs->ErecoBinning();
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
         double prob_surv         = 1;///here
         double eff_reco_m_sk     = fSuperkInputs -> ReconstructionEfficiency(Etrue,m);
         double eff_reco_cc_nd280 = fND280Inputs  -> ReconstructionEfficiency(Etrue); 

         double Ct = fn_ratio * 
                     (xsec_m_sk / xsec_totcc_nd280) * 
                     prob_surv * 
                     (eff_reco_m_sk / eff_reco_cc_nd280);

         double IntgEreco = 0; // integral over E_{reco}

         // E_{reco} bin loop
         for(int irbin = 1; irbin <= ereco_binning.GetNbins(); irbin++) {
            double Ereco = ereco_binning.GetBinCenter(irbin);

            // calc stuff with E_{reco} energy dependency only
            double nobscclike_nd280 = fND280Inputs -> NumOfObsCCLikeEvents (Ereco);
            double nexpnc_nd280     = fND280Inputs -> NumOfEstNCBkgEvents  (Ereco);
            double ncc_nd280        = TMath::Max(0., nobscclike_nd280-nexpnc_nd280);
            double eff_sel_s_m_sk   = fSuperkInputs -> SelectionEfficiency(Etrue,s,m);
            double eff_sel_cc_nd280 = fND280Inputs  -> SelectionEfficiency(Etrue); 
            double Cr = ncc_nd280 * 
                        (eff_sel_s_m_sk / eff_sel_cc_nd280);

            // calc stuff with both E_{true}, E_{reco} energy dependency
            double R_s_m_sk    = fSuperkInputs -> EnergyResolutionMatrix(Etrue,Ereco,s,m);
            double RT_cc_nd280 = fND280Inputs  -> EnergyResolutionMatrix(Etrue,Ereco);
            double Ctr = R_s_m_sk * RT_cc_nd280;

            IntgEreco += (Cr*Ctr);
         }//E_{reco}

         IntgEtrue += (Ct*IntgEreco);
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
  fErBinning       = 0;
  fEtBinning       = 0;
  fObsCCLikeEr     = 0;
  fEstNCBkgEr      = 0;
  fRecoEffCC       = 0;
  fSelEffCC        = 0;
  fResolMtrxCC     = 0;
}
//________________________________________________________________________________________
void SuperkExpectation::ND280Inputs::CleanUp(void)
{
  if ( fErBinning   ) delete fErBinning;
  if ( fEtBinning   ) delete fEtBinning;
  if ( fObsCCLikeEr ) delete fObsCCLikeEr;
  if ( fEstNCBkgEr  ) delete fEstNCBkgEr;
  if ( fRecoEffCC   ) delete fRecoEffCC;
  if ( fSelEffCC    ) delete fSelEffCC;
  if ( fResolMtrxCC ) delete fResolMtrxCC;
}
//________________________________________________________________________________________
bool SuperkExpectation::ND280Inputs::Read(string config_file)
{
  return true;
}
//________________________________________________________________________________________
double SuperkExpectation::ND280Inputs::NumOfObsCCLikeEvents(double Ereco) const
{
// observed number of CC-like events

  int    bin = fObsCCLikeEr->FindBin(Ereco);
  double nev = fObsCCLikeEr->GetBinContent(bin);

  return nev;
}
//________________________________________________________________________________________
double SuperkExpectation::ND280Inputs::NumOfEstNCBkgEvents(double Ereco) const
{
// estimated NC background

  int    bin = fEstNCBkgEr->FindBin(Ereco);
  double nev = fEstNCBkgEr->GetBinContent(bin);

  return nev;
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
double SuperkExpectation::ND280Inputs::ReconstructionEfficiency(double Etrue) const
{
  return 0; 
}
//________________________________________________________________________________________
double SuperkExpectation::ND280Inputs::SelectionEfficiency(double Ereco) const
{
  return 0; 
}
//________________________________________________________________________________________
double SuperkExpectation::ND280Inputs::EnergyResolutionMatrix(
  double Etrue, double Ereco) const
{
  return 0; 
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

}
//________________________________________________________________________________________
void SuperkExpectation::BeamInputs::CleanUp(void)
{

}
//________________________________________________________________________________________
bool SuperkExpectation::BeamInputs::Read(string config_file)
{
  return true;
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
  return 0;
}
//________________________________________________________________________________________
double SuperkExpectation::BeamInputs::ND280Flux(double Etrue) const
{
  return 0;
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


