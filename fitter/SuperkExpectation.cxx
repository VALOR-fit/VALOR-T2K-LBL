//________________________________________________________________________________________
/*
 Changes since v1.0.0:

 @ 16/03/2011 - CA
   Old SuperkInputs, ND280Inputs, BeamInputs and GeneratorInputs nested classes were moved 
   outside SuperkExpectation.h/cxx files and were renamed to SuperkMC, ND280Measurement,
   NuFlux and NuXSec.
 @ 28/03/2011 - CA
   Remove CurrPredictedErecoSpectrum() as it doesn't make sense using the combined spectrum 
   for all the different samples considered. The default was to return the spectrum for the 
   `1-R mu-like' sample but not typing this explicitly could lead to confusion. From now only
   CurrPredictedErecoSpectrum(s) and CurrPredictedErecoSpectrum(s,m) will be used.
 @ 06/06/2011 - CA
   Switched to using SuperkMC::PhysHypothesis() method for setting one of disappearance, 
   decay or decoherence models. Added the `CCQE ring-counting efficiency' systematic err in
   ResetSystematicVariations(), TweakSystematic() and TweakAllSystematics() methods.
   This is an Ereco-dependent error (3 Ereco bins are considered by SK) with Ereco bin-to-bin 
   correlations included/ The old Ereco-independent CCQE efficieny error was kept and renamed 
   to `CCQE other efficiency error' (includes all error sources except ring-counting). 
   There were corresponding modifications in SuperkMC.
 @ 07/06/2011 - CA
   Added systematic params "sk_ccqe_ring_counting_eff_0", "sk_ccqe_ring_counting_eff_1" and 
   "sk_ccqe_ring_counting_eff_2" in TweakSystematic(). They can be used for tweaking the
   CCQE ring-counting efficiency in single Ereco bins neglecting correlations. 
   For validation purposes only.
 @ 05/07/2011 - CA
   Adapted code to use the modified FitParamList and PhysHypothesis classes.
   Adapted code to changes in SuperkMC and NuFlux for including correlated systematics in
   the oscillation fit. Added SetNuisanceParams(const FitParamList &) method.
   There are now 2 interfaces to the SuperkExpectation code: One in which the pdfs can 
   be tweaked using the various TweakXYZ() methods and the extrapolated spectra can be 
   generated via Generate(const PhysHypothesis &). There is also another interface where 
   pdf tweaking and spectrum generation take place by invoking the single 
   Generate(const FitParamList &) method. The inpout FitParamList can contain parameters
   that describe the physics hypothesis and the systematic tweaks to the pdfs. The former
   interface is the one that was installed initially and used in validation. The latter 
   interface is meant to be used in the fits.
 @ 13/07/2011 - CA
   Take into account the statistical error on the measured number of ND280 CC-like events
   and propagate it to SuperK.
 @ 18/07/2011 - CA
   In the 3-flavour analysis, applying systematics and oscillations in the SuperK MC
   templates takes place in separate steps. The SuperkMC interface has changed slightly.
   Call SuperkMC::UpdateTemplates() in SuperkExpectation::Generate(const PhysHypothesis &)
   before updating the Ereco and Etrue spectra.
 @ 26/07/2011 - CA
   Took newly-added `nue signal' templates into account. They should be neglected under all 
   but the 3-flavour oscillation hypothesis. 
 @ 27/09/2011 - CA
   In SetNuisanceParams(const FitParamList & fp) don't call ResetSystematicVariations()
   as this messes things-up in the case where nuisance params are not included in the fit
   but are tweaked externally (via one of the TweakXYZ() methods).
@  06/10/11 - CA
   Added `PrintSystematicVariations()'.

*/
//________________________________________________________________________________________

//#define __debug_level_1__
//#define __debug_level_2__
//#define __debug_level_3__

#include <cassert>
#include <cstdlib>
#include <iostream>

#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TLegend.h>
#include <TSpline.h>
#include <TRandom3.h>
#include <TSystem.h>

#include "SuperkExpectation.h"
#include "ND280Measurement.h"
#include "SuperkMC.h"
#include "NuFlux.h"
#include "NuXSec.h"
#include "Conventions.h"
#include "Utilities.h"

using std::cout;
using std::cerr;
using std::endl;

using namespace t2kana;
using namespace t2kana::conv;

ClassImp(SuperkExpectation)

//________________________________________________________________________________________
SuperkExpectation::SuperkExpectation()
{
  this->Init();
}
//________________________________________________________________________________________
SuperkExpectation::~SuperkExpectation(void)
{
  this->CleanUp();
}
//________________________________________________________________________________________
void SuperkExpectation::Generate(const PhysHypothesis & ph)
{
// Calculate dN^{exp}_{SuperK}(E_{reco}) / dE_{reco} for the input E_{reco} and 
// the input physics hypothesis
//

  if(!fIsConfigured) {
    cerr << "Exp >> Not configured properly! " << endl;
    exit(1);
  }

  // set params driving the oscillation, decoherence or decay model
  fSuperkMC->Hypothesis().Copy(ph);

  // Ignore flux shape tweak and calculate the number of 1Rmu events.
  // Use that number of events to normalise the prediction obtain when the flux shape tweaking 
  // is taken into account.
  double shape_norm = 1.;
  if(this->GetNuFlux()->ShapeTweaked()) {
#ifdef __debug_level_2__
     cout << "Exp >> Flux shape was tweaked! Ignoring temporarily and generating SK spectra..." << endl;
#endif
     this->GetNuFlux()->IgnoreShapeTweaking(true); 
     this->Generate(ph); // generate untweaked predictions
     // get nominal Nsk
     shape_norm = fPredictErecoSpectrumTot[kSKsFCFV1RingMuLike]->Integral(); 
#ifdef __debug_level_2__
     cout << "Exp >> Number of 1Rmu events ignoring flux shape tweak (used for normalisation) = " 
          << shape_norm << endl;
#endif
     this->GetNuFlux()->IgnoreShapeTweaking(false); 
#ifdef __debug_level_2__
     cout << "Exp >> Flux shape tweaking was re-enabled..." << endl;
#endif
  }

  // correction factor from ND280
  double correction = fND280Measurement->CorrectionFactor();

  // update SK MC templates
  fSuperkMC->UpdateTemplates();

  // update Ereco spectrum prediction
  fSuperkMC->UpdateErecoSpectra();

  for(int is = 0; is < kns; is++) {
    SuperkSample_t s = kSKSamples[is];
    if(s == kSKsUndefined) continue;
    int nbins = 0; 
    if(s == kSKsFCFV1RingMuLike) nbins = kNErecoBins1Rmu;
    if(s == kSKsFCFV1RingELike)  nbins = kNErecoBins1Re;
    for(int ibin = 1; ibin <= nbins; ibin++) {
       double ns = correction * fSuperkMC->CurrErecoSpectrum(s)->GetBinContent(ibin);
       fPredictErecoSpectrumTot[s]->SetBinContent(ibin, ns);
       for(int im = 0; im < knm; im++) {
           SuperkMode_t m = kSKReactionModes[is][im];
           if(m == kSKmUndefined) continue;
           double nsm = correction * fSuperkMC->CurrErecoSpectrum(s,m)->GetBinContent(ibin);
           fPredictErecoSpectrum[s][m]->SetBinContent(ibin, nsm);       
       }//m	
    }//ibin
  }//s


  // include beam normalization err (expressed in terms of Nsk/Nnd280)
  if(this->GetNuFlux()->NormTweaked()) {
     double beam_norm = 1;
     if( this->GetNuFlux()->UsingKuboSansCode() ) {
        // input physics params (params included in the fit / systematics params)
        this -> GetSuperkMC() -> UpdateEtrueSpectra();
        this -> GetND280Measurement() -> UpdateEtrueSpectra();
        TH1D*  hnumu1Rmu     = this -> GetSuperkMC() -> Curr1RmuNuMuEtrue    ();
        TH1D*  hnumubar1Rmu  = this -> GetSuperkMC() -> Curr1RmuNuMuBarEtrue ();
        TH1D*  hnue1Rm       = this -> GetSuperkMC() -> Curr1RmuNuEEtrue     ();
        TH1D*  hnumuND       = this -> GetND280Measurement() -> CurrNuMuEtrue    ();
        TH1D*  hnumubarND    = this -> GetND280Measurement() -> CurrNuMuBarEtrue ();
        TH1D*  hnueND        = this -> GetND280Measurement() -> CurrNuEEtrue     ();

        beam_norm = this->GetNuFlux()->SuperkFluxNorm(
                   hnumu1Rmu, hnumubar1Rmu, hnue1Rm, hnumuND, hnumubarND, hnueND);
     } else {
        beam_norm = this->GetNuFlux()->SuperkFluxNorm(0,0,0,0,0,0);
     }

     for(int is = 0; is < kns; is++) {
        SuperkSample_t s = kSKSamples[is];
        if(s == kSKsUndefined) continue;
        fPredictErecoSpectrumTot[s] -> Scale (beam_norm);
        for(int im = 0; im < knm; im++) {
            SuperkMode_t m = kSKReactionModes[is][im];
            if(m == kSKmUndefined) continue;
            fPredictErecoSpectrum[s][m] -> Scale (beam_norm);
        }//m	
     }//s
  }//norm tweaked?

  if(this->GetNuFlux()->ShapeTweaked()) {
     int idx = util::FindSampleSlot(kSKsFCFV1RingMuLike);
     double nsktwk = fPredictErecoSpectrumTot[idx]->Integral();
     shape_norm /= nsktwk;
     for(int is = 0; is < kns; is++) {
        SuperkSample_t s = kSKSamples[is];
        if(s == kSKsUndefined) continue;
        fPredictErecoSpectrumTot[s] -> Scale (shape_norm);
        for(int im = 0; im < knm; im++) {
            SuperkMode_t m = kSKReactionModes[is][im];
            if(m == kSKmUndefined) continue;
            fPredictErecoSpectrum[s][m] -> Scale (shape_norm);
        }//m	
     }//s
  }
}
//________________________________________________________________________________________
void SuperkExpectation::Generate(const FitParamList & fp)
{
// SuperkExpectation interface method for fitters.
// Calculate dN^{exp}_{SuperK}(E_{reco}) / dE_{reco} for the input E_{reco} and 
// the input physics + nuisance parameter list.
//

  // set values of nuisance parameters 
  // (if any nuisance parameter is included in the fit)
  this->SetNuisanceParams(fp);

  // extract the values of params driving the oscillation, decoherence or decay model
  PhysHypothesis ph(fp);

  // call main Generate() method to generate SuperK predictions for the input
  // physics hypothesis using the (potentially) tweaked p.d.fs
  this->Generate(ph);
}
//________________________________________________________________________________________
TH1D* SuperkExpectation::CurrPredictedErecoSpectrum(t2kana::conv::SuperkSample_t s) const 
{
  int sidx = util::FindSampleSlot(s);
  if(sidx == -1) {
    return 0;
  }
  return fPredictErecoSpectrumTot[sidx];  
}
//________________________________________________________________________________________
TH1D* SuperkExpectation::CurrPredictedErecoSpectrum(
  t2kana::conv::SuperkSample_t s, t2kana::conv::SuperkMode_t m) const 
{ 
  int sidx = util::FindSampleSlot(s);
  if(sidx == -1) {
    return 0;
  }
  int midx = util::FindModeSlot(s,m);
  if(midx == -1) {
    return 0;
  }

  return fPredictErecoSpectrum[sidx][midx]; 
}
//________________________________________________________________________________________
void SuperkExpectation::ResetSystematicVariations(void)
{
  cout << "Exp >> Reseting all systematic variations." << endl;

  this->GetSuperkMC()->TweakEnergyScale                (0.);
  this->GetSuperkMC()->TweakOtherEfficiencyCCQE        (0.);
  this->GetSuperkMC()->TweakRingCountingEfficiencyCCQE (true,0.,0.,0.); 
  this->GetSuperkMC()->TweakEfficiencyCCnonQE          (0.);
  this->GetSuperkMC()->TweakEfficiencyNC               (0.);     
  this->GetSuperkMC()->TweakEfficiencyNuECC            (0.);
  this->GetSuperkMC()->TweakEfficiency1ReSig          (0.);
  this->GetSuperkMC()->TweakEfficiency1ReBkg          (0.);

  this->GetND280Measurement()->TweakEfficiency (0.);
  this->GetND280Measurement()->FluctuateNObsCCLike(0.);

  this->GetNuFlux()->TweakNorm (0.);
  this->GetNuFlux()->TweakShape(true,0);  

  this->GetNuXSec()->TweakCCQE        (0.);
  this->GetNuXSec()->TweakCC1piToCCQE (0.);
  this->GetNuXSec()->TweakCCCohToCCQE (0.);
  this->GetNuXSec()->TweakCCOthToCCQE (0.);
  this->GetNuXSec()->TweakNCToCCQE    (0.);
  this->GetNuXSec()->TweakNC1pi0ToCCQE(0.);
  this->GetNuXSec()->TweakNCCohToCCQE (0.);
  this->GetNuXSec()->TweakNCOthToCCQE (0.);
  this->GetNuXSec()->TweakNuEToNuMu   (0.);
  this->GetNuXSec()->TweakSKFSI       (0.);
}
//________________________________________________________________________________________
void SuperkExpectation::TweakSystematic(string systematic, double twkdial)
{
  //
  // SuperK detector systematics
  //
  if(systematic == "sk_energy_scale") {
     this->GetSuperkMC()->TweakEnergyScale(twkdial);
  }
  else
  if(systematic == "sk_ccqe_other_eff") {
     this->GetSuperkMC()->TweakOtherEfficiencyCCQE(twkdial);
  }
  else
  // tweak CCQE energy-dependent ring-counting efficiency generating random
  // efficiency tweaks and taking into account bin-to-bin correlations
  if(systematic == "sk_ccqe_ring_counting_eff") {
     this->GetSuperkMC()->TweakRingCountingEfficiencyCCQE(); //< Cholesky method
  }
  //... or tweak the CCQE ring-counting efficiency in single Ereco bins only
  // neglecting correlations (for validation purposes only)
  else if(systematic == "sk_ccqe_ring_counting_eff_0") {
     this->GetSuperkMC()->TweakRingCountingEfficiencyCCQE(true,twkdial,0.,0.); 
  }
  else if(systematic == "sk_ccqe_ring_counting_eff_1") {
     this->GetSuperkMC()->TweakRingCountingEfficiencyCCQE(true,0.,twkdial,0.); 
  }
  else if(systematic == "sk_ccqe_ring_counting_eff_2") {
     this->GetSuperkMC()->TweakRingCountingEfficiencyCCQE(true,0.,0.,twkdial); 
  }
  else
  if(systematic == "sk_ccnonqe_eff") {
     this->GetSuperkMC()->TweakEfficiencyCCnonQE(twkdial);
  }
  else
  if(systematic == "sk_nc_eff") {
     this->GetSuperkMC()->TweakEfficiencyNC(twkdial);
  }
  else
  if(systematic == "sk_nuecc_eff") {
     this->GetSuperkMC()->TweakEfficiencyNuECC(twkdial);
  }
  else
  if(systematic == "sk_1Resig_eff") {
     this->GetSuperkMC()->TweakEfficiency1ReSig(twkdial);
  }
  else
  if(systematic == "sk_1Rebkg_eff") {
     this->GetSuperkMC()->TweakEfficiency1ReBkg(twkdial);
  }
  //
  // ND280 detector systematics & statistical err on measured CC event count
  //
  else
  if(systematic == "nd280_eff") {
     this->GetND280Measurement()->TweakEfficiency(twkdial);
  }
  else
  if(systematic == "nd280_norm") {
     this->GetND280Measurement()->FluctuateNObsCCLike(twkdial);
  }


  //
  // flux systematics
  //
  else
  if(systematic == "flux_norm") {
     this->GetNuFlux()->TweakNorm(twkdial);
  }
  else
  if(systematic == "flux_shape") {
     this->GetNuFlux()->TweakShape(); //< Cholesky method
  }

  //
  // cross section systematics
  //
  else
  if(systematic == "xsec_ccqe") {
     this->GetNuXSec()->TweakCCQE(twkdial);
  }
  else
  if(systematic == "xsec_cc1pi_to_ccqe") {
     this->GetNuXSec()->TweakCC1piToCCQE(twkdial);
  }
  else
  if(systematic == "xsec_cccoh_to_ccqe") {
     this->GetNuXSec()->TweakCCCohToCCQE(twkdial);
  }
  else
  if(systematic == "xsec_ccoth_to_ccqe") {
     this->GetNuXSec()->TweakCCOthToCCQE(twkdial);
  }
  else
  if(systematic == "xsec_nc_to_ccqe") {
     this->GetNuXSec()->TweakNCToCCQE(twkdial);
  }
  else
  if(systematic == "xsec_nc1pi0_to_ccqe") {
     this->GetNuXSec()->TweakNC1pi0ToCCQE(twkdial);
  }
  else
  if(systematic == "xsec_nccoh_to_ccqe") {
     this->GetNuXSec()->TweakNCCohToCCQE(twkdial);
  }
  else
  if(systematic == "xsec_ncoth_to_ccqe") {
     this->GetNuXSec()->TweakNCOthToCCQE(twkdial);
  }
  else
  if(systematic == "xsec_nue_to_numu") {
     this->GetNuXSec()->TweakNuEToNuMu(twkdial);
  }
  else
  if(systematic == "nuclear_effect") {
     this->GetNuXSec()->TweakSKFSI(twkdial);
  }   
     
  else {
     cerr << "Exp >> Unknown systematic parameter!" << endl;
     exit(1);
  }   
     
}
//________________________________________________________________________________________
void SuperkExpectation::TweakAllSystematics(bool do_tweak)
{
  if(do_tweak) 
  {
     cout << "Exp >> Generating random variations of all systematic params." << endl;

     TRandom3 rand(0);

     // SuperK detector response
     this->GetSuperkMC()->TweakEnergyScale                 (rand.Gaus(0.,1.));
     this->GetSuperkMC()->TweakOtherEfficiencyCCQE         (rand.Gaus(0.,1.));
     this->GetSuperkMC()->TweakRingCountingEfficiencyCCQE  (); //< Cholesky method
     this->GetSuperkMC()->TweakEfficiencyCCnonQE           (rand.Gaus(0.,1.));
     this->GetSuperkMC()->TweakEfficiencyNC                (rand.Gaus(0.,1.));     
     this->GetSuperkMC()->TweakEfficiencyNuECC             (rand.Gaus(0.,1.));
     this->GetSuperkMC()->TweakEfficiency1ReSig           (rand.Gaus(0.,1.));
     this->GetSuperkMC()->TweakEfficiency1ReBkg           (rand.Gaus(0.,1.));

     // ND280 detector response
     this->GetND280Measurement()->TweakEfficiency     (rand.Gaus(0.,1.));
     this->GetND280Measurement()->FluctuateNObsCCLike (rand.Gaus(0.,1.));

     // neutrino flux
     this->GetNuFlux()->TweakNorm (rand.Gaus(0.,1.));
     this->GetNuFlux()->TweakShape(); //< Cholesky method

     // neutrino cross sections
     this->GetNuXSec()->TweakCCQE        (rand.Gaus(0.,1.));
     this->GetNuXSec()->TweakCC1piToCCQE (rand.Gaus(0.,1.));
     this->GetNuXSec()->TweakCCCohToCCQE (rand.Gaus(0.,1.));
     this->GetNuXSec()->TweakCCOthToCCQE (rand.Gaus(0.,1.));
     this->GetNuXSec()->TweakNCToCCQE    (rand.Gaus(0.,1.));
     this->GetNuXSec()->TweakNC1pi0ToCCQE(rand.Gaus(0.,1.));
     this->GetNuXSec()->TweakNCCohToCCQE (rand.Gaus(0.,1.));
     this->GetNuXSec()->TweakNCOthToCCQE (rand.Gaus(0.,1.));
     this->GetNuXSec()->TweakNuEToNuMu   (rand.Gaus(0.,1.));
     this->GetNuXSec()->TweakSKFSI       (rand.Gaus(0.,1.));

  } else {
    cout << "Exp >> No SuperkExpectation tweak" << endl;
  }
}
//________________________________________________________________________________________
void SuperkExpectation::PrintSystematicVariations(void) const
{
  cout << "Exp >> SuperK PDF tweaks: " << endl;
  fND280Measurement -> PrintSystematicVariations();
  fSuperkMC         -> PrintSystematicVariations();
  fNuFlux           -> PrintSystematicVariations();
  fNuXSec           -> PrintSystematicVariations();
}
//________________________________________________________________________________________
bool SuperkExpectation::Read(string config_file)
{
  cout << "Exp >> Reading config file: " << config_file << endl;

  if(fIsConfigured) {
    this->CleanUp();
  } 

  fND280Measurement  = new ND280Measurement;
  fSuperkMC          = new SuperkMC;
  fNuFlux            = new NuFlux; 
  fNuXSec            = new NuXSec; 

  fND280Measurement  -> TellMeAboutCrossSections (fNuXSec);
  fSuperkMC          -> TellMeAboutCrossSections (fNuXSec);
  fSuperkMC          -> TellMeAboutTheFlux       (fNuFlux);

  for(int is = 0; is < kns; is++) {
     SuperkSample_t s = kSKSamples[is];
     if(s == kSKsUndefined) continue;
     if(s == kSKsFCFV1RingMuLike ) 
     {
       fPredictErecoSpectrumTot[s] = new TH1D("","", kNErecoBins1Rmu, kErecoBinEdges1Rmu);
     }
     else
     if(s == kSKsFCFV1RingELike)  
     {
       fPredictErecoSpectrumTot[s] = new TH1D("","", kNErecoBins1Re, kErecoBinEdges1Re);
     }
     fPredictErecoSpectrumTot[s]->SetDirectory(0);
     for(int im = 0; im < knm; im++) {
         SuperkMode_t m = kSKReactionModes[is][im];
         if(m == kSKmUndefined) continue;
         if(s == kSKsFCFV1RingMuLike ) 
         {
           fPredictErecoSpectrum[s][m] = new TH1D("","", kNErecoBins1Rmu, kErecoBinEdges1Rmu);
         }
         else
         if(s == kSKsFCFV1RingELike)  
         { 
           fPredictErecoSpectrum[s][m] = new TH1D("","", kNErecoBins1Re, kErecoBinEdges1Re);
         }    
         fPredictErecoSpectrum[s][m]->SetDirectory(0);
     }//m	
  }//s

  // Parse the input configuration file
  TXMLEngine* xml = new TXMLEngine;   
  XMLDocPointer_t xmldoc = xml->ParseFile(config_file.c_str());
  if (xmldoc==0) {
    cerr << "Exp >> " << "Can not parse input: " << config_file << endl;
    delete xml;
    return false;  
  }
  XMLNodePointer_t root_node = xml->DocGetRootElement(xmldoc);

  if ( ! this->GetND280Measurement()->Read(xml, root_node) ) return false;
  if ( ! this->GetSuperkMC        ()->Read(xml, root_node) ) return false;
  if ( ! this->GetNuFlux          ()->Read(xml, root_node) ) return false;
  if ( ! this->GetNuXSec          ()->Read(xml, root_node) ) return false;

  // Release memory
  xml->FreeDoc(xmldoc);
  delete xml;

  fIsConfigured = true;

  return true;
}
//________________________________________________________________________________________
void SuperkExpectation::Init(void)
{
  fIsConfigured    = false;

  fND280Measurement = 0;
  fSuperkMC         = 0;
  fNuFlux           = 0;
  fNuXSec           = 0;

  for(int is=0; is<kns; is++) {
    fPredictErecoSpectrumTot[is] = 0;
    for(int im=0; im<knm; im++) {
      fPredictErecoSpectrum[is][im] = 0;
    }//im
  }//is
}
//________________________________________________________________________________________
void SuperkExpectation::CleanUp(void)
{
  if (fND280Measurement) delete fND280Measurement;
  if (fSuperkMC        ) delete fSuperkMC;
  if (fNuFlux          ) delete fNuFlux; 
  if (fNuXSec          ) delete fNuXSec;

  for(int is = 0; is < kns; is++) {
     SuperkSample_t s = kSKSamples[is];
     if(s == kSKsUndefined) continue;
     if (fPredictErecoSpectrumTot[s]) delete fPredictErecoSpectrumTot[s];
     for(int im = 0; im < knm; im++) {
       SuperkMode_t m = kSKReactionModes[is][im];
       if(m == kSKmUndefined) continue;
       if (fPredictErecoSpectrum[s][m]) delete fPredictErecoSpectrum[s][m];
     }
  }
}
//________________________________________________________________________________________
void SuperkExpectation::SetNuisanceParams(const FitParamList & fp)
{
   //this->ResetSystematicVariations();

   //Note that this function has NOT been updated with the new systematics errors for the nue analysis

   if(fp.Exists(kNuisSuperkEffCCQEoth)) {
     double twk = fp.CurrValue(kNuisSuperkEffCCQEoth);
     this->GetSuperkMC()->TweakOtherEfficiencyCCQE(twk);
   }
   if(fp.Exists(kNuisSuperkEffCCQErc0u) &&
      fp.Exists(kNuisSuperkEffCCQErc1u) &&
      fp.Exists(kNuisSuperkEffCCQErc2u)) {
     double twk0 = fp.CurrValue(kNuisSuperkEffCCQErc0u);
     double twk1 = fp.CurrValue(kNuisSuperkEffCCQErc1u);
     double twk2 = fp.CurrValue(kNuisSuperkEffCCQErc2u);
     this->GetSuperkMC()->TweakRingCountingEfficiencyCCQE(true,twk0,twk1,twk2,false);
   }
   if(fp.Exists(kNuisSuperkEffCCnonQE)) {
     double twk = fp.CurrValue(kNuisSuperkEffCCnonQE);
     this->GetSuperkMC()->TweakEfficiencyCCnonQE(twk);
   }
   if(fp.Exists(kNuisSuperkEffNC)) {
     double twk = fp.CurrValue(kNuisSuperkEffNC);
     this->GetSuperkMC()->TweakEfficiencyNC(twk);     
   }
   if(fp.Exists(kNuisSuperkEffNueCC)) {
     double twk = fp.CurrValue(kNuisSuperkEffNueCC);
     this->GetSuperkMC()->TweakEfficiencyNuECC(twk);
   }
   if(fp.Exists(kNuisSuperkEnergyScale)) {
     double twk = fp.CurrValue(kNuisSuperkEnergyScale);
     this->GetSuperkMC()->TweakEnergyScale(twk);
   }
   if(fp.Exists(kNuisNd280Eff)) {
     double twk = fp.CurrValue(kNuisNd280Eff);
     this->GetND280Measurement()->TweakEfficiency(twk);
   }
   if(fp.Exists(kNuisFluxNorm)) {
     double twk = fp.CurrValue(kNuisFluxNorm);
     this->GetNuFlux()->TweakNorm(twk);
   }
   if(fp.Exists(kNuisFluxShapeNuMu0u ) &&
      fp.Exists(kNuisFluxShapeNuMu1u ) &&
      fp.Exists(kNuisFluxShapeNuMu2u ) &&
      fp.Exists(kNuisFluxShapeNuMu3u ) &&
      fp.Exists(kNuisFluxShapeNuMu4u ) &&
      fp.Exists(kNuisFluxShapeNuMu5u ) &&
      fp.Exists(kNuisFluxShapeNuMu6u ) &&
      fp.Exists(kNuisFluxShapeNuMu7u ) &&
      fp.Exists(kNuisFluxShapeNuMu8u ) &&
      fp.Exists(kNuisFluxShapeNuMu9u ) &&
      fp.Exists(kNuisFluxShapeNuMu10u) &&
      fp.Exists(kNuisFluxShapeNuMu11u) &&
      fp.Exists(kNuisFluxShapeNuMu12u) &&
      fp.Exists(kNuisFluxShapeNuMu13u) &&
      fp.Exists(kNuisFluxShapeNuMu14u) &&
      fp.Exists(kNuisFluxShapeNuMu15u) &&
      fp.Exists(kNuisFluxShapeNuMu16u) &&
      fp.Exists(kNuisFluxShapeNuMu17u) &&
      fp.Exists(kNuisFluxShapeNuMu18u)) {
     TVectorD twk(19);
     twk[0]  = fp.CurrValue(kNuisFluxShapeNuMu0u);
     twk[1]  = fp.CurrValue(kNuisFluxShapeNuMu1u);
     twk[2]  = fp.CurrValue(kNuisFluxShapeNuMu2u);
     twk[3]  = fp.CurrValue(kNuisFluxShapeNuMu3u);
     twk[4]  = fp.CurrValue(kNuisFluxShapeNuMu4u);
     twk[5]  = fp.CurrValue(kNuisFluxShapeNuMu5u);
     twk[6]  = fp.CurrValue(kNuisFluxShapeNuMu6u);
     twk[7]  = fp.CurrValue(kNuisFluxShapeNuMu7u);
     twk[8]  = fp.CurrValue(kNuisFluxShapeNuMu8u);
     twk[9]  = fp.CurrValue(kNuisFluxShapeNuMu9u);
     twk[10] = fp.CurrValue(kNuisFluxShapeNuMu10u);
     twk[11] = fp.CurrValue(kNuisFluxShapeNuMu11u);
     twk[12] = fp.CurrValue(kNuisFluxShapeNuMu12u);
     twk[13] = fp.CurrValue(kNuisFluxShapeNuMu13u);
     twk[14] = fp.CurrValue(kNuisFluxShapeNuMu14u);
     twk[15] = fp.CurrValue(kNuisFluxShapeNuMu15u);
     twk[16] = fp.CurrValue(kNuisFluxShapeNuMu16u);
     twk[17] = fp.CurrValue(kNuisFluxShapeNuMu17u);
     twk[18] = fp.CurrValue(kNuisFluxShapeNuMu18u);
     bool force=true;
     bool correlated=false;
     this->GetNuFlux()->TweakShape(force,&twk,correlated);
   }
   if(fp.Exists(kNuisXSecCCQE)) {
     double twk = fp.CurrValue(kNuisXSecCCQE);
     this->GetNuXSec()->TweakCCQE(twk);
   }
   if(fp.Exists(kNuisXSecCC1piToCCQE)) {
     double twk = fp.CurrValue(kNuisXSecCC1piToCCQE);
     this->GetNuXSec()->TweakCC1piToCCQE(twk);
   }
   if(fp.Exists(kNuisXSecCCothToCCQE)) {
     double twk = fp.CurrValue(kNuisXSecCCothToCCQE);
     this->GetNuXSec()->TweakCCOthToCCQE(twk);
   }
   if(fp.Exists(kNuisXSecNCToCCQE)) {
     double twk = fp.CurrValue(kNuisXSecNCToCCQE);
     this->GetNuXSec()->TweakNCToCCQE(twk);
   }
   if(fp.Exists(kNuisXSecNuclEffect)) {
     double twk = fp.CurrValue(kNuisXSecNuclEffect);
     this->GetNuXSec()->TweakSKFSI(twk);
   }

}
//________________________________________________________________________________________

