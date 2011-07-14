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
 @ 13/07/2011 - CA
   Take into account the statistical error on the measured number of ND280 CC-like events
   and propagate it to SuperK.
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

  fSamples.push_back(kSKsFCFV1RingMuLike);
}
//________________________________________________________________________________________
SuperkExpectation::SuperkExpectation(const vector<SuperkSample_t> & samples)
{
  this->Init();

  vector<SuperkSample_t>::const_iterator iter = samples.begin();
  for( ; iter != samples.end(); ++iter) {
    fSamples.push_back(*iter);
  }
}
//________________________________________________________________________________________
SuperkExpectation::~SuperkExpectation(void)
{
  this->CleanUp();
}
//________________________________________________________________________________________
void SuperkExpectation::Generate(const PhysParamList & fp)
{
// Calculate dN^{exp}_{SuperK}(E_{reco}) / dE_{reco} for the input E_{reco} and 
// the input physics parameter list.
//

  if(!fIsConfigured) {
    cerr << "Exp >> Not configured properly! " << endl;
  //return -1;
    exit(1);
  }

  const unsigned int ns = fSamples.size();
  const unsigned int nm = 12;  
  SuperkMode_t sk_modes [nm] = { 
     kSKmNuMuCCQE,    kSKmNuMuCC1pi,    kSKmNuMuCCoth,    kSKmNC, 
     kSKmNuMuBarCCQE, kSKmNuMuBarCC1pi, kSKmNuMuBarCCoth, kSKmNuMuBarNC, 
     kSKmNuECCQE,     kSKmNuECC1pi,     kSKmNuECCoth,     kSKmNuENC
  };


  // input physics params (params included in the fit / systematics params)
  if( fp.Exists(kPhypDm23sq) && fp.Exists(kPhypSinsq2Theta23) ) {
    double dmsq23        = fp.CurrValue(kPhypDm23sq);
    double sinsq2theta23 = fp.CurrValue(kPhypSinsq2Theta23);
    fSuperkMC->Hypothesis().EnableOscillations(true);
    fSuperkMC->Hypothesis().EnableDecay(false);
    fSuperkMC->Hypothesis().EnableDecoherence(false);
    fSuperkMC->Hypothesis().UseOscillationParams(sinsq2theta23, dmsq23);
  } else
  if( fp.Exists(kPhypAlphaDecay) && fp.Exists(kPhypSinsq2Theta23) ) {
    double alpha         = fp.CurrValue(kPhypAlphaDecay);
    double sinsq2theta23 = fp.CurrValue(kPhypSinsq2Theta23);
    fSuperkMC->Hypothesis().EnableOscillations(false);
    fSuperkMC->Hypothesis().EnableDecay(true);
    fSuperkMC->Hypothesis().EnableDecoherence(false);
    fSuperkMC->Hypothesis().UseDecayParams(sinsq2theta23, alpha);
  } else
  if( fp.Exists(kPhypMusqDecoherence) && fp.Exists(kPhypSinsq2Theta23) ) {
    double musq          = fp.CurrValue(kPhypMusqDecoherence);
    double sinsq2theta23 = fp.CurrValue(kPhypSinsq2Theta23);
    fSuperkMC->Hypothesis().EnableOscillations(false);
    fSuperkMC->Hypothesis().EnableDecay(false);
    fSuperkMC->Hypothesis().EnableDecoherence(true);
    fSuperkMC->Hypothesis().UseDecoherenceParams(sinsq2theta23, musq);
  }


  double shape_norm = 1.;
  if(this->GetNuFlux()->ShapeTweaked()) {
     this->GetNuFlux()->TweakShape(false); // switch-off beam-shape tweaking
     this->Generate(fp); // generate untweaked predictions
     // get nominal Nsk
     shape_norm = fPredictErecoSpectrumTot[kSKsFCFV1RingMuLike]->Integral(); 
     this->GetNuFlux()->TweakShape(true); // switch-on beam-shape tweaking
  }

  // correction factor from ND280
  double correction = fND280Measurement->CorrectionFactor();

  // update SK MC
  fSuperkMC->UpdateErecoSpectra();

  // update prediction
  for(int ibin = 1; ibin <= kNvErecoBins; ibin++) {
    double n = 0;
    for(unsigned int is = 0; is < ns; is++) {
       SuperkSample_t s = fSamples[is];
       double ns = correction * fSuperkMC->CurrErecoSpectrum(s)->GetBinContent(ibin);
       n += ns;
       fPredictErecoSpectrumTot[s]->SetBinContent(ibin, ns);
       for(unsigned int im = 0; im < nm; im++) {
           SuperkMode_t m = sk_modes[im];
           double nsm = correction * fSuperkMC->CurrErecoSpectrum(s,m)->GetBinContent(ibin);
           fPredictErecoSpectrum[s][m]->SetBinContent(ibin, nsm);       
       }//m	
    }//s
  }//ibin

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
     fPredictErecoSpectrumTot[kSKsFCFV1RingMuLike]                  -> Scale (beam_norm);
     fPredictErecoSpectrum[kSKsFCFV1RingMuLike][kSKmNuMuCCQE]       -> Scale (beam_norm);
     fPredictErecoSpectrum[kSKsFCFV1RingMuLike][kSKmNuMuCC1pi]      -> Scale (beam_norm);
     fPredictErecoSpectrum[kSKsFCFV1RingMuLike][kSKmNuMuCCoth]      -> Scale (beam_norm);
     fPredictErecoSpectrum[kSKsFCFV1RingMuLike][kSKmNC]             -> Scale (beam_norm);
     fPredictErecoSpectrum[kSKsFCFV1RingMuLike][kSKmNuMuBarCCQE]    -> Scale (beam_norm);
     fPredictErecoSpectrum[kSKsFCFV1RingMuLike][kSKmNuMuBarCC1pi]   -> Scale (beam_norm);
     fPredictErecoSpectrum[kSKsFCFV1RingMuLike][kSKmNuMuBarCCoth]   -> Scale (beam_norm);
     fPredictErecoSpectrum[kSKsFCFV1RingMuLike][kSKmNuMuBarNC]      -> Scale (beam_norm);
     fPredictErecoSpectrum[kSKsFCFV1RingMuLike][kSKmNuECCQE]        -> Scale (beam_norm);
     fPredictErecoSpectrum[kSKsFCFV1RingMuLike][kSKmNuECC1pi]       -> Scale (beam_norm);
     fPredictErecoSpectrum[kSKsFCFV1RingMuLike][kSKmNuECCoth]       -> Scale (beam_norm);
     fPredictErecoSpectrum[kSKsFCFV1RingMuLike][kSKmNuENC]          -> Scale (beam_norm);
  }//norm tweaked?

  if(this->GetNuFlux()->ShapeTweaked()) {
     double nsktwk = fPredictErecoSpectrumTot[kSKsFCFV1RingMuLike]->Integral();
     shape_norm /= nsktwk;
     fPredictErecoSpectrumTot[kSKsFCFV1RingMuLike]                  -> Scale (shape_norm);
     fPredictErecoSpectrum[kSKsFCFV1RingMuLike][kSKmNuMuCCQE]       -> Scale (shape_norm);
     fPredictErecoSpectrum[kSKsFCFV1RingMuLike][kSKmNuMuCC1pi]      -> Scale (shape_norm);
     fPredictErecoSpectrum[kSKsFCFV1RingMuLike][kSKmNuMuCCoth]      -> Scale (shape_norm);
     fPredictErecoSpectrum[kSKsFCFV1RingMuLike][kSKmNC]             -> Scale (shape_norm);
     fPredictErecoSpectrum[kSKsFCFV1RingMuLike][kSKmNuMuBarCCQE]    -> Scale (shape_norm);
     fPredictErecoSpectrum[kSKsFCFV1RingMuLike][kSKmNuMuBarCC1pi]   -> Scale (shape_norm);
     fPredictErecoSpectrum[kSKsFCFV1RingMuLike][kSKmNuMuBarCCoth]   -> Scale (shape_norm);
     fPredictErecoSpectrum[kSKsFCFV1RingMuLike][kSKmNuMuBarNC]      -> Scale (shape_norm);
     fPredictErecoSpectrum[kSKsFCFV1RingMuLike][kSKmNuECCQE]        -> Scale (shape_norm);
     fPredictErecoSpectrum[kSKsFCFV1RingMuLike][kSKmNuECC1pi]       -> Scale (shape_norm);
     fPredictErecoSpectrum[kSKsFCFV1RingMuLike][kSKmNuECCoth]       -> Scale (shape_norm);
     fPredictErecoSpectrum[kSKsFCFV1RingMuLike][kSKmNuENC]          -> Scale (shape_norm);
  }

}
//________________________________________________________________________________________
TH1D* SuperkExpectation::CurrPredictedErecoSpectrum(t2kana::conv::SuperkSample_t s) const 
{
  return fPredictErecoSpectrumTot[s];  
}
//________________________________________________________________________________________
TH1D* SuperkExpectation::CurrPredictedErecoSpectrum(
  t2kana::conv::SuperkSample_t s, t2kana::conv::SuperkMode_t m) const 
{ 
  return fPredictErecoSpectrum[s][m]; 
}
//________________________________________________________________________________________
void SuperkExpectation::ResetSystematicVariations(void)
{
  this->GetSuperkMC()->TweakEnergyScale                (0.);
  this->GetSuperkMC()->TweakOtherEfficiencyCCQE        (0.);
  this->GetSuperkMC()->TweakRingCountingEfficiencyCCQE (true,0.,0.,0.); 
  this->GetSuperkMC()->TweakEfficiencyCCnonQE          (0.);
  this->GetSuperkMC()->TweakEfficiencyNC               (0.);     
  this->GetSuperkMC()->TweakEfficiencyNuECC            (0.);

  this->GetND280Measurement()->TweakEfficiency (0.);
  this->GetND280Measurement()->FluctuateNObsCCLike(0.);

  this->GetNuFlux()->TweakNorm (0.);
  this->GetNuFlux()->TweakShape(false);

  this->GetNuXSec()->TweakCCQEXSec             (0.);
  this->GetNuXSec()->TweakCC1piToCCQEXSecRatio (0.);
  this->GetNuXSec()->TweakCCOthToCCQEXSecRatio (0.);
  this->GetNuXSec()->TweakNCToCCQEXSecRatio    (0.);
  this->GetNuXSec()->TweakNuclEffect           (0.);
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

  //
  // ND280 detector systematics
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
     this->GetNuXSec()->TweakCCQEXSec(twkdial);
  }
  else
  if(systematic == "xsec_cc1pi_to_ccqe") {
     this->GetNuXSec()->TweakCC1piToCCQEXSecRatio(twkdial);
  }
  else
  if(systematic == "xsec_ccoth_to_ccqe") {
     this->GetNuXSec()->TweakCCOthToCCQEXSecRatio(twkdial);
  }
  else
  if(systematic == "xsec_nc_to_ccqe") {
     this->GetNuXSec()->TweakNCToCCQEXSecRatio(twkdial);
  }
  else
  if(systematic == "nuclear_effect") {
     this->GetNuXSec()->TweakNuclEffect(twkdial);
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
     TRandom3 rand(0);

     // SuperK detector response
     this->GetSuperkMC()->TweakEnergyScale                 (rand.Gaus(0.,1.));
     this->GetSuperkMC()->TweakOtherEfficiencyCCQE         (rand.Gaus(0.,1.));
     this->GetSuperkMC()->TweakRingCountingEfficiencyCCQE  (); //< Cholesky method
     this->GetSuperkMC()->TweakEfficiencyCCnonQE           (rand.Gaus(0.,1.));
     this->GetSuperkMC()->TweakEfficiencyNC                (rand.Gaus(0.,1.));     
     this->GetSuperkMC()->TweakEfficiencyNuECC             (rand.Gaus(0.,1.));

     // ND280 detector response
     this->GetND280Measurement()->TweakEfficiency     (rand.Gaus(0.,1.));
     this->GetND280Measurement()->FluctuateNObsCCLike (rand.Gaus(0.,1.));

     // neutrino flux
     this->GetNuFlux()->TweakNorm (rand.Gaus(0.,1.));
     this->GetNuFlux()->TweakShape(); //< Cholesky method

     // neutrino cross sections
     this->GetNuXSec()->TweakCCQEXSec             (rand.Gaus(0.,1.));
     this->GetNuXSec()->TweakCC1piToCCQEXSecRatio (rand.Gaus(0.,1.));
     this->GetNuXSec()->TweakCCOthToCCQEXSecRatio (rand.Gaus(0.,1.));
     this->GetNuXSec()->TweakNCToCCQEXSecRatio    (rand.Gaus(0.,1.));
     this->GetNuXSec()->TweakNuclEffect           (rand.Gaus(0.,1.));

  } else {
    cout << "Exp >> No SuperkExpectation tweak" << endl;
  }
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

  const unsigned int ns = fSamples.size();
  const unsigned int nm = 12;  
  SuperkMode_t sk_modes [nm] = { 
     kSKmNuMuCCQE,    kSKmNuMuCC1pi,    kSKmNuMuCCoth,    kSKmNC, 
     kSKmNuMuBarCCQE, kSKmNuMuBarCC1pi, kSKmNuMuBarCCoth, kSKmNuMuBarNC, 
     kSKmNuECCQE,     kSKmNuECC1pi,     kSKmNuECCoth,     kSKmNuENC
  };
  for(unsigned int is = 0; is < ns; is++) {
    SuperkSample_t s = fSamples[is];
    fPredictErecoSpectrumTot[s] = new TH1D("","", kNvErecoBins, kNvErecoBinEdges);
    fPredictErecoSpectrumTot[s]->SetDirectory(0);
    for(unsigned int im = 0; im < nm; im++) {
       SuperkMode_t m = sk_modes[im];
       fPredictErecoSpectrum[s][m] = new TH1D("","", kNvErecoBins, kNvErecoBinEdges);
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

  fSamples.clear();

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

  fSamples.clear();

  for(int is=0; is<kns; is++) {
    if (fPredictErecoSpectrumTot[is]) delete fPredictErecoSpectrumTot[is];
    for(int im=0; im<knm; im++) {
      if (fPredictErecoSpectrum[is][im]) delete fPredictErecoSpectrum[is][im];
    }//im
  }//is
}
//________________________________________________________________________________________

