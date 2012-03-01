//________________________________________________________________________________________
/*

 >> ChangeLog:

 For changes up to version 2010_v3.0.0 (last implementation of the 2010-era analysis)
 see code in the corresponding SVN branch.

 Changes since 2011_v0.9.0:

 @

*/
//________________________________________________________________________________________

#define __show_syst_reset_mesg__
#define __show_syst_tweak_mesg__

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
#include "SuperkMC.h"
#include "NuFlux.h"
#include "NuXSec.h"
#include "BANFF.h"
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

  // update SK MC templates
  fSuperkMC->UpdateTemplates();

  // update Ereco spectrum prediction
  fSuperkMC->UpdateErecoSpectra();

  // copy histograms
  for(int is = 0; is < kns; is++) {
    SuperkSample_t s = kSKSamples[is];
    if(s == kSKsUndefined) continue;
    const int nbins = util::NumErecoBins(s);
    for(int ibin = 1; ibin <= nbins; ibin++) {
       double ns = fSuperkMC->CurrErecoSpectrum(s)->GetBinContent(ibin);
       this->CurrPredictedErecoSpectrum(s)->SetBinContent(ibin, ns);
       for(int im = 0; im < knm; im++) {
           SuperkMode_t m = kSKReactionModes[is][im];
           if(m == kSKmUndefined) continue;
           double nsm = fSuperkMC->CurrErecoSpectrum(s,m)->GetBinContent(ibin);
           this->CurrPredictedErecoSpectrum(s,m)->SetBinContent(ibin, nsm);       
       }//m	
    }//ibin
  }//s

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
#ifdef __show_syst_reset_mesg__
  cout << "Exp >> Reseting all systematic variations." << endl;
#endif

  // SuperK detector systematics
  this->GetSuperkMC()->TweakEnergyScale                    (0.);
  this->GetSuperkMC()->TweakEfficiency1RmuCCQEOther        (0.);
  this->GetSuperkMC()->TweakEfficiency1RmuCCQERingCounting (true,0.,0.,0.); 
  this->GetSuperkMC()->TweakEfficiency1RmuCCnonQE          (0.);
  this->GetSuperkMC()->TweakEfficiency1RmuNC               (0.);     
  this->GetSuperkMC()->TweakEfficiency1RmuNuECC            (0.);

  // BANFF
  this->GetBANFF()->Tweak(true,0);

  // ...
}
//________________________________________________________________________________________
void SuperkExpectation::TweakSystematic(string systematic, double twkdial)
{
#ifdef __show_syst_tweak_mesg__
  cout << "Exp >> Setting systematic parameter: " << systematic 
       << " to value = " << twkdial << endl;
#endif

// Valid systematic parameter names:
// 
// SK detector systematics
//    - sk_energy_scale                     : 
//    - sk_1Rmu_ccqe_ring_counting_eff      : tweak SK 1Rmu/CCQE RC efficiencies in all 3 Ereco bins, correlations taken into account (Cholesky method)
//      or
//      - sk_1Rmu_ccqe_ring_counting_eff_0  : tweak only SK 1Rmu/CCQE RC efficiency in the 1st Ereco bin, correlations with other bins ignored (debug mode)
//      - sk_1Rmu_ccqe_ring_counting_eff_1  : tweak only SK 1Rmu/CCQE RC efficiency in the 2nd Ereco bin, correlations with other bins ignored (debug mode)
//      - sk_1Rmu_ccqe_ring_counting_eff_2  : tweak only SK 1Rmu/CCQE RC efficiency in the 3rd Ereco bin, correlations with other bins ignored (debug mode)
//   - sk_1Rmu_ccqe_other_eff              : 
//   - sk_1Rmu_ccnonqe_eff                  :
//   - sk_1Rmu_nc_eff                       :
//   - sk_1Rmu_nuecc_eff                    :
//
// BANFF
//    - banff                               : tweak all BANFF params, correlations taken into account (Cholesky method)
//      or   
//      - banff_0                           : tweak only BANFF param 00, see Conventions.h for definition, correlations ignored (debug mode)
//      - banff_1                           : tweak only BANFF param 01, see Conventions.h for definition, correlations ignored (debug mode)
//      - ..     
//      - banff_25                          : tweak only BANFF param 25, see Conventions.h for definition, correlations ignored (debug mode)
//

  //
  // SK detector systematics
  //

  if(systematic == "sk_energy_scale") {
     this->GetSuperkMC()->TweakEnergyScale(twkdial);
  }
  else
  if(systematic == "sk_1Rmu_ccqe_other_eff") {
     this->GetSuperkMC()->TweakEfficiency1RmuCCQEOther(twkdial);
  }
  else
  // tweak CCQE energy-dependent ring-counting efficiency generating random
  // efficiency tweaks and taking into account bin-to-bin correlations
  if(systematic == "sk_1Rmu_ccqe_ring_counting_eff") {
     this->GetSuperkMC()->TweakEfficiency1RmuCCQERingCounting(); //< Cholesky method
  }
  //... or tweak the CCQE ring-counting efficiency in single Ereco bins only
  // neglecting correlations (for validation purposes only)
  else if(systematic == "sk_1Rmu_ccqe_ring_counting_eff_0") {
     this->GetSuperkMC()->TweakEfficiency1RmuCCQERingCounting(true,twkdial,0.,0.); 
  }
  else if(systematic == "sk_1Rmu_ccqe_ring_counting_eff_1") {
     this->GetSuperkMC()->TweakEfficiency1RmuCCQERingCounting(true,0.,twkdial,0.); 
  }
  else if(systematic == "sk_1Rmu_ccqe_ring_counting_eff_2") {
     this->GetSuperkMC()->TweakEfficiency1RmuCCQERingCounting(true,0.,0.,twkdial); 
  }
  else
  if(systematic == "sk_1Rmu_ccnonqe_eff") {
     this->GetSuperkMC()->TweakEfficiency1RmuCCnonQE(twkdial);
  }
  else
  if(systematic == "sk_1Rmu_nc_eff") {
     this->GetSuperkMC()->TweakEfficiency1RmuNC(twkdial);
  }
  else
  if(systematic == "sk_1Rmu_nuecc_eff") {
     this->GetSuperkMC()->TweakEfficiency1RmuNuECC(twkdial);
  }
     
  //
  // BANFF 
  //

  // tweak all BANFF parameters taking into account correlations
  else
  if(systematic == "banff") {
     this->GetBANFF()->Tweak(); //< Cholesky method
  }


  else
  if(systematic == "banff_00") {
     TVectorD twkv(26);
     twkv[0]=twkdial;
     this->GetBANFF()->Tweak(true,&twkv); 
  }
  else
  if(systematic == "banff_01") {
     TVectorD twkv(26);
     twkv[1]=twkdial;
     this->GetBANFF()->Tweak(true,&twkv); 
  }
  //
  // etc etc
  // do something more clever here:
  //

  else {
     cerr << "Exp >> Unknown systematic parameter (" << systematic << ")" << endl;
     exit(1);
  }   
     
}
//________________________________________________________________________________________
void SuperkExpectation::TweakAllSystematics(void)
{
#ifdef __show_syst_tweak_mesg__
  cout << "Exp >> Generating random variations of all systematic params." << endl;
#endif

  TRandom3 rand(0);

  // SuperK detector systematics
  this->GetSuperkMC()->TweakEnergyScale              (rand.Gaus(0.,1.));
  this->GetSuperkMC()->TweakEfficiency1RmuCCQEOther  (rand.Gaus(0.,1.));
  this->GetSuperkMC()->TweakEfficiency1RmuCCnonQE    (rand.Gaus(0.,1.));
  this->GetSuperkMC()->TweakEfficiency1RmuNC         (rand.Gaus(0.,1.));     
  this->GetSuperkMC()->TweakEfficiency1RmuNuECC      (rand.Gaus(0.,1.));
  this->GetSuperkMC()->TweakEfficiency1RmuCCQERingCounting  (); // <-- Cholesky method

  // BANFF
  this->GetBANFF()->Tweak(); // <-- Cholesky method

  // ...
}
//________________________________________________________________________________________
void SuperkExpectation::PrintSystematicVariations(void) const
{
  cout << "Exp >> SuperK PDF tweaks: " << endl;
  fSuperkMC -> PrintSystematicVariations();
  fBANFF    -> PrintSystematicVariations();
//fNuXSec   -> PrintSystematicVariations();
}
//________________________________________________________________________________________
bool SuperkExpectation::Read(string config_file)
{
  cout << "Exp >> Reading config file: " << config_file << endl;

  if(fIsConfigured) {
    this->CleanUp();
  } 

  // Instantiate all necessary histograms
  this->Book();

  // Create components used for building-up the p.d.f.
  fBANFF    = new BANFF;
  fNuFlux   = new NuFlux; 
  fNuXSec   = new NuXSec; 
  fSuperkMC = new SuperkMC;

  // Configure inter-dependencies.
  fNuFlux   -> TellMeAboutTheBANFFfit   (fBANFF);
  fNuXSec   -> TellMeAboutTheBANFFfit   (fBANFF);
  fSuperkMC -> TellMeAboutCrossSections (fNuXSec);
  fSuperkMC -> TellMeAboutTheFlux       (fNuFlux);

  // Parse the input XML configuration file
  TXMLEngine* xml = new TXMLEngine;   
  XMLDocPointer_t xmldoc = xml->ParseFile(config_file.c_str());
  if (xmldoc==0) {
    cerr << "Exp >> " << "Can not parse input: " << config_file << endl;
    delete xml;
    return false;  
  }
  XMLNodePointer_t root_node = xml->DocGetRootElement(xmldoc);

  // Configure each of the components 
  if ( ! this->GetBANFF    ()->Read(xml, root_node) ) return false;
  if ( ! this->GetNuFlux   ()->Read(xml, root_node) ) return false;
  if ( ! this->GetNuXSec   ()->Read(xml, root_node) ) return false;
  if ( ! this->GetSuperkMC ()->Read(xml, root_node) ) return false;

  fIsConfigured = true;

  // Release memory
  xml->FreeDoc(xmldoc);
  delete xml;

  cout << "Exp >> Done configuring Superk spectrum predictor... \n\n" << endl;

  return true;
}
//________________________________________________________________________________________
void SuperkExpectation::Init(void)
{
  fIsConfigured = false;

  fSuperkMC = 0;
  fNuFlux   = 0;
  fNuXSec   = 0;
  fBANFF    = 0;

  for(int is=0; is<kns; is++) {
    fPredictErecoSpectrumTot[is] = 0;
    for(int im=0; im<knm; im++) {
      fPredictErecoSpectrum[is][im] = 0;
    }//im
  }//is
}
//________________________________________________________________________________________
void SuperkExpectation::Book(void)
{
  for(int is = 0; is < kns; is++) {
     SuperkSample_t s = kSKSamples[is];
     if(s == kSKsUndefined) continue;
     const int      n_ereco_bins    = util::NumErecoBins (s);
     const double * ereco_bin_edges = util::ErecoBinEdges(s);
     fPredictErecoSpectrumTot[is] = new TH1D("","", n_ereco_bins, ereco_bin_edges);
     fPredictErecoSpectrumTot[is]->SetDirectory(0);
     fPredictErecoSpectrumTot[is]->SetLineWidth(2);
     fPredictErecoSpectrumTot[is]->SetLineColor(kBlack);
     for(int im = 0; im < knm; im++) {
         SuperkMode_t m = kSKReactionModes[is][im];
         if(m == kSKmUndefined) continue;
         fPredictErecoSpectrum[is][im] = new TH1D("","", n_ereco_bins, ereco_bin_edges);
         fPredictErecoSpectrum[is][im]->SetDirectory(0);
         util::SetDefaultSKSpectrumStyle(s,m,fPredictErecoSpectrum[is][im]);
     }//m	
  }//s
}
//________________________________________________________________________________________
void SuperkExpectation::CleanUp(void)
{
  if (fSuperkMC ) delete fSuperkMC;
  if (fNuFlux   ) delete fNuFlux; 
  if (fNuXSec   ) delete fNuXSec;
  if (fBANFF    ) delete fBANFF;

  for(int is = 0; is < kns; is++) {
     if (fPredictErecoSpectrumTot[is]) delete fPredictErecoSpectrumTot[is];
     for(int im = 0; im < knm; im++) {
       if (fPredictErecoSpectrum[is][im]) delete fPredictErecoSpectrum[is][im];
     }
  }
}
//________________________________________________________________________________________
void SuperkExpectation::SetNuisanceParams(const FitParamList & fp)
{
   //this->ResetSystematicVariations();

   // 
   // SuperK detector systematics
   //

   if(fp.Exists(kNuisSuperkEff1RmuCCQEoth)) {
     double twk = fp.CurrValue(kNuisSuperkEff1RmuCCQEoth);
     this->GetSuperkMC()->TweakEfficiency1RmuCCQEOther(twk);
   }
   if(fp.Exists(kNuisSuperkEff1RmuCCQErc0u) &&
      fp.Exists(kNuisSuperkEff1RmuCCQErc1u) &&
      fp.Exists(kNuisSuperkEff1RmuCCQErc2u)) {
     double twk0 = fp.CurrValue(kNuisSuperkEff1RmuCCQErc0u);
     double twk1 = fp.CurrValue(kNuisSuperkEff1RmuCCQErc1u);
     double twk2 = fp.CurrValue(kNuisSuperkEff1RmuCCQErc2u);
     this->GetSuperkMC()->TweakEfficiency1RmuCCQERingCounting(true,twk0,twk1,twk2,false);
   }
   if(fp.Exists(kNuisSuperkEff1RmuCCnonQE)) {
     double twk = fp.CurrValue(kNuisSuperkEff1RmuCCnonQE);
     this->GetSuperkMC()->TweakEfficiency1RmuCCnonQE(twk);
   }
   if(fp.Exists(kNuisSuperkEff1RmuNC)) {
     double twk = fp.CurrValue(kNuisSuperkEff1RmuNC);
     this->GetSuperkMC()->TweakEfficiency1RmuNC(twk);     
   }
   if(fp.Exists(kNuisSuperkEff1RmuNueCC)) {
     double twk = fp.CurrValue(kNuisSuperkEff1RmuNueCC);
     this->GetSuperkMC()->TweakEfficiency1RmuNuECC(twk);
   }


}
//________________________________________________________________________________________

