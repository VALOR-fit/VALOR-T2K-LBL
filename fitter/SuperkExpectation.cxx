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
//#define __show_syst_tweak_debug__

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
#include "SuperkDetector.h"
#include "SuperkMC.h"
#include "NuFlux.h"
#include "NuXSec.h"
#include "BANFF.h"
#include "FSI.h"
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
  int sidx = util::ArrIdx::Instance()->SuperkSampleIdx(s);
  if(sidx == -1) {
    return 0;
  }
  return fPredictErecoSpectrumTot[sidx];  
}
//________________________________________________________________________________________
TH1D* SuperkExpectation::CurrPredictedErecoSpectrum(
  t2kana::conv::SuperkSample_t s, t2kana::conv::SuperkMode_t m) const 
{ 
  int sidx = util::ArrIdx::Instance()->SuperkSampleIdx(s);
  if(sidx == -1) {
    return 0;
  }
  int midx = util::ArrIdx::Instance()->SuperkModeIdx(s,m);
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
  this->GetSuperkDetector()->TweakEnergyScale                    (0.);
  this->GetSuperkDetector()->TweakEfficiency1RmuCCQEOther        (0.);
  this->GetSuperkDetector()->TweakEfficiency1RmuCCQERingCounting (true,0.,0.,0.); 
  this->GetSuperkDetector()->TweakEfficiency1RmuCCnonQE          (0.);
  this->GetSuperkDetector()->TweakEfficiency1RmuNC               (0.);     
  this->GetSuperkDetector()->TweakEfficiency1RmuNuECC            (0.);

  // BANFF
  this->GetBANFF()->Tweak(true,0);

  // FSI
  this->GetFSI()->Tweak(true,0);

  // X Section
  this->GetNuXSec()->TweakFermiMom      (0.);
  this->GetNuXSec()->TweakSpecFunc      (0.);
  this->GetNuXSec()->TweakWshape        (0.);
  this->GetNuXSec()->TweakShapeCCoth    (0.);
  this->GetNuXSec()->TweakPiLessDelta   (0.);
  this->GetNuXSec()->TweakENuShape      (0.);
#if BANFF_BINS != 2
  this->GetNuXSec()->TweakNormCCQEE1    (0.);
  this->GetNuXSec()->TweakNormCCQEE2    (0.);
  this->GetNuXSec()->TweakNormCC1piE1   (0.);
#endif
  this->GetNuXSec()->TweakNormCCcoh     (0.);
  this->GetNuXSec()->TweakNormNC1piC    (0.);
  this->GetNuXSec()->TweakNormNCcoh     (0.);
  this->GetNuXSec()->TweakNormNCoth     (0.);
  this->GetNuXSec()->TweakNormNuEToNuMu (0.);
  this->GetNuXSec()->TweakNormAntiNuToNu(0.);
}
//________________________________________________________________________________________
void SuperkExpectation::TweakSystematic(string systematic, double twkdial, bool best_banff)
{
#ifdef __show_syst_tweak_mesg__
  cout << "Exp >> Setting systematic parameter: " << systematic 
       << " to value = " << twkdial << endl;
#endif

//
// Valid systematic parameter names:
// 
// SK detector systematics
//    - sk_energy_scale                     : 
//    - sk_1Rmu_ccqe_ring_counting_eff      : tweak SK 1Rmu/CCQE RC efficiencies in all 3 Ereco bins, correlations taken into account (Cholesky method)
//      or
//      * sk_1Rmu_ccqe_ring_counting_eff_0  : tweak only SK 1Rmu/CCQE RC efficiency in the 1st Ereco bin, correlations with other bins ignored (debug mode)
//      * sk_1Rmu_ccqe_ring_counting_eff_1  : tweak only SK 1Rmu/CCQE RC efficiency in the 2nd Ereco bin, correlations with other bins ignored (debug mode)
//      * sk_1Rmu_ccqe_ring_counting_eff_2  : tweak only SK 1Rmu/CCQE RC efficiency in the 3rd Ereco bin, correlations with other bins ignored (debug mode)
//   - sk_1Rmu_ccqe_other_eff              : 
//   - sk_1Rmu_ccnonqe_eff                  :
//   - sk_1Rmu_nc_eff                       :
//   - sk_1Rmu_nuecc_eff                    :
//
// BANFF
//    - banff                               : tweak all BANFF params, correlations taken into account (Cholesky method)
//      or   
//      * banff_00                          : tweak only BANFF param 00, see Conventions.h for definition, correlations ignored (debug mode)
//      * banff_01                          : tweak only BANFF param 01, see Conventions.h for definition, correlations ignored (debug mode)
//      * ..     
//      * banff_30(48)                      : tweak only BANFF param 35(48), see Conventions.h for definition, correlations ignored (debug mode)
//      * banff_maqe                        : tweak only BANFF maqe parameter
//      * banff_mares                       : tweak only BANFF mares parameter
//      * banff_ccqeE0                      : tweak only BANFF ccqeE0 parameter
//      * banff_ccqeE1                      : tweak only BANFF ccqeE1 parameter
//      * banff_ccqeE2                      : tweak only BANFF ccqeE2 parameter
//      * banff_cc1pi0E0                    : tweak only BANFF cc1piE0 parameter
//      * banff_cc1pi0E1                    : tweak only BANFF cc1piE1 parameter
//      * banff_nc1pi0                      : tweak only BANFF nc1pi0 parameter
//
// FSI
//    - fsi                                 : tweak all FSI params, correlations taken into account (Cholesky method)
//      or     
//      * fsi_000                           : tweak only FSI param 00, see Conventions.h for definition, correlations ignored (debug mode)
//      * fsi_001                           : tweak only FSI param 01, see Conventions.h for definition, correlations ignored (debug mode)
//      * ..     
//      * fsi_014                           : tweak only FSI param 14, see Conventions.h for definition, correlations ignored (debug mode)
//
// XSec
//    - xsec_pf
//    - xsec_sf
//    - xsec_wshape
//    - xsec_shape_ccoth
//    - xsec_pilessdelta
//    - xsec_enushape
//    - xsec_norm_ccqeE1
//    - xsec_norm_ccqeE2
//    - xsec_norm_cc1piE1
//    - xsec_norm_cccoh
//    - xsec_norm_nc1piC
//    - xsec_norm_nccoh
//    - xsec_norm_ncoth
//    - xsec_norm_nuetonumu
//    - xsec_norm_antinutonu
//

  //
  // SK detector systematics
  //

  if(systematic == "sk_energy_scale") {
     this->GetSuperkDetector()->TweakEnergyScale(twkdial);
  }
  else
  if(systematic == "sk_1Rmu_ccqe_other_eff") {
     this->GetSuperkDetector()->TweakEfficiency1RmuCCQEOther(twkdial);
  }
  else
  // tweak CCQE energy-dependent ring-counting efficiency generating random
  // efficiency tweaks and taking into account bin-to-bin correlations
  if(systematic == "sk_1Rmu_ccqe_ring_counting_eff") {
     this->GetSuperkDetector()->TweakEfficiency1RmuCCQERingCounting(); //< Cholesky method
  }
  //... or tweak the CCQE ring-counting efficiency in single Ereco bins only
  // neglecting correlations (for validation purposes only)
  else if(systematic == "sk_1Rmu_ccqe_ring_counting_eff_0") {
     this->GetSuperkDetector()->TweakEfficiency1RmuCCQERingCounting(true,twkdial,0.,0.); 
  }
  else if(systematic == "sk_1Rmu_ccqe_ring_counting_eff_1") {
     this->GetSuperkDetector()->TweakEfficiency1RmuCCQERingCounting(true,0.,twkdial,0.); 
  }
  else if(systematic == "sk_1Rmu_ccqe_ring_counting_eff_2") {
     this->GetSuperkDetector()->TweakEfficiency1RmuCCQERingCounting(true,0.,0.,twkdial); 
  }
  //special tweak for comparing to KIT analysis (ring counting + other)
  else if(systematic == "sk_1Rmu_ccqe_ring_counting_eff_0+other") {
     this->GetSuperkDetector()->TweakEfficiency1RmuCCQERingCounting(true,twkdial,0.,0.); 
     this->GetSuperkDetector()->TweakEfficiency1RmuCCQEOther(twkdial);
  }
  else if(systematic == "sk_1Rmu_ccqe_ring_counting_eff_1+other") {
     this->GetSuperkDetector()->TweakEfficiency1RmuCCQERingCounting(true,0.,twkdial,0.); 
     this->GetSuperkDetector()->TweakEfficiency1RmuCCQEOther(twkdial);
  }
  else if(systematic == "sk_1Rmu_ccqe_ring_counting_eff_2+other") {
     this->GetSuperkDetector()->TweakEfficiency1RmuCCQERingCounting(true,0.,0.,twkdial); 
     this->GetSuperkDetector()->TweakEfficiency1RmuCCQEOther(twkdial);
  }
  else
  if(systematic == "sk_1Rmu_ccnonqe_eff") {
     this->GetSuperkDetector()->TweakEfficiency1RmuCCnonQE(twkdial);
  }
  else
  if(systematic == "sk_1Rmu_nc_eff") {
     this->GetSuperkDetector()->TweakEfficiency1RmuNC(twkdial);
  }
  else
  if(systematic == "sk_1Rmu_nuecc_eff") {
     this->GetSuperkDetector()->TweakEfficiency1RmuNuECC(twkdial);
  }
     
  //
  // BANFF 
  //

  // tweak all BANFF parameters taking into account correlations
  else
  if(systematic == "banff") {
     this->GetBANFF()->Tweak(); //< Cholesky method
  }
  // tweak an individual BANFF parameter (ignores correlations, debug)
  else
  if(systematic.find("banff_") != string::npos) {
     const int nbanff = conv::kNumOfParamsInBANFFfit;
     TVectorD twkv(this->GetBANFF()->BestFitValues());
     if(!best_banff) {
        for (int ib = 0; ib < nbanff; ib++) { twkv[ib] = 1; }
     }
     PhysParamId_t p = util::PhysParamIdFromString(systematic);
     int ip = util::ArrIdx::Instance()->BANFFSystIdx(p);
     twkv[ip] = twkdial;
#ifdef __show_syst_tweak_debug__
     cout << "Exp >> Tweaking systematic string: " << systematic
	  << ", NuisParam: " << util::AsString(p)
	  << endl;
     util::PrintCompact(twkv);
#endif
     this->GetBANFF()->Tweak(true,&twkv);

     /*
     for (int ib = 0; ib < nbanff; ib++) {
	if(systematic == Form("banff_%02i",ib)) {
	   twkv[ib] = twkdial;
	   this->GetBANFF()->Tweak(true,&twkv);
	   found = true;
	   break;
	}
     }//ib
     if(!found) {
        cerr << "Exp >> Unknown systematic parameter (" << systematic << ")" << endl;
	exit(1);
     }   
     */
  }


  //
  // FSI
  //

  // tweak all FSI parameters taking into account correlations
  else
  if(systematic == "fsi") {
     this->GetFSI()->Tweak(); //< Cholesky method
  }
  // tweak an individual FSI parameter (ignores correlations, debug)
  else
  if(systematic.find("fsi_") != string::npos) {
     bool found = false;
     const int nfsi = conv::kNumOfFSIParams;
     TVectorD twkv(nfsi);
     for (int ip = 0; ip < nfsi; ip++) {
        twkv[ip] = 1;
     }
     for (int ip = 0; ip < nfsi; ip++) {
	if(systematic == Form("fsi_%03i",ip)) {
	   twkv[ip]=twkdial;
	   this->GetFSI()->Tweak(true,&twkv);
	   found = true;
	   break;
	}
     }//ib
     if(!found) {
        cerr << "Exp >> Unknown systematic parameter (" << systematic << ")" << endl;
	exit(1);
     }   
  }


  //
  // XSec
  //

  else if(systematic == "xsec_pf")
    this->GetNuXSec()->TweakFermiMom      (twkdial);
  else if(systematic == "xsec_sf")
    this->GetNuXSec()->TweakSpecFunc      (twkdial);
  else if(systematic == "xsec_wshape")
    this->GetNuXSec()->TweakWshape        (twkdial);
  else if(systematic == "xsec_shape_ccoth")
    this->GetNuXSec()->TweakShapeCCoth    (twkdial);
  else if(systematic == "xsec_pilessdelta")
    this->GetNuXSec()->TweakPiLessDelta   (twkdial);
  else if(systematic == "xsec_enushape")
    this->GetNuXSec()->TweakENuShape      (twkdial);
#if BANFF_BINS != 2
  else if(systematic == "xsec_norm_ccqeE1")
    this->GetNuXSec()->TweakNormCCQEE1    (twkdial);
  else if(systematic == "xsec_norm_ccqeE2")
    this->GetNuXSec()->TweakNormCCQEE2    (twkdial);
  else if(systematic == "xsec_norm_cc1piE1")
    this->GetNuXSec()->TweakNormCC1piE1   (twkdial);
#endif
  else if(systematic == "xsec_norm_cccoh")
    this->GetNuXSec()->TweakNormCCcoh     (twkdial);
  else if(systematic == "xsec_norm_nc1piC")
    this->GetNuXSec()->TweakNormNC1piC    (twkdial);
  else if(systematic == "xsec_norm_nccoh")
    this->GetNuXSec()->TweakNormNCcoh     (twkdial);
  else if(systematic == "xsec_norm_ncoth")
    this->GetNuXSec()->TweakNormNCoth     (twkdial);
  else if(systematic == "xsec_norm_nuetonumu")
    this->GetNuXSec()->TweakNormNuEToNuMu (twkdial);
  else if(systematic == "xsec_norm_antinutonu")
    this->GetNuXSec()->TweakNormAntiNuToNu(twkdial);


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

  this->TweakAllBANFFSystematics();
  this->TweakAllFSISystematics();
  this->TweakAllUncorrelatedXSecSystematics();
  this->TweakAllSKDetectorSystematics();
}
//________________________________________________________________________________________
void SuperkExpectation::TweakAllBANFFSystematics(void)
{
#ifdef __show_syst_tweak_mesg__
  cout << "Exp >> Generating random variations of all BANFF systematic params." << endl;
#endif

  // BANFF
  this->GetBANFF()->Tweak(); // <-- Cholesky method
}
//________________________________________________________________________________________
void SuperkExpectation::TweakAllFSISystematics(void)
{
#ifdef __show_syst_tweak_mesg__
  cout << "Exp >> Generating random variations of all FSI systematic params." << endl;
#endif

  // FSI
  this->GetFSI()->Tweak(); // <-- Cholesky method
}
//________________________________________________________________________________________
void SuperkExpectation::TweakAllUncorrelatedXSecSystematics(void)
{
#ifdef __show_syst_tweak_mesg__
  cout << "Exp >> Generating random variations of all Uncorrelated XSec systematic params." << endl;
#endif

  TRandom3 rand(0);

  // XSec
  this->GetNuXSec()->TweakFermiMom      (rand.Gaus(0.,1.));
  this->GetNuXSec()->TweakSpecFunc      (rand.Gaus(0.,1.));
  this->GetNuXSec()->TweakWshape        (rand.Gaus(0.,1.));
  this->GetNuXSec()->TweakShapeCCoth    (rand.Gaus(0.,1.));
  this->GetNuXSec()->TweakPiLessDelta   (rand.Gaus(0.,1.));
  this->GetNuXSec()->TweakENuShape      (rand.Gaus(0.,1.));
#if BANFF_BINS != 2
  this->GetNuXSec()->TweakNormCCQEE1    (rand.Gaus(0.,1.));
  this->GetNuXSec()->TweakNormCCQEE2    (rand.Gaus(0.,1.));
  this->GetNuXSec()->TweakNormCC1piE1   (rand.Gaus(0.,1.));
#endif
  this->GetNuXSec()->TweakNormCCcoh     (rand.Gaus(0.,1.));
  this->GetNuXSec()->TweakNormNC1piC    (rand.Gaus(0.,1.));
  this->GetNuXSec()->TweakNormNCcoh     (rand.Gaus(0.,1.));
  this->GetNuXSec()->TweakNormNCoth     (rand.Gaus(0.,1.));
  this->GetNuXSec()->TweakNormNuEToNuMu (rand.Gaus(0.,1.));
  this->GetNuXSec()->TweakNormAntiNuToNu(rand.Gaus(0.,1.));
}
//________________________________________________________________________________________
void SuperkExpectation::TweakAllSKDetectorSystematics(void)
{
#ifdef __show_syst_tweak_mesg__
  cout << "Exp >> Generating random variations of all SKDetector systematic params." << endl;
#endif

  TRandom3 rand(0);

  // SuperK detector systematics
  this->GetSuperkDetector()->TweakEnergyScale              (rand.Gaus(0.,1.));
  this->GetSuperkDetector()->TweakEfficiency1RmuCCQEOther  (rand.Gaus(0.,1.));
  this->GetSuperkDetector()->TweakEfficiency1RmuCCnonQE    (rand.Gaus(0.,1.));
  this->GetSuperkDetector()->TweakEfficiency1RmuNC         (rand.Gaus(0.,1.));     
  this->GetSuperkDetector()->TweakEfficiency1RmuNuECC      (rand.Gaus(0.,1.));
  this->GetSuperkDetector()->TweakEfficiency1RmuCCQERingCounting  (); // <-- Cholesky method
}
//________________________________________________________________________________________
void SuperkExpectation::PrintSystematicVariations(ostream & stream) const
{
  cout << "Exp >> SuperK PDF tweaks: " << endl;
  fSuperkDetector -> PrintSystematicVariations(stream);
  fBANFF          -> PrintSystematicVariations(stream);
  fNuXSec         -> PrintSystematicVariations(stream);
  fFSI            -> PrintSystematicVariations(stream);
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

  // Create components used for building-up the p.d.f & implemeting its systematic variation
  fBANFF          = new BANFF;
  fNuFlux         = new NuFlux; 
  fNuXSec         = new NuXSec; 
  fFSI            = new FSI; 
  fSuperkDetector = new SuperkDetector; 
  fSuperkMC       = new SuperkMC;

  // Configure inter-dependencies.
  fNuFlux   -> TellMeAboutTheBANFFfit   (fBANFF);
  fNuXSec   -> TellMeAboutTheBANFFfit   (fBANFF);
  fSuperkMC -> TellMeAboutTheDetector   (fSuperkDetector);
  fSuperkMC -> TellMeAboutCrossSections (fNuXSec);
  fSuperkMC -> TellMeAboutTheFlux       (fNuFlux);
  fSuperkMC -> TellMeAboutFSI           (fFSI);

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
  if ( ! this->GetBANFF          () -> Read (xml, root_node) ) return false;
  if ( ! this->GetNuFlux         () -> Read (xml, root_node) ) return false;
  if ( ! this->GetNuXSec         () -> Read (xml, root_node) ) return false;
  if ( ! this->GetFSI            () -> Read (xml, root_node) ) return false;
  if ( ! this->GetSuperkDetector () -> Read (xml, root_node) ) return false;
  if ( ! this->GetSuperkMC       () -> Read (xml, root_node) ) return false;

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

  fSuperkDetector = 0;
  fSuperkMC       = 0;
  fNuFlux         = 0;
  fNuXSec         = 0;
  fFSI            = 0;
  fBANFF          = 0;

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
  if (fSuperkDetector ) delete fSuperkDetector;
  if (fSuperkMC       ) delete fSuperkMC;
  if (fNuFlux         ) delete fNuFlux; 
  if (fNuXSec         ) delete fNuXSec;
  if (fBANFF          ) delete fBANFF;
  if (fFSI            ) delete fFSI;

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
     this->GetSuperkDetector()->TweakEfficiency1RmuCCQEOther(twk);
   }
   if(fp.Exists(kNuisSuperkEff1RmuCCQErc0u) &&
      fp.Exists(kNuisSuperkEff1RmuCCQErc1u) &&
      fp.Exists(kNuisSuperkEff1RmuCCQErc2u)) {
     double twk0 = fp.CurrValue(kNuisSuperkEff1RmuCCQErc0u);
     double twk1 = fp.CurrValue(kNuisSuperkEff1RmuCCQErc1u);
     double twk2 = fp.CurrValue(kNuisSuperkEff1RmuCCQErc2u);
     this->GetSuperkDetector()->TweakEfficiency1RmuCCQERingCounting(true,twk0,twk1,twk2,false);
   }
   if(fp.Exists(kNuisSuperkEff1RmuCCnonQE)) {
     double twk = fp.CurrValue(kNuisSuperkEff1RmuCCnonQE);
     this->GetSuperkDetector()->TweakEfficiency1RmuCCnonQE(twk);
   }
   if(fp.Exists(kNuisSuperkEff1RmuNC)) {
     double twk = fp.CurrValue(kNuisSuperkEff1RmuNC);
     this->GetSuperkDetector()->TweakEfficiency1RmuNC(twk);     
   }
   if(fp.Exists(kNuisSuperkEff1RmuNueCC)) {
     double twk = fp.CurrValue(kNuisSuperkEff1RmuNueCC);
     this->GetSuperkDetector()->TweakEfficiency1RmuNuECC(twk);
   }

   // 
   // BANFF systematics
   //

   const int nbanff = conv::kNumOfParamsInBANFFfit;
   // BANFF params included in the fit param list?
   TVectorD twkv(nbanff);
   for(int ip = 0; ip < nbanff; ip++) {
     PhysParamId_t p = PhysParamId_t ((int) kNuisBANFF00 + ip);
     bool incl_banff = fp.Exists(p);
     // If BANFF param included, use current value stored in fit param list
     // to tweak the BANFF object
     if(incl_banff) {
       twkv[ip] = fp.CurrValue(p);
     }
     // Otherwise, use whatever is already stored in the BANFF object
     else {
       twkv[ip] = this->GetBANFF()->CurrentValue(p);
     }
   }
   this->GetBANFF()->Tweak(true,&twkv);

   // 
   // Cross-section systematics (other than BANFF ones)
   //
   if(fp.Exists(kNuisXSecFermiMom)) {
     double twk = fp.CurrValue(kNuisXSecFermiMom);
     this->GetNuXSec()->TweakFermiMom(twk);
   }
   if(fp.Exists(kNuisXSecSpecFunc)) {
     double twk = fp.CurrValue(kNuisXSecSpecFunc);
     this->GetNuXSec()->TweakSpecFunc(twk);
   }
   if(fp.Exists(kNuisXSecWshape)) {
     double twk = fp.CurrValue(kNuisXSecWshape);
     this->GetNuXSec()->TweakWshape(twk);
   }
   if(fp.Exists(kNuisXSecShapeCCoth)) {
     double twk = fp.CurrValue(kNuisXSecShapeCCoth);
     this->GetNuXSec()->TweakShapeCCoth(twk);
   }
   if(fp.Exists(kNuisXSecPiLessDelta)) {
     double twk = fp.CurrValue(kNuisXSecPiLessDelta);
     this->GetNuXSec()->TweakPiLessDelta(twk);
   }
   if(fp.Exists(kNuisXSecENuShape)) {
     double twk = fp.CurrValue(kNuisXSecENuShape);
     this->GetNuXSec()->TweakENuShape(twk);
   }
#if BANFF_BINS != 2
   if(fp.Exists(kNuisXSecNormCCQEE1)) {
     double twk = fp.CurrValue(kNuisXSecNormCCQEE1);
     this->GetNuXSec()->TweakNormCCQEE1(twk);
   }
   if(fp.Exists(kNuisXSecNormCCQEE2)) {
     double twk = fp.CurrValue(kNuisXSecNormCCQEE2);
     this->GetNuXSec()->TweakNormCCQEE2(twk);
   }
   if(fp.Exists(kNuisXSecNormCC1piE1)) {
     double twk = fp.CurrValue(kNuisXSecNormCC1piE1);
     this->GetNuXSec()->TweakNormCC1piE1(twk);
   }
#endif
   if(fp.Exists(kNuisXSecNormCCcoh)) {
     double twk = fp.CurrValue(kNuisXSecNormCCcoh);
     this->GetNuXSec()->TweakNormCCcoh(twk);
   }
   if(fp.Exists(kNuisXSecNormNC1piC)) {
     double twk = fp.CurrValue(kNuisXSecNormNC1piC);
     this->GetNuXSec()->TweakNormNC1piC(twk);
   }
   if(fp.Exists(kNuisXSecNormNCcoh)) {
     double twk = fp.CurrValue(kNuisXSecNormNCcoh);
     this->GetNuXSec()->TweakNormNCcoh(twk);
   }
   if(fp.Exists(kNuisXSecNormNCoth)) {
     double twk = fp.CurrValue(kNuisXSecNormNCoth);
     this->GetNuXSec()->TweakNormNCoth(twk);
   }
   if(fp.Exists(kNuisXSecNormNuEToNuMu)) {
     double twk = fp.CurrValue(kNuisXSecNormNuEToNuMu);
     this->GetNuXSec()->TweakNormNuEToNuMu(twk);
   }
   if(fp.Exists(kNuisXSecNormAntiNuToNu)) {
     double twk = fp.CurrValue(kNuisXSecNormAntiNuToNu);
     this->GetNuXSec()->TweakNormAntiNuToNu(twk);
   }

   // 
   // FSI systematics
   //

   const int nfsi = conv::kNumOfFSIParams;
   // FSI params included in the fit param list?
   bool incl_fsi = true;   
   for(int ip = 0; ip < nfsi; ip++) {
      PhysParamId_t p = PhysParamId_t ((int) kNuisFSI000 + ip);
      incl_fsi = incl_fsi && fp.Exists(p);
      if(!incl_fsi) break;
   }
   // If FSI params included, use current values stored in fit param list
   // to tweak the FSI object
   if(incl_fsi) {
     TVectorD twkv(nfsi);
     for(int ip = 0; ip < nfsi; ip++) {
        PhysParamId_t p = PhysParamId_t ((int) kNuisFSI000 + ip);
        twkv[ip] = fp.CurrValue(p);
     }
     this->GetFSI()->Tweak(true,&twkv);    
   }
}
//________________________________________________________________________________________

