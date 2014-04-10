//________________________________________________________________________________________
/*

 \copyright (c) 2010-2013, VaLOR group.

30-12-2013  Nick
Added samples for reverse horn current running

*/
//________________________________________________________________________________________

#define __show_syst_reset_mesg__
//#define __show_syst_tweak_mesg__
//#define __show_syst_tweak_debug__

#include <cassert>
#include <cstdlib>
#include <iostream>

#include <TSystem.h>
#include <TFile.h>
#include <TObjArray.h>
#include <TTree.h>
#include <TString.h>
#include <TLegend.h>
#include <TSpline.h>
#include <TRandom3.h>

#include "SuperkExpectation.h"
#include "SuperkDetectorFSI.h"
#include "SuperkMC.h"
#include "NuFlux.h"
#include "NuXSec.h"
#include "BANFF.h"
#include "Conventions.h"
#include "Utilities.h"
#include "Opt.h"
#include "Rnd.h"

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
void SuperkExpectation::Init(void)
{
  fIsConfigured = false;

  fSuperkMC          = 0;
  fSuperkDetectorFSI = 0;
  fNuFlux            = 0;
  fNuXSec            = 0;
  fBANFF             = 0;

  for(int is=0; is<kns; is++) {
    SuperkSample_t s = kSKSamples[is];
    
    for(int ir=0; ir<kNErecoBins[is]; ir++) {
      if (s==kSKsFCFV1RingMuLike)
        fPredictErecoSpectrumArrayTot1Rmu[ir] = 0;
      else if(s==kSKsFCFV1RingELike)
        fPredictErecoSpectrumArrayTot1Re[ir] = 0;
      else if (s==kSKsRHCFCFV1RingMuLike)
        fPredictRHCErecoSpectrumArrayTot1Rmu[ir] = 0;
      else if(s==kSKsRHCFCFV1RingELike)
        fPredictRHCErecoSpectrumArrayTot1Re[ir] = 0;
      for(int im=0; im<knm; im++) {
        if (s==kSKsFCFV1RingMuLike)
          fPredictErecoSpectrumArray1Rmu[im][ir] = 0;
        else if(s==kSKsFCFV1RingELike)
          fPredictErecoSpectrumArray1Re[im][ir] = 0;
        else if (s==kSKsRHCFCFV1RingMuLike)
          fPredictRHCErecoSpectrumArray1Rmu[im][ir] = 0;
        else if(s==kSKsRHCFCFV1RingELike)
          fPredictRHCErecoSpectrumArray1Re[im][ir] = 0;
      }//im
    }//ir
  }//is
  cout << "done" << endl;
}
//________________________________________________________________________________________
void SuperkExpectation::Book(void)
{
  return;
}
//________________________________________________________________________________________
void SuperkExpectation::CleanUp(void)
{
  if (fSuperkDetectorFSI ) delete fSuperkDetectorFSI;
  if (fSuperkMC       ) delete fSuperkMC;
  if (fNuFlux         ) delete fNuFlux; 
  if (fNuXSec         ) delete fNuXSec;
  if (fBANFF          ) delete fBANFF;
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
  fSuperkDetectorFSI = new SuperkDetectorFSI; 
  fSuperkMC       = new SuperkMC;

  // Configure inter-dependencies.
  fNuFlux   -> TellMeAboutTheBANFFfit   (fBANFF);
  fNuFlux   -> TellMeAboutTheRHCBANFFfit(fBANFF);
  fNuXSec   -> TellMeAboutTheBANFFfit   (fBANFF);
  fNuXSec   -> TellMeAboutTheRHCBANFFfit(fBANFF);
  fSuperkMC -> TellMeAboutCrossSections (fNuXSec);
  fSuperkMC -> TellMeAboutTheFlux       (fNuFlux);
  fSuperkMC -> TellMeAboutTheDetectorFSI(fSuperkDetectorFSI);

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
  if ( ! this->GetBANFF            () -> Read (xml, root_node) ) return false;
  if ( ! this->GetBANFF            () -> ReadRHC (xml, root_node) ) return false;
  if ( ! this->GetNuFlux           () -> Read (xml, root_node) ) return false;
  if ( ! this->GetNuFlux           () -> ReadRHC (xml, root_node) ) return false;
  if ( ! this->GetNuXSec           () -> Read (xml, root_node) ) return false;
  if ( ! this->GetSuperkDetectorFSI() -> Read (xml, root_node) ) return false;
  if ( ! this->GetSuperkDetectorFSI() -> ReadRHC (xml, root_node) ) return false;
  if ( ! this->GetSuperkMC         () -> Read (xml, root_node) ) return false;

  fIsConfigured = true;

  // Release memory
  xml->FreeDoc(xmldoc);
  delete xml;

  cout << "Exp >> Done configuring Superk spectrum predictor... \n\n" << endl;

  return true;
}
//________________________________________________________________________________________
bool SuperkExpectation::Generate(const PhysHypothesis & ph)
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
  bool applied_negevents_condition = fSuperkMC->UpdateErecoSpectra();

  // copy histograms
  for(int is = 0; is < kns; is++) {
    SuperkSample_t s = kSKSamples[is];
    if(s == kSKsUndefined) continue;
    const int nbins = kNErecoBins[is];
    for(int ibin = 1; ibin <= nbins; ibin++) {
      double ns = fSuperkMC->CurrErecoSpectrumArray(s)[ibin-1];
      this->CurrPredictedErecoSpectrumArray(s)[ibin-1] = ns;
      for(int im = 0; im < knm; im++) {
        SuperkMode_t m = kSKReactionModes[is][im];
        if(m == kSKmUndefined) continue;
        double nsm = fSuperkMC->CurrErecoSpectrumArray(s,m)[ibin-1];
        this->CurrPredictedErecoSpectrumArrayMode(s,m)[ibin-1] = nsm;
      }//m
    }//ibin
  }//s
  return applied_negevents_condition;
}
//________________________________________________________________________________________
bool SuperkExpectation::Generate(const FitParamList & fp, bool all_banff_fit, bool all_1Re_eff)
{
// SuperkExpectation interface method for fitters.
// Calculate dN^{exp}_{SuperK}(E_{reco}) / dE_{reco} for the input E_{reco} and 
// the input physics + nuisance parameter list.
//

  // set values of nuisance parameters 
  // (if any nuisance parameter is included in the fit)
  this->SetNuisanceParams(fp, all_banff_fit, all_1Re_eff);

  // extract the values of params driving the oscillation, decoherence or decay model
  PhysHypothesis ph(fp);

  // call main Generate() method to generate SuperK predictions for the input
  // physics hypothesis using the (potentially) tweaked p.d.fs
  bool applied_negevents_condition = this->Generate(ph);
  return applied_negevents_condition;
}
//________________________________________________________________________________________
TH1D* SuperkExpectation::CurrPredictedErecoSpectrum(t2kana::conv::SuperkSample_t s) const 
{
  const double * arr = CurrPredictedErecoSpectrumArray(s);
  if (!arr) return 0;
  
  string name  = Form("curr_ereco_%s"                       ,util::AsStringCompact(s).c_str() );
  string title = Form("Current Predicted Ereco Spectrum: %s",util::AsString(s).c_str() );
  
  int nbins = util::NumErecoBins(s);
  TH1D * h = new TH1D(name.c_str(), title.c_str(), nbins, util::ErecoBinEdges(s) );
  util::FormatErecoHist(h);
  for (int i = 0 ; i < nbins ; i++ ) {
    h->SetBinContent(i+1,arr[i]);
  }
  return h;
}
//________________________________________________________________________________________
TH1D* SuperkExpectation::CurrPredictedErecoSpectrum(
  t2kana::conv::SuperkSample_t s, t2kana::conv::SuperkMode_t m) const 
{ 
  const double * arr = CurrPredictedErecoSpectrumArrayMode(s,m);
  if (!arr) return 0;
  
  int nbins = util::NumErecoBins(s);
  string name  = Form("curr_ereco_%s_%s"                       ,util::AsStringCompact(s).c_str(), util::AsStringCompact(m).c_str());
  string title = Form("Current Predicted Ereco Spectrum: %s %s",util::AsString(s).c_str()       , util::AsString(m).c_str()       );
  TH1D * h = new TH1D(name.c_str(), title.c_str(), nbins, util::ErecoBinEdges(s) );
  util::FormatErecoHist(h);
  for (int i = 0 ; i < nbins ; i++ ) {
    h->SetBinContent(i+1,arr[i]);
  }
  return h;
}
//________________________________________________________________________________________
double* SuperkExpectation::CurrPredictedErecoSpectrumArray(t2kana::conv::SuperkSample_t s) 
{
  int sidx = util::ArrIdx::Instance()->SuperkSampleIdx(s);
  if(sidx == -1) {
    return 0;
  }

  if(s == kSKsRHCFCFV1RingMuLike)
    return fPredictRHCErecoSpectrumArrayTot1Rmu;
  else if(s == kSKsRHCFCFV1RingELike)
    return fPredictRHCErecoSpectrumArrayTot1Re;
  else if(s == kSKsFCFV1RingELike)
    return fPredictErecoSpectrumArrayTot1Re;
  else
    return fPredictErecoSpectrumArrayTot1Rmu;

}
//________________________________________________________________________________________
const double* SuperkExpectation::CurrPredictedErecoSpectrumArray(t2kana::conv::SuperkSample_t s)  const
{
  int sidx = util::ArrIdx::Instance()->SuperkSampleIdx(s);
  if(sidx == -1) {
    return 0;
  }

  if(s == kSKsRHCFCFV1RingMuLike)
    return fPredictRHCErecoSpectrumArrayTot1Rmu;
  else if(s == kSKsRHCFCFV1RingELike)
    return fPredictRHCErecoSpectrumArrayTot1Re;
  else if(s == kSKsFCFV1RingELike)
    return fPredictErecoSpectrumArrayTot1Re;
  else
    return fPredictErecoSpectrumArrayTot1Rmu;

}
//________________________________________________________________________________________
double* SuperkExpectation::CurrPredictedErecoSpectrumArrayMode(
  t2kana::conv::SuperkSample_t s, t2kana::conv::SuperkMode_t m)
{ 
  int sidx = util::ArrIdx::Instance()->SuperkSampleIdx(s);
  if(sidx == -1) {
    return 0;
  }
  int midx = util::ArrIdx::Instance()->SuperkModeIdx(s,m);
  if(midx == -1) {
    return 0;
  }

  if(s == kSKsRHCFCFV1RingMuLike)
    return fPredictRHCErecoSpectrumArray1Rmu[midx];
  else if(s == kSKsRHCFCFV1RingELike)
    return fPredictRHCErecoSpectrumArray1Re[midx];
  else if(s == kSKsFCFV1RingELike)
    return fPredictErecoSpectrumArray1Re[midx];
  else
    return fPredictErecoSpectrumArray1Rmu[midx];

}
//________________________________________________________________________________________
const double* SuperkExpectation::CurrPredictedErecoSpectrumArrayMode(
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

  if(s == kSKsRHCFCFV1RingMuLike)
    return fPredictRHCErecoSpectrumArray1Rmu[midx];
  else if(s == kSKsRHCFCFV1RingELike)
    return fPredictRHCErecoSpectrumArray1Re[midx];
  else if(s == kSKsFCFV1RingELike)
    return fPredictErecoSpectrumArray1Re[midx];
  else
    return fPredictErecoSpectrumArray1Rmu[midx];

}
//________________________________________________________________________________________
double SuperkExpectation::CurrPredictedErecoSpectrumArray(t2kana::conv::SuperkSample_t s, int i) const
{
  const double * arr = CurrPredictedErecoSpectrumArray(s);
  assert(arr);
  return arr[i];
}
//________________________________________________________________________________________
void SuperkExpectation::PrintSystematicVariations(ostream & stream) const
{
  cout << "Exp >> SuperK PDF tweaks: " << endl;
  fBANFF          -> PrintSystematicVariations(stream);
  fBANFF          -> PrintRHCSystematicVariations(stream);
  fNuXSec         -> PrintSystematicVariations(stream);
  fNuXSec         -> PrintRHCSystematicVariations(stream);
  fSuperkDetectorFSI -> PrintSystematicVariations(stream);
  fSuperkDetectorFSI -> PrintRHCSystematicVariations(stream);

}
//________________________________________________________________________________________
bool SuperkExpectation::Write(TObjArray * array, const char * pdf_tag)
{
  if(!array) return false;

  util::SetDefaultStyle();

  const char * append_begin = Form("%s_", pdf_tag);

  for(int is=0; is<kns; is++) {
    SuperkSample_t s = kSKSamples[is];
    if(s == kSKsUndefined) continue;

    TH1D * spectrum = 0;
    TH1D * spectrum_mode[knm] = {0};

     // copy & format Ereco spectra
     spectrum = this->CurrPredictedErecoSpectrum(s);
     spectrum->SetDirectory (0);
     for(int ir = 1; ir <= spectrum->GetNbinsX(); ir++) {
       spectrum->SetBinContent(ir, this->CurrPredictedErecoSpectrumArray(s)[ir-1]);
     }
     for(int im = 0; im < knm; im++) {
        SuperkMode_t m = kSKReactionModes[is][im];
        if(m == kSKmUndefined) continue;
        spectrum_mode[im] = this->CurrPredictedErecoSpectrum(s,m);
        spectrum_mode[im]->SetDirectory (0);
	for(int ir = 1; ir <= spectrum->GetNbinsX(); ir++) {
	  spectrum_mode[im]->SetBinContent(ir, this->CurrPredictedErecoSpectrumArrayMode(s,m)[ir-1]);
	}
     }//imode

     // add spectra to the results array
     spectrum->SetName(util::AsStringCompact(s,kSKmUndefined,append_begin).c_str());
     array->Add(spectrum);
     for(int im = 0; im < knm; im++) {
        SuperkMode_t m = kSKReactionModes[is][im];
        if(m == kSKmUndefined) continue;
        spectrum_mode[im]->SetName(util::AsStringCompact(s,m,append_begin).c_str());
        array->Add(spectrum_mode[im]);
     }//im
  }//is

  return true;
}
//________________________________________________________________________________________
bool SuperkExpectation::Write(
  const char * filename, const char * pdf_tag, const char * option_file)
{
  TObjArray * array = new TObjArray();
  this->Write(array,pdf_tag);

  TFile f(filename, option_file);
  array->Write();
  f.Close();

  return true;
}
//________________________________________________________________________________________
void SuperkExpectation::ResetSystematicVariations(void)
{
#ifdef __show_syst_reset_mesg__
  cout << "Exp >> Reseting all systematic variations." << endl;
#endif

  //
  //FHC
  //

  // SK det + FSI
  this->GetSuperkDetectorFSI()->TweakEnergyScale(0.);
  this->GetSuperkDetectorFSI()->Tweak(true,0);

  // BANFF
  this->GetBANFF()->Tweak(true,0);

  // X Section
  this->GetNuXSec()->Tweak(kNuisXSecFermiMom      , 0.);
  this->GetNuXSec()->Tweak(kNuisXSecSpecFunc      , 0.);
  this->GetNuXSec()->Tweak(kNuisXSecWshape        , 0.);
  this->GetNuXSec()->Tweak(kNuisXSecShapeCCoth    , 0.);
  this->GetNuXSec()->Tweak(kNuisXSecPiLessDelta   , 0.);
#ifndef __NUEONLY__
  this->GetNuXSec()->Tweak(kNuisXSecBindE         , 0.);
#endif
  this->GetNuXSec()->Tweak(kNuisXSecNormCCcoh     , 0.);
#ifndef __NUEONLY__
  this->GetNuXSec()->Tweak(kNuisXSecNormNC1piC    , 0.);
#endif
#ifndef __NUMUONLY__
  this->GetNuXSec()->Tweak(kNuisXSecNormNCcoh     , 0.);
#endif
  this->GetNuXSec()->Tweak(kNuisXSecNormNCoth     , 0.);
  this->GetNuXSec()->Tweak(kNuisXSecNormNuEToNuMu , 0.);
#ifndef __NUEONLY__
  this->GetNuXSec()->Tweak(kNuisXSecNormAntiNuToNu, 0.);
#endif

  //
  //RHC
  //
  this->GetNuXSec()->Tweak(kNuisRHCToFHC, 0.);

  // SK det + FSI
  this->GetSuperkDetectorFSI()->TweakRHCEnergyScale(0.);
  this->GetSuperkDetectorFSI()->TweakRHC(true,0);

  // BANFF
  this->GetBANFF()->TweakRHC(true,0);

  // X Section
  this->GetNuXSec()->Tweak(kRHCNuisXSecFermiMom      , 0.);
  this->GetNuXSec()->Tweak(kRHCNuisXSecSpecFunc      , 0.);
  this->GetNuXSec()->Tweak(kRHCNuisXSecWshape        , 0.);
  this->GetNuXSec()->Tweak(kRHCNuisXSecShapeCCoth    , 0.);
  this->GetNuXSec()->Tweak(kRHCNuisXSecPiLessDelta   , 0.);
#ifndef __NUEONLY__
  this->GetNuXSec()->Tweak(kRHCNuisXSecBindE         , 0.);
#endif
  this->GetNuXSec()->Tweak(kRHCNuisXSecNormCCcoh     , 0.);
#ifndef __NUEONLY__
  this->GetNuXSec()->Tweak(kRHCNuisXSecNormNC1piC    , 0.);
#endif
#ifndef __NUMUONLY__
  this->GetNuXSec()->Tweak(kRHCNuisXSecNormNCcoh     , 0.);
#endif
  this->GetNuXSec()->Tweak(kRHCNuisXSecNormNCoth     , 0.);
  this->GetNuXSec()->Tweak(kRHCNuisXSecNormNuEToNuMu , 0.);
#ifndef __NUEONLY__
  this->GetNuXSec()->Tweak(kRHCNuisXSecNormAntiNuToNu, 0.);
#endif
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
//   - sk_1Re_eff                           :
//   - sk_1Re_eff_...                       :
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
//      * banff_mecnorm                     : tweak only BANFF MEC norm parameter
//      * banff_pdd                         : tweak only BANFF PDD double counting parameter
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
//    - xsec_bindE
//    - xsec_norm_cccoh
//    - xsec_norm_nc1piC
//    - xsec_norm_ncoth
//    - xsec_norm_nuetonumu
//    - xsec_norm_antinutonu
//


  //
  // SK detector systematics
  //

  if(systematic == "sk_energy_scale") {
    this->GetSuperkDetectorFSI()->TweakEnergyScale(twkdial);
  }
  else
  if(systematic == "sk_energy_scale_rhc") {
    this->GetSuperkDetectorFSI()->TweakRHCEnergyScale(twkdial);
  }

  //
  // SKDET + FSI
  //

  // tweak all FSI parameters taking into account correlations
  else
  if(systematic == "skdetfsi") {
    this->GetSuperkDetectorFSI()->Tweak(); //< Cholesky method
  }
  else
  if(systematic == "skdetfsi_rhc") {
    this->GetSuperkDetectorFSI()->TweakRHC(); //< Cholesky method
  }
  // tweak an individual FSI parameter (ignores correlations, debug)
  else
  if(systematic.find("skdetfsi_") != string::npos) {
    bool is_rhc = systematic.find("rhc") != string::npos;
    bool found = false;
    int nskdetfsi = 0;
    if(is_rhc)
      nskdetfsi = conv::kNumOfRHCSKDetFSIParams;
    else
      nskdetfsi = conv::kNumOfSKDetFSIParams;
    TVectorD twkv(nskdetfsi);
    for (int ip = 0; ip < nskdetfsi; ip++) {
      twkv[ip] = 1;
    }
    for (int ip = 0; ip < nskdetfsi; ip++) {
      if(is_rhc) {
	if(systematic == Form("skdetfsi_%03i_rhc",ip)) {
	  twkv[ip]=twkdial;
	  this->GetSuperkDetectorFSI()->TweakRHC(true,&twkv);
	  found = true;
	  break;
	}
      }//is_rhc
      else {
	if(systematic == Form("skdetfsi_%03i",ip)) {
	  twkv[ip]=twkdial;
	  this->GetSuperkDetectorFSI()->Tweak(true,&twkv);
	  found = true;
	  break;
	}
      }//!is_rhc
    }//ip
     if(!found) {
       cerr << "Exp >> Unknown systematic parameter (" << systematic << ")" << endl;
       exit(1);
     }   
  }//skdetfsi
    
  //
  // BANFF & Uncorrelated BANFF
  //

  // tweak all BANFF parameters taking into account correlations
  else
  if(systematic == "banff") {
    this->GetBANFF()->Tweak(false,0); //< Cholesky method
  }
  else
  if(systematic == "banff_rhc") {
    this->GetBANFF()->TweakRHC(false,0); //< Cholesky method
  }

  // tweak an individual BANFF parameter (ignores correlations, debug)
  else
  if(systematic.find("banff_") != string::npos) {
    bool is_rhc = systematic.find("rhc") != string::npos;
    if(is_rhc) {
      const int nbanff = conv::kNumOfRHCParamsInBANFFfit;
      TVectorD twkv(this->GetBANFF()->RHCBestFitValues());
      if(!best_banff) {
	for (int ib = 0; ib < nbanff; ib++) { twkv[ib] = 1; }
      }
      PhysParamId_t p = util::PhysParamIdFromString(systematic);
      int ip = util::ArrIdx::Instance()->RHCBANFFSystIdx(p);
      if(ip < 0) {
	cerr << "Exp >> Unknown systematic parameter (" << systematic << ")" << endl;
	exit(1);
      }
      twkv[ip] = twkdial;
#ifdef __show_syst_tweak_debug__
      cout << "Exp >> Tweaking systematic string: " << systematic
	   << ", NuisParam: " << util::AsString(p)
	   << endl;
      util::PrintCompact(twkv);
#endif
      this->GetBANFF()->TweakRHC(true,&twkv);
    }//is_rhc
    else {
      const int nbanff = conv::kNumOfParamsInBANFFfit;
      TVectorD twkv(this->GetBANFF()->BestFitValues());
      if(!best_banff) {
	for (int ib = 0; ib < nbanff; ib++) { twkv[ib] = 1; }
      }
      PhysParamId_t p = util::PhysParamIdFromString(systematic);
      int ip = util::ArrIdx::Instance()->BANFFSystIdx(p);
      if(ip < 0) {
	cerr << "Exp >> Unknown systematic parameter (" << systematic << ")" << endl;
	exit(1);
      }
      twkv[ip] = twkdial;
#ifdef __show_syst_tweak_debug__
      cout << "Exp >> Tweaking systematic string: " << systematic
	   << ", NuisParam: " << util::AsString(p)
	   << endl;
      util::PrintCompact(twkv);
#endif
      this->GetBANFF()->Tweak(true,&twkv);
    }//!is_rhc
  }//banff

  //
  // XSec
  //
  
  else if(systematic == "xsec_pf")
    this->GetNuXSec()->Tweak(kNuisXSecFermiMom      , twkdial);
  else if(systematic == "xsec_sf")
    this->GetNuXSec()->Tweak(kNuisXSecSpecFunc      , twkdial);
  else if(systematic == "xsec_wshape")
    this->GetNuXSec()->Tweak(kNuisXSecWshape        , twkdial);
  else if(systematic == "xsec_shape_ccoth")
    this->GetNuXSec()->Tweak(kNuisXSecShapeCCoth    , twkdial);
  else if(systematic == "xsec_pilessdelta")
    this->GetNuXSec()->Tweak(kNuisXSecPiLessDelta   , twkdial);
#ifndef __NUEONLY__
  else if(systematic == "xsec_bindE")
    this->GetNuXSec()->Tweak(kNuisXSecBindE         , twkdial);
#endif
  else if(systematic == "xsec_norm_cccoh")
    this->GetNuXSec()->Tweak(kNuisXSecNormCCcoh     , twkdial);
#ifndef __NUMUONLY__
  else if(systematic == "xsec_norm_nccoh")
    this->GetNuXSec()->Tweak(kNuisXSecNormNCcoh     , twkdial);
#endif
#ifndef __NUEONLY__
  else if(systematic == "xsec_norm_nc1piC")
    this->GetNuXSec()->Tweak(kNuisXSecNormNC1piC    , twkdial);
#endif
  else if(systematic == "xsec_norm_ncoth")
    this->GetNuXSec()->Tweak(kNuisXSecNormNCoth     , twkdial);
  else if(systematic == "xsec_norm_nuetonumu")
    this->GetNuXSec()->Tweak(kNuisXSecNormNuEToNuMu , twkdial);
#ifndef __NUEONLY__
  else if(systematic == "xsec_norm_antinutonu")
    this->GetNuXSec()->Tweak(kNuisXSecNormAntiNuToNu, twkdial);
#endif

  else if(systematic == "xsec_pf_rhc")
    this->GetNuXSec()->Tweak(kRHCNuisXSecFermiMom      , twkdial);
  else if(systematic == "xsec_sf_rhc")
    this->GetNuXSec()->Tweak(kRHCNuisXSecSpecFunc      , twkdial);
  else if(systematic == "xsec_wshape_rhc")
    this->GetNuXSec()->Tweak(kRHCNuisXSecWshape        , twkdial);
  else if(systematic == "xsec_shape_ccoth_rhc")
    this->GetNuXSec()->Tweak(kRHCNuisXSecShapeCCoth    , twkdial);
  else if(systematic == "xsec_pilessdelta_rhc")
    this->GetNuXSec()->Tweak(kRHCNuisXSecPiLessDelta   , twkdial);
#ifndef __NUEONLY__
  else if(systematic == "xsec_bindE_rhc")
    this->GetNuXSec()->Tweak(kRHCNuisXSecBindE         , twkdial);
#endif
  else if(systematic == "xsec_norm_cccoh_rhc")
    this->GetNuXSec()->Tweak(kRHCNuisXSecNormCCcoh     , twkdial);
#ifndef __NUMUONLY__
  else if(systematic == "xsec_norm_nccoh_rhc")
    this->GetNuXSec()->Tweak(kRHCNuisXSecNormNCcoh     , twkdial);
#endif
#ifndef __NUEONLY__
  else if(systematic == "xsec_norm_nc1piC_rhc")
    this->GetNuXSec()->Tweak(kRHCNuisXSecNormNC1piC    , twkdial);
#endif
  else if(systematic == "xsec_norm_ncoth_rhc")
    this->GetNuXSec()->Tweak(kRHCNuisXSecNormNCoth     , twkdial);
  else if(systematic == "xsec_norm_nuetonumu_rhc")
    this->GetNuXSec()->Tweak(kRHCNuisXSecNormNuEToNuMu , twkdial);
#ifndef __NUEONLY__
  else if(systematic == "xsec_norm_antinutonu_rhc")
    this->GetNuXSec()->Tweak(kRHCNuisXSecNormAntiNuToNu, twkdial);
#endif

  else if(systematic == "rhcfhcnorm")
    this->GetNuXSec()->Tweak(kNuisRHCToFHC, twkdial);

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

#ifdef HKHACK
  this->GetNuXSec()->Tweak(kNuisXSecNormCCcoh      , Rnd::Instance()->RndToyMCsyst().Gaus(0.,1.));
  this->GetNuXSec()->Tweak(kNuisXSecNormNC1piC     , Rnd::Instance()->RndToyMCsyst().Gaus(0.,1.));
  this->GetNuXSec()->Tweak(kNuisXSecNormNCcoh      , Rnd::Instance()->RndToyMCsyst().Gaus(0.,1.));
  this->GetNuXSec()->Tweak(kNuisXSecNormNCoth      , Rnd::Instance()->RndToyMCsyst().Gaus(0.,1.));
#else

  this->TweakAllBANFFSystematics();
  this->TweakAllUncorrelatedXSecSystematics();
  this->TweakAllSKDetectorFSISystematics();
  this->GetNuXSec()->Tweak(kNuisRHCToFHC, Rnd::Instance()->RndToyMCsyst().Gaus(0.,1.));

#ifdef __FULLY_CORRELATED__
  this->GetNuXSec()->SyncRHCToFHCSystematics();
  this->GetBANFF() ->SyncRHCToFHCSystematics();
  this->GetSuperkDetectorFSI()->SyncRHCToFHCSystematics();
#else
  this->TweakAllRHCBANFFSystematics();
  this->TweakAllRHCUncorrelatedXSecSystematics();
  this->TweakAllRHCSKDetectorFSISystematics();
#endif

#endif//HKHACK
}
//________________________________________________________________________________________
void SuperkExpectation::TweakAllBANFFSystematics(void)
{
#ifdef __show_syst_tweak_mesg__
  cout << "Exp >> Generating random variations of all BANFF systematic params (correlated)" << endl;
#endif

  // BANFF
  this->GetBANFF()->Tweak(false,0); // <-- Cholesky method
}
//________________________________________________________________________________________
void SuperkExpectation::TweakAllSKDetectorFSISystematics(void)
{
#ifdef __show_syst_tweak_mesg__
  cout << "Exp >> Generating random variations of all SKDET + FSI systematic params." << endl;
#endif

  this->GetSuperkDetectorFSI()->TweakEnergyScale(Rnd::Instance()->RndToyMCsyst().Gaus(0.,1.));
  this->GetSuperkDetectorFSI()->Tweak(); // <-- Cholesky method
}
//________________________________________________________________________________________
void SuperkExpectation::TweakAllUncorrelatedXSecSystematics(void)
{
#ifdef __show_syst_tweak_mesg__
  cout << "Exp >> Generating random variations of all Uncorrelated XSec systematic params." << endl;
#endif

  // XSec
  this->GetNuXSec()->Tweak(kNuisXSecFermiMom       , Rnd::Instance()->RndToyMCsyst().Gaus(0.,1.));
  double tweak = 2;
  while(tweak > 1 || tweak < 0)
    tweak = Rnd::Instance()->RndToyMCsyst().Gaus(0.,1.);
  this->GetNuXSec()->Tweak(kNuisXSecSpecFunc       , tweak);
  //this->GetNuXSec()->Tweak(kNuisXSecSpecFunc       , Rnd::Instance()->RndToyMCsyst().Uniform(0.,1.));
  //  this->GetNuXSec()->Tweak(kNuisXSecSpecFunc       , Rnd::Instance()->RndToyMCsyst().Gaus(0.,1.));
  this->GetNuXSec()->Tweak(kNuisXSecWshape         , Rnd::Instance()->RndToyMCsyst().Gaus(0.,1.));
  this->GetNuXSec()->Tweak(kNuisXSecShapeCCoth     , Rnd::Instance()->RndToyMCsyst().Gaus(0.,1.));
  this->GetNuXSec()->Tweak(kNuisXSecPiLessDelta    , Rnd::Instance()->RndToyMCsyst().Gaus(0.,1.));
#ifndef __NUEONLY__
  this->GetNuXSec()->Tweak(kNuisXSecBindE          , Rnd::Instance()->RndToyMCsyst().Gaus(0.,1.));
#endif
  this->GetNuXSec()->Tweak(kNuisXSecNormCCcoh      , Rnd::Instance()->RndToyMCsyst().Gaus(0.,1.));
#ifndef __NUEONLY__
  this->GetNuXSec()->Tweak(kNuisXSecNormNC1piC     , Rnd::Instance()->RndToyMCsyst().Gaus(0.,1.));
#endif
#ifndef __NUMUONLY__
  this->GetNuXSec()->Tweak(kNuisXSecNormNCcoh      , Rnd::Instance()->RndToyMCsyst().Gaus(0.,1.));
#endif
  this->GetNuXSec()->Tweak(kNuisXSecNormNCoth      , Rnd::Instance()->RndToyMCsyst().Gaus(0.,1.));
  this->GetNuXSec()->Tweak(kNuisXSecNormNuEToNuMu  , Rnd::Instance()->RndToyMCsyst().Gaus(0.,1.));
#ifndef __NUEONLY__
  this->GetNuXSec()->Tweak(kNuisXSecNormAntiNuToNu , Rnd::Instance()->RndToyMCsyst().Gaus(0.,1.));
#endif
}
//________________________________________________________________________________________
PhysParamList SuperkExpectation::Tweaks(void) 
{
  PhysParamList pl;

  //
  // BANFF
  //
  const int nbanff = conv::kNumOfParamsInBANFFfit;
  for(int ip = 0; ip < nbanff; ip++) {
     PhysParamId_t p = PhysParamId_t ((int) kNuisBANFF00 + ip);
     double v = this->GetBANFF()->CurrentValue(p);
     pl.SetValue(p,v);
  }

  // 
  //  Uncorrelated xsec
  // 
  pl.SetValue( kNuisXSecFermiMom,        this->GetNuXSec()->CurrentValue(kNuisXSecFermiMom)       );
  pl.SetValue( kNuisXSecSpecFunc,        this->GetNuXSec()->CurrentValue(kNuisXSecSpecFunc)       );
  pl.SetValue( kNuisXSecWshape,          this->GetNuXSec()->CurrentValue(kNuisXSecWshape)         );
  pl.SetValue( kNuisXSecShapeCCoth,      this->GetNuXSec()->CurrentValue(kNuisXSecShapeCCoth)     );
  pl.SetValue( kNuisXSecPiLessDelta,     this->GetNuXSec()->CurrentValue(kNuisXSecPiLessDelta)    );
#ifndef __NUEONLY__
  pl.SetValue( kNuisXSecBindE,           this->GetNuXSec()->CurrentValue(kNuisXSecBindE)          );
#endif
  pl.SetValue( kNuisXSecNormCCcoh,       this->GetNuXSec()->CurrentValue(kNuisXSecNormCCcoh)      );
#ifndef __NUEONLY__
  pl.SetValue( kNuisXSecNormNC1piC,      this->GetNuXSec()->CurrentValue(kNuisXSecNormNC1piC)     );
#endif
#ifndef __NUMUONLY__
  pl.SetValue( kNuisXSecNormNCcoh,       this->GetNuXSec()->CurrentValue(kNuisXSecNormNCcoh)      );
#endif
  pl.SetValue( kNuisXSecNormNCoth,       this->GetNuXSec()->CurrentValue(kNuisXSecNormNCoth)      );
  pl.SetValue( kNuisXSecNormNuEToNuMu,   this->GetNuXSec()->CurrentValue(kNuisXSecNormNuEToNuMu)  );
#ifndef __NUEONLY__
  pl.SetValue( kNuisXSecNormAntiNuToNu,  this->GetNuXSec()->CurrentValue(kNuisXSecNormAntiNuToNu) );
#endif

  // 
  // SuperK detector + FSI systematics
  //
  
  pl.SetValue( kNuisSuperkEnergyScale,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSuperkEnergyScale) );
#if FSI_BINS != 3
  pl.SetValue( kNuisSKDetFSI000,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSKDetFSI000) );
  pl.SetValue( kNuisSKDetFSI001,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSKDetFSI001) );
  pl.SetValue( kNuisSKDetFSI002,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSKDetFSI002) );
  pl.SetValue( kNuisSKDetFSI003,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSKDetFSI003) );
  pl.SetValue( kNuisSKDetFSI004,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSKDetFSI004) );
  pl.SetValue( kNuisSKDetFSI005,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSKDetFSI005) );
#endif
#if FSI_BINS != 2
  pl.SetValue( kNuisSKDetFSI006,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSKDetFSI006) );
  pl.SetValue( kNuisSKDetFSI007,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSKDetFSI007) );
  pl.SetValue( kNuisSKDetFSI008,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSKDetFSI008) );
  pl.SetValue( kNuisSKDetFSI009,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSKDetFSI009) );
  pl.SetValue( kNuisSKDetFSI010,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSKDetFSI010) );
  pl.SetValue( kNuisSKDetFSI011,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSKDetFSI011) );
  pl.SetValue( kNuisSKDetFSI012,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSKDetFSI012) );
  pl.SetValue( kNuisSKDetFSI013,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSKDetFSI013) );
  pl.SetValue( kNuisSKDetFSI014,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSKDetFSI014) );
  pl.SetValue( kNuisSKDetFSI015,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSKDetFSI015) );
  pl.SetValue( kNuisSKDetFSI016,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSKDetFSI016) );
  pl.SetValue( kNuisSKDetFSI017,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSKDetFSI017) );
#ifdef SKCORRELATIONS
  pl.SetValue( kNuisSKDetFSI018,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSKDetFSI018) );
  pl.SetValue( kNuisSKDetFSI019,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSKDetFSI019) );
  pl.SetValue( kNuisSKDetFSI020,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSKDetFSI020) );
  pl.SetValue( kNuisSKDetFSI021,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSKDetFSI021) );
  pl.SetValue( kNuisSKDetFSI022,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSKDetFSI022) );
  pl.SetValue( kNuisSKDetFSI023,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSKDetFSI023) );
  pl.SetValue( kNuisSKDetFSI024,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSKDetFSI024) );
  pl.SetValue( kNuisSKDetFSI025,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSKDetFSI025) );
  pl.SetValue( kNuisSKDetFSI026,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSKDetFSI026) );
  pl.SetValue( kNuisSKDetFSI027,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSKDetFSI027) );
  pl.SetValue( kNuisSKDetFSI028,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSKDetFSI028) );
  pl.SetValue( kNuisSKDetFSI029,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSKDetFSI029) );
  pl.SetValue( kNuisSKDetFSI030,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSKDetFSI030) );
  pl.SetValue( kNuisSKDetFSI031,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSKDetFSI031) );
#endif
#endif // FSI_BINS=1

  return pl;
}
//________________________________________________________________________________________
void SuperkExpectation::SetNuisanceParams(const FitParamList & fp, bool /*all_banff_fit*/, bool /*all_1Re_eff*/)
{
  //this->ResetSystematicVariations();
  

  //below is what it used to do - check that every single systematic parameter exists in the fit,
  // and only then tweak
  //FitParamList has a 'SystMode', which tells you which parameters are being fit - use this instead!

  if(fp.GetSystMode() == 1) {
    //energy scale
    if(fp.Exists(kNuisSuperkEnergyScale)) {
      double twk = fp.CurrValue(kNuisSuperkEnergyScale);
     
      //Todo: fix this so that FHC and RHC can have different tweaks
      this->GetSuperkDetectorFSI()->TweakEnergyScale(twk);
      this->GetSuperkDetectorFSI()->TweakRHCEnergyScale(twk);
    }

    // use current values stored in fit param list
    // to tweak the SKDET + FSI object
#if FSI_BINS == 2
    const int nskdetfsi_1Rmu = conv::kNumOfSKDetFSIParams1Rmu;
    //    if (Opt::Instance()->gOscAnaTypeChoice == conv::kAnaNuMuDisappOnly){
    TVectorD twkv(nskdetfsi_1Rmu);
    for(int ip = 0; ip < nskdetfsi_1Rmu; ip++) {
      PhysParamId_t p = PhysParamId_t ((int) kNuisSKDetFSI000 + ip);
      twkv[ip] = fp.CurrValue(p);
    }

    //Todo: fix this so that FHC and RHC can have different tweaks
    //Also that FHC and RHC can have different numbers of parameters
    this->GetSuperkDetectorFSI()->TweakOnly1Rmu(true,&twkv);
    this->GetSuperkDetectorFSI()->TweakOnlyRHC1Rmu(true,&twkv);
    //}
#else
#if FSI_BINS == 3
    const int nskdetfsi_1Re = conv::kNumOfSKDetFSIParams1Re;
      //    else if (Opt::Instance()->gOscAnaTypeChoice == conv::kAnaNuEAppOnly){
      TVectorD twkv(nskdetfsi_1Re);
      for(int ip = 0; ip < nskdetfsi_1Re; ip++) {
	//	PhysParamId_t p = PhysParamId_t ((int) kNuisSKDetFSI000 + ip + nskdetfsi_1Rmu);
	PhysParamId_t p = PhysParamId_t ((int) kNuisSKDetFSI006 + ip);
	twkv[ip] = fp.CurrValue(p);
      }

      //Todo: fix this so that FHC and RHC can have different tweaks
      //Also that FHC and RHC can have different numbers of parameters
      this->GetSuperkDetectorFSI()->TweakOnly1Re(true,&twkv);
      this->GetSuperkDetectorFSI()->TweakOnlyRHC1Re(true,&twkv);
      //}
#else
      //    else {
      const int nskdetfsi = conv::kNumOfSKDetFSIParams;
      TVectorD twkv(nskdetfsi);
      for(int ip = 0; ip < nskdetfsi; ip++) {
	PhysParamId_t p = PhysParamId_t ((int) kNuisSKDetFSI000 + ip);
	twkv[ip] = fp.CurrValue(p);
      }

      //Todo: fix this so that FHC and RHC can have different tweaks
      //Also that FHC and RHC can have different numbers of parameters
      this->GetSuperkDetectorFSI()->Tweak(true,&twkv);    
      this->GetSuperkDetectorFSI()->TweakRHC(true,&twkv);
      //}
#endif
#endif
   
    // use current values stored in fit param list
    // to tweak the BANFF object
    const int nbanff = conv::kNumOfParamsInBANFFfit;
    TVectorD twkvBANFF(nbanff);
    for(int ip = 0; ip < nbanff; ip++) {
      PhysParamId_t p = PhysParamId_t ((int) kNuisBANFF00 + ip);
      twkvBANFF[ip] = fp.CurrValue(p);
    }

    //Todo: fix this so that FHC and RHC can have different tweaks
    //Also that FHC and RHC can have different numbers of parameters 
    this->GetBANFF()->Tweak(true,&twkvBANFF);
    this->GetBANFF()->TweakRHC(true,&twkvBANFF);

    // Cross-section systematics (other than BANFF ones)
    for(int ip = 0; ip < kNumOfXSecSyst; ip++) {
      PhysParamId_t p = kListOfXSecSyst[ip];
      if( fp.Exists(p) && util::IsXSec(p) ) {
        double twk = fp.CurrValue(p);
        this->GetNuXSec()->Tweak(p,twk);
      }
    }

  }//systmode == 1

  else if(fp.GetSystMode() == 2) {
    cerr << "selected systematics not implemented in SuperkExpectation::SetNuisanceParams(). Exiting..." << endl;
    exit(1);
  }
}
//________________________________________________________________________________________

