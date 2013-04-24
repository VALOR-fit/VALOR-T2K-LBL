//________________________________________________________________________________________
/*

 >> ChangeLog:

 For changes up to version 2010_v3.0.0 (last implementation of the 2010-era analysis)
 see code in the corresponding SVN branch.

 Changes since 2011_v0.9.0:

 @ 18/2/13 LE
 Added tweaking methods for new parameters for the 1Re sample

*/
//________________________________________________________________________________________

#define __show_syst_reset_mesg__
#define __show_syst_tweak_mesg__
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
#include "SuperkDetector.h"
#include "SuperkDetectorFSI.h"
#include "SuperkMC.h"
#include "NuFlux.h"
#include "NuXSec.h"
#include "BANFF.h"
#include "FSI.h"
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

  fSuperkDetector = 0;
  fSuperkDetectorFSI = 0;
  fSuperkMC       = 0;
  fNuFlux         = 0;
  fNuXSec         = 0;
  fFSI            = 0;
  fBANFF          = 0;

  for(int is=0; is<kns; is++) {
    SuperkSample_t s = kSKSamples[is];
    fPredictErecoSpectrumTot[is] = 0;
    for(int im=0; im<knm; im++) {
      fPredictErecoSpectrum[is][im] = 0;
    }//im
    for(int ir=0; ir<util::NumErecoBins(s); ir++) {
      if (s==kSKsFCFV1RingMuLike) fPredictErecoSpectrumArrayTot1Rmu[ir] = 0;
      else fPredictErecoSpectrumArrayTot1Re[ir] = 0;
    }//ir
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
  if (fSuperkDetectorFSI ) delete fSuperkDetectorFSI;
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
  fSuperkDetectorFSI = new SuperkDetectorFSI; 
  fSuperkMC       = new SuperkMC;

  // Configure inter-dependencies.
  fNuFlux   -> TellMeAboutTheBANFFfit   (fBANFF);
  fNuXSec   -> TellMeAboutTheBANFFfit   (fBANFF);
  fSuperkMC -> TellMeAboutCrossSections (fNuXSec);
  fSuperkMC -> TellMeAboutTheFlux       (fNuFlux);
  fSuperkMC -> TellMeAboutTheDetectorFSI(fSuperkDetectorFSI);
  fSuperkMC -> TellMeAboutTheDetector   (fSuperkDetector);
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
  if ( ! this->GetBANFF            () -> Read (xml, root_node) ) return false;
  if ( ! this->GetNuFlux           () -> Read (xml, root_node) ) return false;
  if ( ! this->GetNuXSec           () -> Read (xml, root_node) ) return false;
  if ( ! this->GetFSI              () -> Read (xml, root_node) ) return false;
  if ( ! this->GetSuperkDetectorFSI() -> Read (xml, root_node) ) return false;
  if ( ! this->GetSuperkDetector   () -> Read (xml, root_node) ) return false;
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
    const int nbins = util::NumErecoBins(s);
    for(int ibin = 1; ibin <= nbins; ibin++) {
      //      double ns = fSuperkMC->CurrErecoSpectrum(s)->GetBinContent(ibin);
      double ns = fSuperkMC->CurrErecoSpectrumArray(s)[ibin-1];
      this->CurrPredictedErecoSpectrumArray(s)[ibin-1] = ns;
#ifndef PRODUCTIONMODE
      this->CurrPredictedErecoSpectrum(s)->SetBinContent(ibin, ns);
      for(int im = 0; im < knm; im++) {
	SuperkMode_t m = kSKReactionModes[is][im];
	if(m == kSKmUndefined) continue;
	double nsm = fSuperkMC->CurrErecoSpectrum(s,m)->GetBinContent(ibin);
	this->CurrPredictedErecoSpectrum(s,m)->SetBinContent(ibin, nsm);
      }//m	
#endif
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
double* SuperkExpectation::CurrPredictedErecoSpectrumArray(t2kana::conv::SuperkSample_t s) 
{
  int sidx = util::ArrIdx::Instance()->SuperkSampleIdx(s);
  if(sidx == -1) {
    return 0;
  }
  if (s== kSKsFCFV1RingMuLike )return fPredictErecoSpectrumArrayTot1Rmu;  
  else return fPredictErecoSpectrumArrayTot1Re;
}
//________________________________________________________________________________________
double  SuperkExpectation::CurrPredictedErecoSpectrumArray(t2kana::conv::SuperkSample_t s, int i) 
{
  int sidx = util::ArrIdx::Instance()->SuperkSampleIdx(s);
  if(sidx == -1) {
    return 0;
  }
  if (s== kSKsFCFV1RingMuLike )return fPredictErecoSpectrumArrayTot1Rmu[i];  
  else return fPredictErecoSpectrumArrayTot1Re[i];
}
//________________________________________________________________________________________
void SuperkExpectation::PrintSystematicVariations(ostream & stream) const
{
  cout << "Exp >> SuperK PDF tweaks: " << endl;
  fBANFF          -> PrintSystematicVariations(stream);
  fNuXSec         -> PrintSystematicVariations(stream);
#ifdef SKFSITOGETHER
  fSuperkDetectorFSI -> PrintSystematicVariations(stream);
#endif
#ifndef SKFSITOGETHER
  fSuperkDetector -> PrintSystematicVariations(stream);
  fFSI            -> PrintSystematicVariations(stream);
#endif
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
     spectrum = new TH1D(*this->CurrPredictedErecoSpectrum(s));
     spectrum->SetDirectory (0);
     for(int im = 0; im < knm; im++) {
        SuperkMode_t m = kSKReactionModes[is][im];
        if(m == kSKmUndefined) continue;
        spectrum_mode[im] = new TH1D(*this->CurrPredictedErecoSpectrum(s,m));
        spectrum_mode[im]->SetDirectory (0);
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

#ifndef SKFSITOGETHER
  // SuperK detector systematics
  this->GetSuperkDetector()->TweakEnergyScale                    (0.);
  this->GetSuperkDetector()->TweakEfficiency1RmuCCQEOther        (0.);
  this->GetSuperkDetector()->TweakEfficiency1RmuCCQERingCounting (true,0.,0.,0.); 
  this->GetSuperkDetector()->TweakEfficiency1RmuCCnonQE          (0.);
  this->GetSuperkDetector()->TweakEfficiency1RmuNC               (0.);     
  this->GetSuperkDetector()->TweakEfficiency1RmuNuECC            (0.);

  this->GetSuperkDetector()->TweakEfficiencies1Re                (true,0);

  // FSI
  this->GetFSI()->Tweak(true,0);
#endif
#ifdef SKFSITOGETHER
  this->GetSuperkDetectorFSI()->TweakEnergyScale(0.);
  this->GetSuperkDetectorFSI()->Tweak(true,0);
#endif

  // BANFF
  this->GetBANFF()->Tweak(true,0);


  // X Section
  this->GetNuXSec()->TweakFermiMom      (0.);
  this->GetNuXSec()->TweakSpecFunc      (0.);
  this->GetNuXSec()->TweakWshape        (0.);
  this->GetNuXSec()->TweakShapeCCoth    (0.);
  this->GetNuXSec()->TweakPiLessDelta   (0.);
  this->GetNuXSec()->TweakBindE         (0.);
  this->GetNuXSec()->TweakNormCCcoh     (0.);
  this->GetNuXSec()->TweakNormNCcoh     (0.);
  this->GetNuXSec()->TweakNormNC1piC    (0.);
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
     this->GetSuperkDetector()->TweakEnergyScale(twkdial);
  }

#ifndef SKFSITOGETHER
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
  // tweak 1Re efficiency generating random
  // efficiency tweaks and taking into account bin-to-bin correlations
  else
  if(systematic == "sk_1Re_eff") {
    this->GetSuperkDetector()->TweakEfficiencies1Re();
  }
  else
  if(systematic.find("sk_1Re_eff_") != string::npos) {
    bool found = false;
    const int neff1Re = conv::kNumOfSKDetParams1Re;
    TVectorD twkv(neff1Re);
    for (int ip = 0; ip < neff1Re; ip++) {
      twkv[ip] = 1;
    }
    for (int ip = 0; ip < neff1Re; ip++) {
      if(systematic == Form("sk_1Re_eff_%03i",ip)) {
	twkv[ip]=twkdial;
	this->GetSuperkDetector()->TweakEfficiencies1Re(true,&twkv);
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
#endif

#ifdef SKFSITOGETHER
  //
  // sKDET + FSI
  //

  // tweak all FSI parameters taking into account correlations
  else
  if(systematic == "skdetfsi") {
    this->GetSuperkDetectorFSI()->Tweak(); //< Cholesky method
  }
  // tweak an individual FSI parameter (ignores correlations, debug)
  else
  if(systematic.find("skdetfsi_") != string::npos) {
     bool found = false;
     const int nskdetfsi = conv::kNumOfSKDetFSIParams;
     TVectorD twkv(nskdetfsi);
     for (int ip = 0; ip < nskdetfsi; ip++) {
       twkv[ip] = 1;
     }
     for (int ip = 0; ip < nskdetfsi; ip++) {
       if(systematic == Form("skdetfsi_%03i",ip)) {
	 twkv[ip]=twkdial;
	 this->GetSuperkDetectorFSI()->Tweak(true,&twkv);
	 found = true;
	 break;
       }
     }//ib
     if(!found) {
       cerr << "Exp >> Unknown systematic parameter (" << systematic << ")" << endl;
       exit(1);
     }   
  }
  
  
#endif

  //
  // BANFF & Uncorrelated BANFF
  //

  // tweak all BANFF parameters taking into account correlations
  else
  if(systematic == "banff") {
    this->GetBANFF()->Tweak(false,0); //< Cholesky method
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
  else if(systematic == "xsec_bindE")
    this->GetNuXSec()->TweakBindE         (twkdial);
  else if(systematic == "xsec_norm_cccoh")
    this->GetNuXSec()->TweakNormCCcoh     (twkdial);
  else if(systematic == "xsec_norm_nccoh")
    this->GetNuXSec()->TweakNormNCcoh     (twkdial);
  else if(systematic == "xsec_norm_nc1piC")
    this->GetNuXSec()->TweakNormNC1piC    (twkdial);
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
  this->TweakAllUncorrelatedXSecSystematics();
#ifndef SKFSITOGETHER
  this->TweakAllFSISystematics();
  this->TweakAllSKDetectorSystematics();
#endif
#ifdef SKFSITOGETHER
  this->TweakAllSKDetectorFSISystematics();
#endif
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

  this->GetSuperkDetectorFSI()->TweakEnergyScale              (Rnd::Instance()->RndToyMCsyst().Gaus(0.,1.));
  this->GetSuperkDetectorFSI()->Tweak(); // <-- Cholesky method
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

  // XSec
  this->GetNuXSec()->TweakFermiMom      (Rnd::Instance()->RndToyMCsyst().Gaus(0.,1.));
  this->GetNuXSec()->TweakSpecFunc      (Rnd::Instance()->RndToyMCsyst().Uniform(0.,1.));
  this->GetNuXSec()->TweakWshape        (Rnd::Instance()->RndToyMCsyst().Gaus(0.,1.));
  this->GetNuXSec()->TweakShapeCCoth    (Rnd::Instance()->RndToyMCsyst().Gaus(0.,1.));
  this->GetNuXSec()->TweakPiLessDelta   (Rnd::Instance()->RndToyMCsyst().Gaus(0.,1.));
  this->GetNuXSec()->TweakBindE         (Rnd::Instance()->RndToyMCsyst().Gaus(0.,1.));
  this->GetNuXSec()->TweakNormCCcoh     (Rnd::Instance()->RndToyMCsyst().Gaus(0.,1.));
  this->GetNuXSec()->TweakNormNCcoh     (Rnd::Instance()->RndToyMCsyst().Gaus(0.,1.));
  this->GetNuXSec()->TweakNormNC1piC    (Rnd::Instance()->RndToyMCsyst().Gaus(0.,1.));
  this->GetNuXSec()->TweakNormNCoth     (Rnd::Instance()->RndToyMCsyst().Gaus(0.,1.));
  this->GetNuXSec()->TweakNormNuEToNuMu (Rnd::Instance()->RndToyMCsyst().Gaus(0.,1.));
  this->GetNuXSec()->TweakNormAntiNuToNu(Rnd::Instance()->RndToyMCsyst().Gaus(0.,1.));
}
//________________________________________________________________________________________
void SuperkExpectation::TweakAllSKDetectorSystematics(void)
{
#ifdef __show_syst_tweak_mesg__
  cout << "Exp >> Generating random variations of all SKDetector systematic params." << endl;
#endif


  // SuperK detector systematics
  this->GetSuperkDetector()->TweakEnergyScale              (Rnd::Instance()->RndToyMCsyst().Gaus(0.,1.));
  this->GetSuperkDetector()->TweakEfficiency1RmuCCQEOther  (Rnd::Instance()->RndToyMCsyst().Gaus(0.,1.));
  this->GetSuperkDetector()->TweakEfficiency1RmuCCQERingCounting  (); // <-- Cholesky method
  this->GetSuperkDetector()->TweakEfficiency1RmuCCnonQE    (Rnd::Instance()->RndToyMCsyst().Gaus(0.,1.));
  this->GetSuperkDetector()->TweakEfficiency1RmuNC         (Rnd::Instance()->RndToyMCsyst().Gaus(0.,1.));
  this->GetSuperkDetector()->TweakEfficiency1RmuNuECC      (Rnd::Instance()->RndToyMCsyst().Gaus(0.,1.));
  this->GetSuperkDetector()->TweakEfficiencies1Re          (); // <-- Cholesky method
}
//________________________________________________________________________________________
void SuperkExpectation::SetNuisanceParams(const FitParamList & fp, bool all_banff_fit, bool all_1Re_eff)
{
   //this->ResetSystematicVariations();
  
#ifndef SKFSITOGETHER
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
   if(fp.Exists(kNuisSuperkEnergyScale)) {
     double twk = fp.CurrValue(kNuisSuperkEnergyScale);
     this->GetSuperkDetector()->TweakEnergyScale(twk);
   }

   const int neff1Re = conv::kNumOfSKDetParams1Re;
   TVectorD twkv1Re(neff1Re); 
   // 1Re eff params included in the fit param list?                                                                                                                             
   /*                       
   bool incl_eff1Re = true;
   for(int ip = 0; ip < neff1Re; ip++) {
     PhysParamId_t p = PhysParamId_t ((int) kNuisSuperkEff1ReNueSigCC0 + ip);
     incl_eff1Re = incl_eff1Re && fp.Exists(p);
     if(!incl_eff1Re) break;
   }
   // If 1Re eff params included, use current values stored in fit param list                                                                                                                                
   // to tweak the 1Re eff object 
   if(incl_eff1Re) {
     TVectorD twkv(neff1Re);
     for(int ip = 0; ip < neff1Re; ip++) {
       PhysParamId_t p = PhysParamId_t ((int) kNuisSuperkEff1ReNueSigCC0 + ip);
       twkv[ip] = fp.CurrValue(p);
     }
     this->GetSuperkDetector()->TweakEfficiencies1Re(true,&twkv);
   }
   */
   if(all_1Re_eff){
     for(int ip = 0; ip < neff1Re; ip++) {
       PhysParamId_t p = PhysParamId_t ((int) kNuisSuperkEff1ReNueSigCC0 + ip);
       bool incl_1Re_eff = fp.Exists(p);
       if(incl_1Re_eff)
         twkv1Re[ip] = fp.CurrValue(p);
       else
         twkv1Re[ip] = this->GetSuperkDetector()->CurrentValue(p);
     }
   }
   else{
     const TMatrixD & cov = this->GetSuperkDetector()->Eff1ReCovarRelErr();
     // get the value of things we're actually fitting 
     PhysParamId_t p_signue2 = kNuisSuperkEff1ReNueSigCC2;
     bool incl_signue2 = fp.Exists(p_signue2);
     if(incl_signue2)
       twkv1Re[2] = fp.CurrValue(p_signue2);
     else
       twkv1Re[2] = this->GetSuperkDetector()->CurrentValue(p_signue2);
     double sigma_signue2 = (twkv1Re[2] - 1.0) / cov[2][2];

     PhysParamId_t p_nue2 = kNuisSuperkEff1ReNueCC2;
     bool incl_nue2 = fp.Exists(p_nue2);
     if(incl_nue2)
       twkv1Re[8] = fp.CurrValue(p_nue2);
     else
       twkv1Re[8] = this->GetSuperkDetector()->CurrentValue(p_nue2);
     double sigma_nue2 = (twkv1Re[8] - 1.0) / cov[8][8];

     PhysParamId_t p_nc2 = kNuisSuperkEff1ReNC2;
     bool incl_nc2 = fp.Exists(p_nc2);
     if(incl_nc2)
       twkv1Re[11] = fp.CurrValue(p_nc2);
     else
       twkv1Re[11] = this->GetSuperkDetector()->CurrentValue(p_nc2);
     double sigma_nc2 = (twkv1Re[11] - 1.0) / cov[11][11];

     for(int ip = 0; ip < neff1Re; ip++) {
       PhysParamId_t p = PhysParamId_t ((int) kNuisSuperkEff1ReNueSigCC0 + ip);
       if(p == kNuisSuperkEff1ReNueSigCC2 || p == kNuisSuperkEff1ReNueCC2 || p == kNuisSuperkEff1ReNC2)
         continue;
       twkv1Re[ip] = 1.0;
       twkv1Re[ip] += sigma_signue2  * cov[2] [ip];
       twkv1Re[ip] += sigma_nue2 * cov[8][ip];
       twkv1Re[ip] += sigma_nc2 * cov[11][ip];
     }
   }//!all_1Re_eff
   //And finally apply the tweak 
   this->GetSuperkDetector()->TweakEfficiencies1Re(true,&twkv1Re);
   
   // 
   // FSI systematics
   //
   
   const int nfsi = conv::kNumOfFSIParams;
   const int nfsi_1Rmu = conv::kNumOfFSIParams1Rmu;
   const int nfsi_1Re = conv::kNumOfFSIParams1Re;
   // FSI params included in the fit param list?
   bool incl_fsi = true;

   if (Opt::Instance()->gOscAnaTypeChoice == conv::kAnaNuMuDisappOnly){
     for(int ip = 0; ip < nfsi_1Rmu; ip++) {
       PhysParamId_t p = PhysParamId_t ((int) kNuisFSI000 + ip);
       incl_fsi = incl_fsi && fp.Exists(p);
       if(!incl_fsi) break;
     }
   }
   else if (Opt::Instance()->gOscAnaTypeChoice == conv::kAnaNuEAppOnly){
     for(int ip = nfsi_1Rmu; ip < nfsi; ip++) {
       PhysParamId_t p = PhysParamId_t ((int) kNuisFSI000 + ip);
       incl_fsi = incl_fsi && fp.Exists(p);
       if(!incl_fsi) break;
     }
   }
   else{
     for(int ip = 0; ip < nfsi; ip++) {
       PhysParamId_t p = PhysParamId_t ((int) kNuisFSI000 + ip);
       incl_fsi = incl_fsi && fp.Exists(p);
       if(!incl_fsi) break;
     }
   }
   // If FSI params included, use current values stored in fit param list
   // to tweak the FSI object
   if(incl_fsi) {
     if (Opt::Instance()->gOscAnaTypeChoice == conv::kAnaNuMuDisappOnly){
       TVectorD twkv(nfsi_1Rmu);
       for(int ip = 0; ip < nfsi_1Rmu; ip++) {
         PhysParamId_t p = PhysParamId_t ((int) kNuisFSI000 + ip);
         twkv[ip] = fp.CurrValue(p);
       }
       this->GetFSI()->TweakOnly1Rmu(true,&twkv);
     }
     else if (Opt::Instance()->gOscAnaTypeChoice == conv::kAnaNuEAppOnly){
       TVectorD twkv(nfsi_1Re);
       for(int ip = 0; ip < nfsi_1Re; ip++) {
         PhysParamId_t p = PhysParamId_t ((int) kNuisFSI000 + ip + nfsi_1Rmu);
         twkv[ip] = fp.CurrValue(p);
       }
       this->GetFSI()->TweakOnly1Re(true,&twkv);
     }
     else {
       TVectorD twkv(nfsi);
       for(int ip = 0; ip < nfsi; ip++) {
	 PhysParamId_t p = PhysParamId_t ((int) kNuisFSI000 + ip);
	 twkv[ip] = fp.CurrValue(p);
       }
       this->GetFSI()->Tweak(true,&twkv);    
     }
   }

#endif

#ifdef SKFSITOGETHER

   // 
   // SKDET + FSI systematics
   //

   if(fp.Exists(kNuisSuperkEnergyScale)) {
     double twk = fp.CurrValue(kNuisSuperkEnergyScale);
     this->GetSuperkDetectorFSI()->TweakEnergyScale(twk);
   }
   
   const int nskdetfsi = conv::kNumOfSKDetFSIParams;
   const int nskdetfsi_1Rmu = conv::kNumOfSKDetFSIParams1Rmu;
   const int nskdetfsi_1Re = conv::kNumOfSKDetFSIParams1Re;
   // SKDET + FSI params included in the fit param list?
   bool incl_skdetfsi = true;

   if (Opt::Instance()->gOscAnaTypeChoice == conv::kAnaNuMuDisappOnly){
     for(int ip = 0; ip < nskdetfsi_1Rmu; ip++) {
       PhysParamId_t p = PhysParamId_t ((int) kNuisSKDetFSI000 + ip);
       incl_skdetfsi = incl_skdetfsi && fp.Exists(p);
       if(!incl_skdetfsi) break;
     }
   }
   else if (Opt::Instance()->gOscAnaTypeChoice == conv::kAnaNuEAppOnly){
     for(int ip = nskdetfsi_1Rmu; ip < nskdetfsi; ip++) {
       PhysParamId_t p = PhysParamId_t ((int) kNuisSKDetFSI000 + ip);
       incl_skdetfsi = incl_skdetfsi && fp.Exists(p);
       if(!incl_skdetfsi) break;
     }
   }
   else{
     for(int ip = 0; ip < nskdetfsi; ip++) {
       PhysParamId_t p = PhysParamId_t ((int) kNuisSKDetFSI000 + ip);
       incl_skdetfsi = incl_skdetfsi && fp.Exists(p);
       if(!incl_skdetfsi) break;
     }
   }
   // If SKDET + FSI params included, use current values stored in fit param list
   // to tweak the SKDET + FSI object
   if(incl_skdetfsi) {
     if (Opt::Instance()->gOscAnaTypeChoice == conv::kAnaNuMuDisappOnly){
       TVectorD twkv(nskdetfsi_1Rmu);
       for(int ip = 0; ip < nskdetfsi_1Rmu; ip++) {
         PhysParamId_t p = PhysParamId_t ((int) kNuisSKDetFSI000 + ip);
         twkv[ip] = fp.CurrValue(p);
       }
       this->GetSuperkDetectorFSI()->TweakOnly1Rmu(true,&twkv);
     }
     else if (Opt::Instance()->gOscAnaTypeChoice == conv::kAnaNuEAppOnly){
       TVectorD twkv(nskdetfsi_1Re);
       for(int ip = 0; ip < nskdetfsi_1Re; ip++) {
         PhysParamId_t p = PhysParamId_t ((int) kNuisSKDetFSI000 + ip + nskdetfsi_1Rmu);
         twkv[ip] = fp.CurrValue(p);
       }
       this->GetSuperkDetectorFSI()->TweakOnly1Re(true,&twkv);
     }
     else {
       TVectorD twkv(nskdetfsi);
       for(int ip = 0; ip < nskdetfsi; ip++) {
	 PhysParamId_t p = PhysParamId_t ((int) kNuisSKDetFSI000 + ip);
	 twkv[ip] = fp.CurrValue(p);
       }
       this->GetSuperkDetectorFSI()->Tweak(true,&twkv);    
     }
   }

#endif

   // 
   // BANFF systematics
   //

   const int nbanff = conv::kNumOfParamsInBANFFfit;
   // BANFF params included in the fit param list?
   TVectorD twkv(nbanff);
   if(all_banff_fit) {
     for(int ip = 0; ip < nbanff; ip++) {
       PhysParamId_t p = PhysParamId_t ((int) kNuisBANFF00 + ip);
       bool incl_banff = fp.Exists(p);
       // If BANFF param included, use current value stored in fit param list
       // to tweak the BANFF object
       if(incl_banff)
	 twkv[ip] = fp.CurrValue(p);
       // Otherwise, use whatever is already stored in the BANFF object
       else
	 twkv[ip] = this->GetBANFF()->CurrentValue(p);
     }
   }//all_banff_fit
   else{
     const TMatrixD & cov = this->GetBANFF()->CovMatrix();
     // get the value of things we're actually fitting
     PhysParamId_t p_maqe = kNuisBANFFMAQE;
     int ip_maqe = util::ArrIdx::Instance()->BANFFSystIdx(p_maqe);
     bool incl_maqe = fp.Exists(p_maqe);
     if(incl_maqe)
       twkv[ip_maqe] = fp.CurrValue(p_maqe);
     else
       twkv[ip_maqe] = this->GetBANFF()->CurrentValue(p_maqe);
     double sigma_maqe = (twkv[ip_maqe] - this->GetBANFF()->BestFitValue(p_maqe)) / cov[ip_maqe][ip_maqe];

     PhysParamId_t p_ccqe0 = kNuisBANFFCCQEE0;
     int ip_ccqe0 = util::ArrIdx::Instance()->BANFFSystIdx(p_ccqe0);
     bool incl_ccqe0 = fp.Exists(p_ccqe0);
     if(incl_ccqe0)
       twkv[ip_ccqe0] = fp.CurrValue(p_ccqe0);
     else
       twkv[ip_ccqe0] = this->GetBANFF()->CurrentValue(p_ccqe0);
     double sigma_ccqe0 = (twkv[ip_ccqe0] - this->GetBANFF()->BestFitValue(p_ccqe0)) / cov[ip_ccqe0][ip_ccqe0];
     //new - joint
     double sigma_banff1 = 0;
     int ip_banff1 = -1;
     if (Opt::Instance()->gOscAnaTypeChoice == conv::kAnaCombined){
       PhysParamId_t p_banff1 = kNuisBANFF01;
       ip_banff1 = util::ArrIdx::Instance()->BANFFSystIdx(p_banff1);
       bool incl_banff1 = fp.Exists(p_banff1);
       if(incl_banff1)
	 twkv[ip_banff1] = fp.CurrValue(p_banff1);
       else
	 twkv[ip_banff1] = this->GetBANFF()->CurrentValue(p_banff1);
       sigma_banff1 = (twkv[ip_banff1] - this->GetBANFF()->BestFitValue(p_banff1)) / cov[ip_banff1][ip_banff1];
     }
     // now, tweak the other parameters according to the correlations between them and the fitted parameters
     for(int ip = 0; ip < nbanff; ip++) {
       PhysParamId_t p = PhysParamId_t ((int) kNuisBANFF00 + ip);
       if(p == kNuisBANFFMAQE || p == kNuisBANFFCCQEE0 )
	 continue;
       if (Opt::Instance()->gOscAnaTypeChoice == conv::kAnaCombined) {
	 if (p == kNuisBANFF01) continue;
       }
       twkv[ip] = this->GetBANFF()->BestFitValue(p);
       twkv[ip] += sigma_maqe  * cov[ip_maqe] [ip];
       twkv[ip] += sigma_ccqe0 * cov[ip_ccqe0][ip];
       if (Opt::Instance()->gOscAnaTypeChoice == conv::kAnaCombined) twkv[ip] += sigma_banff1 * cov[ip_banff1][ip];
     }
   }//!all_banff_fit
   //And finally apply the tweak
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
   if(fp.Exists(kNuisXSecBindE)) {
     double twk = fp.CurrValue(kNuisXSecBindE);
     this->GetNuXSec()->TweakBindE(twk);
   }
   if(fp.Exists(kNuisXSecNormCCcoh)) {
     double twk = fp.CurrValue(kNuisXSecNormCCcoh);
     this->GetNuXSec()->TweakNormCCcoh(twk);
   }
   if(fp.Exists(kNuisXSecNormNCcoh)) {
     double twk = fp.CurrValue(kNuisXSecNormNCcoh);
     this->GetNuXSec()->TweakNormNCcoh(twk);
   }
   if(fp.Exists(kNuisXSecNormNC1piC)) {
     double twk = fp.CurrValue(kNuisXSecNormNC1piC);
     this->GetNuXSec()->TweakNormNC1piC(twk);
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


}
//________________________________________________________________________________________

