//________________________________________________________________________________________
/*

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

#include "year1_numu_osc_analysis/SuperkExpectation.h"
#include "year1_numu_osc_analysis/Conventions.h"
#include "year1_numu_osc_analysis/Utilities.h"
#include "year1_numu_osc_analysis/HadErrCalc.h" ///<-- code by Kubo-san

using std::cout;
using std::cerr;
using std::endl;

using namespace t2kana;
using namespace t2kana::yr1numu;
using namespace t2kana::yr1numu::conv;

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
  const unsigned int nm = 6;  
  SuperkMode_t sk_modes [nm] = { 
     kSKmNuMuCCQE, kSKmNuMuCC1pi, kSKmNuMuCCoth, kSKmNC, kSKmNuMuBar, kSKmNuE 
  };

  // input physics params (params included in the fit / systematics params)
  double dmsq23          = fp.CurrValue(kPhypDm23sq);
  double sinsq2theta23   = fp.CurrValue(kPhypSinsq2Theta23);

  // correction factor from ND280
  double correction = fND280Inputs->CorrectionFactor();

  // update SK MC
  fSuperkInputs->UseOscillationParams(sinsq2theta23, dmsq23);
  fSuperkInputs->UpdateErecoSpectra();

  // update prediction
  for(int ibin = 1; ibin <= fPredictErecoSpectrumGrandTot->GetNbinsX(); ibin++) {
    double n = 0;
    for(unsigned int is = 0; is < ns; is++) {
       SuperkSample_t s = fSamples[is];
       double ns = correction * fSuperkInputs->CurrErecoSpectrum(s)->GetBinContent(ibin);
       n += ns;
       fPredictErecoSpectrumTot[s]->SetBinContent(ibin, ns);
       for(unsigned int im = 0; im < nm; im++) {
           SuperkMode_t m = sk_modes[im];
           double nsm = correction * fSuperkInputs->CurrErecoSpectrum(s,m)->GetBinContent(ibin);
           fPredictErecoSpectrum[s][m]->SetBinContent(ibin, nsm);       
       }//m	
    }//s
    double nall = correction * fSuperkInputs->CurrErecoSpectrum()->GetBinContent(ibin);
    fPredictErecoSpectrumGrandTot->SetBinContent(ibin, nall);
  }//ibin

  //
  double beam_norm = 1;
  if( this->GetBeamInputs()->UsingKuboSansCode() ) {
     // input physics params (params included in the fit / systematics params)
     this -> GetSuperkInputs() -> UpdateEtrueSpectra();
     this -> GetND280Inputs()  -> UpdateEtrueSpectra();
     TH1D*  hnumu1Rmu     = this -> GetSuperkInputs() -> Curr1RmuNuMuEtrue    ();
     TH1D*  hnumubar1Rmu  = this -> GetSuperkInputs() -> Curr1RmuNuMuBarEtrue ();
     TH1D*  hnue1Rm       = this -> GetSuperkInputs() -> Curr1RmuNuEEtrue     ();
     TH1D*  hnumuND       = this -> GetND280Inputs()  -> CurrNuMuEtrue        ();
     TH1D*  hnumubarND    = this -> GetND280Inputs()  -> CurrNuMuBarEtrue     ();
     TH1D*  hnueND        = this -> GetND280Inputs()  -> CurrNuEEtrue         ();

     beam_norm = this->GetBeamInputs()->SuperkFluxNorm(
               hnumu1Rmu, hnumubar1Rmu, hnue1Rm, hnumuND, hnumubarND, hnueND);
  } else {
     beam_norm = this->GetBeamInputs()->SuperkFluxNorm(0,0,0,0,0,0);
  }
  if(beam_norm != 1.) {
     fPredictErecoSpectrumGrandTot                             -> Scale (beam_norm);
     fPredictErecoSpectrumTot[kSKsFCFV1RingMuLike]             -> Scale (beam_norm);
     fPredictErecoSpectrum[kSKsFCFV1RingMuLike][kSKmNuMuCCQE]  -> Scale (beam_norm);
     fPredictErecoSpectrum[kSKsFCFV1RingMuLike][kSKmNuMuCC1pi] -> Scale (beam_norm);
     fPredictErecoSpectrum[kSKsFCFV1RingMuLike][kSKmNuMuCCoth] -> Scale (beam_norm);
     fPredictErecoSpectrum[kSKsFCFV1RingMuLike][kSKmNC]        -> Scale (beam_norm);
     fPredictErecoSpectrum[kSKsFCFV1RingMuLike][kSKmNuMuBar]   -> Scale (beam_norm);
     fPredictErecoSpectrum[kSKsFCFV1RingMuLike][kSKmNuE]       -> Scale (beam_norm);
  }

}
//________________________________________________________________________________________
TH1D* SuperkExpectation::CurrPredictedErecoSpectrum (void) const
{
  return fPredictErecoSpectrumGrandTot;  
}
//________________________________________________________________________________________
TH1D* SuperkExpectation::CurrPredictedErecoSpectrum (
  t2kana::yr1numu::conv::SuperkSample_t s) const 
{
  return fPredictErecoSpectrumTot[s];  
}
//________________________________________________________________________________________
TH1D* SuperkExpectation::CurrPredictedErecoSpectrum (
  t2kana::yr1numu::conv::SuperkSample_t s, t2kana::yr1numu::conv::SuperkMode_t m) const 
{ 
  return fPredictErecoSpectrum[s][m]; 
}
//________________________________________________________________________________________
void SuperkExpectation::ResetSystematicVariations(void)
{
  this->GetND280Inputs()->TweakEfficiency( kNDmNuMuCCQE,      0.);
  this->GetND280Inputs()->TweakEfficiency( kNDmNuMuCCRES,     0.);
  this->GetND280Inputs()->TweakEfficiency( kNDmNuMuCCCOH,     0.);
  this->GetND280Inputs()->TweakEfficiency( kNDmNuMuCCOth,     0.);
  this->GetND280Inputs()->TweakEfficiency( kNDmNuMuNCEL,      0.);
  this->GetND280Inputs()->TweakEfficiency( kNDmNuMuNCRES,     0.);
  this->GetND280Inputs()->TweakEfficiency( kNDmNuMuNCCOH,     0.);
  this->GetND280Inputs()->TweakEfficiency( kNDmNuMuNCOth,     0.);
  this->GetND280Inputs()->TweakEfficiency( kNDmNuMuBarCCQE,   0.);
  this->GetND280Inputs()->TweakEfficiency( kNDmNuMuBarCCRES,  0.);
  this->GetND280Inputs()->TweakEfficiency( kNDmNuMuBarCCCOH,  0.);
  this->GetND280Inputs()->TweakEfficiency( kNDmNuMuBarCCOth,  0.);
  this->GetND280Inputs()->TweakEfficiency( kNDmNuMuBarNCEL,   0.);
  this->GetND280Inputs()->TweakEfficiency( kNDmNuMuBarNCRES,  0.);
  this->GetND280Inputs()->TweakEfficiency( kNDmNuMuBarNCCOH,  0.);
  this->GetND280Inputs()->TweakEfficiency( kNDmNuMuBarNCOth,  0.);
  this->GetND280Inputs()->TweakEfficiency( kNDmNuECCQE,       0.);
  this->GetND280Inputs()->TweakEfficiency( kNDmNuECCRES,      0.);
  this->GetND280Inputs()->TweakEfficiency( kNDmNuECCCOH,      0.);
  this->GetND280Inputs()->TweakEfficiency( kNDmNuECCOth,      0.);
  this->GetND280Inputs()->TweakEfficiency( kNDmNuENCEL,       0.);
  this->GetND280Inputs()->TweakEfficiency( kNDmNuENCRES,      0.);
  this->GetND280Inputs()->TweakEfficiency( kNDmNuENCCOH,      0.);
  this->GetND280Inputs()->TweakEfficiency( kNDmNuENCOth,      0.);

  this->GetSuperkInputs()->TweakEnergyScale    (0.);
  this->GetSuperkInputs()->TweakEfficiency (kSKsFCFV1RingMuLike, kSKmNuMuCCQE,    0.);
  this->GetSuperkInputs()->TweakEfficiency (kSKsFCFV1RingMuLike, kSKmNuMuCC1pi,   0.);
  this->GetSuperkInputs()->TweakEfficiency (kSKsFCFV1RingMuLike, kSKmNuMuCCoth,   0.);
  this->GetSuperkInputs()->TweakEfficiency (kSKsFCFV1RingMuLike, kSKmNC,          0.);

  this->GetBeamInputs()->TweakNorm (0.);
  this->GetBeamInputs()->TweakShape(false);

  this->GetGeneratorInputs()->TweakCCQEXSec             (0.);
  this->GetGeneratorInputs()->TweakCC1piToCCQEXSecRatio (0.);
  this->GetGeneratorInputs()->TweakCCOthToCCQEXSecRatio (0.);
  this->GetGeneratorInputs()->TweakNCToCCQEXSecRatio    (0.);
  this->GetGeneratorInputs()->TweakNuclEffect           (0.);
}
//________________________________________________________________________________________
void SuperkExpectation::TweakSystematic(string systematic, double twkdial)
{
  //
  // ND280 detector systematics
  //
  if(systematic == "nd280_eff") {
     this->GetND280Inputs()->TweakEfficiency(twkdial);
  }

  //
  // SuperK detector systematics
  //
  else
  if(systematic == "sk_energy_scale") {
     this->GetSuperkInputs()->TweakEnergyScale(twkdial);
  }
  else
  if(systematic == "sk_ccqe_eff") {
     this->GetSuperkInputs()->TweakEfficiency(kSKsFCFV1RingMuLike, kSKmNuMuCCQE, twkdial);
  }
  else
  if(systematic == "sk_ccnonqe_eff") {
     this->GetSuperkInputs()->TweakEfficiency(kSKsFCFV1RingMuLike, kSKmNuMuCC1pi, twkdial);
  }
  else
  if(systematic == "sk_nc_eff") {
     this->GetSuperkInputs()->TweakEfficiency(kSKsFCFV1RingMuLike, kSKmNC, twkdial);
  }

  //
  // flux systematics
  //
  else
  if(systematic == "flux_norm") {
     this->GetBeamInputs()->TweakNorm(twkdial);
  }
  else
  if(systematic == "flux_shape") {
     this->GetBeamInputs()->TweakShape();
  }

  //
  // cross section systematics
  //
  else
  if(systematic == "xsec_ccqe") {
     this->GetGeneratorInputs()->TweakCCQEXSec(twkdial);
  }
  else
  if(systematic == "xsec_cc1pi_to_ccqe") {
     this->GetGeneratorInputs()->TweakCC1piToCCQEXSecRatio(twkdial);
  }
  else
  if(systematic == "xsec_ccoth_to_ccqe") {
     this->GetGeneratorInputs()->TweakCCOthToCCQEXSecRatio(twkdial);
  }
  else
  if(systematic == "xsec_nc_to_ccqe") {
     this->GetGeneratorInputs()->TweakNCToCCQEXSecRatio(twkdial);
  }
  else
  if(systematic == "nuclear_effect") {
     this->GetGeneratorInputs()->TweakNuclEffect(twkdial);
  }   
     
  else {
     cerr << "Exp >> Unknown systematic parameter!" << endl;
  }   
     
}
//________________________________________________________________________________________
void SuperkExpectation::TweakAllSystematics(bool do_tweak)
{
  if(do_tweak) 
  {
     TRandom3 rand(0);

     this->GetND280Inputs()->TweakEfficiency( kNDmNuMuCCQE,      rand.Gaus(0.,1.));
     this->GetND280Inputs()->TweakEfficiency( kNDmNuMuCCRES,     rand.Gaus(0.,1.));
     this->GetND280Inputs()->TweakEfficiency( kNDmNuMuCCCOH,     rand.Gaus(0.,1.));
     this->GetND280Inputs()->TweakEfficiency( kNDmNuMuCCOth,     rand.Gaus(0.,1.));
     this->GetND280Inputs()->TweakEfficiency( kNDmNuMuNCEL,      rand.Gaus(0.,1.));
     this->GetND280Inputs()->TweakEfficiency( kNDmNuMuNCRES,     rand.Gaus(0.,1.));
     this->GetND280Inputs()->TweakEfficiency( kNDmNuMuNCCOH,     rand.Gaus(0.,1.));
     this->GetND280Inputs()->TweakEfficiency( kNDmNuMuNCOth,     rand.Gaus(0.,1.));
     this->GetND280Inputs()->TweakEfficiency( kNDmNuMuBarCCQE,   rand.Gaus(0.,1.));
     this->GetND280Inputs()->TweakEfficiency( kNDmNuMuBarCCRES,  rand.Gaus(0.,1.));
     this->GetND280Inputs()->TweakEfficiency( kNDmNuMuBarCCCOH,  rand.Gaus(0.,1.));
     this->GetND280Inputs()->TweakEfficiency( kNDmNuMuBarCCOth,  rand.Gaus(0.,1.));
     this->GetND280Inputs()->TweakEfficiency( kNDmNuMuBarNCEL,   rand.Gaus(0.,1.));
     this->GetND280Inputs()->TweakEfficiency( kNDmNuMuBarNCRES,  rand.Gaus(0.,1.));
     this->GetND280Inputs()->TweakEfficiency( kNDmNuMuBarNCCOH,  rand.Gaus(0.,1.));
     this->GetND280Inputs()->TweakEfficiency( kNDmNuMuBarNCOth,  rand.Gaus(0.,1.));
     this->GetND280Inputs()->TweakEfficiency( kNDmNuECCQE,       rand.Gaus(0.,1.));
     this->GetND280Inputs()->TweakEfficiency( kNDmNuECCRES,      rand.Gaus(0.,1.));
     this->GetND280Inputs()->TweakEfficiency( kNDmNuECCCOH,      rand.Gaus(0.,1.));
     this->GetND280Inputs()->TweakEfficiency( kNDmNuECCOth,      rand.Gaus(0.,1.));
     this->GetND280Inputs()->TweakEfficiency( kNDmNuENCEL,       rand.Gaus(0.,1.));
     this->GetND280Inputs()->TweakEfficiency( kNDmNuENCRES,      rand.Gaus(0.,1.));
     this->GetND280Inputs()->TweakEfficiency( kNDmNuENCCOH,      rand.Gaus(0.,1.));
     this->GetND280Inputs()->TweakEfficiency( kNDmNuENCOth,      rand.Gaus(0.,1.));

     this->GetSuperkInputs()->TweakEnergyScale(rand.Gaus(0.,1.));
     this->GetSuperkInputs()->TweakEfficiency (kSKsFCFV1RingMuLike, kSKmNuMuCCQE,    rand.Gaus(0.,1.));
     this->GetSuperkInputs()->TweakEfficiency (kSKsFCFV1RingMuLike, kSKmNuMuCC1pi,   rand.Gaus(0.,1.));
     this->GetSuperkInputs()->TweakEfficiency (kSKsFCFV1RingMuLike, kSKmNuMuCCoth,   rand.Gaus(0.,1.));
     this->GetSuperkInputs()->TweakEfficiency (kSKsFCFV1RingMuLike, kSKmNC,          rand.Gaus(0.,1.));

     this->GetBeamInputs()->TweakNorm (rand.Gaus(0.,1.));
     this->GetBeamInputs()->TweakShape();

     this->GetGeneratorInputs()->TweakCCQEXSec             (rand.Gaus(0.,1.));
     this->GetGeneratorInputs()->TweakCC1piToCCQEXSecRatio (rand.Gaus(0.,1.));
     this->GetGeneratorInputs()->TweakCCOthToCCQEXSecRatio (rand.Gaus(0.,1.));
     this->GetGeneratorInputs()->TweakNCToCCQEXSecRatio    (rand.Gaus(0.,1.));
     this->GetGeneratorInputs()->TweakNuclEffect           (rand.Gaus(0.,1.));

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

  fND280Inputs     = new ND280Inputs;
  fSuperkInputs    = new SuperkInputs;
  fBeamInputs      = new BeamInputs; 
  fGeneratorInputs = new GeneratorInputs; 

  fND280Inputs  -> TellMeAboutCrossSections (fGeneratorInputs);
  fSuperkInputs -> TellMeAboutCrossSections (fGeneratorInputs);
  fSuperkInputs -> TellMeAboutTheFlux       (fBeamInputs);

  const unsigned int ns = fSamples.size();
  const unsigned int nm = 6;  
  SuperkMode_t sk_modes [nm] = { 
     kSKmNuMuCCQE, kSKmNuMuCC1pi, kSKmNuMuCCoth, kSKmNC, kSKmNuMuBar, kSKmNuE 
  };
  fPredictErecoSpectrumGrandTot = new TH1D("","", kNvErecoBins, kNvErecoBinEdges);
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

  if ( ! this->GetND280Inputs    ()->Read(xml, root_node) ) return false;
  if ( ! this->GetSuperkInputs   ()->Read(xml, root_node) ) return false;
  if ( ! this->GetBeamInputs     ()->Read(xml, root_node) ) return false;
  if ( ! this->GetGeneratorInputs()->Read(xml, root_node) ) return false;

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

  fND280Inputs     = 0;
  fSuperkInputs    = 0;
  fBeamInputs      = 0;
  fGeneratorInputs = 0;

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
  if (fND280Inputs    ) delete fND280Inputs;
  if (fSuperkInputs   ) delete fSuperkInputs;
  if (fBeamInputs     ) delete fBeamInputs; 
  if (fGeneratorInputs) delete fGeneratorInputs;

  fSamples.clear();

  if (fPredictErecoSpectrumGrandTot) delete fPredictErecoSpectrumGrandTot;
  for(int is=0; is<kns; is++) {
    if (fPredictErecoSpectrumTot[is]) delete fPredictErecoSpectrumTot[is];
    for(int im=0; im<knm; im++) {
      if (fPredictErecoSpectrum[is][im]) delete fPredictErecoSpectrum[is][im];
    }//im
  }//is
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
  fPOTData       = 0;
  fNumObsCCLike  = 0;
  fPOTMC         = 0;
  for(int i = 0; i < knmnd; i++) {
    fNumSelEt[i] = 0;
    fSelEffEt[i] = 0;
    fSelEffTwkDial[i] = 0;
  }
  fSelEffRelErr  = 0.;
  fGenInputs     = 0; 

  fNuMuEt    = 0;
  fNuMuBarEt = 0;
  fNuEEt     = 0;

  fNomCorr   = 0;
}
//________________________________________________________________________________________
void SuperkExpectation::ND280Inputs::CleanUp(void)
{
  for(int i = 0; i < knmnd; i++) {
    if (fNumSelEt[i]) delete fNumSelEt[i];
    if (fSelEffEt[i]) delete fSelEffEt[i];
  }

  if ( fNuMuEt    ) delete fNuMuEt;
  if ( fNuMuBarEt ) delete fNuMuBarEt;
  if ( fNuEEt     ) delete fNuEEt;
}
//________________________________________________________________________________________
bool SuperkExpectation::ND280Inputs::Read(TXMLEngine* xml, XMLNodePointer_t root_node)
{
  cout << "Exp/nd280 >> Reading nd280 inputs from input XML file." << endl;

  // Get Data POT
  fPOTData = util::XML2Dbl(xml,root_node,"nd280/pot_data");
  cout << "Exp/nd280 >> " << "N_{POT;data} : " << fPOTData << endl;

  // Get observed number of CC-like events
  fNumObsCCLike = util::XML2Dbl(xml,root_node,"nd280/saclay_selection/nobs_cclike");
  cout << "Exp/nd280 >> " << "N^{obs}_{CClike} : " << fNumObsCCLike << endl;

  // Get MC POT
  fPOTMC = util::XML2Dbl(xml,root_node,"nd280/saclay_selection/pot_mc");
  cout << "Exp/nd280 >> " << "N_{POT;mc} : " << fPOTMC << endl;

  // Get NMC_{mode}(Etrue) after cuts
  cout << "Exp/nd280 >> " << "Getting NMC_{mode}(Etrue) histograms" << endl;
  fNumSelEt [kNDmNuMuCCQE]      =
    util::XML2TH1F(xml,root_node,"nd280/saclay_selection/true_energy_spectrum/numu_ccqe");
  fNumSelEt [kNDmNuMuCCRES]     =
    util::XML2TH1F(xml,root_node,"nd280/saclay_selection/true_energy_spectrum/numu_ccres");
  fNumSelEt [kNDmNuMuCCCOH]     =
    util::XML2TH1F(xml,root_node,"nd280/saclay_selection/true_energy_spectrum/numu_cccoh");
  fNumSelEt [kNDmNuMuCCOth]     =
    util::XML2TH1F(xml,root_node,"nd280/saclay_selection/true_energy_spectrum/numu_ccoth");
  fNumSelEt [kNDmNuMuNCEL]      =
    util::XML2TH1F(xml,root_node,"nd280/saclay_selection/true_energy_spectrum/numu_ncel");
  fNumSelEt [kNDmNuMuNCRES]     = 
    util::XML2TH1F(xml,root_node,"nd280/saclay_selection/true_energy_spectrum/numu_ncres");
  fNumSelEt [kNDmNuMuNCCOH]     =
    util::XML2TH1F(xml,root_node,"nd280/saclay_selection/true_energy_spectrum/numu_nccoh");
  fNumSelEt [kNDmNuMuNCOth]     =
    util::XML2TH1F(xml,root_node,"nd280/saclay_selection/true_energy_spectrum/numu_ncoth");
  fNumSelEt [kNDmNuMuBarCCQE]   = 
    util::XML2TH1F(xml,root_node,"nd280/saclay_selection/true_energy_spectrum/numubar_ccqe");
  fNumSelEt [kNDmNuMuBarCCRES]  = 
    util::XML2TH1F(xml,root_node,"nd280/saclay_selection/true_energy_spectrum/numubar_ccres");
  fNumSelEt [kNDmNuMuBarCCCOH]  =
    util::XML2TH1F(xml,root_node,"nd280/saclay_selection/true_energy_spectrum/numubar_cccoh");
  fNumSelEt [kNDmNuMuBarCCOth]  =
    util::XML2TH1F(xml,root_node,"nd280/saclay_selection/true_energy_spectrum/numubar_ccoth");
  fNumSelEt [kNDmNuMuBarNCEL]   =
    util::XML2TH1F(xml,root_node,"nd280/saclay_selection/true_energy_spectrum/numubar_ncel");
  fNumSelEt [kNDmNuMuBarNCRES]  = 
    util::XML2TH1F(xml,root_node,"nd280/saclay_selection/true_energy_spectrum/numubar_ncres");
  fNumSelEt [kNDmNuMuBarNCCOH]  =
    util::XML2TH1F(xml,root_node,"nd280/saclay_selection/true_energy_spectrum/numubar_nccoh");
  fNumSelEt [kNDmNuMuBarNCOth]  =
    util::XML2TH1F(xml,root_node,"nd280/saclay_selection/true_energy_spectrum/numubar_ncoth");
  fNumSelEt [kNDmNuECCQE]       =
    util::XML2TH1F(xml,root_node,"nd280/saclay_selection/true_energy_spectrum/nue_ccqe");
  fNumSelEt [kNDmNuECCRES]      =
    util::XML2TH1F(xml,root_node,"nd280/saclay_selection/true_energy_spectrum/nue_ccres");
  fNumSelEt [kNDmNuECCCOH]      =
    util::XML2TH1F(xml,root_node,"nd280/saclay_selection/true_energy_spectrum/nue_cccoh");
  fNumSelEt [kNDmNuECCOth]      =
    util::XML2TH1F(xml,root_node,"nd280/saclay_selection/true_energy_spectrum/nue_ccoth");
  fNumSelEt [kNDmNuENCEL]       =
    util::XML2TH1F(xml,root_node,"nd280/saclay_selection/true_energy_spectrum/nue_ncel");
  fNumSelEt [kNDmNuENCRES]      =
    util::XML2TH1F(xml,root_node,"nd280/saclay_selection/true_energy_spectrum/nue_ncres");
  fNumSelEt [kNDmNuENCCOH]      =
    util::XML2TH1F(xml,root_node,"nd280/saclay_selection/true_energy_spectrum/nue_nccoh");
  fNumSelEt [kNDmNuENCOth]      =
    util::XML2TH1F(xml,root_node,"nd280/saclay_selection/true_energy_spectrum/nue_ncoth");

  cout << "Exp/nd280 >> Done loading NMC_{mode}(Etrue) histograms" << endl;

  double data_to_mc_pot = fPOTData/fPOTMC;
  cout << "Exp/nd280 >> " << "MC POT scaling factor : " << data_to_mc_pot << endl;

  fNumSelEt [kNDmNuMuCCQE]     -> Scale (data_to_mc_pot);
  fNumSelEt [kNDmNuMuCCRES]    -> Scale (data_to_mc_pot);
  fNumSelEt [kNDmNuMuCCCOH]    -> Scale (data_to_mc_pot);
  fNumSelEt [kNDmNuMuCCOth]    -> Scale (data_to_mc_pot);
  fNumSelEt [kNDmNuMuNCEL]     -> Scale (data_to_mc_pot);
  fNumSelEt [kNDmNuMuNCRES]    -> Scale (data_to_mc_pot);
  fNumSelEt [kNDmNuMuNCCOH]    -> Scale (data_to_mc_pot);
  fNumSelEt [kNDmNuMuNCOth]    -> Scale (data_to_mc_pot);
  fNumSelEt [kNDmNuMuBarCCQE]  -> Scale (data_to_mc_pot);
  fNumSelEt [kNDmNuMuBarCCRES] -> Scale (data_to_mc_pot);
  fNumSelEt [kNDmNuMuBarCCCOH] -> Scale (data_to_mc_pot);
  fNumSelEt [kNDmNuMuBarCCOth] -> Scale (data_to_mc_pot);
  fNumSelEt [kNDmNuMuBarNCEL]  -> Scale (data_to_mc_pot);
  fNumSelEt [kNDmNuMuBarNCRES] -> Scale (data_to_mc_pot);
  fNumSelEt [kNDmNuMuBarNCCOH] -> Scale (data_to_mc_pot);
  fNumSelEt [kNDmNuMuBarNCOth] -> Scale (data_to_mc_pot);
  fNumSelEt [kNDmNuECCQE]      -> Scale (data_to_mc_pot);
  fNumSelEt [kNDmNuECCRES]     -> Scale (data_to_mc_pot);
  fNumSelEt [kNDmNuECCCOH]     -> Scale (data_to_mc_pot);
  fNumSelEt [kNDmNuECCOth]     -> Scale (data_to_mc_pot);
  fNumSelEt [kNDmNuENCEL]      -> Scale (data_to_mc_pot);
  fNumSelEt [kNDmNuENCRES]     -> Scale (data_to_mc_pot);
  fNumSelEt [kNDmNuENCCOH]     -> Scale (data_to_mc_pot);
  fNumSelEt [kNDmNuENCOth]     -> Scale (data_to_mc_pot);

  cout << "Exp/nd280 >> Done rescaling NMC_{mode} histograms to the data POT" << endl;

  // Get efficiency_{mode}(Etrue)
  cout << "Exp/nd280 >> " << "Getting efficiency_{mode}(Etrue) histograms" << endl; 
  fSelEffEt [kNDmNuMuCCQE]      =
    util::XML2TH1F(xml,root_node,"nd280/saclay_selection/efficiency_vs_true_energy/numu_ccqe");
  fSelEffEt [kNDmNuMuCCRES]     =
    util::XML2TH1F(xml,root_node,"nd280/saclay_selection/efficiency_vs_true_energy/numu_ccres");
  fSelEffEt [kNDmNuMuCCCOH]     =
    util::XML2TH1F(xml,root_node,"nd280/saclay_selection/efficiency_vs_true_energy/numu_cccoh");
  fSelEffEt [kNDmNuMuCCOth]     =
    util::XML2TH1F(xml,root_node,"nd280/saclay_selection/efficiency_vs_true_energy/numu_ccoth");
  fSelEffEt [kNDmNuMuNCEL]      =
    util::XML2TH1F(xml,root_node,"nd280/saclay_selection/efficiency_vs_true_energy/numu_ncel");
  fSelEffEt [kNDmNuMuNCRES]     = 
    util::XML2TH1F(xml,root_node,"nd280/saclay_selection/efficiency_vs_true_energy/numu_ncres");
  fSelEffEt [kNDmNuMuNCCOH]     =
    util::XML2TH1F(xml,root_node,"nd280/saclay_selection/efficiency_vs_true_energy/numu_nccoh");
  fSelEffEt [kNDmNuMuNCOth]     =
    util::XML2TH1F(xml,root_node,"nd280/saclay_selection/efficiency_vs_true_energy/numu_ncoth");
  fSelEffEt [kNDmNuMuBarCCQE]   = 
    util::XML2TH1F(xml,root_node,"nd280/saclay_selection/efficiency_vs_true_energy/numubar_ccqe");
  fSelEffEt [kNDmNuMuBarCCRES]  = 
    util::XML2TH1F(xml,root_node,"nd280/saclay_selection/efficiency_vs_true_energy/numubar_ccres");
  fSelEffEt [kNDmNuMuBarCCCOH]  =
    util::XML2TH1F(xml,root_node,"nd280/saclay_selection/efficiency_vs_true_energy/numubar_cccoh");
  fSelEffEt [kNDmNuMuBarCCOth]  =
    util::XML2TH1F(xml,root_node,"nd280/saclay_selection/efficiency_vs_true_energy/numubar_ccoth");
  fSelEffEt [kNDmNuMuBarNCEL]   =
    util::XML2TH1F(xml,root_node,"nd280/saclay_selection/efficiency_vs_true_energy/numubar_ncel");
  fSelEffEt [kNDmNuMuBarNCRES]  = 
    util::XML2TH1F(xml,root_node,"nd280/saclay_selection/efficiency_vs_true_energy/numubar_ncres");
  fSelEffEt [kNDmNuMuBarNCCOH]  =
    util::XML2TH1F(xml,root_node,"nd280/saclay_selection/efficiency_vs_true_energy/numubar_nccoh");
  fSelEffEt [kNDmNuMuBarNCOth]  =
    util::XML2TH1F(xml,root_node,"nd280/saclay_selection/efficiency_vs_true_energy/numubar_ncoth");
  fSelEffEt [kNDmNuECCQE]       =
    util::XML2TH1F(xml,root_node,"nd280/saclay_selection/efficiency_vs_true_energy/nue_ccqe");
  fSelEffEt [kNDmNuECCRES]      =
    util::XML2TH1F(xml,root_node,"nd280/saclay_selection/efficiency_vs_true_energy/nue_ccres");
  fSelEffEt [kNDmNuECCCOH]      =
    util::XML2TH1F(xml,root_node,"nd280/saclay_selection/efficiency_vs_true_energy/nue_cccoh");
  fSelEffEt [kNDmNuECCOth]      =
    util::XML2TH1F(xml,root_node,"nd280/saclay_selection/efficiency_vs_true_energy/nue_ccoth");
  fSelEffEt [kNDmNuENCEL]       =
    util::XML2TH1F(xml,root_node,"nd280/saclay_selection/efficiency_vs_true_energy/nue_ncel");
  fSelEffEt [kNDmNuENCRES]      =
    util::XML2TH1F(xml,root_node,"nd280/saclay_selection/efficiency_vs_true_energy/nue_ncres");
  fSelEffEt [kNDmNuENCCOH]      =
    util::XML2TH1F(xml,root_node,"nd280/saclay_selection/efficiency_vs_true_energy/nue_nccoh");
  fSelEffEt [kNDmNuENCOth]      =
    util::XML2TH1F(xml,root_node,"nd280/saclay_selection/efficiency_vs_true_energy/nue_ncoth");

  // Get error on efficiency
  fSelEffRelErr = util::XML2Dbl(xml,root_node,"nd280/saclay_selection/systematics/efficiency_relative_error");
  cout << "Exp/nd280 >> delta{efficiency} : " << fSelEffRelErr << endl;

  // Set colors and line styles for drawing

  fNumSelEt [kNDmNuMuCCQE]     -> SetLineColor (kBlack);
  fNumSelEt [kNDmNuMuCCRES]    -> SetLineColor (kBlack); 
  fNumSelEt [kNDmNuMuCCCOH]    -> SetLineColor (kBlack); 
  fNumSelEt [kNDmNuMuCCOth]    -> SetLineColor (kBlack); 
  fNumSelEt [kNDmNuMuNCEL]     -> SetLineColor (kRed); 
  fNumSelEt [kNDmNuMuNCRES]    -> SetLineColor (kRed); 
  fNumSelEt [kNDmNuMuNCCOH]    -> SetLineColor (kRed); 
  fNumSelEt [kNDmNuMuNCOth]    -> SetLineColor (kRed); 
  fNumSelEt [kNDmNuMuBarCCQE]  -> SetLineColor (kBlue); 
  fNumSelEt [kNDmNuMuBarCCRES] -> SetLineColor (kBlue); 
  fNumSelEt [kNDmNuMuBarCCCOH] -> SetLineColor (kBlue); 
  fNumSelEt [kNDmNuMuBarCCOth] -> SetLineColor (kBlue); 
  fNumSelEt [kNDmNuMuBarNCEL]  -> SetLineColor (kMagenta); 
  fNumSelEt [kNDmNuMuBarNCRES] -> SetLineColor (kMagenta); 
  fNumSelEt [kNDmNuMuBarNCCOH] -> SetLineColor (kMagenta); 
  fNumSelEt [kNDmNuMuBarNCOth] -> SetLineColor (kMagenta); 
  fNumSelEt [kNDmNuECCQE]      -> SetLineColor (kGreen); 
  fNumSelEt [kNDmNuECCRES]     -> SetLineColor (kGreen); 
  fNumSelEt [kNDmNuECCCOH]     -> SetLineColor (kGreen); 
  fNumSelEt [kNDmNuECCOth]     -> SetLineColor (kGreen); 
  fNumSelEt [kNDmNuENCEL]      -> SetLineColor (kPink); 
  fNumSelEt [kNDmNuENCRES]     -> SetLineColor (kPink); 
  fNumSelEt [kNDmNuENCCOH]     -> SetLineColor (kPink); 
  fNumSelEt [kNDmNuENCOth]     -> SetLineColor (kPink); 

  fNumSelEt [kNDmNuMuCCQE]     -> SetLineStyle (kSolid);
  fNumSelEt [kNDmNuMuCCRES]    -> SetLineStyle (kDashed); 
  fNumSelEt [kNDmNuMuCCCOH]    -> SetLineStyle (kDotted); 
  fNumSelEt [kNDmNuMuCCOth]    -> SetLineStyle (kDashDotted); 
  fNumSelEt [kNDmNuMuNCEL]     -> SetLineStyle (kSolid); 
  fNumSelEt [kNDmNuMuNCRES]    -> SetLineStyle (kDashed); 
  fNumSelEt [kNDmNuMuNCCOH]    -> SetLineStyle (kDotted); 
  fNumSelEt [kNDmNuMuNCOth]    -> SetLineStyle (kDashDotted); 
  fNumSelEt [kNDmNuMuBarCCQE]  -> SetLineStyle (kSolid); 
  fNumSelEt [kNDmNuMuBarCCRES] -> SetLineStyle (kDashed); 
  fNumSelEt [kNDmNuMuBarCCCOH] -> SetLineStyle (kDotted); 
  fNumSelEt [kNDmNuMuBarCCOth] -> SetLineStyle (kDashDotted); 
  fNumSelEt [kNDmNuMuBarNCEL]  -> SetLineStyle (kSolid); 
  fNumSelEt [kNDmNuMuBarNCRES] -> SetLineStyle (kDashed); 
  fNumSelEt [kNDmNuMuBarNCCOH] -> SetLineStyle (kDotted); 
  fNumSelEt [kNDmNuMuBarNCOth] -> SetLineStyle (kDashDotted); 
  fNumSelEt [kNDmNuECCQE]      -> SetLineStyle (kSolid); 
  fNumSelEt [kNDmNuECCRES]     -> SetLineStyle (kDashed); 
  fNumSelEt [kNDmNuECCCOH]     -> SetLineStyle (kDotted); 
  fNumSelEt [kNDmNuECCOth]     -> SetLineStyle (kDashDotted); 
  fNumSelEt [kNDmNuENCEL]      -> SetLineStyle (kSolid); 
  fNumSelEt [kNDmNuENCRES]     -> SetLineStyle (kDashed); 
  fNumSelEt [kNDmNuENCCOH]     -> SetLineStyle (kDotted); 
  fNumSelEt [kNDmNuENCOth]     -> SetLineStyle (kDashDotted); 

  fSelEffEt [kNDmNuMuCCQE]     -> SetLineColor (kBlack);
  fSelEffEt [kNDmNuMuCCRES]    -> SetLineColor (kBlack); 
  fSelEffEt [kNDmNuMuCCCOH]    -> SetLineColor (kBlack); 
  fSelEffEt [kNDmNuMuCCOth]    -> SetLineColor (kBlack); 
  fSelEffEt [kNDmNuMuNCEL]     -> SetLineColor (kRed); 
  fSelEffEt [kNDmNuMuNCRES]    -> SetLineColor (kRed); 
  fSelEffEt [kNDmNuMuNCCOH]    -> SetLineColor (kRed); 
  fSelEffEt [kNDmNuMuNCOth]    -> SetLineColor (kRed); 
  fSelEffEt [kNDmNuMuBarCCQE]  -> SetLineColor (kBlue); 
  fSelEffEt [kNDmNuMuBarCCRES] -> SetLineColor (kBlue); 
  fSelEffEt [kNDmNuMuBarCCCOH] -> SetLineColor (kBlue); 
  fSelEffEt [kNDmNuMuBarCCOth] -> SetLineColor (kBlue); 
  fSelEffEt [kNDmNuMuBarNCEL]  -> SetLineColor (kMagenta); 
  fSelEffEt [kNDmNuMuBarNCRES] -> SetLineColor (kMagenta); 
  fSelEffEt [kNDmNuMuBarNCCOH] -> SetLineColor (kMagenta); 
  fSelEffEt [kNDmNuMuBarNCOth] -> SetLineColor (kMagenta); 
  fSelEffEt [kNDmNuECCQE]      -> SetLineColor (kGreen); 
  fSelEffEt [kNDmNuECCRES]     -> SetLineColor (kGreen); 
  fSelEffEt [kNDmNuECCCOH]     -> SetLineColor (kGreen); 
  fSelEffEt [kNDmNuECCOth]     -> SetLineColor (kGreen); 
  fSelEffEt [kNDmNuENCEL]      -> SetLineColor (kPink); 
  fSelEffEt [kNDmNuENCRES]     -> SetLineColor (kPink); 
  fSelEffEt [kNDmNuENCCOH]     -> SetLineColor (kPink); 
  fSelEffEt [kNDmNuENCOth]     -> SetLineColor (kPink); 

  fSelEffEt [kNDmNuMuCCQE]     -> SetLineStyle (kSolid);
  fSelEffEt [kNDmNuMuCCRES]    -> SetLineStyle (kDashed); 
  fSelEffEt [kNDmNuMuCCCOH]    -> SetLineStyle (kDotted); 
  fSelEffEt [kNDmNuMuCCOth]    -> SetLineStyle (kDashDotted); 
  fSelEffEt [kNDmNuMuNCEL]     -> SetLineStyle (kSolid); 
  fSelEffEt [kNDmNuMuNCRES]    -> SetLineStyle (kDashed); 
  fSelEffEt [kNDmNuMuNCCOH]    -> SetLineStyle (kDotted); 
  fSelEffEt [kNDmNuMuNCOth]    -> SetLineStyle (kDashDotted); 
  fSelEffEt [kNDmNuMuBarCCQE]  -> SetLineStyle (kSolid); 
  fSelEffEt [kNDmNuMuBarCCRES] -> SetLineStyle (kDashed); 
  fSelEffEt [kNDmNuMuBarCCCOH] -> SetLineStyle (kDotted); 
  fSelEffEt [kNDmNuMuBarCCOth] -> SetLineStyle (kDashDotted); 
  fSelEffEt [kNDmNuMuBarNCEL]  -> SetLineStyle (kSolid); 
  fSelEffEt [kNDmNuMuBarNCRES] -> SetLineStyle (kDashed); 
  fSelEffEt [kNDmNuMuBarNCCOH] -> SetLineStyle (kDotted); 
  fSelEffEt [kNDmNuMuBarNCOth] -> SetLineStyle (kDashDotted); 
  fSelEffEt [kNDmNuECCQE]      -> SetLineStyle (kSolid); 
  fSelEffEt [kNDmNuECCRES]     -> SetLineStyle (kDashed); 
  fSelEffEt [kNDmNuECCCOH]     -> SetLineStyle (kDotted); 
  fSelEffEt [kNDmNuECCOth]     -> SetLineStyle (kDashDotted); 
  fSelEffEt [kNDmNuENCEL]      -> SetLineStyle (kSolid); 
  fSelEffEt [kNDmNuENCRES]     -> SetLineStyle (kDashed); 
  fSelEffEt [kNDmNuENCCOH]     -> SetLineStyle (kDotted); 
  fSelEffEt [kNDmNuENCOth]     -> SetLineStyle (kDashDotted); 

  //
  // histograms needed for beam norm err calculation
  //

  fNuMuEt    = new TH1D("fNuMuEt",    "Etrue for numu",    200, 0, 10); /// binning expected by Kubo's code
  fNuMuBarEt = new TH1D("fNuMuBarEt", "Etrue for numubar", 200, 0, 10);
  fNuEEt     = new TH1D("fNuMuEEt",   "Etrue for nue",     200, 0, 10);

  cout << "Exp/nd280 >> Done reading nd280 inputs!" << endl;

  //
  // store nominal correction factor 
  //

  bool nominal = true;
  double nobscclike_nd280  = this -> NObsCCLike();
  double nexpncbkg_nd280   = this -> NExpNCBkg (nominal);
  double nexpcc_nd280      = this -> NExpCC    (nominal);
  fNomCorr = (nobscclike_nd280 - nexpncbkg_nd280) / nexpcc_nd280;
  cout << "Exp/nd280 >> Nominal correction factor from nd280: " << fNomCorr << endl;

  return true;
}
//________________________________________________________________________________________
void SuperkExpectation::ND280Inputs::Draw(TVirtualPad * c, Option_t * option) const
{
  if(!c) return;

  c->cd();

  util::SetDefaultStyle();

  TString opt = option;

  double Emax = 10000;
  double Nmax = 100;
  double Nmin = 0.001;

  if(opt.Contains("Saclay/EtrueSpectrum/NuMu/CC")) { 
    TH1F * hframe = (TH1F*)c->DrawFrame(0,Nmin,Emax,Nmax);
    hframe -> GetXaxis()->SetTitle("E_{#nu} (GeV)");
    hframe -> GetYaxis()->SetTitle(Form("events / %E POT", fPOTData));
    hframe->Draw();
    fNumSelEt [kNDmNuMuCCQE]     -> Draw("LPERRSAME");
    fNumSelEt [kNDmNuMuCCRES]    -> Draw("LPERRSAME");
    fNumSelEt [kNDmNuMuCCCOH]    -> Draw("LPERRSAME");
    fNumSelEt [kNDmNuMuCCOth]    -> Draw("LPERRSAME");
    TLegend * l = new TLegend(0.7,0.6,0.9,0.9);
    l->SetFillColor(0);
    l->AddEntry(fNumSelEt [kNDmNuMuCCQE] , "#nu_{#mu} CCQE",    "LP");
    l->AddEntry(fNumSelEt [kNDmNuMuCCRES], "#nu_{#mu} CCRES",   "LP");
    l->AddEntry(fNumSelEt [kNDmNuMuCCCOH], "#nu_{#mu} CCCOH",   "LP");
    l->AddEntry(fNumSelEt [kNDmNuMuCCOth], "#nu_{#mu} CCOther", "LP");
    l->Draw();
    c->SetLogy();
  }
  else
  if(opt.Contains("Saclay/EtrueSpectrum/NuMu/NC")) { 
    TH1F * hframe = (TH1F*)c->DrawFrame(0,Nmin,Emax,Nmax);
    hframe->Draw();
    fNumSelEt [kNDmNuMuNCEL]     -> Draw("LPERRSAME");
    fNumSelEt [kNDmNuMuNCRES]    -> Draw("LPERRSAME");
    fNumSelEt [kNDmNuMuNCCOH]    -> Draw("LPERRSAME");
    fNumSelEt [kNDmNuMuNCOth]    -> Draw("LPERRSAME");
    fNumSelEt [kNDmNuMuNCEL]     -> GetXaxis()->SetTitle("E_{#nu} (GeV)");
    fNumSelEt [kNDmNuMuNCEL]     -> GetYaxis()->SetTitle(Form("events / %E POT", fPOTData));
    TLegend * l = new TLegend(0.7,0.6,0.9,0.9);
    l->SetFillColor(0);
    l->AddEntry(fNumSelEt [kNDmNuMuNCEL] , "#nu_{#mu} NCEL",    "LP");
    l->AddEntry(fNumSelEt [kNDmNuMuNCRES], "#nu_{#mu} NCRES",   "LP");
    l->AddEntry(fNumSelEt [kNDmNuMuNCCOH], "#nu_{#mu} NCCOH",   "LP");
    l->AddEntry(fNumSelEt [kNDmNuMuNCOth], "#nu_{#mu} NCOther", "LP");
    l->Draw();
    c->SetLogy();
  }
  else
  if(opt.Contains("Saclay/EtrueSpectrum/NuMuBar/CC")) { 
    TH1F * hframe = (TH1F*)c->DrawFrame(0,Nmin,Emax,Nmax);
    hframe->Draw();
    fNumSelEt [kNDmNuMuBarCCQE]     -> Draw("LPERRSAME");
    fNumSelEt [kNDmNuMuBarCCRES]    -> Draw("LPERRSAME");
    fNumSelEt [kNDmNuMuBarCCCOH]    -> Draw("LPERRSAME");
    fNumSelEt [kNDmNuMuBarCCOth]    -> Draw("LPERRSAME");
    fNumSelEt [kNDmNuMuBarCCQE]     -> GetXaxis()->SetTitle("E_{#nu} (GeV)");
    fNumSelEt [kNDmNuMuBarCCQE]     -> GetYaxis()->SetTitle(Form("events / %E POT", fPOTData));
    TLegend * l = new TLegend(0.7,0.6,0.9,0.9);
    l->SetFillColor(0);
    l->AddEntry(fNumSelEt [kNDmNuMuBarCCQE] , "#bar{#nu_{#mu}} CCQE",    "LP");
    l->AddEntry(fNumSelEt [kNDmNuMuBarCCRES], "#bar{#nu_{#mu}} CCRES",   "LP");
    l->AddEntry(fNumSelEt [kNDmNuMuBarCCCOH], "#bar{#nu_{#mu}} CCCOH",   "LP");
    l->AddEntry(fNumSelEt [kNDmNuMuBarCCOth], "#bar{#nu_{#mu}} CCOther", "LP");
    l->Draw();
    c->SetLogy();
  }
  else
  if(opt.Contains("Saclay/EtrueSpectrum/NuMuBar/NC")) { 
    TH1F * hframe = (TH1F*)c->DrawFrame(0,Nmin,Emax,Nmax);
    hframe->Draw();
    fNumSelEt [kNDmNuMuBarNCEL]     -> Draw("LPERRSAME");
    fNumSelEt [kNDmNuMuBarNCRES]    -> Draw("LPERRSAME");
    fNumSelEt [kNDmNuMuBarNCCOH]    -> Draw("LPERRSAME");
    fNumSelEt [kNDmNuMuBarNCOth]    -> Draw("LPERRSAME");
    fNumSelEt [kNDmNuMuBarNCEL]     -> GetXaxis()->SetTitle("E_{#nu} (GeV)");
    fNumSelEt [kNDmNuMuBarNCEL]     -> GetYaxis()->SetTitle(Form("events / %E POT", fPOTData));
    TLegend * l = new TLegend(0.7,0.6,0.9,0.9);
    l->SetFillColor(0);
    l->AddEntry(fNumSelEt [kNDmNuMuBarNCEL] , "#bar{#nu_{#mu}} NCEL",    "LP");
    l->AddEntry(fNumSelEt [kNDmNuMuBarNCRES], "#bar{#nu_{#mu}} NCRES",   "LP");
    l->AddEntry(fNumSelEt [kNDmNuMuBarNCCOH], "#bar{#nu_{#mu}} NCCOH",   "LP");
    l->AddEntry(fNumSelEt [kNDmNuMuBarNCOth], "#bar{#nu_{#mu}} NCOther", "LP");
    l->Draw();
    c->SetLogy();
  }
  else
  if(opt.Contains("Saclay/EtrueSpectrum/NuE/CC")) { 
    TH1F * hframe = (TH1F*)c->DrawFrame(0,Nmin,Emax,Nmax);
    hframe->Draw();
    fNumSelEt [kNDmNuECCQE]     -> Draw("LPERRSAME");
    fNumSelEt [kNDmNuECCRES]    -> Draw("LPERRSAME");
    fNumSelEt [kNDmNuECCCOH]    -> Draw("LPERRSAME");
    fNumSelEt [kNDmNuECCOth]    -> Draw("LPERRSAME");
    fNumSelEt [kNDmNuECCQE]     -> GetXaxis()->SetTitle("E_{#nu} (GeV)");
    fNumSelEt [kNDmNuECCQE]     -> GetYaxis()->SetTitle(Form("events / %E POT", fPOTData));
    TLegend * l = new TLegend(0.7,0.6,0.9,0.9);
    l->SetFillColor(0);
    l->AddEntry(fNumSelEt [kNDmNuECCQE] , "#nu_{e} CCQE",    "LP");
    l->AddEntry(fNumSelEt [kNDmNuECCRES], "#nu_{e} CCRES",   "LP");
    l->AddEntry(fNumSelEt [kNDmNuECCCOH], "#nu_{e} CCCOH",   "LP");
    l->AddEntry(fNumSelEt [kNDmNuECCOth], "#nu_{e} CCOther", "LP");
    l->Draw();
    c->SetLogy();
  }
  else
  if(opt.Contains("Saclay/EtrueSpectrum/NuE/NC")) { 
    TH1F * hframe = (TH1F*)c->DrawFrame(0,Nmin,Emax,Nmax);
    hframe->Draw();
    fNumSelEt [kNDmNuENCEL]     -> Draw("LPERRSAME");
    fNumSelEt [kNDmNuENCRES]    -> Draw("LPERRSAME");
    fNumSelEt [kNDmNuENCCOH]    -> Draw("LPERRSAME");
    fNumSelEt [kNDmNuENCOth]    -> Draw("LPERRSAME");
    fNumSelEt [kNDmNuENCEL]     -> GetXaxis()->SetTitle("E_{#nu} (GeV)");
    fNumSelEt [kNDmNuENCEL]     -> GetYaxis()->SetTitle(Form("events / %E POT", fPOTData));
    TLegend * l = new TLegend(0.7,0.6,0.9,0.9);
    l->SetFillColor(0);
    l->AddEntry(fNumSelEt [kNDmNuENCEL] , "#nu_{e} NCEL",    "LP");
    l->AddEntry(fNumSelEt [kNDmNuENCRES], "#nu_{e} NCRES",   "LP");
    l->AddEntry(fNumSelEt [kNDmNuENCCOH], "#nu_{e} NCCOH",   "LP");
    l->AddEntry(fNumSelEt [kNDmNuENCOth], "#nu_{e} NCOther", "LP");
    l->Draw();
    c->SetLogy();
  }
  else
  if(opt.Contains("Saclay/Efficiency/NuMu/CC")) { 
    TH1F * hframe = (TH1F*)c->DrawFrame(0,0.01,Emax,1);
    hframe->Draw();
    fSelEffEt [kNDmNuMuCCQE]   -> Draw("LPERRSAME");
    fSelEffEt [kNDmNuMuCCRES]  -> Draw("LPERRSAME");
    fSelEffEt [kNDmNuMuCCCOH]  -> Draw("LPERRSAME");
    fSelEffEt [kNDmNuMuCCOth]  -> Draw("LPERRSAME");
    fSelEffEt [kNDmNuMuCCQE]   -> GetXaxis()->SetTitle("E_{#nu} (GeV)");
    fSelEffEt [kNDmNuMuCCQE]   -> GetYaxis()->SetTitle("efficiency");
    TLegend * l = new TLegend(0.7,0.6,0.9,0.9);
    l->SetFillColor(0);
    l->AddEntry(fSelEffEt [kNDmNuMuCCQE] , "#nu_{#mu} CCQE",    "LP");
    l->AddEntry(fSelEffEt [kNDmNuMuCCRES], "#nu_{#mu} CCRES",   "LP");
    l->AddEntry(fSelEffEt [kNDmNuMuCCCOH], "#nu_{#mu} CCCOH",   "LP");
    l->AddEntry(fSelEffEt [kNDmNuMuCCOth], "#nu_{#mu} CCOther", "LP");
    l->Draw();
    c->SetLogy();
  }
  else
  if(opt.Contains("Saclay/Efficiency/NuMu/NC")) { 
    TH1F * hframe = (TH1F*)c->DrawFrame(0,.01,Emax,1);
    hframe->Draw();
    fSelEffEt [kNDmNuMuNCEL]     -> Draw("LPERRSAME");
    fSelEffEt [kNDmNuMuNCRES]    -> Draw("LPERRSAME");
    fSelEffEt [kNDmNuMuNCCOH]    -> Draw("LPERRSAME");
    fSelEffEt [kNDmNuMuNCOth]    -> Draw("LPERRSAME");
    fSelEffEt [kNDmNuMuNCEL]     -> GetXaxis()->SetTitle("E_{#nu} (GeV)");
    fSelEffEt [kNDmNuMuNCEL]     -> GetYaxis()->SetTitle("efficiency");
    TLegend * l = new TLegend(0.7,0.6,0.9,0.9);
    l->SetFillColor(0);
    l->AddEntry(fSelEffEt [kNDmNuMuNCEL] , "#nu_{#mu} NCEL",    "LP");
    l->AddEntry(fSelEffEt [kNDmNuMuNCRES], "#nu_{#mu} NCRES",   "LP");
    l->AddEntry(fSelEffEt [kNDmNuMuNCCOH], "#nu_{#mu} NCCOH",   "LP");
    l->AddEntry(fSelEffEt [kNDmNuMuNCOth], "#nu_{#mu} NCOther", "LP");
    l->Draw();
    c->SetLogy();
  }
  else
  if(opt.Contains("Saclay/Efficiency/NuMuBar/CC")) { 
    TH1F * hframe = (TH1F*)c->DrawFrame(0,0.01,Emax,1);
    hframe->Draw();
    fSelEffEt [kNDmNuMuBarCCQE]     -> Draw("LPERRSAME");
    fSelEffEt [kNDmNuMuBarCCRES]    -> Draw("LPERRSAME");
    fSelEffEt [kNDmNuMuBarCCCOH]    -> Draw("LPERRSAME");
    fSelEffEt [kNDmNuMuBarCCOth]    -> Draw("LPERRSAME");
    fSelEffEt [kNDmNuMuBarCCQE]     -> GetXaxis()->SetTitle("E_{#nu} (GeV)");
    fSelEffEt [kNDmNuMuBarCCQE]     -> GetYaxis()->SetTitle("efficiency");
    TLegend * l = new TLegend(0.7,0.6,0.9,0.9);
    l->SetFillColor(0);
    l->AddEntry(fSelEffEt [kNDmNuMuBarCCQE] , "#bar{#nu_{#mu}} CCQE",    "LP");
    l->AddEntry(fSelEffEt [kNDmNuMuBarCCRES], "#bar{#nu_{#mu}} CCRES",   "LP");
    l->AddEntry(fSelEffEt [kNDmNuMuBarCCCOH], "#bar{#nu_{#mu}} CCCOH",   "LP");
    l->AddEntry(fSelEffEt [kNDmNuMuBarCCOth], "#bar{#nu_{#mu}} CCOther", "LP");
    l->Draw();
    c->SetLogy();
  }
  else
  if(opt.Contains("Saclay/Efficiency/NuMuBar/NC")) { 
    TH1F * hframe = (TH1F*)c->DrawFrame(0,0.01,Emax,1);
    hframe->Draw();
    fSelEffEt [kNDmNuMuBarNCEL]     -> Draw("LPERRSAME");
    fSelEffEt [kNDmNuMuBarNCRES]    -> Draw("LPERRSAME");
    fSelEffEt [kNDmNuMuBarNCCOH]    -> Draw("LPERRSAME");
    fSelEffEt [kNDmNuMuBarNCOth]    -> Draw("LPERRSAME");
    fSelEffEt [kNDmNuMuBarNCEL]     -> GetXaxis()->SetTitle("E_{#nu} (GeV)");
    fSelEffEt [kNDmNuMuBarNCEL]     -> GetYaxis()->SetTitle("efficiency");
    TLegend * l = new TLegend(0.7,0.6,0.9,0.9);
    l->SetFillColor(0);
    l->AddEntry(fSelEffEt [kNDmNuMuBarNCEL] , "#bar{#nu_{#mu}} NCEL",    "LP");
    l->AddEntry(fSelEffEt [kNDmNuMuBarNCRES], "#bar{#nu_{#mu}} NCRES",   "LP");
    l->AddEntry(fSelEffEt [kNDmNuMuBarNCCOH], "#bar{#nu_{#mu}} NCCOH",   "LP");
    l->AddEntry(fSelEffEt [kNDmNuMuBarNCOth], "#bar{#nu_{#mu}} NCOther", "LP");
    l->Draw();
    c->SetLogy();
  }
  else
  if(opt.Contains("Saclay/Efficiency/NuE/CC")) { 
    TH1F * hframe = (TH1F*)c->DrawFrame(0,0.01,Emax,1);
    hframe->Draw();
    fSelEffEt [kNDmNuECCQE]     -> Draw("LPERRSAME");
    fSelEffEt [kNDmNuECCRES]    -> Draw("LPERRSAME");
    fSelEffEt [kNDmNuECCCOH]    -> Draw("LPERRSAME");
    fSelEffEt [kNDmNuECCOth]    -> Draw("LPERRSAME");
    fSelEffEt [kNDmNuECCQE]     -> GetXaxis()->SetTitle("E_{#nu} (GeV)");
    fSelEffEt [kNDmNuECCQE]     -> GetYaxis()->SetTitle("efficiency");
    TLegend * l = new TLegend(0.7,0.6,0.9,0.9);
    l->SetFillColor(0);
    l->AddEntry(fSelEffEt [kNDmNuECCQE] , "#nu_{e} CCQE",    "LP");
    l->AddEntry(fSelEffEt [kNDmNuECCRES], "#nu_{e} CCRES",   "LP");
    l->AddEntry(fSelEffEt [kNDmNuECCCOH], "#nu_{e} CCCOH",   "LP");
    l->AddEntry(fSelEffEt [kNDmNuECCOth], "#nu_{e} CCOther", "LP");
    l->Draw();
    c->SetLogy();
  }
  else
  if(opt.Contains("Saclay/Efficiency/NuE/NC")) { 
    TH1F * hframe = (TH1F*)c->DrawFrame(0,0.01,Emax,1);
    hframe->Draw();
    fSelEffEt [kNDmNuENCEL]     -> Draw("LPERRSAME");
    fSelEffEt [kNDmNuENCRES]    -> Draw("LPERRSAME");
    fSelEffEt [kNDmNuENCCOH]    -> Draw("LPERRSAME");
    fSelEffEt [kNDmNuENCOth]    -> Draw("LPERRSAME");
    fSelEffEt [kNDmNuENCEL]     -> GetXaxis()->SetTitle("E_{#nu} (GeV)");
    fSelEffEt [kNDmNuENCEL]     -> GetYaxis()->SetTitle("efficiency");
    TLegend * l = new TLegend(0.7,0.6,0.9,0.9);
    l->SetFillColor(0);
    l->AddEntry(fSelEffEt [kNDmNuENCEL] , "#nu_{e} NCEL",    "LP");
    l->AddEntry(fSelEffEt [kNDmNuENCRES], "#nu_{e} NCRES",   "LP");
    l->AddEntry(fSelEffEt [kNDmNuENCCOH], "#nu_{e} NCCOH",   "LP");
    l->AddEntry(fSelEffEt [kNDmNuENCOth], "#nu_{e} NCOther", "LP");
    l->Draw();
    c->SetLogy();
  }
	 
  c->Update();
}
//________________________________________________________________________________________
double SuperkExpectation::ND280Inputs::CorrectionFactor(void) const
{
  //
  // avoid looping over all histograms and histogram bins if not needed
  //
  bool xsec_tweaked = fGenInputs->TweakAffectsND280();
  bool eff_tweaked  = false;
  for(int i = 0; i < knmnd; i++) {
    if( TMath::Abs(fSelEffTwkDial[i]) > kASmallNum ) eff_tweaked = true;
  }
  bool tweaked = eff_tweaked || xsec_tweaked;
  if(!tweaked) {
#ifdef __debug_level_2__
    cout << "Exp / nd280 >> nd280 event spectra not tweaked."
         << "Returning nominal value = " << fNomCorr << endl;
#endif
    return fNomCorr;
  }

  //
  // calculate correction factor
  //
  double nobscclike_nd280  = this -> NObsCCLike();
  double nexpncbkg_nd280   = this -> NExpNCBkg();
  double nexpcc_nd280      = this -> NExpCC();

  double correction = (nobscclike_nd280 - nexpncbkg_nd280) / nexpcc_nd280;

#ifdef __debug_level_2__
  cout << "Exp / nd280 >> "
       << "nobscclike_nd280 = "  << nobscclike_nd280
       << ", nexpncbkg_nd280 = " << nexpncbkg_nd280
       << ", nexpcc_nd280 = " << nexpcc_nd280
       << " ==> correction factor = " << correction << endl;
#endif

  return correction;
}
//________________________________________________________________________________________
double SuperkExpectation::ND280Inputs::NObsCCLike(void) const 
{
// observed number of CC-like events

  return fNumObsCCLike;
}
//________________________________________________________________________________________
double SuperkExpectation::ND280Inputs::NExpNCBkg(bool nominal) const
{
// expected NC contamination of CC-like sample

  double nexpncbkg = 
       this->NumOfEvents( kNDmNuMuNCEL,     nominal) +
       this->NumOfEvents( kNDmNuMuNCRES,    nominal) +
       this->NumOfEvents( kNDmNuMuNCCOH,    nominal) +
       this->NumOfEvents( kNDmNuMuNCOth,    nominal) +
       this->NumOfEvents( kNDmNuMuBarNCEL,  nominal) +
       this->NumOfEvents( kNDmNuMuBarNCRES, nominal) +
       this->NumOfEvents( kNDmNuMuBarNCCOH, nominal) +
       this->NumOfEvents( kNDmNuMuBarNCOth, nominal) +
       this->NumOfEvents( kNDmNuENCEL,      nominal) +
       this->NumOfEvents( kNDmNuENCRES,     nominal) +
       this->NumOfEvents( kNDmNuENCCOH,     nominal) +
       this->NumOfEvents( kNDmNuENCOth,     nominal);

  return nexpncbkg;
}
//________________________________________________________________________________________
double SuperkExpectation::ND280Inputs::NExpCC(bool nominal) const
{
// expected number of CC events

  double nexpcc = 
       this->NumOfEvents( kNDmNuMuCCQE,     nominal) +
       this->NumOfEvents( kNDmNuMuCCRES,    nominal) +
       this->NumOfEvents( kNDmNuMuCCCOH,    nominal) +
       this->NumOfEvents( kNDmNuMuCCOth,    nominal) +
       this->NumOfEvents( kNDmNuMuBarCCQE,  nominal) +
       this->NumOfEvents( kNDmNuMuBarCCRES, nominal) +
       this->NumOfEvents( kNDmNuMuBarCCCOH, nominal) +
       this->NumOfEvents( kNDmNuMuBarCCOth, nominal) +
       this->NumOfEvents( kNDmNuECCQE,      nominal) +
       this->NumOfEvents( kNDmNuECCRES,     nominal) +
       this->NumOfEvents( kNDmNuECCCOH,     nominal) +
       this->NumOfEvents( kNDmNuECCOth,     nominal);

  return nexpcc;
}
//________________________________________________________________________________________
double SuperkExpectation::ND280Inputs::NumOfEvents(ND280Mode_t mode, bool nominal) const
{
  TH1F * spectrum = fNumSelEt[mode];
  assert(spectrum);
  double sum_nev = 0;
  for(int ibin = 1; ibin <= spectrum->GetNbinsX(); ibin++) {
     sum_nev += this->NumOfEvents(mode,ibin, nominal);
  }
  return sum_nev;
}
//________________________________________________________________________________________
double SuperkExpectation::ND280Inputs::NumOfEvents(
       t2kana::yr1numu::conv::ND280Mode_t mode, int ibin, bool nominal) const
{
  TH1F * spectrum = fNumSelEt[mode];
  assert(spectrum);
  double mev2gev = 1E-3;
  double Etrue = spectrum->GetBinCenter(ibin) * mev2gev; // GeV

  double nev = spectrum->GetBinContent(ibin);
  if(nominal) return nev;

//double eff        = efficiency->GetBinContent(ibin);
  double efftwkdial = fSelEffTwkDial[mode];
  double efferr     = this->EfficiencyRelErr(Etrue, mode, util::Sign(efftwkdial));
//double twkeff     = eff * (1 + efftwkdial * efferr);
  double effwght    = 1 + efftwkdial * efferr;
  double sigwght    = fGenInputs->ND280XSecRatioWght(Etrue, mode);

  nev *= (effwght*sigwght);
  return nev;
}
//________________________________________________________________________________________
double SuperkExpectation::ND280Inputs::NumOfEvents(
       t2kana::yr1numu::conv::ND280Mode_t mode, double Etrue, bool nominal) const
{
  TH1F * spectrum   = fNumSelEt[mode];
  assert(spectrum);
  double mev2gev = 1E-3;
  int ibin = spectrum->FindBin(Etrue/mev2gev);
  return this->NumOfEvents(mode,ibin, nominal);
}
//________________________________________________________________________________________
void SuperkExpectation::ND280Inputs::UpdateEtrueSpectra(void)
{
  TH1D * spectrum = 0;

  // update numu
  spectrum = fNuMuEt;
  assert(spectrum);
  for(int ibin = 1; ibin <= spectrum->GetNbinsX(); ibin++) {
    double nev = 
         this->NumOfEvents( kNDmNuMuCCQE,  ibin) +
         this->NumOfEvents( kNDmNuMuCCRES, ibin) +
         this->NumOfEvents( kNDmNuMuCCCOH, ibin) +
         this->NumOfEvents( kNDmNuMuCCOth, ibin) +
         this->NumOfEvents( kNDmNuMuNCEL,  ibin) +
         this->NumOfEvents( kNDmNuMuNCRES, ibin) +
         this->NumOfEvents( kNDmNuMuNCCOH, ibin) +
         this->NumOfEvents( kNDmNuMuNCOth, ibin);
    spectrum->SetBinContent(ibin,nev);
  }//bins

  // update numubar
  spectrum = fNuMuBarEt;
  assert(spectrum);
  for(int ibin = 1; ibin <= spectrum->GetNbinsX(); ibin++) {
    double nev = 
         this->NumOfEvents( kNDmNuMuBarCCQE,  ibin) +
         this->NumOfEvents( kNDmNuMuBarCCRES, ibin) +
         this->NumOfEvents( kNDmNuMuBarCCCOH, ibin) +
         this->NumOfEvents( kNDmNuMuBarCCOth, ibin) +
         this->NumOfEvents( kNDmNuMuBarNCEL,  ibin) +
         this->NumOfEvents( kNDmNuMuBarNCRES, ibin) +
         this->NumOfEvents( kNDmNuMuBarNCCOH, ibin) +
         this->NumOfEvents( kNDmNuMuBarNCOth, ibin);
    spectrum->SetBinContent(ibin,nev);
  }//bins

  // update nue
  spectrum = fNuEEt;
  assert(spectrum);
  for(int ibin = 1; ibin <= spectrum->GetNbinsX(); ibin++) {
    double nev = 
         this->NumOfEvents( kNDmNuECCQE,  ibin) +
         this->NumOfEvents( kNDmNuECCRES, ibin) +
         this->NumOfEvents( kNDmNuECCCOH, ibin) +
         this->NumOfEvents( kNDmNuECCOth, ibin) +
         this->NumOfEvents( kNDmNuENCEL,  ibin) +
         this->NumOfEvents( kNDmNuENCRES, ibin) +
         this->NumOfEvents( kNDmNuENCCOH, ibin) +
         this->NumOfEvents( kNDmNuENCOth, ibin);
    spectrum->SetBinContent(ibin,nev);
  }//bins

}
//________________________________________________________________________________________
double SuperkExpectation::ND280Inputs::EfficiencyRelErr(
   double /*Et*/, ND280Mode_t /*mode*/, int /*sign*/) const
{
// relative uncertainty on efficiency 
// (can be asymmetric)

  return fSelEffRelErr;
}
//________________________________________________________________________________________
void SuperkExpectation::ND280Inputs::TweakEfficiency(ND280Mode_t mode, double twkdial)
{
// tweak efficiency for a particular mode
//
  fSelEffTwkDial[mode] = twkdial;
}
//________________________________________________________________________________________
void SuperkExpectation::ND280Inputs::TweakEfficiency(double twkdial)
{
// tweak efficiency of all modes by the same amount
//
  for(int i = 0; i < knmnd; i++) {
     fSelEffTwkDial[i] = twkdial;
  }
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
  fPOT = 0;  

  fMCErecoSpectrumGrandTot = 0;
  for(int is=0; is<kns; is++) {
    fMCErecoSpectrumTot[is] = 0;
    for(int im=0; im<knm; im++) {
      fMCTemplates     [is][im] = 0;
      fMCEtrueSpectrum [is][im] = 0;
      fMCErecoSpectrum [is][im] = 0; 
      fEffRelErr       [is][im] = 0;
      fEffTwkDial      [is][im] = 0;
    }
  }

  fEnergyScaleRelErr  = 0.;
  fEnergyScaleTwkDial = 0.;

  fMCPOT_numu     = 0;
  fMC1Rmu_numu    = 0;
  fMCPOT_numubar  = 0;
  fMC1Rmu_numubar = 0;
  fMCPOT_nue      = 0;
  fMC1Rmu_nue     = 0;
  fBrEreco        = 0;
  fBrEtrue        = 0;
  fBrPdg          = 0;
  fBrGenCode      = 0;

  fBeamInputs     = 0;
  fGenInputs      = 0;

  f1RmNuMuEt      = 0;
  f1RmNuMuBarEt   = 0;
  f1RmNuEEt       = 0;
}
//________________________________________________________________________________________
void SuperkExpectation::SuperkInputs::CleanUp(void)
{
  if(fMC1Rmu_numu   ) delete fMC1Rmu_numu;
  if(fMC1Rmu_numubar) delete fMC1Rmu_numubar;
  if(fMC1Rmu_nue    ) delete fMC1Rmu_nue;

  if( fMCErecoSpectrumGrandTot ) delete fMCErecoSpectrumGrandTot;
  for(int is=0; is<kns; is++) {
    if( fMCErecoSpectrumTot[is] ) delete fMCErecoSpectrumTot[is];
    for(int im=0; im<knm; im++) {
      if( fMCTemplates[is][im]     ) delete fMCTemplates[is][im];
      if( fMCEtrueSpectrum[is][im] ) delete fMCEtrueSpectrum[is][im];
      if( fMCErecoSpectrum[is][im] ) delete fMCErecoSpectrum[is][im];
    }
  }

  if (f1RmNuMuEt   ) delete f1RmNuMuEt;
  if (f1RmNuMuBarEt) delete f1RmNuMuBarEt;
  if (f1RmNuEEt    ) delete f1RmNuEEt;  
}
//________________________________________________________________________________________
bool SuperkExpectation::SuperkInputs::Read(TXMLEngine* xml, XMLNodePointer_t root_node)
{
// Read SuperK inputs to disappearance analysis
//
  cout << "Exp/SK >> Reading SuperK inputs from input XML file." << endl;

  //
  // Get data POT
  //
  fPOT = util::XML2Dbl(xml,root_node,"superk/observation/pot");
  cout << "Exp/SK >> " << "N_{POT}|data : " << fPOT << endl;

  //
  // Get MC trees and MC POT
  //
  fMCPOT_numu = util::XML2Dbl(xml,root_node,"superk/mc/numu/pot");
  fMC1Rmu_numu = util::XML2TTree(xml,root_node,"superk/mc/numu/one_mu_ring");
  cout << "Exp/SK >> " << "N_{POT}|(numu MC) : " << fMCPOT_numu << endl;
  assert(fMC1Rmu_numu);
  fMCPOT_numubar = util::XML2Dbl(xml,root_node,"superk/mc/numubar/pot");
  fMC1Rmu_numubar = util::XML2TTree(xml,root_node,"superk/mc/numubar/one_mu_ring");
  cout << "Exp/SK >> " << "N_{POT}|(numubar MC) : " << fMCPOT_numubar << endl;
  assert(fMC1Rmu_numubar);
  fMCPOT_nue = util::XML2Dbl(xml,root_node,"superk/mc/nuebeam/pot");
  fMC1Rmu_nue = util::XML2TTree(xml,root_node,"superk/mc/nuebeam/one_mu_ring");
  cout << "Exp/SK >> " << "N_{POT}|(beam nue MC) : " << fMCPOT_nue << endl;
  assert(fMC1Rmu_nue);

  // create templates
  this->CreateTemplate (kSKsFCFV1RingMuLike, kSKmNuMuCCQE );
  this->CreateTemplate (kSKsFCFV1RingMuLike, kSKmNuMuCC1pi);
  this->CreateTemplate (kSKsFCFV1RingMuLike, kSKmNuMuCCoth);
  this->CreateTemplate (kSKsFCFV1RingMuLike, kSKmNC       );
  this->CreateTemplate (kSKsFCFV1RingMuLike, kSKmNuMuBar  );
  this->CreateTemplate (kSKsFCFV1RingMuLike, kSKmNuE      );

  //
  // Systematics
  //

  // Get relative error on energy scale
  fEnergyScaleRelErr = util::XML2Dbl(
     xml,root_node,"superk/systematics/energy_scale_relative_error");
  // Get combined relative error on efficiency 
  fEffRelErr[kSKsFCFV1RingMuLike][kSKmNuMuCCQE] =  
    util::XML2Dbl(xml,root_node,
      "superk/systematics/efficiency_relative_error/one_mu_ring/ccqe");
  fEffRelErr[kSKsFCFV1RingMuLike][kSKmNuMuCC1pi] = 
    util::XML2Dbl(xml,root_node,
      "superk/systematics/efficiency_relative_error/one_mu_ring/cc1pi");
  fEffRelErr[kSKsFCFV1RingMuLike][kSKmNuMuCCoth] = 
    util::XML2Dbl(xml,root_node,
      "superk/systematics/efficiency_relative_error/one_mu_ring/ccoth");
  fEffRelErr[kSKsFCFV1RingMuLike][kSKmNC] =  
    util::XML2Dbl(xml,root_node,
      "superk/systematics/efficiency_relative_error/one_mu_ring/nc");

  cout << "Exp/SK >> " << "delta{energy_scale} : " << fEnergyScaleRelErr << endl;
  cout << "Exp/SK >> " << "delta{efficiency/1RingMu/CCQE} : " 
       << fEffRelErr[kSKsFCFV1RingMuLike][kSKmNuMuCCQE] << endl;
  cout << "Exp/SK >> " << "delta{efficiency/1RingMu/CC1pi} : " 
       << fEffRelErr[kSKsFCFV1RingMuLike][kSKmNuMuCC1pi] << endl;
  cout << "Exp/SK >> " << "delta{efficiency/1RingMu/CCoth} : " 
       << fEffRelErr[kSKsFCFV1RingMuLike][kSKmNuMuCCoth] << endl;
  cout << "Exp/SK >> " << "delta{efficiency/1RingMu/NC} : " 
       << fEffRelErr[kSKsFCFV1RingMuLike][kSKmNC] << endl;

  // Reset osc params
  this->UseOscillationParams(0., 0.);

  //
  // Create histograms
  //

  // Etrue histograms needed for beam norm err calculation
  f1RmNuMuEt    = new TH1D("f1RmNuMuEt",    "Etrue for 1-ring mu-like numu",    200, 0, 10); /// binning expected by Kubo's code
  f1RmNuMuBarEt = new TH1D("f1RmNuMuBarEt", "Etrue for 1-ring mu-like numubar", 200, 0, 10);
  f1RmNuEEt     = new TH1D("f1RmNuMuEEt",   "Etrue for 1-ring mu-like nue",     200, 0, 10);
  f1RmNuMuEt     -> SetDirectory(0);
  f1RmNuMuBarEt  -> SetDirectory(0);
  f1RmNuEEt      -> SetDirectory(0);

  // Ereco histograms filled-in by SK spectrum prediction code
  fMCErecoSpectrumGrandTot                              = new TH1D("","", kNvErecoBins, kNvErecoBinEdges);
  fMCErecoSpectrumTot[kSKsFCFV1RingMuLike]              = new TH1D("","", kNvErecoBins, kNvErecoBinEdges);
  fMCErecoSpectrum[kSKsFCFV1RingMuLike][kSKmNuMuCCQE]   = new TH1D("","", kNvErecoBins, kNvErecoBinEdges);
  fMCErecoSpectrum[kSKsFCFV1RingMuLike][kSKmNuMuCC1pi]  = new TH1D("","", kNvErecoBins, kNvErecoBinEdges);
  fMCErecoSpectrum[kSKsFCFV1RingMuLike][kSKmNuMuCCoth]  = new TH1D("","", kNvErecoBins, kNvErecoBinEdges);
  fMCErecoSpectrum[kSKsFCFV1RingMuLike][kSKmNC]         = new TH1D("","", kNvErecoBins, kNvErecoBinEdges);
  fMCErecoSpectrum[kSKsFCFV1RingMuLike][kSKmNuMuBar]    = new TH1D("","", kNvErecoBins, kNvErecoBinEdges);
  fMCErecoSpectrum[kSKsFCFV1RingMuLike][kSKmNuE]        = new TH1D("","", kNvErecoBins, kNvErecoBinEdges);
  fMCErecoSpectrumTot[kSKsFCFV1RingMuLike]              -> SetDirectory(0);
  fMCErecoSpectrum[kSKsFCFV1RingMuLike][kSKmNuMuCCQE]   -> SetDirectory(0);
  fMCErecoSpectrum[kSKsFCFV1RingMuLike][kSKmNuMuCC1pi]  -> SetDirectory(0);
  fMCErecoSpectrum[kSKsFCFV1RingMuLike][kSKmNuMuCCoth]  -> SetDirectory(0);
  fMCErecoSpectrum[kSKsFCFV1RingMuLike][kSKmNC]         -> SetDirectory(0);
  fMCErecoSpectrum[kSKsFCFV1RingMuLike][kSKmNuMuBar]    -> SetDirectory(0);
  fMCErecoSpectrum[kSKsFCFV1RingMuLike][kSKmNuE]        -> SetDirectory(0);

  cout << "Exp/SK >> Done reading SuperK inputs!" << endl;

  return true;
}
//________________________________________________________________________________________
void SuperkExpectation::SuperkInputs::CreateTemplate(SuperkSample_t s, SuperkMode_t m)
{
  cout << "Exp/SK >> Filling event Ereco-Etrue histogram for : " 
       << util::AsString(s) << "/" << util::AsString(m) << endl;

  int is = (int)s;
  int im = (int)m;

  fMCTemplates[is][im] = new TH2D(
     Form("nsk_%d_%d",is,im), util::AsString(m).c_str(), 
       kNvErecoBins, kNvErecoBinEdges, kNEtrueBins, kEtrueMin, kEtrueMax);
  fMCTemplates[is][im]->SetDirectory(0);
  fMCTemplates[is][im]->GetXaxis()->SetTitle("E_{reco} (GeV)");
  fMCTemplates[is][im]->GetYaxis()->SetTitle("E_{true} (GeV)");

  fMCEtrueSpectrum[is][im] = new TH1D(
     Form("nsketrue_%d_%d",is,im), "", kNEtrueBins, kEtrueMin, kEtrueMax);
  fMCEtrueSpectrum[is][im]->SetDirectory(0);
  fMCEtrueSpectrum[is][im]->GetXaxis()->SetTitle("E_{true} (GeV)");

  TTree * mctree = this->MCTree(s,m);
  double  mcpot  = this->MCPOT(s,m);

  double wght = this->AnalyzedSamplePOT() / mcpot;

  cout << "Exp/SK >> Weight (data POT / MC POT) = " << wght << endl;

  TBranch * branchEreco   = mctree->GetBranch("EnuReco");
  TBranch * branchEtrue   = mctree->GetBranch("Enu");
  TBranch * branchPdg     = mctree->GetBranch("Pdg");
  TBranch * branchGenCode = mctree->GetBranch("Mode");
  branchEreco   -> SetAddress(&fBrEreco);
  branchEtrue   -> SetAddress(&fBrEtrue); 
  branchPdg     -> SetAddress(&fBrPdg);
  branchGenCode -> SetAddress(&fBrGenCode);

  for(int iev = 0; iev < mctree->GetEntries(); iev++) 
  {
          mctree->GetEntry(iev);

          bool use_it = false;

          if(m == kSKmNuMuCCQE) {
             bool is_numu = (fBrPdg == kPdgNuMu);
             bool is_ccqe = (fBrGenCode == 1);
             use_it = is_numu && is_ccqe;
          } 
          else
          if(m == kSKmNuMuCC1pi) {
             bool is_numu  = (fBrPdg == kPdgNuMu);
             bool is_cc1pi = (fBrGenCode == 11 || fBrGenCode == 12 || fBrGenCode == 13);
            use_it = is_numu && is_cc1pi;
          } 
          else
          if(m == kSKmNuMuCCoth) {
             bool is_numu  = (fBrPdg == kPdgNuMu);
             bool is_ccoth = (fBrGenCode == 16 || 
                              fBrGenCode == 17 || 
                              fBrGenCode == 21 ||
                              fBrGenCode == 22 || 
                              fBrGenCode == 23 || 
                              fBrGenCode == 26);
             use_it = is_numu && is_ccoth;
          } 
          else
          if(m == kSKmNC) {
             bool is_numu  = (fBrPdg == kPdgNuMu);
             bool is_nc    = (fBrGenCode == 31 ||
                              fBrGenCode == 32 ||
                              fBrGenCode == 36 ||
                              fBrGenCode == 33 ||
                              fBrGenCode == 34 ||
                              fBrGenCode == 38 ||
                              fBrGenCode == 39 ||
                              fBrGenCode == 41 ||
                              fBrGenCode == 42 ||
                              fBrGenCode == 43 ||
                              fBrGenCode == 46 ||
                              fBrGenCode == 51 ||
                              fBrGenCode == 52);
             use_it = is_numu && is_nc;
          } 
          else
          if(m == kSKmNuMuBar) {
             bool is_numubar  = (fBrPdg == kPdgNuMuBar);
             use_it = is_numubar;
          } 
          else
          if(m == kSKmNuE) {
             bool is_nue  = (fBrPdg == kPdgNuE || fBrPdg == kPdgNuEBar);
             use_it = is_nue;
          }
          if(!use_it) continue;

        //cout << "Exp/SK >> Ereco = " << fBrEreco << ", Etrue = " << fBrEtrue << ", wght = " << wght << endl;

          fMCTemplates    [is][im] -> Fill (fBrEreco, fBrEtrue, wght);
          fMCEtrueSpectrum[is][im] -> Fill (fBrEtrue, wght);
  }//iev
}
//________________________________________________________________________________________
void SuperkExpectation::SuperkInputs::Draw(TVirtualPad * c, Option_t * option) const
{
  if(!c) return;

  c->cd();

  TString opt = option;
  util::SetDefaultStyle();

  if(opt.Contains("EtrueSpectra/NuMu")) { 
      fMC1Rmu_numu->Draw("Enu");
      fMC1Rmu_numu->Draw("Enu","Mode==1","same");
      fMC1Rmu_numu->Draw("Enu","Mode==11||Mode==12||Mode==13","same");
      c->Update();
  }

  if(opt.Contains("1RingMuLike/ErecoEtrueTemplates")) { 
    const char * draw_opt = "BOX";
    c->Divide(2,3);
    c->cd(1);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNuMuCCQE]->Draw(draw_opt);
    c->cd(2);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNuMuCC1pi]->Draw(draw_opt);
    c->cd(3);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNuMuCCoth]->Draw(draw_opt);
    c->cd(4);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNC]->Draw(draw_opt);
    c->cd(5);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNuMuBar]->Draw(draw_opt);
    c->cd(6);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNuE]->Draw(draw_opt);
    c->Update();
  }

  if(opt.Contains("1RingMuLike/EtrueSpectra")) { 
    const char * draw_opt = "";
    c->Divide(2,3);
    c->cd(1);
    fMCEtrueSpectrum[kSKsFCFV1RingMuLike][kSKmNuMuCCQE]->Draw(draw_opt);
    c->cd(2);
    fMCEtrueSpectrum[kSKsFCFV1RingMuLike][kSKmNuMuCC1pi]->Draw(draw_opt);
    c->cd(3);
    fMCEtrueSpectrum[kSKsFCFV1RingMuLike][kSKmNuMuCCoth]->Draw(draw_opt);
    c->cd(4);
    fMCEtrueSpectrum[kSKsFCFV1RingMuLike][kSKmNC]->Draw(draw_opt);
    c->cd(5);
    fMCEtrueSpectrum[kSKsFCFV1RingMuLike][kSKmNuMuBar]->Draw(draw_opt);
    c->cd(6);
    fMCEtrueSpectrum[kSKsFCFV1RingMuLike][kSKmNuE]->Draw(draw_opt);
    c->Update();
  }

  c->Update();
}
//________________________________________________________________________________________
double SuperkExpectation::SuperkInputs::NumOfEvents_Ereco(
    double Ereco, SuperkSample_t s, SuperkMode_t m, bool ignore_energy_scale_twk)
{
// Count events
// Apply oscillations, apply flux 10a->10d weights, 
// apply SuperK efficiency tweaks (if any) and SuperK energy scale tweaks (if any),
// apply flux systematic (shape + norm) tweaks.

  TH2D * nsk2d = fMCTemplates[s][m];
  assert(nsk2d);

  double nevtot = 0;

  int reco_bin = nsk2d->GetXaxis()->FindBin(Ereco);
  for(int true_bin = 1; true_bin <= nsk2d->GetYaxis()->GetNbins(); true_bin++) {

      double Etrue     = nsk2d->GetYaxis()->GetBinCenter(true_bin);
      double nev       = nsk2d->GetBinContent(reco_bin, true_bin);
      double prob_surv = util::NuMuSurvProb2D(fCurrSinsq2Theta23, fCurrDm23sq, Etrue, m);

      //
      // Flux weights & systematics
      //
      double flux_10d_weight      = 1.; // 10a->10d weight
      double numu_flux_shape_wght = 1.; // systematic variation (flux shape)
      if(fBeamInputs) {
         switch(m) {
           case (kSKmNuMuCCQE ):
           case (kSKmNuMuCC1pi):
           case (kSKmNuMuCCoth):
           case (kSKmNC       ):
               flux_10d_weight      = fBeamInputs->SuperkNuMuFlux10aTo10dWght(Etrue);
               numu_flux_shape_wght = fBeamInputs->SuperkNuMuFluxShapeWeight (Etrue);
               break;
           case(kSKmNuMuBar):
               flux_10d_weight = fBeamInputs->SuperkNuMuBarFlux10aTo10dWght(Etrue);
               break;
           case(kSKmNuE):
               flux_10d_weight = fBeamInputs->SuperkNuEFlux10aTo10dWght(Etrue);
               break;
           default:
               return 0;
         }//m switch
      }

      //
      // Cross section systematics = f(Etrue)
      //
      double xsec_ccqe_wght  = 1.;
      double xsec_ratio_wght = 1.;
      if(fGenInputs) {
         xsec_ccqe_wght  = fGenInputs->SuperkCCQEXSecWght (Etrue, m);
         xsec_ratio_wght = fGenInputs->SuperkXSecRatioWght(Etrue, m);
      }

      nevtot += ( nev * 
                  prob_surv * 
                  flux_10d_weight * 
                  numu_flux_shape_wght * 
                  xsec_ccqe_wght *
                  xsec_ratio_wght
                );
  }//Etrue

  //
  // Cross section systematics = f(Ereco)
  //
  double nucl_effect_wght = 1;
  if(fGenInputs) {
     nucl_effect_wght = fGenInputs->SuperkNuclEffectWght(Ereco,m);
  }
  nevtot *= nucl_effect_wght;

  //
  // Plug-in "efficiency" uncertainty due to fiducial volume, 
  // reduction chain, ring-counting and PID
  //
  double twkdial = fEffTwkDial[s][m];
  double relerr  = fEffRelErr[s][m];
  double effwght = (1. + twkdial*relerr);
  nevtot *= effwght; 

  //
  // Plug-in energy scale systematic.
  // Let t be the energy scale tweaking dial and sigE the input energy scale err.
  // Estimate the energy shift of the bin edges.
  // If the upper edge of bin i moves higher by dE, and the bin size is EbinSz 
  // then bin i loses dE/EbinSz fraction of its events to bin i+1.
  // But it gains dE'/EbinSz from bin i-1, where dE' is the shift of the lower bin
  // boundary. Vice-versa if the shift is downwards to lower energies.

  if(!ignore_energy_scale_twk) {
     double escale = this->EnergyScale();
     if( TMath::Abs(escale-1) > kASmallNum) {
        double bin_width   = nsk2d->GetXaxis()->GetBinWidth  (reco_bin);
        double bin_max     = nsk2d->GetXaxis()->GetBinUpEdge (reco_bin);
        double bin_min     = nsk2d->GetXaxis()->GetBinLowEdge(reco_bin);
        if(escale > 1) {
           // Losing events to bin i+1 and gaining events from bin i-1
           double lose_to_higher  = (escale - 1) * bin_max / bin_width;
           double gain_from_lower = (escale - 1) * bin_min / bin_width;
           double nev_prev_bin    = 0;
           if(reco_bin > 1) {
             double Ereco_prev_bin = nsk2d->GetXaxis()->GetBinCenter(reco_bin-1);
             nev_prev_bin = this->NumOfEvents_Ereco(Ereco_prev_bin,s,m,true);
           }
           nevtot *= (1-lose_to_higher);
           nevtot += gain_from_lower * nev_prev_bin;
        }//scale>1
        if(escale < 1) {
           // Losing events to bin i-1 and gaining events from bin i+1
           double lose_to_lower    = (1 - escale) * bin_min / bin_width;
           double gain_from_higher = (1 - escale) * bin_max / bin_width;
           double nev_next_bin    = 0;
           if(reco_bin < nsk2d->GetXaxis()->GetNbins()) {
             double Ereco_next_bin = nsk2d->GetXaxis()->GetBinCenter(reco_bin+1);
             nev_next_bin = this->NumOfEvents_Ereco(Ereco_next_bin,s,m,true);
           }
           nevtot *= (1-lose_to_lower);
           nevtot += gain_from_higher * nev_next_bin;
        }//scale<1
     }//energy scale tweaked?
  }//ignore escale

  return nevtot;
}
//________________________________________________________________________________________
double SuperkExpectation::SuperkInputs::NumOfEvents_Etrue(
    double Etrue, SuperkSample_t s, SuperkMode_t m)
{
// Count events in Etrue bin.
// + apply oscillations
// + apply flux 10a->10d weights
// + apply SuperK efficiency tweaks
// + apply flux systematic (shape + norm) tweaks

  int is = (int)s;
  int im = (int)m;

  TH1D * nsk1d = fMCEtrueSpectrum[is][im];
  assert(nsk1d);

  int true_bin = nsk1d->GetXaxis()->FindBin(Etrue);
  double nev = nsk1d->GetBinContent(true_bin);

#ifdef __debug_level_2__
  cout << "Exp/SK/NumOfEvents_Etrue() >> Etrue = " << Etrue << ", bin = " << true_bin << endl;
  cout << "Exp/SK/NumOfEvents_Etrue() >> Nev [unweighted] = " << nev  << endl;
#endif

  // 
  // Oscillations
  //
  double prob_surv = util::NuMuSurvProb2D(fCurrSinsq2Theta23, fCurrDm23sq, Etrue, m);
  nev *= prob_surv;

#ifdef __debug_level_2__
  cout << "Exp/SK/NumOfEvents_Etrue() >> Nev [wght:osc] = " << nev  << endl;
#endif

  //
  // Flux weights & systematics
  //
  double flux_10d_weight      = 1.; // 10a->10d weight
  double numu_flux_shape_wght = 1.; // systematic variation (flux shape)
  if(fBeamInputs) {
         switch(m) {
           case (kSKmNuMuCCQE ):
           case (kSKmNuMuCC1pi):
           case (kSKmNuMuCCoth):
           case (kSKmNC       ):
               flux_10d_weight      = fBeamInputs->SuperkNuMuFlux10aTo10dWght(Etrue);
               numu_flux_shape_wght = fBeamInputs->SuperkNuMuFluxShapeWeight (Etrue);
               break;
           case(kSKmNuMuBar):
               flux_10d_weight = fBeamInputs->SuperkNuMuBarFlux10aTo10dWght(Etrue);
               break;
           case(kSKmNuE):
               flux_10d_weight = fBeamInputs->SuperkNuEFlux10aTo10dWght(Etrue);
               break;
           default:
               return 0;
         }//m switch
  }
  nev *= (flux_10d_weight*numu_flux_shape_wght);

#ifdef __debug_level_2__
  cout << "Exp/SK/NumOfEvents_Etrue() >> Nev [wght:osc+flux] = " << nev  << endl;
#endif

  //
  // Cross section systematics = f(Etrue)
  //
  double xsec_ccqe_wght  = 1.;
  double xsec_ratio_wght = 1.;
  if(fGenInputs) {
      xsec_ccqe_wght  = fGenInputs->SuperkCCQEXSecWght (Etrue, m);
      xsec_ratio_wght = fGenInputs->SuperkXSecRatioWght(Etrue, m);
  }
  nev *= (xsec_ccqe_wght*xsec_ratio_wght);

#ifdef __debug_level_2__
  cout << "Exp/SK/NumOfEvents_Etrue() >> Nev [wght:osc+flux+xsec] = " << nev  << endl;
#endif

  //
  // Plug-in "efficiency" uncertainty due to fiducial volume, 
  // reduction chain, ring-counting and PID
  //
  double twkdial = fEffTwkDial[s][m];
  double relerr  = fEffRelErr[s][m];
  double effwght = (1. + twkdial*relerr);
  nev *= effwght; 

#ifdef __debug_level_2__
  cout << "Exp/SK/NumOfEvents_Etrue() >> Nev [wght:osc+flux+xsec+eff] = " << nev  << endl;
#endif

  return nev;
}
//________________________________________________________________________________________
void SuperkExpectation::SuperkInputs::UpdateErecoSpectra(void)
{
  TH1D * reco_spectrum_grand_tot =  fMCErecoSpectrumGrandTot;
  assert(reco_spectrum_grand_tot);

  for(int ibin = 1; ibin  <= reco_spectrum_grand_tot->GetNbinsX(); ibin++) {
     double Ereco = reco_spectrum_grand_tot->GetBinCenter(ibin);

     double NevGrandTot = 0;    

     for(int is = 0; is < kns; is++) {
       TH1D * reco_spectrum_tot =  fMCErecoSpectrumTot[is];
       if(!reco_spectrum_tot) continue; /// loop only over samples included in osc fit

       double NevTot = 0;

       for(int im = 0; im < knm; im++) { 
            TH1D * reco_spectrum =  fMCErecoSpectrum[is][im];
            assert(reco_spectrum);
            SuperkSample_t s = (SuperkSample_t) is;
            SuperkMode_t   m = (SuperkMode_t  ) im;
            double Nev = this->NumOfEvents_Ereco(Ereco,s,m);
            reco_spectrum->SetBinContent(ibin,Nev);
            NevTot += Nev;
        }//im
        reco_spectrum_tot->SetBinContent(ibin,NevTot);
        NevGrandTot += NevTot;
     }//is

     reco_spectrum_grand_tot->SetBinContent(ibin,NevGrandTot);

  }//bin
}
//________________________________________________________________________________________
void SuperkExpectation::SuperkInputs::UpdateEtrueSpectra(void)
{
  cout << "Exp/SK >> Updating Etrue spectra for calculating beam norm err" << endl;

  TH1D * spectrum = 0;

  // update 1-ring mu-like numu
  spectrum = f1RmNuMuEt;
  for(int ibin = 1; ibin <= spectrum->GetNbinsX(); ibin++) {
    double Etrue = spectrum->GetBinCenter(ibin);
    double nev = 
       this->NumOfEvents_Etrue(Etrue, kSKsFCFV1RingMuLike, kSKmNuMuCCQE ) +
       this->NumOfEvents_Etrue(Etrue, kSKsFCFV1RingMuLike, kSKmNuMuCC1pi) +
       this->NumOfEvents_Etrue(Etrue, kSKsFCFV1RingMuLike, kSKmNuMuCCoth) +
       this->NumOfEvents_Etrue(Etrue, kSKsFCFV1RingMuLike, kSKmNC       );
#ifdef __debug_level_2__
    cout << "Exp/SK/UpdateEtrueSpectra() >> "
         << "1-ring mu / numu : bin = " << ibin << ", nev = " << nev << endl;
#endif
    spectrum->SetBinContent(ibin,nev);
  }//bins

  // update 1-ring mu-like numubar
  spectrum = f1RmNuMuBarEt;
  for(int ibin = 1; ibin <= spectrum->GetNbinsX(); ibin++) {
    double Etrue = spectrum->GetBinCenter(ibin);
    double nev = 
       this->NumOfEvents_Etrue(Etrue, kSKsFCFV1RingMuLike, kSKmNuMuBar);
#ifdef __debug_level_2__
    cout << "Exp/SK/UpdateEtrueSpectra() >> "
         << "1-ring mu / numubar : bin = " << ibin << ", nev = " << nev << endl;
#endif
    spectrum->SetBinContent(ibin,nev);
  }//bins

  // update 1-ring mu-like nue
  spectrum = f1RmNuEEt;
  for(int ibin = 1; ibin <= spectrum->GetNbinsX(); ibin++) {
    double Etrue = spectrum->GetBinCenter(ibin);
    double nev = 
       this->NumOfEvents_Etrue(Etrue, kSKsFCFV1RingMuLike, kSKmNuE);
#ifdef __debug_level_2__
    cout << "Exp/SK/UpdateEtrueSpectra() >> "
         << "1-ring mu / nuebkg : bin = " << ibin << ", nev = " << nev << endl;
#endif
    spectrum->SetBinContent(ibin,nev);
  }//bins
}
//________________________________________________________________________________________
double SuperkExpectation::SuperkInputs::NumOfEventsReadFromTree(
    double Ereco, double ErecoBinSz, SuperkSample_t s, SuperkMode_t m)
{
  double nev = 0;

  TTree * mctree = this->MCTree(s,m);
  double  mcpot  = this->MCPOT(s,m);  

  for(int iev = 0; iev < mctree->GetEntries(); iev++) {
    mctree->GetEntry(iev);

    if(fBrEreco <  Ereco - 0.5*ErecoBinSz) continue;
    if(fBrEreco >= Ereco + 0.5*ErecoBinSz) continue;

    bool use_it = false;

    if(m == kSKmNuMuCCQE) {
      bool is_numu = (fBrPdg == kPdgNuMu);
      bool is_ccqe = (fBrGenCode == 1);
      use_it = is_numu && is_ccqe;
    } 
    else
    if(m == kSKmNuMuCC1pi) {
      bool is_numu  = (fBrPdg == kPdgNuMu);
      bool is_cc1pi = (fBrGenCode == 11 || fBrGenCode == 12 || fBrGenCode == 13);
      use_it = is_numu && is_cc1pi;
    } 
    else
    if(m == kSKmNuMuCCoth) {
      bool is_numu  = (fBrPdg == kPdgNuMu);
      bool is_ccoth = (fBrGenCode == 16 || 
                       fBrGenCode == 17 || 
                       fBrGenCode == 21 ||
                       fBrGenCode == 22 || 
                       fBrGenCode == 23 || 
                       fBrGenCode == 26);
      use_it = is_numu && is_ccoth;
    } 
    else
    if(m == kSKmNC) {
      bool is_numu  = (fBrPdg == kPdgNuMu);
      bool is_nc    = (fBrGenCode == 31 ||
                       fBrGenCode == 32 ||
                       fBrGenCode == 36 ||
                       fBrGenCode == 33 ||
                       fBrGenCode == 34 ||
                       fBrGenCode == 38 ||
                       fBrGenCode == 39 ||
                       fBrGenCode == 41 ||
                       fBrGenCode == 42 ||
                       fBrGenCode == 43 ||
                       fBrGenCode == 46 ||
                       fBrGenCode == 51 ||
                       fBrGenCode == 52);
      use_it = is_numu && is_nc;
    } 
    else
    if(m == kSKmNuMuBar) {
      bool is_numubar  = (fBrPdg == kPdgNuMuBar);
      use_it = is_numubar;
    } 
    else
    if(m == kSKmNuE) {
      bool is_nue  = (fBrPdg == kPdgNuE || fBrPdg == kPdgNuEBar);
      use_it = is_nue;
    }
    if(!use_it) continue;

    double prob_surv = util::NuMuSurvProb2D(fCurrSinsq2Theta23, fCurrDm23sq, fBrEtrue, m);

    nev += prob_surv;
  }//iev

  nev *= (fPOT/mcpot); // MC POT -> DATA POT

  return nev;
}
//________________________________________________________________________________________
double SuperkExpectation::SuperkInputs::AnalyzedSamplePOT(void) const
{
  return fPOT;  
}
//________________________________________________________________________________________
TTree * SuperkExpectation::SuperkInputs::MCTree(SuperkSample_t s, SuperkMode_t m)
{
  assert(s==kSKsFCFV1RingMuLike);

  switch(m) {
     case (kSKmNuMuCCQE ):
     case (kSKmNuMuCC1pi):
     case (kSKmNuMuCCoth):
     case (kSKmNC       ):
        return fMC1Rmu_numu;
        break;
     case(kSKmNuMuBar):
        return fMC1Rmu_numubar;
        break;
     case(kSKmNuE):
        return fMC1Rmu_nue;
        break;
     default:
        return 0;
  }
  return 0;
}
//________________________________________________________________________________________
double SuperkExpectation::SuperkInputs::MCPOT(SuperkSample_t s, SuperkMode_t m)
{
  assert(s==kSKsFCFV1RingMuLike);

  switch(m) {
     case (kSKmNuMuCCQE ):
     case (kSKmNuMuCC1pi):
     case (kSKmNuMuCCoth):
     case (kSKmNC       ):
        return fMCPOT_numu;
        break;
     case(kSKmNuMuBar):
        return fMCPOT_numubar;
        break;
     case(kSKmNuE):
        return fMCPOT_nue;
        break;
     default:
        return 0;
  }
  return 0;
}
//________________________________________________________________________________________
double SuperkExpectation::SuperkInputs::EfficiencyRelErr(
    SuperkSample_t s, SuperkMode_t m) const
{
  double err = fEffRelErr[s][m];
  return err;
}
//________________________________________________________________________________________
double SuperkExpectation::SuperkInputs::EnergyScaleRelErr(void) const
{
// energy scale relative error

  return fEnergyScaleRelErr;
}
//________________________________________________________________________________________
double SuperkExpectation::SuperkInputs::EnergyScale(void) const
{
// energy scale

  double twkdial = fEnergyScaleTwkDial;
  double err     = this->EnergyScaleRelErr();
  double scale   = 1. + twkdial * err;
  return scale;
}
//________________________________________________________________________________________
void SuperkExpectation::SuperkInputs::TweakEnergyScale(double twkdial) 
{ 
  fEnergyScaleTwkDial = twkdial; 
}
//________________________________________________________________________________________
void SuperkExpectation::SuperkInputs::TweakEfficiency(
  SuperkSample_t s, SuperkMode_t m, double twkdial)
{
  fEffTwkDial[s][m] = twkdial;

  // link CC1pi and CCoth as I have a single efficiency for CCnonQE
  if      (m == kSKmNuMuCCoth) { fEffTwkDial[s][kSKmNuMuCC1pi] = twkdial; }
  else if (m == kSKmNuMuCC1pi) { fEffTwkDial[s][kSKmNuMuCCoth] = twkdial; }
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
  fSuperkNuMuFlux              = 0; 
  fSuperkNuMuFluxVarBin        = 0; 
  fSuperkNuMuFluxPOT           = 0;
  fSuperkNuMuFluxCov           = 0;
  fSuperkNuMuWght10aTo10d      = 0;
  fSuperkNuMuBarWght10aTo10d   = 0;
  fSuperkNuEWght10aTo10d       = 0;
  fSuperkNuMuFluxNormTwkDial   = 0;
  fSuperkNuMuFluxCovCholeskyL  = 0;
  fSuperkNuMuFluxHaveCholeskyL = false;
  fSuperkNuMuFluxTweaked       = 0;
  fEnableFluxShapeTweaking     = false;
  fHEC                         = 0;
  fSuperkNormRelErr_simple     = 0;
}
//________________________________________________________________________________________
void SuperkExpectation::BeamInputs::CleanUp(void)
{
  if ( fSuperkNuMuFlux             ) delete fSuperkNuMuFlux;
  if ( fSuperkNuMuFluxVarBin       ) delete fSuperkNuMuFluxVarBin;
  if ( fSuperkNuMuFluxCov          ) delete fSuperkNuMuFluxCov;
  if ( fSuperkNuMuWght10aTo10d     ) delete fSuperkNuMuWght10aTo10d;
  if ( fSuperkNuMuBarWght10aTo10d  ) delete fSuperkNuMuBarWght10aTo10d;
  if ( fSuperkNuEWght10aTo10d      ) delete fSuperkNuEWght10aTo10d;
  if ( fSuperkNuMuFluxTweaked      ) delete fSuperkNuMuFluxTweaked;
  if ( fSuperkNuMuFluxCovCholeskyL ) delete fSuperkNuMuFluxCovCholeskyL;
  if ( fHEC                        ) delete fHEC;
}
//________________________________________________________________________________________
bool SuperkExpectation::BeamInputs::Read(TXMLEngine* xml, XMLNodePointer_t root_node)
{
// Read beam inputs to disappearance analysis
//
  cout << "Exp/Beam >> Reading beam inputs from input XML file." << endl;

  // Get SuperK numu flux
  fSuperkNuMuFlux = util::XML2TH1D(xml,root_node,"beam/superk_flux/v10a/numu");
  if(!fSuperkNuMuFlux) {
    cerr << "Exp/Beam >> No input SuperK numu flux histogram!" << endl;
    exit(1);
  }
  // Get SuperK numu flux (var binning)
  fSuperkNuMuFluxVarBin = util::XML2TH1D(xml,root_node,"beam/superk_flux/v10a/numu_varbin");
  if(!fSuperkNuMuFluxVarBin) {
    cerr << "Exp/Beam >> No input SuperK numu (var binning) flux histogram!" << endl;
    exit(1);
  }
  // Integrated POT for numu flux simulation
  fSuperkNuMuFluxPOT = util::XML2Dbl(xml,root_node,"beam/superk_flux/pot_mc");
  cout << "Exp/Beam >> SK flux simulation POT = " << fSuperkNuMuFluxPOT << endl;

  // Get SuperK numu flux covariance matrix
  fSuperkNuMuFluxCov = util::XML2TMatrixD(
      xml,root_node,"beam/superk_flux/systematics/covariance/numu");
  if(!fSuperkNuMuFluxCov) {
    cerr << "Exp/Beam >> No input SuperK numu flux covariance matrix!" << endl;
    exit(1);
  }
#ifdef __debug_level_3__
  fSuperkNuMuFluxCov->Print();
#endif

  // 10a->10d weights
  fSuperkNuMuWght10aTo10d = 
     util::XML2TH1D(xml,root_node,"beam/superk_flux/weight_v10a_to_v10d/numu");
  if(!fSuperkNuMuWght10aTo10d) {
    cerr << "Exp/Beam >> No input 10a->10d SuperK numu flux weights!" << endl;
    exit(1);
  }
  fSuperkNuMuBarWght10aTo10d = 
     util::XML2TH1D(xml,root_node,"beam/superk_flux/weight_v10a_to_v10d/numubar");
  if(!fSuperkNuMuBarWght10aTo10d) {
    cerr << "Exp/Beam >> No input 10a->10d SuperK numubar flux weights!" << endl;
    exit(1);
  }
  fSuperkNuEWght10aTo10d = 
     util::XML2TH1D(xml,root_node,"beam/superk_flux/weight_v10a_to_v10d/nue");
  if(!fSuperkNuEWght10aTo10d) {
    cerr << "Exp/Beam >> No input 10a->10d SuperK nue flux weights!" << endl;
    exit(1);
  }

  cout << "Exp/Beam >> Loading Kubo-san's norm err inputs!" << endl;

  fSuperkNormRelErr_simple = 0;

  string kubo_data_path =
      util::XML2String(xml,root_node,"beam/superk_flux/systematics/norm_error_inputs");
  if(kubo_data_path != "none") {
     cout << "Exp/Beam >> Using input files from: " << kubo_data_path << endl;
     fHEC = new HadErrCalc(kubo_data_path.c_str());
  } else {
     cout << "Exp/Beam >> Not using Kubo's code / Use a fixed Nsk/Nd280 err instead." << endl;
     fSuperkNormRelErr_simple = util::XML2Dbl(
         xml,root_node,"beam/superk_flux/systematics/norm_error_simple_version");
     cout << "Exp/Beam >> Norm relative err (simple) : " << fSuperkNormRelErr_simple << endl;
  }

  cout << "Exp/Beam >> Done reading beam inputs!" << endl;

  return true;
}
//________________________________________________________________________________________
void SuperkExpectation::BeamInputs::Draw(TVirtualPad * c, Option_t * option) const
{
  if(!c) return;

  c->cd();

  TString opt = option;
  util::SetDefaultStyle();

  if(opt.Contains("SKFlux/Weights/10aTo10d/NuMu")) { 
    fSuperkNuMuWght10aTo10d->Draw();
    fSuperkNuMuWght10aTo10d->GetXaxis()->SetTitle("E_{#nu} (GeV)");
    fSuperkNuMuWght10aTo10d->GetYaxis()->SetTitle("SuperK numu flux weight");
    c->Update();
  }
  else
  if(opt.Contains("SKFlux/10a/NuMu")) { 
    fSuperkNuMuFlux->Draw();
    fSuperkNuMuFlux->GetXaxis()->SetTitle("E_{#nu} (GeV)");
    fSuperkNuMuFlux->GetYaxis()->SetTitle("SuperK numu flux (#neutrinos /cm2 /POT)");
    c->Update();
  } 
  else
  if(opt.Contains("SKFlux/10d/NuMu")) { 
    TH1D * weighted_flux = new TH1D(*fSuperkNuMuFlux);
    weighted_flux->SetDirectory(0);
    for(int bin=1; bin <= weighted_flux->GetNbinsX(); bin++) {
	double Et = weighted_flux->GetBinCenter(bin);
	double f  = this->SuperkNuMuFlux(Et,true);
        weighted_flux->SetBinContent(bin,f);
    }
    weighted_flux->Draw();
    weighted_flux->GetXaxis()->SetTitle("E_{#nu} (GeV)");
    weighted_flux->GetYaxis()->SetTitle("SuperK numu flux (#neutrinos /cm2 /POT)");
    c->Update();
  }
  else
  if(opt.Contains("SKFlux/Both10a10d/NuMu")) { 
    TH1D * weighted_flux = new TH1D(*fSuperkNuMuFlux);
    weighted_flux->SetDirectory(0);
    for(int bin=1; bin <= weighted_flux->GetNbinsX(); bin++) {
	double Et = weighted_flux->GetBinCenter(bin);
	double f  = this->SuperkNuMuFlux(Et,true);
        weighted_flux->SetBinContent(bin,f);
    }
    weighted_flux->SetLineStyle(kDashed);
    fSuperkNuMuFlux->Draw();
    weighted_flux->Draw("same");
    fSuperkNuMuFlux->GetXaxis()->SetTitle("E_{#nu} (GeV)");
    fSuperkNuMuFlux->GetYaxis()->SetTitle("SuperK numu flux (#neutrinos /cm2 /POT)");
    c->Update();
  }
  else
  if(opt.Contains("SKFlux/CovarianceMatrix/NuMu")) { 
    fSuperkNuMuFluxCov->Draw("colz");
    c->Update();
  }
}
//________________________________________________________________________________________
double SuperkExpectation::BeamInputs::SuperkNuMuFlux10aTo10dWght(double Etrue) const
{
  double wght = 1.0;
  if(Etrue < 10/*GeV*/) {
    int wght_bin = fSuperkNuMuWght10aTo10d->FindBin(Etrue);
    wght = fSuperkNuMuWght10aTo10d->GetBinContent(wght_bin);
  }
  return wght;
}
//________________________________________________________________________________________
double SuperkExpectation::BeamInputs::SuperkNuMuBarFlux10aTo10dWght(double Etrue) const
{
  double wght = 1.0;
  if(Etrue < 10/*GeV*/) {
    int wght_bin = fSuperkNuMuBarWght10aTo10d->FindBin(Etrue);
    wght = fSuperkNuMuBarWght10aTo10d->GetBinContent(wght_bin);
  }
  return wght;
}
//________________________________________________________________________________________
double SuperkExpectation::BeamInputs::SuperkNuEFlux10aTo10dWght(double Etrue) const
{
  double wght = 1.0;
  if(Etrue < 10/*GeV*/) {
    int wght_bin = fSuperkNuEWght10aTo10d->FindBin(Etrue);
    wght = fSuperkNuEWght10aTo10d->GetBinContent(wght_bin);
  }
  return wght;
}
//________________________________________________________________________________________
double SuperkExpectation::BeamInputs::SuperkNuMuFlux(double Etrue, bool wght_to_10d) const
{
  assert(fSuperkNuMuFlux);

  int flux_bin = fSuperkNuMuFlux->FindBin(Etrue);
  double flux  = fSuperkNuMuFlux->GetBinContent(flux_bin);
 
  if(wght_to_10d) {
    double wght = this->SuperkNuMuFlux10aTo10dWght(Etrue);
    flux *= wght;
  }

  return flux;
}
//________________________________________________________________________________________
double SuperkExpectation::BeamInputs::SuperkNuMuFluxPOT (void) const
{
  return fSuperkNuMuFluxPOT;
}
//________________________________________________________________________________________
const TMatrixD &  SuperkExpectation::BeamInputs::SuperkNuMuFluxCov(void) const 
{ 
  return *fSuperkNuMuFluxCov; 
} 
//________________________________________________________________________________________
double SuperkExpectation::BeamInputs::SuperkFluxNormErr(
  TH1D * hSkNuMu_1RingMu, TH1D * hSkNuMuBar_1RingMu, TH1D * hSkNuE_1RingMu,
  TH1D * hNdNuMu,         TH1D * hNdNuMuBar,         TH1D * hNdNuE) const
{
// normalization err Nsk/Nd280

  // don't  use Kubo-san's code. 
  if(!fHEC) {
    return fSuperkNormRelErr_simple;
  }

#ifdef __debug_level_3__
  TFile f("./inputs_to_beam_wg_norm_err_calc_code.root","recreate");
  TH1D * h1 = new TH1D(*hSkNuMu_1RingMu);
  TH1D * h2 = new TH1D(*hSkNuMuBar_1RingMu);
  TH1D * h3 = new TH1D(*hSkNuE_1RingMu);
  TH1D * h4 = new TH1D(*hNdNuMu);
  TH1D * h5 = new TH1D(*hNdNuMuBar);
  TH1D * h6 = new TH1D(*hNdNuE);
  h1 -> Write("sk1rmu_numu");
  h2 -> Write("sk1rmu_numubar");
  h3 -> Write("sk1rmu_nue");
  h4 -> Write("nd_numu");
  h5 -> Write("nd_numubar");
  h6 -> Write("nd_nue");
  f.Close();
#endif

  // use code contributed by Kubo-san's on behalf of the Beam WG
  fHEC->Set1RmuHist(hSkNuMu_1RingMu, hSkNuMuBar_1RingMu, hSkNuE_1RingMu);
  fHEC->SetNDHist  (hNdNuMu,         hNdNuMuBar,         hNdNuE);
  fHEC->Calc();

#ifdef __debug_level_3__
  fHEC->ErrorComponents();
  fHEC->ShowDetail();
#endif

  // get relative err on Nsk/Nnd280
  double relative_err = fHEC->GetSumErr(0);

  cout << "Exp/Beam >> Nsk/Nnd280 relative err calculated by Beam WG code = " 
       << relative_err << endl;

  return relative_err;
}
//________________________________________________________________________________________
double SuperkExpectation::BeamInputs::SuperkFluxNorm(
  TH1D * hSkNuMu_1RingMu, TH1D * hSkNuMuBar_1RingMu, TH1D * hSkNuE_1RingMu,
  TH1D * hNdNuMu,         TH1D * hNdNuMuBar,         TH1D * hNdNuE) const
{
// normalization

  double twkdial = fSuperkNuMuFluxNormTwkDial;

  if(TMath::Abs(twkdial) < kASmallNum) return 1; // untweaked - skip calculation

  double relerr  = this->SuperkFluxNormErr(
     hSkNuMu_1RingMu,hSkNuMuBar_1RingMu,hSkNuE_1RingMu,hNdNuMu,hNdNuMuBar,hNdNuE);

#ifdef __debug_level_3__
  cout << "Exp/Beam/SuperkFluxNorm() >> "
       << "Normalization relative error = " << relerr << endl;
#endif
	
  double norm = 1 + twkdial * relerr;
  return norm;
}
//________________________________________________________________________________________
double  SuperkExpectation::BeamInputs::SuperkNuMuFluxShapeWeight(double Etrue)
{
  if(!fEnableFluxShapeTweaking) return 1;

  if(!fSuperkNuMuFluxHaveCholeskyL) {
	exit(1);
  }
  if(!fSuperkNuMuFluxTweaked || !fSuperkNuMuFlux) {
	exit(1);
  }

  int bin = fSuperkNuMuFlux->FindBin(Etrue);

  double flux_nom = fSuperkNuMuFlux        ->GetBinContent(bin);
  double flux_twk = fSuperkNuMuFluxTweaked ->GetBinContent(bin);

  double wght = 1.;
  if(flux_nom > 0) {
	wght = flux_twk/flux_nom;
  }
  return wght;
}
//________________________________________________________________________________________
void SuperkExpectation::BeamInputs::TweakShape(bool enable_tweaking)    
{ 
  fEnableFluxShapeTweaking = enable_tweaking;
  if(!fEnableFluxShapeTweaking) return;

  if(!fSuperkNuMuFluxHaveCholeskyL) {
    cout << "Exp/Beam/TweakShape >> Running Cholesky's method." << endl;
    TMatrixD  sk_flux_cov(this->SuperkNuMuFluxCov());
    TMatrixD L = util::CholeskyDecomposition(sk_flux_cov);
#ifdef __debug_level_3__
    L.Print();
#endif

    fSuperkNuMuFluxCovCholeskyL = new TMatrixD(L);
#ifdef __debug_level_3__
    fSuperkNuMuFluxCovCholeskyL->Print();
#endif
    fSuperkNuMuFluxHaveCholeskyL = true;
    fSuperkNuMuFluxTweaked = new TH1D(*fSuperkNuMuFlux);
  }

#ifdef __debug_level_3__
  cout << "Exp/Beam >> Tweaking SuperK flux shape" << endl;
#endif
 
  //get the variations
  TVectorD tweek_vec = 
     util::CholeskyGenerateCorrelatedParamVariations(*fSuperkNuMuFluxCovCholeskyL);

  // fill in the histo
//int nbins = fSuperkNuMuFlux->GetXaxis()->GetNbins();
//for(int j=1; j <= nbins; j++){
  for(int j=1; j <= 200; j++){ ///TODO
     double rto = util::FluxFractionsFixedBinRelToVarBin(
         j, fSuperkNuMuFlux, fSuperkNuMuFluxVarBin);
     double fixBinFlux = fSuperkNuMuFlux->GetBinContent(j);
     int varBin = util::FluxFix2VarBinNumber(j);
     fSuperkNuMuFluxTweaked->SetBinContent(j, fixBinFlux + tweek_vec[varBin - 1]*rto );
  }//end of bin loop

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
  fCCQEXSecTwkDial             = 0;
  fCC1piToCCQEXSecRatioTwkDial = 0;
  fCCOthToCCQEXSecRatioTwkDial = 0;
  fNCToCCQEXSecRatioTwkDial    = 0;
  fNuclEffectTwkDial           = 0;
  fCCQEXSecRelErr              = 0;
  fCC1piToCCQEXSecRatioRelErr  = 0;
  fCCOthToCCQEXSecRatioRelErr  = 0;
  fNCToCCQEXSecRatioRelErr     = 0;
  fSuperkNuclearEffectRelErr   = 0;
}
//________________________________________________________________________________________
void SuperkExpectation::GeneratorInputs::CleanUp(void)
{
  if ( fCCQEXSecRelErr             ) delete fCCQEXSecRelErr;
  if ( fCC1piToCCQEXSecRatioRelErr ) delete fCC1piToCCQEXSecRatioRelErr;
  if ( fCCOthToCCQEXSecRatioRelErr ) delete fCCOthToCCQEXSecRatioRelErr;
  if ( fNCToCCQEXSecRatioRelErr    ) delete fNCToCCQEXSecRatioRelErr;
  if ( fSuperkNuclearEffectRelErr  ) delete fSuperkNuclearEffectRelErr;
}
//________________________________________________________________________________________
bool SuperkExpectation::GeneratorInputs::Read(TXMLEngine* xml, XMLNodePointer_t root_node)
{
  cout << "Exp/NuGen >> Reading neutrino generator inputs from input XML file." << endl;

  fCCQEXSecRelErr = 
    util::XML2TH1D(xml,root_node,"niwg/systematics/xsec_relative_error/ccqe");
  fCC1piToCCQEXSecRatioRelErr = 
    util::XML2TH1D(xml,root_node,"niwg/systematics/xsec_ratio_relative_error/cc1pi_ccqe");
  fCCOthToCCQEXSecRatioRelErr = 
    util::XML2TH1D(xml,root_node,"niwg/systematics/xsec_ratio_relative_error/ccother_ccqe");
  fNCToCCQEXSecRatioRelErr = 
    util::XML2TH1D(xml,root_node,"niwg/systematics/xsec_ratio_relative_error/nc_ccqe");
  fSuperkNuclearEffectRelErr =
    util::XML2TH1D(xml,root_node,"niwg/systematics/superk_nuclear_effect_relative_error");

  cout << "Exp/NuGen >> Done reading neutrino generator inputs!" << endl;

  return true;  
}
//________________________________________________________________________________________
void SuperkExpectation::GeneratorInputs::Draw(TVirtualPad * c, Option_t * option) const
{
  if(!c) return;

  c->cd();

  TString opt = option;
  util::SetDefaultStyle();

  //
  // ...
  //
  
  c->Update();
}
//________________________________________________________________________________________
double SuperkExpectation::GeneratorInputs::ND280XSecRatioWght(
    double Etrue, t2kana::yr1numu::conv::ND280Mode_t mnd) const
{
   switch(mnd) {
      case ( kNDmNuMuCCQE     ) : return 1.; break;
      case ( kNDmNuMuCCRES    ) : 
        // There is a mismatch between the event categories in the nd280 sample and the 
        // event categories for which errors are given.
        // I need to propagate errors on CCOther/CCQE and CC1pi/CCQE.
        // The CCRES category contains both CC1pi and CCOther.
        // I estimate that 5% of the CCRES at the T2K energies are not CC1pi,
        // but something else multi-pi, N+gamma, etc.
        // Use that number ...to weight the weight... so that in every true energy bin
        // we apply the CC1pi/CCQE weight to 95% of events and the CCOth/CCQE to the
        // other 5%. Although approximate, this helps assigning the proper CCOth/CCQE
        // error on nd280
	return (0.95 * this->CC1piToCCQEXSecRatioWght(Etrue) +
                0.05 * this->CCOthToCCQEXSecRatioWght(Etrue));
	break;
      case ( kNDmNuMuCCCOH    ) : this->CCOthToCCQEXSecRatioWght(Etrue); break;
      case ( kNDmNuMuCCOth    ) : this->CCOthToCCQEXSecRatioWght(Etrue); break;
      case ( kNDmNuMuNCEL     ) : return 1.; break;
      case ( kNDmNuMuNCRES    ) : this->NCToCCQEXSecRatioWght(Etrue);  break;
      case ( kNDmNuMuNCCOH    ) : this->NCToCCQEXSecRatioWght(Etrue);  break;
      case ( kNDmNuMuNCOth    ) : this->NCToCCQEXSecRatioWght(Etrue);  break;
      case ( kNDmNuMuBarCCQE  ) : return 1.; break;
      case ( kNDmNuMuBarCCRES ) : return 1.; break;
      case ( kNDmNuMuBarCCCOH ) : return 1.; break;
      case ( kNDmNuMuBarCCOth ) : return 1.; break;
      case ( kNDmNuMuBarNCEL  ) : return this->NCToCCQEXSecRatioWght(Etrue); break;
      case ( kNDmNuMuBarNCRES ) : return this->NCToCCQEXSecRatioWght(Etrue); break;
      case ( kNDmNuMuBarNCCOH ) : return this->NCToCCQEXSecRatioWght(Etrue); break;
      case ( kNDmNuMuBarNCOth ) : return this->NCToCCQEXSecRatioWght(Etrue); break;
      case ( kNDmNuECCQE      ) : return 1.; break;
      case ( kNDmNuECCRES     ) : return 1.; break;
      case ( kNDmNuECCCOH     ) : return 1.; break;
      case ( kNDmNuECCOth     ) : return 1.; break;
      case ( kNDmNuENCEL      ) : return this->NCToCCQEXSecRatioWght(Etrue); break;
      case ( kNDmNuENCRES     ) : return this->NCToCCQEXSecRatioWght(Etrue); break;
      case ( kNDmNuENCCOH     ) : return this->NCToCCQEXSecRatioWght(Etrue); break;
      case ( kNDmNuENCOth     ) : return this->NCToCCQEXSecRatioWght(Etrue); break;
      default                   : return 1.; break;
   }
   return 1.;

}
//________________________________________________________________________________________
double SuperkExpectation::GeneratorInputs::SuperkXSecRatioWght(
   double Etrue, t2kana::yr1numu::conv::SuperkMode_t m) const
{
   switch(m) {
      case ( kSKmNuMuCCQE  ) : return 1.; break;
      case ( kSKmNuMuCC1pi ) : return this->CC1piToCCQEXSecRatioWght(Etrue); break;
      case ( kSKmNuMuCCoth ) : return this->CCOthToCCQEXSecRatioWght(Etrue); break;
      case ( kSKmNC        ) : return this->NCToCCQEXSecRatioWght(Etrue);    break;
      case ( kSKmNuMuBar   ) : return 1.; break;
      case ( kSKmNuE       ) : return 1.; break;
      default                : return 1.; break;
   }
   return 1.;
}
//________________________________________________________________________________________
double SuperkExpectation::GeneratorInputs::SuperkCCQEXSecWght(
    double Etrue, t2kana::yr1numu::conv::SuperkMode_t m) const
{
   switch(m) {
      case ( kSKmNuMuCCQE  ) : 
      {
          assert(fCCQEXSecRelErr);
          Etrue = TMath::Min( Etrue, fCCQEXSecRelErr->GetXaxis()->GetXmax() );
          int    bin    = fCCQEXSecRelErr->FindBin(Etrue);
          double relerr = fCCQEXSecRelErr->GetBinContent(bin);
          double wght   = (1. + fCCQEXSecTwkDial*relerr);
          return wght;

          break;
      }
      case ( kSKmNuMuCC1pi ) : return 1.; break;
      case ( kSKmNuMuCCoth ) : return 1.; break;
      case ( kSKmNC        ) : return 1.; break;
      case ( kSKmNuMuBar   ) : return 1.; break;
      case ( kSKmNuE       ) : return 1.; break;
      default                : return 1.; break;
   }
   return 1.;
}
//________________________________________________________________________________________
double SuperkExpectation::GeneratorInputs::SuperkNuclEffectWght(
   double Ereco, t2kana::yr1numu::conv::SuperkMode_t /*m*/) const
{
   assert(fSuperkNuclearEffectRelErr); //<-- note: units in MeV, all other NIWG histos have units in GeV
   double mev2gev = 1E-3;
   double E       = TMath::Min( Ereco, fSuperkNuclearEffectRelErr->GetXaxis()->GetXmax() * mev2gev);
   int    bin     = fSuperkNuclearEffectRelErr->FindBin(E/mev2gev);
   double relerr  = fSuperkNuclearEffectRelErr->GetBinContent(bin);
   double twk     = fNuclEffectTwkDial;
   double wght    = (1. + twk*relerr);
   return wght;
}
//________________________________________________________________________________________
double SuperkExpectation::GeneratorInputs::CC1piToCCQEXSecRatioWght(double Etrue) const
{
   return this->XSecRatioWght(
     Etrue, fCC1piToCCQEXSecRatioTwkDial, fCC1piToCCQEXSecRatioRelErr);
}
//________________________________________________________________________________________
double SuperkExpectation::GeneratorInputs::CCOthToCCQEXSecRatioWght(double Etrue) const
{
   return this->XSecRatioWght(
     Etrue, fCCOthToCCQEXSecRatioTwkDial, fCCOthToCCQEXSecRatioRelErr);
}
//________________________________________________________________________________________
double SuperkExpectation::GeneratorInputs::NCToCCQEXSecRatioWght(double Etrue) const
{
   return this->XSecRatioWght(
     Etrue, fNCToCCQEXSecRatioTwkDial, fNCToCCQEXSecRatioRelErr);
}
//________________________________________________________________________________________
double SuperkExpectation::GeneratorInputs::XSecRatioWght(
   double E, double twk, TH1D * hRelErr) const
{
   assert(hRelErr);

   E = TMath::Min( E, hRelErr->GetXaxis()->GetXmax() );

   int    bin   = hRelErr->FindBin(E);
   double relerr = hRelErr->GetBinContent(bin);

   double wght = (1. + twk*relerr);

#ifdef __debug_level_3__
   cout << "Exp/NuGen/XSecRatioWght >> E = " << E 
        << ", relerr = " << relerr << ", twkdial = " << twk << ", wght = " << wght 
        << endl;
#endif

   return wght; 
}
//________________________________________________________________________________________
bool SuperkExpectation::GeneratorInputs::TweakAffectsND280 (void) const
{
  if ( TMath::Abs(fCC1piToCCQEXSecRatioTwkDial) > kASmallNum ) return true;
  if ( TMath::Abs(fCCOthToCCQEXSecRatioTwkDial) > kASmallNum ) return true;
  if ( TMath::Abs(fNCToCCQEXSecRatioTwkDial   ) > kASmallNum ) return true;
  return false;
}
//________________________________________________________________________________________
bool SuperkExpectation::GeneratorInputs::TweakAffectsSuperk(void) const
{
  if ( TMath::Abs(fCCQEXSecTwkDial            ) > kASmallNum ) return true;
  if ( TMath::Abs(fCC1piToCCQEXSecRatioTwkDial) > kASmallNum ) return true;
  if ( TMath::Abs(fCCOthToCCQEXSecRatioTwkDial) > kASmallNum ) return true;
  if ( TMath::Abs(fNCToCCQEXSecRatioTwkDial   ) > kASmallNum ) return true;
  if ( TMath::Abs(fNuclEffectTwkDial          ) > kASmallNum ) return true;
  return false;
}
//________________________________________________________________________________________

