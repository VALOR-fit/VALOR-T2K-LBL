//________________________________________________________________________________________
/*

*/
//________________________________________________________________________________________

#include <cassert>
#include <cstdlib>
#include <iostream>

#include <TFile.h>
#include <TString.h>

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
double SuperkExpectation::Generate(double Ereco, const PhysParamList & fp)
{
// Calculate dN^{exp}_{SuperK}(E_{reco}) / dE_{reco} for the input E_{reco} and 
// the input physics parameter list.
//

  if(!fIsConfigured) {
    cerr << "Exp >> Not configured properly! " << endl;
    return -1;
  }

  double expectation = 0.;

  // input physics params (params included in the fit / systematics params)
  double dmsq23          = fp.CurrValue(kPhypDm23sq);
  double sinsq2theta23   = fp.CurrValue(kPhypSinsq2Theta23);

  // calc stuff with no energy dependency
  double Mfiducial_nd280 = fND280Inputs  -> FiducialMass();
  double Npot_nd280      = fND280Inputs  -> AnalyzedSamplePOT();
  double A_nd280         = fND280Inputs  -> EffectiveTargetMassNumber();
  double Mfiducial_sk    = fSuperkInputs -> FiducialMass();
  double Npot_sk         = fSuperkInputs -> AnalyzedSamplePOT();
  double A_sk            = fSuperkInputs -> EffectiveTargetMassNumber();
  double C0 = (Npot_sk    * Mfiducial_sk    * A_nd280) / 
              (Npot_nd280 * Mfiducial_nd280 * A_sk);

  cout << "C0 = " << C0 << endl;

  const unsigned int ns = fSamples.size();
  const unsigned int nm = 3;  
  SuperkMode_t sk_modes [nm] = { kSKmNuMuCCQE, kSKmNuMuNonCCQE, kSKmNC };

  const TAxis & etrue_binning = fSuperkInputs->EtrueBinning();

  // loop over SK samples & modes/sample
  for(int is = 0; is < ns; is++) {
    SuperkSample_t s = fSamples[is];
    for(int im = 0; im < nm; im++) {
       SuperkMode_t m = sk_modes[im];

       double IntgEtrue = 0; // integral over E_{true}

       // E_{true} bin loop
       for(int itbin = 1; itbin <= etrue_binning.GetNbins(); itbin++) {
         double Etrue = etrue_binning.GetBinCenter(itbin);

         // calc stuff with E_{true} energy dependency only
         double fn_ratio          = fBeamInputs->FNRatio(Etrue);
         double xsec_totcc_nd280  = fGeneratorInputs -> SigmaTotCCND280(Etrue);
         double xsec_tot_sk       = fGeneratorInputs -> SigmaTotSK(Etrue);
         double evtfrac_m_sk      = fGeneratorInputs -> FractModeSK(Etrue,m);
         double xsec_m_sk         = xsec_tot_sk * evtfrac_m_sk;
         double prob_surv         = util::NuMuSurvProb2D(sinsq2theta23, dmsq23, Etrue, m);
         double R_s_m_sk          = fSuperkInputs -> EnergyResolutionMatrix(Etrue,Ereco,s,m);
         double eff_sel_s_m_sk    = fSuperkInputs -> SelectionEfficiency(Ereco,s,m);
         double eff_reco_m_sk     = fSuperkInputs -> ReconstructionEfficiency(Etrue,m);
         double nobscclike_nd280  = fND280Inputs -> NumOfObsCCLikeEvents ();
         double nexpnc_nd280      = fND280Inputs -> NumOfEstNCBkgEvents  ();
         double eff_nd280         = fND280Inputs  -> Efficiency(Etrue); 
/*
         cout << "** " << util::AsString(m) << ", Etrue = " << Etrue 
              << ", F/N = " << fn_ratio 
              << ", Psurv = " << prob_surv 
              << ", SigTot(SK) = " << xsec_tot_sk << ", ModeFrac(SK) = " << evtfrac_m_sk
              << ", SigCCTot(nd280) = " << xsec_totcc_nd280 
              << ", RSK = " << R_s_m_sk 
              << ", EffSelSK = " << eff_sel_s_m_sk << ", EffRecoSK = " << eff_reco_m_sk << endl;
*/         
         double Ct = fn_ratio * 
                     (xsec_m_sk / xsec_totcc_nd280) * 
                     prob_surv * 
                     (nobscclike_nd280 - nexpnc_nd280) *
                     R_s_m_sk * (eff_sel_s_m_sk * eff_reco_m_sk / eff_nd280);     

         IntgEtrue += Ct;
       }//E_{true}

       expectation += (C0*IntgEtrue);

       cout << util::AsString(s) << "/" << util::AsString(m) 
            << " ==> adding to expectation = " << (C0*IntgEtrue) << endl;
    }//m
  }//s

  return expectation;
}
//________________________________________________________________________________________
bool SuperkExpectation::Read(string config_file)
{
  cout << "Exp >> " << "Reading config file: " << config_file << endl;

  fND280Inputs     = new ND280Inputs;
  fSuperkInputs    = new SuperkInputs;
  fBeamInputs      = new BeamInputs; 
  fGeneratorInputs = new GeneratorInputs; 

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
}
//________________________________________________________________________________________
void SuperkExpectation::CleanUp(void)
{
  if (fND280Inputs    ) delete fND280Inputs;
  if (fSuperkInputs   ) delete fSuperkInputs;
  if (fBeamInputs     ) delete fBeamInputs; 
  if (fGeneratorInputs) delete fGeneratorInputs;

  fSamples.clear();
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
  fAeff            = 0;  
  fFidMass         = 0;
  fPOT             = 0;
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
bool SuperkExpectation::ND280Inputs::Read(TXMLEngine* xml, XMLNodePointer_t root_node)
{
  cout << "Exp/nd280 >> Reading nd280 inputs from input XML file." << endl;

  // Get observed number of CC-like events
  fNumObsCCLike = util::XML2Dbl(xml,root_node,"nd280/ncclike");
  // Get nominal estimated number of NC bkg events
  fNumEstNCBkg  = util::XML2Dbl(xml,root_node,"nd280/nncbkg"); 
  // Get fiducial mass
  fFidMass = util::XML2Dbl(xml,root_node,"nd280/fiducial_mass");
  // Get POT
  fPOT = util::XML2Dbl(xml,root_node,"nd280/pot");
  // Get mass number
  fAeff = util::XML2Dbl(xml,root_node,"nd280/mass_number");

  cout << "Exp/nd280 >> " << "N^{obs}_{CClike} : " << fNumObsCCLike << endl;
  cout << "Exp/nd280 >> " << "N^{exp}_{NCbkg}  : " << fNumEstNCBkg  << endl;
  cout << "Exp/nd280 >> " << "M_{fiducial}     : " << fFidMass << " kton" << endl;
  cout << "Exp/nd280 >> " << "N_{POT}          : " << fPOT     << endl;
  cout << "Exp/nd280 >> " << "A_{effective}    : " << fAeff    << endl;

  cout << "Exp/nd280 >> Done reading nd280 inputs!" << endl;

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
double SuperkExpectation::ND280Inputs::EffectiveTargetMassNumber(void) const
{
  return fAeff;
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
  fAeff      = 0;  
  fFidMass   = 0; 
  fPOT       = 0;  
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
bool SuperkExpectation::SuperkInputs::Read(TXMLEngine* xml, XMLNodePointer_t root_node)
{
// Read SuperK inputs to disappearance analysis
//
  cout << "Exp/SK >> Reading SuperK inputs from input XML file." << endl;

  // Get reconstruction efficiencies
  fRecoEff[kSKmNuMuCCQE] = 
    util::XML2TH1D(xml,root_node,"superk/reconstruction_efficiency/ccqe");
  fRecoEff[kSKmNuMuNonCCQE] = 
    util::XML2TH1D(xml,root_node,"superk/reconstruction_efficiency/nonccqe");
  fRecoEff[kSKmNC] = 
    util::XML2TH1D(xml,root_node,"superk/reconstruction_efficiency/nc");

  // Get selection efficiencies
  fSelEff[kSKsFCFV1RingMuLike][kSKmNuMuCCQE] = 
    util::XML2TH1D(xml,root_node,"superk/selection_efficiency/one_mu_ring_tight/ccqe");
  fSelEff[kSKsFCFV1RingMuLike][kSKmNuMuNonCCQE] = 
    util::XML2TH1D(xml,root_node,"superk/selection_efficiency/one_mu_ring_tight/nonccqe");
  fSelEff[kSKsFCFV1RingMuLike][kSKmNC] = 
    util::XML2TH1D(xml,root_node,"superk/selection_efficiency/one_mu_ring_tight/nc");

  // Get energy resolution matrices
  fResolMtrx[kSKsFCFV1RingMuLike][kSKmNuMuCCQE] =
    util::XML2TH2D(xml,root_node,"superk/energy_resolution_matrix/one_mu_ring_tight/ccqe");
  fResolMtrx[kSKsFCFV1RingMuLike][kSKmNuMuNonCCQE] = 
    util::XML2TH2D(xml,root_node,"superk/energy_resolution_matrix/one_mu_ring_tight/nonccqe");
  fResolMtrx[kSKsFCFV1RingMuLike][kSKmNC] = 
    util::XML2TH2D(xml,root_node,"superk/energy_resolution_matrix/one_mu_ring_tight/nc");

//HACK: THEY SHOULD HAVE BEEN NORMALIZED
  fResolMtrx[kSKsFCFV1RingMuLike][kSKmNuMuCCQE] -> Scale(
          1./fResolMtrx[kSKsFCFV1RingMuLike][kSKmNuMuCCQE]->Integral());
  fResolMtrx[kSKsFCFV1RingMuLike][kSKmNuMuNonCCQE] -> Scale(
          1./fResolMtrx[kSKsFCFV1RingMuLike][kSKmNuMuNonCCQE]->Integral());
  fResolMtrx[kSKsFCFV1RingMuLike][kSKmNC] -> Scale(
          1./fResolMtrx[kSKsFCFV1RingMuLike][kSKmNC]->Integral());


  // Get fiducial mass
  fFidMass = util::XML2Dbl(xml,root_node,"superk/fiducial_mass");
  // Get POT
  fPOT = util::XML2Dbl(xml,root_node,"superk/pot");
  // Get mass number
  fAeff = util::XML2Dbl(xml,root_node,"superk/mass_number");
  cout << "Exp/SK >> " << "M_{fiducial}  : " << fFidMass << " kton" << endl;
  cout << "Exp/SK >> " << "N_{POT}       : " << fPOT     << endl;
  cout << "Exp/SK >> " << "A_{effective} : " << fAeff    << endl;

  // Get Ereco and Etrue binning
  fErBinning = util::XML2BinInfo(xml,root_node,"superk/binning/ereco");
  fEtBinning = util::XML2BinInfo(xml,root_node,"superk/binning/etrue");

  cout << "Exp/SK >> Done reading SuperK inputs!" << endl;

  return true;
}
//________________________________________________________________________________________
void SuperkExpectation::SuperkInputs::Draw(TVirtualPad * c, Option_t * option) const
{
  if(!c) return;

  c->cd();

  TString opt = option;
  util::SetDefaultStyle();

  if(opt.Contains("ReconstructionEfficiency")) { 
    fRecoEff[kSKmNuMuCCQE]    -> Draw();
    fRecoEff[kSKmNuMuNonCCQE] -> Draw("same");
    fRecoEff[kSKmNC]          -> Draw("same");
  } 
  else
  if(opt.Contains("SelectionEfficiency")) { 
    if(opt.Contains("OneMuRingTight")) { 
       fSelEff[kSKsFCFV1RingMuLike][kSKmNuMuCCQE]    -> Draw();
       fSelEff[kSKsFCFV1RingMuLike][kSKmNuMuNonCCQE] -> Draw("same");
       fSelEff[kSKsFCFV1RingMuLike][kSKmNC]          -> Draw("same");
    } 
  } 
  else
  if(opt.Contains("EnergyResolutionMatrix")) { 
    if(opt.Contains("OneMuRingTight")) { 
       c->Divide(1,3);
       c->cd(1);
       fResolMtrx[kSKsFCFV1RingMuLike][kSKmNuMuCCQE]-> Draw();
       c->cd(2);
       fResolMtrx[kSKsFCFV1RingMuLike][kSKmNuMuNonCCQE]-> Draw();
       c->cd(3);
       fResolMtrx[kSKsFCFV1RingMuLike][kSKmNC]->Draw();
    }
  }

  c->Update();
}
//________________________________________________________________________________________
double SuperkExpectation::SuperkInputs::EffectiveTargetMassNumber(void) const
{
  return fAeff;
}
//________________________________________________________________________________________
double SuperkExpectation::SuperkInputs::FiducialMass(void) const
{
  return fFidMass; 
}
//________________________________________________________________________________________
double SuperkExpectation::SuperkInputs::AnalyzedSamplePOT(void) const
{
  return fPOT;  
}
//________________________________________________________________________________________
double SuperkExpectation::SuperkInputs::ReconstructionEfficiency(
    double Etrue, SuperkMode_t m) const
{
  TH1D * h = fRecoEff[m];
  assert(h);
  int bin = h->FindBin(Etrue);
  double eff = h->GetBinContent(bin);
  return eff;
}
//________________________________________________________________________________________
double SuperkExpectation::SuperkInputs::SelectionEfficiency(
    double Ereco, SuperkSample_t s, SuperkMode_t m) const
{
  TH1D * h = fSelEff[s][m];
  assert(h);
  int bin = h->FindBin(Ereco);
  double eff = h->GetBinContent(bin);
  return eff;
}
//________________________________________________________________________________________
double SuperkExpectation::SuperkInputs::EnergyResolutionMatrix(
    double Etrue, double Ereco, SuperkSample_t s, SuperkMode_t m) const
{
  TH2D * h2 = fResolMtrx[s][m];
  assert(h2);
  int binx = h2->GetXaxis()->FindBin(Ereco);
  int biny = h2->GetXaxis()->FindBin(Etrue);
  double Rij = h2->GetBinContent(binx,biny);
  return Rij;
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
bool SuperkExpectation::BeamInputs::Read(TXMLEngine* xml, XMLNodePointer_t root_node)
{
// Read beam inputs to disappearance analysis
//
  cout << "Exp/Beam >> Reading beam inputs from input XML file." << endl;

  // Get SuperK flux
  fSuperkFlux = util::XML2TH1D(xml,root_node,"beam/superk_flux");
  if(!fSuperkFlux) {
    cerr << "Exp/Beam >> No input SuperK flux histogram!" << endl;
    exit(1);
  }

  // Get ND280 flux
  fND280Flux  = util::XML2TH1D(xml,root_node,"beam/nd280_flux");
  if(!fND280Flux) {
    cerr << "Exp/Beam >> No input nd280 flux histogram!" << endl;
    exit(1);
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

  if(opt.Contains("ND280Flux")) { 
    fND280Flux->Draw();
    fND280Flux->GetXaxis()->SetTitle("E_{#nu} (GeV)");
    fND280Flux->GetYaxis()->SetTitle("ND280 flux");
  } 
  else
  if(opt.Contains("SKFlux")) { 
    fSuperkFlux->Draw();
    fSuperkFlux->GetXaxis()->SetTitle("E_{#nu} (GeV)");
    fSuperkFlux->GetYaxis()->SetTitle("SuperK flux");
  } 
  else
  if(opt.Contains("OffAxisRatio")) { 
    TH1D * h = new TH1D("h","",100,0,5);
    for(int i=1; i<=h->GetNbinsX(); i++) {
	double e = h->GetBinCenter(i);
	double fn = this->FNRatio(e);
	h->SetBinContent(i,fn);
    }
    h->Draw();
    h->GetXaxis()->SetTitle("E_{#nu} (GeV)");
    h->GetYaxis()->SetTitle("SuperK/ND280 flux ratio");
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
    fFractModeSK[im] = 0;
  }
  fSigmaTotSK      = 0;
  fSigmaTotCCND280 = 0;
}
//________________________________________________________________________________________
void SuperkExpectation::GeneratorInputs::CleanUp(void)
{
  for(int im=0; im<knm; im++) {
    if ( fFractModeSK[im] ) delete fFractModeSK[im];
  }
  if ( fSigmaTotSK      ) delete fSigmaTotSK;
  if ( fSigmaTotCCND280 ) delete fSigmaTotCCND280;
}
//________________________________________________________________________________________
bool SuperkExpectation::GeneratorInputs::Read(TXMLEngine* xml, XMLNodePointer_t root_node)
{
  cout << "Exp/NuGen >> Reading neutrino generator inputs from input XML file." << endl;

  // Get numu+H20 total cross section = f(E_true)
  fSigmaTotSK = util::XML2TH1D(xml,root_node,"genie/sigma_total_H2O");
  if(!fSigmaTotSK) {
    cerr << "Exp/NuGen >> No input numu+H20 total cross section histogram!" << endl;
    exit(1);
  }
  // Get event fractions for different modes in numu+H20 interactions  = f(E_true)
  fFractModeSK[kSKmNuMuCCQE]    = util::XML2TH1D(xml,root_node,"genie/fraction_mode_H2O/ccqe");
  fFractModeSK[kSKmNuMuNonCCQE] = util::XML2TH1D(xml,root_node,"genie/fraction_mode_H2O/nonccqe");
  fFractModeSK[kSKmNC]          = util::XML2TH1D(xml,root_node,"genie/fraction_mode_H2O/nc");
  if ( !fFractModeSK[kSKmNuMuCCQE]    ||
       !fFractModeSK[kSKmNuMuNonCCQE] ||
       !fFractModeSK[kSKmNC] ) {
    cerr << "Exp/NuGen >> No input numu+H20 event fractions for CC,nCCQE,NC modes!" << endl;
    exit(1);
  }

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

  if(opt.Contains("SKSigTotal")) { 
    fSigmaTotSK->Draw();
    fSigmaTotSK->GetXaxis()->SetTitle("E_{#nu} (GeV)");
    fSigmaTotSK->GetYaxis()->SetTitle("#sigma_{tot} (#nu_{#mu}+H20), 1E-38 cm^{2}");
  }
  else
  if(opt.Contains("SKSigMode")) { 
  }
  else
  if(opt.Contains("SKSigAll")) { 
  }
  else
  if(opt.Contains("SKFracMode")) { 
    fFractModeSK[kSKmNuMuCCQE]    -> Draw();
    fFractModeSK[kSKmNuMuNonCCQE] -> Draw("same");
    fFractModeSK[kSKmNC]          -> Draw("same");
    fFractModeSK[kSKmNuMuCCQE] ->GetXaxis()->SetTitle("E_{#nu} (GeV)");
    fFractModeSK[kSKmNuMuCCQE] ->GetYaxis()->SetTitle("event fraction");
  }
  
  c->Update();
}
//________________________________________________________________________________________
double SuperkExpectation::GeneratorInputs::FractModeSK(double Etrue, SuperkMode_t m) const
{
  TH1D * h = fFractModeSK[m];
  assert(h);

  int    bin   = h->FindBin(Etrue);
  double frac  = h->GetBinContent(bin);
  return frac;
}
//________________________________________________________________________________________
double SuperkExpectation::GeneratorInputs::SigmaTotSK(double Etrue) const
{
  TH1D* h = fSigmaTotSK;
  assert(h);

  int    bin  = h->FindBin(Etrue);
  double sig  = h->GetBinContent(bin);
  return sig;
}
//________________________________________________________________________________________
double SuperkExpectation::GeneratorInputs::SigmaModeSK(double Etrue, SuperkMode_t m) const
{
  double sigtot  = this->SigmaTotSK(Etrue);
  double frcmode = this->FractModeSK(Etrue,m);
  double sigmode = sigtot * frcmode; 
  return sigmode;
}
//________________________________________________________________________________________
double SuperkExpectation::GeneratorInputs::SigmaTotCCND280(double Etrue) const
{
//TODO
//CURRENTLY ASSUME A NEAR SK TOO
  return this->SigmaTotSK(Etrue);
}
//________________________________________________________________________________________


