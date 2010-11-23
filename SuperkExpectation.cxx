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

  const unsigned int ns = fSamples.size();
  const unsigned int nm = 3;  
  SuperkMode_t sk_modes [nm] = { kSKmNuMuCCQE, kSKmNuMuNonCCQE, kSKmNC };

  // loop over SK samples & modes/sample
  for(unsigned int is = 0; is < ns; is++) {
    SuperkSample_t s = fSamples[is];
    for(unsigned int im = 0; im < nm; im++) {
       SuperkMode_t m = sk_modes[im];
       double mode_contrib_in_sample = this->Generate(Ereco,fp,s,m);
       expectation += mode_contrib_in_sample;
    }//m
  }//s

#ifdef __debug_level_1__
  cout << "Exp / Generate() >> "
       << "NEvents(all_analysis_samples/all_modes, Ereco = " << Ereco << " GeV) = "
       << expectation << endl;
#endif

  return expectation;
}
//________________________________________________________________________________________
double SuperkExpectation::Generate(
  double Ereco, const PhysParamList & fp, SuperkSample_t s, SuperkMode_t m)
{
  if(!fIsConfigured) {
    cerr << "Exp >> Not configured properly! " << endl;
    return -1;
  }
  
  const double kton2gr        = 1E+9;
  const double xsec_scale     = 1E+38;

  double expectation = 0;

  // input physics params (params included in the fit / systematics params)
  double dmsq23          = fp.CurrValue(kPhypDm23sq);
  double sinsq2theta23   = fp.CurrValue(kPhypSinsq2Theta23);
  double norm            = (fp.Exists(kPhypNorm)) ? fp.CurrValue(kPhypNorm) : 1.;

  // calc stuff with no energy dependency
  double Mfiducial_sk      = fSuperkInputs -> FiducialMass() * kton2gr;
  double Npot_sk           = fSuperkInputs -> AnalyzedSamplePOT();
  double A_sk              = fSuperkInputs -> EffectiveTargetMassNumber();
  double nobscclike_nd280  = fND280Inputs  -> NObsCCLike();
  double nexpncbkg_nd280   = fND280Inputs  -> NExpNCBkg();
  double nexpcc_nd280      = fND280Inputs  -> NExpCC();

  double C0 = Npot_sk * (Mfiducial_sk / A_sk) * kNA;
  double flux_correction = (nobscclike_nd280 - nexpncbkg_nd280) / nexpcc_nd280;

  const TAxis & etrue_binning = fSuperkInputs->EtrueBinning();

  // E_{true} bin loop
  for(int itbin = 1; itbin <= etrue_binning.GetNbins(); itbin++) {
    double Etrue = etrue_binning.GetBinCenter(itbin);

    // calc stuff with E_{true} energy dependency only
    double flux_sk_mc      = fBeamInputs->SuperkFlux(Etrue); // per cm2 per POT
    double flux_sk         = flux_sk_mc * flux_correction;
    double xsec_m_sk       = fGeneratorInputs->SigmaModeSK(Etrue,m) / xsec_scale;
    double prob_surv       = util::NuMuSurvProb2D(sinsq2theta23, dmsq23, Etrue, m);
    double R_s_m_sk        = fSuperkInputs -> EnergyResolutionMatrix(Etrue,Ereco,s,m);
    double eff_sel_s_m_sk  = fSuperkInputs -> SelectionEfficiency(Ereco,s,m);
    double eff_reco_m_sk   = fSuperkInputs -> ReconstructionEfficiency(Etrue,m);

    double expectation_from_true_bin = 
        C0 * flux_sk * xsec_m_sk * prob_surv * R_s_m_sk * eff_sel_s_m_sk * eff_reco_m_sk;     

#ifdef __debug_level_3__
    cout << "** " << util::AsString(m) << ", Etrue = " << Etrue << " GeV"
         << "norm = " << norm << ", Psurv = " << prob_surv 
         << ", F/N = " << fn_ratio << "F(nd280) = " << flux_nd280 << ", F(SK) = " << flux_sk
         << ", SigModeFrac(SK) = " << xsec_m_sk
         << ", RSK = " << R_s_m_sk 
         << ", EffSelSK = " << eff_sel_s_m_sk << ", EffRecoSK = " << eff_reco_m_sk 
         << " ==> expectation(from Etrue bin) = " << expectation_from_true_bin 
         << endl;
#endif

    expectation += expectation_from_true_bin;

  }//E_{true}

  if(fTwkBeamNorm != 0.) {
    double sigma = 0.1;
    expectation *= (1 + fTwkBeamNorm * sigma);
  }

#ifdef __debug_level_2__
  cout << "Exp / Generate() >> "
       << "NEvents(" << util::AsString(s) << "/" << util::AsString(m) 
       << ", Ereco = " << Ereco << " GeV) = " << expectation << endl;
#endif

  return expectation * norm;
}
//________________________________________________________________________________________
TH1D * SuperkExpectation::AsTH1D(const PhysParamList & fp)
{
  int    nbins = fSuperkInputs->fErBinning->GetNbins();
  double emin  = fSuperkInputs->fErBinning->GetXmin();
  double emax  = fSuperkInputs->fErBinning->GetXmax();

  TH1D * hEreco = new TH1D("hEreco", "", nbins, emin, emax);
  hEreco->SetDirectory(0);

  for(int bin=1; bin < hEreco->GetNbinsX(); bin++) {
    double Ereco = hEreco->GetBinCenter(bin);
    double nev   = this->Generate(Ereco,fp);
    hEreco->SetBinContent(bin,nev);
  }

  return hEreco;
}
//________________________________________________________________________________________
TH1D * SuperkExpectation::AsTH1D(const PhysParamList & fp, SuperkSample_t s)
{
  int    nbins = fSuperkInputs->fErBinning->GetNbins();
  double emin  = fSuperkInputs->fErBinning->GetXmin();
  double emax  = fSuperkInputs->fErBinning->GetXmax();

  const unsigned int nm = 3;  
  SuperkMode_t sk_modes [nm] = { kSKmNuMuCCQE, kSKmNuMuNonCCQE, kSKmNC };

  TH1D * hEreco = new TH1D("hEreco", "", nbins, emin, emax);
  hEreco->SetDirectory(0);

  for(int bin=1; bin < hEreco->GetNbinsX(); bin++) {
    double Ereco = hEreco->GetBinCenter(bin);
    double nev = 0;
    for(unsigned int im = 0; im < nm; im++) {
       SuperkMode_t m = sk_modes[im];
       nev  += (this->Generate(Ereco,fp,s,m));
    }//m
    hEreco->SetBinContent(bin,nev);
  }//bin
  return hEreco;
}
//________________________________________________________________________________________
TH1D * SuperkExpectation::AsTH1D(const PhysParamList & fp, SuperkMode_t m)
{
  int    nbins = fSuperkInputs->fErBinning->GetNbins();
  double emin  = fSuperkInputs->fErBinning->GetXmin();
  double emax  = fSuperkInputs->fErBinning->GetXmax();

  TH1D * hEreco = new TH1D("hEreco", "", nbins, emin, emax);
  hEreco->SetDirectory(0);

  const unsigned int ns = fSamples.size();

  for(int bin=1; bin < hEreco->GetNbinsX(); bin++) {
    double Ereco = hEreco->GetBinCenter(bin);
    double nev = 0;
    // loop over SK samples 
    for(unsigned int is = 0; is < ns; is++) {
       SuperkSample_t s = fSamples[is];
       nev  += (this->Generate(Ereco,fp,s,m));
    }//s
    hEreco->SetBinContent(bin,nev);
  }//bin
  return hEreco;
}
//________________________________________________________________________________________
TH1D * SuperkExpectation::AsTH1D(
  const PhysParamList & fp, SuperkSample_t s, SuperkMode_t m)
{
  int    nbins = fSuperkInputs->fErBinning->GetNbins();
  double emin  = fSuperkInputs->fErBinning->GetXmin();
  double emax  = fSuperkInputs->fErBinning->GetXmax();

  TH1D * hEreco = new TH1D("hEreco", "", nbins, emin, emax);
  hEreco->SetDirectory(0);

  for(int bin=1; bin < hEreco->GetNbinsX(); bin++) {
    double Ereco = hEreco->GetBinCenter(bin);
    double nev   = this->Generate(Ereco,fp,s,m);
    hEreco->SetBinContent(bin,nev);
  }

  return hEreco;
}
//________________________________________________________________________________________
bool SuperkExpectation::Read(string config_file)
{
  cout << "Exp >> " << "Reading config file: " << config_file << endl;

  if(fIsConfigured) {
    this->CleanUp();
  } 

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

  fTwkBeamNorm     = 0;
  fTwkSuperkNorm   = 0;
  fTwkSuperkEScale = 0;

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
  fEtBinning      = 0;
  fPOTData        = 0;
  fNumObsCCLike   = 0;
  fPOTMC          = 0;
  for(int i = 0; i < knmnd; i++) {
    fNumSelEt[i] = 0;
    fSelEffEt[i] = 0;
  }
}
//________________________________________________________________________________________
void SuperkExpectation::ND280Inputs::CleanUp(void)
{
  if (fEtBinning) delete fEtBinning;

  for(int i = 0; i < knmnd; i++) {
    if (fNumSelEt[i]) delete fNumSelEt[i];
    if (fSelEffEt[i]) delete fSelEffEt[i];
  }
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

  fNumSelEt [kNDmNuMuCCQE]     -> Sumw2();
  fNumSelEt [kNDmNuMuCCRES]    -> Sumw2();
  fNumSelEt [kNDmNuMuCCCOH]    -> Sumw2();
  fNumSelEt [kNDmNuMuCCOth]    -> Sumw2();
  fNumSelEt [kNDmNuMuNCEL]     -> Sumw2();
  fNumSelEt [kNDmNuMuNCRES]    -> Sumw2();
  fNumSelEt [kNDmNuMuNCCOH]    -> Sumw2();
  fNumSelEt [kNDmNuMuNCOth]    -> Sumw2();
  fNumSelEt [kNDmNuMuBarCCQE]  -> Sumw2();
  fNumSelEt [kNDmNuMuBarCCRES] -> Sumw2();
  fNumSelEt [kNDmNuMuBarCCCOH] -> Sumw2();
  fNumSelEt [kNDmNuMuBarCCOth] -> Sumw2();
  fNumSelEt [kNDmNuMuBarNCEL]  -> Sumw2();
  fNumSelEt [kNDmNuMuBarNCRES] -> Sumw2();
  fNumSelEt [kNDmNuMuBarNCCOH] -> Sumw2();
  fNumSelEt [kNDmNuMuBarNCOth] -> Sumw2();
  fNumSelEt [kNDmNuECCQE]      -> Sumw2();
  fNumSelEt [kNDmNuECCRES]     -> Sumw2();
  fNumSelEt [kNDmNuECCCOH]     -> Sumw2();
  fNumSelEt [kNDmNuECCOth]     -> Sumw2();
  fNumSelEt [kNDmNuENCEL]      -> Sumw2();
  fNumSelEt [kNDmNuENCRES]     -> Sumw2();
  fNumSelEt [kNDmNuENCCOH]     -> Sumw2();
  fNumSelEt [kNDmNuENCOth]     -> Sumw2();

  cout << "Exp/nd280 >> Done rescaling NMC_{mode} histograms to the data POT" << endl;

  cout << "Exp/nd280 >> Done reading nd280 inputs!" << endl;

  return true;
}
//________________________________________________________________________________________
void SuperkExpectation::ND280Inputs::Draw(TVirtualPad * c, Option_t * option) const
{
  if(!c) return;

  c->cd();

  util::SetDefaultStyle();

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

  TString opt = option;

  if(opt.Contains("Saclay/NuMu/CC")) { 
    fNumSelEt [kNDmNuMuCCQE]     -> Draw();
    fNumSelEt [kNDmNuMuCCRES]    -> Draw("SAME");
    fNumSelEt [kNDmNuMuCCCOH]    -> Draw("SAME");
    fNumSelEt [kNDmNuMuCCOth]    -> Draw("SAME");
    fNumSelEt [kNDmNuMuCCQE]     -> GetXaxis()->SetTitle("E_{#nu} (GeV)");
    fNumSelEt [kNDmNuMuCCQE]     -> GetYaxis()->SetTitle(Form("events / %d POT", fPOTData));
  }
  else
  if(opt.Contains("SKSigAll")) { 
  }
 
  c->Update();
}
//________________________________________________________________________________________
double SuperkExpectation::ND280Inputs::NObsCCLike(void) const 
{
// observed number of CC-like events

  return fNumObsCCLike;
}
//________________________________________________________________________________________
double SuperkExpectation::ND280Inputs::NExpNCBkg(void) const
{
// expected NC contamination of CC-like sample

  double nexpncbkg =
       fNumSelEt [kNDmNuMuNCEL]     -> Integral() +
       fNumSelEt [kNDmNuMuNCRES]    -> Integral() +
       fNumSelEt [kNDmNuMuNCCOH]    -> Integral() +
       fNumSelEt [kNDmNuMuNCOth]    -> Integral() +
       fNumSelEt [kNDmNuMuBarNCEL]  -> Integral() +
       fNumSelEt [kNDmNuMuBarNCRES] -> Integral() +
       fNumSelEt [kNDmNuMuBarNCCOH] -> Integral() +
       fNumSelEt [kNDmNuMuBarNCOth] -> Integral() +
       fNumSelEt [kNDmNuENCEL]      -> Integral() +
       fNumSelEt [kNDmNuENCRES]     -> Integral() +
       fNumSelEt [kNDmNuENCCOH]     -> Integral() +
       fNumSelEt [kNDmNuENCOth]     -> Integral();

  return nexpncbkg;
}
//________________________________________________________________________________________
double SuperkExpectation::ND280Inputs::NExpCC(void) const
{
// expected number of CC events

  double nexpcc =
       fNumSelEt [kNDmNuMuCCQE]     -> Integral() +
       fNumSelEt [kNDmNuMuCCRES]    -> Integral() +
       fNumSelEt [kNDmNuMuCCCOH]    -> Integral() +
       fNumSelEt [kNDmNuMuCCOth]    -> Integral() +
       fNumSelEt [kNDmNuMuBarCCQE]  -> Integral() +
       fNumSelEt [kNDmNuMuBarCCRES] -> Integral() +
       fNumSelEt [kNDmNuMuBarCCCOH] -> Integral() +
       fNumSelEt [kNDmNuMuBarCCOth] -> Integral() +
       fNumSelEt [kNDmNuECCQE]      -> Integral() +
       fNumSelEt [kNDmNuECCRES]     -> Integral() +
       fNumSelEt [kNDmNuECCCOH]     -> Integral() +
       fNumSelEt [kNDmNuECCOth]     -> Integral();

  return nexpcc;
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

  this->NormResolutionMtrx( kSKsFCFV1RingMuLike, kSKmNuMuCCQE    );
  this->NormResolutionMtrx( kSKsFCFV1RingMuLike, kSKmNuMuNonCCQE );
  this->NormResolutionMtrx( kSKsFCFV1RingMuLike, kSKmNC          );

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
void SuperkExpectation::SuperkInputs::NormResolutionMtrx(SuperkSample_t s, SuperkMode_t m)
{
// (Input matrices: y->Etrue, x->Ereco)
//
  TH2D * h2 = fResolMtrx[s][m];

  for(int bint=1; bint <= h2->GetYaxis()->GetNbins(); bint++) {

    double sum = 0;
    for(int binr=1; binr <= h2->GetXaxis()->GetNbins(); binr++) {
       double n  = h2->GetBinContent(binr,bint);
       sum += n;
    }//binr
    if(sum<=0) continue;
    for(int binr=1; binr <= h2->GetXaxis()->GetNbins(); binr++) {
       double n = h2->GetBinContent(binr,bint);
       h2->SetBinContent(binr,bint,n/sum);
    }//binr
    
  }//bint

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
       fResolMtrx[kSKsFCFV1RingMuLike][kSKmNuMuCCQE]->Draw("colz");
       c->cd(2);
       fResolMtrx[kSKsFCFV1RingMuLike][kSKmNuMuNonCCQE]->Draw("colz");
       c->cd(3);
       fResolMtrx[kSKsFCFV1RingMuLike][kSKmNC]->Draw("colz");
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
  int bint = h2->GetYaxis()->FindBin(Etrue);
  int binr = h2->GetXaxis()->FindBin(Ereco);
  double Rij = h2->GetBinContent(binr,bint);
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
  fSuperkFluxPOT = 0;
  fSuperkFlux    = 0; 
  fWght10aTo10d  = 0;
}
//________________________________________________________________________________________
void SuperkExpectation::BeamInputs::CleanUp(void)
{
  if (fSuperkFlux  ) delete fSuperkFlux;
  if (fWght10aTo10d) delete fWght10aTo10d;
}
//________________________________________________________________________________________
bool SuperkExpectation::BeamInputs::Read(TXMLEngine* xml, XMLNodePointer_t root_node)
{
// Read beam inputs to disappearance analysis
//
  cout << "Exp/Beam >> Reading beam inputs from input XML file." << endl;

  // Get SuperK flux
  fSuperkFlux = util::XML2TH1D(xml,root_node,"beam/superk_flux/v10a");
  if(!fSuperkFlux) {
    cerr << "Exp/Beam >> No input SuperK flux histogram!" << endl;
    exit(1);
  }
  // Integrated POT for flux simulation
  fSuperkFluxPOT = util::XML2Dbl(xml,root_node,"beam/superk_flux/pot_mc");
  cout << "Exp/Beam >> SK flux simulation POT = " << fSuperkFluxPOT << endl;

  fSuperkFlux->Scale(1./fSuperkFluxPOT);

  // 10a->10d weights
  fWght10aTo10d = util::XML2TH1D(xml,root_node,"beam/superk_flux/weight_v10a_to_v10d");
  if(!fWght10aTo10d) {
    cerr << "Exp/Beam >> No input 10a->10d SuperK flux weights!" << endl;
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

  if(opt.Contains("Weights/10aTo10d")) { 
    fWght10aTo10d->Draw();
    fWght10aTo10d->GetXaxis()->SetTitle("E_{#nu} (GeV)");
    fWght10aTo10d->GetYaxis()->SetTitle("SuperK flux weight");
    c->Update();
  }
  else
  if(opt.Contains("SKFlux/10a")) { 
    fSuperkFlux->Draw();
    fSuperkFlux->GetXaxis()->SetTitle("E_{#nu} (GeV)");
    fSuperkFlux->GetYaxis()->SetTitle("SuperK flux (#neutrinos /cm2 /POT)");
    c->Update();
  } 
  else
  if(opt.Contains("SKFlux/10d")) { 
    fSuperkFlux->Draw();
    fSuperkFlux->GetXaxis()->SetTitle("E_{#nu} (GeV)");
    fSuperkFlux->GetYaxis()->SetTitle("SuperK flux (#neutrinos /cm2 /POT)");
    c->Update();
  }

}
//________________________________________________________________________________________
double SuperkExpectation::BeamInputs::SuperkFlux (double Etrue, bool wght_to_10d) const
{
  assert(fSuperkFlux);

  int flux_bin = fSuperkFlux->FindBin(Etrue);
  double flux  = fSuperkFlux->GetBinContent(flux_bin);
 
  if(wght_to_10d) {
    int wght_bin = fWght10aTo10d->FindBin(Etrue);
    double wght  = fWght10aTo10d->GetBinContent(wght_bin);
    flux *= wght;
  }

  return flux;
}
//________________________________________________________________________________________
double SuperkExpectation::BeamInputs::SuperkFluxPOT (void) const
{
  return fSuperkFluxPOT;
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

  fErrSigmaModeSK[kSKmNuMuCCQE]    = util::XML2TH1DErrEnv(xml,root_node,"genie/err_sigma_mode_H2O/ccqe");
  fErrSigmaModeSK[kSKmNuMuNonCCQE] = util::XML2TH1DErrEnv(xml,root_node,"genie/err_sigma_mode_H2O/nonccqe");
  fErrSigmaModeSK[kSKmNC]          = util::XML2TH1DErrEnv(xml,root_node,"genie/err_sigma_mode_H2O/nc");

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
  else
  if(opt.Contains("SKErrSigMode")) { 
    TH1F * hframe = (TH1F*)c->DrawFrame(0,-.5,4,.5);
    hframe->Draw();
    fErrSigmaModeSK[kSKmNuMuCCQE]    -> fLow  -> Draw("same");
    fErrSigmaModeSK[kSKmNuMuCCQE]    -> fHigh -> Draw("same");
    fErrSigmaModeSK[kSKmNuMuNonCCQE] -> fLow  -> Draw("same");
    fErrSigmaModeSK[kSKmNuMuNonCCQE] -> fHigh -> Draw("same");
    fErrSigmaModeSK[kSKmNC]          -> fLow  -> Draw("same");
    fErrSigmaModeSK[kSKmNC]          -> fHigh -> Draw("same");
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
  double sigtot   = this->SigmaTotSK(Etrue);
  double frcmode  = this->FractModeSK(Etrue,m);
  double sigmode0 = sigtot * frcmode; 
  double twkdial  = fTwkSigmaModeSK[m];
  double relerr   = this->RelErrSigmaModeSK(Etrue,m,util::Sign(twkdial));
  double sigmode  = sigmode0 * (1. + twkdial * relerr);
  return sigmode;
}
//________________________________________________________________________________________
double SuperkExpectation::GeneratorInputs::RelErrSigmaModeSK(
  double Etrue, SuperkMode_t m, int sign) const
{
// Get estimated relative error on sigma_H2O(mode;Etrue)
//
  if(sign==0) return 0;
  TH1D * hrelerr = (sign>0) ? fErrSigmaModeSK[m]->fHigh : fErrSigmaModeSK[m]->fLow;
  assert(hrelerr);
  int    bin    = hrelerr->FindBin(Etrue);
  double relerr = hrelerr->GetBinContent(bin);
  relerr = TMath::Abs(relerr);
  return relerr;
}
//________________________________________________________________________________________
double SuperkExpectation::GeneratorInputs::SigmaTotCCND280(double Etrue) const
{
//TODO
//CURRENTLY ASSUME A NEAR SK TOO
  return this->SigmaTotSK(Etrue);
}
//________________________________________________________________________________________
void SuperkExpectation::GeneratorInputs::TweakSigmaModeSK(
  t2kana::yr1numu::conv::SuperkMode_t m, double twk_dial)
{
  fTwkSigmaModeSK[m] = twk_dial;
}
//________________________________________________________________________________________

