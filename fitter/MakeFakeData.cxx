//____________________________________________________________________________
/*
 
 Copyright (c) 2010-2014, The VALOR group. All rights reserved.
 
 ChangeLog:
 
 */
//____________________________________________________________________________
#include <fstream>

#include <TROOT.h>
#include <TGraph.h>
#include <TF1.h>
#include <TVector.h>
#include <TMath.h>
#include "MakeFakeData.h"
#include "Utilities.h"
#include <iostream>

using std::cout;
using std::cerr;
using std::endl;
using std::ios;
using std::ifstream;

using namespace valor;

double util::UtlMakeFakeData::WeightMECPDD3D()
{
  
  if(util::IsMEC(m))
  {
    int lepton_int = 3; //lepton is at 3 in array of particles
    
    double dir_nu[3]    = {0.669764, -0.742179, 0.024223}; // T2K beam direction
    double costheta = dir_nu[0] * fBrDirNu[lepton_int][0] + dir_nu[1] * fBrDirNu[lepton_int][1] + dir_nu[2] * fBrDirNu[lepton_int][2];

    
    double pdd_weight =1;
    if(util::SuperkMode2SKMCPdg(m) > 0){   //Neutrino PDD rescale

      int binz = PDD_ratio_nu_3D->GetZaxis()->FindBin(Etrue);
      int binx = PDD_ratio_nu_3D->GetXaxis()->FindBin(costheta);
      int biny = PDD_ratio_nu_3D->GetYaxis()->FindBin(fBrPtrue[lepton_int]);
      
      pdd_weight = PDD_ratio_nu_3D->GetBinContent(binx,biny,binz);
      
    }
    else{         //AntiNeutrino PDD rescale
      int binz = PDD_ratio_nubar_3D->GetZaxis()->FindBin(Etrue);
      int binx = PDD_ratio_nubar_3D->GetXaxis()->FindBin(costheta);
      int biny = PDD_ratio_nubar_3D->GetYaxis()->FindBin(fBrPtrue[lepton_int]);
      
      pdd_weight = PDD_ratio_nubar_3D->GetBinContent(binx,biny,binz);
      

    }
    return pdd_weight;
  }
  
  return 1;

};

double util::UtlMakeFakeData::WeightMECPDD2D()
{
  
  if(util::IsMEC(m))
  {
    int lepton_int = 3; //lepton is at 3 in array of particles
    double dir_nu[3]    = {0.669764, -0.742179, 0.024223}; // T2K beam direction
    double costheta = dir_nu[0] * fBrDirNu[lepton_int][0] + dir_nu[1] * fBrDirNu[lepton_int][1] + dir_nu[2] * fBrDirNu[lepton_int][2];

    double pdd_weight =1;
    if(util::SuperkMode2SKMCPdg(m) > 0){   //Neutrino PDD rescale
      int binx = PDD_ratio_nu_2D->GetXaxis()->FindBin(costheta);
      int biny = PDD_ratio_nu_2D->GetYaxis()->FindBin(fBrPtrue[lepton_int]);
      
      pdd_weight = PDD_ratio_nu_2D->GetBinContent(binx,biny);
      
    }
    else{         //AntiNeutrino PDD rescale
      int binx = PDD_ratio_nubar_2D->GetXaxis()->FindBin(costheta);
      int biny = PDD_ratio_nubar_2D->GetYaxis()->FindBin(fBrPtrue[lepton_int]);
      
      pdd_weight = PDD_ratio_nubar_2D->GetBinContent(binx,biny);
    }
    return pdd_weight;
  }
  
  return 1;
};

double util::UtlMakeFakeData::WeightMartini()
{
  if(util::IsMEC(m))
  {
    double martini_nieves_weight =1;
    if(util::SuperkMode2SKMCPdg(m) > 0){   //Neutrino MEC rescale
      martini_nieves_weight = Martini_to_Nieves_MEC_ratio_Nu->Eval(Etrue);
      if(Etrue > 1.5) martini_nieves_weight = 1;
    }
    else{         //AntiNeutrino MEC rescale
      martini_nieves_weight = Martini_to_Nieves_MEC_ratio_NuBar->Eval(Etrue);
      if(Etrue > 1.5) martini_nieves_weight = 1;
    }
    
    return martini_nieves_weight;
  }
  else
    return 1.;
};

double util::UtlMakeFakeData::EReco1P1H()
{
  double dir_nu[3]    = {0.669764, -0.742179, 0.024223}; // T2K beam direction
  double electronMass = 0.511   / 1000.; // GeV
  double muonMass     = 105.658 / 1000.; // GeV
  int lepton_int      = 2;
  if(util::IsMEC(m))  lepton_int = 3;
  double costheta = dir_nu[0] * fBrDirNu[lepton_int][0] + dir_nu[1] * fBrDirNu[lepton_int][1] + dir_nu[2] * fBrDirNu[lepton_int][2];
  double pmax_lep = 0.; // momentum of most energetic lepton
  double mass_lep = 0.;
  int bin_etrue = NievesMomentumOffset->GetXaxis()->FindBin(Etrue*1000);
  int bin_theta = NievesMomentumOffset->GetYaxis()->FindBin(TMath::ACos(costheta)*180/TMath::Pi());
  double ereco_fill = 0;
  
  if( (s == conv::kSKsFCFV1RingELike) || (s == conv::kSKsRHCFCFV1RingELike) ){ // e-like ring
    mass_lep = electronMass; // GeV
    pmax_lep = fBrAmome[0] / 1000.; // GeV
    pmax_lep += NievesMomentumOffset->GetBinContent(bin_etrue,bin_theta)/1000.;
    if(util::IsRHCSample(s))  ereco_fill = this->ErecoFunc_NuBar(pmax_lep,costheta,mass_lep);
    else                      ereco_fill = this->ErecoFunc_Nu(pmax_lep,costheta,mass_lep);
  }
  else if( s == conv::kSKsFCFVMultiRingELike ){ // MultiRingElike
    mass_lep = electronMass; // GeV
    pmax_lep = fBrAmome[0] / 1000.; // GeV
    pmax_lep += NievesMomentumOffset->GetBinContent(bin_etrue,bin_theta)/1000.;
    ereco_fill = this->ErecoFuncCCOnePi_Nu(pmax_lep,costheta,mass_lep);
  }
  else if( (s == conv::kSKsFCFV1RingMuLike) || (s == conv::kSKsRHCFCFV1RingMuLike) ){ // mu-like ring
    mass_lep = muonMass; // GeV
    pmax_lep = fBrAmomm[0] / 1000.; // GeV
    pmax_lep += NievesMomentumOffset->GetBinContent(bin_etrue,bin_theta)/1000.;
    if(util::IsRHCSample(s))  ereco_fill = this->ErecoFunc_NuBar(pmax_lep,costheta,mass_lep);
    else                      ereco_fill = this->ErecoFunc_Nu(pmax_lep,costheta,mass_lep);
  }
  
  return ereco_fill;

};

double util::UtlMakeFakeData::GetEreco()
{
  switch (reweight_option)
  {
    case util::UtlMakeFakeData::ReweightType::Martini : return ERec;
      break;
    case util::UtlMakeFakeData::ReweightType::MECPDD3D : return ERec;
      break;
    case util::UtlMakeFakeData::ReweightType::MECPDD2D : return ERec;
      break;
    case util::UtlMakeFakeData::ReweightType::Nieves1P1H : return this->EReco1P1H();
      break;
    default:
      return ERec;
  }
}

double util::UtlMakeFakeData::GetWeight()
{
 
  switch (reweight_option)
  {
    case util::UtlMakeFakeData::ReweightType::Martini : return this->WeightMartini();
    break;
    case util::UtlMakeFakeData::ReweightType::MECPDD3D : return this->WeightMECPDD3D();
      break;
    case util::UtlMakeFakeData::ReweightType::MECPDD2D : return this->WeightMECPDD2D();
      break;
    case util::UtlMakeFakeData::ReweightType::Nieves1P1H : return 1.;
      break;
    default:
      return 1.;
  }
  
};

void util::UtlMakeFakeData::SetEventInfo(float _fBrPtrue[10],float _fBrDirNu[10][3],float _Etrue,
                                         conv::SuperkMode_t _m, double _ERec, conv::SuperkSample_t _s,
                                         float _fBrAmome[10], float _fBrAmomm[10] )
{
  ERec = _ERec;
  Etrue = _Etrue;
  s = _s;
  m = _m;
  
  for(int i=0;i<10;i++)
  {
    fBrAmome[i] = _fBrAmome[i];
    fBrAmomm[i] = _fBrAmomm[i];
    fBrPtrue[i] = _fBrPtrue[i];
    
    for(int j=0;j<3;j++)
    {
      fBrDirNu[i][j] = _fBrDirNu[i][j];
    }
  }
  
};

void util::UtlMakeFakeData::SetReweightStyle(string _s)
{
  fReweightString = _s;
  if (!_s.compare("Martini")) {
    reweight_option = util::UtlMakeFakeData::ReweightType::Martini;
  }
  else if (!_s.compare("MECPDD3D")) {
    reweight_option = util::UtlMakeFakeData::ReweightType::MECPDD3D;
  }
  else if (!_s.compare("MECPDD2D")) {
    reweight_option = util::UtlMakeFakeData::ReweightType::MECPDD2D;
  }
  else if (!_s.compare("Nieves1P1H")) {
    reweight_option = util::UtlMakeFakeData::ReweightType::Nieves1P1H;
  }
  else{
    cerr << "UtlMakeFakeData/SetReweightStyle >> Can not match reweight option. " << endl;
    cerr << "Options are: " << endl;
    cerr << "\t Martini" << endl;
    cerr << "\t MECPDD3D" << endl;
    cerr << "\t MECPDD2D" << endl;
    cerr << "\t Nieves1P1H" << endl;
    exit(1);
  }

};

//____________________________________________________________________________________
double util::UtlMakeFakeData::ErecoFunc_Nu(double Mom,double costh,double mass_lep){
  
  const double Eb = 27. / 1000;     // oxygen
  const double mp = 938.272 / 1000.; // proton mass (GeV)
  const double mn = 939.565 / 1000.; // neutron mass (GeV)
  
  double Ereco=0.;
  double E_lep=0.;
  
  E_lep = sqrt( Mom*Mom + mass_lep*mass_lep );    // Reco energy muon or electron
  Ereco = mp*mp - (mn - Eb)*(mn - Eb) - mass_lep*mass_lep + 2*(mn - Eb)*E_lep;   // right formula
  Ereco /= 2*(mn - Eb - E_lep + Mom*costh);
  
  return Ereco;  // GeV
}

double util::UtlMakeFakeData::ErecoFunc_NuBar(double Mom,double costh,double mass_lep){
  
  const double Eb = 27. / 1000;     // oxygen
  const double mp = 938.272 / 1000.; // proton mass (GeV)
  const double mn = 939.565 / 1000.; // neutron mass (GeV)
  
  double Ereco=0.;
  double E_lep=0.;
  
  E_lep = sqrt( Mom*Mom + mass_lep*mass_lep );    // Reco energy lepton
  Ereco = mn*mn - (mp - Eb)*(mp - Eb) - mass_lep*mass_lep + 2*(mp - Eb)*E_lep;
  Ereco /= 2*(mp - Eb - E_lep + Mom*costh);
  
  return Ereco;  // MeV
}

//____________________________________________________________________________________
double util::UtlMakeFakeData::ErecoFuncCCOnePi_Nu(double Mom,double costh,double mass_lep){
  
  const double m_start = 938.272 / 1000.; // proton mass (GeV)  //initial nucluon mass
  const double mn_end = 1232.0 / 1000.; // neutron mass (GeV) //final nucleon mass
  
  double Ereco=0.;
  double E_lep=0.;
  
  E_lep = sqrt( Mom*Mom + mass_lep*mass_lep );    // Reco energy muon or electron
  Ereco = mn_end*mn_end - (m_start)*(m_start) - mass_lep*mass_lep + 2*(m_start)*E_lep;   // right formula
  Ereco /= 2*(m_start - E_lep + Mom*costh);
  
  return Ereco;  // GeV
}

//________________________________________________________________________________________
bool util::UtlMakeFakeData::Read(TXMLEngine* xml, XMLNodePointer_t root_node)
{
#ifdef __show_init_mesg_info__
  cout << "UtlMakeFakeData/Read >> Reading neutrino generator inputs from input XML file." << endl;
#endif
  
  switch (reweight_option)
  {
      
    case util::UtlMakeFakeData::ReweightType::Martini :
      Martini_to_Nieves_MEC_ratio_Nu =
      util::XML2TGraph(xml,root_node,"mc/fake_data/martini_to_nieves_nu");
      Martini_to_Nieves_MEC_ratio_NuBar =
      util::XML2TGraph(xml,root_node,"mc/fake_data/martini_to_nieves_nubar");
      return true;
    case util::UtlMakeFakeData::ReweightType::MECPDD3D :
      PDD_ratio_nu_3D =
      util::XML2TH3D(xml,root_node,"mc/fake_data/three_d/pdd_nu");
      PDD_ratio_nubar_3D =
      util::XML2TH3D(xml,root_node,"mc/fake_data/three_d/pdd_nubar");
      return true;
    case util::UtlMakeFakeData::ReweightType::MECPDD2D :
      PDD_ratio_nu_2D =
      util::XML2TH2D(xml,root_node,"mc/fake_data/two_d/pdd_nu");
      PDD_ratio_nubar_2D =
      util::XML2TH2D(xml,root_node,"mc/fake_data/two_d/pdd_nubar");
      return true;
    case util::UtlMakeFakeData::ReweightType::Nieves1P1H :
      NievesMomentumOffset =
      util::XML2TH2D(xml,root_node,"mc/fake_data/nieves_momentum_offset");
      return true;
    default:
      return false;
  }
}
//____________________________________________________________________________
valor::util::UtlMakeFakeData * valor::util::UtlMakeFakeData::fInstance = 0;
//____________________________________________________________________________
util::UtlMakeFakeData::UtlMakeFakeData()
{
  this->Init();
  fInstance = 0;
}
//____________________________________________________________________________
valor::util::UtlMakeFakeData * valor::util::UtlMakeFakeData::Instance()
{
  if(fInstance == 0) {
    static valor::util::UtlMakeFakeData::Cleaner cleaner;
    cleaner.DummyMethodAndSilentCompiler();
    fInstance = new valor::util::UtlMakeFakeData;
  }
  return fInstance;
}
//____________________________________________________________________________
void util::UtlMakeFakeData::Init(void)
{
  
  Martini_to_Nieves_MEC_ratio_Nu      = 0;
  Martini_to_Nieves_MEC_ratio_NuBar   = 0;
  NievesMomentumOffset                = 0;
  PDD_ratio_nu_2D                     = 0;
  PDD_ratio_nubar_2D                  = 0;
  PDD_ratio_nu_3D                     = 0;
  PDD_ratio_nubar_3D                  = 0;
  reweight_option                     = util::UtlMakeFakeData::ReweightType::FDUndefined;
  fActive                             = 0;
  
  
  //Event info
  Etrue     = 0;
  ERec      = 0;
  m         = conv::kSKmUndefined;
  s         = conv::kSKsUndefined;
  
}