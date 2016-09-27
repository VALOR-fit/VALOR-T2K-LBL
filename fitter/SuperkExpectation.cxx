//________________________________________________________________________________________
/*
 
 \copyright (c) 2010-2013, VaLOR group.

30-12-2013  Nick
Added samples for reverse horn current running

21-7-2014 Davide
Added samples for multiring forward horn current

*/
//________________________________________________________________________________________

#define __show_syst_reset_mesg__
//#define __show_syst_tweak_mesg__
#define __show_syst_tweak_debug__
#define __show_physparamid__

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
#include "Utilities.h"
#include "valor/Conventions/Opt.h"
#include "valor/Utilities/Rnd.h"
#include "valor/Utilities/UtlHist.h"
#include "valor/Messenger/Messenger.h"

using std::cout;
using std::cerr;
using std::endl;

using namespace valor;
using namespace valor::conv;

ClassImp(SuperkExpectation)

//________________________________________________________________________________________
SuperkExpectation::SuperkExpectation()
{
  this->Init();
  this->CheckPhysParamId(); 
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
  
  for (int is = 0 ; is < conv::kns ; is++) {
    SuperkSample_t s = conv::kSKSamples[is];
    if (s==kSKsUndefined) {
      fPredictErecoSpectrumArrayTot [is] = 0;
      fPredictErecoSpectrumArray    [is] = 0;
      // NEW 2D TEMPLATES
      fPredictErecoVsQ2SpectrumArrayTot [is] = 0;
      fPredictErecoVsQ2SpectrumArray    [is] = 0;
      //
      continue;
    }
    fPredictErecoSpectrumArrayTot [is]= new VALORArray(util::NumErecoBins(s));
    fPredictErecoSpectrumArray    [is]= new VALORArray(conv::knm, util::NumErecoBins(s));
    
    fPredictErecoSpectrumArrayTot[is] ->SetDimensionName(0, "Ereco");
    fPredictErecoSpectrumArray   [is] ->SetDimensionName(0, "SK Mode");
    fPredictErecoSpectrumArray   [is] ->SetDimensionName(1, "Ereco");
    
    fPredictErecoSpectrumArrayTot[is]->SetTitle(Form("Predicted Ereco Spectrum Array %s total"  ,util::AsStringCompact(s).c_str()));
    fPredictErecoSpectrumArray   [is]->SetTitle(Form("Predicted Ereco Spectrum Array %s by mode",util::AsStringCompact(s).c_str()));
    fPredictErecoSpectrumArrayTot[is]->SetName(Form("ereco_%s_total"  , util::AsStringCompact(s).c_str()));
    fPredictErecoSpectrumArray   [is]->SetName(Form("ereco_%s_mode"   , util::AsStringCompact(s).c_str()));

    // NEW 2D TEMPLATES
    fPredictErecoVsQ2SpectrumArrayTot [is]= new VALORArray(util::NumErecoBins(s),util::NumQ2Bins(s));
    fPredictErecoVsQ2SpectrumArray    [is]= new VALORArray(conv::knm, util::NumErecoBins(s), util::NumQ2Bins(s));
    
    fPredictErecoVsQ2SpectrumArrayTot[is] ->SetDimensionName(0, "Ereco");
    fPredictErecoVsQ2SpectrumArrayTot[is] ->SetDimensionName(1, "Q2");

    fPredictErecoVsQ2SpectrumArray   [is] ->SetDimensionName(0, "SK Mode");
    fPredictErecoVsQ2SpectrumArray   [is] ->SetDimensionName(1, "Ereco");
    fPredictErecoVsQ2SpectrumArray   [is] ->SetDimensionName(2, "Q2");
    
    fPredictErecoVsQ2SpectrumArrayTot[is]->SetTitle(Form("Predicted Ereco Vs Q2 Spectrum Array %s total"  ,util::AsStringCompact(s).c_str()));
    fPredictErecoVsQ2SpectrumArray   [is]->SetTitle(Form("Predicted Ereco Vs Q2 Spectrum Array %s by mode",util::AsStringCompact(s).c_str()));
    fPredictErecoVsQ2SpectrumArrayTot[is]->SetName(Form("ereco_Vs_q2_%s_total"  , util::AsStringCompact(s).c_str()));
    fPredictErecoVsQ2SpectrumArray   [is]->SetName(Form("ereco_Vs_q2_%s_mode"   , util::AsStringCompact(s).c_str()));
    //
  }
  
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
  if (fSuperkMC          ) delete fSuperkMC;
  if (fNuFlux            ) delete fNuFlux; 
  if (fNuXSec            ) delete fNuXSec;
  if (fBANFF             ) delete fBANFF;
}
//________________________________________________________________________________________
void SuperkExpectation::CheckPhysParamId(void)
{
  // Check PhysParamId_t
    
  LOG("SuperkExpectation", pINFO) << "Check FHC BANFF PhysParamId"; 
  for (int i=0; i<conv::kNumOfParamsInBANFFfit;i++){
    PhysParamId_t p = (conv::PhysParamId_t)(i+conv::kNuisBANFF00);
    PhysParamId_t p_list = ( PhysParamId_t)kListOfBANFFParams[i];
#ifdef __show_physparamid__
    std::cout << util::AsString(p) << endl;
#endif
    std::string parname = util::AsString(p);
    if( parname == "--"){
      std::cerr << "SuperkExpectation/CheckPhysParamId >> parameter " << p << " is not defined:" 
        	<< parname << endl;
      util::Die(1);
    }
    if(util::AsString(p_list) == "--"){
      std::cerr << "SuperkExpectation/CheckPhysParamId >> parameter " << p_list << " in kListOfBANFFParams is not defined:" 
    		<< util::AsString(p_list) << endl;
      util::Die(1);
    }    
  }
  
  LOG("SuperkExpectation", pINFO) << "Check RHC BANFF PhysParamId"; 
  for (int i=0; i<kNumOfRHCParamsInBANFFfit;i++){
    PhysParamId_t p = (conv::PhysParamId_t)(i+conv::kRHCNuisBANFF00);
    PhysParamId_t p_list = ( PhysParamId_t)kListOfRHCBANFFParams[i];
#ifdef __show_physparamid__
    std::cout << util::AsString(p) << endl;
#endif
    if(util::AsString(p) == "--"){
      std::cerr << "SuperkExpectation/CheckPhysParamId >> parameter " << p << " is not defined:" 
        	<< util::AsString(p) << endl;
      util::Die(1);
    }
    if(util::AsString(p_list) == "--"){
      std::cerr << "SuperkExpectation/CheckPhysParamId >> parameter " << p_list << " in kListOfRHCBANFFParams is not defined:" 
        	<< util::AsString(p_list) << endl;
      util::Die(1);
    }    
  }
  for (int i=0; i<kNumOfFHCpRHCParamsInBANFFfit;i++){
    PhysParamId_t p_list_rhcpfhc = ( PhysParamId_t)kListOfFHCpRHCBANFFParams[i];
    if(util::AsString(p_list_rhcpfhc) == "--"){
      std::cerr << "SuperkExpectation/CheckPhysParamId >> parameter " << p_list_rhcpfhc << " in kListOfFHCpRHCBANFFParams is not defined:" 
        	<< util::AsString(p_list_rhcpfhc) << endl;
      util::Die(1);
    }    
  }

  LOG("SuperkExpectation", pINFO) << "Check FHC XSec PhysParamId"; 
  for (int i=0; i<kNumOfUncXSecSyst;i++){
    PhysParamId_t p = (conv::PhysParamId_t)(i+conv::kLastBANFFParam+1);
    PhysParamId_t p_list = (conv::PhysParamId_t)kListOfXSecSyst[i];
#ifdef __show_physparamid__
    std::cout << util::AsString(p) << endl;
#endif
    if(util::AsString(p) == "--"){
      std::cerr << "SuperkExpectation/CheckPhysParamId >> parameter " << p << " is not defined:" 
        	<< util::AsString(p) << endl;
      util::Die(1);
    }
    if(util::AsString(p_list) == "--"){
      std::cerr << "SuperkExpectation/CheckPhysParamId >> parameter " << p_list << " in kListOfXSecSyst is not defined:" 
        	<< util::AsString(p_list) << endl;
      util::Die(1);
    }    
  }

  LOG("SuperkExpectation", pINFO) << "Check RHC XSec PhysParamId"; 
  for (int i=0; i<kNumOfRHCUncXSecSyst;i++){
    PhysParamId_t p = (conv::PhysParamId_t)(i+conv::kRHCLastBANFFParam+1);
#ifdef __show_physparamid__
    std::cout << util::AsString(p) << endl;
#endif
    if(util::AsString(p) == "--"){
      std::cerr << "SuperkExpectation/CheckPhysParamId >> parameter " << p << " is not defined:" 
        	<< util::AsString(p) << endl;
      util::Die(1);
    }
  }
  for (int i=0; i<kNumOfRHCXSecSyst;i++){
    PhysParamId_t p = (conv::PhysParamId_t)kListOfRHCXSecSyst[i];
    if(util::AsString(p) == "--"){
      std::cerr << "SuperkExpectation/CheckPhysParamId >> parameter " << p << " in kListOfRHCXSecSyst is not defined:" 
        	<< util::AsString(p) << endl;
      util::Die(1);
    }    
  }

  LOG("SuperkExpectation", pINFO) << "Check FHC SuperkDetectorFSI PhysParamId"; 
  for (int i=0; i<kNumOfSKDetFSIParams;i++){
    PhysParamId_t p = (conv::PhysParamId_t)(i+conv::kNuisSKDetFSI000);
    PhysParamId_t p_list = (conv::PhysParamId_t)kListOfSKDetFSIParams[i];
#ifdef __show_physparamid__
    std::cout << util::AsString(p) << endl;
#endif
    if(util::AsString(p) == "--"){
      std::cerr << "SuperkExpectation/CheckPhysParamId >> parameter " << p << " is not defined:" 
        	<< util::AsString(p) << endl;
      util::Die(1);
    }
    if(util::AsString(p_list) == "--"){
      std::cerr << "SuperkExpectation/CheckPhysParamId >> parameter " << p_list << " in kListOfSKDetFSIParams is not defined:" 
		<< util::AsString(p_list) << endl;
      util::Die(1);
    }    
  }

   LOG("SuperkExpectation", pINFO) << "Check FHC SuperkDetectorFSI MultiRing PhysParamId"; 
  for (int i=0; i<kNumOfSKDetFSIParamsMultiRing;i++){
    // Check for p in MultiRingNuis doesn't currently (2016/04/25) work since there are jumps in the parameters
    // PhysParamId_t p = (PhysParamId_t)(i+conv::kMultiRingNuisSKDetFSI_first);
    PhysParamId_t p_list = (PhysParamId_t)kListOfMultiRingSKDetFSIParams[i];
#ifdef __show_physparamid__
    std::cout << p_list << " " << util::AsString(p_list) << endl;
#endif
    // if(util::AsString(p) == "--"){
    //   std::cerr << "SuperkExpectation/CheckPhysParamId >> parameter " << p << " is not defined:" 
    //     	<< util::AsString(p) << endl;
    //   util::Die(1);
    // }
    if(util::AsString(p_list) == "--"){
      std::cerr << "SuperkExpectation/CheckPhysParamId >> parameter " << p_list << " in kListOfSKDetFSIParams is not defined:" 
        	<< util::AsString(p_list) << endl;
      util::Die(1);
    }    
  }
  
  LOG("SuperkExpectation", pINFO) << "Check RHC SuperkDetectorFSI PhysParamId"; 
  for (int i=0; i<kNumOfRHCSKDetFSIParams;i++){
    PhysParamId_t p = (conv::PhysParamId_t)(i+conv::kRHCNuisSKDetFSI000);
    PhysParamId_t p_list = (conv::PhysParamId_t)kListOfRHCSKDetFSIParams[i];
#ifdef __show_physparamid__
    std::cout << util::AsString(p) << endl;
#endif
    if(util::AsString(p) == "--"){
      std::cerr << "SuperkExpectation/CheckPhysParamId >> parameter " << p << " is not defined:" 
        	<< util::AsString(p) << endl;
      util::Die(1);
    }
    if(util::AsString(p_list) == "--"){
      std::cerr << "SuperkExpectation/CheckPhysParamId >> parameter " << p_list << " in kListOfRHCSKDetFSIParams is not defined:" 
        	<< util::AsString(p_list) << endl;
      util::Die(1);
    }    
  }
 for (int i=0; i<kNumOfFHCpRHCSKDetFSIParams;i++){
    PhysParamId_t p = ( PhysParamId_t)kListOfFHCpRHCSKDetFSIParams[i];
    if(util::AsString(p) == "--"){
      std::cerr << "SuperkExpectation/CheckPhysParamId >> parameter " << p << " in kListOfFHCpRHCSKDetFSIParams is not defined:" 
        	<< util::AsString(p) << endl;
      util::Die(1);
    }    
  }
  
 LOG("SuperkExpectation", pINFO) << " Done check of PhysParamId";
}
//________________________________________________________________________________________
bool SuperkExpectation::Read(string config_file, Option_t * /*opt*/)
{  
  if (T2KOpt::Instance()->ReadXmlPath(config_file)) {
    cerr << "SKExp/Read >> Unable to access config file: " << config_file << endl;
    return false;
  }
  
  cout << "SKExp/Read >> Reading config file: " << config_file << endl;
  
  if(fIsConfigured) {
    this->CleanUp();
  } 

  // Instantiate all necessary histograms
  this->Book();

  // Create components used for building-up the p.d.f & implemeting its systematic variation
  fBANFF             = new BANFF;
  fNuFlux            = new NuFlux; 
  fNuXSec            = new NuXSec; 
  fSuperkDetectorFSI = new SuperkDetectorFSI;
  fSuperkMC          = new SuperkMC;

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
    cerr << "SKExp/Read >> " << "Can not parse input: " << config_file << endl;
    delete xml;
    return false;  
  }
  XMLNodePointer_t root_node = xml->DocGetRootElement(xmldoc);

  // Configure each of the components 
  
  if ( ! this->GetBANFF() -> Read (xml, root_node) ) return false;
  if ( ! this->GetSuperkDetectorFSI() -> Read (xml, root_node) ) return false;

  // NEW FLUX FHC+RHC
#if defined(__FHC_RUNNING__) && defined(__RHC_RUNNING__) && !defined(__FULLY_CORRELATED__)
  
  if ( ! this->GetNuFlux           () -> Read (xml, root_node) ) return false;
  if ( ! this->GetNuFlux           () -> ReadRHC (xml, root_node) ) return false;
  if ( ! this->GetNuXSec           () -> Read (xml, root_node) ) return false;

  // NEW FHC+RHC JOINT
  
#else
  if ( ! this->GetNuFlux           () -> Read (xml, root_node) ) return false;
  if ( ! this->GetNuFlux           () -> ReadRHC (xml, root_node) ) return false;
  if ( ! this->GetNuXSec           () -> Read (xml, root_node) ) return false;
#endif

  if ( ! this->GetSuperkMC         () -> Read (xml, root_node) ) return false;

  fIsConfigured = true;

  // Release memory
  xml->FreeDoc(xmldoc);
  delete xml;

  LOG("SuperkExpectation",pINFO) << "Done configuring Superk spectrum predictor... \n\n";
  
  return true;
}
//________________________________________________________________________________________
//Config file means data file here. Function reads dataset. Copied from SKObs
//________________________________________________________________________________________
bool SuperkExpectation::Read(string config_file, int data_set_id, bool is_xml)
{
  if (T2KOpt::Instance()->ReadXmlPath(config_file)) {
    cerr << "SKExp/Read >> Unable to access data file: " << config_file << endl;
    return false;
  }
  
  cout << "Exp >> " << "Reading data file: " << config_file << endl;
  
  if(is_xml) {
    // Parse the input configuration file
    TXMLEngine* xml = new TXMLEngine;
    XMLDocPointer_t xmldoc = xml->ParseFile(config_file.c_str());
    if (xmldoc==0) {
      cerr << "SKExp/Read >> " << "Can not parse input: " << config_file << endl;
      delete xml;
      exit(1);
    }
    XMLNodePointer_t root_node = xml->DocGetRootElement(xmldoc);
    
    for (int is = 0 ; is < conv::kns ; is++) {
      SuperkSample_t s = conv::kSKSamples[is];
      if (s==kSKsUndefined) continue;
#ifdef __show_read_events_mesg_debug__
      cout << "Exp >> Read: Reading events for sample: " << util::AsString(s) << endl;
#endif
      std::ostringstream datanode;
      datanode << "superk/observation";
      datanode << "/exp_data_" << data_set_id;
      datanode << "/"<<util::GetXMLDatasetName(s);
      
#ifndef __USE_2D_TEMPLATES__ //1D templates
      datanode << "/"<<util::GetXMLDatasetParamsName();
      
      if(util::NodeExists(xml,root_node,datanode.str())) {
        TH1D h("", "", util::NumErecoBins(s), util::ErecoBinEdges(s));
        h.SetDirectory(0);
        vector<double> values = util::XML2DblArray(xml,root_node,datanode.str());
        vector<double>::const_iterator value_iter = values.begin();
        for( ; value_iter != values.end(); ++value_iter) {
          double Ereco = *value_iter;
          h.Fill(Ereco);
        }
        cout << "Exp >> Read: " << h.Integral()
        << " " << util::AsString(s) << " events." << endl;
        for(int ibin = 1; ibin <= util::NumErecoBins(s); ibin++) {
          double nexp = h.GetBinContent(ibin);
          fPredictErecoSpectrumArrayTot[is]->GetArray()[ibin-1] = nexp;
        }
      }
      else {
        cerr << "Exp >> Node " << datanode.str() << " does not exist!" << endl;
      }
#else //2D templates
      //first get the number of events
      int nevents = 0;
      {
        std::ostringstream datanodetmp;
        datanodetmp << datanode.str();
        datanodetmp << "/nevents";
        if(util::NodeExists(xml,root_node,datanodetmp.str())) {
          nevents = util::XML2Int(xml,root_node,datanodetmp.str());
        }
        else {
          cerr << "Exp >> Node " << datanodetmp.str() << " does not exist!" << endl;
        }
      }
      //now get each event
      datanode << "/"<<util::GetXMLDatasetParamsName();
      TH2D h("", "", util::NumErecoBins(s), util::ErecoBinEdges(s), util::NumQ2Bins(s), util::Q2BinEdges(s));
      h.SetDirectory(0);
      for(int iev = 0; iev < nevents; iev++) {
        std::ostringstream datanodetmp;
        datanodetmp << datanode.str();
        datanodetmp << "/ev" << iev;
        if(util::NodeExists(xml,root_node,datanodetmp.str())) {
          pair<double,double> values = util::XML2DblDblPair(xml,root_node,datanodetmp.str());
          double Ereco = values.first;
          double Q2    = values.second;
#ifdef __show_read_events_mesg_debug__
          cout << "Exp >> Read: event " << iev
          << " " << fPredictErecoVsQ2SpectrumArrayTot[is]->GetDimensionName(0) << " = " << Ereco
          << " " << fPredictErecoVsQ2SpectrumArrayTot[is]->GetDimensionName(1) << " = " << Q2
          << endl;
#endif
          h.Fill(Ereco, Q2);
        }
        else {
          cerr << "Exp >> Node " << datanodetmp.str() << " does not exist!" << endl;
        }
      }//iev
      cout << "Exp >> Read: " << h.Integral()
      << " " << util::AsString(s) << " events." << endl;
      for(int ibinx = 1; ibinx <= util::NumErecoBins(s); ibinx++) {
        for(int ibiny = 1; ibiny <= util::NumQ2Bins(s); ibiny++) {
          double nExp = h.GetBinContent(ibinx, ibiny);
          fPredictErecoVsQ2SpectrumArrayTot[is]->SetEntry(nExp, ibinx-1, ibiny-1);
        }//ibinx
      }//ibiny
#endif //__USE_2D_TEMPLATES__
      
    } // is
    
    // Release memory
    xml->FreeDoc(xmldoc);
    delete xml;
  }//is_xml
  
  else {
    
    TFile f(config_file.c_str());
    
#ifndef __USE_2D_TEMPLATES__ //1D templates
    TH1D * hist = 0;
#else
    TH2D * hist = 0;
#endif
    
    if(!f.IsOpen()){
      cerr << endl;
      cerr << "File " << config_file.c_str() << " is not defined" << endl;
      cerr << endl;
      exit(1);
    }
    
    for (int is = 0 ; is < conv::kns ; is++) {
      SuperkSample_t s = conv::kSKSamples[is];
      
      //select which histograms to read in
#if !defined(__FHC_RUNNING__)
      //skip the FHC samples
      if(s == kSKsFCFV1RingMuLike || s == kSKsFCFV1RingELike) continue;
#endif
#if !defined(__RHC_RUNNING__)
      //skip the RHC samples
      if(s == kSKsRHCFCFV1RingMuLike || s == kSKsRHCFCFV1RingELike) continue;
#endif
#if defined(__NUEONLY__)
      //skip the mu-like samples
      if(s == kSKsFCFV1RingMuLike || s == kSKsRHCFCFV1RingMuLike) continue;
#endif
#if defined(__NUMUONLY__)
      //skip the e-like samples
      if(s == kSKsFCFV1RingELike || s == kSKsRHCFCFV1RingELike) continue;
#endif
      
      if (s == conv::kSKsUndefined) continue;
      
#ifndef __USE_2D_TEMPLATES__ //1D templates
      else if(s == conv::kSKsRHCFCFV1RingMuLike) hist = (TH1D*) f.Get(Form("exp_data_RHC1Rmu_%d",data_set_id));
      else if(s == conv::kSKsRHCFCFV1RingELike)	 hist = (TH1D*) f.Get(Form("exp_data_RHC1Re_%d",data_set_id));
      else if(s == conv::kSKsFCFV1RingELike)   	 hist = (TH1D*) f.Get(Form("exp_data_1Re_%d",data_set_id));
      else if(s == conv::kSKsFCFV1RingMuLike)	   hist = (TH1D*) f.Get(Form("exp_data_1Rmu_%d",data_set_id));
      else if(s == conv::kSKsFCFVMultiRingELike) hist = (TH1D*) f.Get(Form("exp_data_MultiRe_%d",data_set_id));
#else //2D templates
      else if(s == conv::kSKsRHCFCFV1RingMuLike) hist = (TH2D*) f.Get(Form("exp_data_RHC1Rmu_%d",data_set_id));
      else if(s == conv::kSKsRHCFCFV1RingELike)	 hist = (TH2D*) f.Get(Form("exp_data_RHC1Re_%d",data_set_id));
      else if(s == conv::kSKsFCFV1RingELike)	   hist = (TH2D*) f.Get(Form("exp_data_1Re_%d",data_set_id));
      else if(s == conv::kSKsFCFV1RingMuLike) 	 hist = (TH2D*) f.Get(Form("exp_data_1Rmu_%d",data_set_id));
      else if(s == conv::kSKsFCFVMultiRingELike) hist = (TH2D*) f.Get(Form("exp_data_MultiRe_%d",data_set_id));
#endif
      else {
        cerr << endl;
        cerr << "Used SK sample " << util::AsString(s) << " not available" << endl;
        cerr << endl;
        exit(1);
      }
      
      if(!hist) {
        cerr << "Exp >> Histogram for " << util::AsString(s) << " does not exist!" << endl;
        if(s == conv::kSKsFCFV1RingMuLike)
        cerr << "Exp >> Need histogram name " << Form("exp_data_1Rmu_%d",data_set_id) << endl;
        else if(s == conv::kSKsRHCFCFV1RingMuLike)
        cerr << "Exp >> Need histogram name " << Form("exp_data_RHC1Rmu_%d",data_set_id) << endl;
        else if(s == conv::kSKsFCFV1RingELike)
        cerr << "Exp >> Need histogram name " << Form("exp_data_1Re_%d",data_set_id) << endl;
        else if(s == conv::kSKsRHCFCFV1RingELike)
        cerr << "Exp >> Need histogram name " << Form("exp_data_RHC1Re_%d",data_set_id) << endl;
        else if(s == conv::kSKsFCFVMultiRingELike)
        cerr << "Exp >> Need histogram name " << Form("mock_data_MultiRe_%d",data_set_id) << endl;
        else{
          cerr << "Exp >> Histogram sample not defined" << endl;
          exit(1);
        }
        exit(1);
      }
      
      
#ifndef __USE_2D_TEMPLATES__ //1D templates
      for(int ibin = 1; ibin <= util::NumErecoBins(s); ibin++) {
        double nExp = hist->GetBinContent(ibin);
        fPredictErecoSpectrumArrayTot[is]->GetArray()[ibin-1] = nExp;
      }
#else //2D templates
#ifdef __show_read_events_mesg_debug__
      cout << "Exp >> Read: " << hist->Integral()
      << " " << util::AsString(s) << " events." << endl;
#endif
      for(int ibinx = 1; ibinx <= util::NumErecoBins(s); ibinx++) {
        for(int ibiny = 1; ibiny <= util::NumQ2Bins(s); ibiny++) {
          // histogram convention (Q2,Ereco) is opposite of VALORArray's one (Ereco,Q2)
          double nExp = hist->GetBinContent(ibiny, ibinx);
          fPredictErecoVsQ2SpectrumArrayTot[is]->SetEntry(nExp, ibinx-1, ibiny-1);
        }//ibiny
      }//ibinx
      
#endif //__USE_2D_TEMPLATES__
      
    }//is
  }//!is_xml
  
  return true;
}
//________________________________________________________________________________________
bool SuperkExpectation::Generate(const PhysHypothesis & ph, Option_t * /*opt*/)
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
  fSuperkMC->UpdateTemplates(); // it already supports both 1D and 2D templates

#ifndef __USE_2D_TEMPLATES__ // 1D templates

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

#else // 2D templates

  // NEW 2D TEMPLATES
  
  // update Ereco spectrum prediction
  bool applied_negevents_condition = fSuperkMC->UpdateErecoVsQ2Spectra();

  // copy histograms
  for(int is = 0; is < kns; is++) {
    SuperkSample_t s = kSKSamples[is];
    if(s == kSKsUndefined) continue;

    const int nbins_q2    = kNQ2Bins[is];
    const int nbins_ereco = kNErecoBins[is];

    typedef double (erq2type)  [nbins_ereco][nbins_q2];
    const erq2type & arr_skmc     = *reinterpret_cast< const erq2type * > ( fSuperkMC->CurrErecoVsQ2SpectrumArray(s) );
    if (!arr_skmc) return 0;

    erq2type & arr_skexp    = *reinterpret_cast< erq2type * > ( this->CurrPredictedErecoVsQ2SpectrumArray(s) );
    if (!arr_skexp) return 0;
    
    for(int ibin_ereco = 1; ibin_ereco <= nbins_ereco; ibin_ereco++) {
      for(int ibin_q2 = 1; ibin_q2 <= nbins_q2; ibin_q2++) {
        
        //double ns = fSuperkMC->CurrErecoVsQ2SpectrumArray(s)[ibin_ereco-1][ibin_q2-1];
        //this->CurrPredictedErecoVsQ2SpectrumArray(s)[ibin_ereco-1][ibin_q2-1] = ns;

        double ns = arr_skmc[ibin_ereco-1][ibin_q2-1];
        arr_skexp[ibin_ereco-1][ibin_q2-1] = ns;

        // now per mode
 
        for(int im = 0; im < knm; im++) {
          SuperkMode_t m = kSKReactionModes[is][im];
          if(m == kSKmUndefined) continue;

          const erq2type & arr_skmc_mode = *reinterpret_cast< const erq2type * > (fSuperkMC->CurrErecoVsQ2SpectrumArray(s,m));  
          if (!arr_skmc_mode) return 0;
          
          erq2type & arr_skexp_mode = *reinterpret_cast< erq2type * > (this->CurrPredictedErecoVsQ2SpectrumArrayMode(s,m));  
          if (!arr_skexp_mode) return 0;
          
          //double nsm = fSuperkMC->CurrErecoVsQ2SpectrumArray(s,m)[ibin_ereco-1][ibin_q2-1]; 
          //this->CurrPredictedErecoVsQ2SpectrumArrayMode(s,m)[ibin_ereco-1][ibin_q2-1] = nsm; 
        
          double nsm = arr_skmc_mode[ibin_ereco-1][ibin_q2-1]; 
          arr_skexp_mode[ibin_ereco-1][ibin_q2-1] = nsm; 

        }//m
      }//ibin q2
    }//ibin ereco
  }//s
  
  //
  
  
  
  
#endif

  return applied_negevents_condition;
}
//________________________________________________________________________________________
bool SuperkExpectation::Generate(const FitParamList & fp, Option_t * /*opt*/)
{
  return Generate(fp, false, false);
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
TH1D* SuperkExpectation::CurrPredictedErecoSpectrum(valor::conv::SuperkSample_t s) const 
{
  const int dimension = 0 ; // 0 = Ereco, 1 = Q2, theta, etc
  const VALORArray * va = 0;
#ifdef __USE_2D_TEMPLATES__
  va = this->CurrPredictedErecoVsQ2VArray(s);
#else
  va = this->CurrPredictedErecoVArray(s);
#endif
  if (!va) {
    std::cerr<<"SuperkExpectation/CurrPredictedErecoSpectrum(s) >> Failed to get VALORArray for "<<util::AsString(s)<<endl;
    return 0;
  }
  
  string name (Form("curr_ereco_%s"                       ,util::AsStringCompact(s).c_str() ) );
  string title(Form("Current Predicted Ereco Spectrum: %s",util::AsString(s).c_str() ) );
  
  const double * bin_edges = util::ErecoBinEdges(s);
  
  TH1D * h = new TH1D(util::Make1DHistogramFromArray(*va,bin_edges,dimension));
  util::FormatErecoHist(h);
  h->SetName(name.c_str());
  h->SetTitle(title.c_str());
  
  return h;
}
//________________________________________________________________________________________
TH1D* SuperkExpectation::CurrPredictedErecoSpectrum(
  valor::conv::SuperkSample_t s, valor::conv::SuperkMode_t m) const 
{
  const int dimension = 0 ; // 0 = Ereco, 1 = Q2, theta, etc
  const VALORArray * va = 0;
#ifdef __USE_2D_TEMPLATES__
  va = this->CurrPredictedErecoVsQ2VArrayMode(s);
#else
  va = this->CurrPredictedErecoVArrayMode(s);
#endif
  if (!va) {
    std::cerr<<"SuperkExpectation/CurrPredictedErecoSpectrum(s,m) >> Failed to get VALORArray for "<<util::AsString(s)<<" - "<<util::AsString(m)<<endl;
    return 0;
  }
  
  int midx = util::ArrIdx::Instance()->SuperkModeIdx(s,m);
  if(midx == -1) {
    std::cerr<<"SuperkExpectation/CurrPredictedErecoSpectrum(s,m) >> Failed to get mode for "<<util::AsString(s)<<" - "<<util::AsString(m)<<endl;
    return 0;
  }
  
  string name  (Form("curr_ereco_%s_%s"                       ,util::AsStringCompact(s).c_str(), util::AsStringCompact(m).c_str()) );
  string title (Form("Current Predicted Ereco Spectrum: %s %s",util::AsString(s).c_str()       , util::AsString(m).c_str()       ) );
  const double * bin_edges = util::ErecoBinEdges(s);
    
  TH1D * h = new TH1D(util::Make1DHistogramFromArray(va->GetSlice(0,midx,midx),bin_edges,dimension));
  util::FormatErecoHist(h);
  util::SetDefaultSKSpectrumStyle(s,m,h);
  h->SetName(name.c_str());
  h->SetTitle(title.c_str());
    
  return h;
}
//________________________________________________________________________________________
VALORArray * SuperkExpectation::CurrPredictedErecoVArray(valor::conv::SuperkSample_t s) 
{
  VALORArray * arr = fPredictErecoSpectrumArrayTot[s];
  if (!arr) {
    std::cerr<<"SKExp/CurrPredictedErecoVArray could not find array for sample "<<util::AsString(s)<<std::endl;
    exit(-1);
  }
  
  return arr;
}
//________________________________________________________________________________________
VALORArray * SuperkExpectation::CurrPredictedErecoVArrayMode(valor::conv::SuperkSample_t s) 
{
  VALORArray * arr = fPredictErecoSpectrumArray[s];
  if (!arr) {
    std::cerr<<"SKExp/CurrPredictedErecoVArrayMode could not find array for sample "<<util::AsString(s)<<std::endl;
    exit(-1);
  }
  
  return arr;
}
//________________________________________________________________________________________
const VALORArray * SuperkExpectation::CurrPredictedErecoVArray(valor::conv::SuperkSample_t s) const
{
  const VALORArray * arr = fPredictErecoSpectrumArrayTot[s];
  if (!arr) {
    std::cerr<<"SKExp/CurrPredictedErecoVArray could not find array for sample "<<util::AsString(s)<<std::endl;
    exit(-1);
  }
  
  return arr;
}
//________________________________________________________________________________________
const VALORArray * SuperkExpectation::CurrPredictedErecoVArrayMode(valor::conv::SuperkSample_t s) const
{
  const VALORArray * arr = fPredictErecoSpectrumArray[s];
  if (!arr) {
    std::cerr<<"SKExp/CurrPredictedErecoVArrayMode could not find array for sample "<<util::AsString(s)<<std::endl;
    exit(-1);
  }
  
  return arr;
}
//________________________________________________________________________________________
double* SuperkExpectation::CurrPredictedErecoSpectrumArray(valor::conv::SuperkSample_t s) 
{
  VALORArray * arr = this->CurrPredictedErecoVArray(s);
  return arr->GetArray();
}
//________________________________________________________________________________________
const double* SuperkExpectation::CurrPredictedErecoSpectrumArray(valor::conv::SuperkSample_t s)  const
{
  const VALORArray * arr = this->CurrPredictedErecoVArray(s);
  return arr->GetArray();
}
//________________________________________________________________________________________
double* SuperkExpectation::CurrPredictedErecoSpectrumArrayMode(
  valor::conv::SuperkSample_t s, valor::conv::SuperkMode_t m)
{ 
  int sidx = util::ArrIdx::Instance()->SuperkSampleIdx(s);
  if(sidx == -1) {
    return 0;
  }
  int midx = util::ArrIdx::Instance()->SuperkModeIdx(s,m);
  if(midx == -1) {
    return 0;
  }
  
  VALORArray * arr = this->CurrPredictedErecoVArrayMode(s);

  return &(arr->GetEntryRef(midx,0));
}
//________________________________________________________________________________________
const double* SuperkExpectation::CurrPredictedErecoSpectrumArrayMode(
  valor::conv::SuperkSample_t s, valor::conv::SuperkMode_t m) const
{ 
  int sidx = util::ArrIdx::Instance()->SuperkSampleIdx(s);
  if(sidx == -1) {
    return 0;
  }
  int midx = util::ArrIdx::Instance()->SuperkModeIdx(s,m);
  if(midx == -1) {
    return 0;
  }
  
  const VALORArray * arr = this->CurrPredictedErecoVArrayMode(s);

  return &(arr->GetEntryRef(midx,0));
}
//________________________________________________________________________________________
double SuperkExpectation::CurrPredictedErecoSpectrumArray(valor::conv::SuperkSample_t s, int i) const
{
  const double * arr = CurrPredictedErecoSpectrumArray(s);
  assert(arr);
  return arr[i];
}
//________________________________________________________________________________________


// NEW 2D TEMPLATES

TH2D* SuperkExpectation::CurrPredictedErecoVsQ2Spectrum(valor::conv::SuperkSample_t s) const 
{
  const int nbins_q2    = util::NumQ2Bins(s);
  const int nbins_ereco = util::NumErecoBins(s);
  typedef double (erq2type)  [nbins_ereco][nbins_q2];  

  const erq2type & arr = *reinterpret_cast< const erq2type * > (CurrPredictedErecoVsQ2SpectrumArray(s));  
  //const double * arr = CurrPredictedErecoSpectrumArray(s);
  if (!arr) return 0;
  
  string name (Form("curr_ereco_Vs_q2_%s"                       ,util::AsStringCompact(s).c_str() ) );
  string title(Form("Current Predicted Ereco Vs Q2 Spectrum: %s",util::AsString(s).c_str() ) );
  
  TH2D * h = new TH2D(name.c_str(), title.c_str(), nbins_q2, util::Q2BinEdges(s), nbins_ereco, util::ErecoBinEdges(s) );
  util::FormatErecoVsQ2Hist(h);
  for (int ir = 0 ; ir < nbins_ereco ; ir++ ) {
    for (int iq = 0 ; iq < nbins_q2 ; iq++ ) {
      h->SetBinContent(iq+1,ir+1,arr[ir][iq]);
    } 
  }
  return h;
}
//________________________________________________________________________________________
TH2D* SuperkExpectation::CurrPredictedErecoVsQ2Spectrum(
  valor::conv::SuperkSample_t s, valor::conv::SuperkMode_t m) const 
{ 
  const int nbins_q2    = util::NumQ2Bins(s);
  const int nbins_ereco = util::NumErecoBins(s);
  typedef double (erq2type)  [nbins_ereco][nbins_q2];  
  const erq2type & arr = *reinterpret_cast< const erq2type * > (CurrPredictedErecoVsQ2SpectrumArrayMode(s,m));  
  //const double * arr = CurrPredictedErecoSpectrumArrayMode(s,m);
  if (!arr) return 0;
  
  string name  (Form("curr_ereco_Vs_q2_%s_%s"                       ,util::AsStringCompact(s).c_str(), util::AsStringCompact(m).c_str()) );
  string title (Form("Current Predicted Ereco Vs Q2 Spectrum: %s %s",util::AsString(s).c_str()       , util::AsString(m).c_str()       ) );
  TH2D * h = new TH2D(name.c_str(), title.c_str(), nbins_q2, util::Q2BinEdges(s), nbins_ereco, util::ErecoBinEdges(s) );

  util::FormatErecoVsQ2Hist(h);
  util::SetDefaultSKSpectrumStyle(s,m,h);
  for (int ir = 0 ; ir < nbins_ereco ; ir++ ) {
    for (int iq = 0 ; iq < nbins_q2 ; iq++ ) {
      h->SetBinContent(iq+1,ir+1,arr[ir][iq]);
    }
  }

  return h;
}
//________________________________________________________________________________________
VALORArray * SuperkExpectation::CurrPredictedErecoVsQ2VArray(valor::conv::SuperkSample_t s) 
{
  VALORArray * arr = fPredictErecoVsQ2SpectrumArrayTot[s];
  if (!arr) {
    std::cerr<<"SKExp/CurrPredictedErecoVsQ2VArray could not find array for sample "<<util::AsString(s)<<std::endl;
    exit(-1);
  }
  
  return arr;
}
//________________________________________________________________________________________
VALORArray * SuperkExpectation::CurrPredictedErecoVsQ2VArrayMode(valor::conv::SuperkSample_t s) 
{
  VALORArray * arr = fPredictErecoVsQ2SpectrumArray[s];
  if (!arr) {
    std::cerr<<"SKExp/CurrPredictedErecoVsQ2VArrayMode could not find array for sample "<<util::AsString(s)<<std::endl;
    exit(-1);
  }
  
  return arr;
}
//________________________________________________________________________________________
const VALORArray * SuperkExpectation::CurrPredictedErecoVsQ2VArray(valor::conv::SuperkSample_t s) const
{
  const VALORArray * arr = fPredictErecoVsQ2SpectrumArrayTot[s];
  if (!arr) {
    std::cerr<<"SKExp/CurrPredictedErecoVsQ2VArray could not find array for sample "<<util::AsString(s)<<std::endl;
    exit(-1);
  }
  
  return arr;
}
//________________________________________________________________________________________
const VALORArray * SuperkExpectation::CurrPredictedErecoVsQ2VArrayMode(valor::conv::SuperkSample_t s) const
{
  const VALORArray * arr = fPredictErecoVsQ2SpectrumArray[s];
  if (!arr) {
    std::cerr<<"SKExp/CurrPredictedErecoVsQ2VArrayMode could not find array for sample "<<util::AsString(s)<<std::endl;
    exit(-1);
  }
  
  return arr;
}
//________________________________________________________________________________________
double* SuperkExpectation::CurrPredictedErecoVsQ2SpectrumArray(valor::conv::SuperkSample_t s) 
{
  VALORArray * arr = this->CurrPredictedErecoVsQ2VArray(s);
  return arr->GetArray();
}
//________________________________________________________________________________________
const double* SuperkExpectation::CurrPredictedErecoVsQ2SpectrumArray(valor::conv::SuperkSample_t s)  const
{
  const VALORArray * arr = this->CurrPredictedErecoVsQ2VArray(s);
  return arr->GetArray();
}
//________________________________________________________________________________________
double* SuperkExpectation::CurrPredictedErecoVsQ2SpectrumArrayMode(
  valor::conv::SuperkSample_t s, valor::conv::SuperkMode_t m)
{ 
  int sidx = util::ArrIdx::Instance()->SuperkSampleIdx(s);
  if(sidx == -1) {
    return 0;
  }
  int midx = util::ArrIdx::Instance()->SuperkModeIdx(s,m);
  if(midx == -1) {
    return 0;
  }
  
  VALORArray * arr = this->CurrPredictedErecoVsQ2VArrayMode(s);

  return &(arr->GetEntryRef(midx,0,0));
}
//________________________________________________________________________________________
const double* SuperkExpectation::CurrPredictedErecoVsQ2SpectrumArrayMode(
  valor::conv::SuperkSample_t s, valor::conv::SuperkMode_t m) const
{ 
  int sidx = util::ArrIdx::Instance()->SuperkSampleIdx(s);
  if(sidx == -1) {
    return 0;
  }
  int midx = util::ArrIdx::Instance()->SuperkModeIdx(s,m);
  if(midx == -1) {
    return 0;
  }
  
  const VALORArray * arr = this->CurrPredictedErecoVsQ2VArrayMode(s);

  return &(arr->GetEntryRef(midx,0,0));
}
//________________________________________________________________________________________
double SuperkExpectation::CurrPredictedErecoVsQ2SpectrumArray(valor::conv::SuperkSample_t s, int ir, int iq) const
{
  const int nbins_q2    = util::NumQ2Bins(s);
  const int nbins_ereco = util::NumErecoBins(s);
  typedef double (erq2type)  [nbins_ereco][nbins_q2];  
  const erq2type & arr = *reinterpret_cast< const erq2type * > (CurrPredictedErecoVsQ2SpectrumArray(s));  
  //const double * arr = CurrPredictedErecoVsQ2SpectrumArray(s);
  if (!arr) return 0;

  assert(arr);
  return arr[ir][iq];
} 
//

//________________________________________________________________________________________
void SuperkExpectation::PrintSystematicVariations(std::ostream & stream) const
{

}
//________________________________________________________________________________________
bool SuperkExpectation::Write(TObjArray * array, const char * pdf_tag, Option_t * write_opt)
{
  if(!array) {
    cout << "SKExp >> Null TObjArray. Not writing data" << endl;
    return false;
  }
  
  bool hist_mode = (write_opt && string(write_opt).find("hist")!=string::npos);
  //bool per_mode  = (write_opt && string(write_opt).find("mode")!=string::npos);
  
  string append_begin(Form("%s_", pdf_tag));
  
  for(int is=0; is<kns; is++) {
    SuperkSample_t s = kSKSamples[is];
    if(s == kSKsUndefined) continue;
    

#ifndef __USE_2D_TEMPLATES__ // 1D templates

    if (hist_mode) {
      util::SetDefaultStyle();
      TH1D * spectrum = 0;
      TH1D * spectrum_mode[knm] = { };
      
      // copy & format Ereco spectra
      spectrum = this->CurrPredictedErecoSpectrum(s);
      spectrum->SetDirectory (0);
      
      for(int im = 0; im < knm; im++) {
         SuperkMode_t m = kSKReactionModes[is][im];
         if(m == kSKmUndefined) continue;
         spectrum_mode[im] = this->CurrPredictedErecoSpectrum(s,m);
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
    } // hist mode
    else {
      VALORArray * total = new VALORArray(*(this->CurrPredictedErecoVArray(s)));
      total->SetName(Form("%s%s",append_begin.c_str(),total->GetName()));
      array->Add(total);
      VALORArray * bymode = new VALORArray(*(this->CurrPredictedErecoVArrayMode(s)));
      bymode->SetName(Form("%s%s",append_begin.c_str(),bymode->GetName()));
      array->Add(bymode);
    }
  }//is

#else // 2D templates
  
  // NEW 2D TEMPLATES
  
  if (hist_mode) {
      util::SetDefaultStyle();
      TH2D * spectrum = 0;
      TH2D * spectrum_mode[knm] = { };
      
      // copy & format Ereco spectra
      spectrum = this->CurrPredictedErecoVsQ2Spectrum(s);
      spectrum->SetDirectory (0);
      
      for(int im = 0; im < knm; im++) {
         SuperkMode_t m = kSKReactionModes[is][im];
         if(m == kSKmUndefined) continue;
         spectrum_mode[im] = this->CurrPredictedErecoVsQ2Spectrum(s,m);
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
    } // hist mode
    else {
      VALORArray * total = new VALORArray(*(this->CurrPredictedErecoVsQ2VArray(s)));
      total->SetName(Form("%s%s",append_begin.c_str(),total->GetName()));
      array->Add(total);
      VALORArray * bymode = new VALORArray(*(this->CurrPredictedErecoVsQ2VArrayMode(s)));
      bymode->SetName(Form("%s%s",append_begin.c_str(),bymode->GetName()));
      array->Add(bymode);
    }
  }//is

 //

#endif

  return true;
}
//________________________________________________________________________________________
bool SuperkExpectation::Write(
  const char * filename, const char * pdf_tag, Option_t * write_opt, Option_t * file_opt)
{
  if (!file_opt) file_opt = "update";
  TObjArray * array = new TObjArray();
  if (write_opt && string(write_opt).find("hist")!=string::npos) array->SetOwner(true);
  this->Write(array, pdf_tag, write_opt);
  TFile f(filename, file_opt);
  array->Write();
  f.Close();

  return true;
}
//________________________________________________________________________________________
void SuperkExpectation::ResetSystematics(Option_t * /*opt*/)
{
#ifdef __show_syst_reset_mesg__
  cout << "Exp >> Resetting all systematic variations." << endl;
#endif

  // SKDetFSI
  this->GetSuperkDetectorFSI()->Tweak(true,0);
  
  // BANFF
  this->GetBANFF()->Tweak(true,0);

  // X Section

  // NEW XSEC --> comment xsec syst that now are in banff or not used any longer
  // --> not commented any longer

  // NEW XSEC NIWG
#if defined(__USE_NIWG2013__)
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
  //RHC
  this->GetNuXSec()->Tweak(kNuisRHCToFHC, 0.);
#elif defined(__USE_NIWG2015__)
  this->GetNuXSec()->Tweak(kNuisXSecRPA , 1.); // NEW XSEC 
#endif

  //
  // RHC
  //

  // X Section

  // NEW XSEC --> comment xsec syst that now are in banff or not used any longer
  // --> not commented any longer

  // NEW XSEC NIWG
#if defined(__USE_NIWG2013__)
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
#elif defined(__USE_NIWG2015__)
  this->GetNuXSec()->Tweak(kRHCNuisXSecRPA           , 1.); // NEW XSEC
#endif

}
//________________________________________________________________________________________
void SuperkExpectation::TweakSystematic(string systematic, double twkdial, Option_t * opt)
{
  PhysParamId_t p = valor::util::PhysParamIdFromString(systematic);
  this->TweakSystematic(p,twkdial,opt);
}
//________________________________________________________________________________________
void SuperkExpectation::TweakSystematic(conv::PhysParamId_t p, double twkdial, Option_t * opt)
{
      
#ifdef __show_syst_tweak_mesg__
  cout << "Exp >> Setting systematic parameter: " << util::AsString(p) 
       << " to value = " << twkdial << endl;
#endif
  
  //
  // SKDET + FSI
  //

  // tweak an individual FSI parameter (ignoring correlations)
  if(util::IsSKDetFSI(p)) {
    this->GetSuperkDetectorFSI()->SetCurrentValue(p,twkdial);
  }//skdetfsi
  
  //
  // BANFF
  //
  
  // tweak an individual FSI parameter (ignoring correlations)
  else if(util::IsBANFF(p)) {
    this->GetBANFF()->SetCurrentValue(p,twkdial);
  }//banff


  //
  // XSec
  //

  // NEW XSEC --> attention to new xsec non-banff systematics. not changed 
  else if (util::IsXSec(p)) {
    this->GetNuXSec()->Tweak(p , twkdial);
  }
  
  else {
     cerr << "Exp/Tweak>> Unhandled systematic parameter (" << util::AsString(p) << ")" << endl;
     exit(1);
  }           
}


//________________________________________________________________________________________
void SuperkExpectation::RandomizeSystematics(Option_t * /*opt*/)
{
#ifdef __show_syst_tweak_mesg__
  cout << "Exp >> Generating random variations of all systematic params." << endl;
#endif

  // common to all the samples
  this->TweakAllUncorrelatedXSecSystematics(); 
  this->GetNuXSec()->SyncRHCToFHCSystematics(); // uncorr xsec are fully correlated rhc/fhc
  this->TweakAllSKDetectorFSISystematics();
  this->TweakAllBANFFSystematics();
  
}
//________________________________________________________________________________________
void SuperkExpectation::TweakAllBANFFSystematics(void)
{
#ifdef __show_syst_tweak_mesg__
    cout << "Exp >> Generating random variations of all BANFF systematic params (correlated)" << endl;
#endif
    // BANFF
    // Additon of util::IsBANFF as tweaking everything includes the energy scale
    this->GetBANFF()->Tweak(false,0,util::IsBANFF); // <-- Cholesky method
}
//________________________________________________________________________________________
void SuperkExpectation::TweakAllSKDetectorFSISystematics(void)
{
#ifdef __show_syst_tweak_mesg__
  cout << "Exp >> Generating random variations of all SKDET + FSI systematic params." << endl;
#endif
  // Energy scale is now included in the parameters
  // Uncorrelated Cholesky tweak should get Gaussian distributed random value
  this->GetSuperkDetectorFSI()->Tweak(false,0); // <-- Cholesky method
}
//________________________________________________________________________________________
void SuperkExpectation::TweakAllUncorrelatedXSecSystematics(void)
{
#ifdef __show_syst_tweak_mesg__
  cout << "Exp >> Generating random variations of all Uncorrelated XSec systematic params." << endl;
#endif

  // NEW XSEC --> comment all old non-banff xsec systematics --> not needed anylonger

  // NEW XSEC NIWG

#if defined(__USE_NIWG2013__)
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

#elif defined(__USE_NIWG2015__)
  // NEW XSEC
  // No longer tweak RPA
  //~ double tweak_discr = Rnd::Instance()->RndToyMCsyst().Uniform(0.,1.);
  //~ if(tweak_discr < 0.5) tweak_discr = 0.;
  //~ else tweak_discr = 1.;
  //~ this->GetNuXSec()->Tweak(kNuisXSecRPA             , tweak_discr);
#endif

}
//________________________________________________________________________________________
PhysParamList SuperkExpectation::Tweaks(Option_t * opt) const
{
  PhysParamList pl;
  this->CopyTweaksTo(pl, opt);
  return pl;
}
//________________________________________________________________________________________
void SuperkExpectation::CopyTweaksTo(PhysParamList & pl, Option_t * /*opt*/) const 
{
  pl.Clear();

  // 
  //  Uncorrelated xsec (common to FHC and RHC)
  // 

  
#if defined(__USE_NIWG2013__)
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

#elif defined(__USE_NIWG2015__)
  pl.SetValue( kNuisXSecRPA,             this->GetNuXSec()->CurrentValue(kNuisXSecRPA)            ); 
#endif

  fBANFF->CopyTweaksTo(pl);
  fSuperkDetectorFSI->CopyTweaksTo(pl);
  
  return;
}


//________________________________________________________________________________________
void SuperkExpectation::SetNuisanceParams(const FitParamList & fp, bool /*all_banff_fit*/, bool /*all_1Re_eff*/)
{
  //////////////////// SKDetFSI ///////////////////
  fSuperkDetectorFSI->UpdateFrom(fp); // EScale included
  
  fBANFF->UpdateFrom(fp);
    
  /////////// UNCORRELATED XSEC //////////
  
  // Cross-section systematics (other than BANFF ones)
  for(int ip = 0; ip < kNumOfXSecSyst; ip++) {
    PhysParamId_t p = kListOfXSecSyst[ip];
    if( fp.Exists(p) && util::IsXSec(p) ) {
      assert(!"CHECK THIS CODE PROPERLY BEFORE USING");
      double twk = fp.CurrValue(p);
      this->GetNuXSec()->Tweak(p,twk);
    }
  }
}
//________________________________________________________________________________________
void SuperkExpectation::PrintSystematics (std::ostream & stream, Option_t * /*opt*/) const
{
  cout << "Exp >> SuperK PDF tweaks: " << endl;
  
  fBANFF             -> PrintSystematicVariations(stream);
  
#ifdef __FHC_RUNNING__ // NEW FLUX
  fNuXSec            -> PrintSystematicVariations(stream);
#endif
#ifdef __RHC_RUNNING__
  fNuXSec            -> PrintRHCSystematicVariations(stream);
#endif

  fSuperkDetectorFSI -> PrintSystematicVariations(stream);
}
//________________________________________________________________________________________
double SuperkExpectation::GetFitStatisticPenalty(Option_t * /*opt*/) const
{
  bool all_syst=true;
  double chisq_penalty = 0.;
  //
  // Add penalty term for BANFF params 
  //
    
  double dchisq_penalty_banff = fBANFF->ChisqPenalty();

  chisq_penalty += dchisq_penalty_banff;
#ifdef __show_contributions_to_penalty_term__
  cout << "FitInputs::GetFitStatisticPenalty >> Current BANFF param variation"
       << " => \\delta\\chi^2_{penalty} = " << dchisq_penalty_banff << endl;
#endif

  //
  // Add penalty term for NuXSec params (other than the ones incl. in BANFF)
  //
  double dchisq_penalty_xsec =  fNuXSec->ChisqPenalty(all_syst);
  chisq_penalty += dchisq_penalty_xsec;
#ifdef __show_contributions_to_penalty_term__
  cout << "FitInputs::GetFitStatisticPenalty >> Current variation of params considered by NuXSec"
       << " => \\delta\\chi^2_{penalty} = " << dchisq_penalty_xsec << endl;
#endif
  
  //
  // Add penalty term for SuperkDetector + FSI params
  //
  
  double dchisq_penalty_skdetfsi = 0.;
  dchisq_penalty_skdetfsi = fSuperkDetectorFSI->ChisqPenalty(); 
  
  chisq_penalty += dchisq_penalty_skdetfsi;
#ifdef __show_contributions_to_penalty_term__
  cout << "FitInputs::GetFitStatisticPenalty >> Current SK+FSI param variation"
       << " => \\delta\\chi^2_{penalty} = " << dchisq_penalty_skdetfsi << endl;
#endif

  return chisq_penalty;
  }
//________________________________________________________________________________________
string SuperkExpectation::Name(void) const
{
  return "SuperkExpectation";
}
//________________________________________________________________________________________
// Should be removed ASAP
void SuperkExpectation::TweakAllRHCSKDetectorFSISystematics       (void)
{
  std::cerr<<"SKExp/TweakAllRHCSKDetectorFSISystematics deprecated - Use TweakAllSKDetectorFSISystematics instead"<<std::endl;
  this->TweakAllSKDetectorFSISystematics();
}
void SuperkExpectation::TweakAllFHCpRHCSKDetectorFSISystematics   (void)
{
  std::cerr<<"SKExp/TweakAllFHCpRHCSKDetectorFSISystematics deprecated - Use TweakAllSKDetectorFSISystematics instead"<<std::endl;
  this->TweakAllSKDetectorFSISystematics();
}
void SuperkExpectation::TweakAllMultiRingSKDetectorFSISystematics (void)
{
  std::cerr<<"SKExp/TweakAllMultiRingSKDetectorFSISystematics deprecated - Use TweakAllSKDetectorFSISystematics instead"<<std::endl;
  this->TweakAllSKDetectorFSISystematics();
}
void SuperkExpectation::TweakAllSampleSKDetectorFSISystematics    (void)
{
  std::cerr<<"SKExp/TweakAllSampleSKDetectorFSISystematics deprecated - Use TweakAllSKDetectorFSISystematics instead"<<std::endl;
  this->TweakAllSKDetectorFSISystematics();
}
void SuperkExpectation::TweakAllRHCBANFFSystematics       (void)
{
  std::cerr<<"SKExp/TweakAllRHCBANFFSystematics deprecated - Use TweakAllBANFFSystematics instead"<<std::endl;
  this->TweakAllBANFFSystematics();
}
void SuperkExpectation::TweakAllFHCpRHCBANFFSystematics   (void)
{
  std::cerr<<"SKExp/TweakAllFHCpRHCBANFFSystematics deprecated - Use TweakAllBANFFSystematics instead"<<std::endl;
  this->TweakAllBANFFSystematics();
}
void SuperkExpectation::TweakAllMultiRingBANFFSystematics (void)
{
  std::cerr<<"SKExp/TweakAllMultiRingBANFFSystematics deprecated - Use TweakAllBANFFSystematics instead"<<std::endl;
  this->TweakAllBANFFSystematics();
}
void SuperkExpectation::TweakAllSampleBANFFSystematics    (void)
{
  std::cerr<<"SKExp/TweakAllSampleBANFFSystematics deprecated - Use TweakAllBANFFSystematics instead"<<std::endl;
  this->TweakAllBANFFSystematics();
}
//________________________________________________________________________________________
