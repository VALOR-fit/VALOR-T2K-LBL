//________________________________________________________________________________________
/*

 >> ChangeLog:

 For changes up to version 2010_v3.0.0 (last implementation of the 2010-era analysis)
 see code in the corresponding SVN branch.

 Changes since 2011_v0.9.0:

 @18/2/13 LE
 Read Ereco spectrum for 1Re.
 Doubled arrays (independent ones for 1Rmu and 1Re)

 31-12-2013  Nick
 Added samples for reverse horn current running

*/
//________________________________________________________________________________________

//#define __show_num_of_events_mesg_debug__

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <sstream>

#include <TFile.h>
#include <TObjArray.h>
#include <TXMLEngine.h>
#include <TString.h>
#include <TAxis.h>
#include <TH1D.h>
#include <TRandom3.h>
#include <TMath.h>

#include "SuperkObservation.h"
#include "Utilities.h"
#include "Rnd.h"

using std::ostringstream;
using std::cout;
using std::cerr;
using std::endl;

using namespace t2kana;
using namespace t2kana::conv;

ClassImp(SuperkObservation)

//________________________________________________________________________________________
SuperkObservation::SuperkObservation(void) 
{
  this->Init();
}
//________________________________________________________________________________________
SuperkObservation::~SuperkObservation(void)
{
  this->CleanUp();
}
//________________________________________________________________________________________
void SuperkObservation::Init(void)
{
// Initialize private data members
//
  for(int is=0; is<kns; is++) {
    SuperkSample_t s = kSKSamples[is];    
    fErDistr[is] = 0;
    for(int ir=0; ir<kNErecoBins[is]; ir++) {
      if (s== kSKsFCFV1RingMuLike) fErDistrArray1Rmu[ir] = 0;
      else fErDistrArray1Re[ir] = 0;
    }
    for(int ir=0; ir<kRHCNErecoBins[is]; ir++) {
      if (s== kSKsRHCFCFV1RingMuLike) fRHCErDistrArray1Rmu[ir] = 0;
      else fRHCErDistrArray1Re[ir] = 0;
    }
  }
  arr_index = util::ArrIdx::Instance();
}
//________________________________________________________________________________________
void SuperkObservation::CleanUp(void)
{
// De-allocate all memory taken up by this object
//
  for(int is=0; is<kns; is++) {
    if(fErDistr[is]) delete fErDistr[is];
  }
}
//________________________________________________________________________________________
bool SuperkObservation::Read(string config_file, int data_set_id, bool is_xml)
{
  cout << "Obs >> " << "Reading config file: " << config_file << endl;

  if(is_xml) { 
    // Parse the input configuration file
    TXMLEngine* xml = new TXMLEngine;
    XMLDocPointer_t xmldoc = xml->ParseFile(config_file.c_str());
    if (xmldoc==0) {
      cerr << "Obs >> " << "Can not parse input: " << config_file << endl;
      delete xml;
      exit(1);
    }
    XMLNodePointer_t root_node = xml->DocGetRootElement(xmldoc);

    //  if(util::NodeExists(xml,root_node,"superk/observation/one_mu_ring/ereco")) {
    //     // Get observed recosntructed energy distribution
    //     fErDistr[kSKsFCFV1RingMuLike] = 
    //        util::XML2TH1D(xml,root_node,"superk/observation/one_mu_ring/ereco");
    //  }
    //  else
#ifdef __FHC_RUNNING__
#ifndef __NUEONLY__ 
    ostringstream datanode1Rmu;
    datanode1Rmu << "superk/observation/one_mu_ring/ereco_values" << "/data_set_" << data_set_id;

    if(util::NodeExists(xml,root_node,datanode1Rmu.str())) {
      int idx = arr_index->SuperkSampleIdx(kSKsFCFV1RingMuLike);
      fErDistr[idx] = new TH1D("", "", kNErecoBins1Rmu,kErecoBinEdges1Rmu);
      fErDistr[idx]->SetDirectory(0);
      vector<double> values = util::XML2DblArray(xml,root_node,datanode1Rmu.str());
      vector<double>::const_iterator value_iter = values.begin();
      for( ; value_iter != values.end(); ++value_iter) {
	double Ereco = *value_iter;
	fErDistr[idx]->Fill(Ereco);
      }
      cout << "Obs >> Read: " << fErDistr[idx]->Integral() 
	   << " " << util::AsString(kSKsFCFV1RingMuLike) << " events." << endl;
      for(int ibin = 1; ibin <= kNErecoBins1Rmu; ibin++) {
	double nobs = fErDistr[idx]->GetBinContent(ibin);
	fErDistrArray1Rmu[ibin-1] = nobs;
      }
    }
    else {
      cerr << "Obs >> Node " << datanode1Rmu.str() << " does not exist!" << endl;
    }
#endif
#ifndef __NUMUONLY__ 
    ostringstream datanode1Re;
    //Ereco spectrum
    datanode1Re << "superk/observation/one_e_ring/ereco_values" << "/data_set_" << data_set_id;
    if(util::NodeExists(xml,root_node,datanode1Re.str())) {
      int idx = arr_index->SuperkSampleIdx(kSKsFCFV1RingELike);
      fErDistr[idx] = new TH1D("", "", kNErecoBins1Re,kErecoBinEdges1Re);
      fErDistr[idx]->SetDirectory(0);
      vector<double> values = util::XML2DblArray(xml,root_node,datanode1Re.str());
      vector<double>::const_iterator value_iter = values.begin();
      for( ; value_iter != values.end(); ++value_iter) {
	double Ereco = *value_iter;
	fErDistr[idx]->Fill(Ereco);
      }
      cout << "Obs >> Read: " << fErDistr[idx]->Integral() 
	   << " " << util::AsString(kSKsFCFV1RingELike) << " events." << endl;
      for(int ibin = 1; ibin <= kNErecoBins1Re; ibin++) {
	double nobs = fErDistr[idx]->GetBinContent(ibin);
	fErDistrArray1Re[ibin-1] = nobs;
      }
    }

    else {
      cerr << "Obs >> Node " << datanode1Re.str() << " does not exist!" << endl;
    }
#endif
#endif

#ifdef __RHC_RUNNING__
#ifndef __NUEONLY__
    ostringstream datanodeRHC1Rmu;
    datanodeRHC1Rmu << "superk/observation/one_rhc_mu_ring/ereco_values" << "/data_set_" << data_set_id;

    if(util::NodeExists(xml,root_node,datanodeRHC1Rmu.str())) {
      int idx = arr_index->SuperkSampleIdx(kSKsRHCFCFV1RingMuLike);
      fErDistr[idx] = new TH1D("", "", kRHCNErecoBins1Rmu,kRHCErecoBinEdges1Rmu);
      fErDistr[idx]->SetDirectory(0);
      vector<double> values = util::XML2DblArray(xml,root_node,datanodeRHC1Rmu.str());
      vector<double>::const_iterator value_iter = values.begin();
      for( ; value_iter != values.end(); ++value_iter) {
        double Ereco = *value_iter;
        fErDistr[idx]->Fill(Ereco);
      }
      cout << "Obs >> Read: " << fErDistr[idx]->Integral()
           << " " << util::AsString(kSKsRHCFCFV1RingMuLike) << " events." << endl;
      for(int ibin = 1; ibin <= kRHCNErecoBins1Rmu; ibin++) {
        double nobs = fErDistr[idx]->GetBinContent(ibin);
        fRHCErDistrArray1Rmu[ibin-1] = nobs;
      }
    }
    else {
      cerr << "Obs >> Node " << datanodeRHC1Rmu.str() << " does not exist!" << endl;
    }
#endif
#ifndef __NUMUONLY__
    ostringstream datanodeRHC1Re;
    //Ereco spectrum                                                                                                                                                                
    datanodeRHC1Re << "superk/observation/one_rhc_e_ring/ereco_values" << "/data_set_" << data_set_id;
    if(util::NodeExists(xml,root_node,datanodeRHC1Re.str())) {
      int idx = arr_index->SuperkSampleIdx(kSKsRHCFCFV1RingELike);
      fErDistr[idx] = new TH1D("", "", kRHCNErecoBins1Re,kRHCErecoBinEdges1Re);
      fErDistr[idx]->SetDirectory(0);
      vector<double> values = util::XML2DblArray(xml,root_node,datanodeRHC1Re.str());
      vector<double>::const_iterator value_iter = values.begin();
      for( ; value_iter != values.end(); ++value_iter) {
        double Ereco = *value_iter;
        fErDistr[idx]->Fill(Ereco);
      }
      cout << "Obs >> Read: " << fErDistr[idx]->Integral()
           << " " << util::AsString(kSKsRHCFCFV1RingELike) << " events." << endl;
      for(int ibin = 1; ibin <= kRHCNErecoBins1Re; ibin++) {
        double nobs = fErDistr[idx]->GetBinContent(ibin);
        fRHCErDistrArray1Re[ibin-1] = nobs;
      }
    }

    else {
      cerr << "Obs >> Node " << datanodeRHC1Re.str() << " does not exist!" << endl;
    }
#endif

#endif

    //Treat as counting experiment
    /*
      datanode1Re << "superk/observation/one_e_ring/event_count" << "/data_set_" << data_set_id;
      if(util::NodeExists(xml,root_node,datanode1Re.str())) {
      int idx = util::ArrIdx::Instance()->SuperkSampleIdx(kSKsFCFV1RingELike);
      fErDistr[idx] = new TH1D("", "", kNErecoBins1Re,kErecoBinEdges1Re);
      fErDistr[idx]->SetDirectory(0);
      double nev = util::XML2Dbl (xml,root_node,datanode1Re.str());
      fErDistr[idx]->SetBinContent(1,nev); ///// just put the number in the 1st bin
      cout << "Obs >> Read: " << fErDistr[idx]->Integral() 
      << " " << util::AsString(kSKsFCFV1RingELike) << " events." << endl;
      for(int ibin = 1; ibin <= kNErecoBins1Re; ibin++) {
      double nobs = fErDistr[idx]->GetBinContent(ibin);
      fErDistrArray1Re[ibin-1] = nobs;
      }
      }
    */

    // Release memory
    xml->FreeDoc(xmldoc);
    delete xml;
  }//is_xml

  else{
#ifndef __NUEONLY__ 
#ifdef __FHC_RUNNING__
    TFile f(config_file.c_str());
    TH1D * hist = 0;
    f.GetObject(Form("mock_data_%d",data_set_id),hist);
    if(!hist) {
      cerr << "Obs >> Histogram for 1Rmu does not exist!" << endl;
    }
    else {
      int idx = arr_index->SuperkSampleIdx(kSKsFCFV1RingMuLike);
      fErDistr[idx] = new TH1D("", "", kNErecoBins1Rmu,kErecoBinEdges1Rmu);
      fErDistr[idx]->SetDirectory(0);
      for(int ix = 1; ix <= hist->GetNbinsX(); ix++) {
	double Ereco = hist->GetBinCenter(ix);
	double nev   = hist->GetBinContent(ix);
	fErDistr[idx]->Fill(Ereco, nev);
      }
      cout << "Obs >> Read: " << fErDistr[idx]->Integral() 
	   << " " << util::AsString(kSKsFCFV1RingMuLike) << " events." << endl;
      for(int ibin = 1; ibin <= kNErecoBins1Rmu; ibin++) {
	double nobs = fErDistr[idx]->GetBinContent(ibin);
	fErDistrArray1Rmu[ibin-1] = nobs;
      }
    }
#endif
#ifdef __RHC_RUNNING__
    TFile fRHC(config_file.c_str());
    TH1D * histRHC = 0;
    fRHC.GetObject(Form("mock_data_%d",data_set_id),histRHC);
    if(!histRHC) {
      cerr << "Obs >> Histogram for RHC 1Rmu does not exist!" << endl;
    }
    else {
      int idx = arr_index->SuperkSampleIdx(kSKsRHCFCFV1RingMuLike);
      fErDistr[idx] = new TH1D("", "", kRHCNErecoBins1Rmu,kRHCErecoBinEdges1Rmu);
      fErDistr[idx]->SetDirectory(0);
      for(int ix = 1; ix <= hist->GetNbinsX(); ix++) {
        double Ereco = hist->GetBinCenter(ix);
        double nev   = hist->GetBinContent(ix);
        fErDistr[idx]->Fill(Ereco, nev);
      }
      cout << "Obs >> Read: " << fErDistr[idx]->Integral()
           << " " << util::AsString(kSKsRHCFCFV1RingMuLike) << " events." << endl;
      for(int ibin = 1; ibin <= kRHCNErecoBins1Rmu; ibin++) {
        double nobs = fErDistr[idx]->GetBinContent(ibin);
        fRHCErDistrArray1Rmu[ibin-1] = nobs;
      }
    }
#endif
#endif
  }//!is_xml



  return true;
}
//________________________________________________________________________________________
double SuperkObservation::NumOfEvents(SuperkSample_t s) const
{
  double nev = 0;

  if(s == kSKsUndefined) {
    return nev;
  }

  int is = arr_index->SuperkSampleIdx(s);
  for(int ir = 0; ir < kNErecoBins[is]; ir++){
    if (s== kSKsFCFV1RingMuLike) nev += fErDistrArray1Rmu[ir];
    else nev += fErDistrArray1Re[ir];
  }
  for(int ir = 0; ir < kRHCNErecoBins[is]; ir++){
    if (s== kSKsRHCFCFV1RingMuLike) nev += fRHCErDistrArray1Rmu[ir];
    else nev += fRHCErDistrArray1Re[ir];
  }
#ifdef __show_num_of_events_mesg_debug__
  cout << "Obs >> NEvents(" 
       << util::AsString(s) << ") = " << nev << endl;
#endif

  return nev;
}
//________________________________________________________________________________________
double SuperkObservation::NumOfEvents(double Er, SuperkSample_t s) const
{
  double nev = 0;

  if(s == kSKsUndefined) {
    return nev;
  }

  int is = arr_index->SuperkSampleIdx(s);
  TH1D * h = fErDistr[is];
  assert(h);
  int bin = h->FindBin(Er);
  nev = h->GetBinContent(bin);
#ifdef __show_num_of_events_mesg_debug__
  cout << "Obs >> NEvents(" 
       << util::AsString(s) << ", Ereco = " << Er << " GeV) = " << nev << endl;
#endif

  return nev;
}
//________________________________________________________________________________________
TH1D* SuperkObservation::ErecoSpectrum (t2kana::conv::SuperkSample_t s) const
{
  if(s == kSKsUndefined) {
    return 0;
  }

  int is = arr_index->SuperkSampleIdx(s);
  return fErDistr[is];
}
//________________________________________________________________________________________
double* SuperkObservation::ErecoSpectrumArray (t2kana::conv::SuperkSample_t s)
{
  if(s == kSKsUndefined) {
    return 0;
  }

  if (s== kSKsRHCFCFV1RingMuLike) return fRHCErDistrArray1Rmu;
  else if (s== kSKsRHCFCFV1RingELike) return fRHCErDistrArray1Re;
  else if (s== kSKsFCFV1RingELike) return fErDistrArray1Re; 
  else return fErDistrArray1Rmu;
}
//________________________________________________________________________________________
bool SuperkObservation::CreateMockData(
      SuperkExpectation & skexp, const FitParamList & fp, bool fluctuate)
{
// Get observed reconstructed energy distribution

  fp.ShowOnlyCurrValuesInNextPrint();
  cout << "Obs >> Generating a toy-MC expt with the following params: " << fp;
  cout << "Obs >> Stat. fluctuations are turned-" << ((fluctuate)? "ON" : "OFF") << endl;

  bool applied_negevents_condition = skexp.Generate(fp); //generate prediction

  int idx = -1;
  idx = arr_index->SuperkSampleIdx(kSKsFCFV1RingMuLike);
  
  if(idx != -1) {
    fErDistr[idx] = new TH1D("", "", kNErecoBins1Rmu,kErecoBinEdges1Rmu);
    fErDistr[idx] -> SetDirectory(0);
    //TH1D* h1Rmu = skexp.CurrPredictedErecoSpectrum(kSKsFCFV1RingMuLike);
    for(int ibin = 1; ibin <= kNErecoBins1Rmu; ibin++) {
       double nexp  = skexp.CurrPredictedErecoSpectrumArray(kSKsFCFV1RingMuLike)[ibin-1];
       double nobs  = (fluctuate) ? Rnd::Instance()->RndToyMCstat().Poisson(nexp) : nexp;
#ifndef PRODUCTIONMODE
       fErDistr[idx]->SetBinContent(ibin,nobs);
#endif
       fErDistrArray1Rmu[ibin-1] = nobs;
    }
  }
  idx = arr_index->SuperkSampleIdx(kSKsFCFV1RingELike);
  if(idx != -1) {
    fErDistr[idx]  = new TH1D("", "", kNErecoBins1Re, kErecoBinEdges1Re);
    fErDistr[idx]  -> SetDirectory(0);
    //    TH1D* h1Re = skexp.CurrPredictedErecoSpectrum(kSKsFCFV1RingELike);
    for(int ibin = 1; ibin <= kNErecoBins1Re; ibin++) {
      double nexp  = skexp.CurrPredictedErecoSpectrumArray(kSKsFCFV1RingELike)[ibin-1];
      double nobs  = (fluctuate) ? Rnd::Instance()->RndToyMCstat().Poisson(nexp) : nexp;
 #ifndef PRODUCTIONMODE
       fErDistr[idx]->SetBinContent(ibin,nobs);
#endif
       fErDistrArray1Re[ibin-1] = nobs;
    }
  }

  idx = arr_index->SuperkSampleIdx(kSKsRHCFCFV1RingMuLike);

  if(idx != -1) {
    fErDistr[idx] = new TH1D("", "", kRHCNErecoBins1Rmu,kRHCErecoBinEdges1Rmu);
    fErDistr[idx] -> SetDirectory(0);
    //TH1D* hRHC1Rmu = skexp.CurrPredictedErecoSpectrum(kSKsRHCFCFV1RingMuLike);
    for(int ibin = 1; ibin <= kRHCNErecoBins1Rmu; ibin++) {
      double nexp  = skexp.CurrPredictedErecoSpectrumArray(kSKsRHCFCFV1RingMuLike)[ibin-1];
      double nobs  = (fluctuate) ? Rnd::Instance()->RndToyMCstat().Poisson(nexp) : nexp;
#ifndef PRODUCTIONMODE
      fErDistr[idx]->SetBinContent(ibin,nobs);
#endif
      fRHCErDistrArray1Rmu[ibin-1] = nobs;
    }
  }
  idx = arr_index->SuperkSampleIdx(kSKsRHCFCFV1RingELike);
  if(idx != -1) {
    fErDistr[idx]  = new TH1D("", "", kRHCNErecoBins1Re, kRHCErecoBinEdges1Re);
    fErDistr[idx]  -> SetDirectory(0);
    //TH1D* hRHC1Re = skexp.CurrPredictedErecoSpectrum(kSKsRHCFCFV1RingELike);
    for(int ibin = 1; ibin <= kRHCNErecoBins1Re; ibin++) {
      double nexp  = skexp.CurrPredictedErecoSpectrumArray(kSKsRHCFCFV1RingELike)[ibin-1];
      double nobs  = (fluctuate) ? Rnd::Instance()->RndToyMCstat().Poisson(nexp) : nexp;
 #ifndef PRODUCTIONMODE
      fErDistr[idx]->SetBinContent(ibin,nobs);
#endif
      fRHCErDistrArray1Re[ibin-1] = nobs;
    }
  }

  return applied_negevents_condition;
}
//________________________________________________________________________________________
bool SuperkObservation::CreateMockData(
     SuperkExpectation & skexp, const PhysHypothesis & ph, bool fluctuate)
{
// Get observed reconstructed energy distribution

  cout << "Obs >> Generating a toy-MC expt with the following params: " << ph;
  cout << "Obs >> Stat. fluctuations are turned-" << ((fluctuate)? "ON" : "OFF") << endl;

  bool applied_negevents_condition = skexp.Generate(ph); //generate prediction

  int idx = -1;
  idx = arr_index->SuperkSampleIdx(kSKsFCFV1RingMuLike);
  if(idx != -1) {
    fErDistr[idx] = new TH1D("", "", kNErecoBins1Rmu,kErecoBinEdges1Rmu);
    fErDistr[idx] -> SetDirectory(0);
    //TH1D* h1Rmu = skexp.CurrPredictedErecoSpectrum(kSKsFCFV1RingMuLike);
    for(int ibin = 1; ibin <= kNErecoBins1Rmu; ibin++) {
       double nexp  = skexp.CurrPredictedErecoSpectrumArray(kSKsFCFV1RingMuLike)[ibin-1];
       double nobs  = (fluctuate) ? Rnd::Instance()->RndToyMCstat().Poisson(nexp) : nexp;
#ifndef PRODUCTIONMODE
       fErDistr[idx]->SetBinContent(ibin,nobs);
#endif
       fErDistrArray1Rmu[ibin-1] = nobs;
    }
  }
  idx = arr_index->SuperkSampleIdx(kSKsFCFV1RingELike);
  if(idx != -1) {
    fErDistr[idx]  = new TH1D("", "", kNErecoBins1Re, kErecoBinEdges1Re);
    fErDistr[idx]  -> SetDirectory(0);
    //    TH1D* h1Re = skexp.CurrPredictedErecoSpectrum(kSKsFCFV1RingELike);
    for(int ibin = 1; ibin <= kNErecoBins1Re; ibin++) {
      double nexp  = skexp.CurrPredictedErecoSpectrumArray(kSKsFCFV1RingELike)[ibin-1];
      double nobs  = (fluctuate) ? Rnd::Instance()->RndToyMCstat().Poisson(nexp) : nexp;
#ifndef PRODUCTIONMODE
       fErDistr[idx]->SetBinContent(ibin,nobs);
#endif
       fErDistrArray1Re[ibin-1] = nobs;
    }
  }

  idx = arr_index->SuperkSampleIdx(kSKsRHCFCFV1RingMuLike);
  if(idx != -1) {
    fErDistr[idx] = new TH1D("", "", kRHCNErecoBins1Rmu,kRHCErecoBinEdges1Rmu);
    fErDistr[idx] -> SetDirectory(0);
    //TH1D* hRHC1Rmu = skexp.CurrPredictedErecoSpectrum(kSKsRHCFCFV1RingMuLike);
    for(int ibin = 1; ibin <= kRHCNErecoBins1Rmu; ibin++) {
      double nexp  = skexp.CurrPredictedErecoSpectrumArray(kSKsRHCFCFV1RingMuLike)[ibin-1];
      double nobs  = (fluctuate) ? Rnd::Instance()->RndToyMCstat().Poisson(nexp) : nexp;
#ifndef PRODUCTIONMODE
      fErDistr[idx]->SetBinContent(ibin,nobs);
#endif
      fRHCErDistrArray1Rmu[ibin-1] = nobs;
    }
  }
  idx = arr_index->SuperkSampleIdx(kSKsRHCFCFV1RingELike);
  if(idx != -1) {
    fErDistr[idx]  = new TH1D("", "", kRHCNErecoBins1Re, kRHCErecoBinEdges1Re);
    fErDistr[idx]  -> SetDirectory(0);
    //TH1D* hRHC1Re = skexp.CurrPredictedErecoSpectrum(kSKsRHCFCFV1RingELike);
    for(int ibin = 1; ibin <= kRHCNErecoBins1Re; ibin++) {
      double nexp  = skexp.CurrPredictedErecoSpectrumArray(kSKsRHCFCFV1RingELike)[ibin-1];
      double nobs  = (fluctuate) ? Rnd::Instance()->RndToyMCstat().Poisson(nexp) : nexp;
#ifndef PRODUCTIONMODE
      fErDistr[idx]->SetBinContent(ibin,nobs);
#endif
      fRHCErDistrArray1Re[ibin-1] = nobs;
    }
  }

  return applied_negevents_condition;
}
//________________________________________________________________________________________
bool SuperkObservation::Write(TObjArray * data_collection, const char * data_tag)
{
  if(!data_collection) {
    cout << "Obs >> Null TObjArray. Not writing data" << endl;
    return false;
  }

  util::SetDefaultStyle();

  for(int is=0; is<kns; is++) {
    SuperkSample_t s = kSKSamples[is];
    if(s == kSKsUndefined) continue;

    TH1D * data = new TH1D(*fErDistr[is]);
    data -> SetDirectory(0);
    data -> SetMarkerStyle(20);
    data -> SetMarkerSize(1);
    data -> GetXaxis()->SetTitle("E_{reco} (GeV)");
    data -> GetYaxis()->SetTitle(Form("Number of %s events",util::AsString(s).c_str()));
    data -> SetName( 
       Form( "%s_%s", 
              data_tag,
              util::AsStringCompact(s).c_str() ) 
            );
#ifdef PRODUCTIONMODE
    for(int ir = 1; ir <= data->GetNbinsX(); ir++) {
      if(s == kSKsFCFV1RingMuLike)
	data->SetBinContent(ir, fErDistrArray1Rmu[ir-1]);
      else if(s == kSKsFCFV1RingELike)
	data->SetBinContent(ir, fErDistrArray1Re[ir-1]);
      if(s == kSKsRHCFCFV1RingMuLike)
        data->SetBinContent(ir, fRHCErDistrArray1Rmu[ir-1]);
      else if(s == kSKsRHCFCFV1RingELike)
        data->SetBinContent(ir, fRHCErDistrArray1Re[ir-1]);
    }//ir
#endif

    data_collection -> Add(data);
  }

  return true;
}
//________________________________________________________________________________________
bool SuperkObservation::Write(
  const char * filename, const char * data_tag, const char * option_file)
{
  TObjArray * array = new TObjArray();
  this->Write(array, data_tag);

  TFile f(filename, option_file);
  array->Write();
  f.Close();

  return true;
}
//________________________________________________________________________________________
