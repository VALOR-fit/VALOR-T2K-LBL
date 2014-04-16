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
#include <TCanvas.h>

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
    if (s==kSKsUndefined) continue;
    double * spectrum = ErecoSpectrumArray(s);
    for(int ir=0; ir<kNErecoBins[is]; ir++) {
      spectrum[ir] = 0;
    }
  }
  arr_index = util::ArrIdx::Instance();
}
//________________________________________________________________________________________
void SuperkObservation::CleanUp(void)
{
  
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
    
#ifdef __FHC_RUNNING__
#ifndef __NUEONLY__ 
    ostringstream datanode1Rmu;
    datanode1Rmu << "superk/observation/one_mu_ring/ereco_values" << "/data_set_" << data_set_id;

    if(util::NodeExists(xml,root_node,datanode1Rmu.str())) {
      TH1D h("", "", kNErecoBins1Rmu,kErecoBinEdges1Rmu);
      h.SetDirectory(0);
      vector<double> values = util::XML2DblArray(xml,root_node,datanode1Rmu.str());
      vector<double>::const_iterator value_iter = values.begin();
      for( ; value_iter != values.end(); ++value_iter) {
	double Ereco = *value_iter;
	h.Fill(Ereco);
      }
      cout << "Obs >> Read: " << h.Integral() 
	   << " " << util::AsString(kSKsFCFV1RingMuLike) << " events." << endl;
      for(int ibin = 1; ibin <= kNErecoBins1Rmu; ibin++) {
	double nobs = h.GetBinContent(ibin);
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
      TH1D h("", "", kNErecoBins1Re,kErecoBinEdges1Re);
      h.SetDirectory(0);
      vector<double> values = util::XML2DblArray(xml,root_node,datanode1Re.str());
      vector<double>::const_iterator value_iter = values.begin();
      for( ; value_iter != values.end(); ++value_iter) {
	double Ereco = *value_iter;
	h.Fill(Ereco);
      }
      cout << "Obs >> Read: " << h.Integral() 
	   << " " << util::AsString(kSKsFCFV1RingELike) << " events." << endl;
      for(int ibin = 1; ibin <= kNErecoBins1Re; ibin++) {
	double nobs = h.GetBinContent(ibin);
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
      TH1D h("", "", kRHCNErecoBins1Rmu,kRHCErecoBinEdges1Rmu);
      h.SetDirectory(0);
      vector<double> values = util::XML2DblArray(xml,root_node,datanodeRHC1Rmu.str());
      vector<double>::const_iterator value_iter = values.begin();
      for( ; value_iter != values.end(); ++value_iter) {
        double Ereco = *value_iter;
        h.Fill(Ereco);
      }
      cout << "Obs >> Read: " << h.Integral()
           << " " << util::AsString(kSKsRHCFCFV1RingMuLike) << " events." << endl;
      for(int ibin = 1; ibin <= kRHCNErecoBins1Rmu; ibin++) {
        double nobs = h.GetBinContent(ibin);
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
      TH1D h("", "", kRHCNErecoBins1Re,kRHCErecoBinEdges1Re);
      h.SetDirectory(0);
      vector<double> values = util::XML2DblArray(xml,root_node,datanodeRHC1Re.str());
      vector<double>::const_iterator value_iter = values.begin();
      for( ; value_iter != values.end(); ++value_iter) {
        double Ereco = *value_iter;
        h.Fill(Ereco);
      }
      cout << "Obs >> Read: " << h.Integral()
           << " " << util::AsString(kSKsRHCFCFV1RingELike) << " events." << endl;
      for(int ibin = 1; ibin <= kRHCNErecoBins1Re; ibin++) {
        double nobs = h.GetBinContent(ibin);
        fRHCErDistrArray1Re[ibin-1] = nobs;
      }
    }

    else {
      cerr << "Obs >> Node " << datanodeRHC1Re.str() << " does not exist!" << endl;
    }
#endif

#endif

    // Release memory
    xml->FreeDoc(xmldoc);
    delete xml;
  }//is_xml

  else { // DEPRECATE?
#ifndef __NUEONLY__ 
#ifdef __FHC_RUNNING__
    TFile f(config_file.c_str());
    TH1D * hist = 0;
    f.GetObject(Form("mock_data_%d",data_set_id),hist);
    if(!hist) {
      cerr << "Obs >> Histogram for 1Rmu does not exist!" << endl;
    }
    else {
      assert(hist->GetNbinsX() == util::NumErecoBins(kSKsFCFV1RingMuLike));
      double integral = 0;
      for(int ix = 1; ix <= hist->GetNbinsX(); ix++) {
        double nev   = hist->GetBinContent(ix);
        integral+=nev;
        fErDistrArray1Rmu[ix-1] = nev;
      }
      cout << "Obs >> Read: " << integral
           << " " << util::AsString(kSKsFCFV1RingMuLike) << " events." << endl;
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
      assert(histRHC->GetNbinsX() == util::NumErecoBins(kSKsRHCFCFV1RingMuLike));
      double integral = 0;
      int idx = arr_index->SuperkSampleIdx(kSKsRHCFCFV1RingMuLike);
      for(int ix = 1; ix <= histRHC->GetNbinsX(); ix++) {
        double Ereco = histRHC->GetBinCenter(ix);
        double nev   = histRHC->GetBinContent(ix);
        integral+=nev;
        fRHCErDistrArray1Rmu[ix-1] = nev;
      }
      cout << "Obs >> Read: " << integral
           << " " << util::AsString(kSKsRHCFCFV1RingMuLike) << " events." << endl;
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
  const double * spectrum = this->ErecoSpectrumArray(s);
  
  for(int ir = 0; ir < kNErecoBins[is]; ir++){
    nev += spectrum[ir];
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
    return 0;
  }
  
  cerr<<"WARNING: SuperkObservation::NumOfEvents(Ereco) needs testing"<<std::endl;
  // NEEDS TESTING
  const double * spectrum = this->ErecoSpectrumArray(s);
  
  int bin = util::ErecoBinContains(s,Er);
  assert(bin!=-1);
  nev = spectrum[util::ErecoBinContains(s,Er)];
  
#ifdef __show_num_of_events_mesg_debug__
  cout << "Obs >> NEvents(" 
       << util::AsString(s) << ", Ereco = " << Er << " GeV) = " << nev << endl;
#endif

  return nev;
}
//________________________________________________________________________________________
int SuperkObservation::NumOfBins() const
{
  int nbins = 0;
  for (int is = 0 ; is < kns ; is++) {
    SuperkSample_t s = kSKSamples[is];
    if (s==kSKsUndefined) continue;
    nbins+=util::NumErecoBins(s);
  }
  return nbins;
}
//________________________________________________________________________________________
TH1D* SuperkObservation::ErecoSpectrum (t2kana::conv::SuperkSample_t s) const
{
  if(s == kSKsUndefined) {
    return 0;
  }
  
  TH1D * data = new TH1D("","",util::NumErecoBins(s),util::ErecoBinEdges(s));
  util::FormatErecoHist(data);
  data -> SetDirectory(0);
  data -> SetMarkerStyle(20);
  data -> SetMarkerSize(1);
  data -> GetXaxis()->SetTitle("E_{reco} (GeV)");
  data -> GetYaxis()->SetTitle(Form("Number of %s events",util::AsString(s).c_str()));
  
  const double * spectrum = this->ErecoSpectrumArray(s);
  
  for(int ir = 1; ir <= data->GetNbinsX(); ir++) {
    data->SetBinContent(ir, spectrum[ir-1]);
  }//ir
  
  return data;
}
//________________________________________________________________________________________
double* SuperkObservation::ErecoSpectrumArray (t2kana::conv::SuperkSample_t s)
{
  /*if(s == kSKsUndefined) {
    return 0;
  } */ // This probably shouldn't fail silently

  if      (s == kSKsRHCFCFV1RingMuLike) return fRHCErDistrArray1Rmu;
  else if (s == kSKsRHCFCFV1RingELike)  return fRHCErDistrArray1Re;
  else if (s == kSKsFCFV1RingELike)     return fErDistrArray1Re; 
  else if (s == kSKsFCFV1RingMuLike)    return fErDistrArray1Rmu;
  
  std::cerr<<"SuperkObservation/ErecoSpectrumArray >> Failed for sample "<<util::AsString(s)<<std::endl;
  exit(-1);
}
//________________________________________________________________________________________
const double* SuperkObservation::ErecoSpectrumArray (t2kana::conv::SuperkSample_t s) const
{
  /*
  if(s == kSKsUndefined) {
    return 0;
  }*/ // This probably shouldn't fail silently

  if      (s == kSKsRHCFCFV1RingMuLike) return fRHCErDistrArray1Rmu;
  else if (s == kSKsRHCFCFV1RingELike)  return fRHCErDistrArray1Re;
  else if (s == kSKsFCFV1RingELike)     return fErDistrArray1Re; 
  else if (s == kSKsFCFV1RingMuLike)    return fErDistrArray1Rmu;
  
  std::cerr<<"SuperkObservation/ErecoSpectrumArray >> Failed for sample "<<util::AsString(s)<<std::endl;
  exit(-1);
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
  
  for (int is = 0 ; is < conv::kns ; is++ ) {
    SuperkSample_t s = kSKSamples[is];
    if (s==kSKsUndefined) continue;
    for(int ibin = 0; ibin < util::NumErecoBins(s); ibin++) {
       double nexp  = skexp.CurrPredictedErecoSpectrumArray(s)[ibin];
       double nobs  = (fluctuate) ? Rnd::Instance()->RndToyMCstat().Poisson(nexp) : nexp;
       this->ErecoSpectrumArray(s)[ibin] = nobs;
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
  
  for (int is = 0 ; is < conv::kns ; is++ ) {
    SuperkSample_t s = kSKSamples[is];
    if (s==kSKsUndefined) continue;
    for(int ibin = 0; ibin < util::NumErecoBins(s); ibin++) {
       double nexp  = skexp.CurrPredictedErecoSpectrumArray(s)[ibin];
       double nobs  = (fluctuate) ? Rnd::Instance()->RndToyMCstat().Poisson(nexp) : nexp;
       this->ErecoSpectrumArray(s)[ibin] = nobs;
    }
  }

  return applied_negevents_condition;
}
//________________________________________________________________________________________
bool SuperkObservation::Write(TObjArray * data_collection, const char * data_tag) const
{
  if(!data_collection) {
    cout << "Obs >> Null TObjArray. Not writing data" << endl;
    return false;
  }

  util::SetDefaultStyle();

  for(int is=0; is<kns; is++) {
    SuperkSample_t s = kSKSamples[is];
    if(s == kSKsUndefined) continue;
    
    // This is probably leaking memory.
    // Shouldn't be too much of an issue.
    TH1D * data = this->ErecoSpectrum(s);
    
    data -> SetName( 
     Form( "%s_%s", 
            data_tag,
            util::AsStringCompact(s).c_str() ) 
    );

    data_collection -> Add(data);
  }

  return true;
}
//________________________________________________________________________________________
bool SuperkObservation::Write(
  const char * filename, const char * data_tag, const char * option_file) const
{
  TObjArray * array = new TObjArray();
  this->Write(array, data_tag);

  TFile f(filename, option_file);
  array->Write();
  f.Close();

  return true;
}
//________________________________________________________________________________________
