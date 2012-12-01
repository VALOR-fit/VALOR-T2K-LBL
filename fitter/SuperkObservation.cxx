//________________________________________________________________________________________
/*

 >> ChangeLog:

 For changes up to version 2010_v3.0.0 (last implementation of the 2010-era analysis)
 see code in the corresponding SVN branch.

 Changes since 2011_v0.9.0:

 @

*/
//________________________________________________________________________________________

//#define __show_num_of_events_mesg_debug__

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <sstream>

#include <TString.h>
#include <TAxis.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TXMLEngine.h>
#include <TRandom3.h>
#include <TMath.h>

#include "SuperkObservation.h"
#include "Utilities.h"

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
double SuperkObservation::NumOfEvents(SuperkSample_t s) const
{
  double nev = 0;

  if(s == kSKsUndefined) {
    return nev;
  }

  int is = util::ArrIdx::Instance()->SuperkSampleIdx(s);
  TH1D * h = fErDistr[is];
  assert(h);
  nev = h->Integral();
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

  int is = util::ArrIdx::Instance()->SuperkSampleIdx(s);
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

  int is = util::ArrIdx::Instance()->SuperkSampleIdx(s);
  return fErDistr[is];
}
//________________________________________________________________________________________
bool SuperkObservation::Read(string config_file, int data_set_id)
{
  cout << "Obs >> " << "Reading config file: " << config_file << endl;
 
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

  ostringstream datanode1Rmu;
  datanode1Rmu << "superk/observation/one_mu_ring/ereco_values" << "/data_set_" << data_set_id;

  if(util::NodeExists(xml,root_node,datanode1Rmu.str())) {
     int idx = util::ArrIdx::Instance()->SuperkSampleIdx(kSKsFCFV1RingMuLike);
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
  }
  else {
     cerr << "Obs >> Node " << datanode1Rmu.str() << " does not exist!" << endl;
  }


  ostringstream datanode1Re;
  datanode1Re << "superk/observation/one_e_ring/event_count" << "/data_set_" << data_set_id;

  if(util::NodeExists(xml,root_node,datanode1Re.str())) {
     int idx = util::ArrIdx::Instance()->SuperkSampleIdx(kSKsFCFV1RingELike);
     fErDistr[idx] = new TH1D("", "", kNErecoBins1Re,kErecoBinEdges1Re);
     fErDistr[idx]->SetDirectory(0);
     double nev = util::XML2Dbl (xml,root_node,datanode1Re.str());
     fErDistr[idx]->SetBinContent(1,nev); ///// just put the number in the 1st bin
     cout << "Obs >> Read: " << fErDistr[idx]->Integral() 
       << " " << util::AsString(kSKsFCFV1RingELike) << " events." << endl;
  }
  else {
     cerr << "Obs >> Node " << datanode1Re.str() << " does not exist!" << endl;
  }

  // Release memory
  xml->FreeDoc(xmldoc);
  delete xml;

  return true;
}
//________________________________________________________________________________________
bool SuperkObservation::CreateMockData(
      SuperkExpectation & skexp, const FitParamList & fp, bool fluctuate)
{
// Get observed reconstructed energy distribution

  fp.ShowOnlyCurrValuesInNextPrint();
  cout << "Obs >> Generating a toy-MC expt with the following params: " << fp;
  cout << "Obs >> Stat. fluctuations are turned-" << ((fluctuate)? "ON" : "OFF") << endl;

  TRandom3 rand(0);
  bool applied_negevents_condition = skexp.Generate(fp); //generate prediction

  int idx = -1;
  idx = util::ArrIdx::Instance()->SuperkSampleIdx(kSKsFCFV1RingMuLike);
  if(idx != -1) {
    fErDistr[idx] = new TH1D("", "", kNErecoBins1Rmu,kErecoBinEdges1Rmu);
    fErDistr[idx] -> SetDirectory(0);
    TH1D* h1Rmu = skexp.CurrPredictedErecoSpectrum(kSKsFCFV1RingMuLike);
    for(int ibin = 1; ibin <= h1Rmu->GetNbinsX(); ibin++) {
       double nexp  = h1Rmu->GetBinContent(ibin);
       double nobs  = (fluctuate) ? rand.Poisson(nexp) : nexp;
       fErDistr[idx]->SetBinContent(ibin,nobs);
    }
  }
  idx = util::ArrIdx::Instance()->SuperkSampleIdx(kSKsFCFV1RingELike);
  if(idx != -1) {
    fErDistr[idx]  = new TH1D("", "", kNErecoBins1Re, kErecoBinEdges1Re);
    fErDistr[idx]  -> SetDirectory(0);
    TH1D* h1Re = skexp.CurrPredictedErecoSpectrum(kSKsFCFV1RingELike);
    for(int ibin = 1; ibin <= h1Re->GetNbinsX(); ibin++) {
       double nexp  = h1Re->GetBinContent(ibin);
       double nobs  = (fluctuate) ? rand.Poisson(nexp) : nexp;
       fErDistr[idx]->SetBinContent(ibin,nobs);
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

  TRandom3 rand(0);
  bool applied_negevents_condition = skexp.Generate(ph); //generate prediction

  int idx = -1;
  idx = util::ArrIdx::Instance()->SuperkSampleIdx(kSKsFCFV1RingMuLike);
  if(idx != -1) {
    fErDistr[idx] = new TH1D("", "", kNErecoBins1Rmu,kErecoBinEdges1Rmu);
    fErDistr[idx] -> SetDirectory(0);
    TH1D* h1Rmu = skexp.CurrPredictedErecoSpectrum(kSKsFCFV1RingMuLike);
    for(int ibin = 1; ibin <= h1Rmu->GetNbinsX(); ibin++) {
       double nexp  = h1Rmu->GetBinContent(ibin);
       double nobs  = (fluctuate) ? rand.Poisson(nexp) : nexp;
       fErDistr[idx]->SetBinContent(ibin,nobs);
    }
  }
  idx = util::ArrIdx::Instance()->SuperkSampleIdx(kSKsFCFV1RingELike);
  if(idx != -1) {
    fErDistr[idx]  = new TH1D("", "", kNErecoBins1Re, kErecoBinEdges1Re);
    fErDistr[idx]  -> SetDirectory(0);
    TH1D* h1Re = skexp.CurrPredictedErecoSpectrum(kSKsFCFV1RingELike);
    for(int ibin = 1; ibin <= h1Re->GetNbinsX(); ibin++) {
       double nexp  = h1Re->GetBinContent(ibin);
       double nobs  = (fluctuate) ? rand.Poisson(nexp) : nexp;
       fErDistr[idx]->SetBinContent(ibin,nobs);
    }
  }
  return applied_negevents_condition;
}
//________________________________________________________________________________________
void SuperkObservation::Draw(TVirtualPad * c, Option_t * option) const
{
  if(!c) return;

  c->cd();

  TString opt = option;
  util::SetDefaultStyle();

  if(opt.Contains("1RmuEreco")) { 
     int idx = util::ArrIdx::Instance()->SuperkSampleIdx(kSKsFCFV1RingMuLike);
     fErDistr[idx]->SetLineColor(kBlack);
     fErDistr[idx]->SetLineWidth(2);
     fErDistr[idx]->SetMarkerColor(kBlack);
     fErDistr[idx]->SetMarkerStyle(20);
     fErDistr[idx]->SetMarkerSize(1);
     fErDistr[idx]->Draw("E");
     fErDistr[idx]->GetXaxis()->SetTitle("E (GeV)");
     fErDistr[idx]->GetYaxis()->SetTitle("Number of 1Rmu events");
  }

  c->Update();
}
//________________________________________________________________________________________
void SuperkObservation::Init(void)
{
// Initialize priv data members
//
  for(int is=0; is<kns; is++) {
    fErDistr[is] = 0;
  }
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
