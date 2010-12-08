//____________________________________________________________________________
/*
Important revisions after version 1.0.0 :

*/
//____________________________________________________________________________

//#define __debug_level_1__
//#define __debug_level_2__
//#define __debug_level_3__

#include <cassert>
#include <cstdlib>
#include <iostream>

#include <TString.h>
#include <TAxis.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TXMLEngine.h>
#include <TRandom3.h>
#include <TMath.h>

#include "year1_numu_osc_analysis/SuperkObservation.h"
#include "year1_numu_osc_analysis/Utilities.h"

using std::cout;
using std::cerr;
using std::endl;

using namespace t2kana;
using namespace t2kana::yr1numu;
using namespace t2kana::yr1numu::conv;

ClassImp(SuperkObservation)

//________________________________________________________________________________________
SuperkObservation::SuperkObservation(void) 
{
  this->Init();

  fSamples.push_back(kSKsFCFV1RingMuLike);
}
//________________________________________________________________________________________
SuperkObservation::SuperkObservation(const vector<SuperkSample_t> & samples)
{
  this->Init();

  vector<SuperkSample_t>::const_iterator iter = samples.begin();
  for( ; iter != samples.end(); ++iter) {
    fSamples.push_back(*iter);
  }
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
    unsigned int nsamples = fSamples.size();
    for(unsigned int isample = 0; isample < nsamples; isample++) {
       SuperkSample_t sp = fSamples[isample];
       TH1D * h = fErDistr[sp];
       assert(h);
       double nevs = h->Integral();
#ifdef __debug_level_3__
       cout << "Obs >> NEvents(" 
          << util::AsString(sp) << ") = " << nevs << endl;
#endif
       nev += nevs;
    }//isample
#ifdef __debug_level_2__
    cout << "Obs >> NEvents(all_analysis_samples) = " << nev << endl;
#endif

  } else {

    TH1D * h = fErDistr[s];
    assert(h);
    nev = h->Integral();
#ifdef __debug_level_2__
    cout << "Obs >> NEvents(" 
         << util::AsString(s) << ") = " << nev << endl;
#endif

  }

  return nev;
}
//________________________________________________________________________________________
double SuperkObservation::NumOfEvents(double Er, SuperkSample_t s) const
{
  double nev = 0;

  if(s == kSKsUndefined) {
    unsigned int nsamples = fSamples.size();
    for(unsigned int isample = 0; isample < nsamples; isample++) {
       SuperkSample_t sp = fSamples[isample];
       TH1D * h = fErDistr[sp];
       assert(h);
       int bin = h->FindBin(Er);
       double nevs = h->GetBinContent(bin);
#ifdef __debug_level_3__
       cout << "Obs >> NEvents(" 
          << util::AsString(sp) << ", Ereco = " << Er << " GeV) = " << nevs << endl;
#endif
       nev += nevs;
    }//isample
#ifdef __debug_level_2__
    cout << "Obs >> NEvents(all_analysis_samples, Ereco = " 
         << Er << ") = " << nev << endl;
#endif

  } else {

    TH1D * h = fErDistr[s];
    assert(h);
    int bin = h->FindBin(Er);
    nev = h->GetBinContent(bin);
#ifdef __debug_level_2__
    cout << "Obs >> NEvents(" 
         << util::AsString(s) << ", Ereco = " << Er << " GeV) = " << nev << endl;
#endif

  }

  return nev;
}
//________________________________________________________________________________________
TH1D* SuperkObservation::AsTH1D(SuperkSample_t s) const
{
  TH1D * hEreco = new TH1D("", "", kNvErecoBins,kNvErecoBinEdges);
  hEreco->SetDirectory(0);

  for(int bin=1; bin <= hEreco->GetNbinsX(); bin++) {
    double Ereco = hEreco->GetBinCenter(bin);
    double nev   = this->NumOfEvents(Ereco,s);
    hEreco->SetBinContent(bin,nev);
  }

  return hEreco;
}
//________________________________________________________________________________________
const TAxis & SuperkObservation::ErecoBinning(void) const 
{ 
  return *(fErDistr[0]->GetXaxis()); 
}
//________________________________________________________________________________________
bool SuperkObservation::Read(string config_file)
{
  cout << "Obs >> " << "Reading config file: " << config_file << endl;
 
  // Parse the input configuration file
  TXMLEngine* xml = new TXMLEngine;
  XMLDocPointer_t xmldoc = xml->ParseFile(config_file.c_str());
  if (xmldoc==0) {
    cerr << "Obs >> " << "Can not parse input: " << config_file << endl;
    delete xml;
    return false;
  }
  XMLNodePointer_t root_node = xml->DocGetRootElement(xmldoc);

  if(util::NodeExists(xml,root_node,"superk/observation/one_mu_ring/ereco")) {
     // Get observed recosntructed energy distribution
     fErDistr[kSKsFCFV1RingMuLike] = 
        util::XML2TH1D(xml,root_node,"superk/observation/one_mu_ring/ereco");
  }
  else
  if(util::NodeExists(xml,root_node,"superk/observation/one_mu_ring/ereco_values")) {
     fErDistr[kSKsFCFV1RingMuLike] = new TH1D("", "", kNvErecoBins,kNvErecoBinEdges);
     fErDistr[kSKsFCFV1RingMuLike]->SetDirectory(0);
     vector<double> values = util::XML2DblArray(
        xml,root_node,"superk/observation/one_mu_ring/ereco_values");
     vector<double>::const_iterator value_iter = values.begin();
     for( ; value_iter != values.end(); ++value_iter) {
       double Ereco = *value_iter;
       fErDistr[kSKsFCFV1RingMuLike]->Fill(Ereco);
     }
  }

  // Release memory
  xml->FreeDoc(xmldoc);
  delete xml;

  return true;
}
//________________________________________________________________________________________
void SuperkObservation::CreateMockData(
      SuperkExpectation & exp, const PhysParamList & fp, bool fluctuate)
{
  // Get observed reconstructed energy distribution
  fErDistr[kSKsFCFV1RingMuLike] = new TH1D("", "", kNvErecoBins,kNvErecoBinEdges);
  fErDistr[kSKsFCFV1RingMuLike]->SetDirectory(0);

  TRandom3 rand(0);
  for(int ibin = 1; ibin <= fErDistr[kSKsFCFV1RingMuLike]->GetNbinsX(); ibin++) {
     double Ereco = fErDistr[kSKsFCFV1RingMuLike]->GetBinCenter(ibin);
     double nexp  = exp.Generate(Ereco, fp);
     double nobs  = (fluctuate) ? rand.Poisson(nexp) : nexp;
     fErDistr[kSKsFCFV1RingMuLike]->SetBinContent(ibin,nobs);
  }
}
//________________________________________________________________________________________
void SuperkObservation::Draw(TVirtualPad * c, Option_t * option) const
{
  if(!c) return;

  c->cd();

  TString opt = option;
  util::SetDefaultStyle();

  if(opt.Contains("Ereco")) { 
     if(opt.Contains("OneMuRingTight")) { 
        fErDistr[kSKsFCFV1RingMuLike]->SetLineColor(kBlack);
        fErDistr[kSKsFCFV1RingMuLike]->SetLineWidth(2);
        fErDistr[kSKsFCFV1RingMuLike]->SetMarkerColor(kBlack);
        fErDistr[kSKsFCFV1RingMuLike]->SetMarkerStyle(20);
        fErDistr[kSKsFCFV1RingMuLike]->SetMarkerSize(1);
        fErDistr[kSKsFCFV1RingMuLike]->Draw("E");
        fErDistr[kSKsFCFV1RingMuLike]->GetXaxis()->SetTitle("E (GeV)");
        fErDistr[kSKsFCFV1RingMuLike]->GetYaxis()->SetTitle("Number 1RingMuTight events");
     }
     else {
        exit(1);
     }
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

  fSamples.clear();
}
//________________________________________________________________________________________
