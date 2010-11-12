//____________________________________________________________________________
/*
Important revisions after version 1.0.0 :

*/
//____________________________________________________________________________

#include <cassert>
#include <cstdlib>
#include <iostream>

#include <TString.h>
#include <TAxis.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TXMLEngine.h>

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
}
//________________________________________________________________________________________
SuperkObservation::~SuperkObservation(void)
{
  this->CleanUp();
}
//________________________________________________________________________________________
double SuperkObservation::NumOfEvents(SuperkSample_t s) const
{
  if(!fErDistr[s]) return -1;

  if(s == t2kana::yr1numu::conv::kSKsUndefined) {
    return ( fErDistr[kSKsFCFV1RingMuLike]->Integral() +
             fErDistr[kSKsFCFVmRingMuLike]->Integral() );
  } else {
    return fErDistr[s]->Integral();
  }
}
//________________________________________________________________________________________
double SuperkObservation::NumOfEvents(double Er, SuperkSample_t s) const
{
  if(!fErDistr[s]) return -1;

  if(s == t2kana::yr1numu::conv::kSKsUndefined) {
    return ( fErDistr[kSKsFCFV1RingMuLike]->GetBinContent(fErDistr[s]->FindBin(Er)) +
             fErDistr[kSKsFCFVmRingMuLike]->GetBinContent(fErDistr[s]->FindBin(Er)) );
  } else {
    return fErDistr[s]->GetBinContent(fErDistr[s]->FindBin(Er));
  }
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

  // Get binning info
  fErBinning = util::XML2BinInfo(xml,root_node,"superk/binning/ereco");

  // Get observed recosntructed energy distribution
  fErDistr[kSKsFCFV1RingMuLike] = 
    util::XML2TH1D(xml,root_node,"superk/observation/one_mu_ring_tight/ereco");

  // Release memory
  xml->FreeDoc(xmldoc);
  delete xml;

  return true;
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
  fErBinning = 0;

  for(int is=0; is<kns; is++) {
    fErDistr[is] = 0;
  }
}
//________________________________________________________________________________________
void SuperkObservation::CleanUp(void)
{
// De-allocate all memory taken up by this object
//
  if(fErBinning) delete fErBinning;

  for(int is=0; is<kns; is++) {
    if(fErDistr[is]) delete fErDistr[is];
  }
}
//________________________________________________________________________________________
