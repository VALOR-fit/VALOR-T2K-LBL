//____________________________________________________________________________
/*
Important revisions after version 1.0.0 :

*/
//____________________________________________________________________________

#include <cassert>

#include <TString.h>
#include <TAxis.h>
#include <TH1D.h>
#include <TCanvas.h>

#include "year1_numu_osc_analysis/SuperkObservation.h"
#include "year1_numu_osc_analysis/Utilities.h"

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
  //
  //...
  //

  return true;
}
//________________________________________________________________________________________
void SuperkObservation::Draw(TVirtualPad * c, Option_t * option) const
{
  if(!c) return;

  c->cd();

  TString opt = option;
  //util::SetDefaultStyle();

  if(opt.Contains("Ereco")) { 
/*
     fErecoCCQElike->SetLineColor(kBlack);
     fErecoCCQElike->SetLineWidth(2);
     fErecoCCQElike->SetMarkerColor(kBlack);
     fErecoCCQElike->SetMarkerStyle(20);
     fErecoCCQElike->SetMarkerSize(1);
     fErecoCCQElike->Draw("E");
     fErecoCCQElike->GetXaxis()->SetTitle("E (GeV)");
     fErecoCCQElike->GetYaxis()->SetTitle("number of CCQE-like events");
*/
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
