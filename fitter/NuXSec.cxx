//________________________________________________________________________________________
/*

 >> ChangeLog:

 For changes up to version 2010_v3.0.0 (last implementation of the 2010-era analysis)
 see code in the corresponding SVN branch.

 Changes since 2011_v0.9.0:

 @

*/
//________________________________________________________________________________________

#define __show_init_mesg_info__
//#define __show_xsec_ratio_mesg_debug__

#include <cassert>
#include <cstdlib>
#include <iostream>

#include <TH1D.h>
#include <TVirtualPad.h>
#include <TLegend.h>
#include <TMath.h>

#include "NuXSec.h"
#include "BANFF.h"
#include "Utilities.h"

using std::cout;
using std::cerr;
using std::endl;

using namespace t2kana;
using namespace t2kana::conv;
using namespace t2kana::util;

ClassImp(NuXSec)

//________________________________________________________________________________________
NuXSec::NuXSec()
{
  this->Init();
}
//________________________________________________________________________________________
NuXSec::~NuXSec()
{
  this->CleanUp();
}
//________________________________________________________________________________________
void NuXSec::Init(void)
{
  fBANFF              = 0;
}
//________________________________________________________________________________________
void NuXSec::CleanUp(void)
{

}
//________________________________________________________________________________________
bool NuXSec::Read(TXMLEngine* xml, XMLNodePointer_t root_node)
{
#ifdef __show_init_mesg_info__
  cout << "Exp/XSec >> Reading neutrino generator inputs from input XML file." << endl;
#endif


#ifdef __show_init_mesg_info__
  cout << "Exp/XSec >> Done reading neutrino cross-section inputs!" << endl;
#endif
  return true;  
}
//________________________________________________________________________________________
void NuXSec::TellMeAboutTheBANFFfit(BANFF * banff)
{
  fBANFF = banff;
}
//________________________________________________________________________________________
double NuXSec::SuperkXSecWght(
   double Etrue, t2kana::conv::SuperkSample_t s, t2kana::conv::SuperkMode_t m) const
{
  // a check to make sure the sample-mode combination exists
  if(s == kSKsUndefined) return 1.;
  if(m == kSKmUndefined) return 1.;
  if(util::FindModeSlot(s, m) < 0) return 1.;

  return 1.;
}
//________________________________________________________________________________________
void NuXSec::PrintSystematicVariations(void) const
{

}
//________________________________________________________________________________________
void NuXSec::Draw(TVirtualPad * c, Option_t * option) const
{
  if(!c) return;

  c->cd();

  TString opt = option;
  util::SetDefaultStyle();

  //
  // ...
  //
  
  c->Update();
}
//________________________________________________________________________________________
