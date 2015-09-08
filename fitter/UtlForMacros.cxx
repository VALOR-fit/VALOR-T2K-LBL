//____________________________________________________________________________
/*

Copyright (c) 2010-2014, The VALOR group. All rights reserved.

ChangeLog:

*/
//____________________________________________________________________________


#include "UtlForMacros.h"

#include "SuperkExpectation.h"
#include "SuperkMC.h"
#include "BANFF.h"
#include "SuperkDetectorFSI.h"

#include "Conventions/Opt.h"
#include "Conventions/PhysParamId.h"
#include "Conventions/PhysHypothesisId.h"
#include "Conventions/Constants.h"
#include "FitCore/PhysHypothesis.h"
#include "Utilities/UtlXML.h"

#include "TSystem.h"
#include "TCanvas.h"

#include <iostream>
#include <string>

using std::cout;
using std::cerr;
using std::endl;
using std::string;

using namespace valor;
using namespace valor::conv;

//____________________________________________________________________________
const char * valor::util::FormFilePath(const char * file, const char * dir)
{
  FileStat_t dummy;
  const char * filename = Form("%s/%s/%s", gSystem->Getenv("T2KOSCANA"), dir, file);
  if(gSystem->GetPathInfo(filename, dummy)) {
    cerr << "File not found: " << filename << endl;
    exit(-1);
  }
  return filename; 
}
//____________________________________________________________________________
string valor::util::SysThrowIdToString(unsigned int tweak_systematics)
{
  string tweak_sys_str;
  switch(tweak_systematics)
    {
    case(0):
      tweak_sys_str = "None"; break;
    case(1):
      tweak_sys_str = "All"; break;
    case(2):
      tweak_sys_str = "Special low energy numu flux increase"; break;
    case(3):
      tweak_sys_str = "Tweak all but SK and FSI"; break;
    case(4):
      tweak_sys_str = "NC1pi norm = 0.5"; break;
    case(5):
      tweak_sys_str = "High energy >1.5 GeV flux increase"; break;
    case(6):
      tweak_sys_str = "BANFF"; break;
    case(7):
      tweak_sys_str = "XSec"; break;
    case(8):
      tweak_sys_str = "SK+FSI"; break;
    default:
      cerr << "SysThrowIdToString() tweak_systematics = " << tweak_systematics
	   << " not implemented" << endl;
      exit(-1);
    }
  return tweak_sys_str;
}
//____________________________________________________________________________
void valor::util::PerformSysThrow(unsigned int tweak_systematics, SuperkExpectation * skexp, bool print)
{
  if(tweak_systematics == 0) {
    //make sure systematics are reset
    skexp->ResetSystematicVariations();
    return;
  }
  else if(tweak_systematics == 1) {
    // tweak all sources of systematic err (gaussian err)
    cout << "- Randomizing all systematic params in pdf..." << endl;
    skexp->TweakAllSystematics();
    // skexp->TweakSystematic("sk_energy_scale",0.0,true);
  }
  else if(tweak_systematics == 2) {
    // special tweak, increasing numu flux parameters <700MeV by 5%
    cout << "- Increasing numu flux parameters <700MeV by 5% in pdf..." << endl;
    TVectorD twkv(skexp->GetBANFF()->BestFitValues());
    twkv[0] += 0.05;
    twkv[1] += 0.05;
    twkv[2] += 0.05;
    twkv[3] += 0.05;
    skexp->GetBANFF()->Tweak(true,&twkv);
  }
  else if(tweak_systematics == 3) {
    cout << "- Randomizing all but SK systematic params in pdf..." << endl;
    skexp->TweakAllBANFFSystematics();
    skexp->TweakAllUncorrelatedXSecSystematics();
  }
  else if(tweak_systematics == 4) {
    cout << "- Setting xsec_norm_nc1piC to -0.5" << endl;
    skexp->TweakSystematic("xsec_norm_nc1piC",-0.5);
  }
  else if(tweak_systematics == 5) {
    cout << "- Increasing flux parameters >1500MeV by 15% in pdf..." << endl;
    TVectorD twkv(skexp->GetBANFF()->BestFitValues());
    twkv[6]  += 0.15;
    twkv[7]  += 0.15;
    twkv[8]  += 0.15;
    twkv[9]  += 0.15;
    twkv[10] += 0.15;
    twkv[14] += 0.15;
    twkv[15] += 0.15;
    skexp->GetBANFF()->Tweak(true,&twkv);
  }
  else if(tweak_systematics == 6) {
    cout << "- Randomizing all BANFF systematic params in pdf..." << endl;
    skexp->TweakAllBANFFSystematics();
  }
  else if(tweak_systematics == 7) {
    cout << "- Randomizing all uncorrelated cross section systematic params in pdf..." << endl;
    skexp->TweakAllUncorrelatedXSecSystematics();
  }
  else if(tweak_systematics == 8) {
    cout << "- Randomizing all SK+FSI systematic params in pdf..." << endl;
    skexp->TweakAllSKDetectorFSISystematics();
  }
  else {
    cerr << "PerformSysThrow() tweak_systematics = " << tweak_systematics
	 << " not implemented" << endl;
    exit(-1);
  }
  if(print)
    skexp->PrintSystematicVariations();
}
//____________________________________________________________________________
void valor::util::SetMassHierarchy(SuperkExpectation * skexp)
{
  //  - iMH=0 --> kNmhNormal
  //  - iMH=1 --> kNmhInverted	
	
  // Don't marginalize Mass Hierarchy. Use the one already set
  if(! OscParams::Instance()->marginalize_mh) {
    OscParams::Instance()->MH = Opt::Instance()->gMassHierarchyAssumed;

    //convert to int
    if(OscParams::Instance()->MH == kNmhNormal)
      OscParams::Instance()->iMH = 0;
    else if(OscParams::Instance()->MH == kNmhInverted)
      OscParams::Instance()->iMH = 1; 
    else {
      cerr << "Mass Hierarchy is not undefined" << endl;
      exit(1);
    }
  }//!marginalize_mh
  // Extract randomly the Mass Hierarchy 
  else {
    OscParams::Instance()->iMH = util::RandomizeMH(OscParams::Instance()->mh_nh_prior);

    //convert to NuMassHierarchy_t
    if(OscParams::Instance()->iMH == 0)
      OscParams::Instance()->MH = kNmhNormal;
    else if(OscParams::Instance()->iMH == 1)
      OscParams::Instance()->MH = kNmhInverted;
    else {
      cerr << "util::SetMassHierarchy() Unknown value of iMH " << OscParams::Instance()->iMH << endl;
      exit(1);
    }
  }//marginalize_mh

  //set the hierarchy
  Opt::Instance()->gMassHierarchyAssumed = OscParams::Instance()->MH;
  skexp->GetSuperkMC()->SetMassHierarchy(OscParams::Instance()->MH);
  skexp->GetSuperkMC()->SetMassHierarchyBar(OscParams::Instance()->MH);

  ////save it
  //if(OscParams::Instance()->hToy_MH)
  //OscParams::Instance()->hToy_MH->Fill(OscParams::Instance()->iMH + conv::kASmallNum);
}
//____________________________________________________________________________
TH2D * valor::util::GetTH2DFromTFile(TFile * f, const char * name)
{
  TH2D * h = 0;
  f->GetObject(name, h);
  if(!h || h->IsZombie()) {
    cerr << "Error: the histogram " << name << " hasn't been found" << endl;
    exit(1);
  }
  return h;
}
//____________________________________________________________________________
//____________________________________________________________________________

