//________________________________________________________________________________________
/*

Copyright (c) 2010-2014, The VALOR group. All rights reserved.

ChangeLog:

*/
//________________________________________________________________________________________

#include <sstream>
#include <TROOT.h>

#include "PhysHypothesisId.h"
#include "Opt.h"

using std::ostringstream;

using namespace t2kana;
using namespace t2kana::conv;

//________________________________________________________________________________________
string t2kana::util::AsString(t2kana::conv::PhysHypothesis_t h)
{
  switch(h) {
    case (kHypDecay) : 
       return "<Neutrino Decay>"; 
       break;
    case (kHypDecoh) : 
       return "<Neutrino Decoherence>"; 
       break;
    case (kHypNuMuDisapp2F) : 
       return "<2-Flavour Neutrino Oscillations: nu_mu disappearance>"; 
       break;
    case (kHypOsc3F) : 
       return "<3-flavour Neutrino Oscillations>"; 
       break;
    case (kHypNSI3F) : 
       return "<3-flavour Neutrino Oscillations with nonstandard matter interactions>"; 
       break;
    case (kHypOscSterile4F) : 
       return "<3+1 Sterile Neutrino Oscillations>"; 
       break;
    case (kHypOscSterile5F) : 
       return "<3+2 Sterile Neutrino Oscillations>"; 
       break;
    default:                
      return "undefined"; 
      break;
  }
  return "";
}
//________________________________________________________________________________________
string t2kana::util::AsString (
  PhysHypothesis_t hyp, double mix12, double mix13, double mix23, 
  double deltacp, double dm12sq, double dm23sq)
{
// The code builds strings which are typically used in building names of
// files containing t2kana results
//

  NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
  NuOscOctant_t   oct = Opt::Instance()->gOscOctantChoice;
  string mixparam = "?";
  switch(omp) {
      case(kNuMixSinsq2Theta):
         switch(oct) {
           case(kNuOctantOne):
             mixparam = "sinsq2theta_octantone";
             break;
           case(kNuOctantTwo):
             mixparam = "sinsq2theta_octanttwo";
             break;
           default:
             return "";
         }
         break;
      case(kNuMixSinsqTheta):
         mixparam = "sinsqtheta";
         break;
      default:
         return "";
  }

  ostringstream out;

  //
  // 3-flavour oscillations
  // Build the string as mix[12_nnn,13_nnn,23_nnn],dcp_nnn,dmsq[12_nnn,23_nnn]eV2
  // where mix is either sinsqtheta or sinsq2theta and nnn is the parameter
  // value with the specified numerical accuracy.
  // Params with value of -999 (randomized during toy-MC generation) will not be
  // included in the output string
  //
  if( hyp == t2kana::conv::kHypOsc3F )
  {
    double rndmz = -900;
    bool show_mix12   = (mix12   > rndmz);
    bool show_mix13   = (mix13   > rndmz);
    bool show_mix23   = (mix23   > rndmz);
    bool show_mix     = show_mix12 || show_mix13 || show_mix23;
    bool show_deltacp = (deltacp > rndmz);
    bool show_dm12sq  = (dm12sq  > rndmz);
    bool show_dm23sq  = (dm23sq  > rndmz);
    bool show_dmsq    = show_dm12sq || show_dm23sq;

    if(show_mix)
    {
       out << mixparam << "[";
       if(show_mix12) {
          out << "12_" << Form("%.3f",mix12);
       }
       if(show_mix13) {
          if(show_mix12) {
             out << ",";
          }
          out << "13_" << Form("%.3f",mix13);
       }
       if(show_mix23) {
          if(show_mix13 || (!show_mix13 && show_mix12)) {
             out << ",";
          }
          out << "23_" << Form("%.3f",mix23);
       }
       out << "]";
    }

    if(show_mix && (show_deltacp || show_dmsq))
    {
      out << ",";
    }

    if(show_deltacp)
    {
      out << "dcp_" << Form("%.3f",deltacp);
    }

    if(show_deltacp && show_dmsq)
    {
      out << ",";
    }

    if(show_dmsq)
    {
      out << "dmsq_" << util::AsString(Opt::Instance()->gMassHierarchyAssumed) << "[";
       if(show_dm12sq) {
          out << "12_" << Form("%.8f",dm12sq);
       }
       if(show_dm23sq) {
          if(show_dm12sq) {
             out << ",";
          }
          out << "23_" << Form("%.6f",dm23sq);
       }
       out << "]eV2";
    }

  }//3f

  // Placeholder for NSI
  if (hyp == t2kana::conv::kHypNSI3F) {
    out<<"nsi3f";
  }
  
  // Placeholder for steriles
  if (hyp == t2kana::conv::kHypOscSterile4F) {
    out<<"sterile_3+1";
  }
  if (hyp == t2kana::conv::kHypOscSterile5F) {
    out<<"sterile_3+2";
  }

  //
  // 2-flavour oscillations, numu disappearance
  // Build the string as mix23_nnn,dm23sq_nnneV2
  // where mix is either sinsqtheta or sinsq2theta and nnn is the parameter
  // value with the specified numerical accuracy.
  //
  if( hyp == t2kana::conv::kHypNuMuDisapp2F)
  {
    out << mixparam << "23_" << Form("%.3f",mix23) << ",dm23sq_" << Form("%.5f",dm23sq) << "eV2";
  }

  return out.str();
}
//________________________________________________________________________________________
PhysHypothesis_t t2kana::util::PhysHypFromString (string sh)
{
  if(sh == "decay"         ) return kHypDecay;
  if(sh == "decoh"         ) return kHypDecoh;
  if(sh == "2f_numu_disap" ) return kHypNuMuDisapp2F;
  if(sh == "3f_osc"        ) return kHypOsc3F;
  if(sh == "3f_nsi"        ) return kHypNSI3F;
  if(sh == "4f_sterile_osc") return kHypOscSterile4F;
  if(sh == "5f_sterile_osc") return kHypOscSterile5F;

  return kHypUndefined;
}
//________________________________________________________________________________________

