//____________________________________________________________________________
/*

\copyright (c) 2010-2013, VaLOR group.

30-12-2013  Nick
Added samples and modes for reverse horn current running

*/
//____________________________________________________________________________

//#define __show_cholesky_decomposition_mesg_debug__
//#define __show_cholesky_gen_corr_param_variations_mesg_debug__
//#define __show_xml_mesg_info__
//#define __show_xml_mesg_debug__
//#define __show_arr_idx_mesg_debug__

#include <cassert>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>

#include <TROOT.h>
#include <TMath.h>
#include <TStyle.h>
#include <TColor.h>
#include <TGraph.h>
#include <TH1.h>
#include <TH2.h>
#include <TAxis.h>
#include <TVectorD.h>
#include <TMatrixD.h>
#include <TRandom3.h>
#include <TXMLEngine.h>
#include <TFile.h>
#include <TSystem.h>

#include "Utilities.h"
#include "Opt.h"
#include "Rnd.h"

using std::cout;
using std::cerr;
using std::endl;
using std::setw;
using std::setfill;
using std::fixed;
using std::scientific;
using std::setprecision;
using std::ofstream;
using std::ostringstream;

using namespace t2kana;
using namespace t2kana::conv;

//____________________________________________________________________________
string t2kana::util::AsString (
  t2kana::conv::PhysHypothesis_t hyp, 
  double mix12, double mix13, double mix23, 
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
//____________________________________________________________________________
string t2kana::util::AsString(t2kana::conv::OscAnaType_t a)
{
  switch(a) {
   case (kAnaCombined) :
     return "Combined nue appearance + numu disappearance analysis";
     break;
   case (kAnaNuMuDisappOnly) :
     return "Numu disappearance analysis only";
     break;
   case (kAnaNuEAppOnly) :
     return "Nue appearance analysis only";
     break;
   default :
     return "Undefined"; 
     break;
  }
  return "";
}
//____________________________________________________________________________
string t2kana::util::AsString(t2kana::conv::NuMassHierarchy_t m)
{
  switch(m) {
   case (kNmhNormal) :
     return "Normal";
     break;
   case (kNmhInverted) :
     return "Inverted";
     break;
   default :
     return "Undefined"; 
     break;
  }
  return "";
}
//____________________________________________________________________________
string t2kana::util::AsString(t2kana::conv::NuOscMixParam_t p)
{
  switch(p) {   
   case (kNuMixSinsq2Theta) :
     return "sin^2(2*theta_{ij})";
     break;
   case (kNuMixSinsqTheta) :
     return "sin^2(theta_{ij})";
     break;
   default :
     return "Undefined"; 
     break;
  }
  return "";
}
//____________________________________________________________________________
string t2kana::util::AsString(t2kana::conv::NuOscOctant_t o)
{
  switch(o) {   
   case (kNuOctantOne) :
     return "sin^2(theta_{23}) <= 0.5";
     break;
   case (kNuOctantTwo) :
     return "sin^2(theta_{23}) >= 0.5";
     break;
   default :
     return "Undefined"; 
     break;
  }
  return "";
}
//____________________________________________________________________________
string t2kana::util::AsString(t2kana::conv::NuOscDmsqScheme_t dmsqs)
{
  switch(dmsqs) {
  case(kNuDmsq12andAbs23) :
    return "dm^2 = Dm^2_{21}, Dm^2 = |Dm^2_{32}|";
    break;
  case(kNuDmsq12andAbs31) :
    return "dm^2 = Dm^2_{21}, Dm^2 = |Dm^2_{31}|";
    break;
  case(kNuDmsqFogliLisi ) :
    return "dm^2 = Dm^2_{21}, Dm^2 = |m^2_{3} - (m^2_{1} + m^2_{1})/2|";
    break;
  default :
    return "Undefined";
    break;
  }
}
//____________________________________________________________________________
string t2kana::util::AsString(t2kana::conv::NuOscDeltaCP_t deltaCP)
{
  switch(deltaCP) {
  case(kNuDeltaCPPhase) :
    return "d_{CP}";
    break;
  case(kNuDeltaCPSin) :
    return "sin(d_{CP})";
    break;
  default :
    return "Undefined";
    break;
  }
}
//____________________________________________________________________________
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
//___________________________________________________________________________
//____________________________________________________________________________
string t2kana::util::AsString(conv::SKMCTemplateType_t t)
{
  switch(t) {
    case ( kSKMCtNom ) : 
      return "nominal";          
      break;
    case ( kSKMCtTun ) : 
      return "tuned"; 
      break;
    case ( kSKMCtSys ) : 
      return "tuned + systematic variation";    
      break;
    case ( kSKMCtFin ) : 
      return "tuned + systematic variation + physics hypothesis";          
      break;
    default : 
      return "--";               
      break;
  }
  return "";
}
//____________________________________________________________________________
string t2kana::util::AsString(t2kana::conv::FluxType_t f)
{
  switch(f) {
    case ( kFlxNominal        ) : return "nominal";          break;
    case ( kFlxExtTuning      ) : return "tuned (external)"; break;
    case ( kFlxFinalTuning    ) : return "tuned (final)";    break;
    case ( kFlxCurrent        ) : return "current";          break;
    case ( kRHCFlxNominal     ) : return "nominal; RHC";          break;
    case ( kRHCFlxExtTuning   ) : return "tuned (external); RHC"; break;
    case ( kRHCFlxFinalTuning ) : return "tuned (final); RHC";    break;
    case ( kRHCFlxCurrent     ) : return "current; RHC";          break;
    default                     : return "--";               break;
  }
  return "";
}
//____________________________________________________________________________
string t2kana::util::AsString(t2kana::conv::OscChannelId_t osch)
{
  switch(osch) {
     case ( kOscNuE2NuE           ) : return "nue -> nue";            break;
     case ( kOscNuE2NuMu          ) : return "nue -> numu";           break;
     case ( kOscNuE2NuTau         ) : return "nue -> nutau";          break;
     case ( kOscNuMu2NuE          ) : return "numu -> nue";           break;  
     case ( kOscNuMu2NuMu         ) : return "numu -> numu";          break;
     case ( kOscNuMu2NuTau        ) : return "numu -> nutau";         break;
     case ( kOscNuTau2NuE         ) : return "nutau -> nue";          break;
     case ( kOscNuTau2NuMu        ) : return "nutau -> numu";         break;
     case ( kOscNuTau2NuTau       ) : return "nutau -> nutau";        break;
     case ( kOscNuEBar2NuEBar     ) : return "nuebar -> nuebar";      break;
     case ( kOscNuEBar2NuMuBar    ) : return "nuebar -> numubar";     break;
     case ( kOscNuEBar2NuTauBar   ) : return "nuebar -> nutaubar";    break;
     case ( kOscNuMuBar2NuEBar    ) : return "numubar -> nuebar";     break;
     case ( kOscNuMuBar2NuMuBar   ) : return "numubar -> numubar";    break;
     case ( kOscNuMuBar2NuTauBar  ) : return "numubar -> nutaubar";   break;
     case ( kOscNuTauBar2NuEBar   ) : return "nutaubar -> nuebar";    break;
     case ( kOscNuTauBar2NuMuBar  ) : return "nutaubar -> numubar";   break;
     case ( kOscNuTauBar2NuTauBar ) : return "nutaubar -> nutaubar";  break;
     default                        : return "--";
  }
  return "--";
}
//____________________________________________________________________________
string t2kana::util::AsStringRooTex (
  t2kana::conv::SuperkSample_t s, t2kana::conv::SuperkMode_t m, 
  string append_begin, string append_end)
{
  ostringstream out;

  if(append_begin.size() > 0) {
    out << append_begin;
  }

  out << t2kana::util::AsStringRooTex(s);

  if(m != kSKmUndefined) {
    out << " (";
    out << t2kana::util::AsStringRooTex(m);
    out << ")";
  }

  if(append_end.size() > 0) {
    out << append_end;
  }

  return out.str();
}
//____________________________________________________________________________
string t2kana::util::AsStringCompact (
  t2kana::conv::SuperkSample_t s, t2kana::conv::SuperkMode_t m, 
  string append_begin, string append_end)
{
// Compact string used for naming output histograms stored ROOT files
//
  ostringstream out;

  if(append_begin.size() > 0) {
    out << append_begin;
  }

  out << t2kana::util::AsStringCompact(s);

  if(m != kSKmUndefined) {
    out << "_";
    out << t2kana::util::AsStringCompact(m);
  }

  if(append_end.size() > 0) {
    out << append_end;
  }

  return out.str();
}
//____________________________________________________________________________
string t2kana::util::AsStringCompact_splinehack (t2kana::conv::SuperkSample_t s)
{
// Compact string used for naming output histograms stored ROOT files
//
  switch(s) {
     case ( kSKsRHCFCFV1RingMuLike ) : return "1Rmu"; break;
     case ( kSKsRHCFCFV1RingELike  ) : return "1Re";  break;

     default                         : return AsStringCompact(s);
  }
  return AsStringCompact(s);
}
//___________________________________________________________________________
t2kana::conv::PhysHypothesis_t t2kana::util::PhysHypFromString(string sh)
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
//____________________________________________________________________________
t2kana::conv::OscAnaType_t  
   t2kana::util::OscAnaFromString(string sh) 
{
  if(sh == "combined"   ) return kAnaCombined;
  if(sh == "numu_disap" ) return kAnaNuMuDisappOnly;
  if(sh == "nue_app"    ) return kAnaNuEAppOnly;
  
  return kAnaUndefined;
}
//____________________________________________________________________________
t2kana::conv::NuOscDmsqScheme_t  
  t2kana::util::DmsqSchemeFromString(string sh) 
{
  if(sh == "abs23"    ) return kNuDmsq12andAbs23;
  if(sh == "abs31"    ) return kNuDmsq12andAbs31;
  if(sh == "fogli"    ) return kNuDmsqFogliLisi;
  
  return kNuDmsqUndefined;
}
//____________________________________________________________________________
t2kana::conv::NuMassHierarchy_t  
  t2kana::util::HierarchyFromString(string sh) 
{
  if(sh == "normal"   ) return kNmhNormal;
  if(sh == "inverted" ) return kNmhInverted;
  
  return kNmhUndefined;
}
//____________________________________________________________________________
t2kana::conv::NuOscOctant_t   
  t2kana::util::OctantFromString(string sh) 
{
  if(sh == "one") return kNuOctantOne;
  if(sh == "two") return kNuOctantTwo;
  
  return kNuOctantUndefined;
}
//____________________________________________________________________________
t2kana::conv::NuOscMixParam_t   
  t2kana::util::MixingTypeFromString(string sh) 
{
  if(sh == "sinsq2theta") return kNuMixSinsq2Theta;
  if(sh == "sinsqtheta" ) return kNuMixSinsqTheta;
  
  return kNuMixUndefined;
}
//____________________________________________________________________________
t2kana::conv::SKMCTemplateType_t   
  t2kana::util::SKMCTemplateTypeFromString(string sh) 
{
  if(sh == "ori")        return kSKMCtOri;
  if(sh == "nom")        return kSKMCtNom;
  if(sh == "tun")        return kSKMCtTun;
  if(sh == "sys")        return kSKMCtSys;
  if(sh == "fin")        return kSKMCtFin;

  return kSKMCtUndefined;
}
//____________________________________________________________________________
t2kana::conv::OscChannelId_t
  t2kana::util::GetOscChannel(int nu_start, int nu_end)
{
  if ( nu_start == kPdgNuE      && nu_end == kPdgNuE      ) return kOscNuE2NuE;
  if ( nu_start == kPdgNuE      && nu_end == kPdgNuMu     ) return kOscNuE2NuMu;
  if ( nu_start == kPdgNuE      && nu_end == kPdgNuTau    ) return kOscNuE2NuTau;
  if ( nu_start == kPdgNuMu     && nu_end == kPdgNuE      ) return kOscNuMu2NuE;
  if ( nu_start == kPdgNuMu     && nu_end == kPdgNuMu     ) return kOscNuMu2NuMu;
  if ( nu_start == kPdgNuMu     && nu_end == kPdgNuTau    ) return kOscNuMu2NuTau;
  if ( nu_start == kPdgNuTau    && nu_end == kPdgNuE      ) return kOscNuTau2NuE;
  if ( nu_start == kPdgNuTau    && nu_end == kPdgNuMu     ) return kOscNuTau2NuMu;
  if ( nu_start == kPdgNuTau    && nu_end == kPdgNuTau    ) return kOscNuTau2NuTau;
     
  if ( nu_start == kPdgNuEBar   && nu_end == kPdgNuEBar   ) return kOscNuEBar2NuEBar;
  if ( nu_start == kPdgNuEBar   && nu_end == kPdgNuMuBar  ) return kOscNuEBar2NuMuBar;
  if ( nu_start == kPdgNuEBar   && nu_end == kPdgNuTauBar ) return kOscNuEBar2NuTauBar;
  if ( nu_start == kPdgNuMuBar  && nu_end == kPdgNuEBar   ) return kOscNuMuBar2NuEBar;
  if ( nu_start == kPdgNuMuBar  && nu_end == kPdgNuMuBar  ) return kOscNuMuBar2NuMuBar;
  if ( nu_start == kPdgNuMuBar  && nu_end == kPdgNuTauBar ) return kOscNuMuBar2NuTauBar;
  if ( nu_start == kPdgNuTauBar && nu_end == kPdgNuEBar   ) return kOscNuTauBar2NuEBar;
  if ( nu_start == kPdgNuTauBar && nu_end == kPdgNuMuBar  ) return kOscNuTauBar2NuMuBar;
  if ( nu_start == kPdgNuTauBar && nu_end == kPdgNuTauBar ) return kOscNuTauBar2NuTauBar;

  //oscillations to or from sterile neutrinos NuS4 or NuS5
  if ( nu_start == kPdgNuE      && nu_end == kPdgNuS4     ) return kOscNuE2NuS4;
  if ( nu_start == kPdgNuE      && nu_end == kPdgNuS5     ) return kOscNuE2NuS5;
  if ( nu_start == kPdgNuMu     && nu_end == kPdgNuS4     ) return kOscNuMu2NuS4;
  if ( nu_start == kPdgNuMu     && nu_end == kPdgNuS5     ) return kOscNuMu2NuS5;
  if ( nu_start == kPdgNuTau    && nu_end == kPdgNuS4     ) return kOscNuTau2NuS4;
  if ( nu_start == kPdgNuTau    && nu_end == kPdgNuS5     ) return kOscNuTau2NuS5;

  if ( nu_start == kPdgNuS4     && nu_end == kPdgNuE      ) return kOscNuS42NuE;
  if ( nu_start == kPdgNuS5     && nu_end == kPdgNuE      ) return kOscNuS52NuE;
  if ( nu_start == kPdgNuS4     && nu_end == kPdgNuMu     ) return kOscNuS42NuMu;
  if ( nu_start == kPdgNuS5     && nu_end == kPdgNuMu     ) return kOscNuS52NuMu;
  if ( nu_start == kPdgNuS4     && nu_end == kPdgNuTau    ) return kOscNuS42NuTau;
  if ( nu_start == kPdgNuS5     && nu_end == kPdgNuTau    ) return kOscNuS52NuTau;

  if ( nu_start == kPdgNuEBar   && nu_end == kPdgNuS4Bar  ) return kOscNuEBar2NuS4Bar;
  if ( nu_start == kPdgNuEBar   && nu_end == kPdgNuS5Bar  ) return kOscNuEBar2NuS5Bar;
  if ( nu_start == kPdgNuMuBar  && nu_end == kPdgNuS4Bar  ) return kOscNuMuBar2NuS4Bar;
  if ( nu_start == kPdgNuMuBar  && nu_end == kPdgNuS5Bar  ) return kOscNuMuBar2NuS5Bar;
  if ( nu_start == kPdgNuTauBar && nu_end == kPdgNuS4Bar  ) return kOscNuTauBar2NuS4Bar;
  if ( nu_start == kPdgNuTauBar && nu_end == kPdgNuS5Bar  ) return kOscNuTauBar2NuS5Bar;

  if ( nu_start == kPdgNuS4Bar  && nu_end == kPdgNuEBar   ) return kOscNuS4Bar2NuEBar;
  if ( nu_start == kPdgNuS5Bar  && nu_end == kPdgNuEBar   ) return kOscNuS5Bar2NuEBar;
  if ( nu_start == kPdgNuS4Bar  && nu_end == kPdgNuMuBar  ) return kOscNuS4Bar2NuMuBar;
  if ( nu_start == kPdgNuS5Bar  && nu_end == kPdgNuMuBar  ) return kOscNuS5Bar2NuMuBar;
  if ( nu_start == kPdgNuS4Bar  && nu_end == kPdgNuTauBar ) return kOscNuS4Bar2NuTauBar;
  if ( nu_start == kPdgNuS5Bar  && nu_end == kPdgNuTauBar ) return kOscNuS5Bar2NuTauBar;

  return kOscUndefined;
}
//____________________________________________________________________________
void t2kana::util::PrintCompact(
        const TVectorD & v, int w, int prec, Option_t * option, ostream & stream)
{
  TString opt = option;

  bool printing_banff_params = opt.Contains("BANFF fit params");
  bool printing_fsi_params   = opt.Contains("SKDetFSI params");
  bool printing_selected_banff_params = opt.Contains("BANFF selected fit params");

  int np = v.GetNoElements();
  stream << fixed;

  for(int ip = 0; ip < np; ip++) {
     stream << "> param = " << setw(3) << ip;
     if(printing_banff_params) {
         stream << " [ ";
         PhysParamId_t p = PhysParamId_t ( (int)kNuisBANFF00 + ip );
         stream << setw(42) << util::AsString(p);
         stream << " ]";
     }
     if(printing_fsi_params) {
         stream << " [ ";
#if FSI_BINS != 3
         PhysParamId_t p = PhysParamId_t ( (int)kNuisSKDetFSI000 + ip );
#else
	 PhysParamId_t p = PhysParamId_t ( (int)kNuisSKDetFSI006 + ip );
#endif
         stream << setw(30) << util::AsString(p);
         stream << " ]";
     }
     if(printing_selected_banff_params) {
         stream << " [ ";
         PhysParamId_t p = kPhypUndefined;
	 if(ip == 0)
	   p = kNuisBANFFMAQE;
	 else if(ip == 1)
	   p = kNuisBANFFCCQEE0;
	 else if(ip == 2)
	   p = kNuisBANFFCC1piE0;
         stream << setw(30) << util::AsString(p);
         stream << " ]";
     }
     stream << " - value = " << setw(w) << setprecision(prec) << v[ip] << endl;
  }
  stream << scientific;
}
//___________________________________________________________________________
void t2kana::util::PrintCompact(
   const TVectorD & current, const TVectorD & nominal, 
   int w, int prec, Option_t * option, ostream & stream)
{
  TString opt = option;

  bool printing_banff_params = opt.Contains("BANFF fit params");
  bool printing_fsi_params   = opt.Contains("FSI params");

  int np = current.GetNoElements();
  stream << fixed;

  for(int ip = 0; ip < np; ip++) {
     stream << "> param = " << setw(3) << ip;
     if(printing_banff_params) {
         stream << " [ ";
         PhysParamId_t p = PhysParamId_t ( (int)kNuisBANFF00 + ip );
         stream << setw(42) << util::AsString(p);
         stream << " ]";
     }
     if(printing_fsi_params) {
         stream << " [ ";
#if FSI_BINS != 3
         PhysParamId_t p = PhysParamId_t ( (int)kNuisFSI000 + ip );
#else
	 PhysParamId_t p = PhysParamId_t ( (int)kNuisFSI006 + ip );
#endif
         stream << setw(30) << util::AsString(p);
         stream << " ]";
     }
     stream << " - current = " 
          << setw(w) << setprecision(prec) << current[ip];
     stream << " - nominal = " 
          << setw(w) << setprecision(prec) << nominal[ip];
     if(nominal[ip] != 0) {
       stream << " - current/nominal = " 
            << setw(w) << setprecision(prec) << current[ip]/nominal[ip];
     }
     stream << endl;
  }
  stream << scientific;
}
//___________________________________________________________________________
void t2kana::util::PrintCompact(
  const TMatrixD & m, int w, int prec, ostream & stream)
{
  int nc = m.GetNcols();
  int nr = m.GetNrows();
  stream << fixed;
  for(int ir=0; ir<nr; ir++) {
    if(ir==0) {
       stream << setw(w) << setfill(' ') << "r\\c:" << " | "; 
       for(int ic=0; ic<nc; ic++) {
         stream << setw(w) << setfill(' ') << setprecision(0) << ic << " | ";
       }
       stream << endl;
    }
    stream << setw(w) << setfill(' ') << setprecision(0) << ir << " | "; 
    for(int ic=0; ic<nc; ic++) {
         stream << setw(w) << setfill(' ') << setprecision(prec) 
              << m[ir][ic] << " | ";
    }//rows
    stream << endl;
  }//columns
  stream << scientific;
}
//___________________________________________________________________________
void t2kana::util::PrintCompactTex(
				   const TMatrixD & m, int w, int prec, ostream & stream, TString opt)
{
  int nc = m.GetNcols();
  int nr = m.GetNrows();
  stream << fixed;
  string columntype = "c";
  if(opt.Contains("banff"))
    columntype = "p{18pt}";
  stream << "{" << endl;
  if(opt.Contains("skdetfsi"))
    stream << "\\small" << endl;
  else if(opt.Contains("banff"))
    stream << "\\tiny" << endl
	   << "\\begin{sideways}" << endl
	   << "\\begin{minipage}{23cm}" << endl
	   << "{" << endl
	   << "\\setlength{\\arraycolsep}{1pt}" << endl;
  stream << "\\begin{equation}" << endl
	 << "C^{";
  if(opt.Contains("skdetfsi"))
    stream << "SK+FSI";
  else if(opt.Contains("banff"))
    stream << "banff";
  stream << "} =";
  if(opt.Contains("banff"))
    stream << "\\end{equation}" << endl
	   << "\\begin{equation*}" << endl;  
  else
    stream << " \\\\" << endl;
  stream << "\\begin{blockarray}{c ";
  for(int ic = 0; ic < nc; ic++)
    stream << columntype;
  stream << "}" << endl;

  for(int ir=0; ir<nr; ir++) {
    if(ir==0) {
      if(!opt.Contains("skdetfsi") && !opt.Contains("banff"))
	stream << setw(w) << setfill(' ') << "r\\c:"; 
      else
	stream << setw(w) << setfill(' ');
      for(int ic=0; ic<nc; ic++) {
	if(opt.Contains("skdetfsi"))
	  stream << " & " << setw(w) << setfill(' ') << setprecision(0) << "{\\color{blue} " << util::AsStringPaperTexNoDollar(kListOfSKDetFSIParams[ic]) << "}";
	else if(opt.Contains("banff"))
	  stream << " & " << setw(w) << setfill(' ') << setprecision(0) << "\\begin{sideways} {\\color{blue} " << util::AsStringPaperTexNoDollar((PhysParamId_t)(ic + (int)kNuisBANFF00)) << "} \\end{sideways}";
	else
	  stream << " & " << setw(w) << setfill(' ') << setprecision(0) << "{\\color{blue} " << ic << "}";
       }
       stream << " \\\\ " << endl;
       stream << "\\begin{block}{c(";
       for(int ic = 0; ic < nc; ic++)
	 stream << columntype;
       if(opt.Contains("banff"))
	 stream << "@{\\hspace*{1pt}})}" << endl;
       else
	 stream << "@{\\hspace*{2pt}})}" << endl;
    }
    if(opt.Contains("skdetfsi"))
      stream << setw(w) << setfill(' ') << setprecision(0) << "{\\color{blue} " << util::AsStringPaperTexNoDollar(kListOfSKDetFSIParams[ir]) << "}  \\vspace{0.2cm}";
    else if(opt.Contains("banff"))
      stream << setw(w) << setfill(' ') << setprecision(0) << "{\\color{blue} " << util::AsStringPaperTexNoDollar((PhysParamId_t)(ir + (int)kNuisBANFF00)) << "}  \\vspace{0.2cm}";
    else
      stream << setw(w) << setfill(' ') << setprecision(0) << "{\\color{blue} " << ir << "}  \\vspace{0.2cm}";
    for(int ic=0; ic<nc; ic++) {
      if(ic == ir)
	stream << " & " << setw(w) << setfill(' ') << setprecision(prec) 
	       << "{\\color{red} " << m[ir][ic] << "}";
      else
	stream << " & " << setw(w) << setfill(' ') << setprecision(prec) 
	       << m[ir][ic];
    }//rows
    stream << " \\\\ " << endl;
  }//columns
  stream << "\\end{block}" << endl;
  if(opt.Contains("skdetfsi"))
    stream << "\\label{eq:Cskfsi}" << endl;
  else if(opt.Contains("banff"))
    stream << "\\label{eq:Cbanff}" << endl;
  stream << "\\end{blockarray}" << endl;
  if(opt.Contains("banff"))
    stream << "\\end{equation*}" << endl;
  else
    stream << "\\end{equation}" << endl;
  if(opt.Contains("banff"))
    stream << "}" << endl
	   << "\\end{minipage}" << endl
	   << "\\end{sideways}" << endl;
  stream << "\\normalsize" << endl
	 << "}" << endl;
  stream << scientific;
}
//___________________________________________________________________________

// A few helper functions for use in UseEventToFillInTemplate
// that make swapping mode definitions a little easier

bool t2kana::util::IsCCQE(int neut_reaction_mode) 
{
  return neut_reaction_mode ==  1;
}
//___________________________________________________________________________
bool t2kana::util::IsCC1pi(int neut_reaction_mode) 
{
  bool is_cc1pi = (neut_reaction_mode == 11 || 
		   neut_reaction_mode == 12 || 
		   neut_reaction_mode == 13);
  return is_cc1pi;
}
//___________________________________________________________________________
bool t2kana::util::IsCCcoh(int neut_reaction_mode) 
{
  return neut_reaction_mode ==  16;
}
//___________________________________________________________________________
bool t2kana::util::IsCCoth(int neut_reaction_mode) 
//this is a catch all for everything CC that is not QE, 1pi, or coh
// needs to be updated if another subcategory is added
{
  bool is_ccoth = (neut_reaction_mode == 17 ||
		     neut_reaction_mode == 21 ||
		     neut_reaction_mode == 22 || 
		     neut_reaction_mode == 23 || 
		     neut_reaction_mode == 26);
  return is_ccoth;
}
//___________________________________________________________________________
//bool t2kana::util::IsCCoth(int neut_reaction_mode) 
//{
//return IsCCoth(neut_reaction_mode);
//}
//___________________________________________________________________________
//bool t2kana::util::IsCCoth(int neut_reaction_mode) 
//{
//return IsCCoth(neut_reaction_mode);
//}
//___________________________________________________________________________
bool t2kana::util::IsNC1pi0(int neut_reaction_mode) 
{
  bool is_cc1pi0 = (neut_reaction_mode == 31 ||
		    neut_reaction_mode == 32);
  return is_cc1pi0;
}
//___________________________________________________________________________
bool t2kana::util::IsNC1piC(int neut_reaction_mode) 
{
  bool is_cc1pi0 = (neut_reaction_mode == 33 ||
		    neut_reaction_mode == 34);
  return is_cc1pi0;
}
//___________________________________________________________________________
bool t2kana::util::IsNCcoh(int neut_reaction_mode) 
{
  return neut_reaction_mode == 36;
}
//___________________________________________________________________________
bool t2kana::util::IsNCoth(int neut_reaction_mode) 
//this is a catch all for everything NC that is not 1pi0, 1piC or coh
// needs to be updated if another subcategory is added
{
  bool is_ncoth = (  neut_reaction_mode == 38 ||
		     neut_reaction_mode == 39 ||
		     neut_reaction_mode == 41 ||
		     neut_reaction_mode == 42 ||
		     neut_reaction_mode == 43 ||
		     neut_reaction_mode == 44 ||
		     neut_reaction_mode == 45 ||
		     neut_reaction_mode == 46 ||
		     neut_reaction_mode == 51 ||
		     neut_reaction_mode == 52);
  return is_ncoth;
}
//___________________________________________________________________________
bool t2kana::util::IsNCoth_mu(int neut_reaction_mode) 
{
  bool is_ncoth_mu = ( IsNCoth (neut_reaction_mode) ||
		       IsNCcoh (neut_reaction_mode) ||
		       IsNC1pi0(neut_reaction_mode) );
  return is_ncoth_mu;
}
//___________________________________________________________________________
bool t2kana::util::IsNCoth_e(int neut_reaction_mode) 
{
  bool is_ncoth_e = ( IsNCoth (neut_reaction_mode) ||
		      IsNC1piC(neut_reaction_mode) );
  return is_ncoth_e;
}
//___________________________________________________________________________
bool t2kana::util::UseEventToFillInTemplate(
  SuperkMode_t m, int nu_pdg_code, int neut_reaction_mode)
{

  bool use_it = false;

  //NuMu templates
  if(m == kSKmNuMuCCQE) 
  {
     bool is_numu = (nu_pdg_code == kPdgNuMu);
     use_it = is_numu && IsCCQE(neut_reaction_mode);
  } 

  else
  if(m == kSKmNuMuCC1pi) 
  {
     bool is_numu  = (nu_pdg_code == kPdgNuMu);
     use_it = is_numu && IsCC1pi(neut_reaction_mode);
  } 

  else
  if(m == kSKmNuMuCCcoh)
  {
    bool is_numu  = (nu_pdg_code == kPdgNuMu);
    use_it = is_numu && IsCCcoh(neut_reaction_mode);
  }

  else
  if(m == kSKmNuMuCCoth) 
  {
     bool is_numu     = (nu_pdg_code == kPdgNuMu);
     use_it = is_numu && IsCCoth(neut_reaction_mode);
  } 

  else
  if(m == kSKmNuMuCCoth)
  {
    bool is_numu    = (nu_pdg_code == kPdgNuMu);
    use_it = is_numu && IsCCoth(neut_reaction_mode);
  }

  else
  if(m == kSKmNuMuNC1piC) 
  {
     bool is_numu   = (nu_pdg_code == kPdgNuMu);
     use_it = is_numu && IsNC1piC(neut_reaction_mode);
  } 

  else
  if(m == kSKmNuMuNC1pi0) 
  {
     bool is_numu   = (nu_pdg_code == kPdgNuMu);
     use_it = is_numu && IsNC1pi0(neut_reaction_mode);
  } 

  else
  if(m == kSKmNuMuNCcoh) 
  {
     bool is_numu  = (nu_pdg_code == kPdgNuMu);
     use_it = is_numu && IsNCcoh(neut_reaction_mode);
  } 

  else
  if(m == kSKmNuMuNCoth_e) 
  {
     bool is_numu    = (nu_pdg_code == kPdgNuMu);
     use_it = is_numu && IsNCoth_e(neut_reaction_mode);
  } 

  else
  if(m == kSKmNuMuNCoth_mu) 
  {
     bool is_numu    = (nu_pdg_code == kPdgNuMu);
     use_it = is_numu && IsNCoth_mu(neut_reaction_mode);
  } 

  //NuMuBar templates
  else
  if(m == kSKmNuMuBarCCQE) 
  {
     bool is_numubar = (nu_pdg_code == kPdgNuMuBar);
     use_it = is_numubar && IsCCQE(-neut_reaction_mode);
  } 

  else
  if(m == kSKmNuMuBarCC1pi) 
  {
     bool is_numubar = (nu_pdg_code == kPdgNuMuBar);
     use_it = is_numubar && IsCC1pi(-neut_reaction_mode);
  } 

  else
  if(m == kSKmNuMuBarCCcoh) 
  {
     bool is_numubar = (nu_pdg_code == kPdgNuMuBar);
     use_it = is_numubar && IsCCcoh(-neut_reaction_mode);
  } 

  else
  if(m == kSKmNuMuBarCCoth) 
  {
     bool is_numubar  = (nu_pdg_code == kPdgNuMuBar);
     use_it = is_numubar && IsCCoth(-neut_reaction_mode);
  } 

  else
  if(m == kSKmNuMuBarCCoth) 
  {
     bool is_numubar = (nu_pdg_code == kPdgNuMuBar);
     use_it = is_numubar && IsCCoth(-neut_reaction_mode);
  } 

  else
  if(m == kSKmNuMuBarNC1piC) 
  {
     bool is_numubar  = (nu_pdg_code == kPdgNuMuBar);
     use_it = is_numubar && IsNC1piC(-neut_reaction_mode);
  } 

  else
  if(m == kSKmNuMuBarNC1pi0) 
  {
     bool is_numubar  = (nu_pdg_code == kPdgNuMuBar);
     use_it = is_numubar && IsNC1pi0(-neut_reaction_mode);
  } 

  else
  if(m == kSKmNuMuBarNCcoh) 
  {
     bool is_numubar  = (nu_pdg_code == kPdgNuMuBar);
     use_it = is_numubar && IsNCcoh(-neut_reaction_mode);
  } 

  else
  if(m == kSKmNuMuBarNCoth_e) 
  {
     bool is_numubar  = (nu_pdg_code == kPdgNuMuBar);
     use_it = is_numubar && IsNCoth_e(-neut_reaction_mode);
  } 

  else
  if(m == kSKmNuMuBarNCoth_mu) 
  {
     bool is_numubar  = (nu_pdg_code == kPdgNuMuBar);
     use_it = is_numubar && IsNCoth_mu(-neut_reaction_mode);
  } 

  //NuE (signal and beam) templates
  else
  if(m == kSKmNuECCQE || m == kSKmOscNuECCQE) 
  {
     bool is_nue  = (nu_pdg_code == kPdgNuE); 
     use_it = is_nue && IsCCQE(neut_reaction_mode);
  } 

  else
  if(m == kSKmNuECC1pi || m == kSKmOscNuECC1pi) 
  {
     bool is_nue   = (nu_pdg_code == kPdgNuE); 
     use_it = is_nue && IsCC1pi(neut_reaction_mode);
  }

  else
  if(m == kSKmNuECCcoh || m == kSKmOscNuECCcoh) 
  {
     bool is_nue   = (nu_pdg_code == kPdgNuE); 
     use_it = is_nue && IsCCcoh(neut_reaction_mode);
  } 

  else
  if(m == kSKmNuECCoth || m == kSKmOscNuECCoth) 
  {
     bool is_nue      = (nu_pdg_code == kPdgNuE); 
     use_it = is_nue && IsCCoth(neut_reaction_mode);
  }

  else
  if(m == kSKmNuECCoth || m == kSKmOscNuECCoth) 
  {
     bool is_nue     = (nu_pdg_code == kPdgNuE); 
     use_it = is_nue && IsCCoth(neut_reaction_mode);
  }

  //this is a little different. nuesig has no CCcoh in 1Re sample, whereas CCcoh is included
  // for the other species
  //  else
  //if(m == kSKmOscNuECCoth) 
  //{
  // bool is_nue     = (nu_pdg_code == kPdgNuE); 
  // use_it = is_nue && (IsCCoth(neut_reaction_mode) || IsCCcoh(neut_reaction_mode));
  //}

  else
  if(m == kSKmNuENC1piC) 
  {
     bool is_nue  = (nu_pdg_code == kPdgNuE); 
     use_it = is_nue && IsNC1piC(neut_reaction_mode);
  }

  else
  if(m == kSKmNuENC1pi0) 
  {
     bool is_nue    = (nu_pdg_code == kPdgNuE); 
     use_it = is_nue && IsNC1pi0(neut_reaction_mode);
  }

  else
  if(m == kSKmNuENCcoh) 
  {
     bool is_nue   = (nu_pdg_code == kPdgNuE); 
     use_it = is_nue && IsNCcoh(neut_reaction_mode);
  }

  else
  if(m == kSKmNuENCoth_e) 
  {
     bool is_nue     = (nu_pdg_code == kPdgNuE); 
     use_it = is_nue && IsNCoth_e(neut_reaction_mode);
  }

  else
  if(m == kSKmNuENCoth_mu) 
  {
     bool is_nue     = (nu_pdg_code == kPdgNuE); 
     use_it = is_nue && IsNCoth_mu(neut_reaction_mode);
  }

  //NuEBar (beam) templates
  else
  if(m == kSKmNuEBarCCQE) 
  {
     bool is_nue  = (nu_pdg_code == kPdgNuEBar); 
     use_it = is_nue && IsCCQE(-neut_reaction_mode);
  } 

  else
  if(m == kSKmNuEBarCC1pi) 
  {
     bool is_nue   = (nu_pdg_code == kPdgNuEBar); 
     use_it = is_nue && IsCC1pi(-neut_reaction_mode);
  }

  else
  if(m == kSKmNuEBarCCcoh) 
  {
     bool is_nue   = (nu_pdg_code == kPdgNuEBar); 
     use_it = is_nue && IsCCcoh(-neut_reaction_mode);
  } 


  else
  if(m == kSKmNuEBarCCoth) 
  {
     bool is_nue      = (nu_pdg_code == kPdgNuEBar); 
     use_it = is_nue && IsCCoth(-neut_reaction_mode);
  }

  else
  if(m == kSKmNuEBarCCoth) 
  {
     bool is_nue     = (nu_pdg_code == kPdgNuEBar); 
     use_it = is_nue && IsCCoth(-neut_reaction_mode);
  }

  else
  if(m == kSKmNuEBarNC1piC) 
  {
     bool is_nue  = (nu_pdg_code == kPdgNuEBar); 
     use_it = is_nue && IsNC1piC(-neut_reaction_mode);
  }

  else
  if(m == kSKmNuEBarNC1pi0) 
  {
     bool is_nue    = (nu_pdg_code == kPdgNuEBar); 
     use_it = is_nue && IsNC1pi0(-neut_reaction_mode);
  }

  else
  if(m == kSKmNuEBarNCcoh) 
  {
     bool is_nue   = (nu_pdg_code == kPdgNuEBar); 
     use_it = is_nue && IsNCcoh(-neut_reaction_mode);
  }

  else
  if(m == kSKmNuEBarNCoth_e) 
  {
     bool is_nue     = (nu_pdg_code == kPdgNuEBar); 
     use_it = is_nue && IsNCoth_e(-neut_reaction_mode);
  }

  else
  if(m == kSKmNuEBarNCoth_mu) 
  {
     bool is_nue     = (nu_pdg_code == kPdgNuEBar); 
     use_it = is_nue && IsNCoth_mu(-neut_reaction_mode);
  }

  return use_it;
}
//___________________________________________________________________________
t2kana::util::ArrIdx * t2kana::util::ArrIdx::fInstance = 0;
//____________________________________________________________________________
t2kana::util::ArrIdx::ArrIdx()
{
  this->Init();
  fInstance = 0;
}
//____________________________________________________________________________
t2kana::util::ArrIdx * t2kana::util::ArrIdx::Instance()
{
  if(fInstance == 0) {
    static t2kana::util::ArrIdx::Cleaner cleaner;
    cleaner.DummyMethodAndSilentCompiler();
    fInstance = new t2kana::util::ArrIdx;
  }
  return fInstance;
}
//____________________________________________________________________________
void t2kana::util::ArrIdx::Init(void)
{
  //
  // sample/mode arrays
  //
  for(int i = 0; i < 2*kns; i++) {
     fSuperkSampleIdx [i] = -1;
     for(int j = 0; j < 2*knm; j++) {
         fSuperkModeIdx [i][j] = -1;
     }
  }
  for(int is = 0; is < kns; is++) {
     SuperkSample_t s = kSKSamples[is];
     if(s == kSKsUndefined) continue;
     fSuperkSampleIdx [s] = is;
#ifdef __show_arr_idx_mesg_debug__
     cout << "Util/ArrIdx >> Sample: " 
          << util::AsString(s)
          << " is at slot: " << fSuperkSampleIdx [s] << endl;
#endif
     for(int im = 0; im < knm; im++) {
         SuperkMode_t m = kSKReactionModes[is][im];
         if(m == kSKmUndefined) continue;
         fSuperkModeIdx [s][m] = im;
#ifdef __show_arr_idx_mesg_debug__
         cout << "Util/ArrIdx >>   Mode: " 
              << util::AsString(s) << "/" << util::AsString(m)
              << " is at slot: " << fSuperkModeIdx [s][m] << endl;
#endif
     }
  } 
  // 
  // phys param arrays
  // 
  //XSEC
  for(int i = 0; i < knppmax; i++) {
    fXSecSystIdx[i] = -1;
  }
  for(int ip = 0; ip < kNumOfXSecSyst; ip++) {
    PhysParamId_t p = kListOfXSecSyst[ip];
    if(p == kPhypUndefined) continue;
    fXSecSystIdx[p] = ip;
#ifdef __show_arr_idx_mesg_debug__
         cout << "Util/ArrIdx >> FHC Systematic: " 
              << util::AsString(p) 
              << " is at slot: " << fXSecSystIdx[p] << endl;
#endif
  }
  for(int i = 0; i < knppmax; i++) {
    fRHCXSecSystIdx[i] = -1;
  }
  for(int ip = 0; ip < kNumOfRHCXSecSyst; ip++) {
    PhysParamId_t p = kListOfRHCXSecSyst[ip];
    if(p == kPhypUndefined) continue;
    fRHCXSecSystIdx[p] = ip;
#ifdef __show_arr_idx_mesg_debug__
         cout << "Util/ArrIdx >> RHC Systematic: " 
              << util::AsString(p) 
              << " is at slot: " << fRHCXSecSystIdx[p] << endl;
#endif
  }

  //BANFF
  for(int i = 0; i < knppmax; i++) {
    fBANFFSystIdx[i] = -1;
  }
  for(int ip = 0; ip < kNumOfParamsInBANFFfit; ip++) {
    PhysParamId_t p = PhysParamId_t ((int) kNuisBANFF00 + ip);
    if(p == kPhypUndefined) continue;
    fBANFFSystIdx[p] = ip;
#ifdef __show_arr_idx_mesg_debug__
     cout << "Util/ArrIdx >> BANFF: " 
          << util::AsString(p)
          << " is at slot: " << fBANFFSystIdx [p] << endl;
#endif
  }
  for(int i = 0; i < knppmax; i++) {
    fRHCBANFFSystIdx[i] = -1;
  }
  for(int ip = 0; ip < kNumOfRHCParamsInBANFFfit; ip++) {
    PhysParamId_t p = PhysParamId_t ((int) kRHCNuisBANFF00 + ip);
    if(p == kPhypUndefined) continue;
    fRHCBANFFSystIdx[p] = ip;
#ifdef __show_arr_idx_mesg_debug__
    cout << "Util/ArrIdx >> RHC BANFF: "
	 << util::AsString(p)
	 << " is at slot: " << fRHCBANFFSystIdx [p] << endl;
#endif
  }

#ifndef SKFSITOGETHER
  //FSI
  for(int i = 0; i < knppmax; i++) {
    fFSISystIdx[i] = -1;
  }
  for(int ip = 0; ip < kNumOfFSIParams; ip++) {
#if FSI_BINS != 3
    PhysParamId_t p = PhysParamId_t ((int) kNuisFSI000 + ip);
#else
    PhysParamId_t p = PhysParamId_t ((int) kNuisFSI006 + ip);
#endif
    if(p == kPhypUndefined) continue;
    fFSISystIdx[p] = ip;
#ifdef __show_arr_idx_mesg_debug__
     cout << "Util/ArrIdx >> FSI: " 
          << util::AsString(p)
          << " is at slot: " << fFSISystIdx [p] << endl;
#endif
  }
#endif

  //SKDET + FSI
#ifdef SKFSITOGETHER
  for(int i = 0; i < knppmax; i++) {
    fSKDetFSISystIdx[i] = -1;
  }
  for(int ip = 0; ip < kNumOfSKDetFSIParams; ip++) {
#if FSI_BINS != 3
    PhysParamId_t p = PhysParamId_t ((int) kNuisSKDetFSI000 + ip);
#else
    PhysParamId_t p = PhysParamId_t ((int) kNuisSKDetFSI006 + ip);
#endif
    if(p == kPhypUndefined) continue;
    fSKDetFSISystIdx[p] = ip;
#ifdef __show_arr_idx_mesg_debug__
     cout << "Util/ArrIdx >> SKDet+FSI: " 
          << util::AsString(p)
          << " is at slot: " << fSKDetFSISystIdx [p] << endl;
#endif
  }
#endif
#ifdef SKFSITOGETHER
  for(int i = 0; i < knppmax; i++) {
    fRHCSKDetFSISystIdx[i] = -1;
  }
  for(int ip = 0; ip < kNumOfRHCSKDetFSIParams; ip++) {
#if FSI_BINS != 3
    PhysParamId_t p = PhysParamId_t ((int) kRHCNuisSKDetFSI000 + ip);
#else
    PhysParamId_t p = PhysParamId_t ((int) kRHCNuisSKDetFSI006 + ip);
#endif
    if(p == kPhypUndefined) continue;
    fRHCSKDetFSISystIdx[p] = ip;
#ifdef __show_arr_idx_mesg_debug__
    cout << "Util/ArrIdx >> RHC SKDet+FSI: "
	 << util::AsString(p)
	 << " is at slot: " << fRHCSKDetFSISystIdx [p] << endl;
#endif
  }
#endif

#ifndef SKFSITOGETHER
#if FSI_BINS != 2
  //1Re SK eff
  for(int i = 0; i < knppmax; i++) {
    fSKDet1ReSystIdx[i] = -1;
  }
  for(int ip = 0; ip < kNumOfSKDetParams1Re; ip++) {
    PhysParamId_t p = PhysParamId_t ((int) kNuisSuperkEff1ReNueSigCC0 + ip);
    if(p == kPhypUndefined) continue;
    fSKDet1ReSystIdx[p] = ip;
#ifdef __show_arr_idx_mesg_debug__
     cout << "Util/ArrIdx >> 1Re SK Det: " 
          << util::AsString(p)
          << " is at slot: " << fSKDet1ReSystIdx [p] << endl;
#endif
  }
#endif
#endif

  //PhysParam
  for(int i = 0; i < knppmax; i++) {
    fPhysParamIdx[i] = -1;
  }
  for(int ip = 0; ip < knppmax; ip++) {
    PhysParamId_t p = PhysParamId_t(ip);
    if(p == kPhypUndefined) continue;
    fPhysParamIdx[p] = ip;
#ifdef __show_arr_idx_mesg_debug__
     cout << "Util/ArrIdx >> PhysParam: " 
          << util::AsString(p)
          << " is at slot: " << fPhysParamIdx [p] << endl;
#endif
  }
}
//____________________________________________________________________________
