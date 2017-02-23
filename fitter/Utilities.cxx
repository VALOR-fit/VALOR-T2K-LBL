//____________________________________________________________________________
/*

\copyright (c) 2010-2013, VaLOR group.

30-12-2013  Nick
Added samples and modes for reverse horn current running

21-7-2014 Davide
Added samples for multiring forward horn current

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

#include "valor/Utilities/Rnd.h"

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

using namespace valor;
using namespace valor::conv;

//____________________________________________________________________________
string valor::util::AsString(valor::conv::OscAnaType_t a)
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
string valor::util::AsString(conv::SKMCTemplateType_t t)
{
  switch(t) {
    case ( kSKMCtOri ) : 
      return "straight from the MC trees";
      break;
    case ( kSKMCtPOT ) : 
      return "flux and POT weighted, no NIWG weights";
      break;
    case ( kSKMCtNom ) : 
      return "nominal with weights from SKMC branches";
      break;
    case ( kSKMCtRPA ) : 
      return "nominal + rel RPA";
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
string valor::util::AsString(valor::conv::FluxType_t f)
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
string valor::util::AsStringRooTex (
  valor::conv::SuperkSample_t s, valor::conv::SuperkMode_t m, 
  string append_begin, string append_end)
{
  ostringstream out;

  if(append_begin.size() > 0) {
    out << append_begin;
  }

  out << valor::util::AsStringRooTex(s);

  if(m != kSKmUndefined) {
    out << " (";
    out << valor::util::AsStringRooTex(m);
    out << ")";
  }

  if(append_end.size() > 0) {
    out << append_end;
  }

  return out.str();
}
//____________________________________________________________________________
string valor::util::AsStringCompact (
  valor::conv::SuperkSample_t s, valor::conv::SuperkMode_t m, 
  string append_begin, string append_end)
{
// Compact string used for naming output histograms stored ROOT files
//
  ostringstream out;

  if(append_begin.size() > 0) {
    out << append_begin;
  }

  out << valor::util::AsStringCompact(s);

  if(m != kSKmUndefined) {
    out << "_";
    out << valor::util::AsStringCompact(m);
  }

  if(append_end.size() > 0) {
    out << append_end;
  }

  return out.str();
}
//____________________________________________________________________________
string valor::util::AsStringCompact_splinehack (valor::conv::SuperkSample_t s)
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
//____________________________________________________________________________
valor::conv::OscAnaType_t  
   valor::util::OscAnaFromString(string sh) 
{
  if(sh == "combined"   ) return kAnaCombined;
  if(sh == "numu_disap" ) return kAnaNuMuDisappOnly;
  if(sh == "nue_app"    ) return kAnaNuEAppOnly;
  
  return kAnaUndefined;
}
//____________________________________________________________________________
valor::conv::SKMCTemplateType_t   
  valor::util::SKMCTemplateTypeFromString(string sh) 
{
  if(sh == "ori")        return kSKMCtOri;
  if(sh == "nom")        return kSKMCtNom;
  if(sh == "pot")        return kSKMCtPOT;
  if(sh == "rpa")        return kSKMCtRPA;
  if(sh == "tun")        return kSKMCtTun;
  if(sh == "sys")        return kSKMCtSys;
  if(sh == "fin")        return kSKMCtFin;

  return kSKMCtUndefined;
}
//____________________________________________________________________________
void valor::util::PrintCompact(
        const TVectorD & v, int w, int prec, Option_t * option, std::ostream & stream)
{
  TString opt = option;

  bool printing_banff_params = opt.Contains("BANFF fit params");
  bool printing_fsi_params   = opt.Contains("SKDetFSI params");
  /// NEW MULTIRING 2
  bool printing_multiring_fsi_params   = opt.Contains("MultiRing SKDetFSI fit params");
  ///
  bool printing_selected_banff_params = opt.Contains("BANFF selected fit params");

  int np = v.GetNoElements();
  stream << fixed;

  for(int ip = 0; ip < np; ip++) {
     stream << "> param = " << setw(3) << ip;
     if(printing_banff_params) {
         stream << " [ ";

         //PhysParamId_t p = kListOfBANFFParams[ip];
         
         // NEW FHC+RHC JOINT
       #if defined(__FHC_RUNNING__) && defined(__RHC_RUNNING__)  
         PhysParamId_t p;
         if      (opt.Contains("FHC+RHC")) p = kListOfFHCpRHCBANFFParams[ip];
         else if (opt.Contains("FHC"))     p = kListOfBANFFParams[ip];   
         else if (opt.Contains("RHC"))     p = kListOfRHCBANFFParams[ip];   
         else {
           cerr << "BANFF options string doesn't contain any options" << endl;
           exit(1);
         }
         //

       #elif defined(__FHC_RUNNING__)
         PhysParamId_t p = kListOfBANFFParams[ip];
       #elif defined(__MULTIRING_FHC_RUNNING__) // ~#{ MultiRing
         PhysParamId_t p = kListOfBANFFParams[ip];
       #elif defined(__RHC_RUNNING__)
         PhysParamId_t p = kListOfRHCBANFFParams[ip];
       #endif
         //
         
         stream << setw(42) << util::AsString(p);
         stream << " ]";
     }
     if(printing_fsi_params) {
         stream << " [ ";


       // #if FSI_BINS != 3
       //   //PhysParamId_t p = PhysParamId_t ( (int)kNuisSKDetFSI000 + ip );
       //        // NEW FLUX
       // #if defined(__FHC_RUNNING__) && defined(__RHC_RUNNING__) && defined(__FULLY_CORRELATED__) // fhc+rhc fully correlated
       //        PhysParamId_t p = PhysParamId_t ( (int)kNuisSKDetFSI000 + ip );
       // #elif defined(__FHC_RUNNING__) // fhc only
       //   PhysParamId_t p = PhysParamId_t ( (int)kNuisSKDetFSI000 + ip );
       // #elif defined(__RHC_RUNNING__) // rhc only
       //   PhysParamId_t p = PhysParamId_t ( (int)kRHCNuisSKDetFSI000 + ip );
       // #endif
       //        //
       // #else
       //        //PhysParamId_t p = PhysParamId_t ( (int)kNuisSKDetFSI006 + ip );
       //        // NEW FLUX
       // #if defined(__FHC_RUNNING__) && defined(__RHC_RUNNING__) && defined(__FULLY_CORRELATED__) // fhc+rhc fully correlated
       //        PhysParamId_t p = PhysParamId_t ( (int)kNuisSKDetFSI006 + ip );
       // #elif defined(__FHC_RUNNING__) // fhc only
       //        PhysParamId_t p = PhysParamId_t ( (int)kNuisSKDetFSI006 + ip );
       // #elif defined(__RHC_RUNNING__) // rhc only
       //        PhysParamId_t p = PhysParamId_t ( (int)kRHCNuisSKDetFSI006 + ip );
       // #endif
       //        
       // #endif
 
                 // NEW FLUX

         // NEW SKDETFSI (use kNuisSKDetFSI_first)

         //PhysParamId_t p = PhysParamId_t ( (int)kNuisSKDetFSI_first + ip );
       #if defined(__FHC_RUNNING__) && defined(__RHC_RUNNING__) && defined(__FULLY_CORRELATED__) // fhc+rhc fully correlated
         PhysParamId_t p = PhysParamId_t ( (int)kNuisSKDetFSI_first + ip );
       #elif defined(__FHC_RUNNING__) // fhc only
         PhysParamId_t p = PhysParamId_t ( (int)kNuisSKDetFSI_first + ip );
       #elif defined(__MULTIRING_FHC_RUNNING__) // ~#{ MultiRing
         PhysParamId_t p = PhysParamId_t ( (int)kNuisSKDetFSI_first + ip );
       #elif defined(__RHC_RUNNING__) // rhc only
         PhysParamId_t p = PhysParamId_t ( (int)kRHCNuisSKDetFSI_first + ip );
       #endif
         //

         stream << setw(30) << util::AsString(p);
         stream << " ]";
     }
 
     /// NEW MULTIRING 2
     if(printing_multiring_fsi_params) {
         stream << " [ ";

         // #if FSI_BINS != 3
         //          PhysParamId_t p = PhysParamId_t ( (int)kMultiRingNuisSKDetFSI000 + ip );
         // #else
         //      PhysParamId_t p = PhysParamId_t ( (int)kMultiRingNuisSKDetFSI006 + ip );
         // #endif

         // NEW SKDETFSI (use kMultiRingNuisSKDetFSI_first)
         PhysParamId_t p = PhysParamId_t ( (int)kMultiRingNuisSKDetFSI_first + ip );

         stream << setw(30) << util::AsString(p);
         stream << " ]";
     }
     ///

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
void valor::util::PrintCompact(
   const TVectorD & current, const TVectorD & nominal, 
   int w, int prec, Option_t * option, std::ostream & stream)
{
  TString opt = option;

  bool printing_banff_params = opt.Contains("BANFF fit params");
  bool printing_fsi_params   = opt.Contains("FSI params");
  /// NEW MULTIRING
  bool printing_multiring_fsi_params   = opt.Contains("MultiRing FSI params");
  ///

  int np = current.GetNoElements();
  stream << fixed;

  for(int ip = 0; ip < np; ip++) {
     stream << "> param = " << setw(3) << ip;
     if(printing_banff_params) {
         stream << " [ ";

         //PhysParamId_t p = kListOfBANFFParams[ip];
         
         // NEW FLUX

         // NEW FHC+RHC JOINT
       #if defined(__FHC_RUNNING__) && defined(__RHC_RUNNING__) 
         PhysParamId_t p;
         if      (opt.Contains("FHC+RHC")) p = kListOfFHCpRHCBANFFParams[ip];
         else if (opt.Contains("FHC"))     p = kListOfBANFFParams[ip];   
         else if (opt.Contains("RHC"))     p = kListOfRHCBANFFParams[ip];   
         else {
           cerr << "BANFF options string doesn't contain any options" << endl;
           exit(1);
         }
         //

       #elif defined(__FHC_RUNNING__)
         PhysParamId_t p = kListOfBANFFParams[ip];
       #elif defined(__MULTIRING_FHC_RUNNING__) // ~#{ MultiRing
         PhysParamId_t p = kListOfBANFFParams[ip];
       #elif defined(__RHC_RUNNING__)
         PhysParamId_t p = kListOfRHCBANFFParams[ip];
       #endif
         //

         stream << setw(42) << util::AsString(p);
         stream << " ]";
     }
     if(printing_fsi_params) {
         stream << " [ ";


         // #if FSI_BINS != 3
         //          //PhysParamId_t p = PhysParamId_t ( (int)kNuisSKDetFSI000 + ip );
         //      // NEW FLUX
         //        #if defined(__FHC_RUNNING__) && defined(__RHC_RUNNING__) && defined(__FULLY_CORRELATED__) // fhc+rhc fully correlated
         //          PhysParamId_t p = PhysParamId_t ( (int)kNuisSKDetFSI000 + ip );
         //        #elif defined(__FHC_RUNNING__)
         //          PhysParamId_t p = PhysParamId_t ( (int)kNuisSKDetFSI000 + ip );
         //        #elif defined(__RHC_RUNNING__)
         //          PhysParamId_t p = PhysParamId_t ( (int)kRHCNuisSKDetFSI000 + ip );
         //        #endif
         //      //
         // #else
         //      //PhysParamId_t p = PhysParamId_t ( (int)kNuisSKDetFSI006 + ip );       
         //      // NEW FLUX
         //        #if defined(__FHC_RUNNING__) && defined(__RHC_RUNNING__) && defined(__FULLY_CORRELATED__) // fhc+rhc fully correlated
         //      PhysParamId_t p = PhysParamId_t ( (int)kNuisSKDetFSI006 + ip );
         //        #elif defined(__FHC_RUNNING__)        
         //      PhysParamId_t p = PhysParamId_t ( (int)kNuisSKDetFSI006 + ip );
         //        #elif defined(__RHC_RUNNING__)
         //      PhysParamId_t p = PhysParamId_t ( (int)kRHCNuisSKDetFSI006 + ip );
         //        #endif
         //      //      
         // #endif



         // NEW FLUX

         //PhysParamId_t p = PhysParamId_t ( (int)kNuisSKDetFSI000 + ip );

         // NEW SKDETFSI (use kNuisSKDetFSI_first)
       #if defined(__FHC_RUNNING__) && defined(__RHC_RUNNING__) && defined(__FULLY_CORRELATED__) // fhc+rhc fully correlated
         PhysParamId_t p = PhysParamId_t ( (int)kNuisSKDetFSI_first + ip );
       #elif defined(__FHC_RUNNING__)
         PhysParamId_t p = PhysParamId_t ( (int)kNuisSKDetFSI_first + ip );
       #elif defined(__MULTIRING_FHC_RUNNING__) // ~#{ MultiRing         
         PhysParamId_t p = PhysParamId_t ( (int)kNuisSKDetFSI_first + ip );
       #elif defined(__RHC_RUNNING__)
         PhysParamId_t p = PhysParamId_t ( (int)kRHCNuisSKDetFSI_first + ip );
       #endif

         stream << setw(30) << util::AsString(p);
         stream << " ]";
     }

     /// NEW MULTIRING 2
     if(printing_multiring_fsi_params) {
       stream << " [ ";

       // #if FSI_BINS != 3
       //        PhysParamId_t p = PhysParamId_t ( (int)kMultiRingNuisSKDetFSI000 + ip );
       // #else
       //        PhysParamId_t p = PhysParamId_t ( (int)kMultiRingNuisSKDetFSI006 + ip );
       // #endif
       
       // NEW SKDETFSI (use kMultiRingSKDetFSI_first)
       PhysParamId_t p = PhysParamId_t ( (int)kMultiRingNuisSKDetFSI_first + ip );

       stream << setw(30) << util::AsString(p);
       stream << " ]";
     }
     ///

     stream << " - current = " 
          << setw(w) << setprecision(prec) << current[ip];
     stream << " - nominal = " 
          << setw(w) << setprecision(prec) << nominal[ip];
     if(fabs(nominal[ip]) > conv::kASmallNum) {
       stream << " - current/nominal = " 
            << setw(w) << setprecision(prec) << current[ip]/nominal[ip];
     }
     stream << endl;
  }
  stream << scientific;
}
//___________________________________________________________________________
// Not in use - commented out April 2016. If we ever need this again, rewrite it generically and put it in the framework.
//~ void valor::util::PrintCompactTex(
                                   //~ const TMatrixD & m, int w, int prec, std::ostream & stream, TString opt)
//~ {
  //~ int nc = m.GetNcols();
  //~ int nr = m.GetNrows();
  //~ stream << fixed;
  //~ string columntype = "c";
  //~ if(opt.Contains("banff"))
    //~ columntype = "p{18pt}";
  //~ stream << "{" << endl;
  //~ if(opt.Contains("skdetfsi"))
    //~ stream << "\\small" << endl;
  //~ else if(opt.Contains("banff"))
    //~ stream << "\\tiny" << endl
           //~ << "\\begin{sideways}" << endl
           //~ << "\\begin{minipage}{23cm}" << endl
           //~ << "{" << endl
           //~ << "\\setlength{\\arraycolsep}{1pt}" << endl;
  //~ stream << "\\begin{equation}" << endl
         //~ << "C^{";
  //~ if(opt.Contains("skdetfsi"))
    //~ stream << "SK+FSI";
  //~ /// NEW MULTIRING 2
  //~ else if(opt.Contains("multiring_skdetfsi"))
    //~ stream << "MultiRing SK+FSI";
  //~ ///
  //~ else if(opt.Contains("banff"))
    //~ stream << "banff";
  //~ stream << "} =";
  //~ if(opt.Contains("banff"))
    //~ stream << "\\end{equation}" << endl
           //~ << "\\begin{equation*}" << endl;  
  //~ else
    //~ stream << " \\\\" << endl;
  //~ stream << "\\begin{blockarray}{c ";
  //~ for(int ic = 0; ic < nc; ic++)
    //~ stream << columntype;
  //~ stream << "}" << endl;
//~ 
  //~ for(int ir=0; ir<nr; ir++) {
    //~ if(ir==0) {
      //~ if(!opt.Contains("skdetfsi") && !opt.Contains("banff"))
        //~ stream << setw(w) << setfill(' ') << "r\\c:"; 
      //~ else
        //~ stream << setw(w) << setfill(' ');
      //~ for(int ic=0; ic<nc; ic++) {
        //~ if(opt.Contains("skdetfsi"))
          //~ stream << " & " << setw(w) << setfill(' ') << setprecision(0) << "{\\color{blue} " << util::AsStringPaperTexNoDollar(kListOfSKDetFSIParams[ic]) << "}";
        //~ /// NEW MULTIRING 2
        //~ else if(opt.Contains("multiring_skdetfsi"))
          //~ stream << " & " << setw(w) << setfill(' ') << setprecision(0) << "{\\color{blue} " << util::AsStringPaperTexNoDollar(kListOfMultiRingSKDetFSIParams[ic]) << "}";
        //~ ///
        //~ else if(opt.Contains("banff"))
          //~ stream << " & " << setw(w) << setfill(' ') << setprecision(0) << "\\begin{sideways} {\\color{blue} " << util::AsStringPaperTexNoDollar((PhysParamId_t)(ic + (int)kNuisBANFF00)) << "} \\end{sideways}";
        //~ else
          //~ stream << " & " << setw(w) << setfill(' ') << setprecision(0) << "{\\color{blue} " << ic << "}";
       //~ }
       //~ stream << " \\\\ " << endl;
       //~ stream << "\\begin{block}{c(";
       //~ for(int ic = 0; ic < nc; ic++)
         //~ stream << columntype;
       //~ if(opt.Contains("banff"))
         //~ stream << "@{\\hspace*{1pt}})}" << endl;
       //~ else
         //~ stream << "@{\\hspace*{2pt}})}" << endl;
    //~ }
    //~ if(opt.Contains("skdetfsi"))
      //~ stream << setw(w) << setfill(' ') << setprecision(0) << "{\\color{blue} " << util::AsStringPaperTexNoDollar(kListOfSKDetFSIParams[ir]) << "}  \\vspace{0.2cm}";
    //~ /// NEW MULTIRING 2
    //~ else if(opt.Contains("multiring_skdetfsi"))
      //~ stream << setw(w) << setfill(' ') << setprecision(0) << "{\\color{blue} " << util::AsStringPaperTexNoDollar(kListOfMultiRingSKDetFSIParams[ir]) << "}  \\vspace{0.2cm}";
    //~ ///
    //~ else if(opt.Contains("banff"))
      //~ stream << setw(w) << setfill(' ') << setprecision(0) << "{\\color{blue} " << util::AsStringPaperTexNoDollar((PhysParamId_t)(ir + (int)kNuisBANFF00)) << "}  \\vspace{0.2cm}";
    //~ else
      //~ stream << setw(w) << setfill(' ') << setprecision(0) << "{\\color{blue} " << ir << "}  \\vspace{0.2cm}";
    //~ for(int ic=0; ic<nc; ic++) {
      //~ if(ic == ir)
        //~ stream << " & " << setw(w) << setfill(' ') << setprecision(prec) 
               //~ << "{\\color{red} " << m[ir][ic] << "}";
      //~ else
        //~ stream << " & " << setw(w) << setfill(' ') << setprecision(prec) 
               //~ << m[ir][ic];
    //~ }//rows
    //~ stream << " \\\\ " << endl;
  //~ }//columns
  //~ stream << "\\end{block}" << endl;
  //~ if(opt.Contains("skdetfsi"))
    //~ stream << "\\label{eq:Cskfsi}" << endl;
  //~ /// NEW MULTIRING 2
  //~ else if(opt.Contains("multiring_skdetfsi"))
    //~ stream << "\\label{eq:Cskfsi}" << endl;
  //~ ///
  //~ else if(opt.Contains("banff"))
    //~ stream << "\\label{eq:Cbanff}" << endl;
  //~ stream << "\\end{blockarray}" << endl;
  //~ if(opt.Contains("banff"))
    //~ stream << "\\end{equation*}" << endl;
  //~ else
    //~ stream << "\\end{equation}" << endl;
  //~ if(opt.Contains("banff"))
    //~ stream << "}" << endl
           //~ << "\\end{minipage}" << endl
           //~ << "\\end{sideways}" << endl;
  //~ stream << "\\normalsize" << endl
         //~ << "}" << endl;
  //~ stream << scientific;
//~ }
//___________________________________________________________________________
void valor::util::ApplyEnergyScaleSystematic(conv::SuperkSample_t s, double * events_per_bin, double escale_tweak)
{
  const int n_ereco_bins = util::NumErecoBins(s);
  
  const double * bin_edges = util::ErecoBinEdges(s);
  
  double new_events_per_bin [n_ereco_bins+1];// With overflow bin
  double scaled_bin_edges   [n_ereco_bins+1];
  double unscaled_bin_edges [n_ereco_bins+2];// With overflow bin
  
  // Calculate bin information
  double total_events = 0;
  for (int ibin = 0 ; ibin <= n_ereco_bins; ibin++) {
    if (ibin  != n_ereco_bins) {
      total_events += events_per_bin[ibin];
    }
    unscaled_bin_edges[ibin] = bin_edges[ibin];
    scaled_bin_edges  [ibin] = bin_edges[ibin] * escale_tweak;
    new_events_per_bin[ibin] = 0;
  }
  
  unscaled_bin_edges[n_ereco_bins+1] = 1E4; //Overflow bin (no underflow unless we start binning from something other than zero).
  
  if (total_events < kASmallNum) return;
  
  double overlap = 0;
      
  for (int ibin_s = 0 ; ibin_s < n_ereco_bins ; ibin_s++) { // Loop over scaled bins
    double bin_min_s = scaled_bin_edges[ibin_s];
    double bin_max_s = scaled_bin_edges[ibin_s+1];
    for (int ibin_u = 0 ; ibin_u < n_ereco_bins+1 ; ibin_u++) {
      double bin_min_u = unscaled_bin_edges[ibin_u];
      double bin_max_u = unscaled_bin_edges[ibin_u+1];
      
      if     (bin_min_u >= bin_min_s && bin_max_u <= bin_max_s) overlap = bin_max_u - bin_min_u; // Entirely subsumed unscaled bin
      else if(bin_min_s >= bin_min_u && bin_max_s <= bin_max_u) overlap = bin_max_s - bin_min_s; // Entirely subsumed scaled bin
      else if(bin_max_u >  bin_max_s && bin_min_u >= bin_min_s) overlap = bin_max_s - bin_min_u; // Scaled bin extends low
      else if(bin_min_u <  bin_min_s && bin_max_u <= bin_max_s) overlap = bin_max_u - bin_min_s; // Scaled bin extends high
      else if(bin_min_u >  bin_max_s) break ; //Unscaled bin entirely larger than scaled bin, stop checking higher unscaled
      else overlap = -999; // No crossover
      
      if (overlap <= 0) continue;
      
      double new_events = events_per_bin[ibin_s] * overlap / (scaled_bin_edges[ibin_s+1]-scaled_bin_edges[ibin_s]);
      new_events_per_bin[ibin_u] += new_events;
    }
  }
  // Move events back into main array.
  for (int i = 0 ; i < n_ereco_bins ; i++) {
    events_per_bin[i] = new_events_per_bin[i];
  }
  events_per_bin[n_ereco_bins-1] += new_events_per_bin[n_ereco_bins]; //Put events from overflow bin back into top bin
  
#ifndef NDEBUG
  double new_total_events = 0;
  for (int i = 0 ; i < n_ereco_bins ; i++) {
    new_total_events+=events_per_bin[i];
  }
  if ( fabs((new_total_events/total_events)-1) > conv::kASmallNum ) {
    cerr<<"Util/ApplyEnergyScaleSystematic - failed normalisation check for sample "<<util::AsString(s)<<endl;
    cerr<<"Events: "<<total_events<<" ==> "<<new_total_events<<"\t"<<new_total_events/total_events
        <<"\t"<<fabs((new_total_events/total_events)-1)<<"\t at scale "<<escale_tweak<<endl;
    assert(!"Broken");
  }
#endif //NDEBUG

  return;
}
//___________________________________________________________________________










// NEW 2D TEMPLATES

void valor::util::ApplyEnergyScaleSystematic2D(conv::SuperkSample_t s, double * events_per_bin_array, double escale_tweak)
{
  const int n_q2_bins = util::NumQ2Bins(s);
  const int n_ereco_bins = util::NumErecoBins(s);

  typedef double (erq2type)  [n_ereco_bins][n_q2_bins];
  erq2type & events_per_bin = *reinterpret_cast< erq2type * > ( events_per_bin_array );
  if (!events_per_bin){
    cerr << endl;
    cerr << "valor::util::ApplyEnergyScaleSystematic2D" << endl;
    cerr << "events_per_bin array doesn't exist!" << endl;
    cerr << endl;
    exit(1);
  }

  // NB: events_par_bin is fine!

  // for (int ibin_q2 = 0 ; ibin_q2 < n_q2_bins; ibin_q2++) { // loop over q2
  //   for (int ibin = 0 ; ibin <= n_ereco_bins; ibin++) {
  //     cout << events_per_bin[ibin][ibin_q2] << ", ";
  //   }
  //   cout << endl;
  // }
  // exit(1);

  const double * bin_edges = util::ErecoBinEdges(s);
  
  double new_events_per_bin [n_ereco_bins+1][n_q2_bins];// With overflow bin
  double scaled_bin_edges   [n_ereco_bins+1][n_q2_bins];
  double unscaled_bin_edges [n_ereco_bins+2][n_q2_bins];// With overflow bin
  
  // Calculate bin information

  for (int ibin_q2 = 0 ; ibin_q2 < n_q2_bins; ibin_q2++) { // loop over q2

    double total_events = 0;
    
    for (int ibin = 0 ; ibin <= n_ereco_bins; ibin++) {
      if (ibin  != n_ereco_bins) {
        total_events += events_per_bin[ibin][ibin_q2];
      }
      unscaled_bin_edges[ibin][ibin_q2] = bin_edges[ibin];
      scaled_bin_edges  [ibin][ibin_q2] = bin_edges[ibin] * escale_tweak;
      new_events_per_bin[ibin][ibin_q2] = 0;
    } // ibin (ereco)
    
    unscaled_bin_edges[n_ereco_bins+1][ibin_q2] = 1E4; //Overflow bin (no underflow unless we start binning from something other than zero).
    
    if (total_events < kASmallNum) return;
    
    double overlap = 0;
      
    for (int ibin_s = 0 ; ibin_s < n_ereco_bins ; ibin_s++) { // Loop over scaled bins
      double bin_min_s = scaled_bin_edges[ibin_s][ibin_q2];
      double bin_max_s = scaled_bin_edges[ibin_s+1][ibin_q2];
      for (int ibin_u = 0 ; ibin_u < n_ereco_bins+1 ; ibin_u++) {
        double bin_min_u = unscaled_bin_edges[ibin_u][ibin_q2];
        double bin_max_u = unscaled_bin_edges[ibin_u+1][ibin_q2];
      
        if     (bin_min_u >= bin_min_s && bin_max_u <= bin_max_s) overlap = bin_max_u - bin_min_u; // Entirely subsumed unscaled bin
        else if(bin_min_s >= bin_min_u && bin_max_s <= bin_max_u) overlap = bin_max_s - bin_min_s; // Entirely subsumed scaled bin
        else if(bin_max_u >  bin_max_s && bin_min_u >= bin_min_s) overlap = bin_max_s - bin_min_u; // Scaled bin extends low
        else if(bin_min_u <  bin_min_s && bin_max_u <= bin_max_s) overlap = bin_max_u - bin_min_s; // Scaled bin extends high
        else if(bin_min_u >  bin_max_s) break ; //Unscaled bin entirely larger than scaled bin, stop checking higher unscaled
        else overlap = -999; // No crossover
        
        if (overlap <= 0) continue;
        
        double new_events = events_per_bin[ibin_s][ibin_q2] * overlap / (scaled_bin_edges[ibin_s+1][ibin_q2] - scaled_bin_edges[ibin_s][ibin_q2]);
        new_events_per_bin[ibin_u][ibin_q2] += new_events;
      }
    }
    // Move events back into main array.
    for (int i = 0 ; i < n_ereco_bins ; i++) {
      events_per_bin[i][ibin_q2] = new_events_per_bin[i][ibin_q2];
    }
    events_per_bin[n_ereco_bins-1][ibin_q2] += new_events_per_bin[n_ereco_bins][ibin_q2]; //Put events from overflow bin back into top bin
    
#ifndef NDEBUG
  double new_total_events = 0;
  for (int i = 0 ; i < n_ereco_bins ; i++) {
    new_total_events+=events_per_bin[i][ibin_q2];
  }
  if ( fabs((new_total_events/total_events)-1) > conv::kASmallNum ) {
    cerr<<"Util/ApplyEnergyScaleSystematic2D - failed normalisation check for sample "<<util::AsString(s)<<endl;
    cerr<<"Q2 bin = " << ibin_q2 << endl;
    cerr<<"Events: "<<total_events<<" ==> "<<new_total_events<<"\t"<<new_total_events/total_events
        <<"\t"<<fabs((new_total_events/total_events)-1)<<"\t at scale "<<escale_tweak<<endl;
    assert(!"Broken");
  }
#endif //NDEBUG

  } // loop over q2 
  
  return;
}
//










//___________________________________________________________________________

// A few helper functions for use in UseEventToFillInTemplate
// that make swapping mode definitions a little easier

bool valor::util::IsCCQE(int neut_reaction_mode) 
{
  return neut_reaction_mode ==  1;
}
//___________________________________________________________________________
bool valor::util::IsCC1pi(int neut_reaction_mode) 
{
  bool is_cc1pi = (neut_reaction_mode == 11 || 
                   neut_reaction_mode == 12 || 
                   neut_reaction_mode == 13);
  return is_cc1pi;
}
//___________________________________________________________________________
bool valor::util::IsCCcoh(int neut_reaction_mode) 
{
  return neut_reaction_mode ==  16;
}
//___________________________________________________________________________
bool valor::util::Is2p2h(int neut_reaction_mode) 
{
  return neut_reaction_mode ==  2;
}
//___________________________________________________________________________
bool valor::util::IsCCoth(int neut_reaction_mode) 
//this is a catch all for everything CC that is not QE, 1pi, or coh
// needs to be updated if another subcategory is added
{
  bool is_ccoth = ( neut_reaction_mode == 17 ||
                    neut_reaction_mode == 21 ||
                    neut_reaction_mode == 22 || 
                    neut_reaction_mode == 23 || 
                    neut_reaction_mode == 26 );
  return is_ccoth;
}
//___________________________________________________________________________
//bool valor::util::IsCCoth(int neut_reaction_mode) 
//{
//return IsCCoth(neut_reaction_mode);
//}
//___________________________________________________________________________
//bool valor::util::IsCCoth(int neut_reaction_mode) 
//{
//return IsCCoth(neut_reaction_mode);
//}
//___________________________________________________________________________
bool valor::util::IsNC1pi0(int neut_reaction_mode) 
{
  bool is_cc1pi0 = (neut_reaction_mode == 31 ||
                    neut_reaction_mode == 32);
  return is_cc1pi0;
}
//___________________________________________________________________________
bool valor::util::IsNC1Gamma(int neut_reaction_mode)
{
    bool is_nc1gamma= (neut_reaction_mode == 38 ||
                       neut_reaction_mode == 39);
    return is_nc1gamma;
}
//___________________________________________________________________________
bool valor::util::IsNC1piC(int neut_reaction_mode) 
{
  bool is_cc1pi0 = (neut_reaction_mode == 33 ||
                    neut_reaction_mode == 34);
  return is_cc1pi0;
}
//___________________________________________________________________________
bool valor::util::IsNCcoh(int neut_reaction_mode) 
{
  return neut_reaction_mode == 36;
}
//___________________________________________________________________________
bool valor::util::IsNCoth(int neut_reaction_mode) 
//this is a catch all for everything NC that is not 1pi0, 1piC or coh
// needs to be updated if another subcategory is added
{
    bool is_ncoth = ( //  neut_reaction_mode == 38 || now in NC 1 gamma!
                     //   neut_reaction_mode == 39 || now in NC 1 gamma!
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
bool valor::util::UseEventToFillInTemplate(
  SuperkMode_t m, int nu_pdg_code, int neut_reaction_mode)
{
  bool is_numu    = (nu_pdg_code == kPdgNuMu);
  bool is_nue     = (nu_pdg_code == kPdgNuE);
  bool is_numubar = (nu_pdg_code == kPdgNuMuBar);
  bool is_nuebar  = (nu_pdg_code == kPdgNuEBar);
  
  int abs_mode = TMath::Abs(neut_reaction_mode);

  switch (m) {
    // Numu
    case kSKmNuMuCCQE         :  return ( IsCCQE    (abs_mode) && is_numu    );
    case kSKmNuMuCC1pi        :  return ( IsCC1pi   (abs_mode) && is_numu    );
    case kSKmNuMuCCcoh        :  return ( IsCCcoh   (abs_mode) && is_numu    );
    case kSKmNuMuCC2p2h        :  return ( Is2p2h     (abs_mode) && is_numu    );
    case kSKmNuMuCCoth        :  return ( IsCCoth   (abs_mode) && is_numu    );
    case kSKmNuMuNC1piC       :  return ( IsNC1piC  (abs_mode) && is_numu    );
    case kSKmNuMuNC1pi0       :  return ( IsNC1pi0  (abs_mode) && is_numu    );
    case kSKmNuMuNCcoh        :  return ( IsNCcoh   (abs_mode) && is_numu    );
    case kSKmNuMuNCoth        :  return ( IsNCoth   (abs_mode) && is_numu    );
      case kSKmNuMuNC1Gamma     :  return ( IsNC1Gamma(abs_mode) && is_numu    );
    
    case kSKmNuMuBarCCQE      :  return ( IsCCQE    (abs_mode) && is_numubar );
    case kSKmNuMuBarCC1pi     :  return ( IsCC1pi   (abs_mode) && is_numubar );
    case kSKmNuMuBarCCcoh     :  return ( IsCCcoh   (abs_mode) && is_numubar );
    case kSKmNuMuBarCC2p2h     :  return ( Is2p2h     (abs_mode) && is_numubar );
    case kSKmNuMuBarCCoth     :  return ( IsCCoth   (abs_mode) && is_numubar );
    case kSKmNuMuBarNC1piC    :  return ( IsNC1piC  (abs_mode) && is_numubar );
    case kSKmNuMuBarNC1pi0    :  return ( IsNC1pi0  (abs_mode) && is_numubar );
    case kSKmNuMuBarNCcoh     :  return ( IsNCcoh   (abs_mode) && is_numubar );
    case kSKmNuMuBarNCoth     :  return ( IsNCoth   (abs_mode) && is_numubar );
      case kSKmNuMuBarNC1Gamma  :  return ( IsNC1Gamma(abs_mode) && is_numubar );

    
    // Osc Numu
    case kSKmOscNuMuCCQE      :  return ( IsCCQE    (abs_mode) && is_numu    );
    case kSKmOscNuMuCC1pi     :  return ( IsCC1pi   (abs_mode) && is_numu    );
    case kSKmOscNuMuCCcoh     :  return ( IsCCcoh   (abs_mode) && is_numu    );
    case kSKmOscNuMuCC2p2h     :  return ( Is2p2h     (abs_mode) && is_numu    );
    case kSKmOscNuMuCCoth     :  return ( IsCCoth   (abs_mode) && is_numu    );
    
    case kSKmOscNuMuBarCCQE   :  return ( IsCCQE    (abs_mode) && is_numubar );
    case kSKmOscNuMuBarCC1pi  :  return ( IsCC1pi   (abs_mode) && is_numubar );
    case kSKmOscNuMuBarCCcoh  :  return ( IsCCcoh   (abs_mode) && is_numubar );
    case kSKmOscNuMuBarCC2p2h  :  return ( Is2p2h     (abs_mode) && is_numubar );
    case kSKmOscNuMuBarCCoth  :  return ( IsCCoth   (abs_mode) && is_numubar );
    
    //NuE
    case kSKmNuECCQE          :  return ( IsCCQE    (abs_mode) && is_nue     );
    case kSKmNuECC1pi         :  return ( IsCC1pi   (abs_mode) && is_nue     );
    case kSKmNuECCcoh         :  return ( IsCCcoh   (abs_mode) && is_nue     );
    case kSKmNuECC2p2h         :  return ( Is2p2h     (abs_mode) && is_nue     );
    case kSKmNuECCoth         :  return ( IsCCoth   (abs_mode) && is_nue     );
    case kSKmNuENC1piC        :  return ( IsNC1piC  (abs_mode) && is_nue     );
    case kSKmNuENC1pi0        :  return ( IsNC1pi0  (abs_mode) && is_nue     );
    case kSKmNuENCcoh         :  return ( IsNCcoh   (abs_mode) && is_nue     );
    case kSKmNuENCoth         :  return ( IsNCoth   (abs_mode) && is_nue     );
      case kSKmNuENC1Gamma      :  return ( IsNC1Gamma(abs_mode) && is_nue     );
          
    case kSKmNuEBarCCQE       :  return ( IsCCQE    (abs_mode) && is_nuebar  );
    case kSKmNuEBarCC1pi      :  return ( IsCC1pi   (abs_mode) && is_nuebar  );
    case kSKmNuEBarCCcoh      :  return ( IsCCcoh   (abs_mode) && is_nuebar  );
    case kSKmNuEBarCC2p2h      :  return ( Is2p2h     (abs_mode) && is_nuebar  );
    case kSKmNuEBarCCoth      :  return ( IsCCoth   (abs_mode) && is_nuebar  );
    case kSKmNuEBarNC1piC     :  return ( IsNC1piC  (abs_mode) && is_nuebar  );
    case kSKmNuEBarNC1pi0     :  return ( IsNC1pi0  (abs_mode) && is_nuebar  );
    case kSKmNuEBarNCcoh      :  return ( IsNCcoh   (abs_mode) && is_nuebar  );
    case kSKmNuEBarNCoth      :  return ( IsNCoth   (abs_mode) && is_nuebar  );
      case kSKmNuEBarNC1Gamma   :  return ( IsNC1Gamma(abs_mode) && is_nuebar  );

    
    // Osc NuE
    case kSKmOscNuECCQE       :  return ( IsCCQE    (abs_mode) && is_nue     );
    case kSKmOscNuECC1pi      :  return ( IsCC1pi   (abs_mode) && is_nue     );
    case kSKmOscNuECCcoh      :  return ( IsCCcoh   (abs_mode) && is_nue     );
    case kSKmOscNuECC2p2h      :  return ( Is2p2h     (abs_mode) && is_nue     );
    case kSKmOscNuECCoth      :  return ( IsCCoth   (abs_mode) && is_nue     );
                                                                             
    case kSKmOscNuEBarCCQE    :  return ( IsCCQE    (abs_mode) && is_nuebar  );
    case kSKmOscNuEBarCC1pi   :  return ( IsCC1pi   (abs_mode) && is_nuebar  );
    case kSKmOscNuEBarCCcoh   :  return ( IsCCcoh   (abs_mode) && is_nuebar  );
    case kSKmOscNuEBarCC2p2h   :  return ( Is2p2h     (abs_mode) && is_nuebar  );
    case kSKmOscNuEBarCCoth   :  return ( IsCCoth   (abs_mode) && is_nuebar  );
    
    
    case kSKmUndefined: break;
  }
  
  std::cerr<<"util/UseEventToFillInTemplate >> Unhandled mode: "<<AsString(m)<<" ("<<m<<")"<<std::endl;
  return false;
}
//___________________________________________________________________________
valor::util::ArrIdx * valor::util::ArrIdx::fInstance = 0;
//____________________________________________________________________________
valor::util::ArrIdx::ArrIdx()
{
  this->Init();
  fInstance = 0;
}
//____________________________________________________________________________
valor::util::ArrIdx * valor::util::ArrIdx::Instance()
{
  if(fInstance == 0) {
    static valor::util::ArrIdx::Cleaner cleaner;
    cleaner.DummyMethodAndSilentCompiler();
    fInstance = new valor::util::ArrIdx;
  }
  return fInstance;
}
//____________________________________________________________________________
void valor::util::ArrIdx::Init(void)
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
     if(s != (SuperkSample_t)is) {
       cerr << "SuperkSample " << (int)s << " " << util::AsString(s)
            << " is at the wrong position " << is << " in kSKSamples"
            << endl;
       exit(-1);
     }
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
              << p << "\t" 
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
              << p << "\t" 
              << util::AsString(p) 
              << " is at slot: " << fRHCXSecSystIdx[p] << endl;
#endif
  }

  //BANFF
  for(int i = 0; i < knppmax; i++) {
    fBANFFSystIdx[i] = -1;
  }
  for(int ip = 0; ip < kNumOfParamsInBANFFfit; ip++) {
    PhysParamId_t p = kListOfBANFFParams[ip];
    if(p == kPhypUndefined) continue;
    fBANFFSystIdx[p] = ip;
#ifdef __show_arr_idx_mesg_debug__
     cout << "Util/ArrIdx >> BANFF: " 
          << p << "\t" 
          << util::AsString(p)
          << " is at slot: " << fBANFFSystIdx [p] << endl;
#endif
  }
  for(int i = 0; i < knppmax; i++) {
    fRHCBANFFSystIdx[i] = -1;
  }
  for(int ip = 0; ip < kNumOfRHCParamsInBANFFfit; ip++) {
    PhysParamId_t p = kListOfRHCBANFFParams[ip];
    if(p == kPhypUndefined) continue;
    fRHCBANFFSystIdx[p] = ip;
#ifdef __show_arr_idx_mesg_debug__
    cout << "Util/ArrIdx >> RHC BANFF: "
         << p << "\t" 
         << util::AsString(p)
         << " is at slot: " << fRHCBANFFSystIdx [p] << endl;
#endif
  }
  
  //BANFF xsec
  for(int i = 0; i < knppmax; i++) {
    fBANFFXSecSystIdx[i] = -1;
  }
  for(int ip = 0; ip < kNumOfParamsInBANFFfit-kNumOfFluxParamsInBANFFfit; ip++) {
    PhysParamId_t p = kListOfBANFFParams[ip+kNumOfFluxParamsInBANFFfit];
    if(p == kPhypUndefined) continue;
    fBANFFXSecSystIdx[p] = ip;
#ifdef __show_arr_idx_mesg_debug__
     cout << "Util/ArrIdx >> BANFF XSec: " 
          << p << "\t" 
          << util::AsString(p)
          << " is at slot: " << fBANFFXSecSystIdx [p] << endl;
#endif
  }
  
  //BANFF xsec
  for(int i = 0; i < knppmax; i++) {
    fRHCBANFFXSecSystIdx[i] = -1;
  }
  for(int ip = 0; ip < kNumOfParamsInBANFFfit-kNumOfFluxParamsInBANFFfit; ip++) {
    PhysParamId_t p = kListOfRHCBANFFParams[ip+kNumOfFluxParamsInBANFFfit];
    if(p == kPhypUndefined) continue;
    fRHCBANFFXSecSystIdx[p] = ip;
#ifdef __show_arr_idx_mesg_debug__
     cout << "Util/ArrIdx >> RHC BANFF XSec: " 
          << p << "\t" 
          << util::AsString(p)
          << " is at slot: " << fRHCBANFFXSecSystIdx [p] << endl;
#endif
  }

//FHC + RHC BANFF
    
    //BANFF
    for(int i = 0; i < knppmax; i++) {
        fFHCpRHCBANFFSystIdx[i] = -1;
    }
    for(int ip = 0; ip < kNumOfFHCpRHCParamsInBANFFfit; ip++) {
        PhysParamId_t p = kListOfFHCpRHCBANFFParams[ip];
        if(p == kPhypUndefined) continue;
        fFHCpRHCBANFFSystIdx[p] = ip;
#ifdef __show_arr_idx_mesg_debug__
        cout << "Util/ArrIdx >> BANFF: "
        << p << "\t"
        << util::AsString(p)
        << " is at slot: " << fFHCpRHCBANFFSystIdx [p] << endl;
#endif
    }

    


  //SKDET + FSI
  // fhc
  for(int i = 0; i < knppmax; i++) {
    fSKDetFSISystIdx[i] = -1;
  }
  for(int ip = 0; ip < kNumOfSKDetFSIParams; ip++) {


    // #if FSI_BINS != 3
    //     PhysParamId_t p = PhysParamId_t ((int) kNuisSKDetFSI000 + ip);
    // #else
    //     PhysParamId_t p = PhysParamId_t ((int) kNuisSKDetFSI006 + ip);
    // #endif

    // NEW SKDETFSI (use kNuisSKDetFSI_first)
    PhysParamId_t p = PhysParamId_t ((int) kNuisSKDetFSI_first + ip);
    
    if(p == kPhypUndefined) continue;
    fSKDetFSISystIdx[p] = ip;
#ifdef __show_arr_idx_mesg_debug__
     cout << "Util/ArrIdx >> SKDet+FSI: "
          << p << "\t"  
          << util::AsString(p)
          << " is at slot: " << fSKDetFSISystIdx [p] << endl;
#endif
  }
  // rhc
  for(int i = 0; i < knppmax; i++) {
    fRHCSKDetFSISystIdx[i] = -1;
  }
  for(int ip = 0; ip < kNumOfRHCSKDetFSIParams; ip++) {

    // #if FSI_BINS != 3
    //     PhysParamId_t p = PhysParamId_t ((int) kRHCNuisSKDetFSI000 + ip);
    // #else
    //     PhysParamId_t p = PhysParamId_t ((int) kRHCNuisSKDetFSI006 + ip);
    // #endif

    // NEW SKDETFSI (use kRHCNuisSKDetFSI_first)
    PhysParamId_t p = PhysParamId_t ((int) kRHCNuisSKDetFSI_first + ip);

    if(p == kPhypUndefined) continue;
    fRHCSKDetFSISystIdx[p] = ip;
#ifdef __show_arr_idx_mesg_debug__
    cout << "Util/ArrIdx >> RHC SKDet+FSI: "
         << p << "\t" 
         << util::AsString(p)
         << " is at slot: " << fRHCSKDetFSISystIdx [p] << endl;
#endif
  }




  // NEW FHC+RHC 
  
  // fhc+rhc
  for(int i = 0; i < knppmax; i++) {
    fFHCpRHCSKDetFSISystIdx[i] = -1;
  }
  for(int ip = 0; ip < kNumOfSKDetFSIParams; ip++) {
    PhysParamId_t p = PhysParamId_t ((int) kNuisSKDetFSI_first + ip);
    if(p == kPhypUndefined) continue;
    fFHCpRHCSKDetFSISystIdx[p] = ip;
#ifdef __show_arr_idx_mesg_debug__
    cout << "Util/ArrIdx >> FHC+RHC SKDet+FSI: "
         << p << "\t"  
         << util::AsString(p)
          << " is at slot: " << fFHCpRHCSKDetFSISystIdx [p] << endl;
#endif
  }
  
  for(int ip = 0; ip < kNumOfRHCSKDetFSIParams; ip++) {
    PhysParamId_t p = PhysParamId_t ((int) kRHCNuisSKDetFSI_first + ip);
    if(p == kPhypUndefined) continue;
    fFHCpRHCSKDetFSISystIdx[p] = ip+kNumOfSKDetFSIParams; // shifted because after the FHC params // NEW FHC+RHC JOINT (NEW)
#ifdef __show_arr_idx_mesg_debug__
    cout << "Util/ArrIdx >> FHC+RHC SKDet+FSI: "
         << p << "\t" 
         << util::AsString(p)
         << " is at slot: " << fFHCpRHCSKDetFSISystIdx [p] << endl;
#endif
  }
  //




  

  /// NEW MULTIRING 2

  // multiring
  for(int i = 0; i < knppmax; i++) {
    fMultiRingSKDetFSISystIdx[i] = -1;
  }
  for(int ip = 0; ip < kNumOfSKDetFSIParamsMultiRing; ip++) {

    // #if FSI_BINS != 3
    //     PhysParamId_t p = PhysParamId_t ((int) kMultiRingNuisSKDetFSI000 + ip);
    // #else
    //     PhysParamId_t p = PhysParamId_t ((int) kMultiRingNuisSKDetFSI006 + ip);
    // #endif
    
    // NEW SKDETFSI
    PhysParamId_t p = (PhysParamId_t)kListOfMultiRingSKDetFSIParams[ip];

    if(p == kPhypUndefined) continue;
    fMultiRingSKDetFSISystIdx[p] = ip;
#ifdef __show_arr_idx_mesg_debug__
    cout << "Util/ArrIdx >> MultiRing SKDet+FSI: "
         << util::AsString(p)
         << " is at slot: " << fMultiRingSKDetFSISystIdx [p] << endl;
#endif
  }

  //FHC+RHC+MR
  for(int i = 0; i < knppmax; i++) {
    fFHCpRHCpMRSKDetFSISystIdx[i] = -1;
  }
  for(int ip = 0; ip < kNumOfSKDetFSIParams; ip++) {
    PhysParamId_t p = PhysParamId_t ((int) kNuisSKDetFSI_first + ip);
    if(p == kPhypUndefined) continue;
    fFHCpRHCpMRSKDetFSISystIdx[p] = ip;
#ifdef __show_arr_idx_mesg_debug__
    cout << "Util/ArrIdx >> FHC+RHC+MR SKDet+FSI: "
         << p << "\t"  
         << util::AsString(p)
         << " is at slot: " << fFHCpRHCpMRSKDetFSISystIdx [p] << endl;
#endif
  }
  
  for(int ip = 0; ip < kNumOfRHCSKDetFSIParams; ip++) {
    PhysParamId_t p = PhysParamId_t ((int) kRHCNuisSKDetFSI_first + ip);
    if(p == kPhypUndefined) continue;
    fFHCpRHCpMRSKDetFSISystIdx[p] = ip+kNumOfSKDetFSIParams; // shifted because after the FHC params // NEW FHC+RHC+MR JOINT (NEW)
#ifdef __show_arr_idx_mesg_debug__
    cout << "Util/ArrIdx >> FHC+RHC+MR SKDet+FSI: "
         << p << "\t" 
         << util::AsString(p)
         << " is at slot: " << fFHCpRHCpMRSKDetFSISystIdx [p] << endl;
#endif
  }
  for(int ip = 0; ip < kNumOfSKDetFSIParamsMultiRing; ip++) {
     
    PhysParamId_t p = (PhysParamId_t)kListOfMultiRingSKDetFSIParams[ip];

    if(p == kPhypUndefined) continue;
    fFHCpRHCpMRSKDetFSISystIdx[p] = ip+kNumOfSKDetFSIParams+kNumOfRHCSKDetFSIParams;
#ifdef __show_arr_idx_mesg_debug__
    cout << "Util/ArrIdx >>  FHC+RHC+MR SKDet+FSI: "
         << util::AsString(p)
         << " is at slot: " << fMultiRingSKDetFSISystIdx [p] << endl;
#endif
  } 
  // Escale in now included in the FSI matrix
  //  PhysParamId_t p = kNuisSuperkEnergyScale;
  fFHCpRHCpMRSKDetFSISystIdx[(PhysParamId_t)kNuisSuperkEnergyScale] = kNumOfSKDetFSIParams+kNumOfRHCSKDetFSIParams+kNumOfSKDetFSIParamsMultiRing;   
    
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
valor::conv::PhysParamId_t valor::util::GetBANFFParam(int FHCpRHC_param_no)
{
  assert(FHCpRHC_param_no >= 0 && FHCpRHC_param_no < kNumOfFHCpRHCParamsInBANFFfit);
  return conv::kListOfFHCpRHCBANFFParams[FHCpRHC_param_no];
}
//____________________________________________________________________________
int valor::util::GetSystNSigma()
{
  int nsig = SYST_NSIGMA;
  return nsig;
}
//____________________________________________________________________________
  










