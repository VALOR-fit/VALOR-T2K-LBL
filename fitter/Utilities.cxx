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

#include "Utilities/Rnd.h"

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
  if(sh == "tun")        return kSKMCtTun;
  if(sh == "sys")        return kSKMCtSys;
  if(sh == "fin")        return kSKMCtFin;

  return kSKMCtUndefined;
}
//____________________________________________________________________________
void valor::util::PrintCompact(
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
         PhysParamId_t p = kListOfBANFFParams[ip];
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
void valor::util::PrintCompact(
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
         PhysParamId_t p = kListOfBANFFParams[ip];
         stream << setw(42) << util::AsString(p);
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
void valor::util::PrintCompact(
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
void valor::util::PrintCompactTex(
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
bool valor::util::IsCCoth(int neut_reaction_mode) 
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
bool valor::util::IsNCoth_mu(int neut_reaction_mode) 
{
  bool is_ncoth_mu = ( IsNCoth (neut_reaction_mode) ||
		       IsNCcoh (neut_reaction_mode) ||
		       IsNC1pi0(neut_reaction_mode) );
  return is_ncoth_mu;
}
//___________________________________________________________________________
bool valor::util::IsNCoth_e(int neut_reaction_mode) 
{
  bool is_ncoth_e = ( IsNCoth (neut_reaction_mode) ||
		      IsNC1piC(neut_reaction_mode) );
  return is_ncoth_e;
}
//___________________________________________________________________________
bool valor::util::UseEventToFillInTemplate(
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
    PhysParamId_t p = kListOfBANFFParams[ip];
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

  //SKDET + FSI
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
