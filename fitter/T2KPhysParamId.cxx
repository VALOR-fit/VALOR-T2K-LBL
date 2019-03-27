//________________________________________________________________________________________
/*

  Copyright (c) 2010-2014, The VALOR group. All rights reserved.

  ChangeLog:

  21-7-2014 Davide
  Added samples for multiring forward horn current

*/
//________________________________________________________________________________________

#include <string>
#include <map>
#include <algorithm>
#include <iostream>

#include <TROOT.h>

#include "T2KPhysParamId.h"
#include "Utilities.h"
#include "Enumerations.h"
#include "valor/Conventions/NameRegistry.h"
#include "valor/Messenger/Messenger.h"
//~ #include "NuOscMixParam.h"


using namespace valor;
using namespace valor::conv;


#if ( defined(__USE_NIWG2015__) || defined (__USE_NIWG2017__) )
#ifndef __USE_ALLFLUX__
#define __USE_ALLFLUX__
#endif
#endif

#if !(BANFF_BINS == 2)
#ifndef __USE_ALLFLUX__
#define __USE_ALLFLUX__
#endif
#endif


//________________________________________________________________________________________
string valor::util::AsString(PhysParamId_t p)
{
  if (!util::IsSystematic(p)) return util::RawParamAsString(p);
  switch(p) {

    /*
     * 2011-analysis nuisance parameters
     */

    //***************************** FORWARD HORN CURRENT **********************************

    //
    // SK detector
    //
  case ( kNuisSuperkEnergyScale     ) : return "SK energy scale";                 break;

    //
    // BANFF flux and cross-section params
    //

  case ( kNuisBANFF00 ) :
  case ( kNuisBANFF01 ) :
  case ( kNuisBANFF02 ) :
  case ( kNuisBANFF03 ) :
  case ( kNuisBANFF04 ) :
  case ( kNuisBANFF05 ) :
  case ( kNuisBANFF06 ) :
  case ( kNuisBANFF07 ) :
  case ( kNuisBANFF08 ) :
  case ( kNuisBANFF09 ) :
  case ( kNuisBANFF10 ) :
  case ( kNuisBANFF11 ) :
  case ( kNuisBANFF12 ) :
  case ( kNuisBANFF13 ) :
  case ( kNuisBANFF14 ) :
  case ( kNuisBANFF15 ) :
#ifdef __USE_ALLFLUX__
  case ( kNuisBANFF16 ) :
  case ( kNuisBANFF17 ) :
  case ( kNuisBANFF18 ) :
  case ( kNuisBANFF19 ) :
  case ( kNuisBANFF20 ) :
  case ( kNuisBANFF21 ) :
  case ( kNuisBANFF22 ) :
  case ( kNuisBANFF23 ) :
  case ( kNuisBANFF24 ) :
#endif
    {
      int ip = (int)p - (int)kNuisBANFF00;
      return std::string(Form("BANFF %2d; SK %s flux; %.1f - %.1f GeV",
			      ip, util::Pdg2String(kNuPdgCodeOfFluxParamInBANFFfit[ip]).c_str(),
			      kEtrueRangeOfFluxParamInBANFFfit[ip][0],
			      kEtrueRangeOfFluxParamInBANFFfit[ip][1]));
      break;
    }
#ifndef ONEOVERMAQE
  case ( kNuisBANFFMAQE    ) : return "BANFF; Ma QE";                           break;
#else
  case ( kNuisBANFFMAQE    ) : return "BANFF; 1 / (Ma QE)^2";                   break;
#endif
  case ( kNuisBANFFMARES   ) : return "BANFF; Ma RES";                          break;
  case ( kNuisBANFFCCQEE0  ) : return "BANFF; Norm; CC QE; 0.0 - 1.5 GeV";      break;
#if BANFF_BINS != 3
  case ( kNuisBANFFCCQEE1  ) : return "BANFF; Norm; CC QE; 1.5 - 3.5 GeV";      break;
  case ( kNuisBANFFCCQEE2  ) : return "BANFF; Norm; CC QE; 3.5 - 30.0 GeV";     break;
#endif
  case ( kNuisBANFFCC1piE0 ) : return "BANFF; Norm; CC 1pi; 0.0 - 2.5 GeV";     break;
#if BANFF_BINS != 3
  case ( kNuisBANFFCC1piE1 ) : return "BANFF; Norm; CC 1pi; 2.5 - 30.0 GeV";    break;
#endif
#if (BANFF_BINS != 2)
  case ( kNuisBANFFNC1pi0  ) : return "BANFF; Norm; NC 1pi0";   break;
#endif

    // NEW XSEC
  case ( kNuisBANFFNorm2p2h       ) : return "BANFF; Norm; 2p2h";                    break;
  case ( kNuisBANFFCA5            ) : return "BANFF; CA5 RES";                      break;
  case ( kNuisBANFFBgRES          ) : return "BANFF; BgRES Isospin 1/2";            break;
    // NEW SPLIT ANA
#if (BANFF_BINS != 2)
  case ( kNuisBANFFSCCA           ) : return "BANFF; Second current class axial";   break;
  case ( kNuisBANFFSCCV           ) : return "BANFF; Second current class vector";  break;
#endif
  case ( kNuisBANFFFermiMom       ) : return "BANFF; Fermi Momentum";               break; 
  case ( kNuisBANFFShapeCCoth     ) : return "BANFF; Shape; CC Oth";                break; 
    //  case ( kNuisBANFFBindE          ) : return "BANFF; Binding energy";               break; 
  case ( kNuisBANFFNormCCcoh      ) : return "BANFF; Norm; CC Coh";                 break; 
  case ( kNuisBANFFNormNCcoh      ) : return "BANFF; Norm; NC Coh";                 break; 
  case ( kNuisBANFFNormNCoth      ) : return "BANFF; Norm; NC Oth";                 break; 
  case ( kNuisBANFFNormNuEToNuMu  ) : return "BANFF; Norm; NuE To NuMu";            break; 
  case ( kNuisBANFFNormNC1Gamma        ) : return "BANFF; Norm; NC 1 Gamma";          break; 
  case ( kNuisBANFFNormNuEBarToNuMuBar ) : return "BANFF; Norm; NuEBar To NuMuBar";   break; 
  case ( kNuisBANFFNorm2p2hBar         ) : return "BANFF; Norm; 2p2hBar";             break;
  case ( kNuisBANFFBeRPA_A       ) : return "BANFF; BeRPA A";            break;
  case ( kNuisBANFFBeRPA_B       ) : return "BANFF; BeRPA B";            break;
  case ( kNuisBANFFBeRPA_D       ) : return "BANFF; BeRPA D";            break;
  case ( kNuisBANFFBeRPA_E       ) : return "BANFF; BeRPA E";            break;
  case ( kNuisBANFFBeRPA_U       ) : return "BANFF; BeRPA U";            break;
  case ( kNuisBANFFShape2p2h     ) : return "BANFF; Shape; 2p2h nu";     break;
  case ( kNuisBANFFNorm2p2hCtoO  ) : return "BANFF; Norm;  2p2h C to O";  break;
  case ( kNuisBANFFEb  ) : return "BANFF; Eb";  break;
  case ( kNuisBANFFMultiReNorm  ) : return "BANFF; FHC MultiRe normalisation";  break;

    //

    //
    // XSec params not correlated with flux
    //

    // NEW XSEC
  case ( kNuisXSecRPA        ) : return "Tuning; Norm; RPA";               break;

  case ( kNuisXSecFermiMom   ) : return "Fermi Momentum";                  break;
  case ( kNuisXSecSpecFunc   ) : return "Spectral Func";                   break;
  case ( kNuisXSecWshape     ) : return "W-shape";                         break;
  case ( kNuisXSecShapeCCoth ) : return "Shape; CC Oth";                   break;
  case ( kNuisXSecPiLessDelta) : return "Pi-less delta decay";             break;
#ifndef __NUEONLY__
  case ( kNuisXSecBindE      ) : return "Binding energy";                  break;
#endif
  case ( kNuisXSecNormCCcoh  ) : return "Norm; CC Coh";                    break;
#ifndef __NUEONLY__
  case ( kNuisXSecNormNC1piC ) : return "Norm; NC 1piC";                   break;
#endif
#ifndef __NUMUONLY__
  case ( kNuisXSecNormNCcoh  ) : return "Norm; NC Coh";                    break;
#endif
  case ( kNuisXSecNormNCoth  ) : return "Norm; NC Oth";                    break;
  case ( kNuisXSecNormNuEToNuMu  ) : return "Norm; NuE To NuMu";           break;
#ifndef __NUEONLY__
  case ( kNuisXSecNormAntiNuToNu ) : return "Norm; AntiNu To Nu";          break;
#endif

    //
    // SKDET+FSI
    //
    //#ifdef SKFSITOGETHER
#ifndef SKCORRELATIONS
#if FSI_BINS != 3
  case ( kNuisSKDetFSI000 ) :
  case ( kNuisSKDetFSI001 ) :
  case ( kNuisSKDetFSI002 ) :
  case ( kNuisSKDetFSI003 ) :
  case ( kNuisSKDetFSI004 ) :
  case ( kNuisSKDetFSI005 ) :

    // NEW SKDETFSI
#if defined(__USE_SKDETFSI2015__)
  case ( kNuisSKDetFSI006 ) :
  case ( kNuisSKDetFSI007 ) :
  case ( kNuisSKDetFSI008 ) :
  case ( kNuisSKDetFSI009 ) :
  case ( kNuisSKDetFSI010 ) :
  case ( kNuisSKDetFSI011 ) :
#endif

    {
      // NEW SKDETFSI
      //int ip = (int)p - (int)kNuisSKDetFSI000;
      int ip = (int)p - (int)kNuisSKDetFSI_first;

      string nutypes = "";
      switch(p) {

        // NEW SKDETFSI
#if defined(__USE_SKDETFSI2013__)
      case (kNuisSKDetFSI000) :
      case (kNuisSKDetFSI001) :
      case (kNuisSKDetFSI002) :
        nutypes = "numu/numubar CCQE (1Rmu)"; break;
      case (kNuisSKDetFSI003) :
        nutypes = "numu/numubar CCnQE (1Rmu)"; break;
      case (kNuisSKDetFSI004) :
        nutypes = "nue/nuebar/signue CC (1Rmu)"; break;
      case (kNuisSKDetFSI005) :
        nutypes = "all NC (1Rmu)"; break;
#elif defined(__USE_SKDETFSI2015__)
      case (kNuisSKDetFSI000) :
      case (kNuisSKDetFSI001) :
      case (kNuisSKDetFSI002) :
        nutypes = "numu/numubar CCQE (1Rmu)"; break;
      case (kNuisSKDetFSI003) :
      case (kNuisSKDetFSI004) :
      case (kNuisSKDetFSI005) :
        nutypes = "numu/numubar CCnQE (1Rmu)"; break;
      case (kNuisSKDetFSI006) :
      case (kNuisSKDetFSI007) :
      case (kNuisSKDetFSI008) :
        nutypes = "nue/nuebar/signue CC (1Rmu)"; break;
      case (kNuisSKDetFSI009) :
      case (kNuisSKDetFSI010) :
      case (kNuisSKDetFSI011) :
        nutypes = "all NC (1Rmu)"; break;
#endif

      default : nutypes = "--";
      }
      return std::string(Form("SKDet + FSI/SI %2d; Ereco range %.2f - %.2f GeV; %s",ip,
			      kSKDetFSIParamsErecoCut[ip][0],kSKDetFSIParamsErecoCut[ip][1],nutypes.c_str()));
      break;
    }//FSI for 1Rmu
#endif
#if FSI_BINS != 2

    // NEW SKDETFSI
#if defined(__USE_SKDETFSI2013__)
  case ( kNuisSKDetFSI006 ) :
  case ( kNuisSKDetFSI007 ) :
  case ( kNuisSKDetFSI008 ) :
  case ( kNuisSKDetFSI009 ) :
  case ( kNuisSKDetFSI010 ) :
  case ( kNuisSKDetFSI011 ) :
  case ( kNuisSKDetFSI012 ) :
  case ( kNuisSKDetFSI013 ) :
  case ( kNuisSKDetFSI014 ) :
  case ( kNuisSKDetFSI015 ) :
  case ( kNuisSKDetFSI016 ) :
  case ( kNuisSKDetFSI017 ) :
#elif defined(__USE_SKDETFSI2015__)
  case ( kNuisSKDetFSI012 ) :
  case ( kNuisSKDetFSI013 ) :
  case ( kNuisSKDetFSI014 ) :
  case ( kNuisSKDetFSI015 ) :
  case ( kNuisSKDetFSI016 ) :
  case ( kNuisSKDetFSI017 ) :
  case ( kNuisSKDetFSI018 ) :
  case ( kNuisSKDetFSI019 ) :
  case ( kNuisSKDetFSI020 ) :
  case ( kNuisSKDetFSI021 ) :
  case ( kNuisSKDetFSI022 ) :
  case ( kNuisSKDetFSI023 ) :
#endif

    {

      // #if FSI_BINS == 3
      //         #if defined(__USE_SKDETFSI2013__)
      //        int ip = (int)p - (int)kNuisSKDetFSI006;
      //       #elif defined(__USE_SKDETFSI2015__)
      //       int ip = (int)p - (int)kNuisSKDetFSI012;
      //       #endif
      // #elif FSI_BINS ==1
      //       int ip = (int)p - (int)kNuisSKDetFSI000;
      // #endif

      // NEW SKDETFSI
      int ip = (int)p - (int)kNuisSKDetFSI_first;

      string nutypes = "";
      switch(p) {

        // NEW SKDETFSI
#if defined(__USE_SKDETFSI2013__)
      case (kNuisSKDetFSI006) :
      case (kNuisSKDetFSI007) :
      case (kNuisSKDetFSI008) :
        nutypes = "oscillated nue CC (1Re)"; break;
      case (kNuisSKDetFSI009) :
      case (kNuisSKDetFSI010) :
      case (kNuisSKDetFSI011) :
        nutypes = "numu/numubar CC (1Re)"; break;
      case (kNuisSKDetFSI012) :
      case (kNuisSKDetFSI013) :
      case (kNuisSKDetFSI014) :
        nutypes = "nue/nuebar CC (1Re)"; break;
      case (kNuisSKDetFSI015) :
      case (kNuisSKDetFSI016) :
      case (kNuisSKDetFSI017) :
        nutypes = "all NC (1Re)"; break;
#elif defined(__USE_SKDETFSI2015__)
      case (kNuisSKDetFSI012) :
      case (kNuisSKDetFSI013) :
      case (kNuisSKDetFSI014) :
        nutypes = "oscillated nue CC (1Re)"; break;
      case (kNuisSKDetFSI015) :
      case (kNuisSKDetFSI016) :
      case (kNuisSKDetFSI017) :
        nutypes = "numu/numubar CC (1Re)"; break;
      case (kNuisSKDetFSI018) :
      case (kNuisSKDetFSI019) :
      case (kNuisSKDetFSI020) :
        nutypes = "nue/nuebar CC (1Re)"; break;
      case (kNuisSKDetFSI021) :
      case (kNuisSKDetFSI022) :
      case (kNuisSKDetFSI023) :
        nutypes = "all NC (1Re)"; break;
#endif

      default : nutypes = "--";
      }
      return std::string(Form("SKDet + FSI/SI %2d; Ereco range %.2f - %.2f GeV; %s",ip,
			      kSKDetFSIParamsErecoCut[ip][0],kSKDetFSIParamsErecoCut[ip][1],nutypes.c_str()));
      break;
    }
#endif
#else

    // OTHER FSI_BINS THAN 1,2,3

  case ( kNuisSKDetFSI000 ) :
  case ( kNuisSKDetFSI001 ) :
  case ( kNuisSKDetFSI002 ) :
  case ( kNuisSKDetFSI003 ) :
  case ( kNuisSKDetFSI004 ) :
  case ( kNuisSKDetFSI005 ) :
  case ( kNuisSKDetFSI006 ) :
  case ( kNuisSKDetFSI007 ) :
  case ( kNuisSKDetFSI008 ) :
  case ( kNuisSKDetFSI009 ) :
  case ( kNuisSKDetFSI010 ) :
  case ( kNuisSKDetFSI011 ) :
  case ( kNuisSKDetFSI012 ) :
  case ( kNuisSKDetFSI013 ) :
  case ( kNuisSKDetFSI014 ) :
  case ( kNuisSKDetFSI015 ) :
  case ( kNuisSKDetFSI016 ) :
  case ( kNuisSKDetFSI017 ) :
  case ( kNuisSKDetFSI018 ) :
  case ( kNuisSKDetFSI019 ) :
    {

      // NEW SKDETFSI
      //int ip = (int)p - (int)kNuisSKDetFSI000;
      int ip = (int)p - (int)kNuisSKDetFSI_first;

      string nutypes = "";
      switch(p) {
      case (kNuisSKDetFSI000) :
      case (kNuisSKDetFSI001) :
      case (kNuisSKDetFSI002) :
      case (kNuisSKDetFSI003) :
      case (kNuisSKDetFSI004) :
        nutypes = "signal nue CC (1Rmu)"; break;
      case (kNuisSKDetFSI005) :
      case (kNuisSKDetFSI006) :
      case (kNuisSKDetFSI007) :
      case (kNuisSKDetFSI008) :
      case (kNuisSKDetFSI009) :
        nutypes = "numu/numubar CC (1Rmu)"; break;
      case (kNuisSKDetFSI010) :
      case (kNuisSKDetFSI011) :
      case (kNuisSKDetFSI012) :
      case (kNuisSKDetFSI013) :
      case (kNuisSKDetFSI014) :
        nutypes = "nue/nuebar  CC (1Rmu)"; break;
      case (kNuisSKDetFSI015) :
      case (kNuisSKDetFSI016) :
      case (kNuisSKDetFSI017) :
      case (kNuisSKDetFSI018) :
      case (kNuisSKDetFSI019) :
        nutypes = "all NC (1Rmu)"; break;
      default : nutypes = "--";
      }
      return std::string(Form("SKDet + FSI/SI %2d; Ereco range %.2f - %.2f GeV; %s",ip,
			      kSKDetFSIParamsErecoCut[ip][0],kSKDetFSIParamsErecoCut[ip][1],nutypes.c_str()));
      break;
    }//FSI for 1Rmu
  case ( kNuisSKDetFSI020 ) :
  case ( kNuisSKDetFSI021 ) :
  case ( kNuisSKDetFSI022 ) :
  case ( kNuisSKDetFSI023 ) :
  case ( kNuisSKDetFSI024 ) :
  case ( kNuisSKDetFSI025 ) :
  case ( kNuisSKDetFSI026 ) :
  case ( kNuisSKDetFSI027 ) :
  case ( kNuisSKDetFSI028 ) :
  case ( kNuisSKDetFSI029 ) :
  case ( kNuisSKDetFSI030 ) :
  case ( kNuisSKDetFSI031 ) :
    {
      // NEW SKDETFSI
      //int ip = (int)p - (int)kNuisSKDetFSI000;
      int ip = (int)p - (int)kNuisSKDetFSI_first;

      string nutypes = "";
      switch(p) {
      case (kNuisSKDetFSI020) :
      case (kNuisSKDetFSI021) :
      case (kNuisSKDetFSI022) :
        nutypes = "oscillated nue CC (1Re)"; break;
      case (kNuisSKDetFSI023) :
      case (kNuisSKDetFSI024) :
      case (kNuisSKDetFSI025) :
        nutypes = "numu/numubar CC (1Re)"; break;
      case (kNuisSKDetFSI026) :
      case (kNuisSKDetFSI027) :
      case (kNuisSKDetFSI028) :
        nutypes = "nue/nuebar CC (1Re)"; break;
      case (kNuisSKDetFSI029) :
      case (kNuisSKDetFSI030) :
      case (kNuisSKDetFSI031) :
        nutypes = "all NC (1Re)"; break;
      default : nutypes = "--";
      }
      return std::string(Form("SKDet + FSI/SI %2d; Ereco range %.2f - %.2f GeV; %s",ip,
			      kSKDetFSIParamsErecoCut[ip][0],kSKDetFSIParamsErecoCut[ip][1],nutypes.c_str()));
      break;
    }
#endif // FSI_BINS
    //#endif //SKDETFSI


  case ( kNuisRHCToFHC ) : return "RHC to FHC normalisation"; break;

    //************************** REVERSE HORN CURRENT *****************************

  case ( kRHCNuisSuperkEnergyScale     ) : return "SK energy scale; RHC";                 break;

    //
    // BANFF flux and cross-section params
    //

  case ( kRHCNuisBANFF00 ) :
  case ( kRHCNuisBANFF01 ) :
  case ( kRHCNuisBANFF02 ) :
  case ( kRHCNuisBANFF03 ) :
  case ( kRHCNuisBANFF04 ) :
  case ( kRHCNuisBANFF05 ) :
  case ( kRHCNuisBANFF06 ) :
  case ( kRHCNuisBANFF07 ) :
  case ( kRHCNuisBANFF08 ) :
  case ( kRHCNuisBANFF09 ) :
  case ( kRHCNuisBANFF10 ) :
  case ( kRHCNuisBANFF11 ) :
  case ( kRHCNuisBANFF12 ) :
  case ( kRHCNuisBANFF13 ) :
  case ( kRHCNuisBANFF14 ) :
  case ( kRHCNuisBANFF15 ) :
#ifdef __USE_ALLFLUX__
  case ( kRHCNuisBANFF16 ) :
  case ( kRHCNuisBANFF17 ) :
  case ( kRHCNuisBANFF18 ) :
  case ( kRHCNuisBANFF19 ) :
  case ( kRHCNuisBANFF20 ) :
  case ( kRHCNuisBANFF21 ) :
  case ( kRHCNuisBANFF22 ) :
  case ( kRHCNuisBANFF23 ) :
  case ( kRHCNuisBANFF24 ) :
#endif
    {
      int ip = (int)p - (int)kRHCNuisBANFF00;
      return std::string(Form("BANFF %2d; SK %s flux; %.1f - %.1f GeV; RHC",
			      ip, util::Pdg2String(kNuPdgCodeOfRHCFluxParamInBANFFfit[ip]).c_str(),
			      kEtrueRangeOfRHCFluxParamInBANFFfit[ip][0],
			      kEtrueRangeOfRHCFluxParamInBANFFfit[ip][1]));
      break;
    }
#ifndef ONEOVERMAQE
  case ( kRHCNuisBANFFMAQE    ) : return "BANFF; Ma QE; RHC";                           break;
#else
  case ( kRHCNuisBANFFMAQE    ) : return "BANFF; 1 / (Ma QE)^2; RHC";                   break;
#endif
  case ( kRHCNuisBANFFMARES   ) : return "BANFF; Ma RES; RHC";                          break;
  case ( kRHCNuisBANFFCCQEE0  ) : return "BANFF; Norm; CC QE; 0.0 - 1.5 GeV; RHC";      break;
#if BANFF_BINS != 3
  case ( kRHCNuisBANFFCCQEE1  ) : return "BANFF; Norm; CC QE; 1.5 - 3.5 GeV; RHC";      break;
  case ( kRHCNuisBANFFCCQEE2  ) : return "BANFF; Norm; CC QE; 3.5 - 30.0 GeV; RHC";     break;
#endif
  case ( kRHCNuisBANFFCC1piE0 ) : return "BANFF; Norm; CC 1pi; 0.0 - 2.5 GeV; RHC";     break;
#if BANFF_BINS != 3
  case ( kRHCNuisBANFFCC1piE1 ) : return "BANFF; Norm; CC 1pi; 2.5 - 30.0 GeV; RHC";    break;
#endif
#if (BANFF_BINS != 2)
  case ( kRHCNuisBANFFNC1pi0  ) : return "BANFF; Norm; NC 1pi0; RHC";   break;
#endif

    // NEW XSEC
  case ( kRHCNuisBANFF2p2h           ) : return "BANFF; Norm; 2p2h; RHC";                  break;
  case ( kRHCNuisBANFFCA5           ) : return "BANFF; CA5 RES; RHC";                      break;
  case ( kRHCNuisBANFFBgRES         ) : return "BANFF; Norm; BgRES Isospin 1/2; RHC";      break;
    // NEW SPLIT ANA
#if (BANFF_BINS != 2)
  case ( kRHCNuisBANFFSCCA          ) : return "BANFF; Second current class axial; RHC";   break;
  case ( kRHCNuisBANFFSCCV          ) : return "BANFF; Second current class vector; RHC";  break;
#endif
  case ( kRHCNuisBANFFFermiMom      ) : return "BANFF; Fermi Momentum; RHC";               break; 
  case ( kRHCNuisBANFFShapeCCoth    ) : return "BANFF; Shape; CC Oth; RHC";                break; 
    // case ( kRHCNuisBANFFBindE         ) : return "BANFF; Binding energy; RHC";               break; 
  case ( kRHCNuisBANFFNormCCcoh     ) : return "BANFF; Norm; CC Coh; RHC";                 break; 
  case ( kRHCNuisBANFFNormNCcoh     ) : return "BANFF; Norm; NC Coh; RHC";                 break; 
  case ( kRHCNuisBANFFNormNCoth     ) : return "BANFF; Norm; NC Oth; RHC";                 break; 
  case ( kRHCNuisBANFFNormNuEToNuMu ) : return "BANFF; Norm; NuE To NuMu; RHC";            break; 
  case ( kRHCNuisBANFFNormNC1Gamma  ) : return "BANFF; Norm; NC 1 Gamma; RHC";             break; 
  case ( kRHCNuisBANFFNormNuEBarToNuMuBar ) : return "BANFF; Norm; NuEBar To NuMuBar; RHC";  break; 
  case ( kRHCNuisBANFF2p2hBar       ) : return "BANFF; Norm; 2p2hBar; RHC"; break;
  case ( kRHCNuisBANFFBeRPA_A       ) : return "BANFF; BeRPA A; RHC";            break;
  case ( kRHCNuisBANFFBeRPA_B       ) : return "BANFF; BeRPA B; RHC";            break;
  case ( kRHCNuisBANFFBeRPA_D       ) : return "BANFF; BeRPA D; RHC";            break;
  case ( kRHCNuisBANFFBeRPA_E       ) : return "BANFF; BeRPA E; RHC";            break;
  case ( kRHCNuisBANFFBeRPA_U       ) : return "BANFF; BeRPA U; RHC";            break;
  case ( kRHCNuisBANFFShape2p2h     ) : return "BANFF; 2p2h nu nubar shape; RHC"; break;
  case ( kRHCNuisBANFFNorm2p2hCtoO  ) : return "BANFF; 2p2h C to O; RHC";
  case ( kRHCNuisBANFFEb  ) : return "BANFF; Eb; RHC";
  case ( kRHCNuisBANFFMultiReNorm   ) : return "BANFF; MultiRe normalisation; RHC; UNUSED";

    //

    //
    // XSec params not correlated with flux
    //

    // NEW XSEC
  case ( kRHCNuisXSecRPA        ) : return "Tuning; RPA; RHC";                     break;

  case ( kRHCNuisXSecFermiMom   ) : return "Fermi Momentum; RHC";                  break;
  case ( kRHCNuisXSecSpecFunc   ) : return "Spectral Func; RHC";                   break;
  case ( kRHCNuisXSecWshape     ) : return "W-shape; RHC";                         break;
  case ( kRHCNuisXSecShapeCCoth ) : return "Shape; CC Oth; RHC";                   break;
  case ( kRHCNuisXSecPiLessDelta) : return "Pi-less delta decay; RHC";             break;
#ifndef __NUEONLY__
  case ( kRHCNuisXSecBindE      ) : return "Binding energy; RHC";                  break;
#endif
  case ( kRHCNuisXSecNormCCcoh  ) : return "Norm; CC Coh; RHC";                    break;
#ifndef __NUEONLY__
  case ( kRHCNuisXSecNormNC1piC ) : return "Norm; NC 1piC; RHC";                   break;
#endif
#ifndef __NUMUONLY__
  case ( kRHCNuisXSecNormNCcoh  ) : return "Norm; NC Coh; RHC";                    break;
#endif
  case ( kRHCNuisXSecNormNCoth  ) : return "Norm; NC Oth; RHC";                    break;
  case ( kRHCNuisXSecNormNuEToNuMu  ) : return "Norm; NuE To NuMu; RHC";           break;
#ifndef __NUEONLY__
  case ( kRHCNuisXSecNormAntiNuToNu ) : return "Norm; AntiNu To Nu; RHC";          break;
#endif

    //
    // SKDET+FSI
    //
#ifndef SKCORRELATIONS
#if FSI_BINS != 3
  case ( kRHCNuisSKDetFSI000 ) :
  case ( kRHCNuisSKDetFSI001 ) :
  case ( kRHCNuisSKDetFSI002 ) :
  case ( kRHCNuisSKDetFSI003 ) :
  case ( kRHCNuisSKDetFSI004 ) :
  case ( kRHCNuisSKDetFSI005 ) :

    // NEW SKDETFSI
#if defined(__USE_SKDETFSI2015__)
  case ( kRHCNuisSKDetFSI006 ) :
  case ( kRHCNuisSKDetFSI007 ) :
  case ( kRHCNuisSKDetFSI008 ) :
  case ( kRHCNuisSKDetFSI009 ) :
  case ( kRHCNuisSKDetFSI010 ) :
  case ( kRHCNuisSKDetFSI011 ) :
#endif

    {
      // NEW SKDETFSI
      //int ip = (int)p - (int)kRHCNuisSKDetFSI000;
      int ip = (int)p - (int)kRHCNuisSKDetFSI_first;

      string nutypes = "";
      switch(p) {

        // NEW SKDETFSI
#if defined(__USE_SKDETFSI2013__)
      case (kRHCNuisSKDetFSI000) :
      case (kRHCNuisSKDetFSI001) :
      case (kRHCNuisSKDetFSI002) :
        nutypes = "numu/numubar CCQE (1Rmu)"; break;
      case (kRHCNuisSKDetFSI003) :
        nutypes = "numu/numubar CCnQE (1Rmu)"; break;
      case (kRHCNuisSKDetFSI004) :
        nutypes = "nue/nuebar/signue CC (1Rmu)"; break;
      case (kRHCNuisSKDetFSI005) :
        nutypes = "all NC (1Rmu)"; break;
#elif defined(__USE_SKDETFSI2015__)
      case (kRHCNuisSKDetFSI000) :
      case (kRHCNuisSKDetFSI001) :
      case (kRHCNuisSKDetFSI002) :
        nutypes = "numu/numubar CCQE (1Rmu)"; break;
      case (kRHCNuisSKDetFSI003) :
      case (kRHCNuisSKDetFSI004) :
      case (kRHCNuisSKDetFSI005) :
        nutypes = "numu/numubar CCnQE (1Rmu)"; break;
      case (kRHCNuisSKDetFSI006) :
      case (kRHCNuisSKDetFSI007) :
      case (kRHCNuisSKDetFSI008) :
        nutypes = "nue/nuebar/signue CC (1Rmu)"; break;
      case (kRHCNuisSKDetFSI009) :
      case (kRHCNuisSKDetFSI010) :
      case (kRHCNuisSKDetFSI011) :
        nutypes = "all NC (1Rmu)"; break;
#endif

      default : nutypes = "--";
      }
      return std::string(Form("SKDet + FSI/SI %2d; Ereco range %.2f - %.2f GeV; %s; RHC",ip,
			      kRHCSKDetFSIParamsErecoCut[ip][0],kRHCSKDetFSIParamsErecoCut[ip][1],nutypes.c_str()));
      break;
    }//FSI for 1Rmu
#endif
#if FSI_BINS != 2

    // NEW SKDETFSI
#if defined(__USE_SKDETFSI2013__)
  case ( kRHCNuisSKDetFSI006 ) :
  case ( kRHCNuisSKDetFSI007 ) :
  case ( kRHCNuisSKDetFSI008 ) :
  case ( kRHCNuisSKDetFSI009 ) :
  case ( kRHCNuisSKDetFSI010 ) :
  case ( kRHCNuisSKDetFSI011 ) :
  case ( kRHCNuisSKDetFSI012 ) :
  case ( kRHCNuisSKDetFSI013 ) :
  case ( kRHCNuisSKDetFSI014 ) :
  case ( kRHCNuisSKDetFSI015 ) :
  case ( kRHCNuisSKDetFSI016 ) :
  case ( kRHCNuisSKDetFSI017 ) :
#elif defined(__USE_SKDETFSI2015__)
  case ( kRHCNuisSKDetFSI012 ) :
  case ( kRHCNuisSKDetFSI013 ) :
  case ( kRHCNuisSKDetFSI014 ) :
  case ( kRHCNuisSKDetFSI015 ) :
  case ( kRHCNuisSKDetFSI016 ) :
  case ( kRHCNuisSKDetFSI017 ) :
  case ( kRHCNuisSKDetFSI018 ) :
  case ( kRHCNuisSKDetFSI019 ) :
  case ( kRHCNuisSKDetFSI020 ) :
  case ( kRHCNuisSKDetFSI021 ) :
  case ( kRHCNuisSKDetFSI022 ) :
  case ( kRHCNuisSKDetFSI023 ) :
#endif

    {
      //#if FSI_BINS == 3
      // #if defined(__USE_SKDETFSI2013__)
      //       int ip = (int)p - (int)kRHCNuisSKDetFSI006;
      // #elif defined(__USE_SKDETFSI2015__)
      //       int ip = (int)p - (int)kRHCNuisSKDetFSI012;
      // #endif
      // #elif FSI_BINS ==1
      //       int ip = (int)p - (int)kRHCNuisSKDetFSI000;
      // #endif

      // NEW SKDETFSI
      int ip = (int)p - (int)kRHCNuisSKDetFSI_first;

      string nutypes = "";
      switch(p) {

        // NEW SKDETFSI
#if defined(__USE_SKDETFSI2013__)
      case (kRHCNuisSKDetFSI006) :
      case (kRHCNuisSKDetFSI007) :
      case (kRHCNuisSKDetFSI008) :
        nutypes = "oscillated nue CC (1Re)"; break;
      case (kRHCNuisSKDetFSI009) :
      case (kRHCNuisSKDetFSI010) :
      case (kRHCNuisSKDetFSI011) :
        nutypes = "numu/numubar CC (1Re)"; break;
      case (kRHCNuisSKDetFSI012) :
      case (kRHCNuisSKDetFSI013) :
      case (kRHCNuisSKDetFSI014) :
        nutypes = "nue/nuebar CC (1Re)"; break;
      case (kRHCNuisSKDetFSI015) :
      case (kRHCNuisSKDetFSI016) :
      case (kRHCNuisSKDetFSI017) :
        nutypes = "all NC (1Re)"; break;
#elif defined(__USE_SKDETFSI2015__)
      case (kRHCNuisSKDetFSI012) :
      case (kRHCNuisSKDetFSI013) :
      case (kRHCNuisSKDetFSI014) :
        nutypes = "oscillated nue CC (1Re)"; break;
      case (kRHCNuisSKDetFSI015) :
      case (kRHCNuisSKDetFSI016) :
      case (kRHCNuisSKDetFSI017) :
        nutypes = "numu/numubar CC (1Re)"; break;
      case (kRHCNuisSKDetFSI018) :
      case (kRHCNuisSKDetFSI019) :
      case (kRHCNuisSKDetFSI020) :
        nutypes = "nue/nuebar CC (1Re)"; break;
      case (kRHCNuisSKDetFSI021) :
      case (kRHCNuisSKDetFSI022) :
      case (kRHCNuisSKDetFSI023) :
        nutypes = "all NC (1Re)"; break;
#endif

      default : nutypes = "--";
      }
      return std::string(Form("SKDet + FSI/SI %2d; Ereco range %.2f - %.2f GeV; %s; RHC",ip,
			      kRHCSKDetFSIParamsErecoCut[ip][0],kRHCSKDetFSIParamsErecoCut[ip][1],nutypes.c_str()));
      break;
    }
#endif
#else

    // OTHER FSI_BINS THAN 1,2,3

  case ( kRHCNuisSKDetFSI000 ) :
  case ( kRHCNuisSKDetFSI001 ) :
  case ( kRHCNuisSKDetFSI002 ) :
  case ( kRHCNuisSKDetFSI003 ) :
  case ( kRHCNuisSKDetFSI004 ) :
  case ( kRHCNuisSKDetFSI005 ) :
  case ( kRHCNuisSKDetFSI006 ) :
  case ( kRHCNuisSKDetFSI007 ) :
  case ( kRHCNuisSKDetFSI008 ) :
  case ( kRHCNuisSKDetFSI009 ) :
  case ( kRHCNuisSKDetFSI010 ) :
  case ( kRHCNuisSKDetFSI011 ) :
  case ( kRHCNuisSKDetFSI012 ) :
  case ( kRHCNuisSKDetFSI013 ) :
  case ( kRHCNuisSKDetFSI014 ) :
  case ( kRHCNuisSKDetFSI015 ) :
  case ( kRHCNuisSKDetFSI016 ) :
  case ( kRHCNuisSKDetFSI017 ) :
  case ( kRHCNuisSKDetFSI018 ) :
  case ( kRHCNuisSKDetFSI019 ) :
    
      {
        // NEW SKDETFSI
        //int ip = (int)p - (int)kRHCNuisSKDetFSI000;
        int ip = (int)p - (int)kRHCNuisSKDetFSI_first;

        string nutypes = "";
        switch(p) {
        case (kRHCNuisSKDetFSI000) :
        case (kRHCNuisSKDetFSI001) :
        case (kRHCNuisSKDetFSI002) :
        case (kRHCNuisSKDetFSI003) :
        case (kRHCNuisSKDetFSI004) :
          nutypes = "signal nue CC (1Rmu)"; break;
        case (kRHCNuisSKDetFSI005) :
        case (kRHCNuisSKDetFSI006) :
        case (kRHCNuisSKDetFSI007) :
        case (kRHCNuisSKDetFSI008) :
        case (kRHCNuisSKDetFSI009) :
          nutypes = "numu/numubar CC (1Rmu)"; break;
        case (kRHCNuisSKDetFSI010) :
        case (kRHCNuisSKDetFSI011) :
        case (kRHCNuisSKDetFSI012) :
        case (kRHCNuisSKDetFSI013) :
        case (kRHCNuisSKDetFSI014) :
          nutypes = "nue/nuebar  CC (1Rmu)"; break;
        case (kRHCNuisSKDetFSI015) :
        case (kRHCNuisSKDetFSI016) :
        case (kRHCNuisSKDetFSI017) :
        case (kRHCNuisSKDetFSI018) :
        case (kRHCNuisSKDetFSI019) :
          nutypes = "all NC (1Rmu)"; break;
        default : nutypes = "--";
        }
        return std::string(Form("SKDet + FSI/SI %2d; Ereco range %.2f - %.2f GeV; %s; RHC",ip,
				kRHCSKDetFSIParamsErecoCut[ip][0],kRHCSKDetFSIParamsErecoCut[ip][1],nutypes.c_str()));
        break;
      }//FSI for 1Rmu
    case ( kRHCNuisSKDetFSI020 ) :
    case ( kRHCNuisSKDetFSI021 ) :
    case ( kRHCNuisSKDetFSI022 ) :
    case ( kRHCNuisSKDetFSI023 ) :
    case ( kRHCNuisSKDetFSI024 ) :
    case ( kRHCNuisSKDetFSI025 ) :
    case ( kRHCNuisSKDetFSI026 ) :
    case ( kRHCNuisSKDetFSI027 ) :
    case ( kRHCNuisSKDetFSI028 ) :
    case ( kRHCNuisSKDetFSI029 ) :
    case ( kRHCNuisSKDetFSI030 ) :
    case ( kRHCNuisSKDetFSI031 ) :
      {
        // NEW SKDETFSI
        //int ip = (int)p - (int)kRHCNuisSKDetFSI000;
        int ip = (int)p - (int)kRHCNuisSKDetFSI_first;

        string nutypes = "";
        switch(p) {
        case (kRHCNuisSKDetFSI020) :
        case (kRHCNuisSKDetFSI021) :
        case (kRHCNuisSKDetFSI022) :
          nutypes = "oscillated nue CC (1Re)"; break;
        case (kRHCNuisSKDetFSI023) :
        case (kRHCNuisSKDetFSI024) :
        case (kRHCNuisSKDetFSI025) :
          nutypes = "numu/numubar CC (1Re)"; break;
        case (kRHCNuisSKDetFSI026) :
        case (kRHCNuisSKDetFSI027) :
        case (kRHCNuisSKDetFSI028) :
          nutypes = "nue/nuebar CC (1Re)"; break;
        case (kRHCNuisSKDetFSI029) :
        case (kRHCNuisSKDetFSI030) :
        case (kRHCNuisSKDetFSI031) :
          nutypes = "all NC (1Re)"; break;
        default : nutypes = "--";
        }
        return std::string(Form("SKDet + FSI/SI %2d; Ereco range %.2f - %.2f GeV; %s; RHC",ip,
				kRHCSKDetFSIParamsErecoCut[ip][0],kRHCSKDetFSIParamsErecoCut[ip][1],nutypes.c_str()));
        break;
      }
#endif


      //************************** MULTIRING FORWARD HORN CURRENT *****************************

      /// NEW MULTIRING 3
    case ( kMultiRingNuisSuperkEnergyScale     ) : return "SK energy scale; MultiRing";  break;
      ///
      
      //// NEW MULTIRING 2

      //
      // SKDET+FSI
      //
#ifndef SKCORRELATIONS
#  if MR_FSI_BINS != 3
    case ( kMultiRingNuisSKDetFSI000 ) :
    case ( kMultiRingNuisSKDetFSI001 ) :
    case ( kMultiRingNuisSKDetFSI002 ) :
    case ( kMultiRingNuisSKDetFSI003 ) :
    case ( kMultiRingNuisSKDetFSI004 ) :
    case ( kMultiRingNuisSKDetFSI005 ) :


      {
        // NEW SKDETFSI
	int ip = ArrIdx::Instance()->MultiRingSKDetFSISystIdx(p);

        string nutypes = "";
        switch(p) {

          // NEW SKDETFSI
#    if defined(__USE_SKDETFSI2013__)
        case (kMultiRingNuisSKDetFSI000) :
        case (kMultiRingNuisSKDetFSI001) :
        case (kMultiRingNuisSKDetFSI002) :
          nutypes = "numu/numubar CCQE (MultiRmu)"; break;
        case (kMultiRingNuisSKDetFSI003) :
          nutypes = "numu/numubar CCnQE (MultiRmu)"; break;
        case (kMultiRingNuisSKDetFSI004) :
          nutypes = "nue/nuebar/signue CC (MultiRmu)"; break;
        case (kMultiRingNuisSKDetFSI005) :
          nutypes = "all NC (MultiRmu)"; break;
#    endif

        default : nutypes = "--";
        }
        return std::string(Form("SKDet + FSI/SI %2d; Ereco range %.2f - %.2f GeV; %s; MultiRing",ip,
				kMultiRingSKDetFSIParamsErecoCut[ip][0],kMultiRingSKDetFSIParamsErecoCut[ip][1],nutypes.c_str()));
        break;
      }//FSI for MultiRmu
#  endif
#  if MR_FSI_BINS != 2

      // NEW SKDETFSI
#    if defined(__USE_SKDETFSI2013__)
    case ( kMultiRingNuisSKDetFSI006 ) :
    case ( kMultiRingNuisSKDetFSI007 ) :
    case ( kMultiRingNuisSKDetFSI008 ) :
    case ( kMultiRingNuisSKDetFSI009 ) :
    case ( kMultiRingNuisSKDetFSI010 ) :
    case ( kMultiRingNuisSKDetFSI011 ) :
    case ( kMultiRingNuisSKDetFSI012 ) :
    case ( kMultiRingNuisSKDetFSI013 ) :
    case ( kMultiRingNuisSKDetFSI014 ) :
    case ( kMultiRingNuisSKDetFSI015 ) :
    case ( kMultiRingNuisSKDetFSI016 ) :
    case ( kMultiRingNuisSKDetFSI017 ) :
#    endif

      {
        // NEW SKDETFSI
	int ip = ArrIdx::Instance()->MultiRingSKDetFSISystIdx(p);

        string nutypes = "";
        switch(p) {

          // NEW SKDETFSI
#    if defined(__USE_SKDETFSI2013__)
        case (kMultiRingNuisSKDetFSI006) :
        case (kMultiRingNuisSKDetFSI007) :
        case (kMultiRingNuisSKDetFSI008) :
          nutypes = "oscillated nue CC (MultiRe)"; break;
        case (kMultiRingNuisSKDetFSI009) :
        case (kMultiRingNuisSKDetFSI010) :
        case (kMultiRingNuisSKDetFSI011) :
          nutypes = "numu/numubar CC (MultiRe)"; break;
        case (kMultiRingNuisSKDetFSI012) :
        case (kMultiRingNuisSKDetFSI013) :
        case (kMultiRingNuisSKDetFSI014) :
          nutypes = "nue/nuebar CC (MultiRe)"; break;
        case (kMultiRingNuisSKDetFSI015) :
        case (kMultiRingNuisSKDetFSI016) :
        case (kMultiRingNuisSKDetFSI017) :
          nutypes = "all NC (MultiRe)"; break;
#    endif

        default : nutypes = "--";
        }
        return std::string(Form("SKDet + FSI/SI %2d; Ereco range %.2f - %.2f GeV; %s; MultiRing",ip,
				kMultiRingSKDetFSIParamsErecoCut[ip][0],kMultiRingSKDetFSIParamsErecoCut[ip][1],nutypes.c_str()));
        break;
      }
#  endif
#  else

      // OTHER FSI_BINS THAN 1,2,3

    case ( kMultiRingNuisSKDetFSI000 ) :
    case ( kMultiRingNuisSKDetFSI001 ) :
    case ( kMultiRingNuisSKDetFSI002 ) :
    case ( kMultiRingNuisSKDetFSI003 ) :
    case ( kMultiRingNuisSKDetFSI004 ) :
    case ( kMultiRingNuisSKDetFSI005 ) :
    case ( kMultiRingNuisSKDetFSI006 ) :
    case ( kMultiRingNuisSKDetFSI007 ) :
    case ( kMultiRingNuisSKDetFSI008 ) :
    case ( kMultiRingNuisSKDetFSI009 ) :
    case ( kMultiRingNuisSKDetFSI010 ) :
    case ( kMultiRingNuisSKDetFSI011 ) :
    case ( kMultiRingNuisSKDetFSI012 ) :
    case ( kMultiRingNuisSKDetFSI013 ) :
    case ( kMultiRingNuisSKDetFSI014 ) :
    case ( kMultiRingNuisSKDetFSI015 ) :
    case ( kMultiRingNuisSKDetFSI016 ) :
    case ( kMultiRingNuisSKDetFSI017 ) :
    case ( kMultiRingNuisSKDetFSI018 ) :
    case ( kMultiRingNuisSKDetFSI019 ) :
      
        {
          // NEW SKDETFSI
	  int ip = ArrIdx::Instance()->MultiRingSKDetFSISystIdx(p);

          string nutypes = "";
          switch(p) {
          case (kMultiRingNuisSKDetFSI000) :
          case (kMultiRingNuisSKDetFSI001) :
          case (kMultiRingNuisSKDetFSI002) :
          case (kMultiRingNuisSKDetFSI003) :
          case (kMultiRingNuisSKDetFSI004) :
            nutypes = "signal nue CC (MultiRmu)"; break;
          case (kMultiRingNuisSKDetFSI005) :
          case (kMultiRingNuisSKDetFSI006) :
          case (kMultiRingNuisSKDetFSI007) :
          case (kMultiRingNuisSKDetFSI008) :
          case (kMultiRingNuisSKDetFSI009) :
            nutypes = "numu/numubar CC (MultiRmu)"; break;
          case (kMultiRingNuisSKDetFSI010) :
          case (kMultiRingNuisSKDetFSI011) :
          case (kMultiRingNuisSKDetFSI012) :
          case (kMultiRingNuisSKDetFSI013) :
          case (kMultiRingNuisSKDetFSI014) :
            nutypes = "nue/nuebar  CC (MultiRmu)"; break;
          case (kMultiRingNuisSKDetFSI015) :
          case (kMultiRingNuisSKDetFSI016) :
          case (kMultiRingNuisSKDetFSI017) :
          case (kMultiRingNuisSKDetFSI018) :
          case (kMultiRingNuisSKDetFSI019) :
            nutypes = "all NC (MultiRmu)"; break;
          default : nutypes = "--";
          }
          return std::string(Form("SKDet + FSI/SI %2d; Ereco range %.2f - %.2f GeV; %s; MultiRing",ip,
				  kMultiRingSKDetFSIParamsErecoCut[ip][0],kMultiRingSKDetFSIParamsErecoCut[ip][1],nutypes.c_str()));
          break;
        }//FSI for MultiRmu
      case ( kMultiRingNuisSKDetFSI020 ) :
      case ( kMultiRingNuisSKDetFSI021 ) :
      case ( kMultiRingNuisSKDetFSI022 ) :
      case ( kMultiRingNuisSKDetFSI023 ) :
      case ( kMultiRingNuisSKDetFSI024 ) :
      case ( kMultiRingNuisSKDetFSI025 ) :
      case ( kMultiRingNuisSKDetFSI026 ) :
      case ( kMultiRingNuisSKDetFSI027 ) :
      case ( kMultiRingNuisSKDetFSI028 ) :
      case ( kMultiRingNuisSKDetFSI029 ) :
      case ( kMultiRingNuisSKDetFSI030 ) :
      case ( kMultiRingNuisSKDetFSI031 ) :
        {
          // NEW SKDETFSI
	  int ip = ArrIdx::Instance()->MultiRingSKDetFSISystIdx(p);

          string nutypes = "";
          switch(p) {
          case (kMultiRingNuisSKDetFSI020) :
          case (kMultiRingNuisSKDetFSI021) :
          case (kMultiRingNuisSKDetFSI022) :
            nutypes = "oscillated nue CC (MultiRe)"; break;
          case (kMultiRingNuisSKDetFSI023) :
          case (kMultiRingNuisSKDetFSI024) :
          case (kMultiRingNuisSKDetFSI025) :
            nutypes = "numu/numubar CC (MultiRe)"; break;
          case (kMultiRingNuisSKDetFSI026) :
          case (kMultiRingNuisSKDetFSI027) :
          case (kMultiRingNuisSKDetFSI028) :
            nutypes = "nue/nuebar CC (MultiRe)"; break;
          case (kMultiRingNuisSKDetFSI029) :
          case (kMultiRingNuisSKDetFSI030) :
          case (kMultiRingNuisSKDetFSI031) :
            nutypes = "all NC (MultiRe)"; break;
          default : nutypes = "--";
          }
          return std::string(Form("SKDet + FSI/SI %2d; Ereco range %.2f - %.2f GeV; %s; MultiRing",ip,
				  kMultiRingSKDetFSIParamsErecoCut[ip][0],kMultiRingSKDetFSIParamsErecoCut[ip][1],nutypes.c_str()));
          break;
        }
#endif
        ////

      default :
        return "--";
      }
      return "--";
}

//________________________________________________________________________________________
string valor::util::AsStringTex(PhysParamId_t p)
{
  if (!util::IsSystematic(p)) return RawParamAsStringTex(p);
  // Specify tex strings all ALL subscripts and superscripts in braces.
  // Otherwise conversion to RooTex won't work.
  switch(p) {
            
    /*
     * 2011-analysis nuisance parameters
     */
            
    //
    // SK detector
    //
  case ( kNuisSuperkEnergyScale     ) : return "SK energy scale";                 break;
    //
    // BANFF flux and cross-section params
    //
  case ( kNuisBANFF00 ) :
  case ( kNuisBANFF01 ) :
  case ( kNuisBANFF02 ) :
  case ( kNuisBANFF03 ) :
  case ( kNuisBANFF04 ) :
  case ( kNuisBANFF05 ) :
  case ( kNuisBANFF06 ) :
  case ( kNuisBANFF07 ) :
  case ( kNuisBANFF08 ) :
  case ( kNuisBANFF09 ) :
  case ( kNuisBANFF10 ) :
  case ( kNuisBANFF11 ) :
  case ( kNuisBANFF12 ) :
  case ( kNuisBANFF13 ) :
  case ( kNuisBANFF14 ) :
  case ( kNuisBANFF15 ) :
#ifdef __USE_ALLFLUX__
  case ( kNuisBANFF16 ) :
  case ( kNuisBANFF17 ) :
  case ( kNuisBANFF18 ) :
  case ( kNuisBANFF19 ) :
  case ( kNuisBANFF20 ) :
  case ( kNuisBANFF21 ) :
  case ( kNuisBANFF22 ) :
  case ( kNuisBANFF23 ) :
  case ( kNuisBANFF24 ) :
#endif
    {
      int ip = (int)p - (int)kNuisBANFF00;
      return std::string(Form("BANFF %2d, SK %s flux, %.1f - %.1f GeV",
                              ip, util::Pdg2String(kNuPdgCodeOfFluxParamInBANFFfit[ip]).c_str(),
                              kEtrueRangeOfFluxParamInBANFFfit[ip][0],
                              kEtrueRangeOfFluxParamInBANFFfit[ip][1]));
      break;
    }
#ifndef ONEOVERMAQE
  case ( kNuisBANFFMAQE    ) : return "BANFF, Ma QE";                           break;
#else
  case ( kNuisBANFFMAQE    ) : return "BANFF, 1 / (Ma QE)$^{2}$";                 break;
#endif
  case ( kNuisBANFFMARES   ) : return "BANFF, Ma RES";                          break;
  case ( kNuisBANFFCCQEE0  ) : return "BANFF; Norm; CC QE; 0.0 - 1.5 GeV";      break;
#if BANFF_BINS != 3
  case ( kNuisBANFFCCQEE1  ) : return "BANFF; Norm; CC QE; 1.5 - 3.5 GeV";      break;
  case ( kNuisBANFFCCQEE2  ) : return "BANFF; Norm; CC QE; 3.5 - 30.0 GeV";     break;
#endif
  case ( kNuisBANFFCC1piE0 ) : return "BANFF; Norm; CC 1pi; 0.0 - 2.5 GeV";     break;
#if BANFF_BINS != 3
  case ( kNuisBANFFCC1piE1 ) : return "BANFF; Norm; CC 1pi; 2.5 - 30.0 GeV";    break;
#endif
#if (BANFF_BINS != 2)
  case ( kNuisBANFFNC1pi0  ) : return "BANFF; Norm; NC 1pi0";    break;
#endif
            
    // NEW XSEC
  case ( kNuisBANFFNorm2p2h      ) : return "BANFF; Norm; 2p2h";                          break;
  case ( kNuisBANFFCA5           ) : return "BANFF; CA5 RES";                             break;
  case ( kNuisBANFFBgRES         ) : return "BANFF; Norm; BgRES Isospin 1/2";             break;
    // NEW SPLIT ANA
#if (BANFF_BINS != 2)
  case ( kNuisBANFFSCCA          ) : return "BANFF; Second current class axial";          break;
  case ( kNuisBANFFSCCV          ) : return "BANFF; Second current class vector";         break;
#endif
  case ( kNuisBANFFFermiMom      ) : return "BANFF; Fermi Momentum";                      break; 
  case ( kNuisBANFFShapeCCoth    ) : return "BANFF; Shape; CC Oth";                       break; 
    // case ( kNuisBANFFBindE         ) : return "BANFF; Binding energy";                      break; 
  case ( kNuisBANFFNormCCcoh     ) : return "BANFF; Norm, CC Coh";                        break; 
  case ( kNuisBANFFNormNCcoh     ) : return "BANFF; Norm, NC Coh";                        break; 
  case ( kNuisBANFFNormNCoth     ) : return "BANFF; Norm, NC Oth";                        break; 
  case ( kNuisBANFFNormNuEToNuMu ) : return "BANFF; Norm, $\\nu_{e}$ To $\\nu_{\\mu}$";  break; 
  case ( kNuisBANFFNormNC1Gamma  ) : return "BANFF; Norm; NC 1 $\\gamma$";                break; 
  case ( kNuisBANFFNormNuEBarToNuMuBar ) : return "BANFF; Norm, $\\bar{\\nu}_{e}$ To $\\bar{\\nu}_{\\mu}$";  break; 
  case ( kNuisBANFFNorm2p2hBar   ) : return "BANFF; Norm; 2p2hBar";              break;
  case ( kNuisBANFFBeRPA_A       ) : return "BANFF; BeRPA A";                    break;
  case ( kNuisBANFFBeRPA_B       ) : return "BANFF; BeRPA B";                    break;
  case ( kNuisBANFFBeRPA_D       ) : return "BANFF; BeRPA D";                    break;
  case ( kNuisBANFFBeRPA_E       ) : return "BANFF; BeRPA E";                    break;
  case ( kNuisBANFFBeRPA_U       ) : return "BANFF; BeRPA U";                    break;
  case ( kNuisBANFFShape2p2h     ) : return "BANFF; Shape; 2p2h";         break;
  case ( kNuisBANFFNorm2p2hCtoO  ) : return "BANFF; Norm; 2p2h C to O";  break;
  case ( kNuisBANFFEb  ) : return "BANFF; Eb";
  case ( kNuisBANFFMultiReNorm  ) : return "BANFF; FHC MultiRe Normalisation";

    //
            
    //
    // XSec params not correlated with flux
    //
            
    // NEW XSEC
  case ( kNuisXSecRPA        ) : return "Tuning; RPA";                     break;
            
  case ( kNuisXSecFermiMom   ) : return "Fermi Momentum";                  break;
  case ( kNuisXSecSpecFunc   ) : return "Spectral Func";                   break;
  case ( kNuisXSecWshape     ) : return "W-shape";                         break;
  case ( kNuisXSecShapeCCoth ) : return "Shape; CC Oth";                   break;
  case ( kNuisXSecPiLessDelta) : return "Pi-less delta decay";             break;
#ifndef __NUEONLY__
  case ( kNuisXSecBindE      ) : return "Binding energy";                  break;
#endif
  case ( kNuisXSecNormCCcoh  ) : return "Norm, CC Coh";                    break;
#ifndef __NUEONLY__
  case ( kNuisXSecNormNC1piC ) : return "Norm, NC 1piC";                   break;
#endif
#ifndef __NUMUONLY__
  case ( kNuisXSecNormNCcoh  ) : return "Norm, NC Coh";                    break;
#endif
  case ( kNuisXSecNormNCoth  ) : return "Norm, NC Oth";                    break;
  case ( kNuisXSecNormNuEToNuMu  ) : return "Norm, $\\nu_{e}$ To $\\nu_{\\mu}$";  break;
#ifndef __NUEONLY__
  case ( kNuisXSecNormAntiNuToNu ) : return "Norm, $\\bar{\\nu}$ To $\\nu$";  break;
#endif
            
    //
    // SKDET+FSI
    //
    //#ifdef SKFSITOGETHER
#ifndef SKCORRELATIONS
#if FSI_BINS != 3
  case ( kNuisSKDetFSI000 ) :
  case ( kNuisSKDetFSI001 ) :
  case ( kNuisSKDetFSI002 ) :
  case ( kNuisSKDetFSI003 ) :
  case ( kNuisSKDetFSI004 ) :
  case ( kNuisSKDetFSI005 ) :
            
    // NEW SKDETFSI
#if defined(__USE_SKDETFSI2015__)
  case ( kNuisSKDetFSI006 ) :
  case ( kNuisSKDetFSI007 ) :
  case ( kNuisSKDetFSI008 ) :
  case ( kNuisSKDetFSI009 ) :
  case ( kNuisSKDetFSI010 ) :
  case ( kNuisSKDetFSI011 ) :
#endif
            
    {
      string nutypes = "";
            
      // NEW SKDETFSI
      //int ip = (int)p - (int)kNuisSKDetFSI000;
      int ip = (int)p - (int)kNuisSKDetFSI_first;
            
      switch(p) {
                    
        // NEW SKDETFSI
#if defined(__USE_SKDETFSI2013__)
      case (kNuisSKDetFSI000) :
      case (kNuisSKDetFSI001) :
      case (kNuisSKDetFSI002) :
        nutypes = "$\\nu_{\\mu}$/$\\bar{\\nu}_{\\mu}$ CCQE ($1R_{\\mu}$)"; break;
      case (kNuisSKDetFSI003) :
        nutypes = "$\\nu_{\\mu}$/$\\bar{\\nu}_{\\mu}$ CCnQE ($1R_{\\mu}$)"; break;
      case (kNuisSKDetFSI004) :
        nutypes = "$\\nu_{e}$/$\\bar{\\nu}_{e}$/sig$\\nu_{e}$ CC ($1R_{\\mu}$)"; break;
      case (kNuisSKDetFSI005) :
        nutypes = "all NC ($1R_{\\mu}$)"; break;
#elif defined(__USE_SKDETFSI2015__)
      case (kNuisSKDetFSI000) :
      case (kNuisSKDetFSI001) :
      case (kNuisSKDetFSI002) :
        nutypes = "$\\nu_{\\mu}$/$\\bar{\\nu}_{\\mu}$ CCQE ($1R_{\\mu}$)"; break;
      case (kNuisSKDetFSI003) :
      case (kNuisSKDetFSI004) :
      case (kNuisSKDetFSI005) :
        nutypes = "$\\nu_{\\mu}$/$\\bar{\\nu}_{\\mu}$ CCnQE ($1R_{\\mu}$)"; break;
      case (kNuisSKDetFSI006) :
      case (kNuisSKDetFSI007) :
      case (kNuisSKDetFSI008) :
        nutypes = "$\\nu_{e}$/$\\bar{\\nu}_{e}$/sig$\\nu_{e}$ CC ($1R_{\\mu}$)"; break;
      case (kNuisSKDetFSI009) :
      case (kNuisSKDetFSI010) :
      case (kNuisSKDetFSI011) :
        nutypes = "all NC ($1R_{\\mu}$)"; break;
#endif
                    
      default : nutypes = "--";
      }
      return std::string(Form("SKDet + FSI/SI %2d; $E_{reco}$( %.2f - %.2f )GeV; %s",ip,
                              kSKDetFSIParamsErecoCut[ip][0],kSKDetFSIParamsErecoCut[ip][1],nutypes.c_str()));
      break;
    }//FSI for 1Rmu
#endif
#if FSI_BINS != 2
            
  // NEW SKDETFSI
#if defined(__USE_SKDETFSI2013__)
  case ( kNuisSKDetFSI006 ) :
  case ( kNuisSKDetFSI007 ) :
  case ( kNuisSKDetFSI008 ) :
  case ( kNuisSKDetFSI009 ) :
  case ( kNuisSKDetFSI010 ) :
  case ( kNuisSKDetFSI011 ) :
  case ( kNuisSKDetFSI012 ) :
  case ( kNuisSKDetFSI013 ) :
  case ( kNuisSKDetFSI014 ) :
  case ( kNuisSKDetFSI015 ) :
  case ( kNuisSKDetFSI016 ) :
  case ( kNuisSKDetFSI017 ) :
#elif defined(__USE_SKDETFSI2015__)
  case ( kNuisSKDetFSI012 ) :
  case ( kNuisSKDetFSI013 ) :
  case ( kNuisSKDetFSI014 ) :
  case ( kNuisSKDetFSI015 ) :
  case ( kNuisSKDetFSI016 ) :
  case ( kNuisSKDetFSI017 ) :
  case ( kNuisSKDetFSI018 ) :
  case ( kNuisSKDetFSI019 ) :
  case ( kNuisSKDetFSI020 ) :
  case ( kNuisSKDetFSI021 ) :
  case ( kNuisSKDetFSI022 ) :
  case ( kNuisSKDetFSI023 ) :
#endif
            
    {
      string nutypes = "";
            
      // #if FSI_BINS != 3
      //       int ip = (int)p - (int)kNuisSKDetFSI000;
      // #else
      //       #if defined(__USE_SKDETFSI2013__)
      //       int ip = (int)p - (int)kNuisSKDetFSI006;
      //       #elif defined(__USE_SKDETFSI2015__)
      //       int ip = (int)p - (int)kNuisSKDetFSI012;
      //       #endif
      // #endif
            
      // NEW SKDETFSI
      int ip = (int)p - (int)kNuisSKDetFSI_first;
            
      switch(p) {
                    
        // NEW SKDETFSI
#if defined(__USE_SKDETFSI2013__)
      case (kNuisSKDetFSI006) :
      case (kNuisSKDetFSI007) :
      case (kNuisSKDetFSI008) :
        nutypes = "Osc. $\\nu_{e}$ CC ($1R_{e}$)"; break;
      case (kNuisSKDetFSI009) :
      case (kNuisSKDetFSI010) :
      case (kNuisSKDetFSI011) :
        nutypes = "$\\nu_{\\mu}$/$\\bar{\\nu}_{\\mu}$ CC ($1R_{e}$)"; break;
      case (kNuisSKDetFSI012) :
      case (kNuisSKDetFSI013) :
      case (kNuisSKDetFSI014) :
        nutypes = "$\\nu_{e}$/$\\bar{\\nu}_{e}$ CC ($1R_{e}$)"; break;
      case (kNuisSKDetFSI015) :
      case (kNuisSKDetFSI016) :
      case (kNuisSKDetFSI017) :
        nutypes = "all NC ($1R_{e}$)"; break;
#elif defined(__USE_SKDETFSI2015__)
      case (kNuisSKDetFSI012) :
      case (kNuisSKDetFSI013) :
      case (kNuisSKDetFSI014) :
        nutypes = "Osc. $\\nu_{e}$ CC ($1R_{e}$)"; break;
      case (kNuisSKDetFSI015) :
      case (kNuisSKDetFSI016) :
      case (kNuisSKDetFSI017) :
        nutypes = "$\\nu_{\\mu}$/$\\bar{\\nu}_{\\mu}$ CC ($1R_{e}$)"; break;
      case (kNuisSKDetFSI018) :
      case (kNuisSKDetFSI019) :
      case (kNuisSKDetFSI020) :
        nutypes = "$\\nu_{e}$/$\\bar{\\nu}_{e}$ CC ($1R_{e}$)"; break;
      case (kNuisSKDetFSI021) :
      case (kNuisSKDetFSI022) :
      case (kNuisSKDetFSI023) :
        nutypes = "all NC ($1R_{e}$)"; break;
#endif
                    
      default : nutypes = "--";
      }
      return std::string(Form("SKDet + FSI/SI %2d; $E_{reco}$( %.2f - %.2f )GeV; %s",ip,
                              kSKDetFSIParamsErecoCut[ip][0],kSKDetFSIParamsErecoCut[ip][1],nutypes.c_str()));
      break;
    }//FSI for 1Re
#endif
#else
            
  // OTHER FSI_BINS THAN 1,2,3
            
  case ( kNuisSKDetFSI000 ) :
  case ( kNuisSKDetFSI001 ) :
  case ( kNuisSKDetFSI002 ) :
  case ( kNuisSKDetFSI003 ) :
  case ( kNuisSKDetFSI004 ) :
  case ( kNuisSKDetFSI005 ) :
  case ( kNuisSKDetFSI006 ) :
  case ( kNuisSKDetFSI007 ) :
  case ( kNuisSKDetFSI008 ) :
  case ( kNuisSKDetFSI009 ) :
  case ( kNuisSKDetFSI010 ) :
  case ( kNuisSKDetFSI011 ) :
  case ( kNuisSKDetFSI012 ) :
  case ( kNuisSKDetFSI013 ) :
  case ( kNuisSKDetFSI014 ) :
  case ( kNuisSKDetFSI015 ) :
  case ( kNuisSKDetFSI016 ) :
  case ( kNuisSKDetFSI017 ) :
  case ( kNuisSKDetFSI018 ) :
  case ( kNuisSKDetFSI019 ) :
    {
      // NEW SKDETFSI
      //int ip = (int)p - (int)kNuisSKDetFSI000;
      int ip = (int)p - (int)kNuisSKDetFSI_first;
            
      string nutypes = "";
      switch(p) {
      case (kNuisSKDetFSI000) :
      case (kNuisSKDetFSI001) :
      case (kNuisSKDetFSI002) :
      case (kNuisSKDetFSI003) :
      case (kNuisSKDetFSI004) :
        nutypes = "Osc. $\\nu_{e}$ CC ($1R_{\\mu}$)"; break;
      case (kNuisSKDetFSI005) :
      case (kNuisSKDetFSI006) :
      case (kNuisSKDetFSI007) :
      case (kNuisSKDetFSI008) :
      case (kNuisSKDetFSI009) :
        nutypes = "$\\nu_{\\mu}$/$\\bar{\\nu}_{\\mu}$ CC ($1R_{\\mu}$)"; break;
      case (kNuisSKDetFSI010) :
      case (kNuisSKDetFSI011) :
      case (kNuisSKDetFSI012) :
      case (kNuisSKDetFSI013) :
      case (kNuisSKDetFSI014) :
        nutypes = "$\\nu_{e}$/$\\bar{\\nu}_{e}$ CC ($1R_{\\mu}$)"; break;
      case (kNuisSKDetFSI015) :
      case (kNuisSKDetFSI016) :
      case (kNuisSKDetFSI017) :
      case (kNuisSKDetFSI018) :
      case (kNuisSKDetFSI019) :
        nutypes = "all NC ($1R_{\\mu}$)"; break;
      default : nutypes = "--";
      }
      return std::string(Form("SKDet + FSI/SI %2d; $E_{reco}$( %.2f - %.2f )GeV; %s",ip,
                              kSKDetFSIParamsErecoCut[ip][0],kSKDetFSIParamsErecoCut[ip][1],nutypes.c_str()));
      break;
    }//FSI for 1Rmu
  case ( kNuisSKDetFSI020 ) :
  case ( kNuisSKDetFSI021 ) :
  case ( kNuisSKDetFSI022 ) :
  case ( kNuisSKDetFSI023 ) :
  case ( kNuisSKDetFSI024 ) :
  case ( kNuisSKDetFSI025 ) :
  case ( kNuisSKDetFSI026 ) :
  case ( kNuisSKDetFSI027 ) :
  case ( kNuisSKDetFSI028 ) :
  case ( kNuisSKDetFSI029 ) :
  case ( kNuisSKDetFSI030 ) :
  case ( kNuisSKDetFSI031 ) :
    {
      // NEW SKDETFSI
      //int ip = (int)p - (int)kNuisSKDetFSI000;
      int ip = (int)p - (int)kNuisSKDetFSI_first;
            
      string nutypes = "";
      switch(p) {
      case (kNuisSKDetFSI020) :
      case (kNuisSKDetFSI021) :
      case (kNuisSKDetFSI022) :
        nutypes = "Osc. $\\nu_{e}$ CC ($1R_{e}$)"; break;
      case (kNuisSKDetFSI023) :
      case (kNuisSKDetFSI024) :
      case (kNuisSKDetFSI025) :
        nutypes = "$\\nu_{\\mu}$/$\\bar{\\nu}_{\\mu}$ CC ($1R_{e}$)"; break;
      case (kNuisSKDetFSI026) :
      case (kNuisSKDetFSI027) :
      case (kNuisSKDetFSI028) :
        nutypes = "$\\nu_{e}$/$\\bar{\\nu}_{e}$ CC ($1R_{e}$)"; break;
      case (kNuisSKDetFSI029) :
      case (kNuisSKDetFSI030) :
      case (kNuisSKDetFSI031) :
        nutypes = "all NC ($1R_{e}$)"; break;
      default : nutypes = "--";
      }
      return std::string(Form("SKDet + FSI/SI %2d; $E_{reco}$( %.2f - %.2f )GeV; %s",ip,
                              kSKDetFSIParamsErecoCut[ip][0],kSKDetFSIParamsErecoCut[ip][1],nutypes.c_str()));
      break;
    }
#endif
  //#endif //SKDETFSITOGETHER
            
            
  /////////////// REVERSE HORN CURRENT /////////////
            
  case ( kNuisRHCToFHC ) : return "RHC to FHC normalisation"; break;
            
            
    //
    // SK detector
    //
  case ( kRHCNuisSuperkEnergyScale     ) : return "SK energy scale";                 break;
    //
    // BANFF flux and cross-section params
    //
  case ( kRHCNuisBANFF00 ) :
  case ( kRHCNuisBANFF01 ) :
  case ( kRHCNuisBANFF02 ) :
  case ( kRHCNuisBANFF03 ) :
  case ( kRHCNuisBANFF04 ) :
  case ( kRHCNuisBANFF05 ) :
  case ( kRHCNuisBANFF06 ) :
  case ( kRHCNuisBANFF07 ) :
  case ( kRHCNuisBANFF08 ) :
  case ( kRHCNuisBANFF09 ) :
  case ( kRHCNuisBANFF10 ) :
  case ( kRHCNuisBANFF11 ) :
  case ( kRHCNuisBANFF12 ) :
  case ( kRHCNuisBANFF13 ) :
  case ( kRHCNuisBANFF14 ) :
  case ( kRHCNuisBANFF15 ) :
#ifdef __USE_ALLFLUX__
  case ( kRHCNuisBANFF16 ) :
  case ( kRHCNuisBANFF17 ) :
  case ( kRHCNuisBANFF18 ) :
  case ( kRHCNuisBANFF19 ) :
  case ( kRHCNuisBANFF20 ) :
  case ( kRHCNuisBANFF21 ) :
  case ( kRHCNuisBANFF22 ) :
  case ( kRHCNuisBANFF23 ) :
  case ( kRHCNuisBANFF24 ) :
#endif
    {
      int ip = (int)p - (int)kRHCNuisBANFF00;
      return std::string(Form("BANFF %2d, SK %s flux, %.1f - %.1f GeV; RHC",
                              ip, util::Pdg2String(kNuPdgCodeOfRHCFluxParamInBANFFfit[ip]).c_str(),
                              kEtrueRangeOfRHCFluxParamInBANFFfit[ip][0],
                              kEtrueRangeOfRHCFluxParamInBANFFfit[ip][1]));
      break;
    }
#ifndef ONEOVERMAQE
  case ( kRHCNuisBANFFMAQE    ) : return "BANFF, Ma QE; RHC";                           break;
#else
  case ( kRHCNuisBANFFMAQE    ) : return "BANFF, 1 / (Ma QE)$^{2}; RHC$";                 break;
#endif
  case ( kRHCNuisBANFFMARES   ) : return "BANFF, Ma RES; RHC";                          break;
  case ( kRHCNuisBANFFCCQEE0  ) : return "BANFF; Norm; CC QE; 0.0 - 1.5 GeV; RHC";      break;
#if BANFF_BINS != 3
  case ( kRHCNuisBANFFCCQEE1  ) : return "BANFF; Norm; CC QE; 1.5 - 3.5 GeV; RHC";      break;
  case ( kRHCNuisBANFFCCQEE2  ) : return "BANFF; Norm; CC QE; 3.5 - 30.0 GeV; RHC";     break;
#endif
  case ( kRHCNuisBANFFCC1piE0 ) : return "BANFF; Norm; CC 1pi; 0.0 - 2.5 GeV; RHC";     break;
#if BANFF_BINS != 3
  case ( kRHCNuisBANFFCC1piE1 ) : return "BANFF; Norm; CC 1pi; 2.5 - 30.0 GeV; RHC";    break;
#endif
#if (BANFF_BINS != 2)
  case ( kRHCNuisBANFFNC1pi0  ) : return "BANFF; Norm; NC 1pi0; RHC";    break;
#endif
            
    // NEW XSEC
  case ( kRHCNuisBANFF2p2h           ) : return "BANFF; Norm; 2p2h; RHC";                           break;
  case ( kRHCNuisBANFFCA5           ) : return "BANFF; CA5 RES; RHC";                             break;
  case ( kRHCNuisBANFFBgRES         ) : return "BANFF; Norm; BgRES Isospin 1/2; RHC";             break;
    // NEW SPLIT ANA
#if (BANFF_BINS != 2)
  case ( kRHCNuisBANFFSCCA          ) : return "BANFF; Second current class axial; RHC";          break;
  case ( kRHCNuisBANFFSCCV          ) : return "BANFF; Second current class vector; RHC";         break;
#endif
  case ( kRHCNuisBANFFFermiMom      ) : return "BANFF; Fermi Momentum; RHC";                      break; // uncorrelated in the past
  case ( kRHCNuisBANFFShapeCCoth    ) : return "BANFF; Shape; CC Oth; RHC";                       break; // uncorrelated in the past
    // case ( kRHCNuisBANFFBindE         ) : return "BANFF; Binding energy; RHC";                      break; // uncorrelated in the past
  case ( kRHCNuisBANFFNormCCcoh     ) : return "BANFF; Norm, CC Coh; RHC";                        break; // uncorrelated in the past
  case ( kRHCNuisBANFFNormNCcoh     ) : return "BANFF; Norm, NC Coh; RHC";                        break; // uncorrelated in the past
  case ( kRHCNuisBANFFNormNCoth     ) : return "BANFF; Norm, NC Oth; RHC";                        break; // uncorrelated in the past
  case ( kRHCNuisBANFFNormNuEToNuMu ) : return "BANFF; Norm, $\\nu_{e}$ To $\\nu_{\\mu}$; RHC";  break; // uncorrelated in the past
  case ( kRHCNuisBANFFNormNC1Gamma      ) : return "BANFF; Norm; NC 1 $\\gamma$; RHC";                break; // uncorrelated in the past
  case ( kRHCNuisBANFFNormNuEBarToNuMuBar ) : return "BANFF; Norm, $\\bar{\\nu}_{e}$ To $\\bar{\\nu}_{\\mu}$; RHC";  break; // uncorrelated in the past
  case ( kRHCNuisBANFF2p2hBar        ) : return "BANFF; Norm; 2p2hBar; RHC";                           break;
            
    //
            
    //
    // XSec params not correlated with flux
    //
            
    // NEW XSEC
  case ( kRHCNuisXSecRPA        ) : return "Tuning; RPA; RHC";                     break;
            
  case ( kRHCNuisXSecFermiMom   ) : return "Fermi Momentum; RHC";                  break;
  case ( kRHCNuisXSecSpecFunc   ) : return "Spectral Func; RHC";                   break;
  case ( kRHCNuisXSecWshape     ) : return "W-shape; RHC";                         break;
  case ( kRHCNuisXSecShapeCCoth ) : return "Shape; CC Oth; RHC";                   break;
  case ( kRHCNuisXSecPiLessDelta) : return "Pi-less delta decay; RHC";             break;
#ifndef __NUEONLY__
  case ( kRHCNuisXSecBindE      ) : return "Binding energy; RHC";                  break;
#endif
  case ( kRHCNuisXSecNormCCcoh  ) : return "Norm, CC Coh; RHC";                    break;
#ifndef __NUEONLY__
  case ( kRHCNuisXSecNormNC1piC ) : return "Norm, NC 1piC; RHC";                   break;
#endif
#ifndef __NUMUONLY__
  case ( kRHCNuisXSecNormNCcoh  ) : return "Norm, NC Coh; RHC";                    break;
#endif
  case ( kRHCNuisXSecNormNCoth  ) : return "Norm, NC Oth; RHC";                    break;
  case ( kRHCNuisXSecNormNuEToNuMu  ) : return "Norm, $\\nu_{e}$ To $\\nu_{\\mu}$; RHC";  break;
#ifndef __NUEONLY__
  case ( kRHCNuisXSecNormAntiNuToNu ) : return "Norm, $\\bar{\\nu}$ To $\\nu$; RHC";  break;
#endif
            
    //
    // SKDET+FSI
    //
    //#ifdef SKFSITOGETHER
#ifndef SKCORRELATIONS
#if FSI_BINS != 3
  case ( kRHCNuisSKDetFSI000 ) :
  case ( kRHCNuisSKDetFSI001 ) :
  case ( kRHCNuisSKDetFSI002 ) :
  case ( kRHCNuisSKDetFSI003 ) :
  case ( kRHCNuisSKDetFSI004 ) :
  case ( kRHCNuisSKDetFSI005 ) :
            
    // NEW SKDETFSI
#if defined(__USE_SKDETFSI2015__)
  case ( kRHCNuisSKDetFSI006 ) :
  case ( kRHCNuisSKDetFSI007 ) :
  case ( kRHCNuisSKDetFSI008 ) :
  case ( kRHCNuisSKDetFSI009 ) :
  case ( kRHCNuisSKDetFSI010 ) :
  case ( kRHCNuisSKDetFSI011 ) :
#endif
            
    {
      string nutypes = "";
            
      // NEW SKDETFSI
      //int ip = (int)p - (int)kRHCNuisSKDetFSI000;
      int ip = (int)p - (int)kRHCNuisSKDetFSI_first;
            
      switch(p) {
                    
        // NEW SKDETFSI
#if defined(__USE_SKDETFSI2013__)
      case (kRHCNuisSKDetFSI000) :
      case (kRHCNuisSKDetFSI001) :
      case (kRHCNuisSKDetFSI002) :
        nutypes = "$\\nu_{\\mu}$/$\\bar{\\nu}_{\\mu}$ CCQE ($1R_{\\mu}$)"; break;
      case (kRHCNuisSKDetFSI003) :
        nutypes = "$\\nu_{\\mu}$/$\\bar{\\nu}_{\\mu}$ CCnQE ($1R_{\\mu}$)"; break;
      case (kRHCNuisSKDetFSI004) :
        nutypes = "$\\nu_{e}$/$\\bar{\\nu}_{e}$/sig$\\nu_{e}$ CC ($1R_{\\mu}$)"; break;
      case (kRHCNuisSKDetFSI005) :
        nutypes = "all NC ($1R_{\\mu}$)"; break;
#elif defined(__USE_SKDETFSI2015__)
      case (kRHCNuisSKDetFSI000) :
      case (kRHCNuisSKDetFSI001) :
      case (kRHCNuisSKDetFSI002) :
        nutypes = "$\\nu_{\\mu}$/$\\bar{\\nu}_{\\mu}$ CCQE ($1R_{\\mu}$)"; break;
      case (kRHCNuisSKDetFSI003) :
      case (kRHCNuisSKDetFSI004) :
      case (kRHCNuisSKDetFSI005) :
        nutypes = "$\\nu_{\\mu}$/$\\bar{\\nu}_{\\mu}$ CCnQE ($1R_{\\mu}$)"; break;
      case (kRHCNuisSKDetFSI006) :
      case (kRHCNuisSKDetFSI007) :
      case (kRHCNuisSKDetFSI008) :
        nutypes = "$\\nu_{e}$/$\\bar{\\nu}_{e}$/sig$\\nu_{e}$ CC ($1R_{\\mu}$)"; break;
      case (kRHCNuisSKDetFSI009) :
      case (kRHCNuisSKDetFSI010) :
      case (kRHCNuisSKDetFSI011) :
        nutypes = "all NC ($1R_{\\mu}$)"; break;
#endif
                    
      default : nutypes = "--";
      }
      return std::string(Form("SKDet + FSI/SI %2d; $E_{reco}$( %.2f - %.2f )GeV; %s; RHC",ip,
                              kRHCSKDetFSIParamsErecoCut[ip][0],kRHCSKDetFSIParamsErecoCut[ip][1],nutypes.c_str()));
      break;
    }//FSI for 1Rmu
#endif
#if FSI_BINS != 2
            
  // NEW SKDETFSI
#if defined(__USE_SKDETFSI2013__)
  case ( kRHCNuisSKDetFSI006 ) :
  case ( kRHCNuisSKDetFSI007 ) :
  case ( kRHCNuisSKDetFSI008 ) :
  case ( kRHCNuisSKDetFSI009 ) :
  case ( kRHCNuisSKDetFSI010 ) :
  case ( kRHCNuisSKDetFSI011 ) :
  case ( kRHCNuisSKDetFSI012 ) :
  case ( kRHCNuisSKDetFSI013 ) :
  case ( kRHCNuisSKDetFSI014 ) :
  case ( kRHCNuisSKDetFSI015 ) :
  case ( kRHCNuisSKDetFSI016 ) :
  case ( kRHCNuisSKDetFSI017 ) :
#elif defined(__USE_SKDETFSI2015__)
  case ( kRHCNuisSKDetFSI012 ) :
  case ( kRHCNuisSKDetFSI013 ) :
  case ( kRHCNuisSKDetFSI014 ) :
  case ( kRHCNuisSKDetFSI015 ) :
  case ( kRHCNuisSKDetFSI016 ) :
  case ( kRHCNuisSKDetFSI017 ) :
  case ( kRHCNuisSKDetFSI018 ) :
  case ( kRHCNuisSKDetFSI019 ) :
  case ( kRHCNuisSKDetFSI020 ) :
  case ( kRHCNuisSKDetFSI021 ) :
  case ( kRHCNuisSKDetFSI022 ) :
  case ( kRHCNuisSKDetFSI023 ) :
#endif
            
    {
      string nutypes = "";
            
      // #if FSI_BINS != 3
      //       int ip = (int)p - (int)kRHCNuisSKDetFSI000;
      // #else
      //       #if defined(__USE_SKDETFSI2013__)
      //       int ip = (int)p - (int)kRHCNuisSKDetFSI006;
      //       #elif defined(__USE_SKDETFSI2015__)
      //       int ip = (int)p - (int)kRHCNuisSKDetFSI012;
      //       #endif
      // #endif
            
      // NEW SKDETFSI
      int ip = (int)p - (int)kRHCNuisSKDetFSI_first;
            
      switch(p) {
                    
        // NEW SKDETFSI
#if defined(__USE_SKDETFSI2013__)
      case (kRHCNuisSKDetFSI006) :
      case (kRHCNuisSKDetFSI007) :
      case (kRHCNuisSKDetFSI008) :
        nutypes = "Osc. $\\nu_{e}$ CC ($1R_{e}$)"; break;
      case (kRHCNuisSKDetFSI009) :
      case (kRHCNuisSKDetFSI010) :
      case (kRHCNuisSKDetFSI011) :
        nutypes = "$\\nu_{\\mu}$/$\\bar{\\nu}_{\\mu}$ CC ($1R_{e}$)"; break;
      case (kRHCNuisSKDetFSI012) :
      case (kRHCNuisSKDetFSI013) :
      case (kRHCNuisSKDetFSI014) :
        nutypes = "$\\nu_{e}$/$\\bar{\\nu}_{e}$ CC ($1R_{e}$)"; break;
      case (kRHCNuisSKDetFSI015) :
      case (kRHCNuisSKDetFSI016) :
      case (kRHCNuisSKDetFSI017) :
        nutypes = "all NC ($1R_{e}$)"; break;
#elif defined(__USE_SKDETFSI2015__)
      case (kRHCNuisSKDetFSI012) :
      case (kRHCNuisSKDetFSI013) :
      case (kRHCNuisSKDetFSI014) :
        nutypes = "Osc. $\\nu_{e}$ CC ($1R_{e}$)"; break;
      case (kRHCNuisSKDetFSI015) :
      case (kRHCNuisSKDetFSI016) :
      case (kRHCNuisSKDetFSI017) :
        nutypes = "$\\nu_{\\mu}$/$\\bar{\\nu}_{\\mu}$ CC ($1R_{e}$)"; break;
      case (kRHCNuisSKDetFSI018) :
      case (kRHCNuisSKDetFSI019) :
      case (kRHCNuisSKDetFSI020) :
        nutypes = "$\\nu_{e}$/$\\bar{\\nu}_{e}$ CC ($1R_{e}$)"; break;
      case (kRHCNuisSKDetFSI021) :
      case (kRHCNuisSKDetFSI022) :
      case (kRHCNuisSKDetFSI023) :
        nutypes = "all NC ($1R_{e}$)"; break;
#endif
                    
      default : nutypes = "--";
      }
      return std::string(Form("SKDet + FSI/SI %2d; $E_{reco}$( %.2f - %.2f )GeV; %s; RHC",ip,
                              kRHCSKDetFSIParamsErecoCut[ip][0],kRHCSKDetFSIParamsErecoCut[ip][1],nutypes.c_str()));
      break;
    }//FSI for 1Re
#endif
#else
            
  // OTHER FSI_BINS THAN 1,2,3
            
  case ( kRHCNuisSKDetFSI000 ) :
  case ( kRHCNuisSKDetFSI001 ) :
  case ( kRHCNuisSKDetFSI002 ) :
  case ( kRHCNuisSKDetFSI003 ) :
  case ( kRHCNuisSKDetFSI004 ) :
  case ( kRHCNuisSKDetFSI005 ) :
  case ( kRHCNuisSKDetFSI006 ) :
  case ( kRHCNuisSKDetFSI007 ) :
  case ( kRHCNuisSKDetFSI008 ) :
  case ( kRHCNuisSKDetFSI009 ) :
  case ( kRHCNuisSKDetFSI010 ) :
  case ( kRHCNuisSKDetFSI011 ) :
  case ( kRHCNuisSKDetFSI012 ) :
  case ( kRHCNuisSKDetFSI013 ) :
  case ( kRHCNuisSKDetFSI014 ) :
  case ( kRHCNuisSKDetFSI015 ) :
  case ( kRHCNuisSKDetFSI016 ) :
  case ( kRHCNuisSKDetFSI017 ) :
  case ( kRHCNuisSKDetFSI018 ) :
  case ( kRHCNuisSKDetFSI019 ) :
    {
      // NEW SKDETFSI
      //int ip = (int)p - (int)kRHCNuisSKDetFSI000;
      int ip = (int)p - (int)kRHCNuisSKDetFSI_first;
            
      string nutypes = "";
      switch(p) {
      case (kRHCNuisSKDetFSI000) :
      case (kRHCNuisSKDetFSI001) :
      case (kRHCNuisSKDetFSI002) :
      case (kRHCNuisSKDetFSI003) :
      case (kRHCNuisSKDetFSI004) :
        nutypes = "Osc. $\\nu_{e}$ CC ($1R_{\\mu}$)"; break;
      case (kRHCNuisSKDetFSI005) :
      case (kRHCNuisSKDetFSI006) :
      case (kRHCNuisSKDetFSI007) :
      case (kRHCNuisSKDetFSI008) :
      case (kRHCNuisSKDetFSI009) :
        nutypes = "$\\nu_{\\mu}$/$\\bar{\\nu}_{\\mu}$ CC ($1R_{\\mu}$)"; break;
      case (kRHCNuisSKDetFSI010) :
      case (kRHCNuisSKDetFSI011) :
      case (kRHCNuisSKDetFSI012) :
      case (kRHCNuisSKDetFSI013) :
      case (kRHCNuisSKDetFSI014) :
        nutypes = "$\\nu_{e}$/$\\bar{\\nu}_{e}$ CC ($1R_{\\mu}$)"; break;
      case (kRHCNuisSKDetFSI015) :
      case (kRHCNuisSKDetFSI016) :
      case (kRHCNuisSKDetFSI017) :
      case (kRHCNuisSKDetFSI018) :
      case (kRHCNuisSKDetFSI019) :
        nutypes = "all NC ($1R_{\\mu}$)"; break;
      default : nutypes = "--";
      }
      return std::string(Form("SKDet + FSI/SI %2d; $E_{reco}$( %.2f - %.2f )GeV; %s; RHC",ip,
                              kRHCSKDetFSIParamsErecoCut[ip][0],kRHCSKDetFSIParamsErecoCut[ip][1],nutypes.c_str()));
      break;
    }//FSI for 1Rmu
  case ( kRHCNuisSKDetFSI020 ) :
  case ( kRHCNuisSKDetFSI021 ) :
  case ( kRHCNuisSKDetFSI022 ) :
  case ( kRHCNuisSKDetFSI023 ) :
  case ( kRHCNuisSKDetFSI024 ) :
  case ( kRHCNuisSKDetFSI025 ) :
  case ( kRHCNuisSKDetFSI026 ) :
  case ( kRHCNuisSKDetFSI027 ) :
  case ( kRHCNuisSKDetFSI028 ) :
  case ( kRHCNuisSKDetFSI029 ) :
  case ( kRHCNuisSKDetFSI030 ) :
  case ( kRHCNuisSKDetFSI031 ) :
    {
      // NEW SKDETFSI
      //int ip = (int)p - (int)kRHCNuisSKDetFSI000;
      int ip = (int)p - (int)kRHCNuisSKDetFSI_first;
            
      string nutypes = "";
      switch(p) {
      case (kRHCNuisSKDetFSI020) :
      case (kRHCNuisSKDetFSI021) :
      case (kRHCNuisSKDetFSI022) :
        nutypes = "Osc. $\\nu_{e}$ CC ($1R_{e}$)"; break;
      case (kRHCNuisSKDetFSI023) :
      case (kRHCNuisSKDetFSI024) :
      case (kRHCNuisSKDetFSI025) :
        nutypes = "$\\nu_{\\mu}$/$\\bar{\\nu}_{\\mu}$ CC ($1R_{e}$)"; break;
      case (kRHCNuisSKDetFSI026) :
      case (kRHCNuisSKDetFSI027) :
      case (kRHCNuisSKDetFSI028) :
        nutypes = "$\\nu_{e}$/$\\bar{\\nu}_{e}$ CC ($1R_{e}$)"; break;
      case (kRHCNuisSKDetFSI029) :
      case (kRHCNuisSKDetFSI030) :
      case (kRHCNuisSKDetFSI031) :
        nutypes = "all NC ($1R_{e}$)"; break;
      default : nutypes = "--";
      }
      return std::string(Form("SKDet + FSI/SI %2d; $E_{reco}$( %.2f - %.2f )GeV; %s; RHC",ip,
                              kRHCSKDetFSIParamsErecoCut[ip][0],kRHCSKDetFSIParamsErecoCut[ip][1],nutypes.c_str()));
      break;
    }
#endif
  //#endif //SKDETFSITOGETHER
            
            
            
  /////////////// MULTIRING FORWARD HORN CURRENT /////////////
            
  /// NEW MULTIRING 3
  case ( kMultiRingNuisSuperkEnergyScale ) : return "SK energy scale; MultiRing"; break;
    ///
            
    /// NEW MULTIRING 2
            
    //
    // SKDET+FSI
    //
    //#ifdef SKFSITOGETHER
#ifndef SKCORRELATIONS
#if MR_FSI_BINS != 3
  case ( kMultiRingNuisSKDetFSI000 ) :
  case ( kMultiRingNuisSKDetFSI001 ) :
  case ( kMultiRingNuisSKDetFSI002 ) :
  case ( kMultiRingNuisSKDetFSI003 ) :
  case ( kMultiRingNuisSKDetFSI004 ) :
  case ( kMultiRingNuisSKDetFSI005 ) :
            
    {
      string nutypes = "";
            
      // NEW SKDETFSI
      int ip = ArrIdx::Instance()->MultiRingSKDetFSISystIdx(p);
            
      switch(p) {
                    
        // NEW SKDETFSI
#if defined(__USE_SKDETFSI2013__)
      case (kMultiRingNuisSKDetFSI000) :
      case (kMultiRingNuisSKDetFSI001) :
      case (kMultiRingNuisSKDetFSI002) :
        nutypes = "$\\nu_{\\mu}$/$\\bar{\\nu}_{\\mu}$ CCQE ($MultiR_{\\mu}$)"; break;
      case (kMultiRingNuisSKDetFSI003) :
        nutypes = "$\\nu_{\\mu}$/$\\bar{\\nu}_{\\mu}$ CCnQE ($MultiR_{\\mu}$)"; break;
      case (kMultiRingNuisSKDetFSI004) :
        nutypes = "$\\nu_{e}$/$\\bar{\\nu}_{e}$/sig$\\nu_{e}$ CC ($MultiR_{\\mu}$)"; break;
      case (kMultiRingNuisSKDetFSI005) :
        nutypes = "all NC ($MultiR_{\\mu}$)"; break;
#elif defined(__USE_SKDETFSI2015__)
      case (kMultiRingNuisSKDetFSI000) :
      case (kMultiRingNuisSKDetFSI001) :
      case (kMultiRingNuisSKDetFSI002) :
        nutypes = "$\\nu_{\\mu}$/$\\bar{\\nu}_{\\mu}$ CCQE ($1R_{\\mu}$)"; break;
      case (kMultiRingNuisSKDetFSI003) :
      case (kMultiRingNuisSKDetFSI004) :
      case (kMultiRingNuisSKDetFSI005) :
        nutypes = "$\\nu_{\\mu}$/$\\bar{\\nu}_{\\mu}$ CCnQE ($1R_{\\mu}$)"; break;
      case (kMultiRingNuisSKDetFSI006) :
      case (kMultiRingNuisSKDetFSI007) :
      case (kMultiRingNuisSKDetFSI008) :
        nutypes = "$\\nu_{e}$/$\\bar{\\nu}_{e}$/sig$\\nu_{e}$ CC ($1R_{\\mu}$)"; break;
      case (kMultiRingNuisSKDetFSI009) :
      case (kMultiRingNuisSKDetFSI010) :
      case (kMultiRingNuisSKDetFSI011) :
        nutypes = "all NC ($1R_{\\mu}$)"; break;
#endif
                    
      default : nutypes = "--";
      }
      return std::string(Form("SKDet + FSI/SI %2d; $E_{reco}$( %.2f - %.2f )GeV; %s; MultiRing",ip,
                              kMultiRingSKDetFSIParamsErecoCut[ip][0],kMultiRingSKDetFSIParamsErecoCut[ip][1],nutypes.c_str()));
      break;
    }//FSI for MultiRmu
#endif
#if FSI_BINS != 2
            
  // NEW SKDETFSI
#if defined(__USE_SKDETFSI2013__)
  case ( kMultiRingNuisSKDetFSI006 ) :
  case ( kMultiRingNuisSKDetFSI007 ) :
  case ( kMultiRingNuisSKDetFSI008 ) :
  case ( kMultiRingNuisSKDetFSI009 ) :
  case ( kMultiRingNuisSKDetFSI010 ) :
  case ( kMultiRingNuisSKDetFSI011 ) :
  case ( kMultiRingNuisSKDetFSI012 ) :
  case ( kMultiRingNuisSKDetFSI013 ) :
  case ( kMultiRingNuisSKDetFSI014 ) :
  case ( kMultiRingNuisSKDetFSI015 ) :
  case ( kMultiRingNuisSKDetFSI016 ) :
  case ( kMultiRingNuisSKDetFSI017 ) :
#elif defined(__USE_SKDETFSI2015__)
  case ( kMultiRingNuisSKDetFSI012 ) :
  case ( kMultiRingNuisSKDetFSI013 ) :
  case ( kMultiRingNuisSKDetFSI014 ) :
  case ( kMultiRingNuisSKDetFSI015 ) :
  case ( kMultiRingNuisSKDetFSI016 ) :
  case ( kMultiRingNuisSKDetFSI017 ) :
  case ( kMultiRingNuisSKDetFSI018 ) :
  case ( kMultiRingNuisSKDetFSI019 ) :
  case ( kMultiRingNuisSKDetFSI020 ) :
  case ( kMultiRingNuisSKDetFSI021 ) :
  case ( kMultiRingNuisSKDetFSI022 ) :
  case ( kMultiRingNuisSKDetFSI023 ) :
#endif
            
    {
      string nutypes = "";
            
      // #if FSI_BINS != 3
      //       int ip = (int)p - (int)kMultiRingNuisSKDetFSI000;
      // #else
      //       #if defined(__USE_SKDETFSI2013__)
      //       int ip = (int)p - (int)kMultiRingNuisSKDetFSI006;
      //       #elif defined(__USE_SKDETFSI2015__)
      //       int ip = (int)p - (int)kMultiRingNuisSKDetFSI012;
      //       #endif
      // #endif
            
      // NEW SKDETFSI
      int ip = ArrIdx::Instance()->MultiRingSKDetFSISystIdx(p);
            
      switch(p) {
                    
        // NEW SKDETFSI
#if defined(__USE_SKDETFSI2013__)
      case (kMultiRingNuisSKDetFSI006) :
      case (kMultiRingNuisSKDetFSI007) :
      case (kMultiRingNuisSKDetFSI008) :
        nutypes = "Osc. $\\nu_{e}$ CC ($MultiR_{e}$)"; break;
      case (kMultiRingNuisSKDetFSI009) :
      case (kMultiRingNuisSKDetFSI010) :
      case (kMultiRingNuisSKDetFSI011) :
        nutypes = "$\\nu_{\\mu}$/$\\bar{\\nu}_{\\mu}$ CC ($MultiR_{e}$)"; break;
      case (kMultiRingNuisSKDetFSI012) :
      case (kMultiRingNuisSKDetFSI013) :
      case (kMultiRingNuisSKDetFSI014) :
        nutypes = "$\\nu_{e}$/$\\bar{\\nu}_{e}$ CC ($MultiR_{e}$)"; break;
      case (kMultiRingNuisSKDetFSI015) :
      case (kMultiRingNuisSKDetFSI016) :
      case (kMultiRingNuisSKDetFSI017) :
        nutypes = "all NC ($MultiR_{e}$)"; break;
#elif defined(__USE_SKDETFSI2015__)
      case (kMultiRingNuisSKDetFSI012) :
      case (kMultiRingNuisSKDetFSI013) :
      case (kMultiRingNuisSKDetFSI014) :
        nutypes = "Osc. $\\nu_{e}$ CC ($1R_{e}$)"; break;
      case (kMultiRingNuisSKDetFSI015) :
      case (kMultiRingNuisSKDetFSI016) :
      case (kMultiRingNuisSKDetFSI017) :
        nutypes = "$\\nu_{\\mu}$/$\\bar{\\nu}_{\\mu}$ CC ($1R_{e}$)"; break;
      case (kMultiRingNuisSKDetFSI018) :
      case (kMultiRingNuisSKDetFSI019) :
      case (kMultiRingNuisSKDetFSI020) :
        nutypes = "$\\nu_{e}$/$\\bar{\\nu}_{e}$ CC ($1R_{e}$)"; break;
      case (kMultiRingNuisSKDetFSI021) :
      case (kMultiRingNuisSKDetFSI022) :
      case (kMultiRingNuisSKDetFSI023) :
        nutypes = "all NC ($1R_{e}$)"; break;
#endif
                    
      default : nutypes = "--";
      }
      return std::string(Form("SKDet + FSI/SI %2d; $E_{reco}$( %.2f - %.2f )GeV; %s; MultiRing",ip,
                              kMultiRingSKDetFSIParamsErecoCut[ip][0],kMultiRingSKDetFSIParamsErecoCut[ip][1],nutypes.c_str()));
      break;
    }//FSI for MultiRe
#endif
#else
            
  // OTHER FSI_BINS THAN 1,2,3
            
  case ( kMultiRingNuisSKDetFSI000 ) :
  case ( kMultiRingNuisSKDetFSI001 ) :
  case ( kMultiRingNuisSKDetFSI002 ) :
  case ( kMultiRingNuisSKDetFSI003 ) :
  case ( kMultiRingNuisSKDetFSI004 ) :
  case ( kMultiRingNuisSKDetFSI005 ) :
  case ( kMultiRingNuisSKDetFSI006 ) :
  case ( kMultiRingNuisSKDetFSI007 ) :
  case ( kMultiRingNuisSKDetFSI008 ) :
  case ( kMultiRingNuisSKDetFSI009 ) :
  case ( kMultiRingNuisSKDetFSI010 ) :
  case ( kMultiRingNuisSKDetFSI011 ) :
  case ( kMultiRingNuisSKDetFSI012 ) :
  case ( kMultiRingNuisSKDetFSI013 ) :
  case ( kMultiRingNuisSKDetFSI014 ) :
  case ( kMultiRingNuisSKDetFSI015 ) :
  case ( kMultiRingNuisSKDetFSI016 ) :
  case ( kMultiRingNuisSKDetFSI017 ) :
  case ( kMultiRingNuisSKDetFSI018 ) :
  case ( kMultiRingNuisSKDetFSI019 ) :
    {
      // NEW SKDETFSI
      int ip = ArrIdx::Instance()->MultiRingSKDetFSISystIdx(p);
            
      string nutypes = "";
      switch(p) {
      case (kMultiRingNuisSKDetFSI000) :
      case (kMultiRingNuisSKDetFSI001) :
      case (kMultiRingNuisSKDetFSI002) :
      case (kMultiRingNuisSKDetFSI003) :
      case (kMultiRingNuisSKDetFSI004) :
        nutypes = "Osc. $\\nu_{e}$ CC ($MultiR_{\\mu}$)"; break;
      case (kMultiRingNuisSKDetFSI005) :
      case (kMultiRingNuisSKDetFSI006) :
      case (kMultiRingNuisSKDetFSI007) :
      case (kMultiRingNuisSKDetFSI008) :
      case (kMultiRingNuisSKDetFSI009) :
        nutypes = "$\\nu_{\\mu}$/$\\bar{\\nu}_{\\mu}$ CC ($MultiR_{\\mu}$)"; break;
      case (kMultiRingNuisSKDetFSI010) :
      case (kMultiRingNuisSKDetFSI011) :
      case (kMultiRingNuisSKDetFSI012) :
      case (kMultiRingNuisSKDetFSI013) :
      case (kMultiRingNuisSKDetFSI014) :
        nutypes = "$\\nu_{e}$/$\\bar{\\nu}_{e}$ CC ($MultiR_{\\mu}$)"; break;
      case (kMultiRingNuisSKDetFSI015) :
      case (kMultiRingNuisSKDetFSI016) :
      case (kMultiRingNuisSKDetFSI017) :
      case (kMultiRingNuisSKDetFSI018) :
      case (kMultiRingNuisSKDetFSI019) :
        nutypes = "all NC ($MultiR_{\\mu}$)"; break;
      default : nutypes = "--";
      }
      return std::string(Form("SKDet + FSI/SI %2d; $E_{reco}$( %.2f - %.2f )GeV; %s; MultiRing",ip,
                              kMultiRingSKDetFSIParamsErecoCut[ip][0],kMultiRingSKDetFSIParamsErecoCut[ip][1],nutypes.c_str()));
      break;
    }//FSI for MultiRmu
  case ( kMultiRingNuisSKDetFSI020 ) :
  case ( kMultiRingNuisSKDetFSI021 ) :
  case ( kMultiRingNuisSKDetFSI022 ) :
  case ( kMultiRingNuisSKDetFSI023 ) :
  case ( kMultiRingNuisSKDetFSI024 ) :
  case ( kMultiRingNuisSKDetFSI025 ) :
  case ( kMultiRingNuisSKDetFSI026 ) :
  case ( kMultiRingNuisSKDetFSI027 ) :
  case ( kMultiRingNuisSKDetFSI028 ) :
  case ( kMultiRingNuisSKDetFSI029 ) :
  case ( kMultiRingNuisSKDetFSI030 ) :
  case ( kMultiRingNuisSKDetFSI031 ) :
    {
      // NEW SKDETFSI
      int ip = ArrIdx::Instance()->MultiRingSKDetFSISystIdx(p);
            
      string nutypes = "";
      switch(p) {
      case (kMultiRingNuisSKDetFSI020) :
      case (kMultiRingNuisSKDetFSI021) :
      case (kMultiRingNuisSKDetFSI022) :
        nutypes = "Osc. $\\nu_{e}$ CC ($MultiR_{e}$)"; break;
      case (kMultiRingNuisSKDetFSI023) :
      case (kMultiRingNuisSKDetFSI024) :
      case (kMultiRingNuisSKDetFSI025) :
        nutypes = "$\\nu_{\\mu}$/$\\bar{\\nu}_{\\mu}$ CC ($MultiR_{e}$)"; break;
      case (kMultiRingNuisSKDetFSI026) :
      case (kMultiRingNuisSKDetFSI027) :
      case (kMultiRingNuisSKDetFSI028) :
        nutypes = "$\\nu_{e}$/$\\bar{\\nu}_{e}$ CC ($MultiR_{e}$)"; break;
      case (kMultiRingNuisSKDetFSI029) :
      case (kMultiRingNuisSKDetFSI030) :
      case (kMultiRingNuisSKDetFSI031) :
        nutypes = "all NC ($MultiR_{e}$)"; break;
      default : nutypes = "--";
      }
      return std::string(Form("SKDet + FSI/SI %2d; $E_{reco}$( %.2f - %.2f )GeV; %s; MultiRing",ip,
                              kMultiRingSKDetFSIParamsErecoCut[ip][0],kMultiRingSKDetFSIParamsErecoCut[ip][1],nutypes.c_str()));
      break;
    }
#endif
  //#endif //SKDETFSITOGETHER
            
  ////
            
  default :
    return "--";
  }
  return "--";
}


//________________________________________________________________________________________
string valor::util::AsStringCompact(PhysParamId_t p)
{
    switch(p) {
            
        case ( kPhypMixing23 ) :
        {
            NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
            switch(omp) {
                case(kNuMixSinsq2Theta) : return "sinsq2theta23";  break;
                case(kNuMixSinsqTheta ) : return "sinsqtheta23";   break;
                default                 : return "--";             break;
            }
        }
            break;
        case ( kPhypMixing13 ) :
        {
            NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
            switch(omp) {
                case(kNuMixSinsq2Theta) : return "sinsq2theta13";  break;
                case(kNuMixSinsqTheta ) : return "sinsqtheta13";   break;
                default                 : return "--";             break;
            }
        }
            break;
        case ( kPhypMixing12 ) :
        {
            NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
            switch(omp) {
                case(kNuMixSinsq2Theta) : return "sinsq2theta12";  break;
                case(kNuMixSinsqTheta ) : return "sinsqtheta12";   break;
                default                 : return "--";             break;
            }
        }
            break;
            
        case ( kPhypCPV                   ) : return "delta_cp";            break;
        case ( kPhypDm23sq                ) :
        {
            NuOscDmsqScheme_t dmscheme = Opt::Instance()->GetEffectiveDmsqScheme();
            switch(dmscheme) {
                case(kNuDmsq12andAbs23) : return "|dm32sq|";                  break;
                case(kNuDmsq12andAbs31) : return "|dm31sq|";                  break;
                case(kNuDmsqFogliLisi ) : return "|Dmsq|";                    break;
                default                 : return "--";                        break;
            }
        }
            break;
        case ( kPhypDm12sq                ) :
        {
            NuOscDmsqScheme_t dmscheme = Opt::Instance()->GetEffectiveDmsqScheme();
            switch(dmscheme) {
                case(kNuDmsq12andAbs23) :
                case(kNuDmsq12andAbs31) :
                    return "dm12sq";                  break;
                case(kNuDmsqFogliLisi ) :
                    return "dmsq";                    break;
                default                 :
                    return "--";                        break;
            }
        }
            break;
            
        case ( kPhypMixing23Bar ) :
        {
            NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
            switch(omp) {
                case(kNuMixSinsq2Theta) : return "sinsq2theta23Bar";  break;
                case(kNuMixSinsqTheta ) : return "sinsqtheta23Bar";   break;
                default                 : return "--";                break;
            }
        }
            break;
        case ( kPhypMixing13Bar ) :
        {
            NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
            switch(omp) {
                case(kNuMixSinsq2Theta) : return "sinsq2theta13Bar";  break;
                case(kNuMixSinsqTheta ) : return "sinsqtheta13Bar";   break;
                default                 : return "--";                break;
            }
        }
            break;
        case ( kPhypMixing12Bar ) :
        {
            NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
            switch(omp) {
                case(kNuMixSinsq2Theta) : return "sinsq2theta12Bar";  break;
                case(kNuMixSinsqTheta ) : return "sinsqtheta12Bar";   break;
                default                 : return "--";                break;
            }
        }
            break;
            
        case ( kPhypCPVBar                ) : return "delta_cpbar";            break;
        case ( kPhypDm23sqBar             ) :
        {
            NuOscDmsqScheme_t dmscheme = Opt::Instance()->GetEffectiveDmsqScheme();
            switch(dmscheme) {
                case(kNuDmsq12andAbs23) : return "|dm32sqbar|";                  break;
                case(kNuDmsq12andAbs31) : return "|dm31sqbar|";                  break;
                case(kNuDmsqFogliLisi ) : return "|Dmsqbar|";                    break;
                default                 : return "--";                           break;
            }
        }
            break;
        case ( kPhypDm12sqBar              ) :
        {
            NuOscDmsqScheme_t dmscheme = Opt::Instance()->GetEffectiveDmsqScheme();
            switch(dmscheme) {
                case(kNuDmsq12andAbs23) :
                case(kNuDmsq12andAbs31) :
                    return "dm12sqbar";                  break;
                case(kNuDmsqFogliLisi ) :
                    return "dmsqbar";                    break;
                default                 :
                    return "--";                         break;
            }
        }
            break;
            
            
            // NEW BETABAR
        case ( kPhypBetaBar            ) : return "beta_bar";                 break;
            //
            
            
            
        case ( kPhypAlphaDecay            ) : return "alpha_decay";                 break;
        case ( kPhypSinsq2ThetaDecay      ) : return "sinsq2theta_decay";           break;
        case ( kPhypMusqDecoh             ) : return "musq_decoh";                  break;
        case ( kPhypSinsq2ThetaDecoh      ) : return "sinsq2theta_decoh";           break;
            
            // Steriles
        case ( kPhypMixing14 ) :
        {
            NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
            switch(omp) {
                case(kNuMixSinsq2Theta) : return "sinsq2theta14";  break;
                case(kNuMixSinsqTheta ) : return "sinsqtheta14";   break;
                default                 : return "--";             break;
            }
        }
            break;
        case ( kPhypMixing15 ) :
        {
            NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
            switch(omp) {
                case(kNuMixSinsq2Theta) : return "sinsq2theta15";  break;
                case(kNuMixSinsqTheta ) : return "sinsqtheta15";   break;
                default                 : return "--";             break;
            }
        }
            break;
        case ( kPhypMixing24 ) :
        {
            NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
            switch(omp) {
                case(kNuMixSinsq2Theta) : return "sinsq2theta24";  break;
                case(kNuMixSinsqTheta ) : return "sinsqtheta24";   break;
                default                 : return "--";             break;
            }
        }
            break;
        case ( kPhypMixing25 ) :
        {
            NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
            switch(omp) {
                case(kNuMixSinsq2Theta) : return "sinsq2theta25";  break;
                case(kNuMixSinsqTheta ) : return "sinsqtheta25";   break;
                default                 : return "--";             break;
            }
        }
            break;
        case ( kPhypMixing34 ) :
        {
            NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
            switch(omp) {
                case(kNuMixSinsq2Theta) : return "sinsq2theta34";  break;
                case(kNuMixSinsqTheta ) : return "sinsqtheta34";   break;
                default                 : return "--";             break;
            }
        }
            break;
        case ( kPhypMixing35 ) :
        {
            NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
            switch(omp) {
                case(kNuMixSinsq2Theta) : return "sinsq2theta35";  break;
                case(kNuMixSinsqTheta ) : return "sinsqtheta35";   break;
                default                 : return "--";             break;
            }
        }
            break;
            
        case ( kPhypDelta12               ) : return "delta_12";                   break;
        case ( kPhypDelta13               ) : return "delta_13";                   break;
        case ( kPhypDelta15               ) : return "delta_15";                   break;
        case ( kPhypDelta34               ) : return "delta_34";                   break;
        case ( kPhypDelta35               ) : return "delta_35";                   break;
            
        case ( kPhypDm41sq                ) : return "|dm41sq|";       break;
        case ( kPhypDm51sq                ) : return "|dm51sq|";       break;
            
            // NSI
        case( kPhypNSIepsMuTau) : return "nsi_epsMuTau" ; break;
            
            /*
             * 2011-analysis nuisance parameters
             */
            
            //**************************** FORWARD HORN CURRENT *******************************
            
        case ( kNuisSuperkEnergyScale     ) : return "f_sk_e_scale";    break;
        case ( kNuisBANFF00               ) : return "f_banff00";      break;
        case ( kNuisBANFF01               ) : return "f_banff01";      break;
        case ( kNuisBANFF02               ) : return "f_banff02";      break;
        case ( kNuisBANFF03               ) : return "f_banff03";      break;
        case ( kNuisBANFF04               ) : return "f_banff04";      break;
        case ( kNuisBANFF05               ) : return "f_banff05";      break;
        case ( kNuisBANFF06               ) : return "f_banff06";      break;
        case ( kNuisBANFF07               ) : return "f_banff07";      break;
        case ( kNuisBANFF08               ) : return "f_banff08";      break;
        case ( kNuisBANFF09               ) : return "f_banff09";      break;
        case ( kNuisBANFF10               ) : return "f_banff10";      break;
        case ( kNuisBANFF11               ) : return "f_banff11";      break;
        case ( kNuisBANFF12               ) : return "f_banff12";      break;
        case ( kNuisBANFF13               ) : return "f_banff13";      break;
        case ( kNuisBANFF14               ) : return "f_banff14";      break;
        case ( kNuisBANFF15               ) : return "f_banff15";      break;
#ifdef __USE_ALLFLUX__
        case ( kNuisBANFF16               ) : return "f_banff16";      break;
        case ( kNuisBANFF17               ) : return "f_banff17";      break;
        case ( kNuisBANFF18               ) : return "f_banff18";      break;
        case ( kNuisBANFF19               ) : return "f_banff19";      break;
        case ( kNuisBANFF20               ) : return "f_banff20";      break;
        case ( kNuisBANFF21               ) : return "f_banff21";      break;
        case ( kNuisBANFF22               ) : return "f_banff22";      break;
        case ( kNuisBANFF23               ) : return "f_banff23";      break;
        case ( kNuisBANFF24               ) : return "f_banff24";      break;
#endif
            
            
            // NEW FHC+RHC JOINT (NEW)
#if defined(__FHC_RUNNING__) && defined(__RHC_RUNNING__)
            
#ifndef ONEOVERMAQE
        case ( kNuisBANFFMAQE             ) : return "f_banffmaqe_fhc_rhc";    break;
#else
        case ( kNuisBANFFMAQE             ) : return "f_banff1overmaqe2_fhc_rhc";  break;
#endif
        case ( kNuisBANFFMARES            ) : return "f_banffmares_fhc_rhc";   break;
        case ( kNuisBANFFCCQEE0           ) : return "f_banffccqee0_fhc_rhc";  break;
#if BANFF_BINS != 3
        case ( kNuisBANFFCCQEE1           ) : return "f_banffccqee1_fhc_rhc";  break;
        case ( kNuisBANFFCCQEE2           ) : return "f_banffccqee2_fhc_rhc";  break;
#endif
        case ( kNuisBANFFCC1piE0          ) : return "f_banffcc1pie0_fhc_rhc"; break;
#if BANFF_BINS != 3
        case ( kNuisBANFFCC1piE1          ) : return "f_banffcc1pie1_fhc_rhc"; break;
#endif
#if (BANFF_BINS != 2)
        case ( kNuisBANFFNC1pi0           ) : return "f_banffnc1pi0_fhc_rhc";  break;
#endif
            
        case ( kNuisBANFFNorm2p2h          ) : return "f_banff2p2h_fhc_rhc";            break;
        case ( kNuisBANFFCA5               ) : return "f_banffca5_fhc_rhc";            break;
        case ( kNuisBANFFBgRES             ) : return "f_banffbgres_fhc_rhc";          break;
#if (BANFF_BINS != 2)
        case ( kNuisBANFFSCCA              ) : return "f_banffscca_fhc_rhc";           break;
        case ( kNuisBANFFSCCV              ) : return "f_banffsccv_fhc_rhc";           break;
#endif
        case ( kNuisBANFFFermiMom          ) : return "f_banffpfermi_fhc_rhc";         break; 
        case ( kNuisBANFFShapeCCoth        ) : return "f_banffshapeCCoth_fhc_rhc";     break; 
          //case ( kNuisBANFFBindE             ) : return "f_banffbindE_fhc_rhc";          break; 
        case ( kNuisBANFFNormCCcoh         ) : return "f_banffnorm_cccoh_fhc_rhc";     break; 
        case ( kNuisBANFFNormNCcoh         ) : return "f_banffnorm_nccoh_fhc_rhc";     break; 
        case ( kNuisBANFFNormNCoth         ) : return "f_banffnorm_ncoth_fhc_rhc";     break; 
        case ( kNuisBANFFNormNuEToNuMu     ) : return "f_banffnorm_nuetonumu_fhc_rhc"; break; 
        case ( kNuisBANFFNormNC1Gamma      ) : return "f_banffnorm_nc1gamma_fhc_rhc";  break; 
        case ( kNuisBANFFNormNuEBarToNuMuBar) : return "f_banffnorm_nuebartonumubar_fhc_rhc"; break; 
        case ( kNuisBANFFNorm2p2hBar       ) : return "f_banff2p2hbar_fhc_rhc";         break;
        case ( kNuisBANFFBeRPA_A           ) : return "f_banffshape_berpa_A_fhc_rhc";    break;
        case ( kNuisBANFFBeRPA_B           ) : return "f_banffshape_berpa_B_fhc_rhc";    break;
        case ( kNuisBANFFBeRPA_D           ) : return "f_banffshape_berpa_D_fhc_rhc";    break;
        case ( kNuisBANFFBeRPA_E           ) : return "f_banffshape_berpa_E_fhc_rhc";    break;
        case ( kNuisBANFFBeRPA_U           ) : return "f_banffshape_berpa_U_fhc_rhc";    break;
        case ( kNuisBANFFShape2p2h         ) : return "f_banffshape_2p2hnu_fhc_rhc";     break;
        case ( kNuisBANFFNorm2p2hCtoO      ) : return "f_banffnorm_2p2hCtoO_fhc_rhc";  break;
        case ( kNuisBANFFEb      ) : return "f_banffshape_eb_fhc_rhc";
        case ( kNuisBANFFMultiReNorm      ) : return "f_banffnorm_MultiRe_fhc";

#else
#ifndef ONEOVERMAQE
        case ( kNuisBANFFMAQE              ) : return "f_banffmaqe";    break;
#else
        case ( kNuisBANFFMAQE              ) : return "f_banff1overmaqe2";  break;
#endif
        case ( kNuisBANFFMARES             ) : return "f_banffmares";   break;
        case ( kNuisBANFFCCQEE0            ) : return "f_banffccqee0";  break;
#if BANFF_BINS != 3
        case ( kNuisBANFFCCQEE1            ) : return "f_banffccqee1";  break;
        case ( kNuisBANFFCCQEE2            ) : return "f_banffccqee2";  break;
#endif
        case ( kNuisBANFFCC1piE0           ) : return "f_banffcc1pie0"; break;
#if BANFF_BINS != 3
        case ( kNuisBANFFCC1piE1           ) : return "f_banffcc1pie1"; break;
#endif
#if (BANFF_BINS != 2)
        case ( kNuisBANFFNC1pi0            ) : return "f_banffnc1pi0";  break;
#endif
        case ( kNuisBANFFNorm2p2h          ) : return "f_banff2p2h";           break;
        case ( kNuisBANFFCA5               ) : return "f_banffca5";            break;
        case ( kNuisBANFFBgRES             ) : return "f_banffbgres";          break;
#if (BANFF_BINS != 2)
        case ( kNuisBANFFSCCA              ) : return "f_banffscca";           break;
        case ( kNuisBANFFSCCV              ) : return "f_banffsccv";           break;
#endif
        case ( kNuisBANFFFermiMom          ) : return "f_banffpfermi";         break; 
        case ( kNuisBANFFShapeCCoth        ) : return "f_banffshapeCCoth";     break; 
          //case ( kNuisBANFFBindE             ) : return "f_banffbindE";          break; 
        case ( kNuisBANFFNormCCcoh         ) : return "f_banffnorm_cccoh";     break; 
        case ( kNuisBANFFNormNCcoh         ) : return "f_banffnorm_nccoh";     break; 
        case ( kNuisBANFFNormNCoth         ) : return "f_banffnorm_ncoth";     break; 
        case ( kNuisBANFFNormNuEToNuMu     ) : return "f_banffnorm_nuetonumu"; break; 
        case ( kNuisBANFFNormNC1Gamma      ) : return "f_banffnorm_nc1gamma";  break; 
        case ( kNuisBANFFNormNuEBarToNuMuBar) : return "f_banffnorm_nuebartonumubar"; break; 
        case ( kNuisBANFFNorm2p2hBar       ) : return "f_banff2p2hbar";          break;
        case ( kNuisBANFFBeRPA_A           ) : return "f_banffshape_berpa_A";    break;
        case ( kNuisBANFFBeRPA_B           ) : return "f_banffshape_berpa_B";    break;
        case ( kNuisBANFFBeRPA_D           ) : return "f_banffshape_berpa_D";    break;
        case ( kNuisBANFFBeRPA_E           ) : return "f_banffshape_berpa_E";    break;
        case ( kNuisBANFFBeRPA_U           ) : return "f_banffshape_berpa_U";    break;
        case ( kNuisBANFFShape2p2h         ) : return "f_banffshape_2p2hnu";     break;
        case ( kNuisBANFFNorm2p2hCtoO      ) : return "f_banffnorm_2p2hCtoO";    break;
        case ( kNuisBANFFEb         ) : return "f_banffshape_eb";
        case ( kNuisBANFFMultiReNorm       ) : return "f_banffnorm_MultiRe";
            //
            
#endif //if defined(__FHC_RUNNING__) && defined(__RHC_RUNNING__)
            
            // NEW XSEC
        case ( kNuisXSecRPA               ) : return "f_rpa";          break;
            
        case ( kNuisXSecFermiMom          ) : return "f_pfermi";       break;
        case ( kNuisXSecSpecFunc          ) : return "f_specfunc";     break;
        case ( kNuisXSecWshape            ) : return "f_wshape";       break;
        case ( kNuisXSecShapeCCoth        ) : return "f_shapeCCoth";   break;
        case ( kNuisXSecPiLessDelta       ) : return "f_pilessdelta";  break;
#ifndef __NUEONLY__
        case ( kNuisXSecBindE             ) : return "f_bindE";        break;
#endif
        case ( kNuisXSecNormCCcoh         ) : return "f_norm_cccoh";   break;
#ifndef __NUEONLY__
        case ( kNuisXSecNormNC1piC        ) : return "f_norm_nc1piC";  break;
#endif
#ifndef __NUMUONLY__
        case ( kNuisXSecNormNCcoh         ) : return "f_norm_nccoh";   break;
#endif
        case ( kNuisXSecNormNCoth         ) : return "f_norm_ncoth";   break;
        case ( kNuisXSecNormNuEToNuMu     ) : return "f_norm_nuetonumu";  break;
#ifndef __NUEONLY__
        case ( kNuisXSecNormAntiNuToNu    ) : return "f_norm_antinutonu"; break;
#endif
            //#ifdef SKFSITOGETHER
#if FSI_BINS != 3
        case ( kNuisSKDetFSI000 ) : return "f_skdetfsi000"; break;
        case ( kNuisSKDetFSI001 ) : return "f_skdetfsi001"; break;
        case ( kNuisSKDetFSI002 ) : return "f_skdetfsi002"; break;
        case ( kNuisSKDetFSI003 ) : return "f_skdetfsi003"; break;
        case ( kNuisSKDetFSI004 ) : return "f_skdetfsi004"; break;
        case ( kNuisSKDetFSI005 ) : return "f_skdetfsi005"; break;
            
            // NEW SKDETFSI
#if defined(__USE_SKDETFSI2015__)
        case ( kNuisSKDetFSI006 ) : return "f_skdetfsi006"; break;
        case ( kNuisSKDetFSI007 ) : return "f_skdetfsi007"; break;
        case ( kNuisSKDetFSI008 ) : return "f_skdetfsi008"; break;
        case ( kNuisSKDetFSI009 ) : return "f_skdetfsi009"; break;
        case ( kNuisSKDetFSI010 ) : return "f_skdetfsi010"; break;
        case ( kNuisSKDetFSI011 ) : return "f_skdetfsi011"; break;
#endif
            
#endif
#if FSI_BINS != 2
            
            // NEW SKDETFSI
#if defined(__USE_SKDETFSI2013__)
        case ( kNuisSKDetFSI006 ) : return "f_skdetfsi006"; break;
        case ( kNuisSKDetFSI007 ) : return "f_skdetfsi007"; break;
        case ( kNuisSKDetFSI008 ) : return "f_skdetfsi008"; break;
        case ( kNuisSKDetFSI009 ) : return "f_skdetfsi009"; break;
        case ( kNuisSKDetFSI010 ) : return "f_skdetfsi010"; break;
        case ( kNuisSKDetFSI011 ) : return "f_skdetfsi011"; break;
        case ( kNuisSKDetFSI012 ) : return "f_skdetfsi012"; break;
        case ( kNuisSKDetFSI013 ) : return "f_skdetfsi013"; break;
        case ( kNuisSKDetFSI014 ) : return "f_skdetfsi014"; break;
        case ( kNuisSKDetFSI015 ) : return "f_skdetfsi015"; break;
        case ( kNuisSKDetFSI016 ) : return "f_skdetfsi016"; break;
        case ( kNuisSKDetFSI017 ) : return "f_skdetfsi017"; break;
#elif defined(__USE_SKDETFSI2013__)
        case ( kNuisSKDetFSI012 ) : return "f_skdetfsi012"; break;
        case ( kNuisSKDetFSI013 ) : return "f_skdetfsi013"; break;
        case ( kNuisSKDetFSI014 ) : return "f_skdetfsi014"; break;
        case ( kNuisSKDetFSI015 ) : return "f_skdetfsi015"; break;
        case ( kNuisSKDetFSI016 ) : return "f_skdetfsi016"; break;
        case ( kNuisSKDetFSI017 ) : return "f_skdetfsi017"; break;
        case ( kNuisSKDetFSI018 ) : return "f_skdetfsi018"; break;
        case ( kNuisSKDetFSI019 ) : return "f_skdetfsi019"; break;
        case ( kNuisSKDetFSI020 ) : return "f_skdetfsi020"; break;
        case ( kNuisSKDetFSI021 ) : return "f_skdetfsi021"; break;
        case ( kNuisSKDetFSI022 ) : return "f_skdetfsi022"; break;
        case ( kNuisSKDetFSI023 ) : return "f_skdetfsi023"; break;
#endif
#endif
            // NEW SKDETFSI (start indices from 24 instead of 18)
#ifdef SKCORRELATIONS
        case ( kNuisSKDetFSI024 ) : return "f_skdetfsi024"; break;
        case ( kNuisSKDetFSI025 ) : return "f_skdetfsi025"; break;
        case ( kNuisSKDetFSI026 ) : return "f_skdetfsi026"; break;
        case ( kNuisSKDetFSI027 ) : return "f_skdetfsi027"; break;
        case ( kNuisSKDetFSI028 ) : return "f_skdetfsi028"; break;
        case ( kNuisSKDetFSI029 ) : return "f_skdetfsi029"; break;
        case ( kNuisSKDetFSI030 ) : return "f_skdetfsi030"; break;
        case ( kNuisSKDetFSI031 ) : return "f_skdetfsi031"; break;
        case ( kNuisSKDetFSI032 ) : return "f_skdetfsi032"; break;
        case ( kNuisSKDetFSI033 ) : return "f_skdetfsi033"; break;
        case ( kNuisSKDetFSI034 ) : return "f_skdetfsi034"; break;
        case ( kNuisSKDetFSI035 ) : return "f_skdetfsi035"; break;
        case ( kNuisSKDetFSI036 ) : return "f_skdetfsi036"; break;
        case ( kNuisSKDetFSI037 ) : return "f_skdetfsi037"; break;
#endif
            //#endif //SKDETFSITOGETHER
            
        case ( kNuisRHCToFHC ) : return "f_rhcfhcnorm"; break;
            
            //***************************** REVERSE HORN CURRENT *************************
            
        case ( kRHCNuisSuperkEnergyScale     ) : return "f_sk_e_scale_rhc";   break;
            
        case ( kRHCNuisBANFF00               ) : return "f_banff00_rhc";      break;
        case ( kRHCNuisBANFF01               ) : return "f_banff01_rhc";      break;
        case ( kRHCNuisBANFF02               ) : return "f_banff02_rhc";      break;
        case ( kRHCNuisBANFF03               ) : return "f_banff03_rhc";      break;
        case ( kRHCNuisBANFF04               ) : return "f_banff04_rhc";      break;
        case ( kRHCNuisBANFF05               ) : return "f_banff05_rhc";      break;
        case ( kRHCNuisBANFF06               ) : return "f_banff06_rhc";      break;
        case ( kRHCNuisBANFF07               ) : return "f_banff07_rhc";      break;
        case ( kRHCNuisBANFF08               ) : return "f_banff08_rhc";      break;
        case ( kRHCNuisBANFF09               ) : return "f_banff09_rhc";      break;
        case ( kRHCNuisBANFF10               ) : return "f_banff10_rhc";      break;
        case ( kRHCNuisBANFF11               ) : return "f_banff11_rhc";      break;
        case ( kRHCNuisBANFF12               ) : return "f_banff12_rhc";      break;
        case ( kRHCNuisBANFF13               ) : return "f_banff13_rhc";      break;
        case ( kRHCNuisBANFF14               ) : return "f_banff14_rhc";      break;
        case ( kRHCNuisBANFF15               ) : return "f_banff15_rhc";      break;
#ifdef __USE_ALLFLUX__
        case ( kRHCNuisBANFF16               ) : return "f_banff16_rhc";      break;
        case ( kRHCNuisBANFF17               ) : return "f_banff17_rhc";      break;
        case ( kRHCNuisBANFF18               ) : return "f_banff18_rhc";      break;
        case ( kRHCNuisBANFF19               ) : return "f_banff19_rhc";      break;
        case ( kRHCNuisBANFF20               ) : return "f_banff20_rhc";      break;
        case ( kRHCNuisBANFF21               ) : return "f_banff21_rhc";      break;
        case ( kRHCNuisBANFF22               ) : return "f_banff22_rhc";      break;
        case ( kRHCNuisBANFF23               ) : return "f_banff23_rhc";      break;
        case ( kRHCNuisBANFF24               ) : return "f_banff24_rhc";      break;
#endif
            
#ifndef ONEOVERMAQE
        case ( kRHCNuisBANFFMAQE             ) : return "f_banffmaqe_rhc";    break;
#else
        case ( kRHCNuisBANFFMAQE             ) : return "f_banff1overmaqe2_rhc";  break;
#endif
        case ( kRHCNuisBANFFMARES            ) : return "f_banffmares_rhc";   break;
        case ( kRHCNuisBANFFCCQEE0           ) : return "f_banffccqee0_rhc";  break;
#if BANFF_BINS != 3
        case ( kRHCNuisBANFFCCQEE1           ) : return "f_banffccqee1_rhc";  break;
        case ( kRHCNuisBANFFCCQEE2           ) : return "f_banffccqee2_rhc";  break;
#endif
        case ( kRHCNuisBANFFCC1piE0          ) : return "f_banffcc1pie0_rhc"; break;
#if BANFF_BINS != 3
        case ( kRHCNuisBANFFCC1piE1          ) : return "f_banffcc1pie1_rhc"; break;
#endif
#if (BANFF_BINS != 2)
        case ( kRHCNuisBANFFNC1pi0           ) : return "f_banffnc1pi0_rhc";  break;
#endif
            
            // NEW XSEC
        case ( kRHCNuisBANFF2p2h               ) : return "f_banff2p2h_rhc";            break;
        case ( kRHCNuisBANFFCA5               ) : return "f_banffca5_rhc";            break;
        case ( kRHCNuisBANFFBgRES             ) : return "f_banffbgres_rhc";          break;
            // NEW SPLIT ANA
#if (BANFF_BINS != 2)
        case ( kRHCNuisBANFFSCCA              ) : return "f_banffscca_rhc";           break;
        case ( kRHCNuisBANFFSCCV              ) : return "f_banffsccv_rhc";           break;
#endif
        case ( kRHCNuisBANFFFermiMom          ) : return "f_banffpfermi_rhc";         break;
        case ( kRHCNuisBANFFShapeCCoth        ) : return "f_banffshapeCCoth_rhc";     break;
          //case ( kRHCNuisBANFFBindE             ) : return "f_banffbindE_rhc";          break;
        case ( kRHCNuisBANFFNormCCcoh         ) : return "f_banffnorm_cccoh_rhc";     break;
        case ( kRHCNuisBANFFNormNCcoh         ) : return "f_banffnorm_nccoh_rhc";     break;
        case ( kRHCNuisBANFFNormNCoth         ) : return "f_banffnorm_ncoth_rhc";     break;
        case ( kRHCNuisBANFFNormNuEToNuMu     ) : return "f_banffnorm_nuetonumu_rhc"; break;
        case ( kRHCNuisBANFFNormNC1Gamma     ) : return "f_banffnorm_nc1gamma_rhc"; break; // uncorrelated in the past
        case ( kRHCNuisBANFFNormNuEBarToNuMuBar) : return "f_banffnorm_nuebartonumubar_rhc"; break; // uncorrelated in the past
        case ( kRHCNuisBANFF2p2hBar            ) : return "f_banff2p2hbar_rhc";            break;
            
            //
            
            // NEW XSEC
        case ( kRHCNuisXSecRPA               ) : return "f_rpa_rhc";          break;
            
        case ( kRHCNuisXSecFermiMom          ) : return "f_pfermi_rhc";       break;
        case ( kRHCNuisXSecSpecFunc          ) : return "f_specfunc_rhc";     break;
        case ( kRHCNuisXSecWshape            ) : return "f_wshape_rhc";       break;
        case ( kRHCNuisXSecShapeCCoth        ) : return "f_shapeCCoth_rhc";   break;
        case ( kRHCNuisXSecPiLessDelta       ) : return "f_pilessdelta_rhc";  break;
#ifndef __NUEONLY__
        case ( kRHCNuisXSecBindE             ) : return "f_bindE_rhc";        break;
#endif
        case ( kRHCNuisXSecNormCCcoh         ) : return "f_norm_cccoh_rhc";   break;
#ifndef __NUEONLY__
        case ( kRHCNuisXSecNormNC1piC        ) : return "f_norm_nc1piC_rhc";  break;
#endif
#ifndef __NUMUONLY__
        case ( kRHCNuisXSecNormNCcoh         ) : return "f_norm_nccoh_rhc";   break;
#endif
        case ( kRHCNuisXSecNormNCoth         ) : return "f_norm_ncoth_rhc";   break;
        case ( kRHCNuisXSecNormNuEToNuMu     ) : return "f_norm_nuetonumu_rhc";  break;
#ifndef __NUEONLY__
        case ( kRHCNuisXSecNormAntiNuToNu    ) : return "f_norm_antinutonu_rhc"; break;
#endif
            //#ifdef SKFSITOGETHER
#if FSI_BINS != 3
        case ( kRHCNuisSKDetFSI000 ) : return "f_skdetfsi000_rhc"; break;
        case ( kRHCNuisSKDetFSI001 ) : return "f_skdetfsi001_rhc"; break;
        case ( kRHCNuisSKDetFSI002 ) : return "f_skdetfsi002_rhc"; break;
        case ( kRHCNuisSKDetFSI003 ) : return "f_skdetfsi003_rhc"; break;
        case ( kRHCNuisSKDetFSI004 ) : return "f_skdetfsi004_rhc"; break;
        case ( kRHCNuisSKDetFSI005 ) : return "f_skdetfsi005_rhc"; break;
            // NEW SKDETFSI
#if defined(__USE_SKDETFSI2015__)
        case ( kRHCNuisSKDetFSI006 ) : return "f_skdetfsi006_rhc"; break;
        case ( kRHCNuisSKDetFSI007 ) : return "f_skdetfsi007_rhc"; break;
        case ( kRHCNuisSKDetFSI008 ) : return "f_skdetfsi008_rhc"; break;
        case ( kRHCNuisSKDetFSI009 ) : return "f_skdetfsi009_rhc"; break;
        case ( kRHCNuisSKDetFSI010 ) : return "f_skdetfsi010_rhc"; break;
        case ( kRHCNuisSKDetFSI011 ) : return "f_skdetfsi011_rhc"; break;
#endif
            
#endif
#if FSI_BINS != 2
            
            // NEW SKDETFSI
#if defined(__USE_SKDETFSI2013__)
        case ( kRHCNuisSKDetFSI006 ) : return "f_skdetfsi006_rhc"; break;
        case ( kRHCNuisSKDetFSI007 ) : return "f_skdetfsi007_rhc"; break;
        case ( kRHCNuisSKDetFSI008 ) : return "f_skdetfsi008_rhc"; break;
        case ( kRHCNuisSKDetFSI009 ) : return "f_skdetfsi009_rhc"; break;
        case ( kRHCNuisSKDetFSI010 ) : return "f_skdetfsi010_rhc"; break;
        case ( kRHCNuisSKDetFSI011 ) : return "f_skdetfsi011_rhc"; break;
        case ( kRHCNuisSKDetFSI012 ) : return "f_skdetfsi012_rhc"; break;
        case ( kRHCNuisSKDetFSI013 ) : return "f_skdetfsi013_rhc"; break;
        case ( kRHCNuisSKDetFSI014 ) : return "f_skdetfsi014_rhc"; break;
        case ( kRHCNuisSKDetFSI015 ) : return "f_skdetfsi015_rhc"; break;
        case ( kRHCNuisSKDetFSI016 ) : return "f_skdetfsi016_rhc"; break;
        case ( kRHCNuisSKDetFSI017 ) : return "f_skdetfsi017_rhc"; break;
#elif defined(__USE_SKDETFSI2013__)
        case ( kRHCNuisSKDetFSI012 ) : return "f_skdetfsi012_rhc"; break;
        case ( kRHCNuisSKDetFSI013 ) : return "f_skdetfsi013_rhc"; break;
        case ( kRHCNuisSKDetFSI014 ) : return "f_skdetfsi014_rhc"; break;
        case ( kRHCNuisSKDetFSI015 ) : return "f_skdetfsi015_rhc"; break;
        case ( kRHCNuisSKDetFSI016 ) : return "f_skdetfsi016_rhc"; break;
        case ( kRHCNuisSKDetFSI017 ) : return "f_skdetfsi017_rhc"; break;
        case ( kRHCNuisSKDetFSI018 ) : return "f_skdetfsi018_rhc"; break;
        case ( kRHCNuisSKDetFSI019 ) : return "f_skdetfsi019_rhc"; break;
        case ( kRHCNuisSKDetFSI020 ) : return "f_skdetfsi020_rhc"; break;
        case ( kRHCNuisSKDetFSI021 ) : return "f_skdetfsi021_rhc"; break;
        case ( kRHCNuisSKDetFSI022 ) : return "f_skdetfsi022_rhc"; break;
        case ( kRHCNuisSKDetFSI023 ) : return "f_skdetfsi023_rhc"; break;
#endif
#endif
            
            // NEW SKDETFSI (start indices from 24 instead of 18)
#ifdef SKCORRELATIONS
        case ( kRHCNuisSKDetFSI024 ) : return "f_skdetfsi024_rhc"; break;
        case ( kRHCNuisSKDetFSI025 ) : return "f_skdetfsi025_rhc"; break;
        case ( kRHCNuisSKDetFSI026 ) : return "f_skdetfsi026_rhc"; break;
        case ( kRHCNuisSKDetFSI027 ) : return "f_skdetfsi027_rhc"; break;
        case ( kRHCNuisSKDetFSI028 ) : return "f_skdetfsi028_rhc"; break;
        case ( kRHCNuisSKDetFSI029 ) : return "f_skdetfsi029_rhc"; break;
        case ( kRHCNuisSKDetFSI030 ) : return "f_skdetfsi030_rhc"; break;
        case ( kRHCNuisSKDetFSI031 ) : return "f_skdetfsi031_rhc"; break;
        case ( kRHCNuisSKDetFSI032 ) : return "f_skdetfsi032_rhc"; break;
        case ( kRHCNuisSKDetFSI033 ) : return "f_skdetfsi033_rhc"; break;
        case ( kRHCNuisSKDetFSI034 ) : return "f_skdetfsi034_rhc"; break;
        case ( kRHCNuisSKDetFSI035 ) : return "f_skdetfsi035_rhc"; break;
        case ( kRHCNuisSKDetFSI036 ) : return "f_skdetfsi036_rhc"; break;
        case ( kRHCNuisSKDetFSI037 ) : return "f_skdetfsi037_rhc"; break;
#endif
            //#endif //SKDETFSITOGETHER
            
            //
            
            
            //***************************** MULTIRING FORWARD HORN CURRENT *************************
            
            /// NEW MULTIRING 3
        case ( kMultiRingNuisSuperkEnergyScale     ) : return "f_sk_e_scale_multiring";      break;
            ///
            
            /// NEW MULTIRING 2
            
            //#ifdef SKFSITOGETHER
#if MR_FSI_BINS != 3
        case ( kMultiRingNuisSKDetFSI000 ) : return "f_skdetfsi000_multiring"; break;
        case ( kMultiRingNuisSKDetFSI001 ) : return "f_skdetfsi001_multiring"; break;
        case ( kMultiRingNuisSKDetFSI002 ) : return "f_skdetfsi002_multiring"; break;
        case ( kMultiRingNuisSKDetFSI003 ) : return "f_skdetfsi003_multiring"; break;
        case ( kMultiRingNuisSKDetFSI004 ) : return "f_skdetfsi004_multiring"; break;
        case ( kMultiRingNuisSKDetFSI005 ) : return "f_skdetfsi005_multiring"; break;            
            
#endif
#if MR_FSI_BINS != 2
            
            // NEW SKDETFSI
#if defined(__USE_SKDETFSI2013__)
        case ( kMultiRingNuisSKDetFSI006 ) : return "f_skdetfsi006_multiring"; break;
        case ( kMultiRingNuisSKDetFSI007 ) : return "f_skdetfsi007_multiring"; break;
        case ( kMultiRingNuisSKDetFSI008 ) : return "f_skdetfsi008_multiring"; break;
        case ( kMultiRingNuisSKDetFSI009 ) : return "f_skdetfsi009_multiring"; break;
        case ( kMultiRingNuisSKDetFSI010 ) : return "f_skdetfsi010_multiring"; break;
        case ( kMultiRingNuisSKDetFSI011 ) : return "f_skdetfsi011_multiring"; break;
        case ( kMultiRingNuisSKDetFSI012 ) : return "f_skdetfsi012_multiring"; break;
        case ( kMultiRingNuisSKDetFSI013 ) : return "f_skdetfsi013_multiring"; break;
        case ( kMultiRingNuisSKDetFSI014 ) : return "f_skdetfsi014_multiring"; break;
        case ( kMultiRingNuisSKDetFSI015 ) : return "f_skdetfsi015_multiring"; break;
        case ( kMultiRingNuisSKDetFSI016 ) : return "f_skdetfsi016_multiring"; break;
        case ( kMultiRingNuisSKDetFSI017 ) : return "f_skdetfsi017_multiring"; break;
#endif
#endif            
            
            // NEW SKDETFSI (start indices from 24 instead of 18)
#ifdef SKCORRELATIONS
        case ( kMultiRingNuisSKDetFSI024 ) : return "f_skdetfsi024_multiring"; break;
        case ( kMultiRingNuisSKDetFSI025 ) : return "f_skdetfsi025_multiring"; break;
        case ( kMultiRingNuisSKDetFSI026 ) : return "f_skdetfsi026_multiring"; break;
        case ( kMultiRingNuisSKDetFSI027 ) : return "f_skdetfsi027_multiring"; break;
        case ( kMultiRingNuisSKDetFSI028 ) : return "f_skdetfsi028_multiring"; break;
        case ( kMultiRingNuisSKDetFSI029 ) : return "f_skdetfsi029_multiring"; break;
        case ( kMultiRingNuisSKDetFSI030 ) : return "f_skdetfsi030_multiring"; break;
        case ( kMultiRingNuisSKDetFSI031 ) : return "f_skdetfsi031_multiring"; break;
        case ( kMultiRingNuisSKDetFSI032 ) : return "f_skdetfsi032_multiring"; break;
        case ( kMultiRingNuisSKDetFSI033 ) : return "f_skdetfsi033_multiring"; break;
        case ( kMultiRingNuisSKDetFSI034 ) : return "f_skdetfsi034_multiring"; break;
        case ( kMultiRingNuisSKDetFSI035 ) : return "f_skdetfsi035_multiring"; break;
        case ( kMultiRingNuisSKDetFSI036 ) : return "f_skdetfsi036_multiring"; break;
        case ( kMultiRingNuisSKDetFSI037 ) : return "f_skdetfsi037_multiring"; break;
#endif
            //#endif //SKDETFSITOGETHER
            
            ///
            
        default :
            return "--";
    }
    return "--";
}

//____________________________________________________________________________


string valor::util::AsStringPaperTex(PhysParamId_t p)
{
    switch(p) {
            
        case ( kPhypMixing23 ) :
        {
            NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
            switch(omp) {
                case(kNuMixSinsq2Theta) : return "$\\sin^{2}(2\\theta_{23})$"; break;
                case(kNuMixSinsqTheta ) : return "$\\sin^{2}(\\theta_{23})$";   break;
                default                 : return "--";                break;
            }
        }
            break;
        case ( kPhypMixing13 ) :
        {
            NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
            switch(omp) {
                case(kNuMixSinsq2Theta) : return "$\\sin^{2}(2\\theta_{13})$";  break;
                case(kNuMixSinsqTheta ) : return "$\\sin^{2}(\\theta_{13})$";    break;
                default                 : return "--";                 break;
            }
        }
            break;
        case ( kPhypMixing12 ) :
        {
            NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
            switch(omp) {
                case(kNuMixSinsq2Theta) : return "$\\sin^{2}(2\\theta_{12})$";  break;
                case(kNuMixSinsqTheta ) : return "$\\sin^{2}(\\theta_{12})$";    break;
                default                 : return "--";                 break;
            }
        }
            break;
            
        case ( kPhypCPV ) :
        {
            conv::NuOscDeltaCP_t cpv = Opt::Instance()->gOscDeltaCP;
            switch(cpv) {
                case(kNuDeltaCPPhase) : return "$\\delta_{CP}$";       break;
                case(kNuDeltaCPSin  ) : return "$\\sin(\\delta_{CP})$";break;
                default               : return "--";                   break;
            }
        }
        case ( kPhypDm23sq                ) :
        {
            NuOscDmsqScheme_t omp = Opt::Instance()->GetEffectiveDmsqScheme();
            switch(omp) {
                case(kNuDmsqFogliLisi) :        return "$|\\Delta m^{2}| = |m^{2}_{3} - (m^{2}_{1} + m^{2}_{1})/2|$";                  break;
                case(kNuDmsq12andAbs23):        return "$|\\Delta m^{2}_{32}|$";                  break;
                case(kNuDmsq12andAbs31):        return "$|\\Delta m^{2}_{31}|$";                  break;
                default                :        return "--";                                   break;
            }
        }
            break;
            
        case ( kPhypDm12sq                ) :
        {
            NuOscDmsqScheme_t omp = Opt::Instance()->GetEffectiveDmsqScheme();
            switch(omp) {
                case(kNuDmsqFogliLisi):
                    return "$\\deltam^{2} = \\Delta m^{2}_{21}$";      break;
                case(kNuDmsq12andAbs23):
                case(kNuDmsq12andAbs31):
                    return "$\\Delta m^{2}_{21}$";                   break;
                default:
                    return "--";                                   break;
            }
        }
        case ( kPhypMixing23Bar ) :
        {
            NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
            switch(omp) {
                case(kNuMixSinsq2Theta) : return "$\\sin^{2}(2\\bar{\\theta_{23}})$"; break;
                case(kNuMixSinsqTheta ) : return "$\\sin^{2}(\\bar{\\theta_{23}})$";   break;
                default                 : return "--";                break;
            }
        }
            break;
        case ( kPhypMixing13Bar ) :
        {
            NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
            switch(omp) {
                case(kNuMixSinsq2Theta) : return "$\\sin^{2}(2\\bar{\\theta_{13}})$";  break;
                case(kNuMixSinsqTheta ) : return "$\\sin^{2}(\\bar{\\theta_{13}})$";    break;
                default                 : return "--";                 break;
            }
        }
            break;
        case ( kPhypMixing12Bar ) :
        {
            NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
            switch(omp) {
                case(kNuMixSinsq2Theta) : return "$\\sin^{2}(2\\bar{\\theta_{12}})$";  break;
                case(kNuMixSinsqTheta ) : return "$\\sin^{2}(\\bar{\\theta_{12}})$";    break;
                default                 : return "--";                 break;
            }
        }
            break;
            
        case ( kPhypCPVBar ) :
        {
            conv::NuOscDeltaCP_t cpv = Opt::Instance()->gOscDeltaCP;
            switch(cpv) {
                case(kNuDeltaCPPhase) : return "$\\bar{\\delta_{CP}}$";       break;
                case(kNuDeltaCPSin  ) : return "$\\sin(\\bar{\\delta_{CP}})$";break;
                default               : return "--";                   break;
            }
        }
        case ( kPhypDm23sqBar                ) :
        {
            NuOscDmsqScheme_t omp = Opt::Instance()->GetEffectiveDmsqScheme();
            switch(omp) {
                case(kNuDmsqFogliLisi) :        return "$|\\bar{\\Delta m^{2}}| = |m^{2}_{3} - (m^{2}_{1} + m^{2}_{1})/2|$";                  break;
                case(kNuDmsq12andAbs23):        return "$|\\bar{\\Delta m^{2}_{32}}|$";                  break;
                case(kNuDmsq12andAbs31):        return "$|\\bar{\\Delta m^{2}_{31}}|$";                  break;
                default                :        return "--";                                   break;
            }
        }
            break;
            
        case ( kPhypDm12sqBar                ) :
        {
            NuOscDmsqScheme_t omp = Opt::Instance()->GetEffectiveDmsqScheme();
            switch(omp) {
                case(kNuDmsqFogliLisi):
                    return "$\\bar{\\deltam^{2}} = \\bar{\\Delta m^{2}_{21}}$";      break;
                case(kNuDmsq12andAbs23):
                case(kNuDmsq12andAbs31):
                    return "$\\bar{\\Delta m^{2}_{21}}$";                   break;
                default:
                    return "--";                                   break;
            }
        }
            
            
            // NEW BETABAR
        case ( kPhypBetaBar            ) : return "$\\beta{\\beta}$";                break;
            //
            
        case ( kPhypAlphaDecay            ) : return "$\\alpha_{decay}$";                break;
        case ( kPhypSinsq2ThetaDecay      ) : return "$\\sin^{2}(2\\theta_{decay})$";      break;
        case ( kPhypMusqDecoh             ) : return "$\\mu^{2}_{decoh}$";                 break;
        case ( kPhypSinsq2ThetaDecoh      ) : return "$\\sin^{2}(2\\theta_{decoh})$";      break;
            
        case ( kPhypMixing14 ) :
        {
            NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
            switch(omp) {
                case(kNuMixSinsq2Theta) : return "$\\sin^{2}(2\\theta_{14})$";  break;
                case(kNuMixSinsqTheta ) : return "$\\sin^{2}(\\theta_{14})$";    break;
                default                 : return "--";                 break;
            }
        }
            break;
        case ( kPhypMixing15 ) :
        {
            NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
            switch(omp) {
                case(kNuMixSinsq2Theta) : return "$\\sin^{2}(2\\theta_{15})$";  break;
                case(kNuMixSinsqTheta ) : return "$\\sin^{2}(\\theta_{15})$";    break;
                default                 : return "--";                 break;
            }
        }
            break;
        case ( kPhypMixing24 ) :
        {
            NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
            switch(omp) {
                case(kNuMixSinsq2Theta) : return "$\\sin^{2}(2\\theta_{24})$";  break;
                case(kNuMixSinsqTheta ) : return "$\\sin^{2}(\\theta_{24})$";    break;
                default                 : return "--";                 break;
            }
        }
            break;
        case ( kPhypMixing25 ) :
        {
            NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
            switch(omp) {
                case(kNuMixSinsq2Theta) : return "$\\sin^{2}(2\\theta_{25})$";  break;
                case(kNuMixSinsqTheta ) : return "$\\sin^{2}(\\theta_{25})$";    break;
                default                 : return "--";                 break;
            }
        }
            break;
        case ( kPhypMixing34 ) :
        {
            NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
            switch(omp) {
                case(kNuMixSinsq2Theta) : return "$\\sin^{2}(2\\theta_{34})$";  break;
                case(kNuMixSinsqTheta ) : return "$\\sin^{2}(\\theta_{34})$";    break;
                default                 : return "--";                 break;
            }
        }
            break;
        case ( kPhypMixing35 ) :
        {
            NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
            switch(omp) {
                case(kNuMixSinsq2Theta) : return "$\\sin^{2}(2\\theta_{35})$";  break;
                case(kNuMixSinsqTheta ) : return "$\\sin^{2}(\\theta_{35})$";    break;
                default                 : return "--";                 break;
            }
        }
            break;
            
        case ( kPhypDelta12               ) : return "$\\delta_{12}$";                   break;
        case ( kPhypDelta13               ) : return "$\\delta_{13}$";                   break;
        case ( kPhypDelta15               ) : return "$\\delta_{15}$";                   break;
        case ( kPhypDelta34               ) : return "$\\delta_{34}$";                   break;
        case ( kPhypDelta35               ) : return "$\\delta_{35}$";                   break;
        case ( kPhypDm41sq                ) : return "$|\\Delta m^{2}_{41}|$";             break;
        case ( kPhypDm51sq                ) : return "$|\\Delta m^{2}_{51}|$";             break;
            
            // Nonstandard matter interaction
        case ( kPhypNSIepsMuTau           ) : return "$nsi_\\epsilon_{\\mu\\tau}$";      break;
            
            
            /*
             * 2011-analysis nuisance parameters
             */
            
            //
            // FHC systematics
            //
            
            //
            // SK detector
            //
            
        case ( kNuisSuperkEnergyScale     ) : return "$f_{E; r}^{SK}$";                 break;
            
            //
            // BANFF flux and cross-section params
            //
        case ( kNuisBANFF00 ) :
        case ( kNuisBANFF01 ) :
        case ( kNuisBANFF02 ) :
        case ( kNuisBANFF03 ) :
        case ( kNuisBANFF04 ) :
        case ( kNuisBANFF05 ) :
        case ( kNuisBANFF06 ) :
        case ( kNuisBANFF07 ) :
        case ( kNuisBANFF08 ) :
        case ( kNuisBANFF09 ) :
        case ( kNuisBANFF10 ) :
        case ( kNuisBANFF11 ) :
        case ( kNuisBANFF12 ) :
        case ( kNuisBANFF13 ) :
        case ( kNuisBANFF14 ) :
        case ( kNuisBANFF15 ) :
#ifdef __USE_ALLFLUX__
        case ( kNuisBANFF16 ) :
        case ( kNuisBANFF17 ) :
        case ( kNuisBANFF18 ) :
        case ( kNuisBANFF19 ) :
        case ( kNuisBANFF20 ) :
        case ( kNuisBANFF21 ) :
        case ( kNuisBANFF22 ) :
        case ( kNuisBANFF23 ) :
        case ( kNuisBANFF24 ) :
#endif
        {
            int ip = (int)p - (int)kNuisBANFF00;
            return std::string(Form("$f_{%d; t, r}^{banff}$", ip));
        }
            
        case ( kNuisBANFFMAQE ) :            return "$f_{M_A^{QE}}^{banff}$";
        case ( kNuisBANFFMARES ) :            return "$f_{M_A^{RES}}^{banff}$";
        case ( kNuisBANFFCCQEE0 ) :
#if BANFF_BINS != 3
        case ( kNuisBANFFCCQEE1 ) :
        case ( kNuisBANFFCCQEE2 ) :
#endif
	  return std::string(Form("$f_{CCQE E%d}^{banff}$", p-kNuisBANFFCCQEE0));
        case ( kNuisBANFFCC1piE0 ) :
#if BANFF_BINS != 3
        case ( kNuisBANFFCC1piE1 ) :
#endif
	  return std::string(Form("$f_{CC1\\pi E%d}^{banff}$", p-kNuisBANFFCC1piE0));
	//return std::string(Form("$f_{CC1\\pi E%d}^{banff}$", p-kNuisBANFFCC1piE1)); // NEW SPLIT ANA
            
#if (BANFF_BINS != 2)
        case ( kNuisBANFFNC1pi0 ) :      return "$f_{NC1\\pi0}^{banff}$";
#endif
            // NEW XSEC
        case ( kNuisBANFFNorm2p2h   ) :  return "$f_{Norm_{2p2h}}^{banff}$";
        case ( kNuisBANFFNorm2p2hBar) :  return "$f_{Norm_{2p2hBar}}^{banff}$";
        case ( kNuisBANFFCA5   ) :       return "$f_{C^{A}_{5}}^{banff}$";
        case ( kNuisBANFFBgRES ):        return "$f_{Bg_{RES}}^{banff}$";
            // NEW SPLIT ANA
#if (BANFF_BINS != 2)
        case ( kNuisBANFFSCCA  ) :       return "$f_{SCCA}^{banff}$";
        case ( kNuisBANFFSCCV  ) :       return "$f_{SCCV}^{banff}$";
#endif
        case ( kNuisBANFFFermiMom   ) : return "$f_{p_{f}}^{banff}$";
        case ( kNuisBANFFShapeCCoth ) : return "$f_{Shape_{CCoth}}^{banff}$";
          //case ( kNuisBANFFBindE      ) : return "$f_{E_{B}}^{banff}$";
        case ( kNuisBANFFNormCCcoh  ) : return "$f_{Norm_{CCcoh}}^{banff}$";
        case ( kNuisBANFFNormNCcoh  ) : return "$f_{Norm_{NCcoh}}^{banff}$";
        case ( kNuisBANFFNormNCoth  ) : return "$f_{Norm_{NCother}}^{banff}$";
        case ( kNuisBANFFNormNuEToNuMu  ) : return "$f_{Norm_{\\nu_e \\rightarrow \\nu_\\mu}}^{banff}$";
        case ( kNuisBANFFNormNC1Gamma   ) : return "$f_{Norm_{NC1\\gamma}}^{banff}$";
        case ( kNuisBANFFNormNuEBarToNuMuBar  ) : return "$f_{Norm_{\\bar{\\nu}_e \\rightarrow \\bar{\\nu}_{\\mu}}}^{banff}$";
        case ( kNuisBANFFBeRPA_A      ) : return "$f_{Shape_{BeRPA_{A}}}^{banff}$";       break;
        case ( kNuisBANFFBeRPA_B      ) : return "$f_{Shape_{BeRPA_{B}}}^{banff}$";       break;
        case ( kNuisBANFFBeRPA_D      ) : return "$f_{Shape_{BeRPA_{D}}}^{banff}$";       break;
        case ( kNuisBANFFBeRPA_E      ) : return "$f_{Shape_{BeRPA_{E}}}^{banff}$";       break;
        case ( kNuisBANFFBeRPA_U      ) : return "$f_{Shape_{BeRPA_{U}}}^{banff}$";       break;
        case ( kNuisBANFFShape2p2h    ) : return "$f_{Shape_{2p2h\\nu}}^{banff}$";        break;
        case ( kNuisBANFFNorm2p2hCtoO ) : return "$f_{Norm_{2p2h C to O}}^{banff}$"; break;
        case ( kNuisBANFFEb    ) : return "$f_{Shape_{Eb}}^{banff}$";
        case ( kNuisBANFFMultiReNorm  ) : return "$f_{Norm_{MultiRe}}^{banff}$";
            //
            
            //
            // XSec params not correlated with flux
            //
            
            // NEW XSEC
        case ( kNuisXSecRPA        ) : return "$f_{RPA; t, r}$";                 break;
            
        case ( kNuisXSecFermiMom   ) : return "$f_{pF; t, r}$";                  break;
        case ( kNuisXSecSpecFunc   ) : return "$f_{SF; t, r}$";                  break;
        case ( kNuisXSecWshape     ) : return "$f_{Wshape; t, r}$";              break;
        case ( kNuisXSecShapeCCoth ) : return "$f_{CCothShape; t, r}$";          break;
        case ( kNuisXSecPiLessDelta) : return "$f_{\\pi-less\\Delta; t, r}$";    break;
#ifndef __NUEONLY__
        case ( kNuisXSecBindE      ) : return "$f_{bindE; t, r}$";               break;
#endif
        case ( kNuisXSecNormCCcoh  ) : return "$f_{CCcoh; t}$";                  break;
#ifndef __NUMUONLY__
        case ( kNuisXSecNormNCcoh  ) : return "$f_{NCcoh; t}$";                  break;
#endif
#ifndef __NUEONLY__
        case ( kNuisXSecNormNC1piC ) : return "$f_{NC1\\pi^\\pm; t}$";           break;
#endif
        case ( kNuisXSecNormNCoth  ) : return "$f_{NCoth; t}$";                  break;
        case ( kNuisXSecNormNuEToNuMu  ) : return "$f_{CC\\nu_{e}; t}$";         break;
#ifndef __NUEONLY__
        case ( kNuisXSecNormAntiNuToNu ) : return "$f_{CC\\bar{\\nu}; t}$";      break;
#endif
            
            
            
#if FSI_BINS != 3
        case ( kNuisSKDetFSI000 ) :
        case ( kNuisSKDetFSI001 ) :
        case ( kNuisSKDetFSI002 ) :
        case ( kNuisSKDetFSI003 ) :
        case ( kNuisSKDetFSI004 ) :
        case ( kNuisSKDetFSI005 ) :
            
            // NEW SKDETFSI
#if defined(__USE_SKDETFSI2015__)
        case ( kNuisSKDetFSI006 ) :
        case ( kNuisSKDetFSI007 ) :
        case ( kNuisSKDetFSI008 ) :
        case ( kNuisSKDetFSI009 ) :
        case ( kNuisSKDetFSI010 ) :
        case ( kNuisSKDetFSI011 ) :
#endif
            
#endif
#if FSI_BINS != 2
            
            // NEW SKDETFSI
#if defined(__USE_SKDETFSI2013__)
        case ( kNuisSKDetFSI006 ) :
        case ( kNuisSKDetFSI007 ) :
        case ( kNuisSKDetFSI008 ) :
        case ( kNuisSKDetFSI009 ) :
        case ( kNuisSKDetFSI010 ) :
        case ( kNuisSKDetFSI011 ) :
        case ( kNuisSKDetFSI012 ) :
        case ( kNuisSKDetFSI013 ) :
        case ( kNuisSKDetFSI014 ) :
        case ( kNuisSKDetFSI015 ) :
        case ( kNuisSKDetFSI016 ) :
        case ( kNuisSKDetFSI017 ) :
#elif defined(__USE_SKDETFSI2013__)
        case ( kNuisSKDetFSI012 ) :
        case ( kNuisSKDetFSI013 ) :
        case ( kNuisSKDetFSI014 ) :
        case ( kNuisSKDetFSI015 ) :
        case ( kNuisSKDetFSI016 ) :
        case ( kNuisSKDetFSI017 ) :
        case ( kNuisSKDetFSI018 ) :
        case ( kNuisSKDetFSI019 ) :
        case ( kNuisSKDetFSI020 ) :
        case ( kNuisSKDetFSI021 ) :
        case ( kNuisSKDetFSI022 ) :
        case ( kNuisSKDetFSI023 ) :
#endif
            
#endif
            // NEW SKDETFSI (start indices from 24 instead of 18)
#ifdef SKCORRELATIONS
        case ( kNuisSKDetFSI024 ) :
        case ( kNuisSKDetFSI025 ) :
        case ( kNuisSKDetFSI026 ) :
        case ( kNuisSKDetFSI027 ) :
        case ( kNuisSKDetFSI028 ) :
        case ( kNuisSKDetFSI029 ) :
        case ( kNuisSKDetFSI030 ) :
        case ( kNuisSKDetFSI031 ) :
        case ( kNuisSKDetFSI032 ) :
        case ( kNuisSKDetFSI033 ) :
        case ( kNuisSKDetFSI034 ) :
        case ( kNuisSKDetFSI035 ) :
        case ( kNuisSKDetFSI036 ) :
        case ( kNuisSKDetFSI037 ) :
#endif
            
        {
            // #if FSI_BINS != 3
            //       int ip = (int)p - (int)kNuisSKDetFSI000;
            // #else
            //       int ip = (int)p - (int)kNuisSKDetFSI006;
            // #endif
            
            // NEW SKDETFSI
            int ip = (int)p - (int)kNuisSKDetFSI_first;
            
            return std::string(Form("$f_{%d; t, r}^{SK+FSI}$", ip));
            break;
        }//FSI
            
        case ( kNuisRHCToFHC ) : return "$f^{RHC/FHC}$"; break;
            
            //
            // RHC systematics
            //
            
            //
            // SK detector
            //
        case ( kRHCNuisSuperkEnergyScale     ) : return "$f_{E; r}^{SK}RHC$";                 break;
            
            //
            // BANFF flux and cross-section params
            //
        case ( kRHCNuisBANFF00 ) :
        case ( kRHCNuisBANFF01 ) :
        case ( kRHCNuisBANFF02 ) :
        case ( kRHCNuisBANFF03 ) :
        case ( kRHCNuisBANFF04 ) :
        case ( kRHCNuisBANFF05 ) :
        case ( kRHCNuisBANFF06 ) :
        case ( kRHCNuisBANFF07 ) :
        case ( kRHCNuisBANFF08 ) :
        case ( kRHCNuisBANFF09 ) :
        case ( kRHCNuisBANFF10 ) :
        case ( kRHCNuisBANFF11 ) :
        case ( kRHCNuisBANFF12 ) :
        case ( kRHCNuisBANFF13 ) :
        case ( kRHCNuisBANFF14 ) :
        case ( kRHCNuisBANFF15 ) :
#ifdef __USE_ALLFLUX__
        case ( kRHCNuisBANFF16 ) :
        case ( kRHCNuisBANFF17 ) :
        case ( kRHCNuisBANFF18 ) :
        case ( kRHCNuisBANFF19 ) :
        case ( kRHCNuisBANFF20 ) :
        case ( kRHCNuisBANFF21 ) :
        case ( kRHCNuisBANFF22 ) :
        case ( kRHCNuisBANFF23 ) :
        case ( kRHCNuisBANFF24 ) :
#endif
        {
            int ip = (int)p - (int)kRHCNuisBANFF00;
            return std::string(Form("$f_{%d; t, r}^{banff} RHC$", ip));
        }
            
        case ( kRHCNuisBANFFMAQE ) :
            return "$f_{M_A^{QE}}^{banff} RHC$";
        case ( kRHCNuisBANFFMARES ) :
            return "$f_{M_A^{RES}}^{banff} RHC$";
        case ( kRHCNuisBANFFCCQEE0 ) :
#if BANFF_BINS != 3
        case ( kRHCNuisBANFFCCQEE1 ) :
        case ( kRHCNuisBANFFCCQEE2 ) :
#endif
	  return std::string(Form("$f_{CCQE E%d}^{banff} RHC$", p-kRHCNuisBANFFCCQEE0));
        case ( kRHCNuisBANFFCC1piE0 ) :
#if BANFF_BINS != 3
        case ( kRHCNuisBANFFCC1piE1 ) :
#endif
	  return std::string(Form("$f_{CC1\\pi E%d}^{banff RHC}$", p-kRHCNuisBANFFCC1piE0));
            //return std::string(Form("$f_{CC1\\pi E%d}^{banff RHC}$", p-kRHCNuisBANFFCC1piE1)); // NEW SPLIT ANA
#if (BANFF_BINS != 2)
        case ( kRHCNuisBANFFNC1pi0 ) :
#endif
            return "$f_{NC1\\pi0}^{banff} RHC$";
            
            // NEW XSEC
        case ( kRHCNuisBANFF2p2h   ) :
            return "$f_{Norm_{2p2h}}^{banff} RHC$";
        case ( kRHCNuisBANFF2p2hBar) :
            return "$f_{Norm_{2p2hBar}}^{banff} RHC$";
        case ( kRHCNuisBANFFCA5   ) :
            return "$f_{C^{A}_{5}}^{banff} RHC$";
        case ( kRHCNuisBANFFBgRES ):
            return "$f_{Bg_{RES}}^{banff} RHC$";
            // NEW SPLIT ANA
#if (BANFF_BINS != 2)
        case ( kRHCNuisBANFFSCCA  ) :
            return "$f_{SCCA}^{banff} RHC$";
        case ( kRHCNuisBANFFSCCV  ) :
            return "$f_{SCCV}^{banff} RHC$";
#endif
        case ( kRHCNuisBANFFFermiMom   ) :     // uncorrelated in the past
            return "$f_{p_{f}}^{banff} RHC$";
        case ( kRHCNuisBANFFShapeCCoth ) :     // uncorrelated in the past
            return "$f_{Shape_{CCoth}}^{banff} RHC$";
            //case ( kRHCNuisBANFFBindE      ) :     // uncorrelated in the past
            //return "$f_{E_{B}}^{banff} RHC$";
        case ( kRHCNuisBANFFNormCCcoh  ) :     // uncorrelated in the past
            return "$f_{Norm_{CCcoh}}^{banff} RHC$";
        case ( kRHCNuisBANFFNormNCcoh  ) :     // uncorrelated in the past
            return "$f_{Norm_{NCcoh}}^{banff} RHC$";
        case ( kRHCNuisBANFFNormNCoth  ) :     // uncorrelated in the past
            return "$f_{Norm_{NCother}}^{banff} RHC$";
        case ( kRHCNuisBANFFNormNuEToNuMu  ) : // uncorrelated in the past
            return "$f_{Norm_{\\nu_e \\rightarrow \\nu_\\mu}}^{banff} RHC$";
        case ( kRHCNuisBANFFNormNC1Gamma  ) :     // uncorrelated in the past
            return "$f_{Norm_{NC1\\gamma}}^{banff} RHC$";
        case ( kRHCNuisBANFFNormNuEBarToNuMuBar  ) : // uncorrelated in the past
            return "$f_{Norm_{\\bar{\\nu}_e \\rightarrow \\bar{\\nu}_{\\mu}}}^{banff} RHC$";
            
            //
            // XSec params not correlated with flux
            //
            
            // NEW XSEC
        case ( kRHCNuisXSecRPA        ) : return "$f_{RPA; t, r}RHC$";                 break;
            
        case ( kRHCNuisXSecFermiMom   ) : return "$f_{pF; t, r}RHC$";                  break;
        case ( kRHCNuisXSecSpecFunc   ) : return "$f_{SF; t, r}RHC$";                  break;
        case ( kRHCNuisXSecWshape     ) : return "$f_{Wshape; t, r}RHC$";              break;
        case ( kRHCNuisXSecShapeCCoth ) : return "$f_{CCothShape; t, r}RHC$";          break;
        case ( kRHCNuisXSecPiLessDelta) : return "$f_{\\pi-less\\Delta; t, r}RHC$";    break;
#ifndef __NUEONLY__
        case ( kRHCNuisXSecBindE      ) : return "$f_{bindE; t, r}RHC$";               break;
#endif
        case ( kRHCNuisXSecNormCCcoh  ) : return "$f_{CCcoh; t}RHC$";                  break;
#ifndef __NUMUONLY__
        case ( kRHCNuisXSecNormNCcoh  ) : return "$f_{NCcoh; t}RHC$";                  break;
#endif
#ifndef __NUEONLY__
        case ( kRHCNuisXSecNormNC1piC ) : return "$f_{NC1\\pi^\\pm; t}RHC$";           break;
#endif
        case ( kRHCNuisXSecNormNCoth  ) : return "$f_{NCoth; t}RHC$";                  break;
        case ( kRHCNuisXSecNormNuEToNuMu  ) : return "$f_{CC\\nu_{e}; t}RHC$";           break;
#ifndef __NUEONLY__
        case ( kRHCNuisXSecNormAntiNuToNu ) : return "$f_{CC\\bar{\\nu}; t}RHC$";      break;
#endif
            
#if FSI_BINS != 3
        case ( kRHCNuisSKDetFSI000 ) :
        case ( kRHCNuisSKDetFSI001 ) :
        case ( kRHCNuisSKDetFSI002 ) :
        case ( kRHCNuisSKDetFSI003 ) :
        case ( kRHCNuisSKDetFSI004 ) :
        case ( kRHCNuisSKDetFSI005 ) :
            
            // NEW SKDETFSI
#if defined(__USE_SKDETFSI2015__)
        case ( kRHCNuisSKDetFSI006 ) :
        case ( kRHCNuisSKDetFSI007 ) :
        case ( kRHCNuisSKDetFSI008 ) :
        case ( kRHCNuisSKDetFSI009 ) :
        case ( kRHCNuisSKDetFSI010 ) :
        case ( kRHCNuisSKDetFSI011 ) :
#endif
            
#endif
#if FSI_BINS != 2
            
            // NEW SKDETFSI
#if defined(__USE_SKDETFSI2013__)
        case ( kRHCNuisSKDetFSI006 ) :
        case ( kRHCNuisSKDetFSI007 ) :
        case ( kRHCNuisSKDetFSI008 ) :
        case ( kRHCNuisSKDetFSI009 ) :
        case ( kRHCNuisSKDetFSI010 ) :
        case ( kRHCNuisSKDetFSI011 ) :
        case ( kRHCNuisSKDetFSI012 ) :
        case ( kRHCNuisSKDetFSI013 ) :
        case ( kRHCNuisSKDetFSI014 ) :
        case ( kRHCNuisSKDetFSI015 ) :
        case ( kRHCNuisSKDetFSI016 ) :
        case ( kRHCNuisSKDetFSI017 ) :
#elif defined(__USE_SKDETFSI2013__)
        case ( kRHCNuisSKDetFSI012 ) :
        case ( kRHCNuisSKDetFSI013 ) :
        case ( kRHCNuisSKDetFSI014 ) :
        case ( kRHCNuisSKDetFSI015 ) :
        case ( kRHCNuisSKDetFSI016 ) :
        case ( kRHCNuisSKDetFSI017 ) :
        case ( kRHCNuisSKDetFSI018 ) :
        case ( kRHCNuisSKDetFSI019 ) :
        case ( kRHCNuisSKDetFSI020 ) :
        case ( kRHCNuisSKDetFSI021 ) :
        case ( kRHCNuisSKDetFSI022 ) :
        case ( kRHCNuisSKDetFSI023 ) :
#endif
            
#endif
            
            // NEW SKDETFSI (start indices from 24 instead of 18)
#ifdef SKCORRELATIONS
        case ( kRHCNuisSKDetFSI024 ) :
        case ( kRHCNuisSKDetFSI025 ) :
        case ( kRHCNuisSKDetFSI026 ) :
        case ( kRHCNuisSKDetFSI027 ) :
        case ( kRHCNuisSKDetFSI028 ) :
        case ( kRHCNuisSKDetFSI029 ) :
        case ( kRHCNuisSKDetFSI030 ) :
        case ( kRHCNuisSKDetFSI031 ) :
        case ( kRHCNuisSKDetFSI032 ) :
        case ( kRHCNuisSKDetFSI033 ) :
        case ( kRHCNuisSKDetFSI034 ) :
        case ( kRHCNuisSKDetFSI035 ) :
        case ( kRHCNuisSKDetFSI036 ) :
        case ( kRHCNuisSKDetFSI037 ) :
#endif
            
        {
            // #if FSI_BINS != 3
            //       int ip = (int)p - (int)kRHCNuisSKDetFSI000;
            // #else
            //       int ip = (int)p - (int)kRHCNuisSKDetFSI006;
            // #endif
            
            // NEW SKDETFSI
            int ip = (int)p - (int)kRHCNuisSKDetFSI_first;
            
            return std::string(Form("$f_{%d; t, r}^{SK+FSI}RHC$", ip));
            break;
        }//FSI
            
            //
            // MultiRing FHC systematics
            //
            
            /// NEW MULTIRING 3
        case ( kMultiRingNuisSuperkEnergyScale     ) : return "$f_{E; r}^{SK}MultiRing$";     break;
            ///
            
            /// NEW MULTIRING 2
            
#if MR_FSI_BINS != 3
        case ( kMultiRingNuisSKDetFSI000 ) :
        case ( kMultiRingNuisSKDetFSI001 ) :
        case ( kMultiRingNuisSKDetFSI002 ) :
        case ( kMultiRingNuisSKDetFSI003 ) :
        case ( kMultiRingNuisSKDetFSI004 ) :
        case ( kMultiRingNuisSKDetFSI005 ) :
#endif
#if MR_FSI_BINS != 2
            
            // NEW SKDETFSI
#if defined(__USE_SKDETFSI2013__)
        case ( kMultiRingNuisSKDetFSI006 ) :
        case ( kMultiRingNuisSKDetFSI007 ) :
        case ( kMultiRingNuisSKDetFSI008 ) :
        case ( kMultiRingNuisSKDetFSI009 ) :
        case ( kMultiRingNuisSKDetFSI010 ) :
        case ( kMultiRingNuisSKDetFSI011 ) :
        case ( kMultiRingNuisSKDetFSI012 ) :
        case ( kMultiRingNuisSKDetFSI013 ) :
        case ( kMultiRingNuisSKDetFSI014 ) :
        case ( kMultiRingNuisSKDetFSI015 ) :
        case ( kMultiRingNuisSKDetFSI016 ) :
        case ( kMultiRingNuisSKDetFSI017 ) :
#endif
#endif            
            // NEW SKDETFSI (start indices from 24 instead of 18)
#ifdef SKCORRELATIONS
        case ( kMultiRingNuisSKDetFSI024 ) :
        case ( kMultiRingNuisSKDetFSI025 ) :
        case ( kMultiRingNuisSKDetFSI026 ) :
        case ( kMultiRingNuisSKDetFSI027 ) :
        case ( kMultiRingNuisSKDetFSI028 ) :
        case ( kMultiRingNuisSKDetFSI029 ) :
        case ( kMultiRingNuisSKDetFSI030 ) :
        case ( kMultiRingNuisSKDetFSI031 ) :
        case ( kMultiRingNuisSKDetFSI032 ) :
        case ( kMultiRingNuisSKDetFSI033 ) :
        case ( kMultiRingNuisSKDetFSI034 ) :
        case ( kMultiRingNuisSKDetFSI035 ) :
        case ( kMultiRingNuisSKDetFSI036 ) :
        case ( kMultiRingNuisSKDetFSI037 ) :
#endif
            
        {
            // #if FSI_BINS != 3
            //       int ip = (int)p - (int)kMultiRingNuisSKDetFSI000;
            // #else
            //       int ip = (int)p - (int)kMultiRingNuisSKDetFSI006;
            // #endif
            
            // NEW SKDETFSI
	    int ip = ArrIdx::Instance()->MultiRingSKDetFSISystIdx(p);
            
            return std::string(Form("$f_{%d; t, r}^{SK+FSI}MultiRing$", ip));
            break;
        }//FSI
            
            ///
        default :
            return "--";
    }
    return "--";
}
//________________________________________________________________________________________
string valor::util::AsStringPaperTable(PhysParamId_t p)
{
  if (util::IsBANFF(p) && util::IsFlux(p)) {
    if (util::IsRHCSystematic(p)) {
      int ip = (int)p - (int)kRHCNuisBANFF00;
      return std::string(Form("RHC %s flux normalisation, E = %.1f - %.1f GeV",
                              util::Pdg2StringTex(kNuPdgCodeOfRHCFluxParamInBANFFfit[ip]).c_str(),
                              kEtrueRangeOfRHCFluxParamInBANFFfit[ip][0],
                              kEtrueRangeOfRHCFluxParamInBANFFfit[ip][1]));
    }
    else {
      int ip = (int)p - (int)kNuisBANFF00;
      return std::string(Form("FHC %s flux normalisation, E = %.1f - %.1f GeV",
                              util::Pdg2StringTex(kNuPdgCodeOfFluxParamInBANFFfit[ip]).c_str(),
                              kEtrueRangeOfFluxParamInBANFFfit[ip][0],
                              kEtrueRangeOfFluxParamInBANFFfit[ip][1]));
    }
  }
  else if (util::IsBANFF(p)) {
    switch (p) {
    case ( kNuisBANFFMAQE              ) :
    case ( kRHCNuisBANFFMAQE           ) :
      return "CCQE axial-mass scaling factor";
    case ( kNuisBANFFMARES             ) :
    case ( kRHCNuisBANFFMARES          ) :
      return "Resonance-production axial-mass scaling factor";
    case ( kNuisBANFFNorm2p2h           ) :
    case ( kRHCNuisBANFF2p2h            ) :
      return "Two particle two hole normalisation for $^{16}$O";
    case ( kNuisBANFFNorm2p2hBar        ) :
    case ( kRHCNuisBANFF2p2hBar         ) :
      return "Antineutrino two particle two hole normalisation for $^{16}$O";
    case ( kNuisBANFFCA5               ) :
    case ( kRHCNuisBANFFCA5            ) :
      return "C$^{A}_{5}$ nucleon to $\\Delta$ transition axial form factor";
    case ( kNuisBANFFBgRES             ) :
    case ( kRHCNuisBANFFBgRES          ) :
      return "Scale of isospin 1/2 nonresonant background";
    case ( kNuisBANFFFermiMom          ) :
    case ( kRHCNuisBANFFFermiMom       ) :
      return "Fermi momentum for $^{16}$O";
    case ( kNuisBANFFShapeCCoth        ) :
    case ( kRHCNuisBANFFShapeCCoth     ) :
      return "CC other shape";
    //case ( kNuisBANFFBindE             ) :
    //case ( kRHCNuisBANFFBindE          ) :
    //return "Binding energy for $^{16}$O";
    case ( kNuisBANFFNormCCcoh         ) :
    case ( kRHCNuisBANFFNormCCcoh      ) :
      return "CC coherent for $^{16}$O normalisation";
    case ( kNuisBANFFNormNCcoh         ) :
    case ( kRHCNuisBANFFNormNCcoh      ) :
      return "NC coherent normalisation";
    case ( kNuisBANFFNormNCoth         ) :
    case ( kRHCNuisBANFFNormNCoth      ) :
      return "NC other normalisation";
    case ( kNuisBANFFNormNuEToNuMu     ) :
    case ( kRHCNuisBANFFNormNuEToNuMu  ) :
      return "CC $\\nu_{e}$ normalisation";
    case ( kNuisBANFFNormNC1Gamma     ) :
    case ( kRHCNuisBANFFNormNC1Gamma  ) :
      return "NC 1$\\gamma$ normalisation";
    case ( kNuisBANFFNormNuEBarToNuMuBar     ) :
    case ( kRHCNuisBANFFNormNuEBarToNuMuBar  ) :
      return "CC $\\bar{\\nu}_{e}$ normalisation";

    case ( kNuisBANFFBeRPA_A         ) : 
    case ( kRHCNuisBANFFBeRPA_A      ) : 
      return "Bernstein Polynomial coefficient A";    break;
    case ( kNuisBANFFBeRPA_B         ) : 
    case ( kRHCNuisBANFFBeRPA_B      ) : 
      return "Bernstein Polynomial coefficient B";    break;
    case ( kNuisBANFFBeRPA_D         ) : 
    case ( kRHCNuisBANFFBeRPA_D      ) : 
      return "Bernstein Polynomial coefficient D";    break;
    case ( kNuisBANFFBeRPA_E         ) : 
    case ( kRHCNuisBANFFBeRPA_E      ) : 
      return "Bernstein Polynomial coefficient E";    break;
    case ( kNuisBANFFBeRPA_U         ) : 
    case ( kRHCNuisBANFFBeRPA_U      ) : 
      return "Bernstein Polynomial coefficient U";    break;
    case ( kNuisBANFFShape2p2h       ) : 
    case ( kRHCNuisBANFFShape2p2h    ) : 
      return "Neutrino two particle two hole $^{16}$O shape";     break;
    case ( kNuisBANFFNorm2p2hCtoO    ) : 
    case ( kRHCNuisBANFFNorm2p2hCtoO ) : 
      return "Two particle two hole $^{12}$C to $^{16}$O normalisation";  break;
    case ( kNuisBANFFEb       ) : 
    case ( kRHCNuisBANFFEb    ) : 
      return "Binding energy shape";     break;
    case ( kNuisBANFFMultiReNorm ) :
      return "FHC MultiRe Normalisation"; break;
    case ( kRHCNuisBANFFMultiReNorm ) :
      return "RHC MultiRe Normalisation. UNUSED"; break;
                
    default:
      return util::AsString(p);
    }
  }
  //~ else if(util::IsSKDetFSI(p) {
  //~ if (util::IsRHCSystematic(p)) {
  //~ int ip = (int)p - (int)kRHCNuisSKDetFSI_first;
  //~ }
  //~ return std::string(Form("SKDet + FSI/SI %2d; Ereco range %.2f - %.2f GeV; %s",ip,
  //~ kSKDetFSIParamsErecoCut[ip][0],kSKDetFSIParamsErecoCut[ip][1],nutypes.c_str());
  //~ break;
  //~ }
  //~ return std::string(Form("RHC"));
  //~ }
  //~ else {
  //~ int ip = (int)p - (int)kNuisSKDetFSI_first;
  //~ return std::string(Form("FHC"));
  //~ }
  //~ }
  return util::AsString(p);
}
//________________________________________________________________________________________
//________________________________________________________________________________________
string valor::util::AsStringPaperRooTex(PhysParamId_t p)
{
  string s = AsStringPaperTexNoDollar(p);
  replace( s.begin(), s.end(), '\\', '#' ); // Make slashes hashes for RooTex
  // Get rid of hashes preceding sin
  size_t ind = s.find("#sin",0);
  if (ind !=  string::npos) (s.replace(ind,4," sin"));
  return util::TrimSpaces(s);
}
//________________________________________________________________________________________
string valor::util::AsStringPaperTexNoDollar(PhysParamId_t p)
{
  string s = AsStringPaperTex(p);
  s.erase( std::remove( s.begin(), s.end(), '$'), s.end() ) ; // Remove dollar signs.
  return s;
}
//____________________________________________________________________________

string valor::util::AsStringCompact_splinehack(PhysParamId_t p)
{
    switch(p) {
            
            //***************************** REVERSE HORN CURRENT *************************
            
        case ( kRHCNuisSuperkEnergyScale     ) : return "f_sk_e_scale";                   break;
            
        case ( kRHCNuisBANFF00               ) : return "f_banff00";      break;
        case ( kRHCNuisBANFF01               ) : return "f_banff01";      break;
        case ( kRHCNuisBANFF02               ) : return "f_banff02";      break;
        case ( kRHCNuisBANFF03               ) : return "f_banff03";      break;
        case ( kRHCNuisBANFF04               ) : return "f_banff04";      break;
        case ( kRHCNuisBANFF05               ) : return "f_banff05";      break;
        case ( kRHCNuisBANFF06               ) : return "f_banff06";      break;
        case ( kRHCNuisBANFF07               ) : return "f_banff07";      break;
        case ( kRHCNuisBANFF08               ) : return "f_banff08";      break;
        case ( kRHCNuisBANFF09               ) : return "f_banff09";      break;
        case ( kRHCNuisBANFF10               ) : return "f_banff10";      break;
        case ( kRHCNuisBANFF11               ) : return "f_banff11";      break;
        case ( kRHCNuisBANFF12               ) : return "f_banff12";      break;
        case ( kRHCNuisBANFF13               ) : return "f_banff13";      break;
        case ( kRHCNuisBANFF14               ) : return "f_banff14";      break;
        case ( kRHCNuisBANFF15               ) : return "f_banff15";      break;
#ifdef __USE_ALLFLUX__
        case ( kRHCNuisBANFF16               ) : return "f_banff16";      break;
        case ( kRHCNuisBANFF17               ) : return "f_banff17";      break;
        case ( kRHCNuisBANFF18               ) : return "f_banff18";      break;
        case ( kRHCNuisBANFF19               ) : return "f_banff19";      break;
        case ( kRHCNuisBANFF20               ) : return "f_banff20";      break;
        case ( kRHCNuisBANFF21               ) : return "f_banff21";      break;
        case ( kRHCNuisBANFF22               ) : return "f_banff22";      break;
        case ( kRHCNuisBANFF23               ) : return "f_banff23";      break;
        case ( kRHCNuisBANFF24               ) : return "f_banff24";      break;
#endif
            
#ifndef ONEOVERMAQE
        case ( kRHCNuisBANFFMAQE             ) : return "f_banffmaqe";    break;
#else
        case ( kRHCNuisBANFFMAQE             ) : return "f_banff1overmaqe2";  break;
#endif
        case ( kRHCNuisBANFFMARES            ) : return "f_banffmares";   break;
        case ( kRHCNuisBANFFCCQEE0           ) : return "f_banffccqee0";  break;
#if BANFF_BINS != 3
        case ( kRHCNuisBANFFCCQEE1           ) : return "f_banffccqee1";  break;
        case ( kRHCNuisBANFFCCQEE2           ) : return "f_banffccqee2";  break;
#endif
        case ( kRHCNuisBANFFCC1piE0          ) : return "f_banffcc1pie0"; break;
#if BANFF_BINS != 3
        case ( kRHCNuisBANFFCC1piE1          ) : return "f_banffcc1pie1"; break;
#endif
#if (BANFF_BINS != 2)
        case ( kRHCNuisBANFFNC1pi0           ) : return "f_banffnc1pi0";  break;
#endif
        case ( kRHCNuisXSecFermiMom          ) : return "f_pfermi";       break;
        case ( kRHCNuisXSecSpecFunc          ) : return "f_specfunc";     break;
        case ( kRHCNuisXSecWshape            ) : return "f_wshape";       break;
        case ( kRHCNuisXSecShapeCCoth        ) : return "f_shapeCCoth";   break;
        case ( kRHCNuisXSecPiLessDelta       ) : return "f_pilessdelta";  break;
#ifndef __NUEONLY__
        case ( kRHCNuisXSecBindE             ) : return "f_bindE";        break;
#endif
        case ( kRHCNuisXSecNormCCcoh         ) : return "f_norm_cccoh";   break;
#ifndef __NUEONLY__
        case ( kRHCNuisXSecNormNC1piC        ) : return "f_norm_nc1piC";  break;
#endif
#ifndef __NUMUONLY__
        case ( kRHCNuisXSecNormNCcoh         ) : return "f_norm_nccoh";   break;
#endif
        case ( kRHCNuisXSecNormNCoth         ) : return "f_norm_ncoth";   break;
        case ( kRHCNuisXSecNormNuEToNuMu     ) : return "f_norm_nuetonumu";  break;
#ifndef __NUEONLY__
        case ( kRHCNuisXSecNormAntiNuToNu    ) : return "f_norm_antinutonu"; break;
#endif
            
            // NEW XSEC
        case ( kRHCNuisXSecRPA            ) : return "f_rpa_rhc";                 break;
            
            // NEW XSEC
        case ( kRHCNuisBANFF2p2h           ) : return "f_banff2p2h_rhc";            break;
        case ( kRHCNuisBANFFCA5           ) : return "f_banffca5_rhc";            break;
        case ( kRHCNuisBANFFBgRES         ) : return "f_banffbgres_rhc";          break;
            // NEW SPLIT ANA
#if (BANFF_BINS != 2)
        case ( kRHCNuisBANFFSCCA          ) : return "f_banffscca_rhc";           break;
        case ( kRHCNuisBANFFSCCV          ) : return "f_banffsccv_rhc";           break;
#endif
        case ( kRHCNuisBANFFFermiMom      ) : return "f_banffpfermi_rhc";         break;    // uncorrelated in the past
        case ( kRHCNuisBANFFShapeCCoth    ) : return "f_banffshapeCCoth_rhc";     break;     // uncorrelated in the past
          //case ( kRHCNuisBANFFBindE         ) : return "f_banffbindE_rhc";          break;    // uncorrelated in the past
        case ( kRHCNuisBANFFNormCCcoh     ) : return "f_banffnorm_cccoh_rhc";     break;    // uncorrelated in the past
        case ( kRHCNuisBANFFNormNCcoh     ) : return "f_banffnorm_nccoh_rhc";     break;    // uncorrelated in the past
        case ( kRHCNuisBANFFNormNCoth     ) : return "f_banffnorm_ncoth_rhc";     break;    // uncorrelated in the past
        case ( kRHCNuisBANFFNormNuEToNuMu ) : return "f_banffnorm_nuetonumu_rhc"; break; // uncorrelated in the past
        case ( kRHCNuisBANFFNormNC1Gamma      ) : return "f_banffnorm_nc1gamma_rhc";  break; // uncorrelated in the past
        case ( kRHCNuisBANFFNormNuEBarToNuMuBar ) : return "f_banffnorm_nuebartonumubar_rhc"; break; // uncorrelated in the past
        case ( kRHCNuisBANFF2p2hBar        ) : return "f_banff2p2hbar_rhc";            break;
            
            //
            
            //#ifdef SKFSITOGETHER
#if FSI_BINS != 3
        case ( kRHCNuisSKDetFSI000 ) : return "f_skdetfsi000"; break;
        case ( kRHCNuisSKDetFSI001 ) : return "f_skdetfsi001"; break;
        case ( kRHCNuisSKDetFSI002 ) : return "f_skdetfsi002"; break;
        case ( kRHCNuisSKDetFSI003 ) : return "f_skdetfsi003"; break;
        case ( kRHCNuisSKDetFSI004 ) : return "f_skdetfsi004"; break;
        case ( kRHCNuisSKDetFSI005 ) : return "f_skdetfsi005"; break;
            
            // NEW SKDETFSI
#if defined(__USE_SKDETFSI2015__)
        case ( kRHCNuisSKDetFSI006 ) : return "f_skdetfsi006"; break;
        case ( kRHCNuisSKDetFSI007 ) : return "f_skdetfsi007"; break;
        case ( kRHCNuisSKDetFSI008 ) : return "f_skdetfsi008"; break;
        case ( kRHCNuisSKDetFSI009 ) : return "f_skdetfsi009"; break;
        case ( kRHCNuisSKDetFSI010 ) : return "f_skdetfsi010"; break;
        case ( kRHCNuisSKDetFSI011 ) : return "f_skdetfsi011"; break;
#endif
            
#endif
#if FSI_BINS != 2
            
            // NEW SKDETFSI
#if defined(__USE_SKDETFSI2013__)
        case ( kRHCNuisSKDetFSI006 ) : return "f_skdetfsi006"; break;
        case ( kRHCNuisSKDetFSI007 ) : return "f_skdetfsi007"; break;
        case ( kRHCNuisSKDetFSI008 ) : return "f_skdetfsi008"; break;
        case ( kRHCNuisSKDetFSI009 ) : return "f_skdetfsi009"; break;
        case ( kRHCNuisSKDetFSI010 ) : return "f_skdetfsi010"; break;
        case ( kRHCNuisSKDetFSI011 ) : return "f_skdetfsi011"; break;
        case ( kRHCNuisSKDetFSI012 ) : return "f_skdetfsi012"; break;
        case ( kRHCNuisSKDetFSI013 ) : return "f_skdetfsi013"; break;
        case ( kRHCNuisSKDetFSI014 ) : return "f_skdetfsi014"; break;
        case ( kRHCNuisSKDetFSI015 ) : return "f_skdetfsi015"; break;
        case ( kRHCNuisSKDetFSI016 ) : return "f_skdetfsi016"; break;
        case ( kRHCNuisSKDetFSI017 ) : return "f_skdetfsi017"; break;
#elif defined(__USE_SKDETFSI2013__)
        case ( kRHCNuisSKDetFSI012 ) : return "f_skdetfsi012"; break;
        case ( kRHCNuisSKDetFSI013 ) : return "f_skdetfsi013"; break;
        case ( kRHCNuisSKDetFSI014 ) : return "f_skdetfsi014"; break;
        case ( kRHCNuisSKDetFSI015 ) : return "f_skdetfsi015"; break;
        case ( kRHCNuisSKDetFSI016 ) : return "f_skdetfsi016"; break;
        case ( kRHCNuisSKDetFSI017 ) : return "f_skdetfsi017"; break;
        case ( kRHCNuisSKDetFSI018 ) : return "f_skdetfsi018"; break;
        case ( kRHCNuisSKDetFSI019 ) : return "f_skdetfsi019"; break;
        case ( kRHCNuisSKDetFSI020 ) : return "f_skdetfsi020"; break;
        case ( kRHCNuisSKDetFSI021 ) : return "f_skdetfsi021"; break;
        case ( kRHCNuisSKDetFSI022 ) : return "f_skdetfsi022"; break;
        case ( kRHCNuisSKDetFSI023 ) : return "f_skdetfsi023"; break;
#endif
            
#endif
            // NEW SKDETFSI (start indices from 24 instead of 18)
#ifdef SKCORRELATIONS
        case ( kRHCNuisSKDetFSI024 ) : return "f_skdetfsi024"; break;
        case ( kRHCNuisSKDetFSI025 ) : return "f_skdetfsi025"; break;
        case ( kRHCNuisSKDetFSI026 ) : return "f_skdetfsi026"; break;
        case ( kRHCNuisSKDetFSI027 ) : return "f_skdetfsi027"; break;
        case ( kRHCNuisSKDetFSI028 ) : return "f_skdetfsi028"; break;
        case ( kRHCNuisSKDetFSI029 ) : return "f_skdetfsi029"; break;
        case ( kRHCNuisSKDetFSI030 ) : return "f_skdetfsi030"; break;
        case ( kRHCNuisSKDetFSI031 ) : return "f_skdetfsi031"; break;
        case ( kRHCNuisSKDetFSI032 ) : return "f_skdetfsi032"; break;
        case ( kRHCNuisSKDetFSI033 ) : return "f_skdetfsi033"; break;
        case ( kRHCNuisSKDetFSI034 ) : return "f_skdetfsi034"; break;
        case ( kRHCNuisSKDetFSI035 ) : return "f_skdetfsi035"; break;
        case ( kRHCNuisSKDetFSI036 ) : return "f_skdetfsi036"; break;
        case ( kRHCNuisSKDetFSI037 ) : return "f_skdetfsi037"; break;
#endif
            //#endif //SKDETFSITOGETHER
            
            
            //***************************** MULTIRING FORWARD HORN CURRENT *************************
            
            /// NEW MULTIRING 3
        case ( kMultiRingNuisSuperkEnergyScale     ) : return "f_sk_e_scale_multiring";        break;
            ///
            
            /// NEW MULTIRING 2
            
            //#ifdef SKFSITOGETHER
#if MR_FSI_BINS != 3
        case ( kMultiRingNuisSKDetFSI000 ) : return "f_skdetfsi000_multiring"; break;
        case ( kMultiRingNuisSKDetFSI001 ) : return "f_skdetfsi001_multiring"; break;
        case ( kMultiRingNuisSKDetFSI002 ) : return "f_skdetfsi002_multiring"; break;
        case ( kMultiRingNuisSKDetFSI003 ) : return "f_skdetfsi003_multiring"; break;
        case ( kMultiRingNuisSKDetFSI004 ) : return "f_skdetfsi004_multiring"; break;
        case ( kMultiRingNuisSKDetFSI005 ) : return "f_skdetfsi005_multiring"; break;
#endif
#if MR_FSI_BINS != 2
            
            // NEW SKDETFSI
#if defined(__USE_SKDETFSI2013__)
        case ( kMultiRingNuisSKDetFSI006 ) : return "f_skdetfsi006_multiring"; break;
        case ( kMultiRingNuisSKDetFSI007 ) : return "f_skdetfsi007_multiring"; break;
        case ( kMultiRingNuisSKDetFSI008 ) : return "f_skdetfsi008_multiring"; break;
        case ( kMultiRingNuisSKDetFSI009 ) : return "f_skdetfsi009_multiring"; break;
        case ( kMultiRingNuisSKDetFSI010 ) : return "f_skdetfsi010_multiring"; break;
        case ( kMultiRingNuisSKDetFSI011 ) : return "f_skdetfsi011_multiring"; break;
        case ( kMultiRingNuisSKDetFSI012 ) : return "f_skdetfsi012_multiring"; break;
        case ( kMultiRingNuisSKDetFSI013 ) : return "f_skdetfsi013_multiring"; break;
        case ( kMultiRingNuisSKDetFSI014 ) : return "f_skdetfsi014_multiring"; break;
        case ( kMultiRingNuisSKDetFSI015 ) : return "f_skdetfsi015_multiring"; break;
        case ( kMultiRingNuisSKDetFSI016 ) : return "f_skdetfsi016_multiring"; break;
        case ( kMultiRingNuisSKDetFSI017 ) : return "f_skdetfsi017_multiring"; break;
#endif
            
#endif
#ifdef SKCORRELATIONS
            // NEW SKDETFSI (start indices from 24 instead of 18)
        case ( kMultiRingNuisSKDetFSI024 ) : return "f_skdetfsi024_multiring"; break;
        case ( kMultiRingNuisSKDetFSI025 ) : return "f_skdetfsi025_multiring"; break;
        case ( kMultiRingNuisSKDetFSI026 ) : return "f_skdetfsi026_multiring"; break;
        case ( kMultiRingNuisSKDetFSI027 ) : return "f_skdetfsi027_multiring"; break;
        case ( kMultiRingNuisSKDetFSI028 ) : return "f_skdetfsi028_multiring"; break;
        case ( kMultiRingNuisSKDetFSI029 ) : return "f_skdetfsi029_multiring"; break;
        case ( kMultiRingNuisSKDetFSI030 ) : return "f_skdetfsi030_multiring"; break;
        case ( kMultiRingNuisSKDetFSI031 ) : return "f_skdetfsi031_multiring"; break;
        case ( kMultiRingNuisSKDetFSI032 ) : return "f_skdetfsi032_multiring"; break;
        case ( kMultiRingNuisSKDetFSI033 ) : return "f_skdetfsi033_multiring"; break;
        case ( kMultiRingNuisSKDetFSI034 ) : return "f_skdetfsi034_multiring"; break;
        case ( kMultiRingNuisSKDetFSI035 ) : return "f_skdetfsi035_multiring"; break;
        case ( kMultiRingNuisSKDetFSI036 ) : return "f_skdetfsi036_multiring"; break;
        case ( kMultiRingNuisSKDetFSI037 ) : return "f_skdetfsi037_multiring"; break;
#endif
            //#endif //SKDETFSITOGETHER
            
            ///
            
        default :
            return AsStringCompact(p);
    }
    return AsStringCompact(p);
}

//________________________________________________________________________________________
bool valor::util::IsMixingAngle(PhysParamId_t p)
{
    switch(p) {
        case ( kPhypMixing23              ) : return true; break;
        case ( kPhypMixing13              ) : return true; break;
        case ( kPhypMixing12              ) : return true; break;
        case ( kPhypSinsq2ThetaDecay      ) : return true; break;
        case ( kPhypSinsq2ThetaDecoh      ) : return true; break;
        case ( kPhypMixing14              ) : return true; break;
        case ( kPhypMixing15              ) : return true; break;
        case ( kPhypMixing24              ) : return true; break;
        case ( kPhypMixing25              ) : return true; break;
        case ( kPhypMixing34              ) : return true; break;
        case ( kPhypMixing35              ) : return true; break;
            
        default :
            return false;
    }
    return false;
}
//________________________________________________________________________________________
bool valor::util::IsSystematic(PhysParamId_t p)
{
    if(p <= kPhypNSIepsMuTau) return false;
    return true;
}
//________________________________________________________________________________________

bool valor::util::IsFHCSystematic(PhysParamId_t p)
{
    //return !IsRHCSystematic(p);
    
    /// NEW MULTIRING 2
    if(p < kNuisRHCToFHC)
        return true;
    else
        return false;
    ///
}
//________________________________________________________________________________________
bool valor::util::IsBANFFUnconstrained(PhysParamId_t p)
{
    switch (p) {
      //case (kNuisBANFFNorm2p2h        )   :
	  //case (kNuisBANFFFermiMom   )   :
	  //case (kNuisBANFFBindE      )   :
	  //case (kNuisBANFFNormCCcoh  )   :
        case (kNuisBANFFNormNuEToNuMu) :
        case (kNuisBANFFNormNuEBarToNuMuBar) :
        case (kNuisBANFFNormNC1Gamma   )   :
        case (kNuisBANFFNormNCoth   )   :
        case (kNuisBANFFEb) :
        case (kNuisBANFFMultiReNorm) :
	  //case (kNuisBANFFNormNuEBarToNuMuBar) :
	  //case (kNuisBANFFNorm2p2hBar     )   :
            
	  //case (kRHCNuisBANFF2p2h        )   :
	  //case (kRHCNuisBANFFFermiMom   )   :
	  //case (kRHCNuisBANFFBindE      )   :
	  //case (kRHCNuisBANFFNormCCcoh  )   :
        case (kRHCNuisBANFFNormNuEToNuMu) :
        case (kRHCNuisBANFFNormNuEBarToNuMuBar) :
        case (kRHCNuisBANFFNormNC1Gamma   )   :
        case (kRHCNuisBANFFNormNCoth   )   :
        case (kRHCNuisBANFFEb ) :
	  //case (kRHCNuisBANFFNormNuEBarToNuMuBar) :
	  //case (kRHCNuisBANFF2p2hBar     )   :
            
            
            return true;
            
        default:
            return false;
    }
    return false;
}
//________________________________________________________________________________________
bool valor::util::IsBANFFNotUnconstrained(PhysParamId_t p)
{
  return !IsBANFFUnconstrained(p);
}
//________________________________________________________________________________________
bool valor::util::IsRPA(PhysParamId_t p)
{
    switch (p) {
        case (kNuisXSecRPA        )   :
        case (kRHCNuisXSecRPA     )   :
            return true;
            
        default:
            return false;
    }
    return false;
}
//________________________________________________________________________________________
bool valor::util::IsBeRPA(PhysParamId_t p)
{
  switch (p) {

    case ( kNuisBANFFBeRPA_A    ) :
    case ( kNuisBANFFBeRPA_B    ) :
    case ( kNuisBANFFBeRPA_D    ) :
    case ( kNuisBANFFBeRPA_E    ) :
    case ( kNuisBANFFBeRPA_U    ) :      

    case ( kRHCNuisBANFFBeRPA_A ) :   
    case ( kRHCNuisBANFFBeRPA_B ) :  
    case ( kRHCNuisBANFFBeRPA_D ) :  
    case ( kRHCNuisBANFFBeRPA_E ) :
    case ( kRHCNuisBANFFBeRPA_U ) :        
      return true;
            
  default:
    return false;
  }
  
}
//________________________________________________________________________________________
bool valor::util::IsRHCSystematic(PhysParamId_t p)
{
    // if(!IsSystematic(p))
    //   return false;
    // if(p < kNuisRHCToFHC)
    //   return false;
    // else if(p < kRHCLastParam)
    //   return true;
    // std::cerr<<"PhysParamId/IsRHCSystematic >> Cannot determine whether or not RHC syst: "<<AsString(p)
    //   <<". Assuming no."<<std::endl;
    // return false;
    

    /// NEW MULTIRING 2
    if(!IsSystematic(p))
        return false;
    if(p < kNuisRHCToFHC || p > kRHCLastParam)
        return false;
    else if(p < kRHCLastParam)
        return true;
    std::cerr<<"PhysParamId/IsRHCSystematic >> Cannot determine whether or not RHC syst: "<<AsString(p)
    <<". Assuming no."<<std::endl;
    return false;
    ///
    
}
//________________________________________________________________________________________

/// NEW MULTIRING 2
bool valor::util::IsMultiRingSystematic(PhysParamId_t p)
{
    if(!IsSystematic(p))
        return false;
    if((p > kNuisXSecNormAntiNuToNu && p < kMultiRingNuisSuperkEnergyScale) || p < kNuisBANFF00)  // BANFF and Uncorr FHC syst are MultiRing as well (no Escale & SKDetFSI
        return false;
    else if(p < kMultiRingLastParam)
        return true;
    std::cerr<<"PhysParamId/IsMultiRingSystematic >> Cannot determine whether or not MultiRing syst: "<<AsString(p)
    <<". Assuming no."<<std::endl;
    return false;
}
//________________________________________________________________________________________
bool valor::util::IsXSec(PhysParamId_t p)
{
  // N.B. Leave the old uncorrelated systematics that will not be used, since now they are all banff (kNuisBANFF*)
    
  //Old uncorr xsec syst (as in the old code)
    
#if ( defined(__USE_NIWG2015__) ||defined(__USE_NIWG2017__) )
  if(p == kNuisXSecRPA || p == kRHCNuisXSecRPA)
    return true;
#endif

  switch(p) {
            
  case ( kNuisBANFFNorm2p2h         ) : return true;
  case ( kNuisBANFFNorm2p2hBar      ) : return true;
  case ( kNuisBANFFCA5              ) : return true;
  case ( kNuisBANFFBgRES            ) : return true;
  case ( kNuisBANFFMAQE             ) : return true;
  case ( kNuisBANFFMARES            ) : return true;
#if (BANFF_BINS != 2)                 
  case ( kNuisBANFFSCCA             ) : return true;
  case ( kNuisBANFFSCCV             ) : return true;
#endif                                
  case ( kNuisBANFFFermiMom         ) : return true;
  case ( kNuisBANFFShapeCCoth       ) : return true;
    //case ( kNuisBANFFBindE            ) : return true;
  case ( kNuisBANFFNormCCcoh        ) : return true;
  case ( kNuisBANFFNormNCcoh        ) : return true;
  case ( kNuisBANFFNormNCoth        ) : return true;
  case ( kNuisBANFFNormNC1Gamma     ) : return true;
         
  case ( kNuisBANFFNormNuEToNuMu      ) : return true;
  case ( kNuisBANFFNormNuEBarToNuMuBar) : return true;
    //case (kNuisBANFFNormNuEToNuMu    ) : return true;
    //case (kNuisBANFFNormNuEBarToNuMuBar) : return true;
  case ( kNuisBANFFBeRPA_A        ) : return true;
  case ( kNuisBANFFBeRPA_B        ) : return true;
  case ( kNuisBANFFBeRPA_D        ) : return true;
  case ( kNuisBANFFBeRPA_E        ) : return true;
  case ( kNuisBANFFBeRPA_U        ) : return true;
  case ( kNuisBANFFShape2p2h      ) : return true;
  case ( kNuisBANFFNorm2p2hCtoO   ) : return true;
  case ( kNuisBANFFEb   ) : return true;
  case ( kNuisBANFFMultiReNorm    ) : return true;

    
    
  case ( kRHCNuisBANFF2p2hBar       ) : return true;
  case ( kRHCNuisBANFF2p2h          ) : return true;
  case ( kRHCNuisBANFFCA5           ) : return true;
  case ( kRHCNuisBANFFBgRES         ) : return true;
  case ( kRHCNuisBANFFMAQE          ) : return true;
  case ( kRHCNuisBANFFMARES         ) : return true;
#if (BANFF_BINS != 2)                 
  case ( kRHCNuisBANFFSCCA          ) : return true;
  case ( kRHCNuisBANFFSCCV          ) : return true;
#endif   
  case ( kRHCNuisBANFFFermiMom      ) : return true;
  case ( kRHCNuisBANFFShapeCCoth    ) : return true;
    //case ( kRHCNuisBANFFBindE         ) : return true;
  case ( kRHCNuisBANFFNormCCcoh     ) : return true;
  case ( kRHCNuisBANFFNormNCcoh     ) : return true;
  case ( kRHCNuisBANFFNormNCoth     ) : return true;
  case ( kRHCNuisBANFFNormNuEToNuMu ) : return true;
  case ( kRHCNuisBANFFNormNuEBarToNuMuBar) : return true;
  case ( kRHCNuisBANFFNormNC1Gamma  ) : return true;

  case ( kRHCNuisBANFFBeRPA_A       ) : return true;
  case ( kRHCNuisBANFFBeRPA_B       ) : return true;
  case ( kRHCNuisBANFFBeRPA_D       ) : return true;
  case ( kRHCNuisBANFFBeRPA_E       ) : return true;
  case ( kRHCNuisBANFFBeRPA_U       ) : return true;
  case ( kRHCNuisBANFFShape2p2h     ) : return true;
  case ( kRHCNuisBANFFNorm2p2hCtoO  ) : return true;
  case ( kRHCNuisBANFFEb  ) : return true;
  case ( kRHCNuisBANFFMultiReNorm   ) : return true;

  default                           : return false;
  }
  return false;
}
//________________________________________________________________________________________
bool valor::util::IsSCC(PhysParamId_t p)
{
#if (defined(__USE_NIWG2015__) && (BANFF_BINS != 2) ) || defined (__USE_NIWG2017__)
  if (p == kNuisBANFFSCCA        ||
      p == kNuisBANFFSCCV        ||
      p == kRHCNuisBANFFSCCA     || 
      p == kRHCNuisBANFFSCCV     ) {
    return true;
  }
#endif
  
  return false;
}
//____________________________________________________________________________
bool valor::util::IsNormXSecParam(conv::PhysParamId_t p)
{
  //Determine if the BANFF nuisance parameter is normalisation
  switch (p) {
    case( kNuisBANFFNorm2p2h            ) : return true; break;
    case( kNuisBANFFNorm2p2hBar         ) : return true; break;

    case( kNuisBANFFNormCCcoh           ) : return true; break;
    case( kNuisBANFFNormNCcoh           ) : return true; break;
    case( kNuisBANFFNormNCoth           ) : return true; break;
    case( kNuisBANFFNormNC1Gamma        ) : return true; break;
    case( kNuisBANFFNormNuEToNuMu       ) : return true; break;
    case( kNuisBANFFNormNuEBarToNuMuBar ) : return true; break;
    case( kNuisBANFFBeRPA_U             ) : return true; break;
    case( kNuisBANFFNorm2p2hCtoO        ) : return true; break;      

    case( kRHCNuisBANFF2p2h               ) : return true; break;
    case( kRHCNuisBANFF2p2hBar            ) : return true; break;

    case( kRHCNuisBANFFNormCCcoh          ) : return true; break;
    case( kRHCNuisBANFFNormNCcoh          ) : return true; break;
    case( kRHCNuisBANFFNormNCoth          ) : return true; break;
    case( kRHCNuisBANFFNormNC1Gamma       ) : return true; break;
    case( kRHCNuisBANFFNormNuEToNuMu      ) : return true; break;
    case( kRHCNuisBANFFNormNuEBarToNuMuBar) : return true; break; 
    case( kRHCNuisBANFFBeRPA_U            ) : return true; break;
    case( kRHCNuisBANFFNorm2p2hCtoO       ) : return true; break;
      
    default                               : return false; break;
   }
  return false;

}
//____________________________________________________________________________
bool valor::util::IsXSecNumuNue(conv::PhysParamId_t p)
{
  switch (p) {
    case ( kNuisBANFFNormNuEToNuMu       ) : return true; break;
    case ( kNuisBANFFNormNuEBarToNuMuBar ) : return true; break;
    default : return false; break;    
  }
  return false;
}
//____________________________________________________________________________
bool valor::util::IsEScale(conv::PhysParamId_t p)
{
  switch (p) {
    case (kNuisSuperkEnergyScale) : return true; break;
    default : return false; break;
  }
  return false;
}
//____________________________________________________________________________
bool valor::util::IsUsingSplines(PhysParamId_t p)
{
  if(p >= kNuisXSecFermiMom && p <= kNuisXSecBindE)
    return true;
  if(p >= kRHCNuisXSecFermiMom && p <= kRHCNuisXSecBindE)
    return true;
    
  switch (p) {
  case kNuisBANFFMAQE:
  case kNuisBANFFMARES:
#if (BANFF_BINS != 2)
  case kNuisBANFFSCCA:
  case kNuisBANFFSCCV:
#endif
  case kNuisBANFFCA5:
  case kNuisBANFFShapeCCoth:
  case kNuisBANFFFermiMom:
    //case kNuisBANFFBindE:
  case kNuisBANFFBgRES:
            
  case kNuisBANFFBeRPA_A:     
  case kNuisBANFFBeRPA_B:    
  case kNuisBANFFBeRPA_D:     
  case kNuisBANFFBeRPA_E:     
  case kNuisBANFFShape2p2h:   
  case kNuisBANFFEb:
    return true;
            
  case kRHCNuisBANFFMAQE:
  case kRHCNuisBANFFMARES:
#if (BANFF_BINS != 2)
  case kRHCNuisBANFFSCCA:
  case kRHCNuisBANFFSCCV:
#endif
  case kRHCNuisBANFFCA5:
  case kRHCNuisBANFFShapeCCoth:
  case kRHCNuisBANFFFermiMom:
    //case kRHCNuisBANFFBindE:
  case kRHCNuisBANFFBgRES:

  case kRHCNuisBANFFBeRPA_A:     
  case kRHCNuisBANFFBeRPA_B:    
  case kRHCNuisBANFFBeRPA_D:     
  case kRHCNuisBANFFBeRPA_E:     
  case kRHCNuisBANFFShape2p2h:   
  case kRHCNuisBANFFEb:
    return true;
            
  default:
    return false;
  }
  return false;
}
//____________________________________________________________________________
bool valor::util::IsEff1Rmu(PhysParamId_t p)
{
    assert(!"IsEff1Re Currently unimplemented");
    switch(p) {
        default :
            return false;
    }
    return false;
}
//____________________________________________________________________________
bool valor::util::IsEff1Re(PhysParamId_t p)
{
    assert(!"IsEff1Re Currently unimplemented");
    switch(p) {
        default :
            return false;
    }
    return false;
}
//____________________________________________________________________________
bool valor::util::IsEff(PhysParamId_t p)
{
    assert(!"IsEff Currently unimplemented");
    switch(p) {
        case ( kNuisSuperkEnergyScale     ) : return true;  break;
        default :
            return false;
    }
    return false;
}
//____________________________________________________________________________
bool valor::util::IsFSI(PhysParamId_t p)
{
    assert(!"IsFSI Currently unimplemented");
    switch(p) {
        default :
            return false;
    }
    return false;
}
//____________________________________________________________________________
bool valor::util::IsSKDetFSI(PhysParamId_t p)
{
    //fhc
#ifdef SKCORRELATIONS
    // NEW SKDETFSI (increase number)
    if(p >= kNuisSKDetFSI000 && p <= kNuisSKDetFSI037)
        return true;
#else
#if FSI_BINS != 3
    
    // NEW SKDETFSI
#if defined(__USE_SKDETFSI2013__)
    if(p >= kNuisSKDetFSI000 && p <= kNuisSKDetFSI005)
        return true;
#elif defined(__USE_SKDETFSI2015__)
    if(p >= kNuisSKDetFSI000 && p <= kNuisSKDetFSI011)
        return true;
#endif
    
#endif
#if FSI_BINS != 2
    
    // NEW SKDETFSI
#if defined(__USE_SKDETFSI2013__)
    if(p >= kNuisSKDetFSI006 && p <= kNuisSKDetFSI017)
        return true;
#elif defined(__USE_SKDETFSI2015__)
    if(p >= kNuisSKDetFSI012 && p <= kNuisSKDetFSI023)
        return true;
#endif
    
#endif
#endif
    if(p == kNuisSuperkEnergyScale)
        return true;
    
    //rhc
#ifdef SKCORRELATIONS
    // NEW SKDETFSI (increase number)
    if(p >= kRHCNuisSKDetFSI000 && p <= kRHCNuisSKDetFSI037)
        return true;
#else
#if FSI_BINS != 3
    
    // NEW SKDETFSI
#if defined(__USE_SKDETFSI2013__)
    if(p >= kRHCNuisSKDetFSI000 && p <= kRHCNuisSKDetFSI005)
        return true;
#elif defined(__USE_SKDETFSI2015__)
    if(p >= kRHCNuisSKDetFSI000 && p <= kRHCNuisSKDetFSI011)
        return true;
#endif
    
#endif
#if FSI_BINS != 2
    
    // NEW SKDETFSI
#if defined(__USE_SKDETFSI2013__)
    if(p >= kRHCNuisSKDetFSI006 && p <= kRHCNuisSKDetFSI017)
        return true;
#elif defined(__USE_SKDETFSI2015__)
    if(p >= kRHCNuisSKDetFSI012 && p <= kRHCNuisSKDetFSI023)
        return true;
#endif
    
#endif
#endif
    if(p == kRHCNuisSuperkEnergyScale)
        return true;
    
    
    /// NEW MULTIRING 2
    //multiring
#ifdef SKCORRELATIONS
    // NEW SKDETFSI (shift number)
    if(p >= kMultiRingNuisSKDetFSI000 && p <= kMultiRingNuisSKDetFSI037)
        return true;
#else
#if FSI_BINS != 3
    
    // NEW SKDETFSI
#if defined(__USE_SKDETFSI2013__)
    if(p >= kMultiRingNuisSKDetFSI_first && p <=  (PhysParamId_t)kListOfMultiRingSKDetFSIParams[conv::kNumOfSKDetFSIParamsMultiRing-1])
        return true;
#elif defined(__USE_SKDETFSI2015__)
    if(p >= kMultiRingNuisSKDetFSI_first && p <=  (PhysParamId_t)kListOfMultiRingSKDetFSIParams[conv::kNumOfSKDetFSIParamsMultiRing-1])
        return true;
#endif
    
#endif
#if FSI_BINS != 2
    
    // NEW SKDETFSI
#if defined(__USE_SKDETFSI2013__)
    if(p >= kMultiRingNuisSKDetFSI_first && p <= (PhysParamId_t)kListOfMultiRingSKDetFSIParams[conv::kNumOfSKDetFSIParamsMultiRing-1])
        return true;
#elif defined(__USE_SKDETFSI2015__)
    if(p >= kMultiRingNuisSKDetFSI_first && p <= (PhysParamId_t)kListOfMultiRingSKDetFSIParams[conv::kNumOfSKDetFSIParamsMultiRing-1])
        return true;
#endif
    
#endif
#endif
    // MultiRingNuisEnergy scale isn't used instead have 
    if(p == kNuisSuperkEnergyScale)   
    return true;
    
    
    /// NEW MULTIRNG 3
    if(p == kMultiRingNuisSuperkEnergyScale)
        return true;
    
    return false;
}
//____________________________________________________________________________
bool valor::util::IsFlux(PhysParamId_t p)
{
    if ( (p >= kNuisBANFF00    && p <= kNuisBANFF24) ||
        (p >= kRHCNuisBANFF00 && p <= kRHCNuisBANFF24) ) {
        return true;
    }
    
    return false;
}
//____________________________________________________________________________
bool valor::util::IsEnergyScale(PhysParamId_t p)
{
    if (p == kNuisSuperkEnergyScale          ||
        p == kRHCNuisSuperkEnergyScale       ||
        p == kMultiRingNuisSuperkEnergyScale)
        return true;
    return false;
}
//____________________________________________________________________________
bool valor::util::IsNormalisation(PhysParamId_t p)
{
    if (util::IsEnergyScale(p)                                              ) return false;
    
    if (p==kNuisRHCToFHC                                                    ) return true;
    
    if (util::IsSKDetFSI(p)                                                 ) return true;
    if (util::IsFlux(p)                                                     ) return true;
    if ( (util::IsXSec(p) || util::IsBANFF(p) ) && !util::IsUsingSplines(p) ) return true;
    return false;
}
//____________________________________________________________________________
bool valor::util::IsBANFF(PhysParamId_t p)
{
    //fhc
    if(p >= kNuisBANFF00 && p < kLastBANFFParam)
        return true;
    
    if(p >= kRHCNuisBANFF00 && p < kRHCLastBANFFParam)
        return true;
    
    return false;
}
//____________________________________________________________________________
bool valor::util::IsRHCBANFF(PhysParamId_t p)
{
  if(p >= kRHCNuisBANFF00 && p < kRHCLastBANFFParam) return true;
  return false;
}
//____________________________________________________________________________
string valor::util::AsStringSpline(PhysParamId_t p)
{
  switch(p) {
  case ( kNuisBANFFMAQE           ) :
  case ( kRHCNuisBANFFMAQE        ) :
    return "maqe";
  case ( kNuisBANFFMARES          ) :
  case ( kRHCNuisBANFFMARES       ) :
    return "manff";
            
  // NEW XSEC
  case ( kNuisBANFFCA5            ) :
  case ( kRHCNuisBANFFCA5         ) :
    return "ca5";
  case ( kNuisBANFFBgRES          ) :
  case ( kRHCNuisBANFFBgRES       ) :
    return "bgscl";
  // NEW SPLIT ANA
#if (BANFF_BINS != 2)
  case ( kNuisBANFFSCCV           ) :
  case ( kRHCNuisBANFFSCCV        ) :
    return "sccvec";
  case ( kNuisBANFFSCCA           ) :
  case ( kRHCNuisBANFFSCCA        ) :
    return "sccaxl";
#endif
  case ( kNuisBANFFFermiMom       ) :
  case ( kRHCNuisBANFFFermiMom    ) :
    return "pfo";
  case ( kNuisBANFFShapeCCoth     ) :
  case ( kRHCNuisBANFFShapeCCoth  ) :
    return "dismpishp";
  case ( kNuisBANFFEb          ) :
  case ( kRHCNuisBANFFEb       ) :
    return "ebo";
  //case ( kNuisBANFFBindE          ) :
  //case ( kRHCNuisBANFFBindE       ) :
  //return "ebo";
  //
            
  // NEW XSEC
  case ( kNuisXSecRPA             ) :
  case ( kRHCNuisXSecRPA          ) :
    return "rpa";
            
  case ( kNuisXSecBindE           ) :
  case ( kRHCNuisXSecBindE        ) :
    return "ebo";
  case ( kNuisXSecFermiMom        ) :
  case ( kRHCNuisXSecFermiMom     ) :
    return "pfo";
  case ( kNuisXSecShapeCCoth      ) :
  case ( kRHCNuisXSecShapeCCoth   ) :
    return "dismpishp";
  case ( kNuisXSecPiLessDelta     ) :
  case ( kRHCNuisXSecPiLessDelta  ) :
    return "piless";
  case ( kNuisXSecSpecFunc        ) :
  case ( kRHCNuisXSecSpecFunc     ) :
    return "sf";
  case ( kNuisXSecWshape          ) :
  case ( kRHCNuisXSecWshape       ) :
    return "wshape";
            
  case ( kNuisBANFFBeRPA_A         ) : 
  case ( kRHCNuisBANFFBeRPA_A      ) : 
    return "erpa_A";    break;
  case ( kNuisBANFFBeRPA_B         ) : 
  case ( kRHCNuisBANFFBeRPA_B      ) : 
    return "erpa_B";    break;
  case ( kNuisBANFFBeRPA_D         ) : 
  case ( kRHCNuisBANFFBeRPA_D      ) : 
    return "erpa_D";    break;
  case ( kNuisBANFFBeRPA_E         ) : 
  case ( kRHCNuisBANFFBeRPA_E      ) : 
    return "erpa_E";    break;
  case ( kNuisBANFFBeRPA_U         ) : 
  case ( kRHCNuisBANFFBeRPA_U      ) : 
    return "erpa_U";    break;
  case ( kNuisBANFFShape2p2h       ) : 
  case ( kRHCNuisBANFFShape2p2h    ) : 
    return "mecpdd";     break;

  default:
    break;
  }
  std::cerr<<"util/AsStringSpline failed for parameter: "<<util::AsString(p)
           <<" ("<<util::RawParamAsString(p)<<")"<<std::endl;
  return "--";
}
//____________________________________________________________________________
string valor::util::AsStringCompact_KDI(PhysParamId_t p)
{
    //Todo: update with new systematics
    
    //Uses the KDI conventions (where different)
    
    switch(p) {
            
        case ( kNuisBANFFMAQE             ) : return "maqe";     break;
        case ( kNuisBANFFMARES            ) : return "mares";    break;
        case ( kNuisXSecFermiMom          ) : return "pf";       break;
        case ( kNuisXSecSpecFunc          ) : return "sf";       break;
        case ( kNuisXSecWshape            ) : return "mdel";     break;
        case ( kNuisXSecShapeCCoth        ) : return "ccothshp"; break;
        case ( kNuisXSecPiLessDelta       ) : return "pdd";      break;
        case ( kNuisXSecBindE             ) : return "eb";       break;
            
        default :
            return AsStringCompact(p);
    }
    return "--";
}
//____________________________________________________________________________
string valor::util::AsStringSimple(PhysParamId_t p)
{
    if      (p==conv::kPhypMixing12   ) { return "mixing12"   ;}
    else if (p==conv::kPhypMixing23   ) { return "mixing23"   ;}
    else if (p==conv::kPhypMixing13   ) { return "mixing13"   ;}
    else if (p==conv::kPhypCPV        ) { return "deltacp"    ;}
    else if (p==conv::kPhypDm12sq     ) { return "dmsq12"     ;}
    else if (p==conv::kPhypDm23sq     ) { return "dmsq23"     ;}
    
    else if (p==conv::kPhypMixing12Bar) { return "mixing12bar";}
    else if (p==conv::kPhypMixing23Bar) { return "mixing23bar";}
    else if (p==conv::kPhypMixing13Bar) { return "mixing13bar";}
    else if (p==conv::kPhypCPVBar     ) { return "deltacpbar" ;}
    else if (p==conv::kPhypDm12sqBar  ) { return "dmsq12bar"  ;}
    else if (p==conv::kPhypDm23sqBar  ) { return "dmsq23bar"  ;}
    
    // NEW BETABAR
    else if (p==conv::kPhypBetaBar    ) { return "betabar"    ;}
    //
    
    return util::AsStringCompact(p);
}
//____________________________________________________________________________
string valor::util::AsStringNIWG(PhysParamId_t p)
{    
  switch(p) {
            
  case (kNuisBANFFNorm2p2h         ): return "2p2h_norm_nu";
  case (kNuisBANFFNorm2p2hBar      ): return "2p2h_norm_nubar";
  case (kNuisBANFFCA5              ): return "CA5";
  case (kNuisBANFFBgRES            ): return "ISO_BKG";
  case (kNuisBANFFMAQE             ): return "MAQE";
  case (kNuisBANFFMARES            ): return "MARES";
#if (BANFF_BINS != 2)
  case (kNuisBANFFSCCA             ): return "SCCA";
  case (kNuisBANFFSCCV             ): return "SCCV";
#endif
  case (kNuisBANFFFermiMom         ): return "pF_O";
  case (kNuisBANFFShapeCCoth       ): return "CC_DIS";
    //case (kNuisBANFFBindE            ): return "Eb_O_nu";
  case (kNuisBANFFNormCCcoh        ): return "CC_Coh_O";
  case (kNuisBANFFNormNCcoh        ): return "NC_Coh";
  case (kNuisBANFFNormNCoth        ): return "NC_other_far";
  case (kNuisBANFFNormNC1Gamma     ): return "NC_1gamma";

  case (kNuisBANFFNormNuEToNuMu      ): return "nue_numu";
  case (kNuisBANFFNormNuEBarToNuMuBar): return "nuebar_numubar";
 
  case (kNuisBANFFBeRPA_A          ) :return "BeRPA_A";    
  case (kNuisBANFFBeRPA_B          ) :return "BeRPA_B";   
  case (kNuisBANFFBeRPA_D          ) :return "BeRPA_D";    
  case (kNuisBANFFBeRPA_E          ) :return "BeRPA_E";    
  case (kNuisBANFFBeRPA_U          ) :return "BeRPA_U";    
  case (kNuisBANFFShape2p2h        ) :return "2p2h_shape_O";  
  case (kNuisBANFFNorm2p2hCtoO     ) :return "2p2h_normCtoO";
  case (kNuisBANFFEb     ) :return "ebo";
  case (kNuisBANFFMultiReNorm      ) :return "FHC_MultiRe_norm";

  case (kRHCNuisBANFF2p2h          ): return "2p2h_norm_O";
  case (kRHCNuisBANFF2p2hBar       ): return "2p2h_norm_nubar_O";
  case (kRHCNuisBANFFCA5           ): return "CA5";
  case (kRHCNuisBANFFBgRES         ): return "ISO_BKG";
  case (kRHCNuisBANFFMAQE          ): return "MAQE";
  case (kRHCNuisBANFFMARES         ): return "MARES";
#if (BANFF_BINS != 2)
  case (kRHCNuisBANFFSCCA          ): return "SCCA";
  case (kRHCNuisBANFFSCCV          ): return "SCCV";
#endif
  case (kRHCNuisBANFFFermiMom      ): return "pF_O";
  case (kRHCNuisBANFFShapeCCoth    ): return "CC_DIS";
    //case (kRHCNuisBANFFBindE         ): return "Eb_O_nu";
  case (kRHCNuisBANFFNormCCcoh     ): return "CC_Coh_O";
  case (kRHCNuisBANFFNormNCcoh     ): return "NC_Coh";
  case (kRHCNuisBANFFNormNCoth     ): return "NC_other_far";
  case (kRHCNuisBANFFNormNuEToNuMu ): return "nue_numu";
  case (kRHCNuisBANFFNormNuEBarToNuMuBar): return "nuebar_numubar";
  case (kRHCNuisBANFFNormNC1Gamma  ): return "NC_1gamma";

  case (kRHCNuisBANFFBeRPA_A       ) :return "BeRPA_A";    
  case (kRHCNuisBANFFBeRPA_B       ) :return "BeRPA_B";   
  case (kRHCNuisBANFFBeRPA_D       ) :return "BeRPA_D";    
  case (kRHCNuisBANFFBeRPA_E       ) :return "BeRPA_E";    
  case (kRHCNuisBANFFBeRPA_U       ) :return "BeRPA_U";    
  case (kRHCNuisBANFFShape2p2h     ) :return "2p2h_shape_O_nu";  
  case (kRHCNuisBANFFNorm2p2hCtoO  ) :return "2p2h_normCtoO";
  case (kRHCNuisBANFFEb     ) :return "ebo";
  case (kRHCNuisBANFFMultiReNorm   ) :return "RHC_MultiRe_norm_UNUSED";

  default                           : return "--";
  }
}
//____________________________________________________________________________
valor::conv::PhysParamId_t
valor::util::PhysParamIdFromString(string sh)
{
    string lsh = util::ToLowerCase(sh);
    
    if (lsh == "dmsq"   ) return conv::kPhypDm23sq;
    if (lsh == "dmsqbar") return conv::kPhypDm23sqBar;
    
    for (int i = 0 ; i < knppmax ; i++) {
        PhysParamId_t k = (PhysParamId_t) i ;
        string s = util::ToLowerCase(AsStringCompact(k));
        string s2;
        if (s.find("f_") == 0) s2 = s.substr(2);
        else s2 = "Bad string";
        string s3 = util::ToLowerCase(AsStringSimple(k));
        string s4 = util::ToLowerCase(AsStringPaperRooTex(k));
        if (k != kPhypUndefined && (lsh == s || lsh == s2 || lsh==s3|| lsh==s4) ) return k;
    }
    std::cerr<<"util::PhysParamIdFromString >> Unmatched string: "<<sh<<std::endl;
    return kPhypUndefined;
}
//____________________________________________________________________________
valor::conv::PhysParamId_t
valor::util::PhysParamIdFromLongString(string sh)
{
    string lsh = util::ToLowerCase(sh);
    for (int i = 0 ; i < knppmax ; i++) {
        PhysParamId_t k = (PhysParamId_t) i ;
        if (k != kPhypUndefined && lsh == util::ToLowerCase(AsString(k)) ) return k;
    }
    std::cerr<<"util::PhysParamIdFromString >> Unmatched string: "<<sh<<std::endl;
    return kPhypUndefined;
}
//____________________________________________________________________________
valor::conv::PhysParamId_t
valor::util::PhysParamIdFromSplineString(string sh)
{
  string lsh = util::ToLowerCase(sh);
  for (int i = 0 ; i < conv::knppmax ; i++) {
    conv::PhysParamId_t k = (conv::PhysParamId_t) i ;
    if (!util::IsUsingSplines(k)) continue;
    if (k != conv::kPhypUndefined && lsh == util::ToLowerCase(util::AsStringSpline(k)) ) return k;
  }
  std::cerr<<"util::PhysParamIdFromString >> Unmatched string: "<<sh<<std::endl;
  return kPhypUndefined;
}
//____________________________________________________________________________
valor::conv::PhysParamId_t
valor::util::GetBarredVersionOfPhysParam(valor::conv::PhysParamId_t p)
{
    switch(p) {
        case conv::kPhypMixing12: return conv::kPhypMixing12Bar;
        case conv::kPhypMixing13: return conv::kPhypMixing13Bar;
        case conv::kPhypMixing23: return conv::kPhypMixing23Bar;
        case conv::kPhypDm12sq  : return conv::kPhypDm12sqBar;
        case conv::kPhypDm23sq  : return conv::kPhypDm23sqBar;
        case conv::kPhypCPV     : return conv::kPhypCPVBar;
            
        default: return kPhypUndefined;
    }
    return kPhypUndefined;
}
//____________________________________________________________________________
string valor::util::ParamAsString                   (conv::PhysParamId_t p) { return util::AsString                  (p);}
string valor::util::ParamAsStringTex                (conv::PhysParamId_t p) { return util::AsStringTex               (p);}
string valor::util::ParamAsStringPaperTexNoDollar   (conv::PhysParamId_t p) { return util::AsStringPaperTexNoDollar  (p);}
string valor::util::ParamAsStringPaperTex           (conv::PhysParamId_t p) { return util::AsStringPaperTex          (p);}
string valor::util::ParamAsStringPaperRooTex        (conv::PhysParamId_t p) { return util::AsStringPaperRooTex       (p);}
string valor::util::ParamAsStringPaperTable         (conv::PhysParamId_t p) { return util::AsStringPaperTable        (p);}
string valor::util::ParamAsStringSimple             (conv::PhysParamId_t p) { return util::AsStringSimple            (p);}
string valor::util::ParamAsStringCompact            (conv::PhysParamId_t p) { return util::AsStringCompact           (p);}
string valor::util::ParamAsStringCompact_splinehack (conv::PhysParamId_t p) { return util::AsStringCompact_splinehack(p);}
string valor::util::ParamAsStringCompact_KDI        (conv::PhysParamId_t p) { return util::AsStringCompact_KDI       (p);}
string valor::util::ParamAsStringSpline             (conv::PhysParamId_t p) { return util::AsStringSpline            (p);}
string valor::util::ParamAsStringNIWG               (conv::PhysParamId_t p) { return util::AsStringNIWG              (p);}
//____________________________________________________________________________
void valor::util::PopulateNameRegistry(void)
{
  NameRegistry * nr = util::NameRegistry::MutableInstance();
  
  for (int ip = -1 ; ip < conv::kLastNuisParam; ip++) {
    conv::PhysParamId_t p = (conv::PhysParamId_t)ip;
    nr->fParamNames      [p] = util::AsString(p);
    nr->fParamNamesSimple[p] = util::AsStringSimple(p);
  }
  return;
}
//____________________________________________________________________________
string
valor::util::AsStringOfficialAxis(valor::conv::PhysParamId_t p)
{
  switch(p) {
    case conv::kPhypMixing12: return "sin^{2}(#theta_{12})";
    case conv::kPhypMixing13: return "sin^{2}(#theta_{13})";
    case conv::kPhypMixing23: return "sin^{2}(#theta_{23})";
    case conv::kPhypDm12sq  : return "|#Delta m^2_{12}| (eV^{2}c^{-4})";
    case conv::kPhypDm23sq  :
    {
      NuOscDmsqScheme_t omp = Opt::Instance()->GetEffectiveDmsqScheme();
      switch(omp) {
        case(kNuDmsqFogliLisi) :        return "|#Delta m^{2}| = |m^{2}_{3} - (m^{2}_{1} + m^{2}_{1})/2| (eV^{2}c^{-4})";                  break;
        case(kNuDmsq12andAbs23):        return "|#Delta m^{2}_{32}| (eV^{2}c^{-4})";                  break;
        case(kNuDmsq12andAbs31):        return "|#Delta m^{2}_{31}| (eV^{2}c^{-4})";                  break;
        default                :        return "--";                                   break;
      }
    }
    case conv::kPhypCPV     : return "#delta_{CP} (Radians)";
      
    case conv::kPhypMixing12Bar: return "sin^{2}(#bar{#theta_{12}})";
    case conv::kPhypMixing13Bar: return "sin^{2}(#bar{#theta_{13}})";
    case conv::kPhypMixing23Bar: return "sin^{2}(#bar{#theta_{23}})";
    case conv::kPhypDm12sqBar  : return "|#Delta #bar{m^2_{12}}| (eV^{2}c^{-4})";
    case conv::kPhypDm23sqBar  :
    {
      NuOscDmsqScheme_t omp = Opt::Instance()->GetEffectiveDmsqScheme();
      switch(omp) {
        case(kNuDmsqFogliLisi) :        return "|#Delta #bar{m^{2}}| = |m^{2}_{3} - (m^{2}_{1} + m^{2}_{1})/2| (eV^{2}c^{-4})";                  break;
        case(kNuDmsq12andAbs23):        return "|#Delta #bar{m^{2}_{32}}| (eV^{2}c^{-4})";                  break;
        case(kNuDmsq12andAbs31):        return "|#Delta #bar{m^{2}_{31}}| (eV^{2}c^{-4})";                  break;
        default                :        return "--";                                   break;
      }
    }
    case conv::kPhypCPVBar     : return "#bar{#delta_{CP}} (Radians)";

    default: return "AsStringOfficialAxisError";
  }
  return "AsStringOfficialAxisError";
}
