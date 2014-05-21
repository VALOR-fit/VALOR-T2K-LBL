//________________________________________________________________________________________
/*

Copyright (c) 2010-2014, The VALOR group. All rights reserved.

ChangeLog:

*/
//________________________________________________________________________________________

#include <string>
#include <algorithm>
#include <iostream>

#include <TROOT.h>

#include "T2KPhysParamId.h"
#include "Utilities.h"
//~ #include "NuOscMixParam.h"
#include "Enumerations.h"

using namespace valor;
using namespace valor::conv;

//________________________________________________________________________________________
string valor::util::AsString(PhysParamId_t p)
{
  switch(p) {

    case ( kPhypMixing23 ) :
    {
        NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
        switch(omp) {
           case(kNuMixSinsq2Theta) : return "sin^2(2*theta_23)"; break;
           case(kNuMixSinsqTheta ) : return "sin^2(theta_23)";   break;
           default                 : return "--";                break;
         }
    }
    break;
    case ( kPhypMixing13 ) :
    {
        NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
        switch(omp) {
           case(kNuMixSinsq2Theta) : return "sin^2(2*theta_13)";  break;
           case(kNuMixSinsqTheta ) : return "sin^2(theta_13)";    break;
           default                 : return "--";                 break;
         }
    }
    break;
    case ( kPhypMixing12 ) :
    {
        NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
        switch(omp) {
           case(kNuMixSinsq2Theta) : return "sin^2(2*theta_12)";  break;
           case(kNuMixSinsqTheta ) : return "sin^2(theta_12)";    break;
           default                 : return "--";                 break;
         }
    }
    break;
    case ( kPhypCPV ) : 
    {
        conv::NuOscDeltaCP_t cpv = Opt::Instance()->gOscDeltaCP;
        switch(cpv) {
           case(kNuDeltaCPPhase) : return "delta_cp";      break;
           case(kNuDeltaCPSin  ) : return "sin(delta_cp)"; break;
           default               : return "--";            break;
        }
    }
    break;

    case ( kPhypDm23sq  ) :
    {
        NuOscDmsqScheme_t dmscheme = Opt::Instance()->gOscDmsqScheme;
        switch(dmscheme) {
          case(kNuDmsq12andAbs23) : return "|Dm^2_32|";                  break;
          case(kNuDmsq12andAbs31) : return "|Dm^2_31|";                  break;
          case(kNuDmsqFogliLisi ) : return "|Dm^2_FL|";                   break;
          default                 : return "--";                         break;
        }
    }
    break;

    case ( kPhypDm12sq ) : 
    {
        NuOscDmsqScheme_t dmscheme = Opt::Instance()->gOscDmsqScheme;
        switch(dmscheme) {
          case(kNuDmsq12andAbs23) : return "Dm^2_12";                  break;
          case(kNuDmsq12andAbs31) : return "Dm^2_12";                  break;
          case(kNuDmsqFogliLisi ) : return "dm^2_FL";                  break;
          default                 : return "--";                       break;
        }
    }
    break;

    // esoteric hypotheses
    case ( kPhypAlphaDecay            ) : return "alpha_decay";                break;
    case ( kPhypSinsq2ThetaDecay      ) : return "sin^2(2*theta_decay)";       break;
    case ( kPhypMusqDecoh             ) : return "mu^2_decoh";                 break;
    case ( kPhypSinsq2ThetaDecoh      ) : return "sin^2(2*theta_decoh)";       break;

    // Steriles
    case ( kPhypMixing14 ) :
    {
        NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
        switch(omp) {
           case(kNuMixSinsq2Theta) : return "sin^2(2*theta_14)";  break;
           case(kNuMixSinsqTheta ) : return "sin^2(theta_14)";    break;
           default                 : return "--";                 break;
         }
    }
    break;
    case ( kPhypMixing15 ) :
    {
        NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
        switch(omp) {
           case(kNuMixSinsq2Theta) : return "sin^2(2*theta_15)";  break;
           case(kNuMixSinsqTheta ) : return "sin^2(theta_15)";    break;
           default                 : return "--";                 break;
         }
    }
    break;
    case ( kPhypMixing24 ) :
    {
        NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
        switch(omp) {
           case(kNuMixSinsq2Theta) : return "sin^2(2*theta_24)";  break;
           case(kNuMixSinsqTheta ) : return "sin^2(theta_24)";    break;
           default                 : return "--";                 break;
         }
    }
    break;
    case ( kPhypMixing25 ) :
    {
        NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
        switch(omp) {
           case(kNuMixSinsq2Theta) : return "sin^2(2*theta_25)";  break;
           case(kNuMixSinsqTheta ) : return "sin^2(theta_25)";    break;
           default                 : return "--";                 break;
         }
    }
    break;
    case ( kPhypMixing34 ) :
    {
        NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
        switch(omp) {
           case(kNuMixSinsq2Theta) : return "sin^2(2*theta_34)";  break;
           case(kNuMixSinsqTheta ) : return "sin^2(theta_34)";    break;
           default                 : return "--";                 break;
         }
    }
    break;
    case ( kPhypMixing35 ) :
    {
        NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
        switch(omp) {
           case(kNuMixSinsq2Theta) : return "sin^2(2*theta_35)";  break;
           case(kNuMixSinsqTheta ) : return "sin^2(theta_35)";    break;
           default                 : return "--";                 break;
         }
    }
    break;

    case ( kPhypDelta12               ) : return "delta_12";        break;
    case ( kPhypDelta13               ) : return "delta_13";        break;
    case ( kPhypDelta15               ) : return "delta_15";        break;
    case ( kPhypDelta34               ) : return "delta_34";        break;
    case ( kPhypDelta35               ) : return "delta_35";        break;

    case ( kPhypDm41sq                ) : return "|Dm^2_41|";       break;
    case ( kPhypDm51sq                ) : return "|Dm^2_51|";       break;

    // Nonstandard matter interaction
    case ( kPhypNSIepsMuTau           ) : return "nsi_epsilonMuTau";    break;

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
#if (BANFF_BINS != 2)
    case ( kNuisBANFF16 ) :
    case ( kNuisBANFF17 ) :
    case ( kNuisBANFF18 ) :
    case ( kNuisBANFF19 ) :
    case ( kNuisBANFF20 ) :
    case ( kNuisBANFF21 ) :
#endif
#if BANFF_BINS == 1
    case ( kNuisBANFF22 ) :
    case ( kNuisBANFF23 ) :
    case ( kNuisBANFF24 ) :
#endif
    {
      int ip = (int)p - (int)kNuisBANFF00;
      return Form("BANFF %2d; SK %s flux; %.1f - %.1f GeV",
                   ip, util::Pdg2String(kNuPdgCodeOfFluxParamInBANFFfit[ip]).c_str(), 
                   kEtrueRangeOfFluxParamInBANFFfit[ip][0],
                   kEtrueRangeOfFluxParamInBANFFfit[ip][1]);     
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

    //
    // XSec params not correlated with flux
    //

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
    {
      int ip = (int)p - (int)kNuisSKDetFSI000;
      string nutypes = "";
      switch(p) {
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
        default : nutypes = "--";
      }
      return Form("SKDet + FSI/SI %2d; Ereco range %.2f - %.2f GeV; %s",ip,
                  kSKDetFSIParamsErecoCut[ip][0],kSKDetFSIParamsErecoCut[ip][1],nutypes.c_str());     
      break;
    }//FSI for 1Rmu
#endif
#if FSI_BINS != 2
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
      {   
#if FSI_BINS == 3
      int ip = (int)p - (int)kNuisSKDetFSI006;
#elif FSI_BINS ==1
      int ip = (int)p - (int)kNuisSKDetFSI000;
#endif
      string nutypes = "";
        switch(p) {
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
        default : nutypes = "--";
        }
      return Form("SKDet + FSI/SI %2d; Ereco range %.2f - %.2f GeV; %s",ip,
                  kSKDetFSIParamsErecoCut[ip][0],kSKDetFSIParamsErecoCut[ip][1],nutypes.c_str());     
      break;
      }
#endif
#else
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
      int ip = (int)p - (int)kNuisSKDetFSI000;
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
      return Form("SKDet + FSI/SI %2d; Ereco range %.2f - %.2f GeV; %s",ip,
                  kSKDetFSIParamsErecoCut[ip][0],kSKDetFSIParamsErecoCut[ip][1],nutypes.c_str());     
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
      int ip = (int)p - (int)kNuisSKDetFSI000;
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
      return Form("SKDet + FSI/SI %2d; Ereco range %.2f - %.2f GeV; %s",ip,
                  kSKDetFSIParamsErecoCut[ip][0],kSKDetFSIParamsErecoCut[ip][1],nutypes.c_str());     
      break;
      }
#endif
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
#if (BANFF_BINS != 2)
  case ( kRHCNuisBANFF16 ) :
  case ( kRHCNuisBANFF17 ) :
  case ( kRHCNuisBANFF18 ) :
  case ( kRHCNuisBANFF19 ) :
  case ( kRHCNuisBANFF20 ) :
  case ( kRHCNuisBANFF21 ) :
#endif
#if BANFF_BINS == 1
  case ( kRHCNuisBANFF22 ) :
  case ( kRHCNuisBANFF23 ) :
  case ( kRHCNuisBANFF24 ) :
#endif
    {
      int ip = (int)p - (int)kRHCNuisBANFF00;
      return Form("BANFF %2d; SK %s flux; %.1f - %.1f GeV; RHC",
                  ip, util::Pdg2String(kNuPdgCodeOfRHCFluxParamInBANFFfit[ip]).c_str(), 
                  kEtrueRangeOfRHCFluxParamInBANFFfit[ip][0],
                  kEtrueRangeOfRHCFluxParamInBANFFfit[ip][1]);     
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

    //
    // XSec params not correlated with flux
    //

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
    {
      int ip = (int)p - (int)kRHCNuisSKDetFSI000;
      string nutypes = "";
      switch(p) {
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
      default : nutypes = "--";
      }
      return Form("SKDet + FSI/SI %2d; Ereco range %.2f - %.2f GeV; %s; RHC",ip,
                  kRHCSKDetFSIParamsErecoCut[ip][0],kRHCSKDetFSIParamsErecoCut[ip][1],nutypes.c_str());     
      break;
    }//FSI for 1Rmu
#endif
#if FSI_BINS != 2
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
    {
#if FSI_BINS == 3
      int ip = (int)p - (int)kRHCNuisSKDetFSI006;
#elif FSI_BINS ==1
      int ip = (int)p - (int)kRHCNuisSKDetFSI000;
#endif
      string nutypes = "";
      switch(p) {
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
      default : nutypes = "--";
      }
      return Form("SKDet + FSI/SI %2d; Ereco range %.2f - %.2f GeV; %s; RHC",ip,
                  kRHCSKDetFSIParamsErecoCut[ip][0],kRHCSKDetFSIParamsErecoCut[ip][1],nutypes.c_str());     
      break;
    }
#endif
#else
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
      {
        int ip = (int)p - (int)kRHCNuisSKDetFSI000;
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
        return Form("SKDet + FSI/SI %2d; Ereco range %.2f - %.2f GeV; %s; RHC",ip,
                    kRHCSKDetFSIParamsErecoCut[ip][0],kRHCSKDetFSIParamsErecoCut[ip][1],nutypes.c_str());     
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
        int ip = (int)p - (int)kRHCNuisSKDetFSI000;
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
        return Form("SKDet + FSI/SI %2d; Ereco range %.2f - %.2f GeV; %s; RHC",ip,
                    kRHCSKDetFSIParamsErecoCut[ip][0],kRHCSKDetFSIParamsErecoCut[ip][1],nutypes.c_str());     
        break;
      }
#endif

    default : 
      return "--";
  }
  return "--";
}
//________________________________________________________________________________________
string valor::util::AsStringTex(PhysParamId_t p)
{
  // Specify tex strings all ALL subscripts and superscripts in braces.
  // Otherwise conversion to RooTex won't work.
  switch(p) {

    case ( kPhypMixing23 ) :
    {
        NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
        switch(omp) {
           case(kNuMixSinsq2Theta) : return "$\\sin^{2}(2\\theta_{23})$"; break;
           case(kNuMixSinsqTheta ) : return "$\\sin^{2}(\\theta_{23})$";  break;
           default                 : return "--";                       break;
         }
    }
    break;

    case ( kPhypMixing13 ) :
    {
        NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
        switch(omp) {
           case(kNuMixSinsq2Theta) : return "$\\sin^{2}(2\\theta_{13})$"; break;
           case(kNuMixSinsqTheta ) : return "$\\sin^{2}(\\theta_{13})$";  break;
           default                 : return "--";                       break;
         }
    }
    break;

    case ( kPhypMixing12 ) :
    {
        NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
        switch(omp) {
           case(kNuMixSinsq2Theta) : return "$\\sin^{2}(2\\theta_{12})$"; break;
           case(kNuMixSinsqTheta ) : return "$\\sin^{2}(\\theta_{12})$";  break;
           default                 : return "--";                       break;
         }
    }
    break;

    case ( kPhypCPV ) : 
    {
        conv::NuOscDeltaCP_t cpv = Opt::Instance()->gOscDeltaCP;
        switch(cpv) {
           case(kNuDeltaCPPhase) : return "$\\delta_{CP}$";      break;
           case(kNuDeltaCPSin  ) : return "$sin(\\delta_{CP})$"; break;
           default               : return "--";                  break;
        }
    }
    break;

    case ( kPhypDm12sq ) : 
    {
        NuOscDmsqScheme_t dmscheme = Opt::Instance()->gOscDmsqScheme;
        switch(dmscheme) {
          case(kNuDmsq12andAbs23) : return "\\Delta m^2_{22}";                  break;
          case(kNuDmsq12andAbs31) : return "\\Delta m^2_{21}";                  break;
          case(kNuDmsqFogliLisi ) : return "\\delta m^2";                       break;
          default                 : return "--";                                  break;
        }
    }
    break;

    case ( kPhypDm23sq ) :
    {
        NuOscDmsqScheme_t dmscheme = Opt::Instance()->gOscDmsqScheme;
        switch(dmscheme) {
          case(kNuDmsq12andAbs23) : return "|\\Delta m^2_{32}|";                  break;
          case(kNuDmsq12andAbs31) : return "|\\Delta m^2_{31}|";                  break;
          case(kNuDmsqFogliLisi ) : return "|\\Delta m^2|";                       break;
          default                 : return "--";                                  break;
        }
    }
    break;

    case ( kPhypAlphaDecay            ) : return "$\\alpha_{decay}$";                break;
    case ( kPhypSinsq2ThetaDecay      ) : return "$\\sin^{2}(2\\theta_{decay})$";      break;
    case ( kPhypMusqDecoh             ) : return "$\\mu^{2}_{decoh}$";                 break;
    case ( kPhypSinsq2ThetaDecoh      ) : return "$\\sin^{2}(2\\theta_{decoh})$";      break;
    
    // Steriles
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
#if (BANFF_BINS != 2)
    case ( kNuisBANFF16 ) :
    case ( kNuisBANFF17 ) :
    case ( kNuisBANFF18 ) :
    case ( kNuisBANFF19 ) :
    case ( kNuisBANFF20 ) :
    case ( kNuisBANFF21 ) :
#endif
#if BANFF_BINS == 1
    case ( kNuisBANFF22 ) :
    case ( kNuisBANFF23 ) :
    case ( kNuisBANFF24 ) :
#endif
    {
      int ip = (int)p - (int)kNuisBANFF00;
      return Form("BANFF %2d, SK %s flux, %.1f - %.1f GeV",
                   ip, util::Pdg2String(kNuPdgCodeOfFluxParamInBANFFfit[ip]).c_str(), 
                   kEtrueRangeOfFluxParamInBANFFfit[ip][0],
                   kEtrueRangeOfFluxParamInBANFFfit[ip][1]);     
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

    //
    // XSec params not correlated with flux
    //

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
    case ( kNuisXSecNormNuEToNuMu  ) : return "Norm, $\\nu_{e}$ To $\\nu_\\{mu}}$";  break;
#ifndef __NUEONLY__
    case ( kNuisXSecNormAntiNuToNu ) : return "Norm, $\\bar{\\nu}$ To $\\nu$";  break;
#endif
//#ifdef SKFSITOGETHER
#ifndef SKCORRELATIONS
#if FSI_BINS != 3
    case ( kNuisSKDetFSI000 ) :
    case ( kNuisSKDetFSI001 ) :
    case ( kNuisSKDetFSI002 ) :
    case ( kNuisSKDetFSI003 ) :
    case ( kNuisSKDetFSI004 ) :
    case ( kNuisSKDetFSI005 ) :
    {
      string nutypes = "";
      int ip = (int)p - (int)kNuisSKDetFSI000;
      switch(p) {
        case (kNuisSKDetFSI000) :
        case (kNuisSKDetFSI001) :
        case (kNuisSKDetFSI002) :
          nutypes = "$\\nu_\\{mu}}$/$\\bar{\\nu_\\{mu}}}$ CCQE ($1R_{\\mu}$)"; break;
        case (kNuisSKDetFSI003) :
          nutypes = "$\\nu_\\{mu}}$/$\\bar{\\nu_\\{mu}}}$ CCnQE ($1R_{\\mu}$)"; break;
        case (kNuisSKDetFSI004) :
          nutypes = "$\\nu_{e}$/$\\bar{\\nu_{e}}$/sig$\\nu_{e}$ CC ($1R_{\\mu}$)"; break;
        case (kNuisSKDetFSI005) :
          nutypes = "all NC ($1R_{\\mu}$)"; break;
        default : nutypes = "--";
      }
      return Form("SKDet + FSI/SI %2d; Ereco range %.2f - %.2f GeV; %s",ip,
                  kSKDetFSIParamsErecoCut[ip][0],kSKDetFSIParamsErecoCut[ip][1],nutypes.c_str());     
      break;
    }//FSI for 1Rmu
#endif
#if FSI_BINS != 2
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

    {
      string nutypes = "";
#if FSI_BINS != 3
      int ip = (int)p - (int)kNuisSKDetFSI000;
#else
      int ip = (int)p - (int)kNuisSKDetFSI006;
#endif
      switch(p) {
        case (kNuisSKDetFSI006) :
        case (kNuisSKDetFSI007) :
        case (kNuisSKDetFSI008) :
          nutypes = "Osc. $\\nu_{e}$ CC ($1R_{e}$)"; break;
        case (kNuisSKDetFSI009) :
        case (kNuisSKDetFSI010) :
        case (kNuisSKDetFSI011) :
          nutypes = "$\\nu_\\{mu}}$/$\\bar{\\nu_\\{mu}}}$ CC ($1R_{e}$)"; break;
        case (kNuisSKDetFSI012) :
        case (kNuisSKDetFSI013) :
        case (kNuisSKDetFSI014) :
          nutypes = "$\\nu_{e}$/$\\bar{\\nu_{e}}$ CC ($1R_{e}$)"; break;
        case (kNuisSKDetFSI015) :
        case (kNuisSKDetFSI016) :
        case (kNuisSKDetFSI017) :
          nutypes = "all NC ($1R_{e}$)"; break;
        default : nutypes = "--";
      }
      return Form("SKDet + FSI/SI %2d; Ereco range %.2f - %.2f GeV; %s",ip,
                  kSKDetFSIParamsErecoCut[ip][0],kSKDetFSIParamsErecoCut[ip][1],nutypes.c_str());     
      break;
    }//FSI for 1Re
#endif
#else
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
      int ip = (int)p - (int)kNuisSKDetFSI000;
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
          nutypes = "$\\nu_\\{mu}}$/$\\bar{\\nu_\\{mu}}}$ CC ($1R_{\\mu}$)"; break;
        case (kNuisSKDetFSI010) :
        case (kNuisSKDetFSI011) :
        case (kNuisSKDetFSI012) :
        case (kNuisSKDetFSI013) :
        case (kNuisSKDetFSI014) :
          nutypes = "$\\nu_{e}$/$\\bar{\\nu_{e}}$ CC ($1R_{\\mu}$)"; break;
        case (kNuisSKDetFSI015) :
        case (kNuisSKDetFSI016) :
        case (kNuisSKDetFSI017) :
        case (kNuisSKDetFSI018) :
        case (kNuisSKDetFSI019) :
          nutypes = "all NC ($1R_{\\mu}$)"; break;
        default : nutypes = "--";
      }
      return Form("SKDet + FSI/SI %2d; Ereco range %.2f - %.2f GeV; %s",ip,
                  kSKDetFSIParamsErecoCut[ip][0],kSKDetFSIParamsErecoCut[ip][1],nutypes.c_str());     
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
      int ip = (int)p - (int)kNuisSKDetFSI000;
      string nutypes = "";
        switch(p) {
        case (kNuisSKDetFSI020) :
        case (kNuisSKDetFSI021) :
        case (kNuisSKDetFSI022) :
          nutypes = "Osc. $\\nu_{e}$ CC ($1R_{e}$)"; break;
        case (kNuisSKDetFSI023) :
        case (kNuisSKDetFSI024) :
        case (kNuisSKDetFSI025) :
          nutypes = "$\\nu_\\{mu}}$/$\\bar{\\nu_\\{mu}}}$ CC ($1R_{e}$)"; break;
        case (kNuisSKDetFSI026) :
        case (kNuisSKDetFSI027) :
        case (kNuisSKDetFSI028) :
          nutypes = "$\\nu_{e}$/$\\bar{\\nu_{e}}$ CC ($1R_{e}$)"; break;
        case (kNuisSKDetFSI029) :
        case (kNuisSKDetFSI030) :
        case (kNuisSKDetFSI031) :
          nutypes = "all NC ($1R_{e}$)"; break;
        default : nutypes = "--";
        }
      return Form("SKDet + FSI/SI %2d; Ereco range %.2f - %.2f GeV; %s",ip,
                  kSKDetFSIParamsErecoCut[ip][0],kSKDetFSIParamsErecoCut[ip][1],nutypes.c_str());     
      break;
      }
#endif
//#endif //SKDETFSITOGETHER

    case ( kNuisRHCToFHC ) : return "RHC to FHC normalisation"; break;

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
        NuOscDmsqScheme_t dmscheme = Opt::Instance()->gOscDmsqScheme;
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
        NuOscDmsqScheme_t dmscheme = Opt::Instance()->gOscDmsqScheme;
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

    case ( kNuisSuperkEnergyScale     ) : return "f_sk_e_scale";                   break;    
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
#if (BANFF_BINS != 2)
    case ( kNuisBANFF16               ) : return "f_banff16";      break;
    case ( kNuisBANFF17               ) : return "f_banff17";      break;
    case ( kNuisBANFF18               ) : return "f_banff18";      break;
    case ( kNuisBANFF19               ) : return "f_banff19";      break;
    case ( kNuisBANFF20               ) : return "f_banff20";      break;
    case ( kNuisBANFF21               ) : return "f_banff21";      break;
#endif
#if BANFF_BINS == 1
    case ( kNuisBANFF22               ) : return "f_banff22";      break;
    case ( kNuisBANFF23               ) : return "f_banff23";      break;
    case ( kNuisBANFF24               ) : return "f_banff24";      break;
#endif

#ifndef ONEOVERMAQE
    case ( kNuisBANFFMAQE             ) : return "f_banffmaqe";    break;
#else
    case ( kNuisBANFFMAQE             ) : return "f_banff1overmaqe2";  break;
#endif
    case ( kNuisBANFFMARES            ) : return "f_banffmares";   break;
    case ( kNuisBANFFCCQEE0           ) : return "f_banffccqee0";  break;
#if BANFF_BINS != 3
    case ( kNuisBANFFCCQEE1           ) : return "f_banffccqee1";  break;
    case ( kNuisBANFFCCQEE2           ) : return "f_banffccqee2";  break;
#endif
    case ( kNuisBANFFCC1piE0          ) : return "f_banffcc1pie0"; break;
#if BANFF_BINS != 3
    case ( kNuisBANFFCC1piE1          ) : return "f_banffcc1pie1"; break;
#endif
#if (BANFF_BINS != 2)
    case ( kNuisBANFFNC1pi0           ) : return "f_banffnc1pi0";  break;
#endif
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
#endif
#if FSI_BINS != 2
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
#endif
#ifdef SKCORRELATIONS
  case ( kNuisSKDetFSI018 ) : return "f_skdetfsi018"; break;
  case ( kNuisSKDetFSI019 ) : return "f_skdetfsi019"; break;
  case ( kNuisSKDetFSI020 ) : return "f_skdetfsi020"; break;
  case ( kNuisSKDetFSI021 ) : return "f_skdetfsi021"; break;
  case ( kNuisSKDetFSI022 ) : return "f_skdetfsi022"; break;
  case ( kNuisSKDetFSI023 ) : return "f_skdetfsi023"; break;
  case ( kNuisSKDetFSI024 ) : return "f_skdetfsi024"; break;
  case ( kNuisSKDetFSI025 ) : return "f_skdetfsi025"; break;
  case ( kNuisSKDetFSI026 ) : return "f_skdetfsi026"; break;
  case ( kNuisSKDetFSI027 ) : return "f_skdetfsi027"; break;
  case ( kNuisSKDetFSI028 ) : return "f_skdetfsi028"; break;
  case ( kNuisSKDetFSI029 ) : return "f_skdetfsi029"; break;
  case ( kNuisSKDetFSI030 ) : return "f_skdetfsi030"; break;
  case ( kNuisSKDetFSI031 ) : return "f_skdetfsi031"; break;
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
#if (BANFF_BINS != 2)
  case ( kRHCNuisBANFF16               ) : return "f_banff16_rhc";      break;
  case ( kRHCNuisBANFF17               ) : return "f_banff17_rhc";      break;
  case ( kRHCNuisBANFF18               ) : return "f_banff18_rhc";      break;
  case ( kRHCNuisBANFF19               ) : return "f_banff19_rhc";      break;
  case ( kRHCNuisBANFF20               ) : return "f_banff20_rhc";      break;
  case ( kRHCNuisBANFF21               ) : return "f_banff21_rhc";      break;
#endif
#if BANFF_BINS == 1
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
#endif
#if FSI_BINS != 2
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
#endif
#ifdef SKCORRELATIONS
  case ( kRHCNuisSKDetFSI018 ) : return "f_skdetfsi018_rhc"; break;
  case ( kRHCNuisSKDetFSI019 ) : return "f_skdetfsi019_rhc"; break;
  case ( kRHCNuisSKDetFSI020 ) : return "f_skdetfsi020_rhc"; break;
  case ( kRHCNuisSKDetFSI021 ) : return "f_skdetfsi021_rhc"; break;
  case ( kRHCNuisSKDetFSI022 ) : return "f_skdetfsi022_rhc"; break;
  case ( kRHCNuisSKDetFSI023 ) : return "f_skdetfsi023_rhc"; break;
  case ( kRHCNuisSKDetFSI024 ) : return "f_skdetfsi024_rhc"; break;
  case ( kRHCNuisSKDetFSI025 ) : return "f_skdetfsi025_rhc"; break;
  case ( kRHCNuisSKDetFSI026 ) : return "f_skdetfsi026_rhc"; break;
  case ( kRHCNuisSKDetFSI027 ) : return "f_skdetfsi027_rhc"; break;
  case ( kRHCNuisSKDetFSI028 ) : return "f_skdetfsi028_rhc"; break;
  case ( kRHCNuisSKDetFSI029 ) : return "f_skdetfsi029_rhc"; break;
  case ( kRHCNuisSKDetFSI030 ) : return "f_skdetfsi030_rhc"; break;
  case ( kRHCNuisSKDetFSI031 ) : return "f_skdetfsi031_rhc"; break;
#endif
//#endif //SKDETFSITOGETHER

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
      NuOscDmsqScheme_t omp = Opt::Instance()->gOscDmsqScheme;
      switch(omp) {
        case(kNuDmsqFogliLisi) :        return "$|\\Deltam^{2}| = |m^{2}_{3} - (m^{2}_{1} + m^{2}_{1})/2|$";                  break;
        case(kNuDmsq12andAbs23):        return "$|\\Deltam^{2}_{32}|$";                  break;
        case(kNuDmsq12andAbs31):        return "$|\\Deltam^{2}_{31}|$";                  break;
        default                :        return "--";                                   break;
      }
    }
    break;

  case ( kPhypDm12sq                ) : 
    {
      NuOscDmsqScheme_t omp = Opt::Instance()->gOscDmsqScheme;
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
#if (BANFF_BINS != 2)
  case ( kNuisBANFF16 ) :
  case ( kNuisBANFF17 ) :
  case ( kNuisBANFF18 ) :
  case ( kNuisBANFF19 ) :
  case ( kNuisBANFF20 ) :
  case ( kNuisBANFF21 ) :
#endif
#if BANFF_BINS == 1
  case ( kNuisBANFF22 ) :
  case ( kNuisBANFF23 ) :
  case ( kNuisBANFF24 ) :
#endif
  case ( kNuisBANFFMAQE ) :
  case ( kNuisBANFFMARES ) :
  case ( kNuisBANFFCCQEE0 ) :
#if BANFF_BINS != 3
  case ( kNuisBANFFCCQEE1 ) :
  case ( kNuisBANFFCCQEE2 ) :
#endif
  case ( kNuisBANFFCC1piE0 ) :
#if BANFF_BINS != 3
  case ( kNuisBANFFCC1piE1 ) :
#endif
#if (BANFF_BINS != 2)
  case ( kNuisBANFFNC1pi0 ) :
#endif
    {
      int ip = (int)p - (int)kNuisBANFF00;
      return Form("$f_{%d; t, r}^{banff}$", ip);
      break;
    }
    
    //
    // XSec params not correlated with flux 
    // 

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
  case ( kNuisXSecNormNuEToNuMu  ) : return "$f_{CC\\nu_{e}; t}$";           break;
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
#endif
#if FSI_BINS != 2
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
#endif
#ifdef SKCORRELATIONS
  case ( kNuisSKDetFSI018 ) :
  case ( kNuisSKDetFSI019 ) :
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
#endif

    {
#if FSI_BINS != 3
      int ip = (int)p - (int)kNuisSKDetFSI000;
#else
      int ip = (int)p - (int)kNuisSKDetFSI006;
#endif
      return Form("$f_{%d; t, r}^{SK+FSI}$", ip);
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
#if (BANFF_BINS != 2)
  case ( kRHCNuisBANFF16 ) :
  case ( kRHCNuisBANFF17 ) :
  case ( kRHCNuisBANFF18 ) :
  case ( kRHCNuisBANFF19 ) :
  case ( kRHCNuisBANFF20 ) :
  case ( kRHCNuisBANFF21 ) :
#endif
#if BANFF_BINS == 1
  case ( kRHCNuisBANFF22 ) :
  case ( kRHCNuisBANFF23 ) :
  case ( kRHCNuisBANFF24 ) :
#endif
  case ( kRHCNuisBANFFMAQE ) :
  case ( kRHCNuisBANFFMARES ) :
  case ( kRHCNuisBANFFCCQEE0 ) :
#if BANFF_BINS != 3
  case ( kRHCNuisBANFFCCQEE1 ) :
  case ( kRHCNuisBANFFCCQEE2 ) :
#endif
  case ( kRHCNuisBANFFCC1piE0 ) :
#if BANFF_BINS != 3
  case ( kRHCNuisBANFFCC1piE1 ) :
#endif
#if (BANFF_BINS != 2)
  case ( kRHCNuisBANFFNC1pi0 ) :
#endif
    {
      int ip = (int)p - (int)kRHCNuisBANFF00;
      return Form("$f_{%d; t, r}^{banff}RHC$", ip);
      break;
    }
    
    //
    // XSec params not correlated with flux 
    // 

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
#endif
#if FSI_BINS != 2
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
#endif
#ifdef SKCORRELATIONS
  case ( kRHCNuisSKDetFSI018 ) :
  case ( kRHCNuisSKDetFSI019 ) :
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
#endif

    {
#if FSI_BINS != 3
      int ip = (int)p - (int)kRHCNuisSKDetFSI000;
#else
      int ip = (int)p - (int)kRHCNuisSKDetFSI006;
#endif
      return Form("$f_{%d; t, r}^{SK+FSI}RHC$", ip);
      break;
    }//FSI

  default :
    return "--";
  }
  return "--";
}
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
#if (BANFF_BINS != 2)
  case ( kRHCNuisBANFF16               ) : return "f_banff16";      break;
  case ( kRHCNuisBANFF17               ) : return "f_banff17";      break;
  case ( kRHCNuisBANFF18               ) : return "f_banff18";      break;
  case ( kRHCNuisBANFF19               ) : return "f_banff19";      break;
  case ( kRHCNuisBANFF20               ) : return "f_banff20";      break;
  case ( kRHCNuisBANFF21               ) : return "f_banff21";      break;
#endif
#if BANFF_BINS == 1
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
//#ifdef SKFSITOGETHER
#if FSI_BINS != 3
  case ( kRHCNuisSKDetFSI000 ) : return "f_skdetfsi000"; break;
  case ( kRHCNuisSKDetFSI001 ) : return "f_skdetfsi001"; break;
  case ( kRHCNuisSKDetFSI002 ) : return "f_skdetfsi002"; break;
  case ( kRHCNuisSKDetFSI003 ) : return "f_skdetfsi003"; break;
  case ( kRHCNuisSKDetFSI004 ) : return "f_skdetfsi004"; break;
  case ( kRHCNuisSKDetFSI005 ) : return "f_skdetfsi005"; break;
#endif
#if FSI_BINS != 2
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
#endif
#ifdef SKCORRELATIONS
  case ( kRHCNuisSKDetFSI018 ) : return "f_skdetfsi018"; break;
  case ( kRHCNuisSKDetFSI019 ) : return "f_skdetfsi019"; break;
  case ( kRHCNuisSKDetFSI020 ) : return "f_skdetfsi020"; break;
  case ( kRHCNuisSKDetFSI021 ) : return "f_skdetfsi021"; break;
  case ( kRHCNuisSKDetFSI022 ) : return "f_skdetfsi022"; break;
  case ( kRHCNuisSKDetFSI023 ) : return "f_skdetfsi023"; break;
  case ( kRHCNuisSKDetFSI024 ) : return "f_skdetfsi024"; break;
  case ( kRHCNuisSKDetFSI025 ) : return "f_skdetfsi025"; break;
  case ( kRHCNuisSKDetFSI026 ) : return "f_skdetfsi026"; break;
  case ( kRHCNuisSKDetFSI027 ) : return "f_skdetfsi027"; break;
  case ( kRHCNuisSKDetFSI028 ) : return "f_skdetfsi028"; break;
  case ( kRHCNuisSKDetFSI029 ) : return "f_skdetfsi029"; break;
  case ( kRHCNuisSKDetFSI030 ) : return "f_skdetfsi030"; break;
  case ( kRHCNuisSKDetFSI031 ) : return "f_skdetfsi031"; break;
#endif
//#endif //SKDETFSITOGETHER

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
  return !IsRHCSystematic(p);
}
//________________________________________________________________________________________
bool valor::util::IsRHCSystematic(PhysParamId_t p)
{
  if(!IsSystematic(p))
    return false;
  if(p < kNuisRHCToFHC)
    return false;
  else if(p < kRHCLastParam)
    return true;
  std::cerr<<"PhysParamId/IsRHCSystematic >> Cannot determine whether or not RHC syst: "<<AsString(p)
    <<". Assuming no."<<std::endl;
  return false;
}
//________________________________________________________________________________________
//____________________________________________________________________________
bool valor::util::IsXSec(PhysParamId_t p)
{
  if(p >= kNuisXSecFermiMom && p <= kNuisXSecNormAntiNuToNu)
    return true;
  if(p >= kRHCNuisXSecFermiMom && p <= kRHCNuisXSecNormAntiNuToNu)
    return true;
  return false;
  
}
//____________________________________________________________________________
bool valor::util::IsUsingSplines(PhysParamId_t p)
{
  assert(!"IsUsingSplines Currently unimplemented");
  if(p >= kNuisXSecFermiMom && p <= kNuisXSecBindE)
    return true;
  if(p >= kRHCNuisXSecFermiMom && p <= kRHCNuisXSecBindE)
    return true;
  if(p == kNuisBANFFMAQE || p == kRHCNuisBANFFMAQE)
    return true;

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
  if(p >= kNuisSKDetFSI000 && p <= kNuisSKDetFSI031)
    return true;
#else
#if FSI_BINS != 3
  if(p >= kNuisSKDetFSI000 && p <= kNuisSKDetFSI005)
    return true;
#endif
#if FSI_BINS != 2
  if(p >= kNuisSKDetFSI006 && p <= kNuisSKDetFSI017)
    return true;
#endif
#endif
  if(p == kNuisSuperkEnergyScale)
    return true;
    
  //rhc
#ifdef SKCORRELATIONS
  if(p >= kRHCNuisSKDetFSI000 && p <= kRHCNuisSKDetFSI031)
    return true;
#else
#if FSI_BINS != 3
  if(p >= kRHCNuisSKDetFSI000 && p <= kRHCNuisSKDetFSI005)
    return true;
#endif
#if FSI_BINS != 2
  if(p >= kRHCNuisSKDetFSI006 && p <= kRHCNuisSKDetFSI017)
    return true;
#endif
#endif
  if(p == kRHCNuisSuperkEnergyScale)
    return true;

  return false;
}
//____________________________________________________________________________   
bool valor::util::IsBANFF(PhysParamId_t p)
{
  //fhc
  if(p >= kNuisBANFF00
#if (BANFF_BINS != 2)
     && p <= kNuisBANFFNC1pi0
#else
     && p <= kNuisBANFFCC1piE1
#endif
     )
    return true;
    
  //rhc
  if(p >= kRHCNuisBANFF00
#if (BANFF_BINS != 2)
     && p <= kRHCNuisBANFFNC1pi0
#else
     && p <= kRHCNuisBANFFCC1piE1
#endif
     )
    return true;
    
  return false;
}
//____________________________________________________________________________

string valor::util::AsStringCompact_KDI(PhysParamId_t p)
{
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
valor::conv::PhysParamId_t
  valor::util::PhysParamIdFromString(string sh) 
{
  string lsh = util::ToLowerCase(sh);
  for (int i = 0 ; i < knppmax ; i++) {
    PhysParamId_t k = (PhysParamId_t) i ;
    string s = util::ToLowerCase(AsStringCompact(k));
    string s2;
    if (s2.find("f_") == 0) s2 = s.substr(2);
    else s2 = "gskfgjalkfjgalfkgjadlfj";
    if (k != kPhypUndefined && (lsh == s || lsh == s2) ) return k;
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
string valor::util::ParamAsString                   (conv::PhysParamId_t p) { return util::AsString                  (p);}
string valor::util::ParamAsStringTex                (conv::PhysParamId_t p) { return util::AsStringTex               (p);}
string valor::util::ParamAsStringPaperTexNoDollar   (conv::PhysParamId_t p) { return util::AsStringPaperTexNoDollar  (p);}
string valor::util::ParamAsStringPaperTex           (conv::PhysParamId_t p) { return util::AsStringPaperTex          (p);}
string valor::util::ParamAsStringPaperRooTex        (conv::PhysParamId_t p) { return util::AsStringPaperRooTex       (p);}
string valor::util::ParamAsStringCompact            (conv::PhysParamId_t p) { return util::AsStringCompact           (p);}
string valor::util::ParamAsStringCompact_splinehack (conv::PhysParamId_t p) { return util::AsStringCompact_splinehack(p);}
string valor::util::ParamAsStringCompact_KDI        (conv::PhysParamId_t p) { return util::AsStringCompact_KDI       (p);}
//____________________________________________________________________________
