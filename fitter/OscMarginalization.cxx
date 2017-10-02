#include "OscMarginalization.h"

#include "Binning.h"
#include "SuperkModes.h"
#include "SuperkSamples.h"
#include "SamplesAndModes.h"

//#include "valor/Conventions/NuOscMixParam.h"
#include "valor/Conventions/Opt.h"
#include "valor/Utilities/UtlNumStat.h"
#include "valor/Utilities/Rnd.h"

#include <cstdlib>
#include <iostream>

#include <TStopwatch.h>

#define __USE_PRIOR_SINSQ2TH13__    // Prior specified in terms of sinsq2th
#define __USE_PRIOR_SINSQ2TH12__    // Prior specified in terms of sinsq2th

namespace valor
{
namespace util
{

void SetOscillationParameters(NuOsc3Flavour& osc, FitParamList* fparm) 
{
   using namespace valor::conv;
   double mixing12 = fparm->CurrValue(conv::kPhypMixing12);
   double mixing13 = fparm->CurrValue(conv::kPhypMixing13);
   double mixing23 = fparm->CurrValue(conv::kPhypMixing23);
   double dm12sq = fparm->CurrValue(conv::kPhypDm12sq);
   double dm23sq = fparm->CurrValue(conv::kPhypDm23sq);
   double deltacp = fparm->CurrValue(conv::kPhypCPV);
    
   NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
   NuOscDmsqScheme_t scheme = Opt::Instance()->GetEffectiveDmsqScheme();
    
   if(scheme == kNuDmsq12andAbs23 && omp == kNuMixSinsq2Theta)
   {
      osc.SetOscParameters(mixing12, mixing13, mixing23, dm12sq, dm23sq, deltacp);
   }
   else if(scheme == kNuDmsq12andAbs23 && omp == kNuMixSinsqTheta)
   {
      osc.SetOscParametersSinSqTAll(mixing12, mixing13, mixing23, dm12sq, dm23sq, deltacp);
   }
   else if(scheme == kNuDmsqFogliLisi && omp == kNuMixSinsq2Theta)
   {
      osc.SetOscParametersSSq2TDmAtmSol(mixing12, mixing13, mixing23, dm12sq, dm23sq, deltacp);
   }
   else if(scheme == kNuDmsqFogliLisi && omp == kNuMixSinsqTheta)
   {
      osc.SetOscParametersDmsqAtmSol(mixing12, mixing13, mixing23, dm12sq, dm23sq, deltacp);
   } 
   else if(scheme == kNuDmsq12andAbs31 && omp == kNuMixSinsq2Theta)
   {
      osc.SetOscParametersSSq2TDmsq13(mixing12, mixing13, mixing23, dm12sq, dm23sq, deltacp);
   }
   else if(scheme == kNuDmsq12andAbs31 && omp == kNuMixSinsqTheta)
   {
      osc.SetOscParametersDmsq13(mixing12, mixing13, mixing23, dm12sq, dm23sq, deltacp);
   }
   else
   {
      std::cerr << "OscMarginalization: Can't handle choice of oscillation mixing params: " <<
           util::AsString(omp) << "\n and independent squared mass-splittings: " <<
           util::AsString(scheme) << std::endl;
      exit(1);
   }
}

void SetOscillationBarParameters(NuOsc3Flavour& osc, FitParamList* fparm) 
{
   using namespace valor::conv;
   double mixing12bar = fparm->CurrValue(conv::kPhypMixing12Bar);
   double mixing13bar = fparm->CurrValue(conv::kPhypMixing13Bar);
   double mixing23bar = fparm->CurrValue(conv::kPhypMixing23Bar);
   double dm12sqbar = fparm->CurrValue(conv::kPhypDm12sqBar);
   double dm23sqbar = fparm->CurrValue(conv::kPhypDm23sqBar);
   double deltacpbar = fparm->CurrValue(conv::kPhypCPVBar);

   NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
   NuOscDmsqScheme_t scheme = Opt::Instance()->GetEffectiveDmsqScheme();
    
   if(scheme == kNuDmsq12andAbs23 && omp == kNuMixSinsq2Theta)
   {
      osc.SetOscParameters(mixing12bar, mixing13bar, mixing23bar, dm12sqbar, dm23sqbar, deltacpbar);
   }
   else if(scheme == kNuDmsq12andAbs23 && omp == kNuMixSinsqTheta)
   {
      osc.SetOscParametersSinSqTAll(mixing12bar, mixing13bar, mixing23bar, dm12sqbar, dm23sqbar, deltacpbar);
   }
   else if(scheme == kNuDmsqFogliLisi && omp == kNuMixSinsq2Theta)
   {
      osc.SetOscParametersSSq2TDmAtmSol(mixing12bar, mixing13bar, mixing23bar, dm12sqbar, dm23sqbar, deltacpbar);
   }
   else if(scheme == kNuDmsqFogliLisi && omp == kNuMixSinsqTheta)
   {
      osc.SetOscParametersDmsqAtmSol(mixing12bar, mixing13bar, mixing23bar, dm12sqbar, dm23sqbar, deltacpbar);
   } 
   else if(scheme == kNuDmsq12andAbs31 && omp == kNuMixSinsq2Theta)
   {
      osc.SetOscParametersSSq2TDmsq13(mixing12bar, mixing13bar, mixing23bar, dm12sqbar, dm23sqbar, deltacpbar);
   }
   else if(scheme == kNuDmsq12andAbs31 && omp == kNuMixSinsqTheta)
   {
      osc.SetOscParametersDmsq13(mixing12bar, mixing13bar, mixing23bar, dm12sqbar, dm23sqbar, deltacpbar);
   }
   else
   {
      std::cerr << "OscMarginalization: Can't handle choice of oscillation mixing params: " <<
           util::AsString(omp) << "\n and independent squared mass-splittings: " <<
           util::AsString(scheme) << std::endl;
      exit(1);
   }
}

void Randomize(FitParamList* fparm)
{
   using namespace valor::conv;
   const int N = fparm->Size();
   for(int i = 0; i < N; ++i)
   {
      PhysParamId_t p = fparm->WhichParam(i);
      if(!fparm->IsMarginalized(p))
      {
         continue;
      }
      double val = conv::kBadResult;
      bool invalid = true;
      do
      {
         if(fparm->GetPriorModeName(p) == "flat")
         {
            double min = fparm->GetMinValueFlat(p);
            double max = fparm->GetMaxValueFlat(p);
            val = Rnd::Instance()->RndToyMCosc().Uniform(min, max);
         }
         else if(fparm->GetPriorModeName(p) == "gaus")
         {
            double mean = fparm->GetMeanGaus(p);
            double sigma = fparm->GetSigmaGaus(p);
            val = Rnd::Instance()->RndToyMCosc().Gaus(mean, sigma);
         }
         else
         {
            std::cerr << "Randomize" << std::endl;
            std::cerr << "Used prior is: " << fparm->GetPriorModeName(p) << std::endl;
            std::cerr << "Prior not well defined for " << util::AsString(p) << std::endl;
            exit(1);
         }
        
         // Convert from sinsq2th to sinsqth or viceversa if mixing angle
         if(util::IsMixingAngle(p))
         {
            if(Opt::Instance()->gOscMixingParamChoice == kNuMixSinsqTheta)
            {
               #ifdef __USE_PRIOR_SINSQ2TH23__
               if(p == conv::kPhypMixing23 || p == conv::kPhypMixing23Bar)
               {
                  val = FromSinsq2thToSinsqth(val);
               }
               #endif
               #ifdef __USE_PRIOR_SINSQ2TH12__
               if(p == conv::kPhypMixing12 || p == conv::kPhypMixing12Bar)
               {
                  val = FromSinsq2thToSinsqth(val);
               }
               #endif
               #ifdef __USE_PRIOR_SINSQ2TH13__
               if(p == conv::kPhypMixing13 || p == conv::kPhypMixing13Bar)
               {
                  val = FromSinsq2thToSinsqth(val);
               }
               #endif
            }
            else if(Opt::Instance()->gOscMixingParamChoice == kNuMixSinsq2Theta)
            {
               #ifndef __USE_PRIOR_SINSQ2TH23__
               if(p == conv::kPhypMixing23 || p == conv::kPhypMixing23Bar)
               {
                  val = FromSinsqthToSinsq2th(val);
               }
               #endif
               #ifndef __USE_PRIOR_SINSQ2TH12__
               if(p == conv::kPhypMixing12 || p == conv::kPhypMixing12Bar)
               {
                  val = FromSinsqthToSinsq2th(val);
               }
               #endif
               #ifndef __USE_PRIOR_SINSQ2TH13__
               if(p == conv::kPhypMixing13 || p == conv::kPhypMixing13Bar)
               {
                  val = FromSinsqthToSinsq2th(val);
               }
               #endif
            }
         }

         if(val >= fparm->MinValue(p) && val <= fparm->MaxValue(p))
         {
            invalid = false;
         }
      }
      while(invalid);
      fparm->Update(p, val);
   }
}

void CalcOscProbabilities(const int Ntoys, FitParamList* fparm)
{
   using namespace valor::conv;
   double p_osc[3][3];
   double p_osc_bar[3][3];
   NuOsc3Flavour osc;

   TStopwatch watch;

   for(int n = 0; n < Ntoys; ++n)
   {
      Randomize(fparm);
      SetOscillationParameters(osc, fparm);

      for(int is = 0; is < kns; ++is)
      {
         SuperkSample_t s = kSKSamples[is];
         if(s == conv::kSKsUndefined)
         {
            continue;
         }
           
         int n_etrue_bins = kNEtrueBins[is];
         for(int it = 1; it <= n_etrue_bins; ++it)
         {
            double Etrue = (util::EtrueBinEdges(s, it - 1) + util::EtrueBinEdges(s, it)) * 0.5;
            osc.CalcMatterOscProbAll(false, Etrue, p_osc);
            osc.CalcMatterOscProbAll(true, Etrue, p_osc_bar);
               
            for(int im=0; im < knm; ++im)
            {
               SuperkMode_t m = kSKReactionModes[is][im];
               if(m == conv::kSKmUndefined)
               {
                  continue;
               }
               double posc = 1.;
               if(kIsCC[is][im])
               {
                  if(util::IsNuMu(m))
                  {
                     posc = p_osc[1][1];
                  }
                  else if(util::IsNuMuBar(m))
                  {
                     posc = p_osc_bar[1][1];
                  }
                  else if(util::IsNuE(m))
                  {
                     posc = p_osc[0][0];
                  }
                  else if(util::IsNuEBar(m))
                  {
                     posc = p_osc_bar[0][0];
                  }
                  else if(util::IsOscNuMu(m))
                  {
                     posc = p_osc[0][1];
                  }
                  else if(util::IsOscNuMuBar(m))
                  {
                     posc = p_osc_bar[0][1];
                  }
                  else if(util::IsOscNuE(m))
                  {
                     posc = p_osc[1][0];
                  }
                  else if(util::IsOscNuEBar(m))
                  {
                     posc = p_osc_bar[1][0];
                  }
               }
            }
         }
      }
      Rnd::Instance()->SetSeeds(
            Rnd::Instance()->RndToyMCstat().GetSeed() + 1,
            Rnd::Instance()->RndToyMCsyst().GetSeed() + 1,
            Rnd::Instance()->RndToyMCosc().GetSeed() + 1,
            Rnd::Instance()->RndFitter().GetSeed() + 1,
            Rnd::Instance()->RndOther().GetSeed()+ 1);
   }

   std::cout << "Real time: " << watch.RealTime() << std::endl;
}


void CalcOscBarProbabilities(const int Ntoys, FitParamList* fparm)
{
   using namespace valor::conv;
   double p_osc[3][3];
   double p_osc_bar[3][3];
   NuOsc3Flavour osc;
   NuOsc3Flavour oscbar;

   TStopwatch watch;

   for(int n = 0; n < Ntoys; ++n)
   {
      Randomize(fparm);
      SetOscillationParameters(osc, fparm);
      SetOscillationBarParameters(oscbar, fparm);
      for(int is = 0; is < kns; ++is)
      {
         SuperkSample_t s = kSKSamples[is];
         if(s == conv::kSKsUndefined)
         {
            continue;
         }
           
         int n_etrue_bins = kNEtrueBins[is];
         for(int it = 1; it <= n_etrue_bins; ++it)
         {
            double Etrue = (util::EtrueBinEdges(s, it - 1) + util::EtrueBinEdges(s, it)) * 0.5;
            osc.CalcMatterOscProbAll(false, Etrue, p_osc);
            oscbar.CalcMatterOscProbAll(true, Etrue, p_osc_bar);
               
            for(int im=0; im < knm; ++im)
            {
               SuperkMode_t m = kSKReactionModes[is][im];
               if(m == conv::kSKmUndefined)
               {
                  continue;
               }
               double posc = 1.;
               if(kIsCC[is][im])
               {
                  if(util::IsNuMu(m))
                  {
                     posc = p_osc[1][1];
                  }
                  else if(util::IsNuMuBar(m))
                  {
                     posc = p_osc_bar[1][1];
                  }
                  else if(util::IsNuE(m))
                  {
                     posc = p_osc[0][0];
                  }
                  else if(util::IsNuEBar(m))
                  {
                     posc = p_osc_bar[0][0];
                  }
                  else if(util::IsOscNuMu(m))
                  {
                     posc = p_osc[0][1];
                  }
                  else if(util::IsOscNuMuBar(m))
                  {
                     posc = p_osc_bar[0][1];
                  }
                  else if(util::IsOscNuE(m))
                  {
                     posc = p_osc[1][0];
                  }
                  else if(util::IsOscNuEBar(m))
                  {
                     posc = p_osc_bar[1][0];
                  }
               }
            }
         }
      }
      Rnd::Instance()->SetSeeds(
            Rnd::Instance()->RndToyMCstat().GetSeed() + 1,
            Rnd::Instance()->RndToyMCsyst().GetSeed() + 1,
            Rnd::Instance()->RndToyMCosc().GetSeed() + 1,
            Rnd::Instance()->RndFitter().GetSeed() + 1,
            Rnd::Instance()->RndOther().GetSeed()+ 1);
   }

   std::cout << "Real time: " << watch.RealTime() << std::endl;
}

}
}

