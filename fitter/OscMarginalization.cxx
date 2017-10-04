#include "OscMarginalization.h"

#include "Binning.h"
#include "Enumerations.h"
#include "SamplesAndModes.h"

//#include "valor/Conventions/NuOscMixParam.h"
#include "valor/Conventions/Opt.h"
#include "valor/Utilities/UtlNumStat.h"
#include "valor/Utilities/Rnd.h"

#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <map>

#include <TStopwatch.h>

#define __USE_PRIOR_SINSQ2TH13__    // Prior specified in terms of sinsq2th
#define __USE_PRIOR_SINSQ2TH12__    // Prior specified in terms of sinsq2th

namespace valor
{

OscProb *OscProb::instance = 0;

OscProb::OscProb(const int Ntoys) :
   N(Ntoys),
   vProb_1Re(0),
   vProb_1Rmu(0),
   vProb_RHC1Re(0),
   vProb_RHC1Rmu(0),
   vProb_MultiRe(0)
{
   InitSampleMap();
   InitModeMaps();
   vProb_1Re = new VALORArrayF(N, util::NumEtrueBins(conv::kSKsFCFV1RingELike), NFlavours);
   vProb_1Rmu = new VALORArrayF(N, util::NumEtrueBins(conv::kSKsFCFV1RingMuLike), NFlavours);
   vProb_RHC1Re = new VALORArrayF(N, util::NumEtrueBins(conv::kSKsRHCFCFV1RingELike), NFlavours);
   vProb_RHC1Rmu = new VALORArrayF(N, util::NumEtrueBins(conv::kSKsRHCFCFV1RingMuLike), NFlavours);
   vProb_MultiRe = new VALORArrayF(N, util::NumEtrueBins(conv::kSKsFCFVMultiRingELike), NFlavours);
}

OscProb::~OscProb()
{
   delete vProb_1Re;
   delete vProb_1Rmu;
   delete vProb_RHC1Re;
   delete vProb_RHC1Rmu;
   delete vProb_MultiRe;
}

OscProb* OscProb::Instance(const int Ntoys)
{
   if(instance == 0)
   {
      instance = new OscProb(Ntoys);
   }

   return instance;
}

OscProb* OscProb::Instance()
{
   if(instance == 0)
   {
      std::cerr << "Error: Attempt to access OscProb prior to initialisation" <<
         " with Instance(const int Ntoys)." << std::endl;
      exit(1);
   }

   return instance;
}

VALORArrayF* OscProb::GetToys(conv::SuperkSample_t s)
{
   switch(s)
   {
      case conv::kSKsFCFV1RingELike:
         return vProb_1Re;
      case conv::kSKsFCFV1RingMuLike:
         return vProb_1Rmu;
      case conv::kSKsRHCFCFV1RingELike:
         return vProb_RHC1Re;
      case conv::kSKsRHCFCFV1RingMuLike:
         return vProb_RHC1Rmu;
      case conv::kSKsFCFVMultiRingELike:
         return vProb_MultiRe;
      default:
         return 0;
   }
}

void OscProb::Write()
{
   for(int is = 0; is < conv::kNSuperkSamples; ++is)
   {
      conv::SuperkSample_t s = conv::kSKSamples[is];
      if(s == conv::kSKsUndefined)
      {
         continue;
      }
      VALORArrayF* vProb = GetToys(s);
      int ndim = vProb->GetNDim();
      int size = vProb->GetSize();
      std::ofstream os_toys(Form("osc_toys_%s.bin", util::AsStringCompact(s).c_str()),
            std::ios::binary | std::ios::out);
      os_toys.write((char*)&ndim, sizeof(int));
      os_toys.write((char*)vProb->GetDimensions(), 6 * sizeof(int));
      os_toys.write((char*)vProb->GetArray(), size * sizeof(int));
      os_toys.close();
   }
}

void OscProb::InitSampleMap()
{
   for(int is = 0; is < conv::kNSuperkSamples; ++is)
   {
      conv::SuperkSample_t s = conv::kSKSamples[is];
      if(s == conv::kSKsUndefined)
      {
         continue;
      }
      sample_map[s] = is;
   }
}

void OscProb::InitModeMaps()
{
   for(int is = 0; is < conv::kNSuperkSamples; ++is)
   {
      conv::SuperkSample_t s = conv::kSKSamples[is];
      if(s == conv::kSKsUndefined)
      {
         continue;
      }
      int idx = 0;
      for(int im = 0; im < conv::knm; ++im)
      {
         conv::SuperkMode_t m = conv::kSKReactionModes[is][im];
         if(m == conv::kSKmUndefined)
         {
            continue;
         }
         if(conv::kIsCC[is][im])
         {
            mode_maps[sample_map[s]][m] = idx;
            ++idx;
         }
      }
   }
}

int OscProb::ModeIndex(conv::SuperkSample_t s, conv::SuperkMode_t m)
{
   return mode_maps[sample_map[s]][m];
}

int OscProb::SampleIndex(conv::SuperkSample_t s)
{
   return sample_map[s];
}

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

void CalcOscProbabilities(const int Ntoys, FitParamList* fparm, bool isbar)
{
   using namespace valor::conv;
   double p_osc[3][3];
   double p_osc_bar[3][3];
   NuOsc3Flavour osc;
   NuOsc3Flavour oscbar;

   TStopwatch watch;
   OscProb* osc_prob = OscProb::Instance(Ntoys);
   for(int n = 0; n < Ntoys; ++n)
   {
      Randomize(fparm);
      SetOscillationParameters(osc, fparm);
      if(isbar)
      {
         SetOscillationParameters(oscbar, fparm);
      }

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
            if(!isbar)
            {
               osc.CalcMatterOscProbAll(true, Etrue, p_osc_bar);
            }
            else
            {
               oscbar.CalcMatterOscProbAll(true, Etrue, p_osc_bar);
            }
               
            for(int im=0; im < knm; ++im)
            {
               SuperkMode_t m = kSKReactionModes[is][im];
               if(m == conv::kSKmUndefined)
               {
                  continue;
               }
               if(kIsCC[is][im])
               {
                  if(util::IsNuMu(m))
                  {
                     osc_prob->GetToys(s)->SetEntry(p_osc[1][1], n, it - 1, kFNuMu);
                  }
                  else if(util::IsNuMuBar(m))
                  {
                     osc_prob->GetToys(s)->SetEntry(p_osc_bar[1][1], n, it - 1, kFNuMuBar);
                  }
                  else if(util::IsNuE(m))
                  {
                     osc_prob->GetToys(s)->SetEntry(p_osc[0][0], n, it - 1, kFNuE);
                  }
                  else if(util::IsNuEBar(m))
                  {
                     osc_prob->GetToys(s)->SetEntry(p_osc_bar[0][0], n, it - 1, kFNuEBar);
                  }
                  else if(util::IsOscNuMu(m))
                  {
                     osc_prob->GetToys(s)->SetEntry(p_osc[0][1], n, it - 1, kFOscNuMu);
                  }
                  else if(util::IsOscNuMuBar(m))
                  {
                     osc_prob->GetToys(s)->SetEntry(p_osc_bar[0][1], n, it - 1, kFOscNuMuBar);
                  }
                  else if(util::IsOscNuE(m))
                  {
                     osc_prob->GetToys(s)->SetEntry(p_osc[1][1], n, it - 1, kFOscNuE);
                  }
                  else if(util::IsOscNuEBar(m))
                  {
                     osc_prob->GetToys(s)->SetEntry(p_osc[1][1], n, it - 1, kFOscNuEBar);
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
//   osc_prob.Write();
}

}
}

