//
// ------------------------------------------------------------------------------------------------           
//
// Run toys for marginalization of oscillation parameters and store into binary files. 
//                                                                                           
// ------------------------------------------------------------------------------------------------

#include "valor/Conventions/Constants.h"

using namespace valor;
using namespace valor::conv;

PhysParamId_t GetFreeParameter(const char* parname)
{
   std::string name(parname);  
   if(name == "mixing23Bar")
   {
      return kPhypMixing23Bar;
   }
   else if(name == "mixing13Bar")
   {
      return kPhypMixing13Bar;
   }
   else if(name == "dmsq23Bar")
   {
      return kPhypDm23sqBar;
   }
   else if(name == "dmsq12Bar")
   {
      return kPhypDm12sqBar;
   }
   else if(name == "mixing12Bar")
   {
      return kPhypMixing12Bar;
   }
   else if(name == "deltacpBar")
   {
      return kPhypCPVBar;
   }
   else if(name == "betaBar")
   {
      return kPhypBetaBar;
   }
   else if(name == "mixing23")
   {
      return kPhypMixing23;
   } 
   else if(name == "mixing13")
   {
      return kPhypMixing13;
   }
   else if(name == "dmsq23")
   {
      return kPhypDm23sq;
   }
   else if(name == "dmsq12")
   {
      return kPhypDm12sq;
   }
   else if(name == "mixing12")
   {
      return kPhypMixing12;
   }
   else if(name == "deltacp")
   {
      return kPhypCPV;
   }
   else
   {
      return kPhypUndefined;
   }
}

void Initialize(FitParamList* fparm)
{
   PhysHypothesis_t hyp = fparm->PhysicsHypothesis();

   double mixing12;
   double mixing13;
   double mixing23;
   double dm12sq;
   double dm23sq;
   double deltacp;
   double beta_bar;

   if(Opt::Instance()->gOscMixingParamChoice == kNuMixSinsq2Theta)
   {
      mixing12 = kGlobFitSinsq2Theta12;
      mixing13 = kGlobFitSinsq2Theta13;
      mixing23 = kGlobFitSinsq2Theta23;
   }
   else
   {
      mixing12 = kGlobFitSinsqTheta12;
      mixing13 = kGlobFitSinsqTheta13;
      mixing23 = kGlobFitSinsqTheta23;
   }

   dm23sq = kGlobFitDm23sq;
   dm12sq = kGlobFitDm12sq;
   deltacp = kGlobFitDeltaCP;
   beta_bar = kGlobFitBetaBar;

   // set oscillation parameter inputs as necessary (used to generate toy-MC expt)
   if(hyp == kHypOsc3F || hyp == kHypOsc3FBar || hyp == kHypOsc3FNormBar ||
         hyp == kHypOsc3FBarNormBar || hyp == kHypOsc3FBarDissap)
   {
      fparm->SetInitValue(kPhypMixing12, mixing12);
      fparm->SetInitValue(kPhypMixing13, mixing13);
      fparm->SetInitValue(kPhypMixing23, mixing23);
      fparm->SetInitValue(kPhypCPV, deltacp);
      fparm->SetInitValue(kPhypDm12sq, dm12sq);
      fparm->SetInitValue(kPhypDm23sq, dm23sq);
   }

   if(hyp == kHypOsc3FBar || hyp == kHypOsc3FBarNormBar)
   {
      fparm->SetInitValue(kPhypMixing12Bar, mixing12);
      fparm->SetInitValue(kPhypMixing13Bar, mixing13);
      fparm->SetInitValue(kPhypMixing23Bar, mixing23);
      fparm->SetInitValue(kPhypCPVBar, deltacp);
      fparm->SetInitValue(kPhypDm12sqBar, dm12sq);
      fparm->SetInitValue(kPhypDm23sqBar, dm23sq);
   }
   
   if(hyp == kHypOsc3FNormBar || hyp == kHypOsc3FBarNormBar)
   {
      fparm->SetInitValue(kPhypBetaBar, beta_bar);
   }
   else if(hyp == kHypOsc3FBarDissap)
   {
      fparm->SetInitValue(kPhypMixing23Bar, mixing23);
      fparm->SetInitValue(kPhypDm23sqBar, dm23sq);
   }
   else if(hyp == kHypNuMuDisapp2F)
   {
      fparm->SetInitValue(kPhypMixing23, mixing23);
      fparm->SetInitValue(kPhypDm23sq, dm23sq);
   }
}

void SetMarginalization(FitParamList* fparm, PhysParamId_t freepar0, PhysParamId_t freepar1)
{
   const int N = fparm->Size();

   for(int i = 0; i < N; ++i)
   {
      PhysParamId_t p = fparm->WhichParam(i);
      if(p != freepar0 && p != freepar1)
      {
         fparm->SetDoMarg(p, true);
      }
   }
}

void SetFreeParameters(FitParamList* fparm, PhysParamId_t freepar0, double val0,
      PhysParamId_t freepar1, double val1)
{
   if(freepar0 != kPhypUndefined)
   {
      fparm->SetInitValue(freepar0, val0);
   }

   if(freepar1 != kPhypUndefined)
   {
      fparm->SetInitValue(freepar1, val1);
   }
}

void SetPrior(FitParamList* fparm, PhysParamId_t p, std::string prior, double a, double b)
{
   fparm->SetPriorMode(p, prior);
  
   if(prior == "flat")
   {
      fparm->SetMinValueFlat(p, a);
      fparm->SetMaxValueFlat(p, b);
      cout << "Setting flat prior " << util::AsString(p) << " Min = " <<
         fparm->GetMinValueFlat(p) << " Max = " << fparm->GetMaxValueFlat(p) << std::endl;
   }
   else if(prior == "gaus")
   {
      fparm->SetMeanGaus(p, a);
      fparm->SetSigmaGaus(p, b);
      cout << "Setting gaussian prior " << util::AsString(p) << " Mean = " <<
         fparm->GetMeanGaus(p) << " Sigma = " << fparm->GetSigmaGaus(p) << std::endl;
   }
}

void SetPriors(FitParamList* fparm, bool reactors)
{
   PhysHypothesis_t hyp = fparm->PhysicsHypothesis();

   // deltacp
   if(fparm->IsMarginalized(kPhypCPV))
   {
      SetPrior(fparm, kPhypCPV, "flat", kPrior_deltaCP_Flat_Low, kPrior_deltaCP_Flat_Up);
   }
   if((hyp == kHypOsc3FBarNormBar || hyp == kHypOsc3FBar) && fparm->IsMarginalized(kPhypCPVBar))
   { 
      SetPrior(fparm, kPhypCPVBar, "flat", kPrior_deltaCP_Flat_Low, kPrior_deltaCP_Flat_Up);
   }
    
   // mixing13
   if(reactors)
   {
      if(fparm->IsMarginalized(kPhypMixing13))
      {
         SetPrior(fparm, kPhypMixing13, "gaus", kPrior_Sinsq2th13_Used_Mean, kPrior_Sinsq2th13_Used_RMS);
      }
      if((hyp == kHypOsc3FBarNormBar || hyp == kHypOsc3FBar) && fparm->IsMarginalized(kPhypMixing13Bar))
      {
         SetPrior(fparm, kPhypMixing13Bar, "gaus", kPrior_Sinsq2th13_Used_Mean, kPrior_Sinsq2th13_Used_RMS);
      }      
   }
   else
   {
      if(fparm->IsMarginalized(kPhypMixing13))
      {
         SetPrior(fparm, kPhypMixing13, "flat", kPrior_Sinsq2th13_Flat_Low, kPrior_Sinsq2th13_Flat_Up);
      }
      if((hyp == kHypOsc3FBarNormBar || hyp == kHypOsc3FBar) && fparm->IsMarginalized(kPhypMixing13Bar))
      {
         SetPrior(fparm, kPhypMixing13Bar, "flat", kPrior_Sinsq2th13_Flat_Low, kPrior_Sinsq2th13_Flat_Up);
      }
   }

   // mixing23
   if(fparm->IsMarginalized(kPhypMixing23))
   {
      SetPrior(fparm, kPhypMixing23, "flat", kPrior_Sinsq2th23_Flat_Low, kPrior_Sinsq2th23_Flat_Up);
   }
   if((hyp == kHypOsc3FBarNormBar || hyp == kHypOsc3FBar || hyp == kHypOsc3FBarDissap) &&
         fparm->IsMarginalized(kPhypMixing23Bar))
   {
      SetPrior(fparm, kPhypMixing23Bar, "flat", kPrior_Sinsq2th23_Flat_Low, kPrior_Sinsq2th23_Flat_Up);
   }
    
   // dmsq23
   if(fparm->IsMarginalized(kPhypDm23sq))
   {
      SetPrior(fparm, kPhypDm23sq, "flat", kPrior_Dm23sq_Flat_Low, kPrior_Dm23sq_Flat_Up);
   }
   if((hyp == kHypOsc3FBarNormBar || hyp == kHypOsc3FBar || hyp == kHypOsc3FBarDissap) &&
         fparm->IsMarginalized(kPhypDm23sqBar))
   {
      SetPrior(fparm, kPhypDm23sqBar, "flat", kPrior_Dm23sq_Flat_Low, kPrior_Dm23sq_Flat_Up);
   }
 
   // mixing12
   if(fparm->IsMarginalized(kPhypMixing12))
   {
      SetPrior(fparm, kPhypMixing12, "gaus", kPrior_Sinsq2th12_Used_Mean, kPrior_Sinsq2th12_Used_RMS);
   }
 
   // dmsq12
   if(fparm->IsMarginalized(kPhypDm12sq))
   {
      SetPrior(fparm, kPhypDm12sq, "gaus", kPrior_Dm12sq_Used_Mean, kPrior_Dm12sq_Used_RMS);
   }
}

void test_osc_marginalization
(
 const int freepar0_bin = -999,
 const int freepar1_bin = -999,
 const char* tag = "osc_toys",
 const char* config_file = "analysis_configuration.xml",
 const char* fit_param_list = "nuosc3f",
 const char* freepar0_name = "mixing13",
 const double freepar0_val = 0., 
 const char* freepar1_name = "deltacp",
 const double freepar1_val = -3.14159265358979312,
 const int NtoysMarg = 10,
 // auto==kNuDmsqAuto abs23==kNuDmsq12andAbs23 abs31==kNuDmsq12andAbs31 fogli==kNuDmsqFogliLisi
 const char* dmsqscheme = "auto",
 const char* hierarchy = "normal",
 // sinsq2theta or sinsqtheta
 const char* mixparamchoice = "sinsqtheta",
 const char* mixparamoctant = "two",
 int seed = 999,
 bool reactors = false
)
{  
   FileStat_t dummy;
   TString config = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);

   if(gSystem->GetPathInfo(config.Data(), dummy))
   {
      config = config_file;
      if(gSystem->GetPathInfo(config.Data(), dummy))
      {
         cerr<<"Config file not found: "<<config_file<<endl;
         exit(-1);
      }
   }
 
   // Global options
   Opt::Instance()->gMassHierarchyAssumed = util::HierarchyFromString(hierarchy);
   Opt::Instance()->gOscDmsqScheme = util::DmsqSchemeFromString(dmsqscheme);
   Opt::Instance()->gOscMixingParamChoice = util::MixingTypeFromString(mixparamchoice);
   Opt::Instance()->gOscOctantChoice = util::OctantFromString(mixparamoctant);
   if(Opt::Instance()->Undefined())
   {
      cerr << "Global options not defined correctly. Exiting..." << endl;
      cerr << *Opt::Instance() << endl;
      exit(1);
   }
   Opt::Instance()->SetDefaultParameterValues();

   // Get superk spectrum predictor & configure
   SuperkExpectation * skexp = new SuperkExpectation;
   skexp->Read(config.Data());

   PhysParamId_t freepar0 = GetFreeParameter(freepar0_name);
   if(freepar0 == kPhypUndefined)
   {  // Parameter 0 must be always defined
      cout << "Parameter 0 is not defined" << endl;
      exit(1);
   }

   PhysParamId_t freepar1 = GetFreeParameter(freepar1_name);
   if(freepar1 == kPhypUndefined)
   {
      cout << "Parameter 1 is not defined" << endl;
      cout << "1D fit" << endl;
   }
   FitParamList* fparm = util::NamedParameterList(fit_param_list);
   // Order is important here, initialize first, then SetPriors should be called
   // after SetMarginalization (SetFreeParameters does not interact with SetPriors
   // or SetMarginalization)
   Initialize(fparm);
   SetFreeParameters(fparm, freepar0, freepar0_val, freepar1, freepar1_val);
   SetMarginalization(fparm, freepar0, freepar1);
   SetPriors(fparm, reactors);

   // stat, syst, osc, fitter, other
   Rnd::Instance()->SetSeeds(seed + 10, seed + 20, seed + 30, seed + 40, seed + 50);

   std::cout << *fparm << std::endl;

   util::CalcOscProbabilities(4000, fparm);
}
