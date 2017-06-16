#include <iomanip.h>
#include <TXMLEngine.h>

void ResetOscParams(valor::PhysHypothesis& ph, double mix12val, double mix13val, double mix23val, double dm12val, double dm23val, double cpvval)
{
   ph.SetCurrValue(conv::kPhypMixing12, mix12val);
   ph.SetCurrValue(conv::kPhypMixing13, mix13val);
   ph.SetCurrValue(conv::kPhypMixing23, mix23val);
   ph.SetCurrValue(conv::kPhypDm12sq, dm12val);
   ph.SetCurrValue(conv::kPhypDm23sq, dm23val);
   ph.SetCurrValue(conv::kPhypCPV, cpvval);
}

bool OscParamsPhysical(double mix12val, double mix13val, double mix23val, double dm12val, double dm23val, double cpvval)
{
   using namespace valor::conv;
   if(mix12val < 0 || mix12val > 1)
   {
      return false;
   }
   if(mix13val < 0 || mix13val > 1)
   {
      return false;
   }
   if(mix23val < 0 || mix23val > 1)
   {
      return false;
   }
   if(dm12val < kMinFitDm12sq || dm12val > kMaxFitDm12sq)
   {
      return false;
   }
   if(dm23val < kMinFitDm23sq || dm23val > kMaxFitDm23sq)
   {
      return false;
   }
   if(cpvval < kPrior_deltaCP_Flat_Low || cpvval > kPrior_deltaCP_Flat_Up)
   {
      return false;
   }

   return true;
}

void RandomizeOsc(valor::PhysHypothesis& ph)
{
   using namespace valor::conv;
   for( ; ; )
   {
      double mix12val = Rnd::Instance()->RndToyMCosc().Gaus(kPrior_Sinsq2th12_Used_Mean, kPrior_Sinsq2th12_Used_RMS);
      double mix13val = Rnd::Instance()->RndToyMCosc().Gaus(kPrior_Sinsq2th13_Used_Mean, kPrior_Sinsq2th13_Used_RMS);
      double mix23val = Rnd::Instance()->RndToyMCosc().Uniform(kPrior_Sinsq2th23_Flat_Low, kPrior_Sinsq2th23_Flat_Up);
      double dm12val = Rnd::Instance()->RndToyMCosc().Gaus(kPrior_Dm12sq_Used_Mean, kPrior_Dm12sq_Used_RMS);
      double dm23val = Rnd::Instance()->RndToyMCosc().Uniform(kPrior_Dm23sq_Flat_Low, kPrior_Dm23sq_Flat_Up);
      double cpvval = Rnd::Instance()->RndToyMCosc().Uniform(kPrior_deltaCP_Flat_Low, kPrior_deltaCP_Flat_Up);

      mix12val = util::FromSinsq2thToSinsqth(mix12val);
      mix13val = util::FromSinsq2thToSinsqth(mix13val);
      mix23val = util::FromSinsq2thToSinsqth(mix23val);

      // Leave the atmospheric parameters fixed
      ph.SetCurrValue(conv::kPhypMixing12, mix12val);
      ph.SetCurrValue(conv::kPhypMixing13, mix13val);
//      ph.SetCurrValue(conv::kPhypMixing23, mix23val);
//      ph.SetCurrValue(conv::kPhypCPV, cpvval);
      ph.SetCurrValue(conv::kPhypDm12sq, dm12val);
//      ph.SetCurrValue(conv::kPhypDm23sq, dm23val);
      if(OscParamsPhysical(mix12val, mix13val, mix23val, dm12val, dm23val, cpvval))
      {
         break;
      }
   }
}

void VaryOsc(valor::SuperkExpectation& skexp, valor::conv::SuperkSample_t s, valor::PhysHypothesis& ph, const int seed, const int N, double* events)
{
   for(int e = 0; e < N; ++e)
   {
      // stat, syst, osc, fitter, other
      Rnd::Instance()->SetSeeds(seed + e + 10, seed + e + 20, seed + e + 30, seed + e + 40, seed + e + 50);

      skexp.ResetSystematics();
      RandomizeOsc(ph);

      skexp.Generate(ph);
      TH1D *hEvents = new TH1D(*skexp.CurrPredictedErecoSpectrum(s));
      events[e] = hEvents->GetEntries();
      delete hEvents;
   }
}

void VaryFlux(valor::SuperkExpectation& skexp, valor::conv::SuperkSample_t s, valor::PhysHypothesis& ph, const int seed, const int N, double* events)
{
   for(int e = 0; e < N; ++e)
   {
      // stat, syst, osc, fitter, other
      Rnd::Instance()->SetSeeds(seed + e + 10, seed + e + 20, seed + e + 30, seed + e + 40, seed + e + 50);

      bool isPhysical = false;
      while(!isPhysical)
      {
         skexp.ResetSystematics();
         skexp.GetBANFF()->TweakSampleFluxOnly(0, 0);

         PhysParamList paramTweaks = skexp.Tweaks();
         isPhysical = valor::T2KPhysicalParameterLimits::Instance()->AllWithinLimits(paramTweaks);
      }

      skexp.Generate(ph);
      TH1D *hEvents = new TH1D(*skexp.CurrPredictedErecoSpectrum(s));
      hEvents->SetName("hEvents");
      hEvents->SetTitle("Flux");
      events[e] = hEvents->GetEntries();
      delete hEvents;
   }
}

void VaryXsec(valor::SuperkExpectation& skexp, valor::conv::SuperkSample_t s, valor::PhysHypothesis& ph, const int seed, const int N, double* events)
{
   for(int e = 0; e < N; ++e)
   {
      // stat, syst, osc, fitter, other
      Rnd::Instance()->SetSeeds(seed + e + 10, seed + e + 20, seed + e + 30, seed + e + 40, seed + e + 50);

      bool isPhysical = false;
      while(!isPhysical)
      {
         skexp.ResetSystematics();
         skexp.TweakAllUncorrelatedXSecSystematics();
         skexp.GetNuXSec()->SyncRHCToFHCSystematics(); // uncorr xsec are fully correlated rhc/fhc
         skexp.GetBANFF()->TweakSampleXSecOnly(0, 0);

         PhysParamList paramTweaks = skexp.Tweaks();
         isPhysical = valor::T2KPhysicalParameterLimits::Instance()->AllWithinLimits(paramTweaks);
      }
      skexp.Generate(ph);
      TH1D* hEvents = new TH1D(*skexp.CurrPredictedErecoSpectrum(s));
      hEvents->SetName("hEvents");
      hEvents->SetTitle("XSec");
      events[e] = hEvents->GetEntries();
      delete hEvents;
   }

}

void VaryFluxAndXsec(valor::SuperkExpectation& skexp, valor::conv::SuperkSample_t s, valor::PhysHypothesis& ph, const int seed, const int N, double* events)
{
   for(int e = 0; e < N; ++e)
   {
      // stat, syst, osc, fitter, other
      Rnd::Instance()->SetSeeds(seed + e + 10, seed + e + 20, seed + e + 30, seed + e + 40, seed + e + 50);

      bool isPhysical = false;
      while(!isPhysical)
      {
         skexp.ResetSystematics();
         skexp.TweakAllUncorrelatedXSecSystematics();
         skexp.GetNuXSec()->SyncRHCToFHCSystematics(); // uncorr xsec are fully correlated rhc/fhc
         skexp.TweakAllBANFFSystematics();

         PhysParamList paramTweaks = skexp.Tweaks();
         isPhysical = valor::T2KPhysicalParameterLimits::Instance()->AllWithinLimits(paramTweaks);
      }

      skexp.Generate(ph);
      TH1D* hEvents = new TH1D(*skexp.CurrPredictedErecoSpectrum(s));
      hEvents->SetName("hEvents");
      hEvents->SetTitle("Flux+XSec");
      events[e] = hEvents->GetEntries();
      delete hEvents;
   }
}

void VaryConstrained(valor::SuperkExpectation& skexp, valor::conv::SuperkSample_t s, valor::PhysHypothesis& ph, const int seed, const int N, double* events)
{
   for(int e = 0; e < N; ++e)
   {
      // stat, syst, osc, fitter, other
      Rnd::Instance()->SetSeeds(seed + e + 10, seed + e + 20, seed + e + 30, seed + e + 40, seed + e + 50);

      bool isPhysical = false;
      while(!isPhysical)
      {
         skexp.ResetSystematics();
         skexp.GetBANFF()->TweakSampleConstrainedOnly();
         PhysParamList paramTweaks = skexp.Tweaks();
         isPhysical = valor::T2KPhysicalParameterLimits::Instance()->AllWithinLimits(paramTweaks);
      }

      skexp.Generate(ph);
      TH1D* hEvents = new TH1D(*skexp.CurrPredictedErecoSpectrum(s));
      events[e] = hEvents->GetEntries();
      delete hEvents;
   }
}

void VaryXsecNumuNue(valor::SuperkExpectation& skexp, valor::conv::SuperkSample_t s, valor::PhysHypothesis& ph, const int seed, const int N, double* events)
{
   for(int e = 0; e < N; ++e)
   {
      // stat, syst, osc, fitter, other
      Rnd::Instance()->SetSeeds(seed + e + 10, seed + e + 20, seed + e + 30, seed + e + 40, seed + e + 50);

      bool isPhysical = false;
      while(!isPhysical)
      {
         skexp.ResetSystematics();
         skexp.GetBANFF()->TweakSampleXSecNumuNueOnly();
         PhysParamList paramTweaks = skexp.Tweaks();
         isPhysical = valor::T2KPhysicalParameterLimits::Instance()->AllWithinLimits(paramTweaks);
      }

      skexp.Generate(ph);
      TH1D* hEvents = new TH1D(*skexp.CurrPredictedErecoSpectrum(s));
      events[e] = hEvents->GetEntries();
      delete hEvents;
   }
}

void VaryXsecNCGamma(valor::SuperkExpectation& skexp, valor::conv::SuperkSample_t s, valor::PhysHypothesis& ph, const int seed, const int N, double* events)
{
   for(int e = 0; e < N; ++e)
   {
      // stat, syst, osc, fitter, other
      Rnd::Instance()->SetSeeds(seed + e + 10, seed + e + 20, seed + e + 30, seed + e + 40, seed + e + 50);

      bool isPhysical = false;
      while(!isPhysical)
      {
         skexp.ResetSystematics();
         skexp.GetBANFF()->TweakSampleNCGamma();
         PhysParamList paramTweaks = skexp.Tweaks();
         isPhysical = valor::T2KPhysicalParameterLimits::Instance()->AllWithinLimits(paramTweaks);
      }

      skexp.Generate(ph);
      TH1D* hEvents = new TH1D(*skexp.CurrPredictedErecoSpectrum(s));
      events[e] = hEvents->GetEntries();
      delete hEvents;
   }
}

void VaryXsecNCoth(valor::SuperkExpectation& skexp, valor::conv::SuperkSample_t s, valor::PhysHypothesis& ph, const int seed, const int N, double* events)
{
   for(int e = 0; e < N; ++e)
   {
      // stat, syst, osc, fitter, other
      Rnd::Instance()->SetSeeds(seed + e + 10, seed + e + 20, seed + e + 30, seed + e + 40, seed + e + 50);

      bool isPhysical = false;
      while(!isPhysical)
      {
         skexp.ResetSystematics();
         skexp.GetBANFF()->TweakSampleNCoth();
         PhysParamList paramTweaks = skexp.Tweaks();
         isPhysical = valor::T2KPhysicalParameterLimits::Instance()->AllWithinLimits(paramTweaks);
      }

      skexp.Generate(ph);
      TH1D* hEvents = new TH1D(*skexp.CurrPredictedErecoSpectrum(s));
      events[e] = hEvents->GetEntries();
      delete hEvents;
   }
}

void VarySKFSI(valor::SuperkExpectation& skexp, valor::conv::SuperkSample_t s, valor::PhysHypothesis& ph, const int seed, const int N, double* events)
{
   for(int e = 0; e < N; ++e)
   {
      // stat, syst, osc, fitter, other
      Rnd::Instance()->SetSeeds(seed + e + 10, seed + e + 20, seed + e + 30, seed + e + 40, seed + e + 50);

      bool isPhysical = false;
      while(!isPhysical)
      {
         skexp.ResetSystematics();
         skexp.TweakAllSKDetectorFSISystematics();

         PhysParamList paramTweaks = skexp.Tweaks();
         isPhysical = valor::T2KPhysicalParameterLimits::Instance()->AllWithinLimits(paramTweaks);
      }
      skexp.Generate(ph);
      TH1D* hEvents = new TH1D(*skexp.CurrPredictedErecoSpectrum(s));
      hEvents->SetName("hEvents");
      hEvents->SetTitle("SK");
      events[e] = hEvents->GetEntries();
      delete hEvents;
   }
}

void VaryAll(valor::SuperkExpectation& skexp, valor::conv::SuperkSample_t s, valor::PhysHypothesis& ph, const int seed, const int N, double* events)
{
   for(int e = 0; e < N; ++e)
   {
      // stat, syst, osc, fitter, other
      Rnd::Instance()->SetSeeds(seed + e + 10, seed + e + 20, seed + e + 30, seed + e + 40, seed + e + 50);

      bool isPhysical = false;
      while(!isPhysical)
      {
         skexp.ResetSystematics();
         skexp.RandomizeSystematics();

         PhysParamList paramTweaks = skexp.Tweaks();
         isPhysical = valor::T2KPhysicalParameterLimits::Instance()->AllWithinLimits(paramTweaks);
      }

      skexp.Generate(ph);
      TH1D* hEvents = new TH1D(*skexp.CurrPredictedErecoSpectrum(s));
      hEvents->SetName("hEvents");
      hEvents->SetTitle("All");
      events[e] = hEvents->GetEntries();
      delete hEvents;
   }
}

void VaryAllWithOsc(valor::SuperkExpectation& skexp, valor::conv::SuperkSample_t s, valor::PhysHypothesis& ph, const int seed, const int N, double* events)
{
   for(int e = 0; e < N; ++e)
   {
      // stat, syst, osc, fitter, other
      Rnd::Instance()->SetSeeds(seed + e + 10, seed + e + 20, seed + e + 30, seed + e + 40, seed + e + 50);

      bool isPhysical = false;
      while(!isPhysical)
      {
         skexp.ResetSystematics();
         skexp.RandomizeSystematics();

         PhysParamList paramTweaks = skexp.Tweaks();
         isPhysical = valor::T2KPhysicalParameterLimits::Instance()->AllWithinLimits(paramTweaks);
      }
      RandomizeOsc(ph);

      skexp.Generate(ph);
      TH1D* hEvents = new TH1D(*skexp.CurrPredictedErecoSpectrum(s));
      hEvents->SetName("hEvents");
      hEvents->SetTitle("All");
      events[e] = hEvents->GetEntries();
      delete hEvents;
   }
}

void make_table_systematic_effects
(
 const char * sk_sample          = "MultiRe",      /* options: "1Rmu", "1Re", "RHC1Rmu", "RHC1Re", "MultiRmu, MultiRe" */
 const char * hypothesis         = "3f_osc",       /* options: "2f_numu_disap", "3f_osc", "3f_oscbar" , "nuosc_3f+nuisance"*/

 double       input_mixing12     = 0.304,//-99,    
 double       input_mixing13     = 0.0217, //-99,    
 double       input_mixing23     = 0.528,    
 double       input_deltacp      = -1.601,   
 double       input_dm12sq       = 7.53E-5,   
 double       input_dm23sq       = 2.509E-3,     

 const char * mixparamchoice     = "sinsqtheta",
 const char * hierarchy          = "normal",

 const char * config_file        = "analysis_configuration_oa2017_sk_error_all",

 const char * group              = "all",
 
 const char*  output_file        = "output"
 ){
   using namespace std;
   using namespace valor;
   using namespace valor::conv;

   std::string group_str(group);

   // Global options
   Opt::Instance()->gMassHierarchyAssumed = util::HierarchyFromString(hierarchy);
   Opt::Instance()->gOscMixingParamChoice = util::MixingTypeFromString(mixparamchoice);
   if(Opt::Instance()->Undefined()) {
      cerr << "Global options not defined correctly. Exiting..." << endl;
      cerr << *Opt::Instance() << endl;
      exit(1);
   }
    
   // Get the sample
   conv::SuperkSample_t s = util::SKSampleFromString(sk_sample); // requested sample
   if(s == conv::kSKsUndefined) {
      cerr << "** Undefined sample!" << endl;
      return;
   }
   int is = util::ArrIdx::Instance()->SuperkSampleIdx(s); // index of specified sample

   cout << endl;
   cout << "Running sample: " << util::AsString(s) << endl;
   cout << endl;

   // Set physics hypothesis
   if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsq2Theta) {
      if(input_mixing13 < -1) input_mixing13 = conv::kGlobFitSinsq2Theta13;
      if(input_mixing23 < -1) input_mixing23 = conv::kGlobFitSinsq2Theta23;
      if(input_mixing12 < -1) input_mixing12 = conv::kGlobFitSinsq2Theta12;
   }
   else {
      if(input_mixing13 < -1) input_mixing13 = conv::kGlobFitSinsqTheta13;
      if(input_mixing23 < -1) input_mixing23 = conv::kGlobFitSinsqTheta23;
      if(input_mixing12 < -1) input_mixing12 = conv::kGlobFitSinsqTheta12;
   }
   if(input_dm23sq < -1)   input_dm23sq  = conv::kGlobFitDm23sq;
   if(input_dm12sq < -1)   input_dm12sq  = conv::kGlobFitDm12sq;
   if(input_deltacp < -90) input_deltacp = conv::kGlobFitDeltaCP;
   PhysHypothesis ph(hypothesis);

   if(ph.Hypothesis() == conv::kHypOsc3F)
   {
      ResetOscParams(ph, input_mixing12, input_mixing13, input_mixing23, input_dm12sq, input_dm23sq, input_deltacp);
   }
   else
   {
      cout << "Unsupported hypothesis: " << hypothesis << endl;
      exit(0);
   }

   std::string banff[] = {"prefit", "postfit"};
   stringstream streamer;
   streamer << fixed << setprecision(2);

   for(int b = 0; b < 2; ++b)
   {
      // Analysis configuration
      const char * config  = Form("%s/config/%s_%s.xml", gSystem->Getenv("T2KOSCANA"), config_file, banff[b].c_str());
      string sconfig(config);
      T2KOpt::Instance()->ReadXmlPath(sconfig);
      
      // Print options
      cout << "\n[Job inputs]" << endl;
      cout << "- Reading analysis configuration from: " << sconfig << endl;
      cout << "- Physics hypothesis: \n" << ph << endl;
      cout << *Opt::Instance() << endl;
       
      // Create & config superk spectrum predictor
      SuperkExpectation skexp;
      skexp.Read(sconfig);

      int N = 10000;
      double* events = new double[N];
    
      if(group_str == "skdet" || group_str == "skfsi")
      {
         // SK variations
         VarySKFSI(skexp, s, ph, 300000, N, events);
         double mean = TMath::Mean(N, events);
         double rms = TMath::RMS(N, events);
         double frac = 100 * rms / mean;
         std::string name;
         if(group_str == "skdet")
         {
            name = "SK Detector";
         }
         else
         {
            name = "SK FSI+SI+PN";
         }
         if(banff[b] == "prefit")
         {
            streamer << name << " & " << mean << " & " << rms << " & " << frac << " & ";
         }
         else
         {
            streamer << mean << " & " << rms << " & " << frac << " \\\\ \\hline\\hline" << endl;
         }
      }
      else if(group_str == "banff_constrained")
      {
         // Constrained variations
         VaryConstrained(skexp, s, ph, 400000, N, events);
         double mean = TMath::Mean(N, events);
         double rms = TMath::RMS(N, events);
         double frac = 100 * rms / mean;
         std::string name = "Flux+Xsec constrained";
         if(banff[b] == "prefit")
         {
            streamer << name << " & " << mean << " & " << rms << " & " << frac << " & ";
         }
         else
         {
            streamer << mean << " & " << rms << " & " << frac << " \\\\ \\hline\\hline" << endl;
         }
      }
      else if(group_str == "xsec_nue")
      {
         // Xsec numu->nue variations
         VaryXsecNumuNue(skexp, s, ph, 600000, N, events);
         double mean = TMath::Mean(N, events);
         double rms = TMath::RMS(N, events);
         double frac = 100 * rms / mean;
         std::string name = "$\\sigma(\\nu_e)/\\sigma(\\bar{\\nu}_e)$";
         if(banff[b] == "prefit")
         {
            streamer << name << " & " << mean << " & " << rms << " & " << frac << " & ";
         }
         else
         {
            streamer << mean << " & " << rms << " & " << frac << " \\\\ \\hline\\hline" << endl;
         }
      }
      else if (group_str == "xsec_ncoth")
      {
         // Xsec NCoth variations
         VaryXsecNCoth(skexp, s, ph, 700000, N, events);
         double mean = TMath::Mean(N, events);
         double rms = TMath::RMS(N, events);
         double frac = 100 * rms / mean;
         std::string name = "NC Other";
         if(banff[b] == "prefit")
         {
            streamer << name << " & " << mean << " & " << rms << " & " << frac << " & ";
         }
         else
         {
            streamer << mean << " & " << rms << " & " << frac << " \\\\ \\hline\\hline" << endl;
         }
      }
      else if (group_str == "xsec_ncgamma")
      {
         // Xsec NCgamma variations
         VaryXsecNCGamma(skexp, s, ph, 800000, N, events);
         double mean = TMath::Mean(N, events);
         double rms = TMath::RMS(N, events);
         double frac = 100 * rms / mean;
         std::string name = "NC1$\\gamma$";
         if(banff[b] == "prefit")
         {
            streamer << name << " & " << mean << " & " << rms << " & " << frac << " & ";
         }
         else
         {
            streamer << mean << " & " << rms << " & " << frac << " \\\\ \\hline\\hline" << endl;
         }
      }
      else if (group_str == "allsyst")
      {
         // All variations
         VaryAll(skexp, s, ph, 500000, N, events);
         double mean = TMath::Mean(N, events);
         double rms = TMath::RMS(N, events);
         double frac = 100 * rms / mean;
         std::string name = "All";
         if(banff[b] == "prefit")
         {
            streamer << name << " & " << mean << " & " << rms << " & " << frac << " & ";
         }
         else
         {
            streamer << mean << " & " << rms << " & " << frac << " \\\\ \\hline\\hline" << endl;
         }
      }
      else if (group_str == "all")
      {
         // All variations
         VaryAllWithOsc(skexp, s, ph, 500000, N, events);
         double mean = TMath::Mean(N, events);
         double rms = TMath::RMS(N, events);
         double frac = 100 * rms / mean;
         std::string name = "All with osc";
         if(banff[b] == "prefit")
         {
            streamer << name << " & " << mean << " & " << rms << " & " << frac << " & ";
         }
         else
         {
            streamer << mean << " & " << rms << " & " << frac << " \\\\ \\hline\\hline" << endl;
         }
      }
      else if (group_str == "osc")
      {
         // Osc variations
         VaryOsc(skexp, s, ph, 500000, N, events);
         double mean = TMath::Mean(N, events);
         double rms = TMath::RMS(N, events);
         double frac = 100 * rms / mean;
         std::string name = "Osc";
         if(banff[b] == "prefit")
         {
            streamer << name << " & " << mean << " & " << rms << " & " << frac << " & ";
         }
         else
         {
            streamer << mean << " & " << rms << " & " << frac << " \\\\ \\hline\\hline" << endl;
         }
      }

      if(banff[b] == "postfit")
      {
         const char* filename = Form("%s_dist_%s.csv", output_file, group);
         stringstream diststreamer;
         diststreamer << fixed << setprecision(2);
         for(int i = 0; i < N; ++i)
         {
            diststreamer << events[i] << endl;
         }
         diststreamer.seekg(0, ifstream::end);
         long size = diststreamer.tellg();
         diststreamer.seekg(0);
         ofstream outfile (filename);
         outfile.write (diststreamer.str().c_str(), size);
         outfile.close();
      }

      delete[] events;

      ResetOscParams(ph, input_mixing12, input_mixing13, input_mixing23, input_dm12sq, input_dm23sq, input_deltacp);
   }

   const char* filename = Form("%s_%s.tex", output_file, group);
   streamer.seekg(0, ifstream::end);
   long size = streamer.tellg();
   streamer.seekg(0);
   ofstream outfile (filename);
   outfile.write (streamer.str().c_str(), size);
   outfile.close();

   return;
}


