#include <iomanip.h>
#include <TXMLEngine.h>

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

 const char * config_file        = "analysis_configuration_oa2017_run1-8_postfit.xml",
 
 const char * mixparamchoice     = "sinsqtheta",
 const char * hierarchy          = "normal",
 
 const char*  output_file        = "output.pdf"
 ){
  using namespace std;
  using namespace valor;
  using namespace valor::conv;

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

  // Analysis configuration
  const char * config  = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);
  //T2KOpt::Instance()->ReadXmlPath(config);
  string sconfig(config);
  T2KOpt::Instance()->ReadXmlPath(sconfig);

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
  if(ph.Hypothesis() == conv::kHypNuMuDisapp2F) {
    ph.SetCurrValue(conv::kPhypMixing23,  input_mixing23);
    ph.SetCurrValue(conv::kPhypDm23sq,    input_dm23sq);
  } else
    if(ph.Hypothesis() == conv::kHypOsc3F) {
      ph.SetCurrValue(conv::kPhypMixing12,  input_mixing12);
      ph.SetCurrValue(conv::kPhypMixing13,  input_mixing13);
      ph.SetCurrValue(conv::kPhypMixing23,  input_mixing23);
      ph.SetCurrValue(conv::kPhypCPV    ,   input_deltacp);
      ph.SetCurrValue(conv::kPhypDm12sq,    input_dm12sq);
      ph.SetCurrValue(conv::kPhypDm23sq,    input_dm23sq);
    } else
      if(ph.Hypothesis() == conv::kHypOsc3FBar){
        ph.SetCurrValue(conv::kPhypMixing12,  input_mixing12);
        ph.SetCurrValue(conv::kPhypMixing13,  input_mixing13);
        ph.SetCurrValue(conv::kPhypMixing23,  input_mixing23);
        ph.SetCurrValue(conv::kPhypCPV    ,   input_deltacp);
        ph.SetCurrValue(conv::kPhypDm12sq,    input_dm12sq);
        ph.SetCurrValue(conv::kPhypDm23sq,    input_dm23sq);
        ph.SetCurrValue(conv::kPhypMixing12Bar,  input_mixing12);
        ph.SetCurrValue(conv::kPhypMixing13Bar,  input_mixing13);
        ph.SetCurrValue(conv::kPhypMixing23Bar,  input_mixing23);
        ph.SetCurrValue(conv::kPhypCPVBar    ,   input_deltacp);
        ph.SetCurrValue(conv::kPhypDm12sqBar,    input_dm12sq);
        ph.SetCurrValue(conv::kPhypDm23sqBar,    input_dm23sq);
      }
    
  // Print options
  cout << "\n[Job inputs]" << endl;
  cout << "- Reading analysis configuration from: " << sconfig << endl;
  cout << "- Physics hypothesis: \n" << ph << endl;
  cout << *Opt::Instance() << endl;
    
  //output file name
  string oscparams =
    util::AsString(ph.Hypothesis(),input_mixing12,input_mixing13,input_mixing23,input_deltacp,input_dm12sq,input_dm23sq);
    
   // Create & config superk spectrum predictor
   SuperkExpectation skexp;
   skexp.Read(sconfig);

   stringstream streamer;
   streamer << fixed << setprecision(2);

   int N = 10;
 
   // SK variations
   double* skEvents = new double[N];
   VarySKFSI(skexp, s, ph, 300000, N, skEvents);
   double meanSK = TMath::Mean(N, skEvents);
   double rmsSK = TMath::RMS(N, skEvents);
   double fracSK = 100 * rmsSK / meanSK;
   streamer << "SK+FSI+SI+PN & " << meanSK << " & " << rmsSK << " & " << fracSK << " & " <<
      meanSK << " & " << rmsSK << " & " << fracSK << " \\\\ \\hline\\hline" << endl;

   // Constrained variations
   double* banffEvents = new double[N];
   VaryConstrained(skexp, s, ph, 400000, N, banffEvents);
   double meanBANFF = TMath::Mean(N, banffEvents);
   double rmsBANFF = TMath::RMS(N, banffEvents);
   double fracBANFF = 100 * rmsBANFF / meanBANFF;
   streamer << "Flux+Xsec constrained & " << meanBANFF << " & " << rmsBANFF << " & " << fracBANFF << " & " <<
      meanBANFF << " & " << rmsBANFF << " & " << fracBANFF << " \\\\ \\hline" << endl;

   // Xsec numu->nue variations
   double* xsecnueEvents = new double[N];
   VaryXsecNumuNue(skexp, s, ph, 600000, N, xsecnueEvents);
   double meanXsecnue = TMath::Mean(N, xsecnueEvents);
   double rmsXsecnue = TMath::RMS(N, xsecnueEvents);
   double fracXsecnue = 100 * rmsXsecnue / meanXsecnue;
   streamer << "$\\sigma(\\nu_e)/\\sigma(\\bar{\\nu}_e)$ & " <<
      meanXsecnue << " & " << rmsXsecnue << " & " << fracXsecnue << " & " <<
      meanXsecnue << " & " << rmsXsecnue << " & " << fracXsecnue << " \\\\ \\hline" << endl;

   // Xsec NCoth variations
   double* xsecncothEvents = new double[N];
   VaryXsecNCoth(skexp, s, ph, 700000, N, xsecncothEvents);
   double meanXsecncoth = TMath::Mean(N, xsecncothEvents);
   double rmsXsecncoth = TMath::RMS(N, xsecncothEvents);
   double fracXsecncoth = 100 * rmsXsecncoth / meanXsecncoth;
   streamer << "NC Other & " << meanXsecncoth << " & " << rmsXsecncoth << " & " << fracXsecncoth << " & " <<
      meanXsecncoth << " & " << rmsXsecncoth << " & " << fracXsecncoth << " \\\\ \\hline" << endl;

   // Xsec NCgamma variations
   double* xsecncgammaEvents = new double[N];
   VaryXsecNCGamma(skexp, s, ph, 800000, N, xsecncgammaEvents);
   double meanXsecncgamma = TMath::Mean(N, xsecncgammaEvents);
   double rmsXsecncgamma = TMath::RMS(N, xsecncgammaEvents);
   double fracXsecncgamma = 100 * rmsXsecncgamma / meanXsecncgamma;
   streamer << "NC1$\\gamma$ & " << meanXsecncgamma << " & " << rmsXsecncgamma << " & " << fracXsecncgamma << " & " <<
      meanXsecncgamma << " & " << rmsXsecncgamma << " & " << fracXsecncgamma << " \\\\ \\hline" << endl;

   // All variations
   double* allEvents = new double[N];
   VaryAll(skexp, s, ph, 500000, N, allEvents);
   double meanAll = TMath::Mean(N, allEvents);
   double rmsAll = TMath::RMS(N, allEvents);
   double fracAll = 100 * rmsAll / meanAll;
   streamer << "All & " << meanAll << " & " << rmsAll << " & " << fracAll << " & " <<
      meanAll << " & " << rmsAll << " & " << fracAll << " \\\\ \\hline" << endl;

   long size;
   streamer.seekg(0, ifstream::end);
   size=streamer.tellg();
   streamer.seekg(0);
   ofstream outfile (output_file);
   outfile.write (streamer.str().c_str(), size);
   outfile.close();

   delete[] skEvents;
   delete[] banffEvents;
   delete[] xsecnueEvents;
   delete[] xsecncothEvents;
   delete[] xsecncgammaEvents;
   delete[] allEvents;

   return;
}


