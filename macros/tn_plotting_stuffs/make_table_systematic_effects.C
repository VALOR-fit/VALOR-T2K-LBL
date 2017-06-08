#include <iomanip.h>
#include <TXMLEngine.h>

void SetRatioNullTo1(TH1D *hratiomy){
  for(int ibin=0;ibin<hratiomy->GetNbinsX();ibin++){
    double cont = hratiomy->GetBinContent(ibin+1);
    if(cont==0)
      hratiomy->SetBinContent(ibin+1,1.);
  }
}

void SetRatio(TH1D* hnom,TH1D* htweak,TH1D* hratio){
  for(int ibin=0;ibin<hnom->GetNbinsX();ibin++){
    double cont_nom   = hnom->GetBinContent(ibin+1);
    double cont_tweak = htweak->GetBinContent(ibin+1);
    double cont_ratio = -999.;
    if(cont_nom>0.){
      cont_ratio = cont_tweak / cont_nom; 
    }
    else{
      cont_ratio = 1.;
    }
    hratio->SetBinContent(ibin+1,cont_ratio);
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

  // Define the tweak histograms
    
  TH1D *hEvents;
   
   int seed = 100000;
   int N = 10;
   // Xsec variations
   double* xsecEvents = new double[N];
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
      hEvents = new TH1D(*skexp.CurrPredictedErecoSpectrum(s));
      hEvents->SetName("hEvents");
      hEvents->SetTitle("XSec");
      xsecEvents[e] = hEvents->GetEntries();
      delete hEvents;
   }
   double meanXsec = TMath::Mean(N, xsecEvents);
   double rmsXsec = TMath::RMS(N, xsecEvents);
   double fracXsec = 100 * rmsXsec / meanXsec;


   // Flux variations
   seed = 200000;
   double* fluxEvents = new double[N];
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
      hEvents = new TH1D(*skexp.CurrPredictedErecoSpectrum(s));
      hEvents->SetName("hEvents");
      hEvents->SetTitle("Flux");
      fluxEvents[e] = hEvents->GetEntries();
      delete hEvents;
   }
   double meanFlux = TMath::Mean(N, fluxEvents);
   double rmsFlux = TMath::RMS(N, fluxEvents);
   double fracFlux = 100 * rmsFlux / meanFlux;
  
   // SK variations
   seed = 300000;
   double* skEvents = new double[N];
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
      hEvents = new TH1D(*skexp.CurrPredictedErecoSpectrum(s));
      hEvents->SetName("hEvents");
      hEvents->SetTitle("SK");
      skEvents[e] = hEvents->GetEntries();
      delete hEvents;
   }
   double meanSK = TMath::Mean(N, skEvents);
   double rmsSK = TMath::RMS(N, skEvents);
   double fracSK = 100 * rmsSK / meanSK;

   // Flux+Xsec variations
   seed = 400000;
   double* banffEvents = new double[N];
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
      hEvents = new TH1D(*skexp.CurrPredictedErecoSpectrum(s));
      hEvents->SetName("hEvents");
      hEvents->SetTitle("Flux+XSec");
      banffEvents[e] = hEvents->GetEntries();
      delete hEvents;
   }
   double meanBANFF = TMath::Mean(N, banffEvents);
   double rmsBANFF = TMath::RMS(N, banffEvents);
   double fracBANFF = 100 * rmsBANFF / meanBANFF;

   // All variations
   seed = 500000;
   double* allEvents = new double[N];
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
      hEvents = new TH1D(*skexp.CurrPredictedErecoSpectrum(s));
      hEvents->SetName("hEvents");
      hEvents->SetTitle("All");
      allEvents[e] = hEvents->GetEntries();
      delete hEvents;
   }
   double meanAll = TMath::Mean(N, allEvents);
   double rmsAll = TMath::RMS(N, allEvents);
   double fracAll = 100 * rmsAll / meanAll;

   stringstream streamer;
   streamer << "\\begin{tabular}{|l||l|l|l||l|l|l|}\n\\hline" << endl;
   streamer << "  & \\multicolumn{3}{|c||}{\\textbf{Pre-BANFF}} & " <<
      "\\multicolumn{3}{|c|}{\\textbf{Post-BANFF}} \\\\ \\hline" << endl;
   streamer << "  Group & Mean & $1\\sigma$ & \\% & Mean & $1\\sigma$ & \\% \\\\ \\hline\\hline" << endl;
   streamer << fixed << setprecision(2);
   streamer << "Flux & " << meanFlux << " & " << rmsFlux << " & " << fracFlux << " & " <<
      meanFlux << " & " << rmsFlux << " & " << fracFlux << " \\\\ \\hline" << endl;
   streamer << "Xsec & " << meanXsec << " & " << rmsXsec << " & " << fracXsec << " & " <<
      meanFlux << " & " << rmsXsec << " & " << fracXsec << " \\\\ \\hline" << endl;
   streamer << "Flux+Xsec & " << meanBANFF << " & " << rmsBANFF << " & " << fracBANFF << " & " <<
      meanBANFF << " & " << rmsBANFF << " & " << fracBANFF << " \\\\ \\hline" << endl;
   streamer << "SK+FSI+SI+PN & " << meanSK << " & " << rmsSK << " & " << fracSK << " & " <<
      meanSK << " & " << rmsSK << " & " << fracSK << " \\\\ \\hline\\hline" << endl;
   streamer << "All & " << meanAll << " & " << rmsAll << " & " << fracAll << " & " <<
      meanAll << " & " << rmsAll << " & " << fracAll << " \\\\ \\hline" << endl;
   streamer << "\\end{tabular}" << endl;

   long size;
   streamer.seekg(0, ifstream::end);
   size=streamer.tellg();
   streamer.seekg(0);
   ofstream outfile (output_file, ios_base::app);
   outfile.write (streamer.str().c_str(), size);
   outfile.close();

   delete[] fluxEvents;
   delete[] xsecEvents;
   delete[] skEvents;
   delete[] banffEvents;
   delete[] allEvents;

  return;
}


