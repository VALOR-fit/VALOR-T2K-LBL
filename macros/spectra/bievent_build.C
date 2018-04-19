//
// ------------------------------------------------------------------------------------------------
// Plot contributions from various event modes in the predicted SuperK 1Re or 1Rmu spectrum.
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

#include <iomanip.h>

using namespace std;
using namespace valor;
using namespace valor::conv; 

void bievent_build
(
 const int ntoys = 100,
 const int seed = 10000,
 const char * tag                  = "ellipse_Run1_9c",
 const char * sk_sample            = "1Re",                         /* options: "1Rmu", "1Re", "RHC1Rmu", "RHC1Re", "MultiRmu" */
 const char * hypothesis           = "3f_osc",                 /* options: "2f_numu_disap", "3f_osc, 3f_oscbar, 3f_osc_normbar, 3f_oscbar_normbar" */ 
 // Best fit values NH
 double       input_mixing12       = 0.304,
 double       input_mixing13       = 0.02235,
 double       input_mixing23       = 0.530,
 double       input_deltacp        = 0,
 double       input_dm12sq         = 7.53E-5,
 double       input_dm23sq         = 2.460E-3,
 double       input_betabar        = -99., //0., //-99,
 bool         disable_banff_tuning = false,                          /* */
 const char * config_file = "analysis_configuration_oa2018_run1-9c_v1_postfit.xml",
 const char * mixparamchoice      = "sinsqtheta",
 const char * hierarchy           = "normal"   ,
 const char * dmsqscheme          = "auto", /* auto==kNuDmsqAuto abs23==kNuDmsq12andAbs23 abs31==kNuDmsq12andAbs31 fogli==kNuDmsqFogliLisi */
 const char * mixparamoctant      = "two",       /* one==kNuOctantOne two==kNuOctantTwo */
 int          iexperiment         = 0,
 const char * rootfile = "official/2018/data/SK_data_Run19c_valor.root", 
 const bool draw_data_spectra = true,
 const bool scale_mc = false
)
{
   using namespace std;
   using namespace valor;

   // Global options
   Opt::Instance()->gMassHierarchyAssumed = util::HierarchyFromString(hierarchy);
   Opt::Instance()->gOscDmsqScheme        = util::DmsqSchemeFromString(dmsqscheme);
   Opt::Instance()->gOscMixingParamChoice = util::MixingTypeFromString(mixparamchoice);
   Opt::Instance()->gOscOctantChoice      = util::OctantFromString(mixparamoctant);
   if(Opt::Instance()->Undefined())
   {
      cerr << "Global options not defined correctly. Exiting..." << endl;
      cerr << *Opt::Instance() << endl;
      exit(1);
   }

   // Set physics hypothesis
   if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsq2Theta)
   {
      if(input_mixing13 < -1) input_mixing13 = conv::kGlobFitSinsq2Theta13;
      if(input_mixing23 < -1) input_mixing23 = conv::kGlobFitSinsq2Theta23;
      if(input_mixing12 < -1) input_mixing12 = conv::kGlobFitSinsq2Theta12;
   }
   else
   {
      if(input_mixing13 < -1) input_mixing13 = conv::kGlobFitSinsqTheta13;
      if(input_mixing23 < -1) input_mixing23 = conv::kGlobFitSinsqTheta23;
      if(input_mixing12 < -1) input_mixing12 = conv::kGlobFitSinsqTheta12;
   }
   if(input_dm23sq < -1)  input_dm23sq = conv::kGlobFitDm23sq;
   if(input_dm12sq < -1)  input_dm12sq = conv::kGlobFitDm12sq;
   if(input_deltacp < -90) input_deltacp = conv::kGlobFitDeltaCP;
   if(input_betabar < -1) input_betabar = conv::kGlobFitBetaBar;
 
   PhysHypothesis ph(hypothesis);
   if(ph.Hypothesis() == conv::kHypNuMuDisapp2F)
   {
      ph.SetCurrValue(conv::kPhypMixing23,  input_mixing23);
      ph.SetCurrValue(conv::kPhypDm23sq,    input_dm23sq);
   }
   else if(ph.Hypothesis() == conv::kHypOsc3F)
   {
      ph.SetCurrValue(conv::kPhypMixing12,  input_mixing12);
      ph.SetCurrValue(conv::kPhypMixing13,  input_mixing13);
      ph.SetCurrValue(conv::kPhypMixing23,  input_mixing23);
      ph.SetCurrValue(conv::kPhypCPV    ,   input_deltacp);
      ph.SetCurrValue(conv::kPhypDm12sq,    input_dm12sq);
      ph.SetCurrValue(conv::kPhypDm23sq,    input_dm23sq);
   }
   else if(ph.Hypothesis() == conv::kHypOsc3FBar)
   {
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
   } else if(ph.Hypothesis() == conv::kHypOsc3FNormBar)
   {
      ph.SetCurrValue(conv::kPhypMixing12,  input_mixing12);
      ph.SetCurrValue(conv::kPhypMixing13,  input_mixing13);
      ph.SetCurrValue(conv::kPhypMixing23,  input_mixing23);
      ph.SetCurrValue(conv::kPhypCPV    ,   input_deltacp);
      ph.SetCurrValue(conv::kPhypDm12sq,    input_dm12sq);
      ph.SetCurrValue(conv::kPhypDm23sq,    input_dm23sq);
      ph.SetCurrValue(conv::kPhypBetaBar,    input_betabar);
   } else if(ph.Hypothesis() == conv::kHypOsc3FBarNormBar)
   {
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
      ph.SetCurrValue(conv::kPhypBetaBar,    input_betabar);
   }
  
   string oscparams = util::AsString(ph.Hypothesis(), input_mixing12,
         input_mixing13, input_mixing23, input_deltacp, input_dm12sq,
         input_dm23sq, input_betabar, -999, -999);

   conv::SuperkSample_t s = util::SKSampleFromString(sk_sample); // requested sample
   if(s == conv::kSKsUndefined)
   {
      cerr << "** Undefined sample!" << endl;
      return;
   }
   int is = util::ArrIdx::Instance()->SuperkSampleIdx(s); // index of specified sample

   // Analysis configuration
   const char * config  = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);
   string sconfig(config);

   // Print options
   cout << "\n[Job inputs]" << endl;
   cout << "- Reading analysis configuration from: " << sconfig << endl;
   cout << "- Disabling BANFF tuning: " << ((disable_banff_tuning) ? "YES" : "NO") << endl;
   cout << "- Physics hypothesis: " << ph << endl;
   cout << *Opt::Instance() << endl;

   // Data files

   TString rootfile2(Form("%s/inputs/%s", gSystem->Getenv("T2KOSCANA"), rootfile));
   SuperkObservation * skobs_fakedata = new SuperkObservation;
   cout << "Use root file" << endl;
   bool draw_data = true;
   TFile *infile = new TFile(rootfile2.Data(),"READ");
   if(!infile->IsOpen())
   {
      cout << endl;
      cout << "The data root file doesn't exist: don't draw data spectra" << endl;
      cout << endl;
      draw_data = false;
   }
   if(!draw_data_spectra)
   {
      draw_data = false;
   }
   infile->Close();

   // Create & config superk spectrum predictor
   SuperkExpectation skexp;
   skexp.Read(sconfig);

   // Disable BANFF tuning if requested
   if(disable_banff_tuning)
   {
      cout << "- Performing special systematic parameter variation which cancels the effect of BANFF tuning..." << endl;
      skexp.GetBANFF()->RestoreNominalMC();
   }

   //
   // - Open loop over values of deltaCP (each octant and mass hierarchy)
   //   and fill TGraph with a point of N_fhc_elike Vs N_rhc_elike
   // - Set lines of the TGraph
   // - Calculate N_fhc_elike Vs N_rhc_elike point from data and plot it on top of the MC one
   // 
    
   double deltaCP_first = -TMath::Pi();
   double deltaCP_last = +TMath::Pi();
   const int dcp_npt = 9;
   double deltaCP_width = (deltaCP_last - deltaCP_first) / (dcp_npt-1);
   const int N23 = 3;
   const int NMH = 2;
   const int NDCP = 4;

   TH2D* hEventsSyst[dcp_npt][N23][NMH];
   TH2D* hEvents[dcp_npt][N23][NMH];
   TH1D* hFHCEvents[dcp_npt][N23][NMH];
   TH1D* hRHCEvents[dcp_npt][N23][NMH];
   for(int idcp=0; idcp < dcp_npt; idcp++)
   {
      for(int iss23 = 0; iss23 < N23; ++iss23)
      {
         for(int imh = 0; imh < NMH; ++imh)
         {
            const char* name_syst = Form("hEventsSyst_%d_%d_%d", idcp, iss23, imh);
            hEventsSyst[idcp][iss23][imh] = new TH2D(name_syst, name_syst, 600, 0, 150, 100, 0, 25);
            const char* name = Form("hEvents_%d_%d_%d", idcp, iss23, imh);
            hEvents[idcp][iss23][imh] = new TH2D(name, name, 150, 0, 150, 25, 0, 25);
            const char* fhcname = Form("hEvents_fhc_%d_%d_%d", idcp, iss23, imh);
            hFHCEvents[idcp][iss23][imh] = new TH1D(fhcname, fhcname, 150, 0, 150);
            const char* rhcname = Form("hEvents_rhc_%d_%d_%d", idcp, iss23, imh);
            hRHCEvents[idcp][iss23][imh] = new TH1D(rhcname, rhcname, 25, 0, 25);
         }
      }
   }

   double ss12[NMH] = {0.304, 0.304};
   double ss13[NMH] = {0.02235, 0.0226125};
   double dm12[NMH] = {7.53E-5, 7.53E-5};
   double dm23[NMH] = {2.460E-3, 2.440E-3};
   double ss23[N23] = {0.450, 0.500, 0.550};

   for(int e = 0; e < ntoys; ++e)
   {
      // stat, syst, osc, fitter, other
      if(e % 10 == 0)
      {
         cout << "Throwing toy " << e << endl;
      }
      Rnd::Instance()->SetSeeds(seed + e + 10, seed + e + 20, seed + e + 30, seed + e + 40, seed + e + 50);

      bool isPhysical = false;
      while(!isPhysical)
      {
         skexp.ResetSystematics();
         skexp.RandomizeSystematics();
         RandomizeOsc(ph);

         PhysParamList paramTweaks = skexp.Tweaks();
         isPhysical = valor::T2KPhysicalParameterLimits::Instance()->AllWithinLimits(paramTweaks);
      }

      for(int idcp=0; idcp < dcp_npt; idcp++)
      {
         for(int iss23 = 0; iss23 < N23; ++iss23)
         {
            for(int imh = 0; imh < NMH; ++imh)
            {
               double Nfhc = -999;
               double Nrhc = -999;
               double Ncc1pi = -999;
    
               double dcp = deltaCP_first + idcp * deltaCP_width;
       
               // Set NH oscillation parameter values
               // Solar parameters and theta13 are now marginalized
               /*
               ph.SetCurrValue(conv::kPhypMixing12, ss12[imh]);
               ph.SetCurrValue(conv::kPhypMixing13, ss13[imh]);
               ph.SetCurrValue(conv::kPhypDm12sq, dm12[imh]);
               */
               ph.SetCurrValue(conv::kPhypDm23sq, dm23[imh]);

               if(imh == 0)
               {
                  Opt::Instance()->gMassHierarchyAssumed = util::HierarchyFromString("normal");
                  skexp.GetSuperkMC()->SetMassHierarchy(conv::kNmhNormal);
               }
               else
               {
                  Opt::Instance()->gMassHierarchyAssumed = util::HierarchyFromString("inverted");
                  skexp.GetSuperkMC()->SetMassHierarchy(conv::kNmhInverted);
               }
               ph.SetCurrValue(conv::kPhypCPV, dcp);
               ph.SetCurrValue(conv::kPhypMixing23, ss23[iss23]);
               // Generate prediction

               skexp.GetSuperkMC()->SetX2p2h(1);
               skexp.Generate(ph);
               GetN_1re(ph,skexp,Nfhc,Nrhc);
               GetN_cc1pi(ph,skexp,Ncc1pi);
               Nfhc += Ncc1pi;
               hEventsSyst[idcp][iss23][imh]->Fill(Nfhc, Nrhc);

               // Apply statistical fluctuations
               Nfhc = Rnd::Instance()->RndToyMCstat().Poisson(Nfhc);
               Nrhc = Rnd::Instance()->RndToyMCstat().Poisson(Nrhc);
               hEvents[idcp][iss23][imh]->Fill(Nfhc, Nrhc);

               hFHCEvents[idcp][iss23][imh]->Fill(Nfhc);
               hRHCEvents[idcp][iss23][imh]->Fill(Nrhc);
            }
         }
      }
   }

   TFile file(Form("events_%d.root", seed), "RECREATE");
   for(int idcp=0; idcp < dcp_npt; ++idcp)
   {
      for(int iss23 = 0; iss23 < N23; ++iss23)
      {
         for(int imh = 0; imh < NMH; ++imh)
         {
            hEvents[idcp][iss23][imh]->Write();
            hEventsSyst[idcp][iss23][imh]->Write();
            hFHCEvents[idcp][iss23][imh]->Write();
            hRHCEvents[idcp][iss23][imh]->Write();
         }
      }
   }
   file.Close();

   return;
}

void GetN_1re(PhysHypothesis ph,SuperkExpectation &skexp,double &Nfhc,double &Nrhc)
{
   // 1Re # of events 
   conv::SuperkSample_t s = util::SKSampleFromString("1Re"); // requested sample
   if(s == conv::kSKsUndefined)
   {
      cerr << "** Undefined sample!" << endl;
      return;
   }
   int is = util::ArrIdx::Instance()->SuperkSampleIdx(s); // index of specified sample
   Nfhc = skexp.CurrPredictedErecoVsQ2Spectrum(s)->Integral();
  
   // RHC1Re # of events 
   s = util::SKSampleFromString("RHC1Re"); // requested sample
   if(s == conv::kSKsUndefined)
   {
      cerr << "** Undefined sample!" << endl;
      return;
   }
   is = util::ArrIdx::Instance()->SuperkSampleIdx(s); // index of specified sample
   Nrhc = skexp.CurrPredictedErecoVsQ2Spectrum(s)->Integral();

   cout << "Events: " << Nfhc << " => ";
}

void GetN_cc1pi(PhysHypothesis ph,SuperkExpectation &skexp,double &Nfhc)
{
   // MultiRe # of events 
   conv::SuperkSample_t s = util::SKSampleFromString("MultiRe"); // requested sample
   if(s == conv::kSKsUndefined)
   {
      cerr << "** Undefined sample!" << endl;
      return;
   }
   int is = util::ArrIdx::Instance()->SuperkSampleIdx(s); // index of specified sample
   Nfhc = skexp.CurrPredictedErecoVsQ2Spectrum(s)->Integral();
}

bool OscParamsPhysical(double mix12val, double mix13val, double dm12val)
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
   if(dm12val < kMinFitDm12sq || dm12val > kMaxFitDm12sq)
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
      double dm12val = Rnd::Instance()->RndToyMCosc().Gaus(kPrior_Dm12sq_Used_Mean, kPrior_Dm12sq_Used_RMS);

      mix12val = util::FromSinsq2thToSinsqth(mix12val);
      mix13val = util::FromSinsq2thToSinsqth(mix13val);

      // Leave the atmospheric parameters fixed
      ph.SetCurrValue(conv::kPhypMixing12, mix12val);
      ph.SetCurrValue(conv::kPhypMixing13, mix13val);
      ph.SetCurrValue(conv::kPhypDm12sq, dm12val);
      if(OscParamsPhysical(mix12val, mix13val, dm12val))
      {
         break;
      }
   }
}

