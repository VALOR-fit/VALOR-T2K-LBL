//
// ------------------------------------------------------------------------------------------------
// Plot contributions from various event modes in the predicted SuperK 1Re or 1Rmu spectrum.
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

#include <iomanip.h>
#include "valor/Utilities/UtlXML.h"

using namespace std;
using namespace valor;
using namespace valor::conv;

using std::string;
using std::cout;
using std::cerr;
using std::endl;

void bievent_build
(
 const int ntoys = 100,
 const int seed = 1,
 
)

{
   using namespace std;
   using namespace valor;



      // Parse the input XML configuration file
   TXMLEngine* xml = new TXMLEngine;
   XMLDocPointer_t doc = xml->ParseFile("bieventplot_params.xml");
   if(!doc)
   {
      util::Die(1);
   }

   string tag = ReadS(xml, doc, "input_parameters/tag");
   string sk_sample = ReadS(xml, doc, "input_parameters/sk_sample");      
   string hypothesis  = ReadS(xml, doc, "input_parameters/hypothesis");

   // Best fit values NH
   double input_mixing12 = ReadD(xml, doc, "bf_values_nh/input_mixing12");
   double input_mixing13 = ReadD(xml, doc, "bf_values_nh/input_mixing13");
   double input_mixing23 = ReadD(xml, doc, "bf_values_nh/input_mixing23");
   double input_deltacp = ReadD(xml, doc, "bf_values_nh/input_deltacp");
   double input_dm12sq = ReadD(xml, doc, "bf_values_nh/input_dm12sq");
   double input_dm23sq = ReadD(xml, doc, "bf_values_nh/input_dm23sq");
   double input_betabar = ReadD(xml, doc, "bf_values_nh/input_betabar");
   bool disable_banff_tuning = ReadB(xml, doc, "bf_values_nh/disable_banff_tuning");
   string config_file = ReadS(xml, doc, "bf_values_nh/config_file");
   string mixparamchoice = ReadS(xml, doc, "bf_values_nh/mixparamchoice");
   string hierarchy = ReadS(xml, doc, "bf_values_nh/hierarchy");
   string dmsqscheme = ReadS(xml, doc, "bf_values_nh/dmsqscheme");
   string mixparamoctant = ReadS(xml, doc, "bf_values_nh/mixparamoctant");

   double deltaCP_first = -TMath::Pi();
   double deltaCP_last = +TMath::Pi();
   const int dcp_npt = ReadI(xml, doc, "osc_params/dcp_npt"); 
   double deltaCP_width = (deltaCP_last - deltaCP_first) / (dcp_npt-1);
   const int N23 = ReadI(xml, doc, "osc_params/N23");
   const int NMH = ReadI(xml, doc, "osc_params/NMH");
   const int NDCP = ReadI(xml, doc, "osc_params/NDCP");

   // PDG 2016 best fits
   vector<double> ss12 = ReadDA(xml,doc,"vectors/ss12");
   vector<double> dm12 = ReadDA(xml,doc,"vectors/dm12");
   // 2018 best fit values from 2D fits (including IH shift for dmsq23)
   vector<double> ss13 = ReadDA(xml,doc,"vectors/ss13");
   vector<double> dm23 = ReadDA(xml,doc,"vectors/dm23");
   dm23[2] = dm23[2] - dm12[2]; //IH shift for dmsq23
   // Alternative values to use for sinsqtheta23
   vector<double> ss23 = ReadDA(xml,doc,"vectors/ss23");

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

   conv::SuperkSample_t s = util::SKSampleFromString(sk_sample); // requested sample
   if(s == conv::kSKsUndefined)
   {
      cerr << "** Undefined sample!" << endl;
      return;
   }
   int is = util::ArrIdx::Instance()->SuperkSampleIdx(s); // index of specified sample

   // Analysis configuration
   const char * config  = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file.c_str());
   string sconfig(config);

   // Print options
   cout << "\n[Job inputs]" << endl;
   cout << "- Reading analysis configuration from: " << sconfig << endl;
   cout << "- Disabling BANFF tuning: " << ((disable_banff_tuning) ? "YES" : "NO") << endl;
   cout << "- Physics hypothesis: " << ph << endl;
   cout << *Opt::Instance() << endl;

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

   for(int e = 0; e < ntoys; ++e)
   {
      // stat, syst, osc, fitter, other
      if(e % 100 == 0)
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
       
               if(imh == 0)
               {
                  Opt::Instance()->gMassHierarchyAssumed = util::HierarchyFromString("normal");
                  skexp.GetSuperkMC()->SetMassHierarchy(conv::kNmhNormal);
                  ph.SetCurrValue(conv::kPhypDm23sq, dm23[imh]);
               }
               else
               {
                  Opt::Instance()->gMassHierarchyAssumed = util::HierarchyFromString("inverted");
                  skexp.GetSuperkMC()->SetMassHierarchy(conv::kNmhInverted);
                  ph.SetCurrValue(conv::kPhypDm23sq, dm23[imh]);
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

   TFile file(Form("%s_%d.root", tag.c_str(), seed), "RECREATE");
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



bool PathExists(TXMLEngine* xml, XMLDocPointer_t doc, string path)
{
   XMLNodePointer_t root_node = xml->DocGetRootElement(doc);
   return root_node && util::NodeExists(xml, root_node, path);
}

////Read int /////
int ReadI(TXMLEngine* xml, XMLDocPointer_t doc, string node)
{
   if(PathExists(xml, doc, node))
   {
      XMLNodePointer_t root_node = xml->DocGetRootElement(doc);
      return util::XML2Int(xml, root_node, node);
   }
   else
   {
      delete xml;
      util::Die(1);
   }

   return 0;
}

////Read string  /////
string ReadS(TXMLEngine* xml, XMLDocPointer_t doc, string node)
{
   if(PathExists(xml, doc, node))
   {
      XMLNodePointer_t root_node = xml->DocGetRootElement(doc);
      return util::XML2String(xml, root_node, node);
   }
   else
   {
      delete xml;
      util::Die(1);
   }

   return 0;
}

////Read double /////
double ReadD(TXMLEngine* xml, XMLDocPointer_t doc, string node)
{ 
   if(PathExists(xml, doc, node))
   {
      XMLNodePointer_t root_node = xml->DocGetRootElement(doc);
      return util::XML2Dbl(xml, root_node, node);
   }
   else
   {
      delete xml;
      util::Die(1);
   }

   return 0;
}

////Read bool /////
bool ReadB(TXMLEngine* xml, XMLDocPointer_t doc, string node)
{ 
   if(PathExists(xml, doc, node))
   {
      XMLNodePointer_t root_node = xml->DocGetRootElement(doc);
      return util::XML2Bool(xml, root_node, node);
   }
   else
   {
      delete xml;
      util::Die(1);
   }

   return 0;
}

////Read double array /////
vector<double> ReadDA(TXMLEngine* xml, XMLDocPointer_t doc, string node)
{ 
   if(PathExists(xml, doc, node))
   {
      XMLNodePointer_t root_node = xml->DocGetRootElement(doc);
      return util::XML2DblArray(xml, root_node, node);
   }
   else
   {
      delete xml;
      util::Die(1);
   }

   return 0;
}

