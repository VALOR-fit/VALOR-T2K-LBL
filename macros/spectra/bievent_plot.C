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

void bievent_plot
(
 bool include_cc1pi                = false,
 const char * tag                  = "ellipse_Run18",
 const string str_preliminary = "T2K Run 1-8 preliminary",
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
 const char * config_file = "analysis_configuration_oa2017_run1-8_postfit.xml",
 bool         store_spectrum      = true,
 bool         print_spectrum      = true,                          /* */
 const char * mixparamchoice      = "sinsqtheta",
 const char * hierarchy           = "normal"   ,
 const char * dmsqscheme          = "auto", /* auto==kNuDmsqAuto abs23==kNuDmsq12andAbs23 abs31==kNuDmsq12andAbs31 fogli==kNuDmsqFogliLisi */
 const char * mixparamoctant      = "two",       /* one==kNuOctantOne two==kNuOctantTwo */
 bool plot_per_50MeV = false,
 double ErecoRangeLow = 0.,
 double ErecoRangeUp = 30.,
 int          iexperiment         = 0,
 const char * rootfile = "data_run1to8/SKdata_Run18.root", 
 const bool draw_data_spectra = true,
 const bool scale_mc = false,
 const bool change_banff_params = false
)
{
   using namespace std;
   using namespace valor;

   char name[300];

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
   double dcp_npt = 9; //201;
   double deltaCP_width = (deltaCP_last - deltaCP_first) / (dcp_npt-1);
   const int N23 = 3;
   const int NMH = 2;
   const int NDCP = 4;

   TGraph gN_events[N23][NMH];
   TGraph gN_events_1pt[N23][NMH][NDCP];

   double ss12[NMH] = {0.304, 0.304};
   double ss13[NMH] = {0.02235, 0.0226125};
   double dm12[NMH] = {7.53E-5, 7.53E-5};
   double dm23[NMH] = {2.460E-3, 2.440E-3};
   double ss23[N23] = {0.450, 0.500, 0.550};

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
            ph.SetCurrValue(conv::kPhypMixing12, ss12[imh]);
            ph.SetCurrValue(conv::kPhypMixing13, ss13[imh]);
            ph.SetCurrValue(conv::kPhypDm12sq, dm12[imh]);
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
            gN_events[iss23][imh].SetPoint(gN_events[iss23][imh].GetN(), Nfhc, Nrhc);

            if(fabs(dcp- (-TMath::Pi() / 2.)) < 1.E-4 ) 
            {  // -pi/2
               gN_events_1pt[iss23][imh][0].SetPoint(gN_events_1pt[iss23][imh][0].GetN(), Nfhc, Nrhc);
            }
            else if(fabs(dcp) < 1.E-4) 
            {  // 0
               gN_events_1pt[iss23][imh][1].SetPoint(gN_events_1pt[iss23][imh][1].GetN(), Nfhc, Nrhc);
            }
            else if(fabs(dcp - TMath::Pi() / 2) < 1.E-4 ) 
            {  // +pi/2
               gN_events_1pt[iss23][imh][2].SetPoint(gN_events_1pt[iss23][imh][2].GetN(), Nfhc, Nrhc);
            }
            else if(fabs(dcp - (-TMath::Pi())) < 1.E-4 )
            {  // +pi
               gN_events_1pt[iss23][imh][3].SetPoint(gN_events_1pt[iss23][imh][3].GetN(), Nfhc, Nrhc);
            }
         }
      }
   }

   // Get data spectra and rescale MC to data

   TGraphErrors *gN_events_data = new TGraphErrors();
   if(draw_data)
   {
      bool using_xml_file = false;
      skobs_fakedata->Read(rootfile2.Data(),iexperiment,using_xml_file);

      // 1Re # of events 
      s = util::SKSampleFromString("1Re"); // requested sample
      if(s == conv::kSKsUndefined)
      {
         cerr << "** Undefined sample!" << endl;
         return;
      }
      is = util::ArrIdx::Instance()->SuperkSampleIdx(s); // index of specified sample
      double nfhc = skobs_fakedata->TwoDSpectrum(s)->Integral();
      s = util::SKSampleFromString("MultiRe"); // requested sample
      if(s == conv::kSKsUndefined)
      {
         cerr << "** Undefined sample!" << endl;
         return;
      }
      is = util::ArrIdx::Instance()->SuperkSampleIdx(s); // index of specified sample
      nfhc += skobs_fakedata->TwoDSpectrum(s)->Integral();

      // RHC1Re # of events 
      s = util::SKSampleFromString("RHC1Re"); // requested sample
      if(s == conv::kSKsUndefined)
      {
         cerr << "** Undefined sample!" << endl;
         return;
      }
      is = util::ArrIdx::Instance()->SuperkSampleIdx(s); // index of specified sample
      double nrhc = skobs_fakedata->TwoDSpectrum(s)->Integral();
    
      double nfhc_err = sqrt(nfhc);
      double nrhc_err = sqrt(nrhc);

      cout << endl;
      cout << "Data: nfhc = " << nfhc << " +- " << nfhc_err << ", nrhc = " << nrhc << " +- " << nrhc_err << endl;
      cout << endl;

      int ipt = gN_events_data->GetN();
      gN_events_data->SetPoint(ipt, nfhc, nrhc);
      gN_events_data->SetPointError(ipt, nfhc_err, nrhc_err);

      cout << "idcp = " << idcp << "" << endl;
   }

   // Init histogram
  
   TH2D *hInit = new TH2D("hInit", "hInit", 500, 0, 100, 500, 0, 100);
   double xmax = 100.0;
   double xmin = 25.0;
   double ymin = 4.0;
   double ymax = 12.0;
   hInit->GetXaxis()->SetRangeUser(xmin,xmax);
   hInit->GetYaxis()->SetRangeUser(ymin,ymax);
   double ttx = (xmax + xmin) * 0.5;
   double ctx = xmin + (xmax + xmin) * 0.1;
   double ty = ymax;
   hInit->GetXaxis()->SetTitle("Neutrino mode 1Re candidates");
   hInit->GetYaxis()->SetTitle("Antineutrino mode 1Re candidates");

   // Draw 
   TCanvas *cEllipse_all = new TCanvas("cEllipse_all","cEllipse_all");
   hInit->Draw();

   int marker_styles[NDCP] = {20, 25, 21, 24};
   for(int it = 0; it < N23; ++it)
   {
      for(int im = 0; im < NMH; ++im)
      {
         if(it == 1)
         {
            gN_events[it][im].SetLineColor(kBlack);
            gN_events[it][im].SetMarkerColor(kBlack);
         }
         else
         {
            gN_events[it][im].SetLineColor(it == 2 ? kBlue : kRed);
            gN_events[it][im].SetMarkerColor(it == 2 ? kBlue : kRed);
         }
         gN_events[it][im].SetLineStyle(im == 0 ? 1 : 2);
         for(int id = 0; id < NDCP; ++id)
         {
            gN_events_1pt[it][im][id].SetMarkerStyle(marker_styles[id]);
            if(it == 1)
            {
               gN_events_1pt[it][im][id].SetMarkerColor(kBlack);
            }
            else
            {
               gN_events_1pt[it][im][id].SetMarkerColor(it == 2 ? kBlue : kRed);
            }
         }
      }
   }

   for(int it = 0; it < N23; ++it)
   {
      for(int im = 0; im < NMH; ++im)
      {
         gN_events[it][im].Draw("l");
         for(int id = 0; id < NDCP; ++id)
         {
            gN_events_1pt[it][im][id].Draw("p");
         }
      }
   }

   gN_events_data->SetMarkerColor(kBlack);
   gN_events_data->SetMarkerSize(1);
   gN_events_data->SetMarkerStyle(4);
   gN_events_data->Draw("PE");

   TLatex* title = new TLatex(ttx, 1.01 * ty, "T2K Run1-8 Preliminary");
   title->SetTextSize(0.05);
   title->Draw();
   TLatex* caveat = new TLatex(ctx, 1.01 * ty, "Final systematics pending");
   caveat->SetTextSize(0.025);
   caveat->Draw();

   TLegend *ldcp_all = new TLegend(0.15,0.19,0.35,0.64);
   ldcp_all->SetFillColor(0);
   ldcp_all->SetLineColor(0);
   ldcp_all->SetTextFont(132);
   ldcp_all->SetTextSize(0.035);
   ldcp_all->AddEntry(&gN_events[1][0],"sin^{2} #theta_{23} = 0.50","l");
   ldcp_all->AddEntry(&gN_events[0][0],"sin^{2} #theta_{23} = 0.45","l");
   ldcp_all->AddEntry(&gN_events[2][0],"sin^{2} #theta_{23} = 0.55","l");
   ldcp_all->AddEntry(&gN_events[1][0],"#Deltam^{2}_{32} = 2.46#times10^{-3} eV^{2}/c^{4}","l");
   ldcp_all->AddEntry(&gN_events[1][1],"#Deltam^{2}_{31} = -2.44#times10^{-3} eV^{2}/c^{4}","l");
   ldcp_all->AddEntry(&gN_events_1pt[1][0][3],"#delta_{CP} = #pi","p");
   ldcp_all->AddEntry(&gN_events_1pt[1][0][2],"#delta_{CP} = +#pi/2","p");
   ldcp_all->AddEntry(&gN_events_1pt[1][0][1],"#delta_{CP} = 0","p");
   ldcp_all->AddEntry(&gN_events_1pt[1][0][0],"#delta_{CP} = -#pi/2","p");
   ldcp_all->AddEntry(gN_events_data,"Data (stat. errors only)","pe");
   ldcp_all->Draw();

   cEllipse_all->SaveAs(Form("%s.pdf", tag));
   cEllipse_all->SaveAs(Form("%s.png", tag));
   cEllipse_all->SaveAs(Form("%s.eps", tag));

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
   hErecoVsQ2Tot = new TH2D(*skexp.CurrPredictedErecoVsQ2Spectrum(s));  
   Nfhc = hErecoVsQ2Tot->Integral();
  
   // RHC1Re # of events 
   s = util::SKSampleFromString("RHC1Re"); // requested sample
   if(s == conv::kSKsUndefined)
   {
      cerr << "** Undefined sample!" << endl;
      return;
   }
   is = util::ArrIdx::Instance()->SuperkSampleIdx(s); // index of specified sample
   hErecoVsQ2Tot = new TH2D(*skexp.CurrPredictedErecoVsQ2Spectrum(s));  
   Nrhc = hErecoVsQ2Tot->Integral();

   return;
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
   hErecoVsQ2Tot = new TH2D(*skexp.CurrPredictedErecoVsQ2Spectrum(s));  
   Nfhc = hErecoVsQ2Tot->Integral();
  
   return;
}

