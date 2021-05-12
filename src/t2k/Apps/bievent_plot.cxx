//
// ------------------------------------------------------------------------------------------------
// Plot contributions from various event modes in the predicted SuperK 
// 1) e-like neutrino and antineutrino spectra  (--compare-by-nu-antinu)
// 2) e-like neutrino events with Erec below or above 500MeV (--compare-by-ereco)
//
// April 2021 
// Maria Antonova < maria.antonova \at ific.uv.es>
// Costas Andreopoulos
//  

// ------------------------------------------------------------------------------------------------
//

#include <iomanip>
#include <vector>
#include <algorithm>
#include <functional>

//ROOT
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TMarker.h"
#include "TMath.h"
#include "Math/QuantFuncMathCore.h"
#include "TGraphAsymmErrors.h"
#include "TStyle.h" 

//valor 
#include "valor/FitCore/CompoundFitElements.h"
#include "valor/FitCore/FitModelABC.h"
#include "valor/FitCore/FitDataABC.h"
#include "valor/CmdLine/CmdLnArgList.h"
#include "valor/Messenger/Messenger.h"

using namespace std;
using namespace valor;
using namespace valor::conv; 


using std::string;
using std::cout;
using std::cerr;
using std::endl; 

void bievent_plot(CmdLnArgList& arg);

double get_contour_level(TH2D* histo);
double get_contour_level(TH1D* histo);
std::pair<double, double> get_interval(TH1D* histo);

int main(int argc, char** argv)
{
   gROOT->SetBatch(1);

   CmdLnArgList arg;
   arg.Read(argc,argv,
      "$VALOR_EXPT/config/apps/bievent_plot.arg");
   LOG("valor_dune", pNOTICE) << arg;
   bievent_plot(arg);
   LOG("valor_dune", pNOTICE) << "Done!";
  
   return 0;
}

void bievent_plot
(CmdLnArgList& arg)
{
   // creating new colours 
   int blue = TColor::GetFreeColorIndex();
   TColor *bi= new TColor(blue, 0.2,0.64,0.80);
   int green = TColor::GetFreeColorIndex();
   TColor *gi= new TColor(green, 0,0.61,0.23);
   int red = TColor::GetFreeColorIndex();
   TColor *ri= new TColor(red, 0.87,0.12,0.18);
   int yellow = TColor::GetFreeColorIndex();
   TColor *yi= new TColor(yellow, 0.98,0.51,0.15);
   //ROOT setup
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);


   Messenger::Instance()->Configure(util::CompleteConfigPath(arg.GetString("mesg-thresholds")));
   // get data, events and nominal files 
   string events_file   = arg.GetString("events-file");
   string nominal_file  = arg.GetString("nominal-file");
   string tag      = arg.GetString("tag");
   bool use_simple_binning = arg.GetBool("use-simple-binning");
   //Reading histograms
   TFile* file = new TFile(events_file.c_str(),"read");
   TFile* filenom = new TFile(nominal_file.c_str(),"read");  
   

   //Read data? 
   const char* read_opt_data = "";
   bool draw_data  = false;
   double bf_dm23sq;  
   double bf_deltacp; 
   double bf_mixing23;
   string bf_hier;    
   if (arg.Exists("draw-data")){
      std::cout<<"draw-data"<<std::endl;
      LOG("bievent", pINFO) << "Will draw data ";
      draw_data = true;
      if (use_simple_binning) {
         read_opt_data = "+SB";
      }

   }
   //Number of the fixed points 
   const int dcp_npt = arg.GetInt("dcp-npt");
   const int N23     = arg.GetInt("N23");
   const int NMH     = arg.GetInt("NMH");
   const int NDCP    = arg.GetInt("NDCP");


   //If plotting stytematic errors contour at BF
   bool display_errors = false;
   if (arg.Exists("display-errors")){
      display_errors = true;
   }
   double deltaCP_first = -TMath::Pi();
   double deltaCP_last = +TMath::Pi();
   double deltaCP_width = (deltaCP_last - deltaCP_first) / (dcp_npt - 1);

   //Define&get graphs and histos 
   TGraph gN_events[N23][NMH];
   TGraph gN_events_1pt[N23][NMH][NDCP];
   TH2D* hEvents;
   TH2D* hEventsSyst;
   TH1D* hFHCEvents;
   TH1D* hRHCEvents;
   double Nfhc = -999.;
   double Nrhc = -999.;
   hEvents = static_cast<TH2D*>(file->Get("hEvents"));
   hEventsSyst = static_cast<TH2D*>(file->Get(Form("hEventsSyst")));
   hFHCEvents = static_cast<TH1D*>(file->Get(Form("hEvents_fhc")));
   hRHCEvents = static_cast<TH1D*>(file->Get(Form("hEvents_rhc")));


   //Draw data? 
   string data_file  = arg.GetString("data-file");
   TGraphAsymmErrors *gN_events_data = new TGraphAsymmErrors();
   if(draw_data && data_file!="-")
   {    

      TFile* dfile = new TFile(data_file.c_str(),"read");

      // FHC # of events 
      double nfhc = static_cast<TGraph*>(dfile->Get("g_kinematics_data_superk_JPARC_FHC_1re"))->GetN();
      nfhc += static_cast<TGraph*>(dfile->Get("g_kinematics_data_superk_JPARC_FHC_multire"))->GetN();

      // RHC1Re # of events 

      double nrhc = static_cast<TGraph*>(dfile->Get("g_kinematics_data_superk_JPARC_RHC_1re"))->GetN();
    
      gN_events_data->SetPoint(gN_events_data->GetN(), nfhc, nrhc);

      // Setting sqrt(N) as the statistical error on a data point is wrong, as the data has no stat error.
      // Instead calculate an error band showing how far away the prediction can be while still keeping
      // the data point within its Poisson stat errors.
      
      const double alpha = 1 - 0.6827;
      double N = nrhc;
      double L =  (N==0) ? 0  : (ROOT::Math::gamma_quantile(alpha/2,N,1.));
      double U =  ROOT::Math::gamma_quantile_c(alpha/2,N+1,1) ;
      gN_events_data->SetPointEYlow(gN_events_data->GetN()-1, N-L);
      gN_events_data->SetPointEYhigh(gN_events_data->GetN()-1, U-N);

      N=nfhc;
      L =  (N==0) ? 0  : (ROOT::Math::gamma_quantile(alpha/2,N,1.));
      U =  ROOT::Math::gamma_quantile_c(alpha/2,N+1,1) ;
      gN_events_data->SetPointEXlow(gN_events_data->GetN()-1, N-L);
      gN_events_data->SetPointEXhigh(gN_events_data->GetN()-1, U-N);
   }

   else if (draw_data && data_file=="-")
   {
      LOG("bievent", pFATAL) << "You want to plot data, but there is no data file...";
       util::Die(-1);
   }


   for(int idcp = 0; idcp < dcp_npt; ++idcp)
   {
      for(int iss23 = 0; iss23 < N23; ++iss23)
      {
         for(int imh = 0; imh < NMH; ++imh)
         {
            double dcp = deltaCP_first + idcp * deltaCP_width;
            // Consider median instead of mean
            TMarker* nom = (TMarker*)((TObjArray*)filenom->Get("Nominals"))->At(idcp * N23 * NMH + iss23 * NMH + imh);
            Nfhc = nom->GetX();
            Nrhc = nom->GetY();

            gN_events[iss23][imh].SetPoint(gN_events[iss23][imh].GetN(), Nfhc, Nrhc);
            // groupinng events to 4 dcp points 
            if(fabs(dcp - (-TMath::Pi() / 2.)) < 1.E-4 ) 
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
   // create error contours 
   if(display_errors)
   {
      double contours[1] = {get_contour_level(hEvents)};
      hEvents->SetContour(1, contours);
      double syst_contours[1] = {get_contour_level(hEventsSyst)};
      hEventsSyst->SetContour(1, syst_contours);
   }
   // Init histogram
  
   TH2D *hInit = new TH2D("hInit", "hInit", 90, 0, 80, 70, 0, 70);
   double xmax = 80.0;//130.0;
   double xmin = 25.0;
   double ymin = 20.0;//4.0;
   double ymax = 70.0;//25.0;
   hInit->GetXaxis()->SetRangeUser(xmin,xmax);
   hInit->GetYaxis()->SetRangeUser(ymin,ymax);
   double ttx = (xmax + xmin) * 0.5;
   double ctx = xmin + (xmax + xmin) * 0.1;
   double ty = ymax;
   hInit->SetLabelFont(132);
   hInit->SetLabelSize(0.04);
   hInit->GetYaxis()->SetLabelSize(0.04);
   hInit->GetYaxis()->SetLabelFont(132);

   hInit->GetXaxis()->SetTitleFont(132);
   hInit->GetYaxis()->SetTitleFont(132);
   hInit->GetXaxis()->SetTitleSize(0.04);
   hInit->GetYaxis()->SetTitleSize(0.04);
    
   hInit->GetXaxis()->SetTitle("Neutrino mode e-like candidates");
   hInit->GetYaxis()->SetTitle("Antineutrino mode e-like candidates");

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
         else if (it == 2)
            {
               gN_events[it][im].SetLineColor(blue);
               gN_events[it][im].SetMarkerColor(blue);
            }
         else
         {
            gN_events[it][im].SetLineColor(it == 3 ? green : red);
            gN_events[it][im].SetMarkerColor(it == 3 ? green : red);
         }
         gN_events[it][im].SetLineStyle(im == 0 ? 1 : 2);
         for(int id = 0; id < NDCP; ++id)
         {
            gN_events_1pt[it][im][id].SetMarkerStyle(marker_styles[id]);
            if(it == 1)
            {
               gN_events_1pt[it][im][id].SetMarkerColor(kBlack);
            }
            else if (it == 2)
            {
               gN_events_1pt[it][im][id].SetMarkerColor(blue);
            }
            else{
               gN_events_1pt[it][im][id].SetMarkerColor(it == 3 ? green : red);
            }
             
         }
      }
   }

   for(int it = 0; it < N23; ++it)
   {
      if (it ==0 || it ==2){
      for(int im = 0; im < NMH; ++im)
      {
      
        gN_events[it][im].Draw("l");
         
         for(int id = 0; id < NDCP; ++id)
         {
           gN_events_1pt[it][im][id].Draw("p");
         }
      }
      }
   }


   TGraph *gN_events_bf = new TGraph();

   if(display_errors)
   {
      // It is possible to plot stat error on MC but for consistency with other simillar plots it was decided to stick to stat error on data
      /*hEvents->SetLineColor(red);
      hEvents->SetLineStyle(3);
      hEvents->Draw("cont3 same");
      */
      //Display syst errors
      hEventsSyst->SetFillColorAlpha(yellow, 0.15);
      hEventsSyst->SetLineColor(yellow);
      hEventsSyst->SetLineStyle(1);
      hEventsSyst->Draw("cont same");
      //hEventsSyst->Draw("cont3 same");
      //draw bf 
      int bf_point = ((TObjArray*)filenom->Get("Nominals"))->GetLast();
      TMarker* nom = (TMarker*)((TObjArray*)filenom->Get("Nominals"))->At(bf_point);
      Nfhc = nom->GetX();
      Nrhc = nom->GetY();
      gN_events_bf->SetPoint(gN_events_bf->GetN(), Nfhc, Nrhc);
      gN_events_bf->SetMarkerColor(yellow);
      gN_events_bf->SetMarkerSize(1);
      gN_events_bf->SetMarkerStyle(22);
     // gN_events_bf->Draw("p");
   }
   if (draw_data)
   {
      gN_events_data->SetMarkerColor(kBlack);
      gN_events_data->SetMarkerSize(1);
      gN_events_data->SetMarkerStyle(23);
      //gN_events_data->Draw("p");
   }
   

   TLatex* title = new TLatex(ttx, 1.01 * ty, "T2K Run1-10 Preliminary");
   title->SetTextSize(0.05);
   title->SetTextFont(132);
   title->Draw();
   TLatex* caveat = new TLatex(ctx, 1.01 * ty, "Final systematics pending");
   caveat->SetTextSize(0.025);
//   caveat->Draw();

   TLegend *ldcp_all = new TLegend(0.15,0.19,0.35,0.48);
   ldcp_all->SetFillColor(0);
   ldcp_all->SetLineColor(0);
   ldcp_all->SetTextFont(132);
   ldcp_all->SetTextSize(0.037);
   ldcp_all->AddEntry(&gN_events[1][0],"sin^{2} #theta_{23} = #color[1181]{0.45}, #color[1]{0.50}, #color[1179]{0.55}, #color[1180]{0.60}","l");
//   ldcp_all->AddEntry(&gN_events[0][0],"sin^{2} #theta_{23} = 0.45","l");
//   ldcp_all->AddEntry(&gN_events[2][0],"sin^{2} #theta_{23} = 0.55","l");
   ldcp_all->AddEntry(&gN_events[1][0],"#Deltam^{2}_{32} = 2.49#times10^{-3} eV^{2}/c^{4}","l");
   ldcp_all->AddEntry(&gN_events[1][1],"#Deltam^{2}_{13} = -2.46#times10^{-3} eV^{2}/c^{4}","l");
   ldcp_all->AddEntry(&gN_events_1pt[1][0][3],"#delta_{CP} = #pi","p");
   ldcp_all->AddEntry(&gN_events_1pt[1][0][2],"#delta_{CP} = +#pi/2","p");
   ldcp_all->AddEntry(&gN_events_1pt[1][0][1],"#delta_{CP} = 0","p");
   ldcp_all->AddEntry(&gN_events_1pt[1][0][0],"#delta_{CP} = -#pi/2","p");

   TLegend *ldcp_pts = new TLegend(0.37,0.19,0.57,0.36);
   ldcp_pts->SetFillColor(0);
   ldcp_pts->SetLineColor(0);
   ldcp_pts->SetTextFont(132);
   ldcp_pts->SetTextSize(0.037);

   if(display_errors)
   {
      ldcp_pts->AddEntry(hEventsSyst, "68% syst. err. at best fit","f");
      //ldcp_pts->AddEntry(hEvents, "stat + syst err","f");
      ldcp_pts->AddEntry(gN_events_bf,"Best fit","p");
   }
   if (draw_data){
    ldcp_pts->AddEntry(gN_events_data,"Data (68% stat. err.)","p");  
   }
   
   //ldcp_pts->Draw();
   //ldcp_all->Draw();

   cEllipse_all->SaveAs(Form("%s.pdf", tag.c_str()));
   cEllipse_all->SaveAs(Form("%s.png", tag.c_str()));
   cEllipse_all->SaveAs(Form("%s.eps", tag.c_str()));

   TCanvas *c_fhc_contour = new TCanvas("c_contour","c_contour");
   hFHCEvents->Draw("cont0");
   std::pair<double, double> interval_fhc = get_interval(hFHCEvents);
   TLine* l_low = new TLine(interval_fhc.first, 0.0, interval_fhc.first, hFHCEvents->GetMaximum());
   TLine* l_high = new TLine(interval_fhc.second, 0.0, interval_fhc.second, hFHCEvents->GetMaximum());
   l_low->Draw("same");
   l_high->Draw("same");
   c_fhc_contour->SaveAs("contour_fhc.pdf");

   TCanvas *c_rhc_contour = new TCanvas("c_contour","c_contour");
   hRHCEvents->Draw("cont0");
   std::pair<double, double> interval_rhc = get_interval(hRHCEvents);
   TLine* l_low2 = new TLine(interval_rhc.first, 0.0, interval_rhc.first, hRHCEvents->GetMaximum());
   TLine* l_high2 = new TLine(interval_rhc.second, 0.0, interval_rhc.second, hRHCEvents->GetMaximum());
   l_low2->Draw("same");
   l_high2->Draw("same");
   c_rhc_contour->SaveAs("contour_rhc.pdf");

   file->Close();

   return;
}

double get_contour_level(TH2D* histo)
{
   int nx = histo->GetNbinsX();
   int ny = histo->GetNbinsY();
   std::vector<double> events;
   double num_events = 0.0;

   for(int x = 1; x <= nx; ++x)
   {
      for(int y = 1; y <= ny; ++y)
      {
         double value = histo->GetBinContent(x, y);
         events.push_back(value);
         num_events += value;
      }
   }
   // Sort descending
   std::sort(events.begin(), events.end());
   double cutoff = num_events * 0.6827 * 0.6827;
   std::cout << "Overall: " << num_events << " " << cutoff << std::endl;
   num_events = 0.0;
   for(int i = events.size() - 1; i >= 0; --i)
   {
      num_events += events[i];

      if(num_events >= cutoff)
      {
         std::cout << "Cut off:" << cutoff << " " << num_events << " " << events[i] << std::endl;
         return events[i];
      }
   }

   return events[events.size() - 1];
}

double get_contour_level(TH1D* histo)
{
   int nx = histo->GetNbinsX();
   std::vector<double> events;
   double num_events = 0.0;

   for(int x = 1; x <= nx; ++x)
   {
      double value = histo->GetBinContent(x);
      events.push_back(value);
      num_events += value;
   }
   // Sort descending
   std::sort(events.begin(), events.end());
   double cutoff = num_events * 0.6827;
   num_events = 0.0;
   for(int i = events.size() - 1; i >= 0; --i)
   {
      num_events += events[i];
      if(num_events >= cutoff)
      {
         return events[i];
      }
   }

   return events[events.size() - 1];
}

std::pair<double, double> get_interval(TH1D* histo)
{
   int nx = histo->GetNbinsX();
   double num_events = 0.0;

   for(int x = 1; x <= nx; ++x)
   {
      double value = histo->GetBinContent(x);
      num_events += value;
   }
   double num_events_low = num_events * 0.15865;
   double num_events_high = num_events * 0.84135;

   num_events = 0.0;
   double low = 0.0;
   double high = 0.0;
   for(int x = 1; x <= nx; ++x)
   {
      double value = histo->GetBinContent(x);
      num_events += value;
      if(low == 0.0 && num_events >= num_events_low)
      {
         low = histo->GetBinLowEdge(x);
      }
      if(high == 0.0 && num_events >= num_events_high)
      {
         high = histo->GetBinLowEdge(x + 1);
      }
   }

   return std::make_pair(low, high);
}
