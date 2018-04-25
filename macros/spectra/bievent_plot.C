//
// ------------------------------------------------------------------------------------------------
// Plot contributions from various event modes in the predicted SuperK 1Re or 1Rmu spectrum.
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

#include <iomanip.h>
#include <vector>
#include <algorithm>
#include <functional>

using namespace std;
using namespace valor;
using namespace valor::conv; 

void bievent_plot
(
 const char * tag                  = "ellipse_Run1_9c",
 int          iexperiment         = 0,
 const char * rootfile = "official/2018/data/SK_data_Run19c_valor.root", 
 const bool draw_data = true
)
{
   using namespace std;
   using namespace valor;

   const int dcp_npt = 9; //201;
   const int N23 = 3;
   const int NMH = 2;
   const int NDCP = 4;

   TGraph gN_events[N23][NMH];
   TGraph gN_events_1pt[N23][NMH][NDCP];
   TH2D* hEvents[dcp_npt][N23][NMH];
   TH2D* hEventsSyst[dcp_npt][N23][NMH];
   TH1D* hFHCEvents[dcp_npt][N23][NMH];
   TH1D* hRHCEvents[dcp_npt][N23][NMH];
   bool display_errors[dcp_npt][N23][NMH];
   display_errors[2][1][0] = true;  // Systematics at dcp = -pi/2, ss23 = 0.5, NH

   double deltaCP_first = -TMath::Pi();
   double deltaCP_last = +TMath::Pi();
   double deltaCP_width = (deltaCP_last - deltaCP_first) / (dcp_npt - 1);

   // PDG 2016 best fits
   double ss12[NMH] = {0.304, 0.304};
   double dm12[NMH] = {7.53E-5, 7.53E-5};
   // 2018 best fit values from 2D fits (including IH shift for dmsq23)
   double ss13[NMH] = {0.02235, 0.02235};
   double dm23[NMH] = {2.440E-3, 2.410E-3 - dm12[2]};
   // Alternative values to use for sinsqtheta23
   double ss23[N23] = {0.450, 0.500, 0.550};

   TString rootfile2(Form("%s/inputs/%s", gSystem->Getenv("T2KOSCANA"), rootfile));
   SuperkObservation * skobs_fakedata = new SuperkObservation;

   TGraph *gN_events_data = new TGraph();
   if(draw_data)
   {
      skobs_fakedata->Read(Form("%s/inputs/%s", gSystem->Getenv("T2KOSCANA"),
               rootfile), iexperiment, false);

      // 1Re # of events 
      conv::SuperkSample_t s = util::SKSampleFromString("1Re");
      int is = util::ArrIdx::Instance()->SuperkSampleIdx(s);
      double nfhc = skobs_fakedata->TwoDSpectrum(s)->Integral();
      s = util::SKSampleFromString("MultiRe");
      is = util::ArrIdx::Instance()->SuperkSampleIdx(s);
      nfhc += skobs_fakedata->TwoDSpectrum(s)->Integral();

      // RHC1Re # of events 
      s = util::SKSampleFromString("RHC1Re"); // requested sample
      is = util::ArrIdx::Instance()->SuperkSampleIdx(s); // index of specified sample
      double nrhc = skobs_fakedata->TwoDSpectrum(s)->Integral();
    
      gN_events_data->SetPoint(gN_events_data->GetN(), nfhc, nrhc);
   }

   TFile* file = new TFile("events.root");
   double Nfhc = 0.0;
   double Nrhc = 0.0;

   for(int idcp = 0; idcp < dcp_npt; ++idcp)
   {
      for(int iss23 = 0; iss23 < N23; ++iss23)
      {
         for(int imh = 0; imh < NMH; ++imh)
         {
            double dcp = deltaCP_first + idcp * deltaCP_width;

            hEvents[idcp][iss23][imh] = static_cast<TH2D*>(file->Get(Form("hEvents_%d_%d_%d", idcp, iss23, imh)));
            hEventsSyst[idcp][iss23][imh] = static_cast<TH2D*>(file->Get(Form("hEventsSyst_%d_%d_%d", idcp, iss23, imh)));
            hFHCEvents[idcp][iss23][imh] = static_cast<TH1D*>(file->Get(Form("hEvents_fhc_%d_%d_%d", idcp, iss23, imh)));
            hRHCEvents[idcp][iss23][imh] = static_cast<TH1D*>(file->Get(Form("hEvents_rhc_%d_%d_%d", idcp, iss23, imh)));
            Nfhc = hEvents[idcp][iss23][imh]->GetMean(1);
            Nrhc = hEvents[idcp][iss23][imh]->GetMean(2);
            //cout << idcp << " " << iss23 << " " << imh << " => " << Nfhc << " " << Nrhc << endl;
            gN_events[iss23][imh].SetPoint(gN_events[iss23][imh].GetN(), Nfhc, Nrhc);

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

            if(display_errors[idcp][iss23][imh])
            {
               double contours[1] = {get_contour_level(hEvents[idcp][iss23][imh])};
               hEvents[idcp][iss23][imh]->SetContour(1, contours);
               double syst_contours[1] = {get_contour_level(hEventsSyst[idcp][iss23][imh])};
               hEventsSyst[idcp][iss23][imh]->SetContour(1, syst_contours);
            }
         }
      }
   }

   // Init histogram
  
   TH2D *hInit = new TH2D("hInit", "hInit", 150, 0, 150, 25, 0, 25);
   double xmax = 100.0;
   double xmin = 25.0;
   double ymin = 4.0;
   double ymax = 18.0;
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

   for(int idcp = 0; idcp < dcp_npt; ++idcp)
   {
      for(int iss23 = 0; iss23 < N23; ++iss23)
      {
         for(int imh = 0; imh < NMH; ++imh)
         {
            if(display_errors[idcp][iss23][imh])
            {
//               hEvents[idcp][iss23][imh]->SetFillColorAlpha(kRed, 0.5);
               hEvents[idcp][iss23][imh]->SetLineColor(kRed);
               hEvents[idcp][iss23][imh]->SetLineStyle(3);
               hEvents[idcp][iss23][imh]->Draw("cont3 same");
//               hEventsSyst[idcp][iss23][imh]->SetFillColorAlpha(kBlue, 0.5);
               hEventsSyst[idcp][iss23][imh]->SetLineColor(kBlue);
               hEventsSyst[idcp][iss23][imh]->SetLineStyle(3);
               hEventsSyst[idcp][iss23][imh]->Draw("cont3 same");
            }
         }
      }
   }

   gN_events_data->SetMarkerColor(kBlack);
   gN_events_data->SetMarkerSize(1);
   gN_events_data->SetMarkerStyle(23);
   gN_events_data->Draw("p");

   TLatex* title = new TLatex(ttx, 1.01 * ty, "T2K Run1-9c Preliminary");
   title->SetTextSize(0.05);
   title->Draw();
   TLatex* caveat = new TLatex(ctx, 1.01 * ty, "Final systematics pending");
   caveat->SetTextSize(0.025);
//   caveat->Draw();

   TLegend *ldcp_all = new TLegend(0.15,0.19,0.35,0.58);
   ldcp_all->SetFillColor(0);
   ldcp_all->SetLineColor(0);
   ldcp_all->SetTextFont(132);
   ldcp_all->SetTextSize(0.035);
   ldcp_all->AddEntry(&gN_events[1][0],"sin^{2} #theta_{23} = #color[1]{0.50}, #color[2]{0.45}, #color[4]{0.55}","l");
//   ldcp_all->AddEntry(&gN_events[0][0],"sin^{2} #theta_{23} = 0.45","l");
//   ldcp_all->AddEntry(&gN_events[2][0],"sin^{2} #theta_{23} = 0.55","l");
   ldcp_all->AddEntry(&gN_events[1][0],"#Deltam^{2}_{32} = 2.44#times10^{-3} eV^{2}/c^{4}","l");
   ldcp_all->AddEntry(&gN_events[1][1],"#Deltam^{2}_{13} = 2.41#times10^{-3} eV^{2}/c^{4}","l");
   ldcp_all->AddEntry(&gN_events_1pt[1][0][3],"#delta_{CP} = #pi","p");
   ldcp_all->AddEntry(&gN_events_1pt[1][0][2],"#delta_{CP} = +#pi/2","p");
   ldcp_all->AddEntry(&gN_events_1pt[1][0][1],"#delta_{CP} = 0","p");
   ldcp_all->AddEntry(&gN_events_1pt[1][0][0],"#delta_{CP} = -#pi/2","p");
   for(int idcp = 0; idcp < dcp_npt; ++idcp)
   {
      for(int iss23 = 0; iss23 < N23; ++iss23)
      {
         for(int imh = 0; imh < NMH; ++imh)
         {
            if(display_errors[idcp][iss23][imh])
            {
               ldcp_all->AddEntry(hEventsSyst[idcp][iss23][imh], "syst err","f");
               ldcp_all->AddEntry(hEvents[idcp][iss23][imh], "stat + syst err","f");
            }
         }
      }
   }
   ldcp_all->AddEntry(gN_events_data,"Data","p");
   ldcp_all->Draw();

   cEllipse_all->SaveAs(Form("%s.pdf", tag));
   cEllipse_all->SaveAs(Form("%s.png", tag));
   cEllipse_all->SaveAs(Form("%s.eps", tag));

   TCanvas *c_fhc_contour = new TCanvas("c_contour","c_contour");
   hFHCEvents[2][1][0]->Draw("cont0");
   std::pair<double, double> interval_fhc = get_interval(hFHCEvents[2][1][0]);
   TLine* l_low = new TLine(interval_fhc.first, 0.0, interval_fhc.first, hFHCEvents[2][1][0]->GetMaximum());
   TLine* l_high = new TLine(interval_fhc.second, 0.0, interval_fhc.second, hFHCEvents[2][1][0]->GetMaximum());
   l_low->Draw("same");
   l_high->Draw("same");
   c_fhc_contour->SaveAs("contour_fhc.pdf");

   TCanvas *c_rhc_contour = new TCanvas("c_contour","c_contour");
   hRHCEvents[2][1][0]->Draw("cont0");
   std::pair<double, double> interval_rhc = get_interval(hRHCEvents[2][1][0]);
   TLine* l_low2 = new TLine(interval_rhc.first, 0.0, interval_rhc.first, hRHCEvents[2][1][0]->GetMaximum());
   TLine* l_high2 = new TLine(interval_rhc.second, 0.0, interval_rhc.second, hRHCEvents[2][1][0]->GetMaximum());
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

