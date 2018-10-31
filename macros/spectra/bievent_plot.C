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
#include "valor/Utilities/UtlXML.h"

using namespace std;
using namespace valor;
using namespace valor::conv; 


using std::string;
using std::cout;
using std::cerr;
using std::endl; 

void bievent_plot
()
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
   
   string tag = ReadS(xml, doc, "plot_data/tag");
   int iexperiment = ReadI(xml,doc,"plot_data/iexperiment");
   string rootfile = ReadS(xml, doc, "plot_data/rootfile");
   bool draw_data = ReadB(xml, doc, "plot_data/draw_data");  

   const int dcp_npt = ReadI(xml, doc, "osc_params/dcp_npt");
   const int N23 = ReadI(xml, doc, "osc_params/N23");
   const int NMH = ReadI(xml, doc, "osc_params/NMH");
   const int NDCP = ReadI(xml, doc, "osc_params/NDCP");

   TGraph gN_events[N23][NMH];
   TGraph gN_events_1pt[N23][NMH][NDCP];
   TH2D* hEvents[dcp_npt][N23][NMH];
   TH2D* hEventsSyst[dcp_npt][N23][NMH];
   TH1D* hFHCEvents[dcp_npt][N23][NMH];
   TH1D* hRHCEvents[dcp_npt][N23][NMH];
   bool display_errors[dcp_npt][N23][NMH];
   int a = ReadI(xml,doc,"plot_data/display_errors/first_param");
   int b = ReadI(xml,doc,"plot_data/display_errors/second_param");
   int c = ReadI(xml,doc,"plot_data/display_errors/third_param");

   for(int idcp = 0; idcp < dcp_npt; ++idcp)
   {
      for(int iss23 = 0; iss23 < N23; ++iss23)
      {
         for(int imh = 0; imh < NMH; ++imh)
         {
            display_errors[idcp][iss23][imh] = false;
         }
      }
   }

   display_errors[a][b][c] = ReadB(xml,doc,"plot_data/display_errors/bool_value");// Systematics at dcp = -pi/2, ss23 = 0.5, NH

   double deltaCP_first = -TMath::Pi();
   double deltaCP_last = +TMath::Pi();
   double deltaCP_width = (deltaCP_last - deltaCP_first) / (dcp_npt - 1);

 // PDG 2016 best fits
   vector<double> ss12 = ReadDA(xml,doc,"vectors/ss12");
   vector<double> dm12 = ReadDA(xml,doc,"vectors/dm12");
// 2018 best fit values from 2D fits (including IH shift for dmsq23)
   vector<double> ss13 = ReadDA(xml,doc,"vectors/ss13");
   vector<double> dm23 = ReadDA(xml,doc,"vectors/dm23");
   dm23[2] = dm23[2] - dm12[2]; //IH shift for dmsq23
// Alternative values to use for sinsqtheta23
   vector<double> ss23 = ReadDA(xml,doc,"vectors/ss23");

   SuperkObservation * skobs_fakedata = new SuperkObservation;

   TGraph *gN_events_data = new TGraph();
   if(draw_data)
   {
      skobs_fakedata->Read(Form("%s/inputs/%s", gSystem->Getenv("T2KOSCANA"),
               rootfile.c_str()), iexperiment, false);

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
   TFile* filenom = new TFile("nominals.root");
   double Nfhc = ReadD(xml,doc,"plot_data/Nfhc");
   double Nrhc = ReadD(xml,doc,"plot_data/Nrhc");

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

            // Consider median instead of mean
            TMarker* nom = (TMarker*)((TObjArray*)filenom->Get("Nominals"))->At(idcp * N23 * NMH + iss23 * NMH + imh);
            Nfhc = nom->GetX();
            Nrhc = nom->GetY();
//            Nfhc = hEventsSyst[idcp][iss23][imh]->GetMean(1);
//            Nrhc = hEventsSyst[idcp][iss23][imh]->GetMean(2);
/*            cout << "Syst" << endl;
            cout << idcp << " " << iss23 << " " << imh << " => " <<
               Nfhc << "(" << hEventsSyst[idcp][iss23][imh]->GetRMS(1) << ") " <<
               Nrhc << "(" << hEventsSyst[idcp][iss23][imh]->GetRMS(2) << ") " << endl;
            cout << "Stat+Syst" << endl;
            cout << idcp << " " << iss23 << " " << imh << " => " <<
               hEvents[idcp][iss23][imh]->GetMean(1) << "(" << hEvents[idcp][iss23][imh]->GetRMS(1) << ") " <<
               hEvents[idcp][iss23][imh]->GetMean(2) << "(" << hEvents[idcp][iss23][imh]->GetRMS(2) << ") " << endl;
*/
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
   double ymax = 25.0;
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
/*
               TEllipse eventsEllipse(hEvents[idcp][iss23][imh]->GetMean(1),
                     hEvents[idcp][iss23][imh]->GetMean(2),
                     hEvents[idcp][iss23][imh]->GetRMS(1),
                     hEvents[idcp][iss23][imh]->GetRMS(2));
               eventsEllipse.SetLineStyle(2);
               eventsEllipse.SetLineColor(kRed);
               eventsEllipse.SetFillStyle(0);
               eventsEllipse.Draw("same");

               TEllipse eventsSystEllipse(hEventsSyst[idcp][iss23][imh]->GetMean(1),
                     hEventsSyst[idcp][iss23][imh]->GetMean(2),
                     hEventsSyst[idcp][iss23][imh]->GetRMS(1),
                     hEventsSyst[idcp][iss23][imh]->GetRMS(2));
               eventsSystEllipse.SetLineStyle(2);
               eventsSystEllipse.SetLineColor(kBlue);
               eventsSystEllipse.SetFillStyle(0);
               eventsSystEllipse.Draw("same");
*/
//               hEventsSyst[idcp][iss23][imh]->Draw("TEXT same");
            }
         }
      }
   }

   gN_events_data->SetMarkerColor(kBlack);
   gN_events_data->SetMarkerSize(1);
   gN_events_data->SetMarkerStyle(23);
   gN_events_data->Draw("p");

   TLatex* title = new TLatex(ttx, 1.01 * ty, "T2K Run1-9 Preliminary");
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

   cEllipse_all->SaveAs(Form("%s.pdf", tag.c_str()));
   cEllipse_all->SaveAs(Form("%s.png", tag.c_str()));
   cEllipse_all->SaveAs(Form("%s.eps", tag.c_str()));

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
