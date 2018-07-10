#include <iomanip.h>
#include <vector>
#include <algorithm>
#include <functional>
#include "valor/Utilities/UtlXML.h"
#include "SetT2KStyle.h"

using namespace valor;
using namespace valor::conv;
using namespace valor::util; 


using std::string;
using std::cout;
using std::cerr;
using std::endl; 

void Plot_dcp_graph_new
()
{
   // Parse the input XML configuration file
   TXMLEngine* xml = new TXMLEngine;
   XMLDocPointer_t doc = xml->ParseFile("bieventplot_params.xml");
   if(!doc)
   {
      util::Die(1);
   }
   XMLNodePointer_t root_node = xml->DocGetRootElement(doc);
 	
   string tag = XML2String(xml, root_node, "plot_data/tag");
   string sk_sample = XML2String(xml, root_node, "input_parameters/sk_sample"); 
   int iexperiment = XML2Int(xml,root_node,"plot_data/iexperiment");
   string rootfile = XML2String(xml, root_node, "plot_data/rootfile");
   bool draw_data = XML2Bool(xml, root_node, "nevents_dcp/draw_data"); 
   string sk_sample = XML2String(xml, root_node, "input_parameters/sk_sample"); 

   const int dcp_npt = ReadI(xml, doc, "nevents_dcp/dcp_npt"); 
   const int N23 =ReadI(xml, doc, "nevents_dcp/N23");
   const int NMH =ReadI(xml, doc, "osc_params/NMH");

   TH1D* hFHCEvents[N23][NMH];
   TH1D* hRHCEvents[N23][NMH];
   TGraphErrors* gFHCEvents[N23][NMH];
   TGraphErrors* gRHCEvents[N23][NMH];

   TH1D* hFHCEventsSystLow[N23][NMH];
   TH1D* hRHCEventsSystLow[N23][NMH];

   TH1D* hFHCEventsSystUp[N23][NMH];
   TH1D* hRHCEventsSystUp[N23][NMH];


   TH1D* hFHCEventsSystLoop[N23][NMH];
   TH1D* hRHCEventsSystLoop[N23][NMH];


 // PDG 2016 best fits
   vector<double> ss12 = XML2DblArray(xml,root_node,"vectors/ss12");
   vector<double> dm12 = XML2DblArray(xml,root_node,"vectors/dm12");
// 2018 best fit values from 2D fits (including IH shift for dmsq23)
   vector<double> ss13 = XML2DblArray(xml,root_node,"vectors/ss13");
   vector<double> dm23 = XML2DblArray(xml,root_node,"vectors/dm23");
   dm23[2] = dm23[2] - dm12[2]; //IH shift for dmsq23
// Alternative values to use for sinsqtheta23
   vector<double> ss23 = XML2DblArray(xml,root_node,"nevents_dcp/ss23");

   TFile* file = new TFile("dcp_syst_grapherror_1_1Re_100_toys.root");
   double Nfhc = XML2Dbl(xml,root_node,"plot_data/Nfhc");
   double Nrhc = XML2Dbl(xml,root_node,"plot_data/Nrhc");
   // Read graphs to plot it 
 	for(int iss23 = 0; iss23 < N23; ++iss23)
   	{
      for(int imh = 0; imh < NMH; ++imh)
      	{
         	gFHCEvents[iss23][imh] = static_cast<TGraphErrors*>(file->Get(Form("gEvents_fhc_%d_%d", iss23, imh)));
         	gRHCEvents[iss23][imh] = static_cast<TGraphErrors*>(file->Get(Form("gEvents_rhc_%d_%d", iss23, imh)));
		}
    }
// ----------------------------------------------------Init lines and styles 
    // Set T2K official style 
    //TStyle* t2kstyle = SetT2KStyle(1,"t2kstyle"); 
    //gROOT->SetStyle(t2kstyle->GetName());
    gStyle->SetOptTitle(1);
   for(int iss23 = 0; iss23 < N23; ++iss23)
   {      
     for(int imh = 0; imh < NMH; ++imh)
      {      	
        gFHCEvents[iss23][imh]->SetLineWidth(4);
        gRHCEvents[iss23][imh]->SetLineWidth(4);
		gFHCEvents[iss23][imh]->SetFillStyle(1001);
		gRHCEvents[iss23][imh]->SetFillStyle(1001);
         
      		//solid lines for sin^{2} #theta_{23} = 0.45 (Solid lines )
        	if (iss23 == 0)
         	{	  
            	gFHCEvents[iss23][imh]->SetLineStyle(1);
            	gRHCEvents[iss23][imh]->SetLineStyle(1);

            	gFHCEvents[iss23][0]->SetLineColor(kBlue);
    			gRHCEvents[iss23][0]->SetLineColor(kBlue);
    			gRHCEvents[iss23][0]->SetFillColorAlpha(kBlue,0.05);
    			gFHCEvents[iss23][0]->SetFillColorAlpha(kBlue,0.05);

    			gFHCEvents[iss23][1]->SetLineColor(kRed);
            	gRHCEvents[iss23][1]->SetLineColor(kRed);
                gRHCEvents[iss23][1]->SetFillColorAlpha(kRed,0.05);
                gFHCEvents[iss23][1]->SetFillColorAlpha(kRed,0.05);

         	}
         	//dash lines for sin^{2} #theta_{23} = 0.55 (dash lines )
        	else
         	{ 
            	gFHCEvents[iss23][imh]->SetLineStyle(9);
            	gRHCEvents[iss23][imh]->SetLineStyle(9); 
            	
            	gFHCEvents[iss23][1]->SetLineColor(kSpring+10);
    			gRHCEvents[iss23][1]->SetLineColor(kSpring+10);
    			gRHCEvents[iss23][1]->SetFillColorAlpha(kSpring+10,0.05);
    			gFHCEvents[iss23][1]->SetFillColorAlpha(kSpring+10,0.05);



    			gFHCEvents[iss23][0]->SetLineColor(kOrange+1);
            	gRHCEvents[iss23][0]->SetLineColor(kOrange+1);
            	gRHCEvents[iss23][0]->SetFillColorAlpha(kOrange+1,0.05);
                gFHCEvents[iss23][0]->SetFillColorAlpha(kOrange+1,0.05);
 
         	}
      	}
   	} 
   	


// Plotting 
   

  // -----------------------------------------------------------------------Draw  FHC NH
  TCanvas *fhc_nh = new TCanvas("cFHC_nh","cFHC_nh");  
  TMultiGraph * gfhc_nh  = new TMultiGraph();
  fhc_nh->SetCanvasSize(3000, 2000);
   double ymax = 140.0;
   double ymin = 40.0;
   double xmin = -3.;
   double xmax =  3.;

   gfhc_nh->Add(gFHCEvents[0][0],"c4");
   gfhc_nh->Add(gFHCEvents[1][0],"c4");
   gfhc_nh->Draw("a");

   gfhc_nh->GetXaxis()->SetLimits(xmin,xmax);
   gfhc_nh->GetYaxis()->SetLimits(ymin,ymax);
   gfhc_nh->GetXaxis()->SetTitle("#delta_{CP} (rad)");
   gfhc_nh->GetYaxis()->SetTitle("Number of events");

   gPad->Modified();
   gPad->Update();

   TLegend *ldcp_fhc = new TLegend(0.20,0.19,0.40,0.48);
   ldcp_fhc->SetFillColor(0);
   ldcp_fhc->SetLineColor(0);
   ldcp_fhc->SetTextFont(132);
   ldcp_fhc->SetTextSize(0.05);
   ldcp_fhc->AddEntry(gFHCEvents[0][0],"sin^{2} #theta_{23} = 0.45","l");
   ldcp_fhc->AddEntry(gFHCEvents[1][0],"sin^{2} #theta_{23} = 0.55","l");
   ldcp_fhc->Draw();

   TLatex* title = new TLatex();
   title->SetTextSize(0.05);
   title->DrawLatex(1.1, 100, Form("%s sample, NH, FHC",sk_sample.c_str()));
   
   gPad->Modified();
   gPad->Update();
   fhc_nh->SaveAs(Form("FHC_graph_nh_%s.pdf", tag.c_str()));
   fhc_nh->SaveAs(Form("FHC_graph_nh_%s.png", tag.c_str()));
   fhc_nh->SaveAs(Form("FHC_graph_nh_%s.eps", tag.c_str()));
  // -----------------------------------------------------------------------Draw  FHC IH
  TCanvas *fhc_ih = new TCanvas("cFHC_ih","cFHC_ih");  
  TMultiGraph * gfhc_ih  = new TMultiGraph();
  fhc_ih->SetCanvasSize(3000, 2000);
   double ymax = 140.0;
   double ymin = 40.0;
   double xmin = -3.;
   double xmax =  3.;

   gfhc_ih->Add(gFHCEvents[0][1],"c4");
   gfhc_ih->Add(gFHCEvents[1][1],"c4");
  
   gfhc_ih->Draw("a");
   gfhc_ih->GetXaxis()->SetLimits(xmin,xmax);
   gfhc_ih->GetYaxis()->SetLimits(ymin,ymax);
   gfhc_ih->GetXaxis()->SetTitle("#delta_{CP} (rad)");
   gfhc_ih->GetYaxis()->SetTitle("Number of events");
   gPad->Modified();
   gPad->Update();

   TLegend *ldcp_fhc_i = new TLegend(0.20,0.19,0.40,0.48);
   ldcp_fhc_i->SetFillColor(0);
   ldcp_fhc_i->SetLineColor(0);
   ldcp_fhc_i->SetTextFont(132);
   ldcp_fhc_i->SetTextSize(0.05);
   ldcp_fhc_i->AddEntry(gFHCEvents[0][1],"sin^{2} #theta_{23} = 0.45","l");
   ldcp_fhc_i->AddEntry(gFHCEvents[1][1],"sin^{2} #theta_{23} = 0.55","l");
   ldcp_fhc_i->Draw();

   TLatex* title = new TLatex();
   title->SetTextSize(0.05);
   title->DrawLatex(1.1, 89.5, Form("%s sample, IH, FHC",sk_sample.c_str()));
   
   gPad->Modified();
   gPad->Update();
   fhc_ih->SaveAs(Form("FHC_graph_ih_%s.pdf", tag.c_str()));
   fhc_ih->SaveAs(Form("FHC_graph_ih_%s.png", tag.c_str()));
   fhc_ih->SaveAs(Form("FHC_graph_ih_%s.eps", tag.c_str()));

 // -----------------------------------------------------------------------Draw  FHC both
  TCanvas *fhc = new TCanvas("cFHC","cFHC");  
  TMultiGraph * gfhc  = new TMultiGraph();
  fhc->SetCanvasSize(3000, 2000);
   double ymax = 140.0;
   double ymin = 40.0;
   double xmin = -3.;
   double xmax =  3.;
	for(int iss23 = 0; iss23 < N23; ++iss23)
   	{      
       for(int imh = 0; imh < NMH; ++imh)
      	{
   			gfhc->Add(gFHCEvents[iss23][imh],"c4");
   			gfhc->Add(gFHCEvents[iss23][imh],"c4");
   		}
	}   
  
     gfhc->Draw("a");

   gfhc->GetXaxis()->SetLimits(xmin,xmax);
   gfhc->GetYaxis()->SetLimits(ymin,ymax);
   gfhc->GetXaxis()->SetTitle("#delta_{CP} (rad)");
   gfhc->GetYaxis()->SetTitle("Number of events");
   gPad->Modified();
   gPad->Update();

   TLegend *ldcp_fhc = new TLegend(0.20,0.19,0.40,0.38);
   ldcp_fhc->SetFillColor(0);
   ldcp_fhc->SetLineColor(0);
   ldcp_fhc->SetTextFont(132);
   ldcp_fhc->SetTextSize(0.05);
   ldcp_fhc->AddEntry(gFHCEvents[0][1],"sin^{2} #theta_{23} = 0.45, IH","l");
   ldcp_fhc->AddEntry(gFHCEvents[1][1],"sin^{2} #theta_{23} = 0.55, IH","l");
   ldcp_fhc->AddEntry(gFHCEvents[0][0],"sin^{2} #theta_{23} = 0.45, NH","l");
   ldcp_fhc->AddEntry(gFHCEvents[1][0],"sin^{2} #theta_{23} = 0.55, NH","l");
   ldcp_fhc->Draw();

   TLatex* title = new TLatex();
   title->SetTextSize(0.05);
   title->DrawLatex(1.1, 100, Form("%s sample, FHC",sk_sample.c_str()));
   
   gPad->Modified();
   gPad->Update();
   fhc->SaveAs(Form("FHC_graph_both_%s.pdf", tag.c_str()));
   fhc->SaveAs(Form("FHC_graph_both_%s.png", tag.c_str()));
   fhc->SaveAs(Form("FHC_graph_both_%s.eps", tag.c_str()));




  // -----------------------------------------------------------------------Draw  RHC NH
  TCanvas *rhc_nh = new TCanvas("cRHC_nh","cRHC_nh");  
  TMultiGraph * grhc_nh  = new TMultiGraph();
  rhc_nh->SetCanvasSize(3000, 2000);
   double ymax = 140.0;
   double ymin = 40.0;
   double xmin = -3.;
   double xmax =  3.;

   grhc_nh->Add(gRHCEvents[0][0],"c4");
   grhc_nh->Add(gRHCEvents[1][0],"c4");
   grhc_nh->Draw("a");
   
   grhc_nh->GetXaxis()->SetLimits(xmin,xmax);
   grhc_nh->GetYaxis()->SetLimits(ymin,ymax);
   grhc_nh->GetXaxis()->SetTitle("#delta_{CP} (rad)");
   grhc_nh->GetYaxis()->SetTitle("Number of events");

   gPad->Modified();
   gPad->Update();

   TLegend *ldcp_rhc = new TLegend(0.60,0.19,0.75,0.29);
   ldcp_rhc->SetFillColor(0);
   ldcp_rhc->SetLineColor(0);
   ldcp_rhc->SetTextFont(132);
   ldcp_rhc->SetTextSize(0.05);
   ldcp_rhc->AddEntry(gRHCEvents[0][0],"sin^{2} #theta_{23} = 0.45","l");
   ldcp_rhc->AddEntry(gRHCEvents[1][0],"sin^{2} #theta_{23} = 0.55","l");
   ldcp_rhc->Draw();

   TLatex* title = new TLatex();
   title->SetTextSize(0.05);
   title->DrawLatex(1., 29., Form("%s sample, NH, RHC",sk_sample.c_str()));
   
   gPad->Modified();
   gPad->Update();
   rhc_nh->SaveAs(Form("RHC_graph_nh_%s.pdf", tag.c_str()));
   rhc_nh->SaveAs(Form("RHC_graph_nh_%s.png", tag.c_str()));
   rhc_nh->SaveAs(Form("RHC_graph_nh_%s.eps", tag.c_str()));
  // -----------------------------------------------------------------------Draw  RHC IH
  TCanvas *rhc_ih = new TCanvas("cRHC_ih","cRHC_ih");  
  TMultiGraph * grhc_ih  = new TMultiGraph();
    rhc_ih->SetCanvasSize(3000, 2000);

   double ymax = 140.0;
   double ymin = 40.0;
   double xmin = -3.;
   double xmax =  3.;

   grhc_ih->Add(gRHCEvents[0][1],"c4");
   grhc_ih->Add(gRHCEvents[1][1],"c4");
   grhc_ih->Draw("a");
   
   grhc_ih->GetXaxis()->SetLimits(xmin,xmax);
   grhc_ih->GetYaxis()->SetLimits(ymin,ymax);
   grhc_ih->GetXaxis()->SetTitle("#delta_{CP} (rad)");
   grhc_ih->GetYaxis()->SetTitle("Number of events");

   gPad->Modified();
   gPad->Update();

   TLegend *ldcp_rhc_i = new TLegend(0.60,0.19,0.75,0.29);
   ldcp_rhc_i->SetFillColor(0);
   ldcp_rhc_i->SetLineColor(0);
   ldcp_rhc_i->SetTextFont(132);
   ldcp_rhc_i->SetTextSize(0.05);
   ldcp_rhc_i->AddEntry(gRHCEvents[0][1],"sin^{2} #theta_{23} = 0.45","l");
   ldcp_rhc_i->AddEntry(gRHCEvents[1][1],"sin^{2} #theta_{23} = 0.55","l");
   ldcp_rhc_i->Draw();

   TLatex* title = new TLatex();
   title->SetTextSize(0.05);
   title->DrawLatex(1., 31.3, Form("%s sample, IH, RHC",sk_sample.c_str()));
   
   gPad->Modified();
   gPad->Update();
   rhc_ih->SaveAs(Form("RHC_graph_ih_%s.pdf", tag.c_str()));
   rhc_ih->SaveAs(Form("RHC_graph_ih_%s.png", tag.c_str()));
   rhc_ih->SaveAs(Form("RHC_graph_ih_%s.eps", tag.c_str()));

// -----------------------------------------------------------------------Draw  RHC both
  TCanvas *rhc = new TCanvas("cRHC","cRHC");  
  TMultiGraph * grhc  = new TMultiGraph();
  rhc->SetCanvasSize(3000, 2000);
   double ymax = 140.0;
   double ymin = 40.0;
   double xmin = -3.;
   double xmax =  3.;
	for(int iss23 = 0; iss23 < N23; ++iss23)
   	{      
       for(int imh = 0; imh < NMH; ++imh)
      	{
   			grhc->Add(gRHCEvents[iss23][imh],"c4");
   			grhc->Add(gRHCEvents[iss23][imh],"c4");
   		}
	}   
  
   grhc->Draw("a");
   grhc->GetXaxis()->SetLimits(xmin,xmax);
   grhc->GetYaxis()->SetLimits(ymin,ymax);
   grhc->GetXaxis()->SetTitle("#delta_{CP} (rad)");
   grhc->GetYaxis()->SetTitle("Number of events");
   gPad->Modified();
   gPad->Update();

   TLegend *ldcp_rhc = new TLegend(0.60,0.19,0.75,0.315);
   ldcp_rhc->SetFillColor(0);
   ldcp_rhc->SetLineColor(0);
   ldcp_rhc->SetTextFont(132);
   ldcp_rhc->SetTextSize(0.04);
   ldcp_rhc->AddEntry(gRHCEvents[0][1],"sin^{2} #theta_{23} = 0.45, IH","l");
   ldcp_rhc->AddEntry(gRHCEvents[1][1],"sin^{2} #theta_{23} = 0.55, IH","l");
   ldcp_rhc->AddEntry(gRHCEvents[0][0],"sin^{2} #theta_{23} = 0.45, NH","l");
   ldcp_rhc->AddEntry(gRHCEvents[1][0],"sin^{2} #theta_{23} = 0.55, NH","l");
   ldcp_rhc->Draw();

   TLatex* title = new TLatex();
   title->SetTextSize(0.05);
   title->DrawLatex(1.1, 31.3, Form("%s sample, RHC",sk_sample.c_str()));
   
   gPad->Modified();
   gPad->Update();
   rhc->SaveAs(Form("RHC_graph_both_%s.pdf", tag.c_str()));
   rhc->SaveAs(Form("RHC_graph_both_%s.png", tag.c_str()));
   rhc->SaveAs(Form("RHC_graph_both_%s.eps", tag.c_str()));

   return;
} 

////Read int /////
int ReadI(TXMLEngine* xml, XMLDocPointer_t doc, string node)
{
   
      XMLNodePointer_t root_node = xml->DocGetRootElement(doc);
      return util::XML2Int(xml, root_node, node);


   return 0;
}
 