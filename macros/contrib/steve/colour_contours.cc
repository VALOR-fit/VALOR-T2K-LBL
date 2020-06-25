////////////////////////////////////////////////////////////
//														                            //
// Scrip to draw sensitivity contours onto a coloured     //
// background giving the chi-squared values.              //
//													                          	  //
// Use the output of the macro "fit_results_to_graphs.py" //
//														                            //
////////////////////////////////////////////////////////////

#include <iostream>
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TPad.h"
#include "TString.h"
#include "TF1.h"
#include "TGraph.h"
#include "TStyle.h"

void plot() {

// Read in output file from fit_results_to_graphs.py
TFile *f = new TFile("contours_as_graphs_sbn_numu.root","read");

// Get contours from file
// Check names in .root file
TGraph *contour_90 		     = (TGraph*) f->Get("contour_90pct_log_0");
TGraph *contour_3_sigma      = (TGraph*) f->Get("contour_3sigma_log_0");
TGraph *contour_5_sigma      = (TGraph*) f->Get("contour_5sigma_log_0");

// Get background 
TH2D *background = (TH2D*) f->Get("contour_log_hist");

// Set colour palette for background
gStyle->SetPalette(56);
// Set log colour scale 
gStyle->SetOptLogz();
// Get the number of different colours - pick a large number to get smooth transitions
background->SetContour(1000);

// Set title and axis (z axis is for the colour)
background->SetTitle("#nu_{#mu} disappearance");
// ~background->SetTitle("#nu_e appearance");

// Range should probably be (-3,0) for numu and (-4,0) for nue
background->GetXaxis()->SetRangeUser(-3,0);
background->GetXaxis()->SetTitle("log_{10}(sin^{2}2#theta_{#mu#mu})");
// ~background->GetXaxis()->SetRangeUser(-4,0);
// ~background->GetXaxis()->SetTitle("log_{10}(sin^{2}2#theta_{#mue})");

background->GetYaxis()->SetRangeUser(-2,2);
background->GetYaxis()->SetTitle("log_{10}(#Deltam^{2}_{41})");
background->GetYaxis()->SetTitleOffset(0.8);

// Have to manually set range so that entire graph is filled
background->GetZaxis()->SetRangeUser(1E-9,1E6);

// Edit canvas so axis titles etc. fit on properly
TCanvas *c1 = new TCanvas("c1","c1");
c1->SetRightMargin(0.2);
c1->SetLeftMargin(0.15);
c1->SetBottomMargin(0.15);

// Draw background 
background->Draw("colz");

// Draw contours
contour_90->SetLineColor(kAzure-4);
contour_90->Draw("same");
contour_3_sigma->SetLineStyle(2);
contour_3_sigma->SetLineColor(kAzure-4);
contour_3_sigma->Draw("same");
contour_5_sigma->SetLineColor(kAzure-4);
contour_5_sigma->SetLineStyle(7);  
contour_5_sigma->Draw("same");  

// Add legend
TLegend *leg = new TLegend(0.65,0.65,0.8,0.85);
leg->SetTextSize(0.05);
leg->SetTextColor(kAzure-4);
leg->SetFillStyle(0);
leg->SetBorderSize(0);
leg->AddEntry(contour_90,"90%","L");
leg->AddEntry(contour_3_sigma,"3#sigma","l");
leg->AddEntry(contour_5_sigma,"5#sigma","l");
leg->Draw("same");

// Add label to colour legend
TLatex *latex = new TLatex();
latex->DrawLatexNDC(0.92,0.88,"#Delta#chi^{2}");




}
