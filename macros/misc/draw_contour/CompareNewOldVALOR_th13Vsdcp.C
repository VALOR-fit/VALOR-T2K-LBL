

// Produce 1D chisq Vs dcp --> plot_contours_dcp.C (both NH and IH)
// Produce 2D chisq Vs dcp Vs sinsqth13 --> draw_contour.C 

using namespace std;
using namespace valor;
using namespace valor::conv;

void CompareNewOldVALOR
( 
 
 char inname_NH_old[200]="", // both 1D and 2D
 char inname_IH_old[200]="", // both 1D and 2D
 
 char inname_2D_NH_new[200]="CONTOURS/fitmarg_sinsqth13Vsdcp_T2Kreact_DataRun14_NH_contour.root",
 char inname_2D_IH_new[200]="CONTOURS/fitmarg_sinsqth13Vsdcp_T2Kreact_DataRun14_IH_contour.root",
 
 char inname_1D_new[200]="CONTOURS/fitmarg_dcp_T2Kreact_DataRun14_contour.root", // contains both NH and IH
 
  ) 
{
  char name[200];
  
  int color_new = kBlack;
  int color_old = kRed;

  // New T2K fit (1D)
  TFile *filein_1D_new = new TFile(inname_1D_new,"READ");
  TH1D *h_deltacp_NH_new = (TH1D*)filein_1D_new->Get("hdeltachisq_NH");
  TH1D *h_deltacp_IH_new = (TH1D*)filein_1D_new->Get("hdeltachisq_IH");
  
  // New T2K fit (2D - NH)
  TFile *filein_2D_NH_new = new TFile(inname_2D_NH_new,"READ");
  TGraph *g90_valor_0_NH_new = (TGraph*)filein_2D_NH_new->Get("g90_valor_0");
  TGraph *g90_valor_1_NH_new = (TGraph*)filein_2D_NH_new->Get("g90_valor_1");
  TGraph *g68_valor_0_NH_new = (TGraph*)filein_2D_NH_new->Get("g68_valor_0");
  TGraph *g68_valor_1_NH_new = (TGraph*)filein_2D_NH_new->Get("g68_valor_1");
  TH2D *chisq_map_2D_NH = (TH2D*)filein_2D_NH_new ->Get("chisq_map_2D");
  chisq_map_2D_NH->Reset();

  // New T2K fit (2D - IH)
  TFile *filein_2D_IH_new = new TFile(inname_2D_IH_new,"READ");
  TGraph *g90_valor_0_IH_new = (TGraph*)filein_2D_IH_new->Get("g90_valor_0");
  TGraph *g90_valor_1_IH_new = (TGraph*)filein_2D_IH_new->Get("g90_valor_1");
  TGraph *g68_valor_0_IH_new = (TGraph*)filein_2D_IH_new->Get("g68_valor_0");
  TGraph *g68_valor_1_IH_new = (TGraph*)filein_2D_IH_new->Get("g68_valor_1");
  TH2D *chisq_map_2D_IH = (TH2D*)filein_2D_IH_new ->Get("chisq_map_2D");
  chisq_map_2D_IH->Reset();
 
  // Old T2K fit (NH)
  TFile *filein_NH_old = new TFile(inname_NH_old,"READ");
  TH2D *h_theta13_deltacp_NH_old = (TH2D*)filein_NH_old->Get("h_theta13_deltacp");
  h_theta13_deltacp_NH_old->SetContour(1);
  h_theta13_deltacp_NH_old->SetContourLevel(0,0.);
  TH2D *h_theta13_deltacp_68_NH_old = (TH2D*)filein_NH_old->Get("h_theta13_deltacp_68");
  h_theta13_deltacp_68_NH_old->SetContour(1);
  h_theta13_deltacp_68_NH_old->SetContourLevel(0,0.);
  TH1D *h_deltacp_NH_old = (TH1D*)filein_NH_old->Get("h_deltacp");

  // Old T2K fit (IH)
  TFile *filein_IH_old = new TFile(inname_IH_old,"READ");
  TH2D *h_theta13_deltacp_IH_old = (TH2D*)filein_IH_old->Get("h_theta13_deltacp");
  h_theta13_deltacp_IH_old->SetContour(1);
  h_theta13_deltacp_IH_old->SetContourLevel(0,0.);
  TH2D *h_theta13_deltacp_68_IH_old = (TH2D*)filein_IH_old->Get("h_theta13_deltacp_68");
  h_theta13_deltacp_68_IH_old->SetContour(1);
  h_theta13_deltacp_68_IH_old->SetContourLevel(0,0.);
  TH1D *h_deltacp_IH_old = (TH1D*)filein_IH_old->Get("h_deltacp");
  

  // From sinsqth13 to sinsq2th13

  ConvertHistoFromSinsqthToSinsq2th(g90_valor_0_NH_new);
  ConvertHistoFromSinsqthToSinsq2th(g90_valor_1_NH_new);
  ConvertHistoFromSinsqthToSinsq2th(g68_valor_0_NH_new);
  ConvertHistoFromSinsqthToSinsq2th(g68_valor_1_NH_new);
  ConvertHistoFromSinsqthToSinsq2th(g90_valor_0_IH_new);
  ConvertHistoFromSinsqthToSinsq2th(g90_valor_1_IH_new);
  ConvertHistoFromSinsqthToSinsq2th(g68_valor_0_IH_new);
  ConvertHistoFromSinsqthToSinsq2th(g68_valor_1_IH_new);

 
 // Draw contours 1D

  TCanvas *ccontours_1d = new TCanvas("ccontours_1d","ccontours_1d");
  h_deltacp_IH_old->SetLineWidth(2);
  h_deltacp_NH_old->SetLineWidth(2);
  h_deltacp_IH_old->SetLineColor(color_old);
  h_deltacp_NH_old->SetLineColor(color_old);
  h_deltacp_IH_old->SetLineStyle(2);

  h_deltacp_IH_new->SetLineWidth(2);
  h_deltacp_NH_new->SetLineWidth(2);
  h_deltacp_IH_new->SetLineColor(color_new);
  h_deltacp_NH_new->SetLineColor(color_new);
  h_deltacp_IH_new->SetLineStyle(2);

  h_deltacp_IH_old->GetXaxis()->SetTitle("#delta_{CP}");
  h_deltacp_IH_old->GetYaxis()->SetTitle("#Delta #chi^{2}");

  h_deltacp_IH_old->Draw();
  h_deltacp_NH_old->Draw("same");
  h_deltacp_IH_new->Draw("same");
  h_deltacp_NH_new->Draw("same");
  h_deltacp_IH_old->GetYaxis()->SetRangeUser(0.,7.);

  TLegend *lcontours_1d = new TLegend(0.5,0.6,0.8,0.8);
  lcontours_1d->SetFillColor(0);
  lcontours_1d->SetLineColor(0);
  lcontours_1d->SetTextFont(132);
  lcontours_1d->SetTextSize(0.05);
  lcontours_1d->AddEntry(h_deltacp_NH_new,"New NH","l");
  lcontours_1d->AddEntry(h_deltacp_IH_new,"New IH","l");
  lcontours_1d->AddEntry(h_deltacp_NH_old,"Old NH","l");
  lcontours_1d->AddEntry(h_deltacp_IH_old,"Old IH","l");
  lcontours_1d->Draw();
  
  TLine *l1sigma = new TLine(h_deltacp_IH_old->GetXaxis()->GetBinLowEdge(1),1.,h_deltacp_IH_old->GetXaxis()->GetBinUpEdge(h_deltacp_IH_old->GetNbinsX()),1.);
  TLine *l90 = new TLine(h_deltacp_IH_old->GetXaxis()->GetBinLowEdge(1),2.71,h_deltacp_IH_old->GetXaxis()->GetBinUpEdge(h_deltacp_IH_old->GetNbinsX()),2.71);
  TLine *l2sigma = new TLine(h_deltacp_IH_old->GetXaxis()->GetBinLowEdge(1),4.,h_deltacp_IH_old->GetXaxis()->GetBinUpEdge(h_deltacp_IH_old->GetNbinsX()),4.);
  l1sigma->SetLineWidth(2);
  l90->SetLineWidth(2);
  l2sigma->SetLineWidth(2);
  l1sigma->SetLineColor(kBlue);
  l90->SetLineColor(kBlue);
  l2sigma->SetLineColor(kBlue);
  l1sigma->Draw();
  l90->Draw();
  l2sigma->Draw();
  
  TLatex *latLine1Sigma = new TLatex(-1.7,1.1,"1#sigma");
  TLatex *latLine90CL = new TLatex(-1.7,2.81,"90% CL");
  TLatex *latLine2Sigma = new TLatex(-1.7,4.1,"2#sigma");
  latLine1Sigma->SetTextColor(kBlue);
  latLine90CL->SetTextColor(kBlue);
  latLine2Sigma->SetTextColor(kBlue);
  latLine1Sigma->SetTextSize(0.07);
  latLine90CL->SetTextSize(0.07);
  latLine2Sigma->SetTextSize(0.07);
  latLine1Sigma->Draw("same");
  latLine90CL->Draw("same");
  latLine2Sigma->Draw("same");


  // Draw contours 2D - NH

  TCanvas *ccontours_NH = new TCanvas("ccontours_NH","ccontours_NH");

  h_theta13_deltacp_NH_old->SetLineColor(color_old);
  h_theta13_deltacp_68_NH_old->SetLineColor(color_old);
  g90_valor_0_NH_new->SetLineColor(color_new);
  g90_valor_1_NH_new->SetLineColor(color_new);
  g68_valor_0_NH_new->SetLineColor(color_new);
  g68_valor_1_NH_new->SetLineColor(color_new);

  h_theta13_deltacp_NH_old->GetXaxis()->SetTitle("sin^{2} 2 #theta_{13}");
  h_theta13_deltacp_NH_old->GetYaxis()->SetTitle("#delta_{CP}");

  h_theta13_deltacp_NH_old->Draw("cont3");
  h_theta13_deltacp_68_NH_old->Draw("cont3 same");
  g90_valor_1_NH_new->Draw("l same");
  g90_valor_0_NH_new->Draw("l same");
  g68_valor_0_NH_new->Draw("l same");
  g68_valor_1_NH_new->Draw("l same");
  ccontours_NH->Update();

  TLegend *lcontours_NH = new TLegend(0.5,0.6,0.8,0.8);
  lcontours_NH->SetFillColor(0);
  lcontours_NH->SetLineColor(0);
  lcontours_NH->SetTextFont(132);
  lcontours_NH->SetTextSize(0.05);
  lcontours_NH->AddEntry(g68_valor_1_NH_new,"New 68% CL","l");
  lcontours_NH->AddEntry(g90_valor_1_NH_new,"New 90% CL","l");
  lcontours_NH->AddEntry(h_theta13_deltacp_68_NH_old,"Old 68% CL","l");
  lcontours_NH->AddEntry(h_theta13_deltacp_NH_old,"Old 90% CL","l");
  lcontours_NH->Draw();
  
  TLatex *latcontours_NH = new TLatex(0.065,1.,"NH");
  latcontours_NH->SetTextSize(0.07);
  latcontours_NH->Draw("same");


  // Draw contours 2D - IH

  TCanvas *ccontours_IH = new TCanvas("ccontours_IH","ccontours_IH");

  h_theta13_deltacp_IH_old->SetLineColor(color_old);
  h_theta13_deltacp_68_IH_old->SetLineColor(color_old);
  g90_valor_0_IH_new->SetLineColor(color_new);
  g90_valor_1_IH_new->SetLineColor(color_new);
  g68_valor_0_IH_new->SetLineColor(color_new);
  g68_valor_1_IH_new->SetLineColor(color_new);

  h_theta13_deltacp_IH_old->GetXaxis()->SetTitle("sin^{2} 2 #theta_{13}");
  h_theta13_deltacp_IH_old->GetYaxis()->SetTitle("#delta_{CP}");

  h_theta13_deltacp_IH_old->Draw("cont3");
  h_theta13_deltacp_68_IH_old->Draw("cont3 same");
  g90_valor_0_IH_new->Draw("l same");
  g90_valor_1_IH_new->Draw("l same");
  g68_valor_0_IH_new->Draw("l same");
  g68_valor_1_IH_new->Draw("l same");
  ccontours_IH->Update();

  TLegend *lcontours_IH = new TLegend(0.5,0.6,0.8,0.8);
  lcontours_IH->SetFillColor(0);
  lcontours_IH->SetLineColor(0);
  lcontours_IH->SetTextFont(132);
  lcontours_IH->SetTextSize(0.05);
  lcontours_IH->AddEntry(g68_valor_1_IH_new,"New 68% CL","l");
  lcontours_IH->AddEntry(g90_valor_1_IH_new,"New 90% CL","l");
  lcontours_IH->AddEntry(h_theta13_deltacp_68_IH_old,"Old 68% CL","l");
  lcontours_IH->AddEntry(h_theta13_deltacp_IH_old,"Old 90% CL","l");
  lcontours_IH->Draw();
  
  TLatex *latcontours_IH = new TLatex(0.065,1.,"IH");
  latcontours_IH->SetTextSize(0.07);
  latcontours_IH->Draw("same");

  return;
}



void ConvertHistoFromSinsqthToSinsq2th(TGraph *graph){
  
  for(int ipt=0;ipt<graph->GetN();ipt++){
    double x = 0.;
    double y = 0.;
    graph->GetPoint(ipt,x,y);
    double sinsq2th = util::FromSinsqthToSinsq2th(x);
    graph->RemovePoint(ipt);
    graph->SetPoint(ipt,sinsq2th,y);
  }
  return;
}
