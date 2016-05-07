
void PlotPValues_nuebar
(
 const char *path = "../TECH_NOTE_14_3_16/BETABAR/PVALUES",

 const double Asimov = 6.6078912
)
{
  char name[300];

  double betabar[6]   = {0,0.5, 1., 1.5, 2., 2.5}; //{0., 0.5, 1., 1.5, 2., 2.5}; //, 3.};
////Run1-6
//  double PValue_r[6]  = {0.503334, 0.242265, 0.0943287, 0.094432, 0.0307989, 0.00860991};
//  double PValue_rs[6] = {0.3896, 0.1853, 0.0797, 0.0316, 0.0112, 0.0037};
//  double nevents[6]   = { 1.7186374, 2.4695405, 3.2204437, 3.9713468, 4.72225, 5.4731531};

  //Run1-7
    double PValue_r[6]  = {0.457, 0.143, 0.0671, 0.0110, 0.00391, 0.000392};
    double PValue_rs[6] = {0.4098, 0.1299, 0.0299, 0.005, 0.0014, 0.0003};
    double nevents[6]   = { 3.5263988, 5.067145, 6.6078912, 8.1486373, 9.6893835, 11.23013};

  TGraph *gPvalue01_Vs_betabar = new TGraph();
  TGraph *gPvalueEvt_Vs_betabar = new TGraph();
  TGraph *gPvalue01_Vs_nevents = new TGraph();
  TGraph *gPvalueEvt_Vs_nevents = new TGraph();

  for(int i=0;i< 6;i++)
  {
    gPvalue01_Vs_betabar->SetPoint(i,betabar[i],PValue_rs[i]);
    gPvalue01_Vs_nevents->SetPoint(i,nevents[i],PValue_rs[i]);

    gPvalueEvt_Vs_betabar->SetPoint(i,betabar[i],PValue_r[i]);
    gPvalueEvt_Vs_nevents->SetPoint(i,nevents[i],PValue_r[i]);
  }
  


  TCanvas *cPvalues_Vs_betabar = new TCanvas("cPvalues_Vs_betabar","cPvalues_Vs_betabar");
  //cPvalues_Vs_betabar->SetLogy();
  gPvalue01_Vs_betabar->SetLineColor(kBlue);
  gPvalue01_Vs_betabar->SetMarkerSize(2.);
  gPvalue01_Vs_betabar->SetMarkerStyle(23);
  gPvalue01_Vs_betabar->SetLineWidth(2);
  gPvalue01_Vs_betabar->GetXaxis()->SetTitle("#bar{#beta}");
  gPvalue01_Vs_betabar->GetXaxis()->SetLabelOffset(0.02);
  gPvalue01_Vs_betabar->GetYaxis()->SetTitle("p-value");
  gPvalue01_Vs_betabar->GetYaxis()->SetRangeUser(0.,0.7);
  gPvalue01_Vs_betabar->GetYaxis()->SetLabelOffset(0.02);
  //gPvalueEvt_Vs_betabar->SetLineColor(kViolet);
  //gPvalueEvt_Vs_betabar->SetMarkerSize(2.);
  //gPvalueEvt_Vs_betabar->SetMarkerStyle(24);
  //gPvalueEvt_Vs_betabar->SetLineWidth(2);
  gPvalue01_Vs_betabar->Draw("alp");
  gPvalueEvt_Vs_betabar->SetLineColor(kGreen+2);
  gPvalueEvt_Vs_betabar->SetMarkerSize(2.);
  gPvalueEvt_Vs_betabar->SetMarkerStyle(22);
  gPvalueEvt_Vs_betabar->SetLineWidth(2);
  gPvalueEvt_Vs_betabar->Draw("lp");
  //gPvalueEvt_Vs_betabar->Draw("lp same");
  TLine *lAsimov_betabar = new TLine(1,0.,1,0.6);
  lAsimov_betabar->SetLineColor(kRed);
  lAsimov_betabar->SetLineWidth(2);
  lAsimov_betabar->Draw();
  TLegend *lPvalues_betabar = new TLegend(0.5,0.6,0.9,0.9);
  lPvalues_betabar->SetFillColor(0);
  lPvalues_betabar->SetLineColor(0);
  lPvalues_betabar->SetTextFont(132);
  lPvalues_betabar->SetTextSize(0.06);
  lPvalues_betabar->AddEntry(gPvalue01_Vs_betabar,"Q = #chi^{2}(#bar{#beta}=0) - #chi^{2}(#bar{#beta}=1)","lp");
  lPvalues_betabar->AddEntry(gPvalueEvt_Vs_nevents,"Q = # of events","lp");
  sprintf(name,"%s","PMNS expectation");
  lPvalues_betabar->AddEntry(lAsimov_betabar,name,"l");
  lPvalues_betabar->Draw();




  TCanvas *cPvalues_Vs_nevents = new TCanvas("cPvalues_Vs_nevents","cPvalues_Vs_nevents");
  //cPvalues_Vs_nevents->SetLogy();
  gPvalue01_Vs_nevents->SetLineColor(kBlue);
  gPvalue01_Vs_nevents->SetMarkerSize(2.);
  gPvalue01_Vs_nevents->SetMarkerStyle(22);
  gPvalue01_Vs_nevents->SetLineWidth(2);
  gPvalue01_Vs_nevents->GetXaxis()->SetTitle("# of observed events");
  gPvalue01_Vs_nevents->GetXaxis()->SetLabelOffset(0.02);
  gPvalue01_Vs_nevents->GetYaxis()->SetTitle("p-value");
  gPvalue01_Vs_nevents->GetYaxis()->SetRangeUser(0.,0.7);
  gPvalue01_Vs_nevents->GetYaxis()->SetLabelOffset(0.02);
  //gPvalueEvt_Vs_nevents->SetLineColor(kViolet);
  //gPvalueEvt_Vs_nevents->SetMarkerSize(2.);
  //gPvalueEvt_Vs_nevents->SetMarkerStyle(24);
  //gPvalueEvt_Vs_nevents->SetLineWidth(2);
  gPvalue01_Vs_nevents->Draw("alp");
  gPvalueEvt_Vs_nevents->SetLineColor(kGreen+2);
  gPvalueEvt_Vs_nevents->SetMarkerSize(2.);
  gPvalueEvt_Vs_nevents->SetMarkerStyle(23);
  gPvalueEvt_Vs_nevents->SetLineWidth(2);
  gPvalueEvt_Vs_nevents->Draw("lp");

  //gPvalueEvt_Vs_nevents->Draw("lp same");
  TLine *lAsimov_nevents = new TLine(Asimov,0.,Asimov,0.6);
  lAsimov_nevents->SetLineColor(kRed);
  lAsimov_nevents->SetLineWidth(2);
  lAsimov_nevents->Draw();  
  TLegend *lPvalues_nevents = new TLegend(0.5,0.6,0.9,0.9);
  lPvalues_nevents->SetFillColor(0);
  lPvalues_nevents->SetLineColor(0);
  lPvalues_nevents->SetTextFont(132);
  lPvalues_nevents->SetTextSize(0.06);
  lPvalues_nevents->AddEntry(gPvalue01_Vs_nevents,"Q = #chi^{2}(#bar{#beta}=0) - #chi^{2}(#bar{#beta}=1)","lp");
  lPvalues_nevents->AddEntry(gPvalueEvt_Vs_nevents,"Q = # of events","lp");
  sprintf(name,"%s","PMNS expectation"); 
  lPvalues_nevents->AddEntry(lAsimov_nevents,name,"l");
  lPvalues_nevents->Draw();


}
