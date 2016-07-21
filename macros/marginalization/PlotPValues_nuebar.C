
void PlotPValues_nuebar
(
 const char *path = "../TECH_NOTE_14_3_16/BETABAR/PVALUES"

)
{
  char name[300];

  double betabar[6]   = {0,0.5, 1., 1.5, 2., 2.5}; //{0., 0.5, 1., 1.5, 2., 2.5}; //, 3.};
//Asimov A
  const int nentries = 6;
//  //Run1-6
    const double Asimov = 3.2204437;

  double PValue_r[6]  = {0.503334,      0.242265,     0.0943287,      0.094432,     0.0307989,    0.00860991};
  double PValue_rs[6] = {0.4158348944, 0.1935544041, 0.08906850643, 0.03436473163, 0.01119689291 , 0.003550691895};
  double nevents[6]   = { 1.7186374, 2.4695405, 3.2204437, 3.9713468, 4.72225, 5.4731531};

//  Run1-7
//  const double Asimov = 6.6078912;
//  double PValue_r[6]  = {0.457, 0.143, 0.0671, 0.0110, 0.00391, 0.000392};
//  double PValue_rs[6] = {0.4355774213, 0.1657780591, 0.05991703725, 0.01054610751, 0.003537013992, 0.0005050986184};
//  double nevents[6]   = { 3.5263988, 5.067145, 6.6078912, 8.1486373, 9.6893835, 11.23013};
//
//  //Asimov B
//  const int nentries = 5;
//  ////Run1-6
//  const double Asimov = 3.35687;
//    double PValue_r[6]  = {0.477893, 0.221406, 0.0827643, 0.0258649, 0.00695958, 0.00163353};
//    double PValue_rs[6] = {0.4605138501, 0.1939057722, 0.07692645339, 0.02575690263, 0.006895707933, 0.0008816364336};
//    double nevents[6]   = { 1.5724182, 2.4646441, 3.35687, 4.2490959, 5.1413218, 6.0335476};

  //Run1-7
//    const double Asimov = 6.8878185;
//      double PValue_r[6]  = {0.420833, 0.122312, 0.0552466, 0.00833692, 0.000892502, 7.22319e-05};
//      double PValue_rs[6] = {0.4550416614, 0.1230357867, 0.01976548254, 0.002040463165, 0, 0};
//      double nevents[6]   = { 3.2263779, 5.0570982, 6.8878185, 8.7185387, 10.549259, 12.379979};

  
  TGraph *gPvalue01_Vs_betabar = new TGraph();
  TGraph *gPvalueEvt_Vs_betabar = new TGraph();
  TGraph *gPvalue01_Vs_nevents = new TGraph();
  TGraph *gPvalueEvt_Vs_nevents = new TGraph();

  for(int i=0;i< nentries;i++)
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
  gPvalue01_Vs_betabar->GetXaxis()->SetTitle("#beta");
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
  lPvalues_betabar->AddEntry(gPvalue01_Vs_betabar,"Q = #chi^{2}(#beta=0) - #chi^{2}(#beta=1)","lp");
  lPvalues_betabar->AddEntry(gPvalueEvt_Vs_betabar,"Q = # of events","lp");
  sprintf(name,"%s","PMNS expectation");
  lPvalues_betabar->AddEntry(lAsimov_betabar,name,"l");
  lPvalues_betabar->Draw();




  TCanvas *cPvalues_Vs_nevents = new TCanvas("cPvalues_Vs_nevents","cPvalues_Vs_nevents");
  //cPvalues_Vs_nevents->SetLogy();
  gPvalue01_Vs_nevents->SetLineColor(kBlue);
  gPvalue01_Vs_nevents->SetMarkerSize(2.);
  gPvalue01_Vs_nevents->SetMarkerStyle(23);
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
  gPvalueEvt_Vs_nevents->SetMarkerStyle(22);
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
  lPvalues_nevents->AddEntry(gPvalue01_Vs_nevents,"Q = #chi^{2}(#beta=0) - #chi^{2}(#beta=1)","lp");
  lPvalues_nevents->AddEntry(gPvalueEvt_Vs_nevents,"Q = # of events","lp");
  sprintf(name,"%s","PMNS expectation"); 
  lPvalues_nevents->AddEntry(lAsimov_nevents,name,"l");
  lPvalues_nevents->Draw();


}
