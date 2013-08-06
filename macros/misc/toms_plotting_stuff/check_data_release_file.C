{
  //Get the histograms from the originals
  TFile fo1("2012results/NEWRESULTS/datafit/dchisq_2all.root");
  TFile fo2("2012results/DATA/octant2_23only/combination/fitresults_3.root");
  TH2D * ho1 = 0;
  fo1.GetObject("hdchisq",ho1);
  TH2D * ho1_90   = new TH2D(*ho1);
  TH2D * ho1_68   = new TH2D(*ho1);
  TH2D * ho2 = 0;
  fo2.GetObject("hdchisq",ho2);
  TH2D * ho2_90   = new TH2D(*ho2);
  TH2D * ho2_68   = new TH2D(*ho2);
  //set up for drawing
  ho1_90->SetContour(1);
  ho1_90->SetContourLevel(0, 4.61);
  ho1_68->SetContour(1);
  ho1_68->SetContourLevel(0, 2.30);
  ho2_90->SetContour(1);
  ho2_90->SetContourLevel(0, 4.61);
  ho2_68->SetContour(1);
  ho2_68->SetContourLevel(0, 2.30);
  ho1_90->SetLineWidth(1);
  ho1_68->SetLineWidth(1);
  ho2_90->SetLineWidth(1);
  ho2_68->SetLineWidth(1);

  //Get the graphs from the copy
  TFile f("t2k_2013paper_likelihoodratio_contour_and_deltachisq.root");
  g_contour_68_o1->SetLineColor(kRed);
  g_contour_68_o1->SetLineStyle(2);
  g_contour_68_o1->SetLineWidth(1);
  g_contour_90_o1->SetLineColor(kRed);
  g_contour_90_o1->SetLineStyle(2);
  g_contour_90_o1->SetLineWidth(1);
  g_contour_68_o2->SetLineColor(kRed);
  g_contour_68_o2->SetLineStyle(2);
  g_contour_68_o2->SetLineWidth(1);
  g_contour_90_o2->SetLineColor(kRed);
  g_contour_90_o2->SetLineStyle(2);
  g_contour_90_o2->SetLineWidth(1);
  //get the dchisq from the copy
  TH2D * hdchisqo1 = 0;
  f.GetObject("h_deltachisq_o1",hdchisqo1);
  TH2D * hdchisqo2 = 0;
  f.GetObject("h_deltachisq_o2",hdchisqo2);
  //get the bf graph from the copy
  g_bestfit_o1.Print();
  g_bestfit_o2.Print();

  //param space
  TH2D * osc_param_space = new TH2D("osc_param_space", "",
       210, 
       0.9 - 0.5*0.005, 
       1.0 + 0.5*0.005,
       250, 
       0.0021 - 0.5*0.00005, 
       0.0033 + 0.5*0.00005);
  osc_param_space->GetXaxis()->SetTitle("sin^{2}(2#theta_{23})"); // adjust if mixing param is sinsqtheta rather than sinsq2theta...
  osc_param_space->GetYaxis()->SetTitle("|#Deltam_{32}^{2}| (eV^{2}/c^{4})");
  osc_param_space->GetYaxis()->SetTitleOffset(1.4);

  //draw
  TCanvas * c = new TCanvas();
  c->Divide(1,3);
  c->cd(1);
  osc_param_space->Draw();
  ho1_90->Draw("cont3same");
  ho1_68->Draw("cont3same");
  ho2_90->Draw("cont3same");
  ho2_68->Draw("cont3same");
  g_contour_68_o1->Draw();
  g_contour_90_o1->Draw();
  g_contour_68_o2->Draw();
  g_contour_90_o2->Draw();
  c->cd(2);
  osc_param_space->Draw();
  hdchisqo1->Draw("colz");
  c->cd(3);
  osc_param_space->Draw();
  hdchisqo2->Draw("colz");

  c->SaveAs("test_contour_graphs.pdf");



}
