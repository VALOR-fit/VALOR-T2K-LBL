TGraph * get_contour(TH2D * h) {
  TCanvas * c = new TCanvas();
  h->Draw("cont list");
  c->Update();
  TObjArray * contours  = gROOT->GetListOfSpecials()->FindObject("contours");
  if(contours->GetSize() != 1)
    exit(1);
  TList     * list      = (TList*)contours->At(0);
  if(list->GetSize() != 1)
    exit(1);
  TGraph * graph = (TGraph*)list->First();
  delete c;
  return graph;
}

void make_data_release_file()
{
  TFile fo1("2012results/NEWRESULTS/datafit/dchisq_2all.root");
  TFile fo2("2012results/DATA/octant2_23only/combination/fitresults_3.root");

  double yoffset = 1.4; //y axis title offset

  //limits for osc param space
  double xmin = 0.7975,   xmax = 1.0025,   nxbins = 41;
  double ymin = 0.000975, ymax = 0.004025, nybins = 61;

  //make an osc param space
  TH2D * h_osc_param_space = new TH2D("h_osc_param_space", "Paramater space to facilitate drawing", nxbins, xmin, xmax, nybins, ymin, ymax);
  h_osc_param_space->GetXaxis()->SetTitle("sin^{2}(2#theta_{23})");
  h_osc_param_space->GetYaxis()->SetTitle("|#Deltam_{32}^{2}| (eV^{2}/c^{4})");
  h_osc_param_space->GetYaxis()->SetTitleOffset(yoffset);

  //limits for chisq surfaces
  xmin = 0.8975,   xmax = 1.0025,   nxbins = 21;
  ymin = 0.001975, ymax = 0.003025, nybins = 21;

  //get the octant one contour
  TH2D * ho1 = 0;
  fo1.GetObject("hdchisq",ho1);
  //need to make sure the range of the 2 octants is the same
  TH2D * ho1_cut = new TH2D("", "", nxbins, xmin, xmax, nybins, ymin, ymax);
  for(int ix = 1; ix <= nxbins; ix++) {
    for(int iy = 1; iy <= nxbins; iy++) {
      double x = ho1_cut->GetXaxis()->GetBinCenter(ix);
      double y = ho1_cut->GetYaxis()->GetBinCenter(iy);
      int bin  = ho1->FindBin(x, y);
      double content = ho1->GetBinContent(bin);
      ho1_cut->SetBinContent(ix, iy, content);
    }//iy
  }//ix
  ho1_cut->SetLineWidth(2);
  //copy
  TH2D * ho1_colz = new TH2D(*ho1_cut);
  TH2D * ho1_90   = new TH2D(*ho1_cut);
  TH2D * ho1_68   = new TH2D(*ho1_cut);
  //set up for drawing
  ho1_colz->SetContour(255);
  ho1_90->SetContour(1);
  ho1_90->SetContourLevel(0, 4.61);
  ho1_68->SetContour(1);
  ho1_68->SetContourLevel(0, 2.30);
  //get the contours
  TGraph * go1_90 = new TGraph(*get_contour(ho1_90));
  go1_90->SetTitle("90%CL contour for #theta_{23} #leq #pi/4, NH (2013 numu paper)");
  go1_90->SetName("g_contour_90_o1");
  TGraph * go1_68 = new TGraph(*get_contour(ho1_68));
  go1_68->SetTitle("68%CL contour for #theta_{23} #leq #pi/4, NH (2013 numu paper)");
  go1_68->SetName("g_contour_68_o1");
  //set axis titles
  ho1_colz->GetXaxis()->SetTitle("sin^{2}(2#theta_{23})");
  ho1_colz->GetYaxis()->SetTitle("|#Deltam_{32}^{2}| (eV^{2}/c^{4})");
  go1_90->  GetXaxis()->SetTitle("sin^{2}(2#theta_{23})");
  go1_90->  GetYaxis()->SetTitle("|#Deltam_{32}^{2}| (eV^{2}/c^{4})");
  go1_68->  GetXaxis()->SetTitle("sin^{2}(2#theta_{23})");
  go1_68->  GetYaxis()->SetTitle("|#Deltam_{32}^{2}| (eV^{2}/c^{4})");
  ho1_colz->GetYaxis()->SetTitleOffset(yoffset);
  go1_90->  GetYaxis()->SetTitleOffset(yoffset);
  go1_68->  GetYaxis()->SetTitleOffset(yoffset);
  //set the dchisq surface name
  ho1_colz->SetTitle("#Delta#chi^{2} surface for #theta_{23} #leq #pi/4, NH (2013 numu paper)");
  ho1_colz->SetName("h_deltachisq_o1");
  //best fit point as graph
  double bfo1_x[1] = {1.000};
  double bfo1_y[1] = {0.00244};
  TGraph * bfo1 = new TGraph(1, bfo1_x, bfo1_y);
  bfo1->SetTitle("best fit point for #theta_{23} #leq #pi/4, NH (2013 numu paper)");
  bfo1->SetName("g_bestfit_o1");
  bfo1->SetMarkerSize(1.8);
  bfo1->SetMarkerColor(kBlack);
  bfo1->SetMarkerStyle(29);
  bfo1->GetXaxis()->SetTitle("sin^{2}(2#theta_{23})");
  bfo1->GetYaxis()->SetTitle("|#Deltam_{32}^{2}| (eV^{2}/c^{4})");
  bfo1->GetYaxis()->SetTitleOffset(yoffset);
  //best fit point as marker
  TMarker * mo1 = new TMarker(1.000,0.00244,29);
  mo1->SetMarkerSize(1.8);
  mo1->SetMarkerColor(kBlack);
  cout << mo1->GetTitle() << endl;
  //mo1->SetTitle("best fit point for #theta_{23} #leq #pi/4, NH (2013 numu paper)");
  //mo1->SetName("m_bestfit_o1");

  //get the octant two contour
  TH2D * ho2 = 0;
  fo2.GetObject("hdchisq",ho2);
  //need to make sure the range of the 2 octants is the same
  TH2D * ho2_cut = new TH2D("", "", nxbins, xmin, xmax, nybins, ymin, ymax);
  for(int ix = 1; ix <= nxbins; ix++) {
    for(int iy = 1; iy <= nxbins; iy++) {
      double x = ho2_cut->GetXaxis()->GetBinCenter(ix);
      double y = ho2_cut->GetYaxis()->GetBinCenter(iy);
      int bin  = ho2->FindBin(x, y);
      double content = ho2->GetBinContent(bin);
      ho2_cut->SetBinContent(ix, iy, content);
    }//iy
  }//ix
  ho2_cut->SetLineWidth(2);
  //copy
  TH2D * ho2_colz = new TH2D(*ho2_cut);
  TH2D * ho2_90   = new TH2D(*ho2_cut);
  TH2D * ho2_68   = new TH2D(*ho2_cut);
  //set up for drawing
  ho2_colz->SetContour(255);
  ho2_90->SetContour(1);
  ho2_90->SetContourLevel(0, 4.61);
  ho2_68->SetContour(1);
  ho2_68->SetContourLevel(0, 2.30);
  //get the contours
  TGraph * go2_90 = new TGraph(*get_contour(ho2_90));
  go2_90->SetTitle("90%CL contour for #theta_{23} #geq #pi/4, NH (2013 numu paper)");
  go2_90->SetName("g_contour_90_o2");
  TGraph * go2_68 = new TGraph(*get_contour(ho2_68));
  go2_68->SetTitle("68%CL contour for #theta_{23} #geq #pi/4, NH (2013 numu paper)");
  go2_68->SetName("g_contour_68_o2");
  //set axis titles
  ho2_colz->GetXaxis()->SetTitle("sin^{2}(2#theta_{23})");
  ho2_colz->GetYaxis()->SetTitle("|#Deltam_{32}^{2}| (eV^{2}/c^{4})");
  go2_90->  GetXaxis()->SetTitle("sin^{2}(2#theta_{23})");
  go2_90->  GetYaxis()->SetTitle("|#Deltam_{32}^{2}| (eV^{2}/c^{4})");
  go2_68->  GetXaxis()->SetTitle("sin^{2}(2#theta_{23})");
  go2_68->  GetYaxis()->SetTitle("|#Deltam_{32}^{2}| (eV^{2}/c^{4})");
  ho2_colz->GetYaxis()->SetTitleOffset(yoffset);
  go2_90->  GetYaxis()->SetTitleOffset(yoffset);
  go2_68->  GetYaxis()->SetTitleOffset(yoffset);
  //set the dchisq surface name
  ho2_colz->SetTitle("#Delta#chi^{2} surface for #theta_{23} #geq #pi/4, NH (2013 numu paper)");
  ho2_colz->SetName("h_deltachisq_o2");
  //best fit point as graph
  double bfo2_x[1] = {0.999};
  double bfo2_y[1] = {0.00244};
  TGraph * bfo2 = new TGraph(1, bfo2_x, bfo2_y);
  bfo2->SetTitle("best fit point for #theta_{23} #leq #pi/4, NH (2013 numu paper)");
  bfo2->SetName("g_bestfit_o2");
  bfo2->SetMarkerSize(1.8);
  bfo2->SetMarkerColor(kBlack);
  bfo2->SetMarkerStyle(29);
  bfo2->GetXaxis()->SetTitle("sin^{2}(2#theta_{23})");
  bfo2->GetYaxis()->SetTitle("|#Deltam_{32}^{2}| (eV^{2}/c^{4})");
  bfo2->GetYaxis()->SetTitleOffset(yoffset);
  //best fit point as marker
  TMarker * mo2 = new TMarker(0.999,0.00244,29);
  mo2->SetMarkerSize(1.8);
  mo2->SetMarkerColor(kBlack);
  cout << mo2->GetTitle() << endl;
  //mo2->SetTitle("best fit point for #theta_{23} #geq #pi/4, NH (2013 numu paper)");
  //mo2->SetName("m_bestfit_o2");

  //save
  TFile fout("t2k_2013paper_likelihoodratio_contour_and_deltachisq.root","RECREATE");
  //octant one
  bfo1->Write();
  go1_68->Write();
  go1_90->Write();
  //mo1->Write("m_bestfit_o1");
  ho1_colz->Write();
  //octant two
  bfo2->Write();
  go2_68->Write();
  go2_90->Write();
  //mo2->Write("m_bestfit_o2");
  ho2_colz->Write();
  //param space
  h_osc_param_space->Write();
  fout.Close();
}
