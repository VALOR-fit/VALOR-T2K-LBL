#include "save_contour.C"
int counter=0;
int format_hist(TH1D* _histogram);

void get_contour_2D(const char* _filename = "",bool draw_asimov_input = true, bool swap_xy = 0)
{
  TGaxis::SetMaxDigits(3); // Sets maximum digits to 4 so we have exponent on axis lables
  
  //   char* filename = "toyA_joint";
  string filename = _filename;
  //removing .root from end of filename
  
  filename.replace(filename.find(".root"), sizeof(".root")-1, "");
  
  TFile f(Form("%s.root",filename.c_str()));

  TTree* t = f.Get("fitnt");
  
  TH2D* h1;
  TH2D* h2;
  TH2D* h3;
  TH2D* h_2sigma;
  TH2D* h_3sigma;
  TH2D* hdchisq;
  TH2D* hlike;
  
  valor::FitParamList* fp = 0;
  valor::PhysParamList* ip = 0;
  double input_x = 0;
  double input_y = 0;
  
  string *trueparname0;
  string *trueparname1;
  string *trueMH;
  
  t->SetBranchAddress("hdeltachisq2D",&hdchisq);
  t->SetBranchAddress("hchisq2D",&h1);
  t->SetBranchAddress("best_fit_params",&fp);
  t->SetBranchAddress("input_params",&ip);
  t->SetBranchAddress("truepar0",&input_x);
  t->SetBranchAddress("truepar1",&input_y);
  t->SetBranchAddress("trueparname0",   &trueparname0   );
  t->SetBranchAddress("trueparname1",   &trueparname1   );
  t->SetBranchAddress("trueMH",   &trueMH   );

  
  t->GetEntry(0);
  string file2 = filename;
  hdchisq->SetDirectory(0);
  h1->SetDirectory(0);
  valor::conv::NuMassHierarchy_t hierarchy = util::HierarchyFromString(*trueMH);

  valor::conv::PhysParamId_t p = valor::util::PhysParamIdFromString(*trueparname0);
  string x_name = valor::util::AsStringOfficialAxis(p);/*,hierarchy);*/
  p = valor::util::PhysParamIdFromString(*trueparname1);
  string y_name = valor::util::AsStringOfficialAxis(p);/*,hierarchy);*/
  h1->GetXaxis()->SetTitle(x_name.c_str());
  h1->GetYaxis()->SetTitle(y_name.c_str());

  if(swap_xy)
  {
    h1 = SwapXY(h1);
    hdchisq = SwapXY(hdchisq);
  }
  h3 = (TH2D*)h1->Clone("h3");
  h2 = (TH2D*)h1->Clone("h2");
  h_2sigma = (TH2D*)h1->Clone("h_2sigma");
  h_3sigma = (TH2D*)h1->Clone("h_3sigma");

  h1->SetDirectory(0);
  h3->SetDirectory(0);
  h2->SetDirectory(0);
  h_2sigma->SetDirectory(0);
  h_3sigma->SetDirectory(0);
  hdchisq->SetDirectory(0);
  double chisq_min = h1->GetMinimum();
  h1->SetContour(1);
  h3->SetContour(1);
  h_2sigma->SetContour(1);
  h_3sigma->SetContour(1);
  
  hlike = (TH2D*)h2->Clone();
  hlike->SetDirectory(0);

  cout << "Start values: " << endl;
  cout << "x: " << input_x << endl;
  cout << "y: " << input_y << endl;
  
  cout << "Chisq min: " << chisq_min << endl;
  
  
  h1->SetContourLevel(0,chisq_min+2.3);
  h3->SetContourLevel(0,chisq_min+4.61);
  h_2sigma->SetContourLevel(0,chisq_min + 6.18);
  h_3sigma->SetContourLevel(0,chisq_min + 11.83);
  h1->SetLineWidth(2);
  h3->SetLineWidth(2);
  h_2sigma->SetLineWidth(2);
  h_3sigma->SetLineWidth(2);
  h1->SetLineColor(kRed);
  h3->SetLineColor(kRed);
  h_2sigma->SetLineColor(kRed);
  h_3sigma->SetLineColor(kRed);
  h1->SetLineStyle(2);
  h_2sigma->SetLineStyle(3);
  
  int minbin = h1->GetMinimumBin();
  int binx;
  int biny;
  h1->GetBinXYZ(minbin,binx,biny,0);
  
  double xmin = h1->GetXaxis()->GetBinCenter(binx);
  double ymin = h1->GetYaxis()->GetBinCenter(biny);
  
  cout << "Best fit values: " << endl;
  cout << "x: " << xmin << endl;
  cout << "y: " << ymin << endl;
 
  double scale = 0.0; // 30
  for(int i = 1;i<((h2->GetNbinsX() +1));i++)
  {
    for(int j = 1;j<(h2->GetNbinsY()+1);j++)
    {
      hlike->SetBinContent(i,j,TMath::Exp(scale - h2->GetBinContent(i,j)/2));
    }
  }
  double max = hlike->Integral();
  cout << "Max likelihood value " << max << endl;

  // Get 1D projections of 2D likelihood surface
  int nbins_x = hlike->GetNbinsX();
  int nbins_y = hlike->GetNbinsY();
  TH1D* hproj_x = new TH1D("hchisq_x", "hchisq_x", nbins_x,
        hlike->GetXaxis()->GetBinLowEdge(1),
        hlike->GetXaxis()->GetBinUpEdge(nbins_x));
  hproj_x->SetDirectory(0);
  hproj_x->GetXaxis()->SetTitle(hlike->GetXaxis()->GetTitle());
  hproj_x->GetYaxis()->SetTitle("#chi^{2}");

  for(int i = 1; i <= hlike->GetNbinsX(); i++)
  {
    double Lavg = 0.0;
    for(int j = 1; j <= hlike->GetNbinsY(); j++)
    {
      Lavg += hlike->GetBinContent(i,j);
    }
    Lavg /= nbins_y;
    hproj_x->SetBinContent(i, 2 * scale - 2 * TMath::Log(Lavg));
  }
  double chisq_min_x = hproj_x->GetMinimum();
  TH1D* hproj_x_delta = hproj_x->Clone("hdchisq_x");
  hproj_x_delta->SetDirectory(0);
  hproj_x_delta->GetYaxis()->SetTitle("#Delta#chi^{2}");
  for(int i = 1; i <= hproj_x->GetNbinsX(); i++)
  {
     hproj_x_delta->SetBinContent(i, hproj_x_delta->GetBinContent(i) - chisq_min_x);
  }

  TH1D* hproj_y = new TH1D("hchisq_y", "hchisq_y", nbins_y,
        hlike->GetYaxis()->GetBinLowEdge(1),
        hlike->GetYaxis()->GetBinUpEdge(nbins_y));
  hproj_y->SetDirectory(0);
  hproj_y->GetXaxis()->SetTitle(hlike->GetYaxis()->GetTitle());
  hproj_y->GetYaxis()->SetTitle("#chi^{2}");
  for(int i = 1; i <= hlike->GetNbinsY(); i++)
  {
    double Lavg = 0.0;
    for(int j = 1; j <= hlike->GetNbinsX(); j++)
    {
      Lavg += hlike->GetBinContent(j, i);
    }
    Lavg /= nbins_x;
    hproj_y->SetBinContent(i, 2 * scale - 2 * TMath::Log(Lavg));
  }
  double chisq_min_y = hproj_y->GetMinimum();
  TH1D* hproj_y_delta = hproj_y->Clone("hdchisq_y");
  hproj_y_delta->SetDirectory(0);
  hproj_y_delta->GetYaxis()->SetTitle("#Delta#chi^{2}");
  for(int i = 1; i <= hproj_y->GetNbinsX(); i++)
  {
     hproj_y_delta->SetBinContent(i, hproj_y_delta->GetBinContent(i) - chisq_min_y);
  }
  for(int i = 1;i<((h2->GetNbinsX() +1));i++)
  {
    for(int j = 1;j<(h2->GetNbinsY()+1);j++)
    {
      hlike->SetBinContent(i,j, hlike->GetBinContent(i,j) / max);
    }
  }
  cout << "Adjusted max likelihood value " << hlike->GetMaximum();
  
  TMarker * bf = new TMarker(xmin,ymin,kFullStar);
  bf->SetMarkerColor(kRed);
  bf->SetMarkerSize(1.5);
  
  TMarker * bf2 = new TMarker(input_x,input_y,kFullStar);
  bf2->SetMarkerColor(kViolet);
  bf2->SetMarkerSize(1.5);
  
  
  format_hist(h2);
  format_hist(hdchisq);
  format_hist(hlike);
 
  TCanvas* c = new TCanvas("c","c",0,0,1200,900);
  c->SetGrid();
 
  TLegend* leg = new TLegend(0.7,0.8,0.9,0.98);
  leg->AddEntry(bf,"Best fit point","p");
  if(draw_asimov_input) leg->AddEntry(bf2,"True asimov point","p");
  leg->AddEntry(h1,"68%","l");
  leg->AddEntry(h3,"90%","l");
  h2->Draw("axis");
  h1->Draw("cont3same");
  h3->Draw("cont3same");
  bf->Draw();
  if(draw_asimov_input)  bf2->Draw();
  
  leg->Draw();
  
  c->SaveAs(Form("%s_contour.pdf",filename.c_str()));
 
  // Plot likelihood scan
  TCanvas* c2 = new TCanvas("c","c",0,0,1200,900);
  c2->SetGrid();
  hlike->Draw("COLZ");
  c2->SaveAs(Form("%s_likelihood.pdf",filename.c_str()));

  // Plot deltachi2 surface
  TCanvas* c3 = new TCanvas("c","c",0,0,1200,900);
  c3->SetGrid();
  hdchisq->SetMinimum(4.0);
  hdchisq->SetMaximum(6.0);
  hdchisq->Draw("COLZ");
  h1->Draw("cont3same");
  h3->Draw("cont3same");
  c3->SaveAs(Form("%s_dchi2surface.pdf",filename.c_str()));

  f.Close();
  
  TFile* f2 = new TFile(Form("%s_contour.root",filename.c_str()),"recreate");
  
  save_contour(h1, "g68_valor");
  save_contour(h3, "g90_valor");
  save_contour(h_2sigma, "g95_valor");
  save_contour(h_3sigma, "g99_valor");
  
  h2->Write("chisq_map_2D");
  hdchisq->Write("dchisq_map_2D");
  hlike->Write("likelihood_2D");
  bf->Write("best_fit_point_marker");
  hproj_x->Write("hchisq_x");
  hproj_y->Write("hchisq_y");
  hproj_x_delta->Write("hdchisq_x");
  hproj_y_delta->Write("hdchisq_y");
  f2->Close();
  
  delete hproj_x;
  delete hproj_y;
  delete hproj_x_delta;
  delete hproj_y_delta;
}


int format_hist(TH2D* _histogram)
{
  _histogram->GetYaxis()->SetTitleOffset(0.77);
  _histogram->GetXaxis()->SetLabelSize(0.03);
  _histogram->GetYaxis()->SetLabelSize(0.03);
  _histogram->GetXaxis()->SetLabelFont(62);
  _histogram->GetYaxis()->SetLabelFont(62);
  _histogram->GetXaxis()->SetTitleSize(0.04);
  _histogram->GetXaxis()->SetTitleFont(62);
  _histogram->GetYaxis()->SetTitleSize(0.04);
  _histogram->GetYaxis()->SetTitleFont(62);
  
}


TGraph *SwapXY(TGraph *graph){
  TGraph *graph_swap = new TGraph();
  for(int ipt=0;ipt<graph->GetN();ipt++){
    double x = 0.;
    double y = 0.;
    graph->GetPoint(ipt,x,y);
    graph_swap->SetPoint(ipt,y,x);
  }
  return graph_swap;
}

void SwapXY(TMarker *graph){
  double x = graph->GetX();
  double y = graph->GetY();
  graph->SetX(y);
  graph->SetY(x);
  
  return;
}

TH2D* SwapXY(TH2D *graph){
  
  int nbinsX = graph->GetNbinsX();
  int nbinsY = graph->GetNbinsY();
  double xmin = graph->GetXaxis()->GetXmin();
  double xmax = graph->GetXaxis()->GetXmax();
  double ymin = graph->GetYaxis()->GetXmin();
  double ymax = graph->GetYaxis()->GetXmax();
  
  double arrayy[500];
  double arrayx[500];
  arrayx[0] = xmin;
  arrayy[0] = ymin;
  
  for(int i=1;i<=nbinsX;i++)
  {
    arrayx[i] = graph->GetXaxis()->GetBinUpEdge(i);
  }
  
  for(int i=1;i<=nbinsY;i++)
  {
    arrayy[i] = graph->GetYaxis()->GetBinUpEdge(i);
  }
  
  TH2D* graph_out = new TH2D(Form("h_temp_%d",counter++),"",nbinsY,arrayy,nbinsX,arrayx);
  
  graph_out->GetXaxis()->SetTitle(graph->GetYaxis()->GetTitle());
  graph_out->GetYaxis()->SetTitle(graph->GetXaxis()->GetTitle());
  
  for(int i=0;i<graph_out->GetNbinsX();i++){
    for(int j=0;j<graph_out->GetNbinsY();j++){
      double bincont = graph->GetBinContent(j+1,i+1);
      graph_out->SetBinContent(i+1,j+1,bincont);
    }
  }
  return graph_out;
  
}
