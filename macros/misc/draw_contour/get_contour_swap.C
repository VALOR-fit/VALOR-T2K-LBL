#include "save_contour.C"

void get_contour_swap(const char* _filename = "",bool draw_asimov_input = true)
{
  TGaxis::SetMaxDigits(3); // Sets maximum digits to 4 so we have exponent on axis lables
  
  //   char* filename = "toyA_joint";
  string filename = _filename;
  //removing .root from end of filename
  
  filename.replace(filename.find(".root"), sizeof(".root")-1, "");
  
  TCanvas* c = new TCanvas("c","c",0,0,1200,900);
  TFile f(Form("%s.root",filename.c_str()));
  TTree* t = f.Get("fitnt");
  
  TH2D* h1;
  TH2D* h2;
  TH2D* h3;
//  TH2D* h11;
//  TH2D* h33;
  TH2D* hdchisq;
  
  valor::FitParamList* fp = 0;
  valor::PhysParamList* ip = 0;
  double input_x = 0;
  double input_y = 0;
  
  t->SetBranchAddress("hdeltachisq2D",&hdchisq);
  t->SetBranchAddress("hchisq2D",&h1);
  t->SetBranchAddress("best_fit_params",&fp);
  t->SetBranchAddress("input_params",&ip);
  t->SetBranchAddress("truepar0",&input_x);
  t->SetBranchAddress("truepar1",&input_y);
  t->GetEntry(0);
  string file2 = filename;
  hdchisq->SetDirectory(0);
  h1->SetDirectory(0);

  SwapXY(h1);
  
  h3 = (TH2D*)h1->Clone();
  h2 = (TH2D*)h1->Clone();
//  h11 = (TH2D*)h1->Clone();
//  h33 = (TH2D*)h1->Clone();
  h3->SetDirectory(0);
  h2->SetDirectory(0);
//  h11->SetDirectory(0);
//  h33->SetDirectory(0);
  double chisq_min = h1->GetMinimum();
  h1->SetContour(1);
  h3->SetContour(1);
  
  cout << "x: " << input_x << endl;
  cout << "y: " << input_y << endl;
  
  cout << "Chisq min: " << chisq_min << endl;
  
  
  h1->SetContourLevel(0,chisq_min+2.3);
  h3->SetContourLevel(0,chisq_min+4.61);
  h1->SetLineWidth(2);
  h3->SetLineWidth(2);
  h1->SetLineColor(kRed);
  h3->SetLineColor(kRed);
  h1->SetLineStyle(2);
//  
//  h11->SetContour(1);
//  h33->SetContour(1);
//  
//  h11->SetContourLevel(0,2.3);
//  h33->SetContourLevel(0,4.61);
//  h11->SetLineWidth(2);
//  h33->SetLineWidth(2);
//  h11->SetLineColor(kRed);
//  h33->SetLineColor(kRed);
//  h11->SetLineStyle(2);
  
  int minbin = h1->GetMinimumBin();
  int binx;
  int biny;
  h1->GetBinXYZ(minbin,binx,biny,0);
  
  double xmin = h1->GetXaxis()->GetBinCenter(binx);
  double ymin = h1->GetYaxis()->GetBinCenter(biny);
  
  TMarker * bf = new TMarker(xmin,ymin,kFullStar);
  bf->SetMarkerColor(kRed);
  bf->SetMarkerSize(1.5);
  
  TMarker * bf2 = new TMarker(input_x,input_y,kFullStar);
  bf2->SetMarkerColor(kViolet);
  bf2->SetMarkerSize(1.5);
  
  format_hist(h2);
  format_hist(hdchisq);
  
  
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
  f.Close();
  
  TFile* f2 = new TFile(Form("%s_contour.root",filename.c_str()),"recreate");
  
  save_contour(h1,h3,"");
//  save_contour(h11,h33,"_dchi0");
  
  
  h2->Write("chisq_map_2D");
  hdchisq->Write("dchisq_map_2D");
  
  bf->Write("best_fit_point_marker");
  f2->Close();
  
//  cout << "xmin: " << xmin << "  ymin: " << ymin << endl;
//  hdchisq->SetBinContent(minbin,999);
//  
//  minbin = hdchisq->GetMinimumBin();
//  hdchisq->GetBinXYZ(minbin,binx,biny,0);
//  
//  xmin = hdchisq->GetXaxis()->GetBinCenter(binx);
//  ymin = hdchisq->GetYaxis()->GetBinCenter(biny);
//  
//  cout << "Second minimum:" << endl;
//  cout << "xmin: " << xmin << "  ymin: " << ymin << "  dchisq: " << hdchisq->GetBinContent(minbin) << endl;
  
}




TH2D* SwapXY(TH2D *graph){
  
  int nbinsX = graph->GetNbinsX();
  int nbinsY = graph->GetNbinsY();
  double xmin = graph->GetXaxis()->GetXmin();
  double xmax = graph->GetXaxis()->GetXmax();
  double ymin = graph->GetYaxis()->GetXmin();
  double ymax = graph->GetYaxis()->GetXmax();
  
  TH2D* graph_out = new TH2D("h_temp","",nbinsY,ymin,ymax,nbinsX,xmin,xmax);
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
  
  return 1;
}

int format_hist(TH1D* _histogram)
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
  
  return 1;
}