
//
// ------------------------------------------------------------------------------------------------
//
// Read and analyze the output from run_fit_marginalization_readbinary_*.C
//
// Davide Sgalaberna
//
// ------------------------------------------------------------------------------------------------

using namespace std;
using namespace valor;
using namespace valor::conv;
TCanvas *ccontours = new TCanvas("ccontours","ccontours");

void CompareContour1D
(
 const char* file_1 = "/Users/rajs/work/studies/final_2015_joint/fit/ccpi/data/dcp_reactor/data_nh_a-_0_0_contour.root",
 const char* file_2 = "/Users/rajs/work/studies/final_2015_joint/fit/ccpi/data/dcp_reactor/data_nh_a-_0_0_contour.root",
 const char* tag    = "",
 const char* tag1   = "",
 const char* tag2   = "",
 
 const bool doGlobalMin  = false,
 const bool draw_preliminary = true,
 
 bool draw_1sigmaCL = true,
 bool draw_90CL     = true,
 bool draw_2sigmaCL = false,
 bool draw_3sigmaCL = true,
 bool draw_4sigmaCL = false,
 bool draw_5sigmaCL = true,
 
 const bool draw_reactors = false,
 const bool paper_style   = false
 )
{
  TGaxis::SetMaxDigits(3);
  gStyle->SetLabelSize(0.05);

  const int Nfiles = 2;
  TString latexname[Nfiles] = {
    "Normal",
    "Inverted"
  };
  
  if(tag1 != "") latexname[0] = tag1;
  if(tag2 != "") latexname[1] = tag2;
  

  TFile f1(file_1,"read");
  TFile f2(file_2,"read");
  
  string histname;
  
  if(doGlobalMin) histname = "h_chisq";
  else histname = "h_deltachisq";
  
  TH1D* h1 = f1.Get(histname.c_str());
  TH1D* h2 = f2.Get(histname.c_str());
  h1->SetLineColor(kBlack);
  h2->SetLineColor(kRed);
  h1->SetLineWidth(2);
  h2->SetLineWidth(2);
  
  cout << "Hist1Min: " << std::setprecision (8) << h1->GetMinimum() << " (" << h1->GetXaxis()->GetBinCenter(h1->GetMinimumBin()) << ") " << endl;
  cout << "Hist2Min: " << std::setprecision (8) << h2->GetMinimum() << " (" << h2->GetXaxis()->GetBinCenter(h2->GetMinimumBin()) << ") " << endl;

  if(doGlobalMin){
  double min = (h1->GetMinimum() < h2->GetMinimum()) ? h1->GetMinimum() : h2->GetMinimum();
    minimise_hist(h1,min);
    minimise_hist(h2,min);
  }
  
  
  valor::RootFinder* find_root = new valor::RootFinder();
  find_root->SetGraph(h1);
  
  find_root->SetDchiSq(1);
  find_root->GetRoots();
  cout << "Hist1: " << endl;
  find_root->PrintRoots();
  delete find_root;
  
  valor::RootFinder* find_root2 = new valor::RootFinder();
  find_root2->SetGraph(h2);
  
  find_root2->SetDchiSq(1);
  find_root2->GetRoots();
  cout << "Hist2: " << endl;
  find_root2->PrintRoots();
  delete find_root2;
  
  double max;
  double max1 = h1->GetMaximum();
  double max2 = h2->GetMaximum();
  double max_hist = (max1>max2)? max1 : max2;
  if(max_hist>17) max = 30;
  else max = max_hist*1.4;
  
  if(max_hist<1.3){
    draw_1sigmaCL=false;
    draw_90CL=false;
    draw_2sigmaCL=false;
    draw_3sigmaCL=false;
    draw_4sigmaCL=false;
    draw_5sigmaCL=false;
  }
  

  
  if(max_hist<3.5){
    draw_90CL=false;
    draw_2sigmaCL=false;
    draw_3sigmaCL=false;
    draw_4sigmaCL=false;
    draw_5sigmaCL=false;
  }
  
  if(max_hist<4.6){
    draw_2sigmaCL=false;
    draw_3sigmaCL=false;
    draw_4sigmaCL=false;
    draw_5sigmaCL=false;
  }

  if(max_hist<13){
    draw_3sigmaCL=false;
    draw_4sigmaCL=false;
    draw_5sigmaCL=false;
  }
  
  if(max_hist<17.6){
    draw_4sigmaCL=false;
    draw_5sigmaCL=false;
  }

  


  h1->GetYaxis()->SetTitle("-2#Deltaln(L)");
//  h1->GetYaxis()->SetLabelSize(0.04);
//  h1->GetXaxis()->SetLabelSize(0.04);

  h1->GetYaxis()->SetRangeUser(0,max);
  if(paper_style) format_hist(h1);
  h1->Draw("l");
  
  const double sinsqth13_nom = 0.0219;
  const double sinsqth13_err = 0.0012;
  double X1 = sinsqth13_nom - sinsqth13_err;
  double X2 = sinsqth13_nom + sinsqth13_err;
  double Y1 = h1->GetMinimum();
  double Y2 = max;
  TBox *breactors = new TBox(X1,Y1,X2,Y2);
  breactors->SetFillColor(kYellow+1);
  
  if(draw_reactors) breactors->Draw("same");
  h1->Draw("lsame");
  h2->SetLineStyle(2);
  h2->Draw("samel");
  
  TLine *l1sigma = new TLine(h1->GetXaxis()->GetBinLowEdge(1),1.,h1->GetXaxis()->GetBinUpEdge(h1->GetNbinsX()),1.);
  TLine *l90 = new TLine(h1->GetXaxis()->GetBinLowEdge(1),2.71,h1->GetXaxis()->GetBinUpEdge(h1->GetNbinsX()),2.71);
  TLine *l2sigma = new TLine(h1->GetXaxis()->GetBinLowEdge(1),4.,h1->GetXaxis()->GetBinUpEdge(h1->GetNbinsX()),4.);
  TLine *l3sigma = new TLine(h1->GetXaxis()->GetBinLowEdge(1),9.,h1->GetXaxis()->GetBinUpEdge(h1->GetNbinsX()),9.);
  TLine *l4sigma = new TLine(h1->GetXaxis()->GetBinLowEdge(1),16.,h1->GetXaxis()->GetBinUpEdge(h1->GetNbinsX()),16.);
  TLine *l5sigma = new TLine(h1->GetXaxis()->GetBinLowEdge(1),25.,h1->GetXaxis()->GetBinUpEdge(h1->GetNbinsX()),25.);
  l1sigma->SetLineWidth(2);
  l90->SetLineWidth(2);
  l2sigma->SetLineWidth(2);
  l3sigma->SetLineWidth(2);
  l4sigma->SetLineWidth(2);
  l5sigma->SetLineWidth(2);
  l1sigma->SetLineColor(kBlue);
  l90->SetLineColor(kBlue);
  l2sigma->SetLineColor(kBlue);
  l3sigma->SetLineColor(kBlue);
  l4sigma->SetLineColor(kBlue);
  l5sigma->SetLineColor(kBlue);
  if(draw_1sigmaCL) l1sigma->Draw();
  if(draw_90CL)     l90->Draw();
  if(draw_2sigmaCL) l2sigma->Draw();
  if(draw_3sigmaCL) l3sigma->Draw();
  if(draw_4sigmaCL) l4sigma->Draw();
  if(draw_5sigmaCL) l5sigma->Draw();
  
  TLatex *latLine1Sigma = new TLatex(h1->GetXaxis()->GetBinLowEdge(5),1.1,"1#sigma");
  TLatex *latLine90CL = new TLatex(h1->GetXaxis()->GetBinLowEdge(5),2.81,"90% CL");
  TLatex *latLine2Sigma = new TLatex(h1->GetXaxis()->GetBinLowEdge(5),4.1,"2#sigma");
  TLatex *latLine3Sigma = new TLatex(h1->GetXaxis()->GetBinLowEdge(5),9.1,"3#sigma");
  TLatex *latLine4Sigma = new TLatex(h1->GetXaxis()->GetBinLowEdge(5),16.1,"4#sigma");
  TLatex *latLine5Sigma = new TLatex(h1->GetXaxis()->GetBinLowEdge(5),25.1,"5#sigma");
  latLine1Sigma->SetTextColor(kBlue);
  latLine90CL->SetTextColor(kBlue);
  latLine2Sigma->SetTextColor(kBlue);
  latLine3Sigma->SetTextColor(kBlue);
  latLine4Sigma->SetTextColor(kBlue);
  latLine5Sigma->SetTextColor(kBlue);
  latLine1Sigma->SetTextSize(0.05);
  latLine90CL->SetTextSize(0.05);
  latLine2Sigma->SetTextSize(0.05);
  latLine3Sigma->SetTextSize(0.05);
  latLine4Sigma->SetTextSize(0.05);
  latLine5Sigma->SetTextSize(0.05);
  if(draw_1sigmaCL) latLine1Sigma->Draw("same");
  if(draw_90CL)     latLine90CL->Draw("same");
  if(draw_2sigmaCL) latLine2Sigma->Draw("same");
  if(draw_3sigmaCL) latLine3Sigma->Draw("same");
  if(draw_4sigmaCL) latLine4Sigma->Draw("same");
  if(draw_5sigmaCL) latLine5Sigma->Draw("same");
  TLegend *lcontours = new TLegend(0.6,0.77,0.75,0.9);
  lcontours->SetFillColor(0);
  lcontours->SetLineColor(0);
  lcontours->SetTextFont(132);
  lcontours->SetTextSize(0.05);
  lcontours->AddEntry(h1,latexname[0]);
  lcontours->AddEntry(h2,latexname[1]);
  if(draw_reactors)
  {
     lcontours->AddEntry(breactors, "PDG 2016", "f");
  }
  lcontours->Draw();
 
  //char* note_text = (doGlobalMin) ? "Global minimum" :"Fixed Mass Ordering";
  char* notetext = "Final systematics pending";
  double tx = h1->GetXaxis()->GetBinLowEdge(h1->GetXaxis()->GetNbins() / 2);
  double nx = h1->GetXaxis()->GetBinLowEdge(h1->GetXaxis()->GetNbins() / 10);
  double ty = max;

  TLatex* title = new TLatex(tx, ty * 1.01, "T2K Run1-8 Preliminary");
  title->SetTextSize(0.05);
  title->Draw("same");
  TLatex* note = new TLatex(nx, ty * 1.01, "Final systematics pending");
  note->SetTextSize(0.025);
//  note->Draw("same");

  if(draw_preliminary)
  {
     title->Draw("same");
  }
  ccontours->SaveAs(Form("%s.pdf",tag));
  ccontours->SaveAs(Form("%s.png",tag));
  ccontours->SaveAs(Form("%s.eps",tag));
  ccontours->SaveAs(Form("%s.root",tag));

  return;
}

void minimise_hist(TH1* h1, double min)
{
  for(int i =1;i<=h1->GetNbinsX();i++)
  {
    h1->SetBinContent(i,h1->GetBinContent(i)-min);
  }
}



int format_hist(TH1D* _histogram)
{
  _histogram->GetYaxis()->SetTitleOffset(0.85);
  _histogram->GetXaxis()->SetTitleOffset(1);
  _histogram->GetXaxis()->SetLabelSize(0.05);
  _histogram->GetYaxis()->SetLabelSize(0.05);
  _histogram->GetXaxis()->SetLabelFont(132);
  _histogram->GetYaxis()->SetLabelFont(132);
  _histogram->GetXaxis()->SetTitleFont(132);
  _histogram->GetYaxis()->SetTitleFont(132);
  _histogram->GetXaxis()->SetTitleSize(0.06);
  _histogram->GetYaxis()->SetTitleSize(0.06);
  _histogram->GetXaxis()->CenterTitle();
  _histogram->GetYaxis()->CenterTitle();
  _histogram->GetZaxis()->CenterTitle();

}
