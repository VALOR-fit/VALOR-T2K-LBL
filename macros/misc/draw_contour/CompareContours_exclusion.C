
// Produce VALOR contours with macros/misc/draw_contours.C
// and read the output files
const int Ncont_68 = 3;
const int Ncont_90 = 3;

const int Nfiles = 2;


TGraph *g68_valor[Ncont_68][Nfiles] = {0};
TGraph *g90_valor[Ncont_90][Nfiles] = {0};


using namespace std;
using namespace valor;
using namespace valor::conv;
void CompareContours_exclusion(char* file_1 = "/Users/rajs/work/studies/final_2015_joint/fit/ccpi/data/cp13_nu/data_fit_nh_a-_0_0_contour.root",
                              char* file_2 = "/Users/rajs/work/studies/final_2015_joint/fit/ccpi/data/cp13_nubar/data_fit_nh_a-_0_0_contour.root",
                              char* tag = "",
                              const char* tag1   = "",
                              const char* tag2   = "",
                              const bool draw_reactors = true,
                              const bool doGlobalMin  = false,
                              const bool draw_preliminary = true,
                              const bool draw_68_contour = true,
                              const bool draw_nu_nubar = true, // Option to draw only the contours, not exclusion bands
                              const bool draw_mass_ordering = true,
                              const char* moretext   = "",
                              const bool paper_style   = false
                              )
{
  
  const bool scalegraphX = false; // scale axis of graphs (e.g. dmsq x 10+3)
  const bool scalegraphY = false; // scale axis of graphs (e.g. dmsq x 10+3)
  const bool draw_best_fit = false;     //Draw best fit
  const bool setinit = true; // true --> set init histo by hand, false --> use already defined set up
                                   // // sinsqth23 Vs dmsq
                                   // const double x_low = 0.3;
                                   // const double x_up = 0.7;
                                   // const double y_low = 1.95;
                                   // const double y_up = 3.05;
                                   // const char *title_x = "sin^{2} #theta_{23}";
                                   // const char *title_y = "|#Delta m^{2}_{32}| (10^{-3} eV^{2}/c^{4})";
  
  // // sinsqth23 Vs dcp (reactors)
  // const double x_low = 0.3;
  // const double x_up = 0.7;
  // const double y_low = -TMath::Pi();
  // const double y_up = +TMath::Pi();
  // const char *title_x = "sin^{2} #theta_{23}";
  // const char *title_y = "#delta_{CP} (radians)";
  
  // // sinsqth13 Vs dcp (reactors)
  // const double x_low = 0.014999;
  // const double x_up = 0.0360001;
  // const double y_low = -TMath::Pi();
  // const double y_up = +TMath::Pi();
  // const char *title_x = "sin^{2} #theta_{13}";
  // const char *title_y = "#delta_{CP} (radians)";
  
  // // sinsqth13 Vs dmsq (reactors)
  // const double x_low = 1.95;
  // const double x_up = 3.05;
  // const double y_low = -TMath::Pi();
  // const double y_up = +TMath::Pi();
  // const char *title_x = "|#Delta m^{2}_{32}| (10^{-3} eV^{2}/c^{4})";
  // const char *title_y = "#delta_{CP} (radians)";
  
  // sinsqth13 Vs dcp (t2k only)
  const double x_low = 0.;
  const double x_up = 0.105;
  const double y_low = -TMath::Pi();
  const double y_up = +TMath::Pi();
  //  const double y_low = -9;
  //  const double y_up = 9;
  
  const char *title_x = "sin^{2} #theta_{13}";
  const char *title_y = "#delta_{CP} (radians)";
  
  TString str_preliminary = "T2K Run1-7 preliminary";

  int Colors[2] = {kRed,kBlue};
  bool draw_reactor_band = false;
  
  TGaxis::SetMaxDigits(4); // Sets maximum digits to 4 so we have exponent on axis lables
  
  //////////////////////////// NEW VALOR contours /////////////////////
  
  // Global scan
  char *inname[Nfiles] = {file_1,file_2};
  const char *contour_68_name[Ncont_68] = { // different graph for each contour (+-)
    "g68_valor_0",
    "g68_valor_1",
    "g68_valor_2"
    
  };
  const int Ncont_90 = 3;
  const char *contour_90_name[Ncont_90] = { // different graph for each contour (+-)
    "g90_valor_0",
    "g90_valor_1",
    "g90_valor_2"
    
  };
  
  char *latexname[Nfiles] = {
    "CC#pi^{+}",
    "1Ring"
  };
  
  if(tag1 != "") latexname[0] = tag1;
  if(tag2 != "") latexname[1] = tag2;
  
  const char *chisq_name = "chisq_map_2D";
  
  const char *bestfitpt_name = "best_fit_point_marker";
  
  ////////////////////////////////////////////////////////////////////////
  
  
  
  // Read contours
  
  TFile *filein[Nfiles];
  TH2D *hchisq_valor[Nfiles];
  TMarker *gbest_fit_global[Nfiles];
  
  TH2D* h2 = 0 ;
  TH2D* h1 = 0 ;
  if(!doGlobalMin)
  {
    for(int ifile=0;ifile<Nfiles;ifile++){
      cout << inname[ifile] << endl;
      filein[ifile] = new TFile(inname[ifile],"READ");
      for(int icont=0;icont<Ncont_68;icont++){
        g68_valor[icont][ifile] = (TGraph*)filein[ifile]->Get(contour_68_name[icont]);
      }
      for(int icont=0;icont<Ncont_90;icont++){
        g90_valor[icont][ifile] = (TGraph*)filein[ifile]->Get(contour_90_name[icont]);
      }
      hchisq_valor[ifile] = (TH2D*)filein[ifile]->Get(chisq_name);
      hchisq_valor[ifile]->Reset(); // used to set the full range though tgraphs
      
      gbest_fit_global[ifile] = (TMarker*)filein[ifile]->Get(bestfitpt_name);
    }
    
    h1 = (TH2D*)filein[0]->Get("chisq_map_2D");
    h2 = (TH2D*)filein[1]->Get("chisq_map_2D");
    assert(h1);
    assert(h2);
    h1->SetDirectory(0);
    h2->SetDirectory(0);
    
  }
  else{
    
    for(int ifile=0;ifile<Nfiles;ifile++){
      filein[ifile] = new TFile(inname[ifile],"READ");
      hchisq_valor[ifile] = (TH2D*)filein[ifile]->Get(chisq_name);
      hchisq_valor[ifile]->Reset(); // used to set the full range though tgraphs
      gbest_fit_global[ifile] = (TMarker*)filein[ifile]->Get(bestfitpt_name);
    }
    
    h1 = (TH2D*) filein[0]->Get("chisq_map_2D");
    h2 = (TH2D*) filein[1]->Get("chisq_map_2D");
    assert(h1);
    assert(h2);
    h1->SetDirectory(0);
    h2->SetDirectory(0);
    
    double min = (h1->GetMinimum() < h2->GetMinimum()) ? h1->GetMinimum() : h2->GetMinimum();
    minimise_hist(h1,min);
    minimise_hist(h2,min);
    
    TH2D* h1_2 = (TH2D*)h1->Clone();
    TH2D* h2_2 = (TH2D*)h2->Clone();
    h1_2->SetDirectory(0);
    h2_2->SetDirectory(0);
    
    h1->SetContour(1);
    h1->SetContourLevel(0,2.3);
    h1_2->SetContour(1);
    h1_2->SetContourLevel(0,4.61);
    h2->SetContour(1);
    h2->SetContourLevel(0,2.3);
    h2_2->SetContour(1);
    h2_2->SetContourLevel(0,4.61);
    
    get_contour_68(h1,0);
    get_contour_68(h2,1);
    get_contour_90(h1_2,0);
    get_contour_90(h2_2,1);
    
  }
  
  
  // Global scan contours
  
  TCanvas *ccontours = new TCanvas("ccontours","ccontours");
  hchisq_valor[0]->GetXaxis()->SetLabelSize(0.04); // to set the range
  hchisq_valor[0]->GetYaxis()->SetLabelSize(0.04); // to set the range
  if(paper_style) format_hist(hchisq_valor[0]);
  hchisq_valor[0]->Draw(); // to set the range
  
  const double sinsqth13_nom = 0.0219;
  const double sinsqth13_err = 0.0012;
  double X1 = sinsqth13_nom + sinsqth13_err;
  double X2 = sinsqth13_nom - sinsqth13_err;
  double Y1 = hchisq_valor[0]->GetYaxis()->GetXmin();
  double Y2 = hchisq_valor[0]->GetYaxis()->GetXmax();
  TBox *breactors = new TBox(X1,Y1+0.05,X2,Y2-0.05);
  breactors->SetFillColor(kYellow+1);
  breactors->SetLineColor(kYellow+1);
  
  if(draw_reactors) breactors->Draw("same");
  
  
  for(int ifile=0;ifile<Nfiles;ifile++){
    for(int icont=0;icont<Ncont_90;icont++){
      if(g90_valor[icont][ifile]){
        g90_valor[icont][ifile]->SetLineWidth(2);
        g90_valor[icont][ifile]->SetLineColor(ifile+1);
        g90_valor[icont][ifile]->SetLineStyle(1);
        g90_valor[icont][ifile]->DrawClone("l same");
      }
    }
    if(draw_68_contour){
      for(int icont=0;icont<Ncont_68;icont++){
        if(g68_valor[icont][ifile]){
          g68_valor[icont][ifile]->SetLineWidth(2);
          g68_valor[icont][ifile]->SetLineColor(ifile+1);
          g68_valor[icont][ifile]->SetLineStyle(2);
          g68_valor[icont][ifile]->DrawClone("l same");
        }
      }
    }
    gbest_fit_global[ifile]->Draw("same");
    gbest_fit_global[ifile]->SetMarkerColor(Colors[ifile]);
  }
  
  TLegend *lcontours = new TLegend(0.6,0.25,0.88,0.48);
  lcontours->SetFillColor(0);
  lcontours->SetLineColor(0);
  lcontours->SetTextFont(132);
  lcontours->SetTextSize(0.04);
  for(int ifile=0;ifile<Nfiles;ifile++){
    char name[200];
    if(draw_68_contour){
      sprintf(name,"%s - 68CL",latexname[ifile]);
      if (g68_valor[0][ifile])
        lcontours->AddEntry(g68_valor[0][ifile],name,"l");
    }
    sprintf(name,"%s - 90CL",latexname[ifile]);
    if (g90_valor[0][ifile])
      lcontours->AddEntry(g90_valor[0][ifile],name,"l");
  }
  
  TLegend *lcontours_2 = new TLegend(0.67,0.52,0.88,0.6);
  lcontours_2->SetFillColor(0);
  lcontours_2->SetLineColor(0);
  lcontours_2->SetTextFont(132);
  lcontours_2->SetTextSize(0.04);
  
  lcontours_2->AddEntry(gbest_fit_global[0],"Best fit","p");
  if(draw_reactors) lcontours_2->AddEntry(breactors,"PDG 2015","f");
  
  
  char* mass_hierarchy_text = (doGlobalMin)? "Global minimum" :"Fixed Mass Ordering" ;
  TLatex *mass_hierarchy_tex = new TLatex(h1->GetXaxis()->GetBinLowEdge((int) (0.15*h1->GetXaxis()->GetNbins())),1.005*Y2,mass_hierarchy_text);
  
  TLatex *latPreliminary = new TLatex(h1->GetXaxis()->GetBinLowEdge((int) (0.56*h1->GetXaxis()->GetNbins())),1.005*Y2,"T2K Run1-7 Preliminary");
  mass_hierarchy_tex->SetTextSize(0.05);
  latPreliminary->SetTextSize(0.05);

  TLatex *latMoreText = new TLatex(h1->GetXaxis()->GetBinLowEdge((int) (0.8*h1->GetXaxis()->GetNbins())),0.8*Y2,moretext);
  latMoreText->SetTextSize(0.05);

  
  if(draw_preliminary) latPreliminary->Draw("same");
  if(draw_mass_ordering)  mass_hierarchy_tex->Draw("same");
  latMoreText->Draw("same");
  //  TLatex *latcontours = new TLatex(0.075,2.,"Global scan");
  //  latcontours->SetTextSize(0.07);
  //  latcontours->Draw("same");
  
  const int fill_style[2][2] =
  {
    {-3005,3005},
    {-3005,3005},
  };
  
  TGraph* g1 = new TGraph();
  int g1_counter=0;
    g1->SetPoint(g1_counter++,0.022,+50);
  
  for(int i = 0;i<g90_valor[0][0]->GetN(); i++)
  {
    double x;
    double y;
    g90_valor[0][0]->GetPoint(i,x,y);
    g1->SetPoint(g1_counter++,x,y);
  }
    g1->SetPoint(g1_counter++,0.022,-50);
  
  if(g90_valor[1][0]){
  for(int i = 0;i<g90_valor[1][0]->GetN(); i++)
  {
    double x;
    double y;
    g90_valor[1][0]->GetPoint(i,x,y);
    g1->SetPoint(g1_counter++,x,y);
  }
  }
  double x;
  double y;
  g1->GetPoint(0,x,y);
  g1->SetPoint(g1_counter++,x,y);
  
  g1->SetLineColor(kRed);
  g1->SetLineWidth(-1002);
  g1->SetFillStyle(3004);
  //  g1->SetFillStyle(3554);
  g1->SetFillColor(kRed);
  if(!draw_nu_nubar)
    g1->Draw("l");
  
  
  TGraph* g2 = new TGraph();
  int g2_counter=0;
  //  g2->SetPoint(g2_counter++,0.02,9);
  for(int i = 0;i<g90_valor[0][1]->GetN(); i++)
  {
    double x;
    double y;
    g90_valor[0][1]->GetPoint(i,x,y);
    g2->SetPoint(g2_counter++,x,y);
  }
  //  g2->SetPoint(g2_counter++,0.02,-9);
  
  g2->SetLineColor(kBlue);
  g2->SetLineWidth(-1002);
  g2->SetFillStyle(3005);
  //  g2->SetFillStyle(3545);
  g2->SetFillColor(kBlue);
  if(!draw_nu_nubar)
    g2->Draw("l");
  
  g90_valor[0][1]->SetLineColor(kBlue);
  g90_valor[0][1]->SetLineWidth(2);
  g90_valor[0][1]->Draw("l");
  
  g90_valor[0][0]->SetLineColor(kRed);
  g90_valor[0][0]->SetLineWidth(2);
  g90_valor[0][0]->Draw("l");
  if(g90_valor[1][0]){
  g90_valor[1][0]->Draw("l");
  g90_valor[1][0]->SetLineColor(kRed);
  g90_valor[1][0]->SetLineWidth(2);
  }
  if(draw_68_contour){
    g68_valor[0][0]->SetLineColor(kRed);
    g68_valor[0][0]->SetLineWidth(2);
    g68_valor[0][0]->Draw("l");
    if(g68_valor[1][0])
    {
      g68_valor[1][0]->SetLineWidth(2);
      g68_valor[1][0]->SetLineColor(kRed);
      g68_valor[1][0]->Draw("l");
    }
    g68_valor[0][1]->SetLineColor(kBlue);
    g68_valor[0][1]->SetLineWidth(2);
    g68_valor[0][1]->Draw("l");
  }
  
  TGraph* gdashed = new TGraph();
  gdashed->SetLineStyle(2);
  gdashed->SetLineWidth(2);
  TGraph* gsolid = new TGraph();
  gsolid->SetLineStyle(2);
  
  
  lcontours->Draw();
  lcontours_2->Draw();
  
  
  
  if(tag != "") ccontours->SaveAs(Form("%s.pdf",tag));
  
  
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

void SwapXY(TGraph *graph){
  for(int ipt=0;ipt<graph->GetN();ipt++){
    double x = 0.;
    double y = 0.;
    graph->GetPoint(ipt,x,y);
    graph->RemovePoint(ipt);
    graph->SetPoint(ipt,y,x);
  }
  return;
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

void SwapXY(TMarker *graph){
  double x = graph->GetX();
  double y = graph->GetY();
  graph->SetX(y);
  graph->SetY(x);
  
  return;
}

void ScaleMarkerAxis(TMarker *marker,double scale_x,double scale_y){
  double x = marker->GetX();
  double y = marker->GetY();
  double x_new = x * scale_x;
  double y_new = y * scale_y;
  marker->SetX(x_new);
  marker->SetY(y_new);
  return;
}
void ScaleGraphAxis(TGraph *graph,double scale_x,double scale_y){
  for(int ipt=0;ipt<graph->GetN();ipt++){
    double x = 0.;
    double y = 0.;
    graph->GetPoint(ipt,x,y);
    double x_new = x * scale_x;
    double y_new = y * scale_y;
    graph->RemovePoint(ipt);
    graph->SetPoint(ipt,x_new,y_new);
  }
  return;
}

int format_hist(TH2D* _histogram)
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