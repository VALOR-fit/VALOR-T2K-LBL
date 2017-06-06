//ToDo add input for legend names
// Produce VALOR contours with macros/misc/draw_contours.C
// and read the output files

using namespace std;
using namespace valor;
using namespace valor::conv;

static const int Ncont_68 = 3;
static const int Ncont_90 = 3;
static const int Nfiles = 2;


TGraph *g68_valor[Ncont_68][Nfiles] = {0};
TGraph *g90_valor[Ncont_90][Nfiles] = {0};


void CompareContour2D(char* file_1 = "test1",
                      char* file_2 = "test2",
                      char* tag = "",
                      const char* tag1   = "",
                      const char* tag2   = "",
                      const bool draw_reactors = false,
                      const bool doGlobalMin  = false,
                      const bool draw_preliminary = true,
                      const bool paper_style   = false

                      )
{
  TGaxis::SetMaxDigits(3);
  gStyle->SetLabelSize(0.05);

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
    cout << "This is done! " << endl;
    
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

  hchisq_valor[0]->GetYaxis()->SetLabelSize(0.04);
  hchisq_valor[0]->GetXaxis()->SetLabelSize(0.04);
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
    for(int icont=0;icont<Ncont_68;icont++){
      if(g68_valor[icont][ifile]){
        g68_valor[icont][ifile]->SetLineWidth(2);
        g68_valor[icont][ifile]->SetLineColor(ifile+1);
        g68_valor[icont][ifile]->SetLineStyle(2);
        g68_valor[icont][ifile]->DrawClone("l same");
      }
    }
    gbest_fit_global[ifile]->SetMarkerColor(ifile+1);
    gbest_fit_global[ifile]->SetMarkerStyle(kFullStar+ifile);
    gbest_fit_global[ifile]->Draw("same");
  }
  
  TLegend *lcontours = new TLegend(0.65,0.73,0.88,0.88);
  lcontours->SetFillColor(0);
  lcontours->SetLineColor(0);
  lcontours->SetTextFont(132);
  lcontours->SetTextSize(0.04);
  for(int ifile=0;ifile<Nfiles;ifile++){
    char name[200];
    sprintf(name,"%s - 68CL",latexname[ifile]);
    if (g68_valor[0][ifile])
    lcontours->AddEntry(g68_valor[0][ifile],name,"l");
    sprintf(name,"%s - 90CL",latexname[ifile]);
    if (g90_valor[0][ifile])
    lcontours->AddEntry(g90_valor[0][ifile],name,"l");
  }
  
  TLegend *lcontours_2 = new TLegend(0.5,0.77,0.65,0.84);
  lcontours_2->SetFillColor(0);
  lcontours_2->SetLineColor(0);
  lcontours_2->SetTextFont(132);
  lcontours_2->SetTextSize(0.04);

  lcontours_2->AddEntry(gbest_fit_global[0],"Best fit","p");
  if(draw_reactors) lcontours_2->AddEntry(breactors,"PDG 2015","f");

  lcontours->Draw();
  lcontours_2->Draw();

  char* mass_hierarchy_text = (doGlobalMin)? "Global minimum" :"Fixed Mass Ordering" ;
  TLatex *mass_hierarchy_tex = new TLatex(h1->GetXaxis()->GetBinLowEdge((int) (0.15*h1->GetXaxis()->GetNbins())),1.005*Y2,mass_hierarchy_text);
  
  TLatex *latPreliminary = new TLatex(h1->GetXaxis()->GetBinLowEdge((int) (0.56*h1->GetXaxis()->GetNbins())),1.005*Y2,"T2K Run1-8");
  mass_hierarchy_tex->SetTextSize(0.05);
  latPreliminary->SetTextSize(0.05);

  if(draw_preliminary) latPreliminary->Draw("same");
  mass_hierarchy_tex->Draw("same");

  //  TLatex *latcontours = new TLatex(0.075,2.,"Global scan");
  //  latcontours->SetTextSize(0.07);
  //  latcontours->Draw("same");
  if(tag != "") ccontours->SaveAs(Form("%s.pdf",tag));
  if(tag != "") ccontours->SaveAs(Form("%s.root",tag));

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


void minimise_hist(TH2* _h1, double min)
{
  for(int i =1;i<=(_h1->GetNbinsX());i++)
  {
     for(int j =1;j<=_h1->GetNbinsY();j++)
    {
    _h1->SetBinContent(i,j,_h1->GetBinContent(i,j)-min);
    }
  }
}

void get_contour_68(TH2* _h1, int file_num)
{
  TCanvas* c2 = new TCanvas("c2","c",0,0,900,600);
  _h1->Draw("cont3same");
  _h1->Draw("CONT LIST");
  c2->Update();
  TObjArray *contours = 0;
  contours = (TObjArray *)gROOT->GetListOfSpecials()->FindObject("contours");
  TList *lister = (TList*)contours->At(0);
  TGraph *g;
  int it=0;
  
  TIter next(lister);
  while ((g = (TGraph*)next()))
  {
    g68_valor[it][file_num] = (TGraph*) g->Clone();
    it++;
  }

  c2->Clear();
  delete c2;
}

void get_contour_90(TH2D* _h1, int file_num)
{
  TCanvas* c2 = new TCanvas("c2","c",0,0,900,600);
  _h1->Draw("cont3same");
  _h1->Draw("CONT LIST");
  c2->Update();
  TObjArray *contours = 0;
  contours = (TObjArray *)gROOT->GetListOfSpecials()->FindObject("contours");
  TList *lister = (TList*)contours->At(0);
  TGraph *g;
  int it=0;
  TIter next(lister);
  while ((g = (TGraph*)next()))
  {
    g90_valor[it][file_num] = (TGraph*) g->Clone();
    it++;
  }
  
  c2->Clear();
  delete c2;
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
