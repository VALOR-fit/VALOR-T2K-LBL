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
                      const bool paper_style   = true

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

  const char* trueMH1;
  const char* trueMH2;
  if(strstr(file_1, "normal") != NULL)
  {
     trueMH1 = "normal";
     std::cout << "File 1:" << file_1 << " " << trueMH1 << std::endl;
  }
  else
  {
     trueMH1 = "inverted";
     std::cout << "File 1:" << file_1 << " " << trueMH1 << std::endl;
  }
  if(strstr(file_2, "normal") != NULL)
  {
     trueMH2 = "normal";
     std::cout << "File 2:" << file_1 << " " << trueMH2 << std::endl;
  }
  else
  {
     trueMH2 = "inverted";
     std::cout << "File 2:" << file_1 << " " << trueMH2 << std::endl;
  }

  // Global scan contours
  
  TCanvas *ccontours = new TCanvas("ccontours","ccontours");
  ccontours->SetTopMargin(0.1);

  const char* xaxis_title = hchisq_valor[0]->GetXaxis()->GetTitle();
  if(strstr(xaxis_title, "Delta m^{2}") != NULL)
  {
     if(strcmp(trueMH1, trueMH2) == 0)
     {  // Same hierarchy
        if(strcmp(trueMH1, "normal") == 0)
        {  // Normal hierarchy
           hchisq_valor[0]->GetXaxis()->SetTitle("#Deltam^{2}_{32} (eV^{2}c^{-4})");
        }
        else
        {  // Inverted hierarchy
           hchisq_valor[0]->GetXaxis()->SetTitle("#Deltam^{2}_{13} (eV^{2}c^{-4})");
        }
     }
     else
     {  // Different hierarchies
        hchisq_valor[0]->GetXaxis()->SetTitle("#Deltam^{2}_{32} (NO), #Deltam^{2}_{13} (IO) (eV^{2}c^{-4})");
     }
  }

  const char* yaxis_title = hchisq_valor[0]->GetYaxis()->GetTitle();
  if(strstr(yaxis_title, "Delta m^{2}") != NULL)
  {
     if(strcmp(trueMH1, trueMH2) == 0)
     {  // Same hierarchy
        if(strcmp(trueMH1, "normal") == 0)
        {  // Normal hierarchy
           hchisq_valor[0]->GetYaxis()->SetTitle("#Deltam^{2}_{32} (eV^{2}c^{-4})");
        }
        else
        {  // Inverted hierarchy
           hchisq_valor[0]->GetYaxis()->SetTitle("#Deltam^{2}_{13} (eV^{2}c^{-4})");
        }
     }
     else
     {  // Different hierarchies
        hchisq_valor[0]->GetYaxis()->SetTitle("#Deltam^{2}_{32} (NO), #Deltam^{2}_{13} (IO) (eV^{2}c^{-4})");
     }
  }

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
  
  if(draw_reactors)
  {
     breactors->Draw("same");
  }

  
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
  if(draw_reactors) lcontours_2->AddEntry(breactors,"PDG 2016","f");

  lcontours->Draw();
  lcontours_2->Draw();

  //char* note_text = (doGlobalMin) ? "Global minimum" :"Fixed Mass Ordering";
  char* notetext = "Final systematics pending";
  double tx = h1->GetXaxis()->GetBinLowEdge(h1->GetXaxis()->GetNbins() * 0.5);
  double nx = h1->GetXaxis()->GetBinLowEdge(h1->GetXaxis()->GetNbins() / 10);
  double ty = Y2;

  TLatex* title = new TLatex(tx, ty * 1.01, "T2K Run 1-9 Preliminary");
  title->SetTextSize(0.05);
  title->Draw("same");
  TLatex* note = new TLatex(nx, ty * 1.01, "Final systematics pending");
  note->SetTextSize(0.025);
//  note->Draw("same");

  ccontours->SaveAs(Form("%s.pdf",tag));
  ccontours->SaveAs(Form("%s.png",tag));
  ccontours->SaveAs(Form("%s.eps",tag));
  ccontours->SaveAs(Form("%s.root",tag));

  // Plot 1D projections
  // x projection 
  TCanvas* canvas_x = new TCanvas("canvas_x", "canvas_x");
  TH1D* hchisq_x1 = (TH1D*)filein[0]->Get("hchisq_x");
  TH1D* hchisq_x2 = (TH1D*)filein[1]->Get("hchisq_x");
  hchisq_x1->SetDirectory(0);
  hchisq_x2->SetDirectory(0);

  valor::RootFinder* find_root = new valor::RootFinder();
  find_root->SetGraph(hchisq_x1);
  find_root->SetDchiSq(1 + hchisq_x1->GetMinimum());
  find_root->GetRoots();
  cout << "X NH: " << endl;
  find_root->PrintRoots();
  delete find_root;
  
  find_root = new valor::RootFinder();
  find_root->SetGraph(hchisq_x2);
  
  find_root->SetDchiSq(1 + hchisq_x2->GetMinimum());
  find_root->GetRoots();
  cout << "X IH: " << endl;
  find_root->PrintRoots();
  delete find_root;

  double min = hchisq_x1->GetMinimum() < hchisq_x2->GetMinimum() ? hchisq_x1->GetMinimum() : hchisq_x2->GetMinimum();
  minimise_hist(hchisq_x1, min);
  minimise_hist(hchisq_x2, min);

  hchisq_x1->GetXaxis()->SetTitle(xaxis_title);
  hchisq_x1->GetYaxis()->SetTitle("-2#Deltaln(L)");
  hchisq_x2->GetXaxis()->SetTitle(xaxis_title);
  hchisq_x2->GetYaxis()->SetTitle("-2#Deltaln(L)");

  double max1 = hchisq_x1->GetMaximum();
  double max2 = hchisq_x2->GetMaximum();
  double max = (max1 > max2)? max1 : max2;
  max = TMath::Ceil(max);
  if(max > 25)
  {  // No need to go past 5 sigma for the plot
     max = 25.0;
  }

  hchisq_x1->GetYaxis()->SetRangeUser(0, max);
  hchisq_x2->GetYaxis()->SetRangeUser(0, max);
  if(paper_style)
  {
     format_hist(hchisq_x1);
     format_hist(hchisq_x2);
  }
  hchisq_x1->SetLineColor(kBlack);
  hchisq_x1->SetLineWidth(1);
  hchisq_x1->SetLineStyle(1);
  hchisq_x2->SetLineColor(kRed);
  hchisq_x2->SetLineWidth(1);
  hchisq_x2->SetLineStyle(2);

  hchisq_x1->Draw("l");
  
  if(draw_reactors && strstr(xaxis_title, "theta_{13}") != NULL)
  {
    delete breactors;
    breactors = new TBox(sinsqth13_nom - sinsqth13_err, hchisq_x1->GetMinimum(), sinsqth13_nom + sinsqth13_err, max);
    breactors->SetFillColor(kYellow+1);
    breactors->Draw("same");
  }

  hchisq_x1->Draw("l same");
  hchisq_x2->Draw("l same");

  double sigmas[4] = {1.0, 4.0, 9.0, 16.0};
  for(int i = 0; i < 4; ++i)
  {
    if(max1 > sigmas[i] || max2 > sigmas[i])
    {
      TLine *lsigma = new TLine(hchisq_x1->GetXaxis()->GetBinLowEdge(1), sigmas[i], hchisq_x1->GetXaxis()->GetBinUpEdge(hchisq_x1->GetNbinsX()), sigmas[i]);
      lsigma->SetLineWidth(2);
      lsigma->SetLineColor(kBlue);
      lsigma->Draw("l same");
  
      TLatex *latLineSigma = new TLatex(hchisq_x1->GetXaxis()->GetBinLowEdge(5), sigmas[i] + 0.1, Form("%d#sigma", (i + 1)));
      latLineSigma->SetTextColor(kBlue);
      latLineSigma->SetTextSize(0.05);
      latLineSigma->Draw("same");
    }
  }
  TLegend *lcontours = new TLegend(0.45,0.77,0.60,0.9);
  lcontours->SetFillStyle(0);
  lcontours->SetBorderSize(0);
  lcontours->SetTextFont(132);
  lcontours->SetTextSize(0.05);
  lcontours->AddEntry(hchisq_x1, latexname[0], "l");
  lcontours->AddEntry(hchisq_x2, latexname[1], "l");
  if(draw_reactors && strstr(xaxis_title, "theta_{13}") != NULL)
  {
     lcontours->AddEntry(breactors, "PDG 2016", "f");
  }
  lcontours->Draw();
 
  double tx = hchisq_x1->GetXaxis()->GetBinLowEdge(hchisq_x1->GetXaxis()->GetNbins() * 0.5);
  double ty = max;

  TLatex* title = new TLatex(tx, ty * 1.01, "T2K Run 1-9 Preliminary");
  title->SetTextSize(0.05);
  title->Draw("same");

  int first = indexof(tag, '_');
  int second = indexof(tag + first + 1, '_') + first + 1;
  char* subtag = new char[strlen(tag) - (second - first) + 1];
  strncpy(subtag, tag, first);
  strncpy(subtag + first, tag + second, strlen(tag) - second);
  subtag[strlen(tag) - (second - first)] = '\0';
  canvas_x->SaveAs(Form("%s_proj.pdf", subtag));
  canvas_x->SaveAs(Form("%s_proj.png", subtag));
  canvas_x->SaveAs(Form("%s_proj.eps", subtag));
  canvas_x->SaveAs(Form("%s_proj.root", subtag));
  delete[] subtag;

  // y projection 
  TCanvas* canvas_y = new TCanvas("canvas_y", "canvas_y");
  TH1D* hchisq_y1 = (TH1D*)filein[0]->Get("hchisq_y");
  TH1D* hchisq_y2 = (TH1D*)filein[1]->Get("hchisq_y");
  hchisq_y1->SetDirectory(0);
  hchisq_y2->SetDirectory(0);

  find_root = new valor::RootFinder();
  find_root->SetGraph(hchisq_y1);
  find_root->SetDchiSq(1 + hchisq_y1->GetMinimum());
  find_root->GetRoots();
  cout << "Y NH: " << endl;
  find_root->PrintRoots();
  delete find_root;
  
  find_root = new valor::RootFinder();
  find_root->SetGraph(hchisq_y2);
  
  find_root->SetDchiSq(1 + hchisq_y2->GetMinimum());
  find_root->GetRoots();
  cout << "Y IH: " << endl;
  find_root->PrintRoots();
  delete find_root;

  min = hchisq_y1->GetMinimum() < hchisq_y2->GetMinimum() ? hchisq_y1->GetMinimum() : hchisq_y2->GetMinimum();
  minimise_hist(hchisq_y1, min);
  minimise_hist(hchisq_y2, min);

  hchisq_y1->GetXaxis()->SetTitle(yaxis_title);
  hchisq_y1->GetYaxis()->SetTitle("-2#Deltaln(L)");
  hchisq_y2->GetXaxis()->SetTitle(xaxis_title);
  hchisq_y2->GetYaxis()->SetTitle("-2#Deltaln(L)");

  max1 = hchisq_y1->GetMaximum();
  max2 = hchisq_y2->GetMaximum();
  max = (max1 > max2)? max1 : max2;
  max = TMath::Ceil(max);
  if(max > 25)
  {  // No need to go past 5 sigma for the plot
     max = 25.0;
  }

  hchisq_y1->GetYaxis()->SetRangeUser(0, max);
  hchisq_y2->GetYaxis()->SetRangeUser(0, max);
  if(paper_style)
  {
     format_hist(hchisq_y1);
     format_hist(hchisq_y2);
  }
  hchisq_y1->SetLineColor(kBlack);
  hchisq_y1->SetLineWidth(1);
  hchisq_y1->SetLineStyle(1);
  hchisq_y2->SetLineColor(kRed);
  hchisq_y2->SetLineWidth(1);
  hchisq_y2->SetLineStyle(2);

  hchisq_y1->Draw("l");
  
  if(draw_reactors && strstr(yaxis_title, "theta_{13}") != NULL)
  {
    delete breactors;
    breactors = new TBox(sinsqth13_nom - sinsqth13_err, hchisq_y1->GetMinimum(), sinsqth13_nom + sinsqth13_err, max);
    breactors->SetFillColor(kYellow+1);
    breactors->Draw("same");
  }

  hchisq_y1->Draw("l same");
  hchisq_y2->Draw("l same");

  for(int i = 0; i < 4; ++i)
  {
    if(max1 > sigmas[i] || max2 > sigmas[i])
    {
      TLine *lsigma = new TLine(hchisq_y1->GetXaxis()->GetBinLowEdge(1), sigmas[i], hchisq_y1->GetXaxis()->GetBinUpEdge(hchisq_y1->GetNbinsX()), sigmas[i]);
      lsigma->SetLineWidth(2);
      lsigma->SetLineColor(kBlue);
      lsigma->Draw("l same");
  
      TLatex *latLineSigma = new TLatex(hchisq_y1->GetXaxis()->GetBinLowEdge(5), sigmas[i] + 0.1, Form("%d#sigma", (i + 1)));
      latLineSigma->SetTextColor(kBlue);
      latLineSigma->SetTextSize(0.05);
      latLineSigma->Draw("same");
    }
  }
  TLegend *lcontours = new TLegend(0.45,0.77,0.60,0.9);
  lcontours->SetFillStyle(0);
  lcontours->SetBorderSize(0);
  lcontours->SetTextFont(132);
  lcontours->SetTextSize(0.05);
  lcontours->AddEntry(hchisq_y1, latexname[0], "l");
  lcontours->AddEntry(hchisq_y2, latexname[1], "l");
  if(draw_reactors && strstr(yaxis_title, "theta_{13}") != NULL)
  {
     lcontours->AddEntry(breactors, "PDG 2016", "f");
  }
  lcontours->Draw();
 
  tx = hchisq_y1->GetXaxis()->GetBinLowEdge(hchisq_y1->GetXaxis()->GetNbins() * 0.5);
  ty = max;

  TLatex* title = new TLatex(tx, ty * 1.01, "T2K Run 1-9 Preliminary");
  title->SetTextSize(0.05);
  title->Draw("same");

  first = indexof(tag, '_');
  subtag = new char[strlen(tag) - (first + 1) + 1];
  strncpy(subtag, tag + first + 1, strlen(tag) - (first + 1));
  subtag[strlen(tag) - (first + 1)] = '\0';
  canvas_y->SaveAs(Form("%s_proj.pdf", subtag));
  canvas_y->SaveAs(Form("%s_proj.png", subtag));
  canvas_y->SaveAs(Form("%s_proj.eps", subtag));
  canvas_y->SaveAs(Form("%s_proj.root", subtag));
  delete[] subtag;

  return;
}

int indexof(const char* str, const char chr)
{
   if(str == 0)
   {
      return -1;
   }
   else
   {
      int len = strlen(str);
      for(int i = 0; i < len; ++i)
      {
         if(str[i] == chr)
         {
            return i;
         }
      }
      return -1;
   }
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

void minimise_hist(TH1* h, double min)
{
   for(int i = 1; i <= h->GetNbinsX(); i++)
   {
      h->SetBinContent(i, h->GetBinContent(i) - min);
   }
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
