
void save_contour(TH2D* h1=0, TH2D* h3=0,const char* label = "")
{
  TCanvas* c2 = new TCanvas(TString::Format("c2%s",label),"c",0,0,900,600);
  h3->Draw("CONT LIST");
  c2->Update();
  TObjArray *contours = 0;
  contours = (TObjArray *)gROOT->GetListOfSpecials()->FindObject("contours");
  TList *list90 = (TList*)contours->At(0);
  TGraph *g;
  int it=0;
  
  TIter next(list90);
  while ((g = (TGraph*)next()))
  {
    TString name;
    name = TString::Format("g90_valor%s_%d",label,it);
    g->SetLineWidth(2);
    g->SetLineStyle(1);
    g->SetLineColor(kRed);
    g->Write(name);
    it++;
  }

  c2->Clear();
  h1->Draw("CONT LIST");
  c2->Update();
  
  contours = (TObjArray*)gROOT->GetListOfSpecials()->FindObject("contours");
  TList *list68 = (TList*)contours->At(0);
  TIter next68(list68);
  it=0;
  while ((g = (TGraph*)next68()))
  {
    TString name;
    name = TString::Format("g68_valor%s_%d",label,it);
    g->SetLineWidth(2);
    g->SetLineStyle(2);
    g->SetLineColor(kRed);
    g->Write(name);
    it++;
  }
  
  delete c2;
  c2=0;
  contours =0;
  list68 =0;
  list90 =0;
  g=0;
  
}
