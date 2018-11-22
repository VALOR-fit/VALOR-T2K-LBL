
void save_contour(TH2D* h, const char* label)
{
  TCanvas* c2 = new TCanvas(TString::Format("c2%s",label),"c2",0,0,900,600);
  h->Draw("CONT LIST");
  c2->Update();
  TObjArray *contours = 0;
  contours = (TObjArray *)gROOT->GetListOfSpecials()->FindObject("contours");
  TList* contour_list = (TList*)contours->At(0);
  TGraph *g;
  int it=0;
  
  TIter next(contour_list);
  while ((g = (TGraph*)next()))
  {
    TString name;
    name = TString::Format("%s_%d",label,it);
    g->SetLineWidth(2);
    g->SetLineStyle(1);
    g->SetLineColor(kRed);
    g->Write(name);
    it++;
  }

  delete c2;
}
