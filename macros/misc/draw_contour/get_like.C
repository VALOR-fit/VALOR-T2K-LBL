//
// Use this macro to draw 60% and 90% credible intervals from an experiment,
// macro assumes 1 experiment per root file and dumps the output
// into a root file and PDF.
//
// Give it the location of the root file, witout the .root at the end eg. $dir/toyA
//

void get_like
(
 const char tag[200] = "fitmarg_sinsqth13Vsdcp_T2Konly_Run1to4_Data_NH",
 const char* filename = "../marginalization/joint_newtoys/fitmarg_sinsqth13Vsdcp_T2Konly_Run1to4_Data_NH-sinsqtheta[12_0.304,13_0.024,23_0.528],dcp_-1.601,dmsq_Normal[12_0.00007530,23_0.002509]eV2_0_0", 
 const char* opt = "cp13",
 const bool SavePDF = false
 )
{
  //    char* filename = "toyA_joint";
  TFile f(Form("%s.root",filename));
  TTree* t = f.Get("fitnt");
  
  TH2D* h1;
  TH2D* h3;
  TH2D* h2;
  TH2D* h4;
  TH2D* h5;
  
  valor::FitParamList* fp = 0;
  valor::PhysParamList* ip = 0;
  
  t->SetBranchAddress("hchisq2D",&h1);
  t->SetBranchAddress("best_fit_params",&fp);
  t->SetBranchAddress("input_params",&ip);
  
  t->GetEntry(0);
  
  f.Close();
  
  h3 = (TH2D*)h1->Clone();
  h2 = (TH2D*)h1->Clone();
  h4 = (TH2D*)h1->Clone();
  h5 = (TH2D*)h1->Clone();
  
  TH1D* hproj_ss23 = 0;
  
  for(int i = 1;i<((h1->GetNbinsX() +1));i++)
    {
      for(int j = 1;j<(h1->GetNbinsY()+1);j++)
        {
	  h3->SetBinContent(i,j,TMath::Exp(30-h1->GetBinContent(i,j)/2));
	  h2->SetBinContent(i,j,TMath::Exp(30-h1->GetBinContent(i,j)/2));
	  h4->SetBinContent(i,j,TMath::Exp(30-h1->GetBinContent(i,j)/2));
	  h5->SetBinContent(i,j,TMath::Exp(30-h1->GetBinContent(i,j)/2));
        }
    }
  
  string file2 = filename;
  
  double input_x,input_y;
  if(file2.find("ss23")!= string::npos)
    {
      input_x = ip->Value(valor::conv::kPhypMixing23);
      input_y = ip->Value(valor::conv::kPhypDm23sq);
    }
  else if (file2.find("cp13")!= string::npos)
    {
      input_x = ip->Value(valor::conv::kPhypMixing13);
      input_y = ip->Value(valor::conv::kPhypCPV);
    }
  else{
    cout << "ERROR: Could not find which parameters are being used" <<endl;
  }
  
  
  hproj_ss23 = h2->ProjectionX();
  
  double total_prob = h3->Integral();
  cout << total_prob << endl;
  double count_prob = 0;
  int bin;
  double cred_68 = 0;
  double cred_90 = 0;
  
  while (count_prob < total_prob*0.68)
    {
      bin =  h3->GetMaximumBin();
      count_prob += h3->GetBinContent(bin);
      cred_68 = h3->GetBinContent(bin);
      h3->SetBinContent(bin,-1);
    }
  
  while (count_prob < total_prob*0.9)
    {
      bin =  h3->GetMaximumBin();
      count_prob += h3->GetBinContent(bin);
      cred_90 = h3->GetBinContent(bin);
      h3->SetBinContent(bin,-1);
    }
  
  h4->SetContour(1);
  h5->SetContour(1);
  
  count_prob = 0;
  
  
  cout << "cred_68:\t" << cred_68 << endl;
  cout << "cred_90:\t" << cred_90 << endl;
  
  h4->Draw("colz");
  
  h4->SetContourLevel(0,cred_68);
  h5->SetContourLevel(0,cred_90);
  h4->SetLineWidth(2);
  h5->SetLineWidth(2);
  h4->SetLineColor(kRed);
  h5->SetLineColor(kRed);
  h4->SetLineStyle(2);
  
  TCanvas* c6 = new TCanvas("c6","c6",0,0,900,600);
  
  h2->Draw("colz");
  h4->Draw("cont3same");
  h5->Draw("cont3same");
  
  TLegend* leg = new TLegend(0.7,0.8,0.9,0.98);
  leg->AddEntry(h4,"68%","l");
  leg->AddEntry(h5,"90%","l");
  leg->Draw();
  
  TString name =  Form("%s_credible.pdf",filename);
  cout << name << endl;
  c6->Print(name.Data());
  
  TFile* f2 = new TFile(Form("%s_contour_credible.root",tag),"recreate");

  save_contour(f2,h4,h5);
  hproj_ss23->Write("ss23");
  h2->Write("marginal_likelihood_2d");
  
  f2->Close();
  
  // h3->Draw("colz");
  
}

void save_contour(TFile* f2, TH2D* h1=0, TH2D* h3=0)
{
  TCanvas* c2 = new TCanvas("c2","c",0,0,900,600);
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
      char name[10];
      sprintf(name, "g90_valor_%d",it);
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
      char name[10];
      sprintf(name, "g68_valor_%d",it);
      g->SetLineWidth(2);
      g->SetLineStyle(2);
      g->SetLineColor(kRed);
      g->Write(name);
      it++;
    }
  
  delete c2;
  delete contours;
  delete list68;
  delete list90;
  delete g;
  
}
