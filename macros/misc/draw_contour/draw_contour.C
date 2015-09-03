//
// Use this macro to draw 60% and 90% confidence regions from an experiment,
// macro assumes 1 experiment per root file and dumps the output
// into a root file and PDF.
//
// Give it the location of the root file, witout the .root at the end eg. $dir/toyA
//
// If the input contains either "ss23" or "cp13" in the name, the bf point will be  added

void draw_contour
(
 const char tag[200] = "fitmarg_sinsqth13Vsdcp_T2Konly_Run1to4_Data_NH",
 const char* filename = "../marginalization/joint_newtoys/fitmarg_sinsqth13Vsdcp_T2Konly_Run1to4_Data_NH-sinsqtheta[12_0.304,13_0.024,23_0.528],dcp_-1.601,dmsq_Normal[12_0.00007530,23_0.002509]eV2_0_0",
 
 //const char tag[200] = "fitmarg_sinsqth13Vsdcp_T2Konly_Run14_Asimov_C_NH",
 //const char* filename = "../marginalization/joint_newtoys/fitmarg_sinsqth13Vsdcp_T2Konly_Run14_Asimov_C_NH-sinsqtheta[12_0.304,13_0.042,23_0.528],dcp_1.910,dmsq_Normal[12_0.00007530,23_0.002509]eV2_0_0",

 // const char tag[200] = "fitmarg_sinsqth13Vsdcp_T2Kreact_DataRun14_NH",
 // const char* filename = "../marginalization/joint_newtoys/fitmarg_sinsqth13Vsdcp_T2Kreact_Run1to4_Data_NH-sinsqtheta[12_0.304,13_0.024,23_0.528],dcp_0.000,dmsq_Normal[12_0.00007530,23_0.002509]eV2_0_0",

 //const char tag[200] = "fitmarg_sinsqth13Vsdcp_T2Kreact_DataRun14_IH",
 //const char* filename = "../marginalization/joint_newtoys/fitmarg_sinsqth13Vsdcp_T2Kreact_Run1to4_Data_IH-sinsqtheta[12_0.304,13_0.024,23_0.528],dcp_0.000,dmsq_Normal[12_0.00007530,23_0.002509]eV2_0_0",

 const char* opt = "cp13",
 const bool SavePDF = false
 )
{
  //   char* filename = "toyA_joint";
  
  TCanvas* c = new TCanvas("c","c",0,0,900,600);
  TFile f(Form("%s.root",filename));
  TTree* t = f.Get("fitnt");
  
  TH2D* h1;
  TH2D* h2;
  TH2D* h3;
  TH2D* h4;
  valor::FitParamList* fp = 0;
  valor::PhysParamList* ip = 0;
  
  t->SetBranchAddress("hchisq2D",&h1);
  t->SetBranchAddress("best_fit_params",&fp);
  t->SetBranchAddress("input_params",&ip);    
  t->GetEntry(0);
  
  f.Close();
  
  
  string file2 = filename;
  string option = opt;
  double input_x,input_y;
  if(file2.find("ss23")!= string::npos || option.find("ss23")!= string::npos )
    {
      input_x = ip->Value(valor::conv::kPhypMixing23);
      input_y = ip->Value(valor::conv::kPhypDm23sq);
    }
  else if ( file2.find("cp13")!= string::npos || option.find("cp13")!= string::npos)
    {
      input_x = ip->Value(valor::conv::kPhypMixing13);
      input_y = ip->Value(valor::conv::kPhypCPV);
    }
  else{
    cout << "\n\n\tERROR: Could not find which parameters are being used\n\n" <<endl;
    input_x = 0;
    input_y = 0;
  }
  
  cout << "x: " << input_x << endl;
  cout << "y: " << input_y << endl;
  
  
  h3 = (TH2D*)h1->Clone();
  h2 = (TH2D*)h1->Clone();
  h4 = (TH2D*)h1->Clone();
  double chisq_min = h1->GetMinimum();

  cout << "chisq_min = " << chisq_min << endl;

  h1->SetContour(1);
  h3->SetContour(1);
  
  h1->SetContourLevel(0,chisq_min+2.3);
  h3->SetContourLevel(0,chisq_min+4.61);
  h1->SetLineWidth(2);
  h3->SetLineWidth(2);
  h1->SetLineColor(kRed);
  h3->SetLineColor(kRed);
  h1->SetLineStyle(2);
  
  int minbin = h1->GetMinimumBin();
  int binx;
  int biny;
  h1->GetBinXYZ(minbin,binx,biny,0);
  
  double xmin = h1->GetXaxis()->GetBinCenter(binx);
  double ymin = h1->GetYaxis()->GetBinCenter(biny);
  
  cout << "best-fit x = " << xmin << endl;
  cout << "best-fit y = " << ymin << endl;

  TMarker * bf = new TMarker(xmin,ymin,kFullStar);
  bf->SetMarkerColor(kRed);
  bf->SetMarkerSize(1.5);
  
  TMarker * bf2 = new TMarker(input_x,input_y,kFullStar);
  bf2->SetMarkerColor(kViolet);
  bf2->SetMarkerSize(1.5);
  
  TLegend* leg = new TLegend(0.7,0.8,0.9,0.98);
  leg->AddEntry(bf,"Best fit point","p");
  leg->AddEntry(bf2,"True asimov point","p");
  leg->AddEntry(h1,"68%","l");
  leg->AddEntry(h3,"90%","l");
  h2->Draw("colz");
  h1->Draw("cont3same");
  h3->Draw("cont3same");
  bf->Draw();
  bf2->Draw();
  
  leg->Draw();
  
  if(SavePDF) c->SaveAs(Form("%s.pdf",tag));
  
  TFile* f2 = new TFile(Form("%s_contour.root",tag),"recreate");
  
  save_contour(f2,h1,h3);
  h2->Write("chisq_map_2D");
  
  bf->Write("best_fit_point_marker");
  f2->Close();
  
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
