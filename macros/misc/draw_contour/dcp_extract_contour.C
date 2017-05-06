
using namespace std;
using namespace valor;
using namespace valor::conv;


dcp_extract_contour(const char* _filename)
{
  //   char* filename = "toyA_joint";
  string filename = _filename;
  filename.replace(filename.find(".root"), sizeof(".root")-1, "");
  
  TCanvas* c = new TCanvas("c","c",0,0,900,600);
  c->SetGrid();

  TFile f(Form("%s.root",filename.c_str()));
  TTree* t = f.Get("fitnt");
  TH1D* h1;
  TH1D* h2;
  
  
  string *trueparname0;
  string *trueMH;

  t->SetBranchAddress("hdeltachisq1D",&h1);
  t->SetBranchAddress("hchisq1D",&h2);
  t->SetBranchAddress("trueparname0",   &trueparname0   );
  t->SetBranchAddress("trueMH",   &trueMH   );

  t->GetEntry(0);
  valor::conv::NuMassHierarchy_t hierarchy = util::HierarchyFromString(*trueMH);
  
  valor::conv::PhysParamId_t p = valor::util::PhysParamIdFromString(*trueparname0);
  string x_name = valor::util::AsStringOfficialAxis(p);//,hierarchy);
  h1->GetXaxis()->SetTitle(x_name.c_str());
  h2->GetXaxis()->SetTitle(x_name.c_str());

  
  format_hist(h1);
  format_hist(h2);
  
  h1->Draw();
  
  c->SaveAs(Form("%s_contour.pdf",filename.c_str()));

  TFile f2(Form("%s_contour.root",filename.c_str()),"UPDATE");
  h1->Write("h_deltachisq");
  h2->Write("h_chisq");
  f2.Close();

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
  
}
