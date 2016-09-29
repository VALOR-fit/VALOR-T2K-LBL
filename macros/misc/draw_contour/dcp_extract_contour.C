
using namespace std;
using namespace valor;
using namespace valor::conv;


dcp_extract_contour(const char* _filename)
{
  //   char* filename = "toyA_joint";
  string filename = _filename;
  filename.replace(filename.find(".root"), sizeof(".root")-1, "");
  
  TCanvas* c = new TCanvas("c","c",0,0,900,600);
  TFile f(Form("%s.root",filename.c_str()));
  TTree* t = f.Get("fitnt");
  TH1D* h1;
  TH1D* h2;
  
  
  
  t->SetBranchAddress("hdeltachisq1D",&h1);
  t->SetBranchAddress("hchisq1D",&h2);
  t->GetEntry(0);
  
  h1->Draw();
  
  c->SaveAs(Form("%s_contour.pdf",filename.c_str()));

  TFile f2(Form("%s_contour.root",filename.c_str()),"UPDATE");
  h1->Write("h_deltachisq");
  h2->Write("h_chisq");
  f2.Close();

}
