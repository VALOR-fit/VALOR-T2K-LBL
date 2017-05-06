#include "get_contour_2D.C"
#include "dcp_extract_contour.C"

void get_contour(const char* _filename = "",bool draw_asimov_input = false, bool swap_xy = 0)
{
  string filename = _filename;
  TFile f(Form("%s",filename.c_str()));
  TTree* t = f.Get("fitnt");
  
  bool isfit2d=0;
  
  t->SetBranchAddress("IsFit2D",&isfit2d);
  t->GetEntry(0);
  
  if(isfit2d) get_contour_2D(_filename,draw_asimov_input,swap_xy);
  else dcp_extract_contour(_filename);
  

}


