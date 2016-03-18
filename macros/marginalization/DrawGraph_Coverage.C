
#include "Build.h"
#include "AnalysisOptions.h"
#include <vector>

using namespace std;
using namespace valor;
using namespace valor::conv;

void DrawGraph_Coverage(){

  TGraph *gcritval_90CL = new TGraph();
  gcritval_90CL->SetPoint(0,-1.601,1.78463);
  gcritval_90CL->SetPoint(1,-1.,1.93531);
  gcritval_90CL->SetPoint(2,0.,2.18557);
  gcritval_90CL->SetPoint(3,+2.,1.94271);
  gcritval_90CL->GetYaxis()->SetRangeUser(0,3);
  gcritval_90CL->GetYaxis()->SetTitle("#Delta #chi^{2}_{crit}");
  gcritval_90CL->GetXaxis()->SetTitle("#delta_{CP}"); 
  TCanvas *ccritval_90CL = new TCanvas("ccritval_90CL","ccritval_90CL");
  gcritval_90CL->Draw("apl");

  TGraph *gcov_90CL = new TGraph();
  gcov_90CL->SetPoint(0,-1.601,0.953);
  gcov_90CL->SetPoint(1,-1.,0.948737);
  gcov_90CL->SetPoint(2,0.,0.9499);
  gcov_90CL->SetPoint(3,+2.,0.949263);
  gcov_90CL->GetYaxis()->SetRangeUser(0.5,1);
  gcov_90CL->GetYaxis()->SetTitle("Coverage");
  gcov_90CL->GetXaxis()->SetTitle("#delta_{CP}"); 

  TGraph *gcov_1sigma = new TGraph();
  gcov_1sigma->SetPoint(0,-1.601,0.796);
  gcov_1sigma->SetPoint(1,-1.,0.778211);
  gcov_1sigma->SetPoint(2,0.,0.6165);
  gcov_1sigma->SetPoint(3,+2.,0.769368);
  gcov_1sigma->GetYaxis()->SetRangeUser(0.5,1);
  gcov_1sigma->GetYaxis()->SetTitle("Coverage");
  gcov_1sigma->GetXaxis()->SetTitle("#delta_{CP}"); 
  TCanvas *ccov = new TCanvas("ccov","ccov");
  gcov_90CL->SetLineColor(kRed);
  gcov_90CL->Draw("apl");
  gcov_1sigma->Draw("pl same");

  return;
}
