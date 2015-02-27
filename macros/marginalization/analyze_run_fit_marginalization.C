
using namespace std;
using namespace valor;
using namespace valor::conv;

void analyze_run_fit_marginalization
(
 char infilename[200] = "ChisqMarg_0.root",
 bool Fit2D           = true
 )
{
  TFile *filein = new TFile(infilename,"READ");
  
  TH1F *hChisqMarg_1D = (TH1F*) filein->Get("hchisq_marg_1D");
  TH2F *hChisqMarg_2D = (TH2F*) filein->Get("hchisq_marg_2D");
  
  TCanvas *cChisqMarg_1D = new TCanvas("cChisqMarg_1D","cChisqMarg_1D");
  hChisqMarg_1D->Draw("");
  
  TCanvas *cChisqMarg_2D = new TCanvas("cChisqMarg_2D","cChisqMarg_2D");
  hChisqMarg_2D->Draw("colz");

  // Print out marginalized chisq
  cout << endl;
  cout << "Marginalized Chisq:" << endl;
  cout << endl;
  int nbinsX = hChisqMarg_2D->GetNbinsX();
  int nbinsY = hChisqMarg_2D->GetNbinsY();
  for(int ibinX=0;ibinX<nbinsX;ibinX++){
    for(int ibinY=0;ibinY<nbinsY;ibinY++){
      cout << hChisqMarg_2D->GetBinContent(ibinX+1,ibinY+1) << ", ";
    }
    cout << endl;
  }
  cout << endl;

  
  // Find best fit

  int minbin       = -999;
  int minbinX      = -999;
  int minbinY      = -999;
  int minbinZ      = -999;
  double chisq_bf  = -999;
  double param0_bf = -999;
  double param1_bf = -999;
  
  cout << endl;
  cout << "Find minimum of marginalized chisq" << endl;
  cout << endl;

  if(!Fit2D){
    chisq_bf_1D = hchisq_marg_1D->GetMinimum();
    minbin = hchisq_marg_1D->GetMinimumBin();
    param0_bf = hchisq_marg_1D->GetBinLowEdge(minbin);
    param1_bf = -999.;
  }
  else{
    chisq_bf = hchisq_marg_2D->GetMinimum();
    hchisq_marg_2D->GetMinimumBin(minbinX,minbinY,minbinZ);
    param0_bf = hchisq_marg_2D->GetXaxis()->GetBinLowEdge(minbinX);
    param1_bf = hchisq_marg_2D->GetYaxis()->GetBinLowEdge(minbinY);
  }

  // Print out best fit values
    cout << "Best-fit marginalized Chisq = " << chisq_bf << endl;
  if(!Fit2D){
    cout << endl;
    cout << "Best-fit param 0 = " << param0_bf << endl;
    cout << endl;
  }
  else{
    cout << endl;
    cout << "Best-fit param 0 = " << param0_bf << endl;
    cout << "Best-fit param 1 = " << param1_bf << endl;
    cout << endl;
  }
  

  return;
}
