
//
// Modify covariance element of parameter sigma_nuenumu from 0.0004 (err=0.02) to 0.0009 (err=0.03)
//

#include <iostream>
#include <fstream>

#include "TROOT.h"
#include "TTree.h"
#include "TSystem.h"
#include "TH1D.h"
#include "TStopwatch.h"
#include "TString.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TVectorD.h"
#include "TRandom3.h"
#include "TSystem.h"
#include "TMath.h"
#include "TMatrixD.h"
#include "TAxis.h"
#include "TXMLEngine.h"

using namespace std;
using namespace valor;
using namespace valor::conv;

void MakeCovMatrix_normnue0p03
(
 const char *filenamein = "/Users/davidesgalaberna/Desktop/ETH_phd/T2K/ANA/VALOR/valor/softw/t2xk/trunk/inputs/official/2015/banff/postfit_banff_2015_data_20150417_valor_order_all_plus_scc.root",
 //const char *filenameout = "postfit_banff_2015_data_20150417_valor_order_rhc_plus_scc_normnue0p03.root",
 //int element = 38, // fhc-only or rhc-only cov matrix
 const char *filenameout = "postfit_banff_2015_data_20150417_valor_order_all_plus_scc_normnue0p03.root",
 int element = 63, // joint fhc+rhc cov matrix
)
{

  TFile *filein = new TFile(filenamein,"READ");
  
  TObjArray *param_list = (TObjArray*) filein->Get("param_list");
  param_list->Print();
  

  TVectorT<double> *prefit_params = (TVectorT<double>*)	filein->Get("prefit_params"); 
  TVectorT<double> *postfit_params = (TVectorT<double>*) filein->Get("postfit_params"); 

  TMatrixT<double> *prefit_cov = (TMatrixT<double>*) filein->Get("prefit_cov");
  cout << prefit_cov(element,element) << endl;
  prefit_cov(element,element) = 0.0009;
  cout << prefit_cov(element,element) << endl;
  //postfit_cov_old->Print();

  TMatrixT<double> *postfit_cov = (TMatrixT<double>*) filein->Get("postfit_cov");
  cout << postfit_cov(element,element) << endl;
  postfit_cov(element,element) = 0.0009;
  cout << postfit_cov(element,element) << endl;
  //postfit_cov_old->Print();

  TFile *fileout = new TFile(filenameout,"RECREATE");
   
  param_list->Write("param_list",TObject::kSingleKey);
  prefit_params->Write("prefit_params");
  prefit_cov->Write("prefit_cov");
  postfit_params->Write("postfit_params");
  postfit_cov->Write("postfit_cov");
  
  fileout->Close();
  
  
}
