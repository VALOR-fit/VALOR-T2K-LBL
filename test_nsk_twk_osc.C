//
// ------------------------------------------------------------------------------------------------
// Plot number of events predicted at SuperK (NSK) vs oscillation params
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

#include "Conventions.h"

void test_nsk_twk_osc()
{
  gSystem->Load("libGpad");  
  gSystem->Load("libPhysics");  
  gSystem->Load("libMatrix");   
  gSystem->Load("libHist");     
  gSystem->Load("libGraf");     
  gSystem->Load("libTree");     
  gSystem->Load("libRIO");   
  gSystem->Load("libXMLIO");
  gSystem->Load("libMinuit");   

  gSystem->Load("$T2KANA/lib/libT2KAnaYr1NuMuOscAnalysis.so");

  using namespace std;
  using namespace t2kana;
  using namespace t2kana::yr1numu;

  double sinsq2theta23_min  = 0.0;
  double sinsq2theta23_max  = 1.0;
  double sinsq2theta23_step = 0.02;
  double dm23sq_min         = 1E-4;
  double dm23sq_max         = 8E-3;
  double dm23sq_step        = 1E-4;

 // make sure that min and max are at the bin centre
 int sinsq2theta23_nbins = 1 + (int) ((sinsq2theta23_max - sinsq2theta23_min) / sinsq2theta23_step);
 int dm23sq_nbins        = 1 + (int) ((dm23sq_max        - dm23sq_min       ) / dm23sq_step       );

 TH2D * nsk_vs_osc_params = 
   new TH2D("nsk_vs_osc_params", "",
       sinsq2theta23_nbins, 
       sinsq2theta23_min - 0.5*sinsq2theta23_step, 
       sinsq2theta23_max + 0.5*sinsq2theta23_step,
       dm23sq_nbins, 
       dm23sq_min - 0.5*dm23sq_step, 
       dm23sq_max + 0.5*dm23sq_step);

  // superk spectrum predictor
  SuperkExpectation skexp;
  skexp.Read("inputs/analysis_configuration.xml");

  // fit params
  PhysParamList * fparm = PhysParamList::NumuDisapFitNoSyst();

  // generate predictions
  for(int i = 1; i <= nsk_vs_osc_params->GetNbinsX(); i++) {
    for(int j = 1; j <= nsk_vs_osc_params->GetNbinsY(); j++) {
       double sinsq2theta23 = nsk_vs_osc_params->GetXaxis()->GetBinCenter(i);
       double dm23sq        = nsk_vs_osc_params->GetYaxis()->GetBinCenter(j);
       fparm->Update(t2kana::yr1numu::conv::kPhypSinsq2Theta23, sinsq2theta23);
       fparm->Update(t2kana::yr1numu::conv::kPhypDm23sq,        dm23sq);

       // generate predictions
       skexp.Generate(*fparm);

       // integrate Ereco distribution
       double nsk = skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike)->Integral();

       nsk_vs_osc_params->SetBinContent(i,j,nsk);
    }//j
  }//i

  nsk_vs_osc_params->Draw("colz");
}
