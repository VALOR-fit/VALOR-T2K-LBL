//
// ------------------------------------------------------------------------------------------------
// Plot predicted SuperK reconstructed neutrino energy spectrum for various oscillation hypotheses.
// The spectra are shown for the input Dm^2_{23} and mixing_{23}.
// The macro produces the 2-flavour hypothesis spectrum and spectra for three different 3-flavour 
// hypotheses. In each 3-flavour hypothesis only one of Dm^2_{12}, mixing_{12}, 
// mixing_{13} or delta_cp is varied while the other parameters remain fixed.
// The 3-flavour / 2-flavour spectrum ratio is also plotted.
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

#include <iomanip.h>

void test_sk_ereco_spectrum_twk_osc_non23sector
(
  int          varied_osc_param   = 0,                              /* 0: vary dmsq12, 1: vary theta12, 2: vary theta13, 3: vary deltacp */
  double       input_mixing23     = 1.00,                           /* */
  double       input_dm23sq       = 2.32E-3,                        /* */
  const char * config_file        = "analysis_configuration.xml"    /* */
)
{
  using namespace std;
  using namespace t2kana;

  // Analysis configuration
  const char * config  = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);

  // Global options
  Opt::Instance()->gOscMixingParamChoice = conv::kNuMixSinsq2Theta;
  Opt::Instance()->gMassHierarchyAssumed = conv::kNmhNormal;

  // Define oscillation hypotheses
  const  int ncases = 4; // one 2-flavour case (0) + three 3-flavour cases (1,2,3)
  const  double pi = TMath::Pi();

  double dm12sq_fix   = 7.58E-5;
  double mixing12_fix = 0.86;
  double mixing13_fix = 0.10;
  double deltacp_fix  = pi/2;
  double dm12sq_var   [ncases-1] = { 6.50E-5, 7.58E-5, 8.50E-5 };
  double mixing12_var [ncases-1] = { 0.700,   0.860,   1.000   };
  double mixing13_var [ncases-1] = { 0.000,   0.100,   0.200   };
  double deltacp_var  [ncases-1] = { 0.000,   pi/2.,   pi      };

  double dm12sq   [ncases] = {-1};
  double mixing12 [ncases] = {-1}; 
  double mixing13 [ncases] = {-1};
  double deltacp  [ncases] = {-1};
  for(int icase = 1; icase < ncases; icase++) { // case 0 is the 2-flavour case and it is skipped
    dm12sq   [icase] = (varied_osc_param == 0) ? dm12sq_var   [icase-1] : dm12sq_fix;
    mixing12 [icase] = (varied_osc_param == 1) ? mixing12_var [icase-1] : mixing12_fix;
    mixing13 [icase] = (varied_osc_param == 2) ? mixing13_var [icase-1] : mixing13_fix;
    deltacp  [icase] = (varied_osc_param == 3) ? deltacp_var  [icase-1] : deltacp_fix;
  }

  // Histogram drawing options
  int hwidth [ncases] = { 2,       1,      1,       1      };
  int hcolor [ncases] = { kBlack,  kRed,   kGreen,  kBlue  };
  int hstyle [ncases] = { kDashed, kSolid, kSolid,  kSolid };

  // Create & config superk spectrum predictor
  SuperkExpectation skexp;
  skexp.Read(config);
  
  // Calculate spectra
  TH1D * hEreco[ncases] = {0};
  for(int icase = 0; icase < ncases; icase++) {
    conv::PhysHypothesis_t phid = 
       (icase == 0) ? (conv::kHypNuMuDisapp2F) : (conv::kHypOsc3F);
    // physics hypothesis 
    PhysHypothesis ph(phid);
    // update osc params
    ph.SetCurrValue(conv::kPhypMixing23,  input_mixing23);
    ph.SetCurrValue(conv::kPhypDm23sq,    input_dm23sq       );
    if(icase > 0) {
      ph.SetCurrValue(conv::kPhypMixing12,  mixing12[icase]);
      ph.SetCurrValue(conv::kPhypMixing13,  mixing13[icase]);
      ph.SetCurrValue(conv::kPhypDeltaCP,   deltacp[icase] );
      ph.SetCurrValue(conv::kPhypDm12sq,    dm12sq[icase]  );
    }
    // generate predictions  
    skexp.Generate(ph);
    // store spectrum
    hEreco[icase] = new TH1D(*skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike));
    hEreco[icase]->SetDirectory(0);
    hEreco[icase]->SetLineWidth(hwidth[icase]);
    hEreco[icase]->SetLineColor(hcolor[icase]);
    hEreco[icase]->SetLineStyle(hstyle[icase]);
  }

  // Form ratios of 3-flavour hypotheses to the 2-flavou hypothesis
  TH1D * hErecoRatio[ncases] = {0};
  for(int icase = 1; icase < ncases; icase++) {
    hErecoRatio[icase] = new TH1D(*hEreco[icase]);
    hErecoRatio[icase]->Divide(hEreco[0]);
    hErecoRatio[icase]->SetDirectory(0);
    hErecoRatio[icase]->SetLineWidth(hwidth[icase]);
    hErecoRatio[icase]->SetLineColor(hcolor[icase]);
    hErecoRatio[icase]->SetLineStyle(hstyle[icase]);
  }

  // Draw
  const char * mixing_param = 0;
  if( Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsq2Theta ) { mixing_param = "sin^{2}2#theta"; }
  if( Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsqTheta  ) { mixing_param = "sin^{2}#theta";  }

  gStyle->SetOptStat(0);
  TCanvas * c = new TCanvas("c","",20,20,700,500);
  TLegend * legend = new TLegend(0.3, 0.7, 0.85, 0.9);
  legend->SetFillColor(0);
  legend->SetBorderSize(0);
  legend->SetHeader(Form("|#Deltam_{23}^{2}| = %.5f eV^{2}/c^{4}, %s_{23} = %.2f",input_dm23sq,mixing_param,input_mixing23));
  for(int icase = 0; icase < ncases; icase++) {
     const char * opt = 
        (icase == 0) ? 
           "" : 
           "same";
     const char * label = 
        (icase == 0) ? 
           "2-flavour" : 
           Form("3-flavour: #Deltam_{12}^{2} = %.5f eV^{2}/c^{4}, %s_{12} = %.2f, %s_{13} = %.2f, #delta_{CP} = %.3f", 
           dm12sq[icase], mixing_param, mixing12[icase], mixing_param, mixing13[icase], deltacp[icase]);
     if(icase==0) {
       hEreco[icase]->GetXaxis()->SetTitle("E_{reco} (GeV)");
       hEreco[icase]->GetYaxis()->SetTitle("Number of events");
     }
     //hEreco[icase]->GetXaxis()->SetRangeUser(0.,2.);
     //hEreco[icase]->GetYaxis()->SetRangeUser(0.,2.);
     hEreco[icase]->Draw(opt);
     legend->AddEntry(hEreco[icase], label, "L");
     cout << " *** " << label << " : " << endl;
     cout << "NumEvents = " 
          <<  setw(10) << setfill(' ') << setprecision(6) << hEreco[icase]->Integral() << endl;
     const int nemax = 7;
     double emax[nemax] = { 0.50, 0.75, 1.00, 1.50, 2.00, 4.00, 10.00 };
     for (int iemax = 0; iemax < nemax; iemax++) {
          cout << "NumEvents (< " 
               << setw(6)  << setfill(' ') << setprecision(3) 
               << emax[iemax] << " GeV): " 
               << setw(10) << setfill(' ') << setprecision(6)
               << hEreco[icase]->Integral(1,hEreco[icase]->FindBin(emax[iemax])) << endl;
     }
  }
  legend->Draw();
  c->Update();
  TCanvas * cr = new TCanvas("cr","",120,120,600,400);
  TLegend * legendr = new TLegend(0.3, 0.7, 0.85, 0.9);
  legendr->SetFillColor(0);
  legendr->SetBorderSize(0);
  legendr->SetHeader(Form("|#Deltam_{23}^{2}| = %.5f eV^{2}/c^{4}, %s_{23} = %.2f",input_dm23sq,mixing_param,input_mixing23));
  for(int icase = 1; icase < ncases; icase++) {
     const char * opt = 
        (icase == 1) ? 
           "" : 
           "same";
     const char * label = 
           Form("#Deltam_{12}^{2} = %.5f eV^{2}/c^{4}, %s_{12} = %.2f, %s_{13} = %.2f, #delta_{CP} = %.3f", 
           dm12sq[icase], mixing_param, mixing12[icase], mixing_param, mixing13[icase], deltacp[icase]);
     if(icase==1) {
       hErecoRatio[icase]->GetXaxis()->SetTitle("E_{reco} (GeV)");
       hErecoRatio[icase]->GetYaxis()->SetTitle("3-flavour / 2-flavour hypothesis ratio");
     }
     //hErecoRatio[icase]->GetXaxis()->SetRangeUser(0.0,2.0);
     //hErecoRatio[icase]->GetYaxis()->SetRangeUser(0.7,1.3);
     hErecoRatio[icase]->Draw(opt);
     legendr->AddEntry(hErecoRatio[icase], label, "L");
  }
  legendr->Draw();
  cr->Update();
  
}

