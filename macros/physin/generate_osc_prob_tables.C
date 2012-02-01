//
// ------------------------------------------------------------------------------------------------
// Generate tables with oscillation probabilites
// to check with the tables used for official analysis
// (generated with Prob3++ package)
//
// Lorena
// ------------------------------------------------------------------------------------------------
//
//

#include <iomanip.h>
//#include <NuOsc3Flavour.h>
//#include "Conventions.h"

void generate_osc_prob_tables
(
  const char * hypothesis                 = "3f_osc",                       /* options: "2f_numu_disap", "3f_osc" */
  double       input_mixing12             = 0.8704,                           /* default values as for official analysis*/
  double       input_mixing13             = 0.10,                           /* */ 
  double       input_mixing23             = 1.0,                           /* */
  double       input_deltacp              = 0.00,                           /* */
  double       input_dm12sq               = 7.6E-5,                         /* */
  double       input_dm23sq               = 2.4E-3,                         /* */
  const char * config_file                = "analysis_configuration.xml"    /* */
)
{
  using namespace std;
  using namespace t2kana;

  // Global options                                                                                                                         
  Opt::Instance()->gOscMixingParamChoice = conv::kNuMixSinsq2Theta;
  Opt::Instance()->gMassHierarchyAssumed = conv::kNmhNormal;


  NuOsc3Flavour calc_osc3f;

  // set physics hypothesis
  PhysHypothesis ph(hypothesis);
  if(ph.Hypothesis() == conv::kHypNuMuDisapp2F) {
    ph.SetCurrValue(conv::kPhypMixing23,  input_mixing23);
    ph.SetCurrValue(conv::kPhypDm23sq,         input_dm23sq);
  } else
  if(ph.Hypothesis() == conv::kHypOsc3F) {
    ph.SetCurrValue(conv::kPhypMixing12,  input_mixing12);
    ph.SetCurrValue(conv::kPhypMixing13,  input_mixing13);
    ph.SetCurrValue(conv::kPhypMixing23,  input_mixing23);
    ph.SetCurrValue(conv::kPhypDeltaCP,        input_deltacp);
    ph.SetCurrValue(conv::kPhypDm12sq,         input_dm12sq);
    ph.SetCurrValue(conv::kPhypDm23sq,         input_dm23sq);
  }
  cout << "Physics hypothesis set: " << endl;
  cout << ph;

  // create & config superk spectrum predictor
  const char * config = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);
  cout << "Reading analysis configuration from: " << config << endl;


  // define the grid
  // deltacp, mixing13 range and step
  double mixing13_min   = 0.0;
  double mixing13_max   = 0.4;
  double mixing13_step  = 0.004;
  double deltacp_min         = -1 * conv::kpi;
  double deltacp_max         = conv::kpi;
  double deltacp_step        = 2*conv::kpi/90.;
  // make sure that min and max are at the bin centre and not at the low/high edge of the first/last bin
  int mixing13_nbins = 1 + (int) ((mixing13_max - mixing13_min) / mixing13_step);
  int deltacp_nbins       = 1 + (int) ((deltacp_max       - deltacp_min      ) / deltacp_step      );
  //energy                                                                                                                            
  double energy_max  = 10.000;
  double energy_min  = 0.250;
  int energy_nbins = 500;
  double energy_step = (energy_max - energy_min)/(energy_nbins-1);
  

  //delta_cp vs energy surfaces                                                                                                                  
  TH2D *NuEToNuE3f_dcp         =  new TH2D("NuEToNuE3f_dcp","3 Flavor P_{#nu_{e}#rightarrow#nu_{e}}",
					   energy_nbins,  energy_min -0.5*energy_step , energy_max + 0.5*energy_step,
					   deltacp_nbins, deltacp_min - 0.5*deltacp_step, deltacp_max + 0.5*deltacp_step);
  TH2D *NuEToNuMu3f_dcp        =  new TH2D("NuEToNuMu3f_dcp","3 Flavor P_{#nu_{e}#rightarrow#nu_{#mu}}",
					   energy_nbins,  energy_min -0.5*energy_step , energy_max + 0.5*energy_step,
					   deltacp_nbins, deltacp_min - 0.5*deltacp_step, deltacp_max + 0.5*deltacp_step);
  TH2D *NuMuToNuE3f_dcp        =  new TH2D("NuMuToNuE3f_dcp","3 Flavor P_{#nu_{#mu}#rightarrow#nu_{e}}",
					   energy_nbins,  energy_min -0.5*energy_step , energy_max + 0.5*energy_step,
					   deltacp_nbins, deltacp_min - 0.5*deltacp_step, deltacp_max + 0.5*deltacp_step);
  TH2D *NuMuToNuMu3f_dcp       =  new TH2D("NuMuToNuMu3f_dcp","3 Flavor P_{#nu_{#mu}#rightarrow#nu_{#mu}}",
					   energy_nbins,  energy_min -0.5*energy_step , energy_max + 0.5*energy_step,
					   deltacp_nbins, deltacp_min - 0.5*deltacp_step, deltacp_max + 0.5*deltacp_step);
  //mix13 vs energy surfaces                                                                                                   
  TH2D *NuEToNuE3f_mix13       =  new TH2D("NuEToNuE3f_mix13","3 Flavor P_{#nu_{e}#rightarrow#nu_{e}}",
					   energy_nbins,  energy_min -0.5*energy_step , energy_max + 0.5*energy_step,
					   mixing13_nbins, mixing13_min - 0.5*mixing13_step, mixing13_max + 0.5*mixing13_step);
  TH2D *NuEToNuMu3f_mix13      =  new TH2D("NuEToNuMu3f_mix13","3 Flavor P_{#nu_{e}#rightarrow#nu_{#mu}}",
					   energy_nbins,  energy_min -0.5*energy_step , energy_max + 0.5*energy_step,
					   mixing13_nbins, mixing13_min - 0.5*mixing13_step, mixing13_max + 0.5*mixing13_step);
  TH2D *NuMuToNuE3f_mix13      =  new TH2D("NuMuToNuE3f_mix13","3 Flavor P_{#nu_{#mu}#rightarrow#nu_{e}}",
					   energy_nbins,  energy_min -0.5*energy_step , energy_max + 0.5*energy_step,
					   mixing13_nbins, mixing13_min - 0.5*mixing13_step, mixing13_max + 0.5*mixing13_step);
  TH2D *NuMuToNuMu3f_mix13     =  new TH2D("NuMuToNuMu3f_mix13","3 Flavor P_{#nu_{#mu}#rightarrow#nu_{#mu}}",
					   energy_nbins,  energy_min -0.5*energy_step , energy_max + 0.5*energy_step,
					   mixing13_nbins, mixing13_min - 0.5*mixing13_step, mixing13_max + 0.5*mixing13_step);


  for(int ibin = 1; ibin <= deltacp_nbins; ibin++) {
    double dcp = deltacp_min + (ibin-1)*deltacp_step;
    for (int jbin = 1; jbin<= energy_nbins; jbin++){
      double energy = energy_min + (jbin-1)*energy_step;
      
      calc_osc3f.SetOscParameters(input_mixing12,input_mixing13,input_mixing23,input_dm12sq,input_dm23sq,dcp);
      
      NuEToNuE3f_dcp  ->Fill(energy,        dcp,       calc_osc3f.CalcMatterOscProb(conv::kPdgNuE,     conv::kPdgNuE,       energy));
      NuEToNuMu3f_dcp ->Fill(energy,        dcp,       calc_osc3f.CalcMatterOscProb(conv::kPdgNuE,     conv::kPdgNuMu,      energy));
      NuMuToNuE3f_dcp ->Fill(energy,        dcp,       calc_osc3f.CalcMatterOscProb(conv::kPdgNuMu,    conv::kPdgNuE,       energy));
      NuMuToNuMu3f_dcp->Fill(energy,        dcp,       calc_osc3f.CalcMatterOscProb(conv::kPdgNuMu,    conv::kPdgNuMu,      energy));
    }
  }
  
  
  for(int ibin = 1; ibin <= mixing13_nbins; ibin++) {
    double mix13 = mixing13_min + (ibin-1)*mixing13_step;
    for (int jbin = 1; jbin<= energy_nbins; jbin++){
      double energy = energy_min + (jbin-1)*energy_step;
   
      calc_osc3f.SetOscParameters(input_mixing12,mix13,input_mixing23,input_dm12sq,input_dm23sq,0.0);
      
      NuEToNuE3f_mix13  ->Fill(energy,        mix13,       calc_osc3f.CalcMatterOscProb(conv::kPdgNuE,     conv::kPdgNuE,       energy));
      NuEToNuMu3f_mix13 ->Fill(energy,        mix13,       calc_osc3f.CalcMatterOscProb(conv::kPdgNuE,     conv::kPdgNuMu,      energy));
      NuMuToNuE3f_mix13 ->Fill(energy,        mix13,       calc_osc3f.CalcMatterOscProb(conv::kPdgNuMu,    conv::kPdgNuE,       energy));
      NuMuToNuMu3f_mix13->Fill(energy,        mix13,       calc_osc3f.CalcMatterOscProb(conv::kPdgNuMu,    conv::kPdgNuMu,      energy));
    }
  }




  //now write the histogram
  TFile f("osc_prob_table.root", "recreate"); 
  f.cd();

  NuEToNuE3f_dcp    ->Write();
  NuEToNuMu3f_dcp   ->Write();
  NuMuToNuE3f_dcp   ->Write();
  NuMuToNuMu3f_dcp  ->Write();

  NuEToNuE3f_mix13    ->Write();
  NuEToNuMu3f_mix13   ->Write();
  NuMuToNuE3f_mix13   ->Write();
  NuMuToNuMu3f_mix13  ->Write();

  f.Write();
  f.Close();
  cout << "Done!" << endl;

}
