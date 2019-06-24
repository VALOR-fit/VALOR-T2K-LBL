#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>

#include "TXMLEngine.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "THStack.h"
#include "TLine.h"
#include "TString.h"
#include "TLegend.h"
#include "TMath.h"
#include "TSystem.h"

#include "SuperkExpectation.h"
#include "BANFF.h"
#include "NuXSec.h"
#include "SuperkDetectorFSI.h"
#include "SuperkMC.h"
#include "SuperkObservation.h"
#include "UtlFitParamListFactory.h"

using namespace std;
using namespace valor;
using namespace valor::conv;

void SetRatioNullTo1(TH1D *hratiomy){
  for(int ibin=0;ibin<hratiomy->GetNbinsX();ibin++){
    double cont = hratiomy->GetBinContent(ibin+1);
    if(cont==0)
      hratiomy->SetBinContent(ibin+1,1.);
  }
}

void SetRatio(TH1D* hnom,TH1D* htweak,TH1D* hratio){
  for(int ibin=0;ibin<hnom->GetNbinsX();ibin++){
    double cont_nom   = hnom->GetBinContent(ibin+1);
    double cont_tweak = htweak->GetBinContent(ibin+1);
    double cont_ratio = -999.;
    if(cont_nom>0.){
      cont_ratio = cont_tweak / cont_nom; 
    }
    else{
      cont_ratio = 1.;
    }
    hratio->SetBinContent(ibin+1,cont_ratio);
  }
}

void DrawForOsc3pp(TH1D * hp1, TH1D * hm1, TH1D * hp3, TH1D * hm3, TH1D * hspecp1, TH1D * hspecm1, conv::PhysParamId_t p, const char * savename, TCanvas & c)
{
  //rangers
  double axis_min = TMath::Min(hp3->GetMinimum(), hm3->GetMinimum());
  double axis_max = TMath::Max(hp3->GetMaximum(), hm3->GetMaximum());
  axis_max *= 1.05;
  axis_min *= 0.95;
  hp3->GetYaxis()->SetRangeUser(axis_min, axis_max);
  hp3->GetYaxis()->SetTitle("Ratio to nominal");
  //styles
  hspecp1->SetLineColor(kRed);
  hp3->SetLineColor(kRed);
  hp1->SetLineColor(kRed);
  hspecm1->SetLineColor(kBlue);
  hm1->SetLineColor(kBlue);
  hm3->SetLineColor(kBlue);
  hp3->SetLineStyle(kDashed);
  hm3->SetLineStyle(kDashed);
  //draws
  c.cd(1);
  hspecp1->Draw();
  hspecm1->Draw("SAME");  
  c.cd(2);
  hp3->Draw();
  hp1->Draw("SAME");
  hm1->Draw("SAME");
  hm3->Draw("SAME");
  //lines
  TLine line(0, 1, hp3->GetXaxis()->GetXmax(), 1);
  line.SetLineColor(kBlack);
  line.Draw();
  //legend
  TLegend l(0.1, 0.9, 1, 1);
  l.SetNColumns(5);
  l.AddEntry((TObject*)0, util::AsString(p).c_str(), "");
  l.AddEntry(hp3, "+3 #sigma", "l");
  l.AddEntry(hp1, "+1 #sigma", "l");
  l.AddEntry(hm1, "-1 #sigma", "l");
  l.AddEntry(hm3, "-3 #sigma", "l");
  l.Draw();
  c.SaveAs(savename);
}

void GetTunedHists(TH1D *& hErecoTun, TH1D *& hRatioTun, TH1D ** hErecoTunMode, TH1D ** hRatioTunMode,
		   TH1D *& hEtrueTun, TH1D *& hTrueRatioTun, TH1D ** hEtrueTunMode, TH1D ** hTrueRatioTunMode,
		   TH1D * hRatioSum, const char * sigma, double tweak, string sysstr, const int ip,
		   SuperkExpectation & skexp, conv::SuperkSample_t s,
		   TH1D * hErecoNom, TH1D ** hErecoNomMode,
		   TH1D * hEtrueNom, TH1D ** hEtrueNomMode)
{
  cout << "hErecoTun before new " << hErecoTun << endl;
  hErecoTun = new TH1D(*skexp.CurrPredictedErecoSpectrum(s));
  cout << "hErecoTun after new " << hErecoTun << endl;
  char name[200];
  sprintf(name,"hErecoTun_%ssigma_%s", sigma, sysstr.c_str());
  char title[200];
  sprintf(title,"E reco, %s, %s sigma = %.3lf",sysstr.c_str(), sigma, tweak);
  hErecoTun->SetName (name);
  hErecoTun->SetTitle(title);
  cout << "Tot # of events " << sigma << " sigma (integral)   = " << hErecoTun->Integral() << endl;
  cout << "Tot # of events " << sigma << " sigma (getentries) = " << hErecoTun->GetEntries() << endl;

  //char name[200];
  sprintf(name,"hRatioTun_%s",sysstr.c_str());
  //char title[200];
  sprintf(title,"Ratio, %s, +1 sigma = %.3lf",sysstr.c_str(),tweak);
  hRatioTun = (TH1D*)hErecoTun->Clone(name);
  SetRatio(hErecoNom,hErecoTun,hRatioTun);

  //true
  hEtrueTun = new TH1D(*skexp.GetSuperkMC()->Get3DHist(s,conv::kSKMCtFin).ProjectionZ());
  //char name[200];
  sprintf(name,"hEtrueTun_%ssigma_%s", sigma, sysstr.c_str());
  //char title[200];
  sprintf(title,"E true, %s, %s sigma = %.3lf",sysstr.c_str(), sigma, tweak);
  hEtrueTun->SetName (name);
  hEtrueTun->SetTitle(title);
  cout << "Tot # of events " << sigma << " sigma (integral)   = " << hEtrueTun->Integral() << endl;
  cout << "Tot # of events " << sigma << " sigma (getentries) = " << hEtrueTun->GetEntries() << endl;

  //char name[200];
  sprintf(name,"hTrueRatioTun_%s",sysstr.c_str());
  //char title[200];
  sprintf(title,"TrueRatio, %s, +1 sigma = %.3lf",sysstr.c_str(),tweak);
  hTrueRatioTun = (TH1D*)hEtrueTun->Clone(name);
  SetRatio(hEtrueNom,hEtrueTun,hTrueRatioTun);

  if(hRatioSum) {
    if(hErecoTun->Integral() > 0)
      hRatioSum->SetBinContent(ip + 1, (hErecoTun->Integral() / hErecoTun->Integral()) - 1);
    else
      hRatioSum->SetBinContent(ip + 1, 0);
  }

  for (int im = 0; im < conv::knm; im++) {
    conv::SuperkMode_t m = conv::kSKReactionModes[s][im];
    if(m == conv::kSKmUndefined) continue;
    if(util::ArrIdx::Instance()->SuperkModeIdx(s,m) < 0) continue;
    string modestr = util::AsStringCompact(m);
    const char * mname = modestr.c_str();
        
    hErecoTunMode[im] = skexp.CurrPredictedErecoSpectrum(s,m);
    hErecoTunMode[im]->SetFillStyle(0);
        
    //char name[200];
    sprintf(name,"hErecoTun_%ssigma_%s_%s",sigma, sysstr.c_str(), modestr.c_str());
    //char name[200];
    sprintf(title,"E reco, %s, %s, %s sigma = %.3lf",sysstr.c_str(),modestr.c_str(), sigma, tweak);
    hErecoTunMode[im]->SetName (name);
    hErecoTunMode[im]->SetTitle(title);
        
    //char name[200];
    sprintf(name,"hRatioTun_%ssigma_%s_%s", sigma, sysstr.c_str(), mname);
    //char title[200];
    sprintf(title,"Ratio, %s, %s, %s sigma = %.3lf",sysstr.c_str(),mname, sigma, tweak);
    hRatioTunMode[im] = (TH1D*)hErecoTunMode[im]->Clone(name);
    SetRatio(hErecoNomMode[im],hErecoTunMode[im],hRatioTunMode[im]);

    //Etrue
    hEtrueTunMode[im] = new TH1D(*skexp.GetSuperkMC()->Get3DHist(s,m,conv::kSKMCtFin).ProjectionZ());
    hEtrueTunMode[im]->SetFillStyle(0);
        
    //char name[200];
    sprintf(name,"hEtrueTun_%ssigma_%s_%s",sigma, sysstr.c_str(), modestr.c_str());
    //char name[200];
    sprintf(title,"E true, %s, %s, %s sigma = %.3lf",sysstr.c_str(),modestr.c_str(), sigma, tweak);
    hEtrueTunMode[im]->SetName (name);
    hEtrueTunMode[im]->SetTitle(title);
        
    //char name[200];
    sprintf(name,"hTrueRatioTun_%ssigma_%s_%s", sigma, sysstr.c_str(), mname);
    //char title[200];
    sprintf(title,"TrueRatio, %s, %s, %s sigma = %.3lf",sysstr.c_str(),mname, sigma, tweak);
    hTrueRatioTunMode[im] = (TH1D*)hEtrueTunMode[im]->Clone(name);
    SetRatio(hEtrueNomMode[im],hEtrueTunMode[im],hTrueRatioTunMode[im]);
  }//im
}


void test_systematic_error_propagation_all_2015_joint_2dtempl
(
 // int          mode               = 1,              /* options: see the macro   */
 bool         update_rootfile    = false,          /* "update" ROOT output file? If false, will "recreate" */
 const char * sk_sample          = "MultiRe",      /* options: "1Rmu", "1Re", "RHC1Rmu", "RHC1Re", "MultiRmu, MultiRe" */
 const char * hypothesis         = "3f_osc",       /* options: "2f_numu_disap", "3f_osc", "3f_oscbar" , "nuosc_3f+nuisance"*/

 double       input_mixing12       = 0.32,//-99,    
 double       input_mixing13       = 0.0256584, //-99,    
 double       input_mixing23       = 0.5,    
 double       input_deltacp        = 0,   
 double       input_dm12sq         = 7.6E-5,   
 double       input_dm23sq         = 2.4E-3,     

 // double       input_mixing12       = 0.,    
 // double       input_mixing13       = 0.,    
 // double       input_mixing23       = 0.,
 // double       input_deltacp        = 0.,
 // double       input_dm12sq         = 7.53E-5,   
 // double       input_dm23sq         = 2.509E-3,     
 
 bool         save_toy_expts     = false,                                       /* save the ereco spectra? */
 
 const char * config_file        = "analysis_configuration_hkfluxsophie_10year_postfit.xml",
 const char * fit_param_list     = "nuosc3f+nuisance",
 
 const char * mixparamchoice     = "sinsqtheta",
 const char * hierarchy          = "normal",
 
 int    tablesyst_idx            = -999,   /* -999 --> full table merged modes, 
                                              0-Nsyst --> mode table, -99 --> Document with all the mode table */
 string tablesyst_name           = "f_tamesu", //"f_banff17", /* follow util::AsStringCompact() convention */
 const bool   drawplots          = false,
 const bool   per_mode           = false,//true,
 const char*  tag                = "postbanff_osc",
 const bool   make_osc3pp_inputs = true
 )
{
  using namespace std;
  using namespace valor;
  using namespace valor::conv;
    
  double spectrum_line_width = 0.8;
  double ratio_line_width    = 1.0;
    
  int NtotParam=0;
  
  bool RHCFlag = false;
  TString str_sk_sample(sk_sample);
  if(str_sk_sample.Contains("RHC"))
    RHCFlag = true;
    
  // Global options
  Opt::Instance()->gMassHierarchyAssumed = util::HierarchyFromString(hierarchy);
  Opt::Instance()->gOscMixingParamChoice = util::MixingTypeFromString(mixparamchoice);
  if(Opt::Instance()->Undefined()) {
    cerr << "Global options not defined correctly. Exiting..." << endl;
    cerr << *Opt::Instance() << endl;
    exit(1);
  }
    
  // Get the sample
  const conv::SuperkSample_t s = util::SKSampleFromString(sk_sample); // requested sample
  if(s == conv::kSKsUndefined) {
    cerr << "** Undefined sample!" << endl;
    return;
  }
  const int is = util::ArrIdx::Instance()->SuperkSampleIdx(s); // index of specified sample
  const string sample_name = util::AsStringCompact(s);
  const char * sname = sample_name.c_str();

  cout << endl;
  cout << "Running sample: " << util::AsString(s) << endl;
  cout << endl;

  // Analysis configuration
  const char * config  = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);
  //T2KOpt::Instance()->ReadXmlPath(config);
  string sconfig(config);
  T2KOpt::Instance()->ReadXmlPath(sconfig);

  // Set physics hypothesis
  if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsq2Theta) {
    if(input_mixing13 < -1) input_mixing13 = conv::kGlobFitSinsq2Theta13;
    if(input_mixing23 < -1) input_mixing23 = conv::kGlobFitSinsq2Theta23;
    if(input_mixing12 < -1) input_mixing12 = conv::kGlobFitSinsq2Theta12;
  }
  else {
    if(input_mixing13 < -1) input_mixing13 = conv::kGlobFitSinsqTheta13;
    if(input_mixing23 < -1) input_mixing23 = conv::kGlobFitSinsqTheta23;
    if(input_mixing12 < -1) input_mixing12 = conv::kGlobFitSinsqTheta12;
  }
  if(input_dm23sq < -1)   input_dm23sq  = conv::kGlobFitDm23sq;
  if(input_dm12sq < -1)   input_dm12sq  = conv::kGlobFitDm12sq;
  if(input_deltacp < -90) input_deltacp = conv::kGlobFitDeltaCP;
  PhysHypothesis ph(hypothesis);
  if(ph.Hypothesis() == conv::kHypNuMuDisapp2F) {
    ph.SetCurrValue(conv::kPhypMixing23,  input_mixing23);
    ph.SetCurrValue(conv::kPhypDm23sq,    input_dm23sq);
  } else
    if(ph.Hypothesis() == conv::kHypOsc3F) {
      ph.SetCurrValue(conv::kPhypMixing12,  input_mixing12);
      ph.SetCurrValue(conv::kPhypMixing13,  input_mixing13);
      ph.SetCurrValue(conv::kPhypMixing23,  input_mixing23);
      ph.SetCurrValue(conv::kPhypCPV    ,   input_deltacp);
      ph.SetCurrValue(conv::kPhypDm12sq,    input_dm12sq);
      ph.SetCurrValue(conv::kPhypDm23sq,    input_dm23sq);
    } else
      if(ph.Hypothesis() == conv::kHypOsc3FBar){
        ph.SetCurrValue(conv::kPhypMixing12,  input_mixing12);
        ph.SetCurrValue(conv::kPhypMixing13,  input_mixing13);
        ph.SetCurrValue(conv::kPhypMixing23,  input_mixing23);
        ph.SetCurrValue(conv::kPhypCPV    ,   input_deltacp);
        ph.SetCurrValue(conv::kPhypDm12sq,    input_dm12sq);
        ph.SetCurrValue(conv::kPhypDm23sq,    input_dm23sq);
        ph.SetCurrValue(conv::kPhypMixing12Bar,  input_mixing12);
        ph.SetCurrValue(conv::kPhypMixing13Bar,  input_mixing13);
        ph.SetCurrValue(conv::kPhypMixing23Bar,  input_mixing23);
        ph.SetCurrValue(conv::kPhypCPVBar    ,   input_deltacp);
        ph.SetCurrValue(conv::kPhypDm12sqBar,    input_dm12sq);
        ph.SetCurrValue(conv::kPhypDm23sqBar,    input_dm23sq);
      }
    
  // Print options
  cout << "\n[Job inputs]" << endl;
  cout << "- Reading analysis configuration from: " << sconfig << endl;
  cout << "- Physics hypothesis: \n" << ph << endl;
  cout << *Opt::Instance() << endl;
    
  //output file name
  string oscparams =
    util::AsString(ph.Hypothesis(),input_mixing12,input_mixing13,input_mixing23,input_deltacp,input_dm12sq,input_dm23sq);
    
  // Precison for printed-out numerical results
  int nprecisionnev = 4; // <-- for number of events
  int nprecisionpct = 4; // <-- for % changes
    
  // Create & config superk spectrum predictor
  SuperkExpectation skexp;
  skexp.Read(sconfig);

  SuperkObservation skobs;
    
  // Define the parameter list

  FitParamList * fparm = util::NamedParameterList(fit_param_list);
  const int np = fparm->GetNumberParameters();
  const int np_loop = np;
  fparm->Print();
  
    
  // Define the tweak histograms
    
  TH1D *hErecoTun;
  TH1D *hErecoTunMode[conv::knm];
  THStack *hErecoTunCheck;

  TH1D *hRatioTun;
  TH1D *hRatioTunMode[conv::knm];
  TH1D *hErecoNom;
  TH1D *hErecoNomMode[conv::knm];
    
  TH1D *hErecoTun_p1sigma_all;
  TH1D *hErecoTun_m1sigma_all;
  TH1D *hErecoTun_p1sigma[np];
  TH1D *hErecoTun_m1sigma[np];
  TH1D *hErecoTunMode_p1sigma[np][conv::knm];
  TH1D *hErecoTunMode_m1sigma[np][conv::knm];
  TH1D *hRatioSum_p1sigma = new TH1D("hRatioSum_p1sigma", "Ratio of total +1 sigma event to nominal for each parameter", np, 0, np);
  TH1D *hRatioSum_m1sigma = new TH1D("hRatioSum_m1sigma", "Ratio of total -1 sigma event to nominal for each parameter", np, 0, np);
    
  TH1D *hErecoTun_p3sigma_all;
  TH1D *hErecoTun_m3sigma_all;
  TH1D *hErecoTun_p3sigma[np];
  TH1D *hErecoTun_m3sigma[np];
  TH1D *hErecoTunMode_p3sigma[np][conv::knm];
  TH1D *hErecoTunMode_m3sigma[np][conv::knm];
    
  TH1D *hRatioTun_p1sigma_all;
  TH1D *hRatioTun_m1sigma_all;
  TH1D *hRatioTun_p1sigma[np];
  TH1D *hRatioTun_m1sigma[np];
  TH1D *hRatioTunMode_p1sigma[np][conv::knm];
  TH1D *hRatioTunMode_m1sigma[np][conv::knm];
    
  TH1D *hRatioTun_p3sigma_all;
  TH1D *hRatioTun_m3sigma_all;
  TH1D *hRatioTun_p3sigma[np];
  TH1D *hRatioTun_m3sigma[np];
  TH1D *hRatioTunMode_p3sigma[np][conv::knm];
  TH1D *hRatioTunMode_m3sigma[np][conv::knm];

  //Etrue
  TH1D *hEtrueTun;
  TH1D *hEtrueTunMode[conv::knm];
  THStack *hEtrueTunCheck;

  TH1D *hTrueRatioTun;
  TH1D *hTrueRatioTunMode[conv::knm];
  TH1D *hEtrueNom;
  TH1D *hEtrueNomMode[conv::knm];
    
  TH1D *hEtrueTun_p1sigma[np];
  TH1D *hEtrueTun_m1sigma[np];
  TH1D *hEtrueTunMode_p1sigma[np][conv::knm];
  TH1D *hEtrueTunMode_m1sigma[np][conv::knm];
    
  TH1D *hEtrueTun_p3sigma[np];
  TH1D *hEtrueTun_m3sigma[np];
  TH1D *hEtrueTunMode_p3sigma[np][conv::knm];
  TH1D *hEtrueTunMode_m3sigma[np][conv::knm];
    
  TH1D *hTrueRatioTun_p1sigma[np];
  TH1D *hTrueRatioTun_m1sigma[np];
  TH1D *hTrueRatioTunMode_p1sigma[np][conv::knm];
  TH1D *hTrueRatioTunMode_m1sigma[np][conv::knm];
    
  TH1D *hTrueRatioTun_p3sigma[np];
  TH1D *hTrueRatioTun_m3sigma[np];
  TH1D *hTrueRatioTunMode_p3sigma[np][conv::knm];
  TH1D *hTrueRatioTunMode_m3sigma[np][conv::knm];
  //
    
  const int nskdetfsi     = conv::kNumOfFHCpRHCSKDetFSIParamsMultiRing; //conv::kNumOfFHCpRHCSKDetFSIParams;
  const int nbanff           = conv::kNumOfFHCpRHCParamsInBANFFfit;
    
  double vp1sigma_skdetfsi[nskdetfsi+1];
  double vp1sigma_banff[nbanff];
  double vp1sigma_skescale;
  double vp1sigma_skescale_rhc;
    
  double vm1sigma_skdetfsi[nskdetfsi+1];
  double vm1sigma_banff[nbanff];
  double vm1sigma_skescale;
  double vm1sigma_skescale_rhc;
    
  double vp3sigma_skdetfsi[nskdetfsi+1];
  double vp3sigma_banff[nbanff];
  double vp3sigma_skescale;
  double vp3sigma_skescale_rhc;
    
  double vm3sigma_skdetfsi[nskdetfsi+1];
  double vm3sigma_banff[nbanff];
  double vm3sigma_skescale;
  double vm3sigma_skescale_rhc;

  //Special values for the sigma tweaks for 2p2h shape parameter kNuisBANFFShape2p2h
  double  sp_2p2h_1sigma = 0.5;
  double  sp_2p2h_3sigma = 1.;    

    
  // central parameter vectors --> to restore the nominal MC for skdetfsi
  TVectorD *vcentral_skdetfsi = new TVectorD(nskdetfsi);
  for(int ip = 0; ip < nskdetfsi; ip++)
    (*vcentral_skdetfsi)[ip] = 1.0;    


  // Get tuned MC - with BANFF tuning
  //
  cout << "Getting nominal spectra" << endl;
  
  skexp.GetBANFF()->Tweak(0,0,"untuned");
  skexp.ResetSystematics();
  skexp.Generate(ph);
  hErecoTun = new TH1D(*skexp.CurrPredictedErecoSpectrum(s));
  //hErecoTun = new TH1D(*skexp.CurrPredictedErecoVsQ2Spectrum(s)->ProjectionZ()); // NEW 2D TEMPLATES
  hErecoTun->SetName("hErecoTun");
  hEtrueTun = new TH1D(*skexp.GetSuperkMC()->Get3DHist(s,conv::kSKMCtFin).ProjectionZ());
  hEtrueTun->SetName("hEtrueTun");  
  cout << endl;
  cout << "Tot # of events tuned reco (integral)   = " << hErecoTun->Integral() << endl;
  cout << "Tot # of events tuned reco (getentries) = " << hErecoTun->GetEntries() << endl;
  cout << "Tot # of events tuned true (integral)   = " << hEtrueTun->Integral() << endl;
  cout << "Tot # of events tuned true (getentries) = " << hEtrueTun->GetEntries() << endl;  
  cout << endl;

  cout << "Nbins in Ereco total sample " << hErecoTun->GetNbinsX() << " max at " << hErecoTun->GetXaxis()->GetXmax() << endl;
  cout << "Nbins in Etrue total sample " << hEtrueTun->GetNbinsX() << " max at " << hEtrueTun->GetXaxis()->GetXmax() << endl;

  hErecoTunCheck = new THStack("recosum", "sum of modes");
  hEtrueTunCheck = new THStack("truesum", "sum of modes");
  
  for (int im = 0; im < conv::knm; im++) {
    conv::SuperkMode_t m = conv::kSKReactionModes[s][im];
    if(m == conv::kSKmUndefined) continue;
    if(util::ArrIdx::Instance()->SuperkModeIdx(s,m) < 0) continue;
    string modestr = util::AsStringCompact(m);
    hErecoTunMode[im] = skexp.CurrPredictedErecoSpectrum(s,m);
    //hErecoTunMode[im] = skexp.CurrPredictedErecoVsQ2Spectrum(s,m)->ProjectionZ(); // NEW 2D TEMPLATES
    hErecoTunMode[im]->SetFillStyle(0);
    hEtrueTunMode[im] = new TH1D(*skexp.GetSuperkMC()->Get3DHist(s,m,conv::kSKMCtFin).ProjectionZ());
    hEtrueTunMode[im]->SetFillStyle(0);
    char name[200];
    sprintf(name,"hErecoTun_%s",modestr.c_str());
    hErecoTunMode[im]->SetName(name);
    hErecoTunCheck->Add(hErecoTunMode[im]);
    sprintf(name,"hEtrueTun_%s",modestr.c_str());
    hEtrueTunMode[im]->SetName(name);
    hEtrueTunCheck->Add(hEtrueTunMode[im]);
  }//im

  cout << "Nbins in Ereco total sample mode0 " << hErecoTunMode[0]->GetNbinsX() << " max at " << hErecoTunMode[0]->GetXaxis()->GetXmax() << endl;
  cout << "Nbins in Etrue total sample mode0 " << hEtrueTunMode[0]->GetNbinsX() << " max at " << hEtrueTunMode[0]->GetXaxis()->GetXmax() << endl;

  hErecoTunCheck->SetName("hErecoTunCheck");    
  hEtrueTunCheck->SetName("hEtrueTunCheck");    
    
  // Get tweaked ratio histograms - with BANFF tuning
  char name[200];
  sprintf(name,"hRatioTun");
  hRatioTun = (TH1D*)hErecoTun->Clone(name);
  SetRatio(hErecoTun,hErecoTun,hRatioTun);

  //char name[200];
  sprintf(name,"hTrueRatioTun");
  hTrueRatioTun = (TH1D*)hEtrueTun->Clone(name);
  SetRatio(hEtrueTun,hEtrueTun,hTrueRatioTun);

  // Get nominal MC - No BANFF tuning
  skexp.ResetSystematics();
  skexp.GetBANFF()->Tweak(0,0,"untuned");
  skexp.Generate(ph);
  
  hErecoNom = new TH1D(*skexp.CurrPredictedErecoSpectrum(s));
  //hErecoNom = new TH1D(*skexp.CurrPredictedErecoVsQ2Spectrum(s)->ProjectionZ()); // NEW 2D TEMPLATES
  hErecoNom->SetName("hErecoNom");
  hEtrueNom = new TH1D(*skexp.GetSuperkMC()->Get3DHist(s,conv::kSKMCtFin).ProjectionZ());
  hEtrueNom->SetName("hEtrueNom");
    
  for (int im = 0; im < conv::knm; im++) {
    conv::SuperkMode_t m = conv::kSKReactionModes[s][im];
    if(m == conv::kSKmUndefined) continue;
    if(util::ArrIdx::Instance()->SuperkModeIdx(s,m) < 0) continue;
    string modestr = util::AsStringCompact(m);
    const char * mname = modestr.c_str();
        
    hErecoNomMode[im] = skexp.CurrPredictedErecoSpectrum(s,m);
    //hErecoNomMode[im] = skexp.CurrPredictedErecoVsQ2Spectrum(s,m)->ProjectionZ(); // NEW 2D TEMPLATES
    hErecoNomMode[im]->SetFillStyle(0);
    //char name[200];
    sprintf(name,"hErecoNomMode_%s",modestr.c_str());
    hErecoNomMode[im]->SetName(name);
        
    // Get tweaked ratio histograms - with BANFF tuning
    //char name[200];
    sprintf(name,"hRatioTunMode_%s",mname);
    hRatioTunMode[im] = (TH1D*)hErecoTunMode[im]->Clone(name);
    SetRatio(hErecoTunMode[im],hErecoTunMode[im],hRatioTunMode[im]);

    //True
    hEtrueNomMode[im] = new TH1D(*skexp.GetSuperkMC()->Get3DHist(s,m,conv::kSKMCtFin).ProjectionZ());
    hEtrueNomMode[im]->SetFillStyle(0);
    //char name[200];
    sprintf(name,"hEtrueNomMode_%s",modestr.c_str());
    hEtrueNomMode[im]->SetName(name);
        
    // Get tweaked ratio histograms - with BANFF tuning
    //char name[200];
    sprintf(name,"hTrueRatioTunMode_%s",mname);
    hTrueRatioTunMode[im] = (TH1D*)hEtrueTunMode[im]->Clone(name);
    SetRatio(hEtrueTunMode[im],hEtrueTunMode[im],hTrueRatioTunMode[im]);
  }//im
    
  for (int im = 0; im < conv::knm; im++) {
    conv::SuperkMode_t m = conv::kSKReactionModes[s][im];
    if(m == conv::kSKmUndefined) continue;
    if(util::ArrIdx::Instance()->SuperkModeIdx(s,m) < 0) continue;
    string modestr = util::AsStringCompact(m);
  }//im

  // Loop over systematics
    
  int ip_banff = 0;
  int ip_skdetfsi = 0;

  cout << "////////////////////----o----////////////////////" << endl;
  cout << "Print list of systematic parameters: " << endl;
  
  for(int ip = 0; ip < np_loop; ip++) {
    conv::PhysParamId_t p = fparm->WhichParam(ip);
    if(!fparm->IsNuisance(p)) continue;
    string sysstr = util::AsStringCompact(p); //if before fparm->IsNuisance(p) gives problems !!!
    if(sysstr=="f_rpa")       continue;  // no tweaks for RPA, which is a discrete parameter

    cout << "Tweak systematic parameter: " << sysstr << endl;
  }//ip

  cout << "////////////////////----o----////////////////////" << endl;
    
  for(int ip = 0; ip < np_loop; ip++) {
    conv::PhysParamId_t p = fparm->WhichParam(ip);
    /*     cout << "Param: " << ip << " -> " << p << " = " << util::AsStringCompact(p) <<
	   " BANFF: " << util::IsBANFF(p) << " RHC BANFF: " << util::IsRHCBANFF(p) <<
	   " Flux:" << util::IsFlux(p) << " xsec: " << util::IsXSec(p) <<
	   " SK: " << util::IsSKDetFSI(p) << " E: " << util::IsEnergyScale(p) << endl;*/

    if(!fparm->IsNuisance(p)) continue;
    string sysstr = util::AsStringCompact(p); //if before fparm->IsNuisance(p) gives problems !!!
    if(sysstr=="f_rpa")       continue;  // no tweaks for RPA, which is a discrete parameter
        
    cout << endl;
    cout << "///////////////////////////////////////////////" << endl;
    cout << endl;
    cout << "Tweak systematic parameter: " << sysstr << endl;

    int idx=-999;
    double tweak_p1=0.;
    double tweak_m1=0.;
    double tweak_p3=0.;
    double tweak_m3=0.;
    double bestfit=0.;
        
    if(util::IsSKDetFSI(p)){
             
      idx = skexp.GetSuperkDetectorFSI()->Param2Idx(p);
      cout << "idx\t" << idx << endl;
      if(p != conv::kNuisSuperkEnergyScale){
        bestfit = 1.;

        tweak_p1 = bestfit + (1 * sqrt(skexp.GetSuperkDetectorFSI()->CovMatrix()(idx,idx)));
        tweak_m1 = bestfit - (1 * sqrt(skexp.GetSuperkDetectorFSI()->CovMatrix()(idx,idx)));
        //tweak_p1 = TMath::Max(tweak_p1,0.);
        //tweak_m1 = TMath::Max(tweak_m1,0.);
        
        vp1sigma_skdetfsi[ip_skdetfsi] = tweak_p1;
        vm1sigma_skdetfsi[ip_skdetfsi] = tweak_m1;
        cout << "idx\t" << idx << endl;
        tweak_p3 = bestfit + (3 * sqrt(skexp.GetSuperkDetectorFSI()->CovMatrix()(idx,idx)));
        tweak_m3 = bestfit - (3 * sqrt(skexp.GetSuperkDetectorFSI()->CovMatrix()(idx,idx)));
        vp3sigma_skdetfsi[ip_skdetfsi] = tweak_p3;
        vm3sigma_skdetfsi[ip_skdetfsi] = tweak_m3;
                
        ip_skdetfsi++;
      }//skdetfsi
      else if(p == conv::kNuisSuperkEnergyScale){
	bestfit = 0.;
	tweak_p1 = bestfit + (1 * sqrt(skexp.GetSuperkDetectorFSI()->CovMatrix()(idx,idx)));
	tweak_m1 = bestfit - (1 * sqrt(skexp.GetSuperkDetectorFSI()->CovMatrix()(idx,idx)));
	//   //tweak_p1 = TMath::Max(tweak_p1,0.);
	//   //tweak_m1 = TMath::Max(tweak_m1,0.);
        
	vp1sigma_skdetfsi[ip_skdetfsi] = tweak_p1;
	vm1sigma_skdetfsi[ip_skdetfsi] = tweak_m1;
	//   cout << "idx\t" << idx << endl;
	tweak_p3 = bestfit + (3 * sqrt(skexp.GetSuperkDetectorFSI()->CovMatrix()(idx,idx)));
	tweak_m3 = bestfit - (3 * sqrt(skexp.GetSuperkDetectorFSI()->CovMatrix()(idx,idx)));

	vp3sigma_skdetfsi[ip_skdetfsi] = tweak_p3;
	vm3sigma_skdetfsi[ip_skdetfsi] = tweak_m3;
                
	ip_skdetfsi++; 
      }//escale
            
    }//skdetfsi
    else if(util::IsBANFF(p)){
      // int idx_old = util::ArrIdx::Instance()->FHCpRHCBANFFSystIdx(p);
      idx = skexp.GetBANFF()->Param2Idx(p);
      // cout << "Index difference " << idx_old - idx << endl;
      if(idx>=0){
        bestfit = skexp.GetBANFF()->InitialValues()(idx);
        cout << "best fit\t" << bestfit << endl;
        
        if (p == kNuisBANFFShape2p2h    ||
            p == kRHCNuisBANFFShape2p2h   ){
          tweak_p1 = bestfit + sp_2p2h_1sigma;
          tweak_m1 = bestfit - sp_2p2h_1sigma;
          vp1sigma_banff[ip_banff] = tweak_p1;
          vm1sigma_banff[ip_banff] = tweak_m1;
                
          tweak_p3 = bestfit + sp_2p2h_3sigma;
          tweak_m3 = bestfit - sp_2p2h_3sigma;
          vp3sigma_banff[ip_banff] = tweak_p3;
          vm3sigma_banff[ip_banff] = tweak_m3;
                
          ip_banff++;
          
        } else {

	  tweak_p1 = bestfit + (1 * sqrt(skexp.GetBANFF()->CovMatrix()(idx,idx)));
	  tweak_m1 = bestfit - (1 * sqrt(skexp.GetBANFF()->CovMatrix()(idx,idx)));
	  vp1sigma_banff[ip_banff] = tweak_p1;
	  vm1sigma_banff[ip_banff] = tweak_m1;
                
	  tweak_p3 = bestfit + (3 * sqrt(skexp.GetBANFF()->CovMatrix()(idx,idx)));
	  tweak_m3 = bestfit - (3 * sqrt(skexp.GetBANFF()->CovMatrix()(idx,idx)));
	  vp3sigma_banff[ip_banff] = tweak_p3;
	  vm3sigma_banff[ip_banff] = tweak_m3;
                
	  ip_banff++;
        }
      }
    }//banff
    
    else if(util::IsXSec(p)){
      if(idx>=0){
        // bestfit = 0;
        // tweak_p1 = skexp.GetNuXSec()->Tweak(kRHCNuisXSecRPA , 1.);
        // tweak_m1 = skexp.GetNuXSec()->Tweak(kRHCNuisXSecRPA , 0.);
      }
    }//xsec
    else{
      cout << endl;     
      cout << "Systematic parameter: " << p << " " << util::AsString(p) << " not found!!!" << endl;
      cout << "IsSKDetFSI?\t " <<  util::IsSKDetFSI(p) << "\n"
           << "IsBANFF?\t "    <<  util::IsBANFF(p)    << "\n"
           << "IsXSec?\t "     <<  util::IsXSec(p)     << endl;
      exit(1);
    }//which class is p in?
    
    
    // Get the tweaked spectrum histograms
    
    // +1 sigma tweak systematic
    cout << "Getting +1 sigma" << endl;
    skexp.GetBANFF()->Tweak(1,0,"");
    skexp.ResetSystematics();
    cout << "Twk p1\t" << sysstr << "\t" << tweak_p1 << endl;
    if(util::IsSKDetFSI(p))   cout << "Nominal = " << skexp.GetSuperkDetectorFSI()->CurrentValue(p)       << endl;
    else if(util::IsBANFF(p)) cout << "Nominal = " << skexp.GetBANFF()            ->CurrentValue(p)       << endl;
    //else if(util::IsXSec(p))  cout << "Nominal = " << skexp.GetNuXSec()           ->FHCpRHCCurrentValue(p)<< endl;
    else{
      cout << "Parameter " << util::AsString(p) << " doesn't belong to any group" << endl;
      exit(1);
    }
    skexp.TweakSystematic(sysstr,tweak_p1);
    if(util::IsSKDetFSI(p))   cout << "+1sigma = " << skexp.GetSuperkDetectorFSI()->CurrentValue(p)       << endl;
    else if(util::IsBANFF(p)) cout << "+1sigma = " << skexp.GetBANFF()            ->CurrentValue(p)       << endl;
    //else if(util::IsXSec(p))  cout << "+1sigma = " << skexp.GetNuXSec()           ->FHCpRHCCurrentValue(p)<< endl;
    else{
      cout << "Parameter " << util::AsString(p) << " doesn't belong to any group" << endl;
      exit(1);
    }
    skexp.Generate(ph);

    cout << "hErecoTun_p1sigma[ip] before GetTunedHists() " << hErecoTun_p1sigma[ip] << endl;
    GetTunedHists(hErecoTun_p1sigma[ip], hRatioTun_p1sigma[ip], hErecoTunMode_p1sigma[ip],  hRatioTunMode_p1sigma[ip],
		  hEtrueTun_p1sigma[ip],  hTrueRatioTun_p1sigma[ip],  hEtrueTunMode_p1sigma[ip],  hTrueRatioTunMode_p1sigma[ip],
		  hRatioSum_p1sigma, "p1", tweak_p1, sysstr, ip,
		  skexp, s,
		  hErecoTun, hErecoTunMode,
		  hEtrueTun, hEtrueTunMode);
    cout << "hErecoTun_p1sigma[ip] after GetTunedHists() " << hErecoTun_p1sigma[ip] << endl;
    skexp.PrintSystematics();

    // -1 sigma tweak systematic
    cout << "Getting -1 sigma" << endl;
    skexp.GetBANFF()->Tweak(1,0,"");
    skexp.ResetSystematics();
    cout << "Twk m1\t" << sysstr << "\t" << tweak_m1 << endl;
    if(util::IsSKDetFSI(p))   cout << "Nominal = " << skexp.GetSuperkDetectorFSI()->CurrentValue(p)       << endl;
    else if(util::IsBANFF(p)) cout << "Nominal = " << skexp.GetBANFF()            ->CurrentValue(p)       << endl;
    //else if(util::IsXSec(p))  cout << "Nominal = " << skexp.GetNuXSec()           ->FHCpRHCCurrentValue(p)<< endl;
    else{
      cout << "Parameter " << util::AsString(p) << " doesn't belong to any group" << endl;
      exit(1);
    }
    skexp.TweakSystematic(sysstr,tweak_m1);
    if(util::IsSKDetFSI(p))   cout << "-1sigma = " << skexp.GetSuperkDetectorFSI()->CurrentValue(p)       << endl;
    else if(util::IsBANFF(p)) cout << "-1sigma = " << skexp.GetBANFF()            ->CurrentValue(p)       << endl;
    //else if(util::IsXSec(p))  cout << "-1sigma = " << skexp.GetNuXSec()           ->FHCpRHCCurrentValue(p)<< endl;
    else{
      cout << "Parameter " << util::AsString(p) << " doesn't belong to any group" << endl;
      exit(1);
    }
    skexp.Generate(ph);
    
    GetTunedHists(hErecoTun_m1sigma[ip], hRatioTun_m1sigma[ip], hErecoTunMode_m1sigma[ip],  hRatioTunMode_m1sigma[ip],
		  hEtrueTun_m1sigma[ip],  hTrueRatioTun_m1sigma[ip],  hEtrueTunMode_m1sigma[ip],  hTrueRatioTunMode_m1sigma[ip],
		  hRatioSum_m1sigma, "m1", tweak_m1, sysstr, ip,
		  skexp, s,
		  hErecoTun, hErecoTunMode,
		  hEtrueTun, hEtrueTunMode);
    
    // +3 sigma tweak systematic
    cout << "Getting +3 sigma" << endl;
    skexp.GetBANFF()->Tweak(1,0,"");
    skexp.ResetSystematics();
    cout << "Twk p3\t" << sysstr << "\t" << tweak_p3 << endl;
    if(util::IsSKDetFSI(p))   cout << "Nominal = " << skexp.GetSuperkDetectorFSI()->CurrentValue(p)       << endl;
    else if(util::IsBANFF(p)) cout << "Nominal = " << skexp.GetBANFF()            ->CurrentValue(p)       << endl;
    //else if(util::IsXSec(p))  cout << "Nominal = " << skexp.GetNuXSec()           ->FHCpRHCCurrentValue(p)<< endl;
    else{
      cout << "Parameter " << util::AsString(p) << " doesn't belong to any group" << endl;
      exit(1);
    }
    skexp.TweakSystematic(sysstr,tweak_p3);
    if(util::IsSKDetFSI(p))   cout << "+3sigma = " << skexp.GetSuperkDetectorFSI()->CurrentValue(p)       << endl;
    else if(util::IsBANFF(p)) cout << "+3sigma = " << skexp.GetBANFF()            ->CurrentValue(p)       << endl;
    //else if(util::IsXSec(p))  cout << "+3sigma = " << skexp.GetNuXSec()           ->FHCpRHCCurrentValue(p)<< endl;
    else{
      cout << "Parameter " << util::AsString(p) << " doesn't belong to any group" << endl;
      exit(1);
    }
    skexp.Generate(ph);
    
    GetTunedHists(hErecoTun_p3sigma[ip], hRatioTun_p3sigma[ip], hErecoTunMode_p3sigma[ip],  hRatioTunMode_p3sigma[ip],
		  hEtrueTun_p3sigma[ip],  hTrueRatioTun_p3sigma[ip],  hEtrueTunMode_p3sigma[ip],  hTrueRatioTunMode_p3sigma[ip],
		  NULL, "p3", tweak_p3, sysstr, ip,
		  skexp, s,
		  hErecoTun, hErecoTunMode,
		  hEtrueTun, hEtrueTunMode);
    
    // -3 sigma tweak systematic
    // skexp.GetBANFF()->RestoreFHCpRHCTunedMC();
    cout << "Getting -3 sigma" << endl;
    skexp.GetBANFF()->Tweak(1,0,"");
    skexp.ResetSystematics();
    cout << "Twk m3\t" << sysstr << "\t" << tweak_m3 << endl;
    if(util::IsSKDetFSI(p))   cout << "Nominal = " << skexp.GetSuperkDetectorFSI()->CurrentValue(p)       << endl;
    else if(util::IsBANFF(p)) cout << "Nominal = " << skexp.GetBANFF()            ->CurrentValue(p)       << endl;
    //else if(util::IsXSec(p))  cout << "Nominal = " << skexp.GetNuXSec()           ->FHCpRHCCurrentValue(p)<< endl;
    else{
      cout << "Parameter " << util::AsString(p) << " doesn't belong to any group" << endl;
      exit(1);
    }
    skexp.TweakSystematic(sysstr,tweak_m3);
    if(util::IsSKDetFSI(p))   cout << "-3sigma = " << skexp.GetSuperkDetectorFSI()->CurrentValue(p)       << endl;
    else if(util::IsBANFF(p)) cout << "-3sigma = " << skexp.GetBANFF()            ->CurrentValue(p)       << endl;
    //else if(util::IsXSec(p))  cout << "-3sigma = " << skexp.GetNuXSec()           ->FHCpRHCCurrentValue(p)<< endl;
    else{
      cout << "Parameter " << util::AsString(p) << " doesn't belong to any group" << endl;
      exit(1);
    }
    skexp.Generate(ph);
    
    GetTunedHists(hErecoTun_m3sigma[ip], hRatioTun_m3sigma[ip], hErecoTunMode_m3sigma[ip],  hRatioTunMode_m3sigma[ip],
		  hEtrueTun_m3sigma[ip],  hTrueRatioTun_m3sigma[ip],  hEtrueTunMode_m3sigma[ip],  hTrueRatioTunMode_m3sigma[ip],
		  NULL, "m3", tweak_m3, sysstr, ip,
		  skexp, s,
		  hErecoTun, hErecoTunMode,
		  hEtrueTun, hEtrueTunMode);
    
    NtotParam++;
    
  } // end loop over systematics ip

  cout << "\n\n#########################################\n" << endl;
  cout << "Ended loop over the systematic parameters" << endl;
  cout << endl;
  cout << endl;
  
  // Get total +1 sigma tweak (all parameters)
    
  skexp.GetBANFF()->Tweak(1,0,"");
  skexp.ResetSystematics();
    
  TVectorD *twkv_banff_p1sigma;
  TVectorD *twkv_skdetfsi_p1sigma;    
    
  twkv_banff_p1sigma    = new TVectorD(nbanff,vp1sigma_banff);
  twkv_skdetfsi_p1sigma = new TVectorD(nskdetfsi+1,vp1sigma_skdetfsi);
  skexp.GetBANFF()->Tweak(1,twkv_banff_p1sigma,"");
  skexp.GetSuperkDetectorFSI()->Tweak(1,twkv_skdetfsi_p1sigma,"");
  skexp.Generate(ph);    
    
  hErecoTun_p1sigma_all = new TH1D(*skexp.CurrPredictedErecoSpectrum(s));
  //char name[200];
  sprintf(name,"hErecoTun_p1sigma_all");
  hErecoTun_p1sigma_all->SetName(name);

  // Get total +3 sigma tweak (all parameters)
    
  skexp.GetBANFF()->Tweak(1,0,"");
  skexp.ResetSystematics();
    
  TVectorD *twkv_banff_p3sigma;
  TVectorD *twkv_skdetfsi_p3sigma;    
    
  twkv_banff_p3sigma    = new TVectorD(nbanff,vp3sigma_banff);
  twkv_skdetfsi_p3sigma = new TVectorD(nskdetfsi+1,vp3sigma_skdetfsi);
  skexp.GetBANFF()->Tweak(1,twkv_banff_p3sigma,"");
  skexp.GetSuperkDetectorFSI()->Tweak(1,twkv_skdetfsi_p3sigma,"");
  skexp.Generate(ph);    
    
  hErecoTun_p3sigma_all = new TH1D(*skexp.CurrPredictedErecoSpectrum(s));
  //char name[200];
  sprintf(name,"hErecoTun_p3sigma_all");
  hErecoTun_p3sigma_all->SetName(name);
  
    
  // Get total -1 sigma tweak (all parameters)
    
  //  skexp.GetBANFF()->RestoreFHCpRHCTunedMC();
  skexp.GetBANFF()->Tweak(1,0,"");
  skexp.ResetSystematics();
    
  TVectorD *twkv_banff_m1sigma;
  TVectorD *twkv_skdetfsi_m1sigma;
    
    
  twkv_banff_m1sigma    = new TVectorD(nbanff,vm1sigma_banff);
  twkv_skdetfsi_m1sigma = new TVectorD(nskdetfsi+1,vm1sigma_skdetfsi);
  skexp.GetBANFF()->Tweak(1,twkv_banff_m1sigma,"");
  skexp.GetSuperkDetectorFSI()->Tweak(1,twkv_skdetfsi_m1sigma,"");
  skexp.Generate(ph);
    
  hErecoTun_m1sigma_all = new TH1D(*skexp.CurrPredictedErecoSpectrum(s));
  //char name[200];
  sprintf(name,"hErecoTun_m1sigma_all");
  hErecoTun_m1sigma_all->SetName(name);
    
  // Get total -3 sigma tweak (all parameters)
    
  skexp.GetBANFF()->Tweak(1,0,"");
  skexp.ResetSystematics();
    
  TVectorD *twkv_banff_m3sigma;
  TVectorD *twkv_skdetfsi_m3sigma;    
    
  twkv_banff_m3sigma    = new TVectorD(nbanff,vm3sigma_banff);
  twkv_skdetfsi_m3sigma = new TVectorD(nskdetfsi+1,vm3sigma_skdetfsi);
  skexp.GetBANFF()->Tweak(1,twkv_banff_m3sigma,"");
  skexp.GetSuperkDetectorFSI()->Tweak(1,twkv_skdetfsi_m3sigma,"");
  skexp.Generate(ph);    
    
  hErecoTun_m3sigma_all = new TH1D(*skexp.CurrPredictedErecoSpectrum(s));
  //char name[200];
  sprintf(name,"hErecoTun_m3sigma_all");
  hErecoTun_m3sigma_all->SetName(name);
  
  if(drawplots){
      
    // Draw +-1 sigma and +-3 sigma plots
    
    for(int ip = 0; ip < np_loop; ip++) {
      conv::PhysParamId_t p = fparm->WhichParam(ip);
      if(!fparm->IsNuisance(p)) continue;
      string sysstr = util::AsStringCompact(p); //if before fparm->IsNuisance(p) gives problems !!!
      if(sysstr=="f_rpa")       continue;  // no tweaks for RPA, which is a discrete parameter
        
      cout << sysstr << endl;
        
      // Spectra
      sprintf(name,"cspectra_%s",sysstr.c_str());
      TCanvas *cspectra = new TCanvas(name,name);
        
      double max_limit = 1.1 * TMath::Max(hErecoTun_p3sigma[ip]->GetMaximum(),hErecoTun_m3sigma[ip]->GetMaximum());
      hErecoTun_m1sigma[ip]->SetMaximum(max_limit);
      hErecoTun_m1sigma[ip]->SetMinimum(0.);
        
      hErecoTun_m1sigma[ip]->Draw();
      hErecoTun_m3sigma[ip]->Draw("same");
      hErecoTun_p1sigma[ip]->Draw("same");
      hErecoTun_p3sigma[ip]->Draw("same");
      hErecoTun->Draw("same");
        
      hErecoTun            ->SetLineWidth(spectrum_line_width);
      hErecoTun_m1sigma[ip]->SetLineWidth(spectrum_line_width);
      hErecoTun_m3sigma[ip]->SetLineWidth(spectrum_line_width);
      hErecoTun_p1sigma[ip]->SetLineWidth(spectrum_line_width);
      hErecoTun_p3sigma[ip]->SetLineWidth(spectrum_line_width);
        
      hErecoTun_m1sigma[ip]->SetLineColor(kViolet);
      hErecoTun_p1sigma[ip]->SetLineColor(kCyan);
      hErecoTun_m3sigma[ip]->SetLineColor(kRed);
      hErecoTun_p3sigma[ip]->SetLineColor(kBlue);
        
      cspectra->SaveAs(Form("Spectrum_%s_%s_nsk_%s.eps",sname,sysstr.c_str(),oscparams.c_str()));
        
      // Ratios
        
      //char name[200];
      sprintf(name,"cratio_%s",sysstr.c_str());
      TCanvas *cratio = new TCanvas(name,name);
        
      max_limit = 1.1 * fabs(TMath::Max(hRatioTun_p3sigma[ip]->GetMaximum(),hRatioTun_m3sigma[ip]->GetMaximum())-1);
      double min_limit = 1.1 * fabs(TMath::Min(hRatioTun_p3sigma[ip]->GetMinimum(),hRatioTun_m3sigma[ip]->GetMinimum())-1);
        
      max_limit = TMath::Max(max_limit,min_limit);
        
      hRatioTun_m1sigma[ip]->SetMaximum(1+max_limit);
      hRatioTun_m1sigma[ip]->SetMinimum(1-max_limit);
        
      hRatioTun_m1sigma[ip]->Draw();
      hRatioTun_m3sigma[ip]->Draw("same");
      hRatioTun_p1sigma[ip]->Draw("same");
      hRatioTun_p3sigma[ip]->Draw("same");
      hRatioTun->Draw("same");
        
      hRatioTun            ->SetLineWidth(ratio_line_width);
      hRatioTun_m1sigma[ip]->SetLineWidth(ratio_line_width);
      hRatioTun_m3sigma[ip]->SetLineWidth(ratio_line_width);
      hRatioTun_p1sigma[ip]->SetLineWidth(ratio_line_width);
      hRatioTun_p3sigma[ip]->SetLineWidth(ratio_line_width);
        
      hRatioTun_m1sigma[ip]->SetLineColor(kViolet);
      hRatioTun_p1sigma[ip]->SetLineColor(kCyan);
      hRatioTun_m3sigma[ip]->SetLineColor(kRed);
      hRatioTun_p3sigma[ip]->SetLineColor(kBlue);
      cratio->SaveAs(Form("Ratio_%s_%s_nsk_%s.eps",sname,sysstr.c_str(),oscparams.c_str()));
        
    } // end loop over systematics ip
    
    if (per_mode) {
      for (int im = 0; im < conv::knm; im++) {
        conv::SuperkMode_t m = conv::kSKReactionModes[s][im];
        if(m == conv::kSKmUndefined) continue;
        if(util::ArrIdx::Instance()->SuperkModeIdx(s,m) < 0) continue;
        string modestr = util::AsStringCompact(m);
        const char * mname = modestr.c_str();
            
        for(int ip = 0; ip < np_loop; ip++) {
          conv::PhysParamId_t p = fparm->WhichParam(ip);
          if(!fparm->IsNuisance(p)) continue;
          string sysstr = util::AsStringCompact(p); //if before fparm->IsNuisance(p) gives problems !!!
          if(sysstr=="f_rpa")       continue;  // no tweaks for RPA, which is a discrete parameter
                
          // Spectra
          sprintf(name,"cspectra_%s_%s",mname,sysstr.c_str());
          TCanvas *cspectra = new TCanvas(name,name);
                
          double max_limit = 1.1 * TMath::Max(hErecoTunMode_p3sigma[ip][im]->GetMaximum(),hErecoTunMode_m3sigma[ip][im]->GetMaximum());
          hErecoTunMode_m1sigma[ip][im]->SetMaximum(max_limit);
          hErecoTunMode_m1sigma[ip][im]->SetMinimum(0.);
                
          hErecoTunMode_m1sigma[ip][im]->Draw();
          hErecoTunMode_m3sigma[ip][im]->Draw("same");
          hErecoTunMode_p1sigma[ip][im]->Draw("same");
          hErecoTunMode_p3sigma[ip][im]->Draw("same");
          hErecoTunMode[im]            ->Draw("same");
                
          hErecoTunMode[im]            ->SetLineWidth(spectrum_line_width);
          hErecoTunMode_m1sigma[ip][im]->SetLineWidth(spectrum_line_width);
          hErecoTunMode_m3sigma[ip][im]->SetLineWidth(spectrum_line_width);
          hErecoTunMode_p1sigma[ip][im]->SetLineWidth(spectrum_line_width);
          hErecoTunMode_p3sigma[ip][im]->SetLineWidth(spectrum_line_width);
                
          hErecoTunMode_m1sigma[ip][im]->SetLineColor(kViolet);
          hErecoTunMode_p1sigma[ip][im]->SetLineColor(kCyan);
          hErecoTunMode_m3sigma[ip][im]->SetLineColor(kRed);
          hErecoTunMode_p3sigma[ip][im]->SetLineColor(kBlue);
                
          cspectra->SaveAs(Form("Spectrum_%s_%s_%s_nsk_%s.eps",sname,mname,sysstr.c_str(),oscparams.c_str()));
                
          // Ratios
                
          //char name[200];
          sprintf(name,"cratio_%s",sysstr.c_str());
          TCanvas *cratio = new TCanvas(name,name);
                
          max_limit = 1.1 * fabs(TMath::Max(hRatioTunMode_p3sigma[ip][im]->GetMaximum(),hRatioTunMode_m3sigma[ip][im]->GetMaximum())-1);
          double min_limit = 1.1 * fabs(TMath::Min(hRatioTunMode_p3sigma[ip][im]->GetMinimum(),hRatioTunMode_m3sigma[ip][im]->GetMinimum())-1);
                
          max_limit = TMath::Max(max_limit,min_limit);
                
          hRatioTunMode_m1sigma[ip][im]->SetMaximum(1+max_limit);
          hRatioTunMode_m1sigma[ip][im]->SetMinimum(1-max_limit);
                
          hRatioTunMode_m1sigma[ip][im]->Draw();
          hRatioTunMode_m3sigma[ip][im]->Draw("same");
          hRatioTunMode_p1sigma[ip][im]->Draw("same");
          hRatioTunMode_p3sigma[ip][im]->Draw("same");
          hRatioTunMode[im]            ->Draw("same");
                
          hRatioTunMode[im]            ->SetLineWidth(ratio_line_width);
          hRatioTunMode_m1sigma[ip][im]->SetLineWidth(ratio_line_width);
          hRatioTunMode_m3sigma[ip][im]->SetLineWidth(ratio_line_width);
          hRatioTunMode_p1sigma[ip][im]->SetLineWidth(ratio_line_width);
          hRatioTunMode_p3sigma[ip][im]->SetLineWidth(ratio_line_width);
                
          hRatioTunMode_m1sigma[ip][im]->SetLineColor(kViolet);
          hRatioTunMode_p1sigma[ip][im]->SetLineColor(kCyan);
          hRatioTunMode_m3sigma[ip][im]->SetLineColor(kRed);
          hRatioTunMode_p3sigma[ip][im]->SetLineColor(kBlue);
          cratio->SaveAs(Form("Ratio_%s_%s_%s_nsk_%s.eps",sname,mname,sysstr.c_str(),oscparams.c_str()));
                
        } // ip
      } //im
    } // per_mode
    
  } // end draw plots


  // Write ROOT file output

  const string file_mode = (update_rootfile ? "UPDATE" : "RECREATE");
  TFile * output_rootfile_osc3pp, * output_rootfile_osc3pp_true, *output_rootfile_osc3pp_mode;
  TCanvas c3pp, c3pp_split;
  if(make_osc3pp_inputs) {
    c3pp.SetRightMargin(0.05);
    c3pp.SaveAs(TString::Format("syst_variations_%s.pdf[", sname));
    c3pp.SaveAs(TString::Format("syst_variations_true_%s.pdf[", sname));
    c3pp_split.Divide(1,2);
    output_rootfile_osc3pp = new TFile(Form("%s.fij.t2k.root", sname), file_mode.c_str());
    output_rootfile_osc3pp_true = new TFile(Form("%s.fij.t2k.true.root", sname), file_mode.c_str());
    output_rootfile_osc3pp_mode = new TFile(Form("%s.fij.t2k.mode.root", sname), file_mode.c_str());
  }//make_osc3pp_inputs
  TFile * output_rootfile = new TFile(Form("systematic_variations%s_%s.root",tag,sname),file_mode.c_str());
  output_rootfile->cd();
  hErecoTun->Write(Form("Spectrum_%s_nominal",sname));
  hErecoTunCheck->Write(Form("Spectrum_%s_nominal_check",sname));
  int systcounter = 0;
  for(int ip = 0; ip < np_loop; ip++) {
    conv::PhysParamId_t p = fparm->WhichParam(ip);
    if(!fparm->IsNuisance(p)) continue;
    string sysstr = util::AsStringCompact(p); //if before fparm->IsNuisance(p) gives problems !!!
    if(sysstr=="f_rpa")       continue;  // no tweaks for RPA, which is a discrete parameter
    
    output_rootfile->cd();

    cout << sysstr << "\t" << hErecoTun_m1sigma[ip] << endl;
    cout << Form("Spectrum_%s_%s_m1",sname,sysstr.c_str()) << endl;
    cout << hErecoTun_m1sigma[ip]->Integral() << endl;
    hErecoTun_m1sigma[ip]->Write(Form("Spectrum_%s_%s_m1",sname,sysstr.c_str()));
    hErecoTun_m3sigma[ip]->Write(Form("Spectrum_%s_%s_m3",sname,sysstr.c_str()));
    hErecoTun_p1sigma[ip]->Write(Form("Spectrum_%s_%s_p1",sname,sysstr.c_str()));
    hErecoTun_p3sigma[ip]->Write(Form("Spectrum_%s_%s_p3",sname,sysstr.c_str()));
    
    hRatioTun_m1sigma[ip]->Write(Form("Ratio_%s_%s_m1",sname,sysstr.c_str()));
    hRatioTun_m3sigma[ip]->Write(Form("Ratio_%s_%s_m3",sname,sysstr.c_str()));
    hRatioTun_p1sigma[ip]->Write(Form("Ratio_%s_%s_p1",sname,sysstr.c_str()));
    hRatioTun_p3sigma[ip]->Write(Form("Ratio_%s_%s_p3",sname,sysstr.c_str()));
    
    if(make_osc3pp_inputs) {
      output_rootfile_osc3pp->cd();
      hRatioTun_p1sigma[ip]->Write(Form("syserre_%d",systcounter));
      if(util::IsUsingSplines(p)) {
	hRatioTun_p1sigma[ip]->Write(Form("syserre_%d_p1",systcounter));
	hRatioTun_p3sigma[ip]->Write(Form("syserre_%d_p3",systcounter));
	hRatioTun_m1sigma[ip]->Write(Form("syserre_%d_m1",systcounter));
	hRatioTun_m3sigma[ip]->Write(Form("syserre_%d_m3",systcounter));
      }
      output_rootfile_osc3pp_true->cd();
      hRatioTun_p1sigma[ip]->Write(Form("syserre_%d_%s",systcounter,sysstr.c_str()));
      hTrueRatioTun_p1sigma[ip]->Write(Form("syserretrue_%d_%s",systcounter,sysstr.c_str()));
      if(util::IsUsingSplines(p)) {
	hTrueRatioTun_p1sigma[ip]->Write(Form("syserretrue_%d_p1",systcounter));
	hTrueRatioTun_p3sigma[ip]->Write(Form("syserretrue_%d_p3",systcounter));
	hTrueRatioTun_m1sigma[ip]->Write(Form("syserretrue_%d_m1",systcounter));
	hTrueRatioTun_m3sigma[ip]->Write(Form("syserretrue_%d_m3",systcounter));
      }
      DrawForOsc3pp(hRatioTun_p1sigma[ip], hRatioTun_m1sigma[ip],
		    hRatioTun_p3sigma[ip], hRatioTun_m3sigma[ip], 
		    hErecoTun_p1sigma[ip], hErecoTun_m1sigma[ip], p,
		    TString::Format("syst_variations_%s.pdf", sname), c3pp_split);      
      DrawForOsc3pp(hTrueRatioTun_p1sigma[ip], hTrueRatioTun_m1sigma[ip],
		    hTrueRatioTun_p3sigma[ip], hTrueRatioTun_m3sigma[ip],
		    hEtrueTun_p1sigma[ip], hEtrueTun_m1sigma[ip], p,
		    TString::Format("syst_variations_true_%s.pdf", sname), c3pp_split);
      if(per_mode) {
	output_rootfile_osc3pp_mode->cd();
	for (int im = 0; im < conv::knm; im++) {
	  conv::SuperkMode_t m = conv::kSKReactionModes[s][im];
	  if(m == conv::kSKmUndefined) continue;
	  if(util::ArrIdx::Instance()->SuperkModeIdx(s,m) < 0) continue;
	  string modestr = util::AsStringCompact(m);
	  const char * mname = modestr.c_str();
	  hRatioTunMode_m1sigma[ip][im]->Write(Form("syserre_%d_%s_%s_m1",systcounter,mname,sysstr.c_str()));
	  hRatioTunMode_m3sigma[ip][im]->Write(Form("syserre_%d_%s_%s_m3",systcounter,mname,sysstr.c_str()));
	  hRatioTunMode_p1sigma[ip][im]->Write(Form("syserre_%d_%s_%s_p1",systcounter,mname,sysstr.c_str()));
	  hRatioTunMode_p3sigma[ip][im]->Write(Form("syserre_%d_%s_%s_p3",systcounter,mname,sysstr.c_str()));
	  hTrueRatioTunMode_m1sigma[ip][im]->Write(Form("syserretrue_%d_%s_%s_m1",systcounter,mname,sysstr.c_str()));
	  hTrueRatioTunMode_m3sigma[ip][im]->Write(Form("syserretrue_%d_%s_%s_m3",systcounter,mname,sysstr.c_str()));
	  hTrueRatioTunMode_p1sigma[ip][im]->Write(Form("syserretrue_%d_%s_%s_p1",systcounter,mname,sysstr.c_str()));
	  hTrueRatioTunMode_p3sigma[ip][im]->Write(Form("syserretrue_%d_%s_%s_p3",systcounter,mname,sysstr.c_str()));
	  if(false) {
	    DrawForOsc3pp(hRatioTunMode_p1sigma[ip][im], hRatioTunMode_m1sigma[ip][im],
			  hRatioTunMode_p3sigma[ip][im], hRatioTunMode_m3sigma[ip][im],
			  hErecoTunMode_p1sigma[ip][im], hErecoTunMode_p1sigma[ip][im], p,
			  TString::Format("syst_variations_%s.pdf", sname), c3pp_split);
	    DrawForOsc3pp(hTrueRatioTunMode_p1sigma[ip][im], hTrueRatioTunMode_m1sigma[ip][im],
			  hTrueRatioTunMode_p3sigma[ip][im], hTrueRatioTunMode_m3sigma[ip][im],
			  hEtrueTunMode_p1sigma[ip][im], hEtrueTunMode_p1sigma[ip][im], p,
			  TString::Format("syst_variations_true_%s.pdf", sname), c3pp_split);
	  }
	  hErecoTunMode_p1sigma[ip][im]->Write(Form("recospectrum_%d_%s_%s_p1",systcounter,mname,sysstr.c_str()));
	  hEtrueTunMode_p1sigma[ip][im]->Write(Form("truespectrum_%d_%s_%s_p1",systcounter,mname,sysstr.c_str()));
	}//im
      }//permode
      systcounter++;
    }//make_osc3pp_inputs
  }//ip

  if (per_mode && false) {
    output_rootfile->cd();
    for (int im = 0; im < conv::knm; im++) {
      conv::SuperkMode_t m = conv::kSKReactionModes[s][im];
      if(m == conv::kSKmUndefined) continue;
      if(util::ArrIdx::Instance()->SuperkModeIdx(s,m) < 0) continue;
      string modestr = util::AsStringCompact(m);
      const char * mname = modestr.c_str();
        
      hErecoTunMode[im]->Write(Form("Spectrum_%s_%s_nominal",sname,mname));
        
      for(int ip = 0; ip < np_loop; ip++) {
	conv::PhysParamId_t p = fparm->WhichParam(ip);
	if(!fparm->IsNuisance(p)) continue;
	string sysstr = util::AsStringCompact(p); //if before fparm->IsNuisance(p) gives problems !!!
	if(sysstr=="f_rpa")       continue;  // no tweaks for RPA, which is a discrete parameter
            
	hErecoTunMode_m1sigma[ip][im]->Write(Form("Spectrum_%s_%s_%s_m1",sname,mname,sysstr.c_str()));
	hErecoTunMode_m3sigma[ip][im]->Write(Form("Spectrum_%s_%s_%s_m3",sname,mname,sysstr.c_str()));
	hErecoTunMode_p1sigma[ip][im]->Write(Form("Spectrum_%s_%s_%s_p1",sname,mname,sysstr.c_str()));
	hErecoTunMode_p3sigma[ip][im]->Write(Form("Spectrum_%s_%s_%s_p3",sname,mname,sysstr.c_str()));
            
	hRatioTunMode_m1sigma[ip][im]->Write(Form("Ratio_%s_%s_%s_m1",sname,mname,sysstr.c_str()));
	hRatioTunMode_m3sigma[ip][im]->Write(Form("Ratio_%s_%s_%s_m3",sname,mname,sysstr.c_str()));
	hRatioTunMode_p1sigma[ip][im]->Write(Form("Ratio_%s_%s_%s_p1",sname,mname,sysstr.c_str()));
	hRatioTunMode_p3sigma[ip][im]->Write(Form("Ratio_%s_%s_%s_p3",sname,mname,sysstr.c_str()));
      } //ip
    } //im
  } // per_mode

  output_rootfile->Close();
  if(make_osc3pp_inputs) {
    output_rootfile_osc3pp->Close();
    output_rootfile_osc3pp_true->Close();
    output_rootfile_osc3pp_mode->Close();
    c3pp.SaveAs(TString::Format("syst_variations_%s.pdf]", sname));
    c3pp.SaveAs(TString::Format("syst_variations_true_%s.pdf]", sname));
    c3pp.Close();
  }//make_osc3pp_inputs


  // Save histo with fractional and total variation of the number of events

  cout << endl;
  cout << "NtotParam = " << NtotParam << endl;
  cout << endl;

  TH1D *hVarFrac_p1s = new TH1D("hVarFrac_p1s","+1sigma fractional variation",NtotParam,0,NtotParam); // store the +1sigma errors
  TH1D *hVarFrac_m1s = new TH1D("hVarFrac_m1s","-1sigma fractional variation",NtotParam,0,NtotParam); // store the -1sigma errors
  TH1D *hVarFrac_p3s = new TH1D("hVarFrac_p3s","+3sigma fractional variation",NtotParam,0,NtotParam); // store the +3sigma errors
  TH1D *hVarFrac_m3s = new TH1D("hVarFrac_m3s","-3sigma fractional variation",NtotParam,0,NtotParam); // store the -3sigma errors

  TH1D *hVarTot_p1s = new TH1D("hVarTot_p1s","+1sigma total variation",NtotParam,0,NtotParam); // store the +1sigma errors
  TH1D *hVarTot_m1s = new TH1D("hVarTot_m1s","-1sigma total variation",NtotParam,0,NtotParam); // store the -1sigma errors
  TH1D *hVarTot_p3s = new TH1D("hVarTot_p3s","+3sigma total variation",NtotParam,0,NtotParam); // store the +3sigma errors
  TH1D *hVarTot_m3s = new TH1D("hVarTot_m3s","-3sigma total variation",NtotParam,0,NtotParam); // store the -3sigma errors

  int iparhisto=0;

  for(int ip = 0; ip < np_loop; ip++) {
    conv::PhysParamId_t p = fparm->WhichParam(ip);
    if(!fparm->IsNuisance(p)) continue;
    string sysstr = util::AsStringCompact(p); // if before fparm->IsNuisance(p) gives problems !!
    if(sysstr=="f_rpa")       continue;  // no tweaks for RPA, which is a discrete parameter
    
    // double totnom = hErecoNom->Integral();
    // double tot    = hErecoTun->Integral();    
    // double tot_p1  = fabs( hErecoTun_p1sigma[ip]->Integral() );
    // double tot_m1  = fabs( hErecoTun_m1sigma[ip]->Integral() );
    // double tot_p3  = fabs( hErecoTun_p3sigma[ip]->Integral() );
    // double tot_m3  = fabs( hErecoTun_m3sigma[ip]->Integral() );
    double totnom = hErecoNom->GetEntries();
    double tot    = hErecoTun->GetEntries();    
    double tot_p1  = fabs( hErecoTun_p1sigma[ip]->GetEntries() );
    double tot_m1  = fabs( hErecoTun_m1sigma[ip]->GetEntries() );
    double tot_p3  = fabs( hErecoTun_p3sigma[ip]->GetEntries() );
    double tot_m3  = fabs( hErecoTun_m3sigma[ip]->GetEntries() );
    double frac_p1 = fabs( (tot_p1 - tot)/tot );
    double frac_m1 = fabs( (tot_m1 - tot)/tot );
    double frac_p3 = fabs( (tot_p3 - tot)/tot );
    double frac_m3 = fabs( (tot_m3 - tot)/tot );
    
    hVarFrac_p1s->SetBinContent(iparhisto+1,frac_p1);
    hVarFrac_m1s->SetBinContent(iparhisto+1,frac_m1);
    hVarTot_p1s->SetBinContent(iparhisto+1,tot_p1);
    hVarTot_m1s->SetBinContent(iparhisto+1,tot_m1);
    
    hVarFrac_p3s->SetBinContent(iparhisto+1,frac_p3);
    hVarFrac_m3s->SetBinContent(iparhisto+1,frac_m3);
    hVarTot_p3s->SetBinContent(iparhisto+1,tot_p3);
    hVarTot_m3s->SetBinContent(iparhisto+1,tot_m3);
    
    iparhisto++;
  }//ip

  string config_file_str(config_file);
  char fname[200];
  if( config_file_str.find("prefit") == string::npos ) sprintf(fname,"%s","_postfit"); // prefit
  else sprintf(fname,"%s","_prefit"); // postfit
  //string file_mode = (update_rootfile ? "UPDATE" : "RECREATE");
  TFile * fouterr = new TFile(Form("fouterr%s_%s.root", fname, sk_sample),file_mode.c_str());
  fouterr->cd();
  hVarFrac_p1s->Write();
  hVarFrac_m1s->Write();
  hVarTot_p1s->Write();
  hVarTot_m1s->Write();
  hVarFrac_p3s->Write();
  hVarFrac_m3s->Write();
  hVarTot_p3s->Write();
  hVarTot_m3s->Write();
  hErecoNom->Write();
  hErecoTun->Write();
  hRatioSum_p1sigma->Write();
  hRatioSum_m1sigma->Write();
  fouterr->Close();


  // Write the table

  //string pctsymbol  = " \\% " ;
  string pctsymbol  = " {(\\si{\\percent})} " ;
  string endcolumn  = " & "   ;
  string endrow     = " \\\\ ";

  int ncol1   = 80;
  int ncolnev = 15;
  int ncolfrc =  8;

  if(tablesyst_idx==-999){  // full table merged modes, all systematics
    
    // write document header to outfile
    
    string output_file (Form("systematic_tables%s_%s.tex",tag,sk_sample));
    
    stringstream streamer;
    // streamer << "\\begin{table}"   << endl;
    // streamer << "\\footnotesize" << endl;
    // streamer << "\\center" << endl;
    // streamer << "\\begin{tabular}{| l | r | r | c | r | r  r | r  r | r |} " << endl;
    streamer << "\\rowcolors{2}{gray!10}{white} " << endl;
    streamer << "\\begin{tabular}{| l | S | S  S | S  S | S  S | S  S | } " << endl; // [input-symbols=()]
    // streamer << "\\hline \\hline" << endl;
    streamer << "\\toprule" << endl;
    streamer << "{\\textbf{ Parameter }} &  "
	     << "\\multicolumn {1}{|c|}{\\textbf{ 1 $\\sigma$ value}} &  "
	     << "\\multicolumn {2}{S|} { {\\textbf{ BANFF tuned -3$\\sigma$}} } & "
	     << "\\multicolumn {2}{S|} { {\\textbf{ BANFF tuned -1$\\sigma$}} } &  "
	     << "\\multicolumn {2}{S|} { {\\textbf{ BANFF tuned +1$\\sigma$}} } &  "
	     << "\\multicolumn {2}{S|} { {\\textbf{ BANFF tuned +3$\\sigma$}} } \\\\  "
	     << endl;
    streamer << "\\midrule" << endl;
    
    // double totnom = hErecoNom->Integral();
    // double tot    = hErecoTun->Integral();    
    // double tot_all_p = hErecoTun_p1sigma_all->Integral();
    // double tot_all_m = hErecoTun_m1sigma_all->Integral();
    double totnom = hErecoNom->GetEntries();
    double tot    = hErecoTun->GetEntries();
    
    double tot_all_p  = hErecoTun_p1sigma_all->GetEntries();
    double tot_all_m  = hErecoTun_m1sigma_all->GetEntries();
    double tot_all_p3 = hErecoTun_p3sigma_all->GetEntries();
    double tot_all_m3 = hErecoTun_m3sigma_all->GetEntries();

    double frac_all_p  = (tot_all_p - tot)/tot * 100;
    double frac_all_m  = (tot_all_m - tot)/tot * 100;
    double frac_all_p3 = (tot_all_p3 - tot)/tot * 100;
    double frac_all_m3 = (tot_all_m3 - tot)/tot * 100;
    
    
    streamer   << setw(ncol1)  << setfill(' ')
	       << "All syst"
	       << endcolumn
	       << setw(ncolnev)  << setfill(' ') << setprecision(5) << ""
	       << endcolumn
	       << setw(ncolnev)  << setfill(' ') << setprecision(6) << tot_all_m3
	       << endcolumn
	       << setw(ncolnev)  << setfill(' ') << setprecision(6) << frac_all_m3  << pctsymbol
	       << endcolumn     
	       << setw(ncolnev)  << setfill(' ') << setprecision(6) << tot_all_m
	       << endcolumn
	       << setw(ncolnev)  << setfill(' ') << setprecision(5) << frac_all_m  << pctsymbol
	       << endcolumn
	       << setw(ncolnev)  << setfill(' ') << setprecision(5) << tot_all_p
	       << endcolumn
	       << setw(ncolnev)  << setfill(' ') << setprecision(5) << frac_all_p << pctsymbol 
	       << endcolumn
	       << setw(ncolnev)  << setfill(' ') << setprecision(6) << tot_all_p3
	       << endcolumn
	       << setw(ncolnev)  << setfill(' ') << setprecision(6) << frac_all_p3 << pctsymbol
	       << endrow
	       << endl;
    streamer << "\\cmidrule(lr){1-1} \\cmidrule(lr){2-2} \\cmidrule(lr){3-4} \\cmidrule(lr){5-6} \\cmidrule(lr){7-8} \\cmidrule(lr){9-10}" << endl;
    
    
    for(int ip = 0; ip < np_loop; ip++) {
      conv::PhysParamId_t p = fparm->WhichParam(ip);
      if(!fparm->IsNuisance(p)) continue;
      string sysstr = util::AsStringCompact(p); // if before fparm->IsNuisance(p) gives problems !!
      if(sysstr=="f_rpa")       continue;  // no tweaks for RPA, which is a discrete parameter

      // double tot_p1  = hErecoTun_p1sigma[ip]->Integral();
      // double tot_m1  = hErecoTun_m1sigma[ip]->Integral();
      // double tot_p3  = hErecoTun_p3sigma[ip]->Integral();
      // double tot_m3  = hErecoTun_m3sigma[ip]->Integral();
      double sigma_val = 0.;
      int idx = -999;
      if ( util::IsSKDetFSI(p) ) {
	idx = skexp.GetSuperkDetectorFSI()->Param2Idx(p);
	sigma_val = sqrt(skexp.GetSuperkDetectorFSI()->CovMatrix()(idx,idx));
      }
      if ( util::IsBANFF(p) )    {
	idx = skexp.GetBANFF()->Param2Idx(p);
	sigma_val = sqrt(skexp.GetBANFF()->CovMatrix()(idx,idx)) ;
      }
      double tot_p1  = hErecoTun_p1sigma[ip]->GetEntries();
      double tot_m1  = hErecoTun_m1sigma[ip]->GetEntries();
      double tot_p3  = hErecoTun_p3sigma[ip]->GetEntries();
      double tot_m3  = hErecoTun_m3sigma[ip]->GetEntries();
      double frac_p1 = (tot_p1 - tot)/tot * 100;
      double frac_m1 = (tot_m1 - tot)/tot * 100;        
      double frac_p3 = (tot_p3 - tot)/tot * 100;
      double frac_m3 = (tot_m3 - tot)/tot * 100;

      if (frac_p1 == 0 && frac_m1 == 0 && frac_p3 == 0 && frac_m3 == 0) continue;
        
      streamer   << setw(ncol1)  << setfill(' ')
	//<< util::AsStringPaperTex(util::PhysParamIdFromString(sysstr))
		 << util::AsStringTex(util::PhysParamIdFromString(sysstr))
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(6) << sigma_val
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(6) << tot_m3
		 << endcolumn 
		 << setw(ncolnev)  << setfill(' ') << setprecision(6) << frac_m3
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(6) << tot_m1
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(5) << frac_m1 
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(6) << tot_p1
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(5) << frac_p1 
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(6) << tot_p3
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(6) << frac_p3
		 << endrow
		 << endl;
        
      if(sysstr=="f_banffnorm_nuetonumu")
	streamer << "\\midrule[0.02em]" << endl;
    }//ip
    streamer << "\\bottomrule" << endl;
    streamer << "\\end{tabular}" << endl;
    // streamer << "\\caption{$\\pm 1,3 \\sigma$ effect of the systematic parameters on the total number of events.}"
    //          << endl;
    // streamer << "\\end{table}" << endl;
    //streamer << "\\clearpage" << endl;
    cout << streamer.str() << endl;
    
    // get size of file
    long size;
    streamer.seekg(0,ifstream::end);
    size=streamer.tellg();
    streamer.seekg(0);
    ofstream outfile (output_file.c_str(),ios_base::app);
    // write to outfile
    outfile.write (streamer.str().c_str(),size);
    outfile.close();
  }//tablesyst_idx==-999


  // Write mode table for a single systematic

  else if(tablesyst_idx==-99 || tablesyst_idx>0){
    
    //output file name
    string output_file_mode (Form("systematic_tables-%s-%s-AllModes.tex",oscparams.c_str(),sk_sample));
    
    // write document header to outfile
    ofstream outfile (output_file_mode.c_str());
    outfile.close();
    
    
    stringstream streamer;
    
    // produce a table for each systematic parameter
    
    for(int ip = 0; ip < np_loop; ip++) {
        
      if( (tablesyst_idx != -99) && (tablesyst_name == "") ){
	cout << endl;
	cout << "You chose the single systematic mode table w/o selecting the paramter" << endl;
	cout << endl;
	exit(1);
      }
        
      conv::PhysParamId_t p = fparm->WhichParam(ip);
      if(!fparm->IsNuisance(p)) continue;
      string sysstr = util::AsStringCompact(p); //if before fparm->IsNuisance(p) gives problems !!!
      if(sysstr=="f_rpa")       continue;  // no tweaks for RPA, which is a discrete parameter
      if( (tablesyst_idx != -99) && (tablesyst_name != sysstr) ) continue;
        
        
      cout << "Make table " << sysstr << endl;
        
      streamer << "\\begin{table}"   << endl;
      streamer << "\\tiny" << endl;
      streamer << "\\center" << endl;
      streamer << "\\begin{tabular}{| l | r | r | r | r  r | r  r | r |}" << endl;
      streamer << "\\hline \\hline" << endl;
      streamer << "{\\textbf{ Component }} &  "
	       << "{\\textbf{ Nominal }}   &  "
	       << "{\\textbf{ BANFF tuned }} &  "
	       << "{\\textbf{ BANFF tuned -3$\\sigma$ }}  &  "
	       << "\\multicolumn {2}{c|} { {\\textbf{ BANFF tuned -1$\\sigma$}} } &  "
	       << "\\multicolumn {2}{c|} { {\\textbf{ BANFF tuned +1$\\sigma$}} } &  "
	       << "{\\textbf{ BANFF tuned +3$\\sigma$}}  \\\\  "
	       << endl;
      streamer << "\\hline" << endl;
        
      for (int im = 0; im < conv::knm; im++) {
	conv::SuperkMode_t m = conv::kSKReactionModes[s][im];
	if(m == conv::kSKmUndefined) continue;
	if(util::ArrIdx::Instance()->SuperkModeIdx(s,m) < 0) continue;
	string modestr = util::AsStringCompact(m);
            
	// double totnom = hErecoNomMode[im]->Integral();
	// double tot = hErecoTunMode[im]->Integral();
	// double tot_p1 = hErecoTunMode_p1sigma[ip][im]->Integral();
	// double tot_m1 = hErecoTunMode_m1sigma[ip][im]->Integral();
	double totnom = hErecoNomMode[im]->GetEntries();
	double tot    = hErecoTunMode[im]->GetEntries();
	double tot_p1 = hErecoTunMode_p1sigma[ip][im]->GetEntries();
	double tot_m1 = hErecoTunMode_m1sigma[ip][im]->GetEntries();
	// double tot_p3 = hErecoTunMode_p3sigma[ip][im]->GetEntries();
	// double tot_m3 = hErecoTunMode_m3sigma[ip][im]->GetEntries();
	double frac_p1 = (tot_p1 - tot)/tot * 100;
	double frac_m1 = (tot_m1 - tot)/tot * 100;
            
	streamer << setw(ncol1)  << setfill(' ')
		 << util::AsStringTex(m)
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(6) << totnom
		 << endcolumn
		 << tot
		 << endcolumn
		 << tot_p1
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(5) << frac_p1 << pctsymbol
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(6) << tot_m1
		 << endcolumn
		 << setw(ncolnev)  << setfill(' ') << setprecision(5) << frac_m1 << pctsymbol
	  // << endcolumn     
	  // << setw(ncolnev)  << setfill(' ') << setprecision(6) << tot_p3
	  // << endcolumn     
	  // << setw(ncolnev)  << setfill(' ') << setprecision(6) << tot_m3 
		 << endrow
		 << endl;
            
	if(modestr=="numuncoth"    ||
	   modestr=="numubarncoth" ||
	   modestr=="nuencoth"     ||
	   modestr=="oscnueccoth"  ||
	   modestr=="nuebarncoth" )
	  streamer << "\\hline" << endl;
            
            
      } // loop over modes //im

      streamer << "\\hline" << endl;
      streamer << "\\end{tabular}" << endl;
      streamer << "\\caption{$\\pm 1 \\sigma$ effect of "
	       << util::AsStringPaperTex(p) 
	       << ".}"
	       << endl;
      streamer << "\\end{table}" << endl;
      streamer << "\\clearpage" << endl;
        
    } // end loop over systematics ip
    
    cout << streamer.str() << endl;
    
    
    // get size of file
    long size;
    streamer.seekg(0,ifstream::end);
    size=streamer.tellg();
    streamer.seekg(0);
    //ofstream outfile (output_file_mode.c_str(),ios_base::app);
    outfile.open(output_file_mode.c_str(),ios_base::app);
    // write to outfile
    outfile.write (streamer.str().c_str(),size);
    outfile.close();
  }//(tablesyst_idx==-99 || tablesyst_idx>0)
  else{
    cout << endl;
    cout << "You didn't choose any table option" << endl;
    cout << "tablesyst_idx = " << tablesyst_idx << endl;
    cout << endl;      
  } // end if table



  return;
}


