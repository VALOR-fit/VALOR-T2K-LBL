#include <iomanip.h>
#include <TXMLEngine.h>

void test_systematic_error_propagation_all_2015
(
  int          mode               = 0,                                         /* options: see above */
  bool         update_rootfile    = false,                                        /* "update" ROOT output file? If false, will "recreate" */
  const char * sk_sample          = "RHC1Rmu",                                   /* options: "1Rmu", "1Re", "RHC1Rmu", "RHC1Re", "MultiRmu" */
  const char * hypothesis         = "3f_osc",                                    /* options: "2f_numu_disap", "3f_osc", "3f_oscbar" */
  double       input_mixing12     = -99,                                         /* */
  double       input_mixing13     = -99,                                         /* */
  double       input_mixing23     = -99,                                         /* */
  double       input_deltacp      = -99,                                         /* */
  double       input_dm12sq       = -99,                                         /* */
  double       input_dm23sq       = -99,                                         /* */
  bool         save_toy_expts     = false,                                       /* save the ereco spectra? */
  const char * config_file        = "analysis_configuration_run1to6_test.xml",   /* */
  //const char * fit_param_list     = "nuosc3fbar_fixnonbar_12fix+nuisance", 
  const char * fit_param_list     = "nuosc3f_12fix+nuisance",
  const char * mixparamchoice     = "sinsqtheta",
  const char * hierarchy          = "normal"
)
{
  using namespace std;
  using namespace valor;
  using namespace valor::conv;

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
  conv::SuperkSample_t s = util::SKSampleFromString(sk_sample); // requested sample
  if(s == conv::kSKsUndefined) {
    cerr << "** Undefined sample!" << endl;
    return;
  }
  int is = util::ArrIdx::Instance()->SuperkSampleIdx(s); // index of specified sample
  
  // Analysis configuration
  const char * config  = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);
  T2KOpt::Instance()->ReadXmlPath(config);
    
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
  if(input_dm23sq < -1)  input_dm23sq = conv::kGlobFitDm23sq;
  if(input_dm12sq < -1)  input_dm12sq = conv::kGlobFitDm12sq;
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
  cout << "- Reading analysis configuration from: " << config << endl;
  cout << "- Physics hypothesis: \n" << ph << endl;
  cout << *Opt::Instance() << endl;
  
 //output file name
  string oscparams = 
    util::AsString(ph.Hypothesis(),input_mixing12,input_mixing13,input_mixing23,input_deltacp,input_dm12sq,input_dm23sq);
  string output_file (Form("systematic_tables-%s-%s.tex",oscparams.c_str(),sk_sample));

  // write document header to outfile
  ofstream outfile (output_file.c_str());
  string f;
  f += "\\documentclass[10pt,a4paper]{article}\n";
  f += ("\\usepackage{amsmath}\n");
  
  f += ("\\usepackage{graphicx}\n");
  f += ("\\usepackage{dcolumn}\n");
  f += ("\\usepackage[cm]{fullpage}\n");
  f += ("\\begin{document}\n");
  outfile.write (f.c_str(),f.size());
  outfile.close();

  // Precison for printed-out numerical results
  int nprecisionnev = 4; // <-- for number of events
  int nprecisionpct = 4; // <-- for % changes

  // Create & config superk spectrum predictor
  SuperkExpectation skexp;
  skexp.Read(config);
  SuperkObservation skobs;

  // Define the parameter list
  FitParamList * fparm = util::NamedParameterList(fit_param_list); 
  const int np = fparm->GetNumberParameters();
  //fparm->Print();

  // Define the tweak histograms
  TH1D *hErecoTun_p1sigma[np];
  TH1D *hErecoTun_m1sigma[np];

  // Get tuned MC - with BANFF tuning
  skexp.GetBANFF()->RestoreTunedMC();
  skexp.GetBANFF()->RestoreRHCTunedMC();
  skexp.Generate(ph);
  TH1D *hErecoTun = new TH1D(*skexp.CurrPredictedErecoSpectrum(s));
  hErecoTun->SetName("hErecoTun");

  // Get nominal MC - No BANFF tuning
  skexp.GetBANFF()->RestoreNominalMC();
  skexp.GetBANFF()->RestoreRHCNominalMC();
  skexp.Generate(ph);
  TH1D *hErecoNom = new TH1D(*skexp.CurrPredictedErecoSpectrum(s));
  hErecoNom->SetName("hErecoNom");

  // Loop over systematics  
  for(int ip = 0; ip < np; ip++) {
    conv::PhysParamId_t p = fparm->WhichParam(ip);
    if(!fparm->IsNuisance(p)) continue; 
    string sysstr = util::AsStringCompact(p); // if before fparm->IsNuisance(p) gives problems !!!

    if(sysstr=="f_rpa")       continue;  // no tweaks for RPA, which is a discrete parameter

    cout << endl;
    cout << "Tweak systematic parameter: " << sysstr << endl;

    int idx=-999;
    double tweak_p=0.;
    double tweak_m=0.;
    double bestfit=0.;

    if(util::IsRHCSystematic(p)){
      if(util::IsSKDetFSI(p)){
	idx = util::ArrIdx::Instance()->RHCSKDetFSISystIdx(p);
      	if(idx>=0){
	  bestfit = 1.;
	  tweak_p = bestfit + (1 * sqrt(skexp.GetSuperkDetectorFSI()->CovMatrix()(idx,idx)));
	  tweak_m = bestfit - (1 * sqrt(skexp.GetSuperkDetectorFSI()->CovMatrix()(idx,idx)));
	}
      }
      else if(util::IsBANFF(p)){
	idx = util::ArrIdx::Instance()->RHCBANFFSystIdx(p);
     	if(idx>=0){
	  bestfit = skexp.GetBANFF()->RHCBestFitValue(p);
	  tweak_p = bestfit + (1 * sqrt(skexp.GetBANFF()->CovMatrix()(idx,idx)));
	  tweak_m = bestfit - (1 * sqrt(skexp.GetBANFF()->CovMatrix()(idx,idx)));
	}
      }
      else if(util::IsXSec(p)){
	if(idx>=0){
	  bestfit = 0;
	  tweak_p = skexp.GetNuXSec()->Tweak(kRHCNuisXSecRPA , 1.);
	  tweak_m = skexp.GetNuXSec()->Tweak(kRHCNuisXSecRPA , 0.);
	}
      }
      else{
	cout << endl;
	cout << "RHC systematic parameter not found!!!" << endl;
	exit(1);
      }
    }
    else if(util::IsFHCSystematic(p)){
      if(util::IsSKDetFSI(p)){
    	idx = util::ArrIdx::Instance()->SKDetFSISystIdx(p);     
    	if(idx>=0){
	  bestfit = 1;
    	  tweak_p = bestfit + (1 * sqrt(skexp.GetSuperkDetectorFSI()->CovMatrix()(idx,idx)));
    	  tweak_m = bestfit - (1 * sqrt(skexp.GetSuperkDetectorFSI()->CovMatrix()(idx,idx)));
    	}
      }
      else if(util::IsBANFF(p)){
    	idx = util::ArrIdx::Instance()->BANFFSystIdx(p);
    	if(idx>=0){
	  bestfit = skexp.GetBANFF()->BestFitValue(p);
    	  tweak_p = bestfit + (1 * sqrt(skexp.GetBANFF()->CovMatrix()(idx,idx)));
    	  tweak_m = bestfit - (1 * sqrt(skexp.GetBANFF()->CovMatrix()(idx,idx)));
    	}   
      }
      else if(util::IsXSec(p)){
      	if(idx>=0){
	  bestfit = 0;
    	  tweak_p = skexp.GetNuXSec()->Tweak(kNuisXSecRPA , 1.);
    	  tweak_m = skexp.GetNuXSec()->Tweak(kNuisXSecRPA , 0.);
    	}
      }
      else{
    	cout << endl;
    	cout << "FHC systematic parameter not found!!!" << endl;
    	exit(1);
      }
    }

    // sk_e_scale is treated in a different way: not in the cov. matrix so no SystIdx() --> idx = -1 
    // RHC/FHC sample is treated inside the function
    if(sysstr.find("sk_e_scale") && idx<0){ 
      tweak_p = 1 + skexp.GetSuperkDetectorFSI()->EnergyScaleRelErr(s,kSKmNuMuCCQE);
      tweak_m = 1 - skexp.GetSuperkDetectorFSI()->EnergyScaleRelErr(s,kSKmNuMuCCQE);
    }

    // +1 sigma tweak systematic
    skexp.GetBANFF()->RestoreTunedMC();
    skexp.GetBANFF()->RestoreRHCTunedMC();
    skexp.TweakSystematic(sysstr,tweak_p);
    skexp.Generate(ph);
    hErecoTun_p1sigma[ip] = new TH1D(*skexp.CurrPredictedErecoSpectrum(s));
    char name[200];
    sprintf(name,"hErecoTun_p1sigma_%s",sysstr.c_str());
    hErecoTun_p1sigma[ip]->SetName(name);

    // +1 sigma tweak systematic
    skexp.GetBANFF()->RestoreTunedMC();
    skexp.GetBANFF()->RestoreRHCTunedMC();
    skexp.TweakSystematic(sysstr,tweak_m);
    skexp.Generate(ph);
    hErecoTun_m1sigma[ip] = new TH1D(*skexp.CurrPredictedErecoSpectrum(s));
    char name[200];
    sprintf(name,"hErecoTun_m1sigma_%s",sysstr.c_str());
    hErecoTun_m1sigma[ip]->SetName(name);
    
    //tot = hErecoTun->Integral();
    //tot_p = hErecoTun_p1sigma[ip]->Integral();
    //tot_m = hErecoTun_m1sigma[ip]->Integral();
  }

  // Write the table 

  string pctsymbol  = " \\% " ;
  string endcolumn  = " & "   ;
  string endrow     = " \\\\ ";
  
  int ncol1   = 80;
  int ncolnev = 15;
  int ncolfrc =  8;

  stringstream streamer;
  streamer << "\\begin{table}"   << endl;
  //streamer << "\\centering"     << endl;
  streamer << "\\tiny" << endl;
  streamer << "\\begin{tabular}{| l | r | r | r  r | r  r |}" << endl;
  streamer << "\\hline \\hline" << endl;
  streamer << "{\\bf Component} &  "
  	   << "{\\bf Nominal}   &  "
  	   << "{\\bf BANFF tuned} &  "
    //<< "{\\bf BANFF tuned +1$\\sigma$} &  "
    //<< "{\\bf BANFF tuned -1$\\sigma$} &  "
	   << "\\multicolumn {2}{c|} { {\\bf BANFF tuned +1$\\sigma$} } &  "
	   << "\\multicolumn {2}{c|} { {\\bf BANFF tuned -1$\\sigma$} } \\\\"
  	   << endl;
  streamer << "\\hline" << endl;

  for(int ip = 0; ip < np; ip++) {
    conv::PhysParamId_t p = fparm->WhichParam(ip);
    if(!fparm->IsNuisance(p)) continue; 
    string sysstr = util::AsStringCompact(p); // if before fparm->IsNuisance(p) gives problems !!
    if(sysstr=="f_rpa")       continue;  // no tweaks for RPA, which is a discrete parameter

    double totnom = hErecoNom->Integral();
    double tot    = hErecoTun->Integral();
    double tot_p  = hErecoTun_p1sigma[ip]->Integral();
    double tot_m  = hErecoTun_m1sigma[ip]->Integral();    
    double frac_p = (tot_p - tot)/tot * 100;
    double frac_m = (tot_m - tot)/tot * 100;

    streamer         << setw(ncol1)  << setfill(' ')
		     << util::AsStringPaperTex(util::PhysParamIdFromString(sysstr))
		     << endcolumn
		     << setw(ncolnev)  << setfill(' ') << setprecision(6) << totnom
		     << endcolumn
		     << tot
		     << endcolumn
		     << tot_p 
		     << endcolumn
		     << setw(ncolnev)  << setfill(' ') << setprecision(5) << frac_p << pctsymbol
		     << endcolumn
		     << setw(ncolnev)  << setfill(' ') << setprecision(6) << tot_m
		     << endcolumn
		     << setw(ncolnev)  << setfill(' ') << setprecision(5) << frac_m << pctsymbol
		     << endrow
		     << endl; 
  
    if(sysstr=="f_banffnorm_nuetonumu_rhc")
      streamer << "\\hline" << endl;    
  }
  streamer << "\\hline" << endl;
  streamer << "\\end{tabular}" << endl;
  streamer << "\\end{table}" << endl;
  streamer << "\\clearpage" << endl;
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
  
  ofstream outfile (output_file.c_str(),ios_base::trunc&&ios_base::out);
  string f;
  f += ("\\end{document}\n");
  outfile.write (f.c_str(),f.size());
  outfile.close();

  return;
}
