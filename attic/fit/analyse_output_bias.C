// 
// Analyse the output of run_fit.C: 
//
// - Used for many toys study
// - Generated toy parameters for a single physics hypothesis
// - Best fit parameters and fit results
// 
// Produce the following output:
// 
// - pulls, best-fit, chisq
//


#define DrawOscFit false
#define SaveOscFit false

#define DrawOscPull false
#define SaveOscPull false

#define DrawOscRes false
#define SaveOscRes false

#define DrawNuisFit true
#define SaveNuisFit false

#define DrawNuisPull false
#define SaveNuisPull false

#define DrawNuisRes false
#define SaveNuisRes false

#define DrawNuisExtr false
#define SaveNuisExtr false

#define DrawOscErr false
#define DrawNuisErr false

#define NBinsFitPar 200
#define NBinsPull 200

#define SaveCombinedPullPlots false

using namespace std;
using namespace valor;
using namespace valor::conv;


// Declare Histograms

TString *NameFitParAll;

TH1F **hExtrParAll;
TH1F **hFitParAll;
TH1F **hErrParAll;
TH1F **hErrPlusParAll;
TH1F **hErrMinusParAll;

TH1F **hPullAll;
TH1F **hPullValAll;
TH1F **hResAll;
TH1F *hPullMean;
TH1F *hPullWidth;
TH1F *hPullValMean;
TH1F *hPullValWidth;

TCanvas **cExtrParAll;
TCanvas **cFitParAll;
TCanvas **cErrParAll;
TCanvas **cErrMinosParAll;
TCanvas **cPullAll;
TCanvas **cResAll;

char name[200];
char canvname[200];

void save_canvas(TCanvas * c, const char * name, const char * tag) {
  c->SaveAs(Form("%s_%s.pdf",tag,name));
  c->SaveAs(Form("%s_%s.eps",tag,name));
  c->SaveAs(Form("%s_%s.C",  tag,name));
  c->SaveAs(Form("%s_%s.root",  tag,name));
}

TLatex *DrawName(char *name,double x,double y){
  TLatex *tname = new TLatex(x,y,name);
  tname->SetTextSize(0.08);
  tname->Draw();
  //tname->DrawClone("");
}

double calc_pull_osc(double centval,double fitval,double fiterr){
  double pull = (fitval-centval)/fiterr; 
  return pull;
}

double calc_res(double centval,double fitval){
  double res = (fitval-centval)/fitval; 
  return res;
}

double calc_pull_nuis(double centval,double centerr,double fitval,double fiterr){
  double pull = (fitval-centval)/sqrt(centerr*centerr - fiterr*fiterr); 
  //double pull = (fitval-centval)/centerr; 
  //double pull = (fitval-centval)/fiterr; 
  //cout << fitval << " " << centval << " " << centerr << endl;
  return pull;
}

double calc_pull_val(double centval,double centerr,double fitval,double fiterr){ 
  double pull = (fitval-centval)/fiterr; 
  return pull;
}
 
inline double get_cent_err(conv::PhysParamId_t p,SuperkExpectation * skexp){ // very slow...

  double centerr = 0.;
  
  if(p == conv::kPhypMixing13) {      
    if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsq2Theta){
      centerr = conv::kGlobFitErrSinsq2Theta13;
    }
    else if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsqTheta){
      centerr = conv::kGlobFitErrSinsqTheta13;
    }
  }
  else if(p == conv::kPhypMixing12) {
    if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsq2Theta){
      centerr = conv::kGlobFitErrSinsq2Theta12;
    }
    else if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsqTheta){
      centerr = conv::kGlobFitErrSinsqTheta12;
    }
  }
  else if(p == conv::kPhypMixing23) {
    if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsq2Theta){
      centerr = conv::kGlobFitErrSinsq2Theta23;
    }
    else if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsqTheta){
      centerr = conv::kGlobFitErrSinsqTheta23;
    }
  }
  else if(p == conv::kPhypDm23sq) {
    centerr = conv::kGlobFitErrDm23sq;
  }
  else if(p == conv::kPhypDm12sq) {
    centerr = conv::kGlobFitErrDm12sq;
  }
  else if(p == conv::kPhypCPV) {
    centerr = conv::kGlobFitErrDeltaCP;
  }
  
  else if(util::IsBANFF(p)){
    //centval = skexp->GetBANFF()->BestFitValue(p);
    centerr = skexp->GetBANFF()->OneSigmaError(p);
  }	  
  else if(util::IsXSec(p)) {
    //"nominal best fit" is 0
    //"1 sigma error" is 1
    //centval = 0;
    centerr = 1;
  }
  else if(p == conv::kNuisSuperkEnergyScale) {
    //"nominal best fit" is 0
    //"1 sigma error" is 1
    //centval = 0;
    centerr = 1;
  }
  else if(util::IsSKDetFSI(p)) {
    //"nominal best fit" is 1
    //centval = 1;
    centerr = skexp->GetSuperkDetectorFSI()->OneSigmaError(p);
  }	
  else {
    cerr << "Unknown pull for parameter " << util::AsStringCompact(p) << endl;
  }
  return centerr;
}




void analyse_output_bias
(
 const char * infile,
 const char * tag,
 const char * hierarchy = "normal",
 const char * dmsqscheme = "abs23",
 const char * mixparamchoice = "sinsqtheta",
 //const char * config_file = "analysis_configuration_test.xml",
 const char * config_file = "analysis_configuration_test_multiring.xml",
 const char * hypothesis = "3f_osc",   // "3f_osc", "4f_sterile_osc", "5f_sterile_osc 
 )
{  
  // Global options
  Opt::Instance()->gMassHierarchyAssumed = util::HierarchyFromString(hierarchy);
  Opt::Instance()->gOscDmsqScheme        = util::DmsqSchemeFromString(dmsqscheme);
  Opt::Instance()->gOscMixingParamChoice = util::MixingTypeFromString(mixparamchoice);
  if(Opt::Instance()->Undefined()) {
    cerr << "Global options not defined correctly. Exiting..." << endl;
    cerr << *Opt::Instance() << endl;
    exit(1);
  }

 Opt::Instance()->SetDefaultParameterValues();
 
  // Analysis configuration
  const char * config = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);
  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read(config);
  PhysHypothesis ph(hypothesis); //initialisation is to no oscillations
  skexp->Generate(ph);

  // Arrays for storing all pieces of fit data and the best-fit pdf (and its components) for each toy-MC.
  
  TObjArray * data = new TObjArray; 
  data->SetOwner(false);
  TObjArray * best_fit_pdf = new TObjArray; 
  best_fit_pdf->SetOwner(false);
  TObjArray * best_fit_pdf_fit1 = new TObjArray;
  best_fit_pdf_fit1->SetOwner(false);
  TObjArray * best_fit_pdf_fit2 = new TObjArray;
  best_fit_pdf_fit2->SetOwner(false);
  int fitstatus1, fitstatus2, fitstatus;
  int fitstatus_err1, fitstatus_err2, fitstatus_err;
  double chisq1, chisq2, chisq;

  // Tree variables
  PhysParamList * input_params; //only used to initialise branch
  FitParamList * best_fit_params; //only used to initialise branch
  FitParamList * best_fit_params_fit1; //first fit (first octant)
  FitParamList * best_fit_params_fit2; //second fit (second octant)

  TFile *filein = new TFile(infile,"READ");
  TTree *fitnt = (TTree*) filein->Get("fitnt");
  
  // Create tree to get outputs
  fitnt->SetBranchAddress("input_params",&input_params);
  fitnt->SetBranchAddress("best_fit_params",&best_fit_params);
  fitnt->SetBranchAddress("best_fit_params_fit1",&best_fit_params_fit1);
  fitnt->SetBranchAddress("best_fit_params_fit2",&best_fit_params_fit2);
  fitnt->SetBranchAddress("fitstatus1",&fitstatus1);
  fitnt->SetBranchAddress("fitstatus2",&fitstatus2);
  fitnt->SetBranchAddress("fitstatus",&fitstatus);
  fitnt->SetBranchAddress("fitstatus_err1",&fitstatus_err1);
  fitnt->SetBranchAddress("fitstatus_err2",&fitstatus_err2);
  fitnt->SetBranchAddress("fitstatus_err",&fitstatus_err);
  fitnt->SetBranchAddress("chisq1",&chisq1);
  fitnt->SetBranchAddress("chisq2",&chisq2);
  fitnt->SetBranchAddress("chisq",&chisq);

  // Initialize Histograms

  cout << endl;
  cout << "Initialize Histograms" << endl;
  cout << endl;
  
  fitnt->GetEntry(0); // Get Tree general informations
  int np = best_fit_params->GetNumberParameters();
  int nph = input_params->NParams();

  NameFitParAll = new TString[np];
  cExtrParAll = new TCanvas* [np];
  cFitParAll = new TCanvas* [np];
  cErrParAll = new TCanvas* [np];
  cErrMinosParAll = new TCanvas* [np];
  cPullAll = new TCanvas* [np];
  cPullOsc = new TCanvas* [np];
  cResAll = new TCanvas* [np];
  cResOsc = new TCanvas* [np];
  hPullAll = new TH1F*[np];
  hPullValAll = new TH1F*[np];
  hPullOsc = new TH1F*[np];
  hPullNuis = new TH1F*[np];
  hResAll = new TH1F*[np];
  hResOsc = new TH1F*[np];
  hResNuis = new TH1F*[np];
  hFitParAll = new TH1F*[np];
  hFitParOsc = new TH1F*[np];
  hFitParNuis = new TH1F*[np];
  hExtrParAll = new TH1F*[np];
  hErrParAll = new TH1F*[np];
  hErrPlusParAll = new TH1F*[np];
  hErrMinusParAll = new TH1F*[np];
  
  hPullMean = new TH1F("hPullMean","hPullMean",np,0,np);
  hPullWidth = new TH1F("hPullWitdh","hPullWidth",np,0,np);
  hPullValMean = new TH1F("hPullValMean","hPullValMean",np,0,np);
  hPullValWidth = new TH1F("hPullValWitdh","hPullValWidth",np,0,np);

  double max=-999.;
  double min=-999.;
  for(int ip = 0; ip < np; ip++) {
    conv::PhysParamId_t p = best_fit_params->WhichParam(ip);

    // Extracted distributions
    max = best_fit_params->MaxValue(p);
    min = best_fit_params->MinValue(p);
    string extrname = "extr_";
    extrname.append(best_fit_params->Name(p));
    sprintf(name,extrname.c_str());
    hExtrParAll[ip] = new TH1F(name,name,NBinsFitPar,min,max);
    
    // Fitted distributions
    max = best_fit_params->MaxValue(p);
    min = best_fit_params->MinValue(p);
    string fitname = "fit_";
    fitname.append(best_fit_params->Name(p));
    sprintf(name,fitname.c_str());
    hFitParAll[ip] = new TH1F(name,name,NBinsFitPar,min,max);

    // Fitted Error distributions
    max = best_fit_params->MaxValue(p);
    min = 0.;
    
    string errname = "err_";
    errname.append(best_fit_params->Name(p));
    sprintf(name,errname.c_str());
    hErrParAll[ip] = new TH1F(name,name,NBinsFitPar,min,max);

    string errname = "err_plus";
    errname.append(best_fit_params->Name(p));
    sprintf(name,errname.c_str());
    hErrPlusParAll[ip] = new TH1F(name,name,NBinsFitPar,min,max);

    string errname = "err_minus";
    errname.append(best_fit_params->Name(p));
    sprintf(name,errname.c_str());
    hErrMinusParAll[ip] = new TH1F(name,name,NBinsFitPar,min,max);
    
    // Pull distributions
    string pullname = "pull_";
    pullname.append(best_fit_params->Name(p));
    sprintf(name,pullname.c_str());
    hPullAll[ip] = new TH1F(name,name,NBinsPull,-5,5);    
    
    // Pull Value distributions
    string pullname = "pullvalue_";
    pullname.append(best_fit_params->Name(p));
    sprintf(name,pullname.c_str());
    hPullValAll[ip] = new TH1F(name,name,NBinsPull,-5,5);    
    
    // Residual distributions
    string resname = "res_";
    resname.append(best_fit_params->Name(p));
    sprintf(name,resname.c_str());
    hResAll[ip] = new TH1F(name,name,NBinsPull,-5,5);        
  }
  
  // Loop on the tree events
  
  int NEntries = fitnt->GetEntries();  
  for(int itoy=0; itoy<NEntries; itoy++) {
    fitnt->GetEntry(itoy);
    
    if(itoy==0) best_fit_params->Print();    
    //if(itoy==0) cout <<*input_params << endl;

    if(!(itoy%100)) cout << "toy " << itoy << endl;
    

    for(int ip = 0; ip < np; ip++) {
      conv::PhysParamId_t p = best_fit_params->WhichParam(ip);
   
      double centval=-999.;
      double centerr=-999.;
      double pull=-999.;

      double bestfitval = best_fit_params->BestFitValue(p);
      double bestfiterr, bestfiterrP, bestfiterrM;
      if(best_fit_params->HesseErrorsSet())
	bestfiterr = best_fit_params->HesseError(p);
      else{ 
	cout << "No HESSE Error --> take MIGRAD Error" << endl;
	bestfiterr = best_fit_params->MigradError(p);
      }

      if(best_fit_params->MinosErrorsSet()){
	bestfiterrM = TMath::Abs(best_fit_params->MinusError(p));
	bestfiterrP = TMath::Abs(best_fit_params->PlusError(p));
	cout << bestfiterrM << " " << bestfiterrP << endl;
      }
      
      centval = input_params->Value(p);  
      
      cout << centval << endl;

 
      // bool EscaleFlag = true;
      // for(int ipescale = 0; ipescale < np; ipescale++) {
      // 	conv::PhysParamId_t pescale = best_fit_params->WhichParam(ipescale);
      // 	//if(pescale == conv::kNuisSuperkEnergyScale){
      // 	if(pescale == conv::kMultiRingNuisSuperkEnergyScale){
      // 	  bestfitval_escale = best_fit_params->BestFitValue(pescale);
      // 	  //if( (bestfitval_escale<-0.0005 || bestfitval_escale>0.0005) ){
      // 	  EscaleFlag = false;
      // 	  //continue;
      // 	  cout << bestfitval_escale << endl;
      // 	  //}
      // 	}
      // }
      // if(!EscaleFlag) continue;
      
      // Extracted parameters
      hExtrParAll[ip]->Fill(centval);
      
      // Fitted parameters
      hFitParAll[ip]->Fill(bestfitval); 

      // Fitted error parameters
      hErrParAll[ip]->Fill(bestfiterr); 

      // Fitted error parameters
      hErrPlusParAll[ip]->Fill(bestfiterrP); 
      hErrMinusParAll[ip]->Fill(bestfiterrM); 

      // Take the parameter names
      
      if(p == conv::kPhypMixing13) {      
	if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsq2Theta){
	  NameFitParAll[ip] = Form("sin^{2} 2 #theta_{13} = %0.4f",centval);	    
	}
	else if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsqTheta){
	  NameFitParAll[ip] = Form("sin^{2} #theta_{13} = %0.4f",centval);	    
	}
      }
      else if(p == conv::kPhypMixing12) {
	if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsq2Theta){
	  NameFitParAll[ip] = Form("sin^{2} 2 #theta_{12} = %0.1f",centval);	    
	}
	else if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsqTheta){
	  NameFitParAll[ip] = Form("sin^{2} #theta_{12} = %0.1f",centval);	    
	}
      }
      else if(p == conv::kPhypMixing23) {
	if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsq2Theta){
	  NameFitParAll[ip] = Form("sin^{2} 2 #theta_{23} = %0.1f",centval);	    
	}
	else if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsqTheta){
	  NameFitParAll[ip] = Form("sin^{2} #theta_{23} = %0.1f",centval);	    
	}
      }
      else if(p == conv::kPhypDm23sq) {
	NameFitParAll[ip] = Form("|#Delta m^{2}_{23}| = %0.4f",centval);	    
      }
      else if(p == conv::kPhypDm12sq) {
	NameFitParAll[ip] = Form("|#Delta m^{2}_{12}| = %0.1f",centval);	    
      }
      else if(p == conv::kPhypCPV) {
	NameFitParAll[ip] = Form("#delta_{CP} = %0.1f",centval);
      }
            
      // Pulls oscillation paramters
      
      if(!best_fit_params->IsNuisance(p)) { // free  params
	pull = calc_pull_osc(centval,bestfitval,bestfiterr);
	hPullAll[ip]->Fill(pull);
      }	

      else { // nuisance  params
	
	// Oscillation params
	
	if(p == conv::kPhypMixing13) {      
	  if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsq2Theta){
	    centerr = conv::kGlobFitErrSinsq2Theta13;
	  }
	  else if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsqTheta){
	    centerr = conv::kGlobFitErrSinsqTheta13;
	  }
	}
	else if(p == conv::kPhypMixing12) {
	  if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsq2Theta){
	    centerr = conv::kGlobFitErrSinsq2Theta12;
	  }
	  else if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsqTheta){
	    centerr = conv::kGlobFitErrSinsqTheta12;
	  }
	}
	else if(p == conv::kPhypMixing23) {
	  if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsq2Theta){
	    centerr = conv::kGlobFitErrSinsq2Theta23;
	  }
	  else if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsqTheta){
	    centerr = conv::kGlobFitErrSinsqTheta23;
	  }
	}
	else if(p == conv::kPhypDm23sq) {
	  centerr = conv::kGlobFitErrDm23sq;
	}
	else if(p == conv::kPhypDm12sq) {
	  centerr = conv::kGlobFitErrDm12sq;
	}
	else if(p == conv::kPhypCPV) {
	  centerr = conv::kGlobFitErrDeltaCP;
	}
	else if(util::IsBANFF(p)){
	  //centval = skexp->GetBANFF()->BestFitValue(p);
	  centerr = skexp->GetBANFF()->OneSigmaError(p);
	}	  
	else if(util::IsXSec(p)) {
	  //"nominal best fit" is 0
	  //"1 sigma error" is 1
	  //centval = 0;
	  centerr = 1;
	}
	else if(p == conv::kNuisSuperkEnergyScale || p == conv::kMultiRingNuisSuperkEnergyScale) {
	  //"nominal best fit" is 0
	  //"1 sigma error" is 1
	  //centval = 0;
	  centerr = 1;
	}      
	else if(util::IsSKDetFSI(p)) {
	  //"nominal best fit" is 1
	  //centval = 1;
	  if(util::IsMultiRingSystematic(p))
	    centerr = skexp->GetSuperkDetectorFSI()->MultiRingOneSigmaError(p); 
	  else
	    centerr = skexp->GetSuperkDetectorFSI()->OneSigmaError(p);
	}	
	else {
	  cerr << "Unknown pull for parameter " << util::AsStringCompact(p) << endl;
	}
	
	pull = calc_pull_nuis(centval,centerr,bestfitval,bestfiterr);	

	//if(centerr < bestfiterr) cout << centerr << " " << bestfiterr << endl;

	hPullAll[ip]->Fill(pull);
		
	//pull = calc_pull_val(centval,centerr,bestfitval,bestfiterr);
	//hPullValAll[ip]->Fill(pull);

      } // end pulls oscillation params

    } // loop np
    
    if(itoy>2000) break;
    
  } // end itoy
  
  // Build pull mean histograms

  double mean=0.;
  double rms=0.;

  for(int ip = 0; ip < np; ip++) {
    conv::PhysParamId_t p = best_fit_params->WhichParam(ip);
    
    if(!best_fit_params->IsNuisance(p)) continue; // no free params

    mean = hPullAll[ip]->GetMean();
    rms  = hPullAll[ip]->GetRMS();

    hPullMean ->SetBinContent(ip+1,mean);
    hPullWidth->SetBinContent(ip+1,rms);
    
    hPullMean ->GetXaxis()->SetBinLabel(ip+1, util::AsStringCompact(p).c_str());
    hPullWidth->GetXaxis()->SetBinLabel(ip+1, util::AsStringCompact(p).c_str());
    
    //mean = hPullValAll[ip]->GetMean();
    //rms  = hPullValAll[ip]->GetRMS();

    //hPullValMean ->SetBinContent(ip+1,mean);
    //hPullValWidth->SetBinContent(ip+1,rms);
  }
  
  // Draw histograms
  //int np = best_fit_params->GetNumberParameters();
  for(int ip = 0; ip < np; ip++) {
    conv::PhysParamId_t p = best_fit_params->WhichParam(ip);
   
    // Fit parameters
    
    if(!best_fit_params->IsNuisance(p) && !best_fit_params->IsFixed(p) && DrawOscFit ){ // Oscillations    

      string fitname = "canv_fit_";
      fitname.append(best_fit_params->Name(p));
      sprintf(canvname,fitname.c_str());
      cFitParAll[ip] = new TCanvas(canvname,canvname);
      cFitParAll[ip]->cd();

      hFitParAll[ip]->Draw();  
      hFitParAll[ip]->SetLineColor(kBlack);   
      double x = best_fit_params->MinValue(p);
      double y = hFitParAll[ip]->GetMaximum()*0.75;
      DrawName(NameFitParAll[ip].Data(),x,y);    
      double mean = hFitParAll[ip]->GetMean();
      double rms = hFitParAll[ip]->GetRMS();
      char *namemean = "Mean = %f";
      sprintf(name,namemean,mean);
      DrawName(name,x,y*0.8);
      char *namerms = "RMS = %f";
      sprintf(name,namerms,rms);
      DrawName(name,x,y*0.7);
     
      if(SaveOscFit) save_canvas(cFitParAll[ip],canvname,tag);
    }
    else if (best_fit_params->IsNuisance(p) && DrawNuisFit){ // Nuisances

      //if(p != conv::kPhypDm12sq && p != kPhypMixing12) continue;
      //if(p == conv::kRHCNuisSuperkEnergyScale) cout << "kRHCNuisSuperkEnergyScale" << endl;
      //if(p != conv::kNuisSuperkEnergyScale) continue;
      
      string fitname = "canv_fit_";
      fitname.append(best_fit_params->Name(p));
      sprintf(canvname,fitname.c_str());
      cFitParAll[ip] = new TCanvas(canvname,canvname);
      cFitParAll[ip]->cd();
      
      hFitParAll[ip]->Draw();
      hFitParAll[ip]->SetLineColor(kBlack);
      double x = best_fit_params->MinValue(p);
      double y = hFitParAll[ip]->GetMaximum()*0.75;
      //DrawName(NameFitParAll[ip].Data(),x,y);    
      double mean = hFitParAll[ip]->GetMean();
      double rms = hFitParAll[ip]->GetRMS();
      char *namemean = "Mean = %f";
      sprintf(name,namemean,mean);
      //DrawName(name,x,y*0.8);
      char *namerms = "RMS = %f";
      sprintf(name,namerms,rms);
      //DrawName(name,x,y*0.7);
      
      if(DrawNuisExtr){
	hExtrParAll[ip]->Draw("same");
	hExtrParAll[ip]->SetLineColor(kRed);
	double x = best_fit_params->MinValue(p);
	double y = hExtrParAll[ip]->GetMaximum()*0.75;
	//DrawName(NameFitParAll[ip].Data(),x,y);    
	double mean = hExtrParAll[ip]->GetMean();
	double rms = hExtrParAll[ip]->GetRMS();
	char *namemean = "Mean = %f";
	sprintf(name,namemean,mean);
	//DrawName(name,x,y*0.8);
	char *namerms = "RMS = %f";
	sprintf(name,namerms,rms);
	//DrawName(name,x,y*0.7);
      }
       
      if(SaveNuisFit) save_canvas(cFitParAll[ip],canvname,tag);
    }
    
    // Fit error parameters
    
    if(!best_fit_params->IsNuisance(p) && !best_fit_params->IsFixed(p) && DrawOscErr ){ // Oscillations    
      string errname = "canv_err_";
      errname.append(best_fit_params->Name(p));
      sprintf(canvname,errname.c_str());
      cErrParAll[ip] = new TCanvas(canvname,canvname);
      cErrParAll[ip]->cd();

      hErrParAll[ip]->Draw();  
      hErrParAll[ip]->SetLineColor(kBlack);   
      //double x = best_fit_params->MinValue(p);
      double x = 0;
      double y = hErrParAll[ip]->GetMaximum()*0.75;
      //DrawName(NameFitParAll[ip].Data(),x,y);    
      double mean = hErrParAll[ip]->GetMean();
      double rms = hErrParAll[ip]->GetRMS();
      char *namemean = "Mean = %f";
      sprintf(name,namemean,mean);
      //DrawName(name,x,y*0.8);
      char *namerms = "RMS = %f";
      sprintf(name,namerms,rms);
      //DrawName(name,x,y*0.7);

      // //

      // string errname = "canv_err_minos";
      // errname.append(best_fit_params->Name(p));
      // sprintf(canvname,errname.c_str());
      // cErrMinosParAll[ip] = new TCanvas(canvname,canvname);
      // cErrMinosParAll[ip]->cd();

      // hErrPlusParAll[ip]->Draw();  
      // hErrPlusParAll[ip]->SetLineColor(kBlack);   
      // //double x = best_fit_params->MinValue(p);
      // double x = 0;
      // double y = hErrPlusParAll[ip]->GetMaximum()*0.75;
      // //DrawName(NameFitParAll[ip].Data(),x,y);    
      // double mean = hErrPlusParAll[ip]->GetMean();
      // double rms = hErrPlusParAll[ip]->GetRMS();
      // char *namemean = "Mean = %f";
      // sprintf(name,namemean,mean);
      // //DrawName(name,x,y*0.8);
      // char *namerms = "RMS = %f";
      // sprintf(name,namerms,rms);
      // //DrawName(name,x,y*0.7);
      
  
     //if(SaveOscFit) save_canvas(cErrParAll[ip],canvname,tag);
    }
    else if (best_fit_params->IsNuisance(p) && !best_fit_params->IsFixed(p) && DrawNuisErr){ // Nuisances


      //if(p != conv::kNuisSuperkEnergyScale) continue;

      TLine *lcenterr = new TLine(centerr,0,centerr,hErrParAll[ip]->GetMaximum()*1.2);
      lcenterr->SetLineColor(kRed);
      //cout << centerr << endl;
      
      // string canverrname = "canv_err_";
      // canverrname.append(best_fit_params->Name(p));
      // sprintf(canvname,canverrname.c_str());
      // cErrParAll[ip] = new TCanvas(canvname,canvname);
      // cErrParAll[ip]->cd();

      // hErrParAll[ip]->Draw();
      // hErrParAll[ip]->SetLineColor(kBlack);
      // double x = best_fit_params->MinValue(p);
      // double y = hErrParAll[ip]->GetMaximum()*0.75;
      // //DrawName(NameFitParAll[ip].Data(),x,y);    
      // double mean = hErrParAll[ip]->GetMean();
      // double rms = hErrParAll[ip]->GetRMS();
      // char *namemean = "Mean = %f";
      // sprintf(name,namemean,mean);
      // //DrawName(name,x,y*0.8);
      // char *namerms = "RMS = %f";
      // sprintf(name,namerms,rms);
      // //DrawName(name,x,y*0.7);
      //lcenterr->Draw();
      // centerr = get_cent_err(p,skexp);
       // //if(SaveNuisFit) save_canvas(cFitParAll[ip],canvname,tag);

      //
      
      string errname = "canv_err_minos";
      errname.append(best_fit_params->Name(p));
      sprintf(canvname,errname.c_str());
      cErrMinosParAll[ip] = new TCanvas(canvname,canvname);
      cErrMinosParAll[ip]->cd();

      hErrParAll[ip]->Draw();  
      hErrParAll[ip]->SetLineColor(kBlue);   
      //hErrPlusParAll[ip]->Draw();  
      //hErrPlusParAll[ip]->SetLineColor(kBlack);   
      //hErrMinusParAll[ip]->Draw("same");  
      //hErrMinusParAll[ip]->SetLineColor(kBlue);   

      //double x = best_fit_params->MinValue(p);
      //double x = 0;
      //double y = hErrPlusParAll[ip]->GetMaximum()*0.75;
      //DrawName(NameFitParAll[ip].Data(),x,y);    
      //double mean = hErrPlusParAll[ip]->GetMean();
      //double rms = hErrPlusParAll[ip]->GetRMS();
      //char *namemean = "Mean = %f";
      //sprintf(name,namemean,mean);
      //DrawName(name,x,y*0.8);
      //char *namerms = "RMS = %f";
      //sprintf(name,namerms,rms);
      //DrawName(name,x,y*0.7);

      lcenterr->Draw();
    }

    
    
    // Pulls

    if(!best_fit_params->IsNuisance(p) && !best_fit_params->IsFixed(p) && DrawOscPull){     
      string pullname = "canv_pull_";
      pullname.append(best_fit_params->Name(p));
      sprintf(canvname,pullname.c_str());
      cPullAll[ip] = new TCanvas(canvname,canvname);
      cPullAll[ip]->cd();
      hPullAll[ip]->Draw();

      double x = -5.;
      double y = hPullAll[ip]->GetMaximum();
      DrawName(NameFitParAll[ip].Data(),x,y);
      double mean = hPullAll[ip]->GetMean();
      double rms = hPullAll[ip]->GetRMS();
      char *namemean = "Mean = %f";
      sprintf(name,namemean,mean);
      DrawName(name,x,y*0.8);
      char *namerms = "RMS = %f";
      sprintf(name,namerms,rms);
      DrawName(name,x,y*0.6);

      if(SaveOscPull) save_canvas(cPullAll[ip],canvname,tag);
    }
    else if (best_fit_params->IsNuisance(p) && DrawNuisPull){

      //if(p != conv::kPhypDm12sq && p != conv::kPhypMixing12) continue;
      //if(!util::IsBANFF(p)) continue;
      //if(p != conv::kNuisSuperkEnergyScale) continue;

      string pullname = "canv_pull_";
      pullname.append(best_fit_params->Name(p));
      sprintf(name,pullname.c_str());
      cPullAll[ip] = new TCanvas(name,name);
      cPullAll[ip]->cd();
      hPullAll[ip]->Draw();
      double x = -5.;
      double y = hPullAll[ip]->GetMaximum();
      //DrawName(name,x,y);
      double mean = hPullAll[ip]->GetMean();
      double rms = hPullAll[ip]->GetRMS();
      char *namemean = "Mean = %f";
      sprintf(name,namemean,mean);
      //DrawName(name,x,y*0.8);
      char *namerms = "RMS = %f";
      sprintf(name,namerms,rms);
      //DrawName(name,x,y*0.6);

      //if(p!=conv::kNuisBANFF00) return;
    }
    
  } // loop over parameters

  // Draw pull means and widths

  //style
  hPullMean->   GetYaxis()->SetRangeUser(-0.5,1.5);
  hPullMean->   GetXaxis()->LabelsOption("v");
  TLine * lpull = new TLine(0,0,hPullMean->GetXaxis()->GetLast(),0);
  TLine * lpull2 = new TLine(0,1,hPullMean->GetXaxis()->GetLast(),1);
  //color
  hPullWidth->SetMarkerColor(kRed);
  hPullWidth->SetLineColor(kRed);
  //title
  hPullMean-> GetYaxis()->SetTitle("Pull = (bestfit - input) / #sigma_{bf}");
  //draw
  TCanvas * cpull = new TCanvas();
  cpull->SetGridx(1);
  cpull->cd();
  hPullMean->Draw("P");
  hPullWidth->Draw("Psame");
  lpull->Draw();
  lpull2->Draw();
  if(SaveNuisPull) save_canvas(cpull,"pulls",tag);
  cpull->cd();

  //count number of parameters with |pull| > 0.1
  int nbigpull = 0;
  for(int i = 1; i <= hPullMean->GetNbinsX(); i++) {
    if(TMath::Abs(hPullMean->GetBinContent(i)) > 0.1) {
      cout << "Param " << i << " " << hPullMean->GetXaxis()->GetBinLabel(i) << " has a large pull: " << hPullMean->GetBinContent(i) << ":" << endl;
      nbigpull++;
    }
  }
  cout << endl << nbigpull << " parameters have |pull| > 0.1" << endl << endl;
  
  //////////////////////////////////////////////////////////////////////
  // Make disappearance style pull plot
  //PULL PLOT
  if (SaveCombinedPullPlots) {
    FitParamList * fparm = best_fit_params;

    //count free parameters
    int np = fparm->GetNumberParameters();
    int npfree = 0;
    for(int ip = 0; ip < np; ip++) {
      conv::PhysParamId_t p = fparm->WhichParam(ip);
      if(fparm->IsNuisance(p))
        npfree++;
    }
    cout << npfree << " free parameters (" << np << " total)" << endl;
    //make the histogram
    TH1D * pulls    = new TH1D("","",npfree,0,npfree);
    TH1D * values   = new TH1D("","",npfree,0,npfree);
    TH1D * errratio = new TH1D("","",npfree,0,npfree); // postfit/prefit error
    //get deviation from best fit (in sigma)
    for(int ip = 0, ipfree = 1; ip < np; ip++) {
      conv::PhysParamId_t p = fparm->WhichParam(ip);
      if(!fparm->IsNuisance(p)) continue;
      TString label;
      if(ipfree%2)
        label = Form("%s           ",util::AsStringPaperRooTex(p).c_str());
      else
        label = Form("%s",util::AsStringPaperRooTex(p).c_str());
      pulls->   GetXaxis()->SetBinLabel(ipfree, label.Data());
      values->  GetXaxis()->SetBinLabel(ipfree, label.Data());
      errratio->GetXaxis()->SetBinLabel(ipfree, label.Data());
      //get best fit info
      double best = fparm->BestFitValue(p);
      double error = 0;
      if(fparm->HesseErrorsSet())
        error = fparm->HesseError(p);
      else
        error = fparm->MigradError(p);

      //osc params
      if(p == conv::kPhypMixing13) {
        if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsq2Theta) {
          pulls-> SetBinContent(ipfree, (best - conv::kGlobFitSinsq2Theta13) / error);
          values->SetBinContent(ipfree, (best - conv::kGlobFitSinsq2Theta13) / (conv::kGlobFitErrSinsq2Theta13 * conv::kGlobFitSinsq2Theta13));
          errratio->SetBinContent(ipfree, error / (conv::kGlobFitErrSinsq2Theta13 * conv::kGlobFitSinsq2Theta13));
        }
        else if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsqTheta) {
          pulls-> SetBinContent(ipfree, (best - conv::kGlobFitSinsqTheta13) / error);
          values->SetBinContent(ipfree, (best - conv::kGlobFitSinsqTheta13) / (conv::kGlobFitErrSinsqTheta13 * conv::kGlobFitSinsqTheta13));
          errratio->SetBinContent(ipfree, error / (conv::kGlobFitErrSinsqTheta13 * conv::kGlobFitSinsqTheta13));
          cout << util::AsString(p) << "        " << error << "        " << conv::kGlobFitErrSinsqTheta13 << "        " << error / conv::kGlobFitErrSinsqTheta13 << endl;
        }
      }
      else if(p == conv::kPhypMixing12) {
        if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsq2Theta) {
          pulls-> SetBinContent(ipfree, (best - conv::kGlobFitSinsq2Theta12) / error);
          values->SetBinContent(ipfree, (best - conv::kGlobFitSinsq2Theta12) / (conv::kGlobFitErrSinsq2Theta12 * conv::kGlobFitSinsq2Theta12));
          errratio->SetBinContent(ipfree, error / (conv::kGlobFitErrSinsq2Theta12 * conv::kGlobFitSinsq2Theta12));
        }
        else if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsqTheta) {
          pulls-> SetBinContent(ipfree, (best - conv::kGlobFitSinsqTheta12) / error);
          values->SetBinContent(ipfree, (best - conv::kGlobFitSinsqTheta12) / (conv::kGlobFitErrSinsqTheta12 * conv::kGlobFitSinsqTheta12));
          errratio->SetBinContent(ipfree, error / (conv::kGlobFitErrSinsqTheta12 * conv::kGlobFitSinsqTheta12));
          cout << util::AsString(p) << "        " << error << "        " << conv::kGlobFitErrSinsqTheta12 << "        " << error / conv::kGlobFitErrSinsqTheta12 << endl;
        }
      }
      else if(p == conv::kPhypCPV    ) {
        pulls-> SetBinContent(ipfree, (best - conv::kGlobFitDeltaCP) / error);
        values->SetBinContent(ipfree, (best - conv::kGlobFitDeltaCP) / conv::kpi);//conv::kGlobFitErrDeltaCP);
        errratio->SetBinContent(ipfree, error / conv::kpi);
        cout << best << endl;
        cout << error << endl;
        cout << 2*conv::kpi << endl;
      }
      else if(p == conv::kPhypDm12sq) {
        pulls-> SetBinContent(ipfree, (best - conv::kGlobFitDm12sq) / error);
        values->SetBinContent(ipfree, (best - conv::kGlobFitDm12sq) / (conv::kGlobFitErrDm12sq * conv::kGlobFitDm12sq));
        errratio->SetBinContent(ipfree, error / (conv::kGlobFitErrDm12sq * conv::kGlobFitDm12sq));
          cout << util::AsString(p) << "        " << error << "        " << conv::kGlobFitErrDm12sq << "        " << error / conv::kGlobFitErrDm12sq << endl;
      }
      else if(util::IsBANFF(p)) {
        pulls-> SetBinContent(ipfree, (best - skexp->GetBANFF()->BestFitValue(p)) / error);
        values->SetBinContent(ipfree, (best - skexp->GetBANFF()->BestFitValue(p)) / skexp->GetBANFF()->OneSigmaError(p));
        errratio->SetBinContent(ipfree, error / skexp->GetBANFF()->OneSigmaError(p));
      }
      else if(util::IsXSec(p)) {
        //"nominal best fit" is 0
        //"1 sigma error" is 1
        pulls-> SetBinContent(ipfree, (best - 0) / error);
        values->SetBinContent(ipfree, (best - 0) / 1);
        errratio->SetBinContent(ipfree, error / 1);
      }
      else if(p == conv::kNuisSuperkEnergyScale || p == conv::kMultiRingNuisSuperkEnergyScale) {
        //"nominal best fit" is 0
        //"1 sigma error" is 1
        pulls-> SetBinContent(ipfree, (best - 0) / error);
        values->SetBinContent(ipfree, (best - 0) / 1);
        errratio->SetBinContent(ipfree, error / 1);
      }
      else if(util::IsSKDetFSI(p)) {
        //"nominal best fit" is 1
        pulls-> SetBinContent(ipfree, (best - 1) / error);
	if(util::IsMultiRingSystematic(p)){
	   values->SetBinContent(ipfree, (best - 1) / skexp->GetSuperkDetectorFSI()->MultiRingOneSigmaError(p));
	   errratio->SetBinContent(ipfree, error / skexp->GetSuperkDetectorFSI()->MultiRingOneSigmaError(p));
	}
	else{
	  values->SetBinContent(ipfree, (best - 1) / skexp->GetSuperkDetectorFSI()->OneSigmaError(p));
	  errratio->SetBinContent(ipfree, error / skexp->GetSuperkDetectorFSI()->OneSigmaError(p));
	}
      }
      else {
        cerr << "Unknown pull for parameter " << util::AsStringCompact(p) << endl;
      }
      cout << util::AsString(p) << "        " << errratio->GetBinContent(ipfree) << endl;
      ipfree++;
    }//ip
    //style
    pulls->   GetYaxis()->SetRangeUser(-1,1);
    pulls->   GetXaxis()->LabelsOption("v");
    values->  GetYaxis()->SetRangeUser(-1,1);
    values->  GetXaxis()->LabelsOption("v");
    errratio->GetYaxis()->SetRangeUser( 0,2);
    errratio->GetXaxis()->LabelsOption("v");
    TLine * lpull = new TLine(0,0,pulls->GetXaxis()->GetLast(),0);
    TLine * lerr  = new TLine(0,1,pulls->GetXaxis()->GetLast(),1);
    //title
    pulls-> GetYaxis()->SetTitle("Pull = (bestfit - input) / #sigma_{bf}");
    values->GetYaxis()->SetTitle("Pull = (bestfit - input) / #sigma_{input}");
    errratio->GetYaxis()->SetTitle("#sigma_{bf} / #sigma_{input}");
    
    // Make plot of pull with unfitted params:
    TH1D * pulls_unfit = new TH1D(*pulls);
    pulls_unfit->Multiply(errratio);
    pulls_unfit->SetMarkerStyle(kOpenSquare);
    pulls_unfit->SetMarkerColor(kRed);
    
    //draw
    TCanvas * cpull = new TCanvas();
    cpull->SetGridx(1);
    cpull->SetBottomMargin(0.4);
    cpull->cd();
    pulls->Draw("P");
    lpull->Draw();
    pulls_unfit->Draw("same P");
    save_canvas(cpull,"combined_pulls",tag);
    /*
    cpull->cd();
    values->Draw("P");
    lpull->Draw();
    save_canvas(cpull,"values",tag);
    cpull->cd();
    errratio->Draw("P");
    lerr->Draw();
    save_canvas(cpull,"err_ratio",tag);*/
  }

  // Draw value pull means and widths

  //style
  hPullValMean->   GetYaxis()->SetRangeUser(-0.5,1.5);
  hPullValMean->   GetXaxis()->LabelsOption("v");
  TLine * lpullval = new TLine(0,0,hPullValMean->GetXaxis()->GetLast(),0);
  //color
  hPullValWidth->SetMarkerColor(kRed);
  hPullValWidth->SetLineColor(kRed);
  //title
  hPullValMean-> GetYaxis()->SetTitle("Pull = (bestfit - input) / #sigma_{input}");
  //draw
  TCanvas * cpullval = new TCanvas();
  cpullval->SetGridx(1);
  cpullval->cd();
  hPullValMean->Draw("P");
  hPullValWidth->Draw("Psame");
  lpullval->Draw();
  if(SaveNuisPull) save_canvas(cpullval,"pullsval",tag);
  cpullval->cd();

  //count number of parameters with |pullval| > 0.1

  // for(int i = 1; i <= hPullValMean->GetNbinsX(); i++) {
  //   if(TMath::Abs(hPullValMean->GetBinContent(i)) > 0.1) {
  //     cout << hPullValMean->GetXaxis()->GetBinLabel(i) << " has a large value: " << hPullValMean->GetBinContent(i) << endl;
  //     nbigpull++;
  //   }
  // }



  return;
}
  
