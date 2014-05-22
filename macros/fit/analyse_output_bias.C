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

#define DrawNuisFit false
#define SaveNuisFit false

#define DrawNuisPull true
#define SaveNuisPull false

#define DrawNuisRes false
#define SaveNuisRes false

#define DrawNuisExtr false
#define SaveNuisExtr false

#define NBinsFitPar 200
#define NBinsPull 200


// Declare Histograms

TString *NameFitParAll;

TH1F **hExtrParAll;
TH1F **hFitParAll;
TH1F **hPullAll;
TH1F **hPullValAll;
TH1F **hResAll;
TH1F *hPullMean;
TH1F *hPullWidth;
TH1F *hPullValMean;
TH1F *hPullValWidth;

TCanvas **cExtrParAll;
TCanvas **cFitParAll;
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
  //double pull = (fitval-centval)/sqrt(centerr*centerr - fiterr*fiterr); 
  //double pull = (fitval-centval)/centerr; 
  double pull = (fitval-centval)/fiterr; 
  //cout << fitval << " " << centval << " " << centerr << endl;
  return pull;
}

double calc_pull_val(double centval,double centerr,double fitval,double fiterr){ 
  double pull = (fitval-centval)/fiterr; 
  return pull;
}


void analyse_output_bias
(
 const char * infile,
 const char * tag,
 const char * hierarchy = "normal",
 const char * dmsqscheme = "abs23",
 const char * mixparamchoice = "sinsqtheta",
 const char * config_file = "analysis_configuration.xml",
 const char * hypothesis = "3f_osc",   // "3f_osc", "4f_sterile_osc", "5f_sterile_osc 
 )
{  
  using namespace std;
  using namespace valor;
  using namespace valor::conv;

  // Global options
  Opt::Instance()->gMassHierarchyAssumed = util::HierarchyFromString(hierarchy);
  Opt::Instance()->gOscDmsqScheme        = util::DmsqSchemeFromString(dmsqscheme);
  Opt::Instance()->gOscMixingParamChoice = util::MixingTypeFromString(mixparamchoice);
  if(Opt::Instance()->Undefined()) {
    cerr << "Global options not defined correctly. Exiting..." << endl;
    cerr << *Opt::Instance() << endl;
    exit(1);
  }
 
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
  //double fittime1, fittime2;
  //int fitstatus1, fitstatus2;
  //double chisq1, chisq2;

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
  //fitnt->SetBranchAddress("data",                 &data);
  //fitnt->SetBranchAddress("best_fit_pdf",         &best_fit_pdf);
  //fitnt->SetBranchAddress("best_fit_pdf_fit1",    &best_fit_pdf_fit1);
  //fitnt->SetBranchAddress("best_fit_pdf_fit2",    &best_fit_pdf_fit2);

  // Initialize Histograms

  cout << endl;
  cout << "Initialize Histograms" << endl;
  cout << endl;
  
  fitnt->GetEntry(0); // Get Tree general informations
  int np = best_fit_params->GetNumberParameters();
  int nph = input_params->NParams();

  // cout << "np = " << np << endl;
  // cout << "nph = " << nph << endl;




  NameFitParAll = new TString[np];
  cExtrParAll = new TCanvas* [np];
  cFitParAll = new TCanvas* [np];
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
    sprintf(name,best_fit_params->Name(p).c_str());
    hExtrParAll[ip] = new TH1F(name,name,NBinsFitPar,min,max);

    // Fitted distributions
    max = best_fit_params->MaxValue(p);
    min = best_fit_params->MinValue(p);
    sprintf(name,best_fit_params->Name(p).c_str());
    hFitParAll[ip] = new TH1F(name,name,NBinsFitPar,min,max);

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
    if(itoy==0) cout <<*input_params << endl;

    if(!(itoy%100)) cout << "toy " << itoy << endl;
    
    for(int ip = 0; ip < np; ip++) {
      conv::PhysParamId_t p = best_fit_params->WhichParam(ip);
   
      double centval=-999.;
      double centerr=-999.;
      double pull=-999.;

      double bestfitval = best_fit_params->BestFitValue(p);
      double bestfiterr = best_fit_params->MigradError(p);
      centval = input_params->Value(p);

      // Extracted parameters

      hExtrParAll[ip]->Fill(centval);
      
      // Fitted parameters

      hFitParAll[ip]->Fill(bestfitval); 
            
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
	
	pull = calc_pull_nuis(centval,centerr,bestfitval,bestfiterr);
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
    
    mean = hPullAll[ip]->GetMean();
    rms = hPullAll[ip]->GetRMS();

    hPullMean->SetBinContent(ip,mean);
    hPullWidth->SetBinContent(ip,rms);
    
    //mean = hPullValAll[ip]->GetMean();
    //rms = hPullValAll[ip]->GetRMS();

    //hPullValMean->SetBinContent(ip,mean);
    //hPullValWidth->SetBinContent(ip,rms);

  }
  
  // Draw histograms
  //int np = best_fit_params->GetNumberParameters();
  for(int ip = 0; ip < np; ip++) {
    conv::PhysParamId_t p = best_fit_params->WhichParam(ip);
    

    // Extracted parameters
    
    if(best_fit_params->IsNuisance(p) && DrawNuisExtr){ // Oscillations    
      sprintf(canvname,best_fit_params->Name(p).c_str());
      cExtrParAll[ip] = new TCanvas(canvname,canvname);
      cExtrParAll[ip]->cd();
      hExtrParAll[ip]->Draw();
      
      double x = best_fit_params->MinValue(p);
      double y = hExtrParAll[ip]->GetMaximum()*0.75;
      DrawName(NameFitParAll[ip].Data(),x,y);    
      double mean = hExtrParAll[ip]->GetMean();
      double rms = hExtrParAll[ip]->GetRMS();
      char *namemean = "Mean = %f";
      sprintf(name,namemean,mean);
      DrawName(name,x,y*0.8);
      char *namerms = "RMS = %f";
      sprintf(name,namerms,rms);
      DrawName(name,x,y*0.7);

      if(SaveNuisExtr) save_canvas(cExtrParAll[ip],canvname,tag);
    }
    
    // Fit parameters
    
    if(!best_fit_params->IsNuisance(p) && !best_fit_params->IsFixed(p) && DrawOscFit){ // Oscillations    
      sprintf(canvname,best_fit_params->Name(p).c_str());
      cFitParAll[ip] = new TCanvas(canvname,canvname);
      cFitParAll[ip]->cd();
      hFitParAll[ip]->Draw();
      
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
      
      sprintf(canvname,best_fit_params->Name(p).c_str());
      cFitParAll[ip] = new TCanvas(canvname,canvname);
      cFitParAll[ip]->cd();
      hFitParAll[ip]->Draw();
      
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
      
      if(SaveNuisFit) save_canvas(cFitParAll[ip],canvname,tag);
    }
    
    // Pulls

    if(!best_fit_params->IsNuisance(p) && !best_fit_params->IsFixed(p) && DrawOscPull){     
      string pullname = "pull_";
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

      string pullname = "pull_";
      pullname.append(best_fit_params->Name(p));
      sprintf(name,pullname.c_str());
      cPullAll[ip] = new TCanvas(name,name);
      cPullAll[ip]->cd();
      hPullAll[ip]->Draw();
      double x = -5.;
      double y = hPullAll[ip]->GetMaximum();
      DrawName(name,x,y);
      double mean = hPullAll[ip]->GetMean();
      double rms = hPullAll[ip]->GetRMS();
      char *namemean = "Mean = %f";
      sprintf(name,namemean,mean);
      DrawName(name,x,y*0.8);
      char *namerms = "RMS = %f";
      sprintf(name,namerms,rms);
      DrawName(name,x,y*0.6);

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
      cout << hPullMean->GetXaxis()->GetBinLabel(i) << " has a large pull: " << hPullMean->GetBinContent(i) << endl;
      nbigpull++;
    }
  }
  cout << endl << nbigpull << " parameters have |pull| > 0.1" << endl << endl;


  // // Draw value pull means and widths

  // //style
  // hPullValMean->   GetYaxis()->SetRangeUser(-0.5,1.5);
  // hPullValMean->   GetXaxis()->LabelsOption("v");
  // TLine * lpullval = new TLine(0,0,hPullValMean->GetXaxis()->GetLast(),0);
  // //color
  // hPullValWidth->SetMarkerColor(kRed);
  // hPullValWidth->SetLineColor(kRed);
  // //title
  // hPullValMean-> GetYaxis()->SetTitle("Pull = (bestfit - input) / #sigma_{input}");
  // //draw
  // TCanvas * cpullval = new TCanvas();
  // cpullval->SetGridx(1);
  // cpullval->cd();
  // hPullValMean->Draw("P");
  // hPullValWidth->Draw("Psame");
  // lpullval->Draw();
  // if(SaveNuisPull) save_canvas(cpullval,"pullsval",tag);
  // cpullval->cd();

  // //count number of parameters with |pullval| > 0.1

  // for(int i = 1; i <= hPullValMean->GetNbinsX(); i++) {
  //   if(TMath::Abs(hPullValMean->GetBinContent(i)) > 0.1) {
  //     cout << hPullValMean->GetXaxis()->GetBinLabel(i) << " has a large value: " << hPullValMean->GetBinContent(i) << endl;
  //     nbigpull++;
  //   }
  // }



  return;
}
  
