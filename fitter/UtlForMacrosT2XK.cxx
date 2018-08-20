//____________________________________________________________________________
/*

Copyright (c) 2010-2014, The VALOR group. All rights reserved.

ChangeLog:

*/
//____________________________________________________________________________


#include "UtlForMacrosT2XK.h"

#include "SuperkExpectation.h"
#include "SuperkObservation.h"
#include "SuperkMC.h"
#include "BANFF.h"
#include "SuperkDetectorFSI.h"
#include "FitMarginalization.h"
#include "Build.h"

#include "valor/Utilities/Env.h"
#include "valor/Conventions/Opt.h"
#include "valor/Conventions/PhysParamId.h"
#include "valor/Conventions/PhysHypothesisId.h"
#include "valor/Conventions/Constants.h"
#include "valor/FitCore/PhysHypothesis.h"
#include "valor/Utilities/UtlXML.h"
#include "valor/Messenger/Messenger.h"

// Deprecated
#include "valor/CmdLine/RunOpt_TO_DELETE.h"

#include "TSystem.h"
#include "TCanvas.h"
#include "TMath.h"

#include <sstream>
#include <fstream>
#include <iomanip>

#include <iostream>
#include <string>

using std::cout;
using std::cerr;
using std::endl;
using std::string;

using namespace valor;
using namespace valor::conv;

//____________________________________________________________________________
void valor::util::PerformSysThrow(unsigned int tweak_systematics, SuperkExpectation * skexp, bool print)
{
  if(tweak_systematics == 0) {
    //make sure systematics are reset
    skexp->ResetSystematicVariations();
    return;
  }
  else if(tweak_systematics == 1) {
    // tweak all sources of systematic err (gaussian err)
    LOG("UtlForMacros",pINFO) << "- Randomizing all systematic params in pdf...";
    skexp->TweakAllSystematics();
    // skexp->TweakSystematic("sk_energy_scale",0.0,true);
  }
  else if(tweak_systematics == 2) {
    // special tweak, increasing numu flux parameters <700MeV by 5%
    LOG("UtlForMacros",pINFO) << "- Increasing numu flux parameters <700MeV by 5% in pdf...";
    TVectorD twkv(skexp->GetBANFF()->InitialValues());
    twkv[0] += 0.05;
    twkv[1] += 0.05;
    twkv[2] += 0.05;
    twkv[3] += 0.05;
    skexp->GetBANFF()->Tweak(true,&twkv);
  }
  else if(tweak_systematics == 3) {
    LOG("UtlForMacros",pINFO) << "- Randomizing all but SK systematic params in pdf...";
    skexp->TweakAllBANFFSystematics();
    skexp->TweakAllUncorrelatedXSecSystematics();
  }
  else if(tweak_systematics == 4) {
    LOG("UtlForMacros",pINFO) << "- Setting xsec_norm_nc1piC to -0.5";
    skexp->TweakSystematic("xsec_norm_nc1piC",-0.5);
  }
  else if(tweak_systematics == 5) {
    LOG("UtlForMacros",pINFO) << "- Increasing flux parameters >1500MeV by 15% in pdf...";
    TVectorD twkv(skexp->GetBANFF()->InitialValues());
    twkv[6]  += 0.15;
    twkv[7]  += 0.15;
    twkv[8]  += 0.15;
    twkv[9]  += 0.15;
    twkv[10] += 0.15;
    twkv[14] += 0.15;
    twkv[15] += 0.15;
    skexp->GetBANFF()->Tweak(true,&twkv);
  }
  else if(tweak_systematics == 6) {
    LOG("UtlForMacros",pINFO) << "- Randomizing all BANFF systematic params in pdf...";
    skexp->TweakAllBANFFSystematics();
  }
  else if(tweak_systematics == 7) {
    LOG("UtlForMacros",pINFO) << "- Randomizing all uncorrelated cross section systematic params in pdf...";
    skexp->TweakAllUncorrelatedXSecSystematics();
  }
  else if(tweak_systematics == 8) {
    LOG("UtlForMacros",pINFO) << "- Randomizing all SK+FSI systematic params in pdf...";
    skexp->TweakAllSKDetectorFSISystematics();
  }
  else {
    LOG("UtlForMacros",pFATAL) << "PerformSysThrow() tweak_systematics = " << tweak_systematics
	 << " not implemented";
    exit(-1);
  }
  if(print)
    skexp->PrintSystematicVariations();
}
//____________________________________________________________________________
#ifndef __USE_2D_TEMPLATES__ //1D templates
void valor::util::AddEventsToXML(ofstream & fxml, ofstream & ftxt, TH1D * h, valor::conv::SuperkSample_t s, string * full_points_to_add, TRandom3 & rand)
#else
void valor::util::AddEventsToXML(ofstream & fxml, ofstream & ftxt, TH2D * h, valor::conv::SuperkSample_t s, string * full_points_to_add, TRandom3 & rand)
#endif
{
  ftxt << util::AsStringCompact(s) << ":" << endl;
  fxml << "        <" << util::GetXMLDatasetName(s) << ">" << endl;
  fxml << "          <" << util::GetXMLDatasetParamsName() << ">" << endl;

  int ievtot = 0;

#ifndef __USE_2D_TEMPLATES__ //1D templates
  const int n_points_per_line_xml = 10; // this should be an option. But 1D templates are going so there's no point...

  // Adding Ereco values as a comma separated list.
  // Assuming uniform distribution within each bin, randomize values within each bin.
  for(int ibin = 1; ibin <= h->GetNbinsX(); ibin++) {
    double Nb   = h->GetBinContent (ibin); // number of events in this Ereco bin
    double Emin = h->GetXaxis()->GetBinLowEdge (ibin); // Ereco bin min
    double Emax = h->GetXaxis()->GetBinUpEdge  (ibin); // Ereco bin max
    //double Ecent = h->GetXaxis()->GetBinCenter (ibin); // Ereco bin center
    for(int iev = 0; iev < Nb; iev++) {
      double E = Emin + rand.Uniform() * (Emax - Emin);
      if(ievtot > 0) {
	fxml << ", "; 
	ftxt << ","; 
	if ( ! (ievtot % n_points_per_line_xml) ) fxml << endl; 
      }
      fxml << E; 
      ftxt << E; 
      ievtot++;
    }//iev
  }//ibin
  // add extra data-points, if any
  if(full_points_to_add) {
    fxml << ",\n" << endl; // should really have an XML comment here but TXMLEngine won't let us.
    fxml << *full_points_to_add;
    ftxt << "," << *full_points_to_add;
  }
  fxml << endl;

#else //2D templates
       
  // Adding 2D template values with each event as its own node containing 2 values
  // Assuming uniform distribution within each bin, randomize values within each bin.
  for(int ibinx = 1; ibinx <= h->GetNbinsX(); ibinx++) {
    double Xmin  = h->GetXaxis()->GetBinLowEdge(ibinx); // X bin min
    double Xmax  = h->GetXaxis()->GetBinUpEdge (ibinx); // X bin max
    //double Xcent = h->GetXaxis()->GetBinCenter (ibinx); // X bin center
    for(int ibiny = 1; ibiny <= h->GetNbinsY(); ibiny++) {
      double Nb    = h->GetBinContent (ibinx, ibiny); // number of events in this bin
      double Ymin  = h->GetYaxis()->GetBinLowEdge(ibiny); // Y bin min
      double Ymax  = h->GetYaxis()->GetBinUpEdge (ibiny); // Y bin max
      //double Ycent = h->GetYaxis()->GetBinCenter (ibiny); // Y bin center
      for(int iev = 0; iev < Nb; iev++) {
	double X = Xmin + rand.Uniform() * (Xmax - Xmin);
	double Y = Ymin + rand.Uniform() * (Ymax - Ymin);
	if(ievtot > 0) {
	  ftxt << ": "; 
	}
	fxml << "            <ev" << ievtot << ">"
	     << X << "," << Y
	     << "</ev" << ievtot << ">"
	     << endl;
	ftxt << X << "," << Y; 
	ievtot++;
      }//iev
    }//ibiny
  }//ibinx
  // add extra data-points, if any
  if(full_points_to_add) {
    fxml << endl; // should really have an XML comment here but TXMLEngine won't let us.
    vector<string> vecpairs = util::Tokens(*full_points_to_add,",");
    assert(vecpairs.size()>0);
    for(unsigned int ipair = 0; ipair < vecpairs.size(); ipair++) {
      vector<string> values = util::Tokens(vecpairs[ipair], ":");
      assert(values.size() == 2);
      fxml << "            <ev" << ievtot << ">"
	   << atof(values[0].c_str()) << "," << atof(values[1].c_str())
	   << "</ev" << ievtot << ">"
	   << endl;
      ievtot++;
    }//ipair
    ftxt << ":" << *full_points_to_add;
  }
  // add the nevents to the XML
  fxml << "            <nevents> " << ievtot << "</nevents>" << endl;
#endif //__USE_2D_TEMPLATES__
	 
  ftxt << endl;
  fxml << "          </" << util::GetXMLDatasetParamsName() << "> " << endl;
  fxml << "        </" << util::GetXMLDatasetName(s) << ">    " << endl;       
}
//____________________________________________________________________________
#ifndef __USE_2D_TEMPLATES__ //1D templates
void valor::util::AddExtraEventsToHist(TH1D * h, string * full_points_to_add)
#else
void valor::util::AddExtraEventsToHist(TH2D * h, string * full_points_to_add)
#endif
{
  if(full_points_to_add) {
#ifndef __USE_2D_TEMPLATES__ //1D templates
    vector<string> vecp = util::Tokens(*full_points_to_add, ",");
    assert(vecp.size() > 0);
    for(unsigned int ip=0; ip < vecp.size(); ip++) {
      h->Fill(atof(vecp[ip].c_str()));
    }//ip
#else //2D templates
    vector<string> vecpairs = util::Tokens(*full_points_to_add,",");
    assert(vecpairs.size()>0);
    for(unsigned int ipair = 0; ipair < vecpairs.size(); ipair++) {
      vector<string> values = util::Tokens(vecpairs[ipair], ":");
      assert(values.size() == 2);
      h->Fill(atof(values[0].c_str()), atof(values[1].c_str()));
    }//ipair
#endif
  }//full_points_to_add
}
//____________________________________________________________________________
string * valor::util::ConstructPointsToAddFromXML(TXMLEngine & xml, XMLNodePointer_t & root_node, valor::conv::SuperkSample_t s, string sextra, int datasetid)
{
  if(util::ArrIdx::Instance()->SuperkSampleIdx(s) == -1)
    return 0;
#ifndef __USE_2D_TEMPLATES__ //1D templates
  string sxml = util::XML2String(&xml, root_node, Form("superk/observation/data_set_%d/%s/ereco_values", datasetid, util::GetXMLMCTreeName(s).c_str()));
  string * sout = new string(sxml);
  if(!sextra.empty()) {
    sout->append(",");
    sout->append(sextra);
  }
#else
  int nevents = util::XML2Int(&xml, root_node, Form("superk/observation/data_set_%d/%s/%s/nevents",
						    datasetid, util::GetXMLMCTreeName(s).c_str(), util::GetXMLDatasetParamsName().c_str()));
  string * sout = new string();
  for(int iev = 0; iev < nevents; iev++) {
    if(iev != 0)
      sout->append(", ");
    string sxml = util::XML2String(&xml, root_node, Form("superk/observation/data_set_%d/%s/%s/ev%d",
							 datasetid, util::GetXMLMCTreeName(s).c_str(), util::GetXMLDatasetParamsName().c_str(), iev));
    std::replace(sxml.begin(), sxml.end(), ',', ':');
    sout->append(sxml);
  }//iev
  if(!sextra.empty()) {
    sout->append(", ");
    sout->append(sextra);
  }
#endif
  return sout;
}
//____________________________________________________________________________
void valor::util::PerformMargFit(int argc, char ** argv,RunOpt *runopt,FitMarginalization *fitmarg, SuperkExpectation *skexp, SuperkObservation *skobs_fakedata, bool /*print*/,string masshier) // NEW FC
{    
  //
  // If mass hierarchy is not set (masshier==""), use the global mass hierarchy.
  // Needed for when Feldman-Cousins is not run
  //

  
  Env::Instance()->AutoSet("t2xk");
  cout << *Env::Instance() << endl;
  //gROOT->SetBatch(true);
  

  //RunOpt *runopt = new RunOpt(argc, argv, kParseCommon | kParseToys | kParseGlobal | kParseFitMarg | kParse3FOscParams);
  //RunOpt runopt(argc, argv, kParseCommon | kParseToys | kParseGlobal | kParseFitMarg | kParse3FOscParams);
  
  //RunOpt::ReadFromCommandLine(argc, argv, parseopt);

  // Set the second parameter (1D fit)
  //const char * freeparname1        = "";  // free parameter
  //const double freeparmin1         = 0.; // minimum value of free parameter
  //const double freeparmax1         = 0.; // maxinum value of free parameter
  //const int    NptScan1            = 1; //201, //21,      // # pts of free parameter to scan
  //

  // // Reset all the variables
  int nfits=-999;
  //int nchisqcalls=-999;
  int niterations=-999;
  int fitstatus=-999;
  double chisqmin=conv::kALargeNum;
  double bestfit0=-999;
  double bestfit1=-999;
  
  // Set free param, ntoys, # pt to scan, min par, max par (if min and max not passed, use default ones for that par)

  // free parameter 

  //PhysParamId_t freepar0 = conv::kPhypUndefined;
  string sfreeparname0(runopt->Freeparname0()); 
  PhysParamId_t freepar1 = conv::kPhypUndefined;
  string sfreeparname1(runopt->Freeparname1());  
     
  PhysParamId_t parfree[2];
  parfree[0] = GetFreePar(sfreeparname0);
  
  if(sfreeparname1=="")  
    parfree[1] = freepar1; // 1D fit
  else 
    parfree[1] = GetFreePar(sfreeparname1); // 2D fit
 
  // # of toys for marginalization
  int ntoys = runopt->NtoysMarg();
  
  // # of free parameter pts to scan 
  int nbinspar[2];
  nbinspar[0] = runopt->NptScan0();
  nbinspar[1] = runopt->NptScan1();
  
  // free parameter range
  double minpar[2];
  double maxpar[2];
  minpar[0] = runopt->Freeparmin0();
  maxpar[0] = runopt->Freeparmax0();
  minpar[1] = runopt->Freeparmin1();
  maxpar[1] = runopt->Freeparmax1();
  
  // free parameter sub-range
  const int    first_bin_freepar0   = -999; // first freepar of the sub-range
  const int    last_bin_freepar0    = -999; // last freepar of the sub-range
  const int    first_bin_freepar1   = -999; // first freepar of the sub-range
  const int    last_bin_freepar1    = -999; // last freepar of the sub-range
  
  int subminbin[2];
  int submaxbin[2];
  subminbin[0] = first_bin_freepar0; submaxbin[0] = last_bin_freepar0;
  subminbin[1] = first_bin_freepar1; submaxbin[1] = last_bin_freepar1;
  
  // Implement the constructor giving the following input parameters:
  //FitMarginalization *fitmarg = new FitMarginalization();   
  fitmarg->SetPhysParFree(parfree,ntoys,nbinspar,minpar,maxpar);
  fitmarg->SetSubRange(subminbin,submaxbin); // By default it scans all the range
  fitmarg->SetTestStatistic(runopt->Teststat()); // Set Test Statistic ("rate-only" or "rate+shape")
  fitmarg->Init(skexp); // Before need SetPhysParFree() to initialize the arrays
  fitmarg->SetStoreLikeWeights(runopt->StoreLikeWeights()); // Set if to store full likelihood for each toy (after fitmarg->Init())
  
  

  // Set Mass Hierarchy marginalization
  Opt::Instance()->gMassHierarchyAssumed = util::HierarchyFromString(runopt->Hierarchy_marg()); // needed to build the right binary file name 
  
  if(masshier==""){ // use global mass hierarchy
    Opt::Instance()->gMassHierarchyAssumed = util::HierarchyFromString(runopt->Hierarchy_marg()); // needed to build the right binary file name 
  }
  else if(masshier=="normal"){
    Opt::Instance()->gMassHierarchyAssumed = util::HierarchyFromString("normal"); // needed to build the right binary file name 
  }
  else if(masshier=="inverted"){
    Opt::Instance()->gMassHierarchyAssumed = util::HierarchyFromString("inverted"); // needed to build the right binary file name 
  }
  else{
    cerr << endl;
    cerr << "You didn't set the marginalization file mode correctly" << endl;
    cerr << "The following options can be chosen:" << endl;
    cerr << " - No mode" << endl;
    cerr << " - normal" << endl;
    cerr << " - inverted" << endl;
    cerr << endl;
    exit(1);
  }

  //
  // MH MUST BE THE SAME AS FOR TOYS!!! TO CHECK  
  //
  fitmarg->SetMargMassHierarchy(runopt->DoMargMH()); 
  //

  // Set fit parameter list
  fitmarg->SetFitParamList(runopt->Fit_param_list().c_str());
  fitmarg->PrintPhysHypothesis();
  
  // Build the map fileID - sub-range
  
  fitmarg->SetMaxIterationSingleFit(2);  // Max # of iterations per fit
  fitmarg->SetMaxFitNumber(3); // Max # of fits max # of iterations
  

  ////// NEW FC
  if(masshier==""){
    fitmarg->SetTagFileName(runopt->Infilename().c_str());
    cout << "No Mass Hierarchy" << endl;
    cout << runopt->Infilename().c_str() << endl;
    cout << endl;
  }
  else if(masshier=="normal"){
    fitmarg->SetTagFileName(runopt->Infilename_nh().c_str());
    cout << "Normal Hierarchy" << endl;
    cout << runopt->Infilename_nh().c_str() << endl;
    cout << endl;
  }
  else if(masshier=="inverted"){
    fitmarg->SetTagFileName(runopt->Infilename_ih().c_str());
    cout << "Inverted Hierarchy" << endl;
    cout << runopt->Infilename_ih().c_str() << endl;
    cout << endl;
  }
  else{
    cerr << endl;
    cerr << "You didn't set the marginalization file mode correctly" << endl;
    cerr << "The following options can be chosen:" << endl;
    cerr << " - No mode" << endl;
    cerr << " - normal" << endl;
    cerr << " - inverted" << endl;
    cerr << endl;
    exit(1);
  }
  /////////

  fitmarg->SetFileID(runopt->Nfiles());



  if(runopt->DoMinimize()) {
    
    // 1D deltaCP fit
    if(runopt->Freeparname0()=="deltacp" && runopt->Freeparname1()==""){
      util::MargFit_1D_deltaCP(runopt,fitmarg,skobs_fakedata,chisqmin,bestfit0,bestfit1,nfits,fitstatus,niterations);
    }
    // 1D betaBar fit (nuebar appearance)
    else if(runopt->Freeparname0()=="betaBar" && runopt->Freeparname1()=="" && runopt->Teststat()=="rate-only"){
      util::MargFit_1D_betaBarRateOnly(runopt,fitmarg,skobs_fakedata,chisqmin,bestfit0,bestfit1,nfits,fitstatus,niterations);
    }
    // 1D betaBar fit (nuebar appearance rate+shape)
    else if(runopt->Freeparname0()=="betaBar" && runopt->Freeparname1()==""){
      util::MargFit_1D_betaBar(runopt,fitmarg,skobs_fakedata,chisqmin,bestfit0,bestfit1,nfits,fitstatus,niterations);
    }
    // 2D sinsqth13 Vs dcp fit
    else if(runopt->Freeparname0()=="mixing13" && runopt->Freeparname1()=="deltacp"){
      
      bool penalty_par0       = runopt->PenaltyTerm0();
      bool penalty_par1       = runopt->PenaltyTerm1();
      bool constrain          = penalty_par0 && penalty_par1;
      double sinsq2th13_mean  = conv::kPrior_Sinsq2th13_Used_Mean; // Used (sinsq2th13)
      double sinsq2th13_sigma = conv::kPrior_Sinsq2th13_Used_RMS; // Used (sinsq2th13)
      
      if(penalty_par0 && !penalty_par1){
	cout << endl;
	cout << "Use penalty term only for free parameter 0:" << endl;
	cout << " sinsq2th13_mean = "  << sinsq2th13_mean  << endl;
	cout << " sinsq2th13_sigma = " << sinsq2th13_sigma << endl;
	cout << endl;
      }
      else if(!penalty_par0 && !penalty_par1){
	cout << endl;
	cout << "No penalty term is used" << endl;
	cout << endl;
      }
      else{
	cerr << endl;
	cerr << "valor::util::PerformMargFit()" << endl;
	cerr << "Penalty term for parameter 1 is not implemented yet" << endl;
	cerr << endl;
	exit(1);
      }
      
      fitmarg->UsePenaltyTerm(penalty_par0,penalty_par1);
      fitmarg->SetPenaltyTerm(sinsq2th13_mean, sinsq2th13_sigma); // to be used after UsePenaltyTerm()
      
      util::MargFit_2D_sinsqth13_deltaCP(runopt,fitmarg,skobs_fakedata,chisqmin,bestfit0,bestfit1,nfits,fitstatus,niterations);
    }
    
    // 2D sinsqth13 Vs dcp fit
    else if(runopt->Freeparname0()=="deltacp" && runopt->Freeparname1()=="mixing13"){
      
      bool penalty_par0       = runopt->PenaltyTerm0();
      bool penalty_par1       = runopt->PenaltyTerm1();
      bool constrain          = penalty_par0 && penalty_par1;
      double sinsq2th13_mean  = conv::kPrior_Sinsq2th13_Used_Mean; // Used (sinsq2th13)
      double sinsq2th13_sigma = conv::kPrior_Sinsq2th13_Used_RMS; // Used (sinsq2th13)
      
      if(penalty_par0 && !penalty_par1){
        cout << endl;
        cout << "Use penalty term only for free parameter 0:" << endl;
        cout << " sinsq2th13_mean = "  << sinsq2th13_mean  << endl;
        cout << " sinsq2th13_sigma = " << sinsq2th13_sigma << endl;
        cout << endl;
      }
      else if(!penalty_par0 && !penalty_par1){
        cout << endl;
        cout << "No penalty term is used" << endl;
        cout << endl;
      }
      else{
        cerr << endl;
        cerr << "valor::util::PerformMargFit()" << endl;
        cerr << "Penalty term for parameter 1 is not implemented yet" << endl;
        cerr << endl;
        exit(1);
      }
      
      fitmarg->UsePenaltyTerm(penalty_par0,penalty_par1);
      fitmarg->SetPenaltyTerm(sinsq2th13_mean, sinsq2th13_sigma); // to be used after UsePenaltyTerm()
      
      util::MargFit_2D_deltaCP_sinsqth13(runopt,fitmarg,skobs_fakedata,chisqmin,bestfit0,bestfit1,nfits,fitstatus,niterations);
    }
    
    // 2D sinsqth23 Vs dcp fit
    else if(runopt->Freeparname0()=="mixing23" && runopt->Freeparname1()=="deltacp"){
      util::MargFit_2D_sinsqth23_deltaCP(runopt,fitmarg,skobs_fakedata,chisqmin,bestfit0,bestfit1,nfits,fitstatus,niterations);
    }
    
    // 2D dmsq32 Vs dcp fit
    else if(runopt->Freeparname0()=="dmsq23" && runopt->Freeparname1()=="deltacp"){
      util::MargFit_2D_dmsq32_deltaCP(runopt,fitmarg,skobs_fakedata,chisqmin,bestfit0,bestfit1,nfits,fitstatus,niterations);
    }
    
    // 2D sinsqth23 Vs dmsq32 fit
    else if(runopt->Freeparname0()=="mixing23" && runopt->Freeparname1()=="dmsq23"){
      util::MargFit_2D_sinsqth23_dmsq32(runopt,fitmarg,skobs_fakedata,chisqmin,bestfit0,bestfit1,nfits,fitstatus,niterations);
    }
    
    // 2D sinsqth23bar Vs dmsq32bar fit
    else if(runopt->Freeparname0()=="mixing23Bar" && runopt->Freeparname1()=="dmsq23Bar"){
      util::MargFit_2D_sinsqth23_dmsq32(runopt,fitmarg,skobs_fakedata,chisqmin,bestfit0,bestfit1,nfits,fitstatus,niterations);
    }
    
    // 2D sinsqth23 Vs sinsqth23bar fit
    else if(runopt->Freeparname0()=="mixing23" && runopt->Freeparname1()=="mixing23Bar"){
      util::MargFit_2D_sinsqth23_sinsqth23bar(runopt,fitmarg,skobs_fakedata,chisqmin,bestfit0,bestfit1,nfits,fitstatus,niterations);
    }
    
    // 2D dmsq32 Vs dmsq32bar fit
    else if(runopt->Freeparname0()=="dmsq23" && runopt->Freeparname1()=="dmsq23Bar"){
      util::MargFit_2D_dmsq_dmsqbar(runopt,fitmarg,skobs_fakedata,chisqmin,bestfit0,bestfit1,nfits,fitstatus,niterations);
    }
    
    // 1D sinsqth13 fit
    else if(runopt->Freeparname0()=="mixing13" && runopt->Freeparname1()==""){
      util::MargFit_1D_sinsqth13(runopt,fitmarg,skobs_fakedata,chisqmin,bestfit0,bestfit1,nfits,fitstatus,niterations);
    }  
    //
    // 1D sinsqth23 fit
    else if(runopt->Freeparname0()=="mixing23" && runopt->Freeparname1()==""){
      util::MargFit_1D_sinsqth23(runopt,fitmarg,skobs_fakedata,chisqmin,bestfit0,bestfit1,nfits,fitstatus,niterations);
    }
    else if(runopt->Freeparname0()=="mixing23Bar" && runopt->Freeparname1()==""){
      util::MargFit_1D_sinsqth23bar(runopt,fitmarg,skobs_fakedata,chisqmin,bestfit0,bestfit1,nfits,fitstatus,niterations);
    }

    //
    // 1D dmsq32 fit
    //
    else if(runopt->Freeparname0()=="dmsq23" && runopt->Freeparname1()==""){
      util::MargFit_1D_dmsq32(runopt,fitmarg,skobs_fakedata,chisqmin,bestfit0,bestfit1,nfits,fitstatus,niterations);
    }
    else if(runopt->Freeparname0()=="dmsq23Bar" && runopt->Freeparname1()==""){
      util::MargFit_1D_dmsq32bar(runopt,fitmarg,skobs_fakedata,chisqmin,bestfit0,bestfit1,nfits,fitstatus,niterations);
    }
    else{
      cout << endl;
      cout << "/////////////////////////////////////////////////////////////" << endl;
      cout << "//                                                         //" << endl;
      cout << "// Minimization of the marginalized chisq is not performed //" << endl;
      cout << "// though it suppoosed to                                  //" << endl;
      cout << "//                                                         //" << endl;
      cout << "/////////////////////////////////////////////////////////////" << endl;
      cout << endl;
      cout << endl;
      

    }
    
    //
    
    cout << endl;
    cout << "Minimization is complete" << endl;
    cout << endl;
    cout << "Set best-fit parameter values" << endl;
    cout << endl;
    
    // Set best-fit parameters in FitMarginalization object
    // Otherwise error computation is done with respect to a wrong minimum
    fitmarg->SetChisqBestFit(chisqmin);
    fitmarg->SetParamBestFit(0,bestfit0);
    if(fitmarg->IsFit2D()) fitmarg->SetParamBestFit(1,bestfit1);
    cout << endl;
    cout << "bestfit0 = " << bestfit0 << endl;
    cout << "bestfit1 = " << bestfit1 << endl;
    cout << "chisqmin = " << chisqmin << endl;
    cout << endl; 
    
    
    // Build the DeltaChisqMarg after the best-fit Chisq has been set
    // Otherwise DeltaChisqMarg becomes negative if shifted wrt a local minimum
    
    cout << "Build deltaChisqMarg distribution" << endl;
    fitmarg->BuildDeltaChisqMarg(); 

  }
  else{
    cout << "Do not minimize the free parameters but setup the penalty terms" << endl;

    bool penalty_par0       = runopt->PenaltyTerm0();
    bool penalty_par1       = runopt->PenaltyTerm1();
    bool constrain          = penalty_par0 && penalty_par1;
    double sinsq2th13_mean  = conv::kPrior_Sinsq2th13_Used_Mean; // Used (sinsq2th13)
    double sinsq2th13_sigma = conv::kPrior_Sinsq2th13_Used_RMS; // Used (sinsq2th13)
    
    if(penalty_par0 && !penalty_par1){
      cout << endl;
      cout << "Use penalty term only for free parameter 0:" << endl;
      cout << " sinsq2th13_mean = "  << sinsq2th13_mean  << endl;
      cout << " sinsq2th13_sigma = " << sinsq2th13_sigma << endl;
      cout << endl;
    }
    else if(!penalty_par0 && !penalty_par1){
      cout << endl;
      cout << "No penalty term is used" << endl;
      cout << endl;
    }
    else{
      cerr << endl;
      cerr << "valor::util::PerformMargFit()" << endl;
      cerr << "Penalty term for parameter 1 is not implemented yet" << endl;
      cerr << endl;
      exit(1);
    }
    
    fitmarg->UsePenaltyTerm(penalty_par0,penalty_par1);
    fitmarg->SetPenaltyTerm(sinsq2th13_mean, sinsq2th13_sigma); // to be used after UsePenaltyTerm()

  }

}
//____________________________________________________________________________
conv::PhysParamId_t valor::util::GetFreePar(string freeparname)
{

  PhysParamId_t freepar;
  string sfreeparname(freeparname);  
  // NuBar mode
  if(sfreeparname=="mixing23Bar")      freepar = conv::kPhypMixing23Bar;
  else if(sfreeparname=="mixing13Bar") freepar = conv::kPhypMixing13Bar;
  else if(sfreeparname=="dmsq23Bar")   freepar = conv::kPhypDm23sqBar;
  else if(sfreeparname=="dmsq12Bar")   freepar = conv::kPhypDm12sqBar;
  else if(sfreeparname=="mixing12Bar") freepar = conv::kPhypMixing12Bar;
  else if(sfreeparname=="deltacpBar")  freepar = conv::kPhypCPVBar;
  else if(sfreeparname=="betaBar")     freepar = conv::kPhypBetaBar;
  // Nu mode
  else if(sfreeparname=="mixing23")    freepar = conv::kPhypMixing23;
  else if(sfreeparname=="mixing13")    freepar = conv::kPhypMixing13;
  else if(sfreeparname=="dmsq23")      freepar = conv::kPhypDm23sq;
  else if(sfreeparname=="dmsq12")      freepar = conv::kPhypDm12sq;
  else if(sfreeparname=="mixing12")    freepar = conv::kPhypMixing12;
  else if(sfreeparname=="deltacp")     freepar = conv::kPhypCPV;
  
  else{
    freepar = conv::kPhypUndefined;
    cout << "UtlForMacrosT2XK/GetFreePar >> Param: "<< freeparname << " is not defined" << endl;
    exit(1);
  }

  return freepar;
}
//____________________________________________________________________________
void valor::util::MargFit_1D_deltaCP(RunOpt *runopt,FitMarginalization *fitmarg,SuperkObservation *skobs_fakedata,double &chisqmin,double &bestfit0, double &bestfit1,int nfits,int fitstatus,int niterations)
{
  
  //double valsecondbin_X = fitmarg->GetMapChisqRangeFileID_2D()->GetXaxis()->GetBinCenter(2);
  //double valsecondbin_Y = fitmarg->GetMapChisqRangeFileID_2D()->GetYaxis()->GetBinCenter(2);
  //int nbinsY = fitmarg->GetMapChisqRangeFileID_2D()->GetNbinsY();
  //double vallastbin_Y = fitmarg->GetMapChisqRangeFileID_2D()->GetYaxis()->GetBinCenter(nbinsY);

  // Fit 1

  cout << endl;
  cout << "################# Fit 1 (dcp=-3.14) #################" << endl;
    
  //fitmarg->SetMinimIntervalRead(valsecondbin_X+conv::kASmallNum,0.3,-TMath::Pi(),+TMath::Pi(),true);
  fitmarg->SetMinimIntervalRead(-TMath::Pi(),+TMath::Pi(),0.,0.,true);
  //fitmarg->SetStartBin(0.1,-TMath::Pi()); 
  fitmarg->SetStartBin(-TMath::Pi(),0.); 
  //fitmarg->SetNumStartBin(10,1);
  fitmarg->SetNumStartBin(1,1);
  //fitmarg->SetGridStartRange(valsecondbin_X+conv::kASmallNum,0.3,-TMath::Pi(),-TMath::Pi()+conv::kASmallNum);
  fitmarg->SetGridStartRange(-TMath::Pi(),-TMath::Pi()+conv::kASmallNum,0.,0.);
  //fitmarg->SetExtractStartPt("GridFix");
  fitmarg->SetExtractStartPt("Single");
  fitmarg->MinimizeChisqStepRead(skobs_fakedata,true);      

  double chisqmin_1 = (double)fitmarg->GetChisqBestFit();
  double bestfit0_1 = fitmarg->GetParamBestFit(0);
  double bestfit1_1 = -999.;
  if(fitmarg->IsFit2D()) bestfit1_1 = fitmarg->GetParamBestFit(1); 
  int nfits_1 = fitmarg->GetNumFit();
  int fitstatus_1 = fitmarg->GetFitStatus();
  int niterations_1 = fitmarg->GetNumIterations();
  if(chisqmin_1<chisqmin){
    chisqmin = chisqmin_1;
    bestfit0 = bestfit0_1; 
    bestfit1 = bestfit1_1; 
    nfits = nfits_1;    
    fitstatus = fitstatus_1;
    niterations = niterations_1;
  }

  // Fit 2

  cout << endl;
  cout << "################# Fit 2 (dcp=-pi/2) #################" << endl;
    
  fitmarg->SetMinimIntervalRead(-TMath::Pi(),+TMath::Pi(),0.,0.,true);
  fitmarg->SetStartBin(-TMath::Pi()/2.,0.); 
  fitmarg->SetNumStartBin(1,1);
  fitmarg->SetGridStartRange(-TMath::Pi()/2.,-TMath::Pi()/2.+conv::kASmallNum,0.,0.);
  //fitmarg->SetExtractStartPt("GridFix");
  fitmarg->SetExtractStartPt("Single");
  fitmarg->MinimizeChisqStepRead(skobs_fakedata,true);      

  double chisqmin_2 = (double)fitmarg->GetChisqBestFit();
  double bestfit0_2 = fitmarg->GetParamBestFit(0);
  double bestfit1_2 = -999.;
  if(fitmarg->IsFit2D()) bestfit1_2 = fitmarg->GetParamBestFit(1); 
  int nfits_2 = fitmarg->GetNumFit();
  int fitstatus_2 = fitmarg->GetFitStatus();
  int niterations_2 = fitmarg->GetNumIterations();
  if(chisqmin_2<chisqmin){
    chisqmin = chisqmin_2;
    bestfit0 = bestfit0_2; 
    bestfit1 = bestfit1_2; 
    nfits = nfits_2;    
    fitstatus = fitstatus_2;
    niterations = niterations_2;
  }

  // Fit 3

  cout << endl;
  cout << "################# Fit 3 (dcp=+pi/2) #################" << endl;
    
  fitmarg->SetMinimIntervalRead(-TMath::Pi(),+TMath::Pi(),0.,0.,true);
  fitmarg->SetStartBin(+TMath::Pi()/2.,0.); 
  fitmarg->SetNumStartBin(1,1);
  fitmarg->SetGridStartRange(+TMath::Pi()/2.,+TMath::Pi()/2.+conv::kASmallNum,0.,0.);
  //fitmarg->SetExtractStartPt("GridFix");
  fitmarg->SetExtractStartPt("Single");
  fitmarg->MinimizeChisqStepRead(skobs_fakedata,true);      

  double chisqmin_3 = (double)fitmarg->GetChisqBestFit();
  double bestfit0_3 = fitmarg->GetParamBestFit(0);
  double bestfit1_3 = -999.;
  if(fitmarg->IsFit2D()) bestfit1_3 = fitmarg->GetParamBestFit(1); 
  int nfits_3 = fitmarg->GetNumFit();
  int fitstatus_3 = fitmarg->GetFitStatus();
  int niterations_3 = fitmarg->GetNumIterations();
  if(chisqmin_3<chisqmin){
    chisqmin = chisqmin_3;
    bestfit0 = bestfit0_3; 
    bestfit1 = bestfit1_3; 
    nfits = nfits_3;    
    fitstatus = fitstatus_3;
    niterations = niterations_3;
  }

  // Fit 4

  cout << endl;
  cout << "################# Fit 4 (dcp=+3.14) #################" << endl;
    
  fitmarg->SetMinimIntervalRead(-TMath::Pi(),+TMath::Pi(),0.,0.,true);
  fitmarg->SetStartBin(+TMath::Pi(),0.); 
  fitmarg->SetNumStartBin(1,1);
  fitmarg->SetGridStartRange(+TMath::Pi()-conv::kASmallNum,+TMath::Pi(),0.,0.);
  //fitmarg->SetExtractStartPt("GridFix");
  fitmarg->SetExtractStartPt("Single");
  fitmarg->MinimizeChisqStepRead(skobs_fakedata,true);      

  double chisqmin_4 = (double)fitmarg->GetChisqBestFit();
  double bestfit0_4 = fitmarg->GetParamBestFit(0);
  double bestfit1_4 = -999.;
  if(fitmarg->IsFit2D()) bestfit1_4 = fitmarg->GetParamBestFit(1); 
  int nfits_4 = fitmarg->GetNumFit();
  int fitstatus_4 = fitmarg->GetFitStatus();
  int niterations_4 = fitmarg->GetNumIterations();
  if(chisqmin_4<chisqmin){
    chisqmin = chisqmin_4;
    bestfit0 = bestfit0_4; 
    bestfit1 = bestfit1_4; 
    nfits = nfits_4;    
    fitstatus = fitstatus_4;
    niterations = niterations_4;
  }
}
//____________________________________________________________________________
void valor::util::MargFit_1D_sinsqth13(RunOpt *runopt,FitMarginalization *fitmarg,SuperkObservation *skobs_fakedata,double &chisqmin,double &bestfit0, double &bestfit1,int nfits,int fitstatus,int niterations)
{

  // Fit 1

  cout << endl;
  cout << "################# Fit 1 (sinsqth13 = conv::kUsed_SinsqTheta13) #################" << endl;
    
  fitmarg->SetMinimIntervalRead(0.007,0.053,0.,0.,true);
  fitmarg->SetStartBin(conv::kUsed_SinsqTheta13,0.); 
  fitmarg->SetNumStartBin(1,1);
  fitmarg->SetGridStartRange(0.007,0.053,0.,0.);
  fitmarg->SetExtractStartPt("Single");
  fitmarg->MinimizeChisqStepRead(skobs_fakedata,true);      

  double chisqmin_1 = (double)fitmarg->GetChisqBestFit();
  double bestfit0_1 = fitmarg->GetParamBestFit(0);
  double bestfit1_1 = -999.;
  if(fitmarg->IsFit2D()) bestfit1_1 = fitmarg->GetParamBestFit(1); 
  int nfits_1 = fitmarg->GetNumFit();
  int fitstatus_1 = fitmarg->GetFitStatus();
  int niterations_1 = fitmarg->GetNumIterations();
  if(chisqmin_1<chisqmin){
    chisqmin = chisqmin_1;
    bestfit0 = bestfit0_1; 
    bestfit1 = bestfit1_1; 
    nfits = nfits_1;    
    fitstatus = fitstatus_1;
    niterations = niterations_1;
  }

}
//____________________________________________________________________________
void valor::util::MargFit_1D_sinsqth23(RunOpt *runopt,FitMarginalization *fitmarg,SuperkObservation *skobs_fakedata,double &chisqmin,double &bestfit0, double &bestfit1,int nfits,int fitstatus,int niterations)
{

  // Fit 1

  cout << endl;
  cout << "################# Fit 1 (sinsqth23 = conv::kUsed_SinsqTheta23) #################" << endl;
    
  fitmarg->SetMinimIntervalRead(0.3,0.7,0.,0.,true);
  fitmarg->SetStartBin(conv::kUsed_SinsqTheta23,0.); 
  fitmarg->SetNumStartBin(1,1);
  fitmarg->SetGridStartRange(0.3,0.7,0.,0.);
  fitmarg->SetExtractStartPt("Single");
  fitmarg->MinimizeChisqStepRead(skobs_fakedata,true);      

  double chisqmin_1 = (double)fitmarg->GetChisqBestFit();
  double bestfit0_1 = fitmarg->GetParamBestFit(0);
  double bestfit1_1 = -999.;
  if(fitmarg->IsFit2D()) bestfit1_1 = fitmarg->GetParamBestFit(1); 
  int nfits_1 = fitmarg->GetNumFit();
  int fitstatus_1 = fitmarg->GetFitStatus();
  int niterations_1 = fitmarg->GetNumIterations();
  if(chisqmin_1<chisqmin){
    chisqmin = chisqmin_1;
    bestfit0 = bestfit0_1; 
    bestfit1 = bestfit1_1; 
    nfits = nfits_1;    
    fitstatus = fitstatus_1;
    niterations = niterations_1;
  }

}
//____________________________________________________________________________
void valor::util::MargFit_1D_sinsqth23bar(RunOpt *runopt,FitMarginalization *fitmarg,SuperkObservation *skobs_fakedata,double &chisqmin,double &bestfit0, double &bestfit1,int nfits,int fitstatus,int niterations)
{
  
  // Fit 1
  
  cout << endl;
  cout << "################# Fit 1 (sinsqth23 = conv::kUsed_SinsqTheta23) #################" << endl;
  
  double valsecondbin_X = fitmarg->GetMapChisqRangeFileID_2D()->GetXaxis()->GetBinCenter(2);
  
  int nbinsX = fitmarg->GetMapChisqRangeFileID_2D()->GetNbinsX();
  double valfirstbin_X = fitmarg->GetMapChisqRangeFileID_2D()->GetXaxis()->GetBinCenter(1);
  double vallastbin_X = fitmarg->GetMapChisqRangeFileID_2D()->GetXaxis()->GetBinCenter(nbinsX);
  
  fitmarg->SetMinimIntervalRead(valsecondbin_X+conv::kASmallNum,vallastbin_X,0.,0.,true);
  fitmarg->SetStartBin(conv::kUsed_SinsqTheta23,0.);
  fitmarg->SetNumStartBin(1,1);
  fitmarg->SetMinimIntervalRead(valsecondbin_X+conv::kASmallNum,vallastbin_X,0.,0.);
  fitmarg->SetExtractStartPt("Single");
  fitmarg->MinimizeChisqStepRead(skobs_fakedata,true);
  
  double chisqmin_1 = (double)fitmarg->GetChisqBestFit();
  double bestfit0_1 = fitmarg->GetParamBestFit(0);
  double bestfit1_1 = -999.;
  if(fitmarg->IsFit2D()) bestfit1_1 = fitmarg->GetParamBestFit(1);
  int nfits_1 = fitmarg->GetNumFit();
  int fitstatus_1 = fitmarg->GetFitStatus();
  int niterations_1 = fitmarg->GetNumIterations();
  if(chisqmin_1<chisqmin){
    chisqmin = chisqmin_1;
    bestfit0 = bestfit0_1;
    bestfit1 = bestfit1_1;
    nfits = nfits_1;
    fitstatus = fitstatus_1;
    niterations = niterations_1;
  }
  
}
//____________________________________________________________________________
void valor::util::MargFit_1D_dmsq32(RunOpt *runopt,FitMarginalization *fitmarg,SuperkObservation *skobs_fakedata,double &chisqmin,double &bestfit0, double &bestfit1,int nfits,int fitstatus,int niterations)
{

  // Fit 1

  cout << endl;
  cout << "################# Fit 1 (dmsq32 = conv::kUsed_Dm23sq) #################" << endl;
    
  fitmarg->SetMinimIntervalRead(2.3E-3,2.7E-3,0.,0.,true);
  fitmarg->SetStartBin(conv::kUsed_Dm23sq,0.); 
  fitmarg->SetNumStartBin(1,1);
  fitmarg->SetGridStartRange(2.3E-3,2.7E-3,0.,0.);
  fitmarg->SetExtractStartPt("Single");
  fitmarg->MinimizeChisqStepRead(skobs_fakedata,true);      

  double chisqmin_1 = (double)fitmarg->GetChisqBestFit();
  double bestfit0_1 = fitmarg->GetParamBestFit(0);
  double bestfit1_1 = -999.;
  if(fitmarg->IsFit2D()) bestfit1_1 = fitmarg->GetParamBestFit(1); 
  int nfits_1 = fitmarg->GetNumFit();
  int fitstatus_1 = fitmarg->GetFitStatus();
  int niterations_1 = fitmarg->GetNumIterations();
  if(chisqmin_1<chisqmin){
    chisqmin = chisqmin_1;
    bestfit0 = bestfit0_1; 
    bestfit1 = bestfit1_1; 
    nfits = nfits_1;    
    fitstatus = fitstatus_1;
    niterations = niterations_1;
  }

}
//____________________________________________________________________________
void valor::util::MargFit_1D_dmsq32bar(RunOpt *runopt,FitMarginalization *fitmarg,SuperkObservation *skobs_fakedata,double &chisqmin,double &bestfit0, double &bestfit1,int nfits,int fitstatus,int niterations)
{
  
  // Fit 1
  
  cout << endl;
  cout << "################# Fit 1 (dmsq32 = conv::kUsed_Dm23sq) #################" << endl;
  
  double valsecondbin_X = fitmarg->GetMapChisqRangeFileID_2D()->GetXaxis()->GetBinCenter(2);
  
  int nbinsX = fitmarg->GetMapChisqRangeFileID_2D()->GetNbinsX();
  double valfirstbin_X = fitmarg->GetMapChisqRangeFileID_2D()->GetXaxis()->GetBinCenter(1);
  double vallastbin_X = fitmarg->GetMapChisqRangeFileID_2D()->GetXaxis()->GetBinCenter(nbinsX);
  
  fitmarg->SetMinimIntervalRead(valsecondbin_X+conv::kASmallNum,vallastbin_X,0.,0.,true);
  fitmarg->SetStartBin(conv::kUsed_Dm23sq,0.);
  fitmarg->SetNumStartBin(1,1);
  fitmarg->SetMinimIntervalRead(valsecondbin_X+conv::kASmallNum,vallastbin_X,0.,0.);
  fitmarg->SetExtractStartPt("Single");
  fitmarg->MinimizeChisqStepRead(skobs_fakedata,true);
  
  
  double chisqmin_1 = (double)fitmarg->GetChisqBestFit();
  double bestfit0_1 = fitmarg->GetParamBestFit(0);
  double bestfit1_1 = -999.;
  if(fitmarg->IsFit2D()) bestfit1_1 = fitmarg->GetParamBestFit(1);
  int nfits_1 = fitmarg->GetNumFit();
  int fitstatus_1 = fitmarg->GetFitStatus();
  int niterations_1 = fitmarg->GetNumIterations();
  if(chisqmin_1<chisqmin){
    chisqmin = chisqmin_1;
    bestfit0 = bestfit0_1;
    bestfit1 = bestfit1_1;
    nfits = nfits_1;
    fitstatus = fitstatus_1;
    niterations = niterations_1;
  }
  
}
//____________________________________________________________________________
void valor::util::MargFit_1D_betaBar(RunOpt *runopt,FitMarginalization *fitmarg,SuperkObservation *skobs_fakedata,double &chisqmin,double &bestfit0, double &bestfit1,int nfits,int fitstatus,int niterations)
{

  // Fit. Just set the fitmarginalization

  cout << endl;
  cout << "################# Fit #################" << endl;
  
  int nbinsX = fitmarg->GetMapChisqRangeFileID_2D()->GetNbinsX();
  double valfirstbin_X = fitmarg->GetMapChisqRangeFileID_2D()->GetXaxis()->GetBinCenter(1);
  double vallastbin_X = fitmarg->GetMapChisqRangeFileID_2D()->GetXaxis()->GetBinCenter(nbinsX);

  fitmarg->SetMinimIntervalRead(valfirstbin_X,vallastbin_X,0.,0.,true);
  fitmarg->SetStartBin(0.,1.);
  fitmarg->SetNumStartBin(1,1);
  fitmarg->SetGridStartRange(valfirstbin_X,vallastbin_X,0.,0.);
  fitmarg->SetExtractStartPt("Single");
  fitmarg->MinimizeChisqStepRead(skobs_fakedata,true);  

  double chisqmin_1 = (double)fitmarg->GetChisqBestFit();
  double bestfit0_1 = fitmarg->GetParamBestFit(0);
  double bestfit1_1 = -999.;
  if(fitmarg->IsFit2D()) bestfit1_1 = fitmarg->GetParamBestFit(1); 
  int nfits_1 = fitmarg->GetNumFit();
  int fitstatus_1 = fitmarg->GetFitStatus();
  int niterations_1 = fitmarg->GetNumIterations();
  if(chisqmin_1<chisqmin){
    chisqmin = chisqmin_1;
    bestfit0 = bestfit0_1; 
    bestfit1 = bestfit1_1; 
    nfits = nfits_1;    
    fitstatus = fitstatus_1;
    niterations = niterations_1;
  }  
}
//____________________________________________________________________________
void valor::util::MargFit_1D_betaBarRateOnly(RunOpt *runopt,FitMarginalization *fitmarg,SuperkObservation *skobs_fakedata,double &chisqmin,double &bestfit0, double &bestfit1,int nfits,int fitstatus,int niterations)
{
  
  // Fit. Just set the fitmarginalization
  
  cout << endl;
  cout << "################# Fit #################" << endl;
  
  fitmarg->SetMinimIntervalRead(0,1.,0.,0.,true);
  fitmarg->SetStartBin(0.,1.);
  fitmarg->SetNumStartBin(1,1);
  fitmarg->SetGridStartRange(0.,1.,0.,0.);
  fitmarg->SetExtractStartPt("Single");
  fitmarg->MinimizeChisqStepRead(skobs_fakedata,true);
  
  double chisqmin_1 = (double)fitmarg->GetChisqBestFit();
  double bestfit0_1 = fitmarg->GetParamBestFit(0);
  double bestfit1_1 = -999.;
  if(fitmarg->IsFit2D()) bestfit1_1 = fitmarg->GetParamBestFit(1);
  int nfits_1 = fitmarg->GetNumFit();
  int fitstatus_1 = fitmarg->GetFitStatus();
  int niterations_1 = fitmarg->GetNumIterations();
  if(chisqmin_1<chisqmin){
    chisqmin = chisqmin_1;
    bestfit0 = bestfit0_1;
    bestfit1 = bestfit1_1;
    nfits = nfits_1;
    fitstatus = fitstatus_1;
    niterations = niterations_1;
  }

}
//____________________________________________________________________________
void valor::util::MargFit_2D_sinsqth13_deltaCP(RunOpt *runopt,FitMarginalization *fitmarg,SuperkObservation *skobs_fakedata,double &chisqmin,double &bestfit0, double &bestfit1,int nfits,int fitstatus,int niterations)
{
  
  double valsecondbin_X = fitmarg->GetMapChisqRangeFileID_2D()->GetXaxis()->GetBinCenter(2);
  double valsecondbin_Y = fitmarg->GetMapChisqRangeFileID_2D()->GetYaxis()->GetBinCenter(2);

  int nbinsX = fitmarg->GetMapChisqRangeFileID_2D()->GetNbinsX();
  double valfirstbin_X = fitmarg->GetMapChisqRangeFileID_2D()->GetXaxis()->GetBinCenter(1);
  double vallastbin_X = fitmarg->GetMapChisqRangeFileID_2D()->GetXaxis()->GetBinCenter(nbinsX);
  
  int nbinsY = fitmarg->GetMapChisqRangeFileID_2D()->GetNbinsY();
  double valfirstbin_Y = fitmarg->GetMapChisqRangeFileID_2D()->GetYaxis()->GetBinCenter(1);
  double vallastbin_Y = fitmarg->GetMapChisqRangeFileID_2D()->GetYaxis()->GetBinCenter(nbinsY);
  

  //double start_sinsqth13_width = (vallastbin_X - valfirstbin_X) / 3.;
  //double start_sinsqth13_first = valfirstbin_X + start_sinsqth13_width;
  //double start_sinsqth13_last  = valfirstbin_X + 2*start_sinsqth13_width;
  
  double start_sinsqth13_first = valfirstbin_X;
  double start_sinsqth13_last  = vallastbin_X;
  
  //cout << endl;
  //cout << "start_sinsqth13_first = " << start_sinsqth13_first << endl;
  //cout << "start_sinsqth13_last  = " << start_sinsqth13_last << endl;
  //cout << endl;


  // Fit 1
  
  cout << endl;
  cout << "/////////////// Fit 1 (mixing13=0) ////////////////" << endl;
  
  //double chisqmin_1 = (double)fitmarg->GetChisqMargTrueRead(0.,0.,skobs_fakedata,true);
  double chisqmin_1 = (double)fitmarg->GetChisqMargTrueRead(valfirstbin_X,valfirstbin_Y,skobs_fakedata,true);
  double bestfit0_1 = 0.;
  double bestfit1_1 = -999.;
  if(fitmarg->IsFit2D()) bestfit1_1 = 0.;     
  int nfits_1 = 1; // it's not a fit but single point
  int fitstatus_1 = 1; // it's not a fit but single point
  int niterations_1 = 1; // it's not a fit but single point
  if(chisqmin_1<chisqmin){
    chisqmin = chisqmin_1;
    bestfit0 = bestfit0_1; 
    bestfit1 = bestfit1_1; 
    nfits = nfits_1;    
    fitstatus = fitstatus_1;
    niterations = niterations_1;
  }
  
  cout << endl;
  cout << "True parameters: " << endl;
  cout << " - par0 = " << bestfit0_1 << endl;
  cout << " - par1 = " << bestfit1_1 << endl;
  cout << "True chisq = " << chisqmin_1;
  cout << endl;

  // Fit 2
  
  cout << endl;
  cout << "/////////////// Fit 2 (mixing13=first, dcp=-3.14) ////////////////" << endl;
  
  //fitmarg->SetMinimIntervalRead(valsecondbin_X+conv::kASmallNum,0.4,-TMath::Pi(),+TMath::Pi(),true);
  fitmarg->SetMinimIntervalRead(valsecondbin_X+conv::kASmallNum,vallastbin_X,valfirstbin_Y,vallastbin_Y,true);
  fitmarg->SetStartBin(start_sinsqth13_first,-TMath::Pi()); 
  fitmarg->SetNumStartBin(1,1);
  //fitmarg->SetGridStartRange(valsecondbin_X+conv::kASmallNum,0.4,-TMath::Pi(),-TMath::Pi()+conv::kASmallNum);
  fitmarg->SetGridStartRange(valsecondbin_X+conv::kASmallNum,vallastbin_X,valfirstbin_Y,vallastbin_Y+conv::kASmallNum);
  fitmarg->SetExtractStartPt("Single");
  fitmarg->MinimizeChisqStepRead(skobs_fakedata,true);      

  double chisqmin_2 = (double)fitmarg->GetChisqBestFit();
  double bestfit0_2 = fitmarg->GetParamBestFit(0);
  double bestfit1_2 = -999.;
  if(fitmarg->IsFit2D()) bestfit1_2 = fitmarg->GetParamBestFit(1); 
  int nfits_2 = fitmarg->GetNumFit();
  int fitstatus_2 = fitmarg->GetFitStatus();
  int niterations_2 = fitmarg->GetNumIterations();
  if(chisqmin_2<chisqmin){
    chisqmin = chisqmin_2;
    bestfit0 = bestfit0_2; 
    bestfit1 = bestfit1_2; 
    nfits = nfits_2;    
    fitstatus = fitstatus_2;
    niterations = niterations_2;
  }
    
  // Fit 3
  
  cout << endl;
  cout << "/////////////// Fit 3 (mixing13=start_sinsqth13_first, dcp=-pi/2.) ////////////////" << endl;
  
  //fitmarg->SetMinimIntervalRead(valsecondbin_X+conv::kASmallNum,0.4,-TMath::Pi(),+TMath::Pi(),true);
  fitmarg->SetMinimIntervalRead(valsecondbin_X+conv::kASmallNum,vallastbin_X,valfirstbin_Y,vallastbin_Y,true);
  fitmarg->SetStartBin(start_sinsqth13_first,-TMath::Pi()/2.); 
  fitmarg->SetNumStartBin(1,1);
  //fitmarg->SetGridStartRange(valsecondbin_X+conv::kASmallNum,0.4,-TMath::Pi(),-TMath::Pi()+conv::kASmallNum);
  fitmarg->SetGridStartRange(valsecondbin_X+conv::kASmallNum,vallastbin_X,valfirstbin_Y,vallastbin_Y+conv::kASmallNum);
  fitmarg->SetExtractStartPt("Single");
  fitmarg->MinimizeChisqStepRead(skobs_fakedata,true);      

  double chisqmin_3 = (double)fitmarg->GetChisqBestFit();
  double bestfit0_3 = fitmarg->GetParamBestFit(0);
  double bestfit1_3 = -999.;
  if(fitmarg->IsFit2D()) bestfit1_3 = fitmarg->GetParamBestFit(1); 
  int nfits_3 = fitmarg->GetNumFit();
  int fitstatus_3 = fitmarg->GetFitStatus();
  int niterations_3 = fitmarg->GetNumIterations();
  if(chisqmin_3<chisqmin){
    chisqmin = chisqmin_3;
    bestfit0 = bestfit0_3; 
    bestfit1 = bestfit1_3; 
    nfits = nfits_3;    
    fitstatus = fitstatus_3;
    niterations = niterations_3;
  }
  
  
  // Fit 4
  
  cout << endl;
  cout << "/////////////// Fit 4 (mixing13=first, dcp=+pi/2.) ////////////////" << endl;
  
  //fitmarg->SetMinimIntervalRead(valsecondbin_X+conv::kASmallNum,0.4,-TMath::Pi(),+TMath::Pi(),true);
  fitmarg->SetMinimIntervalRead(valsecondbin_X+conv::kASmallNum,vallastbin_X,valfirstbin_Y,vallastbin_Y,true);
  fitmarg->SetStartBin(start_sinsqth13_first,+TMath::Pi()/2.); 
  fitmarg->SetNumStartBin(1,1);
  //fitmarg->SetGridStartRange(valsecondbin_X+conv::kASmallNum,0.4,-TMath::Pi(),-TMath::Pi()+conv::kASmallNum);
  fitmarg->SetGridStartRange(valsecondbin_X+conv::kASmallNum,vallastbin_X,valfirstbin_Y,vallastbin_Y+conv::kASmallNum);
  fitmarg->SetExtractStartPt("Single");
  fitmarg->MinimizeChisqStepRead(skobs_fakedata,true);      
  
  double chisqmin_4 = (double)fitmarg->GetChisqBestFit();
  double bestfit0_4 = fitmarg->GetParamBestFit(0);
  double bestfit1_4 = -999.;
  if(fitmarg->IsFit2D()) bestfit1_4 = fitmarg->GetParamBestFit(1); 
  int nfits_4 = fitmarg->GetNumFit();
  int fitstatus_4 = fitmarg->GetFitStatus();
  int niterations_4 = fitmarg->GetNumIterations();
  if(chisqmin_4<chisqmin){
    chisqmin = chisqmin_4;
    bestfit0 = bestfit0_4; 
    bestfit1 = bestfit1_4; 
    nfits = nfits_4;    
    fitstatus = fitstatus_4;
    niterations = niterations_4;
  }
  
  
  // Fit 5
  
  cout << endl;
  cout << "/////////////// Fit 5 (mixing13=first, dcp=+3.14) ////////////////" << endl;
  
  //fitmarg->SetMinimIntervalRead(valsecondbin_X+conv::kASmallNum,0.4,-TMath::Pi(),+TMath::Pi(),true);
  fitmarg->SetMinimIntervalRead(valsecondbin_X+conv::kASmallNum,vallastbin_X,valfirstbin_Y,vallastbin_Y,true);
  fitmarg->SetStartBin(start_sinsqth13_first,+TMath::Pi()); 
  fitmarg->SetNumStartBin(1,1);
  //fitmarg->SetGridStartRange(valsecondbin_X+conv::kASmallNum,0.4,-TMath::Pi(),-TMath::Pi()+conv::kASmallNum);
  fitmarg->SetGridStartRange(valsecondbin_X+conv::kASmallNum,vallastbin_X,valfirstbin_Y,vallastbin_Y+conv::kASmallNum);
  fitmarg->SetExtractStartPt("Single");
  fitmarg->MinimizeChisqStepRead(skobs_fakedata,true);      
  
  double chisqmin_5 = (double)fitmarg->GetChisqBestFit();
  double bestfit0_5 = fitmarg->GetParamBestFit(0);
  double bestfit1_5 = -999.;
  if(fitmarg->IsFit2D()) bestfit1_5 = fitmarg->GetParamBestFit(1); 
  int nfits_5 = fitmarg->GetNumFit();
  int fitstatus_5 = fitmarg->GetFitStatus();
  int niterations_5 = fitmarg->GetNumIterations();
  if(chisqmin_5<chisqmin){
    chisqmin = chisqmin_5;
    bestfit0 = bestfit0_5; 
    bestfit1 = bestfit1_5; 
    nfits = nfits_5;    
    fitstatus = fitstatus_5;
    niterations = niterations_5;
  }
  

  // Fit 6
  
  cout << endl;
  cout << "/////////////// Fit 6 (mixing13=last, dcp=-3.14) ////////////////" << endl;
  
  //fitmarg->SetMinimIntervalRead(valsecondbin_X+conv::kASmallNum,0.4,-TMath::Pi(),+TMath::Pi(),true);
  fitmarg->SetMinimIntervalRead(valsecondbin_X+conv::kASmallNum,vallastbin_X,valfirstbin_Y,vallastbin_Y,true);
  fitmarg->SetStartBin(start_sinsqth13_last,-TMath::Pi()); 
  fitmarg->SetNumStartBin(1,1);
  //fitmarg->SetGridStartRange(valsecondbin_X+conv::kASmallNum,0.4,-TMath::Pi(),-TMath::Pi()+conv::kASmallNum);
  fitmarg->SetGridStartRange(valsecondbin_X+conv::kASmallNum,vallastbin_X,valfirstbin_Y,vallastbin_Y+conv::kASmallNum);
  fitmarg->SetExtractStartPt("Single");
  fitmarg->MinimizeChisqStepRead(skobs_fakedata,true);      
  
  double chisqmin_6 = (double)fitmarg->GetChisqBestFit();
  double bestfit0_6 = fitmarg->GetParamBestFit(0);
  double bestfit1_6 = -999.;
  if(fitmarg->IsFit2D()) bestfit1_6 = fitmarg->GetParamBestFit(1); 
  int nfits_6 = fitmarg->GetNumFit();
  int fitstatus_6 = fitmarg->GetFitStatus();
  int niterations_6 = fitmarg->GetNumIterations();
  if(chisqmin_6<chisqmin){
    chisqmin = chisqmin_6;
    bestfit0 = bestfit0_6; 
    bestfit1 = bestfit1_6; 
    nfits = nfits_6;    
    fitstatus = fitstatus_6;
    niterations = niterations_6;
  }
    
  // Fit 7
  
  cout << endl;
  cout << "/////////////// Fit 3 (mixing13=last, dcp=-pi/2.) ////////////////" << endl;
  
  //fitmarg->SetMinimIntervalRead(valsecondbin_X+conv::kASmallNum,0.4,-TMath::Pi(),+TMath::Pi(),true);
  fitmarg->SetMinimIntervalRead(valsecondbin_X+conv::kASmallNum,vallastbin_X,valfirstbin_Y,vallastbin_Y,true);
  fitmarg->SetStartBin(start_sinsqth13_last,-TMath::Pi()/2.); 
  fitmarg->SetNumStartBin(1,1);
  //fitmarg->SetGridStartRange(valsecondbin_X+conv::kASmallNum,0.4,-TMath::Pi(),-TMath::Pi()+conv::kASmallNum);
  fitmarg->SetGridStartRange(valsecondbin_X+conv::kASmallNum,vallastbin_X,valfirstbin_Y,vallastbin_Y+conv::kASmallNum);
  fitmarg->SetExtractStartPt("Single");
  fitmarg->MinimizeChisqStepRead(skobs_fakedata,true);      
  
  double chisqmin_7 = (double)fitmarg->GetChisqBestFit();
  double bestfit0_7 = fitmarg->GetParamBestFit(0);
  double bestfit1_7 = -999.;
  if(fitmarg->IsFit2D()) bestfit1_7 = fitmarg->GetParamBestFit(1); 
  int nfits_7 = fitmarg->GetNumFit();
  int fitstatus_7 = fitmarg->GetFitStatus();
  int niterations_7 = fitmarg->GetNumIterations();
  if(chisqmin_7<chisqmin){
    chisqmin = chisqmin_7;
    bestfit0 = bestfit0_7; 
    bestfit1 = bestfit1_7; 
    nfits = nfits_7;    
    fitstatus = fitstatus_7;
    niterations = niterations_7;
  }
  
  
  // Fit 8
  
  cout << endl;
  cout << "/////////////// Fit 8 (mixing13=last, dcp=+pi/2.) ////////////////" << endl;
  
  //fitmarg->SetMinimIntervalRead(valsecondbin_X+conv::kASmallNum,0.4,-TMath::Pi(),+TMath::Pi(),true);
  fitmarg->SetMinimIntervalRead(valsecondbin_X+conv::kASmallNum,vallastbin_X,valfirstbin_Y,vallastbin_Y,true);
  fitmarg->SetStartBin(start_sinsqth13_last,+TMath::Pi()/2.); 
  fitmarg->SetNumStartBin(1,1);
  //fitmarg->SetGridStartRange(valsecondbin_X+conv::kASmallNum,0.4,-TMath::Pi(),-TMath::Pi()+conv::kASmallNum);
  fitmarg->SetGridStartRange(valsecondbin_X+conv::kASmallNum,vallastbin_X,valfirstbin_Y,vallastbin_Y+conv::kASmallNum);
  fitmarg->SetExtractStartPt("Single");
  fitmarg->MinimizeChisqStepRead(skobs_fakedata,true);      
  
  double chisqmin_8 = (double)fitmarg->GetChisqBestFit();
  double bestfit0_8 = fitmarg->GetParamBestFit(0);
  double bestfit1_8 = -999.;
  if(fitmarg->IsFit2D()) bestfit1_8 = fitmarg->GetParamBestFit(1); 
  int nfits_8 = fitmarg->GetNumFit();
  int fitstatus_8 = fitmarg->GetFitStatus();
  int niterations_8 = fitmarg->GetNumIterations();
  if(chisqmin_8<chisqmin){
    chisqmin = chisqmin_8;
    bestfit0 = bestfit0_8; 
    bestfit1 = bestfit1_8; 
    nfits = nfits_8;    
    fitstatus = fitstatus_8;
    niterations = niterations_8;
  }
  
  
  // Fit 9
  
  cout << endl;
  cout << "/////////////// Fit 9 (mixing13=last, dcp=+3.14) ////////////////" << endl;
  
  //fitmarg->SetMinimIntervalRead(valsecondbin_X+conv::kASmallNum,0.4,-TMath::Pi(),+TMath::Pi(),true);
  fitmarg->SetMinimIntervalRead(valsecondbin_X+conv::kASmallNum,vallastbin_X,valfirstbin_Y,vallastbin_Y,true);
  fitmarg->SetStartBin(start_sinsqth13_last,+TMath::Pi()); 
  fitmarg->SetNumStartBin(1,1);
  //fitmarg->SetGridStartRange(valsecondbin_X+conv::kASmallNum,0.4,-TMath::Pi(),-TMath::Pi()+conv::kASmallNum);
  fitmarg->SetGridStartRange(valsecondbin_X+conv::kASmallNum,vallastbin_X,valfirstbin_Y,vallastbin_Y+conv::kASmallNum);
  fitmarg->SetExtractStartPt("Single");
  fitmarg->MinimizeChisqStepRead(skobs_fakedata,true);      
  
  double chisqmin_9 = (double)fitmarg->GetChisqBestFit();
  double bestfit0_9 = fitmarg->GetParamBestFit(0);
  double bestfit1_9 = -999.;
  if(fitmarg->IsFit2D()) bestfit1_9 = fitmarg->GetParamBestFit(1); 
  int nfits_9 = fitmarg->GetNumFit();
  int fitstatus_9 = fitmarg->GetFitStatus();
  int niterations_9 = fitmarg->GetNumIterations();
  if(chisqmin_9<chisqmin){
    chisqmin = chisqmin_9;
    bestfit0 = bestfit0_9; 
    bestfit1 = bestfit1_9; 
    nfits = nfits_9;    
    fitstatus = fitstatus_9;
    niterations = niterations_9;
  }


  
  
  
  
  
  // // Fit 2
  // cout << endl;
  // cout << "/////////////// Fit 2 (mixing13>0, dcp=-3.14) ////////////////" << endl;
  // fitmarg->SetMinimIntervalRead(valsecondbin_X+conv::kASmallNum,0.3,-TMath::Pi(),+TMath::Pi(),true);
  // fitmarg->SetStartBin(0.1,-TMath::Pi()); 
  // fitmarg->SetNumStartBin(10,1);
  // fitmarg->SetGridStartRange(valsecondbin_X+conv::kASmallNum,0.3,-TMath::Pi(),-TMath::Pi()+conv::kASmallNum);
  // fitmarg->SetExtractStartPt("GridFix");
  // fitmarg->MinimizeChisqStepRead(skobs_fakedata,true);      
  // double chisqmin_2 = (double)fitmarg->GetChisqBestFit();
  // double bestfit0_2 = fitmarg->GetParamBestFit(0);
  // double bestfit1_2 = -999.;
  // if(fitmarg->IsFit2D()) bestfit1_2 = fitmarg->GetParamBestFit(1); 
  // int nfits_2 = fitmarg->GetNumFit();
  // int fitstatus_2 = fitmarg->GetFitStatus();
  // int niterations_2 = fitmarg->GetNumIterations();
  // if(chisqmin_2<chisqmin){
  //   chisqmin = chisqmin_2;
  //   bestfit0 = bestfit0_2; 
  //   bestfit1 = bestfit1_2; 
  //   nfits = nfits_2;    
  //   fitstatus = fitstatus_2;
  //   niterations = niterations_2;
  // }
  
  // // Fit 3
  // cout << endl;
  // cout << "/////////////// Fit 3 (mixing13>0, start dcp = 0.) ////////////////" << endl;
  // int bindcp0_Y = fitmarg->GetMapChisqRangeFileID_2D()->GetYaxis()->FindBin(0.);
  // double valdcp0_Y = fitmarg->GetMapChisqRangeFileID_2D()->GetYaxis()->GetBinCenter(bindcp0_Y);
  // double binwidth = fitmarg->GetMapChisqRangeFileID_2D()->GetYaxis()->GetBinWidth(1);
  // fitmarg->SetMinimIntervalRead(valsecondbin_X+conv::kASmallNum,0.3,-TMath::Pi(),+TMath::Pi(),true);
  // fitmarg->SetStartBin(0.1,-TMath::Pi()); 
  // fitmarg->SetNumStartBin(10,1);
  // fitmarg->SetGridStartRange(valsecondbin_X+conv::kASmallNum,0.3,valdcp0_Y-binwidth,valdcp0_Y+binwidth);
  // fitmarg->SetExtractStartPt("GridFix");
  // fitmarg->MinimizeChisqStepRead(skobs_fakedata,true);          
  // double chisqmin_3 = (double)fitmarg->GetChisqBestFit();
  // double bestfit0_3 = fitmarg->GetParamBestFit(0);
  // double bestfit1_3 = -999.;
  // if(fitmarg->IsFit2D()) bestfit1_3 = fitmarg->GetParamBestFit(1); 
  // int nfits_3 = fitmarg->GetNumFit();
  // int fitstatus_3 = fitmarg->GetFitStatus();
  // int niterations_3 = fitmarg->GetNumIterations();
  // if(chisqmin_3<chisqmin){
  //   chisqmin = chisqmin_3;
  //   bestfit0 = bestfit0_3; 
  //   bestfit1 = bestfit1_3; 
  //   nfits = nfits_3;    
  //   fitstatus = fitstatus_3;
  //   niterations = niterations_3;
  // }
  
  // // Fit 4
  // cout << endl;
  // cout << "/////////////// Fit 4 (mixing13>0, start dcp=+3.14) ////////////////" << endl;
  // fitmarg->SetMinimIntervalRead(valsecondbin_X+conv::kASmallNum,0.3,-TMath::Pi(),+TMath::Pi(),true);
  // fitmarg->SetStartBin(0.1,-TMath::Pi()); 
  // fitmarg->SetNumStartBin(10,1);
  // fitmarg->SetGridStartRange(valsecondbin_X+conv::kASmallNum,0.3,vallastbin_Y-conv::kASmallNum,+TMath::Pi());
  // fitmarg->SetExtractStartPt("GridFix");
  // fitmarg->MinimizeChisqStepRead(skobs_fakedata,true);      
  // double chisqmin_4 = (double)fitmarg->GetChisqBestFit();
  // double bestfit0_4 = fitmarg->GetParamBestFit(0);
  // double bestfit1_4 = -999.;
  // if(fitmarg->IsFit2D()) bestfit1_4 = fitmarg->GetParamBestFit(1); 
  // int nfits_4 = fitmarg->GetNumFit();
  // int fitstatus_4 = fitmarg->GetFitStatus();
  // int niterations_4 = fitmarg->GetNumIterations();
  // if(chisqmin_4<chisqmin){
  //   chisqmin = chisqmin_4;
  //   bestfit0 = bestfit0_4; 
  //   bestfit1 = bestfit1_4; 
  //   nfits = nfits_4;    
  //   fitstatus = fitstatus_4;
  //   niterations = niterations_4;
  // }
}
//____________________________________________________________________________
void valor::util::MargFit_2D_deltaCP_sinsqth13(RunOpt *runopt,FitMarginalization *fitmarg,SuperkObservation *skobs_fakedata,double &chisqmin,double &bestfit0, double &bestfit1,int nfits,int fitstatus,int niterations)
{
  
  double valsecondbin_X = fitmarg->GetMapChisqRangeFileID_2D()->GetXaxis()->GetBinCenter(2);
  double valsecondbin_Y = fitmarg->GetMapChisqRangeFileID_2D()->GetYaxis()->GetBinCenter(2);
  
  int nbinsX = fitmarg->GetMapChisqRangeFileID_2D()->GetNbinsX();
  double valfirstbin_X = fitmarg->GetMapChisqRangeFileID_2D()->GetXaxis()->GetBinCenter(1);
  double vallastbin_X = fitmarg->GetMapChisqRangeFileID_2D()->GetXaxis()->GetBinCenter(nbinsX);
  
  int nbinsY = fitmarg->GetMapChisqRangeFileID_2D()->GetNbinsY();
  double valfirstbin_Y = fitmarg->GetMapChisqRangeFileID_2D()->GetYaxis()->GetBinCenter(1);
  double vallastbin_Y = fitmarg->GetMapChisqRangeFileID_2D()->GetYaxis()->GetBinCenter(nbinsY);
  
  
  //double start_sinsqth13_width = (vallastbin_X - valfirstbin_X) / 3.;
  //double start_sinsqth13_first = valfirstbin_X + start_sinsqth13_width;
  //double start_sinsqth13_last  = valfirstbin_X + 2*start_sinsqth13_width;
  
  double start_sinsqth13_first = valfirstbin_Y;
  double start_sinsqth13_last  = vallastbin_Y;
  
  //cout << endl;
  //cout << "start_sinsqth13_first = " << start_sinsqth13_first << endl;
  //cout << "start_sinsqth13_last  = " << start_sinsqth13_last << endl;
  //cout << endl;
  
  
  // Fit 1
  
  cout << endl;
  cout << "/////////////// Fit 1 (mixing13=0) ////////////////" << endl;
  
  //double chisqmin_1 = (double)fitmarg->GetChisqMargTrueRead(0.,0.,skobs_fakedata,true);
  double chisqmin_1 = (double)fitmarg->GetChisqMargTrueRead(valfirstbin_X,valfirstbin_Y,skobs_fakedata,true);
  double bestfit0_1 = 0.;
  double bestfit1_1 = -999.;
  if(fitmarg->IsFit2D()) bestfit1_1 = 0.;
  int nfits_1 = 1; // it's not a fit but single point
  int fitstatus_1 = 1; // it's not a fit but single point
  int niterations_1 = 1; // it's not a fit but single point
  if(chisqmin_1<chisqmin){
    chisqmin = chisqmin_1;
    bestfit0 = bestfit0_1;
    bestfit1 = bestfit1_1;
    nfits = nfits_1;
    fitstatus = fitstatus_1;
    niterations = niterations_1;
  }
  
  cout << endl;
  cout << "True parameters: " << endl;
  cout << " - par0 = " << bestfit0_1 << endl;
  cout << " - par1 = " << bestfit1_1 << endl;
  cout << "True chisq = " << chisqmin_1;
  cout << endl;
  
  // Fit 2
  
  cout << endl;
  cout << "/////////////// Fit 2 (mixing13=first, dcp=-3.14) ////////////////" << endl;
  
  //fitmarg->SetMinimIntervalRead(valsecondbin_X+conv::kASmallNum,0.4,-TMath::Pi(),+TMath::Pi(),true);
  fitmarg->SetMinimIntervalRead(valsecondbin_X+conv::kASmallNum,vallastbin_X,valfirstbin_Y,vallastbin_Y,true);
  fitmarg->SetStartBin(-TMath::Pi(),start_sinsqth13_first);
  fitmarg->SetNumStartBin(1,1);
  //fitmarg->SetGridStartRange(valsecondbin_X+conv::kASmallNum,0.4,-TMath::Pi(),-TMath::Pi()+conv::kASmallNum);
  fitmarg->SetGridStartRange(valsecondbin_X+conv::kASmallNum,vallastbin_X,valfirstbin_Y,vallastbin_Y+conv::kASmallNum);
  fitmarg->SetExtractStartPt("Single");
  fitmarg->MinimizeChisqStepRead(skobs_fakedata,true);
  
  double chisqmin_2 = (double)fitmarg->GetChisqBestFit();
  double bestfit0_2 = fitmarg->GetParamBestFit(0);
  double bestfit1_2 = -999.;
  if(fitmarg->IsFit2D()) bestfit1_2 = fitmarg->GetParamBestFit(1);
  int nfits_2 = fitmarg->GetNumFit();
  int fitstatus_2 = fitmarg->GetFitStatus();
  int niterations_2 = fitmarg->GetNumIterations();
  if(chisqmin_2<chisqmin){
    chisqmin = chisqmin_2;
    bestfit0 = bestfit0_2;
    bestfit1 = bestfit1_2;
    nfits = nfits_2;
    fitstatus = fitstatus_2;
    niterations = niterations_2;
  }
  
  // Fit 3
  
  cout << endl;
  cout << "/////////////// Fit 3 (mixing13=start_sinsqth13_first, dcp=-pi/2.) ////////////////" << endl;
  
  //fitmarg->SetMinimIntervalRead(valsecondbin_X+conv::kASmallNum,0.4,-TMath::Pi(),+TMath::Pi(),true);
  fitmarg->SetMinimIntervalRead(valsecondbin_X+conv::kASmallNum,vallastbin_X,valfirstbin_Y,vallastbin_Y,true);
  fitmarg->SetStartBin(-TMath::Pi()/2.,start_sinsqth13_first);
  fitmarg->SetNumStartBin(1,1);
  //fitmarg->SetGridStartRange(valsecondbin_X+conv::kASmallNum,0.4,-TMath::Pi(),-TMath::Pi()+conv::kASmallNum);
  fitmarg->SetGridStartRange(valsecondbin_X+conv::kASmallNum,vallastbin_X,valfirstbin_Y,vallastbin_Y+conv::kASmallNum);
  fitmarg->SetExtractStartPt("Single");
  fitmarg->MinimizeChisqStepRead(skobs_fakedata,true);
  
  double chisqmin_3 = (double)fitmarg->GetChisqBestFit();
  double bestfit0_3 = fitmarg->GetParamBestFit(0);
  double bestfit1_3 = -999.;
  if(fitmarg->IsFit2D()) bestfit1_3 = fitmarg->GetParamBestFit(1);
  int nfits_3 = fitmarg->GetNumFit();
  int fitstatus_3 = fitmarg->GetFitStatus();
  int niterations_3 = fitmarg->GetNumIterations();
  if(chisqmin_3<chisqmin){
    chisqmin = chisqmin_3;
    bestfit0 = bestfit0_3;
    bestfit1 = bestfit1_3;
    nfits = nfits_3;
    fitstatus = fitstatus_3;
    niterations = niterations_3;
  }
  
  
  // Fit 4
  
  cout << endl;
  cout << "/////////////// Fit 4 (mixing13=first, dcp=+pi/2.) ////////////////" << endl;
  
  //fitmarg->SetMinimIntervalRead(valsecondbin_X+conv::kASmallNum,0.4,-TMath::Pi(),+TMath::Pi(),true);
  fitmarg->SetMinimIntervalRead(valsecondbin_X+conv::kASmallNum,vallastbin_X,valfirstbin_Y,vallastbin_Y,true);
  fitmarg->SetStartBin(+TMath::Pi()/2., start_sinsqth13_first);
  fitmarg->SetNumStartBin(1,1);
  //fitmarg->SetGridStartRange(valsecondbin_X+conv::kASmallNum,0.4,-TMath::Pi(),-TMath::Pi()+conv::kASmallNum);
  fitmarg->SetGridStartRange(valsecondbin_X+conv::kASmallNum,vallastbin_X,valfirstbin_Y,vallastbin_Y+conv::kASmallNum);
  fitmarg->SetExtractStartPt("Single");
  fitmarg->MinimizeChisqStepRead(skobs_fakedata,true);
  
  double chisqmin_4 = (double)fitmarg->GetChisqBestFit();
  double bestfit0_4 = fitmarg->GetParamBestFit(0);
  double bestfit1_4 = -999.;
  if(fitmarg->IsFit2D()) bestfit1_4 = fitmarg->GetParamBestFit(1);
  int nfits_4 = fitmarg->GetNumFit();
  int fitstatus_4 = fitmarg->GetFitStatus();
  int niterations_4 = fitmarg->GetNumIterations();
  if(chisqmin_4<chisqmin){
    chisqmin = chisqmin_4;
    bestfit0 = bestfit0_4;
    bestfit1 = bestfit1_4;
    nfits = nfits_4;
    fitstatus = fitstatus_4;
    niterations = niterations_4;
  }
  
  
  // Fit 5
  
  cout << endl;
  cout << "/////////////// Fit 5 (mixing13=first, dcp=+3.14) ////////////////" << endl;
  
  //fitmarg->SetMinimIntervalRead(valsecondbin_X+conv::kASmallNum,0.4,-TMath::Pi(),+TMath::Pi(),true);
  fitmarg->SetMinimIntervalRead(valsecondbin_X+conv::kASmallNum,vallastbin_X,valfirstbin_Y,vallastbin_Y,true);
  fitmarg->SetStartBin(+TMath::Pi(),start_sinsqth13_first);
  fitmarg->SetNumStartBin(1,1);
  //fitmarg->SetGridStartRange(valsecondbin_X+conv::kASmallNum,0.4,-TMath::Pi(),-TMath::Pi()+conv::kASmallNum);
  fitmarg->SetGridStartRange(valsecondbin_X+conv::kASmallNum,vallastbin_X,valfirstbin_Y,vallastbin_Y+conv::kASmallNum);
  fitmarg->SetExtractStartPt("Single");
  fitmarg->MinimizeChisqStepRead(skobs_fakedata,true);
  
  double chisqmin_5 = (double)fitmarg->GetChisqBestFit();
  double bestfit0_5 = fitmarg->GetParamBestFit(0);
  double bestfit1_5 = -999.;
  if(fitmarg->IsFit2D()) bestfit1_5 = fitmarg->GetParamBestFit(1);
  int nfits_5 = fitmarg->GetNumFit();
  int fitstatus_5 = fitmarg->GetFitStatus();
  int niterations_5 = fitmarg->GetNumIterations();
  if(chisqmin_5<chisqmin){
    chisqmin = chisqmin_5;
    bestfit0 = bestfit0_5;
    bestfit1 = bestfit1_5;
    nfits = nfits_5;
    fitstatus = fitstatus_5;
    niterations = niterations_5;
  }
  
  
  // Fit 6
  
  cout << endl;
  cout << "/////////////// Fit 6 (mixing13=last, dcp=-3.14) ////////////////" << endl;
  
  //fitmarg->SetMinimIntervalRead(valsecondbin_X+conv::kASmallNum,0.4,-TMath::Pi(),+TMath::Pi(),true);
  fitmarg->SetMinimIntervalRead(valsecondbin_X+conv::kASmallNum,vallastbin_X,valfirstbin_Y,vallastbin_Y,true);
  fitmarg->SetStartBin(-TMath::Pi(),start_sinsqth13_last);
  fitmarg->SetNumStartBin(1,1);
  //fitmarg->SetGridStartRange(valsecondbin_X+conv::kASmallNum,0.4,-TMath::Pi(),-TMath::Pi()+conv::kASmallNum);
  fitmarg->SetGridStartRange(valsecondbin_X+conv::kASmallNum,vallastbin_X,valfirstbin_Y,vallastbin_Y+conv::kASmallNum);
  fitmarg->SetExtractStartPt("Single");
  fitmarg->MinimizeChisqStepRead(skobs_fakedata,true);
  
  double chisqmin_6 = (double)fitmarg->GetChisqBestFit();
  double bestfit0_6 = fitmarg->GetParamBestFit(0);
  double bestfit1_6 = -999.;
  if(fitmarg->IsFit2D()) bestfit1_6 = fitmarg->GetParamBestFit(1);
  int nfits_6 = fitmarg->GetNumFit();
  int fitstatus_6 = fitmarg->GetFitStatus();
  int niterations_6 = fitmarg->GetNumIterations();
  if(chisqmin_6<chisqmin){
    chisqmin = chisqmin_6;
    bestfit0 = bestfit0_6;
    bestfit1 = bestfit1_6;
    nfits = nfits_6;
    fitstatus = fitstatus_6;
    niterations = niterations_6;
  }
  
  // Fit 7
  
  cout << endl;
  cout << "/////////////// Fit 3 (mixing13=last, dcp=-pi/2.) ////////////////" << endl;
  
  //fitmarg->SetMinimIntervalRead(valsecondbin_X+conv::kASmallNum,0.4,-TMath::Pi(),+TMath::Pi(),true);
  fitmarg->SetMinimIntervalRead(valsecondbin_X+conv::kASmallNum,vallastbin_X,valfirstbin_Y,vallastbin_Y,true);
  fitmarg->SetStartBin(-TMath::Pi()/2.,start_sinsqth13_last);
  fitmarg->SetNumStartBin(1,1);
  //fitmarg->SetGridStartRange(valsecondbin_X+conv::kASmallNum,0.4,-TMath::Pi(),-TMath::Pi()+conv::kASmallNum);
  fitmarg->SetGridStartRange(valsecondbin_X+conv::kASmallNum,vallastbin_X,valfirstbin_Y,vallastbin_Y+conv::kASmallNum);
  fitmarg->SetExtractStartPt("Single");
  fitmarg->MinimizeChisqStepRead(skobs_fakedata,true);
  
  double chisqmin_7 = (double)fitmarg->GetChisqBestFit();
  double bestfit0_7 = fitmarg->GetParamBestFit(0);
  double bestfit1_7 = -999.;
  if(fitmarg->IsFit2D()) bestfit1_7 = fitmarg->GetParamBestFit(1);
  int nfits_7 = fitmarg->GetNumFit();
  int fitstatus_7 = fitmarg->GetFitStatus();
  int niterations_7 = fitmarg->GetNumIterations();
  if(chisqmin_7<chisqmin){
    chisqmin = chisqmin_7;
    bestfit0 = bestfit0_7;
    bestfit1 = bestfit1_7;
    nfits = nfits_7;
    fitstatus = fitstatus_7;
    niterations = niterations_7;
  }
  
  
  // Fit 8
  
  cout << endl;
  cout << "/////////////// Fit 8 (mixing13=last, dcp=+pi/2.) ////////////////" << endl;
  
  //fitmarg->SetMinimIntervalRead(valsecondbin_X+conv::kASmallNum,0.4,-TMath::Pi(),+TMath::Pi(),true);
  fitmarg->SetMinimIntervalRead(valsecondbin_X+conv::kASmallNum,vallastbin_X,valfirstbin_Y,vallastbin_Y,true);
  fitmarg->SetStartBin(+TMath::Pi()/2., start_sinsqth13_last);
  fitmarg->SetNumStartBin(1,1);
  //fitmarg->SetGridStartRange(valsecondbin_X+conv::kASmallNum,0.4,-TMath::Pi(),-TMath::Pi()+conv::kASmallNum);
  fitmarg->SetGridStartRange(valsecondbin_X+conv::kASmallNum,vallastbin_X,valfirstbin_Y,vallastbin_Y+conv::kASmallNum);
  fitmarg->SetExtractStartPt("Single");
  fitmarg->MinimizeChisqStepRead(skobs_fakedata,true);
  
  double chisqmin_8 = (double)fitmarg->GetChisqBestFit();
  double bestfit0_8 = fitmarg->GetParamBestFit(0);
  double bestfit1_8 = -999.;
  if(fitmarg->IsFit2D()) bestfit1_8 = fitmarg->GetParamBestFit(1);
  int nfits_8 = fitmarg->GetNumFit();
  int fitstatus_8 = fitmarg->GetFitStatus();
  int niterations_8 = fitmarg->GetNumIterations();
  if(chisqmin_8<chisqmin){
    chisqmin = chisqmin_8;
    bestfit0 = bestfit0_8;
    bestfit1 = bestfit1_8;
    nfits = nfits_8;
    fitstatus = fitstatus_8;
    niterations = niterations_8;
  }
  
  
  // Fit 9
  
  cout << endl;
  cout << "/////////////// Fit 9 (mixing13=last, dcp=+3.14) ////////////////" << endl;
  
  //fitmarg->SetMinimIntervalRead(valsecondbin_X+conv::kASmallNum,0.4,-TMath::Pi(),+TMath::Pi(),true);
  fitmarg->SetMinimIntervalRead(valsecondbin_X+conv::kASmallNum,vallastbin_X,valfirstbin_Y,vallastbin_Y,true);
  fitmarg->SetStartBin(+TMath::Pi(), start_sinsqth13_last);
  fitmarg->SetNumStartBin(1,1);
  //fitmarg->SetGridStartRange(valsecondbin_X+conv::kASmallNum,0.4,-TMath::Pi(),-TMath::Pi()+conv::kASmallNum);
  fitmarg->SetGridStartRange(valsecondbin_X+conv::kASmallNum,vallastbin_X,valfirstbin_Y,vallastbin_Y+conv::kASmallNum);
  fitmarg->SetExtractStartPt("Single");
  fitmarg->MinimizeChisqStepRead(skobs_fakedata,true);
  
  double chisqmin_9 = (double)fitmarg->GetChisqBestFit();
  double bestfit0_9 = fitmarg->GetParamBestFit(0);
  double bestfit1_9 = -999.;
  if(fitmarg->IsFit2D()) bestfit1_9 = fitmarg->GetParamBestFit(1);
  int nfits_9 = fitmarg->GetNumFit();
  int fitstatus_9 = fitmarg->GetFitStatus();
  int niterations_9 = fitmarg->GetNumIterations();
  if(chisqmin_9<chisqmin){
    chisqmin = chisqmin_9;
    bestfit0 = bestfit0_9;
    bestfit1 = bestfit1_9;
    nfits = nfits_9;
    fitstatus = fitstatus_9;
    niterations = niterations_9;
  }
  
  
  
  
  
  
  
  // // Fit 2
  // cout << endl;
  // cout << "/////////////// Fit 2 (mixing13>0, dcp=-3.14) ////////////////" << endl;
  // fitmarg->SetMinimIntervalRead(valsecondbin_X+conv::kASmallNum,0.3,-TMath::Pi(),+TMath::Pi(),true);
  // fitmarg->SetStartBin(0.1,-TMath::Pi());
  // fitmarg->SetNumStartBin(10,1);
  // fitmarg->SetGridStartRange(valsecondbin_X+conv::kASmallNum,0.3,-TMath::Pi(),-TMath::Pi()+conv::kASmallNum);
  // fitmarg->SetExtractStartPt("GridFix");
  // fitmarg->MinimizeChisqStepRead(skobs_fakedata,true);
  // double chisqmin_2 = (double)fitmarg->GetChisqBestFit();
  // double bestfit0_2 = fitmarg->GetParamBestFit(0);
  // double bestfit1_2 = -999.;
  // if(fitmarg->IsFit2D()) bestfit1_2 = fitmarg->GetParamBestFit(1);
  // int nfits_2 = fitmarg->GetNumFit();
  // int fitstatus_2 = fitmarg->GetFitStatus();
  // int niterations_2 = fitmarg->GetNumIterations();
  // if(chisqmin_2<chisqmin){
  //   chisqmin = chisqmin_2;
  //   bestfit0 = bestfit0_2;
  //   bestfit1 = bestfit1_2;
  //   nfits = nfits_2;
  //   fitstatus = fitstatus_2;
  //   niterations = niterations_2;
  // }
  
  // // Fit 3
  // cout << endl;
  // cout << "/////////////// Fit 3 (mixing13>0, start dcp = 0.) ////////////////" << endl;
  // int bindcp0_Y = fitmarg->GetMapChisqRangeFileID_2D()->GetYaxis()->FindBin(0.);
  // double valdcp0_Y = fitmarg->GetMapChisqRangeFileID_2D()->GetYaxis()->GetBinCenter(bindcp0_Y);
  // double binwidth = fitmarg->GetMapChisqRangeFileID_2D()->GetYaxis()->GetBinWidth(1);
  // fitmarg->SetMinimIntervalRead(valsecondbin_X+conv::kASmallNum,0.3,-TMath::Pi(),+TMath::Pi(),true);
  // fitmarg->SetStartBin(0.1,-TMath::Pi());
  // fitmarg->SetNumStartBin(10,1);
  // fitmarg->SetGridStartRange(valsecondbin_X+conv::kASmallNum,0.3,valdcp0_Y-binwidth,valdcp0_Y+binwidth);
  // fitmarg->SetExtractStartPt("GridFix");
  // fitmarg->MinimizeChisqStepRead(skobs_fakedata,true);
  // double chisqmin_3 = (double)fitmarg->GetChisqBestFit();
  // double bestfit0_3 = fitmarg->GetParamBestFit(0);
  // double bestfit1_3 = -999.;
  // if(fitmarg->IsFit2D()) bestfit1_3 = fitmarg->GetParamBestFit(1);
  // int nfits_3 = fitmarg->GetNumFit();
  // int fitstatus_3 = fitmarg->GetFitStatus();
  // int niterations_3 = fitmarg->GetNumIterations();
  // if(chisqmin_3<chisqmin){
  //   chisqmin = chisqmin_3;
  //   bestfit0 = bestfit0_3;
  //   bestfit1 = bestfit1_3;
  //   nfits = nfits_3;
  //   fitstatus = fitstatus_3;
  //   niterations = niterations_3;
  // }
  
  // // Fit 4
  // cout << endl;
  // cout << "/////////////// Fit 4 (mixing13>0, start dcp=+3.14) ////////////////" << endl;
  // fitmarg->SetMinimIntervalRead(valsecondbin_X+conv::kASmallNum,0.3,-TMath::Pi(),+TMath::Pi(),true);
  // fitmarg->SetStartBin(0.1,-TMath::Pi());
  // fitmarg->SetNumStartBin(10,1);
  // fitmarg->SetGridStartRange(valsecondbin_X+conv::kASmallNum,0.3,vallastbin_Y-conv::kASmallNum,+TMath::Pi());
  // fitmarg->SetExtractStartPt("GridFix");
  // fitmarg->MinimizeChisqStepRead(skobs_fakedata,true);
  // double chisqmin_4 = (double)fitmarg->GetChisqBestFit();
  // double bestfit0_4 = fitmarg->GetParamBestFit(0);
  // double bestfit1_4 = -999.;
  // if(fitmarg->IsFit2D()) bestfit1_4 = fitmarg->GetParamBestFit(1);
  // int nfits_4 = fitmarg->GetNumFit();
  // int fitstatus_4 = fitmarg->GetFitStatus();
  // int niterations_4 = fitmarg->GetNumIterations();
  // if(chisqmin_4<chisqmin){
  //   chisqmin = chisqmin_4;
  //   bestfit0 = bestfit0_4;
  //   bestfit1 = bestfit1_4;
  //   nfits = nfits_4;
  //   fitstatus = fitstatus_4;
  //   niterations = niterations_4;
  // }
}
//____________________________________________________________________________
void valor::util::MargFit_2D_sinsqth23_deltaCP(RunOpt *runopt,FitMarginalization *fitmarg,SuperkObservation *skobs_fakedata,double &chisqmin,double &bestfit0, double &bestfit1,int nfits,int fitstatus,int niterations)
{
  
  double valsecondbin_X = fitmarg->GetMapChisqRangeFileID_2D()->GetXaxis()->GetBinCenter(2);
  double valsecondbin_Y = fitmarg->GetMapChisqRangeFileID_2D()->GetYaxis()->GetBinCenter(2);

  int nbinsX = fitmarg->GetMapChisqRangeFileID_2D()->GetNbinsX();
  double valfirstbin_X = fitmarg->GetMapChisqRangeFileID_2D()->GetXaxis()->GetBinCenter(1);
  double vallastbin_X = fitmarg->GetMapChisqRangeFileID_2D()->GetXaxis()->GetBinCenter(nbinsX);
  
  int nbinsY = fitmarg->GetMapChisqRangeFileID_2D()->GetNbinsY();
  double valfirstbin_Y = fitmarg->GetMapChisqRangeFileID_2D()->GetYaxis()->GetBinCenter(1);
  double vallastbin_Y = fitmarg->GetMapChisqRangeFileID_2D()->GetYaxis()->GetBinCenter(nbinsY);
    
  double start_sinsqth23_first = valfirstbin_X;
  double start_sinsqth23_last  = vallastbin_X;
  
  // Fit 1
   
  cout << endl;
  cout << "/////////////// Fit 1 (mixing23=first_bin, dcp=-3.14) ////////////////" << endl;
  
  fitmarg->SetMinimIntervalRead(valfirstbin_X+conv::kASmallNum,vallastbin_X,valfirstbin_Y,vallastbin_Y,true);
  fitmarg->SetStartBin(start_sinsqth23_first,-TMath::Pi());
  fitmarg->SetNumStartBin(1,1);
  fitmarg->SetGridStartRange(valfirstbin_X+conv::kASmallNum,vallastbin_X,valfirstbin_Y,vallastbin_Y+conv::kASmallNum);
  fitmarg->SetExtractStartPt("Single");
  fitmarg->MinimizeChisqStepRead(skobs_fakedata,true);      

  double chisqmin_1 = (double)fitmarg->GetChisqBestFit();
  double bestfit0_1 = fitmarg->GetParamBestFit(0);
  double bestfit1_1 = -999.;
  if(fitmarg->IsFit2D()) bestfit1_1 = fitmarg->GetParamBestFit(1); 
  int nfits_1 = fitmarg->GetNumFit();
  int fitstatus_1 = fitmarg->GetFitStatus();
  int niterations_1 = fitmarg->GetNumIterations();
  if(chisqmin_1<chisqmin){
    chisqmin = chisqmin_1;
    bestfit0 = bestfit0_1; 
    bestfit1 = bestfit1_1; 
    nfits = nfits_1;    
    fitstatus = fitstatus_1;
    niterations = niterations_1;
  }
    
  // Fit 2
  
  cout << endl;
  cout << "/////////////// Fit 2 (mixing23=first_bin, dcp=-pi/2.) ////////////////" << endl;
  
  fitmarg->SetMinimIntervalRead(valfirstbin_X+conv::kASmallNum,vallastbin_X,valfirstbin_Y,vallastbin_Y,true);
  fitmarg->SetStartBin(start_sinsqth23_first,-TMath::Pi()/2.); 
  fitmarg->SetNumStartBin(1,1);
  fitmarg->SetGridStartRange(valfirstbin_X+conv::kASmallNum,vallastbin_X,valfirstbin_Y,vallastbin_Y+conv::kASmallNum);
  fitmarg->SetExtractStartPt("Single");
  fitmarg->MinimizeChisqStepRead(skobs_fakedata,true);      

  double chisqmin_2 = (double)fitmarg->GetChisqBestFit();
  double bestfit0_2 = fitmarg->GetParamBestFit(0);
  double bestfit1_2 = -999.;
  if(fitmarg->IsFit2D()) bestfit1_2 = fitmarg->GetParamBestFit(1); 
  int nfits_2 = fitmarg->GetNumFit();
  int fitstatus_2 = fitmarg->GetFitStatus();
  int niterations_2 = fitmarg->GetNumIterations();
  if(chisqmin_2<chisqmin){
    chisqmin = chisqmin_2;
    bestfit0 = bestfit0_2; 
    bestfit1 = bestfit1_2; 
    nfits = nfits_2;    
    fitstatus = fitstatus_2;
    niterations = niterations_2;
  }
  
  
  // Fit 3
  
  cout << endl;
  cout << "/////////////// Fit 3 (mixing23=first, dcp=+pi/2.) ////////////////" << endl;
  
  fitmarg->SetMinimIntervalRead(valfirstbin_X+conv::kASmallNum,vallastbin_X,valfirstbin_Y,vallastbin_Y,true);
  fitmarg->SetStartBin(start_sinsqth23_first,+TMath::Pi()/2.); 
  fitmarg->SetNumStartBin(1,1);
  fitmarg->SetGridStartRange(valfirstbin_X+conv::kASmallNum,vallastbin_X,valfirstbin_Y,vallastbin_Y+conv::kASmallNum);
  fitmarg->SetExtractStartPt("Single");
  fitmarg->MinimizeChisqStepRead(skobs_fakedata,true);      
  
  double chisqmin_3 = (double)fitmarg->GetChisqBestFit();
  double bestfit0_3 = fitmarg->GetParamBestFit(0);
  double bestfit1_3 = -999.;
  if(fitmarg->IsFit2D()) bestfit1_3 = fitmarg->GetParamBestFit(1); 
  int nfits_3 = fitmarg->GetNumFit();
  int fitstatus_3 = fitmarg->GetFitStatus();
  int niterations_3 = fitmarg->GetNumIterations();
  if(chisqmin_3<chisqmin){
    chisqmin = chisqmin_3;
    bestfit0 = bestfit0_3; 
    bestfit1 = bestfit1_3; 
    nfits = nfits_3;    
    fitstatus = fitstatus_3;
    niterations = niterations_3;
  }
  
  
  // Fit 4
  
  cout << endl;
  cout << "/////////////// Fit 4 (mixing23=first, dcp=+3.14) ////////////////" << endl;
  
  fitmarg->SetMinimIntervalRead(valfirstbin_X+conv::kASmallNum,vallastbin_X,valfirstbin_Y,vallastbin_Y,true);
  fitmarg->SetStartBin(start_sinsqth23_first,+TMath::Pi()); 
  fitmarg->SetNumStartBin(1,1);
  fitmarg->SetGridStartRange(valfirstbin_X+conv::kASmallNum,vallastbin_X,valfirstbin_Y,vallastbin_Y+conv::kASmallNum);
  fitmarg->SetExtractStartPt("Single");
  fitmarg->MinimizeChisqStepRead(skobs_fakedata,true);      
  
  double chisqmin_4 = (double)fitmarg->GetChisqBestFit();
  double bestfit0_4 = fitmarg->GetParamBestFit(0);
  double bestfit1_4 = -999.;
  if(fitmarg->IsFit2D()) bestfit1_4 = fitmarg->GetParamBestFit(1); 
  int nfits_4 = fitmarg->GetNumFit();
  int fitstatus_4 = fitmarg->GetFitStatus();
  int niterations_4 = fitmarg->GetNumIterations();
  if(chisqmin_4<chisqmin){
    chisqmin = chisqmin_4;
    bestfit0 = bestfit0_4; 
    bestfit1 = bestfit1_4; 
    nfits = nfits_4;    
    fitstatus = fitstatus_4;
    niterations = niterations_4;
  }
  

  // Fit 5
  
  cout << endl;
  cout << "/////////////// Fit 5 (mixing23=last_bin, dcp=-3.14) ////////////////" << endl;
  
  fitmarg->SetMinimIntervalRead(valfirstbin_X+conv::kASmallNum,vallastbin_X,valfirstbin_Y,vallastbin_Y,true);
  fitmarg->SetStartBin(start_sinsqth23_last,-TMath::Pi()); 
  fitmarg->SetNumStartBin(1,1);
  fitmarg->SetGridStartRange(valfirstbin_X+conv::kASmallNum,vallastbin_X,valfirstbin_Y,vallastbin_Y+conv::kASmallNum);
  fitmarg->SetExtractStartPt("Single");
  fitmarg->MinimizeChisqStepRead(skobs_fakedata,true);      
  
  double chisqmin_5 = (double)fitmarg->GetChisqBestFit();
  double bestfit0_5 = fitmarg->GetParamBestFit(0);
  double bestfit1_5 = -999.;
  if(fitmarg->IsFit2D()) bestfit1_5 = fitmarg->GetParamBestFit(1); 
  int nfits_5 = fitmarg->GetNumFit();
  int fitstatus_5 = fitmarg->GetFitStatus();
  int niterations_5 = fitmarg->GetNumIterations();
  if(chisqmin_5<chisqmin){
    chisqmin = chisqmin_5;
    bestfit0 = bestfit0_5; 
    bestfit1 = bestfit1_5; 
    nfits = nfits_5;    
    fitstatus = fitstatus_5;
    niterations = niterations_5;
  }
    
  // Fit 6
  
  cout << endl;
  cout << "/////////////// Fit 6 (mixing23=last_bin, dcp=-pi/2.) ////////////////" << endl;
  
  fitmarg->SetMinimIntervalRead(valfirstbin_X+conv::kASmallNum,vallastbin_X,valfirstbin_Y,vallastbin_Y,true);
  fitmarg->SetStartBin(start_sinsqth23_last,-TMath::Pi()/2.); 
  fitmarg->SetNumStartBin(1,1);
  fitmarg->SetGridStartRange(valfirstbin_X+conv::kASmallNum,vallastbin_X,valfirstbin_Y,vallastbin_Y+conv::kASmallNum);
  fitmarg->SetExtractStartPt("Single");
  fitmarg->MinimizeChisqStepRead(skobs_fakedata,true);      
  
  double chisqmin_6 = (double)fitmarg->GetChisqBestFit();
  double bestfit0_6 = fitmarg->GetParamBestFit(0);
  double bestfit1_6 = -999.;
  if(fitmarg->IsFit2D()) bestfit1_6 = fitmarg->GetParamBestFit(1); 
  int nfits_6 = fitmarg->GetNumFit();
  int fitstatus_6 = fitmarg->GetFitStatus();
  int niterations_6 = fitmarg->GetNumIterations();
  if(chisqmin_6<chisqmin){
    chisqmin = chisqmin_6;
    bestfit0 = bestfit0_6; 
    bestfit1 = bestfit1_6; 
    nfits = nfits_6;    
    fitstatus = fitstatus_6;
    niterations = niterations_6;
  }
  
  
  // Fit 7
  
  cout << endl;
  cout << "/////////////// Fit 7 (mixing23=last_bin, dcp=+pi/2.) ////////////////" << endl;
  
  fitmarg->SetMinimIntervalRead(valfirstbin_X+conv::kASmallNum,vallastbin_X,valfirstbin_Y,vallastbin_Y,true);
  fitmarg->SetStartBin(start_sinsqth23_last,+TMath::Pi()/2.); 
  fitmarg->SetNumStartBin(1,1);
  fitmarg->SetGridStartRange(valfirstbin_X+conv::kASmallNum,vallastbin_X,valfirstbin_Y,vallastbin_Y+conv::kASmallNum);
  fitmarg->SetExtractStartPt("Single");
  fitmarg->MinimizeChisqStepRead(skobs_fakedata,true);      
  
  double chisqmin_7 = (double)fitmarg->GetChisqBestFit();
  double bestfit0_7 = fitmarg->GetParamBestFit(0);
  double bestfit1_7 = -999.;
  if(fitmarg->IsFit2D()) bestfit1_7 = fitmarg->GetParamBestFit(1); 
  int nfits_7 = fitmarg->GetNumFit();
  int fitstatus_7 = fitmarg->GetFitStatus();
  int niterations_7 = fitmarg->GetNumIterations();
  if(chisqmin_7<chisqmin){
    chisqmin = chisqmin_7;
    bestfit0 = bestfit0_7; 
    bestfit1 = bestfit1_7; 
    nfits = nfits_7;    
    fitstatus = fitstatus_7;
    niterations = niterations_7;
  }
  
  
  // Fit 8
  
  cout << endl;
  cout << "/////////////// Fit 8 (mixing23=last_bin, dcp=+3.14) ////////////////" << endl;
  
  fitmarg->SetMinimIntervalRead(valfirstbin_X+conv::kASmallNum,vallastbin_X,valfirstbin_Y,vallastbin_Y,true);
  fitmarg->SetStartBin(start_sinsqth23_last,+TMath::Pi()); 
  fitmarg->SetNumStartBin(1,1);
  fitmarg->SetGridStartRange(valfirstbin_X+conv::kASmallNum,vallastbin_X,valfirstbin_Y,vallastbin_Y+conv::kASmallNum);
  fitmarg->SetExtractStartPt("Single");
  fitmarg->MinimizeChisqStepRead(skobs_fakedata,true);      
  
  double chisqmin_8 = (double)fitmarg->GetChisqBestFit();
  double bestfit0_8 = fitmarg->GetParamBestFit(0);
  double bestfit1_8 = -999.;
  if(fitmarg->IsFit2D()) bestfit1_8 = fitmarg->GetParamBestFit(1); 
  int nfits_8 = fitmarg->GetNumFit();
  int fitstatus_8 = fitmarg->GetFitStatus();
  int niterations_8 = fitmarg->GetNumIterations();
  if(chisqmin_8<chisqmin){
    chisqmin = chisqmin_8;
    bestfit0 = bestfit0_8; 
    bestfit1 = bestfit1_8; 
    nfits = nfits_8;    
    fitstatus = fitstatus_8;
    niterations = niterations_8;
  }
}
//____________________________________________________________________________


















void valor::util::MargFit_2D_dmsq32_deltaCP(RunOpt *runopt,FitMarginalization *fitmarg,SuperkObservation *skobs_fakedata,double &chisqmin,double &bestfit0, double &bestfit1,int nfits,int fitstatus,int niterations)
{
  
  double valsecondbin_X = fitmarg->GetMapChisqRangeFileID_2D()->GetXaxis()->GetBinCenter(2);
  double valsecondbin_Y = fitmarg->GetMapChisqRangeFileID_2D()->GetYaxis()->GetBinCenter(2);

  int nbinsX = fitmarg->GetMapChisqRangeFileID_2D()->GetNbinsX();
  double valfirstbin_X = fitmarg->GetMapChisqRangeFileID_2D()->GetXaxis()->GetBinCenter(1);
  double vallastbin_X = fitmarg->GetMapChisqRangeFileID_2D()->GetXaxis()->GetBinCenter(nbinsX);
  
  int nbinsY = fitmarg->GetMapChisqRangeFileID_2D()->GetNbinsY();
  double valfirstbin_Y = fitmarg->GetMapChisqRangeFileID_2D()->GetYaxis()->GetBinCenter(1);
  double vallastbin_Y = fitmarg->GetMapChisqRangeFileID_2D()->GetYaxis()->GetBinCenter(nbinsY);
    
  double start_dmsq32_first = valfirstbin_X;
  double start_dmsq32_last  = vallastbin_X;
  
  // Fit 1
   
  cout << endl;
  cout << "/////////////// Fit 1 (dmsq32=first_bin, dcp=-3.14) ////////////////" << endl;
  
  fitmarg->SetMinimIntervalRead(valfirstbin_X+conv::kASmallNum,vallastbin_X,valfirstbin_Y,vallastbin_Y,true);
  fitmarg->SetStartBin(start_dmsq32_first,-TMath::Pi()); 
  fitmarg->SetNumStartBin(1,1);
  fitmarg->SetGridStartRange(valfirstbin_X+conv::kASmallNum,vallastbin_X,valfirstbin_Y,vallastbin_Y+conv::kASmallNum);
  fitmarg->SetExtractStartPt("Single");
  fitmarg->MinimizeChisqStepRead(skobs_fakedata,true);      

  double chisqmin_1 = (double)fitmarg->GetChisqBestFit();
  double bestfit0_1 = fitmarg->GetParamBestFit(0);
  double bestfit1_1 = -999.;
  if(fitmarg->IsFit2D()) bestfit1_1 = fitmarg->GetParamBestFit(1); 
  int nfits_1 = fitmarg->GetNumFit();
  int fitstatus_1 = fitmarg->GetFitStatus();
  int niterations_1 = fitmarg->GetNumIterations();
  if(chisqmin_1<chisqmin){
    chisqmin = chisqmin_1;
    bestfit0 = bestfit0_1; 
    bestfit1 = bestfit1_1; 
    nfits = nfits_1;    
    fitstatus = fitstatus_1;
    niterations = niterations_1;
  }
    
  // Fit 2
  
  cout << endl;
  cout << "/////////////// Fit 2 (dmsq32=first_bin, dcp=-pi/2.) ////////////////" << endl;
  
  fitmarg->SetMinimIntervalRead(valfirstbin_X+conv::kASmallNum,vallastbin_X,valfirstbin_Y,vallastbin_Y,true);
  fitmarg->SetStartBin(start_dmsq32_first,-TMath::Pi()/2.); 
  fitmarg->SetNumStartBin(1,1);
  fitmarg->SetGridStartRange(valfirstbin_X+conv::kASmallNum,vallastbin_X,valfirstbin_Y,vallastbin_Y+conv::kASmallNum);
  fitmarg->SetExtractStartPt("Single");
  fitmarg->MinimizeChisqStepRead(skobs_fakedata,true);      

  double chisqmin_2 = (double)fitmarg->GetChisqBestFit();
  double bestfit0_2 = fitmarg->GetParamBestFit(0);
  double bestfit1_2 = -999.;
  if(fitmarg->IsFit2D()) bestfit1_2 = fitmarg->GetParamBestFit(1); 
  int nfits_2 = fitmarg->GetNumFit();
  int fitstatus_2 = fitmarg->GetFitStatus();
  int niterations_2 = fitmarg->GetNumIterations();
  if(chisqmin_2<chisqmin){
    chisqmin = chisqmin_2;
    bestfit0 = bestfit0_2; 
    bestfit1 = bestfit1_2; 
    nfits = nfits_2;    
    fitstatus = fitstatus_2;
    niterations = niterations_2;
  }
  
  
  // Fit 3
  
  cout << endl;
  cout << "/////////////// Fit 3 (dmsq32=first, dcp=+pi/2.) ////////////////" << endl;
  
  fitmarg->SetMinimIntervalRead(valfirstbin_X+conv::kASmallNum,vallastbin_X,valfirstbin_Y,vallastbin_Y,true);
  fitmarg->SetStartBin(start_dmsq32_first,+TMath::Pi()/2.); 
  fitmarg->SetNumStartBin(1,1);
  fitmarg->SetGridStartRange(valfirstbin_X+conv::kASmallNum,vallastbin_X,valfirstbin_Y,vallastbin_Y+conv::kASmallNum);
  fitmarg->SetExtractStartPt("Single");
  fitmarg->MinimizeChisqStepRead(skobs_fakedata,true);      
  
  double chisqmin_3 = (double)fitmarg->GetChisqBestFit();
  double bestfit0_3 = fitmarg->GetParamBestFit(0);
  double bestfit1_3 = -999.;
  if(fitmarg->IsFit2D()) bestfit1_3 = fitmarg->GetParamBestFit(1); 
  int nfits_3 = fitmarg->GetNumFit();
  int fitstatus_3 = fitmarg->GetFitStatus();
  int niterations_3 = fitmarg->GetNumIterations();
  if(chisqmin_3<chisqmin){
    chisqmin = chisqmin_3;
    bestfit0 = bestfit0_3; 
    bestfit1 = bestfit1_3; 
    nfits = nfits_3;    
    fitstatus = fitstatus_3;
    niterations = niterations_3;
  }
  
  
  // Fit 4
  
  cout << endl;
  cout << "/////////////// Fit 4 (dmsq32=first, dcp=+3.14) ////////////////" << endl;
  
  fitmarg->SetMinimIntervalRead(valfirstbin_X+conv::kASmallNum,vallastbin_X,valfirstbin_Y,vallastbin_Y,true);
  fitmarg->SetStartBin(start_dmsq32_first,+TMath::Pi()); 
  fitmarg->SetNumStartBin(1,1);
  fitmarg->SetGridStartRange(valfirstbin_X+conv::kASmallNum,vallastbin_X,valfirstbin_Y,vallastbin_Y+conv::kASmallNum);
  fitmarg->SetExtractStartPt("Single");
  fitmarg->MinimizeChisqStepRead(skobs_fakedata,true);      
  
  double chisqmin_4 = (double)fitmarg->GetChisqBestFit();
  double bestfit0_4 = fitmarg->GetParamBestFit(0);
  double bestfit1_4 = -999.;
  if(fitmarg->IsFit2D()) bestfit1_4 = fitmarg->GetParamBestFit(1); 
  int nfits_4 = fitmarg->GetNumFit();
  int fitstatus_4 = fitmarg->GetFitStatus();
  int niterations_4 = fitmarg->GetNumIterations();
  if(chisqmin_4<chisqmin){
    chisqmin = chisqmin_4;
    bestfit0 = bestfit0_4; 
    bestfit1 = bestfit1_4; 
    nfits = nfits_4;    
    fitstatus = fitstatus_4;
    niterations = niterations_4;
  }
  

  // Fit 5
  
  cout << endl;
  cout << "/////////////// Fit 5 (dmsq32=last_bin, dcp=-3.14) ////////////////" << endl;
  
  fitmarg->SetMinimIntervalRead(valfirstbin_X+conv::kASmallNum,vallastbin_X,valfirstbin_Y,vallastbin_Y,true);
  fitmarg->SetStartBin(start_dmsq32_last,-TMath::Pi()); 
  fitmarg->SetNumStartBin(1,1);
  fitmarg->SetGridStartRange(valfirstbin_X+conv::kASmallNum,vallastbin_X,valfirstbin_Y,vallastbin_Y+conv::kASmallNum);
  fitmarg->SetExtractStartPt("Single");
  fitmarg->MinimizeChisqStepRead(skobs_fakedata,true);      
  
  double chisqmin_5 = (double)fitmarg->GetChisqBestFit();
  double bestfit0_5 = fitmarg->GetParamBestFit(0);
  double bestfit1_5 = -999.;
  if(fitmarg->IsFit2D()) bestfit1_5 = fitmarg->GetParamBestFit(1); 
  int nfits_5 = fitmarg->GetNumFit();
  int fitstatus_5 = fitmarg->GetFitStatus();
  int niterations_5 = fitmarg->GetNumIterations();
  if(chisqmin_5<chisqmin){
    chisqmin = chisqmin_5;
    bestfit0 = bestfit0_5; 
    bestfit1 = bestfit1_5; 
    nfits = nfits_5;    
    fitstatus = fitstatus_5;
    niterations = niterations_5;
  }
    
  // Fit 6
  
  cout << endl;
  cout << "/////////////// Fit 6 (dmsq32=last_bin, dcp=-pi/2.) ////////////////" << endl;
  
  fitmarg->SetMinimIntervalRead(valfirstbin_X+conv::kASmallNum,vallastbin_X,valfirstbin_Y,vallastbin_Y,true);
  fitmarg->SetStartBin(start_dmsq32_last,-TMath::Pi()/2.); 
  fitmarg->SetNumStartBin(1,1);
  fitmarg->SetGridStartRange(valfirstbin_X+conv::kASmallNum,vallastbin_X,valfirstbin_Y,vallastbin_Y+conv::kASmallNum);
  fitmarg->SetExtractStartPt("Single");
  fitmarg->MinimizeChisqStepRead(skobs_fakedata,true);      
  
  double chisqmin_6 = (double)fitmarg->GetChisqBestFit();
  double bestfit0_6 = fitmarg->GetParamBestFit(0);
  double bestfit1_6 = -999.;
  if(fitmarg->IsFit2D()) bestfit1_6 = fitmarg->GetParamBestFit(1); 
  int nfits_6 = fitmarg->GetNumFit();
  int fitstatus_6 = fitmarg->GetFitStatus();
  int niterations_6 = fitmarg->GetNumIterations();
  if(chisqmin_6<chisqmin){
    chisqmin = chisqmin_6;
    bestfit0 = bestfit0_6; 
    bestfit1 = bestfit1_6; 
    nfits = nfits_6;    
    fitstatus = fitstatus_6;
    niterations = niterations_6;
  }
  
  
  // Fit 7
  
  cout << endl;
  cout << "/////////////// Fit 7 (dmsq32=last_bin, dcp=+pi/2.) ////////////////" << endl;
  
  fitmarg->SetMinimIntervalRead(valfirstbin_X+conv::kASmallNum,vallastbin_X,valfirstbin_Y,vallastbin_Y,true);
  fitmarg->SetStartBin(start_dmsq32_last,+TMath::Pi()/2.); 
  fitmarg->SetNumStartBin(1,1);
  fitmarg->SetGridStartRange(valfirstbin_X+conv::kASmallNum,vallastbin_X,valfirstbin_Y,vallastbin_Y+conv::kASmallNum);
  fitmarg->SetExtractStartPt("Single");
  fitmarg->MinimizeChisqStepRead(skobs_fakedata,true);      
  
  double chisqmin_7 = (double)fitmarg->GetChisqBestFit();
  double bestfit0_7 = fitmarg->GetParamBestFit(0);
  double bestfit1_7 = -999.;
  if(fitmarg->IsFit2D()) bestfit1_7 = fitmarg->GetParamBestFit(1); 
  int nfits_7 = fitmarg->GetNumFit();
  int fitstatus_7 = fitmarg->GetFitStatus();
  int niterations_7 = fitmarg->GetNumIterations();
  if(chisqmin_7<chisqmin){
    chisqmin = chisqmin_7;
    bestfit0 = bestfit0_7; 
    bestfit1 = bestfit1_7; 
    nfits = nfits_7;    
    fitstatus = fitstatus_7;
    niterations = niterations_7;
  }
  
  
  // Fit 8
  
  cout << endl;
  cout << "/////////////// Fit 8 (dmsq32=last_bin, dcp=+3.14) ////////////////" << endl;
  
  fitmarg->SetMinimIntervalRead(valfirstbin_X+conv::kASmallNum,vallastbin_X,valfirstbin_Y,vallastbin_Y,true);
  fitmarg->SetStartBin(start_dmsq32_last,+TMath::Pi()); 
  fitmarg->SetNumStartBin(1,1);
  fitmarg->SetGridStartRange(valfirstbin_X+conv::kASmallNum,vallastbin_X,valfirstbin_Y,vallastbin_Y+conv::kASmallNum);
  fitmarg->SetExtractStartPt("Single");
  fitmarg->MinimizeChisqStepRead(skobs_fakedata,true);      
  
  double chisqmin_8 = (double)fitmarg->GetChisqBestFit();
  double bestfit0_8 = fitmarg->GetParamBestFit(0);
  double bestfit1_8 = -999.;
  if(fitmarg->IsFit2D()) bestfit1_8 = fitmarg->GetParamBestFit(1); 
  int nfits_8 = fitmarg->GetNumFit();
  int fitstatus_8 = fitmarg->GetFitStatus();
  int niterations_8 = fitmarg->GetNumIterations();
  if(chisqmin_8<chisqmin){
    chisqmin = chisqmin_8;
    bestfit0 = bestfit0_8; 
    bestfit1 = bestfit1_8; 
    nfits = nfits_8;    
    fitstatus = fitstatus_8;
    niterations = niterations_8;
  }
}
//____________________________________________________________________________










void valor::util::MargFit_2D_sinsqth23_dmsq32(RunOpt *runopt,FitMarginalization *fitmarg,SuperkObservation *skobs_fakedata,double &chisqmin,double &bestfit0, double &bestfit1,int nfits,int fitstatus,int niterations)
{
    
  // Fit 1

  cout << endl;
  cout << "/////////////// Fit 1 (start dmsq = 0.0025, first octant) ////////////////" << endl;
  
  fitmarg->SetMinimIntervalRead(runopt->Freeparmin0(),runopt->Freeparmax0(),runopt->Freeparmin1(),runopt->Freeparmax1(),true);
  fitmarg->SetStartBin(0.4,0.0025); 
  fitmarg->SetNumStartBin(5,2);
  fitmarg->SetGridStartRange(runopt->Freeparmin0(),runopt->Freeparmax0(),runopt->Freeparmin1(),runopt->Freeparmax1()); // sin2th23<0.4 sure in the 1st octant
  fitmarg->SetExtractStartPt("Single");
  fitmarg->MinimizeChisqStepRead(skobs_fakedata,true);      
  
  double chisqmin_1 = (double)fitmarg->GetChisqBestFit();
  double bestfit0_1 = fitmarg->GetParamBestFit(0);
  double bestfit1_1 = -999.;
  if(fitmarg->IsFit2D()) bestfit1_1 = fitmarg->GetParamBestFit(1); 
  int nfits_1 = fitmarg->GetNumFit();
  int fitstatus_1 = fitmarg->GetFitStatus();
  int niterations_1 = fitmarg->GetNumIterations();
  if(chisqmin_1<chisqmin){
    chisqmin = chisqmin_1;
    bestfit0 = bestfit0_1; 
    bestfit1 = bestfit1_1; 
    nfits = nfits_1;    
    fitstatus = fitstatus_1;
    niterations = niterations_1;
  }
  
  // Fit 2
  
  cout << endl;
  cout << "/////////////// Fit 2 (start dmsq = 0.0025, second octant) ////////////////" << endl;
  
  fitmarg->SetMinimIntervalRead(runopt->Freeparmin0(),runopt->Freeparmax0(),runopt->Freeparmin1(),runopt->Freeparmax1(),true);
  fitmarg->SetStartBin(0.6,0.0025); 
  fitmarg->SetNumStartBin(5,2);
  fitmarg->SetGridStartRange(runopt->Freeparmin0(),runopt->Freeparmax0(),runopt->Freeparmin1(),runopt->Freeparmax1()); // start value sin2th23<0.4 sure in the 1st octant
  fitmarg->SetExtractStartPt("GridFix");
  fitmarg->MinimizeChisqStepRead(skobs_fakedata,true);      
  
  double chisqmin_2 = (double)fitmarg->GetChisqBestFit();
  double bestfit0_2 = fitmarg->GetParamBestFit(0);
  double bestfit1_2 = -999.;
  if(fitmarg->IsFit2D()) bestfit1_2 = fitmarg->GetParamBestFit(1); 
  int nfits_2 = fitmarg->GetNumFit();
  int fitstatus_2 = fitmarg->GetFitStatus();
  int niterations_2 = fitmarg->GetNumIterations();
  if(chisqmin_2<chisqmin){
    chisqmin = chisqmin_2;
    bestfit0 = bestfit0_2; 
    bestfit1 = bestfit1_2; 
    nfits = nfits_2;    
    fitstatus = fitstatus_2;
    niterations = niterations_2;
  }





  // cout << endl;
  // cout << "/////////////// Fit 1 (start dmsq = 0.002, first octant) ////////////////" << endl;  
  // fitmarg->SetMinimIntervalRead(runopt->Freeparmin0(),runopt->Freeparmax0(),runopt->Freeparmin1(),runopt->Freeparmax1(),true);
  // fitmarg->SetStartBin(0.4,0.0025); 
  // fitmarg->SetNumStartBin(5,2);
  // fitmarg->SetGridStartRange(runopt->Freeparmin0(),0.4,runopt->Freeparmin1(),0.004); // sin2th23<0.4 sure in the 1st octant
  // fitmarg->SetExtractStartPt("GridFix");
  // fitmarg->MinimizeChisqStepRead(skobs_fakedata,true);        
  // double chisqmin_1 = (double)fitmarg->GetChisqBestFit();
  // double bestfit0_1 = fitmarg->GetParamBestFit(0);
  // double bestfit1_1 = -999.;
  // if(fitmarg->IsFit2D()) bestfit1_1 = fitmarg->GetParamBestFit(1); 
  // int nfits_1 = fitmarg->GetNumFit();
  // int fitstatus_1 = fitmarg->GetFitStatus();
  // int niterations_1 = fitmarg->GetNumIterations();
  // if(chisqmin_1<chisqmin){
  //   chisqmin = chisqmin_1;
  //   bestfit0 = bestfit0_1; 
  //   bestfit1 = bestfit1_1; 
  //   nfits = nfits_1;    
  //   fitstatus = fitstatus_1;
  //   niterations = niterations_1;
  // }  
  // // Fit 2  
  // cout << endl;
  // cout << "/////////////// Fit 2 (start dmsq = 0.002, second octant) ////////////////" << endl;  
  // fitmarg->SetMinimIntervalRead(runopt->Freeparmin0(),runopt->Freeparmax0(),runopt->Freeparmin1(),runopt->Freeparmax1(),true);
  // fitmarg->SetStartBin(0.1,0.0025); 
  // fitmarg->SetNumStartBin(5,2);
  // fitmarg->SetGridStartRange(0.6,runopt->Freeparmax0(),runopt->Freeparmin1(),0.004); // start value sin2th23<0.4 sure in the 1st octant
  // fitmarg->SetExtractStartPt("GridFix");
  // fitmarg->MinimizeChisqStepRead(skobs_fakedata,true);        
  // double chisqmin_2 = (double)fitmarg->GetChisqBestFit();
  // double bestfit0_2 = fitmarg->GetParamBestFit(0);
  // double bestfit1_2 = -999.;
  // if(fitmarg->IsFit2D()) bestfit1_2 = fitmarg->GetParamBestFit(1); 
  // int nfits_2 = fitmarg->GetNumFit();
  // int fitstatus_2 = fitmarg->GetFitStatus();
  // int niterations_2 = fitmarg->GetNumIterations();
  // if(chisqmin_2<chisqmin){
  //   chisqmin = chisqmin_2;
  //   bestfit0 = bestfit0_2; 
  //   bestfit1 = bestfit1_2; 
  //   nfits = nfits_2;    
  //   fitstatus = fitstatus_2;
  //   niterations = niterations_2;
  // }

}
//____________________________________________________________________________
void valor::util::MargFit_2D_sinsqth23_sinsqth23bar(RunOpt *runopt,FitMarginalization *fitmarg,SuperkObservation *skobs_fakedata,double &chisqmin,double &bestfit0, double &bestfit1,int nfits,int fitstatus,int niterations)
{
    
    // Fit 1
    
    cout << endl;
    cout << "/////////////// Fit 1 (start dmsq = 0.0025, first octant) ////////////////" << endl;
    
    fitmarg->SetMinimIntervalRead(runopt->Freeparmin0(),runopt->Freeparmax0(),runopt->Freeparmin1(),runopt->Freeparmax1(),true);
    fitmarg->SetStartBin(0.4,0.4);
    fitmarg->SetNumStartBin(5,5);
    fitmarg->SetGridStartRange(runopt->Freeparmin0(),runopt->Freeparmax0(),runopt->Freeparmin1(),runopt->Freeparmax1()); // sin2th23<0.4 sure in the 1st octant
    fitmarg->SetExtractStartPt("Single");
    fitmarg->MinimizeChisqStepRead(skobs_fakedata,true);
    
    double chisqmin_1 = (double)fitmarg->GetChisqBestFit();
    double bestfit0_1 = fitmarg->GetParamBestFit(0);
    double bestfit1_1 = -999.;
    if(fitmarg->IsFit2D()) bestfit1_1 = fitmarg->GetParamBestFit(1);
    int nfits_1 = fitmarg->GetNumFit();
    int fitstatus_1 = fitmarg->GetFitStatus();
    int niterations_1 = fitmarg->GetNumIterations();
    if(chisqmin_1<chisqmin){
        chisqmin = chisqmin_1;
        bestfit0 = bestfit0_1;
        bestfit1 = bestfit1_1;
        nfits = nfits_1;
        fitstatus = fitstatus_1;
        niterations = niterations_1;
    }
    
    // Fit 2
    
    cout << endl;
    cout << "/////////////// Fit 2 (start ss23 = 0.6, ss23bar = 0.6, second octant) ////////////////" << endl;
    
    fitmarg->SetMinimIntervalRead(runopt->Freeparmin0(),runopt->Freeparmax0(),runopt->Freeparmin1(),runopt->Freeparmax1(),true);
    fitmarg->SetStartBin(0.6,0.6);
    fitmarg->SetNumStartBin(5,5);
    fitmarg->SetGridStartRange(runopt->Freeparmin0(),runopt->Freeparmax0(),runopt->Freeparmin1(),runopt->Freeparmax1()); // start value sin2th23 = 0.6 sure in the 2nd octant
    fitmarg->SetExtractStartPt("GridFix");
    fitmarg->MinimizeChisqStepRead(skobs_fakedata,true);
    
    double chisqmin_2 = (double)fitmarg->GetChisqBestFit();
    double bestfit0_2 = fitmarg->GetParamBestFit(0);
    double bestfit1_2 = -999.;
    if(fitmarg->IsFit2D()) bestfit1_2 = fitmarg->GetParamBestFit(1);
    int nfits_2 = fitmarg->GetNumFit();
    int fitstatus_2 = fitmarg->GetFitStatus();
    int niterations_2 = fitmarg->GetNumIterations();
    if(chisqmin_2<chisqmin){
        chisqmin = chisqmin_2;
        bestfit0 = bestfit0_2;
        bestfit1 = bestfit1_2;
        nfits = nfits_2;
        fitstatus = fitstatus_2;
        niterations = niterations_2;
    }
    
}
//____________________________________________________________________________
void valor::util::MargFit_2D_dmsq_dmsqbar(RunOpt *runopt,FitMarginalization *fitmarg,SuperkObservation *skobs_fakedata,double &chisqmin,double &bestfit0, double &bestfit1,int nfits,int fitstatus,int niterations)
{
    
    // Fit 1
    
    cout << endl;
    cout << "/////////////// Fit 1 (start dmsq = 0.0025, first octant) ////////////////" << endl;
    
    fitmarg->SetMinimIntervalRead(runopt->Freeparmin0(),runopt->Freeparmax0(),runopt->Freeparmin1(),runopt->Freeparmax1(),true);
    fitmarg->SetStartBin(0.0025,0.0025);
    fitmarg->SetNumStartBin(2,2);
    fitmarg->SetGridStartRange(runopt->Freeparmin0(),runopt->Freeparmax0(),runopt->Freeparmin1(),runopt->Freeparmax1());
    fitmarg->SetExtractStartPt("Single");
    fitmarg->MinimizeChisqStepRead(skobs_fakedata,true);
    
    double chisqmin_1 = (double)fitmarg->GetChisqBestFit();
    double bestfit0_1 = fitmarg->GetParamBestFit(0);
    double bestfit1_1 = -999.;
    if(fitmarg->IsFit2D()) bestfit1_1 = fitmarg->GetParamBestFit(1);
    int nfits_1 = fitmarg->GetNumFit();
    int fitstatus_1 = fitmarg->GetFitStatus();
    int niterations_1 = fitmarg->GetNumIterations();
    if(chisqmin_1<chisqmin){
        chisqmin = chisqmin_1;
        bestfit0 = bestfit0_1;
        bestfit1 = bestfit1_1;
        nfits = nfits_1;
        fitstatus = fitstatus_1;
        niterations = niterations_1;
    }
    
}
//____________________________________________________________________________
void valor::util::SetMassHierarchy(SuperkExpectation * skexp)
{
  //  - iMH=0 --> kNmhNormal
  //  - iMH=1 --> kNmhInverted	
	
  //convert to NuMassHierarchy_t
  if(OscParams::Instance()->iMH == 0)
    OscParams::Instance()->MH = kNmhNormal;
  else if(OscParams::Instance()->iMH == 1)
    OscParams::Instance()->MH = kNmhInverted;
  else {
    cerr << "util::SetMassHierarchy() Unknown value of iMH " << OscParams::Instance()->iMH << endl;
    exit(1);
  }
  
  //set the hierarchy
  Opt::Instance()->gMassHierarchyAssumed = OscParams::Instance()->MH;
  skexp->GetSuperkMC()->SetMassHierarchy(OscParams::Instance()->MH);
  skexp->GetSuperkMC()->SetMassHierarchyBar(OscParams::Instance()->MH);

  ////save it
  //if(OscParams::Instance()->hToy_MH)
  //OscParams::Instance()->hToy_MH->Fill(OscParams::Instance()->iMH + conv::kASmallNum);
}
//____________________________________________________________________________
//____________________________________________________________________________

