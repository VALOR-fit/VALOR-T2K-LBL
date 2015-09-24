//____________________________________________________________________________
/*

Copyright (c) 2010-2014, The VALOR group. All rights reserved.

ChangeLog:

*/
//____________________________________________________________________________


#include "UtlForMacros.h"

#include "SuperkExpectation.h"
#include "SuperkObservation.h"
#include "SuperkMC.h"
#include "BANFF.h"
#include "SuperkDetectorFSI.h"
#include "FitMarginalization.h"
#include "Build.h"
#include "UtlForMacros.h"
#include "OscParams.h"
#include "RunOpt.h"
#include "Utilities/Env.h"

#include "Conventions/Opt.h"
#include "Conventions/PhysParamId.h"
#include "Conventions/PhysHypothesisId.h"
#include "Conventions/Constants.h"
#include "FitCore/PhysHypothesis.h"
#include "Utilities/UtlXML.h"

#include "TSystem.h"
#include "TCanvas.h"
#include "TMath.h"

#include <iostream>
#include <string>

using std::cout;
using std::cerr;
using std::endl;
using std::string;

using namespace valor;
using namespace valor::conv;

//____________________________________________________________________________
const char * valor::util::FormFilePath(const char * file, const char * dir)
{
  FileStat_t dummy;
  const char * filename = Form("%s/%s/%s", gSystem->Getenv("T2KOSCANA"), dir, file);
  if(gSystem->GetPathInfo(filename, dummy)) {
    cerr << "File not found: " << filename << endl;
    exit(-1);
  }
  return filename; 
}
//____________________________________________________________________________
string valor::util::SysThrowIdToString(unsigned int tweak_systematics)
{
  string tweak_sys_str;
  switch(tweak_systematics)
    {
    case(0):
      tweak_sys_str = "None"; break;
    case(1):
      tweak_sys_str = "All"; break;
    case(2):
      tweak_sys_str = "Special low energy numu flux increase"; break;
    case(3):
      tweak_sys_str = "Tweak all but SK and FSI"; break;
    case(4):
      tweak_sys_str = "NC1pi norm = 0.5"; break;
    case(5):
      tweak_sys_str = "High energy >1.5 GeV flux increase"; break;
    case(6):
      tweak_sys_str = "BANFF"; break;
    case(7):
      tweak_sys_str = "XSec"; break;
    case(8):
      tweak_sys_str = "SK+FSI"; break;
    default:
      cerr << "SysThrowIdToString() tweak_systematics = " << tweak_systematics
	   << " not implemented" << endl;
      exit(-1);
    }
  return tweak_sys_str;
}
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
    cout << "- Randomizing all systematic params in pdf..." << endl;
    skexp->TweakAllSystematics();
    // skexp->TweakSystematic("sk_energy_scale",0.0,true);
  }
  else if(tweak_systematics == 2) {
    // special tweak, increasing numu flux parameters <700MeV by 5%
    cout << "- Increasing numu flux parameters <700MeV by 5% in pdf..." << endl;
    TVectorD twkv(skexp->GetBANFF()->BestFitValues());
    twkv[0] += 0.05;
    twkv[1] += 0.05;
    twkv[2] += 0.05;
    twkv[3] += 0.05;
    skexp->GetBANFF()->Tweak(true,&twkv);
  }
  else if(tweak_systematics == 3) {
    cout << "- Randomizing all but SK systematic params in pdf..." << endl;
    skexp->TweakAllBANFFSystematics();
    skexp->TweakAllUncorrelatedXSecSystematics();
  }
  else if(tweak_systematics == 4) {
    cout << "- Setting xsec_norm_nc1piC to -0.5" << endl;
    skexp->TweakSystematic("xsec_norm_nc1piC",-0.5);
  }
  else if(tweak_systematics == 5) {
    cout << "- Increasing flux parameters >1500MeV by 15% in pdf..." << endl;
    TVectorD twkv(skexp->GetBANFF()->BestFitValues());
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
    cout << "- Randomizing all BANFF systematic params in pdf..." << endl;
    skexp->TweakAllBANFFSystematics();
  }
  else if(tweak_systematics == 7) {
    cout << "- Randomizing all uncorrelated cross section systematic params in pdf..." << endl;
    skexp->TweakAllUncorrelatedXSecSystematics();
  }
  else if(tweak_systematics == 8) {
    cout << "- Randomizing all SK+FSI systematic params in pdf..." << endl;
    skexp->TweakAllSKDetectorFSISystematics();
  }
  else {
    cerr << "PerformSysThrow() tweak_systematics = " << tweak_systematics
	 << " not implemented" << endl;
    exit(-1);
  }
  if(print)
    skexp->PrintSystematicVariations();
}
//____________________________________________________________________________
void valor::util::PerformMargFit(int argc, char ** argv, FitMarginalization *fitmarg, SuperkExpectation *skexp, SuperkObservation *skobs_fakedata, bool /*print*/)
{    
  Env::Instance()->AutoSet("t2xk");
  cout << *Env::Instance() << endl;
  //gROOT->SetBatch(true);
  

  RunOpt runopt(argc, argv, kParseCommon | kParseToys | kParseGlobal | kParseFitMarg | kParse3FOscParams);

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
  string sfreeparname0(runopt.Freeparname0()); 
  PhysParamId_t freepar1 = conv::kPhypUndefined;
  string sfreeparname1(runopt.Freeparname1());  
     
  PhysParamId_t parfree[2];
  parfree[0] = GetFreePar(sfreeparname0);
  parfree[1] = freepar1;

  // # of toys for marginalization
  int ntoys = runopt.NtoysMarg();
  
  // # of free parameter pts to scan 
  int nbinspar[2];
  nbinspar[0] = runopt.NptScan0();
  nbinspar[1] = runopt.NptScan1();
  
  // free parameter range
  double minpar[2];
  double maxpar[2];
  minpar[0] = runopt.Freeparmin0();
  maxpar[0] = runopt.Freeparmax0();
  minpar[1] = runopt.Freeparmin1();
  maxpar[1] = runopt.Freeparmax1();
  
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
  fitmarg->SetTestStatistic(runopt.Teststat()); // Set Test Statistic ("rate-only" or "rate+shape")
  fitmarg->Init(skexp); // Before need SetPhysParFree() to initialize the arrays
  
  // Set Mass Hierarchy marginalization
  Opt::Instance()->gMassHierarchyAssumed = util::HierarchyFromString(runopt.Hierarchy_marg()); // needed to build the right binary file name 
  
  //
  // MH MUST BE THE SAME AS FOR TOYS!!! TO CHECK  
  //
  fitmarg->SetMargMassHierarchy(runopt.DoMargMH()); 
  //
  
  // Set fit parameter list
  fitmarg->SetFitParamList(runopt.Fit_param_list().c_str());
  fitmarg->PrintPhysHypothesis();
  
  
  // Build the map fileID - sub-range
  
  fitmarg->SetMaxIterationSingleFit(2);  // Max # of iterations per fit
  fitmarg->SetMaxFitNumber(3); // Max # of fits max # of iterations
  
  fitmarg->SetTagFileName(runopt.Infilename().c_str());
  fitmarg->SetFileID(runopt.Nfiles());


  // 1D deltaCP fit
  if(runopt.Freeparname0()=="deltacp" && runopt.Freeparname1()==""){
    util::MargFit_1D_deltaCP(fitmarg,skobs_fakedata,chisqmin,bestfit0,bestfit1,nfits,fitstatus,niterations);
  }
  
  // 1D sinsqth13 fit
  //
  // 1D sinsqth23 fit
  //
  // 1D dmsq32 fit
  // 2D sinsqth13 Vs dcp fit
  //
  // 2D sinsqth23 Vs dmsq32 fit
  //
  
  // Set best-fit parameters in FitMarginalization object
  // Otherwise error computation is done with respect to a wrong minimum
  fitmarg->SetChisqBestFit(chisqmin);
  fitmarg->SetParamBestFit(0,bestfit0);
  if(fitmarg->IsFit2D()) fitmarg->SetParamBestFit(1,bestfit1);
  
  // Build the DeltaChisqMarg after the best-fit Chisq has been set
  // Otherwise DeltaChisqMarg becomes negative if shifted wrt a local minimum
  fitmarg->BuildDeltaChisqMarg();
  
  //
  // Get fit errors
  //
  // Since the chisq was minimized reading binary files at each step you 
  // must adopt the same approach for error computations as well.
  //
  
  if(runopt.TypeError().Contains("MINOS")) fitmarg->CalcBestFitErrors(skobs_fakedata,"MINOS",true);
  if(runopt.TypeError().Contains("HESSE")) fitmarg->CalcBestFitErrors(skobs_fakedata,"HESSE",true);

  // Draw Full Chisq
  if(runopt.DrawFullChisq()){
    fitmarg->FillChisqMarg(skobs_fakedata,true);
    fitmarg->MinimizeHisto();
    //bestfitFullChisq0 = fitmarg->GetParamBestFit(0);
    //if(fitmarg->IsFit2D()) bestfitFullChisq1 = fitmarg->GetParamBestFit(1); 
    //chisqminFullChisq = (double)fitmarg->GetChisqBestFit();
    fitmarg->BuildDeltaChisqMarg();
  }

  // Store toy likelihoods (needed to draw best-fit spectra w/ marginalization)
  fitmarg->WriteLikeToys2Binary(runopt.FileTag().c_str());
  
}
//____________________________________________________________________________
conv::PhysParamId_t valor::util::GetFreePar(string freeparname){

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
    cout << "Param is not defined" << endl;
    exit(1);
  }

  return freepar;
}
//____________________________________________________________________________
void valor::util::MargFit_1D_deltaCP(FitMarginalization *fitmarg,SuperkObservation *skobs_fakedata,double &chisqmin,double &bestfit0, double &bestfit1,int nfits,int fitstatus,int niterations){
  
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
  cout << "################# Fit 2 (dcp=0.) #################" << endl;
    
  fitmarg->SetMinimIntervalRead(-TMath::Pi(),+TMath::Pi(),0.,0.,true);
  fitmarg->SetStartBin(0.,0.); 
  fitmarg->SetNumStartBin(1,1);
  fitmarg->SetGridStartRange(0.,0.+conv::kASmallNum,0.,0.);
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
  cout << "################# Fit 3 (dcp=+3.14) #################" << endl;
    
  fitmarg->SetMinimIntervalRead(-TMath::Pi(),+TMath::Pi(),0.,0.,true);
  fitmarg->SetStartBin(+TMath::Pi(),0.); 
  fitmarg->SetNumStartBin(1,1);
  fitmarg->SetGridStartRange(+TMath::Pi()-conv::kASmallNum,+TMath::Pi(),0.,0.);
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
}
//____________________________________________________________________________
void valor::util::MargFit_1D_sinsqth13(FitMarginalization * /*fitmarg*/,SuperkObservation * /*skobs_fakedata*/,double & /*chisqmin*/,double & /*bestfit0*/, double & /*bestfit1*/,int /*nfits*/,int /*fitstatus*/,int /*niterations*/){
}
//____________________________________________________________________________
void valor::util::MargFit_2D_sinsqth13_deltaCP(FitMarginalization * /*fitmarg*/,SuperkObservation * /*skobs_fakedata*/,double & /*chisqmin*/,double & /*bestfit0*/, double & /*bestfit1*/,int /*nfits*/,int /*fitstatus*/,int /*niterations*/){
}
//____________________________________________________________________________
void valor::util::MargFit_1D_sinsqth23(FitMarginalization * /*fitmarg*/,SuperkObservation * /*skobs_fakedata*/,double & /*chisqmin*/,double & /*bestfit0*/, double & /*bestfit1*/,int /*nfits*/,int /*fitstatus*/,int /*niterations*/){
}
//____________________________________________________________________________
void valor::util::MargFit_2D_sinsqth23_dmsq32(FitMarginalization * /*fitmarg*/,SuperkObservation * /*skobs_fakedata*/,double & /*chisqmin*/,double & /*bestfit0*/, double & /*bestfit1*/,int /*nfits*/,int /*fitstatus*/,int /*niterations*/){
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
TH2D * valor::util::GetTH2DFromTFile(TFile * f, const char * name)
{
  TH2D * h = 0;
  f->GetObject(name, h);
  if(!h || h->IsZombie()) {
    cerr << "Error: the histogram " << name << " hasn't been found" << endl;
    exit(1);
  }
  return h;
}
//____________________________________________________________________________
//____________________________________________________________________________

