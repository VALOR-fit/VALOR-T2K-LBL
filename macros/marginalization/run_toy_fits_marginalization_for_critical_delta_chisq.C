
// ------------------------------------------------------------------------------------------------           
//
// Intended to calculate a critical delta Chi-square surface for the marginalized toys.
// Fits N random fakedata sets at each grid point.
// Saves both toy spectra and fit results.
// 
// Steve Dennis, March 2015
//                                                                                           
// ------------------------------------------------------------------------------------------------


#ifndef __CINT__
#include "TROOT.h"
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
#include "TAxis.h"
#include "TXMLEngine.h"
#include "TTree.h"
#include "TBranch.h"

#include "Conventions/Opt.h"
#include "Conventions.h"
#include "Utilities.h"
#include "Utilities/Rnd.h"
#include "Utilities/Env.h"
#include "FitCore/PhysHypothesis.h"
#include "FitCore/FitInputs.h"
#include "FitCore/FitterMinuit.h"
#include "FitCore/UtlFit.h"

#include "AnalysisOptions.h"
#include "BANFF.h"
#include "BANFFSystConv.h"
#include "Binning.h"
#include "Build.h"
#include "Conventions.h"
#include "Enumerations.h"
#include "FitMarginalization.h"
#include "Formatting.h"
#include "NuFlux.h"
#include "NuXSec.h"
#include "SamplesAndModes.h"
#include "SKFSISystConv.h"
#include "SuperkDetectorFSI.h"
#include "SuperkExpectation.h"
#include "SuperkMC.h"
#include "SuperkModes.h"
#include "SuperkObservation.h"
#include "SuperkSamples.h"
#include "T2KOpt.h"
#include "T2KPhysicalParameterLimits.h"
#include "T2KPhysParamId.h"
#include "Utilities.h"
#include "UtlFitParamListFactory.h"
#include "UtlT2KStyle.h"
#include "XSecSystConv.h"

#endif

void run_toy_fits_marginalization_for_critical_delta_chisq
(
  const bool   readtoys             = true,
  const char * outfile_tag          = "marginalization_toys",
  
  const int    Nfakedata           = 10000,
  
  const int    IDoutfile            = 0,
  const int    first_bin_freepar0   = -999,   // first freepar of the sub-range
  const int    last_bin_freepar0    = -999, // last freepar of the sub-range
  const int    first_bin_freepar1   = -999,   // first freepar of the sub-range
  const int    last_bin_freepar1    = -999, // last freepar of the sub-range
  
  const char * infilename           = "toys",
  
  const int    nfiles               = 1,         // # of files to read
  
  const char * config_file          = "analysis_configuration_2.5E20.xml",  
  
  const char * fit_param_list       = "nuosc3fbar_numudisapp_fixnonbar_fix4+nuisance",
  
  const int    NtoysMarg            = 10000,   // # toys marginalization 
  
  const char * freeparname0        = "mixing23Bar",  // free parameter
  const double freeparmin0         = 0.,             // minimum value of free parameter
  const double freeparmax0         = 1.,             // maxinum value of free parameter
  const int    NptScan0            = 21,             // # pts of free parameter to scan
  
  const char * freeparname1        = "dmsq23Bar",  // free parameter
  const double freeparmin1         = 0.0002, // minimum value of free parameter
  const double freeparmax1         = 0.02, // maxinum value of free parameter
  const int    NptScan1            = 21,      // # pts of free parameter to scan
  
  int         tweak_systematics    = 1,       /* generate toy MC with randomized systematics? */
  
  const bool  doMargMH             = true,    // marginalize Mass Hierarchy
  const bool  doRethrowUnphysical  = true,    // rethrow unphysical systematic toys
  
  const bool  doSaveFakedataTree   = true,    // keep spectra and parameters for fakedata
  
  const char * dmsqscheme          = "abs23",      /* abs23==kNuDmsq12andAbs23 abs31==kNuDmsq12andAbs31 fogli==kNuDmsqFogliLisi */
  const char * hierarchy           = "normal",     /* normal or inverted */
  const char * mixparamchoice      = "sinsqtheta", /* sinsq2theta or sinsqtheta */
  const char * mixparamoctant      = "one"       /* one==kNuOctantOne two==kNuOctantTwo */
)
{
  
  using namespace std;
  using namespace valor;
  using namespace valor::conv;
  
  TStopwatch w1,w2,w3;
  w1.Start();
  w2.Start();
  w3.Start();

  // Global options
  Opt::Instance()->gMassHierarchyAssumed = util::HierarchyFromString(hierarchy);
  Opt::Instance()->gOscDmsqScheme        = util::DmsqSchemeFromString(dmsqscheme);
  Opt::Instance()->gOscMixingParamChoice = util::MixingTypeFromString(mixparamchoice);
  Opt::Instance()->gOscOctantChoice      = util::OctantFromString(mixparamoctant);
  if(Opt::Instance()->Undefined()) {
    cerr << "Global options not defined correctly. Exiting..." << endl;
    cerr << *Opt::Instance() << endl;
    exit(1);
  }
  Opt::Instance()->SetDefaultParameterValues();
  
  FileStat_t dummy;

  TString config    = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);

  // Check for config file, and XML data file.
  // If not found in inputs or config, see if it works as a normal path
  if (gSystem->GetPathInfo(config.Data(), dummy)) 
  {
    config = config_file;
    if (gSystem->GetPathInfo(config.Data(), dummy)) {
      cerr<<"Config file not found: "<<config_file<<endl;
      exit(-1);
    }
  }
  
  // Get superk spectrum predictor & configure
  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read(config.Data());
  
  // Work out the parameters
  bool fitIs1D = false;
  PhysParamId_t freepar0 = util::PhysParamIdFromString(freeparname0);
  if (freepar0 < 0 || util::IsSystematic(freepar0))
  {
    cout << "Cannot handle parameter 0 "<< freeparname0 <<", p= "<<freepar0 << "\t"<<util::AsString(freepar0) << endl;
    exit(1);
  }

  PhysParamId_t freepar1;
  string sfreeparname1(freeparname1);  
  if ( sfreeparname1=="" || sfreeparname1=="none" )
  {
    cout << "1D fit" << endl;
    freepar1 = conv::kPhypUndefined;
    fitIs1D = true;
  }
  else {
    freepar1 = util::PhysParamIdFromString(freeparname1);
    
    if (freepar1 < 1 || util::IsSystematic(freepar1))
    {
      cout << "Cannot handle parameter 1 "<< freeparname1 <<", p= "<<freepar1 << "\t"<<util::AsString(freepar1) << endl;
      exit(1);
    }
  }
  
  // Set free param, ntoys, # pt to scan, min par, max par (if min and max not passed, use default ones for that par)

  // free parameter 
  PhysParamId_t parfree[2];
  parfree[0] = freepar0;
  parfree[1] = freepar1;
  
  // # of toys for marginalization
  int ntoys = NtoysMarg;

  // # of free parameter pts to scan 
  int nbinspar[2];
  nbinspar[0] = NptScan0;
  nbinspar[1] = NptScan1;

  // free parameter range
  double minpar[2];
  double maxpar[2];
  minpar[0] = freeparmin0;
  maxpar[0] = freeparmax0;
  minpar[1] = freeparmin1;
  maxpar[1] = freeparmax1;

  // free parameter sub-range 
  int subminbin[2];
  int submaxbin[2];
  subminbin[0] = first_bin_freepar0; submaxbin[0] = last_bin_freepar0;
  subminbin[1] = first_bin_freepar1; submaxbin[1] = last_bin_freepar1;
  
  // Implement the constructor giving the following input parameters:
  FitMarginalization *fitmarg = new FitMarginalization();   
  fitmarg->SetPhysParFree(parfree,ntoys,nbinspar,minpar,maxpar);
  fitmarg->SetSubRange(subminbin,submaxbin); // By default it scans all the range
  fitmarg->Init(); // Before need SetPhysParFree() to initialize the arrays
  
  // Set fit parameter list
  fitmarg->SetFitParamList(fit_param_list);
  fitmarg->PrintPhysHypothesis();
      
  // Set Mass Hierarchy marginalization
  fitmarg->SetMargMassHierarchy(doMargMH); // false by default 
  
  fitmarg->BuildDefaultPriors();
  
  const FitParamList * fparm = fitmarg->GetFitParamList();
  
  // Update the variables with the actual bin numbers
  int real_first_bin_freepar0 = fitmarg->GetFirstBinInSubRange(0);
  int real_first_bin_freepar1 = fitmarg->GetFirstBinInSubRange(1);
  int real_last_bin_freepar0  = fitmarg->GetLastBinInSubRange (0);
  int real_last_bin_freepar1  = fitmarg->GetLastBinInSubRange (1);
  
  int n_bins_this_job_par0 = real_last_bin_freepar0-real_first_bin_freepar0+1;
  int n_bins_this_job_par1 = 1;
  
  if (freepar1!=conv::kPhypUndefined) {
    n_bins_this_job_par1 = real_last_bin_freepar1-real_first_bin_freepar1+1;
  }
  
  int n_bins_this_job = n_bins_this_job_par0*n_bins_this_job_par1;
  
  // Set up output file for toys
  TFile * outfile = new TFile(Form("%s_%d.root",outfile_tag,IDoutfile),"RECREATE");
  
  
  TTree * fakedata_tree = (doSaveFakedataTree ? new TTree("fakedata","fakedata") : 0);
  
  VALORArrayF   fillva[conv::kns]; // Used for filling branch
  PhysParamList fillpp   ; // Used for initialising branch
  
  TBranch * spectrum_branches[conv::kns] = {};
  TBranch * params_branch = 0;
  
  int ibin,ipar0,ipar1;
  double value=-999;
  
  if (doSaveFakedataTree) {
    params_branch  = fakedata_tree->Branch("parameters"   ,"valor::PhysParamList",&fillpp);
    for (int is = 0 ; is < conv::kns ; is++) {
      SuperkSample_t s = conv::kSKSamples[is];
      if (s==conv::kSKsUndefined) {continue;}
      string branch_name(Form("ErecoSpectrum%s",util::AsStringCompact(s).c_str()));
      spectrum_branches[is] = fakedata_tree->Branch(branch_name.c_str(),"valor::VALORArrayF"   ,&fillva[is]);
    }
    
    fakedata_tree->Branch("ibin_par0" ,&ipar0   ,"ibin_par_0/I");
    fakedata_tree->Branch("ibin_par1" ,&ipar1   ,"ibin_par_1/I");
    fakedata_tree->Branch("freepar0"  ,&freepar0,"freepar0/I");
    fakedata_tree->Branch("freepar1"  ,&freepar1,"freepar1/I");
  }
  
  // Start generating fake datasets
  SuperkObservation * skobs = new SuperkObservation;
  
  VALORArrayF   ** fFakeDataEreco = (VALORArrayF  **) calloc(n_bins_this_job*Nfakedata*conv::kns, sizeof(VALORArrayF   *));
  PhysParamList ** fTrueValues    = (PhysParamList**) calloc(n_bins_this_job*Nfakedata          , sizeof(PhysParamList *));
  
  PhysHypothesis ph(*fparm);
  
  std::cout<<"Setup is complete, CPU time is "<<Form("%.2lf",w1.CpuTime())<<" seconds, walltime is "<<Form("%.2lf",w1.RealTime())<<" seconds"<<std::endl;
  
  std::cout<<"Generating "<<Nfakedata<<" fake datasets each for "<<n_bins_this_job<<" bins"<<std::endl;
  
  ibin = 0;
  for (ipar0 = real_first_bin_freepar0 ; ipar0 <= real_last_bin_freepar0 ; ipar0++) {
    for (ipar1 = real_first_bin_freepar1 ; ipar1 <= real_last_bin_freepar1 ; ipar1++) {
      for(int ifakedata = 0 ; ifakedata < Nfakedata; ifakedata++) {
        bool got_good_tweaks = false;
        while (true) {
          skexp->RandomizeSystematics();
          PhysParamList true_nuisance_params = skexp->Tweaks();
          got_good_tweaks = (!doRethrowUnphysical) || T2KPhysicalParameterLimits::Instance()->AllWithinLimits(true_nuisance_params);
          if (got_good_tweaks) {
            fTrueValues[ifakedata] = new PhysParamList(true_nuisance_params);
            break;
          }
        }
        
        for (int ip = 0 ; ip < conv::kNuis000; ip++) {// loop over physics parameters
          PhysParamId_t p = (PhysParamId_t)ip;
          
          if (!(fparm->Exists(p))) {continue;}
          
          if (p==freepar0) {
            value = fitmarg->GetBinParameterValue(0,ipar0);
          }
          else if (p==freepar1){
            value = fitmarg->GetBinParameterValue(1,ipar1);
          }
          else if (!fparm->IsFixed(p)) {
            value = fitmarg->Randomize(p);
          }
          else {
            value = fparm->CurrValue(p);
          }
          
          fTrueValues[ifakedata]->SetValue(p,value);
          
          ph.SetCurrValue(p,value);
        }
        
        skobs->Generate(skexp,ph,true);
        
        for (int is = 0; is < conv::kns ; is++) {
          int index = ibin*Nfakedata*conv::kns + conv::kns*ifakedata+is;
          SuperkSample_t s = conv::kSKSamples[is];
          if (s==conv::kSKsUndefined) {
            fFakeDataEreco[index] = new VALORArrayF(1);
          }
          else {
            fFakeDataEreco[index] = new VALORArrayF(*(skobs->ErecoSpectrumVArray(s)));
          }
        }
        
        if (doSaveFakedataTree) {
          for (int is = 0; is < conv::kns ; is++) {
            int index = ibin*Nfakedata*conv::kns + conv::kns*ifakedata+is;
            if (spectrum_branches[is]) {
              fillva[is].CopyFrom(*fFakeDataEreco[index]);
            }
          }
          fillpp.CopyFrom(*fTrueValues[ifakedata]);
          fakedata_tree->Fill();
        } //save
        
        //~ cout<<"Filled fakedata "<<ifakedata<<" for bin "<<ipar0<<","<<ipar1<<endl;
      } //ifakedata
      ibin++;
      cout<<"Made all fakedata for bin "<<ipar0<<","<<ipar1<<endl;
    } //ipar1
  } //ipar0
  
  outfile->cd();
  if (doSaveFakedataTree) {fakedata_tree->Write();}
  
  std::cout<<"Fake data generation is complete, CPU time is "<<Form("%.2lf",w2.CpuTime())<<" seconds, walltime is "<<Form("%.2lf",w2.RealTime())<<" seconds"<<std::endl;
  // 
  
  
  //Set up fit output tree.
  double min_chisq=-999,this_chisq=-999,freepar0_val=-999,freepar1_val=-999;
  outfile->cd();
  TTree * fit_result_tree = new TTree("fit_results","fit_results");
  
  fit_result_tree->Branch("ibin_par0"    ,&ipar0        ,"ibin_par_0/I");
  fit_result_tree->Branch("ibin_par1"    ,&ipar1        ,"ibin_par_1/I");
  fit_result_tree->Branch("freepar0"     ,&freepar0     ,"freepar0/I");
  fit_result_tree->Branch("freepar1"     ,&freepar1     ,"freepar1/I");
  //fit_result_tree->Branch("min_chisq"    ,&min_chisq    ,"min_chisq/D");
  fit_result_tree->Branch("chisq"        ,&this_chisq   ,"chisq/D    ");
  fit_result_tree->Branch("freepar0_val" ,&freepar0_val ,"freepar0_val/D");
  fit_result_tree->Branch("freepar1_val" ,&freepar1_val ,"freepar1_val/D");
  
  // Read toys from binary files and fill marginalized chisq
  for(int ifile=0;ifile<nfiles;ifile++){
    
  
    cout << endl;
    cout << "////////////////////////////////" << endl;
    cout << endl;
    cout << "Read file ID " << ifile << endl; 
    cout << endl;

    fitmarg->ReadToysFromBinary(ifile,infilename);
    
    for (int ifakedata = 0 ; ifakedata < Nfakedata; ifakedata++) {
      
      ibin=0;
      if (ifakedata % 500 == 0) {
        cout<<"Calculating chi-square for fake data set "<<ifakedata<<std::endl;
      }
      
      if (fitIs1D) {
        for (ipar0 = real_first_bin_freepar0 ; ipar0 <= real_last_bin_freepar0 ; ipar0++) {
          
          // Copy the right fakedata into SKobs
          for (int is = 0; is < conv::kns ; is++) {
            int index = ibin*Nfakedata*conv::kns + conv::kns*ifakedata+is;
            SuperkSample_t s = conv::kSKSamples[is];
            if (s==conv::kSKsUndefined) {continue;}
            skobs->ErecoSpectrumVArray(s)->CopyFrom(*fFakeDataEreco[index]);
          }//is
          
          // Get the fit point
          freepar0_val = fitmarg->GetBinParameterValue(0,ipar0);
          this_chisq   = fitmarg->CalcChisqMargPt(skobs,ipar0,0,false);
          
          ibin++;
          fit_result_tree->Fill();
        }//ipar0
      }//fitIs1D
      else {
        for (ipar0 = real_first_bin_freepar0 ; ipar0 <= real_last_bin_freepar0 ; ipar0++) {
          for (ipar1 = real_first_bin_freepar1 ; ipar1 <= real_last_bin_freepar1 ; ipar1++) {
            
            // Copy the right fakedata into SKobs
            for (int is = 0; is < conv::kns ; is++) {
              int index = ibin*Nfakedata*conv::kns + conv::kns*ifakedata+is;
              SuperkSample_t s = conv::kSKSamples[is];
              if (s==conv::kSKsUndefined) {continue;}
              skobs->ErecoSpectrumVArray(s)->CopyFrom(*fFakeDataEreco[index]);
            }//is
          
            freepar0_val = fitmarg->GetBinParameterValue(0,ipar0);
            freepar1_val = fitmarg->GetBinParameterValue(1,ipar1);
            this_chisq = fitmarg->CalcChisqMargPt(skobs,ipar0,ipar1,false);
            
            ibin++;
            fit_result_tree->Fill();
          }//ipar1
        }//ipar0
      }//2d fit
      //fit_result_tree->AutoSave();
    }//ifakedata
  }  //ifile
  
  outfile->cd();
  fit_result_tree->Write();
  
  std::cout<<"Job finished, CPU time is "<<Form("%.2lf",w3.CpuTime())<<" seconds, walltime is "<<Form("%.2lf",w3.RealTime())<<" seconds"<<std::endl;
  
  return;
}



#ifndef __CINT__
int main()
{
  valor::Env::Instance()->AutoSet("t2xk");
  TString conf_dir  = valor::Env::Instance()->TopLevelConfDir();
  TString verbosity = conf_dir + "/verbosity.xml";
  valor::util::Verbosity(verbosity.Data());

  run_toy_fits_marginalization_for_critical_delta_chisq(true,"marginalization_toys",1000,0,-999,-999,-999,-999,"toys",1,"analysis_configuration_2.5E20.xml","nuosc3fbar_numudisapp_fixnonbar_fix4+nuisance",10000,"mixing23Bar",0,1,21,"dmsq23bar",2E-4,2E-2,21,1);
  
  return 0;
}
#endif
