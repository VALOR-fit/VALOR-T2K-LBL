//
// ------------------------------------------------------------------------------------------------
// Test best-fit point shift as a result of tweaking a single systematic parameter by +/- 1 sigma.
// Adapted from test_systematic_error_propagation.C
// The best-fit parameter shift is written out in a text file.
// One can assemble all such text files and produce the `star plot'.
//
// The input systematic should be given as a PhysParamId_t number
//      
// Steve Dennis, adapted from Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//



#ifndef __CINT__
#include "iostream"
#include "fstream"

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
#include "TMatrixD.h"
#include "TAxis.h"
#include "TXMLEngine.h"

#include "Conventions/Opt.h"
#include "Conventions.h"
#include "Utilities.h"
#include "Utilities/Rnd.h"
#include "Utilities/Env.h"
#include "FitCore/PhysHypothesis.h"
#include "SuperkExpectation.h"
#include "SuperkObservation.h"
#include "FitCore/FitInputs.h"
#include "FitCore/FitterMinuit.h"
#include "FitCore/UtlFit.h"
#include "BANFF.h"
#include "SuperkDetectorFSI.h"
#include "NuXSec.h"
#include "UtlFitParamListFactory.h"
#endif


void test_bestfit_point_shift_due_to_syst_twk
(
 const char * systematic                                                ,  
 const char * sksample                   = "1Rmu"                       ,  /* options: "1Rmu", "1Re" */
 const char * fit_param_list             = "nuosc3f_numudisapp_fix4"    ,  /* */
 const char * fit_option                 = "minuit"                     ,  /* oscillation fitter (see run_fitter.C) */
 double       input_mixing12             = -99                          ,  /* */
 double       input_mixing13             = -99                          ,  /* */
 double       input_mixing23             = -99                          ,  /* */
 double       input_deltacp              = -99                          ,  /* */
 double       input_dm12sq               = -99                          ,  /* */
 double       input_dm23sq               = -99                          ,  /* */
 const char * config_file                = "analysis_configuration.xml" ,  /* */
 bool         correlated_parameters      = false                        ,  /* tweak matrix parameters respecting correlations? */
 const int    throw_n_toys               = 0                               /* if <=0, will use +1 and -1 sigma tweaks. Otherwise, throw n toys & find standard deviation */
 )
{
  using namespace std;
  using namespace valor;
  using namespace valor::conv;
  
  string sys(systematic);
  bool do_nominal = (sys == "nominal");
  
  // Analysis configuration
  const char * config  = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);
  
  // Create & config superk spectrum predictor
  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read(config);
  
  // Global options
  Opt::Instance()->gOscMixingParamChoice = kNuMixSinsqTheta;
  Opt::Instance()->gMassHierarchyAssumed = kNmhNormal;
  
  Opt::Instance()->SetDefaultParameterValues();

  // Get requested sample
  SuperkSample_t s = util::SKSampleFromString(sksample);
  assert(s != kSKsUndefined);
  
  // Set physics hypothesis
  if(Opt::Instance()->gOscMixingParamChoice == kNuMixSinsq2Theta) {
    if(input_mixing13 < -1) input_mixing13 = kGlobFitSinsq2Theta13;
    if(input_mixing23 < -1) input_mixing23 = kGlobFitSinsq2Theta23;
    if(input_mixing12 < -1) input_mixing12 = kGlobFitSinsq2Theta12;
  }
  else {
    if(input_mixing13 < -1) input_mixing13 = kGlobFitSinsqTheta13;
    if(input_mixing23 < -1) input_mixing23 = kGlobFitSinsqTheta23;
    if(input_mixing12 < -1) input_mixing12 = kGlobFitSinsqTheta12;
  }
  if(input_dm23sq < -1)  input_dm23sq = kGlobFitDm23sq;
  if(input_dm12sq < -1)  input_dm12sq = kGlobFitDm12sq;
  if(input_deltacp < -90) input_deltacp = kGlobFitDeltaCP;

  // Instantiate fit params & set specified inputs as necessary
  FitParamList * fparm = util::NamedParameterList(fit_param_list);
  PhysHypothesis_t hyp = fparm->PhysicsHypothesis();
  cout << "- The following physics hypothesis is used in the fit: " << util::AsString(hyp) << endl;
  if(hyp == kHypOsc3F) {
    fparm->SetInitValue(kPhypMixing12, input_mixing12);
    fparm->SetInitValue(kPhypMixing13, input_mixing13);
    fparm->SetInitValue(kPhypMixing23, input_mixing23);
    fparm->SetInitValue(kPhypCPV    ,  input_deltacp);
    fparm->SetInitValue(kPhypDm12sq,   input_dm12sq);
    fparm->SetInitValue(kPhypDm23sq,   input_dm23sq);
  }
  if(hyp == kHypOsc3FBar) {
    fparm->SetInitValue(kPhypMixing12   , input_mixing12);
    fparm->SetInitValue(kPhypMixing13   , input_mixing13);
    fparm->SetInitValue(kPhypMixing23   , input_mixing23);
    fparm->SetInitValue(kPhypCPV        ,  input_deltacp);
    fparm->SetInitValue(kPhypDm12sq     ,   input_dm12sq);
    fparm->SetInitValue(kPhypDm23sq     ,   input_dm23sq);
    fparm->SetInitValue(kPhypMixing12Bar, input_mixing12);
    fparm->SetInitValue(kPhypMixing13Bar, input_mixing13);
    fparm->SetInitValue(kPhypMixing23Bar, input_mixing23);
    fparm->SetInitValue(kPhypCPVBar     ,  input_deltacp);
    fparm->SetInitValue(kPhypDm12sqBar  ,   input_dm12sq);
    fparm->SetInitValue(kPhypDm23sqBar  ,   input_dm23sq);
  }
  else if(hyp == kHypNuMuDisapp2F) {
    fparm->SetInitValue(kPhypMixing23, input_mixing23);
    fparm->SetInitValue(kPhypDm23sq,   input_dm23sq);
  }  
  else {
    return;
  }
  
  PhysParamId_t p = util::PhysParamIdFromString(systematic);
  
  int first_banff = (int) (util::IsRHCSample(s) ? kRHCNuisBANFF00 : kNuisBANFF00 );
  int first_skdet_fsi  = (util::IsRHCSample(s) ? kRHCNuisSKDetFSI000 : kNuisSKDetFSI000 );
  
  const TMatrixD & cov_skdetfsi = (util::IsRHCSample(s) ? skexp->GetSuperkDetectorFSI()->RHCCovMatrix() : skexp->GetSuperkDetectorFSI()->CovMatrix() );
  
  if (!do_nominal) {
    if (!fparm->Exists(p)) {
      fparm->Print(std::cerr);
      cerr<<"\n\nFitParamList does not contain parameter "<<util::AsString(p) << " (p="<<p<<")"<<endl;
      exit(-1);
    }
  }

  bool fluctuate_stat = false; // <- don't fluctuate, best-fit point shift is purely due to systematic tweak

  // Print options
  cout << "\n[Job inputs]" << endl;
  cout << " - Reading analysis configuration from: " << config << endl;
  cout << " - SuperK sample: " << util::AsString(s) << endl;
  cout << " - Toy MC expt info: " << endl;
  cout << "    - Varied systematic parameter (by +/- 1 sigma): " << util::AsString(p) << endl;
  cout << "    - Including statistical fluctuations? : " << fluctuate_stat << endl;
  cout << "    - True oscillation parameters: " << endl;
  cout << "       > " << util::AsString(kPhypMixing12) << " = " << input_mixing12 << endl;
  cout << "       > " << util::AsString(kPhypMixing13) << " = " << input_mixing13 << endl;
  cout << "       > " << util::AsString(kPhypMixing23) << " = " << input_mixing23 << endl;
  cout << "       > " << util::AsString(kPhypCPV    )  << " = " << input_deltacp << endl;
  cout << "       > " << util::AsString(kPhypDm12sq)   << " = " << input_dm12sq << " eV^2/c^4" << endl;
  cout << "       > " << util::AsString(kPhypDm23sq)   << " = " << input_dm23sq << " eV^2/c^4" << endl;
  cout << " - Fit info: " << endl;
  cout << "    - Fit parameter list: " << fit_param_list << endl;
  cout << "    - Fitter option: " << fit_option << endl;
  cout << *Opt::Instance() << endl;

  SuperkObservation * skobs = 0;
  FitInputs * fit_inputs = 0;
  FitterABC * fitter = 0;

  //initialise output
  double bf_variations[kNuis000][2] = {}; // Array big enough to store values of all physics parameters, +1 and -1 sigma
  
  if(throw_n_toys <= 0) {
    
    for (int isigma = 0 ; isigma < 2 ; isigma++) {
      //
      //
      // +/- 1 sigma (relative to tuned BANFF)
      //
      //
      
      double tweak_size = +1;
      if (isigma) {
        tweak_size = -1;
      }
      
      cout << "\n\nApplying "<<tweak_size<<" sigma tweak" << endl;
      if (util::IsBANFF(p)) {
        //find parameter to tweak
        int ip = -1;
        if (util::IsRHCSample(s)) { ip = util::ArrIdx::Instance()->RHCBANFFSystIdx(p);}
        else                      { ip = util::ArrIdx::Instance()->BANFFSystIdx   (p);}
        assert(ip!=-1);
        //initialise tweak vector
        const TMatrixD & cov = skexp->GetBANFF()->SampleCovMatrix();
        const int nbanff = cov.GetNrows();
        TVectorD twkv(nbanff);
        //tweak the parameter to tweak
        twkv[ip] = skexp->GetBANFF()->BestFitValue(p) + (tweak_size * sqrt(cov[ip][ip]));
        //if we're tweaking the other parameters in a correlated way too
        double sigma = (twkv[ip] - skexp->GetBANFF()->BestFitValue(p)) / cov[ip][ip];
        for(int ipp = 0; ipp < nbanff; ipp++) {
          int position = first_banff + ipp;
          PhysParamId_t pp = (PhysParamId_t)position;
          if(pp == p)
            continue;
          //set non-tweaked values to best fit
          twkv[ipp] = skexp->GetBANFF()->BestFitValue(pp);
          //now, tweak the other parameters according to the correlations between them and the fitted parameters
          if(correlated_parameters)
            twkv[ipp] += sigma * cov[ip][ipp];
        }
        skexp->GetBANFF()->TweakSample(true,&twkv);
      } //isBANFF
      else if (util::IsEnergyScale(p)) {
        skexp->TweakSystematic(p, tweak_size);
      }
      else if (util::IsSKDetFSI(p)) {
        int ip = -1;
        if (util::IsRHCSample(s)) { ip = util::ArrIdx::Instance()->RHCSKDetFSISystIdx(p);}
        else                      { ip = util::ArrIdx::Instance()->SKDetFSISystIdx   (p);}
        assert(ip!=-1);
        //initialise tweak vector
        const TMatrixD & cov = cov_skdetfsi;
        const int nskdetfsi = cov.GetNrows();
        TVectorD twkv(nskdetfsi);
        //tweak the parameter to tweak
        twkv[ip] = 1 + (tweak_size * sqrt(cov[ip][ip]));
        //if we're tweaking the other parameters in a correlated way too
        double sigma = (twkv[ip] - 1) / cov[ip][ip];
        for(int ipp = 0; ipp < nskdetfsi; ipp++) {
          int position = first_skdet_fsi + ipp;
          PhysParamId_t pp = (PhysParamId_t)position;
          if(pp == p)
            continue;
          //set non-tweaked values to best fit
          twkv[ipp] = 1;
          //now, tweak the other parameters according to the correlations between them and the fitted parameters
          if(correlated_parameters)
            twkv[ipp] += sigma * cov[ip][ipp];
        }
        skexp->GetSuperkDetectorFSI()->Tweak(true,&twkv);
      } //isSKDetFSI
      else if (!util::IsSystematic(p)) {
        double input_val = fparm->InitValue(p);
        double tweak  = input_val * (1 + tweak_size*Opt::Instance()->GetGlobalFitError(p));
        if(tweak >  fparm->MaxValue(p)) { tweak = fparm->MaxValue(p); }
        if(tweak <  fparm->MinValue(p)) { tweak = fparm->MinValue(p); }
        fparm->SetInitValue(p, tweak);
      }
      else if(do_nominal) {
        //do nothing
      }
      else {
        cerr<<"test_bestfit_point_shift_due_to_syst_twk >> Unrecognised systematic: "<<util::AsString(p)<<std::endl;
      }

      // Generate toy MC
      cout << "Generating toy-MC..." << endl;
      cout << "p.d.f parameters: " << endl;
      skexp->PrintSystematicVariations();
      skobs = new SuperkObservation;
      skobs->CreateMockData(*skexp, PhysHypothesis(*fparm), fluctuate_stat);
      
      // Reset oscillation parameters before fitting
      if(hyp == kHypOsc3F) {
        fparm->SetInitValue(kPhypMixing12, input_mixing12);
        fparm->SetInitValue(kPhypMixing13, input_mixing13);
        fparm->SetInitValue(kPhypMixing23, input_mixing23);
        fparm->SetInitValue(kPhypCPV    ,  input_deltacp);
        fparm->SetInitValue(kPhypDm12sq,   input_dm12sq);
        fparm->SetInitValue(kPhypDm23sq,   input_dm23sq);
      }
      else if(hyp == kHypOsc3FBar) {
        fparm->SetInitValue(kPhypMixing12   , input_mixing12);
        fparm->SetInitValue(kPhypMixing13   , input_mixing13);
        fparm->SetInitValue(kPhypMixing23   , input_mixing23);
        fparm->SetInitValue(kPhypCPV        ,  input_deltacp);
        fparm->SetInitValue(kPhypDm12sq     ,   input_dm12sq);
        fparm->SetInitValue(kPhypDm23sq     ,   input_dm23sq);
        fparm->SetInitValue(kPhypMixing12Bar, input_mixing12);
        fparm->SetInitValue(kPhypMixing13Bar, input_mixing13);
        fparm->SetInitValue(kPhypMixing23Bar, input_mixing23);
        fparm->SetInitValue(kPhypCPVBar     ,  input_deltacp);
        fparm->SetInitValue(kPhypDm12sqBar  ,   input_dm12sq);
        fparm->SetInitValue(kPhypDm23sqBar  ,   input_dm23sq);
      }
      else if(hyp == kHypNuMuDisapp2F) {
        fparm->SetInitValue(kPhypMixing23, input_mixing23);
        fparm->SetInitValue(kPhypDm23sq,   input_dm23sq);
      }  
      else {
        return;
      }

      // Reset systematic variation
      skexp->ResetSystematicVariations();

      // Get fitter
      fit_inputs = new FitInputs(skobs,skexp,fparm);
      fitter = new FitterMinuit(fit_inputs,fit_option);

      // Do the fit
      fitter->DoTheFit();
      
      // Save the bf variations
      for (int ip = 0 ; ip < kNuis000 ; ip++) {
        PhysParamId_t pp = (PhysParamId_t) ip;
        if (fparm->Exists(pp)) {
          bf_variations[ip][isigma] = fparm->BestFitValue(pp);
          //~ cerr<<util::AsString(pp)<<"\t"<<fparm->BestFitValue(pp)<<endl;
        }
      }//ip

      // Clean-up in preparation for the next fit
      delete skobs;
      delete fit_inputs;
      delete fitter;
      skexp->ResetSystematicVariations();
      fparm->ResetCurrToInitValues();
      //need to do this, just incase InitValues were changed inside FitterKernel::DoTheFit()
      if(hyp == kHypOsc3F) {
        fparm->SetInitValue(kPhypMixing12, input_mixing12);
        fparm->SetInitValue(kPhypMixing13, input_mixing13);
        fparm->SetInitValue(kPhypMixing23, input_mixing23);
        fparm->SetInitValue(kPhypCPV    ,  input_deltacp);
        fparm->SetInitValue(kPhypDm12sq,   input_dm12sq);
        fparm->SetInitValue(kPhypDm23sq,   input_dm23sq);
      }
      else if(hyp == kHypOsc3FBar) {
        fparm->SetInitValue(kPhypMixing12   , input_mixing12);
        fparm->SetInitValue(kPhypMixing13   , input_mixing13);
        fparm->SetInitValue(kPhypMixing23   , input_mixing23);
        fparm->SetInitValue(kPhypCPV        ,  input_deltacp);
        fparm->SetInitValue(kPhypDm12sq     ,   input_dm12sq);
        fparm->SetInitValue(kPhypDm23sq     ,   input_dm23sq);
        fparm->SetInitValue(kPhypMixing12Bar, input_mixing12);
        fparm->SetInitValue(kPhypMixing13Bar, input_mixing13);
        fparm->SetInitValue(kPhypMixing23Bar, input_mixing23);
        fparm->SetInitValue(kPhypCPVBar     ,  input_deltacp);
        fparm->SetInitValue(kPhypDm12sqBar  ,   input_dm12sq);
        fparm->SetInitValue(kPhypDm23sqBar  ,   input_dm23sq);
      }
      else if(hyp == kHypNuMuDisapp2F) {
        fparm->SetInitValue(kPhypMixing23, input_mixing23);
        fparm->SetInitValue(kPhypDm23sq,   input_dm23sq);
      }  
    } //isigma
    
  }//throw_n_toys <= 0
  else {

    //initialise
    double bf_toys[kNuis000] = {}; // big enough to hold any physics parameters

    for(int itoy = 0; itoy < throw_n_toys; itoy++) {
      cout << "\nToy " << itoy << endl;

      // Reset systematic variation
      skexp->ResetSystematicVariations();

      skobs = new SuperkObservation;

      //apply the tweak
      bool applied_negevents_condition = true;
      while(applied_negevents_condition) {
        if(util::IsBANFF(p)) {
          int ip = -1;
          if (util::IsRHCSample(s)) { ip = util::ArrIdx::Instance()->RHCBANFFSystIdx(p);}
          else                      { ip = util::ArrIdx::Instance()->BANFFSystIdx   (p);}
          assert(ip!=-1);
          //initialise tweak vector
          const TMatrixD & cov = skexp->GetBANFF()->SampleCovMatrix();
          const int nbanff = cov.GetNrows();
          TVectorD twkv(nbanff);
          
          //tweak the parameter to tweak
          cout << *Rnd::Instance() << endl;
          int seed = -1;
          twkv[ip] = Rnd::Instance()->RndOther().Gaus(skexp->GetBANFF()->BestFitValue(p),sqrt(cov[ip][ip]));
          //if we're tweaking the other parameters in a correlated way too
          double sigma = (twkv[ip] - skexp->GetBANFF()->BestFitValue(p)) / cov[ip][ip];
          for(int ipp = 0; ipp < nbanff; ipp++) {
            int position = first_banff + ipp;
            PhysParamId_t pp = (PhysParamId_t)position;
            if(pp == p)
              continue;
            //set non-tweaked values to best fit
            twkv[ipp] = skexp->GetBANFF()->BestFitValue(pp);
            //now, tweak the other parameters according to the correlations between them and the fitted parameters
            if(correlated_parameters)
              twkv[ipp] += sigma * cov[ip][ipp];
          }
          skexp->GetBANFF()->TweakSample(true,&twkv);
        }
        else if (util::IsEnergyScale(p)) {
          double tweak_size = Rnd::Instance()->RndOther().Gaus(0,1);
          skexp->TweakSystematic(p, tweak_size);
        }
        else if(util::IsSKDetFSI(p)) {
          int ip = -1;
          if (util::IsRHCSample(s)) { ip = util::ArrIdx::Instance()->RHCSKDetFSISystIdx(p);}
          else                      { ip = util::ArrIdx::Instance()->SKDetFSISystIdx   (p);}
          assert(ip!=-1);
          //initialise tweak vector
          const TMatrixD & cov = cov_skdetfsi;
          const int nskdetfsi = cov.GetNrows();
          TVectorD twkv(nskdetfsi);
          //tweak the parameter to tweak
          twkv[ip] = Rnd::Instance()->RndOther().Gaus(1,sqrt(cov[ip][ip]));
          //if we're tweaking the other parameters in a correlated way too
          double sigma = (twkv[ip] - 1) / cov[ip][ip];
          for(int ipp = 0; ipp < nskdetfsi; ipp++) {
            int position = first_skdet_fsi + ipp;
            PhysParamId_t pp = (PhysParamId_t)position;
            if(pp == p)
              continue;
            //set non-tweaked values to best fit
            twkv[ipp] = 1;
            //now, tweak the other parameters according to the correlations between them and the fitted parameters
            if(correlated_parameters)
              twkv[ipp] += sigma * cov[ip][ipp];
          }
          skexp->GetSuperkDetectorFSI()->Tweak(true,&twkv);
        }
        else if (!util::IsSystematic(p)) {
          double tweak  = util::RandomizeParameter(p);
          if(tweak >  fparm->MaxValue(p)) { tweak = fparm->MaxValue(p); }
          if(tweak <  fparm->MinValue(p)) { tweak = fparm->MinValue(p); }
          fparm->SetInitValue(p, tweak);
        }
        else if(do_nominal) {
          //do nothing
        }
        else {
          cerr<<"test_bestfit_point_shift_due_to_syst_twk >> Unrecognised systematic: "<<util::AsString(p)<<std::endl;
        }
        // Generate toy MC
        cout << "Generating toy-MC..." << endl;
        cout << "p.d.f parameters: " << endl;
        skexp->PrintSystematicVariations();

        applied_negevents_condition = skobs->CreateMockData(*skexp, PhysHypothesis(*fparm), fluctuate_stat);
      }//applied_negevents_condition
     
      // Reset oscillation parameters before fitting
      if(hyp == kHypOsc3F) {
        fparm->SetInitValue(kPhypMixing12, input_mixing12);
        fparm->SetInitValue(kPhypMixing13, input_mixing13);
        fparm->SetInitValue(kPhypMixing23, input_mixing23);
        fparm->SetInitValue(kPhypCPV    ,  input_deltacp);
        fparm->SetInitValue(kPhypDm12sq,   input_dm12sq);
        fparm->SetInitValue(kPhypDm23sq,   input_dm23sq);
      }
      else if(hyp == kHypOsc3FBar) {
        fparm->SetInitValue(kPhypMixing12   , input_mixing12);
        fparm->SetInitValue(kPhypMixing13   , input_mixing13);
        fparm->SetInitValue(kPhypMixing23   , input_mixing23);
        fparm->SetInitValue(kPhypCPV        ,  input_deltacp);
        fparm->SetInitValue(kPhypDm12sq     ,   input_dm12sq);
        fparm->SetInitValue(kPhypDm23sq     ,   input_dm23sq);
        fparm->SetInitValue(kPhypMixing12Bar, input_mixing12);
        fparm->SetInitValue(kPhypMixing13Bar, input_mixing13);
        fparm->SetInitValue(kPhypMixing23Bar, input_mixing23);
        fparm->SetInitValue(kPhypCPVBar     ,  input_deltacp);
        fparm->SetInitValue(kPhypDm12sqBar  ,   input_dm12sq);
        fparm->SetInitValue(kPhypDm23sqBar  ,   input_dm23sq);
      }
      else if(hyp == kHypNuMuDisapp2F) {
        fparm->SetInitValue(kPhypMixing23, input_mixing23);
        fparm->SetInitValue(kPhypDm23sq,   input_dm23sq);
      }  
      else {
        return;
      }


      // Reset systematic variation
      skexp->ResetSystematicVariations();
      
      // Get fitter
      fit_inputs = new FitInputs(skobs,skexp,fparm);
      fitter = new FitterMinuit(fit_inputs,fit_option);

      // Do the fit
      fitter->DoTheFit();

      // Save the bf variations
      for (int ip = 0 ; ip < kNuis000 ; ip++) {
        PhysParamId_t pp = (PhysParamId_t) ip;
        if (fparm->Exists(pp)) {
          bf_toys[ip] += fparm->BestFitValue(pp);
        }
      }//ip

      //clean
      delete skobs;
      delete fit_inputs;
      delete fitter;

      //reset osc params to initial
      fparm->ResetCurrToInitValues();
      //need to do this, just incase InitValues were changed inside FitterKernel::DoTheFit()
      if(hyp == kHypOsc3F) {
        fparm->SetInitValue(kPhypMixing12, input_mixing12);
        fparm->SetInitValue(kPhypMixing13, input_mixing13);
        fparm->SetInitValue(kPhypMixing23, input_mixing23);
        fparm->SetInitValue(kPhypCPV    ,  input_deltacp);
        fparm->SetInitValue(kPhypDm12sq,   input_dm12sq);
        fparm->SetInitValue(kPhypDm23sq,   input_dm23sq);
      }
      else if(hyp == kHypOsc3FBar) {
        fparm->SetInitValue(kPhypMixing12   , input_mixing12);
        fparm->SetInitValue(kPhypMixing13   , input_mixing13);
        fparm->SetInitValue(kPhypMixing23   , input_mixing23);
        fparm->SetInitValue(kPhypCPV        ,  input_deltacp);
        fparm->SetInitValue(kPhypDm12sq     ,   input_dm12sq);
        fparm->SetInitValue(kPhypDm23sq     ,   input_dm23sq);
        fparm->SetInitValue(kPhypMixing12Bar, input_mixing12);
        fparm->SetInitValue(kPhypMixing13Bar, input_mixing13);
        fparm->SetInitValue(kPhypMixing23Bar, input_mixing23);
        fparm->SetInitValue(kPhypCPVBar     ,  input_deltacp);
        fparm->SetInitValue(kPhypDm12sqBar  ,   input_dm12sq);
        fparm->SetInitValue(kPhypDm23sqBar  ,   input_dm23sq);
      }
      else if(hyp == kHypNuMuDisapp2F) {
        fparm->SetInitValue(kPhypMixing23, input_mixing23);
        fparm->SetInitValue(kPhypDm23sq,   input_dm23sq);
      }  
    }//end loop over itoy
    
    // Save the bf variations
    for (int ip = 0 ; ip < kNuis000 ; ip++) {
      PhysParamId_t pp = (PhysParamId_t) ip;
      if (fparm->Exists(pp)) {
        for (int isigma = 0 ; isigma < 2 ; isigma++) {
          bf_variations[ip][isigma] = bf_toys[ip] / throw_n_toys;
        }
      }
    }//ip
    
  }//throw_n_toys > 0

  delete skexp;
  
  string simple_p_name = (do_nominal ? "nominal" : util::AsStringSimple(p));
  string tex_p_name    = (do_nominal ? "nominal" : util::AsStringPaperTex(p));
  
  // Write-out results
  ostringstream outfilename;
  string cor = "";
  if(correlated_parameters)
    cor = "_correlated";
  outfilename << "star_plot_line-" << cor << simple_p_name << ".data";
  ofstream outfile;
  outfile.open(outfilename.str().c_str());
  outfile << tex_p_name;
  for (int isigma = 0 ; isigma < 2 ; isigma++) {
    for (int ip = 0 ; ip < kNuis000 ; ip++) {
      PhysParamId_t pp = (PhysParamId_t) ip;
      if (fparm->Exists(pp)) {
        outfile << "\t" << bf_variations[ip][isigma];
      }
    }
  }
  outfile << endl;
  outfile.close();
  
  cout<<"\n\n========================\n\n"<<endl;
  for (int ip = 0 ; ip < kNuis000 ; ip++) {
    PhysParamId_t pp = (PhysParamId_t) ip;
    if (fparm->Exists(pp) && (!fparm->IsFixed(pp))) {
      cout << util::PadToLength(util::AsString(pp),20) << "\t" << bf_variations[ip][0] << "\t"<<bf_variations[ip][1]<<endl;
    }
  }
  cout << "\nDone!" << endl;
}










#ifndef __CINT__
int main()
{
  valor::Env::Instance()->AutoSet("t2xk");
  TString conf_dir  = valor::Env::Instance()->TopLevelConfDir();
  TString verbosity = conf_dir + "/verbosity.xml";
  valor::util::Verbosity(verbosity.Data());
  
  //~ test_bestfit_point_shift_due_to_syst_twk(valor::conv::kRHCNuisBANFFMAQE,"RHC1Rmu",
    //~ "nuosc3fbar_numudisapp_fixnonbar_fix4+nuisance","minuit",
    //~ -99,-99,-99,-99,-99,-99,
    //~ "analysis_configuration_2.5E20.xml",true,5);
  test_bestfit_point_shift_due_to_syst_twk("f_skdetfsi000_rhc","RHC1Rmu",
    "nuosc3fbar_numudisapp_fixnonbar_fix4+nuisance","minuit",
    -99,-99,-99,-99,-99,-99,
    "analysis_configuration_2.5E20.xml",false,02  );
  return 0;
}
#endif
