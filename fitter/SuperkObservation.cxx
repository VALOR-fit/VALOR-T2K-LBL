//________________________________________________________________________________________
/*

 >> ChangeLog:

 For changes up to version 2010_v3.0.0 (last implementation of the 2010-era analysis)
 see code in the corresponding SVN branch.

 Changes since 2011_v0.9.0:

 @18/2/13 LE
 Read Ereco spectrum for 1Re.
 Doubled arrays (independent ones for 1Rmu and 1Re)

 31-12-2013  Nick
 Added samples for reverse horn current running

 21-7-2014 Davide
 Added samples for multiring forward horn current

*/
//________________________________________________________________________________________

//#define __show_num_of_events_mesg_debug__
//#define __show_bin_contributions_to_chisq__
//#define __announce_generate__

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <sstream>

#include <TFile.h>
#include <TObjArray.h>
#include <TXMLEngine.h>
#include <TString.h>
#include <TAxis.h>
#include <TH1D.h>
#include <TRandom3.h>
#include <TMath.h>
#include <TCanvas.h>

#include "Enumerations.h"
#include "AnalysisOptions.h"

#include "SuperkObservation.h"
#include "Utilities.h"
#include "Utilities/Rnd.h"

using std::istringstream;
using std::ostringstream;
using std::cout;
using std::cerr;
using std::endl;

using namespace valor;
using namespace valor::conv;

ClassImp(SuperkObservation)

//________________________________________________________________________________________
SuperkObservation::SuperkObservation(void) 
{
  this->Init();
}
//________________________________________________________________________________________
SuperkObservation::~SuperkObservation(void)
{
  this->CleanUp();
}
//________________________________________________________________________________________
void SuperkObservation::Init(void)
{
// Initialize private data members
//
  for (int is = 0 ; is < conv::kns ;  is ++ ) {
    SuperkSample_t s = kSKSamples[is];
    if (s==kSKsUndefined) {
      fErDistrArray[is]=0;
    }
    else {
      fErDistrArray[is] = new VALORArray(util::NumErecoBins(s));
      fErDistrArray[is]->SetTitle(Form("Ereco distribution: %s",util::AsString(s).c_str()));
      fErDistrArray[is]->SetName(Form("%s",util::AsStringCompact(s).c_str()));
      fErDistrArray[is]->SetDimensionName(0,"Reconstructed energy /GeV");
    }
  }
  
  arr_index = util::ArrIdx::Instance();
}
//________________________________________________________________________________________
void SuperkObservation::CleanUp(void)
{
  for (int is = 0 ; is < conv::kns ; is++) {
    if (fErDistrArray[is]) delete fErDistrArray[is];
  }
}
//________________________________________________________________________________________
bool SuperkObservation::Read(string config_file, Option_t * opt)
{
  int data_set_id = 0;
  bool is_xml = true;
  if ( opt ) {
    string stopt(opt);
    size_t pos =stopt.find("data_set=");
    if (pos != string::npos) {
      istringstream iss(stopt.c_str()+pos+9);
      iss >> data_set_id;
    }
    if (stopt.find("root")!=string::npos) {
      is_xml = false;
    }
  }
  return Read(config_file,data_set_id,is_xml);
}
//________________________________________________________________________________________
bool SuperkObservation::Read(string config_file, int data_set_id, bool is_xml)
{
  cout << "Obs >> " << "Reading data file: " << config_file << endl;

  if(is_xml) { 
    // Parse the input configuration file
    TXMLEngine* xml = new TXMLEngine;
    XMLDocPointer_t xmldoc = xml->ParseFile(config_file.c_str());
    if (xmldoc==0) {
      cerr << "SKObs/Read >> " << "Can not parse input: " << config_file << endl;
      delete xml;
      std::exit(1);
    }
    XMLNodePointer_t root_node = xml->DocGetRootElement(xmldoc);
    
    for (int is = 0 ; is < conv::kns ; is++) {
      SuperkSample_t s = conv::kSKSamples[is];
      if (s==kSKsUndefined) continue;
      ostringstream datanode;
      datanode << "superk/observation";
      datanode << "/data_set_" << data_set_id;
      datanode << "/"<<util::GetXMLDatasetName(s);
      datanode << "/ereco_values"; 

      if(util::NodeExists(xml,root_node,datanode.str())) {
        TH1D h("", "", util::NumErecoBins(s), util::ErecoBinEdges(s));
        h.SetDirectory(0);
        vector<double> values = util::XML2DblArray(xml,root_node,datanode.str());
        vector<double>::const_iterator value_iter = values.begin();
        for( ; value_iter != values.end(); ++value_iter) {
          double Ereco = *value_iter;
          h.Fill(Ereco);
        }
        cout << "Obs >> Read: " << h.Integral() 
             << " " << util::AsString(s) << " events." << endl;
        for(int ibin = 1; ibin <= util::NumErecoBins(s); ibin++) {
          double nobs = h.GetBinContent(ibin);
          fErDistrArray[is]->GetArray()[ibin-1] = nobs;
        }
      }
      else {
        cerr << "Obs >> Node " << datanode.str() << " does not exist!" << endl;
      }
    } // is
    
    // Release memory
    xml->FreeDoc(xmldoc);
    delete xml;
  }//is_xml

  else { // DEPRECATE?
   
    //assert(!"This code branch is deprecated");


    // NEW FHC+RHC JOINT

#if defined(__RHC_RUNNING__) && defined(__RHC_RUNNING__) // joint fit
    
    TFile f(config_file.c_str());
    TH1D * hist = 0;
 
    for (int is = 0 ; is < conv::kns ; is++) {
      SuperkSample_t s = conv::kSKSamples[is];
      
      if (s == conv::kSKsUndefined) continue;
      else if(s == conv::kSKsRHCFCFV1RingMuLike) f.GetObject(Form("mock_data_1Rmu_RHC_%d",data_set_id),hist);
      else if(s == conv::kSKsRHCFCFV1RingELike) f.GetObject(Form("mock_data_1Re_RHC_%d",data_set_id),hist);
      else if(s == conv::kSKsFCFV1RingELike) f.GetObject(Form("mock_data_1Re_%d",data_set_id),hist);
      else if(s == conv::kSKsFCFV1RingMuLike) f.GetObject(Form("mock_data_1Rmu_%d",data_set_id),hist);
      else {
	cerr << endl;
	cerr << "Used SK sample " << util::AsString(s) << " not available" << endl;
	cerr << endl;
	std::exit(1);
      }
      
      if(!hist) {
	cerr << "Obs >> Histogram for " << util::AsString(s) << " does not exist!" << endl;
	std::exit(1);
      }

      for(int ibin = 1; ibin <= util::NumErecoBins(s); ibin++) {
	double nobs = hist->GetBinContent(ibin);
	fErDistrArray[is]->GetArray()[ibin-1] = nobs;
      }
      cout << "Obs >> Read: " << hist->Integral()
     	   << " " << util::AsString(s) << " events." << endl;
    }

#elif defined(__RHC_RUNNING__) // rhc only analyses (no joint fit)
    
    // nuebar app only
#if defined(__NUEONLY__) 
    TFile f(config_file.c_str());
    TH1D * hist = 0;
    f.GetObject(Form("mock_data_1Re_RHC_%d",data_set_id),hist);
    if(!hist) {
      cerr << "Obs >> Histogram for RHC1Re does not exist!" << endl;
      std::exit(1);
    }
 
    for (int is = 0 ; is < conv::kns ; is++) {
      SuperkSample_t s = conv::kSKSamples[is];
      if (s!=kSKsRHCFCFV1RingELike) continue;
      
      for(int ibin = 1; ibin <= util::NumErecoBins(s); ibin++) {
	double nobs = hist->GetBinContent(ibin);
	fErDistrArray[is]->GetArray()[ibin-1] = nobs;
      }
      cout << "Obs >> Read: " << hist->Integral()
     	   << " " << util::AsString(s) << " events." << endl;
    }

    // numubar disap only
#elif defined(__NUMUONLY__) 

    TFile f(config_file.c_str());
    TH1D * hist = 0;
    f.GetObject(Form("mock_data_1Rmu_RHC_%d",data_set_id),hist);
    if(!hist) {
      cerr << "Obs >> Histogram for RHC1Rmu does not exist!" << endl;
      std::exit(1);
    }
 
    for (int is = 0 ; is < conv::kns ; is++) {
      SuperkSample_t s = conv::kSKSamples[is];
      if (s!=kSKsRHCFCFV1RingMuLike) continue;
      
      for(int ibin = 1; ibin <= util::NumErecoBins(s); ibin++) {
	double nobs = hist->GetBinContent(ibin);
	fErDistrArray[is]->GetArray()[ibin-1] = nobs;
      }
      cout << "Obs >> Read: " << hist->Integral()
     	   << " " << util::AsString(s) << " events." << endl;
    }

#endif // __NUEONLY__ - __NUMUONLY__

#else
    cerr << endl;
    cerr << "SuperkObservation::Read()" << endl;
    cerr << "Analysis mode not defined" << endl;
    cerr << endl;
    std::exit(1);
#endif // __RHC_RUNNING__ - __FHC_RUNNING__

    //





    
    // else {
    //   assert(hist->GetNbinsX() == util::NumErecoBins(kSKsRHCFCFV1RingELike));
    //   double integral = 0;
    //   for(int ix = 1; ix <= hist->GetNbinsX(); ix++) {
    // 	double nev   = hist->GetBinContent(ix);
    // 	integral+=nev;
    // 	fErDistrArray[is]->GetArray()[ix-1] = nev;
    //   }
    //   cout << "Obs >> Read: " << integral
    // 	   << " " << util::AsString(kSKsRHCFCFV1RingELike) << " events." << endl;
    // }
    //#endif
    
  }//!is_xml

  return true;
}
//________________________________________________________________________________________
double SuperkObservation::NumOfEvents(SuperkSample_t s) const
{
  double nev = 0;

  if(s == kSKsUndefined) {
    return nev;
  }

  int is = arr_index->SuperkSampleIdx(s);
  const double * spectrum = this->ErecoSpectrumArray(s);
  
  for(int ir = 0; ir < kNErecoBins[is]; ir++){
    nev += spectrum[ir];
  }

#ifdef __show_num_of_events_mesg_debug__
  cout << "Obs >> NEvents(" 
       << util::AsString(s) << ") = " << nev << endl;
#endif

  return nev;
}
//________________________________________________________________________________________
double SuperkObservation::NumOfEvents(double Er, SuperkSample_t s) const
{
  double nev = 0;

  if(s == kSKsUndefined) {
    return 0;
  }
  
  cerr<<"WARNING: SuperkObservation::NumOfEvents(Ereco) needs testing"<<std::endl;
  // NEEDS TESTING
  const double * spectrum = this->ErecoSpectrumArray(s);
  
  int bin = util::ErecoBinContains(s,Er);
  assert(bin!=-1);
  nev = spectrum[util::ErecoBinContains(s,Er)];
  
#ifdef __show_num_of_events_mesg_debug__
  cout << "Obs >> NEvents(" 
       << util::AsString(s) << ", Ereco = " << Er << " GeV) = " << nev << endl;
#endif

  return nev;
}
//________________________________________________________________________________________
int SuperkObservation::NumOfBins()
{
  int nbins = 0;
  for (int is = 0 ; is < kns ; is++) {
    SuperkSample_t s = kSKSamples[is];
    if (s==kSKsUndefined) continue;
    nbins+=util::NumErecoBins(s);
  }
  return nbins;
}
//________________________________________________________________________________________
TH1D* SuperkObservation::ErecoSpectrum (valor::conv::SuperkSample_t s) const
{
  if(s == kSKsUndefined) {
    return 0;
  }
  
  TH1D * data = new TH1D("","",util::NumErecoBins(s),util::ErecoBinEdges(s));
  util::FormatErecoHist(data);
  data -> SetDirectory(0);
  data -> SetMarkerStyle(20);
  data -> SetMarkerSize(1);
  data -> GetXaxis()->SetTitle("E_{reco} (GeV)");
  data -> GetYaxis()->SetTitle(Form("Number of %s events",util::AsString(s).c_str()));
  
  const double * spectrum = this->ErecoSpectrumArray(s);
  
  for(int ir = 1; ir <= data->GetNbinsX(); ir++) {
    data->SetBinContent(ir, spectrum[ir-1]);
  }//ir
  
  return data;
}
//________________________________________________________________________________________
VALORArray * SuperkObservation::ErecoSpectrumVArray (valor::conv::SuperkSample_t s)
{
  return fErDistrArray[s];
}
//________________________________________________________________________________________
const VALORArray * SuperkObservation::ErecoSpectrumVArray (valor::conv::SuperkSample_t s) const
{
  return fErDistrArray[s];
}
//________________________________________________________________________________________
double* SuperkObservation::ErecoSpectrumArray (valor::conv::SuperkSample_t s)
{
  VALORArray * va = ErecoSpectrumVArray(s);
  if (not va) return 0;
  return ErecoSpectrumVArray(s)->GetArray();
}
//________________________________________________________________________________________
const double* SuperkObservation::ErecoSpectrumArray (valor::conv::SuperkSample_t s) const
{
  const VALORArray * va = ErecoSpectrumVArray(s);
  if (not va) return 0;
  return ErecoSpectrumVArray(s)->GetArray();
}
//________________________________________________________________________________________
bool SuperkObservation::Generate(
     FitModelABC * model, bool fluctuate)
{
  SuperkExpectation * skexp = dynamic_cast<SuperkExpectation*>(model);
  
  assert(skexp);
 
  for (int is = 0 ; is < conv::kns ; is++ ) {
    SuperkSample_t s = kSKSamples[is];
    if (s==kSKsUndefined) continue;
    for(int ibin = 0; ibin < util::NumErecoBins(s); ibin++) {
       double nexp  = skexp->CurrPredictedErecoSpectrumArray(s)[ibin];
       double nobs  = (fluctuate) ? Rnd::Instance()->RndToyMCstat().Poisson(nexp) : nexp;
       this->ErecoSpectrumArray(s)[ibin] = nobs;
    }
  }

  return false;
}
//________________________________________________________________________________________
bool SuperkObservation::Generate(
      FitModelABC * model, const FitParamList & fp, bool fluctuate)
{
// Get observed reconstructed energy distribution
  
  fp.ShowOnlyCurrValuesInNextPrint();
#ifdef __announce_generate__
  cout << "Obs >> Generating a toy-MC expt with the following fitparamlist: " << fp;
  cout << "Obs >> Stat. fluctuations are turned-" << ((fluctuate)? "ON" : "OFF") << endl;
#endif
  
  // Should do this properly.
  SuperkExpectation * skexp = dynamic_cast<SuperkExpectation*>(model);
  
  assert(skexp);
  
  bool applied_negevents_condition = skexp->Generate(fp); //generate prediction
  
  applied_negevents_condition = this->Generate(model,fluctuate) || applied_negevents_condition;
  return applied_negevents_condition;
}
//________________________________________________________________________________________
bool SuperkObservation::Generate(
     FitModelABC * model, const PhysHypothesis & ph, bool fluctuate)
{

// Get observed reconstructed energy distribution
#ifdef __announce_generate__
  cout << "Obs >> Generating a toy-MC expt with the following physics hypothesis: " << ph;
  cout << "Obs >> Stat. fluctuations are turned-" << ((fluctuate)? "ON" : "OFF") << endl;
#endif
  
  // Should do this properly.
  SuperkExpectation * skexp = dynamic_cast<SuperkExpectation*>(model);

  assert(skexp);
  
  bool applied_negevents_condition = skexp->Generate(ph); //generate prediction
  applied_negevents_condition = this->Generate(model,fluctuate) || applied_negevents_condition;

  return applied_negevents_condition;
}
//________________________________________________________________________________________
bool SuperkObservation::Write(TObjArray * data_collection, const char * data_tag, Option_t * write_opt)
{
  if(!data_collection) {
    cout << "Obs >> Null TObjArray. Not writing data" << endl;
    return false;
  }
  
  bool hist_mode = (write_opt && string(write_opt).find("hist")!=string::npos);

  for(int is=0; is<kns; is++) {
    SuperkSample_t s = kSKSamples[is];
    if(s == kSKsUndefined) continue;
    
    if (hist_mode) {
      util::SetDefaultStyle();
      // This is a new object, leaks memory if data_collection has not been set as owner
      TH1D * data = this->ErecoSpectrum(s);
      
      data -> SetName( 
       Form( "%s_%s", 
              data_tag,
              util::AsStringCompact(s).c_str() ) 
      );
      data_collection -> Add(data);
    }
    else {
      VALORArray * data = this->ErecoSpectrumVArray(s);
      data -> SetName( Form( "%s_%s", data_tag, util::AsStringCompact(s).c_str() ) );
      data_collection -> Add(data);
    }

  }
  
  return true;
}
//________________________________________________________________________________________
bool SuperkObservation::Write(
  const char * filename, const char * data_tag, Option_t * write_opt, Option_t * file_opt)
{
  if (!file_opt) file_opt="update";
  TObjArray * array = new TObjArray();
  if (write_opt && string(write_opt).find("hist")!=string::npos) array->SetOwner(true);
  this->Write(array, data_tag, write_opt);

  TFile f(filename, file_opt);
  array->Write();
  f.Close();

  return true;
}
//________________________________________________________________________________________
double SuperkObservation::Statistic (FitModelABC * model, const FitParamList & fp,   Option_t * stat_opt)
{
  bool applied_negevents_condition = model->Generate(fp);
  double total_chisq = this->Statistic(model,stat_opt);
  double phys_penalty = fp.HypothesisFitStatisticPenalty(stat_opt);
  
  if (applied_negevents_condition) cerr<<"SuperkObservation/Statistic >> Warning: applied negative events condition"<<std::endl;
  
  //~ std::cerr<<"Statistic: "<<total_chisq<<"\t"<<phys_penalty<<endl;
  total_chisq += phys_penalty;
  return total_chisq;
}
//________________________________________________________________________________________
double SuperkObservation::Statistic (FitModelABC * model, const PhysHypothesis & ph, Option_t * stat_opt)
{
  bool applied_negevents_condition = model->Generate(ph);
  double total_chisq = this->Statistic(model,stat_opt);
  if (applied_negevents_condition) cerr<<"SuperkObservation/Statistic >> Warning: applied negative events condition"<<std::endl;
  
  return total_chisq;
}
//________________________________________________________________________________________
double SuperkObservation::Statistic (FitModelABC * model, Option_t * stat_opt)
{
  // We will try casting to either of the SK models (single expectation, or many toys)
  SuperkExpectation * skexp = dynamic_cast<SuperkExpectation*>(model);
  
  //std::cerr<<"Tried casting "<<model->Name()<<" to SuperkExpectation, got "<<skexp<<std::endl;
  if (skexp) {
  
    bool gauss   = false;
    bool gof     = false; //goodness of fit
    if (stat_opt) {
      string stopt(stat_opt);
      if (stopt.find("gauss") != string::npos) gauss = true;
      if (stopt.find("gof")   != string::npos) gof   = true;
    }
    
    double total_chisq = 0;
    
    double nobs,nexp,delta_chisq;
    
    for (int is = 0 ; is < kns ; is++) {
      SuperkSample_t s = kSKSamples[is];
      // current sample is undefined: skip
      if (s == conv::kSKsUndefined) continue;
      int n_ereco_bins = kNErecoBins[is];
      
      const VALORArray * obs_v = this->ErecoSpectrumVArray(s);
      const VALORArray * exp_v = skexp->CurrPredictedErecoVArray(s);
      const double * obs_arr = obs_v->GetArray();
      const double * exp_arr = exp_v->GetArray();
      
      if (not gof) {
          
        for (int ibin = 0 ; ibin < n_ereco_bins ; ibin++) {
          nobs = obs_arr[ibin];
          nexp = exp_arr[ibin];
          
          if    (gauss) { delta_chisq = util::ChisqGauss  (nobs,nexp); }
          else          { delta_chisq = util::ChisqPoisson(nobs,nexp); }
          
          total_chisq += delta_chisq;
        
  #ifdef __show_bin_contributions_to_chisq__
          double xlo = util::ErecoBinEdges(s,ibin);
          double xup = util::ErecoBinEdges(s,ibin+1);
            cout << "SuperkObservation::Statistic >> delta_chisq (s = "<<util::AsString(s)<<", bin = " << ibin 
                 << "; [" << xlo << ", " << xup << "] GeV) = " << delta_chisq << endl;
            cout << "\tNobs: " << nobs << "\tNexp:" << nexp << endl;
  #endif
        }//ibin
      }// not goodness of fit
      else {
        
        const int nbins = kNbinsGoodnessOfFit;
        const double * binedges = kBinEdgesGoodnessOfFit;
      
        //loop over Ereco bins
        for(int ibin = 0; ibin < nbins; ibin++) {

          int ErecoMinBin = TMath::BinarySearch(n_ereco_bins,util::ErecoBinEdges(s),binedges[ibin]);
          int ErecoMaxBin = TMath::BinarySearch(n_ereco_bins,util::ErecoBinEdges(s),binedges[ibin+1]);
          nobs = 0;
          nexp = 0;
          for(int ir = ErecoMinBin; ir <= ErecoMaxBin; ir++) {
            nobs += obs_arr[ir];
            nexp += exp_arr[ir];
          }
          
          if    (gauss) { delta_chisq = util::ChisqGauss  (nobs,nexp); }
          else          { delta_chisq = util::ChisqPoisson(nobs,nexp); }
          
          total_chisq += delta_chisq;
            
  #ifdef __show_bin_contributions_to_chisq__
          double xlo = util::ErecoBinEdges(s,ibin);
          double xup = util::ErecoBinEdges(s,ibin+1);
          cout << "FitInputs::Chisquare()> delta_chisq (s = "<<util::AsString(s)<<", bin = " << ibin 
             << "; [" << xlo << ", " << xup << "] GeV) = " << delta_chisq << endl;
          cout << "\tNobs: " << nobs << "\tNexp:" << nexp << endl;
  #endif
        }//ibin
      } //goodness of fit
    }//is
    
    double syst_chisq_pen = model->GetFitStatisticPenalty();
    total_chisq += syst_chisq_pen;
    
    return total_chisq;
  } // end of SKExp cast
  
  // If we get here, the cast has failed.
  std::cerr<<"SuperkOberservation/Generate >> Failed to interpret model: "<<model->Name()<<std::endl;
  std::exit(-1);
}
//________________________________________________________________________________________
