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

*/
//________________________________________________________________________________________

//#define __show_num_of_events_mesg_debug__
//#define __show_bin_contributions_to_chisq__

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
  fErDistrArray1Rmu    = new VALORArray (conv::kNErecoBins1Rmu    );
  fErDistrArray1Re     = new VALORArray (conv::kNErecoBins1Re     );
  fRHCErDistrArray1Rmu = new VALORArray (conv::kRHCNErecoBins1Rmu );
  fRHCErDistrArray1Re  = new VALORArray (conv::kRHCNErecoBins1Re  );
  
  fErDistrArray1Rmu    ->SetAll(0.);
  fErDistrArray1Re     ->SetAll(0.);
  fRHCErDistrArray1Rmu ->SetAll(0.);
  fRHCErDistrArray1Re  ->SetAll(0.);
  
  arr_index = util::ArrIdx::Instance();
}
//________________________________________________________________________________________
void SuperkObservation::CleanUp(void)
{
  delete fErDistrArray1Rmu    ;
  delete fErDistrArray1Re     ;
  delete fRHCErDistrArray1Rmu ;
  delete fRHCErDistrArray1Re  ;
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
  cout << "Obs >> " << "Reading config file: " << config_file << endl;

  if(is_xml) { 
    // Parse the input configuration file
    TXMLEngine* xml = new TXMLEngine;
    XMLDocPointer_t xmldoc = xml->ParseFile(config_file.c_str());
    if (xmldoc==0) {
      cerr << "Obs >> " << "Can not parse input: " << config_file << endl;
      delete xml;
      exit(1);
    }
    XMLNodePointer_t root_node = xml->DocGetRootElement(xmldoc);
    
#ifdef __FHC_RUNNING__
#ifndef __NUEONLY__ 
    ostringstream datanode1Rmu;
    datanode1Rmu << "superk/observation";
    datanode1Rmu << "/data_set_" << data_set_id;
    datanode1Rmu << "/one_mu_ring/ereco_values"; 

    if(util::NodeExists(xml,root_node,datanode1Rmu.str())) {
      TH1D h("", "", kNErecoBins1Rmu,kErecoBinEdges1Rmu);
      h.SetDirectory(0);
      vector<double> values = util::XML2DblArray(xml,root_node,datanode1Rmu.str());
      vector<double>::const_iterator value_iter = values.begin();
      for( ; value_iter != values.end(); ++value_iter) {
        double Ereco = *value_iter;
        h.Fill(Ereco);
      }
      cout << "Obs >> Read: " << h.Integral() 
           << " " << util::AsString(kSKsFCFV1RingMuLike) << " events." << endl;
      for(int ibin = 1; ibin <= kNErecoBins1Rmu; ibin++) {
        double nobs = h.GetBinContent(ibin);
        fErDistrArray1Rmu->GetArray()[ibin-1] = nobs;
      }
    }
    else {
      cerr << "Obs >> Node " << datanode1Rmu.str() << " does not exist!" << endl;
    }
#endif
#ifndef __NUMUONLY__ 
    ostringstream datanode1Re;
    //Ereco spectrum
    datanode1Re << "superk/observation";
    datanode1Re << "/data_set_" << data_set_id;
    datanode1Re << "/one_e_ring/ereco_values";
    
    if(util::NodeExists(xml,root_node,datanode1Re.str())) {
      TH1D h("", "", kNErecoBins1Re,kErecoBinEdges1Re);
      h.SetDirectory(0);
      vector<double> values = util::XML2DblArray(xml,root_node,datanode1Re.str());
      vector<double>::const_iterator value_iter = values.begin();
      for( ; value_iter != values.end(); ++value_iter) {
        double Ereco = *value_iter;
        h.Fill(Ereco);
      }
      cout << "Obs >> Read: " << h.Integral() 
           << " " << util::AsString(kSKsFCFV1RingELike) << " events." << endl;
      for(int ibin = 1; ibin <= kNErecoBins1Re; ibin++) {
        double nobs = h.GetBinContent(ibin);
        fErDistrArray1Re->GetArray()[ibin-1] = nobs;
      }
    }

    else {
      cerr << "Obs >> Node " << datanode1Re.str() << " does not exist!" << endl;
    }
#endif
#endif

#ifdef __RHC_RUNNING__
#ifndef __NUEONLY__
    ostringstream datanodeRHC1Rmu;
    datanodeRHC1Rmu << "superk/observation";
    datanodeRHC1Rmu << "/data_set_" << data_set_id;
    datanodeRHC1Rmu << "/one_mu_ring_rhc/ereco_values"; 

    if(util::NodeExists(xml,root_node,datanodeRHC1Rmu.str())) {
      int idx = arr_index->SuperkSampleIdx(kSKsRHCFCFV1RingMuLike);
      TH1D h("", "", kRHCNErecoBins1Rmu,kRHCErecoBinEdges1Rmu);
      h.SetDirectory(0);
      vector<double> values = util::XML2DblArray(xml,root_node,datanodeRHC1Rmu.str());
      vector<double>::const_iterator value_iter = values.begin();
      for( ; value_iter != values.end(); ++value_iter) {
        double Ereco = *value_iter;
        h.Fill(Ereco);
      }
      cout << "Obs >> Read: " << h.Integral()
           << " " << util::AsString(kSKsRHCFCFV1RingMuLike) << " events." << endl;
      for(int ibin = 1; ibin <= kRHCNErecoBins1Rmu; ibin++) {
        double nobs = h.GetBinContent(ibin);
        fRHCErDistrArray1Rmu->GetArray()[ibin-1] = nobs;
      }
    }
    else {
      cerr << "Obs >> Node " << datanodeRHC1Rmu.str() << " does not exist!" << endl;
    }
#endif
#ifndef __NUMUONLY__
    ostringstream datanodeRHC1Re;
    //Ereco spectrum                                                                            
    datanodeRHC1Re << "superk/observation";
    datanodeRHC1Re << "/data_set_" << data_set_id;
    datanodeRHC1Re << "/one_e_ring_rhc/ereco_values"; 

    if(util::NodeExists(xml,root_node,datanodeRHC1Re.str())) {
      int idx = arr_index->SuperkSampleIdx(kSKsRHCFCFV1RingELike);
      TH1D h("", "", kRHCNErecoBins1Re,kRHCErecoBinEdges1Re);
      h.SetDirectory(0);
      vector<double> values = util::XML2DblArray(xml,root_node,datanodeRHC1Re.str());
      vector<double>::const_iterator value_iter = values.begin();
      for( ; value_iter != values.end(); ++value_iter) {
        double Ereco = *value_iter;
        h.Fill(Ereco);
      }
      cout << "Obs >> Read: " << h.Integral()
           << " " << util::AsString(kSKsRHCFCFV1RingELike) << " events." << endl;
      for(int ibin = 1; ibin <= kRHCNErecoBins1Re; ibin++) {
        double nobs = h.GetBinContent(ibin);
        fRHCErDistrArray1Re->GetArray()[ibin-1] = nobs;
      }
    }

    else {
      cerr << "Obs >> Node " << datanodeRHC1Re.str() << " does not exist!" << endl;
    }
#endif

#endif

    // Release memory
    xml->FreeDoc(xmldoc);
    delete xml;
  }//is_xml

  else { // DEPRECATE?
#ifndef __NUEONLY__ 
#ifdef __FHC_RUNNING__
    TFile f(config_file.c_str());
    TH1D * hist = 0;
    f.GetObject(Form("mock_data_%d",data_set_id),hist);
    if(!hist) {
      cerr << "Obs >> Histogram for 1Rmu does not exist!" << endl;
    }
    else {
      assert(hist->GetNbinsX() == util::NumErecoBins(kSKsFCFV1RingMuLike));
      double integral = 0;
      for(int ix = 1; ix <= hist->GetNbinsX(); ix++) {
        double nev   = hist->GetBinContent(ix);
        integral+=nev;
        fErDistrArray1Rmu->GetArray()[ix-1] = nev;
      }
      cout << "Obs >> Read: " << integral
           << " " << util::AsString(kSKsFCFV1RingMuLike) << " events." << endl;
    }
#endif
#ifdef __RHC_RUNNING__
    TFile fRHC(config_file.c_str());
    TH1D * histRHC = 0;
    fRHC.GetObject(Form("mock_data_%d",data_set_id),histRHC);
    if(!histRHC) {
      cerr << "Obs >> Histogram for RHC 1Rmu does not exist!" << endl;
    }
    else {
      assert(histRHC->GetNbinsX() == util::NumErecoBins(kSKsRHCFCFV1RingMuLike));
      double integral = 0;
      int idx = arr_index->SuperkSampleIdx(kSKsRHCFCFV1RingMuLike);
      for(int ix = 1; ix <= histRHC->GetNbinsX(); ix++) {
        double Ereco = histRHC->GetBinCenter(ix);
        double nev   = histRHC->GetBinContent(ix);
        integral+=nev;
        fRHCErDistrArray1Rmu->GetArray()[ix-1] = nev;
      }
      cout << "Obs >> Read: " << integral
           << " " << util::AsString(kSKsRHCFCFV1RingMuLike) << " events." << endl;
    }
#endif
#endif
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
  if      (s == kSKsRHCFCFV1RingMuLike) return fRHCErDistrArray1Rmu;
  else if (s == kSKsRHCFCFV1RingELike)  return fRHCErDistrArray1Re;
  else if (s == kSKsFCFV1RingELike)     return fErDistrArray1Re; 
  else if (s == kSKsFCFV1RingMuLike)    return fErDistrArray1Rmu; 
  std::cerr<<"SuperkObservation/ErecoSpectrumValorArray >> Failed for sample "<<util::AsString(s)<<std::endl;
  exit(-1);
}
//________________________________________________________________________________________
const VALORArray * SuperkObservation::ErecoSpectrumVArray (valor::conv::SuperkSample_t s) const
{
  if      (s == kSKsRHCFCFV1RingMuLike) return fRHCErDistrArray1Rmu;
  else if (s == kSKsRHCFCFV1RingELike)  return fRHCErDistrArray1Re;
  else if (s == kSKsFCFV1RingELike)     return fErDistrArray1Re; 
  else if (s == kSKsFCFV1RingMuLike)    return fErDistrArray1Rmu; 
  std::cerr<<"SuperkObservation/ErecoSpectrumValorArray >> Failed for sample "<<util::AsString(s)<<std::endl;
  exit(-1);
}
//________________________________________________________________________________________
double* SuperkObservation::ErecoSpectrumArray (valor::conv::SuperkSample_t s)
{
  return ErecoSpectrumVArray(s)->GetArray();
}
//________________________________________________________________________________________
const double* SuperkObservation::ErecoSpectrumArray (valor::conv::SuperkSample_t s) const
{
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
  cout << "Obs >> Generating a toy-MC expt with the following params: " << fp;
  cout << "Obs >> Stat. fluctuations are turned-" << ((fluctuate)? "ON" : "OFF") << endl;
  
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
  
  cout << "Obs >> Generating a toy-MC expt with the following params: " << ph;
  cout << "Obs >> Stat. fluctuations are turned-" << ((fluctuate)? "ON" : "OFF") << endl;
  
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

  util::SetDefaultStyle();

  for(int is=0; is<kns; is++) {
    SuperkSample_t s = kSKSamples[is];
    if(s == kSKsUndefined) continue;
    
    // This is a new object, leaks memory if data_collection has not been set as owner
    TH1D * data = this->ErecoSpectrum(s);
    
    data -> SetName( 
     Form( "%s_%s", 
            data_tag,
            util::AsStringCompact(s).c_str() ) 
    );

    data_collection -> Add(data);
  }
  
  return true;
}
//________________________________________________________________________________________
bool SuperkObservation::Write(
  const char * filename, const char * data_tag, Option_t * write_opt, Option_t * file_opt)
{
  if (!file_opt) file_opt="update";
  TObjArray * array = new TObjArray();
  array->SetOwner(true);
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
  
  //~ std::cerr<<"Statistic: "<<total_chisq<<"\t"<<phys_penalty<<endl;
  total_chisq += phys_penalty;
  return total_chisq;
}
//________________________________________________________________________________________
double SuperkObservation::Statistic (FitModelABC * model, const PhysHypothesis & ph, Option_t * stat_opt)
{
  bool applied_negevents_condition = model->Generate(ph);
  double total_chisq = this->Statistic(model,stat_opt);
  
  return total_chisq;
}
//________________________________________________________________________________________
double SuperkObservation::Statistic (FitModelABC * model, Option_t * stat_opt)
{
  SuperkExpectation * skexp = dynamic_cast<SuperkExpectation*>(model);
  assert(skexp);
  
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
      const int nbins = 5;
      double binedges[nbins+1] = {0.0,0.4,0.7,1.0,2.0,30.0}; 
    
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
}
//________________________________________________________________________________________
