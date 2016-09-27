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

//#define __show_read_events_mesg_debug__
//#define __show_num_of_events_mesg_debug__
//#define __show_num_of_bins_mesg_debug__
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
#include "valor/Utilities/Rnd.h"

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
#ifndef __USE_2D_TEMPLATES__ //1D templates
      fErDistrArray[is]=0;
#else
      f2DDistrArray[is]=0;
#endif
    }
    else {
#ifndef __USE_2D_TEMPLATES__ //1D templates
      fErDistrArray[is] = new VALORArray(util::NumErecoBins(s));
      fErDistrArray[is]->SetTitle(Form("Ereco distribution: %s",util::AsString(s).c_str()));
      fErDistrArray[is]->SetName(Form("%s",util::AsStringCompact(s).c_str()));
      fErDistrArray[is]->SetDimensionName(0, util::GetAxisNameCompact(0));
#else //2D templates
      f2DDistrArray[is]= new VALORArray(util::NumErecoBins(s),util::NumQ2Bins(s));   
      f2DDistrArray[is]->SetDimensionName(0, util::GetAxisNameCompact(0));
      f2DDistrArray[is]->SetDimensionName(1, util::GetAxisNameCompact(1));
      f2DDistrArray[is]->SetTitle(Form("2D Observed Spectrum Array %s", util::AsStringCompact(s).c_str()));
      f2DDistrArray[is]->SetName (Form("obs_2D_arr_%s", util::AsStringCompact(s).c_str()));
#endif //__USE_2D_TEMPLATES__
    }
  }
  
  arr_index = util::ArrIdx::Instance();
}
//________________________________________________________________________________________
void SuperkObservation::CleanUp(void)
{
  for (int is = 0 ; is < conv::kns ; is++) {
#ifndef __USE_2D_TEMPLATES__ //1D templates
    if (fErDistrArray[is]) delete fErDistrArray[is];
#else
    if (f2DDistrArray[is]) delete f2DDistrArray[is];
#endif
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
      exit(1);
    }
    XMLNodePointer_t root_node = xml->DocGetRootElement(xmldoc);
    
    for (int is = 0 ; is < conv::kns ; is++) {
      SuperkSample_t s = conv::kSKSamples[is];
      if (s==kSKsUndefined) continue;
#ifdef __show_read_events_mesg_debug__
      cout << "Obs >> Read: Reading events for sample: " << util::AsString(s) << endl;
#endif
      ostringstream datanode;
      datanode << "superk/observation";
      datanode << "/data_set_" << data_set_id;
      datanode << "/"<<util::GetXMLDatasetName(s);

#ifndef __USE_2D_TEMPLATES__ //1D templates
      datanode << "/"<<util::GetXMLDatasetParamsName();

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
#else //2D templates
      //first get the number of events
      int nevents = 0;
      {
        ostringstream datanodetmp;
        datanodetmp << datanode.str();
        datanodetmp << "/nevents";
        if(util::NodeExists(xml,root_node,datanodetmp.str())) {
          nevents = util::XML2Int(xml,root_node,datanodetmp.str());
        }
        else {
          cerr << "Obs >> Node " << datanodetmp.str() << " does not exist!" << endl;
        }
      }
      //now get each event
      datanode << "/"<<util::GetXMLDatasetParamsName();
      TH2D h("", "", util::NumErecoBins(s), util::ErecoBinEdges(s), util::NumQ2Bins(s), util::Q2BinEdges(s));
      h.SetDirectory(0);
      for(int iev = 0; iev < nevents; iev++) {
        ostringstream datanodetmp;
        datanodetmp << datanode.str();
        datanodetmp << "/ev" << iev;
        if(util::NodeExists(xml,root_node,datanodetmp.str())) {
          pair<double,double> values = util::XML2DblDblPair(xml,root_node,datanodetmp.str());
          double Ereco = values.first;
          double Q2    = values.second;
#ifdef __show_read_events_mesg_debug__
          cout << "Obs >> Read: event " << iev
               << " " << f2DDistrArray[is]->GetDimensionName(0) << " = " << Ereco
               << " " << f2DDistrArray[is]->GetDimensionName(1) << " = " << Q2
               << endl;
#endif
          h.Fill(Ereco, Q2);
        }
        else {
          cerr << "Obs >> Node " << datanodetmp.str() << " does not exist!" << endl;
        }
      }//iev
      cout << "Obs >> Read: " << h.Integral()
           << " " << util::AsString(s) << " events." << endl;
      for(int ibinx = 1; ibinx <= util::NumErecoBins(s); ibinx++) {
        for(int ibiny = 1; ibiny <= util::NumQ2Bins(s); ibiny++) {
          double nobs = h.GetBinContent(ibinx, ibiny);
          f2DDistrArray[is]->SetEntry(nobs, ibinx-1, ibiny-1);
        }//ibinx
      }//ibiny
#endif //__USE_2D_TEMPLATES__

    } // is
    
    // Release memory
    xml->FreeDoc(xmldoc);
    delete xml;
  }//is_xml

  else {
   
    TFile f(config_file.c_str());
      
#ifndef __USE_2D_TEMPLATES__ //1D templates
    TH1D * hist = 0;
#else
    TH2D * hist = 0;
#endif

    if(!f.IsOpen()){
      cerr << endl;
      cerr << "File " << config_file.c_str() << " is not defined" << endl;
      cerr << endl;
      exit(1);
    }

    for (int is = 0 ; is < conv::kns ; is++) {
      SuperkSample_t s = conv::kSKSamples[is];

      //select which histograms to read in
#if !defined(__FHC_RUNNING__)
      //skip the FHC samples
      if(s == kSKsFCFV1RingMuLike || s == kSKsFCFV1RingELike) continue;
#endif
#if !defined(__RHC_RUNNING__)
      //skip the RHC samples
      if(s == kSKsRHCFCFV1RingMuLike || s == kSKsRHCFCFV1RingELike) continue;
#endif
#if defined(__NUEONLY__)
      //skip the mu-like samples
      if(s == kSKsFCFV1RingMuLike || s == kSKsRHCFCFV1RingMuLike) continue;
#endif
#if defined(__NUMUONLY__)
      //skip the e-like samples
      if(s == kSKsFCFV1RingELike || s == kSKsRHCFCFV1RingELike) continue;
#endif
      
      if (s == conv::kSKsUndefined) continue;

#ifndef __USE_2D_TEMPLATES__ //1D templates
      else if(s == conv::kSKsRHCFCFV1RingMuLike) hist = (TH1D*) f.Get(Form("mock_data_RHC1Rmu_%d",data_set_id));
      else if(s == conv::kSKsRHCFCFV1RingELike)  hist = (TH1D*) f.Get(Form("mock_data_RHC1Re_%d",data_set_id));
      else if(s == conv::kSKsFCFV1RingELike)     hist = (TH1D*) f.Get(Form("mock_data_1Re_%d",data_set_id));
      else if(s == conv::kSKsFCFV1RingMuLike)    hist = (TH1D*) f.Get(Form("mock_data_1Rmu_%d",data_set_id));
      else if(s == conv::kSKsFCFVMultiRingELike) hist = (TH1D*) f.Get(Form("mock_data_MultiRe_%d",data_set_id));
#else //2D templates
      else if(s == conv::kSKsRHCFCFV1RingMuLike) hist = (TH2D*) f.Get(Form("mock_data_RHC1Rmu_%d",data_set_id));
      else if(s == conv::kSKsRHCFCFV1RingELike)  hist = (TH2D*) f.Get(Form("mock_data_RHC1Re_%d",data_set_id));
      else if(s == conv::kSKsFCFV1RingELike)     hist = (TH2D*) f.Get(Form("mock_data_1Re_%d",data_set_id));
      else if(s == conv::kSKsFCFV1RingMuLike)    hist = (TH2D*) f.Get(Form("mock_data_1Rmu_%d",data_set_id));
      else if(s == conv::kSKsFCFVMultiRingELike) hist = (TH2D*) f.Get(Form("mock_data_MultiRe_%d",data_set_id));
#endif
      else {
        cerr << endl;
        cerr << "Used SK sample " << util::AsString(s) << " not available" << endl;
        cerr << endl;
        exit(1);
      }
      
      if(!hist) {
        cerr << "Obs >> Histogram for " << util::AsString(s) << " does not exist!" << endl;
        if(s == conv::kSKsFCFV1RingMuLike)
          cerr << "Obs >> Need histogram name " << Form("mock_data_1Rmu_%d",data_set_id) << endl;
        else if(s == conv::kSKsRHCFCFV1RingMuLike)
          cerr << "Obs >> Need histogram name " << Form("mock_data_RHC1Rmu_%d",data_set_id) << endl;
        else if(s == conv::kSKsFCFV1RingELike)
          cerr << "Obs >> Need histogram name " << Form("mock_data_1Re_%d",data_set_id) << endl;
        else if(s == conv::kSKsRHCFCFV1RingELike)
          cerr << "Obs >> Need histogram name " << Form("mock_data_RHC1Re_%d",data_set_id) << endl;
        else{
          cerr << "Obs >> Histogram sample not defined" << endl;
          exit(1);
        }
        exit(1);
      }


#ifndef __USE_2D_TEMPLATES__ //1D templates
      for(int ibin = 1; ibin <= util::NumErecoBins(s); ibin++) {
        double nobs = hist->GetBinContent(ibin);
        fErDistrArray[is]->GetArray()[ibin-1] = nobs;
      }
#else //2D templates
      
      for(int ibinx = 1; ibinx <= util::NumErecoBins(s); ibinx++) {
        for(int ibiny = 1; ibiny <= util::NumQ2Bins(s); ibiny++) {
          // histogram convention (Q2,Ereco) is opposite of VALORArray's one (Ereco,Q2)
          double nobs = hist->GetBinContent(ibiny, ibinx);
          f2DDistrArray[is]->SetEntry(nobs, ibinx-1, ibiny-1);
        }//ibiny
      }//ibinx
      
#endif //__USE_2D_TEMPLATES__
      
    }//is
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

#ifndef __USE_2D_TEMPLATES__ //1D templates
  nev = ErecoSpectrumVArray(s)->Integral();
#else
  nev = TwoDSpectrumVArray(s)->Integral();
#endif

#ifdef __show_num_of_events_mesg_debug__
  cout << "Obs >> NEvents(" 
       << util::AsString(s) << ") = " << nev << endl;
#endif

  return nev;
}
//________________________________________________________________________________________
#ifndef __USE_2D_TEMPLATES__ //1D templates 
double SuperkObservation::NumOfEvents(double Er, SuperkSample_t s) const
#else
double SuperkObservation::NumOfEvents(double xvalue, double yvalue, SuperkSample_t s) const
#endif
{
  double nev = 0;

  if(s == kSKsUndefined) {
    return 0;
  }
  
#ifndef __USE_2D_TEMPLATES__ //1D templates 
  cerr<<"WARNING: SuperkObservation::NumOfEvents(Ereco) needs testing"<<std::endl;
  // NEEDS TESTING

  const double * spectrum = this->ErecoSpectrumArray(s);
  
  int bin = util::ErecoBinContains(s,Er);
  assert(bin!=-1);
  nev = spectrum[util::ErecoBinContains(s,Er)];
#else //2D templates
  const VALORArray * spectrum = this->TwoDSpectrumVArray(s);

  int binx = -1, biny = -1;
  if(xvalue >= 0)
    binx = util::ErecoBinContains(s, xvalue);
  if(yvalue >= 0)
    biny = util::Q2BinContains(s, yvalue);

  nev = spectrum->Integral(binx,binx, biny,biny);
#endif //__USE_2D_TEMPLATES__

#ifdef __show_num_of_events_mesg_debug__
  cout << "Obs >> NEvents(" 
       << util::AsString(s) << ", ";
#ifndef __USE_2D_TEMPLATES__
  cout << "Ereco = " << Er << " GeV";
#else
  cout << spectrum->GetDimensionName(0) << " = ";
  if(xvalue < 0) cout << "integrated";
  else cout << xvalue;
  cout << ", "
       << spectrum->GetDimensionName(1) << " = ";
  if(yvalue < 0) cout << "integrated";
  else cout << yvalue;
#endif
  cout << ") = " << nev << endl;
#endif

  return nev;
}
//________________________________________________________________________________________
int SuperkObservation::NumOfBins(Option_t * /*opt*/) const
{
  int nbins = 0;
#ifdef __show_num_of_bins_mesg_debug__
  cout << "SuperkObservation::NumOfBins() ";
#endif
  for (int is = 0 ; is < kns ; is++) {
    SuperkSample_t s = kSKSamples[is];
    if (s==kSKsUndefined) continue;
#ifndef __USE_2D_TEMPLATES__ //1D templates
    nbins+=util::NumErecoBins(s);
#else
    nbins += util::NumErecoBins(s) * util::NumQ2Bins(s);
#ifdef __show_num_of_bins_mesg_debug__
    cout << util::AsStringCompact(s) << ": " 
         << util::NumErecoBins(s) * util::NumQ2Bins(s) << " + ";
#endif
#endif //__USE_2D_TEMPLATES__
  }
#ifdef __show_num_of_bins_mesg_debug__
  cout << " = " << nbins << std::endl;
#endif
  return nbins;
}
//________________________________________________________________________________________


//________________________________________________________________________________________
#ifndef __USE_2D_TEMPLATES__ //1D templates
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
//#endif //__USE_2D_TEMPLATES__
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
#else
//#ifdef __USE_2D_TEMPLATES__ //2D templates
//________________________________________________________________________________________
TH2D* SuperkObservation::TwoDSpectrum (valor::conv::SuperkSample_t s) const
{
  if(s == kSKsUndefined) {
    return 0;
  }

  const int nbins_q2    = util::NumQ2Bins(s);
  const int nbins_ereco = util::NumErecoBins(s);
  typedef double (erq2type)  [nbins_ereco][nbins_q2];  

  const erq2type & arr = *reinterpret_cast< const erq2type * > (TwoDSpectrumArray(s));
  if (!arr) return 0;

  string name (Form("obs_2D_%s"               ,util::AsStringCompact(s).c_str() ) );
  string title(Form("2D Observed Spectrum: %s",util::AsString(s).c_str() ) );

  TH2D * data = new TH2D(name.c_str(), title.c_str(), nbins_q2, util::Q2BinEdges(s), nbins_ereco, util::ErecoBinEdges(s) );
  util::FormatErecoVsQ2Hist(data);
  data -> SetDirectory(0);
  data -> SetMarkerStyle(20);
  data -> SetMarkerSize(1);
  
  for (int ir = 0 ; ir < nbins_ereco ; ir++ ) {
    for (int iq = 0 ; iq < nbins_q2 ; iq++ ) {
      data->SetBinContent(iq+1,ir+1,arr[ir][iq]);
    }//iq
  }//ir
  
  return data;
}
//________________________________________________________________________________________
TH1D* SuperkObservation::XSpectrum (valor::conv::SuperkSample_t s) const
{
  if(s == kSKsUndefined) {
    return 0;
  }

  const int nbins = util::NumErecoBins(s);
  
  const VALORArray * varray = TwoDSpectrumVArray(s);
  if (!varray) return 0;

  string name (Form("obs_%s_%s", varray->GetDimensionName(0).c_str(), util::AsStringCompact(s).c_str() ) );
  string title(Form("1D Observed Spectrum: %s", util::AsString(s).c_str() ) );

  TH1D * data = new TH1D(name.c_str(), title.c_str(), nbins, util::ErecoBinEdges(s) );
  util::FormatErecoHist(data);
  data -> SetDirectory(0);
  data -> SetMarkerStyle(20);
  data -> SetMarkerSize(1);
  
  for (int ib = 0 ; ib < nbins ; ib++ ) {
    data->SetBinContent(ib+1, varray->Integral(ib,ib, -1,-1));
  }//ib
  
  return data;
}
//________________________________________________________________________________________
TH1D* SuperkObservation::YSpectrum (valor::conv::SuperkSample_t s) const
{
  if(s == kSKsUndefined) {
    return 0;
  }

  const int nbins = util::NumQ2Bins(s);

  const VALORArray * varray = TwoDSpectrumVArray(s);
  if (!varray) return 0;

  string name (Form("obs_%s_%s", varray->GetDimensionName(1).c_str(), util::AsStringCompact(s).c_str() ) );
  string title(Form("1D Observed Spectrum: %s", util::AsString(s).c_str() ) );

  TH1D * data = new TH1D(name.c_str(), title.c_str(), nbins, util::Q2BinEdges(s) );
  util::FormatQ2Hist(data);
  data -> SetDirectory(0);
  data -> SetMarkerStyle(20);
  data -> SetMarkerSize(1);
  
  for (int ib = 0 ; ib < nbins ; ib++ ) {
    data->SetBinContent(ib+1, varray->Integral(-1,-1, ib,ib));
  }//ib
  
  return data;
}
//________________________________________________________________________________________
VALORArray * SuperkObservation::TwoDSpectrumVArray (valor::conv::SuperkSample_t s)
{
  return f2DDistrArray[s];
}
//________________________________________________________________________________________
const VALORArray * SuperkObservation::TwoDSpectrumVArray (valor::conv::SuperkSample_t s) const
{
  return f2DDistrArray[s];
}
//________________________________________________________________________________________
double* SuperkObservation::TwoDSpectrumArray (valor::conv::SuperkSample_t s)
{
  VALORArray * va = TwoDSpectrumVArray(s);
  if (not va) return 0;
  return TwoDSpectrumVArray(s)->GetArray();
}
//________________________________________________________________________________________
const double* SuperkObservation::TwoDSpectrumArray (valor::conv::SuperkSample_t s) const
{
  const VALORArray * va = TwoDSpectrumVArray(s);
  if (not va) return 0;
  return TwoDSpectrumVArray(s)->GetArray();
}
#endif //__USE_2D_TEMPLATES__
//________________________________________________________________________________________

//________________________________________________________________________________________
bool SuperkObservation::Generate(
     FitModelABC * model, bool fluctuate)
{
  SuperkExpectation * skexp = dynamic_cast<SuperkExpectation*>(model);
  
  assert(skexp);
 
  for (int is = 0 ; is < conv::kns ; is++ ) {
    SuperkSample_t s = kSKSamples[is];
    if (s==kSKsUndefined) continue;
#ifndef __USE_2D_TEMPLATES__ //1D templates
    for(int ibin = 0; ibin < util::NumErecoBins(s); ibin++) {
       double nexp  = skexp->CurrPredictedErecoSpectrumArray(s)[ibin];
       double nobs  = (fluctuate) ? Rnd::Instance()->RndToyMCstat().Poisson(nexp) : nexp;
       this->ErecoSpectrumArray(s)[ibin] = nobs;
    }//ibin
#else //2D templates
    for(int ibinx = 0; ibinx < util::NumErecoBins(s); ibinx++) {
      for(int ibiny = 0; ibiny < util::NumQ2Bins(s); ibiny++) {
        double nexp = skexp->CurrPredictedErecoVsQ2VArray(s)->GetEntry(ibinx, ibiny);
        double nobs = (fluctuate) ? Rnd::Instance()->RndToyMCstat().Poisson(nexp) : nexp;
        this->TwoDSpectrumVArray(s)->SetEntry(nobs, ibinx, ibiny);
      }//ibiny
    }//ibinx
#endif //__USE_2D_TEMPLATES__
  }//is

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

#ifndef __USE_2D_TEMPLATES__ //1D templates 
      // This is a new object, leaks memory if data_collection has not been set as owner
      TH1D * data = this->ErecoSpectrum(s);      
#else //2D templates
      // This is a new object, leaks memory if data_collection has not been set as owner
      TH2D * data = this->TwoDSpectrum(s);

      TH1D * datax = this->XSpectrum(s);
      datax->SetName(Form("x_%s_%s", data_tag, util::AsStringCompact(s).c_str()));
      data_collection->Add(datax);
      TH1D * datay = this->YSpectrum(s);
      datay->SetName(Form("y_%s_%s", data_tag, util::AsStringCompact(s).c_str()));
      data_collection->Add(datay);
#endif //__USE_2D_TEMPLATES__

      data->SetName(Form("%s_%s", data_tag, util::AsStringCompact(s).c_str()));
      data_collection -> Add(data);
    }
    else {
#ifndef __USE_2D_TEMPLATES__ //1D templates  
      VALORArray * data = this->ErecoSpectrumVArray(s);
#else //2D templates
      VALORArray * data = this->TwoDSpectrumVArray(s);
#endif //__USE_2D_TEMPLATES__
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
double SuperkObservation::GetFitStatistic (FitModelABC * model, const FitParamList & fp,   Option_t * stat_opt)
{
  bool applied_negevents_condition = model->Generate(fp);
  double total_chisq = this->GetFitStatistic(model,stat_opt);
  
  if (applied_negevents_condition) cerr<<"SuperkObservation/Statistic >> Warning: applied negative events condition"<<std::endl;
  
  return total_chisq;
}
//________________________________________________________________________________________
double SuperkObservation::GetFitStatistic (FitModelABC * model, const PhysHypothesis & ph, Option_t * stat_opt)
{
  bool applied_negevents_condition = model->Generate(ph);
  double total_chisq = this->GetFitStatistic(model,stat_opt);
  if (applied_negevents_condition) cerr<<"SuperkObservation/Statistic >> Warning: applied negative events condition"<<std::endl;
  
  return total_chisq;
}
//________________________________________________________________________________________
double SuperkObservation::GetFitStatistic (FitModelABC * model, Option_t * stat_opt)
{
  // We will try casting to an SK model
  SuperkExpectation * skexp = dynamic_cast<SuperkExpectation*>(model);
  
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
      int n_ereco_bins = util::NumErecoBins(s);

#ifndef __USE_2D_TEMPLATES__ //1D templates
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
#ifndef __RHC_RUNNING__
        const int nbins = 5;
        double binedges[nbins+1] = {0.0,0.4,0.7,1.0,2.0,30.0};
#else
        const int nbins = 2;
        double binedges[nbins+1] = {0.0,2.5,30.0};
#endif
      
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
          cout << "SuperkObservation::Statistic()> delta_chisq (s = "<<util::AsString(s)<<", bin = " << ibin 
             << "; [" << xlo << ", " << xup << "] GeV) = " << delta_chisq << endl;
          cout << "\tNobs: " << nobs << "\tNexp:" << nexp << endl;
  #endif
        }//ibin
      } //goodness of fit

#else //2D templates

      const VALORArray * obs_v = this->TwoDSpectrumVArray(s);
      const VALORArray * exp_v = skexp->CurrPredictedErecoVsQ2VArray(s);
      int n_q2_bins = util::NumQ2Bins(s);
      
      if (not gof) {
          
        for (int ibinx = 0 ; ibinx < n_ereco_bins ; ibinx++) {
          for (int ibiny = 0 ; ibiny < n_q2_bins ; ibiny++) {
            nobs = obs_v->GetEntry(ibinx, ibiny);
            nexp = exp_v->GetEntry(ibinx, ibiny);
          
            if    (gauss) { delta_chisq = util::ChisqGauss  (nobs,nexp); }
            else          { delta_chisq = util::ChisqPoisson(nobs,nexp); }
          
            total_chisq += delta_chisq;
        
#ifdef __show_bin_contributions_to_chisq__
            double xlo = util::ErecoBinEdges(s,ibinx);
            double xup = util::ErecoBinEdges(s,ibinx+1);
            double ylo = util::ErecoBinEdges(s,ibiny);
            double yup = util::ErecoBinEdges(s,ibiny+1);
            cout << "SuperkObservation::Statistic >> delta_chisq (s = "<<util::AsString(s)
                 << ", bin = " << ibinx << "," << ibiny 
                 << "; " << obs_v->GetDimensionName(0) <<  " = [" << xlo << ", " << xup << "]"
                 << "; " << obs_v->GetDimensionName(1) <<  " = [" << ylo << ", " << yup << "]"
                 << ") = " << delta_chisq << endl;
            cout << "\tNobs: " << nobs << "\tNexp:" << nexp << endl;
#endif
          }//ibiny
        }//ibinx

      }// not goodness of fit
      else {
        const int nbinsx = 2;
        double binedgesx[nbinsx+1] = {0.0,2.5,30.0};
        const int nbinsy = 2;
        double binedgesy[nbinsy+1] = {0.0,1.0,30.0};

        cerr << "Warning** SuperKObservation goodness of fit arrays for 2D templates have arbitrary values set. Please update!" << endl;
        cerr << "SuperKObservation goodness of fit arrays should be moved to Binning.h" << endl;
      
        //loop over Ereco bins
        for(int ibinx = 0; ibinx < nbinsx; ibinx++) {
          int ErecoMinBin = TMath::BinarySearch(n_ereco_bins,util::ErecoBinEdges(s),binedgesx[ibinx]);
          int ErecoMaxBin = TMath::BinarySearch(n_ereco_bins,util::ErecoBinEdges(s),binedgesx[ibinx+1]);
          for(int ibiny = 0; ibiny < nbinsy; ibiny++) {
            int Q2MinBin = TMath::BinarySearch(n_q2_bins,util::Q2BinEdges(s),binedgesy[ibiny]);
            int Q2MaxBin = TMath::BinarySearch(n_q2_bins,util::Q2BinEdges(s),binedgesy[ibiny+1]);

            nobs = obs_v->Integral(ErecoMinBin, ErecoMaxBin, Q2MinBin, Q2MaxBin);
            nexp = exp_v->Integral(ErecoMinBin, ErecoMaxBin, Q2MinBin, Q2MaxBin);
          
            if    (gauss) { delta_chisq = util::ChisqGauss  (nobs,nexp); }
            else          { delta_chisq = util::ChisqPoisson(nobs,nexp); }
          
            total_chisq += delta_chisq;
            
#ifdef __show_bin_contributions_to_chisq__
            double xlo = util::ErecoBinEdges(s,ibinx);
            double xup = util::ErecoBinEdges(s,ibinx+1);
            double ylo = util::ErecoBinEdges(s,ibiny);
            double yup = util::ErecoBinEdges(s,ibiny+1);
            cout << "SuperkObservation::Statistic >> delta_chisq (s = "<<util::AsString(s)
                 << ", bin = " << ibinx << "," << ibiny 
                 << "; " << obs_v->GetDimensionName(0) <<  " = [" << xlo << ", " << xup << "]"
                 << "; " << obs_v->GetDimensionName(1) <<  " = [" << ylo << ", " << yup << "]"
                 << ") = " << delta_chisq << endl;
            cout << "\tNobs: " << nobs << "\tNexp:" << nexp << endl;
  #endif
          }//ibiny
        }//ibinx
      } //goodness of fit

#endif //__USE_2D_TEMPLATES__

    }//is
    
    return total_chisq;
  } // end of SKExp cast
  
  // If we get here, the cast has failed.
  std::cerr<<"SuperkOberservation/Generate >> Failed to interpret model: "<<model->Name()<<std::endl;
  exit(-1);
}
//________________________________________________________________________________________
