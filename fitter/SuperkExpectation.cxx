//________________________________________________________________________________________
/*
 
 \copyright (c) 2010-2013, VaLOR group.

30-12-2013  Nick
Added samples for reverse horn current running

21-7-2014 Davide
Added samples for multiring forward horn current

*/
//________________________________________________________________________________________

//#define __show_syst_reset_mesg__
//#define __show_syst_tweak_mesg__
//#define __show_syst_tweak_debug__

#include <cassert>
#include <cstdlib>
#include <iostream>

#include <TSystem.h>
#include <TFile.h>
#include <TObjArray.h>
#include <TTree.h>
#include <TString.h>
#include <TLegend.h>
#include <TSpline.h>
#include <TRandom3.h>

#include "SuperkExpectation.h"
#include "SuperkDetectorFSI.h"
#include "SuperkMC.h"
#include "NuFlux.h"
#include "NuXSec.h"
#include "BANFF.h"
#include "Utilities.h"
#include "Conventions/Opt.h"
#include "Utilities/Rnd.h"

using std::cout;
using std::cerr;
using std::endl;

using namespace valor;
using namespace valor::conv;

ClassImp(SuperkExpectation)

//________________________________________________________________________________________
SuperkExpectation::SuperkExpectation()
{
  this->Init();
}
//________________________________________________________________________________________
SuperkExpectation::~SuperkExpectation(void)
{
  this->CleanUp();
}
//________________________________________________________________________________________
void SuperkExpectation::Init(void)
{
  fIsConfigured = false;

  fSuperkMC          = 0;
  fSuperkDetectorFSI = 0;
  fNuFlux            = 0;
  fNuXSec            = 0;
  fBANFF             = 0;
  
  for (int is = 0 ; is < conv::kns ; is++) {
    SuperkSample_t s = conv::kSKSamples[is];
    if (s==kSKsUndefined) {
      fPredictErecoSpectrumArrayTot [is] = 0;
      fPredictErecoSpectrumArray    [is] = 0;
      // NEW 2D TEMPLATES
      fPredictErecoVsQ2SpectrumArrayTot [is] = 0;
      fPredictErecoVsQ2SpectrumArray    [is] = 0;
      //
      continue;
    }
    fPredictErecoSpectrumArrayTot [is]= new VALORArray(util::NumErecoBins(s));
    fPredictErecoSpectrumArray    [is]= new VALORArray(conv::knm, util::NumErecoBins(s));
    
    fPredictErecoSpectrumArrayTot[is] ->SetDimensionName(0, "Ereco");
    fPredictErecoSpectrumArray   [is] ->SetDimensionName(0, "SK Mode");
    fPredictErecoSpectrumArray   [is] ->SetDimensionName(1, "Ereco");
    
    fPredictErecoSpectrumArrayTot[is]->SetTitle(Form("Predicted Ereco Spectrum Array %s total"  ,util::AsStringCompact(s).c_str()));
    fPredictErecoSpectrumArray   [is]->SetTitle(Form("Predicted Ereco Spectrum Array %s by mode",util::AsStringCompact(s).c_str()));
    fPredictErecoSpectrumArrayTot[is]->SetName(Form("ereco_%s_total"  , util::AsStringCompact(s).c_str()));
    fPredictErecoSpectrumArray   [is]->SetName(Form("ereco_%s_mode"   , util::AsStringCompact(s).c_str()));

    // NEW 2D TEMPLATES
    fPredictErecoVsQ2SpectrumArrayTot [is]= new VALORArray(util::NumErecoBins(s),util::NumQ2Bins(s));
    fPredictErecoVsQ2SpectrumArray    [is]= new VALORArray(conv::knm, util::NumErecoBins(s), util::NumQ2Bins(s));
    
    fPredictErecoVsQ2SpectrumArrayTot[is] ->SetDimensionName(0, "Ereco");
    fPredictErecoVsQ2SpectrumArrayTot[is] ->SetDimensionName(1, "Q2");

    fPredictErecoVsQ2SpectrumArray   [is] ->SetDimensionName(0, "SK Mode");
    fPredictErecoVsQ2SpectrumArray   [is] ->SetDimensionName(1, "Ereco");
    fPredictErecoVsQ2SpectrumArray   [is] ->SetDimensionName(2, "Q2");
    
    fPredictErecoVsQ2SpectrumArrayTot[is]->SetTitle(Form("Predicted Ereco Vs Q2 Spectrum Array %s total"  ,util::AsStringCompact(s).c_str()));
    fPredictErecoVsQ2SpectrumArray   [is]->SetTitle(Form("Predicted Ereco Vs Q2 Spectrum Array %s by mode",util::AsStringCompact(s).c_str()));
    fPredictErecoVsQ2SpectrumArrayTot[is]->SetName(Form("ereco_Vs_q2_%s_total"  , util::AsStringCompact(s).c_str()));
    fPredictErecoVsQ2SpectrumArray   [is]->SetName(Form("ereco_Vs_q2_%s_mode"   , util::AsStringCompact(s).c_str()));
    //
  }
  
}
//________________________________________________________________________________________
void SuperkExpectation::Book(void)
{
  return;
}
//________________________________________________________________________________________
void SuperkExpectation::CleanUp(void)
{
  if (fSuperkDetectorFSI ) delete fSuperkDetectorFSI;
  if (fSuperkMC          ) delete fSuperkMC;
  if (fNuFlux            ) delete fNuFlux; 
  if (fNuXSec            ) delete fNuXSec;
  if (fBANFF             ) delete fBANFF;
}

//________________________________________________________________________________________
bool SuperkExpectation::Read(string config_file, Option_t * /*opt*/)
{  
  T2KOpt::Instance()->ReadXmlPath(config_file); 
  
  cout << "SKExp/Read >> Reading config file: " << config_file << endl;

  if(fIsConfigured) {
    this->CleanUp();
  } 

  // Instantiate all necessary histograms
  this->Book();

  // Create components used for building-up the p.d.f & implemeting its systematic variation
  fBANFF             = new BANFF;
  fNuFlux            = new NuFlux; 
  fNuXSec            = new NuXSec; 
  fSuperkDetectorFSI = new SuperkDetectorFSI;
  fSuperkMC          = new SuperkMC;

  // Configure inter-dependencies.
  fNuFlux   -> TellMeAboutTheBANFFfit   (fBANFF);
  fNuFlux   -> TellMeAboutTheRHCBANFFfit(fBANFF);
  fNuXSec   -> TellMeAboutTheBANFFfit   (fBANFF);
  fNuXSec   -> TellMeAboutTheRHCBANFFfit(fBANFF);
  fSuperkMC -> TellMeAboutCrossSections (fNuXSec);
  fSuperkMC -> TellMeAboutTheFlux       (fNuFlux);
  fSuperkMC -> TellMeAboutTheDetectorFSI(fSuperkDetectorFSI);


  // Parse the input XML configuration file
  TXMLEngine* xml = new TXMLEngine;   
  XMLDocPointer_t xmldoc = xml->ParseFile(config_file.c_str());
  if (xmldoc==0) {
    cerr << "SKExp/Read >> " << "Can not parse input: " << config_file << endl;
    delete xml;
    return false;  
  }
  XMLNodePointer_t root_node = xml->DocGetRootElement(xmldoc);

  // Configure each of the components 

  // NEW FLUX FHC+RHC
#if defined(__FHC_RUNNING__) && defined(__RHC_RUNNING__) && !defined(__FULLY_CORRELATED__)
  
  if ( ! this->GetBANFF            () -> ReadFHCpRHC (xml, root_node) ) return false;
  
  if ( ! this->GetNuFlux           () -> Read (xml, root_node) ) return false;
  if ( ! this->GetNuFlux           () -> ReadRHC (xml, root_node) ) return false;

  if ( ! this->GetNuXSec           () -> Read (xml, root_node) ) return false;

  // NEW FHC+RHC JOINT
  //if ( ! this->GetSuperkDetectorFSI() -> Read (xml, root_node) ) return false;
  //if ( ! this->GetSuperkDetectorFSI() -> ReadRHC (xml, root_node) ) return false;
  if ( ! this->GetSuperkDetectorFSI() -> ReadFHCpRHC (xml, root_node) ) return false;
  //
#else
  //

  if ( ! this->GetBANFF            () -> Read (xml, root_node) ) return false;
  if ( ! this->GetBANFF            () -> ReadRHC (xml, root_node) ) return false;
  if ( ! this->GetNuFlux           () -> Read (xml, root_node) ) return false;
  if ( ! this->GetNuFlux           () -> ReadRHC (xml, root_node) ) return false;
  if ( ! this->GetNuXSec           () -> Read (xml, root_node) ) return false;
  if ( ! this->GetSuperkDetectorFSI() -> Read (xml, root_node) ) return false;
  if ( ! this->GetSuperkDetectorFSI() -> ReadRHC (xml, root_node) ) return false;
#endif

  //

  /// NEW MULTIRING 2
#ifdef __MULTIRING_FHC_RUNNING__
  if ( ! this->GetSuperkDetectorFSI() -> ReadMultiRing (xml, root_node) ) return false; /// Take same cov. matrix of FHC sample
#endif
  ///

  if ( ! this->GetSuperkMC         () -> Read (xml, root_node) ) return false;

  fIsConfigured = true;

  // Release memory
  xml->FreeDoc(xmldoc);
  delete xml;

  cout << "Exp >> Done configuring Superk spectrum predictor... \n\n" << endl;
  
  return true;
}
//________________________________________________________________________________________
bool SuperkExpectation::Generate(const PhysHypothesis & ph, Option_t * /*opt*/)
{

// Calculate dN^{exp}_{SuperK}(E_{reco}) / dE_{reco} for the input E_{reco} and 
// the input physics hypothesis
//
  if(!fIsConfigured) {
    cerr << "Exp >> Not configured properly! " << endl;
    exit(1);
  }

  // set params driving the oscillation, decoherence or decay model
  fSuperkMC->Hypothesis().Copy(ph);

  // update SK MC templates
  fSuperkMC->UpdateTemplates(); // it already supports both 1D and 2D templates

#ifndef __USE_2D_TEMPLATES__ // 1D templates

  // update Ereco spectrum prediction
  bool applied_negevents_condition = fSuperkMC->UpdateErecoSpectra();

  // copy histograms
  for(int is = 0; is < kns; is++) {
    SuperkSample_t s = kSKSamples[is];
    if(s == kSKsUndefined) continue;
    const int nbins = kNErecoBins[is];
    for(int ibin = 1; ibin <= nbins; ibin++) {
      double ns = fSuperkMC->CurrErecoSpectrumArray(s)[ibin-1];
      this->CurrPredictedErecoSpectrumArray(s)[ibin-1] = ns;
      for(int im = 0; im < knm; im++) {
        SuperkMode_t m = kSKReactionModes[is][im];
        if(m == kSKmUndefined) continue;
        double nsm = fSuperkMC->CurrErecoSpectrumArray(s,m)[ibin-1];
        this->CurrPredictedErecoSpectrumArrayMode(s,m)[ibin-1] = nsm;
      }//m
    }//ibin
  }//s

#else // 2D templates

  // NEW 2D TEMPLATES
  
  // update Ereco spectrum prediction
  bool applied_negevents_condition = fSuperkMC->UpdateErecoVsQ2Spectra();

  // copy histograms
  for(int is = 0; is < kns; is++) {
    SuperkSample_t s = kSKSamples[is];
    if(s == kSKsUndefined) continue;

    const int nbins_q2    = kNQ2Bins[is];
    const int nbins_ereco = kNErecoBins[is];

    typedef double (erq2type)  [nbins_ereco][nbins_q2];
    const erq2type & arr_skmc     = *reinterpret_cast< const erq2type * > ( fSuperkMC->CurrErecoVsQ2SpectrumArray(s) );
    if (!arr_skmc) return 0;

    erq2type & arr_skexp    = *reinterpret_cast< erq2type * > ( this->CurrPredictedErecoVsQ2SpectrumArray(s) );
    if (!arr_skexp) return 0;
    
    for(int ibin_ereco = 1; ibin_ereco <= nbins_ereco; ibin_ereco++) {
      for(int ibin_q2 = 1; ibin_q2 <= nbins_q2; ibin_q2++) {
	
	//double ns = fSuperkMC->CurrErecoVsQ2SpectrumArray(s)[ibin_ereco-1][ibin_q2-1];
	//this->CurrPredictedErecoVsQ2SpectrumArray(s)[ibin_ereco-1][ibin_q2-1] = ns;

	double ns = arr_skmc[ibin_ereco-1][ibin_q2-1];
	arr_skexp[ibin_ereco-1][ibin_q2-1] = ns;

	// now per mode
 
	for(int im = 0; im < knm; im++) {
	  SuperkMode_t m = kSKReactionModes[is][im];
	  if(m == kSKmUndefined) continue;

	  const erq2type & arr_skmc_mode = *reinterpret_cast< const erq2type * > (fSuperkMC->CurrErecoVsQ2SpectrumArray(s,m));  
	  if (!arr_skmc_mode) return 0;
	  
	  erq2type & arr_skexp_mode = *reinterpret_cast< erq2type * > (this->CurrPredictedErecoVsQ2SpectrumArrayMode(s,m));  
	  if (!arr_skexp_mode) return 0;
	  
	  //double nsm = fSuperkMC->CurrErecoVsQ2SpectrumArray(s,m)[ibin_ereco-1][ibin_q2-1]; 
	  //this->CurrPredictedErecoVsQ2SpectrumArrayMode(s,m)[ibin_ereco-1][ibin_q2-1] = nsm; 
	
	  double nsm = arr_skmc_mode[ibin_ereco-1][ibin_q2-1]; 
	  arr_skexp_mode[ibin_ereco-1][ibin_q2-1] = nsm; 

	}//m
      }//ibin q2
    }//ibin ereco
  }//s
  
  //
  
  
  
  
#endif

  return applied_negevents_condition;
}
//________________________________________________________________________________________
bool SuperkExpectation::Generate(const FitParamList & fp, Option_t * /*opt*/)
{
  return Generate(fp, false, false);
}
//________________________________________________________________________________________
bool SuperkExpectation::Generate(const FitParamList & fp, bool all_banff_fit, bool all_1Re_eff)
{

// SuperkExpectation interface method for fitters.
// Calculate dN^{exp}_{SuperK}(E_{reco}) / dE_{reco} for the input E_{reco} and 
// the input physics + nuisance parameter list.
//

  // set values of nuisance parameters 
  // (if any nuisance parameter is included in the fit)
  this->SetNuisanceParams(fp, all_banff_fit, all_1Re_eff);

  // extract the values of params driving the oscillation, decoherence or decay model
  PhysHypothesis ph(fp);

  // call main Generate() method to generate SuperK predictions for the input
  // physics hypothesis using the (potentially) tweaked p.d.fs
  bool applied_negevents_condition = this->Generate(ph);
  return applied_negevents_condition;
}
//________________________________________________________________________________________
TH1D* SuperkExpectation::CurrPredictedErecoSpectrum(valor::conv::SuperkSample_t s) const 
{
  const double * arr = CurrPredictedErecoSpectrumArray(s);
  if (!arr) return 0;
  
  string name (Form("curr_ereco_%s"                       ,util::AsStringCompact(s).c_str() ) );
  string title(Form("Current Predicted Ereco Spectrum: %s",util::AsString(s).c_str() ) );
  
  int nbins = util::NumErecoBins(s);
  TH1D * h = new TH1D(name.c_str(), title.c_str(), nbins, util::ErecoBinEdges(s) );
  util::FormatErecoHist(h);
  for (int i = 0 ; i < nbins ; i++ ) {
    h->SetBinContent(i+1,arr[i]);
  }
  return h;
}
//________________________________________________________________________________________
TH1D* SuperkExpectation::CurrPredictedErecoSpectrum(
  valor::conv::SuperkSample_t s, valor::conv::SuperkMode_t m) const 
{ 
  const double * arr = CurrPredictedErecoSpectrumArrayMode(s,m);
  if (!arr) return 0;
  
  int nbins = util::NumErecoBins(s);
  string name  (Form("curr_ereco_%s_%s"                       ,util::AsStringCompact(s).c_str(), util::AsStringCompact(m).c_str()) );
  string title (Form("Current Predicted Ereco Spectrum: %s %s",util::AsString(s).c_str()       , util::AsString(m).c_str()       ) );
  TH1D * h = new TH1D(name.c_str(), title.c_str(), nbins, util::ErecoBinEdges(s) );
  util::FormatErecoHist(h);
  util::SetDefaultSKSpectrumStyle(s,m,h);
  for (int i = 0 ; i < nbins ; i++ ) {
    h->SetBinContent(i+1,arr[i]);
  }
  return h;
}
//________________________________________________________________________________________
VALORArray * SuperkExpectation::CurrPredictedErecoVArray(valor::conv::SuperkSample_t s) 
{
  VALORArray * arr = fPredictErecoSpectrumArrayTot[s];
  if (!arr) {
    std::cerr<<"SKExp/CurrPredictedErecoVArray could not find array for sample "<<util::AsString(s)<<std::endl;
    exit(-1);
  }
  
  return arr;
}
//________________________________________________________________________________________
VALORArray * SuperkExpectation::CurrPredictedErecoVArrayMode(valor::conv::SuperkSample_t s) 
{
  VALORArray * arr = fPredictErecoSpectrumArray[s];
  if (!arr) {
    std::cerr<<"SKExp/CurrPredictedErecoVArrayMode could not find array for sample "<<util::AsString(s)<<std::endl;
    exit(-1);
  }
  
  return arr;
}
//________________________________________________________________________________________
const VALORArray * SuperkExpectation::CurrPredictedErecoVArray(valor::conv::SuperkSample_t s) const
{
  const VALORArray * arr = fPredictErecoSpectrumArrayTot[s];
  if (!arr) {
    std::cerr<<"SKExp/CurrPredictedErecoVArray could not find array for sample "<<util::AsString(s)<<std::endl;
    exit(-1);
  }
  
  return arr;
}
//________________________________________________________________________________________
const VALORArray * SuperkExpectation::CurrPredictedErecoVArrayMode(valor::conv::SuperkSample_t s) const
{
  const VALORArray * arr = fPredictErecoSpectrumArray[s];
  if (!arr) {
    std::cerr<<"SKExp/CurrPredictedErecoVArrayMode could not find array for sample "<<util::AsString(s)<<std::endl;
    exit(-1);
  }
  
  return arr;
}
//________________________________________________________________________________________
double* SuperkExpectation::CurrPredictedErecoSpectrumArray(valor::conv::SuperkSample_t s) 
{
  VALORArray * arr = this->CurrPredictedErecoVArray(s);
  return arr->GetArray();
}
//________________________________________________________________________________________
const double* SuperkExpectation::CurrPredictedErecoSpectrumArray(valor::conv::SuperkSample_t s)  const
{
  const VALORArray * arr = this->CurrPredictedErecoVArray(s);
  return arr->GetArray();
}
//________________________________________________________________________________________
double* SuperkExpectation::CurrPredictedErecoSpectrumArrayMode(
  valor::conv::SuperkSample_t s, valor::conv::SuperkMode_t m)
{ 
  int sidx = util::ArrIdx::Instance()->SuperkSampleIdx(s);
  if(sidx == -1) {
    return 0;
  }
  int midx = util::ArrIdx::Instance()->SuperkModeIdx(s,m);
  if(midx == -1) {
    return 0;
  }
  
  VALORArray * arr = this->CurrPredictedErecoVArrayMode(s);

  return &(arr->GetEntryRef(midx,0));
}
//________________________________________________________________________________________
const double* SuperkExpectation::CurrPredictedErecoSpectrumArrayMode(
  valor::conv::SuperkSample_t s, valor::conv::SuperkMode_t m) const
{ 
  int sidx = util::ArrIdx::Instance()->SuperkSampleIdx(s);
  if(sidx == -1) {
    return 0;
  }
  int midx = util::ArrIdx::Instance()->SuperkModeIdx(s,m);
  if(midx == -1) {
    return 0;
  }
  
  const VALORArray * arr = this->CurrPredictedErecoVArrayMode(s);

  return &(arr->GetEntryRef(midx,0));
}
//________________________________________________________________________________________
double SuperkExpectation::CurrPredictedErecoSpectrumArray(valor::conv::SuperkSample_t s, int i) const
{
  const double * arr = CurrPredictedErecoSpectrumArray(s);
  assert(arr);
  return arr[i];
}
//________________________________________________________________________________________












// NEW 2D TEMPLATES

TH2D* SuperkExpectation::CurrPredictedErecoVsQ2Spectrum(valor::conv::SuperkSample_t s) const 
{
  const int nbins_q2    = util::NumQ2Bins(s);
  const int nbins_ereco = util::NumErecoBins(s);
  typedef double (erq2type)  [nbins_ereco][nbins_q2];  

  const erq2type & arr = *reinterpret_cast< const erq2type * > (CurrPredictedErecoVsQ2SpectrumArray(s));  
  //const double * arr = CurrPredictedErecoSpectrumArray(s);
  if (!arr) return 0;
  
  string name (Form("curr_ereco_Vs_q2_%s"                       ,util::AsStringCompact(s).c_str() ) );
  string title(Form("Current Predicted Ereco Vs Q2 Spectrum: %s",util::AsString(s).c_str() ) );
  
  TH2D * h = new TH2D(name.c_str(), title.c_str(), nbins_q2, util::Q2BinEdges(s), nbins_ereco, util::ErecoBinEdges(s) );
  util::FormatErecoVsQ2Hist(h);
  for (int ir = 0 ; ir < nbins_ereco ; ir++ ) {
    for (int iq = 0 ; iq < nbins_q2 ; iq++ ) {
      h->SetBinContent(iq+1,ir+1,arr[ir][iq]);
    } 
  }
  return h;
}
//________________________________________________________________________________________
TH2D* SuperkExpectation::CurrPredictedErecoVsQ2Spectrum(
  valor::conv::SuperkSample_t s, valor::conv::SuperkMode_t m) const 
{ 
  const int nbins_q2    = util::NumQ2Bins(s);
  const int nbins_ereco = util::NumErecoBins(s);
  typedef double (erq2type)  [nbins_ereco][nbins_q2];  
  const erq2type & arr = *reinterpret_cast< const erq2type * > (CurrPredictedErecoVsQ2SpectrumArrayMode(s,m));  
  //const double * arr = CurrPredictedErecoSpectrumArrayMode(s,m);
  if (!arr) return 0;
  
  string name  (Form("curr_ereco_Vs_q2_%s_%s"                       ,util::AsStringCompact(s).c_str(), util::AsStringCompact(m).c_str()) );
  string title (Form("Current Predicted Ereco Vs Q2 Spectrum: %s %s",util::AsString(s).c_str()       , util::AsString(m).c_str()       ) );
  TH2D * h = new TH2D(name.c_str(), title.c_str(), nbins_q2, util::Q2BinEdges(s), nbins_ereco, util::ErecoBinEdges(s) );

  util::FormatErecoVsQ2Hist(h);
  util::SetDefaultSKSpectrumStyle(s,m,h);
  for (int ir = 0 ; ir < nbins_ereco ; ir++ ) {
    for (int iq = 0 ; iq < nbins_q2 ; iq++ ) {
      h->SetBinContent(iq+1,ir+1,arr[ir][iq]);
    }
  }

  return h;
}
//________________________________________________________________________________________
VALORArray * SuperkExpectation::CurrPredictedErecoVsQ2VArray(valor::conv::SuperkSample_t s) 
{
  VALORArray * arr = fPredictErecoVsQ2SpectrumArrayTot[s];
  if (!arr) {
    std::cerr<<"SKExp/CurrPredictedErecoVsQ2VArray could not find array for sample "<<util::AsString(s)<<std::endl;
    exit(-1);
  }
  
  return arr;
}
//________________________________________________________________________________________
VALORArray * SuperkExpectation::CurrPredictedErecoVsQ2VArrayMode(valor::conv::SuperkSample_t s) 
{
  VALORArray * arr = fPredictErecoVsQ2SpectrumArray[s];
  if (!arr) {
    std::cerr<<"SKExp/CurrPredictedErecoVsQ2VArrayMode could not find array for sample "<<util::AsString(s)<<std::endl;
    exit(-1);
  }
  
  return arr;
}
//________________________________________________________________________________________
const VALORArray * SuperkExpectation::CurrPredictedErecoVsQ2VArray(valor::conv::SuperkSample_t s) const
{
  const VALORArray * arr = fPredictErecoVsQ2SpectrumArrayTot[s];
  if (!arr) {
    std::cerr<<"SKExp/CurrPredictedErecoVsQ2VArray could not find array for sample "<<util::AsString(s)<<std::endl;
    exit(-1);
  }
  
  return arr;
}
//________________________________________________________________________________________
const VALORArray * SuperkExpectation::CurrPredictedErecoVsQ2VArrayMode(valor::conv::SuperkSample_t s) const
{
  const VALORArray * arr = fPredictErecoVsQ2SpectrumArray[s];
  if (!arr) {
    std::cerr<<"SKExp/CurrPredictedErecoVsQ2VArrayMode could not find array for sample "<<util::AsString(s)<<std::endl;
    exit(-1);
  }
  
  return arr;
}
//________________________________________________________________________________________
double* SuperkExpectation::CurrPredictedErecoVsQ2SpectrumArray(valor::conv::SuperkSample_t s) 
{
  VALORArray * arr = this->CurrPredictedErecoVsQ2VArray(s);
  return arr->GetArray();
}
//________________________________________________________________________________________
const double* SuperkExpectation::CurrPredictedErecoVsQ2SpectrumArray(valor::conv::SuperkSample_t s)  const
{
  const VALORArray * arr = this->CurrPredictedErecoVsQ2VArray(s);
  return arr->GetArray();
}
//________________________________________________________________________________________
double* SuperkExpectation::CurrPredictedErecoVsQ2SpectrumArrayMode(
  valor::conv::SuperkSample_t s, valor::conv::SuperkMode_t m)
{ 
  int sidx = util::ArrIdx::Instance()->SuperkSampleIdx(s);
  if(sidx == -1) {
    return 0;
  }
  int midx = util::ArrIdx::Instance()->SuperkModeIdx(s,m);
  if(midx == -1) {
    return 0;
  }
  
  VALORArray * arr = this->CurrPredictedErecoVsQ2VArrayMode(s);

  return &(arr->GetEntryRef(midx,0,0));
}
//________________________________________________________________________________________
const double* SuperkExpectation::CurrPredictedErecoVsQ2SpectrumArrayMode(
  valor::conv::SuperkSample_t s, valor::conv::SuperkMode_t m) const
{ 
  int sidx = util::ArrIdx::Instance()->SuperkSampleIdx(s);
  if(sidx == -1) {
    return 0;
  }
  int midx = util::ArrIdx::Instance()->SuperkModeIdx(s,m);
  if(midx == -1) {
    return 0;
  }
  
  const VALORArray * arr = this->CurrPredictedErecoVsQ2VArrayMode(s);

  return &(arr->GetEntryRef(midx,0,0));
}
//________________________________________________________________________________________
double SuperkExpectation::CurrPredictedErecoVsQ2SpectrumArray(valor::conv::SuperkSample_t s, int ir, int iq) const
{
  const int nbins_q2    = util::NumQ2Bins(s);
  const int nbins_ereco = util::NumErecoBins(s);
  typedef double (erq2type)  [nbins_ereco][nbins_q2];  
  const erq2type & arr = *reinterpret_cast< const erq2type * > (CurrPredictedErecoVsQ2SpectrumArray(s));  
  //const double * arr = CurrPredictedErecoVsQ2SpectrumArray(s);
  if (!arr) return 0;

  assert(arr);
  return arr[ir][iq];
} 
//














//________________________________________________________________________________________
void SuperkExpectation::PrintSystematicVariations(std::ostream & stream) const
{
  cout << "Exp >> SuperK PDF tweaks: " << endl;
  
#ifdef __FHC_RUNNING__ // NEW FLUX
  fBANFF             -> PrintSystematicVariations(stream);
#endif
#ifdef __RHC_RUNNING__
  fBANFF             -> PrintRHCSystematicVariations(stream);
#endif
#ifdef __FHC_RUNNING__ // NEW FLUX
  fNuXSec            -> PrintSystematicVariations(stream);
#endif
#ifdef __RHC_RUNNING__
  fNuXSec            -> PrintRHCSystematicVariations(stream);
#endif
#ifdef __FHC_RUNNING__
  fSuperkDetectorFSI -> PrintSystematicVariations(stream);
#endif
#ifdef __RHC_RUNNING__
  fSuperkDetectorFSI -> PrintRHCSystematicVariations(stream);
#endif
  /// NEW MULTIRING 2
#ifdef __MULTIRING_FHC_RUNNING__
  fSuperkDetectorFSI -> PrintMultiRingSystematicVariations(stream);
#endif  
  ///
}
//________________________________________________________________________________________
bool SuperkExpectation::Write(TObjArray * array, const char * pdf_tag, Option_t * write_opt)
{
  if(!array) {
    cout << "SKExp >> Null TObjArray. Not writing data" << endl;
    return false;
  }
  
  bool hist_mode = (write_opt && string(write_opt).find("hist")!=string::npos);
  //bool per_mode  = (write_opt && string(write_opt).find("mode")!=string::npos);
  
  string append_begin(Form("%s_", pdf_tag));
  
  for(int is=0; is<kns; is++) {
    SuperkSample_t s = kSKSamples[is];
    if(s == kSKsUndefined) continue;
    

#ifndef __USE_2D_TEMPLATES__ // 1D templates

    if (hist_mode) {
      util::SetDefaultStyle();
      TH1D * spectrum = 0;
      TH1D * spectrum_mode[knm] = { };
      
      // copy & format Ereco spectra
      spectrum = this->CurrPredictedErecoSpectrum(s);
      spectrum->SetDirectory (0);
      
      for(int im = 0; im < knm; im++) {
         SuperkMode_t m = kSKReactionModes[is][im];
         if(m == kSKmUndefined) continue;
         spectrum_mode[im] = this->CurrPredictedErecoSpectrum(s,m);
         spectrum_mode[im]->SetDirectory (0);
      }//imode
      
      // add spectra to the results array
      spectrum->SetName(util::AsStringCompact(s,kSKmUndefined,append_begin).c_str());
      array->Add(spectrum);
      for(int im = 0; im < knm; im++) {
        SuperkMode_t m = kSKReactionModes[is][im];
         if(m == kSKmUndefined) continue;
         spectrum_mode[im]->SetName(util::AsStringCompact(s,m,append_begin).c_str());
         array->Add(spectrum_mode[im]);
      }//im
    } // hist mode
    else {
      VALORArray * total = new VALORArray(*(this->CurrPredictedErecoVArray(s)));
      total->SetName(Form("%s%s",append_begin.c_str(),total->GetName()));
      array->Add(total);
      VALORArray * bymode = new VALORArray(*(this->CurrPredictedErecoVArrayMode(s)));
      bymode->SetName(Form("%s%s",append_begin.c_str(),bymode->GetName()));
      array->Add(bymode);
    }
  }//is

#else // 2D templates
  
  // NEW 2D TEMPLATES
  
  if (hist_mode) {
      util::SetDefaultStyle();
      TH2D * spectrum = 0;
      TH2D * spectrum_mode[knm] = { };
      
      // copy & format Ereco spectra
      spectrum = this->CurrPredictedErecoVsQ2Spectrum(s);
      spectrum->SetDirectory (0);
      
      for(int im = 0; im < knm; im++) {
         SuperkMode_t m = kSKReactionModes[is][im];
         if(m == kSKmUndefined) continue;
         spectrum_mode[im] = this->CurrPredictedErecoVsQ2Spectrum(s,m);
         spectrum_mode[im]->SetDirectory (0);
      }//imode
      
      // add spectra to the results array
      spectrum->SetName(util::AsStringCompact(s,kSKmUndefined,append_begin).c_str());
      array->Add(spectrum);
      for(int im = 0; im < knm; im++) {
         SuperkMode_t m = kSKReactionModes[is][im];
         if(m == kSKmUndefined) continue;
         spectrum_mode[im]->SetName(util::AsStringCompact(s,m,append_begin).c_str());
         array->Add(spectrum_mode[im]);
      }//im
    } // hist mode
    else {
      VALORArray * total = new VALORArray(*(this->CurrPredictedErecoVsQ2VArray(s)));
      total->SetName(Form("%s%s",append_begin.c_str(),total->GetName()));
      array->Add(total);
      VALORArray * bymode = new VALORArray(*(this->CurrPredictedErecoVsQ2VArrayMode(s)));
      bymode->SetName(Form("%s%s",append_begin.c_str(),bymode->GetName()));
      array->Add(bymode);
    }
  }//is

 //

#endif

  return true;
}
//________________________________________________________________________________________
bool SuperkExpectation::Write(
  const char * filename, const char * pdf_tag, Option_t * write_opt, Option_t * file_opt)
{
  if (!file_opt) file_opt = "update";
  TObjArray * array = new TObjArray();
  if (write_opt && string(write_opt).find("hist")!=string::npos) array->SetOwner(true);
  this->Write(array, pdf_tag, write_opt);
  TFile f(filename, file_opt);
  array->Write();
  f.Close();

  return true;
}
//________________________________________________________________________________________
void SuperkExpectation::ResetSystematics(void)
{
#ifdef __show_syst_reset_mesg__
  cout << "Exp >> Reseting all systematic variations." << endl;
#endif

  // NEW FHC+RHC JOINT

  //
  // FHC+RHC 
  //
  
  // BANFF
#if defined(__FHC_RUNNING__) && defined(__RHC_RUNNING__) 
  this->GetBANFF()->TweakFHCpRHC(true,0);
  // SK det + FSI
  this->GetSuperkDetectorFSI()->TweakFHCpRHC(true,0);
  this->GetSuperkDetectorFSI()->TweakEnergyScale(0.);
  this->GetSuperkDetectorFSI()->TweakRHCEnergyScale(0.);
#endif
  //


  
  //
  //FHC
  //
  

#if defined(__FHC_RUNNING__) && !defined(__RHC_RUNNING__) // NEW FHC+RHC JOINT
  // BANFF
  this->GetBANFF()->Tweak(true,0);
  // SK det + FSI
  this->GetSuperkDetectorFSI()->TweakEnergyScale(0.);
  this->GetSuperkDetectorFSI()->Tweak(true,0);
#endif

  // X Section

  // NEW XSEC --> comment xsec syst that now are in banff or not used any longer
  // --> not commented any longer

  // NEW XSEC NIWG
#if defined(__USE_NIWG2013__)
  this->GetNuXSec()->Tweak(kNuisXSecFermiMom      , 0.);
  this->GetNuXSec()->Tweak(kNuisXSecSpecFunc      , 0.);
  this->GetNuXSec()->Tweak(kNuisXSecWshape        , 0.);
  this->GetNuXSec()->Tweak(kNuisXSecShapeCCoth    , 0.);
  this->GetNuXSec()->Tweak(kNuisXSecPiLessDelta   , 0.);
  #ifndef __NUEONLY__
  this->GetNuXSec()->Tweak(kNuisXSecBindE         , 0.);
  #endif
  this->GetNuXSec()->Tweak(kNuisXSecNormCCcoh     , 0.);
  #ifndef __NUEONLY__
  this->GetNuXSec()->Tweak(kNuisXSecNormNC1piC    , 0.);
  #endif
  #ifndef __NUMUONLY__
  this->GetNuXSec()->Tweak(kNuisXSecNormNCcoh     , 0.);
  #endif
  this->GetNuXSec()->Tweak(kNuisXSecNormNCoth     , 0.);
  this->GetNuXSec()->Tweak(kNuisXSecNormNuEToNuMu , 0.);
  #ifndef __NUEONLY__
  this->GetNuXSec()->Tweak(kNuisXSecNormAntiNuToNu, 0.);
  #endif
  //RHC
  this->GetNuXSec()->Tweak(kNuisRHCToFHC, 0.);
#elif defined(__USE_NIWG2015__)
  this->GetNuXSec()->Tweak(kNuisXSecRPA , 1.); // NEW XSEC 
#endif



  //
  // RHC
  //

#if !defined(__FHC_RUNNING__) && defined(__RHC_RUNNING__) // NEW FHC+RHC JOINT
  // BANFF
  this->GetBANFF()->TweakRHC(true,0);
  // SK det + FSI
  this->GetSuperkDetectorFSI()->TweakRHCEnergyScale(0.);
  this->GetSuperkDetectorFSI()->TweakRHC(true,0);
#endif

  // X Section

  // NEW XSEC --> comment xsec syst that now are in banff or not used any longer
  // --> not commented any longer

  // NEW XSEC NIWG
#if defined(__USE_NIWG2013__)
  this->GetNuXSec()->Tweak(kRHCNuisXSecFermiMom      , 0.);
  this->GetNuXSec()->Tweak(kRHCNuisXSecSpecFunc      , 0.);
  this->GetNuXSec()->Tweak(kRHCNuisXSecWshape        , 0.);
  this->GetNuXSec()->Tweak(kRHCNuisXSecShapeCCoth    , 0.);
  this->GetNuXSec()->Tweak(kRHCNuisXSecPiLessDelta   , 0.);
  #ifndef __NUEONLY__
  this->GetNuXSec()->Tweak(kRHCNuisXSecBindE         , 0.);
  #endif
  this->GetNuXSec()->Tweak(kRHCNuisXSecNormCCcoh     , 0.);
  #ifndef __NUEONLY__
  this->GetNuXSec()->Tweak(kRHCNuisXSecNormNC1piC    , 0.);
  #endif
  #ifndef __NUMUONLY__
  this->GetNuXSec()->Tweak(kRHCNuisXSecNormNCcoh     , 0.);
  #endif
  this->GetNuXSec()->Tweak(kRHCNuisXSecNormNCoth     , 0.);
  this->GetNuXSec()->Tweak(kRHCNuisXSecNormNuEToNuMu , 0.);
  #ifndef __NUEONLY__
  this->GetNuXSec()->Tweak(kRHCNuisXSecNormAntiNuToNu, 0.);
  #endif
#elif defined(__USE_NIWG2015__)
  this->GetNuXSec()->Tweak(kRHCNuisXSecRPA           , 1.); // NEW XSEC
#endif

  /// NEW MULTIRING 2
  //
  //FHC MULTIRING
  //
  // SK det + FSI
#ifdef __MULTIRING_FHC_RUNNING__
  /// NEW MULTIRING 3
  this->GetSuperkDetectorFSI()->TweakMultiRingEnergyScale(0.);   
  ///
  this->GetSuperkDetectorFSI()->TweakMultiRing(true,0);
#endif
  ///

}
//________________________________________________________________________________________
void SuperkExpectation::TweakSystematic(string systematic, double twkdial, Option_t * opt)
{
  bool best_banff = !opt || string(opt).find("best_banff")!=string::npos;
  TweakSystematicLegacy(systematic, twkdial, best_banff);
}
//________________________________________________________________________________________
void SuperkExpectation::TweakSystematicLegacy(string systematic, double twkdial, bool best_banff)
{

#ifdef __show_syst_tweak_mesg__
  cout << "Exp >> Setting systematic parameter: " << systematic 
       << " to value = " << twkdial << endl;
#endif

 {
    PhysParamId_t p = util::PhysParamIdFromString(systematic);
    if (p!=kPhypUndefined) {

      //cout << util::AsString(p) << endl;
      if(best_banff) this->TweakSystematic(p, twkdial,"best_banff");
      else this->TweakSystematic(p, twkdial);
      return;
    }
 }
  
//
// Valid systematic parameter names:
// 
// SK detector systematics
//    - sk_energy_scale                     : 
//    - sk_1Rmu_ccqe_ring_counting_eff      : tweak SK 1Rmu/CCQE RC efficiencies in all 3 Ereco bins, correlations taken into account (Cholesky method)
//      or
//      * sk_1Rmu_ccqe_ring_counting_eff_0  : tweak only SK 1Rmu/CCQE RC efficiency in the 1st Ereco bin, correlations with other bins ignored (debug mode)
//      * sk_1Rmu_ccqe_ring_counting_eff_1  : tweak only SK 1Rmu/CCQE RC efficiency in the 2nd Ereco bin, correlations with other bins ignored (debug mode)
//      * sk_1Rmu_ccqe_ring_counting_eff_2  : tweak only SK 1Rmu/CCQE RC efficiency in the 3rd Ereco bin, correlations with other bins ignored (debug mode)
//   - sk_1Rmu_ccqe_other_eff              : 
//   - sk_1Rmu_ccnonqe_eff                  :
//   - sk_1Rmu_nc_eff                       :
//   - sk_1Rmu_nuecc_eff                    :
//   - sk_1Re_eff                           :
//   - sk_1Re_eff_...                       :
//
// BANFF
//    - banff                               : tweak all BANFF params, correlations taken into account (Cholesky method)
//      or   
//      * banff_00                          : tweak only BANFF param 00, see Conventions.h for definition, correlations ignored (debug mode)
//      * banff_01                          : tweak only BANFF param 01, see Conventions.h for definition, correlations ignored (debug mode)
//      * ..     
//      * banff_30(48)                      : tweak only BANFF param 35(48), see Conventions.h for definition, correlations ignored (debug mode)
//      * banff_maqe                        : tweak only BANFF maqe parameter
//      * banff_mares                       : tweak only BANFF mares parameter
//      * banff_ccqeE0                      : tweak only BANFF ccqeE0 parameter
//      * banff_ccqeE1                      : tweak only BANFF ccqeE1 parameter
//      * banff_ccqeE2                      : tweak only BANFF ccqeE2 parameter
//      * banff_cc1pi0E0                    : tweak only BANFF cc1piE0 parameter
//      * banff_cc1pi0E1                    : tweak only BANFF cc1piE1 parameter
//      * banff_nc1pi0                      : tweak only BANFF nc1pi0 parameter
//      * banff_mecnorm                     : tweak only BANFF MEC norm parameter
//      * banff_pdd                         : tweak only BANFF PDD double counting parameter
//
// FSI
//    - fsi                                 : tweak all FSI params, correlations taken into account (Cholesky method)
//      or     
//      * fsi_000                           : tweak only FSI param 00, see Conventions.h for definition, correlations ignored (debug mode)
//      * fsi_001                           : tweak only FSI param 01, see Conventions.h for definition, correlations ignored (debug mode)
//      * ..     
//      * fsi_014                           : tweak only FSI param 14, see Conventions.h for definition, correlations ignored (debug mode)
//
// XSec
//    - xsec_pf
//    - xsec_sf
//    - xsec_wshape
//    - xsec_shape_ccoth
//    - xsec_pilessdelta
//    - xsec_bindE
//    - norm_cccoh
//    - norm_nc1piC
//    - norm_ncoth
//    - norm_nuetonumu
//    - norm_antinutonu
//


  //
  // SK detector systematics
  //
  
  if(systematic == "sk_e_scale") {
    this->GetSuperkDetectorFSI()->TweakEnergyScale(twkdial);
  }
  else if(systematic == "sk_e_scale_rhc") {
    this->GetSuperkDetectorFSI()->TweakRHCEnergyScale(twkdial);
  }

  // NEW FHC+RHC JOINT
  else if(systematic == "sk_e_scale_fhcprhc") {
    this->GetSuperkDetectorFSI()->TweakEnergyScale(twkdial);
    this->GetSuperkDetectorFSI()->TweakRHCEnergyScale(twkdial);
  }
  //

  /// NEW MULTIRING 3
  else if(systematic == "sk_e_scale_multiring") {
    this->GetSuperkDetectorFSI()->TweakMultiRingEnergyScale(twkdial);
  }
  ///

  //
  // SKDET + FSI
  //

  // tweak all FSI parameters taking into account correlations
  else if(systematic == "skdetfsi") {
    this->GetSuperkDetectorFSI()->Tweak(); //< Cholesky method
  }
  else if(systematic == "skdetfsi_rhc") {
    this->GetSuperkDetectorFSI()->TweakRHC(); //< Cholesky method
  }

  // NEW FHC+RHC 
  else if(systematic == "skdetfsi_fhcprhc") {
    this->GetSuperkDetectorFSI()->TweakFHCpRHC(); //< Cholesky method
  }
  //

  /// NEW MULTIRING 2
  else if(systematic == "skdetfsi_multiring") {
    this->GetSuperkDetectorFSI()->TweakMultiRing(); //< Cholesky method
  }
  ///
  // tweak an individual FSI parameter (ignores correlations, debug)
  else if(systematic.find("skdetfsi") != string::npos) {

    // NEW FHC+RHC JOINT
    if(systematic.find("rhc")!=string::npos && systematic.find("fhc")!=string::npos){
      cerr << endl;
      cerr << "SuperkExpectation::TweakSystematicLegacy" << endl;
      cerr << "Tweak FHC and RHC SKDetFSI systematic separately, if w/o correlations" << endl;
      exit(1);
    }
    //
    
    bool is_rhc = systematic.find("rhc") != string::npos;
    /// NEW MULTIRING 2
    bool is_multiring = systematic.find("multiring") != string::npos;
    ///
    bool found = false;
    int nskdetfsi = 0;
    if(is_rhc)
      nskdetfsi = conv::kNumOfRHCSKDetFSIParams;
    /// NEW MULTIRING 2
    else if(is_multiring)
      nskdetfsi = conv::kNumOfSKDetFSIParamsMultiRing;
    ///
    else
      nskdetfsi = conv::kNumOfSKDetFSIParams;
    TVectorD twkv(nskdetfsi);
    for (int ip = 0; ip < nskdetfsi; ip++) {
      twkv[ip] = 1;
    }
    for (int ip = 0; ip < nskdetfsi; ip++) {
      if(is_rhc) {
        if(systematic == Form("skdetfsi%03i_rhc",ip)) {
          twkv[ip]=twkdial;
          this->GetSuperkDetectorFSI()->TweakRHC(true,&twkv);
          found = true;
          break;
        }
      }//is_rhc


      /// NEW MULTIRING 2
      if(is_multiring) {
        if(systematic == Form("skdetfsi%03i_multiring",ip)) {
          twkv[ip]=twkdial;
          this->GetSuperkDetectorFSI()->TweakMultiRing(true,&twkv);
          found = true;
          break;
        }
      }//is_multiring
      ///


      else {
        if(systematic == Form("skdetfsi%03i",ip)) {
          twkv[ip]=twkdial;
          this->GetSuperkDetectorFSI()->Tweak(true,&twkv);
          found = true;
          break;
        }
      }//!is_rhc
    }//ip
     if(!found) {
       cerr << "Exp >> Unknown systematic parameter (" << systematic << ")" << endl;
       exit(1);
     }   
  }//skdetfsi
    
  //
  // BANFF & Uncorrelated BANFF
  //

  // tweak all BANFF parameters taking into account correlations
  else
  if(systematic == "banff") {
    this->GetBANFF()->Tweak(false,0); //< Cholesky method
  }
  else
  if(systematic == "banff_rhc") {
    this->GetBANFF()->TweakRHC(false,0); //< Cholesky method
  }
  
  // NEW FHC+RHC JOINT
  else
  if(systematic == "banff_fhc_rhc") {
    this->GetBANFF()->TweakFHCpRHC(false,0); //< Cholesky method
  }
  //

  // tweak an individual BANFF parameter (ignores correlations, debug)
  else
  if(systematic.find("banff") != string::npos) {

    // NEW FHC+RHC JOINT 
    if(systematic.find("rhc")!=string::npos && systematic.find("fhc")!=string::npos){
      cerr << endl;
      cerr << "SuperkExpectation::TweakSystematicLegacy" << endl;
      cerr << "Tweak FHC and RHC BANFF systematic separately, if w/o correlations" << endl;
      exit(1);
    }
    
    // NEW FHC+RHC JOINT (NEW)
    
    else if(systematic.find("fhc_rhc")!=string::npos){
      
      const int nbanff = conv::kNumOfFHCpRHCParamsInBANFFfit;
      TVectorD twkv(this->GetBANFF()->FHCpRHCBestFitValues());
      if(!best_banff) {
        for (int ib = 0; ib < nbanff; ib++) { twkv[ib] = 1; }
      }
      PhysParamId_t p = util::PhysParamIdFromString(systematic);
      int ip = util::ArrIdx::Instance()->FHCpRHCBANFFSystIdx(p);
      if(ip < 0) {
        cerr << "Exp >> Unknown systematic parameter (" << systematic << ")" << endl;
        exit(1);
      }
      twkv[ip] = twkdial;
#ifdef __show_syst_tweak_debug__
      cout << "Exp >> Tweaking systematic string: " << systematic
           << ", NuisParam: " << util::AsString(p)
           << endl;
      util::PrintCompact(twkv);
#endif
      this->GetBANFF()->TweakFHCpRHC(true,&twkv);

      //

    }
    
    else {

      bool is_rhc = systematic.find("rhc") != string::npos;
      if(is_rhc) {
	const int nbanff = conv::kNumOfRHCParamsInBANFFfit;
	TVectorD twkv(this->GetBANFF()->RHCBestFitValues());
	if(!best_banff) {
	  for (int ib = 0; ib < nbanff; ib++) { twkv[ib] = 1; }
	}
	PhysParamId_t p = util::PhysParamIdFromString(systematic);
	int ip = util::ArrIdx::Instance()->RHCBANFFSystIdx(p);
	if(ip < 0) {
	  cerr << "Exp >> Unknown systematic parameter (" << systematic << ")" << endl;
	  exit(1);
	}
	twkv[ip] = twkdial;
#ifdef __show_syst_tweak_debug__
	cout << "Exp >> Tweaking systematic string: " << systematic
	     << ", NuisParam: " << util::AsString(p)
	     << endl;
	util::PrintCompact(twkv);
#endif
	this->GetBANFF()->TweakRHC(true,&twkv);
      }//is_rhc
      else {
	const int nbanff = conv::kNumOfParamsInBANFFfit;
	TVectorD twkv(this->GetBANFF()->BestFitValues());
	if(!best_banff) {
	  for (int ib = 0; ib < nbanff; ib++) { twkv[ib] = 1; }
	}
	
	PhysParamId_t p = util::PhysParamIdFromString(systematic);
	int ip = util::ArrIdx::Instance()->BANFFSystIdx(p);
	if(ip < 0) {
	  cerr << "Exp >> Unknown systematic parameter (" << systematic << ")" << endl;
	  exit(1);
	}
	twkv[ip] = twkdial;
	
#ifdef __show_syst_tweak_debug__
	cout << "Exp >> Tweaking systematic string: " << systematic
	     << ", NuisParam: " << util::AsString(p)
	     << endl;
	util::PrintCompact(twkv);
#endif
	this->GetBANFF()->Tweak(true,&twkv);
      }//!is_rhc

    }
  }//banff
  
  

  //
  // XSec
  //

  // NEW XSEC --> comment xsec syst that now are in banff or not used any longer
  // --> not commented any longer

  // NEW XSEC NIWG

#if defined(__USE_NIWG2013__)

  else if(systematic == "pfermi")
    this->GetNuXSec()->Tweak(kNuisXSecFermiMom      , twkdial);  
  else if(systematic == "specfunc")
    this->GetNuXSec()->Tweak(kNuisXSecSpecFunc      , twkdial);
  else if(systematic == "wshape")
      this->GetNuXSec()->Tweak(kNuisXSecWshape        , twkdial);
    else if(systematic == "shapeccoth")
      this->GetNuXSec()->Tweak(kNuisXSecShapeCCoth    , twkdial);
    else if(systematic == "pilessdelta")
      this->GetNuXSec()->Tweak(kNuisXSecPiLessDelta   , twkdial);
    #ifndef __NUEONLY__
    else if(systematic == "bindE")
      this->GetNuXSec()->Tweak(kNuisXSecBindE         , twkdial);
    #endif
    else if(systematic == "norm_cccoh")
      this->GetNuXSec()->Tweak(kNuisXSecNormCCcoh     , twkdial);
    #ifndef __NUMUONLY__
    else if(systematic == "norm_nccoh")
      this->GetNuXSec()->Tweak(kNuisXSecNormNCcoh     , twkdial);
    #endif
    #ifndef __NUEONLY__
    else if(systematic == "norm_nc1piC")
      this->GetNuXSec()->Tweak(kNuisXSecNormNC1piC    , twkdial);
    #endif
    else if(systematic == "norm_ncoth")
      this->GetNuXSec()->Tweak(kNuisXSecNormNCoth     , twkdial);
    else if(systematic == "norm_nuetonumu")
      this->GetNuXSec()->Tweak(kNuisXSecNormNuEToNuMu , twkdial);
    #ifndef __NUEONLY__
    else if(systematic == "norm_antinutonu")
      this->GetNuXSec()->Tweak(kNuisXSecNormAntiNuToNu, twkdial);
    #endif
    else if(systematic == "pfermi_rhc")
      this->GetNuXSec()->Tweak(kRHCNuisXSecFermiMom      , twkdial);
    else if(systematic == "specfunc_rhc")
      this->GetNuXSec()->Tweak(kRHCNuisXSecSpecFunc      , twkdial);
    else if(systematic == "wshape_rhc")
      this->GetNuXSec()->Tweak(kRHCNuisXSecWshape        , twkdial);
    else if(systematic == "shapeccoth_rhc")
      this->GetNuXSec()->Tweak(kRHCNuisXSecShapeCCoth    , twkdial);
    else if(systematic == "pilessdelta_rhc")
      this->GetNuXSec()->Tweak(kRHCNuisXSecPiLessDelta   , twkdial);
    #ifndef __NUEONLY__
    else if(systematic == "bindE_rhc")
      this->GetNuXSec()->Tweak(kRHCNuisXSecBindE         , twkdial);
    #endif
    else if(systematic == "norm_cccoh_rhc")
      this->GetNuXSec()->Tweak(kRHCNuisXSecNormCCcoh     , twkdial);
    #ifndef __NUMUONLY__
    else if(systematic == "norm_nccoh_rhc")
      this->GetNuXSec()->Tweak(kRHCNuisXSecNormNCcoh     , twkdial);
    #endif
    #ifndef __NUEONLY__
    else if(systematic == "norm_nc1piC_rhc")
      this->GetNuXSec()->Tweak(kRHCNuisXSecNormNC1piC    , twkdial);
    #endif
    else if(systematic == "norm_ncoth_rhc")
      this->GetNuXSec()->Tweak(kRHCNuisXSecNormNCoth     , twkdial);
    else if(systematic == "norm_nuetonumu_rhc")
      this->GetNuXSec()->Tweak(kRHCNuisXSecNormNuEToNuMu , twkdial);
    #ifndef __NUEONLY__
    else if(systematic == "norm_antinutonu_rhc")
      this->GetNuXSec()->Tweak(kRHCNuisXSecNormAntiNuToNu, twkdial);
    #endif
    else if(systematic == "rhcfhcnorm")
      this->GetNuXSec()->Tweak(kNuisRHCToFHC, twkdial);
  
#elif defined(__USE_NIWG2015__)
    // NEW XSEC
    else if(systematic == "rpa")
      this->GetNuXSec()->Tweak(kNuisXSecRPA      , twkdial);
    else if(systematic == "rpa_rhc")
      this->GetNuXSec()->Tweak(kRHCNuisXSecRPA   , twkdial); 
#endif
  
  else {
     cerr << "Exp >> Unknown systematic parameter (" << systematic << ")" << endl;
     exit(1);
  }   
}
//________________________________________________________________________________________
void SuperkExpectation::TweakSystematic(conv::PhysParamId_t p, double twkdial, Option_t * opt)
{

  cout << "SuperkExpectation::TweakSystematic(p) 1" << endl;

      
#ifdef __show_syst_tweak_mesg__
  cout << "Exp >> Setting systematic parameter: " << util::AsString(p) 
       << " to value = " << twkdial << endl;
#endif
    
  if(p == kNuisSuperkEnergyScale) {
    this->GetSuperkDetectorFSI()->TweakEnergyScale(twkdial);
  }
  /// NEW MULTIRING 3
  else if(p == kMultiRingNuisSuperkEnergyScale) {
    this->GetSuperkDetectorFSI()->TweakMultiRingEnergyScale(twkdial);
  }
  // ///
  else if(p == kRHCNuisSuperkEnergyScale) {
    this->GetSuperkDetectorFSI()->TweakRHCEnergyScale(twkdial);
  }

  //
  // SKDET + FSI
  //

  // tweak an individual FSI parameter (ignores correlations, debug)
  else if(util::IsSKDetFSI(p)) {

    bool is_rhc = util::IsRHCSystematic(p);
    bool is_fhc = util::IsFHCSystematic(p);
    bool is_multiring = util::IsMultiRingSystematic(p);
   

#if defined(__FHC_RUNNING__) && defined(__RHC_RUNNING__) 

    // NEW 2D TEMPLATES (NEW)

    const int nskdetfsi = conv::kNumOfFHCpRHCSKDetFSIParams;
    TVectorD twkv(nskdetfsi);
    for (int ib = 0; ib < nskdetfsi; ib++) { twkv[ib] = 1; }
  
    int ip = util::ArrIdx::Instance()->FHCpRHCSKDetFSISystIdx(p);
    if(ip < 0) {
      cerr << "Exp >> Unknown systematic parameter (" << util::AsString(p) << ")" << endl;
      exit(1);
    }
    twkv[ip] = twkdial;

    cout << "ip = " << ip << " - twkdial = " << twkdial << " - twkv[ip] = " << twkv[ip] << endl;
    this->GetSuperkDetectorFSI()->TweakFHCpRHC(true,&twkv);  

    //
#else

    int nskdetfsi = 0;
    if(is_rhc)
      nskdetfsi = conv::kNumOfRHCSKDetFSIParams;
    else if(is_fhc)
      nskdetfsi = conv::kNumOfSKDetFSIParams;
    else
      nskdetfsi = conv::kNumOfSKDetFSIParamsMultiRing;
    ///

    TVectorD twkv(nskdetfsi);
    for (int ip = 0; ip < nskdetfsi; ip++) {
      twkv[ip] = 1;
    }
    int ip = 0;
    
    if (is_rhc)      ip = p-kListOfRHCSKDetFSIParams[0];
    else if (is_fhc) ip = p-kListOfSKDetFSIParams[0];
    else             ip = p-kListOfMultiRingSKDetFSIParams[0];
    
    assert(ip >=0 && ip<twkv.GetNoElements());
    twkv[ip]=twkdial;

    if (is_rhc)      this->GetSuperkDetectorFSI()->TweakRHC(true,&twkv);
    else if (is_fhc) this->GetSuperkDetectorFSI()->Tweak(true,&twkv);
    else             this->GetSuperkDetectorFSI()->TweakMultiRing(true,&twkv);
    
#endif


  }//skdetfsi
    
  //
  // BANFF & Uncorrelated BANFF
  //
  
  // tweak an individual BANFF parameter (ignores correlations, debug)
  else if(util::IsBANFF(p)) {

#if defined(__FHC_RUNNING__) && defined(__RHC_RUNNING__) 
    
    // NEW FHC+RHC JOINT (NEW)

    bool best_banff = opt && string(opt).find("best_banff")!=string::npos;
    const int nbanff = conv::kNumOfFHCpRHCParamsInBANFFfit;
    TVectorD twkv(this->GetBANFF()->FHCpRHCBestFitValues());
    if(!best_banff) {
      for (int ib = 0; ib < nbanff; ib++) { twkv[ib] = 1; }
    }
    int ip = util::ArrIdx::Instance()->FHCpRHCBANFFSystIdx(p);
    if(ip < 0) {
      cerr << "Exp >> Unknown systematic parameter (" << util::AsString(p) << ")" << endl;
      exit(1);
    }
    twkv[ip] = twkdial;
    this->GetBANFF()->TweakFHCpRHC(true,&twkv);  
    //
    
#else
  
    bool best_banff = opt && string(opt).find("best_banff")!=string::npos;
    bool is_rhc = util::IsRHCSystematic(p);
    if(is_rhc) {
      const int nbanff = conv::kNumOfRHCParamsInBANFFfit;
      TVectorD twkv(this->GetBANFF()->RHCBestFitValues());
      if(!best_banff) {
        for (int ib = 0; ib < nbanff; ib++) { twkv[ib] = 1; }
      }
      int ip = util::ArrIdx::Instance()->RHCBANFFSystIdx(p);
      if(ip < 0) {
        cerr << "Exp >> Unknown systematic parameter (" << util::AsString(p) << ")" << endl;
        exit(1);
      }
      twkv[ip] = twkdial;
      this->GetBANFF()->TweakRHC(true,&twkv);
    }//is_rhc
    else {
      const int nbanff = conv::kNumOfParamsInBANFFfit;
      TVectorD twkv(this->GetBANFF()->BestFitValues());
      if(!best_banff) {
        for (int ib = 0; ib < nbanff; ib++) { twkv[ib] = 1; }
      }
      int ip = util::ArrIdx::Instance()->BANFFSystIdx(p);
      if(ip < 0) {
        cerr << "Exp >> Unknown systematic parameter (" << util::AsString(p) << ")" << endl;
        exit(1);
      }
      twkv[ip] = twkdial;
      this->GetBANFF()->Tweak(true,&twkv);
    }//!is_rhc

#endif // if defined(__FHC_RUNNING__) && defined(__RHC_RUNNING__) 

  }//banff


  //
  // XSec
  //

  // NEW XSEC --> attention to new xsec non-banff systematics. not changed 
  else if (util::IsXSec(p)) {
    this->GetNuXSec()->Tweak(p , twkdial);
  }
  
  else {
     cerr << "Exp/Tweak>> Unhandled systematic parameter (" << util::AsString(p) << ")" << endl;
     exit(1);
  }           
}


//________________________________________________________________________________________
void SuperkExpectation::RandomizeSystematics(void)
{
#ifdef __show_syst_tweak_mesg__
  cout << "Exp >> Generating random variations of all systematic params." << endl;
#endif

//   // common to all the samples
//   this->TweakAllBANFFSystematics();
//   this->TweakAllUncorrelatedXSecSystematics();

// #ifdef __FHC_RUNNING__
//   this->TweakAllSKDetectorFSISystematics();
// #endif
  
// #ifdef __RHC_RUNNING__   /// NEW MULTIRING 3
//   this->GetNuXSec()->Tweak(kNuisRHCToFHC, Rnd::Instance()->RndToyMCsyst().Gaus(0.,1.));
//   #ifdef __FULLY_CORRELATED__
//   this->GetNuXSec()->SyncRHCToFHCSystematics();
//   this->GetBANFF() ->SyncRHCToFHCSystematics();
//   this->GetSuperkDetectorFSI()->SyncRHCToFHCSystematics();   
//   #else
//   #ifdef __BANFF_UNCORRELATED__
//   this->TweakAllRHCBANFFSystematics();
//   this->GetNuXSec()->SyncRHCToFHCSystematics();
//   this->GetSuperkDetectorFSI()->SyncRHCToFHCSystematics();
//   #else
//   this->TweakAllRHCBANFFSystematics();
//   this->TweakAllRHCUncorrelatedXSecSystematics();
//   this->TweakAllRHCSKDetectorFSISystematics();
//   #endif
//   #endif
// #endif


  // Don't need to comment TweakAllUncorrelatedXSecSystematics() because there could be 
  // future NIWG uncorrelated xsec systematics

  // NEW FLUX
  
  // common to all the samples
  this->TweakAllUncorrelatedXSecSystematics();  // NEW XSEC NIWG
  this->GetNuXSec()->SyncRHCToFHCSystematics(); // uncorr xsec are fully correlated rhc/fhc

#if  defined(__FHC_RUNNING__) && !defined(__RHC_RUNNING__)
  this->TweakAllBANFFSystematics();
  this->TweakAllSKDetectorFSISystematics();
  this->TweakAllUncorrelatedXSecSystematics();  // NEW XSEC NIWG
#endif
  
#ifdef __RHC_RUNNING__   /// NEW MULTIRING 3
 
  #if defined(__FULLY_CORRELATED__) && defined(__FHC_RUNNING__) // RHC && FHC fully correlated
   this->TweakAllBANFFSystematics();
   this->GetBANFF()->SyncRHCToFHCSystematics();  
   this->TweakAllSKDetectorFSISystematics();
   this->GetSuperkDetectorFSI()->SyncRHCToFHCSystematics();  
 
   // NEW FHC+RHC JOINT 
  #elif defined(__FHC_RUNNING__)   // RHC && FHC with correlations
   this->TweakAllFHCpRHCBANFFSystematics();
   this->TweakAllFHCpRHCSKDetectorFSISystematics();
   //
   
  #else
   this->TweakAllRHCBANFFSystematics(); // RHC only 
   this->TweakAllRHCSKDetectorFSISystematics();
   
  #endif
#endif
  
//
  
  /// NEW MULTIRING 2
#ifdef __MULTIRING_FHC_RUNNING__
  #ifdef __FULLY_CORRELATED_MULTIRING__
  this->GetSuperkDetectorFSI()->SyncMultiRingToFHCSystematics();
  #else
  this->TweakAllMultiRingSKDetectorFSISystematics();
  #endif
#endif
  ///
}








//________________________________________________________________________________________
void SuperkExpectation::TweakAllSampleBANFFSystematics(void)
{
  // NEW FHC+RHC JOINT
#if defined(__FHC_RUNNING__) && defined(__RHC_RUNNING__)
  this->TweakAllFHCpRHCBANFFSystematics();
#elif defined(__FHC_RUNNING__)
  this->TweakAllBANFFSystematics();
#elif defined(__RHC_RUNNING__)
  this->TweakAllRHCBANFFSystematics();
#endif
#ifdef __MULTIRING_RUNNING__
  this->TweakAllMultiRingBANFFSystematics();
#endif
  //

// #ifdef __FHC_RUNNING__
//   this->TweakAllBANFFSystematics();
// #endif
// #ifdef __RHC_RUNNING__
//   this->TweakAllRHCBANFFSystematics();
// #endif
// #ifdef __MULTIRING_RUNNING__
//   this->TweakAllMultiRingBANFFSystematics();
// #endif
}

//________________________________________________________________________________________


// NEW FHC+RHC JOINT
void SuperkExpectation::TweakAllFHCpRHCBANFFSystematics(void)
{
#ifdef __show_syst_tweak_mesg__
  cout << "Exp >> Generating random variations of all FHC+RHC BANFF systematic params (correlated)" << endl;
#endif

  // BANFF
  this->GetBANFF()->TweakFHCpRHC(false,0); // <-- Cholesky method
}
//


//________________________________________________________________________________________
void SuperkExpectation::TweakAllRHCBANFFSystematics(void)
{
#ifdef __show_syst_tweak_mesg__
  cout << "Exp >> Generating random variations of all RHC BANFF systematic params (correlated)" << endl;
#endif

  // BANFF
  this->GetBANFF()->TweakRHC(false,0); // <-- Cholesky method
}
//________________________________________________________________________________________
void SuperkExpectation::TweakAllBANFFSystematics(void)
{
#ifdef __show_syst_tweak_mesg__
    cout << "Exp >> Generating random variations of all BANFF systematic params (correlated)" << endl;
#endif
    // BANFF
    this->GetBANFF()->Tweak(false,0); // <-- Cholesky method
}

//________________________________________________________________________________________
void SuperkExpectation::TweakAllSampleSKDetectorFSISystematics(void)
{
  // NEW FHC+RHC JOINT
#if defined(__FHC_RUNNING__) && defined(__RHC_RUNNING__)
  this->TweakAllFHCpRHCSKDetectorFSISystematics();
#elif defined(__FHC_RUNNING__)
  this->TweakAllSKDetectorFSISystematics();
#elif defined(__RHC_RUNNING__)
  this->TweakAllRHCSKDetectorFSISystematics();
#endif
#ifdef __MULTIRING_RUNNING__
  this->TweakAllMultiRingSKDetectorFSISystematics();
#endif
  //
  
// #ifdef __FHC_RUNNING__
//   this->TweakAllSKDetectorFSISystematics();
// #endif
// #ifdef __RHC_RUNNING__
//   this->TweakAllRHCSKDetectorFSISystematics();
// #endif
// #ifdef __MULTIRING_RUNNING__
//   this->TweakAllMultiRingSKDetectorFSISystematics();
// #endif
}
//________________________________________________________________________________________
void SuperkExpectation::TweakAllSKDetectorFSISystematics(void)
{
#ifdef __show_syst_tweak_mesg__
  cout << "Exp >> Generating random variations of all SKDET + FSI systematic params." << endl;
#endif

  this->GetSuperkDetectorFSI()->TweakEnergyScale(Rnd::Instance()->RndToyMCsyst().Gaus(0.,1.));
  this->GetSuperkDetectorFSI()->Tweak(); // <-- Cholesky method
}
//________________________________________________________________________________________

// NEW FLUX

void SuperkExpectation::TweakAllRHCSKDetectorFSISystematics(void)
{
#ifdef __show_syst_tweak_mesg__
  cout << "Exp >> Generating random variations of all RHC SKDET + FSI systematic params." << endl;
#endif

  this->GetSuperkDetectorFSI()->TweakRHCEnergyScale(Rnd::Instance()->RndToyMCsyst().Gaus(0.,1.));
  this->GetSuperkDetectorFSI()->TweakRHC(); // <-- Cholesky method
}
//


//________________________________________________________________________________________

// NEW FHC+RHC JOINT

void SuperkExpectation::TweakAllFHCpRHCSKDetectorFSISystematics(void)
{
#ifdef __show_syst_tweak_mesg__
  cout << "Exp >> Generating random variations of all FHC+RHC SKDET + FSI systematic params." << endl;
#endif

  this->GetSuperkDetectorFSI()->TweakEnergyScale(Rnd::Instance()->RndToyMCsyst().Gaus(0.,1.));
  this->GetSuperkDetectorFSI()->TweakRHCEnergyScale(Rnd::Instance()->RndToyMCsyst().Gaus(0.,1.));
  this->GetSuperkDetectorFSI()->TweakFHCpRHC(); // <-- Cholesky method
}
//

//________________________________________________________________________________________

/// NEW MULTIRING 2
void SuperkExpectation::TweakAllMultiRingSKDetectorFSISystematics(void)
{
#ifdef __show_syst_tweak_mesg__
  cout << "Exp >> Generating random variations of all MultiRing SKDET + FSI systematic params." << endl;
#endif
  /// NEW MULTIRING 3
  this->GetSuperkDetectorFSI()->TweakMultiRingEnergyScale(Rnd::Instance()->RndToyMCsyst().Gaus(0.,1.));
  ///
  this->GetSuperkDetectorFSI()->TweakMultiRing(); // <-- Cholesky method
}
//________________________________________________________________________________________
void SuperkExpectation::TweakAllUncorrelatedXSecSystematics(void)
{
#ifdef __show_syst_tweak_mesg__
  cout << "Exp >> Generating random variations of all Uncorrelated XSec systematic params." << endl;
#endif

  // NEW XSEC --> comment all old non-banff xsec systematics --> not needed anylonger

  // NEW XSEC NIWG

#if defined(__USE_NIWG2013__)
  this->GetNuXSec()->Tweak(kNuisXSecFermiMom       , Rnd::Instance()->RndToyMCsyst().Gaus(0.,1.));
  double tweak = 2;
  while(tweak > 1 || tweak < 0)
    tweak = Rnd::Instance()->RndToyMCsyst().Gaus(0.,1.);
  this->GetNuXSec()->Tweak(kNuisXSecSpecFunc       , tweak);
  //this->GetNuXSec()->Tweak(kNuisXSecSpecFunc       , Rnd::Instance()->RndToyMCsyst().Uniform(0.,1.));
  //  this->GetNuXSec()->Tweak(kNuisXSecSpecFunc       , Rnd::Instance()->RndToyMCsyst().Gaus(0.,1.));
  this->GetNuXSec()->Tweak(kNuisXSecWshape         , Rnd::Instance()->RndToyMCsyst().Gaus(0.,1.));
  this->GetNuXSec()->Tweak(kNuisXSecShapeCCoth     , Rnd::Instance()->RndToyMCsyst().Gaus(0.,1.));
  this->GetNuXSec()->Tweak(kNuisXSecPiLessDelta    , Rnd::Instance()->RndToyMCsyst().Gaus(0.,1.));
  #ifndef __NUEONLY__
  this->GetNuXSec()->Tweak(kNuisXSecBindE          , Rnd::Instance()->RndToyMCsyst().Gaus(0.,1.));
  #endif
  this->GetNuXSec()->Tweak(kNuisXSecNormCCcoh      , Rnd::Instance()->RndToyMCsyst().Gaus(0.,1.));
  #ifndef __NUEONLY__
  this->GetNuXSec()->Tweak(kNuisXSecNormNC1piC     , Rnd::Instance()->RndToyMCsyst().Gaus(0.,1.));
  #endif
  #ifndef __NUMUONLY__
  this->GetNuXSec()->Tweak(kNuisXSecNormNCcoh      , Rnd::Instance()->RndToyMCsyst().Gaus(0.,1.));
  #endif
  this->GetNuXSec()->Tweak(kNuisXSecNormNCoth      , Rnd::Instance()->RndToyMCsyst().Gaus(0.,1.));
  this->GetNuXSec()->Tweak(kNuisXSecNormNuEToNuMu  , Rnd::Instance()->RndToyMCsyst().Gaus(0.,1.));
  #ifndef __NUEONLY__
  this->GetNuXSec()->Tweak(kNuisXSecNormAntiNuToNu , Rnd::Instance()->RndToyMCsyst().Gaus(0.,1.));
  #endif

#elif defined(__USE_NIWG2015__)
  // NEW XSEC
  // No longer tweak RPA
  //~ double tweak_discr = Rnd::Instance()->RndToyMCsyst().Uniform(0.,1.);
  //~ if(tweak_discr < 0.5) tweak_discr = 0.;
  //~ else tweak_discr = 1.;
  //~ this->GetNuXSec()->Tweak(kNuisXSecRPA             , tweak_discr);
#endif

}
//________________________________________________________________________________________
PhysParamList SuperkExpectation::Tweaks(void) 
{
  PhysParamList pl;
  this->CopyTweaksTo(pl);
  return pl;
}
//________________________________________________________________________________________
void SuperkExpectation::CopyTweaksTo(PhysParamList & pl) 
{
  pl.Clear();
  //
  // BANFF
  //
#ifdef __FHC_RUNNING__  
  const int nbanff = conv::kNumOfParamsInBANFFfit;
  for(int ip = 0; ip < nbanff; ip++) {
    PhysParamId_t p = kListOfBANFFParams[ip];
    double v = this->GetBANFF()->CurrentValue(p);
    pl.SetValue(p,v);
  }
#endif
  
  //
  // RHC BANFF
  //
#ifdef __RHC_RUNNING__   
  const int nbanff_rhc = conv::kNumOfRHCParamsInBANFFfit;
  for(int ip = 0; ip < nbanff_rhc; ip++) {
     PhysParamId_t p = kListOfRHCBANFFParams[ip];
     double v = this->GetBANFF()->RHCCurrentValue(p);
     pl.SetValue(p,v);
  }
#endif

  // 
  //  Uncorrelated xsec (common to FHC and RHC)
  // 

  
#if defined(__USE_NIWG2013__)
  pl.SetValue( kNuisXSecFermiMom,        this->GetNuXSec()->CurrentValue(kNuisXSecFermiMom)       );
  pl.SetValue( kNuisXSecSpecFunc,        this->GetNuXSec()->CurrentValue(kNuisXSecSpecFunc)       );
  pl.SetValue( kNuisXSecWshape,          this->GetNuXSec()->CurrentValue(kNuisXSecWshape)         );
  pl.SetValue( kNuisXSecShapeCCoth,      this->GetNuXSec()->CurrentValue(kNuisXSecShapeCCoth)     );
  pl.SetValue( kNuisXSecPiLessDelta,     this->GetNuXSec()->CurrentValue(kNuisXSecPiLessDelta)    );
  #ifndef __NUEONLY__
  pl.SetValue( kNuisXSecBindE,           this->GetNuXSec()->CurrentValue(kNuisXSecBindE)          );
  #endif
  pl.SetValue( kNuisXSecNormCCcoh,       this->GetNuXSec()->CurrentValue(kNuisXSecNormCCcoh)      );
  #ifndef __NUEONLY__
  pl.SetValue( kNuisXSecNormNC1piC,      this->GetNuXSec()->CurrentValue(kNuisXSecNormNC1piC)     );
  #endif
  #ifndef __NUMUONLY__
  pl.SetValue( kNuisXSecNormNCcoh,       this->GetNuXSec()->CurrentValue(kNuisXSecNormNCcoh)      );
  #endif
  pl.SetValue( kNuisXSecNormNCoth,       this->GetNuXSec()->CurrentValue(kNuisXSecNormNCoth)      );
  pl.SetValue( kNuisXSecNormNuEToNuMu,   this->GetNuXSec()->CurrentValue(kNuisXSecNormNuEToNuMu)  );
  #ifndef __NUEONLY__
  pl.SetValue( kNuisXSecNormAntiNuToNu,  this->GetNuXSec()->CurrentValue(kNuisXSecNormAntiNuToNu) );
  #endif

#elif defined(__USE_NIWG2015__)
  pl.SetValue( kNuisXSecRPA,             this->GetNuXSec()->CurrentValue(kNuisXSecRPA)            ); 
#endif




  // 
  // SuperK detector + FSI systematics
  //
  
  // FHC
#ifdef __FHC_RUNNING__  
  pl.SetValue( kNuisSuperkEnergyScale,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSuperkEnergyScale) );
#if FSI_BINS != 3
  pl.SetValue( kNuisSKDetFSI000,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSKDetFSI000) );
  pl.SetValue( kNuisSKDetFSI001,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSKDetFSI001) );
  pl.SetValue( kNuisSKDetFSI002,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSKDetFSI002) );
  pl.SetValue( kNuisSKDetFSI003,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSKDetFSI003) );
  pl.SetValue( kNuisSKDetFSI004,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSKDetFSI004) );
  pl.SetValue( kNuisSKDetFSI005,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSKDetFSI005) );
  
  #if defined(__USE_SKDETFSI2015__)
  pl.SetValue( kNuisSKDetFSI006,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSKDetFSI006) );
  pl.SetValue( kNuisSKDetFSI007,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSKDetFSI007) );
  pl.SetValue( kNuisSKDetFSI008,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSKDetFSI008) );
  pl.SetValue( kNuisSKDetFSI009,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSKDetFSI009) );
  pl.SetValue( kNuisSKDetFSI010,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSKDetFSI010) );
  pl.SetValue( kNuisSKDetFSI011,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSKDetFSI011) );
  #endif

#endif
#if FSI_BINS != 2
  #if defined(__USE_SKDETFSI2013__)
  pl.SetValue( kNuisSKDetFSI006,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSKDetFSI006) );
  pl.SetValue( kNuisSKDetFSI007,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSKDetFSI007) );
  pl.SetValue( kNuisSKDetFSI008,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSKDetFSI008) );
  pl.SetValue( kNuisSKDetFSI009,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSKDetFSI009) );
  pl.SetValue( kNuisSKDetFSI010,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSKDetFSI010) );
  pl.SetValue( kNuisSKDetFSI011,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSKDetFSI011) );
  pl.SetValue( kNuisSKDetFSI012,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSKDetFSI012) );
  pl.SetValue( kNuisSKDetFSI013,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSKDetFSI013) );
  pl.SetValue( kNuisSKDetFSI014,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSKDetFSI014) );
  pl.SetValue( kNuisSKDetFSI015,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSKDetFSI015) );
  pl.SetValue( kNuisSKDetFSI016,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSKDetFSI016) );
  pl.SetValue( kNuisSKDetFSI017,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSKDetFSI017) );
  #elif defined(__USE_SKDETFSI2015__)
  pl.SetValue( kNuisSKDetFSI012,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSKDetFSI012) );
  pl.SetValue( kNuisSKDetFSI013,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSKDetFSI013) );
  pl.SetValue( kNuisSKDetFSI014,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSKDetFSI014) );
  pl.SetValue( kNuisSKDetFSI015,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSKDetFSI015) );
  pl.SetValue( kNuisSKDetFSI016,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSKDetFSI016) );
  pl.SetValue( kNuisSKDetFSI017,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSKDetFSI017) );
  pl.SetValue( kNuisSKDetFSI018,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSKDetFSI018) );
  pl.SetValue( kNuisSKDetFSI019,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSKDetFSI019) );
  pl.SetValue( kNuisSKDetFSI020,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSKDetFSI020) );
  pl.SetValue( kNuisSKDetFSI021,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSKDetFSI021) );
  pl.SetValue( kNuisSKDetFSI022,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSKDetFSI022) );
  pl.SetValue( kNuisSKDetFSI023,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSKDetFSI023) );
  #endif

#ifdef SKCORRELATIONS

  #if defined(__USE_SKDETFSI2013__)
  pl.SetValue( kNuisSKDetFSI018,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSKDetFSI018) );
  pl.SetValue( kNuisSKDetFSI019,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSKDetFSI019) );
  pl.SetValue( kNuisSKDetFSI020,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSKDetFSI020) );
  pl.SetValue( kNuisSKDetFSI021,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSKDetFSI021) );
  pl.SetValue( kNuisSKDetFSI022,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSKDetFSI022) );
  pl.SetValue( kNuisSKDetFSI023,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSKDetFSI023) );
  pl.SetValue( kNuisSKDetFSI024,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSKDetFSI024) );
  pl.SetValue( kNuisSKDetFSI025,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSKDetFSI025) );
  pl.SetValue( kNuisSKDetFSI026,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSKDetFSI026) );
  pl.SetValue( kNuisSKDetFSI027,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSKDetFSI027) );
  pl.SetValue( kNuisSKDetFSI028,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSKDetFSI028) );
  pl.SetValue( kNuisSKDetFSI029,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSKDetFSI029) );
  pl.SetValue( kNuisSKDetFSI030,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSKDetFSI030) );
  pl.SetValue( kNuisSKDetFSI031,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSKDetFSI031) );
  #elif defined(__USE_SKDETFSI2015__)
  pl.SetValue( kNuisSKDetFSI024,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSKDetFSI024) );
  pl.SetValue( kNuisSKDetFSI025,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSKDetFSI025) );
  pl.SetValue( kNuisSKDetFSI026,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSKDetFSI026) );
  pl.SetValue( kNuisSKDetFSI027,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSKDetFSI027) );
  pl.SetValue( kNuisSKDetFSI028,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSKDetFSI028) );
  pl.SetValue( kNuisSKDetFSI029,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSKDetFSI029) );
  pl.SetValue( kNuisSKDetFSI030,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSKDetFSI030) );
  pl.SetValue( kNuisSKDetFSI031,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSKDetFSI031) );
  pl.SetValue( kNuisSKDetFSI032,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSKDetFSI032) );
  pl.SetValue( kNuisSKDetFSI033,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSKDetFSI033) );
  pl.SetValue( kNuisSKDetFSI034,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSKDetFSI034) );
  pl.SetValue( kNuisSKDetFSI035,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSKDetFSI035) );
  pl.SetValue( kNuisSKDetFSI036,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSKDetFSI036) );
  pl.SetValue( kNuisSKDetFSI037,  this->GetSuperkDetectorFSI()->CurrentValue(kNuisSKDetFSI037) );
  #endif

#endif
#endif // FSI_BINS=1
#endif // FHC


  // RHC
#ifdef __RHC_RUNNING__  
  pl.SetValue( kRHCNuisSuperkEnergyScale,  this->GetSuperkDetectorFSI()->RHCCurrentValue(kRHCNuisSuperkEnergyScale) );
#if FSI_BINS != 3
  pl.SetValue( kRHCNuisSKDetFSI000,  this->GetSuperkDetectorFSI()->RHCCurrentValue(kRHCNuisSKDetFSI000) );
  pl.SetValue( kRHCNuisSKDetFSI001,  this->GetSuperkDetectorFSI()->RHCCurrentValue(kRHCNuisSKDetFSI001) );
  pl.SetValue( kRHCNuisSKDetFSI002,  this->GetSuperkDetectorFSI()->RHCCurrentValue(kRHCNuisSKDetFSI002) );
  pl.SetValue( kRHCNuisSKDetFSI003,  this->GetSuperkDetectorFSI()->RHCCurrentValue(kRHCNuisSKDetFSI003) );
  pl.SetValue( kRHCNuisSKDetFSI004,  this->GetSuperkDetectorFSI()->RHCCurrentValue(kRHCNuisSKDetFSI004) );
  pl.SetValue( kRHCNuisSKDetFSI005,  this->GetSuperkDetectorFSI()->RHCCurrentValue(kRHCNuisSKDetFSI005) );
  
  #if defined(__USE_SKDETFSI2015__)
  pl.SetValue( kRHCNuisSKDetFSI006,  this->GetSuperkDetectorFSI()->RHCCurrentValue(kRHCNuisSKDetFSI006) );
  pl.SetValue( kRHCNuisSKDetFSI007,  this->GetSuperkDetectorFSI()->RHCCurrentValue(kRHCNuisSKDetFSI007) );
  pl.SetValue( kRHCNuisSKDetFSI008,  this->GetSuperkDetectorFSI()->RHCCurrentValue(kRHCNuisSKDetFSI008) );
  pl.SetValue( kRHCNuisSKDetFSI009,  this->GetSuperkDetectorFSI()->RHCCurrentValue(kRHCNuisSKDetFSI009) );
  pl.SetValue( kRHCNuisSKDetFSI010,  this->GetSuperkDetectorFSI()->RHCCurrentValue(kRHCNuisSKDetFSI010) );
  pl.SetValue( kRHCNuisSKDetFSI011,  this->GetSuperkDetectorFSI()->RHCCurrentValue(kRHCNuisSKDetFSI011) );
  #endif

#endif
#if FSI_BINS != 2
  #if defined(__USE_SKDETFSI2013__)
  pl.SetValue( kRHCNuisSKDetFSI006,  this->GetSuperkDetectorFSI()->RHCCurrentValue(kRHCNuisSKDetFSI006) );
  pl.SetValue( kRHCNuisSKDetFSI007,  this->GetSuperkDetectorFSI()->RHCCurrentValue(kRHCNuisSKDetFSI007) );
  pl.SetValue( kRHCNuisSKDetFSI008,  this->GetSuperkDetectorFSI()->RHCCurrentValue(kRHCNuisSKDetFSI008) );
  pl.SetValue( kRHCNuisSKDetFSI009,  this->GetSuperkDetectorFSI()->RHCCurrentValue(kRHCNuisSKDetFSI009) );
  pl.SetValue( kRHCNuisSKDetFSI010,  this->GetSuperkDetectorFSI()->RHCCurrentValue(kRHCNuisSKDetFSI010) );
  pl.SetValue( kRHCNuisSKDetFSI011,  this->GetSuperkDetectorFSI()->RHCCurrentValue(kRHCNuisSKDetFSI011) );
  pl.SetValue( kRHCNuisSKDetFSI012,  this->GetSuperkDetectorFSI()->RHCCurrentValue(kRHCNuisSKDetFSI012) );
  pl.SetValue( kRHCNuisSKDetFSI013,  this->GetSuperkDetectorFSI()->RHCCurrentValue(kRHCNuisSKDetFSI013) );
  pl.SetValue( kRHCNuisSKDetFSI014,  this->GetSuperkDetectorFSI()->RHCCurrentValue(kRHCNuisSKDetFSI014) );
  pl.SetValue( kRHCNuisSKDetFSI015,  this->GetSuperkDetectorFSI()->RHCCurrentValue(kRHCNuisSKDetFSI015) );
  pl.SetValue( kRHCNuisSKDetFSI016,  this->GetSuperkDetectorFSI()->RHCCurrentValue(kRHCNuisSKDetFSI016) );
  pl.SetValue( kRHCNuisSKDetFSI017,  this->GetSuperkDetectorFSI()->RHCCurrentValue(kRHCNuisSKDetFSI017) );
  #elif defined(__USE_SKDETFSI2015__)
  pl.SetValue( kRHCNuisSKDetFSI012,  this->GetSuperkDetectorFSI()->RHCCurrentValue(kRHCNuisSKDetFSI012) );
  pl.SetValue( kRHCNuisSKDetFSI013,  this->GetSuperkDetectorFSI()->RHCCurrentValue(kRHCNuisSKDetFSI013) );
  pl.SetValue( kRHCNuisSKDetFSI014,  this->GetSuperkDetectorFSI()->RHCCurrentValue(kRHCNuisSKDetFSI014) );
  pl.SetValue( kRHCNuisSKDetFSI015,  this->GetSuperkDetectorFSI()->RHCCurrentValue(kRHCNuisSKDetFSI015) );
  pl.SetValue( kRHCNuisSKDetFSI016,  this->GetSuperkDetectorFSI()->RHCCurrentValue(kRHCNuisSKDetFSI016) );
  pl.SetValue( kRHCNuisSKDetFSI017,  this->GetSuperkDetectorFSI()->RHCCurrentValue(kRHCNuisSKDetFSI017) );
  pl.SetValue( kRHCNuisSKDetFSI018,  this->GetSuperkDetectorFSI()->RHCCurrentValue(kRHCNuisSKDetFSI018) );
  pl.SetValue( kRHCNuisSKDetFSI019,  this->GetSuperkDetectorFSI()->RHCCurrentValue(kRHCNuisSKDetFSI019) );
  pl.SetValue( kRHCNuisSKDetFSI020,  this->GetSuperkDetectorFSI()->RHCCurrentValue(kRHCNuisSKDetFSI020) );
  pl.SetValue( kRHCNuisSKDetFSI021,  this->GetSuperkDetectorFSI()->RHCCurrentValue(kRHCNuisSKDetFSI021) );
  pl.SetValue( kRHCNuisSKDetFSI022,  this->GetSuperkDetectorFSI()->RHCCurrentValue(kRHCNuisSKDetFSI022) );
  pl.SetValue( kRHCNuisSKDetFSI023,  this->GetSuperkDetectorFSI()->RHCCurrentValue(kRHCNuisSKDetFSI023) );
  #endif

#ifdef SKCORRELATIONS

  #if defined(__USE_SKDETFSI2013__)
  pl.SetValue( kRHCNuisSKDetFSI018,  this->GetSuperkDetectorFSI()->RHCCurrentValue(kRHCNuisSKDetFSI018) );
  pl.SetValue( kRHCNuisSKDetFSI019,  this->GetSuperkDetectorFSI()->RHCCurrentValue(kRHCNuisSKDetFSI019) );
  pl.SetValue( kRHCNuisSKDetFSI020,  this->GetSuperkDetectorFSI()->RHCCurrentValue(kRHCNuisSKDetFSI020) );
  pl.SetValue( kRHCNuisSKDetFSI021,  this->GetSuperkDetectorFSI()->RHCCurrentValue(kRHCNuisSKDetFSI021) );
  pl.SetValue( kRHCNuisSKDetFSI022,  this->GetSuperkDetectorFSI()->RHCCurrentValue(kRHCNuisSKDetFSI022) );
  pl.SetValue( kRHCNuisSKDetFSI023,  this->GetSuperkDetectorFSI()->RHCCurrentValue(kRHCNuisSKDetFSI023) );
  pl.SetValue( kRHCNuisSKDetFSI024,  this->GetSuperkDetectorFSI()->RHCCurrentValue(kRHCNuisSKDetFSI024) );
  pl.SetValue( kRHCNuisSKDetFSI025,  this->GetSuperkDetectorFSI()->RHCCurrentValue(kRHCNuisSKDetFSI025) );
  pl.SetValue( kRHCNuisSKDetFSI026,  this->GetSuperkDetectorFSI()->RHCCurrentValue(kRHCNuisSKDetFSI026) );
  pl.SetValue( kRHCNuisSKDetFSI027,  this->GetSuperkDetectorFSI()->RHCCurrentValue(kRHCNuisSKDetFSI027) );
  pl.SetValue( kRHCNuisSKDetFSI028,  this->GetSuperkDetectorFSI()->RHCCurrentValue(kRHCNuisSKDetFSI028) );
  pl.SetValue( kRHCNuisSKDetFSI029,  this->GetSuperkDetectorFSI()->RHCCurrentValue(kRHCNuisSKDetFSI029) );
  pl.SetValue( kRHCNuisSKDetFSI030,  this->GetSuperkDetectorFSI()->RHCCurrentValue(kRHCNuisSKDetFSI030) );
  pl.SetValue( kRHCNuisSKDetFSI031,  this->GetSuperkDetectorFSI()->RHCCurrentValue(kRHCNuisSKDetFSI031) );
  #elif defined(__USE_SKDETFSI2015__)
  pl.SetValue( kRHCNuisSKDetFSI024,  this->GetSuperkDetectorFSI()->RHCCurrentValue(kRHCNuisSKDetFSI024) );
  pl.SetValue( kRHCNuisSKDetFSI025,  this->GetSuperkDetectorFSI()->RHCCurrentValue(kRHCNuisSKDetFSI025) );
  pl.SetValue( kRHCNuisSKDetFSI026,  this->GetSuperkDetectorFSI()->RHCCurrentValue(kRHCNuisSKDetFSI026) );
  pl.SetValue( kRHCNuisSKDetFSI027,  this->GetSuperkDetectorFSI()->RHCCurrentValue(kRHCNuisSKDetFSI027) );
  pl.SetValue( kRHCNuisSKDetFSI028,  this->GetSuperkDetectorFSI()->RHCCurrentValue(kRHCNuisSKDetFSI028) );
  pl.SetValue( kRHCNuisSKDetFSI029,  this->GetSuperkDetectorFSI()->RHCCurrentValue(kRHCNuisSKDetFSI029) );
  pl.SetValue( kRHCNuisSKDetFSI030,  this->GetSuperkDetectorFSI()->RHCCurrentValue(kRHCNuisSKDetFSI030) );
  pl.SetValue( kRHCNuisSKDetFSI031,  this->GetSuperkDetectorFSI()->RHCCurrentValue(kRHCNuisSKDetFSI031) );
  pl.SetValue( kRHCNuisSKDetFSI032,  this->GetSuperkDetectorFSI()->RHCCurrentValue(kRHCNuisSKDetFSI032) );
  pl.SetValue( kRHCNuisSKDetFSI033,  this->GetSuperkDetectorFSI()->RHCCurrentValue(kRHCNuisSKDetFSI033) );
  pl.SetValue( kRHCNuisSKDetFSI034,  this->GetSuperkDetectorFSI()->RHCCurrentValue(kRHCNuisSKDetFSI034) );
  pl.SetValue( kRHCNuisSKDetFSI035,  this->GetSuperkDetectorFSI()->RHCCurrentValue(kRHCNuisSKDetFSI035) );
  pl.SetValue( kRHCNuisSKDetFSI036,  this->GetSuperkDetectorFSI()->RHCCurrentValue(kRHCNuisSKDetFSI036) );
  pl.SetValue( kRHCNuisSKDetFSI037,  this->GetSuperkDetectorFSI()->RHCCurrentValue(kRHCNuisSKDetFSI037) );
  #endif

#endif
#endif // FSI_BINS=1
#endif // RHC






#ifdef __MULTIRING_FHC_RUNNING__

  pl.SetValue( kMultiRingNuisSuperkEnergyScale,  this->GetSuperkDetectorFSI()->MultiRingCurrentValue(kMultiRingNuisSuperkEnergyScale) ); 
  ///
#if FSI_BINS != 3
  pl.SetValue( kMultiRingNuisSKDetFSI000,  this->GetSuperkDetectorFSI()->MultiRingCurrentValue(kMultiRingNuisSKDetFSI000) );
  pl.SetValue( kMultiRingNuisSKDetFSI001,  this->GetSuperkDetectorFSI()->MultiRingCurrentValue(kMultiRingNuisSKDetFSI001) );
  pl.SetValue( kMultiRingNuisSKDetFSI002,  this->GetSuperkDetectorFSI()->MultiRingCurrentValue(kMultiRingNuisSKDetFSI002) );
  pl.SetValue( kMultiRingNuisSKDetFSI003,  this->GetSuperkDetectorFSI()->MultiRingCurrentValue(kMultiRingNuisSKDetFSI003) );
  pl.SetValue( kMultiRingNuisSKDetFSI004,  this->GetSuperkDetectorFSI()->MultiRingCurrentValue(kMultiRingNuisSKDetFSI004) );
  pl.SetValue( kMultiRingNuisSKDetFSI005,  this->GetSuperkDetectorFSI()->MultiRingCurrentValue(kMultiRingNuisSKDetFSI005) );

  // NEW SKDETFSI
  #if defined(__USE_SKDETFSI2015__)
  pl.SetValue( kMultiRingNuisSKDetFSI006,  this->GetSuperkDetectorFSI()->MultiRingCurrentValue(kMultiRingNuisSKDetFSI006) );
  pl.SetValue( kMultiRingNuisSKDetFSI007,  this->GetSuperkDetectorFSI()->MultiRingCurrentValue(kMultiRingNuisSKDetFSI007) );
  pl.SetValue( kMultiRingNuisSKDetFSI008,  this->GetSuperkDetectorFSI()->MultiRingCurrentValue(kMultiRingNuisSKDetFSI008) );
  pl.SetValue( kMultiRingNuisSKDetFSI009,  this->GetSuperkDetectorFSI()->MultiRingCurrentValue(kMultiRingNuisSKDetFSI009) );
  pl.SetValue( kMultiRingNuisSKDetFSI010,  this->GetSuperkDetectorFSI()->MultiRingCurrentValue(kMultiRingNuisSKDetFSI010) );
  pl.SetValue( kMultiRingNuisSKDetFSI011,  this->GetSuperkDetectorFSI()->MultiRingCurrentValue(kMultiRingNuisSKDetFSI011) );
  #endif

#endif
#if FSI_BINS != 2

  // NEW SKDETFSI
  #if defined(__USE_SKDETFSI2013__)
  pl.SetValue( kMultiRingNuisSKDetFSI006,  this->GetSuperkDetectorFSI()->MultiRingCurrentValue(kMultiRingNuisSKDetFSI006) );
  pl.SetValue( kMultiRingNuisSKDetFSI007,  this->GetSuperkDetectorFSI()->MultiRingCurrentValue(kMultiRingNuisSKDetFSI007) );
  pl.SetValue( kMultiRingNuisSKDetFSI008,  this->GetSuperkDetectorFSI()->MultiRingCurrentValue(kMultiRingNuisSKDetFSI008) );
  pl.SetValue( kMultiRingNuisSKDetFSI009,  this->GetSuperkDetectorFSI()->MultiRingCurrentValue(kMultiRingNuisSKDetFSI009) );
  pl.SetValue( kMultiRingNuisSKDetFSI010,  this->GetSuperkDetectorFSI()->MultiRingCurrentValue(kMultiRingNuisSKDetFSI010) );
  pl.SetValue( kMultiRingNuisSKDetFSI011,  this->GetSuperkDetectorFSI()->MultiRingCurrentValue(kMultiRingNuisSKDetFSI011) );
  pl.SetValue( kMultiRingNuisSKDetFSI012,  this->GetSuperkDetectorFSI()->MultiRingCurrentValue(kMultiRingNuisSKDetFSI012) );
  pl.SetValue( kMultiRingNuisSKDetFSI013,  this->GetSuperkDetectorFSI()->MultiRingCurrentValue(kMultiRingNuisSKDetFSI013) );
  pl.SetValue( kMultiRingNuisSKDetFSI014,  this->GetSuperkDetectorFSI()->MultiRingCurrentValue(kMultiRingNuisSKDetFSI014) );
  pl.SetValue( kMultiRingNuisSKDetFSI015,  this->GetSuperkDetectorFSI()->MultiRingCurrentValue(kMultiRingNuisSKDetFSI015) );
  pl.SetValue( kMultiRingNuisSKDetFSI016,  this->GetSuperkDetectorFSI()->MultiRingCurrentValue(kMultiRingNuisSKDetFSI016) );
  pl.SetValue( kMultiRingNuisSKDetFSI017,  this->GetSuperkDetectorFSI()->MultiRingCurrentValue(kMultiRingNuisSKDetFSI017) );
  #elif defined(__USE_SKDETFSI2015__)
  pl.SetValue( kMultiRingNuisSKDetFSI012,  this->GetSuperkDetectorFSI()->MultiRingCurrentValue(kMultiRingNuisSKDetFSI012) );
  pl.SetValue( kMultiRingNuisSKDetFSI013,  this->GetSuperkDetectorFSI()->MultiRingCurrentValue(kMultiRingNuisSKDetFSI013) );
  pl.SetValue( kMultiRingNuisSKDetFSI014,  this->GetSuperkDetectorFSI()->MultiRingCurrentValue(kMultiRingNuisSKDetFSI014) );
  pl.SetValue( kMultiRingNuisSKDetFSI015,  this->GetSuperkDetectorFSI()->MultiRingCurrentValue(kMultiRingNuisSKDetFSI015) );
  pl.SetValue( kMultiRingNuisSKDetFSI016,  this->GetSuperkDetectorFSI()->MultiRingCurrentValue(kMultiRingNuisSKDetFSI016) );
  pl.SetValue( kMultiRingNuisSKDetFSI017,  this->GetSuperkDetectorFSI()->MultiRingCurrentValue(kMultiRingNuisSKDetFSI017) );
  pl.SetValue( kMultiRingNuisSKDetFSI018,  this->GetSuperkDetectorFSI()->MultiRingCurrentValue(kMultiRingNuisSKDetFSI018) );
  pl.SetValue( kMultiRingNuisSKDetFSI019,  this->GetSuperkDetectorFSI()->MultiRingCurrentValue(kMultiRingNuisSKDetFSI019) );
  pl.SetValue( kMultiRingNuisSKDetFSI020,  this->GetSuperkDetectorFSI()->MultiRingCurrentValue(kMultiRingNuisSKDetFSI020) );
  pl.SetValue( kMultiRingNuisSKDetFSI021,  this->GetSuperkDetectorFSI()->MultiRingCurrentValue(kMultiRingNuisSKDetFSI021) );
  pl.SetValue( kMultiRingNuisSKDetFSI022,  this->GetSuperkDetectorFSI()->MultiRingCurrentValue(kMultiRingNuisSKDetFSI022) );
  pl.SetValue( kMultiRingNuisSKDetFSI023,  this->GetSuperkDetectorFSI()->MultiRingCurrentValue(kMultiRingNuisSKDetFSI023) );
  #endif

#ifdef SKCORRELATIONS
  // NEW SKDETFSI
  #if defined(__USE_SKDETFSI2013__)
  pl.SetValue( kMultiRingNuisSKDetFSI018,  this->GetSuperkDetectorFSI()->MultiRingCurrentValue(kMultiRingNuisSKDetFSI018) );
  pl.SetValue( kMultiRingNuisSKDetFSI019,  this->GetSuperkDetectorFSI()->MultiRingCurrentValue(kMultiRingNuisSKDetFSI019) );
  pl.SetValue( kMultiRingNuisSKDetFSI020,  this->GetSuperkDetectorFSI()->MultiRingCurrentValue(kMultiRingNuisSKDetFSI020) );
  pl.SetValue( kMultiRingNuisSKDetFSI021,  this->GetSuperkDetectorFSI()->MultiRingCurrentValue(kMultiRingNuisSKDetFSI021) );
  pl.SetValue( kMultiRingNuisSKDetFSI022,  this->GetSuperkDetectorFSI()->MultiRingCurrentValue(kMultiRingNuisSKDetFSI022) );
  pl.SetValue( kMultiRingNuisSKDetFSI023,  this->GetSuperkDetectorFSI()->MultiRingCurrentValue(kMultiRingNuisSKDetFSI023) );
  pl.SetValue( kMultiRingNuisSKDetFSI024,  this->GetSuperkDetectorFSI()->MultiRingCurrentValue(kMultiRingNuisSKDetFSI024) );
  pl.SetValue( kMultiRingNuisSKDetFSI025,  this->GetSuperkDetectorFSI()->MultiRingCurrentValue(kMultiRingNuisSKDetFSI025) );
  pl.SetValue( kMultiRingNuisSKDetFSI026,  this->GetSuperkDetectorFSI()->MultiRingCurrentValue(kMultiRingNuisSKDetFSI026) );
  pl.SetValue( kMultiRingNuisSKDetFSI027,  this->GetSuperkDetectorFSI()->MultiRingCurrentValue(kMultiRingNuisSKDetFSI027) );
  pl.SetValue( kMultiRingNuisSKDetFSI028,  this->GetSuperkDetectorFSI()->MultiRingCurrentValue(kMultiRingNuisSKDetFSI028) );
  pl.SetValue( kMultiRingNuisSKDetFSI029,  this->GetSuperkDetectorFSI()->MultiRingCurrentValue(kMultiRingNuisSKDetFSI029) );
  pl.SetValue( kMultiRingNuisSKDetFSI030,  this->GetSuperkDetectorFSI()->MultiRingCurrentValue(kMultiRingNuisSKDetFSI030) );
  pl.SetValue( kMultiRingNuisSKDetFSI031,  this->GetSuperkDetectorFSI()->MultiRingCurrentValue(kMultiRingNuisSKDetFSI031) );
  #elif defined(__USE_SKDETFSI2015__)
  pl.SetValue( kMultiRingNuisSKDetFSI024,  this->GetSuperkDetectorFSI()->MultiRingCurrentValue(kMultiRingNuisSKDetFSI024) );
  pl.SetValue( kMultiRingNuisSKDetFSI025,  this->GetSuperkDetectorFSI()->MultiRingCurrentValue(kMultiRingNuisSKDetFSI025) );
  pl.SetValue( kMultiRingNuisSKDetFSI026,  this->GetSuperkDetectorFSI()->MultiRingCurrentValue(kMultiRingNuisSKDetFSI026) );
  pl.SetValue( kMultiRingNuisSKDetFSI027,  this->GetSuperkDetectorFSI()->MultiRingCurrentValue(kMultiRingNuisSKDetFSI027) );
  pl.SetValue( kMultiRingNuisSKDetFSI028,  this->GetSuperkDetectorFSI()->MultiRingCurrentValue(kMultiRingNuisSKDetFSI028) );
  pl.SetValue( kMultiRingNuisSKDetFSI029,  this->GetSuperkDetectorFSI()->MultiRingCurrentValue(kMultiRingNuisSKDetFSI029) );
  pl.SetValue( kMultiRingNuisSKDetFSI030,  this->GetSuperkDetectorFSI()->MultiRingCurrentValue(kMultiRingNuisSKDetFSI030) );
  pl.SetValue( kMultiRingNuisSKDetFSI031,  this->GetSuperkDetectorFSI()->MultiRingCurrentValue(kMultiRingNuisSKDetFSI031) );
  pl.SetValue( kMultiRingNuisSKDetFSI032,  this->GetSuperkDetectorFSI()->MultiRingCurrentValue(kMultiRingNuisSKDetFSI032) );
  pl.SetValue( kMultiRingNuisSKDetFSI033,  this->GetSuperkDetectorFSI()->MultiRingCurrentValue(kMultiRingNuisSKDetFSI033) );
  pl.SetValue( kMultiRingNuisSKDetFSI034,  this->GetSuperkDetectorFSI()->MultiRingCurrentValue(kMultiRingNuisSKDetFSI034) );
  pl.SetValue( kMultiRingNuisSKDetFSI035,  this->GetSuperkDetectorFSI()->MultiRingCurrentValue(kMultiRingNuisSKDetFSI035) );
  pl.SetValue( kMultiRingNuisSKDetFSI036,  this->GetSuperkDetectorFSI()->MultiRingCurrentValue(kMultiRingNuisSKDetFSI036) );
  pl.SetValue( kMultiRingNuisSKDetFSI037,  this->GetSuperkDetectorFSI()->MultiRingCurrentValue(kMultiRingNuisSKDetFSI037) );
  #endif

#endif
#endif // FSI_BINS=1 
#endif // MULTIRING FHC

  return;
}


//________________________________________________________________________________________
void SuperkExpectation::SetNuisanceParams(const FitParamList & fp, bool /*all_banff_fit*/, bool /*all_1Re_eff*/)
{
  //this->ResetSystematicVariations();
  
  //below is what it used to do - check that every single systematic parameter exists in the fit,
  // and only then tweak
  //FitParamList has a 'SystMode', which tells you which parameters are being fit - use this instead!
  if(fp.GetSystMode() == 1 || fp.GetSystMode() == 3) {


    /////////////// ENERGY SCALE //////////////

    // if(fp.Exists(kNuisSuperkEnergyScale)) {
    //   double twk = fp.CurrValue(kNuisSuperkEnergyScale);
    //   //Todo: fix this so that FHC and RHC can have different tweaks
    //   this->GetSuperkDetectorFSI()->TweakEnergyScale(twk);
    //   this->GetSuperkDetectorFSI()->TweakRHCEnergyScale(twk);
    // }

#if defined(__FHC_RUNNING__) && defined(__RHC_RUNNING__)
  #if defined(__FULLY_CORRELATED__) // same as before
    if(fp.Exists(kNuisSuperkEnergyScale)) {
      double twk = fp.CurrValue(kNuisSuperkEnergyScale);
      this->GetSuperkDetectorFSI()->TweakEnergyScale(twk);
      this->GetSuperkDetectorFSI()->TweakRHCEnergyScale(twk);
    }
  #else

    // NEW FHC+RHC JOINT
    if(fp.Exists(kNuisSuperkEnergyScale)) {
      double twk = fp.CurrValue(kNuisSuperkEnergyScale);
      this->GetSuperkDetectorFSI()->TweakEnergyScale(twk);
    }
    if(fp.Exists(kRHCNuisSuperkEnergyScale)) {
      double twk = fp.CurrValue(kRHCNuisSuperkEnergyScale);
      this->GetSuperkDetectorFSI()->TweakRHCEnergyScale(twk);
    }
    //
    
  #endif
#else 
    // fhc only
    if(fp.Exists(kNuisSuperkEnergyScale)) {
      double twk = fp.CurrValue(kNuisSuperkEnergyScale);
      this->GetSuperkDetectorFSI()->TweakEnergyScale(twk);
    }
    // rhc
    if(fp.Exists(kRHCNuisSuperkEnergyScale)) {
      double twk = fp.CurrValue(kRHCNuisSuperkEnergyScale);
      this->GetSuperkDetectorFSI()->TweakRHCEnergyScale(twk);
    }
#endif
    
    if(fp.Exists(kMultiRingNuisSuperkEnergyScale)) {
      double twk = fp.CurrValue(kMultiRingNuisSuperkEnergyScale);
      this->GetSuperkDetectorFSI()->TweakMultiRingEnergyScale(twk);
    }



    //////////////////// SKDetFSI ///////////////////

    
    // use current values stored in fit param list
    // to tweak the SKDET + FSI object

    //
    // N.B.
    // Don't need to use TweakFHCpRHC for correlations, unless Cholesky is performed.
    //

#if FSI_BINS == 2 // numu disappearance only

    // #ifdef __RHC_RUNNING__
    //     //Todo: fix this so that FHC and RHC can have different tweaks
    //     //Also that FHC and RHC can have different numbers of parameters
    //     const int nskdetfsi_1Rmu_rhc = conv::kNumOfSKDetFSIParams1Rmu; // same FHC
    //     TVectorD twkv_rhc(nskdetfsi_1Rmu_rhc);
    //     for(int ip = 0; ip < nskdetfsi_1Rmu_rhc; ip++) {
    //       PhysParamId_t p = PhysParamId_t ((int) kNuisSKDetFSI000 + ip); // same FHC
    //       twkv_rhc[ip] = fp.CurrValue(p);
    //     }
    //     this->GetSuperkDetectorFSI()->TweakOnlyRHC1Rmu(true,&twkv_rhc); // same FHC
    // #endif

    // fhc+rhc 
#if defined(__FHC_RUNNING__) && defined(__RHC_RUNNING__) 
    
    //fully correlated
    
    #if defined(__FULLY_CORRELATED__)
  
      const int nskdetfsi_1Rmu = conv::kNumOfSKDetFSIParams1Rmu;
      TVectorD twkv(nskdetfsi_1Rmu);
      for(int ip = 0; ip < nskdetfsi_1Rmu; ip++) {
	// Be sure kNuisSKDetFSI_first = kNuisSKDetFSI000 
	// kNuisSKDetFSI_first defined in T2KPhysParamId.h for each different analysis mode  
	PhysParamId_t p = PhysParamId_t ((int) kNuisSKDetFSI_first + ip);	
	twkv[ip] = fp.CurrValue(p);
      }
      this->GetSuperkDetectorFSI()->TweakOnly1Rmu(true,&twkv);
      this->GetSuperkDetectorFSI()->TweakOnlyRHC1Rmu(true,&twkv); // same FHC

      // //Todo: fix this so that FHC and RHC can have different tweaks
      // //Also that FHC and RHC can have different numbers of parameters
      // const int nskdetfsi_1Rmu_rhc = conv::kNumOfSKDetFSIParams1Rmu; // same FHC
      // TVectorD twkv_rhc(nskdetfsi_1Rmu_rhc);
      // for(int ip = 0; ip < nskdetfsi_1Rmu_rhc; ip++) {
      // 	// Be sure kNuisSKDetFSI_first = kRHCNuisSKDetFSI000 
      // 	// kNuisSKDetFSI_first defined in T2KPhysParamId.h for each different analysis mode  
      // 	PhysParamId_t p = PhysParamId_t ((int) kNuisSKDetFSI_first + ip); // same FHC
      // 	twkv_rhc[ip] = fp.CurrValue(p);
      // }
      // this->GetSuperkDetectorFSI()->TweakOnlyRHC1Rmu(true,&twkv_rhc); // same FHC
      

      // fhc/rhc correlations
    #else 

      // NEW FHC+RHC JOINT
      const int nskdetfsi_1Rmu = conv::kNumOfSKDetFSIParams1Rmu;
      TVectorD twkv(nskdetfsi_1Rmu);
      for(int ip = 0; ip < nskdetfsi_1Rmu; ip++) {
	// Be sure kNuisSKDetFSI_first = kNuisSKDetFSI000
	// kNuisSKDetFSI_first defined in T2KPhysParamId.h for each different analysis mode  
	PhysParamId_t p = PhysParamId_t ((int) kNuisSKDetFSI_first + ip);
	twkv[ip] = fp.CurrValue(p);
      }
      this->GetSuperkDetectorFSI()->TweakOnly1Rmu(true,&twkv);

      const int nskdetfsi_1Rmu_rhc = conv::kNumOfRHCSKDetFSIParams1Rmu; 
      TVectorD twkv_rhc(nskdetfsi_1Rmu_rhc);
      for(int ip = 0; ip < nskdetfsi_1Rmu_rhc; ip++) {
	// Be sure kRHCNuisSKDetFSI_first = kRHCNuisSKDetFSI000
	// kNuisSKDetFSI_first defined in T2KPhysParamId.h for each different analysis mode  
        PhysParamId_t p = PhysParamId_t ((int) kRHCNuisSKDetFSI_first + ip); 
	twkv_rhc[ip] = fp.CurrValue(p);
      }
      this->GetSuperkDetectorFSI()->TweakOnlyRHC1Rmu(true,&twkv_rhc); 
      //

    #endif
      
      // fhc only
#elif defined(__FHC_RUNNING__)
    const int nskdetfsi_1Rmu = conv::kNumOfSKDetFSIParams1Rmu;
    TVectorD twkv(nskdetfsi_1Rmu);
    for(int ip = 0; ip < nskdetfsi_1Rmu; ip++) {
      // Be sure kNuisSKDetFSI_first = kNuisSKDetFSI000 
      // kNuisSKDetFSI_first defined in T2KPhysParamId.h for each different analysis mode  
      PhysParamId_t p = PhysParamId_t ((int) kNuisSKDetFSI_first + ip);
      twkv[ip] = fp.CurrValue(p);
    }
    this->GetSuperkDetectorFSI()->TweakOnly1Rmu(true,&twkv);

    // rhc only
#elif defined(__RHC_RUNNING__)
    const int nskdetfsi_1Rmu_rhc = conv::kNumOfRHCSKDetFSIParams1Rmu; 
    TVectorD twkv_rhc(nskdetfsi_1Rmu_rhc);
    for(int ip = 0; ip < nskdetfsi_1Rmu_rhc; ip++) {
      // Be sure kRHCNuisSKDetFSI_first = kRHCNuisSKDetFSI000 
      // kNuisSKDetFSI_first defined in T2KPhysParamId.h for each different analysis mode  
      PhysParamId_t p = PhysParamId_t ((int) kRHCNuisSKDetFSI_first + ip); 
      twkv_rhc[ip] = fp.CurrValue(p);
    }
    this->GetSuperkDetectorFSI()->TweakOnlyRHC1Rmu(true,&twkv_rhc); 
#endif
    //

    // multiring
#ifdef __MULTIRING_FHC_RUNNING__
    const int nskdetfsi_MultiRmu = conv::kNumOfSKDetFSIParamsMultiRmu;
    TVectorD twkv_multiring(nskdetfsi_MultiRmu);
    for(int ip = 0; ip < nskdetfsi_MultiRmu; ip++) {   
      // Be sure kMultiRingNuisSKDetFSI_first = kMultiRingNuisSKDetFSI000 
      // kNuisSKDetFSI_first defined in T2KPhysParamId.h for each different analysis mode  
      PhysParamId_t p = PhysParamId_t ((int) kMultiRingNuisSKDetFSI_first + ip);
      twkv_multiring[ip] = fp.CurrValue(p);
    }
    this->GetSuperkDetectorFSI()->TweakOnlyMultiRmu(true,&twkv_multiring);
#endif
    ///

#else
#if FSI_BINS == 3 // nuebar appearance only

    // #ifdef __RHC_RUNNING__
    //       //Todo: fix this so that FHC and RHC can have different tweaks
    //       //Also that FHC and RHC can have different numbers of parameters
    //       const int nskdetfsi_1Re_rhc = conv::kNumOfSKDetFSIParams1Re; // same as FHC
    //       TVectorD twkv_rhc(nskdetfsi_1Re_rhc);
    //       for(int ip = 0; ip < nskdetfsi_1Re_rhc; ip++) {
    //       PhysParamId_t p = PhysParamId_t ((int) kNuisSKDetFSI006 + ip); // same as FHC
    //       twkv_rhc[ip] = fp.CurrValue(p);
    //       }
    //       this->GetSuperkDetectorFSI()->TweakOnlyRHC1Re(true,&twkv_rhc); // same as FHC
    // #endif

    // fhc+rhc 
#if defined(__FHC_RUNNING__) && defined(__RHC_RUNNING__) 

    //fully correlated
    #if defined(__FULLY_CORRELATED__)
    
    const int nskdetfsi_1Re = conv::kNumOfSKDetFSIParams1Re;
      TVectorD twkv(nskdetfsi_1Re);
      for(int ip = 0; ip < nskdetfsi_1Re; ip++) {
	// Be sure kNuisSKDetFSI_first = kNuisSKDetFSI006 
	// kNuisSKDetFSI_first defined in T2KPhysParamId.h for each different analysis mode  
	PhysParamId_t p = PhysParamId_t ((int) kNuisSKDetFSI_first + ip); 
        twkv[ip] = fp.CurrValue(p);
      }
      this->GetSuperkDetectorFSI()->TweakOnly1Re(true,&twkv);
      this->GetSuperkDetectorFSI()->TweakOnlyRHC1Re(true,&twkv); // same as FHC
      
      // //Todo: fix this so that FHC and RHC can have different tweaks
      // //Also that FHC and RHC can have different numbers of parameters
      // const int nskdetfsi_1Re_rhc = conv::kNumOfSKDetFSIParams1Re; // same as FHC
      // TVectorD twkv_rhc(nskdetfsi_1Re_rhc);
      // for(int ip = 0; ip < nskdetfsi_1Re_rhc; ip++) {
      // 	// Be sure kNuisSKDetFSI_first = kNuisSKDetFSI006 
      // 	// kNuisSKDetFSI_first defined in T2KPhysParamId.h for each different analysis mode  
      // 	PhysParamId_t p = PhysParamId_t ((int) kNuisSKDetFSI_first + ip); // same as FHC 
      // 	twkv_rhc[ip] = fp.CurrValue(p);
      // }
      // this->GetSuperkDetectorFSI()->TweakOnlyRHC1Re(true,&twkv_rhc); // same as FHC

      
      // fhc/rhc correlations
     #else

      // NEW FHC+RHC JOINT
      const int nskdetfsi_1Re = conv::kNumOfSKDetFSIParams1Re;
      TVectorD twkv(nskdetfsi_1Re);
      for(int ip = 0; ip < nskdetfsi_1Re; ip++) {
	// Be sure kNuisSKDetFSI_first = kNuisSKDetFSI006 
	// kNuisSKDetFSI_first defined in T2KPhysParamId.h for each different analysis mode  
        PhysParamId_t p = PhysParamId_t ((int) kNuisSKDetFSI_first + ip);
        twkv[ip] = fp.CurrValue(p);
      }
      this->GetSuperkDetectorFSI()->TweakOnly1Re(true,&twkv);

      const int nskdetfsi_1Re_rhc = conv::kNumOfRHCSKDetFSIParams1Re; 
      TVectorD twkv_rhc(nskdetfsi_1Re_rhc);
      for(int ip = 0; ip < nskdetfsi_1Re_rhc; ip++) {
	// Be sure kRHCNuisSKDetFSI_first = kRHCNuisSKDetFSI006 
	// kNuisSKDetFSI_first defined in T2KPhysParamId.h for each different analysis mode  
	PhysParamId_t p = PhysParamId_t ((int) kRHCNuisSKDetFSI_first + ip); 
	twkv_rhc[ip] = fp.CurrValue(p);
      }
      this->GetSuperkDetectorFSI()->TweakOnlyRHC1Re(true,&twkv_rhc); 
      //

   #endif
      //
      
      // fhc only
#elif defined(__FHC_RUNNING__)
    const int nskdetfsi_1Re = conv::kNumOfSKDetFSIParams1Re;
      TVectorD twkv(nskdetfsi_1Re);
      for(int ip = 0; ip < nskdetfsi_1Re; ip++) {
	// Be sure kNuisSKDetFSI_first = kNuisSKDetFSI006 
	// kNuisSKDetFSI_first defined in T2KPhysParamId.h for each different analysis mode  
	PhysParamId_t p = PhysParamId_t ((int) kNuisSKDetFSI_first + ip);
        twkv[ip] = fp.CurrValue(p);
      }
      this->GetSuperkDetectorFSI()->TweakOnly1Re(true,&twkv);
      
      // rhc only
#elif defined(__RHC_RUNNING__)
      const int nskdetfsi_1Re_rhc = conv::kNumOfRHCSKDetFSIParams1Re; 
      TVectorD twkv_rhc(nskdetfsi_1Re_rhc);
      for(int ip = 0; ip < nskdetfsi_1Re_rhc; ip++) {
	// Be sure kRHCNuisSKDetFSI_first = kRHCNuisSKDetFSI006 
	// kNuisSKDetFSI_first defined in T2KPhysParamId.h for each different analysis mode  
	PhysParamId_t p = PhysParamId_t ((int) kRHCNuisSKDetFSI_first + ip); 
      twkv_rhc[ip] = fp.CurrValue(p);
      }
      this->GetSuperkDetectorFSI()->TweakOnlyRHC1Re(true,&twkv_rhc); 
#endif
      
      // multiring
#ifdef __MULTIRING_FHC_RUNNING__
      const int nskdetfsi_MultiRe = conv::kNumOfSKDetFSIParamsMultiRe;
      TVectorD twkv_multiring(nskdetfsi_MultiRe);
      for(int ip = 0; ip < nskdetfsi_MultiRe; ip++) {
	// Be sure kMultiRingNuisSKDetFSI_first = kMultiRingNuisSKDetFSI006 
	// kNuisSKDetFSI_first defined in T2KPhysParamId.h for each different analysis mode  
	PhysParamId_t p = PhysParamId_t ((int) kMultiRingNuisSKDetFSI_first + ip);	
	twkv_multiring[ip] = fp.CurrValue(p);
      }
      this->GetSuperkDetectorFSI()->TweakOnlyMultiRe(true,&twkv_multiring);
#endif

#else // JOINT FIT
            
      // #ifdef __RHC_RUNNING__
      //       //Todo: fix this so that FHC and RHC can have different tweaks
      //       //Also that FHC and RHC can have different numbers of parameters
      //       const int nskdetfsi_rhc = conv::kNumOfSKDetFSIParams; // same as FHC
      //       TVectorD twkv_rhc(nskdetfsi_rhc);
      //       for(int ip = 0; ip < nskdetfsi_rhc; ip++) {
      //       PhysParamId_t p = PhysParamId_t ((int) kNuisSKDetFSI000 + ip); // same as FHC
      //       twkv_rhc[ip] = fp.CurrValue(p); 
      //       }
      //       // N.B. You have to tweak the FHC syst since RHC SKDetFSI doesn't exist!!!
      //       this->GetSuperkDetectorFSI()->Tweak(true,&twkv_rhc); // same as FHC    
      //       this->GetSuperkDetectorFSI()->TweakRHC(true,&twkv_rhc); // same as FHC    
      // #endif 
    
      // fhc+rhc 
#if defined(__FHC_RUNNING__) && defined(__RHC_RUNNING__) 
   
      //fully correlated
    #if defined(__FULLY_CORRELATED__)
      
      const int nskdetfsi = conv::kNumOfSKDetFSIParams;
      TVectorD twkv(nskdetfsi);
      for(int ip = 0; ip < nskdetfsi; ip++) {
	// Be sure kNuisSKDetFSI_first = kNuisSKDetFSI000
	// kNuisSKDetFSI_first defined in T2KPhysParamId.h for each different analysis mode  
	PhysParamId_t p = PhysParamId_t ((int) kNuisSKDetFSI_first + ip);
	twkv[ip] = fp.CurrValue(p); 
      }
      this->GetSuperkDetectorFSI()->Tweak(true,&twkv);    
      this->GetSuperkDetectorFSI()->TweakRHC(true,&twkv); // same as FHC    

      // //Todo: fix this so that FHC and RHC can have different tweaks
      // //Also that FHC and RHC can have different numbers of parameters
      // const int nskdetfsi_rhc = conv::kNumOfSKDetFSIParams; // same as FHC
      // TVectorD twkv_rhc(nskdetfsi_rhc);
      // for(int ip = 0; ip < nskdetfsi_rhc; ip++) {
      // 	// Be sure kNuisSKDetFSI_first = kRHCNuisSKDetFSI000
      // 	// kNuisSKDetFSI_first defined in T2KPhysParamId.h for each different analysis mode  
      // 	PhysParamId_t p = PhysParamId_t ((int) kNuisSKDetFSI_first + ip); // same as FHC	
      // 	twkv_rhc[ip] = fp.CurrValue(p); 
      // }
      // // N.B. You have to tweak the FHC syst since RHC SKDetFSI doesn't exist!!!
      // this->GetSuperkDetectorFSI()->Tweak(true,&twkv_rhc); // same as FHC    
      // this->GetSuperkDetectorFSI()->TweakRHC(true,&twkv_rhc); // same as FHC    


      // fhc/rhc with correlations
    #else 

      // NEW FHC+RHC JOINT

      const int nskdetfsi = conv::kNumOfSKDetFSIParams;
      TVectorD twkv(nskdetfsi);
      for(int ip = 0; ip < nskdetfsi; ip++) {
	// Be sure kNuisSKDetFSI_first = kNuisSKDetFSI000
	// kNuisSKDetFSI_first defined in T2KPhysParamId.h for each different analysis mode  
	PhysParamId_t p = PhysParamId_t ((int) kNuisSKDetFSI_first + ip);
	twkv[ip] = fp.CurrValue(p); 
      }
      this->GetSuperkDetectorFSI()->Tweak(true,&twkv);    
      const int nskdetfsi_rhc = conv::kNumOfRHCSKDetFSIParams; 
      TVectorD twkv_rhc(nskdetfsi_rhc);
      for(int ip = 0; ip < nskdetfsi_rhc; ip++) {
	// Be sure kRHCNuisSKDetFSI_first = kRHCNuisSKDetFSI000
	// kNuisSKDetFSI_first defined in T2KPhysParamId.h for each different analysis mode  
	PhysParamId_t p = PhysParamId_t ((int) kRHCNuisSKDetFSI_first + ip); 
	twkv_rhc[ip] = fp.CurrValue(p); 
      }
      this->GetSuperkDetectorFSI()->TweakRHC(true,&twkv_rhc); 
      // Update the merged FHC+RHC current values
      this->GetSuperkDetectorFSI()->UpdateFHCpRHCCurrentValues();

      //
      
    #endif

      // fhc only
#elif defined(__FHC_RUNNING__)
      const int nskdetfsi = conv::kNumOfSKDetFSIParams;
      TVectorD twkv(nskdetfsi);
      for(int ip = 0; ip < nskdetfsi; ip++) {
	// Be sure kNuisSKDetFSI_first = kNuisSKDetFSI000
	// kNuisSKDetFSI_first defined in T2KPhysParamId.h for each different analysis mode  
	PhysParamId_t p = PhysParamId_t ((int) kNuisSKDetFSI_first + ip);
       	twkv[ip] = fp.CurrValue(p); 
      }
      this->GetSuperkDetectorFSI()->Tweak(true,&twkv);    

      // rhc only
#elif defined(__RHC_RUNNING__)
      const int nskdetfsi_rhc = conv::kNumOfRHCSKDetFSIParams; 
      TVectorD twkv_rhc(nskdetfsi_rhc);
      for(int ip = 0; ip < nskdetfsi_rhc; ip++) {
	// Be sure kRHCNuisSKDetFSI_first = kRHCNuisSKDetFSI000
	// kNuisSKDetFSI_first defined in T2KPhysParamId.h for each different analysis mode  
	PhysParamId_t p = PhysParamId_t ((int) kRHCNuisSKDetFSI_first + ip); 
	twkv_rhc[ip] = fp.CurrValue(p); 
      }
      this->GetSuperkDetectorFSI()->TweakRHC(true,&twkv_rhc); // same as FHC    
#endif 
      //

      // multiring
#ifdef __MULTIRING_FHC_RUNNING__
      const int nskdetfsi_multiring = conv::kNumOfSKDetFSIParamsMultiRing;
      TVectorD twkv_multiring(nskdetfsi_multiring);
      for(int ip = 0; ip < nskdetfsi_multiring; ip++) {
	// Be sure kMultiRingNuisSKDetFSI_first = kMultiRingNuisSKDetFSI000 
	// kNuisSKDetFSI_first defined in T2KPhysParamId.h for each different analysis mode  
	PhysParamId_t p = PhysParamId_t ((int) kMultiRingNuisSKDetFSI_first + ip);
	twkv_multiring[ip] = fp.CurrValue(p); 
      }
      this->GetSuperkDetectorFSI()->TweakMultiRing(true,&twkv_multiring);    
#endif
      ///

#endif //FSI_BINS == 3
#endif //FSI_BINS == 2 
         

      ////////////////////// BANFF //////////////////////////
            
      // // use current values stored in fit param list
      // // to tweak the BANFF object
      // const int nbanff = conv::kNumOfParamsInBANFFfit;
      // TVectorD twkvBANFF(nbanff);
      // for(int ip = 0; ip < nbanff; ip++) {
      //   PhysParamId_t p = kListOfBANFFParams[ip];
      //   twkvBANFF[ip] = fp.CurrValue(p);
      // }
      // //Todo: fix this so that FHC and RHC can have different tweaks
      // //Also that FHC and RHC can have different numbers of parameters 
      // this->GetBANFF()->Tweak(true,&twkvBANFF);
      //   if(fp.GetSystMode() == 3){
      //   TVectorD twkvBANFF(nbanff);
      //   for(int ip = 0; ip < nbanff; ip++) {
      //       PhysParamId_t p = PhysParamId_t ((int) kRHCNuisBANFF00 + ip);
      //       twkvBANFF[ip] = fp.CurrValue(p);
      //   }
      //   this->GetBANFF()->TweakRHC(true,&twkvBANFF);
      //   }
      //   else this->GetBANFF()->TweakRHC(true,&twkvBANFF);    
   
   // fhc+rhc 
#if defined(__FHC_RUNNING__) && defined(__RHC_RUNNING__) 

   // fully correlated
   #if defined(__FULLY_CORRELATED__)

    // use current values stored in fit param list
    // to tweak the BANFF object
    const int nbanff = conv::kNumOfParamsInBANFFfit;
    TVectorD twkvBANFF(nbanff);
    for(int ip = 0; ip < nbanff; ip++) {
      PhysParamId_t p = kListOfBANFFParams[ip];
      twkvBANFF[ip] = fp.CurrValue(p);
    }
   //Todo: fix this so that FHC and RHC can have different tweaks
   //Also that FHC and RHC can have different numbers of parameters 
   this->GetBANFF()->Tweak(true,&twkvBANFF);
   if(fp.GetSystMode() == 3){
      TVectorD twkvBANFF(nbanff);
      for(int ip = 0; ip < nbanff; ip++) {
          PhysParamId_t p = PhysParamId_t ((int) kRHCNuisBANFF00 + ip);
          twkvBANFF[ip] = fp.CurrValue(p);
      }
      this->GetBANFF()->TweakRHC(true,&twkvBANFF);
      }
      else this->GetBANFF()->TweakRHC(true,&twkvBANFF);
   

   // fhc/rhc correlations
   #else

   // NEW FHC+RHC JOINT 

   //
   // Since the BANFF-XSec parameters are shared by FHC and RHC,
   // we use only the FHC ones parameters in the FHC+RHC vector.
   // e.g. RHC-MEC is not included (use FHC-enum) --> need to use TweakFHCpRHC()!!!
   // If Tweak() and TweakRHC() are used separately, it crashes.
   //
   // Use TweakFHCpRHC   
   //
   this->GetBANFF()->UpdateFHCpRHCCurrentValues();  
   const int nbanff = conv::kNumOfFHCpRHCParamsInBANFFfit;
   TVectorD twkvBANFF(nbanff);
   for(int ip = 0; ip < nbanff; ip++) {
    PhysParamId_t p = kListOfFHCpRHCBANFFParams[ip];
    twkvBANFF[ip] = fp.CurrValue(p);
   }
   this->GetBANFF()->TweakFHCpRHC(true,&twkvBANFF);
   // Split FHC+RHC into FHC-only and RHC-only
   this->GetBANFF()->UpdateFHCandRHConlyCurrentValues();   
   //
   #endif
   
   // fhc only
#elif defined(__FHC_RUNNING__)
    const int nbanff = conv::kNumOfParamsInBANFFfit;
    TVectorD twkvBANFF(nbanff);
    for(int ip = 0; ip < nbanff; ip++) {
      PhysParamId_t p = kListOfBANFFParams[ip];
      twkvBANFF[ip] = fp.CurrValue(p);
    }
    this->GetBANFF()->Tweak(true,&twkvBANFF);
    // rhc only
#elif defined(__RHC_RUNNING__)
    const int nbanff = conv::kNumOfRHCParamsInBANFFfit;
    TVectorD twkvBANFF(nbanff);
    for(int ip = 0; ip < nbanff; ip++) {
      PhysParamId_t p = kListOfRHCBANFFParams[ip];	
      twkvBANFF[ip] = fp.CurrValue(p);
    }      
    this->GetBANFF()->TweakRHC(true,&twkvBANFF);
#endif
    
    
    /////////// UNCORRELATED XSEC //////////
    
    // Cross-section systematics (other than BANFF ones)
    for(int ip = 0; ip < kNumOfXSecSyst; ip++) {
      PhysParamId_t p = kListOfXSecSyst[ip];
      if( fp.Exists(p) && util::IsXSec(p) ) {
	double twk = fp.CurrValue(p);
	this->GetNuXSec()->Tweak(p,twk);
      }
    }
    
    this->GetNuXSec()->SyncRHCToFHCSystematics(); // Uncorr Xsec syst are common to fhc and rhc
    
  }//systmode == 1 || systmode == 3
 
 else if(fp.GetSystMode() == 2) {
    cerr << "selected systematics not implemented in SuperkExpectation::SetNuisanceParams(). Exiting..." << endl;
    exit(1);
  }
}

//________________________________________________________________________________________
void SuperkExpectation::PrintSystematics (std::ostream & stream) const
{
  stream<<"SuperkExpectation/PrintSystematics not yet implemented"<<endl;
}
//________________________________________________________________________________________
double SuperkExpectation::GetFitStatisticPenalty(Option_t * /*opt*/) const
{
  bool all_syst=true;
  double chisq_penalty = 0.;
  //
  // Add penalty term for BANFF params 
  //
  
  //double dchisq_penalty_banff = fBANFF->ChisqPenalty() + fBANFF->RHCChisqPenalty();
    
  double dchisq_penalty_banff = 0.;

  // NEW FHC+RHC JOINT
#if defined(__FHC_RUNNING__) && defined(__RHC_RUNNING__) && !defined(__FULLY_CORRELATED__)
  dchisq_penalty_banff += fBANFF->FHCpRHCChisqPenalty();
#else
  //
  #if defined(__FHC_RUNNING__)
  dchisq_penalty_banff += fBANFF->ChisqPenalty();
  #endif
  #if defined(__RHC_RUNNING__)
  dchisq_penalty_banff += fBANFF->RHCChisqPenalty();
  #endif
  #if defined(__FHC_RUNNING__) && defined(__RHC_RUNNING__) && defined(__FULLY_CORRELATED__)
  dchisq_penalty_banff /= 2.;
  #endif
#endif

  chisq_penalty += dchisq_penalty_banff;
#ifdef __show_contributions_to_penalty_term__
  cout << "FitInputs::GetFitStatisticPenalty >> Current BANFF param variation"
       << " => \\delta\\chi^2_{penalty} = " << dchisq_penalty_banff << endl;
#endif

  //
  // Add penalty term for NuXSec params (other than the ones incl. in BANFF)
  //
  double dchisq_penalty_xsec =  fNuXSec->ChisqPenalty(all_syst);
  chisq_penalty += dchisq_penalty_xsec;
#ifdef __show_contributions_to_penalty_term__
  cout << "FitInputs::GetFitStatisticPenalty >> Current variation of params considered by NuXSec"
       << " => \\delta\\chi^2_{penalty} = " << dchisq_penalty_xsec << endl;
#endif
  
  //
  // Add penalty term for SuperkDetector + FSI params
  //

  //
  // OLD
  //
  // #if defined(__FHC_RUNNING__) || defined(__RHC_RUNNING__)
  //   // if the FHC is not defined, the penatly term returns 0 (same for RHC)
  //   double dchisq_penalty_skdetfsi = fSuperkDetectorFSI->ChisqPenalty() + fSuperkDetectorFSI->RHCChisqPenalty(); 
  //  
  // #if (defined(__FULLY_CORRELATED__) || defined(__BANFF_UNCORRELATED__)) && defined(__FHC_RUNNING__) && defined(__RHC_RUNNING__)
  //   dchisq_penalty_skdetfsi /= 2.;
  // #endif

  // NEW FHC+RHC JOINT
  double dchisq_penalty_skdetfsi = 0.;
#if defined(__FHC_RUNNING__) && defined(__RHC_RUNNING__)
  #if defined(__FULLY_CORRELATED__)
  dchisq_penalty_skdetfsi = fSuperkDetectorFSI->ChisqPenalty() + fSuperkDetectorFSI->RHCChisqPenalty(); 
  dchisq_penalty_skdetfsi /= 2.;
  #else
  dchisq_penalty_skdetfsi = fSuperkDetectorFSI->FHCpRHCChisqPenalty(); 
  #endif
#elif defined(__FHC_RUNNING__)
  dchisq_penalty_skdetfsi = fSuperkDetectorFSI->ChisqPenalty(); 
#elif defined(__RHC_RUNNING__)
  dchisq_penalty_skdetfsi = fSuperkDetectorFSI->RHCChisqPenalty(); 
#else
  cerr << endl;
  cerr << "No samples are defined" << endl;
  cerr << endl;
  exit(1);
#endif
  //
 
  chisq_penalty += dchisq_penalty_skdetfsi;
#ifdef __show_contributions_to_penalty_term__
  cout << "FitInputs::GetFitStatisticPenalty >> Current SK+FSI param variation"
       << " => \\delta\\chi^2_{penalty} = " << dchisq_penalty_skdetfsi << endl;
#endif

  // 
  // OLD
  //
  //#endif //defined(__FHC_RUNNING__) || defined(__RHC_RUNNING__)
  //

  // Add penalty term for SuperkDetector + FSI params MultiRing
#ifdef __MULTIRING_FHC_RUNNING__
  double dchisq_penalty_skdetfsi_multiring = fSuperkDetectorFSI->MultiRingChisqPenalty();
  chisq_penalty += dchisq_penalty_skdetfsi_multiring;
#ifdef __show_contributions_to_penalty_term__
  cout << "FitInputs::GetFitStatisticPenalty >> Current SK+FSI param variation MultiRing"
       << " => \\delta\\chi^2_{penalty} = " << dchisq_penalty_skdetfsi_multiring << endl;
#endif
#endif
  ///

  return chisq_penalty;
  }
//________________________________________________________________________________________
string valor::SuperkExpectation::Name(void)
{
  return "SuperkExpectation";
}
//________________________________________________________________________________________
