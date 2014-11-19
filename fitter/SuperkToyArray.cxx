//________________________________________________________________________________________
/*
 
 \copyright (c) 2010-2014, VaLOR group.

*/
//________________________________________________________________________________________

#define __show_syst_reset_mesg__
#define __show_syst_tweak_mesg__
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

#include "SuperkToyArray.h"
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

ClassImp(SuperkToyArray)

//________________________________________________________________________________________
SuperkToyArray::SuperkToyArray(SuperkExpectation * master_skexp)
{
  for (int itoy = 0 ; itoy < NUMOFTOYSTOMARGINALIZE ; itoy++) {
    for (int is = 0 ; is < conv::kns ; is++) {
      fPredictErecoSpectrumArrayMode  [itoy][is] = 0;
      fPredictErecoSpectrumArrayTot   [itoy][is] = 0;
      fPredictSpectrumArrayOriginal   [itoy][is] = 0;
      fPredictSpectrumArrayWithPhysics[itoy][is] = 0;
    }
  }
    
  this->Init(master_skexp);
}
//________________________________________________________________________________________
SuperkToyArray::~SuperkToyArray(void)
{
  this->CleanUp();
}
//________________________________________________________________________________________
void SuperkToyArray::Init(SuperkExpectation * master_skexp)
{
  
  fMasterExpectation = master_skexp;
  
  SuperkMC * skmc = fMasterExpectation->GetSuperkMC();
  
  // Apply no physics hypothesis.
  PhysHypothesis ph(kHypVoid);
  
  for (int itoy = 0 ; itoy < NUMOFTOYSTOMARGINALIZE ; itoy++) {
    
    // Make a new toy.
    fMasterExpectation->RandomizeSystematics();
    fMasterExpectation->Generate(ph);
    
    for (int is = 0 ; is < conv::kns ; is++) {
      SuperkSample_t s = kSKSamples[is];
      if (s == conv::kSKsUndefined) continue;
        
      VALORArray * arr = new VALORArray(*skmc->CurrSpectrumVArray(s));
      
      fPredictSpectrumArrayOriginal   [itoy][is] = arr;
      fPredictSpectrumArrayWithPhysics[itoy][is] = new VALORArray(*arr);
      fPredictSpectrumArrayWithPhysics[itoy][is]->SetAll(0);
      
      int n_ereco_bins = kNErecoBins[is];
      
      fPredictErecoSpectrumArrayTot  [itoy][is] = new VALORArray(n_ereco_bins);
      fPredictErecoSpectrumArrayMode [itoy][is] = new VALORArray(conv::knm, n_ereco_bins);
      
      fPredictErecoSpectrumArrayTot  [itoy][is]->SetDimensionName(0,"Ereco");
      fPredictErecoSpectrumArrayMode [itoy][is]->SetDimensionName(0,"Mode");
      fPredictErecoSpectrumArrayMode [itoy][is]->SetDimensionName(1,"Ereco");
      
      fPredictSpectrumArrayOriginal   [itoy][is]->SetName(Form("spectrum_original_%s_%d"   ,util::AsStringCompact(s).c_str(),itoy));
      fPredictSpectrumArrayWithPhysics[itoy][is]->SetName(Form("spectrum_physics_%s_%d"    ,util::AsStringCompact(s).c_str(),itoy));
      fPredictErecoSpectrumArrayTot   [itoy][is]->SetName(Form("spectrum_ereco_total_%s_%d",util::AsStringCompact(s).c_str(),itoy));
      fPredictErecoSpectrumArrayMode  [itoy][is]->SetName(Form("spectrum_ereco_mode_%s_%d" ,util::AsStringCompact(s).c_str(),itoy));
    }
  }
  
  fMasterExpectation->ResetSystematics();
  
  return;
}
//________________________________________________________________________________________
void SuperkToyArray::Book(void)
{
  return;
}
//________________________________________________________________________________________
void SuperkToyArray::CleanUp(void)
{
  for (int itoy = 0 ; itoy < NUMOFTOYSTOMARGINALIZE ; itoy++) {
    for (int is = 0 ; is < conv::kns ; is++) {
      if (fPredictSpectrumArrayWithPhysics[itoy][is]) delete fPredictSpectrumArrayWithPhysics[itoy][is];
      if (fPredictSpectrumArrayOriginal   [itoy][is]) delete fPredictSpectrumArrayOriginal   [itoy][is];
      if (fPredictErecoSpectrumArrayTot   [itoy][is]) delete fPredictErecoSpectrumArrayTot   [itoy][is];
      if (fPredictErecoSpectrumArrayMode  [itoy][is]) delete fPredictErecoSpectrumArrayMode  [itoy][is];
    }
  }
}

//________________________________________________________________________________________
bool SuperkToyArray::Read(string /*config_file*/, Option_t * /*opt*/)
{
  // SuperkToyArray has no parameters to read for now.
  return true;
}
//________________________________________________________________________________________
bool SuperkToyArray::Generate(const PhysHypothesis & ph, Option_t * /*opt*/)
{
  // Calculate oscillation probabilities
  
  bool applied_negevents_condition = false;
  
  PhysHypothesis_t hyp = ph.Hypothesis();
      
  double p_osc[3][3];
  double p_osc_bar[3][3];

  //
  // 3-flavour oscillations
  //
  if(hyp == kHypOsc3F)
  {
    NuOsc3Flavour osc;
    // set oscillation parameters to the oscillation probability calculation code
    double mixing12 = ph.CurrValue(kPhypMixing12);
    double mixing13 = ph.CurrValue(kPhypMixing13);
    double mixing23 = ph.CurrValue(kPhypMixing23);
    //    double dm12sq   = ph.CurrValue(kPhypDm12sq);
    //double dm23sq   = ph.CurrValue(kPhypDm23sq);
    double dmsq     = ph.CurrValue(kPhypDm12sq);
    double Dmsq     = ph.CurrValue(kPhypDm23sq);
    double delta_cp = ph.CurrValue(kPhypCPV   );
#ifdef __show_hypothesis__
    cout << "**  Exp/SKMC >> " << fPhysHypothesis;
#endif
    NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
    NuOscDmsqScheme_t  dmsqscheme = Opt::Instance()->gOscDmsqScheme;
  
    // mixing params are: sin^2(2*theta)
    // independent squared mass splittings are: dm^2 = Dm^2_{12}, Dm^2 = |Dm^2_{23}|                                                                                             
    if(dmsqscheme == kNuDmsq12andAbs23 && omp == kNuMixSinsq2Theta) {
      osc.SetOscParameters(mixing12,mixing13,mixing23,dmsq,Dmsq,delta_cp);
    }
    // mixing params are: sin^2(theta)                                                                                                                                         
    // independent squared mass splittings are: dm^2 = Dm^2_{12}, Dm^2 = |Dm^2_{23}|                                                                                              
    else if(dmsqscheme == kNuDmsq12andAbs23 && omp == kNuMixSinsqTheta) {
      osc.SetOscParametersSinSqTAll(mixing12,mixing13,mixing23,dmsq,Dmsq,delta_cp);
    }
    // mixing params are: sin^2(2*theta)                                                                                                                                        
    // independent squared mass splittings are: dm^2 = Dm^2_{12}, Dm^2 = |m^2_{3} - (m^2_{1} + m^2_{1})/2|                                                                        
    else if(dmsqscheme == kNuDmsqFogliLisi && omp == kNuMixSinsq2Theta) {
      osc.SetOscParametersSSq2TDmAtmSol(mixing12,mixing13,mixing23,dmsq,Dmsq,delta_cp);
    }
    // mixing params are: sin^2(theta)                                                                                                                                             
    // independent squared mass splittings are: dm^2 = Dm^2_{12}, Dm^2 = |m^2_{3} - (m^2_{1} + m^2_{1})/2|                                                                         
    else if(dmsqscheme == kNuDmsqFogliLisi && omp == kNuMixSinsqTheta) {
      osc.SetOscParametersDmsqAtmSol(mixing12,mixing13,mixing23,dmsq,Dmsq,delta_cp);
    }
    // mixing paras are: sin^2(2*theta)
    // independent squared mass splittings are: dm^2 = Dm^2_{12}, Dm^2 = |Dm^2_{31}|
    else if(dmsqscheme == kNuDmsq12andAbs31 && omp == kNuMixSinsq2Theta) {
      osc.SetOscParametersSSq2TDmsq13(mixing12,mixing13,mixing23,dmsq,Dmsq,delta_cp);
    }
    // mixing paras are: sin^2(theta)
    // independent squared mass splittings are: dm^2 = Dm^2_{12}, Dm^2 = |Dm^2_{31}|
    else if(dmsqscheme == kNuDmsq12andAbs31 && omp == kNuMixSinsqTheta) {
      osc.SetOscParametersDmsq13(mixing12,mixing13,mixing23,dmsq,Dmsq,delta_cp);
    }
    //unknown
    else {
      cerr << "**  Exp/SKMC >> Can't handle choice of "
           << "\n oscillation mixing params: " << util::AsString(omp)
           << "\n and independent squared mass-splittings: " << util::AsString(dmsqscheme)
           << endl;
      exit(1);
    }
  
    // loop over samples
    for(int is = 0; is < kns; is++) {
      SuperkSample_t s = kSKSamples[is];
      if(s == kSKsUndefined) continue;
      
      int n_ereco_bins = kNErecoBins[is];
      int n_etrue_bins = kNEtrueBins[is];
        
      // loop over Etrue
      for(int it = 1; it <= n_etrue_bins; it++) {
         
        // calculate oscillation probabilities (incl. matter effects) for the current Etrue value
        double Etrue = (util::EtrueBinEdges(s,it-1) + util::EtrueBinEdges(s,it) ) / 2;
         
        osc.CalcMatterOscProbAll(false, Etrue, p_osc    );
        osc.CalcMatterOscProbAll(true , Etrue, p_osc_bar);
          
        double p_numu_numu       = p_osc    [1][1];
        double p_numubar_numubar = p_osc_bar[1][1];
        double p_nue_nue         = p_osc    [0][0];
        double p_nuebar_nuebar   = p_osc_bar[0][0];
        double p_numu_nue        = p_osc    [1][0];
        double p_numubar_nuebar  = p_osc_bar[1][0];
        double p_nue_numu        = p_osc    [0][1];
        double p_nuebar_numubar  = p_osc_bar[0][1];

        // Loop over the toy experiments
        for (int itoy = 0 ; itoy < NUMOFTOYSTOMARGINALIZE ; itoy++) {
          
          // loop over modes
          for(int im=0; im < knm; im++) {
            SuperkMode_t m = kSKReactionModes[is][im];
            if(m == kSKmUndefined) continue;
            
            // Get a pointer to the beginning of the ereco array for this toy/sample/mode/truebin
            const double * orig_array = &(fPredictSpectrumArrayOriginal   [itoy][is]->GetEntryRef(im,it-1,0));
            double * physics_array    = &(fPredictSpectrumArrayWithPhysics[itoy][is]->GetEntryRef(im,it-1,0));
            // osc. probability
            double posc = 1.;
            if(kIsCC[is][im])
            {
              if     ( util::IsNuMu(m)      ) { posc = p_numu_numu;       }
              else if( util::IsNuMuBar(m)   ) { posc = p_numubar_numubar; }
              else if( util::IsNuE(m)       ) { posc = p_nue_nue;         } 
              else if( util::IsNuEBar(m)    ) { posc = p_nuebar_nuebar;   }
              else if( util::IsOscNuMu(m)   ) { posc = p_nue_numu;        }
              else if( util::IsOscNuMuBar(m)) { posc = p_nuebar_numubar;  }
              else if( util::IsOscNuE(m)    ) { posc = p_numu_nue;        }
              else if( util::IsOscNuEBar(m) ) { posc = p_numubar_nuebar;  }
            }
            // loop over Ereco
            for(int ir = 1; ir <= n_ereco_bins ; ir++) {
              // number of events in (Ereco,Etrue) bin (systematics applied)
              double nev = 0;
              
              nev = orig_array[ir-1];
              
#ifdef PHYSHYPWEIGHTSAREONE
              posc = 1;
#ifdef PHYSHYPWEIGHTAREONEANDNUEISZERO
              if(util::IsOscNuE(m) || util::IsOscNuEBar(m) || util::IsOscNuMuBar(m) || util::IsOscNuMu(m) ) posc = 0;
#endif
#endif
              // number of events in (Ereco,Etrue) bin (oscillations + systematics applied)
              double nev_osc = nev * posc;
              
              physics_array[ir-1] = nev_osc;
            } // ireco
          }//imode
        }//itoy
      }//itrue
    }//isample
  }//hyp==NuOsc3f
  else {
    std::cerr<<"Unimplemented hypothesis: "<<util::AsString(hyp);
    exit(-1);
  }
  
  this->UpdateErecoSpectra();
  
  return applied_negevents_condition;
}
//________________________________________________________________________________________
bool SuperkToyArray::Generate(const FitParamList & fp, Option_t * /*opt*/)
{
  // extract the values of params driving the oscillation, decoherence or decay model
  PhysHypothesis ph(fp);
  
  return Generate(ph);
}
//________________________________________________________________________________________
void SuperkToyArray::UpdateErecoSpectra(void)
{
  for (int itoy = 0 ; itoy < this->GetNToys() ; itoy++) {
    // loop over samples
    for(int is = 0; is < kns; is++) {
      SuperkSample_t s = kSKSamples[is];
      if(s == kSKsUndefined) continue;

      int n_ereco_bins = kNErecoBins[is];
      int n_etrue_bins = kNEtrueBins[is];
      
      // Recast VALORArray blocks so they can be used addressed like normal arrays (fast).
      // Hard to read, so should only be used in code for which speed is essential.
      typedef double (templmtype)[knm][n_etrue_bins][n_ereco_bins];
      typedef double (ermtype)   [knm][n_ereco_bins];
      const templmtype & templFin = *reinterpret_cast< const templmtype * > (fPredictSpectrumArrayWithPhysics [itoy][is]->GetArray());
      ermtype & erecoMode         = *reinterpret_cast< ermtype          * > (fPredictErecoSpectrumArrayMode   [itoy][is]->GetArray());
      
      double * erecoTot = fPredictErecoSpectrumArrayTot[itoy][is]->GetArray();

      // loop over Ereco
      for(int ir = 1; ir <= n_ereco_bins; ir++) {
        double nev_tot = 0;
        // loop over modes
        for(int im = 0; im < knm; im++) { 
          SuperkMode_t m = kSKReactionModes[is][im];
          if(m == kSKmUndefined) continue;
          
          // loop over Etrue
          double nev = 0;
          for(int it = 1; it <= n_etrue_bins; it++) {
            nev += templFin[im][it-1][ir-1];
          }//it
          nev_tot += nev;
#ifdef __apply_minevents_to_individual_modes__
          if(nev < 0) {
#ifdef __show_negative_events_mesg_debug__
          cout << "Exp/SuperKToyArray/UpdateErecoSpectra >> set negative number of events " << nev
               << " to zero " << 0 << "\t" << " Ereco bin: " << ir << endl;
#endif
            nev = 0;
          } // if nev < 0
#endif //__apply_minevents_to_individual_modes__

          erecoMode[im][ir-1] = nev;   
        }//im
        
        erecoTot[ir-1] = nev_tot;
        
      }//ir
    }//is
    //applied_negevents_condition = this->ApplyEnergyScaleSystematicToTotalSpectrum();
  }
  
  return;
}
//________________________________________________________________________________________
TH1D* SuperkToyArray::CurrPredictedErecoSpectrum(int itoy, valor::conv::SuperkSample_t s) const 
{
  const double * arr = this->CurrPredictedErecoSpectrumArray(itoy, s);
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
TH1D* SuperkToyArray::CurrPredictedErecoSpectrum(
  int itoy, valor::conv::SuperkSample_t s, valor::conv::SuperkMode_t m) const 
{ 
  const double * arr = this->CurrPredictedErecoSpectrumArrayMode(itoy, s,m);
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
VALORArray * SuperkToyArray::CurrPredictedErecoVArray(int itoy, valor::conv::SuperkSample_t s) 
{
  assert(itoy >= 0 && itoy < NUMOFTOYSTOMARGINALIZE);
  VALORArray * arr = fPredictErecoSpectrumArrayTot[itoy][s];
  if (!arr) {
    std::cerr<<"SKToyArray/CurrPredictedErecoVArray could not find array for toy "<<itoy <<" and sample "<<util::AsString(s)<<std::endl;
    exit(-1);
  }
  
  return arr;
}
//________________________________________________________________________________________
VALORArray * SuperkToyArray::CurrPredictedErecoVArrayMode(int itoy, valor::conv::SuperkSample_t s) 
{
  assert(itoy >= 0 && itoy < NUMOFTOYSTOMARGINALIZE);
  VALORArray * arr = fPredictSpectrumArrayWithPhysics[itoy][s];
  if (!arr) {
    std::cerr<<"SKToyArray/CurrPredictedErecoVArraMode could not find array for toy "<<itoy <<" and sample "<<util::AsString(s)<<std::endl;
    exit(-1);
  }
  
  return arr;
}
//________________________________________________________________________________________
const VALORArray * SuperkToyArray::CurrPredictedErecoVArray(int itoy, valor::conv::SuperkSample_t s) const
{
  assert(itoy >= 0 && itoy < NUMOFTOYSTOMARGINALIZE);
  const VALORArray * arr = fPredictErecoSpectrumArrayTot[itoy][s];
  if (!arr) {
    std::cerr<<"SKToyArray/CurrPredictedErecoVArray const could not find array for toy "<<itoy <<" and sample "<<util::AsString(s)<<std::endl;
    exit(-1);
  }
  
  return arr;
}
//________________________________________________________________________________________
const VALORArray * SuperkToyArray::CurrPredictedErecoVArrayMode(int itoy, valor::conv::SuperkSample_t s) const
{
  assert(itoy >= 0 && itoy < NUMOFTOYSTOMARGINALIZE);
  const VALORArray * arr = fPredictErecoSpectrumArrayMode[itoy][s];
  if (!arr) {
    std::cerr<<"SKToyArray/CurrPredictedErecoVArrayMode const could not find array for toy "<<itoy <<" and sample "<<util::AsString(s)<<std::endl;
    exit(-1);
  }
  
  return arr;
}
//________________________________________________________________________________________
double* SuperkToyArray::CurrPredictedErecoSpectrumArray(int itoy, valor::conv::SuperkSample_t s) 
{
  VALORArray * arr = this->CurrPredictedErecoVArray(itoy,s);
  return arr->GetArray();
}
//________________________________________________________________________________________
const double* SuperkToyArray::CurrPredictedErecoSpectrumArray(int itoy, valor::conv::SuperkSample_t s)  const
{
  const VALORArray * arr = this->CurrPredictedErecoVArray(itoy, s);
  return arr->GetArray();
}
//________________________________________________________________________________________
double* SuperkToyArray::CurrPredictedErecoSpectrumArrayMode( int itoy, 
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
  
  VALORArray * arr = this->CurrPredictedErecoVArrayMode(itoy, s);

  return &(arr->GetEntryRef(midx,0));
}
//________________________________________________________________________________________
const double* SuperkToyArray::CurrPredictedErecoSpectrumArrayMode(int itoy,
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
  
  const VALORArray * arr = this->CurrPredictedErecoVArrayMode(itoy, s);

  return &(arr->GetEntryRef(midx,0));
}
//________________________________________________________________________________________
const VALORArray * SuperkToyArray::CurrPredictedVArrayMode (int itoy, valor::conv::SuperkSample_t s) const
{
  int sidx = util::ArrIdx::Instance()->SuperkSampleIdx(s);
  if(sidx == -1) {
    return 0;
  }
  
  return fPredictSpectrumArrayWithPhysics[itoy][sidx];
}
//________________________________________________________________________________________
VALORArray *       SuperkToyArray::CurrPredictedVArrayMode (int itoy, valor::conv::SuperkSample_t s)
{
  int sidx = util::ArrIdx::Instance()->SuperkSampleIdx(s);
  if(sidx == -1) {
    return 0;
  }
  
  return fPredictSpectrumArrayWithPhysics[itoy][sidx];
}
//________________________________________________________________________________________

void SuperkToyArray::PrintSystematicVariations(std::ostream & stream) const
{
  stream<<"SuperkToyArray/PrintSystematicVariations not implemented"<<std::endl;
}
//________________________________________________________________________________________
bool SuperkToyArray::Write(TObjArray * array, const char * pdf_tag, Option_t * write_opt)
{
  if(!array) {
    cout << "SKExp >> Null TObjArray. Not writing data" << endl;
    return false;
  }
  
  bool hist_mode = (write_opt && string(write_opt).find("hist")!=string::npos);
  
  string append_begin(Form("%s_", pdf_tag));
  
  for (int itoy = 0 ; itoy < NUMOFTOYSTOMARGINALIZE ; itoy++ ) {
    for(int is=0; is<kns; is++) {
      SuperkSample_t s = kSKSamples[is];
      if(s == kSKsUndefined) continue;

      if (hist_mode) {
        util::SetDefaultStyle();
        TH1D * spectrum = 0;
        TH1D * spectrum_mode[knm] = { };
        
        // copy & format Ereco spectra
        spectrum = this->CurrPredictedErecoSpectrum(itoy,s);
        spectrum->SetDirectory (0);
        
        for(int im = 0; im < knm; im++) {
           SuperkMode_t m = kSKReactionModes[is][im];
           if(m == kSKmUndefined) continue;
           spectrum_mode[im] = this->CurrPredictedErecoSpectrum(itoy,s,m);
           spectrum_mode[im]->SetDirectory (0);
        }//imode
        
        // add spectra to the results array
        spectrum->SetName( Form("%s_%d",util::AsStringCompact(s,kSKmUndefined,append_begin).c_str(),itoy));
        array->Add(spectrum);
        for(int im = 0; im < knm; im++) {
           SuperkMode_t m = kSKReactionModes[is][im];
           if(m == kSKmUndefined) continue;
           spectrum_mode[im]->SetName(Form("%s_%d",util::AsStringCompact(s,m,append_begin).c_str(),itoy));
           array->Add(spectrum_mode[im]);
        }//im
      } // hist mode
      else {
        VALORArray * total = this->CurrPredictedErecoVArray(itoy,s);
        array->Add(total);
        VALORArray * bymode = this->CurrPredictedErecoVArrayMode(itoy,s);
        array->Add(bymode);
      }
    }//is
  } //itoy

  return true;
}
//________________________________________________________________________________________
bool SuperkToyArray::Write(
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
void SuperkToyArray::ResetSystematics(void)
{
  cout<<"SuperkToyArray/ResetSystematics not yet implemented"<<endl;
}
//________________________________________________________________________________________
void SuperkToyArray::RandomizeSystematics(void)
{
  cout<<"SuperkToyArray/RandomizeSystematics not yet implemented"<<endl;
  return;
}
//________________________________________________________________________________________
void SuperkToyArray::PrintSystematics (std::ostream & stream) const
{
  stream<<"SuperkToyArray/PrintSystematics not yet implemented"<<endl;
}
//________________________________________________________________________________________
void SuperkToyArray::TweakSystematic      (string /*syst_name*/, double /*twkdial*/, Option_t * /*opt*/)
{
  cout<<"SuperkToyArray/TweakSystematic not yet implemented"<<endl;
  return;
}
//________________________________________________________________________________________
void SuperkToyArray::TweakSystematic      (conv::PhysParamId_t /*syst*/, double /*twkdial*/, Option_t * /*opt*/)
{
  cout<<"SuperkToyArray/TweakSystematic not yet implemented"<<endl;
  return;
}
//________________________________________________________________________________________
PhysParamList SuperkToyArray::Tweaks (void)
{
  cout<<"SuperkToyArray/Tweaks not yet implemented"<<endl;
  return PhysParamList();
}
//________________________________________________________________________________________
  
double SuperkToyArray::GetFitStatisticPenalty(Option_t * /*opt*/) const
{
  double chisq_penalty = 0.;
  cout<<"SuperkToyArray/GetFitStatisticPenalty not implemented, returning zero."<<endl;
  return chisq_penalty;
}
//________________________________________________________________________________________
string valor::SuperkToyArray::Name(void)
{
  return "SuperkToyArray";
}
//________________________________________________________________________________________
int valor::SuperkToyArray::GetNToys(void) const
{
  return NUMOFTOYSTOMARGINALIZE;
}
//________________________________________________________________________________________
SuperkExpectation * valor::SuperkToyArray::GetMasterExpectation(void){
  return fMasterExpectation;
}
//________________________________________________________________________________________
const SuperkExpectation * valor::SuperkToyArray::GetMasterExpectation(void) const
{
  return fMasterExpectation;
}
//________________________________________________________________________________________
