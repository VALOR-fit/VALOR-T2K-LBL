//________________________________________________________________________________________
/*

*/
//________________________________________________________________________________________

#include <cassert>

#include <TAxis.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TSpline.h>

#include "year1_numu_osc_analysis/SuperkExpectation.h"

using namespace t2kana;
using namespace t2kana::year1_numuosc;

ClassImp(SuperkExpectation)

//________________________________________________________________________________________
SuperkExpectation::SuperkExpectation(void)
{
  this->Init();
}
//________________________________________________________________________________________
SuperkExpectation::~SuperkExpectation(void)
{
  this->CleanUp();
}
//________________________________________________________________________________________
bool SuperkExpectation::Generate(double Er)
{
//
//
//

/*
   for(int it = 0; it < fSuperkInputs->NEnergyBins(); it++) {



   //
   // calculate get dN^{total}_{nd280} (Et) / dEt
   //
   
   ntot_nd280_t = 0;

   for(int ir = 0; ir < fND280Inputs->NEnergyBins(); ir++) {

     double Erp  = fND280Inputs->BinCentre(ir);
     double dErp = fND280Inputs->BinWidth(ir);

     // get dN^{observed}_{nd280} (Er') / dEr'
     double nobs_nd280_r  = fND280Inputs -> NEvents(Erp);

     // get nd280 efficiency (Er')
     double eff_nd280  = fND280Inputs -> Efficiency(Erp);

     // calculate get dN^{total}_{nd280} (Er') / dEr'
     double ntot_nd280_r = nobs_nd280_r / eff_nd280;

     // get nd280 energy resolution matrix
     double Rnd280  = fND280Inputs->EnergyResolutionMatrix(Et, Erp);

     ntot_nd280_t += (dErp * Rnd280 * ntot_nd280_r);
  }

   //
   // calculate dFlux_{nd280}(Et)/dEt
   //

   double mass_nd280    = fND280Inputs -> FiducialMass();
   double pot_nd280     = fND280Inputs -> POT();
   double xsec_nd280    = fGeneratorInputs -> XXX();
   double flux_nd280    = ntot_nd280_t / (mass_nd280 * pot_nd280 * xsec_nd280);

   //
   // calculate dFlux_{SK}(Et)/dEt
   //

   double fn_ratio = fBeamInputs->FNRatio(Et);
   double flux_sk  = flux_nd280 * fn_ratio;

   //
   // calculate dN^{total;unoscillated}_{SK}(Et)/dEt
   //

   double mass_sk    = fSuperkInputs -> FiducialMass();
   double pot_sk     = fSuperkInputs -> POT();
   double xsec_water = fGeneratorInputs -> XXX();
   double evt_frac   = fGeneratorInputs -> XXX();

   double ntot_unosc_sk_t = flux_sk * xsec_water * evt_frac * mass_sk * pot_sk;

   //
   // calculate dN^{total;oscillated}_{SK}(Et)/dEt
   //
   double prob_survival = XXX;
   ntot_sk_t = ntot_unosc_sk_t * prob_survival;
*/
  return 0;
}
//________________________________________________________________________________________
void SuperkExpectation::Update(PhysParamList & fp)
{
// update SK energy spectrum expectation

}
//________________________________________________________________________________________
double SuperkExpectation::Chisquare(SuperkObservation & skobs, Option_t * opt)
{
// judge the agreement between the input data and the current prediction

  return 0;
}
//________________________________________________________________________________________
void SuperkExpectation::Init(void)
{
  fND280Inputs     = new ND280Inputs;
  fSuperkInputs    = new SuperkInputs;
  fBeamInputs      = new BeamInputs; 
  fGeneratorInputs = new GeneratorInputs; 
}
//________________________________________________________________________________________
void SuperkExpectation::CleanUp(void)
{
  delete fND280Inputs;
  delete fSuperkInputs;
  delete fBeamInputs; 
  delete fGeneratorInputs;
}
//________________________________________________________________________________________
//........................................................................................
// nested class ND280Inputs
//........................................................................................
//________________________________________________________________________________________
SuperkExpectation::ND280Inputs::ND280Inputs()
{
  this->Init();
}
//________________________________________________________________________________________
SuperkExpectation::ND280Inputs::~ND280Inputs()
{
  this->CleanUp();
}
//________________________________________________________________________________________
void SuperkExpectation::ND280Inputs::Init(void)
{
  fPOT             = 0;
  fFidMass         = 0;
  fErBinning       = 0;
  fEtBinning       = 0;
  fObsCCLikeEr     = 0;
  fEstNCBkgEr      = 0;
  fRecoEffCC       = 0;
  fSelEffCC        = 0;
  fResolMtrxCC     = 0;
}
//________________________________________________________________________________________
void SuperkExpectation::ND280Inputs::CleanUp(void)
{
  if ( fErBinning   ) delete fErBinning;
  if ( fEtBinning   ) delete fEtBinning;
  if ( fObsCCLikeEr ) delete fObsCCLikeEr;
  if ( fEstNCBkgEr  ) delete fEstNCBkgEr;
  if ( fRecoEffCC   ) delete fRecoEffCC;
  if ( fSelEffCC    ) delete fSelEffCC;
  if ( fResolMtrxCC ) delete fResolMtrxCC;
}
//________________________________________________________________________________________
bool SuperkExpectation::ND280Inputs::Read(string config_file)
{
  return true;
}
//________________________________________________________________________________________
//........................................................................................
// nested class SuperkInputs
//........................................................................................
//________________________________________________________________________________________
SuperkExpectation::SuperkInputs::SuperkInputs()
{
  this->Init();
}
//________________________________________________________________________________________
SuperkExpectation::SuperkInputs::~SuperkInputs()
{
  this->CleanUp();
}
//________________________________________________________________________________________
void SuperkExpectation::SuperkInputs::Init(void)
{
  fPOT       = 0;  
  fFidMass   = 0; 
  fErBinning = 0;
  fEtBinning = 0;

  for(int im=0; im<knm; im++) {
    fRecoEff[im] = 0;
  }
  for(int is=0; is<kns; is++) {
    for(int im=0; im<knm; im++) {
      fSelEff   [is][im] = 0; 
      fResolMtrx[is][im] = 0; 
    }
  }
}
//________________________________________________________________________________________
void SuperkExpectation::SuperkInputs::CleanUp(void)
{
  if (fErBinning ) delete fErBinning;
  if (fEtBinning ) delete fEtBinning;

  for(int im=0; im<knm; im++) {
    if(fRecoEff[im]) delete fRecoEff[im];
  }
  for(int is=0; is<kns; is++) {
    for(int im=0; im<knm; im++) {
      if ( fSelEff   [is][im] ) delete  fSelEff   [is][im];
      if ( fResolMtrx[is][im] ) delete  fResolMtrx[is][im];
    }
  }
}
//________________________________________________________________________________________
bool SuperkExpectation::SuperkInputs::Read(string config_file)
{
  return true;
}
//________________________________________________________________________________________
//........................................................................................
// nested class BeamInputs
//........................................................................................
//________________________________________________________________________________________
SuperkExpectation::BeamInputs::BeamInputs()
{

}
//________________________________________________________________________________________
SuperkExpectation::BeamInputs::~BeamInputs()
{

}
//________________________________________________________________________________________
//........................................................................................
// nested class GeneratorInputs
//........................................................................................
//________________________________________________________________________________________
SuperkExpectation::GeneratorInputs::GeneratorInputs()
{
  this->Init();
}
//________________________________________________________________________________________
SuperkExpectation::GeneratorInputs::~GeneratorInputs()
{
  this->CleanUp();
}
//________________________________________________________________________________________
void SuperkExpectation::GeneratorInputs::Init(void)
{
  for(int im=0; im<knm; im++) {
    fModeFraction[im] = 0;
  }
  fSigmaTotWater   = 0;
  fSigmaTotCCND280 = 0;
}
//________________________________________________________________________________________
void SuperkExpectation::GeneratorInputs::CleanUp(void)
{
  for(int im=0; im<knm; im++) {
    if ( fModeFraction[im] ) delete fModeFraction[im];
  }
  if ( fSigmaTotWater   ) delete fSigmaTotWater;
  if ( fSigmaTotCCND280 ) delete fSigmaTotCCND280;
}
//________________________________________________________________________________________
bool SuperkExpectation::GeneratorInputs::Read(string config_file)
{
  return true;
}
//________________________________________________________________________________________


