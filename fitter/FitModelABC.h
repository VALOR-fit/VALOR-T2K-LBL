//____________________________________________________________________________
/*!

\class     t2kana::FitModelABC

\brief     Pure abstract base class (pABC) for fit models.

\author    Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
           University of Liverpool & STFC Rutherford Appleton Laboratory

\created   March 30, 2014

\copyright (c) 2010-2014, The VALOR group

*/
//____________________________________________________________________________

#ifndef _t2kana_fit_model_abc_h_
#define _t2kana_fit_model_abc_h_

#include <string>

//~ #include "Conventions/PhysParamId.h"
//~ #include "FitCore/PhysHypothesis.h"
//~ #include "FitCore/PhysParamList.h"
#include "PhysParamId.h"
#include "PhysHypothesis.h"
#include "PhysParamList.h"

class TObjArray;

namespace t2kana {

class FitModelABC
{
public:
  virtual ~FitModelABC();

  // Identify fit model with a name (to be propagated in output files, etc)
  virtual string Name(void) const = 0;

  // Read configuration data from file
  virtual bool Read(string file, Option_t * opt = 0) = 0;

  // Below are a set methods to generate model predictions for a given physics hypothesis and a set of nuisance parameters.
  // The model is first tweaked appropriately by setting systematics to particular values or by randomizing them.
  // Then one uses the tweaked model to generate predictions for the given physics hypothesis.
  // The input option in TweakSystematics() can be used as needed in a particular context to interpret the input tweak
  // (for example, in T2K context, to specify whether the input tweak is wrt to the nominal and BANFF-tuned MC parameters).
  virtual bool Generate             (const PhysHypothesis & ph) = 0;
  virtual void RandomizeSystematics (void) = 0; 
  virtual void TweakSystematic      (string syst_name, double twkdial=0., Option_t * opt = 0) = 0; 
  virtual void TweakSystematic      (conv::PhysParamId_t syst, double twkdial=0., Option_t * opt = 0) = 0; 
  // Additional related methods to reset the systematics and enquire their values.
  virtual void ResetSystematics (void) = 0;  
  virtual void PrintSystematics (ostream & stream = std::cout) const = 0;
  virtual PhysParamList Tweaks(void) = 0;

  // This is a second interface for generating predictions for a given physics hypothesis and systematics.
  // This interface method was developed for fitters. 
  // The fit parameter list contains all relevant parameters, describing both a physics hypothesis and systematic variations.
  virtual bool Generate (const FitParamList & fp) = 0;

  // Write all predictions to the input container or to the input file
  virtual bool Write (TObjArray* collection, const char * tag = "model", Option_t * write_opt = 0) const = 0;
  virtual bool Write (const char * filename, const char * tag = "model", Option_t * write_opt = 0, Option_t * file_opt = 0) const = 0;

protected:
   FitModelABC();
};

} // t2kana

#endif

