//____________________________________________________________________________
/*!

\class    t2kana::SuperkObservation

\brief    

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  October 20, 2010

*/
//____________________________________________________________________________

#ifndef _t2kana_sk_observation_h_
#define _t2kana_sk_observation_h_

#include <string>

class TH1D;
class TVirtualPad;

#include "Conventions.h"
#include "SuperkExpectation.h"
#include "FitParamList.h"

using std::string;
using t2kana::conv::kns;

namespace t2kana {

class SuperkObservation: public TObject
{
public:
  SuperkObservation();
 ~SuperkObservation();

  // Read data
  bool Read(string config_file, int data_set_id=0);

  // Create mock data
  void CreateMockData(SuperkExpectation & exp, const FitParamList & fp,   bool fluctuate=true);
  void CreateMockData(SuperkExpectation & exp, const PhysHypothesis & ph, bool fluctuate=true);

  // Get data
  double NumOfEvents (t2kana::conv::SuperkSample_t s = t2kana::conv::kSKsUndefined) const;             
  double NumOfEvents (double Er, t2kana::conv::SuperkSample_t s = t2kana::conv::kSKsUndefined) const;  

  TH1D*  ErecoSpectrum (t2kana::conv::SuperkSample_t s) const; 

  // Draw.
  // Set option to "1RmuEreco"  to show the reconstructed Ereco distributions
  // Set option to ... to show ...
  void Draw(TVirtualPad * c, Option_t * opt = "Ereco") const;

private:

  void Init    (void);
  void CleanUp (void);

  TH1D *  fErDistr[kns];  ///< Reconstructed energy spectrum for each sample

ClassDef(SuperkObservation, 3)
};

} // t2kana

#endif 

