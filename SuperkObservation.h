//____________________________________________________________________________
/*!

\class    t2kana::year1_numuosc::SuperkObservation

\brief    

\author   STFC, Rutherford Appleton Laboratory / Oxford T2K groups

          Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          Jelena Ilic         <Jelena.Ilic \at stfc.ac.uk>
          Nick Grant          <nick.grant \at stfc.ac.uk>
          Geoff Pearce        <geoff.pearce \at stfc.ac.uk>
          Alfons Weber        <alfons.weber \at stfc.ac.uk>

\created  October 20, 2010

*/
//____________________________________________________________________________

#ifndef _t2kana_year1_numu_analysis_sk_observation_h_
#define _t2kana_year1_numu_analysis_sk_observation_h_

#include <string>

class TAxis;
class TH1D;
class TVirtualPad;

#include "year1_numu_osc_analysis/Conventions.h"

using std::string;
//using t2kana::year1_numuosc::conventions::SuperkSample_t;
using t2kana::year1_numuosc::conventions::kns;

namespace t2kana {
namespace year1_numuosc {

class SuperkObservation: public TObject
{
public:
  SuperkObservation();
 ~SuperkObservation();

  // Read data
  bool Read(string file);

  // Get data
  double NumOfEvents (t2kana::year1_numuosc::conventions::SuperkSample_t s) const;             ///< all events in input sample
  double NumOfEvents (double Er, t2kana::year1_numuosc::conventions::SuperkSample_t s) const;  ///< all events in given bin of input sample

  // Binning info
  const TAxis& ErBinning (void) const { return *fErBinning; }

  // Draw.
  // Set option to "Ereco"  to show the reconstructed Ereco distributions
  // Set option to ... to show ...
  void Draw(TVirtualPad * c, Option_t * opt = "Ereco") const;
          
private:

  void Init    (void);
  void CleanUp (void);

  TAxis * fErBinning;     ///< E_{reco} binning (both fixed and variable bin widths supported)
  TH1D *  fErDistr[kns];  ///< Reconstructed energy spectrum for each sample

ClassDef(SuperkObservation, 1)
};

} // year1_numuosc
} // t2kana

#endif // _t2kana_year1_numu_analysis_sk_expectation_h_

