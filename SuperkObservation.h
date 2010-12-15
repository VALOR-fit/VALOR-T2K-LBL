//____________________________________________________________________________
/*!

\class    t2kana::yr1numu::SuperkObservation

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
#include <vector>

class TH1D;
class TVirtualPad;

#include "year1_numu_osc_analysis/Conventions.h"
#include "year1_numu_osc_analysis/SuperkExpectation.h"
#include "year1_numu_osc_analysis/PhysParamList.h"

using std::string;
using std::vector;
using t2kana::yr1numu::conv::kns;

namespace t2kana {
namespace yr1numu {

class SuperkObservation: public TObject
{
public:
  SuperkObservation();
  SuperkObservation(const vector<t2kana::yr1numu::conv::SuperkSample_t> & samples);
 ~SuperkObservation();

  // Read data
  bool Read(string config_file);

  // Create mock data
  void CreateMockData(SuperkExpectation & exp, const PhysParamList & fp, bool fluctuate=true);

  // Get data
  double NumOfEvents (t2kana::yr1numu::conv::SuperkSample_t s = t2kana::yr1numu::conv::kSKsUndefined) const;             
  double NumOfEvents (double Er, t2kana::yr1numu::conv::SuperkSample_t s = t2kana::yr1numu::conv::kSKsUndefined) const;  

  TH1D*  ErecoSpectrum (void) const; 
  TH1D*  ErecoSpectrum (t2kana::yr1numu::conv::SuperkSample_t s) const; 

///  TH1D*  AsTH1D (t2kana::yr1numu::conv::SuperkSample_t s = t2kana::yr1numu::conv::kSKsUndefined) const; 

  // Draw.
  // Set option to "Ereco"  to show the reconstructed Ereco distributions
  // Set option to ... to show ...
  void Draw(TVirtualPad * c, Option_t * opt = "Ereco") const;

private:

  void Init    (void);
  void CleanUp (void);

  TH1D *  fErDistr[kns];  ///< Reconstructed energy spectrum for each sample
  TH1D *  fErDistrTot;    ///< Reconstructed energy spectrum for all samples

  vector<t2kana::yr1numu::conv::SuperkSample_t> fSamples;  

ClassDef(SuperkObservation, 1)
};

} // yr1numu
} // t2kana

#endif // _t2kana_year1_numu_analysis_sk_observation_h_

