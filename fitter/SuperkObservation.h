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
class TObjArray;

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
  bool Read(string config_file, int data_set_id=0, bool is_xml = true);

  // Create mock data
  bool CreateMockData(SuperkExpectation & exp, const FitParamList & fp,   bool fluctuate=true);
  bool CreateMockData(SuperkExpectation & exp, const PhysHypothesis & ph, bool fluctuate=true);

  // Get data
  double NumOfEvents (t2kana::conv::SuperkSample_t s = t2kana::conv::kSKsUndefined) const;             
  double NumOfEvents (double Er, t2kana::conv::SuperkSample_t s = t2kana::conv::kSKsUndefined) const;  

  TH1D*  ErecoSpectrum (t2kana::conv::SuperkSample_t s) const; // Returns new object, mut be deleted
  double*  ErecoSpectrumArray (t2kana::conv::SuperkSample_t s); 
  const double*  ErecoSpectrumArray (t2kana::conv::SuperkSample_t s) const; 

  // Write datasets to the input container or to the input file
  bool Write(TObjArray * data_collection, const char * tag = "data");
  bool Write(const char * filename, const char * tag = "data", const char * option_file="update");
  
private:
  
  double  fErDistrArray1Rmu   [conv::kNErecoBins1Rmu];
  double  fErDistrArray1Re    [conv::kNErecoBins1Re];
  double  fRHCErDistrArray1Rmu[conv::kRHCNErecoBins1Rmu];
  double  fRHCErDistrArray1Re [conv::kRHCNErecoBins1Re];

  void Init    (void);
  void CleanUp (void);
  
  util::ArrIdx * arr_index; //! Pointer to singleton, do not save in TObject

ClassDef(SuperkObservation, 4)
};

} // t2kana

#endif 

