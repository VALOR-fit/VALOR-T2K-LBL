//____________________________________________________________________________
/*!

\class    valor::SuperkObservation

\brief    

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  October 20, 2010

*/
//____________________________________________________________________________

#ifndef _valor_sk_observation_h_
#define _valor_sk_observation_h_

#include <string>

class TH1D;
class TObjArray;

#include "FitCore/FitDataABC.h"

#include "FitCore/FitParamList.h"

#include "SuperkExpectation.h"

using std::string;
using valor::conv::kns;

namespace valor {

class SuperkObservation: public FitDataABC
{
public:
  SuperkObservation();
 ~SuperkObservation();
 
  // FitDataABC interface
  string Name(void) {return string("SuperkObservation");};

  // Read data
  bool Read(string config_file, Option_t * opt = 0); // Valid options: dataset=n

  // Create mock data
  bool Generate(FitModelABC * exp, const FitParamList & fp,   bool fluctuate=true);
  bool Generate(FitModelABC * exp, const PhysHypothesis & ph, bool fluctuate=true);
  bool Generate(FitModelABC * exp, bool fluctuate=true);
  
  // Write datasets to the input container or to the input file
  bool Write (TObjArray * collection, const char * tag = "data", Option_t * write_opt = 0);
  bool Write (const char * filename,  const char * tag = "data", Option_t * write_opt = 0, Option_t * file_opt = 0);
  
  // Number of bins in fit data, or number of events for unbinned fits
  int NumOfBins(void);
  
  // Calculate the specified statistic
  // The corresponding prediction is taken from the input model.
  // Similar considerations as in Write() method, see above, apply for how this prediction is generated.
  // This method calculates the statistical term. 
  // The penalty terms is requested from the input model, which knows of all the systematic tweaks applied to it [TO BE IMPLEMENTED]
  double Statistic (FitModelABC * pred, const FitParamList & fp,   Option_t * stat_opt = "likelihood-ratio");
  double Statistic (FitModelABC * pred, const PhysHypothesis & ph, Option_t * stat_opt = "likelihood-ratio");
  double Statistic (FitModelABC * pred, Option_t * stat_opt = "likelihood-ratio");
  
  //  *** T2K LEGACY INTERFACE ***
  bool Read(string config_file, int data_set_id, bool is_xml = true);

  // Create mock data
  bool CreateMockData(SuperkExpectation & exp, const FitParamList & fp,   bool fluctuate=true) {return Generate(&exp,fp,fluctuate);};
  bool CreateMockData(SuperkExpectation & exp, const PhysHypothesis & ph, bool fluctuate=true) {return Generate(&exp,ph,fluctuate);};

  // Get data
  double NumOfEvents (valor::conv::SuperkSample_t s = valor::conv::kSKsUndefined) const;             
  double NumOfEvents (double Er, valor::conv::SuperkSample_t s = valor::conv::kSKsUndefined) const;  

  TH1D*  ErecoSpectrum (valor::conv::SuperkSample_t s) const; // Returns new object, mut be deleted
  VALORArray       *  ErecoSpectrumVArray (valor::conv::SuperkSample_t s);
  const VALORArray *  ErecoSpectrumVArray (valor::conv::SuperkSample_t s) const;
  double           *  ErecoSpectrumArray  (valor::conv::SuperkSample_t s);
  const double     *  ErecoSpectrumArray  (valor::conv::SuperkSample_t s) const;
  
private:
  
  VALORArray * fErDistrArray1Rmu   ;
  VALORArray * fErDistrArray1Re    ;
  VALORArray * fRHCErDistrArray1Rmu;
  VALORArray * fRHCErDistrArray1Re ;

  void Init    (void);
  void CleanUp (void);
  
  util::ArrIdx * arr_index; //! Pointer to singleton, do not save in TObject

ClassDef(SuperkObservation, 4)
};

} // valor

#endif 

