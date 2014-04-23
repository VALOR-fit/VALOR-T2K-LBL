//____________________________________________________________________________
/*!

\class     t2kana::FitDataABC

\brief     Pure abstract base class (pABC) for fit data.

\author    Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
           University of Liverpool & STFC Rutherford Appleton Laboratory

\created   March 30, 2014

\copyright (c) 2010-2014, The VALOR group

*/
//____________________________________________________________________________

#ifndef _t2kana_fit_data_abc_h_
#define _t2kana_fit_data_abc_h_

#include <string>
#include <TROOT.h>

class TObjArray;

using std::string;

namespace t2kana {

class FitModelABC;
class FitParamList;
class PhysHypothesis;

class FitDataABC
{
public:
  virtual ~FitDataABC();

  // Identify fit data with a name (to be propagated in output files, etc)
  virtual string Name(void) const = 0;  

  // Read data from the input file (ROOT, XML, other)
  // Use input option as needed to specify what exactly to read from the input file.
  // Return true if data were succesfully read out
  virtual bool Read(string file, Option_t * opt = 0) = 0;  

  // Generate mock data (toy-MC)
  // Toy-MC is generated from using the input model 
  // If a FitParamList is specified, all physics parameters and systematics are set from its values
  // If a PhysHypothesis is specified, all physics parameters are set from its values and it is
  // assumed that systematics were tweaked independently
  // If no parameters are specified then the pre-generated model will be used.
  virtual bool Generate (FitModelABC * model, const FitParamList & fp,   bool fluctuate=true) = 0;
  virtual bool Generate (FitModelABC * model, const PhysHypothesis & ph, bool fluctuate=true) = 0;
  virtual bool Generate (FitModelABC * model, bool fluctuate=true) = 0;

  // Write datasets to the input container or to the input file
  // The input tag will be appended to the return value of Name() and the name of each stored object.
  // For example, if Name() returns "t2k" and tag is "run1to4", a histogram with name "1ringmu" will be
  // saved out "t2k_1ringmu_run1to4".
  // The input write_opt option can be used to specify what exactly to write out. For example, if data are held
  // in arrays (as primary source, for CPU efficiency) and also in histograms (for their graphical representation)
  // one can choose to write out one or both. Depending on how each concrete implementation works and what
  // is needed by each analysis step, this may save time and/or space.
  // If null, everything should be written out.
  // The input file_opt option can be used to write out a number of results into the same file (e.g. by calling
  // the Write() method with file_opt="recreate" in the first instance and with file_opt="update" subsequently).
  // If null, the default is "recreate".
  virtual bool Write (TObjArray * collection, const char * tag = "data", Option_t * write_opt = 0) const = 0;
  virtual bool Write (const char * filename,  const char * tag = "data", Option_t * write_opt = 0, Option_t * file_opt = 0) const = 0;

  // Number of bins in fit data, or number of events for unbinned fits
  virtual int NumOfBins (void) const = 0;

  // Calculate the specified statistic
  // The corresponding prediction is taken from the input model.
  // Similar considerations as in Write() method, see above, apply for how this prediction is generated.
  // This method calculates the statistical term. 
  // The penalty terms is requested from the input model, which knows of all the systematic tweaks applied to it [TO BE IMPLEMENTED]
  virtual double Statistic (FitModelABC * pred, const FitParamList & fp,   Option_t * stat_opt = "likelihood-ratio") = 0;
  virtual double Statistic (FitModelABC * pred, const PhysHypothesis & ph, Option_t * stat_opt = "likelihood-ratio") = 0;
  virtual double Statistic (FitModelABC * pred, Option_t * stat_opt = "likelihood-ratio") = 0;

protected:
   FitDataABC();
};

} // t2kana

#endif

