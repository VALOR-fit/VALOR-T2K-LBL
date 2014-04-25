//____________________________________________________________________________
/*!

\class    t2kana::VALORArray

\brief    N-dimensional array (for N<=6)

\author   Steve Dennis <stephen.dennis \at warwick.ac.uk>
          University of Warwick, Rutherford Appleton Laboratory

\created  April 2014

*/
//____________________________________________________________________________

#ifndef _t2kana_array_h_
#define _t2kana_array_h_

#include <TObject.h>
#include <vector>

namespace t2kana {

class VALORArray: public TObject
{
public:
  // Constructors
  VALORArray();
  VALORArray(const std::vector<int> dimensions);   // List of dimension sizes
  VALORArray(int nd, const int * dimensions);      // List of dimension sizes
  VALORArray(int dim1, int dim2=-1, int dim3=-1, int dim4=-1, int dim5=-1, int dim6=-1); // Dimension sizes individually
  VALORArray(const VALORArray & arr);
 ~VALORArray();
  
  // Retrives size of total array, or size of dimension
  int GetSize(int dim = -1) const;
  
  // Get value at specific point
  double GetEntry           (int ind1, int ind2=-1, int ind3=-1, int ind4=-1, int ind5=-1, int ind6=-1) const;
  double & GetEntryRef      (int ind1, int ind2=-1, int ind3=-1, int ind4=-1, int ind5=-1, int ind6=-1);
  const double & GetEntryRef(int ind1, int ind2=-1, int ind3=-1, int ind4=-1, int ind5=-1, int ind6=-1) const;
  
  // Set value at specific point
  void SetEntry(double value, int ind1, int ind2=-1, int ind3=-1, int ind4=-1, int ind5=-1, int ind6=-1);
  //double GetEntry(int dummy,...) const;
  //void SetEntry(double value, ...);
  
  int GetNDim() const;
  const double * GetArray(void) const;
  double * GetArray(void);
  
  const int * GetDimensions(void) const;
  
  void SetAll(double value);
  
private:
  int fNDim;
  int fSize;
  
  int * fDimensions ; //[fNDim] Stores size of each dimension
  int * fIncrement;   //[fNDim] Stores number to increment in main array to get +1 in this dimension
  
  double * fArray;    //[fSize] Stores actual contents of array
  
  int Address(int ind1, int ind2=-1, int ind3=-1, int ind4=-1, int ind5=-1, int ind6=-1) const;
  //int ArgsToAddress(va_list args) const;
  void Construct(int nd, const int * dimensions);
  
ClassDef(VALORArray, 1)
};

} //t2kana

#endif
