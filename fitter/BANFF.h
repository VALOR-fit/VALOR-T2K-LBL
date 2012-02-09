//____________________________________________________________________________
/*!

\class    t2kana::BANFF

\brief    Class to hold analysis inputs from BANFF (Beam and ND280 Flux
          Extrapolation Task Force)

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  February 09, 2012

*/
//____________________________________________________________________________

#ifndef _t2kana_banff_h_
#define _t2kana_banff_h_

#include <TObject.h>
#include <TMatrixD.h>
#include <TMatrixTSym.h>
#include <TXMLEngine.h>

class TH1D;
class TH2D;
class TVirtualPad;

namespace t2kana {

class BANFF: public TObject
{
public:
  BANFF();   
 ~BANFF();

  // Get BANFF fit param values (nominal + tweaked)
  double NominalParamValue (PhysParamId_t p) const;
  double CurrentParamValue (PhysParamId_t p) const;
   
  // Get BANFF covariance matrix
  const TMatrixD & CovMatrix(void) const { return *fCovMtrx; }

  // Systematic variations 
  void Tweak (bool force_input_twkdials = false, TVectorD * input_twkv = 0, bool correlated = true, bool fractional = false);
   
  // Read inputs from XML file
  bool Read (TXMLEngine* xml, XMLNodePointer_t root_node);

  // Drawing options:
  // -
  void Draw(TVirtualPad * c, Option_t * opt = "") const;

  // Print
  void PrintSystematicVariations (void) const;

private:

  void Init    (void);
  void CleanUp (void);

  TVectorD *  fBestFitValues;  ///<
  TMatrixD *  fCovMtrx;        ///< BANFF covariance matrix
  TMatrixD *  fCovMtrxChoL;    ///< Triangular matrix obtained from Cholesky decomposition of BANFF covariance matrix
  bool        fIsTweaked;      ///< Was there a flux shape tweaking?
      
ClassDef(BANFF, 1)
};

} //t2kana

#endif
