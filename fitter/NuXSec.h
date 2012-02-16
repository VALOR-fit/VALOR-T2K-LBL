//____________________________________________________________________________
/*!
   
\class    t2kana::NuXSec
      
\brief    Neutrino generator / cross section info used in the oscillation fit
   
\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  October 20, 2010

*/
//____________________________________________________________________________
   
#ifndef _t2kana_nu_xsec_h_
#define _t2kana_nu_xsec_h_

#include <TObject.h>
#include <TXMLEngine.h>

#include "Conventions.h"

class TH1D;
class TVirtualPad;

namespace t2kana {

class BANFF;

class NuXSec: public TObject
{
public:
   NuXSec();
  ~NuXSec();

   // Get weight
   double SuperkXSecWght (double Etrue, t2kana::conv::SuperkSample_t s, t2kana::conv::SuperkMode_t m) const;

   // Systematic variations
   // tweak...

   // Print
   void PrintSystematicVariations (void) const;

   // Read inputs from XML file
   bool Read (TXMLEngine* xml, XMLNodePointer_t root_node);

   // BANFF fit info
   void TellMeAboutTheBANFFfit(BANFF * banff);

   // Draw inputs:
   //
   void Draw(TVirtualPad * c, Option_t * opt = "") const;
   
private:

   void Init    (void);
   void CleanUp (void);


   BANFF *  fBANFF;  ///< BANFF fit info

ClassDef(NuXSec, 4);
};

} // t2kana

#endif
