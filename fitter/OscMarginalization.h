//____________________________________________________________________________
/*!

\class    valor::OscMarginalization

\brief    Code to perform marginalization of oscillation parameters

\author   Andy Chappell <A.Chappell \at warwick.ac.uk>
          
\created  Sep 28, 2017

\copyright (c) 2010-2017, VALOR group.

*/
//____________________________________________________________________________

#ifndef _valor_oscmarginalization__h_
#define _valor_oscmarginalization__h_


#include <string>

#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>

#include "valor/VArray/VALORArrayD.h"
#include "valor/FitCore/FitParamList.h"
#include "valor/Oscillation/NuOsc3Flavour.h"

namespace valor
{
namespace util
{
   // This should probably just use the FitParamList to work out what values to generate.
   // 6D VALORArray should hold the values, 2D for the scan point, further 3D for sample,
   // Etrue and mode. Possibly 1D for the number of toys (might just store these as
   // separate VALORArrays and load from file for each toy).
   /**
      This function computes oscillation probabilities for a given free parameter
      point whilst marginalizing the remaining oscillation parameters.

      \param Ntoys   The number of marginalization toys.
      \param fparm   The FitParamList that determines which oscillation parameters to use.
   */
   void CalcOscProbabilities(const int Ntoys, FitParamList* fparm);
   void CalcOscBarProbabilities(const int Ntoys, FitParamList* fparm);

//   void SetOscillationParameters(NuOsc3Flavour& osc, FitParamList* fparm);
//   void SetOscillationBarParameters(NuOsc3Flavour& osc, FitParamList* fparm);

}
}

#endif
