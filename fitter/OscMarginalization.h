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

#include "valor/VArray/VALORArrayF.h"
#include "valor/FitCore/FitParamList.h"
#include "valor/Oscillation/NuOsc3Flavour.h"
#include "SuperkSamples.h"
#include "SuperkModes.h"

namespace valor
{

const int NFlavours = 8;
typedef enum EFlavours 
{
   kFNuMu,
   kFNuMuBar,
   kFNuE,
   kFNuEBar,
   kFOscNuMu,
   kFOscNuMuBar,
   kFOscNuE,
   kFOscNuEBar
}
Flavours_t;

class OscProb
{
public:
   OscProb(const int Ntoys);
   virtual ~OscProb();
   VALORArrayF* GetToys(conv::SuperkSample_t s);
   void Write();
   int ModeIndex(conv::SuperkSample_t s, conv::SuperkMode_t m);
   int SampleIndex(conv::SuperkSample_t s);

private:
   void InitSampleMap();
   void InitModeMaps();

private:
   const int N;
   VALORArrayF* vProb_1Re;
   VALORArrayF* vProb_1Rmu;
   VALORArrayF* vProb_RHC1Re;
   VALORArrayF* vProb_RHC1Rmu;
   VALORArrayF* vProb_MultiRe;

   std::map<conv::SuperkMode_t, int> mode_maps[conv::kNSuperkSamples];
   std::map<conv::SuperkSample_t, int> sample_map;

   ClassDef(OscProb,99)
};

namespace util
{
   /**
      This function computes oscillation probabilities for a given free parameter
      point whilst marginalizing the remaining oscillation parameters.

      \param Ntoys   The number of marginalization toys.
      \param fparm   The FitParamList that determines which oscillation parameters to use.
   */
   void CalcOscProbabilities(const int Ntoys, FitParamList* fparm, bool isbar = false);

//   void SetOscillationParameters(NuOsc3Flavour& osc, FitParamList* fparm);
//   void SetOscillationBarParameters(NuOsc3Flavour& osc, FitParamList* fparm);

}
}

#endif
